"""RERUN Stage A — re-derive the age axis WITHIN batch (A1-A4) and evaluate STOP A.

A1  per-gene partition under  expression ~ C(cell_type) + C(pool) + log(mean UMI/cell) + log(mean genes/cell) + age
    (unique_age is now within-pool, depth-adjusted age variance); distribution vs the prior pooled run; a within-pool
    donor-level permutation null for unique_age (the nominal F test is anti-conservative because a donor's 16 pseudobulks
    share its age), giving empirical FDR per threshold and per-gene q-values.
A2  gene sets on the new partition with the UNCHANGED primary thresholds (t_age=0.005, k_ct=0.20, t_ct=0.50, k_age=0.001)
    and the unchanged sweep grid. A small A set is the expected, correct result and is not grounds to relax anything.
A3  known positives CDKN2A (p16, up with age) and LRRN3 (down with age): recovered? direction? class?
A4  the prior run's confound filters (same lists, same fixed thresholds) applied to the NEW A set; survival.
STOP A: fewer than 10 A-genes survive everything -> halt.

Outputs: results/rerun/stageA_report.txt, results/rerun/stageA_STOP_verdict.txt, results/tables/rerun_stageA_*.csv,
results/figures/rerun_stageA_*.png
"""
import sys
import json
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from scipy import stats
import statsmodels.formula.api as smf

sys.path.insert(0, str(__import__("pathlib").Path(__file__).resolve().parent))
from rerun_common import (RERUN_DIR, RERUN_SEED, tab, fig, Logger, load_data, decompose_corrected, build_design_blocks,  # noqa
                          residualize, permute_age_within_pool, onehot, bh_q)
from decomp import classify, SWEEP, PRIMARY_THRESH  # noqa  (thresholds and grid reused UNCHANGED)
from confounds import build_flags  # noqa            (confound gene lists / flags reused)
from ensembl_chrom import fetch_chromosomes  # noqa

N_PERM = 1000
STOP_A_MIN_SURVIVING = 10
# A3 gating positives (expected sign of the age slope). Literature context rows are reported, not gating.
KNOWN_POSITIVES = {"CDKN2A": +1, "LRRN3": -1}
LITERATURE_CONTEXT = {  # bulk-blood / T-cell age genes (Peters et al. 2015 Nat Commun meta-analysis; Liu et al. 2009 for p16)
    "CDKN2A": +1, "LRRN3": -1, "CD248": -1, "NELL2": -1, "LEF1": -1, "CCR7": -1, "GZMH": +1, "B3GAT1": +1, "KLRG1": +1,
}
FILTER_STEPS = [  # identical to run_1b_confounds.survival_table (order and thresholds unchanged)
    ("cell-cycle (Tirosh S/G2M + proliferation list)", "flag_cc_list"),
    ("cell-cycle (data: |r|>=0.3 with cycling covariates)", "flag_cc_data"),
    ("ribosomal (RPL/RPS/MRP)", "flag_ribo"),
    ("mitochondrial (MT-)", "flag_mito"),
    ("sex-linked (chrX/chrY)", "flag_sex_chrom"),
    ("sex-DE (|d|>=0.5, p<1e-6, within cell type)", "flag_sex_data_only"),
    ("depth-correlated (|r|>=0.3 with mean counts/genes per cell)", "flag_depth"),
]


def survival_table(genes_A, F):
    fa = F.loc[genes_A].copy()
    fa["flag_sex_data_only"] = fa.flag_sex_data & ~fa.flag_sex_chrom
    n0 = len(fa)
    rows, removed = [], pd.Series(False, index=fa.index)
    for name, col in FILTER_STEPS:
        flag = fa[col].astype(bool)
        rows.append(dict(filter=name, flagged_marginal=int(flag.sum()), removed_incremental=int((flag & ~removed).sum()),
                         surviving=n0 - int((removed | flag).sum())))
        removed |= flag
    return pd.DataFrame(rows), removed


def fmt(x):
    return f"{x:.4f}"


def main():
    log = Logger(RERUN_DIR / "stageA_report.txt")
    rng = np.random.default_rng(RERUN_SEED)
    log("=" * 100)
    log("RERUN STAGE A — age axis within batch.  seed =", RERUN_SEED, "| N_PERM =", N_PERM)
    log("=" * 100)
    Y, genes, obs = load_data(log)
    gid, sym = genes.gene_id.to_numpy(), genes.symbol.to_numpy()
    ct, pool, age = obs.celltype.to_numpy(), obs.pool_number.to_numpy(), obs.age.to_numpy(float)

    # ------------------------------------------------------------------ A0. design audit (re-stated for the rerun)
    don = obs.groupby("donor").agg(age=("age", "first"), pool=("pool_number", "first"))
    f_pool = smf.ols("age ~ C(pool)", data=don.assign(pool=don.pool.astype(str))).fit()
    log(f"\n[A0] donors={len(don)}, pools={don.pool.nunique()} (donors/pool: min {don.pool.value_counts().min()}, "
        f"median {don.pool.value_counts().median():.0f}, max {don.pool.value_counts().max()})")
    log(f"[A0] age ~ C(pool): R2={f_pool.rsquared:.3f}, adj R2={f_pool.rsquared_adj:.3f}, p={f_pool.f_pvalue:.2e}  "
        f"-> within-pool age variance available to the model = {1-f_pool.rsquared:.1%} of total donor age variance")
    B = build_design_blocks(obs)
    for j, nm in enumerate(("log mean UMI/cell", "log mean genes/cell")):
        v = B["depth"][:, j]
        r2 = 1 - ((residualize(v[:, None], np.hstack([B['ones'], B['pool']]))) ** 2).sum() / ((v - v.mean()) ** 2).sum()
        log(f"[A0] R2( {nm} [within-type-centred] ~ C(pool) ) = {r2:.3f}")
    n_par = B["ct"].shape[1] + B["pool"].shape[1] + B["depth"].shape[1] + 1
    log(f"[A0] corrected model columns: 16 cell types + 74 pool contrasts + 2 depth + 1 age = {n_par} parameters on {len(obs):,} pseudobulks")

    # ------------------------------------------------------------------ A1. decomposition
    log("\n" + "=" * 100 + "\n[A1] per-gene variance partition, corrected model\n" + "=" * 100)
    df = decompose_corrected(Y, obs, gene_ids=gid, symbols=sym)
    q = [0.5, 0.9, 0.95, 0.99, 0.999, 1.0]
    log("\nCorrected model — distribution over %d genes (fraction of total variance):" % len(df))
    log(df[["unique_ct", "unique_age", "shared", "nuisance", "interaction", "resid", "partial_age"]].quantile(q).T.to_string(float_format=fmt))
    log("\nunique_age: prior pooled model vs pool-only vs corrected (pool + depth):")
    cmp = pd.DataFrame({"pooled (prior run: ct + age)": df.unique_age_pooled.quantile(q),
                        "pool-only (ct + pool + age)": df.unique_age_pool_only.quantile(q),
                        "corrected (ct + pool + depth + age)": df.unique_age.quantile(q)}).T
    log(cmp.to_string(float_format=lambda x: f"{x:.5f}"))
    cmp.to_csv(tab("stageA_unique_age_quantiles_prior_vs_corrected.csv"))
    for t in (0.001, 0.002, 0.005, 0.01, 0.02):
        log(f"   genes with unique_age >= {t:<6}: pooled={int((df.unique_age_pooled >= t).sum()):>5}  "
            f"pool-only={int((df.unique_age_pool_only >= t).sum()):>4}  corrected={int((df.unique_age >= t).sum()):>4}  "
            f"(in both pooled & corrected: {int(((df.unique_age_pooled >= t) & (df.unique_age >= t)).sum())})")
    log(f"\nSpearman(unique_age pooled, corrected) = {df.unique_age_pooled.corr(df.unique_age, method='spearman'):.3f}; "
        f"Spearman(unique_age pool-only, corrected) = {df.unique_age_pool_only.corr(df.unique_age, method='spearman'):.3f}")
    log(f"Spearman(unique_ct pooled, corrected) = {df.unique_ct_pooled.corr(df.unique_ct, method='spearman'):.4f}; "
        f"max |diff| = {(df.unique_ct_pooled - df.unique_ct).abs().max():.4f} (identity axis unchanged by the correction)")
    means = df[["unique_ct", "unique_age", "shared", "nuisance", "interaction", "resid"]].mean()
    log("Mean partition per gene: " + ", ".join(f"{k}={v:.4f}" for k, v in means.items()))
    ret = (df.unique_age / df.unique_age_pooled.replace(0, np.nan))
    log(f"Fraction of pooled age variance retained under the corrected model, genes with pooled unique_age >= 0.005 "
        f"(n={int((df.unique_age_pooled >= .005).sum())}): median {ret[df.unique_age_pooled >= .005].median():.3f}")
    log("\nTop 30 genes by CORRECTED unique_age:")
    top = df.sort_values("unique_age", ascending=False).head(30)
    log(top[["symbol", "unique_age", "unique_age_pool_only", "unique_age_pooled", "unique_ct", "interaction", "age_slope_per_year", "mean_logcpm"]]
        .to_string(float_format=fmt))

    # ------------------------------------------------------------------ A1-null. within-pool donor-level permutation
    log("\n" + "-" * 100 + f"\n[A1-null] permutation null for unique_age: donor ages permuted WITHIN pool, {N_PERM} draws\n" + "-" * 100)
    X0 = np.hstack([B["ct"], B["pool"], B["depth"]])
    Yc = Y - Y.mean(0, keepdims=True)
    ss_tot = (Yc ** 2).sum(0)
    Yr = residualize(Yc, X0)                       # Y with ct + pool + depth projected out
    age_r = residualize(B["age"], X0)[:, 0]
    ua_obs = (Yr.T @ age_r) ** 2 / ((age_r ** 2).sum() * ss_tot)
    assert np.allclose(ua_obs, df.unique_age.to_numpy(), atol=1e-9), "closed-form unique_age disagrees with QR partition"
    null = np.empty((N_PERM, Y.shape[1]), dtype=np.float32)
    for b in range(N_PERM):
        a_p = permute_age_within_pool(obs, rng)
        a_pr = residualize((a_p - a_p.mean())[:, None], X0)[:, 0]
        null[b] = (Yr.T @ a_pr) ** 2 / ((a_pr ** 2).sum() * ss_tot)
    null_max = null.max(1)
    rows = []
    for t in (0.001, 0.002, 0.005, 0.01):
        n_obs = int((ua_obs >= t).sum())
        n_null = (null >= t).sum(1)
        rows.append(dict(t_age=t, observed=n_obs, null_mean=float(n_null.mean()), null_p95=float(np.percentile(n_null, 95)),
                         null_max=int(n_null.max()), empirical_FDR=float(n_null.mean() / max(n_obs, 1))))
    nulltab = pd.DataFrame(rows)
    log(nulltab.to_string(index=False, float_format=lambda x: f"{x:.4f}"))
    nulltab.to_csv(tab("stageA_permutation_null_counts.csv"), index=False)
    fwer95 = float(np.percentile(null_max, 95))
    log(f"null max unique_age over genes: median {np.median(null_max):.5f}, 95th pct {fwer95:.5f}, max {null_max.max():.5f}  "
        f"-> genes exceeding the FWER-5% threshold: {int((ua_obs > fwer95).sum())}")
    # per-gene empirical p (with +1 correction) and Benjamini-Hochberg q
    p_emp = (1 + (null >= ua_obs[None, :]).sum(0)) / (N_PERM + 1)
    q_bh = bh_q(p_emp)
    df["age_p_perm"] = p_emp
    df["age_q_BH_perm"] = q_bh
    df["age_q_BH_nominal"] = bh_q(df.age_p_nominal.to_numpy())
    for a in (0.05, 0.01):
        log(f"genes with permutation BH q < {a}: {int((q_bh < a).sum())}   (nominal-F BH q < {a}: {int((df.age_q_BH_nominal < a).sum())} — "
            f"anti-conservative, orientation only)")
    log(f"genes at the permutation floor p = 1/{N_PERM+1} : {int((p_emp <= 1/(N_PERM+1)).sum())}")
    np.save(RERUN_DIR / "stageA_null_max_unique_age.npy", null_max)

    # ------------------------------------------------------------------ A2. gene sets (thresholds UNCHANGED)
    log("\n" + "=" * 100 + f"\n[A2] gene sets on the corrected partition. PRIMARY thresholds (unchanged): {PRIMARY_THRESH}\n" + "=" * 100)
    Pth = PRIMARY_THRESH
    lab = classify(df, **Pth)
    df["gene_class_rerun"] = lab
    vc = lab.value_counts()
    log("class counts (primary): " + ", ".join(f"{k}={int(vc.get(k, 0))}" for k in ("A", "I", "MIXED", "OTHER")))
    A = df.index[lab == "A"]
    log(f"A-genes (n={len(A)}): " + (", ".join(df.loc[A].sort_values('unique_age', ascending=False).symbol) or "(none)"))
    M = df.index[lab == "MIXED"]
    log(f"MIXED (n={len(M)}): " + (", ".join(df.loc[M].sort_values('unique_age', ascending=False).symbol) or "(none)"))
    I = df.index[lab == "I"]
    log(f"I-genes: n={len(I)}")
    # prior sets for comparison (recomputed from the pooled columns with the same thresholds; no dependency on prior files)
    lab_prior = classify(df, **Pth, age_col="unique_age_pooled", ct_col="unique_ct_pooled")
    df["gene_class_prior_model"] = lab_prior
    A_prior = df.index[lab_prior == "A"]
    log(f"\nprior-model A-genes recomputed here: n={len(A_prior)} (FINDINGS.md: 140); overlap with new A: {len(set(A) & set(A_prior))}; "
        f"prior A-genes now classed as: {df.loc[A_prior, 'gene_class_rerun'].value_counts().to_dict()}")
    log(f"prior-model I-genes: n={int((lab_prior=='I').sum())}; new I-genes: {len(I)}; overlap: {len(set(I) & set(df.index[lab_prior=='I']))}")
    log("What happened to the genes that pass t_age in the corrected model but are not A:")
    hi = df[df.unique_age >= Pth["t_age"]]
    log(f"   {len(hi)} genes have corrected unique_age >= {Pth['t_age']}; of these unique_ct <= {Pth['k_ct']} (A): {int((hi.unique_ct <= Pth['k_ct']).sum())}, "
        f"unique_ct in ({Pth['k_ct']}, {Pth['t_ct']}) (neither): {int(((hi.unique_ct > Pth['k_ct']) & (hi.unique_ct < Pth['t_ct'])).sum())}, "
        f"unique_ct >= {Pth['t_ct']} (MIXED): {int((hi.unique_ct >= Pth['t_ct']).sum())}")
    log(hi.sort_values("unique_age", ascending=False)[["symbol", "unique_age", "unique_ct", "age_slope_per_year", "gene_class_rerun", "age_q_BH_perm"]]
        .to_string(float_format=fmt))
    # where on the identity axis do the genes with a detectable within-pool age association sit? (marginal of the A/I/MIXED grid)
    bins = [-1, Pth["k_ct"], Pth["t_ct"], 1.01]
    labels_ = [f"unique_ct <= {Pth['k_ct']} (A-eligible)", f"{Pth['k_ct']} < unique_ct < {Pth['t_ct']} (neither)", f"unique_ct >= {Pth['t_ct']} (MIXED-eligible)"]
    rows = []
    for name, mask in (("all expressed genes", np.ones(len(df), bool)), ("permutation q<0.05", (df.age_q_BH_perm < 0.05).to_numpy()),
                       ("unique_age >= 0.002", (df.unique_age >= 0.002).to_numpy()), ("unique_age >= 0.005", (df.unique_age >= 0.005).to_numpy())):
        c = pd.cut(df.unique_ct[mask], bins=bins, labels=labels_).value_counts().reindex(labels_)
        rows.append(dict(gene_set=name, n=int(mask.sum()), **{k: int(v) for k, v in c.items()},
                         **{f"frac_{i}": float(v / max(mask.sum(), 1)) for i, v in enumerate(c.values)}))
    idax = pd.DataFrame(rows)
    log("\nIdentity-axis position of genes with a within-pool age association (counts; fractions in the csv):")
    log(idax[["gene_set", "n"] + labels_].to_string(index=False))
    idax.to_csv(tab("stageA_age_genes_by_identity_axis.csv"), index=False)
    # sweep
    rows = []
    for t_age in SWEEP["t_age"]:
        for k_ct in SWEEP["k_ct"]:
            for t_ct in SWEEP["t_ct"]:
                for k_age in SWEEP["k_age"]:
                    l_ = classify(df, t_age, k_ct, t_ct, k_age)
                    v_ = l_.value_counts()
                    rows.append(dict(t_age=t_age, k_ct=k_ct, t_ct=t_ct, k_age=k_age, n_A=int(v_.get("A", 0)), n_I=int(v_.get("I", 0)),
                                     n_MIXED=int(v_.get("MIXED", 0)), n_OTHER=int(v_.get("OTHER", 0))))
    sweep = pd.DataFrame(rows)
    sweep.to_csv(tab("stageA_threshold_sweep_set_sizes.csv"), index=False)
    log(f"\nA-gene set size by (t_age, k_ct) [corrected model; I thresholds fixed t_ct={Pth['t_ct']}, k_age={Pth['k_age']}]:")
    log(sweep[(sweep.t_ct == Pth["t_ct"]) & (sweep.k_age == Pth["k_age"])].pivot(index="t_age", columns="k_ct", values="n_A").to_string())
    log("   (prior pooled model, same grid: 328/700/872/1096 | 65/140/170/242 | 13/30/39/64 | 6/10/11/19)")
    log(f"MIXED by (t_age, t_ct) [k_ct={Pth['k_ct']}, k_age={Pth['k_age']}]:")
    log(sweep[(sweep.k_ct == Pth["k_ct"]) & (sweep.k_age == Pth["k_age"])].pivot(index="t_age", columns="t_ct", values="n_MIXED").to_string())
    log(f"I-genes by (t_ct, k_age) [t_age={Pth['t_age']}, k_ct={Pth['k_ct']}]:")
    log(sweep[(sweep.t_age == Pth["t_age"]) & (sweep.k_ct == Pth["k_ct"])].pivot(index="t_ct", columns="k_age", values="n_I").to_string())

    # ------------------------------------------------------------------ A3. known positives
    log("\n" + "=" * 100 + "\n[A3] known positives: recovered with the expected direction?\n" + "=" * 100)
    df["rank_unique_age"] = df.unique_age.rank(ascending=False, method="min").astype(int)
    df["rank_unique_age_pooled"] = df.unique_age_pooled.rank(ascending=False, method="min").astype(int)
    rows = []
    for s, sign in LITERATURE_CONTEXT.items():
        if s not in set(sym):
            rows.append(dict(symbol=s, expected_sign=sign, gating=s in KNOWN_POSITIVES, direction_ok=False, note="not in expressed gene set")); continue
        r = df[df.symbol == s].iloc[0]
        rows.append(dict(symbol=s, expected_sign=sign, slope_corrected=r.age_slope_per_year, slope_pooled=r.age_slope_per_year_pooled,
                         direction_ok=bool(np.sign(r.age_slope_per_year) == sign), unique_age_corrected=r.unique_age,
                         rank_corrected=int(r.rank_unique_age), unique_age_pooled=r.unique_age_pooled, rank_pooled=int(r.rank_unique_age_pooled),
                         unique_ct=r.unique_ct, q_perm=r.age_q_BH_perm, gene_class=r.gene_class_rerun, gating=s in KNOWN_POSITIVES))
    pos = pd.DataFrame(rows)
    log(pos.to_string(index=False, float_format=fmt))
    pos.to_csv(tab("stageA_known_positives.csv"), index=False)
    # per-cell-type within-pool correlation for the two gating positives
    log("\nPer-cell-type within-pool Pearson r(expression, age) [both pool-centred within the cell type] for the gating positives:")
    pct_rows = []
    for s in KNOWN_POSITIVES:
        j = int(np.flatnonzero(sym == s)[0])
        for c in sorted(np.unique(ct)):
            m_ = ct == c
            Xp = np.hstack([np.ones((m_.sum(), 1)), onehot(pool[m_])[:, 1:]])
            yr = residualize(Y[m_, j][:, None], Xp)[:, 0]
            ar = residualize(age[m_][:, None], Xp)[:, 0]
            r_, p_ = stats.pearsonr(yr, ar)
            pct_rows.append(dict(symbol=s, celltype=c, n=int(m_.sum()), r_within_pool=r_, p=p_))
    pct = pd.DataFrame(pct_rows)
    log(pct.pivot(index="celltype", columns="symbol", values="r_within_pool").to_string(float_format=lambda x: f"{x:+.3f}"))
    pct.to_csv(tab("stageA_known_positives_per_celltype.csv"), index=False)
    a3_ok = all(pos[pos.gating].direction_ok) and all(pos[pos.gating].q_perm < 0.05)
    a3_top = pos[pos.gating].rank_corrected.max()
    log(f"\nA3 verdict: gating positives recovered with expected direction and permutation q<0.05: {a3_ok}; worst rank among them = {a3_top} of {len(df)}. "
        f"Classes: {dict(zip(pos[pos.gating].symbol, pos[pos.gating].gene_class))}")
    if not a3_ok:
        log("   !! at least one gating positive lost -> correction may be over-aggressive; investigate before proceeding (reported in FINDINGS_RERUN.md)")

    # ------------------------------------------------------------------ A4. confound filters (lists/thresholds unchanged)
    log("\n" + "=" * 100 + "\n[A4] prior-run confound filters applied to the NEW A-gene set\n" + "=" * 100)
    chrom = fetch_chromosomes(list(gid), verbose=False).set_index("gene_id").chrom
    obs_f = obs.copy()
    F = build_flags(Y, obs_f, gid, sym, chrom)
    F.to_csv(tab("stageA_gene_flags.csv"))
    log("genome-wide flag counts (all expressed genes): " + ", ".join(
        f"{c}={int(F[c].sum())}" for c in ("flag_cc_list", "flag_cc_data", "flag_ribo", "flag_mito", "flag_sex_chrom", "flag_sex_data", "flag_depth")))
    if len(A):
        surv_tab, removed = survival_table(A, F)
        log(f"\nsequential survival of the {len(A)} primary A-genes:")
        log(surv_tab.to_string(index=False))
        surv = df.loc[A][~removed.to_numpy()]
        flagged = df.loc[A][removed.to_numpy()]
        flag_cols = ("flag_cc_list", "flag_cc_data", "flag_ribo", "flag_mito", "flag_sex_chrom", "flag_sex_data", "flag_depth")
        log(f"\nremoved ({len(flagged)}): " + (", ".join(
            f"{r.symbol}[{','.join(c.replace('flag_', '') for c in flag_cols if bool(F.loc[g, c]))}]" for g, r in flagged.iterrows()) or "(none)"))
        log(f"SURVIVING A-genes ({len(surv)}): " + (", ".join(f"{r.symbol}({'+' if r.age_slope_per_year>0 else '-'}; ua={r.unique_age:.4f}; ct={r.unique_ct:.3f}; q={r.age_q_BH_perm:.3f})"
                                                        for _, r in surv.sort_values('unique_age', ascending=False).iterrows()) or "(none)"))
    else:
        surv_tab = pd.DataFrame(columns=["filter", "flagged_marginal", "removed_incremental", "surviving"])
        removed = pd.Series(dtype=bool)
        surv = df.iloc[0:0]
        log("primary A-gene set is EMPTY; nothing to filter.")
    surv_tab.to_csv(tab("stageA_Agene_survival_primary.csv"), index=False)
    n_surv = len(surv)
    # sweep survival
    rows = []
    for t_age in SWEEP["t_age"]:
        for k_ct in SWEEP["k_ct"]:
            l_ = classify(df, t_age, k_ct, Pth["t_ct"], Pth["k_age"])
            A_ = l_.index[l_ == "A"]
            if len(A_) == 0:
                rows.append(dict(t_age=t_age, k_ct=k_ct, n_A=0, n_removed=0, n_surviving=0, n_surviving_q05=0)); continue
            _, rem_ = survival_table(A_, F)
            sv = df.loc[A_][~rem_.to_numpy()]
            rows.append(dict(t_age=t_age, k_ct=k_ct, n_A=len(A_), n_removed=int(rem_.sum()), n_surviving=len(sv),
                             n_surviving_q05=int((sv.age_q_BH_perm < 0.05).sum())))
    sw = pd.DataFrame(rows)
    sw.to_csv(tab("stageA_Agene_survival_sweep.csv"), index=False)
    log("\nA-genes surviving all confound filters, across the sweep (t_age rows x k_ct columns):")
    log(sw.pivot(index="t_age", columns="k_ct", values="n_surviving").to_string())
    log("   of which permutation q<0.05:")
    log(sw.pivot(index="t_age", columns="k_ct", values="n_surviving_q05").to_string())

    # per-gene table for the primary sets
    df["A_survives_filters"] = False
    if len(A):
        df.loc[surv.index, "A_survives_filters"] = True
    out_cols = ["symbol", "gene_class_rerun", "gene_class_prior_model", "unique_age", "unique_age_pool_only", "unique_age_pooled", "unique_ct",
                "unique_ct_pooled", "shared", "nuisance", "interaction", "resid", "partial_age", "age_slope_per_year", "age_slope_per_year_pooled",
                "age_F", "age_p_nominal", "age_q_BH_nominal", "age_p_perm", "age_q_BH_perm", "rank_unique_age", "rank_unique_age_pooled",
                "mean_logcpm", "A_survives_filters"]
    df[out_cols].to_csv(tab("stageA_variance_decomposition.csv"))
    df.loc[df.gene_class_rerun.isin(["A", "MIXED"]), out_cols].join(
        F[["flag_cc_list", "flag_cc_data", "flag_ribo", "flag_mito", "flag_sex_chrom", "flag_sex_data", "flag_depth"]]).to_csv(tab("stageA_A_and_MIXED_genes_primary.csv"))

    # ------------------------------------------------------------------ STOP A
    fired = n_surv < STOP_A_MIN_SURVIVING
    lines = ["=" * 100, f"STOP CONDITION A — VERDICT  (primary thresholds {Pth}; filters = prior run's, unchanged)", "=" * 100,
             f"corrected-model A-genes at primary thresholds: {len(A)}",
             f"A-genes surviving all confound filters:       {n_surv}   (STOP A fires if < {STOP_A_MIN_SURVIVING})",
             f"of these with within-pool permutation q<0.05:  {int((surv.age_q_BH_perm < 0.05).sum()) if len(surv) else 0}",
             f"genes with unique_age >= t_age={Pth['t_age']} regardless of cell-type variance: {int((df.unique_age >= Pth['t_age']).sum())} "
             f"(permutation null mean {nulltab.loc[nulltab.t_age == Pth['t_age'], 'null_mean'].iloc[0]:.2f})",
             f"genes with ANY within-pool age association (permutation BH q<0.05): {int((q_bh < 0.05).sum())}",
             f"known positives: " + "; ".join(f"{r.symbol} slope {r.slope_corrected:+.4f}/yr (expected {'+' if r.expected_sign>0 else '-'}), "
                                             f"rank {r.rank_corrected}, unique_ct {r.unique_ct:.3f}, class {r.gene_class}" for _, r in pos[pos.gating].iterrows()),
             "", f"VERDICT: STOP CONDITION A {'FIRED' if fired else 'NOT fired'}."]
    if fired:
        lines += ["  Fewer than 10 A-genes survive. Under the corrected within-pool model the within-cell-type transcriptomic age signal in PBMC is",
                  "  too sparse to build a gene-set-restricted age score from at the pre-registered thresholds. Per the operating mode the rerun halts",
                  "  here: Stages B, C and D are NOT run. Thresholds were not relaxed. See FINDINGS_RERUN.md."]
    else:
        lines += ["  >= 10 A-genes survive; proceeding to Stage B."]
    lines.append("=" * 100)
    txt = "\n".join(lines)
    log("\n" + txt)
    (RERUN_DIR / "stageA_STOP_verdict.txt").write_text(txt, encoding="utf-8")
    json.dump(dict(fired=bool(fired), n_A_primary=int(len(A)), n_surviving=int(n_surv), thresholds=Pth, seed=RERUN_SEED, n_perm=N_PERM,
                   surviving_symbols=surv.symbol.tolist(), A_symbols=df.loc[A].symbol.tolist(), a3_ok=bool(a3_ok)),
              open(RERUN_DIR / "stageA_summary.json", "w"), indent=1)

    # ------------------------------------------------------------------ figures
    make_figures(df, A, A_prior, M, surv, null, null_max, nulltab, sw, surv_tab, pos, Pth)
    log.close()
    return df, fired


def make_figures(df, A, A_prior, M, surv, null, null_max, nulltab, sw, surv_tab, pos, Pth):
    # 1. pooled vs corrected unique_age
    figu, axes = plt.subplots(1, 2, figsize=(13, 5.2))
    ax = axes[0]
    x = df.unique_age_pooled.clip(lower=1e-6); y = df.unique_age.clip(lower=1e-6)
    ax.scatter(x, y, s=3, alpha=0.25, c="grey", label="all genes")
    ax.scatter(x.loc[A_prior], y.loc[A_prior], s=9, c="tab:red", label=f"prior-model A-genes (n={len(A_prior)})")
    if len(A):
        ax.scatter(x.loc[A], y.loc[A], s=22, c="tab:blue", marker="s", label=f"corrected-model A-genes (n={len(A)})")
    for s in ("CDKN2A", "LRRN3"):
        r = df[df.symbol == s].iloc[0]
        ax.annotate(s, (max(r.unique_age_pooled, 1e-6), max(r.unique_age, 1e-6)), fontsize=8, xytext=(4, 4), textcoords="offset points")
    lim = [1e-6, 0.05]
    ax.plot(lim, lim, "k--", lw=0.7, label="y = x"); ax.plot(lim, [0.77 * l for l in lim], "b:", lw=0.8, label="y = 0.77 x (real-effect expectation)")
    ax.axhline(Pth["t_age"], color="tab:red", lw=0.6, ls="--"); ax.axvline(Pth["t_age"], color="tab:red", lw=0.6, ls="--")
    ax.set_xscale("log"); ax.set_yscale("log"); ax.set_xlim(lim); ax.set_ylim(lim)
    ax.set_xlabel("unique age variance — prior pooled model (ct + age)"); ax.set_ylabel("unique age variance — corrected (ct + pool + depth + age)")
    ax.set_title("per-gene age variance before vs after the batch/depth correction"); ax.legend(fontsize=7, loc="upper left")
    ax = axes[1]
    for col, lab_, c in (("unique_age_pooled", "pooled (prior run)", "tab:red"), ("unique_age_pool_only", "pool only", "tab:orange"),
                         ("unique_age", "pool + depth (corrected)", "tab:blue")):
        v = np.sort(df[col].clip(lower=1e-7).to_numpy())[::-1]
        ax.plot(np.arange(1, len(v) + 1), v, label=lab_, color=c, lw=1.2)
    nm = np.sort(null.mean(0))[::-1]
    ax.plot(np.arange(1, len(nm) + 1), np.clip(nm, 1e-7, None), color="k", lw=0.9, ls=":", label="permutation null (mean over draws)")
    ax.axhline(Pth["t_age"], color="grey", lw=0.6, ls="--"); ax.text(1.2, Pth["t_age"] * 1.1, f"t_age = {Pth['t_age']}", fontsize=7)
    ax.set_xscale("log"); ax.set_yscale("log"); ax.set_ylim(1e-5, 0.05); ax.set_xlim(1, len(df))
    ax.set_xlabel("gene rank"); ax.set_ylabel("unique age variance"); ax.set_title("rank profile of per-gene age variance"); ax.legend(fontsize=7)
    figu.tight_layout(); figu.savefig(fig("stageA_pooled_vs_corrected_unique_age.png"), dpi=140); plt.close(figu)

    # 2. joint distribution
    figu, ax = plt.subplots(figsize=(7.5, 5.2))
    ax.scatter(df.unique_ct, df.unique_age.clip(lower=1e-6), s=3, alpha=0.25, c="grey")
    if len(A):
        ax.scatter(df.loc[A].unique_ct, df.loc[A].unique_age, s=22, c="tab:blue", marker="s", label=f"A (n={len(A)})")
    if len(M):
        ax.scatter(df.loc[M].unique_ct, df.loc[M].unique_age, s=14, c="tab:purple", label=f"MIXED (n={len(M)})")
    for _, r in pos.dropna(subset=["unique_age_corrected"]).iterrows():
        ax.scatter([r.unique_ct], [max(r.unique_age_corrected, 1e-6)], s=30, facecolors="none", edgecolors="k")
        ax.annotate(r.symbol, (r.unique_ct, max(r.unique_age_corrected, 1e-6)), fontsize=7, xytext=(3, 3), textcoords="offset points")
    for t in SWEEP["t_age"]:
        ax.axhline(t, color="tab:red", lw=0.5, ls="--")
    for k in SWEEP["k_ct"]:
        ax.axvline(k, color="tab:blue", lw=0.5, ls=":")
    ax.axvline(Pth["t_ct"], color="tab:purple", lw=0.6, ls="-.")
    ax.set_yscale("log"); ax.set_ylim(1e-5, 0.05)
    ax.set_xlabel("unique cell-type variance"); ax.set_ylabel("unique WITHIN-POOL age variance (log)")
    ax.set_title("corrected model: identity vs age axis per gene (circled = literature age genes)"); ax.legend(fontsize=8)
    figu.tight_layout(); figu.savefig(fig("stageA_joint_variance_distribution.png"), dpi=140); plt.close(figu)

    # 3. permutation null
    figu, axes = plt.subplots(1, 2, figsize=(12, 4.6))
    ax = axes[0]
    ts = nulltab.t_age.tolist()
    counts_null = [(null >= t).sum(1) for t in ts]
    ax.boxplot(counts_null, positions=range(len(ts)), widths=0.5, showfliers=True)
    ax.scatter(range(len(ts)), nulltab.observed, c="tab:red", zorder=5, label="observed")
    ax.set_xticks(range(len(ts))); ax.set_xticklabels([str(t) for t in ts]); ax.set_yscale("symlog", linthresh=1); ax.set_ylim(bottom=0)
    ax.set_xlabel("t_age"); ax.set_ylabel("# genes with unique_age >= t_age"); ax.set_title("observed vs within-pool permutation null"); ax.legend()
    ax = axes[1]
    ax.hist(null_max, bins=40, color="lightgrey", label="null: max over genes per draw")
    ax.axvline(df.unique_age.max(), color="tab:red", label=f"observed max = {df.unique_age.max():.4f} ({df.sort_values('unique_age').symbol.iloc[-1]})")
    ax.axvline(np.percentile(null_max, 95), color="k", ls="--", label="null 95th pct (FWER 5%)")
    ax.set_xlabel("unique age variance"); ax.set_title("family-wise null of the maximum"); ax.legend(fontsize=8)
    figu.tight_layout(); figu.savefig(fig("stageA_permutation_null.png"), dpi=140); plt.close(figu)

    # 4. survival
    figu, axes = plt.subplots(1, 2, figsize=(12, 4.6))
    ax = axes[0]
    if len(surv_tab):
        ax.bar(range(len(surv_tab)), surv_tab.surviving, color="steelblue")
        ax.set_xticks(range(len(surv_tab))); ax.set_xticklabels([s.split(" (")[0] for s in surv_tab["filter"]], rotation=35, ha="right", fontsize=8)
        ax.axhline(len(A), color="k", ls="--", lw=0.8)
    ax.axhline(STOP_A_MIN_SURVIVING, color="tab:red", ls=":", lw=1.0, label="STOP A line (10)")
    ax.set_ylabel("A-genes surviving (cumulative)"); ax.set_title(f"primary corrected-model A-genes n={len(A)} -> surviving {len(surv)}"); ax.legend()
    ax = axes[1]
    piv = sw.pivot(index="t_age", columns="k_ct", values="n_surviving")
    im = ax.imshow(piv.to_numpy(), cmap="Blues", aspect="auto")
    for i in range(piv.shape[0]):
        for j in range(piv.shape[1]):
            ax.text(j, i, str(piv.iat[i, j]), ha="center", va="center", fontsize=9, color="k")
    ax.set_xticks(range(piv.shape[1])); ax.set_xticklabels(piv.columns); ax.set_yticks(range(piv.shape[0])); ax.set_yticklabels(piv.index)
    ax.set_xlabel("k_ct (max cell-type variance)"); ax.set_ylabel("t_age (min within-pool age variance)")
    ax.set_title("A-genes surviving all filters across the sweep (primary = 0.005 / 0.2)")
    figu.tight_layout(); figu.savefig(fig("stageA_Agene_survival.png"), dpi=140); plt.close(figu)


if __name__ == "__main__":
    main()
