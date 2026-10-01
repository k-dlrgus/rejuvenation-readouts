"""STAGE P1 — per-gene variance decomposition on 3,323 DLPFC pseudobulks.

Model (P0 extras filled in at runtime):
  expression ~ C(cell_type) + C(Source) + log(mean UMI/nucleus) + log(mean genes/nucleus)
               + [P0 covariates] + age

Thresholds are QUANTILES of this corrected distribution, declared a priori
(PRIMARY_QUANTILES in brain_phase1_common.py). Numeric cutoffs and set sizes are
written to FINDINGS_BRAIN_PHASE1.md BEFORE any score is built.

Usage: python src/brain_phase1_p1.py
"""
from __future__ import annotations

import sys
from pathlib import Path

import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

sys.path.insert(0, str(Path(__file__).resolve().parent))
from brain_phase1_common import (  # noqa: E402
    PHASE1_DIR, PHASE1_FIG, PHASE1_SEED, Logger, HASH_POOL_LIMITATION,
    PRIMARY_QUANTILES, SWEEP_HIGH_Q, SWEEP_LOW_Q, N_PERM_P1,
    BLOOD_UA_MAX, BLOOD_UA_P99, BLOOD_UA_P50, BLOOD_UCT_P50,
    ensure_logcpm, merge_p0_into_obs, load_p0_extras, decompose_brain,
    quantiles_to_thresholds, apply_thresholds, dump_json,
    build_design,
)
from brain_phase1_findings import write_findings  # noqa: E402
from rerun_common import residualize as _resid  # noqa: E402


def permute_age_within_site(obs, rng):
    don = obs.groupby("donor").agg(age=("age", "first"), site=("Source", "first"))
    perm = don.age.copy()
    for _, idx in don.groupby("site").groups.items():
        perm.loc[idx] = rng.permutation(don.age.loc[idx].to_numpy())
    return obs.donor.map(perm).to_numpy(float)


def run():
    rng = np.random.default_rng(PHASE1_SEED)
    log = Logger(PHASE1_DIR / "p1_report.txt")
    try:
        return _run(log, rng)
    finally:
        log.close()


def _run(log, rng):
    log("=" * 100)
    log("BRAIN PHASE 1  STAGE P1 — variance decomposition.  seed =", PHASE1_SEED)
    log("=" * 100)
    log(HASH_POOL_LIMITATION)
    log(f"DECLARED quantile rule (a priori, before looking at set sizes or scores): {PRIMARY_QUANTILES}")
    log("  A: unique_age >= Q{t_age_q} AND unique_ct <= Q{k_ct_q}".format(**PRIMARY_QUANTILES))
    log("  I: unique_ct  >= Q{t_ct_q} AND unique_age <= Q{k_age_q}".format(**PRIMARY_QUANTILES))
    log("  MIXED: unique_age >= Q{t_age_q} AND unique_ct >= Q{t_ct_q}  (excluded from both scores)")
    log("Sweep is reported; selection is NOT on downstream R².")

    Y, genes, obs = ensure_logcpm(log)
    obs = merge_p0_into_obs(obs)
    extra = load_p0_extras()
    log(f"[P1] pseudobulks={len(obs):,} genes={Y.shape[1]:,} donors={obs.donor.nunique()} "
        f"types={obs.celltype.nunique()} Source={obs.Source.value_counts().to_dict()}")
    log(f"[P1] P0 extra covariates in design: {extra or '(none; depth + Source only)'}")
    n_par = 1 + (obs.celltype.nunique()) + (obs.Source.nunique() - 1) + 2 + len(extra)
    log(f"[P1] design: {obs.celltype.nunique()} cell types + {obs.Source.nunique() - 1} Source contrast "
        f"+ 2 depth + {len(extra)} P0 extras + 1 age  (~{n_par} columns) on {len(obs):,} rows")

    # donor nesting
    n_src = obs.groupby("donor").Source.nunique()
    assert (n_src == 1).all(), "a donor appears in more than one Source"
    log("[P1] every donor in exactly one Source: True")

    df = decompose_brain(Y, obs, extra_cols=extra,
                         gene_ids=genes.gene_id.to_numpy(), symbols=genes.symbol.to_numpy())
    q = [0.5, 0.8, 0.9, 0.95, 0.99, 0.999, 1.0]
    cols = ["unique_ct", "unique_age", "shared", "nuisance", "interaction", "resid", "partial_age"]
    log("\nCorrected model — distribution over %d genes (fraction of total variance):" % len(df))
    log(df[cols].quantile(q).T.to_string(float_format=lambda x: f"{x:.5f}"))
    means = df[cols].mean()
    log("Mean partition: " + ", ".join(f"{k}={v:.4f}" for k, v in means.items()))
    log(f"\nBlood within-pool unique_age (FINDINGS_RERUN.md): max={BLOOD_UA_MAX:.4f} ({BLOOD_UA_MAX:.2%}), "
        f"99th pct={BLOOD_UA_P99:.4f} ({BLOOD_UA_P99:.2%}), median={BLOOD_UA_P50:.5f}")
    log(f"Brain unique_age: max={df.unique_age.max():.5f} ({df.unique_age.max():.2%}), "
        f"99th={df.unique_age.quantile(0.99):.5f} ({df.unique_age.quantile(0.99):.2%}), "
        f"median={df.unique_age.median():.5f}")
    log(f"Brain unique_ct: max={df.unique_ct.max():.4f}, median={df.unique_ct.median():.4f} "
        f"(blood median unique_ct={BLOOD_UCT_P50:.4f})")

    log("\nunique_age: pooled (ct+age) vs Source-only vs corrected (Source+depth+P0):")
    cmp = pd.DataFrame({
        "pooled (ct+age)": df.unique_age_pooled.quantile(q),
        "source-only (ct+Source+age)": df.unique_age_source_only.quantile(q),
        "corrected (ct+Source+depth+P0+age)": df.unique_age.quantile(q),
    }).T
    log(cmp.to_string(float_format=lambda x: f"{x:.5f}"))
    cmp.to_csv(PHASE1_DIR / "p1_unique_age_quantiles_models.csv")
    log(f"Spearman(pooled, corrected unique_age) = {df.unique_age_pooled.corr(df.unique_age, method='spearman'):.3f}")
    log(f"Spearman(source-only, corrected) = {df.unique_age_source_only.corr(df.unique_age, method='spearman'):.3f}")

    log("\nTop 25 genes by CORRECTED unique_age:")
    top = df.sort_values("unique_age", ascending=False).head(25)
    log(top[["symbol", "unique_age", "unique_ct", "unique_age_pooled", "age_slope_per_year", "mean_logcpm"]]
        .to_string(float_format=lambda x: f"{x:.4f}"))

    # permutation null (orientation; thresholds are quantile-based, not FDR-based)
    log("\n" + "-" * 80)
    log(f"[P1-null] donor ages permuted WITHIN site, {N_PERM_P1} draws (orientation only)")
    B = build_design(obs, extra_cols=extra)
    X0 = _hs_safe(B)
    Yc = Y - Y.mean(0, keepdims=True)
    ss_tot = (Yc ** 2).sum(0)
    Yr = _resid(Yc, X0)
    age_r = _resid(B["age"], X0)[:, 0]
    ua_obs = (Yr.T @ age_r) ** 2 / ((age_r ** 2).sum() * ss_tot)
    assert np.allclose(ua_obs, df.unique_age.to_numpy(), atol=1e-8), "unique_age closed form mismatch"
    null_max = np.empty(N_PERM_P1)
    # store only max and count-above-t; full null is large
    t_grid = [0.001, 0.002, 0.005, 0.01, 0.02]
    counts_null = {t: [] for t in t_grid}
    for b in range(N_PERM_P1):
        a_p = permute_age_within_site(obs, rng)
        a_pr = _resid((a_p - a_p.mean())[:, None], X0)[:, 0]
        ua_n = (Yr.T @ a_pr) ** 2 / ((a_pr ** 2).sum() * ss_tot)
        null_max[b] = ua_n.max()
        for t in t_grid:
            counts_null[t].append(int((ua_n >= t).sum()))
        if (b + 1) % 50 == 0:
            log(f"   perm {b+1}/{N_PERM_P1}")
    rows = []
    for t in t_grid:
        n_obs = int((ua_obs >= t).sum())
        nn = np.array(counts_null[t], float)
        rows.append(dict(t_age=t, observed=n_obs, null_mean=float(nn.mean()),
                         null_p95=float(np.percentile(nn, 95)),
                         empirical_FDR=float(nn.mean() / max(n_obs, 1))))
    nulltab = pd.DataFrame(rows)
    log(nulltab.to_string(index=False, float_format=lambda x: f"{x:.4f}"))
    nulltab.to_csv(PHASE1_DIR / "p1_permutation_null_counts.csv", index=False)
    log(f"null max unique_age: median {np.median(null_max):.5f}, 95th {np.percentile(null_max, 95):.5f}, "
        f"max {null_max.max():.5f}")
    np.save(PHASE1_DIR / "p1_null_max_unique_age.npy", null_max)

    # ---- DECLARE THRESHOLDS from quantiles (before any score) ----
    thresh = quantiles_to_thresholds(df, PRIMARY_QUANTILES)
    log("\n" + "=" * 80)
    log("[P1] DECLARED THRESHOLDS (quantile rule applied to the corrected distribution)")
    log(f"  t_age = Q{PRIMARY_QUANTILES['t_age_q']:.2f}(unique_age) = {thresh['t_age']:.6f}")
    log(f"  k_ct  = Q{PRIMARY_QUANTILES['k_ct_q']:.2f}(unique_ct)  = {thresh['k_ct']:.6f}")
    log(f"  t_ct  = Q{PRIMARY_QUANTILES['t_ct_q']:.2f}(unique_ct)  = {thresh['t_ct']:.6f}")
    log(f"  k_age = Q{PRIMARY_QUANTILES['k_age_q']:.2f}(unique_age) = {thresh['k_age']:.6f}")
    lab = apply_thresholds(df, thresh)
    df["gene_class"] = lab
    vc = lab.value_counts()
    nA, nI, nM, nO = (int(vc.get(k, 0)) for k in ("A", "I", "MIXED", "OTHER"))
    log(f"  set sizes: A={nA}  I={nI}  MIXED={nM}  OTHER={nO}  (n genes={len(df)})")
    log("  These sizes are written to FINDINGS_BRAIN_PHASE1.md NOW, before P2 filters and before any score.")
    dump_json(PHASE1_DIR / "p1_thresholds.json", dict(
        **thresh, n_A=nA, n_I=nI, n_MIXED=nM, n_OTHER=nO, n_genes=int(len(df)),
        extra_covariates=extra, seed=PHASE1_SEED, declared_before_scores=True,
        quantile_rule=PRIMARY_QUANTILES,
        unique_age_max=float(df.unique_age.max()),
        unique_age_p99=float(df.unique_age.quantile(0.99)),
        unique_age_median=float(df.unique_age.median()),
        unique_ct_max=float(df.unique_ct.max()),
        unique_ct_median=float(df.unique_ct.median()),
        blood_unique_age_max=BLOOD_UA_MAX, blood_unique_age_p99=BLOOD_UA_P99,
    ))
    df.to_csv(PHASE1_DIR / "p1_variance_decomposition.csv")
    df.loc[lab == "A"].sort_values("unique_age", ascending=False).to_csv(PHASE1_DIR / "p1_A_genes_prefilter.csv")
    df.loc[lab == "I"].sort_values("unique_ct", ascending=False).to_csv(PHASE1_DIR / "p1_I_genes_prefilter.csv")
    df.loc[lab == "MIXED"].sort_values("unique_age", ascending=False).to_csv(PHASE1_DIR / "p1_MIXED_genes.csv")

    # sweep (not used for selection)
    sweep_rows = []
    for hq in SWEEP_HIGH_Q:
        for lq in SWEEP_LOW_Q:
            t = dict(
                t_age=float(np.quantile(df.unique_age, hq)),
                k_ct=float(np.quantile(df.unique_ct, lq)),
                t_ct=float(np.quantile(df.unique_ct, hq)),
                k_age=float(np.quantile(df.unique_age, lq)),
            )
            lab_s = apply_thresholds(df, t)
            vc_s = lab_s.value_counts()
            sweep_rows.append(dict(
                high_q=hq, low_q=lq, **t,
                n_A=int(vc_s.get("A", 0)), n_I=int(vc_s.get("I", 0)),
                n_MIXED=int(vc_s.get("MIXED", 0)), n_OTHER=int(vc_s.get("OTHER", 0)),
                is_primary=bool(hq == PRIMARY_QUANTILES["t_age_q"] and lq == PRIMARY_QUANTILES["k_ct_q"]),
            ))
    sweep = pd.DataFrame(sweep_rows)
    sweep.to_csv(PHASE1_DIR / "p1_threshold_sweep_set_sizes.csv", index=False)
    log("\nSweep set sizes (high_q x low_q) — NOT selected on R²:")
    log(sweep.pivot(index="high_q", columns="low_q", values="n_A").to_string())
    log("I-gene sizes:")
    log(sweep.pivot(index="high_q", columns="low_q", values="n_I").to_string())

    # figures
    fig, ax = plt.subplots(figsize=(6.4, 5.2))
    ax.scatter(df.unique_ct, df.unique_age, s=4, c="0.6", alpha=0.35, linewidths=0)
    ax.axhline(thresh["t_age"], c="tab:red", ls="--", lw=1, label=f"t_age Q90={thresh['t_age']:.4f}")
    ax.axhline(thresh["k_age"], c="tab:red", ls=":", lw=1, label=f"k_age Q50={thresh['k_age']:.4f}")
    ax.axvline(thresh["k_ct"], c="tab:blue", ls="--", lw=1, label=f"k_ct Q50={thresh['k_ct']:.3f}")
    ax.axvline(thresh["t_ct"], c="tab:blue", ls=":", lw=1, label=f"t_ct Q90={thresh['t_ct']:.3f}")
    ax.scatter(df.loc[lab == "A", "unique_ct"], df.loc[lab == "A", "unique_age"],
               s=10, c="tab:red", label=f"A n={nA}", zorder=3)
    ax.scatter(df.loc[lab == "I", "unique_ct"], df.loc[lab == "I", "unique_age"],
               s=8, c="tab:blue", label=f"I n={nI}", zorder=3)
    ax.scatter(df.loc[lab == "MIXED", "unique_ct"], df.loc[lab == "MIXED", "unique_age"],
               s=12, c="tab:purple", label=f"MIXED n={nM}", zorder=3)
    ax.set_xlabel("unique_ct")
    ax.set_ylabel("unique_age")
    ax.set_title("P1 corrected partition (thresholds = quantiles)")
    ax.legend(fontsize=7, loc="upper right")
    fig.tight_layout()
    fig.savefig(PHASE1_FIG / "p1_unique_age_vs_unique_ct.png", dpi=140)
    plt.close(fig)

    fig, ax = plt.subplots(figsize=(6.4, 3.8))
    ax.hist(np.log10(np.clip(df.unique_age, 1e-8, None)), bins=60, color="steelblue",
            alpha=0.85, label="DLPFC corrected")
    ax.axvline(np.log10(BLOOD_UA_MAX), c="0.3", ls="--", label=f"blood max {BLOOD_UA_MAX:.3%}")
    ax.axvline(np.log10(BLOOD_UA_P99), c="0.3", ls=":", label=f"blood 99th {BLOOD_UA_P99:.2%}")
    ax.axvline(np.log10(max(thresh["t_age"], 1e-8)), c="tab:red", ls="--", label="brain t_age")
    ax.set_xlabel("log10 unique_age")
    ax.set_ylabel("n genes")
    ax.legend(fontsize=7)
    ax.set_title("unique_age: DLPFC vs blood within-pool")
    fig.tight_layout()
    fig.savefig(PHASE1_FIG / "p1_unique_age_hist_vs_blood.png", dpi=140)
    plt.close(fig)

    # write findings NOW with declared thresholds (before P2/P3)
    write_findings(log=log)
    stamp = PHASE1_DIR / "p1_DECLARED_BEFORE_SCORES.flag"
    stamp.write_text(
        f"thresholds and set sizes written {pd.Timestamp.now(tz='UTC').isoformat()}\n"
        f"A={nA} I={nI} MIXED={nM}\n{thresh}\n",
        encoding="utf-8",
    )
    log(f"\n[P1] FINDINGS_BRAIN_PHASE1.md updated with declared thresholds. flag={stamp}")
    log("[P1] done.")
    return dict(n_A=nA, n_I=nI, n_MIXED=nM, thresh=thresh)


def _hs_safe(B):
    from brain_phase1_common import _hs
    return _hs([B["ct"], B["source"], B["depth"], B["extra"]])


if __name__ == "__main__":
    run()
