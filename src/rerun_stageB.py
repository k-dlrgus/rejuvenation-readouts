"""RERUN Stage B — the compositional question. Measurement only; no stop condition; gene sets and thresholds untouched.

B1  per-donor cell-type proportions (ALL cells of the raw h5ad, not the >=20-cell pseudobulk groups) -> CLR -> age ~ C(pool) + CLR.
    How much of the WITHIN-POOL donor age variance does composition explain (in-sample, with a within-pool permutation null)?
B2  pool-grouped 5-fold CV (HARD RULE 1 asserted), 10 repeats: age predicted from (a) composition alone, (b) within-cell-type
    expression, all 13,904 genes, unrestricted kernel ridge — per cell type and all types jointly, (c) 2-df combination of the
    out-of-fold predictions of (a) and (b). Primary metric: WITHIN-POOL R^2 (features pool-centred, target pool-centred, predictions
    re-centred within each held-out pool). Secondary: the naive pool-grouped R^2 on raw age (includes whatever between-pool
    batch/depth structure transfers to unseen pools). CIs: pool-cluster bootstrap of OOF predictions + spread over CV repeats. One
    within-pool age shuffle per model as a pipeline sanity check.
B3  of the summed within-pool per-gene age variance (Stage A unique_age), what fraction sits in A vs MIXED vs I vs OTHER (Stage A
    primary classes) and in the three identity bins; also in absolute variance units and restricted to permutation q<0.05 genes.
Outputs: results/rerun/stageB_report.txt, results/tables/rerun_stageB_*.csv, results/figures/rerun_stageB_*.png
"""
import sys
import json
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from scipy import stats
import statsmodels.api as sm

sys.path.insert(0, str(__import__("pathlib").Path(__file__).resolve().parent))
from rerun_common import RERUN_DIR, RERUN_SEED, tab, fig, Logger, load_data, permute_age_within_pool, assert_disjoint  # noqa
from config import PRIMARY  # noqa
from run_pseudobulk_onek1k import CELLTYPE_MAP  # noqa

N_REPEATS = 10
N_OUTER = 5
N_INNER = 3
N_BOOT = 2000
N_PERM_B1 = 1000
ALPHAS = np.array([0.01, 0.03, 0.1, 0.3, 1, 3, 10, 30, 100, 300, 1000])
NON_PBMC = {"Doublet", "Platelet", "Eryth", "HSPC"}  # excluded from the composition denominator (logged)


# ----------------------------------------------------------------------------- helpers
def pool_center(M, pool):
    """Subtract, per column, the mean over rows sharing a pool (unsupervised; uses no labels)."""
    M = np.asarray(M, float)
    out = M.copy()
    for p in np.unique(pool):
        m = pool == p
        out[m] -= M[m].mean(0, keepdims=True)
    return out


def within_pool_r2_mae(y_c, pred, pool):
    """y_c is already pool-centred; predictions are re-centred within each pool (unsupervised) before scoring."""
    pred_c = pool_center(pred[:, None], pool)[:, 0]
    ss_res = ((y_c - pred_c) ** 2).sum()
    ss_tot = (y_c ** 2).sum()
    return 1 - ss_res / ss_tot, np.abs(y_c - pred_c).mean()


def raw_r2_mae(y, pred):
    return 1 - ((y - pred) ** 2).sum() / ((y - y.mean()) ** 2).sum(), np.abs(y - pred).mean()


def make_pool_folds(pools_all, n_folds, rng):
    pools = np.array(sorted(pools_all))
    rng.shuffle(pools)
    return {p: i % n_folds for i, p in enumerate(pools)}


def kernel_ridge_oof(K, y, pool, donor, fold_of_pool, alphas, rng, n_inner=N_INNER, obs_for_assert=None):
    """Out-of-fold kernel-ridge predictions with pool-grouped outer folds and pool-grouped inner alpha selection.
    K: n x n kernel over rows; y: target (n,); pool: (n,); returns (pred_oof, chosen_alpha_per_fold)."""
    n = len(y)
    fold = np.array([fold_of_pool[p] for p in pool])
    pred = np.full(n, np.nan)
    chosen = []
    for k in np.unique(fold):
        te = np.flatnonzero(fold == k)
        tr = np.flatnonzero(fold != k)
        # HARD RULE 1 (donor AND pool disjoint)
        assert not (set(donor[tr]) & set(donor[te])), "donor leakage"
        assert not (set(pool[tr]) & set(pool[te])), "pool leakage"
        # inner pool-grouped CV over alphas
        tr_pools = np.array(sorted(set(pool[tr])))
        rng.shuffle(tr_pools)
        inner_of_pool = {p: i % n_inner for i, p in enumerate(tr_pools)}
        inner_fold = np.array([inner_of_pool[p] for p in pool[tr]])
        score = np.zeros(len(alphas))
        for ki in range(n_inner):
            ite = tr[inner_fold == ki]
            itr = tr[inner_fold != ki]
            mu_i = y[itr].mean()                      # intercept = training mean (kernel ridge is fit without a bias term)
            S, U = np.linalg.eigh(K[np.ix_(itr, itr)])
            Uy = U.T @ (y[itr] - mu_i)
            Kte = K[np.ix_(ite, itr)]
            for a_i, a in enumerate(alphas):
                dual = U @ (Uy / (S + a))
                p_ = mu_i + Kte @ dual
                score[a_i] += -((y[ite] - p_) ** 2).sum()
        a_best = alphas[int(np.argmax(score))]
        chosen.append(float(a_best))
        mu = y[tr].mean()
        S, U = np.linalg.eigh(K[np.ix_(tr, tr)])
        dual = U @ ((U.T @ (y[tr] - mu)) / (S + a_best))
        pred[te] = mu + K[np.ix_(te, tr)] @ dual
    assert not np.isnan(pred).any()
    return pred, chosen


def bootstrap_ci(y, pred, pool, rng, n_boot=N_BOOT, within=True):
    """Pool-cluster bootstrap (resample POOLS with replacement; donors within a pool share CV fold and batch, so the pool is the
    exchangeable unit). Predictions are re-centred within pool ONCE on the full data (within-pool version), not per resample."""
    pools = np.unique(pool)
    idx_by_pool = [np.flatnonzero(pool == p) for p in pools]
    pred_use = pool_center(pred[:, None], pool)[:, 0] if within else pred
    vals = np.empty(n_boot)
    for b in range(n_boot):
        sel = rng.integers(0, len(pools), len(pools))
        idx = np.concatenate([idx_by_pool[s] for s in sel])
        yy, pp = y[idx], pred_use[idx]
        if within:
            vals[b] = 1 - ((yy - pp) ** 2).sum() / (yy ** 2).sum()
        else:
            vals[b] = 1 - ((yy - pp) ** 2).sum() / ((yy - yy.mean()) ** 2).sum()
    return float(np.percentile(vals, 2.5)), float(np.percentile(vals, 97.5))


# ----------------------------------------------------------------------------- main
def main():
    log = Logger(RERUN_DIR / "stageB_report.txt")
    rng = np.random.default_rng(RERUN_SEED)
    log("=" * 100)
    log("RERUN STAGE B — the compositional question.  seed =", RERUN_SEED, f"| {N_REPEATS} repeats x {N_OUTER}-fold pool-grouped CV, {N_INNER}-fold inner")
    log("=" * 100)
    Y, genes, obs = load_data(log)
    sym = genes.symbol.to_numpy()
    donors = obs.groupby("donor").agg(age=("age", "first"), pool=("pool_number", "first"), sex=("sex", "first")).sort_index()
    donors["age_pool_mean"] = donors.groupby("pool").age.transform("mean")
    donors["age_c"] = donors.age - donors.age_pool_mean
    r2_pool = 1 - (donors.age_c ** 2).sum() / ((donors.age - donors.age.mean()) ** 2).sum()
    log(f"[B0] {len(donors)} donors, {donors.pool.nunique()} pools; R2(age ~ pool) = {r2_pool:.3f}; within-pool age SD = {donors.age_c.std():.2f} y "
        f"(total SD {donors.age.std():.2f} y). All B2 metrics are on within-pool age unless labelled 'naive'.")

    # ================================================================== B1. composition
    log("\n" + "=" * 100 + "\n[B1] per-donor cell-type composition (all cells of the raw h5ad) -> within-pool age\n" + "=" * 100)
    import anndata as ad
    a = ad.read_h5ad(PRIMARY["local"], backed="r")
    o = a.obs[["donor_id", "predicted.celltype.l2"]].copy()
    o["donor"] = o.donor_id.astype(str)
    o["l2"] = o["predicted.celltype.l2"].astype(str)
    n0 = len(o)
    o = o[o.donor.isin(donors.index)]
    log(f"[filter] cells of the 981 analysis donors                     cells {n0:,} -> {len(o):,} (removed {n0-len(o):,})")
    n1 = len(o)
    o = o[~o.l2.isin(NON_PBMC)]
    log(f"[filter] drop non-PBMC / artefact labels {sorted(NON_PBMC)}   cells {n1:,} -> {len(o):,} (removed {n1-len(o):,})")
    analysis_types = sorted(set(v for v in CELLTYPE_MAP.values() if v is not None) & set(obs.celltype.unique()))
    o["cat"] = o.l2.map(lambda s: CELLTYPE_MAP.get(s) if CELLTYPE_MAP.get(s) in analysis_types else "other_PBMC")
    other_labels = sorted(set(o.l2[o.cat == "other_PBMC"]))
    log(f"composition categories: {len(analysis_types)} analysis types + 'other_PBMC' = {other_labels} "
        f"({int((o.cat=='other_PBMC').sum()):,} cells, {100*(o.cat=='other_PBMC').mean():.2f}%)")
    counts = o.groupby(["donor", "cat"], observed=True).size().unstack(fill_value=0).reindex(donors.index).fillna(0)
    cats = analysis_types + ["other_PBMC"]
    counts = counts.reindex(columns=cats, fill_value=0)
    Ncell = counts.sum(1)
    prop = counts.div(Ncell, axis=0)
    log(f"cells per donor used for composition: median {Ncell.median():.0f} (min {Ncell.min():.0f}, max {Ncell.max():.0f})")
    K_ = len(cats)
    clr = np.log((counts + 0.5).div(Ncell + 0.5 * K_, axis=0))
    clr = clr.sub(clr.mean(1), axis=0)
    comp_out = pd.concat([donors[["age", "pool", "sex", "age_c"]], Ncell.rename("n_cells"), prop.add_prefix("prop_"), clr.add_prefix("clr_")], axis=1)
    comp_out.to_csv(tab("stageB_donor_composition.csv"))
    log("mean proportion per category (%): " + ", ".join(f"{c} {100*prop[c].mean():.1f}" for c in cats))

    # per-category within-pool association with age
    pool_d = donors.pool.to_numpy()
    prop_c = pool_center(prop.to_numpy(), pool_d)
    age_c = donors.age_c.to_numpy()
    rows = []
    for j, c in enumerate(cats):
        r, p = stats.pearsonr(prop_c[:, j], age_c)
        rs, ps = stats.spearmanr(prop_c[:, j], age_c)
        rows.append(dict(category=c, mean_prop=prop[c].mean(), r_within_pool=r, p=p, spearman_within_pool=rs,
                         r_pooled=stats.pearsonr(prop[c], donors.age)[0]))
    percat = pd.DataFrame(rows).sort_values("r_within_pool")
    log("\nper-category proportion vs age, within pool (Pearson on pool-centred values; 'pooled' = uncorrected for reference):")
    log(percat.to_string(index=False, float_format=lambda x: f"{x:+.3f}"))
    percat.to_csv(tab("stageB_composition_per_category.csv"), index=False)
    # summary axes
    T_cols = [c for c in cats if c.startswith(("CD4_", "CD8_", "Treg", "gdT", "MAIT"))]
    axes_ = pd.DataFrame({
        "naive share of CD4 T": counts.CD4_Naive / counts[[c for c in cats if c.startswith("CD4_") or c == "Treg"]].sum(1),
        "naive share of CD8 T": counts.CD8_Naive / counts[[c for c in cats if c.startswith("CD8_")]].sum(1),
        "CD4_CTL share of CD4 T": counts.CD4_CTL / counts[[c for c in cats if c.startswith("CD4_") or c == "Treg"]].sum(1),
        "naive share of B": counts.B_naive / counts[["B_naive", "B_memory", "B_intermediate"]].sum(1),
        "NK share of PBMC": prop.NK, "T share of PBMC": prop[T_cols].sum(1), "monocyte share of PBMC": prop.CD14_Mono + prop.CD16_Mono,
    }, index=donors.index)
    log("\nsummary axes vs age, within pool:")
    for c in axes_.columns:
        v = pool_center(axes_[c].to_numpy()[:, None], pool_d)[:, 0]
        r, p = stats.pearsonr(v, age_c)
        log(f"   {c:<26} r = {r:+.3f}  (p = {p:.1e});  mean {axes_[c].mean():.3f}")
    axes_.to_csv(tab("stageB_composition_summary_axes.csv"))

    # in-sample partition: age ~ C(pool) + CLR
    Xp = pd.get_dummies(donors.pool.astype("category"), dtype=float).to_numpy()
    Xc = clr.to_numpy()[:, :-1]  # CLR rows sum to 0 -> drop one column for full rank
    y = donors.age.to_numpy(float)
    f0 = sm.OLS(y, Xp).fit()
    f1 = sm.OLS(y, np.hstack([Xp, Xc])).fit()
    inc = f1.rsquared - f0.rsquared
    frac_within = inc / (1 - f0.rsquared)
    ftest = f1.compare_f_test(f0)
    # permutation null for the increment (age shuffled within pool at the donor level)
    null_inc = np.empty(N_PERM_B1)
    Q0, _ = np.linalg.qr(Xp)
    Q1, _ = np.linalg.qr(np.hstack([Xp, Xc]))
    for b in range(N_PERM_B1):
        yp = permute_age_within_pool(obs, rng)  # row-level -> collapse to donor level
        yp_d = pd.Series(yp, index=obs.donor.to_numpy()).groupby(level=0).first().reindex(donors.index).to_numpy()
        ypc = yp_d - yp_d.mean()
        null_inc[b] = (((Q1.T @ ypc) ** 2).sum() - ((Q0.T @ ypc) ** 2).sum()) / (ypc ** 2).sum()
    log(f"\n[B1] in-sample: R2(age ~ pool) = {f0.rsquared:.3f}; R2(age ~ pool + CLR composition[{Xc.shape[1]} df]) = {f1.rsquared:.3f}; "
        f"increment = {inc:.3f} of total age variance = {frac_within:.1%} OF WITHIN-POOL AGE VARIANCE (adj R2 {f0.rsquared_adj:.3f} -> {f1.rsquared_adj:.3f}); "
        f"F p = {ftest[1]:.1e}")
    log(f"     within-pool permutation null for the increment ({N_PERM_B1} draws): mean {null_inc.mean():.4f} (= {null_inc.mean()/(1-f0.rsquared):.1%} of within-pool variance), "
        f"95th pct {np.percentile(null_inc, 95):.4f}, max {null_inc.max():.4f}; observed {inc:.4f} -> p < {1/N_PERM_B1}" if inc > null_inc.max()
        else f"     null mean {null_inc.mean():.4f}; observed {inc:.4f}; p = {(null_inc >= inc).mean():.3f}")
    b1 = dict(r2_pool=f0.rsquared, r2_pool_comp=f1.rsquared, increment=inc, frac_within_pool_insample=frac_within,
              adj_r2_pool=f0.rsquared_adj, adj_r2_pool_comp=f1.rsquared_adj, F_p=float(ftest[1]), null_increment_mean=float(null_inc.mean()),
              null_increment_p95=float(np.percentile(null_inc, 95)), n_comp_df=int(Xc.shape[1]))

    # ================================================================== B2. pool-grouped CV
    log("\n" + "=" * 100 + "\n[B2] pool-grouped CV: age from composition vs from within-cell-type expression (all genes)\n" + "=" * 100)
    ct = obs.celltype.to_numpy(); pool = obs.pool_number.to_numpy(); donor = obs.donor.to_numpy()
    donor_idx = pd.Index(donors.index)
    # --- kernels (features pool-centred within cell type, then standardised; unsupervised) ---
    types = sorted(obs.celltype.unique())
    kern = {}      # type -> (rows idx into obs, K within-pool version, K naive version)
    log("building kernels (all 13,904 genes, standardised within cell type):")
    for t in types:
        rows_t = np.flatnonzero(ct == t)
        X = Y[rows_t]
        Xc_ = pool_center(X, pool[rows_t]); sd = Xc_.std(0); sd[sd == 0] = 1; Xc_ = Xc_ / sd
        Xn = X - X.mean(0); sdn = Xn.std(0); sdn[sdn == 0] = 1; Xn = Xn / sdn
        kern[t] = (rows_t, Xc_ @ Xc_.T / X.shape[1], Xn @ Xn.T / X.shape[1])
        log(f"   {t:<15} n={len(rows_t):>4} pools={len(np.unique(pool[rows_t]))}")
    # all-type donor-level kernel = sum over types (zero for missing donor x type)
    nD = len(donors)
    K_all_w = np.zeros((nD, nD)); K_all_n = np.zeros((nD, nD))
    for t in types:
        rows_t, Kw, Kn = kern[t]
        di = donor_idx.get_indexer(donor[rows_t])
        K_all_w[np.ix_(di, di)] += Kw; K_all_n[np.ix_(di, di)] += Kn
    K_all_w /= len(types); K_all_n /= len(types)
    # composition kernels
    Cw = pool_center(clr.to_numpy(), pool_d); Cw = Cw / np.where(Cw.std(0) == 0, 1, Cw.std(0))
    Cn = clr.to_numpy() - clr.to_numpy().mean(0); Cn = Cn / np.where(Cn.std(0) == 0, 1, Cn.std(0))
    K_comp_w = Cw @ Cw.T / Cw.shape[1]; K_comp_n = Cn @ Cn.T / Cn.shape[1]
    y_d = donors.age.to_numpy(float); yc_d = donors.age_c.to_numpy(float)

    results = []      # long table: repeat, model, metric values
    oof_store = {}    # (model) -> OOF preds of repeat 0 (within version), for bootstrap and for the 2-df combination
    oof_store_all_repeats = {}
    for rep in range(N_REPEATS):
        rrng = np.random.default_rng([RERUN_SEED, rep])
        fold_of_pool = make_pool_folds(donors.pool.unique(), N_OUTER, rrng)
        # (a) composition
        for label, Kw_, Kn_ in (("composition (CLR, 17 categories)", K_comp_w, K_comp_n), ("expression: all 16 cell types jointly", K_all_w, K_all_n)):
            pw, aw = kernel_ridge_oof(Kw_, yc_d, pool_d, donors.index.to_numpy(), fold_of_pool, ALPHAS, rrng)
            pn, an = kernel_ridge_oof(Kn_, y_d, pool_d, donors.index.to_numpy(), fold_of_pool, ALPHAS, rrng)
            r2w, maew = within_pool_r2_mae(yc_d, pw, pool_d); r2n, maen = raw_r2_mae(y_d, pn)
            results.append(dict(repeat=rep, model=label, celltype="(donor level)", n=nD, r2_within_pool=r2w, mae_within_pool=maew,
                                r2_naive=r2n, mae_naive=maen, alpha_within=np.median(aw), alpha_naive=np.median(an)))
            oof_store_all_repeats.setdefault(label, []).append(pw)
            if rep == 0:
                oof_store[label] = (pw, pn)
        # (b) per cell type
        for t in types:
            rows_t, Kw, Kn = kern[t]
            yt = obs.age.to_numpy(float)[rows_t]; ytc = yt - donors.age_pool_mean.reindex(donor[rows_t]).to_numpy()
            pw, aw = kernel_ridge_oof(Kw, ytc, pool[rows_t], donor[rows_t], fold_of_pool, ALPHAS, rrng)
            pn, an = kernel_ridge_oof(Kn, yt, pool[rows_t], donor[rows_t], fold_of_pool, ALPHAS, rrng)
            r2w, maew = within_pool_r2_mae(ytc, pw, pool[rows_t]); r2n, maen = raw_r2_mae(yt, pn)
            results.append(dict(repeat=rep, model="expression: within cell type", celltype=t, n=len(rows_t), r2_within_pool=r2w, mae_within_pool=maew,
                                r2_naive=r2n, mae_naive=maen, alpha_within=np.median(aw), alpha_naive=np.median(an)))
            if rep == 0:
                oof_store[("type", t)] = (pw, pn, rows_t, ytc, yt)
        log(f"   repeat {rep}: comp R2w={results[-1-len(types)-1]['r2_within_pool']:.3f}  all-type expr R2w={results[-1-len(types)]['r2_within_pool']:.3f}  "
            f"median per-type R2w={np.median([r['r2_within_pool'] for r in results[-len(types):]]):.3f}")
    res = pd.DataFrame(results)
    res.to_csv(tab("stageB_cv_results_long.csv"), index=False)

    # summaries with CIs
    def summarise(sub, key):
        return dict(mean=sub[key].mean(), min=sub[key].min(), max=sub[key].max())
    summ_rows = []
    for label in ("composition (CLR, 17 categories)", "expression: all 16 cell types jointly"):
        sub = res[res.model == label]
        pw, pn = oof_store[label]
        lo, hi = bootstrap_ci(yc_d, pw, pool_d, np.random.default_rng([RERUN_SEED, 99]), within=True)
        lon, hin = bootstrap_ci(y_d, pn, pool_d, np.random.default_rng([RERUN_SEED, 98]), within=False)
        summ_rows.append(dict(model=label, celltype="(donor level)", n=nD, r2_within_mean=sub.r2_within_pool.mean(), r2_within_repeat_min=sub.r2_within_pool.min(),
                              r2_within_repeat_max=sub.r2_within_pool.max(), r2_within_boot_lo=lo, r2_within_boot_hi=hi, mae_within=sub.mae_within_pool.mean(),
                              r2_naive_mean=sub.r2_naive.mean(), r2_naive_boot_lo=lon, r2_naive_boot_hi=hin, mae_naive=sub.mae_naive.mean()))
    for t in types:
        sub = res[(res.model == "expression: within cell type") & (res.celltype == t)]
        pw, pn, rows_t, ytc, yt = oof_store[("type", t)]
        lo, hi = bootstrap_ci(ytc, pw, pool[rows_t], np.random.default_rng([RERUN_SEED, 97, types.index(t)]), within=True)
        lon, hin = bootstrap_ci(yt, pn, pool[rows_t], np.random.default_rng([RERUN_SEED, 96, types.index(t)]), within=False)
        summ_rows.append(dict(model="expression: within cell type", celltype=t, n=len(rows_t), r2_within_mean=sub.r2_within_pool.mean(),
                              r2_within_repeat_min=sub.r2_within_pool.min(), r2_within_repeat_max=sub.r2_within_pool.max(), r2_within_boot_lo=lo, r2_within_boot_hi=hi,
                              mae_within=sub.mae_within_pool.mean(), r2_naive_mean=sub.r2_naive.mean(), r2_naive_boot_lo=lon, r2_naive_boot_hi=hin, mae_naive=sub.mae_naive.mean()))
    summ = pd.DataFrame(summ_rows)
    summ.to_csv(tab("stageB_cv_summary.csv"), index=False)
    log("\n[B2] pool-grouped CV summary (R2 mean over 10 repeats; [pool-cluster bootstrap 95% CI, repeat 0]; {repeat min-max}):")
    log(f"{'model':<42}{'n':>5}  {'within-pool R2':<34}{'MAE_w':>7}   {'naive R2 (raw age)':<30}{'MAE_n':>7}")
    for _, r in summ.iterrows():
        name = r.model if r.celltype == "(donor level)" else f"   expression: {r.celltype}"
        log(f"{name:<42}{r.n:>5}  {r.r2_within_mean:+.3f} [{r.r2_within_boot_lo:+.3f}, {r.r2_within_boot_hi:+.3f}] {{{r.r2_within_repeat_min:+.3f},{r.r2_within_repeat_max:+.3f}}}"
            f"{r.mae_within:>7.2f}   {r.r2_naive_mean:+.3f} [{r.r2_naive_boot_lo:+.3f}, {r.r2_naive_boot_hi:+.3f}]      {r.mae_naive:>7.2f}")
    per_t = summ[summ.model == "expression: within cell type"]
    log(f"\nper-cell-type within-pool R2: median {per_t.r2_within_mean.median():+.3f}, range {per_t.r2_within_mean.min():+.3f} .. {per_t.r2_within_mean.max():+.3f} "
        f"(best: {per_t.sort_values('r2_within_mean').celltype.iloc[-1]}); naive median {per_t.r2_naive_mean.median():+.3f}")
    log(f"within-pool age SD = {donors.age_c.std():.2f} y -> a within-pool MAE of {donors.age_c.abs().mean():.2f} y is the 'predict the pool mean' baseline")

    # (c) 2-df combination of OOF predictions (repeat 0) and paired comparison
    pc, _ = oof_store["composition (CLR, 17 categories)"]; pe, _ = oof_store["expression: all 16 cell types jointly"]
    pc_c = pool_center(pc[:, None], pool_d)[:, 0]; pe_c = pool_center(pe[:, None], pool_d)[:, 0]
    def r2_of(Xcols):
        X_ = np.column_stack(Xcols)
        b, *_ = np.linalg.lstsq(X_, yc_d, rcond=None)
        return 1 - ((yc_d - X_ @ b) ** 2).sum() / (yc_d ** 2).sum()
    r2_c, r2_e, r2_ce = r2_of([pc_c]), r2_of([pe_c]), r2_of([pc_c, pe_c])
    r_ce = stats.pearsonr(pc_c, pe_c)[0]
    log(f"\n[B2c] linear recombination of OOF predictions (within-pool; 1-2 in-sample coefficients on 981 donors):")
    log(f"   age_c ~ comp_pred: R2 = {r2_c:.3f};  age_c ~ expr_pred: R2 = {r2_e:.3f};  age_c ~ comp_pred + expr_pred: R2 = {r2_ce:.3f}  "
        f"(gain over composition alone {r2_ce - r2_c:+.3f}; over expression alone {r2_ce - r2_e:+.3f}); corr(comp_pred, expr_pred) = {r_ce:+.3f}")
    # paired repeat-level difference
    d_ = (res[res.model == "expression: all 16 cell types jointly"].sort_values("repeat").r2_within_pool.to_numpy()
          - res[res.model == "composition (CLR, 17 categories)"].sort_values("repeat").r2_within_pool.to_numpy())
    log(f"   paired difference (expression all-type − composition) over {N_REPEATS} repeats: mean {d_.mean():+.3f}, min {d_.min():+.3f}, max {d_.max():+.3f}")

    # sanity shuffle (within-pool, donor-level) — one draw per model, repeat-0 folds
    log("\n[B2-shuffle] within-pool donor-level age shuffle (one draw; pipeline sanity, not the Stage D control):")
    srng = np.random.default_rng([RERUN_SEED, 7])
    yp = permute_age_within_pool(obs, srng)
    yp_d = pd.Series(yp, index=obs.donor.to_numpy()).groupby(level=0).first().reindex(donors.index).to_numpy()
    ypc_d = yp_d - donors.age_pool_mean.to_numpy()
    fold_of_pool0 = make_pool_folds(donors.pool.unique(), N_OUTER, np.random.default_rng([RERUN_SEED, 0]))
    shuf = {}
    for label, Kw_ in (("composition", K_comp_w), ("expression all-type", K_all_w)):
        pw, _ = kernel_ridge_oof(Kw_, ypc_d, pool_d, donors.index.to_numpy(), fold_of_pool0, ALPHAS, srng)
        shuf[label] = within_pool_r2_mae(ypc_d, pw, pool_d)[0]
    per_t_shuf = []
    for t in types:
        rows_t, Kw, Kn = kern[t]
        ytc_p = pd.Series(yp, index=obs.donor.to_numpy()).groupby(level=0).first().reindex(donor[rows_t]).to_numpy() - donors.age_pool_mean.reindex(donor[rows_t]).to_numpy()
        pw, _ = kernel_ridge_oof(Kw, ytc_p, pool[rows_t], donor[rows_t], fold_of_pool0, ALPHAS, srng)
        per_t_shuf.append(within_pool_r2_mae(ytc_p, pw, pool[rows_t])[0])
    log(f"   shuffled within-pool R2: composition {shuf['composition']:+.3f}; expression all-type {shuf['expression all-type']:+.3f}; "
        f"per-type median {np.median(per_t_shuf):+.3f} (range {min(per_t_shuf):+.3f} .. {max(per_t_shuf):+.3f})")

    # ================================================================== B3. where the per-gene age signal sits
    log("\n" + "=" * 100 + "\n[B3] share of the summed within-pool per-gene age variance by Stage A class (classes UNCHANGED)\n" + "=" * 100)
    d = pd.read_csv(tab("stageA_variance_decomposition.csv"), index_col=0)
    assert (d.index == genes.gene_id).all()
    var_tot = Y.var(0, ddof=1)
    d["abs_age_var"] = d.unique_age * var_tot           # absolute log2-CPM variance explained by within-pool age
    d["abs_age_var_pooled"] = d.unique_age_pooled * var_tot
    d["id_bin"] = pd.cut(d.unique_ct, [-1, 0.2, 0.5, 1.01], labels=["unique_ct<=0.2", "0.2<unique_ct<0.5", "unique_ct>=0.5"])
    rows = []
    for name, mask in (("all genes", np.ones(len(d), bool)), ("permutation q<0.05", (d.age_q_BH_perm < 0.05).to_numpy()),
                       ("unique_age>=0.002", (d.unique_age >= 0.002).to_numpy())):
        sub = d[mask]
        for grp_col in ("gene_class_rerun", "id_bin"):
            g = sub.groupby(grp_col, observed=False).agg(n=("unique_age", "size"), sum_unique_age=("unique_age", "sum"), sum_abs_var=("abs_age_var", "sum"),
                                                          sum_abs_var_pooled=("abs_age_var_pooled", "sum"))
            g["share_unique_age"] = g.sum_unique_age / sub.unique_age.sum()
            g["share_abs_var"] = g.sum_abs_var / sub.abs_age_var.sum()
            g["share_abs_var_pooled_model"] = g.sum_abs_var_pooled / sub.abs_age_var_pooled.sum()
            g["gene_set"] = name; g["grouping"] = grp_col
            rows.append(g.reset_index().rename(columns={grp_col: "group"}))
    b3 = pd.concat(rows, ignore_index=True)
    b3.to_csv(tab("stageB_age_signal_by_class.csv"), index=False)
    for name in ("all genes", "permutation q<0.05", "unique_age>=0.002"):
        log(f"\n{name}:")
        for grp_col in ("gene_class_rerun", "id_bin"):
            s = b3[(b3.gene_set == name) & (b3.grouping == grp_col)]
            log(f"   by {grp_col}:")
            for _, r in s.iterrows():
                log(f"      {str(r.group):<20} n={int(r.n):>6}  share of Σ unique_age = {r.share_unique_age:6.1%}   share of Σ absolute age variance = {r.share_abs_var:6.1%}"
                    f"   (prior pooled model, abs: {r.share_abs_var_pooled_model:6.1%})")
    A_share = float(b3[(b3.gene_set == "all genes") & (b3.grouping == "gene_class_rerun") & (b3.group == "A")].share_unique_age.iloc[0])
    MI_share = float(b3[(b3.gene_set == "all genes") & (b3.grouping == "gene_class_rerun") & (b3.group.isin(["MIXED", "I"]))].share_unique_age.sum())
    hi_share = float(b3[(b3.gene_set == "all genes") & (b3.grouping == "id_bin") & (b3.group == "unique_ct>=0.5")].share_unique_age.iloc[0])
    log(f"\n[B3] clean A-genes (2 genes) hold {A_share:.1%} of the summed within-pool per-gene age variance; MIXED+I genes hold {MI_share:.1%}; "
        f"genes with unique_ct>=0.5 hold {hi_share:.1%} (they are {(d.unique_ct>=0.5).mean():.1%} of genes).")
    log("   Caveat: summed per-gene unique_age is an accounting of marginal effects, not multivariate predictability (that is B2 / Stage C).")

    # ================================================================== summary json + figures
    comp_r = summ[summ.model.str.startswith("composition")].iloc[0]; expr_r = summ[summ.model.str.startswith("expression: all")].iloc[0]
    summary = dict(seed=RERUN_SEED, n_donors=nD, r2_age_pool=r2_pool, B1=b1,
                   B2=dict(composition_r2_within=float(comp_r.r2_within_mean), composition_ci=[float(comp_r.r2_within_boot_lo), float(comp_r.r2_within_boot_hi)],
                           composition_mae_within=float(comp_r.mae_within), composition_r2_naive=float(comp_r.r2_naive_mean),
                           expression_alltype_r2_within=float(expr_r.r2_within_mean), expression_alltype_ci=[float(expr_r.r2_within_boot_lo), float(expr_r.r2_within_boot_hi)],
                           expression_alltype_mae_within=float(expr_r.mae_within), expression_alltype_r2_naive=float(expr_r.r2_naive_mean),
                           expression_per_type_r2_within_median=float(per_t.r2_within_mean.median()), expression_per_type_r2_within_max=float(per_t.r2_within_mean.max()),
                           expression_per_type_best=str(per_t.sort_values('r2_within_mean').celltype.iloc[-1]),
                           combined_2df_r2=float(r2_ce), corr_comp_expr_pred=float(r_ce), shuffle=dict(composition=float(shuf["composition"]),
                           expression_alltype=float(shuf["expression all-type"]), per_type_median=float(np.median(per_t_shuf)))),
                   B3=dict(A_share_sum_unique_age=A_share, MIXED_plus_I_share=MI_share, unique_ct_ge_0p5_share=hi_share))
    json.dump(summary, open(RERUN_DIR / "stageB_summary.json", "w"), indent=1)
    make_figures(percat, axes_, donors, prop, summ, res, b3, pc_c, pe_c, yc_d)
    log("\nStage B complete (measurement only; no stop condition).")
    log.close()
    return summary


def make_figures(percat, axes_, donors, prop, summ, res, b3, pc_c, pe_c, yc_d):
    pool_d = donors.pool.to_numpy(); age_c = donors.age_c.to_numpy()
    # B1
    figu, axs = plt.subplots(1, 2, figsize=(13, 4.8))
    ax = axs[0]
    cols = ["tab:red" if p < 0.001 else ("salmon" if p < 0.05 else "lightgrey") for p in percat.p]
    ax.barh(range(len(percat)), percat.r_within_pool, color=cols)
    ax.set_yticks(range(len(percat))); ax.set_yticklabels(percat.category, fontsize=8); ax.axvline(0, c="k", lw=0.6)
    ax.set_xlabel("within-pool Pearson r (proportion, age)"); ax.set_title("B1: cell-type proportion vs age within pool (red p<0.001, salmon p<0.05)")
    ax = axs[1]
    v = pool_center(axes_["naive share of CD4 T"].to_numpy()[:, None], pool_d)[:, 0]
    ax.scatter(age_c, v, s=6, alpha=0.5, c="tab:blue")
    b = np.polyfit(age_c, v, 1); xs = np.linspace(age_c.min(), age_c.max(), 10); ax.plot(xs, np.polyval(b, xs), "k-", lw=1)
    ax.set_xlabel("donor age − pool mean age (years)"); ax.set_ylabel("naive share of CD4 T − pool mean")
    ax.set_title(f"naive CD4 share vs within-pool age: r = {stats.pearsonr(v, age_c)[0]:+.3f}")
    figu.tight_layout(); figu.savefig(fig("stageB_composition_vs_age.png"), dpi=140); plt.close(figu)
    # B2
    figu, axs = plt.subplots(1, 2, figsize=(14, 5), gridspec_kw=dict(width_ratios=[1, 1.6]))
    ax = axs[0]
    top = summ[summ.celltype == "(donor level)"]
    per_t = summ[summ.model == "expression: within cell type"]
    labels = ["composition\n(17 categories)", "expression\nall 16 types", f"expression\nper type (median)"]
    vals = [top.iloc[0].r2_within_mean, top.iloc[1].r2_within_mean, per_t.r2_within_mean.median()]
    err = np.clip([[top.iloc[0].r2_within_mean - top.iloc[0].r2_within_boot_lo, top.iloc[1].r2_within_mean - top.iloc[1].r2_within_boot_lo, per_t.r2_within_mean.median() - per_t.r2_within_mean.min()],
                   [top.iloc[0].r2_within_boot_hi - top.iloc[0].r2_within_mean, top.iloc[1].r2_within_boot_hi - top.iloc[1].r2_within_mean, per_t.r2_within_mean.max() - per_t.r2_within_mean.median()]], 0, None)
    ax.bar(range(3), vals, yerr=err, capsize=4, color=["tab:orange", "tab:blue", "lightsteelblue"])
    naive = [top.iloc[0].r2_naive_mean, top.iloc[1].r2_naive_mean, per_t.r2_naive_mean.median()]
    ax.scatter(range(3), naive, marker="_", s=300, c="k", label="naive pool-grouped R² (raw age)")
    ax.set_xticks(range(3)); ax.set_xticklabels(labels, fontsize=8); ax.set_ylabel("held-out R² (within-pool age)"); ax.axhline(0, c="k", lw=0.5)
    ax.set_title("B2: age from composition vs within-type expression\n(bars: within-pool R², 95% pool-cluster bootstrap CI; per-type bar: range over types)", fontsize=9); ax.legend(fontsize=7)
    ax = axs[1]
    pt = per_t.sort_values("r2_within_mean")
    ax.barh(range(len(pt)), pt.r2_within_mean, xerr=np.clip([pt.r2_within_mean - pt.r2_within_boot_lo, pt.r2_within_boot_hi - pt.r2_within_mean], 0, None), color="lightsteelblue", capsize=2)
    ax.scatter(pt.r2_naive_mean, range(len(pt)), marker="|", s=120, c="k", label="naive R²")
    ax.axvline(top.iloc[0].r2_within_mean, color="tab:orange", ls="--", label="composition alone (within-pool)")
    ax.set_yticks(range(len(pt))); ax.set_yticklabels([f"{c} (n={n})" for c, n in zip(pt.celltype, pt.n)], fontsize=8); ax.axvline(0, c="k", lw=0.5)
    ax.set_xlabel("held-out within-pool R²"); ax.set_title("expression, all genes, per cell type (pool-grouped CV)"); ax.legend(fontsize=7)
    figu.tight_layout(); figu.savefig(fig("stageB_cv_composition_vs_expression.png"), dpi=140); plt.close(figu)
    # B2c scatter
    figu, axs = plt.subplots(1, 3, figsize=(14, 4.4))
    for ax, (x, lab_) in zip(axs[:2], ((pc_c, "composition OOF prediction"), (pe_c, "expression (all types) OOF prediction"))):
        ax.scatter(x, yc_d, s=5, alpha=0.4); ax.set_xlabel(lab_ + " (pool-centred)"); ax.set_ylabel("age − pool mean")
        ax.set_title(f"r = {stats.pearsonr(x, yc_d)[0]:+.3f}")
    axs[2].scatter(pc_c, pe_c, s=5, alpha=0.4, c="tab:green"); axs[2].set_xlabel("composition prediction"); axs[2].set_ylabel("expression prediction")
    axs[2].set_title(f"agreement of the two predictors: r = {stats.pearsonr(pc_c, pe_c)[0]:+.3f}")
    figu.tight_layout(); figu.savefig(fig("stageB_oof_predictions.png"), dpi=140); plt.close(figu)
    # B3
    figu, axs = plt.subplots(1, 2, figsize=(12, 4.4))
    for ax, grp_col, title in zip(axs, ("gene_class_rerun", "id_bin"), ("by Stage A primary class", "by identity-variance bin")):
        sets = ["all genes", "permutation q<0.05", "unique_age>=0.002"]
        groups = list(b3[(b3.grouping == grp_col) & (b3.gene_set == "all genes")].group)
        bottom = np.zeros(len(sets))
        for g in groups:
            vals = [float(b3[(b3.grouping == grp_col) & (b3.gene_set == s) & (b3.group == g)].share_unique_age.sum()) for s in sets]
            ax.bar(range(len(sets)), vals, bottom=bottom, label=str(g)); bottom += vals
        ax.set_xticks(range(len(sets))); ax.set_xticklabels(sets, fontsize=8); ax.set_ylabel("share of Σ within-pool unique_age"); ax.set_title(f"B3: {title}"); ax.legend(fontsize=7)
    figu.tight_layout(); figu.savefig(fig("stageB_age_signal_by_class.png"), dpi=140); plt.close(figu)


if __name__ == "__main__":
    main()
