"""Verify the 1b batch finding: (a) simulation of a true age effect -> expected within-pool retention;
(b) depth ~ pool confounding; (c) genome-wide within-pool age variance and how many genes pass the
primary A-gene threshold on within-pool age variance; (d) direct within-pool vs pooled age correlations."""
import sys
import numpy as np
import pandas as pd
import statsmodels.formula.api as smf
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

sys.path.insert(0, str(__import__("pathlib").Path(__file__).resolve().parent))
from config import DATA_PROC, TAB, FIG, SEED  # noqa
from decomp import decompose, PRIMARY_THRESH  # noqa
from confounds import within_pool_age_r2, residualize_on_celltype, corr_with_covariate  # noqa

rng = np.random.default_rng(SEED)


def main():
    Y = np.load(DATA_PROC / "onek1k_logcpm.npy")
    genes = pd.read_csv(DATA_PROC / "onek1k_logcpm_genes.csv")
    obs = pd.read_csv(DATA_PROC / "onek1k_pseudobulk_obs.csv", index_col=0)
    F = pd.read_csv(TAB / "phase1b_gene_flags.csv", index_col=0)
    ct, pool, age = obs.celltype.to_numpy(), obs.pool_number.to_numpy(), obs.age.to_numpy(float)

    # (a) simulation: true additive age effect + pool batch effect (uncorrelated with age) + noise
    print("=== (a) simulation: expected within-pool retention for a REAL age effect ===")
    pools = pd.Categorical(pool)
    pool_eff = rng.normal(0, 0.3, len(pools.categories))[pools.codes]
    ct_eff = rng.normal(0, 1.0, len(pd.Categorical(ct).categories))[pd.Categorical(ct).codes]
    rows = []
    for beta in (0.005, 0.01, 0.02):
        for rep in range(20):
            y = beta * age + ct_eff + pool_eff + rng.normal(0, 0.5, len(age))
            d = decompose(y[:, None], ct, age)
            wp = within_pool_age_r2(y[:, None], ct, pool, age)[0]
            rows.append(dict(beta=beta, unique_age=d.unique_age[0], within_pool=wp, retained=wp / d.unique_age[0]))
    sim = pd.DataFrame(rows).groupby("beta").agg(unique_age=("unique_age", "mean"), within_pool=("within_pool", "mean"),
                                                retained_mean=("retained", "mean"), retained_min=("retained", "min"))
    print(sim.to_string(float_format=lambda x: f"{x:.4f}"))
    don = obs.groupby("donor").agg(age=("age", "first"), pool=("pool_number", "first"))
    r2_pool = smf.ols("age ~ C(pool)", data=don.assign(pool=don.pool.astype(str))).fit().rsquared
    print(f"expected retention ~ 1 - R2(age~pool) = {1-r2_pool:.3f}   (observed median for primary A-genes: "
          f"{F.loc[F.gene_class_primary=='A','age_signal_retained_within_pool'].median():.3f})")

    # (b) depth ~ pool
    print("\n=== (b) is per-cell depth a pool-level (batch) property? ===")
    for c in ("mean_counts_per_cell", "mean_genes_per_cell"):
        v = residualize_on_celltype(np.log(obs[c].to_numpy(float))[:, None], ct)[:, 0]
        d = pd.DataFrame(dict(v=v, pool=pool.astype(str)))
        r2 = smf.ols("v ~ C(pool)", data=d).fit().rsquared
        print(f"   R2( log {c} [cell-type-residualized] ~ C(pool) ) = {r2:.3f}")
    # age vs depth at donor level, within vs between pool
    dd = obs.groupby("donor").agg(age=("age", "first"), pool=("pool_number", "first"),
                                  depth=("mean_counts_per_cell", "mean"))
    dd["depth_pool_mean"] = dd.groupby("pool").depth.transform("mean")
    dd["age_pool_mean"] = dd.groupby("pool").age.transform("mean")
    r_between = np.corrcoef(dd.groupby("pool").age.mean(), dd.groupby("pool").depth.mean())[0, 1]
    r_within = np.corrcoef(dd.age - dd.age_pool_mean, dd.depth - dd.depth_pool_mean)[0, 1]
    print(f"   donor-level corr(age, depth): between-pool means r={r_between:+.3f} (n={dd.pool.nunique()} pools); within-pool r={r_within:+.3f}")

    # (c) genome-wide within-pool age variance
    print("\n=== (c) genome-wide: within-pool unique age variance ===")
    F["unique_age_within_pool"] = within_pool_age_r2(Y, ct, pool, age)
    q = [0.5, 0.9, 0.95, 0.99, 0.999, 1.0]
    print(pd.DataFrame({"pooled unique_age": F.unique_age.quantile(q), "within-pool unique_age": F.unique_age_within_pool.quantile(q)}).T.to_string(float_format=lambda x: f"{x:.5f}"))
    for t in (0.002, 0.005, 0.01):
        n_pooled = int((F.unique_age >= t).sum()); n_wp = int((F.unique_age_within_pool >= t).sum())
        n_both = int(((F.unique_age >= t) & (F.unique_age_within_pool >= t)).sum())
        print(f"   genes with unique_age >= {t}: pooled={n_pooled}, within-pool={n_wp}, both={n_both}")
    print("Spearman(pooled unique_age, within-pool unique_age) over all genes:", f"{F.unique_age.corr(F.unique_age_within_pool, method='spearman'):.3f}")
    print("\nTop 25 genes by WITHIN-POOL unique age variance:")
    cols = ["symbol", "unique_age_within_pool", "unique_age", "unique_ct", "r_mean_counts", "flag_cell_cycle", "flag_technical", "flag_sex"]
    print(F.sort_values("unique_age_within_pool", ascending=False).head(25)[cols].to_string(float_format=lambda x: f"{x:.4f}"))

    # (d) direct correlations for primary A-genes: pooled vs within-pool
    print("\n=== (d) primary A-genes: age correlation pooled vs within-pool (cell-type-residualized expression) ===")
    A = F.index[F.gene_class_primary == "A"]
    Yr = residualize_on_celltype(Y, ct)
    age_r = residualize_on_celltype(age[:, None], ct)[:, 0]
    D = pd.get_dummies(pd.Categorical(ct), dtype=float).to_numpy()
    Dp = pd.get_dummies(pd.Categorical(pool), dtype=float).to_numpy()
    Q, _ = np.linalg.qr(np.hstack([D, Dp[:, 1:]]))
    Yrp = Y - Q @ (Q.T @ Y); age_rp = (age[:, None] - Q @ (Q.T @ age[:, None]))[:, 0]
    r_pooled = corr_with_covariate(Yr, age_r); r_wp = corr_with_covariate(Yrp, age_rp)
    F["r_age_pooled"] = r_pooled; F["r_age_within_pool"] = r_wp
    print(F.loc[A].sort_values("unique_age", ascending=False).head(20)[["symbol", "r_age_pooled", "r_age_within_pool", "r_mean_counts"]].to_string(float_format=lambda x: f"{x:+.3f}"))
    print(f"\n   primary A-genes: median |r_age| pooled={F.loc[A,'r_age_pooled'].abs().median():.3f}, within-pool={F.loc[A,'r_age_within_pool'].abs().median():.3f}")
    print(f"   all genes: 99th pct |r_age| pooled={F.r_age_pooled.abs().quantile(.99):.3f}, within-pool={F.r_age_within_pool.abs().quantile(.99):.3f}")
    # null for within-pool r: shuffle age within pool
    age_perm = age.copy()
    for p in np.unique(pool):
        idx = np.flatnonzero(pool == p)
        # shuffle donors' ages within pool (keep donor structure: permute at donor level)
        dsub = obs.iloc[idx]
        dmap = dict(zip(dsub.donor.unique(), rng.permutation(dsub.groupby("donor").age.first().values)))
        age_perm[idx] = dsub.donor.map(dmap).values
    age_perm_rp = (age_perm[:, None] - Q @ (Q.T @ age_perm[:, None]))[:, 0]
    r_null = corr_with_covariate(Yrp, age_perm_rp)
    print(f"   null (age permuted within pool at donor level): 99th pct |r| = {np.nanquantile(np.abs(r_null), .99):.3f}, max={np.nanmax(np.abs(r_null)):.3f}")
    F.to_csv(TAB / "phase1b_gene_flags.csv")

    fig, axes = plt.subplots(1, 2, figsize=(12, 5))
    ax = axes[0]
    ax.scatter(F.unique_age, F.unique_age_within_pool, s=3, alpha=0.3, c="grey", label="all genes")
    ax.scatter(F.loc[A].unique_age, F.loc[A].unique_age_within_pool, s=10, c="tab:red", label="primary A-genes")
    lim = max(F.unique_age.max(), F.unique_age_within_pool.max()) * 1.05
    ax.plot([0, lim], [0, lim], "k--", lw=0.7, label="y = x"); ax.plot([0, lim], [0, (1 - r2_pool) * lim], "b:", lw=0.9, label=f"y = {1-r2_pool:.2f}x (real-effect expectation)")
    ax.set_xlabel("unique age variance, pooled model (ct + age)"); ax.set_ylabel("unique age variance, within-pool (ct + pool + age)")
    ax.legend(fontsize=8); ax.set_title("age signal before vs after batch (pool) adjustment")
    ax = axes[1]
    ax.scatter(F.r_age_pooled, F.r_age_within_pool, s=3, alpha=0.3, c="grey")
    ax.scatter(F.loc[A].r_age_pooled, F.loc[A].r_age_within_pool, s=10, c="tab:red")
    ax.axhline(0, c="k", lw=0.5); ax.axvline(0, c="k", lw=0.5)
    ax.set_xlabel("r(expression, age) pooled"); ax.set_ylabel("r(expression, age) within pool")
    ax.set_title("per-gene age correlation: pooled vs within-pool")
    fig.tight_layout(); fig.savefig(FIG / "phase1b_batch_adjustment_age_signal.png", dpi=140)


if __name__ == "__main__":
    main()
