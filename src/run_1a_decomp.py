"""Phase 1a runner: normalize pseudobulk, decompose variance per gene, report distributions, sweep set sizes, plot."""
import sys
import json
import numpy as np
import pandas as pd
import anndata as ad
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

sys.path.insert(0, str(__import__("pathlib").Path(__file__).resolve().parent))
from config import DATA_PROC, TAB, FIG, SEED  # noqa
from decomp import normalize_logcpm, decompose, classify, SWEEP, PRIMARY_THRESH  # noqa

np.random.seed(SEED)


def main():
    pb = ad.read_h5ad(DATA_PROC / "onek1k_pseudobulk.h5ad")
    print(f"pseudobulk: {pb.shape[0]:,} samples x {pb.shape[1]:,} genes; donors={pb.obs.donor.nunique()}; cell types={pb.obs.celltype.nunique()}")
    print(pb.obs.groupby("celltype", observed=True).agg(n_pseudobulks=("donor", "size"), n_donors=("donor", "nunique"),
                                                        median_cells=("n_cells", "median")).to_string())
    Y, gid, sym = normalize_logcpm(pb)
    np.save(DATA_PROC / "onek1k_logcpm.npy", Y)
    pd.DataFrame({"gene_id": gid, "symbol": sym}).to_csv(DATA_PROC / "onek1k_logcpm_genes.csv", index=False)
    pb.obs.to_csv(DATA_PROC / "onek1k_pseudobulk_obs.csv")

    df = decompose(Y, pb.obs.celltype.to_numpy(), pb.obs.age.to_numpy(), gene_ids=gid, symbols=sym)
    df.to_csv(TAB / "phase1a_variance_decomposition.csv")
    print("\n=== variance component distribution over genes (fraction of total variance) ===")
    q = [0.5, 0.9, 0.95, 0.99, 0.999, 1.0]
    print(df[["unique_ct", "unique_age", "shared", "resid", "interaction", "partial_age"]].quantile(q).T.to_string(float_format=lambda x: f"{x:.4f}"))
    print(f"\nshared component: min={df.shared.min():.4f} max={df.shared.max():.4f} (near 0 => age and cell type ~orthogonal in design)")
    print("\nTop 25 genes by unique age variance:")
    print(df.sort_values("unique_age", ascending=False).head(25)[["symbol", "unique_age", "unique_ct", "interaction", "age_slope_per_year", "mean_logcpm"]].to_string(float_format=lambda x: f"{x:.4f}"))
    print("\nTop 15 genes by unique cell-type variance:")
    print(df.sort_values("unique_ct", ascending=False).head(15)[["symbol", "unique_age", "unique_ct", "interaction", "mean_logcpm"]].to_string(float_format=lambda x: f"{x:.4f}"))

    # ---- sweep set sizes ----
    rows = []
    for t_age in SWEEP["t_age"]:
        for k_ct in SWEEP["k_ct"]:
            for t_ct in SWEEP["t_ct"]:
                for k_age in SWEEP["k_age"]:
                    lab = classify(df, t_age, k_ct, t_ct, k_age)
                    vc = lab.value_counts()
                    rows.append(dict(t_age=t_age, k_ct=k_ct, t_ct=t_ct, k_age=k_age, n_A=int(vc.get("A", 0)),
                                     n_I=int(vc.get("I", 0)), n_MIXED=int(vc.get("MIXED", 0)), n_OTHER=int(vc.get("OTHER", 0))))
    sweep = pd.DataFrame(rows)
    sweep.to_csv(TAB / "phase1a_threshold_sweep_set_sizes.csv", index=False)
    P = PRIMARY_THRESH
    print(f"\n=== set sizes: A-genes by (t_age, k_ct) [I thresholds fixed t_ct={P['t_ct']},k_age={P['k_age']}] ===")
    print(sweep[(sweep.t_ct == P["t_ct"]) & (sweep.k_age == P["k_age"])].pivot(index="t_age", columns="k_ct", values="n_A").to_string())
    print(f"\n=== set sizes: MIXED by (t_age, t_ct) [k_ct={P['k_ct']},k_age={P['k_age']}] ===")
    print(sweep[(sweep.k_ct == P["k_ct"]) & (sweep.k_age == P["k_age"])].pivot(index="t_age", columns="t_ct", values="n_MIXED").to_string())
    print(f"\n=== set sizes: I-genes by (t_ct, k_age) [t_age={P['t_age']},k_ct={P['k_ct']}] ===")
    print(sweep[(sweep.t_age == P["t_age"]) & (sweep.k_ct == P["k_ct"])].pivot(index="t_ct", columns="k_age", values="n_I").to_string())

    lab = classify(df, **P)
    df["gene_class_primary"] = lab
    df.to_csv(TAB / "phase1a_variance_decomposition.csv")
    print(f"\n=== PRIMARY thresholds {P} ===")
    print(lab.value_counts().to_string())
    print("A-genes (primary):", ", ".join(df.loc[lab == "A"].sort_values("unique_age", ascending=False).symbol.tolist()))
    print("MIXED (primary):", ", ".join(df.loc[lab == "MIXED"].sort_values("unique_age", ascending=False).symbol.tolist()))

    # ---- sensitivity: only well-sampled pseudobulks (>=100 cells) to gauge sampling-noise deflation ----
    big = (pb.obs.n_cells >= 100).to_numpy()
    cts_big = pb.obs.celltype[big].value_counts()
    cts_big = cts_big[cts_big >= 100].index
    sel = big & pb.obs.celltype.isin(cts_big).to_numpy()
    df_big = decompose(Y[sel], pb.obs.celltype.to_numpy()[sel], pb.obs.age.to_numpy()[sel], gene_ids=gid, symbols=sym)
    df_big.to_csv(TAB / "phase1a_variance_decomposition_ge100cells.csv")
    print(f"\n=== sensitivity: pseudobulks with >=100 cells only (n={sel.sum():,}; cell types={list(cts_big)}) ===")
    print(df_big[["unique_ct", "unique_age", "resid", "partial_age"]].quantile(q).T.to_string(float_format=lambda x: f"{x:.4f}"))
    print("Spearman(unique_age all vs >=100-cell):", f"{df.unique_age.corr(df_big.unique_age, method='spearman'):.3f}")

    # ---- batch (pool) vs age confound check at the donor level ----
    if "pool_number" in pb.obs:
        don = pb.obs.groupby("donor", observed=True).agg(age=("age", "first"), pool=("pool_number", "first"),
                                                          n_pools=("pool_number", "nunique"))
        assert (don.n_pools == 1).all(), "a donor appears in >1 pool"
        import statsmodels.formula.api as smf
        fit = smf.ols("age ~ C(pool)", data=don.assign(pool=don.pool.astype(str))).fit()
        print(f"\nBatch check: age ~ C(pool) across {len(don)} donors, {don.pool.nunique()} pools: R2={fit.rsquared:.3f}, "
              f"adj R2={fit.rsquared_adj:.3f}, F p={fit.f_pvalue:.3g}")

    # ---- joint distribution plot ----
    fig, axes = plt.subplots(1, 2, figsize=(12, 5))
    ax = axes[0]
    ax.scatter(df.unique_ct, df.unique_age, s=3, alpha=0.3, c="grey")
    ax.set_xlabel("unique cell-type variance (fraction of total)"); ax.set_ylabel("unique age variance (fraction of total)")
    ax.set_title(f"OneK1K pseudobulk: {len(df):,} genes")
    ax = axes[1]
    ax.scatter(df.unique_ct, df.unique_age, s=3, alpha=0.3, c="grey")
    ax.set_yscale("log"); ax.set_ylim(1e-5, 1); ax.set_xlabel("unique cell-type variance"); ax.set_ylabel("unique age variance (log)")
    for t in SWEEP["t_age"]:
        ax.axhline(t, color="tab:red", lw=0.5, ls="--")
    for k in SWEEP["k_ct"]:
        ax.axvline(k, color="tab:blue", lw=0.5, ls=":")
    ax.set_title("sweep grid: red = t_age, blue = k_ct")
    fig.tight_layout(); fig.savefig(FIG / "phase1a_joint_variance_distribution.png", dpi=140)
    # stacked bar of mean components
    comp = df[["unique_ct", "unique_age", "shared", "interaction", "resid"]].clip(lower=0)
    comp["resid"] = comp["resid"] - comp["interaction"]
    fig, ax = plt.subplots(figsize=(6, 3.5))
    means = comp.mean()
    ax.bar(means.index, means.values, color=["tab:blue", "tab:red", "tab:purple", "tab:orange", "lightgrey"])
    ax.set_ylabel("mean fraction of variance (all genes)"); ax.set_title("average variance partition per gene")
    fig.tight_layout(); fig.savefig(FIG / "phase1a_mean_variance_partition.png", dpi=140)
    print("\nMean partition over genes:", means.round(4).to_dict())
    return df


if __name__ == "__main__":
    main()
