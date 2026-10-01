"""Is sequencing depth the mediator? Compare age-variance retention after adjusting for
(i) pool (74 df) vs (ii) only log mean UMI/cell + log mean genes/cell (2 df) vs (iii) both."""
import sys
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

sys.path.insert(0, str(__import__("pathlib").Path(__file__).resolve().parent))
from config import DATA_PROC, TAB, FIG  # noqa


def unique_age_given(Y, ct, covs, age):
    Yc = Y - Y.mean(0, keepdims=True); ss = (Yc ** 2).sum(0)
    D = pd.get_dummies(pd.Categorical(ct), dtype=float).to_numpy()
    X0 = np.hstack([D] + covs) if covs else D
    X1 = np.hstack([X0, (age - age.mean())[:, None]])
    Q0, _ = np.linalg.qr(X0); Q1, _ = np.linalg.qr(X1)
    return (((Q1.T @ Yc) ** 2).sum(0) - ((Q0.T @ Yc) ** 2).sum(0)) / ss


def main():
    Y = np.load(DATA_PROC / "onek1k_logcpm.npy").astype(np.float64)
    obs = pd.read_csv(DATA_PROC / "onek1k_pseudobulk_obs.csv", index_col=0)
    F = pd.read_csv(TAB / "phase1b_gene_flags.csv", index_col=0)
    ct, age = obs.celltype.to_numpy(), obs.age.to_numpy(float)
    pool = pd.get_dummies(pd.Categorical(obs.pool_number), dtype=float).to_numpy()[:, 1:]
    depth = np.column_stack([np.log(obs.mean_counts_per_cell), np.log(obs.mean_genes_per_cell)])
    depth = depth - depth.mean(0)
    base = unique_age_given(Y, ct, [], age)
    F["ua_given_depth"] = unique_age_given(Y, ct, [depth], age)
    F["ua_given_pool"] = unique_age_given(Y, ct, [pool], age)
    F["ua_given_pool_depth"] = unique_age_given(Y, ct, [pool, depth], age)
    A = F.index[F.gene_class_primary == "A"]
    print("Primary A-genes (n=%d): fraction of pooled age variance retained after adjusting for..." % len(A))
    for c, lab in (("ua_given_depth", "depth covariates only (2 df)"), ("ua_given_pool", "pool only (74 df)"), ("ua_given_pool_depth", "pool + depth")):
        r = (F.loc[A, c] / base[F.index.get_indexer(A)])
        print(f"   {lab:<32} median={r.median():.3f}  mean={r.mean():.3f}  n(<0.25)={int((r<0.25).sum())}  n(<0.5)={int((r<0.5).sum())}")
    print("\nAll genes with pooled unique_age >= 0.002 (n=%d):" % int((base >= 0.002).sum()))
    sel = base >= 0.002
    for c, lab in (("ua_given_depth", "depth only"), ("ua_given_pool", "pool only"), ("ua_given_pool_depth", "pool + depth")):
        r = F[c].to_numpy()[sel] / base[sel]
        print(f"   {lab:<14} median retention={np.median(r):.3f}; genes still >=0.002: {int((F[c].to_numpy()[sel] >= 0.002).sum())}")
    F.to_csv(TAB / "phase1b_gene_flags.csv")

    # figure: pool mean age vs pool mean depth
    dd = obs.groupby("donor").agg(age=("age", "first"), pool=("pool_number", "first"), depth=("mean_counts_per_cell", "mean"))
    pm = dd.groupby("pool").agg(age=("age", "mean"), depth=("depth", "mean"), n=("age", "size"))
    fig, axes = plt.subplots(1, 2, figsize=(11, 4.2))
    axes[0].scatter(pm.age, pm.depth, s=pm.n * 3, alpha=0.7)
    axes[0].set_xlabel("pool mean donor age (years)"); axes[0].set_ylabel("pool mean UMI per cell")
    axes[0].set_title(f"75 pools: r = {np.corrcoef(pm.age, pm.depth)[0,1]:+.2f}  (size = n donors)")
    order = pm.sort_values("age").index
    axes[1].boxplot([dd.age[dd.pool == p] for p in order], positions=range(len(order)), widths=0.6, showfliers=False)
    axes[1].set_xticks([]); axes[1].set_xlabel("pool (sorted by mean age)"); axes[1].set_ylabel("donor age")
    axes[1].set_title("donor age by 10x pool: R2(age ~ pool) = 0.23, p = 4e-21")
    fig.tight_layout(); fig.savefig(FIG / "phase1b_pool_age_depth_confound.png", dpi=140)


if __name__ == "__main__":
    main()
