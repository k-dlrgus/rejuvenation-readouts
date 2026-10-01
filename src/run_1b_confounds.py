"""Phase 1b runner: flag confounds, report A-gene survival per filter, evaluate STOP CONDITION 1b."""
import sys
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

sys.path.insert(0, str(__import__("pathlib").Path(__file__).resolve().parent))
from config import DATA_PROC, TAB, FIG, SEED  # noqa
from decomp import classify, SWEEP, PRIMARY_THRESH  # noqa
from confounds import build_flags, covariate_age_report  # noqa
from ensembl_chrom import fetch_chromosomes  # noqa

np.random.seed(SEED)


def survival_table(genes_A: pd.Index, F: pd.DataFrame):
    """Sequential and marginal removal counts for an A-gene set."""
    fa = F.loc[genes_A]
    n0 = len(fa)
    steps = [("cell-cycle (Tirosh S/G2M + proliferation list)", fa.flag_cc_list),
             ("cell-cycle (data: |r|>=0.3 with cycling covariates)", fa.flag_cc_data),
             ("ribosomal (RPL/RPS/MRP)", fa.flag_ribo),
             ("mitochondrial (MT-)", fa.flag_mito),
             ("sex-linked (chrX/chrY)", fa.flag_sex_chrom),
             ("sex-DE (|d|>=0.5, p<1e-6, within cell type)", fa.flag_sex_data & ~fa.flag_sex_chrom),
             ("depth-correlated (|r|>=0.3 with mean counts/genes per cell)", fa.flag_depth)]
    rows, removed = [], pd.Series(False, index=fa.index)
    for name, flag in steps:
        marg = int(flag.sum())
        new = int((flag & ~removed).sum())
        removed |= flag
        rows.append(dict(filter=name, flagged_marginal=marg, removed_incremental=new, surviving=n0 - int(removed.sum())))
    t = pd.DataFrame(rows)
    return t, removed


def main():
    Y = np.load(DATA_PROC / "onek1k_logcpm.npy")
    genes = pd.read_csv(DATA_PROC / "onek1k_logcpm_genes.csv")
    obs = pd.read_csv(DATA_PROC / "onek1k_pseudobulk_obs.csv", index_col=0)
    df = pd.read_csv(TAB / "phase1a_variance_decomposition.csv", index_col=0)
    assert (df.index == genes.gene_id).all()
    chrom = fetch_chromosomes(genes.gene_id.tolist(), verbose=False).set_index("gene_id").chrom

    print("=== covariates vs age (is proliferation / depth itself age-associated?) ===")
    cov_rep = covariate_age_report(obs)
    cov_rep.to_csv(TAB / "phase1b_covariate_age_correlations.csv", index=False)
    don = obs.groupby("donor").agg(age=("age", "first"), sex=("sex", "first"))
    from scipy import stats
    t = stats.ttest_ind(don.age[don.sex == "male"], don.age[don.sex == "female"], equal_var=False)
    print(f"\nsex vs age across donors: male {don.age[don.sex=='male'].mean():.1f}y (n={int((don.sex=='male').sum())}), "
          f"female {don.age[don.sex=='female'].mean():.1f}y (n={int((don.sex=='female').sum())}), Welch p={t.pvalue:.3f} -> "
          f"{'NOT confounded' if t.pvalue > 0.05 else 'CONFOUNDED'}")

    F = build_flags(Y, obs, genes.gene_id.to_numpy(), genes.symbol.to_numpy(), chrom)
    F = F.join(df[["unique_age", "unique_ct", "interaction", "age_slope_per_year", "gene_class_primary"]])
    F["age_signal_retained_within_pool"] = F.unique_age_within_pool / F.unique_age.replace(0, np.nan)
    F.to_csv(TAB / "phase1b_gene_flags.csv")
    print("\n=== genome-wide flag counts (all 13,904 expressed genes) ===")
    for c in ["flag_cc_list", "flag_cc_data", "flag_cell_cycle", "flag_ribo", "flag_mito", "flag_sex_chrom", "flag_sex_data", "flag_sex", "flag_depth", "flag_technical"]:
        print(f"  {c:<18} {int(F[c].sum()):>6}")

    # ---- PRIMARY A-genes ----
    A = F.index[F.gene_class_primary == "A"]
    print(f"\n=== PRIMARY A-genes: n={len(A)} (thresholds {PRIMARY_THRESH}) ===")
    t, removed = survival_table(A, F)
    print(t.to_string(index=False))
    t.to_csv(TAB / "phase1b_Agene_survival_primary.csv", index=False)
    n_cc = int(F.loc[A, "flag_cell_cycle"].sum())
    n_tech = int((F.loc[A, "flag_technical"] & ~F.loc[A, "flag_cell_cycle"]).sum())
    n_cc_or_tech = int((F.loc[A, "flag_cell_cycle"] | F.loc[A, "flag_technical"]).sum())
    n_sex = int((F.loc[A, "flag_sex"] & ~(F.loc[A, "flag_cell_cycle"] | F.loc[A, "flag_technical"])).sum())
    n_surv = len(A) - int(removed.sum())
    print(f"\nA-genes removed as cell-cycle: {n_cc}; additionally as technical (ribo/mito/depth): {n_tech}; "
          f"cell-cycle OR technical: {n_cc_or_tech}/{len(A)} = {100*n_cc_or_tech/len(A):.1f}%")
    print(f"A-genes removed as sex-linked only (biological covariate, reported separately): {n_sex}")
    print(f"A-genes surviving ALL filters: {n_surv}/{len(A)} = {100*n_surv/len(A):.1f}%")
    print("Flagged A-genes:")
    fl = F.loc[A][F.loc[A, ["flag_cell_cycle", "flag_technical", "flag_sex"]].any(axis=1)]
    for gid, r in fl.sort_values("unique_age", ascending=False).iterrows():
        tags = [k.replace("flag_", "") for k in ["flag_cc_list", "flag_cc_data", "flag_ribo", "flag_mito", "flag_sex_chrom", "flag_sex_data", "flag_depth"] if r[k]]
        print(f"   {r.symbol:<14} age={r.unique_age:.4f} ct={r.unique_ct:.3f} r_cc={r.r_cc_umi_frac:+.2f} r_prolif={r.r_prolif_cell_frac:+.2f} "
              f"r_depth={r.r_mean_counts:+.2f}/{r.r_mean_genes:+.2f} d_sex={r.sex_cohen_d:+.2f} chr={r.chrom} [{', '.join(tags)}]")
    surv = F.loc[A][~removed]
    print("\nSurviving A-genes:", ", ".join(surv.sort_values("unique_age", ascending=False).symbol))
    print("\nBatch check on A-genes: fraction of age signal retained after adjusting for pool:")
    print(F.loc[A, "age_signal_retained_within_pool"].describe().to_string(float_format=lambda x: f"{x:.3f}"))
    print(f"   A-genes with <50% of age signal retained within pool: {int((F.loc[A,'age_signal_retained_within_pool']<0.5).sum())}")
    print(f"   A-genes with <25% of age signal retained within pool: {int((F.loc[A,'age_signal_retained_within_pool']<0.25).sum())}")

    stop_1b_literal = n_cc_or_tech > len(A) / 2
    print(f"\nSTOP 1b, literal pre-specified flag criterion only (cell-cycle | ribo | mito | depth |r|>=0.3): "
          f"{'would fire' if stop_1b_literal else 'would NOT fire on its own'} ({n_cc_or_tech}/{len(A)} = {100*n_cc_or_tech/len(A):.1f}%)")
    print("   -> final STOP 1b verdict is issued by run_1b_verdict.py after the batch/depth-mediation analysis.")

    # ---- sweep: survival across all A-gene definitions ----
    rows = []
    for t_age in SWEEP["t_age"]:
        for k_ct in SWEEP["k_ct"]:
            lab = classify(df, t_age, k_ct, PRIMARY_THRESH["t_ct"], PRIMARY_THRESH["k_age"])
            A_ = lab.index[lab == "A"]
            if len(A_) == 0:
                rows.append(dict(t_age=t_age, k_ct=k_ct, n_A=0)); continue
            fa = F.loc[A_]
            rows.append(dict(t_age=t_age, k_ct=k_ct, n_A=len(A_),
                             n_cell_cycle=int(fa.flag_cell_cycle.sum()), n_cc_list=int(fa.flag_cc_list.sum()), n_cc_data=int(fa.flag_cc_data.sum()),
                             n_ribo=int(fa.flag_ribo.sum()), n_mito=int(fa.flag_mito.sum()), n_depth=int(fa.flag_depth.sum()),
                             n_sex=int(fa.flag_sex.sum()),
                             n_cc_or_tech=int((fa.flag_cell_cycle | fa.flag_technical).sum()),
                             frac_cc_or_tech=float((fa.flag_cell_cycle | fa.flag_technical).mean()),
                             n_surviving=int((~(fa.flag_cell_cycle | fa.flag_technical | fa.flag_sex)).sum())))
    sw = pd.DataFrame(rows)
    sw.to_csv(TAB / "phase1b_Agene_survival_sweep.csv", index=False)
    print("\n=== sweep: A-gene confound removal across thresholds ===")
    print(sw.to_string(index=False, float_format=lambda x: f"{x:.3f}"))

    # I-genes flags (for completeness; I-genes are not filtered for cell-cycle per spec, but report)
    I = F.index[F.gene_class_primary == "I"]
    fi = F.loc[I]
    print(f"\nPRIMARY I-genes n={len(I)}: cell-cycle {int(fi.flag_cell_cycle.sum())}, ribo {int(fi.flag_ribo.sum())}, mito {int(fi.flag_mito.sum())}, "
          f"sex {int(fi.flag_sex.sum())}, depth {int(fi.flag_depth.sum())}")

    # ---- figure ----
    fig, axes = plt.subplots(1, 2, figsize=(12, 4.5))
    ax = axes[0]
    ax.bar(range(len(t)), t.surviving, color="steelblue")
    ax.set_xticks(range(len(t))); ax.set_xticklabels([s.split(" (")[0] for s in t["filter"]], rotation=35, ha="right", fontsize=8)
    ax.axhline(len(A), color="k", ls="--", lw=0.8); ax.axhline(len(A) / 2, color="tab:red", ls=":", lw=0.8)
    ax.set_ylabel("A-genes surviving (cumulative)"); ax.set_title(f"primary A-genes n={len(A)}; red = STOP-1b majority line")
    ax = axes[1]
    fa = F.loc[A]
    ax.scatter(fa.r_cc_umi_frac, fa.r_mean_counts, c=np.where(fa.flag_cell_cycle | fa.flag_technical, "tab:red", "grey"), s=14)
    ax.axvline(0.3, ls=":", c="k"); ax.axvline(-0.3, ls=":", c="k"); ax.axhline(0.3, ls=":", c="k"); ax.axhline(-0.3, ls=":", c="k")
    ax.set_xlabel("r(gene, cell-cycle UMI fraction) within cell type"); ax.set_ylabel("r(gene, mean UMI per cell) within cell type")
    ax.set_title("A-genes: proliferation vs depth correlation")
    fig.tight_layout(); fig.savefig(FIG / "phase1b_Agene_confound_filters.png", dpi=140)
    return F


if __name__ == "__main__":
    main()
