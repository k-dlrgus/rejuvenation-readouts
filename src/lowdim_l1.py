"""STAGE L1 — build the low-dimensional spaces and flag technical components.

Per cell type, on the 233-donor DLPFC pseudobulks. k in {20, 50, 100, 150}.

  A  unsupervised PCA, descriptive full-cohort fit (evaluation in L2 is
     nested inside training folds).
  B  curated Hallmark/GO/Reactome module scores; k modules ranked by
     mean within-type variance (Hallmark first).
  C  WGCNA-style co-expression modules (Ward clustering of 2000 HVGs
     in sample-space; average-linkage on correlation chained to
     singletons and was not used).

CRITICAL CONTROL: any component with |Spearman r| > 0.5 vs site, log depth
(mean UMI/nucleus, mean genes/nucleus), or neuronal fraction is flagged.
L2 reports results with and without those components.

This script's PCA/modules/WGCNA fits are DESCRIPTIVE (full cohort) for the
technical-correlation audit. They are NOT used as frozen features in L2.

Usage: python src/lowdim_l1.py
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
from lowdim_common import (  # noqa: E402
    LD_DIR, LD_FIG, LD_SEED, DATASET_ID, KS, KMAX, VARIANTS, VARIANT_LABEL,
    TECH_R_FLAG, TECH_NAMES, MIN_SET_GENES, MIN_N_TYPE,
    Logger, dump_json, ld_log_banner, load_bundle, load_gene_sets,
    zscore_train, fit_spaces_by_type, raw_scores, flag_components, FilterLog,
)


def _run(log, rng):
    ld_log_banner(log, "L1 — low-dimensional spaces + technical control")
    flog = FilterLog()
    data = load_bundle(log)
    Y, genes, obs, tech, types = data["Y"], data["genes"], data["obs"], data["tech"], data["types"]
    flog.add("L1", "pseudobulks", len(obs), n_genes=int(Y.shape[1]), n_donors=int(data["n_donors"]))
    flog.add("L1", "cell_types", len(types))
    sets, set_filt = load_gene_sets(genes, log)
    n_usable = sum(1 for s in sets if s["n_mapped"] >= MIN_SET_GENES)
    flog.add("L1", "gene_sets_listed", len(sets), usable=n_usable, min_genes=MIN_SET_GENES)

    Xs, mu, sd = zscore_train(Y)
    log(f"[L1] global z-score on full cohort (DESCRIPTIVE only): mu finite={np.isfinite(mu).mean():.3f}")
    ct = obs.celltype.astype(str).to_numpy()

    long_rows = []
    sum_rows = []
    for variant in VARIANTS:
        log(f"\n[L1] variant={variant} ({VARIANT_LABEL[variant]})  kmax={KMAX}  DESCRIPTIVE")
        spaces = fit_spaces_by_type(Xs, Y, ct, types, variant, KMAX, sets)
        for t in types:
            sp = spaces.get(t)
            n = int((ct == t).sum())
            if sp is None or n < MIN_N_TYPE:
                log(f"   {t}: SKIP n={n}")
                continue
            X_t = Xs[ct == t]
            tech_t = tech[ct == t]
            for k in KS:
                spec = sp.spec(k)
                k_req = len(spec["names"])
                if k_req < 1:
                    sum_rows.append(dict(
                        variant=variant, celltype=t, k=k, k_eff=0, n=n,
                        n_flagged=0, frac_flagged=np.nan, n_flag_site=0,
                        n_flag_depth=0, n_flag_nfrac=0, reverted=False,
                    ))
                    continue
                S = raw_scores(X_t, spec)
                R, flags = flag_components(S, tech_t)
                n_flag = int(flags.sum())
                n_site = int((np.abs(R[:, 0]) > TECH_R_FLAG).sum())
                n_depth = int(((np.abs(R[:, 1]) > TECH_R_FLAG) | (np.abs(R[:, 2]) > TECH_R_FLAG)).sum())
                n_nfrac = int((np.abs(R[:, 3]) > TECH_R_FLAG).sum())
                log(f"   {t} k={k} k_eff={k_req}: flagged={n_flag}/{k_req}  "
                    f"site={n_site} depth={n_depth} nfrac={n_nfrac}")
                sum_rows.append(dict(
                    variant=variant, celltype=t, k=int(k), k_eff=int(k_req), n=n,
                    n_flagged=n_flag, frac_flagged=n_flag / max(k_req, 1),
                    n_flag_site=n_site, n_flag_depth=n_depth, n_flag_nfrac=n_nfrac,
                    median_abs_r_site=float(np.nanmedian(np.abs(R[:, 0]))),
                    median_abs_r_logumi=float(np.nanmedian(np.abs(R[:, 1]))),
                    median_abs_r_loggenes=float(np.nanmedian(np.abs(R[:, 2]))),
                    median_abs_r_nfrac=float(np.nanmedian(np.abs(R[:, 3]))),
                    max_abs_r_any=float(np.nanmax(np.abs(R))) if np.isfinite(R).any() else np.nan,
                    n_hvg=int(sp.extra.get("n_hvg", np.nan)) if variant == "wgcna" else np.nan,
                    n_valid_sets=int(sp.extra.get("n_valid_sets", np.nan)) if variant == "curated" else np.nan,
                ))
                for j, name in enumerate(spec["names"]):
                    rec = dict(
                        variant=variant, celltype=t, k=int(k), component=j + 1,
                        name=name, flagged=bool(flags[j]),
                    )
                    for ti, tn in enumerate(TECH_NAMES):
                        rec[f"r_{tn}"] = float(R[j, ti]) if np.isfinite(R[j, ti]) else np.nan
                        rec[f"flag_{tn}"] = bool(np.isfinite(R[j, ti]) and abs(R[j, ti]) > TECH_R_FLAG)
                    long_rows.append(rec)

    long_df = pd.DataFrame(long_rows)
    sum_df = pd.DataFrame(sum_rows)
    long_df.to_csv(LD_DIR / "l1_component_tech.csv", index=False)
    sum_df.to_csv(LD_DIR / "l1_flag_summary.csv", index=False)
    flog.add("L1", "component_rows", len(long_df))
    flog.to_csv(LD_DIR / "l1_filter_log.csv")

    # median-over-types fraction flagged
    log("\n[L1] median-over-types fraction of components flagged (|r|>0.5):")
    overview = []
    for variant in VARIANTS:
        for k in KS:
            sub = sum_df[(sum_df.variant == variant) & (sum_df.k == k)]
            if not len(sub):
                continue
            rec = dict(
                variant=variant, label=VARIANT_LABEL[variant], k=int(k),
                median_frac_flagged=float(sub.frac_flagged.median()),
                median_n_flagged=float(sub.n_flagged.median()),
                median_k_eff=float(sub.k_eff.median()),
                types_any_flagged=int((sub.n_flagged > 0).sum()),
                n_types=int(len(sub)),
                median_abs_r_site=float(sub.median_abs_r_site.median()),
                median_abs_r_logumi=float(sub.median_abs_r_logumi.median()),
                median_abs_r_nfrac=float(sub.median_abs_r_nfrac.median()),
            )
            overview.append(rec)
            log(f"   {VARIANT_LABEL[variant]} k={k}: frac_flagged={rec['median_frac_flagged']:.3f}  "
                f"|r_site|={rec['median_abs_r_site']:.3f}  "
                f"|r_umi|={rec['median_abs_r_logumi']:.3f}  "
                f"|r_nfrac|={rec['median_abs_r_nfrac']:.3f}  "
                f"types_with_flag={rec['types_any_flagged']}/{rec['n_types']}")
    ov = pd.DataFrame(overview)
    ov.to_csv(LD_DIR / "l1_overview.csv", index=False)

    _figures(sum_df, long_df)
    summary = dict(
        seed=LD_SEED, dataset_id=DATASET_ID,
        n_donors=int(data["n_donors"]), n_types=len(types), n_genes=int(Y.shape[1]),
        n_rows=int(len(obs)), k=list(KS), variants=list(VARIANTS),
        tech_r_flag=TECH_R_FLAG, n_gene_sets=len(sets), n_gene_sets_usable=n_usable,
        gene_set_filter=set_filt, overview=ov.to_dict(orient="records"),
        note="L1 fits are descriptive (full cohort) for the technical audit; L2 refits inside training folds.",
    )
    dump_json(LD_DIR / "l1_summary.json", summary)
    log("[L1] wrote results/lowdim/l1_*")
    log("[L1] done.")
    return summary


def _figures(sum_df, long_df):
    fig, ax = plt.subplots(figsize=(7.2, 4.2))
    colors = dict(pca="tab:blue", curated="tab:green", wgcna="tab:orange")
    for variant in VARIANTS:
        sub = sum_df[sum_df.variant == variant]
        if not len(sub):
            continue
        g = sub.groupby("k").frac_flagged.median()
        ax.plot(g.index, g.values, "-o", color=colors[variant], label=VARIANT_LABEL[variant])
    ax.set_xlabel("k")
    ax.set_ylabel("median-over-types fraction of components flagged")
    ax.set_title("L1  technical components (|r|>0.5 with site / depth / nfrac)")
    ax.legend(fontsize=8)
    ax.set_ylim(-0.02, 1.02)
    fig.tight_layout()
    fig.savefig(LD_FIG / "l1_frac_flagged.png", dpi=140)
    plt.close(fig)

    fig, axes = plt.subplots(1, 3, figsize=(11.4, 3.6), sharey=True)
    for ax, tn, title in zip(axes, ("site", "log_umi", "frac_neuronal"),
                             ("site (MSSM dummy)", "log UMI/nucleus", "neuronal fraction")):
        col = f"r_{tn}"
        for variant, c in colors.items():
            sub = long_df[(long_df.variant == variant) & (long_df.k == 50)]
            if not len(sub) or col not in sub:
                continue
            vals = sub[col].abs().dropna().to_numpy()
            if len(vals):
                ax.hist(vals, bins=25, histtype="step", color=c, label=VARIANT_LABEL[variant], density=True)
        ax.axvline(TECH_R_FLAG, c="k", ls="--", lw=1)
        ax.set_title(title, fontsize=9)
        ax.set_xlabel("|Spearman r|")
    axes[0].set_ylabel("density (k=50 components)")
    axes[2].legend(fontsize=7)
    fig.suptitle("L1  component–technical correlations (k=50, all types)", fontsize=10)
    fig.tight_layout()
    fig.savefig(LD_FIG / "l1_tech_r_hist_k50.png", dpi=140)
    plt.close(fig)

    # per-type flagged count at k=50
    fig, ax = plt.subplots(figsize=(8.4, 5.6))
    k50 = sum_df[sum_df.k == 50]
    types = sorted(k50.celltype.unique())
    y_pos = np.arange(len(types))
    width = 0.25
    for i, variant in enumerate(VARIANTS):
        sub = k50[k50.variant == variant].set_index("celltype").reindex(types)
        ax.barh(y_pos + (i - 1) * width, sub.frac_flagged.fillna(0), height=width,
                color=colors[variant], label=VARIANT_LABEL[variant])
    ax.set_yticks(y_pos)
    ax.set_yticklabels(types, fontsize=6)
    ax.set_xlabel("fraction of k=50 components flagged")
    ax.set_title("L1  per-type technical fraction (k=50)")
    ax.legend(fontsize=7)
    ax.set_xlim(0, 1.05)
    fig.tight_layout()
    fig.savefig(LD_FIG / "l1_per_type_frac_k50.png", dpi=140)
    plt.close(fig)


def run():
    rng = np.random.default_rng(LD_SEED)
    log = Logger(LD_DIR / "l1_report.txt")
    try:
        return _run(log, rng)
    finally:
        log.close()


if __name__ == "__main__":
    run()
