"""STAGE P0 — Brain-specific confounders on the DLPFC adult-control cohort.

Audits obs (+ uns) for PMI, RIN/RNA quality, ambient-RNA / soup proxies, neuronal
nucleus fraction (snRNA dissociation bias), cause of death, and tissue pH.

Any donor-level expression-affecting variable with within-site Spearman |r| vs age
>= 0.20 (declared a priori) enters the P1 design alongside depth.

The 6-plex hashing pools are not in obs and are not recovered from barcodekey.

Usage: python src/brain_phase1_p0.py
"""
from __future__ import annotations

import json
import sys
from pathlib import Path

import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from scipy import stats

sys.path.insert(0, str(Path(__file__).resolve().parent))
from brain_phase1_common import (  # noqa: E402
    PHASE1_DIR, PHASE1_FIG, PHASE1_SEED, Logger, HASH_POOL_LIMITATION,
    P0_AGE_CORR_ABS, NEURONAL_CLASSES, LOW_N_GENES_CUTOFF, OBS_ABSENT_KEYWORDS,
    ADULT_MIN_AGE, DLPFC_NAME, PB_PATH, dump_json, ensure_logcpm,
)
from brain_common import BRAIN_DIR  # noqa: E402
from tissue_download_top import local_path  # noqa: E402
from tissue_peek_obs import read_obs_column, obs_column_names  # noqa: E402
from search_cxg import parse_age_years  # noqa: E402


def _decode(v):
    if isinstance(v, bytes):
        return v.decode("utf-8", "replace")
    return v


def _keyword_hits(name: str) -> list[str]:
    n = name.lower().replace(" ", "_")
    hits = []
    for kw in OBS_ABSENT_KEYWORDS:
        if kw in n:
            hits.append(kw)
    # 'ph' as a token (avoid matching Schizophrenia etc.)
    tokens = set(n.replace("-", "_").split("_"))
    if "ph" in tokens or n in {"ph", "p_h"}:
        hits.append("ph")
    return sorted(set(hits))


def spearman_safe(x, y):
    x = np.asarray(x, float)
    y = np.asarray(y, float)
    m = np.isfinite(x) & np.isfinite(y)
    if m.sum() < 8:
        return np.nan, np.nan, int(m.sum())
    if np.nanstd(x[m]) < 1e-12 or np.nanstd(y[m]) < 1e-12:
        return np.nan, np.nan, int(m.sum())
    r, p = stats.spearmanr(x[m], y[m])
    return float(r), float(p), int(m.sum())


def audit_raw_obs(h5ad_path: Path, log) -> dict:
    import h5py
    rec = dict(path=str(h5ad_path), size_bytes=int(h5ad_path.stat().st_size))
    with h5py.File(h5ad_path, "r") as f:
        cols = obs_column_names(f)
        rec["n_obs_columns"] = len(cols)
        rec["obs_columns"] = cols
        rec["uns_keys"] = list(f["uns"].keys()) if "uns" in f else []
        rec["obsm_keys"] = list(f["obsm"].keys()) if "obsm" in f else []
        rec["layers"] = list(f["layers"].keys()) if "layers" in f else []
        if "batch_condition" in f["uns"]:
            bc = f["uns"]["batch_condition"][()]
            rec["uns_batch_condition"] = [_decode(x) for x in np.atleast_1d(bc)]
        hits = []
        for c in cols:
            kw = _keyword_hits(c)
            if kw:
                hits.append(dict(column=c, keywords=kw))
        rec["keyword_hits_in_obs"] = hits
        rec["has_pmi"] = any("pmi" in c.lower() or "postmortem" in c.lower().replace("-", "")
                             or "post_mortem" in c.lower() for c in cols)
        rec["has_rin"] = any("rin" == c.lower() or "rna_integrity" in c.lower() or c.lower().endswith("_rin")
                             for c in cols)
        rec["has_ph"] = any(c.lower() in {"ph", "tissue_ph", "ph_tissue"} or "tissue_ph" in c.lower()
                            for c in cols)
        rec["has_cod"] = any("death" in c.lower() or c.lower() in {"cod", "cause_of_death"}
                             for c in cols)
        rec["has_soup"] = any(any(k in c.lower() for k in ("soup", "ambient", "decontx", "soupx"))
                              for c in cols)
        rec["has_hash_pool"] = any(any(k in c.lower() for k in ("hash_pool", "hashtag", "hashing", "multiplex"))
                                   for c in cols)
        log(f"[P0] obs columns ({len(cols)}): {cols}")
        log(f"[P0] uns keys: {rec['uns_keys']}")
        log(f"[P0] uns batch_condition: {rec.get('uns_batch_condition')}")
        log(f"[P0] PMI in obs: {rec['has_pmi']}  RIN: {rec['has_rin']}  pH: {rec['has_ph']}  "
            f"cause of death: {rec['has_cod']}  soup/ambient: {rec['has_soup']}  hash pool: {rec['has_hash_pool']}")
        log(f"[P0] keyword hits in obs column names: {hits or 'none'}")
        log("[P0] " + HASH_POOL_LIMITATION)
    return rec


def donor_level_from_raw(h5ad_path: Path, analysis_donors: pd.DataFrame, log) -> pd.DataFrame:
    """Neuronal fraction, depth, low-gene debris proxy — from ALL adult nuclei of analysis donors.

    Does not load the 12 GB count matrix.
    """
    import h5py
    log(f"[P0] reading donor-level obs from {h5ad_path.name} (no matrix)")
    with h5py.File(h5ad_path, "r") as f:
        donor = np.array([_decode(x) for x in read_obs_column(f, "donor_id")])
        klass = np.array([_decode(x) for x in read_obs_column(f, "class")])
        ct = np.array([_decode(x) for x in read_obs_column(f, "cell_type")])
        stage = np.array([_decode(x) for x in read_obs_column(f, "development_stage")])
        source = np.array([_decode(x) for x in read_obs_column(f, "Source")])
        n_counts = np.asarray(read_obs_column(f, "n_counts"), float)
        n_genes = np.asarray(read_obs_column(f, "n_genes"), float)
        sex = np.array([_decode(x) for x in read_obs_column(f, "sex")])
    age = np.array([parse_age_years(s) for s in stage], dtype=object)
    age_f = np.array([np.nan if a is None else float(a) for a in age])
    keep_d = set(analysis_donors.donor.astype(str))
    adult = np.isfinite(age_f) & (age_f >= ADULT_MIN_AGE) & np.isin(donor, list(keep_d))
    log(f"[P0] nuclei total={len(donor):,}  adult analysis-donor nuclei={int(adult.sum()):,}  "
        f"donors in mask={pd.Series(donor[adult]).nunique()}")

    d = pd.DataFrame({
        "donor": donor[adult], "class": klass[adult], "cell_type": ct[adult],
        "Source": source[adult], "n_counts": n_counts[adult], "n_genes": n_genes[adult],
        "age": age_f[adult], "sex": sex[adult],
    })
    is_neu_class = d["class"].isin(NEURONAL_CLASSES)
    is_neu_name = d["cell_type"].str.contains("neuron", case=False, na=False)
    d["neuronal_class"] = is_neu_class
    d["neuronal_name"] = is_neu_name
    d["low_n_genes"] = d.n_genes < LOW_N_GENES_CUTOFF

    rows = []
    for donor, g in d.groupby("donor", observed=True):
        n = len(g)
        rows.append(dict(
            donor=str(donor),
            n_nuclei=n,
            frac_neuronal=float(g.neuronal_class.mean()),
            frac_neuronal_by_name=float(g.neuronal_name.mean()),
            frac_EN=float((g["class"] == "EN").mean()),
            frac_IN=float((g["class"] == "IN").mean()),
            frac_Oligo=float((g["class"] == "Oligo").mean()),
            frac_Immune=float((g["class"] == "Immune").mean()),
            mean_n_counts=float(g.n_counts.mean()),
            median_n_counts=float(g.n_counts.median()),
            mean_n_genes=float(g.n_genes.mean()),
            median_n_genes=float(g.n_genes.median()),
            frac_low_n_genes=float(g.low_n_genes.mean()),
            log_mean_umi=float(np.log(max(g.n_counts.mean(), 1e-12))),
            log_mean_genes=float(np.log(max(g.n_genes.mean(), 1e-12))),
            Source=str(g.Source.iloc[0]),
            age=float(g.age.iloc[0]),
            sex=str(g.sex.iloc[0]),
        ))
    don = pd.DataFrame(rows)
    # class vs name neuronal definition should agree closely
    r_agree, p_agree, _ = spearman_safe(don.frac_neuronal, don.frac_neuronal_by_name)
    log(f"[P0] frac_neuronal (class EN+IN) vs name-contains-'neuron': Spearman r={r_agree:.3f} "
        f"(mean |diff|={(don.frac_neuronal - don.frac_neuronal_by_name).abs().mean():.4f})")
    disagree = int(((don.frac_neuronal - don.frac_neuronal_by_name).abs() > 0.05).sum())
    log(f"[P0] donors with |class−name neuronal frac| > 0.05: {disagree}/{len(don)}")
    return don


def pb_quality_proxies(log) -> pd.DataFrame:
    """MALAT1 and MT UMI fractions from the existing pseudobulk (no raw-matrix reread).

    These are soup / RNA-quality proxies computable without PMI/RIN. Weighted to donor
    level by n_cells.
    """
    import anndata as ad
    pb = ad.read_h5ad(PB_PATH)
    X = np.asarray(pb.X, dtype=np.float64)
    lib = X.sum(1)
    lib[lib == 0] = np.nan
    sym = pb.var["feature_name"].astype(str).to_numpy() if "feature_name" in pb.var else pb.var.index.astype(str).to_numpy()
    malat = np.isin(sym, ["MALAT1", "Malat1"])
    mt = np.array([s.startswith("MT-") or s.startswith("mt-") for s in sym])
    log(f"[P0] MALAT1 columns in pb var: {int(malat.sum())}  MT- columns: {int(mt.sum())}")
    malat_frac = (X[:, malat].sum(1) / lib) if malat.any() else np.full(len(pb), np.nan)
    mt_frac = (X[:, mt].sum(1) / lib) if mt.any() else np.full(len(pb), np.nan)
    tmp = pb.obs[["donor", "n_cells"]].copy()
    tmp["donor"] = tmp.donor.astype(str)
    tmp["malat1_frac"] = malat_frac
    tmp["mt_frac"] = mt_frac
    tmp["cc_umi_frac"] = pb.obs["cc_umi_frac"].to_numpy(float) if "cc_umi_frac" in pb.obs else np.nan
    tmp["n_cells"] = tmp.n_cells.to_numpy(float)

    def wavg(g, col):
        w = g.n_cells.to_numpy(float)
        v = g[col].to_numpy(float)
        m = np.isfinite(v) & np.isfinite(w)
        if not m.any() or w[m].sum() <= 0:
            return np.nan
        return float(np.average(v[m], weights=w[m]))

    rows = []
    for don, g in tmp.groupby("donor"):
        rows.append(dict(
            donor=str(don),
            malat1_frac=wavg(g, "malat1_frac"),
            mt_frac=wavg(g, "mt_frac"),
            cc_umi_frac=wavg(g, "cc_umi_frac"),
            n_pb=len(g),
        ))
    return pd.DataFrame(rows)


def correlation_table(don: pd.DataFrame, covariates: list[str], log) -> pd.DataFrame:
    rows = []
    sites = sorted(don.Source.astype(str).unique())
    for c in covariates:
        if c not in don.columns:
            continue
        for site in sites + ["within_site_pooled"]:
            if site == "within_site_pooled":
                x = don[c].to_numpy(float)
                y = don.age.to_numpy(float)
                # residualize on site
                D = pd.get_dummies(don.Source.astype(str), dtype=float).to_numpy()
                Q, _ = np.linalg.qr(D)
                x = x - Q @ (Q.T @ np.nan_to_num(x, nan=np.nanmean(x)))
                y = y - Q @ (Q.T @ y)
                r, p, n = spearman_safe(x, y)
            else:
                g = don[don.Source.astype(str) == site]
                r, p, n = spearman_safe(g[c], g.age)
            rows.append(dict(covariate=c, site=site, spearman_r=r, p=p, n=n,
                             abs_r=abs(r) if r == r else np.nan))
    tab = pd.DataFrame(rows)
    log("\n[P0] within-site Spearman(age, covariate):")
    piv = tab.pivot(index="covariate", columns="site", values="spearman_r")
    log(piv.to_string(float_format=lambda x: f"{x:+.3f}"))
    return tab


def decide_included(corr: pd.DataFrame, candidates: list[str], always: list[str], log) -> list[str]:
    included = list(always)
    reasons = {c: "always in design (specified: log depth)" for c in always}
    for c in candidates:
        sub = corr[corr.covariate == c]
        hit = sub[sub.abs_r >= P0_AGE_CORR_ABS]
        if len(hit):
            if c not in included:
                included.append(c)
            reasons[c] = (
                f"within-site |Spearman r| >= {P0_AGE_CORR_ABS} "
                f"({', '.join(f'{r.site} r={r.spearman_r:+.3f}' for r in hit.itertuples())})"
            )
        else:
            reasons.setdefault(c, f"not included: all within-site |r| < {P0_AGE_CORR_ABS}")
    log(f"\n[P0] inclusion threshold |r| >= {P0_AGE_CORR_ABS} (declared a priori)")
    for c, why in reasons.items():
        log(f"   {c}: {why}")
    return included, reasons


def figures(don: pd.DataFrame, covariates: list[str]):
    covs = [c for c in covariates if c in don.columns]
    n = len(covs)
    if not n:
        return
    fig, axes = plt.subplots(int(np.ceil(n / 3)), 3, figsize=(11, 3.2 * np.ceil(n / 3)), squeeze=False)
    colors = {"H": "tab:blue", "M": "tab:orange"}
    for i, c in enumerate(covs):
        ax = axes.flat[i]
        for src, g in don.groupby(don.Source.astype(str)):
            ax.scatter(g.age, g[c], s=18, alpha=0.7, c=colors.get(src, "grey"), label=f"{src} n={len(g)}")
        ax.set_xlabel("donor age (y)")
        ax.set_ylabel(c)
        ax.legend(fontsize=7, loc="best")
    for j in range(i + 1, axes.size):
        axes.flat[j].axis("off")
    fig.suptitle("P0 donor-level confounders vs age (H=HBCC, M=MSSM)", fontsize=11)
    fig.tight_layout()
    fig.savefig(PHASE1_FIG / "p0_covariates_vs_age.png", dpi=140)
    plt.close(fig)


def run():
    rng = np.random.default_rng(PHASE1_SEED)
    log = Logger(PHASE1_DIR / "p0_report.txt")
    try:
        return _run(log, rng)
    finally:
        log.close()


def _run(log, rng):
    log("=" * 100)
    log("BRAIN PHASE 1  STAGE P0 — confounders.  seed =", PHASE1_SEED)
    log("=" * 100)
    log(HASH_POOL_LIMITATION)
    log(f"Adult cutoff age >= {ADULT_MIN_AGE:.0f}. Inclusion |r| >= {P0_AGE_CORR_ABS} (a priori).")
    log(f"Neuronal classes (a priori): {NEURONAL_CLASSES}. Low-n_genes debris cutoff: {LOW_N_GENES_CUTOFF}.")

    donors_r1 = pd.read_csv(BRAIN_DIR / "r1_donors.csv")
    donors_r1["donor"] = donors_r1["donor"].astype(str)
    log(f"[P0] R1 analysis donors: {len(donors_r1)}  ages {donors_r1.age.min():.0f}–{donors_r1.age.max():.0f}  "
        f"Source={donors_r1.Source.astype(str).value_counts().to_dict()}")

    h5ad = local_path(DLPFC_NAME)
    if not h5ad.exists():
        raise FileNotFoundError(h5ad)
    audit = audit_raw_obs(h5ad, log)
    don = donor_level_from_raw(h5ad, donors_r1, log)
    qprox = pb_quality_proxies(log)
    don = don.merge(qprox, on="donor", how="left")
    # attach R2 depth already computed on this cohort
    r2d = pd.read_csv(BRAIN_DIR / "r2_donors.csv")
    r2d["donor"] = r2d["donor"].astype(str)
    if "mean_counts_per_cell" in r2d:
        don = don.merge(r2d[["donor", "mean_counts_per_cell", "mean_genes_per_cell"]], on="donor", how="left")

    log(f"[P0] donor table n={len(don)}  missing malat1={int(don.malat1_frac.isna().sum())}  "
        f"missing mt={int(don.mt_frac.isna().sum())}")
    log(f"[P0] frac_neuronal: mean={don.frac_neuronal.mean():.3f}  "
        f"range {don.frac_neuronal.min():.3f}–{don.frac_neuronal.max():.3f}")
    log(f"[P0] malat1_frac: mean={don.malat1_frac.mean():.4f}  mt_frac: mean={don.mt_frac.mean():.5f}")
    log(f"[P0] frac_low_n_genes (n_genes<{LOW_N_GENES_CUTOFF}): mean={don.frac_low_n_genes.mean():.3f}")

    # Covariates to audit. Depth stays in the model regardless.
    candidates = [
        "frac_neuronal", "frac_low_n_genes", "malat1_frac", "mt_frac",
        "mean_n_counts", "mean_n_genes", "n_nuclei", "cc_umi_frac",
        "frac_Oligo", "frac_Immune",
    ]
    always_in_model = ["log_mean_umi", "log_mean_genes"]  # represented as log(mean_counts/genes) in decomp
    corr = correlation_table(don, candidates + always_in_model, log)
    corr.to_csv(PHASE1_DIR / "p0_covariate_age_correlations.csv", index=False)

    # Always-in-design depth is implemented as log(mean_counts_per_cell), log(mean_genes_per_cell)
    # on the pseudobulk — not re-added as extra. Extra = P0 hits among candidates.
    included, reasons = decide_included(corr, candidates, always=[], log=log)
    # depth always noted
    log("[P0] depth covariates ALWAYS in P1 design: log(mean UMI/nucleus), log(mean genes/nucleus) "
        "(on each pseudobulk; specified, not selected on |r|)")

    dump_json(PHASE1_DIR / "p0_included_covariates.json", dict(
        included=included,
        always_in_design=["log_mean_umi_per_nucleus", "log_mean_genes_per_nucleus"],
        inclusion_threshold_abs_r=P0_AGE_CORR_ABS,
        reasons=reasons,
        absent=dict(PMI=not audit["has_pmi"], RIN=not audit["has_rin"],
                    tissue_pH=not audit["has_ph"], cause_of_death=not audit["has_cod"],
                    soupX_decontX=not audit["has_soup"], hashing_pools=not audit["has_hash_pool"]),
        hash_pool_limitation=HASH_POOL_LIMITATION,
        seed=PHASE1_SEED,
        n_donors=int(len(don)),
        neuronal_classes=list(NEURONAL_CLASSES),
        low_n_genes_cutoff=LOW_N_GENES_CUTOFF,
    ))
    don.to_csv(PHASE1_DIR / "p0_donor_covariates.csv", index=False)
    dump_json(PHASE1_DIR / "p0_obs_audit.json", audit)
    figures(don, ["frac_neuronal", "malat1_frac", "mt_frac", "frac_low_n_genes",
                  "mean_n_genes", "n_nuclei", "cc_umi_frac"])

    # filter log for this stage (no nuclei dropped; audit only)
    filt = [
        dict(step="P0 audit (no additional cell filter)", unit="donors",
             before=int(len(donors_r1)), after=int(len(don)), removed=int(len(donors_r1) - len(don)),
             note="P0 does not drop donors; neuronal fraction uses adult nuclei of R1 donors"),
    ]
    pd.DataFrame(filt).to_csv(PHASE1_DIR / "p0_filter_log.csv", index=False)

    log(f"\n[P0] PMI present in obs: {audit['has_pmi']}  — ABSENT is a real limitation.")
    log(f"[P0] RIN present in obs: {audit['has_rin']}  — ABSENT is a real limitation.")
    log(f"[P0] tissue pH present: {audit['has_ph']}  cause of death: {audit['has_cod']}")
    log(f"[P0] included extra covariates for P1: {included or '(none beyond specified depth)'}")
    log(f"[P0] wrote {PHASE1_DIR}")
    log(f"[P0] rng consumed (seed {PHASE1_SEED}); no stochastic draws in P0.")
    dump_json(PHASE1_DIR / "p0_summary.json", dict(
        n_donors=int(len(don)),
        included_extra=included,
        has_pmi=audit["has_pmi"], has_rin=audit["has_rin"], has_ph=audit["has_ph"],
        has_cod=audit["has_cod"], has_soup=audit["has_soup"], has_hash_pool=audit["has_hash_pool"],
        frac_neuronal_mean=float(don.frac_neuronal.mean()),
        seed=PHASE1_SEED,
    ))
    from brain_phase1_findings import write_findings
    write_findings(log=log)
    return dict(included=included, n_donors=len(don))


if __name__ == "__main__":
    run()
