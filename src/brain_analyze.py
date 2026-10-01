"""DLPFC Aging_Cohort R1–R3: adult-control cohort, batch audit, turnover test.

Reuses Stage B / T3 metric code from src/tissue_analyze.py (kernel ridge, batch-grouped
CV, within-batch R², CLR composition, pool-cluster bootstrap). Seed 20260914.

Usage: python src/brain_analyze.py
"""
import hashlib
import json
import sys
from pathlib import Path

import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

sys.path.insert(0, str(Path(__file__).resolve().parent))
from brain_common import (  # noqa: E402
    BRAIN_DIR, BRAIN_FIG, BRAIN_SEED, Logger, DLPFC_NAME, ADULT_MIN_AGE,
    STOP_R1_MIN_DONORS, STOP_R1_MIN_SPAN, BLOOD_COMP_R2, BLOOD_EXPR_R2,
    WM_COMP_R2, WM_EXPR_R2,
)
from tissue_download_top import resolve, local_path  # noqa: E402
from tissue_analyze import (  # noqa: E402
    pool_center, within_pool_r2_mae, raw_r2_mae, make_pool_folds,
    kernel_ridge_oof, bootstrap_ci, _safe_ols_r2, build_composition,
    N_REPEATS, N_OUTER, N_INNER, N_BOOT, MIN_CELLS, ALPHAS,
)
from search_cxg import parse_age_years  # noqa: E402
from pseudobulk import build, FILTER_LOG  # noqa: E402
from decomp import normalize_logcpm  # noqa: E402

# Forced batch: R0 showed no sequencing pool/library/run. Source (H=HBCC, M=MSSM)
# is the only shared-donor site factor. Do NOT let pick_batch choose ancestry/sex.
FORCE_BATCH = "Source"
NOT_BATCH = {
    "sex", "sex_ontology_term_id", "genetic_ancestry", "self_reported_ethnicity",
    "self_reported_ethnicity_ontology_term_id", "disease", "disease_ontology_term_id",
    "development_stage", "development_stage_ontology_term_id", "age", "donor",
    "n_cells", "mean_counts_per_cell", "mean_genes_per_cell", "n_cells_comp",
    "class", "subclass", "subtype", "cell_type", "cell_type_ontology_term_id",
    "is_primary_data",
}
BRAIN_DISEASE_COLS = (
    "Schizophrenia", "Bipolar_Disorder", "Parkinson_disease",
    "DLBD_status", "FTD_status", "Tauopathy_status", "Vascular_status",
)


def _name_salt(name: str) -> int:
    """Stable namespace for default_rng; do not use hash() (PYTHONHASHSEED randomizes it)."""
    return int(hashlib.md5(name.encode("utf-8")).hexdigest()[:8], 16) % (2**31)


def kernel_ridge_oof_with_donor_inner(K, y, pool, donor, fold_of_pool, alphas, rng, n_inner=N_INNER):
    """Same outer batch-grouped CV as tissue_analyze.kernel_ridge_oof.

    When the training set has <2 batches (leave-one-site-out with 2 sites), batch-grouped
    inner CV cannot score alphas (all-zero scores → accidental alpha=0.01). Fall back to
    donor-grouped inner CV *inside the training site* so alpha is selected. Outer split
    remains batch-grouped; donor and batch leakage asserts still fire.
    """
    n = len(y)
    fold = np.array([fold_of_pool[p] for p in pool])
    pred = np.full(n, np.nan)
    chosen = []
    inner_mode = []
    for k in np.unique(fold):
        te = np.flatnonzero(fold == k)
        tr = np.flatnonzero(fold != k)
        assert not (set(donor[tr]) & set(donor[te])), "donor leakage"
        assert not (set(pool[tr]) & set(pool[te])), "pool leakage"
        tr_pools = np.array(sorted(set(pool[tr])))
        score = np.zeros(len(alphas))
        used_donor_inner = False
        if len(tr_pools) >= 2:
            rng.shuffle(tr_pools)
            n_in = max(2, min(n_inner, len(tr_pools)))
            inner_of_pool = {p: i % n_in for i, p in enumerate(tr_pools)}
            inner_fold = np.array([inner_of_pool[p] for p in pool[tr]])
            for ki in range(n_in):
                ite = tr[inner_fold == ki]
                itr = tr[inner_fold != ki]
                if len(itr) < 3 or len(ite) < 1:
                    continue
                mu_i = y[itr].mean()
                S, U = np.linalg.eigh(K[np.ix_(itr, itr)])
                Uy = U.T @ (y[itr] - mu_i)
                Kte = K[np.ix_(ite, itr)]
                for a_i, a in enumerate(alphas):
                    dual = U @ (Uy / (S + a))
                    p_ = mu_i + Kte @ dual
                    score[a_i] += -((y[ite] - p_) ** 2).sum()
        if not np.any(score != 0):
            used_donor_inner = True
            tr_don = np.array(sorted(set(donor[tr])))
            rng.shuffle(tr_don)
            n_in = max(2, min(n_inner, len(tr_don)))
            inner_of_don = {d: i % n_in for i, d in enumerate(tr_don)}
            inner_fold = np.array([inner_of_don[d] for d in donor[tr]])
            for ki in range(n_in):
                ite = tr[inner_fold == ki]
                itr = tr[inner_fold != ki]
                if len(itr) < 3 or len(ite) < 1:
                    continue
                mu_i = y[itr].mean()
                S, U = np.linalg.eigh(K[np.ix_(itr, itr)])
                Uy = U.T @ (y[itr] - mu_i)
                Kte = K[np.ix_(ite, itr)]
                for a_i, a in enumerate(alphas):
                    dual = U @ (Uy / (S + a))
                    p_ = mu_i + Kte @ dual
                    score[a_i] += -((y[ite] - p_) ** 2).sum()
        a_best = alphas[int(np.argmax(score))]
        chosen.append(float(a_best))
        inner_mode.append("donor" if used_donor_inner else "batch")
        mu = y[tr].mean()
        S, U = np.linalg.eigh(K[np.ix_(tr, tr)])
        dual = U @ ((U.T @ (y[tr] - mu)) / (S + a_best))
        pred[te] = mu + K[np.ix_(te, tr)] @ dual
    assert not np.isnan(pred).any()
    return pred, chosen, inner_mode


def _log_filter(rows, log, step, before, after, unit, note=""):
    removed = before - after
    rows.append(dict(step=step, unit=unit, before=int(before), after=int(after),
                     removed=int(removed), note=note))
    log(f"[filter] {step:<55} {unit:<10} {before:>10,} -> {after:>10,}  "
        f"(removed {removed:,}) {note}")


def run(name=DLPFC_NAME):
    rng = np.random.default_rng([BRAIN_SEED, _name_salt(name)])
    info = resolve(name)
    path = local_path(name)
    log = Logger(BRAIN_DIR / "r1_r3_analyze.txt")
    try:
        return _run(name, info, path, rng, log)
    finally:
        log.close()


def _run(name, info, path, rng, log):
    if not path.exists():
        raise FileNotFoundError(f"{path} missing — run tissue_download_top.py {name}")
    log("=" * 100)
    log(f"BRAIN ANALYZE  name={name}  seed={BRAIN_SEED}  dataset_id={info['dataset_id']}")
    log(f"  title={info['title']}  tissue={info['tissue']}")
    log(f"  local={path}  size={path.stat().st_size:,}")
    log(f"  adult cutoff age>={ADULT_MIN_AGE:.0f} (paper: adolescence 12–19 is developmental; "
        f"transcriptome stabilizes after 20).  FORCE_BATCH={FORCE_BATCH}")
    log("  NOTE: this is single-NUCLEUS RNA (10x 3' v3). Fewer transcripts/cell and different "
        "gene coverage than blood scRNA. Metrics are comparable; absolute sensitivity is not.")
    log("=" * 100)

    import anndata as ad
    a = ad.read_h5ad(str(path), backed="r")
    obs = a.obs.copy()
    obs["_pos"] = np.arange(len(obs))
    n0 = len(obs)
    log(f"[load] {n0:,} nuclei; obs cols={list(obs.columns)}")
    filt = []

    # is_primary_data is False for every nucleus (duplicated into the parent PsychAD collection).
    if "is_primary_data" in obs:
        n_pri = int(obs.is_primary_data.astype(bool).sum())
        log(f"[filter] is_primary_data==True would keep {n_pri:,}/{n0:,} "
            f"({100*n_pri/max(n0,1):.1f}%) — NOT applied (Aging_Cohort cells are marked "
            f"non-primary because they also sit in the parent cross-disorder collection)")
        _log_filter(filt, log, "is_primary_data (NOT applied; all False)", n0, n0, "cells",
                    "kept all; flag is collection-duplication not double-counting inside this file")

    # R1: neurotypical controls. Psychiatric/neurodegeneration flags are constant in this file
    # (it is already the Aging_Cohort neurotypical subset). Log them; drop only if a positive
    # brain-disease label appears.
    n = len(obs)
    brain_disease_pos = np.zeros(n, dtype=bool)
    for c in BRAIN_DISEASE_COLS:
        if c not in obs.columns:
            log(f"[R1] brain-disease column {c} not in obs")
            continue
        vc = obs[c].astype(str).value_counts(dropna=False)
        log(f"[R1] {c}: {list(vc.items())}")
        pos = ~obs[c].astype(str).str.lower().isin(
            {"nan", "none", "no", "false", "0", "control", "negative", "normal", "na", ""}
        )
        # constant "unknown"/single category that isn't an obvious negative: if nunique==1, keep all
        if obs[c].nunique(dropna=False) <= 1:
            log(f"[R1] {c} is constant — file is the neurotypical subset for this flag; not used to drop")
            continue
        brain_disease_pos |= pos.to_numpy()
    if brain_disease_pos.any():
        keep = ~brain_disease_pos
        _log_filter(filt, log, "drop nuclei with a brain-disease flag", n, int(keep.sum()), "cells")
        obs = obs[keep]
        n = len(obs)
    else:
        _log_filter(filt, log, "neurotypical (no brain-disease flag to drop)", n, n, "cells",
                    "Aging_Cohort is already the neurotypical subset of PsychAD; CXG `disease` is "
                    "systemic comorbidity (atherosclerosis/diabetes), kept")

    if "disease" in obs.columns:
        log(f"[R1] CXG disease (systemic, kept): {list(obs.disease.astype(str).value_counts().items())}")

    obs["age"] = obs["development_stage"].map(parse_age_years).astype(float)
    n = len(obs)
    obs = obs[obs.age.notna()]
    _log_filter(filt, log, "nuclei with numeric donor age", n, len(obs), "cells",
                "non-numeric development_stage dropped")
    n = len(obs)
    n_dev = int((obs.age < ADULT_MIN_AGE).sum())
    obs = obs[obs.age >= ADULT_MIN_AGE]
    _log_filter(filt, log, f"adult nuclei (age >= {ADULT_MIN_AGE:.0f})", n, len(obs), "cells",
                f"developmental/adolescent nuclei removed={n_dev:,}")

    obs["donor"] = obs["donor_id"].astype(str)
    obs["celltype"] = obs["cell_type"].astype(str)
    if "sex" in obs.columns:
        obs["sex"] = obs["sex"].astype(str)

    extra = [c for c in (
        "assay", "suspension_type", "tissue", "disease", "Source", "genetic_ancestry",
        "class", "subclass", "self_reported_ethnicity", *BRAIN_DISEASE_COLS,
    ) if c in obs.columns]

    def _mode(s):
        s = s.dropna().astype(str)
        if s.empty:
            return "NA"
        return s.value_counts().index[0]

    agg = dict(age=("age", "first"), n_cells=("age", "size"), age_nunique=("age", "nunique"))
    if "sex" in obs:
        agg["sex"] = ("sex", "first")
    for c in extra:
        agg[c] = (c, _mode)
    donors = obs.groupby("donor", observed=True).agg(**agg)
    n_bad_age = int((donors.age_nunique != 1).sum())
    if n_bad_age:
        log(f"[warn] {n_bad_age} donors have >1 age value; using first")
    span = float(donors.age.max() - donors.age.min())
    log(f"[R1] {len(donors)} adult control donors; age {donors.age.min():.0f}–{donors.age.max():.0f} "
        f"(nunique={donors.age.nunique()}, span={span:.0f} y); "
        f"nuclei/donor median {donors.n_cells.median():.0f}; nuclei total {int(donors.n_cells.sum()):,}")
    if "Source" in donors:
        log(f"[R1] Source (brain bank) donor counts: {donors.Source.astype(str).value_counts().to_dict()}")
        log("      " + donors.groupby(donors.Source.astype(str)).age.agg(["min", "max", "median", "count"]).to_string().replace("\n", "\n      "))

    r1 = dict(n_donors=int(len(donors)), n_nuclei=int(donors.n_cells.sum()),
              age_min=float(donors.age.min()), age_max=float(donors.age.max()),
              age_span=span, adult_min_age=ADULT_MIN_AGE, n_batches_source=int(donors.Source.nunique())
              if "Source" in donors else None)
    pd.DataFrame([r1]).to_csv(BRAIN_DIR / "r1_cohort.csv", index=False)
    donors.reset_index().to_csv(BRAIN_DIR / "r1_donors.csv", index=False)

    if len(donors) < STOP_R1_MIN_DONORS or span < STOP_R1_MIN_SPAN:
        log(f"[STOP R1] adult control donors={len(donors)} (need >={STOP_R1_MIN_DONORS}) "
            f"span={span:.0f} y (need >={STOP_R1_MIN_SPAN:.0f}). Halt — file adds little over "
            f"the 20-donor white matter; run R0–R3 on retina snRNA instead.")
        json.dump(dict(stop_r1=True, r1=r1, seed=BRAIN_SEED), open(BRAIN_DIR / "summary.json", "w"), indent=1)
        pd.DataFrame(filt).to_csv(BRAIN_DIR / "filter_log.csv", index=False)
        return dict(stop_r1=True, r1=r1)

    # ---- R2 batch column: force Source (R0 decision) ----
    if FORCE_BATCH not in donors.columns:
        raise RuntimeError(f"FORCE_BATCH {FORCE_BATCH} not in donor table")
    donors["batch"] = donors[FORCE_BATCH].astype(str)
    n_lev = donors.batch.nunique()
    sizes = donors.batch.value_counts()
    aliased = n_lev == len(donors)
    log(f"[R2] batch column FORCED to {FORCE_BATCH}: n_levels={n_lev}  "
        f"donors/level min={sizes.min()} max={sizes.max()} median={sizes.median():.1f}  "
        f"aliased_with_donor={aliased}")
    log("     No sequencing-pool / library / run column exists in obs (R0). Hashing pools of 6 "
        "described in the methods were not uploaded. Within-batch = within brain bank (H vs M).")
    if aliased or n_lev < 2 or int((sizes >= 2).sum()) < 2:
        raise RuntimeError("forced batch is not a shared-donor factor")
    batch_col = FORCE_BATCH
    batch_usable = True
    r2_age_batch = _safe_ols_r2(donors.age, donors.batch)
    log(f"[R2] R2(age ~ C(batch={batch_col})) = {r2_age_batch:.3f}   n_batches={n_lev}")

    keep_pos = obs["_pos"].to_numpy()
    keep_mask = np.zeros(n0, dtype=bool)
    keep_mask[keep_pos] = True
    min_don_ct = max(8, min(10, len(donors) // 3))
    try:
        a.file.close()
    except Exception:
        pass
    del a

    log(f"\n[pseudobulk] min_cells={MIN_CELLS}  min_donors_per_type={min_don_ct}")
    FILTER_LOG.clear()
    pb = build(
        str(path), f"tissue_{name}",
        donor_col="donor_id", ct_col="cell_type", age_col="development_stage",
        sex_col="sex" if "sex" in obs.columns else "sex",
        extra_obs=tuple(extra),
        min_cells=MIN_CELLS, min_donors=min_don_ct,
        require_primary=False,
        keep_mask=keep_mask,
    )
    for row in FILTER_LOG:
        filt.append(row)
        log(f"[filter] {row.get('step','')}")
    n_pb = len(pb)
    keep_d = set(donors.index)
    pb = pb[pb.obs.donor.isin(keep_d)].copy()
    _log_filter(filt, log, "pseudobulks of adult analysis donors", n_pb, pb.n_obs, "groups")
    ct_n = pb.obs.groupby("celltype").donor.nunique()
    keep_ct = ct_n[ct_n >= min_don_ct].index
    n1 = pb.n_obs
    pb = pb[pb.obs.celltype.isin(keep_ct)].copy()
    _log_filter(filt, log, f"cell types with >= {min_don_ct} adult donors", n1, pb.n_obs, "groups",
                f"types {len(ct_n)} -> {len(keep_ct)} dropped={dict(ct_n[ct_n < min_don_ct])}")

    Y, gid, sym = normalize_logcpm(pb, log=log)
    pb.obs["donor"] = pb.obs["donor"].astype(str)
    pb.obs["celltype"] = pb.obs["celltype"].astype(str)
    pb.obs["batch"] = pb.obs.donor.map(donors.batch)

    depth_d = pb.obs.groupby("donor").apply(
        lambda g: np.average(g.mean_counts_per_cell, weights=g.n_cells), include_groups=False
    )
    genes_d = pb.obs.groupby("donor").apply(
        lambda g: np.average(g.mean_genes_per_cell, weights=g.n_cells), include_groups=False
    )
    donors = donors.join(depth_d.rename("mean_counts_per_cell"), how="inner")
    donors = donors.join(genes_d.rename("mean_genes_per_cell"), how="left")
    _log_filter(filt, log, "donors with >=1 retained pseudobulk", len(keep_d), len(donors), "donors")
    counts, prop, clr, cats, Ncell = build_composition(obs, donors, log, min_donors_type=min_don_ct)
    donors = donors.join(Ncell.rename("n_cells_comp"), how="left")

    r2_depth_batch = _safe_ols_r2(np.log(donors.mean_counts_per_cell), donors.batch)
    log(f"[R2] R2(log mean UMI/nucleus ~ C(batch)) = {r2_depth_batch:.3f}")
    donors["age_b"] = donors.groupby("batch").age.transform("mean")
    donors["dep_b"] = donors.groupby("batch").mean_counts_per_cell.transform("mean")
    bm = donors.groupby("batch").agg(age=("age", "mean"), depth=("mean_counts_per_cell", "mean"))
    r_between = float(np.corrcoef(bm.age, bm.depth)[0, 1]) if len(bm) > 2 else np.nan
    r_within = float(np.corrcoef(donors.age - donors.age_b, donors.mean_counts_per_cell - donors.dep_b)[0, 1])
    r_pooled = float(np.corrcoef(donors.age, donors.mean_counts_per_cell)[0, 1])
    log(f"[R2] corr(age, depth): pooled r={r_pooled:+.3f}; between-batch r={r_between:+.3f} "
        f"(n={len(bm)} batches; NaN if only 2 batches); within-batch r={r_within:+.3f}")
    log("[R2] reference: blood R2(age~batch)=0.23  R2(depth~batch)=0.86  r_between/within=+0.43/−0.03")
    log("[R2] reference: white matter R2(age~batch)=0.272  R2(depth~batch)=0.626  r=+0.625/+0.006")

    t2 = dict(name=name, dataset_id=info["dataset_id"], tissue="brain (DLPFC)",
              n_donors=int(len(donors)), n_batches=int(donors.batch.nunique()),
              batch_column=batch_col, batch_usable=True,
              batch_kind="brain_bank_not_sequencing_pool",
              r2_age_batch=float(r2_age_batch), r2_depth_batch=float(r2_depth_batch),
              r_age_depth_pooled=float(r_pooled),
              r_age_depth_between=None if not np.isfinite(r_between) else float(r_between),
              r_age_depth_within=float(r_within),
              age_min=float(donors.age.min()), age_max=float(donors.age.max()),
              n_cell_types=int(pb.obs.celltype.nunique()), n_pseudobulks=int(pb.n_obs),
              mean_umi_per_nucleus=float(donors.mean_counts_per_cell.mean()),
              blood_r2_age_batch=0.23, blood_r2_depth_batch=0.86,
              wm_r2_age_batch=0.272, wm_r2_depth_batch=0.626)
    pd.DataFrame([t2]).to_csv(BRAIN_DIR / "r2_batch_audit.csv", index=False)
    donors.reset_index().to_csv(BRAIN_DIR / "r2_donors.csv", index=False)

    # ---- R3 turnover test (T3 / Stage B path) ----
    log("\n" + "=" * 100)
    log("[R3] batch-grouped CV: composition vs within-cell-type expression")
    log("=" * 100)
    n_batch = donors.batch.nunique()
    n_outer = int(min(N_OUTER, n_batch))
    n_rep = N_REPEATS if n_batch >= 5 else max(3, n_batch)
    log(f"  n_donors={len(donors)} n_batches={n_batch} n_outer={n_outer} n_repeats={n_rep}")
    if n_batch < 3:
        log("  [note] n_batches<3: outer CV is leave-one-site-out. Repeats are not independent "
            "(fold assignment is determined by the 2 sites). Inner alpha uses donor-grouped CV "
            "inside the training site because batch-grouped inner CV is degenerate. "
            "A literal tissue_analyze kernel_ridge_oof (accidental alpha=0.01) is also run as a check.")

    types = sorted(pb.obs.celltype.unique())
    donor_idx = pd.Index(donors.index)
    obs_pb = pb.obs
    ct = obs_pb.celltype.to_numpy()
    pool_pb = obs_pb.batch.to_numpy()
    donor_pb = obs_pb.donor.to_numpy()
    kern = {}
    log("building per-type kernels:")
    for t in types:
        rows_t = np.flatnonzero(ct == t)
        X = Y[rows_t]
        Xc = pool_center(X, pool_pb[rows_t])
        sd = Xc.std(0); sd[sd == 0] = 1; Xc = Xc / sd
        Xn = X - X.mean(0); sdn = Xn.std(0); sdn[sdn == 0] = 1; Xn = Xn / sdn
        kern[t] = (rows_t, Xc @ Xc.T / X.shape[1], Xn @ Xn.T / X.shape[1])
        log(f"   {t:<55} n={len(rows_t):>4} batches={len(np.unique(pool_pb[rows_t]))}")
    nD = len(donors)
    K_all_w = np.zeros((nD, nD)); K_all_n = np.zeros((nD, nD))
    for t in types:
        rows_t, Kw, Kn = kern[t]
        di = donor_idx.get_indexer(donor_pb[rows_t])
        K_all_w[np.ix_(di, di)] += Kw
        K_all_n[np.ix_(di, di)] += Kn
    K_all_w /= max(len(types), 1); K_all_n /= max(len(types), 1)
    pool_d = donors.batch.to_numpy()
    Cw = pool_center(clr.to_numpy(), pool_d)
    sd = Cw.std(0); sd[sd == 0] = 1; Cw = Cw / sd
    Cn = clr.to_numpy() - clr.to_numpy().mean(0); sdn = Cn.std(0); sdn[sdn == 0] = 1; Cn = Cn / sdn
    K_comp_w = Cw @ Cw.T / max(Cw.shape[1], 1)
    K_comp_n = Cn @ Cn.T / max(Cn.shape[1], 1)
    y_d = donors.age.to_numpy(float)
    donors["age_c"] = donors.age - donors.groupby("batch").age.transform("mean")
    yc_d = donors.age_c.to_numpy(float)
    donor_ids = donors.index.to_numpy()

    results = []
    oof_store = {}
    oof_literal = {}
    per_type_oof = {}
    for rep in range(n_rep):
        rrng = np.random.default_rng([BRAIN_SEED, _name_salt(name), rep])
        fold_of_pool = make_pool_folds(donors.batch.unique(), n_outer, rrng)
        models = (
            (f"composition (CLR, {clr.shape[1]} categories)", K_comp_w, K_comp_n, None),
            (f"expression: all {len(types)} cell types jointly", K_all_w, K_all_n, None),
        )
        for label, Kw_, Kn_, _ in models:
            pw, aw, mode = kernel_ridge_oof_with_donor_inner(
                Kw_, yc_d, pool_d, donor_ids, fold_of_pool, ALPHAS, rrng, n_inner=min(N_INNER, n_outer))
            pn, an, _ = kernel_ridge_oof_with_donor_inner(
                Kn_, y_d, pool_d, donor_ids, fold_of_pool, ALPHAS, rrng, n_inner=min(N_INNER, n_outer))
            r2w, maew = within_pool_r2_mae(yc_d, pw, pool_d)
            r2n, maen = raw_r2_mae(y_d, pn)
            results.append(dict(repeat=rep, model=label, celltype="(donor level)", n=nD,
                                r2_within_batch=r2w, mae_within=maew, r2_naive=r2n, mae_naive=maen,
                                alpha_mean=float(np.mean(aw)), inner_mode=";".join(mode),
                                path="donor_inner_fallback"))
            if rep == 0:
                oof_store[label] = (pw, pn)
            # literal tissue_analyze path (may pick alpha=0.01 when inner is degenerate)
            pw0, aw0 = kernel_ridge_oof(
                Kw_, yc_d, pool_d, donor_ids, fold_of_pool, ALPHAS, rrng, n_inner=min(N_INNER, n_outer))
            pn0, an0 = kernel_ridge_oof(
                Kn_, y_d, pool_d, donor_ids, fold_of_pool, ALPHAS, rrng, n_inner=min(N_INNER, n_outer))
            r2w0, maew0 = within_pool_r2_mae(yc_d, pw0, pool_d)
            r2n0, maen0 = raw_r2_mae(y_d, pn0)
            results.append(dict(repeat=rep, model=label, celltype="(donor level)", n=nD,
                                r2_within_batch=r2w0, mae_within=maew0, r2_naive=r2n0, mae_naive=maen0,
                                alpha_mean=float(np.mean(aw0)), inner_mode="literal_tissue_analyze",
                                path="literal_T3"))
            if rep == 0:
                oof_literal[label] = (pw0, pn0)
        # per-cell-type expression (primary path = donor-inner fallback)
        for t in types:
            rows_t, Kw_t, Kn_t = kern[t]
            y_t = yc_d[donor_idx.get_indexer(donor_pb[rows_t])]
            yraw_t = y_d[donor_idx.get_indexer(donor_pb[rows_t])]
            pool_t = pool_pb[rows_t]
            don_t = donor_pb[rows_t]
            if len(np.unique(pool_t)) < 2:
                continue
            rrng_t = np.random.default_rng([BRAIN_SEED, _name_salt(name), rep, _name_salt(t)])
            pw, aw, mode = kernel_ridge_oof_with_donor_inner(
                Kw_t, y_t, pool_t, don_t, fold_of_pool, ALPHAS, rrng_t, n_inner=min(N_INNER, n_outer))
            pn, an, _ = kernel_ridge_oof_with_donor_inner(
                Kn_t, yraw_t, pool_t, don_t, fold_of_pool, ALPHAS, rrng_t, n_inner=min(N_INNER, n_outer))
            r2w, maew = within_pool_r2_mae(y_t, pw, pool_t)
            r2n, maen = raw_r2_mae(yraw_t, pn)
            results.append(dict(repeat=rep, model="expression: within cell type", celltype=t, n=len(rows_t),
                                r2_within_batch=r2w, mae_within=maew, r2_naive=r2n, mae_naive=maen,
                                alpha_mean=float(np.mean(aw)), inner_mode=";".join(mode),
                                path="donor_inner_fallback"))
            if rep == 0:
                per_type_oof[t] = (pw, pn, y_t, yraw_t, pool_t)
        rec = [r for r in results if r["repeat"] == rep and r["path"] == "donor_inner_fallback"
               and r["celltype"] == "(donor level)"]
        log(f"   repeat {rep}: comp R2w={rec[0]['r2_within_batch']:+.3f}  "
            f"expr R2w={rec[1]['r2_within_batch']:+.3f}  "
            f"comp naive={rec[0]['r2_naive']:+.3f}  expr naive={rec[1]['r2_naive']:+.3f}  "
            f"alpha_comp={rec[0]['alpha_mean']:.3g} alpha_expr={rec[1]['alpha_mean']:.3g}")

    res = pd.DataFrame(results)
    res.to_csv(BRAIN_DIR / "r3_cv_repeats.csv", index=False)

    def _summ(sub, pw, pn, y, yc, pool):
        lo, hi = bootstrap_ci(yc, pw, pool, np.random.default_rng([BRAIN_SEED, 99]), within=True, n_boot=N_BOOT)
        lon, hin = bootstrap_ci(y, pn, pool, np.random.default_rng([BRAIN_SEED, 98]), within=False, n_boot=N_BOOT)
        return dict(r2_within_mean=float(sub.r2_within_batch.mean()),
                    r2_within_min=float(sub.r2_within_batch.min()),
                    r2_within_max=float(sub.r2_within_batch.max()),
                    r2_within_lo=float(lo), r2_within_hi=float(hi),
                    r2_naive_mean=float(sub.r2_naive.mean()),
                    r2_naive_lo=float(lon), r2_naive_hi=float(hin),
                    mae_within=float(sub.mae_within.mean()), mae_naive=float(sub.mae_naive.mean()))

    prim = res[(res.path == "donor_inner_fallback") & (res.celltype == "(donor level)")]
    labels = list(dict.fromkeys(prim.model))
    summ_rows = []
    for lab in labels:
        sub = prim[prim.model == lab]
        pw, pn = oof_store[lab]
        s = _summ(sub, pw, pn, y_d, yc_d, pool_d)
        s.update(model=lab, name=name, path="donor_inner_fallback", celltype="(donor level)")
        summ_rows.append(s)
        log(f"  {lab}: within-batch R2={s['r2_within_mean']:+.3f} "
            f"[{s['r2_within_lo']:+.3f}, {s['r2_within_hi']:+.3f}] "
            f"{{{s['r2_within_min']:+.3f},{s['r2_within_max']:+.3f}}}  "
            f"naive R2={s['r2_naive_mean']:+.3f}")
    lit = res[(res.path == "literal_T3") & (res.celltype == "(donor level)")]
    for lab in labels:
        sub = lit[lit.model == lab]
        pw, pn = oof_literal[lab]
        s = _summ(sub, pw, pn, y_d, yc_d, pool_d)
        s.update(model=lab, name=name, path="literal_T3", celltype="(donor level)")
        summ_rows.append(s)
        log(f"  [literal T3] {lab}: within-batch R2={s['r2_within_mean']:+.3f} "
            f"{{{s['r2_within_min']:+.3f},{s['r2_within_max']:+.3f}}}")

    per_t_rows = []
    per_t = res[(res.path == "donor_inner_fallback") & (res.model == "expression: within cell type")]
    for t, sub in per_t.groupby("celltype"):
        pw, pn, yt, yraw, pt = per_type_oof[t]
        s = _summ(sub, pw, pn, yraw, yt, pt)
        s.update(model="expression: within cell type", name=name, path="donor_inner_fallback", celltype=t, n=int(sub.n.iloc[0]))
        per_t_rows.append(s)
        log(f"  per-type {t}: within-batch R2={s['r2_within_mean']:+.3f} "
            f"{{{s['r2_within_min']:+.3f},{s['r2_within_max']:+.3f}}}")
    summ = pd.concat([pd.DataFrame(summ_rows), pd.DataFrame(per_t_rows)], ignore_index=True)
    summ.to_csv(BRAIN_DIR / "r3_cv_summary.csv", index=False)

    comp = [s for s in summ_rows if s["path"] == "donor_inner_fallback" and s["model"].startswith("composition")][0]
    expr = [s for s in summ_rows if s["path"] == "donor_inner_fallback" and s["model"].startswith("expression")][0]
    c_r2, e_r2 = comp["r2_within_mean"], expr["r2_within_mean"]
    delta = e_r2 - c_r2
    prim_comp = prim[prim.model.str.startswith("composition")].sort_values("repeat")
    prim_expr = prim[prim.model.str.startswith("expression")].sort_values("repeat")
    n_expr_gt = int((prim_expr.r2_within_batch.to_numpy() > prim_comp.r2_within_batch.to_numpy()).sum())
    log(f"\n[R3] DLPFC: composition R2={c_r2:+.3f}  expression R2={e_r2:+.3f}  expr−comp={delta:+.3f}")
    log(f"     expression > composition in {n_expr_gt}/{len(prim_comp)} repeats")
    log(f"     blood: composition {BLOOD_COMP_R2:.3f}  expression {BLOOD_EXPR_R2:.3f}")
    log(f"     white matter: composition {WM_COMP_R2:.3f}  expression {WM_EXPR_R2:.3f}")
    if len(per_t_rows):
        pt = pd.DataFrame(per_t_rows).sort_values("r2_within_mean")
        log(f"     per-type expression R2: median {pt.r2_within_mean.median():+.3f}  "
            f"range {pt.r2_within_mean.min():+.3f} .. {pt.r2_within_mean.max():+.3f}  "
            f"best={pt.celltype.iloc[-1]}")
        n_pos = int((pt.r2_within_mean > 0).sum())
        log(f"     per-type: {n_pos}/{len(pt)} types have mean within-batch R2 > 0")

    t3 = dict(name=name, dataset_id=info["dataset_id"], tissue="brain (DLPFC)",
              n_donors=int(nD), n_cell_types=len(types), n_comp_categories=int(clr.shape[1]),
              n_batches=int(n_batch), batch_column=batch_col, metric="within_batch",
              composition_r2=float(c_r2), expression_r2=float(e_r2),
              delta_expr_minus_comp=float(delta),
              n_repeats_expr_gt_comp=n_expr_gt, n_repeats=int(len(prim_comp)),
              composition_r2_lo=comp["r2_within_lo"], composition_r2_hi=comp["r2_within_hi"],
              expression_r2_lo=expr["r2_within_lo"], expression_r2_hi=expr["r2_within_hi"],
              composition_r2_min=comp["r2_within_min"], composition_r2_max=comp["r2_within_max"],
              expression_r2_min=expr["r2_within_min"], expression_r2_max=expr["r2_within_max"],
              composition_r2_naive=comp["r2_naive_mean"], expression_r2_naive=expr["r2_naive_mean"],
              blood_composition_r2=BLOOD_COMP_R2, blood_expression_r2=BLOOD_EXPR_R2,
              wm_composition_r2=WM_COMP_R2, wm_expression_r2=WM_EXPR_R2,
              seed=BRAIN_SEED)
    pd.DataFrame([t3]).to_csv(BRAIN_DIR / "r3_summary.csv", index=False)
    pd.DataFrame(filt).to_csv(BRAIN_DIR / "filter_log.csv", index=False)

    # figures
    fig, ax = plt.subplots(figsize=(7.6, 4.4))
    labs = ["blood\ncomp", "blood\nexpr", "WM\ncomp", "WM\nexpr", "DLPFC\ncomp", "DLPFC\nexpr"]
    vals = [BLOOD_COMP_R2, BLOOD_EXPR_R2, WM_COMP_R2, WM_EXPR_R2, c_r2, e_r2]
    cols = ["0.75", "0.55", "tab:orange", "tab:blue", "darkorange", "royalblue"]
    ax.bar(range(6), vals, color=cols)
    ax.set_xticks(range(6)); ax.set_xticklabels(labs, fontsize=8)
    ax.set_ylabel("held-out within-batch R²")
    ax.axhline(0, c="k", lw=0.5)
    ax.set_title(f"DLPFC adult controls (n={nD}): expr−comp = {delta:+.3f}")
    fig.tight_layout(); fig.savefig(BRAIN_FIG / "r3_comp_vs_expr.png", dpi=140); plt.close(fig)

    fig, ax = plt.subplots(figsize=(6.2, 3.4))
    ax.hist(donors.age, bins=20, color="steelblue")
    ax.set_xlabel("donor age (years)"); ax.set_ylabel("n donors")
    ax.set_title(f"DLPFC adult controls: {len(donors)} donors, {donors.age.min():.0f}–{donors.age.max():.0f} y")
    fig.tight_layout(); fig.savefig(BRAIN_FIG / "r1_age_hist.png", dpi=130); plt.close(fig)

    if "Source" in donors.columns:
        fig, ax = plt.subplots(figsize=(6.2, 3.6))
        for src, g in donors.groupby(donors.Source.astype(str)):
            ax.hist(g.age, bins=15, alpha=0.55, label=f"{src} n={len(g)}")
        ax.set_xlabel("donor age (years)"); ax.set_ylabel("n donors"); ax.legend(fontsize=8)
        ax.set_title("Age by Source (H=HBCC, M=MSSM)")
        fig.tight_layout(); fig.savefig(BRAIN_FIG / "r2_age_by_source.png", dpi=130); plt.close(fig)

    if per_t_rows:
        pt = pd.DataFrame(per_t_rows).sort_values("r2_within_mean")
        fig, ax = plt.subplots(figsize=(8.0, max(3.5, 0.28 * len(pt))))
        ax.barh(pt.celltype, pt.r2_within_mean, color="tab:blue")
        ax.axvline(0, c="k", lw=0.5)
        ax.set_xlabel("held-out within-batch R² (expression)")
        ax.set_title("Per-cell-type age prediction (DLPFC adult controls)")
        fig.tight_layout(); fig.savefig(BRAIN_FIG / "r3_per_type.png", dpi=140); plt.close(fig)

    json.dump(dict(r1=r1, t2=t2, t3=t3, types=types, composition_categories=cats,
                   seed=BRAIN_SEED, adult_min_age=ADULT_MIN_AGE, force_batch=FORCE_BATCH),
              open(BRAIN_DIR / "summary.json", "w"), indent=1)

    cmp = pd.DataFrame([
        dict(tissue="blood PBMC (OneK1K, closed)", n_donors=981, n_batches=75, batch="pool (75 10x)",
             composition_r2=BLOOD_COMP_R2, expression_r2=BLOOD_EXPR_R2,
             expr_minus_comp=BLOOD_EXPR_R2 - BLOOD_COMP_R2, repeats_expr_gt_comp="n/a (Stage B mean)",
             data="scRNA"),
        dict(tissue="brain white matter (c05e6940)", n_donors=20, n_batches=7, batch="SequencingPool",
             composition_r2=WM_COMP_R2, expression_r2=WM_EXPR_R2,
             expr_minus_comp=WM_EXPR_R2 - WM_COMP_R2, repeats_expr_gt_comp="10/10", data="snRNA"),
        dict(tissue="brain DLPFC Aging_Cohort (adult controls)", n_donors=int(nD),
             n_batches=int(n_batch), batch="Source (H/M brain bank)",
             composition_r2=float(c_r2), expression_r2=float(e_r2),
             expr_minus_comp=float(delta),
             repeats_expr_gt_comp=f"{n_expr_gt}/{len(prim_comp)}", data="snRNA"),
    ])
    cmp.to_csv(BRAIN_DIR / "r4_comparison.csv", index=False)
    if per_t_rows:
        pd.DataFrame(per_t_rows)[["celltype", "n", "r2_within_mean", "r2_within_min", "r2_within_max"]].sort_values(
            "r2_within_mean", ascending=False
        ).to_csv(BRAIN_DIR / "r3_per_type.csv", index=False)

    log(f"\n[brain_aging] done. wrote {BRAIN_DIR}")
    return dict(r1=r1, t2=t2, t3=t3, stop_r1=False)


if __name__ == "__main__":
    run()
