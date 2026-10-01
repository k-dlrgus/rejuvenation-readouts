"""Shared helpers for the leakage audit / gene-count control / geometric reframe.

Namespaced under results/geometry/. Does not modify prior FINDINGS*.md or FALSIFICATION.md.
Does not derive new A/I gene-set definitions; uses the frozen P2 sets and the declared
P1 quantile rule (nested inside folds only as a leakage test).
Seed 20260914.
"""
from __future__ import annotations

import json
import sys
from pathlib import Path

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parent))
from config import RESULTS  # noqa: E402
from brain_common import DLPFC_DATASET_ID  # noqa: E402
from brain_phase1_common import (  # noqa: E402
    PHASE1_DIR, PHASE1_SEED, HASH_POOL_LIMITATION, ALPHAS, N_SITE_FOLDS,
    ensure_logcpm, merge_p0_into_obs, load_p0_extras, load_json, dump_json,
    loso_folds, within_site_donor_folds, r2_mae, Logger,
)
from brain_phase1_models import (  # noqa: E402
    gene_index, age_oof_for_genes, permute_age_within_site_obs, permute_celltype,
    select_genes_on_train, load_chrom_series, ridge_predict,
)
from brain_phase1_p3 import _scheme_folds, _assert_no_donor_leak  # noqa: E402

GEO_DIR = RESULTS / "geometry"
GEO_FIG = GEO_DIR / "figures"
GEO_TAB = GEO_DIR
for _p in (GEO_DIR, GEO_FIG):
    _p.mkdir(parents=True, exist_ok=True)

GEO_SEED = PHASE1_SEED  # 20260914
N_PERM_A = 20
N_PERM_C = 100
N_RANDOM_DRAWS = 50
N_SWEEP_DRAWS = 20
N_MATCH_DRAWS = 20
SIZE_SWEEP = (10, 50, 100, 500, 1000, 1500)
NULL_USABLE = 0.05  # STOP A if corrected shuffle R² still above this
DATASET_ID = DLPFC_DATASET_ID  # fetched in T1/R0; not invented


def geo_log_banner(log, part: str):
    log("=" * 100)
    log(f"GEOMETRY {part}  seed={GEO_SEED}  dataset_id={DATASET_ID}")
    log("=" * 100)
    log(HASH_POOL_LIMITATION)
    log("FALSIFICATION.md and prior FINDINGS*.md are not modified.")


def load_phase1_matrix(log=print):
    """Load the cached Phase-1 log-CPM matrix and frozen A/I sets. Log every filter count."""
    Y, genes, obs = ensure_logcpm(log)
    obs = merge_p0_into_obs(obs)
    n0 = len(obs)
    n_don = obs.donor.nunique()
    n_ct = obs.celltype.nunique()
    site_counts = obs.groupby("Source").donor.nunique().to_dict()
    log(f"[filter] log-CPM matrix: {Y.shape[0]:,} pseudobulks x {Y.shape[1]:,} genes")
    log(f"[filter] donors={n_don}  cell_types={n_ct}  sites={site_counts}  rows={n0}")
    frozen = load_json(PHASE1_DIR / "p2_frozen_gene_sets.json")
    A_idx = gene_index(genes, frozen["A_genes"])
    I_idx = gene_index(genes, frozen["I_genes"])
    log(f"[filter] frozen A-genes: requested={frozen['n_A']} mapped={len(A_idx)}")
    log(f"[filter] frozen I-genes: requested={frozen['n_I']} mapped={len(I_idx)}")
    if len(A_idx) != int(frozen["n_A"]) or len(I_idx) != int(frozen["n_I"]):
        log("[filter] WARNING: frozen-set mapping size mismatch")
    extra = load_p0_extras()
    log(f"[filter] P0 extras on obs: {extra}")
    return dict(Y=Y, genes=genes, obs=obs, frozen=frozen, A_idx=A_idx, I_idx=I_idx,
                extra=extra, n_donors=n_don, n_types=n_ct, site_donor_counts=site_counts)


def make_folds(obs, rng):
    return _scheme_folds(obs, rng)


def assert_donor_folds(obs, folds, log, name):
    _assert_no_donor_leak(obs, folds, log, name)
    leaked = []
    for i, fold in enumerate(folds):
        tr, te = fold[0], fold[1]
        ov = set(obs.iloc[tr].donor.astype(str)) & set(obs.iloc[te].donor.astype(str))
        if ov:
            leaked.append((i, sorted(ov)[:8]))
    return dict(n_folds=len(folds), n_leaked_folds=len(leaked), examples=leaked)


def permute_age_per_pseudobulk(obs, rng):
    """WRONG shuffle: permute the age column at row level (breaks donor-constant age)."""
    out = obs.copy()
    out["age"] = rng.permutation(out.age.to_numpy(float))
    return out


def permute_age_within_donor_rows(obs, rng):
    """WRONG shuffle if ages were heterogeneous within donor; a no-op if age is donor-constant."""
    out = obs.copy()
    age = out.age.to_numpy(float).copy()
    donors = out.donor.astype(str).to_numpy()
    for d in pd.unique(donors):
        idx = np.flatnonzero(donors == d)
        age[idx] = rng.permutation(age[idx])
    out["age"] = age
    return out


def donor_age_is_constant(obs) -> bool:
    g = obs.groupby("donor")["age"].nunique()
    return bool((g == 1).all())


def intercept_only_oof(obs, folds):
    """Per-cell-type training-mean predictor (the intercept a ridge always has)."""
    pred = np.full(len(obs), np.nan)
    y = obs.age.to_numpy(float)
    ct = obs.celltype.astype(str).to_numpy()
    types = sorted(pd.unique(ct))
    for fold in folds:
        tr, te = fold[0], fold[1]
        for t in types:
            tr_t = tr[ct[tr] == t]
            te_t = te[ct[te] == t]
            if len(tr_t) < 1 or len(te_t) < 1:
                continue
            pred[te_t] = float(np.mean(y[tr_t]))
    return pred


def ridge_svd_predict(Xtr, ytr, Xte, alphas=ALPHAS):
    """Ridge with closed-form LOO alpha selection (PRESS). Training intercept only."""
    Xtr = np.asarray(Xtr, np.float64)
    Xte = np.asarray(Xte, np.float64)
    ytr = np.asarray(ytr, np.float64)
    n, p = Xtr.shape
    if n < 5 or p == 0:
        return np.full(len(Xte), np.nan), np.nan
    mu = Xtr.mean(0)
    sd = Xtr.std(0)
    sd[sd < 1e-12] = 1.0
    Xtr_s = (Xtr - mu) / sd
    Xte_s = (Xte - mu) / sd
    ymu = float(ytr.mean())
    yc = ytr - ymu
    U, S, Vt = np.linalg.svd(Xtr_s, full_matrices=False)
    keep = S > 1e-10
    if not np.any(keep):
        return np.full(len(Xte), ymu), np.nan
    U, S, Vt = U[:, keep], S[keep], Vt[keep]
    Uy = U.T @ yc
    best_a, best_sse = float(alphas[0]), np.inf
    for a in alphas:
        lam = (S ** 2) / (S ** 2 + float(a))
        h = (U ** 2) @ lam
        fit = U @ (lam * Uy)
        denom = 1.0 - h
        denom[np.abs(denom) < 1e-12] = 1e-12
        sse = float(np.sum(((yc - fit) / denom) ** 2))
        if sse < best_sse:
            best_sse = sse
            best_a = float(a)
    d = S / (S ** 2 + best_a)
    coef = Vt.T @ (d * Uy)
    pred = ymu + Xte_s @ coef
    return pred, best_a


def ridge_svd_fit(X, y, alphas=ALPHAS):
    """Return (coef_on_standardized_X, intercept, mu, sd, alpha)."""
    X = np.asarray(X, np.float64)
    y = np.asarray(y, np.float64)
    n, p = X.shape
    if n < 5 or p == 0:
        return np.zeros(p), float(y.mean()) if n else np.nan, np.zeros(p), np.ones(p), np.nan
    mu = X.mean(0)
    sd = X.std(0)
    sd[sd < 1e-12] = 1.0
    Xs = (X - mu) / sd
    ymu = float(y.mean())
    yc = y - ymu
    U, S, Vt = np.linalg.svd(Xs, full_matrices=False)
    keep = S > 1e-10
    if not np.any(keep):
        return np.zeros(p), ymu, mu, sd, np.nan
    U, S, Vt = U[:, keep], S[keep], Vt[keep]
    Uy = U.T @ yc
    best_a, best_sse = float(alphas[0]), np.inf
    for a in alphas:
        lam = (S ** 2) / (S ** 2 + float(a))
        h = (U ** 2) @ lam
        fit = U @ (lam * Uy)
        denom = 1.0 - h
        denom[np.abs(denom) < 1e-12] = 1e-12
        sse = float(np.sum(((yc - fit) / denom) ** 2))
        if sse < best_sse:
            best_sse = sse
            best_a = float(a)
    d = S / (S ** 2 + best_a)
    coef = Vt.T @ (d * Uy)
    return coef, ymu, mu, sd, best_a


def age_oof_svd(Y, obs, gene_idx, folds, y_override=None):
    """Within-type OOF age predictions with SVD-LOO ridge. y_override lets us score shuffled ages."""
    pred = np.full(len(obs), np.nan)
    y = obs.age.to_numpy(float) if y_override is None else np.asarray(y_override, float)
    ct = obs.celltype.astype(str).to_numpy()
    types = sorted(pd.unique(ct))
    gi = gene_idx
    alphas = []
    for fi, fold in enumerate(folds):
        tr, te = fold[0], fold[1]
        g = gi[fi] if isinstance(gi, (list, tuple)) else gi
        g = np.asarray(g, dtype=int)
        for t in types:
            tr_t = tr[ct[tr] == t]
            te_t = te[ct[te] == t]
            if len(tr_t) < 8 or len(te_t) < 1:
                continue
            Xtr = Y[np.ix_(tr_t, g)] if g.size else np.zeros((len(tr_t), 0))
            Xte = Y[np.ix_(te_t, g)] if g.size else np.zeros((len(te_t), 0))
            p, a = ridge_svd_predict(Xtr, y[tr_t], Xte)
            pred[te_t] = p
            alphas.append(a)
    return dict(pred=pred, alpha_mean=float(np.nanmean(alphas) if alphas else np.nan), y=y)


def _finite_mask(y, pred):
    y = np.asarray(y, float)
    pred = np.asarray(pred, float)
    return np.isfinite(y) & np.isfinite(pred)


def score_age_predictions(obs, pred, y=None):
    """Both aggregations: original pooled-across-sites, and site-stratified.

    Original P3/P4 median-over-types R² pools a cell type's OOF rows from HBCC and MSSM
    before computing R². Within-site CV trains a separate intercept in each bank, so that
    pooled R² credits recovery of the between-bank age gap.

    Site-stratified R² computes R² inside each (cell type, site) cell, averages the two
    banks (equal weight), then takes the median over cell types. Intercept-only is ~0.
    """
    y = obs.age.to_numpy(float) if y is None else np.asarray(y, float)
    pred = np.asarray(pred, float)
    ct = obs.celltype.astype(str).to_numpy()
    site = obs.Source.astype(str).to_numpy()
    donors = obs.donor.astype(str).to_numpy()
    types = sorted(pd.unique(ct))
    sites = sorted(pd.unique(site))
    rows = []
    type_pooled, type_strat = [], []
    for t in types:
        m = (ct == t) & _finite_mask(y, pred)
        met_p = r2_mae(y[m], pred[m])
        rows.append(dict(celltype=t, site="pooled", **met_p))
        if np.isfinite(met_p["r2"]):
            type_pooled.append(met_p["r2"])
        site_r2, site_mae = [], []
        for s in sites:
            ms = m & (site == s)
            met_s = r2_mae(y[ms], pred[ms])
            rows.append(dict(celltype=t, site=s, **met_s))
            if np.isfinite(met_s["r2"]):
                site_r2.append(met_s["r2"])
                site_mae.append(met_s["mae"])
        if site_r2:
            type_strat.append(float(np.mean(site_r2)))
            rows.append(dict(celltype=t, site="mean_of_sites",
                             r2=float(np.mean(site_r2)),
                             mae=float(np.mean(site_mae)) if site_mae else np.nan,
                             n=int(m.sum())))
    # donor-level joint
    tmp = pd.DataFrame({"donor": donors, "y": y, "pred": pred, "site": site})
    tmp = tmp[_finite_mask(tmp.y, tmp.pred)]
    don = tmp.groupby("donor", sort=False).agg(y=("y", "first"), pred=("pred", "mean"),
                                               site=("site", "first"))
    joint_pooled = r2_mae(don.y, don.pred)
    site_joint = []
    for s in sites:
        sub = don[don.site == s]
        site_joint.append(r2_mae(sub.y, sub.pred)["r2"])
    per_type = pd.DataFrame(rows)
    cells = per_type[per_type.site.isin(sites)]
    return dict(
        pooled_median_type_r2=_nanmedian(type_pooled),
        site_strat_median_type_r2=_nanmedian(type_strat),
        site_cell_median_r2=float(cells.r2.median()) if len(cells) else np.nan,
        joint_pooled_r2=joint_pooled["r2"],
        joint_site_strat_r2=_nanmean(site_joint),
        pooled_median_type_mae=_nanmedian([r["mae"] for r in rows if r["site"] == "pooled"]),
        n_types_pooled=int(len(type_pooled)),
        n_types_strat=int(len(type_strat)),
        per_type=per_type,
        donor_pred=don,
    )


def _nanmedian(xs):
    xs = np.asarray(list(xs), float)
    xs = xs[np.isfinite(xs)]
    return float(np.median(xs)) if len(xs) else np.nan


def _nanmean(xs):
    xs = np.asarray(list(xs), float)
    xs = xs[np.isfinite(xs)]
    return float(np.mean(xs)) if len(xs) else np.nan


def primary_r2(scheme, scores):
    """The number that should be interpreted for this CV scheme."""
    if scheme == "within_site":
        return scores["site_strat_median_type_r2"]
    return scores["pooled_median_type_r2"]


def r2_age_on_site(obs):
    """How much donor-age variance is between the two banks."""
    don = obs.groupby("donor").agg(age=("age", "first"), site=("Source", "first"))
    y = don.age.to_numpy(float)
    # one-hot site (drop first)
    s = (don.site.astype(str) == sorted(don.site.astype(str).unique())[-1]).to_numpy(float)
    X = np.column_stack([np.ones(len(y)), s])
    pred = X @ np.linalg.lstsq(X, y, rcond=None)[0]
    return r2_mae(y, pred)["r2"]


def centering_audit(obs, pred, y=None):
    """Does R² change if ss_tot uses the full-cohort mean vs the evaluation-set mean?"""
    y = obs.age.to_numpy(float) if y is None else np.asarray(y, float)
    pred = np.asarray(pred, float)
    m = _finite_mask(y, pred)
    yv, pv = y[m], pred[m]
    ss_res = float(np.sum((yv - pv) ** 2))
    ss_tot_eval = float(np.sum((yv - yv.mean()) ** 2))
    ss_tot_full = float(np.sum((yv - y.mean()) ** 2))  # full-cohort mean, including rows without pred
    ss_tot_all_eval_fullmean = float(np.sum((yv - y[m].mean()) ** 2))  # same as eval if m is all scored
    r2_eval = np.nan if ss_tot_eval <= 0 else 1.0 - ss_res / ss_tot_eval
    r2_fullmean = np.nan if ss_tot_full <= 0 else 1.0 - ss_res / ss_tot_full
    return dict(r2_eval_mean=r2_eval, r2_full_cohort_mean_ss_tot=r2_fullmean,
                ss_tot_eval=ss_tot_eval, ss_tot_full=ss_tot_full, n=int(m.sum()))


def load_p3_oof(scheme):
    path = PHASE1_DIR / f"p3_oof_age_{scheme}.csv"
    return pd.read_csv(path)


def fmt(x, d=3):
    try:
        x = float(x)
    except (TypeError, ValueError):
        return "NA"
    if not np.isfinite(x):
        return "NA"
    return f"{x:+.{d}f}" if d else str(x)
