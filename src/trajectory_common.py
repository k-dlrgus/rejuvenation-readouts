"""Geometric age / identity scores — no gene selection.

Every gene in the Phase-1 log-CPM matrix enters every model with a weight.
Namespaced under results/trajectory/. Does not modify prior FINDINGS*.md or FALSIFICATION.md.
Seed 20260914 (same fold construction as Phase 1 / geometry).
"""
from __future__ import annotations

import sys
from pathlib import Path

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parent))
from config import RESULTS  # noqa: E402
from geometry_common import (  # noqa: E402
    GEO_SEED, DATASET_ID, NULL_USABLE, ALPHAS,
    Logger, dump_json, load_json, load_phase1_matrix, make_folds,
    permute_age_within_site_obs, score_age_predictions, primary_r2,
    ridge_svd_fit, HASH_POOL_LIMITATION,
)
from geometry_partC import (  # noqa: E402
    identity_basis, pls1_direction, proj_r2, angle_to_subspace_deg, _zscore,
)
from brain_phase1_models import permute_celltype  # noqa: E402
from brain_phase1_common import r2_mae  # noqa: E402

TRAJ_DIR = RESULTS / "trajectory"
TRAJ_FIG = TRAJ_DIR / "figures"
TRAJ_TAB = TRAJ_DIR
for _p in (TRAJ_DIR, TRAJ_FIG):
    _p.mkdir(parents=True, exist_ok=True)

TRAJ_SEED = GEO_SEED  # 20260914 — folds identical to geometry / Phase 1
N_PERM_G1 = 20
N_PERM_G2 = 50
N_PERM_ID = 20
# STOP G1: corrected (site-stratified within-site) shuffle R² still above this,
# OR observed age R² near zero under BOTH CV schemes.
STOP_NULL = NULL_USABLE  # 0.05
STOP_R2_COLLAPSE = 0.05
# G2b redundancy (mirrors the spirit of the gene-set cross-test, not a rewrite of it)
G2B_RATIO = 0.80
G2B_DELTA = 0.10
# G2c collapse: residual age R² loses most of the age-direction signal
G2C_RATIO = 0.50


def traj_log_banner(log, stage: str):
    log("=" * 100)
    log(f"TRAJECTORY {stage}  seed={TRAJ_SEED}  dataset_id={DATASET_ID}")
    log("=" * 100)
    log(HASH_POOL_LIMITATION)
    log("NO GENE SELECTION: every gene in the Phase-1 matrix enters every model with a weight.")
    log("FALSIFICATION.md and prior FINDINGS*.md are not modified.")


def all_gene_idx(n_genes: int) -> np.ndarray:
    return np.arange(int(n_genes), dtype=int)


def zscore_train(Xtr, Xte=None):
    """Train-only mean/sd. Xte optional."""
    Xtr = np.asarray(Xtr, np.float64)
    mu = Xtr.mean(0)
    sd = Xtr.std(0)
    sd = np.where(sd < 1e-12, 1.0, sd)
    Xtr_s = (Xtr - mu) / sd
    if Xte is None:
        return Xtr_s, mu, sd
    Xte = np.asarray(Xte, np.float64)
    return Xtr_s, (Xte - mu) / sd, mu, sd


def calibrate_1d(s_tr, y_tr, s_te):
    """OLS y ~ 1 + score on train; apply to test. Returns pred_te, intercept, slope."""
    s_tr = np.asarray(s_tr, np.float64)
    y_tr = np.asarray(y_tr, np.float64)
    s_te = np.asarray(s_te, np.float64)
    if len(y_tr) < 3 or np.nanstd(s_tr) < 1e-12:
        mu = float(np.nanmean(y_tr)) if len(y_tr) else np.nan
        return np.full(len(s_te), mu), mu, 0.0
    A = np.column_stack([np.ones(len(s_tr)), s_tr])
    coef, *_ = np.linalg.lstsq(A, y_tr, rcond=None)
    pred = coef[0] + coef[1] * s_te
    return pred, float(coef[0]), float(coef[1])


def pls1_predict(Xtr, ytr, Xte):
    """Projection onto the PLS-1 covariance direction, then 1D OLS to years.

    Xtr/Xte are already standardized on train. Every gene gets weight X'y.
    """
    w, nrm = pls1_direction(Xtr, ytr)
    if nrm < 1e-12 or not np.isfinite(w).any():
        mu = float(np.mean(ytr)) if len(ytr) else np.nan
        return np.full(len(Xte), mu), w, np.full(len(Xte), 0.0), (mu, 0.0)
    s_tr = Xtr @ w
    s_te = Xte @ w
    pred, a, b = calibrate_1d(s_tr, ytr, s_te)
    return pred, w, s_te, (a, b)


def ridge_predict_dir(Xtr, ytr, Xte, alphas=ALPHAS):
    """Ridge on all genes (SVD-LOO alpha). Score = X_s @ coef (years, with intercept)."""
    coef, ymu, mu, sd, a = ridge_svd_fit(Xtr, ytr, alphas=alphas)
    # Xtr/Xte already standardized: mu~0, sd~1. Re-standardizing is a no-op if already z-scored
    # but ridge_svd_fit re-z-scores internally. Apply the same transform to Xte.
    Xte = np.asarray(Xte, np.float64)
    Xte_s = (Xte - mu) / sd
    pred = ymu + Xte_s @ coef
    nrm = float(np.linalg.norm(coef))
    w = coef / nrm if nrm > 1e-12 else coef
    return pred, w, pred - ymu, (float(ymu), a, nrm)


def ridge_from_svd_cache(U, S, Vt, ytr, Xte_s, alphas=ALPHAS):
    """Reuse train SVD when only y changes (permutation null). Xte_s already train-standardized."""
    ytr = np.asarray(ytr, np.float64)
    n = len(ytr)
    p = Vt.shape[1]
    if n < 5 or S.size == 0:
        mu = float(ytr.mean()) if n else np.nan
        return np.full(len(Xte_s), mu), np.zeros(p), mu, np.nan
    ymu = float(ytr.mean())
    yc = ytr - ymu
    Uy = U.T @ yc
    best_a, best_sse = float(alphas[0]), np.inf
    U2 = U ** 2
    for a in alphas:
        aa = float(a)
        lam = (S ** 2) / (S ** 2 + aa)
        h = U2 @ lam
        fit = U @ (lam * Uy)
        denom = 1.0 - h
        denom[np.abs(denom) < 1e-12] = 1e-12
        sse = float(np.sum(((yc - fit) / denom) ** 2))
        if sse < best_sse:
            best_sse = sse
            best_a = aa
    d = S / (S ** 2 + best_a)
    coef = Vt.T @ (d * Uy)
    pred = ymu + np.asarray(Xte_s, np.float64) @ coef
    return pred, coef, ymu, best_a


def svd_train(Xtr):
    Xtr = np.asarray(Xtr, np.float64)
    U, S, Vt = np.linalg.svd(Xtr, full_matrices=False)
    keep = S > 1e-10
    if not np.any(keep):
        return U[:, :0], S[:0], Vt[:0]
    return U[:, keep], S[keep], Vt[keep]


def type_masks(ct, idx, t):
    return idx[ct[idx] == t]


def nearest_centroid_predict(Ztr, ytr, Zte, yte=None):
    """Nearest type-centroid in the identity subspace. ytr are train cell-type labels."""
    ytr = np.asarray(ytr).astype(str)
    types = sorted(pd.unique(ytr))
    cents, used = [], []
    for t in types:
        m = ytr == t
        if int(m.sum()) < 1:
            continue
        cents.append(Ztr[m].mean(0))
        used.append(t)
    if not used:
        n = len(Zte)
        return np.array(["__none__"] * n, dtype=object), np.full(n, np.nan), np.zeros((0, Ztr.shape[1])), []
    C = np.vstack(cents)
    d2 = ((Zte[:, None, :] - C[None, :, :]) ** 2).sum(-1)
    j = np.argmin(d2, axis=1)
    pred = np.array([used[i] for i in j], dtype=object)
    dist_pred = np.sqrt(d2[np.arange(len(Zte)), j])
    dist_own = np.full(len(Zte), np.nan)
    if yte is not None:
        yte = np.asarray(yte).astype(str)
        col = {t: i for i, t in enumerate(used)}
        for i, lab in enumerate(yte):
            k = col.get(lab)
            if k is not None:
                dist_own[i] = float(np.sqrt(d2[i, k]))
    return pred, dist_own, C, used


def identity_project(X, basis, mu_row=None):
    """Coordinates in the identity subspace. X is already standardized in the same space as basis."""
    X = np.asarray(X, np.float64)
    if mu_row is not None:
        X = X - mu_row
    if basis.size == 0:
        return np.zeros((len(X), 0))
    return X @ basis


def residualize_identity(X, basis):
    """X - P_id X, without forming the p x p projector."""
    X = np.asarray(X, np.float64)
    if basis.size == 0:
        return X
    Z = X @ basis
    return X - Z @ basis.T


def permutation_p(obs_v, null_v, greater=True):
    null_v = np.asarray(null_v, float)
    null_v = null_v[np.isfinite(null_v)]
    if (not np.isfinite(obs_v)) or len(null_v) == 0:
        return np.nan
    if greater:
        return float((np.sum(null_v >= obs_v) + 1) / (len(null_v) + 1))
    return float((np.sum(null_v <= obs_v) + 1) / (len(null_v) + 1))


def summarize_null_col(v):
    v = np.asarray(v, float)
    v = v[np.isfinite(v)]
    if len(v) == 0:
        return dict(mean=np.nan, median=np.nan, p05=np.nan, p95=np.nan, min=np.nan, max=np.nan, n=0)
    return dict(
        mean=float(np.mean(v)), median=float(np.median(v)),
        p05=float(np.percentile(v, 5)), p95=float(np.percentile(v, 95)),
        min=float(np.min(v)), max=float(np.max(v)), n=int(len(v)),
    )


def scores_json(sc):
    skip = {"per_type", "donor_pred", "pred", "y"}
    out = {}
    for k, v in sc.items():
        if k in skip:
            continue
        if isinstance(v, (np.floating, float)):
            out[k] = float(v) if np.isfinite(v) else None
        else:
            out[k] = v
    return out


def primary_mae(scheme, scores):
    pt = scores["per_type"]
    if scheme == "within_site":
        sub = pt[pt.site == "mean_of_sites"]
        v = sub.mae.to_numpy(float) if len(sub) else np.array([])
        v = v[np.isfinite(v)]
        return float(np.median(v)) if len(v) else np.nan
    return scores.get("pooled_median_type_mae", np.nan)


def identity_metrics(y_true, y_pred, chance):
    y_true = np.asarray(y_true).astype(str)
    y_pred = np.asarray(y_pred).astype(str)
    m = np.array([p != "None" and p != "__none__" for p in y_pred])
    if m.sum() < 1:
        return dict(acc=np.nan, macro_f1=np.nan, n=0, chance=chance)
    yt, yp = y_true[m], y_pred[m]
    acc = float((yt == yp).mean())
    types = sorted(set(yt) | set(yp))
    f1s = []
    for t in types:
        tp = int(((yt == t) & (yp == t)).sum())
        fp = int(((yt != t) & (yp == t)).sum())
        fn = int(((yt == t) & (yp != t)).sum())
        prec = tp / (tp + fp) if (tp + fp) else 0.0
        rec = tp / (tp + fn) if (tp + fn) else 0.0
        f1s.append(0.0 if (prec + rec) == 0 else 2 * prec * rec / (prec + rec))
    return dict(acc=acc, macro_f1=float(np.mean(f1s)) if f1s else np.nan, n=int(m.sum()), chance=float(chance))


def chance_identity(obs):
    vc = obs.celltype.astype(str).value_counts(normalize=True)
    n = max(obs.celltype.nunique(), 1)
    return float(max(1.0 / n, vc.max() if len(vc) else 0.0))


def fmt(x, d=3):
    try:
        x = float(x)
    except (TypeError, ValueError):
        return "NA"
    if not np.isfinite(x):
        return "NA"
    return f"{x:+.{d}f}"


def fmt_u(x, d=3):
    s = fmt(x, d)
    return s.lstrip("+") if s != "NA" else s
