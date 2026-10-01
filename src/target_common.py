"""Target vector — expression space only.

Direction that opposes aging while displacing identity as little as possible.
No gene selection: every gene in the Phase-1 log-CPM matrix carries a weight.
Namespaced under results/target/. Does not modify prior FINDINGS*.md or FALSIFICATION.md.
Seed 20260914 (same fold construction as Phase 1 / geometry / trajectory).

Working space: training-fold (or analysis-set) GLOBAL z-score of all genes.
Identity subspace and age directions are both fit in that space so TARGET
= (I - P_id) w_age is a well-defined projection. Ridge does NOT re-z-score
per cell type (that would put coefficients in a different metric than P_id).
"""
from __future__ import annotations

import sys
from pathlib import Path

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parent))
from config import RESULTS  # noqa: E402
from geometry_common import (  # noqa: E402
    GEO_SEED, DATASET_ID, NULL_USABLE, ALPHAS, HASH_POOL_LIMITATION,
    Logger, dump_json, load_json, load_phase1_matrix, make_folds,
    permute_age_within_site_obs, score_age_predictions, primary_r2,
)
from geometry_partC import (  # noqa: E402
    identity_basis, pls1_direction, proj_r2, angle_to_subspace_deg,
    pairwise_angles_deg, consensus_direction,
)
from trajectory_common import (  # noqa: E402
    zscore_train, calibrate_1d, permutation_p, summarize_null_col,
    fmt, fmt_u, scores_json, svd_train, primary_mae,
)
from brain_phase1_common import r2_mae  # noqa: E402

TARGET_DIR = RESULTS / "target"
TARGET_FIG = TARGET_DIR / "figures"
for _p in (TARGET_DIR, TARGET_FIG):
    _p.mkdir(parents=True, exist_ok=True)

TARGET_SEED = GEO_SEED  # 20260914
N_BOOT_RIDGE = 20
N_BOOT_PLS = 40
N_PERM = 20
N_PERM_PAIR = 50
STOP_BOOT_ANGLE = 60.0  # degrees; median pairwise bootstrap angle
STOP_NULL = NULL_USABLE  # 0.05
MIN_N_TYPE = 8

# V3: consensus is "near" per-type if it keeps most of the held-out R²
V3_RATIO_NEAR = 0.80
V3_LOSS_NEAR = 0.10
V3_RATIO_LOST = 0.50


def tgt_log_banner(log, stage: str):
    log("=" * 100)
    log(f"TARGET {stage}  seed={TARGET_SEED}  dataset_id={DATASET_ID}")
    log("=" * 100)
    log(HASH_POOL_LIMITATION)
    log("NO GENE SELECTION: every gene in the Phase-1 matrix enters every model with a weight.")
    log("NO perturbation data. No TF Atlas. Expression space only.")
    log("FALSIFICATION.md and prior FINDINGS*.md are not modified.")
    log("Working space: global train z-score; ridge coefficients live in that space (no second per-type z-score).")


def unsigned_angle_deg(a, b):
    a = np.asarray(a, np.float64).ravel()
    b = np.asarray(b, np.float64).ravel()
    na, nb = float(np.linalg.norm(a)), float(np.linalg.norm(b))
    if na < 1e-12 or nb < 1e-12 or not np.isfinite(a).all() or not np.isfinite(b).all():
        return np.nan
    d = float(np.clip(np.abs(np.dot(a, b) / (na * nb)), 0.0, 1.0))
    return float(np.degrees(np.arccos(d)))


def unit(w):
    w = np.asarray(w, np.float64).ravel()
    n = float(np.linalg.norm(w))
    if n < 1e-12 or not np.isfinite(n):
        return w, 0.0
    return w / n, n


def align_age_sign(w, X, y):
    """Flip w so cov(X @ w, y) >= 0 (older is the positive direction)."""
    w = np.asarray(w, np.float64).ravel()
    if not np.isfinite(w).all() or float(np.linalg.norm(w)) < 1e-12:
        return w
    s = np.asarray(X, np.float64) @ w
    y = np.asarray(y, np.float64)
    if len(y) < 3 or float(np.std(s)) < 1e-12:
        return w
    if float(np.cov(s, y, ddof=0)[0, 1]) < 0:
        return -w
    return w


def target_from_age(w, basis):
    """Remove the identity-subspace component of unit age direction w; renormalize.

    Returns (target_unit, norm_retained, angle_to_age_deg, proj_r2_on_id).
    TARGET is aligned with aging (positive = older). A rejuvenating move is -TARGET.
    """
    w, nrm = unit(w)
    if nrm < 1e-12:
        return w, 0.0, np.nan, np.nan
    pr = proj_r2(w, basis)
    if basis.size == 0:
        return w, 1.0, 0.0, 0.0
    z = basis.T @ w
    resid = w - basis @ z
    retained = float(np.linalg.norm(resid))
    if retained < 1e-12:
        return np.zeros_like(w), 0.0, 90.0, float(pr) if np.isfinite(pr) else 1.0
    tgt = resid / retained
    ang = unsigned_angle_deg(tgt, w)
    return tgt, retained, ang, float(pr) if np.isfinite(pr) else np.nan


def ridge_prestd(X, y, alphas=ALPHAS, svd=None):
    """Ridge on already-standardized X via the n×n Gram (no p×n Vᵀ).

    Cache is (U, S2) from eigh(X Xᵀ), tiny compared with SVD's Vᵀ.
    Returns coef, intercept, alpha, cache.
    """
    X = np.asarray(X, np.float64)
    y = np.asarray(y, np.float64)
    n, p = X.shape
    if n < 5 or p == 0:
        return np.zeros(p), float(y.mean()) if n else np.nan, np.nan, None
    ymu = float(y.mean())
    yc = y - ymu
    if svd is None:
        G = X @ X.T
        G = 0.5 * (G + G.T)
        evals, U = np.linalg.eigh(G)
        keep = evals > 1e-10
        if not np.any(keep):
            return np.zeros(p), ymu, np.nan, None
        U, S2 = U[:, keep], np.clip(evals[keep], 0.0, None)
    else:
        U, S2 = svd
        if S2.size == 0:
            return np.zeros(p), ymu, np.nan, (U, S2)
    Uy = U.T @ yc
    best_a, best_sse = float(alphas[0]), np.inf
    U2 = U ** 2
    for a in alphas:
        aa = float(a)
        lam = S2 / (S2 + aa)
        h = U2 @ lam
        fit = U @ (lam * Uy)
        denom = 1.0 - h
        denom[np.abs(denom) < 1e-12] = 1e-12
        sse = float(np.sum(((yc - fit) / denom) ** 2))
        if sse < best_sse:
            best_sse = sse
            best_a = aa
    # coef = X' U (Uy / (S2 + α))  — equivalent to V' (S/(S²+α) U' yc)
    coef = X.T @ (U @ (Uy / (S2 + best_a)))
    return coef, ymu, best_a, (U, S2)


def age_direction(X, y, method="ridge", svd=None):
    """Unit age direction in the same space as X. method in {ridge, pls1}."""
    X = np.asarray(X, np.float64)
    y = np.asarray(y, np.float64)
    p = X.shape[1]
    if len(y) < MIN_N_TYPE or p == 0:
        return np.full(p, np.nan), dict(ok=False, n=int(len(y)), method=method)
    if method == "pls1":
        w, nrm = pls1_direction(X, y)
        w = align_age_sign(w, X, y)
        w, _ = unit(w)
        return w, dict(ok=True, n=int(len(y)), method="pls1", norm=float(nrm), alpha=np.nan)
    coef, ymu, a, svd_out = ridge_prestd(X, y, svd=svd)
    w, nrm = unit(coef)
    w = align_age_sign(w, X, y)
    w, _ = unit(w)
    return w, dict(ok=True, n=int(len(y)), method="ridge", norm=float(nrm),
                   alpha=float(a) if np.isfinite(a) else np.nan, intercept=float(ymu),
                   svd=svd_out, coef=coef)


def per_type_directions(X, y, ct, types, method="ridge", svd_cache=None):
    """W: p x T. svd_cache: optional dict type -> (U,S,Vt) for ridge perms."""
    p = X.shape[1]
    W = np.full((p, len(types)), np.nan)
    meta = []
    built = {} if svd_cache is None else dict(svd_cache)
    for j, t in enumerate(types):
        m = ct == t
        n = int(m.sum())
        if n < MIN_N_TYPE:
            meta.append(dict(celltype=t, n=n, ok=False, method=method))
            continue
        svd = None if svd_cache is None else svd_cache.get(t)
        w, rec = age_direction(X[m], y[m], method=method, svd=svd)
        W[:, j] = w
        if method == "ridge" and svd_cache is None and rec.get("ok") and rec.get("svd") is not None:
            built[t] = rec["svd"]
        rec = {k: v for k, v in rec.items() if k not in ("svd", "coef")}
        rec["celltype"] = t
        meta.append(rec)
    return W, pd.DataFrame(meta), built


def per_type_targets(W_age, basis):
    p, T = W_age.shape
    W_tgt = np.full((p, T), np.nan)
    rows = []
    for j in range(T):
        w = W_age[:, j]
        if not np.isfinite(w).all() or float(np.linalg.norm(w)) < 1e-12:
            rows.append(dict(ok=False, norm_retained=np.nan, angle_to_age_deg=np.nan,
                             proj_r2_id=np.nan, angle_to_id_deg=np.nan))
            continue
        tgt, retained, ang, pr = target_from_age(w, basis)
        W_tgt[:, j] = tgt
        rows.append(dict(ok=True, norm_retained=retained, angle_to_age_deg=ang,
                         proj_r2_id=pr, angle_to_id_deg=angle_to_subspace_deg(w, basis)))
    return W_tgt, pd.DataFrame(rows)


def donor_bootstrap_index(obs, rng):
    """Resample donors with replacement; concatenate all of each draw's rows (multiplicity)."""
    donors = obs.donor.astype(str).to_numpy()
    uniq = pd.unique(donors)
    draw = rng.choice(uniq, size=len(uniq), replace=True)
    # map donor -> row indices
    groups = {}
    for i, d in enumerate(donors):
        groups.setdefault(d, []).append(i)
    parts = [np.asarray(groups[d], dtype=int) for d in draw]
    return np.concatenate(parts)


def pairwise_median_deg(W):
    """Median unsigned pairwise angle among finite unit columns of W (p x T)."""
    ok = []
    for j in range(W.shape[1]):
        w = W[:, j]
        if np.isfinite(w).all() and float(np.linalg.norm(w)) > 1e-12:
            ok.append(j)
    if len(ok) < 2:
        return np.nan, None
    Wc = np.column_stack([W[:, j] / np.linalg.norm(W[:, j]) for j in ok])
    ang = pairwise_angles_deg(Wc)
    iu = np.triu_indices(ang.shape[0], 1)
    return float(np.median(ang[iu])), ang


def oof_project_both(Y, obs, folds, method, y_override=None, svd_cache=None):
    """OOF predictions onto raw age direction AND TARGET from one train fit."""
    pred_age = np.full(len(obs), np.nan)
    pred_tgt = np.full(len(obs), np.nan)
    y = obs.age.to_numpy(float) if y_override is None else np.asarray(y_override, float)
    ct = obs.celltype.astype(str).to_numpy()
    types = sorted(pd.unique(ct))
    new_cache = {} if svd_cache is None and method == "ridge" else None
    for fi, fold in enumerate(folds):
        tr, te = fold[0], fold[1]
        Xtr, Xte, _, _ = zscore_train(Y[tr], Y[te])
        basis, _, _ = identity_basis(Xtr, ct[tr], types)
        ytr = y[tr]
        cache_f = None if svd_cache is None else svd_cache.get(fi)
        W_age, _, built = per_type_directions(
            Xtr, ytr, ct[tr], types, method=method, svd_cache=cache_f)
        W_tgt, _ = per_type_targets(W_age, basis)
        if new_cache is not None:
            new_cache[fi] = built
        for j, t in enumerate(types):
            tr_m = ct[tr] == t
            te_m = ct[te] == t
            if int(tr_m.sum()) < MIN_N_TYPE or int(te_m.sum()) < 1:
                continue
            for W, pred in ((W_age, pred_age), (W_tgt, pred_tgt)):
                w = W[:, j]
                if not np.isfinite(w).all() or float(np.linalg.norm(w)) < 1e-12:
                    continue
                pte, _, _ = calibrate_1d(Xtr[tr_m] @ w, ytr[tr_m], Xte[te_m] @ w)
                pred[te[te_m]] = pte
    out = dict(pred_age=pred_age, pred_tgt=pred_tgt, y=y, types=types)
    if new_cache is not None:
        out["svd_cache"] = new_cache
    return out


def oof_project(Y, obs, folds, method, use_target=True, y_override=None,
                svd_cache=None, collect_W=False):
    """OOF 1D OLS of age ~ projection onto per-type (TARGET or raw age) direction.

    Directions and identity basis fit on the training fold only, in train z-score space.
    svd_cache: optional {fold_i: {type: (U,S,Vt)}} for ridge permutations.
    """
    pred = np.full(len(obs), np.nan)
    y = obs.age.to_numpy(float) if y_override is None else np.asarray(y_override, float)
    ct = obs.celltype.astype(str).to_numpy()
    types = sorted(pd.unique(ct))
    new_cache = {} if svd_cache is None and method == "ridge" else None
    fold_pack = [] if collect_W else None
    for fi, fold in enumerate(folds):
        tr, te = fold[0], fold[1]
        Xtr, Xte, mu, sd = zscore_train(Y[tr], Y[te])
        basis, _, _ = identity_basis(Xtr, ct[tr], types)
        ytr = y[tr]
        cache_f = None if svd_cache is None else svd_cache.get(fi)
        W_age, meta, built = per_type_directions(
            Xtr, ytr, ct[tr], types, method=method, svd_cache=cache_f)
        if use_target:
            W, geom = per_type_targets(W_age, basis)
        else:
            W, geom = W_age, None
        if new_cache is not None:
            new_cache[fi] = built
        if collect_W:
            fold_pack.append(dict(
                fold=fi, tag=str(fold[2] if len(fold) > 2 else fi),
                W_age=W_age, W_tgt=(W if use_target else W_age),
                basis=basis, mu=mu, sd=sd, types=types,
            ))
        for j, t in enumerate(types):
            tr_m = ct[tr] == t
            te_m = ct[te] == t
            w = W[:, j]
            if int(tr_m.sum()) < MIN_N_TYPE or int(te_m.sum()) < 1:
                continue
            if not np.isfinite(w).all() or float(np.linalg.norm(w)) < 1e-12:
                continue
            s_tr = Xtr[tr_m] @ w
            s_te = Xte[te_m] @ w
            pte, _, _ = calibrate_1d(s_tr, ytr[tr_m], s_te)
            pred[te[te_m]] = pte
    out = dict(pred=pred, y=y, types=types)
    if new_cache is not None:
        out["svd_cache"] = new_cache
    if collect_W:
        out["folds"] = fold_pack
    return out


def oof_project_fixed_W(Y, obs, folds, W, y_override=None):
    """Project onto a FROZEN per-type matrix W (p x T), still z-scoring on train only.

    Used for the consensus vector (same w for every type: pass a p x T with repeated columns)
    and for site-transfer / age-split evaluations that apply a pre-fit direction.
    """
    pred = np.full(len(obs), np.nan)
    y = obs.age.to_numpy(float) if y_override is None else np.asarray(y_override, float)
    ct = obs.celltype.astype(str).to_numpy()
    types = sorted(pd.unique(ct))
    assert W.shape[1] == len(types)
    for fold in folds:
        tr, te = fold[0], fold[1]
        Xtr, Xte, _, _ = zscore_train(Y[tr], Y[te])
        ytr = y[tr]
        for j, t in enumerate(types):
            tr_m = ct[tr] == t
            te_m = ct[te] == t
            w = W[:, j]
            if int(tr_m.sum()) < MIN_N_TYPE or int(te_m.sum()) < 1:
                continue
            if not np.isfinite(w).all() or float(np.linalg.norm(w)) < 1e-12:
                continue
            s_tr = Xtr[tr_m] @ w
            s_te = Xte[te_m] @ w
            pte, _, _ = calibrate_1d(s_tr, ytr[tr_m], s_te)
            pred[te[te_m]] = pte
    return pred


def oof_consensus_vs_own(Y, obs, folds, method="ridge", y_override=None, svd_cache=None):
    """Per fold: fit per-type TARGETs on train, consensus = PC1, predict test with both."""
    pred_own = np.full(len(obs), np.nan)
    pred_cons = np.full(len(obs), np.nan)
    y = obs.age.to_numpy(float) if y_override is None else np.asarray(y_override, float)
    ct = obs.celltype.astype(str).to_numpy()
    types = sorted(pd.unique(ct))
    new_cache = {} if svd_cache is None and method == "ridge" else None
    pc1_fracs = []
    for fi, fold in enumerate(folds):
        tr, te = fold[0], fold[1]
        Xtr, Xte, _, _ = zscore_train(Y[tr], Y[te])
        basis, _, _ = identity_basis(Xtr, ct[tr], types)
        ytr = y[tr]
        cache_f = None if svd_cache is None else svd_cache.get(fi)
        W_age, _, built = per_type_directions(
            Xtr, ytr, ct[tr], types, method=method, svd_cache=cache_f)
        if new_cache is not None:
            new_cache[fi] = built
        W_tgt, _ = per_type_targets(W_age, basis)
        ok = [j for j in range(len(types))
              if np.isfinite(W_tgt[:, j]).all() and float(np.linalg.norm(W_tgt[:, j])) > 1e-12]
        if len(ok) < 2:
            w_cons, pc1 = None, np.nan
        else:
            Wc = np.column_stack([W_tgt[:, j] / np.linalg.norm(W_tgt[:, j]) for j in ok])
            w_cons, pc1 = consensus_direction(Wc)
            # sign: positive mean overlap with per-type targets
            dots = Wc.T @ w_cons
            if float(np.nanmean(dots)) < 0:
                w_cons = -w_cons
        pc1_fracs.append(pc1)
        for j, t in enumerate(types):
            tr_m = ct[tr] == t
            te_m = ct[te] == t
            if int(tr_m.sum()) < MIN_N_TYPE or int(te_m.sum()) < 1:
                continue
            w = W_tgt[:, j]
            if np.isfinite(w).all() and float(np.linalg.norm(w)) > 1e-12:
                pte, _, _ = calibrate_1d(Xtr[tr_m] @ w, ytr[tr_m], Xte[te_m] @ w)
                pred_own[te[te_m]] = pte
            if w_cons is not None:
                pte, _, _ = calibrate_1d(Xtr[tr_m] @ w_cons, ytr[tr_m], Xte[te_m] @ w_cons)
                pred_cons[te[te_m]] = pte
    out = dict(pred_own=pred_own, pred_cons=pred_cons, y=y,
               pc1_frac=float(np.nanmean(pc1_fracs)) if pc1_fracs else np.nan)
    if new_cache is not None:
        out["svd_cache"] = new_cache
    return out


def per_type_primary_r2(obs, pred, scheme, y=None):
    """One R² per cell type under the scheme's primary metric."""
    sc = score_age_predictions(obs, pred, y=y)
    pt = sc["per_type"]
    if scheme == "within_site":
        sub = pt[pt.site == "mean_of_sites"]
    else:
        sub = pt[pt.site == "pooled"]
    return sub.set_index("celltype")["r2"].to_dict(), sc


def shuffle_null_oof(Y, obs, folds, method, use_target, rng, n_perm, svd_cache, scheme):
    rows, per_type_rows = [], []
    for i in range(n_perm):
        obs_p = permute_age_within_site_obs(obs, rng)
        yp = obs_p.age.to_numpy(float)
        out = oof_project(Y, obs, folds, method, use_target=use_target, y_override=yp,
                          svd_cache=svd_cache)
        rmap, sc = per_type_primary_r2(obs_p, out["pred"], scheme, y=yp)
        rows.append(dict(
            perm=i, pooled=sc["pooled_median_type_r2"],
            site_strat=sc["site_strat_median_type_r2"],
            primary=primary_r2(scheme, sc),
        ))
        for t, r in rmap.items():
            per_type_rows.append(dict(perm=i, celltype=t, r2=r))
        if (i + 1) % 5 == 0:
            print(f"      perm {i+1}/{n_perm}", flush=True)
    return pd.DataFrame(rows), pd.DataFrame(per_type_rows)


def type_loses_age(obs_r2, null_r2s):
    """True if held-out R² is at or below the permutation null mean."""
    null_r2s = np.asarray(null_r2s, float)
    null_r2s = null_r2s[np.isfinite(null_r2s)]
    if not np.isfinite(obs_r2) or len(null_r2s) == 0:
        return True
    return bool(obs_r2 <= float(np.mean(null_r2s)))


def fit_full(Y, obs, method="ridge"):
    """Full-cohort descriptive fit (not for evaluation)."""
    ct = obs.celltype.astype(str).to_numpy()
    types = sorted(pd.unique(ct))
    y = obs.age.to_numpy(float)
    Xs, mu, sd = zscore_train(Y)
    basis, S_id, used = identity_basis(Xs, ct, types)
    W_age, meta, built = per_type_directions(Xs, y, ct, types, method=method)
    W_tgt, geom = per_type_targets(W_age, basis)
    geom = geom.drop(columns=[c for c in ("ok",) if c in geom.columns])
    geom = pd.concat([meta.reset_index(drop=True), geom.reset_index(drop=True)], axis=1)
    return dict(
        X=Xs, mu=mu, sd=sd, y=y, ct=ct, types=types,
        basis=basis, S_id=S_id, used=used,
        W_age=W_age, W_tgt=W_tgt, meta=geom, svd=built,
    )


def joint_within_type_direction(X, y, ct, types, method="ridge"):
    """One age direction on type-mean-centered X (and type-centered y). Already ~orthogonal to identity."""
    X = np.asarray(X, np.float64).copy()
    y = np.asarray(y, np.float64).copy()
    for t in types:
        m = ct == t
        if int(m.sum()) >= 3:
            X[m] = X[m] - X[m].mean(0)
            y[m] = y[m] - y[m].mean()
    w, rec = age_direction(X, y, method=method)
    return w, rec


def save_matrix_bundle(path, genes, types, **arrays):
    payload = dict(
        gene_id=genes.gene_id.astype(str).to_numpy(),
        symbol=genes.symbol.astype(str).to_numpy(),
        types=np.array(types, dtype=object),
        **arrays,
    )
    np.savez_compressed(path, **payload)


def md_table(df, cols=None):
    if cols is not None:
        df = df[cols]
    lines = ["| " + " | ".join(df.columns) + " |", "|" + "|".join("---" for _ in df.columns) + "|"]
    for _, row in df.iterrows():
        cells = []
        for c in df.columns:
            v = row[c]
            if isinstance(v, (float, np.floating)):
                cells.append("NA" if not np.isfinite(v) else f"{float(v):+.3f}")
            elif isinstance(v, (bool, np.bool_)):
                cells.append(str(bool(v)))
            elif v is None:
                cells.append("NA")
            else:
                cells.append(str(v))
        lines.append("| " + " | ".join(cells) + " |")
    return "\n".join(lines)
