"""PART C — geometric reframe (age and identity as directions, not gene buckets).

C1. Per-cell-type age direction (PLS-1 / ridge coefficients).
C2. Identity subspace (type-centroid SVD).
C3. Angle of the age direction with the identity subspace, plus a donor-age permutation null.
C4. Pairwise angles among per-type age directions (shared vs type-specific).

Both CV schemes: directions estimated on each training fold; full-data descriptive
reported alongside. Does not derive new gene sets — all genes in the Phase-1 matrix.

Usage: python src/geometry_partC.py
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
from geometry_common import (  # noqa: E402
    GEO_DIR, GEO_FIG, GEO_SEED, N_PERM_C, NULL_USABLE, DATASET_ID,
    Logger, dump_json, load_json, geo_log_banner, load_phase1_matrix, make_folds,
    permute_age_within_site_obs, ridge_svd_fit,
)


def _zscore(X, mu=None, sd=None):
    X = np.asarray(X, np.float64)
    if mu is None:
        mu = X.mean(0)
        sd = X.std(0)
        sd[sd < 1e-12] = 1.0
    return (X - mu) / sd, mu, sd


def identity_basis(X, ct, types):
    """Orthonormal basis of the between-type centroid subspace. X is standardized."""
    cents = []
    used = []
    for t in types:
        m = ct == t
        if int(m.sum()) < 3:
            continue
        cents.append(X[m].mean(0))
        used.append(t)
    C = np.vstack(cents)
    C = C - C.mean(0, keepdims=True)
    # C is n_types x p; right singular vectors span the centroid subspace
    _, S, Vt = np.linalg.svd(C, full_matrices=False)
    rank = int(np.sum(S > 1e-8))
    rank = max(rank, 0)
    if rank == 0:
        return np.zeros((X.shape[1], 0)), S, used
    return Vt[:rank].T, S[:rank], used


def pls1_direction(X, y):
    """Unit covariance direction (univariate PLS / age-weighted mean gene)."""
    X = np.asarray(X, np.float64)
    y = np.asarray(y, np.float64)
    if len(y) < 5 or X.shape[1] == 0:
        return np.zeros(X.shape[1]), 0.0
    Xc = X - X.mean(0)
    yc = y - y.mean()
    w = Xc.T @ yc
    nrm = float(np.linalg.norm(w))
    if nrm < 1e-12:
        return w, 0.0
    return w / nrm, nrm


def ridge_direction(X, y):
    coef, _, _, _, a = ridge_svd_fit(X, y)
    nrm = float(np.linalg.norm(coef))
    if nrm < 1e-12:
        return coef, a, 0.0
    return coef / nrm, a, nrm


def proj_r2(w, basis):
    if basis.size == 0 or not np.isfinite(w).any():
        return np.nan
    nrm = float(np.linalg.norm(w))
    if nrm < 1e-12:
        return np.nan
    w = w / nrm
    return float(np.sum((basis.T @ w) ** 2))


def angle_to_subspace_deg(w, basis):
    r = proj_r2(w, basis)
    if not np.isfinite(r):
        return np.nan
    r = float(np.clip(r, 0.0, 1.0))
    return float(np.degrees(np.arccos(np.sqrt(r))))


def pairwise_angles_deg(W):
    """W: p x T unit columns. Return T x T matrix of unsigned angles in [0, 90]."""
    T = W.shape[1]
    dots = np.clip(np.abs(W.T @ W), 0.0, 1.0)
    ang = np.degrees(np.arccos(dots))
    np.fill_diagonal(ang, 0.0)
    return ang


def consensus_direction(W):
    """PC1 of the per-type age directions (columns of W, unit)."""
    if W.shape[1] < 2:
        w = W[:, 0]
        n = float(np.linalg.norm(w))
        return (w / n) if n > 0 else w, 1.0
    # SVD of T x p
    _, S, Vt = np.linalg.svd(W.T, full_matrices=False)
    ev = (S ** 2)
    ev = ev / ev.sum() if ev.sum() > 0 else ev
    w = Vt[0]
    n = float(np.linalg.norm(w))
    return (w / n) if n > 0 else w, float(ev[0]) if len(ev) else np.nan


def directions_on_mask(X, y, ct, types, mask, ridge=False):
    """Age directions for each type on a boolean/index mask of rows.

    PLS-1 is the primary direction (cheap, well-defined for p >> n). Ridge is
    optional sensitivity — do not enable it inside permutation/CV loops.
    """
    W_pls, W_ridge = [], []
    rows = []
    idx = np.flatnonzero(mask) if mask.dtype == bool else np.asarray(mask)
    X_m, y_m, ct_m = X[idx], y[idx], ct[idx]
    p = X.shape[1]
    for t in types:
        m = ct_m == t
        if int(m.sum()) < 8:
            W_pls.append(np.full(p, np.nan))
            if ridge:
                W_ridge.append(np.full(p, np.nan))
            rows.append(dict(celltype=t, n=int(m.sum()), ok=False))
            continue
        w_p, nrm_p = pls1_direction(X_m[m], y_m[m])
        W_pls.append(w_p)
        rec = dict(celltype=t, n=int(m.sum()), ok=True, pls_norm=nrm_p)
        if ridge:
            w_r, a, nrm_r = ridge_direction(X_m[m], y_m[m])
            W_ridge.append(w_r)
            rec.update(ridge_norm=nrm_r, ridge_alpha=a)
        rows.append(rec)
    Wp = np.column_stack(W_pls)
    Wr = np.column_stack(W_ridge) if ridge else None
    return Wp, Wr, pd.DataFrame(rows)


def score_geometry(W, basis, types):
    """Per-type projection R² / angle, plus consensus."""
    rows = []
    ok_cols = []
    for j, t in enumerate(types):
        w = W[:, j]
        if not np.isfinite(w).all() or np.linalg.norm(w) < 1e-12:
            rows.append(dict(celltype=t, proj_r2=np.nan, angle_deg=np.nan, ok=False))
            continue
        pr = proj_r2(w, basis)
        rows.append(dict(celltype=t, proj_r2=pr, angle_deg=angle_to_subspace_deg(w, basis), ok=True))
        ok_cols.append(j)
    tab = pd.DataFrame(rows)
    if ok_cols:
        Wc = W[:, ok_cols]
        # re-unit (ridge/pls already unit, but NaN cols dropped)
        for j in range(Wc.shape[1]):
            n = np.linalg.norm(Wc[:, j])
            if n > 0:
                Wc[:, j] = Wc[:, j] / n
        w_cons, pc1 = consensus_direction(Wc)
        ang = pairwise_angles_deg(Wc)
        iu = np.triu_indices(ang.shape[0], 1)
        cons_r2 = proj_r2(w_cons, basis)
        cons_ang = angle_to_subspace_deg(w_cons, basis)
        pair_med = float(np.median(ang[iu])) if len(iu[0]) else np.nan
        pair_mean = float(np.mean(ang[iu])) if len(iu[0]) else np.nan
    else:
        w_cons, pc1 = None, np.nan
        ang = None
        cons_r2 = cons_ang = pair_med = pair_mean = np.nan
    return dict(
        per_type=tab,
        median_proj_r2=float(tab.proj_r2.median()) if tab.ok.any() else np.nan,
        median_angle_deg=float(tab.angle_deg.median()) if tab.ok.any() else np.nan,
        mean_proj_r2=float(tab.proj_r2.mean()) if tab.ok.any() else np.nan,
        consensus_proj_r2=cons_r2,
        consensus_angle_deg=cons_ang,
        pairwise_median_deg=pair_med,
        pairwise_mean_deg=pair_mean,
        shared_pc1_frac=pc1,
        pairwise=ang,
        n_ok=int(tab.ok.sum()),
    )


def run():
    rng = np.random.default_rng(GEO_SEED)
    log = Logger(GEO_DIR / "partC_report.txt")
    try:
        return _run(log, rng)
    finally:
        log.close()


def _run(log, rng):
    geo_log_banner(log, "PART C — geometric reframe")
    a_sum = GEO_DIR / "partA_summary.json"
    if a_sum.exists() and load_json(a_sum).get("stop_A"):
        log("[C] STOP A is set — refusing to run Part C on an unsound pipeline.")
        dump_json(GEO_DIR / "partC_summary.json", dict(skipped=True, reason="STOP A"))
        return dict(skipped=True)

    data = load_phase1_matrix(log)
    Y, obs = data["Y"], data["obs"]
    log(f"[filter] geometry uses ALL {Y.shape[1]:,} genes; no A/I bucketing.")
    ct = obs.celltype.astype(str).to_numpy()
    types = sorted(pd.unique(ct))
    y = obs.age.to_numpy(float)
    donors = obs.donor.astype(str).to_numpy()
    X, mu, sd = _zscore(Y)
    log(f"[C] standardized {X.shape[0]} x {X.shape[1]} (global mean/sd on full matrix for descriptive)")

    basis, S_id, used = identity_basis(X, ct, types)
    log(f"[C2] identity subspace: {basis.shape[1]} axes from {len(used)} type centroids; "
        f"singular values head={np.round(S_id[:5], 3).tolist()}")

    # ------------------------------------------------------------------ C1/C3 full-data descriptive
    W_pls, W_ridge, meta = directions_on_mask(
        X, y, ct, types, np.ones(len(obs), bool), ridge=True)
    geo_pls = score_geometry(W_pls, basis, types)
    geo_ridge = score_geometry(W_ridge, basis, types)
    geo_pls["per_type"].assign(method="pls1").merge(meta, on="celltype").to_csv(
        GEO_DIR / "partC_per_type_full_pls.csv", index=False)
    geo_ridge["per_type"].assign(method="ridge").merge(meta, on="celltype").to_csv(
        GEO_DIR / "partC_per_type_full_ridge.csv", index=False)
    log(f"[C3] FULL PLS-1  median angle={geo_pls['median_angle_deg']:.1f}°  "
        f"median proj R²={geo_pls['median_proj_r2']:.3f}  "
        f"consensus angle={geo_pls['consensus_angle_deg']:.1f}°  "
        f"PC1 shared={geo_pls['shared_pc1_frac']:.3f}")
    log(f"[C3] FULL ridge  median angle={geo_ridge['median_angle_deg']:.1f}°  "
        f"median proj R²={geo_ridge['median_proj_r2']:.3f}  "
        f"consensus angle={geo_ridge['consensus_angle_deg']:.1f}°  "
        f"PC1 shared={geo_ridge['shared_pc1_frac']:.3f}")

    if geo_pls["pairwise"] is not None:
        pair = pd.DataFrame(geo_pls["pairwise"], index=types, columns=types)
        pair.to_csv(GEO_DIR / "partC_pairwise_angles_pls.csv")
    if geo_ridge["pairwise"] is not None:
        pd.DataFrame(geo_ridge["pairwise"], index=types, columns=types).to_csv(
            GEO_DIR / "partC_pairwise_angles_ridge.csv")

    # ------------------------------------------------------------------ permutation null (full-data PLS-1; identity basis fixed)
    log(f"\n[C3] permutation null: {N_PERM_C} donor-level within-site age shuffles (PLS-1)")
    null_rows = []
    for i in range(N_PERM_C):
        obs_p = permute_age_within_site_obs(obs, rng)
        yp = obs_p.age.to_numpy(float)
        Wp, _, _ = directions_on_mask(X, yp, ct, types, np.ones(len(obs), bool), ridge=False)
        g = score_geometry(Wp, basis, types)
        null_rows.append(dict(
            perm=i,
            median_angle_deg=g["median_angle_deg"],
            median_proj_r2=g["median_proj_r2"],
            consensus_angle_deg=g["consensus_angle_deg"],
            consensus_proj_r2=g["consensus_proj_r2"],
            pairwise_median_deg=g["pairwise_median_deg"],
            shared_pc1_frac=g["shared_pc1_frac"],
        ))
        if (i + 1) % 20 == 0:
            log(f"   perm {i+1}/{N_PERM_C}")
    null_tab = pd.DataFrame(null_rows)
    null_tab.to_csv(GEO_DIR / "partC_perm_null.csv", index=False)

    def _p(obs_v, null_v, greater=True):
        null_v = np.asarray(null_v, float)
        null_v = null_v[np.isfinite(null_v)]
        if not np.isfinite(obs_v) or len(null_v) == 0:
            return np.nan
        if greater:
            return float((np.sum(null_v >= obs_v) + 1) / (len(null_v) + 1))
        return float((np.sum(null_v <= obs_v) + 1) / (len(null_v) + 1))

    c3 = dict(
        method="pls1",
        median_angle_deg=geo_pls["median_angle_deg"],
        median_proj_r2=geo_pls["median_proj_r2"],
        consensus_angle_deg=geo_pls["consensus_angle_deg"],
        consensus_proj_r2=geo_pls["consensus_proj_r2"],
        pairwise_median_deg=geo_pls["pairwise_median_deg"],
        shared_pc1_frac=geo_pls["shared_pc1_frac"],
        null_median_angle_mean=float(null_tab.median_angle_deg.mean()),
        null_median_proj_r2_mean=float(null_tab.median_proj_r2.mean()),
        null_consensus_angle_mean=float(null_tab.consensus_angle_deg.mean()),
        null_pairwise_median_mean=float(null_tab.pairwise_median_deg.mean()),
        null_shared_pc1_mean=float(null_tab.shared_pc1_frac.mean()),
        p_proj_r2_ge_obs=_p(geo_pls["median_proj_r2"], null_tab.median_proj_r2, greater=True),
        p_angle_le_obs=_p(geo_pls["median_angle_deg"], null_tab.median_angle_deg, greater=False),
        p_shared_pc1_ge_obs=_p(geo_pls["shared_pc1_frac"], null_tab.shared_pc1_frac, greater=True),
        n_perm=N_PERM_C,
        n_identity_axes=int(basis.shape[1]),
        p_genes=int(X.shape[1]),
    )
    # ridge as sensitivity
    c3["ridge_median_angle_deg"] = geo_ridge["median_angle_deg"]
    c3["ridge_median_proj_r2"] = geo_ridge["median_proj_r2"]
    c3["ridge_consensus_angle_deg"] = geo_ridge["consensus_angle_deg"]
    c3["ridge_shared_pc1_frac"] = geo_ridge["shared_pc1_frac"]
    log(f"[C3] obs median angle {c3['median_angle_deg']:.1f}° vs null {c3['null_median_angle_mean']:.1f}° "
        f"(p_angle≤obs={c3['p_angle_le_obs']:.3f})")
    log(f"     obs median proj R² {c3['median_proj_r2']:.3f} vs null {c3['null_median_proj_r2_mean']:.3f} "
        f"(p={c3['p_proj_r2_ge_obs']:.3f})")
    log(f"     shared PC1 {c3['shared_pc1_frac']:.3f} vs null {c3['null_shared_pc1_mean']:.3f} "
        f"(p={c3['p_shared_pc1_ge_obs']:.3f})")

    # ------------------------------------------------------------------ CV schemes
    rng_folds = np.random.default_rng(GEO_SEED)
    schemes = make_folds(obs, rng_folds)
    cv_rec = {}
    for name, folds in schemes.items():
        log(f"\n[C] CV scheme={name}: age directions on each training fold")
        fold_rows = []
        # z-score using training rows only, then identity basis on train
        W_acc_pls = []
        for fi, fold in enumerate(folds):
            tr = fold[0]
            Xtr, mu_tr, sd_tr = _zscore(Y[tr])
            basis_tr, _, _ = identity_basis(Xtr, ct[tr], types)
            ytr = y[tr]
            # map back: directions_on_mask expects full X aligned to obs; pass train-standardized via a slice
            Wp, Wr, meta_f = directions_on_mask(
                Xtr, ytr, ct[tr], types, np.ones(len(tr), bool), ridge=False)
            g = score_geometry(Wp, basis_tr, types)
            fold_rows.append(dict(
                fold=fi, tag=str(fold[2] if len(fold) > 2 else fi),
                median_angle_deg=g["median_angle_deg"],
                median_proj_r2=g["median_proj_r2"],
                consensus_angle_deg=g["consensus_angle_deg"],
                pairwise_median_deg=g["pairwise_median_deg"],
                shared_pc1_frac=g["shared_pc1_frac"],
            ))
            W_acc_pls.append(Wp)
        ft = pd.DataFrame(fold_rows)
        ft.to_csv(GEO_DIR / f"partC_cv_{name}_folds.csv", index=False)
        cv_rec[name] = dict(
            median_angle_deg=float(ft.median_angle_deg.median()),
            median_proj_r2=float(ft.median_proj_r2.median()),
            consensus_angle_deg=float(ft.consensus_angle_deg.median()),
            pairwise_median_deg=float(ft.pairwise_median_deg.median()),
            shared_pc1_frac=float(ft.shared_pc1_frac.median()),
            n_folds=len(ft),
        )
        log(f"   {name} median-over-folds angle={cv_rec[name]['median_angle_deg']:.1f}°  "
            f"proj R²={cv_rec[name]['median_proj_r2']:.3f}  "
            f"pairwise={cv_rec[name]['pairwise_median_deg']:.1f}°  "
            f"PC1={cv_rec[name]['shared_pc1_frac']:.3f}")

        # permutation null on CV: 20 draws, median-over-folds statistic
        n_cv_perm = 20
        log(f"   {name} permutation null ({n_cv_perm} draws of median-over-folds)...")
        cv_null = []
        for i in range(n_cv_perm):
            obs_p = permute_age_within_site_obs(obs, rng)
            yp = obs_p.age.to_numpy(float)
            angs, projs = [], []
            for fold in folds:
                tr = fold[0]
                Xtr, _, _ = _zscore(Y[tr])
                basis_tr, _, _ = identity_basis(Xtr, ct[tr], types)
                Wp, _, _ = directions_on_mask(
                    Xtr, yp[tr], ct[tr], types, np.ones(len(tr), bool), ridge=False)
                g = score_geometry(Wp, basis_tr, types)
                angs.append(g["median_angle_deg"])
                projs.append(g["median_proj_r2"])
            cv_null.append(dict(perm=i, median_angle_deg=float(np.nanmedian(angs)),
                                median_proj_r2=float(np.nanmedian(projs))))
        cvn = pd.DataFrame(cv_null)
        cvn.to_csv(GEO_DIR / f"partC_cv_{name}_perm.csv", index=False)
        cv_rec[name]["null_median_angle_mean"] = float(cvn.median_angle_deg.mean())
        cv_rec[name]["null_median_proj_r2_mean"] = float(cvn.median_proj_r2.mean())
        cv_rec[name]["p_angle_le_obs"] = _p(cv_rec[name]["median_angle_deg"],
                                            cvn.median_angle_deg, greater=False)
        cv_rec[name]["p_proj_r2_ge_obs"] = _p(cv_rec[name]["median_proj_r2"],
                                              cvn.median_proj_r2, greater=True)
        cv_rec[name]["n_perm"] = n_cv_perm
        log(f"   {name} null angle={cv_rec[name]['null_median_angle_mean']:.1f}°  "
            f"p={cv_rec[name]['p_angle_le_obs']:.3f}")

    # ------------------------------------------------------------------ interpretation
    # Orthogonal: high angle (~90), proj R² near the null (which for a random vector in high-p
    # into a 19-dim subspace is ~19/p ≈ 0). Fusion: proj R² >> null, angle << 90.
    p = X.shape[1]
    k = basis.shape[1]
    random_proj = k / p
    obs_r = c3["median_proj_r2"]
    null_r = c3["null_median_proj_r2_mean"]
    # "substantially within" if observed proj is much larger than null and angle well below 90
    if np.isfinite(obs_r) and obs_r > max(5 * max(null_r, random_proj), 0.20) and c3["median_angle_deg"] < 60:
        geom = "fused"
        geom_text = (
            "the age direction lies substantially WITHIN the identity subspace — fusion is geometric "
            "and the earlier verdict is confirmed at a deeper level."
        )
    elif np.isfinite(obs_r) and (obs_r <= 2 * max(null_r, random_proj) or c3["median_angle_deg"] >= 75):
        geom = "orthogonal"
        geom_text = (
            "the age direction is largely ORTHOGONAL to identity directions — the axes are separable "
            "in geometry even though genes cannot be partitioned."
        )
    else:
        geom = "partial"
        geom_text = (
            "partial overlap: the age direction is neither clearly orthogonal nor clearly inside the "
            "identity subspace. Report the angle and projection R²; do not overclaim fusion or separability."
        )
    log(f"\n[C3 claim] {geom}: {geom_text}")
    log(f"   (random vector into {k}-dim subspace of R^{p} has E[proj R²]={random_proj:.4f})")

    shared = (
        "shared across cell types" if (c3["shared_pc1_frac"] > 0.35
                                       and c3["p_shared_pc1_ge_obs"] < 0.05)
        else "largely type-specific"
    )
    log(f"[C4] per-type age directions are {shared} "
        f"(pairwise median {c3['pairwise_median_deg']:.1f}°; PC1 {c3['shared_pc1_frac']:.3f}; "
        f"null PC1 {c3['null_shared_pc1_mean']:.3f})")

    # ------------------------------------------------------------------ figures
    fig, axes = plt.subplots(1, 2, figsize=(9.6, 3.8))
    ax = axes[0]
    ax.hist(null_tab.median_angle_deg, bins=20, color="0.75", label="perm null")
    ax.axvline(c3["median_angle_deg"], c="tab:red", lw=2, label=f"obs {c3['median_angle_deg']:.1f}°")
    ax.axvline(c3["null_median_angle_mean"], c="k", ls=":", label="null mean")
    ax.set_xlabel("median-over-types angle to identity subspace (deg)")
    ax.set_ylabel("permutations")
    ax.set_title("C3  PLS-1 age vs identity")
    ax.legend(fontsize=7)
    ax = axes[1]
    ax.hist(null_tab.median_proj_r2, bins=20, color="0.75", label="perm null")
    ax.axvline(c3["median_proj_r2"], c="tab:red", lw=2, label=f"obs {c3['median_proj_r2']:.3f}")
    ax.axvline(c3["null_median_proj_r2_mean"], c="k", ls=":")
    ax.axvline(random_proj, c="tab:blue", ls="--", lw=1, label=f"k/p={random_proj:.4f}")
    ax.set_xlabel("median-over-types proj R² on identity subspace")
    ax.set_title("C3  variance of age dir. in identity subspace")
    ax.legend(fontsize=7)
    fig.tight_layout()
    fig.savefig(GEO_FIG / "partC_angle_vs_null.png", dpi=140)
    plt.close(fig)

    if geo_pls["pairwise"] is not None:
        fig, ax = plt.subplots(figsize=(7.6, 6.6))
        im = ax.imshow(geo_pls["pairwise"], cmap="viridis", vmin=0, vmax=90)
        ax.set_xticks(range(len(types)))
        ax.set_yticks(range(len(types)))
        short = [t if len(t) < 28 else t[:25] + "…" for t in types]
        ax.set_xticklabels(short, rotation=90, fontsize=6)
        ax.set_yticklabels(short, fontsize=6)
        fig.colorbar(im, ax=ax, fraction=0.046, pad=0.04, label="unsigned angle (deg)")
        ax.set_title(f"C4  pairwise angles of per-type PLS-1 age directions\n"
                     f"median {c3['pairwise_median_deg']:.1f}°  PC1 {c3['shared_pc1_frac']:.2f}")
        fig.tight_layout()
        fig.savefig(GEO_FIG / "partC_pairwise_angles.png", dpi=140)
        plt.close(fig)

    fig, ax = plt.subplots(figsize=(7.2, 4.4))
    pt = geo_pls["per_type"].sort_values("angle_deg")
    y_pos = np.arange(len(pt))
    ax.barh(y_pos, pt.angle_deg, color="tab:red", height=0.7)
    ax.axvline(c3["null_median_angle_mean"], c="k", ls=":", label="perm-null median angle")
    ax.axvline(90, c="0.7", lw=0.6)
    ax.set_yticks(y_pos)
    ax.set_yticklabels(pt.celltype, fontsize=6)
    ax.set_xlabel("angle of age direction to identity subspace (deg)")
    ax.set_title("C3  per-type PLS-1")
    ax.legend(fontsize=7)
    fig.tight_layout()
    fig.savefig(GEO_FIG / "partC_per_type_angles.png", dpi=140)
    plt.close(fig)

    fig, ax = plt.subplots(figsize=(6.4, 3.6))
    ax.hist(null_tab.shared_pc1_frac, bins=20, color="0.75", label="perm null")
    ax.axvline(c3["shared_pc1_frac"], c="tab:red", lw=2, label=f"obs {c3['shared_pc1_frac']:.3f}")
    ax.axvline(c3["null_shared_pc1_mean"], c="k", ls=":")
    ax.set_xlabel("PC1 fraction of per-type age directions")
    ax.set_title("C4  is the age direction shared?")
    ax.legend(fontsize=7)
    fig.tight_layout()
    fig.savefig(GEO_FIG / "partC_shared_pc1.png", dpi=140)
    plt.close(fig)

    summary = dict(
        seed=GEO_SEED, dataset_id=DATASET_ID,
        n_genes=int(Y.shape[1]), n_types=len(types), n_identity_axes=int(basis.shape[1]),
        random_proj_expect=float(random_proj),
        c3=c3, c3_ridge=dict(
            median_angle_deg=geo_ridge["median_angle_deg"],
            median_proj_r2=geo_ridge["median_proj_r2"],
            consensus_angle_deg=geo_ridge["consensus_angle_deg"],
            shared_pc1_frac=geo_ridge["shared_pc1_frac"],
            pairwise_median_deg=geo_ridge["pairwise_median_deg"],
        ),
        cv=cv_rec,
        geom=geom, geom_text=geom_text, shared=shared,
        n_perm_c=N_PERM_C,
    )
    dump_json(GEO_DIR / "partC_summary.json", summary)
    (GEO_DIR / "partC_claim.txt").write_text(
        f"{geom}\n{geom_text}\nshared={shared}\n"
        f"angle={c3['median_angle_deg']:.2f} null={c3['null_median_angle_mean']:.2f}\n"
        f"proj_r2={c3['median_proj_r2']:.4f} null={c3['null_median_proj_r2_mean']:.4f}\n",
        encoding="utf-8",
    )
    log("\n[C] wrote results/geometry/partC_*")
    log("[C] done.")
    return summary


if __name__ == "__main__":
    run()
