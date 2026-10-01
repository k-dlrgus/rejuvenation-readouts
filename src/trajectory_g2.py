"""STAGE G2 — orthogonality, properly tested.

G2a. Angle of the age direction to the identity subspace; directions and the
     identity basis fit inside training folds only.
G2b. Can the identity subspace predict age? Age R² from identity coordinates
     vs from the age direction (same CV, same metric, same null).
G2c. Residual test: regress out the identity subspace, refit the age direction
     on the residual. Does age R² survive?

Halt before G3 if G2b shows redundancy or G2c shows collapse.

Usage: python src/trajectory_g2.py
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
from trajectory_common import (  # noqa: E402
    TRAJ_DIR, TRAJ_FIG, TRAJ_SEED, DATASET_ID, N_PERM_G2, ALPHAS,
    G2B_RATIO, G2B_DELTA, G2C_RATIO, STOP_R2_COLLAPSE,
    Logger, dump_json, load_json, traj_log_banner, load_phase1_matrix, make_folds,
    permute_age_within_site_obs, score_age_predictions, primary_r2, primary_mae,
    scores_json, pls1_predict, ridge_predict_dir, ridge_from_svd_cache, svd_train,
    identity_basis, identity_project, residualize_identity, zscore_train,
    permutation_p, summarize_null_col, proj_r2, angle_to_subspace_deg,
    fmt,
)
from geometry_partC import score_geometry  # noqa: E402
from trajectory_g1 import _age_oof  # noqa: E402


def _fold_pack(Y, tr, te, ct, types):
    Xtr, Xte, mu, sd = zscore_train(Y[tr], Y[te])
    basis, S_id, used = identity_basis(Xtr, ct[tr], types)
    return dict(Xtr=Xtr, Xte=Xte, mu=mu, sd=sd, basis=basis, S_id=S_id, used=used)


def _age_dirs_on_train(Xtr, ytr, ct_tr, types, method, svd_cache=None):
    p = Xtr.shape[1]
    W = []
    built = {} if svd_cache is None and method == "ridge" else None
    for t in types:
        m = ct_tr == t
        if int(m.sum()) < 8:
            W.append(np.full(p, np.nan))
            continue
        if method == "pls1":
            _, w, _, _ = pls1_predict(Xtr[m], ytr[m], Xtr[m])
        else:
            if svd_cache is not None and t in svd_cache:
                U, S, Vt = svd_cache[t]
                _, coef, _, _ = ridge_from_svd_cache(U, S, Vt, ytr[m], Xtr[m], alphas=ALPHAS)
                nrm = float(np.linalg.norm(coef))
                w = coef / nrm if nrm > 1e-12 else coef
            else:
                _, w, _, _ = ridge_predict_dir(Xtr[m], ytr[m], Xtr[m])
                if built is not None:
                    built[t] = svd_train(Xtr[m])
        W.append(w)
    return np.column_stack(W), built


def _g2a_scheme(Y, obs, folds, types, method, rng, n_perm, log, name):
    ct = obs.celltype.astype(str).to_numpy()
    y = obs.age.to_numpy(float)
    fold_rows = []
    fold_packs = []
    fold_svds = []
    for fi, fold in enumerate(folds):
        tr, te = fold[0], fold[1]
        pack = _fold_pack(Y, tr, te, ct, types)
        W, built = _age_dirs_on_train(pack["Xtr"], y[tr], ct[tr], types, method)
        fold_packs.append(pack)
        fold_svds.append(built or {})
        g = score_geometry(W, pack["basis"], types)
        fold_rows.append(dict(
            fold=fi, tag=str(fold[2] if len(fold) > 2 else fi),
            median_angle_deg=g["median_angle_deg"],
            median_proj_r2=g["median_proj_r2"],
            consensus_angle_deg=g["consensus_angle_deg"],
            pairwise_median_deg=g["pairwise_median_deg"],
            n_axes=int(pack["basis"].shape[1]),
        ))
    ft = pd.DataFrame(fold_rows)
    obs_ang = float(ft.median_angle_deg.median())
    obs_pr = float(ft.median_proj_r2.median())
    log(f"   {name} {method}  median-over-folds angle={obs_ang:.1f}°  proj R²={obs_pr:.3f}")

    log(f"   {name} {method} permutation null ({n_perm} draws)...")
    null_rows = []
    for i in range(n_perm):
        obs_p = permute_age_within_site_obs(obs, rng)
        yp = obs_p.age.to_numpy(float)
        angs, projs = [], []
        for fi, fold in enumerate(folds):
            tr = fold[0]
            pack = fold_packs[fi]
            W, _ = _age_dirs_on_train(
                pack["Xtr"], yp[tr], ct[tr], types, method, svd_cache=fold_svds[fi] or None)
            g = score_geometry(W, pack["basis"], types)
            angs.append(g["median_angle_deg"])
            projs.append(g["median_proj_r2"])
        null_rows.append(dict(perm=i, median_angle_deg=float(np.nanmedian(angs)),
                              median_proj_r2=float(np.nanmedian(projs))))
        if (i + 1) % 10 == 0:
            log(f"      perm {i+1}/{n_perm}")
    nt = pd.DataFrame(null_rows)
    p_ang = permutation_p(obs_ang, nt.median_angle_deg, greater=False)
    p_pr = permutation_p(obs_pr, nt.median_proj_r2, greater=True)
    log(f"   {name} {method} null angle={nt.median_angle_deg.mean():.1f}°  "
        f"p_angle={p_ang:.3f}  p_proj={p_pr:.3f}")
    return dict(
        median_angle_deg=obs_ang, median_proj_r2=obs_pr,
        null_angle_mean=float(nt.median_angle_deg.mean()),
        null_proj_r2_mean=float(nt.median_proj_r2.mean()),
        p_angle=p_ang, p_proj=p_pr, n_perm=n_perm,
        folds=ft, null=nt,
    )


def _predict_from_Z(Ztr, ytr, Zte, method="ridge"):
    """Age from identity coordinates (k << n). All k axes used; no gene selection."""
    Ztr = np.asarray(Ztr, np.float64)
    Zte = np.asarray(Zte, np.float64)
    ytr = np.asarray(ytr, np.float64)
    if Ztr.shape[1] == 0 or len(ytr) < 8:
        return np.full(len(Zte), float(np.mean(ytr)) if len(ytr) else np.nan)
    if method == "pls1":
        p, _, _, _ = pls1_predict(Ztr, ytr, Zte)
        return p
    p, _, _, _ = ridge_predict_dir(Ztr, ytr, Zte)
    return p


def _g2b_oof(Y, obs, folds, y_override=None):
    """Within-type age from identity-subspace coordinates (train basis, train z-score)."""
    pred = np.full(len(obs), np.nan)
    y = obs.age.to_numpy(float) if y_override is None else np.asarray(y_override, float)
    ct = obs.celltype.astype(str).to_numpy()
    types = sorted(pd.unique(ct))
    for fold in folds:
        tr, te = fold[0], fold[1]
        pack = _fold_pack(Y, tr, te, ct, types)
        Ztr = identity_project(pack["Xtr"], pack["basis"])
        Zte = identity_project(pack["Xte"], pack["basis"])
        for t in types:
            tr_t = ct[tr] == t
            te_t = ct[te] == t
            if int(tr_t.sum()) < 8 or int(te_t.sum()) < 1:
                continue
            pred[te[te_t]] = _predict_from_Z(Ztr[tr_t], y[tr][tr_t], Zte[te_t])
    return pred


def _g2c_oof(Y, obs, folds, method, y_override=None, svd_cache=None):
    """Age direction fit on identity-residual expression, train-only projector."""
    pred = np.full(len(obs), np.nan)
    y = obs.age.to_numpy(float) if y_override is None else np.asarray(y_override, float)
    ct = obs.celltype.astype(str).to_numpy()
    types = sorted(pd.unique(ct))
    new_cache = {} if svd_cache is None and method == "ridge" else None
    for fi, fold in enumerate(folds):
        tr, te = fold[0], fold[1]
        pack = _fold_pack(Y, tr, te, ct, types)
        Rtr = residualize_identity(pack["Xtr"], pack["basis"])
        Rte = residualize_identity(pack["Xte"], pack["basis"])
        built = {}
        fold_cache = (svd_cache or {}).get(fi, {})
        for t in types:
            tr_t = ct[tr] == t
            te_t = ct[te] == t
            if int(tr_t.sum()) < 8 or int(te_t.sum()) < 1:
                continue
            Xtr_t, Xte_t = Rtr[tr_t], Rte[te_t]
            ytr = y[tr][tr_t]
            if method == "pls1":
                p, _, _, _ = pls1_predict(Xtr_t, ytr, Xte_t)
            else:
                if svd_cache is not None and t in fold_cache:
                    U, S, Vt = fold_cache[t]
                    p, _, _, _ = ridge_from_svd_cache(U, S, Vt, ytr, Xte_t, alphas=ALPHAS)
                else:
                    p, _, _, _ = ridge_predict_dir(Xtr_t, ytr, Xte_t)
                    if new_cache is not None:
                        built[t] = svd_train(Xtr_t)
            pred[te[te_t]] = p
        if new_cache is not None:
            new_cache[fi] = built
    return pred, new_cache


def _null_from_pred_fn(obs, folds, rng, n_perm, pred_fn, log_every=None, log=None):
    rows = []
    for i in range(n_perm):
        obs_p = permute_age_within_site_obs(obs, rng)
        yp = obs_p.age.to_numpy(float)
        pred = pred_fn(yp)
        sc = score_age_predictions(obs_p, pred, y=yp)
        rows.append(dict(perm=i, pooled=sc["pooled_median_type_r2"],
                         site_strat=sc["site_strat_median_type_r2"]))
        if log is not None and log_every and (i + 1) % log_every == 0:
            log(f"      perm {i+1}/{n_perm}")
    return pd.DataFrame(rows)


def run():
    rng = np.random.default_rng(TRAJ_SEED)
    log = Logger(TRAJ_DIR / "g2_report.txt")
    try:
        return _run(log, rng)
    finally:
        log.close()


def _run(log, rng):
    traj_log_banner(log, "G2 — orthogonality controls")
    g1_path = TRAJ_DIR / "g1_summary.json"
    if not g1_path.exists():
        raise FileNotFoundError("results/trajectory/g1_summary.json missing — run G1 first")
    g1 = load_json(g1_path)
    if g1.get("stop_G1"):
        log("[G2] STOP G1 is set — refusing to run G2.")
        dump_json(TRAJ_DIR / "g2_summary.json", dict(skipped=True, reason="STOP G1"))
        (TRAJ_DIR / "g2_STOP.txt").write_text("stop=True\nreason=STOP G1\n", encoding="utf-8")
        return dict(skipped=True, reason="STOP G1")

    data = load_phase1_matrix(log)
    Y, obs = data["Y"], data["obs"]
    p = Y.shape[1]
    log(f"[filter] G2 uses ALL {p:,} genes; no threshold.")
    ct = obs.celltype.astype(str).to_numpy()
    types = sorted(pd.unique(ct))
    rng_folds = np.random.default_rng(TRAJ_SEED)
    schemes = make_folds(obs, rng_folds)

    # ------------------------------------------------------------------ G2a
    g2a = {}
    for name, folds in schemes.items():
        g2a[name] = {}
        for method in ("pls1", "ridge"):
            log(f"\n[G2a] scheme={name} method={method}")
            rec = _g2a_scheme(Y, obs, folds, types, method, rng, N_PERM_G2, log, name)
            rec["folds"].to_csv(TRAJ_DIR / f"g2a_folds_{name}_{method}.csv", index=False)
            rec["null"].to_csv(TRAJ_DIR / f"g2a_null_{name}_{method}.csv", index=False)
            g2a[name][method] = {k: v for k, v in rec.items() if k not in ("folds", "null")}

    # ------------------------------------------------------------------ G2b / G2c  (ridge primary; PLS-1 if G1 said they differ)
    methods = ("pls1", "ridge")
    g2b, g2c = {}, {}
    for name, folds in schemes.items():
        log(f"\n[G2b] identity-subspace coordinates → age   scheme={name}")
        pred_b = _g2b_oof(Y, obs, folds)
        sc_b = score_age_predictions(obs, pred_b)
        r2_b = primary_r2(name, sc_b)
        log(f"   identity-coords age R²={r2_b:+.3f}  "
            f"(pooled={sc_b['pooled_median_type_r2']:+.3f}  "
            f"site_strat={sc_b['site_strat_median_type_r2']:+.3f})")
        log(f"   permutation null ({N_PERM_G2} draws)...")
        nb = _null_from_pred_fn(
            obs, folds, rng, N_PERM_G2,
            lambda yp, folds=folds: _g2b_oof(Y, obs, folds, y_override=yp),
            log_every=10, log=log,
        )
        nb.to_csv(TRAJ_DIR / f"g2b_null_{name}.csv", index=False)
        key = "site_strat" if name == "within_site" else "pooled"
        null_b = summarize_null_col(nb[key])
        p_b = permutation_p(r2_b, nb[key], greater=True)
        log(f"   null {key} mean={null_b['mean']:+.3f}  p={p_b:.3f}")
        sc_b["per_type"].assign(scheme=name).to_csv(
            TRAJ_DIR / f"g2b_per_type_{name}.csv", index=False)

        age_dir = {}
        resid = {}
        for method in methods:
            # age-direction R² from G1 (same folds, same metric) — recompute here so G2 is self-contained
            log(f"[G2b/c] age-direction {method} scheme={name} (recompute, train-only)")
            out_a = _age_oof(Y, obs, folds, method)
            sc_a = score_age_predictions(obs, out_a["pred"])
            r2_a = primary_r2(name, sc_a)
            log(f"   age-dir {method} R²={r2_a:+.3f}")
            age_dir[method] = dict(r2=r2_a, mae=primary_mae(name, sc_a),
                                   scores=scores_json(sc_a))

            log(f"[G2c] residual age direction {method} scheme={name}")
            pred_c, cache = _g2c_oof(Y, obs, folds, method)
            sc_c = score_age_predictions(obs, pred_c)
            r2_c = primary_r2(name, sc_c)
            log(f"   residual {method} R²={r2_c:+.3f}")
            log(f"   permutation null ({N_PERM_G2} draws)...")
            nc = _null_from_pred_fn(
                obs, folds, rng, N_PERM_G2,
                lambda yp, folds=folds, method=method, cache=cache: _g2c_oof(
                    Y, obs, folds, method, y_override=yp, svd_cache=cache)[0],
                log_every=10, log=log,
            )
            nc.to_csv(TRAJ_DIR / f"g2c_null_{name}_{method}.csv", index=False)
            null_c = summarize_null_col(nc[key])
            p_c = permutation_p(r2_c, nc[key], greater=True)
            log(f"   residual null {key} mean={null_c['mean']:+.3f}  p={p_c:.3f}")
            sc_c["per_type"].assign(scheme=name, method=method).to_csv(
                TRAJ_DIR / f"g2c_per_type_{name}_{method}.csv", index=False)
            resid[method] = dict(
                r2=r2_c, mae=primary_mae(name, sc_c),
                scores=scores_json(sc_c), null=null_c, p=p_c, n_perm=N_PERM_G2,
            )

        # redundancy / collapse vs primary ridge age-dir
        r2_age = age_dir["ridge"]["r2"]
        r2_id = r2_b
        r2_res = resid["ridge"]["r2"]
        redundant = (
            np.isfinite(r2_age) and np.isfinite(r2_id) and r2_age > STOP_R2_COLLAPSE and (
                (r2_id >= G2B_RATIO * r2_age) or ((r2_age - r2_id) < G2B_DELTA)
            )
        )
        collapsed = (
            (not np.isfinite(r2_res) or r2_res <= STOP_R2_COLLAPSE) or (
                np.isfinite(r2_age) and r2_age > STOP_R2_COLLAPSE and r2_res < G2C_RATIO * r2_age
            )
        )
        log(f"   G2b ridge: R²_id={r2_id:+.3f} vs R²_age={r2_age:+.3f}  "
            f"ratio={r2_id / r2_age if r2_age else np.nan:.3f}  redundant={redundant}")
        log(f"   G2c ridge: R²_resid={r2_res:+.3f} vs R²_age={r2_age:+.3f}  collapsed={collapsed}")

        g2b[name] = dict(
            r2_identity_coords=r2_id, mae=primary_mae(name, sc_b),
            scores=scores_json(sc_b), null=null_b, p=p_b, n_perm=N_PERM_G2,
            age_dir=age_dir, redundant=bool(redundant),
            ratio=(float(r2_id / r2_age) if (r2_age and np.isfinite(r2_age) and r2_age != 0) else None),
            delta=float(r2_age - r2_id) if (np.isfinite(r2_age) and np.isfinite(r2_id)) else None,
        )
        g2c[name] = dict(residual=resid, collapsed=bool(collapsed))

    # Halt rule: G2b redundancy OR G2c collapse on the PRIMARY scheme (within_site),
    # with LOSO reported. If within_site is clean and LOSO is not, report both and proceed
    # (LOSO is 1 df of bank transfer). Halt if within_site fails.
    redundant_ws = bool(g2b["within_site"]["redundant"])
    collapsed_ws = bool(g2c["within_site"]["collapsed"])
    stop = bool(redundant_ws or collapsed_ws)
    if redundant_ws and collapsed_ws:
        reason = "G2b redundant AND G2c collapsed on within_site — scores are not independent"
    elif redundant_ws:
        reason = "G2b: identity coordinates predict age nearly as well as the age direction (within_site)"
    elif collapsed_ws:
        reason = "G2c: age R² collapses after removing the identity subspace (within_site)"
    else:
        reason = ""
    log(f"\n[STOP G2] stop={stop}  redundant_ws={redundant_ws}  collapsed_ws={collapsed_ws}")
    if reason:
        log(f"   {reason}")
        log("   Halt before G3 — a trajectory on redundant axes is a figure that means nothing.")
    else:
        log("   G2b not redundant; G2c age signal survives. Proceed to G3.")
    (TRAJ_DIR / "g2_STOP.txt").write_text(
        f"stop={stop}\nredundant_ws={redundant_ws}\ncollapsed_ws={collapsed_ws}\nreason={reason}\n",
        encoding="utf-8",
    )

    _figures(g2a, g2b, g2c, log)

    summary = dict(
        seed=TRAJ_SEED, dataset_id=DATASET_ID, n_genes=int(p), n_perm=N_PERM_G2,
        g2a=g2a, g2b=g2b,         g2c={
            s: dict(
                collapsed=g2c[s]["collapsed"],
                residual={m: dict(
                    r2=g2c[s]["residual"][m]["r2"],
                    mae=g2c[s]["residual"][m]["mae"],
                    null=g2c[s]["residual"][m]["null"],
                    p=g2c[s]["residual"][m]["p"],
                    n_perm=g2c[s]["residual"][m]["n_perm"],
                    scores=g2c[s]["residual"][m]["scores"],
                ) for m in g2c[s]["residual"]},
            ) for s in g2c
        },
        stop_G2=stop, stop_reason=reason,
        redundant_within_site=redundant_ws, collapsed_within_site=collapsed_ws,
        thresholds=dict(g2b_ratio=G2B_RATIO, g2b_delta=G2B_DELTA, g2c_ratio=G2C_RATIO),
        note="Directions and identity basis fit on training folds only. All genes.",
    )
    dump_json(TRAJ_DIR / "g2_summary.json", summary)
    _write_table(g2a, g2b, g2c)
    log("[G2] wrote results/trajectory/g2_*")
    log("[G2] done.")
    return summary


def _write_table(g2a, g2b, g2c):
    rows = []
    for scheme in g2a:
        for method, rec in g2a[scheme].items():
            rows.append(dict(control="G2a_angle", scheme=scheme, method=method,
                             value=rec["median_angle_deg"], null=rec["null_angle_mean"],
                             p=rec["p_angle"]))
            rows.append(dict(control="G2a_proj_r2", scheme=scheme, method=method,
                             value=rec["median_proj_r2"], null=rec["null_proj_r2_mean"],
                             p=rec["p_proj"]))
        rows.append(dict(control="G2b_id_coords_age_r2", scheme=scheme, method="ridge_on_Z",
                         value=g2b[scheme]["r2_identity_coords"],
                         null=g2b[scheme]["null"]["mean"], p=g2b[scheme]["p"]))
        for method in ("pls1", "ridge"):
            rows.append(dict(control="G2b_age_dir_r2", scheme=scheme, method=method,
                             value=g2b[scheme]["age_dir"][method]["r2"],
                             null=np.nan, p=np.nan))
            rc = g2c[scheme]["residual"][method]
            rows.append(dict(control="G2c_residual_age_r2", scheme=scheme, method=method,
                             value=rc["r2"], null=rc["null"]["mean"], p=rc["p"]))
    pd.DataFrame(rows).to_csv(TRAJ_DIR / "g2_summary_table.csv", index=False)


def _figures(g2a, g2b, g2c, log):
    # G2a null histograms (within_site PLS-1)
    rec = g2a.get("within_site", {}).get("pls1")
    path = TRAJ_DIR / "g2a_null_within_site_pls1.csv"
    if path.exists() and rec:
        nt = pd.read_csv(path)
        fig, axes = plt.subplots(1, 2, figsize=(9.6, 3.6))
        ax = axes[0]
        ax.hist(nt.median_angle_deg, bins=15, color="0.75")
        ax.axvline(rec["median_angle_deg"], c="tab:red", lw=2,
                   label=f"obs {rec['median_angle_deg']:.1f}°")
        ax.axvline(rec["null_angle_mean"], c="k", ls=":")
        ax.set_xlabel("median-over-folds angle (deg)")
        ax.set_title("G2a  PLS-1 angle (within-site CV)")
        ax.legend(fontsize=7)
        ax = axes[1]
        ax.hist(nt.median_proj_r2, bins=15, color="0.75")
        ax.axvline(rec["median_proj_r2"], c="tab:red", lw=2,
                   label=f"obs {rec['median_proj_r2']:.3f}")
        ax.axvline(rec["null_proj_r2_mean"], c="k", ls=":")
        ax.set_xlabel("median proj R²")
        ax.set_title("G2a  projection onto identity")
        ax.legend(fontsize=7)
        fig.tight_layout()
        fig.savefig(TRAJ_FIG / "g2a_angle_vs_null.png", dpi=140)
        plt.close(fig)

    fig, axes = plt.subplots(1, 2, figsize=(9.8, 4.0))
    for ax, scheme in zip(axes, ("within_site", "loso")):
        r_age = g2b[scheme]["age_dir"]["ridge"]["r2"]
        r_pls = g2b[scheme]["age_dir"]["pls1"]["r2"]
        r_id = g2b[scheme]["r2_identity_coords"]
        r_res = g2c[scheme]["residual"]["ridge"]["r2"]
        r_res_p = g2c[scheme]["residual"]["pls1"]["r2"]
        vals = [r_age, r_pls, r_id, r_res, r_res_p]
        labs = ["age-dir ridge", "age-dir PLS-1", "identity coords",
                "residual ridge", "residual PLS-1"]
        colors = ["tab:red", "tab:orange", "tab:blue", "tab:green", "tab:cyan"]
        ax.bar(np.arange(len(vals)), vals, color=colors)
        ax.axhline(0, c="0.4", lw=0.6)
        ax.set_xticks(np.arange(len(vals)))
        ax.set_xticklabels(labs, rotation=25, ha="right", fontsize=7)
        ax.set_ylabel("primary age R²")
        ax.set_title(f"G2b / G2c  {scheme}")
    fig.tight_layout()
    fig.savefig(TRAJ_FIG / "g2b_g2c_r2.png", dpi=140)
    plt.close(fig)
    log("[G2] figures written")


if __name__ == "__main__":
    run()
