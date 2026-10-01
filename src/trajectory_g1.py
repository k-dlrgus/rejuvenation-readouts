"""STAGE G1 — geometric age and identity scores, no gene selection.

Age score: within-type PLS-1 (primary geometric direction) and ridge (predictive
sensitivity). Score = projection onto the direction (ridge: the linear predictor).
Identity score: cell-type centroid SVD subspace; position in that subspace plus
distance to the sample's own type centroid. Classification = nearest centroid
in the subspace.

Directions and the identity basis are fit on TRAINING FOLDS ONLY.
Both CV schemes; site-stratified R² for within-site; permutation null on every
predictive number.

STOP G1: corrected shuffle null > ~0.05, or age R² collapses under both schemes.

Usage: python src/trajectory_g1.py
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
    TRAJ_DIR, TRAJ_FIG, TRAJ_SEED, DATASET_ID, N_PERM_G1, N_PERM_ID,
    STOP_NULL, STOP_R2_COLLAPSE, ALPHAS,
    Logger, dump_json, traj_log_banner, load_phase1_matrix, make_folds,
    permute_age_within_site_obs, permute_celltype,
    score_age_predictions, primary_r2, primary_mae, scores_json,
    pls1_predict, ridge_predict_dir, ridge_from_svd_cache, svd_train,
    identity_basis, identity_project, nearest_centroid_predict,
    zscore_train, chance_identity, identity_metrics,
    permutation_p, summarize_null_col, fmt, fmt_u,
)


def _age_oof(Y, obs, folds, method, y_override=None, cache=None):
    """Within-type OOF age predictions. Global train z-score, then a per-type direction.

    cache: optional dict fold_i -> {t -> (U,S,Vt, tr_local, te_local)} built on true y's X
    (X does not depend on y). Used to speed ridge permutations.
    """
    pred = np.full(len(obs), np.nan)
    raw = np.full(len(obs), np.nan)
    y = obs.age.to_numpy(float) if y_override is None else np.asarray(y_override, float)
    ct = obs.celltype.astype(str).to_numpy()
    types = sorted(pd.unique(ct))
    alphas = []
    new_cache = {} if cache is None and method == "ridge" else None
    for fi, fold in enumerate(folds):
        tr, te = fold[0], fold[1]
        Xtr, Xte, _, _ = zscore_train(Y[tr], Y[te])
        ct_tr, ct_te = ct[tr], ct[te]
        y_tr = y[tr]
        fold_cache = (cache or {}).get(fi, {}) if cache is not None else {}
        built = {}
        for t in types:
            tr_m = ct_tr == t
            te_m = ct_te == t
            if int(tr_m.sum()) < 8 or int(te_m.sum()) < 1:
                continue
            Xtr_t, ytr_t = Xtr[tr_m], y_tr[tr_m]
            Xte_t = Xte[te_m]
            if method == "pls1":
                p, _, s, _ = pls1_predict(Xtr_t, ytr_t, Xte_t)
                pred[te[te_m]] = p
                raw[te[te_m]] = s
            elif method == "ridge":
                if fi in (cache or {}) and t in fold_cache:
                    U, S, Vt = fold_cache[t]
                    p, coef, ymu, a = ridge_from_svd_cache(U, S, Vt, ytr_t, Xte_t, alphas=ALPHAS)
                else:
                    # First pass: fit (also caches SVD of already-zscored Xtr_t)
                    p, w, s, extra = ridge_predict_dir(Xtr_t, ytr_t, Xte_t)
                    a = extra[1]
                    if new_cache is not None:
                        U, S, Vt = svd_train(Xtr_t)
                        built[t] = (U, S, Vt)
                pred[te[te_m]] = p
                raw[te[te_m]] = p
                alphas.append(a)
            else:
                raise ValueError(method)
        if new_cache is not None:
            new_cache[fi] = built if built else {
                t: svd_train(Xtr[ct_tr == t])
                for t in types
                if int((ct_tr == t).sum()) >= 8 and int((ct_te == t).sum()) >= 1
            }
    out = dict(pred=pred, raw=raw, y=y,
               alpha_mean=float(np.nanmean(alphas) if alphas else np.nan))
    if new_cache is not None:
        out["svd_cache"] = new_cache
    return out


def _identity_oof(Y, obs, folds, ct_override=None):
    """Identity subspace on train (all genes), nearest-centroid on test projections."""
    ct = obs.celltype.astype(str).to_numpy() if ct_override is None else np.asarray(ct_override).astype(str)
    types = sorted(pd.unique(ct))
    pred = np.array([None] * len(obs), dtype=object)
    dist_own = np.full(len(obs), np.nan)
    n_axes = []
    for fold in folds:
        tr, te = fold[0], fold[1]
        Xtr, Xte, _, _ = zscore_train(Y[tr], Y[te])
        basis, S_id, used = identity_basis(Xtr, ct[tr], types)
        n_axes.append(int(basis.shape[1]))
        Ztr = identity_project(Xtr, basis)
        Zte = identity_project(Xte, basis)
        yhat, d_own, _, _ = nearest_centroid_predict(Ztr, ct[tr], Zte, yte=ct[te])
        pred[te] = yhat
        dist_own[te] = d_own
    chance = chance_identity(obs if ct_override is None else obs.assign(celltype=ct))
    met = identity_metrics(ct, np.array([str(p) if p is not None else "__none__" for p in pred]), chance)
    met.update(dist_own=dist_own, pred=pred, n_axes_mean=float(np.mean(n_axes) if n_axes else np.nan),
               y=ct)
    return met


def _shuffle_age_null(Y, obs, folds, method, rng, n_perm, svd_cache=None):
    rows = []
    for i in range(n_perm):
        obs_p = permute_age_within_site_obs(obs, rng)
        yp = obs_p.age.to_numpy(float)
        out = _age_oof(Y, obs, folds, method, y_override=yp, cache=svd_cache)
        sc = score_age_predictions(obs_p, out["pred"], y=yp)
        rows.append(dict(
            perm=i, method=method,
            pooled=sc["pooled_median_type_r2"],
            site_strat=sc["site_strat_median_type_r2"],
            pooled_mae=sc["pooled_median_type_mae"],
            site_strat_mae=primary_mae("within_site", sc),
        ))
    return pd.DataFrame(rows)


def _shuffle_ct_null(Y, obs, folds, rng, n_perm):
    rows = []
    for i in range(n_perm):
        obs_p = permute_celltype(obs, rng)
        met = _identity_oof(Y, obs, folds, ct_override=obs_p.celltype.astype(str).to_numpy())
        rows.append(dict(perm=i, acc=met["acc"], macro_f1=met["macro_f1"]))
    return pd.DataFrame(rows)


def _freeze_full(Y, obs, genes, log):
    """Full-cohort frozen scores for G3 (not used for G1 evaluation)."""
    ct = obs.celltype.astype(str).to_numpy()
    types = sorted(pd.unique(ct))
    y = obs.age.to_numpy(float)
    Xs, mu, sd = zscore_train(Y)
    basis, S_id, used = identity_basis(Xs, ct, types)
    Z = identity_project(Xs, basis)
    cents, used_c = [], []
    for t in used:
        m = ct == t
        cents.append(Z[m].mean(0))
        used_c.append(t)
    C = np.vstack(cents) if cents else np.zeros((0, basis.shape[1]))
    W_pls, W_ridge, intercepts, alphas, nrms = [], [], [], [], []
    meta_rows = []
    for t in types:
        m = ct == t
        n = int(m.sum())
        if n < 8:
            W_pls.append(np.full(Y.shape[1], np.nan))
            W_ridge.append(np.full(Y.shape[1], np.nan))
            intercepts.append(dict(pls=(np.nan, np.nan), ridge=np.nan))
            alphas.append(np.nan)
            nrms.append(np.nan)
            meta_rows.append(dict(celltype=t, n=n, ok=False))
            continue
        p_pls, w_pls, _, cal = pls1_predict(Xs[m], y[m], Xs[m])
        p_r, w_r, _, extra = ridge_predict_dir(Xs[m], y[m], Xs[m])
        W_pls.append(w_pls)
        W_ridge.append(w_r)
        intercepts.append(dict(pls=cal, ridge=extra[0]))
        alphas.append(extra[1])
        nrms.append(extra[2])
        meta_rows.append(dict(celltype=t, n=n, ok=True,
                              ridge_alpha=extra[1], ridge_norm=extra[2],
                              pls_intercept=cal[0], pls_slope=cal[1]))
    np.savez_compressed(
        TRAJ_DIR / "frozen_scores.npz",
        mu=mu, sd=sd, basis=basis, centroids=C,
        W_pls=np.column_stack(W_pls), W_ridge=np.column_stack(W_ridge),
        gene_id=genes.gene_id.astype(str).to_numpy(),
        symbol=genes.symbol.astype(str).to_numpy(),
        types=np.array(types, dtype=object),
        centroid_types=np.array(used_c, dtype=object),
        S_id=S_id,
    )
    dump_json(TRAJ_DIR / "frozen_meta.json", dict(
        n_genes=int(Y.shape[1]), n_types=len(types), n_identity_axes=int(basis.shape[1]),
        types=types, centroid_types=used_c,
        per_type=meta_rows, seed=TRAJ_SEED, dataset_id=DATASET_ID,
        note="Frozen on the full 233-donor Phase-1 cohort AFTER CV evaluation. Not used to score G1.",
    ))
    log(f"[freeze] wrote frozen_scores.npz  genes={Y.shape[1]}  types={len(types)}  "
        f"identity_axes={basis.shape[1]}")
    pd.DataFrame(meta_rows).to_csv(TRAJ_DIR / "g1_frozen_per_type.csv", index=False)


def run():
    rng = np.random.default_rng(TRAJ_SEED)
    log = Logger(TRAJ_DIR / "g1_report.txt")
    try:
        return _run(log, rng)
    finally:
        log.close()


def _run(log, rng):
    traj_log_banner(log, "G1 — geometric scores")
    data = load_phase1_matrix(log)
    Y, genes, obs = data["Y"], data["genes"], data["obs"]
    n, p = Y.shape
    log(f"[filter] G1 uses ALL {p:,} genes; gene_idx = arange({p}). No threshold, no A/I bucket.")
    log(f"[filter] cohort: {data['n_donors']} donors, {data['n_types']} cell types, "
        f"sites={data['site_donor_counts']}, rows={n}")
    ct = obs.celltype.astype(str).to_numpy()
    types = sorted(pd.unique(ct))
    chance = chance_identity(obs)
    log(f"[filter] identity chance = {chance:.4f}  (max of 1/{len(types)} and majority class)")

    rng_folds = np.random.default_rng(TRAJ_SEED)
    schemes = make_folds(obs, rng_folds)
    log(f"[folds] seed={TRAJ_SEED}  loso={len(schemes['loso'])}  "
        f"within_site={len(schemes['within_site'])}")

    methods = ("pls1", "ridge")
    age_rec = {}
    svd_caches = {}
    oof_frames = []

    for name, folds in schemes.items():
        log(f"\n[G1 age] scheme={name}")
        age_rec[name] = {}
        for method in methods:
            log(f"   fitting {method} on all {p} genes, training folds only...")
            out = _age_oof(Y, obs, folds, method)
            if method == "ridge" and "svd_cache" in out:
                svd_caches[(name, method)] = out["svd_cache"]
            sc = score_age_predictions(obs, out["pred"])
            r2 = primary_r2(name, sc)
            mae = primary_mae(name, sc)
            log(f"   {method}  primary R²={r2:+.3f}  MAE={mae:.2f} y  "
                f"(pooled R²={sc['pooled_median_type_r2']:+.3f}  "
                f"site_strat R²={sc['site_strat_median_type_r2']:+.3f})")
            pt = sc["per_type"].copy()
            pt["scheme"] = name
            pt["method"] = method
            pt.to_csv(TRAJ_DIR / f"g1_age_per_type_{name}_{method}.csv", index=False)
            oof = pd.DataFrame({
                "row": np.arange(len(obs)),
                "donor": obs.donor.astype(str).to_numpy(),
                "celltype": ct,
                "site": obs.Source.astype(str).to_numpy(),
                "age": obs.age.to_numpy(float),
                "pred": out["pred"],
                "scheme": name,
                "method": method,
            })
            oof_frames.append(oof)
            log(f"   permutation null ({N_PERM_G1} donor-level within-site shuffles)...")
            cache = svd_caches.get((name, method)) if method == "ridge" else None
            if method == "ridge" and cache is None:
                # build cache on a dedicated pass
                cache = _age_oof(Y, obs, folds, "ridge")["svd_cache"]
                svd_caches[(name, method)] = cache
            null_tab = _shuffle_age_null(Y, obs, folds, method, rng, N_PERM_G1, svd_cache=cache)
            null_tab.to_csv(TRAJ_DIR / f"g1_age_null_{name}_{method}.csv", index=False)
            key = "site_strat" if name == "within_site" else "pooled"
            null_s = summarize_null_col(null_tab[key])
            p_gt = permutation_p(r2, null_tab[key], greater=True)
            usable = bool(np.isfinite(null_s["mean"]) and null_s["mean"] <= STOP_NULL)
            log(f"   {method} null {key} mean={null_s['mean']:+.3f}  p95={null_s['p95']:+.3f}  "
                f"p(null>=obs)={p_gt:.3f}  usable={usable}")
            age_rec[name][method] = dict(
                r2=float(r2) if np.isfinite(r2) else None,
                mae=float(mae) if np.isfinite(mae) else None,
                pooled_r2=sc["pooled_median_type_r2"],
                site_strat_r2=sc["site_strat_median_type_r2"],
                pooled_mae=sc["pooled_median_type_mae"],
                scores=scores_json(sc),
                null=null_s,
                p_r2=p_gt,
                null_usable=usable,
                n_perm=N_PERM_G1,
                alpha_mean=out.get("alpha_mean"),
            )

        # identity
        log(f"[G1 identity] scheme={name}  nearest centroid in train identity subspace")
        idm = _identity_oof(Y, obs, folds)
        log(f"   acc={idm['acc']:.4f}  macro_F1={idm['macro_f1']:.4f}  "
            f"chance={idm['chance']:.4f}  mean_dist_own={np.nanmean(idm['dist_own']):.3f}  "
            f"n_axes={idm['n_axes_mean']:.1f}")
        log(f"   permutation null ({N_PERM_ID} cell-type shuffles)...")
        id_null = _shuffle_ct_null(Y, obs, folds, rng, N_PERM_ID)
        id_null.to_csv(TRAJ_DIR / f"g1_id_null_{name}.csv", index=False)
        id_null_s = summarize_null_col(id_null.acc)
        p_acc = permutation_p(idm["acc"], id_null.acc, greater=True)
        log(f"   acc null mean={id_null_s['mean']:.4f}  p={p_acc:.3f}")
        age_rec[name]["identity"] = dict(
            acc=idm["acc"], macro_f1=idm["macro_f1"], chance=idm["chance"],
            mean_dist_own=float(np.nanmean(idm["dist_own"])),
            n_axes_mean=idm["n_axes_mean"], n=idm["n"],
            null_acc=id_null_s, p_acc=p_acc, n_perm=N_PERM_ID,
        )
        pd.DataFrame({
            "donor": obs.donor.astype(str).to_numpy(),
            "celltype": ct,
            "site": obs.Source.astype(str).to_numpy(),
            "pred": [str(p) if p is not None else "" for p in idm["pred"]],
            "dist_own": idm["dist_own"],
            "scheme": name,
        }).to_csv(TRAJ_DIR / f"g1_id_oof_{name}.csv", index=False)

    if oof_frames:
        pd.concat(oof_frames, ignore_index=True).to_csv(TRAJ_DIR / "g1_age_oof.csv", index=False)

    # compare PLS vs ridge
    diffs = {}
    for name in schemes:
        r_pls = age_rec[name]["pls1"]["r2"]
        r_rid = age_rec[name]["ridge"]["r2"]
        d = abs((r_pls or np.nan) - (r_rid or np.nan))
        diffs[name] = dict(abs_delta=float(d) if np.isfinite(d) else None,
                           material=bool(np.isfinite(d) and d > 0.05))
        log(f"[G1] {name} |R²_ridge − R²_pls1|={d:.3f}  "
            f"material={'YES' if diffs[name]['material'] else 'no'}")

    # STOP G1
    # Null: within-site site-stratified shuffle mean for the primary age method.
    # Collapse: observed primary R² ~ 0 under BOTH schemes.
    primary_method = "ridge"  # predictive clock; PLS-1 is the geometric companion
    ws_null = age_rec["within_site"][primary_method]["null"]["mean"]
    loso_null = age_rec["loso"][primary_method]["null"]["mean"]
    ws_r2 = age_rec["within_site"][primary_method]["r2"] or np.nan
    loso_r2 = age_rec["loso"][primary_method]["r2"] or np.nan
    # Also check PLS-1; halt only if BOTH methods fail.
    ws_pls = age_rec["within_site"]["pls1"]["r2"] or np.nan
    loso_pls = age_rec["loso"]["pls1"]["r2"] or np.nan

    null_bad = (np.isfinite(ws_null) and ws_null > STOP_NULL)
    collapse = (
        (not np.isfinite(ws_r2) or ws_r2 <= STOP_R2_COLLAPSE) and
        (not np.isfinite(loso_r2) or loso_r2 <= STOP_R2_COLLAPSE) and
        (not np.isfinite(ws_pls) or ws_pls <= STOP_R2_COLLAPSE) and
        (not np.isfinite(loso_pls) or loso_pls <= STOP_R2_COLLAPSE)
    )
    stop = bool(null_bad or collapse)
    if null_bad:
        reason = (f"within-site {primary_method} shuffle null {ws_null:+.3f} > {STOP_NULL} "
                  "(corrected metric still unusable)")
    elif collapse:
        reason = (f"age R² collapsed under both CV schemes "
                  f"(ridge within_site={ws_r2:+.3f} loso={loso_r2:+.3f}; "
                  f"PLS-1 within_site={ws_pls:+.3f} loso={loso_pls:+.3f})")
    else:
        reason = ""
    log(f"\n[STOP G1] stop={stop}  null_bad={null_bad}  collapse={collapse}")
    if reason:
        log(f"   {reason}")
    else:
        log("   age direction produces usable held-out R²; null is not inflated. Proceed to G2.")
    (TRAJ_DIR / "g1_STOP.txt").write_text(
        f"stop={stop}\nnull_bad={null_bad}\ncollapse={collapse}\nreason={reason}\n",
        encoding="utf-8",
    )

    # freeze full-cohort scores only if we did not halt (G3 needs them; freeze is not evaluation)
    if not stop:
        _freeze_full(Y, obs, genes, log)
    else:
        log("[freeze] skipped because STOP G1 is set")

    # figures
    _figures(age_rec, schemes, log)

    summary = dict(
        seed=TRAJ_SEED, dataset_id=DATASET_ID,
        n_genes=int(p), n_donors=int(data["n_donors"]), n_types=len(types),
        n_rows=int(n), site_donor_counts=data["site_donor_counts"],
        methods=list(methods), primary_method=primary_method,
        chance_identity=chance,
        schemes=age_rec, pls_vs_ridge=diffs,
        stop_G1=stop, stop_reason=reason, null_bad=null_bad, collapse=collapse,
        n_perm_age=N_PERM_G1, n_perm_id=N_PERM_ID,
        note="No gene selection. Directions fit on training folds only.",
    )
    dump_json(TRAJ_DIR / "g1_summary.json", summary)
    _write_tables(age_rec)
    log("[G1] wrote results/trajectory/g1_*")
    log("[G1] done.")
    return summary


def _write_tables(age_rec):
    rows = []
    for scheme, rec in age_rec.items():
        for method in ("pls1", "ridge"):
            m = rec[method]
            rows.append(dict(
                scheme=scheme, model=f"age_{method}",
                r2=m["r2"], mae=m["mae"],
                pooled_r2=m["pooled_r2"], site_strat_r2=m["site_strat_r2"],
                null_mean=m["null"]["mean"], null_p95=m["null"]["p95"],
                p=m["p_r2"], null_usable=m["null_usable"], n_perm=m["n_perm"],
            ))
        idm = rec["identity"]
        rows.append(dict(
            scheme=scheme, model="identity_nearest_centroid",
            r2=idm["acc"], mae=idm["macro_f1"],
            pooled_r2=idm["acc"], site_strat_r2=np.nan,
            null_mean=idm["null_acc"]["mean"], null_p95=idm["null_acc"]["p95"],
            p=idm["p_acc"], null_usable=bool(idm["null_acc"]["mean"] <= 2 * idm["chance"]
                                            if np.isfinite(idm["null_acc"]["mean"]) else False),
            n_perm=idm["n_perm"],
        ))
    pd.DataFrame(rows).to_csv(TRAJ_DIR / "g1_summary_table.csv", index=False)


def _figures(age_rec, schemes, log):
    # per-type R² for within_site ridge vs PLS
    fig, axes = plt.subplots(1, 2, figsize=(11.2, 5.4), sharey=False)
    for ax, method in zip(axes, ("pls1", "ridge")):
        path = TRAJ_DIR / f"g1_age_per_type_within_site_{method}.csv"
        if not path.exists():
            continue
        pt = pd.read_csv(path)
        sub = pt[pt.site == "mean_of_sites"].sort_values("r2")
        y = np.arange(len(sub))
        ax.barh(y, sub.r2, color="tab:red" if method == "pls1" else "tab:blue", height=0.7)
        ax.axvline(0, c="0.5", lw=0.6)
        null_m = age_rec["within_site"][method]["null"]["mean"]
        ax.axvline(null_m, c="k", ls=":", label=f"shuffle mean {null_m:+.2f}")
        ax.set_yticks(y)
        ax.set_yticklabels(sub.celltype, fontsize=6)
        ax.set_xlabel("site-stratified R²")
        ax.set_title(f"G1 within-site {method}")
        ax.legend(fontsize=7)
    fig.tight_layout()
    fig.savefig(TRAJ_FIG / "g1_age_per_type.png", dpi=140)
    plt.close(fig)

    fig, ax = plt.subplots(figsize=(6.4, 3.8))
    labels, obs_r, nul_r = [], [], []
    for scheme in ("within_site", "loso"):
        for method in ("pls1", "ridge"):
            labels.append(f"{scheme}\n{method}")
            obs_r.append(age_rec[scheme][method]["r2"] or np.nan)
            nul_r.append(age_rec[scheme][method]["null"]["mean"])
    x = np.arange(len(labels))
    ax.bar(x - 0.18, obs_r, 0.36, color="tab:red", label="observed")
    ax.bar(x + 0.18, nul_r, 0.36, color="0.75", label="shuffle null")
    ax.axhline(0, c="0.4", lw=0.6)
    ax.set_xticks(x)
    ax.set_xticklabels(labels, fontsize=7)
    ax.set_ylabel("primary age R²")
    ax.set_title("G1  age direction vs permutation null")
    ax.legend(fontsize=7)
    fig.tight_layout()
    fig.savefig(TRAJ_FIG / "g1_age_vs_null.png", dpi=140)
    plt.close(fig)

    fig, ax = plt.subplots(figsize=(5.6, 3.4))
    labs, accs, nulls, chances = [], [], [], []
    for scheme in ("within_site", "loso"):
        labs.append(scheme)
        accs.append(age_rec[scheme]["identity"]["acc"])
        nulls.append(age_rec[scheme]["identity"]["null_acc"]["mean"])
        chances.append(age_rec[scheme]["identity"]["chance"])
    x = np.arange(len(labs))
    ax.bar(x - 0.2, accs, 0.4, color="tab:green", label="subspace nearest centroid")
    ax.bar(x + 0.2, nulls, 0.4, color="0.75", label="label-shuffle null")
    ax.axhline(chances[0], c="k", ls=":", label=f"chance {chances[0]:.3f}")
    ax.set_xticks(x)
    ax.set_xticklabels(labs)
    ax.set_ylim(0, 1.05)
    ax.set_ylabel("accuracy")
    ax.set_title("G1  identity from subspace projection")
    ax.legend(fontsize=7)
    fig.tight_layout()
    fig.savefig(TRAJ_FIG / "g1_identity_acc.png", dpi=140)
    plt.close(fig)
    log("[G1] figures written")


if __name__ == "__main__":
    run()
