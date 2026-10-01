"""STAGE V2 — the target vector.

Per cell type:
  TARGET = age direction with its identity-subspace component removed
           (project onto the orthogonal complement of the identity subspace,
           then renormalize).

This is the "move along age, not along identity" direction.

Report per cell type:
  - how much of the age direction survived the projection (norm retained)
  - the angle between TARGET and the raw age direction
  - whether TARGET still predicts age (project held-out samples onto TARGET,
    regress on age, both CV schemes, with permutation null)

STOP V2: TARGET loses age predictivity (held-out R² at or below null) in a
majority of cell types. Halt — it would mean the age signal lives entirely
in the identity subspace after all, contradicting G2c.

Usage: python src/target_v2.py
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
from target_common import (  # noqa: E402
    TARGET_DIR, TARGET_FIG, TARGET_SEED, DATASET_ID, N_PERM, STOP_NULL,
    Logger, dump_json, tgt_log_banner, load_phase1_matrix, make_folds,
    fit_full, oof_project, oof_project_both, shuffle_null_oof, per_type_primary_r2,
    type_loses_age, permute_age_within_site_obs, primary_r2, primary_mae,
    permutation_p, summarize_null_col, save_matrix_bundle,
)


def _pack_scores(obs, pred, scheme, method, use_target, nt, nt_pt, log):
    rmap, sc = per_type_primary_r2(obs, pred, scheme)
    prim = primary_r2(scheme, sc)
    mae = primary_mae(scheme, sc)
    null_prim = nt.primary.to_numpy(float)
    pval = permutation_p(prim, null_prim, greater=True)
    null_sum = summarize_null_col(null_prim)
    lose = {}
    pt_null_mean = nt_pt.groupby("celltype").r2.mean() if len(nt_pt) else pd.Series(dtype=float)
    for t, r in rmap.items():
        lose[t] = type_loses_age(r, nt_pt.loc[nt_pt.celltype == t, "r2"] if len(nt_pt) else [])
    n_lose, n_t = int(sum(lose.values())), int(len(lose))
    tag = "TARGET" if use_target else "age"
    log(f"      {method} {tag}: primary R²={prim:+.3f}  MAE={mae:.2f} y  "
        f"null={null_sum['mean']:+.3f}  p={pval:.3f}  types≤null={n_lose}/{n_t}")
    stem = f"v2_{scheme}_{method}_{'tgt' if use_target else 'age'}"
    sc["per_type"].assign(method=method, scheme=scheme, use_target=use_target).to_csv(
        TARGET_DIR / f"{stem}_per_type.csv", index=False)
    nt.assign(method=method, scheme=scheme, use_target=use_target).to_csv(
        TARGET_DIR / f"{stem}_null.csv", index=False)
    nt_pt.assign(method=method, scheme=scheme, use_target=use_target).to_csv(
        TARGET_DIR / f"{stem}_null_per_type.csv", index=False)
    return dict(
        scheme=scheme, method=method, use_target=use_target,
        primary_r2=float(prim), mae=float(mae),
        pooled_r2=float(sc["pooled_median_type_r2"]),
        site_strat_r2=float(sc["site_strat_median_type_r2"]),
        null_mean=null_sum["mean"], null_p95=null_sum["p95"],
        p=float(pval), n_perm=N_PERM,
        null_usable=bool(np.isfinite(null_sum["mean"]) and null_sum["mean"] <= STOP_NULL),
        n_types_lose=n_lose, n_types=n_t,
        per_type_r2=rmap, per_type_lose=lose,
        per_type_null_mean=pt_null_mean.to_dict(),
    )


def _shuffle_both(Y, obs, folds, method, rng, n_perm, svd_cache, scheme, log):
    rows_age, rows_tgt, pt_age, pt_tgt = [], [], [], []
    log(f"      permutation null ({n_perm} donor-level within-site shuffles)...")
    for i in range(n_perm):
        obs_p = permute_age_within_site_obs(obs, rng)
        yp = obs_p.age.to_numpy(float)
        out = oof_project_both(Y, obs, folds, method, y_override=yp, svd_cache=svd_cache)
        for pred, rows, pt in (
            (out["pred_age"], rows_age, pt_age),
            (out["pred_tgt"], rows_tgt, pt_tgt),
        ):
            rmap, sc = per_type_primary_r2(obs_p, pred, scheme, y=yp)
            rows.append(dict(
                perm=i, pooled=sc["pooled_median_type_r2"],
                site_strat=sc["site_strat_median_type_r2"],
                primary=primary_r2(scheme, sc),
            ))
            for t, r in rmap.items():
                pt.append(dict(perm=i, celltype=t, r2=r))
        if (i + 1) % 5 == 0:
            log(f"         perm {i+1}/{n_perm}")
    return (pd.DataFrame(rows_age), pd.DataFrame(pt_age),
            pd.DataFrame(rows_tgt), pd.DataFrame(pt_tgt))


def _scheme_ridge(Y, obs, folds, scheme, rng, log):
    log(f"   OOF ridge (age + TARGET in one pass)  scheme={scheme}")
    out = oof_project_both(Y, obs, folds, "ridge")
    cache = out.get("svd_cache")
    nt_age, pt_age, nt_tgt, pt_tgt = _shuffle_both(
        Y, obs, folds, "ridge", rng, N_PERM, cache, scheme, log)
    rec_age = _pack_scores(obs, out["pred_age"], scheme, "ridge", False, nt_age, pt_age, log)
    rec_tgt = _pack_scores(obs, out["pred_tgt"], scheme, "ridge", True, nt_tgt, pt_tgt, log)
    return rec_age, rec_tgt


def _scheme_block(Y, obs, folds, scheme, method, rng, log, use_target=True):
    log(f"   OOF {method}  target={use_target}  scheme={scheme}")
    out = oof_project(Y, obs, folds, method, use_target=use_target, collect_W=False)
    cache = out.get("svd_cache")
    rmap, sc = per_type_primary_r2(obs, out["pred"], scheme)
    prim = primary_r2(scheme, sc)
    mae = primary_mae(scheme, sc)
    log(f"      primary R²={prim:+.3f}  MAE={mae:.2f} y  "
        f"pooled={sc['pooled_median_type_r2']:+.3f}  strat={sc['site_strat_median_type_r2']:+.3f}")
    log(f"      permutation null ({N_PERM} donor-level within-site shuffles)...")
    nt, nt_pt = shuffle_null_oof(
        Y, obs, folds, method, use_target, rng, N_PERM, cache, scheme)
    null_prim = nt.primary.to_numpy(float)
    pval = permutation_p(prim, null_prim, greater=True)
    null_sum = summarize_null_col(null_prim)
    log(f"      null mean={null_sum['mean']:+.3f}  p95={null_sum['p95']:+.3f}  p={pval:.3f}")
    # per-type lose?
    lose = {}
    pt_null_mean = nt_pt.groupby("celltype").r2.mean() if len(nt_pt) else pd.Series(dtype=float)
    for t, r in rmap.items():
        lose[t] = type_loses_age(r, nt_pt.loc[nt_pt.celltype == t, "r2"] if len(nt_pt) else [])
    n_lose = int(sum(lose.values()))
    n_t = int(len(lose))
    log(f"      types at-or-below null: {n_lose}/{n_t}")
    sc["per_type"].assign(method=method, scheme=scheme, use_target=use_target).to_csv(
        TARGET_DIR / f"v2_{scheme}_{method}_{'tgt' if use_target else 'age'}_per_type.csv",
        index=False)
    nt.assign(method=method, scheme=scheme, use_target=use_target).to_csv(
        TARGET_DIR / f"v2_{scheme}_{method}_{'tgt' if use_target else 'age'}_null.csv",
        index=False)
    nt_pt.assign(method=method, scheme=scheme, use_target=use_target).to_csv(
        TARGET_DIR / f"v2_{scheme}_{method}_{'tgt' if use_target else 'age'}_null_per_type.csv",
        index=False)
    return dict(
        scheme=scheme, method=method, use_target=use_target,
        primary_r2=float(prim), mae=float(mae),
        pooled_r2=float(sc["pooled_median_type_r2"]),
        site_strat_r2=float(sc["site_strat_median_type_r2"]),
        null_mean=null_sum["mean"], null_p95=null_sum["p95"],
        p=float(pval), n_perm=N_PERM,
        null_usable=bool(np.isfinite(null_sum["mean"]) and null_sum["mean"] <= STOP_NULL),
        n_types_lose=n_lose, n_types=n_t,
        per_type_r2=rmap, per_type_lose=lose,
        per_type_null_mean=pt_null_mean.to_dict(),
    )


def run():
    rng = np.random.default_rng(TARGET_SEED)
    log = Logger(TARGET_DIR / "v2_report.txt")
    try:
        return _run(log, rng)
    finally:
        log.close()


def _run(log, rng):
    tgt_log_banner(log, "V2 — target = age dir minus identity subspace")
    v1 = TARGET_DIR / "v1_STOP.txt"
    if v1.exists() and "stop_V1=True" in v1.read_text(encoding="utf-8"):
        log("[V2] STOP V1 is set — refusing to derive a target from an unstable age direction.")
        dump_json(TARGET_DIR / "v2_summary.json", dict(skipped=True, reason="STOP V1"))
        (TARGET_DIR / "v2_STOP.txt").write_text("stop_V2=True\nreason=STOP V1\n", encoding="utf-8")
        return dict(skipped=True, stop_V2=True)

    data = load_phase1_matrix(log)
    Y, genes, obs = data["Y"], data["genes"], data["obs"]
    log(f"[filter] V2 uses ALL {Y.shape[1]:,} genes.")

    pack_r = fit_full(Y, obs, method="ridge")
    pack_p = fit_full(Y, obs, method="pls1")
    types = pack_r["types"]

    # descriptive geometry of TARGET vs age
    geom_r = pack_r["meta"].copy()
    geom_p = pack_p["meta"].copy()
    geom_r["method"] = "ridge"
    geom_p["method"] = "pls1"
    geom = pd.concat([geom_r, geom_p], ignore_index=True)
    geom.to_csv(TARGET_DIR / "v2_geom_full.csv", index=False)
    log("[V2] full-cohort TARGET vs raw age direction (ridge):")
    log(f"   median norm retained = {geom_r.norm_retained.median():.3f}")
    log(f"   median angle TARGET vs age = {geom_r.angle_to_age_deg.median():.2f}°")
    log(f"   median angle age vs identity subspace = {geom_r.angle_to_id_deg.median():.2f}°")
    for _, row in geom_r.iterrows():
        log(f"   {row.celltype}: retained={row.norm_retained:.3f}  "
            f"∠(T,age)={row.angle_to_age_deg:.2f}°  ∠(age,id)={row.angle_to_id_deg:.1f}°")

    save_matrix_bundle(
        TARGET_DIR / "target_W.npz",
        genes, types,
        W_target_ridge=pack_r["W_tgt"],
        W_age_ridge=pack_r["W_age"],
        W_target_pls=pack_p["W_tgt"],
        W_age_pls=pack_p["W_age"],
        basis=pack_r["basis"],
        mu=pack_r["mu"], sd=pack_r["sd"],
        S_id=pack_r["S_id"],
    )
    log("[V2] wrote target_W.npz  (per-type TARGET matrix for later use)")

    rng_folds = np.random.default_rng(TARGET_SEED)
    schemes = make_folds(obs, rng_folds)
    log(f"[folds] seed={TARGET_SEED}  loso={len(schemes['loso'])}  "
        f"within_site={len(schemes['within_site'])}")

    rec = {}
    for name, folds in schemes.items():
        log(f"\n[V2] scheme={name}")
        rec[name] = {}
        rec_age, rec_tgt = _scheme_ridge(Y, obs, folds, name, rng, log)
        rec[name]["ridge_age"] = rec_age
        rec[name]["ridge_target"] = rec_tgt
        rec[name]["pls1_target"] = _scheme_block(
            Y, obs, folds, name, "pls1", rng, log, use_target=True)

    # STOP: majority of cell types lose TARGET age predictivity on within_site (the usable scheme)
    ws = rec["within_site"]["ridge_target"]
    n_lose, n_t = ws["n_types_lose"], ws["n_types"]
    stop = bool(n_lose > n_t / 2)
    (TARGET_DIR / "v2_STOP.txt").write_text(
        f"stop_V2={stop}\n"
        f"within_site_ridge_target_primary_r2={ws['primary_r2']}\n"
        f"null_mean={ws['null_mean']}\n"
        f"types_lose={n_lose}/{n_t}\n"
        f"reason={'TARGET loses age predictivity vs null in a majority of cell types — age lives in identity (contradicts G2c)' if stop else 'TARGET retains held-out age predictivity in a majority of cell types'}\n",
        encoding="utf-8",
    )
    log(f"\n[STOP V2] {stop}  within-site ridge TARGET types at/below null: {n_lose}/{n_t}")
    log(f"          primary R²={ws['primary_r2']:+.3f} vs null {ws['null_mean']:+.3f}")

    # figures
    fig, ax = plt.subplots(figsize=(8.2, 5.6))
    g = geom_r.sort_values("norm_retained")
    y_pos = np.arange(len(g))
    ax.barh(y_pos, g.norm_retained, color="tab:purple", height=0.7)
    ax.set_yticks(y_pos)
    ax.set_yticklabels(g.celltype, fontsize=6)
    ax.set_xlabel("‖(I − P_id) w_age‖  (norm retained)")
    ax.set_xlim(0, 1.02)
    ax.set_title("V2  how much of the ridge age direction survives")
    fig.tight_layout()
    fig.savefig(TARGET_FIG / "v2_norm_retained.png", dpi=140)
    plt.close(fig)

    fig, ax = plt.subplots(figsize=(8.2, 5.6))
    # per-type R² target vs null, within_site ridge
    pt = pd.read_csv(TARGET_DIR / "v2_within_site_ridge_tgt_per_type.csv")
    sub = pt[pt.site == "mean_of_sites"].copy()
    nm = pd.Series(ws["per_type_null_mean"])
    sub["null_mean"] = sub.celltype.map(nm)
    sub = sub.sort_values("r2")
    y_pos = np.arange(len(sub))
    ax.barh(y_pos, sub.r2, color="tab:purple", height=0.65, label="TARGET")
    ax.scatter(sub.null_mean, y_pos, c="k", s=12, zorder=3, label="perm null mean")
    ax.axvline(0, c="0.7", lw=0.6)
    ax.set_yticks(y_pos)
    ax.set_yticklabels(sub.celltype, fontsize=6)
    ax.set_xlabel("site-stratified R² (within-site CV)")
    ax.set_title("V2  TARGET still predicts age?")
    ax.legend(fontsize=7)
    fig.tight_layout()
    fig.savefig(TARGET_FIG / "v2_per_type_r2.png", dpi=140)
    plt.close(fig)

    fig, ax = plt.subplots(figsize=(6.4, 3.8))
    labels, obs_v, null_v = [], [], []
    for sch in ("within_site", "loso"):
        for key, lab in (("ridge_target", "ridge TARGET"), ("ridge_age", "ridge age"),
                         ("pls1_target", "PLS-1 TARGET")):
            recs = rec[sch][key]
            labels.append(f"{sch}\n{lab}")
            obs_v.append(recs["primary_r2"])
            null_v.append(recs["null_mean"])
    x = np.arange(len(labels))
    ax.bar(x, obs_v, color="tab:purple", label="observed")
    ax.scatter(x, null_v, c="k", zorder=3, label="null mean")
    ax.axhline(0, c="0.7", lw=0.6)
    ax.set_xticks(x)
    ax.set_xticklabels(labels, fontsize=7)
    ax.set_ylabel("primary R²")
    ax.set_title("V2  TARGET vs raw age (both CV schemes)")
    ax.legend(fontsize=7)
    fig.tight_layout()
    fig.savefig(TARGET_FIG / "v2_r2_vs_age.png", dpi=140)
    plt.close(fig)

    summary = dict(
        seed=TARGET_SEED, dataset_id=DATASET_ID,
        n_genes=int(Y.shape[1]), n_types=len(types),
        n_identity_axes=int(pack_r["basis"].shape[1]),
        ridge_median_norm_retained=float(geom_r.norm_retained.median()),
        ridge_median_angle_to_age_deg=float(geom_r.angle_to_age_deg.median()),
        ridge_median_angle_to_id_deg=float(geom_r.angle_to_id_deg.median()),
        pls_median_norm_retained=float(geom_p.norm_retained.median()),
        pls_median_angle_to_age_deg=float(geom_p.angle_to_age_deg.median()),
        cv=rec, stop_V2=stop, n_perm=N_PERM,
    )
    dump_json(TARGET_DIR / "v2_summary.json", summary)
    log("[V2] wrote results/target/v2_*")
    log("[V2] done.")
    return summary


if __name__ == "__main__":
    run()
