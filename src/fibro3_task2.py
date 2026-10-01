"""Task 2 — extrapolation distance vs frozen age score. Report-only, no gate."""
from __future__ import annotations

import sys
from pathlib import Path

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parent))
from fibro3_common import (  # noqa: E402
    FIBRO3_DIR, FIBRO3_SEED, FIBRO3_BOOT, N_PERM, N_BOOT, MIN_CELLS, PCA_MAHAL_K,
    PREREG_TASK2, FROZEN_RULER,
    StopStep, Logger, dump_json, jsonable, fibro3_log_banner,
    load_manifest, save_manifest, record_failure, log_columns,
    load_frozen_ruler, progress_snapshot,
)
from fibro2_extrap import _pca_mahal_fit, _nn_and_mahal, _gtex_z  # noqa: E402
from fibro_stage1 import load_pack  # noqa: E402
from gtex_common import spearman_safe, pearson_safe, pred_scores  # noqa: E402
from trajectory_common import permutation_p  # noqa: E402


def _row_bootstrap_ci(x, y, n_boot=N_BOOT, seed=FIBRO3_BOOT):
    x = np.asarray(x, float)
    y = np.asarray(y, float)
    m = np.isfinite(x) & np.isfinite(y)
    x, y = x[m], y[m]
    n = int(len(x))
    rng = np.random.default_rng(int(seed))
    vals = []
    if n < 8:
        return dict(p025=np.nan, p975=np.nan, n=n, n_ok=0)
    for _ in range(int(n_boot)):
        idx = rng.integers(0, n, size=n)
        r = spearman_safe(x[idx], y[idx])
        if np.isfinite(r):
            vals.append(float(r))
    v = np.asarray(vals, float)
    if v.size < 2:
        return dict(p025=np.nan, p975=np.nan, n=n, n_ok=int(v.size))
    return dict(
        p025=float(np.percentile(v, 2.5)),
        p975=float(np.percentile(v, 97.5)),
        median=float(np.median(v)),
        n=n, n_ok=int(v.size),
    )


def _rho_perm(x, y, n_perm=N_PERM, seed=FIBRO3_SEED):
    x = np.asarray(x, float)
    y = np.asarray(y, float)
    m = np.isfinite(x) & np.isfinite(y)
    x, y = x[m], y[m]
    rho = spearman_safe(x, y)
    rng = np.random.default_rng(int(seed))
    nulls = []
    for _ in range(int(n_perm)):
        yp = y.copy()
        rng.shuffle(yp)
        nulls.append(spearman_safe(x, yp))
    nulls = np.asarray(nulls, float)
    rho_null = float(np.nanmedian(nulls)) if np.isfinite(nulls).any() else np.nan
    p = permutation_p(rho, nulls, greater=True)
    return rho, rho_null, p, nulls


def _corr_block(dist, age, tag, donor=None):
    rho, rho_null, p, nulls = _rho_perm(dist, age)
    ci = _row_bootstrap_ci(dist, age)
    sc = pred_scores(np.asarray(age, float), np.asarray(dist, float))
    rec = dict(
        distance=tag, n=int(np.isfinite(dist).sum()),
        rho=rho, rho_null=rho_null, rho_p=p,
        rho_ci_lo=ci.get("p025"), rho_ci_hi=ci.get("p975"),
        r=sc["r"], r2=sc["r2"], cal_r2=sc["cal_r2"],
        n_perm=N_PERM, n_boot=N_BOOT, seed=FIBRO3_SEED, boot_seed=FIBRO3_BOOT,
        age_rises_with_distance=bool(np.isfinite(rho) and rho > 0),
    )
    if donor is not None:
        rec["cell_line"] = donor
    return rec, nulls


def run_task2(log=None):
    if not (FIBRO3_DIR / "PREREG_TASK2.flag").exists():
        raise StopStep("prereg", "PREREG_TASK2.flag missing")
    close_log = False
    if log is None:
        log = Logger(FIBRO3_DIR / "t2_report.txt")
        close_log = True
    fibro3_log_banner(log, "TASK2")
    log(PREREG_TASK2)
    zpath = FIBRO3_DIR / "t1_cluster_z.npz"
    tpath = FIBRO3_DIR / "t1_cluster_table.csv"
    if not zpath.exists():
        raise StopStep("task2", f"missing {zpath} — run Task 1")
    if not tpath.exists():
        raise StopStep("task2", f"missing {tpath}")
    if not FROZEN_RULER.exists():
        raise StopStep("frozen_ruler", f"missing {FROZEN_RULER}")
    frozen = load_frozen_ruler()
    pack = load_pack()
    log(f"[task2] GTEx pack X shape={np.asarray(pack['X']).shape} obs columns actually read: {list(pack['obs'].columns)}")
    log_columns("t2_gtex_obs", list(pack["obs"].columns), str(FIBRO3_DIR.parent / "fibro" / "stage1_obs.csv"))
    Xz = _gtex_z(pack, frozen)
    log(f"[task2] GTEx train z n={Xz.shape[0]} p={Xz.shape[1]}")
    pca = _pca_mahal_fit(Xz, PCA_MAHAL_K, FIBRO3_SEED)
    log(f"[task2] PCA k={pca['k']}")

    z = np.load(zpath, allow_pickle=True)
    log(f"[task2] t1_cluster_z keys actually read: {list(z.files)}")
    log_columns("t2_t1_z", list(z.files), str(zpath))
    if "Z" not in z.files:
        raise StopStep("task2", f"{zpath.name} missing Z")
    Z = np.asarray(z["Z"], np.float64)
    meta = pd.read_csv(tpath)
    log(f"[task2] t1_cluster_table columns actually read: {list(meta.columns)}")
    log_columns("t2_t1_table", list(meta.columns), str(tpath))
    if len(meta) != len(Z):
        raise StopStep("task2", f"table n={len(meta)} Z n={len(Z)}")
    nn, nn_idx, mahal = _nn_and_mahal(Z, Xz, pca)
    out = meta.copy()
    out["nn_euclidean"] = nn
    out["nn_train_index"] = nn_idx
    out["mahalanobis_pca"] = mahal
    out["pca_k"] = int(pca["k"])
    out.to_csv(FIBRO3_DIR / "t2_extrap.csv", index=False)

    scored = out[~out.below_min_cells].copy() if "below_min_cells" in out.columns else out
    scored = scored[np.isfinite(scored.age_score.to_numpy(float))]
    if len(scored) < 8:
        raise StopStep(
            "task2",
            f"n={len(scored)} donor×cluster×timepoint rows with n_cells≥{MIN_CELLS}; "
            "spearman_safe requires n≥8. Not lowering the bar.",
        )
    age = scored.age_score.to_numpy(float)
    rows = []
    nulls_store = {}
    for tag, dist in (
        ("nn_euclidean", scored.nn_euclidean.to_numpy(float)),
        ("mahalanobis_pca", scored.mahalanobis_pca.to_numpy(float)),
    ):
        rec, nulls = _corr_block(dist, age, tag)
        rec["scope"] = "all_donor_cluster_timepoint_n_ge_20"
        rows.append(rec)
        nulls_store[tag] = nulls
        log(f"[task2] {tag} ρ={rec['rho']:+.3f} null={rec['rho_null']:+.3f} p={rec['rho_p']:+.3f} "
            f"CI=[{rec['rho_ci_lo']:+.3f}, {rec['rho_ci_hi']:+.3f}] "
            f"age_rises_with_distance={rec['age_rises_with_distance']} n={rec['n']}")
        for line, g in scored.groupby("cell_line"):
            if len(g) < 8:
                rows.append(dict(
                    distance=tag, scope=f"per_donor_{line}", cell_line=line,
                    n=int(len(g)), rho=np.nan, note="n<8; spearman_safe not applied",
                ))
                continue
            rec_d, _ = _corr_block(g[tag].to_numpy(float), g.age_score.to_numpy(float), tag, donor=line)
            rec_d["scope"] = f"per_donor_{line}"
            rows.append(rec_d)
            log(f"[task2] {line} {tag} ρ={rec_d['rho']:+.3f} CI=[{rec_d['rho_ci_lo']:+.3f}, {rec_d['rho_ci_hi']:+.3f}] n={rec_d['n']}")

    corr = pd.DataFrame(rows)
    corr.to_csv(FIBRO3_DIR / "t2_correlation.csv", index=False)
    np.savez_compressed(FIBRO3_DIR / "t2_null_rho.npz", **{k: np.asarray(v, float) for k, v in nulls_store.items()})

    pooled = corr[corr.scope == "all_donor_cluster_timepoint_n_ge_20"]
    obs_notes = []
    for _, r in pooled.iterrows():
        if bool(r.age_rises_with_distance):
            obs_notes.append(
                f"{r.distance}: age score rises with distance ρ={float(r.rho):+.3f} "
                f"CI=[{float(r.rho_ci_lo):+.3f}, {float(r.rho_ci_hi):+.3f}] n={int(r.n)}. "
                "Observation only: the ruler is being read outside its calibrated region "
                "and the d10 value is an extrapolation artifact rather than a measurement. "
                "Not a gate and does not overturn any verdict."
            )
        else:
            obs_notes.append(
                f"{r.distance}: age score does not rise with distance ρ={float(r.rho):+.3f} "
                f"CI=[{float(r.rho_ci_lo):+.3f}, {float(r.rho_ci_hi):+.3f}] n={int(r.n)}. "
                "Not a gate."
            )
    summary = dict(
        pca_k=int(pca["k"]), gtex_n=int(Xz.shape[0]), gtex_p=int(Xz.shape[1]),
        n_rows_all=int(len(out)), n_rows_scored=int(len(scored)),
        min_cells=MIN_CELLS, note="report-only; no pass/fail; not a gate",
        observations=obs_notes,
        pooled=pooled.to_dict(orient="records"),
    )
    dump_json(FIBRO3_DIR / "t2_summary.json", jsonable(summary))
    man = load_manifest()
    man["status"] = "TASK2_DONE"
    man["t2_observations"] = obs_notes
    save_manifest(man)
    progress_snapshot("Task 3 per-cluster random-direction null", stop="TASK2_DONE")
    if close_log:
        log.close()
    return summary


if __name__ == "__main__":
    try:
        run_task2()
    except StopStep as e:
        record_failure(e.step, e.message, e.details)
        progress_snapshot(f"fix STOP [{e.step}]", stop=f"STOP {e.step}: {e.message}")
        raise
