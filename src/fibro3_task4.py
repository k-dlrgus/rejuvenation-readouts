"""Task 4 — supersede T-A primary to CELLxGENE 10x. Option (a) is not run."""
from __future__ import annotations

import sys
from pathlib import Path

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parent))
from fibro3_common import (  # noqa: E402
    FIBRO3_DIR, FIBRO3_SEED, FIBRO3_BOOT, N_PERM, N_BOOT,
    PREREG_TASK4, PREREG_TASK4_FLAG, CXG_DATASET_ID, CXG_TAG, CXG_SCORES,
    CXG_RESULT, CXG_PROJ, CXG_H5AD, FROZEN_RULER,
    StopStep, Logger, dump_json, jsonable, fibro3_log_banner,
    load_manifest, save_manifest, record_failure, log_columns,
    load_frozen_ruler, progress_snapshot, pass_rho_bar, rho_null_donor,
)
from fibro_common import unit  # noqa: E402
from gtex_common import pred_scores, spearman_safe  # noqa: E402
from fibro2_common import ADULT_MIN_AGE  # noqa: E402


def _task4_reading(out):
    passed = bool(out.get("pass_rho")) and float(out.get("rho") or 0) > 0
    if passed:
        text = (
            "the frozen fibroblast ruler is a validated instrument outside GTEx and "
            "the Stage 2 trajectory is interpretable as a statement about the cells."
        )
        return dict(key="instrument_validated", text=text, interpretable=True, pass_rho=True)
    text = (
        "Stage 2 remains uninterpretable and every number in Tasks 1–3 is descriptive only."
    )
    return dict(key="instrument_failed", text=text, interpretable=False, pass_rho=False)


def _bootstrap_donor_rho(y, pred, donor, n_boot=N_BOOT, seed=FIBRO3_BOOT):
    y = np.asarray(y, float)
    pred = np.asarray(pred, float)
    donor = np.asarray(donor).astype(str)
    rng = np.random.default_rng(int(seed))
    d_u = pd.unique(donor)
    n_d = int(len(d_u))
    by = {d: np.flatnonzero(donor == d) for d in d_u}
    vals = []
    for _ in range(int(n_boot)):
        take_d = rng.choice(d_u, size=n_d, replace=True)
        yy, pp = [], []
        for d in take_d:
            rows = by[d]
            yy.extend(y[rows].tolist())
            pp.extend(np.asarray(pred)[rows].tolist())
        r = spearman_safe(np.asarray(yy, float), np.asarray(pp, float))
        if np.isfinite(r):
            vals.append(float(r))
    if len(vals) < 2:
        return np.nan, np.nan, int(len(vals))
    return float(np.percentile(vals, 2.5)), float(np.percentile(vals, 97.5)), int(len(vals))


def _score_from_proj(frozen, log):
    if not CXG_PROJ.exists():
        raise StopStep("task4", f"missing frozen projection {CXG_PROJ}")
    z = np.load(CXG_PROJ, allow_pickle=True)
    keys = list(z.files)
    log(f"[task4] proj npz keys actually read: {keys} path={CXG_PROJ}")
    log_columns("t4_proj_npz", keys, str(CXG_PROJ))
    need = {"pred", "y", "donor", "Z"}
    missing = sorted(need - set(keys))
    if missing:
        raise StopStep("task4", f"{CXG_PROJ.name} missing fields {missing}. Not substituting.")
    pred = np.asarray(z["pred"], float)
    y = np.asarray(z["y"], float)
    donor = np.asarray(z["donor"]).astype(str)
    Z = np.asarray(z["Z"], float)
    w = np.asarray(frozen["w"], float)
    w, _ = unit(w)
    pred_chk = Z @ w
    if pred_chk.shape != pred.shape:
        raise StopStep("task4", f"Z@w shape {pred_chk.shape} != pred {pred.shape}")
    max_diff = float(np.nanmax(np.abs(pred_chk - pred)))
    log(f"[task4] freeze-and-project check: max |Z@w − pred|={max_diff:.6g}")
    if max_diff > 1e-4:
        raise StopStep(
            "task4",
            f"frozen ruler · Z does not match stored pred (max abs {max_diff}). Not using a rescaled score.",
        )
    return pred, y, donor, dict(n=int(len(y)), max_zw_diff=max_diff, npz_keys=keys)


def _obs_from_scores(log):
    if not CXG_SCORES.exists():
        raise StopStep("task4", f"missing {CXG_SCORES}")
    obs = pd.read_csv(CXG_SCORES)
    log(f"[task4] scores csv columns actually read: {list(obs.columns)} path={CXG_SCORES} n={len(obs)}")
    log_columns("t4_scores_csv", list(obs.columns), str(CXG_SCORES))
    need = {"donor", "age_years", "age_score"}
    missing = sorted(need - set(obs.columns))
    if missing:
        raise StopStep("task4", f"{CXG_SCORES.name} missing columns {missing}. Not substituting.")
    return obs


def run_task4(log=None):
    if not PREREG_TASK4_FLAG.exists():
        raise StopStep("prereg", "PREREG_TASK4.flag missing — write supersession before any ρ")
    close_log = False
    if log is None:
        log = Logger(FIBRO3_DIR / "t4_report.txt")
        close_log = True
    fibro3_log_banner(log, "TASK4")
    log(PREREG_TASK4)
    log("Option (a) is not run. GSE113957 is not scored. GSE226189 is not the primary.")
    if not FROZEN_RULER.exists():
        raise StopStep("frozen_ruler", f"missing {FROZEN_RULER}")
    frozen = load_frozen_ruler()
    log(f"[frozen] n_genes={len(frozen['w'])} n_donors={int(frozen['n_donors'])} "
        f"method={frozen['method']} regime={frozen['regime']} path={FROZEN_RULER}")

    if not CXG_SCORES.exists() and not CXG_PROJ.exists():
        details = dict(scores=str(CXG_SCORES), proj=str(CXG_PROJ), h5ad=str(CXG_H5AD),
                       h5ad_exists=CXG_H5AD.exists())
        raise StopStep(
            "task4",
            "superseded primary CELLxGENE 10x cohort has no frozen scores and no projection npz. "
            "Not opening a substitute cohort. Not running option (a).",
            details,
        )

    obs = _obs_from_scores(log)
    pred_npz, y_npz, donor_npz, chk = _score_from_proj(frozen, log)

    y = obs.age_years.to_numpy(float)
    pred_csv = obs.age_score.to_numpy(float)
    donor = obs.donor.astype(str).to_numpy()
    if len(y) != len(y_npz):
        raise StopStep("task4", f"scores csv n={len(y)} != proj n={len(y_npz)}")
    if float(np.nanmax(np.abs(y - y_npz))) > 1e-9:
        raise StopStep("task4", "age_years in csv does not match proj y. Not repairing.")
    if float(np.nanmax(np.abs(pred_csv - pred_npz))) > 1e-4:
        raise StopStep("task4", "age_score in csv does not match proj pred. Not repairing.")
    if list(donor) != list(donor_npz):
        raise StopStep("task4", "donor order in csv does not match proj donor. Not repairing.")

    n_child = int((y < ADULT_MIN_AGE).sum())
    if n_child:
        log(f"[task4] excluding {n_child} samples with stated age < {ADULT_MIN_AGE}")
        keep = y >= ADULT_MIN_AGE
        y, pred_csv, donor = y[keep], pred_csv[keep], donor[keep]
        obs = obs.loc[keep].reset_index(drop=True)
    if len(y) < 8:
        raise StopStep("task4", f"n={len(y)} after adult filter; spearman_safe requires n≥8. Not lowering the bar.")

    # one row per donor already in T-A scores; still run donor-level null
    sc, rho_null, p, nulls = rho_null_donor(y, pred_csv, donor, n_perm=N_PERM, seed=FIBRO3_SEED)
    ci_lo, ci_hi, n_boot_ok = _bootstrap_donor_rho(y, pred_csv, donor)
    passed = pass_rho_bar(sc["rho"], rho_null)
    age_src = None
    if "age_source" in obs.columns:
        age_src = sorted(set(obs.age_source.astype(str)))
    assay = None
    if "assay" in obs.columns:
        assay = obs.assay.astype(str).value_counts().to_dict()
    prev = {}
    if CXG_RESULT.exists():
        from brain_phase1_common import load_json
        prev = load_json(CXG_RESULT)
        log(f"[task4] T-A cached result.json keys actually read: {sorted(prev.keys())}")
        log_columns("t4_ta_result_json", sorted(prev.keys()), str(CXG_RESULT))

    out = dict(
        accession=CXG_DATASET_ID, role="superseded_primary_external",
        platform_kind="single_cell", platform_detail="10x (restricted)",
        tag=CXG_TAG, choice="b", option_a_run=False,
        n_samples=int(len(obs)), n_donors=int(pd.Series(donor).nunique()),
        n_excluded_child=n_child,
        n_overlap=prev.get("n_overlap"), n_missing_z0=prev.get("n_missing_z0"),
        n_ruler=prev.get("n_ruler"),
        rho=sc["rho"], rho_null=rho_null, rho_p=p,
        rho_ci_lo=ci_lo, rho_ci_hi=ci_hi,
        r=sc["r"], r2=sc["r2"], cal_r2=sc["cal_r2"],
        pass_rho=bool(passed),
        n_perm=N_PERM, n_boot=N_BOOT, n_boot_ok=n_boot_ok,
        seed=FIBRO3_SEED, boot_seed=FIBRO3_BOOT,
        age_min=float(np.nanmin(y)), age_max=float(np.nanmax(y)),
        n_unique_age=int(pd.Series(y).nunique()),
        age_source=age_src, assay_counts=assay,
        frozen_ruler=str(FROZEN_RULER),
        scores_path=str(CXG_SCORES), proj_path=str(CXG_PROJ),
        zw_check=chk,
        used_10x_only=True,
        adult_min_age=ADULT_MIN_AGE,
        note="freeze-and-project on stored T-A Z; ruler not refit; option (a) not run",
    )
    reading = _task4_reading(out)
    out["reading_key"] = reading["key"]
    out["interpretable"] = reading["interpretable"]
    log(f"[task4] ρ={out['rho']:+.3f} null={rho_null:+.3f} p={p:+.3f} "
        f"CI=[{ci_lo:+.3f}, {ci_hi:+.3f}] calR²={sc['cal_r2']:+.3f} R²={sc['r2']:+.3f} "
        f"pass={passed} n={out['n_samples']} n_donors={out['n_donors']}")
    log(f"[task4 reading] key={reading['key']} interpretable={reading['interpretable']} {reading['text']}")

    obs_out = obs.copy()
    obs_out["age_score"] = pred_csv if n_child == 0 else obs.age_score.to_numpy(float)
    obs_out.to_csv(FIBRO3_DIR / "t4_scores.csv", index=False)
    pd.DataFrame([out]).to_csv(FIBRO3_DIR / "t4_result.csv", index=False)
    np.save(FIBRO3_DIR / "t4_null_rho.npy", np.asarray(nulls, float))
    dump_json(FIBRO3_DIR / "t4_result.json", jsonable(out))
    dump_json(FIBRO3_DIR / "t4_reading.json", jsonable(reading))
    dump_json(FIBRO3_DIR / "t4_summary.json", jsonable(dict(result=out, reading=reading)))
    man = load_manifest()
    man["status"] = "TASK4_DONE"
    man["t4_reading"] = reading
    man["t4_pass"] = bool(passed)
    save_manifest(man)
    progress_snapshot("Task 1 joint clustering (settings already frozen)", stop="TASK4_DONE")
    if close_log:
        log.close()
    return dict(result=out, reading=reading)


if __name__ == "__main__":
    try:
        run_task4()
    except StopStep as e:
        record_failure(e.step, e.message, e.details)
        progress_snapshot(f"fix STOP [{e.step}]", stop=f"STOP {e.step}: {e.message}")
        raise
