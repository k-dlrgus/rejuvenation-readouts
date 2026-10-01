"""C1 — planted-signal control (semi-synthetic, primary).

Permute age within bank → age*. Plant β (age*−mean) wᵀ into the real log-CPM matrix.
Calibrate β so within-site site-stratified ridge R² is 0.30 ± 0.03 (also 0.15 / 0.50
for dense_shared). Run the gene-space L2 battery on (X_planted, age*).

Usage: python src/posctrl_c1.py
       python src/posctrl_c1.py dense_shared 0.30
"""
from __future__ import annotations

import sys
from pathlib import Path

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parent))
from posctrl_common import (  # noqa: E402
    PC_DIR, PC_SEED, W_SEED, CAL_SEED, DATASET_ID, REGIMES, N_SUPPORT,
    N_BOOT_FULL, N_BOOT_HALF, N_PERM_FULL, N_PERM_HALF, TARGET_R2_PRIMARY,
    CAL_R2_TOL, N_TYPES_CPM, CPM_CUT, SEX_ANGLE_REF,
    Logger, pc_log_banner, load_phase1_matrix, permute_age_within_site_obs,
    eligible_sparse_genes, draw_planted_w, save_planted_w, plant, donor_mean,
    calibrate_beta, full_within_site_r2, cheap_within_site_r2,
    run_l2_gene_battery, save_cell, load_cell, cell_done, cell_path,
    dump_json, load_json, jsonable, progress_snapshot,
    p1_blocks_c1_c3, pass_transfer_primary, refresh_progress, finite,
)


def cell_name(regime, target_r2):
    rtag = f"r{int(round(float(target_r2)*100)):03d}"
    return f"c1_{regime}_{rtag}"


def _ensure_age_star(obs, log):
    path = PC_DIR / "age_star.npy"
    meta_p = PC_DIR / "age_star.json"
    if path.exists() and meta_p.exists():
        y = np.load(path)
        log(f"[C1] loaded frozen age* from {path.name}  mean={float(np.mean(y)):.2f}")
        return y
    rng = np.random.default_rng(PC_SEED)
    obs_p = permute_age_within_site_obs(obs, rng)
    y = obs_p.age.to_numpy(float)
    np.save(path, y)
    dump_json(meta_p, dict(
        seed=PC_SEED, n=int(len(y)), mean=float(np.mean(y)),
        donor_mean=donor_mean(obs, y),
        site_means=obs_p.groupby("Source").age.mean().to_dict(),
        note="donor-level within-bank permutation of chronological age",
    ))
    log(f"[C1] froze age*  donor-mean={donor_mean(obs, y):.2f}  "
        f"site means={obs_p.groupby('Source').age.mean().to_dict()}")
    return y


def _ensure_w(Y, genes, obs, types, log):
    flag = PC_DIR / "planted_w_ready.json"
    elig, hits = eligible_sparse_genes(Y, obs, types)
    log(f"[C1] sparse eligible genes (mean CPM>{CPM_CUT} in ≥{N_TYPES_CPM}/20 types): {len(elig)} / {Y.shape[1]}")
    if flag.exists() and all((PC_DIR / f"planted_w_{r}.npz").exists() for r in REGIMES):
        log("[C1] planted w already frozen")
        packs = {}
        for r in REGIMES:
            z = np.load(PC_DIR / f"planted_w_{r}.npz", allow_pickle=True)
            kind = str(np.asarray(z["kind"]).reshape(-1)[0])
            if kind == "shared":
                sup = z["support"] if z["support"].size else None
                packs[r] = dict(
                    kind=kind, w=z["w"],
                    support=None if sup is None or len(sup) == 0 else np.asarray(sup, int),
                )
            else:
                W = z["w"]
                if int(z["support_per_type"]):
                    off = z["support_offsets"]
                    flat = z["support"]
                    sups = [np.asarray(flat[int(off[i]):int(off[i + 1])], int) for i in range(len(types))]
                else:
                    sups = None
                packs[r] = dict(kind=kind, w=W, support=sups)
        return packs, elig
    rng = np.random.default_rng(W_SEED)
    packs = draw_planted_w(Y.shape[1], types, elig, rng)
    for r in REGIMES:
        save_planted_w(r, packs[r], genes, types)
        log(f"[C1] wrote planted_w_{r}.npz  kind={packs[r]['kind']}")
    dump_json(flag, dict(
        seed=W_SEED, n_eligible=int(len(elig)), n_support=N_SUPPORT,
        n_types_cpm=N_TYPES_CPM, regimes=list(REGIMES),
    ))
    pd.DataFrame(dict(n_types_mean_cpm_gt1=hits)).to_csv(PC_DIR / "sparse_eligible_hits.csv", index=False)
    return packs, elig


def _ensure_beta(Y, obs, ystar, types, packs, regime, target_r2, log):
    key = f"{regime}_r{int(round(target_r2*100)):03d}"
    path = PC_DIR / "frozen_beta.json"
    rec = load_json(path) if path.exists() else {}
    if key in rec and rec[key].get("beta") is not None:
        prev = rec[key]
        prev_r2 = prev.get("cheap_then_full_r2")
        if finite(prev_r2) and abs(float(prev_r2) - target_r2) <= CAL_R2_TOL:
            log(f"[C1] frozen β[{key}]={prev['beta']}  full-fold R²={prev_r2}")
            return float(prev["beta"]), prev
        log(f"[C1] frozen β[{key}] R²={prev_r2} outside {target_r2:.2f}±{CAL_R2_TOL}; recalibrating")
    log(f"[C1] calibrating β for {key}  target R²={target_r2:.2f} ± {CAL_R2_TOL:.2f}")
    beta, grid = calibrate_beta(Y, obs, ystar, types, packs[regime], target_r2, log)
    grid.to_csv(PC_DIR / f"cal_{key}_grid.csv", index=False)
    history = []

    def _measure(b):
        Yp = plant(Y, ystar, obs, types, packs[regime], b)
        r2 = full_within_site_r2(Yp, obs, ystar)
        history.append(dict(beta=float(b), full_fold_r2=float(r2) if finite(r2) else None))
        log(f"[C1] measured full-fold within-site R²={r2:+.3f} at β={b:.6g}")
        return r2

    def _next_beta(b, r2):
        # Prefer log-β interpolation between measured points that straddle the target.
        ok = [(h["beta"], h["full_fold_r2"]) for h in history
              if h.get("full_fold_r2") is not None and np.isfinite(h["full_fold_r2"])]
        below = [x for x in ok if x[1] <= target_r2]
        above = [x for x in ok if x[1] >= target_r2]
        if below and above:
            b0, r0 = max(below, key=lambda t: t[0])
            b1, r1 = min(above, key=lambda t: t[0])
            if abs(r1 - r0) > 1e-9 and b0 > 0 and b1 > 0:
                t = (target_r2 - r0) / (r1 - r0)
                return float(np.exp(np.log(b0) + t * (np.log(b1) - np.log(b0))))
        if not (finite(r2) and 0.01 < r2 < 0.95):
            return float(b * 1.5) if (not finite(r2) or r2 <= target_r2) else float(b * 0.5)
        s_obs = r2 / max(1.0 - r2, 1e-6)
        s_tgt = target_r2 / max(1.0 - target_r2, 1e-6)
        return float(b * np.sqrt(s_tgt / max(s_obs, 1e-12)))

    achieved = _measure(beta)
    n_extra = 0
    while (not finite(achieved) or abs(achieved - target_r2) > CAL_R2_TOL) and n_extra < 6:
        beta2 = _next_beta(beta, achieved)
        if not np.isfinite(beta2) or beta2 <= 0:
            break
        if any(abs(h["beta"] - beta2) / max(abs(beta2), 1e-12) < 1e-4 for h in history):
            okh = [(h["beta"], h["full_fold_r2"]) for h in history
                   if h.get("full_fold_r2") is not None and np.isfinite(h["full_fold_r2"])]
            blo = [x[0] for x in okh if x[1] <= target_r2]
            abv = [x[0] for x in okh if x[1] >= target_r2]
            if blo and abv:
                beta2 = float(np.sqrt(max(blo) * min(abv)))
            if any(abs(h["beta"] - beta2) / max(abs(beta2), 1e-12) < 1e-4 for h in history):
                log("[C1] next β already measured; stop iterate")
                break
        log(f"[C1] rescale β {beta:.6g} → {beta2:.6g} (R² {achieved:+.3f} → want {target_r2:.2f})")
        achieved = _measure(beta2)
        beta = beta2
        n_extra += 1
    # pick closest measured
    ok = [h for h in history if h.get("full_fold_r2") is not None and np.isfinite(h["full_fold_r2"])]
    if ok:
        best = min(ok, key=lambda h: abs(h["full_fold_r2"] - target_r2))
        beta, achieved = float(best["beta"]), float(best["full_fold_r2"])
    in_band = bool(finite(achieved) and abs(achieved - target_r2) <= CAL_R2_TOL)
    log(f"[C1] freeze β={beta:.6g}  full-fold R²={achieved:+.3f}  in_band={in_band}  "
        f"target {target_r2:.2f}±{CAL_R2_TOL:.2f}")
    info = dict(
        beta=float(beta), cheap_then_full_r2=float(achieved), target_r2=float(target_r2),
        cal_seed=CAL_SEED, in_band=in_band, history=history, n_extra=n_extra,
    )
    rec[key] = info
    dump_json(path, rec)
    return float(beta), info


def _run_cell(regime, target_r2, Y, obs, ystar, types, packs, log):
    name = cell_name(regime, target_r2)
    if cell_done(name):
        log(f"[C1] {name} already complete")
        return load_cell(name)
    reduced = not (abs(target_r2 - TARGET_R2_PRIMARY) < 1e-9)
    n_boot = N_BOOT_HALF if reduced else N_BOOT_FULL
    n_perm = N_PERM_HALF if reduced else N_PERM_FULL
    beta, cal = _ensure_beta(Y, obs, ystar, types, packs, regime, target_r2, log)
    log(f"[C1] plant {regime} target_r2={target_r2:.2f} β={beta:.6g}  "
        f"n_boot={n_boot}{' (B/2)' if reduced else ''} n_perm={n_perm}")
    Yp = plant(Y, ystar, obs, types, packs[regime], beta)
    rng = np.random.default_rng(PC_SEED)
    bat = run_l2_gene_battery(
        Yp, obs, ystar, rng, log,
        n_boot=n_boot, n_perm=n_perm,
        planted_pack=packs[regime], y_is_binary=False,
        split_on_chronological=False, tag=name,
    )
    rec = dict(
        cell=name, regime=regime, target_r2=float(target_r2),
        beta=float(beta), cal=cal, b_reduced=bool(reduced),
        n_boot=n_boot, n_perm=n_perm,
        achieved_within_site_r2=bat.get("within_site_r2"),
        within_site_null=bat.get("within_site_null"),
        **{k: bat.get(k) for k in (
            "boot_median_angle", "recovery_angle", "loso_r2", "loso_null",
            "transfer_H_to_M", "transfer_M_to_H", "transfer_H_to_M_null", "transfer_M_to_H_null",
            "young_old_balanced", "young_old_within_H", "young_old_within_M",
            "pass_angle", "pass_transfer", "pass_transfer_any", "pass_transfer_primary",
            "angle_vs_sex_ref", "pass_both",
            "support_frac", "elapsed_s",
        )},
        path=str(cell_path(name)),
        cal_achieved_before_battery=cal.get("cheap_then_full_r2"),
        sex_angle_ref=SEX_ANGLE_REF,
    )
    rec["pass_transfer"] = pass_transfer_primary(
        rec.get("transfer_H_to_M"), rec.get("transfer_M_to_H"),
        rec.get("transfer_H_to_M_null"), rec.get("transfer_M_to_H_null"),
    )
    rec["pass_transfer_primary"] = rec["pass_transfer"]
    rec["pass_both"] = rec["pass_transfer"]  # angle is not a gate under supersession
    if finite(rec.get("boot_median_angle")):
        rec["angle_vs_sex_ref"] = float(rec["boot_median_angle"]) - SEX_ANGLE_REF
    save_cell(name, jsonable(rec))
    log(f"[C1] wrote {cell_path(name)}")
    refresh_progress()
    from posctrl_findings import write_findings
    write_findings()
    return rec


def jobs_from_args(argv):
    if len(argv) >= 3:
        return [(argv[1], float(argv[2]))]
    if len(argv) == 2:
        return [(argv[1], TARGET_R2_PRIMARY)]
    jobs = [(r, TARGET_R2_PRIMARY) for r in REGIMES]
    jobs += [("dense_shared", 0.15), ("dense_shared", 0.50)]
    return jobs


def run(argv=None):
    argv = list(sys.argv if argv is None else argv)
    log = Logger(PC_DIR / "c1_report.txt", mode="a" if (PC_DIR / "c1_report.txt").exists() else "w")
    try:
        pc_log_banner(log, "C1 — planted-signal control")
        if p1_blocks_c1_c3():
            log("[C1] STOP P1 is set — refusing to run planted controls.")
            return dict(skipped=True, reason="STOP P1")
        log("[C1] P1 superseded 2026-09-16; 35° bar retired; running planted controls.")
        refresh_progress()
        data = load_phase1_matrix(log)
        Y, genes, obs = data["Y"], data["genes"], data["obs"]
        types = sorted(pd.unique(obs.celltype.astype(str)))
        ystar = _ensure_age_star(obs, log)
        packs, elig = _ensure_w(Y, genes, obs, types, log)
        dump_json(PC_DIR / "c1_setup.json", dict(
            n_donors=int(obs.donor.nunique()), n_types=len(types), n_genes=int(Y.shape[1]),
            n_eligible_sparse=int(len(elig)), seed=PC_SEED, w_seed=W_SEED, cal_seed=CAL_SEED,
        ))
        recs = {}
        for regime, t_r2 in jobs_from_args(argv):
            recs[cell_name(regime, t_r2)] = _run_cell(regime, t_r2, Y, obs, ystar, types, packs, log)
        return recs
    finally:
        log.close()


if __name__ == "__main__":
    run()
