"""C2 — real-data sex control.

C2a: all genes. Must pass (XIST, Y-linked). Gate for C1/C3.
C2b: chrX, chrY, PAR removed.

Target = sex coded ±1 (male +1, female −1). Ridge, plus held-out AUC.
Nulls: donor-level sex permutation within bank (same function as age perm).

Usage: python src/posctrl_c2.py            # both
       python src/posctrl_c2.py c2a
       python src/posctrl_c2.py c2b
"""
from __future__ import annotations

import sys
from pathlib import Path

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parent))
from posctrl_common import (  # noqa: E402
    PC_DIR, PC_SEED, DATASET_ID, N_BOOT_FULL, N_PERM_FULL, PREREG_BLOCK,
    Logger, pc_log_banner, load_phase1_matrix, sex_pm1, sex_chrom_mask,
    run_l2_gene_battery, save_cell, load_cell, cell_done, cell_path,
    pass_angle, pass_transfer_any, progress_snapshot, jsonable,
)


def _run_one(which, Y, genes, obs, types, log):
    tag = which
    if cell_done(tag):
        rec = load_cell(tag)
        log(f"[{which}] already complete → {cell_path(tag)}")
        return rec
    y, n_m, n_f = sex_pm1(obs)
    n_nan = int(np.isnan(y).sum())
    log(f"[{which}] sex ±1  male_rows={n_m} female_rows={n_f} nan={n_nan}  "
        f"donors_male={obs.loc[obs.sex.astype(str).str.lower()=='male','donor'].nunique()}  "
        f"donors_female={obs.loc[obs.sex.astype(str).str.lower()=='female','donor'].nunique()}")
    Y_use, genes_use = Y, genes
    chrom_info = None
    if which == "c2b":
        keep, chrom_info = sex_chrom_mask(genes)
        log(f"[c2b] drop chrX/chrY/PAR: n_drop={chrom_info['n_drop']} n_keep={chrom_info['n_keep']}  "
            f"chrom_xy={chrom_info['n_chrom_xy']} par_symbol={chrom_info['n_par_symbol']}")
        log(f"[c2b] chrom counts: {chrom_info['chrom_counts']}")
        Y_use = np.ascontiguousarray(Y[:, keep])
        genes_use = genes.iloc[np.flatnonzero(keep)].reset_index(drop=True)
        pd.Series(keep.astype(int), index=genes.gene_id.astype(str)).to_csv(
            PC_DIR / "c2b_gene_keep.csv", header=["keep"])
    rng = np.random.default_rng(PC_SEED)
    bat = run_l2_gene_battery(
        Y_use, obs, y, rng, log,
        n_boot=N_BOOT_FULL, n_perm=N_PERM_FULL,
        planted_pack=None, y_is_binary=True,
        split_on_chronological=True, tag=tag,
    )
    rec = dict(
        cell=which, regime=which, target="sex_pm1",
        n_male_rows=n_m, n_female_rows=n_f,
        n_genes=int(Y_use.shape[1]),
        n_boot=N_BOOT_FULL, n_perm=N_PERM_FULL, b_reduced=False,
        chrom=chrom_info,
        **{k: bat.get(k) for k in (
            "boot_median_angle", "recovery_angle", "within_site_r2", "within_site_null",
            "within_site_p", "loso_r2", "loso_null", "loso_p",
            "transfer_H_to_M", "transfer_M_to_H", "transfer_H_to_M_null", "transfer_M_to_H_null",
            "transfer_H_to_M_p", "transfer_M_to_H_p",
            "young_old_balanced", "young_old_within_H", "young_old_within_M",
            "pass_angle", "pass_transfer", "pass_transfer_any", "pass_both",
            "site_strat_auc", "pooled_auc", "elapsed_s",
        )},
        achieved_within_site_r2=bat.get("within_site_r2"),
        path=str(cell_path(tag)),
    )
    rec["p1_fail"] = bool(which == "c2a" and (
        not pass_angle(rec.get("boot_median_angle")) or
        not pass_transfer_any(rec.get("transfer_H_to_M"), rec.get("transfer_M_to_H"))
    ))
    save_cell(tag, jsonable(rec))
    log(f"[{which}] wrote {cell_path(tag)}")
    log(f"[{which}] boot={rec['boot_median_angle']:.1f}°  "
        f"H→M={rec['transfer_H_to_M']:+.3f}  M→H={rec['transfer_M_to_H']:+.3f}  "
        f"ws R²={rec['achieved_within_site_r2']:+.3f}  AUC={rec.get('site_strat_auc')}")
    return rec


def run(which=None):
    log = Logger(PC_DIR / "c2_report.txt", mode="a" if (PC_DIR / "c2_report.txt").exists() else "w")
    try:
        pc_log_banner(log, "C2 — sex control")
        data = load_phase1_matrix(log)
        Y, genes, obs = data["Y"], data["genes"], data["obs"]
        types = sorted(pd.unique(obs.celltype.astype(str)))
        jobs = [which] if which in ("c2a", "c2b") else ["c2a", "c2b"]
        recs = {}
        for job in jobs:
            recs[job] = _run_one(job, Y, genes, obs, types, log)
        if "c2a" in recs and recs["c2a"].get("p1_fail"):
            log("[STOP P1] C2a failed bootstrap ≤35° or transfer >0. Do not proceed to C1/C3.")
            (PC_DIR / "STOP_P1.txt").write_text(
                "stop_P1=True\n"
                f"boot_median_angle={recs['c2a'].get('boot_median_angle')}\n"
                f"transfer_H_to_M={recs['c2a'].get('transfer_H_to_M')}\n"
                f"transfer_M_to_H={recs['c2a'].get('transfer_M_to_H')}\n"
                "reason=C2a sex control failed L2 bars; pipeline suspect\n",
                encoding="utf-8",
            )
        return recs
    finally:
        log.close()


if __name__ == "__main__":
    arg = sys.argv[1] if len(sys.argv) > 1 else None
    run(arg)
