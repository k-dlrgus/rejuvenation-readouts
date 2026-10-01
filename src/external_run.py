"""Run the SEA-AD MTG external-transfer pipeline.

Order is fixed: Stage 0 → pre-register (before any fit) → freeze DLPFC directions
→ project → FINDINGS_EXTERNAL.md. Does not modify any other file.

Usage: python src/external_run.py
"""
from __future__ import annotations

import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
from external_common import (  # noqa: E402
    EXT_DIR, EXT_SEED, BOOT_SEED, N_PERM, N_BOOT, TRANSFER_NULL_BAR, MIN_N_R,
    PREREG_DATE, PREREG_FLAG, PREREG_BARS, StopStep, load_json, refresh_progress,
)
from external_stage0 import stage0  # noqa: E402
from external_freeze import freeze  # noqa: E402
from external_project import project  # noqa: E402
from external_findings import write_findings  # noqa: E402


def write_prereg():
    """Write the pre-registration flag from Stage 0 fields. Never after a fit."""
    if PREREG_FLAG.exists():
        return str(PREREG_FLAG)
    if (EXT_DIR / "frozen_directions.npz").exists() or (EXT_DIR / "primary_raw.json").exists():
        raise StopStep(
            "prereg",
            "refusing to write the pre-registration flag after a fit already exists",
            dict(frozen=(EXT_DIR / "frozen_directions.npz").exists()),
        )
    summ = load_json(EXT_DIR / "stage0_summary.json")
    text = PREREG_BARS.format(
        date=PREREG_DATE,
        n_perm=N_PERM,
        perm_seed=EXT_SEED,
        n_boot=N_BOOT,
        boot_seed=BOOT_SEED,
        min_n=MIN_N_R,
        null_bar=TRANSFER_NULL_BAR,
        path_field=summ.get("pathology_field"),
        path_file=summ.get("pathology_source_file") or summ.get("donor_metadata_path"),
        low_levels=summ.get("pathology_low_levels"),
        encoding=f"{summ.get('pathology_family')} declared ordinal table",
    )
    PREREG_FLAG.write_text(text.strip() + "\n", encoding="utf-8")
    refresh_progress(next_action="freeze DLPFC raw + identity-residual ridge on 233 donors")
    return str(PREREG_FLAG)


def main():
    try:
        need_stage0 = True
        if (EXT_DIR / "stage0_summary.json").exists() and (EXT_DIR / "manifest.json").exists():
            if load_json(EXT_DIR / "manifest.json").get("status") == "STAGE0_OK":
                need_stage0 = False
        if need_stage0:
            stage0()
        man = load_json(EXT_DIR / "manifest.json")
        if man.get("status") == "STOP":
            write_findings()
            refresh_progress(next_action="STOP recorded in FINDINGS_EXTERNAL.md", stop="STOP")
            return
        write_prereg()
        if not (EXT_DIR / "frozen_directions.npz").exists():
            freeze()
        if not (EXT_DIR / "project_summary.json").exists():
            project()
        write_findings()
        refresh_progress(next_action="done", stop="none")
    except StopStep as e:
        print(f"[STOP] {e.step}: {e.message}", flush=True)
        try:
            write_findings()
        except Exception as werr:
            print(f"[warn] FINDINGS_EXTERNAL.md incomplete: {werr}", flush=True)
        refresh_progress(next_action=f"STOP {e.step}", stop=f"STOP {e.step}")
        raise SystemExit(2)


if __name__ == "__main__":
    main()
