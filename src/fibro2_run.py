"""Run Prompt F: T-A, T-B, T-C, extrapolation. Write FINDINGS_FIBRO2.md.

Usage: python src/fibro2_run.py [--cell prereg|ta|tb|tc|extrap|findings|all]

Writes T-A/T-B/T-C pre-registration flags and FINDINGS_FIBRO2.md skeleton
BEFORE any scoring. Does not modify existing FINDINGS*.md or FALSIFICATION.md.
Does not open GSE325735.
"""
from __future__ import annotations

import argparse
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
from fibro2_common import (  # noqa: E402
    FIBRO2_DIR, FROZEN_RULER, PREREG_TA, PREREG_TB, PREREG_TC,
    write_all_prereg_flags, progress_snapshot, load_manifest, save_manifest,
    fibro2_log_banner, Logger, StopStep, record_failure,
)
from fibro2_findings import write_findings  # noqa: E402


def _prereg_and_skeleton():
    write_all_prereg_flags()
    man = load_manifest()
    man["status"] = "PREREG_WRITTEN"
    man["frozen_ruler"] = str(FROZEN_RULER)
    save_manifest(man)
    progress_snapshot(
        "T-A resolve GSE113957 (pre-registration is on disk)",
        stop="PREREG_WRITTEN",
        extra="T-A/T-B/T-C pre-registered blocks written before any task run.",
    )
    write_findings()
    return "prereg written; FINDINGS_FIBRO2.md skeleton on disk"


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--cell", default="all",
                    help="prereg | ta | tb | tc | extrap | findings | all")
    args = ap.parse_args()
    log = Logger(FIBRO2_DIR / "run_report.txt", mode="a" if args.cell != "all" else "w")
    fibro2_log_banner(log, f"RUN cell={args.cell}")
    try:
        if args.cell in ("all", "prereg"):
            log(_prereg_and_skeleton())
        if args.cell in ("all", "ta"):
            write_all_prereg_flags()
            from fibro2_ta import run_ta
            run_ta(log=log)
            write_findings()
            progress_snapshot("T-B QC + d0→d7 null", stop="TA_DONE")
        if args.cell in ("all", "tb"):
            write_all_prereg_flags()
            from fibro2_tb import run_tb
            run_tb(log=log)
            write_findings()
            progress_snapshot("T-C HistGradientBoosting vs ridge", stop="TB_DONE")
        if args.cell in ("all", "tc"):
            write_all_prereg_flags()
            from fibro2_tc import run_tc
            run_tc(log=log)
            write_findings()
            progress_snapshot("extrapolation distances", stop="TC_DONE")
        if args.cell in ("all", "extrap"):
            from fibro2_extrap import run_extrap
            run_extrap(log=log)
            write_findings()
        if args.cell in ("all", "findings", "prereg", "ta", "tb", "tc", "extrap"):
            p = write_findings()
            log(f"[findings] {p}")
        progress_snapshot("done", stop="FIBRO2_DONE" if args.cell == "all" else load_manifest().get("status"))
    except StopStep as e:
        log(f"STOP [{e.step}] {e.message}")
        record_failure(e.step, e.message, e.details)
        write_findings()
        progress_snapshot(f"fix STOP [{e.step}]", stop=f"STOP {e.step}: {e.message}")
        raise
    except Exception as e:
        log(f"FAIL {type(e).__name__}: {e}")
        record_failure("fail", f"{type(e).__name__}: {e}")
        write_findings()
        progress_snapshot(f"fix {type(e).__name__}", stop=f"FAIL {type(e).__name__}: {e}")
        raise
    finally:
        log.close()


if __name__ == "__main__":
    main()
