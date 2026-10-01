"""Run Prompt G: Task 4 supersession, joint-cluster trajectory, extrap, per-cluster null.

Usage: python src/fibro3_run.py [--cell prereg|t4|t1|t2|t3|findings|all]

Writes pre-registration flags and FINDINGS_FIBRO3.md skeleton BEFORE any scoring.
Does not modify existing FINDINGS*.md or FALSIFICATION.md.
Does not open GSE325735. Does not promote d7. Does not run Task 4 option (a).
"""
from __future__ import annotations

import argparse
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
from fibro3_common import (  # noqa: E402
    FIBRO3_DIR, FROZEN_RULER, STANDING_CONSTRAINT,
    write_all_prereg_flags, progress_snapshot, load_manifest, save_manifest,
    fibro3_log_banner, Logger, StopStep, record_failure,
)
from fibro3_findings import write_findings  # noqa: E402


def _prereg_and_skeleton():
    write_all_prereg_flags()
    man = load_manifest()
    man["status"] = "PREREG_WRITTEN"
    man["frozen_ruler"] = str(FROZEN_RULER)
    man["task4_choice"] = "b"
    man["option_a_run"] = False
    save_manifest(man)
    progress_snapshot(
        "Task 4 score superseded primary CELLxGENE 10x (option a not run)",
        stop="PREREG_WRITTEN",
        extra="Task 4/1/2/3 pre-registered blocks and CLUSTER_SETTINGS written before any score. "
              + STANDING_CONSTRAINT,
    )
    write_findings()
    return "prereg written; FINDINGS_FIBRO3.md skeleton on disk"


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--cell", default="all",
                    help="prereg | t4 | t1 | t2 | t3 | findings | all")
    args = ap.parse_args()
    log = Logger(FIBRO3_DIR / "run_report.txt", mode="a" if args.cell != "all" else "w")
    fibro3_log_banner(log, f"RUN cell={args.cell}")
    try:
        if args.cell in ("all", "prereg"):
            log(_prereg_and_skeleton())
        if args.cell in ("all", "t4"):
            write_all_prereg_flags()
            from fibro3_task4 import run_task4
            run_task4(log=log)
            write_findings()
            progress_snapshot("Task 1 joint clustering (settings already frozen)", stop="TASK4_DONE")
        if args.cell in ("all", "t1"):
            write_all_prereg_flags()
            from fibro3_task1 import run_task1
            run_task1(log=log)
            write_findings()
            progress_snapshot("Task 2 extrapolation distances", stop="TASK1_DONE")
        if args.cell in ("all", "t2"):
            write_all_prereg_flags()
            from fibro3_task2 import run_task2
            run_task2(log=log)
            write_findings()
            progress_snapshot("Task 3 per-cluster random-direction null", stop="TASK2_DONE")
        if args.cell in ("all", "t3"):
            write_all_prereg_flags()
            from fibro3_task3 import run_task3
            run_task3(log=log)
            write_findings()
        if args.cell in ("all", "findings", "prereg", "t4", "t1", "t2", "t3"):
            p = write_findings()
            log(f"[findings] {p}")
        man = load_manifest()
        next_map = {
            "prereg": ("Task 4 score superseded primary CELLxGENE 10x (option a not run)", "PREREG_WRITTEN"),
            "t4": ("Task 1 joint clustering (settings already frozen)", "TASK4_DONE"),
            "t1": ("Task 2 extrapolation distances", "TASK1_DONE"),
            "t2": ("Task 3 per-cluster random-direction null", "TASK2_DONE"),
            "t3": ("write FINDINGS_FIBRO3.md from disk", "TASK3_DONE"),
            "findings": ("done", man.get("status")),
        }
        if args.cell == "all" and man.get("status") != "STOP":
            man["status"] = "FIBRO3_DONE"
            save_manifest(man)
            progress_snapshot("done", stop="FIBRO3_DONE")
        elif args.cell in next_map:
            nxt, st = next_map[args.cell]
            progress_snapshot(nxt, stop=st)
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
