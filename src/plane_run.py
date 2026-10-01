"""Orchestrate Prompt M step 2 (age–identity plane). Does not modify existing FINDINGS*.md or FALSIFICATION.md."""
from __future__ import annotations

import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
from plane_common import (  # noqa: E402
    PLANE_DIR, StopStep, Logger, plane_log_banner, write_all_prereg_flags,
    record_failure, progress_snapshot, load_manifest, save_manifest,
)
from plane_task1 import run_task1  # noqa: E402
from plane_task2 import run_task2  # noqa: E402
from plane_task3 import run_task3  # noqa: E402
from plane_findings import write_findings  # noqa: E402


def main():
    log = Logger(PLANE_DIR / "run_report.txt")
    plane_log_banner(log, "RUN")
    write_all_prereg_flags()
    progress_snapshot(
        "Task 1 occupancy from md2 labels, sanity check, then DECLARED_BEFORE_SCORES, then medians/nulls",
        stop="PREREG_WRITTEN",
    )
    try:
        run_task1(log=log)
        run_task2(log=log)
        run_task3(log=log)
        write_findings()
        man = load_manifest()
        man["status"] = "DONE"
        save_manifest(man)
        progress_snapshot("done", stop="DONE")
        log("[run] DONE")
    except StopStep as e:
        record_failure(e.step, e.message, e.details)
        progress_snapshot(f"fix STOP [{e.step}]", stop=f"STOP {e.step}: {e.message}")
        log(f"[run] STOP [{e.step}] {e.message}")
        try:
            write_findings()
        except Exception as fe:
            log(f"[run] findings write after STOP failed: {fe}")
        raise
    finally:
        log.close()


if __name__ == "__main__":
    main()
