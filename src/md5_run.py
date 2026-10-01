"""Orchestrate Prompt L (MD5). Does not modify existing FINDINGS*.md or FALSIFICATION.md."""
from __future__ import annotations

import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
from md5_common import (  # noqa: E402
    MD5_DIR, StopStep, Logger, md5_log_banner, write_all_prereg_flags,
    record_failure, progress_snapshot, load_manifest, save_manifest,
)
from md5_task1 import run_task1  # noqa: E402
from md5_task2 import run_task2  # noqa: E402
from md5_findings import write_findings  # noqa: E402


def main():
    log = Logger(MD5_DIR / "run_report.txt")
    md5_log_banner(log, "RUN")
    write_all_prereg_flags()
    progress_snapshot(
        "Task 1 occupancy from md2 labels, sanity check, then md4 cell_scores contrasts",
        stop="PREREG_WRITTEN",
    )
    try:
        run_task1(log=log)
        run_task2(log=log)
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
