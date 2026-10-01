"""Run Prompt H: reconstruct the authors' MD claim from public sources, then (if gated) test it.

Usage: python src/md_run.py [--cell prereg|stage1|stage2|findings|all]

Writes PREREG_STAGE2.flag BEFORE any MD score.
Does not modify existing FINDINGS*.md or FALSIFICATION.md.
Does not open GSE325735. Does not circumvent a paywall. Does not approximate an MD score.
"""
from __future__ import annotations

import argparse
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
from md_common import (  # noqa: E402
    MD_DIR, FROZEN_RULER, write_prereg_stage2, progress_snapshot,
    load_manifest, save_manifest, md_log_banner, Logger, StopStep, record_failure,
)
from md_findings import write_findings  # noqa: E402


def _prereg():
    write_prereg_stage2()
    man = load_manifest()
    man["status"] = "PREREG_WRITTEN"
    man["frozen_ruler"] = str(FROZEN_RULER)
    save_manifest(man)
    progress_snapshot(
        "Stage 1 public-source retrieval (GEO, PubMed, PMC, preprint, code, open supplements)",
        stop="PREREG_WRITTEN",
        extra="PREREG_STAGE2 written before any MD score.",
    )
    write_findings()
    return "prereg written"


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--cell", default="all", help="prereg | stage1 | stage2 | findings | all")
    args = ap.parse_args()
    log = Logger(MD_DIR / "run_report.txt", mode="a" if args.cell != "all" else "w")
    md_log_banner(log, f"RUN cell={args.cell}")
    try:
        if args.cell in ("all", "prereg"):
            log(_prereg())
        if args.cell in ("all", "stage1"):
            write_prereg_stage2()
            from md_stage1 import run_stage1
            run_stage1(log=log)
            write_findings()
        if args.cell in ("all", "stage2"):
            write_prereg_stage2()
            from md_stage2 import run_stage2
            run_stage2(log=log)
            write_findings()
        if args.cell in ("all", "findings", "prereg", "stage1", "stage2"):
            p = write_findings()
            log(f"[findings] {p}")
    except StopStep as e:
        record_failure(e.step, e.message)
        log(f"STOP [{e.step}] {e.message}")
        progress_snapshot(f"stopped at {e.step}", stop="STOP")
        write_findings()
        raise SystemExit(2)
    finally:
        log.close()


if __name__ == "__main__":
    main()
