"""Run the fibroblast age-ruler + GSE297234 pipeline.

Usage: python src/fibro_run.py [--stage 1|2|findings|all]
Stage 2 runs only if Stage 1 passed. Does not modify existing FINDINGS*.md
or FALSIFICATION.md. Does not open GSE325735.
"""
from __future__ import annotations

import argparse
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
from fibro_common import (  # noqa: E402
    FIBRO_DIR, PREREG_STAGE1_FLAG, write_prereg_stage1, progress_snapshot,
    load_json, fibro_log_banner, Logger,
)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--stage", default="all", help="1 | 2 | findings | all")
    args = ap.parse_args()
    write_prereg_stage1()
    progress_snapshot("Stage 1 cohort" if args.stage in ("all", "1") else "as requested",
                      stop="starting")
    if args.stage in ("all", "1"):
        from fibro_stage1 import main as s1
        s1()
    gate_p = FIBRO_DIR / "stage1_gate.json"
    passed = bool(gate_p.exists() and load_json(gate_p).get("pass"))
    if args.stage in ("all", "2"):
        if not passed:
            progress_snapshot(
                "write FINDINGS_FIBRO.md; do not open GSE297234",
                stop="STAGE1_FAIL",
            )
            if args.stage == "2":
                raise SystemExit("Stage 1 did not pass; Stage 2 not run")
        else:
            from fibro_stage2 import main as s2
            s2()
    if args.stage in ("all", "findings", "1", "2"):
        from fibro_findings import write_findings
        p = write_findings()
        print(p)


if __name__ == "__main__":
    main()
