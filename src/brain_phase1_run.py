"""Run DLPFC Phase 1 stages in order. Halt on STOP P2 / STOP P3.

Usage:
  python src/brain_phase1_run.py
  python src/brain_phase1_run.py --stage p0
"""
from __future__ import annotations

import argparse
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--stage", default="all",
                    choices=["all", "p0", "p1", "p2", "p3", "p4", "p5"])
    args = ap.parse_args()
    order = ["p0", "p1", "p2", "p3", "p4", "p5"] if args.stage == "all" else [args.stage]
    for st in order:
        print("\n" + "#" * 80 + f"\n# PHASE 1  {st.upper()}\n" + "#" * 80, flush=True)
        if st == "p0":
            from brain_phase1_p0 import run
            run()
        elif st == "p1":
            from brain_phase1_p1 import run
            run()
        elif st == "p2":
            from brain_phase1_p2 import run
            rec = run()
            if rec.get("stop") and args.stage == "all":
                print("[run] STOP P2 — not continuing to P3–P5", flush=True)
                return 2
        elif st == "p3":
            from brain_phase1_p3 import run
            rec = run()
            if rec.get("stop") and args.stage == "all":
                print("[run] STOP P3 — not continuing to P4–P5", flush=True)
                return 3
        elif st == "p4":
            from brain_phase1_p4 import run
            rec = run()
            # P4 verdict of 'not separable' is a result, not a crash; still run P5.
        elif st == "p5":
            from brain_phase1_p5 import run
            run()
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
