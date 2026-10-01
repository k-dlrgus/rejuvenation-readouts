"""Run the geometry audit: A, then B/C only if STOP A is clear; write FINDINGS_GEOMETRY.md.

Usage: python src/geometry_run.py
"""
from __future__ import annotations

import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))


def main():
    import geometry_make_notebooks  # writes notebooks/geometry_part*.ipynb
    from geometry_partA import run as run_A
    from geometry_findings import write_findings

    print("[run] Part A...")
    a = run_A()
    p = write_findings()
    if a.get("stop_A"):
        print("[run] STOP A - not running B or C. See FINDINGS_GEOMETRY.md")
        return
    print("[run] Part B...")
    from geometry_partB import run as run_B
    run_B()
    write_findings()
    print("[run] Part C...")
    from geometry_partC import run as run_C
    run_C()
    p = write_findings()
    print("[run] wrote FINDINGS_GEOMETRY.md")


if __name__ == "__main__":
    main()
