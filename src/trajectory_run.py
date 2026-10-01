"""Run G1 → G2 → G3 search → G3 → G4 findings, honouring STOP files.

Usage: python src/trajectory_run.py
"""
from __future__ import annotations

import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))


def main():
    import trajectory_make_notebooks  # noqa: F401  writes notebooks/trajectory_g*.ipynb
    from trajectory_g1 import run as run_g1
    from trajectory_findings import write_findings

    print("[run] G1...")
    g1 = run_g1()
    write_findings()
    if g1.get("stop_G1"):
        print("[run] STOP G1 — not running G2/G3. See FINDINGS_TRAJECTORY.md")
        return

    print("[run] G2...")
    from trajectory_g2 import run as run_g2
    g2 = run_g2()
    write_findings()
    if g2.get("stop_G2") or g2.get("skipped"):
        print("[run] STOP G2 — not running G3 trajectory. See FINDINGS_TRAJECTORY.md")
        # still run the search so the limitation/needed-dataset record exists
        print("[run] G3 search (documentation only)...")
        from trajectory_g3_search import run as run_search
        run_search()
        write_findings()
        return

    print("[run] G3 search...")
    from trajectory_g3_search import run as run_search
    run_search()
    print("[run] G3 transfer / trajectory...")
    from trajectory_g3 import run as run_g3
    run_g3()
    p = write_findings()
    print("[run] wrote", p)


if __name__ == "__main__":
    main()
