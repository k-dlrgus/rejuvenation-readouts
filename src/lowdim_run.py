"""Run L1 → L2 → L3 → findings, honouring STOP L2.

Usage: python src/lowdim_run.py
"""
from __future__ import annotations

import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))


def main():
    import lowdim_make_notebooks  # noqa: F401  writes notebooks/lowdim_l*.ipynb
    from lowdim_l1 import run as run_l1
    from lowdim_findings import write_findings

    print("[run] L1...")
    run_l1()
    write_findings()

    print("[run] L2...")
    from lowdim_l2 import run as run_l2
    l2 = run_l2()
    write_findings()
    if l2.get("stop_L2"):
        print("[run] STOP L2 — not running L3. See FINDINGS_LOWDIM.md")
        return

    print("[run] L3...")
    from lowdim_l3 import run as run_l3
    run_l3()
    p = write_findings()
    print("[run] wrote", p)


if __name__ == "__main__":
    main()
