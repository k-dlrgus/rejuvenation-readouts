"""Run V1 → V2 → V3 → V4 → findings, honouring STOP files.

Usage: python src/target_run.py
"""
from __future__ import annotations

import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))


def main():
    import target_make_notebooks  # noqa: F401  writes notebooks/target_v*.ipynb
    from target_v1 import run as run_v1
    from target_findings import write_findings

    print("[run] V1...")
    v1 = run_v1()
    write_findings()
    if v1.get("stop_V1"):
        print("[run] STOP V1 — not running V2–V4. See FINDINGS_TARGET.md")
        return

    print("[run] V2...")
    from target_v2 import run as run_v2
    v2 = run_v2()
    write_findings()
    if v2.get("stop_V2") or v2.get("skipped"):
        print("[run] STOP V2 — not running V3–V4. See FINDINGS_TARGET.md")
        return

    print("[run] V3...")
    from target_v3 import run as run_v3
    run_v3()
    write_findings()

    print("[run] V4...")
    from target_v4 import run as run_v4
    run_v4()
    p = write_findings()
    print("[run] wrote", p)


if __name__ == "__main__":
    main()
