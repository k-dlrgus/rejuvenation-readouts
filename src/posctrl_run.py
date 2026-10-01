"""Orchestrate POSCTRL: C2a gate, then C2b, C1, C3, findings.

Usage: python src/posctrl_run.py
Resumes from results/posctrl/*.json (complete=true skipped).
"""
from __future__ import annotations

import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
from posctrl_common import (  # noqa: E402
    PC_DIR, PC_SEED, progress_snapshot, load_json, cell_path, PREREG_BLOCK,
    TARGET_R2_PRIMARY, REGIMES, SUPERSESSION_FLAG, refresh_progress,
)
from posctrl_c2 import run as run_c2
from posctrl_c1 import run as run_c1, cell_name
from posctrl_c3 import run as run_c3
from posctrl_findings import write_findings, _load_if


def _collect():
    cells = {}
    for p in sorted(PC_DIR.glob("*.json")):
        if p.name.endswith(("_battery.json", "_bootstrap.json", "_within_site.json",
                            "_loso.json", "_youngold.json", "_setup.json")):
            continue
        try:
            rec = load_json(p)
        except Exception:
            continue
        if rec.get("complete") or rec.get("boot_median_angle") is not None:
            cells[p.stem] = rec
            rec["path"] = str(p)
    return cells


def _betas():
    p = PC_DIR / "frozen_beta.json"
    if not p.exists():
        return None
    rec = load_json(p)
    return {k: v.get("beta") if isinstance(v, dict) else v for k, v in rec.items()}


def main():
    PC_DIR.mkdir(parents=True, exist_ok=True)
    progress_snapshot(_collect(), "run C2a (gate)", betas=_betas())
    print("[posctrl] C2a", flush=True)
    c2 = run_c2("c2a")
    c2a = (c2 or {}).get("c2a") or _load_if("c2a")
    if c2a and c2a.get("p1_fail") and not SUPERSESSION_FLAG.exists():
        progress_snapshot(_collect(), "STOP P1 — write findings only",
                          stop="P1", betas=_betas())
        write_findings()
        print("[posctrl] STOP P1", flush=True)
        return 2
    if c2a and c2a.get("p1_fail") and SUPERSESSION_FLAG.exists():
        print("[posctrl] P1 superseded 2026-09-16; continuing to C1/C3", flush=True)
        refresh_progress()
    progress_snapshot(_collect(), "run C2b", betas=_betas())
    print("[posctrl] C2b", flush=True)
    run_c2("c2b")
    progress_snapshot(_collect(), "run C1 dense_shared 0.30 (primary)", betas=_betas())
    print("[posctrl] C1 primary dense_shared 0.30", flush=True)
    run_c1(["posctrl_c1.py", "dense_shared", "0.30"])
    progress_snapshot(_collect(), "run remaining C1 regimes @ 0.30", betas=_betas())
    for r in REGIMES:
        if r == "dense_shared":
            continue
        print(f"[posctrl] C1 {r} 0.30", flush=True)
        run_c1(["posctrl_c1.py", r, "0.30"])
        progress_snapshot(_collect(), f"next C1 after {r}", betas=_betas())
    print("[posctrl] C1 dense_shared 0.15 / 0.50 (B/2)", flush=True)
    run_c1(["posctrl_c1.py", "dense_shared", "0.15"])
    run_c1(["posctrl_c1.py", "dense_shared", "0.50"])
    progress_snapshot(_collect(), "run C3 real-age then planted", betas=_betas())
    print("[posctrl] C3", flush=True)
    run_c3()
    progress_snapshot(_collect(), "write FINDINGS_POSCTRL.md", betas=_betas())
    fired = write_findings()
    progress_snapshot(_collect(), "done", stop=fired, betas=_betas(),
                      extra=f"Fired: {fired}")
    print(f"[posctrl] done fired={fired}", flush=True)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
