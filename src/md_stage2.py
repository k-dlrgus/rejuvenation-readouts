"""Stage 2 — gated. Reproduce the authors' MD score only if Stage 1 recovered it exactly.

Does not approximate. Does not refit the frozen ruler. Does not open GSE325735.
"""
from __future__ import annotations

import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
from md_common import (  # noqa: E402
    MD_DIR, GATE_PATH, PREREG_STAGE2_FLAG, PREREG_STAGE2, FROZEN_RULER,
    StopStep, Logger, load_json, dump_json, jsonable, load_manifest, save_manifest,
    md_log_banner, progress_snapshot, record_failure, write_prereg_stage2,
)


def run_stage2(log=None):
    write_prereg_stage2()
    close_log = False
    if log is None:
        log = Logger(MD_DIR / "s2_report.txt")
        close_log = True
    md_log_banner(log, "STAGE2")
    log(PREREG_STAGE2)
    if not PREREG_STAGE2_FLAG.exists():
        raise StopStep("prereg", "PREREG_STAGE2.flag missing")
    if not GATE_PATH.exists():
        raise StopStep("stage2_gate", "gate.json missing — run Stage 1 first")
    gate = load_json(GATE_PATH)
    log(f"[gate] open={gate.get('open')} reason={gate.get('reason')}")
    if not gate.get("open"):
        rec = dict(
            ran=False,
            reason=gate.get("reason"),
            missing=gate.get("missing"),
        )
        dump_json(MD_DIR / "s2_summary.json", jsonable(rec))
        man = load_manifest()
        man["status"] = "STAGE2_NOT_RUN_GATE_CLOSED"
        man["s2_summary"] = rec
        save_manifest(man)
        progress_snapshot(
            "STOP. Stage 1 gate closed; Stage 2 does not run. Write FINDINGS_MD.md Stage 1 only.",
            stop="STAGE2_NOT_RUN_GATE_CLOSED",
        )
        log("[stage2] not run — gate closed")
        if close_log and hasattr(log, "close"):
            log.close()
        return rec
    if not FROZEN_RULER.exists():
        raise StopStep("frozen_ruler", f"missing {FROZEN_RULER}")
    # If the gate opens, the recovered gene list and computation must be loaded from
    # results/md/ (not reconstructed from memory) and implemented here. That branch
    # is not written until Stage 1 actually recovers both.
    raise StopStep(
        "stage2_impl",
        "Gate is open but Stage 2 implementation was not generated from a recovered "
        "gene list/computation in this session. Refusing to approximate.",
    )


if __name__ == "__main__":
    run_stage2()
