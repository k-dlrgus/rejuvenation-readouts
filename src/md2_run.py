"""Orchestrate Prompt I (MD2). Does not modify existing FINDINGS*.md or FALSIFICATION.md."""
from __future__ import annotations

import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
from md2_common import (  # noqa: E402
    MD2_DIR, StopStep, Logger, md2_log_banner, write_all_prereg_flags,
    record_failure, progress_snapshot, load_manifest, save_manifest,
)
from md2_genesets import build_gene_sets  # noqa: E402
from md2_cluster import run_clustering  # noqa: E402
from md2_task1 import run_task1  # noqa: E402
from md2_task2 import run_task2  # noqa: E402
from md2_task3 import run_task3  # noqa: E402
from md2_findings import write_findings  # noqa: E402


def _clean_copyrighted_extracts():
    for p in (
        MD2_DIR / "paper_extract.txt",
        Path(__file__).resolve().parent / "_tmp_extract_pdf.py",
    ):
        if p.exists():
            p.unlink()


def main():
    log = Logger(MD2_DIR / "run_report.txt")
    md2_log_banner(log, "RUN")
    _clean_copyrighted_extracts()
    write_all_prereg_flags()
    progress_snapshot("build gene sets from Hallmark + mmc3 + STAR Methods", stop="PREREG_WRITTEN")
    try:
        build_gene_sets(log=log)
        run_clustering(log=log)
        run_task1(log=log)
        run_task2(log=log)
        run_task3(log=log)
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
        raise
    finally:
        log.close()


if __name__ == "__main__":
    main()
