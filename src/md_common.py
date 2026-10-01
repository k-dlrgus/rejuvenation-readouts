"""Prompt H — reconstruct Lu et al. MD claim from public sources, then (if gated) test it.

Namespaced under results/md/. Does not modify any existing FINDINGS_*.md or FALSIFICATION.md.
Does not open GSE325735. Does not circumvent a paywall.
Seed 20260914. Bootstrap seed 20260918. Spearman ρ is the reported statistic.
Refit nothing. Frozen ruler used as-is.
"""
from __future__ import annotations

import json
import sys
from datetime import date
from pathlib import Path

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parent))
from config import RESULTS, ROOT, SEED  # noqa: E402
from brain_phase1_common import Logger, dump_json, load_json  # noqa: E402
from gtex_common import StopStep  # noqa: E402
from fibro_common import (  # noqa: E402
    FROZEN_RULER, FIBRO_SEED, BOOT_SEED, N_PERM, N_BOOT, N_RANDOM_DIR, jsonable,
)
from fibro3_common import MIN_CELLS, AGED_LINE, YOUNG_LINE  # noqa: E402

MD_DIR = RESULTS / "md"
MD_RAW = MD_DIR / "raw"
MD_FIG = MD_DIR / "figures"
for _p in (MD_DIR, MD_RAW, MD_FIG):
    _p.mkdir(parents=True, exist_ok=True)

MD_SEED = int(SEED)  # 20260914
MD_BOOT = int(BOOT_SEED)  # 20260918
RETRIEVAL_DATE = date(2026, 9, 19).isoformat()
PMID = "40816266"
DOI = "10.1016/j.cell.2025.07.031"
PII = "S0092-8674(25)00853-0"
PII_SD = "S0092867425008530"
GSE_SC = "GSE297234"
GSE_BULK = "GSE297233"
TITLE = "Prevalent mesenchymal drift in aging and disease is reversed by partial reprogramming"

MANIFEST_PATH = MD_DIR / "manifest.json"
PROGRESS_PATH = ROOT / "PROGRESS_MD.md"
FINDINGS_PATH = ROOT / "FINDINGS_MD.md"
PREREG_STAGE2_FLAG = MD_DIR / "PREREG_STAGE2.flag"
GATE_PATH = MD_DIR / "gate.json"
SOURCES_CSV = MD_DIR / "sources.csv"
SOURCES_JSON = MD_DIR / "sources.json"

# Public-source allowlist. Third-party mirrors of the paywalled PDF are excluded.
DISALLOWED_HOST_SUBSTR = (
    "rapamycin.news",
    "sci-hub",
    "scihub",
    "libgen",
    "annas-archive",
    "reddit.com",
)

PREREG_STAGE2 = """Written **before any MD score or frozen-ruler comparison**, 2026-09-19. No threshold in this block is re-tuned after numbers exist. Stage 2 runs only if Stage 1 recovered the authors' MD gene list and computation from public sources (GEO, PubMed, PMC, preprint servers, code repositories, open supplements). Exact or not at all. Refit nothing. Recalibrate nothing. Rescale nothing. Frozen ruler `results/fibro/frozen_ruler_ridge_raw.npz` used as-is. Clustering and cell filters already frozen in `results/fibro3/` used as-is. GSE325735 is out of scope.

On the GSE297234 matrices already used by FINDINGS_FIBRO / FINDINGS_FIBRO3 (Cell Ranger `filtered_feature_bc_matrix.h5`), using the joint-cluster labels and MIN_CELLS=20 already frozen in `results/fibro3/`:

1. **Reproduce their metric.** Compute the MD score as published, per cell and per (donor × cluster × timepoint). Report its trajectory. State plainly whether our reproduction reproduces their reported direction and magnitude. If it does not, say so and stop — a failed reproduction of their metric means Stage 2's comparison is not valid.
2. **Run their comparison both ways.** For the exact timepoint comparison they made, report side by side: their MD score, and our frozen age ruler. Both on whole-sample averages, and both per cluster. Four cells. Never averaged together.
3. **Cross-check the metrics against each other.** Across all donor × cluster × timepoint rows with ≥20 cells, report Spearman ρ between MD score and frozen age score, with permutation null and bootstrap CI. Two age-related metrics on the same cells should agree if they measure the same thing.
4. **Composition test on their metric.** Apply the FINDINGS_FIBRO3.md per-cluster decomposition to *their* score: does the MD change survive when populations are analysed separately, or is it a mix shift? Report with the same 200-direction null where applicable.

Seed `20260914`. boot `20260918`. n_perm=200. n_boot=200. n_random=200. Spearman ρ with permutation null is the reported statistic. Nothing averaged across donors, clusters, or regimes.

**Pre-registered reading (only the outcome that fired):**
- MD moves per cluster, our ruler does not → the two instruments disagree on the same cells. Report as an instrument disagreement, describe both, and do not declare a winner. Name the test that would settle it.
- Neither moves per cluster, both move on whole-sample averages → the published claim rests on composition, as ours did before decomposition. This is a methodological finding about how reprogramming time courses are read, and it must be stated carefully and without overreach: one dataset, two donors, our reproduction of their metric.
- Both move per cluster → their claim holds and our ruler is insensitive to what OSKM does in these cells. Report it as a limitation of our ruler. This is a real possible outcome and must not be argued away.
- MD and our ruler correlate strongly across rows but disagree on the trajectory → report the discrepancy as unexplained.
"""

STAGE2_VERDICT_FIBRO = (
    "age score does not decline → the ruler reads a static donor property, not a "
    "modifiable state. Step 3 is not supported by this data."
)
FIBRO3_TASK1_SENTENCE = (
    "Cell-intrinsic. The same cells' clusters individually reverse from d7 to d10. "
    "Then reprogramming genuinely pushes cells back along the age axis after d7, "
    "and that is a biological claim worth its own test. clusters=[0, 1, 2]."
)
FIBRO3_TASK3_AGED_D0D7 = (
    "Aged GM00731 d0→d7 (the all-cell T-B cell had p=+0.005): cluster 0 decline=-0.582 "
    "p=+0.652 n_ge=130/200 n_d0=4941 n_d7=2719; cluster 1 decline=+2.215 p=+0.144 "
    "n_ge=28/200 n_d0=76 n_d7=850; cluster 2 skipped (n_cells<20 at an endpoint)."
)


def load_manifest():
    if MANIFEST_PATH.exists():
        return load_json(MANIFEST_PATH)
    return dict(status="INIT", failures=[], columns={})


def save_manifest(man):
    dump_json(MANIFEST_PATH, jsonable(man))
    return MANIFEST_PATH


def record_failure(step, message, extra=None):
    man = load_manifest()
    rec = dict(step=step, message=str(message))
    if extra:
        rec["extra"] = jsonable(extra)
    fails = man.setdefault("failures", [])
    if not any(f.get("step") == step and f.get("message") == rec["message"] for f in fails):
        fails.append(rec)
    man["status"] = "STOP"
    save_manifest(man)
    return rec


def write_prereg_stage2():
    if PREREG_STAGE2_FLAG.exists():
        existing = PREREG_STAGE2_FLAG.read_text(encoding="utf-8")
        if existing.strip() != PREREG_STAGE2.strip():
            raise StopStep(
                "prereg_stage2",
                "PREREG_STAGE2.flag already exists and does not match the frozen block. Not rewriting.",
            )
        return PREREG_STAGE2_FLAG
    PREREG_STAGE2_FLAG.write_text(PREREG_STAGE2, encoding="utf-8")
    return PREREG_STAGE2_FLAG


def md_log_banner(log, stage: str):
    log("=" * 100)
    log(f"MD {stage}  seed={MD_SEED}  boot_seed={MD_BOOT}  retrieval_date={RETRIEVAL_DATE}")
    log("=" * 100)
    log("FALSIFICATION.md and existing FINDINGS*.md are not modified.")
    log("Public sources only: GEO, PubMed, PMC, preprint servers, code repositories, open supplements.")
    log("Do not circumvent a paywall. Do not approximate an MD score. Do not open GSE325735.")
    log("Refit nothing. Frozen ruler used as-is.")
    log(f"frozen_ruler exists={FROZEN_RULER.exists()} path={FROZEN_RULER}")


def _read_csv_if(path: Path):
    if path.exists():
        return pd.read_csv(path)
    return None


def _disk_row_lines(title, path, max_cols=14, max_rows=40):
    df = _read_csv_if(path)
    out = [f"### {title} (`{path.name}`)"]
    if df is None or not len(df):
        out.append("- missing or empty")
        out.append("")
        return out
    cols = list(df.columns)[:max_cols]
    for _, r in df.head(max_rows).iterrows():
        bits = [f"{c}={r[c]}" for c in cols]
        out.append("- " + " ".join(bits))
    if len(df) > max_rows:
        out.append(f"- … {len(df) - max_rows} more rows")
    out.append("")
    return out


def progress_snapshot(next_action, stop=None, extra=""):
    man = load_manifest() if MANIFEST_PATH.exists() else {}
    stop_s = stop if stop else man.get("status", "running")
    searches = man.get("searches_run") or []
    lines = [
        "# PROGRESS_MD",
        "",
        f"**STOP status:** {stop_s}",
        f"**Next action:** {next_action}",
        "",
        "## Seeds / gates",
        "",
        f"- seed `{MD_SEED}`  boot `{MD_BOOT}`  n_perm={N_PERM}  n_boot={N_BOOT}  n_random={N_RANDOM_DIR}",
        "- primary metric: Spearman ρ with permutation null; uncalibrated R² is never a gate",
        f"- frozen ruler: `{FROZEN_RULER}` exists={FROZEN_RULER.exists()}",
        f"- PREREG_STAGE2.flag exists={PREREG_STAGE2_FLAG.exists()}",
        f"- gate.json exists={GATE_PATH.exists()}",
        "- GSE325735: out of scope",
        "- Refit nothing. Do not circumvent a paywall.",
        "",
        "## Gate status",
        "",
    ]
    if GATE_PATH.exists():
        gate = load_json(GATE_PATH)
        lines.append(f"- open={gate.get('open')}  reason={gate.get('reason')}")
        lines.append(f"- gene_list_recovered={gate.get('gene_list_recovered')}")
        lines.append(f"- computation_recovered={gate.get('computation_recovered')}")
        lines.append(f"- missing={gate.get('missing')}")
        lines.append("")
    else:
        lines.append("- not evaluated")
        lines.append("")
    lines += [
        "## Searches already run (do not repeat)",
        "",
    ]
    if searches:
        for s in searches:
            lines.append(
                f"- `{s.get('query')}` via {s.get('via')} → {s.get('result', '')[:400]}"
            )
        lines.append("")
    else:
        lines.append("- (none recorded in manifest yet)")
        lines.append("")
    lines += ["## Finished cells (from disk)", ""]
    for name, path in (
        ("sources", SOURCES_CSV),
        ("Stage 2 MD trajectory all-cell", MD_DIR / "s2_md_allcell.csv"),
        ("Stage 2 MD trajectory cluster", MD_DIR / "s2_md_cluster.csv"),
        ("Stage 2 four-cell comparison", MD_DIR / "s2_four_cell.csv"),
        ("Stage 2 MD vs age ρ", MD_DIR / "s2_md_vs_age_rho.csv"),
        ("Stage 2 composition", MD_DIR / "s2_composition.csv"),
        ("Stage 2 random-direction null", MD_DIR / "s2_null.csv"),
    ):
        lines += _disk_row_lines(name, path)
    for js_name in (
        "gate.json", "s1_summary.json", "s2_summary.json", "s2_reading.json",
        "mmc3_MD_signatures.json", "geo_gds_two_series.json",
    ):
        p = MD_DIR / js_name
        if p.exists():
            rec = load_json(p)
            if js_name == "mmc3_MD_signatures.json" and isinstance(rec, dict):
                rec = {
                    "url": rec.get("url"),
                    "file": rec.get("file"),
                    "sheets": rec.get("sheets"),
                    "n_MD_score": rec.get("n_MD_score"),
                    "n_TGFB_score": rec.get("n_TGFB_score"),
                    "gene_lists_omitted_from_progress": True,
                }
            lines.append(f"### {js_name}")
            lines.append(f"- {json.dumps(jsonable(rec), ensure_ascii=False)[:4000]}")
            lines.append("")
    if extra:
        lines += ["## Note", "", extra, ""]
    if man.get("failures"):
        lines += ["## Failures (manifest)", ""]
        for f in man["failures"]:
            lines.append(f"- **{f.get('step')}:** {f.get('message')}")
        lines.append("")
    lines += [
        "## Files",
        "",
        "- `src/md_common.py`, `md_stage1.py`, `md_stage2.py`, `md_findings.py`, `md_run.py`, `md_finalize_stage1.py`",
        "- `results/md/`",
        "- `FINDINGS_MD.md`",
        "- `PROGRESS_MD.md`",
        "",
    ]
    PROGRESS_PATH.write_text("\n".join(lines) + "\n", encoding="utf-8")
    return PROGRESS_PATH
