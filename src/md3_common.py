"""Prompt J — test the ruler on the claim population; compare instruments vs donor age.

Namespaced under results/md3/. Does not modify any existing FINDINGS_*.md or FALSIFICATION.md.
Does not open GSE325735. Refit nothing. Frozen ruler used as-is.
Seed 20260914. Bootstrap seed 20260918. Spearman ρ is the reported statistic.
Reuse src/md2_*.py Louvain + AddModuleScore exactly; do not re-tune.
"""
from __future__ import annotations

import json
import sys
from pathlib import Path

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parent))
from config import RESULTS, DATA_RAW, DATA_PROC, ROOT, SEED  # noqa: E402
from brain_phase1_common import Logger, dump_json, load_json  # noqa: E402
from gtex_common import StopStep, spearman_safe  # noqa: E402
from fibro_common import (  # noqa: E402
    FROZEN_RULER, BOOT_SEED, N_PERM, N_BOOT, N_RANDOM_DIR,
    jsonable, bootstrap_rho_ci, PLURI_ENDOGENOUS, FIBRO_IDENTITY, OSKM_FAMILY,
)
from fibro3_common import MIN_CELLS, AGED_LINE, YOUNG_LINE, DAYS  # noqa: E402
from fibro2_common import load_frozen_ruler, rho_null_donor, ADULT_MIN_AGE  # noqa: E402
from md2_common import (  # noqa: E402
    MD2_DIR, MD2_PROC, MMC3_XLSX, GSE113957_FPKM, CXG_H5AD, CXG_DATASET_ID,
    GSE226189_TAR, AMS_NBIN, AMS_CTRL, SEURAT_LOGNORMALIZE_SCALE,
    LOUVAIN_RES, CLUSTER_SETTINGS,
)

MD3_DIR = RESULTS / "md3"
MD3_FIG = MD3_DIR / "figures"
MD3_PROC = DATA_PROC / "md3"
for _p in (MD3_DIR, MD3_FIG, MD3_PROC):
    _p.mkdir(parents=True, exist_ok=True)

MD3_SEED = int(SEED)  # 20260914
MD3_BOOT = int(BOOT_SEED)  # 20260918
PREREG_DATE = "2026-09-19"
FIBRO2_DIR = RESULTS / "fibro2"

MANIFEST_PATH = MD3_DIR / "manifest.json"
PROGRESS_PATH = ROOT / "PROGRESS_MD3.md"
FINDINGS_PATH = ROOT / "FINDINGS_MD3.md"
PREREG_TASK1_FLAG = MD3_DIR / "PREREG_TASK1.flag"
PREREG_TASK2_FLAG = MD3_DIR / "PREREG_TASK2.flag"
PREREG_TASK3_FLAG = MD3_DIR / "PREREG_TASK3.flag"
DECLARED_BEFORE_SCORES_FLAG = MD3_DIR / "DECLARED_BEFORE_SCORES.flag"

MIN_CELLS_20 = int(MIN_CELLS)  # 20; FINDINGS_MD2.md threshold, unchanged
MIN_CELLS_10 = 10  # sensitivity only; never primary
MIN_CELLS_POOLED = 1  # pooled state has no extra threshold; n is reported

STATE_NAMES = ("Fibroblast", "PartialReprog", "EarlyPluripotency", "Pluripotency", "NonReprog")
# Prompt's "EarlyPluri" is mmc3 column EarlyPluripotency (md2 name, reused).
STATE_ORDINAL = {
    "Fibroblast": 0,
    "PartialReprog": 1,
    "EarlyPluripotency": 2,
    "Pluripotency": 3,
}
# NonReprog is off the OSKM progression axis and is excluded from the ordinal ρ.

TABLE_S3_REASONING = """## Why Table S3 age-up / age-down lists are exact, not an approximation

FINDINGS_MD2.md did not run the authors' aging signature because GSE113957 is FPKM
and the STAR Methods DESeq2 old-versus-young split (adjusted p < 0.05, log2FC > 0.5,
old/young from PCA separation, DESeq2 v1.40.2) could not be rebuilt from integer
counts. That refusal still stands: this prompt does not rebuild the DE analysis.

mmc3.xlsx Table S3 (sheet `Aging_signatures`, columns `Age up` and `Age down`) is
the authors' own published output of that DESeq2 analysis. Scoring those lists is
using the published result, not reconstructing the test that produced it. That is
therefore permitted here, where approximating the DE from FPKM was not.

Mismatch versus a rebuilt DE, if one were later released, would be reported; it is
not a reason to withhold the published lists.
"""

PREREG_TASK1 = """Written **before any MD3 age / MD / TGF-β / Table S3 / pluripotency score**, 2026-09-19. No threshold in this block is re-tuned after numbers exist.

Using the src/md2_*.py Louvain clustering already computed (19 clusters per donor; not re-run, not re-tuned). Cell-state labels are the md2 assignment: argmax of mean AddModuleScore of mmc3 Reprog_cell_state_signatures columns Fibroblast / PartialReprog / EarlyPluripotency / Pluripotency / NonReprog. Labels are read from results/md2/t2_cluster_labels_*.csv. Not re-labelled.

Three aggregation levels, each reported separately and never averaged together:
- (a) pooled by cell state — all PartialReprog clusters pooled into one pseudobulk per timepoint, likewise each other state. This is the level their Figure 3G operates at and is the **primary cell**. No MIN_CELLS cutoff is applied to the pool (pooling is the reason this prompt exists); n cells at each timepoint is reported. A state with n=0 at an endpoint is not pairable.
- (b) per cluster at ≥20 cells (the FINDINGS_MD2.md threshold, reproduced unchanged).
- (c) per cluster at ≥10 cells, reported as a relaxed-threshold sensitivity check and labelled as such. Not primary.

Per donor, never averaged: frozen ruler score, MD score, TGF-β score, Table S3 age-up score, Table S3 age-down score, and pluripotency−fibroblast score (frozen FINDINGS_FIBRO.md lists, same differentiation_score as md2), at d0/d3/d7/d10, with d0→d7 and d0→d10 endpoint changes.

Nulls (not interchangeable):
- Frozen ruler: 200 permuted-ruler-weight directions, seed 20260914. Empirical p = (n_random ≥ real + 1) / (n_random + 1). Endpoint decline = score(day 0) − score(day T). Positive = younger on the ruler.
- Every gene-list score (MD, TGF-β, Age up, Age down): 200 size-matched, expression-bin-matched random gene sets, AddModuleScore nbin=24 ctrl=100, seed 20260914 with a documented per-list offset so MD uses the same stream as md2. Not permuted ruler weights.
- Pluripotency−fibroblast is the frozen Stage 2 mean-z difference on the same TMM/z as the ruler, not AddModuleScore. Its null is 200 size-matched random gene sets in frozen-ruler gene space on that Z (not AMS bins). The difference is stated.

Donors GM00731 and GM23815 reported separately. d7 is not the Stage 2 endpoint. Frozen ruler results/fibro/frozen_ruler_ridge_raw.npz used as-is. Refit nothing. Louvain resolution, AddModuleScore parameters, and cell filters are inherited from md2 and are not re-tuned.

**Pre-registered reading (only the outcome that fired):**
- Ruler declines in the pooled PartialReprog state with p ≤ 0.05 → the FINDINGS_MD2.md `negative_holds` reading was limited by cell counts, not by the ruler. Say so plainly and quote the sentence superseded.
- Ruler does not decline in the pooled PartialReprog state while MD does, both with their own nulls → the instruments disagree on exactly the population the published claim concerns. This is the sharpest form of the disagreement; report both numbers side by side and declare no winner.
- Neither moves in the pooled state → our reproduction does not recover their Figure 3G at the pooled level; report the reproduction failure as the finding, not as a refutation.
- Pooled and per-cluster (≥20) levels disagree → report both, say which is which, do not pick.

“Declines” / “moves” = endpoint decline d0→d7 or d0→d10 with p ≤ 0.05 on that instrument's own null, on the aged donor GM00731 pooled PartialReprog state for (a), or on at least one aged-donor Louvain cluster for (b). GM23815 is tabulated as a contrast, never pooled, never the reading.
"""

PREREG_TASK2 = """Written **before any MD3 three-way ρ vs donor age or paired Δρ**, 2026-09-19. No threshold in this block is re-tuned after numbers exist.

Three instruments, all applied unchanged to the same donors: (1) frozen GTEx ruler from results/fibro/frozen_ruler_ridge_raw.npz, used as-is, not refit; (2) MD score, AddModuleScore as in md2; (3) their age-up minus age-down score (mmc3 Table S3 Age up and Age down, AddModuleScore, same implementation as md2). Instrument (3) is AMS(Age up) − AMS(Age down).

Cohorts:
- CELLxGENE a19d1667-a7b5-4556-9e5f-f9bfa690c0f1 (10x, n=93) — already in results/fibro2/ and results/md2/. Frozen-ruler and MD numbers are read from those files and recomputed from the saved per-donor scores as a check. Age-up minus age-down is newly scored on the cached donor-count matrix.
- GSE226189 (bulk, n=82) — same: ruler and MD from disk; age-up minus age-down newly scored. geneCOUNT Tracking_ID is Ensembl; ID type is checked before mapping.
- GSE113957 (bulk, n=143, ages 1–92). FPKM blocks DESeq2 and TMM/log2-CPM but does not block rank-based scoring. Rank-transform genes within sample and score all three instruments on ranks. This is a different transform from the GTEx spec and is therefore a weaker test. Restrict to adults ≥18 to match the ruler's training range; report the excluded count. Transcript IDs are RefSeq; gene symbols are parsed from Annotation/Divergence. Unmappable IDs stop this cohort, recorded, not substituted.

Per instrument × cohort: Spearman ρ vs donor age, permutation null (shuffle ages, hold scores; n_perm=200, seed 20260914), p, bootstrap CI (B=200, seed 20260918), n genes mapped out of the list's total, and the ID type of the matrix. Uncalibrated R² is reported and is never a gate.

Then the statistic that matters: per cohort, the paired bootstrap difference in ρ between the frozen ruler and each of the other two instruments — resample donors, recompute both ρ on the same resampled donors, report Δρ = ρ_ruler − ρ_other with its 95% percentile CI. A Δρ CI that excludes zero is the claim; two separate ρ values with overlapping CIs are not. Report Δρ even where it favours the other instrument.

**Pre-registered reading:**
- Δρ CI excludes zero in favour of the ruler (CI lower bound > 0) in ≥2 of 3 named cohorts → MD and/or their aging signature predict donor age materially worse than the ruler on independent fibroblast cohorts.
- Δρ CI includes zero (in the cohorts that distinguish the instruments under test) → the instruments are not distinguishable on age prediction at these n, and the Task 1 disagreement cannot be attributed to one being a weaker age instrument. Say which.
- Mixed (one cohort excludes zero in favour of the ruler, another in favour of the other instrument, or fewer than two ruler-favouring cohorts) → report per cohort, do not average, do not declare a winner.

An unscored cohort does not count as excluding zero. Do not report two separate ρ values as evidence of a difference without the paired Δρ CI. Do not declare a winner between the instruments.
"""

PREREG_TASK3 = """Written **before any MD3 MD co-variation ρ**, 2026-09-19. Report-only. No gate.

Across the GSE297234 Louvain cluster × timepoint rows with ≥20 cells (the FINDINGS_MD2.md unit, so the MD-vs-ruler cell is a check), report Spearman ρ between MD score and each of: the frozen ruler (reproduce the FINDINGS_MD2.md values as a check), the pluripotency−fibroblast score, the TGF-β score, and cell-state label treated as an ordinal along their trajectory.

Ordinal, frozen before this ρ: Fibroblast=0, PartialReprog=1, EarlyPluripotency=2, Pluripotency=3. NonReprog is not on that axis and is excluded from the ordinal ρ; n excluded is reported. This is our coding of their OSKM progression, not theirs.

Permutation null n_perm=200 seed 20260914 (shuffle MD, hold the other variable). Bootstrap CI B=200 seed 20260918. One table, both donors separate, never averaged. This says what MD co-varies with in this dataset without claiming what it “really” measures.
"""

MD2_NEGATIVE_HOLDS = (
    "No cluster at any resolution shows p ≤ 0.05. The per-cluster negative holds and is "
    "not a resolution artifact. The frozen ruler does not move in these cells at any "
    "granularity tested."
)
MD2_INSTRUMENTS_DISAGREE = (
    "MD moves per cluster, the frozen ruler does not. The instruments disagree on the "
    "same cells. What would settle it: an independent fibroblast age instrument scored "
    "on these same Louvain clusters, or the authors' Seurat/sctransform object with the "
    "frozen ruler projected onto it. No winner is declared."
)
MD2_PARTIAL_ZERO = (
    "Aged-donor Louvain clusters pairable at both endpoints (ok=True): n=1 (10:NonReprog). "
    "PartialReprog among them: n=0."
)
STAGE2_VERDICT_FIBRO = (
    "age score does not decline → the ruler reads a static donor property, not a "
    "modifiable state. Step 3 is not supported by this data."
)
FIBRO3_TASK3_AGED_D0D7 = (
    "Aged GM00731 d0→d7 (the all-cell T-B cell had p=+0.005): cluster 0 decline=-0.582 "
    "p=+0.652 n_ge=130/200 n_d0=4941 n_d7=2719; cluster 1 decline=+2.215 p=+0.144 "
    "n_ge=28/200 n_d0=76 n_d7=850; cluster 2 skipped (n_cells<20 at an endpoint)."
)

GENELIST_NULL_SEEDS = {
    "md_score": MD3_SEED,
    "tgfb_score": MD3_SEED + 101,
    "age_up": MD3_SEED + 202,
    "age_down": MD3_SEED + 303,
    "pluri_primary": MD3_SEED + 404,
}


def md3_log_banner(log, stage: str):
    log("=" * 100)
    log(f"MD3 {stage}  seed={MD3_SEED}  boot_seed={MD3_BOOT}  "
        f"n_perm={N_PERM}  n_boot={N_BOOT}  n_random={N_RANDOM_DIR}")
    log("=" * 100)
    log("FALSIFICATION.md and existing FINDINGS*.md are not modified.")
    log(f"Frozen ruler: {FROZEN_RULER} exists={FROZEN_RULER.exists()}")
    log("Primary metric: Spearman ρ (n_perm=200); uncalibrated R² is never a gate.")
    log("Refit nothing. Do not open GSE325735. Do not promote d7 to the Stage 2 endpoint.")
    log("Fail loudly: no substitute columns, no inferred ages, no silent coerce.")
    log("Louvain / AddModuleScore inherited from md2; not re-tuned.")


def load_manifest():
    if MANIFEST_PATH.exists():
        return load_json(MANIFEST_PATH)
    return dict(
        seed=MD3_SEED, boot_seed=MD3_BOOT, files={}, failures=[],
        columns={}, unusable=[], unverified=[], gene_sets_read=[],
        id_types={}, status="INIT",
    )


def save_manifest(man):
    dump_json(MANIFEST_PATH, jsonable(man))
    return MANIFEST_PATH


def _clip_msg(message, n=400):
    s = str(message)
    return s if len(s) <= n else s[:n] + "…"


def record_failure(step, message, details=None):
    man = load_manifest()
    rec = dict(step=step, message=_clip_msg(message), details=jsonable(details or {}))
    fails = man.setdefault("failures", [])
    if not any(f.get("step") == step and f.get("message") == rec["message"] for f in fails):
        fails.append(rec)
    save_manifest(man)
    return rec


def log_columns(man_key, columns, source):
    man = load_manifest()
    man.setdefault("columns", {})[man_key] = dict(
        source=str(source), columns=[str(c) for c in list(columns)],
    )
    save_manifest(man)


def log_geneset(man_key, name, genes, source, id_type=None):
    man = load_manifest()
    man.setdefault("gene_sets_read", [])
    rec = dict(
        key=man_key, name=str(name), n=int(len(genes)),
        source=str(source), genes=[str(g) for g in list(genes)],
        id_type=id_type,
    )
    existing = man["gene_sets_read"]
    if not any(x.get("key") == man_key and x.get("name") == rec["name"] for x in existing):
        existing.append(rec)
    save_manifest(man)


def log_id_type(man_key, rec):
    man = load_manifest()
    man.setdefault("id_types", {})[man_key] = jsonable(rec)
    save_manifest(man)


def write_prereg_flag(path: Path, text: str):
    if path.exists():
        existing = path.read_text(encoding="utf-8")
        if existing.strip() != text.strip():
            raise StopStep(
                "prereg",
                f"{path.name} already exists and does not match the frozen block. Not rewriting.",
            )
        return path
    path.write_text(text, encoding="utf-8")
    return path


def write_all_prereg_flags():
    write_prereg_flag(PREREG_TASK1_FLAG, PREREG_TASK1)
    write_prereg_flag(PREREG_TASK2_FLAG, PREREG_TASK2)
    write_prereg_flag(PREREG_TASK3_FLAG, PREREG_TASK3)
    return PREREG_TASK1_FLAG, PREREG_TASK2_FLAG, PREREG_TASK3_FLAG


def _fmt_num(x, d=3):
    if x is None:
        return "NA"
    if isinstance(x, (bool, np.bool_)):
        return str(bool(x))
    try:
        v = float(x)
    except (TypeError, ValueError):
        return str(x)
    if not np.isfinite(v):
        return "NA"
    return f"{v:+.{d}f}"


def _read_csv_if(path: Path):
    if path.exists():
        return pd.read_csv(path)
    return None


def _disk_row_lines(title, path, max_cols=16, max_rows=60):
    df = _read_csv_if(path)
    out = [f"### {title} (`{path.name}`)"]
    if df is None or not len(df):
        out.append("- missing or empty")
        out.append("")
        return out
    for i, (_, r) in enumerate(df.iterrows()):
        if i >= max_rows:
            out.append(f"- … {len(df) - max_rows} more rows")
            break
        bits = [f"{c}={r[c]}" for c in list(df.columns)[:max_cols]]
        out.append("- " + " ".join(bits))
    out.append("")
    return out


def progress_snapshot(next_action, stop=None, extra=""):
    man = load_manifest() if MANIFEST_PATH.exists() else {}
    stop_s = stop if stop else man.get("status", "running")
    lines = [
        "# PROGRESS_MD3",
        "",
        f"**STOP status:** {stop_s}",
        f"**Next action:** {next_action}",
        "",
        "## Seeds / gates",
        "",
        f"- seed `{MD3_SEED}`  boot `{MD3_BOOT}`  n_perm={N_PERM}  n_boot={N_BOOT}  n_random={N_RANDOM_DIR}",
        "- primary metric: Spearman ρ; uncalibrated R² is never a gate",
        f"- frozen ruler: `{FROZEN_RULER}` exists={FROZEN_RULER.exists()}  **not refit**",
        f"- MIN_CELLS_20={MIN_CELLS_20} (primary per-cluster)  MIN_CELLS_10={MIN_CELLS_10} (sensitivity only)  pooled n reported, no extra cutoff",
        f"- Louvain resolution={LOUVAIN_RES} inherited from md2; AddModuleScore nbin={AMS_NBIN} ctrl={AMS_CTRL}",
        f"- PREREG_TASK1.flag exists={PREREG_TASK1_FLAG.exists()}",
        f"- PREREG_TASK2.flag exists={PREREG_TASK2_FLAG.exists()}",
        f"- PREREG_TASK3.flag exists={PREREG_TASK3_FLAG.exists()}",
        f"- DECLARED_BEFORE_SCORES.flag exists={DECLARED_BEFORE_SCORES_FLAG.exists()}",
        "- GSE325735: out of scope",
        "- d7 is not the Stage 2 endpoint",
        "- Refit nothing. Do not re-tune Louvain or AddModuleScore.",
        "",
        "## Frozen ruler path",
        "",
        f"- `{FROZEN_RULER}` exists={FROZEN_RULER.exists()}",
        "",
        "## Clustering and scoring settings inherited from md2 (not re-tuned)",
        "",
        CLUSTER_SETTINGS,
        "",
        "## Pre-registered blocks (verbatim)",
        "",
        "### Task 1",
        "",
        PREREG_TASK1,
        "",
        "### Task 2",
        "",
        PREREG_TASK2,
        "",
        "### Task 3",
        "",
        PREREG_TASK3,
        "",
        "## Table S3 reasoning",
        "",
        TABLE_S3_REASONING,
        "",
        "## Gene-list versions and mapped counts (from disk; never from memory)",
        "",
    ]
    gs_path = MD3_DIR / "genesets_summary.json"
    if gs_path.exists():
        rec = load_json(gs_path)
        lines.append(f"- {json.dumps(jsonable(rec), ensure_ascii=False)[:8000]}")
        lines.append("")
    else:
        lines.append("- genesets_summary.json not yet written")
        lines.append("")
    idt = man.get("id_types") or {}
    if idt:
        lines += ["## ID type per matrix (from disk)", ""]
        for k, v in idt.items():
            lines.append(f"- **{k}** {json.dumps(jsonable(v), ensure_ascii=False)[:1500]}")
        lines.append("")
    lines += ["## Finished cells (from disk)", ""]
    for name, path in (
        ("Task 1 occupancy", MD3_DIR / "t1_occupancy.csv"),
        ("Task 1 pooled scores", MD3_DIR / "t1_pooled_scores.csv"),
        ("Task 1 cluster ≥20 scores", MD3_DIR / "t1_cluster20_scores.csv"),
        ("Task 1 cluster ≥10 scores", MD3_DIR / "t1_cluster10_scores.csv"),
        ("Task 1 pooled null", MD3_DIR / "t1_pooled_null.csv"),
        ("Task 1 cluster ≥20 null", MD3_DIR / "t1_cluster20_null.csv"),
        ("Task 1 cluster ≥10 null", MD3_DIR / "t1_cluster10_null.csv"),
        ("Task 1 reading", MD3_DIR / "t1_reading.json"),
        ("Task 2 instruments", MD3_DIR / "t2_instruments.csv"),
        ("Task 2 paired Δρ", MD3_DIR / "t2_delta_rho.csv"),
        ("Task 2 reading", MD3_DIR / "t2_reading.json"),
        ("Task 3 MD covariation", MD3_DIR / "t3_md_covariation.csv"),
    ):
        if str(path).endswith(".json"):
            if path.exists():
                rec = load_json(path)
                lines.append(f"### {name} (`{path.name}`)")
                lines.append(f"- {json.dumps(jsonable(rec), ensure_ascii=False)[:4000]}")
                lines.append("")
            else:
                lines.append(f"### {name} (`{path.name}`)")
                lines.append("- missing")
                lines.append("")
        else:
            lines += _disk_row_lines(name, path)
    for js_name in (
        "genesets_summary.json", "t1_summary.json", "t2_summary.json", "t3_summary.json",
    ):
        p = MD3_DIR / js_name
        if p.exists():
            rec = load_json(p)
            lines.append(f"### {js_name}")
            lines.append(f"- {json.dumps(jsonable(rec), ensure_ascii=False)[:4000]}")
            lines.append("")
    cols = man.get("columns") or {}
    if cols:
        lines += ["## Columns actually read", ""]
        for k, v in cols.items():
            lines.append(f"- **{k}** source=`{v.get('source')}` columns={v.get('columns')}")
        lines.append("")
    gsets = man.get("gene_sets_read") or []
    if gsets:
        lines += ["## Gene-set names actually read", ""]
        for g in gsets:
            lines.append(
                f"- **{g.get('key')}** name=`{g.get('name')}` n={g.get('n')} "
                f"id_type={g.get('id_type')} source=`{g.get('source')}`"
            )
        lines.append("")
    if extra:
        lines += ["## Note", "", extra, ""]
    if man.get("failures"):
        lines += ["## Failures (manifest)", ""]
        for f in man["failures"]:
            lines.append(f"- **{f.get('step')}:** {_clip_msg(f.get('message'), 240)}")
        lines.append("")
    lines += [
        "## Files",
        "",
        "- `src/md3_common.py`, `md3_idtype.py`, `md3_genesets.py`, `md3_task1.py`, "
        "`md3_task2.py`, `md3_task3.py`, `md3_findings.py`, `md3_run.py`",
        "- `results/md3/`",
        "- `FINDINGS_MD3.md`",
        "- `PROGRESS_MD3.md`",
        "",
    ]
    PROGRESS_PATH.write_text("\n".join(lines) + "\n", encoding="utf-8")
    return PROGRESS_PATH
