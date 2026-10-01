"""Fibroblast d7–d10 trajectory (Prompt G). Namespaced under results/fibro3/.

Does not modify any existing FINDINGS_*.md or FALSIFICATION.md. Does not open GSE325735.
Seed 20260914. Bootstrap seed 20260918. Spearman ρ is the reported statistic.
Uncalibrated R² is reported and is never a gate. Refit nothing.
d7 is not the primary Stage 2 endpoint.
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
    FROZEN_RULER, FIBRO_DIR, FIBRO_SEED, BOOT_SEED, N_PERM, N_BOOT,
    N_RANDOM_DIR, jsonable, fmt_ci,
)
from fibro2_common import (  # noqa: E402
    PCA_MAHAL_K, CXG_FIBRO_ATLAS_ID, pass_rho_bar, load_frozen_ruler,
    rho_null_donor, bootstrap_rho_ci,
)

FIBRO3_DIR = RESULTS / "fibro3"
FIBRO3_FIG = FIBRO3_DIR / "figures"
FIBRO3_PROC = DATA_PROC / "fibro3"
for _p in (FIBRO3_DIR, FIBRO3_FIG, FIBRO3_PROC):
    _p.mkdir(parents=True, exist_ok=True)

FIBRO3_SEED = int(SEED)  # 20260914
FIBRO3_BOOT = int(BOOT_SEED)  # 20260918
PREREG_DATE = "2026-09-18"
MIN_CELLS = 20
CLUSTER_K = 3
CLUSTER_NPC = 20
CLUSTER_BATCH = 2048
CLUSTER_N_INIT = 10
TRANSFER_NULL_BAR = 0.05

MANIFEST_PATH = FIBRO3_DIR / "manifest.json"
PROGRESS_PATH = ROOT / "PROGRESS_FIBRO3.md"
FINDINGS_PATH = ROOT / "FINDINGS_FIBRO3.md"
PREREG_TASK1_FLAG = FIBRO3_DIR / "PREREG_TASK1.flag"
PREREG_TASK4_FLAG = FIBRO3_DIR / "PREREG_TASK4.flag"
CLUSTER_SETTINGS_FLAG = FIBRO3_DIR / "CLUSTER_SETTINGS.flag"
DECLARED_BEFORE_SCORES_FLAG = FIBRO3_DIR / "DECLARED_BEFORE_SCORES.flag"

AGED_LINE = "GM00731"
YOUNG_LINE = "GM23815"
DAYS = (0, 3, 7, 10)
CXG_TAG = "cxg_a19d1667"
CXG_SCORES = FIBRO_DIR.parent / "fibro2" / f"ta_{CXG_TAG}_scores.csv"
CXG_RESULT = FIBRO_DIR.parent / "fibro2" / f"ta_{CXG_TAG}_result.json"
CXG_PROJ = DATA_PROC / "fibro2" / f"ta_{CXG_TAG}_proj.npz"
CXG_H5AD = DATA_RAW / "fibro2" / "cxg" / "8ae4471a-2e70-47b5-bc42-4f0df26a9996.h5ad"
CXG_DATASET_ID = CXG_FIBRO_ATLAS_ID

STANDING_CONSTRAINT = (
    "d7 was not the pre-registered Stage 2 endpoint. Nothing in this prompt may "
    "promote d7 to the primary result. The pre-registered d0→d10 verdict stands as "
    "written in `FINDINGS_FIBRO.md`. This prompt characterizes the trajectory; it "
    "does not re-decide it."
)

STAGE2_VERDICT_SENTENCE = (
    "age score does not decline → the ruler reads a static donor property, not a "
    "modifiable state. Step 3 is not supported by this data."
)
FIBRO2_TA_SENTENCE = (
    "Stage 2 no_decline must not be reported as a finding about reprogramming "
    "until the primary external cell is scored and passes."
)
FIBRO2_TB_SENTENCE = (
    "d0→d7 gives p ≤ 0.05 with no technical difference at d10 → the discrepancy "
    "is unexplained."
)

CLUSTER_SETTINGS = """Written **before any joint-cluster age score**, 2026-09-18. Not swept after numbers exist.

JOINT CLUSTERING (code-listed; chosen before looking at any age score):
- Scope: all cells of one donor, all timepoints concatenated. Donors clustered separately. Never pooled across donors.
- Genes: frozen-ruler overlap genes with count-sum > 0 (same gene filter as Stage 2 cluster spec in FINDINGS_FIBRO.md).
- Transform: log1p(CP10k) on those genes. Library size = row sum of the kept overlap genes (Stage 2).
- Dimensionality: sklearn.decomposition.PCA(n_components=20, or n_cells−1 if smaller, svd_solver="randomized", random_state=20260914) on float32 log1p(CP10k). Cell order: libraries concatenated in day order 0, 3, 7, 10.
- Cluster: sklearn.cluster.KMeans(n_clusters=3, random_state=20260914, n_init=10).
- k=3 is the Stage 2 frozen k, applied jointly so cluster IDs are comparable across days of the same donor. Resolution is k (not a Leiden resolution). Not swept.
- Seed `20260914`.
- MIN_CELLS = 20 (prompt minimum) to score age / pluri / extrapolation / random-direction null at a cluster×timepoint. QC (n, fraction, median UMI, median genes, mito fraction) is reported at all n, including n<20.

Amendment, still before any age score: the first listed dimensionality was IncrementalPCA(n_components=20, batch_size=2048). A run reached IncrementalPCA partial_fit on 2048×20354 and did not finish in usable time (~3 min/batch). No cluster labels and no age scores were written. Dimensionality was replaced by randomized PCA on the same matrix, same n_components, same KMeans k. This is not a k/resolution sweep.
"""

PREREG_TASK1 = """Written **before any joint-cluster age score**, 2026-09-18. No threshold in this block is re-tuned after numbers exist.

Using the GSE297234 Cell Ranger `filtered_feature_bc_matrix.h5` matrices already on disk. Cluster cells jointly across all timepoints within each donor (not within timepoint). Clustering settings are in `results/fibro3/CLUSTER_SETTINGS.flag` and are listed before any age score.

Per donor × cluster × timepoint report: n cells, fraction of that timepoint's cells, median UMI, median genes, mito fraction, frozen age score, pluripotency−fibroblast score (frozen lists from FINDINGS_FIBRO.md, with and without OSKM-family genes). Frozen ruler `results/fibro/frozen_ruler_ridge_raw.npz` used as-is. Refit nothing. Missing overlap genes left at z=0.

MIN_CELLS = 20. Do not read a d7–d10 trend off a cluster×timepoint with n < 20.

Operational definitions (frozen before any joint-cluster age score):
- C_ok(donor) = clusters with n ≥ 20 at both d7 and d10.
- Cluster reverse = age_score(d10) > age_score(d7) (score increased; older on the ruler).
- Cluster stable_or_decline = age_score(d10) ≤ age_score(d7).
- High-age cluster at d7 = cluster with n ≥ 20 at d7 and age_score(d7) equal to the maximum among those.
- Low-age cluster at d7 = cluster with n ≥ 20 at d7 and age_score(d7) equal to the minimum among those.
- Fraction = n_cells(cluster, day) / n_cells(donor, day). A missing cluster has fraction 0.
- Mix shift (reversal-direction) = frac(high-age, d10) > frac(high-age, d7) OR frac(low-age, d10) < frac(low-age, d7).
- Primary reading is the aged donor GM00731. GM23815 is tabulated as a contrast, never pooled, never the reading.

**Pre-registered reading** (only the outcome that fired; d7 is not the Stage 2 endpoint):
- **Compositional.** Individual clusters' age scores are stable or keep declining from d7 to d10, but the *mix* shifts — a high-age-score cluster grows, or a low-scoring one shrinks or drops out. Then the all-cell reversal is a composition artifact, the same failure mode as the blood result in `FINDINGS_BRAIN_PHASE1.md` §2. Say so in those words.
- **Cell-intrinsic.** The same cells' clusters individually reverse from d7 to d10. Then reprogramming genuinely pushes cells back along the age axis after d7, and that is a biological claim worth its own test.
- **Mixed.** Report which clusters do which. Do not average them into one number.
- **Unresolvable.** Too few cells per cluster × timepoint (threshold 20 cells) to read either way. Say so and stop rather than reading a trend off thin cells.

Firing rule (aged donor only):
- Unresolvable if C_ok is empty.
- Cell-intrinsic if C_ok is non-empty AND every cluster in C_ok reverses.
- Compositional if C_ok is non-empty AND every cluster in C_ok is stable_or_decline AND mix shift is True.
- Mixed if C_ok is non-empty AND at least one cluster reverses AND at least one does not.
- If all C_ok are stable_or_decline and mix shift is False: key=`no_listed_bullet`; report the tables; do not invent a fifth story.

Clusters are not lineage-tracked. "The same cells" is an inference from joint cluster identity, not from barcodes.
"""

PREREG_TASK4 = """Written **before any Task 4 ρ**, 2026-09-18. No threshold in this block is re-tuned after numbers exist. Option (a) is not run.

Choice: **(b)**. The primary external cell is superseded from GSE113957 to CELLxGENE dataset `a19d1667-a7b5-4556-9e5f-f9bfa690c0f1` (Human dermal fibroblast atlas; 10x-restricted adult donors).

Ground: FINDINGS_FIBRO2.md T-A pre-registration already stated that a 10x cohort is worth more than a larger bulk one because the platform shift is the thing under test. Stage 2 projects a GTEx bulk polyA (RNASeQCv2.4.2) frozen ruler onto 10x 3' scRNA-seq pseudobulks. GSE113957 is bulk FPKM; an FPKM rank-transform would be a different transform from the GTEx TMM/log2-CPM spec and would remain bulk-to-bulk, so it would not test the platform shift. The CELLxGENE 10x cohort tests that shift directly.

Take the Stage 1 ridge direction exactly as frozen in `results/fibro/frozen_ruler_ridge_raw.npz`. Refit nothing. Recalibrate nothing. Rescale nothing. Score Spearman ρ vs donor age on the superseded primary.

Freeze-and-project path (FINDINGS_EXTERNAL.md / T-A): TMM/edgeR log2-CPM (`prior.count=2`); genes restricted to the overlap with the frozen ruler; missing overlap genes left at z=0; z uses the frozen train mean and sd. Adult restriction ≥18. Ages from the record column actually used. Do not infer a donor age that the record does not state.

Metric: Spearman ρ vs stated age. Donor-level permutation null n_perm=200 (scores held fixed; shuffle donor ages; fresh Generator on seed `20260914`). Donor-bootstrap 95% percentile CI, B=200, seed `20260918`. Calibrated R² and uncalibrated R² reported; uncalibrated R² is never a gate. Pearson r reported, not a gate.

**Pre-registered reading:** if the superseded primary passes (ρ > 0, null ≤ 0.05), the frozen fibroblast ruler is a validated instrument outside GTEx and the Stage 2 trajectory is interpretable as a statement about the cells. If it fails, Stage 2 remains uninterpretable and every number in Tasks 1–3 is descriptive only.

GSE113957 is not scored in this prompt. Option (a) is not run. GSE226189 is not the superseded primary. GSE325735 is out of scope. Do not report both (a) and (b) and pick the better number.
"""

PREREG_TASK2 = """Written **before any Task 2 distance or correlation**, 2026-09-18.

Report-only, no gate. Per donor × cluster × timepoint with n ≥ 20, report extrapolation distance to the GTEx fibroblast training distribution: PCA k=50 nearest-neighbour Euclidean and Mahalanobis, same measures as FINDINGS_FIBRO2.md. Missing overlap genes at z=0. Frozen ruler used as-is.

Then Spearman ρ, across all donor × cluster × timepoint rows with n ≥ 20, between each distance and frozen age score. Donor-level permutation is not the unit here (the unit is a cluster×timepoint pseudobulk). Permutation null: shuffle age scores, n_perm=200, seed `20260914`. Bootstrap 95% percentile CI on the rows, B=200, seed `20260918`. Per-donor ρ is reported and is not averaged. Uncalibrated R² reported, never a gate.

If age score rises with distance from the training distribution, the ruler is being read outside its calibrated region and the d10 value is an extrapolation artifact rather than a measurement. State this as an observation with its correlation and CI. It is not a gate and does not overturn any verdict.
"""

PREREG_TASK3 = """Written **before any per-cluster random-direction p**, 2026-09-18.

For each donor × cluster with ≥20 cells at both endpoints, rerun the pre-registered random-direction null (200 permuted-weight directions, seed `20260914`) on d0→d7 and d0→d10 separately. Endpoint decline = age_score(day 0) − age_score(day T). Positive = younger on the ruler. Empirical p = (n_random ≥ real + 1) / (n_random + 1). Joint cluster IDs pair the same cluster across days (not size-rank). Do not average clusters or donors. d7 is not the Stage 2 endpoint.
"""


def fibro3_log_banner(log, stage: str):
    log("=" * 100)
    log(f"FIBRO3 {stage}  seed={FIBRO3_SEED}  boot_seed={FIBRO3_BOOT}  "
        f"n_perm={N_PERM}  n_boot={N_BOOT}  n_random={N_RANDOM_DIR}")
    log("=" * 100)
    log("FALSIFICATION.md and existing FINDINGS*.md are not modified.")
    log(f"Frozen ruler: {FROZEN_RULER} exists={FROZEN_RULER.exists()}")
    log("Primary metric: Spearman ρ (n_perm=200); uncalibrated R² is never a gate.")
    log("Refit nothing. Do not open GSE325735. Do not promote d7 to the Stage 2 endpoint.")
    log(STANDING_CONSTRAINT)
    log("Fail loudly: no substitute columns, no inferred ages, no silent coerce.")


def load_manifest():
    if MANIFEST_PATH.exists():
        return load_json(MANIFEST_PATH)
    return dict(
        seed=FIBRO3_SEED, boot_seed=FIBRO3_BOOT, files={}, failures=[],
        columns={}, unusable=[], unverified=[], status="INIT",
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
    man["status"] = "STOP"
    save_manifest(man)
    return rec


def log_columns(man_key, columns, source):
    man = load_manifest()
    man.setdefault("columns", {})[man_key] = dict(
        source=str(source), columns=[str(c) for c in list(columns)],
    )
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
    write_prereg_flag(PREREG_TASK4_FLAG, PREREG_TASK4)
    write_prereg_flag(PREREG_TASK1_FLAG, PREREG_TASK1)
    write_prereg_flag(CLUSTER_SETTINGS_FLAG, CLUSTER_SETTINGS)
    write_prereg_flag(FIBRO3_DIR / "PREREG_TASK2.flag", PREREG_TASK2)
    write_prereg_flag(FIBRO3_DIR / "PREREG_TASK3.flag", PREREG_TASK3)
    return PREREG_TASK4_FLAG, PREREG_TASK1_FLAG, CLUSTER_SETTINGS_FLAG


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


def _disk_row_lines(title, path, max_cols=18):
    df = _read_csv_if(path)
    out = [f"### {title} (`{path.name}`)"]
    if df is None or not len(df):
        out.append("- missing or empty")
        out.append("")
        return out
    for _, r in df.iterrows():
        bits = [f"{c}={r[c]}" for c in list(df.columns)[:max_cols]]
        out.append("- " + " ".join(bits))
    out.append("")
    return out


def progress_snapshot(next_action, stop=None, extra=""):
    man = load_manifest() if MANIFEST_PATH.exists() else {}
    stop_s = stop if stop else man.get("status", "running")
    lines = [
        "# PROGRESS_FIBRO3",
        "",
        f"**STOP status:** {stop_s}",
        f"**Next action:** {next_action}",
        "",
        "## Seeds / gates",
        "",
        f"- seed `{FIBRO3_SEED}`  boot `{FIBRO3_BOOT}`  n_perm={N_PERM}  n_boot={N_BOOT}  n_random={N_RANDOM_DIR}",
        f"- primary metric: Spearman ρ; uncalibrated R² is never a gate",
        f"- frozen ruler: `{FROZEN_RULER}` exists={FROZEN_RULER.exists()}",
        f"- MIN_CELLS={MIN_CELLS}  CLUSTER_K={CLUSTER_K}  CLUSTER_NPC={CLUSTER_NPC}  "
        f"CLUSTER_BATCH={CLUSTER_BATCH}  n_init={CLUSTER_N_INIT}",
        f"- PREREG_TASK4.flag exists={PREREG_TASK4_FLAG.exists()}",
        f"- PREREG_TASK1.flag exists={PREREG_TASK1_FLAG.exists()}",
        f"- CLUSTER_SETTINGS.flag exists={CLUSTER_SETTINGS_FLAG.exists()}",
        f"- DECLARED_BEFORE_SCORES.flag exists={DECLARED_BEFORE_SCORES_FLAG.exists()}",
        f"- Task 4 choice: (b) CELLxGENE `{CXG_DATASET_ID}` 10x; option (a) not run",
        f"- GSE325735: out of scope",
        f"- d7 is not the Stage 2 endpoint",
        "",
        "## Standing constraint",
        "",
        STANDING_CONSTRAINT,
        "",
        "## Frozen ruler path",
        "",
        f"- `{FROZEN_RULER}` exists={FROZEN_RULER.exists()}",
        "",
        "## Clustering settings (verbatim, frozen before any age score)",
        "",
        CLUSTER_SETTINGS,
        "",
        "## Pre-registered blocks (verbatim)",
        "",
        "### Task 4 supersession",
        "",
        PREREG_TASK4,
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
        "## Finished cells (from disk)",
        "",
    ]
    for name, path in (
        ("Task 4 result", FIBRO3_DIR / "t4_result.csv"),
        ("Task 1 cluster table", FIBRO3_DIR / "t1_cluster_table.csv"),
        ("Task 1 composition d7–d10", FIBRO3_DIR / "t1_composition_d7_d10.csv"),
        ("Task 2 extrap", FIBRO3_DIR / "t2_extrap.csv"),
        ("Task 2 correlation", FIBRO3_DIR / "t2_correlation.csv"),
        ("Task 3 null", FIBRO3_DIR / "t3_null.csv"),
    ):
        lines += _disk_row_lines(name, path)
    for js_name in (
        "t4_summary.json", "t1_summary.json", "t1_reading.json",
        "t2_summary.json", "t3_summary.json",
    ):
        p = FIBRO3_DIR / js_name
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
    if extra:
        lines += ["## Note", "", extra, ""]
    if man.get("unusable"):
        lines += ["## unusable", ""]
        for u in man["unusable"]:
            lines.append(f"- **{u.get('accession', u.get('step'))}:** {u.get('reason', u.get('message'))}")
        lines.append("")
    if man.get("failures"):
        lines += ["## Failures (manifest)", ""]
        for f in man["failures"]:
            lines.append(f"- **{f.get('step')}:** {_clip_msg(f.get('message'), 240)}")
        lines.append("")
    lines += [
        "## Files",
        "",
        "- `src/fibro3_common.py`, `fibro3_task4.py`, `fibro3_task1.py`, `fibro3_task2.py`, "
        "`fibro3_task3.py`, `fibro3_findings.py`, `fibro3_run.py`",
        "- `results/fibro3/`",
        "- `FINDINGS_FIBRO3.md`",
        "- `PROGRESS_FIBRO3.md`",
        "",
    ]
    PROGRESS_PATH.write_text("\n".join(lines) + "\n", encoding="utf-8")
    return PROGRESS_PATH
