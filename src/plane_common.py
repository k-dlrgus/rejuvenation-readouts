"""Prompt M step 2 — age–identity plane: does the curve bend?

Namespaced under results/plane/. Does not modify any existing FINDINGS_*.md or FALSIFICATION.md.
Does not open GSE325735. Refit nothing. Frozen ruler used as-is.
Seed 20260914. Bootstrap seed 20260918. n_perm=200. n_boot=200. n_random=200.
Reuse src/md2_*.py Louvain labels + AddModuleScore exactly; do not re-cluster, re-label, or re-tune.
Reuse md4 cell_scores for observed frozen-ruler / Table S3 / pluripotency−fibroblast vectors.
"""
from __future__ import annotations

import json
import re
import sys
from pathlib import Path

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parent))
from config import RESULTS, DATA_PROC, ROOT, SEED  # noqa: E402
from brain_phase1_common import Logger, dump_json, load_json  # noqa: E402
from gtex_common import StopStep  # noqa: E402
from fibro_common import (  # noqa: E402
    FROZEN_RULER, BOOT_SEED, N_PERM, N_BOOT, N_RANDOM_DIR, jsonable,
    PLURI_ENDOGENOUS, FIBRO_IDENTITY, OSKM_FAMILY,
)
from fibro3_common import AGED_LINE, YOUNG_LINE, DAYS  # noqa: E402
from fibro2_common import load_frozen_ruler, PCA_MAHAL_K  # noqa: E402
from md2_common import (  # noqa: E402
    MD2_DIR, MD2_PROC, MMC3_XLSX, AMS_NBIN, AMS_CTRL, SEURAT_LOGNORMALIZE_SCALE,
    LOUVAIN_RES, CLUSTER_SETTINGS,
)
from md3_common import MD3_DIR, GENELIST_NULL_SEEDS  # noqa: E402
from md4_common import (  # noqa: E402
    MD4_DIR, MD4_PROC, STATE_NAMES, DONOR_INDEX,
    MD2_POOLED_MD_PR, MD2_POOLED_MD_NR, MD2_N_PR, MD2_N_NR,
    MD3_TASK1_STATUS, MD3_TASK1_READING, MD3_TASK2_READING,
    MD2_NEGATIVE_HOLDS, STAGE2_VERDICT_FIBRO,
)
from md5_common import (  # noqa: E402
    MD5_DIR,
    MD4_TASK1_YOUNGER_CLAIM, MD4_TASK1_STATUS_FIRED,
)

PLANE_DIR = RESULTS / "plane"
PLANE_FIG = PLANE_DIR / "figures"
PLANE_PROC = DATA_PROC / "plane"
for _p in (PLANE_DIR, PLANE_FIG, PLANE_PROC):
    _p.mkdir(parents=True, exist_ok=True)

PLANE_SEED = int(SEED)  # 20260914
PLANE_BOOT = int(BOOT_SEED)  # 20260918
PREREG_DATE = "2026-09-21"

MANIFEST_PATH = PLANE_DIR / "manifest.json"
PROGRESS_PATH = ROOT / "PROGRESS_PLANE.md"
FINDINGS_PATH = ROOT / "FINDINGS_PLANE.md"
PREREG_TASK1_FLAG = PLANE_DIR / "PREREG_TASK1.flag"
PREREG_TASK2_FLAG = PLANE_DIR / "PREREG_TASK2.flag"
PREREG_TASK3_FLAG = PLANE_DIR / "PREREG_TASK3.flag"
DECLARED_BEFORE_SCORES_FLAG = PLANE_DIR / "DECLARED_BEFORE_SCORES.flag"

MIN_CELLS = 30
TRAJECTORY = ("Fibroblast", "PartialReprog", "EarlyPluripotency", "Pluripotency")
OFF_TRAJECTORY = "NonReprog"
SEGMENT_A = ("Fibroblast", "PartialReprog")
SEGMENT_B_PREFERRED_END = "Pluripotency"
SEGMENT_B_FALLBACK_END = "EarlyPluripotency"
IDENTITY_KEY = "pluri_primary"
AGE_INSTRUMENTS = ("frozen_ruler", "md_score", "age_up_minus_age_down")
INST_INDEX = {k: i for i, k in enumerate(AGE_INSTRUMENTS)}
INST_LABELS = {
    "frozen_ruler": "frozen ruler",
    "md_score": "MD",
    "age_up_minus_age_down": "age-up − age-down",
}
DX_EPS = 1e-12
AMS_CTRL_STREAM_OFFSET = 17
RULER_NULL_CHUNK = 512

FIBRO3_EXTRAP_RHO_NN = -0.625
FIBRO3_EXTRAP_CI = (-0.834, -0.320)
FIBRO3_EXTRAP_N = 22
MD4_EXTRAP_RHO_NN = -0.280
MD3_MD_VS_PLURI_GM00731 = -0.938
MD3_MD_VS_PLURI_GM23815 = -0.960

MD5_FIBRO_STAGE2_QUOTE = (
    "age score does not decline → the ruler reads a static donor property, not a "
    "modifiable state. Step 3 is not supported by this data."
)
MD5_MD3_T1_QUOTE = (
    "Ruler does not decline in the pooled PartialReprog state while MD does, both with "
    "their own nulls. The instruments disagree on exactly the population the published "
    "claim concerns. No winner is declared."
)
MD5_MD3_T2_QUOTE = (
    "Δρ CI excludes zero in favour of the ruler in 2 of 3 named cohorts. MD and/or their "
    "aging signature predict donor age materially worse than the ruler on independent "
    "fibroblast cohorts."
)
MD5_MD4_T1_STATUS_QUOTE = (
    "Task 1: GM00731 d3: `both_lower`; d7: `both_lower`; GM23815 d10: `both_lower`. "
    "Fired: GM00731 `both_lower`; GM23815 `both_lower`."
)
MD5_MD4_T1_READING_QUOTE = (
    "Both instruments agree that partially reprogrammed cells read younger. Their claim "
    "is supported and our earlier negatives reflected the wrong test (a within-state "
    "trajectory with no baseline)."
)
MD5_FIRED_GM00731 = (
    "PartialReprog below Fibroblast on the frozen ruler, beating its null. Responding "
    "cells read younger than their own starting population. The rejuvenation reading "
    "holds and is not an artifact of NonReprog reading old."
)
MD5_STATUS_LINE = (
    "Task 1: GM00731 scheme=a d10: `pr_below_fibroblast`; GM23815 scheme=a d3: "
    "`pr_below_fibroblast`; d7: `pr_below_fibroblast`; d10: `pr_below_fibroblast`. "
    "Fired: GM00731 `pr_below_fibroblast` scheme=a; GM23815 `pr_below_fibroblast` scheme=a."
)

PREREG_TASK1 = """Written **before any PLANE state-median, slope, bend statistic, or null p-value**, 2026-09-21. No threshold in this block is re-tuned after numbers exist.

Using the src/md2_*.py Louvain clustering already computed (not re-run, not re-tuned). Cell-state labels are the md2 assignment: argmax of mean AddModuleScore of mmc3 Reprog_cell_state_signatures columns Fibroblast / PartialReprog / EarlyPluripotency / Pluripotency / NonReprog. Labels are read from results/md2/t2_cluster_labels_*.csv. Not re-labelled. AddModuleScore for MD is the per-cell vector already stored in data/processed/md2/louvain_ams_*.npz (the md2 computation). Table S3 age-up minus age-down and the frozen-ruler / pluripotency−fibroblast per-cell scores are the vectors already stored in data/processed/md4/cell_scores_*.npz (the md4 computation). Frozen ruler: results/fibro/frozen_ruler_ridge_raw.npz used as-is. Refit nothing. Do not re-cluster, re-label, re-tune, rescore observed instruments, or re-order the trajectory after seeing scores.

Trajectory order, fixed before scoring: Fibroblast → PartialReprog → EarlyPluri → Pluripotency. NonReprog is off-trajectory: score and plot it, but exclude it from every trajectory statistic.

Unit of analysis: the cell state, per donor. Pool timepoints within each state (the state, not the day, is the trajectory position here). Report the timepoint composition of every state. Require ≥30 cells per state per donor; a state below that is reported as not scored, not dropped silently.

Two coordinates per cell: x = identity = pluripotency−fibroblast score (frozen lists from FINDINGS_FIBRO.md; primary version as stored in md4 cell_scores, drop_oskm=False). y = age. Three age instruments, plotted and tested separately, never combined: frozen ruler; MD; Table S3 age-up minus age-down.

Per donor, per instrument:

1. State medians and bootstrap 95% CIs (B=200, cells resampled within state) on both coordinates.
2. Segment slopes. Δage / Δidentity for segment A (Fibroblast → PartialReprog) and segment B (PartialReprog → Pluripotency; if Pluripotency is not scored, use EarlyPluri and say so).
3. Bend statistic = slope_A − slope_B, with bootstrap CI.
4. Null: for the frozen ruler, the same statistic on 200 permuted-weight random directions (seed 20260914). For each gene-list score (MD, age-up−age-down), 200 size-matched, expression-bin-matched random gene sets. Report the empirical p for the observed bend being at least as large as the null's.

Young direction = destination state median age below origin (Δage < 0). If Δidentity > 0, that is a negative slope; if Δidentity < 0, a positive slope. Slope CI excludes 0 in the young direction iff the bootstrap CI of the slope lies entirely on that side of 0. Segment B is flat if its slope CI includes 0; reversed (old direction) if the slope CI excludes 0 on the age-rises side (Δage > 0).

Empirical p for the bend: p_more_negative = (n_null ≤ obs + 1) / (n_null + 1). That is the one-sided test that the observed bend is at least as large as the null's in the rejuvenation-separable direction (slope_A more negative than slope_B when identity increases). Also report p_more_positive. “Beats its null” = p_more_negative ≤ 0.05. Seed 20260914 for random directions and random gene sets. Bootstrap seed 20260918 with a documented per-(donor, instrument) offset. Donors never pooled. Instruments never combined. NonReprog never in a trajectory statistic.

Sanity, before any new median or bend: (1) reproduce FINDINGS_MD2.md pooled MD means on the aged donor (PartialReprog +0.150, NonReprog +0.439; t2_reproduction.json); (2) occupancy matches results/md4/t1_cell_counts.csv summed within state. If either fails, stop. Cell-state labels have drifted.

**Pre-registered reading (per donor, per instrument; only what fired):**
- **Bend:** age falls in segment A (slope_A CI excludes 0 in the young direction) **and** segment B's age change is flat or reversed (slope_B CI includes 0 or is in the old direction) **and** the bend statistic beats its null at p ≤ 0.05 → age moves before identity is lost; rejuvenation is separable from dedifferentiation on this instrument.
- **Lockstep:** age falls in both segments with overlapping slopes and the bend statistic does not beat its null → this instrument cannot separate rejuvenation from identity loss.
- **Mixed / other:** report exactly what fired.

Instrument comparison, report-only: state which instruments bend and which run in lockstep. FINDINGS_MD3.md Task 3 found MD correlates ρ ≈ −0.94 with the pluripotency−fibroblast score, so MD is expected to run in lockstep; say whether it does. If the frozen ruler bends and MD does not, that is the separation between the instruments, stated as observed, with no winner declared.
"""

PREREG_TASK2 = """Written **before any PLANE extrapolation distance**, 2026-09-21. Report-only. No gate.

Per state per donor, report the extrapolation distance to the GTEx fibroblast training distribution (same PCA k=50 nearest-neighbour Euclidean and Mahalanobis measures as FINDINGS_FIBRO2.md / src/fibro2_extrap.py). Sum Louvain-kept UMIs to a pseudobulk on the frozen-ruler gene space (md4 louvain_rulerY aligned). TMM among that donor's scored-state rows, then frozen μ/σ/w. Do not re-TMM among cells. Do not refit PCA. Do not rescore the per-cell instruments of Task 1.

State the direction extrapolation pushes the result, using the evidence already in the repo: FINDINGS_FIBRO3.md Task 2 (the same PCA k=50 NN / Mahalanobis measures as FINDINGS_FIBRO2.md) found the ruler score falls with distance from training (ρ = −0.625); FINDINGS_MD4.md found a weaker version (ρ = −0.280). Pluripotent cells sit furthest from fibroblast training. If extrapolation pushes readings down with distance, it manufactures lockstep, not a bend — so a bend observed in the frozen ruler would be conservative with respect to extrapolation, and lockstep would be uninterpretable. State whichever applies from the numbers, and say so plainly. Do not use this to rescue a lockstep result.
"""

PREREG_TASK3 = """Written **before any PLANE figure**, 2026-09-21.

One figure per donor, three panels (frozen ruler, MD, age-up−age-down), x = identity, y = age, the four trajectory states as points with 2D bootstrap CIs joined in order, NonReprog as a separate marker, the random-direction or random-set envelope shaded behind. Same identity-axis scaling within a donor. results/plane/figures/plane_<donor>.png.
"""


def plane_log_banner(log, stage: str):
    log("=" * 100)
    log(f"PLANE {stage}  seed={PLANE_SEED}  boot_seed={PLANE_BOOT}  "
        f"n_perm={N_PERM}  n_boot={N_BOOT}  n_random={N_RANDOM_DIR}")
    log("=" * 100)
    log("FALSIFICATION.md and existing FINDINGS*.md are not modified.")
    log(f"Frozen ruler: {FROZEN_RULER} exists={FROZEN_RULER.exists()}")
    log("Primary metric: Spearman ρ (n_perm=200); uncalibrated R² is never a gate.")
    log("Refit nothing. Do not open GSE325735. Do not promote d7 to the Stage 2 endpoint.")
    log("Fail loudly: no substitute columns, no inferred ages, no silent coerce.")
    log("Louvain / labels / AddModuleScore inherited from md2; not re-tuned.")
    log("md4 cell_scores reused for observed instruments; cells are not rescored.")
    log(f"MIN_CELLS={MIN_CELLS}. Trajectory={list(TRAJECTORY)}. "
        f"NonReprog off-trajectory. Donors never pooled. Instruments never combined.")


def load_manifest():
    if MANIFEST_PATH.exists():
        return load_json(MANIFEST_PATH)
    return dict(
        seed=PLANE_SEED, boot_seed=PLANE_BOOT, files={}, failures=[],
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


def record_unusable(step, message, details=None):
    man = load_manifest()
    rec = dict(step=step, message=_clip_msg(message), details=jsonable(details or {}))
    rows = man.setdefault("unusable", [])
    if not any(f.get("step") == step and f.get("message") == rec["message"] for f in rows):
        rows.append(rec)
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


_ENS_RE = re.compile(r"^ENS[A-Z]*[GT]\d+", re.I)
_NM_RE = re.compile(r"^(NM_|NR_|XM_|XR_)\d+", re.I)
_SYM_RE = re.compile(r"^[A-Z][A-Z0-9\-.]{1,18}$")


def detect_id_type(ids, log, tag, source):
    ids = [str(x) for x in np.asarray(ids).astype(str)]
    n = int(len(ids))
    if n == 0:
        rec = dict(tag=tag, source=str(source), n=0, kind="empty",
                   n_ensembl=0, n_refseq=0, n_symbol_like=0, examples=[])
        log_id_type(tag, rec)
        if log is not None:
            log(f"[id_type {tag}] EMPTY source={source}")
        return rec
    n_ens = int(sum(1 for x in ids if _ENS_RE.match(x.split(".")[0])))
    n_ref = int(sum(1 for x in ids if _NM_RE.match(x)))
    n_sym = int(sum(1 for x in ids if _SYM_RE.match(x.upper()) and not _ENS_RE.match(x) and not _NM_RE.match(x)))
    frac_ens = n_ens / n
    frac_ref = n_ref / n
    frac_sym = n_sym / n
    if frac_ens >= 0.8:
        kind = "ensembl"
    elif frac_ref >= 0.5:
        kind = "refseq_transcript"
    elif frac_sym >= 0.5:
        kind = "symbol"
    else:
        kind = "unknown"
    rec = dict(
        tag=tag, source=str(source), n=n, kind=kind,
        n_ensembl=n_ens, n_refseq=n_ref, n_symbol_like=n_sym,
        frac_ensembl=frac_ens, frac_refseq=frac_ref, frac_symbol_like=frac_sym,
        examples=ids[:12],
    )
    log_id_type(tag, rec)
    if log is not None:
        log(f"[id_type {tag}] kind={kind} n={n} n_ensembl={n_ens} n_refseq={n_ref} "
            f"n_symbol_like={n_sym} source={source} examples={ids[:8]}")
    return rec


def require_mappable(rec, allowed, log, step):
    if rec.get("kind") not in allowed:
        msg = (
            f"{rec.get('tag')}: ID type {rec.get('kind')!r} not in {allowed}. "
            f"n={rec.get('n')} n_ensembl={rec.get('n_ensembl')} n_refseq={rec.get('n_refseq')} "
            f"n_symbol_like={rec.get('n_symbol_like')} examples={rec.get('examples')}. "
            "Not mapping. This step stops."
        )
        if log is not None:
            log(f"[id_type STOP] {msg}")
        raise StopStep(step, msg, details=rec)
    return rec


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


def _slope(dy, dx):
    if not (np.isfinite(dy) and np.isfinite(dx)):
        return np.nan
    if abs(float(dx)) < DX_EPS:
        return np.nan
    return float(dy / dx)


def segment_B_end(qrow):
    end = (qrow or {}).get("segment_B_end")
    if end in (None, "NA") or (isinstance(end, float) and not np.isfinite(end)):
        return None, False
    return str(end), bool((qrow or {}).get("segment_B_end_is_fallback"))


def flag_true(v):
    if v is None:
        return False
    if isinstance(v, (float, np.floating)) and not np.isfinite(v):
        return False
    if isinstance(v, str) and v.strip().lower() in ("", "nan", "none", "nat", "na"):
        return False
    return bool(v)


def boot_seed(line, instrument):
    return int(PLANE_BOOT + 1000 * DONOR_INDEX[line] + 10 * INST_INDEX[instrument])


def ci_includes_0(lo, hi):
    if lo is None or hi is None:
        return None
    try:
        lo = float(lo)
        hi = float(hi)
    except (TypeError, ValueError):
        return None
    if not (np.isfinite(lo) and np.isfinite(hi)):
        return None
    return bool(lo <= 0.0 <= hi)


def cis_overlap(a_lo, a_hi, b_lo, b_hi):
    vals = [a_lo, a_hi, b_lo, b_hi]
    try:
        vals = [float(v) for v in vals]
    except (TypeError, ValueError):
        return None
    if not all(np.isfinite(v) for v in vals):
        return None
    return bool(not (vals[1] < vals[2] or vals[3] < vals[0]))


def slope_young_side(dx):
    """Sign of slope that corresponds to Δage < 0 given observed Δidentity."""
    if dx is None or (not np.isfinite(dx)) or abs(float(dx)) < DX_EPS:
        return None
    return float(-np.sign(dx))


def ci_excludes_0_on_side(lo, hi, side):
    """True iff CI is entirely on `side` of 0. side is -1 (negative) or +1 (positive)."""
    if side is None or lo is None or hi is None:
        return None
    try:
        lo = float(lo)
        hi = float(hi)
        side = float(side)
    except (TypeError, ValueError):
        return None
    if not (np.isfinite(lo) and np.isfinite(hi) and np.isfinite(side) and side != 0):
        return None
    if side < 0:
        return bool(hi < 0.0)
    return bool(lo > 0.0)


def _read_csv_if(path: Path):
    if path.exists():
        return pd.read_csv(path)
    return None


def _disk_row_lines(title, path, max_cols=18, max_rows=80):
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
        "# PROGRESS_PLANE",
        "",
        f"**STOP status:** {stop_s}",
        f"**Next action:** {next_action}",
        "",
        "## Seeds / gates",
        "",
        f"- seed `{PLANE_SEED}`  boot `{PLANE_BOOT}`  n_perm={N_PERM}  n_boot={N_BOOT}  n_random={N_RANDOM_DIR}",
        "- primary metric: Spearman ρ; uncalibrated R² is never a gate",
        f"- frozen ruler: `{FROZEN_RULER}` exists={FROZEN_RULER.exists()}  **not refit**",
        f"- MIN_CELLS={MIN_CELLS} (per state per donor; timepoints pooled within state; never pooled across donors)",
        f"- trajectory (fixed): {' → '.join(TRAJECTORY)}",
        f"- NonReprog off-trajectory; scored and plotted; excluded from every trajectory statistic",
        f"- age instruments (never combined): {list(AGE_INSTRUMENTS)}",
        f"- identity: `{IDENTITY_KEY}` = pluripotency−fibroblast (FINDINGS_FIBRO.md frozen lists; md4 primary, drop_oskm=False)",
        f"- Louvain resolution={LOUVAIN_RES} inherited from md2; AddModuleScore nbin={AMS_NBIN} ctrl={AMS_CTRL}",
        f"- PCA_MAHAL_K={PCA_MAHAL_K} (FINDINGS_FIBRO2.md measures; Task 2 state-pooled TMM, not the md4 state×timepoint rows)",
        f"- PREREG_TASK1.flag exists={PREREG_TASK1_FLAG.exists()}",
        f"- PREREG_TASK2.flag exists={PREREG_TASK2_FLAG.exists()}",
        f"- PREREG_TASK3.flag exists={PREREG_TASK3_FLAG.exists()}",
        f"- DECLARED_BEFORE_SCORES.flag exists={DECLARED_BEFORE_SCORES_FLAG.exists()}",
        "- GSE325735: out of scope",
        "- d7 is not the Stage 2 endpoint",
        "- Refit nothing. Do not re-cluster, re-label, or re-tune Louvain or AddModuleScore.",
        "- Do not rescore observed cells. md4 cell_scores / md2 louvain_ams reused.",
        "- Donors never pooled. Instruments never combined.",
        "",
        "## Frozen ruler path",
        "",
        f"- `{FROZEN_RULER}` exists={FROZEN_RULER.exists()}",
        "",
        "## Frozen gene lists (FINDINGS_FIBRO.md; not reconstructed)",
        "",
        f"- pluripotency endogenous: {list(PLURI_ENDOGENOUS)}",
        f"- fibroblast identity: {list(FIBRO_IDENTITY)}",
        f"- OSKM-family (dropped only for the without-OSKM series; md4 primary keeps POU5F1): {list(OSKM_FAMILY)}",
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
        "## Gene-list versions and mapped counts (from disk; never from memory)",
        "",
    ]
    gs_path = PLANE_DIR / "genesets_summary.json"
    if gs_path.exists():
        rec = load_json(gs_path)
        lines.append(f"- {json.dumps(jsonable(rec), ensure_ascii=False)[:8000]}")
        lines.append("")
    else:
        md4_gs = MD4_DIR / "genesets_summary.json"
        if md4_gs.exists():
            rec = load_json(md4_gs)
            lines.append(f"- inherited from md4 `{md4_gs}`: "
                         f"{json.dumps(jsonable(rec), ensure_ascii=False)[:8000]}")
            lines.append("")
        else:
            lines.append("- genesets_summary.json not yet written")
            lines.append("")
    for line in (AGED_LINE, YOUNG_LINE):
        p = MD4_DIR / f"t1_s3_ams_meta_{line}.json"
        if p.exists():
            rec = load_json(p)
            lines.append(f"- md4 S3 AMS meta `{p.name}`: {json.dumps(jsonable(rec), ensure_ascii=False)[:4000]}")
            lines.append("")
    idt = man.get("id_types") or {}
    if idt:
        lines += ["## ID type per matrix (from disk)", ""]
        for k, v in idt.items():
            lines.append(f"- **{k}** {json.dumps(jsonable(v), ensure_ascii=False)[:1500]}")
        lines.append("")
    lines += ["## Finished cells (from disk)", ""]
    for name, path in (
        ("Task 1 cell counts", PLANE_DIR / "t1_cell_counts.csv"),
        ("Task 1 timepoint composition", PLANE_DIR / "t1_timepoint_composition.csv"),
        ("Task 1 qualify", PLANE_DIR / "t1_qualify.csv"),
        ("Task 1 sanity", PLANE_DIR / "t1_sanity.json"),
        ("Task 1 state medians", PLANE_DIR / "t1_state_medians.csv"),
        ("Task 1 slopes", PLANE_DIR / "t1_slopes.csv"),
        ("Task 1 bend", PLANE_DIR / "t1_bend.csv"),
        ("Task 1 null summary", PLANE_DIR / "t1_null_summary.csv"),
        ("Task 1 reading", PLANE_DIR / "t1_reading.json"),
        ("Task 2 extrap", PLANE_DIR / "t2_extrap.csv"),
        ("Task 2 summary", PLANE_DIR / "t2_summary.json"),
        ("Task 3 figures", PLANE_DIR / "t3_figures.json"),
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
    if man.get("unusable"):
        lines += ["## Unusable / not computable (manifest)", ""]
        for f in man["unusable"]:
            lines.append(f"- **{f.get('step')}:** {_clip_msg(f.get('message'), 240)}")
        lines.append("")
    lines += [
        "## Files",
        "",
        "- `src/plane_common.py`, `plane_task1.py`, `plane_nulls.py`, "
        "`plane_task2.py`, `plane_task3.py`, `plane_findings.py`, `plane_run.py`",
        "- `results/plane/`",
        "- `FINDINGS_PLANE.md`",
        "- `PROGRESS_PLANE.md`",
        "",
    ]
    PROGRESS_PATH.write_text("\n".join(lines) + "\n", encoding="utf-8")
    return PROGRESS_PATH
