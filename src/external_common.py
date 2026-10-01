"""SEA-AD MTG external transfer of the frozen DLPFC age direction.

Does the DLPFC age direction (raw ridge and identity-residualized ridge, the
0.317→0.300 fit) predict chronological age in SEA-AD MTG, a cohort it was
never trained on? Namespaced under results/external/. Does not modify any
prior FINDINGS*.md or FALSIFICATION.md. Seed 20260914.

Fail loudly: no substitute columns, no invented cell-type synonyms, no silent
coerce. Dataset ID is the fetched CELLxGENE record, not invented.
"""
from __future__ import annotations

import json
import re
import sys
from pathlib import Path

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parent))
from config import RESULTS, DATA_RAW, DATA_PROC, ROOT, SEED  # noqa: E402
from brain_phase1_common import (  # noqa: E402
    Logger, dump_json, load_json, strip_ensembl, r2_mae, LOGCPM_OBS, HASH_POOL_LIMITATION,
)
from gtex_common import (  # noqa: E402
    StopStep, pearson_safe, spearman_safe, pred_scores, pass_r_bar, fmt,
)
from target_common import md_table, MIN_N_TYPE  # noqa: E402
from search_cxg import API as CXG_API  # noqa: E402

EXT_DIR = RESULTS / "external"
EXT_FIG = EXT_DIR / "figures"
EXT_RAW = DATA_RAW / "seaad_mtg"
EXT_PROC = DATA_PROC / "seaad_mtg"
for _p in (EXT_DIR, EXT_FIG, EXT_RAW, EXT_PROC):
    _p.mkdir(parents=True, exist_ok=True)

EXT_SEED = int(SEED)  # 20260914
BOOT_SEED = 20260918  # POSCTRL S3 donor-bootstrap seed
N_PERM = 200
N_BOOT = 200
TRANSFER_NULL_BAR = 0.05
MIN_DONORS_MAP = 40
MIN_CELLS_PB = 20  # project constant MIN_CELLS_PER_PSEUDOBULK; not retuned
MIN_N_R = int(MIN_N_TYPE)  # 8; Pearson r needs this many donors in a type
PREREG_DATE = "2026-09-17"

# Fetched in T1 / phase0 (results/tissue/t1_cxg_rescored.csv, results/phase0_cxg_search.txt).
# Live API must return this id; we do not invent a replacement.
DATASET_ID = "c2876b1b-06d8-4d96-a56b-5304f815b99a"
COLLECTION_ID = "1ca90a2d-2943-483d-b678-b809bf464c30"
TITLE_REQUIRED = "Whole Taxonomy - MTG: Seattle Alzheimer's Disease Atlas (SEA-AD)"
TISSUE_REQUIRED = "middle temporal gyrus"
OUR_DATA_URL = "https://brain-map.org/consortia/sea-ad/our-data"
CXG_TABLE = RESULTS / "tissue" / "t1_cxg_rescored.csv"
DLPFC_OBS_CSV = LOGCPM_OBS

MANIFEST_PATH = EXT_DIR / "manifest.json"
PROGRESS_PATH = ROOT / "PROGRESS_EXTERNAL.md"
FINDINGS_PATH = ROOT / "FINDINGS_EXTERNAL.md"
PREREG_FLAG = EXT_DIR / "PREREG_20260917.flag"
STOP_PATH = EXT_DIR / "STOP.txt"

# Column-name detectors. A hit must be an actual column in a file we read.
# "disease" is not a pathology score and is not used as a substitute.
PATHOLOGY_COL_RE = re.compile(
    r"(cerad|braak|\badnc\b|adnc_|_adnc|"
    r"overall\s*ad\s*neuropath|"
    r"ad\s*neuropatholog|"
    r"nia[-_ ]?aa|"
    r"neuropathologic\s*change)",
    re.I,
)
DONOR_ID_NAME_RE = re.compile(r"^(donor|donor_id|donor id|individual_id|individual id|death_id)$", re.I)

PREREG_BARS = """Written **before any DLPFC fit or SEA-AD projection**, {date}. No threshold in this block is re-tuned after numbers exist.

FINDINGS_POSCTRL.md Supersession 2–3: the DLPFC age direction transfers between HBCC and MSSM as Pearson r (gene_ridge r +0.568 / +0.578; identity-residualized r +0.573 / +0.576; n_perm=200). Uncalibrated transfer R² is a calibration artifact. FINDINGS_GTEX.md: a stable bulk cortex age direction exists in GTEx BA9 (held-out r) and is not driven by recorded logistics.

This cell asks a different question: does the **frozen** DLPFC direction, trained on all 233 DLPFC donors and never fit on SEA-AD, predict chronological age in SEA-AD single-nucleus MTG (Allen Institute; a third bank, region, and lab)?

Two directions are frozen on all 233 DLPFC donors, per cell type, before any SEA-AD score is computed:
- **raw ridge** — global z-score, per-type SVD-LOO ridge age direction (same fitter as POSCTRL gene_ridge / TARGET age companion).
- **identity-residualized ridge** — the G2c / FINDINGS_GEOMETRY Part C 0.317→0.300 fit: type-centroid identity basis on the globally z-scored matrix, residualize, refit ridge age per type.

SEA-AD donor × mapped-type pseudobulks use the same log2(CPM+1) normalization as the DLPFC Phase-1 matrix (`decomp.normalize_logcpm`). Gene overlap is reported. Missing DLPFC genes are left at the DLPFC mean (z=0) so the frozen vector is not rewritten.

**Primary cell:** donors with no or low AD pathology, using the pathology field frozen below from Stage 0 (CERAD, Braak, or ADNC/overall pathology score; not `disease`).
**Secondary cell:** all mapped donors, pathology score as a covariate (partial Pearson r: residualize score and age on the pathology score).

Metric: Pearson r vs chronological age. Donor-level permutation null n_perm={n_perm} (X held fixed; each cell starts a fresh Generator on seed {perm_seed}). Donor-bootstrap 95% percentile CI on the median-over-types r, B={n_boot}, seed {boot_seed}. Also report Spearman ρ, calibrated R² (intercept and slope refit on SEA-AD), and uncalibrated R² (DLPFC 1-D calibration applied to SEA-AD). Aggregation: median over mapped cell types with n≥{min_n} donors (POSCTRL S2/S3). Angle is not a gate.

**Pass = r > 0 with permutation-null r ≤ {null_bar} in the primary cell.**

Pathology field frozen from Stage 0 (not chosen after seeing r):
- field: `{path_field}`
- source file: `{path_file}`
- no/low levels (primary): {low_levels}
- numeric encoding for the covariate cell: {encoding}
"""

PREREG_READING = """- **pass** (primary r > 0 with null ≤ 0.05) → the direction transfers to a third bank, region, and lab
- **fail** on the primary cell, **pass** on the pathology-covariate cell → transfer is masked by pathology
- **both fail** → does not transfer beyond the two DLPFC banks

Only the outcome that fired is the reading. Do not average primary and secondary.
"""

LOWDIM_L4_QUOTE = (
    "**The negative result stands.** What would be needed to change it: more donors per "
    "cell type (so that even gene space is not p ≫ n, or so that a k ≪ n space is estimated "
    "from hundreds of independent brains rather than ~200); **more sites** (two banks are one "
    "degree of freedom of transfer, and that transfer failed); or a **different tissue** with "
    "less bank structure and more independent donors. Re-fitting more unsupervised variants on "
    "this same 233-donor DLPFC matrix will not make an underdetermined direction unique."
)


def ext_log_banner(log, stage: str):
    log("=" * 100)
    log(f"EXTERNAL {stage}  seed={EXT_SEED}  seaad_dataset_id={DATASET_ID}")
    log("=" * 100)
    log(HASH_POOL_LIMITATION)
    log("NO perturbation data. No TF Atlas. No candidate interventions.")
    log("FALSIFICATION.md and prior FINDINGS*.md are not modified.")
    log("Fail loudly: no substitute columns, no invented cell-type synonyms.")
    log(f"Cohort: {TITLE_REQUIRED}")
    log(f"Pass bar: r > 0 with permutation-null r ≤ {TRANSFER_NULL_BAR} in the primary (no/low pathology) cell.")
    log(f"n_perm={N_PERM}  n_boot={N_BOOT}  bootstrap seed={BOOT_SEED}")


def jsonable(x):
    if isinstance(x, dict):
        return {str(k): jsonable(v) for k, v in x.items()}
    if isinstance(x, (list, tuple)):
        return [jsonable(v) for v in x]
    if isinstance(x, (np.integer,)):
        return int(x)
    if isinstance(x, (np.floating,)):
        return float(x) if np.isfinite(x) else None
    if isinstance(x, np.ndarray):
        return jsonable(x.tolist())
    if isinstance(x, (pd.Timestamp, Path)):
        return str(x)
    if x is None or isinstance(x, (str, int, float, bool)):
        if isinstance(x, float) and not np.isfinite(x):
            return None
        return x
    return str(x)


def write_progress(text: str):
    PROGRESS_PATH.write_text(text, encoding="utf-8")


def finite(x):
    try:
        return bool(np.isfinite(float(x)))
    except (TypeError, ValueError):
        return False


def load_dlpfc_types(log=print):
    """The 20 Phase-1 cell types. Missing file is fatal; the list is not invented."""
    if not DLPFC_OBS_CSV.exists():
        raise StopStep(
            "dlpfc_types",
            f"DLPFC Phase-1 obs missing: {DLPFC_OBS_CSV}. Cannot map SEA-AD types.",
            dict(path=str(DLPFC_OBS_CSV)),
        )
    obs = pd.read_csv(DLPFC_OBS_CSV)
    if "celltype" not in obs.columns:
        raise StopStep(
            "dlpfc_types",
            f"{DLPFC_OBS_CSV.name} has no celltype column. columns={list(obs.columns)}",
            dict(columns=list(obs.columns)),
        )
    types = sorted(obs.celltype.astype(str).unique().tolist())
    n_don = int(obs.donor.nunique()) if "donor" in obs.columns else None
    log(f"[dlpfc] {len(types)} cell types from {DLPFC_OBS_CSV.name}  donors={n_don}")
    if len(types) != 20:
        raise StopStep(
            "dlpfc_types",
            f"expected 20 DLPFC Phase-1 types, found {len(types)}: {types}",
            dict(n=len(types), types=types),
        )
    if n_don is not None and n_don != 233:
        raise StopStep(
            "dlpfc_types",
            f"expected 233 DLPFC donors, found {n_don}",
            dict(n_donors=n_don),
        )
    return types, obs


def require_columns(names, required, *, file_tag: str, step: str):
    have = list(names)
    missing = [c for c in required if c not in have]
    if missing:
        raise StopStep(
            step,
            f"{file_tag}: required column(s) missing or renamed: {missing}. "
            f"columns actually read={have}",
            dict(file=file_tag, missing=missing, columns_read=have),
        )
    return have


def resolve_column(names, options, *, field: str, step: str, file_tag: str):
    """Pick the unique present name from a declared list of spellings of ONE field.

    Used for Allen/CXG capitalization (subclass vs Subclass). Does not allow a
    different variable (e.g. cell_type or Supertype) to stand in.
    """
    have = list(names)
    present = [c for c in options if c in have]
    if len(present) == 0:
        raise StopStep(
            step,
            f"{file_tag}: {field} not found. looked for {list(options)}. "
            f"columns actually read={have}. Not substituting another column.",
            dict(file=file_tag, field=field, looked_for=list(options), columns_read=have),
        )
    if len(present) > 1:
        raise StopStep(
            step,
            f"{file_tag}: {field} is present under multiple names {present}. Not picking.",
            dict(file=file_tag, field=field, present=present),
        )
    return present[0]


def pathology_columns(columns):
    hits = []
    for c in columns:
        if PATHOLOGY_COL_RE.search(str(c)):
            hits.append(c)
    return hits


def pick_pathology_field(columns, log):
    """Prefer ADNC / overall neuropathologic change, then Braak, then CERAD.

    Priority is the prompt order of allowed fields, not a substitute for a missing
    field. If none of the three families exist, STOP.
    """
    hits = pathology_columns(columns)
    log(f"[pathology] columns matching CERAD/Braak/ADNC/overall: {hits or '(none)'}")
    if not hits:
        raise StopStep(
            "pathology",
            "no CERAD / Braak / ADNC / overall AD neuropathology column in the files read. "
            "Not substituting disease=normal/dementia. "
            f"columns={list(columns)}",
            dict(columns=list(columns)),
        )

    def _fam(c):
        s = str(c).lower()
        if "adnc" in s or "neuropathologic" in s or "neuropathological change" in s or "overall ad" in s:
            return 0
        if "braak" in s:
            return 1
        if "cerad" in s:
            return 2
        return 9

    ranked = sorted(hits, key=lambda c: (_fam(c), str(c)))
    chosen = ranked[0]
    fam = _fam(chosen)
    if fam == 9:
        raise StopStep(
            "pathology",
            f"matched columns {hits} but none is CERAD, Braak, or ADNC/overall. Not substituting.",
            dict(hits=hits),
        )
    family = {0: "ADNC_or_overall", 1: "Braak", 2: "CERAD"}[fam]
    log(f"[pathology] using {chosen!r} (family={family}); all hits kept in the manifest")
    return chosen, family, hits


# No/low pathology level sets. Compared to the file's actual strings after strip.
# If the file uses a level not in ANY of these sets, STOP rather than guess.
ADNC_LOW = frozenset({
    "not ad", "not_ad", "none", "absent", "no ad", "not", "0", "low",
    "not ad/low", "not ad / low", "reference", "no/low", "not ad, low",
})
BRAAK_LOW = frozenset({
    "0", "braak 0", "braak0", "i", "ii", "1", "2",
    "braak i", "braak ii", "braak 1", "braak 2",
    "stage 0", "stage i", "stage ii", "stage 1", "stage 2",
})
CERAD_LOW = frozenset({
    "0", "absent", "none", "negative", "sparse", "1",
    "cerad 0", "cerad 1", "cerad absent", "cerad none", "cerad sparse",
})
ADNC_ORD = {
    "not ad": 0, "not_ad": 0, "none": 0, "absent": 0, "no ad": 0, "not": 0, "0": 0,
    "reference": 0,  # declared in ADNC_LOW; SEA-AD neurotypical-reference level
    "low": 1, "1": 1,
    "not ad/low": 1, "not ad / low": 1, "no/low": 1, "not ad, low": 1,
    "intermediate": 2, "2": 2, "moderate": 2,
    "high": 3, "3": 3, "severe": 3, "frequent": 3,
}
BRAAK_ORD = {
    "0": 0, "braak 0": 0, "braak0": 0, "stage 0": 0,
    "i": 1, "1": 1, "braak i": 1, "braak 1": 1, "stage i": 1, "stage 1": 1,
    "ii": 2, "2": 2, "braak ii": 2, "braak 2": 2, "stage ii": 2, "stage 2": 2,
    "iii": 3, "3": 3, "braak iii": 3, "braak 3": 3, "stage iii": 3, "stage 3": 3,
    "iv": 4, "4": 4, "braak iv": 4, "braak 4": 4, "stage iv": 4, "stage 4": 4,
    "v": 5, "5": 5, "braak v": 5, "braak 5": 5, "stage v": 5, "stage 5": 5,
    "vi": 6, "6": 6, "braak vi": 6, "braak 6": 6, "stage vi": 6, "stage 6": 6,
}
CERAD_ORD = {
    "0": 0, "absent": 0, "none": 0, "negative": 0, "cerad 0": 0, "cerad absent": 0, "cerad none": 0,
    "1": 1, "sparse": 1, "cerad 1": 1, "cerad sparse": 1,
    "2": 2, "moderate": 2, "cerad 2": 2, "cerad moderate": 2,
    "3": 3, "frequent": 3, "cerad 3": 3, "cerad frequent": 3,
}


def _norm_level(v):
    try:
        if pd.isna(v):
            return None
    except (ValueError, TypeError):
        pass
    if v is None or (isinstance(v, float) and not np.isfinite(v)):
        return None
    s = str(v).strip()
    # "None" is a real CERAD level; do not treat it as missing.
    if s == "" or s.lower() in ("nan", "<na>", "."):
        return None
    return re.sub(r"\s+", " ", s.lower())


def encode_pathology(series, family: str, *, step="pathology_encode"):
    """Map file levels onto a numeric score using the declared tables. Unknown → STOP."""
    raw = pd.Series(series)
    norm = raw.map(_norm_level)
    if family == "ADNC_or_overall":
        table, low = ADNC_ORD, ADNC_LOW
    elif family == "Braak":
        table, low = BRAAK_ORD, BRAAK_LOW
    elif family == "CERAD":
        table, low = CERAD_ORD, CERAD_LOW
    else:
        raise StopStep(step, f"unknown pathology family {family!r}")
    observed = sorted(set(x for x in norm.dropna().unique()))
    unknown = [x for x in observed if x not in table]
    if unknown:
        raise StopStep(
            step,
            f"{family}: levels not in the declared encoding table: {unknown}. "
            f"observed={observed}  table={sorted(table)}. Not guessing.",
            dict(family=family, unknown=unknown, observed=observed),
        )
    score = norm.map(lambda x: np.nan if x is None else table[x])
    is_low = norm.map(lambda x: False if x is None else x in low)
    used_low = sorted(x for x in observed if x in low)
    if not used_low:
        raise StopStep(
            step,
            f"{family}: no observed level is in the declared no/low set {sorted(low)}. "
            f"observed={observed}. Not expanding the set after seeing the file.",
            dict(family=family, observed=observed, low=sorted(low)),
        )
    return score.astype(float), is_low.astype(bool), dict(
        family=family, observed=observed, unknown=[],
        n_low=int(is_low.sum()), n_missing=int(norm.isna().sum()),
        low_levels_used=used_low, encoding=table,
    )


def median_over_types(values):
    v = np.asarray(values, float)
    v = v[np.isfinite(v)]
    if v.size == 0:
        return np.nan
    return float(np.median(v))


def load_manifest():
    if not MANIFEST_PATH.exists():
        return dict(columns_read={}, files=[], failures=[], status="EMPTY")
    return load_json(MANIFEST_PATH)


def save_manifest(man: dict):
    dump_json(MANIFEST_PATH, jsonable(man))


def refresh_progress(next_action=None, stop=None, extra=""):
    man = load_manifest() if MANIFEST_PATH.exists() else {}
    stop = stop if stop is not None else man.get("status", "none")
    if next_action is None:
        if man.get("status") == "STOP":
            next_action = "STOP — write FINDINGS_EXTERNAL.md from Stage 0"
        elif not PREREG_FLAG.exists():
            next_action = "write pre-registration flag (before any fit)"
        elif not (EXT_DIR / "frozen_directions.npz").exists():
            next_action = "freeze DLPFC raw + identity-residual ridge on 233 donors"
        elif not (EXT_DIR / "primary.json").exists():
            next_action = "project SEA-AD pseudobulks; primary + secondary cells"
        else:
            next_action = "write FINDINGS_EXTERNAL.md"
    lines = [
        "# PROGRESS_EXTERNAL",
        "",
        f"**STOP status:** {stop}.",
        f"**Next action:** {next_action}",
        "",
        "## Prompt (re-read after any context summary)",
        "",
        "Goal: does the frozen DLPFC age direction predict age in an independent cohort it was never trained on? "
        "Cohort: SEA-AD single-nucleus MTG (Allen Institute), open access.",
        "",
        "Stage 0: download donor metadata and per-nucleus cell-type labels; map SEA-AD subclass → 20 DLPFC types; "
        "report unmapped types and n donors per type; report age range, sex, and AD pathology fields. "
        "Log every column name read into `results/external/manifest.json`. "
        "Stop if not downloadable or fewer than 40 donors map.",
        "",
        "Pre-register before any fit. Primary = no/low AD pathology. Secondary = all donors, pathology covariate. "
        "Pass = r > 0 with null ≤ 0.05 in the primary cell. "
        "Reading, only the outcome that fired: pass → third bank/region/lab; "
        "fail + pathology-covariate pass → masked by pathology; both fail → does not transfer beyond the two DLPFC banks.",
        "",
        "## Pre-registered bars",
        "",
        f"- pass: r > 0 with permutation-null r ≤ {TRANSFER_NULL_BAR} in the primary cell",
        f"- n_perm={N_PERM}  n_boot={N_BOOT}  perm seed `{EXT_SEED}`  boot seed `{BOOT_SEED}`",
        "- directions: raw ridge and identity-residualized ridge (0.317→0.300 fit), trained on all 233 DLPFC donors",
        "- aggregation: median over mapped types with n≥8 donors",
        "",
        "## Seeds",
        "",
        f"- main / permutation: `{EXT_SEED}`",
        f"- donor-bootstrap: `{BOOT_SEED}`",
        f"- SEA-AD dataset: `{DATASET_ID}`",
        f"- DLPFC dataset: `4442d412-91cb-4261-acca-8adf5fa04c11`",
        "",
        "## Finished cells",
        "",
    ]
    cells = []
    for p in sorted(EXT_DIR.glob("*.json")):
        if p.name in ("manifest.json",):
            continue
        cells.append(p.name)
    if not cells:
        lines.append("None.")
    else:
        for name in cells:
            lines.append(f"- `{name}`")
    lines += ["", "## Paths written", ""]
    if MANIFEST_PATH.exists():
        lines.append("- `results/external/manifest.json`")
    if PREREG_FLAG.exists():
        lines.append("- `results/external/PREREG_20260917.flag`")
    if STOP_PATH.exists():
        lines.append("- `results/external/STOP.txt`")
    lines.append("- `PROGRESS_EXTERNAL.md`")
    if FINDINGS_PATH.exists():
        lines.append("- `FINDINGS_EXTERNAL.md`")
    if extra:
        lines += ["", extra, ""]
    write_progress("\n".join(lines) + "\n")
