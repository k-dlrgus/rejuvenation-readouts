"""Fibroblast age ruler (GTEx V10) and OSKM projection (GSE297234).

Namespaced under results/fibro/. Does not modify any existing FINDINGS*.md
or FALSIFICATION.md. Seed 20260914.

Stage 1 pass = Spearman ρ > 0 in both SMCENTER transfer directions
(B1→C1 and C1→B1) with permutation-null ρ ≤ 0.05, ridge, raw.
Uncalibrated R² is reported and is never a gate. PLS-1 is reported and is
never averaged with ridge. Stage 2 runs only if Stage 1 passes.
"""
from __future__ import annotations

import json
import sys
from pathlib import Path

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parent))
from config import RESULTS, DATA_RAW, DATA_PROC, ROOT, SEED  # noqa: E402
from brain_phase1_common import Logger, dump_json, load_json, strip_ensembl  # noqa: E402
from gtex_common import (  # noqa: E402
    StopStep, pearson_safe, spearman_safe, pred_scores, pass_r_bar,
    fmt, GTEX_RAW, GTEX_PROC, FILES, SMAFRZE_RNASEQ, SAMPLE_REQUIRED,
    PHENO_REQUIRED, EXPECTED_SEX, EXPECTED_DTHHRDY, EXPECTED_AGE_BINS,
    AGE_MID, RELEASE, RNASEQ_PIPELINE, DBGAP, MIN_COUNT, MIN_FRAC,
    EDGE_R_PRIOR, sampid_to_subj, require_columns, require_levels,
    age_midpoint, tmm_norm_factors, log2_cpm_edger, gene_filter_mask,
    load_dlpfc_gene_universe, assert_disjoint, grouped_kfold,
    residualize_train_test, onehot_train, apply_onehot, dist_summary,
    DLPFC_GENES_CSV,
)
from target_common import (  # noqa: E402
    age_direction, zscore_train, unit, align_age_sign, ridge_prestd, md_table,
)
from trajectory_common import permutation_p, summarize_null_col, calibrate_1d  # noqa: E402

FIBRO_DIR = RESULTS / "fibro"
FIBRO_FIG = FIBRO_DIR / "figures"
FIBRO_RAW = DATA_RAW / "gse297234"
FIBRO_PROC = DATA_PROC / "fibro"
for _p in (FIBRO_DIR, FIBRO_FIG, FIBRO_RAW, FIBRO_PROC):
    _p.mkdir(parents=True, exist_ok=True)

FIBRO_SEED = int(SEED)  # 20260914
BOOT_SEED = 20260918
N_PERM = 200
N_BOOT = 200
N_FOLDS = 5
TRANSFER_NULL_BAR = 0.05
N_RANDOM_DIR = 200
PREREG_DATE = "2026-09-17"

SMTSD_FIBRO = "Cells - Cultured fibroblasts"
TRANSFER_SITES = ("B1", "C1")
EXCLUDE_SITE = "D1"
EXCLUDE_SITE_REASON = (
    "D1 n is too small for a transfer cell; excluded from B1↔C1 transfer and from "
    "within-site CV, stated, not silently dropped from the cohort table."
)
SITE_CV_MIN_N = 10
TRANSFER_MIN_N = 25

MANIFEST_PATH = FIBRO_DIR / "manifest.json"
PROGRESS_PATH = ROOT / "PROGRESS_FIBRO.md"
FINDINGS_PATH = ROOT / "FINDINGS_FIBRO.md"
PREREG_STAGE1_FLAG = FIBRO_DIR / "PREREG_STAGE1.flag"
PREREG_STAGE2_FLAG = FIBRO_DIR / "PREREG_STAGE2.flag"
FROZEN_RULER = FIBRO_DIR / "frozen_ruler_ridge_raw.npz"

# Frozen Stage 2 gene lists (written into the Stage 2 flag before any score).
PLURI_ENDOGENOUS = ("POU5F1", "NANOG", "LIN28A", "SALL4", "DPPA4", "ZFP42", "DNMT3B")
FIBRO_IDENTITY = ("COL1A1", "COL1A2", "THY1", "S100A4", "FN1", "VIM", "POSTN")
OSKM_FAMILY = ("POU5F1", "OCT4", "OCT3", "SOX2", "KLF4", "MYC", "MYCL", "MYCN")
GSE_ACCESSION = "GSE297234"
GSE_EXPECTED_DONORS = (
    dict(cell_line="GM00731", age_years=96, role="aged"),
    dict(cell_line="GM23815", age_years=22, role="young"),
)
GSE_EXPECTED_DAYS = (0, 3, 7, 10)
GSE_SUPPL = (
    "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE297nnn/GSE297234/suppl/GSE297234_RAW.tar",
    "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE297nnn/GSE297234/suppl/GSE297234_GM00731_SEVOSKM.rds",
    "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE297nnn/GSE297234/suppl/GSE297234_HFIB_COMBINED_SEVOSKM.rds",
)

CAVEAT_AGE_BIN = (
    "GTEx public `AGE` is a 10-year bin, so ρ is against a coarse ordinal target "
    "and year-level error is not meaningful."
)

PREREG_STAGE1 = """Written **before any Stage 1 fit**, 2026-09-17. No threshold in this block is re-tuned after numbers exist.

Cohort: GTEx Analysis V10 RNASEQ freeze (`SMAFRZE==RNASEQ`), `SMTSD == "Cells - Cultured fibroblasts"`.
One sample per donor. Target = public `AGE` bin midpoint. Model = SVD-LOO ridge (same fitter as
POSCTRL gene_ridge / TARGET age companion). PLS-1 is fit and reported; it is **not** averaged with
ridge (GTEx T3: PLS-1 has an inflated null).

Preprocessing (FINDINGS_GTEX.md Stage 0, verbatim): gene filter ≥6 counts in ≥20% of samples;
TMM/edgeR log2-CPM (`prior.count=2`); per-gene z-score on training rows only.

Covariate regimes `raw` and `resid`. `resid` columns: `SMRIN`, `SMTSISCH`, `DTHHRDY` one-hot, `SEX`,
residualized on training rows only. `resid` is reported alongside `raw`; it is not a gate.

Within-site CV grouped by `SMNABTCH` (batch-held-out). Donor-held-out is asserted in the same folds
(one sample per donor). Metrics: Spearman ρ, donor-level permutation null (n_perm=200, seed `20260914`),
donor-bootstrap 95% percentile CI (B=200, seed `20260918`), calibrated R², uncalibrated R².
Pearson r is reported and is not a gate.

Transfer: B1→C1 and C1→B1 (`SMCENTER`). D1 is too small and is excluded, stated, not silently dropped.

**Stage 1 pass = ρ > 0 in both transfer directions with permutation-null ρ ≤ 0.05, in `raw`, ridge.**
If Stage 1 fails, STOP. Do not open GSE297234. The fibroblast ruler does not transfer between GTEx
sites, so a reprogramming projection would be uninterpretable.

GTEx public `AGE` is a 10-year bin, so ρ is against a coarse ordinal target and year-level error is
not meaningful. Gene-universe overlap with the DLPFC 25526-gene matrix (Ensembl, version-stripped)
is reported before Stage 2.
"""

PREREG_STAGE2 = """Written **before the first GSE297234 projection**, 2026-09-17. No threshold in this block is re-tuned after numbers exist. Stage 2 runs only if Stage 1 passed. Nothing is fitted, refitted, or calibrated on GSE297234.

Frozen age ruler: SVD-LOO ridge, `raw`, fit on **all** GTEx V10 RNASEQ cultured-fibroblast donors (including D1) after Stage 1 passed. Weights, train mean, and train sd are frozen. Missing overlap genes are left at z=0 (FINDINGS_EXTERNAL.md rule). Same TMM/edgeR log2-CPM (`prior.count=2`) as the GTEx matrix, genes restricted to the overlap.

Dataset: GSE297234. Verify the sample table from the record before use. Expected (scoping): Sendai OSKM, primary dermal fibroblasts, days 0 / 3 / 7 / 10, donors GM23815 age 22 and GM00731 age 96. If the record differs, report the discrepancy and use the record. Do not infer a donor age that is not in the record.

Two coordinates:
- Age: pseudobulk per (donor × timepoint). (a) all-cell pseudobulk is primary. (b) unsupervised cluster within each donor × timepoint (PCA 20 PCs on log1p-CP10k of overlap genes, k-means k=3, seed `20260914`; skip a group if n_cells < 50). (b) is reported and never averaged into (a).
- Differentiation / pluripotency: mean z of frozen endogenous set (POU5F1, NANOG, LIN28A, SALL4, DPPA4, ZFP42, DNMT3B) minus mean z of frozen fibroblast-identity set (COL1A1, COL1A2, THY1, S100A4, FN1, VIM, POSTN). z is per-gene across the scored (a) pseudobulks (and separately for (b)). Lists frozen before scoring.

Sendai transgene contamination is a hard requirement. Before any score: check whether the processed matrix contains Sendai/vector features; report what is found. Report the pluripotency score both with and without POU5F1 and any other OSKM-family gene in the endogenous set. If the two disagree, the without-OSKM version is primary.

Day-0 sanity check: with no reprogramming yet, the 96-year donor must score higher on the age ruler than the 22-year donor. Two donors is one comparison and proves nothing statistically — but if it fails, report it and treat every downstream number as suspect.

Primary question: in the aged donor (GM00731), does the age score decline monotonically across days 0 → 3 → 7 → 10?

Primary null — random directions: draw 200 unit vectors in the same gene space, matched to the frozen ruler's weight distribution (permute the ruler's weights across genes; seed `20260914`). Score every (a) pseudobulk on each. Report the fraction of random directions whose day-0→day-10 decline on the aged donor is at least as large as the real ruler's. **Pass = empirical p ≤ 0.05.**

Secondary null: the same 200 permuted-weight directions against the monotonicity of the aged-donor trajectory (Spearman ρ of score vs day), not just the endpoint drop.

Reading (only the outcome that fired; do not average donors, regimes, or pseudobulk variants):
- age score declines, beats the random-direction null, **and** the pluripotency score rises later than the age drop → the curve bends; there is a measurable window where age moves before identity does.
- age score declines but does **not** beat the random-direction null → uninformative; wholesale transcriptome change moves every direction, including this one. Do not report a bend.
- age score declines and pluripotency rises in lockstep → no window; the field's core assumption is challenged. Report it as such, not as a failure.
- age score does not decline → the ruler reads a static donor property, not a modifiable state. Step 3 is not supported by this data; say so plainly and stop.

Operational definitions (frozen before any GSE297234 score):
- Age endpoint decline = age_score(day 0) − age_score(day 10) on the aged donor, all-cell pseudobulk (a). Positive = younger on the ruler.
- Age drop time = earliest d ∈ {3, 7, 10} with age_score(d) < age_score(0). If none, no drop.
- Pluripotency rise time = earliest d ∈ {3, 7, 10} with pluri(d) > pluri(0). If none, no rise.
- "rises later than the age drop" = age drop time exists AND pluri rise time exists AND pluri rise time > age drop time.
- "lockstep" = both times exist AND they are equal.
- "age score does not decline" = endpoint decline ≤ 0.
- Two series of pluripotency "disagree" if the day-0→day-10 sign differs or the rise times differ; then without-OSKM is primary.
- Unsupervised clusters: within each donor × timepoint independently, log1p(CP10k) on genes overlapping the frozen ruler, PCA 20 (or n_cells−1 if smaller), k-means k=3, seed `20260914`. Skip a group if n_cells < 50.

The young donor (GM23815) is reported alongside as a contrast, never pooled with the aged donor.
Do not open GSE325735.
"""


def fibro_log_banner(log, stage: str):
    log("=" * 100)
    log(f"FIBRO {stage}  seed={FIBRO_SEED}  boot_seed={BOOT_SEED}  "
        f"release={RELEASE}  pipeline={RNASEQ_PIPELINE}  {DBGAP}")
    log("=" * 100)
    log("FALSIFICATION.md and existing FINDINGS*.md are not modified.")
    log(f"SMTSD exact: {SMTSD_FIBRO!r}  freeze={SMAFRZE_RNASEQ}")
    log("Primary transfer metric: Spearman ρ (n_perm=200); uncalibrated R² is never a gate.")
    log("PLS-1 is reported and is not averaged with ridge.")
    log("Fail loudly: no substitute columns, no silent coerce, no dropping rows to pass an assertion.")
    log("Do not infer a donor age that is not in the record.")
    log("Do not open GSE325735.")


def pass_rho_bar(rho, rho_null, bar=TRANSFER_NULL_BAR):
    """Pass = ρ > 0 with permutation-null ρ ≤ 0.05."""
    try:
        rv, nv = float(rho), float(rho_null)
    except (TypeError, ValueError):
        return False
    return bool(np.isfinite(rv) and rv > 0 and np.isfinite(nv) and nv <= float(bar))


def fmt_ci(lo, hi, d=3):
    try:
        a, b = float(lo), float(hi)
    except (TypeError, ValueError):
        return "[NA, NA]"
    if not (np.isfinite(a) and np.isfinite(b)):
        return "[NA, NA]"
    return f"[{a:+.{d}f}, {b:+.{d}f}]"


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


def jsonable(obj):
    if isinstance(obj, dict):
        return {str(k): jsonable(v) for k, v in obj.items()}
    if isinstance(obj, (list, tuple)):
        return [jsonable(v) for v in obj]
    if isinstance(obj, Path):
        return str(obj)
    if isinstance(obj, (np.integer,)):
        return int(obj)
    if isinstance(obj, (np.floating,)):
        return float(obj) if np.isfinite(obj) else None
    if isinstance(obj, (np.bool_,)):
        return bool(obj)
    if isinstance(obj, np.ndarray):
        return jsonable(obj.tolist())
    if obj is None or isinstance(obj, (str, int, float, bool)):
        if isinstance(obj, float) and not np.isfinite(obj):
            return None
        return obj
    return str(obj)


def load_manifest():
    if MANIFEST_PATH.exists():
        return load_json(MANIFEST_PATH)
    return dict(seed=FIBRO_SEED, boot_seed=BOOT_SEED, files={}, failures=[],
                columns={}, status="INIT")


def save_manifest(man):
    dump_json(MANIFEST_PATH, jsonable(man))
    return MANIFEST_PATH


def record_failure(step, message, details=None):
    man = load_manifest()
    rec = dict(step=step, message=message, details=jsonable(details or {}))
    man.setdefault("failures", []).append(rec)
    man["status"] = "STOP"
    save_manifest(man)
    return rec


def bootstrap_rho_ci(y, pred, rng, n_boot=N_BOOT):
    """Donor-level percentile CI on Spearman ρ (one row per donor)."""
    y = np.asarray(y, float)
    pred = np.asarray(pred, float)
    m = np.isfinite(y) & np.isfinite(pred)
    y, pred = y[m], pred[m]
    n = int(len(y))
    vals = []
    if n < 3:
        return dict(p025=np.nan, p975=np.nan, n_ok=0, n=n)
    for _ in range(int(n_boot)):
        idx = rng.integers(0, n, size=n)
        r = spearman_safe(y[idx], pred[idx])
        if np.isfinite(r):
            vals.append(float(r))
    v = np.asarray(vals, float)
    if v.size < 2:
        return dict(p025=np.nan, p975=np.nan, n_ok=int(v.size), n=n)
    return dict(
        p025=float(np.percentile(v, 2.5)),
        p975=float(np.percentile(v, 97.5)),
        median=float(np.median(v)),
        n_ok=int(v.size),
        n=n,
    )


def write_prereg_stage1():
    if PREREG_STAGE1_FLAG.exists():
        return PREREG_STAGE1_FLAG
    PREREG_STAGE1_FLAG.write_text(PREREG_STAGE1, encoding="utf-8")
    return PREREG_STAGE1_FLAG


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


def _read_csv_if(path: Path):
    if path.exists():
        return pd.read_csv(path)
    return None


def progress_snapshot(next_action, stop=None, extra=""):
    """Rewrite PROGRESS_FIBRO.md from results/fibro on disk. Called after every cell."""
    man = load_manifest() if MANIFEST_PATH.exists() else {}
    gate = man.get("stage1_gate") or {}
    if (FIBRO_DIR / "stage1_gate.json").exists():
        gate = load_json(FIBRO_DIR / "stage1_gate.json")
    stop_s = stop if stop else man.get("status", "running")
    lines = [
        "# PROGRESS_FIBRO",
        "",
        f"**STOP status:** {stop_s}",
        f"**Next action:** {next_action}",
        "",
        "## Stage gate",
        "",
        f"- Stage 1 prereg flag: `{PREREG_STAGE1_FLAG.name}` exists={PREREG_STAGE1_FLAG.exists()}",
        f"- Stage 2 prereg flag: `{PREREG_STAGE2_FLAG.name}` exists={PREREG_STAGE2_FLAG.exists()}",
        f"- Stage 1 pass (ridge raw both transfer directions): {gate.get('pass')}",
        f"- Stage 1 reason: {gate.get('reason', 'not scored')}",
        f"- Stage 2 opened: {bool(man.get('stage2_opened'))}",
        f"- frozen ruler: `{FROZEN_RULER.name}` exists={FROZEN_RULER.exists()}",
        f"- seed `{FIBRO_SEED}`  boot seed `{BOOT_SEED}`  n_perm={N_PERM}  n_boot={N_BOOT}",
        f"- primary metric: Spearman ρ; uncalibrated R² is never a gate",
        f"- D1 exclusion: {EXCLUDE_SITE_REASON}",
        "",
        "## Frozen preprocessing spec (FINDINGS_GTEX.md Stage 0, verbatim)",
        "",
        f"- `SMAFRZE=={SMAFRZE_RNASEQ}` and `SMTSD == {SMTSD_FIBRO!r}`",
        f"- Gene filter: ≥{MIN_COUNT} counts in ≥{MIN_FRAC:.0%} of samples",
        f"- TMM (edgeR defaults) → log2-CPM `prior.count={EDGE_R_PRIOR:g}`",
        "- z-score per gene on training rows only inside every fold (Stage 1)",
        "- Map to DLPFC gene universe, Ensembl version-stripped",
        "",
        "## Frozen gene lists (Stage 2; not scored until Stage 1 passes)",
        "",
        f"- Pluripotency endogenous: {', '.join(PLURI_ENDOGENOUS)}",
        f"- Fibroblast identity: {', '.join(FIBRO_IDENTITY)}",
        f"- OSKM-family genes dropped for the without-OSKM pluripotency score: {', '.join(OSKM_FAMILY)}",
        "",
        "## Frozen ruler path",
        "",
        f"- `{FROZEN_RULER}` exists={FROZEN_RULER.exists()}",
        "",
        "## Stage 1 pre-registration (verbatim)",
        "",
        PREREG_STAGE1,
        "",
        "## Finished cells (from disk)",
        "",
    ]
    cohort = _read_csv_if(FIBRO_DIR / "stage1_cohort.csv")
    if cohort is not None and len(cohort):
        r = cohort.iloc[0]
        lines.append("### Cohort")
        lines.append(
            f"- SMTSD={SMTSD_FIBRO!r} n_samples={int(r.get('n_samples', np.nan))} "
            f"n_donors={int(r.get('n_donors', np.nan))} "
            f"SMCENTER={r.get('SMCENTER_counts', '')} "
            f"n_SMNABTCH={int(r.get('n_SMNABTCH', np.nan)) if pd.notna(r.get('n_SMNABTCH', np.nan)) else 'NA'}"
        )
        lines.append("")
    overlap_p = FIBRO_DIR / "stage1_gene_overlap.json"
    if overlap_p.exists():
        ov = load_json(overlap_p)
        lines.append("### Gene overlap")
        lines.append(
            f"- n_genes_raw={ov.get('n_genes_raw')} n_genes_filter={ov.get('n_genes_filter')} "
            f"n_overlap_dlpfc={ov.get('n_overlap_dlpfc')} "
            f"dlpfc_n={ov.get('n_dlpfc')}"
        )
        lines.append("")
    cv = _read_csv_if(FIBRO_DIR / "stage1_cv.csv")
    if cv is not None and len(cv):
        lines.append("### Within-site CV (SMNABTCH-grouped)")
        for _, r in cv.iterrows():
            lines.append(
                f"- site={r.get('site')} regime={r.get('regime')} method={r.get('method')} "
                f"ρ={_fmt_num(r.get('rho'))} (null {_fmt_num(r.get('rho_null'))} p={_fmt_num(r.get('rho_p'))}) "
                f"CI={fmt_ci(r.get('rho_ci_lo'), r.get('rho_ci_hi'))} "
                f"calR²={_fmt_num(r.get('cal_r2'))} R²={_fmt_num(r.get('r2'))} "
                f"n={r.get('n')} n_perm={r.get('n_perm')}"
            )
        lines.append("")
    tr = _read_csv_if(FIBRO_DIR / "stage1_transfer.csv")
    if tr is not None and len(tr):
        lines.append("### Transfer B1↔C1")
        for _, r in tr.iterrows():
            lines.append(
                f"- {r.get('direction')} regime={r.get('regime')} method={r.get('method')} "
                f"ρ={_fmt_num(r.get('rho'))} (null {_fmt_num(r.get('rho_null'))} p={_fmt_num(r.get('rho_p'))}) "
                f"CI={fmt_ci(r.get('rho_ci_lo'), r.get('rho_ci_hi'))} "
                f"calR²={_fmt_num(r.get('cal_r2'))} R²={_fmt_num(r.get('r2'))} "
                f"pass_rho={r.get('pass_rho')} n_train={r.get('n_train')} n_test={r.get('n_test')} "
                f"n_perm={r.get('n_perm')}"
            )
        lines.append("")
    if gate:
        lines.append("### Stage 1 gate")
        lines.append(
            f"- pass={gate.get('pass')} method={gate.get('method')} regime={gate.get('regime')} "
            f"B1→C1 ρ={_fmt_num(gate.get('rho_B1_to_C1'))} (null {_fmt_num(gate.get('rho_null_B1_to_C1'))}) "
            f"C1→B1 ρ={_fmt_num(gate.get('rho_C1_to_B1'))} (null {_fmt_num(gate.get('rho_null_C1_to_B1'))}) "
            f"reason={gate.get('reason')}"
        )
        lines.append("")
    s2s = FIBRO_DIR / "stage2_summary.json"
    if s2s.exists():
        s2 = load_json(s2s)
        lines.append("### Stage 2 (from disk)")
        lines.append(f"- {json.dumps(jsonable(s2), ensure_ascii=False)[:4000]}")
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
        "- `src/fibro_common.py`, `fibro_stage1.py`, `fibro_stage2.py`, `fibro_findings.py`, `fibro_run.py`",
        "- `results/fibro/`",
        "- `FINDINGS_FIBRO.md`",
        "- `PROGRESS_FIBRO.md`",
        "",
    ]
    PROGRESS_PATH.write_text("\n".join(lines) + "\n", encoding="utf-8")
    return PROGRESS_PATH
