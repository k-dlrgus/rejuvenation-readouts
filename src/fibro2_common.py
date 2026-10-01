"""Fibroblast ruler interpretability (Prompt F). Namespaced under results/fibro2/.

Does not modify any existing FINDINGS_*.md or FALSIFICATION.md. Does not open GSE325735.
Seed 20260914. Bootstrap seed 20260918. Spearman ρ is the reported statistic.
Uncalibrated R² is reported and is never a gate. Refit nothing on external or OSKM data.
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
from brain_phase1_common import Logger, dump_json, load_json, strip_ensembl  # noqa: E402
from gtex_common import (  # noqa: E402
    StopStep, pearson_safe, spearman_safe, pred_scores, EDGE_R_PRIOR,
    tmm_norm_factors, log2_cpm_edger,
)
from target_common import unit, md_table, zscore_train  # noqa: E402
from trajectory_common import permutation_p, calibrate_1d  # noqa: E402
from fibro_common import (  # noqa: E402
    FROZEN_RULER, FIBRO_DIR, FIBRO_RAW, FIBRO_SEED, BOOT_SEED, N_PERM, N_BOOT,
    N_RANDOM_DIR, jsonable, bootstrap_rho_ci, pass_rho_bar, fmt_ci,
    SMTSD_FIBRO,
)

FIBRO2_DIR = RESULTS / "fibro2"
FIBRO2_FIG = FIBRO2_DIR / "figures"
FIBRO2_RAW = DATA_RAW / "fibro2"
FIBRO2_PROC = DATA_PROC / "fibro2"
for _p in (FIBRO2_DIR, FIBRO2_FIG, FIBRO2_RAW, FIBRO2_PROC):
    _p.mkdir(parents=True, exist_ok=True)

FIBRO2_SEED = int(SEED)  # 20260914
FIBRO2_BOOT = int(BOOT_SEED)  # 20260918
PREREG_DATE = "2026-09-18"
ADULT_MIN_AGE = 18
TRANSFER_NULL_BAR = 0.05
PCA_MAHAL_K = 50

MANIFEST_PATH = FIBRO2_DIR / "manifest.json"
PROGRESS_PATH = ROOT / "PROGRESS_FIBRO2.md"
FINDINGS_PATH = ROOT / "FINDINGS_FIBRO2.md"
PREREG_TA_FLAG = FIBRO2_DIR / "PREREG_TA.flag"
PREREG_TB_FLAG = FIBRO2_DIR / "PREREG_TB.flag"
PREREG_TC_FLAG = FIBRO2_DIR / "PREREG_TC.flag"

GSE_PRIMARY = "GSE113957"
GSE_NEAR_BULK = "GSE226189"
GSE_NEAR_N9 = "GSE307377"
CXG_FIBRO_ATLAS_ID = "a19d1667-a7b5-4556-9e5f-f9bfa690c0f1"
OUT_OF_SCOPE = ("GSE325735",)

PREREG_TA = """Written **before any T-A download, projection, or ρ**, 2026-09-18. No threshold in this block is re-tuned after numbers exist.

Take the Stage 1 ridge direction **exactly as frozen in `results/fibro/frozen_ruler_ridge_raw.npz`**. Refit nothing. Recalibrate nothing. Rescale nothing. Project, score, report Spearman ρ vs donor age.

Same normalization as the frozen GTEx spec (FINDINGS_GTEX.md Stage 0 / FINDINGS_FIBRO.md Stage 1): TMM/edgeR log2-CPM (`prior.count=2`). Genes restricted to the overlap with the frozen ruler. Missing overlap genes left at z=0 (FINDINGS_EXTERNAL.md rule). Gene-overlap count reported. z uses the frozen train mean and sd.

Metric: Spearman ρ vs the age the record states. Donor-level permutation null n_perm=200 (scores held fixed; shuffle donor ages; each cell starts a fresh Generator on seed `20260914`). Donor-bootstrap 95% percentile CI, B=200, seed `20260918`. Calibrated R² and uncalibrated R² reported; uncalibrated R² is never a gate. Pearson r reported, not a gate. Nothing averaged across cohorts.

Adult restriction on the primary cell: GSE113957 samples with stated age < 18 are excluded to match the GTEx training range; the excluded count is reported. Do not infer a donor age that the record does not state. A record that does not state per-sample donor age is unusable.

Cohorts, in priority order, each resolved from the record this session:
1. **GSE113957** — primary external cell. Human dermal fibroblast. Ages stated in titles. Bulk. Adults ≥18.
2. **Any single-cell fibroblast/skin aging cohort with ≥10 donors and stated ages** resolved this session (CELLxGENE, GEO, or SCOPING_TISSUE.md `near_miss` rows GSE226189 n=82 bulk, GSE307377 n=9). A 10x cohort is worth more than a larger bulk one because the platform shift is the thing under test. Say which platform each is.
3. If no single-cell fibroblast cohort resolves, say so plainly. Do not substitute a bulk cohort and call the platform question answered.

GSE325735 is out of scope.

**Pre-registered reading** (only the outcome that fired; T-A cannot overturn the Stage 2 `no_decline` verdict by itself):
- ρ > 0 with null ≤ 0.05 in the primary external cell → the ruler is a working instrument outside its training cohort; the Stage 2 `no_decline` verdict is interpretable as a statement about OSKM.
- ρ ≤ 0 or null > 0.05 → the ruler does not validate externally; Stage 2 is uninterpretable and `no_decline` must not be reported as a finding about reprogramming.
- bulk external passes, single-cell external fails (or none resolves) → the platform shift is unresolved; the Stage 2 projection rests on an unvalidated cross-platform step; do not resolve it by argument.
"""

PREREG_TB = """Written **before any T-B QC table or d0→d7 null**, 2026-09-18. No threshold in this block is re-tuned after numbers exist.

From the GSE297234 Cell Ranger `filtered_feature_bc_matrix.h5` matrices already on disk under `data/raw/gse297234/RAW/`, report **per donor × timepoint**: cells recovered, median UMI/cell, median genes/cell, mitochondrial read fraction, and n cells per unsupervised cluster. Table only — this is a QC audit, not a model. Unsupervised clusters are the Stage 2 spec: within each donor × timepoint independently, log1p(CP10k) on genes overlapping the frozen ruler, PCA 20 (or n_cells−1 if smaller), k-means k=3, seed `20260914`. Skip a group if n_cells < 50.

Then rerun the Stage 2 **random-direction null exactly as pre-registered** (200 permuted-weight directions, seed `20260914`) on the **d0→d7 endpoint** in addition to the d0→d10 endpoint already reported, for both donors, both pseudobulk variants (all-cell primary, per-cluster secondary). Cluster labels are assigned independently per donor × timepoint; the secondary pairing is by **within-timepoint cell-count rank** (largest / middle / smallest), not by k-means ID. Do not average cluster ranks into the all-cell primary. Do not average donors.

Endpoint decline = age_score(day 0) − age_score(day T). Positive = younger on the ruler. Empirical p = (n_random ≥ real + 1) / (n_random + 1). Pass bar remains p ≤ 0.05. A post-hoc d7 endpoint is not permitted to become the primary Stage 2 result.

**Pre-registered reading** (caveat on the verdict, not a reversal of it):
- d0→d7 also gives p > 0.05 → the `no_decline` verdict does not depend on d10; d10's behaviour is a side observation.
- d0→d7 gives p ≤ 0.05 while d0→d10 does not, **and** the QC table shows d10 differing materially from d0/d3/d7 on cells recovered, depth, or mitochondrial fraction → the endpoint choice is confounded by a technical difference at d10. Report as a caveat. State what would settle it (a cohort with more donors and a pre-registered endpoint).
- d0→d7 gives p ≤ 0.05 with no technical difference at d10 → report the discrepancy and say it is unexplained.
"""

PREREG_TC = """Written **before any T-C boosting fit**, 2026-09-18. No threshold in this block is re-tuned after numbers exist.

On GTEx cultured fibroblasts only (the Stage 1 cohort and preprocessing, unchanged: `SMAFRZE==RNASEQ`, `SMTSD == "Cells - Cultured fibroblasts"`, gene filter ≥6 counts in ≥20% of samples, TMM/edgeR log2-CPM `prior.count=2`, per-gene z-score on training rows only). D1 n is too small; excluded from B1↔C1 transfer and from within-site CV, stated, not silently dropped.

Fit `sklearn.ensemble.HistGradientBoostingRegressor`. Hyperparameters by inner grouped CV on training folds only (outer folds unchanged: within-site `SMNABTCH`-grouped; transfer is B1→C1 and C1→B1). Inner-CV selection metric is Spearman ρ on the inner held-out batches; the selected hyperparameter tuple is frozen for that train split before any permutation of that split.

Report within-site ρ and **B1→C1 / C1→B1 transfer ρ** with the same donor-level permutation null (n_perm=200, seed `20260914`) and donor-bootstrap 95% CI (B=200, seed `20260918`), side by side with the ridge numbers from `results/fibro/stage1_cv.csv` and `results/fibro/stage1_transfer.csv`. Uncalibrated R² reported, never a gate.

This is a **diagnostic, not a replacement**. A nonlinear model yields a score, not a direction; the project's plane, angles, and "age movement per unit identity loss" criterion all require a direction, so a nonlinear model cannot be substituted into Stage 2 without discarding the geometry. Do not rescore GSE297234 with boosting.

**Pre-registered reading:**
- boosting transfer ρ within the ridge bootstrap CI (both directions) → linearity is not the limiting factor; ridge stays.
- boosting transfer ρ materially above the ridge CI (upper bound exceeded in both directions) → the linear direction is leaving signal on the table. Report the gap. Do not rerun Stage 2 with boosting. A direction-preserving nonlinear extension (e.g. a kernel or autoencoder latent with a linear age axis inside it) would be the next design question.
- boosting transfer ρ below ridge → note it and move on.
"""

HGB_GRID = (
    dict(max_depth=3, learning_rate=0.05, max_iter=80, min_samples_leaf=20, l2_regularization=0.0),
    dict(max_depth=3, learning_rate=0.10, max_iter=80, min_samples_leaf=20, l2_regularization=0.0),
    dict(max_depth=6, learning_rate=0.05, max_iter=80, min_samples_leaf=20, l2_regularization=1.0),
    dict(max_depth=6, learning_rate=0.10, max_iter=80, min_samples_leaf=20, l2_regularization=1.0),
)

STAGE2_VERDICT_SENTENCE = (
    "age score does not decline → the ruler reads a static donor property, not a "
    "modifiable state. Step 3 is not supported by this data."
)


def fibro2_log_banner(log, stage: str):
    log("=" * 100)
    log(f"FIBRO2 {stage}  seed={FIBRO2_SEED}  boot_seed={FIBRO2_BOOT}  "
        f"n_perm={N_PERM}  n_boot={N_BOOT}  n_random={N_RANDOM_DIR}")
    log("=" * 100)
    log("FALSIFICATION.md and existing FINDINGS*.md are not modified.")
    log(f"Frozen ruler: {FROZEN_RULER} exists={FROZEN_RULER.exists()}")
    log("Primary metric: Spearman ρ (n_perm=200); uncalibrated R² is never a gate.")
    log("Refit nothing on external or GSE297234. Do not open GSE325735.")
    log("Fail loudly: no substitute columns, no inferred ages, no dropping rows to pass an assertion.")


def load_manifest():
    if MANIFEST_PATH.exists():
        return load_json(MANIFEST_PATH)
    return dict(
        seed=FIBRO2_SEED, boot_seed=FIBRO2_BOOT, files={}, failures=[],
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


def record_unusable(accession, reason, details=None):
    man = load_manifest()
    rec = dict(accession=accession, reason=_clip_msg(reason, 800), details=jsonable(details or {}))
    existing = man.setdefault("unusable", [])
    if not any(u.get("accession") == accession and u.get("reason") == rec["reason"] for u in existing):
        existing.append(rec)
    save_manifest(man)
    return rec


def record_unverified(accession, reason, details=None):
    man = load_manifest()
    rec = dict(accession=accession, reason=_clip_msg(reason, 800), details=jsonable(details or {}))
    existing = man.setdefault("unverified", [])
    if not any(u.get("accession") == accession and u.get("reason") == rec["reason"] for u in existing):
        existing.append(rec)
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
    write_prereg_flag(PREREG_TA_FLAG, PREREG_TA)
    write_prereg_flag(PREREG_TB_FLAG, PREREG_TB)
    write_prereg_flag(PREREG_TC_FLAG, PREREG_TC)
    return PREREG_TA_FLAG, PREREG_TB_FLAG, PREREG_TC_FLAG


def load_frozen_ruler():
    if not FROZEN_RULER.exists():
        raise StopStep("frozen_ruler", f"missing {FROZEN_RULER}")
    return np.load(FROZEN_RULER, allow_pickle=True)


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


def geo_nnn(acc: str) -> str:
    acc = str(acc).upper().strip()
    if not acc.startswith("GSE"):
        raise StopStep("geo", f"not a GSE accession: {acc!r}")
    num = int(acc[3:])
    return f"GSE{num // 1000}nnn"


def strict_age_years(text):
    """Return an integer year iff the record text states a numeric age. Else None.

    Accepts '72-year-old stage', 'age: 72', 'AGE22', '22yr', '22 y', 'donor age 79'.
    Rejects bare integers, 'adult', '80 and over' as a substitute for a year.
    """
    if text is None or (isinstance(text, float) and not np.isfinite(text)):
        return None
    s = str(text).strip()
    if not s or s.lower() in {"nan", "none", "na", ""}:
        return None
    if re.search(r"and over|or older|over stage|adult stage|elderly|octogenarian", s, re.I):
        if re.search(r"and over|or older|over stage", s, re.I):
            return None
    patterns = (
        r"(?:donor\s+)?age\s*(?:\(years\))?\s*[:_\s]\s*(\d{1,3})\b",
        r"\bAGE\s*[_:-]?\s*(\d{1,3})\b",
        r"\b(\d{1,3})\s*[-_]year[-_]old\b",
        r"\b(\d{1,3})\s*year[- ]old\b",
        r"\b(\d{1,3})\s*(?:years?|yrs?|yr)\b",
        r"\((\d{1,3})\s*yr\)",
        r"\b(\d{1,3})\s*y(?:o|o\.|r old)\b",
    )
    for pat in patterns:
        m = re.search(pat, s, re.I)
        if m:
            v = int(m.group(1))
            if 0 <= v <= 120:
                return v
    if re.fullmatch(r"\d{1,3}", s):
        v = int(s)
        if 0 <= v <= 120:
            return v
    return None


def align_counts_to_ruler(C, gene_ids, symbols, frozen, log, tag=""):
    """C is genes × samples (dense). Return samples × n_ruler counts; missing genes stay 0."""
    C = np.asarray(C, np.float64)
    if C.ndim != 2:
        raise StopStep("align", f"{tag}: counts not 2-d")
    ruler_ens = np.asarray(frozen["ensembl"]).astype(str)
    ruler_sym = np.array([str(s).upper() for s in np.asarray(frozen["symbol"])], dtype=object)
    gid = np.asarray(gene_ids).astype(str)
    ens = np.array([strip_ensembl(g) for g in gid])
    if symbols is None:
        sym_u = np.array([str(g).upper() for g in gid], dtype=object)
    else:
        sym_u = np.array([str(s).upper() for s in np.asarray(symbols)], dtype=object)
    is_ens = np.array([e.startswith("ENS") for e in ens])
    pos_e, pos_s = {}, {}
    for i, e in enumerate(ens):
        if is_ens[i]:
            pos_e.setdefault(e, i)
        pos_s.setdefault(sym_u[i], i)
    n_r = len(ruler_ens)
    idx = np.full(n_r, -1, dtype=int)
    n_ens = n_sym = 0
    for j, e in enumerate(ruler_ens):
        i = pos_e.get(e)
        if i is not None:
            idx[j] = i
            n_ens += 1
            continue
        i = pos_s.get(ruler_sym[j])
        if i is not None:
            idx[j] = i
            n_sym += 1
    n_ov = int((idx >= 0).sum())
    log(f"[align {tag}] ruler_genes={n_r} overlap={n_ov} ensembl={n_ens} "
        f"symbol_fallback={n_sym} missing={n_r - n_ov}")
    n_samp = int(C.shape[1])
    Y = np.zeros((n_samp, n_r), dtype=np.float64)
    take_r = np.flatnonzero(idx >= 0)
    if take_r.size:
        Y[:, take_r] = C[idx[take_r], :].T
    rec = dict(
        n_ruler=n_r, n_overlap=n_ov, n_ensembl=n_ens, n_symbol_fallback=n_sym,
        n_missing=int(n_r - n_ov), tag=tag,
    )
    return Y, rec


def tmm_logcpm_rows(Y, log, tag=""):
    """Y is samples × genes counts. TMM among these samples, then log2-CPM."""
    Y = np.asarray(Y, np.float64)
    nf = tmm_norm_factors(Y.T)
    logcpm = log2_cpm_edger(Y.T, nf).T
    log(f"[tmm {tag}] n_samples={Y.shape[0]} n_genes={Y.shape[1]} "
        f"nf_median={float(np.median(nf)):.4f} prior.count={EDGE_R_PRIOR:g}")
    return logcpm, nf


def score_frozen(logcpm, frozen, counts=None):
    mu = np.asarray(frozen["mu"], float)
    sd = np.asarray(frozen["sd"], float)
    sd = np.where(sd < 1e-12, 1.0, sd)
    w = np.asarray(frozen["w"], float)
    w, _ = unit(w)
    if counts is not None:
        missing = (np.asarray(counts) == 0).all(0)
    else:
        missing = (np.asarray(logcpm) == 0).all(0)
    Z = (np.asarray(logcpm, float) - mu) / sd
    Z[:, missing] = 0.0
    return Z @ w, Z, missing


def rho_null_fixed_scores(y, pred, n_perm=N_PERM, seed=FIBRO2_SEED):
    """Donor-level permutation of y; pred held fixed. Returns rho, rho_null, p, nulls."""
    y = np.asarray(y, float)
    pred = np.asarray(pred, float)
    sc = pred_scores(y, pred)
    rng = np.random.default_rng(int(seed))
    nulls = []
    n = int(len(y))
    for _ in range(int(n_perm)):
        yp = y.copy()
        rng.shuffle(yp)
        nulls.append(spearman_safe(yp, pred))
    nulls = np.asarray(nulls, float)
    rho_null = float(np.nanmedian(nulls)) if np.isfinite(nulls).any() else np.nan
    p = permutation_p(sc["rho"], nulls, greater=True)
    return sc, rho_null, p, nulls


def donor_permute_y(y, donor, rng):
    """Shuffle unique donor ages and map back. Age must be constant within donor."""
    y = np.asarray(y, float).copy()
    donor = np.asarray(donor).astype(str)
    d_u, inv = np.unique(donor, return_inverse=True)
    age_d = np.empty(len(d_u), dtype=np.float64)
    for i, d in enumerate(d_u):
        vals = y[donor == d]
        if float(np.max(vals) - np.min(vals)) > 1e-9:
            raise StopStep("age_perm", f"donor {d} has non-constant age")
        age_d[i] = float(vals[0])
    rng.shuffle(age_d)
    return age_d[inv]


def rho_null_donor(y, pred, donor, n_perm=N_PERM, seed=FIBRO2_SEED):
    y = np.asarray(y, float)
    pred = np.asarray(pred, float)
    sc = pred_scores(y, pred)
    rng = np.random.default_rng(int(seed))
    nulls = []
    for _ in range(int(n_perm)):
        yp = donor_permute_y(y, donor, rng)
        nulls.append(spearman_safe(yp, pred))
    nulls = np.asarray(nulls, float)
    rho_null = float(np.nanmedian(nulls)) if np.isfinite(nulls).any() else np.nan
    p = permutation_p(sc["rho"], nulls, greater=True)
    return sc, rho_null, p, nulls


def looks_like_counts(arr, frac=0.999, max_check=2_000_000):
    v = np.asarray(arr, float).ravel()
    v = v[np.isfinite(v)]
    if v.size == 0:
        return False
    if v.size > max_check:
        rng = np.random.default_rng(FIBRO2_SEED)
        v = rng.choice(v, size=max_check, replace=False)
    frac_int = float(np.mean(np.abs(v - np.round(v)) < 1e-6))
    return bool(frac_int >= frac and float(np.nanmin(v)) >= -1e-6)


def _read_csv_if(path: Path):
    if path.exists():
        return pd.read_csv(path)
    return None


def progress_snapshot(next_action, stop=None, extra=""):
    man = load_manifest() if MANIFEST_PATH.exists() else {}
    stop_s = stop if stop else man.get("status", "running")
    lines = [
        "# PROGRESS_FIBRO2",
        "",
        f"**STOP status:** {stop_s}",
        f"**Next action:** {next_action}",
        "",
        "## Seeds / gates",
        "",
        f"- seed `{FIBRO2_SEED}`  boot `{FIBRO2_BOOT}`  n_perm={N_PERM}  n_boot={N_BOOT}  n_random={N_RANDOM_DIR}",
        f"- primary metric: Spearman ρ; uncalibrated R² is never a gate",
        f"- frozen ruler: `{FROZEN_RULER}` exists={FROZEN_RULER.exists()}",
        f"- PREREG_TA.flag exists={PREREG_TA_FLAG.exists()}",
        f"- PREREG_TB.flag exists={PREREG_TB_FLAG.exists()}",
        f"- PREREG_TC.flag exists={PREREG_TC_FLAG.exists()}",
        f"- GSE325735: out of scope",
        "",
        "## Frozen ruler path",
        "",
        f"- `{FROZEN_RULER}` exists={FROZEN_RULER.exists()}",
        "",
        "## Pre-registered blocks (verbatim)",
        "",
        "### T-A",
        "",
        PREREG_TA,
        "",
        "### T-B",
        "",
        PREREG_TB,
        "",
        "### T-C",
        "",
        PREREG_TC,
        "",
        "## Resolved accessions",
        "",
    ]
    resolved = man.get("resolved") or []
    if not resolved:
        lines.append("(none yet)")
        lines.append("")
    else:
        for r in resolved:
            lines.append(f"- {json.dumps(jsonable(r), ensure_ascii=False)}")
        lines.append("")
    lines += ["## Finished cells (from disk)", ""]
    for name, path in (
        ("T-A cohort table", FIBRO2_DIR / "ta_cohorts.csv"),
        ("T-A results", FIBRO2_DIR / "ta_results.csv"),
        ("T-B QC", FIBRO2_DIR / "tb_qc.csv"),
        ("T-B null", FIBRO2_DIR / "tb_null.csv"),
        ("T-C within-site", FIBRO2_DIR / "tc_cv.csv"),
        ("T-C transfer", FIBRO2_DIR / "tc_transfer.csv"),
        ("extrapolation", FIBRO2_DIR / "extrap_distances.csv"),
    ):
        df = _read_csv_if(path)
        if df is None or not len(df):
            continue
        lines.append(f"### {name} (`{path.name}`)")
        for _, r in df.iterrows():
            bits = [f"{c}={r[c]}" for c in df.columns[:18]]
            lines.append("- " + " ".join(bits))
        lines.append("")
    for js_name in ("ta_summary.json", "tb_summary.json", "tc_summary.json", "extrap_summary.json"):
        p = FIBRO2_DIR / js_name
        if p.exists():
            rec = load_json(p)
            lines.append(f"### {js_name}")
            lines.append(f"- {json.dumps(jsonable(rec), ensure_ascii=False)[:4000]}")
            lines.append("")
    if extra:
        lines += ["## Note", "", extra, ""]
    if man.get("unusable"):
        lines += ["## unusable", ""]
        for u in man["unusable"]:
            lines.append(f"- **{u.get('accession')}:** {u.get('reason')}")
        lines.append("")
    if man.get("unverified"):
        lines += ["## unverified", ""]
        for u in man["unverified"]:
            lines.append(f"- **{u.get('accession')}:** {u.get('reason')}")
        lines.append("")
    if man.get("failures"):
        lines += ["## Failures (manifest)", ""]
        for f in man["failures"]:
            lines.append(f"- **{f.get('step')}:** {_clip_msg(f.get('message'), 240)}")
        lines.append("")
    lines += [
        "## Files",
        "",
        "- `src/fibro2_common.py`, `fibro2_ta.py`, `fibro2_tb.py`, `fibro2_tc.py`, `fibro2_extrap.py`, `fibro2_findings.py`, `fibro2_run.py`",
        "- `results/fibro2/`",
        "- `FINDINGS_FIBRO2.md`",
        "- `PROGRESS_FIBRO2.md`",
        "",
    ]
    PROGRESS_PATH.write_text("\n".join(lines) + "\n", encoding="utf-8")
    return PROGRESS_PATH
