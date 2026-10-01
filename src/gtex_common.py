"""GTEx v10 cortex logistics audit (FINDINGS_GTEX.md).

Stage 0: download, QC, covariate audit. Stage 1: identifiability battery,
gated on FINDINGS_POSCTRL.md Supersession 2 outcome (i) or (ii) (planted r transfers).
Does not modify prior FINDINGS*.md or FALSIFICATION.md. Seed 20260914.
"""
from __future__ import annotations

import hashlib
import json
import sys
from pathlib import Path

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parent))
from config import RESULTS, DATA_RAW, DATA_PROC, ROOT, SEED  # noqa: E402
from brain_phase1_common import Logger, dump_json, load_json, strip_ensembl  # noqa: E402
from target_common import (  # noqa: E402
    ALPHAS, TARGET_DIR, N_BOOT_RIDGE, N_PERM, HASH_POOL_LIMITATION,
    age_direction, zscore_train, unsigned_angle_deg, pairwise_angles_deg,
    md_table, unit, align_age_sign, ridge_prestd,
)
from trajectory_common import (  # noqa: E402
    permutation_p, summarize_null_col, fmt, fmt_u, calibrate_1d,
)
from geometry_partC import pls1_direction, identity_basis, angle_to_subspace_deg, proj_r2  # noqa: E402
from geometry_common import ridge_svd_fit, ridge_svd_predict  # noqa: E402

GTEX_DIR = RESULTS / "gtex"
GTEX_FIG = GTEX_DIR / "figures"
GTEX_RAW = DATA_RAW / "gtex"
GTEX_PROC = DATA_PROC / "gtex"
for _p in (GTEX_DIR, GTEX_FIG, GTEX_RAW, GTEX_PROC):
    _p.mkdir(parents=True, exist_ok=True)

GTEX_SEED = int(SEED)  # 20260914
RELEASE = "GTEx Analysis V10"
RNASEQ_PIPELINE = "RNASeQCv2.4.2"
DBGAP = "phs000424.v10"
PORTAL_DATE = "2026-09-16"

GCS_BUCKET = "https://storage.googleapis.com/adult-gtex"
FILES = {
    "gene_reads": dict(
        url=f"{GCS_BUCKET}/bulk-gex/v10/rna-seq/GTEx_Analysis_v10_RNASeQCv2.4.2_gene_reads.gct.gz",
        name="GTEx_Analysis_v10_RNASeQCv2.4.2_gene_reads.gct.gz",
        gcs_size=941668791,
        gcs_md5_b64="+EWiWJNdzeBD9VR/GWlBtw==",
        kind="gct_gz",
    ),
    "sample_attributes": dict(
        url=f"{GCS_BUCKET}/annotations/v10/metadata-files/GTEx_Analysis_v10_Annotations_SampleAttributesDS.txt",
        name="GTEx_Analysis_v10_Annotations_SampleAttributesDS.txt",
        gcs_size=38497961,
        gcs_md5_b64="TEXiXD/F7JKpfaDUvLoSBw==",
        kind="tsv",
    ),
    "subject_phenotypes": dict(
        url=f"{GCS_BUCKET}/annotations/v10/metadata-files/GTEx_Analysis_v10_Annotations_SubjectPhenotypesDS.txt",
        name="GTEx_Analysis_v10_Annotations_SubjectPhenotypesDS.txt",
        gcs_size=20292,
        gcs_md5_b64="/6FaaAhVxshPJVG/uoAXug==",
        kind="tsv",
    ),
}

# Six tissues, no more. Exact SMTSD strings; do not substitute.
TISSUES = (
    "Brain - Frontal Cortex (BA9)",
    "Brain - Cortex",
    "Brain - Anterior cingulate cortex (BA24)",
    "Brain - Hippocampus",
    "Heart - Left Ventricle",
    "Muscle - Skeletal",
)
PRIMARY_TISSUE = "Brain - Frontal Cortex (BA9)"
PAIRED_TISSUE = "Brain - Cortex"
CROSS_REGION = (
    "Brain - Anterior cingulate cortex (BA24)",
    "Brain - Hippocampus",
)
NONBRAIN = (
    "Heart - Left Ventricle",
    "Muscle - Skeletal",
)
BRAIN_TISSUES = (
    "Brain - Frontal Cortex (BA9)",
    "Brain - Cortex",
    "Brain - Anterior cingulate cortex (BA24)",
    "Brain - Hippocampus",
)
TISSUE_SLUG = {
    "Brain - Frontal Cortex (BA9)": "ba9",
    "Brain - Cortex": "cortex",
    "Brain - Anterior cingulate cortex (BA24)": "ba24",
    "Brain - Hippocampus": "hippocampus",
    "Heart - Left Ventricle": "heart_lv",
    "Muscle - Skeletal": "muscle",
}

SMAFRZE_RNASEQ = "RNASEQ"
MIN_COUNT = 6
MIN_FRAC = 0.20
TMM_LOGRATIO_TRIM = 0.30
TMM_SUM_TRIM = 0.05
EDGE_R_PRIOR = 2.0
N_BOOT = int(N_BOOT_RIDGE)  # 20
N_PERM_GTEX = 200  # pre-registered 2026-09-17; transfer Pearson r null
N_FOLDS = 5
SEX_ANGLE_REF = 41.3
PLANTED_ANGLE_N233 = 41.191  # C3 planted dense_shared n=233 bootstrap (nearest to BA9 n=269)
TRANSFER_NULL_BAR = 0.05
BA9_N_UNDERPOWERED = 150
CENTER_TRANSFER_MIN_N = 25
CENTER_TRANSFER_MIN_LEVELS = 3
T6_MIN_OVERLAP = 15000
PREREG_DATE = "2026-09-17"
# Per-cell RNGs so T-cells are independently resumable (n_perm=200 must not shift later cells).
CELL_SEEDS = dict(T1=GTEX_SEED + 1, T2=GTEX_SEED + 2, T3=GTEX_SEED + 3,
                  T4=GTEX_SEED + 4, T5=GTEX_SEED + 5, T6=GTEX_SEED + 6,
                  COMP=GTEX_SEED + 7)

SAMPLE_REQUIRED = (
    "SAMPID", "SMTSD", "SMAFRZE", "SMCENTER", "SMNABTCH", "SMGEBTCH",
    "SMRIN", "SMTSISCH",
)
PHENO_REQUIRED = ("SUBJID", "SEX", "AGE", "DTHHRDY")
EXPECTED_SEX = frozenset({"1", "2"})
EXPECTED_DTHHRDY = frozenset({"0", "1", "2", "3", "4"})
EXPECTED_AGE_BINS = ("20-29", "30-39", "40-49", "50-59", "60-69", "70-79")
AGE_MID = {
    "20-29": 24.5,
    "30-39": 34.5,
    "40-49": 44.5,
    "50-59": 54.5,
    "60-69": 64.5,
    "70-79": 74.5,
}

# Frozen marker lists. Neuronal set is the prompt list. Glial set frozen a priori
# (astrocyte / oligodendrocyte / microglia), listed before any score is computed.
NEURONAL_MARKERS = ("RBFOX3", "SNAP25", "SYT1", "GAD1", "SLC17A7")
GLIAL_MARKERS = ("GFAP", "AQP4", "MBP", "PLP1", "OLIG2", "CX3CR1", "AIF1")

DLPFC_GENES_CSV = DATA_PROC / "brain_aging_phase1_logcpm_genes.csv"
MANIFEST_PATH = GTEX_DIR / "manifest.json"
PROGRESS_PATH = ROOT / "PROGRESS_GTEX.md"
FINDINGS_PATH = ROOT / "FINDINGS_GTEX.md"
POSCTRL_FINDINGS = ROOT / "FINDINGS_POSCTRL.md"
PREREG_FLAG = GTEX_DIR / "PREREG_20260917.flag"
PC_DIR = RESULTS / "posctrl"

CAVEATS = """1. Public `AGE` is a 10-year bin. Use bin midpoints; report Spearman with the bin as an ordinal alongside R². Year-level MAE is
   not meaningful. Exact ages need dbGaP (`phs000424.v10`) — out of scope here.
2. Bulk tissue. A stable bulk age direction can be compositional (the blood lesson). Composition is handled explicitly below.
3. Per the GTEx FAQ, `Brain - Cortex` was sampled at the collection site and preserved in PAXgene; `Brain - Frontal Cortex (BA9)`
   and the other brain sub-regions were sampled later at the Miami Brain Endowment Bank and snap-frozen (longer ischemic time).
   So `SMCENTER` is the collection site, not the dissection site, for BA9. The BA9/Cortex pair on the same donors is a
   same-region, same-donor, different-preservation contrast — use it as such."""

PREREG_BARS_ORIGINAL = """primary bar is **transfer R² > 0 with permutation null ≤ 0.05**; bootstrap median
pairwise angle is **reported against the POSCTRL reference values** (sex control 41.3°, and the planted dense_shared angle at the
nearest n from the C3 power curve), not against 35°. The 35° bar is retired."""

PREREG_READING_ORIGINAL = """- BA9 raw passes transfer AND resid does not change it materially (sign unchanged, Δangle < 5°) → a stable bulk cortex age direction
  exists and is not driven by recorded logistics; the DLPFC instability is then attributable to something GTEx does not share:
  two-bank structure, snRNA/pseudobulk noise, or unrecorded covariates. **Say all three; do not pick one.**
- BA9 raw fails, resid passes → recorded logistics (RIN / ischemic time / Hardy / composition) were hiding a stable direction. This
  makes unrecorded PMI/RIN a plausible, not proven, cause of the DLPFC failure.
- Both fail at BA9 n ≥ the n where POSCTRL's planted direction transferred → cortex age directions fail transfer in a second,
  differently-structured cohort. The negative generalizes. Say so.
- Both fail at n below that → uninformative; report as such.
- T3 fails while T2 passes → the preservation pathway defeats transfer even when collection center does not. Report as a distinct
  finding."""

PREREG_BARS = """Written **before any Stage 1 fit**, 2026-09-17. No threshold in this block is re-tuned after numbers exist.

FINDINGS_POSCTRL.md Supersession 2 fired outcome (i): the DLPFC age direction transfers between banks as a correlation
(gene_ridge r +0.568 / +0.578, nulls ≈ 0); uncalibrated transfer R² is a calibration artifact. Stage 1 gate is open.

Primary transfer metric is held-out Pearson r between predicted and true age-bin midpoint, with a donor-level permutation
null within `SMCENTER` (`n_perm=200`). **Pass = r > 0 with null ≤ 0.05.** Also report uncalibrated R², calibrated R²
(intercept and slope refit on the test set), and Spearman ρ. Bootstrap median pairwise angle is **reported** against the
POSCTRL reference values (sex control 41.3°, planted dense_shared at C3 n=233: 41.191° — nearest n to BA9 n=269), not
against 35°. The 35° bar is retired. **Angle is not a pass/fail gate.**

POSCTRL planted `dense_shared` transferred as r at n=150 (10/10 draws) and n=233 (1/1). BA9 n=269 is at or above that n.

The original Prompt B bar (transfer R²) was never evaluated on GTEx numbers; it is superseded by this block before any
T-cell is fit."""

PREREG_READING = """- BA9 raw passes transfer (held-out Pearson r > 0 with permutation null ≤ 0.05) AND resid does not change it materially
  (sign of r unchanged, Δangle < 5°) → a stable bulk cortex age direction exists and is not driven by recorded logistics;
  the DLPFC instability is then attributable to something GTEx does not share: two-bank structure, snRNA/pseudobulk noise,
  or unrecorded covariates. **Say all three; do not pick one.**
- BA9 raw fails, resid passes → recorded logistics (RIN / ischemic time / Hardy / composition) were hiding a stable
  direction. This makes unrecorded PMI/RIN a plausible, not proven, cause of the DLPFC failure.
- Both fail at BA9 n ≥ the n where POSCTRL's planted direction transferred as r → cortex age directions fail transfer in a
  second, differently-structured cohort. The negative generalizes. Say so.
- Both fail at n below that → uninformative; report as such.
- T3 fails while T2 passes → the preservation pathway defeats transfer even when collection center does not. Report as a
  distinct finding."""

LOWDIM_L4_QUOTE = (
    "**The negative result stands.** What would be needed to change it: more donors per "
    "cell type (so that even gene space is not p ≫ n, or so that a k ≪ n space is estimated "
    "from hundreds of independent brains rather than ~200); **more sites** (two banks are one "
    "degree of freedom of transfer, and that transfer failed); or a **different tissue** with "
    "less bank structure and more independent donors. Re-fitting more unsupervised variants on "
    "this same 233-donor DLPFC matrix will not make an underdetermined direction unique."
)


class StopStep(Exception):
    """Fail loudly: missing/renamed column, unexpected levels, failed assertion."""

    def __init__(self, step: str, message: str, details=None):
        super().__init__(f"[{step}] {message}")
        self.step = step
        self.message = message
        self.details = details or {}


def gtex_log_banner(log, stage: str):
    log("=" * 100)
    log(f"GTEX {stage}  seed={GTEX_SEED}  release={RELEASE}  pipeline={RNASEQ_PIPELINE}  {DBGAP}")
    log("=" * 100)
    log("NO perturbation data. No TF Atlas. No candidate interventions.")
    log("FALSIFICATION.md and prior FINDINGS*.md are not modified.")
    log(f"Six tissues only: {list(TISSUES)}")
    log(f"Neuronal markers (frozen): {list(NEURONAL_MARKERS)}")
    log(f"Glial markers (frozen a priori): {list(GLIAL_MARKERS)}")
    log("z-score per gene on training rows only inside every fold.")
    log("Primary transfer metric: held-out Pearson r (n_perm=200); R²/calR²/ρ reported, not gates.")
    log("Within-batch R² stratified by grouping variable, never pooled.")
    log("Fail loudly: no substitute columns, no silent coerce, no dropping rows to pass an assertion.")


def file_md5(path: Path, chunk: int = 1 << 20) -> str:
    h = hashlib.md5()
    with open(path, "rb") as f:
        while True:
            b = f.read(chunk)
            if not b:
                break
            h.update(b)
    return h.hexdigest()


def gcs_md5_hex(b64: str) -> str:
    import base64
    return base64.b64decode(b64).hex()


def sampid_to_subj(sampid: str) -> str:
    parts = str(sampid).split("-")
    if len(parts) < 2 or parts[0] != "GTEX":
        raise StopStep("donor_id", f"SAMPID does not match GTEX-<donor>-... : {sampid!r}")
    return f"{parts[0]}-{parts[1]}"


def require_columns(df: pd.DataFrame, required, *, file_tag: str, step: str):
    have = list(df.columns)
    missing = [c for c in required if c not in df.columns]
    if missing:
        raise StopStep(
            step,
            f"{file_tag}: required column(s) missing or renamed: {missing}. "
            f"columns actually read={have}",
            dict(file=file_tag, missing=missing, columns_read=have),
        )
    return have


def require_levels(series: pd.Series, expected: frozenset, *, field: str, step: str, allow_na=True):
    vals = series.astype("string")
    if allow_na:
        present = set(v for v in vals.dropna().unique() if str(v) != "<NA>")
    else:
        present = set(str(v) for v in vals.unique())
    extra = sorted(present - set(expected))
    missing = sorted(set(expected) - present)
    if extra:
        raise StopStep(
            step,
            f"{field}: unexpected levels {extra}. observed={sorted(present)} expected={sorted(expected)}",
            dict(field=field, observed=sorted(present), expected=sorted(expected), extra=extra),
        )
    return dict(field=field, observed=sorted(present), expected=sorted(expected),
                missing_expected_not_in_data=missing, n_na=int(vals.isna().sum()))


def age_midpoint(age_bin: str) -> float:
    s = str(age_bin)
    if s not in AGE_MID:
        raise StopStep("age_bin", f"AGE bin {s!r} has no pre-declared midpoint. "
                       f"declared={list(AGE_MID)}")
    return float(AGE_MID[s])


def tmm_norm_factors(counts: np.ndarray, logratio_trim=TMM_LOGRATIO_TRIM,
                     sum_trim=TMM_SUM_TRIM, acutoff=-1e10) -> np.ndarray:
    """edgeR-style TMM factors. counts: genes × samples. Geom-mean scaled to 1."""
    counts = np.asarray(counts, np.float64)
    if counts.ndim != 2:
        raise ValueError("counts must be genes x samples")
    n_genes, n_samp = counts.shape
    lib = counts.sum(0)
    lib = np.where(lib <= 0, np.nan, lib)
    q75 = np.quantile(counts, 0.75, axis=0)
    f75 = q75 / lib
    finite = np.isfinite(f75)
    if not finite.any():
        return np.ones(n_samp)
    ref = int(np.nanargmin(np.abs(f75 - np.nanmean(f75))))
    factors = np.ones(n_samp, dtype=np.float64)
    ref_c = counts[:, ref]
    n_r = float(lib[ref])
    for j in range(n_samp):
        if j == ref or not np.isfinite(lib[j]):
            continue
        obs = counts[:, j]
        n_o = float(lib[j])
        with np.errstate(divide="ignore", invalid="ignore"):
            log_r = np.log2((obs / n_o) / (ref_c / n_r))
            abs_e = 0.5 * (np.log2(obs / n_o) + np.log2(ref_c / n_r))
            var = (n_o - obs) / n_o / obs + (n_r - ref_c) / n_r / ref_c
        keep = np.isfinite(log_r) & np.isfinite(abs_e) & (abs_e > acutoff) & np.isfinite(var) & (var > 0)
        n = int(keep.sum())
        if n < 10:
            factors[j] = 1.0
            continue
        lr, ae, v = log_r[keep], abs_e[keep], var[keep]
        r_m = pd.Series(lr).rank(method="average").to_numpy()
        r_a = pd.Series(ae).rank(method="average").to_numpy()
        lo_l = n * logratio_trim + 1.0
        hi_l = n + 1.0 - lo_l
        lo_s = n * sum_trim + 1.0
        hi_s = n + 1.0 - lo_s
        trim = (r_m >= lo_l) & (r_m <= hi_l) & (r_a >= lo_s) & (r_a <= hi_s)
        if int(trim.sum()) < 2:
            factors[j] = 1.0
            continue
        w = 1.0 / v[trim]
        factors[j] = float(2 ** (np.sum(lr[trim] * w) / np.sum(w)))
    pos = factors[np.isfinite(factors) & (factors > 0)]
    if pos.size:
        factors = factors / float(np.exp(np.mean(np.log(pos))))
    factors[~np.isfinite(factors)] = 1.0
    return factors


def log2_cpm_edger(counts: np.ndarray, norm_factors: np.ndarray, prior=EDGE_R_PRIOR) -> np.ndarray:
    """edgeR cpm(..., log=TRUE, prior.count=prior). counts genes × samples."""
    counts = np.asarray(counts, np.float64)
    nf = np.asarray(norm_factors, np.float64).ravel()
    lib = counts.sum(0) * nf
    ave = float(np.mean(lib)) if np.isfinite(lib).any() else 1.0
    if ave <= 0:
        ave = 1.0
    adj_prior = prior * lib / ave
    adj = counts + adj_prior
    den = lib + 2.0 * adj_prior
    den = np.where(den <= 0, 1.0, den)
    cpm = adj / den * 1e6
    return np.log2(np.clip(cpm, 1e-12, None))


def gene_filter_mask(counts: np.ndarray, min_count=MIN_COUNT, min_frac=MIN_FRAC) -> np.ndarray:
    """counts genes × samples. Keep genes with ≥ min_count in ≥ min_frac of samples."""
    counts = np.asarray(counts)
    n = counts.shape[1]
    if n == 0:
        return np.zeros(counts.shape[0], dtype=bool)
    return (counts >= min_count).sum(1) >= (min_frac * n)


def ols_r2(y, X, add_intercept=True):
    """Descriptive OLS R². Drops incomplete rows for this regression only; reports n."""
    y = np.asarray(y, np.float64).ravel()
    X = np.asarray(X, np.float64)
    if X.ndim == 1:
        X = X[:, None]
    if add_intercept:
        X = np.column_stack([np.ones(len(y)), X])
    ok = np.isfinite(y) & np.isfinite(X).all(1)
    n = int(ok.sum())
    p = int(X.shape[1])
    if n < p + 2 or float(np.std(y[ok])) < 1e-12:
        return dict(r2=np.nan, n=n, n_dropped=int((~ok).sum()), p=p, ok=False)
    yy, XX = y[ok], X[ok]
    beta, *_ = np.linalg.lstsq(XX, yy, rcond=None)
    pred = XX @ beta
    ss_res = float(np.sum((yy - pred) ** 2))
    ss_tot = float(np.sum((yy - yy.mean()) ** 2))
    r2 = 1.0 - ss_res / ss_tot if ss_tot > 0 else np.nan
    return dict(r2=float(r2) if np.isfinite(r2) else np.nan, n=n,
                n_dropped=int((~ok).sum()), p=p, ok=True)


def onehot_train(values, *, drop_first=False):
    """One-hot with explicit level list. Returns matrix, level names."""
    s = pd.Series(values).astype("string")
    levels = [lv for lv in s.dropna().unique().tolist() if str(lv) != "<NA>"]
    levels = sorted(str(lv) for lv in levels)
    if drop_first and len(levels) > 1:
        levels = levels[1:]
    n = len(s)
    X = np.zeros((n, len(levels)), dtype=np.float64)
    for j, lv in enumerate(levels):
        X[:, j] = (s.astype(str).to_numpy() == lv).astype(np.float64)
        X[s.isna().to_numpy(), j] = np.nan
    return X, levels


def apply_onehot(values, levels):
    s = pd.Series(values).astype("string")
    n = len(s)
    X = np.zeros((n, len(levels)), dtype=np.float64)
    raw = s.astype(str).to_numpy()
    for j, lv in enumerate(levels):
        X[:, j] = (raw == lv).astype(np.float64)
        X[s.isna().to_numpy(), j] = np.nan
    unseen = sorted(set(raw[~s.isna().to_numpy()]) - set(levels) - {"<NA>"})
    return X, unseen


def residualize_train_test(Xtr, Xte, Ztr, Zte):
    """OLS residualize each column of X on Z (with intercept), train-only fit."""
    Xtr = np.asarray(Xtr, np.float64)
    Xte = np.asarray(Xte, np.float64)
    Ztr = np.asarray(Ztr, np.float64)
    Zte = np.asarray(Zte, np.float64)
    if Ztr.ndim == 1:
        Ztr = Ztr[:, None]
        Zte = Zte[:, None]
    Ztr_i = np.column_stack([np.ones(len(Ztr)), Ztr])
    Zte_i = np.column_stack([np.ones(len(Zte)), Zte])
    ok = np.isfinite(Ztr_i).all(1)
    if int(ok.sum()) < Ztr_i.shape[1] + 2:
        return Xtr.copy(), Xte.copy(), dict(ok=False, n_train=int(ok.sum()))
    B, *_ = np.linalg.lstsq(Ztr_i[ok], Xtr[ok], rcond=None)
    return Xtr - Ztr_i @ B, Xte - Zte_i @ B, dict(ok=True, n_train=int(ok.sum()), p_z=int(Ztr_i.shape[1]))


def spearman_safe(a, b):
    from scipy import stats
    a = np.asarray(a, float)
    b = np.asarray(b, float)
    m = np.isfinite(a) & np.isfinite(b)
    if int(m.sum()) < 8:
        return np.nan
    a, b = a[m], b[m]
    if float(np.std(a)) < 1e-12 or float(np.std(b)) < 1e-12:
        return np.nan
    r, _ = stats.spearmanr(a, b)
    return float(r) if np.isfinite(r) else np.nan


def pearson_safe(a, b):
    from scipy import stats
    a = np.asarray(a, float)
    b = np.asarray(b, float)
    m = np.isfinite(a) & np.isfinite(b)
    if int(m.sum()) < 8:
        return np.nan
    a, b = a[m], b[m]
    if float(np.std(a)) < 1e-12 or float(np.std(b)) < 1e-12:
        return np.nan
    r, _ = stats.pearsonr(a, b)
    return float(r) if np.isfinite(r) else np.nan


def calibrated_r2(y, pred):
    """R² after OLS y ~ 1 + pred on the same points (test-set intercept and slope)."""
    from brain_phase1_common import r2_mae
    y = np.asarray(y, float)
    pred = np.asarray(pred, float)
    m = np.isfinite(y) & np.isfinite(pred)
    if int(m.sum()) < 3:
        return np.nan
    pcal, _, _ = calibrate_1d(pred[m], y[m], pred[m])
    return r2_mae(y[m], pcal)["r2"]


def pred_scores(y, pred):
    """Uncalibrated R², Pearson r, Spearman ρ, calibrated R². Primary transfer metric is r."""
    from brain_phase1_common import r2_mae
    sc = r2_mae(y, pred)
    return dict(
        r2=sc["r2"],
        r=pearson_safe(y, pred),
        rho=spearman_safe(y, pred),
        cal_r2=calibrated_r2(y, pred),
        n=sc["n"],
    )


def pass_r_bar(r, r_null, bar=TRANSFER_NULL_BAR):
    """Pass = r > 0 with permutation-null r ≤ 0.05."""
    try:
        rv, nv = float(r), float(r_null)
    except (TypeError, ValueError):
        return False
    return bool(np.isfinite(rv) and rv > 0 and np.isfinite(nv) and nv <= float(bar))


def dist_summary(x, numeric=True):
    x = pd.Series(x)
    if not numeric:
        vc = x.astype("string").value_counts(dropna=False).to_dict()
        return dict(n=int(len(x)), n_missing=int(x.isna().sum()), n_levels=int(x.nunique(dropna=True)),
                    counts={str(k): int(v) for k, v in vc.items()})
    v = pd.to_numeric(x, errors="coerce")
    return dict(
        n=int(len(v)), n_missing=int(v.isna().sum()),
        min=float(v.min()) if v.notna().any() else np.nan,
        p25=float(v.quantile(0.25)) if v.notna().any() else np.nan,
        median=float(v.median()) if v.notna().any() else np.nan,
        p75=float(v.quantile(0.75)) if v.notna().any() else np.nan,
        max=float(v.max()) if v.notna().any() else np.nan,
        mean=float(v.mean()) if v.notna().any() else np.nan,
    )


def load_dlpfc_gene_universe():
    if not DLPFC_GENES_CSV.exists():
        raise StopStep(
            "dlpfc_universe",
            f"DLPFC gene universe missing: {DLPFC_GENES_CSV}",
            dict(path=str(DLPFC_GENES_CSV)),
        )
    genes = pd.read_csv(DLPFC_GENES_CSV)
    if "gene_id" not in genes.columns:
        raise StopStep(
            "dlpfc_universe",
            f"{DLPFC_GENES_CSV.name} has no gene_id column. columns={list(genes.columns)}",
            dict(columns=list(genes.columns)),
        )
    ids = [strip_ensembl(g) for g in genes.gene_id.astype(str)]
    return set(ids), genes


def stage1_gate():
    """Stage 1 runs if POSCTRL S2 planted r transfers (outcomes i or ii).

    S2 2026-09-17 restated the transfer bar as held-out Pearson r (uncalibrated R² is
    a calibration artifact). Outcome (i) = planted r>0 and real-age r>0 both ways;
    (ii) planted r>0, real age r≤0; (iii) planted r≤0 → gate closed.

    If S2 files are missing, fall back to the 2026-09-16 R² transfer bar (P2/Mixed).
    """
    out = dict(open=False, reason="", posctrl_exists=POSCTRL_FINDINGS.exists(),
               c1_ran=False, status=None, source="", transfer=None, s2=None)
    if not POSCTRL_FINDINGS.exists():
        out["reason"] = "FINDINGS_POSCTRL.md does not exist"
        return out
    planted_s2 = PC_DIR / "s2_c1_dense_shared_r030.json"
    real_s2 = PC_DIR / "s2_gene_ridge.json"
    if planted_s2.exists() and real_s2.exists():
        planted = load_json(planted_s2)
        real = load_json(real_s2)
        out["c1_ran"] = True
        out["source"] = f"{planted_s2.name}+{real_s2.name}"
        def _side(rec, side):
            return dict(
                r=rec.get(f"r_{side}"),
                r_null=rec.get(f"r_{side}_null"),
                r2=rec.get(f"r2_{side}"),
                cal_r2=rec.get(f"cal_r2_{side}"),
                rho=rec.get(f"rho_{side}"),
            )
        s2 = dict(
            planted_pass=bool(planted.get("pass_r")),
            real_pass=bool(real.get("pass_r")),
            planted=dict(H_to_M=_side(planted, "H_to_M"), M_to_H=_side(planted, "M_to_H"),
                         pass_r=planted.get("pass_r")),
            real_age=dict(H_to_M=_side(real, "H_to_M"), M_to_H=_side(real, "M_to_H"),
                          pass_r=real.get("pass_r")),
        )
        out["s2"] = s2
        out["transfer"] = dict(
            H_to_M=real.get("r_H_to_M"), M_to_H=real.get("r_M_to_H"),
            H_to_M_null=real.get("r_H_to_M_null"), M_to_H_null=real.get("r_M_to_H_null"),
            pass_transfer=real.get("pass_r"),
            planted_H_to_M=planted.get("r_H_to_M"), planted_M_to_H=planted.get("r_M_to_H"),
            planted_H_to_M_null=planted.get("r_H_to_M_null"),
            planted_M_to_H_null=planted.get("r_M_to_H_null"),
            planted_pass_r=planted.get("pass_r"),
            metric="pearson_r",
            uncal_R2_H_to_M=real.get("r2_H_to_M"), uncal_R2_M_to_H=real.get("r2_M_to_H"),
        )
        if s2["planted_pass"] and s2["real_pass"]:
            out["status"] = "S2_i"
            out["open"] = True
            out["reason"] = (
                "gate open: S2 outcome (i) — planted and real-age r transfer both ways "
                f"(gene_ridge r HBCC→MSSM {float(real.get('r_H_to_M')):+.3f} / "
                f"MSSM→HBCC {float(real.get('r_M_to_H')):+.3f}; uncalibrated R² is a calibration artifact)"
            )
            return out
        if s2["planted_pass"] and not s2["real_pass"]:
            out["status"] = "S2_ii"
            out["open"] = True
            out["reason"] = "gate open: S2 outcome (ii) — planted r transfers; real age does not"
            return out
        out["status"] = "S2_iii"
        out["reason"] = (
            "S2 outcome (iii): planted dense_shared did not transfer as r. "
            "GTEx n would not be interpretable."
        )
        return out

    rec, src = _load_c1_dense_shared(PC_DIR)
    out["source"] = src
    if rec is None:
        out["c1_ran"] = False
        out["status"] = "C1_not_run"
        out["reason"] = "C1 dense_shared @ R²≈0.30 has not run (no c1_dense_shared_r030.json) and S2 files missing"
        return out
    out["c1_ran"] = True
    h = rec.get("transfer_H_to_M")
    m = rec.get("transfer_M_to_H")
    if h is None or m is None:
        out["status"] = "C1_incomplete"
        out["reason"] = f"{src} missing transfer_H_to_M / transfer_M_to_H"
        return out
    h, m = float(h), float(m)
    out["transfer"] = dict(
        H_to_M=h, M_to_H=m,
        H_to_M_null=rec.get("transfer_H_to_M_null"),
        M_to_H_null=rec.get("transfer_M_to_H_null"),
        pass_transfer=rec.get("pass_transfer"),
        achieved_within_site_r2=rec.get("achieved_within_site_r2"),
        target_r2=rec.get("target_r2"),
        boot_median_angle=rec.get("boot_median_angle"),
        metric="r2_legacy",
    )
    if h > 0 and m > 0:
        status = "P2"
    elif h > 0 or m > 0:
        status = "Mixed"
    else:
        status = "P3"
    out["status"] = status
    if status == "P3":
        out["reason"] = (
            f"P3: C1 dense_shared @ R²≈0.30 did not transfer "
            f"(HBCC→MSSM {h:+.3f}, MSSM→HBCC {m:+.3f}). "
            "GTEx n would not be interpretable."
        )
        return out
    out["open"] = True
    out["reason"] = f"gate open: post-supersession status={status} (legacy R² bar; S2 files missing)"
    return out


def _load_c1_dense_shared(pc_dir: Path):
    path = pc_dir / "c1_dense_shared_r030.json"
    if path.exists():
        return load_json(path), path.name
    for p in sorted(pc_dir.glob("c1_dense_shared*.json")):
        if p.name.endswith("_battery.json") or p.name.endswith("_bootstrap.json"):
            continue
        rec = load_json(p)
        tr = rec.get("target_r2")
        if tr is not None and abs(float(tr) - 0.30) < 0.04:
            return rec, p.name
    return None, None


def assert_disjoint(ids_tr, ids_te, *, what: str, where: str):
    overlap = set(map(str, ids_tr)) & set(map(str, ids_te))
    if overlap:
        raise StopStep(
            "held_out_assert",
            f"{where}: {what} overlap between train and test ({len(overlap)} ids). "
            "Not repairing the fold.",
            dict(where=where, what=what, n_overlap=len(overlap),
                 example=sorted(overlap)[:10]),
        )


def grouped_kfold(groups, n_splits=N_FOLDS, min_test=1):
    """Leave-group-out style k-fold. Groups with a single level → cannot split; caller must STOP."""
    groups = np.asarray(groups).astype(str)
    levels = pd.unique(groups)
    if len(levels) < 2:
        raise StopStep(
            "batch_held_out",
            f"grouped k-fold needs ≥2 group levels, got {len(levels)} {levels.tolist()}",
            dict(n_levels=int(len(levels)), levels=levels.tolist()),
        )
    n_splits = int(min(n_splits, len(levels)))
    order = np.array(levels)
    folds = []
    for k in range(n_splits):
        te_lv = set(order[k::n_splits].tolist())
        te = np.flatnonzero(np.isin(groups, list(te_lv)))
        tr = np.flatnonzero(~np.isin(groups, list(te_lv)))
        if len(te) < min_test or len(tr) < 5:
            raise StopStep(
                "batch_held_out",
                f"fold {k}: train n={len(tr)} test n={len(te)} after grouping; not repairing.",
                dict(fold=k, n_train=int(len(tr)), n_test=int(len(te)), test_levels=sorted(te_lv)),
            )
        folds.append((tr, te, sorted(te_lv)))
    return folds


def write_progress(text: str):
    PROGRESS_PATH.write_text(text, encoding="utf-8")


def _fmt_num(x):
    if x is None:
        return "NA"
    try:
        v = float(x)
    except (TypeError, ValueError):
        return str(x)
    if not np.isfinite(v):
        return "NA"
    return f"{v:+.3f}" if abs(v) < 1000 else str(v)


def progress_snapshot(next_action, stop=None, extra=""):
    """Rewrite PROGRESS_GTEX.md from results/gtex on disk. Called after every T-cell."""
    gate = stage1_gate()
    stage0 = load_json(GTEX_DIR / "stage0_summary.json") if (GTEX_DIR / "stage0_summary.json").exists() else {}
    flags = stage0.get("stop_flags") or {}
    cells = []
    for name in ("t1.json", "t2.json", "t3_t4.json", "t3.json", "t4.json", "t5.json", "t6.json",
                 "t1_comp.json", "stage1_summary.json"):
        p = GTEX_DIR / name
        if p.exists():
            cells.append(name)
    lines = [
        "# PROGRESS_GTEX",
        "",
        f"**STOP status:** {stop if stop else ('Stage 1 running' if gate.get('open') else 'Stage 1 gate closed')}.",
        f"**Next action:** {next_action}",
        "",
        "## Stage gate",
        "",
        f"- `FINDINGS_POSCTRL.md` exists: {bool(gate.get('posctrl_exists'))}",
        f"- status: `{gate.get('status')}`  open={gate.get('open')}",
        f"- source: `{gate.get('source')}`",
        f"- reason: {gate.get('reason')}",
        f"- primary transfer metric: held-out Pearson r, n_perm={N_PERM_GTEX}, pass = r>0 with null ≤ {TRANSFER_NULL_BAR}",
        f"- bootstrap angle reported only (sex ref {SEX_ANGLE_REF}°, planted n=233 {PLANTED_ANGLE_N233}°)",
        f"- pre-registration {PREREG_DATE} written before any T-cell: {PREREG_FLAG.exists()}",
        "",
        "## Pre-registered bars (verbatim, 2026-09-17)",
        "",
        PREREG_BARS,
        "",
        "## Pre-registered reading (verbatim, r not R²)",
        "",
        PREREG_READING,
        "",
        "## Seeds",
        "",
        f"- main: `{GTEX_SEED}`",
        f"- T-cell RNGs: {CELL_SEEDS}",
        "",
        "## Frozen preprocessing spec (Stage 0)",
        "",
        "- Tissues (SMTSD), six only, `SMAFRZE == RNASEQ`",
        "- Gene filter: ≥6 counts in ≥20% of samples per tissue",
        "- TMM (edgeR defaults) → log2-CPM `prior.count=2`",
        "- z-score per gene on training rows only inside every fold (Stage 1)",
        "- Map to DLPFC gene universe, Ensembl version-stripped",
        "",
        "## Frozen marker lists",
        "",
        f"- Neuronal: {', '.join(NEURONAL_MARKERS)}",
        f"- Glial: {', '.join(GLIAL_MARKERS)}",
        "",
        "## Stage 0 STOP checks",
        "",
        f"- BA9 n={flags.get('ba9_n')}  t2_not_testable={flags.get('t2_not_testable')}  "
        f"n_SMCENTER_ge25={flags.get('n_SMCENTER_ge25')}",
        f"- BA9 SMCENTER: `{flags.get('center_counts')}`",
        "",
        "## Finished Stage 1 files",
        "",
    ]
    if not cells:
        lines.append("None yet (pre-registration only).")
    else:
        for name in cells:
            lines.append(f"- `results/gtex/{name}`")
    lines += ["", "## Finished cell numbers (from disk)", ""]
    # T1
    t1p = GTEX_DIR / "t1.csv"
    if t1p.exists():
        df = pd.read_csv(t1p)
        lines.append("### T1")
        keep = [c for c in ("tissue", "regime", "method", "boot_median_angle",
                            "within_batch_r_median_over_centers", "within_batch_r_null",
                            "within_batch_r2_median_over_centers", "within_batch_r2_null",
                            "pass_r") if c in df.columns]
        for _, r in df.iterrows():
            bits = [str(r.get("tissue", "")), str(r.get("regime", "")), str(r.get("method", ""))]
            if "boot_median_angle" in df.columns:
                bits.append(f"angle={_fmt_num(r.boot_median_angle)}°")
            if "within_batch_r_median_over_centers" in df.columns:
                bits.append(f"r={_fmt_num(r.within_batch_r_median_over_centers)} (null {_fmt_num(r.within_batch_r_null) if 'within_batch_r_null' in df.columns else 'NA'})")
            if "within_batch_r2_median_over_centers" in df.columns:
                bits.append(f"R²={_fmt_num(r.within_batch_r2_median_over_centers)}")
            lines.append("- " + " | ".join(bits))
        lines.append("")
    t2p = GTEX_DIR / "t2.json"
    if t2p.exists():
        t2 = load_json(t2p)
        lines.append("### T2")
        if isinstance(t2, dict) and t2.get("skipped"):
            lines.append(f"- SKIPPED: {t2.get('reason')}")
        elif isinstance(t2, list):
            for rec in t2:
                lines.append(
                    f"- {rec.get('test')} {rec.get('regime')} {rec.get('method')} "
                    f"r={_fmt_num(rec.get('r'))} (null {_fmt_num(rec.get('r_null'))}) "
                    f"R²={_fmt_num(rec.get('r2'))} calR²={_fmt_num(rec.get('cal_r2'))} "
                    f"ρ={_fmt_num(rec.get('rho'))} pass={rec.get('pass_r')}"
                )
        lines.append("")
    t34p = GTEX_DIR / "t3_t4.csv"
    if t34p.exists():
        df = pd.read_csv(t34p)
        lines.append("### T3 / T4")
        for _, r in df.iterrows():
            lines.append(
                f"- {r.get('test')} {r.get('label')} {r.get('regime')} {r.get('method')} "
                f"r={_fmt_num(r.get('r'))} (null {_fmt_num(r.get('r_null') if 'r_null' in df.columns else r.get('null_median'))}) "
                f"R²={_fmt_num(r.get('r2'))} calR²={_fmt_num(r.get('cal_r2') if 'cal_r2' in df.columns else np.nan)} "
                f"ρ={_fmt_num(r.get('rho') if 'rho' in df.columns else r.get('spearman'))} "
                f"pass={r.get('pass_r') if 'pass_r' in df.columns else 'NA'}"
            )
        lines.append("")
    if (GTEX_DIR / "t5.json").exists():
        t5 = load_json(GTEX_DIR / "t5.json")
        lines.append("### T5")
        lines.append(f"- {t5}")
        lines.append("")
    if (GTEX_DIR / "t6.json").exists():
        t6 = load_json(GTEX_DIR / "t6.json")
        lines.append("### T6")
        lines.append(f"- {t6}")
        lines.append("")
    if extra:
        lines += ["## Note", "", extra, ""]
    lines += [
        "## Files",
        "",
        "- `src/gtex_common.py`, `gtex_stage0.py`, `gtex_stage1.py`, `gtex_findings.py`",
        "- `notebooks/gtex_stage0.ipynb`, `gtex_stage1.ipynb`",
        "- `FINDINGS_GTEX.md`",
        "- `PROGRESS_GTEX.md`",
        "",
    ]
    write_progress("\n".join(lines) + "\n")
    return PROGRESS_PATH
