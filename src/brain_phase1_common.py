"""Shared constants / data / design / CV helpers for DLPFC Phase 1 (FINDINGS_BRAIN_PHASE1.md).

Namespaced under results/brain_phase1/. Does not modify prior FINDINGS*.md or FALSIFICATION.md.
Seed 20260914. Cohort is the R1 adult-control DLPFC set (233 donors, 3,323 pseudobulks).

DECLARED BEFORE ANY SCORE (not selected on downstream R²):
  A-genes: unique_age >= Q90(unique_age) AND unique_ct <= Q50(unique_ct)
  I-genes: unique_ct  >= Q90(unique_ct)  AND unique_age <= Q50(unique_age)
  MIXED:   unique_age >= Q90(unique_age) AND unique_ct >= Q90(unique_ct)  — excluded from both scores

P0 covariate inclusion (declared before looking at expression models):
  a donor-level, expression-affecting variable enters the P1 design iff its within-site
  Spearman |r| with age is >= 0.20 in either site or in the site-residualized pool.
  Depth (log mean UMI/nucleus, log mean genes/nucleus) is ALWAYS in the design.
"""
from __future__ import annotations

import json
import sys
from pathlib import Path

import numpy as np
import pandas as pd

if not hasattr(np, "unicode_"):
    np.unicode_ = np.str_  # type: ignore[attr-defined]

sys.path.insert(0, str(Path(__file__).resolve().parent))
from brain_common import BRAIN_DIR, BRAIN_SEED, DLPFC_NAME, ADULT_MIN_AGE, Logger  # noqa: E402
from config import DATA_PROC, RESULTS, ROOT  # noqa: E402
from rerun_common import onehot, _q, r2_matrix, residualize, bh_q  # noqa: E402
from decomp import normalize_logcpm, classify  # noqa: E402
from tissue_download_top import resolve, local_path  # noqa: E402

PHASE1_DIR = RESULTS / "brain_phase1"
PHASE1_FIG = PHASE1_DIR / "figures"
PHASE1_TAB = PHASE1_DIR
for p in (PHASE1_DIR, PHASE1_FIG):
    p.mkdir(parents=True, exist_ok=True)

PHASE1_SEED = BRAIN_SEED  # 20260914
N_SITE_FOLDS = 5
N_INNER = 3
N_BOOT = 1000
N_MATCH_DRAWS = 20
N_PERM_P1 = 200
STOP_P2_MIN_A = 30
STOP_P3_AGE_R2 = 0.15
STOP_P3_UNREST_R2 = 0.35
IDENTITY_NEAR_CHANCE_MULT = 2.0
ALPHAS = np.array([0.01, 0.03, 0.1, 0.3, 1.0, 3.0, 10.0, 30.0, 100.0, 300.0, 1000.0, 3000.0, 10000.0])
LOGREG_CS = np.array([0.01, 0.1, 1.0, 10.0])

# Quantile LEVELS declared a priori (the numeric cutoffs are the corresponding
# percentiles of the corrected unique_age / unique_ct distribution, written to
# FINDINGS_BRAIN_PHASE1.md after P1 and before any score is built).
PRIMARY_QUANTILES = dict(t_age_q=0.90, k_ct_q=0.50, t_ct_q=0.90, k_age_q=0.50)
SWEEP_HIGH_Q = (0.80, 0.90, 0.95, 0.99)
SWEEP_LOW_Q = (0.25, 0.50, 0.75)

# P0: within-site |Spearman r| with age required to add a covariate.
P0_AGE_CORR_ABS = 0.20
NEURONAL_CLASSES = ("EN", "IN")
LOW_N_GENES_CUTOFF = 500  # snRNA debris / empty-ish nuclei proxy

# Blood within-pool corrected unique_age (FINDINGS_RERUN.md). Do not recompute.
BLOOD_UA_MAX = 0.0126
BLOOD_UA_P99 = 0.0019
BLOOD_UA_P50 = 0.00011
BLOOD_UCT_P50 = 0.1248
BLOOD_UCT_MAX = 0.9610
BLOOD_N_A_SURVIVING = 2

HASH_POOL_LIMITATION = (
    "LIMITATION: the source paper describes 6-sample hashing pools; they are not in obs "
    "and cannot be corrected for. Source (HBCC/MSSM) is a 2-level brain-bank factor (1 df), "
    "not equivalent to OneK1K's 75 multiplexed 10x pools. Barcodekey looks like "
    "Donor-{n}-{10x_barcode}-{gem_group} but is not used as a pool ID."
)

OBS_ABSENT_KEYWORDS = (
    "pmi", "postmortem", "post-mortem", "post_mortem", "rin", "rna_integrity",
    "tissue_ph", "ph_", "cause_of_death", "causeofdeath", "manner_of_death",
    "cod", "soup", "ambient", "decontx", "soupx", "hash_pool", "hashing",
    "hashtag", "multiplex_pool", "tissue_handling",
)

LOGCPM_NPY = DATA_PROC / "brain_aging_phase1_logcpm.npy"
LOGCPM_GENES = DATA_PROC / "brain_aging_phase1_logcpm_genes.csv"
LOGCPM_OBS = DATA_PROC / "brain_aging_phase1_obs.csv"
PB_PATH = DATA_PROC / "tissue_brain_aging_pseudobulk.h5ad"


def dump_json(path: Path, obj) -> None:
    path = Path(path)
    path.write_text(json.dumps(obj, indent=1, default=_json_default), encoding="utf-8")


def _json_default(o):
    if isinstance(o, (np.integer,)):
        return int(o)
    if isinstance(o, (np.floating,)):
        return float(o) if np.isfinite(o) else None
    if isinstance(o, np.ndarray):
        return o.tolist()
    if isinstance(o, Path):
        return str(o)
    raise TypeError(type(o))


def load_json(path: Path):
    return json.loads(Path(path).read_text(encoding="utf-8"))


def ensure_logcpm(log=print):
    """Load (or rebuild) the 3,323-pseudobulk log-CPM matrix used in R3 / Phase 1."""
    import anndata as ad

    if LOGCPM_NPY.exists() and LOGCPM_GENES.exists() and LOGCPM_OBS.exists():
        Y = np.load(LOGCPM_NPY).astype(np.float64)
        genes = pd.read_csv(LOGCPM_GENES)
        obs = pd.read_csv(LOGCPM_OBS)
        if Y.shape == (len(obs), len(genes)):
            _cast_obs(obs)
            log(f"[data] cached log-CPM {Y.shape[0]:,} x {Y.shape[1]:,} from {LOGCPM_NPY.name}")
            return Y, genes, obs
        log("[data] cache shape mismatch — rebuilding log-CPM")

    if not PB_PATH.exists():
        raise FileNotFoundError(
            f"{PB_PATH} missing — run src/brain_analyze.py first (R1–R3 pseudobulk)"
        )
    log(f"[data] reading {PB_PATH}")
    pb = ad.read_h5ad(PB_PATH)
    Y, gid, sym = normalize_logcpm(pb, log=log)
    genes = pd.DataFrame({"gene_id": gid, "symbol": np.asarray(sym).astype(str)})
    obs = pb.obs.copy().reset_index(drop=True)
    _cast_obs(obs)
    np.save(LOGCPM_NPY, Y.astype(np.float32))
    genes.to_csv(LOGCPM_GENES, index=False)
    obs.to_csv(LOGCPM_OBS, index=False)
    log(f"[data] wrote {LOGCPM_NPY} ({Y.shape[0]:,} x {Y.shape[1]:,})")
    return Y.astype(np.float64), genes, obs


def _cast_obs(obs: pd.DataFrame) -> None:
    obs["donor"] = obs["donor"].astype(str)
    obs["celltype"] = obs["celltype"].astype(str)
    if "Source" in obs:
        obs["Source"] = obs["Source"].astype(str)
    if "sex" in obs:
        obs["sex"] = obs["sex"].astype(str)
    obs["age"] = obs["age"].astype(float)
    obs["batch"] = obs["Source"].astype(str) if "Source" in obs else obs.get("batch", "NA").astype(str)


def assert_donors_disjoint(train_obs: pd.DataFrame, test_obs: pd.DataFrame, extra_name="site"):
    d_overlap = set(train_obs.donor.astype(str)) & set(test_obs.donor.astype(str))
    assert not d_overlap, f"donor leakage: {sorted(d_overlap)[:8]}"
    if extra_name == "site" and "Source" in train_obs and "Source" in test_obs:
        # LOSO: sites must also be disjoint. Within-site k-fold: they must match.
        pass
    return True


def loso_folds(obs: pd.DataFrame):
    """Leave-one-site-out: 2 folds. Returns list of (train_idx, test_idx, test_site)."""
    sites = sorted(obs.Source.astype(str).unique())
    assert len(sites) == 2, f"expected 2 sites, got {sites}"
    folds = []
    idx = np.arange(len(obs))
    for te_site in sites:
        te = idx[obs.Source.astype(str).to_numpy() == te_site]
        tr = idx[obs.Source.astype(str).to_numpy() != te_site]
        tr_obs, te_obs = obs.iloc[tr], obs.iloc[te]
        assert_donors_disjoint(tr_obs, te_obs)
        s_overlap = set(tr_obs.Source.astype(str)) & set(te_obs.Source.astype(str))
        assert not s_overlap, f"site leakage in LOSO: {s_overlap}"
        folds.append((tr, te, te_site))
    return folds


def within_site_donor_folds(obs: pd.DataFrame, n_folds: int, rng: np.random.Generator):
    """Donor-grouped k-fold NESTED WITHIN each site. Concatenated OOF covers every donor.

    Fold ids are per-site; a training fold never contains a donor from the test fold,
    and never mixes the other site into test.
    Returns list of (train_idx, test_idx, site, fold_k).
    """
    folds = []
    idx = np.arange(len(obs))
    for site in sorted(obs.Source.astype(str).unique()):
        site_mask = obs.Source.astype(str).to_numpy() == site
        site_idx = idx[site_mask]
        donors = np.array(sorted(obs.iloc[site_idx].donor.astype(str).unique()))
        rng.shuffle(donors)
        n_f = int(min(n_folds, max(2, len(donors))))
        fold_of = {d: i % n_f for i, d in enumerate(donors)}
        don_arr = obs.iloc[site_idx].donor.astype(str).to_numpy()
        fold = np.array([fold_of[d] for d in don_arr])
        for k in range(n_f):
            te_local = np.flatnonzero(fold == k)
            tr_local = np.flatnonzero(fold != k)
            te = site_idx[te_local]
            tr = site_idx[tr_local]
            assert_donors_disjoint(obs.iloc[tr], obs.iloc[te])
            te_sites = set(obs.iloc[te].Source.astype(str))
            tr_sites = set(obs.iloc[tr].Source.astype(str))
            assert te_sites == {site} and tr_sites == {site}, (te_sites, tr_sites, site)
            folds.append((tr, te, site, k))
    return folds


def donor_grouped_inner_splits(groups, rng, n_inner=N_INNER):
    groups = np.asarray(groups)
    uniq = np.array(sorted(set(groups)))
    rng.shuffle(uniq)
    n_in = int(max(2, min(n_inner, len(uniq))))
    fold_of = {g: i % n_in for i, g in enumerate(uniq)}
    fold = np.array([fold_of[g] for g in groups])
    out = []
    for k in range(n_in):
        te = np.flatnonzero(fold == k)
        tr = np.flatnonzero(fold != k)
        if len(tr) < 4 or len(te) < 1:
            continue
        assert not (set(groups[tr]) & set(groups[te]))
        out.append((tr, te))
    return out


def build_design(obs: pd.DataFrame, extra_cols=()):
    """Corrected Phase-1 design.

    Nuisance N = C(Source) + within-type-centred log depth (2) + within-type-centred extra P0 covs.
    unique_age = R²(N+ct+age) − R²(N+ct)   (within-site, depth/P0-adjusted age variance)
    unique_ct  = R²(N+ct+age) − R²(N+age)
    """
    D_ct = onehot(obs.celltype.to_numpy())
    src = obs.Source.astype(str).to_numpy() if "Source" in obs else obs.batch.astype(str).to_numpy()
    D_src = onehot(src)
    if D_src.shape[1] >= 2:
        D_src = D_src[:, 1:]
    else:
        D_src = np.zeros((len(obs), 0))
    umi = np.log(np.clip(obs.mean_counts_per_cell.to_numpy(float), 1e-12, None))
    ngenes = np.log(np.clip(obs.mean_genes_per_cell.to_numpy(float), 1e-12, None))
    depth = np.column_stack([umi, ngenes])
    Qct = _q(D_ct)
    depth_c = depth - Qct @ (Qct.T @ depth)
    extras, names = [], []
    for c in extra_cols:
        if c not in obs.columns:
            continue
        v = obs[c].to_numpy(float)
        if not np.isfinite(v).any():
            continue
        mu = np.nanmean(v)
        v = np.where(np.isfinite(v), v, mu)
        if np.nanstd(v) < 1e-12:
            continue
        vc = v - Qct @ (Qct.T @ v[:, None])[:, 0]
        extras.append(vc[:, None])
        names.append(c)
    extra = np.hstack(extras) if extras else np.zeros((len(obs), 0))
    age = obs.age.to_numpy(float)
    age_c = (age - np.mean(age))[:, None]
    ones = np.ones((len(obs), 1))
    return dict(ct=D_ct, source=D_src, depth=depth_c, extra=extra, extra_names=names,
                age=age_c, ones=ones, age_raw=age)


def _hs(parts):
    parts = [p for p in parts if p is not None and np.asarray(p).ndim == 2 and p.shape[1] > 0]
    if not parts:
        raise ValueError("empty design")
    return np.hstack(parts)


def nuisance_block(B):
    return _hs([B["ones"], B["source"], B["depth"], B["extra"]])


def decompose_brain(Y, obs, extra_cols=(), gene_ids=None, symbols=None):
    """Per-gene commonality partition under the brain Phase-1 model."""
    Y = np.asarray(Y, dtype=np.float64)
    Yc = Y - Y.mean(0, keepdims=True)
    ss_tot = (Yc ** 2).sum(0)
    ss_tot[ss_tot == 0] = np.nan
    B = build_design(obs, extra_cols=extra_cols)
    N = nuisance_block(B)
    R2 = {}
    R2["N"] = r2_matrix(N, Yc, ss_tot)
    R2["Nct"] = r2_matrix(_hs([B["ct"], B["source"], B["depth"], B["extra"]]), Yc, ss_tot)
    R2["Nage"] = r2_matrix(_hs([N, B["age"]]), Yc, ss_tot)
    X_full = _hs([B["ct"], B["source"], B["depth"], B["extra"], B["age"]])
    R2["full"] = r2_matrix(X_full, Yc, ss_tot)
    R2["int"] = r2_matrix(_hs([X_full[:, :-1], B["ct"] * B["age"]]), Yc, ss_tot)
    R2["ct"] = r2_matrix(B["ct"], Yc, ss_tot)
    R2["age"] = r2_matrix(_hs([B["ones"], B["age"]]), Yc, ss_tot)
    R2["ct_age"] = r2_matrix(_hs([B["ct"], B["age"]]), Yc, ss_tot)
    R2["ct_src"] = r2_matrix(_hs([B["ct"], B["source"]]), Yc, ss_tot)
    R2["ct_src_age"] = r2_matrix(_hs([B["ct"], B["source"], B["age"]]), Yc, ss_tot)

    df = pd.DataFrame({f"R2_{k}": v for k, v in R2.items()})
    df["unique_ct"] = df.R2_full - df.R2_Nage
    df["unique_age"] = df.R2_full - df.R2_Nct
    df["shared"] = df.R2_Nct + df.R2_Nage - df.R2_full - df.R2_N
    df["nuisance"] = df.R2_N
    df["resid"] = 1.0 - df.R2_full
    df["interaction"] = df.R2_int - df.R2_full
    df["partial_age"] = df.unique_age / (1.0 - df.R2_Nct)
    df["unique_age_pooled"] = df.R2_ct_age - df.R2_ct
    df["unique_ct_pooled"] = df.R2_ct_age - df.R2_age
    df["unique_age_source_only"] = df.R2_ct_src_age - df.R2_ct_src
    beta, *_ = np.linalg.lstsq(X_full, Yc, rcond=None)
    df["age_slope_per_year"] = beta[-1]
    n, p = X_full.shape
    F = df.unique_age * (n - p) / df.resid.replace(0, np.nan)
    from scipy import stats
    df["age_F"] = F
    df["age_p_nominal"] = stats.f.sf(F, 1, max(n - p, 1))
    df["mean_logcpm"] = Y.mean(0)
    if gene_ids is not None:
        df.index = pd.Index(gene_ids, name="gene_id")
    if symbols is not None:
        df.insert(0, "symbol", np.asarray(symbols))
    return df


def quantiles_to_thresholds(df: pd.DataFrame, qdict=None):
    qdict = qdict or PRIMARY_QUANTILES
    ua = df["unique_age"].to_numpy(float)
    uc = df["unique_ct"].to_numpy(float)
    t = dict(
        t_age=float(np.quantile(ua, qdict["t_age_q"])),
        k_ct=float(np.quantile(uc, qdict["k_ct_q"])),
        t_ct=float(np.quantile(uc, qdict["t_ct_q"])),
        k_age=float(np.quantile(ua, qdict["k_age_q"])),
        t_age_q=qdict["t_age_q"],
        k_ct_q=qdict["k_ct_q"],
        t_ct_q=qdict["t_ct_q"],
        k_age_q=qdict["k_age_q"],
    )
    return t


def apply_thresholds(df: pd.DataFrame, thresh: dict) -> pd.Series:
    return classify(df, t_age=thresh["t_age"], k_ct=thresh["k_ct"],
                    t_ct=thresh["t_ct"], k_age=thresh["k_age"])


def r2_mae(y, pred):
    y = np.asarray(y, float)
    pred = np.asarray(pred, float)
    m = np.isfinite(y) & np.isfinite(pred)
    y, pred = y[m], pred[m]
    if len(y) < 3:
        return dict(r2=np.nan, mae=np.nan, n=int(len(y)))
    ss_tot = ((y - y.mean()) ** 2).sum()
    r2 = np.nan if ss_tot <= 0 else 1.0 - ((y - pred) ** 2).sum() / ss_tot
    return dict(r2=float(r2), mae=float(np.abs(y - pred).mean()), n=int(len(y)))


def bootstrap_r2_mae(y, pred, rng, n_boot=N_BOOT):
    y = np.asarray(y, float)
    pred = np.asarray(pred, float)
    n = len(y)
    if n < 8:
        m = r2_mae(y, pred)
        return dict(r2=m["r2"], mae=m["mae"], r2_lo=np.nan, r2_hi=np.nan, mae_lo=np.nan, mae_hi=np.nan, n=n)
    r2s, maes = np.empty(n_boot), np.empty(n_boot)
    for b in range(n_boot):
        i = rng.integers(0, n, n)
        m = r2_mae(y[i], pred[i])
        r2s[b], maes[b] = m["r2"], m["mae"]
    m = r2_mae(y, pred)
    return dict(r2=m["r2"], mae=m["mae"], n=n,
                r2_lo=float(np.nanpercentile(r2s, 2.5)), r2_hi=float(np.nanpercentile(r2s, 97.5)),
                mae_lo=float(np.nanpercentile(maes, 2.5)), mae_hi=float(np.nanpercentile(maes, 97.5)))


def strip_ensembl(gid: str) -> str:
    s = str(gid)
    return s.split(".")[0]


def load_p0_extras():
    """P0 covariates that were included in the P1 design (may be empty).

    Donor-level columns are merged as p0_<name> so they do not collide with
    per-pseudobulk fields (cc_umi_frac, n_cells, ...).
    """
    p = PHASE1_DIR / "p0_included_covariates.json"
    if not p.exists():
        return []
    rec = load_json(p)
    out = []
    for c in rec.get("included", []):
        out.append(c if str(c).startswith("p0_") else f"p0_{c}")
    return out


def merge_p0_into_obs(obs: pd.DataFrame) -> pd.DataFrame:
    """Attach donor-level P0 covariates onto pseudobulk obs as p0_* columns."""
    path = PHASE1_DIR / "p0_donor_covariates.csv"
    if not path.exists():
        return obs
    don = pd.read_csv(path)
    don["donor"] = don["donor"].astype(str)
    extra_cols = [c for c in don.columns if c not in ("donor", "Source", "age", "sex")]
    don2 = don[["donor"] + extra_cols].copy()
    don2 = don2.rename(columns={c: f"p0_{c}" for c in extra_cols})
    return obs.merge(don2, on="donor", how="left")
