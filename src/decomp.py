"""Phase 1a: per-gene variance decomposition on pseudobulk log-CPM.

Model (as specified): expression ~ cell_type + donor_age  (additive OLS, one row per donor x cell type).
Partition of total variance for each gene via commonality analysis on R^2:
    R2_ct     : R^2 of  y ~ C(cell_type)
    R2_age    : R^2 of  y ~ age
    R2_full   : R^2 of  y ~ C(cell_type) + age
    unique_ct  = R2_full - R2_age        (cell-type-explained, unique)
    unique_age = R2_full - R2_ct         (age-explained, unique)
    shared     = R2_ct + R2_age - R2_full (shared; ~0 when design is balanced)
    resid      = 1 - R2_full
Extras reported: R2_int (y ~ C(cell_type) * age) and interaction = R2_int - R2_full, i.e. the variance explained
by cell-type-specific age slopes; partial_age = unique_age / (1 - R2_ct) = fraction of within-cell-type variance
explained by age.
"""
import sys
import numpy as np
import pandas as pd
import anndata as ad

sys.path.insert(0, str(__import__("pathlib").Path(__file__).resolve().parent))
from config import TAB  # noqa


def normalize_logcpm(pb: ad.AnnData, min_frac_detected=0.10, min_mean_logcpm=None, log=print):
    X = np.asarray(pb.X, dtype=np.float64)
    lib = X.sum(1, keepdims=True)
    cpm = X / lib * 1e6
    logcpm = np.log2(cpm + 1.0)
    det = (cpm > 1.0).mean(0)
    keep = det >= min_frac_detected
    n0 = X.shape[1]
    log(f"[filter] genes with CPM>1 in >= {min_frac_detected:.0%} of pseudobulks   genes  {n0:,} -> {keep.sum():,}  (removed {n0-keep.sum():,})")
    # drop non-protein-coding? keep all biotypes but record. Drop feature_is_filtered genes (CELLxGENE flag).
    if "feature_is_filtered" in pb.var:
        ff = pb.var["feature_is_filtered"].astype(bool).to_numpy()
        n1 = keep.sum()
        keep &= ~ff
        log(f"[filter] drop CELLxGENE feature_is_filtered genes                genes  {n1:,} -> {keep.sum():,}  (removed {n1-keep.sum():,})")
    return logcpm[:, keep].astype(np.float32), pb.var.index[keep].to_numpy(), pb.var["feature_name"].astype(str).to_numpy()[keep]


def _qr_resid_ss(Q, Y):
    """Residual sum of squares of Y projected off span(Q); Q has orthonormal columns."""
    coef = Q.T @ Y
    fitted_ss = (coef ** 2).sum(0)
    return fitted_ss  # explained SS (Y assumed centered)


def design_matrices(ct, age):
    ct = pd.Categorical(ct)
    D_ct = pd.get_dummies(ct, dtype=float).to_numpy()  # n x k (includes implicit intercept via full one-hot)
    age_c = (np.asarray(age, float) - np.mean(age))[:, None]
    ones = np.ones((len(age_c), 1))
    X_ct = D_ct
    X_age = np.hstack([ones, age_c])
    X_full = np.hstack([D_ct, age_c])
    X_int = np.hstack([D_ct, D_ct * age_c])
    return X_ct, X_age, X_full, X_int


def decompose(Y, ct, age, gene_ids=None, symbols=None):
    """Y: (n_pseudobulk x n_genes) log-CPM. Returns per-gene DataFrame of variance components."""
    Y = np.asarray(Y, dtype=np.float64)
    Yc = Y - Y.mean(0, keepdims=True)
    ss_tot = (Yc ** 2).sum(0)
    ss_tot[ss_tot == 0] = np.nan
    out = {}
    for name, X in zip(("ct", "age", "full", "int"), design_matrices(ct, age)):
        Q, _ = np.linalg.qr(X - 0.0)  # QR of design (intercept is in the column space of one-hots)
        out[f"R2_{name}"] = _qr_resid_ss(Q, Yc) / ss_tot
    df = pd.DataFrame(out)
    df["unique_ct"] = df.R2_full - df.R2_age
    df["unique_age"] = df.R2_full - df.R2_ct
    df["shared"] = df.R2_ct + df.R2_age - df.R2_full
    df["resid"] = 1.0 - df.R2_full
    df["interaction"] = df.R2_int - df.R2_full
    df["partial_age"] = df.unique_age / (1.0 - df.R2_ct)
    # signed age slope (pooled additive model) for direction
    X_full = design_matrices(ct, age)[2]
    beta, *_ = np.linalg.lstsq(X_full, Yc, rcond=None)
    df["age_slope_per_year"] = beta[-1]
    df["mean_logcpm"] = Y.mean(0)
    if gene_ids is not None:
        df.index = gene_ids
    if symbols is not None:
        df.insert(0, "symbol", symbols)
    return df


def classify(df, t_age, k_ct, t_ct, k_age, age_col="unique_age", ct_col="unique_ct"):
    """A: age>=t_age & ct<=k_ct ; I: ct>=t_ct & age<=k_age ; MIXED: age>=t_age & ct>=t_ct ; else OTHER."""
    a = df[age_col].to_numpy(); c = df[ct_col].to_numpy()
    lab = np.full(len(df), "OTHER", dtype=object)
    lab[(a >= t_age) & (c <= k_ct)] = "A"
    lab[(c >= t_ct) & (a <= k_age)] = "I"
    lab[(a >= t_age) & (c >= t_ct)] = "MIXED"
    return pd.Series(lab, index=df.index, name="gene_class")


# Sweep grid. HISTORY (recorded for transparency): the first grid written before seeing any data was
# t_age in {0.01,0.02,0.03,0.05,0.08,0.12}; the observed per-gene unique-age variance turned out to have
# max 0.035 / 99th pct 0.0074, so that grid produced EMPTY A-gene sets above 0.03. The grid below is
# re-centred on the observed distribution (first run of run_1a_decomp.py, results/phase1a_decomp_initial.txt).
# This recalibration happened BEFORE any downstream (1c/1d) model was fitted and is not a STOP-condition threshold.
SWEEP = dict(
    t_age=[0.002, 0.005, 0.01, 0.02],              # min unique age variance for A-genes (~top 9%, 2.5%, 0.6%, 0.15%)
    k_ct=[0.10, 0.20, 0.30, 0.50],                 # max cell-type variance for A-genes
    t_ct=[0.30, 0.50, 0.70],                       # min cell-type variance for I-genes
    k_age=[0.0005, 0.001, 0.002, 0.005],           # max age variance for I-genes
)
# Primary configuration, chosen from the variance distribution (set sizes), before any downstream model:
PRIMARY_THRESH = dict(t_age=0.005, k_ct=0.20, t_ct=0.50, k_age=0.001)
