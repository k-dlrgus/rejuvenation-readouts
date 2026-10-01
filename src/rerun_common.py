"""Shared code for the Phase 1 RERUN (corrected within-pool model).

The prior run (FINDINGS.md) halted at STOP 1b because `expression ~ cell_type + age` selects depth-artifact genes on OneK1K
(donors not age-randomised across the 75 multiplexed 10x pools; per-cell depth is a pool property). This rerun re-derives
Phase 1 from 1a under the corrected model

    expression ~ C(cell_type) + C(pool) + log(mean UMI/cell) + log(mean genes/cell) + age

and every output is namespaced (results/rerun/, results/tables/rerun_*.csv, results/figures/rerun_*.png) so nothing from the
prior run is overwritten. Reused from the prior run: cached raw data, pseudobulk construction (src/pseudobulk.py,
src/run_pseudobulk_onek1k.py), log-CPM normalisation (decomp.normalize_logcpm), confound gene lists/flags (src/confounds.py,
src/gene_lists.py). The prior A-gene set is NOT reused.
"""
import sys
from pathlib import Path

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parent))
from config import DATA_PROC, RESULTS, FIG, TAB, SEED, PRIMARY  # noqa: E402

RERUN_DIR = RESULTS / "rerun"
RERUN_DIR.mkdir(parents=True, exist_ok=True)
RERUN_SEED = SEED  # 20260914 — same project seed; every rng in the rerun is numpy.random.default_rng(RERUN_SEED)
P = "rerun_"       # file-name prefix for tables/figures


def tab(name):
    return TAB / f"{P}{name}"


def fig(name):
    return FIG / f"{P}{name}"


class Logger:
    """print() to console and append to a text file (so the notebook and results/rerun/*.txt hold the same log)."""

    def __init__(self, path: Path, mode="w"):
        self.path = Path(path)
        self.fh = open(self.path, mode, encoding="utf-8")

    def __call__(self, *args, **kw):
        s = " ".join(str(a) for a in args)
        try:
            print(s, flush=True, **kw)
        except UnicodeEncodeError:
            enc = getattr(getattr(sys.stdout, "encoding", None), "lower", lambda: "")() or "ascii"
            print(s.encode(enc, "replace").decode(enc), flush=True, **kw)
        self.fh.write(s + "\n")
        self.fh.flush()

    def close(self):
        self.fh.close()


# --------------------------------------------------------------------------- data
def ensure_logcpm(log=print):
    """Guarantee data/processed/onek1k_logcpm.npy (+genes, +obs) exist; rebuild from the pseudobulk if missing.
    Runnable from a clean checkout: run_pseudobulk_onek1k.py builds the pseudobulk from the cached raw h5ad."""
    import subprocess
    if not (DATA_PROC / "onek1k_pseudobulk.h5ad").exists():
        log("[rerun] pseudobulk missing -> building via src/run_pseudobulk_onek1k.py (reused construction)")
        r = subprocess.run([sys.executable, str(Path(__file__).parent / "run_pseudobulk_onek1k.py")], capture_output=True,
                           text=True, encoding="utf-8", errors="replace")
        log(r.stdout)
        if r.returncode:
            raise RuntimeError(r.stderr)
    if not (DATA_PROC / "onek1k_logcpm.npy").exists():
        import anndata as ad
        from decomp import normalize_logcpm
        log("[rerun] log-CPM matrix missing -> regenerating with decomp.normalize_logcpm (reused normalisation)")
        pb = ad.read_h5ad(DATA_PROC / "onek1k_pseudobulk.h5ad")
        Y, gid, sym = normalize_logcpm(pb, log=log)
        np.save(DATA_PROC / "onek1k_logcpm.npy", Y)
        pd.DataFrame({"gene_id": gid, "symbol": sym}).to_csv(DATA_PROC / "onek1k_logcpm_genes.csv", index=False)
        pb.obs.to_csv(DATA_PROC / "onek1k_pseudobulk_obs.csv")


def load_data(log=print):
    ensure_logcpm(log)
    Y = np.load(DATA_PROC / "onek1k_logcpm.npy").astype(np.float64)
    genes = pd.read_csv(DATA_PROC / "onek1k_logcpm_genes.csv")
    obs = pd.read_csv(DATA_PROC / "onek1k_pseudobulk_obs.csv", index_col=0)
    obs["donor"] = obs["donor"].astype(str)
    obs["celltype"] = obs["celltype"].astype(str)
    obs["pool_number"] = obs["pool_number"].astype(int)
    assert Y.shape == (len(obs), len(genes)), (Y.shape, len(obs), len(genes))
    # donors nest in pools (HARD RULE 1 relies on this)
    n_pools_per_donor = obs.groupby("donor").pool_number.nunique()
    assert (n_pools_per_donor == 1).all(), "a donor appears in more than one pool"
    log(f"[data] log-CPM pseudobulks: {Y.shape[0]:,} x {Y.shape[1]:,} genes | donors={obs.donor.nunique()} | "
        f"pools={obs.pool_number.nunique()} | cell types={obs.celltype.nunique()} | every donor in exactly one pool: True")
    return Y, genes, obs


# --------------------------------------------------------------------------- design matrices
def onehot(x):
    return pd.get_dummies(pd.Categorical(x), dtype=float).to_numpy()


def _q(X):
    Q, _ = np.linalg.qr(X)
    return Q


def build_design_blocks(obs):
    """Blocks of the corrected model. Depth covariates are centred WITHIN cell type (residualised on the cell-type one-hots) so
    that they capture between-donor/between-pool depth variation only; this leaves the full-model fit unchanged (same column
    space as raw log depth + cell type) but keeps cell-type variance attributed to cell type in the partition."""
    n = len(obs)
    D_ct = onehot(obs.celltype.to_numpy())                    # n x 16, spans the intercept
    D_pool = onehot(obs.pool_number.to_numpy())[:, 1:]        # n x 74 (one level dropped; intercept in D_ct / ones)
    depth = np.column_stack([np.log(obs.mean_counts_per_cell.to_numpy(float)),
                             np.log(obs.mean_genes_per_cell.to_numpy(float))])
    Qct = _q(D_ct)
    depth_c = depth - Qct @ (Qct.T @ depth)                   # within-cell-type-centred log depth (n x 2)
    age = obs.age.to_numpy(float)
    age_c = (age - age.mean())[:, None]
    ones = np.ones((n, 1))
    return dict(ct=D_ct, pool=D_pool, depth=depth_c, age=age_c, ones=ones, age_raw=age)


def r2_matrix(X, Yc, ss_tot):
    """R^2 of every column of (centred) Yc on the column space of X."""
    Q = _q(X)
    return ((Q.T @ Yc) ** 2).sum(0) / ss_tot


def decompose_corrected(Y, obs, gene_ids=None, symbols=None):
    """Per-gene commonality partition under the corrected model, with the prior (pooled) model alongside for comparison.

    Nuisance N = C(pool) + within-type-centred log depth (2 cols).
      R2_N, R2_Nct = R2(N + ct), R2_Nage = R2(N + age), R2_full = R2(N + ct + age), R2_int = R2(N + ct + ct:age)
      unique_ct  = R2_full - R2_Nage        (cell-type variance not explained by pool/depth/age)
      unique_age = R2_full - R2_Nct         (WITHIN-POOL, depth-adjusted age variance)  <- the A-gene axis
      shared     = R2_Nct + R2_Nage - R2_full - R2_N
      nuisance   = R2_N ; resid = 1 - R2_full ; interaction = R2_int - R2_full ; partial_age = unique_age / (1 - R2_Nct)
    Prior-model columns for side-by-side comparison:
      unique_age_pooled = R2(ct+age) - R2(ct) ; unique_ct_pooled = R2(ct+age) - R2(age) ; unique_age_pool_only = R2(ct+pool+age) - R2(ct+pool)
    """
    Y = np.asarray(Y, dtype=np.float64)
    Yc = Y - Y.mean(0, keepdims=True)
    ss_tot = (Yc ** 2).sum(0)
    ss_tot[ss_tot == 0] = np.nan
    B = build_design_blocks(obs)
    hs = np.hstack
    R2 = {}
    R2["N"] = r2_matrix(hs([B["ones"], B["pool"], B["depth"]]), Yc, ss_tot)
    R2["Nct"] = r2_matrix(hs([B["ct"], B["pool"], B["depth"]]), Yc, ss_tot)
    R2["Nage"] = r2_matrix(hs([B["ones"], B["pool"], B["depth"], B["age"]]), Yc, ss_tot)
    X_full = hs([B["ct"], B["pool"], B["depth"], B["age"]])
    R2["full"] = r2_matrix(X_full, Yc, ss_tot)
    R2["int"] = r2_matrix(hs([B["ct"], B["pool"], B["depth"], B["ct"] * B["age"]]), Yc, ss_tot)
    # prior model (pooled) and pool-only (no depth) for comparison
    R2["ct"] = r2_matrix(B["ct"], Yc, ss_tot)
    R2["age"] = r2_matrix(hs([B["ones"], B["age"]]), Yc, ss_tot)
    R2["ct_age"] = r2_matrix(hs([B["ct"], B["age"]]), Yc, ss_tot)
    R2["ct_pool"] = r2_matrix(hs([B["ct"], B["pool"]]), Yc, ss_tot)
    R2["ct_pool_age"] = r2_matrix(hs([B["ct"], B["pool"], B["age"]]), Yc, ss_tot)

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
    df["unique_age_pool_only"] = df.R2_ct_pool_age - df.R2_ct_pool
    # signed age slope (log2-CPM per year) in the corrected model
    beta, *_ = np.linalg.lstsq(X_full, Yc, rcond=None)
    df["age_slope_per_year"] = beta[-1]
    beta_p, *_ = np.linalg.lstsq(hs([B["ct"], B["age"]]), Yc, rcond=None)
    df["age_slope_per_year_pooled"] = beta_p[-1]
    # nominal F test for the age term (1 df), n - p residual df. Anti-conservative (pseudobulks of one donor are not
    # independent) -> the permutation null in Stage A is the primary inference; this is reported for orientation only.
    n, p = X_full.shape
    F = df.unique_age * (n - p) / df.resid
    from scipy import stats
    df["age_F"] = F
    df["age_p_nominal"] = stats.f.sf(F, 1, n - p)
    df["mean_logcpm"] = Y.mean(0)
    if gene_ids is not None:
        df.index = pd.Index(gene_ids, name="gene_id")
    if symbols is not None:
        df.insert(0, "symbol", np.asarray(symbols))
    return df


def residualize(Y, X):
    Q = _q(X)
    return Y - Q @ (Q.T @ Y)


def bh_q(p):
    """Benjamini-Hochberg q-values for a vector of p-values."""
    p = np.asarray(p, float)
    m = len(p)
    order = np.argsort(p)
    q = np.empty(m)
    q[order] = np.minimum.accumulate((p[order] * m / np.arange(1, m + 1))[::-1])[::-1]
    return np.minimum(q, 1.0)


def permute_age_within_pool(obs, rng):
    """Permute donor ages at the DONOR level, within each pool (keeps pool age means and the donor x cell-type structure)."""
    don = obs.groupby("donor").agg(age=("age", "first"), pool=("pool_number", "first"))
    perm = don.age.copy()
    for _, idx in don.groupby("pool").groups.items():
        perm.loc[idx] = rng.permutation(don.age.loc[idx].to_numpy())
    return obs.donor.map(perm).to_numpy(float)


# --------------------------------------------------------------------------- pool-grouped CV (HARD RULE 1)
def pool_grouped_folds(obs_sub, n_folds, rng):
    """Assign POOLS (not donors) to folds; return list of (train_idx, test_idx) over rows of obs_sub, and assert that no donor
    and no pool appears in both train and test of any fold."""
    pools = np.array(sorted(obs_sub.pool_number.unique()))
    rng.shuffle(pools)
    fold_of_pool = {p: i % n_folds for i, p in enumerate(pools)}
    fold = obs_sub.pool_number.map(fold_of_pool).to_numpy()
    folds = []
    for k in range(n_folds):
        te = np.flatnonzero(fold == k)
        tr = np.flatnonzero(fold != k)
        assert_disjoint(obs_sub.iloc[tr], obs_sub.iloc[te])
        folds.append((tr, te))
    return folds


def assert_disjoint(train_obs, test_obs):
    """HARD RULE 1: no donor AND no pool in both train and test."""
    d_overlap = set(train_obs.donor) & set(test_obs.donor)
    p_overlap = set(train_obs.pool_number) & set(test_obs.pool_number)
    assert not d_overlap, f"donor leakage between train and test: {sorted(d_overlap)[:5]}"
    assert not p_overlap, f"pool leakage between train and test: {sorted(p_overlap)[:5]}"
    return True
