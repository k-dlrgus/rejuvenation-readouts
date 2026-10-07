"""SENG2 shared layer.

Reads the frozen instruments (the GTEx fibroblast ruler, the SENG split, the
SENG proliferation quartiles, the SENG axes) and rebuilds SENG's analysis panel
in SENG's own coordinate space. Nothing is refit. Writes only results/seng2/.

`.X` is streamed: with ~2.6 GB of RAM free the 148M-nnz matrix cannot be held,
so group pseudobulks are accumulated in row blocks and only the few thousand
cells that need per-cell resampling are materialised as a submatrix.
"""
from __future__ import annotations

import sys
import time
from pathlib import Path

import h5py
import numpy as np
import pandas as pd
from scipy import sparse

SRC = Path(__file__).resolve().parent
if str(SRC) not in sys.path:
    sys.path.insert(0, str(SRC))

import seng_run as S  # noqa: E402  the frozen SENG instrument, imported not modified
from fibro2_common import load_frozen_ruler  # noqa: E402
from same_run import axis_stats, ci_excludes_zero, cosine_full  # noqa: E402
from toward_run import percentile_ci, ruler_score, tmm_logcpm_quiet, zscore_frozen  # noqa: E402

ROOT = SRC.parent
SENG = ROOT / "results" / "seng"
SENG_WHY = ROOT / "results" / "seng_why"
OUT = ROOT / "results" / "seng2"
H5 = ROOT / "data" / "raw" / "seng" / "exp3_merged_not_normalized.h5ad"
MD5_EXPECT = "bfccaa94a916c22e4b0793913a04128b"

SEED = 20260914
BOOT_SEED = 20260918
N_NULL = 200
N_PERM = 200
N_BOOT = 200
MIN_CELLS = 30
P_BAR = 0.05
BLOCK = 2000
NAMED = ["CRA_E2F3", "CRA_EZH2", "CRI_STAT3", "CRI_ZFX"]
IDENTITY_PRESENT = ["COL1A1", "COL1A2", "FN1", "LUM", "PDGFRA", "PDGFRB",
                    "POSTN", "PRRX1", "SERPINH1", "VIM"]


class Log:
    def __init__(self, path):
        self.fh = open(path, "w", encoding="utf-8")
        self.t0 = time.time()

    def __call__(self, *a):
        msg = " ".join(str(x) for x in a)
        line = f"[{time.time() - self.t0:8.1f}s] {msg}"
        print(line, flush=True)
        self.fh.write(line + "\n")
        self.fh.flush()

    def close(self):
        self.fh.close()


# --------------------------------------------------------------------------
# alignment and streaming
# --------------------------------------------------------------------------
def build_expander(src, n_file_cols):
    """Join map from a file column to every local position it occupies.

    `src` may repeat a file column (SENG_WHY recorded 10 such duplicate maps),
    and `mat[:, src]` keeps each repeat as its own column. The expansion below
    reproduces that exactly instead of de-duplicating.
    """
    src = np.asarray(src, np.int64)
    order = np.argsort(src, kind="stable")
    ssort = src[order]
    left = np.searchsorted(ssort, np.arange(n_file_cols), "left")
    right = np.searchsorted(ssort, np.arange(n_file_cols), "right")
    return dict(order=order, left=left, right=right, n_local=src.size)


def expand_block(ind, dat, rnnz, exp):
    """Map one row block's (file column, value) pairs into local columns."""
    cnt = (exp["right"][ind] - exp["left"][ind]).astype(np.int64)
    keep_rows = np.add.reduceat(cnt, np.r_[0, np.cumsum(rnnz)[:-1]]) if rnnz.size else cnt[:0]
    keep_rows = np.where(rnnz > 0, keep_rows, 0)
    tot = int(cnt.sum())
    if tot == 0:
        return (np.zeros(0, np.int32), np.zeros(0, np.float64),
                np.zeros(rnnz.size + 1, np.int64))
    starts = np.repeat(exp["left"][ind], cnt)
    within = np.arange(tot, dtype=np.int64) - np.repeat(np.r_[0, np.cumsum(cnt)[:-1]], cnt)
    loc = exp["order"][starts + within].astype(np.int32)
    val = np.repeat(np.asarray(dat, np.float64), cnt)
    indptr = np.r_[0, np.cumsum(keep_rows)].astype(np.int64)
    return loc, val, indptr


def stream_pseudobulk(f, idx_list, exp, log, block=BLOCK):
    """Sum counts over each group's cells, in local (ruler-present) columns."""
    X = f["X"]
    n_cells = int(X.attrs["shape"][0])
    indptr = X["indptr"][:]
    idx_ds, dat_ds = X["indices"], X["data"]
    n_g, n_loc = len(idx_list), exp["n_local"]
    rows, cols = [], []
    for gi, idx in enumerate(idx_list):
        rows.append(np.full(len(idx), gi, np.int32))
        cols.append(np.asarray(idx, np.int32))
    Sind = sparse.csr_matrix(
        (np.ones(sum(len(c) for c in cols), np.float64),
         (np.concatenate(rows), np.concatenate(cols))),
        shape=(n_g, n_cells))
    acc = np.zeros((n_g, n_loc), np.float64)
    for a in range(0, n_cells, block):
        b = min(a + block, n_cells)
        lo, hi = int(indptr[a]), int(indptr[b])
        rnnz = np.diff(indptr[a:b + 1]).astype(np.int64)
        if hi == lo:
            continue
        loc, val, ip = expand_block(idx_ds[lo:hi], dat_ds[lo:hi], rnnz, exp)
        B = sparse.csr_matrix((val, loc, ip), shape=(b - a, n_loc))
        sub = Sind[:, a:b]
        if sub.nnz:
            acc += np.asarray((sub @ B).todense())
        if (a // block) % 5 == 0:
            log(f"[stream] rows {a}-{b}")
    return acc


def submatrix(f, rows, exp):
    """Materialise the given cells as a CSR in local columns (few thousand rows)."""
    X = f["X"]
    indptr = X["indptr"][:]
    idx_ds, dat_ds = X["indices"], X["data"]
    rows = np.asarray(rows, int)
    data_l, loc_l, ip = [], [], [0]
    for r in rows:
        lo, hi = int(indptr[r]), int(indptr[r + 1])
        if hi > lo:
            l_, v_, _ = expand_block(idx_ds[lo:hi], dat_ds[lo:hi],
                                     np.array([hi - lo], np.int64), exp)
            loc_l.append(l_)
            data_l.append(v_)
            ip.append(ip[-1] + l_.size)
        else:
            ip.append(ip[-1])
    loc = np.concatenate(loc_l) if loc_l else np.zeros(0, np.int32)
    dat = np.concatenate(data_l) if data_l else np.zeros(0, np.float64)
    return sparse.csr_matrix((dat, loc, np.asarray(ip, np.int64)),
                             shape=(rows.size, exp["n_local"]))


# --------------------------------------------------------------------------
# the panel
# --------------------------------------------------------------------------
def to_ruler(local_sum, dst, n_ruler):
    out = np.zeros(n_ruler, np.float64)
    out[dst] = np.asarray(local_sum, np.float64).ravel()
    return out


def load_context(log):
    """Everything frozen: ruler, obs, assignment, split, quartiles, axes."""
    frozen = load_frozen_ruler()
    n_ruler = int(np.asarray(frozen["mu"]).shape[0])
    f = h5py.File(H5, "r")
    obs = S.read_obs_var(f["obs"])
    var = S.read_obs_var(f["var"])
    symbols = np.asarray(var.index, dtype=object)
    gene_ids = symbols
    wide = S.parse_wide_obs(obs, log)
    ruler_idx = S.ruler_column_index(gene_ids, symbols, frozen, log)
    dst = np.flatnonzero(ruler_idx >= 0)
    src = ruler_idx[dst]
    exp = build_expander(src, len(symbols))
    log(f"[align] ruler={n_ruler} present={dst.size} dup_file_cols="
        f"{int(src.size - np.unique(src).size)}")

    split = pd.read_json(SENG / "split.json", typ="series")
    sidx = np.load(SENG / "split_idx.npz")
    qidx = np.load(SENG / "prolif_quartile_idx.npz")
    axes = np.load(SENG / "axes.npz")

    kind = wide["kind"].to_numpy()
    labels = wide["labels"].to_numpy()
    modality = wide["modality"]
    pd_labels = wide["pd_labels"].to_numpy()
    is_wt = wide["is_wt"]
    nt_mask = kind == "NT"
    nt_mod = {m: np.flatnonzero(nt_mask & (modality == m)) for m in ("CRA", "CRI")}

    groups = []

    def add(name, idx, kindname):
        groups.append(dict(name=name, idx=np.asarray(idx, int), kind=kindname))

    add("EARLY_A", sidx["early_a"], "anchor")
    add("LATE_A", sidx["late_a"], "anchor")
    for m in ("CRA", "CRI"):
        add(f"NT_{m}", nt_mod[m], "nt")
        add(f"TOP_{m}", qidx[f"{m}_top"], "quartile")
        add(f"BOT_{m}", qidx[f"{m}_bottom"], "quartile")
    qual = []
    for lab, pos in pd.Series(np.arange(len(labels))).groupby(labels, sort=True).groups.items():
        lab = str(lab)
        if not (S.CRA_RE.fullmatch(lab) or S.CRI_RE.fullmatch(lab)):
            continue
        idx = np.asarray(pos, int)
        m = "CRA" if lab.startswith("CRA_") else "CRI"
        if not np.all(modality[idx] == m):
            log(f"[panel] {lab}: modality disagrees with prefix, skipped")
            continue
        if idx.size < MIN_CELLS:
            continue
        qual.append((m, lab, idx))
    qual.sort(key=lambda t: (t[0], t[1]))
    for m, lab, idx in qual:
        add(lab, idx, "perturbation")
    log(f"[panel] qualifying perturbations n>={MIN_CELLS}: {len(qual)}")

    return dict(f=f, obs=obs, var=var, symbols=symbols, wide=wide, frozen=frozen,
                n_ruler=n_ruler, dst=dst, src=src, exp=exp, groups=groups,
                split=split, sidx=sidx, qidx=qidx, axes=axes, kind=kind,
                labels=labels, modality=modality, pd_labels=pd_labels,
                is_wt=is_wt, nt_mod=nt_mod, qual=qual,
                n_counts=obs["n_counts"].to_numpy(np.float64))


def panel_z_from_counts(count_rows, frozen):
    """SENG's panel path, read frozen: joint TMM, log2-CPM, frozen GTEx z."""
    C = np.vstack([np.asarray(r, np.float64).ravel() for r in count_rows])
    logcpm, nf = tmm_logcpm_quiet(C)
    Z, missing = zscore_frozen(logcpm, frozen, C)
    return Z, missing, nf, C


# --------------------------------------------------------------------------
# thinning (amendment A3)
# --------------------------------------------------------------------------
def thin_counts(vec, p, rng):
    """Binomial thinning of a pseudobulk by a common per-group factor.

    Independent binomials with a common p sum to a binomial of the sum, so
    thinning the group pseudobulk is exactly thinning every cell by p.
    """
    if p >= 1.0:
        return np.asarray(vec, np.float64).copy()
    v = np.asarray(vec, np.float64)
    out = np.zeros_like(v)
    nz = v > 0
    out[nz] = rng.binomial(v[nz].astype(np.int64), p).astype(np.float64)
    return out


def thin_factors(medians, names):
    """target = lowest group median; factor = target / that group's median."""
    target = float(min(medians[n] for n in names))
    return target, {n: min(1.0, target / float(medians[n])) for n in names}


def fmt(x, nd=6):
    if x is None or (isinstance(x, float) and not np.isfinite(x)):
        return "n/a"
    return f"{x:.{nd}f}"


def binom_ci(k, n):
    """Wilson interval, used for the chance baseline beside every count."""
    if n == 0:
        return (np.nan, np.nan)
    from scipy.stats import beta
    lo = 0.0 if k == 0 else float(beta.ppf(0.025, k, n - k + 1))
    hi = 1.0 if k == n else float(beta.ppf(0.975, k + 1, n - k))
    return lo, hi
