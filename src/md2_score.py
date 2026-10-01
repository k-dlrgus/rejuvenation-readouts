"""Seurat AddModuleScore (exact control-bin algorithm) and frozen-ruler scoring.

AddModuleScore is implemented from Seurat v5 AddModuleScore:
  gene means on the log-normalized matrix; nbin quantile bins with infinitesimal
  Gaussian noise to break ties (ggplot2::cut_number); ctrl genes sampled without
  replacement from the same bin as each feature gene; unique control pool;
  score = mean(features) − mean(controls) per cell.
Seed for sampling is house seed 20260914 (Seurat default seed=1 is not used).
"""
from __future__ import annotations

import sys
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import sparse

sys.path.insert(0, str(Path(__file__).resolve().parent))
from md2_common import (  # noqa: E402
    AMS_CTRL, AMS_NBIN, MD2_SEED, SEURAT_LOGNORMALIZE_SCALE, StopStep,
)
from fibro_common import unit  # noqa: E402
from fibro2_common import tmm_logcpm_rows, score_frozen  # noqa: E402
from fibro_stage2 import differentiation_score  # noqa: E402


def lognormalize_csr(X, scale=SEURAT_LOGNORMALIZE_SCALE):
    """Seurat NormalizeData LogNormalize: ln(1 + scale * count / tot). Sparse-preserving."""
    X = X.tocsr().copy()
    if not np.issubdtype(X.dtype, np.floating):
        X = X.astype(np.float64)
    lib = np.asarray(X.sum(axis=1)).ravel()
    lib = np.where(lib <= 0, 1.0, lib)
    factor = float(scale) / lib
    X.data = X.data * np.repeat(factor, np.diff(X.indptr))
    X.data = np.log1p(X.data)
    return X


def lognormalize_dense(C, scale=SEURAT_LOGNORMALIZE_SCALE):
    """Same LogNormalize on a dense samples × genes count matrix."""
    C = np.asarray(C, np.float64)
    lib = C.sum(axis=1, keepdims=True)
    lib = np.where(lib <= 0, 1.0, lib)
    return np.log1p(C / lib * float(scale))


def _cut_number(x, nbin, rng):
    """ggplot2::cut_number analogue: quantile bins, right=FALSE, nbin labels 0..nbin-1."""
    x = np.asarray(x, float)
    noise = rng.normal(loc=0.0, scale=1.0, size=x.shape[0]) / 1e30
    z = x + noise
    qs = np.quantile(z, np.linspace(0.0, 1.0, int(nbin) + 1))
    # Identical breaks (all-zero genes) would collapse bins. Noise of 1e-30
    # breaks exact ties; if still collapsed, fail loudly.
    for i in range(1, len(qs)):
        if qs[i] <= qs[i - 1]:
            qs[i] = np.nextafter(qs[i - 1], np.inf)
    # cut(..., right=FALSE, include.lowest=TRUE) → digitize with right=False,
    # plus a left-closed first bin.
    bins = np.searchsorted(qs[1:-1], z, side="right")
    if int(bins.min()) < 0 or int(bins.max()) >= int(nbin):
        raise StopStep(
            "addmodulescore",
            f"cut_number produced bin ids outside 0..{nbin-1}: "
            f"min={int(bins.min())} max={int(bins.max())}",
        )
    n_unique = int(np.unique(bins).size)
    if n_unique != int(nbin):
        raise StopStep(
            "addmodulescore",
            f"cut_number did not yield {nbin} bins (got {n_unique}). "
            "Not dropping bins.",
        )
    return bins.astype(int)


def _column_means_csr(X):
    return np.asarray(X.mean(axis=0)).ravel()


def _mean_cols_csr(X, idx):
    if idx is None or len(idx) == 0:
        return np.zeros(X.shape[0], dtype=np.float64)
    idx = np.asarray(idx, dtype=int)
    return np.asarray(X[:, idx].mean(axis=1)).ravel()


def symbol_index(symbols):
    """Uppercase symbol → first column index. Duplicates keep the first."""
    pos = {}
    for i, s in enumerate(symbols):
        u = str(s).upper()
        if u not in pos:
            pos[u] = i
    return pos


def map_features(symbols, genes, log, tag):
    pos = symbol_index(symbols)
    idx, missing, mapped = [], [], []
    seen = set()
    for g in genes:
        u = str(g).upper()
        j = pos.get(u)
        if j is None:
            missing.append(u)
            continue
        if j in seen:
            continue
        seen.add(j)
        idx.append(int(j))
        mapped.append(u)
    if log is not None:
        log(f"[AMS {tag}] requested={len(list(genes))} mapped={len(idx)} missing={len(missing)} "
            f"missing_symbols={missing[:20]}")
    return np.asarray(idx, dtype=int), mapped, missing


def add_module_score(
    data, symbols, gene_sets, log=None, nbin=AMS_NBIN, ctrl=AMS_CTRL,
    seed=MD2_SEED, tag="", bins=None, gene_mean=None,
):
    """AddModuleScore on a cells × genes matrix (CSR or dense).

    gene_sets: dict name -> list of symbols.
    Returns scores dict name->(n_cells,), plus meta (bins, means, mapped).
    """
    dense = not sparse.issparse(data)
    if dense:
        X = np.asarray(data, np.float64)
        n_cells, n_genes = X.shape
        if gene_mean is None:
            gene_mean = X.mean(axis=0)
    else:
        X = data.tocsr()
        n_cells, n_genes = X.shape
        if gene_mean is None:
            gene_mean = _column_means_csr(X)
    rng = np.random.default_rng(int(seed))
    if bins is None:
        bins = _cut_number(gene_mean, nbin, rng)
    # Fresh RNG stream for control sampling (after bin noise draw).
    rng_ctrl = np.random.default_rng(int(seed) + 1)
    out, meta = {}, {}
    for name, genes in gene_sets.items():
        feat_idx, mapped, missing = map_features(symbols, genes, log, f"{tag}:{name}")
        if feat_idx.size == 0:
            raise StopStep(
                "addmodulescore",
                f"{tag}:{name}: 0 feature genes mapped of {len(list(genes))}. Not scoring.",
            )
        ctrl_pool = []
        for j in feat_idx:
            b = int(bins[j])
            pool = np.flatnonzero(bins == b)
            if pool.size < int(ctrl):
                raise StopStep(
                    "addmodulescore",
                    f"{tag}:{name}: bin {b} has {pool.size} genes < ctrl={ctrl}. "
                    "Seurat samples without replacement; not reducing ctrl.",
                )
            ctrl_pool.extend(rng_ctrl.choice(pool, size=int(ctrl), replace=False).tolist())
        ctrl_idx = np.unique(np.asarray(ctrl_pool, dtype=int))
        if dense:
            feat_mean = X[:, feat_idx].mean(axis=1)
            ctrl_mean = X[:, ctrl_idx].mean(axis=1)
        else:
            feat_mean = _mean_cols_csr(X, feat_idx)
            ctrl_mean = _mean_cols_csr(X, ctrl_idx)
        out[name] = np.asarray(feat_mean - ctrl_mean, np.float64)
        meta[name] = dict(
            n_requested=int(len(list(genes))), n_mapped=int(feat_idx.size),
            n_missing=int(len(missing)), missing=list(missing), mapped=list(mapped),
            n_ctrl=int(ctrl_idx.size), nbin=int(nbin), ctrl_per_gene=int(ctrl),
            seed=int(seed),
        )
        if log is not None:
            log(f"[AMS {tag}:{name}] n_mapped={feat_idx.size} n_ctrl={ctrl_idx.size} "
                f"score_mean={float(np.mean(out[name])):+.4f}")
    return out, dict(
        nbin=int(nbin), ctrl=int(ctrl), seed=int(seed), n_cells=int(n_cells),
        n_genes=int(n_genes), gene_mean=np.asarray(gene_mean, np.float64),
        bins=np.asarray(bins, np.int32), per_set=meta,
        algorithm="Seurat_AddModuleScore_nbin24_ctrl100_LogNormalize",
    )


def draw_size_matched_indices(feature_idx, bins, rng, n_draws, exclude_extra=None):
    """200 draws of the same n genes, each gene replaced by another from its mean-expression bin.

    Without replacement within a draw. Excludes the original feature indices (and optional extras).
    """
    feature_idx = np.asarray(feature_idx, dtype=int)
    bins = np.asarray(bins, dtype=int)
    n = int(feature_idx.size)
    by_bin = {}
    for b in np.unique(bins):
        by_bin[int(b)] = np.flatnonzero(bins == int(b))
    blocked = set(int(i) for i in feature_idx)
    if exclude_extra is not None:
        blocked.update(int(i) for i in exclude_extra)
    draws = np.full((n_draws, n), -1, dtype=int)
    for d in range(int(n_draws)):
        used = set(blocked)
        row = np.empty(n, dtype=int)
        for k, gi in enumerate(feature_idx):
            pool = by_bin[int(bins[gi])]
            avail = pool[np.array([int(p) not in used for p in pool], dtype=bool)]
            if avail.size == 0:
                avail = pool[pool != int(gi)]
            if avail.size == 0:
                raise StopStep(
                    "md_null",
                    f"bin {int(bins[gi])} has no genes to draw for size-matched null",
                )
            pick = int(rng.choice(avail))
            row[k] = pick
            used.add(pick)
        draws[d] = row
    return draws


def ams_from_indices(data, feat_idx, bins, rng, ctrl=AMS_CTRL):
    """AddModuleScore given integer feature indices (for random gene-set draws)."""
    feat_idx = np.asarray(feat_idx, dtype=int)
    dense = not sparse.issparse(data)
    ctrl_pool = []
    for j in feat_idx:
        b = int(bins[j])
        pool = np.flatnonzero(bins == b)
        if pool.size < int(ctrl):
            raise StopStep("addmodulescore", f"bin {b} has {pool.size} < ctrl={ctrl}")
        ctrl_pool.extend(rng.choice(pool, size=int(ctrl), replace=False).tolist())
    ctrl_idx = np.unique(np.asarray(ctrl_pool, dtype=int))
    if dense:
        X = np.asarray(data, np.float64)
        feat_mean = X[:, feat_idx].mean(axis=1)
        ctrl_mean = X[:, ctrl_idx].mean(axis=1)
    else:
        feat_mean = _mean_cols_csr(data, feat_idx)
        ctrl_mean = _mean_cols_csr(data, ctrl_idx)
    return np.asarray(feat_mean - ctrl_mean, np.float64)


def score_frozen_pseudobulk(counts, frozen, log, tag):
    """counts: samples × ruler genes. TMM among these samples, frozen z, missing→0."""
    logcpm, nf = tmm_logcpm_rows(counts, log, tag=tag)
    age, Z, missing = score_frozen(logcpm, frozen, counts=counts)
    pluri_w, meta_w = differentiation_score(Z, frozen, drop_oskm=False, log=log)
    pluri_wo, meta_wo = differentiation_score(Z, frozen, drop_oskm=True, log=log)
    return dict(
        logcpm=logcpm, nf=nf, age=age, Z=Z, missing=missing,
        pluri_with=pluri_w, pluri_without=pluri_wo,
        pluri_meta_with=meta_w, pluri_meta_without=meta_wo,
    )
