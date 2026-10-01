"""Low-dimensional identifiable target vector.

Fit the age direction in a space with fewer dimensions than donors, then
re-test the identifiability failures from FINDINGS_TARGET.md (bootstrap
angles, site transfer, age-range split, junk gene weights).

Dimensionality reduction is fit on TRAINING FOLDS ONLY. Seed 20260914.
Namespaced under results/lowdim/. Does not modify prior FINDINGS*.md or
FALSIFICATION.md. No perturbation data, no TF Atlas, no candidate interventions.
"""
from __future__ import annotations

import json
import sys
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.cluster.hierarchy import fcluster, linkage
from scipy.spatial.distance import pdist
from scipy import stats

sys.path.insert(0, str(Path(__file__).resolve().parent))
from config import RESULTS, DATA_PROC  # noqa: E402
from target_common import (  # noqa: E402
    TARGET_DIR, TARGET_SEED, DATASET_ID, MIN_N_TYPE, N_BOOT_RIDGE, N_PERM,
    N_PERM_PAIR, STOP_NULL, ALPHAS, HASH_POOL_LIMITATION,
    Logger, dump_json, load_json, load_phase1_matrix, make_folds,
    permute_age_within_site_obs, score_age_predictions, primary_r2,
    unsigned_angle_deg, unit, align_age_sign, ridge_prestd, age_direction,
    donor_bootstrap_index, pairwise_median_deg, pairwise_angles_deg,
    calibrate_1d, zscore_train, identity_basis, target_from_age,
    consensus_direction, md_table, per_type_primary_r2,
)
from trajectory_common import permutation_p, summarize_null_col, fmt, fmt_u, primary_mae  # noqa: E402
from brain_phase1_common import r2_mae  # noqa: E402
from geometry_partC import proj_r2, angle_to_subspace_deg  # noqa: E402

LD_DIR = RESULTS / "lowdim"
LD_FIG = LD_DIR / "figures"
LD_CACHE = LD_DIR / "cache"
for _p in (LD_DIR, LD_FIG, LD_CACHE):
    _p.mkdir(parents=True, exist_ok=True)

LD_SEED = TARGET_SEED  # 20260914
KS = (20, 50, 100, 150)
KMAX = max(KS)
VARIANTS = ("pca", "curated", "wgcna")
VARIANT_LABEL = {
    "pca": "A_pca",
    "curated": "B_curated",
    "wgcna": "C_wgcna",
}
TECH_R_FLAG = 0.5
N_HVG = 2000
MIN_SET_GENES = 10
MIN_KEEP = 3
N_BOOT = N_BOOT_RIDGE  # 20
# Materially below gene-space 45.1°: a drop of at least 10°.
ANGLE_MATERIAL = 35.0
GENE_BOOT_ANGLE = 45.1
GENE_WITHIN_V2 = 0.231
GENE_WITHIN_G1 = 0.317
GENE_TRANSFER = {"H_to_M": -0.213, "M_to_H": -0.499}
GENE_YOUNG_OLD = 86.1
SITE_NAME = {"H": "HBCC", "M": "MSSM"}
ENRICHR_BASE = "https://maayanlab.cloud/Enrichr/geneSetLibrary"
ENRICHR_LIBS = (
    ("Hallmark", "MSigDB_Hallmark_2020"),
    ("Reactome", "Reactome_2022"),
    ("GO_BP", "GO_Biological_Process_2023"),
)
TECH_NAMES = ("site", "log_umi", "log_genes", "frac_neuronal")


def ld_log_banner(log, stage: str):
    log("=" * 100)
    log(f"LOWDIM {stage}  seed={LD_SEED}  dataset_id={DATASET_ID}")
    log("=" * 100)
    log(HASH_POOL_LIMITATION)
    log("NO perturbation data. No TF Atlas. No candidate interventions.")
    log("Dimensionality reduction is fit on training folds only (no full-cohort PCA/NMF then held-out).")
    log("FALSIFICATION.md and prior FINDINGS*.md are not modified.")
    log(f"k in {list(KS)}; variants={list(VARIANTS)}; tech |r| flag={TECH_R_FLAG}")


def attach_tech(obs: pd.DataFrame) -> pd.DataFrame:
    """Add per-row technical covariates used to flag components."""
    out = obs.copy()
    out["log_umi"] = np.log(np.clip(out.mean_counts_per_cell.to_numpy(float), 1e-12, None))
    out["log_genes"] = np.log(np.clip(out.mean_genes_per_cell.to_numpy(float), 1e-12, None))
    src = out.Source.astype(str).to_numpy()
    out["site_code"] = (src == "M").astype(np.float64)
    if "p0_frac_neuronal" in out.columns:
        out["frac_neuronal"] = out["p0_frac_neuronal"].to_numpy(float)
    else:
        out["frac_neuronal"] = np.nan
    return out


def tech_matrix(obs: pd.DataFrame) -> np.ndarray:
    """n x 4: site dummy, log UMI/nucleus, log genes/nucleus, donor neuronal fraction."""
    return np.column_stack([
        obs.site_code.to_numpy(float),
        obs.log_umi.to_numpy(float),
        obs.log_genes.to_numpy(float),
        obs.frac_neuronal.to_numpy(float),
    ])


def cap_k(k, n, p=None):
    k = int(k)
    hi = max(1, int(n) - 2)
    if p is not None:
        hi = min(hi, int(p))
    return int(max(1, min(k, hi)))


def spearman_safe(a, b):
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


def flag_components(scores, tech, thresh=TECH_R_FLAG):
    """scores n x k, tech n x 4. Vectorized Spearman via ranks + z-score multiply."""
    scores = np.asarray(scores, float)
    tech = np.asarray(tech, float)
    n, k = scores.shape
    nt = tech.shape[1]
    if n < 8 or k == 0:
        return np.full((k, nt), np.nan), np.zeros(k, dtype=bool)

    def _rank_z(X):
        R = np.empty(X.shape, dtype=np.float64)
        ok = np.ones(X.shape[1], dtype=bool)
        for j in range(X.shape[1]):
            col = X[:, j]
            m = np.isfinite(col)
            if int(m.sum()) < 8 or float(np.std(col[m])) <= 1e-12:
                R[:, j] = np.nan
                ok[j] = False
                continue
            r = np.full(len(col), np.nan)
            r[m] = stats.rankdata(col[m])
            mu = float(np.nanmean(r))
            sd = float(np.nanstd(r))
            if not np.isfinite(sd) or sd < 1e-12:
                R[:, j] = np.nan
                ok[j] = False
            else:
                R[:, j] = (r - mu) / sd
        return R, ok

    A, okA = _rank_z(scores)
    B, okB = _rank_z(tech)
    row = np.isfinite(B).all(1)
    if int(row.sum()) < 8:
        return np.full((k, nt), np.nan), np.zeros(k, dtype=bool)
    A2 = np.nan_to_num(A[row], nan=0.0)
    B2 = np.nan_to_num(B[row], nan=0.0)
    Rmat = (A2.T @ B2) / float(row.sum())
    Rmat[~okA] = np.nan
    Rmat[:, ~okB] = np.nan
    flagged = np.any(np.abs(Rmat) > thresh, axis=1)
    flagged = np.where(np.isfinite(Rmat).any(1), flagged, False)
    return Rmat, flagged


def zscore_cols(X):
    X = np.asarray(X, np.float64)
    if X.ndim == 1:
        X = X[:, None]
    mu = X.mean(0)
    sd = X.std(0)
    sd = np.where(sd < 1e-12, 1.0, sd)
    return (X - mu) / sd, mu, sd


def raw_scores(Xz, spec):
    Xz = np.asarray(Xz, np.float64)
    if spec["kind"] == "dense":
        return Xz @ spec["L"]
    members = spec["members"]
    S = np.empty((len(Xz), len(members)), dtype=np.float64)
    for j, idx in enumerate(members):
        idx = np.asarray(idx, dtype=int)
        if idx.size == 0:
            S[:, j] = 0.0
        else:
            S[:, j] = Xz[:, idx].mean(1)
    return S


def gene_map(w, spec, sd_comp, p):
    """Gradient of (z-scored component scores @ w) w.r.t. gene-z X."""
    w = np.asarray(w, np.float64).ravel()
    sd = np.asarray(sd_comp, np.float64).ravel()
    sd = np.where(np.abs(sd) < 1e-12, 1.0, sd)
    scale = w / sd
    v = np.zeros(int(p), dtype=np.float64)
    if spec["kind"] == "dense":
        L = np.asarray(spec["L"], np.float64)
        if L.size == 0 or scale.size == 0:
            return v
        return L @ scale
    for j, idx in enumerate(spec["members"]):
        idx = np.asarray(idx, dtype=int)
        if idx.size == 0 or j >= scale.size:
            continue
        v[idx] += scale[j] / float(idx.size)
    return v


class SpaceFit:
    """Max-k space fit on a training matrix Xz (n x p, already gene-z on train)."""

    def __init__(self, variant, n, p, kmax):
        self.variant = variant
        self.n = int(n)
        self.p = int(p)
        self.kmax = int(kmax)
        self.L = None
        self.names = []
        self.members = None          # curated: ranked list of index arrays
        self.Z_link = None           # wgcna linkage
        self.hvg_idx = None
        self.evals = None
        self.extra = {}

    def k_eff(self, k):
        return cap_k(k, self.n, self.p)

    def spec(self, k):
        ke = self.k_eff(k)
        if self.variant == "pca":
            ke = min(ke, 0 if self.L is None else self.L.shape[1])
            if ke < 1:
                return dict(kind="dense", L=np.zeros((self.p, 0)), names=[])
            return dict(kind="dense", L=self.L[:, :ke], names=self.names[:ke])
        if self.variant == "curated":
            ke = min(ke, 0 if self.members is None else len(self.members))
            if ke < 1:
                return dict(kind="members", members=[], names=[])
            return dict(kind="members", members=self.members[:ke], names=self.names[:ke])
        # wgcna: cut tree at ke
        return self._wgcna_spec(ke)

    def _wgcna_spec(self, ke):
        if self.Z_link is None or self.hvg_idx is None or ke < 1:
            return dict(kind="members", members=[], names=[])
        n_hvg = len(self.hvg_idx)
        ke = min(ke, n_hvg)
        labels = fcluster(self.Z_link, ke, criterion="maxclust")
        members, names = [], []
        for m in range(1, ke + 1):
            idx = self.hvg_idx[np.flatnonzero(labels == m)]
            members.append(np.asarray(idx, dtype=int))
            names.append(f"M{m}_n{len(idx)}")
        return dict(kind="members", members=members, names=names)


def fit_pca(Xz, kmax) -> SpaceFit:
    Xz = np.ascontiguousarray(Xz, np.float64)
    n, p = Xz.shape
    kmax = cap_k(kmax, n, p)
    sp = SpaceFit("pca", n, p, kmax)
    if n < 5 or p < 2 or kmax < 1:
        sp.L = np.zeros((p, 0))
        return sp
    # n << p: eigendecompose the n×n Gram instead of a thin SVD that materialises Vᵀ.
    G = Xz @ Xz.T
    G = 0.5 * (G + G.T)
    evals, evecs = np.linalg.eigh(G)
    evals = evals[::-1]
    evecs = evecs[:, ::-1]
    keep = evals > 1e-10
    if not np.any(keep):
        sp.L = np.zeros((p, 0))
        return sp
    evals, evecs = evals[keep], evecs[:, keep]
    kmax = min(kmax, int(evals.size))
    evals, evecs = evals[:kmax], evecs[:, :kmax]
    scale = 1.0 / np.sqrt(np.clip(evals, 1e-12, None))
    sp.kmax = kmax
    sp.L = (Xz.T @ evecs) * scale
    sp.evals = evals / max(n - 1, 1)
    sp.names = [f"PC{i+1}" for i in range(kmax)]
    tot = float(np.sum(np.clip(evals, 0.0, None)))  # not total variance if truncated
    sp.extra = dict(n_sv=int(kmax), gram=True)
    return sp


_SET_MAT = {}


def set_score_matrix(sets, p, min_genes=MIN_SET_GENES):
    """Sparse p x n_valid mean-operator for curated module scores (cached)."""
    from scipy.sparse import csr_matrix
    key = (id(sets), int(p), int(min_genes))
    if key in _SET_MAT:
        return _SET_MAT[key]
    valid = [s for s in sets if int(s["n_mapped"]) >= min_genes]
    rows, cols, data = [], [], []
    for j, s in enumerate(valid):
        idx = np.asarray(s["idx"], dtype=int)
        if idx.size == 0:
            continue
        w = 1.0 / float(idx.size)
        rows.extend(idx.tolist())
        cols.extend([j] * int(idx.size))
        data.extend([w] * int(idx.size))
    M = csr_matrix((data, (rows, cols)), shape=(int(p), len(valid)), dtype=np.float64)
    _SET_MAT[key] = (valid, M)
    return valid, M


def fit_curated(Xz, kmax, sets, ct, types, min_genes=MIN_SET_GENES) -> SpaceFit:
    """Shared gene-set basis. Rank by mean within-type variance of the module score."""
    Xz = np.asarray(Xz, np.float64)
    n, p = Xz.shape
    sp = SpaceFit("curated", n, p, kmax)
    valid, M = set_score_matrix(sets or [], p, min_genes=min_genes)
    if not valid or n < 5 or M.shape[1] == 0:
        sp.members, sp.names = [], []
        sp.kmax = 0
        return sp
    scores = np.ascontiguousarray(Xz @ M)
    acc = np.zeros(len(valid), dtype=np.float64)
    n_t = 0
    ct = np.asarray(ct).astype(str)
    for t in types:
        m = ct == t
        if int(m.sum()) < MIN_N_TYPE:
            continue
        acc += scores[m].var(0)
        n_t += 1
    rank_var = acc / max(n_t, 1)
    hall_i = [i for i, s in enumerate(valid) if s["library"] == "Hallmark"]
    other_i = [i for i, s in enumerate(valid) if s["library"] != "Hallmark"]
    hall_i = sorted(hall_i, key=lambda i: -rank_var[i])
    other_i = sorted(other_i, key=lambda i: -rank_var[i])
    kmax = cap_k(kmax, n, len(valid))
    if kmax <= len(hall_i):
        order = hall_i[:kmax]
    else:
        order = hall_i + other_i[: max(0, kmax - len(hall_i))]
        order = order[:kmax]
    sp.kmax = len(order)
    sp.members = [np.asarray(valid[i]["idx"], dtype=int) for i in order]
    sp.names = [str(valid[i]["name"]) for i in order]
    sp.extra = dict(
        n_valid_sets=int(len(valid)),
        n_hallmark=int(len(hall_i)),
        n_other=int(len(other_i)),
        n_types_for_rank=int(n_t),
        rank_var=[float(rank_var[i]) for i in order],
        libraries=[str(valid[i]["library"]) for i in order],
    )
    return sp


def fit_wgcna(Xz, kmax, var_raw, n_hvg=N_HVG) -> SpaceFit:
    """Ward clustering of HVGs in sample-space (WGCNA-style modules).

    Not full WGCNA (no soft-threshold / TOM); tree is cut at each k. Fit on train only.
    Average-linkage on correlation distance is the fallback if Ward fails; it is
    not the primary method (it chained to singletons on this matrix).
    """
    Xz = np.asarray(Xz, np.float64)
    n, p = Xz.shape
    sp = SpaceFit("wgcna", n, p, kmax)
    var_raw = np.asarray(var_raw, float) if var_raw is not None else Xz.var(0)
    var_raw = np.where(np.isfinite(var_raw), var_raw, -1.0)
    n_hvg = int(min(n_hvg, p, max(n_hvg, 1)))
    order = np.argsort(var_raw)[::-1]
    picked = []
    for i in order:
        if var_raw[i] <= 1e-12:
            break
        if float(np.std(Xz[:, i])) < 1e-12:
            continue
        picked.append(int(i))
        if len(picked) >= n_hvg:
            break
    hvg = np.asarray(picked, dtype=int)
    kmax = cap_k(kmax, n, len(hvg))
    sp.kmax = kmax
    sp.hvg_idx = hvg
    sp.extra = dict(n_hvg_requested=int(n_hvg), n_hvg=int(len(hvg)))
    if len(hvg) < 4 or n < 5 or kmax < 1:
        return sp
    Xh = np.ascontiguousarray(Xz[:, hvg])
    # Ward on genes in sample-space (already gene-z). Average-linkage on
    # correlation distance chained into singletons on this matrix.
    try:
        sp.Z_link = linkage(Xh.T, method="ward")
    except Exception:
        y = pdist(Xh.T, metric="correlation")
        y = np.clip(np.where(np.isfinite(y), y, 1.0), 0.0, 2.0)
        sp.Z_link = linkage(y, method="average")
    sizes = []
    spec = sp._wgcna_spec(kmax)
    sizes = [len(m) for m in spec["members"]]
    sp.extra["module_size_kmax_median"] = float(np.median(sizes)) if sizes else np.nan
    sp.extra["module_size_kmax_min"] = int(min(sizes)) if sizes else 0
    return sp


def fit_space(Xz, variant, kmax, *, var_raw=None, sets=None, ct=None, types=None):
    if variant == "pca":
        return fit_pca(Xz, kmax)
    if variant == "curated":
        return fit_curated(Xz, kmax, sets or [], ct, types or [])
    if variant == "wgcna":
        return fit_wgcna(Xz, kmax, var_raw)
    raise ValueError(variant)


def apply_slice(space, Xtr, Xte, k, tech_tr, drop_tech):
    """Train-only component z-score; optional drop of |r|>0.5 technical components.

    Returns dict with Ztr, Zte, w-ready matrices, spec_kept, sd, keep, flags, reverted.
    """
    k = cap_k(k, len(Xtr), space.p)
    spec = space.spec(k)
    names = spec["names"]
    k_req = len(names)
    empty = dict(
        Ztr=np.zeros((len(Xtr), 0)), Zte=np.zeros((len(Xte), 0)),
        spec=spec, spec_keep=spec, sd=np.ones(0), mu=np.zeros(0),
        keep=np.zeros(k_req, dtype=bool), flags=np.zeros(k_req, dtype=bool),
        Rtech=np.zeros((k_req, 4)), reverted=False, names=names, k_req=k_req,
        k_keep=0,
    )
    if k_req < 1 or len(Xtr) < MIN_N_TYPE:
        return empty
    raw_tr = raw_scores(Xtr, spec)
    raw_te = raw_scores(Xte, spec)
    Rtech, flags = flag_components(raw_tr, tech_tr)
    keep = np.ones(k_req, dtype=bool)
    reverted = False
    if drop_tech:
        keep = ~flags
        if int(keep.sum()) < MIN_KEEP:
            keep = np.ones(k_req, dtype=bool)
            reverted = True
    Ztr_k = raw_tr[:, keep]
    Zte_k = raw_te[:, keep]
    Ztr, mu, sd = zscore_cols(Ztr_k)
    Zte = (Zte_k - mu) / sd
    spec_keep = _take_keep(spec, keep)
    return dict(
        Ztr=Ztr, Zte=Zte, spec=spec, spec_keep=spec_keep, sd=sd, mu=mu,
        keep=keep, flags=flags, Rtech=Rtech, reverted=reverted,
        names=names, k_req=k_req, k_keep=int(keep.sum()),
    )


def _take_keep(spec, keep):
    keep = np.asarray(keep, dtype=bool)
    names = [n for n, k in zip(spec["names"], keep) if k]
    if spec["kind"] == "dense":
        return dict(kind="dense", L=spec["L"][:, keep], names=names)
    members = [m for m, k in zip(spec["members"], keep) if k]
    return dict(kind="members", members=members, names=names)


class TypeProjector:
    """Cache component scores and tech flags so each (type, boot/fold) pays SVD once."""

    def __init__(self, space, Xtr, tech_tr, Xte=None):
        self.space = space
        self.Xtr = np.asarray(Xtr, np.float64)
        self.Xte = None if Xte is None else np.asarray(Xte, np.float64)
        self.tech_tr = np.asarray(tech_tr, float)
        self.nested = space.variant in ("pca", "curated")
        self._max = None
        self._by_k = {}

    def _empty(self, spec):
        names = spec["names"]
        k_req = len(names)
        nte = 0 if self.Xte is None else len(self.Xte)
        return dict(
            Ztr=np.zeros((len(self.Xtr), 0)), Zte=np.zeros((nte, 0)),
            spec=spec, spec_keep=spec, sd=np.ones(0), mu=np.zeros(0),
            keep=np.zeros(k_req, dtype=bool), flags=np.zeros(k_req, dtype=bool),
            Rtech=np.zeros((k_req, 4)), reverted=False, names=names, k_req=k_req,
            k_keep=0,
        )

    def _compute(self, ke):
        spec = self.space.spec(ke)
        raw_tr = raw_scores(self.Xtr, spec)
        if self.Xte is None:
            raw_te = np.zeros((0, len(spec["names"])))
        else:
            raw_te = raw_scores(self.Xte, spec)
        R, flags = flag_components(raw_tr, self.tech_tr)
        return dict(spec=spec, raw_tr=raw_tr, raw_te=raw_te, R=R, flags=flags)

    def pack(self, k, drop_tech):
        ke = cap_k(k, len(self.Xtr), self.space.p)
        if ke < 1 or len(self.Xtr) < MIN_N_TYPE:
            return self._empty(self.space.spec(max(ke, 0)))
        if self.nested:
            if self._max is None:
                self._max = self._compute(self.space.kmax if self.space.kmax else ke)
            base = self._max
            ke = min(ke, len(base["spec"]["names"]))
            if ke < 1:
                return self._empty(base["spec"])
            mask = np.zeros(len(base["spec"]["names"]), dtype=bool)
            mask[:ke] = True
            spec = _take_keep(base["spec"], mask)
            raw_tr = base["raw_tr"][:, :ke]
            raw_te = base["raw_te"][:, :ke]
            R = base["R"][:ke]
            flags = base["flags"][:ke]
        else:
            if ke not in self._by_k:
                self._by_k[ke] = self._compute(ke)
            base = self._by_k[ke]
            spec, raw_tr, raw_te, R, flags = base["spec"], base["raw_tr"], base["raw_te"], base["R"], base["flags"]
        k_req = len(spec["names"])
        if k_req < 1:
            return self._empty(spec)
        keep = np.ones(k_req, dtype=bool)
        reverted = False
        if drop_tech:
            keep = ~flags
            if int(keep.sum()) < MIN_KEEP:
                keep = np.ones(k_req, dtype=bool)
                reverted = True
        Ztr, mu, sd = zscore_cols(raw_tr[:, keep])
        if raw_te.shape[0] == 0:
            Zte = np.zeros((0, int(keep.sum())))
        else:
            Zte = (raw_te[:, keep] - mu) / sd
        return dict(
            Ztr=Ztr, Zte=Zte, spec=spec, spec_keep=_take_keep(spec, keep),
            sd=sd, mu=mu, keep=keep, flags=flags, Rtech=R, reverted=reverted,
            names=spec["names"], k_req=k_req, k_keep=int(keep.sum()),
        )


def direction_and_pred(Ztr, ytr, Zte):
    w, rec = age_direction(Ztr, ytr, method="ridge")
    if not rec.get("ok") or not np.isfinite(w).all() or float(np.linalg.norm(w)) < 1e-12:
        mu = float(np.mean(ytr)) if len(ytr) else np.nan
        return w, rec, np.full(len(Zte), mu)
    pred, _, _ = calibrate_1d(Ztr @ w, ytr, Zte @ w)
    return w, rec, pred


def gene_direction_from_pack(pack, w, p):
    if pack["k_keep"] < 1 or not np.isfinite(w).all():
        return np.full(p, np.nan)
    return gene_map(w, pack["spec_keep"], pack["sd"], p)


def fit_pack(Y, obs, tech, idx, types, variant, sets, kmax=KMAX):
    """Train-only global z-score + per-type (or shared curated) space on row subset idx."""
    idx = np.asarray(idx)
    Xz, mu, sd = zscore_train(Y[idx])
    ct = obs.celltype.astype(str).to_numpy()[idx]
    y = obs.age.to_numpy(float)[idx]
    spaces = fit_spaces_by_type(Xz, Y[idx], ct, types, variant, kmax, sets)
    return dict(Xz=Xz, mu=mu, sd=sd, ct=ct, y=y, tech=tech[idx], spaces=spaces, idx=idx)


def fit_spaces_by_type(Xz, Yraw, ct, types, variant, kmax, sets, var_raw_all=None):
    """If variant is curated, one shared basis ranked on all rows; else per-type."""
    ct = np.asarray(ct).astype(str)
    spaces = {}
    if variant == "curated":
        shared = fit_space(Xz, "curated", kmax, sets=sets, ct=ct, types=types)
        for t in types:
            spaces[t] = shared
        return spaces
    for t in types:
        m = ct == t
        n = int(m.sum())
        if n < MIN_N_TYPE:
            spaces[t] = None
            continue
        var_raw = np.asarray(Yraw[m], np.float64).var(0) if variant == "wgcna" else None
        spaces[t] = fit_space(np.ascontiguousarray(Xz[m]), variant, kmax, var_raw=var_raw)
    return spaces


def load_bundle(log=print):
    data = load_phase1_matrix(log)
    Y, genes, obs = data["Y"], data["genes"], data["obs"]
    obs = attach_tech(obs)
    tech = tech_matrix(obs)
    n_nan = int(np.isnan(tech).sum())
    log(f"[filter] tech matrix {tech.shape}  NaNs={n_nan}  "
        f"site H/M={int((obs.Source.astype(str)=='H').sum())}/{int((obs.Source.astype(str)=='M').sum())} rows")
    log(f"[filter] log_umi finite={int(np.isfinite(obs.log_umi).sum())}  "
        f"frac_neuronal finite={int(np.isfinite(obs.frac_neuronal).sum())}")
    n0 = len(obs)
    assert n0 == Y.shape[0]
    types = sorted(pd.unique(obs.celltype.astype(str)))
    n_per = obs.groupby("celltype").donor.nunique().to_dict()
    log(f"[filter] types={len(types)}  donors/type min/median/max="
        f"{min(n_per.values())}/{float(np.median(list(n_per.values()))):.0f}/{max(n_per.values())}")
    for t, n in sorted(n_per.items(), key=lambda kv: kv[1]):
        k_caps = {k: cap_k(k, int((obs.celltype.astype(str) == t).sum())) for k in KS}
        log(f"[filter]   {t}: n_rows={(obs.celltype.astype(str)==t).sum()} n_donors={n} k_cap={k_caps}")
    data.update(dict(Y=Y, genes=genes, obs=obs, tech=tech, types=types, n_per_type=n_per))
    return data


def symbol_index(genes):
    sy = genes.symbol.astype(str).str.upper().to_numpy()
    mp = {}
    dups = 0
    for i, s in enumerate(sy):
        if s in mp:
            dups += 1
            continue
        mp[s] = i
    return mp, int(dups)


def _looks_like_gene(tok: str) -> bool:
    tok = tok.strip()
    if not tok or " " in tok or len(tok) > 25:
        return False
    return True


def parse_enrichr_library(text: str, library: str, sym_map: dict):
    rows = []
    for line in text.splitlines():
        line = line.strip()
        if not line:
            continue
        parts = line.split("\t")
        if len(parts) < 2:
            continue
        name = parts[0].strip()
        genes = []
        for p in parts[1:]:
            if "," in p and " " not in p.split(",")[0]:
                genes.extend(g.strip() for g in p.split(",") if g.strip())
            else:
                genes.append(p.strip())
        genes = [g for g in genes if _looks_like_gene(g)]
        idx = []
        seen = set()
        for g in genes:
            j = sym_map.get(g.upper())
            if j is None or j in seen:
                continue
            seen.add(j)
            idx.append(j)
        rows.append(dict(
            name=name, library=library, n_listed=len(genes),
            n_mapped=len(idx), idx=np.asarray(idx, dtype=int),
        ))
    return rows


def load_gene_sets(genes, log=print, force=False):
    """Download Hallmark / Reactome / GO_BP from Enrichr; cache under results/lowdim."""
    import requests

    cache_json = LD_DIR / "genesets_index.json"
    cache_dir = LD_DIR / "genesets"
    cache_dir.mkdir(parents=True, exist_ok=True)
    sym_map, n_dups = symbol_index(genes)
    log(f"[filter] gene symbols: {len(genes)}  unique-upper={len(sym_map)}  dup_skipped={n_dups}")
    if cache_json.exists() and not force:
        rec = load_json(cache_json)
        sets = []
        for row in rec["sets"]:
            sets.append(dict(
                name=row["name"], library=row["library"],
                n_listed=row["n_listed"], n_mapped=row["n_mapped"],
                idx=np.asarray(row["idx"], dtype=int),
            ))
        log(f"[filter] cached gene sets: {len(sets)} from {cache_json.name}")
        for lib, n in rec.get("n_per_library", {}).items():
            log(f"[filter]   {lib}: {n}")
        log(f"[filter] sets with >= {MIN_SET_GENES} mapped genes: "
            f"{sum(1 for s in sets if s['n_mapped'] >= MIN_SET_GENES)}")
        return sets, rec.get("filter", {})

    headers = {"User-Agent": "AgeIdentityLowdim/1.0 (research; seed 20260914)"}
    sets = []
    filt = dict(n_symbol_dups=n_dups, libraries={})
    for lib_name, lib_id in ENRICHR_LIBS:
        url = f"{ENRICHR_BASE}?mode=text&libraryName={lib_id}"
        path = cache_dir / f"{lib_id}.txt"
        text = None
        if path.exists() and not force:
            text = path.read_text(encoding="utf-8", errors="replace")
            log(f"[filter] gene-set library {lib_name}: read cache {path.name} ({len(text):,} chars)")
        else:
            log(f"[filter] gene-set library {lib_name}: GET {lib_id}")
            try:
                r = requests.get(url, headers=headers, timeout=120)
                r.raise_for_status()
                text = r.text
                path.write_text(text, encoding="utf-8")
                log(f"[filter]   wrote {path.name} ({len(text):,} chars)")
            except Exception as e:
                log(f"[filter]   DOWNLOAD FAILED {lib_name}: {e}")
                filt["libraries"][lib_name] = dict(ok=False, error=str(e))
                continue
        parsed = parse_enrichr_library(text, lib_name, sym_map)
        n_ok = sum(1 for s in parsed if s["n_mapped"] >= MIN_SET_GENES)
        log(f"[filter]   {lib_name}: listed={len(parsed)}  mapped>={MIN_SET_GENES}: {n_ok}")
        filt["libraries"][lib_name] = dict(
            ok=True, n_listed=len(parsed), n_mapped_ok=n_ok,
            median_mapped=float(np.median([s["n_mapped"] for s in parsed])) if parsed else 0,
        )
        sets.extend(parsed)
    serial = []
    for s in sets:
        serial.append(dict(
            name=s["name"], library=s["library"],
            n_listed=int(s["n_listed"]), n_mapped=int(s["n_mapped"]),
            idx=s["idx"].tolist(),
        ))
    n_per = {lib: int(sum(1 for s in sets if s["library"] == lib)) for lib, _ in ENRICHR_LIBS}
    payload = dict(
        seed=LD_SEED, min_set_genes=MIN_SET_GENES, n_sets=len(sets),
        n_per_library=n_per, filter=filt, sets=serial,
    )
    dump_json(cache_json, payload)
    log(f"[filter] wrote {cache_json}  total_sets={len(sets)}  "
        f"usable={sum(1 for s in sets if s['n_mapped'] >= MIN_SET_GENES)}")
    return sets, filt


def config_id(variant, k, clean):
    lab = VARIANT_LABEL[variant]
    tag = f"{lab}_k{int(k):03d}"
    if clean:
        tag += "_clean"
    return tag


def gene_space_baseline():
    """Numbers from FINDINGS_TARGET / G1, filled from saved json when present."""
    out = dict(
        boot_median_angle=GENE_BOOT_ANGLE,
        within_site_r2_v2=GENE_WITHIN_V2,
        within_site_r2_g1=GENE_WITHIN_G1,
        transfer_H_to_M=GENE_TRANSFER["H_to_M"],
        transfer_M_to_H=GENE_TRANSFER["M_to_H"],
        young_old_angle=GENE_YOUNG_OLD,
        source="FINDINGS_TARGET.md constants",
    )
    try:
        v1 = load_json(TARGET_DIR / "v1_summary.json")
        out["boot_median_angle"] = float(v1["ridge_median_pairwise_deg"])
        out["source"] = "results/target json"
    except Exception:
        pass
    try:
        v2 = load_json(TARGET_DIR / "v2_summary.json")
        rt = (v2.get("within_site") or {}).get("ridge_target") or {}
        if rt.get("primary_r2") is not None:
            out["within_site_r2_v2"] = float(rt["primary_r2"])
    except Exception:
        pass
    try:
        v4 = load_json(TARGET_DIR / "v4_summary.json")
        for row in v4.get("v4a", []):
            tr, te = str(row.get("train_site")), str(row.get("test_site"))
            if tr == "H" and te == "M":
                out["transfer_H_to_M"] = float(row["median_r2"])
            if tr == "M" and te == "H":
                out["transfer_M_to_H"] = float(row["median_r2"])
        v4b = v4.get("v4b") or {}
        if v4b.get("median_angle_deg") is not None:
            out["young_old_angle"] = float(v4b["median_angle_deg"])
    except Exception:
        pass
    try:
        g1 = load_json(RESULTS / "trajectory" / "g1_summary.json")
        ridge = ((g1.get("schemes") or {}).get("within_site") or {}).get("ridge") or {}
        if ridge.get("site_strat_r2") is not None:
            out["within_site_r2_g1"] = float(ridge["site_strat_r2"])
        elif ridge.get("r2") is not None:
            out["within_site_r2_g1"] = float(ridge["r2"])
    except Exception:
        pass
    return out


def junk_symbol(sym: str) -> bool:
    import re
    s = str(sym)
    if s.startswith("ENSG") or s.startswith("LINC") or s.startswith("OR"):
        return True
    if "-AS" in s or s.endswith("-DT") or re.search(r"P\d+$", s):
        return True
    return False


class FilterLog:
    def __init__(self):
        self.rows = []

    def add(self, stage, step, n, **kw):
        rec = dict(stage=stage, step=step, n=n)
        rec.update(kw)
        self.rows.append(rec)

    def to_csv(self, path):
        pd.DataFrame(self.rows).to_csv(path, index=False)


def median_pairwise_unit_cols(cols):
    cols = [np.asarray(c, float) for c in cols
            if np.isfinite(c).all() and float(np.linalg.norm(c)) > 1e-12]
    if len(cols) < 2:
        return np.nan, 0
    Wc = np.column_stack([c / np.linalg.norm(c) for c in cols])
    ang = pairwise_angles_deg(Wc)
    iu = np.triu_indices(ang.shape[0], 1)
    return float(np.median(ang[iu])), int(len(iu[0]))
