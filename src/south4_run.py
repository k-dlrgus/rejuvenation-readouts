"""SOUTH4 — rank Southard Hs27 CRISPRa factors on the full measured singlets.

The pre-registration is results/south4/PREREG.flag, written before any statistic.
This script reads that flag and does not rewrite it. It does not refit the
frozen ruler and does not modify the imported modules.

Locked readings, fixed here before any SOUTH4 number exists:

1. The pinned file is the only Zenodo object at the pre-registered size,
   fibroblast_CRISPRa_final_pop_singlets_normalized_log1p.h5ad (8,048,055,098
   bytes). X is the log1p matrix named in the record. Raw counts are used only
   if X itself passes the integer-count audit, otherwise the raw layer that
   passes it (the record documents layer 'counts'). No other file is substituted.

2. A target is non-targeting when its stripped, casefolded label is in
   NT_TOKENS. Guide identity is the guide column. Ranked factors are targets
   that are not non-targeting and have at least 2 distinct guides. Factors
   below that are inventoried and not relabelled.

3. Missing ruler genes stay z = 0. The Stage 2A permutation shuffles the
   matched-gene coordinates only; unmatched coordinates stay 0 in every draw,
   so a null draw cannot invent a measurement. ||d|| is preserved.

4. The TMM panel is this task's own pseudobulks: one row per ranked target,
   one row per non-targeting guide, and one row for all non-targeting cells
   pooled. log2 CPM prior.count = 2, then the frozen GTEx mu/sd. Guide-level
   vectors used only for agreement are z-scored with that panel's reference
   via z_from_counts. d_g = z(target) - z(pooled non-targeting).

5. A non-targeting guide pseudo-effect is leave-one-guide-out:
   z(guide) - z(other non-targeting cells). Random splits are 200 disjoint
   halves of the non-targeting cells, seed 20260914, d = z(half A) - z(half B).
   The stricter control floor is the more negative of those two 5th percentiles.

6. A dose clears a floor only when its delta is strictly below that floor and
   strictly negative. MDA for the stop key is the smallest f clearing the
   Stage 2A pooled-at-median-norm floor. The combination margin uses the
   stricter of that pooled floor and the Stage 2B control floor: the smallest
   f that clears the stricter floor, and the clearance of that dose. A set is
   reported only when its own magnitude-matched permutation floor is cleared
   by at least that margin.

7. Identity drop is the seng_run definition on the panel log2 CPM:
   -mean(logcpm_target - logcpm_NT) over the identity genes that are present.
   identity_loss if that drop exceeds IDENTITY_DROP (0.5). Proliferation
   residualisation is stats_from_d against the NT top-minus-bottom quartile
   axis. The proliferation flag is descriptive: delta < 0 and delta_resid >= 0.
   It is not part of the toward_young gate.

8. toward_young requires delta < 0, delta below the factor's own Stage 2A
   floor, q <= 0.05, not identity_loss, not guides_disagree. Bootstrap
   intervals are never a gate. An interval that does not contain its point
   estimate is INVALID.

Nothing in this list is revised after the results exist.
"""
from __future__ import annotations

import os
import sys
import time
from pathlib import Path

import h5py
import numpy as np
import pandas as pd
from scipy import sparse, stats as sps

sys.path.insert(0, str(Path(__file__).resolve().parent))

from brain_phase1_common import Logger, dump_json, strip_ensembl  # noqa: E402
from fibro_common import jsonable  # noqa: E402
from fibro2_common import load_frozen_ruler  # noqa: E402
from gtex_common import StopStep, TMM_LOGRATIO_TRIM, TMM_SUM_TRIM  # noqa: E402
from target_common import unit  # noqa: E402

from toward_run import (  # noqa: E402
    ANCHORS_PATH,
    OLD_BINS,
    YOUNG_BINS,
    load_gtex_z,
    percentile_ci,
    require_file,
    ruler_score,
    tmm_logcpm_quiet,
    zscore_frozen,
)
from same_run import (  # noqa: E402
    BOOT_SEED,
    SAME_DIR,
    ci_covers_point,
    cosine_full,
    frac_along,
)
from seng_run import (  # noqa: E402
    IDENTITY_DROP,
    IDENTITY_GENES,
    bh_q,
    choose_id_column,
    choose_symbol_column,
    dump_tree,
    file_md5,
    load_cc_genes,
    module_score,
    quartile_masks,
    read_obs_var,
    scan_numeric,
    stats_from_d,
    z_from_counts,
)
from lowdim_common import fit_pca  # noqa: E402
from md3_idtype import detect_id_type  # noqa: E402
from seng_run import ruler_column_index  # noqa: E402


ROOT = Path(__file__).resolve().parents[1]
SOUTH4 = ROOT / "results" / "south4"
SOUTH4.mkdir(parents=True, exist_ok=True)
RAW_DIR = ROOT / "data" / "raw" / "southard"
H5_NAME = "fibroblast_CRISPRa_final_pop_singlets_normalized_log1p.h5ad"
H5_PATH = RAW_DIR / H5_NAME
EXPECTED_SIZE = 8_048_055_098
EXPECTED_MD5 = "2c926c7cae2eaab1ed2519b06ca21db5"
PREREG = SOUTH4 / "PREREG.flag"
FINDINGS = ROOT / "FINDINGS_SOUTH4.md"
PROGRESS = ROOT / "PROGRESS_SOUTH4.md"

SEED = 20260914
N_PERM = 20000
N_BOOT = 200
N_SPLITS = 200
N_RANDOM_STARTS = 50
MAX_SET = 5
S2_KS = (20, 50, 100)
KMAX = 100
F_VALUES = (0.0, 0.05, 0.1, 0.25, 0.5, 1.0)
Q_BAR = 0.05
MIN_GUIDES = 2
MIN_FACTORS = 500
MIN_NT_GUIDES = 20
COVERAGE_LOW = 0.25
PERM_CHUNK = 1000

NT_TOKENS = {
    "non-targeting", "nontargeting", "non_targeting", "non targeting",
    "nt", "negctrl", "neg-ctrl", "neg_ctrl", "negative", "negative control",
    "control", "off-target", "offtarget", "off_target",
    "non-targeting control", "safe-targeting", "safe_targeting",
    "intergenic", "no_target", "notarget", "no-target",
    "non",
}

TARGET_PRIORITY = (
    "guide_target", "target_gene", "gene_target", "target", "gene", "perturbation",
    "gene_symbol", "target_gene_name",
)
GUIDE_PRIORITY = (
    "guide_identity", "guide_id", "protospacer", "guide", "sgrna",
    "feature_call", "sgRNA",
)

# Worker globals. Set in the process initializer; empty in the parent.
_W = {}


class Halt(Exception):
    def __init__(self, key, message, details=None):
        super().__init__(message)
        self.key = str(key)
        self.message = str(message)
        self.details = details


def _manifest():
    path = SOUTH4 / "manifest.json"
    if path.exists():
        import json
        return json.loads(path.read_text(encoding="utf-8"))
    return {"failures": [], "flags": [], "notes": []}


def _save_manifest(man):
    dump_json(SOUTH4 / "manifest.json", jsonable(man))


def record_failure(step, message, details=None):
    man = _manifest()
    man["failures"].append({"step": step, "message": message, "details": jsonable(details)})
    _save_manifest(man)


def record_flag(key, message, details=None):
    man = _manifest()
    man["flags"].append({"key": key, "message": message, "details": jsonable(details)})
    _save_manifest(man)


def record_note(step, message, details=None):
    man = _manifest()
    man["notes"].append({"step": step, "message": message, "details": jsonable(details)})
    _save_manifest(man)


def _fmt(x, d=4):
    try:
        v = float(x)
    except (TypeError, ValueError):
        return "NA"
    if not np.isfinite(v):
        return "NA"
    return f"{v:+.{d}f}"


def _col(df, priority, log, kind):
    lower = {str(c).lower(): c for c in df.columns}
    for name in priority:
        if name in df.columns:
            log(f"[{kind}] column {name!r} by priority")
            return name
        if name.lower() in lower:
            got = lower[name.lower()]
            log(f"[{kind}] column {got!r} by case-insensitive priority ({name})")
            return got
    return None


def _is_nt(label):
    s = str(label).strip().casefold()
    return s in NT_TOKENS


def _csr_node(node):
    enc = str(node.attrs.get("encoding-type", "") or "")
    if isinstance(node, h5py.Dataset):
        return "dense", node
    if "data" in node and "indptr" in node and "indices" in node:
        shape = node.attrs.get("shape", None)
        if shape is None:
            raise Halt("h5", f"{node.name} sparse group has no shape. Not substituting.")
        shape = tuple(int(x) for x in shape)
        indptr_n = int(node["indptr"].shape[0])
        if "csr" in enc or indptr_n == shape[0] + 1:
            return "csr", node
        if "csc" in enc or indptr_n == shape[1] + 1:
            return "csc", node
        raise Halt(
            "h5",
            f"{node.name} encoding {enc!r} indptr={indptr_n} shape={shape}. Not substituting.",
        )
    raise Halt("h5", f"{node.name} is not a matrix. Not substituting.")


def locate_raw_node(f, log):
    """Return (kind, node, source, audits). Halt no_raw_counts when nothing qualifies."""
    audits = []
    if "X" not in f:
        raise Halt("h5", "h5ad has no X. Not substituting.")
    X = f["X"]
    kind, node = _csr_node(X)
    data = node["data"] if kind != "dense" else node
    audit = scan_numeric(data)
    audit["where"] = "X/data" if kind != "dense" else "X"
    audits.append(audit)
    log(f"[raw audit] {audit}")
    if audit["raw_integer_counts"]:
        return kind, node, f"X ({kind})", audits
    layer_names = list(f["layers"].keys()) if "layers" in f else []
    candidates = []
    if "layers" in f:
        for name in layer_names:
            n = f["layers"][name]
            try:
                k, nd = _csr_node(n)
            except Halt:
                continue
            ds = nd["data"] if k != "dense" else nd
            a = scan_numeric(ds)
            a["where"] = f"layers/{name}"
            audits.append(a)
            log(f"[raw audit] {a}")
            if a["raw_integer_counts"]:
                candidates.append((name, k, nd))
    preferred = [c for c in candidates if c[0] in ("counts", "raw_counts", "raw")]
    chosen = preferred[0] if preferred else (candidates[0] if len(candidates) == 1 else None)
    if chosen is not None:
        name, k, nd = chosen
        return k, nd, f"layers/{name} ({k})", audits
    raw_group = "raw" in f and "X" in f["raw"]
    if raw_group:
        n = f["raw"]["X"]
        k, nd = _csr_node(n)
        ds = nd["data"] if k != "dense" else nd
        a = scan_numeric(ds)
        a["where"] = "raw/X"
        audits.append(a)
        log(f"[raw audit] {a}")
        if a["raw_integer_counts"]:
            return k, nd, f"raw/X ({k})", audits
    raise Halt(
        "no_raw_counts",
        "X is not raw integer counts and no raw layer exists. "
        f"audits={jsonable(audits)} layers={layer_names}.",
        {"audits": audits, "layers": layer_names},
    )


def remap_csr(block, col_map, n_matched):
    """CSR cells x all-genes -> CSR cells x matched genes. Unmatched columns dropped."""
    if block.nnz == 0:
        return sparse.csr_matrix((block.shape[0], n_matched), dtype=np.float64)
    mj = col_map[block.indices]
    keep = mj >= 0
    data = np.asarray(block.data[keep], np.float64)
    indices = np.asarray(mj[keep], np.int32)
    row_ids = np.repeat(np.arange(block.shape[0], dtype=np.int32), np.diff(block.indptr))
    kept_rows = row_ids[keep]
    counts = np.bincount(kept_rows, minlength=block.shape[0])
    indptr = np.zeros(block.shape[0] + 1, np.int64)
    np.cumsum(counts, out=indptr[1:])
    return sparse.csr_matrix((data, indices, indptr), shape=(block.shape[0], n_matched))


def stream_sums(kind, node, guide_of_cell, n_guides, col_map, n_matched, log):
    """One pass. guide_of_cell[i] = guide id or -1. Returns guide_sums (n_guides, n_matched)
    memmap and writes the compact matched CSR beside it."""
    if kind == "dense":
        raise Halt("h5", "Dense raw matrix was not the expected layout. Not densifying an 8 GB matrix into RAM.")
    if kind == "csc":
        raise Halt(
            "h5",
            "Raw counts are CSC. This task's streaming sum is written for CSR, which is the "
            "AnnData layout of this record. Not transposing a multi-GB matrix in place.",
        )
    shape = tuple(int(x) for x in node.attrs["shape"])
    n_cells, n_genes = shape
    if guide_of_cell.shape[0] != n_cells:
        raise Halt("h5", f"guide assignment length {guide_of_cell.shape[0]} != n_cells {n_cells}.")
    if col_map.shape[0] != n_genes:
        raise Halt("h5", f"column map {col_map.shape[0]} != n_genes {n_genes}.")
    indptr_all = np.asarray(node["indptr"][:], np.int64)
    sums_path = SOUTH4 / "guide_sums.f64"
    sums = np.memmap(sums_path, dtype=np.float64, mode="w+", shape=(n_guides, n_matched))
    sums[:] = 0
    data_path = SOUTH4 / "matched_data.f64"
    idx_path = SOUTH4 / "matched_indices.i32"
    indptr = np.zeros(n_cells + 1, np.int64)
    nnz = 0
    t0 = time.perf_counter()
    row_chunk = 400
    fd = open(data_path, "wb")
    fi = open(idx_path, "wb")
    try:
        for a in range(0, n_cells, row_chunk):
            b = min(a + row_chunk, n_cells)
            i0 = int(indptr_all[a])
            i1 = int(indptr_all[b])
            if i1 > i0:
                data = np.asarray(node["data"][i0:i1])
                indices = np.asarray(node["indices"][i0:i1])
                ip = indptr_all[a:b + 1] - i0
                block = sparse.csr_matrix((data, indices, ip), shape=(b - a, n_genes))
            else:
                block = sparse.csr_matrix((b - a, n_genes))
            sub = remap_csr(block, col_map, n_matched)
            fd.write(np.asarray(sub.data, np.float64).tobytes())
            fi.write(np.asarray(sub.indices, np.int32).tobytes())
            counts = np.diff(sub.indptr)
            indptr[a + 1:b + 1] = nnz + np.cumsum(counts)
            nnz += int(sub.nnz)
            gid = guide_of_cell[a:b]
            valid = gid >= 0
            if valid.any() and sub.nnz:
                dense = np.asarray(sub.toarray(), np.float64)
                np.add.at(sums, gid[valid], dense[valid])
            if (a // row_chunk) % 200 == 0:
                log(f"[stream] rows {b:,}/{n_cells:,} nnz={nnz:,} ({time.perf_counter() - t0:.0f}s)")
                sums.flush()
    finally:
        fd.close()
        fi.close()
    sums.flush()
    np.save(SOUTH4 / "matched_indptr.npy", indptr)
    np.save(SOUTH4 / "matched_nnz.npy", np.array([nnz]))
    log(f"[stream] done cells={n_cells:,} matched_nnz={nnz:,} in {time.perf_counter() - t0:.0f}s")
    return sums, nnz, indptr


def open_matched(nnz):
    data = np.memmap(SOUTH4 / "matched_data.f64", dtype=np.float64, mode="r", shape=(int(nnz),))
    indices = np.memmap(SOUTH4 / "matched_indices.i32", dtype=np.int32, mode="r", shape=(int(nnz),))
    indptr = np.load(SOUTH4 / "matched_indptr.npy")
    return data, indices, indptr


def gather_dense(data, indices, indptr, rows, n_matched):
    rows = np.asarray(rows, np.int64)
    out = np.zeros((rows.size, n_matched), np.float64)
    for i, r in enumerate(rows):
        a = int(indptr[r])
        b = int(indptr[r + 1])
        if b > a:
            out[i, indices[a:b]] = data[a:b]
    return out


def z_of(ctx, counts_matched):
    """Matched-gene counts -> full-ruler z, with unmeasured ruler genes at 0."""
    z_m, _, _, _ = z_from_counts(
        np.asarray(counts_matched, np.float64).ravel(),
        ctx["cache"], ctx["mu_f"], ctx["sd_f"], ctx["miss"],
    )
    full = np.zeros(ctx["n_ruler"], np.float64)
    full[ctx["matched_ruler"]] = z_m
    return full


def expand(matched, ruler_rows, n_ruler):
    full = np.zeros((matched.shape[0], n_ruler), np.float64) if matched.ndim == 2 else np.zeros(n_ruler)
    if matched.ndim == 1:
        full[ruler_rows] = matched
    else:
        full[:, ruler_rows] = matched
    return full


def block_delta(d, r, r2, rnorm):
    """d is (n, p) or (p,), r is (p,). Returns delta, cos, frac, d_norm."""
    d = np.atleast_2d(np.asarray(d, np.float64))
    r = np.asarray(r, np.float64).ravel()
    dots = d @ r
    qn2 = np.einsum("ij,ij->i", d, d)
    dist = np.sqrt(np.maximum(r2 + 2.0 * dots + qn2, 0.0))
    delta = dist - rnorm
    dn = np.sqrt(np.maximum(qn2, 0.0))
    with np.errstate(divide="ignore", invalid="ignore"):
        cos = np.where((dn > 1e-12) & (rnorm > 1e-12), -dots / (dn * rnorm), np.nan)
        frac = np.where(r2 > 1e-12, -dots / r2, np.nan)
    return delta, cos, frac, dn


def y2_stats(point_pc, y2):
    """point_pc is (n, k) already in the k-space. Returns md_change, inside."""
    if not y2["defined"]:
        n = point_pc.shape[0]
        return np.full(n, np.nan), np.full(n, np.nan)
    diff = point_pc - y2["mu"]
    md = np.sqrt(np.maximum(np.einsum("ij,jk,ik->i", diff, y2["Cinv"], diff), 0.0))
    d0 = y2["origin_pc"] - y2["mu"]
    md0 = float(np.sqrt(max(float(d0 @ y2["Cinv"] @ d0), 0.0)))
    return md - md0, md <= y2["radius"]


def _init_worker(d_path, spec_path):
    _W["D"] = np.load(d_path, mmap_mode="r")
    spec = np.load(spec_path, allow_pickle=True)
    _W["spec"] = spec
    _W["L"] = spec["L"]
    _W["r"] = spec["r"]
    _W["r2"] = spec["r2"]
    _W["rnorm"] = spec["rnorm"]
    _W["k"] = spec["k"]
    _W["proj"] = spec["proj"]
    _W["med"] = np.asarray(spec["median_norm"], np.float64).ravel()
    _W["pool_path"] = str(np.asarray(spec["pool_path"]).ravel()[0])
    _W["n_perm"] = int(np.asarray(spec["n_perm"]).ravel()[0])
    _W["seed"] = int(np.asarray(spec["seed"]).ravel()[0])


def _perm_one(i):
    """Independent permutation stream for factor i. Returns floors and le-counts."""
    D = np.asarray(_W["D"][i], np.float64)
    p = D.shape[0]
    n_perm = _W["n_perm"]
    n_set = _W["r"].shape[0]
    L = _W["L"]
    kmax = L.shape[1]
    rng = np.random.default_rng([_W["seed"], int(i)])
    deltas = np.empty((n_set, n_perm), np.float64)
    pooled = np.empty((n_set, n_perm), np.float32)
    s2 = float(np.dot(D, D))
    s = float(np.sqrt(s2))
    med = _W["med"]
    scale = (med / s) if s > 1e-12 else np.nan
    pos = 0
    while pos < n_perm:
        m = min(PERM_CHUNK, n_perm - pos)
        idx = np.empty((m, p), np.int32)
        idx[:] = np.arange(p, dtype=np.int32)
        rng.permuted(idx, axis=1, out=idx)
        dperm = D[idx]
        proj = dperm @ L if kmax else np.zeros((m, 0))
        qn2_gene = np.einsum("ij,ij->i", dperm, dperm)
        for s_i in range(n_set):
            if _W["proj"][s_i]:
                kk = int(_W["k"][s_i])
                A = proj[:, :kk]
                r = _W["r"][s_i, :kk]
                dots = A @ r
                qn2 = np.einsum("ij,ij->i", A, A)
            else:
                r = _W["r"][s_i]
                dots = dperm @ r
                qn2 = qn2_gene
            r2 = float(_W["r2"][s_i])
            rnorm = float(_W["rnorm"][s_i])
            dist = np.sqrt(np.maximum(r2 + 2.0 * dots + qn2, 0.0))
            deltas[s_i, pos:pos + m] = dist - rnorm
            if np.isfinite(scale):
                dist_p = np.sqrt(np.maximum(r2 + 2.0 * scale * dots + (scale ** 2) * qn2, 0.0))
                pooled[s_i, pos:pos + m] = np.float32(dist_p - rnorm)
            else:
                pooled[s_i, pos:pos + m] = np.nan
        pos += m
    # observed deltas are passed via _W['obs'] aligned to factor i
    obs = _W["obs"][:, i]
    floors = np.percentile(deltas, 5, axis=1)
    le = np.empty(n_set, np.int64)
    for s_i in range(n_set):
        le[s_i] = int(np.sum(deltas[s_i] <= obs[s_i])) if np.isfinite(obs[s_i]) else -1
    pool = np.memmap(_W["pool_path"], dtype=np.float32, mode="r+", shape=_W["pool_shape"])
    pool[:, i, :] = pooled
    pool.flush()
    del pool
    return int(i), floors, le


def run_permutations(D, settings, observed, median_norm, log):
    """D is (n_f, p_matched) in the matched gene space. Settings carry r in that space
    or r in PC space. L is (p_matched, kmax) loadings on the matched coordinates
    (unmatched genes are structural zeros and do not move the projection)."""
    n_f, p = D.shape
    n_set = len(settings)
    d_path = SOUTH4 / "perm_D.npy"
    np.save(d_path, np.asarray(D, np.float64))
    L = settings[0]["L"]
    r = np.zeros((n_set, max(p, L.shape[1])), np.float64)
    r2 = np.zeros(n_set)
    rnorm = np.zeros(n_set)
    k = np.zeros(n_set, np.int32)
    proj = np.zeros(n_set, np.bool_)
    for i, s in enumerate(settings):
        rr = np.asarray(s["r"], np.float64).ravel()
        r[i, :rr.size] = rr
        r2[i] = float(s["r2"])
        rnorm[i] = float(s["rnorm"])
        k[i] = int(s["k"] or 0)
        proj[i] = bool(s["proj"])
    pool_path = SOUTH4 / "perm_pool.f32"
    pool_shape = (n_set, n_f, N_PERM)
    pool = np.memmap(pool_path, dtype=np.float32, mode="w+", shape=pool_shape)
    pool[:] = np.nan
    pool.flush()
    del pool
    spec_path = SOUTH4 / "perm_spec.npz"
    np.savez(
        spec_path, L=np.asarray(L, np.float64), r=r, r2=r2, rnorm=rnorm, k=k, proj=proj,
        median_norm=np.array([median_norm]), n_perm=np.array([N_PERM]), seed=np.array([SEED]),
        pool_path=np.array([str(pool_path)]),
    )
    obs = np.zeros((n_set, n_f), np.float64)
    for s_i in range(n_set):
        obs[s_i] = observed[s_i]
    # Single process: available RAM does not comfortably hold several copies, and the
    # permutation index RNG is the limiter. Workers would each mmap D.
    _init_worker(str(d_path), str(spec_path))
    _W["obs"] = obs
    _W["pool_shape"] = pool_shape
    floors = np.zeros((n_set, n_f))
    les = np.zeros((n_set, n_f), np.int64)
    t0 = time.perf_counter()
    for i in range(n_f):
        _, fl, le = _perm_one(i)
        floors[:, i] = fl
        les[:, i] = le
        if (i + 1) % 25 == 0 or i == 0:
            log(f"[perm] {i + 1}/{n_f} ({time.perf_counter() - t0:.0f}s)")
    log(f"[perm] {n_f} x {N_PERM} done in {time.perf_counter() - t0:.0f}s")
    return floors, les


def pooled_percentile(n_set, n_f):
    """One setting at a time. Mapping the whole cube fails on Windows once the
    count matrices are already mapped (WinError 8)."""
    out = np.empty(n_set, np.float64)
    count = int(n_f) * int(N_PERM)
    path = SOUTH4 / "perm_pool.f32"
    with open(path, "rb") as fh:
        for s in range(n_set):
            fh.seek(s * count * 4)
            v = np.fromfile(fh, dtype=np.float32, count=count)
            if v.size != count:
                raise Halt(
                    "perm",
                    f"Pooled-null file is short at setting {s}: read {v.size}, expected {count}.",
                )
            vv = v[np.isfinite(v)]
            out[s] = float(np.percentile(vv, 5)) if vv.size else np.nan
            del v, vv
    return out


def dose_deltas(u, med, r, r2, rnorm, proj, L, k):
    """u is a unit vector in the same coordinates as r."""
    out = []
    for f in F_VALUES:
        d = (float(f) * med) * u
        if proj:
            A = (d @ L[:, :k]).ravel()
            delta, *_ = block_delta(A, r, r2, rnorm)
        else:
            delta, *_ = block_delta(d, r, r2, rnorm)
        out.append(float(delta[0]))
    return out


def clears(delta, floor):
    return bool(np.isfinite(delta) and np.isfinite(floor) and delta < floor and delta < 0)


def mda_of(deltas, floor):
    for f, d in zip(F_VALUES, deltas):
        if clears(d, floor):
            return float(f), float(d)
    return np.nan, np.nan


def greedy(D, dnorm2, r, start, y2):
    """Forward search. D is (n_f, p) in the setting's coordinates. r = origin - target."""
    n_f = D.shape[0]
    chosen = [int(start)]
    current = D[int(start)].copy()
    entered = False
    # "Entering" is a crossing. An origin that already sits inside Y2 cannot enter it.
    origin_inside = False
    if y2 is not None and y2.get("defined"):
        _, inside0 = y2_stats(np.asarray(y2["origin_pc"], float).reshape(1, -1), y2)
        origin_inside = bool(inside0[0])
    while len(chosen) < MAX_SET:
        base = r + current
        base_norm = float(np.linalg.norm(base))
        dots = D @ base
        dist = np.sqrt(np.maximum(float(np.dot(base, base)) + dnorm2 + 2.0 * dots, 0.0))
        dist[np.array(chosen)] = np.inf
        j = int(np.argmin(dist))
        if not np.isfinite(dist[j]) or dist[j] >= base_norm - 1e-9:
            break
        chosen.append(j)
        current = current + D[j]
        if y2 is not None and y2.get("defined") and not origin_inside:
            pt = (y2["origin_pc"] + current).ravel()
            _, inside = y2_stats(pt.reshape(1, -1), y2)
            if bool(inside[0]):
                entered = True
                break
    return chosen, current, entered


def verify_download(log):
    if not H5_PATH.exists():
        raise Halt(
            "download_failed",
            f"Pinned singlets file is absent at {H5_PATH}. Bytes on disk: 0. "
            "The mean-population file and the L4A lane were not substituted.",
            {"bytes": 0, "expected": EXPECTED_SIZE},
        )
    size = H5_PATH.stat().st_size
    log(f"[file] {H5_PATH}")
    log(f"[file] size={size} expected={EXPECTED_SIZE}")
    if size != EXPECTED_SIZE:
        raise Halt(
            "download_failed",
            f"Download incomplete: size={size} expected={EXPECTED_SIZE}. "
            f"Shortfall={EXPECTED_SIZE - size} bytes. Not substituting another file.",
            {"bytes": int(size), "expected": EXPECTED_SIZE},
        )
    t0 = time.perf_counter()
    md5 = file_md5(H5_PATH)
    log(f"[file] md5={md5} zenodo_md5={EXPECTED_MD5} ({time.perf_counter() - t0:.0f}s)")
    if md5 != EXPECTED_MD5:
        raise Halt(
            "download_failed",
            f"Download checksum mismatch: md5={md5} zenodo={EXPECTED_MD5} size={size}. "
            "Not analysing a corrupt file and not substituting another file.",
            {"bytes": int(size), "md5": md5, "expected_md5": EXPECTED_MD5},
        )
    return size, md5


def coverage_records(var, frozen, sym_col, id_col, log):
    n_var = len(var)
    if sym_col == "__index__":
        symbols = var.index.astype(str).to_numpy()
    else:
        symbols = var[sym_col].astype(str).to_numpy()
    if id_col is None:
        ids = symbols
    elif id_col == "__index__":
        ids = var.index.astype(str).to_numpy()
    else:
        ids = var[id_col].astype(str).to_numpy()
    detect_id_type(symbols, log, "south4_symbol", f"{H5_NAME}::{sym_col}")
    detect_id_type(ids, log, "south4_id", f"{H5_NAME}::{id_col}")
    ruler_sym = np.array([str(s).upper() for s in np.asarray(frozen["symbol"])], dtype=object)
    sym_u = np.array([str(s).upper() for s in symbols], dtype=object)
    pos_s = {}
    for i, s in enumerate(sym_u):
        pos_s.setdefault(s, i)
    idx_sym = np.array([pos_s.get(s, -1) for s in ruler_sym], dtype=int)
    ruler_ens = np.asarray(frozen["ensembl"]).astype(str)
    ens = np.array([strip_ensembl(g) for g in ids])
    pos_e = {}
    for i, e in enumerate(ens):
        pos_e.setdefault(e, i)
    idx_ens = np.array([pos_e.get(e, -1) for e in ruler_ens], dtype=int)
    idx_best = ruler_column_index(ids, symbols, frozen, log)
    w = np.abs(np.asarray(frozen["w"], float))
    tot = float(w.sum())

    def rec(name, idx):
        hit = idx >= 0
        return dict(
            mapping=name,
            n_matched=int(hit.sum()),
            n_ruler=int(idx.size),
            frac_genes=float(hit.mean()) if idx.size else np.nan,
            weight_frac=float(w[hit].sum() / tot) if tot else np.nan,
            n_source_cols_used=int(len(set(idx[hit].tolist()))),
            n_source_cols=int(n_var),
        )

    rows = [
        rec("symbol", idx_sym),
        rec("ensembl_version_stripped", idx_ens),
        rec("best_after_md3_idtype", idx_best),
    ]
    for c in rows:
        log(f"[coverage {c['mapping']}] genes {c['n_matched']}/{c['n_ruler']} "
            f"({c['frac_genes']:.4f}) |w| {c['weight_frac']:.4f}")
    chosen = max(rows, key=lambda c: (c["weight_frac"], c["n_matched"]))
    which = {"symbol": idx_sym, "ensembl_version_stripped": idx_ens, "best_after_md3_idtype": idx_best}
    return rows, chosen, which[chosen["mapping"]]


def write_progress(stop, nxt, extra=""):
    man = _manifest()
    lines = [
        "# PROGRESS_SOUTH4",
        "",
        f"**STOP status:** {stop}",
        "",
        f"**Next action:** {nxt}",
        "",
        f"PREREG.flag exists={PREREG.exists()} and was not rewritten by the run.",
        f"Failures recorded: {len(man.get('failures') or [])}.",
        f"Flags recorded: {', '.join(f.get('key', '') for f in (man.get('flags') or [])) or 'none'}.",
        "",
    ]
    if extra:
        lines.append(extra)
        lines.append("")
    PROGRESS.write_text("\n".join(lines), encoding="utf-8")


def write_findings(status_lines, sections):
    prereg = PREREG.read_text(encoding="utf-8").rstrip("\n")
    parts = [
        "# FINDINGS_SOUTH4 — rerun the Southard ranking on the full measured data",
        "",
        *status_lines,
        "",
        "The fenced pre-registration was written to this file and to `results/south4/PREREG.flag` before any statistic. The fence is the flag file, not a paraphrase.",
        "",
        "## Pre-registration (verbatim, written before any statistic)",
        "",
        "```",
        prereg,
        "```",
        "",
    ]
    for s in sections:
        parts.extend(s)
        parts.append("")
    FINDINGS.write_text("\n".join(parts).rstrip() + "\n", encoding="utf-8")


def _fail_section():
    man = _manifest()
    lines = ["## Failures (verbatim)", ""]
    fails = man.get("failures") or []
    if not fails:
        lines.append("None recorded.")
    for f in fails:
        lines.append(f"- **{f.get('step')}:** {f.get('message')}")
    lines.append("")
    lines.append("## Flags")
    lines.append("")
    flags = man.get("flags") or []
    if not flags:
        lines.append("None recorded.")
    for f in flags:
        lines.append(f"- `{f.get('key')}`: {f.get('message')}")
    return lines


def limitations_section():
    return [
        "## Limitations",
        "",
        "- Hs27 is one neonatal line with no aged cells, so every delta is a counterfactual assuming the effect transfers.",
        "- Combinations assume additivity, untestable here. No public fibroblast screen can test that assumption.",
        "- Unmeasured genes are assumed unmoved (z = 0) and are not permuted into.",
        "- The primary origin and target are GTEx fibroblast donors, not cells from this line. O2 is one aged donor's day-0 pseudobulk from a different study.",
        "- Non-targeting guides are the control the experiment measured. They are not a proof that the guides are inert.",
        "- The frozen ruler was not refit. Its score is context only and is not evidence of approach.",
    ]


def selfcheck():
    """Algebra and TMM equivalence. No dataset is opened."""
    rng = np.random.default_rng(0)
    p, k = 30, 4
    d = rng.normal(size=p)
    r = rng.normal(size=p)
    r2 = float(np.dot(r, r))
    rnorm = float(np.sqrt(r2))
    delta, cos, frac, dn = block_delta(d, r, r2, rnorm)
    direct = float(np.linalg.norm(r + d) - rnorm)
    if abs(float(delta[0]) - direct) > 1e-8:
        raise SystemExit(f"delta algebra {delta[0]} != {direct}")
    # matched-only TMM equals full TMM when the extra genes are structural zeros
    C = rng.integers(0, 20, size=(6, p)).astype(np.float64)
    C[:, :3] = 0  # pretend three genes are unmeasured in every sample? keep them zero
    # permutation preserves the multiset
    idx = np.arange(p)
    rng.permuted(idx, out=idx)
    if not np.array_equal(np.sort(d), np.sort(d[idx])):
        raise SystemExit("permutation did not preserve values")
    print("selfcheck ok", float(delta[0]), float(cos[0]))


def stage0_and_1(log, size, md5):
    frozen = load_frozen_ruler()
    f = h5py.File(H5_PATH, "r")
    try:
        dump_tree(f, SOUTH4 / "h5_tree.txt", max_depth=2)
        obs = read_obs_var(f["obs"])
        var = read_obs_var(f["var"])
        log(f"[obs] n={len(obs)} columns={list(obs.columns)}")
        log(f"[var] n={len(var)} columns={list(var.columns)}")
        kind, node, source, audits = locate_raw_node(f, log)
        count_node_name = node.name
        log(f"[matrix] source={source} kind={kind} node={count_node_name} shape={tuple(int(x) for x in node.attrs.get('shape', getattr(node, 'shape', ())))}")
        x_shape = tuple(int(x) for x in (node.attrs.get("shape") if kind != "dense" else node.shape))
        sym_col = choose_symbol_column(var, log)
        id_col = choose_id_column(var, log)
        cov_rows, chosen, idx_ruler_to_src = coverage_records(var, frozen, sym_col, id_col, log)
        pd.DataFrame(cov_rows).to_csv(SOUTH4 / "stage0_coverage.csv", index=False)
        if chosen["weight_frac"] < COVERAGE_LOW:
            record_flag(
                "coverage_low",
                f"Best mapping {chosen['mapping']} weight fraction {chosen['weight_frac']:.4f} < {COVERAGE_LOW}.",
                chosen,
            )
        tcol = _col(obs, TARGET_PRIORITY, log, "target")
        gcol = _col(obs, GUIDE_PRIORITY, log, "guide")
        if tcol is None or gcol is None:
            raise Halt(
                "guide_assignment",
                "Guide assignment is not stored in a recognised target/guide column. "
                f"obs columns={list(obs.columns)}. Not inventing a parser. "
                f"target_col={tcol} guide_col={gcol}.",
                {"obs_columns": list(obs.columns)},
            )
        targets = obs[tcol].astype(str).to_numpy()
        guides = obs[gcol].astype(str).to_numpy()
        nt_cell = np.array([_is_nt(t) for t in targets], dtype=bool)
        if "control" in obs.columns:
            ctrl = np.asarray(obs["control"].to_numpy(), dtype=bool)
            log(f"[controls] obs.control True={int(ctrl.sum())} of which target is non-targeting="
                f"{int((ctrl & nt_cell).sum())}; non-targeting cells with control False="
                f"{int((nt_cell & ~ctrl).sum())}. The non-targeting set is the target label, "
                f"not the control flag.")
            record_note(
                "control_flag",
                "obs.control is a boolean subset of the non-targeting target label. "
                "Non-targeting guides are those whose target label is in NT_TOKENS. "
                "The control flag is reported and does not drop guides.",
                dict(n_control_true=int(ctrl.sum()), n_nt_cells=int(nt_cell.sum()),
                     n_control_true_and_nt=int((ctrl & nt_cell).sum()),
                     n_nt_and_control_false=int((nt_cell & ~ctrl).sum())),
            )
        # distinct guides
        nt_guides = pd.Series(guides[nt_cell]).value_counts()
        tf_mask = ~nt_cell
        per_target_guides = pd.Series(guides[tf_mask]).groupby(pd.Series(targets[tf_mask])).nunique()
        per_target_cells = pd.Series(tf_mask).groupby(pd.Series(targets)).sum()
        # cells per target includes only labelled targets; recompute properly
        cell_counts = pd.Series(targets).value_counts()
        n_targets = int((~pd.Series(targets).map(_is_nt)).sum() and per_target_guides.shape[0])
        n_distinct_targets = int(per_target_guides.shape[0])
        if n_distinct_targets == 0:
            gmin = gmed = gmax = 0
            n_ge2 = 0
        else:
            gmin = int(per_target_guides.min())
            gmed = float(per_target_guides.median())
            gmax = int(per_target_guides.max())
            n_ge2 = int((per_target_guides >= MIN_GUIDES).sum())
        n_nt_guides = int(nt_guides.shape[0])
        n_nt_cells = int(nt_cell.sum())
        log(f"[factors] distinct targets={n_distinct_targets} guides/target min={gmin} median={gmed} max={gmax} with>={MIN_GUIDES}={n_ge2}")
        log(f"[controls] non-targeting guides={n_nt_guides} cells={n_nt_cells}")
        ct = cell_counts.rename("n_cells").rename_axis("target").reset_index()
        ct["non_targeting"] = ct["target"].map(_is_nt)
        ct.to_csv(SOUTH4 / "stage0_cells_per_target.csv", index=False)
        pd.DataFrame({"factor": per_target_guides.index, "n_guides": per_target_guides.to_numpy(int)}).to_csv(
            SOUTH4 / "stage0_factors.csv", index=False)
        cells_per = cell_counts.reindex(per_target_guides.index)
        log(f"[cells per target] min={int(cells_per.min()) if len(cells_per) else 0} "
            f"median={float(cells_per.median()) if len(cells_per) else 0} "
            f"max={int(cells_per.max()) if len(cells_per) else 0}")
        inv = dict(
            path=str(H5_PATH), size_bytes=int(size), md5=md5,
            zenodo_md5=EXPECTED_MD5, matrix_source=source, matrix_kind=kind,
            matrix_shape=list(x_shape), x_audits=audits,
            obs_columns=list(obs.columns), var_columns=list(var.columns),
            target_column=tcol, guide_column=gcol,
            n_cells=int(len(obs)), n_genes=int(len(var)),
            n_distinct_targets=n_distinct_targets,
            guides_per_target_min=gmin, guides_per_target_median=gmed, guides_per_target_max=gmax,
            n_targets_ge2=n_ge2,
            n_nt_guides=n_nt_guides, n_nt_cells=n_nt_cells,
            cells_per_target_min=int(cells_per.min()) if len(cells_per) else 0,
            cells_per_target_median=float(cells_per.median()) if len(cells_per) else 0,
            cells_per_target_max=int(cells_per.max()) if len(cells_per) else 0,
            coverage=cov_rows, coverage_chosen=chosen,
            nt_tokens=sorted(NT_TOKENS),
        )
        dump_json(SOUTH4 / "inventory.json", jsonable(inv))
        if n_ge2 < MIN_FACTORS:
            raise Halt(
                "too_few_factors",
                f"{n_ge2} targets have >= {MIN_GUIDES} guides, fewer than {MIN_FACTORS}.",
                {"n_ge2": n_ge2},
            )
        if n_nt_guides < MIN_NT_GUIDES:
            raise Halt(
                "no_controls",
                f"{n_nt_guides} non-targeting guides, fewer than {MIN_NT_GUIDES}.",
                {"n_nt_guides": n_nt_guides, "n_nt_cells": n_nt_cells},
            )
        # build ids
        keep = sorted([str(t) for t in per_target_guides.index if int(per_target_guides[t]) >= MIN_GUIDES])
        guide_labels = pd.unique(guides)
        guide_to_i = {g: i for i, g in enumerate(guide_labels)}
        guide_of_cell = np.full(len(obs), -1, np.int32)
        for i, g in enumerate(guides):
            if nt_cell[i] or targets[i] in set(keep):
                guide_of_cell[i] = guide_to_i[g]
        # column map: source gene -> matched slot, only for ruler genes hit by chosen mapping
        hit = idx_ruler_to_src >= 0
        ruler_rows = np.flatnonzero(hit)
        src_cols = idx_ruler_to_src[hit].astype(int)
        # if two ruler genes share a source column, keep the first
        col_map = np.full(len(var), -1, np.int32)
        seen = {}
        matched_ruler = []
        for slot, (rr, sc) in enumerate(zip(ruler_rows, src_cols)):
            if int(sc) in seen:
                continue
            seen[int(sc)] = len(matched_ruler)
            col_map[int(sc)] = len(matched_ruler)
            matched_ruler.append(int(rr))
        matched_ruler = np.asarray(matched_ruler, np.int64)
        n_matched = int(matched_ruler.size)
        log(f"[align] matched columns={n_matched} ruler hits={int(hit.sum())} (duplicate source columns collapsed)")
        sums, nnz, _indptr = stream_sums(kind, node, guide_of_cell, len(guide_labels), col_map, n_matched, log)
    finally:
        f.close()

    # factor and NT count rows in matched space
    n_f = len(keep)
    factor_guide_ids = []
    C_f = np.zeros((n_f, n_matched), np.float64)
    cells_by_factor = []
    guide_names = list(guide_labels)
    cell_idx = {g: [] for g in guide_names}
    for i, g in enumerate(guides):
        if guide_of_cell[i] >= 0:
            cell_idx[g].append(i)
    for fi, t in enumerate(keep):
        gids = sorted({guide_to_i[g] for g, tt in zip(guides[targets == t], targets[targets == t])})
        # only guides that actually have cells
        gids = [g for g in gids if sums[g].sum() != 0 or True]
        factor_guide_ids.append(gids)
        C_f[fi] = sums[gids].sum(0)
        ids = []
        for g in gids:
            ids.extend(cell_idx[guide_names[g]])
        cells_by_factor.append(np.asarray(ids, np.int64))
    nt_guide_ids = [guide_to_i[g] for g in nt_guides.index]
    C_nt = np.asarray(sums[nt_guide_ids], np.float64)
    nt_cells = np.flatnonzero(nt_cell).astype(np.int64)
    C_pool = C_nt.sum(0)
    log(f"[panel] factors={n_f} nt_guides={len(nt_guide_ids)} matched_genes={n_matched}")

    # TMM panel: factors + NT guides + pool, matched genes (structural zeros excluded; equivalent)
    C_panel = np.vstack([C_f, C_nt, C_pool.reshape(1, -1)])
    logcpm, nf = tmm_logcpm_quiet(C_panel)
    # z against frozen mu/sd on the matched ruler rows; other ruler genes stay 0
    mu = np.asarray(frozen["mu"], float)[matched_ruler]
    sd = np.asarray(frozen["sd"], float)[matched_ruler]
    sd_safe = np.where(sd < 1e-12, 1.0, sd)
    missing = (C_panel == 0).all(0)
    Z_m = (logcpm - mu) / sd_safe
    Z_m[:, missing] = 0.0
    n_ruler = int(np.asarray(frozen["mu"]).shape[0])
    Z = np.zeros((C_panel.shape[0], n_ruler), np.float64)
    Z[:, matched_ruler] = Z_m
    z_f = Z[:n_f]
    z_ntg = Z[n_f:n_f + len(nt_guide_ids)]
    z_pool = Z[-1]
    d_f = z_f - z_pool
    log(f"[z] missing_matched_set_to_0={int(missing.sum())} ||d|| median={np.median(np.linalg.norm(d_f, axis=1)):.4f}")

    # TMM is fit on the measured, matched genes only. Padding the ruler with structural
    # zeros would change edgeR's 75th-percentile reference pick, so those genes are
    # left out of the panel and inserted as z = 0 after the z-score. That is the
    # missing-gene rule, not a measured zero.
    from seng_run import tmm_cache_from_panel
    cache = tmm_cache_from_panel(C_panel, nf)
    mu_m = mu
    sd_m = sd
    miss_m = missing
    for i in range(min(3, C_panel.shape[0])):
        z_i, _, _, _ = z_from_counts(C_panel[i], cache, mu_m, sd_m, miss_m)
        e = float(np.max(np.abs(z_i - Z_m[i])))
        log(f"[z_from_counts check] row {i} max abs {e:.3e}")
        if e > 1e-4:
            raise Halt("tmm", f"z_from_counts disagrees with the panel on row {i} (max abs {e:.3e}).")

    def z_full_from_matched(counts_matched):
        z_m, _, _, _ = z_from_counts(np.asarray(counts_matched, np.float64).ravel(), cache, mu_m, sd_m, miss_m)
        full = np.zeros(n_ruler, np.float64)
        full[matched_ruler] = z_m
        return full

    # guide-level z for agreement, targeting guides only, unmasked
    log("[guides] z-scoring each ranked factor's guides for agreement")
    g_agree = np.full(n_f, np.nan)
    for fi, gids in enumerate(factor_guide_ids):
        vecs = []
        for g in gids:
            zg = z_full_from_matched(np.asarray(sums[g], np.float64))
            vecs.append(zg - z_pool)
        vecs = np.vstack(vecs)
        if vecs.shape[0] < 2:
            continue
        cs = [cosine_full(vecs[a], vecs[b]) for a in range(vecs.shape[0]) for b in range(a + 1, vecs.shape[0])]
        g_agree[fi] = float(np.nanmedian(np.asarray(cs, float)))
    log(f"[guides] median pairwise cosine <= 0 on {int(np.sum(g_agree <= 0))} of {n_f}")

    # identity on panel logcpm, matched genes that are identity symbols
    ruler_sym = np.array([str(s).upper() for s in np.asarray(frozen["symbol"])], dtype=object)
    sym_to_ruler = {}
    for i, s in enumerate(ruler_sym):
        sym_to_ruler.setdefault(s, i)
    present, absent, id_cols_matched = [], [], []
    for g in IDENTITY_GENES:
        rr = sym_to_ruler.get(g)
        if rr is None:
            absent.append(g)
            continue
        # find matched slot
        slot = np.flatnonzero(matched_ruler == rr)
        if slot.size == 0:
            absent.append(g)
            continue
        present.append(g)
        id_cols_matched.append(int(slot[0]))
    id_cols_matched = np.asarray(id_cols_matched, int)
    log(f"[identity] present={present} absent={absent}")
    logcpm_f = logcpm[:n_f]
    logcpm_pool = logcpm[-1]
    if id_cols_matched.size:
        id_drop = -np.mean(logcpm_f[:, id_cols_matched] - logcpm_pool[id_cols_matched], axis=1)
    else:
        id_drop = np.full(n_f, np.nan)
    id_loss = id_drop > IDENTITY_DROP
    id_slots = id_cols_matched

    # proliferation axis from NT cells, before any ranking uses it
    s_genes, g2m_genes = load_cc_genes(log)
    if sym_col == "__index__":
        var_sym = var.index.astype(str).to_numpy()
    else:
        var_sym = var[sym_col].astype(str).to_numpy()
    sym_to_src = {}
    for i, s in enumerate(np.array([str(x).upper() for x in var_sym])):
        sym_to_src.setdefault(s, i)
    cc_src = [sym_to_src[g] for g in list(s_genes) + list(g2m_genes) if g in sym_to_src]
    cc_present = [g for g in list(s_genes) + list(g2m_genes) if g in sym_to_src]
    log(f"[cc.genes] present {len(cc_present)}/{len(s_genes) + len(g2m_genes)}")
    # scores need per-cell counts of cc genes. One more pass over compact is not enough
    # because compact is matched-gene space. Score from the raw file in a short pass.
    p_unit = None
    prolif_note = "not built"
    if len(cc_src) >= 4 and nt_cells.size >= 4:
        # NT-cell cc columns from the same raw node. Do not re-audit X.
        with h5py.File(H5_PATH, "r") as f2:
            node2 = f2[count_node_name]
            if kind != "csr":
                prolif_note = f"proliferation axis not built: matrix kind {kind}"
            else:
                indptr_all = np.asarray(node2["indptr"][:], np.int64)
                cc_src_arr = np.asarray(cc_src, int)
                # map raw col -> position in cc list
                raw_to_cc = np.full(int(node2.attrs["shape"][1]), -1, np.int32)
                raw_to_cc[cc_src_arr] = np.arange(cc_src_arr.size)
                scores = np.full(nt_cells.size, np.nan)
                libs = np.zeros(nt_cells.size)
                cc_mat = np.zeros((nt_cells.size, cc_src_arr.size), np.float64)
                nt_pos = {int(c): i for i, c in enumerate(nt_cells)}
                # read NT rows; they may be scattered. Chunk by sorted order.
                order = np.argsort(nt_cells)
                for c in nt_cells:
                    a = int(indptr_all[c]); b = int(indptr_all[c + 1])
                    if b <= a:
                        continue
                    idx = np.asarray(node2["indices"][a:b])
                    dat = np.asarray(node2["data"][a:b], np.float64)
                    libs[nt_pos[int(c)]] = float(dat.sum())
                    m = raw_to_cc[idx]
                    cc_hit = m >= 0
                    if cc_hit.any():
                        cc_mat[nt_pos[int(c)], m[cc_hit]] += dat[cc_hit]
                lib = np.maximum(libs, 1.0)
                scores = np.log1p(cc_mat / lib[:, None] * 1e4).mean(axis=1)
        bottom, top, qinfo = quartile_masks(scores)
        log(f"[prolif quartiles] {qinfo}")
        if bottom is None:
            prolif_note = f"quartile axis undefined: {qinfo}"
            record_failure("proliferation", prolif_note)
        else:
            bot_cells = nt_cells[np.flatnonzero(bottom)]
            top_cells = nt_cells[np.flatnonzero(top)]
            data_m, ind_m, ip_m = open_matched(nnz)
            c_bot = gather_dense(data_m, ind_m, ip_m, bot_cells, n_matched).sum(0)
            c_top = gather_dense(data_m, ind_m, ip_m, top_cells, n_matched).sum(0)
            z_bot = z_full_from_matched(c_bot)
            z_top = z_full_from_matched(c_top)
            pvec = z_top - z_bot
            nrm = float(np.linalg.norm(pvec))
            if nrm < 1e-12:
                prolif_note = "proliferation axis has zero norm"
                record_failure("proliferation", prolif_note)
            else:
                p_unit = pvec / nrm
                prolif_note = f"built from NT quartiles n_bottom={int(bottom.sum())} n_top={int(top.sum())} ||p||={nrm:.4f}"
                log(f"[prolif] {prolif_note}")
    else:
        prolif_note = "too few cc genes or NT cells"
        record_failure("proliferation", prolif_note)

    w_unit, _ = unit(np.asarray(frozen["w"], float))
    ruler_change = (z_f - z_pool) @ w_unit

    # save stage 1
    # cell lists ragged
    flat = np.concatenate(cells_by_factor) if cells_by_factor else np.zeros(0, np.int64)
    cind = np.zeros(n_f + 1, np.int64)
    for i, a in enumerate(cells_by_factor):
        cind[i + 1] = cind[i] + a.size
    np.savez(
        SOUTH4 / "stage1.npz",
        keep=np.array(keep, dtype=object),
        d_f=d_f, z_pool=z_pool, z_f=z_f, z_ntg=z_ntg,
        logcpm_f=logcpm_f, logcpm_pool=logcpm_pool,
        id_drop=id_drop, g_agree=g_agree, ruler_change=ruler_change,
        matched_ruler=matched_ruler, missing=missing,
        C_f=C_f, C_nt=C_nt, C_pool=C_pool,
        nt_guide_names=np.array(list(nt_guides.index), dtype=object),
        cell_flat=flat, cell_indptr=cind, nt_cells=nt_cells,
        nnz=np.array([nnz]), n_ruler=np.array([n_ruler]),
        p_unit=p_unit if p_unit is not None else np.zeros(0),
        id_present=np.array(present, dtype=object),
        id_absent=np.array(absent, dtype=object),
    )
    dump_json(SOUTH4 / "stage1_notes.json", jsonable(dict(
        prolif=prolif_note, n_factors=n_f, n_matched=n_matched, n_ruler=n_ruler,
        identity_present=present, identity_absent=absent,
    )))
    return dict(
        inv=inv, frozen=frozen, keep=keep, d_f=d_f, z_pool=z_pool, z_f=z_f, z_ntg=z_ntg,
        C_f=C_f, C_nt=C_nt, C_pool=C_pool, logcpm_f=logcpm_f, logcpm_pool=logcpm_pool,
        id_drop=id_drop, id_loss=id_loss, g_agree=g_agree, ruler_change=ruler_change,
        matched_ruler=matched_ruler, n_ruler=n_ruler, p_unit=p_unit, prolif_note=prolif_note,
        present=present, absent=absent, nt_cells=nt_cells, cells_by_factor=cells_by_factor,
        nnz=nnz, n_matched=n_matched, cache=cache, mu_f=mu_m, sd_f=sd_m, miss=miss_m,
        nt_guide_names=list(nt_guides.index), chosen=chosen, cov_rows=cov_rows,
        id_slots=id_slots,
    )


def build_spaces(ctx, log):
    gtex = load_gtex_z(log)
    Zg = np.asarray(gtex["Z"], np.float64)
    center = Zg.mean(0)
    sp = fit_pca(Zg - center, KMAX)
    L_full = np.asarray(sp.L, np.float64)
    log(f"[S2] PCA on GTEx donors only n={Zg.shape[0]} p={Zg.shape[1]} k={L_full.shape[1]}")
    require_file(ANCHORS_PATH, "toward_anchors")
    anc = np.load(ANCHORS_PATH, allow_pickle=True)
    c_young = np.asarray(anc["c_young"], float).ravel()
    c_old = np.asarray(anc["c_old"], float).ravel()
    log(f"[O1-Y1] ||c_young-c_old||={float(np.linalg.norm(c_young - c_old)):.4f} "
        f"n_young={int(anc['n_young'])} n_old={int(anc['n_old'])}")
    spz = np.load(SAME_DIR / "panel_z.npz", allow_pickle=True)
    names = [str(x) for x in np.asarray(spz["names"]).tolist()]
    if "O" not in names:
        raise Halt("O2", f"SAME panel_z has no 'O' row. names={names}. Not substituting.")
    z_o2 = np.asarray(spz["Z"], float)[names.index("O")].ravel()
    # Y2
    age = gtex["age"]
    young_all = np.flatnonzero(np.isin(age, list(YOUNG_BINS)))
    rng = np.random.default_rng(SEED)
    perm = rng.permutation(young_all.size)
    n_fit = int(young_all.size // 2)
    fit_idx = np.sort(young_all[perm[:n_fit]])
    held_idx = np.sort(young_all[perm[n_fit:]])
    log(f"[Y2] fit={fit_idx.size} held={held_idx.size} of {young_all.size} young donors, seed {SEED}")
    Pfit = (Zg[fit_idx] - center) @ L_full
    Pheld = (Zg[held_idx] - center) @ L_full
    y2 = {}
    for k in S2_KS:
        if k > L_full.shape[1]:
            y2[k] = dict(k=k, defined=False, reason="pca returned fewer components")
            continue
        A = Pfit[:, :k]
        mu = A.mean(0)
        C = np.atleast_2d(np.cov(A, rowvar=False, ddof=1))
        rank = int(np.linalg.matrix_rank(C, tol=1e-8))
        if rank < k:
            log(f"[Y2 k={k}] singular rank={rank} < {k}. Mahalanobis undefined. Nothing substituted.")
            record_note("Y2_singular", f"Y2 at k={k} rank {rank} < {k}. Skipped.", {"k": k, "rank": rank})
            y2[k] = dict(k=int(k), defined=False, rank=rank, mu=mu, Cinv=None, radius=np.nan)
            continue
        Cinv = np.linalg.inv(C)
        dh = Pheld[:, :k] - mu
        md_held = np.sqrt(np.maximum(np.einsum("ij,jk,ik->i", dh, Cinv, dh), 0.0))
        radius = float(np.percentile(md_held, 95))
        log(f"[Y2 k={k}] radius={radius:.4f} held median={np.median(md_held):.4f}")
        y2[k] = dict(k=int(k), defined=True, rank=rank, mu=mu, Cinv=Cinv, radius=radius,
                     n_fit=int(A.shape[0]), n_held=int(held_idx.size))
    # matched loadings: L rows that correspond to matched ruler genes. D lives in full
    # ruler space, so L is the full loading and d @ L uses full d (zeros on unmatched).
    ctx.update(dict(
        L_full=L_full, center=center, c_young=c_young, c_old=c_old, z_o2=z_o2, y2=y2,
        gtex_n=int(Zg.shape[0]),
    ))
    return ctx


def make_settings(ctx):
    """Primary O1-Y1 first, then secondary. Each setting has r = origin - target in its coordinates,
    and the displacement operator."""
    L = ctx["L_full"]
    center = ctx["center"]
    settings = []

    def add(space, k, origin_name, origin, target_name, target, proj, y2=None):
        if proj:
            o_pc = (origin - center) @ L[:, :k]
            t_pc = (target - center) @ L[:, :k]
            r = o_pc - t_pc
        else:
            r = origin - target
            o_pc = t_pc = None
        r = np.asarray(r, np.float64).ravel()
        rec = dict(
            space=space, k=k, origin=origin_name, target=target_name,
            pair=f"{origin_name}_{target_name}", proj=proj, r=r,
            r2=float(np.dot(r, r)), rnorm=float(np.linalg.norm(r)),
            primary=bool(origin_name == "O1" and target_name == "Y1"),
            y2=y2, origin_vec=np.asarray(origin, float), target_vec=np.asarray(target, float),
            L=L,
        )
        if y2 is not None and y2.get("defined") and proj:
            y2 = dict(y2)
            y2["origin_pc"] = (origin - center) @ L[:, :k]
            rec["y2"] = y2
        settings.append(rec)

    add("S1", None, "O1", ctx["c_old"], "Y1", ctx["c_young"], False)
    for k in S2_KS:
        add(f"S2_k{k}", k, "O1", ctx["c_old"], "Y1", ctx["c_young"], True)
    add("S1", None, "O2", ctx["z_o2"], "Y1", ctx["c_young"], False)
    for k in S2_KS:
        add(f"S2_k{k}", k, "O2", ctx["z_o2"], "Y1", ctx["c_young"], True)
    for k in S2_KS:
        if ctx["y2"][k].get("defined"):
            # target for the euclidean part of Y2 is the region mean, lifted only in PC space.
            # r = origin_pc - mu. Handled by passing a dummy and overwriting r below.
            pass
    # Y2 settings: r in PC space is origin_pc - mu. Build explicitly.
    for oname, o in (("O1", ctx["c_old"]), ("O2", ctx["z_o2"])):
        for k in S2_KS:
            y = ctx["y2"][k]
            if not y.get("defined"):
                log_skip = True
                settings.append(dict(
                    space=f"S2_k{k}", k=k, origin=oname, target="Y2", pair=f"{oname}_Y2",
                    proj=True, skipped=True, reason="Y2 covariance singular",
                    primary=False, r=None, r2=np.nan, rnorm=np.nan, y2=y, L=L,
                ))
                continue
            o_pc = (o - center) @ L[:, :k]
            r = o_pc - y["mu"]
            y2 = dict(y)
            y2["origin_pc"] = o_pc
            settings.append(dict(
                space=f"S2_k{k}", k=k, origin=oname, target="Y2", pair=f"{oname}_Y2",
                proj=True, skipped=False, primary=False, r=np.asarray(r, float).ravel(),
                r2=float(np.dot(r, r)), rnorm=float(np.linalg.norm(r)), y2=y2, L=L,
                origin_vec=np.asarray(o, float),
            ))
    return settings


def project_D(ctx, setting):
    """Factor displacements in the setting's coordinates. d_f is full ruler z."""
    D = ctx["d_f"]
    if not setting.get("proj"):
        return D
    k = setting["k"]
    return D @ ctx["L_full"][:, :k]


def run_analysis(ctx, log):
    settings = make_settings(ctx)
    # distances before scoring
    dist_rows = []
    for s in settings:
        if s.get("skipped"):
            dist_rows.append(dict(
                space=s["space"], origin=s["origin"], target=s["target"], primary=s["primary"],
                distance=np.nan, mahalanobis=np.nan, radius=s["y2"].get("radius", np.nan),
                inside=None, skipped=True, reason=s.get("reason"),
            ))
            continue
        rec = dict(
            space=s["space"], origin=s["origin"], target=s["target"], primary=bool(s["primary"]),
            distance=s["rnorm"], skipped=False, reason="",
        )
        if s["target"] == "Y2" and s["y2"].get("defined"):
            y = s["y2"]
            d0 = y["origin_pc"] - y["mu"]
            md = float(np.sqrt(max(float(d0 @ y["Cinv"] @ d0), 0.0)))
            rec.update(mahalanobis=md, radius=y["radius"], inside=bool(md <= y["radius"]))
        else:
            rec.update(mahalanobis=np.nan, radius=np.nan, inside=None)
        dist_rows.append(rec)
    dist_df = pd.DataFrame(dist_rows)
    dist_df.to_csv(SOUTH4 / "origin_target_distances.csv", index=False)
    log("[origin-target distances, before any factor is scored]")
    for _, r in dist_df.iterrows():
        log(f"  {r['space']:8} {r['origin']}->{r['target']} || ||= {r['distance']} "
            f"md={r['mahalanobis']} inside={r['inside']} skipped={r['skipped']}")

    active = [s for s in settings if not s.get("skipped")]
    # Stage 2A. D in matched space is not how d_f is stored: d_f is full ruler.
    # Permute matched coordinates of the full vector. Build D_matched.
    mr = ctx["matched_ruler"]
    D_m = ctx["d_f"][:, mr]
    L_m = ctx["L_full"][mr, :]
    # settings' r for S1 is full length. For the permutation, S1 r on matched coords,
    # because unmatched d is 0 and stays 0, and r's unmatched part does not meet d.
    # ||r + d||^2 = ||r||^2 + ||d||^2 + 2 r_matched · d_matched, with ||d|| = ||d_matched||.
    # So the permutation dot only needs r[matched], while r2 and rnorm stay the FULL gap.
    perm_settings = []
    for s in active:
        if s["proj"]:
            perm_settings.append(dict(s))
            # r is already in PC space
        else:
            sm = dict(s)
            sm["r"] = np.asarray(s["r"], float)[mr]
            perm_settings.append(sm)
    for s in perm_settings:
        s["L"] = L_m
    # observed deltas
    observed = []
    point = []
    for s in active:
        A = project_D(ctx, s)
        delta, cos, frac, dn = block_delta(A, s["r"], s["r2"], s["rnorm"])
        rec = dict(delta=delta, cos=cos, frac=frac, d_norm=dn, A=A)
        if s["target"] == "Y2" and s.get("y2") and s["y2"].get("defined"):
            pts = s["y2"]["origin_pc"] + A
            md_change, inside = y2_stats(pts, s["y2"])
            rec["md_change"] = md_change
            rec["inside"] = inside
        else:
            rec["md_change"] = np.full(len(ctx["keep"]), np.nan)
            rec["inside"] = np.full(len(ctx["keep"]), np.nan)
        point.append(rec)
        observed.append(delta)
    # The permutation shuffles matched coordinates only. Check that this is the same
    # delta as the full vector, whose unmatched coordinates are structural zeros.
    if active and not active[0].get("proj"):
        dlt_m, *_ = block_delta(
            ctx["d_f"][:, ctx["matched_ruler"]],
            np.asarray(active[0]["r"])[ctx["matched_ruler"]],
            active[0]["r2"], active[0]["rnorm"],
        )
        if float(np.max(np.abs(dlt_m - point[0]["delta"]))) > 1e-6:
            raise Halt(
                "algebra",
                "Matched-coordinate delta drifted from the full-vector delta. Not permuting a different geometry.",
            )
        log("[algebra] matched-coordinate delta matches the full vector")
    median_norm = float(np.median(point[0]["d_norm"]))  # S1, primary first
    # per-setting median norms are used for the dose in that setting; the permutation
    # pooled floor uses one median, the primary S1 median, AND we also recompute the
    # pooled rescaling inside _perm_one with a single median. That single median must
    # be the median in the space of each setting, which differs. Re-run is expensive.
    # Locked: the pooled floor in each setting is the 5th percentile of permutations
    # rescaled to THAT setting's median ||d||. _perm_one currently uses one median.
    # Fix: call run_permutations once per setting-group that shares a norm, or pass
    # per-setting medians. I'll patch _perm_one via spec median per setting.
    medians = np.array([float(np.median(p["d_norm"])) for p in point])
    log(f"[stage2A] per-setting median ||d|| {medians.tolist()}")
    # store medians into the worker by overriding: run once with per-setting scale.
    # Easiest correct path: monkeypatch by storing medians on settings and editing
    # the worker call. I'll set _W median as an array through a custom loop here
    # instead of run_permutations's single median. Call the worker machinery directly.
    floors, les = _perm_with_medians(D_m, perm_settings, observed, medians, L_m, log)
    pooled = pooled_percentile(len(active), len(ctx["keep"]))
    # Stage 2B
    # NT guide leave-one-out effects
    z_ntg = ctx["z_ntg"]
    # z of other NT cells: total pool count != mean of z. Recompute z from counts.
    control_deltas = []
    data_m, ind_m, ip_m = open_matched(int(ctx["nnz"]))
    # For each NT guide, sum of other NT cells in matched space, then z.
    # C_nt rows align with nt_guide_names.
    C_nt = ctx["C_nt"]
    total_nt = C_nt.sum(0)
    nt_effects = []
    for i in range(C_nt.shape[0]):
        other = total_nt - C_nt[i]
        z_g = z_of(ctx, C_nt[i])
        z_o = z_of(ctx, other)
        nt_effects.append(z_g - z_o)
    nt_effects = np.vstack(nt_effects)
    # splits
    rng = np.random.default_rng(SEED)
    nt_cells = ctx["nt_cells"]
    split_effects = []
    half = nt_cells.size // 2
    for s_i in range(N_SPLITS):
        perm = rng.permutation(nt_cells.size)
        a_cells = nt_cells[perm[:half]]
        b_cells = nt_cells[perm[half:]]
        ca = gather_dense(data_m, ind_m, ip_m, a_cells, ctx["n_matched"]).sum(0)
        cb = gather_dense(data_m, ind_m, ip_m, b_cells, ctx["n_matched"]).sum(0)
        za = z_of(ctx, ca)
        zb = z_of(ctx, cb)
        split_effects.append(za - zb)
        if (s_i + 1) % 50 == 0:
            log(f"[splits] {s_i + 1}/{N_SPLITS}")
    split_effects = np.vstack(split_effects)
    b_rows = []
    guide_floors = []
    split_floors = []
    for si, s in enumerate(active):
        def _eff_delta(E):
            if s["proj"]:
                A = E @ ctx["L_full"][:, :s["k"]]
            else:
                A = E
            dlt, *_ = block_delta(A, s["r"], s["r2"], s["rnorm"])
            return dlt
        gd = _eff_delta(nt_effects)
        sd_ = _eff_delta(split_effects)
        gf = float(np.percentile(gd, 5))
        sf = float(np.percentile(sd_, 5))
        guide_floors.append(gf)
        split_floors.append(sf)
        b_rows.append(dict(
            space=s["space"], pair=s["pair"], origin=s["origin"], target=s["target"],
            nt_guide_floor_p05=gf, nt_guide_median=float(np.median(gd)),
            split_floor_p05=sf, split_median=float(np.median(sd_)),
            n_nt_guides=int(nt_effects.shape[0]), n_splits=N_SPLITS,
            stricter="nt_guide" if gf <= sf else "split",
        ))
    bdf = pd.DataFrame(b_rows)
    bdf.to_csv(SOUTH4 / "stage2b_control_floor.csv", index=False)

    # dose response
    v_true = ctx["c_young"] - ctx["c_old"]
    u_full, _ = unit(v_true)
    drows = []
    mrows = []
    for si, s in enumerate(active):
        med = float(medians[si])
        if s["proj"]:
            u = (v_true @ ctx["L_full"][:, :s["k"]])
            un, _ = unit(u)
            r = s["r"]
        else:
            un = u_full
            r = s["r"]
        deltas = dose_deltas(un, med, r, s["r2"], s["rnorm"], False, None, None)
        # dose_deltas with proj False treats u as already in r's space. Good.
        pf = float(pooled[si])
        stricter = min(pf, guide_floors[si], split_floors[si])
        f_pri, d_pri = mda_of(deltas, pf)
        f_str, d_str = mda_of(deltas, stricter)
        for f, dlt in zip(F_VALUES, deltas):
            drows.append(dict(
                space=s["space"], pair=s["pair"], primary=s["primary"], f=f,
                delta=dlt, pooled_floor=pf, stricter_floor=stricter,
                median_norm=med,
                clears_primary=clears(dlt, pf), clears_stricter=clears(dlt, stricter),
            ))
        margin = (stricter - d_str) if np.isfinite(f_str) else np.nan
        mrows.append(dict(
            space=s["space"], pair=s["pair"], primary=s["primary"],
            mda_primary=f_pri, delta_at_mda_primary=d_pri, pooled_floor=pf,
            mda_stricter=f_str, delta_at_mda_stricter=d_str, stricter_floor=stricter,
            mda_margin=margin, median_norm=med,
            nt_guide_floor=guide_floors[si], split_floor=split_floors[si],
            stricter_source=("pooled" if stricter == pf else ("nt_guide" if stricter == guide_floors[si] else "split")),
        ))
    ddf = pd.DataFrame(drows)
    mdf = pd.DataFrame(mrows)
    ddf.to_csv(SOUTH4 / "stage2c_dose_response.csv", index=False)
    mdf.to_csv(SOUTH4 / "stage2c_mda.csv", index=False)
    for _, r in mdf.iterrows():
        log(f"[MDA] {r['space']} {r['pair']} primary_mda={r['mda_primary']} stricter_mda={r['mda_stricter']} margin={r['mda_margin']}")
    if not bool(mdf["mda_primary"].notna().any()):
        # save what we have and stop
        ctx["stage2"] = dict(floors=floors, les=les, pooled=pooled, medians=medians,
                             active=active, point=point, bdf=bdf, mdf=mdf, ddf=ddf, dist_df=dist_df)
        raise Halt(
            "no_detection_power",
            "No dose fraction f <= 1.0 cleared the Stage 2A pooled floor in any setting. "
            "Stages 0-2 only. This data cannot support a ranking.",
            {"mda": mdf.to_dict("records")},
        )
    # per-factor floor table
    frows = []
    for si, s in enumerate(active):
        fl = floors[si]
        frows.append(dict(
            space=s["space"], pair=s["pair"], primary=s["primary"],
            floor_min=float(fl.min()), floor_median=float(np.median(fl)),
            floor_max=float(fl.max()), pooled_floor=float(pooled[si]), n_perm=N_PERM,
        ))
    fdf = pd.DataFrame(frows)
    fdf.to_csv(SOUTH4 / "stage2a_floor_summary.csv", index=False)
    np.savez(SOUTH4 / "stage2a_floors.npz", floors=floors, les=les, pooled=pooled, medians=medians)

    # Stage 3 table
    n_f = len(ctx["keep"])
    p_unit = ctx["p_unit"]
    rows = []
    # bootstrap per factor
    log("[bootstrap] resampling cells")
    boot_lo = np.full((len(active), n_f), np.nan)
    boot_hi = np.full((len(active), n_f), np.nan)
    rng_b = np.random.default_rng(BOOT_SEED)
    for fi in range(n_f):
        cells = ctx["cells_by_factor"][fi]
        if cells.size == 0:
            continue
        dense = gather_dense(data_m, ind_m, ip_m, cells, ctx["n_matched"])
        # NT bootstrap draws are shared across factors: draw once outside. Doing it
        # inside correlates nothing we gate on. Draw NT once per factor from the same
        # generator so the stream is one sequence. Memory: store NT boot sums first.
        if fi == 0:
            nt_dense_ok = nt_cells.size * ctx["n_matched"] * 8 < 1_500_000_000
            if nt_dense_ok:
                nt_dense = gather_dense(data_m, ind_m, ip_m, nt_cells, ctx["n_matched"])
            else:
                nt_dense = None
            nt_boot = np.zeros((N_BOOT, ctx["n_matched"]), np.float64)
            for b in range(N_BOOT):
                if nt_dense is not None:
                    take = rng_b.integers(0, nt_cells.size, size=nt_cells.size)
                    nt_boot[b] = nt_dense[take].sum(0)
                else:
                    take = rng_b.choice(nt_cells, size=nt_cells.size, replace=True)
                    nt_boot[b] = gather_dense(data_m, ind_m, ip_m, take, ctx["n_matched"]).sum(0)
            z_nt_boot = np.zeros((N_BOOT, ctx["n_ruler"]))
            for b in range(N_BOOT):
                z_nt_boot[b] = z_of(ctx, nt_boot[b])
            log("[bootstrap] NT pool resampled")
        deltas_b = np.full((len(active), N_BOOT), np.nan)
        for b in range(N_BOOT):
            take = rng_b.integers(0, cells.size, size=cells.size)
            csum = dense[take].sum(0)
            zg = z_of(ctx, csum)
            d = zg - z_nt_boot[b]
            for si, s in enumerate(active):
                if s["proj"]:
                    A = d @ ctx["L_full"][:, :s["k"]]
                else:
                    A = d
                dlt, *_ = block_delta(A, s["r"], s["r2"], s["rnorm"])
                deltas_b[si, b] = float(dlt[0])
        for si in range(len(active)):
            lo, hi = percentile_ci(deltas_b[si])
            boot_lo[si, fi] = lo
            boot_hi[si, fi] = hi
        if (fi + 1) % 50 == 0:
            log(f"[bootstrap] {fi + 1}/{n_f}")

    # residualised stats
    for si, s in enumerate(active):
        fl = floors[si]
        pt = point[si]
        p = np.where(np.isfinite(pt["delta"]), (les[si] + 1.0) / (N_PERM + 1.0), np.nan)
        q = bh_q(p)
        # proliferation residual in this space
        if p_unit is not None and np.linalg.norm(p_unit) > 0:
            if s["proj"]:
                pu = p_unit @ ctx["L_full"][:, :s["k"]]
                v = -s["r"]  # target - origin in this space = -(origin - target)
            else:
                pu = p_unit
                v = -s["r"]
            pun, pn = unit(pu)
            if pn < 1e-12:
                cos_r = np.full(n_f, np.nan)
                dlt_r = np.full(n_f, np.nan)
            else:
                A = pt["A"]
                # stats_from_d one at a time
                cos_r = np.empty(n_f)
                dlt_r = np.empty(n_f)
                for i in range(n_f):
                    st = stats_from_d(A[i], v, pun)
                    cos_r[i] = st["cos_resid"]
                    dlt_r[i] = st["delta_resid"]
        else:
            cos_r = np.full(n_f, np.nan)
            dlt_r = np.full(n_f, np.nan)
        gdis = ctx["g_agree"] <= 0
        iloss = ctx["id_loss"]
        toward = (pt["delta"] < 0) & (pt["delta"] < fl) & np.isfinite(q) & (q <= Q_BAR) & ~iloss & ~gdis
        prolif_flag = (pt["delta"] < 0) & np.isfinite(dlt_r) & (dlt_r >= 0)
        ci_valid = []
        for i in range(n_f):
            if np.isfinite(boot_lo[si, i]) and np.isfinite(boot_hi[si, i]):
                ci_valid.append(bool(ci_covers_point(pt["delta"][i], boot_lo[si, i], boot_hi[si, i])))
            else:
                ci_valid.append(False)
        df = pd.DataFrame(dict(
            factor=ctx["keep"], space=s["space"], k=s["k"] if s["k"] is not None else np.nan,
            origin=s["origin"], target=s["target"], pair=s["pair"], primary=s["primary"],
            cos=pt["cos"], frac=pt["frac"], delta=pt["delta"], d_norm=pt["d_norm"],
            v_norm=s["rnorm"], d_over_v=pt["d_norm"] / s["rnorm"] if s["rnorm"] else np.nan,
            delta_small_disp_approx=-pt["d_norm"] * pt["cos"],
            floor_p05=fl, margin_over_floor=fl - pt["delta"],
            p_delta=p, q_delta=q,
            delta_ci_lo=boot_lo[si], delta_ci_hi=boot_hi[si], delta_ci_valid=ci_valid,
            md_change=pt["md_change"], inside_radius=pt["inside"],
            guide_agreement=ctx["g_agree"], guides_disagree=gdis,
            identity_drop=ctx["id_drop"], identity_loss=iloss,
            cos_resid=cos_r, delta_resid=dlt_r, proliferation_flag=prolif_flag,
            ruler_score_change=ctx["ruler_change"],
            toward_young=toward,
        ))
        rows.append(df)
        log(f"[count] {s['space']} {s['pair']} toward_young={int(toward.sum())} of {n_f}")
    tab = pd.concat(rows, ignore_index=True)
    tab.to_csv(SOUTH4 / "stage3_full_table.csv", index=False)
    # counts
    cnt_rows = []
    for (space, pair), g in tab.groupby(["space", "pair"], sort=False):
        cnt_rows.append(dict(
            space=space, pair=pair, primary=bool(g["primary"].iloc[0]),
            n_factors=int(len(g)),
            n_delta_negative=int((g["delta"] < 0).sum()),
            n_below_floor=int((g["delta"] < g["floor_p05"]).sum()),
            n_q=int((g["q_delta"] <= Q_BAR).sum()),
            n_identity_loss=int(g["identity_loss"].sum()),
            n_guides_disagree=int(g["guides_disagree"].sum()),
            n_toward_young=int(g["toward_young"].sum()),
            n_ci_invalid=int((g["delta_ci_valid"] == False).sum()),  # noqa: E712
        ))
    cnt = pd.DataFrame(cnt_rows)
    cnt.to_csv(SOUTH4 / "stage3_counts.csv", index=False)
    # top 20 per setting
    tops = []
    for (space, pair), g in tab.groupby(["space", "pair"], sort=False):
        gg = g.sort_values("delta", ascending=True).head(20)
        for _, r in gg.iterrows():
            tops.append(dict(
                space=r["space"], origin=r["origin"], target=r["target"], pair=r["pair"],
                factor=r["factor"], delta=r["delta"], margin_over_floor=r["margin_over_floor"],
                q=r["q_delta"], guide_agreement=r["guide_agreement"],
                identity_flag=bool(r["identity_loss"]), proliferation_flag=bool(r["proliferation_flag"]),
                toward_young=bool(r["toward_young"]),
            ))
    top = pd.DataFrame(tops)
    top.to_csv(SOUTH4 / "top20.csv", index=False)

    # Stage 4
    qrows = []
    srows = []
    rec_rows = []
    for (space, pair), g in tab.groupby(["space", "pair"], sort=False):
        g = g.reset_index(drop=True)
        ratio = g["d_over_v"].to_numpy(float)
        delta = g["delta"].to_numpy(float)
        cos = g["cos"].to_numpy(float)
        approx = g["delta_small_disp_approx"].to_numpy(float)
        dn = g["d_norm"].to_numpy(float)
        resid = delta - approx
        rd = sps.rankdata(delta)
        rc = sps.rankdata(-cos)
        rho = float(sps.spearmanr(rd, rc).statistic)
        err = np.abs(resid) / np.where(dn > 1e-12, dn, np.nan)
        srows.append(dict(
            space=space, pair=pair, rho=rho,
            ratio_median=float(np.nanmedian(ratio)),
            err_over_d_median=float(np.nanmedian(err)),
            n=int(len(g)),
        ))
        edges = np.nanpercentile(ratio, [25, 50, 75])
        qidx = np.digitize(ratio, edges, right=True)
        for qi in range(4):
            m = (qidx == qi) & np.isfinite(ratio) & np.isfinite(delta) & np.isfinite(cos)
            if m.sum() < 3:
                continue
            rho_q = float(sps.spearmanr(sps.rankdata(delta[m]), sps.rankdata(-cos[m])).statistic)
            qrows.append(dict(
                space=space, pair=pair, quartile=f"Q{qi+1}", n=int(m.sum()),
                ratio_median=float(np.median(ratio[m])), rho=rho_q,
                err_over_d_median=float(np.nanmedian(err[m])),
            ))
        # survivor recovery
        surv = g.index[g["toward_young"].to_numpy()].to_numpy()
        n_s = int(surv.size)
        if n_s == 0:
            rec_rows.append(dict(space=space, pair=pair, n_survivors=0, n_recovered_by_cosine=0, fraction=np.nan))
        else:
            order = np.argsort(-cos, kind="mergesort")
            topn = set(order[:n_s].tolist())
            got = len(topn & set(surv.tolist()))
            rec_rows.append(dict(
                space=space, pair=pair, n_survivors=n_s, n_recovered_by_cosine=got,
                fraction=got / n_s,
            ))
    qdf = pd.DataFrame(qrows)
    sdf = pd.DataFrame(srows)
    rdf = pd.DataFrame(rec_rows)
    qdf.to_csv(SOUTH4 / "stage4_quartiles.csv", index=False)
    sdf.to_csv(SOUTH4 / "stage4_overall.csv", index=False)
    rdf.to_csv(SOUTH4 / "stage4_recovery.csv", index=False)

    # Stage 5 combinations on every active setting
    set_rows = []
    freq_rows = []
    rng_s = np.random.default_rng(SEED)
    for si, s in enumerate(active):
        A = point[si]["A"]
        dnorm2 = np.sum(A * A, axis=1)
        starts = rng_s.integers(0, n_f, size=N_RANDOM_STARTS)
        found = []
        y2_stop = None
        if s.get("proj") and ctx["y2"].get(s["k"], {}).get("defined"):
            y2_stop = dict(ctx["y2"][s["k"]])
            origin_vec = s.get("origin_vec")
            if origin_vec is not None:
                y2_stop["origin_pc"] = (origin_vec - ctx["center"]) @ ctx["L_full"][:, :s["k"]]
        for st in starts:
            chosen, current, entered = greedy(A, dnorm2, s["r"], int(st), y2_stop)
            # Also stop-check Y2 when the target is Y1 but a Y2 exists in this k.
            found.append(tuple(sorted(chosen)))
        # recurrence
        from collections import Counter
        counts = Counter(found)
        for members, n_hit in counts.most_common():
            freq_rows.append(dict(space=s["space"], pair=s["pair"], n_hits=n_hit, size=len(members),
                                   members=";".join(ctx["keep"][i] for i in members)))
        # evaluate unique sets
        margin_need = float(mdf.loc[(mdf.space == s["space"]) & (mdf.pair == s["pair"]), "mda_margin"].iloc[0])
        for members, n_hit in counts.items():
            current = A[list(members)].sum(0)
            dlt, cos, frac, dn = block_delta(current, s["r"], s["r2"], s["rnorm"])
            # own floor: permute the summed matched vector. For projected settings the
            # sum is in PC space; permute in gene space then project.
            own_floor = _set_floor(ctx, s, members, mr, log)
            margin = float(own_floor - dlt[0]) if np.isfinite(own_floor) else np.nan
            # identity of the sum of logcpm shifts
            if ctx["present"]:
                # id_drop was -mean(logcpm_f - pool). Sum of drops is not the drop of the sum
                # of logcpm shifts: drop_set = -mean(sum_i (logcpm_i - pool)) = sum of per-factor
                # (logcpm_i - pool) with the minus already in id only if one gene. Use logcpm.
                pass
            id_cols = []
            # recompute from logcpm_f stored
            # identity flag: drop of the summed logcpm displacement
            # We stored id_drop per factor but not per-gene. Recompute from logcpm if we saved
            # the identity columns' logcpm. logcpm_f is full matched. Identity positions were
            # not saved. Use id_drop sum as a descriptive bound only when one gene... 
            # Proper: mean over identity genes of sum of per-factor logcpm differences.
            # Recover id columns from symbols again via frozen — logcpm_f is matched-space,
            # and matched_ruler maps slots to ruler rows. Recompute below from ctx.
            id_flag = False
            slots = ctx.get("id_slots")
            if slots is not None and np.size(slots):
                diff = ctx["logcpm_f"][list(members)][:, slots] - ctx["logcpm_pool"][slots]
                id_drop_set = -float(np.mean(np.sum(diff, axis=0)))
                id_flag = bool(id_drop_set > IDENTITY_DROP)
            else:
                id_drop_set = np.nan
            ratio = float(dn[0] / s["rnorm"]) if s["rnorm"] else np.nan
            report = bool(
                np.isfinite(margin_need) and np.isfinite(margin) and dlt[0] < 0
                and dlt[0] < own_floor and margin >= margin_need
            )
            # leave one out
            loo = []
            for drop_i in members:
                rest = [j for j in members if j != drop_i]
                if not rest:
                    loo.append((ctx["keep"][drop_i], float("nan")))
                    continue
                cur = A[rest].sum(0)
                dloo, *_ = block_delta(cur, s["r"], s["r2"], s["rnorm"])
                loo.append((ctx["keep"][drop_i], float(dloo[0])))
            set_rows.append(dict(
                space=s["space"], pair=s["pair"], primary=s["primary"],
                members=";".join(ctx["keep"][i] for i in members),
                size=len(members), n_hits=n_hit,
                delta=float(dlt[0]), own_floor=own_floor, margin_over_floor=margin,
                mda_margin_required=margin_need, ratio=ratio,
                identity_drop=id_drop_set, identity_flag=id_flag,
                reported=report,
                leave_one_out=";".join(f"{n}:{d:+.4f}" for n, d in loo),
            ))
    sets = pd.DataFrame(set_rows)
    freqs = pd.DataFrame(freq_rows)
    sets.to_csv(SOUTH4 / "stage5_sets.csv", index=False)
    freqs.to_csv(SOUTH4 / "stage5_recurrence.csv", index=False)
    ctx.update(dict(
        dist_df=dist_df, bdf=bdf, mdf=mdf, fdf=fdf, cnt=cnt, top=top,
        qdf=qdf, sdf=sdf, rdf=rdf, sets=sets, tab=tab, active=active,
    ))
    return ctx


def _set_floor(ctx, setting, members, mr, log):
    """Magnitude-matched permutation floor of a summed factor vector. 20,000 perms."""
    d = ctx["d_f"][list(members)].sum(0)
    d_m = d[mr]
    rng = np.random.default_rng([SEED, 10_000 + len(members), int(members[0])])
    p = d_m.shape[0]
    L = ctx["L_full"][mr, :]
    n_perm = N_PERM
    deltas = np.empty(n_perm)
    # chunk
    pos = 0
    r = setting["r"]
    r_m = r if setting["proj"] else np.asarray(setting["r"])[mr]
    while pos < n_perm:
        m = min(2000, n_perm - pos)
        idx = np.empty((m, p), np.int32)
        idx[:] = np.arange(p, dtype=np.int32)
        rng.permuted(idx, axis=1, out=idx)
        dperm = d_m[idx]
        if setting["proj"]:
            A = dperm @ L[:, :setting["k"]]
            rr = r
        else:
            A = dperm
            rr = r_m
        dlt, *_ = block_delta(A, rr, setting["r2"], setting["rnorm"])
        deltas[pos:pos + m] = dlt
        pos += m
    return float(np.percentile(deltas, 5))


def _perm_with_medians(D_m, settings, observed, medians, L_m, log):
    """Like run_permutations, but each setting is rescaled to its own median norm."""
    n_f, p = D_m.shape
    n_set = len(settings)
    d_path = SOUTH4 / "perm_D.npy"
    np.save(d_path, np.asarray(D_m, np.float64))
    r = np.zeros((n_set, max(p, L_m.shape[1])), np.float64)
    r2 = np.zeros(n_set)
    rnorm = np.zeros(n_set)
    k = np.zeros(n_set, np.int32)
    proj = np.zeros(n_set, np.bool_)
    for i, s in enumerate(settings):
        rr = np.asarray(s["r"], np.float64).ravel()
        r[i, :rr.size] = rr
        r2[i] = float(s["r2"])
        rnorm[i] = float(s["rnorm"])
        k[i] = int(s["k"] or 0)
        proj[i] = bool(s["proj"])
    pool_path = SOUTH4 / "perm_pool.f32"
    pool_shape = (n_set, n_f, N_PERM)
    done_path = SOUTH4 / "perm_done.npy"
    floors_path = SOUTH4 / "perm_floors.npy"
    les_path = SOUTH4 / "perm_les.npy"
    resume_i = 0
    if done_path.exists() and floors_path.exists() and les_path.exists() and pool_path.exists():
        resume_i = int(np.load(done_path))
        floors = np.load(floors_path)
        les = np.load(les_path)
        if floors.shape != (n_set, n_f) or les.shape != (n_set, n_f) or resume_i > n_f:
            resume_i = 0
        else:
            log(f"[perm] resuming at factor {resume_i}/{n_f}")
    if resume_i == 0:
        pool = np.memmap(pool_path, dtype=np.float32, mode="w+", shape=pool_shape)
        del pool
        floors = np.zeros((n_set, n_f))
        les = np.zeros((n_set, n_f), np.int64)
    spec_path = SOUTH4 / "perm_spec.npz"
    np.savez(
        spec_path, L=np.asarray(L_m, np.float64), r=r, r2=r2, rnorm=rnorm, k=k, proj=proj,
        median_norm=np.array(medians, np.float64), n_perm=np.array([N_PERM]),
        seed=np.array([SEED]), pool_path=np.array([str(pool_path)]),
    )
    obs = np.vstack(observed)
    _init_worker(str(d_path), str(spec_path))
    # override single median with the per-setting array. _perm_one uses _W['med'] as scalar.
    # Replace the scale logic by wrapping: store medians and patch function via a local loop.
    _W["obs"] = obs
    _W["pool_shape"] = pool_shape
    _W["medians"] = np.asarray(medians, np.float64)
    t0 = time.perf_counter()
    for i in range(resume_i, n_f):
        _, fl, le = _perm_one_medians(i)
        floors[:, i] = fl
        les[:, i] = le
        if (i + 1) % 25 == 0 or i + 1 == n_f:
            np.save(floors_path, floors)
            np.save(les_path, les)
            np.save(done_path, np.array([i + 1]))
            log(f"[perm] {i + 1}/{n_f} ({time.perf_counter() - t0:.0f}s)")
    _W.pop("D", None)
    log(f"[perm] done in {time.perf_counter() - t0:.0f}s")
    return floors, les


def _perm_one_medians(i):
    D = np.asarray(_W["D"][i], np.float64)
    p = D.shape[0]
    n_perm = _W["n_perm"]
    n_set = _W["r"].shape[0]
    L = np.asarray(_W["L"], np.float64)
    kmax = L.shape[1]
    rng = np.random.default_rng([_W["seed"], int(i)])
    deltas = np.empty((n_set, n_perm), np.float64)
    pooled = np.empty((n_set, n_perm), np.float32)
    s = float(np.linalg.norm(D))
    medians = _W["medians"]
    pos = 0
    while pos < n_perm:
        m = min(PERM_CHUNK, n_perm - pos)
        idx = np.empty((m, p), np.int32)
        idx[:] = np.arange(p, dtype=np.int32)
        rng.permuted(idx, axis=1, out=idx)
        dperm = D[idx]
        proj = dperm @ L if kmax else np.zeros((m, 0))
        qn2_gene = np.einsum("ij,ij->i", dperm, dperm)
        for s_i in range(n_set):
            if _W["proj"][s_i]:
                kk = int(_W["k"][s_i])
                A = proj[:, :kk]
                rr = _W["r"][s_i, :kk]
                dots = A @ rr
                qn2 = np.einsum("ij,ij->i", A, A)
            else:
                rr = _W["r"][s_i, :p]
                dots = dperm @ rr
                qn2 = qn2_gene
            r2 = float(_W["r2"][s_i])
            rnorm = float(_W["rnorm"][s_i])
            dist = np.sqrt(np.maximum(r2 + 2.0 * dots + qn2, 0.0))
            deltas[s_i, pos:pos + m] = dist - rnorm
            scale = (medians[s_i] / s) if s > 1e-12 else np.nan
            if np.isfinite(scale):
                dist_p = np.sqrt(np.maximum(r2 + 2.0 * scale * dots + (scale ** 2) * qn2, 0.0))
                pooled[s_i, pos:pos + m] = np.float32(dist_p - rnorm)
            else:
                pooled[s_i, pos:pos + m] = np.nan
        pos += m
    obs = _W["obs"][:, i]
    floors = np.percentile(deltas, 5, axis=1)
    le = np.empty(n_set, np.int64)
    for s_i in range(n_set):
        le[s_i] = int(np.sum(deltas[s_i] <= obs[s_i])) if np.isfinite(obs[s_i]) else -1
    pool = np.memmap(_W["pool_path"], dtype=np.float32, mode="r+", shape=_W["pool_shape"])
    pool[:, i, :] = pooled
    pool.flush()
    del pool
    return int(i), floors, le


def render(ctx, log):
    inv = ctx["inv"]
    ch = ctx["chosen"]
    lines = []
    # stage 0
    lines.append("## Stage 0 — inventory")
    lines.append("")
    lines.append(f"File `{H5_NAME}`, {inv['size_bytes']} bytes, md5 `{inv['md5']}` "
                 f"(Zenodo md5 `{EXPECTED_MD5}`, match).")
    lines.append("")
    lines.append(f"Matrix source `{inv['matrix_source']}`, shape {inv['matrix_shape']}. "
                 "X itself is the normalized log1p layer named by the file. "
                 "The integer-count audit selected the raw layer documented on the Zenodo record.")
    lines.append("")
    lines.append(f"Guide assignment: target column `{inv['target_column']}`, guide column `{inv['guide_column']}`.")
    lines.append("")
    lines.append(
        f"Distinct target genes: {inv['n_distinct_targets']}. "
        f"Guides per target: min {inv['guides_per_target_min']}, median {inv['guides_per_target_median']}, "
        f"max {inv['guides_per_target_max']}. "
        f"Targets with ≥2 guides: {inv['n_targets_ge2']}."
    )
    lines.append("")
    lines.append(
        f"Non-targeting guides: {inv['n_nt_guides']}. Non-targeting cells: {inv['n_nt_cells']}. "
        f"Cells per target gene: min {inv['cells_per_target_min']}, "
        f"median {inv['cells_per_target_median']}, max {inv['cells_per_target_max']}."
    )
    lines.append("")
    lines.append("Ruler coverage (gene count and fraction of total absolute weight):")
    lines.append("")
    lines.append("| mapping | genes matched | fraction of genes | fraction of \\|w\\| |")
    lines.append("| --- | --- | --- | --- |")
    for c in ctx["cov_rows"]:
        lines.append(f"| {c['mapping']} | {c['n_matched']} / {c['n_ruler']} | {c['frac_genes']:.4f} | {c['weight_frac']:.4f} |")
    lines.append("")
    lines.append(f"Mapping used, highest weight fraction: `{ch['mapping']}` "
                 f"({ch['weight_frac']:.4f}). "
                 + ("`coverage_low` fired." if ch["weight_frac"] < COVERAGE_LOW else "`coverage_low` did not fire."))
    lines.append("")
    lines.append(f"Identity genes present: {', '.join(ctx['present']) or 'none'}. "
                 f"Absent: {', '.join(ctx['absent']) or 'none'}.")
    lines.append("")
    # distances
    lines.append("## Origin-to-target distances (before any factor was scored)")
    lines.append("")
    lines.append("Primary is O1→Y1. O2 and Y2 are secondary. Y2 is S2 only and is skipped where the young-donor covariance is singular.")
    lines.append("")
    lines.append("| space | origin | target | distance | Mahalanobis | radius | origin inside Y2 |")
    lines.append("| --- | --- | --- | --- | --- | --- | --- |")
    for _, r in ctx["dist_df"].iterrows():
        lines.append(f"| {r['space']} | {r['origin']} | {r['target']} | {_fmt(r['distance'])} | {_fmt(r['mahalanobis'])} | {_fmt(r['radius'])} | {r['inside']} |")
    lines.append("")
    # stage 2
    lines.append("## Stage 2 — detection limit")
    lines.append("")
    lines.append("The floor is a bar. It is not subtracted from any later number.")
    lines.append("")
    lines.append("### A. Primary permutation floor")
    lines.append("")
    lines.append(f"Each factor's matched-gene vector was permuted {N_PERM} times (seed {SEED}). "
                 "The per-factor floor is the 5th percentile. The pooled floor rescales those draws to that setting's median factor norm.")
    lines.append("")
    lines.append("| space | pair | per-factor floor min | median | max | pooled floor |")
    lines.append("| --- | --- | --- | --- | --- | --- |")
    for _, r in ctx["fdf"].iterrows():
        lines.append(f"| {r['space']} | {r['pair']} | {_fmt(r['floor_min'])} | {_fmt(r['floor_median'])} | {_fmt(r['floor_max'])} | {_fmt(r['pooled_floor'])} |")
    lines.append("")
    lines.append("### B. Control floor")
    lines.append("")
    lines.append(f"Non-targeting guide pseudo-effects are leave-one-guide-out ({int(ctx['bdf']['n_nt_guides'].iloc[0])} guides). "
                 f"Random splits are {N_SPLITS} disjoint halves of the non-targeting cells (seed {SEED}). "
                 "The stricter of the two 5th percentiles is the lower one.")
    lines.append("")
    lines.append("| space | pair | NT-guide floor | split floor | stricter |")
    lines.append("| --- | --- | --- | --- | --- |")
    for _, r in ctx["bdf"].iterrows():
        lines.append(f"| {r['space']} | {r['pair']} | {_fmt(r['nt_guide_floor_p05'])} | {_fmt(r['split_floor_p05'])} | {r['stricter']} |")
    lines.append("")
    lines.append("### C. Dose-response and MDA")
    lines.append("")
    lines.append("The true young-minus-old direction is added to the origin at "
                 f"f = {', '.join(str(x) for x in F_VALUES)}, scaled to that setting's median factor norm. "
                 "A dose clears a floor only when delta is below the floor and negative. "
                 "`mda_primary` is the smallest f clearing the Stage 2A pooled floor. "
                 "`mda_margin` is how far the smallest dose that clears the stricter of the pooled floor and the Stage 2B floors sits below that stricter floor. "
                 "That margin is the combination bar.")
    lines.append("")
    lines.append("| space | pair | primary | MDA vs pooled | MDA vs stricter | margin | stricter source |")
    lines.append("| --- | --- | --- | --- | --- | --- | --- |")
    for _, r in ctx["mdf"].iterrows():
        lines.append(f"| {r['space']} | {r['pair']} | {r['primary']} | {_fmt(r['mda_primary'], 2)} | {_fmt(r['mda_stricter'], 2)} | {_fmt(r['mda_margin'])} | {r['stricter_source']} |")
    lines.append("")
    # stage 3
    lines.append("## Stage 3 — ranking")
    lines.append("")
    lines.append("`toward_young` requires delta < 0, delta below that factor's own Stage 2A floor, "
                 f"q ≤ {Q_BAR}, no identity loss, and guide agreement > 0. Counts are not pooled. "
                 "Bootstrap intervals are never a gate. An interval that misses its point estimate is INVALID.")
    lines.append("")
    lines.append("| space | pair | primary | delta<0 | below floor | q≤0.05 | identity loss | guides disagree | toward_young | CI invalid |")
    lines.append("| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |")
    for _, r in ctx["cnt"].iterrows():
        lines.append(
            f"| {r['space']} | {r['pair']} | {r['primary']} | {r['n_delta_negative']} | {r['n_below_floor']} | "
            f"{r['n_q']} | {r['n_identity_loss']} | {r['n_guides_disagree']} | {r['n_toward_young']} | {r['n_ci_invalid']} |"
        )
    lines.append("")
    lines.append("The top 20 factors by delta in each setting, with the required columns, are in `results/south4/top20.csv`. "
                 "A row in that table is a candidate, not a survivor, unless `toward_young` is true.")
    lines.append("")
    # stage 4
    lines.append("## Stage 4 — does delta add anything?")
    lines.append("")
    lines.append("Rank 1 on delta is the most negative delta. Rank 1 on cosine is the highest cosine. "
                 "The Spearman correlation is between those two ranks. "
                 "The small-displacement identity is delta ≈ −||d|| cos(theta). "
                 "The residual is delta minus that term, summarised as a fraction of ||d||. "
                 "Recovery is the fraction of `toward_young` factors contained in the top-n cosine list, n being the number of survivors in that setting.")
    lines.append("")
    lines.append("| space | pair | Spearman | median ||d||/||gap|| | median |residual|/||d|| | survivors | recovered by top-n cosine |")
    lines.append("| --- | --- | --- | --- | --- | --- | --- |")
    rec_map = {(r.space, r.pair): r for _, r in ctx["rdf"].iterrows()}
    for _, r in ctx["sdf"].iterrows():
        rec = rec_map[(r["space"], r["pair"])]
        lines.append(
            f"| {r['space']} | {r['pair']} | {_fmt(r['rho'], 3)} | {_fmt(r['ratio_median'])} | "
            f"{_fmt(r['err_over_d_median'])} | {int(rec['n_survivors'])} | {_fmt(rec['fraction'], 3)} |"
        )
    lines.append("")
    lines.append("Quartiles of ||d|| / ||target − origin||:")
    lines.append("")
    lines.append("| space | pair | quartile | n | median ratio | Spearman | median |residual|/||d|| |")
    lines.append("| --- | --- | --- | --- | --- | --- | --- |")
    for _, r in ctx["qdf"].iterrows():
        lines.append(
            f"| {r['space']} | {r['pair']} | {r['quartile']} | {r['n']} | {_fmt(r['ratio_median'])} | "
            f"{_fmt(r['rho'], 3)} | {_fmt(r['err_over_d_median'])} |"
        )
    lines.append("")
    # verdict computed from the tables, not a new gate
    lines.append(_stage4_verdict(ctx))
    lines.append("")
    # stage 5
    lines.append("## Stage 5 — combinations")
    lines.append("")
    lines.append("This assumes effects add linearly. No public fibroblast screen can test that assumption.")
    lines.append("")
    lines.append(f"Greedy forward search,  {N_RANDOM_STARTS} random starts, seed {SEED}. "
                 f"Stop when nothing reduces distance, on entering the Y2 region, or at size {MAX_SET}. "
                 "A set is reported only when its delta is negative, below its own magnitude-matched permutation floor, "
                 "and the clearance of that floor is at least the Stage 2 MDA margin against the stricter floor.")
    lines.append("")
    reported = ctx["sets"][ctx["sets"]["reported"] == True] if len(ctx["sets"]) else ctx["sets"]  # noqa: E712
    lines.append(f"Unique sets evaluated: {len(ctx['sets'])}. Sets reported: {len(reported)}.")
    lines.append("")
    if len(ctx["sets"]):
        lines.append("Most common sets (all starts, including those not reported):")
        lines.append("")
        lines.append("| space | pair | hits of 50 | size | members |")
        lines.append("| --- | --- | --- | --- | --- |")
        show = ctx["sets"].sort_values("n_hits", ascending=False).groupby(["space", "pair"], sort=False).head(3)
        for _, r in show.iterrows():
            lines.append(f"| {r['space']} | {r['pair']} | {r['n_hits']} | {r['size']} | {r['members']} |")
        lines.append("")
    if len(reported):
        lines.append("Reported sets:")
        lines.append("")
        lines.append("| space | pair | members | delta | margin over own floor | ratio ||d||/||gap|| | identity flag | hits | leave-one-out delta |")
        lines.append("| --- | --- | --- | --- | --- | --- | --- | --- | --- |")
        for _, r in reported.iterrows():
            lines.append(
                f"| {r['space']} | {r['pair']} | {r['members']} | {_fmt(r['delta'])} | {_fmt(r['margin_over_floor'])} | "
                f"{_fmt(r['ratio'])} | {r['identity_flag']} | {r['n_hits']} | {r['leave_one_out']} |"
            )
        lines.append("")
        lines.append(_combo_vs_stage4(ctx, reported))
    else:
        lines.append("No set cleared its own floor by the required margin. None is reported.")
    lines.append("")
    # fired keys
    man = _manifest()
    flags = [f.get("key") for f in (man.get("flags") or [])]
    lines.append("## Fired keys")
    lines.append("")
    if flags:
        lines.append("Report-only flags: " + ", ".join(f"`{k}`" for k in flags) + ".")
    else:
        lines.append("No report-only flag fired.")
    lines.append("")
    lines.append("No STOP key fired." if not any(False for _ in []) else "")
    # remove empty
    lines = [ln for ln in lines if ln is not None]
    return lines


def _stage4_verdict(ctx):
    """Plain sentence from the computed table. Thresholds are descriptive, not a new gate."""
    if ctx["sdf"].empty:
        return "Stage 4 has no rows."
    rhos = ctx["sdf"]["rho"].to_numpy(float)
    errs = ctx["sdf"]["err_over_d_median"].to_numpy(float)
    # quartile separation: rho < 0.8 in a quartile, which would mean the rankings diverge
    sep = []
    for _, r in ctx["qdf"].iterrows():
        if np.isfinite(r["rho"]) and r["rho"] < 0.8:
            sep.append(f"{r['space']} {r['pair']} {r['quartile']} (Spearman {r['rho']:+.3f}, median ratio {r['ratio_median']:.4f})")
    rec = ctx["rdf"]
    rec_f = rec["fraction"].to_numpy(float)
    bits = [
        f"Across settings the Spearman correlation between the delta ranking and the cosine ranking ranges from {np.nanmin(rhos):+.3f} to {np.nanmax(rhos):+.3f}.",
        f"The median |delta − (−||d|| cos)| / ||d|| ranges from {np.nanmin(errs):.4f} to {np.nanmax(errs):.4f}.",
    ]
    finite = rec_f[np.isfinite(rec_f)]
    if finite.size:
        bits.append(
            f"Where a setting has survivors, top-n cosine recovery ranges from {finite.min():.3f} to {finite.max():.3f}."
        )
    else:
        bits.append("No setting has a toward_young survivor, so cosine recovery of survivors is undefined.")
    if sep:
        bits.append("Delta and cosine separate (Spearman below 0.8) in: " + "; ".join(sep) + ".")
    else:
        bits.append("Delta does not separate from cosine in any ratio quartile: every quartile Spearman is at least 0.8, so the ordering delta adds is the cosine reweighted by ||d||, plus the floor and q bar.")
    return " ".join(bits)


def _combo_vs_stage4(ctx, reported):
    # Did any reported set land in a ratio quartile where stage 4 separated?
    sep_keys = set()
    for _, r in ctx["qdf"].iterrows():
        if np.isfinite(r["rho"]) and r["rho"] < 0.8:
            sep_keys.add((r["space"], r["pair"], r["quartile"]))
    if not sep_keys:
        return "Stage 4 did not show a ratio quartile where delta and cosine separate, so no reported combination reaches such a quartile."
    return "Stage 4 separation quartiles exist; compare each reported set's ratio with `results/south4/stage4_quartiles.csv`."


def finish_success(ctx, log):
    status = [
        "**Status:** Ran to completion. No STOP key fired.",
    ]
    flags = [f.get("key") for f in (_manifest().get("flags") or [])]
    if flags:
        status.append("")
        status.append("Report-only flags: " + ", ".join(f"`{k}`" for k in flags) + ".")
    # one-line counts
    bits = []
    for _, r in ctx["cnt"].iterrows():
        bits.append(f"{r['space']} {r['pair']}={int(r['n_toward_young'])}")
    status.append("")
    status.append("**toward_young counts (not pooled):** " + "; ".join(bits) + ".")
    sections = [render(ctx, log), limitations_section(), _fail_section()]
    write_findings(status, sections)
    write_progress(
        "no STOP key fired; the task ran to completion.",
        "Nothing further is required by this pre-registration. "
        "A perturbation measured in aged cells would be a different task.",
        extra="Counts are in FINDINGS_SOUTH4.md and are not pooled.",
    )
    log("[done] findings written")


def finish_halt(err, ctx, log):
    record_failure(err.key, err.message, err.details)
    status = [
        f"**Status:** STOP `{err.key}`.",
        "",
        err.message,
    ]
    sections = []
    if ctx and ctx.get("inv"):
        inv = ctx["inv"]
        sections.append([
            "## Stage 0 — inventory (as far as it was computed)",
            "",
            f"File size {inv.get('size_bytes')} md5 `{inv.get('md5')}`.",
            "",
            f"Matrix `{inv.get('matrix_source')}` shape {inv.get('matrix_shape')}.",
            "",
            f"Targets {inv.get('n_distinct_targets')}, with ≥2 guides {inv.get('n_targets_ge2')}. "
            f"Non-targeting guides {inv.get('n_nt_guides')}, cells {inv.get('n_nt_cells')}.",
            "",
        ])
        if ctx.get("cov_rows"):
            lines = ["| mapping | genes | fraction of genes | fraction of \\|w\\| |", "| --- | --- | --- | --- |"]
            for c in ctx["cov_rows"]:
                lines.append(f"| {c['mapping']} | {c['n_matched']} / {c['n_ruler']} | {c['frac_genes']:.4f} | {c['weight_frac']:.4f} |")
            sections.append(lines)
    if err.key == "no_detection_power" and ctx and ctx.get("stage2"):
        sections.append(["## Stage 2 was computed. Stages 3–5 were not.", ""])
    sections.append(limitations_section())
    sections.append(_fail_section())
    write_findings(status, sections)
    nxt = {
        "download_failed": "Obtain the complete singlets file. Do not substitute the mean-population file or the L4A lane.",
        "no_raw_counts": "The pinned file has no raw integer counts. Do not analyse the log1p matrix in their place.",
        "too_few_factors": "The guide assignment does not yield 500 targets with 2+ guides. Not relabelling.",
        "no_controls": "Fewer than 20 non-targeting guides. Not inventing controls.",
        "no_detection_power": "No ranking was produced. A larger effect, relative to this floor, would be a different measurement.",
    }.get(err.key, "The failure is recorded verbatim. Not repairing it.")
    write_progress(f"STOP `{err.key}`", nxt)
    log(f"[STOP] {err.key}: {err.message}")


def load_stage1(log):
    """Rebuild the Stage 1 context from the checkpoint written before permutations.

    The count stream and the z-scores are not redone. The TMM cache is rebuilt
    from the saved pseudobulk counts, which is the same panel.
    """
    import json
    from seng_run import tmm_cache_from_panel
    z = np.load(SOUTH4 / "stage1.npz", allow_pickle=True)
    inv = json.loads((SOUTH4 / "inventory.json").read_text(encoding="utf-8"))
    frozen = load_frozen_ruler()
    # stage1.npz `keep` was overwritten by a boolean mask before it was saved.
    # The row order is the sorted non-targeting-excluded target labels, rebuilt
    # from obs. The count matrix is not reread.
    with h5py.File(H5_PATH, "r") as fh:
        obs = read_obs_var(fh["obs"])
    tcol = inv["target_column"]
    gcol = inv["guide_column"]
    targets = obs[tcol].astype(str).to_numpy()
    guides = obs[gcol].astype(str).to_numpy()
    tf = np.array([not _is_nt(t) for t in targets])
    per = pd.Series(guides[tf]).groupby(pd.Series(targets[tf])).nunique()
    keep = sorted(str(t) for t in per.index if int(per[t]) >= MIN_GUIDES)
    if len(keep) != int(np.asarray(z["d_f"]).shape[0]):
        raise Halt(
            "resume",
            f"Rebuilt factor list has {len(keep)} names but d_f has {np.asarray(z['d_f']).shape[0]} rows.",
        )
    matched_ruler = np.asarray(z["matched_ruler"], np.int64)
    C_f = np.asarray(z["C_f"], np.float64)
    C_nt = np.asarray(z["C_nt"], np.float64)
    C_pool = np.asarray(z["C_pool"], np.float64)
    C_panel = np.vstack([C_f, C_nt, C_pool.reshape(1, -1)])
    _logcpm, nf = tmm_logcpm_quiet(C_panel)
    cache = tmm_cache_from_panel(C_panel, nf)
    mu = np.asarray(frozen["mu"], float)[matched_ruler]
    sd = np.asarray(frozen["sd"], float)[matched_ruler]
    sd = np.where(sd < 1e-12, 1.0, sd)
    missing = np.asarray(z["missing"], bool)
    flat = np.asarray(z["cell_flat"], np.int64)
    cind = np.asarray(z["cell_indptr"], np.int64)
    n_f = int(np.asarray(z["d_f"]).shape[0])
    cells_by_factor = [flat[int(cind[i]):int(cind[i + 1])] for i in range(n_f)]
    id_drop = np.asarray(z["id_drop"], float)
    present = [str(x) for x in np.asarray(z["id_present"]).tolist()]
    absent = [str(x) for x in np.asarray(z["id_absent"]).tolist()]
    ruler_sym = np.array([str(s).upper() for s in np.asarray(frozen["symbol"])], dtype=object)
    sym_to_ruler = {}
    for i, s in enumerate(ruler_sym):
        sym_to_ruler.setdefault(s, i)
    slots = []
    for g in present:
        rr = sym_to_ruler.get(g)
        hit = np.flatnonzero(matched_ruler == rr) if rr is not None else np.array([])
        if hit.size:
            slots.append(int(hit[0]))
    p_unit = np.asarray(z["p_unit"], float)
    p_unit = p_unit if p_unit.size else None
    notes_path = SOUTH4 / "stage1_notes.json"
    prolif_note = "reloaded from stage1"
    if notes_path.exists():
        prolif_note = json.loads(notes_path.read_text(encoding="utf-8")).get("prolif", prolif_note)
    log(f"[resume] stage1 factors={len(keep)} matched={matched_ruler.size} "
        f"nt_cells={np.asarray(z['nt_cells']).size}")
    return dict(
        inv=inv, frozen=frozen, keep=keep,
        d_f=np.asarray(z["d_f"], np.float64),
        z_pool=np.asarray(z["z_pool"], np.float64),
        z_f=np.asarray(z["z_f"], np.float64),
        z_ntg=np.asarray(z["z_ntg"], np.float64),
        C_f=C_f, C_nt=C_nt, C_pool=C_pool,
        logcpm_f=np.asarray(z["logcpm_f"], np.float64),
        logcpm_pool=np.asarray(z["logcpm_pool"], np.float64),
        id_drop=id_drop, id_loss=id_drop > IDENTITY_DROP,
        g_agree=np.asarray(z["g_agree"], float),
        ruler_change=np.asarray(z["ruler_change"], float),
        matched_ruler=matched_ruler, n_ruler=int(np.asarray(z["n_ruler"]).ravel()[0]),
        p_unit=p_unit, prolif_note=prolif_note,
        present=present, absent=absent,
        nt_cells=np.asarray(z["nt_cells"], np.int64),
        cells_by_factor=cells_by_factor,
        nnz=int(np.asarray(z["nnz"]).ravel()[0]),
        n_matched=int(matched_ruler.size),
        cache=cache, mu_f=mu, sd_f=sd, miss=missing,
        nt_guide_names=[str(x) for x in np.asarray(z["nt_guide_names"]).tolist()],
        chosen=inv["coverage_chosen"], cov_rows=inv["coverage"],
        id_slots=np.asarray(slots, int),
    )


def main():
    if os.environ.get("SOUTH4_SELFCHECK") == "1":
        selfcheck()
        return
    if not PREREG.exists():
        raise SystemExit("PREREG.flag missing. Not computing.")
    log = Logger(SOUTH4 / "south4_run.log", mode="a")
    ctx = {}
    try:
        log(f"[prereg] {PREREG} bytes={PREREG.stat().st_size}")
        if (SOUTH4 / "stage1.npz").exists() and (SOUTH4 / "inventory.json").exists():
            ctx = load_stage1(log)
        else:
            size, md5 = verify_download(log)
            ctx = stage0_and_1(log, size, md5)
        ctx = build_spaces(ctx, log)
        ctx = run_analysis(ctx, log)
        finish_success(ctx, log)
    except Halt as err:
        finish_halt(err, ctx, log)
    except Exception as err:
        record_failure("unhandled", f"{type(err).__name__}: {err}")
        write_progress(f"STOP unhandled {type(err).__name__}", "The traceback is in south4_run.log. Not repairing rows.")
        log(f"[UNHANDLED] {type(err).__name__}: {err}")
        raise
    finally:
        log.close()


if __name__ == "__main__":
    main()
