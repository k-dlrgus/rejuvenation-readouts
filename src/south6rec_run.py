"""SOUTH6REC — recovery control for the SOUTH6 readout, on a depth grid.

The pre-registration is results/south6rec/PREREG.flag, written before any
number. This script reads it and does not rewrite it. It reads results/south6/
as stored and writes nothing there. It imports src/south6_run.py,
src/toward_run.py, src/seng_run.py and src/lowdim_common.py unchanged and does
not refit the frozen ruler or the TOWARD anchors.

Run from the repository root:   python src/south6rec_run.py

Stage 0 rebuilds the SOUTH6 panel, settings, split floor and dose response from
the matched-gene CSR cache and checks them against the stored files. Any
mismatch is STOP `south6_not_reproduced` and no Stage 1 number is computed.
Stage 1 plants f * g_m into null NT pseudobulks of n cells and reads delta in
S1 and S2_k20/50/100 (O1 -> Y1).
"""
from __future__ import annotations

import sys
import time
from pathlib import Path

import h5py
import numpy as np
import pandas as pd
from scipy import sparse, stats as sps

sys.path.insert(0, str(Path(__file__).resolve().parent))

from brain_phase1_common import Logger  # noqa: E402
from fibro2_common import load_frozen_ruler  # noqa: E402
from target_common import unit  # noqa: E402
from toward_run import ANCHORS_PATH, load_gtex_z, tmm_logcpm_quiet  # noqa: E402
from seng_run import file_md5, read_obs_var, tmm_cache_from_panel, z_from_counts  # noqa: E402
from lowdim_common import fit_pca  # noqa: E402
from south6_run import (  # noqa: E402
    EXPECTED_MD5, EXPECTED_SIZE, F_VALUES, GUIDE_COL, H5_PATH, KMAX, N_SPLITS,
    NON_LABEL, S2_KS, SEED, SOUTH6, TARGET_COL, block_delta,
)

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "results" / "south6rec"
PREREG = OUT / "PREREG.flag"
FINDINGS = ROOT / "FINDINGS_SOUTH6REC.md"

BOOT_SEED = 20260918
DEPTHS = (106, 390, 2500)
DOSES = (0.0, 0.25, 0.50, 1.00)
R = 200
N_BOOT = 2000
POWER_BAR = 0.80
TOL = 1e-4
MIN_GUIDES_R1 = 6
SPACES = ("S1",) + tuple(f"S2_k{k}" for k in S2_KS)
VARIANTS = ("primary", "R1_guide", "R2_pool", "R3_fullgap")


class Stop(Exception):
    def __init__(self, key, message):
        super().__init__(message)
        self.key = key
        self.message = message


# ---------------------------------------------------------------- helpers

def check_num(rows, block, item, new, stored, tol=TOL, gating=True):
    new = np.asarray(new, np.float64).ravel()
    stored = np.asarray(stored, np.float64).ravel()
    if new.shape != stored.shape:
        rows.append(dict(block=block, item=item, kind="numeric", n=int(stored.size),
                         max_abs_diff=np.nan, max_scaled_diff=np.nan, tol=tol,
                         stored_example=np.nan, new_example=np.nan,
                         n_mismatch=int(max(new.size, stored.size)), passed=False, gating=gating,
                         note=f"shape new {new.shape} stored {stored.shape}"))
        return
    # Chunked, so large matrices are compared without full-size temporaries.
    max_abs = 0.0
    max_sc = 0.0
    n_bad = 0
    step = 2_000_000
    for a in range(0, stored.size, step):
        x = new[a:a + step]
        y = stored[a:a + step]
        both_nan = np.isnan(x) & np.isnan(y)
        diff = np.where(both_nan, 0.0, np.abs(x - y))
        sc = diff / np.maximum(1.0, np.abs(np.nan_to_num(y)))
        bad = ~(sc <= tol) & ~both_nan
        n_bad += int(bad.sum())
        if diff.size:
            max_abs = max(max_abs, float(np.nanmax(diff)) if np.isfinite(diff).any() else np.inf)
            max_sc = max(max_sc, float(np.nanmax(sc)) if np.isfinite(sc).any() else np.inf)
    rows.append(dict(block=block, item=item, kind="numeric", n=int(stored.size),
                     max_abs_diff=max_abs, max_scaled_diff=max_sc, tol=tol,
                     stored_example=float(stored[0]) if stored.size else np.nan,
                     new_example=float(new[0]) if new.size else np.nan,
                     n_mismatch=n_bad, passed=bool(n_bad == 0), gating=gating, note=""))


def check_cat(rows, block, item, new, stored, gating=True):
    new = np.asarray(new).ravel()
    stored = np.asarray(stored).ravel()
    same = new.shape == stored.shape and bool(np.all(new == stored))
    n_mis = int(np.sum(new != stored)) if new.shape == stored.shape else int(max(new.size, stored.size))
    rows.append(dict(block=block, item=item, kind="categorical", n=int(stored.size),
                     max_abs_diff=np.nan, max_scaled_diff=np.nan, tol=0.0,
                     stored_example=str(stored[0]) if stored.size else "",
                     new_example=str(new[0]) if new.size else "",
                     n_mismatch=n_mis, passed=same, gating=gating, note=""))


def indicator(rows_list, n_cols):
    """Sparse 0/1 matrix, one row per list of column positions."""
    indptr = np.zeros(len(rows_list) + 1, np.int64)
    for i, r in enumerate(rows_list):
        indptr[i + 1] = indptr[i] + len(r)
    idx = np.concatenate([np.asarray(r, np.int64) for r in rows_list]) if rows_list else np.zeros(0, np.int64)
    return sparse.csr_matrix((np.ones(idx.size), idx, indptr), shape=(len(rows_list), n_cols))


def wilson(k, n, z=1.959963984540054):
    if n == 0:
        return np.nan, np.nan
    p = k / n
    den = 1 + z * z / n
    c = (p + z * z / (2 * n)) / den
    h = z * np.sqrt(p * (1 - p) / n + z * z / (4 * n * n)) / den
    return c - h, c + h


def fmt(x, d=3):
    if x is None or (isinstance(x, float) and not np.isfinite(x)):
        return "NA"
    return f"{x:.{d}f}"


# ---------------------------------------------------------------- stage 0

def stage0(ctx, log):
    rows = []
    ctx["rep_rows"] = rows
    t0 = time.time()

    size = H5_PATH.stat().st_size
    md5 = file_md5(H5_PATH)
    log(f"[file] size={size} md5={md5} ({time.time() - t0:.0f}s)")
    if size != EXPECTED_SIZE or md5 != EXPECTED_MD5:
        ctx["rep"] = pd.DataFrame(rows)
        raise Stop("file_mismatch", f"size {size} md5 {md5}; expected {EXPECTED_SIZE} {EXPECTED_MD5}.")
    ctx.update(md5=md5, size=size)

    with h5py.File(H5_PATH, "r") as f:
        obs = read_obs_var(f["obs"])
    targets = obs[TARGET_COL].astype(str).to_numpy()
    guides = obs[GUIDE_COL].astype(str).to_numpy()
    nt_obs = np.flatnonzero(targets == NON_LABEL).astype(np.int64)
    del obs

    s1 = np.load(SOUTH6 / "stage1.npz", allow_pickle=True)
    nt_cells = np.asarray(s1["nt_cells"], np.int64)
    if not np.array_equal(nt_obs, nt_cells):
        ctx["rep"] = pd.DataFrame(rows)
        raise Stop("nt_cells_mismatch", f"obs NT cells {nt_obs.size} vs stored {nt_cells.size}; "
                   f"sets equal={set(nt_obs.tolist()) == set(nt_cells.tolist())}.")
    check_cat(rows, "obs", "NT cells (guide_target == 'non') equal stage1.npz nt_cells", nt_obs, nt_cells)
    nt_guide_names = [str(x) for x in s1["nt_guide_names"]]
    nt_guide_of = guides[nt_cells]
    gpos = {g: i for i, g in enumerate(nt_guide_names)}
    if set(nt_guide_of) != set(nt_guide_names):
        ctx["rep"] = pd.DataFrame(rows)
        raise Stop("nt_cells_mismatch", "NT guide labels in obs differ from stage1.npz nt_guide_names.")
    guide_idx = np.array([gpos[g] for g in nt_guide_of], np.int64)
    matched_ruler = np.asarray(s1["matched_ruler"], np.int64)
    n_matched = matched_ruler.size
    cell_flat = np.asarray(s1["cell_flat"], np.int64)
    cell_indptr = np.asarray(s1["cell_indptr"], np.int64)
    keep = [str(x) for x in s1["keep"]]
    n_f = len(keep)
    nnz = int(np.asarray(s1["nnz"]).ravel()[0])
    n_ruler = int(np.asarray(s1["n_ruler"]).ravel()[0])

    data = np.memmap(SOUTH6 / "matched_data.f64", dtype=np.float64, mode="r", shape=(nnz,))
    indices = np.memmap(SOUTH6 / "matched_indices.i32", dtype=np.int32, mode="r", shape=(nnz,))
    indptr = np.load(SOUTH6 / "matched_indptr.npy")

    def rows_sum(cells):
        parts_d, parts_i = [], []
        for c in cells:
            a, b = int(indptr[c]), int(indptr[c + 1])
            if b > a:
                parts_d.append(np.asarray(data[a:b]))
                parts_i.append(np.asarray(indices[a:b]))
        if not parts_d:
            return np.zeros(n_matched)
        return np.bincount(np.concatenate(parts_i), weights=np.concatenate(parts_d), minlength=n_matched)

    # In-memory CSR of the NT cells, rows in nt_cells order.
    t0 = time.time()
    lens = np.array([int(indptr[c + 1] - indptr[c]) for c in nt_cells], np.int64)
    ip = np.zeros(nt_cells.size + 1, np.int64)
    ip[1:] = np.cumsum(lens)
    dat = np.empty(int(ip[-1]), np.float64)
    ind = np.empty(int(ip[-1]), np.int32)
    for j, c in enumerate(nt_cells):
        a, b = int(indptr[c]), int(indptr[c + 1])
        dat[ip[j]:ip[j + 1]] = data[a:b]
        ind[ip[j]:ip[j + 1]] = indices[a:b]
    X_nt = sparse.csr_matrix((dat, ind, ip), shape=(nt_cells.size, n_matched))
    log(f"[nt csr] cells={nt_cells.size} nnz={X_nt.nnz} ({time.time() - t0:.0f}s)")

    # (a) counts
    t0 = time.time()
    C_f = np.zeros((n_f, n_matched), np.float64)
    for fi in range(n_f):
        C_f[fi] = rows_sum(cell_flat[cell_indptr[fi]:cell_indptr[fi + 1]])
    log(f"[C_f] rebuilt from cells ({time.time() - t0:.0f}s)")
    del data, indices
    C_nt = np.zeros((len(nt_guide_names), n_matched))
    for gi in range(len(nt_guide_names)):
        C_nt[gi] = np.asarray(X_nt[np.flatnonzero(guide_idx == gi)].sum(0)).ravel()
    C_pool = np.asarray(X_nt.sum(0)).ravel()
    st = s1["C_f"]
    check_num(rows, "(a) counts", "C_f, 1836 factors x 20264 genes", C_f, st, tol=0.0)
    del st
    check_num(rows, "(a) counts", "C_nt, 78 guides x 20264 genes", C_nt, s1["C_nt"], tol=0.0)
    check_num(rows, "(a) counts", "C_pool", C_pool, s1["C_pool"], tol=0.0)
    check_num(rows, "(a) counts", "C_pool equals sum of C_nt", C_pool, C_nt.sum(0), tol=0.0)

    # (b) panel transform
    t0 = time.time()
    C_panel = np.vstack([C_f, C_nt, C_pool.reshape(1, -1)])
    del C_f
    logcpm, nf = tmm_logcpm_quiet(C_panel)
    frozen = load_frozen_ruler()
    mu = np.asarray(frozen["mu"], float)[matched_ruler]
    sd = np.asarray(frozen["sd"], float)[matched_ruler]
    sd_safe = np.where(sd < 1e-12, 1.0, sd)
    missing = (C_panel == 0).all(0)
    logcpm -= mu
    logcpm /= sd_safe
    logcpm[:, missing] = 0.0
    Z_m = logcpm
    z_pool_m = Z_m[-1].copy()
    d_m = Z_m[:n_f] - z_pool_m
    cache = tmm_cache_from_panel(C_panel, nf)
    del C_panel, Z_m, logcpm
    log(f"[panel] transform rebuilt ({time.time() - t0:.0f}s)")
    check_cat(rows, "(b) transform", "missing-gene mask", missing, s1["missing"])
    zp = np.asarray(s1["z_pool"], float)
    check_num(rows, "(b) transform", "z_pool (matched genes)", z_pool_m, zp[matched_ruler])
    unm = np.ones(n_ruler, bool)
    unm[matched_ruler] = False
    check_num(rows, "(b) transform", "z_pool unmatched genes are 0", zp[unm], np.zeros(int(unm.sum())), tol=0.0)
    st = s1["d_f"]
    check_num(rows, "(b) transform", "d_f, matched genes", d_m, st[:, matched_ruler])
    check_num(rows, "(b) transform", "d_f, unmatched genes are 0", st[:, unm].ravel(),
              np.zeros(int(unm.sum()) * n_f), tol=0.0)
    del st, zp
    s1.close()

    # (c) spaces
    gtex = load_gtex_z(log)
    Zg = np.asarray(gtex["Z"], np.float64)
    age = np.asarray(gtex["age"]).astype(str)
    center = Zg.mean(0)
    sp = fit_pca(Zg - center, KMAX)
    L_full = np.asarray(sp.L, np.float64)
    log(f"[S2] PCA on GTEx donors n={Zg.shape[0]} p={Zg.shape[1]} k={L_full.shape[1]}")
    anc = np.load(ANCHORS_PATH, allow_pickle=True)
    c_young = np.asarray(anc["c_young"], float).ravel()
    c_old = np.asarray(anc["c_old"], float).ravel()
    g = c_young - c_old
    otd = pd.read_csv(SOUTH6 / "origin_target_distances.csv")
    sp_r = {}
    for s in SPACES:
        if s == "S1":
            r = c_old - c_young
        else:
            k = int(s.split("k")[1])
            r = (c_old - center) @ L_full[:, :k] - (c_young - center) @ L_full[:, :k]
        sp_r[s] = dict(r=r, r2=float(r @ r), rnorm=float(np.linalg.norm(r)))
        st = otd[(otd.space == s) & (otd.origin == "O1") & (otd.target == "Y1")]["distance"].to_numpy()
        check_num(rows, "(c) spaces", f"||O1 - Y1|| {s}", sp_r[s]["rnorm"], st)
    Lm = L_full[matched_ruler]

    def proj(E_m, s):
        """E_m rows in matched coordinates (unmatched = 0) -> coordinates in space s."""
        if s == "S1":
            full = np.zeros((E_m.shape[0], n_ruler))
            full[:, matched_ruler] = E_m
            return full
        return E_m @ Lm[:, :int(s.split("k")[1])]

    # (d) per-factor stats
    t0 = time.time()
    tab = pd.read_csv(SOUTH6 / "stage3_full_table.csv",
                      usecols=["factor", "space", "pair", "delta", "cos", "frac", "d_norm"])
    tab = tab[tab.pair == "O1_Y1"]
    medians = {}
    for s in SPACES:
        if s == "S1":
            # d_f is 0 off the matched genes, so d . r only needs r on them; r2 and ||r|| stay full.
            D, r_use = d_m, sp_r[s]["r"][matched_ruler]
        else:
            D, r_use = proj(d_m, s), sp_r[s]["r"]
        dl, cs, fr, dn = block_delta(D, r_use, sp_r[s]["r2"], sp_r[s]["rnorm"])
        medians[s] = float(np.median(dn))
        del D
        t = tab[tab.space == s].set_index("factor").reindex(keep)
        for name, new in (("delta", dl), ("cos", cs), ("frac", fr), ("d_norm", dn)):
            check_num(rows, f"(d) stage3 {s} O1_Y1", name, new, t[name].to_numpy())
    del tab
    log(f"[stage3] per-factor stats compared ({time.time() - t0:.0f}s)")

    def zm(counts):
        z, _, _, _ = z_from_counts(np.asarray(counts, np.float64).ravel(), cache, mu, sd, missing)
        return z

    # (e) split floor
    t0 = time.time()
    rng = np.random.default_rng(SEED)
    half = nt_cells.size // 2
    A, B = [], []
    for _ in range(N_SPLITS):
        perm = rng.permutation(nt_cells.size)
        A.append(perm[:half])
        B.append(perm[half:])
    CA = np.asarray((indicator(A, nt_cells.size) @ X_nt).todense())
    CB = np.asarray((indicator(B, nt_cells.size) @ X_nt).todense())
    E = np.vstack([zm(CA[i]) - zm(CB[i]) for i in range(N_SPLITS)])
    del CA, CB
    b2 = pd.read_csv(SOUTH6 / "stage2b_control_floor.csv")
    for s in SPACES:
        dl, *_ = block_delta(proj(E, s), sp_r[s]["r"], sp_r[s]["r2"], sp_r[s]["rnorm"])
        row = b2[(b2.space == s) & (b2.pair == "O1_Y1") & (b2.control_set == "all78")]
        check_num(rows, f"(e) stage2b {s} O1_Y1 all78", "split_floor_p05", np.percentile(dl, 5),
                  row["split_floor_p05"].to_numpy())
        check_num(rows, f"(e) stage2b {s} O1_Y1 all78", "split_median", np.median(dl),
                  row["split_median"].to_numpy())
    del E
    log(f"[stage2b] split floor recomputed ({time.time() - t0:.0f}s)")

    # (f) dose response
    dr = pd.read_csv(SOUTH6 / "stage2c_dose_response.csv")
    for s in SPACES:
        if s == "S1":
            un, _ = unit(g)
        else:
            un, _ = unit(g @ L_full[:, :int(s.split("k")[1])])
        new = [float(block_delta((f * medians[s]) * un, sp_r[s]["r"], sp_r[s]["r2"], sp_r[s]["rnorm"])[0][0])
               for f in F_VALUES]
        st = dr[(dr.space == s) & (dr.pair == "O1_Y1")].sort_values("f")["delta"].to_numpy()
        check_num(rows, f"(f) stage2c {s} O1_Y1", "delta at f = 0..1", new, st)

    rep = pd.DataFrame(rows)
    rep.to_csv(OUT / "stage0_reproduction.csv", index=False)
    failed = rep[(~rep.passed) & rep.gating]
    log(f"[stage0] checks={len(rep)} failed={len(failed)}")

    ctx.update(
        rep=rep, X_nt=X_nt, nt_cells=nt_cells, guide_idx=guide_idx, nt_guide_names=nt_guide_names,
        C_nt=C_nt, C_pool=C_pool, z_pool_m=z_pool_m, zm=zm, proj=proj, sp_r=sp_r, g=g,
        matched_ruler=matched_ruler, missing=missing, n_ruler=n_ruler, L_full=L_full, Lm=Lm,
        Zg=Zg, age=age, n_anchor_young=int(anc["n_young"]), n_anchor_old=int(anc["n_old"]),
    )
    ctx["rep"] = rep
    if len(failed):
        ctx["fail_items"] = failed
        raise Stop("south6_not_reproduced",
                   "; ".join(f"{r.block} / {r.item}: {r.n_mismatch} of {r.n} outside tol "
                             f"(max scaled {r.max_scaled_diff})" for r in failed.itertuples()))
    return ctx


# ---------------------------------------------------------------- stage 1

def primary_draws(n, N):
    """Disjoint-block draws: positions into the NT pool, plus their permutation id."""
    rng = np.random.default_rng([SEED, n])
    m = N // n
    per = m // 4
    out = {di: [] for di in range(len(DOSES))}
    perm_of = {di: [] for di in range(len(DOSES))}
    n_perm = 0
    while min(len(v) for v in out.values()) < R:
        perm = rng.permutation(N)
        for j in range(4 * per):
            di = j % 4
            if len(out[di]) < R:
                out[di].append(np.sort(perm[j * n:(j + 1) * n]))
                perm_of[di].append(n_perm)
        n_perm += 1
    return out, perm_of, dict(m=m, blocks_used=4 * per, per_dose_per_perm=per, n_perm=n_perm)


def guide_draws(n, guide_idx, n_guides):
    rng = np.random.default_rng([SEED, n, 1])
    cells_of = [np.flatnonzero(guide_idx == gi) for gi in range(n_guides)]
    out = {di: [] for di in range(len(DOSES))}
    chosen = {di: [] for di in range(len(DOSES))}
    for di in range(len(DOSES)):
        for _ in range(R):
            order = rng.permutation(n_guides)
            pick, tot = [], 0
            for gi in order:
                pick.append(int(gi))
                tot += cells_of[gi].size
                if len(pick) >= MIN_GUIDES_R1 and tot >= n:
                    break
            if tot < n:
                raise Stop("r1_undefined", f"guide-structured draw at n={n} cannot reach n cells.")
            pool = np.concatenate([cells_of[gi] for gi in pick])
            out[di].append(np.sort(rng.choice(pool, size=n, replace=False)))
            chosen[di].append(pick)
    return out, chosen


def overlap_stats(draws, n, N):
    S = indicator(draws, N)
    O = (S @ S.T).toarray()
    iu = np.triu_indices(O.shape[0], 1)
    v = O[iu] / n
    return float(v.mean()), float(v.max()), float(np.mean(v == 0))


def stage1(ctx, log):
    X_nt = ctx["X_nt"]
    N = X_nt.shape[0]
    zm = ctx["zm"]
    proj = ctx["proj"]
    sp_r = ctx["sp_r"]
    g = ctx["g"]
    mr = ctx["matched_ruler"]
    meas = np.zeros(ctx["n_ruler"], bool)
    meas[mr[~ctx["missing"]]] = True
    g_m = np.where(meas, g, 0.0)
    C_pool = ctx["C_pool"]
    C_nt = ctx["C_nt"]
    z_pool_m = ctx["z_pool_m"]
    n_guides = len(ctx["nt_guide_names"])

    # Plant vectors per space.
    plant = {}
    for s in SPACES:
        if s == "S1":
            plant[s] = dict(meas=g_m, full=g)
        else:
            L = ctx["L_full"][:, :int(s.split("k")[1])]
            plant[s] = dict(meas=g_m @ L, full=g @ L)
    gm_ratio = {s: float(np.linalg.norm(plant[s]["meas"]) / np.linalg.norm(plant[s]["full"])) for s in SPACES}
    for s in SPACES:
        log(f"[plant] {s} ||g_m||={np.linalg.norm(plant[s]['meas']):.4f} ||g||={np.linalg.norm(plant[s]['full']):.4f} "
            f"ratio={gm_ratio[s]:.4f}")

    long_rows = []
    ov_rows = []
    design_rows = []
    for n in DEPTHS:
        t0 = time.time()
        pd_draws, perm_of, info = primary_draws(n, N)
        gd_draws, gchosen = guide_draws(n, ctx["guide_idx"], n_guides)
        design_rows.append(dict(depth=n, pool=N, m_disjoint_blocks=info["m"],
                                blocks_used_per_perm=info["blocks_used"],
                                draws_per_dose_per_perm=info["per_dose_per_perm"],
                                permutations_used=info["n_perm"],
                                cells_used_per_perm=info["blocks_used"] * n))
        for di, f in enumerate(DOSES):
            pids = np.asarray(perm_of[di])
            max_disjoint = int(np.max(np.bincount(pids)))
            mo, mx, fz = overlap_stats(pd_draws[di], n, N)
            ov_rows.append(dict(depth=n, dose=f, variant="primary", R=R, mutually_disjoint_draws=max_disjoint,
                                n_permutations_spanned=int(np.unique(pids).size),
                                mean_pairwise_overlap_frac=mo, max_pairwise_overlap_frac=mx,
                                frac_pairs_disjoint=fz, expected_overlap_independent=n / N,
                                draw_in_reference_frac=0.0))
            mo, mx, fz = overlap_stats(gd_draws[di], n, N)
            ng = np.array([len(c) for c in gchosen[di]])
            ov_rows.append(dict(depth=n, dose=f, variant="R1_guide", R=R, mutually_disjoint_draws=np.nan,
                                n_permutations_spanned=np.nan,
                                mean_pairwise_overlap_frac=mo, max_pairwise_overlap_frac=mx,
                                frac_pairs_disjoint=fz, expected_overlap_independent=n / N,
                                draw_in_reference_frac=0.0, guides_per_draw_median=float(np.median(ng)),
                                guides_per_draw_min=int(ng.min()), guides_per_draw_max=int(ng.max())))

        flat_p = [d for di in range(len(DOSES)) for d in pd_draws[di]]
        CD = np.asarray((indicator(flat_p, N) @ X_nt).todense())
        E_prim = np.empty((CD.shape[0], CD.shape[1]))
        E_pool = np.empty_like(E_prim)
        for i in range(CD.shape[0]):
            zd = zm(CD[i])
            E_prim[i] = zd - zm(C_pool - CD[i])
            E_pool[i] = zd - z_pool_m
        del CD
        flat_g = [d for di in range(len(DOSES)) for d in gd_draws[di]]
        flat_gc = [c for di in range(len(DOSES)) for c in gchosen[di]]
        CG = np.asarray((indicator(flat_g, N) @ X_nt).todense())
        E_guide = np.empty((CG.shape[0], CG.shape[1]))
        for i in range(CG.shape[0]):
            ref = C_pool - C_nt[flat_gc[i]].sum(0)
            E_guide[i] = zm(CG[i]) - zm(ref)
        del CG
        dose_of = np.repeat(np.arange(len(DOSES)), R)
        draw_of = np.tile(np.arange(R), len(DOSES))
        for s in SPACES:
            r, r2, rn = sp_r[s]["r"], sp_r[s]["r2"], sp_r[s]["rnorm"]
            P = {"primary": proj(E_prim, s), "R2_pool": proj(E_pool, s), "R1_guide": proj(E_guide, s)}
            P["R3_fullgap"] = P["primary"]
            for v in VARIANTS:
                pv = plant[s]["full"] if v == "R3_fullgap" else plant[s]["meas"]
                fvec = np.asarray(DOSES)[dose_of][:, None]
                dl, *_ = block_delta(P[v] + fvec * pv, r, r2, rn)
                for i in range(dl.size):
                    long_rows.append((n, DOSES[dose_of[i]], int(draw_of[i]), v, s, float(dl[i])))
            del P
        del E_prim, E_pool, E_guide
        log(f"[depth {n}] draws done perms={info['n_perm']} ({time.time() - t0:.0f}s)")

    long = pd.DataFrame(long_rows, columns=["depth", "dose", "draw", "variant", "space", "delta"])
    long.to_csv(OUT / "stage1_draws.csv", index=False)
    ov = pd.DataFrame(ov_rows)
    ov.to_csv(OUT / "overlap.csv", index=False)
    design = pd.DataFrame(design_rows)
    design.to_csv(OUT / "draw_design.csv", index=False)

    # Summaries.
    srows = []
    grows = []
    for vi, v in enumerate(VARIANTS):
        for si, s in enumerate(SPACES):
            pv = plant[s]["full"] if v == "R3_fullgap" else plant[s]["meas"]
            r, r2, rn = sp_r[s]["r"], sp_r[s]["r2"], sp_r[s]["rnorm"]
            for n in DEPTHS:
                sub = long[(long.variant == v) & (long.space == s) & (long.depth == n)]
                x0 = sub[sub.dose == 0.0].sort_values("draw")["delta"].to_numpy()
                q05_0 = float(np.percentile(x0, 5))
                meds = []
                powers = {}
                for di, f in enumerate(DOSES):
                    x = sub[sub.dose == f].sort_values("draw")["delta"].to_numpy()
                    rng_m = np.random.default_rng([BOOT_SEED, n, di, si, vi])
                    bi = rng_m.integers(0, x.size, (N_BOOT, x.size))
                    bm = np.median(x[bi], axis=1)
                    nfree = float(block_delta(f * pv, r, r2, rn)[0][0])
                    row = dict(variant=v, space=s, depth=n, dose=f, R=int(x.size),
                               median_delta=float(np.median(x)),
                               median_ci_lo=float(np.percentile(bm, 2.5)),
                               median_ci_hi=float(np.percentile(bm, 97.5)),
                               mean_delta=float(x.mean()), sd_delta=float(x.std(ddof=1)),
                               p05=float(np.percentile(x, 5)), p95=float(np.percentile(x, 95)),
                               noise_free_delta=nfree, q05_dose0=q05_0)
                    if f > 0:
                        k = int(np.sum(x < q05_0))
                        pw = k / x.size
                        rng_p = np.random.default_rng([BOOT_SEED, n, di, si, vi, 1])
                        i0 = rng_p.integers(0, x0.size, (N_BOOT, x0.size))
                        i1 = rng_p.integers(0, x.size, (N_BOOT, x.size))
                        q = np.percentile(x0[i0], 5, axis=1)
                        bp = np.mean(x[i1] < q[:, None], axis=1)
                        u = sps.mannwhitneyu(x, x0, alternative="less")
                        row.update(power=pw, power_ci_lo=float(np.percentile(bp, 2.5)),
                                   power_ci_hi=float(np.percentile(bp, 97.5)),
                                   power_wilson_lo=wilson(k, x.size)[0], power_wilson_hi=wilson(k, x.size)[1],
                                   separable=bool(pw >= POWER_BAR),
                                   mwu_p_less=float(u.pvalue), auc=float(u.statistic / (x.size * x0.size)))
                        powers[f] = pw
                    srows.append(row)
                    meds.append(float(np.median(x)))
                mono = bool(all(meds[i + 1] < meds[i] for i in range(len(meds) - 1)))
                sep25 = powers[0.25] >= POWER_BAR
                smallest = next((f for f in (0.25, 0.50, 1.00) if powers[f] >= POWER_BAR), None)
                grows.append(dict(variant=v, space=s, depth=n, monotone=mono,
                                  power_25=powers[0.25], power_50=powers[0.50], power_100=powers[1.00],
                                  separable_25=sep25,
                                  smallest_separable=("none" if smallest is None else f"{int(smallest * 100)}%"),
                                  qualifies=bool(mono and sep25),
                                  gated=(v == "primary")))
    summ = pd.DataFrame(srows)
    summ.to_csv(OUT / "stage1_summary.csv", index=False)
    gate = pd.DataFrame(grows)
    gate.to_csv(OUT / "stage1_gate.csv", index=False)
    ctx.update(long=long, summ=summ, gate=gate, ov=ov, design=design, gm_ratio=gm_ratio,
               plant_norms={s: (float(np.linalg.norm(plant[s]["meas"])), float(np.linalg.norm(plant[s]["full"])))
                            for s in SPACES})
    return ctx


def r5_retention(ctx, log):
    Zg, age, g = ctx["Zg"], ctx["age"], ctx["g"]
    mid = np.isin(age, ["40-49", "50-59"])
    Zm = Zg[mid]
    sp = fit_pca(Zm - Zm.mean(0), KMAX)
    Lmid = np.asarray(sp.L, np.float64)
    gn = float(np.linalg.norm(g))
    rows = []
    for k in S2_KS:
        rows.append(dict(k=k, n_donors_south6_basis=int(Zg.shape[0]),
                         retention_south6_basis=float(np.linalg.norm(g @ ctx["L_full"][:, :k]) / gn),
                         n_donors_mid_basis=int(mid.sum()),
                         retention_mid_only_basis=float(np.linalg.norm(g @ Lmid[:, :k]) / gn) if Lmid.shape[1] >= k else np.nan))
    df = pd.DataFrame(rows)
    df.to_csv(OUT / "r5_basis_retention.csv", index=False)
    log(f"[R5] {df.to_dict('records')}")
    ctx["r5"] = df
    return ctx


# ---------------------------------------------------------------- findings

def md_table(df, cols, heads=None, d=3):
    heads = heads or cols
    out = ["| " + " | ".join(heads) + " |", "|" + "---|" * len(cols)]
    for r in df.itertuples(index=False):
        vals = []
        for c in cols:
            x = getattr(r, c)
            if isinstance(x, (float, np.floating)):
                vals.append(fmt(float(x), d))
            else:
                vals.append(str(x))
        out.append("| " + " | ".join(vals) + " |")
    return "\n".join(out)


def recon_section(ctx):
    age = ctx.get("age")
    if age is None:
        age = np.asarray(load_gtex_z(lambda *a, **k: None)["age"]).astype(str)
    counts = pd.Series(age).value_counts().to_dict() if age is not None else {}
    ny = int(counts.get("20-29", 0) + counts.get("30-39", 0))
    no = int(counts.get("60-69", 0) + counts.get("70-79", 0))
    nm = int(counts.get("40-49", 0) + counts.get("50-59", 0))
    ntot = int(sum(counts.values()))
    L = []
    L.append("## 1. Reconnaissance: what the SOUTH6 settings are (read-only)")
    L.append("")
    L.append(f"> **The anchors entered the S2 basis.** The PCA behind S2_k20, S2_k50 and S2_k100 was fitted on all "
             f"{ntot} GTEx fibroblast donors. That set includes the {ny} donors aged 20-39 whose mean is Y1 and the "
             f"{no} donors aged 60-79 whose mean is O1: {ny + no} of the {ntot} rows ({100 * (ny + no) / max(ntot, 1):.1f}%) "
             f"are anchor donors. Only the {nm} donors aged 40-59 are not. Gap retention in the S2 spaces is therefore "
             "partly circular. The basis was built from rows that contain the young-old difference it is later asked "
             "to retain. R5 below measures how much retention survives a basis that excludes every anchor donor.")
    L.append("")
    L.append("Sources: `src/south6_run.py` `build_spaces` (lines 1293-1346) and `make_settings` (lines 1349-1402); "
             "`src/toward_run.py` `load_gtex_z` and `freeze_anchors`; `src/lowdim_common.py` `fit_pca`; "
             "`results/south6/south6_run.log`.")
    L.append("")
    L.append("| label | what it is |")
    L.append("|---|---|")
    L.append("| S1 | Gene space. The 23,485 genes of the frozen ruler, in frozen-GTEx z units (log2 CPM minus the ruler's mu, over its sd). "
             "20,265 ruler genes match the Southard panel (20,264 distinct source columns). The other 3,220, and 164 matched genes with zero counts across the panel, sit at z = 0 in every screen vector. |")
    L.append("| S2_k20, S2_k50, S2_k100 | Projection of S1 vectors onto the first k principal components of the GTEx donor matrix. "
             "Fit: `fit_pca(Zg - Zg.mean(0), 100)` on Zg = all 652 GTEx fibroblast donors x 23,485 genes in frozen z, all age bins "
             f"({', '.join(f'{k}: {v}' for k, v in sorted(counts.items()))}). Gram eigendecomposition, loadings scaled to unit norm. "
             "One fit at k = 100; k = 20 and 50 are its leading columns. No Southard cell enters the fit. Origins, targets and screen vectors are projected as (x - Zg.mean(0)) @ L[:, :k]. |")
    L.append("| O1 (primary origin) | `c_old` in `results/toward/anchors.npz`: the mean frozen-z profile of the 231 GTEx donors aged 60-79 (60-69: 210, 70-79: 21). Frozen, not refit. |")
    L.append("| Y1 (primary target) | `c_young` in the same file: the mean of the 113 donors aged 20-39 (20-29: 60, 30-39: 53). Frozen, not refit. |")
    L.append("| O2 (secondary origin) | The Lu aged donor GM00731 day-0 pseudobulk, row 'O' of `results/same/panel_z.npz`. |")
    L.append("| Y2 (secondary target) | A young region in S2 only. Its mean and covariance come from 56 of the 113 young donors (default_rng(20260914) permutation). Its radius is the 95th-percentile Mahalanobis distance of the other 57. Defined at k = 20 and 50; skipped at k = 100 (covariance rank 55 < 100). O1 already sits inside Y2 at both defined k (Mahalanobis 1.06 against radius 11.88 at k = 20; 6.01 against 59.30 at k = 50). |")
    L.append("")
    L.append("Gap retention in the stored settings: ||O1 - Y1|| = 32.284 in S1, 28.645 in S2_k20 (88.7%), 30.864 in S2_k50 (95.6%), "
             "31.434 in S2_k100 (97.4%) (`origin_target_distances.csv`).")
    L.append("")
    L.append("Also shared with the anchors, and not part of the basis question: the frozen mu and sd that define z in every space, S1 included, "
             "come from the same 652 GTEx donors (`frozen_ruler_ridge_raw.npz`, note 'trained on all GTEx fibroblast donors').")
    L.append("")
    L.append("The stored Stage 2C dose response, for contrast with Stage 1 here, is noise-free. It adds f x (median factor norm) along the unit gap direction to O1, "
             "with no sampling noise and including gap coordinates the panel cannot measure. Its MDA of 0.05 therefore says nothing about recovery from real pseudobulks.")
    return "\n".join(L)


def write_findings(ctx, status, stop=None):
    pre = PREREG.read_text(encoding="utf-8")
    L = [f"# FINDINGS_SOUTH6REC — recovery control for the SOUTH6 readout, depth grid", ""]
    L.append(f"Status: {status}")
    L.append("")
    L.append("## Pre-registration (verbatim, `results/south6rec/PREREG.flag`)")
    L.append("")
    L.append("```")
    L.append(pre.rstrip("\n"))
    L.append("```")
    L.append("")
    L.append(recon_section(ctx))
    L.append("")
    L.append("## 2. Stage 0: SOUTH6 reproduced with this task's code (gate)")
    L.append("")
    rep = ctx.get("rep")
    if rep is not None:
        n_g = int(rep.gating.sum())
        n_fail = int(((~rep.passed) & rep.gating).sum())
        L.append(f"{len(rep)} checks, {n_g} gating, {n_fail} gating failures. Tolerance |new - stored| <= 1e-4 x max(1, |stored|); counts exact. "
                 "The full table is `results/south6rec/stage0_reproduction.csv`.")
        L.append("")
        L.append(md_table(rep, ["block", "item", "n", "max_abs_diff", "n_mismatch", "passed"], d=3))
    else:
        L.append("Stage 0 did not complete.")
    L.append("")
    if stop is not None:
        L.append(f"## STOP `{stop.key}`")
        L.append("")
        L.append(stop.message)
        L.append("")
        L.append("No Stage 1 number was computed. Per the pre-registration, the run exits here and is not resumed.")
        FINDINGS.write_text("\n".join(L) + "\n", encoding="utf-8")
        return
    gate, summ, ov, design = ctx["gate"], ctx["summ"], ctx["ov"], ctx["design"]
    L.append("## 3. Stage 1: the gate, per space and depth (primary null, measured-gene plant)")
    L.append("")
    L.append("Plant f x g_m, where g_m is the young-minus-old gap restricted to the coordinates the panel measures. "
             "Its norm relative to the full gap in each space:")
    L.append("")
    L.append("| space | ||g_m|| | ||g|| | ratio |")
    L.append("|---|---|---|---|")
    for s in SPACES:
        a, b = ctx["plant_norms"][s]
        L.append(f"| {s} | {a:.3f} | {b:.3f} | {a / b:.3f} |")
    L.append("")
    gp = gate[gate.variant == "primary"]
    L.append(md_table(gp, ["space", "depth", "monotone", "power_25", "power_50", "power_100", "smallest_separable", "qualifies"],
                      ["space", "cells", "monotone", "power 25%", "power 50%", "power 100%", "smallest separable", "QUALIFIES"]))
    L.append("")
    L.append("Power = fraction of the 200 dose draws whose delta is below the 5th percentile of the 200 dose-0 draws (bar 0.80).")
    L.append("")
    L.append("### Readings that fired (pre-registered wording)")
    L.append("")
    for r in gp.itertuples():
        if r.qualifies:
            L.append(f"- {r.space}, {r.depth} cells: recovered. Median delta falls monotonically with dose and a 25% plant "
                     f"is separated from the 0% draws (power {r.power_25:.3f}).")
        else:
            L.append(f"- {r.space}, {r.depth} cells: not recovered. Monotone: {'yes' if r.monotone else 'no'}; 25% plant "
                     f"power {r.power_25:.3f} against the 0.80 bar. At this depth and in this space, the readout cannot "
                     f"detect an approach of 25% of the gap.")
    L.append("")
    L.append("## 4. Measured change in distance per dose (primary)")
    L.append("")
    L.append("Median delta over 200 draws with its 95% bootstrap CI. Noise-free = delta with no sampling noise (e = 0). "
             "Units are each space's own; compare within a row block, not across spaces.")
    for s in SPACES:
        L.append("")
        L.append(f"### {s}")
        L.append("")
        L.append("| cells | dose | median delta [95% CI] | mean | SD | 5th-95th pct | noise-free | power [95% CI] | MWU p (R4) | AUC (R4) |")
        L.append("|---|---|---|---|---|---|---|---|---|---|")
        sub = summ[(summ.variant == "primary") & (summ.space == s)]
        for r in sub.itertuples():
            pw = "—" if r.dose == 0 else f"{r.power:.3f} [{r.power_ci_lo:.3f}, {r.power_ci_hi:.3f}]"
            mp = "—" if r.dose == 0 else f"{r.mwu_p_less:.2e}"
            au = "—" if r.dose == 0 else f"{r.auc:.3f}"
            L.append(f"| {r.depth} | {int(r.dose * 100)}% | {r.median_delta:.3f} [{r.median_ci_lo:.3f}, {r.median_ci_hi:.3f}] | "
                     f"{r.mean_delta:.3f} | {r.sd_delta:.3f} | {r.p05:.3f} to {r.p95:.3f} | {r.noise_free_delta:.3f} | {pw} | {mp} | {au} |")
    L.append("")
    L.append("## 5. Overlap and variance (stated explicitly)")
    L.append("")
    L.append("Draw design per depth (14,013 NT cells):")
    L.append("")
    L.append(md_table(design, ["depth", "m_disjoint_blocks", "blocks_used_per_perm", "draws_per_dose_per_perm", "permutations_used"],
                      ["cells", "disjoint blocks per permutation", "blocks used", "draws per dose per permutation", "permutations"]))
    L.append("")
    op = ov[ov.variant == "primary"]
    L.append(md_table(op, ["depth", "dose", "R", "mutually_disjoint_draws", "n_permutations_spanned", "mean_pairwise_overlap_frac",
                           "max_pairwise_overlap_frac", "frac_pairs_disjoint", "expected_overlap_independent"],
                      ["cells", "dose", "draws", "mutually disjoint (largest one-permutation set)", "permutations spanned",
                       "mean pairwise overlap / n", "max pairwise overlap / n", "share of pairs disjoint", "n / 14013"]))
    L.append("")
    L.append("Draws within one permutation share no cell, across doses too. Draws in the same condition from different permutations overlap at random, "
             "at about n/14,013 of their cells per pair. At 2,500 cells only one draw per dose fits in a permutation. Each of the 200 draws per dose there "
             "comes from its own permutation, and any two share about 18% of their cells. The primary reference is the complement of the draw, so a draw never shares a cell with its own reference.")
    L.append("")
    L.append("What the replicate spread does and does not include:")
    L.append("")
    L.append("1. Every draw resamples the same 14,013 cells. The spread is the spread of resampling this pool, not of new cells from the line. Given the pool, draws are independent, so overlap does not shrink the spread around the pool mean. But no draw can show variation the pool lacks. At 2,500 cells, 200 draws cover each cell about 36 times per dose.")
    L.append("2. Random cells mix all 78 guides, so the primary null carries no guide-to-guide variation. A factor's cells come from about 6 guides. R1 below restores that structure; compare its spread with the primary.")
    L.append("3. A draw's reference is the other NT cells. A real factor is compared with all 14,013. The primary reference is smaller by n cells, which adds a little noise at 2,500 cells (11,513 reference cells).")
    L.append("")
    L.append("## 6. Report-only readings (pre-registered, no decision rule, never a gate)")
    L.append("")
    for v, title in (("R1_guide", "R1 guide-structured null"), ("R2_pool", "R2 pool reference (draw inside its own reference)"),
                     ("R3_fullgap", "R3 full-gap plant (includes unmeasurable coordinates)")):
        L.append(f"### {title}")
        L.append("")
        gv = gate[gate.variant == v]
        L.append(md_table(gv, ["space", "depth", "monotone", "power_25", "power_50", "power_100", "smallest_separable"],
                          ["space", "cells", "monotone", "power 25%", "power 50%", "power 100%", "smallest separable"]))
        L.append("")
        sv = summ[(summ.variant == v) & (summ.dose == 0.0)][["space", "depth", "median_delta", "sd_delta"]].rename(
            columns={"median_delta": "median_delta_0", "sd_delta": "sd_delta_0"})
        sp_ = summ[(summ.variant == "primary") & (summ.dose == 0.0)][["space", "depth", "sd_delta"]].rename(columns={"sd_delta": "sd_primary_0"})
        mm = sv.merge(sp_, on=["space", "depth"])
        mm["sd_ratio"] = mm.sd_delta_0 / mm.sd_primary_0
        L.append("Dose-0 spread against the primary:")
        L.append("")
        L.append(md_table(mm, ["space", "depth", "median_delta_0", "sd_delta_0", "sd_primary_0", "sd_ratio"],
                          ["space", "cells", "median delta at 0%", "SD at 0%", "primary SD at 0%", "SD ratio"]))
        L.append("")
        if v == "R1_guide":
            o1 = ov[ov.variant == "R1_guide"]
            L.append(md_table(o1, ["depth", "dose", "guides_per_draw_median", "guides_per_draw_min", "guides_per_draw_max",
                                   "mean_pairwise_overlap_frac", "frac_pairs_disjoint"],
                              ["cells", "dose", "guides per draw (median)", "min", "max", "mean pairwise overlap / n", "share of pairs disjoint"]))
            L.append("")
    L.append("### R4 Mann-Whitney")
    L.append("")
    L.append("Per condition in the Section 4 tables (one-sided p that dose draws sit below dose-0 draws, and AUC). With 200 against 200 draws, "
             "this tests a shift of the whole distribution, not whether one factor would be called.")
    L.append("")
    L.append("### R5 Basis circularity: gap retention ||proj_k g|| / ||g||")
    L.append("")
    L.append(md_table(ctx["r5"], ["k", "n_donors_south6_basis", "retention_south6_basis", "n_donors_mid_basis", "retention_mid_only_basis"],
                      ["k", "donors in SOUTH6 basis", "retention, SOUTH6 basis", "donors in 40-59 basis", "retention, 40-59-only basis"]))
    L.append("")
    L.append("## 7. Limitations")
    L.append("")
    L.append("- The plant is added in z after the transform. A real effect would pass through counts, TMM and the log, and would add its own cell-to-cell variance. Here only the null carries sampling noise.")
    L.append("- The plant points exactly along the gap. A real factor's displacement also carries components off the gap axis. In S1, the screen's median ||d|| (123.4) is about 3.8 times the gap (32.3), so off-axis movement dominates real vectors in a way this control does not include.")
    L.append("- The plant is restricted to measured genes. At 100% it closes only the measured part of the gap, so even a noise-free 100% plant leaves the distance above zero (see 'noise-free').")
    L.append("- One NT pool of 14,013 cells from one neonatal line, with no aged cells. Every delta is a counterfactual against GTEx anchors from tissue donors.")
    L.append("- The S2 basis shares donors with the anchors (Section 1, R5). The gate uses that basis as stored.")
    L.append("- Separability is called against the 5th percentile of 200 dose-0 draws. That percentile is itself estimated, and its uncertainty is carried only in the bootstrap CI on power, which is reported but not gated on.")
    L.append("")
    L.append("## 8. Files")
    L.append("")
    L.append("`src/south6rec_run.py`; `results/south6rec/`: `PREREG.flag`, `stage0_reproduction.csv`, `draw_design.csv`, `overlap.csv`, "
             "`stage1_draws.csv` (every draw x space x variant), `stage1_summary.csv`, `stage1_gate.csv`, `r5_basis_retention.csv`, `south6rec_run.log`.")
    FINDINGS.write_text("\n".join(L) + "\n", encoding="utf-8")


def main():
    if not PREREG.exists():
        raise SystemExit(f"PREREG.flag missing at {PREREG}. Not computing.")
    log = Logger(OUT / "south6rec_run.log", mode="a")
    ctx = {}
    try:
        log(f"[prereg] {PREREG} bytes={PREREG.stat().st_size}")
        ctx = stage0(ctx, log)
        ctx = stage1(ctx, log)
        ctx = r5_retention(ctx, log)
        write_findings(ctx, "completed. Stage 0 reproduced SOUTH6; Stage 1 and the report-only readings ran.")
        log("[done] FINDINGS_SOUTH6REC.md written")
    except Stop as err:
        log(f"[STOP] {err.key}: {err.message}")
        if ctx.get("rep") is None and ctx.get("rep_rows") is not None:
            ctx["rep"] = pd.DataFrame(ctx["rep_rows"])
        write_findings(ctx, f"STOP `{err.key}`", stop=err)
    finally:
        log.close()


if __name__ == "__main__":
    main()
