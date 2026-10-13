"""TRANSFER — does a factor's displacement in Southard point the same way as
the same factor's displacement in Sengstack's CRA arm?

The pre-registration is results/transfer/PREREG.flag (2026-10-13T18:44:34Z) and
the author's answers are results/transfer/AMENDMENTS.md. This script reads both
and rewrites neither. Every frozen instrument is read, never refit: the frozen
ruler, the frozen GTEx mu/sd, the 20-component S2A basis of results/south7/
basis.npz, the SENG split and quartiles, and both proliferation axes.

Stages, run from the repository root:
    python src/transfer_run.py stage0     reproduction gate + F + G
    python src/transfer_run.py south      Southard displacement vectors
    python src/transfer_run.py seng       Sengstack displacement vectors
    python src/transfer_run.py read       gate, primary readout, diagnostics
    python src/transfer_run.py findings   FINDINGS_TRANSFER.md
"""
from __future__ import annotations

import json
import sys
import time
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import sparse, stats as sps

SRC = Path(__file__).resolve().parent
if str(SRC) not in sys.path:
    sys.path.insert(0, str(SRC))

import seng_run as S                      # noqa: E402  frozen SENG instrument
import seng2_common as K                  # noqa: E402  frozen SENG2 panel layer
import south6_run as S6                   # noqa: E402  frozen SOUTH6 pipeline
from fibro2_common import load_frozen_ruler  # noqa: E402
from same_run import ci_covers_point      # noqa: E402
from toward_run import percentile_ci, tmm_logcpm_quiet  # noqa: E402

ROOT = SRC.parent
OUT = ROOT / "results" / "transfer"
PREREG = OUT / "PREREG.flag"
AMEND = OUT / "AMENDMENTS.md"
FINDINGS = ROOT / "FINDINGS_TRANSFER.md"
SOUTH6 = ROOT / "results" / "south6"
SOUTH7 = ROOT / "results" / "south7"
SENG = ROOT / "results" / "seng"
SENG2 = ROOT / "results" / "seng2"
SENG3 = ROOT / "results" / "seng3"

SEED = 20260914
N_SPLIT = 100
N_PERM_SHUFFLE = 1000
N_BOOT = 10_000
KDIM = 20
R_FLOOR = 0.3
M_SAME_BAR = 0.85
GATE_P = 0.01
MIN_FQ = 20
TOL = 1e-4                 # |new - stored| <= TOL * max(1, |stored|), SENG3 Stage 0
CRA_PREFIX = "CRA_"
SCREEN_ID = {"south": 1, "seng": 2}
NAMED_N = {"CRA_E2F3": 128, "CRA_EZH2": 107, "CRI_STAT3": 96, "CRI_ZFX": 157}
SENG_CRA_COUNT = 167


def require_prereg():
    for p in (PREREG, AMEND):
        if not p.exists():
            raise SystemExit(f"{p} missing. Not computing.")


def dump(name, obj):
    (OUT / name).write_text(json.dumps(_js(obj), indent=2), encoding="utf-8")


def _js(o):
    if isinstance(o, dict):
        return {str(k): _js(v) for k, v in o.items()}
    if isinstance(o, (list, tuple)):
        return [_js(v) for v in o]
    if isinstance(o, np.integer):
        return int(o)
    if isinstance(o, np.floating):
        return None if not np.isfinite(o) else float(o)
    if isinstance(o, np.bool_):
        return bool(o)
    if isinstance(o, np.ndarray):
        return _js(o.tolist())
    if isinstance(o, float) and not np.isfinite(o):
        return None
    return o


def fmt(x, d=4):
    if x is None:
        return "n/a"
    x = float(x)
    if not np.isfinite(x):
        return "n/a"
    if x != 0.0 and abs(x) < 10 ** (-d):
        return f"{x:.2e}"
    return f"{x:.{d}f}"


class Stop(Exception):
    def __init__(self, key, message):
        super().__init__(message)
        self.key = key
        self.message = message


def write_stop(key, message, extra=None):
    dump("STOP.json", dict(key=key, message=message, extra=extra or {}))


# --------------------------------------------------------------------------
# geometry helpers
# --------------------------------------------------------------------------
def load_basis(log):
    z = np.load(SOUTH7 / "basis.npz", allow_pickle=True)
    L = np.asarray(z["L_mid_full"], np.float64)[:, :KDIM]
    err = float(np.max(np.abs(L.T @ L - np.eye(KDIM))))
    n_mid = int(np.asarray(z["n_mid"]).ravel()[0])
    log(f"[basis] S2A_k{KDIM} read frozen from results/south7/basis.npz L={L.shape} "
        f"n_mid={n_mid} max|LtL-I|={err:.3e}")
    if err > 1e-8:
        raise Stop("basis_not_orthonormal",
                   f"frozen basis columns are not orthonormal (max|LtL-I| = {err:.3e}).")
    return L, n_mid, err


def unit_rows(A):
    A = np.atleast_2d(np.asarray(A, np.float64))
    n = np.linalg.norm(A, axis=1)
    out = np.zeros_like(A)
    ok = n > 0
    out[ok] = A[ok] / n[ok, None]
    return out


def cos_rows(A, B):
    return np.einsum("ij,ij->i", unit_rows(A), unit_rows(B))


def spearman_brown(c_half):
    c = np.asarray(c_half, float)
    r = np.full(c.shape, np.nan)
    ok = c > 0
    r[ok] = 2.0 * c[ok] / (1.0 + c[ok])
    return r


def halves(rng, n):
    """Disjoint halves of range(n); the extra cell goes to half A."""
    perm = rng.permutation(n)
    k = (n + 1) // 2
    return perm[:k], perm[k:]


# --------------------------------------------------------------------------
# Southard side
# --------------------------------------------------------------------------
def south_ctx(log, need_cache=True):
    s1 = np.load(SOUTH6 / "stage1.npz", allow_pickle=True)
    keep = [str(x) for x in s1["keep"]]
    matched_ruler = np.asarray(s1["matched_ruler"], np.int64)
    missing = np.asarray(s1["missing"], bool)
    nnz = int(np.asarray(s1["nnz"]).ravel()[0])
    n_ruler = int(np.asarray(s1["n_ruler"]).ravel()[0])
    cell_flat = np.asarray(s1["cell_flat"], np.int64)
    cell_indptr = np.asarray(s1["cell_indptr"], np.int64)
    nt_cells = np.asarray(s1["nt_cells"], np.int64)
    data, indices, indptr = S6.open_matched(nnz)
    frozen = load_frozen_ruler()
    mu = np.asarray(frozen["mu"], float)[matched_ruler]
    sd = np.asarray(frozen["sd"], float)[matched_ruler]
    ctx = dict(s1=s1, keep=keep, matched_ruler=matched_ruler, missing=missing,
               nnz=nnz, n_ruler=n_ruler, cell_flat=cell_flat, cell_indptr=cell_indptr,
               nt_cells=nt_cells, data=data, indices=indices, indptr=indptr,
               frozen=frozen, mu=mu, sd=sd, n_matched=int(matched_ruler.size))
    if need_cache:
        c = np.load(OUT / "south_cache.npz", allow_pickle=True)
        ctx["cache"] = dict(ref=int(np.asarray(c["ref"]).ravel()[0]),
                            ref_counts=np.asarray(c["ref_counts"], float),
                            ref_lib=float(np.asarray(c["ref_lib"]).ravel()[0]),
                            geom=float(np.asarray(c["geom"]).ravel()[0]),
                            ave=float(np.asarray(c["ave"]).ravel()[0]),
                            lib=np.asarray(c["lib"], float), nf=np.asarray(c["nf"], float))
        ctx["C_pool"] = np.asarray(c["C_pool"], float)
    log(f"[south] factors={len(keep)} matched={ctx['n_matched']} nt_cells={nt_cells.size} "
        f"nnz={nnz:,}")
    return ctx


def south_sum(ctx, cells):
    return S6.gather_sum(ctx["data"], ctx["indices"], ctx["indptr"],
                         np.asarray(cells, np.int64), ctx["n_matched"])


def south_z(ctx, counts_matched):
    z, _, _, _ = S.z_from_counts(np.asarray(counts_matched, np.float64).ravel(),
                                 ctx["cache"], ctx["mu"], ctx["sd"], ctx["missing"])
    return z


def south_cells(ctx, i):
    a, b = int(ctx["cell_indptr"][i]), int(ctx["cell_indptr"][i + 1])
    return ctx["cell_flat"][a:b]


# --------------------------------------------------------------------------
# Sengstack side
# --------------------------------------------------------------------------
def seng_ctx(log):
    size = K.H5.stat().st_size
    md5 = S.file_md5(K.H5)
    log(f"[seng file] size={size} md5={md5}")
    if md5 != K.MD5_EXPECT:
        raise Stop("file_mismatch",
                   f"{K.H5} md5 {md5}, expected {K.MD5_EXPECT}.")
    ctx = K.load_context(log)
    z = np.load(SENG2 / "panel_counts.npz", allow_pickle=True)
    names = [str(x) for x in z["names"]]
    ni = {n: i for i, n in enumerate(names)}
    C = np.asarray(z["C"], np.float64)
    frozen = ctx["frozen"]
    Z, missing, nf, Cc = K.panel_z_from_counts(C, frozen)
    cache = S.tmm_cache_from_panel(Cc, nf)
    mu = np.asarray(frozen["mu"], float)
    sd = np.asarray(frozen["sd"], float)
    ctx.update(names=names, ni=ni, C=C, Z=Z, missing=np.asarray(missing, bool),
               cache=cache, mu=mu, sd=sd, md5=md5, size=size)
    log(f"[seng] panel rows={len(names)} ruler={ctx['n_ruler']} present={ctx['dst'].size} "
        f"qual={len(ctx['qual'])}")
    return ctx


def seng_z(ctx, counts_ruler):
    z, _, _, _ = S.z_from_counts(np.asarray(counts_ruler, np.float64).ravel(),
                                 ctx["cache"], ctx["mu"], ctx["sd"], ctx["missing"])
    return z


# ==========================================================================
# STAGE 0
# ==========================================================================
def stage0():
    log = K.Log(OUT / "stage0.log")
    rows = []
    t_start = time.time()

    def num(block, item, stored, got, tol=TOL):
        stored = float(stored)
        got = float(got)
        ok = bool(np.isfinite(got) and abs(got - stored) <= tol * max(1.0, abs(stored)))
        rows.append(dict(block=block, item=item, kind="numeric", n=1, n_mismatch=int(not ok),
                         stored=stored, got=got,
                         max_scaled_diff=abs(got - stored) / max(1.0, abs(stored)), passed=ok))
        return ok

    def arr(block, item, stored, got, tol=TOL):
        a = np.asarray(stored, float).ravel()
        b = np.asarray(got, float).ravel()
        if a.shape != b.shape:
            rows.append(dict(block=block, item=item, kind="array", n=a.size,
                             n_mismatch=a.size, stored=np.nan, got=np.nan,
                             max_scaled_diff=np.inf, passed=False))
            return False
        d = np.abs(b - a) / np.maximum(1.0, np.abs(a))
        nm = int((d > tol).sum())
        rows.append(dict(block=block, item=item, kind="array", n=a.size, n_mismatch=nm,
                         stored=np.nan, got=np.nan, max_scaled_diff=float(d.max() if d.size else 0.0),
                         passed=bool(nm == 0)))
        return nm == 0

    def cnt(block, item, stored, got):
        ok = bool(int(stored) == int(got))
        rows.append(dict(block=block, item=item, kind="count", n=1, n_mismatch=int(not ok),
                         stored=float(stored), got=float(got),
                         max_scaled_diff=0.0 if ok else np.inf, passed=ok))
        return ok

    L, n_mid, berr = load_basis(log)

    # ---------------- (a) Southard: re-accumulate and rebuild the transform
    sc = south_ctx(log, need_cache=False)
    s1 = sc["s1"]
    n_f = len(sc["keep"])
    n_m = sc["n_matched"]
    log("[s0a] re-accumulating C_f from the matched CSR cache")
    C_f = np.zeros((n_f, n_m), np.float64)
    t0 = time.time()
    for i in range(n_f):
        C_f[i] = south_sum(sc, south_cells(sc, i))
        if (i + 1) % 400 == 0:
            log(f"[s0a] C_f {i + 1}/{n_f} ({time.time() - t0:.0f}s)")
    C_pool = south_sum(sc, sc["nt_cells"])
    stored_Cf = np.asarray(s1["C_f"], np.float64)
    bad = int((C_f != stored_Cf).sum())
    rows.append(dict(block="(a) Southard counts", item="C_f re-accumulated from cells",
                     kind="array", n=int(C_f.size), n_mismatch=bad, stored=np.nan, got=np.nan,
                     max_scaled_diff=0.0 if bad == 0 else float(np.abs(C_f - stored_Cf).max()),
                     passed=bool(bad == 0)))
    del stored_Cf
    C_nt = np.asarray(s1["C_nt"], np.float64)
    stored_pool = np.asarray(s1["C_pool"], np.float64)
    cnt("(a) Southard counts", "C_pool equals stored C_pool exactly",
        0, int((C_pool != stored_pool).sum()))
    cnt("(a) Southard counts", "C_pool equals the sum of the 78 stored C_nt rows",
        0, int((C_pool != C_nt.sum(0)).sum()))

    log("[s0a] rebuilding the panel transform")
    C_panel = np.vstack([C_f, C_nt, C_pool.reshape(1, -1)])
    del C_f
    logcpm, nf = tmm_logcpm_quiet(C_panel)
    mu, sd = sc["mu"], sc["sd"]
    sd_safe = np.where(sd < 1e-12, 1.0, sd)
    missing = (C_panel == 0).all(0)
    cnt("(b) Southard transform", "missing-gene mask size",
        int(np.asarray(s1["missing"]).sum()), int(missing.sum()))
    cnt("(b) Southard transform", "missing-gene mask identical",
        0, int((missing != np.asarray(s1["missing"], bool)).sum()))
    cache = S.tmm_cache_from_panel(C_panel, nf)
    del C_panel
    Z = (logcpm - mu) / sd_safe
    Z[:, missing] = 0.0
    del logcpm
    z_pool_m = Z[-1].copy()
    stored_zpool = np.asarray(s1["z_pool"], float)
    arr("(b) Southard transform", "z_pool on the matched genes",
        stored_zpool[sc["matched_ruler"]], z_pool_m)
    unm = np.ones(sc["n_ruler"], bool)
    unm[sc["matched_ruler"]] = False
    cnt("(b) Southard transform", "z_pool is 0 on the unmatched ruler genes",
        0, int((stored_zpool[unm] != 0).sum()))
    stored_df = np.asarray(s1["d_f"], float)
    nm_d, max_d = 0, 0.0
    for a in range(0, n_f, 200):
        b = min(a + 200, n_f)
        blk = np.abs((Z[a:b] - z_pool_m) - stored_df[a:b][:, sc["matched_ruler"]])
        blk /= np.maximum(1.0, np.abs(stored_df[a:b][:, sc["matched_ruler"]]))
        nm_d += int((blk > TOL).sum())
        max_d = max(max_d, float(blk.max()))
        del blk
    rows.append(dict(block="(b) Southard transform", item="d_f on the matched genes",
                     kind="array", n=int(n_f * n_m), n_mismatch=nm_d,
                     stored=np.nan, got=np.nan, max_scaled_diff=max_d, passed=bool(nm_d == 0)))
    cnt("(b) Southard transform", "d_f is 0 on the unmatched ruler genes",
        0, int((stored_df[:, unm] != 0).sum()))
    del stored_df, Z
    np.savez(OUT / "south_cache.npz", ref=np.array([cache["ref"]]),
             ref_counts=cache["ref_counts"], ref_lib=np.array([cache["ref_lib"]]),
             geom=np.array([cache["geom"]]), ave=np.array([cache["ave"]]),
             lib=cache["lib"], nf=cache["nf"], C_pool=C_pool,
             matched_ruler=sc["matched_ruler"], missing=missing)
    log(f"[s0a] panel cache written ({time.time() - t0:.0f}s)")

    # ---------------- (b) Sengstack: count, n_cells, and the compressed-space values
    jc = seng_ctx(log)
    Zj, ni = jc["Z"], jc["ni"]
    qual = jc["qual"]
    cra = [(lab, int(idx.size)) for m, lab, idx in qual if m == "CRA"]
    cnt("(c) Sengstack", f"qualifying CRA perturbations (n >= {K.MIN_CELLS})",
        SENG_CRA_COUNT, len(cra))
    d4 = pd.read_csv(SENG3 / "d4_table.csv")
    d4c = d4[d4.modality == "CRA"].set_index("perturbation")
    labs = [l for l, _ in cra]
    miss_lab = [l for l in labs if l not in d4c.index]
    cnt("(c) Sengstack", "every CRA perturbation is a row of seng3/d4_table.csv",
        0, len(miss_lab))
    if not miss_lab:
        arr("(c) Sengstack", "n_cells per CRA perturbation",
            d4c.loc[labs, "n_cells"].to_numpy(float),
            np.array([n for _, n in cra], float), tol=0.0)
        r_i = (Zj[ni["NT_CRA"]] - Zj[ni["EARLY_A"]]) @ L
        r_ii = (Zj[ni["LATE_A"]] - Zj[ni["EARLY_A"]]) @ L
        D = np.vstack([(Zj[ni[l]] - Zj[ni["NT_CRA"]]) for l in labs]) @ L
        for tag, r in (("i", r_i), ("ii", r_ii)):
            r2 = float(r @ r)
            rn = float(np.sqrt(r2))
            dots = D @ r
            qn2 = np.einsum("ij,ij->i", D, D)
            delta = np.sqrt(np.maximum(r2 + 2.0 * dots + qn2, 0.0)) - rn
            frac = -dots / r2
            arr("(c) Sengstack", f"delta_{tag} in S2A_k20, CRA rows",
                d4c.loc[labs, f"delta_{tag}"].to_numpy(float), delta)
            arr("(c) Sengstack", f"frac_{tag} in S2A_k20, CRA rows",
                d4c.loc[labs, f"frac_{tag}"].to_numpy(float), frac)

    # ---------------- (c) the four named cell counts
    sizes = {lab: int(idx.size) for m, lab, idx in qual}
    for lab, n in NAMED_N.items():
        cnt("(d) assignment rule", f"n_cells[{lab}]", n, sizes.get(lab, -1))

    # ---------------- F and G
    keep_set = {k.upper(): k for k in sc["keep"]}
    cra_sym = {l[len(CRA_PREFIX):].upper(): l for l, _ in cra}
    shared = sorted(set(keep_set) & set(cra_sym))
    s_only_n = len(set(keep_set) - set(cra_sym))
    j_only = sorted(set(cra_sym) - set(keep_set))
    G = np.intersect1d(sc["matched_ruler"], np.asarray(jc["dst"], np.int64))
    log(f"[F] shared={len(shared)} southard_only={s_only_n} sengstack_cra_only={len(j_only)}")
    log(f"[G] |G|={G.size} of ruler {sc['n_ruler']}")
    if not shared:
        raise Stop("no_shared_factors", "no gene symbol matched exactly between the screens.")

    fs7 = pd.read_csv(SOUTH7 / "factor_set.csv")
    in_active = set(fs7.loc[fs7.in_active_212.astype(str).str.upper() == "TRUE", "factor"])
    in_exp = set(fs7.loc[fs7.in_expanded_102.astype(str).str.upper() == "TRUE", "factor"])
    s_i = {k.upper(): i for i, k in enumerate(sc["keep"])}
    cell_n = np.diff(sc["cell_indptr"])
    jn = {l: n for l, n in cra}
    frame = pd.DataFrame([dict(
        symbol=s, south_label=keep_set[s], seng_label=cra_sym[s],
        n_S=int(cell_n[s_i[s]]), n_J=int(jn[cra_sym[s]]),
        in_active_212=bool(keep_set[s] in in_active),
        in_expanded_102=bool(keep_set[s] in in_exp)) for s in shared])
    frame.to_csv(OUT / "factor_set.csv", index=False)
    np.savez(OUT / "shared_genes.npz", G=G,
             pos_in_matched=np.searchsorted(sc["matched_ruler"], G))
    s_missing_ruler = sc["matched_ruler"][np.asarray(missing, bool)]
    j_missing_ruler = np.flatnonzero(jc["missing"])
    dump("stage0_sets.json", dict(
        n_shared=len(shared), n_southard_only=s_only_n, n_seng_cra_only=len(j_only),
        seng_cra_only_names=j_only, shared_symbols=shared,
        n_G=int(G.size), n_ruler=int(sc["n_ruler"]),
        n_southard_matched=int(sc["n_matched"]), n_seng_present=int(jc["dst"].size),
        n_G_in_southard_missing=int(np.intersect1d(G, s_missing_ruler).size),
        n_G_in_seng_missing=int(np.intersect1d(G, j_missing_ruler).size),
        basis_n_mid=n_mid, basis_orthonormality_err=berr,
        seng_md5=jc["md5"], seng_size=jc["size"],
        southard_h5ad_opened=False,
        n_S_min=int(frame.n_S.min()), n_S_median=float(frame.n_S.median()),
        n_S_max=int(frame.n_S.max()),
        n_J_min=int(frame.n_J.min()), n_J_median=float(frame.n_J.median()),
        n_J_max=int(frame.n_J.max())))

    rep = pd.DataFrame(rows)
    rep.to_csv(OUT / "stage0_reproduction.csv", index=False)
    failed = rep[~rep.passed]
    dump("stage0_summary.json", dict(passed=bool(len(failed) == 0), n_checks=len(rep),
                                     n_failed=len(failed),
                                     failed=[f"{r.block} / {r.item}" for r in failed.itertuples()],
                                     runtime_s=time.time() - t_start))
    log(f"[stage0] checks={len(rep)} failed={len(failed)} ({time.time() - t_start:.0f}s)")
    jc["f"].close()
    if len(failed):
        msg = "; ".join(f"{r.block} / {r.item}: {r.n_mismatch} of {r.n} outside tol "
                        f"(max scaled {r.max_scaled_diff})" for r in failed.itertuples())
        write_stop("prior_results_not_reproduced", msg)
        log(f"[stage0] STOP prior_results_not_reproduced :: {msg}")
        log.close()
        raise SystemExit(3)
    log.close()


# ==========================================================================
# DISPLACEMENT VECTORS, ONE SCREEN AT A TIME
# ==========================================================================
def build(screen):
    log = K.Log(OUT / f"{screen}.log")
    t_start = time.time()
    L, _, _ = load_basis(log)
    g = np.load(OUT / "shared_genes.npz")
    G = np.asarray(g["G"], np.int64)
    LG = L[G, :]
    frame = pd.read_csv(OUT / "factor_set.csv")
    sid = SCREEN_ID[screen]
    nf = len(frame)

    if screen == "south":
        ctx = south_ctx(log)
        pos = np.asarray(g["pos_in_matched"], np.int64)
        keep_i = {k.upper(): i for i, k in enumerate(ctx["keep"])}
        fac_cells = [south_cells(ctx, keep_i[s.upper()]) for s in frame.symbol]
        ctrl = ctx["nt_cells"]
        C_ctrl = np.asarray(ctx["C_pool"], float)
        C_fac = [south_sum(ctx, c) for c in fac_cells]
        zc = lambda v: south_z(ctx, v)[pos]
        zfull = lambda v: south_z(ctx, v)
        sum_cells = lambda cells: south_sum(ctx, cells)
        ctrl_n = ctrl.size
        ctrl_take = lambda take: south_sum(ctx, ctrl[take])
        p_axis_full = np.asarray(ctx["s1"]["p_unit"], float).ravel()
        p_axis_G = p_axis_full[G] if p_axis_full.size == ctx["n_ruler"] else np.zeros(G.size)
        p_src = "results/south6/stage1.npz key p_unit"
    else:
        ctx = seng_ctx(log)
        n_ruler = ctx["n_ruler"]
        dst = np.asarray(ctx["dst"], np.int64)
        lab_idx = {lab: idx for m, lab, idx in ctx["qual"]}
        fac_cells = [np.asarray(lab_idx[l], int) for l in frame.seng_label]
        ctrl_cells = np.asarray(ctx["nt_mod"]["CRA"], int)
        log(f"[seng] materialising NT_CRA cells={ctrl_cells.size} and "
            f"{int(sum(c.size for c in fac_cells))} factor cells")
        X_ctrl = K.submatrix(ctx["f"], ctrl_cells, ctx["exp"]).tocsr()
        flat = np.concatenate(fac_cells)
        X_fac = K.submatrix(ctx["f"], flat, ctx["exp"]).tocsr()
        off = np.r_[0, np.cumsum([c.size for c in fac_cells])]
        to_r = lambda loc: K.to_ruler(loc, dst, n_ruler)
        C_ctrl = to_r(np.asarray(X_ctrl.sum(0)).ravel())
        C_fac = [to_r(np.asarray(X_fac[off[i]:off[i + 1]].sum(0)).ravel()) for i in range(nf)]
        zc = lambda v: seng_z(ctx, v)[G]
        zfull = lambda v: seng_z(ctx, v)
        ctrl_n = ctrl_cells.size
        ctrl_take = lambda take: to_r(np.asarray(X_ctrl[np.asarray(take, int)].sum(0)).ravel())
        def fac_take(i, take):
            rows = off[i] + np.asarray(take, int)
            return to_r(np.asarray(X_fac[rows].sum(0)).ravel())
        p_axis_full = np.asarray(np.load(SENG / "axes.npz")["p_CRA"], float).ravel()
        p_axis_G = p_axis_full[G]
        p_src = "results/seng/axes.npz key p_CRA"

    if screen == "south":
        def fac_take(i, take):
            return south_sum(ctx, fac_cells[i][np.asarray(take, int)])

    n_cells = np.array([c.size for c in fac_cells], int)

    # ---- full-depth displacement
    zctrl_G = zc(C_ctrl)
    zctrl_full = zfull(C_ctrl)
    dG = np.vstack([zc(C_fac[i]) - zctrl_G for i in range(nf)])
    dfull_norm = np.array([float(np.linalg.norm(zfull(C_fac[i]) - zctrl_full)) for i in range(nf)])
    dP = dG @ LG
    log(f"[{screen}] full-depth displacements built ({time.time() - t_start:.0f}s)")

    # ---- reliability splits (100), and the gate's quarter splits
    hA = np.zeros((nf, N_SPLIT, KDIM))
    hB = np.zeros((nf, N_SPLIT, KDIM))
    gq = np.zeros((nf, N_SPLIT, 4, KDIM))
    # gate contexts: one fixed split of factor and control cells
    rng_g = np.random.default_rng([SEED, sid, 7])
    cgA, cgB = halves(rng_g, ctrl_n)
    Cg = [ctrl_take(cgA), None]
    Cg[1] = C_ctrl - Cg[0]
    zCg = [zc(Cg[0]), zc(Cg[1])]
    g1 = np.zeros((nf, KDIM))
    g2 = np.zeros((nf, KDIM))
    fac_gh = []
    for i in range(nf):
        r = np.random.default_rng([SEED, sid, 7, i])
        a, b = halves(r, n_cells[i])
        fac_gh.append((a, b))
        ca = fac_take(i, a)
        cb = C_fac[i] - ca
        g1[i] = (zc(ca) - zCg[0]) @ LG
        g2[i] = (zc(cb) - zCg[1]) @ LG

    for s in range(N_SPLIT):
        rc = np.random.default_rng([SEED, sid, s])
        cA, _ = halves(rc, ctrl_n)
        CcA = ctrl_take(cA)
        CcB = C_ctrl - CcA
        zA, zB = zc(CcA), zc(CcB)
        # quarter splits of the two gate control halves
        rq = np.random.default_rng([SEED, sid, 8, s])
        q = []
        for h in (0, 1):
            cq = Cg[h]
            members = cgA if h == 0 else cgB
            take, _ = halves(rq, int(members.size))
            part = ctrl_take(members[take])
            q.append((zc(part), zc(cq - part)))
        for i in range(nf):
            r = np.random.default_rng([SEED, sid, i, s])
            a, _ = halves(r, n_cells[i])
            ca = fac_take(i, a)
            cb = C_fac[i] - ca
            hA[i, s] = (zc(ca) - zA) @ LG
            hB[i, s] = (zc(cb) - zB) @ LG
            rqf = np.random.default_rng([SEED, sid, 8, i, s])
            for h in (0, 1):
                cells_h = fac_gh[i][h]
                ch = fac_take(i, cells_h)
                t, _ = halves(rqf, cells_h.size)
                p1 = fac_take(i, cells_h[t])
                p2 = ch - p1
                gq[i, s, 2 * h] = (zc(p1) - q[h][0]) @ LG
                gq[i, s, 2 * h + 1] = (zc(p2) - q[h][1]) @ LG
        if (s + 1) % 10 == 0:
            log(f"[{screen}] splits {s + 1}/{N_SPLIT} ({time.time() - t_start:.0f}s)")

    pG = p_axis_G @ LG
    np.savez(OUT / f"vec_{screen}.npz", symbol=np.array(frame.symbol, dtype=object),
             n_cells=n_cells, dG=dG, dP=dP, dfull_norm=dfull_norm,
             hA=hA, hB=hB, g1=g1, g2=g2, gq=gq, p_proj=pG,
             p_axis_norm_G=np.array([float(np.linalg.norm(p_axis_G))]))
    dump(f"{screen}_build.json", dict(screen=screen, n_factors=nf, n_G=int(G.size),
                                      ctrl_cells=int(ctrl_n), prolif_axis_source=p_src,
                                      prolif_axis_norm_in_space=float(np.linalg.norm(pG)),
                                      runtime_s=time.time() - t_start))
    log(f"[{screen}] done ({time.time() - t_start:.0f}s)")
    if screen == "seng":
        ctx["f"].close()
    log.close()


# ==========================================================================
# READOUT
# ==========================================================================
def readout(DS, DJ, hA_S, hB_S, hA_J, hB_J):
    """The whole readout on two sets of projected displacements."""
    c_obs = cos_rows(DS, DJ)
    cS = np.array([np.mean(cos_rows(hA_S[i], hB_S[i])) for i in range(DS.shape[0])])
    cJ = np.array([np.mean(cos_rows(hA_J[i], hB_J[i])) for i in range(DS.shape[0])])
    rS, rJ = spearman_brown(cS), spearman_brown(cJ)
    with np.errstate(invalid="ignore"):
        c_corr = c_obs / np.sqrt(rS * rJ)
    q = np.isfinite(rS) & np.isfinite(rJ) & (rS >= R_FLOOR) & (rJ >= R_FLOOR)
    out = dict(c_obs=c_obs, c_half_S=cS, c_half_J=cJ, r_S=rS, r_J=rJ,
               c_corr=c_corr, qual=q)
    out.update(summarise(c_obs, c_corr, q, DS, DJ))
    return out


def summarise(c_obs, c_corr, q, DS, DJ):
    idx = np.flatnonzero(q)
    nq = idx.size
    res = dict(n_q=int(nq))
    spec = np.zeros(c_obs.shape, bool)
    p95 = np.full(c_obs.shape, np.nan)
    if nq >= 2:
        A = unit_rows(DS[idx])
        B = unit_rows(DJ[idx])
        Cm = A @ B.T
        for j in range(nq):
            off = np.delete(Cm[j], j)
            p95[idx[j]] = float(np.percentile(off, 95))
            spec[idx[j]] = bool(Cm[j, j] > p95[idx[j]])
        k = int(spec[idx].sum())
        bt = sps.binomtest(k, nq, 0.05, alternative="greater")
        res.update(k_specific=k, chance=0.05 * nq, spec_p=float(bt.pvalue))
    else:
        res.update(k_specific=0, chance=0.05 * nq, spec_p=np.nan)
    vals = c_corr[idx]
    fin = np.isfinite(vals)
    res["n_nonfinite_ccorr"] = int((~fin).sum())
    if fin.sum() >= 1:
        v = vals[fin]
        M = float(np.median(v))
        rng = np.random.default_rng(SEED)
        bs = np.array([np.median(v[rng.integers(0, v.size, v.size)]) for _ in range(N_BOOT)])
        lo, hi = percentile_ci(bs)
        res.update(M=M, ci_lo=lo, ci_hi=hi, ci_valid=bool(ci_covers_point(M, lo, hi)),
                   angle=float(np.degrees(np.arccos(np.clip(M, -1, 1)))),
                   median_c_obs=float(np.median(c_obs[idx][fin])))
    else:
        res.update(M=np.nan, ci_lo=np.nan, ci_hi=np.nan, ci_valid=False, angle=np.nan,
                   median_c_obs=np.nan)
    res["specific"] = spec
    res["p95"] = p95
    return res


def read():
    log = K.Log(OUT / "read.log")
    frame = pd.read_csv(OUT / "factor_set.csv")
    vs = np.load(OUT / "vec_south.npz", allow_pickle=True)
    vj = np.load(OUT / "vec_seng.npz", allow_pickle=True)
    if [str(x) for x in vs["symbol"]] != [str(x) for x in vj["symbol"]]:
        raise Stop("prior_results_not_reproduced", "factor order differs between the two vector files.")
    g = np.load(OUT / "shared_genes.npz")
    nG = int(np.asarray(g["G"]).size)

    # ---------------- Stage 4 gate, each screen separately
    gate = {}
    for tag, v in (("south", vs), ("seng", vj)):
        g1, g2, gq = v["g1"], v["g2"], v["gq"]
        r = readout(g1, g2, gq[:, :, 0, :], gq[:, :, 1, :], gq[:, :, 2, :], gq[:, :, 3, :])
        gate[tag] = dict(M_same=r["M"], ci_lo=r["ci_lo"], ci_hi=r["ci_hi"],
                         ci_valid=r["ci_valid"], n_q=r["n_q"], k_specific=r["k_specific"],
                         chance=r["chance"], spec_p=r["spec_p"],
                         median_c_obs=r["median_c_obs"],
                         passed=bool(np.isfinite(r["M"]) and r["M"] >= M_SAME_BAR
                                     and np.isfinite(r["spec_p"]) and r["spec_p"] < GATE_P
                                     and r["n_q"] > 0))
        pd.DataFrame(dict(symbol=frame.symbol, n_cells=v["n_cells"], c_obs=r["c_obs"],
                          r_1=r["r_S"], r_2=r["r_J"], c_corr=r["c_corr"],
                          qualifying=r["qual"], specific=r["specific"])).to_csv(
            OUT / f"gate_{tag}.csv", index=False)
        log(f"[gate {tag}] M_same={fmt(r['M'])} CI=({fmt(r['ci_lo'])},{fmt(r['ci_hi'])}) "
            f"n_q={r['n_q']} specific={r['k_specific']} chance={fmt(r['chance'],2)} "
            f"p={fmt(r['spec_p'],6)} passed={gate[tag]['passed']}")
    gate_passed = bool(gate["south"]["passed"] and gate["seng"]["passed"])
    dump("gate_summary.json", dict(gate=gate, passed=gate_passed, bar_M=M_SAME_BAR,
                                   bar_p=GATE_P))

    # The flag, section 4: if the gate fails, STOP `unmeasurable`, the primary
    # statistic is not computed and no cross-screen number is reported.
    if not gate_passed:
        failed = [k for k, v in gate.items() if not v["passed"]]
        msg = ("the same-context control did not recover zero-drift agreement in "
               f"{' and '.join(failed)} (bars: M_same >= {M_SAME_BAR}, specificity "
               f"p < {GATE_P})")
        write_stop("unmeasurable", msg)
        dump("summary.json", dict(reading="UNMEASURABLE", stop_key="unmeasurable",
                                  stop_message=msg, gate=gate, gate_passed=False,
                                  gate_failed_in=failed, primary=None, diagnostics=None,
                                  n_F=int(len(frame)), n_G=nG, min_fq=MIN_FQ,
                                  bar_M=M_SAME_BAR, bar_p=GATE_P))
        log(f"[reading] UNMEASURABLE :: STOP unmeasurable :: {msg}")
        log("[stop] the primary statistic is not computed and no cross-screen number "
            "is reported, per flag section 4")
        log.close()
        return

    # ---------------- primary readout
    prim = readout(vs["dP"], vj["dP"], vs["hA"], vs["hB"], vj["hA"], vj["hB"])
    tab = pd.DataFrame(dict(
        symbol=frame.symbol, south_label=frame.south_label, seng_label=frame.seng_label,
        n_S=frame.n_S, n_J=frame.n_J, in_active_212=frame.in_active_212,
        in_expanded_102=frame.in_expanded_102,
        c_half_S=prim["c_half_S"], c_half_J=prim["c_half_J"],
        r_S=prim["r_S"], r_J=prim["r_J"], c_obs=prim["c_obs"],
        angle_deg=np.degrees(np.arccos(np.clip(prim["c_obs"], -1, 1))),
        c_corr=prim["c_corr"], qualifying=prim["qual"],
        mismatch_p95=prim["p95"], specific=prim["specific"],
        d_norm_S_full=vs["dfull_norm"], d_norm_J_full=vj["dfull_norm"],
        d_norm_S_masked=np.linalg.norm(vs["dG"], axis=1),
        d_norm_J_masked=np.linalg.norm(vj["dG"], axis=1),
        d_norm_S_proj=np.linalg.norm(vs["dP"], axis=1),
        d_norm_J_proj=np.linalg.norm(vj["dP"], axis=1)))
    tab["share_norm2_kept_by_mask_S"] = (tab.d_norm_S_masked / tab.d_norm_S_full) ** 2
    tab["share_norm2_kept_by_mask_J"] = (tab.d_norm_J_masked / tab.d_norm_J_full) ** 2
    tab["share_norm2_kept_by_proj_S"] = (tab.d_norm_S_proj / tab.d_norm_S_masked) ** 2
    tab["share_norm2_kept_by_proj_J"] = (tab.d_norm_J_proj / tab.d_norm_J_masked) ** 2
    tab["c_obs_full_genes"] = cos_rows(vs["dG"], vj["dG"])
    tab["angle_full_genes_deg"] = np.degrees(np.arccos(np.clip(tab.c_obs_full_genes, -1, 1)))
    tab = tab.sort_values("c_corr", ascending=False, na_position="last")
    tab["rank_by_c_corr"] = np.arange(1, len(tab) + 1)
    tab.to_csv(OUT / "per_factor.csv", index=False)
    log(f"[primary] n_q={prim['n_q']} M={fmt(prim['M'])} "
        f"CI=({fmt(prim['ci_lo'])},{fmt(prim['ci_hi'])}) valid={prim['ci_valid']} "
        f"specific={prim['k_specific']} chance={fmt(prim['chance'],2)} p={fmt(prim['spec_p'],6)}")

    # ---------------- the reading
    nq = prim["n_q"]
    if not gate_passed or nq < MIN_FQ:
        reading = "UNMEASURABLE"
    elif (prim["ci_valid"] and np.isfinite(prim["ci_hi"]) and prim["ci_hi"] < 0.5) or \
         (not np.isfinite(prim["spec_p"])) or prim["spec_p"] >= 0.05:
        reading = "LARGE DRIFT"
    elif prim["ci_valid"] and np.isfinite(prim["ci_lo"]) and prim["ci_lo"] >= 0.7 \
            and prim["spec_p"] < GATE_P:
        reading = "SMALL DRIFT"
    else:
        reading = "INTERMEDIATE"
    log(f"[reading] {reading}")

    # ---------------- diagnostics
    diag = {}
    idx = np.flatnonzero(prim["qual"])
    DS, DJ = np.asarray(vs["dP"]), np.asarray(vj["dP"])
    if idx.size >= 2:
        A, B = unit_rows(DS[idx]), unit_rows(DJ[idx])
        Cm = A @ B.T
        rS, rJ = prim["r_S"][idx], prim["r_J"][idx]
        rng = np.random.default_rng(SEED)
        meds = []
        for _ in range(N_PERM_SHUFFLE):
            p = rng.permutation(idx.size)
            ok = p != np.arange(idx.size)
            v = Cm[np.arange(idx.size)[ok], p[ok]] / np.sqrt(rS[ok] * rJ[p[ok]])
            v = v[np.isfinite(v)]
            if v.size:
                meds.append(float(np.median(v)))
        meds = np.asarray(meds)
        diag["D1_shuffled_pairing"] = dict(
            n_perm=int(meds.size), median=float(np.median(meds)),
            lo=float(np.percentile(meds, 2.5)), hi=float(np.percentile(meds, 97.5)),
            frac_ge_observed=float(np.mean(meds >= prim["M"])) if np.isfinite(prim["M"]) else None)

    def rerun(DS2, DJ2, hAS, hBS, hAJ, hBJ):
        r = readout(DS2, DJ2, hAS, hBS, hAJ, hBJ)
        return {k: r[k] for k in ("n_q", "M", "ci_lo", "ci_hi", "ci_valid", "angle",
                                  "k_specific", "chance", "spec_p")}

    def resid(X, p):
        p = np.asarray(p, float).ravel()
        n = float(np.linalg.norm(p))
        if n == 0:
            return np.asarray(X, float)
        ph = p / n
        X = np.asarray(X, float)
        if X.ndim == 3:
            return X - np.einsum("isk,k->is", X, ph)[..., None] * ph
        return X - np.outer(X @ ph, ph)

    pS, pJ = vs["p_proj"], vj["p_proj"]
    diag["D2_proliferation_removed"] = rerun(
        resid(DS, pS), resid(DJ, pJ), resid(vs["hA"], pS), resid(vs["hB"], pS),
        resid(vj["hA"], pJ), resid(vj["hB"], pJ))
    diag["D2_note"] = ("Both axes are estimated from the same control cells that define each "
                       "screen's control pool, so this diagnostic is not independent of the "
                       "displacements it corrects.")

    def demean(X, m):
        X = np.asarray(X, float)
        return X - m if X.ndim == 2 else X - m[None, None, :]
    mS = DS[idx].mean(0) if idx.size else np.zeros(KDIM)
    mJ = DJ[idx].mean(0) if idx.size else np.zeros(KDIM)
    diag["D3_mean_displacement_removed"] = rerun(
        demean(DS, mS), demean(DJ, mJ), demean(vs["hA"], mS), demean(vs["hB"], mS),
        demean(vj["hA"], mJ), demean(vj["hB"], mJ))

    fg = tab.loc[tab.qualifying, "angle_full_genes_deg"]
    diag["D4_full_gene_angles"] = dict(
        n_G=nG, median_deg=float(fg.median()) if len(fg) else None,
        lo_deg=float(np.percentile(fg, 2.5)) if len(fg) else None,
        hi_deg=float(np.percentile(fg, 97.5)) if len(fg) else None,
        median_cos=float(tab.loc[tab.qualifying, "c_obs_full_genes"].median()) if len(fg) else None)
    diag["D5_norm_shares"] = {
        c: dict(median=float(tab[c].median()), lo=float(np.percentile(tab[c], 2.5)),
                hi=float(np.percentile(tab[c], 97.5)))
        for c in ("share_norm2_kept_by_mask_S", "share_norm2_kept_by_mask_J",
                  "share_norm2_kept_by_proj_S", "share_norm2_kept_by_proj_J")}
    sub = tab[tab.qualifying & np.isfinite(tab.c_corr)]
    if len(sub) >= 3:
        rs = sps.spearmanr(sub.c_corr, sub.n_S)
        rj = sps.spearmanr(sub.c_corr, sub.n_J)
        diag["D6_cell_count"] = dict(n=len(sub), rho_c_corr_vs_n_S=float(rs.statistic),
                                     p_n_S=float(rs.pvalue), rho_c_corr_vs_n_J=float(rj.statistic),
                                     p_n_J=float(rj.pvalue))
    dump("diagnostics.json", diag)

    dump("summary.json", dict(
        reading=reading, gate=gate, gate_passed=gate_passed,
        primary={k: prim[k] for k in ("n_q", "M", "ci_lo", "ci_hi", "ci_valid", "angle",
                                      "k_specific", "chance", "spec_p", "median_c_obs",
                                      "n_nonfinite_ccorr")},
        n_F=int(len(frame)), n_G=nG, min_fq=MIN_FQ,
        n_excluded=int((~prim["qual"]).sum()),
        excluded=[dict(symbol=s, n_S=int(a), n_J=int(b), r_S=None if not np.isfinite(x) else float(x),
                       r_J=None if not np.isfinite(y) else float(y))
                  for s, a, b, x, y in zip(frame.symbol, frame.n_S, frame.n_J,
                                           prim["r_S"], prim["r_J"])
                  if not (np.isfinite(x) and np.isfinite(y) and x >= R_FLOOR and y >= R_FLOOR)]))
    if reading == "UNMEASURABLE" and not gate_passed:
        write_stop("unmeasurable", "the same-context gate did not pass in both screens")
    log.close()


# ==========================================================================
# FINDINGS
# ==========================================================================
def md_table(df, cols, heads=None, d=4):
    heads = heads or cols
    out = ["| " + " | ".join(heads) + " |", "|" + "|".join(["---"] * len(cols)) + "|"]
    for _, r in df.iterrows():
        cells = []
        for c in cols:
            v = r[c]
            if isinstance(v, (bool, np.bool_)):
                cells.append("yes" if v else "no")
            elif isinstance(v, (int, np.integer)):
                cells.append(str(int(v)))
            elif isinstance(v, (float, np.floating)):
                cells.append(fmt(v, d))
            else:
                cells.append(str(v))
        out.append("| " + " | ".join(cells) + " |")
    return out


def findings():
    s = json.loads((OUT / "summary.json").read_text(encoding="utf-8"))
    sets = json.loads((OUT / "stage0_sets.json").read_text(encoding="utf-8"))
    s0 = json.loads((OUT / "stage0_summary.json").read_text(encoding="utf-8"))
    rep = pd.read_csv(OUT / "stage0_reproduction.csv")
    sb = {t: json.loads((OUT / f"{t}_build.json").read_text(encoding="utf-8"))
          for t in ("south", "seng")}
    p = s.get("primary")
    diag = (json.loads((OUT / "diagnostics.json").read_text(encoding="utf-8"))
            if p and (OUT / "diagnostics.json").exists() else None)
    tab = pd.read_csv(OUT / "per_factor.csv") if p and (OUT / "per_factor.csv").exists() else None
    L = []
    A = L.append

    A("# FINDINGS — TRANSFER")
    A("")
    if p is None:
        A(f"Reading: **{s['reading']}**. STOP key `{s['stop_key']}`: {s['stop_message']}.")
        A("")
        A(f"The pre-registered qualification gate failed, so **the primary statistic was not "
          f"computed and no cross-screen number is reported** (flag section 4). |F| = "
          f"{s['n_F']} shared factors, |G| = {s['n_G']} shared genes.")
    else:
        A(f"Reading: **{s['reading']}**. Gate: "
          f"{'passed' if s['gate_passed'] else 'FAILED'}. "
          f"|F| = {s['n_F']}, |F_q| = {p['n_q']}, |G| = {s['n_G']}.")
    A("")
    A("## The pre-registration, verbatim")
    A("")
    A("```")
    A(PREREG.read_text(encoding="utf-8").rstrip())
    A("```")
    A("")
    A("## The amendments, verbatim")
    A("")
    A(AMEND.read_text(encoding="utf-8").rstrip())
    A("")
    A("## Deviations and choices")
    A("")
    for line in DEVIATIONS:
        A(line)
    A("")
    A("## Stage 0 — reproduction of the stored numbers")
    A("")
    A(f"{s0['n_checks']} checks, {s0['n_failed']} failed. Tolerance: counts and the "
      "re-accumulated Southard counts exact; numerics "
      "|new - stored| <= 1e-4 * max(1, |stored|), the tolerance SENG3's Stage 0 used.")
    A("")
    L.extend(md_table(rep, ["block", "item", "kind", "n", "n_mismatch", "max_scaled_diff", "passed"],
                      ["block", "item", "kind", "n", "mismatches", "max scaled diff", "passed"], 6))
    A("")
    A("## The factor set and the gene set")
    A("")
    A(f"- |F| = {s['n_F']} shared factors by exact gene symbol. Southard-only: "
      f"{sets['n_southard_only']}. Sengstack-CRA-only: {sets['n_seng_cra_only']}"
      + (f" ({', '.join(sets['seng_cra_only_names'])})" if sets["n_seng_cra_only"] else "") + ".")
    A(f"- Southard cell counts over F: min {sets['n_S_min']}, median {sets['n_S_median']:.0f}, "
      f"max {sets['n_S_max']}. Sengstack: min {sets['n_J_min']}, median "
      f"{sets['n_J_median']:.0f}, max {sets['n_J_max']}.")
    A(f"- |G| = {sets['n_G']} shared ruler genes, computed at run time from the arrays "
      f"(Southard matched {sets['n_southard_matched']} columns, Sengstack present "
      f"{sets['n_seng_present']} of {sets['n_ruler']} ruler genes). On the 20,265/20,264 "
      "question see the deviation headed \"matched coordinates\" in FINDINGS_COORDSPACE.md: "
      "20,265 ruler genes match the Southard file, two map to one source column, and SOUTH6 "
      "collapses them, so the panel has 20,264 columns.")
    A(f"- G genes inside a screen's all-zero mask: Southard {sets['n_G_in_southard_missing']}, "
      f"Sengstack {sets['n_G_in_seng_missing']}. Those contribute exactly 0 to any difference "
      "either way.")
    A(f"- Frozen basis: results/south7/basis.npz L_mid_full[:, :20], PCA of "
      f"{sets['basis_n_mid']} GTEx donors aged 40-59, max|L'L - I| = "
      f"{sets['basis_orthonormality_err']:.2e}. Read, not refit.")
    A("")
    A("## Stage 4 — the same-context gate, and the result")
    A("")
    A("True drift is zero by construction: the two \"contexts\" are two disjoint halves of the "
      "same cells of the same factor in the same screen, each measured against its own half of "
      "that screen's control pool. The whole readout — observed cosine, split-half "
      "reliabilities from quarter splits, the r floor, the corrected cosine, specificity and "
      "the median with its bootstrap CI — is run on those two contexts. This is also the "
      "recovery check this project requires beside any null: it asks whether the statistic can "
      "see agreement that is certainly there.")
    A("")
    gt = pd.DataFrame([dict(screen=k, **{kk: vv for kk, vv in v.items()}) for k, v in s["gate"].items()])
    L.extend(md_table(gt, ["screen", "n_q", "median_c_obs", "M_same", "ci_lo", "ci_hi",
                           "ci_valid", "k_specific", "chance", "spec_p", "passed"],
                      ["screen", "n_q", "median raw cosine", "M_same", "CI lo", "CI hi",
                       "CI valid", "specific", "chance", "binomial p", "passed"], 4))
    A("")
    A(f"Bars, fixed in advance: **M_same >= {M_SAME_BAR}** and specificity above chance at "
      f"**p < {GATE_P}**, in both screens. Gate "
      f"{'passed' if s['gate_passed'] else '**FAILED**'}"
      + (f" — failed in: {', '.join(s['gate_failed_in'])}." if not s["gate_passed"] else "."))
    A("")
    A("Within-screen context for the gate, from the gate's own tables "
      "(`gate_south.csv`, `gate_seng.csv`) and the build logs:")
    A("")
    for t, label in (("south", "Southard"), ("seng", "Sengstack CRA")):
        d = pd.read_csv(OUT / f"gate_{t}.csv")
        A(f"- **{label}**: control pool {sb[t]['ctrl_cells']} cells; factor cells over F "
          f"min {int(d.n_cells.min())}, median {d.n_cells.median():.0f}, max "
          f"{int(d.n_cells.max())}. Median split-half reliability of context 1 "
          f"{fmt(d.r_1.median())} and of context 2 {fmt(d.r_2.median())}; "
          f"{int((d.r_1 < R_FLOOR).sum() + d.r_1.isna().sum())} factors fell below the "
          f"r >= {R_FLOOR} floor on context 1, {int((d.r_2 < R_FLOOR).sum() + d.r_2.isna().sum())} "
          f"on context 2; {int(d.qualifying.sum())} of {len(d)} qualified.")
    A("")
    A("What the gate says, and what it does not. Specificity is comfortably above chance in "
      "both screens: within a screen, a factor's half-1 displacement is closer to its own "
      "half-2 displacement than to other factors' — the readout can tell factors apart. What "
      "fails is **magnitude recovery in the Sengstack arm**: on data whose true agreement is 1 "
      "by construction, the reliability-corrected cosine recovers only "
      f"{fmt(s['gate']['seng']['M_same'])} (CI {fmt(s['gate']['seng']['ci_lo'])} to "
      f"{fmt(s['gate']['seng']['ci_hi'])}), against {fmt(s['gate']['south']['M_same'])} in "
      "Southard. A corrected cosine that cannot reach 1 where the answer is 1 cannot be read "
      "as a measure of how far short of 1 two screens fall.")
    A("")
    A("A likely mechanism, stated as an observation about the instrument and **not** as a "
      f"tested claim: the Sengstack CRA control pool is {sb['seng']['ctrl_cells']} cells, "
      "against a median of "
      f"{pd.read_csv(OUT / 'gate_seng.csv').n_cells.median():.0f} cells per factor. The "
      "correction assumes the noise a split can see is the noise that limits the displacement, "
      "but control-side noise is shared between the two halves of a split and is not measured "
      "by it, and at this control depth that noise is not negligible. Southard's control pool "
      f"is {sb['south']['ctrl_cells']} cells and its gate recovers "
      f"{fmt(s['gate']['south']['M_same'])}. Testing that mechanism would need a new "
      "pre-registration; nothing here does it.")
    A("")
    A("## The pre-registered reading that fired")
    A("")
    A(f"**{s['reading']}**")
    A("")
    A(reading_text(s, diag))
    A("")
    if p is None:
        A("Per flag section 4 the primary statistic was not computed and no cross-screen number "
          "is reported. The pre-registered diagnostics D1-D6 are cross-screen and are not "
          "reported either. What was computed before the gate's outcome was converted into a "
          "stop is quarantined under `results/transfer/NOT_A_RESULT_gate_failed/`, with its own "
          "README; see deviation 1.")
        A("")
    else:
        A("## Primary readout")
        A("")
        nq, chance = p["n_q"], p["chance"]
        A(f"- |F_q| = {nq} of {s['n_F']} factors had r_S >= {R_FLOOR} and r_J >= {R_FLOOR}; "
          f"{s['n_excluded']} were excluded.")
        A(f"- M = median corrected cosine over F_q = **{fmt(p['M'])}** "
          f"(95% bootstrap CI over factors {fmt(p['ci_lo'])} to {fmt(p['ci_hi'])}, "
          f"{'valid' if p['ci_valid'] else 'INVALID'}), angle {fmt(p['angle'], 1)} degrees.")
        A(f"- Median uncorrected cosine over F_q: {fmt(p['median_c_obs'])}.")
        A(f"- Specific factors: **{p['k_specific']} of {nq}**, against a chance baseline of "
          f"{fmt(chance, 2)} by construction (one-sided exact binomial p = "
          f"{fmt(p['spec_p'], 6)}; approximate, because the mismatched cosines are not "
          "independent across factors).")
        A(f"- Non-finite corrected cosines inside F_q: {p['n_nonfinite_ccorr']}.")
        A("")
        A("## Diagnostics (not gated results)")
        A("")
        for line in diag_lines(diag, s):
            A(line)
        A("")
        A("## Per-factor table")
        A("")
        cols = ["symbol", "n_S", "n_J", "r_S", "r_J", "c_obs", "angle_deg", "c_corr",
                "qualifying", "specific", "mismatch_p95", "rank_by_c_corr", "in_active_212"]
        heads = ["factor", "n_S", "n_J", "r_S", "r_J", "c_obs", "angle", "c_corr",
                 "qualifying", "specific", "mismatch 95th", "rank", "in SOUTH7 active set"]
        L.extend(md_table(tab, cols, heads, 4))
        A("")
        A("Full columns, including the norm bookkeeping and the full-gene angles, are in "
          "`results/transfer/per_factor.csv`.")
        A("")
    A("## What this does not cover")
    A("")
    for line in NOT_COVERED:
        A(line)
    FINDINGS.write_text("\n".join(L) + "\n", encoding="utf-8")
    print(f"wrote {FINDINGS}")


def reading_text(s, diag):
    p = s.get("primary")
    if p is None:
        return ("\"The method could not see transfer. The same-context control did not recover "
                "zero-drift agreement in both screens. No statement is made about whether "
                "factor directions transfer between the Southard and Sengstack screens. This "
                "is a statement about the measurement, not about the biology.\"")
    nq, chance = p["n_q"], p["chance"]
    M, lo, hi = fmt(p["M"]), fmt(p["ci_lo"]), fmt(p["ci_hi"])
    ang = fmt(p["angle"], 1)
    k, pv = p["k_specific"], fmt(p["spec_p"], 6)
    if s["reading"] == "UNMEASURABLE":
        why = ("The same-context control did not recover zero-drift agreement in both screens"
               if not s["gate_passed"] else
               f"Only {nq} factors had sufficient within-screen reliability, fewer than the 20 "
               "fixed in advance")
        return ("\"The method could not see transfer. " + why + ". No statement is made about "
                "whether factor directions transfer between the Southard and Sengstack screens. "
                "This is a statement about the measurement, not about the biology.\"")
    if s["reading"] == "LARGE DRIFT":
        return ("\"Factor directions measured in Southard do not carry over to another "
                "fibroblast screen even for the factors most likely to transfer. The median "
                f"reliability-corrected cosine between the same factor's displacement in the two "
                f"screens is {M} (95% CI {lo} to {hi}, {ang} degrees) over {nq} qualifying "
                f"factors, and {k} of {nq} factors were specific against a chance baseline of "
                f"{fmt(chance, 2)} (binomial p = {pv}). A forward model cannot be initialised "
                "from Southard geometry. The combination screen must measure single-factor "
                "effects in its own cells.\n\nThe two screens differ in lab, in guide library "
                "and in platform as well as in cell state, so the drift measured here is not "
                "attributed to cell state alone. That does not weaken the kill, because any "
                "future screen differs from Southard in the same ways.\"")
    if s["reading"] == "SMALL DRIFT":
        return (f"\"For these {nq} factors, displacement directions measured in Southard point "
                f"the same way as in Sengstack's CRA arm: median reliability-corrected cosine "
                f"{M} (95% CI {lo} to {hi}, {ang} degrees), with {k} of {nq} factors specific "
                f"against a chance baseline of {fmt(chance, 2)} (binomial p = {pv}).\n\nThis "
                f"holds only for these {nq} factors, on the most favourable comparison "
                "available: the shared factors were chosen by Sengstack from passage- and "
                "proliferation-linked modules and are not a random sample of Southard's 1,836 "
                "targets. It says nothing about the untested factors. It does NOT justify "
                "initialising a forward model from Southard geometry.\n\nA SMALL DRIFT result "
                "must not be cited as support for initialising the forward model from "
                "Southard.\"")
    return ("\"The measurement is between the two pre-registered thresholds: median "
            f"reliability-corrected cosine {M} (95% CI {lo} to {hi}, {ang} degrees) over {nq} "
            f"qualifying factors, with {k} of {nq} specific against a chance baseline of "
            f"{fmt(chance, 2)} (binomial p = {pv}). This is not evidence for initialising a "
            "forward model from Southard. A result in this band leaves the kill-check "
            "unresolved and does not license using Southard geometry as a starting point.\"")


def diag_lines(diag, s):
    out = []
    d1 = diag.get("D1_shuffled_pairing")
    if d1:
        out.append(f"- **D1, shuffled-identity pairing** (diagnostic, not a gated result). "
                   f"{d1['n_perm']} permutations of the Sengstack pairing: median of the median "
                   f"corrected cosine {fmt(d1['median'])} (2.5-97.5% {fmt(d1['lo'])} to "
                   f"{fmt(d1['hi'])}). Share of permutations at or above the observed M: "
                   f"{fmt(d1['frac_ge_observed'])}. This is the \"nothing factor-specific "
                   f"transfers\" reference.")
    for key, label in (("D2_proliferation_removed",
                        "D2, each screen's own frozen proliferation axis projected out"),
                       ("D3_mean_displacement_removed",
                        "D3, each screen's mean displacement over F_q removed")):
        d = diag.get(key)
        if d:
            out.append(f"- **{label}** (diagnostic, not a gated result). n_q {d['n_q']}, "
                       f"M {fmt(d['M'])} (CI {fmt(d['ci_lo'])} to {fmt(d['ci_hi'])}), "
                       f"{d['k_specific']} specific against chance {fmt(d['chance'], 2)}, "
                       f"binomial p {fmt(d['spec_p'], 6)}.")
    if diag.get("D2_note"):
        out.append(f"  - {diag['D2_note']}")
    d4 = diag.get("D4_full_gene_angles")
    if d4 and d4.get("median_deg") is not None:
        out.append(f"- **D4, full-gene-space angles on G** (diagnostic, not a gated result). "
                   f"Median {fmt(d4['median_deg'], 1)} degrees (2.5-97.5% "
                   f"{fmt(d4['lo_deg'], 1)} to {fmt(d4['hi_deg'], 1)}), median cosine "
                   f"{fmt(d4['median_cos'])}, over {d4['n_G']} shared genes. Expected near 90 "
                   "degrees: in a space of that dimension two displacement vectors that share "
                   "no measurement noise are nearly orthogonal whatever their shared signal, "
                   "which is the property that made the full-gene space blind in SENG2 and "
                   "SOUTH6REC.")
    d5 = diag.get("D5_norm_shares")
    if d5:
        out.append("- **D5, norm bookkeeping** (diagnostic, not a gated result), as shares of "
                   "squared norm, median [2.5-97.5%]:")
        nice = {"share_norm2_kept_by_mask_S": "Southard, kept by the shared-gene mask",
                "share_norm2_kept_by_mask_J": "Sengstack, kept by the shared-gene mask",
                "share_norm2_kept_by_proj_S": "Southard, kept by the 20-component projection",
                "share_norm2_kept_by_proj_J": "Sengstack, kept by the 20-component projection"}
        for k2, v in d5.items():
            out.append(f"  - {nice[k2]}: {fmt(v['median'])} [{fmt(v['lo'])}, {fmt(v['hi'])}]")
    d6 = diag.get("D6_cell_count")
    if d6:
        out.append(f"- **D6, cell-count check** (diagnostic, not a gated result). Spearman of "
                   f"c_corr with n_S: rho {fmt(d6['rho_c_corr_vs_n_S'])} (p "
                   f"{fmt(d6['p_n_S'], 4)}); with n_J: rho {fmt(d6['rho_c_corr_vs_n_J'])} "
                   f"(p {fmt(d6['p_n_J'], 4)}), over {d6['n']} factors. Group size in these "
                   "screens is an accident of the experiment, so this check is reported "
                   "whatever it shows.")
    return out


DEVIATIONS = [
    "Everything here is a deviation from the flag, an addition to it, or a choice the flag did "
    "not fix. Housekeeping is included.",
    "",
    "1. **The flag's stop ordering was violated on the first pass, and this is the most "
    "serious deviation in this task.** Flag section 4 says that if the gate fails, the primary "
    "statistic is not computed and no cross-screen number is reported. The gate failed, but the "
    "first pass of `src/transfer_run.py read` evaluated the gate and then ran the primary "
    "readout and the cross-screen diagnostics unconditionally, so those numbers were computed "
    "and the session that ran them has seen them. They are quarantined, not deleted, under "
    "`results/transfer/NOT_A_RESULT_gate_failed/` with a README stating what they are; no "
    "number in them is quoted anywhere in this findings file and no reading is drawn from them. "
    "`src/transfer_run.py` was then corrected so that the `read` stage stops at a failed gate, "
    "and the corrected stage was re-run: the gate is deterministic at seed 20260914 and its "
    "outcome, its numbers and the stop key are unchanged. That re-run changed no threshold and "
    "is not a resumption of a stopped run. A reader who wants the gate re-verified from tracked "
    "files alone can run `stage0`, `south`, `seng`, `read` in order and will get this same stop.",
    "2. **Stage 0(b) narrowed** (amendment A2). The SENG3 floors and the floor-dependent "
    "boolean columns are NOT recomputed. Stage 0(b) reproduces the CRA count, every CRA "
    "`n_cells`, and `delta_i`, `frac_i`, `delta_ii`, `frac_ii` for the CRA rows of "
    "`results/seng3/d4_table.csv`. Stage 0 therefore does not re-verify SENG3's floors; "
    "nothing in TRANSFER depends on them.",
    "3. **The Southard h5ad is never opened** (amendment, implementation definition 4). Counts "
    "come from the SOUTH6 matched CSR cache, so its md5 is not recomputed and the "
    "`file_mismatch` stop key can fire only on the Sengstack file, whose md5 is verified.",
    "4. **Stage 0(a) does not re-accumulate the 78 per-guide NT rows.** Rebuilding `C_nt` per "
    "guide needs the per-cell guide labels, which live in the Southard h5ad obs, excluded by "
    "deviation 3. `C_f` (every factor, from its own cells) and `C_pool` (the pooled control, "
    "from `nt_cells`) are re-accumulated and compared exactly; `C_pool` is also checked against "
    "the sum of the stored `C_nt` rows; and the stored `C_nt` rows are used as panel input when "
    "rebuilding the TMM transform. This is weaker than a full panel rebuild and is reported as "
    "such.",
    "5. **The half-split displacement is factor-half minus an independent control-half.** This "
    "differs from SENG3 amendment B5's within-set paired construction "
    "(`z(draw) - z(pool minus draw)`), which is not the construction `d_S`/`d_J` uses. The "
    "complement half's pseudobulk is obtained by subtracting the drawn half's counts from the "
    "group total, which is exact in float64 for integer counts.",
    "6. **Seed 20260914 is used for the factor bootstrap**, per the flag, departing from this "
    "project's usual 20260918 bootstrap seed.",
    "7. **Control split streams** are keyed `default_rng([20260914, screen_id, split_index])`, "
    "a three-element key that cannot collide with the flag's four-element per-factor keys. The "
    "gate's fixed contexts use key element 7 and its quarter splits key element 8.",
    "8. **Per-stage caching.** `stage0` writes the Southard panel TMM cache to "
    "`results/transfer/south_cache.npz` and the later stages read it instead of rebuilding the "
    "panel. Each stage's random streams run uninterrupted from their seeds, so this is not "
    "checkpoint/resume of a draw stream.",
    "9. **No PROGRESS file**, per the flag; the stop status is carried here.",
    "10. **The Southard factor universe is all 1,836 ranked targets**, not SOUTH7's published "
    "304-gene primary set; membership in that set is carried as a descriptive column only.",
    "11. **No depth-matched control variant** (amendment A4) and **no feasibility probe** "
    "(amendment A5).",
]

NOT_COVERED = [
    "- **The kill-check is unresolved.** UNMEASURABLE is not a weak version of either outcome: "
    "LARGE DRIFT did not fire, and neither did SMALL DRIFT. This task establishes neither that "
    "Southard's factor directions fail to transfer nor that they transfer. Step 1 of the "
    "blueprint is not answered, and nothing here licenses initialising a forward model from "
    "Southard geometry or declines to. A further attempt needs a statistic whose zero-drift "
    "recovery passes in both screens, pre-registered separately, before any cross-screen number "
    "is computed.",
    "- **Neither screen is aged donor cells.** Southard is Hs27, a neonatal foreskin fibroblast "
    "line. Sengstack's \"old\" is passage number, not donor age, and this project has already "
    "found that passage does not resemble donor aging (cos(v,u) = -0.040, permutation p = 1.0, "
    "`results/seng/axes.json`). Nothing here transfers to aged human cells, in either direction.",
    "- Transfer is judged **only inside the 20-component space** the search objective is scored "
    "in. A factor pair could agree in some other space and disagree here, or the reverse.",
    "- The 20-component basis was selected, in SOUTH7 and SOUTH9, because it was the space in "
    "which the **Southard** screen's statistic recovered planted effects. The space is chosen on "
    "Southard's behalf, which makes it more favourable to Southard than to Sengstack. Disclosed, "
    "not corrected.",
    "- The shared factors are **not a random sample** of Southard's 1,836 targets. Sengstack "
    "selected their perturbation set from passage- and proliferation-linked modules.",
    "- The specificity test treats factors as independent. They are not: every mismatched cosine "
    "is drawn from the same two displacement pools, so the binomial p is approximate.",
    "- The reliability correction measures only the noise a split can see. Control-side noise is "
    "shared between the two halves and is not measured or corrected (amendment A4).",
    "- This task measures **direction agreement only**. It says nothing about effect magnitude, "
    "about whether either displacement is in a useful direction, or about whether any factor "
    "approaches young in either screen.",
]


STAGES = dict(stage0=stage0, south=lambda: build("south"), seng=lambda: build("seng"),
              read=read, findings=findings)

if __name__ == "__main__":
    require_prereg()
    stage = sys.argv[1]
    try:
        STAGES[stage]()
    except Stop as e:
        write_stop(e.key, e.message)
        print(f"STOP {e.key}: {e.message}", file=sys.stderr)
        raise SystemExit(3)
