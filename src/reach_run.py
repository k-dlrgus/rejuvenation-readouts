"""REACH — the reachability bound. How close to the young anchor can any
weighted sum of Southard's measured factor displacements get?

Pre-registration: results/reach/PREREG.flag (2026-10-17T22:22Z), written
before any REACH number. This script reads it and never rewrites it. It is
"sizing information, never evidence" (TF search model blueprint).

Instruments are read, never refit: SOUTH6 d_f and counts (results/south6/
stage1.npz and the matched-gene CSR cache), the S2A_k20 basis and the O1/Y1
anchors (results/south7/basis.npz), the frozen GTEx mu/sd. src/south6_run.py,
src/south6rec_run.py, src/south7_run.py, src/seng_run.py and src/toward_run.py
are imported unchanged. Nothing is written outside results/reach/ and
FINDINGS_REACH.md.

Inputs: results/south6/stage1.npz, matched_data.f64, matched_indices.i32 and
matched_indptr.npy are git-ignored; they are regenerated from the Southard
download by src/south6_run.py (REPO_INVENTORY.md). The readings are
re-verifiable from tracked files alone (results/reach/inputs_k20.npz):
    python src/reach_run.py verify

Usage, from the repository root:
    python src/reach_run.py           Stage 0, draws, sets A-D, readings,
                                      diagnostics, findings
    python src/reach_run.py verify    re-solve sets A-D from inputs_k20.npz and
                                      re-derive the reading
    python src/reach_run.py findings  rewrite FINDINGS_REACH.md from outputs
"""
from __future__ import annotations

import json
import os
import sys
import time
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.optimize import lsq_linear, minimize

SRC = Path(__file__).resolve().parent
if str(SRC) not in sys.path:
    sys.path.insert(0, str(SRC))

ROOT = SRC.parent
OUT = ROOT / "results" / "reach"
PREREG = OUT / "PREREG.flag"
FINDINGS = ROOT / "FINDINGS_REACH.md"
SOUTH6 = ROOT / "results" / "south6"
SOUTH7 = ROOT / "results" / "south7"

SEED = 20260914
KDIM = 20
KS = (1, 3, 5, 10, 50)
K_ALL = ("1", "3", "5", "10", "50", "unlimited")
N_TARGETS = 200
N_NT = 14013
N_FACTORS = 1836
N_RULER = 23485
TOL = 1e-4            # |new - stored| <= TOL * max(1, |stored|)
ORTHO_TOL = 1e-8
GAP_REL = 1e-8        # duality gap <= GAP_REL * max(1, 0.5 ||b||^2)
NNZ_EPS = 1e-6
TIE = 1e-6
N_BISECT = 60
LAM_LO = 1e-10
N_RESTART = 20
N_SAMPLE = 20
EIG_REL = 1e-10
N_WORKERS = max(1, min(16, (os.cpu_count() or 2) - 2))

SET_LABEL = {
    "A": "A raw",
    "B": "B shrunk",
    "C": "C null bank",
    "Cs": "C null bank, shrunk (diagnostic)",
}


class Stop(Exception):
    def __init__(self, key, message):
        super().__init__(message)
        self.key = key
        self.message = message


class Log:
    def __init__(self, path):
        self.fh = open(path, "a", encoding="utf-8")
        self.t0 = time.time()

    def __call__(self, msg):
        line = f"[{time.time() - self.t0:7.0f}s] {msg}"
        print(line, flush=True)
        self.fh.write(line + "\n")
        self.fh.flush()


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


def dump(name, obj):
    (OUT / name).write_text(json.dumps(_js(obj), indent=2), encoding="utf-8")


def fmt(x, d=3):
    if x is None:
        return "n/a"
    x = float(x)
    if not np.isfinite(x):
        return "n/a"
    if x != 0.0 and abs(x) < 10 ** (-d):
        return f"{x:.2e}"
    return f"{x:.{d}f}"


# --------------------------------------------------------------------------
# solvers (module level and import-light, so pool workers can use them)
# --------------------------------------------------------------------------
class Problem:
    """min ||A w - b||^2 + c0 over w in [0, 1]^n; frac = (g - sqrt(.)) / g.

    c0 is the squared residual no column can touch (0 in S2A_k20; the part of
    b outside the column span in the Gram-reduced full-space form)."""

    def __init__(self, A, b, g, c0=0.0, tol_scale=None):
        self.A = np.ascontiguousarray(A, np.float64)
        self.b = np.asarray(b, np.float64).ravel()
        self.g = float(g)
        self.c0 = float(c0)
        half = 0.5 * (self.g ** 2) if tol_scale is None else float(tol_scale)
        self.tol = GAP_REL * max(1.0, half)

    def frac(self, w):
        e = self.A @ w - self.b
        return (self.g - np.sqrt(max(float(e @ e) + self.c0, 0.0))) / self.g


def gap_of(A, b, w, lam):
    y = A @ w - b
    P = 0.5 * float(y @ y) + lam * float(w.sum())
    D = -0.5 * float(y @ y) - float(b @ y) + float(np.minimum(0.0, A.T @ y + lam).sum())
    return P - D


def bvls(A, b, tol):
    n = A.shape[1]
    if n == 0:
        return np.zeros(0), True, dict(status=1, gap=0.0)
    r = lsq_linear(A, b, bounds=(0.0, 1.0), method="bvls", tol=1e-12, max_iter=10000)
    w = np.clip(r.x, 0.0, 1.0)
    gp = gap_of(A, b, w, 0.0)
    ok = bool(r.status >= 1 and gp <= tol)
    return w, ok, dict(status=int(r.status), gap=float(gp))


def _clean(w):
    w = np.asarray(w, np.float64).copy()
    w[w <= NNZ_EPS] = 0.0
    return w


def solve_unlimited(P):
    w, ok, info = bvls(P.A, P.b, P.tol)
    w = _clean(w)
    return dict(w=w, nnz=int((w > 0).sum()), frac=P.frac(w), ok=ok, gap=info["gap"])


def greedy(P, kmax=max(KS)):
    A, b = P.A, P.b
    n = A.shape[1]
    colsq = np.einsum("ij,ij->j", A, A)
    usable = colsq > 0
    in_S = np.zeros(n, bool)
    S = []
    w = np.zeros(n)
    e = -b.copy()
    best = {k: dict(w=np.zeros(n), nnz=0, frac=0.0) for k in KS}
    ok_all, max_gap, steps = True, 0.0, 0
    while True:
        dots = A.T @ e
        with np.errstate(divide="ignore", invalid="ignore"):
            t = np.clip(np.where(usable, -dots / np.where(usable, colsq, 1.0), 0.0), 0.0, 1.0)
        gain = -(2.0 * t * dots + t * t * colsq)
        gain[in_S | ~usable] = -np.inf
        j = int(np.argmax(gain))
        if not np.isfinite(gain[j]) or gain[j] <= 0:
            break
        S.append(j)
        in_S[j] = True
        wS, ok, info = bvls(A[:, S], b, P.tol)
        ok_all &= ok
        max_gap = max(max_gap, info["gap"])
        w = np.zeros(n)
        w[S] = wS
        w = _clean(w)
        e = A @ w - b
        nnz = int((w > 0).sum())
        fr = P.frac(w)
        steps += 1
        for k in KS:
            if nnz <= k and fr > best[k]["frac"]:
                best[k] = dict(w=w.copy(), nnz=nnz, frac=fr)
        if nnz >= kmax or len(S) >= n or float(e @ e) <= 1e-24:
            break
    return best, ok_all, dict(max_gap=max_gap, steps=steps)


def l1_solve(A, b, lam, w0, tol):
    n = A.shape[1]
    bounds = [(0.0, 1.0)] * n

    def f(w):
        e = A @ w - b
        return 0.5 * float(e @ e) + lam * float(w.sum()), A.T @ e + lam

    w = w0.copy()
    gp = np.inf
    for _ in range(N_RESTART + 1):
        res = minimize(f, w, jac=True, method="L-BFGS-B", bounds=bounds,
                       options=dict(maxiter=20000, ftol=1e-16, gtol=1e-12, maxcor=30))
        w = np.clip(res.x, 0.0, 1.0)
        gp = gap_of(A, b, w, lam)
        if gp <= tol:
            return w, True, gp
    return w, False, gp


def l1(P):
    A, b = P.A, P.b
    n = A.shape[1]
    c = A.T @ b
    lam_max = float(c.max())
    out = {}
    if lam_max <= 0:
        for k in KS:
            out[k] = dict(w=np.zeros(n), nnz=0, frac=0.0, nnz_hit=0, hit=False,
                          lam=np.nan, ok=True)
        return out, True, dict(max_gap=0.0, n_solves=0)
    solved = {}
    ok_all, max_gap, n_solves = True, 0.0, 0
    lo, hi = np.log(lam_max * LAM_LO), np.log(lam_max)

    def solve_at(t):
        nonlocal ok_all, max_gap, n_solves
        if t in solved:
            return solved[t]
        if solved:
            near = min(solved, key=lambda s: abs(s - t))
            w0 = solved[near]
        else:
            w0 = np.zeros(n)
        w, ok, gp = l1_solve(A, b, float(np.exp(t)), w0, P.tol)
        ok_all &= ok
        max_gap = max(max_gap, gp)
        n_solves += 1
        solved[t] = w
        return w

    for k in KS:
        a, z = lo, hi
        cand = None
        for _ in range(N_BISECT):
            t = 0.5 * (a + z)
            w = solve_at(t)
            nnz = int((w > NNZ_EPS).sum())
            if nnz <= k and (cand is None or nnz > cand[2] or (nnz == cand[2] and t < cand[0])):
                cand = (t, w, nnz)
            if nnz == k:
                break
            if nnz > k:
                a = t
            else:
                z = t
        if cand is None:
            out[k] = dict(w=np.zeros(n), nnz=0, frac=0.0, nnz_hit=0, hit=False, lam=np.nan, ok=True)
            continue
        t, wl, nnz_hit = cand
        supp = np.flatnonzero(wl > NNZ_EPS)
        wS, ok, info = bvls(A[:, supp], b, P.tol)
        ok_all &= ok
        max_gap = max(max_gap, info["gap"])
        w = np.zeros(n)
        w[supp] = wS
        w = _clean(w)
        out[k] = dict(w=w, nnz=int((w > 0).sum()), frac=P.frac(w), nnz_hit=nnz_hit,
                      hit=bool(nnz_hit == k), lam=float(np.exp(t)), ok=ok)
    return out, ok_all, dict(max_gap=max_gap, n_solves=n_solves)


def jaccard(w1, w2):
    a, b = set(np.flatnonzero(w1 > 0)), set(np.flatnonzero(w2 > 0))
    if not a and not b:
        return 1.0
    return len(a & b) / len(a | b)


def run_set(P):
    """Every k for one problem. Returns per-k dicts and an overall ok flag."""
    gr, ok_g, ig = greedy(P)
    lr, ok_l, il = l1(P)
    un = solve_unlimited(P)
    res = {}
    for k in KS:
        G, L = gr[k], lr[k]
        res[str(k)] = dict(
            greedy=G, l1=L, gated=max(G["frac"], L["frac"]),
            abs_diff=abs(G["frac"] - L["frac"]), jaccard=jaccard(G["w"], L["w"]),
        )
    res["unlimited"] = dict(bvls=un, gated=un["frac"])
    ok = bool(ok_g and ok_l and un["ok"])
    info = dict(greedy=ig, l1=il, unlimited_gap=un["gap"], ok_greedy=ok_g, ok_l1=ok_l,
                ok_bvls=un["ok"])
    return res, ok, info


_W = {}


def _worker_init(A, g):
    _W["A"] = np.asarray(A, np.float64)
    _W["g"] = float(g)


def _worker_target(args):
    t, b = args
    P = Problem(_W["A"], b, _W["g"])
    res, ok, info = run_set(P)
    row = dict(target=int(t), ok=ok)
    for k in KS:
        r = res[str(k)]
        row[f"greedy_k{k}"] = r["greedy"]["frac"]
        row[f"l1_k{k}"] = r["l1"]["frac"]
        row[f"l1_hit_k{k}"] = r["l1"]["hit"]
        row[f"gated_k{k}"] = r["gated"]
    row["gated_kunlimited"] = res["unlimited"]["gated"]
    row["nnz_unlimited"] = res["unlimited"]["bvls"]["nnz"]
    row["max_gap"] = max(info["greedy"]["max_gap"], info["l1"]["max_gap"], info["unlimited_gap"])
    return row


def run_targets(A, g, targets, log):
    """Set (D): (B)'s vectors against each random target, in a process pool."""
    import multiprocessing as mp
    for v in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS"):
        os.environ[v] = "1"
    rows = []
    t0 = time.time()
    with mp.get_context("spawn").Pool(N_WORKERS, initializer=_worker_init, initargs=(A, g)) as pool:
        for row in pool.imap(_worker_target, list(enumerate(targets))):
            rows.append(row)
            if not row["ok"]:
                pool.terminate()
                raise Stop("solver_not_converged",
                           f"set D target {row['target']}: a solve failed the duality-gap test "
                           f"(max gap {row['max_gap']:.3e}).")
            if len(rows) % 20 == 0:
                log(f"[D] {len(rows)}/{len(targets)} targets ({time.time() - t0:.0f}s)")
    return pd.DataFrame(rows)


# --------------------------------------------------------------------------
# Stage 0
# --------------------------------------------------------------------------
def check(rows, item, stored, got, kind="tol", gated=True, extra=False, note=""):
    s = np.asarray(stored).ravel()
    r = np.asarray(got).ravel()
    if s.shape != r.shape:
        rows.append(dict(item=item, kind=kind, gated=gated, extra=extra, n=int(max(s.size, r.size)),
                         stored=str(s[:1]), reproduced=str(r[:1]), max_abs_diff=np.nan,
                         max_scaled_diff=np.nan, n_mismatch=int(max(s.size, r.size)), ok=False,
                         note=f"shape stored {s.shape} reproduced {r.shape}"))
        return False
    if kind == "exact":
        bad = int(np.sum(s != r))
        mad = msd = np.nan
    else:
        s = s.astype(np.float64)
        r = r.astype(np.float64)
        both_nan = np.isnan(s) & np.isnan(r)
        diff = np.where(both_nan, 0.0, np.abs(s - r))
        sc = diff / np.maximum(1.0, np.abs(np.nan_to_num(s)))
        bad = int(np.sum(~(sc <= TOL)))
        mad = float(np.nanmax(diff)) if diff.size else 0.0
        msd = float(np.nanmax(sc)) if sc.size else 0.0
    one = s.size == 1
    rows.append(dict(item=item, kind=kind, gated=gated, extra=extra, n=int(s.size),
                     stored=(s[0].item() if one and hasattr(s[0], "item") else (s[0] if one else "")),
                     reproduced=(r[0].item() if one and hasattr(r[0], "item") else (r[0] if one else "")),
                     max_abs_diff=mad, max_scaled_diff=msd, n_mismatch=bad, ok=bool(bad == 0),
                     note=note))
    return bad == 0


def stage0_static(log, rows):
    """(a)-(e): counts, basis, distances, per-factor norms. Cheap."""
    s1 = np.load(SOUTH6 / "stage1.npz", allow_pickle=True)
    d_f = np.asarray(s1["d_f"], np.float64)
    keep = [str(x) for x in s1["keep"]]
    cell_indptr = np.asarray(s1["cell_indptr"], np.int64)
    n_cells = np.diff(cell_indptr)
    log(f"[stage1] d_f {d_f.shape} keep {len(keep)}")

    cpt = pd.read_csv(SOUTH6 / "stage0_cells_per_target.csv")
    cpt_f = cpt[cpt.target != "non"]
    check(rows, "(a) factor count: d_f rows", N_FACTORS, d_f.shape[0], "exact")
    check(rows, "(a) factor count: len(keep)", N_FACTORS, len(keep), "exact")
    check(rows, "(a) factor count: stage0_cells_per_target.csv factor rows", N_FACTORS, len(cpt_f), "exact")
    check(rows, "(a) d_f columns (ruler coordinates)", N_RULER, d_f.shape[1], "exact")
    stored_n = cpt_f.set_index("target").reindex(keep)["n_cells"].to_numpy()
    check(rows, "(a) cells per factor equal stage0_cells_per_target.csv", stored_n, n_cells, "exact")

    bz = np.load(SOUTH7 / "basis.npz", allow_pickle=True)
    L = np.asarray(bz["L_mid_full"], np.float64)[:, :KDIM]
    mr = np.asarray(bz["matched_ruler"], np.int64)
    ortho = float(np.max(np.abs(L.T @ L - np.eye(KDIM))))
    rows.append(dict(item="(b) basis orthonormality max|L20^T L20 - I| <= 1e-8", kind="threshold",
                     gated=True, extra=False, n=KDIM * KDIM, stored=ORTHO_TOL, reproduced=ortho,
                     max_abs_diff=ortho, max_scaled_diff=np.nan, n_mismatch=int(ortho > ORTHO_TOL),
                     ok=bool(ortho <= ORTHO_TOL), note="threshold, not a stored value"))
    log(f"[basis] L20 {L.shape} max|LtL-I| = {ortho:.3e}")
    c_old = np.asarray(bz["c_old"], np.float64).ravel()
    c_young = np.asarray(bz["c_young"], np.float64).ravel()
    r_full = c_old - c_young
    r20 = r_full @ L
    otd = pd.read_csv(SOUTH6 / "origin_target_distances.csv")
    st_s1 = otd[(otd.space == "S1") & (otd.origin == "O1") & (otd.target == "Y1")]["distance"].to_numpy()
    check(rows, "(c) ||O1 - Y1|| S1 vs south6/origin_target_distances.csv", st_s1, np.linalg.norm(r_full))
    ret = pd.read_csv(SOUTH7 / "basis_retention.csv")
    st_20 = ret[ret.k == 20]["origin_target_distance"].to_numpy()
    check(rows, "(c) ||O1 - Y1|| S2A_k20 vs south7/basis_retention.csv (deviation: not in the south6 file)",
          st_20, np.linalg.norm(r20))
    spaces = [str(x) for x in bz["spaces"]]
    check(rows, "(c) ||O1 - Y1|| S2A_k20 vs south7/basis.npz rnorm", bz["rnorm"][spaces.index("S2A_k20")],
          np.linalg.norm(r20), extra=True)
    check(rows, "(c) r_S2A_k20 vector vs south7/basis.npz", bz["r_S2A_k20"], r20, extra=True)
    check(rows, "(c) r_S1 on matched genes vs south7/basis.npz", bz["r_S1_matched"], r_full[mr], extra=True)

    D20 = d_f @ L
    dn_s1 = np.linalg.norm(d_f, axis=1)
    dn_20 = np.linalg.norm(D20, axis=1)
    tab = pd.read_csv(SOUTH6 / "stage3_full_table.csv", usecols=["factor", "space", "pair", "d_norm"])
    t6 = tab[(tab.space == "S1") & (tab.pair == "O1_Y1")].set_index("factor").reindex(keep)["d_norm"].to_numpy()
    obs7 = pd.read_csv(SOUTH7 / "stage_a1_observed.csv", usecols=["factor", "dnorm_S2A_k20"])
    t7 = obs7.set_index("factor").reindex(keep)["dnorm_S2A_k20"].to_numpy()
    pick = np.sort(np.random.default_rng([SEED, 0]).choice(N_FACTORS, N_SAMPLE, replace=False))
    check(rows, "(d) ||d_f|| S1, 20 sampled factors, vs south6/stage3_full_table.csv d_norm",
          t6[pick], dn_s1[pick], note="factors " + ",".join(keep[i] for i in pick))
    check(rows, "(d) ||d_f L20|| S2A_k20, 20 sampled factors, vs south7/stage_a1_observed.csv",
          t7[pick], dn_20[pick])
    check(rows, "(e) ||d_f|| S1, all 1836 factors, vs south6/stage3_full_table.csv", t6, dn_s1, extra=True)
    check(rows, "(e) ||d_f L20|| S2A_k20, all 1836 factors, vs south7/stage_a1_observed.csv", t7, dn_20,
          extra=True)
    sample = pd.DataFrame(dict(factor=[keep[i] for i in pick], n_cells=n_cells[pick],
                               dnorm_S1=dn_s1[pick], stored_S1=t6[pick],
                               dnorm_S2A_k20=dn_20[pick], stored_S2A_k20=t7[pick]))
    ctx = dict(s1=s1, d_f=d_f, keep=keep, n_cells=n_cells, L=L, mr=mr, c_old=c_old,
               c_young=c_young, r_full=r_full, r20=r20, D20=D20, sample=sample, ortho=ortho)
    return ctx


def stage0_panel(ctx, log, rows):
    """(f): the SOUTH6 panel and NT counts rebuilt from stored inputs."""
    from scipy import sparse
    from fibro2_common import load_frozen_ruler
    from seng_run import tmm_cache_from_panel
    from toward_run import tmm_logcpm_quiet
    import south6_run as S6

    s1 = ctx["s1"]
    mr = np.asarray(s1["matched_ruler"], np.int64)
    check(rows, "(f) stage1.npz matched_ruler equals basis.npz matched_ruler", ctx["mr"], mr, "exact",
          extra=True)
    C_f = np.asarray(s1["C_f"], np.float64)
    C_nt = np.asarray(s1["C_nt"], np.float64)
    C_pool = np.asarray(s1["C_pool"], np.float64)
    C_panel = np.vstack([C_f, C_nt, C_pool.reshape(1, -1)])
    n_f = C_f.shape[0]
    del C_f
    t0 = time.time()
    logcpm, nf = tmm_logcpm_quiet(C_panel)
    frozen = load_frozen_ruler()
    mu = np.asarray(frozen["mu"], float)[mr]
    sd = np.asarray(frozen["sd"], float)[mr]
    sd_safe = np.where(sd < 1e-12, 1.0, sd)
    missing = (C_panel == 0).all(0)
    cache = tmm_cache_from_panel(C_panel, nf)
    del C_panel
    logcpm -= mu
    logcpm /= sd_safe
    logcpm[:, missing] = 0.0
    d_m = logcpm[:n_f] - logcpm[-1]
    del logcpm
    log(f"[panel] rebuilt from stored counts ({time.time() - t0:.0f}s)")
    check(rows, "(f) missing-gene mask identical to stage1.npz", np.asarray(s1["missing"], bool), missing,
          "exact", extra=True)
    check(rows, "(f) rebuilt d_f on matched genes vs stage1.npz d_f", ctx["d_f"][:, mr], d_m, extra=True)
    del d_m

    nt_cells = np.asarray(s1["nt_cells"], np.int64)
    nnz = int(np.asarray(s1["nnz"]).ravel()[0])
    data, indices, indptr = S6.open_matched(nnz)
    t0 = time.time()
    lens = (indptr[nt_cells + 1] - indptr[nt_cells]).astype(np.int64)
    ip = np.zeros(nt_cells.size + 1, np.int64)
    ip[1:] = np.cumsum(lens)
    dat = np.empty(int(ip[-1]), np.float64)
    ind = np.empty(int(ip[-1]), np.int32)
    for j, c in enumerate(nt_cells):
        a, b = int(indptr[c]), int(indptr[c + 1])
        dat[ip[j]:ip[j + 1]] = data[a:b]
        ind[ip[j]:ip[j + 1]] = indices[a:b]
    X_nt = sparse.csr_matrix((dat, ind, ip), shape=(nt_cells.size, mr.size))
    del data, indices
    log(f"[nt csr] cells={nt_cells.size} nnz={X_nt.nnz} ({time.time() - t0:.0f}s)")
    check(rows, "(f) X_nt rows (non-targeting cells)", N_NT, X_nt.shape[0], "exact", extra=True)
    check(rows, "(f) X_nt column sums equal stage1.npz C_pool", C_pool, np.asarray(X_nt.sum(0)).ravel(),
          "exact", extra=True)
    ctx.update(X_nt=X_nt, C_pool=C_pool, cache=cache, mu=mu, sd=sd, missing=missing, mr_s1=mr)


def null_draws_all(ctx, log, rows):
    """(g) and the draws themselves: SOUTH7's depth-matched NT nulls at every
    distinct factor cell count. Keeps projected draws, full squared norms and
    the full null-bank rows; checks SOUTH7's stored deltas."""
    from seng_run import z_from_counts
    from south6_run import block_delta
    from south6rec_run import indicator
    from south7_run import null_draws

    X_nt, C_pool, cache = ctx["X_nt"], ctx["C_pool"], ctx["cache"]
    mu, sd, missing, mr = ctx["mu"], ctx["sd"], ctx["missing"], ctx["mr"]
    L = ctx["L"]
    Lm = L[mr]
    r20 = ctx["r20"]
    r_full = ctx["r_full"]
    r2_20, rn_20 = float(r20 @ r20), float(np.linalg.norm(r20))
    r2_1, rn_1 = float(r_full @ r_full), float(np.linalg.norm(r_full))
    r1_m = r_full[mr]
    n_cells = ctx["n_cells"]
    N = X_nt.shape[0]

    occ = np.zeros(N_FACTORS, np.int64)          # m_i: draw index for factor i
    seen = {}
    for i, n in enumerate(n_cells):
        occ[i] = seen.get(int(n), 0)
        seen[int(n)] = occ[i] + 1
    ns = sorted(seen)
    log(f"[draws] distinct n = {len(ns)} ({ns[0]}-{ns[-1]}); max multiplicity {int(occ.max()) + 1}")

    def zm(counts):
        z, _, _, _ = z_from_counts(np.asarray(counts, np.float64).ravel(), cache, mu, sd, missing)
        return z

    stored_csv = pd.read_csv(SOUTH7 / "stage_a1_depth_nulls.csv")
    stored_npz = np.load(SOUTH7 / "stage_a1_null_draws.npz")
    E20_all = {}
    sq20 = {}
    sqfull = {}
    C_full = np.zeros((N_FACTORS, mr.size), np.float64)
    C20 = np.zeros((N_FACTORS, KDIM))
    sum_rows = []
    per_draw_fail = {"S2A_k20": [0, 0.0], "S1": [0, 0.0]}
    t0 = time.time()
    for j, n in enumerate(ns):
        draws, perm_of, info = null_draws(n, N)
        CD = np.asarray((indicator(draws, N) @ X_nt).todense())
        E = np.empty_like(CD)
        for i in range(CD.shape[0]):
            E[i] = zm(CD[i]) - zm(C_pool - CD[i])
        del CD
        E20 = E @ Lm
        E20_all[n] = E20
        sq20[n] = np.einsum("ij,ij->i", E20, E20)
        sqfull[n] = np.einsum("ij,ij->i", E, E)
        fac = np.flatnonzero(n_cells == n)
        C_full[fac] = E[occ[fac]]
        C20[fac] = E20[occ[fac]]
        for space, A, rr, r2, rn in (("S2A_k20", E20, r20, r2_20, rn_20), ("S1", E, r1_m, r2_1, rn_1)):
            dl, *_ = block_delta(A, rr, r2, rn)
            st = stored_npz[f"{n}|{space}"]
            diff = np.abs(dl - st) / np.maximum(1.0, np.abs(st))
            per_draw_fail[space][0] += int(np.sum(~(diff <= TOL)))
            per_draw_fail[space][1] = max(per_draw_fail[space][1], float(diff.max()))
            sum_rows.append(dict(n_cells=n, space=space, floor_p05=float(np.percentile(dl, 5)),
                                 median=float(np.median(dl)), sd=float(dl.std(ddof=1)),
                                 p95=float(np.percentile(dl, 95))))
        del E
        if (j + 1) % 25 == 0 or j + 1 == len(ns):
            log(f"[draws] {j + 1}/{len(ns)} depths ({time.time() - t0:.0f}s)")
    stored_npz.close()
    new = pd.DataFrame(sum_rows)
    for space in ("S2A_k20", "S1"):
        st = stored_csv[stored_csv.space == space].set_index("n_cells").reindex(ns)
        nw = new[new.space == space].set_index("n_cells").reindex(ns)
        for col in ("floor_p05", "median", "sd", "p95"):
            check(rows, f"(g) SOUTH7 depth nulls {space} {col}, {len(ns)} depths "
                  "vs south7/stage_a1_depth_nulls.csv", st[col].to_numpy(), nw[col].to_numpy(), extra=True)
        nb, mx = per_draw_fail[space]
        rows.append(dict(item=f"(g) SOUTH7 per-draw deltas {space}, {len(ns)} x 200 draws vs "
                         "south7/stage_a1_null_draws.npz (n outside tol)", kind="tol", gated=True,
                         extra=True, n=len(ns) * 200, stored=0, reproduced=nb, max_abs_diff=np.nan,
                         max_scaled_diff=mx, n_mismatch=nb, ok=bool(nb == 0), note=""))
    ctx.update(E20_all=E20_all, sq20=sq20, sqfull=sqfull, C_full=C_full, C20=C20, occ=occ, ns=ns)


def gate_stage0(rows, log):
    rep = pd.DataFrame(rows)
    rep.to_csv(OUT / "stage0_reproduction.csv", index=False)
    failed = rep[(~rep.ok) & rep.gated]
    log(f"[stage0] checks={len(rep)} gated={int(rep.gated.sum())} failed={len(failed)}")
    if len(failed):
        msg = "; ".join(f"{r.item}: {r.n_mismatch} of {r.n} outside tol (max scaled {r.max_scaled_diff})"
                        for r in failed.itertuples())
        if failed["item"].str.startswith("(b)").any():
            raise Stop("basis_not_orthonormal", msg)
        raise Stop("prior_results_not_reproduced", msg)


# --------------------------------------------------------------------------
# main run
# --------------------------------------------------------------------------
def frac_rows(space, set_key, res, gated_set):
    rows = []
    for k in KS:
        r = res[str(k)]
        for m in ("greedy", "l1"):
            x = r[m]
            rows.append(dict(space=space, set=set_key, set_label=SET_LABEL[set_key], k=str(k), method=m,
                             frac=x["frac"], nnz=x["nnz"],
                             l1_nnz_hit=x.get("nnz_hit", np.nan) if m == "l1" else np.nan,
                             l1_hit_exact=x.get("hit") if m == "l1" else None,
                             l1_lam=x.get("lam", np.nan) if m == "l1" else np.nan,
                             gated_set=gated_set, abs_diff_greedy_l1=np.nan, jaccard=np.nan))
        rows.append(dict(space=space, set=set_key, set_label=SET_LABEL[set_key], k=str(k), method="gated",
                         frac=r["gated"], nnz=np.nan, l1_nnz_hit=np.nan, l1_hit_exact=None, l1_lam=np.nan,
                         gated_set=gated_set, abs_diff_greedy_l1=r["abs_diff"], jaccard=r["jaccard"]))
    u = res["unlimited"]["bvls"]
    rows.append(dict(space=space, set=set_key, set_label=SET_LABEL[set_key], k="unlimited", method="bvls",
                     frac=u["frac"], nnz=u["nnz"], l1_nnz_hit=np.nan, l1_hit_exact=None, l1_lam=np.nan,
                     gated_set=gated_set, abs_diff_greedy_l1=np.nan, jaccard=np.nan))
    return rows


def weight_rows(space, set_key, res, keep, n_cells, pct):
    rows = []
    for k in K_ALL:
        methods = (("greedy", res[k]["greedy"]), ("l1", res[k]["l1"])) if k != "unlimited" else \
            (("bvls", res[k]["bvls"]),)
        for m, x in methods:
            for j in np.flatnonzero(x["w"] > 0):
                rows.append(dict(space=space, set=set_key, k=k, method=m, factor=keep[j], weight=x["w"][j],
                                 n_cells=int(n_cells[j]), n_cells_pctile=float(pct[j])))
    return rows


def spearman_rows(space, set_key, res, n_cells):
    from scipy import stats as sps
    rows = []
    for k in K_ALL:
        methods = (("greedy", res[k]["greedy"]), ("l1", res[k]["l1"])) if k != "unlimited" else \
            (("bvls", res[k]["bvls"]),)
        for m, x in methods:
            sel = np.flatnonzero(x["w"] > 0)
            if sel.size >= 3:
                rho, p = sps.spearmanr(x["w"][sel], n_cells[sel])
            else:
                rho, p = np.nan, np.nan
            rows.append(dict(space=space, set=set_key, k=k, method=m, over="selected", n=int(sel.size),
                             spearman_rho=float(rho), p=float(p),
                             median_n_cells_selected=float(np.median(n_cells[sel])) if sel.size else np.nan,
                             median_n_cells_all=float(np.median(n_cells))))
            if k == "unlimited":
                rho, p = sps.spearmanr(x["w"], n_cells)
                rows.append(dict(space=space, set=set_key, k=k, method=m, over="all 1836 (zeros included)",
                                 n=int(x["w"].size), spearman_rho=float(rho), p=float(p),
                                 median_n_cells_selected=np.nan, median_n_cells_all=float(np.median(n_cells))))
    return rows


def gram_problem(Dm, b_full):
    """Full-space problem in Gram-reduced form. Dm rows = vectors (matched coords
    already scattered into b's coordinates by the caller)."""
    G = Dm @ Dm.T
    c = Dm @ b_full
    lam, U = np.linalg.eigh(G)
    keep = lam > EIG_REL * lam.max()
    lam, U = lam[keep], U[:, keep]
    sq = np.sqrt(lam)
    A = sq[:, None] * U.T
    bp = (U.T @ c) / sq
    bb = float(b_full @ b_full)
    c0 = bb - float(bp @ bp)
    return Problem(A, bp, np.sqrt(bb), c0=max(c0, 0.0), tol_scale=0.5 * bb), dict(rank=int(keep.sum()),
                                                                                    c0_raw=c0)


def readings_of(res_young, dtab):
    B = {k: float(res_young["B"][k]["gated"]) for k in K_ALL}
    C = {k: float(res_young["C"][k]["gated"]) for k in K_ALL}
    Dp5, Dp95, Dmed = {}, {}, {}
    for k in K_ALL:
        v = dtab[f"gated_k{k}"].to_numpy(np.float64)
        Dp5[k] = float(np.percentile(v, 5))
        Dp95[k] = float(np.percentile(v, 95))
        Dmed[k] = float(np.median(v))
    c_clause = [k for k in K_ALL if C[k] >= B[k] - TIE]
    d_clause = [k for k in K_ALL if Dp5[k] - TIE <= B[k] <= Dp95[k] + TIE]
    fired = "R0" if (c_clause or d_clause) else "R1"
    return dict(B=B, C=C, Dp5=Dp5, Dp95=Dp95, Dmed=Dmed, c_clause=c_clause, d_clause=d_clause,
                fired=fired)


def reading_text(rd):
    if rd["fired"] == "R0":
        parts = []
        if rd["c_clause"]:
            parts.append("the null bank (C) reaches at least as far as the shrunken real factors (B) at k = "
                         + ", ".join(rd["c_clause"]))
        if rd["d_clause"]:
            parts.append("the young target sits inside the middle 90% of the random-target distribution (D) "
                         "at k = " + ", ".join(rd["d_clause"]))
        return ("UNINTERPRETABLE. The bound cannot distinguish the measured factors from noise vectors of the "
                "same depth, so it carries no information about whether young is reachable, and no number "
                "from (A) or (B) may be quoted as a bound. Fired because: " + "; ".join(parts) + ".")
    B, C, p5, p95 = rd["B"], rd["C"], rd["Dp5"], rd["Dp95"]
    above = [k for k in K_ALL if B[k] > p95[k] + TIE]
    below = [k for k in K_ALL if B[k] < p5[k] - TIE]
    if len(above) == len(K_ALL):
        side = "above that range at every k"
    elif len(below) == len(K_ALL):
        side = "below that range at every k"
    else:
        side = ("above that range at k = " + ", ".join(above) + " and below it at k = " + ", ".join(below))
    return ("BOUND REPORTED. In S2A_k20, with Southard's factor displacements shrunk for depth noise (set B), "
            "the largest fraction of the O1-Y1 gap closed by any combination with doses between 0 and the "
            "delivered dose is " + ", ".join(f"k={k}: {fmt(B[k])}" for k in K_ALL) + ". Beside it: the null "
            "factor bank (set C) closes " + ", ".join(f"k={k}: {fmt(C[k])}" for k in K_ALL) + "; against 200 "
            "random targets of the same length (set D) the 5th-95th percentile range is "
            + ", ".join(f"k={k}: {fmt(p5[k])}-{fmt(p95[k])}" for k in K_ALL) + " and the young target sits "
            + side + ". This is an upper bound computed under assumptions known to be optimistic: that these "
            "directions are measured rather than estimated, that they add, and that they transfer to aged "
            "cells. It is sizing information for the combination screen. It is not evidence that any "
            "combination rejuvenates aged cells, and it must not be cited as such.")


def solve_gated(P, label, log):
    t0 = time.time()
    res, ok, info = run_set(P)
    log(f"[{label}] " + ", ".join(f"k={k}: {res[k]['gated']:.4f}" for k in K_ALL)
        + f" ok={ok} ({time.time() - t0:.0f}s; L1 solves {info['l1']['n_solves']})")
    if not ok:
        raise Stop("solver_not_converged",
                   f"{label}: greedy ok={info['ok_greedy']}, L1 ok={info['ok_l1']}, BVLS ok={info['ok_bvls']} "
                   f"(max gaps greedy {info['greedy']['max_gap']:.3e}, L1 {info['l1']['max_gap']:.3e}, "
                   f"BVLS {info['unlimited_gap']:.3e}; tolerance {P.tol:.3e}).")
    return res, info


def run(log):
    if not PREREG.exists():
        raise SystemExit(f"{PREREG} missing. Not computing.")
    rows = []
    ctx = stage0_static(log, rows)
    gate_stage0(rows, log)
    stage0_panel(ctx, log, rows)
    gate_stage0(rows, log)
    null_draws_all(ctx, log, rows)
    gate_stage0(rows, log)
    log("[stage0] all gated checks pass")

    keep, n_cells, L = ctx["keep"], ctx["n_cells"], ctx["L"]
    D20 = ctx["D20"]
    b20 = -ctx["r20"]
    g = float(np.linalg.norm(b20))
    sigma2_n = {n: float(ctx["sq20"][n].mean()) for n in ctx["ns"]}
    sigma2 = np.array([sigma2_n[int(n)] for n in n_cells])
    dn2 = np.einsum("ij,ij->i", D20, D20)
    s_f = np.maximum(0.0, 1.0 - sigma2 / dn2)
    B20 = D20 * s_f[:, None]
    C20 = ctx["C20"]
    cn2 = np.einsum("ij,ij->i", C20, C20)
    s_c = np.maximum(0.0, 1.0 - sigma2 / cn2)
    Cs20 = C20 * s_c[:, None]
    z = np.random.default_rng([SEED, 4]).standard_normal((N_TARGETS, KDIM))
    targets = z / np.linalg.norm(z, axis=1, keepdims=True) * g

    pd.DataFrame(dict(n_cells=ctx["ns"], sigma2_S2A_k20=[sigma2_n[n] for n in ctx["ns"]],
                      sigma2_S1=[float(ctx["sqfull"][n].mean()) for n in ctx["ns"]], R=200)
                 ).to_csv(OUT / "null_depths.csv", index=False)
    pct = pd.Series(n_cells).rank(pct=True, method="average").to_numpy() * 100.0
    fac = pd.DataFrame(dict(factor=keep, n_cells=n_cells, n_cells_pctile=pct, dnorm_S2A_k20=np.sqrt(dn2),
                            sigma2_S2A_k20=sigma2, s_f=s_f, null_draw_index=ctx["occ"],
                            null_dnorm_S2A_k20=np.sqrt(cn2), null_s=s_c))
    fac.to_csv(OUT / "factors.csv", index=False)
    np.savez(OUT / "inputs_k20.npz", keep=np.array(keep, dtype=object), n_cells=n_cells, D20=D20, b20=b20,
             g=np.array([g]), sigma2=sigma2, s_f=s_f, C20=C20, targets=targets)
    ctx["sample"].to_csv(OUT / "stage0_norm_sample.csv", index=False)
    log(f"[inputs] g = {g:.6f}; s_f = 0 for {int((s_f == 0).sum())} of {N_FACTORS} factors; "
        f"null bank s = 0 for {int((s_c == 0).sum())}")

    # Solver check: (B) plus a planted column equal to b.
    Pp = Problem(np.hstack([B20.T, b20[:, None]]), b20, g)
    gr, okg, _ = greedy(Pp)
    lr, okl, _ = l1(Pp)
    sel_g = np.flatnonzero(gr[1]["w"] > 0).tolist()
    sel_l = np.flatnonzero(lr[1]["w"] > 0).tolist()
    plant = dict(greedy_frac=gr[1]["frac"], greedy_selected=sel_g, l1_frac=lr[1]["frac"], l1_selected=sel_l,
                 ok_greedy=okg, ok_l1=okl)
    plant_ok = bool(okg and okl and gr[1]["frac"] >= 1 - TIE and lr[1]["frac"] >= 1 - TIE
                    and sel_g == [N_FACTORS] and sel_l == [N_FACTORS])
    plant["passed"] = plant_ok
    dump("solver_check.json", plant)
    log(f"[solver check] {plant}")
    if not plant_ok:
        raise Stop("solver_not_converged", f"planted-column solver check failed: {plant}")

    res = {}
    infos = {}
    for key, M in (("A", D20), ("B", B20), ("C", C20)):
        res[key], infos[key] = solve_gated(Problem(M.T, b20, g), f"{key} young", log)
    dtab = run_targets(B20.T, g, targets, log)
    dtab.to_csv(OUT / "random_targets.csv", index=False)
    rd = readings_of(res, dtab)
    rd["text"] = reading_text(rd)
    rd["g_S2A_k20"] = g
    rd["solver_check_passed"] = plant_ok
    rd["l1_solves"] = {k: v["l1"]["n_solves"] for k, v in infos.items()}
    dump("readings.json", rd)
    log(f"[reading] {rd['fired']} c_clause={rd['c_clause']} d_clause={rd['d_clause']}")

    # ---------------- diagnostics (not gated) ----------------
    log("[diag] C shrunk")
    res_cs, ok_cs, info_cs = run_set(Problem(Cs20.T, b20, g))
    ft = []
    for key in ("A", "B", "C"):
        ft += frac_rows("S2A_k20", key, res[key], True)
    ft += frac_rows("S2A_k20", "Cs", res_cs, False)
    wr, sr = [], []
    for key in ("A", "B"):
        wr += weight_rows("S2A_k20", key, res[key], keep, n_cells, pct)
        sr += spearman_rows("S2A_k20", key, res[key], n_cells)

    log("[diag] full gene space S1")
    mr = ctx["mr"]
    b_full = (ctx["c_young"] - ctx["c_old"])
    full_info = {}
    sigma2_full = np.array([float(ctx["sqfull"][int(n)].mean()) for n in n_cells])
    Dm = ctx["d_f"][:, mr]
    # b's off-matched coordinates are untouchable; fold them into c0 by
    # working on matched coordinates and adding their squared norm back.
    off = np.ones(N_RULER, bool)
    off[mr] = False
    b_m = b_full[mr]
    c_off = float(b_full[off] @ b_full[off])
    dn2_full = np.einsum("ij,ij->i", Dm, Dm)
    s_full = np.maximum(0.0, 1.0 - sigma2_full / dn2_full)
    C_full = ctx["C_full"]
    for key, M in (("A", Dm), ("B", Dm * s_full[:, None]), ("C", C_full)):
        t0 = time.time()
        P, gi = gram_problem(M, b_m)
        P.c0 += c_off
        P.g = float(np.linalg.norm(b_full))
        P.tol = GAP_REL * max(1.0, 0.5 * P.g ** 2)
        r_, ok_, inf_ = run_set(P)
        full_info[key] = dict(ok=ok_, rank=gi["rank"], c0_raw_matched=gi["c0_raw"], c_off=c_off,
                              s_zero=int((s_full == 0).sum()) if key == "B" else None,
                              seconds=time.time() - t0)
        ft += [dict(r, space="S1") for r in frac_rows("S1", key, r_, False)]
        if key in ("A", "B"):
            sr += spearman_rows("S1", key, r_, n_cells)
        log(f"[diag S1 {key}] " + ", ".join(f"k={k}: {r_[k]['gated']:.4f}" for k in K_ALL)
            + f" ok={ok_} rank={gi['rank']} ({time.time() - t0:.0f}s)")
    pd.DataFrame(ft).to_csv(OUT / "frac_table.csv", index=False)
    pd.DataFrame(wr).to_csv(OUT / "weights.csv", index=False)
    pd.DataFrame(sr).to_csv(OUT / "spearman.csv", index=False)
    dump("diagnostics.json", dict(
        c_shrunk_ok=ok_cs, full_space=full_info, s_f_zero_B=int((s_f == 0).sum()),
        s_zero_null_bank=int((s_c == 0).sum()), n_factors=N_FACTORS,
        s_f_zero_B_full=int((s_full == 0).sum()),
        gaps={k: dict(greedy=v["greedy"]["max_gap"], l1=v["l1"]["max_gap"], bvls=v["unlimited_gap"])
              for k, v in infos.items()},
        d_max_gap=float(dtab["max_gap"].max()),
        n_workers=N_WORKERS))
    log("[run] DONE")


# --------------------------------------------------------------------------
# verify (tracked files only)
# --------------------------------------------------------------------------
def verify(log):
    z = np.load(OUT / "inputs_k20.npz", allow_pickle=True)
    D20, b20, g = z["D20"], z["b20"], float(z["g"][0])
    B20 = D20 * z["s_f"][:, None]
    res = {}
    for key, M in (("A", D20), ("B", B20), ("C", z["C20"])):
        res[key], _ = solve_gated(Problem(M.T, b20, g), f"verify {key}", log)
    dtab = run_targets(B20.T, g, z["targets"], log)
    rd = readings_of(res, dtab)
    st = json.loads((OUT / "readings.json").read_text(encoding="utf-8"))
    ft = pd.read_csv(OUT / "frac_table.csv")
    bad = []
    for key in ("A", "B", "C"):
        for k in K_ALL:
            m = "bvls" if k == "unlimited" else "gated"
            sv = ft[(ft.space == "S2A_k20") & (ft.set == key) & (ft.k == k) & (ft.method == m)]["frac"].iloc[0]
            if abs(sv - res[key][k]["gated"]) > TOL * max(1, abs(sv)):
                bad.append(f"{key} k={k}: stored {sv} recomputed {res[key][k]['gated']}")
    dst = pd.read_csv(OUT / "random_targets.csv")
    for k in K_ALL:
        d = np.abs(dst[f"gated_k{k}"].to_numpy() - dtab[f"gated_k{k}"].to_numpy())
        if np.any(d > TOL):
            bad.append(f"D k={k}: {int((d > TOL).sum())} targets differ")
    if rd["fired"] != st["fired"]:
        bad.append(f"reading stored {st['fired']} recomputed {rd['fired']}")
    log(f"[verify] reading {rd['fired']} (stored {st['fired']}); mismatches: {bad if bad else 'none'}")
    return not bad


# --------------------------------------------------------------------------
# findings
# --------------------------------------------------------------------------
def md_table(df, cols, heads=None):
    heads = heads or cols
    out = ["| " + " | ".join(heads) + " |", "|" + "---|" * len(cols)]
    for _, r in df.iterrows():
        cells = []
        for c in cols:
            v = r[c]
            if isinstance(v, (float, np.floating)):
                cells.append(fmt(v, 4) if np.isfinite(v) else "")
            else:
                cells.append(str(v))
        out.append("| " + " | ".join(cells) + " |")
    return "\n".join(out)


DEVIATIONS = [
    "**S2A_k20 O1-Y1 distance compared against `results/south7/basis_retention.csv`, not "
    "`results/south6/origin_target_distances.csv`.** The task named the SOUTH6 file for both S1 and S2A_k20; "
    "that file has no S2A row (its k20 row is SOUTH6's S2_k20 basis, 28.645). Fixed in the flag before any "
    "number; S1 is compared against the SOUTH6 file as asked.",
    "**sigma_f^2 built with this task's code.** `results/south6/stage2a_floors.npz` holds gene-label permutation "
    "delta floors (keys `floors`, `les`, `pooled`, `medians`), not indexed by cell count and not squared norms, "
    "so it does not supply sigma_f^2. sigma^2(n) is the mean of ||E L20||^2 over the 200 depth-matched "
    "non-targeting draws of `south7_run.null_draws(n, 14013)` at each factor's own n; no interpolation.",
    "**\"The SOUTH6 control-draw construction\" read as the SOUTH6REC/SOUTH7 depth-matched construction.** "
    "`src/south6_run.py` has no draw at a factor's cell count (its control draws are 7,006/7,007-cell half "
    "splits and leave-one-guide-out). The construction used is SOUTH6REC's primary scheme, which SOUTH7 used "
    "for its depth-matched floor: E = z(draw) - z(C_pool - draw) with the SOUTH6 panel TMM cache. Its "
    "reproduction of SOUTH7's stored null deltas is gated in Stage 0 (g).",
    "**(C) is not shrunk.** The task did not say; the flag fixed (C) as the null analogue of (A), so R0 compares "
    "unshrunk noise with shrunk real factors (the stricter comparison for the bound). (C) shrunk by the same "
    "rule is a labelled diagnostic.",
    "**(C) reuses the sigma^2 draws.** Factor i with count n takes draw index m_i (its rank among factors with "
    "that n, keep order; max 32 of 200). The same draws estimate sigma^2(n), which enters only (B).",
    "**Decisions the task did not specify, fixed in the flag:** gated value at finite k = max(greedy, L1-refit); "
    "L1 support refit by BVLS without penalty; L1 bisection rule and fallback when nnz never equals k; "
    "nnz counts weights > 1e-6 (weights at or below 1e-6 are zeroed after every BVLS fit); convergence "
    "certified by duality gap <= 1e-8 * max(1, 0.5 g^2); R0 tie tolerance 1e-6 with a tie (including both "
    "sets reaching the target exactly) counting as \"at least as far\"; \"at any k\" applied to both R0 "
    "clauses; the D clause inclusive with 1e-6 slack; random targets from default_rng([20260914, 4]); "
    "Stage 0 sample from default_rng([20260914, 0]).",
    "**Added beyond the task text:** the planted-column solver check (the task defines no positive control); "
    "gated extra Stage 0 checks (e)-(g) (all-factor norms, panel rebuild from stored counts, reproduction of "
    "SOUTH7's stored null deltas); `stage0_norm_sample.csv`, `null_depths.csv`, `factors.csv`, "
    "`inputs_k20.npz`, `solver_check.json`, `diagnostics.json`, `verify` mode.",
    "**Imported unchanged:** `south6_run` (block_delta, open_matched), `south6rec_run.indicator`, "
    "`south7_run.null_draws`, `seng_run` (tmm_cache_from_panel, z_from_counts), `toward_run.tmm_logcpm_quiet`, "
    "`fibro2_common.load_frozen_ruler`. Importing `south7_run` runs its `OUT.mkdir` on `results/south7/` "
    "(already present); nothing is written there.",
    "**Where written:** only `results/reach/`, `src/reach_run.py`, `FINDINGS_REACH.md`. "
    "`results/transfer/NOT_A_RESULT_gate_failed/` was not opened. The pool for set (D) runs 16 single-threaded "
    "workers; BLAS thread count can change floating-point summation order, not results beyond ~1e-12.",
    "**Script edited while the run was in progress.** The findings renderer (`render`, `_stop_section`, the "
    "lists in this section) was added to `src/reach_run.py` after the run started, at Stage 0 (f)/(g). The "
    "main process had already loaded the earlier file; set (D)'s spawned workers re-import the script from "
    "disk and so loaded the edited file. No solver, draw or Stage 0 function was changed by those edits; "
    "the committed file is the edited one.",
    "**Flag wording slip, not edited:** the flag says S1 includes \"the 3,221 ruler genes the factors cannot "
    "move\". 3,221 is the unmatched count (23,485 - 20,264); genes in the missing-gene mask (zero in every panel "
    "row) are also z = 0 in every d_f, so the untouchable set is at least that large.",
]

WHAT_NOT = [
    "Southard is Hs27 neonatal foreskin fibroblasts, not aged donors. Nothing here transfers to aged cells.",
    "The transfer check (FINDINGS_TRANSFER.md, commit 7230271) could not measure whether these directions carry "
    "to another fibroblast context.",
    "Doses above the delivered one are not modelled; any w in (0, 1) assumes a partial dose acts as a scaled "
    "displacement, itself an assumption. Doses above 1 would be extrapolation.",
    "Additivity of displacements is assumed, not tested.",
    "Nothing here models toxicity, identity loss or proliferation.",
    "Circularity, disclosed: the null draws, sigma^2 and (C) all come from the same 14,013 non-targeting cells "
    "that form d_f's reference pool; (C)'s vectors are among the draws that estimate sigma^2. The S2A basis "
    "was fit on GTEx donors aged 40-59 with no anchor donor in the fit (SOUTH7), so the space is not fit to "
    "the factors or the anchors.",
]


def _stop_section(stop):
    tol = GAP_REL * max(1.0, 0.5 * 27.929154451474936 ** 2)
    sc = json.loads((OUT / "solver_check.json").read_text(encoding="utf-8")) if (OUT / "solver_check.json").exists() else {}
    return [
        "## What happened, in order", "",
        "1. Stage 0 (a)-(g): 30 of 30 gated checks passed (table above), including the reproduction of "
        "SOUTH7's stored depth-matched null deltas for all 212 cell counts.",
        "2. Inputs were built and written before any solve: `null_depths.csv` (sigma^2(n)), `factors.csv` "
        "(per-factor ||D20_f||, sigma^2, s_f, null-bank draw index and norm), `inputs_k20.npz`.",
        f"3. Planted-column solver check passed (greedy and L1 at k=1 both selected the planted column; "
        f"passed = {sc.get('passed')}). `solver_check.json`.",
        "4. Sets (A), (B), (C) on the young target were solved and each passed the convergence test. Their "
        "values were printed to `run.log` as the code logs every set. **They are not results.** No reading "
        "was evaluated: both R0 clauses and R1 need set (D), and only a pre-registered reading whose "
        "conditions fired is a conclusion. They are not tabulated here and are not a bound.",
        "5. Set (D): 200 random targets in a pool of 16 workers, results consumed in target order. Target 11 "
        "(0-based) returned a failed convergence test: the largest duality gap over that target's solves was "
        f"4.268e-06 against the pre-registered tolerance 1e-8 x max(1, 0.5 g^2) = {tol:.3e} (ratio 1.09). "
        "The pool was terminated and STOP `solver_not_converged` written. Which of that target's solves failed "
        "(a greedy BVLS refit, an L1 L-BFGS-B solve after 20 restarts, the L1 refit, or the unlimited BVLS) is "
        "not recorded: the worker returns only the maximum gap. Targets 0-10 had completed; results of targets "
        "in flight on other workers were discarded unrecorded; `random_targets.csv` was not written.",
        "6. After the STOP nothing was computed: no rerun, no change of tolerance or solver, no "
        "re-solve of target 11 to find which step failed.", "",
        "## What the STOP does and does not say", "",
        "- It is a failure of the convergence certificate fixed in the flag, which was chosen from solver "
        "behaviour on synthetic matrices. The gap exceeded it by 9%. It is not evidence about reachability "
        "in either direction, and it leaves no reading on record.",
        "- No frac_closed table across the four sets is reported, because (D) was not completed and the "
        "readings were never evaluated.",
        "- Whether a new task should re-pose the bound with a different, pre-registered convergence rule is a "
        "decision for outside this session. This run is not to be resumed or adjusted.",
        "- The STOP is reproducible from tracked files: `python src/reach_run.py verify` re-solves from "
        "`inputs_k20.npz` and runs the same convergence test (not run in this session, as it would compute "
        "after the STOP).", "",
    ]


def _prereg_block():
    return "```\n" + PREREG.read_text(encoding="utf-8").rstrip("\n") + "\n```"


def _stage0_section():
    rep = pd.read_csv(OUT / "stage0_reproduction.csv")
    L = ["## Stage 0 — reproduction of stored numbers", "",
         "Counts exact; numerics within 1e-4 x max(1, |stored|). Every row is gated; rows marked `extra` were "
         "added beyond the task text. For multi-value rows `stored`/`reproduced` are blank and `n_mismatch` "
         "counts values outside tolerance. File: `results/reach/stage0_reproduction.csv`.", ""]
    rep = rep.copy()
    rep["extra"] = rep["extra"].map(lambda x: "extra" if str(x) == "True" else "")
    L.append(md_table(rep.fillna(""), ["item", "kind", "extra", "n", "stored", "reproduced", "max_scaled_diff",
                                        "n_mismatch", "ok"]))
    n_fail = int(((rep.ok.astype(str) != "True") & (rep.gated.astype(str) == "True")).sum())
    L += ["", f"{len(rep)} checks, {int((rep.gated.astype(str) == 'True').sum())} gated, {n_fail} gated failures."]
    p = OUT / "stage0_norm_sample.csv"
    if p.exists():
        s = pd.read_csv(p)
        L += ["", "The 20 sampled factors for (d) (default_rng([20260914, 0])): "
              + ", ".join(f"{a} ({n})" for a, n in zip(s.factor, s.n_cells)) + " — name (cells). "
              "Values: `results/reach/stage0_norm_sample.csv`."]
    return L


def render(stop=None):
    rd = json.loads((OUT / "readings.json").read_text(encoding="utf-8")) if (OUT / "readings.json").exists() else None
    L = ["# FINDINGS_REACH — the reachability bound", ""]
    if stop:
        L += [f"**Status: STOP `{stop['key']}`.** {stop['message']}", "",
              "Per CLAUDE.md, nothing was computed after the STOP and the run was not restarted.", ""]
    else:
        L += [f"**Status: DONE.** Reading: `{'UNINTERPRETABLE' if rd['fired'] == 'R0' else 'BOUND REPORTED'}` "
              f"({rd['fired']}).", "",
              "This is sizing information, never evidence (TF search model blueprint). No reading of this task "
              "licenses a claim about aged cells. A positive bound sizes an experiment; it does not support one.",
              ""]
    L += ["Reproduced by `python src/reach_run.py` then `python src/reach_run.py findings`. "
          + ("`python src/reach_run.py verify` re-solves sets A-D from tracked `results/reach/inputs_k20.npz` "
             "and applies the same convergence test, so the STOP is re-checkable from tracked files." if stop else
             "The reading is re-derived from tracked files alone by `python src/reach_run.py verify` (re-solves "
             "sets A-D from `results/reach/inputs_k20.npz`)."), "",
          "## Pre-registration (verbatim from `results/reach/PREREG.flag`)", "", _prereg_block(), ""]
    L += ["## Deviations and decisions, disclosed", ""] + [f"{i + 1}. {d}" for i, d in enumerate(DEVIATIONS)] + [""]
    L += _stage0_section() + [""]
    if stop or rd is None:
        L += _stop_section(stop)
        L += ["## What this does not cover", ""] + [f"- {w}" for w in WHAT_NOT] + [""]
        present = sorted(p.name for p in OUT.iterdir() if p.is_file())
        L += ["## Files", "", "`src/reach_run.py`; `results/reach/`: " + ", ".join(f"`{p}`" for p in present)
              + ". Not written, because the STOP fired first: `random_targets.csv`, `readings.json`, "
              "`frac_table.csv`, `weights.csv`, `spearman.csv`, `diagnostics.json`."]
        return "\n".join(L) + "\n"

    sc = json.loads((OUT / "solver_check.json").read_text(encoding="utf-8"))
    L += ["## Solver check (added; not a positive control for a test)", "",
          f"Set (B)'s 1,836 vectors plus one planted column equal to the gap b. Greedy k=1: frac "
          f"{fmt(sc['greedy_frac'], 6)}, selected column {sc['greedy_selected']}; L1 k=1: frac "
          f"{fmt(sc['l1_frac'], 6)}, selected {sc['l1_selected']} (the planted column is index 1836). "
          f"Passed: {sc['passed']}. This shows only that the optimiser finds a direction that is present; "
          "the bound is an optimisation, not a test, so it has no recovery power in the usual sense. Sets "
          "(C) and (D) are its comparisons.", ""]

    B, C, p5, p95, med = rd["B"], rd["C"], rd["Dp5"], rd["Dp95"], rd["Dmed"]
    fired = rd["fired"]
    L += [f"## Pre-registered reading: `{'UNINTERPRETABLE' if fired == 'R0' else 'BOUND REPORTED'}` ({fired})",
          "", "> " + rd["text"], "",
          "No reading of this task licenses a claim about aged cells. A positive bound sizes an experiment; it "
          "does not support one.", "",
          "Conditions, per k (gated values; finite k = max(greedy, L1-refit), unlimited = exact BVLS):", ""]
    cond = pd.DataFrame(dict(
        k=list(K_ALL), B=[B[k] for k in K_ALL], C=[C[k] for k in K_ALL],
        C_ge_B=[C[k] >= B[k] - TIE for k in K_ALL],
        D_p5=[p5[k] for k in K_ALL], D_p95=[p95[k] for k in K_ALL],
        B_inside_D=[p5[k] - TIE <= B[k] <= p95[k] + TIE for k in K_ALL]))
    L += [md_table(cond, ["k", "B", "C", "C_ge_B", "D_p5", "D_p95", "B_inside_D"],
                   ["k", "B shrunk", "C null bank", "C >= B - 1e-6", "D 5th pct", "D 95th pct",
                    "B inside D middle 90%"]), ""]
    if fired == "R0":
        L += ["Under R0 the numbers from (A) and (B) are not a bound and must not be quoted as one. They appear "
              "below only because the hand-off requires the table across sets; they are labelled as such.", ""]

    ft = pd.read_csv(OUT / "frac_table.csv")
    dt = pd.read_csv(OUT / "random_targets.csv")
    k20 = ft[ft.space == "S2A_k20"]

    def gv(s, k, space="S2A_k20"):
        m = "bvls" if k == "unlimited" else "gated"
        x = ft[(ft.space == space) & (ft.set == s) & (ft.k == k) & (ft.method == m)]
        return float(x["frac"].iloc[0]) if len(x) else np.nan

    tab = pd.DataFrame(dict(
        k=list(K_ALL), A=[gv("A", k) for k in K_ALL], B=[gv("B", k) for k in K_ALL],
        C=[gv("C", k) for k in K_ALL], Dmed=[med[k] for k in K_ALL], Dp5=[p5[k] for k in K_ALL],
        Dp95=[p95[k] for k in K_ALL],
        Dmin=[float(dt[f"gated_k{k}"].min()) for k in K_ALL],
        Dmax=[float(dt[f"gated_k{k}"].max()) for k in K_ALL]))
    L += ["## frac_closed across k, all four sets (S2A_k20)", "",
          f"g = ||c_old - c_young|| in S2A_k20 = {fmt(rd['g_S2A_k20'], 6)}. Gated value per set and k. (D) is "
          "(B)'s vectors against 200 random targets of length g; its columns summarise those 200 values. "
          + ("Under the R0 reading that fired, none of these numbers is a bound." if fired == "R0" else
             "(A) is the unshrunk diagnostic (not gated)."), "",
          md_table(tab, ["k", "A", "B", "C", "Dmed", "Dp5", "Dp95", "Dmin", "Dmax"],
                   ["k", "A raw", "B shrunk", "C null bank", "D median", "D 5th", "D 95th", "D min", "D max"]),
          "", "Files: `results/reach/frac_table.csv` (every set, k and method), `results/reach/random_targets.csv` "
          "(every target, k and method).", ""]

    rows = []
    for s in ("A", "B", "C", "Cs"):
        for k in KS:
            g_ = k20[(k20.set == s) & (k20.k == str(k))]
            G = g_[g_.method == "greedy"].iloc[0]
            Lr = g_[g_.method == "l1"].iloc[0]
            Q = g_[g_.method == "gated"].iloc[0]
            rows.append(dict(set=SET_LABEL[s], k=k, greedy=G.frac, greedy_nnz=int(G.nnz), l1=Lr.frac,
                             l1_nnz=int(Lr.nnz), l1_hit=f"{int(Lr.l1_nnz_hit)}{'' if str(Lr.l1_hit_exact) == 'True' else ' (not exact)'}",
                             diff=Q.abs_diff_greedy_l1, jac=Q.jaccard))
    mt = pd.DataFrame(rows)
    dis = []
    for k in KS:
        d = np.abs(dt[f"greedy_k{k}"] - dt[f"l1_k{k}"])
        dis.append(dict(k=k, med=float(d.median()), mx=float(d.max()),
                        l1_better=int((dt[f"l1_k{k}"] > dt[f"greedy_k{k}"] + TIE).sum()),
                        gr_better=int((dt[f"greedy_k{k}"] > dt[f"l1_k{k}"] + TIE).sum()),
                        miss=int((dt[f"l1_hit_k{k}"].astype(str) != "True").sum())))
    L += ["## Greedy versus L1 at finite k", "",
          "Both methods are feasible solutions; neither is guaranteed optimal for k > 1 (greedy is exhaustive "
          "at k = 1). `l1 nnz hit` is the support size the bisection reached before the BVLS refit; "
          "\"not exact\" means no lam in 60 steps gave exactly k and the fallback rule was used.", "",
          md_table(mt, ["set", "k", "greedy", "greedy_nnz", "l1", "l1_nnz", "l1_hit", "diff", "jac"],
                   ["set", "k", "greedy frac", "greedy nnz", "L1 frac", "L1 nnz", "L1 nnz hit",
                    "abs diff", "support Jaccard"]), "",
          "Set (D), over the 200 targets:", "",
          md_table(pd.DataFrame(dis), ["k", "med", "mx", "l1_better", "gr_better", "miss"],
                   ["k", "median abs diff", "max abs diff", "L1 better (targets)", "greedy better (targets)",
                    "L1 did not hit k exactly"]), ""]

    dg = json.loads((OUT / "diagnostics.json").read_text(encoding="utf-8"))
    fac = pd.read_csv(OUT / "factors.csv")
    w = pd.read_csv(OUT / "weights.csv")
    sp = pd.read_csv(OUT / "spearman.csv")
    L += ["## Diagnostics (not gated)", "",
          "### What shrinkage costs: (A) beside (B)", "",
          md_table(tab, ["k", "A", "B"], ["k", "A raw", "B shrunk"]), "",
          f"Factors shrunk to zero in (B) (||D20_f||^2 <= sigma^2(n_f)): {dg['s_f_zero_B']} of {dg['n_factors']}. "
          f"Chance baseline: the same rule shrinks {dg['s_zero_null_bank']} of the {dg['n_factors']} null-bank "
          "vectors (C) to zero. Median s_f "
          f"{fmt(float(fac.s_f.median()))}; median ||D20_f|| {fmt(float(fac.dnorm_S2A_k20.median()))} against "
          f"median sigma(n_f) {fmt(float(np.sqrt(fac.sigma2_S2A_k20).median()))}.", "",
          "### (C) shrunk by the same rule", "",
          md_table(pd.DataFrame(dict(k=list(K_ALL), Cs=[gv("Cs", k) for k in K_ALL], C=[gv("C", k) for k in K_ALL],
                                     B=[gv("B", k) for k in K_ALL])),
                   ["k", "Cs", "C", "B"], ["k", "C shrunk (diag)", "C null bank", "B shrunk"]), "",
          "### Which factors carry the weight, and cell count", "",
          "This project's shortlist-versus-depth check. Cell-count percentile is among the 1,836 factors "
          f"(median cell count {int(fac.n_cells.median())}).", ""]
    for s in ("B", "A"):
        for k in ("1", "3", "5", "10"):
            for m in ("greedy", "l1"):
                x = w[(w.space == "S2A_k20") & (w.set == s) & (w.k.astype(str) == k) & (w.method == m)]
                x = x.sort_values("weight", ascending=False)
                L.append(f"- {SET_LABEL[s]}, k={k}, {m}: " + (", ".join(
                    f"{r.factor} (w {r.weight:.2f}, {r.n_cells} cells, p{r.n_cells_pctile:.0f})"
                    for r in x.itertuples()) if len(x) else "none"))
    for s in ("B", "A"):
        for k, m in (("50", "greedy"), ("50", "l1"), ("unlimited", "bvls")):
            x = w[(w.space == "S2A_k20") & (w.set == s) & (w.k.astype(str) == k) & (w.method == m)]
            if len(x):
                L.append(f"- {SET_LABEL[s]}, k={k}, {m}: {len(x)} factors; median cell count "
                         f"{float(x.n_cells.median()):.0f} (all factors {float(fac.n_cells.median()):.0f}); "
                         f"{int((x.weight >= 1 - 1e-9).sum())} at the dose cap w = 1.")
    L += ["", "Full lists: `results/reach/weights.csv`.", "",
          md_table(sp[sp.space == "S2A_k20"], ["set", "k", "method", "over", "n", "spearman_rho", "p",
                                                "median_n_cells_selected"],
                   ["set", "k", "method", "over", "n", "Spearman rho (weight vs cells)", "p",
                    "median cells, selected"]), ""]
    fs = ft[ft.space == "S1"]
    fi = dg["full_space"]
    L += ["### Full gene space S1 (caveat)", "",
          "All 23,485 ruler coordinates, b = c_young - c_old, Gram-reduced form. S1 includes coordinates no "
          "factor can move (unmatched and missing-mask genes, z = 0 in every d_f) and many noise dimensions; "
          "the numbers are not comparable to S2A_k20 and are not gated. (D) is not run here. (B) uses "
          "sigma^2 from the same draws as full-space squared norms; "
          f"{dg['s_f_zero_B_full']} of 1,836 factors shrink to zero there.", "",
          md_table(pd.DataFrame(dict(k=list(K_ALL), A=[gv("A", k, "S1") for k in K_ALL],
                                     B=[gv("B", k, "S1") for k in K_ALL], C=[gv("C", k, "S1") for k in K_ALL])),
                   ["k", "A", "B", "C"], ["k", "A raw", "B shrunk", "C null bank"]), "",
          "Solver status in S1: " + "; ".join(f"{k}: converged {v['ok']}, Gram rank {v['rank']}"
                                              for k, v in fi.items()) + ".", ""]
    L += ["## What this does not cover", ""] + [f"- {x}" for x in WHAT_NOT] + [""]
    L += ["## Files", "",
          "`src/reach_run.py`; `results/reach/`: `PREREG.flag`, `stage0_reproduction.csv`, "
          "`stage0_norm_sample.csv`, `null_depths.csv`, `factors.csv`, `inputs_k20.npz`, `solver_check.json`, "
          "`frac_table.csv`, `random_targets.csv`, `weights.csv`, `spearman.csv`, `readings.json`, "
          "`diagnostics.json`, `run.log`, `run_stdout.txt`, `run_stderr.txt`."]
    return "\n".join(L) + "\n"


def write_findings(stop=None):
    FINDINGS.write_text(render(stop), encoding="utf-8")


def main():
    mode = sys.argv[1] if len(sys.argv) > 1 else "run"
    OUT.mkdir(parents=True, exist_ok=True)
    if mode == "findings":
        stop = json.loads((OUT / "STOP.json").read_text(encoding="utf-8")) if (OUT / "STOP.json").exists() else None
        write_findings(stop)
        return
    log = Log(OUT / ("run.log" if mode == "run" else "verify.log"))
    log("=" * 80)
    log(f"REACH {mode} seed={SEED} k={K_ALL} targets={N_TARGETS} workers={N_WORKERS}")
    try:
        if mode == "run":
            run(log)
        elif mode == "verify":
            ok = verify(log)
            sys.exit(0 if ok else 1)
        else:
            raise SystemExit(f"unknown mode {mode}")
    except Stop as e:
        log(f"[STOP] {e.key}: {e.message}")
        if mode == "run":
            dump("STOP.json", dict(key=e.key, message=e.message))
        sys.exit(2)


if __name__ == "__main__":
    main()
