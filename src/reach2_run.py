"""REACH2 — the REACH reachability bound rerun with a workable convergence rule.

Pre-registration: results/reach2/PREREG.flag (2026-10-17T22:59Z), written
before any REACH2 number. It carries REACH's flag (results/reach/PREREG.flag,
lines 15-217) verbatim and changes only three things: the convergence test is
relative (gap <= 1e-6 * max(P, 1)), a failed solve gets one retry with a
stricter solver setting, and every solve is recorded. This script reads the
flag and never rewrites it. It is "sizing information, never evidence".

REACH's Stage 0, readings and table helpers are imported unchanged from
src/reach_run.py; its solvers are copied here with only the convergence test,
the retry and the per-solve record changed. Nothing is written outside
results/reach2/ and FINDINGS_REACH2.md. results/reach/run.log is never opened.

Inputs: as REACH (results/south6/stage1.npz and the matched-gene CSR cache are
git-ignored; regenerated from the Southard download by src/south6_run.py).
The gated solves are re-verifiable from tracked files alone:
    python src/reach2_run.py verify

Usage, from the repository root:
    python src/reach2_run.py           Stage 0, inputs, solver check, sets A-D,
                                       readings, diagnostics, findings
    python src/reach2_run.py verify    re-solve sets A-D from inputs_k20.npz
    python src/reach2_run.py findings  rewrite FINDINGS_REACH2.md from outputs
    python src/reach2_run.py selftest DIR   the whole solve/readings/findings
                                       path on synthetic matrices, into DIR
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

import reach_run as R  # noqa: E402  (unchanged; defines functions and constants only)

ROOT = SRC.parent
OUT = Path(os.environ.get("REACH2_OUT", ROOT / "results" / "reach2"))
PREREG = ROOT / "results" / "reach2" / "PREREG.flag"
FINDINGS = Path(os.environ.get("REACH2_FINDINGS", ROOT / "FINDINGS_REACH2.md"))
REACH_OUT = ROOT / "results" / "reach"

SEED = R.SEED
KDIM = R.KDIM
KS = R.KS
K_ALL = R.K_ALL
N_TARGETS = R.N_TARGETS
N_FACTORS = R.N_FACTORS
TOL = R.TOL
TIE = R.TIE
NNZ_EPS = R.NNZ_EPS
N_BISECT = R.N_BISECT
LAM_LO = R.LAM_LO
N_RESTART = R.N_RESTART
EIG_REL = R.EIG_REL
N_WORKERS = R.N_WORKERS

GAP_REL = 1e-6           # O1: gap <= 1e-6 * max(P, 1)
# selftest-only switches (refused in run and verify): force first attempts to fail, or shrink the tolerance
_ST_FAIL_FIRST = os.environ.get("REACH2_SELFTEST_FAIL_FIRST") == "1"
if os.environ.get("REACH2_SELFTEST_GAP_REL"):
    GAP_REL = float(os.environ["REACH2_SELFTEST_GAP_REL"])
OLD_GAP_REL = 1e-8       # REACH: gap <= 1e-8 * max(1, 0.5 ||b||^2), counterfactual only
BVLS_BASE = dict(tol=1e-12, max_iter=10000)
BVLS_STRICT = dict(tol=1e-15, max_iter=100000)
LB_BASE = dict(maxiter=20000, ftol=1e-16, gtol=1e-12, maxcor=30, maxfun=15000, maxls=20)
LB_STRICT = dict(maxiter=200000, ftol=1e-20, gtol=1e-15, maxcor=30, maxfun=200000, maxls=20)

Stop = R.Stop
Log = R.Log


def dump(name, obj):
    (OUT / name).write_text(json.dumps(R._js(obj), indent=2), encoding="utf-8")


def tol_of(P):
    return GAP_REL * max(float(P), 1.0)


# --------------------------------------------------------------------------
# solvers: REACH's algorithms; O1 test, O2 retry, O3 record
# --------------------------------------------------------------------------
class SolveFail(Exception):
    """Second failure of one solve in a gated problem (O2)."""

    def __init__(self, call):
        super().__init__(f"{call}")
        self.call = call


class Problem:
    """min 0.5||A w - b||^2 (+ 0.5 c0) over w in [0, 1]^n; frac = (g - sqrt(.)) / g."""

    def __init__(self, A, b, g, c0=0.0):
        self.A = np.ascontiguousarray(A, np.float64)
        self.b = np.asarray(b, np.float64).ravel()
        self.g = float(g)
        self.c0 = float(c0)

    def frac(self, w):
        e = self.A @ w - self.b
        return (self.g - np.sqrt(max(float(e @ e) + self.c0, 0.0))) / self.g


class Rec:
    """Per-solve record for one problem (one space, set, target)."""

    def __init__(self, space, set_key, target, gated, g):
        self.base = dict(space=space, set=set_key, target=str(target))
        self.gated = bool(gated)
        self.old_tol = OLD_GAP_REL * max(1.0, 0.5 * float(g) ** 2)
        self.calls = []

    def add(self, **kw):
        row = dict(self.base)
        row.update(kw)
        row["ratio"] = row["gap"] / row["tol"] if row["tol"] > 0 else np.inf
        row["old_tol"] = self.old_tol
        self.calls.append(row)
        if not row["passed"] and self.gated:
            raise SolveFail(row)
        return row


def obj_gap(A, b, w, lam, c0):
    """P (with 0.5 c0, the constant a Gram reduction removes) and P - D."""
    y = A @ w - b
    yy = float(y @ y)
    P = 0.5 * yy + lam * float(w.sum())
    D = -0.5 * yy - float(b @ y) + float(np.minimum(0.0, A.T @ y + lam).sum())
    return P + 0.5 * c0, P - D


def _bvls_once(A, b, c0, setting):
    r = lsq_linear(A, b, bounds=(0.0, 1.0), method="bvls", **setting)
    w = np.clip(r.x, 0.0, 1.0)
    P, gp = obj_gap(A, b, w, 0.0, c0)
    t = tol_of(P)
    return w, dict(P=P, gap=gp, tol=t, status=int(r.status), passed=bool(r.status >= 1 and gp <= t))


def _score(info):
    return info["gap"] / info["tol"] if info["status"] >= 1 else np.inf


def bvls(A, b, c0, rec, **meta):
    n = A.shape[1]
    if n == 0:
        return np.zeros(0), True
    w, info = _bvls_once(A, b, c0, BVLS_BASE)
    if _ST_FAIL_FIRST:
        info["passed"] = False
    first = dict(info)
    retry = False
    if not info["passed"]:
        retry = True
        w2, info2 = _bvls_once(A, b, c0, BVLS_STRICT)
        if info2["passed"] or _score(info2) < _score(info):
            w, info = w2, info2
    rec.add(**meta, lam=0.0, P=info["P"], gap=info["gap"], tol=info["tol"], restarts=0,
            status=info["status"], retry_used=retry, gap_first=first["gap"], tol_first=first["tol"],
            status_first=first["status"], passed=info["passed"])
    return w, info["passed"]


def _clean(w):
    w = np.asarray(w, np.float64).copy()
    w[w <= NNZ_EPS] = 0.0
    return w


def solve_unlimited(P, rec):
    w, ok = bvls(P.A, P.b, P.c0, rec, method="bvls", k="unlimited", kind="bvls_unlimited", step=0)
    w = _clean(w)
    return dict(w=w, nnz=int((w > 0).sum()), frac=P.frac(w), ok=ok)


def _k_row_of_step(s):
    for k in KS:
        if s <= k:
            return str(k)
    return str(KS[-1])


def greedy(P, rec, kmax=max(KS)):
    A, b = P.A, P.b
    n = A.shape[1]
    colsq = np.einsum("ij,ij->j", A, A)
    usable = colsq > 0
    in_S = np.zeros(n, bool)
    S = []
    w = np.zeros(n)
    e = -b.copy()
    best = {k: dict(w=np.zeros(n), nnz=0, frac=0.0) for k in KS}
    ok_all, steps = True, 0
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
        steps += 1
        wS, ok = bvls(A[:, S], b, P.c0, rec, method="greedy", k=_k_row_of_step(steps), kind="greedy_refit",
                      step=steps)
        ok_all &= ok
        w = np.zeros(n)
        w[S] = wS
        w = _clean(w)
        e = A @ w - b
        nnz = int((w > 0).sum())
        fr = P.frac(w)
        for k in KS:
            if nnz <= k and fr > best[k]["frac"]:
                best[k] = dict(w=w.copy(), nnz=nnz, frac=fr)
        if nnz >= kmax or len(S) >= n or float(e @ e) <= 1e-24:
            break
    return best, ok_all


def _lbfgsb_once(A, b, lam, w0, c0, opts):
    n = A.shape[1]
    bounds = [(0.0, 1.0)] * n

    def f(w):
        e = A @ w - b
        return 0.5 * float(e @ e) + lam * float(w.sum()), A.T @ e + lam

    w = w0.copy()
    for i in range(N_RESTART + 1):
        res = minimize(f, w, jac=True, method="L-BFGS-B", bounds=bounds, options=opts)
        w = np.clip(res.x, 0.0, 1.0)
        P, gp = obj_gap(A, b, w, lam, c0)
        t = tol_of(P)
        if gp <= t:
            return w, dict(P=P, gap=gp, tol=t, restarts=i, status=int(res.status), passed=True)
    return w, dict(P=P, gap=gp, tol=t, restarts=N_RESTART, status=int(res.status), passed=False)


def l1_solve(A, b, lam, w0, c0, rec, **meta):
    w, info = _lbfgsb_once(A, b, lam, w0, c0, LB_BASE)
    if _ST_FAIL_FIRST:
        info["passed"] = False
    first = dict(info)
    retry = False
    if not info["passed"]:
        retry = True
        w2, info2 = _lbfgsb_once(A, b, lam, w0, c0, LB_STRICT)
        if info2["passed"] or info2["gap"] / info2["tol"] < info["gap"] / info["tol"]:
            w, info = w2, info2
    rec.add(**meta, lam=lam, P=info["P"], gap=info["gap"], tol=info["tol"], restarts=info["restarts"],
            status=info["status"], retry_used=retry, gap_first=first["gap"], tol_first=first["tol"],
            status_first=first["status"], passed=info["passed"])
    return w, info["passed"]


def l1(P, rec):
    A, b = P.A, P.b
    n = A.shape[1]
    c = A.T @ b
    lam_max = float(c.max())
    out = {}
    if lam_max <= 0:
        for k in KS:
            out[k] = dict(w=np.zeros(n), nnz=0, frac=0.0, nnz_hit=0, hit=False, lam=np.nan, ok=True)
        return out, True
    solved = {}
    ok_all = True
    lo, hi = np.log(lam_max * LAM_LO), np.log(lam_max)
    state = dict(k=None, i=0)

    def solve_at(t):
        nonlocal ok_all
        if t in solved:
            return solved[t]
        if solved:
            near = min(solved, key=lambda s: abs(s - t))
            w0 = solved[near]
        else:
            w0 = np.zeros(n)
        w, ok = l1_solve(A, b, float(np.exp(t)), w0, P.c0, rec, method="l1", k=str(state["k"]),
                         kind="lbfgsb", step=state["i"])
        ok_all &= ok
        solved[t] = w
        return w

    for k in KS:
        state["k"] = k
        a, z = lo, hi
        cand = None
        for i in range(N_BISECT):
            state["i"] = i + 1
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
        wS, ok = bvls(A[:, supp], b, P.c0, rec, method="l1", k=str(k), kind="l1_refit", step=0)
        ok_all &= ok
        w = np.zeros(n)
        w[supp] = wS
        w = _clean(w)
        out[k] = dict(w=w, nnz=int((w > 0).sum()), frac=P.frac(w), nnz_hit=nnz_hit,
                      hit=bool(nnz_hit == k), lam=float(np.exp(t)), ok=ok)
    return out, ok_all


def run_set(P, rec):
    """Every k for one problem, in REACH's order (greedy, L1, unlimited)."""
    gr, ok_g = greedy(P, rec)
    lr, ok_l = l1(P, rec)
    un = solve_unlimited(P, rec)
    res = {}
    for k in KS:
        G, L = gr[k], lr[k]
        res[str(k)] = dict(greedy=G, l1=L, gated=max(G["frac"], L["frac"]),
                           abs_diff=abs(G["frac"] - L["frac"]), jaccard=R.jaccard(G["w"], L["w"]))
    res["unlimited"] = dict(bvls=un, gated=un["frac"])
    return res, bool(ok_g and ok_l and un["ok"])


# --------------------------------------------------------------------------
# the gated pool: A, B, C young and D targets 0-199
# --------------------------------------------------------------------------
_W = {}


def _worker_init(mats, b20, targets, g):
    _W["mats"] = {k: np.asarray(v, np.float64) for k, v in mats.items()}
    _W["b20"] = np.asarray(b20, np.float64)
    _W["targets"] = np.asarray(targets, np.float64)
    _W["g"] = float(g)


def _d_row(t, res):
    row = dict(target=int(t))
    for k in KS:
        r = res[str(k)]
        row[f"greedy_k{k}"] = r["greedy"]["frac"]
        row[f"l1_k{k}"] = r["l1"]["frac"]
        row[f"l1_hit_k{k}"] = r["l1"]["hit"]
        row[f"gated_k{k}"] = r["gated"]
    row["gated_kunlimited"] = res["unlimited"]["gated"]
    row["nnz_unlimited"] = res["unlimited"]["bvls"]["nnz"]
    return row


def _worker(task):
    key, t = task
    if key == "D":
        A, b, target = _W["mats"]["B"], _W["targets"][t], int(t)
    else:
        A, b, target = _W["mats"][key], _W["b20"], "young"
    rec = Rec("S2A_k20", key, target, True, _W["g"])
    t0 = time.time()
    try:
        res, ok = run_set(Problem(A, b, _W["g"]), rec)
    except SolveFail as e:
        return dict(key=key, target=target, stop=True, call=e.call, calls=rec.calls, seconds=time.time() - t0)
    out = dict(key=key, target=target, stop=False, ok=ok, calls=rec.calls, seconds=time.time() - t0)
    out["res"] = res if key != "D" else None
    out["drow"] = _d_row(t, res) if key == "D" else None
    return out


class GatedStop(Stop):
    """STOP inside the gated pool, carrying what completed."""

    def __init__(self, key, message, calls, status):
        super().__init__(key, message)
        self.calls = calls
        self.status = status


def _pname(k, t):
    return ("S2A_k20", k, "young" if t < 0 else str(t))


def run_gated(mats, b20, targets, g, log):
    """Returns (res for A, B, C; D table; calls). Raises GatedStop on a second failure."""
    import multiprocessing as mp
    for v in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS"):
        os.environ[v] = "1"
    tasks = [("A", -1), ("B", -1), ("C", -1)] + [("D", t) for t in range(len(targets))]
    res, drows, calls, status = {}, [], [], {}
    t0 = time.time()
    pool = mp.get_context("spawn").Pool(N_WORKERS, initializer=_worker_init, initargs=(mats, b20, targets, g))
    try:
        for out in pool.imap_unordered(_worker, tasks):
            calls += out["calls"]
            name = _pname(out["key"], -1 if out["target"] == "young" else int(out["target"]))
            if out["stop"]:
                pool.terminate()
                status[name] = "stopped here"
                pending = [_pname(k, t) for k, t in tasks if _pname(k, t) not in status]
                for p in pending:
                    status[p] = "not completed (pool terminated)"
                c = out["call"]
                msg = (f"set {out['key']} target {out['target']}: second failure of one solve "
                       f"(method {c['method']}, k {c['k']}, kind {c['kind']}, step {c['step']}, lam {c['lam']:.6g}; "
                       f"gap {c['gap']:.3e} vs tolerance {c['tol']:.3e}, objective P {c['P']:.6g}, status "
                       f"{c['status']}; first attempt gap {c['gap_first']:.3e} vs {c['tol_first']:.3e}). "
                       f"Problems completed before it: {len(status) - 1 - len(pending)}; not completed when the "
                       f"pool was terminated: {len(pending)}.")
                raise GatedStop("solver_not_converged", msg, calls, status)
            status[name] = "completed"
            if out["key"] == "D":
                drows.append(out["drow"])
            else:
                res[out["key"]] = out["res"]
            n = len(status)
            if n % 20 == 0 or n == len(tasks):
                cc = pd.DataFrame(calls)
                log(f"[gated] {n}/{len(tasks)} problems ({time.time() - t0:.0f}s); solves {len(cc)}, "
                    f"retries {int(cc.retry_used.sum())}, max gap/tol {cc.ratio.max():.3g}")
    finally:
        pool.terminate()
        pool.join()
    dtab = pd.DataFrame(drows).sort_values("target").reset_index(drop=True)
    return res, dtab, calls


# --------------------------------------------------------------------------
# solver diagnostics tables (O3)
# --------------------------------------------------------------------------
CALL_COLS = ["space", "set", "target", "method", "k", "kind", "step", "lam", "P", "gap", "tol", "ratio",
             "restarts", "status", "retry_used", "gap_first", "tol_first", "status_first", "passed", "old_tol"]


def row_skeleton(space, set_key, target):
    rows = [dict(space=space, set=set_key, target=str(target), k=str(k), method=m)
            for k in KS for m in ("greedy", "l1")]
    rows.append(dict(space=space, set=set_key, target=str(target), k="unlimited", method="bvls"))
    return rows


def diag_rows(calls, problems, status):
    """One row per (space, set, target, k, method). problems: list of (space, set, target, methods);
    status: (space, set, target) -> completed / stopped here / not completed (pool terminated) / not run."""
    cc = pd.DataFrame(calls, columns=CALL_COLS)
    sk = []
    for space, set_key, target, only in problems:
        for r in row_skeleton(space, set_key, target):
            if only is None or r["method"] in only:
                r["problem_status"] = status.get((space, set_key, str(target)), "not run")
                sk.append(r)
    sk = pd.DataFrame(sk)
    if len(cc):
        cc = cc.copy()
        cc["ratio"] = cc["ratio"].astype(float)
        idx = cc.groupby(["space", "set", "target", "k", "method"])["gap"].idxmax()
        worst = cc.loc[idx, ["space", "set", "target", "k", "method", "gap", "P", "tol"]].rename(
            columns=dict(gap="final_gap", P="P_at_final_gap", tol="tol_at_final_gap"))
        agg = cc.groupby(["space", "set", "target", "k", "method"]).agg(
            n_solves=("gap", "size"), max_ratio=("ratio", "max"), n_retries=("retry_used", "sum"),
            all_passed=("passed", "all")).reset_index()
        agg = agg.merge(worst, on=["space", "set", "target", "k", "method"])
        out = sk.merge(agg, on=["space", "set", "target", "k", "method"], how="left")
    else:
        out = sk.assign(n_solves=np.nan, max_ratio=np.nan, n_retries=np.nan, all_passed=np.nan,
                        final_gap=np.nan, P_at_final_gap=np.nan, tol_at_final_gap=np.nan)
    out["n_solves"] = out["n_solves"].fillna(0).astype(int)
    out["n_retries"] = out["n_retries"].fillna(0).astype(int)
    out["retry_used"] = out["n_retries"] > 0
    # a row with solves: did all pass; a row without: nothing to certify if its problem completed
    # (greedy stopped early, L1 reused cached solves), not reached otherwise
    conv = []
    for r in out.itertuples():
        if r.n_solves > 0:
            conv.append(str(bool(r.all_passed)))
        else:
            conv.append("True" if r.problem_status == "completed" else "not reached")
    out["converged"] = conv
    return out[["space", "set", "target", "k", "method", "problem_status", "n_solves", "final_gap",
                "P_at_final_gap", "tol_at_final_gap", "max_ratio", "retry_used", "n_retries", "converged"]]


GATED_SETS = ("A", "B", "C", "D", "B+planted")


def write_solver_tables(calls, problems, status):
    cc = pd.DataFrame(calls, columns=CALL_COLS)
    cc.to_csv(OUT / "solver_calls.csv.gz", index=False, compression="gzip")
    dg = diag_rows(calls, problems, status)
    dg.to_csv(OUT / "solver_diagnostics.csv", index=False)
    gc = cc[cc.set.isin(GATED_SETS) & (cc.space == "S2A_k20")]
    w1 = gc.sort_values("gap", ascending=False).head(20).assign(ranked_by="gap")
    w2 = gc.sort_values("ratio", ascending=False).head(20).assign(ranked_by="gap/tol")
    pd.concat([w1, w2]).to_csv(OUT / "worst_gaps.csv", index=False)
    old_fail_new_pass = int((gc.passed & (gc.gap > gc.old_tol)).sum())
    is_bvls = gc.kind != "lbfgsb"
    first_fail = ~((gc.gap_first <= gc.tol_first) & (~is_bvls | (gc.status_first >= 1)))
    old_ok_first = (gc.gap_first <= gc.old_tol) & (~is_bvls | (gc.status_first >= 1))
    summ = dict(
        n_solves_gated=int(len(gc)), n_retries_gated=int(gc.retry_used.sum()),
        n_first_attempt_fail_gated=int(first_fail.sum()),
        n_second_fail_gated=int((~gc.passed).sum()),
        max_gap_gated=float(gc.gap.max()) if len(gc) else None,
        max_ratio_gated=float(gc.ratio.max()) if len(gc) else None,
        old_tol=float(gc.old_tol.iloc[0]) if len(gc) else None,
        pass_new_fail_old=old_fail_new_pass,
        fail_new_first_pass_old=int((first_fail & old_ok_first).sum()),
        n_solves_all=int(len(cc)), n_retries_all=int(cc.retry_used.sum()),
        n_second_fail_all=int((~cc.passed).sum()) if len(cc) else 0,
        by_kind={k: dict(n=int(len(x)), max_gap=float(x.gap.max()), max_ratio=float(x.ratio.max()),
                         retries=int(x.retry_used.sum()), median_P=float(x.P.median()),
                         max_restarts=int(x.restarts.max()))
                 for k, x in gc.groupby("kind")},
    )
    dump("solver_summary.json", summ)
    return dg, summ


# --------------------------------------------------------------------------
# Stage 0 additions: row-for-row against REACH, (h) against REACH's inputs
# --------------------------------------------------------------------------
NUM_COLS = ("stored", "reproduced", "max_abs_diff", "max_scaled_diff")
EXACT_COLS = ("item", "kind", "gated", "extra", "n", "n_mismatch", "ok")


def _num(x):
    try:
        v = float(x)
    except (TypeError, ValueError):
        return None
    return v


def compare_rowwise(mine, theirs):
    """mine/theirs: DataFrames read back from CSV (same dtypes as REACH's file)."""
    out = []
    n = max(len(mine), len(theirs))
    for i in range(n):
        a = mine.iloc[i] if i < len(mine) else None
        b = theirs.iloc[i] if i < len(theirs) else None
        if a is None or b is None:
            out.append(dict(row=i, item=(b if b is not None else a)["item"], column="(row)", reach="present" if b
                            is not None else "absent", reach2="present" if a is not None else "absent", ok=False))
            continue
        for c in EXACT_COLS:
            ok = str(a[c]) == str(b[c])
            out.append(dict(row=i, item=b["item"], column=c, reach=str(b[c]), reach2=str(a[c]), ok=ok))
        for c in NUM_COLS:
            rv, mv = b[c], a[c]
            rblank = pd.isna(rv) or str(rv) == ""
            mblank = pd.isna(mv) or str(mv) == ""
            if rblank:
                ok = mblank
            else:
                rn, mn = _num(rv), _num(mv)
                ok = (rn is not None and mn is not None and abs(mn - rn) <= TOL * max(1.0, abs(rn)))
            out.append(dict(row=i, item=b["item"], column=c, reach="" if rblank else str(rv),
                            reach2="" if mblank else str(mv), ok=bool(ok)))
        out.append(dict(row=i, item=b["item"], column="note (not gated)", reach=str(b["note"]),
                        reach2=str(a["note"]), ok=str(a["note"]) == str(b["note"])))
    df = pd.DataFrame(out)
    df["gated"] = df["column"] != "note (not gated)"
    return df


def stage0_h(rows, inp):
    """(h): this task's inputs against REACH's stored inputs (written before REACH's first solve)."""
    z = np.load(REACH_OUT / "inputs_k20.npz", allow_pickle=True)
    R.check(rows, "(h) inputs_k20.npz keep", np.array([str(x) for x in z["keep"]]), np.array(inp["keep"]),
            "exact", extra=True)
    R.check(rows, "(h) inputs_k20.npz n_cells", z["n_cells"], inp["n_cells"], "exact", extra=True)
    for key in ("D20", "b20", "g", "sigma2", "s_f", "C20", "targets"):
        R.check(rows, f"(h) inputs_k20.npz {key}", z[key], inp[key], extra=True)
    nd = pd.read_csv(REACH_OUT / "null_depths.csv")
    mine = inp["null_depths"]
    for c in nd.columns:
        kind = "exact" if c in ("n_cells", "R") else "tol"
        R.check(rows, f"(h) null_depths.csv {c}", nd[c].to_numpy(), mine[c].to_numpy(), kind, extra=True)
    fa = pd.read_csv(REACH_OUT / "factors.csv")
    mine = inp["factors"]
    for c in fa.columns:
        kind = "exact" if c in ("factor", "n_cells", "null_draw_index") else "tol"
        R.check(rows, f"(h) factors.csv {c}", fa[c].to_numpy(), mine[c].to_numpy(), kind, extra=True)


def gate(rows, log):
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


def build_inputs(ctx):
    """Exactly REACH's run() input construction."""
    keep, n_cells = ctx["keep"], ctx["n_cells"]
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
    null_depths = pd.DataFrame(dict(n_cells=ctx["ns"], sigma2_S2A_k20=[sigma2_n[n] for n in ctx["ns"]],
                                    sigma2_S1=[float(ctx["sqfull"][n].mean()) for n in ctx["ns"]], R=200))
    pct = pd.Series(n_cells).rank(pct=True, method="average").to_numpy() * 100.0
    factors = pd.DataFrame(dict(factor=keep, n_cells=n_cells, n_cells_pctile=pct, dnorm_S2A_k20=np.sqrt(dn2),
                                sigma2_S2A_k20=sigma2, s_f=s_f, null_draw_index=ctx["occ"],
                                null_dnorm_S2A_k20=np.sqrt(cn2), null_s=s_c))
    return dict(keep=keep, n_cells=n_cells, D20=D20, b20=b20, g=np.array([g]), sigma2=sigma2, s_f=s_f, C20=C20,
                targets=targets, B20=B20, Cs20=Cs20, s_c=s_c, pct=pct, null_depths=null_depths, factors=factors)


def write_inputs(inp, ctx):
    inp["null_depths"].to_csv(OUT / "null_depths.csv", index=False)
    inp["factors"].to_csv(OUT / "factors.csv", index=False)
    np.savez(OUT / "inputs_k20.npz", keep=np.array(inp["keep"], dtype=object), n_cells=inp["n_cells"],
             D20=inp["D20"], b20=inp["b20"], g=inp["g"], sigma2=inp["sigma2"], s_f=inp["s_f"], C20=inp["C20"],
             targets=inp["targets"])
    ctx["sample"].to_csv(OUT / "stage0_norm_sample.csv", index=False)


def solver_check(B20, b20, g, log):
    """Carried: (B) plus a planted column equal to b; greedy and L1 at k = 1."""
    Pp = Problem(np.hstack([B20.T, b20[:, None]]), b20, g)
    rec = Rec("S2A_k20", "B+planted", "planted", True, g)
    try:
        gr, okg = greedy(Pp, rec)
        lr, okl = l1(Pp, rec)
    except SolveFail as e:
        return None, rec.calls, e.call
    sel_g = np.flatnonzero(gr[1]["w"] > 0).tolist()
    sel_l = np.flatnonzero(lr[1]["w"] > 0).tolist()
    nplant = B20.shape[0]
    plant = dict(greedy_frac=gr[1]["frac"], greedy_selected=sel_g, l1_frac=lr[1]["frac"], l1_selected=sel_l,
                 ok_greedy=okg, ok_l1=okl, planted_index=nplant)
    plant["passed"] = bool(okg and okl and gr[1]["frac"] >= 1 - TIE and lr[1]["frac"] >= 1 - TIE
                           and sel_g == [nplant] and sel_l == [nplant])
    return plant, rec.calls, None


# --------------------------------------------------------------------------
# diagnostics (not gated): C shrunk, full gene space
# --------------------------------------------------------------------------
def gram_problem(Dm, b_full):
    """REACH's Gram reduction (reach_run.gram_problem) for this task's Problem."""
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
    return Problem(A, bp, np.sqrt(bb), c0=max(c0, 0.0)), dict(rank=int(keep.sum()), c0_raw=c0)


def run_diagnostics(ctx, inp, b20, g, log):
    calls, info = [], {}
    rec = Rec("S2A_k20", "Cs", "young", False, g)
    res_cs, ok_cs = run_set(Problem(inp["Cs20"].T, b20, g), rec)
    calls += rec.calls
    info["Cs"] = dict(ok=ok_cs)
    log(f"[diag] C shrunk: converged {ok_cs}")
    mr = ctx["mr"]
    b_full = ctx["c_young"] - ctx["c_old"]
    n_cells = inp["n_cells"]
    sigma2_full = np.array([float(ctx["sqfull"][int(n)].mean()) for n in n_cells])
    Dm = ctx["d_f"][:, mr]
    off = np.ones(R.N_RULER, bool)
    off[mr] = False
    b_m = b_full[mr]
    c_off = float(b_full[off] @ b_full[off])
    dn2_full = np.einsum("ij,ij->i", Dm, Dm)
    s_full = np.maximum(0.0, 1.0 - sigma2_full / dn2_full)
    res_full = {}
    for key, M in (("A", Dm), ("B", Dm * s_full[:, None]), ("C", ctx["C_full"])):
        t0 = time.time()
        P, gi = gram_problem(M, b_m)
        P.c0 += c_off
        P.g = float(np.linalg.norm(b_full))
        rec = Rec("S1", key, "young", False, P.g)
        r_, ok_ = run_set(P, rec)
        calls += rec.calls
        res_full[key] = r_
        info[f"S1_{key}"] = dict(ok=ok_, rank=gi["rank"], c0_raw_matched=gi["c0_raw"], c_off=c_off,
                                 seconds=time.time() - t0)
        log(f"[diag S1 {key}] converged {ok_} rank={gi['rank']} ({time.time() - t0:.0f}s)")
    info["s_f_zero_B_full"] = int((s_full == 0).sum())
    return res_cs, res_full, calls, info


# --------------------------------------------------------------------------
# main run
# --------------------------------------------------------------------------
def stop_out(e, calls, problems, status, log):
    write_solver_tables(calls, problems, status)
    dump("STOP.json", dict(key=e.key, message=e.message))
    log(f"[STOP] {e.key}: {e.message}")


def all_problems(n_targets):
    p = [("S2A_k20", "B+planted", "planted", ("greedy", "l1"))]
    p += [("S2A_k20", s, "young", None) for s in ("A", "B", "C")]
    p += [("S2A_k20", "D", t, None) for t in range(n_targets)]
    return p


def solve_and_read(inp, log, n_targets=N_TARGETS):
    """Solver check, then the gated pool. Shared by run and selftest."""
    b20 = inp["b20"]
    g = float(inp["g"][0])
    problems = all_problems(n_targets)
    status = {}
    plant, calls, fail = solver_check(inp["B20"], b20, g, log)
    if fail is not None:
        status[("S2A_k20", "B+planted", "planted")] = "stopped here"
        e = Stop("solver_not_converged", f"solver check: second failure of one solve ({fail}).")
        stop_out(e, calls, problems, status, log)
        raise e
    status[("S2A_k20", "B+planted", "planted")] = "completed"
    dump("solver_check.json", plant)
    log(f"[solver check] passed={plant['passed']} greedy_selected={plant['greedy_selected']} "
        f"l1_selected={plant['l1_selected']}")
    if not plant["passed"]:
        e = Stop("solver_not_converged", f"planted-column solver check failed: {plant}")
        stop_out(e, calls, problems, status, log)
        raise e
    mats = dict(A=inp["D20"].T, B=inp["B20"].T, C=inp["C20"].T)
    try:
        res, dtab, gcalls = run_gated(mats, b20, inp["targets"][:n_targets], g, log)
    except GatedStop as e:
        status.update(e.status)
        stop_out(e, calls + e.calls, problems, status, log)
        raise
    calls += gcalls
    for p in problems:
        status[(p[0], p[1], str(p[2]))] = "completed"
    return plant, res, dtab, calls, problems, status


def run(log):
    if not PREREG.exists():
        raise SystemExit(f"{PREREG} missing. Not computing.")
    rows = []
    ctx = R.stage0_static(log, rows)
    gate(rows, log)
    R.stage0_panel(ctx, log, rows)
    gate(rows, log)
    R.null_draws_all(ctx, log, rows)
    gate(rows, log)
    mine = pd.read_csv(OUT / "stage0_reproduction.csv")
    theirs = pd.read_csv(REACH_OUT / "stage0_reproduction.csv")
    cmp_ = compare_rowwise(mine, theirs)
    cmp_.to_csv(OUT / "stage0_vs_reach.csv", index=False)
    bad = cmp_[cmp_.gated & ~cmp_.ok]
    log(f"[stage0] row-for-row vs results/reach/stage0_reproduction.csv: {len(theirs)} REACH rows, "
        f"{len(mine)} REACH2 rows, {int(cmp_.gated.sum())} gated cells, {len(bad)} mismatches")
    if len(bad):
        raise Stop("prior_results_not_reproduced", "row-for-row against results/reach/stage0_reproduction.csv: "
                   + "; ".join(f"row {r.row} {r.item} [{r.column}] REACH {r.reach!r} REACH2 {r.reach2!r}"
                               for r in bad.itertuples()))
    inp = build_inputs(ctx)
    stage0_h(rows, inp)
    gate(rows, log)
    log("[stage0] all gated checks pass, (a)-(h) and row-for-row")
    write_inputs(inp, ctx)
    log(f"[inputs] written; s_f = 0 for {int((inp['s_f'] == 0).sum())} of {N_FACTORS} factors; "
        f"null bank s = 0 for {int((inp['s_c'] == 0).sum())}")

    plant, res, dtab, calls, problems, status = solve_and_read(inp, log)
    finish(inp, res, dtab, calls, problems, status, log, ctx=ctx)


def finish(inp, res, dtab, calls, problems, status, log, ctx=None):
    """All 203 gated problems converged: readings first, then diagnostics."""
    g = float(inp["g"][0])
    b20 = inp["b20"]
    dtab.to_csv(OUT / "random_targets.csv", index=False)
    rd = R.readings_of(res, dtab)
    rd["text"] = R.reading_text(rd)
    rd["g_S2A_k20"] = g
    rd["solver_check_passed"] = True
    dump("readings.json", rd)
    log(f"[reading] {rd['fired']} c_clause={rd['c_clause']} d_clause={rd['d_clause']}")

    keep, n_cells, pct = inp["keep"], inp["n_cells"], inp["pct"]
    ft, wr, sr = [], [], []
    for key in ("A", "B", "C"):
        ft += R.frac_rows("S2A_k20", key, res[key], True)
    for key in ("A", "B"):
        wr += R.weight_rows("S2A_k20", key, res[key], keep, n_cells, pct)
        sr += R.spearman_rows("S2A_k20", key, res[key], n_cells)
    info = {}
    if ctx is not None:
        res_cs, res_full, dcalls, info = run_diagnostics(ctx, inp, b20, g, log)
        calls = calls + dcalls
        ft += R.frac_rows("S2A_k20", "Cs", res_cs, False)
        for key, r_ in res_full.items():
            ft += [dict(r, space="S1") for r in R.frac_rows("S1", key, r_, False)]
            if key in ("A", "B"):
                sr += R.spearman_rows("S1", key, r_, n_cells)
        problems = problems + [("S2A_k20", "Cs", "young", None)] + [("S1", k, "young", None) for k in "ABC"]
        status = dict(status)
        for p in [("S2A_k20", "Cs", "young")] + [("S1", k, "young") for k in "ABC"]:
            status[p] = "completed"
    pd.DataFrame(ft).to_csv(OUT / "frac_table.csv", index=False)
    pd.DataFrame(wr).to_csv(OUT / "weights.csv", index=False)
    pd.DataFrame(sr).to_csv(OUT / "spearman.csv", index=False)
    dg, summ = write_solver_tables(calls, problems, status)
    dump("diagnostics.json", dict(
        diag_solver=info, s_f_zero_B=int((inp["s_f"] == 0).sum()), s_zero_null_bank=int((inp["s_c"] == 0).sum()),
        n_factors=int(len(keep)), s_f_zero_B_full=info.get("s_f_zero_B_full"), n_workers=N_WORKERS,
        diag_rows_not_converged=int(((dg.converged != "True") & ~dg.set.isin(list(GATED_SETS))).sum())))
    log(f"[solver] gated solves {summ['n_solves_gated']}, retries {summ['n_retries_gated']}, max gap "
        f"{summ['max_gap_gated']:.3e}, max gap/tol {summ['max_ratio_gated']:.3g}")
    log("[run] DONE")


# --------------------------------------------------------------------------
# verify (tracked files only)
# --------------------------------------------------------------------------
def verify(log):
    z = np.load(OUT / "inputs_k20.npz", allow_pickle=True)
    D20, b20, g = z["D20"], z["b20"], float(z["g"][0])
    B20 = D20 * z["s_f"][:, None]
    mats = dict(A=D20.T, B=B20.T, C=z["C20"].T)
    res, dtab, calls = run_gated(mats, b20, z["targets"], g, log)
    rd = R.readings_of(res, dtab)
    st = json.loads((OUT / "readings.json").read_text(encoding="utf-8"))
    ft = pd.read_csv(OUT / "frac_table.csv")
    bad = []
    for key in ("A", "B", "C"):
        for k in K_ALL:
            m = "bvls" if k == "unlimited" else "gated"
            sv = ft[(ft.space == "S2A_k20") & (ft.set == key) & (ft.k.astype(str) == k) & (ft.method == m)]["frac"]
            sv = float(sv.iloc[0])
            if abs(sv - res[key][k]["gated"]) > TOL * max(1, abs(sv)):
                bad.append(f"{key} k={k}: stored {sv} recomputed {res[key][k]['gated']}")
    dst = pd.read_csv(OUT / "random_targets.csv")
    maxd = 0.0
    for k in K_ALL:
        d = np.abs(dst[f"gated_k{k}"].to_numpy() - dtab[f"gated_k{k}"].to_numpy())
        maxd = max(maxd, float(d.max()))
        if np.any(d > TOL):
            bad.append(f"D k={k}: {int((d > TOL).sum())} targets differ")
    if rd["fired"] != st["fired"]:
        bad.append(f"reading stored {st['fired']} recomputed {rd['fired']}")
    cc = pd.DataFrame(calls)
    out = dict(fired=rd["fired"], stored=st["fired"], mismatches=bad, max_abs_diff_D=maxd,
               solves=int(len(cc)), retries=int(cc.retry_used.sum()), max_ratio=float(cc.ratio.max()))
    dump("verify.json", out)
    log(f"[verify] reading {rd['fired']} (stored {st['fired']}); mismatches: {bad if bad else 'none'}")
    return not bad


# --------------------------------------------------------------------------
# findings
# --------------------------------------------------------------------------
fmt = R.fmt
md_table = R.md_table
SET_LABEL = dict(R.SET_LABEL, D="D random targets", **{"B+planted": "B + planted column (solver check)"})

CARRIED = [
    ("THE OPTIMISATION: origin c_old, target c_young, b = (c_young - c_old) @ L20, g = ||b||, w in [0, 1]^1836 "
     "(dose cap), frac_closed(w) = (g - ||D20^T w - b||) / g; sparsity caps k = 1, 3, 5, 10, 50 and unlimited; "
     "unlimited by BVLS; greedy and L1 (L-BFGS-B, log-space bisection on lam, BVLS refit of the support) both "
     "run at finite k; gated value max(greedy, L1-refit)", "results/reach/PREREG.flag lines 54-81 (flag text, copied)"),
    ("The duality-gap definition (P, D, lam = 0 for BVLS, restricted to the active columns for greedy re-fits), "
     "BVLS status >= 1, and up to 20 L-BFGS-B restarts from the last point",
     "results/reach/PREREG.flag lines 82-87; only the threshold is replaced (O1)"),
    ("INSTRUMENTS: d_f and keep from results/south6/stage1.npz, L20 = south7 basis L_mid_full[:, :20], c_old and "
     "c_young from south7 basis.npz, projection of differences only, cell counts from cell_indptr, the "
     "non-targeting pool and the panel TMM cache rebuilt from stored counts", "results/reach/PREREG.flag lines 36-52"),
    ("Set (A) raw D20", "results/reach/PREREG.flag line 90"),
    ("Set (B) shrinkage s_f = max(0, 1 - sigma_f^2 / ||D20_f||^2), sigma^2(n) the mean of ||E L20||^2 over the 200 "
     "depth-matched draws of south7_run.null_draws(n, 14013), exact per n", "results/reach/PREREG.flag lines 91-101"),
    ("Set (C) null factor bank: the depth-matched construction (SOUTH6REC primary scheme, as SOUTH7), draw index "
     "m_i, not shrunk", "results/reach/PREREG.flag lines 102-113"),
    ("Set (D) 200 random targets u_t * g from default_rng([20260914, 4]).standard_normal((200, 20)), against (B)'s "
     "vectors, every k, both methods", "results/reach/PREREG.flag lines 114-117"),
    ("Stage 0 (a)-(g), tolerance 1e-4 * max(1, |stored|), Stage 0 sample default_rng([20260914, 0])",
     "results/reach/PREREG.flag lines 119-146"),
    ("The planted-column solver check (B's matrix plus a column equal to b; greedy and L1 at k = 1 must reach "
     "frac >= 1 - 1e-6 and select that column)", "results/reach/PREREG.flag lines 148-153"),
    ("Readings R0 (UNINTERPRETABLE) and R1 (BOUND REPORTED), their wording, their evaluation order, the tie "
     "tolerance 1e-6, \"at any k\" on both R0 clauses, numpy default percentiles",
     "results/reach/PREREG.flag lines 155-184"),
    ("The sentence \"No reading of this task licenses a claim about aged cells. A positive bound sizes an "
     "experiment; it does not support one.\"", "results/reach/PREREG.flag lines 34 and 186"),
    ("The R1 slot when B sits above the D range at some k and below at others",
     "src/reach_run.py reading_text (REACH's code, not its flag); stated in the REACH2 flag"),
    ("DIAGNOSTICS: (A) beside (B); weights and cell counts with Spearman; s_f zero count with the null-bank "
     "baseline; (C) shrunk; the full gene space S1 in Gram-reduced form (D not run there; failures recorded, "
     "not a STOP)", "results/reach/PREREG.flag lines 188-206"),
    ("STOP keys prior_results_not_reproduced and basis_not_orthonormal; \"no computation after a STOP, no "
     "restart with an adjustment\"", "results/reach/PREREG.flag lines 208-214"),
    ("WHAT THIS DOES NOT COVER (all six items)", "results/reach/PREREG.flag lines 216-228"),
    ("Seed 20260914 and every generator named", "results/reach/PREREG.flag line 52 and the set definitions"),
    ("L-BFGS-B maxfun 15000 and maxls 20 (scipy defaults)", "src/reach_run.py l1_solve (REACH's flag did not name "
     "them)"),
    ("16 single-threaded spawned workers for the pool; weights at or below 1e-6 zeroed after every BVLS fit",
     "src/reach_run.py (N_WORKERS, _clean)"),
]

DEVIATIONS = [
    "**Only the three pre-registered changes alter the computation** (O1 relative convergence test, O2 one "
    "retry with a stricter setting, O3 the per-solve record). Everything else is REACH's, listed under "
    "\"Carried over\".",
    "**\"Its own objective value\" fixed as the primal objective P** of that solve at the w it returns; in the "
    "S1 Gram-reduced diagnostic P includes the constant 0.5 c0 the reduction removes. Fixed in the flag; the task "
    "text did not say primal or dual.",
    "**\"A solve\" fixed as one solver call** (one BVLS call, or one L-BFGS-B solve at one lam with its up to 20 "
    "restarts); the retry re-solves that call from the same starting point. The per-row table (set, target, k, "
    "method) aggregates the calls of a row; the rule assigning calls to rows is in the flag (O3).",
    "**Strict solver setting chosen by this session** (the task asked only for a tighter internal tolerance and "
    "more iterations): BVLS tol 1e-15, max_iter 100000; L-BFGS-B ftol 1e-20, gtol 1e-15, maxiter 200000, "
    "maxfun 200000. Not tuned on any data.",
    "**Added beyond the task text:** Stage 0 (h), this task's inputs against REACH's stored `inputs_k20.npz`, "
    "`null_depths.csv` and `factors.csv` (gated); `stage0_vs_reach.csv`; `solver_calls.csv.gz`, "
    "`worst_gaps.csv`, `solver_summary.json`; the old-rule counterfactual count (not gated); `selftest` mode; "
    "`verify` mode as in REACH.",
    "**Files of REACH's directory read by this task's code:** `stage0_reproduction.csv` (row-for-row), "
    "`inputs_k20.npz`, `null_depths.csv`, `factors.csv` (Stage 0 (h)). Per `src/reach_run.py` these were all "
    "written before REACH's first solve and hold no solve output. Not opened by this session or its code: "
    "`run.log`, `run_stdout.txt`, `run_stderr.txt`, `STOP.json`, `solver_check.json`.",
    "**Execution changed, not the computation:** sets A, B, C now run in the same pool of 16 spawned "
    "single-threaded workers as set D (REACH ran them in the main process with multi-threaded BLAS), results are "
    "taken as they complete, and no frac value is logged before all 203 gated problems have converged. BLAS "
    "thread count can change floating-point summation order, not results beyond ~1e-12.",
    "**Imported unchanged from `src/reach_run.py`:** `stage0_static`, `stage0_panel`, `null_draws_all`, `check`, "
    "`readings_of`, `reading_text`, `frac_rows`, `weight_rows`, `spearman_rows`, `jaccard`, `fmt`, `md_table`, "
    "`Stop`, `Log`, `_js` and its constants. Through them REACH's imports stand (`south6_run`, "
    "`south6rec_run.indicator`, `south7_run.null_draws`, `seng_run`, `toward_run`, `fibro2_common`); importing "
    "`south7_run`, `seng_run`, `toward_run` and `south6_run` runs their `mkdir` on directories that already exist; "
    "nothing is written there. REACH's `gate_stage0` is not called (it writes into `results/reach/`); this task "
    "has its own gate. REACH's solvers and Gram reduction are copied into `src/reach2_run.py` with only O1-O3 "
    "changed.",
    "**Where written:** `results/reach2/`, `src/reach2_run.py`, `FINDINGS_REACH2.md`. Nothing in "
    "`results/reach/` or any other earlier results directory. `results/transfer/NOT_A_RESULT_gate_failed/` was "
    "not opened.",
    "**The first flag-writing attempt failed** (a shell here-document was cut off before the flag file was "
    "created; no file was written and nothing had been computed). The flag was then assembled by a script that "
    "copies REACH's lines 15-217 byte for byte and asserts the copy.",
    "**Synthetic tests before the real run, no project data.** `selftest` (added) runs the solver check, pool, "
    "readings and findings on random matrices (120 vectors, 24 targets) into the session scratchpad. Four runs: "
    "clean (2,357 solves, no retry); every first attempt forced to fail (all 2,357 retried, all passed); "
    "tolerance 1e-30 (STOP at the solver check); tolerance 1e-11 (STOP inside the pool, diagnostics table "
    "written). The forcing switches are environment variables (`REACH2_SELFTEST_FAIL_FIRST`, "
    "`REACH2_SELFTEST_GAP_REL`, `REACH2_OUT`, `REACH2_FINDINGS`); `run`, `verify` and `findings` refuse to start "
    "when any is set. Seen before the real run and recorded here: in the 1e-11 test the strict L-BFGS-B retry "
    "ended at exactly the gap of the first attempt, both ending in a line-search stop (scipy status 2). Where "
    "L-BFGS-B stops on its line search, a tighter tolerance and more iterations do not help, so O2's retry "
    "mainly helps BVLS and iteration-limited solves. The flag was not changed.",
    "**A stray file, removed.** Testing that `verify` refuses an empty directory created "
    "`results/reach2/verify.log` with only its header line (no number). It was deleted before the real run.",
    "**Script edited after the run, renderer only.** The run and the verify re-solve used `src/reach2_run.py` "
    "with SHA-256 a46a5d2b082bee8b27e3a42f3661ddee9fd20ffa88af6120f98d533f387ad154 (hashed at launch, "
    "unchanged at the end of the run). After verify finished, two edits were made: the hand-written notes "
    "(`NOTES`, labelled as notes in the findings) and string formatting of the counts in the set-D greedy/L1 "
    "table. No solver, Stage 0, reading or output-writing code was changed. The findings were then rewritten "
    "by `findings` mode, which computes nothing. The committed file is the edited one.",
    "**Console output kept outside the repository.** The run's stdout (the same lines as `run.log`) and stderr "
    "went to the session scratchpad and are not committed. stderr was empty. "
    "`solver_calls.csv.gz` is 553 KB, under the 1 MB line, so it is tracked and no `.gitignore` line was "
    "added.",
]

PRIOR_PARTIAL = (
    "**Prior partial output, disclosed.** REACH solved sets A, B and C against the young target before its STOP "
    "and printed those values to `results/reach/run.log`. They are not a result. This session did not open "
    "`run.log` or `run_stdout.txt` (the same log stream), and this task's code does not read them. No number from "
    "REACH's partial solve reached this session by any route. The REACH numbers this session saw are its Stage 0 "
    "table, the STOP's gap (4.268e-06) and tolerance (3.900e-06), and g = 27.929154.")

WHAT_NOT = list(R.WHAT_NOT) + [
    "The relative convergence test certifies each solve against its own optimum to the stated accuracy. It does "
    "not make greedy or L1 optimal for k > 1: both remain lower estimates of the k-sparse optimum, as REACH said.",
]

# Written by hand after the run and labelled as notes, not readings. Numbers are from the outputs named.
NOTES = [
    "Both R0 clauses fired at every k, and either alone fires R0. The D clause does not depend on the flag's "
    "choice to compare unshrunk (C) with shrunk (B). At unlimited k, B (0.174) sits next to the D median "
    "(0.168), inside a middle 90% of 0.030-0.395. The C clause, by contrast, rests on that choice: the "
    "diagnostic (C) shrunk by B's rule closes less than B at every k (0.096 at unlimited, "
    "against 0.174); that comparison is not gated and does not change the reading.",
    "R0 is not a finding that young is unreachable. It says this bound carries no information either way. The "
    "planted-column check shows the optimiser finds a direction that is present, so the bound could register "
    "such a direction. A positive R1 would still have been sizing information only.",
    "For sets A, B and C the unlimited optimum is sparse: 9 factors for A, 10 for B and 6 for C "
    "(`frac_table.csv`). So k = 10, 50 and unlimited give the same value. In 20 dimensions, a handful of "
    "vectors is enough to reach the best point of the cone the dose cap allows.",
    "Shortlist against cell count (this project's known trap). Nine of the 10 factors in B's unlimited support "
    "sit above the 75th cell-count percentile (chance about 2.5 of 10); CDX2 at p21 is the exception. The "
    "support's median cell count is 158 against 106 overall, and Spearman rho between weight and cell count "
    "over all 1,836 factors is 0.080 (p = 0.0006). Part of this is mechanical: shrinkage cuts low-count "
    "factors hardest. Under R0 no factor list here is a shortlist.",
    "What the convergence change did. 531 gated solves end with a gap above REACH's 3.900e-06 and pass O1 "
    "(largest gap/tolerance 0.199). All 531 are intermediate L-BFGS-B solves of the L1 bisection: 526 in "
    "set D (122 of the 200 targets), 3 in B, 2 in A. Under REACH's rule each would have been restarted. They "
    "decide which support the bisection finds; they do not decide the reported frac. Every reported frac comes "
    "from a BVLS fit (greedy re-fit, L1 refit or unlimited), and the largest BVLS gap relative to max(P, 1) "
    "is 9.6e-15. No solve needed a retry. Three L-BFGS-B solves used a restart. On set D target 11, where "
    "REACH stopped, the largest gap is 1.01e-05: above REACH's tolerance, and 0.045 of that solve's O1 "
    "tolerance.",
    "Because REACH's partial values for sets A, B and C were never seen, this run's A-C values were not "
    "compared with them, deliberately. This run is not a check on REACH's partial solve.",
]


def _prereg_block():
    return "```\n" + PREREG.read_text(encoding="utf-8").rstrip("\n") + "\n```"


def _stage0_section():
    rep = pd.read_csv(OUT / "stage0_reproduction.csv")
    L = ["## Stage 0 — reproduction of stored numbers", "",
         "Counts exact; numerics within 1e-4 x max(1, |stored|). Every row is gated; rows marked `extra` were "
         "added beyond the task text (REACH's or this task's). Rows (a)-(g) are REACH's 30 checks re-run by this "
         "task's script; rows (h) are added by REACH2. File: `results/reach2/stage0_reproduction.csv`.", ""]
    rep = rep.copy()
    rep["extra"] = rep["extra"].map(lambda x: "extra" if str(x) == "True" else "")
    L.append(md_table(rep.fillna(""), ["item", "kind", "extra", "n", "stored", "reproduced", "max_scaled_diff",
                                        "n_mismatch", "ok"]))
    n_fail = int(((rep.ok.astype(str) != "True") & (rep.gated.astype(str) == "True")).sum())
    L += ["", f"{len(rep)} checks, {int((rep.gated.astype(str) == 'True').sum())} gated, {n_fail} gated failures."]
    p = OUT / "stage0_vs_reach.csv"
    if p.exists():
        c = pd.read_csv(p)
        g_ = c[c.gated.astype(str) == "True"]
        nb = int((g_.ok.astype(str) != "True").sum())
        ng = c[c.gated.astype(str) != "True"]
        L += ["", "### Row for row against `results/reach/stage0_reproduction.csv`", "",
              f"{c.row.nunique()} rows compared, {len(g_)} gated cells (item, kind, gated, extra, n, n_mismatch, ok "
              f"identical; stored, reproduced, max_abs_diff, max_scaled_diff within tolerance or blank where "
              f"REACH's is blank): **{nb} mismatches**. The `note` column (not gated) is identical in "
              f"{int((ng.ok.astype(str) == 'True').sum())} of {len(ng)} rows. File: "
              "`results/reach2/stage0_vs_reach.csv`."]
        if nb:
            L += ["", md_table(g_[g_.ok.astype(str) != "True"], ["row", "item", "column", "reach", "reach2"])]
    return L


def _solver_section(summ, dg):
    L = ["## Solver diagnostics (O3)", "",
         "Every solve was recorded whatever the outcome. A solve is one solver call; a row of "
         "`solver_diagnostics.csv` is one (space, set, target, k, method) and aggregates the solves assigned to it "
         "by the flag's rule. `final gap` is the largest gap any solve of the row ended with (after any retry); "
         "`gap/tol` is gap over that solve's own O1 tolerance 1e-6 x max(P, 1).", ""]
    if summ:
        L += [f"- Gated solves (solver check, A, B, C young, D targets): {summ['n_solves_gated']}. Failed O1 at "
              f"the first attempt: {summ['n_first_attempt_fail_gated']}. Retries used: {summ['n_retries_gated']}. "
              f"Second failures: {summ['n_second_fail_gated']}.",
              f"- Largest gated gap: {summ['max_gap_gated']:.3e}; largest gap/tolerance ratio: "
              f"{summ['max_ratio_gated']:.3g}.",
              f"- Counterfactual against REACH's rule (not gated; evaluated on this run's gaps, not a replay of "
              f"REACH, whose restart loop stopped at different points): {summ['pass_new_fail_old']} gated solves "
              f"pass O1 but have a gap above REACH's {summ['old_tol']:.3e}; {summ['fail_new_first_pass_old']} "
              f"failed O1 at the "
              f"first attempt but would have passed REACH's test.", ""]
        bk = pd.DataFrame([dict(kind=k, **v) for k, v in summ["by_kind"].items()])
        if len(bk):
            L += ["Gated solves by kind:", "",
                  md_table(bk, ["kind", "n", "max_gap", "max_ratio", "retries", "median_P", "max_restarts"],
                           ["kind", "solves", "max gap", "max gap/tol", "retries", "median objective P",
                            "max L-BFGS-B restarts"]), ""]
    wg = OUT / "worst_gaps.csv"
    if wg.exists():
        w = pd.read_csv(wg)
        cols = ["set", "target", "k", "method", "kind", "step", "lam", "P", "gap", "tol", "ratio", "retry_used",
                "passed"]
        for by in ("gap", "gap/tol"):
            x = w[w.ranked_by == by]
            L += [f"The 20 largest gated gaps, ranked by {by}:", "", md_table(x, cols), ""]
    if dg is not None and len(dg):
        x = dg[dg.space == "S2A_k20"].copy()
        x["setk"] = x["set"]
        agg = x.groupby(["set", "k", "method"], sort=False).agg(
            rows=("target", "size"), solves=("n_solves", "sum"), max_final_gap=("final_gap", "max"),
            max_ratio=("max_ratio", "max"), retries=("n_retries", "sum"),
            not_converged=("converged", lambda s: int((s.astype(str) == "False").sum())),
            not_reached=("converged", lambda s: int((s.astype(str) == "not reached").sum()))).reset_index()
        L += ["Per (set, k, method), over targets (set D has 200 rows per k and method; full per-target table: "
              "`results/reach2/solver_diagnostics.csv`):", "",
              md_table(agg, ["set", "k", "method", "rows", "solves", "max_final_gap", "max_ratio", "retries",
                             "not_converged", "not_reached"],
                       ["set", "k", "method", "rows", "solves", "max final gap", "max gap/tol", "retries",
                        "rows not converged", "rows not reached"]), ""]
        s1 = dg[dg.space == "S1"]
        if len(s1):
            a1 = s1.groupby(["set"]).agg(solves=("n_solves", "sum"), max_final_gap=("final_gap", "max"),
                                         max_ratio=("max_ratio", "max"), retries=("n_retries", "sum"),
                                         not_converged=("converged", lambda s: int((s.astype(str) == "False").sum()))
                                         ).reset_index()
            L += ["Full gene space S1 (diagnostic, not gated):", "",
                  md_table(a1, ["set", "solves", "max_final_gap", "max_ratio", "retries", "not_converged"]), ""]
    return L


def render(stop=None):
    rd = json.loads((OUT / "readings.json").read_text(encoding="utf-8")) if (OUT / "readings.json").exists() else None
    summ = json.loads((OUT / "solver_summary.json").read_text(encoding="utf-8")) \
        if (OUT / "solver_summary.json").exists() else None
    dg = pd.read_csv(OUT / "solver_diagnostics.csv", dtype={"target": str}) \
        if (OUT / "solver_diagnostics.csv").exists() else None
    L = ["# FINDINGS_REACH2 — the reachability bound, rerun with a relative convergence rule", ""]
    if stop:
        L += [f"**Status: STOP `{stop['key']}`.** {stop['message']}", "",
              "Per CLAUDE.md, nothing was computed after the STOP and the run was not restarted. No frac value of "
              "any set is written.", ""]
    else:
        L += [f"**Status: DONE.** Reading: `{'UNINTERPRETABLE' if rd['fired'] == 'R0' else 'BOUND REPORTED'}` "
              f"({rd['fired']}).", "",
              "This is sizing information, never evidence (TF search model blueprint). No reading of this task "
              "licenses a claim about aged cells. A positive bound sizes an experiment; it does not support one.",
              ""]
    L += [PRIOR_PARTIAL, "",
          "Reproduced by `python src/reach2_run.py` then `python src/reach2_run.py findings`. "
          "`python src/reach2_run.py verify` re-solves sets A-D from tracked `results/reach2/inputs_k20.npz` "
          "under the same convergence rule.", "",
          "## Pre-registration (verbatim from `results/reach2/PREREG.flag`)", "", _prereg_block(), ""]
    L += ["## Carried over from REACH, item by item", "",
          md_table(pd.DataFrame(CARRIED, columns=["item", "source"]), ["item", "source"]), ""]
    L += ["## Deviations and decisions, disclosed", ""] + [f"{i + 1}. {d}" for i, d in enumerate(DEVIATIONS)] + [""]
    if (OUT / "stage0_reproduction.csv").exists():
        L += _stage0_section() + [""]
    sc = json.loads((OUT / "solver_check.json").read_text(encoding="utf-8")) \
        if (OUT / "solver_check.json").exists() else None
    if sc:
        L += ["## Solver check (carried; the planted problem with a known answer)", "",
              f"Set (B)'s 1,836 vectors plus one planted column equal to the gap b (index {sc['planted_index']}). "
              f"Greedy k=1: frac {fmt(sc['greedy_frac'], 6)}, selected {sc['greedy_selected']}; L1 k=1: frac "
              f"{fmt(sc['l1_frac'], 6)}, selected {sc['l1_selected']}. Passed: {sc['passed']}. This shows the "
              "optimiser finds a direction that is present; the bound is an optimisation, not a test, so it has no "
              "recovery power in the usual sense. Sets (C) and (D) are its comparisons.", ""]
    if stop or rd is None:
        L += _solver_section(summ, dg)
        L += ["## What this does not cover", ""] + [f"- {w}" for w in WHAT_NOT] + [""]
        present = sorted(p.name for p in OUT.iterdir() if p.is_file())
        L += ["## Files", "", "`src/reach2_run.py`; `results/reach2/`: " + ", ".join(f"`{p}`" for p in present) + "."]
        return "\n".join(L) + "\n"

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
        x = ft[(ft.space == space) & (ft.set == s) & (ft.k.astype(str) == k) & (ft.method == m)]
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
          "", "Files: `results/reach2/frac_table.csv` (every set, k and method), "
          "`results/reach2/random_targets.csv` (every target, k and method).", ""]
    for n in NOTES:
        L += [f"**Note (not a pre-registered reading):** {n}", ""]
    L += _solver_section(summ, dg)

    rows = []
    for s in ("A", "B", "C", "Cs"):
        for k in KS:
            g_ = k20[(k20.set == s) & (k20.k.astype(str) == str(k))]
            if not len(g_):
                continue
            G = g_[g_.method == "greedy"].iloc[0]
            Lr = g_[g_.method == "l1"].iloc[0]
            Q = g_[g_.method == "gated"].iloc[0]
            rows.append(dict(set=SET_LABEL[s], k=k, greedy=G.frac, greedy_nnz=int(G.nnz), l1=Lr.frac,
                             l1_nnz=int(Lr.nnz),
                             l1_hit=f"{int(Lr.l1_nnz_hit)}{'' if str(Lr.l1_hit_exact) == 'True' else ' (not exact)'}",
                             diff=Q.abs_diff_greedy_l1, jac=Q.jaccard))
    mt = pd.DataFrame(rows)
    dis = []
    for k in KS:
        d = np.abs(dt[f"greedy_k{k}"] - dt[f"l1_k{k}"])
        # counts as strings so md_table does not print them as floats
        dis.append(dict(k=str(k), med=float(d.median()), mx=float(d.max()),
                        l1_better=str(int((dt[f"l1_k{k}"] > dt[f"greedy_k{k}"] + TIE).sum())),
                        gr_better=str(int((dt[f"greedy_k{k}"] > dt[f"l1_k{k}"] + TIE).sum())),
                        miss=str(int((dt[f"l1_hit_k{k}"].astype(str) != "True").sum()))))
    L += ["## Greedy versus L1 at finite k", "",
          "Both methods are feasible solutions; neither is guaranteed optimal for k > 1 (greedy is exhaustive "
          "at k = 1). `L1 nnz hit` is the support size the bisection reached before the BVLS refit; "
          "\"not exact\" means no lam in 60 steps gave exactly k and the fallback rule was used.", "",
          md_table(mt, ["set", "k", "greedy", "greedy_nnz", "l1", "l1_nnz", "l1_hit", "diff", "jac"],
                   ["set", "k", "greedy frac", "greedy nnz", "L1 frac", "L1 nnz", "L1 nnz hit",
                    "abs diff", "support Jaccard"]), "",
          "Set (D), over the 200 targets:", "",
          md_table(pd.DataFrame(dis), ["k", "med", "mx", "l1_better", "gr_better", "miss"],
                   ["k", "median abs diff", "max abs diff", "L1 better (targets)", "greedy better (targets)",
                    "L1 did not hit k exactly"]), ""]

    dgj = json.loads((OUT / "diagnostics.json").read_text(encoding="utf-8"))
    fac = pd.read_csv(OUT / "factors.csv")
    w = pd.read_csv(OUT / "weights.csv")
    sp = pd.read_csv(OUT / "spearman.csv")
    L += ["## Diagnostics (not gated)", "",
          "### What shrinkage costs: (A) beside (B)", "",
          md_table(tab, ["k", "A", "B"], ["k", "A raw", "B shrunk"]), "",
          f"Factors shrunk to zero in (B) (||D20_f||^2 <= sigma^2(n_f)): {dgj['s_f_zero_B']} of {dgj['n_factors']}. "
          f"Chance baseline: the same rule shrinks {dgj['s_zero_null_bank']} of the {dgj['n_factors']} null-bank "
          "vectors (C) to zero. Median s_f "
          f"{fmt(float(fac.s_f.median()))}; median ||D20_f|| {fmt(float(fac.dnorm_S2A_k20.median()))} against "
          f"median sigma(n_f) {fmt(float(np.sqrt(fac.sigma2_S2A_k20).median()))}.", ""]
    if "Cs" in set(ft.set):
        L += ["### (C) shrunk by the same rule", "",
              md_table(pd.DataFrame(dict(k=list(K_ALL), Cs=[gv("Cs", k) for k in K_ALL],
                                         C=[gv("C", k) for k in K_ALL], B=[gv("B", k) for k in K_ALL])),
                       ["k", "Cs", "C", "B"], ["k", "C shrunk (diag)", "C null bank", "B shrunk"]), ""]
    L += ["### Which factors carry the weight, and cell count", "",
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
    L += ["", "Full lists: `results/reach2/weights.csv`.", "",
          md_table(sp[sp.space == "S2A_k20"], ["set", "k", "method", "over", "n", "spearman_rho", "p",
                                                "median_n_cells_selected"],
                   ["set", "k", "method", "over", "n", "Spearman rho (weight vs cells)", "p",
                    "median cells, selected"]), ""]
    if "S1" in set(ft.space):
        fi = dgj["diag_solver"]
        L += ["### Full gene space S1 (caveat)", "",
              "All 23,485 ruler coordinates, b = c_young - c_old, Gram-reduced form. S1 includes coordinates no "
              "factor can move (unmatched and missing-mask genes, z = 0 in every d_f) and many noise dimensions; "
              "the numbers are not comparable to S2A_k20 and are not gated. (D) is not run here. (B) uses "
              "sigma^2 from the same draws as full-space squared norms; "
              f"{dgj['s_f_zero_B_full']} of 1,836 factors shrink to zero there.", "",
              md_table(pd.DataFrame(dict(k=list(K_ALL), A=[gv("A", k, "S1") for k in K_ALL],
                                         B=[gv("B", k, "S1") for k in K_ALL], C=[gv("C", k, "S1") for k in K_ALL])),
                       ["k", "A", "B", "C"], ["k", "A raw", "B shrunk", "C null bank"]), "",
              "Solver status in S1: " + "; ".join(f"{k[3:]}: converged {v['ok']}, Gram rank {v['rank']}"
                                                  for k, v in fi.items() if k.startswith("S1_")) + ".", ""]
    vj = OUT / "verify.json"
    if vj.exists():
        v = json.loads(vj.read_text(encoding="utf-8"))
        L += ["## Verify (re-solve from tracked inputs)", "",
              f"`python src/reach2_run.py verify` re-solved sets A-D from `results/reach2/inputs_k20.npz`: reading "
              f"{v['fired']} (stored {v['stored']}); mismatches beyond 1e-4: "
              f"{'; '.join(v['mismatches']) if v['mismatches'] else 'none'}; largest |difference| over the D table "
              f"{v['max_abs_diff_D']:.2e}; {v['solves']} solves, {v['retries']} retries, max gap/tol "
              f"{v['max_ratio']:.3g}.", ""]
    L += ["## What this does not cover", ""] + [f"- {x}" for x in WHAT_NOT] + [""]
    present = sorted(p.name for p in OUT.iterdir() if p.is_file())
    L += ["## Files", "", "`src/reach2_run.py`; `results/reach2/`: " + ", ".join(f"`{p}`" for p in present) + "."]
    return "\n".join(L) + "\n"


def write_findings(stop=None):
    FINDINGS.write_text(render(stop), encoding="utf-8")


# --------------------------------------------------------------------------
# selftest: synthetic matrices only, no project data
# --------------------------------------------------------------------------
def selftest(log, n_f=120, n_t=24):
    rng = np.random.default_rng([SEED, 99])
    D20 = rng.standard_normal((n_f, KDIM)) * 0.8
    C20 = rng.standard_normal((n_f, KDIM)) * 0.8
    b20 = rng.standard_normal(KDIM) * 3
    g = float(np.linalg.norm(b20))
    s_f = np.clip(rng.uniform(-0.2, 1.0, n_f), 0, 1)
    z = rng.standard_normal((n_t, KDIM))
    targets = z / np.linalg.norm(z, axis=1, keepdims=True) * g
    n_cells = rng.integers(13, 342, n_f)
    keep = [f"F{i}" for i in range(n_f)]
    inp = dict(keep=keep, n_cells=n_cells, D20=D20, b20=b20, g=np.array([g]), sigma2=np.ones(n_f), s_f=s_f,
               C20=C20, targets=targets, B20=D20 * s_f[:, None], Cs20=C20, s_c=np.ones(n_f),
               pct=pd.Series(n_cells).rank(pct=True).to_numpy() * 100)
    pd.DataFrame(dict(factor=keep, n_cells=n_cells, s_f=s_f, dnorm_S2A_k20=np.linalg.norm(D20, axis=1),
                      sigma2_S2A_k20=np.ones(n_f))).to_csv(OUT / "factors.csv", index=False)
    np.savez(OUT / "inputs_k20.npz", keep=np.array(keep, dtype=object), n_cells=n_cells, D20=D20, b20=b20,
             g=np.array([g]), sigma2=np.ones(n_f), s_f=s_f, C20=C20, targets=targets)
    plant, res, dtab, calls, problems, status = solve_and_read(inp, log, n_targets=n_t)
    finish(inp, res, dtab, calls, problems, status, log, ctx=None)


def main():
    mode = sys.argv[1] if len(sys.argv) > 1 else "run"
    if mode != "selftest" and (_ST_FAIL_FIRST or os.environ.get("REACH2_SELFTEST_GAP_REL")
                               or os.environ.get("REACH2_OUT") or os.environ.get("REACH2_FINDINGS")):
        raise SystemExit("selftest switches are set in the environment; refusing to run or verify.")
    global OUT, FINDINGS
    if mode == "selftest":
        OUT = Path(sys.argv[2]).resolve()
        FINDINGS = OUT / "FINDINGS_selftest.md"
        os.environ["REACH2_OUT"] = str(OUT)
        os.environ["REACH2_FINDINGS"] = str(FINDINGS)
    OUT.mkdir(parents=True, exist_ok=True)
    if mode == "findings":
        stop = json.loads((OUT / "STOP.json").read_text(encoding="utf-8")) if (OUT / "STOP.json").exists() else None
        write_findings(stop)
        return
    log = Log(OUT / {"run": "run.log", "verify": "verify.log", "selftest": "selftest.log"}[mode])
    log("=" * 80)
    log(f"REACH2 {mode} seed={SEED} k={K_ALL} targets={N_TARGETS} workers={N_WORKERS} "
        f"gap <= {GAP_REL:g} * max(P, 1); retry BVLS {BVLS_STRICT} L-BFGS-B {LB_STRICT}")
    try:
        if mode == "run":
            run(log)
            write_findings(None)
        elif mode == "verify":
            ok = verify(log)
            sys.exit(0 if ok else 1)
        elif mode == "selftest":
            selftest(log)
            write_findings(None)
        else:
            raise SystemExit(f"unknown mode {mode}")
    except Stop as e:
        if mode == "verify":
            dump("verify.json", dict(stop=e.key, message=e.message))
            log(f"[verify STOP] {e.key}: {e.message}")
            sys.exit(2)
        if not (OUT / "STOP.json").exists():
            dump("STOP.json", dict(key=e.key, message=e.message))
            log(f"[STOP] {e.key}: {e.message}")
        write_findings(dict(key=e.key, message=e.message))
        sys.exit(2)


if __name__ == "__main__":
    main()
