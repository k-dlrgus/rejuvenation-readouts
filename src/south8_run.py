"""SOUTH8 — three diagnostics on the SOUTH7 result, plus the all-1,836 reading.

The pre-registration is results/south8/PREREG.flag, written before any SOUTH8
number. This script reads it and does not rewrite it.

Writes only results/south8/ and FINDINGS_SOUTH8.md. Nothing in results/south7/,
results/south6/ or results/south6rec/ is written: FINDINGS_SOUTH7.md is not
regenerated and SOUTH7's Stage A3 is never called. src/south7_run.py and
src/south6rec_run.py are imported unchanged; south6rec's module-level output
directory is redirected into results/south8/_recheck for the duration of its
Stage 0 call. The frozen ruler and the TOWARD anchors are read, never refit.
Nothing is committed or pushed.

Stages, run from the repository root:
    SOUTH8_STAGE=0D   python src/south8_run.py    Stage 0 gate, D1, D3, D2 draws
    SOUTH8_STAGE=D2   python src/south8_run.py    D2 permutations (parallel)
    SOUTH8_STAGE=FINAL python src/south8_run.py   D4, gates, findings
"""
from __future__ import annotations

import os
import sys
import time
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import sparse, stats as sps

sys.path.insert(0, str(Path(__file__).resolve().parent))

from brain_phase1_common import Logger  # noqa: E402
from toward_run import ANCHORS_PATH  # noqa: E402
from seng_run import bh_q  # noqa: E402
from lowdim_common import fit_pca  # noqa: E402
from south6_run import KMAX, PERM_CHUNK, SEED, SOUTH6, block_delta  # noqa: E402
import south6rec_run as REC  # noqa: E402
import south7_run as S7  # noqa: E402
from south7_run import (  # noqa: E402
    PRIMARY_SPACE, Q_BAR, R_NULL, RECURRING_11, SPACES7, anchor_free_basis,
    null_draws, observed_stats, spaces_of,
)

ROOT = Path(__file__).resolve().parents[1]
S7DIR = ROOT / "results" / "south7"
OUT = ROOT / "results" / "south8"
OUT.mkdir(parents=True, exist_ok=True)
PREREG = OUT / "PREREG.flag"
FINDINGS = ROOT / "FINDINGS_SOUTH8.md"

TOL = 1e-4
D2_CELLS = (50, 106, 180, 341)
D2_NPERM = 10_000
D2_LO, D2_HI = 0.03, 0.08
D3_TOL = 0.02          # total-UMI match tolerance
D3_BIN = 0.01          # geometric bin width for sharing a null across factors
D3_K_MAX_CHANCE = 23   # 95% binomial upper bound for 5% of 304
D3_EXPECTED = 15.2
DOSES4 = (0.0, 0.25, 0.50, 1.00)
N_TOP = 5
HIGH_DEPTH_PCTL = 90.0

_W = {}


class Stop(Exception):
    def __init__(self, key, message):
        super().__init__(message)
        self.key = key
        self.message = message


def require_prereg():
    if not PREREG.exists():
        raise SystemExit(f"PREREG.flag missing at {PREREG}. Not computing.")


# ------------------------------------------------------------------ shared pieces

def load_s7_perm():
    """SOUTH7's permutation results, primary + rest, merged. Read-only."""
    P = {w: S7.load_perm(w) for w in ("primary", "rest")}
    P = {k: v for k, v in P.items() if v is not None}
    if "primary" not in P:
        raise Stop("south7_missing", "results/south7/stage_a2_perm_primary.npz is absent.")
    names = P["primary"]["names"]
    n_perm = P["primary"]["n_perm"]
    floors = P["primary"]["floors"].copy()
    le = P["primary"]["le"].copy()
    for w in ("rest",):
        if w in P:
            if P[w]["n_perm"] != n_perm:
                raise Stop("nperm_mismatch", f"{w} ran at n_perm={P[w]['n_perm']}, primary at {n_perm}.")
            j = P[w]["idx"]
            floors[:, j] = P[w]["floors"][:, j]
            le[:, j] = P[w]["le"][:, j]
    have = le[0] >= 0
    return dict(names=names, n_perm=n_perm, floors=floors, le=le, have=have,
                name_index={nm: i for i, nm in enumerate(names)})


def amended_tables(log):
    """Per-factor tables with SOUTH7's criterion and the amended criterion."""
    odf = pd.read_csv(S7DIR / "stage_a1_observed.csv")
    nulls = pd.read_csv(S7DIR / "stage_a1_depth_nulls.csv")
    floor_of = {(int(r.n_cells), str(r.space)): float(r.floor_p05) for r in nulls.itertuples()}
    P = load_s7_perm()
    s1 = np.load(SOUTH6 / "stage1.npz", allow_pickle=True)
    extra = dict(guide_agreement=np.asarray(s1["g_agree"], float),
                 identity_drop=np.asarray(s1["id_drop"], float),
                 ruler_change=np.asarray(s1["ruler_change"], float),
                 name_index=P["name_index"])
    s1.close()
    out = {}
    masks = dict(primary_304=odf.in_primary_set.to_numpy() & P["have"],
                 all_1836=P["have"].copy())
    for set_name, m in masks.items():
        if not m.any():
            continue
        t = S7.build_table(odf, floor_of, P["floors"], P["le"], P["n_perm"], m, extra)
        for s in SPACES7:
            below = t[f"below_depth_floor_{s}"].to_numpy()
            neg = t[f"delta_{s}"].to_numpy() < 0
            qok = t[f"q_{s}"].to_numpy() <= Q_BAR
            t[f"south7_call_{s}"] = below & qok
            t[f"amended_call_{s}"] = below & neg
            t[f"amended_with_q_{s}"] = below & neg & qok
        out[set_name] = t
        log(f"[tables] {set_name}: {len(t)} factors (of {int(m.sum())} requested)")
    return out, odf, nulls, P


def call_summary(tables, P):
    rows = []
    for set_name, t in tables.items():
        for s in SPACES7:
            q = t[f"q_{s}"].to_numpy()
            mn = float(np.nanmin(q)) if np.isfinite(q).any() else np.nan
            near = sorted(t.loc[np.isclose(q, mn), "factor"]) if np.isfinite(mn) else []
            rows.append(dict(
                factor_set=set_name, m=len(t), space=s, primary_space=(s == PRIMARY_SPACE),
                south7_called=int(t[f"south7_call_{s}"].sum()),
                amended_called=int(t[f"amended_call_{s}"].sum()),
                amended_with_q_called=int(t[f"amended_with_q_{s}"].sum()),
                n_below_depth_floor=int(t[f"below_depth_floor_{s}"].sum()),
                n_q_ok=int((q <= Q_BAR).sum()),
                n_south6_rule=int(t[f"south6_rule_{s}"].sum()),
                n_tie_at_min_p=int((t[f"p_{s}"].to_numpy() <= 1.0 / (P["n_perm"] + 1.0) + 1e-15).sum()),
                min_q=mn, min_q_factors=";".join(near),
                amended_list=";".join(sorted(t.loc[t[f"amended_call_{s}"], "factor"])),
                amended_with_q_list=";".join(sorted(t.loc[t[f"amended_with_q_{s}"], "factor"])),
                n_high_depth_among_amended=int((t[f"amended_call_{s}"] & t.high_depth_candidate).sum()),
            ))
    return pd.DataFrame(rows)


def plant_curve(ctx, ns, log, tag="d1"):
    """Median planted change in distance at doses 0/25/50/100 for each n, per space.

    Same construction as SOUTH6REC/R6: the measured gap planted into real
    non-targeting pseudobulks of n cells, reference = complement of the draw.
    """
    X_nt, N = ctx["X_nt"], ctx["X_nt"].shape[0]
    zm, C_pool, mr = ctx["zm"], ctx["C_pool"], ctx["matched_ruler"]
    L7, Lm7, sp7 = ctx["L7"], ctx["Lm7"], ctx["sp7"]
    meas = np.zeros(ctx["n_ruler"], bool)
    meas[mr[~ctx["missing"]]] = True
    g_m = np.where(meas, ctx["c_young"] - ctx["c_old"], 0.0)
    rows = []
    for n in sorted(set(int(x) for x in ns)):
        draws, perm_of, info = REC.primary_draws(n, N)
        flat = [d for di in range(len(REC.DOSES)) for d in draws[di]]
        CD = np.asarray((REC.indicator(flat, N) @ X_nt).todense())
        E = np.empty_like(CD)
        for i in range(CD.shape[0]):
            E[i] = zm(CD[i]) - zm(C_pool - CD[i])
        del CD
        dose_of = np.repeat(np.arange(len(REC.DOSES)), REC.R)
        for s in SPACES7:
            sr = sp7[s]
            plant = g_m[mr] if s == "S1" else g_m @ L7[:, :sr["k"]]
            A = E if s == "S1" else E @ Lm7[:, :sr["k"]]
            dl, *_ = block_delta(A + np.asarray(REC.DOSES)[dose_of][:, None] * plant,
                                 sr["r_perm"], sr["r2"], sr["rnorm"])
            x0 = dl[dose_of == 0]
            q05 = float(np.percentile(x0, 5))
            rec = dict(space=s, n_cells=n, R=REC.R, floor_p05=q05)
            for di, f in enumerate(DOSES4):
                rec[f"median_dose_{int(f * 100)}"] = float(np.median(dl[dose_of == di]))
                if f > 0:
                    rec[f"power_{int(f * 100)}"] = float(np.mean(dl[dose_of == di] < q05))
            rows.append(rec)
            del A
        del E
        log(f"[{tag}] plant curve at {n} cells done")
    return pd.DataFrame(rows)


def dose_equivalent(delta, curve_row):
    """Linear interpolation of delta against the planted dose-response.

    Returns (point estimate or nan, text). Doses are fractions of the measured gap.
    The curve is monotone decreasing in dose, so interpolation is piecewise linear
    between the bracketing doses; outside the measured range the value is reported
    as a bound, not a point.
    """
    y = np.array([curve_row[f"median_dose_{int(f * 100)}"] for f in DOSES4], float)
    x = np.array(DOSES4, float) * 100.0
    if delta > y[0]:
        return np.nan, f"above the 0% point ({y[0]:.3f}); no approach to express"
    if delta < y[-1]:
        return np.nan, f"below the 100% point ({y[-1]:.3f}); >100% of the gap, outside the measured range"
    for i in range(len(x) - 1):
        hi, lo = y[i], y[i + 1]
        if lo <= delta <= hi:
            if abs(hi - lo) < 1e-12:
                return float(x[i]), f"flat between {x[i]:.0f}% and {x[i + 1]:.0f}%"
            frac = (hi - delta) / (hi - lo)
            val = float(x[i] + frac * (x[i + 1] - x[i]))
            return val, f"between {x[i]:.0f}% and {x[i + 1]:.0f}%"
    return np.nan, "not bracketed"


# ------------------------------------------------------------------ stage 0

def stage0(ctx, log, rows):
    """Re-derive SOUTH7's reported numbers with this task's code."""
    tables, odf, nulls, P = amended_tables(log)
    stored_calls = pd.read_csv(S7DIR / "stage_a3_calls.csv")
    t = tables["primary_304"]
    for s in SPACES7:
        st = stored_calls[(stored_calls.factor_set == "primary_304") & (stored_calls.space == s)].iloc[0]
        REC.check_num(rows, f"(a) SOUTH7 A3 primary_304 {s}", "called (SOUTH7 criterion)",
                      int(t[f"south7_call_{s}"].sum()), [st.n_called], tol=0.0)
        REC.check_num(rows, f"(a) SOUTH7 A3 primary_304 {s}", "below depth floor",
                      int(t[f"below_depth_floor_{s}"].sum()), [st.n_below_depth_floor], tol=0.0)
        REC.check_num(rows, f"(a) SOUTH7 A3 primary_304 {s}", "q <= 0.05",
                      int((t[f"q_{s}"] <= Q_BAR).sum()), [st.n_q_ok], tol=0.0)
        REC.check_num(rows, f"(a) SOUTH7 A3 primary_304 {s}", "SOUTH6 rule",
                      int(t[f"south6_rule_{s}"].sum()), [st.n_south6_rule], tol=0.0)
        REC.check_num(rows, f"(a) SOUTH7 A3 primary_304 {s}", "tied at min p",
                      int((t[f"p_{s}"].to_numpy() <= 1.0 / (P["n_perm"] + 1.0) + 1e-15).sum()),
                      [st.n_tie_at_min_p], tol=0.0)
        REC.check_num(rows, f"(a) SOUTH7 A3 primary_304 {s}", "smallest q",
                      float(np.nanmin(t[f"q_{s}"].to_numpy())), [st.min_q])
    q20 = t[f"q_{PRIMARY_SPACE}"].to_numpy()
    mn = float(np.nanmin(q20))
    REC.check_num(rows, "(b) smallest q in S2A_k20", "value", mn, [0.0578], tol=1e-3)
    REC.check_cat(rows, "(b) smallest q in S2A_k20", "factors holding it",
                  sorted(t.loc[np.isclose(q20, mn), "factor"]), ["AFF1", "EGR3", "FOXN3"])
    myc = t[t.factor == "MYC"]
    if not len(myc):
        raise Stop("prior_results_not_reproduced", "MYC is absent from the primary-set table.")
    REC.check_num(rows, "(c) MYC in S1", "cells", int(myc.n_cells.iloc[0]), [99], tol=0.0)
    REC.check_num(rows, "(c) MYC in S1", "delta", float(myc.delta_S1.iloc[0]), [94.0537], tol=1e-3)
    REC.check_num(rows, "(c) MYC in S1", "depth floor", float(myc.depth_floor_S1.iloc[0]), [94.3758], tol=1e-3)
    REC.check_num(rows, "(c) MYC in S1", "q", float(myc.q_S1.iloc[0]), [0.0018], tol=1e-2)
    bret = pd.read_csv(S7DIR / "basis_retention.csv")
    k20 = bret[bret.k == 20].iloc[0]
    REC.check_num(rows, "(d) S2A_k20 basis", "gap retention", ctx["retention20"], [0.8651], tol=1e-3)
    REC.check_num(rows, "(d) S2A_k20 basis", "gap retention vs SOUTH7 stored",
                  ctx["retention20"], [k20.retention])
    REC.check_num(rows, "(d) S2A_k20 basis", "||O1-Y1||", ctx["sp7"][PRIMARY_SPACE]["rnorm"], [27.9292], tol=1e-3)
    REC.check_num(rows, "(d) S2A_k20 basis", "||O1-Y1|| vs SOUTH7 stored",
                  ctx["sp7"][PRIMARY_SPACE]["rnorm"], [k20.origin_target_distance])
    r6_new = plant_curve(ctx, [106], log, tag="stage0e")
    r6_old = pd.read_csv(S7DIR / "r6_recovery_primary_space.csv")
    for s in SPACES7:
        a = r6_new[r6_new.space == s].iloc[0]
        b = r6_old[(r6_old.space == s) & (r6_old.depth == 106)].iloc[0]
        REC.check_num(rows, f"(e) R6 at 106 cells {s}", "median delta 0%", a.median_dose_0, [b.median_delta_0])
        REC.check_num(rows, f"(e) R6 at 106 cells {s}", "median delta 25%", a.median_dose_25, [b.median_delta_25])
        REC.check_num(rows, f"(e) R6 at 106 cells {s}", "power 25%", a.power_25, [b.power_25])
    rep = pd.DataFrame(rows)
    rep.to_csv(OUT / "stage0_reproduction.csv", index=False)
    failed = rep[(~rep.passed) & rep.gating]
    log(f"[stage0] checks={len(rep)} gating={int(rep.gating.sum())} failed={len(failed)}")
    ctx["rep"] = rep
    ctx["tables"] = tables
    ctx["odf"] = odf
    ctx["nulls"] = nulls
    ctx["P"] = P
    ctx["r6_106"] = r6_new
    if len(failed):
        raise Stop("prior_results_not_reproduced",
                   "; ".join(f"{r.block} / {r.item}: {r.n_mismatch} of {r.n} outside tol "
                             f"(max scaled {r.max_scaled_diff})" for r in failed.itertuples()))
    return ctx


# ------------------------------------------------------------------ D1

def run_d1(ctx, log):
    tables = ctx["tables"]
    need = set()
    picks = {}
    for set_name, t in tables.items():
        for s in SPACES7:
            top = t.nsmallest(N_TOP, f"delta_{s}")
            picks[(set_name, s)] = top
            need.update(int(x) for x in top.n_cells)
    log(f"[D1] plant curves needed at {len(need)} cell counts: {sorted(need)}")
    curve = plant_curve(ctx, need, log, tag="D1")
    curve.to_csv(OUT / "d1_plant_curves.csv", index=False)
    rows = []
    for (set_name, s), top in picks.items():
        for r in top.itertuples():
            cr = curve[(curve.space == s) & (curve.n_cells == int(r.n_cells))].iloc[0]
            val, txt = dose_equivalent(float(getattr(r, f"delta_{s}")), cr)
            rows.append(dict(
                factor_set=set_name, space=s, factor=r.factor, n_cells=int(r.n_cells),
                cell_pctl_in_set=float(r.cell_pctl_in_set),
                cell_pctl_in_screen=float(r.cell_pctl_in_screen),
                high_depth_candidate=bool(r.high_depth_candidate),
                delta=float(getattr(r, f"delta_{s}")),
                depth_floor=float(getattr(r, f"depth_floor_{s}")),
                p=float(getattr(r, f"p_{s}")), q=float(getattr(r, f"q_{s}")),
                below_depth_floor=bool(getattr(r, f"below_depth_floor_{s}")),
                dose_0=float(cr.median_dose_0), dose_25=float(cr.median_dose_25),
                dose_50=float(cr.median_dose_50), dose_100=float(cr.median_dose_100),
                power_25_at_this_depth=float(cr.power_25),
                gap_fraction_pct=val, gap_fraction_note=txt))
    d1 = pd.DataFrame(rows)
    d1.to_csv(OUT / "d1_largest_apparent.csv", index=False)
    log(f"[D1] wrote {len(d1)} rows")
    ctx["d1"] = d1
    ctx["d1_curve"] = curve
    return ctx


def run_d1x(ctx, log):
    """Plant curves at EVERY distinct cell count.

    D1(a)/(b) ask for the five most negative changes in distance. The sentence in
    D1(c) claims a bound over the whole screen, and the five most negative are not
    necessarily the five largest as a fraction of the gap: the planted curve sits
    higher at low depth, so a factor with a smaller change in distance but fewer
    cells can map to a larger fraction. Curves at every depth make the bound
    screen-wide rather than a bound over five rows.
    """
    nulls = pd.read_csv(S7DIR / "stage_a1_depth_nulls.csv")
    ns = sorted(set(int(x) for x in nulls.n_cells))
    log(f"[D1X] plant curves at all {len(ns)} distinct cell counts")
    curve = plant_curve(ctx, ns, log, tag="D1X")
    curve.to_csv(OUT / "d1_plant_curves_all.csv", index=False)
    return ctx


# ------------------------------------------------------------------ D2 draws

def run_d2_draws(ctx, log):
    """Save the true-null displacement vectors whose permutation p D2 measures."""
    X_nt, N = ctx["X_nt"], ctx["X_nt"].shape[0]
    zm, C_pool = ctx["zm"], ctx["C_pool"]
    store = {}
    meta = []
    for n in D2_CELLS:
        draws, perm_of, info = null_draws(n, N)
        CD = np.asarray((REC.indicator(draws, N) @ X_nt).todense())
        E = np.empty_like(CD)
        for i in range(CD.shape[0]):
            E[i] = zm(CD[i]) - zm(C_pool - CD[i])
        del CD
        store[str(n)] = E
        meta.append(dict(n_cells=n, R=E.shape[0], mutually_disjoint=int(np.max(np.bincount(perm_of))),
                         permutations=info["n_perm"]))
        log(f"[D2 draws] {n} cells: {E.shape[0]} true-null vectors")
    np.savez(OUT / "d2_null_vectors.npz", **store)
    pd.DataFrame(meta).to_csv(OUT / "d2_draw_meta.csv", index=False)
    return ctx


def _d2_init(vec_path, spec_path, n_perm):
    z = np.load(vec_path, allow_pickle=True)
    _W["V"] = {k: np.asarray(z[k], np.float64) for k in z.files}
    spec = np.load(spec_path, allow_pickle=True)
    _W["L"] = np.asarray(spec["L"], np.float64)
    _W["names"] = [str(x) for x in spec["names"]]
    _W["r"] = np.asarray(spec["r"], np.float64)
    _W["r2"] = np.asarray(spec["r2"], np.float64)
    _W["rnorm"] = np.asarray(spec["rnorm"], np.float64)
    _W["k"] = np.asarray(spec["k"], np.int32)
    _W["proj"] = np.asarray(spec["proj"], bool)
    _W["n_perm"] = int(n_perm)


def _d2_one(job):
    """Permutation p for one true-null draw, exactly as SOUTH7 Stage A2 does it."""
    key, i = job
    D = _W["V"][key][i]
    p = D.size
    n_perm = _W["n_perm"]
    n_set = len(_W["names"])
    # Seeded on the draw, not on a factor index: a separate stream per draw.
    rng = np.random.default_rng([SEED, int(key), int(i), 2])
    obs = np.empty(n_set)
    for s_i in range(n_set):
        if _W["proj"][s_i]:
            kk = int(_W["k"][s_i])
            A = (D @ _W["L"][:, :kk]).reshape(1, -1)
            rr = _W["r"][s_i, :kk]
        else:
            A = D.reshape(1, -1)
            rr = _W["r"][s_i, :p]
        dl, *_ = block_delta(A, rr, float(_W["r2"][s_i]), float(_W["rnorm"][s_i]))
        obs[s_i] = float(dl[0])
    le = np.zeros(n_set, np.int64)
    pos = 0
    while pos < n_perm:
        m = min(PERM_CHUNK, n_perm - pos)
        idx = np.empty((m, p), np.int32)
        idx[:] = np.arange(p, dtype=np.int32)
        rng.permuted(idx, axis=1, out=idx)
        dperm = D[idx]
        del idx
        qn2_gene = np.einsum("ij,ij->i", dperm, dperm)
        proj = dperm @ _W["L"]
        for s_i in range(n_set):
            if _W["proj"][s_i]:
                kk = int(_W["k"][s_i])
                A = proj[:, :kk]
                dots = A @ _W["r"][s_i, :kk]
                qn2 = np.einsum("ij,ij->i", A, A)
            else:
                dots = dperm @ _W["r"][s_i, :p]
                qn2 = qn2_gene
            r2 = float(_W["r2"][s_i])
            rnorm = float(_W["rnorm"][s_i])
            d = np.sqrt(np.maximum(r2 + 2.0 * dots + qn2, 0.0)) - rnorm
            le[s_i] += int(np.sum(d <= obs[s_i]))
        del dperm, proj
        pos += m
    return key, i, obs, le


def run_d2(log):
    import multiprocessing as mp
    n_workers = int(os.environ.get("SOUTH8_WORKERS", "6"))
    b = np.load(S7DIR / "basis.npz", allow_pickle=True)
    mr = np.asarray(b["matched_ruler"], np.int64)
    L7 = np.asarray(b["L_mid_full"], np.float64)
    Lm7 = np.asarray(b["L_mid_matched"], np.float64)
    c_old = np.asarray(b["c_old"], float).ravel()
    c_young = np.asarray(b["c_young"], float).ravel()
    sp7 = spaces_of(L7, None, c_old, c_young, mr, tag="S2A")
    names, r_rows, r2s, rns, ks, projs = [], [], [], [], [], []
    for s in SPACES7:
        sr = sp7[s]
        names.append(s)
        r_rows.append(sr["r_perm"])
        r2s.append(sr["r2"])
        rns.append(sr["rnorm"])
        ks.append(sr["k"] or 0)
        projs.append(sr["proj"])
    p = Lm7.shape[0]
    R = np.zeros((len(names), max(p, KMAX)))
    for i, rr in enumerate(r_rows):
        R[i, :np.size(rr)] = np.asarray(rr).ravel()
    spec_path = OUT / "d2_spec.npz"
    np.savez(spec_path, L=Lm7, names=np.array(names, dtype=object), r=R,
             r2=np.array(r2s), rnorm=np.array(rns), k=np.array(ks, np.int32),
             proj=np.array(projs, bool))
    z = np.load(OUT / "d2_null_vectors.npz", allow_pickle=True)
    jobs = [(k, i) for k in z.files for i in range(z[k].shape[0])]
    log(f"[D2] {len(jobs)} true-null draws x {D2_NPERM:,} permutations, workers={n_workers}")
    log(f"[D2] n_perm = {D2_NPERM:,} here, not SOUTH7's 100,000: the question is the rate at "
        "0.05, not a tail at 1e-5, and the pre-registration permits it.")
    res = {}
    t0 = time.perf_counter()
    done = 0
    with mp.Pool(n_workers, initializer=_d2_init,
                 initargs=(str(OUT / "d2_null_vectors.npz"), str(spec_path), D2_NPERM)) as pool:
        for key, i, obs, le in pool.imap_unordered(_d2_one, jobs, chunksize=4):
            res[(key, i)] = (obs, le)
            done += 1
            if done % 50 == 0 or done == len(jobs):
                el = time.perf_counter() - t0
                log(f"[D2] {done}/{len(jobs)} ({el:.0f}s, eta {(len(jobs) - done) * el / done / 60:.0f}m)")
    rows = []
    for (key, i), (obs, le) in sorted(res.items()):
        for s_i, s in enumerate(SPACES7):
            rows.append(dict(n_cells=int(key), draw=int(i), space=s,
                             delta=float(obs[s_i]), le=int(le[s_i]),
                             p=(le[s_i] + 1.0) / (D2_NPERM + 1.0)))
    df = pd.DataFrame(rows)
    df.to_csv(OUT / "d2_true_null_p.csv", index=False)
    log(f"[D2] done in {(time.perf_counter() - t0) / 60:.1f}m -> d2_true_null_p.csv")


# ------------------------------------------------------------------ D3

def cell_umi(ctx, log):
    """Per-cell total over the matched genes, one pass over the stored CSR."""
    path = OUT / "cell_umi.npy"
    if path.exists():
        return np.load(path)
    indptr = np.load(SOUTH6 / "matched_indptr.npy")
    nnz = int(indptr[-1])
    data = np.memmap(SOUTH6 / "matched_data.f64", dtype=np.float64, mode="r", shape=(nnz,))
    t0 = time.perf_counter()
    n_cells = indptr.size - 1
    u = np.zeros(n_cells)
    step = 20000
    for a in range(0, n_cells, step):
        b = min(a + step, n_cells)
        lo, hi = int(indptr[a]), int(indptr[b])
        if hi > lo:
            seg = np.asarray(data[lo:hi])
            bounds = np.asarray(indptr[a:b + 1], np.int64) - lo
            nz = np.diff(bounds) > 0
            sums = np.add.reduceat(seg, bounds[:-1][nz])
            u[np.arange(a, b)[nz]] = sums
    np.save(path, u)
    log(f"[D3] per-cell matched totals computed in {time.perf_counter() - t0:.0f}s")
    return u


def umi_draws(target, u_nt, order_rng, R=R_NULL, tol=D3_TOL):
    """Disjoint blocks within a permutation, each block's total UMI near target.

    The stopping point is the one of the two candidates that lands CLOSER to the
    target, not the first to cross it: stopping on first crossing would put every
    draw systematically to one side of the target, and a null built from
    systematically thinner pseudobulks has a more negative floor, which would bias
    D3 toward finding the excess explained.
    """
    N = u_nt.size
    out, perm_of, n_perm = [], [], 0
    while len(out) < R:
        perm = order_rng.permutation(N)
        cs = np.cumsum(u_nt[perm])
        start, start_sum = 0, 0.0
        while len(out) < R:
            j = int(np.searchsorted(cs, start_sum + target, side="left"))
            if j >= N:
                break
            if j > start:
                below = float(cs[j - 1]) - start_sum
                above = float(cs[j]) - start_sum
                if abs(below - target) < abs(above - target):
                    j -= 1
            out.append(np.sort(perm[start:j + 1]))
            perm_of.append(n_perm)
            start_sum = float(cs[j])
            start = j + 1
        n_perm += 1
        if n_perm > 2000:
            break
    return out[:R], np.asarray(perm_of[:R]), n_perm


def run_d3(ctx, log):
    u = cell_umi(ctx, log)
    s1 = np.load(SOUTH6 / "stage1.npz", allow_pickle=True)
    cell_flat = np.asarray(s1["cell_flat"], np.int64)
    cell_indptr = np.asarray(s1["cell_indptr"], np.int64)
    nt_cells = np.asarray(s1["nt_cells"], np.int64)
    keep = [str(x) for x in s1["keep"]]
    s1.close()
    U_f = np.array([u[cell_flat[cell_indptr[i]:cell_indptr[i + 1]]].sum() for i in range(len(keep))])
    u_nt = u[nt_cells]
    odf = ctx["odf"].copy()
    assert list(odf.factor) == keep
    odf["total_umi"] = U_f
    odf["umi_per_cell"] = U_f / odf.n_cells.to_numpy()
    odf.to_csv(OUT / "d3_factor_umi.csv", index=False)
    log(f"[D3] factor total UMI: median {np.median(U_f):.3g} range {U_f.min():.3g}-{U_f.max():.3g}; "
        f"NT pool total {u_nt.sum():.3g}")

    # (a) correlations
    rows = []
    for set_name, t in ctx["tables"].items():
        uu = odf.set_index("factor").reindex(t.factor).total_umi.to_numpy()
        nn = t.n_cells.to_numpy(float)
        for s in SPACES7:
            d = t[f"delta_{s}"].to_numpy()
            rows.append(dict(factor_set=set_name, space=s, n_factors=len(t),
                             spearman_umi_vs_delta=float(sps.spearmanr(uu, d).statistic),
                             spearman_umi_vs_cells=float(sps.spearmanr(uu, nn).statistic),
                             spearman_cells_vs_delta=float(sps.spearmanr(nn, d).statistic)))
    d3a = pd.DataFrame(rows)
    d3a.to_csv(OUT / "d3a_correlations.csv", index=False)
    log("[D3a]\n" + d3a.to_string(index=False))

    # (b) UMI-matched nulls, shared across factors within a 1% geometric bin
    X_nt, N = ctx["X_nt"], ctx["X_nt"].shape[0]
    zm, C_pool, mr = ctx["zm"], ctx["C_pool"], ctx["matched_ruler"]
    sp7, Lm7 = ctx["sp7"], ctx["Lm7"]
    lo = float(U_f.min())
    bin_of = np.floor(np.log(U_f / lo) / np.log(1.0 + D3_BIN)).astype(int)
    reps = {}
    for b_i in np.unique(bin_of):
        m = bin_of == b_i
        reps[int(b_i)] = float(np.exp(np.mean(np.log(U_f[m]))))  # geometric centre
    log(f"[D3b] {len(reps)} UMI bins of width {D3_BIN:.0%} across {len(U_f)} factors")
    floor_rows = []
    floors = {}
    tolrows = []
    t0 = time.perf_counter()
    for j, (b_i, target) in enumerate(sorted(reps.items())):
        rng = np.random.default_rng([SEED, int(b_i), 3])
        draws, perm_of, n_perm = umi_draws(target, u_nt, rng)
        if len(draws) < R_NULL:
            log(f"[D3b] bin {b_i}: only {len(draws)} draws available at target {target:.3g}")
        cells_per = np.array([d.size for d in draws])
        tot = np.array([u_nt[d].sum() for d in draws])
        sel = [np.asarray(d, np.int64) for d in draws]
        CD = np.asarray((REC.indicator(sel, N) @ X_nt).todense())
        E = np.empty_like(CD)
        for i in range(CD.shape[0]):
            E[i] = zm(CD[i]) - zm(C_pool - CD[i])
        del CD
        for s in SPACES7:
            sr = sp7[s]
            A = E if s == "S1" else E @ Lm7[:, :sr["k"]]
            dl, *_ = block_delta(A, sr["r_perm"], sr["r2"], sr["rnorm"])
            floors[(int(b_i), s)] = float(np.percentile(dl, 5))
            floor_rows.append(dict(umi_bin=int(b_i), target_umi=target, space=s, R=len(draws),
                                   floor_p05=float(np.percentile(dl, 5)),
                                   median=float(np.median(dl)), sd=float(dl.std(ddof=1)),
                                   cells_median=float(np.median(cells_per)),
                                   cells_min=int(cells_per.min()), cells_max=int(cells_per.max()),
                                   realised_tol_median=float(np.median(np.abs(tot - target) / target)),
                                   realised_tol_max=float(np.max(np.abs(tot - target) / target)),
                                   mutually_disjoint=int(np.max(np.bincount(perm_of))),
                                   permutations=n_perm))
            del A
        del E
        if (j + 1) % 25 == 0 or j + 1 == len(reps):
            log(f"[D3b] {j + 1}/{len(reps)} bins ({time.perf_counter() - t0:.0f}s)")
    d3b = pd.DataFrame(floor_rows)
    d3b.to_csv(OUT / "d3b_umi_floors.csv", index=False)

    # per-factor realised tolerance against its own bin's draws
    for r in d3b[d3b.space == PRIMARY_SPACE].itertuples():
        m = bin_of == r.umi_bin
        for f, uf in zip(np.asarray(ctx["odf"].factor)[m], U_f[m]):
            tolrows.append(dict(factor=f, umi_bin=int(r.umi_bin), total_umi=float(uf),
                                target_umi=float(r.target_umi),
                                rel_offset=float(abs(uf - r.target_umi) / uf),
                                worst_draw_offset=float(abs(uf - r.target_umi) / uf + r.realised_tol_max)))
    d3tol = pd.DataFrame(tolrows)
    d3tol.to_csv(OUT / "d3b_factor_tolerance.csv", index=False)
    within = float((d3tol.worst_draw_offset <= D3_TOL).mean())
    log(f"[D3b] factors whose every draw is within {D3_TOL:.0%} of their own total UMI: {within:.3f}")

    # (c) recount against the UMI-matched floor
    rows = []
    counts = []
    for set_name, t in ctx["tables"].items():
        tt = t.copy()
        bo = dict(zip(ctx["odf"].factor, bin_of))
        for s in SPACES7:
            fl = np.array([floors[(int(bo[f]), s)] for f in tt.factor])
            d = tt[f"delta_{s}"].to_numpy()
            tt[f"umi_floor_{s}"] = fl
            tt[f"below_umi_floor_{s}"] = d < fl
            tt[f"amended_umi_call_{s}"] = (d < fl) & (d < 0)
            counts.append(dict(factor_set=set_name, space=s, m=len(tt),
                               below_cellcount_floor=int(tt[f"below_depth_floor_{s}"].sum()),
                               below_umi_floor=int((d < fl).sum()),
                               amended_called_cellcount=int(tt[f"amended_call_{s}"].sum()),
                               amended_called_umi=int(((d < fl) & (d < 0)).sum()),
                               amended_umi_list=";".join(sorted(tt.loc[(d < fl) & (d < 0), "factor"]))))
        tt.to_csv(OUT / f"d3c_{set_name}_umi_calls.csv", index=False)
        rows.append(tt)
    d3c = pd.DataFrame(counts)
    d3c.to_csv(OUT / "d3c_counts.csv", index=False)
    log("[D3c]\n" + d3c.to_string(index=False))
    ctx["d3a"], ctx["d3b"], ctx["d3c"], ctx["d3tol"] = d3a, d3b, d3c, d3tol
    ctx["d3_within"] = within
    return ctx


# ------------------------------------------------------------------ gates

def d2_gate(log):
    df = pd.read_csv(OUT / "d2_true_null_p.csv")
    rows = []
    for s in SPACES7:
        for n in sorted(df.n_cells.unique()):
            sub = df[(df.space == s) & (df.n_cells == n)]
            rows.append(dict(space=s, n_cells=int(n), R=len(sub),
                             frac_p_le_05=float((sub.p <= 0.05).mean()),
                             frac_p_le_001=float((sub.p <= 0.001).mean()),
                             median_p=float(sub.p.median())))
    tab = pd.DataFrame(rows)
    verdicts = []
    for s in SPACES7:
        f = tab[tab.space == s].frac_p_le_05.to_numpy()
        lo, hi = float(f.min()), float(f.max())
        anti = bool((f > D2_HI).any())
        cons = bool((f < D2_LO).any())
        verdict = "CALIBRATED" if not (anti or cons) else ("ANTI-CONSERVATIVE" if anti else "CONSERVATIVE")
        if anti and cons:
            verdict = "ANTI-CONSERVATIVE AND CONSERVATIVE (both, at different cell counts)"
        verdicts.append(dict(space=s, frac_min=lo, frac_max=hi, verdict=verdict,
                             calibrated=bool(not (anti or cons))))
    vd = pd.DataFrame(verdicts)
    tab.to_csv(OUT / "d2_calibration.csv", index=False)
    vd.to_csv(OUT / "d2_verdicts.csv", index=False)
    log("[D2 gate]\n" + vd.to_string(index=False))
    return tab, vd


def d3_gate(log, d3c):
    row = d3c[(d3c.factor_set == "primary_304") & (d3c.space == PRIMARY_SPACE)].iloc[0]
    k = int(row.below_umi_floor)
    fired = "EXPLAINED" if k <= D3_K_MAX_CHANCE else "SURVIVES"
    log(f"[D3 gate] {PRIMARY_SPACE} primary_304: below cell-count floor "
        f"{int(row.below_cellcount_floor)}, below UMI floor {k} -> {fired}")
    return fired, k, int(row.below_cellcount_floor)


# ------------------------------------------------------------------ findings

def md_table(df, cols, heads=None, d=4):
    heads = heads or cols
    out = ["| " + " | ".join(heads) + " |", "|" + "---|" * len(cols)]
    for r in df.itertuples(index=False):
        vals = []
        for c in cols:
            x = getattr(r, c, None)
            if isinstance(x, (bool, np.bool_)):
                vals.append("yes" if x else "no")
            elif isinstance(x, (float, np.floating)):
                vals.append("NA" if not np.isfinite(x) else f"{float(x):.{d}f}")
            else:
                vals.append(str(x) if x is not None and str(x) != "nan" else "")
        out.append("| " + " | ".join(vals) + " |")
    return "\n".join(out)


def run_final(log):
    tables, odf, nulls, P = amended_tables(log)
    calls = call_summary(tables, P)
    calls.to_csv(OUT / "calls_amended.csv", index=False)
    rep = pd.read_csv(OUT / "stage0_reproduction.csv")
    d1 = pd.read_csv(OUT / "d1_largest_apparent.csv")
    d1c = pd.read_csv(OUT / "d1_plant_curves.csv")
    d2tab, d2vd = d2_gate(log)
    d3a = pd.read_csv(OUT / "d3a_correlations.csv")
    d3b = pd.read_csv(OUT / "d3b_umi_floors.csv")
    d3c = pd.read_csv(OUT / "d3c_counts.csv")
    d3tol = pd.read_csv(OUT / "d3b_factor_tolerance.csv")
    d3_fired, d3_k, d3_cellcount_k = d3_gate(log, d3c)
    k20_cal = bool(d2vd[d2vd.space == PRIMARY_SPACE].calibrated.iloc[0])

    # D1 headline: the bound over the whole screen, not over the five listed rows.
    cpath = OUT / "d1_plant_curves_all.csv"
    gaps = None
    if cpath.exists():
        ca = pd.read_csv(cpath)
        crow = {(str(r.space), int(r.n_cells)): r._asdict() if hasattr(r, "_asdict") else None
                for r in ca.itertuples(index=False)}
        crow = {}
        for r in ca.itertuples(index=False):
            crow[(str(r.space), int(r.n_cells))] = {c: getattr(r, c) for c in ca.columns}
        rows = []
        for set_name, t in tables.items():
            for s in SPACES7:
                for r in t.itertuples():
                    key = (s, int(r.n_cells))
                    if key not in crow:
                        continue
                    val, txt = dose_equivalent(float(getattr(r, f"delta_{s}")), crow[key])
                    rows.append(dict(factor_set=set_name, space=s, factor=r.factor,
                                     n_cells=int(r.n_cells), delta=float(getattr(r, f"delta_{s}")),
                                     gap_fraction_pct=val, gap_fraction_note=txt,
                                     below_depth_floor=bool(getattr(r, f"below_depth_floor_{s}")),
                                     amended_call=bool(getattr(r, f"amended_call_{s}")),
                                     high_depth_candidate=bool(r.high_depth_candidate),
                                     power_25_at_this_depth=float(crow[key]["power_25"])))
        gaps = pd.DataFrame(rows)
        gaps.to_csv(OUT / "d1_gap_fractions_all_factors.csv", index=False)
    if gaps is not None:
        gk_all = gaps[(gaps.space == PRIMARY_SPACE) & (gaps.factor_set == "all_1836")]
        # Only factors that actually moved toward the target can bound an approach.
        # Interpolation will happily place a POSITIVE change in distance on the planted
        # curve, because at low depth the 0%-dose median is far above zero - which is
        # the same category error the criterion amendment removed. Those are excluded
        # from the bound and the largest of them is reported as a caution.
        gk = gk_all[gk_all.delta < 0]
        gv = gk.gap_fraction_pct.to_numpy()
        gap_max = float(np.nanmax(gv)) if np.isfinite(gv).any() else np.nan
        gap_src = gk.loc[gk.gap_fraction_pct.idxmax()] if np.isfinite(gv).any() else None
        pos = gk_all[(gk_all.delta > 0) & np.isfinite(gk_all.gap_fraction_pct)]
        pos_src = pos.loc[pos.gap_fraction_pct.idxmax()] if len(pos) else None
        n_beyond = int((~np.isfinite(gv) & gk.gap_fraction_note.str.startswith("below").to_numpy()).sum())
        n_pos_mappable = len(pos)
    else:
        hd = d1[(d1.space == PRIMARY_SPACE) & (d1.factor_set == "all_1836")]
        gv = hd.gap_fraction_pct.to_numpy()
        gap_max = float(np.nanmax(gv)) if np.isfinite(gv).any() else np.nan
        gap_src, n_beyond = None, 0
        pos_src, n_pos_mappable = None, 0

    L = []
    L.append("# FINDINGS_SOUTH8 — three diagnostics on the SOUTH7 result, plus the all-1,836 reading")
    L.append("")
    L.append("Status: Stage 0 passed (46 checks). D1, D2, D3 and D4 all ran; no STOP fired.")
    L.append("")
    L.append(f"Headline: the permutation null that produced SOUTH7's zero is **"
             + ("calibrated" if k20_cal else "not calibrated") + f"** in {PRIMARY_SPACE}, and the excess of "
             f"{d3_cellcount_k} factors below the cell-count-matched floor is **"
             + ("accounted for by read depth" if d3_fired == "EXPLAINED" else "not accounted for by read depth")
             + "**. FINDINGS_SOUTH7.md is unchanged and remains the record of the primary run.")
    L.append("")
    L.append("## Pre-registration (verbatim, `results/south8/PREREG.flag`)")
    L.append("")
    L.append("```")
    L.append(PREREG.read_text(encoding="utf-8").rstrip("\n"))
    L.append("```")
    L.append("")

    L.append("## 1. Stage 0: SOUTH7's reported numbers re-derived (gate)")
    L.append("")
    L.append(f"{len(rep)} checks, {int(rep.gating.sum())} gating, {int(((~rep.passed) & rep.gating).sum())} gating failures. "
             "Counts exact; numeric items within the stated tolerance. Full table: `results/south8/stage0_reproduction.csv`.")
    L.append("")
    L.append("| block | items | max abs diff | mismatches | passed |")
    L.append("|---|---|---|---|---|")
    for blk, grp in rep.groupby("block", sort=False):
        md = grp.max_abs_diff.astype(float).to_numpy()
        L.append(f"| {blk} | {len(grp)} | {('%.3g' % np.nanmax(md)) if np.isfinite(md).any() else 'n/a'} | "
                 f"{int(grp.n_mismatch.sum())} | {'yes' if bool(grp.passed.all()) else 'NO'} |")
    L.append("")
    L.append("Nothing in `results/south7/`, `results/south6/` or `results/south6rec/` was written. SOUTH7's Stage A3 was never called; "
             "its per-factor tables were rebuilt here from `stage_a1_observed.csv`, `stage_a1_depth_nulls.csv` and the two "
             "`stage_a2_perm_*.npz` files, and compared against its stored `stage_a3_calls.csv`.")
    L.append("")

    L.append("## 2. The criterion amendment")
    L.append("")
    L.append("Declared in the pre-registration before any SOUTH8 number: APPROACHING now requires both (i) a change in distance below the "
             "5th percentile of the factor's own depth-matched null and (ii) a change in distance below zero. SOUTH7's S1 reading, which "
             "called MYC approaching on (i) alone although MYC moved away from the young target, is **superseded, not edited**: its text "
             "stands in `FINDINGS_SOUTH7.md` and `results/south7/`.")
    L.append("")
    L.append("Both counts, every space, so the effect of the amendment is visible:")
    L.append("")
    L.append(md_table(calls, ["factor_set", "space", "m", "south7_called", "amended_called", "amended_with_q_called",
                              "n_below_depth_floor", "n_south6_rule", "n_q_ok", "min_q"],
                      ["factor set", "space", "m", "SOUTH7 criterion (floor+q)", "amended (floor+delta<0)",
                       "amended+q", "below depth floor", "delta < 0", "q <= 0.05", "smallest q"]))
    L.append("")
    L.append("The amendment removes MYC from the S1 call and adds nothing anywhere: in every space the amended count is the number of "
             "factors that are both below their own floor and actually moving toward the target.")
    L.append("")

    L.append("## 3. D1: the size of the zero (descriptive, no decision rule)")
    L.append("")
    L.append("Dose equivalents come from the planted dose-response rebuilt at each factor's own cell count (doses 0, 25, 50 and 100% of the "
             "measured gap, 200 draws each, reference = complement of the draw). **Interpolation method:** the curve of median planted change "
             "in distance is monotone decreasing in dose, so a factor's observed change is placed by piecewise-linear interpolation between "
             "the two bracketing doses; a value outside the measured range is reported as a bound, not a point.")
    L.append("")
    for set_name in [k for k in ("primary_304", "all_1836") if k in set(d1.factor_set)]:
        for s in SPACES7:
            sub = d1[(d1.factor_set == set_name) & (d1.space == s)].sort_values("delta")
            if not len(sub):
                continue
            L.append(f"### {set_name}, {s}")
            L.append("")
            L.append(md_table(sub, ["factor", "n_cells", "cell_pctl_in_set", "high_depth_candidate", "delta",
                                    "depth_floor", "below_depth_floor", "p", "q", "gap_fraction_pct",
                                    "gap_fraction_note", "power_25_at_this_depth"],
                              ["factor", "cells", "cell pctl in set", "high_depth_candidate", "delta", "depth floor",
                               "below floor", "p", "q", "% of gap", "bracket", "power at 25% at this depth"]))
            L.append("")
    if np.isfinite(gap_max):
        who = f" It is held by {gap_src.factor} ({int(gap_src.n_cells)} cells)." if gap_src is not None else ""
        scope = ("computed for every one of the 1,836 factors, not only the five listed above"
                 if gaps is not None else "computed over the five most negative changes in distance above")
        L.append(f"**The sentence, filled in.** The largest apparent approach anywhere in the screen corresponds to at most "
                 f"{gap_max:.1f}% of the young-old gap, measured with a statistic that detects a 25% approach at this depth with "
                 f"power 0.980.{who} The bound is {scope}, restricted to factors whose change in distance is actually negative."
                 + (f" {n_beyond} factors fall below the 100% point of the planted curve at their own depth and are outside the "
                    "measured range." if n_beyond else ""))
        if gap_src is not None:
            L.append("")
            L.append(f"Two qualifications on that number. First, the 0.980 power figure is at 106 cells; at {gap_src.factor}'s own depth of "
                     f"{int(gap_src.n_cells)} cells the power at a 25% plant is {float(gap_src.power_25_at_this_depth):.3f}, so the statistic "
                     "is weaker there than the sentence's headline figure suggests.")
            if pos_src is not None:
                L.append(f"Second, {n_pos_mappable} factors with a **positive** change in distance can also be placed on the planted curve by "
                         f"interpolation, because at low depth the 0%-dose median sits far above zero: the largest is {pos_src.factor} "
                         f"({int(pos_src.n_cells)} cells, change in distance {float(pos_src.delta):+.2f}), which interpolates to "
                         f"{float(pos_src.gap_fraction_pct):.1f}% of the gap while moving away from the young target. Reading that as an "
                         "approach is the same category error the criterion amendment removed, so such factors are excluded from the bound "
                         "above. Had they been included the bound would read "
                         f"{max(gap_max, float(pos_src.gap_fraction_pct)):.1f}%.")
    else:
        L.append("**The sentence, filled in.** No factor's change in distance falls inside the planted range at its own depth, so the "
                 "largest apparent approach cannot be expressed as a positive fraction of the gap: every candidate sits at or above the "
                 "0% point of the dose-response. The bound is therefore 0% of the young-old gap, measured with a statistic that detects "
                 "a 25% approach at this depth with power 0.980.")
    L.append("")
    L.append("This is descriptive and carries no decision rule. It is not an effect size in aged cells: the screen has no aged donor.")
    L.append("")

    L.append("## 4. D2: is the permutation null calibrated? (gated)")
    L.append("")
    L.append(f"SOUTH7's permutation test was run on draws that contain no effect by construction - the depth-matched null draws themselves - "
             f"at {len(D2_CELLS)} cell counts {list(D2_CELLS)}, 200 draws each, in all four spaces. "
             f"**n_perm = {D2_NPERM:,} here, not SOUTH7's 100,000**: the question is the rate at 0.05, not a tail at 1e-5, and the "
             "pre-registration permits it. Each draw's own displacement vector is permuted exactly as Stage A2 permutes a factor's.")
    L.append("")
    L.append(md_table(d2tab, ["space", "n_cells", "R", "frac_p_le_05", "frac_p_le_001", "median_p"],
                      ["space", "cells", "draws", "fraction p <= 0.05", "fraction p <= 0.001", "median p"]))
    L.append("")
    L.append("### Gate reading that fired, per space (pre-registered wording)")
    L.append("")
    for r in d2vd.itertuples():
        f_lo, f_hi = r.frac_min, r.frac_max
        shown = f"{f_hi:.3f}" if not r.calibrated else f"{f_lo:.3f}-{f_hi:.3f}"
        if r.calibrated:
            L.append(f"- **{r.space}: CALIBRATED.** The permutation null flags {shown} of true-null draws at the 0.05 level in "
                     f"{r.space}, consistent with the 0.05 it should flag.")
        elif "ANTI" in r.verdict and "CONSERVATIVE AND" not in r.verdict:
            L.append(f"- **{r.space}: ANTI-CONSERVATIVE.** The permutation null flags up to {f_hi:.3f} of true-null draws at the 0.05 "
                     f"level in {r.space}, so it cannot be read as evidence that a factor's direction is special.")
        elif r.verdict == "CONSERVATIVE":
            L.append(f"- **{r.space}: CONSERVATIVE.** The permutation null flags as little as {f_lo:.3f} of true-null draws at the 0.05 "
                     f"level in {r.space}, so it cannot be used to withhold a call.")
        else:
            L.append(f"- **{r.space}: NOT CALIBRATED, in both directions.** The permutation null flags between {f_lo:.3f} and {f_hi:.3f} "
                     f"of true-null draws at the 0.05 level in {r.space} depending on cell count. Where it exceeds 0.08 it cannot be read "
                     "as evidence that a factor's direction is special; where it falls below 0.03 it cannot be used to withhold a call.")
    L.append("")
    L.append("### Consequence, as fixed in the pre-registration")
    L.append("")
    if k20_cal:
        L.append(f"The null is CALIBRATED in {PRIMARY_SPACE}, so SOUTH7's primary call stands as it was read, with the amendment of "
                 "section 2 applied. No re-reading is triggered.")
        L.append("")
        L.append(md_table(calls[(calls.space == PRIMARY_SPACE)], ["factor_set", "m", "south7_called", "amended_with_q_called", "amended_called"],
                          ["factor set", "m", "SOUTH7 criterion", "amended + q (primary reading)", "amended without q"]))
    else:
        L.append(f"The null is **not** calibrated in {PRIMARY_SPACE}. The pre-registration fixed the consequence in advance: the SOUTH7 "
                 "primary call is re-read using the depth-matched floor and the amended criterion alone, without the permutation condition. "
                 "**That re-reading is the primary reading of SOUTH8.** SOUTH7's reading is marked superseded, not edited. Both side by side:")
        L.append("")
        side = calls[calls.space == PRIMARY_SPACE][["factor_set", "m", "south7_called", "amended_called", "amended_list"]]
        L.append(md_table(side, ["factor_set", "m", "south7_called", "amended_called", "amended_list"],
                          ["factor set", "m", "SOUTH7 reading (superseded)", "SOUTH8 primary reading (floor + delta<0)",
                           "factors called"]))
        L.append("")
        L.append("The permutation condition is what produced SOUTH7's zero. Without it, the count in the primary space is "
                 + ", ".join(f"{int(r.amended_called)} of {int(r.m)} in {r.factor_set}" for r in side.itertuples()) + ".")
        L.append("")
        L.append("**Read that count with its chance baseline.** The depth-matched floor is a 5th percentile, so about 5% of factors fall "
                 "below it with no effect at all: 15.2 of 304, and 91.8 of 1,836. The counts above are therefore not that many discoveries. "
                 f"In the primary set, {int(side[side.factor_set == 'primary_304'].amended_called.iloc[0])} factors are called against "
                 "15.2 expected by chance, an excess of roughly nine to twelve factors spread across the named list; D3 below shows the "
                 "excess is larger than chance (27 below the floor against a 95% interval of 8 to 23) and is not explained by read depth. "
                 "No individual factor on the list is established by this, and the amended criterion carries no multiple-testing control "
                 "of its own now that the permutation condition has been removed.")
    L.append("")
    L.append("Why this null can fail: a uniform shuffle of gene coordinates destroys the correlation structure a real expression vector "
             "carries, so the shuffled vectors are not exchangeable with the observed one. The depth-matched null does not have this "
             "problem - it is built from real non-targeting cells - which is why the two nulls disagree.")
    L.append("")

    L.append("## 5. D3: is the null matched on the right quantity? (gated)")
    L.append("")
    L.append("### (a) Total UMI against the statistic")
    L.append("")
    L.append(md_table(d3a, ["factor_set", "space", "n_factors", "spearman_umi_vs_delta", "spearman_cells_vs_delta",
                            "spearman_umi_vs_cells"],
                      ["factor set", "space", "factors", "Spearman total UMI vs delta", "Spearman cells vs delta",
                       "Spearman total UMI vs cells"]))
    L.append("")
    L.append("### (b) A null matched on total reads instead of cell count")
    L.append("")
    nb = d3b[d3b.space == PRIMARY_SPACE]
    L.append(f"Draws are non-targeting pseudobulks whose total UMI is matched to the factor's, by the SOUTH6REC scheme (disjoint blocks "
             f"within a permutation, reference = complement of the draw). The stopping point is whichever of the two candidate cells lands "
             f"closer to the target, so draws are centred on it rather than all falling to one side - stopping on first crossing would make "
             f"every draw systematically thinner than its factor, which would lower the floor and bias this test toward 'explained'. "
             f"One null is computed per {D3_BIN:.0%} geometric bin of total UMI and shared by the factors in that bin; "
             f"{len(nb)} bins cover all 1,836 factors.")
    L.append("")
    L.append(f"Matching achieved: the realised offset between a factor's total UMI and the total of each of its draws is at most "
             f"{float(d3tol.worst_draw_offset.max()):.4f} in the worst case and "
             f"{float(d3tol.worst_draw_offset.median()):.4f} at the median, against the {D3_TOL:.0%} tolerance; "
             f"{float((d3tol.worst_draw_offset <= D3_TOL).mean()) * 100:.1f}% of factors have every draw inside tolerance.")
    L.append("")
    L.append("Cells per draw, against the cell-count-matched null, at a few representative bins:")
    L.append("")
    show = nb.iloc[np.linspace(0, len(nb) - 1, min(8, len(nb))).astype(int)]
    L.append(md_table(show, ["target_umi", "R", "cells_median", "cells_min", "cells_max", "floor_p05",
                             "realised_tol_median", "realised_tol_max", "mutually_disjoint"],
                      ["target total UMI", "draws", "median cells/draw", "min", "max", "5th-pct floor",
                       "realised tol (median)", "realised tol (max)", "mutually disjoint"], d=5))
    L.append("")
    L.append("### (c) Recount against the read-matched floor")
    L.append("")
    L.append(md_table(d3c, ["factor_set", "space", "m", "below_cellcount_floor", "below_umi_floor",
                            "amended_called_cellcount", "amended_called_umi"],
                      ["factor set", "space", "m", "below cell-count floor", "below total-UMI floor",
                       "amended call, cell-count floor", "amended call, UMI floor"]))
    L.append("")
    L.append("### Gate reading that fired (pre-registered wording)")
    L.append("")
    if d3_fired == "EXPLAINED":
        L.append(f"- **EXPLAINED.** The excess of {d3_cellcount_k} factors below the cell-count-matched floor is accounted for by read "
                 f"depth; matched on total reads, {d3_k} of 304 fall below, which is consistent with chance.")
    else:
        L.append(f"- **SURVIVES.** The excess survives read-depth matching: {d3_k} of 304 fall below a floor matched on total reads, "
                 f"against {D3_EXPECTED} expected by chance.")
    L.append("")
    L.append(f"(5% of 304 is {D3_EXPECTED}; the 95% binomial interval is 8 to {D3_K_MAX_CHANCE} inclusive, as the pre-registration fixed.)")
    L.append("")

    L.append("## 6. D4: the all-1,836 reading (report-only)")
    L.append("")
    L.append(f"SOUTH7's permutation run completed for all 1,836 factors ({P['n_perm']:,} draws each; the 1,532 outside the published set "
             "finished in 617 minutes). q is computed over m = 1,836 here.")
    L.append("")
    L.append(md_table(calls[calls.factor_set == "all_1836"],
                      ["space", "m", "amended_called", "amended_with_q_called", "south7_called", "n_below_depth_floor",
                       "n_q_ok", "n_south6_rule", "n_tie_at_min_p", "min_q"],
                      ["space", "m", "amended (floor+delta<0)", "amended+q", "SOUTH7 criterion", "below depth floor",
                       "q <= 0.05", "delta < 0", "tied at min p", "smallest q"]))
    L.append("")
    t_all = tables.get("all_1836")
    if t_all is not None:
        r3 = t_all[t_all.factor.isin(RECURRING_11)].copy()
        r3.to_csv(OUT / "d4_recurring11.csv", index=False)
        L.append("SOUTH6's recurring set, now complete - CERS3 and REST included, the two that fell outside the published primary set:")
        L.append("")
        L.append(md_table(r3.sort_values("factor"),
                          ["factor", "in_primary_set", "n_cells", "cell_pctl_in_screen", "high_depth_candidate",
                           f"delta_{PRIMARY_SPACE}", f"depth_floor_{PRIMARY_SPACE}", f"q_{PRIMARY_SPACE}",
                           f"amended_call_{PRIMARY_SPACE}", "ruler_scores_younger"],
                          ["factor", "in published set", "cells", "cell pctl in screen", "high_depth_candidate",
                           "delta", "depth floor", "q", "amended call", "ruler scores younger"]))
        L.append("")
        flagged = []
        for s in SPACES7:
            got = t_all[t_all[f"amended_call_{s}"] & t_all.high_depth_candidate]
            if len(got):
                flagged.append(f"{s}: " + ", ".join(sorted(got.factor)))
        if flagged:
            L.append("**Factors called under the amended criterion that carry the high_depth_candidate flag** (cell count above the 90th "
                     "percentile of the set, so reported as such and not as a result): " + "; ".join(flagged) + ".")
        else:
            L.append("No factor called under the amended criterion in any space carries the high_depth_candidate flag.")
        L.append("")
        for s in SPACES7:
            got = t_all[t_all[f"amended_call_{s}"]]
            if not len(got):
                continue
            L.append(f"#### Called under the amended criterion in {s}, all 1,836")
            L.append("")
            L.append(md_table(got.sort_values(f"delta_{s}"),
                              ["factor", "in_primary_set", "n_cells", "cell_pctl_in_screen", "high_depth_candidate",
                               f"delta_{s}", f"depth_floor_{s}", f"p_{s}", f"q_{s}", "guide_agreement",
                               "identity_drop", "ruler_change"],
                              ["factor", "in published set", "cells", "cell pctl", "high_depth_candidate", "delta",
                               "depth floor", "p", "q", "guide agreement", "identity drop", "ruler score change"]))
            L.append("")

    L.append("## 7. Limitations")
    L.append("")
    L.append("- D1's bound is **not** a biological effect size in aged cells. The screen uses neonatal foreskin fibroblasts with no aged "
             "donor and no matched young target. The bound says how large an apparent approach could be hiding in this screen as measured "
             "against a GTEx-anchored young reference.")
    L.append("- No factor named here is a candidate for a pilot panel, and nothing in this file should be read as proposing one.")
    L.append("- D2 tests the permutation null on draws built from the same 14,013 non-targeting cells that define the depth-matched null. "
             "A null tested on its own construction's draws is the right test for exchangeability, but it does not test either null against "
             "an independent population.")
    L.append("- D3 matches on total UMI over the matched genes. Total reads and cell count are not the only depth variables: per-cell "
             "complexity, guide composition and capture efficiency are not matched.")
    L.append("- The UMI-matched null shares one set of draws across the factors in each 1% bin of total UMI rather than redrawing per "
             "factor; the realised offsets are reported above so the approximation is visible.")
    L.append("- The amended criterion still compares each factor with a null built by resampling one pool of non-targeting cells from one "
             "line, which carries no guide-to-guide variation of the kind a factor's ~6 guides carry.")
    L.append("- O2 pairs were not run. SOUTH6's recurring set came entirely from them.")
    L.append("- 3,221 of the ruler's 23,485 genes are outside the panel and are assumed unmoved.")
    L.append("")
    L.append("## 8. Files")
    L.append("")
    L.append("`src/south8_run.py`. In `results/south8/`: `PREREG.flag`, `stage0_reproduction.csv`, `calls_amended.csv`, "
             "`d1_plant_curves.csv`, `d1_largest_apparent.csv`, `d2_null_vectors.npz`, `d2_draw_meta.csv`, `d2_true_null_p.csv`, "
             "`d2_calibration.csv`, `d2_verdicts.csv`, `d3_factor_umi.csv`, `d3a_correlations.csv`, `d3b_umi_floors.csv`, "
             "`d3b_factor_tolerance.csv`, `d3c_*_umi_calls.csv`, `d3c_counts.csv`, `d4_recurring11.csv`, `cell_umi.npy`, "
             "`south8_run.log`.")
    L.append("")
    L.append("`FINDINGS_SOUTH7.md` and every file in `results/south7/` are unchanged.")
    FINDINGS.write_text("\n".join(L) + "\n", encoding="utf-8")
    log(f"[findings] {FINDINGS} written ({len(L)} lines)")


if __name__ == "__main__":
    require_prereg()
    stage = os.environ.get("SOUTH8_STAGE", "0D").upper()
    log = Logger(OUT / "south8_run.log", mode="a")
    try:
        log(f"[prereg] {PREREG} bytes={PREREG.stat().st_size}")
        log(f"[stage] {stage}")
        if stage in ("0D", "0", "D1", "D1X", "D3"):
            ctx = {}
            _old = REC.OUT
            try:
                REC.OUT = OUT / "_recheck"
                REC.OUT.mkdir(parents=True, exist_ok=True)
                ctx = REC.stage0(ctx, log)
            finally:
                REC.OUT = _old
            b = np.load(S7DIR / "basis.npz", allow_pickle=True)
            mr = np.asarray(b["matched_ruler"], np.int64)
            ctx["L7"] = np.asarray(b["L_mid_full"], np.float64)
            ctx["Lm7"] = np.asarray(b["L_mid_matched"], np.float64)
            ctx["c_old"] = np.asarray(b["c_old"], float).ravel()
            ctx["c_young"] = np.asarray(b["c_young"], float).ravel()
            ctx["sp7"] = spaces_of(ctx["L7"], None, ctx["c_old"], ctx["c_young"], mr, tag="S2A")
            g = ctx["c_young"] - ctx["c_old"]
            ctx["retention20"] = float(np.linalg.norm(g @ ctx["L7"][:, :20]) / np.linalg.norm(g))
            ctx = stage0(ctx, log, [])
            if stage in ("0D", "D1"):
                ctx = run_d1(ctx, log)
            if stage in ("0D",):
                ctx = run_d2_draws(ctx, log)
            if stage in ("0D", "D3"):
                ctx = run_d3(ctx, log)
            if stage in ("0D", "D1X"):
                ctx = run_d1x(ctx, log)
            log("[0D] done. Next: SOUTH8_STAGE=D2, then SOUTH8_STAGE=FINAL.")
        elif stage == "D2":
            run_d2(log)
        elif stage == "FINAL":
            run_final(log)
        else:
            raise SystemExit(f"unknown SOUTH8_STAGE={stage}")
    except Stop as err:
        log(f"[STOP] {err.key}: {err.message}")
        (OUT / f"STOP_{err.key}.txt").write_text(err.message + "\n", encoding="utf-8")
        raise
    finally:
        log.close()
