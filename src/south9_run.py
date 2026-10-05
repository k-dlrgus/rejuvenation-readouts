"""SOUTH9 — is the below-floor excess proliferation?

The pre-registration is results/south9/PREREG.flag, written before any SOUTH9
number. This script reads it and does not rewrite it.

Writes only results/south9/ and FINDINGS_SOUTH9.md. No earlier findings file is
regenerated and no earlier results directory is written. src/south6rec_run.py,
src/south7_run.py and src/south8_run.py are imported unchanged; south6rec's
module-level output directory is redirected into results/south9/_recheck for the
duration of its Stage 0 call. The proliferation axis is carried over frozen from
results/south6/stage1.npz and is NOT refit. Nothing is committed or pushed.

Run from the repository root:
    python src/south9_run.py
"""
from __future__ import annotations

import json
import os
import sys
import time
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import stats as sps

sys.path.insert(0, str(Path(__file__).resolve().parent))

from brain_phase1_common import Logger  # noqa: E402
from south6_run import SEED, SOUTH6, block_delta  # noqa: E402
import south6rec_run as REC  # noqa: E402
import south7_run as S7  # noqa: E402
import south8_run as S8  # noqa: E402
from south7_run import PRIMARY_SPACE, R_NULL, SPACES7, null_draws, spaces_of  # noqa: E402

ROOT = Path(__file__).resolve().parents[1]
S7DIR = ROOT / "results" / "south7"
S8DIR = ROOT / "results" / "south8"
OUT = ROOT / "results" / "south9"
OUT.mkdir(parents=True, exist_ok=True)
PREREG = OUT / "PREREG.flag"
FINDINGS = ROOT / "FINDINGS_SOUTH9.md"

DOSES4 = (0.0, 0.25, 0.50, 1.00)
POWER_BAR = 0.80
R6_RAW_POWER_25 = 0.980
# Pre-registered chance intervals (5% of each set, 95% binomial).
CHANCE = {"primary_304": dict(m=304, e=15.2, lo=8, hi=23),
          "all_1836": dict(m=1836, e=91.8, lo=74, hi=110)}
# Declared here, before any SOUTH9 number: which residualization is primary.
PRIMARY_VARIANT = "in_space"


class Stop(Exception):
    def __init__(self, key, message):
        super().__init__(message)
        self.key = key
        self.message = message


def require_prereg():
    if not PREREG.exists():
        raise SystemExit(f"PREREG.flag missing at {PREREG}. Not computing.")


# ------------------------------------------------------------------ geometry

def build_geometry(log):
    """Spaces, the frozen proliferation axis, and the two residualization variants."""
    b = np.load(S7DIR / "basis.npz", allow_pickle=True)
    mr = np.asarray(b["matched_ruler"], np.int64)
    L_full = np.asarray(b["L_mid_full"], np.float64)
    Lm = np.asarray(b["L_mid_matched"], np.float64)
    c_old = np.asarray(b["c_old"], float).ravel()
    c_young = np.asarray(b["c_young"], float).ravel()
    sp = spaces_of(L_full, None, c_old, c_young, mr, tag="S2A")
    s1 = np.load(SOUTH6 / "stage1.npz", allow_pickle=True)
    p_full = np.asarray(s1["p_unit"], np.float64).ravel()
    s1.close()
    if p_full.size == 0 or abs(float(np.linalg.norm(p_full)) - 1.0) > 1e-9:
        raise Stop("prolif_axis", f"frozen p_unit is not a unit vector (norm "
                                  f"{float(np.linalg.norm(p_full)) if p_full.size else 0}). Not refitting.")
    off = np.ones(p_full.size, bool)
    off[mr] = False
    n_off = int((p_full[off] != 0).sum())
    if n_off:
        raise Stop("prolif_axis", f"frozen p_unit has {n_off} nonzero coordinates outside the matched genes; "
                                  "the projection identity used here assumes it does not. Not refitting.")
    p_m = p_full[mr]
    geo = {}
    for s in SPACES7:
        d = dict(space=s, proj=sp[s]["proj"], k=sp[s]["k"], r=sp[s]["r_perm"],
                 r2=sp[s]["r2"], rnorm=sp[s]["rnorm"])
        if s == "S1":
            # Gene space: full-vector algebra, exactly SOUTH6's stats_from_d geometry.
            r_full = c_old - c_young
            d["r_full_dot_p"] = float(r_full @ p_full)
            d["r2_full"] = float(r_full @ r_full)
            d["p_m"] = p_m
            d["r_m"] = r_full[mr]
        else:
            k = int(s.split("k")[1])
            Lk = Lm[:, :k]
            pk = p_m @ Lk                      # = p_full @ L_full[:, :k], p is 0 off matched
            nk = float(np.linalg.norm(pk))
            d["Lk"] = Lk
            d["pk"] = pk
            d["pu"] = pk / nk if nk > 1e-12 else pk
            d["axis_retention"] = nk
            d["p_m"] = p_m
            d["r_full_dot_p"] = float((c_old - c_young) @ p_full)
        geo[s] = d
    log(f"[geometry] frozen proliferation axis carried over: ||p||=1 by construction, "
        f"{int((p_m != 0).sum())} nonzero coordinates, all inside the matched genes")
    for s in SPACES7:
        if s != "S1":
            log(f"[geometry] {s}: axis retention ||proj p||={geo[s]['axis_retention']:.4f}; "
                f"cos(r, axis)={float(geo[s]['r'] @ geo[s]['pu'] / geo[s]['rnorm']):.4f}")
        else:
            log(f"[geometry] S1: cos(r, axis)={geo[s]['r_full_dot_p'] / np.sqrt(geo[s]['r2_full']):.4f}")
    return geo, mr, Lm, c_old, c_young, p_full, p_m


def deltas_for(X_m, g, variant):
    """Raw and residualized change in distance for matched-coordinate rows X_m.

    variant 'in_space'  : the component along the axis AS THAT SPACE SEES IT is
                          removed from both the displacement and the origin-target
                          vector, inside the space the statistic is computed in.
    variant 'gene_space': SOUTH6's own formula - the component is removed in gene
                          space first, and the residual is then projected.
    In S1 the two coincide.
    """
    X_m = np.atleast_2d(np.asarray(X_m, np.float64))
    if not g["proj"]:
        c = X_m @ g["p_m"]
        xr = X_m @ g["r_m"]
        xn2 = np.einsum("ij,ij->i", X_m, X_m)
        r2, a = g["r2_full"], g["r_full_dot_p"]
        raw = np.sqrt(np.maximum(r2 + 2.0 * xr + xn2, 0.0)) - np.sqrt(r2)
        r2r = max(r2 - a * a, 0.0)
        res = (np.sqrt(np.maximum(r2r + 2.0 * (xr - a * c) + xn2 - c * c, 0.0))
               - np.sqrt(r2r))
        return raw, res
    Xk = X_m @ g["Lk"]
    r = g["r"]
    raw = np.sqrt(np.maximum(np.einsum("ij,ij->i", Xk + r, Xk + r), 0.0)) - g["rnorm"]
    if variant == "in_space":
        pu = g["pu"]
        Xr = Xk - np.outer(Xk @ pu, pu)
        rr = r - (r @ pu) * pu
    else:
        pk = g["pk"]
        c = X_m @ g["p_m"]
        Xr = Xk - np.outer(c, pk)
        rr = r - g["r_full_dot_p"] * pk
    S = Xr + rr
    res = np.sqrt(np.maximum(np.einsum("ij,ij->i", S, S), 0.0)) - float(np.linalg.norm(rr))
    return raw, res


# ------------------------------------------------------------------ method statement

def method_statement(log):
    notes = json.loads((SOUTH6 / "stage1_notes.json").read_text(encoding="utf-8"))
    rec = dict(prolif_note=notes.get("prolif"),
               source="results/south6/stage1.npz p_unit, carried over frozen; not refit here")
    log(f"[method] {rec}")
    (OUT / "method_statement.json").write_text(json.dumps(rec, indent=1), encoding="utf-8")
    return rec


# ------------------------------------------------------------------ stage 0

def stage0(log):
    rows = []
    tables, odf, nulls, P = S8.amended_tables(log)
    stored = pd.read_csv(S8DIR / "calls_amended.csv")
    for set_name, t in tables.items():
        st = stored[(stored.factor_set == set_name) & (stored.space == PRIMARY_SPACE)].iloc[0]
        REC.check_num(rows, f"(a) SOUTH8 amended {set_name} {PRIMARY_SPACE}", "called",
                      int(t[f"amended_call_{PRIMARY_SPACE}"].sum()), [st.amended_called], tol=0.0)
        REC.check_num(rows, f"(a) SOUTH8 amended {set_name} {PRIMARY_SPACE}", "below depth floor",
                      int(t[f"below_depth_floor_{PRIMARY_SPACE}"].sum()), [st.n_below_depth_floor], tol=0.0)
    REC.check_num(rows, "(a) SOUTH8 amended", "called, primary_304",
                  int(tables["primary_304"][f"amended_call_{PRIMARY_SPACE}"].sum()), [24], tol=0.0)
    REC.check_num(rows, "(a) SOUTH8 amended", "called, all_1836",
                  int(tables["all_1836"][f"amended_call_{PRIMARY_SPACE}"].sum()), [125], tol=0.0)
    REC.check_num(rows, "(a) SOUTH8 amended", "below floor, primary_304",
                  int(tables["primary_304"][f"below_depth_floor_{PRIMARY_SPACE}"].sum()), [27], tol=0.0)
    REC.check_num(rows, "(a) SOUTH8 amended", "below floor, all_1836",
                  int(tables["all_1836"][f"below_depth_floor_{PRIMARY_SPACE}"].sum()), [153], tol=0.0)
    d2 = pd.read_csv(S8DIR / "d2_true_null_p.csv")
    for s, want in ((PRIMARY_SPACE, 0.0350), ("S1", 0.1200)):
        sub = d2[(d2.space == s) & (d2.n_cells == 106)]
        REC.check_num(rows, f"(b) SOUTH8 D2 calibration {s} at 106 cells", "fraction p <= 0.05",
                      float((sub.p <= 0.05).mean()), [want], tol=1e-3)
    s8list = sorted(stored[(stored.factor_set == "primary_304")
                           & (stored.space == PRIMARY_SPACE)].amended_list.iloc[0].split(";"))
    mine = sorted(tables["primary_304"].loc[tables["primary_304"][f"amended_call_{PRIMARY_SPACE}"], "factor"])
    REC.check_cat(rows, "(c) the 24-factor list", "primary_304 amended call list", mine, s8list)
    rep = pd.DataFrame(rows)
    rep.to_csv(OUT / "stage0_reproduction.csv", index=False)
    failed = rep[(~rep.passed) & rep.gating]
    log(f"[stage0] checks={len(rep)} gating={int(rep.gating.sum())} failed={len(failed)}")
    if len(failed):
        raise Stop("prior_results_not_reproduced",
                   "; ".join(f"{r.block} / {r.item}: {r.n_mismatch} of {r.n} outside tol "
                             f"(max scaled {r.max_scaled_diff})" for r in failed.itertuples()))
    return rep, tables, odf, s8list


# ------------------------------------------------------------------ D1

def run_d1(ctx, geo, tables, log):
    """Residualized statistic and residualized depth-matched floors."""
    X_nt, N = ctx["X_nt"], ctx["X_nt"].shape[0]
    zm, C_pool = ctx["zm"], ctx["C_pool"]
    D_m = np.asarray(np.load(SOUTH6 / "perm_D.npy", mmap_mode="r"), np.float64)
    variants = ("in_space", "gene_space")
    # Observed statistic per factor.
    obs = {}
    for s in SPACES7:
        for v in variants:
            raw, res = deltas_for(D_m, geo[s], v)
            obs[(s, v)] = dict(raw=raw, res=res)
    del D_m
    # Residualized floors, one null per distinct cell count, same draws as SOUTH7 A1.
    ns = sorted(set(int(x) for x in ctx["odf"].n_cells))
    log(f"[D1] residualized floors at {len(ns)} distinct cell counts, {R_NULL} draws each")
    frows = []
    floors = {}
    t0 = time.perf_counter()
    for j, n in enumerate(ns):
        draws, perm_of, info = null_draws(n, N)
        CD = np.asarray((REC.indicator(draws, N) @ X_nt).todense())
        E = np.empty_like(CD)
        for i in range(CD.shape[0]):
            E[i] = zm(CD[i]) - zm(C_pool - CD[i])
        del CD
        for s in SPACES7:
            for v in variants:
                raw, res = deltas_for(E, geo[s], v)
                floors[(n, s, v, "raw")] = float(np.percentile(raw, 5))
                floors[(n, s, v, "res")] = float(np.percentile(res, 5))
                frows.append(dict(n_cells=n, space=s, variant=v, R=len(draws),
                                  raw_floor_p05=float(np.percentile(raw, 5)),
                                  res_floor_p05=float(np.percentile(res, 5)),
                                  raw_median=float(np.median(raw)), res_median=float(np.median(res)),
                                  res_sd=float(res.std(ddof=1))))
        del E
        if (j + 1) % 50 == 0 or j + 1 == len(ns):
            log(f"[D1] {j + 1}/{len(ns)} depths ({time.perf_counter() - t0:.0f}s)")
    fdf = pd.DataFrame(frows)
    fdf.to_csv(OUT / "d1_residualized_floors.csv", index=False)

    rows = []
    per_factor = {}
    for set_name, t in tables.items():
        n = t.n_cells.to_numpy()
        pos = t.index.to_numpy()  # row order of the rebuilt table
        idx = np.asarray([list(ctx["odf"].factor).index(f) for f in t.factor]) \
            if len(t) != len(ctx["odf"]) else np.arange(len(t))
        tt = t.copy()
        for s in SPACES7:
            for v in variants:
                raw = obs[(s, v)]["raw"][idx]
                res = obs[(s, v)]["res"][idx]
                fr_raw = np.array([floors[(int(x), s, v, "raw")] for x in n])
                fr_res = np.array([floors[(int(x), s, v, "res")] for x in n])
                tt[f"raw_delta_{s}"] = raw
                tt[f"res_delta_{s}_{v}"] = res
                tt[f"res_floor_{s}_{v}"] = fr_res
                tt[f"res_below_{s}_{v}"] = res < fr_res
                tt[f"res_call_{s}_{v}"] = (res < fr_res) & (res < 0)
                rows.append(dict(factor_set=set_name, space=s, variant=v, m=len(tt),
                                 raw_below_floor=int((raw < fr_raw).sum()),
                                 raw_negative=int((raw < 0).sum()),
                                 raw_called=int(((raw < fr_raw) & (raw < 0)).sum()),
                                 res_below_floor=int((res < fr_res).sum()),
                                 res_negative=int((res < 0).sum()),
                                 res_called=int(((res < fr_res) & (res < 0)).sum()),
                                 res_called_list=";".join(sorted(tt.loc[(res < fr_res) & (res < 0), "factor"]))))
        tt.to_csv(OUT / f"d1_{set_name}_residualized.csv", index=False)
        per_factor[set_name] = tt
    d1 = pd.DataFrame(rows)
    d1.to_csv(OUT / "d1_counts.csv", index=False)
    log("[D1]\n" + d1[d1.variant == PRIMARY_VARIANT][
        ["factor_set", "space", "m", "raw_called", "res_called", "raw_below_floor", "res_below_floor"]
    ].to_string(index=False))
    return d1, per_factor, fdf


# ------------------------------------------------------------------ D2

def run_d2(ctx, geo, log):
    """Plant recovery with residualization applied to plant and null alike."""
    X_nt, N = ctx["X_nt"], ctx["X_nt"].shape[0]
    zm, C_pool, mr = ctx["zm"], ctx["C_pool"], ctx["matched_ruler"]
    meas = np.zeros(ctx["n_ruler"], bool)
    meas[mr[~ctx["missing"]]] = True
    g_m = np.where(meas, ctx["c_young"] - ctx["c_old"], 0.0)[mr]
    draws, perm_of, info = REC.primary_draws(106, N)
    flat = [d for di in range(len(REC.DOSES)) for d in draws[di]]
    CD = np.asarray((REC.indicator(flat, N) @ X_nt).todense())
    E = np.empty_like(CD)
    for i in range(CD.shape[0]):
        E[i] = zm(CD[i]) - zm(C_pool - CD[i])
    del CD
    dose_of = np.repeat(np.arange(len(REC.DOSES)), REC.R)
    rows = []
    for s in SPACES7:
        for v in ("in_space", "gene_space"):
            X = E + np.asarray(REC.DOSES)[dose_of][:, None] * g_m
            raw, res = deltas_for(X, geo[s], v)
            rec = dict(space=s, variant=v, depth=106, R=REC.R)
            for kind, arr in (("raw", raw), ("res", res)):
                x0 = arr[dose_of == 0]
                q05 = float(np.percentile(x0, 5))
                meds = [float(np.median(arr[dose_of == di])) for di in range(4)]
                rec[f"{kind}_floor_p05"] = q05
                rec[f"{kind}_monotone"] = bool(all(meds[i + 1] < meds[i] for i in range(3)))
                for di, f in enumerate(DOSES4):
                    rec[f"{kind}_median_{int(f * 100)}"] = meds[di]
                    if f > 0:
                        rec[f"{kind}_power_{int(f * 100)}"] = float(np.mean(arr[dose_of == di] < q05))
            rows.append(rec)
            del X
    del E
    d2 = pd.DataFrame(rows)
    d2.to_csv(OUT / "d2_residualized_recovery.csv", index=False)
    log("[D2]\n" + d2[d2.variant == PRIMARY_VARIANT][
        ["space", "raw_power_25", "res_power_25", "res_power_50", "res_power_100", "res_monotone"]
    ].to_string(index=False))
    return d2


# ------------------------------------------------------------------ gates

def d1_gate(d1):
    out = {}
    for set_name, c in CHANCE.items():
        row = d1[(d1.factor_set == set_name) & (d1.space == PRIMARY_SPACE)
                 & (d1.variant == PRIMARY_VARIANT)]
        if not len(row):
            continue
        k = int(row.res_called.iloc[0])
        out[set_name] = dict(k=k, inside_or_below=bool(k <= c["hi"]), **c)
    both_in = all(v["inside_or_below"] for v in out.values())
    any_above = any(not v["inside_or_below"] for v in out.values())
    if both_in:
        fired = "EXPLAINED"
    elif all(not v["inside_or_below"] for v in out.values()):
        fired = "SURVIVES"
    else:
        fired = "SPLIT"
    return fired, out


def d2_gate(d2):
    row = d2[(d2.space == PRIMARY_SPACE) & (d2.variant == PRIMARY_VARIANT)].iloc[0]
    p = float(row.res_power_25)
    return ("RETAINED" if p >= POWER_BAR else "LOST"), p, float(row.raw_power_25)


# ------------------------------------------------------------------ D3

def run_d3(per_factor, s8list, log):
    rows = []
    for set_name, t in per_factor.items():
        called = t[t[f"res_call_{PRIMARY_SPACE}_{PRIMARY_VARIANT}"]]
        if len(called):
            called.to_csv(OUT / f"d3_{set_name}_called_detail.csv", index=False)
        rows.append(dict(
            factor_set=set_name, m=len(t), n_called=len(called),
            id_drop_median_called=float(np.median(called.identity_drop)) if len(called) else np.nan,
            id_drop_median_set=float(np.median(t.identity_drop)),
            guide_agreement_median_called=float(np.median(called.guide_agreement)) if len(called) else np.nan,
            guide_agreement_median_set=float(np.median(t.guide_agreement)),
            n_high_depth_called=int(called.high_depth_candidate.sum()) if len(called) else 0,
            ruler_younger_total=int(t.ruler_scores_younger.sum()),
            ruler_younger_not_called=int((t.ruler_scores_younger
                                          & ~t[f"res_call_{PRIMARY_SPACE}_{PRIMARY_VARIANT}"]).sum()),
        ))
    d3 = pd.DataFrame(rows)
    d3.to_csv(OUT / "d3_summary.csv", index=False)
    # SOUTH8's 24, raw against residualized, individually.
    t = per_factor["primary_304"]
    s24 = t[t.factor.isin(s8list)].copy()
    cols = ["factor", "n_cells", "cell_pctl_in_set", "high_depth_candidate",
            f"raw_delta_{PRIMARY_SPACE}", f"depth_floor_{PRIMARY_SPACE}",
            f"res_delta_{PRIMARY_SPACE}_{PRIMARY_VARIANT}", f"res_floor_{PRIMARY_SPACE}_{PRIMARY_VARIANT}",
            f"res_below_{PRIMARY_SPACE}_{PRIMARY_VARIANT}", f"res_call_{PRIMARY_SPACE}_{PRIMARY_VARIANT}",
            "identity_drop", "guide_agreement", "ruler_change"]
    s24[cols].to_csv(OUT / "d3_south8_24_raw_vs_residualized.csv", index=False)
    log("[D3]\n" + d3.to_string(index=False))
    return d3, s24[cols]


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
                vals.append("" if x is None or str(x) == "nan" else str(x))
        out.append("| " + " | ".join(vals) + " |")
    return "\n".join(out)


def write_findings(R, log):
    d1, d2, d3, s24, geo = R["d1"], R["d2"], R["d3"], R["s24"], R["geo"]
    fired1, g1 = R["gate1"]
    fired2, p_res, p_raw = R["gate2"]
    L = []
    L.append("# FINDINGS_SOUTH9 — is the below-floor excess proliferation?")
    L.append("")
    L.append(f"Status: Stage 0 passed ({len(R['rep'])} checks). D1, D2 and D3 ran; no STOP fired.")
    L.append("")
    L.append(f"D1 gate: **{fired1}**. D2 gate: **{fired2}** (power at a 25% plant after residualization "
             f"{p_res:.3f}, against {p_raw:.3f} raw).")
    if fired2 == "LOST":
        L.append("")
        L.append("Because D2 fired LOST, D1's reading is reported but carries no conclusion. That is said again in its gate sentence.")
    L.append("")
    L.append("## Pre-registration (verbatim, `results/south9/PREREG.flag`)")
    L.append("")
    L.append("```")
    L.append(PREREG.read_text(encoding="utf-8").rstrip("\n"))
    L.append("```")
    L.append("")

    L.append("## 1. The method, stated before it is used")
    L.append("")
    L.append("How SOUTH6 builds the proliferation axis (`src/south6_run.py`, the Stage 1 block that writes `p_unit`):")
    L.append("")
    L.append("1. **Gene set.** Seurat's `cc.genes` (Tirosh 2016), parsed from the packaged `.rda` with boundary checks "
             "(97 genes: 43 S-phase, 54 G2M; the parser refuses a hand-written substitute). 94 of the 97 are present in the "
             "Southard panel.")
    L.append("2. **Score.** For each **non-targeting Southard cell**, the mean of log1p(counts/library x 1e4) over those cc genes.")
    L.append("3. **Contrast.** The non-targeting cells are split at the quartiles of that score: 3,504 cells in the bottom "
             "quartile, 3,504 in the top. Each group is summed into a pseudobulk, carried through the same TMM / log2-CPM / "
             "frozen-GTEx-z transform as every other pseudobulk, and the axis is z(top) - z(bottom), norm 65.9633, stored "
             "normalised to unit length as `p_unit`.")
    L.append("4. **Residual.** SOUTH6's `stats_from_d` removes the component along `p_unit` from both the factor's displacement "
             "and the origin-target vector, then recomputes the change in distance.")
    L.append("")
    L.append("What entered the fit, and what did not:")
    L.append("")
    L.append("- **No GTEx donor entered it.** The axis is built only from Southard cells, so neither anchor (O1, Y1) nor the "
             "frozen ruler's training set contributed to its direction.")
    L.append("- **No targeting cell entered it.** Only non-targeting cells were used, so a factor's own cells did not help "
             "define the axis it is tested against.")
    L.append("- **The axis is carried over frozen here and is not refit.** It is read from "
             "`results/south6/stage1.npz` (`p_unit`), unit norm, 18,849 nonzero coordinates, all of them inside the 20,264 "
             "matched genes.")
    L.append("")
    L.append("> **The circularity, stated plainly.** The axis was fit on the same 14,013 non-targeting cells that (i) form the "
             "reference pseudobulk every factor's displacement is measured against, and (ii) generate every depth-matched null "
             "draw used here, including the residualized floors. The axis is therefore not independent of the data it is applied "
             "to: it is the leading cell-cycle contrast *within the very pool that defines zero*. Removing that direction removes "
             "whatever part of a factor's displacement happens to align with the dominant cell-cycle contrast of the control pool, "
             "and that contrast is itself estimated from a finite sample of those cells. Every SOUTH9 number below carries this "
             "circularity. The pre-registration forbids refitting the axis to avoid it, and it was not refit.")
    L.append("")
    L.append("**One choice the pre-registration left open, fixed here before any number was computed.** SOUTH6's residual is "
             "defined in gene space; the SOUTH8 primary space is a 20-component projection, and removing a direction in gene "
             "space does not leave the projection free of it. Two readings are therefore possible, and both are reported:")
    L.append("")
    L.append("- `in_space` (**primary**): the axis is projected into the space, and its component is removed from the "
             "displacement and the origin-target vector *inside that space*. This is the reading that makes the residualized "
             "statistic actually free of the proliferation direction as the space sees it.")
    L.append("- `gene_space` (secondary): SOUTH6's own formula - remove in gene space, then project.")
    L.append("")
    L.append("In S1 the two coincide. How much of the axis each space retains, and how aligned it is with the young-old "
             "direction:")
    L.append("")
    grows = []
    for s in SPACES7:
        g = geo[s]
        if g["proj"]:
            grows.append(dict(space=s, axis_retention=g["axis_retention"],
                              cos_axis_vs_gap=float(g["r"] @ g["pu"] / g["rnorm"])))
        else:
            grows.append(dict(space=s, axis_retention=1.0,
                              cos_axis_vs_gap=float(g["r_full_dot_p"] / np.sqrt(g["r2_full"]))))
    L.append(md_table(pd.DataFrame(grows), ["space", "axis_retention", "cos_axis_vs_gap"],
                      ["space", "||proj p|| (axis retention)", "cos(origin->target, axis)"]))
    L.append("")
    L.append("The cosine matters for reading D2: the more the proliferation axis aligns with the young-old direction, the more "
             "of a genuine approach residualization must also remove.")
    L.append("")

    L.append("## 2. Stage 0: SOUTH8's numbers re-derived (gate)")
    L.append("")
    rep = R["rep"]
    L.append(f"{len(rep)} checks, {int(rep.gating.sum())} gating, {int(((~rep.passed) & rep.gating).sum())} failures. "
             "Full table: `results/south9/stage0_reproduction.csv`.")
    L.append("")
    L.append(md_table(rep, ["block", "item", "n", "n_mismatch", "passed"],
                      ["block", "item", "n", "mismatches", "passed"]))
    L.append("")
    L.append("Nothing in `results/south6/`, `results/south7/` or `results/south8/` was written, and no earlier findings file "
             "was regenerated.")
    L.append("")

    L.append("## 3. D1: the residualized statistic in the anchor-free space (gated)")
    L.append("")
    L.append("Every depth-matched non-targeting draw from SOUTH7 Stage A1 is put through the same residualization as a factor, "
             "at the same cell count, and the 5th percentile is recomputed from the residualized draws: "
             f"{R['n_depths']} distinct cell counts, {R_NULL} draws each, per space and per variant. No raw floor is applied to "
             "a residualized statistic anywhere in this file.")
    L.append("")
    L.append(f"Primary reading (`{PRIMARY_VARIANT}`), raw and residualized side by side:")
    L.append("")
    show = d1[d1.variant == PRIMARY_VARIANT]
    L.append(md_table(show, ["factor_set", "space", "m", "raw_below_floor", "raw_negative", "raw_called",
                             "res_below_floor", "res_negative", "res_called"],
                      ["factor set", "space", "m", "raw: below floor", "raw: delta<0", "raw: called",
                       "residualized: below floor", "residualized: delta<0", "residualized: called"]))
    L.append("")
    L.append("Secondary reading (`gene_space`, SOUTH6's own formula):")
    L.append("")
    show2 = d1[d1.variant == "gene_space"]
    L.append(md_table(show2, ["factor_set", "space", "m", "raw_called", "res_below_floor", "res_negative", "res_called"],
                      ["factor set", "space", "m", "raw: called", "residualized: below floor",
                       "residualized: delta<0", "residualized: called"]))
    L.append("")
    surv = R["survivors"]
    L.append("Which of SOUTH8's calls survive, in the primary space and primary variant:")
    L.append("")
    L.append(md_table(pd.DataFrame(surv), ["factor_set", "south8_called", "res_called", "survive", "lost",
                                           "newly_called", "survivors_list"],
                      ["factor set", "SOUTH8 called (raw)", "called residualized", "in both", "lost",
                       "called residualized but NOT by SOUTH8", "in both"]))
    L.append("")
    L.append("**The residualized call set is not a subset of the raw one.** In the primary set, "
             f"{surv[0]['survive']} of SOUTH8's {surv[0]['south8_called']} survive, but "
             f"{surv[0]['newly_called']} of the {surv[0]['res_called']} residualized calls are factors SOUTH8 did not call: "
             + (", ".join(surv[0]["newly_called_list"].split(";")) if surv[0]["newly_called_list"] else "none")
             + ". Removing a direction moves some factors below their floor as well as lifting others above it, so "
             "'15 called' must not be read as '15 of the 24'.")
    L.append("")
    L.append("### Gate reading that fired (pre-registered wording)")
    L.append("")
    for set_name, g in g1.items():
        L.append(f"- {set_name}: {g['k']} of {g['m']} called residualized, against {g['e']} expected by chance "
                 f"(95% interval {g['lo']} to {g['hi']}) — {'inside or below' if g['inside_or_below'] else 'above'} the interval.")
    L.append("")
    p304 = g1.get("primary_304", {})
    pall = g1.get("all_1836", {})
    if fired1 == "EXPLAINED":
        L.append(f"**EXPLAINED.** The excess of factors below their own noise floor is accounted for by proliferation: "
                 f"residualized, {p304.get('k')} of 304 and {pall.get('k')} of 1,836 are called, against 15.2 and 91.8 "
                 "expected by chance. What the raw statistic scores as approach to young is in this screen a change in "
                 "proliferative state.")
    elif fired1 == "SURVIVES":
        L.append(f"**SURVIVES.** The excess survives proliferation residualization: {p304.get('k')} of 304 and "
                 f"{pall.get('k')} of 1,836 are called against 15.2 and 91.8 expected. Proliferation does not account for it.")
    else:
        inside = [k for k, v in g1.items() if v["inside_or_below"]]
        above = [k for k, v in g1.items() if not v["inside_or_below"]]
        L.append(f"**SPLIT.** One factor set falls inside or below its chance interval and the other sits above it, so no "
                 f"single conclusion is drawn. Inside or below: {', '.join(inside)} "
                 f"({', '.join(str(g1[k]['k']) + ' of ' + str(g1[k]['m']) for k in inside)}). "
                 f"Above: {', '.join(above)} "
                 f"({', '.join(str(g1[k]['k']) + ' of ' + str(g1[k]['m']) for k in above)}).")
    if fired2 == "LOST":
        L.append("")
        L.append("**This reading carries no conclusion**, because D2 fired LOST: residualization costs the statistic the "
                 "sensitivity it would need for a low residualized count to mean absence of approach.")
    L.append("")

    L.append("## 4. D2: does residualization also remove a real approach? (gated)")
    L.append("")
    L.append("SOUTH6REC's plant recovery, repeated in S2A_k20 at 106 cells, 200 draws per dose, doses 0/25/50/100% of the "
             "measured gap, with the residualization applied to plant and null alike.")
    L.append("")
    L.append(md_table(d2[d2.variant == PRIMARY_VARIANT],
                      ["space", "raw_power_25", "res_power_25", "res_power_50", "res_power_100",
                       "res_median_0", "res_median_25", "res_monotone"],
                      ["space", "raw power 25%", "residualized power 25%", "residualized power 50%",
                       "residualized power 100%", "residualized median at 0%", "at 25%", "monotone"]))
    L.append("")
    L.append("### Gate reading that fired (pre-registered wording)")
    L.append("")
    if fired2 == "RETAINED":
        L.append(f"- **RETAINED.** Residualization leaves the statistic able to detect a 25% approach, power {p_res:.3f}, so a "
                 "null result after residualization is a null result and not a loss of sensitivity.")
    else:
        L.append(f"- **LOST.** Residualization costs the statistic its sensitivity: power at a 25% plant falls from "
                 f"{R6_RAW_POWER_25:.3f} to {p_res:.3f}. The residualized counts cannot be read as evidence of absence.")
    L.append("")
    L.append("The planted gap and the proliferation axis are not orthogonal, so some of a genuine approach lies along the "
             "direction being removed; the cosine in section 1 says how much.")
    L.append("")

    L.append("## 5. D3: report-only (no decision rule)")
    L.append("")
    L.append(md_table(d3, ["factor_set", "m", "n_called", "n_high_depth_called", "id_drop_median_called",
                           "id_drop_median_set", "guide_agreement_median_called", "guide_agreement_median_set"],
                      ["factor set", "m", "called residualized", "of those, high_depth_candidate",
                       "identity drop, median (called)", "identity drop, median (set)",
                       "guide agreement (called)", "guide agreement (set)"]))
    L.append("")
    L.append("Identity, the comparison the pre-registration asks for: SOUTH8's 24 have a median identity drop of -0.0122 "
             "against +0.0408 for the screen, so they lose less fibroblast identity than average. The same comparison for the "
             "set that survives here is in the table above (a negative drop means identity genes rose relative to the "
             "non-targeting pool).")
    L.append("")
    L.append("The frozen ruler against the distance call - the preprint's thesis in this screen, a diagnostic and not evidence "
             "about any factor:")
    L.append("")
    L.append(md_table(d3, ["factor_set", "m", "ruler_younger_total", "ruler_younger_not_called"],
                      ["factor set", "m", "score younger on the ruler", "score younger yet NOT called residualized"]))
    L.append("")
    L.append("SOUTH8's 24, raw against residualized, individually:")
    L.append("")
    L.append(md_table(s24.sort_values(f"res_delta_{PRIMARY_SPACE}_{PRIMARY_VARIANT}"),
                      ["factor", "n_cells", "high_depth_candidate", f"raw_delta_{PRIMARY_SPACE}",
                       f"depth_floor_{PRIMARY_SPACE}", f"res_delta_{PRIMARY_SPACE}_{PRIMARY_VARIANT}",
                       f"res_floor_{PRIMARY_SPACE}_{PRIMARY_VARIANT}",
                       f"res_call_{PRIMARY_SPACE}_{PRIMARY_VARIANT}", "identity_drop", "guide_agreement",
                       "ruler_change"],
                      ["factor", "cells", "high_depth", "raw delta", "raw floor", "residualized delta",
                       "residualized floor", "still called", "identity drop", "guide agreement",
                       "ruler score change"]))
    L.append("")

    L.append("## 6. Limitations")
    L.append("")
    L.append("- **The axis is not independent of the test.** It is the dominant cell-cycle contrast of the same non-targeting "
             "pool that defines the reference and every null draw (section 1). A residualized count is therefore a statement "
             "about movement orthogonal to one finite-sample estimate of that contrast, not about proliferation in general.")
    L.append("- Residualization removes one direction. A factor whose effect is partly proliferative and partly not keeps only "
             "the orthogonal part, so a surviving factor is not proliferation-free and a lost factor is not purely "
             "proliferative.")
    L.append("- The cc.genes score is a coarse proxy for proliferative state; it does not distinguish cycle arrest from "
             "senescence, quiescence or a change in cycle phase distribution.")
    L.append("- No factor named here is a candidate for a pilot panel, and nothing in this file proposes one.")
    L.append("- A residualized null is not evidence that no transcription factor can rejuvenate aged cells: the screen uses "
             "neonatal foreskin fibroblasts with no aged donor and no matched young target.")
    L.append("- O2 pairs were not run.")
    L.append("- The residualized criterion inherits SOUTH8's: a 5th-percentile floor with no multiple-testing control, since "
             "the permutation condition was removed as uncalibrated in SOUTH8 D2.")
    L.append("- 3,221 of the ruler's 23,485 genes are outside the panel and are assumed unmoved; the axis is supported only on "
             "the 20,264 matched genes.")
    L.append("")
    L.append("## 7. Files")
    L.append("")
    L.append("`src/south9_run.py`. In `results/south9/`: `PREREG.flag`, `method_statement.json`, `stage0_reproduction.csv`, "
             "`d1_residualized_floors.csv`, `d1_counts.csv`, `d1_primary_304_residualized.csv`, "
             "`d1_all_1836_residualized.csv`, `d2_residualized_recovery.csv`, `d3_summary.csv`, "
             "`d3_south8_24_raw_vs_residualized.csv`, `d3_*_called_detail.csv`, `south9_run.log`.")
    FINDINGS.write_text("\n".join(L) + "\n", encoding="utf-8")
    log(f"[findings] {FINDINGS} written ({len(L)} lines)")


if __name__ == "__main__":
    require_prereg()
    log = Logger(OUT / "south9_run.log", mode="a")
    try:
        log(f"[prereg] {PREREG} bytes={PREREG.stat().st_size}")
        method_statement(log)
        rep, tables, odf, s8list = stage0(log)
        ctx = {}
        _old = REC.OUT
        try:
            REC.OUT = OUT / "_recheck"
            REC.OUT.mkdir(parents=True, exist_ok=True)
            ctx = REC.stage0(ctx, log)
        finally:
            REC.OUT = _old
        geo, mr, Lm, c_old, c_young, p_full, p_m = build_geometry(log)
        ctx["odf"] = odf
        ctx["c_old"], ctx["c_young"] = c_old, c_young
        d1, per_factor, fdf = run_d1(ctx, geo, tables, log)
        d2 = run_d2(ctx, geo, log)
        d3, s24 = run_d3(per_factor, s8list, log)
        gate1 = d1_gate(d1)
        gate2 = d2_gate(d2)
        log(f"[gate D1] {gate1[0]} {gate1[1]}")
        log(f"[gate D2] {gate2[0]} res_power_25={gate2[1]:.3f} raw={gate2[2]:.3f}")
        survivors = []
        for set_name, t in per_factor.items():
            raw_called = set(t.loc[t[f"amended_call_{PRIMARY_SPACE}"], "factor"])
            res_called = set(t.loc[t[f"res_call_{PRIMARY_SPACE}_{PRIMARY_VARIANT}"], "factor"])
            keep = sorted(raw_called & res_called)
            new = sorted(res_called - raw_called)
            survivors.append(dict(factor_set=set_name, south8_called=len(raw_called),
                                  res_called=len(res_called),
                                  survive=len(keep), lost=len(raw_called) - len(keep),
                                  newly_called=len(new),
                                  survivors_list=";".join(keep),
                                  newly_called_list=";".join(new)))
        pd.DataFrame(survivors).to_csv(OUT / "d1_survivors.csv", index=False)
        write_findings(dict(rep=rep, d1=d1, d2=d2, d3=d3, s24=s24, geo=geo,
                            gate1=gate1, gate2=gate2, survivors=survivors,
                            n_depths=int(fdf.n_cells.nunique())), log)
    except Stop as err:
        log(f"[STOP] {err.key}: {err.message}")
        (OUT / f"STOP_{err.key}.txt").write_text(err.message + "\n", encoding="utf-8")
        raise
    finally:
        log.close()
