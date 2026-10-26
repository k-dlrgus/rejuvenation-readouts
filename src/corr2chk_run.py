"""CORR2CHK — do paper 2's Southard and young-donor "move away" clauses survive a
size-matched no-effect baseline?

Pre-registration: results/corr2chk/PREREG.flag, written before any number here existed.

Two statements are checked, both of the form "these units move away from young donors":

  A1b  paper 2 para 31, "and move away from young donors in the Southard screen".
       Units: the four gene names the sentence is about, in the Southard screen.
       Nothing is drawn: the size-matched baseline is already stored, keyed to each
       factor's own cell count, in results/south7/stage_a1_depth_nulls.csv.

  C5   paper 2 para 30, "and all 80 move away from them".
       Units: the 80 perturbations called closer to early-passage cells.
       The baseline is re-derived by re-running AXANG's own D3b floor loop with one
       change: the median and the 95th percentile are recorded alongside the 5th
       percentile that AXANG kept.

The decision statistic for both is the pair (d - M(n), d >= P95(n)) where d is the
stored change in distance and M, P95 come from no-effect draws at the unit's own cell
count. FLOORCHK's quadrature length correction is deliberately NOT used: the statistic
is already a change in distance taken against draws of the same size, so subtracting a
noise length again would double-count it (flag, D1).

Writes results/corr2chk/ only. Reads, and does not modify, results/south6/,
results/south7/, results/south8/, results/axang/, results/seng2/, results/seng3/,
results/toward/ and the Sengstack h5ad.

Run from the repository root with the Anaconda interpreter:
    python -I src/corr2chk_run.py run
    python -I src/corr2chk_run.py verify
"""
from __future__ import annotations

import json
import sys
from pathlib import Path

import numpy as np
import pandas as pd

SRC = Path(__file__).resolve().parent
if str(SRC) not in sys.path:
    sys.path.insert(0, str(SRC))

import seng2_common as K  # noqa: E402
import seng3_run as R3  # noqa: E402
from south6_run import block_delta  # noqa: E402

ROOT = SRC.parent
OUT = ROOT / "results" / "corr2chk"
OUT.mkdir(parents=True, exist_ok=True)
SOUTH6, SOUTH7, SOUTH8 = (ROOT / "results" / x for x in ("south6", "south7", "south8"))
AXANG, SENG2, SENG3, TOWARD = (ROOT / "results" / x for x in
                               ("axang", "seng2", "seng3", "toward"))

KDIM = 20
SPACE = "S2A_k20"
R_REPS = 200                 # AXANG's K.N_NULL, and SOUTH7's R_NULL
PRIMARY_SEED = K.SEED        # 20260914, AXANG's own; see flag D7
STABILITY_SEED = 20261030    # fresh, gates nothing; see flag D7
NAMED4 = ("E2F3", "EZH2", "STAT3", "ZFX")

# Tolerances, flag D11.
TOL_CROSS_FILE = 1e-12       # same value stored in two files
TOL_REBUILD = 1e-9           # recomputed from stored vectors, or a reseeded loop at its seed


class Stop(Exception):
    def __init__(self, key, msg):
        super().__init__(f"STOP {key}: {msg}")
        self.key, self.msg = key, msg


# --------------------------------------------------------------------------- repro table

class Repro:
    def __init__(self):
        self.rows = []

    def num(self, check, item, expected, observed, tol=TOL_REBUILD, gating=True):
        e, o = float(expected), float(observed)
        scale = max(1.0, abs(e))
        ok = bool(np.isfinite(e) and np.isfinite(o) and abs(e - o) <= tol * scale)
        self.rows.append(dict(check=check, item=item, expected=f"{e:.10g}",
                              observed=f"{o:.10g}", tol=f"{tol:g}", match=ok,
                              gating=gating))
        return ok

    def cnt(self, check, item, expected, observed, gating=True):
        ok = bool(int(expected) == int(observed))
        self.rows.append(dict(check=check, item=item, expected=str(int(expected)),
                              observed=str(int(observed)), tol="exact", match=ok,
                              gating=gating))
        return ok

    def cat(self, check, item, expected, observed, gating=True):
        ok = bool(str(expected) == str(observed))
        self.rows.append(dict(check=check, item=item, expected=str(expected),
                              observed=str(observed), tol="exact", match=ok,
                              gating=gating))
        return ok

    def frame(self):
        return pd.DataFrame(self.rows)

    def failures(self):
        return [r for r in self.rows if r["gating"] and not r["match"]]


# --------------------------------------------------------------------------- the statistic

def verdict(d, med, p95):
    """Flag D2. SUPPORTS / INSIDE NOISE / CONTRADICTS, for a claim that d moves away."""
    if d >= p95:
        return "SUPPORTS"
    if d < med:
        return "CONTRADICTS"
    return "INSIDE_NOISE"


def tally(df):
    v = df["verdict"].value_counts()
    return dict(m=int(len(df)),
                supports=int(v.get("SUPPORTS", 0)),
                inside_noise=int(v.get("INSIDE_NOISE", 0)),
                contradicts=int(v.get("CONTRADICTS", 0)),
                exp_supports_chance=round(0.05 * len(df), 2),
                exp_contradicts_chance=round(0.50 * len(df), 2))


# --------------------------------------------------------------------------- A1b

def a1b(log, rp):
    """Southard, four named factors, S2A_k20. Nothing is drawn (flag D5)."""
    stored = pd.read_csv(SOUTH8 / "d3c_all_1836_umi_calls.csv")
    nulls = pd.read_csv(SOUTH7 / "stage_a1_depth_nulls.csv")
    nulls = nulls[nulls.space == SPACE]
    obs = pd.read_csv(SOUTH7 / "stage_a1_observed.csv")
    cpt = pd.read_csv(SOUTH6 / "stage0_cells_per_target.csv")
    curves = pd.read_csv(SOUTH8 / "d1_plant_curves_all.csv")
    curves = curves[curves.space == SPACE]

    # R4 -- each name appears exactly once
    for g in NAMED4:
        rp.cnt("R4 Southard table", f"rows for {g}", 1, int((stored.factor == g).sum()))

    # R1 -- recompute the change in distance from the stored displacement and basis
    bas = np.load(SOUTH7 / "basis.npz", allow_pickle=True)
    spaces = [str(x) for x in bas["spaces"]]
    si = spaces.index(SPACE)
    r20 = np.asarray(bas[f"r_{SPACE}"], np.float64)
    r2 = float(np.asarray(bas["r2"])[si])
    rnorm = float(np.asarray(bas["rnorm"])[si])
    rp.num("R1 basis", "||O1-Y1|| in S2A_k20 (stored rnorm vs ||r||)",
           rnorm, float(np.linalg.norm(r20)))
    L = np.asarray(bas["L_mid_full"], np.float64)[:, :KDIM]

    s1 = np.load(SOUTH6 / "stage1.npz", allow_pickle=True)
    keep = [str(x) for x in s1["keep"]]
    pos = {g: keep.index(g) for g in NAMED4}
    d_f = s1["d_f"]
    d20 = np.vstack([np.asarray(d_f[pos[g]], np.float64) for g in NAMED4]) @ L
    del d_f
    rebuilt, _, _, _ = block_delta(d20, r20, r2, rnorm)

    rows, cell_n = [], {}
    for i, g in enumerate(NAMED4):
        srow = stored[stored.factor == g].iloc[0]
        n = int(srow.n_cells)
        cell_n[g] = n
        rp.num(f"R1 {g}", "change in distance, rebuilt vs stored delta_S2A_k20",
               float(srow[f"delta_{SPACE}"]), float(rebuilt[i]))
        orow = obs[obs.factor == g]
        if len(orow):
            rp.num(f"R1 {g}", "stored delta in stage_a1_observed vs d3c table",
                   float(srow[f"delta_{SPACE}"]), float(orow.iloc[0][f"delta_{SPACE}"]),
                   tol=TOL_CROSS_FILE)
        # R3 -- cell count against SOUTH6 stage 0
        c0 = cpt[cpt.target.astype(str) == g]
        rp.cnt(f"R3 {g}", "rows in SOUTH6 stage 0 cells-per-target", 1, len(c0))
        if len(c0):
            rp.cnt(f"R3 {g}", "cells per target vs SOUTH6 stage 0", n,
                   int(c0.iloc[0].n_cells))
        # R2 -- the floor read here is the floor the published call was made against
        nrow = nulls[nulls.n_cells == n]
        if not len(nrow):
            raise Stop("reproduction_failed",
                       f"no stored depth null at {n} cells in {SPACE} for {g}")
        nrow = nrow.iloc[0]
        rp.num(f"R2 {g}", "depth floor, SOUTH7 nulls vs SOUTH8 table",
               float(srow[f"depth_floor_{SPACE}"]), float(nrow.floor_p05),
               tol=TOL_CROSS_FILE)
        rp.cnt(f"R2 {g}", "draws behind the baseline", R_REPS, int(nrow.R))
        crow = curves[curves.n_cells == n]
        pw = float(crow.iloc[0].power_25) if len(crow) else float("nan")
        d = float(srow[f"delta_{SPACE}"])
        rows.append(dict(statement="A1b", unit=g, n_cells=n, change_in_distance=d,
                         floor_p05=float(nrow.floor_p05),
                         null_median=float(nrow["median"]),
                         null_sd=float(nrow["sd"]), null_p95=float(nrow.p95),
                         d_minus_median=d - float(nrow["median"]),
                         below_floor=bool(d < float(nrow.floor_p05)),
                         power_25_at_this_depth=pw,
                         verdict=verdict(d, float(nrow["median"]), float(nrow.p95))))
    df = pd.DataFrame(rows)
    log(f"[a1b] {tally(df)}")
    return df, cell_n


# --------------------------------------------------------------------------- C5

def axang_prologue(log):
    """AXANG d3's prologue, repeated line for line from src/axang_run.py:449-486.

    axang_run.d3() is not called because it writes into results/axang/ (flag D6).
    R5 is the proof that this repetition is faithful.
    """
    axz = np.load(AXANG / "axes.npz", allow_pickle=True)
    L = np.asarray(axz["L"], float)
    anch = np.load(TOWARD / "anchors.npz", allow_pickle=True)
    cy, co = np.asarray(anch["c_young"], float), np.asarray(anch["c_old"], float)
    r_full = co - cy                      # origin minus target, SOUTH's convention
    r20 = r_full @ L
    gap20 = float(np.linalg.norm((cy - co) @ L))
    sp3 = R3.space3(log, need_extra=False)
    ctx, zf, Z3, ni3 = sp3["ctx"], sp3["zf"], sp3["Z"], sp3["ni"]
    pc = np.load(SENG2 / "panel_counts.npz")
    dst, nr = np.asarray(pc["dst"], int), sp3["n_ruler"]
    C = pc["C"]
    blk, pool = {}, {}
    for m in ("CRA", "CRI"):
        blk[m] = np.asarray(K.submatrix(ctx["f"], ctx["nt_mod"][m], ctx["exp"]).todense(),
                            np.float64)
        pool[m] = C[ni3[f"NT_{m}"]]
        log(f"[c5] NT {m} {blk[m].shape}")

    def draws(m, n, reps=R_REPS, seed=PRIMARY_SEED):
        rng = np.random.default_rng(seed)
        D = []
        for _ in range(reps):
            t = rng.choice(blk[m].shape[0], size=n, replace=False)
            sub = K.to_ruler(blk[m][t].sum(axis=0), dst, nr)
            D.append(zf(sub) - zf(pool[m] - sub))
        return np.vstack(D) @ L

    def delta(d20):
        r2 = float(r20 @ r20)
        rn = float(np.sqrt(r2))
        dots = d20 @ r20
        qn2 = np.einsum("ij,ij->i", d20, d20)
        return np.sqrt(np.maximum(r2 + 2.0 * dots + qn2, 0.0)) - rn

    return dict(L=L, ctx=ctx, Z3=Z3, ni3=ni3, draws=draws, delta=delta, gap20=gap20)


def c5(log, rp):
    P = axang_prologue(log)
    ctx, Z3, ni3, L = P["ctx"], P["Z3"], P["ni3"], P["L"]
    draws, delta = P["draws"], P["delta"]
    stored = pd.read_csv(AXANG / "d3b_table.csv")
    d4 = pd.read_csv(SENG3 / "d4_table.csv")

    # R7 -- the units
    rp.cnt("R7 units", "rows in seng3 d4_table", 340, len(d4))
    rp.cnt("R7 units", "rows in axang d3b_table", 340, len(stored))
    rp.cnt("R7 units", "called closer to early-passage cells (own-control start)",
           80, int(d4.called_delta_i.sum()))
    j = d4[["modality", "perturbation", "n_cells", "called_delta_i", "delta_i"]].merge(
        stored[["modality", "perturbation", "n_cells", "delta", "floor",
                "below_floor", "decreased", "called"]],
        on=["modality", "perturbation", "n_cells"], how="inner")
    rp.cnt("R7 units", "rows surviving the join", 340, len(j))

    # R8 -- the published counts and the four named changes in distance
    rp.cnt("R8 axang", "called closer to young donors, of 340", 5, int(stored.called.sum()))
    rp.cnt("R8 axang", "distance decreased, of 340", 17, int(stored.decreased.sum()))
    for lab, want in (("CRA_E2F3", 11.448), ("CRA_EZH2", 9.022),
                      ("CRI_STAT3", 5.327), ("CRI_ZFX", 15.336)):
        rp.num(f"R8 {lab}", "stored young-donor change in distance (3 dp)", want,
               round(float(stored[stored.perturbation == lab].iloc[0].delta), 3),
               tol=1e-9)

    # R6 -- recompute every stored young-donor change in distance
    worst = 0.0
    for m, lab, idx in ctx["qual"]:
        d_full = Z3[ni3[lab]] - Z3[ni3[f"NT_{m}"]]
        dl = float(delta((d_full @ L)[None, :])[0])
        srow = stored[(stored.modality == m) & (stored.perturbation == lab)].iloc[0]
        worst = max(worst, abs(dl - float(srow.delta)) / max(1.0, abs(float(srow.delta))))
    rp.num("R6 axang", "max relative error, 340 recomputed changes in distance",
           0.0, worst)

    # The re-derived baseline: AXANG's loop, plus the median and P95 (flag D6)
    by = {"CRA": [], "CRI": []}
    for m, lab, idx in ctx["qual"]:
        by[m].append((lab, int(idx.size)))
    base, stab = {}, {}
    for m in ("CRA", "CRI"):
        ns = sorted({n for _, n in by[m]})
        log(f"[c5] {m}: {len(ns)} distinct cell counts {ns[0]}-{ns[-1]}")
        for i, n in enumerate(ns):
            dl = delta(draws(m, n))
            base[(m, n)] = (float(np.percentile(dl, 5)), float(np.median(dl)),
                            float(np.percentile(dl, 95)), float(np.std(dl, ddof=1)))
            ds = delta(draws(m, n, seed=STABILITY_SEED))
            stab[(m, n)] = (float(np.percentile(ds, 5)), float(np.median(ds)),
                            float(np.percentile(ds, 95)))
            if (i + 1) % 20 == 0:
                log(f"[c5 baseline] {m} {i + 1}/{len(ns)}")

    # R5 -- the gate the ruling requires
    frows, worst5 = [], 0.0
    for (m, n), (f5, med, p95, sd) in sorted(base.items()):
        srow = stored[(stored.modality == m) & (stored.n_cells == n)]
        sf = float(srow.iloc[0].floor)
        worst5 = max(worst5, abs(f5 - sf) / max(1.0, abs(sf)))
        s5, smed, sp95 = stab[(m, n)]
        frows.append(dict(modality=m, n_cells=n, R=R_REPS,
                          floor_p05=f5, null_median=med, null_p95=p95, null_sd=sd,
                          stored_floor=sf, floor_abs_err=abs(f5 - sf),
                          stability_floor_p05=s5, stability_median=smed,
                          stability_p95=sp95))
    fdf = pd.DataFrame(frows)
    ok5 = rp.num("R5 C5 baseline", "max relative error, re-derived vs stored 5th percentile",
                 0.0, worst5)
    if not ok5:
        raise Stop("axang_floor_not_reproduced",
                   f"re-derived 5th percentiles differ from results/axang/d3b_table.csv "
                   f"by up to {worst5:.3g} relative; C5 is PAPER2-UNRESOLVED and no "
                   f"re-derived median is used.")

    # R9 -- descriptive: AXANG's own dose-0 medians at its four gate depths
    plant = pd.read_csv(AXANG / "d3a_plant.csv")
    for _, pr in plant[plant.dose == 0.0].iterrows():
        key = (str(pr.library), int(pr.n_cells))
        if key in base:
            rp.num(f"R9 {key[0]} {key[1]}", "dose-0 median, re-derived vs AXANG d3a",
                   float(pr.delta_median), base[key][1], tol=1e-9, gating=False)

    rows = []
    for r in j[j.called_delta_i].itertuples():
        f5, med, p95, sd = base[(r.modality, int(r.n_cells))]
        d = float(r.delta)
        rows.append(dict(statement="C5", unit=r.perturbation, modality=r.modality,
                         n_cells=int(r.n_cells), change_in_distance=d,
                         floor_p05=f5, null_median=med, null_sd=sd, null_p95=p95,
                         d_minus_median=d - med, below_floor=bool(d < f5),
                         passage_change_in_distance=float(r.delta_i),
                         verdict=verdict(d, med, p95)))
    df = pd.DataFrame(rows).sort_values("n_cells").reset_index(drop=True)
    log(f"[c5] {tally(df)}")
    return df, fdf


# --------------------------------------------------------------------------- readings

def readings_of(a1b_df, c5_df):
    L = []
    sp = "k20 (paper 2's 20-component space)"
    for sid, df, src in (("A1b", a1b_df, "paper 2 paragraph 31"),
                         ("C5", c5_df, "paper 2 paragraph 30")):
        t = tally(df)
        L.append(f"== {sid} ({src}) ==")
        nums = (f"{t['supports']} of {t['m']} units move away further than the 95th "
                f"percentile of no-effect draws at their own cell count "
                f"(about {t['exp_supports_chance']} expected by chance), "
                f"{t['inside_noise']} are inside the no-effect spread, and "
                f"{t['contradicts']} move away less far than the median no-effect draw "
                f"(about {t['exp_contradicts_chance']} expected by chance)")
        if t["supports"] == t["m"]:
            L.append(f"CONCLUSION  PAPER2-STANDS: {sid}, {sp}: the statement survives the "
                     f"size-matched floor; {nums}.")
        else:
            L.append(f"CONCLUSION  PAPER2-DEPTH: {sid}, {sp}: the statement was reached "
                     f"partly by sampling noise; {nums}. Corrected wording: see "
                     f"FINDINGS_CORR2CHK.md, section \"Corrected wording\".")
    L.append("== the 17 statements not checked ==")
    for sid in ("A1a", "A2", "A3", "A4", "B4", "B5", "C1", "C2", "C3", "C4", "C6",
                "C7", "C8", "C9", "C10", "C11", "C12"):
        L.append(f"CONCLUSION  PAPER2-NOT-AFFECTED: {sid}: its threshold was already "
                 f"matched to group size, so FLOORCHK's correction does not apply.")
    L.append("            (qualification, flag Part 3) C7 has no threshold of any kind; it "
             "is not affected because there is nothing to size-match.")
    L.append("            (qualification, flag Part 3) A2's quoted numbers are themselves "
             "size-matched positive-control powers.")
    L.append("== B1, B2, B3 ==")
    L.append("            no size-matched floor exists for these three length ratios; per "
             "the ruling they are recorded with the direction of their noise bias and no "
             "key fires for them.")
    return L


# --------------------------------------------------------------------------- drivers

def run():
    log = K.Log(OUT / "corr2chk_run.log")
    rp = Repro()
    try:
        log("[stage] A1b, from stored values only")
        a1b_df, _ = a1b(log, rp)
        log("[stage] C5, re-deriving AXANG's discarded medians")
        c5_df, c5_floor = c5(log, rp)
    except Stop as e:
        rp.frame().to_csv(OUT / "reproduction.csv", index=False)
        (OUT / "readings.txt").write_text(f"{e}\n", encoding="utf-8")
        log(str(e))
        log.close()
        raise

    fails = rp.failures()
    rp.frame().to_csv(OUT / "reproduction.csv", index=False)
    log(f"[repro] {len(rp.rows)} checks, {len(fails)} gating mismatches")
    if fails:
        names = "; ".join(f"{f['check']} / {f['item']}: expected {f['expected']}, "
                          f"observed {f['observed']}" for f in fails)
        e = Stop("reproduction_failed", names)
        (OUT / "readings.txt").write_text(f"{e}\n", encoding="utf-8")
        log(str(e))
        log.close()
        raise e

    a1b_df.to_csv(OUT / "a1b_units.csv", index=False)
    c5_df.to_csv(OUT / "c5_units.csv", index=False)
    c5_floor.to_csv(OUT / "c5_baseline.csv", index=False)
    lines = readings_of(a1b_df, c5_df)
    (OUT / "readings.txt").write_text("\n".join(lines) + "\n", encoding="utf-8")
    summ = dict(a1b=tally(a1b_df), c5=tally(c5_df),
                repro_checks=len(rp.rows), repro_mismatches=0,
                primary_seed=PRIMARY_SEED, stability_seed=STABILITY_SEED, reps=R_REPS,
                a1b_contradicting_units=sorted(
                    a1b_df[a1b_df.verdict == "CONTRADICTS"].unit.tolist()),
                c5_contradicting_units=sorted(
                    c5_df[c5_df.verdict == "CONTRADICTS"].unit.tolist()),
                a1b_power_25_range=[float(a1b_df.power_25_at_this_depth.min()),
                                    float(a1b_df.power_25_at_this_depth.max())])
    (OUT / "summary.json").write_text(json.dumps(summ, indent=2), encoding="utf-8")
    for ln in lines:
        log("  " + ln)
    log.close()


def verify():
    """Re-derive every reading from the tracked CSVs alone."""
    a1b_df = pd.read_csv(OUT / "a1b_units.csv")
    c5_df = pd.read_csv(OUT / "c5_units.csv")
    want = (OUT / "readings.txt").read_text(encoding="utf-8").splitlines()
    got = readings_of(a1b_df, c5_df)
    bad = [(i, a, b) for i, (a, b) in enumerate(zip(want, got)) if a != b]
    print(f"{len(got)} reading lines re-derived from the tracked CSVs; "
          f"{len(bad)} mismatches; line counts {len(want)} vs {len(got)}")
    for i, a, b in bad:
        print(f"  line {i}:\n    stored: {a}\n    rebuilt: {b}")
    rp = pd.read_csv(OUT / "reproduction.csv")
    g = rp[rp.gating]
    print(f"reproduction.csv: {len(rp)} checks, {len(g)} gating, "
          f"{int((~g.match).sum())} gating mismatches")
    for df, sid in ((a1b_df, "A1b"), (c5_df, "C5")):
        t = tally(df)
        n_lt = int((df.change_in_distance < df.null_median).sum())
        print(f"{sid}: {t}; units below their own no-effect median: {n_lt}")
    return 0 if not bad and int((~g.match).sum()) == 0 else 1


MODES = dict(run=run, verify=verify)

if __name__ == "__main__":
    if len(sys.argv) != 2 or sys.argv[1] not in MODES:
        print(f"usage: python -I {Path(__file__).name} {'|'.join(MODES)}")
        raise SystemExit(2)
    raise SystemExit(MODES[sys.argv[1]]() or 0)
