"""DEPTH — at what per-sample cell count does the mixing-series positive control
still detect a young spike-in?

SAME's positive control P (src/same_run.py run_posctrl) mixed aged and young
half-B day-0 Fibroblast cells at f = 0, 0.10, 0.25, 0.50 with n_total = 2,470.
This task repeats P's mixtures at n_total = 400, 1200 and 2470, 200 independent
draws each, and asks whether a 10% / 25% young spike-in is still detected
against the f = 0 draws.

Pre-registration: results/depth/PREREG.flag (written before any DEPTH number).
This script never rewrites it. Instruments are read, not refit. src/same_run.py
and src/toward_run.py are imported, not modified; same_run.SAME_DIR is
redirected to results/depth/stage0/ for the Stage 0 re-run of P only.

Inputs: data/processed/md4/louvain_rulerY_<line>.npz and cell_scores_<line>.npz
are git-ignored; they are regenerated from the GSE297234 download by the MD4
scripts (REPO_INVENTORY.md). The gates are re-verifiable from the tracked
results/depth/draws.csv alone: python src/depth_run.py verify

Usage:
  python src/depth_run.py          Stage 0, draws, tables, readings, findings
  python src/depth_run.py verify   recompute tables and readings from draws.csv
"""
from __future__ import annotations

import sys
from pathlib import Path

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parent))
from config import RESULTS, ROOT  # noqa: E402
from brain_phase1_common import Logger, dump_json, load_json  # noqa: E402
from gtex_common import StopStep  # noqa: E402
from fibro_common import FROZEN_RULER, jsonable  # noqa: E402
from fibro2_common import load_frozen_ruler  # noqa: E402
from fibro3_common import AGED_LINE, YOUNG_LINE  # noqa: E402
from target_common import md_table  # noqa: E402
import same_run  # noqa: E402
from same_run import (  # noqa: E402
    F_VALUES,
    SAME_SEED,
    axis_stats,
    d0_fibroblast_idx,
    mix_sum,
    monotonic_increasing,
    panel_z,
    run_posctrl,
    split_half,
)
from toward_run import load_donor, sum_rows  # noqa: E402


DEPTH_DIR = RESULTS / "depth"
STAGE0_DIR = DEPTH_DIR / "stage0"
SAME_STORED = RESULTS / "same"
FINDINGS_PATH = ROOT / "FINDINGS_DEPTH.md"
PREREG_FLAG = DEPTH_DIR / "PREREG.flag"

DEPTH_SEED = 20260914
N_GRID = (400, 1200, 2470)
GATED_DEPTHS = (400, 1200)
N_DRAWS = 200
FLOOR_PCT = 5.0
RATE_BAR = 0.80
TOL_REL = 1e-4

READING_TEXT = {
    "USABLE": (
        "At {N} cells per sample, this depth preserves the sensitivity the mixing "
        "control was validated at."
    ),
    "MARGINAL": (
        "At {N} cells per sample, the mixing control would detect a quarter-gap "
        "approach (a 25% young spike-in) but not a tenth-gap one (a 10% young "
        "spike-in). An effect smaller than a quarter of the donor gap would be "
        "undetectable by construction at this depth."
    ),
    "INSUFFICIENT": (
        "At {N} cells per sample, the mixing control does not work: even a 25% young "
        "spike-in is not detected in 80% of draws. A null result from an experiment "
        "at this depth could not be read."
    ),
    "NO_READING": (
        "At {N} cells per sample, no pre-registered reading fired: both spike-in rungs "
        "were detected in at least 80% of draws, but frac was monotonic in f in fewer "
        "than 80% of draws. No conclusion about this depth's sensitivity is drawn."
    ),
}
RESCUE_SENTENCE = (
    "A depth that fails here cannot be rescued by analysis. The reading is about "
    "what the experiment could see, not about biology."
)


def load_prereg_text():
    if not PREREG_FLAG.exists():
        raise StopStep("prereg", "results/depth/PREREG.flag missing. Not computing.")
    return PREREG_FLAG.read_text(encoding="utf-8").rstrip("\n")


def n_split(n_total, f):
    n_young = int(np.round(float(f) * int(n_total)))
    return int(n_total) - n_young, n_young


def close_enough(new, stored):
    new, stored = float(new), float(stored)
    if not (np.isfinite(new) and np.isfinite(stored)):
        return False
    return abs(new - stored) <= TOL_REL * max(1.0, abs(stored))


# --------------------------------------------------------------------------- stage 0
def stage0(Y_a, Y_y, obs_a, obs_y, frozen, log):
    rows = []

    def add(item, stored, new, kind, gated=True):
        if kind == "exact":
            ok = bool(stored == new)
        else:
            ok = close_enough(new, stored)
        rows.append(dict(item=item, kind=kind, gated=gated, stored=stored, reproduced=new,
                         abs_diff=(abs(float(new) - float(stored)) if kind == "tol" else np.nan),
                         ok=ok))

    # (a) split
    split_tab = pd.read_csv(SAME_STORED / "split_table.csv")
    stored_idx = np.load(SAME_STORED / "split_idx.npz")
    d0_a = d0_fibroblast_idx(obs_a)
    d0_y = d0_fibroblast_idx(obs_y)
    o_new, ab_new = split_half(d0_a, SAME_SEED)
    y_new, yb_new = split_half(d0_y, SAME_SEED)
    new_sizes = {(AGED_LINE, "A"): o_new.size, (AGED_LINE, "B"): ab_new.size,
                 (YOUNG_LINE, "A"): y_new.size, (YOUNG_LINE, "B"): yb_new.size}
    for _, r in split_tab.iterrows():
        add(f"split_table {r.donor} half {r.half} n_cells", int(r.n_cells),
            int(new_sizes[(str(r.donor), str(r.half))]), "exact")
    for key, arr in (("O", o_new), ("aged_B", ab_new), ("Y", y_new), ("young_B", yb_new)):
        st = np.asarray(stored_idx[key], int)
        add(f"split_idx.npz {key} indices identical (extra check)", True,
            bool(st.shape == arr.shape and np.array_equal(st, arr)), "exact")

    O = np.asarray(stored_idx["O"], int)
    Yi = np.asarray(stored_idx["Y"], int)
    aged_B = np.asarray(stored_idx["aged_B"], int)
    young_B = np.asarray(stored_idx["young_B"], int)

    # (b) P re-run, output redirected away from results/same/
    STAGE0_DIR.mkdir(parents=True, exist_ok=True)
    saved_dir = same_run.SAME_DIR
    same_run.SAME_DIR = STAGE0_DIR
    try:
        log(f"[stage0] same_run.SAME_DIR redirected {saved_dir} -> {STAGE0_DIR}")
        pdf_new, psum_new = run_posctrl(Y_a, Y_y, O, Yi, aged_B, young_B, frozen, log)
    finally:
        same_run.SAME_DIR = saved_dir
    pdf_old = pd.read_csv(SAME_STORED / "posctrl.csv")
    psum_old = load_json(SAME_STORED / "posctrl_summary.json")
    for _, r_old in pdf_old.iterrows():
        m = np.isclose(pdf_new["f"].to_numpy(float), float(r_old.f), atol=1e-12)
        if m.sum() != 1:
            add(f"posctrl f={r_old.f} row present", True, False, "exact")
            continue
        r_new = pdf_new[m].iloc[0]
        tag = f"posctrl f={float(r_old.f):.2f}"
        for c in ("n_aged", "n_young", "n_total"):
            add(f"{tag} {c}", int(r_old[c]), int(r_new[c]), "exact")
        for c in ("frac", "delta", "cos"):
            add(f"{tag} {c}", float(r_old[c]), float(r_new[c]), "tol")
        for c in ("frac_ci_lo", "frac_ci_hi", "delta_ci_lo", "delta_ci_hi"):
            add(f"{tag} {c} (diagnostic)", float(r_old[c]), float(r_new[c]), "tol", gated=False)
    for c in ("pass", "monotonic", "frac50_ci_excludes_0", "delta50_neg"):
        add(f"posctrl_summary {c}", bool(psum_old[c]), bool(psum_new[c]), "exact")
    add("posctrl_summary n_total", int(psum_old["n_total"]), int(psum_new["n_total"]), "exact")

    df = pd.DataFrame(rows)
    df.to_csv(DEPTH_DIR / "stage0_reproduction.csv", index=False)
    bad = df[df.gated & ~df.ok]
    log(f"[stage0] checks={len(df)} gated={int(df.gated.sum())} gated_failures={len(bad)} "
        f"diagnostic_failures={int((~df.gated & ~df.ok).sum())}")
    if len(bad):
        names = "; ".join(f"{r.item}: stored={r.stored} reproduced={r.reproduced}"
                          for r in bad.itertuples())
        raise StopStep("prior_results_not_reproduced", names,
                       dict(mismatches=bad.item.tolist()))
    return dict(O=O, Y=Yi, aged_B=aged_B, young_B=young_B)


# --------------------------------------------------------------------------- draws
def check_pools(aged_B, young_B):
    for n in N_GRID:
        for f in F_VALUES:
            n_aged, n_young = n_split(n, f)
            if n_aged > aged_B.size or n_young > young_B.size:
                raise StopStep(
                    "pool_too_small",
                    f"n_total={n} f={f} needs n_aged={n_aged} n_young={n_young}; "
                    f"half-B pools are {aged_B.size}/{young_B.size}.",
                )


def draw_cells(n_index, f_index, draw_index, n_total, f, aged_B, young_B):
    rng = np.random.default_rng([DEPTH_SEED, n_index, f_index, draw_index])
    n_aged, n_young = n_split(n_total, f)
    take_a = rng.choice(aged_B, size=n_aged, replace=False) if n_aged > 0 else np.array([], int)
    take_y = rng.choice(young_B, size=n_young, replace=False) if n_young > 0 else np.array([], int)
    return n_aged, n_young, np.asarray(take_a, int), np.asarray(take_y, int)


def run_draws(Y_a, Y_y, idx, frozen, log):
    check_pools(idx["aged_B"], idx["young_B"])
    c_O = sum_rows(Y_a, idx["O"])
    c_Y = sum_rows(Y_y, idx["Y"])
    rows = []
    for n_index, n_total in enumerate(N_GRID):
        for draw in range(N_DRAWS):
            counts = [c_O, c_Y]
            meta = []
            for f_index, f in enumerate(F_VALUES):
                n_aged, n_young, ta, ty = draw_cells(
                    n_index, f_index, draw, n_total, f, idx["aged_B"], idx["young_B"])
                if ta.size != n_aged or ty.size != n_young or np.unique(ta).size != n_aged \
                        or np.unique(ty).size != n_young:
                    raise StopStep("draw", f"n={n_total} f={f} draw={draw}: bad sample sizes.")
                counts.append(mix_sum(Y_a, Y_y, ta, ty))
                meta.append((f_index, f, n_aged, n_young))
            Z, missing, _, _ = panel_z(counts, frozen)
            for i, (f_index, f, n_aged, n_young) in enumerate(meta):
                st = axis_stats(Z[0], Z[1], Z[2 + i])
                rows.append(dict(
                    n_index=n_index, n_total=n_total, f_index=f_index, f=f, draw=draw,
                    n_aged=n_aged, n_young=n_young,
                    frac=st["frac"], delta=st["delta"], cos=st["cos"],
                    d_norm=st["d_norm"], v_norm=st["v_norm"],
                    n_missing_z0=int(np.asarray(missing).sum()),
                ))
            if (draw + 1) % 50 == 0:
                log(f"[draws] n_total={n_total} {draw + 1}/{N_DRAWS}")
    df = pd.DataFrame(rows)
    df.to_csv(DEPTH_DIR / "draws.csv", index=False)
    return df


# --------------------------------------------------------------------------- tables
def tables_from_draws(df):
    det_rows, mono_rows, readings = [], [], {}
    for n_total in N_GRID:
        dn = df[df.n_total == n_total]
        d0 = dn[dn.f_index == 0]["delta"].to_numpy(float)
        if d0.size != N_DRAWS or not np.isfinite(d0).all():
            raise StopStep("floor", f"n_total={n_total}: f=0 deltas n={d0.size}, "
                                    f"finite={int(np.isfinite(d0).sum())}.")
        p5 = float(np.percentile(d0, FLOOR_PCT))
        rate = {}
        for f_index, f in enumerate(F_VALUES):
            g = dn[dn.f_index == f_index].sort_values("draw")
            dl = g["delta"].to_numpy(float)
            fr = g["frac"].to_numpy(float)
            det = np.isfinite(dl) & (dl < p5) & (dl < 0)
            rate[f_index] = float(det.sum()) / N_DRAWS
            det_rows.append(dict(
                n_total=n_total, f=f,
                n_aged=int(g.n_aged.iloc[0]), n_young=int(g.n_young.iloc[0]),
                n_draws=int(len(g)),
                n_detected=int(det.sum()), detection_rate=rate[f_index],
                role=("chance baseline (f=0 against its own floor)" if f_index == 0
                      else "spike-in"),
                floor_p5_delta=p5,
                median_frac=float(np.median(fr)), median_delta=float(np.median(dl)),
                diag_p5_frac=float(np.percentile(fr, 5)), diag_p95_frac=float(np.percentile(fr, 95)),
                diag_p5_delta=float(np.percentile(dl, 5)), diag_p95_delta=float(np.percentile(dl, 95)),
                diag_median_cos=float(np.median(g["cos"].to_numpy(float))),
                diag_share_delta_lt_0=float(np.mean(dl < 0)),
                # Added after the first run's readings existed; not gated.
                diag_share_below_floor=float(np.mean(np.isfinite(dl) & (dl < p5))),
                diag_p95_frac_f0=float(np.percentile(
                    dn[dn.f_index == 0]["frac"].to_numpy(float), 95)),
                diag_share_frac_gt_f0_p95=float(np.mean(fr > np.percentile(
                    dn[dn.f_index == 0]["frac"].to_numpy(float), 95))),
            ))
        piv = dn.pivot(index="draw", columns="f_index", values="frac").sort_index()
        mono = np.array([monotonic_increasing(piv.loc[d, list(range(len(F_VALUES)))].to_numpy(float))
                         for d in piv.index])
        mono_share = float(mono.mean())
        mono_rows.append(dict(n_total=n_total, n_draws=int(mono.size),
                              n_monotonic=int(mono.sum()), monotonic_share=mono_share))
        det10, det25 = rate[1], rate[2]
        cond = dict(det10=det10, det25=det25, mono=mono_share,
                    det10_ge_080=det10 >= RATE_BAR, det25_ge_080=det25 >= RATE_BAR,
                    mono_ge_080=mono_share >= RATE_BAR)
        if det10 >= RATE_BAR and det25 >= RATE_BAR and mono_share >= RATE_BAR:
            key = "USABLE"
        elif det25 >= RATE_BAR and det10 < RATE_BAR:
            key = "MARGINAL"
        elif det25 < RATE_BAR:
            key = "INSUFFICIENT"
        else:
            key = "NO_READING"
        readings[str(n_total)] = dict(
            key=key if n_total in GATED_DEPTHS else None,
            rule_outcome=key,
            gated=n_total in GATED_DEPTHS,
            text=(READING_TEXT[key].format(N=f"{n_total:,}") if n_total in GATED_DEPTHS else
                  "Reference rung; no gated reading (PREREG.flag DESIGN)."),
            conditions=cond, floor_p5_delta=p5,
        )
    return pd.DataFrame(det_rows), pd.DataFrame(mono_rows), readings


# --------------------------------------------------------------------------- findings
def _num(x, d=3):
    return "NA" if x is None or not np.isfinite(float(x)) else f"{float(x):+.{d}f}"


def _pct(x):
    return f"{100.0 * float(x):.1f}%"


def write_findings(prereg, stop=None, s0=None, det=None, mono=None, readings=None):
    L = [
        "# FINDINGS_DEPTH — cell depth for the mixing-series positive control",
        "",
    ]
    if stop:
        L += [f"**Status: STOP `{stop['step']}`.** {stop['message']}", ""]
    else:
        keys = ", ".join(f"n_total={n}: `{readings[str(n)]['key']}`" for n in GATED_DEPTHS)
        L += [f"**Status: DONE.** Readings: {keys}.", ""]
    L += [
        "Reproduced by `python src/depth_run.py`; gates re-checked from tracked files by "
        "`python src/depth_run.py verify`.",
        "",
        "## Pre-registration (verbatim from `results/depth/PREREG.flag`)",
        "",
        "```",
        prereg,
        "```",
        "",
    ]
    if s0 is not None:
        L += [
            "## Stage 0 — reproduction of SAME's stored numbers",
            "",
            "Counts exact; numerics within 1e-4 x max(1, |stored|). CI-bound rows are "
            "diagnostic (not gated). Full table: `results/depth/stage0_reproduction.csv`.",
            "",
        ]
        show = s0.copy()
        show["stored"] = [str(v) for v in show["stored"]]
        show["reproduced"] = [str(v) for v in show["reproduced"]]
        show["abs_diff"] = ["" if not np.isfinite(v) else f"{v:.2e}" for v in show["abs_diff"]]
        L += [md_table(show[["item", "kind", "gated", "stored", "reproduced", "abs_diff", "ok"]]), ""]
        n_g = int(s0.gated.sum())
        n_gf = int((s0.gated & ~s0.ok).sum())
        n_df = int((~s0.gated & ~s0.ok).sum())
        L += [f"{len(s0)} checks, {n_g} gated, {n_gf} gated failures, "
              f"{n_df} diagnostic failures.", ""]
    if stop:
        L += _not_covered() + _deviations() + _files(stop=True)
        FINDINGS_PATH.write_text("\n".join(L) + "\n", encoding="utf-8")
        return
    L += ["## Pre-registered readings (only the reading that fired, per depth)", ""]
    for n in GATED_DEPTHS:
        r = readings[str(n)]
        c = r["conditions"]
        L += [
            f"### n_total = {n:,}: `{r['key']}`",
            "",
            f"> {r['text']}",
            ">",
            f"> {RESCUE_SENTENCE}",
            "",
            f"Conditions: detection at f=0.10 = {_pct(c['det10'])} (>= 80%: {c['det10_ge_080']}); "
            f"detection at f=0.25 = {_pct(c['det25'])} (>= 80%: {c['det25_ge_080']}); "
            f"monotonic share = {_pct(c['mono'])} (>= 80%: {c['mono_ge_080']}). "
            f"Chance baseline (f=0 draws against their own floor): "
            f"{_pct(det[(det.n_total == n) & (det.f == 0.0)].detection_rate.iloc[0])}.",
            "",
        ]
    L += [
        "## Detection table",
        "",
        "Per depth and f, over 200 draws. `floor_p5_delta` is the 5th percentile of that "
        "depth's f=0 deltas. A draw at f>0 is detected iff delta < floor_p5_delta and "
        "delta < 0. The f=0 row is the chance baseline: the same rule applied to the floor's "
        "own draws (in-sample, at most 5% by construction). The 2470 rung is a reference "
        "rung with no gated reading; its f=0 draws are all the same 2,470 cells, so its "
        "floor is not a sampling floor. File: `results/depth/detection.csv`.",
        "",
    ]
    show = det[["n_total", "f", "n_aged", "n_young", "n_detected", "detection_rate",
                "floor_p5_delta", "median_frac", "median_delta"]].copy()
    show["detection_rate"] = [_pct(v) for v in show["detection_rate"]]
    show["f"] = [f"{v:.2f}" for v in show["f"]]
    L += [md_table(show), ""]
    L += [
        "Monotonicity: share of the 200 draw indices in which frac(f=0) <= frac(0.10) <= "
        "frac(0.25) <= frac(0.50) (`same_run.monotonic_increasing`). Cells are not shared "
        "across f within a draw index. File: `results/depth/monotonic.csv`.",
        "",
    ]
    ms = mono.copy()
    ms["monotonic_share"] = [_pct(v) for v in ms["monotonic_share"]]
    L += [md_table(ms), ""]
    ref = readings["2470"]
    L += [
        "## Diagnostics (not gated)",
        "",
        f"- 2470 reference rung: the rule, if applied, would give `{ref['rule_outcome']}` "
        f"(det10={_pct(ref['conditions']['det10'])}, det25={_pct(ref['conditions']['det25'])}, "
        f"mono={_pct(ref['conditions']['mono'])}). Not a reading: the f=0 floor at this depth "
        "has no sampling spread (every draw is the full aged pool), so the floor is "
        "narrower than a real replicate would be.",
        "- Spread of frac and delta over draws, median cos and the share of draws with "
        "delta < 0, per depth and f:",
        "",
    ]
    dg = det[["n_total", "f", "diag_p5_frac", "diag_p95_frac", "diag_p5_delta",
              "diag_p95_delta", "diag_median_cos", "diag_share_delta_lt_0"]].copy()
    dg["f"] = [f"{v:.2f}" for v in dg["f"]]
    L += [md_table(dg), ""]
    L += [
        "## Which clause of the detection rule binds (diagnostic; added after the "
        "readings were computed)",
        "",
        "This section was added after the first run's readings existed. It changes no "
        "reading and is not a gate. It splits the detection rule into its two clauses: "
        "below the f=0 floor, and delta < 0. It also gives a frac-based separation that the "
        "rule does not use.",
        "",
    ]
    cl = det[det.f > 0][["n_total", "f", "detection_rate", "diag_share_below_floor",
                         "diag_share_delta_lt_0", "diag_p95_frac_f0",
                         "diag_share_frac_gt_f0_p95"]].copy()
    cl["f"] = [f"{v:.2f}" for v in cl["f"]]
    for c in ("detection_rate", "diag_share_below_floor", "diag_share_delta_lt_0",
              "diag_share_frac_gt_f0_p95"):
        cl[c] = [_pct(v) for v in cl[c]]
    L += [md_table(cl), ""]
    for n in GATED_DEPTHS:
        g = det[det.n_total == n].set_index("f")
        p5 = float(g.loc[0.0, "floor_p5_delta"])
        L.append(
            f"- n_total={n:,}: the f=0 floor's 5th-percentile delta is {_num(p5)}, and "
            f"{_pct(g.loc[0.0, 'diag_share_delta_lt_0'])} of f=0 draws have delta < 0. At "
            f"f=0.10, {_pct(g.loc[0.10, 'diag_share_below_floor'])} of draws fall below the "
            f"floor and {_pct(g.loc[0.10, 'diag_share_delta_lt_0'])} have delta < 0. At "
            f"f=0.25 the shares are {_pct(g.loc[0.25, 'diag_share_below_floor'])} and "
            f"{_pct(g.loc[0.25, 'diag_share_delta_lt_0'])}. frac exceeds the f=0 draws' "
            f"95th percentile ({_num(g.loc[0.0, 'diag_p95_frac_f0'])}) in "
            f"{_pct(g.loc[0.10, 'diag_share_frac_gt_f0_p95'])} of f=0.10 draws."
        )
    L += [
        "",
        "Plainly: wherever the floor clause holds but detection fails, the failure is "
        "the delta < 0 clause. A shallow S is farther from Y than the deep half-A origin O "
        "is, through sampling noise alone. The spike-in moves S toward Y, but not past O's "
        "distance. frac is a projection, so noise does not add to it the way it adds to a "
        "distance. But frac is not depth-invariant either: its f=0 median is "
        + ", ".join(f"{_num(det[(det.n_total == n) & (det.f == 0.0)].median_frac.iloc[0])} "
                    f"at {n:,}" for n in N_GRID)
        + ". This describes how the "
        "pre-registered rule behaves. It is not a re-reading. Swapping to a frac-based rule "
        "now would be the after-the-fact rescue the flag forbids. A future flag that wanted "
        "a depth-robust rule would have to fix it in advance, and would have to model "
        "anchors at the experiment's own depth, which this design keeps deep.",
        "",
    ]
    L += _not_covered() + _deviations() + _files(stop=False)
    FINDINGS_PATH.write_text("\n".join(L) + "\n", encoding="utf-8")


def _not_covered():
    return [
        "## What this does not cover",
        "",
        "- Two donors, one young (GM23815) and one aged (GM00731), from a reprogramming "
        "dataset, not the six donors of the planned pilot. This measures depth sensitivity "
        "only; it says nothing about donor-to-donor variation, which more donors and not "
        "more cells would address.",
        "- Cells here are subsampled from a deeply sequenced run. A CRO's 400 recovered "
        "cells may also differ in reads per cell, which this cannot model. The result is "
        "therefore optimistic for the low-depth rungs.",
        "- The spike-in is cells of a genuinely younger donor, which is the most favourable "
        "possible positive signal. A partial approach to young by an intervention is a "
        "weaker and less clean signal than mixing in real young cells.",
        "- Only the sample S is subsampled. The anchors O and Y stay at half-A depth "
        "(2,469 and 3,851 cells). In a real experiment at 400 or 1,200 cells per sample the "
        "references would also be shallow, so this is a second reason the low-depth rungs "
        "are optimistic.",
        "- Circularity, disclosed: the f=0 floor and the spike-in draws come from the same "
        "half-B pools and share cells across draws, and the floor is the 5th percentile of "
        "the same draws the chance baseline is computed on.",
        "- Day-0 fibroblasts of two cell lines only. Nothing here speaks to other cell types, "
        "to partially reprogrammed cells, or to donors without an aged counterpart.",
        "",
    ]


def _deviations():
    return [
        "## Deviations and decisions beyond the task text",
        "",
        "- **Reading overlap and gap.** The task text's USABLE and INSUFFICIENT could both "
        "hold (10% rung and monotonicity pass, 25% rung fails), and no reading held when both "
        "rungs passed but monotonicity failed. Asked in chat before the flag was written. The "
        "author chose to make USABLE also require the 25% rung and to add `NO_READING` for "
        "the gap. The USABLE condition is therefore stricter than the task text.",
        "- **Panel composition** (not fixed by the task text): each (depth, draw index) is "
        "one 6-row TMM panel [O, Y, four mixes], the shape P used. TMM coupling means the "
        "O and Y z-scores move slightly between panels, as in P.",
        "- **Cells across f:** not shared. Each (depth, f, draw) has its own generator.",
        "- **Stage 0 re-runs P's bootstrap** (B=200, seed 20260918), to reproduce the stored "
        "`frac50_ci_excludes_0` flag. No DEPTH reading uses a bootstrap.",
        "- **Output redirect:** `same_run.SAME_DIR` was set to `results/depth/stage0/` for "
        "the Stage 0 call of `same_run.run_posctrl`, which writes `posctrl.csv`, "
        "`posctrl_summary.json` and `posctrl_z.npz` there. Nothing was written into "
        "`results/same/`.",
        "- **Import side effects:** importing `same_run` and `toward_run` runs "
        "`mkdir(exist_ok=True)` on `results/same/` and `results/toward/`. Both already "
        "existed, so nothing was created and no file in them was written.",
        "- **Extra Stage 0 check:** the recomputed split index arrays were compared for "
        "exact equality with `results/same/split_idx.npz`, beyond the half sizes.",
        "- **Fixed in the flag, not in the task text:** numpy-default (linear) percentile "
        "for the floor; strict inequalities; a non-finite delta counts as not detected; "
        "monotonic means non-decreasing (P's `monotonic_increasing`); rates compared as "
        ">= 0.80; the 2470 rung gets no gated reading.",
        "- **Added diagnostics:** cos, d_norm, v_norm, n_missing_z0 per draw; 5th/95th "
        "percentiles of frac and delta; share of draws with delta < 0. None of these are "
        "gated.",
        "- **Added after the readings existed:** the columns `diag_share_below_floor`, "
        "`diag_p95_frac_f0` and `diag_share_frac_gt_f0_p95`, and the section on which "
        "clause binds. They were added after the first complete run had produced the "
        "readings. Then the whole script was re-run from scratch, so the committed code "
        "produces every committed file. The readings and every gated number of the re-run "
        "are identical to the first run.",
        "- **Not opened:** `results/transfer/NOT_A_RESULT_gate_failed/`.",
        "",
    ]


def _files(stop):
    L = [
        "## Files",
        "",
        "- `src/depth_run.py`",
        "- `results/depth/PREREG.flag`",
        "- `results/depth/stage0_reproduction.csv`",
        "- `results/depth/stage0/` (P re-run: posctrl.csv, posctrl_summary.json, posctrl_z.npz)",
        "- `results/depth/run.log`",
    ]
    if stop:
        L += ["- `results/depth/STOP.json`"]
    else:
        L += [
            "- `results/depth/draws.csv` (one row per draw: frac, delta, cos, ...)",
            "- `results/depth/detection.csv`",
            "- `results/depth/monotonic.csv`",
            "- `results/depth/readings.json`",
        ]
    L += ["- `FINDINGS_DEPTH.md`", ""]
    return L


# --------------------------------------------------------------------------- main
def run(log):
    prereg = load_prereg_text()
    log("=" * 90)
    log(f"DEPTH seed={DEPTH_SEED} N_GRID={N_GRID} f={F_VALUES} n_draws={N_DRAWS}")
    log("=" * 90)
    frozen = load_frozen_ruler()
    log(f"[frozen] {FROZEN_RULER} n_genes={len(np.asarray(frozen['mu']))}")
    obs_a, Y_a = load_donor(AGED_LINE, log)
    obs_y, Y_y = load_donor(YOUNG_LINE, log)
    same_run.check_gene_space(Y_a, frozen, AGED_LINE)
    same_run.check_gene_space(Y_y, frozen, YOUNG_LINE)
    s0 = None
    try:
        idx = stage0(Y_a, Y_y, obs_a, obs_y, frozen, log)
        s0 = pd.read_csv(DEPTH_DIR / "stage0_reproduction.csv")
        df = run_draws(Y_a, Y_y, idx, frozen, log)
        det, mono, readings = tables_from_draws(df)
    except StopStep as e:
        stop = dict(step=e.step, message=e.message, details=jsonable(e.details))
        dump_json(DEPTH_DIR / "STOP.json", stop)
        p = DEPTH_DIR / "stage0_reproduction.csv"
        s0 = pd.read_csv(p) if p.exists() else None
        write_findings(prereg, stop=stop, s0=s0)
        log(f"[STOP] {e.step}: {e.message}")
        return 1
    det.to_csv(DEPTH_DIR / "detection.csv", index=False)
    mono.to_csv(DEPTH_DIR / "monotonic.csv", index=False)
    dump_json(DEPTH_DIR / "readings.json", jsonable(readings))
    for n in N_GRID:
        r = readings[str(n)]
        log(f"[reading] n_total={n} key={r['key']} rule_outcome={r['rule_outcome']} {r['conditions']}")
    write_findings(prereg, s0=s0, det=det, mono=mono, readings=readings)
    log("[run] DONE")
    return 0


def verify():
    df = pd.read_csv(DEPTH_DIR / "draws.csv")
    det, mono, readings = tables_from_draws(df)
    det_s = pd.read_csv(DEPTH_DIR / "detection.csv")
    mono_s = pd.read_csv(DEPTH_DIR / "monotonic.csv")
    read_s = load_json(DEPTH_DIR / "readings.json")
    ok = True
    for c in ("n_detected", "n_aged", "n_young"):
        ok &= bool((det[c].to_numpy() == det_s[c].to_numpy()).all())
    for c in ("detection_rate", "floor_p5_delta", "median_frac", "median_delta"):
        ok &= bool(np.allclose(det[c].to_numpy(float), det_s[c].to_numpy(float), rtol=1e-9, atol=1e-12))
    ok &= bool((mono["n_monotonic"].to_numpy() == mono_s["n_monotonic"].to_numpy()).all())
    for n in N_GRID:
        ok &= readings[str(n)]["rule_outcome"] == read_s[str(n)]["rule_outcome"]
        ok &= readings[str(n)]["key"] == read_s[str(n)]["key"]
        print(f"n_total={n}: key={readings[str(n)]['key']} rule_outcome={readings[str(n)]['rule_outcome']}")
    print("verify:", "OK" if ok else "MISMATCH")
    return 0 if ok else 1


def main():
    if len(sys.argv) > 1 and sys.argv[1] == "verify":
        return verify()
    DEPTH_DIR.mkdir(parents=True, exist_ok=True)
    log = Logger(DEPTH_DIR / "run.log")
    try:
        return run(log)
    finally:
        log.close()


if __name__ == "__main__":
    sys.exit(main())
