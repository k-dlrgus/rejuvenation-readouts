"""DEPTH2 — the DEPTH test with references at the experiment's own depth.

DEPTH (src/depth_run.py, commit 5255f87) subsampled only the sample S and kept
the anchors O and Y at half-A depth. DEPTH2 builds O, Y and S all at n_total
cells in every draw (O from GM00731 half A, Y from GM23815 half A, S from the
half-B pools exactly as DEPTH did) and applies DEPTH's detection rule and
readings unchanged, at n_total = 400 and 1,200.

Pre-registration: results/depth2/PREREG.flag (written before any DEPTH2
number). This script never rewrites it. Instruments are read, not refit.
src/depth_run.py, src/same_run.py and src/toward_run.py are imported, not
modified. Nothing is written into results/depth/ or results/same/.

Inputs: data/processed/md4/louvain_rulerY_<line>.npz and cell_scores_<line>.npz
are git-ignored; they are regenerated from the GSE297234 download by the MD4
scripts (REPO_INVENTORY.md). The gates are re-verifiable from tracked files
alone (results/depth2/draws.csv, results/depth2/stage0/draws_deep.csv and
results/depth/): python src/depth2_run.py verify

Usage:
  python src/depth2_run.py          Stage 0, draws, tables, readings, findings
  python src/depth2_run.py verify   recompute Stage 0 comparison, tables and
                                    readings from the tracked draw files
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
    split_half,
)
from toward_run import load_donor, sum_rows  # noqa: E402
import depth_run  # noqa: E402
from depth_run import (  # noqa: E402
    FLOOR_PCT,
    N_DRAWS,
    RATE_BAR,
    READING_TEXT,
    RESCUE_SENTENCE,
    TOL_REL,
    close_enough,
    n_split,
)


DEPTH2_DIR = RESULTS / "depth2"
STAGE0_DIR = DEPTH2_DIR / "stage0"
DEPTH_STORED = RESULTS / "depth"
SAME_STORED = RESULTS / "same"
FINDINGS_PATH = ROOT / "FINDINGS_DEPTH2.md"
PREREG_FLAG = DEPTH2_DIR / "PREREG.flag"

DEPTH2_SEED = 20260915
N_GRID = (400, 1200)  # n_index 0, 1 (DEPTH's mapping)
ROLE_O, ROLE_Y, ROLE_S = 0, 1, 2

# Columns of results/depth/detection.csv compared exactly in Stage 0.
EXACT_COLS = ("n_aged", "n_young", "n_draws", "n_detected", "role")


def load_prereg_text():
    if not PREREG_FLAG.exists():
        raise StopStep("prereg", "results/depth2/PREREG.flag missing. Not computing.")
    return PREREG_FLAG.read_text(encoding="utf-8").rstrip("\n")


def classify(det10, det25, mono):
    """DEPTH's nesting (results/depth/PREREG.flag READINGS)."""
    if det10 >= RATE_BAR and det25 >= RATE_BAR and mono >= RATE_BAR:
        return "USABLE"
    if det25 >= RATE_BAR and det10 < RATE_BAR:
        return "MARGINAL"
    if det25 < RATE_BAR:
        return "INSUFFICIENT"
    return "NO_READING"


def score_row(Z, n_index, n_total, f_index, f, draw, n_aged, n_young, missing):
    st = axis_stats(Z[0], Z[1], Z[2])
    return dict(
        n_index=n_index, n_total=n_total, f_index=f_index, f=f, draw=draw,
        n_aged=n_aged, n_young=n_young,
        frac=st["frac"], delta=st["delta"], cos=st["cos"],
        d_norm=st["d_norm"], v_norm=st["v_norm"],
        n_missing_z0=int(np.asarray(missing).sum()),
    )


# --------------------------------------------------------------------------- tables
def tables_from_draws(df, depths):
    """DEPTH's detection rule and table columns, for the given depths only."""
    det_rows, mono_rows, readings = [], [], {}
    for n_total in depths:
        dn = df[df.n_total == n_total]
        d0 = dn[dn.f_index == 0]["delta"].to_numpy(float)
        if d0.size != N_DRAWS or not np.isfinite(d0).all():
            raise StopStep("floor", f"n_total={n_total}: f=0 deltas n={d0.size}, "
                                    f"finite={int(np.isfinite(d0).sum())}.")
        p5 = float(np.percentile(d0, FLOOR_PCT))
        fr0 = dn[dn.f_index == 0]["frac"].to_numpy(float)
        p95_fr0 = float(np.percentile(fr0, 95))
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
                diag_share_below_floor=float(np.mean(np.isfinite(dl) & (dl < p5))),
                diag_p95_frac_f0=p95_fr0,
                diag_share_frac_gt_f0_p95=float(np.mean(fr > p95_fr0)),
                diag_median_v_norm=float(np.median(g["v_norm"].to_numpy(float))),
            ))
        piv = dn.pivot(index="draw", columns="f_index", values="frac").sort_index()
        mono = np.array([monotonic_increasing(piv.loc[d, list(range(len(F_VALUES)))].to_numpy(float))
                         for d in piv.index])
        mono_share = float(mono.mean())
        mono_rows.append(dict(n_total=n_total, n_draws=int(mono.size),
                              n_monotonic=int(mono.sum()), monotonic_share=mono_share))
        det10, det25 = rate[1], rate[2]
        key = classify(det10, det25, mono_share)
        readings[str(n_total)] = dict(
            key=key, rule_outcome=key, gated=True,
            text=READING_TEXT[key].format(N=f"{n_total:,}"),
            conditions=dict(det10=det10, det25=det25, mono=mono_share,
                            det10_ge_080=det10 >= RATE_BAR, det25_ge_080=det25 >= RATE_BAR,
                            mono_ge_080=mono_share >= RATE_BAR),
            floor_p5_delta=p5,
        )
    return pd.DataFrame(det_rows), pd.DataFrame(mono_rows), readings


# --------------------------------------------------------------------------- stage 0
def read_pools():
    stored = np.load(SAME_STORED / "split_idx.npz")
    return {k: np.asarray(stored[k], int) for k in ("O", "Y", "aged_B", "young_B")}


def draws_deep(Y_a, Y_y, idx, frozen, log):
    """DEPTH's configuration: deep half-A anchors, 6-row panel, DEPTH's seeds."""
    c_O = sum_rows(Y_a, idx["O"])
    c_Y = sum_rows(Y_y, idx["Y"])
    rows = []
    for n_index, n_total in enumerate(N_GRID):
        for draw in range(N_DRAWS):
            counts = [c_O, c_Y]
            meta = []
            for f_index, f in enumerate(F_VALUES):
                n_aged, n_young, ta, ty = depth_run.draw_cells(
                    n_index, f_index, draw, n_total, f, idx["aged_B"], idx["young_B"])
                counts.append(mix_sum(Y_a, Y_y, ta, ty))
                meta.append((f_index, f, n_aged, n_young))
            Z, missing, _, _ = panel_z(counts, frozen)
            for i, (f_index, f, n_aged, n_young) in enumerate(meta):
                Zi = Z[[0, 1, 2 + i]]
                rows.append(score_row(Zi, n_index, n_total, f_index, f, draw,
                                      n_aged, n_young, missing))
            if (draw + 1) % 50 == 0:
                log(f"[stage0 draws] n_total={n_total} {draw + 1}/{N_DRAWS}")
    return pd.DataFrame(rows)


def compare_stage0(deep_df, split_rows):
    """Gated comparison of DEPTH-configuration output with results/depth/."""
    rows = list(split_rows)

    def add(item, stored, new, kind, gated=True):
        ok = bool(stored == new) if kind == "exact" else close_enough(new, stored)
        rows.append(dict(item=item, kind=kind, gated=gated, stored=stored, reproduced=new,
                         abs_diff=(abs(float(new) - float(stored)) if kind == "tol" else np.nan),
                         ok=ok))

    det_new, mono_new, read_new = tables_from_draws(deep_df, N_GRID)
    det_old = pd.read_csv(DEPTH_STORED / "detection.csv")
    det_old = det_old[det_old.n_total.isin(N_GRID)]
    for _, r_old in det_old.iterrows():
        m = (det_new.n_total == r_old.n_total) & np.isclose(det_new.f, r_old.f, atol=1e-12)
        if m.sum() != 1:
            add(f"detection n={r_old.n_total} f={r_old.f} row present", True, False, "exact")
            continue
        r_new = det_new[m].iloc[0]
        tag = f"detection n={int(r_old.n_total)} f={float(r_old.f):.2f}"
        for c in det_old.columns:
            if c in ("n_total", "f"):
                continue
            if c in EXACT_COLS:
                cast = str if c == "role" else int
                add(f"{tag} {c}", cast(r_old[c]), cast(r_new[c]), "exact")
            else:
                add(f"{tag} {c}", float(r_old[c]), float(r_new[c]), "tol")
    # (c) per-draw frac and delta
    dr_old = pd.read_csv(DEPTH_STORED / "draws.csv")
    dr_old = dr_old[dr_old.n_total.isin(N_GRID)]
    mg = dr_old.merge(deep_df, on=["n_total", "f_index", "draw"], suffixes=("_old", "_new"),
                      how="outer", indicator=True)
    add("draws.csv 400/1200 rows matched (extra check)", int(len(dr_old)),
        int((mg["_merge"] == "both").sum()), "exact")
    for c in ("n_aged", "n_young"):
        add(f"draws.csv per-draw {c} all equal (extra check)", True,
            bool((mg[f"{c}_old"] == mg[f"{c}_new"]).all()), "exact")
    for c in ("frac", "delta"):
        old = mg[f"{c}_old"].to_numpy(float)
        new = mg[f"{c}_new"].to_numpy(float)
        ok = np.array([close_enough(a, b) for a, b in zip(new, old)])
        rows.append(dict(item=f"draws.csv per-draw {c} within tol, n_failed (extra check)",
                         kind="exact", gated=True, stored=0, reproduced=int((~ok).sum()),
                         abs_diff=float(np.nanmax(np.abs(new - old))), ok=bool(ok.all())))
    # (d) monotonic, (e) readings
    mono_old = pd.read_csv(DEPTH_STORED / "monotonic.csv")
    read_old = load_json(DEPTH_STORED / "readings.json")
    for n in N_GRID:
        add(f"monotonic n={n} n_monotonic (extra check)",
            int(mono_old[mono_old.n_total == n].n_monotonic.iloc[0]),
            int(mono_new[mono_new.n_total == n].n_monotonic.iloc[0]), "exact")
        add(f"readings n={n} rule_outcome (extra check)",
            str(read_old[str(n)]["rule_outcome"]), str(read_new[str(n)]["rule_outcome"]), "exact")
    return pd.DataFrame(rows), det_new, mono_new, read_new


def split_checks(obs_a, obs_y, idx):
    rows = []

    def add(item, stored, new):
        rows.append(dict(item=item, kind="exact", gated=True, stored=stored, reproduced=new,
                         abs_diff=np.nan, ok=bool(stored == new)))

    o_new, ab_new = split_half(d0_fibroblast_idx(obs_a), SAME_SEED)
    y_new, yb_new = split_half(d0_fibroblast_idx(obs_y), SAME_SEED)
    for key, arr in (("O", o_new), ("aged_B", ab_new), ("Y", y_new), ("young_B", yb_new)):
        st = idx[key]
        add(f"split_idx.npz {key} indices identical (extra check)", True,
            bool(st.shape == arr.shape and np.array_equal(st, arr)))
    add("O pool (GM00731 half A) and aged_B disjoint (extra check)", 0,
        int(np.intersect1d(idx["O"], idx["aged_B"]).size))
    add("Y pool (GM23815 half A) and young_B disjoint (extra check)", 0,
        int(np.intersect1d(idx["Y"], idx["young_B"]).size))
    return rows


def stage0(Y_a, Y_y, obs_a, obs_y, frozen, log):
    STAGE0_DIR.mkdir(parents=True, exist_ok=True)
    idx = read_pools()
    srows = split_checks(obs_a, obs_y, idx)
    deep = draws_deep(Y_a, Y_y, idx, frozen, log)
    deep.to_csv(STAGE0_DIR / "draws_deep.csv", index=False)
    s0, det, mono, read = compare_stage0(deep, srows)
    det.to_csv(STAGE0_DIR / "detection_deep.csv", index=False)
    mono.to_csv(STAGE0_DIR / "monotonic_deep.csv", index=False)
    dump_json(STAGE0_DIR / "readings_deep.json", jsonable(read))
    s0.to_csv(DEPTH2_DIR / "stage0_reproduction.csv", index=False)
    bad = s0[s0.gated & ~s0.ok]
    log(f"[stage0] checks={len(s0)} gated={int(s0.gated.sum())} gated_failures={len(bad)}")
    if len(bad):
        names = "; ".join(f"{r.item}: stored={r.stored} reproduced={r.reproduced}"
                          for r in bad.itertuples())
        raise StopStep("prior_results_not_reproduced", names,
                       dict(mismatches=bad.item.tolist()))
    return idx


# --------------------------------------------------------------------------- draws
def check_pools(idx):
    for n in N_GRID:
        if n > idx["O"].size or n > idx["Y"].size:
            raise StopStep("pool_too_small",
                           f"n_total={n} needs {n} reference cells; half-A pools are "
                           f"{idx['O'].size}/{idx['Y'].size}.")
        for f in F_VALUES:
            n_aged, n_young = n_split(n, f)
            if n_aged > idx["aged_B"].size or n_young > idx["young_B"].size:
                raise StopStep("pool_too_small",
                               f"n_total={n} f={f} needs n_aged={n_aged} n_young={n_young}; "
                               f"half-B pools are {idx['aged_B'].size}/{idx['young_B'].size}.")


def role_rng(n_index, f_index, draw, role):
    return np.random.default_rng([DEPTH2_SEED, n_index, f_index, draw, role])


def draw_matched(n_index, f_index, draw, n_total, f, idx):
    o = role_rng(n_index, f_index, draw, ROLE_O).choice(idx["O"], size=n_total, replace=False)
    y = role_rng(n_index, f_index, draw, ROLE_Y).choice(idx["Y"], size=n_total, replace=False)
    rng = role_rng(n_index, f_index, draw, ROLE_S)
    n_aged, n_young = n_split(n_total, f)
    ta = rng.choice(idx["aged_B"], size=n_aged, replace=False) if n_aged > 0 else np.array([], int)
    ty = rng.choice(idx["young_B"], size=n_young, replace=False) if n_young > 0 else np.array([], int)
    out = [np.asarray(a, int) for a in (o, y, ta, ty)]
    for a, n in zip(out, (n_total, n_total, n_aged, n_young)):
        if a.size != n or np.unique(a).size != n:
            raise StopStep("draw", f"n={n_total} f={f} draw={draw}: bad sample sizes.")
    return n_aged, n_young, out


def run_draws(Y_a, Y_y, idx, frozen, log):
    check_pools(idx)
    rows = []
    for n_index, n_total in enumerate(N_GRID):
        for f_index, f in enumerate(F_VALUES):
            for draw in range(N_DRAWS):
                n_aged, n_young, (o, y, ta, ty) = draw_matched(
                    n_index, f_index, draw, n_total, f, idx)
                counts = [sum_rows(Y_a, o), sum_rows(Y_y, y), mix_sum(Y_a, Y_y, ta, ty)]
                Z, missing, _, _ = panel_z(counts, frozen)
                rows.append(score_row(Z, n_index, n_total, f_index, f, draw,
                                      n_aged, n_young, missing))
            log(f"[draws] n_total={n_total} f={f} {N_DRAWS}/{N_DRAWS}")
    df = pd.DataFrame(rows)
    df.to_csv(DEPTH2_DIR / "draws.csv", index=False)
    return df


# --------------------------------------------------------------------------- findings
def _num(x, d=3):
    return "NA" if x is None or not np.isfinite(float(x)) else f"{float(x):+.{d}f}"


def _pct(x):
    return f"{100.0 * float(x):.1f}%"


def write_findings(prereg, stop=None, s0=None, det=None, mono=None, readings=None):
    L = ["# FINDINGS_DEPTH2 — the DEPTH test with references at the experiment's own depth", ""]
    if stop:
        L += [f"**Status: STOP `{stop['step']}`.** {stop['message']}", ""]
    else:
        keys = ", ".join(f"n_total={n}: `{readings[str(n)]['key']}`" for n in N_GRID)
        L += [f"**Status: DONE.** Readings: {keys}.", ""]
    L += [
        "Reproduced by `python src/depth2_run.py`; Stage 0 comparison, tables and readings "
        "re-checked from tracked files by `python src/depth2_run.py verify`. DEPTH's readings "
        "(400 `INSUFFICIENT`, 1,200 `MARGINAL`, FINDINGS_DEPTH.md) stand on record and are "
        "not revised here.",
        "",
        "## Pre-registration (verbatim from `results/depth2/PREREG.flag`)",
        "",
        "```",
        prereg,
        "```",
        "",
    ]
    if s0 is not None:
        L += [
            "## Stage 0 — reproduction of DEPTH's stored numbers (DEPTH configuration, deep anchors)",
            "",
            "This task's code, run with DEPTH's deep anchors, 6-row panels and seeds, compared "
            "with `results/depth/` for n_total = 400 and 1,200. Counts exact; numerics within "
            "1e-4 x max(1, |stored|). All rows are gated; rows marked \"extra check\" were added "
            "beyond the task text. For the two per-draw rows, `stored` is the required number "
            "of failures (0), `reproduced` the number found, and `abs_diff` the largest absolute "
            "difference over the 1,600 draws. Full table: `results/depth2/stage0_reproduction.csv`.",
            "",
        ]
        show = s0.copy()
        show["stored"] = [str(v) for v in show["stored"]]
        show["reproduced"] = [str(v) for v in show["reproduced"]]
        show["abs_diff"] = ["" if not np.isfinite(v) else f"{v:.2e}" for v in show["abs_diff"]]
        L += [md_table(show[["item", "kind", "gated", "stored", "reproduced", "abs_diff", "ok"]]), ""]
        n_gf = int((s0.gated & ~s0.ok).sum())
        L += [f"{len(s0)} checks, {int(s0.gated.sum())} gated, {n_gf} gated failures.", ""]
    if stop:
        L += _not_covered() + _deviations() + _files(stop=True)
        FINDINGS_PATH.write_text("\n".join(L) + "\n", encoding="utf-8")
        return
    L += ["## Pre-registered readings (only the reading that fired, per depth)", ""]
    for n in N_GRID:
        r = readings[str(n)]
        c = r["conditions"]
        base = det[(det.n_total == n) & (det.f == 0.0)].detection_rate.iloc[0]
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
            f"Chance baseline (f=0 draws against their own floor, in-sample, at most 5% by "
            f"construction): {_pct(base)}.",
            "",
        ]
    L += [
        "## Detection table",
        "",
        "Per depth and f, over 200 draws; O, Y and S all at n_total cells, each (n_total, f, "
        "draw) its own 3-row panel. `floor_p5_delta` is the 5th percentile of that depth's "
        "f=0 deltas. A draw is detected iff delta < floor_p5_delta and delta < 0. The f=0 row "
        "is the chance baseline: the same rule applied to the floor's own draws. "
        "File: `results/depth2/detection.csv`.",
        "",
    ]
    show = det[["n_total", "f", "n_aged", "n_young", "n_detected", "detection_rate",
                "floor_p5_delta", "median_frac", "median_delta"]].copy()
    show["detection_rate"] = [_pct(v) for v in show["detection_rate"]]
    show["f"] = [f"{v:.2f}" for v in show["f"]]
    L += [md_table(show), ""]
    L += [
        "Monotonicity: share of the 200 draw indices in which frac(f=0) <= frac(0.10) <= "
        "frac(0.25) <= frac(0.50). Within a draw index each f has its own independently drawn "
        "O, Y and S. File: `results/depth2/monotonic.csv`.",
        "",
    ]
    ms = mono.copy()
    ms["monotonic_share"] = [_pct(v) for v in ms["monotonic_share"]]
    L += [md_table(ms), ""]
    L += _notes(det)
    L += _diagnostics(det, mono)
    L += _not_covered() + _deviations() + _files(stop=False)
    FINDINGS_PATH.write_text("\n".join(L) + "\n", encoding="utf-8")


def _notes(det):
    """Added after the first run's readings existed. Changes no reading; not gated."""
    need = int(np.ceil(RATE_BAR * N_DRAWS))
    L = [
        "## Notes on the readings (not gated; added after the first run's numbers existed)",
        "",
        "This section changes no reading. It records what the reading keys alone do not say.",
        "",
    ]
    for n in N_GRID:
        g = det[(det.n_total == n) & np.isclose(det.f, 0.10)].iloc[0]
        p = float(g.detection_rate)
        se = (f"binomial standard error of the rate at this value over {N_DRAWS} draws: "
              f"{100.0 * np.sqrt(p * (1 - p) / N_DRAWS):.1f} points" if 0.0 < p < 1.0 else
              f"rate at the boundary, so a binomial standard error is not informative; the "
              f"rule-of-three 95% bound is {_pct(3.0 / N_DRAWS if p == 0.0 else 1.0 - 3.0 / N_DRAWS)}")
        L.append(
            f"- Margin at n_total={n:,}, f=0.10: {int(g.n_detected)} of {N_DRAWS} detected, "
            f"{need} needed for the 80% bar ({se})."
        )
    L += [
        "- `USABLE` is the pre-registered key for meeting the 80% bars. Its wording says the "
        "depth \"preserves the sensitivity the mixing control was validated at\"; it does not "
        "mean detection is unchanged. At 400 cells the 10% rung is detected in "
        + _pct(det[(det.n_total == 400) & np.isclose(det.f, 0.10)].detection_rate.iloc[0])
        + " of draws, against 100% in DEPTH's 2,470 reference rung (deep anchors).",
        "- frac at f=0 is not zero: its median is "
        + ", ".join(f"{_num(det[(det.n_total == n) & (det.f == 0.0)].median_frac.iloc[0])} at {n:,}"
                    for n in N_GRID)
        + ". O's sampling noise enters both d = z_S - z_O and v = z_Y - z_O, which adds a "
        "positive term to dot(d, v) and grows as depth falls; this is a property of the "
        "statistic, not tested separately here. frac is therefore only readable against the "
        "same depth's f=0 draws, never as an absolute number across depths.",
        "- The panel-shape change (6 rows in DEPTH, 3 here) is confounded with the anchor-depth "
        "change. This task cannot say how much of the difference from DEPTH's 400 "
        "`INSUFFICIENT` is due to each.",
        "- The f=0 floor contains cell-sampling noise of one aged line only. A pilot's aged "
        "reference from different donors would add between-donor spread, widening the floor; "
        "nothing here measures that.",
        "",
    ]
    return L


def _diagnostics(det, mono):
    L = [
        "## Diagnostics (not gated)",
        "",
        "### Which clause of the detection rule binds",
        "",
        "The same breakdown DEPTH reported: share of draws below the f=0 floor, share with "
        "delta < 0, and a frac-based separation that the rule does not use (share of draws "
        "whose frac exceeds the f=0 draws' 95th percentile). Not gated.",
        "",
    ]
    cl = det[["n_total", "f", "detection_rate", "diag_share_below_floor",
              "diag_share_delta_lt_0", "diag_p95_frac_f0", "diag_share_frac_gt_f0_p95"]].copy()
    cl["f"] = [f"{v:.2f}" for v in cl["f"]]
    for c in ("detection_rate", "diag_share_below_floor", "diag_share_delta_lt_0",
              "diag_share_frac_gt_f0_p95"):
        cl[c] = [_pct(v) for v in cl[c]]
    L += [md_table(cl), ""]
    for n in N_GRID:
        g = det[det.n_total == n].set_index("f")
        L.append(
            f"- n_total={n:,}: the f=0 floor's 5th-percentile delta is "
            f"{_num(g.loc[0.0, 'floor_p5_delta'])}, and {_pct(g.loc[0.0, 'diag_share_delta_lt_0'])} "
            f"of f=0 draws have delta < 0. At f=0.10, {_pct(g.loc[0.10, 'diag_share_below_floor'])} "
            f"of draws fall below the floor and {_pct(g.loc[0.10, 'diag_share_delta_lt_0'])} have "
            f"delta < 0. At f=0.25 the shares are {_pct(g.loc[0.25, 'diag_share_below_floor'])} and "
            f"{_pct(g.loc[0.25, 'diag_share_delta_lt_0'])}. frac exceeds the f=0 draws' 95th "
            f"percentile ({_num(g.loc[0.0, 'diag_p95_frac_f0'])}) in "
            f"{_pct(g.loc[0.10, 'diag_share_frac_gt_f0_p95'])} of f=0.10 draws and "
            f"{_pct(g.loc[0.25, 'diag_share_frac_gt_f0_p95'])} of f=0.25 draws."
        )
    L += [
        "",
        "### Spread",
        "",
        "5th/95th percentiles of frac and delta over draws, median cos, median ||v|| "
        "(v = z(Y) - z(O) is recomputed per draw). Not gated.",
        "",
    ]
    dg = det[["n_total", "f", "diag_p5_frac", "diag_p95_frac", "diag_p5_delta",
              "diag_p95_delta", "diag_median_cos", "diag_median_v_norm"]].copy()
    dg["f"] = [f"{v:.2f}" for v in dg["f"]]
    L += [md_table(dg), ""]
    old = pd.read_csv(DEPTH_STORED / "detection.csv")
    old_m = pd.read_csv(DEPTH_STORED / "monotonic.csv")
    rows = []
    for n in N_GRID:
        for f in F_VALUES:
            o = old[(old.n_total == n) & np.isclose(old.f, f)].iloc[0]
            d = det[(det.n_total == n) & np.isclose(det.f, f)].iloc[0]
            rows.append(dict(n_total=n, f=f"{f:.2f}",
                             DEPTH_detection=_pct(o.detection_rate),
                             DEPTH2_detection=_pct(d.detection_rate),
                             DEPTH_floor_p5=float(o.floor_p5_delta),
                             DEPTH2_floor_p5=float(d.floor_p5_delta),
                             DEPTH_median_frac=float(o.median_frac),
                             DEPTH2_median_frac=float(d.median_frac)))
    L += [
        "### DEPTH beside DEPTH2 (descriptive)",
        "",
        "DEPTH's stored numbers (deep anchors, one 6-row panel per draw index) beside "
        "DEPTH2's (anchors at n_total, one 3-row panel per draw). Two things changed between "
        "the columns: anchor depth, and the TMM panel shape. A difference is therefore not "
        "attributable to anchor depth alone. Descriptive only; DEPTH's readings are not "
        "revised.",
        "",
        md_table(pd.DataFrame(rows)),
        "",
        "Monotonic share, DEPTH vs DEPTH2: " + "; ".join(
            f"{n:,}: {_pct(old_m[old_m.n_total == n].monotonic_share.iloc[0])} vs "
            f"{_pct(mono[mono.n_total == n].monotonic_share.iloc[0])}" for n in N_GRID) + ".",
        "",
    ]
    return L


def _not_covered():
    return [
        "## What this does not cover",
        "",
        "- Two donors, one young (GM23815) and one aged (GM00731), from a reprogramming "
        "dataset (GSE297234), not the six donors of the planned pilot.",
        "- Reads per cell are not modelled. Cells are subsampled from a deeply sequenced run; "
        "a CRO's recovered cells may be shallower per cell. The result is optimistic in that "
        "respect.",
        "- The spike-in is genuinely young cells, the most favourable possible signal. A "
        "partial approach to young by an intervention is weaker and less clean.",
        "- Nothing here about donor-to-donor variation. Each reference is a resample of one "
        "donor line, so the f=0 floor contains cell-sampling noise only, not the between-donor "
        "spread a real aged reference would carry.",
        "- Circularity, disclosed: the f=0 floor and the spike-in draws come from the same "
        "pools and share cells across draws, and the floor is the 5th percentile of the same "
        "draws the chance baseline is computed on.",
        "- Day-0 fibroblasts of two cell lines only. Nothing here speaks to other cell types, "
        "to partially reprogrammed cells, or to donors without an aged counterpart.",
        "",
    ]


def _deviations():
    return [
        "## Deviations and decisions beyond the task text",
        "",
        "- **Scoring panel** (not fixed by the task text; fixed in the flag before computing): "
        "the seed key includes f_index for O and Y, so the references are independent per "
        "(n_total, f, draw) and DEPTH's shared 6-row panel cannot be formed. Each draw is its "
        "own 3-row TMM panel [O, Y, S]. This changes TMM's normalization context relative to "
        "DEPTH, so DEPTH-vs-DEPTH2 differences are not attributable to anchor depth alone.",
        "- **Monotonicity under independent references:** within a draw index, the four f "
        "values have different O and Y. Kept as the definition, stated in the flag.",
        "- **n_index mapping** 0 = 400, 1 = 1200 (DEPTH's), stated in the flag; the task "
        "text gave the key but not the mapping.",
        "- **Detection-table code** is this task's own function (`tables_from_draws(df, "
        "depths)`), because DEPTH's version iterates over its module-level 3-depth grid. "
        "It applies DEPTH's rule and produces DEPTH's columns; Stage 0 checks it against "
        "every column of DEPTH's stored table. One diagnostic column was added: "
        "`diag_median_v_norm`. Readings text, the rescue sentence, the floor percentile, "
        "the 0.80 bar, the tolerance and `n_split` are imported from `src/depth_run.py`, "
        "and the Stage 0 S draws use `depth_run.draw_cells`.",
        "- **Extra Stage 0 checks** (all gated): split indices recomputed and equal to "
        "`results/same/split_idx.npz`; disjointness of the O pool vs aged_B and the Y pool vs "
        "young_B; all 1,600 per-draw frac and delta values of `results/depth/draws.csv` "
        "(400, 1,200); `monotonic.csv` counts; `readings.json` rule outcomes.",
        "- **Stage 0 outputs** are written to `results/depth2/stage0/` (draws_deep.csv, "
        "detection_deep.csv, monotonic_deep.csv, readings_deep.json). Nothing was written "
        "into `results/depth/` or `results/same/`.",
        "- **Import side effects:** importing `same_run` and `toward_run` (directly and via "
        "`depth_run`) runs `mkdir(exist_ok=True)` on `results/same/` and `results/toward/`. "
        "Both already existed; no file in them was written. Importing `depth_run` writes "
        "nothing.",
        "- **Added diagnostics** beyond the task text: spread percentiles of frac and delta, "
        "median cos, median ||v||, and DEPTH's stored rates printed beside DEPTH2's. All "
        "declared in the flag, none gated.",
        "- **Added after the readings existed:** the section \"Notes on the readings\". It was "
        "added after the first complete run had produced the readings. Then the whole script "
        "was re-run from scratch, so the committed code produces every committed file. The "
        "readings and every gated number of the re-run are identical to the first run.",
        "- **Not opened:** `results/transfer/NOT_A_RESULT_gate_failed/`.",
        "",
    ]


def _files(stop):
    L = [
        "## Files",
        "",
        "- `src/depth2_run.py`",
        "- `results/depth2/PREREG.flag`",
        "- `results/depth2/stage0_reproduction.csv`",
        "- `results/depth2/stage0/` (DEPTH-configuration draws_deep.csv, detection_deep.csv, "
        "monotonic_deep.csv, readings_deep.json)",
        "- `results/depth2/run.log`",
    ]
    if stop:
        L += ["- `results/depth2/STOP.json`"]
    else:
        L += [
            "- `results/depth2/draws.csv` (one row per draw: frac, delta, cos, d_norm, v_norm, ...)",
            "- `results/depth2/detection.csv`",
            "- `results/depth2/monotonic.csv`",
            "- `results/depth2/readings.json`",
        ]
    L += ["- `FINDINGS_DEPTH2.md`", ""]
    return L


# --------------------------------------------------------------------------- main
def run(log):
    prereg = load_prereg_text()
    log("=" * 90)
    log(f"DEPTH2 seed={DEPTH2_SEED} N_GRID={N_GRID} f={F_VALUES} n_draws={N_DRAWS}")
    log("=" * 90)
    frozen = load_frozen_ruler()
    log(f"[frozen] {FROZEN_RULER} n_genes={len(np.asarray(frozen['mu']))}")
    obs_a, Y_a = load_donor(AGED_LINE, log)
    obs_y, Y_y = load_donor(YOUNG_LINE, log)
    same_run.check_gene_space(Y_a, frozen, AGED_LINE)
    same_run.check_gene_space(Y_y, frozen, YOUNG_LINE)
    try:
        idx = stage0(Y_a, Y_y, obs_a, obs_y, frozen, log)
        s0 = pd.read_csv(DEPTH2_DIR / "stage0_reproduction.csv")
        df = run_draws(Y_a, Y_y, idx, frozen, log)
        det, mono, readings = tables_from_draws(df, N_GRID)
    except StopStep as e:
        stop = dict(step=e.step, message=e.message, details=jsonable(e.details))
        dump_json(DEPTH2_DIR / "STOP.json", stop)
        p = DEPTH2_DIR / "stage0_reproduction.csv"
        write_findings(prereg, stop=stop, s0=pd.read_csv(p) if p.exists() else None)
        log(f"[STOP] {e.step}: {e.message}")
        return 1
    det.to_csv(DEPTH2_DIR / "detection.csv", index=False)
    mono.to_csv(DEPTH2_DIR / "monotonic.csv", index=False)
    dump_json(DEPTH2_DIR / "readings.json", jsonable(readings))
    for n in N_GRID:
        r = readings[str(n)]
        log(f"[reading] n_total={n} key={r['key']} {r['conditions']}")
    write_findings(prereg, s0=s0, det=det, mono=mono, readings=readings)
    log("[run] DONE")
    return 0


def verify():
    ok = True
    idx = read_pools()
    srows = [dict(item="O/aged_B disjoint", kind="exact", gated=True, stored=0,
                  reproduced=int(np.intersect1d(idx["O"], idx["aged_B"]).size), abs_diff=np.nan,
                  ok=np.intersect1d(idx["O"], idx["aged_B"]).size == 0),
             dict(item="Y/young_B disjoint", kind="exact", gated=True, stored=0,
                  reproduced=int(np.intersect1d(idx["Y"], idx["young_B"]).size), abs_diff=np.nan,
                  ok=np.intersect1d(idx["Y"], idx["young_B"]).size == 0)]
    s0, _, _, _ = compare_stage0(pd.read_csv(STAGE0_DIR / "draws_deep.csv"), srows)
    n_bad = int((s0.gated & ~s0.ok).sum())
    print(f"stage0 (from tracked draws_deep.csv; split recomputation needs raw inputs and is "
          f"not re-run here): checks={len(s0)} gated_failures={n_bad}")
    ok &= n_bad == 0
    df = pd.read_csv(DEPTH2_DIR / "draws.csv")
    det, mono, readings = tables_from_draws(df, N_GRID)
    det_s = pd.read_csv(DEPTH2_DIR / "detection.csv")
    mono_s = pd.read_csv(DEPTH2_DIR / "monotonic.csv")
    read_s = load_json(DEPTH2_DIR / "readings.json")
    for c in ("n_detected", "n_aged", "n_young"):
        ok &= bool((det[c].to_numpy() == det_s[c].to_numpy()).all())
    for c in ("detection_rate", "floor_p5_delta", "median_frac", "median_delta"):
        ok &= bool(np.allclose(det[c].to_numpy(float), det_s[c].to_numpy(float), rtol=1e-9, atol=1e-12))
    ok &= bool((mono["n_monotonic"].to_numpy() == mono_s["n_monotonic"].to_numpy()).all())
    for n in N_GRID:
        ok &= readings[str(n)]["key"] == read_s[str(n)]["key"]
        print(f"n_total={n}: key={readings[str(n)]['key']}")
    print("verify:", "OK" if ok else "MISMATCH")
    return 0 if ok else 1


def main():
    if len(sys.argv) > 1 and sys.argv[1] == "verify":
        return verify()
    DEPTH2_DIR.mkdir(parents=True, exist_ok=True)
    log = Logger(DEPTH2_DIR / "run.log")
    try:
        return run(log)
    finally:
        log.close()


if __name__ == "__main__":
    sys.exit(main())
