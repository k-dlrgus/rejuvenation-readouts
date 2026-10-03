"""COORDSPACE — coordinate-space robustness check (SECONDARY).

The pre-registration is results/coordspace/PREREG.flag, written before any
statistic. This script reads that flag, refuses to compute if it is missing,
and never rewrites it.

Three spaces. Everything upstream of the last transform step is shared: the
pseudobulk panels, the summed counts, the TMM factors (estimated on each
analysis's own panel exactly as SAME and SOUTH6 estimate them), the SAME
half-split, the mixture draws, and the bootstrap and permutation streams.

  A  TMM CPM, no log, no z          x = 1e6 * count / (nf * library size)
  B  TMM log2 CPM, no z             x = log2_cpm_edger(count, nf, prior=2)
  C  log2 CPM + frozen GTEx z       x = (B - mu) / sd   (the current pipeline)

Usage:
  python src/coordspace_run.py                       full run, Stages 0-4 and the decision table
  python src/coordspace_run.py --through-stage N     stop after Stage N (0-3)

Every run starts at Stage 0. There is no resume path. A run that stops early is
repeated from the start, and every random draw comes from its own fresh
generator seeded exactly as SAME and SOUTH6 seed it, so a repeated stage draws
the same numbers.

src/same_run.py, src/south6_run.py, src/toward_run.py, src/seng_run.py and
src/gtex_common.py are imported and not modified. The frozen ruler is not refit.
"""
from __future__ import annotations

import argparse
import multiprocessing as mp
import os
import shutil
import sys
import tempfile
import time
from pathlib import Path

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parent))
from config import RESULTS, ROOT  # noqa: E402
from brain_phase1_common import Logger, dump_json  # noqa: E402
from fibro_common import jsonable  # noqa: E402
from fibro2_common import load_frozen_ruler  # noqa: E402
from fibro3_common import AGED_LINE, YOUNG_LINE  # noqa: E402
from target_common import md_table, unit  # noqa: E402
import same_run as SR  # noqa: E402
import toward_run as TR  # noqa: E402
import south6_run as S6R  # noqa: E402
from seng_run import (  # noqa: E402
    bh_q,
    logcpm_one,
    tmm_cache_from_panel,
    unscaled_tmm_factor,
    z_from_counts,
)


COORD_DIR = RESULTS / "coordspace"
COORD_DIR.mkdir(parents=True, exist_ok=True)
PREREG = COORD_DIR / "PREREG.flag"
MANIFEST = COORD_DIR / "manifest.json"
LOG_PATH = COORD_DIR / "coordspace_run.log"
FINDINGS = ROOT / "FINDINGS_COORDSPACE.md"
PROGRESS = ROOT / "PROGRESS_COORDSPACE.md"
SOUTH6_DIR = RESULTS / "south6"

SPACES = ("A", "B", "C")
TOL_SAME = 1e-6
TOL_SOUTH6 = 1e-4
S6_SEED = S6R.SEED
S6_N_PERM = S6R.N_PERM
S6_PERM_CHUNK = S6R.PERM_CHUNK
S6_N_SPLITS = S6R.N_SPLITS
S6_Q_BAR = S6R.Q_BAR
S6_PRIMARY_SETTING_INDEX = 0  # S1 O1_Y1 is the first active setting in south6_run.make_settings

PREREG_TEXT_NOTES = [
    dict(
        item="matched coordinates",
        flag_says="20,265 matched genes; permutation across the 20,265 matched coordinates",
        actual="20,264. SOUTH6's coverage table counts 20,265 ruler genes matched by "
               "best_after_md3_idtype; two of them map to one source column, which SOUTH6 "
               "collapses, so its panel, its d_g and its Stage 2A permutation run over 20,264 "
               "coordinates (results/south6/stage1.npz, matched_ruler).",
        effect="None on any computation. The pre-registered rule is 'the same matched genes, "
               "unchanged', and this task uses SOUTH6's own matched set.",
    ),
    dict(
        item="identity_loss units",
        flag_says="the 0.5 identity_loss threshold is defined in frozen-z units",
        actual="In src/south6_run.py identity_drop is computed from the panel's log2-CPM "
               "(logcpm_f, logcpm_pool), so the 0.5 threshold is in log2-CPM units, the units of "
               "space B, not z units.",
        effect="None on any computation. The flags are carried unchanged from SOUTH6 space C as "
               "pre-registered. In space B a recomputed identity_drop would equal the carried one "
               "exactly; space A still has no counterpart.",
    ),
]


class Halt(Exception):
    def __init__(self, step, message, details=None):
        super().__init__(f"[{step}] {message}")
        self.step = step
        self.message = message
        self.details = details or {}


# ----------------------------------------------------------------------------- bookkeeping
def load_prereg():
    if not PREREG.exists():
        raise SystemExit(f"PREREG.flag missing at {PREREG}. Not computing.")
    return PREREG.read_text(encoding="utf-8").strip()


def new_manifest():
    return dict(
        task="COORDSPACE", status="RUNNING", failures=[], notes=[],
        prereg_text_notes=PREREG_TEXT_NOTES, stages={},
        seeds=dict(same=SR.SAME_SEED, same_boot=SR.BOOT_SEED, south6=S6_SEED),
    )


def save_manifest(man):
    dump_json(MANIFEST, jsonable(man))


def scaled_diff(new, stored):
    new = np.asarray(new, np.float64).ravel()
    stored = np.asarray(stored, np.float64).ravel()
    if new.shape != stored.shape:
        return np.inf, np.inf, False
    both_nan = np.isnan(new) & np.isnan(stored)
    d = np.abs(new - stored)
    d[both_nan] = 0.0
    sc = d / np.maximum(1.0, np.abs(stored))
    sc[both_nan] = 0.0
    bad = ~np.isfinite(sc)
    if bad.any():
        return np.inf, np.inf, False
    return float(np.max(d)) if d.size else 0.0, float(np.max(sc)) if sc.size else 0.0, True


def cmp_num(rows, block, item, new, stored, tol, gating=True):
    max_abs, max_sc, ok_shape = scaled_diff(new, stored)
    ok = bool(ok_shape and max_sc <= tol)
    rows.append(dict(
        block=block, item=item, kind="numeric", n=int(np.asarray(stored).size),
        max_abs_diff=max_abs, max_scaled_diff=max_sc, tol=tol,
        stored_example=_example(stored), new_example=_example(new),
        n_mismatch=_n_mismatch(new, stored, tol), pass_=ok, gating=gating,
    ))
    return ok


def cmp_cat(rows, block, item, new, stored, gating=True):
    a = np.asarray(new).ravel()
    b = np.asarray(stored).ravel()
    ok = bool(a.shape == b.shape and np.array_equal(a, b))
    nm = int(np.sum(a != b)) if a.shape == b.shape else int(max(a.size, b.size))
    rows.append(dict(
        block=block, item=item, kind="categorical", n=int(b.size),
        max_abs_diff=np.nan, max_scaled_diff=np.nan, tol=np.nan,
        stored_example=_example(b), new_example=_example(a),
        n_mismatch=nm, pass_=ok, gating=gating,
    ))
    return ok


def _example(x):
    a = np.asarray(x).ravel()
    if a.size == 0:
        return ""
    v = a[0]
    if isinstance(v, (float, np.floating)):
        return f"{float(v):.10g}"
    return str(v)


def _n_mismatch(new, stored, tol):
    a = np.asarray(new, np.float64).ravel()
    b = np.asarray(stored, np.float64).ravel()
    if a.shape != b.shape:
        return int(max(a.size, b.size))
    both_nan = np.isnan(a) & np.isnan(b)
    sc = np.abs(a - b) / np.maximum(1.0, np.abs(b))
    sc[both_nan] = 0.0
    return int(np.sum(~(sc <= tol)))


# ----------------------------------------------------------------------------- transforms
def transform_panel(C, mu, sd, space):
    """C is rows x genes summed counts. TMM on this panel exactly as SAME/SOUTH6 do it.
    Genes with no counts in any row are set to one constant per gene: 0 in A, mu in B, 0 in C."""
    C = np.asarray(C, np.float64)
    mu = np.asarray(mu, np.float64)
    sd = np.asarray(sd, np.float64)
    logcpm, nf = TR.tmm_logcpm_quiet(C)
    missing = (C == 0).all(0)
    if space == "C":
        sd_safe = np.where(sd < 1e-12, 1.0, sd)
        X = (logcpm - mu) / sd_safe
        X[:, missing] = 0.0
    elif space == "B":
        X = np.array(logcpm, np.float64, copy=True)
        X[:, missing] = mu[missing]
    elif space == "A":
        lib = C.sum(1) * np.asarray(nf, np.float64)
        X = C / lib[:, None] * 1e6
        X[:, missing] = 0.0
    else:
        raise ValueError(space)
    return X, missing, nf, logcpm


def same_panel(count_rows, frozen, space):
    if space == "C":
        Z, missing, nf, C = SR.panel_z(count_rows, frozen)
        return Z, missing, nf
    C = np.vstack([np.asarray(r, np.float64).ravel() for r in count_rows])
    X, missing, nf, _ = transform_panel(C, frozen["mu"], frozen["sd"], space)
    return X, missing, nf


def x_one(counts, cache, mu, sd, missing, space):
    """One new row normalised against a cached panel reference, as seng_run.z_from_counts does."""
    counts = np.asarray(counts, np.float64).ravel()
    if space == "C":
        z, _, _, _ = z_from_counts(counts, cache, mu, sd, missing)
        return z
    fac = unscaled_tmm_factor(counts, float(np.sum(counts)), cache["ref_counts"], cache["ref_lib"])
    nf = fac / cache["geom"]
    if space == "B":
        lc, _ = logcpm_one(counts, nf, cache["ave"])
        x = np.array(lc, np.float64, copy=True)
        x[missing] = np.asarray(mu, np.float64)[missing]
        return x
    lib = float(np.sum(counts)) * float(nf)
    x = counts / lib * 1e6
    x[missing] = 0.0
    return x


def share10(x, symbols=None):
    """Share of the squared Euclidean length carried by the 10 largest squared coordinates."""
    x = np.asarray(x, np.float64).ravel()
    sq = x * x
    tot = float(sq.sum())
    if not np.isfinite(tot) or tot <= 0:
        return np.nan, []
    top = np.argpartition(sq, -10)[-10:]
    top = top[np.argsort(-sq[top])]
    genes = [str(symbols[j]) for j in top] if symbols is not None else [int(j) for j in top]
    return float(sq[top].sum() / tot), genes


def share10_rows(M, offset=None, chunk=200):
    out = np.full(M.shape[0], np.nan)
    for a in range(0, M.shape[0], chunk):
        v = np.asarray(M[a:a + chunk], np.float64)
        if offset is not None:
            v = v + np.asarray(offset, np.float64)[None, :]
        sq = v * v
        tot = sq.sum(1)
        top = -np.partition(-sq, 9, axis=1)[:, :10]
        with np.errstate(divide="ignore", invalid="ignore"):
            out[a:a + chunk] = np.where(tot > 0, top.sum(1) / tot, np.nan)
    return out


def gtex_centroids(space, log):
    """O1 and Y1 in one space: the arithmetic mean of the member donors' coordinates.
    C is read from results/toward/anchors.npz as stored; B and A are built from the GTEx pack."""
    anc = np.load(TR.ANCHORS_PATH, allow_pickle=True)
    if space == "C":
        return np.asarray(anc["c_old"], np.float64), np.asarray(anc["c_young"], np.float64)
    gtex = TR.load_gtex_z(log)
    age = gtex["age"]
    young_m = np.isin(age, TR.YOUNG_BINS)
    old_m = np.isin(age, TR.OLD_BINS)
    if int(young_m.sum()) != int(anc["n_young"]) or int(old_m.sum()) != int(anc["n_old"]):
        raise Halt("gtex_anchor_mismatch", "GTEx YOUNG/OLD donor counts differ from anchors.npz.")
    if space == "B":
        X = np.asarray(gtex["pack"]["X"], np.float64)
    else:
        counts = np.asarray(gtex["pack"]["counts"], np.float64).T
        nf = np.asarray(gtex["pack"]["tmm_factors"], np.float64).ravel()
        lib = counts.sum(1) * nf
        X = counts / lib[:, None] * 1e6
    return X[old_m].mean(0), X[young_m].mean(0)


# ----------------------------------------------------------------------------- SAME
def same_load(log):
    frozen = load_frozen_ruler()
    obs_a, Y_a = TR.load_donor(AGED_LINE, log)
    obs_y, Y_y = TR.load_donor(YOUNG_LINE, log)
    SR.check_gene_space(Y_a, frozen, AGED_LINE)
    SR.check_gene_space(Y_y, frozen, YOUNG_LINE)
    sp = np.load(SR.SAME_DIR / "split_idx.npz")
    o_idx = np.asarray(sp["O"], int)
    y_idx = np.asarray(sp["Y"], int)
    aged_b = np.asarray(sp["aged_B"], int)
    young_b = np.asarray(sp["young_B"], int)
    oa, ab = SR.split_half(SR.d0_fibroblast_idx(obs_a), SR.SAME_SEED)
    ya, yb = SR.split_half(SR.d0_fibroblast_idx(obs_y), SR.SAME_SEED)
    split_matches = bool(
        np.array_equal(oa, o_idx) and np.array_equal(ab, aged_b)
        and np.array_equal(ya, y_idx) and np.array_equal(yb, young_b)
    )
    log(f"[same] stored split O={o_idx.size} aged_B={aged_b.size} Y={y_idx.size} young_B={young_b.size}; "
        f"equals a fresh split_half draw: {split_matches}")
    return dict(
        frozen=frozen, obs_a=obs_a, Y_a=Y_a, obs_y=obs_y, Y_y=Y_y,
        o_idx=o_idx, y_idx=y_idx, aged_b=aged_b, young_b=young_b,
        split_matches=split_matches,
    )


def same_posctrl(space, S, log):
    """SAME positive control P (src/same_run.py run_posctrl), with the transform swapped."""
    frozen, Y_a, Y_y = S["frozen"], S["Y_a"], S["Y_y"]
    mixes, n_tot = SR.mix_indices(S["aged_b"], S["young_b"], SR.SAME_SEED)
    counts = [SR.sum_rows(Y_a, S["o_idx"]), SR.sum_rows(Y_y, S["y_idx"])]
    for m in mixes:
        counts.append(SR.mix_sum(Y_a, Y_y, m["aged_idx"], m["young_idx"]))
    X, missing, nf = same_panel(counts, frozen, space)
    obs_stats = []
    for i, m in enumerate(mixes):
        obs_stats.append(SR.axis_stats(X[0], X[1], X[2 + i]))
    ruler = TR.ruler_score(X, frozen) if space == "C" else np.full(X.shape[0], np.nan)

    rng_b = np.random.default_rng(SR.BOOT_SEED)
    boot = {m["f"]: dict(frac=[], delta=[], cos=[]) for m in mixes}
    c_O_obs, c_Y_obs = counts[0], counts[1]
    for b in range(SR.N_BOOT):
        b_counts = [c_O_obs, c_Y_obs]
        for m in mixes:
            b_counts.append(SR.mix_sum(
                Y_a, Y_y,
                SR.resample_idx(m["aged_idx"], rng_b),
                SR.resample_idx(m["young_idx"], rng_b),
            ))
        Xb, _, _ = same_panel(b_counts, frozen, space)
        for i, m in enumerate(mixes):
            st = SR.axis_stats(Xb[0], Xb[1], Xb[2 + i])
            boot[m["f"]]["frac"].append(st["frac"])
            boot[m["f"]]["delta"].append(st["delta"])
            boot[m["f"]]["cos"].append(st["cos"])
    rows = []
    for i, m in enumerate(mixes):
        st = obs_stats[i]
        flo, fhi = TR.percentile_ci(boot[m["f"]]["frac"])
        dlo, dhi = TR.percentile_ci(boot[m["f"]]["delta"])
        rows.append(dict(
            space=space, f=m["f"], n_aged=m["n_aged"], n_young=m["n_young"], n_total=m["n_total"],
            frac=st["frac"], delta=st["delta"], cos=st["cos"],
            d_norm=st["d_norm"], v_norm=st["v_norm"],
            frac_ci_lo=flo, frac_ci_hi=fhi, delta_ci_lo=dlo, delta_ci_hi=dhi,
            ruler_score_context=float(ruler[2 + i]),
            noise_floor=bool(m["f"] == 0.0),
            frac_ci_covers_point=SR.ci_covers_point(st["frac"], flo, fhi),
            delta_ci_covers_point=SR.ci_covers_point(st["delta"], dlo, dhi),
        ))
    pdf = pd.DataFrame(rows)
    fracs = [float(r["frac"]) for r in rows]
    row50 = next(r for r in rows if abs(r["f"] - 0.50) < 1e-12)
    mono = SR.monotonic_increasing(fracs)
    frac50_ex0 = SR.ci_excludes_zero(row50["frac_ci_lo"], row50["frac_ci_hi"])
    frac50_valid = bool(row50["frac_ci_covers_point"])
    delta50_neg = bool(np.isfinite(row50["delta"]) and row50["delta"] < 0)
    pass_same_rule = bool(mono and frac50_ex0 and delta50_neg)
    pass_coord = bool(mono and frac50_ex0 and frac50_valid and delta50_neg)
    summary = dict(
        space=space, pass_same_rule=pass_same_rule, pass_coordspace=pass_coord,
        monotonic=mono, frac50_ci_excludes_0=frac50_ex0, frac50_ci_valid=frac50_valid,
        delta50_neg=delta50_neg, n_total=n_tot, n_missing=int(np.asarray(missing).sum()),
    )
    log(f"[P {space}] pass(SAME rule)={pass_same_rule} pass(COORDSPACE rule)={pass_coord} "
        f"monotonic={mono} frac50_ci_excludes_0={frac50_ex0} frac50_ci_valid={frac50_valid} "
        f"delta50<0={delta50_neg}")
    return pdf, summary


def same_main(space, S, u_vec, log):
    """SAME forward, reverse, asymmetry and context (src/same_run.py run), transform swapped."""
    frozen = S["frozen"]
    obs_a, obs_y, Y_a, Y_y = S["obs_a"], S["obs_y"], S["Y_a"], S["Y_y"]
    s_idx = SR.state_idx(obs_a, SR.DEST_STATE)
    srev_idx = SR.state_idx(obs_y, SR.DEST_STATE)
    pluri_idx = SR.state_idx(obs_a, SR.PLURI_STATE)
    nr_idx = SR.state_idx(obs_a, SR.NONREPROG_STATE)
    panel_spec = [
        ("O", Y_a, S["o_idx"]), ("Y", Y_y, S["y_idx"]),
        ("S", Y_a, s_idx), ("S_rev", Y_y, srev_idx),
        ("Pluri", Y_a, pluri_idx), ("NonReprog", Y_a, nr_idx),
    ]
    names = [p[0] for p in panel_spec]
    pos = {nm: i for i, nm in enumerate(names)}
    counts = [SR.sum_rows(Y, ix) for _, Y, ix in panel_spec]
    X, missing, nf = same_panel(counts, frozen, space)
    x_O, x_Y = X[pos["O"]], X[pos["Y"]]
    ruler = TR.ruler_score(X, frozen) if space == "C" else np.full(X.shape[0], np.nan)
    fwd = SR.axis_stats(x_O, x_Y, X[pos["S"]])
    rev = SR.axis_stats(x_Y, x_O, X[pos["S_rev"]])
    pluri = SR.axis_stats(x_O, x_Y, X[pos["Pluri"]])
    nr = SR.axis_stats(x_O, x_Y, X[pos["NonReprog"]])
    asym = float(fwd["frac"] - rev["frac"])
    p1, _, _ = SR.n1_cos_p(fwd["d"], fwd["v"], SR.SAME_SEED, SR.N_PERM)
    p1r, _, _ = SR.n1_cos_p(rev["d"], rev["v"], SR.SAME_SEED, SR.N_PERM)
    cos_v_u = SR.cosine_full(fwd["v"], u_vec)

    rng_b = np.random.default_rng(SR.BOOT_SEED)
    bs = dict(frac_S=[], frac_rev=[], delta_S=[], delta_rev=[], asymmetry=[],
              frac_pluri=[], delta_pluri=[], frac_nr=[], delta_nr=[], cos_S=[], cos_rev=[])
    for b in range(SR.N_BOOT):
        b_counts = []
        for i, (_, Y, ix) in enumerate(panel_spec):
            if names[i] in ("O", "Y"):
                b_counts.append(counts[i])
            else:
                b_counts.append(SR.sum_rows(Y, SR.resample_idx(ix, rng_b)))
        Xb, _, _ = same_panel(b_counts, frozen, space)
        xO, xY = Xb[pos["O"]], Xb[pos["Y"]]
        st_f = SR.axis_stats(xO, xY, Xb[pos["S"]])
        st_r = SR.axis_stats(xY, xO, Xb[pos["S_rev"]])
        st_p = SR.axis_stats(xO, xY, Xb[pos["Pluri"]])
        st_n = SR.axis_stats(xO, xY, Xb[pos["NonReprog"]])
        bs["frac_S"].append(st_f["frac"])
        bs["frac_rev"].append(st_r["frac"])
        bs["delta_S"].append(st_f["delta"])
        bs["delta_rev"].append(st_r["delta"])
        bs["asymmetry"].append(st_f["frac"] - st_r["frac"])
        bs["frac_pluri"].append(st_p["frac"])
        bs["delta_pluri"].append(st_p["delta"])
        bs["frac_nr"].append(st_n["frac"])
        bs["delta_nr"].append(st_n["delta"])
        bs["cos_S"].append(st_f["cos"])
        bs["cos_rev"].append(st_r["cos"])
    ci = {k: TR.percentile_ci(v) for k, v in bs.items()}
    fwd["frac_ci_lo"], fwd["frac_ci_hi"] = ci["frac_S"]
    fwd["delta_ci_lo"], fwd["delta_ci_hi"] = ci["delta_S"]
    rev["frac_ci_lo"], rev["frac_ci_hi"] = ci["frac_rev"]
    rev["delta_ci_lo"], rev["delta_ci_hi"] = ci["delta_rev"]
    pluri["frac_ci_lo"], pluri["frac_ci_hi"] = ci["frac_pluri"]
    pluri["delta_ci_lo"], pluri["delta_ci_hi"] = ci["delta_pluri"]
    nr["frac_ci_lo"], nr["frac_ci_hi"] = ci["frac_nr"]
    nr["delta_ci_lo"], nr["delta_ci_hi"] = ci["delta_nr"]
    a_lo, a_hi = ci["asymmetry"]

    stats_df = pd.DataFrame([
        dict(space=space, test="forward", n_cells_S=int(s_idx.size), cos=fwd["cos"], frac=fwd["frac"],
             delta=fwd["delta"], p_N1=p1, frac_ci_lo=fwd["frac_ci_lo"], frac_ci_hi=fwd["frac_ci_hi"],
             delta_ci_lo=fwd["delta_ci_lo"], delta_ci_hi=fwd["delta_ci_hi"],
             d_norm=fwd["d_norm"], v_norm=fwd["v_norm"], ruler_score_S=float(ruler[pos["S"]]),
             dist_start=fwd["dist_origin_to_target"], dist_final=fwd["dist_to_target"]),
        dict(space=space, test="reverse", n_cells_S=int(srev_idx.size), cos=rev["cos"], frac=rev["frac"],
             delta=rev["delta"], p_N1=p1r, frac_ci_lo=rev["frac_ci_lo"], frac_ci_hi=rev["frac_ci_hi"],
             delta_ci_lo=rev["delta_ci_lo"], delta_ci_hi=rev["delta_ci_hi"],
             d_norm=rev["d_norm"], v_norm=rev["v_norm"], ruler_score_S=float(ruler[pos["S_rev"]]),
             dist_start=rev["dist_origin_to_target"], dist_final=rev["dist_to_target"]),
    ])
    ctx_rows = []
    for nm, rec, n in (("S", fwd, s_idx.size), ("S_rev", rev, srev_idx.size),
                       ("Pluri", pluri, pluri_idx.size), ("NonReprog", nr, nr_idx.size)):
        ctx_rows.append(dict(
            space=space, pseudobulk=nm, n_cells=int(n), ruler_score_context=float(ruler[pos[nm]]),
            frac=rec["frac"], delta=rec["delta"], cos=rec["cos"],
            frac_ci_lo=rec["frac_ci_lo"], frac_ci_hi=rec["frac_ci_hi"],
            delta_ci_lo=rec["delta_ci_lo"], delta_ci_hi=rec["delta_ci_hi"],
        ))
    ctx_df = pd.DataFrame(ctx_rows)
    reading_same_rule = SR.fire_reading(True, fwd, rev, asym, (a_lo, a_hi), p1, int(s_idx.size), pluri)
    out = dict(
        space=space, X=X, names=names, missing=missing, fwd=fwd, rev=rev, pluri=pluri, nr=nr,
        asym=asym, asym_ci=(a_lo, a_hi), p1=p1, p1r=p1r, cos_v_u=cos_v_u,
        stats_df=stats_df, ctx_df=ctx_df, reading_same_rule=reading_same_rule,
        ruler=dict(zip(names, [float(x) for x in ruler])),
        n=dict(S=int(s_idx.size), S_rev=int(srev_idx.size), Pluri=int(pluri_idx.size),
               NonReprog=int(nr_idx.size)),
    )
    log(f"[SAME {space}] delta_S={fwd['delta']:.6f} frac_S={fwd['frac']:.6f} "
        f"asym={asym:.6f} [{a_lo:.6f}, {a_hi:.6f}] key(SAME rule)={reading_same_rule['key']}")
    return out


def coord_reading(main, qualified):
    """SAME decision order 1-5 and overshoot flag, with the COORDSPACE interval check:
    a gating interval (frac_S, Asymmetry, Pluripotency frac) that does not contain its own
    point estimate makes its condition `ci_invalid`, and any key whose firing depends on it
    becomes `mixed`, naming the interval."""
    fwd, rev, pluri = main["fwd"], main["rev"], main["pluri"]
    a_lo, a_hi = main["asym_ci"]
    asym = main["asym"]
    dlt = float(fwd["delta"])
    invalid = []
    frac_valid = SR.ci_covers_point(fwd["frac"], fwd["frac_ci_lo"], fwd["frac_ci_hi"])
    asym_valid = SR.ci_covers_point(asym, a_lo, a_hi)
    pluri_valid = SR.ci_covers_point(pluri["frac"], pluri["frac_ci_lo"], pluri["frac_ci_hi"])
    if not frac_valid:
        invalid.append("frac_S CI")
    if not asym_valid:
        invalid.append("Asymmetry CI")
    if not pluri_valid:
        invalid.append("Pluripotency frac CI")
    frac_pos = frac_valid and SR.ci_positive(fwd["frac_ci_lo"], fwd["frac_ci_hi"])
    frac_incl0 = frac_valid and SR.ci_includes_zero(fwd["frac_ci_lo"], fwd["frac_ci_hi"])
    delta_neg = bool(np.isfinite(dlt) and dlt < 0)
    delta_nonneg = bool(np.isfinite(dlt) and dlt >= 0)
    n1_ok = bool(np.isfinite(main["p1"]) and main["p1"] <= SR.P_BAR)
    asym_pos = asym_valid and SR.ci_positive(a_lo, a_hi)
    asym_incl0_or_neg = asym_valid and (SR.ci_includes_zero(a_lo, a_hi) or SR.ci_negative(a_lo, a_hi))
    flags = []
    if not qualified:
        key = "pos_control_broken"
        why = "mixture control failed in this space; numbers reported as DISQUALIFIED, no reading"
    elif delta_nonneg:
        key = "no_approach_same_platform"
        why = "delta_S >= 0"
    elif not frac_valid:
        key = "mixed"
        why = "delta_S < 0 and the frac_S CI does not contain its point estimate (ci_invalid)"
    elif frac_pos and delta_neg and n1_ok and asym_pos:
        key = "toward_young_same_platform"
        why = "frac_S CI > 0, delta_S < 0, N1 p <= 0.05, Asymmetry CI > 0"
    elif frac_pos and delta_neg and n1_ok and not asym_valid:
        key = "mixed"
        why = "forward conditions met but the Asymmetry CI does not contain its point estimate (ci_invalid)"
    elif frac_pos and delta_neg and n1_ok and asym_incl0_or_neg:
        key = "convergence_not_rejuvenation"
        why = "forward conditions met, Asymmetry CI includes 0 or is negative"
    elif frac_incl0:
        key = "no_approach_same_platform"
        why = "frac_S CI includes 0"
    else:
        key = "mixed"
        why = (f"frac_S CI > 0={frac_pos}; delta_S < 0={delta_neg}; N1 p <= 0.05={n1_ok} "
               f"(p={main['p1']:.4f}); Asymmetry CI > 0={asym_pos}")
    pdl = pluri["delta"]
    if np.isfinite(pdl) and pdl < 0:
        if not pluri_valid:
            flags.append("overshoot_not_evaluable_ci_invalid")
        elif SR.ci_positive(pluri["frac_ci_lo"], pluri["frac_ci_hi"]):
            flags.append("overshoot_also_closer")
    sign = "HOLDS" if delta_nonneg else ("FLIPS" if delta_neg else "NA")
    R_S = fwd["dist_to_target"] / fwd["dist_origin_to_target"] if fwd["dist_origin_to_target"] > 0 else np.nan
    return dict(
        key=key, why=why, flags=flags, sign=sign, R_S=float(R_S), ci_invalid=invalid,
        conditions=dict(
            frac_S_CI_positive=frac_pos, frac_S_CI_includes_0=frac_incl0, delta_S_lt_0=delta_neg,
            N1_p_le_0_05=n1_ok, Asymmetry_CI_positive=asym_pos,
            Asymmetry_CI_includes_0_or_negative=asym_incl0_or_neg,
            frac_S_CI_valid=frac_valid, Asymmetry_CI_valid=asym_valid, Pluri_frac_CI_valid=pluri_valid,
        ),
    )


def same_share10(main, symbols):
    X, names = main["X"], main["names"]
    pos = {nm: i for i, nm in enumerate(names)}
    xO, xY, xS, xR = X[pos["O"]], X[pos["Y"]], X[pos["S"]], X[pos["S_rev"]]
    rows = []
    for vec, label in ((xO - xY, "start distance x(O) - x(Y)"),
                       (xS - xY, "final distance x(S) - x(Y)"),
                       (xS - xO, "displacement d_S"),
                       (xR - xO, "reverse final distance x(S_rev) - x(O)"),
                       (xR - xY, "reverse displacement d_rev")):
        sh, genes = share10(vec, symbols)
        rows.append(dict(space=main["space"], half="SAME", vector=label, share10=sh,
                         top10_genes=", ".join(genes)))
    return rows


# ----------------------------------------------------------------------------- SOUTH6
def s6_load(log):
    z = np.load(SOUTH6_DIR / "stage1.npz", allow_pickle=True)
    keep = [str(x) for x in z["keep"].tolist()]
    C_nt = np.asarray(z["C_nt"], np.float64)
    C_pool = np.asarray(z["C_pool"], np.float64)
    C_panel = np.vstack([np.asarray(z["C_f"], np.float64), C_nt, C_pool.reshape(1, -1)])
    matched_ruler = np.asarray(z["matched_ruler"], np.int64)
    nt_cells = np.asarray(z["nt_cells"], np.int64)
    nnz = int(np.asarray(z["nnz"]).ravel()[0])
    n_ruler = int(np.asarray(z["n_ruler"]).ravel()[0])
    stored_missing = np.asarray(z["missing"], bool)
    del z
    frozen = load_frozen_ruler()
    mu = np.asarray(frozen["mu"], np.float64)
    sd = np.asarray(frozen["sd"], np.float64)
    if mu.shape[0] != n_ruler:
        raise Halt("genes", f"frozen mu n={mu.shape[0]} but SOUTH6 n_ruler={n_ruler}. Not aligning.")
    tab = pd.read_csv(SOUTH6_DIR / "stage3_full_table.csv")
    t1 = tab[(tab["space"] == "S1") & (tab["pair"] == "O1_Y1")].reset_index(drop=True)
    if t1["factor"].astype(str).tolist() != keep:
        raise Halt("factors", "stage3_full_table S1 O1_Y1 factor order differs from stage1 keep. Not reordering.")
    gdis = t1["guides_disagree"].astype(bool).to_numpy()
    iloss = t1["identity_loss"].astype(bool).to_numpy()
    cpt = pd.read_csv(SOUTH6_DIR / "stage0_cells_per_target.csv").set_index("target")["n_cells"]
    n_cells = np.array([int(cpt.get(f, -1)) for f in keep], np.int64)
    symbols = np.array([str(s) for s in np.asarray(frozen["symbol"])], dtype=object)
    log(f"[south6] factors={len(keep)} matched={matched_ruler.size} ruler={n_ruler} nt_guides={C_nt.shape[0]} "
        f"nt_cells={nt_cells.size}; carried flags: guides_disagree={int(gdis.sum())} identity_loss={int(iloss.sum())}")
    return dict(
        keep=keep, C_nt=C_nt, C_panel=C_panel,
        matched_ruler=matched_ruler, nt_cells=nt_cells, nnz=nnz, n_ruler=n_ruler,
        mu=mu, sd=sd, mu_m=mu[matched_ruler], sd_m=sd[matched_ruler],
        stored_missing=stored_missing, table=t1,
        gdis=gdis, iloss=iloss, n_f=len(keep), n_nt=C_nt.shape[0],
        n_cells=n_cells, symbols=symbols,
    )


def s6_shared(S6, log):
    """Everything the three spaces share: the panel TMM, its cache, and the Stage 2B count sums."""
    logcpm, nf = TR.tmm_logcpm_quiet(S6["C_panel"])
    cache = tmm_cache_from_panel(S6["C_panel"], nf)
    missing = (S6["C_panel"] == 0).all(0)
    data_m, ind_m, ip_m = S6R.open_matched(S6["nnz"])
    n_matched = S6["matched_ruler"].size
    rng = np.random.default_rng(S6_SEED)
    cells = S6["nt_cells"]
    half = cells.size // 2
    split_counts = np.zeros((S6_N_SPLITS, 2, n_matched), np.float64)
    t0 = time.perf_counter()
    for s_i in range(S6_N_SPLITS):
        perm = rng.permutation(cells.size)
        split_counts[s_i, 0] = S6R.gather_sum(data_m, ind_m, ip_m, cells[perm[:half]], n_matched)
        split_counts[s_i, 1] = S6R.gather_sum(data_m, ind_m, ip_m, cells[perm[half:]], n_matched)
        if (s_i + 1) % 50 == 0:
            log(f"[south6 splits] {s_i + 1}/{S6_N_SPLITS} ({time.perf_counter() - t0:.0f}s)")
    del data_m, ind_m, ip_m
    return dict(logcpm=logcpm, nf=nf, cache=cache, missing=missing, split_counts=split_counts)


_PW = {}


def _perm_init(d_path, r_m, r2, rnorm, median, obs, n_perm, seed, chunk):
    _PW.update(D=np.load(d_path, mmap_mode="r"), r=np.asarray(r_m, np.float64), r2=float(r2),
               rnorm=float(rnorm), med=float(median), obs=np.asarray(obs, np.float64),
               n_perm=int(n_perm), seed=int(seed), chunk=int(chunk))


def _perm_factor(i):
    """south6_run._perm_one_medians restricted to one S1 setting. The index stream is
    default_rng([seed, i]) in chunks of PERM_CHUNK, exactly as SOUTH6 drew it."""
    D = np.asarray(_PW["D"][i], np.float64)
    p = D.shape[0]
    n_perm = _PW["n_perm"]
    rr = _PW["r"]
    r2 = _PW["r2"]
    rnorm = _PW["rnorm"]
    rng = np.random.default_rng([_PW["seed"], int(i)])
    deltas = np.empty(n_perm, np.float64)
    pooled = np.empty(n_perm, np.float32)
    s = float(np.linalg.norm(D))
    scale = (_PW["med"] / s) if s > 1e-12 else np.nan
    pos = 0
    while pos < n_perm:
        m = min(_PW["chunk"], n_perm - pos)
        idx = np.empty((m, p), np.int32)
        idx[:] = np.arange(p, dtype=np.int32)
        rng.permuted(idx, axis=1, out=idx)
        dperm = D[idx]
        qn2 = np.einsum("ij,ij->i", dperm, dperm)
        dots = dperm @ rr
        dist = np.sqrt(np.maximum(r2 + 2.0 * dots + qn2, 0.0))
        deltas[pos:pos + m] = dist - rnorm
        if np.isfinite(scale):
            dist_p = np.sqrt(np.maximum(r2 + 2.0 * scale * dots + (scale ** 2) * qn2, 0.0))
            pooled[pos:pos + m] = np.float32(dist_p - rnorm)
        else:
            pooled[pos:pos + m] = np.nan
        pos += m
    floor = float(np.percentile(deltas, 5))
    o = _PW["obs"][i]
    le = int(np.sum(deltas <= o)) if np.isfinite(o) else -1
    return int(i), floor, le, pooled


def s6_permutations(space, D_m, r_m, r2, rnorm, median, obs, workers, log):
    tmp = Path(tempfile.mkdtemp(prefix=f"coordspace_{space}_"))
    try:
        d_path = tmp / "D.npy"
        np.save(d_path, np.asarray(D_m, np.float64))
        n_f = D_m.shape[0]
        floors = np.full(n_f, np.nan)
        les = np.full(n_f, -1, np.int64)
        pooled = np.full((n_f, S6_N_PERM), np.nan, np.float32)
        t0 = time.perf_counter()
        ctx = mp.get_context("spawn")
        with ctx.Pool(
            processes=int(workers), initializer=_perm_init,
            initargs=(str(d_path), r_m, r2, rnorm, median, obs, S6_N_PERM, S6_SEED, S6_PERM_CHUNK),
        ) as pool:
            done = 0
            for i, fl, le, pv in pool.imap_unordered(_perm_factor, range(n_f), chunksize=2):
                floors[i] = fl
                les[i] = le
                pooled[i] = pv
                done += 1
                if done % 100 == 0 or done == n_f:
                    log(f"[perm {space}] {done}/{n_f} ({time.perf_counter() - t0:.0f}s)")
        vv = pooled[np.isfinite(pooled)]
        pooled_floor = float(np.percentile(vv, 5)) if vv.size else np.nan
        return floors, les, pooled_floor
    finally:
        shutil.rmtree(tmp, ignore_errors=True)


def s6_space(space, S6, SH, c_old, c_young, workers, log):
    """SOUTH6 S1 O1->Y1 in one space: Stage 1 vectors, 2A, 2B (78 guides), 2C, Stage 3 labels."""
    mr = S6["matched_ruler"]
    n_ruler = S6["n_ruler"]
    n_f, n_nt = S6["n_f"], S6["n_nt"]
    missing = SH["missing"]
    logcpm, nf = SH["logcpm"], SH["nf"]
    mu_m, sd_m = S6["mu_m"], S6["sd_m"]
    if space == "C":
        sd_safe = np.where(sd_m < 1e-12, 1.0, sd_m)
        X_m = (logcpm - mu_m) / sd_safe
        X_m[:, missing] = 0.0
        fill = np.zeros(n_ruler)
    elif space == "B":
        X_m = np.array(logcpm, np.float64, copy=True)
        X_m[:, missing] = mu_m[missing]
        fill = S6["mu"].copy()
    else:
        lib = S6["C_panel"].sum(1) * np.asarray(nf, np.float64)
        X_m = S6["C_panel"] / lib[:, None] * 1e6
        X_m[:, missing] = 0.0
        fill = np.zeros(n_ruler)
    for i in range(3):
        xi = x_one(S6["C_panel"][i], SH["cache"], mu_m, sd_m, missing, space)
        e = float(np.max(np.abs(xi - X_m[i])) / max(1.0, float(np.max(np.abs(X_m[i])))))
        if e > 1e-4:
            raise Halt("tmm", f"space {space}: single-row transform disagrees with the panel on row {i} ({e:.3e}).")

    def full(xm_rows):
        xm_rows = np.atleast_2d(xm_rows)
        out = np.tile(fill, (xm_rows.shape[0], 1))
        out[:, mr] = xm_rows
        return out

    x_pool = full(X_m[-1])[0]
    d_f = full(X_m[:n_f])
    d_f -= x_pool
    del X_m
    r = np.asarray(c_old, np.float64) - np.asarray(c_young, np.float64)
    r2 = float(np.dot(r, r))
    rnorm = float(np.linalg.norm(r))
    delta, cos, frac, dn = S6R.block_delta(d_f, r, r2, rnorm)
    dlt_m, *_ = S6R.block_delta(d_f[:, mr], r[mr], r2, rnorm)
    if float(np.max(np.abs(dlt_m - delta))) > 1e-6 * max(1.0, float(np.max(np.abs(delta)))):
        raise Halt("algebra", f"space {space}: matched-coordinate delta drifted from the full-vector delta.")
    median = float(np.median(dn))
    log(f"[south6 {space}] ||O1-Y1||={rnorm:.6f} median ||d||={median:.6f} n_delta<0={int((delta < 0).sum())}")
    floors, les, pooled_floor = s6_permutations(space, d_f[:, mr], r[mr], r2, rnorm, median, delta, workers, log)

    def xfull_one(counts):
        return full(x_one(counts, SH["cache"], mu_m, sd_m, missing, space))[0]

    C_nt = S6["C_nt"]
    tot = C_nt.sum(0)
    loo = np.vstack([xfull_one(C_nt[i]) - xfull_one(tot - C_nt[i]) for i in range(n_nt)])
    spl = np.vstack([xfull_one(SH["split_counts"][s, 0]) - xfull_one(SH["split_counts"][s, 1])
                     for s in range(S6_N_SPLITS)])
    gd, *_ = S6R.block_delta(loo, r, r2, rnorm)
    sdl, *_ = S6R.block_delta(spl, r, r2, rnorm)
    gf78 = float(np.percentile(gd, 5))
    sf78 = float(np.percentile(sdl, 5))
    u_full, _ = unit(np.asarray(c_young, np.float64) - np.asarray(c_old, np.float64))
    doses = S6R.dose_deltas(u_full, median, r, r2, rnorm)
    stricter = min(pooled_floor, gf78, sf78)
    f_pri, d_pri = S6R.mda_of(doses, pooled_floor)
    f_str, d_str = S6R.mda_of(doses, stricter)
    p = np.where(np.isfinite(delta), (les + 1.0) / (S6_N_PERM + 1.0), np.nan)
    q = bh_q(p)
    toward = (delta < 0) & (delta < floors) & np.isfinite(q) & (q <= S6_Q_BAR) & ~S6["iloss"] & ~S6["gdis"]
    n_neg_coords = None
    n_neg_factors = None
    if space == "A":
        n_neg_coords = 0
        n_neg_factors = 0
        co = np.asarray(c_old, np.float64)[None, :]
        for a in range(0, n_f, 200):
            neg = (co + d_f[a:a + 200]) < 0
            n_neg_coords += int(neg.sum())
            n_neg_factors += int(neg.any(axis=1).sum())
    sh_gap, sh_gap_genes = share10(-r, S6["symbols"])
    sh_d = share10_rows(d_f)
    sh_final = share10_rows(d_f, offset=r)
    out = dict(
        space=space, d_f=d_f, r=r, rnorm=rnorm, delta=delta, cos=cos, frac=frac, d_norm=dn,
        median=median, floors=floors, les=les, pooled_floor=pooled_floor,
        nt_guide_floor=gf78, split_floor=sf78, stricter_floor=stricter,
        stricter_source=("pooled" if stricter == pooled_floor else ("nt_guide" if stricter == gf78 else "split")),
        doses=doses, mda_primary=f_pri, delta_at_mda_primary=d_pri,
        mda_stricter=f_str, delta_at_mda_stricter=d_str,
        mda_margin=(stricter - d_str) if np.isfinite(f_str) else np.nan,
        p=p, q=q, toward=toward, n_neg_coords=n_neg_coords, n_neg_factors=n_neg_factors,
        keep=list(S6["keep"]), n_cells=S6["n_cells"],
        share10_gap=sh_gap, share10_gap_genes=sh_gap_genes,
        share10_d=sh_d, share10_final=sh_final,
        counts=dict(
            n_delta_negative=int((delta < 0).sum()), n_below_floor=int((delta < floors).sum()),
            n_q=int((q <= S6_Q_BAR).sum()), n_toward_young=int(toward.sum()),
        ),
    )
    log(f"[south6 {space}] pooled floor={pooled_floor:.6f} MDA primary={f_pri} stricter={f_str} "
        f"toward_young={int(toward.sum())}")
    return out


# ----------------------------------------------------------------------------- Stage 0
def stage0(man, workers, log):
    rows = []
    log("=" * 90)
    log("STAGE 0 — reproduce space C before any A or B number")
    log("=" * 90)

    # (a) SAME
    S = same_load(log)
    cmp_cat(rows, "SAME", "stored split equals a fresh split_half draw", [S["split_matches"]], [True], gating=False)
    pdf, psum = same_posctrl("C", S, log)
    st_p = pd.read_csv(SR.SAME_DIR / "posctrl.csv")
    for col in ("f", "n_aged", "n_young", "n_total", "frac", "delta", "cos", "d_norm", "v_norm",
                "frac_ci_lo", "frac_ci_hi", "delta_ci_lo", "delta_ci_hi", "ruler_score_context"):
        cmp_num(rows, "SAME posctrl.csv", col, pdf[col].to_numpy(), st_p[col].to_numpy(), TOL_SAME)
    st_ps = SR.load_json(SR.SAME_DIR / "posctrl_summary.json")
    cmp_cat(rows, "SAME posctrl_summary.json", "P pass", [psum["pass_same_rule"]], [bool(st_ps["pass"])])
    pdf_C = pdf
    anc = np.load(TR.ANCHORS_PATH, allow_pickle=True)
    main = same_main("C", S, np.asarray(anc["u"], np.float64), log)
    st_s = pd.read_csv(SR.SAME_DIR / "stats.csv")
    new_s = main["stats_df"]
    for col in ("n_cells_S", "cos", "frac", "delta", "p_N1", "frac_ci_lo", "frac_ci_hi",
                "delta_ci_lo", "delta_ci_hi", "d_norm", "v_norm", "ruler_score_S"):
        cmp_num(rows, "SAME stats.csv", col, new_s[col].to_numpy(), st_s[col].to_numpy(), TOL_SAME)
    st_sum = SR.load_json(SR.SAME_DIR / "summary.json")
    cmp_num(rows, "SAME summary.json", "asymmetry", [main["asym"]], [st_sum["asymmetry"]], TOL_SAME)
    cmp_num(rows, "SAME summary.json", "asymmetry_ci_lo", [main["asym_ci"][0]], [st_sum["asymmetry_ci_lo"]], TOL_SAME)
    cmp_num(rows, "SAME summary.json", "asymmetry_ci_hi", [main["asym_ci"][1]], [st_sum["asymmetry_ci_hi"]], TOL_SAME)
    st_c = pd.read_csv(SR.SAME_DIR / "context.csv").set_index("pseudobulk")
    new_c = main["ctx_df"].set_index("pseudobulk")
    for col in ("n_cells", "ruler_score_context", "frac", "delta", "cos",
                "frac_ci_lo", "frac_ci_hi", "delta_ci_lo", "delta_ci_hi"):
        cmp_num(rows, "SAME context.csv", col, new_c.loc[list(new_c.index), col].to_numpy(),
                st_c.loc[list(new_c.index), col].to_numpy(), TOL_SAME)
    cmp_num(rows, "SAME context.csv", "ruler_score_context (O, Y)",
            [main["ruler"]["O"], main["ruler"]["Y"]],
            [st_c.loc["O", "ruler_score_context"], st_c.loc["Y", "ruler_score_context"]], TOL_SAME)
    cmp_cat(rows, "SAME reading", "key", [main["reading_same_rule"]["key"]], [st_sum["key"]])
    cmp_cat(rows, "SAME reading", "key is no_approach_same_platform",
            [main["reading_same_rule"]["key"]], ["no_approach_same_platform"])
    cmp_num(rows, "SAME reading", "delta_S (+130.894)", [round(main["fwd"]["delta"], 3)], [130.894], TOL_SAME)
    same_c = dict(
        delta_S=main["fwd"]["delta"], frac_S=main["fwd"]["frac"], key=main["reading_same_rule"]["key"],
        asym=main["asym"], asym_ci=list(main["asym_ci"]), p_pass=psum["pass_same_rule"],
        p_pass_coordspace=psum["pass_coordspace"],
    )
    main_C = main
    del S

    # (c) GTEx centroids
    gtex = TR.load_gtex_z(log)
    age = gtex["age"]
    young_m = np.isin(age, TR.YOUNG_BINS)
    old_m = np.isin(age, TR.OLD_BINS)
    cmp_num(rows, "GTEx centroids", "n_young (113), n_old (231)",
            [young_m.sum(), old_m.sum()], [int(anc["n_young"]), int(anc["n_old"])], 0.0)
    frozen = gtex["frozen"]
    mu = np.asarray(frozen["mu"], np.float64)
    sd_safe = np.where(np.asarray(frozen["sd"], np.float64) < 1e-12, 1.0, np.asarray(frozen["sd"], np.float64))
    Xg = np.asarray(gtex["pack"]["X"], np.float64)
    cB_y = Xg[young_m].mean(0)
    cB_o = Xg[old_m].mean(0)
    c_y = np.asarray(anc["c_young"], np.float64)
    c_o = np.asarray(anc["c_old"], np.float64)
    cmp_num(rows, "GTEx centroids", "B young vs mu + sd * C young", cB_y, mu + sd_safe * c_y, TOL_SAME)
    cmp_num(rows, "GTEx centroids", "B old vs mu + sd * C old", cB_o, mu + sd_safe * c_o, TOL_SAME)
    Zg = np.asarray(gtex["Z"], np.float64)
    cmp_num(rows, "GTEx centroids", "C young recomputed vs anchors.npz (diagnostic)",
            Zg[young_m].mean(0), c_y, TOL_SAME, gating=False)
    cmp_num(rows, "GTEx centroids", "C old recomputed vs anchors.npz (diagnostic)",
            Zg[old_m].mean(0), c_o, TOL_SAME, gating=False)
    del gtex, Xg, Zg, cB_y, cB_o

    # (b) SOUTH6, S1 O1_Y1
    S6 = s6_load(log)
    SH = s6_shared(S6, log)
    cmp_cat(rows, "SOUTH6 stage1.npz", "panel missing-gene mask", SH["missing"], S6["stored_missing"], gating=False)
    res = s6_space("C", S6, SH, c_o, c_y, workers, log)
    stored_d_f = np.asarray(np.load(SOUTH6_DIR / "stage1.npz", allow_pickle=True)["d_f"], np.float64)
    cmp_num(rows, "SOUTH6 stage1.npz", "d_g, all factors and genes (diagnostic)",
            res["d_f"], stored_d_f, TOL_SOUTH6, gating=False)
    del stored_d_f
    otd = pd.read_csv(SOUTH6_DIR / "origin_target_distances.csv")
    row = otd[(otd["space"] == "S1") & (otd["origin"] == "O1") & (otd["target"] == "Y1")].iloc[0]
    cmp_num(rows, "SOUTH6 origin_target_distances.csv", "||O1 - Y1||", [res["rnorm"]], [row["distance"]], TOL_SOUTH6)
    fl = np.load(SOUTH6_DIR / "stage2a_floors.npz")
    fsum = pd.read_csv(SOUTH6_DIR / "stage2a_floor_summary.csv")
    if not (fsum.iloc[S6_PRIMARY_SETTING_INDEX]["space"] == "S1" and fsum.iloc[S6_PRIMARY_SETTING_INDEX]["pair"] == "O1_Y1"):
        raise Halt("setting_order", "stage2a_floor_summary row 0 is not S1 O1_Y1. Not guessing the index.")
    k = S6_PRIMARY_SETTING_INDEX
    cmp_num(rows, "SOUTH6 stage2a_floors.npz", "per-factor floor", res["floors"], fl["floors"][k], TOL_SOUTH6)
    cmp_num(rows, "SOUTH6 stage2a_floors.npz", "per-factor permutation count le", res["les"], fl["les"][k], 0.0)
    cmp_num(rows, "SOUTH6 stage2a_floors.npz", "pooled floor", [res["pooled_floor"]], [fl["pooled"][k]], TOL_SOUTH6)
    cmp_num(rows, "SOUTH6 stage2a_floors.npz", "median ||d||", [res["median"]], [fl["medians"][k]], TOL_SOUTH6)
    mda = pd.read_csv(SOUTH6_DIR / "stage2c_mda.csv")
    mrow = mda[(mda["space"] == "S1") & (mda["pair"] == "O1_Y1")].iloc[0]
    for col, val in (("mda_primary", res["mda_primary"]), ("delta_at_mda_primary", res["delta_at_mda_primary"]),
                     ("pooled_floor", res["pooled_floor"]), ("mda_stricter", res["mda_stricter"]),
                     ("delta_at_mda_stricter", res["delta_at_mda_stricter"]),
                     ("stricter_floor", res["stricter_floor"]), ("mda_margin", res["mda_margin"]),
                     ("median_norm", res["median"]), ("nt_guide_floor", res["nt_guide_floor"]),
                     ("split_floor", res["split_floor"])):
        cmp_num(rows, "SOUTH6 stage2c_mda.csv", col, [val], [mrow[col]], TOL_SOUTH6)
    cmp_cat(rows, "SOUTH6 stage2c_mda.csv", "stricter_source", [res["stricter_source"]], [mrow["stricter_source"]])
    t1 = S6["table"]
    for col, val in (("delta", res["delta"]), ("cos", res["cos"]), ("frac", res["frac"]),
                     ("d_norm", res["d_norm"]), ("floor_p05", res["floors"]),
                     ("p_delta", res["p"]), ("q_delta", res["q"])):
        cmp_num(rows, "SOUTH6 stage3_full_table.csv (S1 O1_Y1)", col, val, t1[col].to_numpy(), TOL_SOUTH6)
    cmp_cat(rows, "SOUTH6 stage3_full_table.csv (S1 O1_Y1)", "below floor (per factor)",
            res["delta"] < res["floors"], (t1["delta"] < t1["floor_p05"]).to_numpy())
    cmp_cat(rows, "SOUTH6 stage3_full_table.csv (S1 O1_Y1)", "q <= 0.05 (per factor)",
            res["q"] <= S6_Q_BAR, (t1["q_delta"] <= S6_Q_BAR).to_numpy())
    cmp_cat(rows, "SOUTH6 stage3_full_table.csv (S1 O1_Y1)", "toward_young (per factor)",
            res["toward"], t1["toward_young"].astype(bool).to_numpy())
    south6_c = dict(
        rnorm=res["rnorm"], median=res["median"], pooled_floor=res["pooled_floor"],
        mda_primary=res["mda_primary"], mda_stricter=res["mda_stricter"], counts=res["counts"],
        floor_min=float(np.min(res["floors"])), floor_median=float(np.median(res["floors"])),
        floor_max=float(np.max(res["floors"])),
    )
    res.pop("d_f", None)
    s6_C = res
    del S6, SH

    rdf = pd.DataFrame(rows).rename(columns={"pass_": "pass"})
    rdf.to_csv(COORD_DIR / "stage0_reproduction.csv", index=False)
    gating = rdf[rdf["gating"]]
    failed = gating[~gating["pass"]]
    diag_failed = rdf[(~rdf["gating"]) & (~rdf["pass"])]
    passed = bool(len(failed) == 0)
    summary = dict(
        passed=passed, n_checks=int(len(rdf)), n_gating=int(len(gating)),
        n_gating_failed=int(len(failed)), n_diagnostic_failed=int(len(diag_failed)),
        tol_same=TOL_SAME, tol_south6=TOL_SOUTH6, same_C=same_c, south6_C=south6_c,
        failed_items=[f"{r.block} :: {r.item}" for r in failed.itertuples()],
        diagnostic_failed_items=[f"{r.block} :: {r.item}" for r in diag_failed.itertuples()],
    )
    dump_json(COORD_DIR / "stage0_summary.json", jsonable(summary))
    man["stages"]["stage0"] = summary
    save_manifest(man)
    log(f"[stage0] gating checks {len(gating) - len(failed)}/{len(gating)} passed; "
        f"diagnostic failures {len(diag_failed)}")
    if not passed:
        raise Halt(
            "space_c_not_reproduced",
            "Space C did not reproduce the stored SAME/SOUTH6 results: " + "; ".join(summary["failed_items"]),
            summary,
        )
    return dict(summary=summary, pdf_C=pdf_C, psum_C=psum, main_C=main_C, s6_C=s6_C)


# ----------------------------------------------------------------------------- Stage 1
def stage1(man, ctx, prereg, log):
    log("=" * 90)
    log("STAGE 1 — mixture control in every space, before any PartialReprog vector in A or B")
    log("=" * 90)
    S = same_load(log)
    tables = {"C": ctx["pdf_C"]}
    sums = {"C": ctx["psum_C"]}
    for sp in ("A", "B"):
        pdf, ps = same_posctrl(sp, S, log)
        tables[sp] = pdf
        sums[sp] = ps
    ptab = pd.concat([tables[sp].assign(space=sp) for sp in SPACES], ignore_index=True)
    ptab.to_csv(COORD_DIR / "stage1_posctrl.csv", index=False)
    qual = {sp: bool(sums[sp]["pass_coordspace"]) for sp in SPACES}
    qrows = []
    for sp in SPACES:
        s = sums[sp]
        qrows.append(dict(
            space=sp, qualified=qual[sp], monotonic=s["monotonic"],
            frac50_ci_excludes_0=s["frac50_ci_excludes_0"], frac50_ci_valid=s["frac50_ci_valid"],
            delta50_lt_0=s["delta50_neg"], pass_same_rule=s["pass_same_rule"],
        ))
    qdf = pd.DataFrame(qrows)
    qdf.to_csv(COORD_DIR / "stage1_qualification.csv", index=False)
    block = ["", "--- QUALIFICATION, per space (mixture control P) ---"]
    for r in qrows:
        block.append(
            f"space {r['space']}: {'QUALIFIED' if r['qualified'] else 'DISQUALIFIED'}; "
            f"frac monotonic in f={r['monotonic']}; frac CI at f=0.50 excludes 0={r['frac50_ci_excludes_0']}; "
            f"that CI contains its point estimate={r['frac50_ci_valid']}; delta at f=0.50 < 0={r['delta50_lt_0']}"
        )
    qpath = COORD_DIR / "QUALIFIED_BEFORE_PRIMARY.flag"
    body = "\n".join(block) + "\n"
    if qpath.exists():
        old = qpath.read_text(encoding="utf-8")
        if body.strip() not in old:
            raise Halt("qualified_flag_mismatch",
                       "QUALIFIED_BEFORE_PRIMARY.flag exists with different qualification. Not rewriting it.")
        log("[stage1] QUALIFIED_BEFORE_PRIMARY.flag exists with the same qualification; not rewritten")
    else:
        stamp = time.strftime("%Y-%m-%dT%H:%M:%S%z")
        stamp = stamp[:-2] + ":" + stamp[-2:]
        head = (f"Written after Stage 1 and before any PartialReprog vector in space A or B, {stamp}.\n"
                "Not rewritten afterwards.\n\n")
        qpath.write_text(head + prereg + "\n" + body, encoding="utf-8")
        log(f"[stage1] wrote {qpath}")
    man["stages"]["stage1"] = dict(qualified=qual, conditions=qrows)
    save_manifest(man)
    ctx.update(S=S, p_tables=tables, p_sums=sums, qualified=qual)
    return ctx


# ----------------------------------------------------------------------------- Stage 2
def stage2(man, ctx, symbols, log):
    log("=" * 90)
    log("STAGE 2 — SAME in every space")
    log("=" * 90)
    S = ctx["S"]
    mains = {"C": ctx["main_C"]}
    for sp in ("A", "B"):
        c_old, c_young = ctx["centroids"][sp]
        mains[sp] = same_main(sp, S, c_young - c_old, log)
    readings = {sp: coord_reading(mains[sp], ctx["qualified"][sp]) for sp in SPACES}
    srows, share_rows = [], []
    for sp in SPACES:
        m, rd = mains[sp], readings[sp]
        srows.append(dict(
            space=sp, qualified=ctx["qualified"][sp], key=rd["key"], why=rd["why"],
            flags=";".join(rd["flags"]), sign=rd["sign"], ci_invalid=";".join(rd["ci_invalid"]),
            n_cells_S=m["n"]["S"], delta_S=m["fwd"]["delta"], frac_S=m["fwd"]["frac"], cos_S=m["fwd"]["cos"],
            R_S=rd["R_S"], dist_start=m["fwd"]["dist_origin_to_target"], dist_final=m["fwd"]["dist_to_target"],
            p_N1=m["p1"], frac_S_ci_lo=m["fwd"]["frac_ci_lo"], frac_S_ci_hi=m["fwd"]["frac_ci_hi"],
            delta_S_ci_lo=m["fwd"]["delta_ci_lo"], delta_S_ci_hi=m["fwd"]["delta_ci_hi"],
            frac_rev=m["rev"]["frac"], delta_rev=m["rev"]["delta"], p_N1_rev=m["p1r"],
            frac_rev_ci_lo=m["rev"]["frac_ci_lo"], frac_rev_ci_hi=m["rev"]["frac_ci_hi"],
            asymmetry=m["asym"], asymmetry_ci_lo=m["asym_ci"][0], asymmetry_ci_hi=m["asym_ci"][1],
            pluri_frac=m["pluri"]["frac"], pluri_delta=m["pluri"]["delta"],
            nonreprog_frac=m["nr"]["frac"], nonreprog_delta=m["nr"]["delta"],
            cos_v_u_toward=m["cos_v_u"], **{f"cond_{k}": v for k, v in rd["conditions"].items()},
        ))
        share_rows += same_share10(m, symbols)
    sdf = pd.DataFrame(srows)
    sdf.to_csv(COORD_DIR / "stage2_same_summary.csv", index=False)
    pd.concat([mains[sp]["stats_df"] for sp in SPACES], ignore_index=True).to_csv(
        COORD_DIR / "stage2_same_stats.csv", index=False)
    pd.concat([mains[sp]["ctx_df"] for sp in SPACES], ignore_index=True).to_csv(
        COORD_DIR / "stage2_same_context.csv", index=False)
    man["stages"]["stage2"] = dict(
        readings={sp: dict(key=readings[sp]["key"], sign=readings[sp]["sign"], R_S=readings[sp]["R_S"],
                           flags=readings[sp]["flags"], ci_invalid=readings[sp]["ci_invalid"])
                  for sp in SPACES})
    save_manifest(man)
    for sp in SPACES:
        log(f"[stage2 {sp}] qualified={ctx['qualified'][sp]} key={readings[sp]['key']} "
            f"sign={readings[sp]['sign']} R_S={readings[sp]['R_S']:.4f}")
    ctx.update(mains=mains, readings=readings, same_share_rows=share_rows)
    del ctx["S"]
    return ctx


# ----------------------------------------------------------------------------- Stage 3
def stage3(man, ctx, workers, log):
    log("=" * 90)
    log("STAGE 3 — SOUTH6 S1 O1->Y1 in every space")
    log("=" * 90)
    S6 = s6_load(log)
    SH = s6_shared(S6, log)
    res = {"C": ctx["s6_C"]}
    for sp in ("A", "B"):
        c_old, c_young = ctx["centroids"][sp]
        r = s6_space(sp, S6, SH, c_old, c_young, workers, log)
        r.pop("d_f", None)
        res[sp] = r
    del S6, SH
    keep = res["C"]["keep"]
    toward_C = set(np.asarray(keep)[res["C"]["toward"]].tolist())
    sum_rows, fac_rows, out_rows = [], [], []
    for sp in SPACES:
        r = res[sp]
        ndp = not np.isfinite(r["mda_primary"])
        tw = sorted(np.asarray(keep)[r["toward"]].tolist())
        if not ctx["qualified"][sp]:
            outcome = "DISQUALIFIED"
        elif ndp:
            outcome = "no_detection_power"
        elif sp == "C":
            outcome = "reference"
        else:
            outcome = "UNCHANGED" if set(tw) == toward_C else "CHANGED"
        sum_rows.append(dict(
            space=sp, qualified=ctx["qualified"][sp], outcome=outcome, gap=r["rnorm"],
            median_d_norm=r["median"], floor_min=float(np.min(r["floors"])),
            floor_median=float(np.median(r["floors"])), floor_max=float(np.max(r["floors"])),
            pooled_floor=r["pooled_floor"], median_floor_over_gap=float(np.median(r["floors"])) / r["rnorm"],
            nt_guide_floor=r["nt_guide_floor"], split_floor=r["split_floor"],
            stricter_floor=r["stricter_floor"], stricter_source=r["stricter_source"],
            mda_primary=r["mda_primary"], mda_stricter=r["mda_stricter"], mda_margin=r["mda_margin"],
            no_detection_power=ndp, **r["counts"],
            n_tie_at_min_p=int(np.sum(np.abs(r["p"] - 1.0 / (S6_N_PERM + 1.0)) < 1e-12)),
            n_negative_coords_A=r["n_neg_coords"], n_factors_with_negative_coords_A=r["n_neg_factors"],
            toward_young_factors=";".join(tw),
        ))
        for i, f in enumerate(keep):
            fac_rows.append(dict(
                space=sp, factor=f, n_cells=int(r["n_cells"][i]), delta=r["delta"][i], cos=r["cos"][i],
                frac=r["frac"][i], d_norm=r["d_norm"][i], floor_p05=r["floors"][i],
                p_delta=r["p"][i], q_delta=r["q"][i], toward_young=bool(r["toward"][i]),
                share10_d=r["share10_d"][i], share10_final=r["share10_final"][i],
            ))
        if tw:
            for f in tw:
                i = keep.index(f)
                out_rows.append(dict(space=sp, factor=f, delta=r["delta"][i], floor_p05=r["floors"][i],
                                     q_delta=r["q"][i], n_cells=int(r["n_cells"][i])))
    sdf = pd.DataFrame(sum_rows)
    sdf.to_csv(COORD_DIR / "stage3_south6_summary.csv", index=False)
    pd.DataFrame(fac_rows).to_csv(COORD_DIR / "stage3_south6_factors.csv", index=False)
    pd.DataFrame(out_rows, columns=["space", "factor", "delta", "floor_p05", "q_delta", "n_cells"]).to_csv(
        COORD_DIR / "stage3_toward_young.csv", index=False)
    from scipy.stats import spearmanr
    rk = []
    for a, b in (("A", "C"), ("B", "C"), ("A", "B")):
        rho = float(spearmanr(res[a]["delta"], res[b]["delta"]).correlation)
        ta = set(np.asarray(keep)[np.argsort(res[a]["delta"], kind="stable")[:20]].tolist())
        tb = set(np.asarray(keep)[np.argsort(res[b]["delta"], kind="stable")[:20]].tolist())
        wa = set(np.asarray(keep)[res[a]["toward"]].tolist())
        wb = set(np.asarray(keep)[res[b]["toward"]].tolist())
        rk.append(dict(pair=f"{a}-{b}", spearman_delta=rho, top20_overlap=len(ta & tb),
                       toward_young_overlap=len(wa & wb), toward_young_a=len(wa), toward_young_b=len(wb)))
    pd.DataFrame(rk).to_csv(COORD_DIR / "stage3_rank_agreement.csv", index=False)
    man["stages"]["stage3"] = dict(summary=jsonable(sum_rows), rank_agreement=rk)
    save_manifest(man)
    for r in sum_rows:
        log(f"[stage3 {r['space']}] qualified={r['qualified']} outcome={r['outcome']} "
            f"toward_young={r['n_toward_young']} MDA={r['mda_primary']}")
    ctx.update(s6=res)
    return ctx


# ----------------------------------------------------------------------------- Stage 4
def stage4(man, ctx, log):
    rows = list(ctx["same_share_rows"])
    for sp in SPACES:
        r = ctx["s6"][sp]
        rows.append(dict(space=sp, half="SOUTH6", vector="O1 - Y1 (gap)", share10=r["share10_gap"],
                         top10_genes=", ".join(r["share10_gap_genes"])))
        for lab, arr in (("d_g across factors", r["share10_d"]),
                         ("(O1 + d_g) - Y1 across factors", r["share10_final"])):
            a = np.asarray(arr, float)
            a = a[np.isfinite(a)]
            rows.append(dict(space=sp, half="SOUTH6", vector=lab,
                             share10=float(np.median(a)) if a.size else np.nan,
                             p05=float(np.percentile(a, 5)) if a.size else np.nan,
                             p95=float(np.percentile(a, 95)) if a.size else np.nan, top10_genes=""))
    df = pd.DataFrame(rows)
    df.to_csv(COORD_DIR / "stage4_share10.csv", index=False)
    man["stages"]["stage4"] = dict(n_rows=int(len(df)))
    save_manifest(man)
    return ctx


# ----------------------------------------------------------------------------- decision
def _share_final(ctx, sp):
    for r in ctx["same_share_rows"]:
        if r["space"] == sp and r["vector"].startswith("final distance"):
            return r["share10"], r["top10_genes"]
    return np.nan, ""


def decide(man, ctx, log):
    qual, rd, mains = ctx["qualified"], ctx["readings"], ctx["mains"]
    flips = {sp: rd[sp]["sign"] == "FLIPS" for sp in SPACES}
    row3 = [sp for sp in SPACES if qual[sp] and flips[sp]]
    row2 = [sp for sp in SPACES if (not qual[sp]) and flips[sp]]
    surviving = [sp for sp in SPACES if qual[sp]]
    disq = [sp for sp in SPACES if not qual[sp]]
    fired = [3] if row3 else ([1] + ([2] if row2 else []))
    trivial = not (qual["A"] or qual["B"])
    words = []
    if 1 in fired:
        parts = []
        for sp in surviving:
            sh, genes = _share_final(ctx, sp)
            parts.append(f"{sp}: R_S = {rd[sp]['R_S']:.3f}, top-10 genes carry {sh:.1%} of the squared final distance")
        w = ("In every coordinate space that detected known young-cell mixtures ("
             + ", ".join(surviving) + "), aged partially reprogrammed cells end farther from the young donor "
             "than they started (" + "; ".join(parts) + "). The same-platform result does not depend on the "
             "log or the z-score.")
        if disq:
            w += " Disqualified by the mixture control: " + ", ".join(disq) + "."
        if trivial:
            w += (" No alternative space could detect a known approach, so coordinate robustness was not "
                  "tested.")
        words.append(dict(row=1, text=w))
    if 2 in fired:
        for sp in row2:
            words.append(dict(row=2, text=(
                f"Space {sp} ends with delta_S < 0 but failed its mixture control. Its numbers are reported "
                "as disqualified and support no conclusion.")))
    if 3 in fired:
        for sp in row3:
            key = rd[sp]["key"]
            sh, _ = _share_final(ctx, sp)
            if key == "toward_young_same_platform":
                t = (f"In space {sp}, which detected known young-cell mixtures, aged partially reprogrammed "
                     "cells end closer to the young donor than they started and approach it more than the "
                     "reverse. The same-platform result depends on the choice of coordinates.")
            elif key == "convergence_not_rejuvenation":
                t = (f"In space {sp}, which detected known young-cell mixtures, aged partially reprogrammed "
                     "cells end closer to the young donor than they started, but the reverse test does not "
                     "separate this from mutual convergence of the two donors. The statement that the cells "
                     "end farther from young depends on the choice of coordinates.")
            else:
                t = (f"In space {sp}, which detected known young-cell mixtures, the distance to the young "
                     f"donor falls (delta_S < 0), but {rd[sp]['why']}. The sign of the change in distance "
                     "depends on the choice of coordinates.")
            t += f" (R_S = {rd[sp]['R_S']:.3f}; top-10 genes carry {sh:.1%} of the squared final distance.)"
            words.append(dict(row=3, text=t))
        stamp = time.strftime("%Y-%m-%dT%H:%M:%S%z")
        (COORD_DIR / "REVISION_TRIGGERED.flag").write_text(
            f"Decision-table row 3 fired, {stamp[:-2] + ':' + stamp[-2:]}.\n"
            + "\n".join(f"space {sp}: key {rd[sp]['key']}, delta_S {mains[sp]['fwd']['delta']:.6f}" for sp in row3)
            + "\nSpace C stays primary. Report to the investigator before anything else is run.\n",
            encoding="utf-8")
    s6 = man["stages"].get("stage3", {}).get("summary", [])
    screen = []
    for r in s6:
        if r["space"] == "C":
            continue
        if r["outcome"] == "CHANGED":
            screen.append(f"Screen, space {r['space']}: CHANGED. toward_young = {r['n_toward_young']} "
                          f"({r['toward_young_factors']}) against 0 in space C. Reported immediately; "
                          "SOUTH6 space C remains the screen's primary.")
        elif r["outcome"] == "UNCHANGED":
            screen.append(f"Screen, space {r['space']}: UNCHANGED. toward_young = {r['n_toward_young']}, "
                          "the same set as space C.")
        else:
            screen.append(f"Screen, space {r['space']}: {r['outcome']}.")
    dec = dict(rows_fired=fired, row3_spaces=row3, row2_spaces=row2, surviving=surviving,
               disqualified=disq, trivial=trivial, wording=words, screen=screen)
    dump_json(COORD_DIR / "decision.json", jsonable(dec))
    man["stages"]["decision"] = dec
    save_manifest(man)
    log(f"[decision] rows fired={fired} surviving={surviving} disqualified={disq}")
    return dec


# ----------------------------------------------------------------------------- reports
def _tbl(path, cols=None):
    if not path.exists():
        return None
    df = pd.read_csv(path)
    if cols:
        df = df[[c for c in cols if c in df.columns]]
    return md_table(df)


def write_reports(man, prereg, status, next_action):
    st = man["stages"]
    s0 = st.get("stage0")
    lines = [
        "# FINDINGS_COORDSPACE — coordinate-space robustness check (SECONDARY)",
        "",
        f"**Status:** {status}",
        "",
        "Reproduced by `src/coordspace_run.py`. Imports `src/same_run.py`, `src/south6_run.py`, "
        "`src/toward_run.py`, `src/seng_run.py` and `src/gtex_common.py` without modifying them. "
        "The frozen ruler is not refit. Space C is the pre-registered primary and is not replaced by "
        "anything below.",
        "",
    ]
    dec = st.get("decision")
    if dec:
        if 3 in dec["rows_fired"]:
            lines += ["**REVISION_TRIGGERED** (`results/coordspace/REVISION_TRIGGERED.flag`).", ""]
        lines += ["## Decision-table outcome", ""]
        lines.append(f"Rows fired: {', '.join(str(r) for r in dec['rows_fired'])}. "
                     f"Surviving (qualified) spaces: {', '.join(dec['surviving']) or 'none'}. "
                     f"Disqualified: {', '.join(dec['disqualified']) or 'none'}.")
        lines.append("")
        for w in dec["wording"]:
            lines.append(f"- Row {w['row']}: {w['text']}")
        for s in dec["screen"]:
            lines.append(f"- {s}")
        lines.append("")
    lines += [
        "## Pre-registration (verbatim, written before any COORDSPACE statistic)",
        "",
        "```",
        prereg,
        "```",
        "",
        "## Notes on the pre-registration text",
        "",
        "The flag is a timestamped record and was not edited. Two statements in it are factually "
        "wrong about SOUTH6. Neither changes any computation.",
        "",
    ]
    for n in man.get("prereg_text_notes", []):
        lines.append(f"- **{n['item']}.** The flag says: {n['flag_says']}. Actual: {n['actual']} "
                     f"Effect: {n['effect']}")
    lines.append("")
    lines += ["## Stage 0 — reproduction of space C", ""]
    if s0:
        lines.append(
            f"Gate: every categorical outcome identical, and every compared number within "
            f"|new − stored| ≤ tol × max(1, |stored|), tol = {TOL_SAME:g} for SAME and the GTEx "
            f"centroids, {TOL_SOUTH6:g} for SOUTH6. Checks marked diagnostic are not in the "
            "pre-registered list and do not gate.")
        lines.append("")
        lines.append(
            f"**Result: {'PASS' if s0['passed'] else 'FAIL'}.** Gating checks passed: "
            f"{s0['n_gating'] - s0['n_gating_failed']} of {s0['n_gating']}. "
            f"Diagnostic checks failed: {s0['n_diagnostic_failed']}.")
        lines.append("")
        t = _tbl(COORD_DIR / "stage0_reproduction.csv",
                 ["block", "item", "kind", "n", "max_abs_diff", "max_scaled_diff", "tol", "n_mismatch",
                  "pass", "gating"])
        if t:
            lines += [t, ""]
    if st.get("stage1"):
        lines += ["## Stage 1 — mixture control in every space", "",
                  "SAME positive control P, unchanged, in each space. A space is QUALIFIED only if frac rises "
                  "monotonically with f, the frac CI at f = 0.50 excludes 0 and contains its own point estimate, "
                  "and delta at f = 0.50 < 0. The f = 0 row is the noise floor and is not subtracted.", ""]
        t = _tbl(COORD_DIR / "stage1_qualification.csv")
        if t:
            lines += [t, ""]
        t = _tbl(COORD_DIR / "stage1_posctrl.csv",
                 ["space", "f", "n_total", "frac", "frac_ci_lo", "frac_ci_hi", "delta", "delta_ci_lo",
                  "delta_ci_hi", "cos", "frac_ci_covers_point", "delta_ci_covers_point", "noise_floor"])
        if t:
            lines += [t, ""]
    if st.get("stage2"):
        lines += ["## Stage 2 — SAME in every space", "",
                  "Raw distances are in different units in each space and are not compared across spaces. "
                  "Cross-space statements use the sign of delta_S, the key, frac_S and R_S = final / start "
                  "distance.", ""]
        p = COORD_DIR / "stage2_same_summary.csv"
        df = pd.read_csv(p)
        cols = ["space", "key", "sign", "flags", "ci_invalid", "n_cells_S", "delta_S", "R_S", "frac_S",
                "frac_S_ci_lo", "frac_S_ci_hi", "p_N1", "frac_rev", "asymmetry", "asymmetry_ci_lo",
                "asymmetry_ci_hi", "pluri_delta", "cos_v_u_toward"]
        q = df[df["qualified"]]
        dq = df[~df["qualified"]]
        lines += ["### Qualified spaces", "", md_table(q[cols]) if len(q) else "None.", ""]
        for _, r in q.iterrows():
            sh, genes = np.nan, ""
            for rr in pd.read_csv(COORD_DIR / "stage4_share10.csv").itertuples() if (COORD_DIR / "stage4_share10.csv").exists() else []:
                if rr.space == r["space"] and str(rr.vector).startswith("final distance"):
                    sh, genes = rr.share10, rr.top10_genes
            lines.append(
                f"- Space {r['space']}: GM00731 PartialReprog (n_cells = {int(r['n_cells_S'])}) "
                f"{'ends farther from' if r['sign'] == 'HOLDS' else 'ends closer to'} the young donor than it "
                f"started (delta_S = {r['delta_S']:+.4g} in space {r['space']} units, R_S = {r['R_S']:.3f}); "
                f"the top 10 genes carry {sh:.1%} of the squared final distance ({genes}). "
                f"Key `{r['key']}` ({r['why']}).")
        lines.append("")
        lines += ["### DISQUALIFIED spaces (no reading is drawn from these numbers)", "",
                  md_table(dq[cols]) if len(dq) else "None.", ""]
        t = _tbl(COORD_DIR / "stage2_same_context.csv")
        if t:
            lines += ["Context, reported only (frozen-ruler score in space C only):", "", t, ""]
    if st.get("stage3"):
        lines += ["## Stage 3 — SOUTH6 S1 O1→Y1 in every space", "",
                  "identity_loss and guides_disagree are carried from SOUTH6 space C. Floors are recomputed "
                  "in each space. Raw floors and deltas are in each space's units; the unit-free comparison "
                  "is median floor / ||O1 − Y1||.", ""]
        df = pd.read_csv(COORD_DIR / "stage3_south6_summary.csv")
        cols = ["space", "qualified", "outcome", "gap", "median_d_norm", "floor_min", "floor_median",
                "floor_max", "pooled_floor", "median_floor_over_gap", "nt_guide_floor", "split_floor",
                "mda_primary", "mda_stricter", "n_delta_negative", "n_below_floor", "n_q", "n_toward_young",
                "n_tie_at_min_p", "n_negative_coords_A", "n_factors_with_negative_coords_A"]
        q = df[df["qualified"]]
        dq = df[~df["qualified"]]
        lines += ["### Qualified spaces", "", md_table(q[cols]) if len(q) else "None.", ""]
        t = _tbl(COORD_DIR / "stage3_toward_young.csv")
        lines += ["toward_young factors, every space:", "", t if t else "None.", ""]
        lines += ["### DISQUALIFIED spaces (no conclusion is drawn from these numbers)", "",
                  md_table(dq[cols]) if len(dq) else "None.", ""]
        t = _tbl(COORD_DIR / "stage3_rank_agreement.csv")
        if t:
            lines += ["Rank agreement, reported only (all 1,836 factors, S1 O1→Y1):", "", t, ""]
    if st.get("stage4"):
        lines += ["## Stage 4 — diagnostic: share of squared distance in the top 10 genes", "",
                  "share10 = (sum of the 10 largest squared coordinates) / (squared length). Report only, "
                  "never a gate.", ""]
        t = _tbl(COORD_DIR / "stage4_share10.csv")
        if t:
            lines += [t, ""]
    if st.get("stage2"):
        lines += [
            "## Limitations",
            "",
            "- SAME rests on two donors. No choice of coordinates changes that.",
            "- B and C share the prior count and A has none, so A and B differ in the log and the prior together.",
            "- TMM factors are shared across spaces and are themselves estimated from log-ratios. This task does not test TMM.",
            "- In A, the screen's displacements are Hs27 CPM differences added to a GTEx CPM origin; negative "
            "coordinates are kept and counted above.",
            "- identity_loss and guides_disagree are not recomputed in A or B.",
            "- Passing P shows that a space can see a known approach between these two donors' day-0 cells, not "
            "across the bulk-to-single-cell gap where the screen's GTEx origin and target sit. Stage 2C is the "
            "screen's own check on that.",
            "- SOUTH6 Stage 1 count sums are read from results/south6/stage1.npz and results/south6/matched_* as "
            "SOUTH6 wrote them; Stage 0 shows they reproduce SOUTH6's space C vectors.",
            "",
        ]
    fails = man.get("failures") or []
    lines += ["## STOP / failures", ""]
    if fails:
        for f in fails:
            lines.append(f"- **{f['step']}:** {f['message']}")
    else:
        lines.append("None recorded.")
    lines += ["", "## Files", "", "- `src/coordspace_run.py`",
              "- `results/coordspace/PREREG.flag` (not rewritten)"]
    for f in sorted(COORD_DIR.iterdir()):
        if f.name != "PREREG.flag":
            lines.append(f"- `results/coordspace/{f.name}`")
    lines += ["- `FINDINGS_COORDSPACE.md`, `PROGRESS_COORDSPACE.md`", ""]
    FINDINGS.write_text("\n".join(lines) + "\n", encoding="utf-8")
    prog = [
        "# PROGRESS_COORDSPACE",
        "",
        f"**STOP status:** {status}",
        "",
        f"**Next action:** {next_action}",
        "",
        "**Terminal command:**",
        "",
        "```",
        "python src/coordspace_run.py",
        "```",
        "",
    ]
    if man.get("runtime_s") is not None:
        prog += [f"**Runtime of this run:** {man['runtime_s']:.0f} s", ""]
    prog += [f"PREREG.flag exists={PREREG.exists()} and was not rewritten by the run.",
             f"Failures recorded: {len(fails)}.", ""]
    PROGRESS.write_text("\n".join(prog) + "\n", encoding="utf-8")


# ----------------------------------------------------------------------------- main
def _pause(man, prereg, t0, stage, log):
    man["status"] = f"PAUSED_AFTER_STAGE{stage}"
    man["runtime_s"] = time.perf_counter() - t0
    save_manifest(man)
    write_reports(man, prereg,
                  f"Stopped after Stage {stage} on --through-stage {stage}. Not a STOP key.",
                  "Run again without --through-stage. Every run starts at Stage 0.")
    log(f"[run] paused after Stage {stage}")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--through-stage", type=int, choices=[0, 1, 2, 3, 4], default=4,
                    help="Last stage to run. Every run starts at Stage 0.")
    ap.add_argument("--workers", type=int, default=int(os.environ.get("COORDSPACE_WORKERS", "5")))
    args = ap.parse_args()
    prereg = load_prereg()
    log = Logger(LOG_PATH)
    man = new_manifest()
    save_manifest(man)
    t0 = time.perf_counter()
    try:
        log(f"[prereg] {PREREG} bytes={PREREG.stat().st_size}")
        log(f"[run] through-stage={args.through_stage} workers={args.workers}")
        ctx = stage0(man, args.workers, log)
        if args.through_stage == 0:
            return _pause(man, prereg, t0, 0, log)
        ctx = stage1(man, ctx, prereg, log)
        if args.through_stage == 1:
            return _pause(man, prereg, t0, 1, log)
        ctx["centroids"] = {sp: gtex_centroids(sp, log) for sp in SPACES}
        symbols = np.array([str(s) for s in np.asarray(load_frozen_ruler()["symbol"])], dtype=object)
        ctx = stage2(man, ctx, symbols, log)
        if args.through_stage == 2:
            return _pause(man, prereg, t0, 2, log)
        ctx = stage3(man, ctx, args.workers, log)
        if args.through_stage == 3:
            return _pause(man, prereg, t0, 3, log)
        ctx = stage4(man, ctx, log)
        dec = decide(man, ctx, log)
        man["status"] = "DONE_REVISION_TRIGGERED" if 3 in dec["rows_fired"] else "DONE"
        man["runtime_s"] = time.perf_counter() - t0
        save_manifest(man)
        dump_json(COORD_DIR / "summary.json", jsonable(dict(
            status=man["status"], decision=dec, stages=man["stages"])))
        write_reports(
            man, prereg,
            "DONE. " + ("Decision-table row 3 fired: REVISION_TRIGGERED. Report to the investigator before "
                        "anything else is run." if 3 in dec["rows_fired"]
                        else f"Decision-table rows fired: {dec['rows_fired']}."),
            "Report to the investigator." if 3 in dec["rows_fired"] else "Nothing further is required by this pre-registration.",
        )
        log("[run] DONE")
    except Halt as e:
        man["status"] = "STOP"
        man["failures"].append(dict(step=e.step, message=e.message, details=jsonable(e.details)))
        man["runtime_s"] = time.perf_counter() - t0
        save_manifest(man)
        log(f"[run] STOP [{e.step}] {e.message}")
        write_reports(man, prereg, f"STOP `{e.step}`: {e.message}",
                      "None. Never compute after a STOP and never resume a stopped run.")
        raise SystemExit(1)
    finally:
        log.close()


if __name__ == "__main__":
    main()
