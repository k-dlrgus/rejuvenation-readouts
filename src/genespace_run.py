"""GENESPACE — do the SAME and TOWARD displacement tests give the same answer
inside smaller gene spaces (MD genes; published age-up ∪ age-down genes)?

Pre-registration: results/genespace/PREREG.flag (incl. Amendment 1). This script
does not rewrite it. The full frozen-ruler space is rerun first as a
reproduction check; if it misses, nothing in MD or AGE is computed.

Imports pure functions from src/same_run.py and src/toward_run.py; calls none of
their writers. Reference numbers come only from results/same/ and
results/toward/. paper/ and paper_package/ are not read.

Usage: python src/genespace_run.py
"""
from __future__ import annotations

import hashlib
import json
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
from md4_common import STATE_NAMES  # noqa: E402
from target_common import md_table, unit  # noqa: E402
from toward_run import (  # noqa: E402
    tmm_logcpm_quiet,
    zscore_frozen,
    percentile_ci,
    sum_rows,
    load_donor,
    load_gtex_z,
    donor_panel,
    cosine_align,
    null_p,
    meets_toward,
    reading_for_donor,
    _fmt,
    _fmt_p,
    _fmt_ci,
    YOUNG_BINS,
    OLD_BINS,
    MIDDLE_BINS,
    DEST_STATES,
    MIN_CELLS,
    N_FOLDS,
    N_RANDOM,
    C3_COS_BAR,
    C3_MIN_YOUNG_TRAIN,
    C3_MIN_OLD_TRAIN,
    C3_MIN_OLD_TEST,
    DONOR_BOOT_OFFSET,
    P_BAR,
    ORIGIN_STATE,
)
from same_run import (  # noqa: E402
    cosine_full,
    axis_stats,
    ci_excludes_zero,
    ci_covers_point,
    check_gene_space,
    split_half,
    d0_fibroblast_idx,
    state_idx,
    mix_indices,
    mix_sum,
    resample_idx,
    panel_z,
    n1_cos_p,
    monotonic_increasing,
    F_VALUES,
    DEST_STATE,
    PLURI_STATE,
    NONREPROG_STATE,
)


GS_DIR = RESULTS / "genespace"
PREREG_FLAG = GS_DIR / "PREREG.flag"
FINDINGS_PATH = ROOT / "FINDINGS_GENESPACE.md"
PROGRESS_PATH = ROOT / "PROGRESS_GENESPACE.md"
SAME_DIR = RESULTS / "same"
TOWARD_DIR = RESULTS / "toward"
MD3_DIR = RESULTS / "md3"
GENESETS_PATH = MD3_DIR / "genesets.json"
SRC = ROOT / "src"

SEED = 20260914
BOOT_SEED = 20260918
N_PERM = 200
N_BOOT = 200
N_BASELINE = 200
N_BINS = 24
REPRO_TOL = 1e-6
IDENTITY_TOL = 1e-8
SPACES = ("FULL", "MD", "AGE")
LISTS = ("MD", "age_up", "age_down")
DONORS = (AGED_LINE, YOUNG_LINE)

EXPECTED_LISTS = {
    "MD": dict(n_entries=205, n_unique=205, n_mapped=205, n_unmapped=0, dups=[]),
    "age_up": dict(n_entries=1533, n_unique=1532, n_mapped=1312, n_unmapped=220, dups=["GOLGA8M"]),
    "age_down": dict(n_entries=2007, n_unique=2006, n_mapped=1777, n_unmapped=229, dups=["TMEM191A"]),
}
EXPECTED_SPACE_COLS = {"FULL": 23485, "MD": 205, "AGE": 3089}
EXPECTED_AGE_UNIQUE = 3538
EXPECTED_MD_AGE = dict(total=43, age_up=18, age_down=25)
EXPECTED_MULTI_COL_SYMBOLS = 13

MAIN_NAMES = ("O", "Y", "S", "S_rev", "Pluri", "NonReprog")
# test -> (origin row, target row, S row) in the main panel
MAIN_TESTS = {
    "forward": (0, 1, 2),
    "reverse": (1, 0, 3),
    "Pluri": (0, 1, 4),
    "NonReprog": (0, 1, 5),
}
MAIN_TEST_META = {
    "forward": (AGED_LINE, "PartialReprog", f"{AGED_LINE} d0 half A (O)", f"{YOUNG_LINE} d0 half A (Y)"),
    "reverse": (YOUNG_LINE, "PartialReprog", f"{YOUNG_LINE} d0 half A (Y)", f"{AGED_LINE} d0 half A (O)"),
    "Pluri": (AGED_LINE, "Pluripotency", f"{AGED_LINE} d0 half A (O)", f"{YOUNG_LINE} d0 half A (Y)"),
    "NonReprog": (AGED_LINE, "NonReprog", f"{AGED_LINE} d0 half A (O)", f"{YOUNG_LINE} d0 half A (Y)"),
}
SAME_STATS = ("cos", "progress", "delta", "perp", "rel_delta")
GTEX_STATS = ("cos_S", "progress", "perp", "delta_young", "rel_delta_young", "delta_old", "rel_delta_old")
BASELINES = ("sec6_mu", "A1_mu_zero")

GSE325735_BAN = "GSE325735 is out of scope and is not opened."


# --------------------------------------------------------------------------- integrity
def protected_paths():
    paths = []
    for d in (SAME_DIR, TOWARD_DIR, MD3_DIR):
        paths += sorted(p for p in d.rglob("*") if p.is_file())
    paths += [
        FROZEN_RULER, SRC / "same_run.py", SRC / "toward_run.py",
        ROOT / "FINDINGS_SAME.md", ROOT / "FINDINGS_TOWARD.md",
    ]
    return paths


def sha256_file(path):
    h = hashlib.sha256()
    with open(path, "rb") as f:
        for chunk in iter(lambda: f.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def hash_protected():
    return {str(Path(p).relative_to(ROOT)).replace("\\", "/"): sha256_file(p) for p in protected_paths()}


# --------------------------------------------------------------------------- statistics
def geometry_label(delta, progress, perp):
    if not (np.isfinite(delta) and np.isfinite(progress) and np.isfinite(perp)):
        return "NA"
    if delta < 0:
        return "APPROACH"
    if delta > 0 and progress > 1 and (progress - 1) ** 2 > perp ** 2:
        return "OVERSHOOT"
    if delta > 0:
        return "SIDEWAYS"
    return "TIE"


def same_axis(z_origin, z_target, z_S, tag):
    """Section-4 statistics, same-platform frame. STOP on identity mismatch."""
    st = axis_stats(z_origin, z_target, z_S)
    vn = st["v_norm"]
    if (not np.isfinite(vn)) or vn < 1e-12:
        raise StopStep("v_norm", f"{tag}: ||v||={vn} < 1e-12. STOP for this space.")
    prog = st["frac"]
    perp = float(np.linalg.norm(st["d"] - prog * st["v"]) / vn)
    start = st["dist_origin_to_target"]
    final = st["dist_to_target"]
    lhs = final ** 2 / vn ** 2
    rhs = (1.0 - prog) ** 2 + perp ** 2
    mism = abs(lhs - rhs) / max(abs(lhs), 1e-300)
    if not (mism <= IDENTITY_TOL):
        raise StopStep(
            "identity",
            f"{tag}: final²/||v||²={lhs!r} vs (1-progress)²+perp²={rhs!r}; "
            f"relative mismatch {mism:.3e} > {IDENTITY_TOL:g}. Code bug. STOP.",
        )
    return dict(
        cos=st["cos"], progress=prog, delta=st["delta"], perp=perp,
        rel_delta=float(st["delta"] / start) if start > 0 else np.nan,
        dist_final=final, dist_start=start, v_norm=vn, d_norm=st["d_norm"],
        identity_mismatch=float(mism),
        label=geometry_label(st["delta"], prog, perp),
        v=st["v"], d=st["d"],
    )


def gtex_axis(z0, zS, c_y, c_o):
    """Section-4 statistics against w(X) = c_young(X) - c_old(X). Origin = donor d0."""
    w = np.asarray(c_y, float) - np.asarray(c_o, float)
    u, wn = unit(w)
    d = np.asarray(zS, float) - np.asarray(z0, float)
    if wn < 1e-12:
        nan = np.nan
        return dict(cos_S=nan, progress=nan, perp=nan, d=d, u=u, w_norm=wn)
    prog = float(np.dot(d, w) / wn ** 2)
    perp = float(np.linalg.norm(d - prog * w) / wn)
    dy = float(np.linalg.norm(zS - c_y))
    dy0 = float(np.linalg.norm(z0 - c_y))
    do = float(np.linalg.norm(zS - c_o))
    do0 = float(np.linalg.norm(z0 - c_o))
    return dict(
        cos_S=cosine_align(d, u), progress=prog, perp=perp,
        dist_young=dy, dist_young_origin=dy0, delta_young=dy - dy0,
        rel_delta_young=(dy - dy0) / dy0 if dy0 > 0 else np.nan,
        dist_old=do, dist_old_origin=do0, delta_old=do - do0,
        rel_delta_old=(do - do0) / do0 if do0 > 0 else np.nan,
        d_norm=float(np.linalg.norm(d)), w_norm=float(wn), d=d, u=u,
    )


def ci_with_flag(vals, point):
    lo, hi = percentile_ci(vals)
    return lo, hi, bool(ci_covers_point(point, lo, hi))


def valid_str(v):
    return "valid" if v else "invalid"


# --------------------------------------------------------------------------- gene spaces
def map_gene_spaces(frozen, genesets, log):
    ruler_sym_u = np.array([str(s).upper() for s in np.asarray(frozen["symbol"])], dtype=object)
    sym_to_cols = {}
    for j, s in enumerate(ruler_sym_u):
        sym_to_cols.setdefault(s, []).append(j)
    multi = {s: c for s, c in sym_to_cols.items() if len(c) > 1}
    lists = {}
    entry_rows = []
    for key in LISTS:
        raw = list(genesets[key])
        seen = {}
        uniq, dups = [], []
        for pos, e in enumerate(raw):
            n = str(e).strip().upper()
            if n in seen:
                dups.append(n)
                entry_rows.append(dict(record="duplicate_entry", list=key, list_pos=pos,
                                       entry=str(e), entry_norm=n, first_pos=seen[n]))
                continue
            seen[n] = pos
            uniq.append(n)
        mapped = [n for n in uniq if n in sym_to_cols]
        unmapped = [n for n in uniq if n not in sym_to_cols]
        for n in unmapped:
            entry_rows.append(dict(record="unmapped_entry", list=key, list_pos=seen[n],
                                   entry=str(raw[seen[n]]), entry_norm=n, first_pos=seen[n]))
        cols = sorted({c for n in mapped for c in sym_to_cols[n]})
        lists[key] = dict(
            n_entries=len(raw), n_unique=len(uniq), n_mapped=len(mapped),
            n_unmapped=len(unmapped), dups=dups, cols=np.asarray(cols, int),
            n_cols=len(cols), unmapped=unmapped, mapped=set(mapped),
            multi_in_list=sorted(set(mapped) & set(multi)),
        )
        log(f"[map] {key}: entries={len(raw)} unique={len(uniq)} mapped={len(mapped)} "
            f"unmapped={len(unmapped)} dups={dups} cols={len(cols)}")

    up_set, down_set = set(lists["age_up"]["mapped"]), set(lists["age_down"]["mapped"])
    up_list = {str(e).strip().upper() for e in genesets["age_up"]}
    down_list = {str(e).strip().upper() for e in genesets["age_down"]}
    md_cols = set(lists["MD"]["cols"].tolist())
    up_cols = set(lists["age_up"]["cols"].tolist())
    down_cols = set(lists["age_down"]["cols"].tolist())
    age_cols = sorted(up_cols | down_cols)
    n_ruler = len(ruler_sym_u)
    spaces = {
        "FULL": np.arange(n_ruler, dtype=int),
        "MD": np.asarray(sorted(md_cols), int),
        "AGE": np.asarray(age_cols, int),
    }
    checks = dict(
        up_down_overlap_lists=len(up_list & down_list),
        up_down_overlap_mapped=len(up_set & down_set),
        age_unique=len(up_list | down_list),
        md_age_total=len(md_cols & set(age_cols)),
        md_age_up=len(md_cols & up_cols),
        md_age_down=len(md_cols & down_cols),
        multi_col_symbols=len(multi),
        multi_col_in_lists=sorted(set().union(*[set(lists[k]["multi_in_list"]) for k in LISTS])),
        space_cols={k: int(v.size) for k, v in spaces.items()},
    )
    log(f"[map] checks {checks}")

    mismatches = []
    for key, exp in EXPECTED_LISTS.items():
        got = lists[key]
        for f in ("n_entries", "n_unique", "n_mapped", "n_unmapped"):
            if got[f] != exp[f]:
                mismatches.append(f"{key}.{f}: got {got[f]} expected {exp[f]}")
        if sorted(got["dups"]) != sorted(exp["dups"]):
            mismatches.append(f"{key}.dups: got {got['dups']} expected {exp['dups']}")
        if got["n_cols"] != got["n_mapped"]:
            mismatches.append(f"{key}: n_cols {got['n_cols']} != n_mapped {got['n_mapped']}")
    for k, n in EXPECTED_SPACE_COLS.items():
        if checks["space_cols"][k] != n:
            mismatches.append(f"space {k}: got {checks['space_cols'][k]} cols expected {n}")
    if checks["up_down_overlap_lists"] != 0 or checks["up_down_overlap_mapped"] != 0:
        mismatches.append(f"age_up ∩ age_down: lists={checks['up_down_overlap_lists']} "
                          f"mapped={checks['up_down_overlap_mapped']} expected 0/0")
    if checks["age_unique"] != EXPECTED_AGE_UNIQUE:
        mismatches.append(f"AGE unique entries: got {checks['age_unique']} expected {EXPECTED_AGE_UNIQUE}")
    for f, exp in (("md_age_total", "total"), ("md_age_up", "age_up"), ("md_age_down", "age_down")):
        if checks[f] != EXPECTED_MD_AGE[exp]:
            mismatches.append(f"MD ∩ AGE {exp}: got {checks[f]} expected {EXPECTED_MD_AGE[exp]}")
    if checks["multi_col_symbols"] != EXPECTED_MULTI_COL_SYMBOLS:
        mismatches.append(f"ruler symbols spanning >1 column: got {checks['multi_col_symbols']} "
                          f"expected {EXPECTED_MULTI_COL_SYMBOLS}")
    if checks["multi_col_in_lists"]:
        mismatches.append(f"multi-column symbols in lists: {checks['multi_col_in_lists']} (expected none)")
    if mismatches:
        raise StopStep("gene_counts", "mapping counts differ from PREREG.flag: " + "; ".join(mismatches),
                       dict(checks=checks, mismatches=mismatches))
    membership = dict(
        MD=np.isin(np.arange(n_ruler), lists["MD"]["cols"]),
        age_up=np.isin(np.arange(n_ruler), lists["age_up"]["cols"]),
        age_down=np.isin(np.arange(n_ruler), lists["age_down"]["cols"]),
    )
    return dict(lists=lists, spaces=spaces, checks=checks, entry_rows=entry_rows,
                membership=membership, ruler_sym_u=ruler_sym_u)


def mu_bins(mu):
    mu = np.asarray(mu, float)
    edges = np.quantile(mu, np.linspace(0.0, 1.0, N_BINS + 1))
    b = np.clip(np.searchsorted(edges, mu, side="right") - 1, 0, N_BINS - 1)
    return b.astype(int), edges


def strata_counts(cols, bins, zero):
    cols = np.asarray(cols, int)
    out = {}
    for b in range(N_BINS):
        in_b = bins[cols] == b
        out[f"bin{b:02d}_zero"] = int((in_b & zero[cols]).sum())
        out[f"bin{b:02d}_nonzero"] = int((in_b & ~zero[cols]).sum())
    return out


def bin_counts(cols, bins):
    cols = np.asarray(cols, int)
    return {f"bin{b:02d}": int((bins[cols] == b).sum()) for b in range(N_BINS)}


def draw_baseline_sets(cols, bins, zero, mode):
    """One Generator(SEED) per call. Strata in order; space genes not excluded."""
    rng = np.random.default_rng(SEED)
    cols = np.asarray(cols, int)
    if mode == "sec6_mu":
        strata = [(bins == b) for b in range(N_BINS)]
    else:
        strata = []
        for b in range(N_BINS):
            strata.append((bins == b) & zero)
            strata.append((bins == b) & ~zero)
    pools = [np.flatnonzero(m) for m in strata]
    need = [int(np.isin(cols, p).sum()) for p in pools]
    sets = np.empty((N_BASELINE, cols.size), dtype=int)
    for i in range(N_BASELINE):
        parts = []
        for pool, n in zip(pools, need):
            if n > 0:
                parts.append(rng.choice(pool, size=n, replace=False))
        sets[i] = np.concatenate(parts)
    return sets, need


# --------------------------------------------------------------------------- per-space point blocks
def same_point_block(Zmain, Zmix, cols, mixes, sp):
    zm = Zmain[:, cols]
    zx = Zmix[:, cols]
    main = {}
    for test, (o, t, s) in MAIN_TESTS.items():
        main[test] = same_axis(zm[o], zm[t], zm[s], f"{sp} {test}")
    mix = [same_axis(zx[0], zx[1], zx[2 + i], f"{sp} mix f={m['f']:.2f}") for i, m in enumerate(mixes)]
    return main, mix


def gtex_point_block(donor_z, anchors_full, cols):
    c_y = anchors_full["c_young"][cols]
    c_o = anchors_full["c_old"][cols]
    out = {}
    for line in DONORS:
        dz = donor_z[line]
        z0 = dz["Z"][0, cols]
        out[line] = {}
        for j, nm in enumerate(dz["names"]):
            if nm == ORIGIN_STATE:
                continue
            out[line][nm] = gtex_axis(z0, dz["Z"][j, cols], c_y, c_o)
    return out


def flat_point_stats(cols, ctx, include_gtex, sp):
    main, mix = same_point_block(ctx["Zmain"], ctx["Zmix"], cols, ctx["mixes"], sp)
    out = {}
    for test, st in main.items():
        for stat in SAME_STATS:
            out[("same", "", test, stat)] = st[stat]
    out[("same", "", "forward", "geometry_label")] = main["forward"]["label"]
    for m, st in zip(ctx["mixes"], mix):
        out[("mixture", "", f"f={m['f']:.2f}", "progress")] = st["progress"]
        out[("mixture", "", f"f={m['f']:.2f}", "delta")] = st["delta"]
    if include_gtex:
        g = gtex_point_block(ctx["donor_z"], ctx["anchors"], cols)
        for line in DONORS:
            for nm, rec in g[line].items():
                for stat in GTEX_STATS:
                    out[("gtex", line, nm, stat)] = rec[stat]
    return out


# --------------------------------------------------------------------------- C3
def c3_prepare(gtex, log):
    Z, age = gtex["Z"], gtex["age"]
    n = int(Z.shape[0])
    rng = np.random.default_rng(SEED)
    order = rng.permutation(n)
    folds = np.empty(n, dtype=int)
    folds[order] = np.arange(n) % N_FOLDS
    out = []
    for k in range(N_FOLDS):
        tr = folds != k
        te = folds == k
        age_tr, age_te = age[tr], age[te]
        young_tr = np.isin(age_tr, YOUNG_BINS)
        old_tr = np.isin(age_tr, OLD_BINS)
        old_te = np.isin(age_te, OLD_BINS)
        rec = dict(fold=k, n_train=int(tr.sum()), n_test=int(te.sum()),
                   n_young_train=int(young_tr.sum()), n_old_train=int(old_tr.sum()),
                   n_old_test=int(old_te.sum()), construction_ok=True, reason="")
        if (rec["n_young_train"] < C3_MIN_YOUNG_TRAIN or rec["n_old_train"] < C3_MIN_OLD_TRAIN
                or rec["n_old_test"] < C3_MIN_OLD_TEST):
            rec["construction_ok"] = False
            rec["reason"] = (f"construction minimum failed: n_young_train={rec['n_young_train']} "
                             f"n_old_train={rec['n_old_train']} n_old_test={rec['n_old_test']}")
            out.append(rec)
            continue
        Ztr, Zte = Z[tr], Z[te]
        rec["c_y"] = Ztr[young_tr].mean(0)
        rec["c_o"] = Ztr[old_tr].mean(0)
        rec["z_old_te"] = Zte[old_te]
        out.append(rec)
    log(f"[C3 prep] folds={[(r['fold'], r['n_young_train'], r['n_old_train'], r['n_old_test']) for r in out]}")
    return out


def c3_space(folds, cols, sp, log):
    rng_boot = np.random.default_rng(BOOT_SEED)
    rows, boots = [], {}
    construction_ok = True
    for fr in folds:
        rec = {k: fr[k] for k in ("fold", "n_train", "n_test", "n_young_train", "n_old_train",
                                  "n_old_test", "construction_ok", "reason")}
        rec.update(space=sp, cos=np.nan, d_norm=np.nan, cos_ci_lo=np.nan, cos_ci_hi=np.nan,
                   cos_ci_valid=False, pass_bar=False)
        if not fr["construction_ok"]:
            construction_ok = False
            rows.append(rec)
            continue
        c_y = fr["c_y"][cols]
        u_tr, nrm = unit(c_y - fr["c_o"][cols])
        if nrm < 1e-12:
            rec["construction_ok"] = False
            rec["reason"] = "train ||c_young-c_old|| degenerate"
            construction_ok = False
            rows.append(rec)
            continue
        zte = fr["z_old_te"][:, cols]
        d = c_y - zte.mean(0)
        rec["d_norm"] = float(np.linalg.norm(d))
        rec["cos"] = cosine_align(d, u_tr)
        rec["axis_norm_train"] = float(nrm)
        n_ote = int(zte.shape[0])
        b = []
        for _ in range(N_BOOT):
            take = rng_boot.integers(0, n_ote, size=n_ote)
            b.append(cosine_align(c_y - zte[take].mean(0), u_tr))
        boots[fr["fold"]] = np.asarray(b, float)
        rec["cos_ci_lo"], rec["cos_ci_hi"], rec["cos_ci_valid"] = ci_with_flag(b, rec["cos"])
        rec["pass_bar"] = bool(np.isfinite(rec["cos"]) and rec["cos"] > C3_COS_BAR)
        rows.append(rec)
        log(f"[C3 {sp}] fold {fr['fold']} cos={_fmt(rec['cos'], 4)} "
            f"CI={_fmt_ci(rec['cos_ci_lo'], rec['cos_ci_hi'])} {valid_str(rec['cos_ci_valid'])} "
            f"pass_bar={rec['pass_bar']}")
    cosv = np.array([r["cos"] for r in rows], float)
    passed = bool(len(rows) == N_FOLDS and construction_ok and np.isfinite(cosv).all()
                  and (cosv > C3_COS_BAR).all())
    log(f"[C3 {sp}] pass={passed} min_cos={_fmt(np.nanmin(cosv) if np.isfinite(cosv).any() else np.nan)}")
    return rows, boots, passed


# --------------------------------------------------------------------------- reading
def mixture_gate(mix_rows):
    prog = [float(r["progress"]) for r in mix_rows]
    r50 = next(r for r in mix_rows if abs(r["f"] - 0.50) < 1e-12)
    mono = monotonic_increasing(prog)
    lo, hi = r50["progress_ci_lo"], r50["progress_ci_hi"]
    excl = ci_excludes_zero(lo, hi)
    valid = bool(r50["progress_ci_valid"])
    d50 = r50["delta"]
    dneg = bool(np.isfinite(d50) and d50 < 0)
    failed = []
    if not mono:
        failed.append("progress not non-decreasing in f (f=0, 0.10, 0.25, 0.50: "
                      + ", ".join(_fmt(p) for p in prog) + ")")
    if not excl:
        failed.append(f"progress CI at f=0.50 {_fmt_ci(lo, hi)} does not exclude 0")
    if not valid:
        failed.append(f"progress CI at f=0.50 {_fmt_ci(lo, hi)} is invalid "
                      f"(does not contain the point {_fmt(r50['progress'])})")
    if not dneg:
        failed.append(f"delta at f=0.50 = {_fmt(d50)} is not < 0")
    return dict(
        passed=bool(mono and excl and valid and dneg), monotonic=mono,
        progress50_ci_excludes_0=excl, progress50_ci_valid=valid, delta50_lt_0=dneg,
        progress=prog, progress50=r50["progress"], progress50_ci=[lo, hi], delta50=d50,
        failed=failed,
    )


def fire_space_reading(sp, n_genes, n_S, fwd, rev, gate):
    if not gate["passed"]:
        key = "NO READING"
        text = (f"In the {sp} gene space ({n_genes} ruler genes), the mixture positive control failed "
                f"({'; '.join(gate['failed'])}); the statistic cannot detect a known approach in this "
                "space, and its forward numbers are not interpreted.")
    else:
        lab = fwd["label"]
        a = (fwd["progress"] - 1.0) ** 2
        b = fwd["perp"] ** 2
        if lab == "APPROACH":
            key = "APPROACH"
            text = (f"In the {sp} gene space ({n_genes} ruler genes), aged partially-reprogrammed cells "
                    f"(GM00731 PartialReprog, n_cells={n_S}) end closer to the young donor's day-0 cells "
                    f"than they started: change in distance {_fmt(fwd['delta'])} ({_fmt(fwd['rel_delta'])} "
                    f"of the starting gap), progress {_fmt(fwd['progress'])}, sideways {_fmt(fwd['perp'])}. "
                    f"The reverse test in the same space gives progress {_fmt(rev['progress'])} and change "
                    f"in distance {_fmt(rev['delta'])}; approach alone does not separate rejuvenation from "
                    "donor convergence. One young donor, so 'young' is confounded with that donor.")
        elif lab == "OVERSHOOT":
            key = "OVERSHOOT"
            text = (f"In the {sp} gene space ({n_genes} ruler genes), aged partially-reprogrammed cells "
                    f"(GM00731 PartialReprog, n_cells={n_S}) move past the young donor's day-0 cells along "
                    f"the aged-to-young axis (progress {_fmt(fwd['progress'])} > 1) and end farther from them "
                    f"than they started (change in distance {_fmt(fwd['delta'])}); the overshoot term "
                    f"(progress − 1)² = {_fmt(a)} is larger than the sideways term perp² = {_fmt(b)}. "
                    "This is not approach to the young donor.")
        elif lab == "SIDEWAYS":
            key = "SIDEWAYS"
            text = (f"In the {sp} gene space ({n_genes} ruler genes), aged partially-reprogrammed cells "
                    f"(GM00731 PartialReprog, n_cells={n_S}) end farther from the young donor's day-0 cells "
                    f"than they started (change in distance {_fmt(fwd['delta'])}), mainly through movement "
                    f"off the aged-to-young axis: sideways term perp² = {_fmt(b)}, overshoot term "
                    f"(progress − 1)² = {_fmt(a)}, progress {_fmt(fwd['progress'])}. "
                    "This is not approach to the young donor.")
        elif lab == "TIE":
            key = "NO READING (tie)"
            text = (f"In the {sp} gene space ({n_genes} ruler genes), change in distance is exactly 0; "
                    "no reading.")
        else:
            raise StopStep("reading", f"{sp}: forward label {lab!r} is not defined (non-finite statistic).")
    ci_line = (f"delta CI {_fmt_ci(fwd['delta_ci_lo'], fwd['delta_ci_hi'])} "
               f"({valid_str(fwd['delta_ci_valid'])}); progress CI "
               f"{_fmt_ci(fwd['progress_ci_lo'], fwd['progress_ci_hi'])} "
               f"({valid_str(fwd['progress_ci_valid'])}).")
    return dict(space=sp, key=key, text=text, ci_line=ci_line, gate=gate,
                n_genes=int(n_genes), n_S=int(n_S))


# --------------------------------------------------------------------------- run
class RunState:
    def __init__(self):
        self.man = dict(
            task="GENESPACE", status="INIT", failures=[],
            seeds=dict(split_mixture_perm_donor_splits_c3_folds_baseline=SEED, bootstrap=BOOT_SEED),
            n_perm=N_PERM, n_random=N_RANDOM, n_boot=N_BOOT, n_baseline=N_BASELINE,
            repro_tol=REPRO_TOL, identity_tol=IDENTITY_TOL, n_mu_bins=N_BINS,
        )
        self.reading = {}
        self.summary = {}
        self.findings_ctx = {}


def run(log, rs):
    prereg_text = PREREG_FLAG.read_text(encoding="utf-8")
    rs.prereg_text = prereg_text
    log("=" * 100)
    log(f"GENESPACE  seed={SEED}  boot={BOOT_SEED}  n_perm={N_PERM}  n_boot={N_BOOT}  "
        f"n_random={N_RANDOM}  n_baseline={N_BASELINE}")
    log("=" * 100)
    log(GSE325735_BAN)
    log("PREREG.flag read; not rewritten. paper/ and paper_package/ are not read.")
    rs.man["status"] = "RUNNING"
    rs.man["prereg_sha256"] = sha256_file(PREREG_FLAG)
    write_progress(rs, "hash protected files; map gene spaces")

    rs.man["hashes_before"] = hash_protected()
    log(f"[hash] protected files before run: n={len(rs.man['hashes_before'])}")

    # ---- section 1: gene spaces
    frozen_npz = load_frozen_ruler()
    frozen = {k: np.asarray(frozen_npz[k]) for k in ("mu", "sd", "w", "symbol", "ensembl")}
    genesets = load_json(GENESETS_PATH)
    gm = map_gene_spaces(frozen, genesets, log)
    spaces = gm["spaces"]
    rs.gm = gm
    rs.man["gene_counts"] = dict(
        lists={k: {f: gm["lists"][k][f] for f in ("n_entries", "n_unique", "n_mapped", "n_unmapped", "dups", "n_cols")}
               for k in LISTS},
        spaces={k: int(v.size) for k, v in spaces.items()},
        checks=gm["checks"],
    )

    # ---- section 2 inputs: donors, split, panels
    obs_a, Y_a = load_donor(AGED_LINE, log)
    obs_y, Y_y = load_donor(YOUNG_LINE, log)
    check_gene_space(Y_a, frozen, AGED_LINE)
    check_gene_space(Y_y, frozen, YOUNG_LINE)
    sp_file = np.load(SAME_DIR / "split_idx.npz")
    o_idx, y_idx = np.asarray(sp_file["O"], int), np.asarray(sp_file["Y"], int)
    aged_b, young_b = np.asarray(sp_file["aged_B"], int), np.asarray(sp_file["young_B"], int)
    ra, rb = split_half(d0_fibroblast_idx(obs_a), SEED)
    ya, yb = split_half(d0_fibroblast_idx(obs_y), SEED)
    split_ok = bool(np.array_equal(ra, o_idx) and np.array_equal(rb, aged_b)
                    and np.array_equal(ya, y_idx) and np.array_equal(yb, young_b))
    log(f"[split] O={o_idx.size} aged_B={aged_b.size} Y={y_idx.size} young_B={young_b.size} "
        f"rerun split_half equals split_idx.npz: {split_ok}")
    if not split_ok:
        raise StopStep("split", "split_half(d0, 20260914) does not equal results/same/split_idx.npz. STOP.")

    s_idx = state_idx(obs_a, DEST_STATE)
    srev_idx = state_idx(obs_y, DEST_STATE)
    pluri_idx = state_idx(obs_a, PLURI_STATE)
    nr_idx = state_idx(obs_a, NONREPROG_STATE)
    main_spec = [
        ("O", Y_a, o_idx), ("Y", Y_y, y_idx), ("S", Y_a, s_idx),
        ("S_rev", Y_y, srev_idx), ("Pluri", Y_a, pluri_idx), ("NonReprog", Y_a, nr_idx),
    ]
    n_cells = {nm: int(ix.size) for nm, _, ix in main_spec}
    log(f"[main panel] n_cells={n_cells}")
    main_counts = [sum_rows(Ym, ix) for _, Ym, ix in main_spec]
    zero = (np.vstack(main_counts) == 0).all(0)
    bins, edges = mu_bins(frozen["mu"])
    log(f"[A1] zero genes in main SAME panel (fixed before any statistic): {int(zero.sum())}")

    mixes, n_tot = mix_indices(aged_b, young_b, SEED)
    mix_counts = [sum_rows(Y_a, o_idx), sum_rows(Y_y, y_idx)]
    for m in mixes:
        mix_counts.append(mix_sum(Y_a, Y_y, m["aged_idx"], m["young_idx"]))
    log(f"[mixtures] n_total={n_tot} " + " ".join(
        f"f={m['f']:.2f}:{m['n_aged']}+{m['n_young']}" for m in mixes))

    write_gene_map(gm, frozen, bins, zero)
    space_info = build_space_info(gm, bins, zero, edges)
    rs.space_info = space_info
    rs.n_cells = n_cells
    rs.mixes_meta = [dict(f=m["f"], n_aged=m["n_aged"], n_young=m["n_young"], n_total=m["n_total"]) for m in mixes]

    # ---- full-gene z matrices (stats are sliced from these)
    Zmain, miss_main, _, _ = panel_z(main_counts, frozen)
    Zmix, miss_mix, _, _ = panel_z(mix_counts, frozen)
    if not np.array_equal(miss_main, zero):
        raise StopStep("zero_mask", "main-panel z=0 mask differs from the A1 zero mask. STOP.")

    gtex = load_gtex_z(log)
    age = gtex["age"]
    young_m, old_m, mid_m = np.isin(age, YOUNG_BINS), np.isin(age, OLD_BINS), np.isin(age, MIDDLE_BINS)
    anchors = dict(c_young=gtex["Z"][young_m].mean(0), c_old=gtex["Z"][old_m].mean(0))
    rs.anchor_n = dict(n_young=int(young_m.sum()), n_old=int(old_m.sum()), n_middle=int(mid_m.sum()),
                       n_donors=int(len(age)))
    log(f"[anchors] {rs.anchor_n}")

    donors = {AGED_LINE: (obs_a, Y_a), YOUNG_LINE: (obs_y, Y_y)}
    donor_z = {}
    for line in DONORS:
        obs, Ym = donors[line]
        panel = donor_panel(obs, Ym, log, line)
        Zd, miss_d = zscore_frozen(panel["logcpm"], frozen, panel["C"])
        n_by = {st: int((obs["label"].astype(str) == st).sum()) for st in STATE_NAMES}
        n_by[ORIGIN_STATE] = int(panel["n_d0"])
        donor_z[line] = dict(Z=Zd, names=list(panel["names"]), idx=panel["idx"], missing=miss_d,
                             skipped=panel["skipped"], n_by=n_by, Y=Ym)
    rs.donor_meta = {line: dict(names=donor_z[line]["names"], skipped=donor_z[line]["skipped"],
                                n_by=donor_z[line]["n_by"]) for line in DONORS}

    zero_by_panel = {}
    for sp in SPACES:
        cols = spaces[sp]
        zero_by_panel[sp] = dict(
            main_same=int(miss_main[cols].sum()), mixture=int(miss_mix[cols].sum()),
            **{f"gtex_{line}": int(donor_z[line]["missing"][cols].sum()) for line in DONORS},
        )
    rs.zero_by_panel = zero_by_panel
    log(f"[z=0 per space and panel] {zero_by_panel}")

    ctx = dict(Zmain=Zmain, Zmix=Zmix, mixes=mixes, donor_z=donor_z, anchors=anchors)

    # ---- section 8: stop rule on FULL
    write_progress(rs, "FULL reproduction check (section 8)")
    full = spaces["FULL"]
    f_main, f_mix = same_point_block(Zmain, Zmix, full, mixes, "FULL")
    f_gtex = gtex_point_block(donor_z, anchors, full)
    ref_same = load_json(SAME_DIR / "summary.json")
    ref_tow = pd.read_csv(TOWARD_DIR / "stats.csv")
    trow = ref_tow[(ref_tow.cell_line.astype(str) == AGED_LINE) & (ref_tow.state.astype(str) == "PartialReprog")]
    if len(trow) != 1:
        raise StopStep("repro", f"results/toward/stats.csv has {len(trow)} GM00731 PartialReprog rows. STOP.")
    trow = trow.iloc[0]
    pairs = [
        ("same/summary.json frac_S", f_main["forward"]["progress"], ref_same["frac_S"]),
        ("same/summary.json delta_S", f_main["forward"]["delta"], ref_same["delta_S"]),
        ("same/summary.json frac_rev", f_main["reverse"]["progress"], ref_same["frac_rev"]),
    ]
    for i, m in enumerate(mixes):
        pairs.append((f"same/summary.json p.fracs[f={m['f']:.2f}]", f_mix[i]["progress"], ref_same["p"]["fracs"][i]))
    pairs.append(("toward/stats.csv GM00731 PartialReprog cos_S",
                  f_gtex[AGED_LINE]["PartialReprog"]["cos_S"], float(trow["cos_S"])))
    pairs.append(("toward/stats.csv GM00731 PartialReprog delta_young",
                  f_gtex[AGED_LINE]["PartialReprog"]["delta_young"], float(trow["delta_young"])))
    point_rows = []
    for name, new, ref in pairs:
        diff = float(new) - float(ref) if ref is not None else np.nan
        ok = bool(np.isfinite(diff) and abs(diff) <= REPRO_TOL)
        point_rows.append(dict(quantity=name, rerun=float(new), reference=ref, diff=diff, pass_=ok))
        log(f"[repro] {name}: rerun={float(new)!r} ref={ref!r} diff={diff:.3e} {'OK' if ok else 'MISS'}")
    repro_pass = all(r["pass_"] for r in point_rows)
    rs.repro = dict(tol=REPRO_TOL, passed=repro_pass,
                    point=[{**{k: v for k, v in r.items() if k != "pass_"}, "pass": r["pass_"]} for r in point_rows])
    dump_json(GS_DIR / "repro_check.json", jsonable(rs.repro))
    if not repro_pass:
        rs.reading = dict(key="reproduction_failed",
                          text="FULL rerun missed the section-8 stop rule. No MD or AGE statistic is computed.")
        dump_json(GS_DIR / "reading.json", jsonable(rs.reading))
        raise StopStep("reproduction_failed", "FULL rerun misses the stop rule (|new - reference| > 1e-6). "
                       "No MD or AGE statistic computed.", rs.repro)
    log("[repro] stop rule PASSED; MD and AGE may be computed.")

    # ---- section 2: same-platform, per space
    write_progress(rs, "same-platform point estimates, N1, bootstrap")
    same_pt, mix_pt = {}, {}
    for sp in SPACES:
        same_pt[sp], mix_pt[sp] = same_point_block(Zmain, Zmix, spaces[sp], mixes, sp)
        f = same_pt[sp]["forward"]
        log(f"[same {sp}] n_genes={spaces[sp].size} fwd progress={_fmt(f['progress'])} delta={_fmt(f['delta'])} "
            f"perp={_fmt(f['perp'])} label={f['label']} | rev progress={_fmt(same_pt[sp]['reverse']['progress'])}")
        for m, st in zip(mixes, mix_pt[sp]):
            log(f"[mix {sp}] f={m['f']:.2f} progress={_fmt(st['progress'])} delta={_fmt(st['delta'])} "
                f"perp={_fmt(st['perp'])} cos={_fmt(st['cos'])}")

    nulls = {}
    n1 = {}
    for sp in SPACES:
        for test in ("forward", "reverse"):
            st = same_pt[sp][test]
            p, cos_obs, null = n1_cos_p(st["d"], st["v"], SEED, N_PERM)
            n1[(sp, test)] = dict(p=p, null_mean=float(np.nanmean(null)), null_p95=float(np.nanpercentile(null, 95)))
            nulls[f"null_same_N1__{sp}__{test}"] = null
            log(f"[N1 same {sp} {test}] cos={_fmt(cos_obs)} p={_fmt_p(p)}")

    boot = {}

    def bput(key, val):
        boot.setdefault(key, []).append(float(val) if val is not None else np.nan)

    rng_b = np.random.default_rng(BOOT_SEED)
    for b in range(N_BOOT):
        b_counts = [mix_counts[0], mix_counts[1]]
        for m in mixes:
            b_counts.append(mix_sum(Y_a, Y_y, resample_idx(m["aged_idx"], rng_b),
                                    resample_idx(m["young_idx"], rng_b)))
        Zb, _, _, _ = panel_z(b_counts, frozen)
        for sp in SPACES:
            zb = Zb[:, spaces[sp]]
            for i, m in enumerate(mixes):
                st = same_axis(zb[0], zb[1], zb[2 + i], f"boot{b} {sp} mix f={m['f']:.2f}")
                for stat in SAME_STATS:
                    bput(f"same_mix__{sp}__f{m['f']:.2f}__{stat}", st[stat])
        if (b + 1) % 50 == 0:
            log(f"[mix boot] {b + 1}/{N_BOOT}")

    rng_b = np.random.default_rng(BOOT_SEED)
    for b in range(N_BOOT):
        b_counts = []
        for i, (nm, Ym, ix) in enumerate(main_spec):
            if nm in ("O", "Y"):
                b_counts.append(main_counts[i])
            else:
                b_counts.append(sum_rows(Ym, resample_idx(ix, rng_b)))
        Zb, _, _, _ = panel_z(b_counts, frozen)
        for sp in SPACES:
            zb = Zb[:, spaces[sp]]
            res = {}
            for test, (o, t, s) in MAIN_TESTS.items():
                res[test] = same_axis(zb[o], zb[t], zb[s], f"boot{b} {sp} {test}")
                for stat in SAME_STATS:
                    bput(f"same_main__{sp}__{test}__{stat}", res[test][stat])
            bput(f"same_main__{sp}__asymmetry", res["forward"]["progress"] - res["reverse"]["progress"])
        if (b + 1) % 50 == 0:
            log(f"[main boot] {b + 1}/{N_BOOT}")

    same_rows, mix_rows = [], []
    for sp in SPACES:
        asym = same_pt[sp]["forward"]["progress"] - same_pt[sp]["reverse"]["progress"]
        for test in MAIN_TESTS:
            st = same_pt[sp][test]
            line, state, origin, target = MAIN_TEST_META[test]
            n_c = n_cells[MAIN_NAMES[MAIN_TESTS[test][2]]]
            row = dict(space=sp, n_genes=int(spaces[sp].size), test=test,
                       role="primary" if test == "forward" else ("reverse" if test == "reverse" else "context"),
                       cell_line=line, state=state, origin=origin, target=target, n_cells=n_c)
            for stat in SAME_STATS:
                row[stat] = st[stat]
            row.update(dist_final=st["dist_final"], dist_start=st["dist_start"], v_norm=st["v_norm"],
                       d_norm=st["d_norm"], identity_mismatch=st["identity_mismatch"],
                       geometry_label=st["label"],
                       label_use="context" if (test != "forward" or sp == "FULL") else "reading (section 7)")
            nk = n1.get((sp, test))
            row.update(p_N1=nk["p"] if nk else np.nan, N1_null_mean=nk["null_mean"] if nk else np.nan,
                       N1_null_p95=nk["null_p95"] if nk else np.nan)
            for stat in SAME_STATS:
                lo, hi, ok = ci_with_flag(boot[f"same_main__{sp}__{test}__{stat}"], st[stat])
                row[f"{stat}_ci_lo"], row[f"{stat}_ci_hi"], row[f"{stat}_ci_valid"] = lo, hi, ok
                st[f"{stat}_ci_lo"], st[f"{stat}_ci_hi"], st[f"{stat}_ci_valid"] = lo, hi, ok
            if test == "forward":
                lo, hi, ok = ci_with_flag(boot[f"same_main__{sp}__asymmetry"], asym)
                row.update(asymmetry=asym, asymmetry_ci_lo=lo, asymmetry_ci_hi=hi, asymmetry_ci_valid=ok)
            same_rows.append(row)
        for m, st in zip(mixes, mix_pt[sp]):
            row = dict(space=sp, n_genes=int(spaces[sp].size), f=m["f"], n_aged=m["n_aged"],
                       n_young=m["n_young"], n_total=m["n_total"], noise_floor=bool(m["f"] == 0.0))
            for stat in SAME_STATS:
                row[stat] = st[stat]
            row.update(identity_mismatch=st["identity_mismatch"], geometry_label=st["label"])
            for stat in SAME_STATS:
                lo, hi, ok = ci_with_flag(boot[f"same_mix__{sp}__f{m['f']:.2f}__{stat}"], st[stat])
                row[f"{stat}_ci_lo"], row[f"{stat}_ci_hi"], row[f"{stat}_ci_valid"] = lo, hi, ok
            mix_rows.append(row)
    same_df = pd.DataFrame(same_rows)
    mix_df = pd.DataFrame(mix_rows)
    same_df.to_csv(GS_DIR / "same_stats.csv", index=False)
    mix_df.to_csv(GS_DIR / "same_posctrl.csv", index=False)

    gates = {sp: mixture_gate([r for r in mix_rows if r["space"] == sp]) for sp in SPACES}
    for sp in SPACES:
        log(f"[gate {sp}] pass={gates[sp]['passed']} failed={gates[sp]['failed']}")

    # ---- A3: FULL CI comparison (reported only)
    ref_pc = pd.read_csv(SAME_DIR / "posctrl.csv")
    ffw, frv = same_pt["FULL"]["forward"], same_pt["FULL"]["reverse"]
    ci_pairs = [
        ("frac_S_ci_lo", ffw["progress_ci_lo"], ref_same.get("frac_S_ci_lo")),
        ("frac_S_ci_hi", ffw["progress_ci_hi"], ref_same.get("frac_S_ci_hi")),
        ("delta_S_ci_lo", ffw["delta_ci_lo"], ref_same.get("delta_S_ci_lo")),
        ("delta_S_ci_hi", ffw["delta_ci_hi"], ref_same.get("delta_S_ci_hi")),
        ("frac_rev_ci_lo", frv["progress_ci_lo"], ref_same.get("frac_rev_ci_lo")),
        ("frac_rev_ci_hi", frv["progress_ci_hi"], ref_same.get("frac_rev_ci_hi")),
    ]
    full_mix = {r["f"]: r for r in mix_rows if r["space"] == "FULL"}
    ci_pairs.append(("p.frac50_ci_lo", full_mix[0.50]["progress_ci_lo"], ref_same["p"].get("frac50_ci_lo")))
    ci_pairs.append(("p.frac50_ci_hi", full_mix[0.50]["progress_ci_hi"], ref_same["p"].get("frac50_ci_hi")))
    for fval in (0.0, 0.10, 0.25):
        rr = ref_pc[np.isclose(ref_pc["f"].astype(float), fval)]
        ref_lo = float(rr["frac_ci_lo"].iloc[0]) if len(rr) else None
        ref_hi = float(rr["frac_ci_hi"].iloc[0]) if len(rr) else None
        ci_pairs.append((f"posctrl.csv f={fval:.2f} frac_ci_lo", full_mix[fval]["progress_ci_lo"], ref_lo))
        ci_pairs.append((f"posctrl.csv f={fval:.2f} frac_ci_hi", full_mix[fval]["progress_ci_hi"], ref_hi))
    ci_rows = []
    for name, new, ref in ci_pairs:
        diff = float(new) - float(ref) if (ref is not None and np.isfinite(new)) else np.nan
        ci_rows.append(dict(bound=name, rerun=float(new), reference=ref, diff=diff))
        log(f"[A3] {name}: rerun={float(new)!r} ref={ref!r} diff={diff:.3e}")
    rs.repro["ci_comparison_A3"] = dict(note="reported only; never a stop rule, gate or reading", rows=ci_rows)
    dump_json(GS_DIR / "repro_check.json", jsonable(rs.repro))

    # ---- section 3: GTEx-anchored test, per space
    write_progress(rs, "GTEx C3 per space, then GTEx state statistics")
    folds = c3_prepare(gtex, log)
    c3_rows, c3_pass = [], {}
    for sp in SPACES:
        rows, cboots, passed = c3_space(folds, spaces[sp], sp, log)
        c3_rows += rows
        c3_pass[sp] = passed
        for k, v in cboots.items():
            boot[f"c3__{sp}__fold{k}__cos"] = list(v)
    pd.DataFrame(c3_rows).to_csv(GS_DIR / "gtex_c3.csv", index=False)

    gtex_pt = {}
    gtex_null = {}
    for sp in SPACES:
        if not c3_pass[sp]:
            log(f"[gtex {sp}] c3_broken: nothing about GSE297234 is reported for {sp} in this test.")
            continue
        cols = spaces[sp]
        gtex_pt[sp] = gtex_point_block(donor_z, anchors, cols)
        u_x, _ = unit(anchors["c_young"][cols] - anchors["c_old"][cols])
        rng1 = np.random.default_rng(SEED)
        U1 = np.empty((N_RANDOM, cols.size))
        for i in range(N_RANDOM):
            U1[i] = rng1.permutation(u_x)
        rng2 = np.random.default_rng(SEED)
        Zx = gtex["Z"][:, cols]
        n = int(Zx.shape[0])
        n_half = n // 2
        U2 = np.empty((N_PERM, cols.size))
        for i in range(N_PERM):
            idx = rng2.permutation(n)
            ur, nrm = unit(Zx[idx[:n_half]].mean(0) - Zx[idx[n_half:2 * n_half]].mean(0))
            U2[i] = ur if nrm >= 1e-12 else np.nan
        del Zx
        for line in DONORS:
            for nm, rec in gtex_pt[sp][line].items():
                p1, nc1 = null_p(rec["cos_S"], U1, rec["d"])
                p2, nc2 = null_p(rec["cos_S"], U2, rec["d"])
                rec["p_N1"], rec["p_N2"] = p1, p2
                rec["pass_nulls"] = bool(np.isfinite(p1) and np.isfinite(p2) and p1 <= P_BAR and p2 <= P_BAR)
                rec["skipped"] = False
                rec["meets_toward"] = meets_toward(rec)
                gtex_null[(sp, line, nm)] = (float(np.nanmean(nc1)), float(np.nanpercentile(nc1, 95)),
                                             float(np.nanmean(nc2)), float(np.nanpercentile(nc2, 95)))
                nulls[f"null_gtex_N1__{sp}__{line}__{nm}"] = nc1
                nulls[f"null_gtex_N2__{sp}__{line}__{nm}"] = nc2
                log(f"[gtex {sp} {line} {nm}] cos_S={_fmt(rec['cos_S'])} progress={_fmt(rec['progress'])} "
                    f"perp={_fmt(rec['perp'])} delta_young={_fmt(rec['delta_young'])} "
                    f"delta_old={_fmt(rec['delta_old'])} p_N1={_fmt_p(p1)} p_N2={_fmt_p(p2)}")

    gtex_sp = [sp for sp in SPACES if c3_pass[sp]]
    if gtex_sp:
        for line in DONORS:
            dz = donor_z[line]
            rng = np.random.default_rng(BOOT_SEED + DONOR_BOOT_OFFSET[line])
            for b in range(N_BOOT):
                mats = []
                for nm in dz["names"]:
                    ix = dz["idx"][nm]
                    mats.append(sum_rows(dz["Y"], rng.choice(ix, size=int(ix.size), replace=True)))
                C = np.vstack(mats)
                logcpm, _ = tmm_logcpm_quiet(C)
                Zb, _ = zscore_frozen(logcpm, frozen, C)
                for sp in gtex_sp:
                    cols = spaces[sp]
                    c_y, c_o = anchors["c_young"][cols], anchors["c_old"][cols]
                    z0 = Zb[0, cols]
                    for j, nm in enumerate(dz["names"]):
                        if nm == ORIGIN_STATE:
                            continue
                        g = gtex_axis(z0, Zb[j, cols], c_y, c_o)
                        for stat in GTEX_STATS:
                            bput(f"gtex__{sp}__{line}__{nm}__{stat}", g[stat])
                if (b + 1) % 50 == 0:
                    log(f"[gtex boot {line}] {b + 1}/{N_BOOT}")

    gtex_rows = []
    gtex_keys = {}
    for sp in gtex_sp:
        gtex_keys[sp] = {}
        for line in DONORS:
            dz = donor_z[line]
            recs = {}
            for st_name in DEST_STATES:
                rec = gtex_pt[sp][line].get(st_name)
                if rec is None:
                    sk = next((s for s in dz["skipped"] if s["state"] == st_name), None)
                    recs[st_name] = dict(skipped=True, reason=sk["reason"] if sk else "",
                                         cos_S=np.nan, delta_young=np.nan, p_N1=np.nan, p_N2=np.nan)
                    gtex_rows.append(dict(space=sp, n_genes=int(spaces[sp].size), cell_line=line,
                                          state=st_name, n_cells=dz["n_by"].get(st_name, 0), skipped=True,
                                          skipped_reason=sk["reason"] if sk else ""))
                    continue
                recs[st_name] = rec
                row = dict(space=sp, n_genes=int(spaces[sp].size), cell_line=line, state=st_name,
                           role="contrast" if line == YOUNG_LINE else "claim_donor",
                           n_cells=int(dz["idx"][st_name].size), skipped=False, skipped_reason="")
                for k in ("cos_S", "progress", "perp", "dist_young", "dist_young_origin", "delta_young",
                          "rel_delta_young", "dist_old", "dist_old_origin", "delta_old", "rel_delta_old",
                          "d_norm", "w_norm", "p_N1", "p_N2", "pass_nulls", "meets_toward"):
                    row[k] = rec[k]
                nm1, p951, nm2, p952 = gtex_null[(sp, line, st_name)]
                row.update(N1_null_mean=nm1, N1_null_p95=p951, N2_null_mean=nm2, N2_null_p95=p952)
                for stat in GTEX_STATS:
                    lo, hi, ok = ci_with_flag(boot[f"gtex__{sp}__{line}__{st_name}__{stat}"], rec[stat])
                    row[f"{stat}_ci_lo"], row[f"{stat}_ci_hi"], row[f"{stat}_ci_valid"] = lo, hi, ok
                gtex_rows.append(row)
            key = reading_for_donor(line, recs, dz["n_by"], False, c3_pass=True)
            gtex_keys[sp][line] = key
            log(f"[gtex key {sp} {line}] {key['key']}")
    for sp in SPACES:
        if not c3_pass[sp]:
            gtex_keys[sp] = {"key": "c3_broken",
                             "text": f"C3 failed in {sp}; nothing about GSE297234 is reported for {sp} in the GTEx test."}
    pd.DataFrame(gtex_rows).to_csv(GS_DIR / "gtex_stats.csv", index=False)

    # ---- section 5: readout-vector cosine
    mem = gm["membership"]
    readout_rows = []
    for sp in SPACES:
        cols = spaces[sp]
        r_md = np.where(mem["MD"][cols], -1.0, 0.0)
        r_age = np.where(mem["age_up"][cols], -1.0, 0.0) + np.where(mem["age_down"][cols], 1.0, 0.0)
        rvecs = {"r_MD": r_md, "r_AGE": r_age}
        vecs = []
        for test in MAIN_TESTS:
            line, state, _, _ = MAIN_TEST_META[test]
            vecs.append(("same", test, line, state, same_pt[sp][test]["d"]))
        for m, st in zip(mixes, mix_pt[sp]):
            vecs.append(("same", "mixture", "", f"f={m['f']:.2f}", st["d"]))
        vecs.append(("same", "axis", "", "v = z(Y) - z(O)", same_pt[sp]["forward"]["v"]))
        u_x, _ = unit(anchors["c_young"][cols] - anchors["c_old"][cols])
        vecs.append(("gtex", "axis", "", "u(X)", u_x))
        if c3_pass[sp]:
            for line in DONORS:
                for nm, rec in gtex_pt[sp][line].items():
                    vecs.append(("gtex", "d0_to_state", line, nm, rec["d"]))
        for rname, r in rvecs.items():
            nnz = int(np.count_nonzero(r))
            if nnz < 1:
                continue
            for test_block, test, line, state, vec in vecs:
                readout_rows.append(dict(space=sp, readout=rname, n_nonzero=nnz, test=test_block,
                                         vector=test, cell_line=line, state=state, cos=cosine_full(vec, r)))
    pd.DataFrame(readout_rows).to_csv(GS_DIR / "readout_cos.csv", index=False)

    # ---- section 7: reading
    readings = {}
    for sp in ("MD", "AGE"):
        readings[sp] = fire_space_reading(sp, spaces[sp].size, n_cells["S"], same_pt[sp]["forward"],
                                          same_pt[sp]["reverse"], gates[sp])
        log(f"[reading {sp}] key={readings[sp]['key']}")
        log(readings[sp]["text"])
        log(readings[sp]["ci_line"])
    context_labels = {sp: {t: same_pt[sp][t]["label"] for t in ("reverse", "Pluri", "NonReprog")} for sp in SPACES}
    context_labels["FULL"]["forward"] = same_pt["FULL"]["forward"]["label"]
    rs.reading = dict(
        readings=readings, context_labels=context_labels, gates=gates,
        gtex_secondary_keys=gtex_keys,
        note="Only readings[MD] and readings[AGE] are the pre-registered reading. Context labels, the "
             "FULL gate and GTEx secondary keys never change it.",
    )
    dump_json(GS_DIR / "reading.json", jsonable(rs.reading))

    # ---- section 6 + A1: size baselines
    write_progress(rs, "size baselines (section 6 and Amendment 1)")
    base_rows = []
    base_sets = {}
    for mode in BASELINES:
        for sp in ("MD", "AGE"):
            cols = spaces[sp]
            sets, need = draw_baseline_sets(cols, bins, zero, mode)
            base_sets[f"{mode}__{sp}"] = sets
            inc_g = bool(c3_pass[sp])
            real = flat_point_stats(cols, ctx, inc_g, sp)
            rand = [flat_point_stats(np.sort(s), ctx, inc_g, f"{mode} {sp} draw{i}") for i, s in enumerate(sets)]
            overlap = np.array([np.isin(s, cols).sum() for s in sets], float)
            log(f"[baseline {mode} {sp}] draws={len(sets)} mean_overlap={overlap.mean():.2f} "
                f"({overlap.mean() / cols.size:.4f})")
            for key, rv in real.items():
                blk, line, state, stat = key
                base = dict(baseline=mode, space=sp, n_genes=int(cols.size), test=blk, cell_line=line,
                            state=state, stat=stat, mean_overlap=float(overlap.mean()),
                            mean_overlap_frac=float(overlap.mean() / cols.size))
                if stat == "geometry_label":
                    labs = [r[key] for r in rand]
                    counts = {lab: int(sum(1 for x in labs if x == lab))
                              for lab in ("APPROACH", "OVERSHOOT", "SIDEWAYS", "TIE", "NA")}
                    base.update(real=rv, label_counts=json.dumps(counts))
                else:
                    vals = np.array([r[key] for r in rand], float)
                    fin = vals[np.isfinite(vals)]
                    pos = (float(np.sum(fin < rv)) + 0.5 * float(np.sum(fin == rv))) / N_BASELINE
                    base.update(real=rv, position=pos, rand_p2_5=float(np.percentile(fin, 2.5)),
                                rand_p50=float(np.percentile(fin, 50)), rand_p97_5=float(np.percentile(fin, 97.5)),
                                n_finite=int(fin.size))
                base_rows.append(base)
    pd.DataFrame(base_rows).to_csv(GS_DIR / "size_baseline.csv", index=False)
    np.savez_compressed(GS_DIR / "size_baseline_sets.npz", **base_sets)

    np.savez_compressed(GS_DIR / "boot.npz",
                        **{k: np.asarray(v, float) for k, v in boot.items()},
                        **{k: np.asarray(v, float) for k, v in nulls.items()})

    rs.tables = dict(same=same_df, mix=mix_df, c3=pd.DataFrame(c3_rows), gtex=pd.DataFrame(gtex_rows),
                     readout=pd.DataFrame(readout_rows), base=pd.DataFrame(base_rows))
    rs.c3_pass = c3_pass
    rs.gates = gates
    rs.summary = dict(
        repro_pass=repro_pass,
        spaces={sp: int(spaces[sp].size) for sp in SPACES},
        n_cells=n_cells, anchors=rs.anchor_n, c3_pass=c3_pass,
        same={sp: {t: {s: same_pt[sp][t][s] for s in SAME_STATS + ("label",)} for t in MAIN_TESTS} for sp in SPACES},
        asymmetry={sp: same_pt[sp]["forward"]["progress"] - same_pt[sp]["reverse"]["progress"] for sp in SPACES},
        mixture={sp: [{s: st[s] for s in SAME_STATS} | {"f": m["f"]} for m, st in zip(mixes, mix_pt[sp])]
                 for sp in SPACES},
        gates={sp: {k: v for k, v in gates[sp].items()} for sp in SPACES},
        readings={sp: readings[sp]["key"] for sp in readings},
        gtex_claim_PartialReprog={sp: {k: gtex_pt[sp][AGED_LINE]["PartialReprog"][k] for k in GTEX_STATS}
                                  for sp in gtex_sp if "PartialReprog" in gtex_pt[sp][AGED_LINE]},
        zero_by_panel=zero_by_panel,
        seed=SEED, boot_seed=BOOT_SEED, n_perm=N_PERM, n_boot=N_BOOT, n_random=N_RANDOM, n_baseline=N_BASELINE,
    )
    dump_json(GS_DIR / "summary.json", jsonable(rs.summary))
    rs.man["status"] = "DONE"
    log("[run] DONE")


# --------------------------------------------------------------------------- files
def write_gene_map(gm, frozen, bins, zero):
    mem = gm["membership"]
    n = len(frozen["mu"])
    df = pd.DataFrame(dict(
        record="ruler_column", col=np.arange(n), ensembl=np.asarray(frozen["ensembl"]).astype(str),
        symbol=np.asarray(frozen["symbol"]).astype(str), mu=np.asarray(frozen["mu"], float),
        mu_bin=bins, zero_main_same_panel=zero, in_MD=mem["MD"], in_age_up=mem["age_up"],
        in_age_down=mem["age_down"],
    ))
    df["in_AGE"] = df.in_age_up | df.in_age_down
    ent = pd.DataFrame(gm["entry_rows"])
    out = pd.concat([df, ent], ignore_index=True, sort=False)
    out.to_csv(GS_DIR / "gene_map.csv", index=False)


def build_space_info(gm, bins, zero, edges):
    info = dict(mu_bin_edges=edges.tolist(), n_zero_ruler=int(zero.sum()), spaces={}, lists={})
    sets = dict(gm["spaces"])
    for name, cols in list(sets.items()) + [(k, gm["lists"][k]["cols"]) for k in ("age_up", "age_down")]:
        rec = dict(n_mapped=int(cols.size), n_nonzero_main_same=int((~zero[cols]).sum()),
                   n_zero_main_same=int(zero[cols].sum()),
                   per_bin=bin_counts(cols, bins), per_stratum_A1=strata_counts(cols, bins, zero))
        if name in SPACES:
            rec["cols"] = cols.tolist() if name != "FULL" else "all 23,485 ruler columns"
            info["spaces"][name] = rec
        else:
            info["lists"][name] = rec
    for k in LISTS:
        L = gm["lists"][k]
        info["lists"].setdefault(k, {})
        info["lists"][k].update(n_entries=L["n_entries"], n_unique=L["n_unique"], n_unmapped=L["n_unmapped"],
                                duplicates=L["dups"])
        if k == "MD":
            info["lists"][k].update(n_mapped=int(L["cols"].size),
                                    n_nonzero_main_same=int((~zero[L["cols"]]).sum()),
                                    n_zero_main_same=int(zero[L["cols"]].sum()))
    info["checks"] = gm["checks"]
    dump_json(GS_DIR / "spaces.json", jsonable(info))
    return info


def write_progress(rs, next_action):
    lines = [
        "# PROGRESS_GENESPACE", "",
        f"**Status:** {rs.man.get('status')}",
        f"**Next action:** {next_action}", "",
        f"- seeds {SEED} / {BOOT_SEED}; n_perm={N_PERM} n_random={N_RANDOM} n_boot={N_BOOT} n_baseline={N_BASELINE}",
        "- PREREG.flag (incl. Amendment 1) not rewritten; paper/ and paper_package/ not read.",
        "- Imports pure functions from src/same_run.py and src/toward_run.py; calls none of their writers.",
    ]
    if getattr(rs, "repro", None):
        lines.append(f"- FULL reproduction check (section 8): passed={rs.repro.get('passed')}")
    if rs.reading.get("readings"):
        for sp, r in rs.reading["readings"].items():
            lines.append(f"- reading {sp}: `{r['key']}`")
    elif rs.reading.get("key"):
        lines.append(f"- reading: `{rs.reading['key']}`")
    if rs.man.get("hashes_unchanged") is not None:
        lines.append(f"- protected files unchanged: {rs.man['hashes_unchanged']}")
    fails = rs.man.get("failures") or []
    if fails:
        lines += ["", "## Failures", ""]
        lines += [f"- **{f['step']}:** {f['message']}" for f in fails]
    lines += ["", "## Files", "", "- `src/genespace_run.py`", "- `results/genespace/`",
              "- `FINDINGS_GENESPACE.md`", "- `PROGRESS_GENESPACE.md`", ""]
    PROGRESS_PATH.write_text("\n".join(lines) + "\n", encoding="utf-8")


def _tbl(df, cols):
    cols = [c for c in cols if c in df.columns]
    return md_table(df[cols])


def write_findings(rs):
    prereg = getattr(rs, "prereg_text", None) or PREREG_FLAG.read_text(encoding="utf-8")
    man = rs.man
    L = [
        "# FINDINGS_GENESPACE — the SAME and TOWARD tests inside smaller gene spaces", "",
        f"**Status:** {man.get('status')}. Run valid (protected files unchanged): "
        f"{man.get('hashes_unchanged')}. Seeds {SEED} / {BOOT_SEED}. n_perm={N_PERM}, n_random={N_RANDOM}, "
        f"n_boot={N_BOOT}, n_baseline={N_BASELINE}.", "",
        "Reproduced by `src/genespace_run.py`. Only the section-7 reading at the end is the conclusion; "
        "every other number is **diagnostic**. No manuscript is edited; paper/ and paper_package/ are not read.", "",
        "## Pre-registration (verbatim copy of results/genespace/PREREG.flag)", "",
        "~~~~text", prereg.rstrip("\n"), "~~~~", "",
        "## Results", "",
        "### Run integrity", "",
    ]
    fails = man.get("failures") or []
    if fails:
        L += [f"- **STOP [{f['step']}]:** {f['message']}" for f in fails]
    else:
        L.append("- No STOP.")
    L.append(f"- SHA-256 of {len(man.get('hashes_before') or {})} protected files before and after: "
             f"unchanged = {man.get('hashes_unchanged')}.")
    if man.get("hashes_changed"):
        L.append(f"- **Changed files (run flagged invalid):** {man['hashes_changed']}")
    L.append("")

    si = getattr(rs, "space_info", None)
    if si:
        L += ["### Gene spaces (section 1, Amendment 1)", "",
              "Counts recomputed by the script match PREREG.flag (otherwise the run would have stopped). "
              "Zero = summed counts 0 in every row of the full-gene main SAME panel [O, Y, S, S_rev, Pluri, NonReprog].", ""]
        rows = []
        for kind, name, rec in ([("space", k, v) for k, v in si["spaces"].items()]
                                + [("list", k, v) for k, v in si["lists"].items()]):
            rows.append(dict(set=f"{kind} {name}", n_mapped_to_ruler=rec.get("n_mapped"),
                             n_nonzero_main_SAME=rec.get("n_nonzero_main_same"),
                             n_zero_main_SAME=rec.get("n_zero_main_same")))
        L += [md_table(pd.DataFrame(rows)), ""]
        zb = getattr(rs, "zero_by_panel", None)
        if zb:
            L += ["Genes at z=0 per space and panel (zero counts in every row of that panel):", "",
                  md_table(pd.DataFrame([dict(space=k, **v) for k, v in zb.items()])), ""]
        L += ["Per-stratum counts (24 GTEx-mu bins × zero/nonzero) are in `results/genespace/spaces.json`.", ""]

    rep = getattr(rs, "repro", None)
    if rep:
        L += ["### Reproduction check, FULL space (section 8; stop rule |new − reference| ≤ 1e-6)", "",
              f"Passed: **{rep['passed']}**.", "",
              md_table(pd.DataFrame(rep["point"]).assign(
                  rerun=lambda d: d.rerun.map(lambda x: f"{x:.9g}"),
                  reference=lambda d: d.reference.map(lambda x: f"{float(x):.9g}"),
                  diff=lambda d: d["diff"].map(lambda x: f"{x:.2e}"))), ""]
        a3 = rep.get("ci_comparison_A3")
        if a3:
            L += ["### FULL bootstrap CI comparison (Amendment 1, A3; reported only, not a stop rule)", "",
                  md_table(pd.DataFrame(a3["rows"]).assign(
                      rerun=lambda d: d.rerun.map(lambda x: f"{x:.9g}"),
                      reference=lambda d: d.reference.map(lambda x: "NA" if x is None else f"{float(x):.9g}"),
                      diff=lambda d: d["diff"].map(lambda x: f"{x:.2e}"))), ""]

    T = getattr(rs, "tables", None)
    if T is not None:
        L += ["### Same-platform test per space (diagnostic except where used by the section-7 reading)", "",
              "Forward: v = z(Y) − z(O), S = GM00731 PartialReprog. Reverse: v_rev = z(O) − z(Y), S_rev = GM23815 "
              "PartialReprog. Pluri / NonReprog: GM00731, forward frame (context). N1 p for forward and reverse. "
              "Geometry labels are context except MD/AGE forward, which the section-7 reading uses only after "
              "the gate. Raw delta is not compared across spaces.", ""]
        L += [_tbl(T["same"], ["space", "n_genes", "test", "n_cells", "cos", "progress", "delta", "perp",
                               "rel_delta", "geometry_label", "p_N1", "asymmetry"]), ""]
        ci = T["same"].copy()
        for s in ("progress", "delta", "perp", "cos", "rel_delta"):
            ci[f"{s}_CI"] = [f"{_fmt_ci(a, b)} {valid_str(v)}" for a, b, v in
                             zip(ci[f"{s}_ci_lo"], ci[f"{s}_ci_hi"], ci[f"{s}_ci_valid"])]
        ci["asymmetry_CI"] = [f"{_fmt_ci(a, b)} {valid_str(v)}" if t == "forward" else "" for a, b, v, t in
                              zip(ci.get("asymmetry_ci_lo"), ci.get("asymmetry_ci_hi"),
                                  ci.get("asymmetry_ci_valid"), ci["test"])]
        L += ["95% percentile bootstrap CIs (B=200; O and Y frozen; valid iff lo ≤ point ≤ hi):", "",
              _tbl(ci, ["space", "test", "progress_CI", "delta_CI", "perp_CI", "cos_CI", "rel_delta_CI",
                        "asymmetry_CI"]), ""]

        L += ["### Mixture positive control per space (f = 0 row = noise floor, Amendment 1 A2; not subtracted)", ""]
        mx = T["mix"].copy()
        for s in ("progress", "delta"):
            mx[f"{s}_CI"] = [f"{_fmt_ci(a, b)} {valid_str(v)}" for a, b, v in
                             zip(mx[f"{s}_ci_lo"], mx[f"{s}_ci_hi"], mx[f"{s}_ci_valid"])]
        L += [_tbl(mx, ["space", "f", "n_aged", "n_young", "progress", "delta", "perp", "cos",
                        "progress_CI", "delta_CI", "geometry_label"]), ""]
        grow = []
        for sp in SPACES:
            g = rs.gates[sp]
            grow.append(dict(space=sp, gate_pass=g["passed"], non_decreasing=g["monotonic"],
                             ci50_excludes_0=g["progress50_ci_excludes_0"], ci50_valid=g["progress50_ci_valid"],
                             delta50_lt_0=g["delta50_lt_0"],
                             use="reproduction check only" if sp == "FULL" else "section-7 gate"))
        L += ["Section-7 gate per space:", "", md_table(pd.DataFrame(grow)), ""]

        L += ["### GTEx-anchored test per space (secondary; diagnostic)", "",
              f"Anchors: YOUNG (20-29 ∪ 30-39) n={rs.anchor_n['n_young']}, OLD (60-69 ∪ 70-79) "
              f"n={rs.anchor_n['n_old']}, middle (40-49, 50-59) n={rs.anchor_n['n_middle']} (reported, not in anchors), "
              f"all donors n={rs.anchor_n['n_donors']}.", "",
              "C3 five-fold positive control (pass = all folds finite and cos_fold > 0.50):", "",
              _tbl(T["c3"].assign(cos_CI=[f"{_fmt_ci(a, b)} {valid_str(v)}" for a, b, v in
                                          zip(T["c3"].cos_ci_lo, T["c3"].cos_ci_hi, T["c3"].cos_ci_valid)]),
                   ["space", "fold", "n_young_train", "n_old_train", "n_old_test", "cos", "cos_CI", "pass_bar",
                    "construction_ok"]), "",
              "C3 pass per space: " + ", ".join(f"{sp}={rs.c3_pass[sp]}" for sp in SPACES) + ".", ""]
        g = T["gtex"]
        if len(g):
            gs = g[~g.skipped.astype(bool)].copy()
            L += ["State statistics (d = z(S) − z(d0 Fibroblast); v = w(X) = c_young(X) − c_old(X); "
                  "delta_young/delta_old computed directly):", "",
                  _tbl(gs, ["space", "cell_line", "state", "n_cells", "cos_S", "progress", "perp", "delta_young",
                            "rel_delta_young", "delta_old", "rel_delta_old", "p_N1", "p_N2", "meets_toward"]), ""]
            for s in ("cos_S", "delta_young", "delta_old", "progress"):
                gs[f"{s}_CI"] = [f"{_fmt_ci(a, b)} {valid_str(v)}" for a, b, v in
                                 zip(gs[f"{s}_ci_lo"], gs[f"{s}_ci_hi"], gs[f"{s}_ci_valid"])]
            L += ["Bootstrap CIs (Generator(20260918 + donor offset); origin and destinations resampled):", "",
                  _tbl(gs, ["space", "cell_line", "state", "cos_S_CI", "progress_CI", "delta_young_CI",
                            "delta_old_CI"]), ""]
            sk = g[g.skipped.astype(bool)]
            if len(sk):
                L += ["Skipped destination states (n < 30):", "",
                      _tbl(sk, ["space", "cell_line", "state", "n_cells", "skipped_reason"]), ""]
        L += ["Secondary key per donor (results/toward/PREREG.flag decision order in X; GM23815 is a contrast; "
              "never changes the section-7 reading):", ""]
        for sp in SPACES:
            k = rs.reading["gtex_secondary_keys"][sp]
            if "key" in k:
                L.append(f"- {sp}: `{k['key']}` — {k['text']}")
            else:
                for line in DONORS:
                    L.append(f"- {sp} {line}: `{k[line]['key']}` — {k[line]['text']}")
        L.append("")

        ro = T["readout"]
        L += ["### Readout-vector cosine (section 5; context only)", "",
              "Unweighted sign vectors (r_MD: −1 on MD genes; r_AGE: −1 age_up, +1 age_down) dotted with frozen-GTEx "
              "z of TMM pseudobulk displacements. This approximates the published AddModuleScore and does not "
              "reproduce it (AddModuleScore is per cell on log-normalized counts, subtracts expression-bin-matched "
              "control genes, and is summarized per cluster).", "",
              _tbl(ro, ["space", "readout", "n_nonzero", "test", "vector", "cell_line", "state", "cos"]), ""]

        b = T["base"]
        L += ["### Size baselines (section 6 and Amendment 1 A1; context, never changes a reading)", "",
              "`sec6_mu`: section 6 as written (24 equal-count GTEx-mu bins). `A1_mu_zero`: Amendment 1 "
              "(24 bins × zero/nonzero in the main SAME panel). Position = (#{random < real} + 0.5·#{random = real}) / 200. "
              "Selected rows; the full table is `results/genespace/size_baseline.csv`.", ""]
        sel = b[
            ((b.test == "same") & b.state.isin(["forward", "reverse"]) & b.stat.isin(["progress", "perp", "rel_delta", "cos"]))
            | ((b.test == "mixture") & (b.state == "f=0.50"))
            | ((b.test == "gtex") & (b.cell_line == AGED_LINE) & (b.state == "PartialReprog")
               & b.stat.isin(["cos_S", "progress", "perp", "rel_delta_young"]))
        ].copy()
        sel = sel[sel.stat != "geometry_label"]
        sel["real"] = sel["real"].astype(float)
        L += [_tbl(sel, ["baseline", "space", "test", "cell_line", "state", "stat", "real", "position",
                         "rand_p2_5", "rand_p50", "rand_p97_5"]), ""]
        lab = b[b.stat == "geometry_label"]
        L += ["Forward geometry label: real space vs counts among 200 random sets (no gate applied to random sets):", "",
              _tbl(lab, ["baseline", "space", "real", "label_counts", "mean_overlap", "mean_overlap_frac"]), ""]

    L += ["## Pre-registered reading (section 7; only what fired)", ""]
    rd = rs.reading
    if rd.get("key") == "reproduction_failed":
        L += [f"**key: `reproduction_failed`.** {rd.get('text')}", ""]
    elif rd.get("readings"):
        for sp in ("MD", "AGE"):
            r = rd["readings"][sp]
            L += [f"### {sp}: `{r['key']}`", "", r["text"], "", r["ci_line"], ""]
        L += ["MD and AGE are reported side by side; they are not averaged, ranked or combined. "
              "FULL is a reproduction check and gets no reading.", "",
              "Context geometry labels (never a reading): " + "; ".join(
                  f"{sp}: " + ", ".join(f"{t}={lab}" for t, lab in rd["context_labels"][sp].items())
                  for sp in SPACES) + ".", ""]
    else:
        L += ["No reading fired (run stopped before section 7).", ""]
    L += ["## Files", "", "- `src/genespace_run.py`", "- `results/genespace/` (manifest.json, gene_map.csv, spaces.json, "
          "repro_check.json, same_stats.csv, same_posctrl.csv, gtex_stats.csv, gtex_c3.csv, readout_cos.csv, "
          "size_baseline.csv, size_baseline_sets.npz, boot.npz [bootstrap draws and null cosines], reading.json, "
          "summary.json, run_report.txt)", "- `FINDINGS_GENESPACE.md`", "- `PROGRESS_GENESPACE.md`", ""]
    FINDINGS_PATH.write_text("\n".join(L) + "\n", encoding="utf-8")


def finish(rs, log):
    try:
        after = hash_protected()
        before = rs.man.get("hashes_before") or {}
        changed = sorted(k for k in set(before) | set(after) if before.get(k) != after.get(k))
        rs.man["hashes_after"] = after
        rs.man["hashes_unchanged"] = bool(before) and not changed
        rs.man["hashes_changed"] = changed
        if changed:
            rs.man["run_valid"] = False
            log(f"[hash] CHANGED protected files: {changed}. Run flagged invalid.")
        else:
            rs.man["run_valid"] = rs.man.get("status") == "DONE"
            log(f"[hash] protected files unchanged (n={len(after)}).")
    except Exception as e:  # noqa: BLE001
        rs.man["hashes_unchanged"] = None
        log(f"[hash] after-run hashing failed: {e}")
    dump_json(GS_DIR / "manifest.json", jsonable(rs.man))
    write_findings(rs)
    write_progress(rs, "done" if rs.man.get("status") == "DONE" else "run stopped; see failures")


def main():
    log = Logger(GS_DIR / "run_report.txt")
    rs = RunState()
    try:
        run(log, rs)
    except StopStep as e:
        rs.man["status"] = "STOP"
        rs.man["failures"].append(dict(step=e.step, message=e.message))
        log(f"[run] STOP [{e.step}] {e.message}")
        finish(rs, log)
        log.close()
        raise
    finish(rs, log)
    log.close()


if __name__ == "__main__":
    main()
