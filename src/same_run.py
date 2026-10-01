"""SAME — same-platform young target.

Does aged partially-reprogrammed tissue move toward the YOUNG donor's untreated
day-0 fibroblasts from the same experiment? And is that rejuvenation, or just
both donors converging?

Read-only on existing matrices, labels, and the frozen ruler. Refit nothing.
Re-cluster nothing. Rescore nothing. Re-label nothing. GSE325735 out of scope.
Reuse src/toward_run.py functions; do not change them.

Pre-registration must already exist at results/same/PREREG.flag (written before
any statistic). This script does not rewrite that block.

Usage: python src/same_run.py
"""
from __future__ import annotations

import sys
from pathlib import Path

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parent))
from config import RESULTS, ROOT  # noqa: E402
from brain_phase1_common import Logger, dump_json, load_json  # noqa: E402
from gtex_common import StopStep, EDGE_R_PRIOR  # noqa: E402
from fibro_common import FROZEN_RULER, jsonable  # noqa: E402
from fibro2_common import load_frozen_ruler  # noqa: E402
from fibro3_common import AGED_LINE, YOUNG_LINE, DAYS  # noqa: E402
from target_common import md_table, unit  # noqa: E402
from trajectory_common import permutation_p  # noqa: E402
from toward_run import (  # noqa: E402
    tmm_logcpm_quiet,
    zscore_frozen,
    ruler_score,
    percentile_ci,
    sum_rows,
    load_donor,
    require_file,
    _fmt,
    _fmt_p,
    _fmt_ci,
    cosine_align,
    log_columns,
)


SAME_DIR = RESULTS / "same"
SAME_DIR.mkdir(parents=True, exist_ok=True)
FINDINGS_PATH = ROOT / "FINDINGS_SAME.md"
PROGRESS_PATH = ROOT / "PROGRESS_SAME.md"
PREREG_FLAG = SAME_DIR / "PREREG.flag"
MANIFEST_PATH = SAME_DIR / "manifest.json"
TOWARD_ANCHORS = RESULTS / "toward" / "anchors.npz"

SAME_SEED = 20260914
BOOT_SEED = 20260918
N_PERM = 200
N_BOOT = 200
MIN_HALF = 500
P_BAR = 0.05
F_VALUES = (0.0, 0.10, 0.25, 0.50)
ORIGIN_STATE = "Fibroblast"
DEST_STATE = "PartialReprog"
PLURI_STATE = "Pluripotency"
NONREPROG_STATE = "NonReprog"

GSE325735_BAN = "GSE325735 is out of scope and is not opened."


def load_prereg_text():
    if not PREREG_FLAG.exists():
        raise StopStep(
            "prereg",
            "PREREG.flag missing — write the pre-registration before any statistic. Not computing.",
        )
    return PREREG_FLAG.read_text(encoding="utf-8").strip()


def load_manifest():
    if MANIFEST_PATH.exists():
        return load_json(MANIFEST_PATH)
    return dict(seed=SAME_SEED, boot_seed=BOOT_SEED, failures=[], status="INIT")


def save_manifest(man):
    dump_json(MANIFEST_PATH, jsonable(man))
    return MANIFEST_PATH


def record_failure(step, message, details=None):
    man = load_manifest()
    rec = dict(step=step, message=message, details=jsonable(details or {}))
    man.setdefault("failures", []).append(rec)
    man["status"] = "STOP"
    save_manifest(man)
    return rec


def cosine_full(a, b):
    """cos = dot(a, b) / (||a|| * ||b||). Matches the SAME pre-registration."""
    a = np.asarray(a, float).ravel()
    b = np.asarray(b, float).ravel()
    na = float(np.linalg.norm(a))
    nb = float(np.linalg.norm(b))
    if (not np.isfinite(na)) or (not np.isfinite(nb)) or na < 1e-12 or nb < 1e-12:
        return np.nan
    if (not np.isfinite(a).all()) or (not np.isfinite(b).all()):
        return np.nan
    return float(np.dot(a, b) / (na * nb))


def frac_along(d, v):
    d = np.asarray(d, float).ravel()
    v = np.asarray(v, float).ravel()
    vv = float(np.dot(v, v))
    if (not np.isfinite(vv)) or vv < 1e-12 or (not np.isfinite(d).all()) or (not np.isfinite(v).all()):
        return np.nan
    return float(np.dot(d, v) / vv)


def axis_stats(z_origin, z_target, z_S):
    z_origin = np.asarray(z_origin, float).ravel()
    z_target = np.asarray(z_target, float).ravel()
    z_S = np.asarray(z_S, float).ravel()
    v = z_target - z_origin
    d = z_S - z_origin
    rec = dict(
        v=v,
        d=d,
        v_norm=float(np.linalg.norm(v)),
        d_norm=float(np.linalg.norm(d)),
        cos=cosine_full(d, v),
        frac=frac_along(d, v),
        delta=float(np.linalg.norm(z_S - z_target) - np.linalg.norm(z_origin - z_target)),
        dist_to_target=float(np.linalg.norm(z_S - z_target)),
        dist_origin_to_target=float(np.linalg.norm(z_origin - z_target)),
    )
    return rec


def ci_excludes_zero(lo, hi):
    if not (np.isfinite(lo) and np.isfinite(hi)):
        return False
    return bool(lo > 0 or hi < 0)


def ci_positive(lo, hi):
    return bool(np.isfinite(lo) and np.isfinite(hi) and lo > 0)


def ci_negative(lo, hi):
    return bool(np.isfinite(lo) and np.isfinite(hi) and hi < 0)


def ci_includes_zero(lo, hi):
    if not (np.isfinite(lo) and np.isfinite(hi)):
        return False
    return bool(lo <= 0 <= hi)


def ci_covers_point(point, lo, hi):
    if not (np.isfinite(point) and np.isfinite(lo) and np.isfinite(hi)):
        return False
    return bool(lo <= point <= hi)


def check_gene_space(Y, frozen, tag):
    mu = np.asarray(frozen["mu"], float)
    if int(Y.shape[1]) != int(mu.shape[0]):
        raise StopStep(
            "genes",
            f"{tag}: rulerY n_genes={Y.shape[1]} frozen mu n={mu.shape[0]}. Not aligning ad hoc.",
        )


def split_half(idx, seed):
    idx = np.asarray(idx, int)
    rng = np.random.default_rng(int(seed))
    perm = rng.permutation(idx.size)
    n_a = int(idx.size // 2)
    half_a = np.sort(idx[perm[:n_a]])
    half_b = np.sort(idx[perm[n_a:]])
    return half_a, half_b


def d0_fibroblast_idx(obs):
    return np.flatnonzero(
        (obs["label"].astype(str) == ORIGIN_STATE) & (obs["day"].astype(int) == 0)
    )


def state_idx(obs, state):
    return np.flatnonzero(obs["label"].astype(str) == str(state))


def day_composition(obs, state, line):
    m = obs["label"].astype(str) == str(state)
    rec = dict(cell_line=line, state=state, n_cells=int(m.sum()))
    for d in DAYS:
        rec[f"n_d{d}"] = int(((obs["day"].astype(int) == int(d)) & m).sum())
    return rec


def mix_indices(aged_b, young_b, seed):
    """Cell-count mixtures. Total N = min(n_aged_B, n_young_B). Independent Generator(seed)."""
    aged_b = np.asarray(aged_b, int)
    young_b = np.asarray(young_b, int)
    n_tot = int(min(aged_b.size, young_b.size))
    rng = np.random.default_rng(int(seed))
    mixes = []
    for f in F_VALUES:
        n_young = int(np.round(float(f) * n_tot))
        n_aged = int(n_tot - n_young)
        if n_aged > aged_b.size or n_young > young_b.size:
            raise StopStep(
                "posctrl",
                f"mixture f={f} needs n_aged={n_aged} n_young={n_young} "
                f"but half-B has {aged_b.size}/{young_b.size}. Not relaxing N.",
            )
        take_a = (
            rng.choice(aged_b, size=n_aged, replace=False) if n_aged > 0
            else np.array([], dtype=int)
        )
        take_y = (
            rng.choice(young_b, size=n_young, replace=False) if n_young > 0
            else np.array([], dtype=int)
        )
        mixes.append(dict(
            f=float(f), n_total=n_tot, n_aged=n_aged, n_young=n_young,
            aged_idx=np.asarray(take_a, int), young_idx=np.asarray(take_y, int),
        ))
    return mixes, n_tot


def mix_sum(Y_aged, Y_young, aged_idx, young_idx):
    return sum_rows(Y_aged, aged_idx) + sum_rows(Y_young, young_idx)


def resample_idx(idx, rng):
    idx = np.asarray(idx, int)
    n = int(idx.size)
    if n == 0:
        return idx
    return rng.choice(idx, size=n, replace=True)


def panel_z(count_rows, frozen):
    C = np.vstack([np.asarray(r, np.float64).ravel() for r in count_rows])
    logcpm, nf = tmm_logcpm_quiet(C)
    Z, missing = zscore_frozen(logcpm, frozen, C)
    return Z, missing, nf, C


def n1_cos_p(d, v, seed, n_perm=N_PERM):
    """Permute v's entries across genes. Independent Generator(seed)."""
    d = np.asarray(d, float).ravel()
    v = np.asarray(v, float).ravel()
    cos_obs = cosine_full(d, v)
    rng = np.random.default_rng(int(seed))
    null = np.empty(int(n_perm), dtype=float)
    for i in range(int(n_perm)):
        null[i] = cosine_full(d, rng.permutation(v))
    p = permutation_p(cos_obs, null, greater=True)
    return float(p) if np.isfinite(p) else np.nan, cos_obs, null


def monotonic_increasing(vals):
    v = np.asarray(vals, float)
    if v.size < 2 or (not np.isfinite(v).all()):
        return False
    return bool(np.all(v[1:] >= v[:-1]))


def write_progress(next_action, stop, extra=""):
    man = load_manifest()
    split = load_json(SAME_DIR / "split.json") if (SAME_DIR / "split.json").exists() else {}
    psum = load_json(SAME_DIR / "posctrl_summary.json") if (SAME_DIR / "posctrl_summary.json").exists() else {}
    reading = load_json(SAME_DIR / "reading.json") if (SAME_DIR / "reading.json").exists() else {}
    lines = [
        "# PROGRESS_SAME",
        "",
        f"**STOP status:** {stop}",
        f"**Next action:** {next_action}",
        "",
        "## Seeds / gates",
        "",
        f"- seed `{SAME_SEED}`  boot `{BOOT_SEED}`  n_perm={N_PERM}  n_boot={N_BOOT}",
        "- primary: frac_S, delta_S, Asymmetry = frac_S − frac_rev; uncalibrated R² is never a gate",
        f"- frozen ruler: context column only; not the statistic; not refit; exists={FROZEN_RULER.exists()}",
        f"- MIN_HALF={MIN_HALF}  P_BAR={P_BAR}  f={list(F_VALUES)}",
        f"- GSE325735: out of scope",
        "- Refit nothing. Do not re-cluster, re-label, or re-tune Louvain or AddModuleScore.",
        "- Donors never pooled. States never averaged. PartialReprog pooled across timepoints only.",
        f"- PREREG.flag exists={PREREG_FLAG.exists()} (not rewritten)",
        "",
        "## Split (from disk)",
        "",
        f"- GM00731 d0 Fibroblast half A n={split.get('n_O')}  half B n={split.get('n_aged_B')}",
        f"- GM23815 d0 Fibroblast half A n={split.get('n_Y')}  half B n={split.get('n_young_B')}",
        f"- split_ok={split.get('ok')}",
        "",
        "## Positive control P",
        "",
        f"- pass={psum.get('pass')}  monotonic={psum.get('monotonic')}  "
        f"frac_f0.50_ci_excludes_0={psum.get('frac50_ci_excludes_0')}  "
        f"delta_f0.50<0={psum.get('delta50_neg')}",
        "",
        "## Reading (from disk)",
        "",
        f"- key=`{reading.get('key', 'not fired')}`  flags={reading.get('flags', [])}",
        "",
    ]
    if extra:
        lines += ["## Note", "", extra, ""]
    fails = man.get("failures") or []
    if fails:
        lines += ["## Failures (manifest)", ""]
        for f in fails:
            lines.append(f"- **{f.get('step')}:** {f.get('message')}")
        lines.append("")
    lines += [
        "## Files",
        "",
        "- `src/same_run.py`",
        "- `src/toward_run.py` (functions reused, not modified)",
        "- `results/same/`",
        "- `FINDINGS_SAME.md`",
        "- `PROGRESS_SAME.md`",
        "",
    ]
    PROGRESS_PATH.write_text("\n".join(lines) + "\n", encoding="utf-8")
    return PROGRESS_PATH


def write_findings():
    prereg = load_prereg_text()
    man = load_manifest()
    status = man.get("status", "unknown")
    split = load_json(SAME_DIR / "split.json") if (SAME_DIR / "split.json").exists() else {}
    psum = load_json(SAME_DIR / "posctrl_summary.json") if (SAME_DIR / "posctrl_summary.json").exists() else {}
    reading = load_json(SAME_DIR / "reading.json") if (SAME_DIR / "reading.json").exists() else {}
    summary = load_json(SAME_DIR / "summary.json") if (SAME_DIR / "summary.json").exists() else {}
    key = reading.get("key", "not fired")
    flags = reading.get("flags") or []
    flag_s = (", flags: " + ", ".join(f"`{x}`" for x in flags)) if flags else ""
    status_line = (
        f"**Status:** {status}. key=`{key}`{flag_s}. "
        f"Seed `{SAME_SEED}`. boot `{BOOT_SEED}`. n_perm={N_PERM}. n_boot={N_BOOT}. "
        f"P pass={psum.get('pass')}. Frozen ruler exists={FROZEN_RULER.exists()}."
    )
    lines = [
        "# FINDINGS_SAME — same-platform young target",
        "",
        status_line,
        "",
        "Does aged partially-reprogrammed tissue move toward the YOUNG donor's untreated "
        "day-0 fibroblasts from the same experiment? And is that rejuvenation, or just "
        "both donors converging?",
        "",
        "Does not modify any existing `FINDINGS_*.md`, `PROGRESS_*.md`, or `FALSIFICATION.md`. "
        "No GSE325735. Uncalibrated R² is never a gate. Frozen ruler is a context column only. "
        "Donors never pooled. States never averaged.",
        "",
        "Reproduced by `src/same_run.py`. Reuses functions from `src/toward_run.py` "
        f"(not modified). Frozen ruler: `{FROZEN_RULER}` exists={FROZEN_RULER.exists()}.",
        "",
        "## Pre-registration (verbatim, written before any SAME statistic)",
        "",
        prereg,
        "",
        f"Flag: `results/same/PREREG.flag` exists={PREREG_FLAG.exists()}.",
        "",
        "## STOP / failures",
        "",
        "Not substituting columns or repairing rows.",
        "",
    ]
    fails = man.get("failures") or []
    if fails:
        for f in fails:
            lines.append(f"- **{f.get('step')}:** {f.get('message')}")
        lines.append("")
    else:
        lines += ["None recorded.", ""]

    lines += ["## Split sizes (frozen before any statistic)", ""]
    lines.append(
        "Each donor's d0 Fibroblast cells split 50/50 with an independent "
        f"`Generator({SAME_SEED})`. Half A = anchors (O, Y). Half B = positive control P only. "
        f"STOP if any half has <{MIN_HALF} cells."
    )
    lines.append("")
    sp = SAME_DIR / "split_table.csv"
    if sp.exists():
        sdf = pd.read_csv(sp)
        lines += [md_table(sdf), ""]
    elif split:
        lines.append(
            f"- GM00731 half A (O) n={split.get('n_O')}; half B n={split.get('n_aged_B')}."
        )
        lines.append(
            f"- GM23815 half A (Y) n={split.get('n_Y')}; half B n={split.get('n_young_B')}."
        )
        lines.append("")
    if split.get("ok") is False:
        lines.append(
            f"**Split STOP:** a half has n<{MIN_HALF}. Not relaxing the minimum. "
            "P and downstream statistics are not reported."
        )
        lines.append("")
        lines += _files_touched()
        FINDINGS_PATH.write_text("\n".join(lines) + "\n", encoding="utf-8")
        return FINDINGS_PATH

    lines += ["## Positive control P (run BEFORE any PartialReprog vector)", ""]
    lines.append(
        "Mixtures of aged half-B d0 Fibroblast cells with young half-B d0 Fibroblast cells "
        f"at young fractions f={list(F_VALUES)}; cell counts, not weights; "
        f"total cells fixed = min(n_aged_B, n_young_B); seed `{SAME_SEED}`. "
        "Scored as S in the forward test against half-A anchors. "
        "P PASSES only if frac is monotonically increasing in f, AND frac at f=0.50 has "
        "bootstrap CI excluding 0, AND delta at f=0.50 < 0."
    )
    lines.append("")
    if psum:
        lines.append(
            f"- P pass={psum.get('pass')}. monotonic={psum.get('monotonic')}. "
            f"frac(f=0.50) CI excludes 0 = {psum.get('frac50_ci_excludes_0')} "
            f"{_fmt_ci(psum.get('frac50_ci_lo'), psum.get('frac50_ci_hi'))}. "
            f"delta(f=0.50)={_fmt(psum.get('delta50'))} < 0 = {psum.get('delta50_neg')}."
        )
        lines.append(
            f"- Noise floor f=0: frac={_fmt(psum.get('frac0'))} "
            f"delta={_fmt(psum.get('delta0'))} (should be ~0). "
            f"n_total={psum.get('n_total')}."
        )
        if psum.get("reason"):
            lines.append(f"- fail reason: {psum.get('reason')}")
        lines.append("")
    pp = SAME_DIR / "posctrl.csv"
    if pp.exists():
        pdf = pd.read_csv(pp)
        show = pdf.copy()
        for c in ("p_N1",):
            if c in show.columns:
                show[c] = [_fmt_p(v) for v in show[c]]
        keep = [c for c in (
            "f", "n_aged", "n_young", "n_total", "frac", "delta", "cos",
            "frac_ci_lo", "frac_ci_hi", "delta_ci_lo", "delta_ci_hi",
            "noise_floor",
        ) if c in show.columns]
        lines += [md_table(show[keep]), ""]
    if psum.get("pass") is False:
        lines.append(
            "**P failed. Statistic cannot detect a known approach. STOP. "
            "Nothing downstream is reported.**"
        )
        lines.append("")

    if psum.get("pass") is False:
        lines += ["## Pre-registered reading (only the outcome that fired)", ""]
        lines.append(f"**key: `{key}`.**")
        lines.append("")
        lines.append(reading.get("text") or "Positive control P failed. STOP.")
        lines.append("")
        lines += _limitations()
        lines += _files_touched()
        FINDINGS_PATH.write_text("\n".join(lines) + "\n", encoding="utf-8")
        return FINDINGS_PATH

    occ_p = SAME_DIR / "occupancy.csv"
    if occ_p.exists():
        lines += ["## Timepoint composition (PartialReprog pooled; not split)", ""]
        lines.append(
            "Labels: md2 argmax of mean AddModuleScore (mmc3 Reprog_cell_state_signatures), "
            "exactly as used in TOWARD. Not re-labelled. PartialReprog is pooled across "
            "timepoints; composition is reported and is not used to slice a better result. "
            "Half-B d0 Fibroblast cells are not in any destination vector."
        )
        lines.append("")
        odf = pd.read_csv(occ_p)
        lines += [md_table(odf), ""]

    stats_p = SAME_DIR / "stats.csv"
    if stats_p.exists():
        lines += ["## Forward, reverse, and asymmetry", ""]
        lines.append(
            "Forward: v = z(Y)−z(O), S = GM00731 PartialReprog (timepoints pooled). "
            "Reverse: donors swapped; S_rev = GM23815 PartialReprog (timepoints pooled). "
            "Asymmetry = frac_S − frac_rev. N1: 200 permutations of v (resp. v_rev) entries "
            f"(seed {SAME_SEED}); p = (#{'{cos ≥ observed}'} + 1)/201. "
            f"Bootstrap B={N_BOOT} seed {BOOT_SEED}; 95% percentile CIs. "
            "Uncalibrated R² is never a gate."
        )
        lines.append("")
        sdf = pd.read_csv(stats_p)
        show = sdf.copy()
        if "p_N1" in show.columns:
            show["p_N1"] = [_fmt_p(v) for v in show["p_N1"]]
        keep = [c for c in (
            "test", "origin", "target", "S", "n_cells_S",
            "cos", "frac", "delta",
            "p_N1", "frac_ci_lo", "frac_ci_hi", "delta_ci_lo", "delta_ci_hi",
            "d_norm", "v_norm",
        ) if c in show.columns]
        lines += [md_table(show[keep]), ""]
        if "asymmetry" in sdf.columns or True:
            lines.append(
                f"- Asymmetry = frac_S − frac_rev = {_fmt(summary.get('asymmetry'))} "
                f"CI {_fmt_ci(summary.get('asymmetry_ci_lo'), summary.get('asymmetry_ci_hi'))}."
            )
            lines.append("")

    ctx_p = SAME_DIR / "context.csv"
    if ctx_p.exists():
        lines += ["## Context (reported, never the reading)", ""]
        lines.append(
            "GM00731 Pluripotency and NonReprog scored as S in the forward test. "
            "Frozen-ruler score of every pseudobulk (Z @ w). "
            "cos(v, u_TOWARD) from `results/toward/anchors.npz` — how much of the "
            "same-platform donor gap looks like age as GTEx defines it. Report only."
        )
        lines.append("")
        cdf = pd.read_csv(ctx_p)
        show = cdf.copy()
        if "p_N1" in show.columns:
            show["p_N1"] = [_fmt_p(v) for v in show["p_N1"]]
        lines += [md_table(show), ""]
        if summary.get("cos_v_u_toward") is not None:
            lines.append(
                f"- cos(v, u_TOWARD)={_fmt(summary.get('cos_v_u_toward'))} "
                f"(report only; not a gate)."
            )
            lines.append("")

    lines += ["## Pre-registered reading (only the outcome that fired)", ""]
    lines.append(f"**key: `{key}`.**" + (f" Flag: `{flags[0]}`." if flags else ""))
    lines.append("")
    if flags:
        lines.append("Flags: " + ", ".join(f"`{x}`" for x in flags) + ".")
        lines.append("")
    lines.append(reading.get("text") or "")
    lines.append("")
    cond = reading.get("conditions") or {}
    if cond and key == "mixed":
        lines.append("Conditions that held:")
        for ck, cv in cond.items():
            lines.append(f"- {ck}: {cv}")
        lines.append("")

    lines += _limitations()
    lines += _files_touched()
    FINDINGS_PATH.write_text("\n".join(lines) + "\n", encoding="utf-8")
    return FINDINGS_PATH


def _limitations():
    return [
        "## Limitations",
        "",
        "- One young donor, so \"young\" is confounded with that donor.",
        "- Day-0 culture/passage differences between donors are part of v.",
        "- Two donors total. Donors are not averaged.",
        "- State labels are Louvain-cluster argmax of signature scores, not per-cell argmax.",
        "- TMM is among a handful of pseudobulk rows (this task's panel).",
        "- Bootstrap CIs re-TMM the panel, so they are CIs on the TMM-coupled statistic "
        "and need not contain the observed point estimate. Not a gate.",
        "- Half-A anchors are a random 50% of d0 Fibroblast cells; the complementary half "
        "is used only in P.",
        "- Uncalibrated R² is never a gate and is not computed here.",
        f"- {GSE325735_BAN}",
        "",
    ]


def _files_touched():
    return [
        "## Files touched",
        "",
        "- `src/same_run.py`",
        "- `FINDINGS_SAME.md`",
        "- `PROGRESS_SAME.md`",
        "- `results/same/PREREG.flag` (not rewritten if it already existed)",
        "- `results/same/split.json`",
        "- `results/same/split_table.csv`",
        "- `results/same/split_idx.npz`",
        "- `results/same/posctrl.csv`",
        "- `results/same/posctrl_summary.json`",
        "- `results/same/occupancy.csv`",
        "- `results/same/stats.csv`",
        "- `results/same/context.csv`",
        "- `results/same/reading.json`",
        "- `results/same/summary.json`",
        "- `results/same/manifest.json`",
        "- `results/same/run_report.txt`",
        "- Existing `FINDINGS_*.md`, `PROGRESS_*.md`, and `FALSIFICATION.md` were not "
        "modified (except FINDINGS_SAME.md and PROGRESS_SAME.md).",
        "",
    ]


def fire_reading(p_pass, fwd, rev, asym, asym_ci, p_n1, n_cells_S, pluri, conditions_extra=None):
    flags = []
    frac_lo, frac_hi = fwd.get("frac_ci_lo", np.nan), fwd.get("frac_ci_hi", np.nan)
    dlt = float(fwd["delta"]) if fwd.get("delta") is not None and np.isfinite(fwd.get("delta", np.nan)) else np.nan
    frac = float(fwd["frac"]) if fwd.get("frac") is not None and np.isfinite(fwd.get("frac", np.nan)) else np.nan
    if isinstance(asym_ci, dict):
        a_lo, a_hi = asym_ci.get("lo", np.nan), asym_ci.get("hi", np.nan)
    else:
        a_lo, a_hi = (np.nan, np.nan) if asym_ci is None else (float(asym_ci[0]), float(asym_ci[1]))

    frac_ci_pos = ci_positive(frac_lo, frac_hi)
    frac_ci_incl0 = ci_includes_zero(frac_lo, frac_hi)
    delta_neg = bool(np.isfinite(dlt) and dlt < 0)
    delta_nonneg = bool(np.isfinite(dlt) and dlt >= 0)
    n1_ok = bool(np.isfinite(p_n1) and p_n1 <= P_BAR)
    asym_ci_pos = ci_positive(a_lo, a_hi)
    asym_ci_incl0 = ci_includes_zero(a_lo, a_hi)
    asym_ci_neg = ci_negative(a_lo, a_hi)

    pluri_overshoot = False
    if pluri:
        plo, phi = pluri.get("frac_ci_lo", np.nan), pluri.get("frac_ci_hi", np.nan)
        pdlt = pluri.get("delta", np.nan)
        pluri_overshoot = bool(
            np.isfinite(pdlt) and pdlt < 0 and ci_positive(plo, phi)
        )

    conditions = dict(
        P_pass=bool(p_pass),
        frac_S_CI_excludes_0_positive=frac_ci_pos,
        delta_S_lt_0=delta_neg,
        N1_p_le_0_05=n1_ok,
        Asymmetry_CI_excludes_0_positive=asym_ci_pos,
        Asymmetry_CI_includes_0=asym_ci_incl0,
        Asymmetry_CI_negative=asym_ci_neg,
        delta_S_ge_0=delta_nonneg,
        frac_S_CI_includes_0=frac_ci_incl0,
        overshoot_pluri_delta_lt_0_and_frac_CI_gt_0=pluri_overshoot,
        frac_S=frac,
        delta_S=dlt,
        p_N1=p_n1,
        n_cells_S=int(n_cells_S) if n_cells_S is not None else None,
        asymmetry=asym,
        frac_S_CI=[frac_lo, frac_hi],
        asymmetry_CI=[a_lo, a_hi],
    )
    if conditions_extra:
        conditions.update(conditions_extra)

    overshoot_note = (
        " The forward result may reflect generic de-differentiation, not approach "
        "to young fibroblast."
    )

    if not p_pass:
        key = "pos_control_broken"
        text = (
            "Positive control P failed. The statistic cannot detect a known approach. STOP. "
            "Nothing downstream is reported."
        )
        if pluri_overshoot:
            flags.append("overshoot_also_closer")
            text = text + overshoot_note
        return dict(key=key, flags=flags, text=text, conditions=conditions)

    forward_core = bool(frac_ci_pos and delta_neg and n1_ok)

    if forward_core and asym_ci_pos:
        key = "toward_young_same_platform"
        text = (
            f"Aged partially-reprogrammed cells (n_cells={int(n_cells_S)}, frac_S={_fmt(frac)}) "
            "approach the young donor's untreated cells more than the reverse. "
            "One young donor, so \"young\" is confounded with that donor."
        )
    elif forward_core and (asym_ci_incl0 or asym_ci_neg):
        key = "convergence_not_rejuvenation"
        text = (
            "Both donors drift toward each other; reprogramming erases donor differences "
            "rather than age."
        )
    elif delta_nonneg or frac_ci_incl0:
        key = "no_approach_same_platform"
        text = (
            "The TOWARD negative holds with the platform gap removed. "
            f"GM00731 PartialReprog n_cells={int(n_cells_S) if n_cells_S is not None else 'NA'}, "
            f"frac_S={_fmt(frac)}, delta_S={_fmt(dlt)} "
            "(delta_S >= 0: not closer to the young donor)."
        )
    else:
        key = "mixed"
        text = (
            "Mixed: the pre-registered toward-young / convergence / no-approach conditions "
            "did not fire as a block. "
            f"frac_S CI excludes 0 (positive)={frac_ci_pos} {_fmt_ci(frac_lo, frac_hi)}; "
            f"delta_S<0={delta_neg} (delta_S={_fmt(dlt)}); "
            f"N1 p<=0.05={n1_ok} (p={_fmt_p(p_n1)}); "
            f"Asymmetry CI excludes 0 (positive)={asym_ci_pos} {_fmt_ci(a_lo, a_hi)}; "
            f"frac_S CI includes 0={frac_ci_incl0}; delta_S>=0={delta_nonneg}."
        )

    if pluri_overshoot:
        flags.append("overshoot_also_closer")
        text = text + overshoot_note
    return dict(key=key, flags=flags, text=text, conditions=conditions)


def run_posctrl(Y_aged, Y_young, o_idx, y_idx, aged_b, young_b, frozen, log):
    mixes, n_tot = mix_indices(aged_b, young_b, SAME_SEED)
    log(f"[P] n_total=min(n_aged_B, n_young_B)={n_tot}  n_mixes={len(mixes)} seed={SAME_SEED}")
    names = ["O", "Y"] + [f"mix_f{m['f']:.2f}" for m in mixes]
    counts = [
        sum_rows(Y_aged, o_idx),
        sum_rows(Y_young, y_idx),
    ]
    for m in mixes:
        counts.append(mix_sum(Y_aged, Y_young, m["aged_idx"], m["young_idx"]))
        log(f"[P] f={m['f']:.2f} n_aged={m['n_aged']} n_young={m['n_young']} n_total={m['n_total']}")
    Z, missing, nf, C = panel_z(counts, frozen)
    log(f"[P tmm] n_rows={C.shape[0]} names={names} nf_median={float(np.median(nf)):.4f} "
        f"prior.count={EDGE_R_PRIOR:g} n_missing_z0={int(np.asarray(missing).sum())}")
    z_O, z_Y = Z[0], Z[1]
    rows = []
    obs_stats = []
    for i, m in enumerate(mixes):
        st = axis_stats(z_O, z_Y, Z[2 + i])
        st.update(f=m["f"], n_aged=m["n_aged"], n_young=m["n_young"], n_total=m["n_total"])
        obs_stats.append(st)
        log(f"[P stat] f={m['f']:.2f} frac={_fmt(st['frac'])} delta={_fmt(st['delta'])} "
            f"cos={_fmt(st['cos'])}")
    ruler = ruler_score(Z, frozen)

    rng_b = np.random.default_rng(BOOT_SEED)
    boot = {m["f"]: dict(frac=[], delta=[], cos=[]) for m in mixes}
    # Frozen half-A anchors (same count vectors as the point estimate).
    # Resampling O/Y here rebuilds v and the origin of d from a different
    # high-p draw than the point estimate and shifts the percentile CI off
    # the observed frac (shared origin noise). Mix rows are still resampled.
    c_O_obs = counts[0]
    c_Y_obs = counts[1]
    for b in range(N_BOOT):
        b_counts = [c_O_obs, c_Y_obs]
        for m in mixes:
            b_counts.append(mix_sum(
                Y_aged, Y_young,
                resample_idx(m["aged_idx"], rng_b),
                resample_idx(m["young_idx"], rng_b),
            ))
        Zb, _, _, _ = panel_z(b_counts, frozen)
        zO, zY = Zb[0], Zb[1]
        for i, m in enumerate(mixes):
            st = axis_stats(zO, zY, Zb[2 + i])
            boot[m["f"]]["frac"].append(st["frac"])
            boot[m["f"]]["delta"].append(st["delta"])
            boot[m["f"]]["cos"].append(st["cos"])
        if (b + 1) % 50 == 0:
            log(f"[P boot] {b+1}/{N_BOOT}")

    for i, m in enumerate(mixes):
        st = obs_stats[i]
        flo, fhi = percentile_ci(boot[m["f"]]["frac"])
        dlo, dhi = percentile_ci(boot[m["f"]]["delta"])
        rows.append(dict(
            f=m["f"], n_aged=m["n_aged"], n_young=m["n_young"], n_total=m["n_total"],
            frac=st["frac"], delta=st["delta"], cos=st["cos"],
            d_norm=st["d_norm"], v_norm=st["v_norm"],
            frac_ci_lo=flo, frac_ci_hi=fhi,
            delta_ci_lo=dlo, delta_ci_hi=dhi,
            ruler_score_context=float(ruler[2 + i]),
            noise_floor=bool(m["f"] == 0.0),
        ))
    pdf = pd.DataFrame(rows)
    pdf.to_csv(SAME_DIR / "posctrl.csv", index=False)
    for r in rows:
        for stat in ("frac", "delta"):
            ok = ci_covers_point(r[stat], r[f"{stat}_ci_lo"], r[f"{stat}_ci_hi"])
            log(
                f"[P CI sanity] f={r['f']:.2f} {stat}={_fmt(r[stat])} "
                f"{_fmt_ci(r[f'{stat}_ci_lo'], r[f'{stat}_ci_hi'])} "
                f"{'OK' if ok else 'INVALID'}"
            )

    fracs = [float(r["frac"]) for r in rows]
    row50 = next(r for r in rows if abs(r["f"] - 0.50) < 1e-12)
    row0 = next(r for r in rows if abs(r["f"] - 0.0) < 1e-12)
    mono = monotonic_increasing(fracs)
    frac50_ex0 = ci_excludes_zero(row50["frac_ci_lo"], row50["frac_ci_hi"])
    delta50_neg = bool(np.isfinite(row50["delta"]) and row50["delta"] < 0)
    passed = bool(mono and frac50_ex0 and delta50_neg)
    reasons = []
    if not mono:
        reasons.append(f"frac not monotonically increasing in f: {[_fmt(x) for x in fracs]}")
    if not frac50_ex0:
        reasons.append(
            f"frac(f=0.50) CI does not exclude 0: {_fmt_ci(row50['frac_ci_lo'], row50['frac_ci_hi'])}"
        )
    if not delta50_neg:
        reasons.append(f"delta(f=0.50)={_fmt(row50['delta'])} is not < 0")
    summary = {
        "pass": passed,
        "monotonic": mono,
        "frac50_ci_excludes_0": frac50_ex0,
        "delta50_neg": delta50_neg,
        "frac50": row50["frac"],
        "delta50": row50["delta"],
        "frac50_ci_lo": row50["frac_ci_lo"],
        "frac50_ci_hi": row50["frac_ci_hi"],
        "frac0": row0["frac"],
        "delta0": row0["delta"],
        "n_total": n_tot,
        "fracs": fracs,
        "f_values": list(F_VALUES),
        "reason": "; ".join(reasons) if reasons else "",
        "n_boot": N_BOOT,
        "seed": SAME_SEED,
        "boot_seed": BOOT_SEED,
        "n_missing_z0": int(np.asarray(missing).sum()),
    }
    dump_json(SAME_DIR / "posctrl_summary.json", jsonable(summary))
    np.savez_compressed(
        SAME_DIR / "posctrl_z.npz",
        Z=Z, names=np.array(names), missing=np.asarray(missing),
    )
    log(f"[P] pass={passed} monotonic={mono} frac50_ci_excludes_0={frac50_ex0} "
        f"delta50_neg={delta50_neg} reason={summary['reason'] or 'NA'}")
    return pdf, summary


def run(log):
    prereg = load_prereg_text()
    log("=" * 100)
    log("SAME  seed=%s  boot=%s  n_perm=%s  n_boot=%s" % (
        SAME_SEED, BOOT_SEED, N_PERM, N_BOOT))
    log("=" * 100)
    log("FALSIFICATION.md and existing FINDINGS*.md / PROGRESS*.md are not modified "
        "(except FINDINGS_SAME.md and PROGRESS_SAME.md).")
    log(GSE325735_BAN)
    log("Primary statistics: frac_S, delta_S, Asymmetry. Uncalibrated R² is never a gate.")
    log("Frozen ruler is a context column only. Not refit.")
    log("PREREG.flag exists; not rewriting the block.")
    log(prereg[:400] + " ...")
    man = load_manifest()
    man["status"] = "RUNNING"
    man["failures"] = []
    save_manifest(man)
    write_progress("load donors and freeze 50/50 split", stop="RUNNING")

    require_file(FROZEN_RULER, "frozen_ruler")
    frozen = load_frozen_ruler()
    log(f"[frozen] {FROZEN_RULER} n_genes={len(np.asarray(frozen['mu']))} keys={list(frozen.files)}")

    obs_a, Y_a = load_donor(AGED_LINE, log)
    obs_y, Y_y = load_donor(YOUNG_LINE, log)
    check_gene_space(Y_a, frozen, AGED_LINE)
    check_gene_space(Y_y, frozen, YOUNG_LINE)
    if int(Y_a.shape[1]) != int(Y_y.shape[1]):
        raise StopStep("genes", f"aged n_genes={Y_a.shape[1]} young n_genes={Y_y.shape[1]}. Not aligning.")

    d0_a = d0_fibroblast_idx(obs_a)
    d0_y = d0_fibroblast_idx(obs_y)
    log(f"[d0 Fibroblast] {AGED_LINE} n={d0_a.size}  {YOUNG_LINE} n={d0_y.size}")
    o_idx, aged_b = split_half(d0_a, SAME_SEED)
    y_idx, young_b = split_half(d0_y, SAME_SEED)
    if set(o_idx) & set(aged_b):
        raise StopStep("split", f"{AGED_LINE}: half A and half B overlap. Not using overlapping cells.")
    if set(y_idx) & set(young_b):
        raise StopStep("split", f"{YOUNG_LINE}: half A and half B overlap. Not using overlapping cells.")

    split_rows = [
        dict(donor=AGED_LINE, half="A", n_cells=int(o_idx.size),
             role="O (aged origin; forward-test anchor)", used_in="anchors"),
        dict(donor=AGED_LINE, half="B", n_cells=int(aged_b.size),
             role="held out; P mixtures only", used_in="P"),
        dict(donor=YOUNG_LINE, half="A", n_cells=int(y_idx.size),
             role="Y (young target; forward-test anchor)", used_in="anchors"),
        dict(donor=YOUNG_LINE, half="B", n_cells=int(young_b.size),
             role="held out; P mixtures only", used_in="P"),
    ]
    split_df = pd.DataFrame(split_rows)
    split_df.to_csv(SAME_DIR / "split_table.csv", index=False)
    halves = [int(o_idx.size), int(aged_b.size), int(y_idx.size), int(young_b.size)]
    split_ok = bool(all(n >= MIN_HALF for n in halves))
    split_meta = dict(
        n_O=int(o_idx.size), n_Y=int(y_idx.size),
        n_aged_B=int(aged_b.size), n_young_B=int(young_b.size),
        n_d0_aged=int(d0_a.size), n_d0_young=int(d0_y.size),
        min_half=MIN_HALF, ok=split_ok, seed=SAME_SEED,
        note="independent Generator(seed) per donor; half A = n//2",
    )
    dump_json(SAME_DIR / "split.json", jsonable(split_meta))
    np.savez_compressed(
        SAME_DIR / "split_idx.npz",
        O=o_idx, Y=y_idx, aged_B=aged_b, young_B=young_b,
    )
    log(f"[split] O={o_idx.size} aged_B={aged_b.size} Y={y_idx.size} young_B={young_b.size} ok={split_ok}")
    if not split_ok:
        msg = (
            f"a half has n<{MIN_HALF}: O={o_idx.size} aged_B={aged_b.size} "
            f"Y={y_idx.size} young_B={young_b.size}. STOP. Not relaxing the minimum."
        )
        raise StopStep("split", msg, split_meta)

    write_progress("positive control P (no PartialReprog vector yet)", stop="SPLIT_FROZEN")
    pdf, psum = run_posctrl(Y_a, Y_y, o_idx, y_idx, aged_b, young_b, frozen, log)
    if not psum["pass"]:
        reading = fire_reading(
            False, dict(frac=np.nan, delta=np.nan, frac_ci_lo=np.nan, frac_ci_hi=np.nan),
            {}, np.nan, (np.nan, np.nan), np.nan, None, None,
        )
        dump_json(SAME_DIR / "reading.json", jsonable(reading))
        dump_json(SAME_DIR / "summary.json", jsonable(dict(
            p=psum, split=split_meta, c3_stop=False, p_stop=True,
            seed=SAME_SEED, boot_seed=BOOT_SEED, n_perm=N_PERM, n_boot=N_BOOT,
        )))
        raise StopStep(
            "P",
            "P failed. Statistic cannot detect a known approach. STOP. "
            "Nothing downstream is reported. " + (psum.get("reason") or ""),
            psum,
        )

    write_progress("PartialReprog forward/reverse after P pass", stop="P_PASS")

    # Destination indices — first time any PartialReprog vector is built.
    s_idx = state_idx(obs_a, DEST_STATE)
    srev_idx = state_idx(obs_y, DEST_STATE)
    pluri_idx = state_idx(obs_a, PLURI_STATE)
    nr_idx = state_idx(obs_a, NONREPROG_STATE)
    if int(s_idx.size) < 1:
        raise StopStep("S", f"{AGED_LINE} {DEST_STATE} n={s_idx.size}. Cannot build S.")
    if int(srev_idx.size) < 1:
        raise StopStep("S_rev", f"{YOUNG_LINE} {DEST_STATE} n={srev_idx.size}. Cannot build S_rev.")

    occ_rows = [
        day_composition(obs_a, ORIGIN_STATE, AGED_LINE),
        day_composition(obs_a, DEST_STATE, AGED_LINE),
        day_composition(obs_a, PLURI_STATE, AGED_LINE),
        day_composition(obs_a, NONREPROG_STATE, AGED_LINE),
        day_composition(obs_y, ORIGIN_STATE, YOUNG_LINE),
        day_composition(obs_y, DEST_STATE, YOUNG_LINE),
    ]
    occ_df = pd.DataFrame(occ_rows)
    occ_df.to_csv(SAME_DIR / "occupancy.csv", index=False)
    log("[occupancy] " + occ_df.to_string(index=False).replace("\n", " | "))

    panel_spec = [
        ("O", Y_a, o_idx, AGED_LINE, "anchor_halfA_d0_Fibroblast"),
        ("Y", Y_y, y_idx, YOUNG_LINE, "anchor_halfA_d0_Fibroblast"),
        ("S", Y_a, s_idx, AGED_LINE, "PartialReprog_pooled"),
        ("S_rev", Y_y, srev_idx, YOUNG_LINE, "PartialReprog_pooled"),
        ("Pluri", Y_a, pluri_idx, AGED_LINE, "Pluripotency_pooled"),
        ("NonReprog", Y_a, nr_idx, AGED_LINE, "NonReprog_pooled"),
    ]
    for name, Y, ix, line, role in panel_spec:
        log(f"[panel] {name} {line} n={ix.size} role={role}")
    # Indices are per-donor row positions. Only same-donor halves can leak.
    if set(o_idx) & set(aged_b):
        raise StopStep("leak", "O (GM00731 half A) overlaps aged half B. Half-B is P only.")
    if set(y_idx) & set(young_b):
        raise StopStep("leak", "Y (GM23815 half A) overlaps young half B. Half-B is P only.")

    counts = [sum_rows(Y, ix) for _, Y, ix, _, _ in panel_spec]
    names = [p[0] for p in panel_spec]
    Z, missing, nf, C = panel_z(counts, frozen)
    log(f"[tmm main] n_rows={C.shape[0]} names={names} nf_median={float(np.median(nf)):.4f} "
        f"prior.count={EDGE_R_PRIOR:g} n_missing_z0={int(np.asarray(missing).sum())}")
    pos = {nm: i for i, nm in enumerate(names)}
    z_O, z_Y = Z[pos["O"]], Z[pos["Y"]]
    z_S, z_Srev = Z[pos["S"]], Z[pos["S_rev"]]
    z_pluri, z_nr = Z[pos["Pluri"]], Z[pos["NonReprog"]]
    ruler = ruler_score(Z, frozen)

    fwd = axis_stats(z_O, z_Y, z_S)
    rev = axis_stats(z_Y, z_O, z_Srev)
    pluri = axis_stats(z_O, z_Y, z_pluri)
    nr = axis_stats(z_O, z_Y, z_nr)
    asym = float(fwd["frac"] - rev["frac"]) if (
        np.isfinite(fwd["frac"]) and np.isfinite(rev["frac"])
    ) else np.nan
    log(f"[forward] n_S={s_idx.size} cos={_fmt(fwd['cos'])} frac={_fmt(fwd['frac'])} "
        f"delta={_fmt(fwd['delta'])} ||d||={_fmt(fwd['d_norm'])} ||v||={_fmt(fwd['v_norm'])}")
    log(f"[reverse] n_Srev={srev_idx.size} cos={_fmt(rev['cos'])} frac={_fmt(rev['frac'])} "
        f"delta={_fmt(rev['delta'])}")
    log(f"[asymmetry] frac_S - frac_rev = {_fmt(asym)}")
    log(f"[context Pluri] n={pluri_idx.size} frac={_fmt(pluri['frac'])} delta={_fmt(pluri['delta'])}")
    log(f"[context NonReprog] n={nr_idx.size} frac={_fmt(nr['frac'])} delta={_fmt(nr['delta'])}")

    p1, _, null1 = n1_cos_p(fwd["d"], fwd["v"], SAME_SEED, N_PERM)
    p1r, _, null1r = n1_cos_p(rev["d"], rev["v"], SAME_SEED, N_PERM)
    log(f"[N1] forward p={_fmt_p(p1)} reverse p={_fmt_p(p1r)} n_perm={N_PERM} seed={SAME_SEED}")
    np.savez_compressed(SAME_DIR / "null_cos_forward.npz", null=null1, cos_obs=np.array(fwd["cos"]))
    np.savez_compressed(SAME_DIR / "null_cos_reverse.npz", null=null1r, cos_obs=np.array(rev["cos"]))

    require_file(TOWARD_ANCHORS, "toward_anchors")
    anc = np.load(TOWARD_ANCHORS)
    u_toward = np.asarray(anc["u"], float).ravel()
    log_columns(log, "toward_anchors", list(anc.files), str(TOWARD_ANCHORS))
    if u_toward.size != z_O.size:
        raise StopStep(
            "toward_u",
            f"u_TOWARD n={u_toward.size} z n={z_O.size}. Not aligning ad hoc.",
        )
    # cosine_align assumes unit u; u_TOWARD is stored unit. Also report full cosine.
    u_unit, u_n = unit(u_toward)
    cos_v_u = cosine_full(fwd["v"], u_toward)
    cos_v_u_align = cosine_align(fwd["v"], u_unit) if u_n >= 1e-12 else np.nan
    log(f"[context] cos(v, u_TOWARD)={_fmt(cos_v_u)} (unit-align={_fmt(cos_v_u_align)}; report only)")

    rng_b = np.random.default_rng(BOOT_SEED)
    boot_store = dict(
        frac_S=[], frac_rev=[], delta_S=[], delta_rev=[], asymmetry=[],
        frac_pluri=[], delta_pluri=[], frac_nr=[], delta_nr=[],
        cos_S=[], cos_rev=[],
    )
    # Frozen O/Y count vectors (same as the point-estimate panel). Destinations
    # (S, S_rev, Pluri, NonReprog) are resampled within the frozen cell sets.
    # Rebuilding O/Y from resampled cells was shifting every percentile CI off
    # the observed point (high-p shared origin noise in frac = dot(d,v)/||v||^2).
    for b in range(N_BOOT):
        b_counts = []
        for i, (_, Y, ix, _, _) in enumerate(panel_spec):
            if names[i] in ("O", "Y"):
                b_counts.append(counts[i])
            else:
                b_counts.append(sum_rows(Y, resample_idx(ix, rng_b)))
        Zb, _, _, _ = panel_z(b_counts, frozen)
        zO, zY = Zb[pos["O"]], Zb[pos["Y"]]
        st_f = axis_stats(zO, zY, Zb[pos["S"]])
        st_r = axis_stats(zY, zO, Zb[pos["S_rev"]])
        st_p = axis_stats(zO, zY, Zb[pos["Pluri"]])
        st_n = axis_stats(zO, zY, Zb[pos["NonReprog"]])
        boot_store["frac_S"].append(st_f["frac"])
        boot_store["frac_rev"].append(st_r["frac"])
        boot_store["delta_S"].append(st_f["delta"])
        boot_store["delta_rev"].append(st_r["delta"])
        boot_store["asymmetry"].append(
            (st_f["frac"] - st_r["frac"]) if (np.isfinite(st_f["frac"]) and np.isfinite(st_r["frac"])) else np.nan
        )
        boot_store["frac_pluri"].append(st_p["frac"])
        boot_store["delta_pluri"].append(st_p["delta"])
        boot_store["frac_nr"].append(st_n["frac"])
        boot_store["delta_nr"].append(st_n["delta"])
        boot_store["cos_S"].append(st_f["cos"])
        boot_store["cos_rev"].append(st_r["cos"])
        if (b + 1) % 50 == 0:
            log(f"[boot] {b+1}/{N_BOOT}")
    np.savez_compressed(SAME_DIR / "boot.npz", **{k: np.asarray(v, float) for k, v in boot_store.items()})

    def _ci(key):
        return percentile_ci(boot_store[key])

    fwd["frac_ci_lo"], fwd["frac_ci_hi"] = _ci("frac_S")
    fwd["delta_ci_lo"], fwd["delta_ci_hi"] = _ci("delta_S")
    fwd["cos_ci_lo"], fwd["cos_ci_hi"] = _ci("cos_S")
    rev["frac_ci_lo"], rev["frac_ci_hi"] = _ci("frac_rev")
    rev["delta_ci_lo"], rev["delta_ci_hi"] = _ci("delta_rev")
    rev["cos_ci_lo"], rev["cos_ci_hi"] = _ci("cos_rev")
    pluri["frac_ci_lo"], pluri["frac_ci_hi"] = _ci("frac_pluri")
    pluri["delta_ci_lo"], pluri["delta_ci_hi"] = _ci("delta_pluri")
    nr["frac_ci_lo"], nr["frac_ci_hi"] = _ci("frac_nr")
    nr["delta_ci_lo"], nr["delta_ci_hi"] = _ci("delta_nr")
    a_lo, a_hi = _ci("asymmetry")

    sanity_rows = [
        ("forward frac_S", fwd["frac"], fwd["frac_ci_lo"], fwd["frac_ci_hi"]),
        ("forward delta_S", fwd["delta"], fwd["delta_ci_lo"], fwd["delta_ci_hi"]),
        ("reverse frac_rev", rev["frac"], rev["frac_ci_lo"], rev["frac_ci_hi"]),
        ("reverse delta_rev", rev["delta"], rev["delta_ci_lo"], rev["delta_ci_hi"]),
        ("Asymmetry", asym, a_lo, a_hi),
        ("Pluri frac", pluri["frac"], pluri["frac_ci_lo"], pluri["frac_ci_hi"]),
        ("Pluri delta", pluri["delta"], pluri["delta_ci_lo"], pluri["delta_ci_hi"]),
        ("NonReprog frac", nr["frac"], nr["frac_ci_lo"], nr["frac_ci_hi"]),
        ("NonReprog delta", nr["delta"], nr["delta_ci_lo"], nr["delta_ci_hi"]),
    ]
    for tag, pt, lo, hi in sanity_rows:
        ok = ci_covers_point(pt, lo, hi)
        log(f"[CI sanity] {tag}={_fmt(pt)} {_fmt_ci(lo, hi)} {'OK' if ok else 'INVALID'}")

    stats_df = pd.DataFrame([
        dict(
            test="forward", origin=f"{AGED_LINE} d0 Fibroblast half A",
            target=f"{YOUNG_LINE} d0 Fibroblast half A",
            S=f"{AGED_LINE} PartialReprog pooled", n_cells_S=int(s_idx.size),
            cos=fwd["cos"], frac=fwd["frac"], delta=fwd["delta"],
            p_N1=p1, frac_ci_lo=fwd["frac_ci_lo"], frac_ci_hi=fwd["frac_ci_hi"],
            delta_ci_lo=fwd["delta_ci_lo"], delta_ci_hi=fwd["delta_ci_hi"],
            d_norm=fwd["d_norm"], v_norm=fwd["v_norm"],
            ruler_score_S=float(ruler[pos["S"]]),
        ),
        dict(
            test="reverse", origin=f"{YOUNG_LINE} d0 Fibroblast half A",
            target=f"{AGED_LINE} d0 Fibroblast half A",
            S=f"{YOUNG_LINE} PartialReprog pooled", n_cells_S=int(srev_idx.size),
            cos=rev["cos"], frac=rev["frac"], delta=rev["delta"],
            p_N1=p1r, frac_ci_lo=rev["frac_ci_lo"], frac_ci_hi=rev["frac_ci_hi"],
            delta_ci_lo=rev["delta_ci_lo"], delta_ci_hi=rev["delta_ci_hi"],
            d_norm=rev["d_norm"], v_norm=rev["v_norm"],
            ruler_score_S=float(ruler[pos["S_rev"]]),
        ),
    ])
    stats_df.to_csv(SAME_DIR / "stats.csv", index=False)

    ctx_df = pd.DataFrame([
        dict(pseudobulk="O", cell_line=AGED_LINE, state="Fibroblast_d0_halfA",
             n_cells=int(o_idx.size), role="origin_forward",
             ruler_score_context=float(ruler[pos["O"]]),
             frac=np.nan, delta=np.nan, cos=np.nan,
             frac_ci_lo=np.nan, frac_ci_hi=np.nan, delta_ci_lo=np.nan, delta_ci_hi=np.nan),
        dict(pseudobulk="Y", cell_line=YOUNG_LINE, state="Fibroblast_d0_halfA",
             n_cells=int(y_idx.size), role="target_forward",
             ruler_score_context=float(ruler[pos["Y"]]),
             frac=np.nan, delta=np.nan, cos=np.nan,
             frac_ci_lo=np.nan, frac_ci_hi=np.nan, delta_ci_lo=np.nan, delta_ci_hi=np.nan),
        dict(pseudobulk="S", cell_line=AGED_LINE, state="PartialReprog",
             n_cells=int(s_idx.size), role="forward_S",
             ruler_score_context=float(ruler[pos["S"]]),
             frac=fwd["frac"], delta=fwd["delta"], cos=fwd["cos"],
             frac_ci_lo=fwd["frac_ci_lo"], frac_ci_hi=fwd["frac_ci_hi"],
             delta_ci_lo=fwd["delta_ci_lo"], delta_ci_hi=fwd["delta_ci_hi"]),
        dict(pseudobulk="S_rev", cell_line=YOUNG_LINE, state="PartialReprog",
             n_cells=int(srev_idx.size), role="reverse_S",
             ruler_score_context=float(ruler[pos["S_rev"]]),
             frac=rev["frac"], delta=rev["delta"], cos=rev["cos"],
             frac_ci_lo=rev["frac_ci_lo"], frac_ci_hi=rev["frac_ci_hi"],
             delta_ci_lo=rev["delta_ci_lo"], delta_ci_hi=rev["delta_ci_hi"]),
        dict(pseudobulk="Pluri", cell_line=AGED_LINE, state="Pluripotency",
             n_cells=int(pluri_idx.size), role="overshoot_reference_forward",
             ruler_score_context=float(ruler[pos["Pluri"]]),
             frac=pluri["frac"], delta=pluri["delta"], cos=pluri["cos"],
             frac_ci_lo=pluri["frac_ci_lo"], frac_ci_hi=pluri["frac_ci_hi"],
             delta_ci_lo=pluri["delta_ci_lo"], delta_ci_hi=pluri["delta_ci_hi"]),
        dict(pseudobulk="NonReprog", cell_line=AGED_LINE, state="NonReprog",
             n_cells=int(nr_idx.size), role="off_trajectory_forward",
             ruler_score_context=float(ruler[pos["NonReprog"]]),
             frac=nr["frac"], delta=nr["delta"], cos=nr["cos"],
             frac_ci_lo=nr["frac_ci_lo"], frac_ci_hi=nr["frac_ci_hi"],
             delta_ci_lo=nr["delta_ci_lo"], delta_ci_hi=nr["delta_ci_hi"]),
    ])
    ctx_df.to_csv(SAME_DIR / "context.csv", index=False)

    np.savez_compressed(
        SAME_DIR / "panel_z.npz",
        Z=Z, names=np.array(names), missing=np.asarray(missing),
        v=fwd["v"], v_rev=rev["v"], u_toward=u_toward,
        cos_v_u_toward=np.array(cos_v_u),
    )

    reading = fire_reading(
        True, fwd, rev, asym, (a_lo, a_hi), p1, int(s_idx.size), pluri,
    )
    dump_json(SAME_DIR / "reading.json", jsonable(reading))
    log(f"[reading] key={reading['key']} flags={reading['flags']}")
    log(reading["text"])

    summary = dict(
        p=psum, split=split_meta,
        n_S=int(s_idx.size), n_Srev=int(srev_idx.size),
        n_pluri=int(pluri_idx.size), n_nonreprog=int(nr_idx.size),
        frac_S=fwd["frac"], delta_S=fwd["delta"], cos_S=fwd["cos"], p_N1=p1,
        frac_rev=rev["frac"], delta_rev=rev["delta"], cos_rev=rev["cos"], p_N1_rev=p1r,
        asymmetry=asym, asymmetry_ci_lo=a_lo, asymmetry_ci_hi=a_hi,
        frac_S_ci_lo=fwd["frac_ci_lo"], frac_S_ci_hi=fwd["frac_ci_hi"],
        delta_S_ci_lo=fwd["delta_ci_lo"], delta_S_ci_hi=fwd["delta_ci_hi"],
        frac_rev_ci_lo=rev["frac_ci_lo"], frac_rev_ci_hi=rev["frac_ci_hi"],
        delta_rev_ci_lo=rev["delta_ci_lo"], delta_rev_ci_hi=rev["delta_ci_hi"],
        pluri_frac=pluri["frac"], pluri_delta=pluri["delta"],
        pluri_frac_ci_lo=pluri["frac_ci_lo"], pluri_frac_ci_hi=pluri["frac_ci_hi"],
        nr_frac=nr["frac"], nr_delta=nr["delta"],
        cos_v_u_toward=cos_v_u,
        ruler={nm: float(ruler[pos[nm]]) for nm in names},
        n_missing_z0=int(np.asarray(missing).sum()),
        n_genes=int(Z.shape[1]),
        key=reading["key"], flags=reading["flags"],
        seed=SAME_SEED, boot_seed=BOOT_SEED, n_perm=N_PERM, n_boot=N_BOOT,
        min_half=MIN_HALF, p_bar=P_BAR,
    )
    dump_json(SAME_DIR / "summary.json", jsonable(summary))
    man = load_manifest()
    man["status"] = "DONE"
    man["key"] = reading["key"]
    man["flags"] = reading["flags"]
    save_manifest(man)
    write_findings()
    write_progress("done", stop="DONE")
    log("[run] DONE")
    return dict(reading=reading, p=psum)


def main():
    log = Logger(SAME_DIR / "run_report.txt")
    try:
        run(log)
    except StopStep as e:
        record_failure(e.step, e.message, e.details)
        log(f"[run] STOP [{e.step}] {e.message}")
        try:
            write_findings()
        except Exception as fe:
            log(f"[run] findings write after STOP failed: {fe}")
        try:
            write_progress(f"fix STOP [{e.step}]", stop=f"STOP {e.step}: {e.message}")
        except Exception:
            pass
        raise
    finally:
        log.close()


if __name__ == "__main__":
    main()
