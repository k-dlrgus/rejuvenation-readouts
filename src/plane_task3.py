"""Task 3 — one age–identity plane figure per donor, three instruments."""
from __future__ import annotations

import sys
from pathlib import Path

import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

sys.path.insert(0, str(Path(__file__).resolve().parent))
from plane_common import (  # noqa: E402
    PLANE_DIR, PLANE_FIG, AGED_LINE, YOUNG_LINE, TRAJECTORY, OFF_TRAJECTORY,
    AGE_INSTRUMENTS, INST_LABELS, IDENTITY_KEY, PREREG_TASK3, PREREG_TASK3_FLAG,
    StopStep, Logger, dump_json, jsonable, plane_log_banner,
    load_manifest, save_manifest, progress_snapshot,
)

STATE_COLOR = {
    "Fibroblast": "#4C78A8",
    "PartialReprog": "#F58518",
    "EarlyPluripotency": "#54A24B",
    "Pluripotency": "#E45756",
    "NonReprog": "#9D755D",
}
SHORT = {
    "Fibroblast": "Fib",
    "PartialReprog": "PR",
    "EarlyPluripotency": "EP",
    "Pluripotency": "Pluri",
    "NonReprog": "NR",
}


def _is_scored(r):
    v = r.get("scored") if hasattr(r, "get") else r["scored"]
    if isinstance(v, str):
        return v.strip().lower() in ("true", "1", "yes")
    if isinstance(v, (float, np.floating)) and not np.isfinite(v):
        return False
    return bool(v)


def _row(df, line, inst, st):
    hit = df[(df.cell_line == line) & (df.instrument == inst) & (df.label == st)]
    if not len(hit):
        return None
    return hit.iloc[0]


def _xyerr(r):
    x = float(r.identity_median)
    y = float(r.age_median)
    xlo, xhi = float(r.identity_ci_lo), float(r.identity_ci_hi)
    ylo, yhi = float(r.age_ci_lo), float(r.age_ci_hi)
    xerr = np.array([[x - xlo], [xhi - x]], dtype=float)
    yerr = np.array([[y - ylo], [yhi - y]], dtype=float)
    xerr = np.clip(xerr, 0, None)
    yerr = np.clip(yerr, 0, None)
    return x, y, xerr, yerr


def _xlim_for_donor(mdf, line):
    g = mdf[(mdf.cell_line == line) & (mdf.scored.astype(str).str.lower().isin(("true", "1")))]
    if not len(g):
        g = mdf[mdf.cell_line == line]
    xs = []
    for _, r in g.iterrows():
        for v in (r.get("identity_median"), r.get("identity_ci_lo"), r.get("identity_ci_hi")):
            try:
                fv = float(v)
            except (TypeError, ValueError):
                continue
            if np.isfinite(fv):
                xs.append(fv)
    if not xs:
        return None
    lo, hi = float(min(xs)), float(max(xs))
    pad = 0.08 * (hi - lo if hi > lo else 1.0)
    return lo - pad, hi + pad


def _ylim_observed(mdf, line, inst):
    g = mdf[(mdf.cell_line == line) & (mdf.instrument == inst)]
    ys = []
    for _, r in g.iterrows():
        if not _is_scored(r):
            continue
        for v in (r.get("age_median"), r.get("age_ci_lo"), r.get("age_ci_hi")):
            try:
                fv = float(v)
            except (TypeError, ValueError):
                continue
            if np.isfinite(fv):
                ys.append(fv)
    if not ys:
        return None
    lo, hi = float(min(ys)), float(max(ys))
    pad = 0.12 * (hi - lo if hi > lo else 1.0)
    return lo - pad, hi + pad


def plot_donor(line, mdf, env, log):
    xlim = _xlim_for_donor(mdf, line)
    fig, axes = plt.subplots(1, 3, figsize=(11.2, 3.6), sharex=True)
    fig.suptitle(f"{line}  x = {IDENTITY_KEY} (pluripotency−fibroblast)  NonReprog off-trajectory",
                 fontsize=10)
    for ax, inst in zip(axes, AGE_INSTRUMENTS):
        ax.set_title(INST_LABELS.get(inst, inst), fontsize=10)
        ax.set_xlabel("identity")
        ax.set_ylabel("age")
        if xlim is not None:
            ax.set_xlim(*xlim)
        ylim = _ylim_observed(mdf, line, inst)
        if ylim is not None:
            ax.set_ylim(*ylim)
        # envelope: vertical 2.5–97.5 at each scored state's observed identity
        env_label_used = False
        xs, elo, ehi = [], [], []
        for st in list(TRAJECTORY) + [OFF_TRAJECTORY]:
            e = _row(env, line, inst, st)
            r = _row(mdf, line, inst, st)
            if e is None or r is None or not _is_scored(r):
                continue
            x = float(r.identity_median)
            a = float(e.age_p025)
            b = float(e.age_p975)
            if not (np.isfinite(x) and np.isfinite(a) and np.isfinite(b)):
                continue
            ax.plot(
                [x, x], [a, b], color="#9ecae1", alpha=0.55, linewidth=7,
                zorder=1, solid_capstyle="round",
                label="null envelope 2.5–97.5" if not env_label_used else None,
            )
            env_label_used = True
            if st in TRAJECTORY:
                xs.append(x)
                elo.append(a)
                ehi.append(b)
        if len(xs) >= 2 and ylim is not None:
            obs_span = abs(ylim[1] - ylim[0])
            env_span = float(max(ehi) - min(elo))
            if obs_span > 0 and env_span <= 3.0 * obs_span:
                ax.fill_between(
                    xs, elo, ehi, color="#9ecae1", alpha=0.25, linewidth=0, zorder=0,
                )
        # trajectory states
        tx, ty = [], []
        for st in TRAJECTORY:
            r = _row(mdf, line, inst, st)
            if r is None or not _is_scored(r):
                ax.plot([], [], "o", color=STATE_COLOR[st], label=f"{SHORT[st]} not scored")
                continue
            x, y, xerr, yerr = _xyerr(r)
            if not (np.isfinite(x) and np.isfinite(y)):
                continue
            ax.errorbar(
                x, y, xerr=xerr, yerr=yerr, fmt="o", color=STATE_COLOR[st],
                ecolor=STATE_COLOR[st], elinewidth=1.0, capsize=2.5, markersize=6,
                zorder=3, label=SHORT[st],
            )
            tx.append(x)
            ty.append(y)
        if len(tx) >= 2:
            ax.plot(tx, ty, "-", color="#333333", linewidth=1.0, zorder=2)
        # NonReprog
        r = _row(mdf, line, inst, OFF_TRAJECTORY)
        if r is not None and _is_scored(r):
            x, y, xerr, yerr = _xyerr(r)
            if np.isfinite(x) and np.isfinite(y):
                ax.errorbar(
                    x, y, xerr=xerr, yerr=yerr, fmt="s", color=STATE_COLOR[OFF_TRAJECTORY],
                    ecolor=STATE_COLOR[OFF_TRAJECTORY], elinewidth=1.0, capsize=2.5,
                    markersize=6, zorder=3, label=SHORT[OFF_TRAJECTORY],
                )
        ax.legend(fontsize=7, loc="best", frameon=False, ncol=2)
        ax.axhline(0.0, color="#bbbbbb", linewidth=0.6, zorder=0)
        ax.axvline(0.0, color="#bbbbbb", linewidth=0.6, zorder=0)
    fig.tight_layout(rect=[0, 0, 1, 0.93])
    out = PLANE_FIG / f"plane_{line}.png"
    fig.savefig(out, dpi=160)
    plt.close(fig)
    log(f"[t3] wrote {out}")
    return out


def run_task3(log=None):
    if not PREREG_TASK3_FLAG.exists():
        raise StopStep("prereg", "PREREG_TASK3.flag missing")
    close_log = False
    if log is None:
        log = Logger(PLANE_DIR / "t3_report.txt")
        close_log = True
    plane_log_banner(log, "TASK3")
    log(PREREG_TASK3)
    med_p = PLANE_DIR / "t1_state_medians.csv"
    env_p = PLANE_DIR / "t1_null_envelope.csv"
    if not med_p.exists():
        raise StopStep("task3", f"missing {med_p}")
    if not env_p.exists():
        raise StopStep("task3", f"missing {env_p}")
    mdf = pd.read_csv(med_p)
    env = pd.read_csv(env_p)
    paths = {}
    for line in (AGED_LINE, YOUNG_LINE):
        paths[line] = str(plot_donor(line, mdf, env, log))
    rec = dict(
        figures=paths,
        identity=IDENTITY_KEY,
        instruments=list(AGE_INSTRUMENTS),
        trajectory=list(TRAJECTORY),
        NonReprog=OFF_TRAJECTORY,
        xlim="shared identity axis within donor",
        ylim="per instrument",
    )
    dump_json(PLANE_DIR / "t3_figures.json", jsonable(rec))
    man = load_manifest()
    man["status"] = "TASK3_DONE"
    save_manifest(man)
    progress_snapshot("write FINDINGS_PLANE.md from disk", stop="TASK3_DONE")
    if close_log:
        log.close()
    return rec


if __name__ == "__main__":
    try:
        run_task3()
    except StopStep as e:
        from plane_common import record_failure  # noqa: E402
        record_failure(e.step, e.message, e.details)
        progress_snapshot(f"fix STOP [{e.step}]", stop=f"STOP {e.step}: {e.message}")
        raise
