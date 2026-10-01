"""C3 — power curve.

dense_shared at achieved R²≈0.30 and the real age target, subsample donors
n ∈ {60, 100, 150, 233} stratified by bank, 10 draws per n (seed 20260914+3).
Each draw: bootstrap median angle, HBCC→MSSM, MSSM→HBCC.

Planted C3 uses B/2. Real-age C3 uses full B. n=233 is the full cohort (1 unique draw).

Usage: python src/posctrl_c3.py
       python src/posctrl_c3.py real
       python src/posctrl_c3.py planted
"""
from __future__ import annotations

import sys
from pathlib import Path

import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

sys.path.insert(0, str(Path(__file__).resolve().parent))
from posctrl_common import (  # noqa: E402
    PC_DIR, PC_FIG, PC_SEED, C3_SEED, DATASET_ID, SEX_ANGLE_REF, TRANSFER_NULL_BAR,
    N_BOOT_FULL, N_BOOT_HALF, N_PERM_FULL, N_PERM_HALF,
    C3_NS, N_C3_DRAWS, TARGET_R2_PRIMARY,
    Logger, pc_log_banner, load_phase1_matrix, load_json, dump_json,
    plant, run_l2_gene_battery, subsample_donors, restrict,
    save_cell, load_cell, cell_done, cell_path, jsonable, finite,
    p1_blocks_c1_c3, pass_transfer_primary, refresh_progress,
)


def n_cross_threshold(ns, values, threshold, decreasing=True):
    """Fit y = a + b/sqrt(n) (decreasing) or y = a + b*sqrt(n) (increasing).

    n where the curve equals `threshold`. Extrapolation if n>233.
    """
    ns = np.asarray(ns, float)
    val = np.asarray(values, float)
    m = np.isfinite(ns) & np.isfinite(val)
    ns, val = ns[m], val[m]
    form = "a+b/sqrt(n)" if decreasing else "a+b*sqrt(n)"
    empty = dict(n_cross=np.nan, a=np.nan, b=np.nan, form=form, extrapolation=None,
                 threshold=threshold, n_cross_observed=np.nan, n_cross_fitted=np.nan)
    if len(ns) < 2:
        return empty
    if decreasing:
        x = 1.0 / np.sqrt(ns)
        A = np.column_stack([np.ones(len(x)), x])
        coef, *_ = np.linalg.lstsq(A, val, rcond=None)
        a, b = float(coef[0]), float(coef[1])
        if b < 0:
            b = 0.0
            a = float(np.mean(val))
        below = ns[val <= threshold]
        obs_n = float(np.min(below)) if len(below) else np.nan
        if b < 1e-12:
            n_hat = np.nan if a > threshold else (float(np.min(ns)) if a <= threshold else np.nan)
        else:
            n_hat = np.nan if a >= threshold else float((b / (threshold - a)) ** 2)
    else:
        x = np.sqrt(ns)
        A = np.column_stack([np.ones(len(x)), x])
        coef, *_ = np.linalg.lstsq(A, val, rcond=None)
        a, b = float(coef[0]), float(coef[1])
        above = ns[val > threshold]
        obs_n = float(np.min(above)) if len(above) else np.nan
        if abs(b) < 1e-12:
            n_hat = float(np.min(ns)) if a > threshold else np.nan
        else:
            # a + b*sqrt(n) = threshold  => n = ((threshold-a)/b)^2 if b>0 and threshold>a
            n_hat = float(((threshold - a) / b) ** 2) if b > 0 and threshold >= a else np.nan
    if finite(obs_n):
        n_cross, extra = obs_n, False
    else:
        n_cross, extra = n_hat, bool(finite(n_hat) and n_hat > 233)
    return dict(
        n_cross=n_cross, a=a, b=b, form=form, threshold=threshold,
        extrapolation=extra, n_cross_observed=obs_n, n_cross_fitted=n_hat,
    )


def n_cross_35(ns, med_angles):
    """Kept for the P1 record; C3 reporting uses n_cross_threshold vs 41.3° / transfer 0."""
    return n_cross_threshold(ns, med_angles, SEX_ANGLE_REF, decreasing=True)


def _draw_cell_name(kind, n, draw):
    return f"c3_{kind}_n{int(n):03d}_d{int(draw):02d}"


def _run_draw(kind, n, draw, Ys, obs_s, ys, types, pack, n_boot, n_perm, log):
    name = _draw_cell_name(kind, n, draw)
    if cell_done(name):
        return load_cell(name)
    log(f"[C3] {name} n_donors={obs_s.donor.nunique()}  "
        f"H={int((obs_s.groupby('donor').Source.first()=='H').sum())} "
        f"M={int((obs_s.groupby('donor').Source.first()=='M').sum())}")
    rng = np.random.default_rng(PC_SEED)
    bat = run_l2_gene_battery(
        Ys, obs_s, ys, rng, log,
        n_boot=n_boot, n_perm=n_perm,
        planted_pack=pack if kind == "planted" else None,
        parts=("boot", "transfer"),
        y_is_binary=False, split_on_chronological=False, tag=name,
    )
    rec = dict(
        cell=name, kind=kind, n_requested=int(n), n_donors=int(obs_s.donor.nunique()),
        draw=int(draw), n_boot=n_boot, n_perm=n_perm, b_reduced=bool(n_boot < N_BOOT_FULL),
        n_picked=int(obs_s.donor.nunique()),
        boot_median_angle=bat.get("boot_median_angle"),
        transfer_H_to_M=bat.get("transfer_H_to_M"),
        transfer_M_to_H=bat.get("transfer_M_to_H"),
        transfer_H_to_M_null=bat.get("transfer_H_to_M_null"),
        transfer_M_to_H_null=bat.get("transfer_M_to_H_null"),
        recovery_angle=bat.get("recovery_angle"),
        elapsed_s=bat.get("elapsed_s"),
        sex_angle_ref=SEX_ANGLE_REF,
        path=str(cell_path(name)),
    )
    rec["pass_transfer"] = pass_transfer_primary(
        rec.get("transfer_H_to_M"), rec.get("transfer_M_to_H"),
        rec.get("transfer_H_to_M_null"), rec.get("transfer_M_to_H_null"),
    )
    rec["pass_transfer_primary"] = rec["pass_transfer"]
    if finite(rec.get("boot_median_angle")):
        rec["angle_vs_sex_ref"] = float(rec["boot_median_angle"]) - SEX_ANGLE_REF
    save_cell(name, jsonable(rec))
    refresh_progress()
    return rec


def _figure(tab):
    """One figure: bootstrap angle and both transfers vs n, with draw spread."""
    fig, axes = plt.subplots(3, 1, figsize=(7.2, 9.2), sharex=True)
    kinds = (
        ("planted", "tab:blue", "planted dense_shared R²≈0.30"),
        ("real", "tab:red", "real age"),
    )

    ax = axes[0]
    for kind, col, lab in kinds:
        sub = tab[tab.kind == kind]
        if not len(sub):
            continue
        g = sub.groupby("n_requested")
        ns = np.array(sorted(g.groups))
        med = np.array([g.get_group(n).boot_median_angle.median() for n in ns])
        ax.plot(ns, med, "-o", color=col, label=lab, lw=1.8)
        for n in ns:
            ys = g.get_group(n).boot_median_angle.to_numpy(float)
            if np.isfinite(ys).any():
                ax.vlines(n, np.nanmin(ys), np.nanmax(ys), color=col, lw=0.8, alpha=0.5)
                ax.scatter(np.full(len(ys), n, float), ys, color=col, s=18, alpha=0.45, zorder=3)
    ax.axhline(SEX_ANGLE_REF, c="k", ls="--", lw=1, label=f"sex reference {SEX_ANGLE_REF:.1f}°")
    ax.set_ylabel("bootstrap median pairwise angle (deg)")
    ax.set_title("C3 power curve")
    ax.legend(fontsize=7)
    ax.set_xticks(list(C3_NS))

    for ax, coln, title in zip(
        axes[1:],
        ("transfer_H_to_M", "transfer_M_to_H"),
        ("HBCC → MSSM", "MSSM → HBCC"),
    ):
        for kind, col, lab in kinds:
            sub = tab[tab.kind == kind]
            if not len(sub):
                continue
            g = sub.groupby("n_requested")
            ns = np.array(sorted(g.groups))
            med = np.array([g.get_group(n)[coln].median() for n in ns])
            ax.plot(ns, med, "-o", color=col, label=lab)
            for n in ns:
                ys = g.get_group(n)[coln].to_numpy(float)
                if np.isfinite(ys).any():
                    ax.vlines(n, np.nanmin(ys), np.nanmax(ys), color=col, lw=0.8, alpha=0.5)
                    ax.scatter(np.full(len(ys), n, float), ys, color=col, s=18, alpha=0.45, zorder=3)
        ax.axhline(0, c="0.6", lw=0.8)
        ax.set_title(title)
        ax.set_ylabel("median-over-types R²")
        ax.legend(fontsize=7)
        ax.set_xticks(list(C3_NS))
    axes[-1].set_xlabel("n donors (bank-stratified subsample)")
    fig.tight_layout()
    fig.savefig(PC_FIG / "c3_power_vs_n.png", dpi=140)
    plt.close(fig)


def _summarize(tab, log):
    rows = []
    for kind in ("planted", "real"):
        sub = tab[tab.kind == kind]
        for n in C3_NS:
            s = sub[sub.n_requested == n]
            if not len(s):
                continue
            rows.append(dict(
                kind=kind, n=int(n), n_draws=int(len(s)),
                boot_median_angle=float(s.boot_median_angle.median()),
                boot_p05=float(np.nanpercentile(s.boot_median_angle, 5)),
                boot_p95=float(np.nanpercentile(s.boot_median_angle, 95)),
                transfer_H_to_M=float(s.transfer_H_to_M.median()),
                transfer_M_to_H=float(s.transfer_M_to_H.median()),
                n_pass_transfer_primary=int(s["pass_transfer"].fillna(False).sum()) if "pass_transfer" in s.columns else np.nan,
                n_boot=int(s.n_boot.iloc[0]),
                b_reduced=bool(s.b_reduced.iloc[0]),
            ))
    summ = pd.DataFrame(rows)
    summ.to_csv(PC_DIR / "c3_summary.csv", index=False)
    planted = tab[tab.kind == "planted"]
    real = tab[tab.kind == "real"]
    fit = {}
    if len(planted):
        g_ang = planted.groupby("n_requested").boot_median_angle.median()
        g_hm = planted.groupby("n_requested").transfer_H_to_M.median()
        g_mh = planted.groupby("n_requested").transfer_M_to_H.median()
        fit_ang = n_cross_threshold(g_ang.index.to_numpy(), g_ang.to_numpy(), SEX_ANGLE_REF, decreasing=True)
        fit_hm = n_cross_threshold(g_hm.index.to_numpy(), g_hm.to_numpy(), 0.0, decreasing=False)
        fit_mh = n_cross_threshold(g_mh.index.to_numpy(), g_mh.to_numpy(), 0.0, decreasing=False)
        fit = dict(angle_vs_sex_ref=fit_ang, transfer_H_to_M=fit_hm, transfer_M_to_H=fit_mh)
        for lab, f in (("angle vs 41.3°", fit_ang), ("H→M vs 0", fit_hm), ("M→H vs 0", fit_mh)):
            extra = f.get("extrapolation")
            log(f"[C3] planted {lab}: n_cross={f.get('n_cross')}  "
                f"extrapolation={extra}"
                + ("  (n>233 labeled extrapolation)" if extra else ""))
    # track vs planted
    track = None
    if len(planted) and len(real):
        notes = []
        for n in C3_NS:
            p = planted[planted.n_requested == n].boot_median_angle.to_numpy(float)
            r = real[real.n_requested == n].boot_median_angle.to_numpy(float)
            if len(p) < 1 or len(r) < 1:
                continue
            rmed, pmed = float(np.nanmedian(r)), float(np.nanmedian(p))
            plo, phi = float(np.nanmin(p)), float(np.nanmax(p))
            inside = plo <= rmed <= phi
            notes.append(dict(n=int(n), real_med=rmed, planted_med=pmed, real_inside_planted_spread=bool(inside)))
        # planted improves if Spearman n vs angle is negative
        def _slope(sub):
            g = sub.groupby("n_requested").boot_median_angle.median()
            if len(g) < 2:
                return np.nan
            return float(np.corrcoef(g.index.astype(float), g.to_numpy())[0, 1])
        track = dict(
            planted_corr_n_angle=_slope(planted),
            real_corr_n_angle=_slope(real),
            per_n=notes,
        )
        log(f"[C3] corr(n, angle) planted={track['planted_corr_n_angle']:+.3f}  "
            f"real={track['real_corr_n_angle']:+.3f}")
    dump_json(PC_DIR / "c3_fit.json", dict(planted_curve=fit, track=track))
    return summ, fit, track


def run(kind=None):
    log = Logger(PC_DIR / "c3_report.txt", mode="a" if (PC_DIR / "c3_report.txt").exists() else "w")
    try:
        pc_log_banner(log, "C3 — power curve")
        if p1_blocks_c1_c3():
            log("[C3] STOP P1 is set — refusing.")
            return dict(skipped=True, reason="STOP P1")
        log("[C3] P1 superseded 2026-09-16; 35° bar retired; running power curve.")
        data = load_phase1_matrix(log)
        Y, genes, obs = data["Y"], data["genes"], data["obs"]
        types = sorted(pd.unique(obs.celltype.astype(str)))
        kinds = [kind] if kind in ("real", "planted") else ["real", "planted"]
        # planted X
        y_real = obs.age.to_numpy(float)
        y_plant, pack, beta = None, None, None
        if "planted" in kinds:
            from posctrl_c1 import _ensure_age_star, _ensure_w, _ensure_beta
            y_plant = _ensure_age_star(obs, log)
            packs, _ = _ensure_w(Y, genes, obs, types, log)
            pack = packs["dense_shared"]
            beta, _ = _ensure_beta(Y, obs, y_plant, types, packs, "dense_shared", TARGET_R2_PRIMARY, log)
            Y_plant = plant(Y, y_plant, obs, types, pack, beta)
            log(f"[C3] planted dense_shared β={beta:.6g}")
        else:
            Y_plant = None
        rng_sub = np.random.default_rng(C3_SEED)
        recs = []
        n_don = int(obs.donor.nunique())
        for kd in ("real", "planted"):
            Y_use = Y if kd == "real" else Y_plant
            y_use = y_real if kd == "real" else y_plant
            pk = None if kd == "real" else pack
            n_boot = N_BOOT_FULL if kd == "real" else N_BOOT_HALF
            n_perm = N_PERM_FULL if kd == "real" else N_PERM_HALF
            for n in C3_NS:
                n_draws = 1 if int(n) >= n_don else N_C3_DRAWS
                for d in range(N_C3_DRAWS):
                    if int(n) >= n_don:
                        Ys, obs_s, ys = Y_use, obs, y_use
                    else:
                        row, _picked = subsample_donors(obs, int(n), rng_sub)
                        if kd not in kinds or d >= n_draws or Y_use is None:
                            continue
                        Ys, obs_s = restrict(Y_use, obs, row)
                        ys = np.asarray(y_use, float)[row]
                    if kd not in kinds or d >= n_draws or Y_use is None:
                        continue
                    rec = _run_draw(kd, n, d, Ys, obs_s, ys, types, pk, n_boot, n_perm, log)
                    recs.append(rec)
        tab = pd.DataFrame(recs)
        tab.to_csv(PC_DIR / "c3_draws.csv", index=False)
        summ, fit, track = _summarize(tab, log)
        _figure(tab)
        save_cell("c3_summary", jsonable(dict(
            n_draws_table=summ.to_dict(orient="records"),
            fit=fit, track=track, n_cells=int(len(tab)),
        )))
        log("[C3] wrote c3_summary.json and figures/c3_power_vs_n.png")
        refresh_progress()
        from posctrl_findings import write_findings
        write_findings()
        return dict(table=tab, summary=summ, fit=fit)
    finally:
        log.close()


if __name__ == "__main__":
    arg = sys.argv[1] if len(sys.argv) > 1 else None
    run(arg)
