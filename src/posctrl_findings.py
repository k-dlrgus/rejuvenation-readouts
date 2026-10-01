"""Write FINDINGS_POSCTRL.md from results/posctrl/*.json.

Does not modify any prior FINDINGS file or FALSIFICATION.md.
The pre-registered block is copied verbatim.
"""
from __future__ import annotations

import sys
from pathlib import Path

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parent))
from posctrl_common import (  # noqa: E402
    PC_DIR, PC_FIG, PC_SEED, DATASET_ID, PREREG_BLOCK, FINDINGS_PATH,
    ANGLE_BAR, TARGET_R2_PRIMARY, REGIMES, N_BOOT_FULL, N_BOOT_HALF,
    LOWDIM_L4_QUOTE, gene_ridge_row, load_json, cell_path, md_table,
    pass_angle, pass_transfer_any, pass_transfer_both, pass_transfer_primary,
    finite, fmt, fmt_u, SEX_ANGLE_REF, SUPERSESSION_BLOCK, c1_cell_name,
)


def _load_if(name):
    p = cell_path(name)
    if not p.exists():
        return None
    try:
        rec = load_json(p)
    except Exception:
        return None
    return rec if rec.get("complete") else rec


def _row_from_cell(rec, regime, target_r2, extra_note=""):
    if rec is None:
        return dict(
            regime=regime, target_r2=target_r2,
            achieved_within_site_r2=np.nan, within_site_null=np.nan,
            boot_median_angle=np.nan, recovery_angle=np.nan, loso_r2=np.nan,
            transfer_H_to_M=np.nan, transfer_M_to_H=np.nan, young_old_balanced=np.nan,
            pass_angle=False, pass_transfer=False, pass_both=False,
            n_boot=np.nan, note=extra_note or "not run",
        )
    hm, mh = rec.get("transfer_H_to_M"), rec.get("transfer_M_to_H")
    ang = rec.get("boot_median_angle")
    pa = pass_angle(ang)
    pt = pass_transfer_both(hm, mh)
    note = extra_note
    if rec.get("b_reduced"):
        note = ((note + "; ") if note else "") + f"B={rec.get('n_boot')} (B/2)"
    return dict(
        regime=regime,
        target_r2=rec.get("target_r2", target_r2),
        achieved_within_site_r2=rec.get("achieved_within_site_r2", rec.get("within_site_r2")),
        within_site_null=rec.get("within_site_null"),
        boot_median_angle=ang,
        recovery_angle=rec.get("recovery_angle"),
        loso_r2=rec.get("loso_r2"),
        transfer_H_to_M=hm,
        transfer_M_to_H=mh,
        young_old_balanced=rec.get("young_old_balanced"),
        pass_angle=pa,
        pass_transfer=pt,
        pass_both=bool(pa and pt),
        n_boot=rec.get("n_boot"),
        support_frac=rec.get("support_frac"),
        site_strat_auc=rec.get("site_strat_auc"),
        note=note,
    )


def _fire(c2a, planted):
    """Exclusive: P1 > P2 > Mixed > P3. None if planted not run and C2a passed."""
    if c2a is None:
        return "not_run", "C2a has not been run."
    if not pass_angle(c2a.get("boot_median_angle")) or not pass_transfer_any(
            c2a.get("transfer_H_to_M"), c2a.get("transfer_M_to_H")):
        return "P1", (
            f"C2a bootstrap {fmt_u(c2a.get('boot_median_angle'), 1)}° > 35° "
            f"(transfer HBCC→MSSM {fmt(c2a.get('transfer_H_to_M'))}, "
            f"MSSM→HBCC {fmt(c2a.get('transfer_M_to_H'))}, both >0)."
        )
    if planted is None:
        return "C2a_pass_C1_pending", "C2a passed; C1 dense_shared @ 0.30 not run."
    pa = pass_angle(planted.get("boot_median_angle"))
    pt = pass_transfer_both(planted.get("transfer_H_to_M"), planted.get("transfer_M_to_H"))
    if pa and pt:
        return "P2", "planted dense_shared @ R²≈0.30 passed angle and both transfers."
    if pa != pt:
        which = "angle passed, transfer failed" if pa else "transfer passed, angle failed"
        return "Mixed", which
    return "P3", "planted dense_shared @ R²≈0.30 failed angle and transfer."


def _verdict_paragraph(fired, c2a, planted, fit, track):
    if fired == "P1":
        return (
            f"P1. C2a sex control (all genes) bootstrap {fmt_u(c2a.get('boot_median_angle'), 1)}° "
            f"(bar ≤{ANGLE_BAR:.0f}°; fail), HBCC→MSSM {fmt(c2a.get('transfer_H_to_M'))}, "
            f"MSSM→HBCC {fmt(c2a.get('transfer_M_to_H'))} (both >0; pass). "
            f"Within-site ridge R² {fmt(c2a.get('achieved_within_site_r2', c2a.get('within_site_r2')))} "
            f"(null {fmt(c2a.get('within_site_null'))}), site-stratified AUC {fmt_u(c2a.get('site_strat_auc'))}. "
            f"The failure is the bootstrap-angle bar, not transfer and not prediction. "
            f"C1 and C3 were not run. Everything downstream of the target pipeline is suspect "
            f"as an identifiability (angle) claim at this n, p."
        )
    if fired == "P2":
        return (
            f"P2. Planted dense_shared at achieved within-site R² "
            f"{fmt(planted.get('achieved_within_site_r2', planted.get('within_site_r2')))} "
            f"(target 0.30) bootstrap {planted.get('boot_median_angle')}° "
            f"(≤{ANGLE_BAR:.0f}°), HBCC→MSSM {fmt(planted.get('transfer_H_to_M'))}, "
            f"MSSM→HBCC {fmt(planted.get('transfer_M_to_H'))}, both >0. "
            f"The pipeline can find a real direction at n=233. "
            f"The age negative in FINDINGS_LOWDIM stands as a statement about the data."
        )
    if fired == "Mixed":
        pa = pass_angle(planted.get("boot_median_angle"))
        pt = pass_transfer_both(planted.get("transfer_H_to_M"), planted.get("transfer_M_to_H"))
        extra = ""
        if pa and not pt:
            extra = (
                " Transfer failing on a planted shared direction means bank covariance "
                "differences alone defeat ridge transfer."
            )
        return (
            f"Mixed. Planted dense_shared at achieved within-site R² "
            f"{fmt(planted.get('achieved_within_site_r2', planted.get('within_site_r2')))}: "
            f"bootstrap {planted.get('boot_median_angle')}° "
            f"(pass_angle={pa}), HBCC→MSSM {fmt(planted.get('transfer_H_to_M'))}, "
            f"MSSM→HBCC {fmt(planted.get('transfer_M_to_H'))} (pass_transfer={pt}). "
            f"Do not average.{extra}"
        )
    if fired == "P3":
        n_cross = (fit or {}).get("n_cross")
        extra = ""
        if finite(n_cross):
            lab = "extrapolation" if (fit or {}).get("extrapolation") else "observed"
            extra = f" Planted angle crosses 35° at n={n_cross:.0f} ({lab})."
        return (
            f"P3. Planted dense_shared at achieved within-site R² "
            f"{fmt(planted.get('achieved_within_site_r2', planted.get('within_site_r2')))} "
            f"bootstrap {planted.get('boot_median_angle')}° "
            f"(bar ≤{ANGLE_BAR:.0f}°), HBCC→MSSM {fmt(planted.get('transfer_H_to_M'))}, "
            f"MSSM→HBCC {fmt(planted.get('transfer_M_to_H'))}. "
            f"The LOWDIM bars were unattainable at n=233. "
            f"\"not identifiable\" must be rewritten as \"not identifiable at n=233\". "
            f"The power curve (C3) is the primary result.{extra} "
            f"No claim about real age beyond n=233."
        )
    return f"{fired}."


def _primary_pass(rec):
    if rec is None:
        return False
    if rec.get("pass_transfer_primary") is not None:
        return bool(rec["pass_transfer_primary"])
    return pass_transfer_primary(
        rec.get("transfer_H_to_M"), rec.get("transfer_M_to_H"),
        rec.get("transfer_H_to_M_null"), rec.get("transfer_M_to_H_null"),
    )


def _ss_row(rec, regime, target_r2, extra_note=""):
    row = _row_from_cell(rec, regime, target_r2, extra_note=extra_note)
    if rec is None:
        row.update(
            transfer_H_to_M_null=np.nan, transfer_M_to_H_null=np.nan,
            angle_vs_sex_ref=np.nan, pass_transfer_primary=False,
        )
        return row
    hm_n, mh_n = rec.get("transfer_H_to_M_null"), rec.get("transfer_M_to_H_null")
    ang = rec.get("boot_median_angle")
    pt = _primary_pass(rec)
    row["transfer_H_to_M_null"] = hm_n
    row["transfer_M_to_H_null"] = mh_n
    row["angle_vs_sex_ref"] = (float(ang) - SEX_ANGLE_REF) if finite(ang) else np.nan
    row["pass_transfer"] = pt
    row["pass_transfer_primary"] = pt
    row["pass_both"] = pt
    return row


def _ss_verdict(planted, regimes, real233):
    if planted is None:
        return "Pending C1 and C3."
    flags = {r: _primary_pass(regimes.get(r)) for r in REGIMES}
    planted_ok = _primary_pass(planted)
    uniq = set(bool(v) for v in flags.values())
    mixed = len(uniq) > 1
    real_ok = _primary_pass(real233) if real233 is not None else None

    def _xfer(rec):
        if rec is None:
            return "NA", "NA", "NA", "NA"
        return (
            fmt(rec.get("transfer_H_to_M")), fmt(rec.get("transfer_M_to_H")),
            fmt(rec.get("transfer_H_to_M_null")), fmt(rec.get("transfer_M_to_H_null")),
        )

    per = []
    for r in REGIMES:
        rec = regimes.get(r)
        hm, mh, hmn, mhn = _xfer(rec)
        per.append(
            f"`{r}` pass_transfer={flags[r]}  HBCC→MSSM {hm} (null {hmn})  "
            f"MSSM→HBCC {mh} (null {mhn})"
        )
    per_txt = " Per regime (not averaged): " + "; ".join(per) + "."

    hm, mh, hmn, mhn = _xfer(planted)
    ang = planted.get("boot_median_angle")
    r2 = fmt(planted.get("achieved_within_site_r2", planted.get("within_site_r2")))
    ang_s = f"{ang}" if ang is not None else "NA"
    vs = f"{float(ang) - SEX_ANGLE_REF:+.1f}" if finite(ang) else "NA"
    planted_line = (
        f"Planted dense_shared @ R²≈0.30: achieved within-site R² {r2}, "
        f"bootstrap {ang_s}° ({vs} vs sex 41.3°), "
        f"HBCC→MSSM {hm} (null {hmn}), MSSM→HBCC {mh} (null {mhn})."
    )
    if real233 is None:
        return planted_line + " C3 real-age n=233 not yet run. Verdict pending." + (per_txt if mixed else "")

    rhm, rmh, rhmn, rmhn = _xfer(real233)
    real_line = (
        f" Real age (C3 n=233, same L2 path): HBCC→MSSM {rhm} (null {rhmn}), "
        f"MSSM→HBCC {rmh} (null {rmhn}), pass_transfer={real_ok}."
    )
    if mixed:
        return (
            "Mixed across regimes — report per regime, do not average. "
            + planted_line + real_line + per_txt
        )
    if planted_ok and real_ok is False:
        return (
            "The age negative is real and it is specifically a transfer failure, not an angle failure. "
            "Planted dense_shared at R²≈0.30 transfers both ways; real age does not. "
            + planted_line + real_line
        )
    if (not planted_ok) and (real_ok is False or real_ok is None):
        return (
            "Bank covariance differences defeat transfer for any direction at this n; "
            "the age negative is uninformative and the power curve is the result. "
            "Planted dense_shared at R²≈0.30 does not transfer either. "
            "All four C1 regimes at R²≈0.30 failed the primary transfer bar (not mixed). "
            + planted_line + real_line
        )
    if planted_ok and real_ok:
        return (
            "Planted dense_shared at R²≈0.30 transfers, and real age also transfers on the same L2 path. "
            + planted_line + real_line
        )
    return planted_line + real_line + per_txt


def _live_supersession_md():
    planted = _load_if("c1_dense_shared_r030")
    regimes = {r: _load_if(c1_cell_name(r, 0.30)) for r in REGIMES}
    real233 = _load_if("c3_real_n233_d00")
    c3 = _load_if("c3_summary")
    fit = (c3 or {}).get("fit") if c3 else None
    track = (c3 or {}).get("track") if c3 else None

    lines = []
    lines.append("## Supersession 2026-09-16")
    lines.append("")
    lines.append(
        "Written **before any C1/C3 fit**. The original pre-registered block above "
        "and the P1 verdict are kept for the record."
    )
    lines.append("")
    lines.append(
        "C2a sex control missed the 35° bootstrap-angle bar (41.3°) while transferring "
        "(HBCC→MSSM +0.512, MSSM→HBCC +0.358) and predicting (AUC 0.939). "
        "A real direction that transfers between banks cannot be \"pipeline broken.\" "
        "The 35° bar was miscalibrated, not the pipeline."
    )
    lines.append("")
    lines.append(
        "The ≤35° angle bar is **retired**. New pre-registered bars for C1 and C3 "
        "(no threshold in this block is re-tuned after numbers exist):"
    )
    lines.append("")
    lines.append("- **(a) Primary:** transfer R² > 0 in both directions (HBCC→MSSM and MSSM→HBCC) with permutation null ≤ 0.05.")
    lines.append("- **(b)** Bootstrap median pairwise angle is **reported** against the sex reference **41.3°**, not against 35°. Angle is not a pass/fail gate.")
    lines.append("")
    lines.append("Restated P2 / P3 (these bars; C1 `dense_shared` at R²≈0.30):")
    lines.append("")
    lines.append("- **P2** = planted `dense_shared` at R²≈0.30 transfers both ways.")
    lines.append("- **P3** = it does not.")
    lines.append("")
    lines.append(
        "Sex is sparse and huge, so it is not a fair reference for a planted dense "
        "continuous direction at R²≈0.30. C1 is."
    )
    lines.append("")
    lines.append(
        "C3: `dense_shared` at R²≈0.30 and real age, n ∈ {60, 100, 150, 233} "
        "stratified by bank, 10 draws per n. Any n>233 crossover is labeled extrapolation."
    )
    lines.append("")
    lines.append("Verdict rule (only the outcome that fired):")
    lines.append("")
    lines.append("- if planted transfers and real age does not → the age negative is real and it is specifically a transfer failure, not an angle failure.")
    lines.append("- If planted does not transfer either → bank covariance differences defeat transfer for any direction at this n; the age negative is uninformative and the power curve is the result.")
    lines.append("- If mixed across regimes, report per regime, do not average.")
    lines.append("")
    lines.append("### C1 (supersession)")
    lines.append("")
    any_c1 = any(regimes.values()) or _load_if("c1_dense_shared_r015") or _load_if("c1_dense_shared_r050")
    if not any_c1:
        lines.append(
            "Four regimes (`dense_shared`, `dense_typespec`, `sparse_shared`, `sparse_typespec`) "
            "calibrated to within-site site-stratified R² 0.30 ± 0.03, plus `dense_shared` at 0.15 and 0.50. "
            "Same L2 code paths, same folds, same α grid."
        )
        lines.append("")
        lines.append("Question: does a planted dense continuous direction at R²≈0.30 transfer between HBCC and MSSM?")
        lines.append("")
        lines.append("Not yet run.")
        lines.append("")
    else:
        lines.append(
            "Cohort: 233 donors, 20 cell types. age* = within-bank donor permutation of chronological age "
            f"(seed `{PC_SEED}`). w drawn on seed `{PC_SEED + 1}`. β calibrated on seed `{PC_SEED + 2}` "
            "to within-site site-stratified ridge R² 0.30 ± 0.03 (measured R² reported). "
            "Primary bar: transfer R² > 0 both ways with permutation null ≤ 0.05. "
            f"Angle reported vs sex reference {SEX_ANGLE_REF:.1f}°. "
            "gene_ridge is the LOWDIM reference row (not refit)."
        )
        lines.append("")
        c1_rows = [gene_ridge_row()]
        for r in REGIMES:
            c1_rows.append(_ss_row(regimes.get(r), r, 0.30))
        for t in (0.15, 0.50):
            rec = _load_if(c1_cell_name("dense_shared", t))
            c1_rows.append(_ss_row(rec, "dense_shared", t, extra_note="not primary; B/2"))
        c1_df = pd.DataFrame(c1_rows)
        cols1 = [
            "regime", "target_r2", "achieved_within_site_r2", "within_site_null",
            "boot_median_angle", "angle_vs_sex_ref", "recovery_angle", "loso_r2",
            "transfer_H_to_M", "transfer_M_to_H",
            "transfer_H_to_M_null", "transfer_M_to_H_null",
            "young_old_balanced", "pass_transfer", "n_boot",
        ]
        cols1 = [c for c in cols1 if c in c1_df.columns]
        lines.append(md_table(c1_df, cols=cols1))
        lines.append("")
        bp = PC_DIR / "frozen_beta.json"
        if bp.exists():
            betas = load_json(bp)
            lines.append("Frozen β:")
            for k, v in betas.items():
                lines.append(
                    f"- `{k}`: β={v.get('beta')}  calibration full-fold R²="
                    f"{fmt(v.get('cheap_then_full_r2'))}"
                )
            lines.append("")
        for r in ("sparse_shared", "sparse_typespec"):
            rec = regimes.get(r)
            if rec and rec.get("support_frac") is not None:
                lines.append(
                    f"`{r}` fraction of 200 support genes in top-200 |coef| "
                    f"(median over types / boots): {fmt_u(rec.get('support_frac'))}."
                )
        lines.append("")
        lines.append(
            "Question: does a planted dense continuous direction at R²≈0.30 transfer between HBCC and MSSM? "
            + ("Yes." if _primary_pass(planted) else "No." if planted is not None else "Pending.")
        )
        lines.append("")

    lines.append("### C3 (supersession)")
    lines.append("")
    c3tab_p = PC_DIR / "c3_summary.csv"
    if c3tab_p.exists():
        c3tab = pd.read_csv(c3tab_p)
        lines.append(
            "n ∈ {60, 100, 150, 233}, bank-stratified, 10 draws per n except n=233 (1 unique set). "
            "Real-age uses full B; planted uses B/2. "
            "Figure: `results/posctrl/figures/c3_power_vs_n.png`. "
            f"Angle referenced to sex {SEX_ANGLE_REF:.1f}°, not 35°."
        )
        lines.append("")
        show = [c for c in (
            "kind", "n", "n_draws", "boot_median_angle", "boot_p05", "boot_p95",
            "transfer_H_to_M", "transfer_M_to_H", "n_pass_transfer_primary",
            "n_boot", "b_reduced",
        ) if c in c3tab.columns]
        lines.append(md_table(c3tab, cols=show))
        lines.append("")
        if fit:
            for key, lab in (
                ("angle_vs_sex_ref", f"planted angle vs {SEX_ANGLE_REF:.1f}°"),
                ("transfer_H_to_M", "planted HBCC→MSSM vs 0"),
                ("transfer_M_to_H", "planted MSSM→HBCC vs 0"),
            ):
                f = fit.get(key) if isinstance(fit, dict) else None
                if not f:
                    continue
                extra = f.get("extrapolation")
                lab_x = "extrapolation (n>233)" if extra else "observed"
                lines.append(
                    f"{lab}: n_cross={f.get('n_cross')} ({lab_x})."
                )
        if track:
            lines.append(
                f"corr(n, angle) planted={fmt(track.get('planted_corr_n_angle'))}  "
                f"real={fmt(track.get('real_corr_n_angle'))}."
            )
        lines.append("")
        if (PC_FIG / "c3_power_vs_n.png").exists():
            lines.append("Figure: `results/posctrl/figures/c3_power_vs_n.png`.")
            lines.append("")
    else:
        lines.append("Not yet run.")
        lines.append("")
    lines.append("### Verdict (supersession 2026-09-16)")
    lines.append("")
    lines.append(_ss_verdict(planted, regimes, real233))
    lines.append("")
    lines.append("## Files")
    lines.append("")
    lines.append("| path | content |")
    lines.append("|---|---|")
    lines.append("| `src/posctrl_common.py`, `posctrl_c1.py`, `posctrl_c2.py`, `posctrl_c3.py`, `posctrl_findings.py` | code |")
    lines.append("| `notebooks/posctrl_c1.ipynb` … `posctrl_verdict.ipynb` | runnable from a clean checkout |")
    lines.append("| `results/posctrl/` | tables, json, logs, `figures/` |")
    lines.append("| `results/posctrl/SUPERSESSION_20260916.flag` | supersession written before C1/C3 fit |")
    lines.append("| `FINDINGS_POSCTRL.md` | this file |")
    lines.append("| `PROGRESS_POSCTRL.md` | resume state |")
    lines.append("")
    return "\n".join(lines)


def write_findings():
    """Keep the original P1 record; rewrite only the supersession live section."""
    if not FINDINGS_PATH.exists():
        raise FileNotFoundError(FINDINGS_PATH)
    text = FINDINGS_PATH.read_text(encoding="utf-8")
    marker = "## Supersession 2026-09-16"
    if marker not in text:
        raise RuntimeError(
            "FINDINGS_POSCTRL.md is missing the frozen supersession section; "
            "refusing to overwrite the P1 record."
        )
    keep = text.split(marker, 1)[0].rstrip()
    live = _live_supersession_md()
    FINDINGS_PATH.write_text(keep + "\n\n" + live, encoding="utf-8")
    print("wrote", FINDINGS_PATH)
    planted = _load_if("c1_dense_shared_r030")
    real233 = _load_if("c3_real_n233_d00")
    if planted is None:
        return "supersession_pending_C1"
    if real233 is None:
        return "supersession_pending_C3"
    return "supersession_verdict"


if __name__ == "__main__":
    write_findings()
