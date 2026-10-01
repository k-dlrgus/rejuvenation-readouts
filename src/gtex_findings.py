"""Write FINDINGS_GTEX.md from results/gtex/. Does not modify other FINDINGS*.md."""
from __future__ import annotations

import json
import sys
from pathlib import Path

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parent))
from gtex_common import (  # noqa: E402
    GTEX_DIR, GTEX_SEED, RELEASE, RNASEQ_PIPELINE, DBGAP, TISSUES, PRIMARY_TISSUE,
    NEURONAL_MARKERS, GLIAL_MARKERS, CAVEATS, PREREG_BARS, PREREG_READING,
    PREREG_BARS_ORIGINAL, PREREG_READING_ORIGINAL, PREREG_DATE, PREREG_FLAG,
    LOWDIM_L4_QUOTE, SMAFRZE_RNASEQ, EXPECTED_AGE_BINS, FINDINGS_PATH,
    load_json, dump_json, stage1_gate, fmt, md_table, TISSUE_SLUG,
    N_PERM_GTEX, SEX_ANGLE_REF, PLANTED_ANGLE_N233,
)


def _signed(x, d=3):
    try:
        v = float(x)
    except (TypeError, ValueError):
        return "NA"
    if not np.isfinite(v):
        return "NA"
    return f"{v:+.{d}f}"


def _audit_table(audit: pd.DataFrame) -> str:
    cov = ["SMRIN", "SMTSISCH", "DTHHRDY", "SEX", "SMCENTER", "SMNABTCH", "JOINT"]
    lines = [
        "| tissue | " + " | ".join(cov) + " |",
        "|" + "|".join(["---"] * (1 + len(cov))) + "|",
    ]
    for t in TISSUES:
        d = audit[audit.tissue == t]
        cells = []
        for c in cov:
            row = d[d.covariate == c]
            if row.empty:
                cells.append("NA")
            else:
                r = row.iloc[0]
                extra = ""
                if c in ("SMNABTCH", "SMCENTER", "DTHHRDY", "SEX", "JOINT") and "n_levels" in r.index:
                    extra = f", L={int(r.n_levels)}"
                cells.append(f"{_signed(r.r2)} (n={int(r.n)}{extra})")
        lines.append(f"| {t} | " + " | ".join(cells) + " |")
    return "\n".join(lines)


def _cohort_table(summary: pd.DataFrame) -> str:
    cols = ["tissue", "n_samples", "n_donors", "n_SMCENTER", "n_SMNABTCH",
            "n_SMGEBTCH", "n_male", "n_female"]
    lines = [
        "| " + " | ".join(cols) + " |",
        "|" + "|".join(["---"] + ["---:"] * (len(cols) - 1)) + "|",
    ]
    for _, r in summary.iterrows():
        lines.append(
            f"| {r.tissue} | {int(r.n_samples)} | {int(r.n_donors)} | {int(r.n_SMCENTER)} | "
            f"{int(r.n_SMNABTCH)} | {int(r.n_SMGEBTCH)} | {int(r.n_male)} | {int(r.n_female)} |"
        )
    return "\n".join(lines)


def _age_table(summary: pd.DataFrame) -> str:
    lines = [
        "| tissue | " + " | ".join(EXPECTED_AGE_BINS) + " |",
        "|" + "|".join(["---"] * (1 + len(EXPECTED_AGE_BINS))) + "|",
    ]
    for _, r in summary.iterrows():
        cells = " | ".join(str(int(r[f"n_AGE_{b}"])) for b in EXPECTED_AGE_BINS)
        lines.append(f"| {r.tissue} | {cells} |")
    return "\n".join(lines)


def write_findings():
    man_path = GTEX_DIR / "manifest.json"
    if not man_path.exists():
        raise SystemExit("results/gtex/manifest.json missing — run Stage 0")
    man = load_json(man_path)
    gate = stage1_gate()
    man["stage1_gate"] = gate
    dump_json(man_path, man)
    if (GTEX_DIR / "stage0_summary.json").exists():
        summ = load_json(GTEX_DIR / "stage0_summary.json")
        summ["stage1_gate"] = gate
        dump_json(GTEX_DIR / "stage0_summary.json", summ)
    failures = man.get("failures") or []
    status_file = man.get("status", "UNKNOWN")

    if status_file == "STOP" or failures:
        status_line = (
            f"**Status:** STOP Stage 0. Seed `{GTEX_SEED}`. Release {RELEASE} "
            f"({RNASEQ_PIPELINE}; {DBGAP}). Failures recorded in `results/gtex/manifest.json`."
        )
    elif not gate.get("open"):
        status_line = (
            f"**Status:** Stage 0 done; Stage 1 gated. Seed `{GTEX_SEED}`. "
            f"Release {RELEASE} ({RNASEQ_PIPELINE}; {DBGAP}). "
            f"Stage 1 not run: {str(gate.get('reason')).rstrip('.')}."
        )
    elif not (GTEX_DIR / "stage1_summary.json").exists():
        status_line = (
            f"**Status:** Stage 0 done; Stage 1 open ({gate.get('status')}). "
            f"Pre-registration {PREREG_DATE} written **before any T-cell fit**. "
            f"Seed `{GTEX_SEED}`. Release {RELEASE} ({RNASEQ_PIPELINE}; {DBGAP})."
        )
    else:
        status_line = (
            f"**Status:** Stage 1 complete ({gate.get('status')}). "
            f"Seed `{GTEX_SEED}`. Release {RELEASE} ({RNASEQ_PIPELINE}; {DBGAP}). "
            f"Primary transfer metric: held-out Pearson r (n_perm=200)."
        )

    lines = [
        "# FINDINGS_GTEX — GTEx cortex: is the age direction stable when logistics are recorded?",
        "",
        status_line,
        "",
        "Does not modify `FINDINGS_LOWDIM.md`, `FINDINGS_GEOMETRY.md`, `FINDINGS_POSCTRL.md`, "
        "`FINDINGS_TARGET.md`, `FINDINGS_TRAJECTORY.md`, `FINDINGS_BRAIN_PHASE1.md`, or "
        "`FALSIFICATION.md`. No perturbation data. No TF Atlas. No candidate interventions.",
        "",
        f"Reproduced by `notebooks/gtex_stage0.ipynb`, `gtex_stage1.ipynb`. "
        f"RNA-seq analysis freeze: `SMAFRZE=={SMAFRZE_RNASEQ}`. Six tissues, no more.",
        "",
        "## Caveats (GTEx, before Stage 1)",
        "",
        CAVEATS,
        "",
    ]

    if failures:
        lines += ["## STOP", "", "Stage 0 halted. Not substituting columns or repairing rows.", ""]
        for f in failures:
            lines.append(f"- **{f.get('step')}:** {f.get('message')}")
        lines += ["", "See `results/gtex/manifest.json`.", ""]

    # Stage 0 tables if present
    cohort_csv = GTEX_DIR / "stage0_cohort.csv"
    audit_csv = GTEX_DIR / "stage0_covariate_audit.csv"
    overlap_csv = GTEX_DIR / "stage0_gene_overlap.csv"
    comp_csv = GTEX_DIR / "stage0_composition.csv"
    summary = load_json(GTEX_DIR / "stage0_summary.json") if (GTEX_DIR / "stage0_summary.json").exists() else {}

    if cohort_csv.exists():
        coh = pd.read_csv(cohort_csv)
        lines += ["## Stage 0 — cohort", "", _cohort_table(coh), "", "AGE bins (n samples):", "", _age_table(coh), ""]
        rin_lines = ["| tissue | SMRIN median (min–max, miss) | SMTSISCH median (min–max, miss) |",
                     "|---|---|---|"]
        for _, r in coh.iterrows():
            rin_lines.append(
                f"| {r.tissue} | {r.SMRIN_median:.2f} ({r.SMRIN_min:.2f}–{r.SMRIN_max:.2f}, "
                f"miss={int(r.SMRIN_n_missing)}) | {r.SMTSISCH_median:.1f} "
                f"({r.SMTSISCH_min:.1f}–{r.SMTSISCH_max:.1f}, miss={int(r.SMTSISCH_n_missing)}) |"
            )
        lines += rin_lines + [""]
        lines.append("Full SMCENTER / SMNABTCH / DTHHRDY / SEX counts: `results/gtex/stage0_cohort.md`.")
        lines.append("")

    flags = (summary.get("stop_flags") or man.get("stop_flags") or {})
    if flags:
        lines += [
            "## Stage 0 STOP checks",
            "",
            f"- BA9 n after freeze filter: **{flags.get('ba9_n')}** "
            f"(underpowered vs POSCTRL n=150 bar if n<150: {flags.get('ba9_n_underpowered')}).",
            f"- SMCENTER levels with ≥25 BA9 samples: **{flags.get('n_SMCENTER_ge25')}** "
            f"(T2 not testable if <3: {flags.get('t2_not_testable')}).",
            f"- BA9 SMCENTER counts: `{flags.get('center_counts')}`. "
            "Values are the file strings; BA9/BA24/hippocampus are comma-joined (`B1, A1`); they are not split.",
            "",
        ]

    if audit_csv.exists():
        audit = pd.read_csv(audit_csv)
        lines += [
            "## Stage 0 — covariate audit",
            "",
            "OLS R² of age-bin midpoint on each covariate, descriptive, no CV. "
            "Categorical covariates one-hot (drop first). Joint = all six together. "
            "n is complete cases for that regression; incomplete rows are not dropped from the cohort. "
            "L = number of levels (SMNABTCH L is large relative to n; that R² is a saturated-batch fit, not a logistics slope).",
            "",
            _audit_table(audit),
            "",
        ]

    if overlap_csv.exists():
        ov = pd.read_csv(overlap_csv)
        lines += ["## Stage 0 — gene filter and DLPFC overlap", "",
                  "| tissue | n_samples | n_genes_raw | n_genes_filter | n_overlap_dlpfc |",
                  "|---|---:|---:|---:|---:|"]
        for _, r in ov.iterrows():
            lines.append(
                f"| {r.tissue} | {int(r.n_samples)} | {int(r.n_genes_raw)} | "
                f"{int(r.n_genes_filter)} | {int(r.n_overlap_dlpfc)} |"
            )
        lines += ["",
                  "Filter (frozen): ≥6 counts in ≥20% of samples per tissue. "
                  "Then TMM (edgeR defaults) → log2-CPM (`prior.count=2`). "
                  "z-score per gene on training rows only inside every fold (Stage 1).",
                  ""]

    if comp_csv.exists():
        comp = pd.read_csv(comp_csv)
        lines += [
            "## Stage 0 — composition proxies (brain tissues)",
            "",
            f"Neuronal markers (frozen): {', '.join(NEURONAL_MARKERS)}. "
            f"Glial markers (frozen a priori): {', '.join(GLIAL_MARKERS)}. "
            "Score = sum of gene-z log-CPM over the set (full-tissue z, Stage 0 descriptive).",
            "",
            "| tissue | score | Spearman age_mid | Spearman age ordinal | Spearman SMTSISCH | R² ~ C(SMCENTER) |",
            "|---|---|---:|---:|---:|---:|",
        ]
        for _, r in comp.iterrows():
            lines.append(
                f"| {r.tissue} | {r.score} | {_signed(r.spearman_age_mid)} | "
                f"{_signed(r.spearman_age_ordinal)} | {_signed(r.spearman_SMTSISCH)} | "
                f"{_signed(r.r2_SMCENTER)} |"
            )
        lines.append("")

    lines += [
        "## Stage 1 gate",
        "",
        f"- FINDINGS_POSCTRL.md exists: {bool(gate.get('posctrl_exists'))}",
        f"- C1 ran: {gate.get('c1_ran')} (source={gate.get('source')})",
        f"- post-supersession status: {gate.get('status')}",
        f"- open: {gate.get('open')}",
        f"- reason: {gate.get('reason')}",
        "",
    ]
    tr = gate.get("transfer") or {}
    if tr:
        metric = tr.get("metric", "r2_legacy")
        if metric == "pearson_r":
            lines.append(
                f"S2 transfer used for the gate (Pearson r): real-age gene_ridge "
                f"HBCC→MSSM {_signed(tr.get('H_to_M'))} (null {_signed(tr.get('H_to_M_null'))}), "
                f"MSSM→HBCC {_signed(tr.get('M_to_H'))} (null {_signed(tr.get('M_to_H_null'))}); "
                f"uncalibrated R² {_signed(tr.get('uncal_R2_H_to_M'))} / {_signed(tr.get('uncal_R2_M_to_H'))} "
                f"(calibration artifact). Planted dense_shared r "
                f"{_signed(tr.get('planted_H_to_M'))} / {_signed(tr.get('planted_M_to_H'))} "
                f"(null {_signed(tr.get('planted_H_to_M_null'))} / {_signed(tr.get('planted_M_to_H_null'))})."
            )
        else:
            lines.append(
                f"C1 dense_shared transfer used for the gate: HBCC→MSSM {_signed(tr.get('H_to_M'))} "
                f"(null {_signed(tr.get('H_to_M_null'))}), "
                f"MSSM→HBCC {_signed(tr.get('M_to_H'))} "
                f"(null {_signed(tr.get('M_to_H_null'))}); "
                f"achieved within-site R² {_signed(tr.get('achieved_within_site_r2'))} "
                f"(target {_signed(tr.get('target_r2'))}); "
                f"bootstrap {tr.get('boot_median_angle')}° vs sex ref 41.3°."
            )
        lines.append("")

    if gate.get("open"):
        lines += [
            f"## Pre-registration {PREREG_DATE} (verbatim, written before any Stage 1 fit)",
            "",
            PREREG_BARS,
            "",
            "## Pre-registered reading (verbatim, r in place of R²)",
            "",
            PREREG_READING,
            "",
            "## Original Prompt B bars (superseded before any GTEx Stage 1 number)",
            "",
            PREREG_BARS_ORIGINAL,
            "",
            PREREG_READING_ORIGINAL,
            "",
        ]
        if not PREREG_FLAG.exists():
            PREREG_FLAG.write_text(
                f"pre-registration {PREREG_DATE} written before any Stage 1 T-cell fit\n"
                f"n_perm={N_PERM_GTEX}\nprimary=held-out Pearson r\n"
                f"pass = r > 0 with null <= 0.05\nangle reported only\n",
                encoding="utf-8",
            )

    stage1_ran = (GTEX_DIR / "stage1_summary.json").exists()
    if not gate.get("open"):
        lines += [
            "Stage 1 T-cells were not fit. Pre-registered bars and reading are not evaluated. "
            "They will be copied into this file **before** any T-cell is fit if the gate opens.",
            "",
        ]
    elif not stage1_ran:
        lines += [
            "Stage 1 T-cells have **not yet been fit**. The pre-registration above is frozen. "
            f"n_perm={N_PERM_GTEX}. Angle is reported only (sex ref {SEX_ANGLE_REF}°, "
            f"planted n=233 {PLANTED_ANGLE_N233}°).",
            "",
        ]
    else:
        t1_path = GTEX_DIR / "t1.csv" if (GTEX_DIR / "t1.csv").exists() else GTEX_DIR / "t1.json"
        if t1_path.exists():
            t1 = pd.read_csv(t1_path) if t1_path.suffix == ".csv" else pd.read_json(t1_path)
            cols = [c for c in t1.columns if c not in ("per_center", "folds")]
            lines += ["## T1", "", md_table(t1[cols]), ""]
        if (GTEX_DIR / "t2.json").exists():
            t2j = load_json(GTEX_DIR / "t2.json")
            if isinstance(t2j, dict) and t2j.get("skipped"):
                lines += ["## T2", "", f"SKIPPED: {t2j.get('reason')}", ""]
            elif (GTEX_DIR / "t2.csv").exists():
                t2 = pd.read_csv(GTEX_DIR / "t2.csv")
                lines += ["## T2", "", md_table(t2), ""]
        elif (GTEX_DIR / "t2.csv").exists():
            t2 = pd.read_csv(GTEX_DIR / "t2.csv")
            lines += ["## T2", "", md_table(t2), ""]
        if (GTEX_DIR / "t3_t4.csv").exists():
            t34 = pd.read_csv(GTEX_DIR / "t3_t4.csv")
            lines += ["## T3 / T4", "", md_table(t34), ""]
        if (GTEX_DIR / "t5.json").exists():
            t5 = load_json(GTEX_DIR / "t5.json")
            lines += ["## T5", "", f"```\n{json.dumps(t5, indent=1, default=str)}\n```", ""]
        if (GTEX_DIR / "t6.json").exists():
            t6 = load_json(GTEX_DIR / "t6.json")
            lines += ["## T6", "", f"```\n{json.dumps(t6, indent=1, default=str)}\n```", ""]
        if (GTEX_DIR / "t1_comp.csv").exists():
            tc = pd.read_csv(GTEX_DIR / "t1_comp.csv")
            lines += ["## Composition-only rerun (BA9 raw passed)", "", md_table(tc), ""]
        verdict = ""
        vp = GTEX_DIR / "stage1_verdict.txt"
        if vp.exists():
            verdict = vp.read_text(encoding="utf-8").strip()
        lines += [
            "## Verdict",
            "",
            verdict if verdict else "Only the pre-registered outcome that fired is stated here. See T1–T6 tables.",
            "",
            "## What this changes in FINDINGS_LOWDIM L4 (quote, not edited)",
            "",
            LOWDIM_L4_QUOTE,
            "",
        ]

    lines += [
        "## Limitations",
        "",
        CAVEATS,
        "",
        "4. Two covariate regimes (`raw`, `resid`) are both primary; `resid` is not selected after seeing numbers.",
        "5. Seed `20260914`. SVD-LOO ridge α grid and PLS-1 are the geometry/target methods called on a different X.",
        "6. Public AGE is a bin; Spearman with the ordinal bin is reported alongside R². Year-level MAE is not used.",
        "7. `SMCENTER` is collection site, not BA9 dissection site (Miami Brain Endowment Bank).",
        "8. Exact ages require dbGaP `phs000424.v10` and are out of scope.",
        "",
        "## Files",
        "",
        "| path | content |",
        "|---|---|",
        "| `src/gtex_common.py`, `gtex_stage0.py`, `gtex_stage1.py`, `gtex_findings.py` | code |",
        "| `notebooks/gtex_stage0.ipynb`, `gtex_stage1.ipynb` | runnable from a clean checkout |",
        "| `results/gtex/` | manifest, cohort, audit, composition, figures |",
        "| `FINDINGS_GTEX.md` | this file |",
        "| `PROGRESS_GTEX.md` | resume state |",
        "",
    ]

    FINDINGS_PATH.write_text("\n".join(lines), encoding="utf-8")
    print("wrote", FINDINGS_PATH)


if __name__ == "__main__":
    write_findings()
