"""Write FINDINGS_FIBRO.md from results/fibro/. Does not modify other FINDINGS*.md."""
from __future__ import annotations

import sys
from pathlib import Path

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parent))
from fibro_common import (  # noqa: E402
    FIBRO_DIR, FIBRO_FIG, FIBRO_SEED, BOOT_SEED, N_PERM, N_BOOT, N_RANDOM_DIR,
    SMTSD_FIBRO, EXCLUDE_SITE, EXCLUDE_SITE_REASON, CAVEAT_AGE_BIN,
    PREREG_STAGE1, PREREG_STAGE2, PREREG_STAGE1_FLAG, PREREG_STAGE2_FLAG,
    FROZEN_RULER, FINDINGS_PATH, PLURI_ENDOGENOUS, FIBRO_IDENTITY, OSKM_FAMILY,
    GSE_ACCESSION, RELEASE, RNASEQ_PIPELINE, DBGAP, load_json, md_table,
    fmt, fmt_ci, jsonable, progress_snapshot,
)


def _signed(x, d=3):
    try:
        v = float(x)
    except (TypeError, ValueError):
        return "NA"
    if not np.isfinite(v):
        return "NA"
    return f"{v:+.{d}f}"


def _md_from_csv(path, cols):
    if not path.exists():
        return "_file missing._"
    df = pd.read_csv(path)
    keep = [c for c in cols if c in df.columns]
    return md_table(df[keep]) if keep else "_no matching columns._"


def write_findings():
    man = load_json(FIBRO_DIR / "manifest.json") if (FIBRO_DIR / "manifest.json").exists() else {}
    gate = load_json(FIBRO_DIR / "stage1_gate.json") if (FIBRO_DIR / "stage1_gate.json").exists() else man.get("stage1_gate") or {}
    failures = man.get("failures") or []
    status = man.get("status", "UNKNOWN")
    passed = bool(gate.get("pass")) if gate else False

    if status == "STOP" or (failures and not (FIBRO_DIR / "stage1_gate.json").exists()):
        status_line = (
            f"**Status:** STOP. Seed `{FIBRO_SEED}`. {RELEASE} ({RNASEQ_PIPELINE}; {DBGAP}). "
            f"`SMTSD == {SMTSD_FIBRO!r}`. Failures in `results/fibro/manifest.json`."
        )
    elif gate and not passed:
        status_line = (
            f"**Status:** Stage 1 fail. Seed `{FIBRO_SEED}`. boot `{BOOT_SEED}`. "
            f"{RELEASE} ({RNASEQ_PIPELINE}; {DBGAP}). Primary metric: Spearman ρ "
            f"(n_perm={N_PERM}). GSE297234 not opened."
        )
    elif gate and passed and not (FIBRO_DIR / "stage2_summary.json").exists():
        status_line = (
            f"**Status:** Stage 1 pass; Stage 2 not complete. Seed `{FIBRO_SEED}`. "
            f"Primary metric: Spearman ρ (n_perm={N_PERM})."
        )
    elif (FIBRO_DIR / "stage2_summary.json").exists():
        s2 = load_json(FIBRO_DIR / "stage2_summary.json")
        key = (s2.get("outcome") or {}).get("key", "")
        status_line = (
            f"**Status:** Stage 2 complete (outcome `{key}`). Seed `{FIBRO_SEED}`. "
            f"boot `{BOOT_SEED}`. {RELEASE} ({RNASEQ_PIPELINE}; {DBGAP}). "
            f"Primary Stage 1 metric: Spearman ρ (n_perm={N_PERM}). "
            f"Stage 2 random-direction n={N_RANDOM_DIR}."
        )
    else:
        status_line = (
            f"**Status:** in progress (`{status}`). Seed `{FIBRO_SEED}`. "
            f"{RELEASE} ({RNASEQ_PIPELINE}; {DBGAP})."
        )

    lines = [
        "# FINDINGS_FIBRO — fibroblast age ruler, then the reprogramming curve",
        "",
        status_line,
        "",
        "Does not modify any existing `FINDINGS_*.md` or `FALSIFICATION.md`. "
        "No GSE325735. Uncalibrated R² is reported and is never a gate. "
        "PLS-1 is reported and is not averaged with ridge.",
        "",
        f"Reproduced by `src/fibro_run.py`. Freeze: `SMAFRZE==RNASEQ`. Tissue: `{SMTSD_FIBRO}`.",
        "",
        "## Stage 1 pre-registration (verbatim, written before any Stage 1 fit)",
        "",
        PREREG_STAGE1,
        "",
        f"Flag: `results/fibro/{PREREG_STAGE1_FLAG.name}` exists={PREREG_STAGE1_FLAG.exists()}.",
        "",
        f"**Caveat:** {CAVEAT_AGE_BIN}",
        "",
    ]

    if failures:
        lines += ["## STOP / failures", "", "Not substituting columns or repairing rows.", ""]
        for f in failures:
            lines.append(f"- **{f.get('step')}:** {f.get('message')}")
        if (FIBRO_DIR / "stage2_summary.json").exists():
            lines.append(
                "The earlier MTX-triplet stop was a format mismatch; Cell Ranger "
                "`filtered_feature_bc_matrix.h5` in `GSE297234_RAW.tar` was used. "
                "Seurat `.rds` files were listed on GEO and not loaded."
            )
            lines.append("")

    coh_p = FIBRO_DIR / "stage1_cohort.csv"
    if coh_p.exists():
        coh = pd.read_csv(coh_p).iloc[0]
        cj = load_json(FIBRO_DIR / "stage1_cohort.json") if (FIBRO_DIR / "stage1_cohort.json").exists() else {}
        counts = (cj.get("SMCENTER_counts") if isinstance(cj, dict) else None) or {}
        lines += [
            "## Stage 1 — cohort",
            "",
            f"- SMTSD `{SMTSD_FIBRO}`  n_samples={int(coh.n_samples)}  n_donors={int(coh.n_donors)}  "
            f"one sample per donor.",
            f"- SMCENTER (file strings): `{counts}`  n_SMNABTCH={int(coh.n_SMNABTCH)}  n_SMGEBTCH={int(coh.n_SMGEBTCH)}",
            f"- SEX 1/2 (GTEx coding) = {int(coh.n_male)}/{int(coh.n_female)}",
            f"- SMRIN median {float(coh.SMRIN_median):.2f} "
            f"({float(coh.SMRIN_min):.2f}–{float(coh.SMRIN_max):.2f}, miss={int(coh.SMRIN_n_missing)})",
            f"- SMTSISCH median {float(coh.SMTSISCH_median):.1f} "
            f"({float(coh.SMTSISCH_min):.1f}–{float(coh.SMTSISCH_max):.1f}, miss={int(coh.SMTSISCH_n_missing)})",
            f"- DTHHRDY miss={int(coh.DTHHRDY_n_missing)}",
            f"- {EXCLUDE_SITE} n={int(coh.D1_n)}. {EXCLUDE_SITE_REASON}",
            "",
            "AGE bins (n samples):",
            "",
        ]
        bins = ["20-29", "30-39", "40-49", "50-59", "60-69", "70-79"]
        lines.append("| " + " | ".join(bins) + " |")
        lines.append("|" + "|".join(["---:"] * len(bins)) + "|")
        lines.append("| " + " | ".join(str(int(coh[f"n_AGE_{b}"])) for b in bins) + " |")
        lines.append("")
        cols_read = []
        if (FIBRO_DIR / "stage1_cohort.json").exists():
            blob = load_json(FIBRO_DIR / "stage1_cohort.json")
            cols_read.append(f"- sample_attributes: {blob.get('sample_columns_read')}")
            cols_read.append(f"- subject_phenotypes: {blob.get('phenotype_columns_read')}")
        if cols_read:
            lines += ["Columns actually read:", ""] + cols_read + [""]

    ov_p = FIBRO_DIR / "stage1_gene_overlap.json"
    if ov_p.exists():
        ov = load_json(ov_p)
        lines += [
            "## Stage 1 — gene filter and DLPFC overlap",
            "",
            f"- n_genes_raw={ov.get('n_genes_raw')}  n_genes_filter={ov.get('n_genes_filter')}  "
            f"n_overlap_dlpfc={ov.get('n_overlap_dlpfc')}  DLPFC n={ov.get('n_dlpfc')}",
            f"- filter: {ov.get('gene_filter')}  normalization: {ov.get('normalization')}",
            f"- DLPFC columns actually read: {ov.get('dlpfc_columns_read')}",
            "",
        ]

    cv_p = FIBRO_DIR / "stage1_cv.csv"
    if cv_p.exists():
        lines += [
            "## Stage 1 — within-site CV (SMNABTCH-grouped; batch-held-out and donor-held-out asserted)",
            "",
            "Every cell has its permutation null. Uncalibrated R² is not a gate.",
            "",
            _md_from_csv(cv_p, [
                "site", "regime", "method", "rho", "rho_null", "rho_p",
                "rho_ci_lo", "rho_ci_hi", "cal_r2", "r2", "r", "n", "n_perm", "pass_rho",
            ]),
            "",
        ]
        skipped = (load_json(FIBRO_DIR / "stage1_cv.json") or {}).get("skipped") if (FIBRO_DIR / "stage1_cv.json").exists() else None
        if skipped:
            lines.append(f"Skipped sites: `{skipped}`")
            lines.append("")

    tr_p = FIBRO_DIR / "stage1_transfer.csv"
    if tr_p.exists():
        lines += [
            "## Stage 1 — B1↔C1 transfer",
            "",
            f"{EXCLUDE_SITE} excluded (stated). PLS-1 is not averaged with ridge.",
            "",
            _md_from_csv(tr_p, [
                "direction", "regime", "method", "rho", "rho_null", "rho_p",
                "rho_ci_lo", "rho_ci_hi", "cal_r2", "r2", "r",
                "n_train", "n_test", "n_perm", "pass_rho",
            ]),
            "",
        ]

    lines += ["## Stage 1 gate outcome", ""]
    if gate:
        lines.append(
            f"- pass={gate.get('pass')}  method={gate.get('method')}  regime={gate.get('regime')}"
        )
        lines.append(
            f"- B1→C1 ρ={_signed(gate.get('rho_B1_to_C1'))} "
            f"(null {_signed(gate.get('rho_null_B1_to_C1'))} p={_signed(gate.get('rho_p_B1_to_C1'))} "
            f"CI={fmt_ci(*(gate.get('rho_ci_B1_to_C1') or [None, None]))}) "
            f"calR²={_signed(gate.get('cal_r2_B1_to_C1'))} R²={_signed(gate.get('r2_B1_to_C1'))} "
            f"pass={gate.get('pass_B1_to_C1')}"
        )
        lines.append(
            f"- C1→B1 ρ={_signed(gate.get('rho_C1_to_B1'))} "
            f"(null {_signed(gate.get('rho_null_C1_to_B1'))} p={_signed(gate.get('rho_p_C1_to_B1'))} "
            f"CI={fmt_ci(*(gate.get('rho_ci_C1_to_B1') or [None, None]))}) "
            f"calR²={_signed(gate.get('cal_r2_C1_to_B1'))} R²={_signed(gate.get('r2_C1_to_B1'))} "
            f"pass={gate.get('pass_C1_to_B1')}"
        )
        lines.append(f"- reason: {gate.get('reason')}")
        lines.append("")
        if not passed:
            lines.append(
                "The fibroblast ruler does not transfer between GTEx sites so the "
                "reprogramming projection would be uninterpretable. GSE297234 was not opened."
            )
            lines.append("")
    else:
        lines += ["Gate not scored.", ""]

    lines += [
        "## Stage 2 pre-registration (verbatim, written before the first projection; gated on Stage 1)",
        "",
        PREREG_STAGE2,
        "",
        f"Flag: `results/fibro/{PREREG_STAGE2_FLAG.name}` exists={PREREG_STAGE2_FLAG.exists()}.",
        f"Frozen ruler: `{FROZEN_RULER}` exists={FROZEN_RULER.exists()}.",
        "",
        "Frozen lists:",
        f"- pluripotency endogenous: {', '.join(PLURI_ENDOGENOUS)}",
        f"- fibroblast identity: {', '.join(FIBRO_IDENTITY)}",
        f"- OSKM-family (dropped for without-OSKM score): {', '.join(OSKM_FAMILY)}",
        "",
    ]

    s2p = FIBRO_DIR / "stage2_summary.json"
    if not s2p.exists():
        if passed:
            lines += ["Stage 2 did not complete.", ""]
        else:
            lines += ["Stage 2 not run (Stage 1 gate closed).", ""]
    else:
        s2 = load_json(s2p)
        sendai = s2.get("sendai") or {}
        lines += [
            "## Sendai contamination check",
            "",
            f"- n_hits={sendai.get('n_hits')}  n_features_scanned={sendai.get('n_features_scanned')}",
            f"- n_ensg={sendai.get('n_ensg')}  n_non_ensg={sendai.get('n_non_ensg')}  "
            f"genomes={sendai.get('genomes')}  feature_types={sendai.get('feature_types')}",
            f"- columns actually read: {sendai.get('columns_actually_read')}",
            f"- matrix format: {s2.get('matrix_format') or man.get('stage2_matrix_format')}",
            f"- names: `{sendai.get('names')}`",
            f"- {sendai.get('note')}",
            "",
        ]
        ver = s2.get("verify") or {}
        if ver.get("discrepancies"):
            lines += ["## GSE297234 record vs scoping description", "",
                      f"Discrepancies (record used): `{ver.get('discrepancies')}`", ""]
        lines += [
            "## Stage 2 — all-cell pseudobulk (a), primary",
            "",
            _md_from_csv(FIBRO_DIR / "stage2_allcell_scores.csv", [
                "cell_line", "day", "age_years", "n_cells", "age_score",
                "pluri_with_OSKM", "pluri_without_OSKM", "pluri_primary", "gsm",
            ]),
            "",
        ]
        aged = s2.get("aged") or {}
        young = s2.get("young") or {}
        san = s2.get("sanity") or {}
        lines += [
            "### Day-0 sanity and aged-donor nulls",
            "",
            f"- day0 aged={_signed(san.get('day0_aged'), 4)} young={_signed(san.get('day0_young'), 4)} "
            f"aged>young={san.get('aged_gt_young')} (sanity test only)",
            f"- aged GM00731 days={aged.get('days')} age_score={aged.get('age_score')}",
            f"- aged decline 0→10={_signed(aged.get('decline_0_to_10'), 4)}  "
            f"Spearman score vs day={_signed(aged.get('spearman_vs_day'))}",
            f"- primary null endpoint p={_signed(aged.get('p_endpoint'))}  "
            f"n_random≥real={aged.get('n_null_ge_real')}/{aged.get('n_random')}",
            f"- secondary null monotonicity p={_signed(aged.get('p_monotonicity'))}",
            f"- young GM23815 days={young.get('days')} age_score={young.get('age_score')} "
            f"decline 0→10={_signed(young.get('decline_0_to_10'), 4)} "
            f"Spearman vs day={_signed(young.get('spearman_vs_day'))} (contrast; not pooled)",
            f"- pluripotency primary={s2.get('pluri_primary')}  with/without OSKM disagree={s2.get('pluri_disagree')}",
            "",
        ]
        pm_p = FIBRO_DIR / "stage2_pluri_meta.json"
        if pm_p.exists():
            pm = load_json(pm_p)
            w = (pm.get("with_OSKM") or {}).get("pluri") or {}
            o = (pm.get("without_OSKM") or {}).get("pluri") or {}
            f = (pm.get("with_OSKM") or {}).get("fibro") or {}
            lines += [
                "Frozen pluripotency genes actually mapped onto the GTEx fibroblast ruler:",
                f"- with OSKM: mapped={w.get('mapped')} missing={w.get('missing')}",
                f"- without OSKM: mapped={o.get('mapped')} missing={o.get('missing')}",
                f"- fibroblast identity: mapped={f.get('mapped')} missing={f.get('missing')}",
                "",
            ]
        if (FIBRO_DIR / "stage2_cluster_scores.csv").exists():
            lines += [
                "## Stage 2 — cluster pseudobulks (b), not averaged into (a)",
                "",
                _md_from_csv(FIBRO_DIR / "stage2_cluster_scores.csv", [
                    "cell_line", "day", "cluster", "n_cells", "age_score", "pluri_primary",
                ]),
                "",
            ]
        fig = FIBRO_FIG / "plane.png"
        if fig.exists():
            lines += ["## Figure", "", f"![]({fig.as_posix().split('Age–Identity Separability Benchmark/')[-1] if False else 'results/fibro/figures/plane.png'})", ""]
        out = s2.get("outcome") or {}
        lines += ["## Verdict", "", out.get("text", "missing"), "",
                  f"Fired key: `{out.get('key')}`. drop_time={out.get('drop_time')} "
                  f"rise_time={out.get('rise_time')} monotonic={out.get('monotonic')} "
                  f"beats_null={out.get('beats_null')} sanity_ok={out.get('sanity_ok')}.",
                  ""]

    lines += [
        "## Limitations",
        "",
        "1. Two donors in GSE297234. The time course can show a bend; it cannot estimate donor-to-donor variance.",
        "2. " + CAVEAT_AGE_BIN,
        "3. Cross-platform shift from GTEx bulk polyA (RNASeQCv2.4.2) to 10x 3' scRNA-seq pseudobulk.",
        "4. One cell type (cultured fibroblast), so no identity(type) ruler. Differentiation uses a frozen pluripotency − fibroblast-identity score.",
        "5. Sendai contamination status is in the Sendai section above; vector reads can map onto endogenous OSKM genes.",
        "6. D1 (n=7) is excluded from transfer and within-site CV, stated, and is included in the frozen all-donor ruler if Stage 1 passed.",
        "7. Seed `20260914` (permutation / random directions); donor-bootstrap seed `20260918`. n_perm=200, n_boot=200, n_random=200.",
        "8. Missing overlap genes left at z=0. Nothing fitted on GSE297234. GSE325735 not opened.",
        "9. Frozen GTEx fibroblast gene set (n=23485 after Stage 0 filter) did not contain NANOG, LIN28A, DPPA4, ZFP42. "
        "The with-OSKM pluripotency mean used POU5F1, SALL4, DNMT3B; the without-OSKM mean used SALL4, DNMT3B.",
        "10. Donor ages 22 and 96 are from `Series_overall_design`, not `Sample_characteristics_ch1` (those fields had no numeric age).",
        "",
        "## Files",
        "",
        "| path | content |",
        "|---|---|",
        "| `src/fibro_common.py`, `fibro_stage1.py`, `fibro_stage2.py`, `fibro_findings.py`, `fibro_run.py` | code |",
        "| `results/fibro/` | manifest, tables, frozen ruler, logs, `figures/plane.png` |",
        "| `FINDINGS_FIBRO.md` | this file |",
        "| `PROGRESS_FIBRO.md` | resume state |",
        "",
    ]
    FINDINGS_PATH.write_text("\n".join(lines) + "\n", encoding="utf-8")
    return FINDINGS_PATH


if __name__ == "__main__":
    p = write_findings()
    progress_snapshot("done" if (FIBRO_DIR / "stage2_summary.json").exists() or
                      ((FIBRO_DIR / "stage1_gate.json").exists() and not load_json(FIBRO_DIR / "stage1_gate.json").get("pass"))
                      else "continue Stage 1/2",
                      extra=f"FINDINGS_FIBRO.md written from disk → {p}")
    print(p)
