"""Write FINDINGS_FIBRO2.md from results/fibro2/. Does not modify other FINDINGS*.md."""
from __future__ import annotations

import sys
from pathlib import Path

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parent))
from fibro2_common import (  # noqa: E402
    FIBRO2_DIR, FIBRO2_SEED, FIBRO2_BOOT, N_PERM, N_BOOT, N_RANDOM_DIR,
    PREREG_TA, PREREG_TB, PREREG_TC, PREREG_TA_FLAG, PREREG_TB_FLAG, PREREG_TC_FLAG,
    FINDINGS_PATH, FROZEN_RULER, STAGE2_VERDICT_SENTENCE,
    load_json, load_manifest, jsonable, progress_snapshot,
)
from target_common import md_table  # noqa: E402
from fibro_common import fmt_ci  # noqa: E402


def _signed(x, d=3):
    try:
        v = float(x)
    except (TypeError, ValueError):
        return "NA"
    if not np.isfinite(v):
        return "NA"
    return f"{v:+.{d}f}"


def _md(path, cols=None):
    if not path.exists():
        return "_file missing._"
    df = pd.read_csv(path)
    if cols:
        keep = [c for c in cols if c in df.columns]
        df = df[keep] if keep else df
    if df.empty:
        return "_empty._"
    return md_table(df)


def _status_line(man):
    ta = (man.get("ta_reading") or {})
    tb = (man.get("tb_reading") or {})
    tc = (man.get("tc_reading") or {})
    ta_s = ta.get("key", "T-A not scored")
    if ta.get("interpretable") is True:
        ta_txt = "T-A: interpretable (primary ρ pass)."
    elif ta.get("interpretable") is False:
        ta_txt = f"T-A: not interpretable (`{ta_s}`)."
    else:
        ta_txt = f"T-A: `{ta_s}`."
    tb_key = tb.get("key", "T-B not scored")
    if tb.get("d10_technical") is True:
        tb_txt = f"T-B: d10 technical (`{tb_key}`)."
    elif tb.get("d10_technical") is False:
        tb_txt = f"T-B: d10 not a material technical difference (`{tb_key}`)."
    else:
        tb_txt = f"T-B: `{tb_key}`."
    if tc.get("linearity_limiting") is True:
        tc_txt = f"T-C: linearity limiting (`{tc.get('key')}`)."
    elif tc.get("linearity_limiting") is False:
        tc_txt = f"T-C: linearity not the limiting factor (`{tc.get('key')}`)."
    else:
        tc_txt = f"T-C: `{tc.get('key', 'not scored')}`."
    return (
        f"**Status:** {ta_txt} {tb_txt} {tc_txt} "
        f"Seed `{FIBRO2_SEED}`. boot `{FIBRO2_BOOT}`. n_perm={N_PERM}. n_boot={N_BOOT}. "
        f"n_random={N_RANDOM_DIR}. Frozen ruler `{FROZEN_RULER.name}` exists={FROZEN_RULER.exists()}."
    )


def write_findings():
    man = load_manifest() if (FIBRO2_DIR / "manifest.json").exists() else {}
    ta_sum = load_json(FIBRO2_DIR / "ta_summary.json") if (FIBRO2_DIR / "ta_summary.json").exists() else {}
    tb_sum = load_json(FIBRO2_DIR / "tb_summary.json") if (FIBRO2_DIR / "tb_summary.json").exists() else {}
    tc_sum = load_json(FIBRO2_DIR / "tc_summary.json") if (FIBRO2_DIR / "tc_summary.json").exists() else {}
    ex_sum = load_json(FIBRO2_DIR / "extrap_summary.json") if (FIBRO2_DIR / "extrap_summary.json").exists() else {}

    lines = [
        "# FINDINGS_FIBRO2 — is the fibroblast ruler the instrument, or is the result real?",
        "",
        _status_line(man),
        "",
        "Does not modify any existing `FINDINGS_*.md` or `FALSIFICATION.md`. No GSE325735. "
        "Uncalibrated R² is reported and is never a gate. Nothing averaged across cohorts, donors, or regimes. "
        "T-A / T-B / T-C cannot overturn the Stage 2 `no_decline` verdict by themselves; they establish whether it is interpretable.",
        "",
        f"Reproduced by `src/fibro2_run.py`. Frozen ruler: `{FROZEN_RULER}` exists={FROZEN_RULER.exists()}.",
        "",
        f"Flag: `results/fibro2/PREREG_TA.flag` exists={PREREG_TA_FLAG.exists()}. "
        f"`PREREG_TB.flag` exists={PREREG_TB_FLAG.exists()}. "
        f"`PREREG_TC.flag` exists={PREREG_TC_FLAG.exists()}.",
        "",
        "## T-A pre-registration (verbatim, written before any T-A download or ρ)",
        "",
        PREREG_TA,
        "",
        "## T-A — cohort table",
        "",
        _md(FIBRO2_DIR / "ta_cohorts.csv", [
            "accession", "role", "title", "platform_kind", "types", "platforms",
            "n_gsm", "n_with_age", "n_missing_age", "n_adult", "n_child",
            "n_donors_with_age", "usable", "unusable_reason",
        ]),
        "",
        "## T-A — results table",
        "",
        _md(FIBRO2_DIR / "ta_results.csv", [
            "accession", "role", "platform_kind", "tag", "n_samples", "n_donors",
            "n_excluded_no_age", "n_excluded_child", "n_overlap", "n_missing_z0",
            "rho", "rho_null", "rho_p", "rho_ci_lo", "rho_ci_hi",
            "cal_r2", "r2", "r", "pass_rho", "n_perm", "n_boot",
            "used_10x_only", "platform_detail",
        ]),
        "",
        "## T-A — reading",
        "",
    ]
    reading_ta = (ta_sum.get("reading") or man.get("ta_reading") or {})
    if reading_ta:
        lines += [
            f"Fired key: `{reading_ta.get('key')}`. interpretable={reading_ta.get('interpretable')}.",
            "",
            reading_ta.get("text", ""),
            "",
        ]
    else:
        lines += ["_T-A not scored._", ""]
    if ta_sum.get("sc_resolved") is False:
        lines += [
            f"Single-cell fibroblast cohort resolved this session: **False**. "
            f"reason: {ta_sum.get('sc_reason')}",
            "",
            "Not substituting a bulk cohort and declaring the platform question closed.",
            "",
        ]
    elif ta_sum.get("sc_resolved") is True:
        lines += ["Single-cell fibroblast cohort resolved this session: **True**.", ""]

    lines += [
        "## T-B pre-registration (verbatim, written before any T-B QC table or d0→d7 null)",
        "",
        PREREG_TB,
        "",
        "## T-B — QC table (per donor × timepoint)",
        "",
        _md(FIBRO2_DIR / "tb_qc.csv", [
            "cell_line", "day", "gsm", "age_years", "n_cells",
            "median_UMI_per_cell", "median_genes_per_cell",
            "mito_read_fraction_library", "mito_read_fraction_median_cell",
            "n_mito_genes", "n_cluster_0", "n_cluster_1", "n_cluster_2",
            "chemistry", "software",
        ]),
        "",
    ]
    mito_p = FIBRO2_DIR / "tb_mito_rule.json"
    if mito_p.exists():
        mito = load_json(mito_p)
        lines += [
            f"Mitochondrial rule: `{mito.get('rule')}` n_mito_genes={mito.get('n_mito')} "
            f"symbols={mito.get('mito_symbols')}.",
            "",
        ]
    lines += [
        "## T-B — random-direction null table",
        "",
        _md(FIBRO2_DIR / "tb_null.csv", [
            "variant", "cell_line", "endpoint", "primary", "rank",
            "n_cells_d0", "n_cells_end", "score_d0", "score_end",
            "decline", "p", "n_random_ge_real", "n_random", "pass_p", "ok",
        ]),
        "",
        "## T-B — reading",
        "",
    ]
    reading_tb = (tb_sum.get("reading") or man.get("tb_reading") or {})
    if reading_tb:
        lines += [
            f"Fired key: `{reading_tb.get('key')}`. d10_technical={reading_tb.get('d10_technical')}. "
            f"p_d0d7={reading_tb.get('p7')} p_d0d10={reading_tb.get('p10')}.",
            "",
            reading_tb.get("text", ""),
            "",
            "A post-hoc d7 endpoint is not the primary Stage 2 result.",
            "",
        ]
    else:
        lines += ["_T-B not scored._", ""]
    qc_diff = tb_sum.get("qc_diff") or {}
    if qc_diff.get("notes"):
        lines += ["d10 vs d0/d3/d7 (from disk):", ""]
        for n in qc_diff["notes"]:
            lines.append(f"- {n}")
        lines.append("")

    lines += [
        "## T-C pre-registration (verbatim, written before any boosting fit)",
        "",
        PREREG_TC,
        "",
        "## T-C — ridge vs HistGradientBoosting",
        "",
        _md(FIBRO2_DIR / "tc_side_by_side.csv", [
            "cell", "site", "ridge_rho", "ridge_rho_null", "ridge_p",
            "ridge_ci_lo", "ridge_ci_hi", "ridge_cal_r2", "ridge_r2",
            "hgb_rho", "hgb_rho_null", "hgb_p", "hgb_ci_lo", "hgb_ci_hi",
            "hgb_cal_r2", "hgb_r2", "n",
        ]),
        "",
        "## T-C — reading",
        "",
    ]
    reading_tc = (tc_sum.get("reading") or man.get("tc_reading") or {})
    if reading_tc:
        lines += [
            f"Fired key: `{reading_tc.get('key')}`. linearity_limiting={reading_tc.get('linearity_limiting')}.",
            "",
            reading_tc.get("text", ""),
            "",
            "HGB null: donor-level shuffle of ages; fitted scores held fixed (n_perm=200, seed `20260914`). "
            "Ridge Stage 1 p-values on disk used refit-on-permuted-y. T-C reading compares HGB transfer ρ "
            "to the ridge bootstrap CI, not to ridge p-values. D1 n=7 skipped (n<10), stated.",
            "",
        ]
        if reading_tc.get("comps"):
            lines.append("Per-direction comparison (from disk):")
            lines.append("")
            for c in reading_tc["comps"]:
                lines.append(
                    f"- {c.get('direction')} hgb ρ={_signed(c.get('hgb_rho'))} "
                    f"ridge ρ={_signed(c.get('ridge_rho'))} "
                    f"ridge CI=[{_signed(c.get('ridge_ci_lo'))}, {_signed(c.get('ridge_ci_hi'))}] "
                    f"inside={c.get('inside')} above={c.get('above')}"
                )
            lines.append("")
    else:
        lines += ["_T-C not scored._", ""]
    lines += [
        "## Direction vs score",
        "",
        tc_sum.get("direction_vs_score") or (
            "A nonlinear model yields a score, not a direction; the project's plane, angles, "
            "and \"age movement per unit identity loss\" criterion all require a direction, so a "
            "nonlinear model cannot be substituted into Stage 2 without discarding the geometry."
        ),
        "",
        "## Extrapolation (report-only, no gate)",
        "",
        "Mahalanobis is in the leading PCA of the GTEx fibroblast training z-space. "
        "Nearest neighbour is Euclidean in that z-space. Missing overlap genes at z=0.",
        "",
        "GSE297234 all-cell rows:",
        "",
    ]
    extrap_p = FIBRO2_DIR / "extrap_distances.csv"
    if extrap_p.exists():
        ed = pd.read_csv(extrap_p)
        ac = ed[(ed.cohort == "GSE297234") & (ed.variant == "all_cell")].copy()
        cols = [c for c in (
            "cohort", "variant", "cell_line", "day", "n_cells", "age_score",
            "nn_euclidean", "mahalanobis_pca", "pca_k",
        ) if c in ac.columns]
        lines.append(md_table(ac[cols]) if len(ac) and cols else "_no all-cell rows._")
        lines += ["", "Cohort summaries:", "", _md(FIBRO2_DIR / "extrap_summary_by_cohort.csv")]
        cl = ed[(ed.cohort == "GSE297234") & (ed.variant == "cluster")].copy()
        if len(cl):
            lines += ["", "GSE297234 cluster rows:", ""]
            ccols = [c for c in (
                "cell_line", "day", "cluster", "n_cells", "age_score",
                "nn_euclidean", "mahalanobis_pca",
            ) if c in cl.columns]
            lines.append(md_table(cl[ccols]) if ccols else "_cluster columns missing._")
        lines.append("")
    else:
        lines += ["_extrap_distances.csv missing._", ""]

    lines += [
        "## unusable",
        "",
    ]
    unusable = man.get("unusable") or []
    if not unusable:
        lines.append("None recorded.")
    else:
        seen = set()
        for u in unusable:
            key = (u.get("accession"), u.get("reason"))
            if key in seen:
                continue
            seen.add(key)
            lines.append(f"- **{u.get('accession')}:** {u.get('reason')}")
    if (FIBRO2_DIR / "ta_GSE113957_samples.csv").exists():
        lines += [
            "",
            "GSE113957 record (from SOFT + files this session): n_gsm=143, n_with_stated_age=135, "
            "n_missing_age=8, n_stated_age_<18=24, n_adult≥18=111. Ages from Sample_title "
            "(e.g. years in the title). Series_type is high-throughput sequencing (GPL16791, GPL18573). "
            "Supplementary `GSE113957_fpkm.txt.gz` is FPKM, not integer counts "
            "(min=0 max=249213345 mean=1009452.0862567581). Series matrix HEAD 404. "
            "TMM/log2-CPM not applied. Primary external cell not scored.",
        ]
    if (FIBRO2_DIR / "ta_GSE307377_samples.csv").exists():
        lines += [
            "",
            "GSE307377 record (from SOFT + files this session): n_gsm=9, ages stated in "
            "`Sample_characteristics_ch1` key `age` (23–72 y). `GSE307377_raw.txt.gz` is HOMER "
            "`analyzeRepeats.pl` output: gene IDs are `NM_`, columns include `start`/`end`/`Length`/`Copies` "
            "and `tags/*` with half-integer values. Not integer counts; TMM/log2-CPM not applied. "
            "`GSE307377_tpm.txt.gz` is TPM. Series matrix table has 0 genes. n=9 is the SCOPING_TISSUE "
            "`near_miss` row. SOFT `cell line: Primary` is not a donor id; ages are per GSM. Not scored.",
        ]
    lines += ["", "## unverified", ""]
    unverified = man.get("unverified") or []
    if not unverified:
        lines.append("None recorded.")
    else:
        for u in unverified:
            lines.append(f"- **{u.get('accession')}:** {u.get('reason')}")
    lines += [
        "",
        "## What this changes in FINDINGS_FIBRO.md (quote, not edited)",
        "",
        "Stage 2 verdict sentence, quoted:",
        "",
        f"> {STAGE2_VERDICT_SENTENCE}",
        "",
        "This file does not edit `FINDINGS_FIBRO.md`. T-A/T-B/T-C do not replace that sentence. "
        "They say whether it can be read as a statement about reprogramming.",
        "",
    ]
    if reading_ta.get("interpretable") is False:
        lines.append(
            "T-A: Stage 2 `no_decline` must not be reported as a finding about reprogramming "
            f"(`{reading_ta.get('key')}`)."
        )
    elif reading_ta.get("interpretable") is True:
        lines.append(
            "T-A: the Stage 2 `no_decline` verdict is interpretable as a statement about OSKM "
            "on the primary external cell. The platform caveat is in the T-A reading if single-cell failed."
        )
    else:
        lines.append("T-A reading not yet on disk.")
    lines.append("")
    if reading_tb.get("text"):
        lines.append(f"T-B: {reading_tb.get('text')}")
        lines.append("")
    if reading_tc.get("text"):
        lines.append(f"T-C: {reading_tc.get('text')}")
        lines.append("")

    lines += [
        "## Limitations",
        "",
        "1. Frozen ruler trained on GTEx V10 cultured fibroblasts, public `AGE` 10-year bins; ρ is against a coarse ordinal on GTEx and against stated years on external records.",
        "2. Cross-platform: GTEx bulk polyA RNASeQCv2.4.2 vs 10x 3' scRNA-seq pseudobulk (GSE297234) vs whatever T-A resolved.",
        "3. GSE297234 has two donors. T-B cannot estimate donor-to-donor variance of a technical effect.",
        "4. Cluster labels are independent per donor × timepoint; the secondary null pairs by cell-count rank, not by biological identity of a cluster.",
        "5. Missing overlap genes at z=0. Nothing fitted on GSE297234 or on T-A matrices.",
        "6. HistGradientBoosting is a score, not a direction; T-C is not a Stage 2 replacement.",
        "7. Seed `20260914` (permutation / random directions / HGB). Donor-bootstrap seed `20260918`. n_perm=200, n_boot=200, n_random=200.",
        "8. GSE325735 not opened.",
        "9. Extrapolation distances have no pass/fail.",
        "10. A record that does not state per-sample donor age is unusable; ages are not inferred.",
        "11. GSE113957 supplementary is FPKM (transcript-level IDs), not integer counts; series matrix 404. Primary cell unscored.",
        "12. HGB permutation null holds scores fixed. Ridge Stage 1 null refit on permuted y.",
        "13. CXG atlas scored on `development_stage` (stated year-old stages); `Age` is also in obs and was not the column used.",
        "14. GSE226189 is bulk (GPL24676), not a substitute for the primary cell and not a 10x platform check.",
        "",
        "## Files",
        "",
        "| path | content |",
        "|---|---|",
        "| `src/fibro2_common.py`, `fibro2_ta.py`, `fibro2_tb.py`, `fibro2_tc.py`, `fibro2_extrap.py`, `fibro2_findings.py`, `fibro2_run.py` | code |",
        "| `results/fibro2/` | manifest, tables, logs, frozen-projection scores |",
        "| `FINDINGS_FIBRO2.md` | this file |",
        "| `PROGRESS_FIBRO2.md` | resume state |",
        "",
    ]
    FINDINGS_PATH.write_text("\n".join(lines) + "\n", encoding="utf-8")
    return FINDINGS_PATH


if __name__ == "__main__":
    p = write_findings()
    print(p)
