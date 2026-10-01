"""Write FINDINGS_MD5.md from results/md5/. Does not modify other FINDINGS*.md or FALSIFICATION.md."""
from __future__ import annotations

import sys
from pathlib import Path

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parent))
from md5_common import (  # noqa: E402
    MD5_DIR, MD4_DIR, MD5_SEED, MD5_BOOT, N_PERM, N_BOOT, N_RANDOM_DIR,
    MIN_CELLS_A, PREREG_TASK1, PREREG_TASK2, PREREG_TASK3,
    PREREG_TASK1_FLAG, PREREG_TASK2_FLAG, PREREG_TASK3_FLAG, DECLARED_BEFORE_SCORES_FLAG,
    FINDINGS_PATH, FROZEN_RULER, AMS_NBIN, AMS_CTRL, LOUVAIN_RES, PCA_MAHAL_K,
    MD3_TASK1_STATUS, MD3_TASK1_READING, MD3_TASK2_READING,
    STAGE2_VERDICT_FIBRO, MD4_TASK1_YOUNGER_CLAIM, MD4_TASK1_STATUS_FIRED,
    MD4_CONTEXT_GM23815_D10_FIB_NR, PAPER_NONREPROG_QUOTE, FIG3G_QUOTE,
    WHAT_WOULD_SETTLE, AGED_LINE, YOUNG_LINE, SCHEME_B_FIB_DAYS, SCHEME_B_PR_DAYS,
    load_json, load_manifest, jsonable, progress_snapshot,
)
from target_common import md_table  # noqa: E402


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


def _failures_block(man):
    fails = man.get("failures") or []
    unus = man.get("unusable") or []
    bits = []
    if fails:
        for f in fails:
            bits.append(f"- **{f.get('step')}:** {f.get('message')}")
    if unus:
        for f in unus:
            bits.append(f"- **unusable {f.get('step')}:** {f.get('message')}")
    if not bits:
        return "None recorded.\n"
    return "\n".join(bits) + "\n"


def _status_line(t1r, t2s):
    by = (t1r or {}).get("by_donor") or {}
    bits = []
    for line in (AGED_LINE, YOUNG_LINE):
        rec = by.get(line) or {}
        key = rec.get("key", "not scored")
        scheme = rec.get("scheme_used")
        if rec.get("per_timepoint"):
            inner = "; ".join(
                f"d{p['day']}: `{p['key']}`" for p in rec["per_timepoint"]
            )
            bits.append(f"{line} scheme={scheme} {inner}")
        else:
            bits.append(f"{line} scheme={scheme}: `{key}`")
    t1_txt = "Task 1: " + "; ".join(bits) + "."
    t2_txt = "Task 2: report-only QC / extrap / NonReprog − Fibroblast tabulated; no gate."
    t3_txt = "Task 3: ledger of standing/withdrawn readings."
    fired = "; ".join(
        f"{line} `{((by.get(line) or {}).get('key', 'NA'))}` scheme={((by.get(line) or {}).get('scheme_used'))}"
        for line in (AGED_LINE, YOUNG_LINE)
    )
    return (
        f"**Status:** {t1_txt} {t2_txt} {t3_txt} "
        f"Fired: {fired}. "
        f"Seed `{MD5_SEED}`. boot `{MD5_BOOT}`. n_perm={N_PERM}. n_boot={N_BOOT}. "
        f"n_random={N_RANDOM_DIR}. Frozen ruler `{FROZEN_RULER.name}` exists={FROZEN_RULER.exists()}."
    )


def _contrast_cols():
    return [
        "cell_line", "day", "state_a", "state_b", "instrument", "qualifies",
        "n_a", "n_b", "mean_a", "mean_b", "delta_mean", "p_mean_lower", "p_mean_two_sided",
        "delta_mean_ci_lo", "delta_mean_ci_hi", "median_a", "median_b", "delta_median",
        "p_median_lower", "cohens_d", "beats_null_lower", "ok", "reason",
    ]


def write_findings():
    man = load_manifest() if (MD5_DIR / "manifest.json").exists() else {}
    gs = load_json(MD5_DIR / "genesets_summary.json") if (MD5_DIR / "genesets_summary.json").exists() else {}
    t1s = load_json(MD5_DIR / "t1_summary.json") if (MD5_DIR / "t1_summary.json").exists() else {}
    t1r = load_json(MD5_DIR / "t1_reading.json") if (MD5_DIR / "t1_reading.json").exists() else {}
    t1san = load_json(MD5_DIR / "t1_sanity.json") if (MD5_DIR / "t1_sanity.json").exists() else {}
    t2s = load_json(MD5_DIR / "t2_summary.json") if (MD5_DIR / "t2_summary.json").exists() else {}

    lines = [
        "# FINDINGS_MD5 — PartialReprog vs Fibroblast (did responding cells get younger than where they started?)",
        "",
        _status_line(t1r, t2s),
        "",
        "Does not modify any existing `FINDINGS_*.md` or `FALSIFICATION.md`. No GSE325735. "
        "Uncalibrated R² is reported and is never a gate. Nothing averaged across donors, "
        "states, timepoints, or instruments. Refit nothing. d7 is not the primary Stage 2 "
        "endpoint. Task 2 is report-only and is not used to discount Task 1. "
        "Scheme (b) is the weaker design and is not primary where scheme (a) qualified.",
        "",
        f"Reproduced by `src/md5_run.py`. Frozen ruler: `{FROZEN_RULER}` exists={FROZEN_RULER.exists()}.",
        "",
        f"Flag: `PREREG_TASK1.flag` exists={PREREG_TASK1_FLAG.exists()}. "
        f"`PREREG_TASK2.flag` exists={PREREG_TASK2_FLAG.exists()}. "
        f"`PREREG_TASK3.flag` exists={PREREG_TASK3_FLAG.exists()}. "
        f"`DECLARED_BEFORE_SCORES.flag` exists={DECLARED_BEFORE_SCORES_FLAG.exists()}.",
        "",
        "## Why this contrast is the decisive one",
        "",
        "FINDINGS_MD4.md Task 1 found PartialReprog lower than NonReprog on both the frozen "
        "ruler and MD, at every qualifying timepoint, in both donors (p at the permutation floor). "
        f"Its own context rows raise the alternative named in the prompt: {MD4_CONTEXT_GM23815_D10_FIB_NR}. "
        "If starting fibroblasts already read younger than NonReprog, then PartialReprog < NonReprog "
        "may reflect NonReprog reading old rather than PartialReprog reading young.",
        "",
        "Lu et al., Cell 2025, 188:5895–5911:",
        "",
        f"> {PAPER_NONREPROG_QUOTE}",
        "",
        "The contrast that distinguishes the two readings is PartialReprog vs Fibroblast. "
        "On aged GM00731 that pair never reached md4's n≥50 threshold at any timepoint. "
        f"This file runs it at n≥{MIN_CELLS_A} (scheme a) and, as the weaker design, pooled "
        f"Fibroblast d0+d3 vs PartialReprog d3+d7 with timepoint-stratified permutation (scheme b).",
        "",
        "Paper quote (Lu et al., Figure 3G):",
        "",
        f"> {FIG3G_QUOTE}",
        "",
        "## Failures recorded",
        "",
        _failures_block(man),
        "## Task 1 pre-registration (verbatim, written before any MD5 contrast p-value)",
        "",
        PREREG_TASK1,
        "",
        "## Task 1 — cell counts by state × timepoint (md2 Louvain labels, not re-labelled)",
        "",
        "Assignment basis: argmax of mean AddModuleScore of mmc3 `Reprog_cell_state_signatures` "
        "(Fibroblast, PartialReprog, EarlyPluripotency, Pluripotency, NonReprog). "
        "Labels from `results/md2/t2_cluster_labels_*.csv`. Not re-labelled. "
        "Occupancy matches `results/md4/t1_cell_counts.csv` or this file stops.",
        "",
        "Fibroblast and PartialReprog at every timepoint (the pair under test):",
        "",
        _md(MD5_DIR / "t1_fib_pr_counts.csv"),
        "",
        "Pooled n cells by state (timepoints not pooled for scheme a):",
        "",
        _md(MD5_DIR / "t1_cell_counts_pooled.csv"),
        "",
        "Per state × timepoint:",
        "",
        _md(MD5_DIR / "t1_cell_counts.csv", [
            "cell_line", "day", "label", "n_cells", "n_timepoint", "frac_timepoint", "ge30", "ge50",
        ]),
        "",
        f"Scheme (a): a timepoint qualifies iff both Fibroblast and PartialReprog have ≥{MIN_CELLS_A} cells. "
        "Timepoints are not pooled to reach the threshold. Primary if it qualifies.",
        "",
        _md(MD5_DIR / "t1_qualify_a.csv", [
            "cell_line", "day", "state_a", "state_b", "n_a", "n_b", "min_cells",
            "qualifies", "primary", "reason",
        ]),
        "",
        "Scheme (b), weaker design: Fibroblast d0+d3 vs PartialReprog d3+d7, pooled within donor. "
        "Permutation shuffles the state label within timepoint strata. Timepoint composition of both arms:",
        "",
        _md(MD5_DIR / "t1_scheme_b_composition.csv"),
        "",
        _md(MD5_DIR / "t1_qualify_b.csv", [
            "cell_line", "state_a", "state_b", "n_a", "n_b", "min_cells", "qualifies",
            "primary", "scheme_a_qualifies_on_this_donor",
            "n_Fibroblast_d0", "n_Fibroblast_d3", "n_PartialReprog_d3", "n_PartialReprog_d7",
            "reason", "note",
        ]),
        "",
        "## Task 1 — sanity check (FINDINGS_MD2.md pooled MD means; occupancy vs md4)",
        "",
    ]
    if t1san:
        lines += [
            f"- n_PartialReprog={t1san.get('n_PartialReprog')} "
            f"(md2 json {t1san.get('md2_json_n_PartialReprog')})",
            f"- n_NonReprog={t1san.get('n_NonReprog')} "
            f"(md2 json {t1san.get('md2_json_n_NonReprog')})",
            f"- md_mean_PartialReprog={_signed(t1san.get('md_mean_PartialReprog'), 6)} "
            f"(md2 json {_signed(t1san.get('md2_json_md_mean_PartialReprog'), 6)}; "
            f"FINDINGS_MD2.md printed {t1san.get('findings_printed_PartialReprog')})",
            f"- md_mean_NonReprog={_signed(t1san.get('md_mean_NonReprog'), 6)} "
            f"(md2 json {_signed(t1san.get('md2_json_md_mean_NonReprog'), 6)}; "
            f"FINDINGS_MD2.md printed {t1san.get('findings_printed_NonReprog')})",
            f"- n_match={t1san.get('n_match')} md_match_1e12={t1san.get('md_match_1e12')} "
            f"printed_3dp_match={t1san.get('printed_3dp_match')} ok={t1san.get('ok')}",
            "- occupancy vs `results/md4/t1_cell_counts.csv`: match required or STOP",
            "",
        ]
    else:
        lines += ["_t1_sanity.json missing._", ""]

    lines += [
        "## Task 1 — scheme (a) contrasts (primary where it qualifies)",
        "",
        "Δ = mean(PartialReprog) − mean(Fibroblast). "
        "Permutation: md4_task1._one_contrast unchanged — 200 shuffles of the two-state "
        f"label among cells in those two states at that donor × timepoint. Bootstrap: B={N_BOOT} "
        f"cells within state, seed `{MD5_BOOT}` with md4 pair-offset for PartialReprog vs Fibroblast. "
        "Effect size: Cohen's d (pooled sd, n−1), same sign as Δmean. "
        f"Rows with qualifies=False are not tested. MIN_CELLS={MIN_CELLS_A}.",
        "",
        _md(MD5_DIR / "t1_contrasts_a.csv", _contrast_cols()),
        "",
        "Qualifying scheme-(a) rows only:",
        "",
    ]
    apath = MD5_DIR / "t1_contrasts_a.csv"
    if apath.exists():
        adf = pd.read_csv(apath)
        hit = adf[adf["qualifies"] == True].copy()  # noqa: E712
        keepc = [c for c in [
            "cell_line", "day", "instrument", "n_a", "n_b",
            "mean_a", "mean_b", "delta_mean", "p_mean_lower", "p_mean_two_sided",
            "delta_mean_ci_lo", "delta_mean_ci_hi", "delta_median", "p_median_lower",
            "cohens_d", "beats_null_lower",
        ] if c in hit.columns]
        lines += [md_table(hit[keepc]) if len(hit) else "_no qualifying scheme-(a) rows._", ""]
    else:
        lines += ["_t1_contrasts_a.csv missing._", ""]

    lines += [
        "## Task 1 — scheme (b) contrasts (weaker design; not primary where scheme (a) qualified)",
        "",
        f"Fibroblast cells from days {list(SCHEME_B_FIB_DAYS)} vs PartialReprog cells from days "
        f"{list(SCHEME_B_PR_DAYS)}, pooled within donor. "
        "Permutation shuffles the state label **within timepoint strata** (200 shuffles). "
        "Bootstrap resamples cells within each state independently (B=200), as in md4. "
        "This design lets timepoint composition of the two arms persist in every permutation. "
        "It is the weaker design. It is not the primary reading on any donor where scheme (a) qualified.",
        "",
        _md(MD5_DIR / "t1_contrasts_b.csv", [
            "cell_line", "instrument", "qualifies", "primary", "n_a", "n_b",
            "n_Fibroblast_d0", "n_Fibroblast_d3", "n_PartialReprog_d3", "n_PartialReprog_d7",
            "mean_a", "mean_b", "delta_mean", "p_mean_lower", "p_mean_two_sided",
            "delta_mean_ci_lo", "delta_mean_ci_hi", "delta_median", "p_median_lower",
            "cohens_d", "beats_null_lower", "n_strata_mixed", "ok", "reason",
        ]),
        "",
        "## Task 1 reading",
        "",
        "Scheme (a) is primary where it qualifies. Scheme (b) is tabulated above and is used "
        "for the fired bullet only where scheme (a) does not qualify. Donors are not pooled. "
        "Timepoints are not averaged. No winner is declared between instruments.",
        "",
    ]
    by = (t1r or {}).get("by_donor") or {}
    if by:
        for line in (AGED_LINE, YOUNG_LINE):
            rec = by.get(line) or {}
            lines.append(f"### {line}")
            lines.append("")
            lines.append(
                f"Fired key: `{rec.get('key', 'NA')}`. scheme_used=`{rec.get('scheme_used')}`. "
                f"scheme_a_is_primary={rec.get('scheme_a_is_primary')}. "
                f"n_qualifying_a={rec.get('n_qualifying_a')} n_qualifying_b={rec.get('n_qualifying_b')}."
            )
            lines.append("")
            pts = rec.get("per_timepoint") or []
            if pts:
                for p in pts:
                    lines.append(f"scheme={p.get('scheme')} d{p.get('day')}: `{p.get('key')}`.")
                    lines.append("")
                    lines.append(p.get("text") or "_no text._")
                    lines.append("")
            else:
                lines.append(rec.get("text") or "_no text._")
                lines.append("")
            if rec.get("scheme_b_tabulated_not_primary"):
                lines.append(
                    "Scheme (b) also qualifies on this donor and is tabulated; it is not the primary reading."
                )
                lines.append("")
            keys_here = {p.get("key") for p in pts} if pts else {rec.get("key")}
            if "instruments_disagree" in keys_here:
                lines.append(f"What would settle it: {WHAT_WOULD_SETTLE}")
                lines.append("")
    else:
        lines += ["_t1_reading.json missing._", ""]

    lines += [
        "Δmean < 0 with p_lower ≤ 0.05 on the frozen ruler is “PartialReprog below Fibroblast / beating the null”. "
        "Donors are not pooled. Scheme (a) timepoints are not pooled.",
        "",
        "## Task 2 pre-registration (verbatim, written before any MD5 QC / extrap / NR − Fib p-value)",
        "",
        PREREG_TASK2,
        "",
        "## Task 2 — QC (report-only, no gate)",
        "",
        f"Per donor × timepoint × {{NonReprog, Fibroblast, PartialReprog}} with n≥{MIN_CELLS_A}: "
        "median UMI, median genes, mitochondrial fraction from `results/md2/louvain_obs_*.csv`. "
        "Cell-cycle phase fractions if computable from the md2 object.",
        "",
        _md(MD5_DIR / "t2_qc.csv", [
            "cell_line", "day", "label", "n_cells", "ok",
            "median_umi", "median_genes", "median_mito_frac",
            "cell_cycle_computable", "reason",
        ]),
        "",
    ]
    cc_path = MD5_DIR / "t2_cellcycle.json"
    if cc_path.exists():
        cc = load_json(cc_path)
        lines.append("Cell-cycle from the md2 object:")
        lines.append("")
        for line in (AGED_LINE, YOUNG_LINE):
            rec = cc.get(line) or {}
            lines.append(
                f"- {line}: computable={rec.get('computable')} "
                f"obs_hits={rec.get('obs_cell_cycle_hits')} ams_hits={rec.get('ams_cell_cycle_hits')}"
            )
            if rec.get("reason"):
                lines.append(f"  - {rec.get('reason')}")
        lines.append("")
    lines += [
        "## Task 2 — extrapolation (report-only, no gate; read from md4, not recomputed)",
        "",
        f"PCA k={PCA_MAHAL_K}. Nearest neighbour is Euclidean in GTEx fibroblast training z-space. "
        "Mahalanobis is in the leading PCA of that z-space. Rows read from "
        "`results/md4/t2_extrap.csv` (TMM among that donor's state×timepoint rows, then frozen μ/σ/w). "
        "Same measures as FINDINGS_FIBRO2.md / `src/fibro2_extrap.py`. Not refit. Not re-TMM'd.",
        "",
        _md(MD5_DIR / "t2_extrap.csv", [
            "cell_line", "day", "label", "n_cells", "age_score",
            "nn_euclidean", "mahalanobis_pca", "pca_k",
        ]),
        "",
        "## Task 2 — NonReprog − Fibroblast contrasts (report-only, no gate)",
        "",
        "Δ = mean(NonReprog) − mean(Fibroblast). Same nulls and CIs as Task 1 scheme (a): "
        "md4_task1._one_contrast unchanged, within timepoint, n≥30 both states. "
        "`nr_older` = Δmean > 0 and p_two_sided ≤ 0.05.",
        "",
        _md(MD5_DIR / "t2_contrasts_nr_fib.csv", [
            "cell_line", "day", "instrument", "qualifies", "n_a", "n_b",
            "mean_a", "mean_b", "delta_mean", "p_mean_lower", "p_mean_two_sided",
            "delta_mean_ci_lo", "delta_mean_ci_hi", "delta_median", "cohens_d",
            "beats_null_lower", "nr_older", "ok", "reason",
        ]),
        "",
        _md(MD5_DIR / "t2_nr_older_summary.csv"),
        "",
    ]
    if t2s:
        lines += [
            t2s.get("note") or "_t2 note missing._",
            "",
            f"- report_only={t2s.get('report_only')} used_to_discount_task1={t2s.get('used_to_discount_task1')}",
            "",
        ]
    else:
        lines += ["_t2_summary.json missing._", ""]

    lines += [
        "## Task 3 — one-line ledger (quote, do not edit the source files)",
        "",
        PREREG_TASK3,
        "",
        "### FINDINGS_FIBRO.md Stage 2 d0→d10 (stands, all-cell)",
        "",
        "File: `FINDINGS_FIBRO.md`. Stands as written. Nothing here promotes d7 to the primary Stage 2 result.",
        "",
        f"> {STAGE2_VERDICT_FIBRO}",
        "",
        "### FINDINGS_MD3.md Task 1 (withdrawn in md4)",
        "",
        "File: `FINDINGS_MD3.md`. Withdrawn as unsupported in `FINDINGS_MD4.md` (n₀=3 PartialReprog cells at day 0). Not edited.",
        "",
        f"> {MD3_TASK1_READING}",
        "",
        f"Status quoted: Task 1: `{MD3_TASK1_STATUS}`.",
        "",
        "### FINDINGS_MD3.md Task 2 Δρ (stands)",
        "",
        "File: `FINDINGS_MD3.md`. Unaffected by md4 or by this file.",
        "",
        f"> {MD3_TASK2_READING}",
        "",
        "### FINDINGS_MD4.md Task 1 (stands, as qualified by Task 1 here)",
        "",
        "File: `FINDINGS_MD4.md`. Not edited. Qualification is this file's Task 1 reading, per donor, above.",
        "",
        f"> {MD4_TASK1_STATUS_FIRED}",
        "",
        f"> {MD4_TASK1_YOUNGER_CLAIM}",
        "",
        "## What this changes in FINDINGS_MD4.md (quote, do not edit)",
        "",
        "This file does not edit `FINDINGS_MD4.md`, `FINDINGS_MD3.md`, `FINDINGS_MD2.md`, "
        "`FINDINGS_FIBRO.md`, `FINDINGS_FIBRO2.md`, `FINDINGS_FIBRO3.md`, or `FALSIFICATION.md`.",
        "",
        "FINDINGS_MD4.md Task 1 reading, quoted:",
        "",
        f"> {MD4_TASK1_YOUNGER_CLAIM}",
        "",
    ]
    aged = (by.get(AGED_LINE) or {}) if by else {}
    young = (by.get(YOUNG_LINE) or {}) if by else {}
    aged_key = aged.get("key")
    if aged_key == "pr_below_fibroblast":
        lines += [
            "Qualification from this file, aged donor, primary scheme: PartialReprog below "
            "Fibroblast on the frozen ruler, beating its null. The md4 PartialReprog < NonReprog "
            "reading is not an artifact of NonReprog reading old on that donor. The quoted md4 "
            "sentence stands with that qualification.",
            "",
        ]
    elif aged_key == "pr_not_below_fibroblast":
        lines += [
            "Qualification from this file, aged donor, primary scheme: PartialReprog is not below "
            "Fibroblast. The md4 result is driven by NonReprog reading old, not by PartialReprog "
            "reading young. The rejuvenation reading is not supported by this data. The quoted "
            "md4 sentence must be restated in those terms. It is not edited in FINDINGS_MD4.md.",
            "",
        ]
    elif aged_key == "instruments_disagree":
        lines += [
            "Qualification from this file, aged donor, primary scheme: instruments disagree on "
            "PartialReprog vs Fibroblast. No winner is declared. The quoted md4 sentence is not "
            "restated as a resolved rejuvenation claim.",
            "",
        ]
    elif aged_key == "neither_scheme_qualifies":
        lines += [
            "Qualification from this file: neither scheme qualifies on the aged donor. The "
            "decisive test cannot be run in this dataset. The quoted md4 sentence is not "
            "resolved by a PartialReprog vs Fibroblast test on GM00731. The young donor does "
            "not stand in for the aged one.",
            "",
        ]
    elif aged_key == "timepoints_disagree_do_not_pick":
        lines += [
            "Qualification from this file, aged donor: qualifying timepoints fire different "
            "bullets; reported separately, not picked, not averaged. The quoted md4 sentence "
            "is not collapsed into a single restatement.",
            "",
        ]
    else:
        lines += [
            f"Qualification from this file, aged donor fired key=`{aged_key}`. See Task 1 reading.",
            "",
        ]
    lines += [
        f"Young donor fired key=`{young.get('key')}` scheme=`{young.get('scheme_used')}`. "
        "Not pooled with the aged donor. Not a substitute for the aged-donor test.",
        "",
        "## Limitations",
        "",
        "1. Two donors in GSE297234.",
        "2. Clusters are not lineage-tracked. PartialReprog / Fibroblast / NonReprog labels are "
        "argmax of mmc3 Reprog_cell_state_signatures AddModuleScore, not Slingshot trajectories "
        "from day-0 fibroblasts. States are inferred from expression rather than lineage.",
        "3. A state contrast is not a trajectory. This file does not revive the MD3 d0→d7 / "
        "d0→d10 within-PartialReprog test.",
        "4. Scheme (b) mixes timepoints. Permutation is stratified within timepoint, so "
        "timepoint composition of the two arms persists in every shuffle. Scheme (b) is the "
        "weaker design.",
        "5. Louvain is a Python deviation from their R/sctransform pipeline, as recorded in "
        "FINDINGS_MD2.md (LogNormalize + quadratic HVG + percent.mt residualization + PCA 50 + "
        f"SNN + networkx louvain, resolution={LOUVAIN_RES} Seurat default; STAR Methods omit resolution).",
        "6. Our reproduction of their metrics is ours, not theirs. AddModuleScore is a Python "
        "reimplementation of the published algorithm, not Seurat's C++/R object.",
        "7. Frozen ruler trained on GTEx V10 cultured fibroblasts, public `AGE` 10-year bins. "
        "Per-cell scores are the md4 log2-CPM (TMM nf=1) vectors; not refit.",
        "8. Cross-platform shift from GTEx bulk polyA (RNASeQCv2.4.2) to 10x 3' scRNA-seq.",
        "9. Missing overlap genes at z=0. Nothing fitted on GSE297234.",
        f"10. Scheme (a) MIN_CELLS={MIN_CELLS_A}. The threshold is not lowered further. Donors are not pooled.",
        f"11. Seed `{MD5_SEED}`. Bootstrap seed `{MD5_BOOT}`. n_perm={N_PERM}, n_boot={N_BOOT}, n_random={N_RANDOM_DIR}.",
        "12. GSE325735 not opened. d7 is not the Stage 2 endpoint.",
        "13. Task 2 distances are on TMM state×timepoint pseudobulks from md4; Task 1 ruler scores "
        "are per-cell log2-CPM. Those two ruler numbers are not interchangeable.",
        "14. Scheme (a) on the aged donor, if it qualifies, is at the timepoint(s) that reach n≥30 "
        "in both states. Fibroblast cells at that timepoint are expression-labelled fibroblasts, "
        "not lineage-tracked day-0 cells.",
        "15. Task 2 QC and extrapolation are not used to explain away a Task 1 result.",
        "",
        "## Open questions",
        "",
        "1. Would the authors' sctransform v2 object change which Louvain clusters are "
        "PartialReprog vs Fibroblast at these timepoints?",
        "2. If cells were lineage-tracked from day-0 fibroblasts, would the same contrast hold?",
        "3. An independent fibroblast age instrument on these same labelled cells would settle "
        "an instrument disagreement if one remains.",
        "",
        "## Gene-list versions",
        "",
    ]
    s3 = (gs.get("s3_mapped") or {}) if gs else {}
    md4sum = (gs.get("md4_genesets_summary") or {}) if gs else {}
    if md4sum:
        lines += [
            f"- MD n={md4sum.get('n_MD')} TGFB n={md4sum.get('n_TGFB')} "
            f"Age up n={md4sum.get('n_age_up')} Age down n={md4sum.get('n_age_down')}",
            f"- mmc3 equals md2 copies: up={md4sum.get('mmc3_equals_md2_age_up')} "
            f"down={md4sum.get('mmc3_equals_md2_age_down')}",
            f"- used_as_deseq2_rebuild={md4sum.get('used_as_deseq2_rebuild')} "
            f"used_as_published_lists={md4sum.get('used_as_published_lists')}",
            f"- AddModuleScore nbin={AMS_NBIN} ctrl={AMS_CTRL}",
            "",
        ]
    for line in (AGED_LINE, YOUNG_LINE):
        rec = s3.get(line) or {}
        up = rec.get("age_up") or {}
        dn = rec.get("age_down") or {}
        if up or dn:
            lines.append(
                f"- {line} Age up mapped {up.get('n_mapped')}/{up.get('n_requested')} "
                f"missing={up.get('n_missing')}; Age down mapped {dn.get('n_mapped')}/{dn.get('n_requested')} "
                f"missing={dn.get('n_missing')} (from md4 `t1_s3_ams_meta_{line}.json`)"
            )
    lines += ["", "## ID types actually read", ""]
    idt = man.get("id_types") or {}
    if idt:
        for k, v in idt.items():
            lines.append(
                f"- **{k}** kind=`{v.get('kind')}` n={v.get('n')} n_ensembl={v.get('n_ensembl')} "
                f"n_refseq={v.get('n_refseq')} n_symbol_like={v.get('n_symbol_like')} "
                f"source=`{v.get('source')}` examples={v.get('examples')}"
            )
    else:
        lines.append("- none recorded")
    lines += ["", "## Columns actually read", ""]
    cols = man.get("columns") or {}
    if cols:
        for k, v in cols.items():
            lines.append(f"- **{k}** source=`{v.get('source')}` columns={v.get('columns')}")
    else:
        lines.append("- none recorded")
    lines += ["", "## Gene-set names actually read", ""]
    gsets = man.get("gene_sets_read") or []
    if gsets:
        for g in gsets:
            lines.append(
                f"- **{g.get('key')}** name=`{g.get('name')}` n={g.get('n')} "
                f"id_type={g.get('id_type')} source=`{g.get('source')}`"
            )
    else:
        lines.append("- none recorded")
    lines += [
        "",
        "## Files",
        "",
        "- `src/md5_common.py`, `md5_task1.py`, `md5_task2.py`, `md5_findings.py`, `md5_run.py`",
        "- `results/md5/`",
        "- `FINDINGS_MD5.md`",
        "- `PROGRESS_MD5.md`",
        "",
    ]
    FINDINGS_PATH.write_text("\n".join(lines) + "\n", encoding="utf-8")
    progress_snapshot("done", stop="FINDINGS_WRITTEN")
    return FINDINGS_PATH
