"""Write FINDINGS_MD4.md from results/md4/. Does not modify other FINDINGS*.md or FALSIFICATION.md."""
from __future__ import annotations

import sys
from pathlib import Path

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parent))
from md4_common import (  # noqa: E402
    MD4_DIR, MD4_SEED, MD4_BOOT, N_PERM, N_BOOT, N_RANDOM_DIR, MIN_CELLS_CONTRAST,
    PREREG_TASK1, PREREG_TASK2, PREREG_TASK3, FIG3G_QUOTE,
    PREREG_TASK1_FLAG, PREREG_TASK2_FLAG, PREREG_TASK3_FLAG, DECLARED_BEFORE_SCORES_FLAG,
    FINDINGS_PATH, FROZEN_RULER, AMS_NBIN, AMS_CTRL, LOUVAIN_RES, PCA_MAHAL_K,
    MD3_TASK1_STATUS, MD3_TASK1_READING, MD3_TASK1_NUMBERS, MD3_TASK2_READING,
    MD2_NEGATIVE_HOLDS, MD2_INSTRUMENTS_DISAGREE, STAGE2_VERDICT_FIBRO,
    FIBRO3_EXTRAP_RHO, WHAT_WOULD_SETTLE, AGED_LINE, YOUNG_LINE,
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
    if not fails:
        return "None recorded.\n"
    lines = []
    for f in fails:
        lines.append(f"- **{f.get('step')}:** {f.get('message')}")
    return "\n".join(lines) + "\n"


def _status_line(t1r, t2s):
    by = (t1r or {}).get("by_donor") or {}
    bits = []
    for line in (AGED_LINE, YOUNG_LINE):
        rec = by.get(line) or {}
        key = rec.get("key", "not scored")
        tps = rec.get("timepoints") or rec.get("per_timepoint") or []
        if rec.get("per_timepoint"):
            inner = "; ".join(f"d{p['day']}: `{p['key']}`" for p in rec["per_timepoint"])
            bits.append(f"{line} {inner}")
        else:
            bits.append(f"{line}: `{key}`")
    t1_txt = "Task 1: " + "; ".join(bits) + "."
    t2_txt = "Task 2: report-only extrapolation tabulated; no gate."
    t3_txt = (
        f"Task 3: FINDINGS_MD3.md Task 1 `{MD3_TASK1_STATUS}` withdrawn as unsupported "
        "(n₀=3); superseded here. MD3 Task 2 Δρ stands."
    )
    fired = "; ".join(
        f"{line} `{((by.get(line) or {}).get('key', 'NA'))}`" for line in (AGED_LINE, YOUNG_LINE)
    )
    return (
        f"**Status:** {t1_txt} {t2_txt} {t3_txt} "
        f"Fired: {fired}. "
        f"Seed `{MD4_SEED}`. boot `{MD4_BOOT}`. n_perm={N_PERM}. n_boot={N_BOOT}. "
        f"n_random={N_RANDOM_DIR}. Frozen ruler `{FROZEN_RULER.name}` exists={FROZEN_RULER.exists()}."
    )


def write_findings():
    man = load_manifest() if (MD4_DIR / "manifest.json").exists() else {}
    gs = load_json(MD4_DIR / "genesets_summary.json") if (MD4_DIR / "genesets_summary.json").exists() else {}
    t1s = load_json(MD4_DIR / "t1_summary.json") if (MD4_DIR / "t1_summary.json").exists() else {}
    t1r = load_json(MD4_DIR / "t1_reading.json") if (MD4_DIR / "t1_reading.json").exists() else {}
    t1san = load_json(MD4_DIR / "t1_sanity.json") if (MD4_DIR / "t1_sanity.json").exists() else {}
    t2s = load_json(MD4_DIR / "t2_summary.json") if (MD4_DIR / "t2_summary.json").exists() else {}

    lines = [
        "# FINDINGS_MD4 — cell-state contrast at a shared timepoint (Figure 3G)",
        "",
        _status_line(t1r, t2s),
        "",
        "Does not modify any existing `FINDINGS_*.md` or `FALSIFICATION.md`. No GSE325735. "
        "Uncalibrated R² is reported and is never a gate. Nothing averaged across donors, "
        "states, timepoints, or instruments. Refit nothing. d7 is not the primary Stage 2 "
        "endpoint. Task 2 extrapolation is report-only and is not used to discount Task 1.",
        "",
        f"Reproduced by `src/md4_run.py`. Frozen ruler: `{FROZEN_RULER}` exists={FROZEN_RULER.exists()}.",
        "",
        f"Flag: `PREREG_TASK1.flag` exists={PREREG_TASK1_FLAG.exists()}. "
        f"`PREREG_TASK2.flag` exists={PREREG_TASK2_FLAG.exists()}. "
        f"`PREREG_TASK3.flag` exists={PREREG_TASK3_FLAG.exists()}. "
        f"`DECLARED_BEFORE_SCORES.flag` exists={DECLARED_BEFORE_SCORES_FLAG.exists()}.",
        "",
        "## Why the MD3 test was the wrong one",
        "",
        "FINDINGS_MD3.md Task 1 tested a within-state trajectory (d0→d7, d0→d10) on the "
        "pooled PartialReprog population. At day 0 that pool has n₀=3 cells on aged GM00731 "
        "(FINDINGS_MD3.md occupancy: PartialReprog n_d0=3, n_d3=5220, n_d7=578, n_d10=31). "
        "The state essentially does not exist before induction. Every number resting on that "
        "baseline (ruler p=+0.403, MD p=+0.015) is uninterpretable, and the recorded status "
        f"`{MD3_TASK1_STATUS}` is not supported by it.",
        "",
        "Lu et al. Figure 3G compares cell states against each other at timepoints where both "
        "exist: partially reprogrammed vs non-reprogrammed. FINDINGS_MD2.md already had the "
        "counts for that contrast on the aged donor (PartialReprog n=5,832, NonReprog n=5,459, "
        "pooled across timepoints). This file runs that contrast per donor, per timepoint, "
        "with no pooling.",
        "",
        "Paper quote (Lu et al., Cell 2025, 188:5895–5911):",
        "",
        f"> {FIG3G_QUOTE}",
        "",
        "## Failures recorded",
        "",
        _failures_block(man),
        "## Task 1 pre-registration (verbatim, written before any MD4 cell-level instrument score)",
        "",
        PREREG_TASK1,
        "",
        "## Task 1 — cell counts by state × timepoint (md2 Louvain labels, not re-labelled)",
        "",
        "Assignment basis: argmax of mean AddModuleScore of mmc3 `Reprog_cell_state_signatures` "
        "(Fibroblast, PartialReprog, EarlyPluripotency, Pluripotency, NonReprog). "
        "Labels from `results/md2/t2_cluster_labels_*.csv`. Not re-labelled.",
        "",
        "Pooled n cells by state (sum of clusters with that label; timepoints not pooled for testing):",
        "",
        _md(MD4_DIR / "t1_cell_counts_pooled.csv"),
        "",
        "Per state × timepoint:",
        "",
        _md(MD4_DIR / "t1_cell_counts.csv", [
            "cell_line", "day", "label", "n_cells", "n_timepoint", "frac_timepoint", "ge50",
        ]),
        "",
        f"A timepoint qualifies for a two-state contrast iff both states have ≥{MIN_CELLS_CONTRAST} cells. "
        "Timepoints are not pooled to reach the threshold.",
        "",
        _md(MD4_DIR / "t1_qualify.csv", [
            "cell_line", "day", "state_a", "state_b", "n_a", "n_b", "min_cells",
            "qualifies", "primary", "reason",
        ]),
        "",
        "## Task 1 — sanity check (FINDINGS_MD2.md pooled MD means)",
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
            "",
        ]
    else:
        lines += ["_t1_sanity.json missing._", ""]

    lines += [
        "## Task 1 — cell-level means by state × timepoint",
        "",
        "These are means of per-cell scores, not TMM among pseudobulk rows. "
        "Frozen ruler is per-cell log2-CPM (TMM nf=1). MD and TGF-β are md2 AddModuleScore.",
        "",
        _md(MD4_DIR / "t1_state_timepoint_means.csv"),
        "",
        "## Task 1 — contrasts (per donor, per timepoint, per instrument)",
        "",
        "Δ = mean(state_A) − mean(state_B) and median(state_A) − median(state_B). "
        "Permutation: 200 shuffles of the two-state label among cells in those two states "
        f"at that donor × timepoint. Bootstrap: B={N_BOOT} cells within state, seed `{MD4_BOOT}`. "
        "Effect size: Cohen's d (pooled sd, n−1), same sign as Δmean. "
        f"Primary pair: PartialReprog vs NonReprog. Context rows: Fibroblast vs NonReprog, "
        "PartialReprog vs Fibroblast. Rows with qualifies=False are not tested.",
        "",
        _md(MD4_DIR / "t1_contrasts.csv", [
            "cell_line", "day", "state_a", "state_b", "instrument", "primary", "qualifies",
            "n_a", "n_b", "mean_a", "mean_b", "delta_mean", "p_mean_lower", "p_mean_two_sided",
            "delta_mean_ci_lo", "delta_mean_ci_hi", "median_a", "median_b", "delta_median",
            "p_median_lower", "cohens_d", "beats_null_lower", "ok", "reason",
        ]),
        "",
        "Primary pair, qualifying rows only (excerpt of the table above):",
        "",
    ]
    cpath = MD4_DIR / "t1_contrasts.csv"
    if cpath.exists():
        cdf = pd.read_csv(cpath)
        hit = cdf[(cdf["primary"] == True) & (cdf["qualifies"] == True)].copy()  # noqa: E712
        keepc = [c for c in [
            "cell_line", "day", "state_a", "state_b", "instrument", "n_a", "n_b",
            "mean_a", "mean_b", "delta_mean", "p_mean_lower", "p_mean_two_sided",
            "delta_mean_ci_lo", "delta_mean_ci_hi", "delta_median", "p_median_lower",
            "cohens_d", "beats_null_lower",
        ] if c in hit.columns]
        lines += [md_table(hit[keepc]) if len(hit) else "_no qualifying primary rows._", ""]
    else:
        lines += ["_t1_contrasts.csv missing._", ""]

    lines += [
        "## Task 1 reading",
        "",
    ]
    by = (t1r or {}).get("by_donor") or {}
    if by:
        for line in (AGED_LINE, YOUNG_LINE):
            rec = by.get(line) or {}
            lines.append(f"### {line}")
            lines.append("")
            lines.append(f"Fired key: `{rec.get('key', 'NA')}`. n_qualifying={rec.get('n_qualifying')}.")
            lines.append("")
            pts = rec.get("per_timepoint") or []
            if pts:
                for p in pts:
                    lines.append(f"d{p.get('day')}: `{p.get('key')}`.")
                    lines.append("")
                    lines.append(p.get("text") or "_no text._")
                    lines.append("")
            else:
                lines.append(rec.get("text") or "_no text._")
                lines.append("")
            keys_here = {p.get("key") for p in pts} if pts else {rec.get("key")}
            if "md_lower_ruler_not" in keys_here:
                lines.append(f"What would settle it: {WHAT_WOULD_SETTLE}")
                lines.append("")
            if "both_lower" in keys_here:
                lines.append("Sentences superseded (quoted, not edited in the source files):")
                lines.append("")
                lines.append(f"> {MD3_TASK1_READING}")
                lines.append("")
                lines.append(f"> {MD2_INSTRUMENTS_DISAGREE}")
                lines.append("")
    else:
        lines += ["_t1_reading.json missing._", ""]

    lines += [
        "Δmean < 0 with p_lower ≤ 0.05 is “lower / beating the null”. Donors are not pooled. "
        "Timepoints are not pooled. No winner is declared between instruments unless the "
        "pre-registered both_lower bullet fired, in which case they agree on direction; "
        "that is not a declaration that one instrument is better.",
        "",
        "## Task 2 pre-registration (verbatim, written before any MD4 extrapolation distance or ρ)",
        "",
        PREREG_TASK2,
        "",
        "## Task 2 — extrapolation (report-only, no gate)",
        "",
        f"PCA k={PCA_MAHAL_K}. Nearest neighbour is Euclidean in GTEx fibroblast training z-space. "
        "Mahalanobis is in the leading PCA of that z-space. Missing overlap genes at z=0. "
        "TMM among that donor's state×timepoint rows, then frozen μ/σ/w. Same measures as "
        "FINDINGS_FIBRO2.md / `src/fibro2_extrap.py`.",
        "",
        f"FINDINGS_FIBRO3.md Task 2 (same measures), quoted: `{FIBRO3_EXTRAP_RHO}`",
        "",
        _md(MD4_DIR / "t2_extrap.csv", [
            "cell_line", "day", "label", "n_cells", "age_score",
            "nn_euclidean", "mahalanobis_pca", "pca_k",
        ]),
        "",
        "Spearman ρ across state × timepoint rows (shuffle ruler score, hold distance). "
        "`rho_p` is p_greater = (n_null ≥ ρ_obs + 1)/(n_null+1), not p_lower. "
        "A negative ρ with large p_greater is on the 'ruler falls with distance' side; "
        "the CI is the bootstrap on ρ.",
        "",
        _md(MD4_DIR / "t2_correlation.csv", [
            "scope", "cell_line", "distance", "n", "rho", "rho_null", "rho_p",
            "rho_ci_lo", "rho_ci_hi", "note",
        ]),
        "",
        "Direction this would push Task 1 (not used to explain Task 1 away):",
        "",
        _md(MD4_DIR / "t2_direction.csv"),
        "",
        "At Task 1 qualifying primary timepoints (from `t2_direction.csv`, not reconstructed):",
        "",
    ]
    dpath = MD4_DIR / "t2_direction.csv"
    qpath = MD4_DIR / "t1_qualify.csv"
    if dpath.exists() and qpath.exists():
        ddf = pd.read_csv(dpath)
        qdf = pd.read_csv(qpath)
        qpri = qdf[(qdf["primary"] == True) & (qdf["qualifies"] == True)]  # noqa: E712
        keys = set(zip(qpri["cell_line"].astype(str), qpri["day"].astype(int)))
        if keys:
            dhit = ddf[ddf.apply(lambda r: (str(r["cell_line"]), int(r["day"])) in keys, axis=1)]
            keepd = [c for c in [
                "cell_line", "day", "n_PR", "n_NR", "nn_PR", "nn_NR", "nn_PR_minus_NR",
                "mahal_PR", "mahal_NR", "mahal_PR_minus_NR", "age_PR", "age_NR",
                "age_PR_minus_NR", "PR_farther_nn", "would_push_ruler_PR_lower_via_nn",
            ] if c in dhit.columns]
            lines += [md_table(dhit[keepd]) if len(dhit) else "_no matching rows._", ""]
        else:
            lines += ["_no qualifying primary timepoints._", ""]
    else:
        lines += ["_t2_direction.csv or t1_qualify.csv missing._", ""]

    if t2s:
        lines += [
            f"- n_rows={t2s.get('n_rows')} pca_k={t2s.get('pca_k')} "
            f"ρ_nn_all={_signed(t2s.get('rho_nn_euclidean_all'))}",
            f"- {t2s.get('note')}",
            "",
        ]
    lines += [
        "If PartialReprog sits farther from the GTEx training distribution than NonReprog, "
        "and ruler score falls with distance, that would push the Task 1 frozen-ruler contrast "
        "toward PartialReprog reading younger. Reported alongside Task 1. Not a gate.",
        "",
        "## Task 3 — withdrawal of FINDINGS_MD3.md Task 1 (dated 2026-09-19)",
        "",
        PREREG_TASK3,
        "",
        "Quoted from FINDINGS_MD3.md, not edited there:",
        "",
        f"> {MD3_TASK1_READING}",
        "",
        f"Status line quoted: Task 1: `{MD3_TASK1_STATUS}`.",
        "",
        f"Numbers quoted: {MD3_TASK1_NUMBERS}",
        "",
        "Those numbers rest on n₀=3 PartialReprog cells at day 0. The status is withdrawn as "
        "unsupported. It is superseded by Task 1 of this file, which compares PartialReprog vs "
        "NonReprog at timepoints where both exist with ≥50 cells.",
        "",
        "FINDINGS_MD3.md Task 2 (paired Δρ vs donor age on independent fibroblast cohorts) is "
        "unaffected and stands. Quoted:",
        "",
        f"> {MD3_TASK2_READING}",
        "",
        "## What this changes in FINDINGS_MD3.md / FINDINGS_MD2.md / FINDINGS_FIBRO.md (quote, not edited)",
        "",
        "This file does not edit `FINDINGS_MD3.md`, `FINDINGS_MD2.md`, `FINDINGS_FIBRO.md`, "
        "`FINDINGS_FIBRO2.md`, `FINDINGS_FIBRO3.md`, or `FALSIFICATION.md`.",
        "",
        "FINDINGS_MD3.md Task 1 status and reading, quoted, withdrawn as unsupported by this file:",
        "",
        f"> {MD3_TASK1_READING}",
        "",
        "FINDINGS_MD3.md Task 2 Δρ reading, quoted, stands:",
        "",
        f"> {MD3_TASK2_READING}",
        "",
        "FINDINGS_MD2.md Task 1 reading, quoted, is a within-cluster d0→dT test and is not this contrast:",
        "",
        f"> {MD2_NEGATIVE_HOLDS}",
        "",
        "FINDINGS_MD2.md Task 2 reading, quoted, was also a within-cluster trajectory test "
        "(pairable Louvain clusters at both endpoints; PartialReprog among them: n=0):",
        "",
        f"> {MD2_INSTRUMENTS_DISAGREE}",
        "",
        "FINDINGS_FIBRO.md Stage 2 verdict sentence, quoted. That d0→d10 all-cell verdict stands. "
        "Nothing here promotes d7 to the primary Stage 2 result:",
        "",
        f"> {STAGE2_VERDICT_FIBRO}",
        "",
        "## Limitations",
        "",
        "1. Two donors in GSE297234.",
        "2. Clusters are not lineage-tracked. PartialReprog / NonReprog labels are argmax of "
        "mmc3 Reprog_cell_state_signatures AddModuleScore, not Slingshot trajectories from day-0 "
        "fibroblasts. States are inferred from expression, not lineage.",
        "3. A state contrast at one timepoint is not a trajectory. This file does not revive "
        "the MD3 d0→d7 / d0→d10 within-PartialReprog test.",
        "4. Louvain is a Python deviation from their R/sctransform pipeline, as recorded in "
        "FINDINGS_MD2.md (LogNormalize + quadratic HVG + percent.mt residualization + PCA 50 + "
        f"SNN + networkx louvain, resolution={LOUVAIN_RES} Seurat default; STAR Methods omit resolution).",
        "5. Our reproduction of their metrics is ours, not theirs. AddModuleScore is a Python "
        "reimplementation of the published algorithm, not Seurat's C++/R object.",
        "6. Frozen ruler trained on GTEx V10 cultured fibroblasts, public `AGE` 10-year bins. "
        "Per-cell scores use log2-CPM with TMM nf=1; TMM among 20k cells was not fit. "
        "That is a different size-factor treatment from the cluster×timepoint TMM in md2/md3.",
        "7. Cross-platform shift from GTEx bulk polyA (RNASeQCv2.4.2) to 10x 3' scRNA-seq.",
        "8. Missing overlap genes at z=0. Nothing fitted on GSE297234.",
        f"9. MIN_CELLS_CONTRAST={MIN_CELLS_CONTRAST}. Timepoints below that on either side are not tested.",
        f"10. Seed `{MD4_SEED}`. Bootstrap seed `{MD4_BOOT}`. n_perm={N_PERM}, n_boot={N_BOOT}, n_random={N_RANDOM_DIR}.",
        "11. GSE325735 not opened. d7 is not the Stage 2 endpoint.",
        "12. Table S3 Age up/down are the published lists, not a rebuilt DESeq2 analysis.",
        "13. Task 2 distances are on TMM state×timepoint pseudobulks; Task 1 ruler scores are per-cell log2-CPM. Those two ruler numbers are not interchangeable.",
        "",
        "## Open questions",
        "",
        "1. Would the authors' sctransform v2 object change which Louvain clusters are PartialReprog vs NonReprog at these timepoints?",
        "2. If cells were lineage-tracked from day-0 fibroblasts, would the same contrast hold?",
        "3. An independent fibroblast age instrument on these same labelled cells at d3/d7 (aged) and d10 (young) would settle an instrument disagreement if one remains.",
        "",
        "## Gene-list versions",
        "",
        f"- MD n={gs.get('n_MD')} TGFB n={gs.get('n_TGFB')} "
        f"Age up n={gs.get('n_age_up')} Age down n={gs.get('n_age_down')}",
        f"- mmc3 equals md2 copies: up={gs.get('mmc3_equals_md2_age_up')} down={gs.get('mmc3_equals_md2_age_down')}",
        f"- used_as_deseq2_rebuild={gs.get('used_as_deseq2_rebuild')} used_as_published_lists={gs.get('used_as_published_lists')}",
        f"- AddModuleScore nbin={AMS_NBIN} ctrl={AMS_CTRL}",
        "",
    ]
    for line in (AGED_LINE, YOUNG_LINE):
        meta_p = MD4_DIR / f"t1_s3_ams_meta_{line}.json"
        if meta_p.exists():
            meta = load_json(meta_p)
            au = meta.get("age_up") or {}
            ad = meta.get("age_down") or {}
            lines.append(
                f"- {line} Age up mapped {au.get('n_mapped')}/{au.get('n_requested')} "
                f"missing={au.get('n_missing')}; Age down mapped {ad.get('n_mapped')}/{ad.get('n_requested')} "
                f"missing={ad.get('n_missing')} (from `{meta_p.name}`)"
            )
    lines += [
        "",
        "## ID types actually read",
        "",
    ]
    idt = man.get("id_types") or {}
    if idt:
        for k, v in idt.items():
            lines.append(
                f"- **{k}** kind=`{v.get('kind')}` n={v.get('n')} "
                f"n_ensembl={v.get('n_ensembl')} n_refseq={v.get('n_refseq')} "
                f"n_symbol_like={v.get('n_symbol_like')} source=`{v.get('source')}` "
                f"examples={v.get('examples')}"
            )
        lines.append("")
    else:
        lines += ["_none yet._", ""]
    lines += [
        "## Columns actually read",
        "",
    ]
    cols = man.get("columns") or {}
    if cols:
        for k, v in cols.items():
            lines.append(f"- **{k}** source=`{v.get('source')}` columns={v.get('columns')}")
        lines.append("")
    else:
        lines += ["_none yet._", ""]
    lines += [
        "## Gene-set names actually read",
        "",
    ]
    gsets = man.get("gene_sets_read") or []
    if gsets:
        for g in gsets:
            lines.append(
                f"- **{g.get('key')}** name=`{g.get('name')}` n={g.get('n')} "
                f"id_type={g.get('id_type')} source=`{g.get('source')}`"
            )
        lines.append("")
    else:
        lines += ["_none yet._", ""]
    lines += [
        "## Files",
        "",
        "- `src/md4_common.py`, `md4_task1.py`, `md4_task2.py`, `md4_findings.py`, `md4_run.py`",
        "- `results/md4/`",
        "- `FINDINGS_MD4.md`",
        "- `PROGRESS_MD4.md`",
        "",
    ]
    FINDINGS_PATH.write_text("\n".join(lines) + "\n", encoding="utf-8")
    progress_snapshot("done", stop=man.get("status", "FINDINGS_WRITTEN"))
    return FINDINGS_PATH
