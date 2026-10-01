"""Write FINDINGS_MD3.md from results/md3/. Does not modify other FINDINGS*.md or FALSIFICATION.md."""
from __future__ import annotations

import sys
from pathlib import Path

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parent))
from md3_common import (  # noqa: E402
    MD3_DIR, MD3_SEED, MD3_BOOT, N_PERM, N_BOOT, N_RANDOM_DIR, MIN_CELLS_20, MIN_CELLS_10,
    PREREG_TASK1, PREREG_TASK2, PREREG_TASK3, TABLE_S3_REASONING,
    PREREG_TASK1_FLAG, PREREG_TASK2_FLAG, PREREG_TASK3_FLAG, DECLARED_BEFORE_SCORES_FLAG,
    FINDINGS_PATH, FROZEN_RULER,
    MD2_NEGATIVE_HOLDS, MD2_INSTRUMENTS_DISAGREE, MD2_PARTIAL_ZERO,
    STAGE2_VERDICT_FIBRO, FIBRO3_TASK3_AGED_D0D7, STATE_ORDINAL,
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


def write_findings():
    man = load_manifest() if (MD3_DIR / "manifest.json").exists() else {}
    gs = load_json(MD3_DIR / "genesets_summary.json") if (MD3_DIR / "genesets_summary.json").exists() else {}
    t1s = load_json(MD3_DIR / "t1_summary.json") if (MD3_DIR / "t1_summary.json").exists() else {}
    t1r = load_json(MD3_DIR / "t1_reading.json") if (MD3_DIR / "t1_reading.json").exists() else {}
    t2s = load_json(MD3_DIR / "t2_summary.json") if (MD3_DIR / "t2_summary.json").exists() else {}
    t2r = load_json(MD3_DIR / "t2_reading.json") if (MD3_DIR / "t2_reading.json").exists() else {}
    t3s = load_json(MD3_DIR / "t3_summary.json") if (MD3_DIR / "t3_summary.json").exists() else {}

    t1_txt = f"Task 1: `{t1r.get('key', 'not scored')}`."
    t2_txt = f"Task 2: `{t2r.get('key', 'not scored')}`."
    t3_txt = "Task 3: report-only MD co-variation tabulated; no gate."
    extra = ""
    if t1r.get("levels_disagree"):
        extra = " Task 1 also fired pooled vs per-cluster disagreement (both reported, not picked)."
    status = (
        f"**Status:** {t1_txt} {t2_txt} {t3_txt}{extra} "
        f"Fired: Task 1 `{t1r.get('key', 'NA')}`; Task 2 `{t2r.get('key', 'NA')}`. "
        f"Seed `{MD3_SEED}`. boot `{MD3_BOOT}`. n_perm={N_PERM}. n_boot={N_BOOT}. "
        f"n_random={N_RANDOM_DIR}. Frozen ruler `{FROZEN_RULER.name}` exists={FROZEN_RULER.exists()}."
    )

    lines = [
        "# FINDINGS_MD3 — test the ruler where their claim lives; compare instruments against donor age",
        "",
        status,
        "",
        "Does not modify any existing `FINDINGS_*.md` or `FALSIFICATION.md`. No GSE325735. "
        "Uncalibrated R² is reported and is never a gate. Nothing averaged across donors, clusters, "
        "or instruments. Refit nothing. d7 is not the primary Stage 2 endpoint. "
        "≥10-cell threshold is a sensitivity check, not primary. Paired Δρ is the age-instrument "
        "difference; overlapping separate CIs are not.",
        "",
        f"Reproduced by `src/md3_run.py`. Frozen ruler: `{FROZEN_RULER}` exists={FROZEN_RULER.exists()}.",
        "",
        f"Flag: `PREREG_TASK1.flag` exists={PREREG_TASK1_FLAG.exists()}. "
        f"`PREREG_TASK2.flag` exists={PREREG_TASK2_FLAG.exists()}. "
        f"`PREREG_TASK3.flag` exists={PREREG_TASK3_FLAG.exists()}. "
        f"`DECLARED_BEFORE_SCORES.flag` exists={DECLARED_BEFORE_SCORES_FLAG.exists()}.",
        "",
        TABLE_S3_REASONING,
        "",
        f"- mmc3 Age up n={gs.get('n_age_up')} id_type={gs.get('age_up_id_type')} "
        f"Age down n={gs.get('n_age_down')} id_type={gs.get('age_down_id_type')}",
        f"- mmc3 equals md2 genesets.json copies: up={gs.get('mmc3_equals_md2_age_up')} "
        f"down={gs.get('mmc3_equals_md2_age_down')}",
        f"- used_as_deseq2_rebuild={gs.get('used_as_deseq2_rebuild')} "
        f"used_as_published_lists={gs.get('used_as_published_lists')}",
        "",
        "## Failures recorded",
        "",
        _failures_block(man),
        "## Task 1 pre-registration (verbatim, written before any MD3 instrument score)",
        "",
        PREREG_TASK1,
        "",
        "## Task 1 — cluster labels (md2 assignment, reused)",
        "",
        "Assignment basis: argmax of mean AddModuleScore of mmc3 `Reprog_cell_state_signatures` "
        "(Fibroblast, PartialReprog, EarlyPluripotency, Pluripotency, NonReprog). "
        "Labels from `results/md2/t2_cluster_labels_*.csv`. Not re-labelled. Prompt's EarlyPluri "
        "is this EarlyPluripotency column.",
        "",
        _md(MD3_DIR / "t1_occupancy.csv", [
            "cell_line", "cluster", "label", "n_cells_all", "n_d0", "n_d3", "n_d7", "n_d10",
            "mean_Fibroblast", "mean_PartialReprog", "mean_EarlyPluripotency",
            "mean_Pluripotency", "mean_NonReprog",
        ]),
        "",
        "Pooled n cells by state (sum of clusters with that label):",
        "",
        _md(MD3_DIR / "t1_occupancy_pooled.csv"),
        "",
        "## Task 1 — three aggregation levels (not averaged)",
        "",
        "### Primary cell excerpt — aged GM00731 pooled PartialReprog (not a substitute for the tables below)",
        "",
    ]
    ppath = MD3_DIR / "t1_pooled_scores.csv"
    npath = MD3_DIR / "t1_pooled_null.csv"
    if ppath.exists():
        pdf = pd.read_csv(ppath)
        hit = pdf[(pdf.cell_line.astype(str) == "GM00731") & (pdf.label.astype(str) == "PartialReprog")]
        keep = [c for c in [
            "cell_line", "label", "day", "n_cells", "age_score", "md_score", "tgfb_score",
            "age_up", "age_down", "pluri_primary",
        ] if c in hit.columns]
        lines += [md_table(hit[keep]) if len(hit) else "_no pooled PartialReprog rows._", ""]
    else:
        lines += ["_t1_pooled_scores.csv missing._", ""]
    if npath.exists():
        ndf = pd.read_csv(npath)
        hitn = ndf[(ndf.cell_line.astype(str) == "GM00731") & (ndf.label.astype(str) == "PartialReprog")
                   & (ndf.instrument.isin(["frozen_ruler", "md_score"]))]
        keepn = [c for c in [
            "instrument", "endpoint", "n_cells_d0", "n_cells_end", "score_d0", "score_end",
            "decline", "p", "n_random_ge_real", "pass_p", "ok", "reason", "null_kind",
        ] if c in hitn.columns]
        lines += [
            "Aged pooled PartialReprog ruler and MD nulls:",
            "",
            md_table(hitn[keepn]) if len(hitn) else "_no pooled PartialReprog null rows._",
            "",
        ]
    lines += [
        "### (a) pooled by cell state — primary cell (Figure 3G level)",
        "",
        _md(MD3_DIR / "t1_pooled_scores.csv", [
            "cell_line", "label", "day", "n_cells", "age_score", "md_score", "tgfb_score",
            "age_up", "age_down", "pluri_primary", "below_min_cells",
        ]),
        "",
        "Pooled endpoint nulls:",
        "",
        _md(MD3_DIR / "t1_pooled_null.csv", [
            "cell_line", "label", "instrument", "endpoint", "n_cells_d0", "n_cells_end",
            "score_d0", "score_end", "decline", "p", "n_random_ge_real", "pass_p", "ok", "reason", "null_kind",
        ]),
        "",
        "### (b) per cluster at ≥20 cells (FINDINGS_MD2.md threshold, unchanged)",
        "",
        _md(MD3_DIR / "t1_cluster20_scores.csv", [
            "cell_line", "cluster", "label", "day", "n_cells", "below_min_cells",
            "age_score", "md_score", "tgfb_score", "age_up", "age_down", "pluri_primary",
        ]),
        "",
        _md(MD3_DIR / "t1_cluster20_null.csv", [
            "cell_line", "cluster", "label", "instrument", "endpoint", "n_cells_d0", "n_cells_end",
            "score_d0", "score_end", "decline", "p", "n_random_ge_real", "pass_p", "ok", "reason", "null_kind",
        ]),
        "",
        "### (c) per cluster at ≥10 cells — sensitivity check, not primary",
        "",
        _md(MD3_DIR / "t1_cluster10_scores.csv", [
            "cell_line", "cluster", "label", "day", "n_cells", "below_min_cells",
            "age_score", "md_score", "tgfb_score", "age_up", "age_down", "pluri_primary",
        ]),
        "",
        _md(MD3_DIR / "t1_cluster10_null.csv", [
            "cell_line", "cluster", "label", "instrument", "endpoint", "n_cells_d0", "n_cells_end",
            "score_d0", "score_end", "decline", "p", "n_random_ge_real", "pass_p", "ok", "reason", "null_kind",
        ]),
        "",
        "MD cluster means vs `results/md2/t2_cluster_table.csv` (check, not a gate):",
        "",
        _md(MD3_DIR / "t1_md_check_vs_md2.csv"),
        "",
        "## Task 1 reading",
        "",
        f"Fired key: `{t1r.get('key', 'NA')}`.",
        "",
        t1r.get("text", "not scored"),
        "",
    ]
    if t1r.get("extra"):
        lines += [t1r["extra"], ""]
    pr = t1r.get("pooled_PR_ruler") or {}
    pm = t1r.get("pooled_PR_md") or {}
    lines += [
        f"Aged GM00731 pooled PartialReprog ruler: n_ok={pr.get('n_ok')} n_pass={pr.get('n_pass')} hits={pr.get('hits')}.",
        f"Aged GM00731 pooled PartialReprog MD: n_ok={pm.get('n_ok')} n_pass={pm.get('n_pass')} hits={pm.get('hits')}.",
        "",
        "d7 is not the Stage 2 endpoint. Rows are not averaged across donors, clusters, or instruments. "
        "≥10 is not primary.",
        "",
        "## Task 2 pre-registration (verbatim, written before any three-way ρ or Δρ)",
        "",
        PREREG_TASK2,
        "",
        "## Task 2 — instrument × cohort (ρ vs donor age)",
        "",
        _md(MD3_DIR / "t2_instruments.csv", [
            "instrument", "cohort", "accession", "n", "n_donors", "rho", "rho_null", "rho_p",
            "rho_ci_lo", "rho_ci_hi", "r", "cal_r2", "r2", "platform", "id_type", "transform",
            "n_MD_mapped", "n_age_up_mapped", "n_age_down_mapped", "n_overlap", "n_excluded_child",
            "weaker_test", "reason",
        ]),
        "",
        "## Task 2 — paired bootstrap Δρ (ρ_ruler − ρ_other; this is the claim)",
        "",
        _md(MD3_DIR / "t2_delta_rho.csv", [
            "cohort", "other", "n_donors", "rho_ruler", "rho_other", "delta_rho",
            "delta_ci_lo", "delta_ci_hi", "excludes_zero", "favours", "note",
        ]),
        "",
        "## Task 2 reading",
        "",
        f"Fired key: `{t2r.get('key', 'NA')}`.",
        "",
        t2r.get("text", "not scored"),
        "",
        "GSE113957 MD and Table S3 AddModuleScore on ranks: not run if `t2_instruments.csv` "
        "records a reason for those instruments. nbin=24 was not reduced. The frozen-ruler "
        "rank ρ is still reported. That cohort does not contribute paired Δρ.",
        "",
        "## Task 3 pre-registration (verbatim, written before any MD co-variation ρ)",
        "",
        PREREG_TASK3,
        "",
        f"Ordinal used: {STATE_ORDINAL}. NonReprog excluded.",
        "",
        "## Task 3 — MD co-variation (report-only, no gate)",
        "",
        _md(MD3_DIR / "t3_md_covariation.csv", [
            "cell_line", "other", "n", "rho", "rho_null", "rho_p", "rho_ci_lo", "rho_ci_hi",
            "n_excluded_NonReprog", "note",
        ]),
        "",
        "MD vs frozen ruler vs FINDINGS_MD2.md (check):",
        "",
        _md(MD3_DIR / "t3_md_vs_ruler_check.csv"),
        "",
        "## What this changes in FINDINGS_MD2.md / FINDINGS_FIBRO3.md / FINDINGS_FIBRO.md (quote, not edited)",
        "",
        "This file does not edit `FINDINGS_MD2.md`, `FINDINGS_FIBRO.md`, `FINDINGS_FIBRO2.md`, "
        "`FINDINGS_FIBRO3.md`, or `FALSIFICATION.md`.",
        "",
        "FINDINGS_MD2.md Task 1 reading, quoted:",
        "",
        f"> {MD2_NEGATIVE_HOLDS}",
        "",
        "FINDINGS_MD2.md Task 2 reading, quoted:",
        "",
        f"> {MD2_INSTRUMENTS_DISAGREE}",
        "",
        "FINDINGS_MD2.md pairable-cluster sentence, quoted:",
        "",
        f"> {MD2_PARTIAL_ZERO}",
        "",
        "FINDINGS_FIBRO.md Stage 2 verdict sentence, quoted:",
        "",
        f"> {STAGE2_VERDICT_FIBRO}",
        "",
        "That d0→d10 all-cell verdict stands. Nothing here promotes d7 to the primary Stage 2 result.",
        "",
        "FINDINGS_FIBRO3.md per-cluster d0→d7 sentence, quoted:",
        "",
        f"> {FIBRO3_TASK3_AGED_D0D7}",
        "",
        "Task 1 of this file tests the pooled PartialReprog state, which those sentences did not. "
        "The quoted sentences are not edited there.",
        "",
        "## Limitations",
        "",
        "1. Two donors in GSE297234.",
        "2. Clusters are not lineage-tracked. PartialReprog / NonReprog labels are argmax of mmc3 "
        "Reprog_cell_state_signatures AddModuleScore, not Slingshot trajectories from day-0 fibroblasts. "
        "Pooled PartialReprog at d0 is the cells that sit in clusters whose all-timepoint argmax is PartialReprog, "
        "not a lineage from a d0 fibroblast cluster.",
        "3. Our reproduction of their metrics is ours, not theirs. AddModuleScore is a Python reimplementation "
        "of the published algorithm, not Seurat's C++/R object.",
        "4. Louvain is a Python deviation from their R/sctransform pipeline, as recorded in FINDINGS_MD2.md "
        "(LogNormalize + quadratic HVG + percent.mt residualization + PCA 50 + SNN + networkx louvain, "
        "resolution=0.8 Seurat default; STAR Methods omit resolution).",
        "5. GSE113957 uses a rank transform. That is a different transform from the GTEx TMM/log2-CPM spec "
        "and is a weaker test. Frozen mu/sd (log2-CPM units) are applied as-is to ranks; the ruler is not refit. "
        "AddModuleScore nbin=24 on those ranks collapsed to 22 bins; nbin was not reduced; MD and Table S3 "
        "AMS on this cohort were not scored. Annotation/Divergence first tokens were unique per transcript "
        "(27142 transcripts → 27142 symbols); isoforms were not summed under a shared gene symbol. "
        "Adult ≥18 n=111; excluded missing age n=8; excluded age<18 n=24.",
        "6. Frozen ruler trained on GTEx V10 cultured fibroblasts, public `AGE` 10-year bins.",
        "7. Cross-platform shift from GTEx bulk polyA (RNASeQCv2.4.2) to 10x 3' scRNA-seq pseudobulk and to bulk FPKM ranks.",
        "8. Missing overlap genes at z=0. Nothing fitted on GSE297234 or on the external matrices.",
        "9. MIN_CELLS=20 is primary for per-cluster. ≥10 is a labelled sensitivity check. Pooled state has no extra cutoff; n is reported.",
        "10. Seed `20260914`. Bootstrap seed `20260918`. n_perm=200, n_boot=200, n_random=200.",
        "11. GSE325735 not opened. d7 is not the Stage 2 endpoint.",
        "12. Table S3 Age up/down are the published lists, not a rebuilt DESeq2 analysis.",
        "13. Pluripotency−fibroblast null is size-matched in frozen-ruler Z space, not AMS expression-bin matching.",
        "",
        "## Open questions",
        "",
        "1. Would the authors' sctransform v2 object put more PartialReprog cells at d0, so the pooled d0→d7 ruler test is not a 3-cell versus 578-cell comparison?",
        "2. If GSE113957 integer counts and the PCA old/young split were released, how would a rebuilt DE list compare to Table S3 on these cohorts?",
        "3. Slingshot from day-0 fibroblast clusters was not run.",
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
        "- `src/md3_common.py`, `md3_idtype.py`, `md3_genesets.py`, `md3_task1.py`, "
        "`md3_task2.py`, `md3_task3.py`, `md3_findings.py`, `md3_run.py`",
        "- `results/md3/`",
        "- `FINDINGS_MD3.md`",
        "- `PROGRESS_MD3.md`",
        "",
    ]
    FINDINGS_PATH.write_text("\n".join(lines) + "\n", encoding="utf-8")
    progress_snapshot("done", stop="FINDINGS_WRITTEN")
    return FINDINGS_PATH
