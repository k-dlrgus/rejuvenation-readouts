"""Write FINDINGS_MD2.md from results/md2/. Does not modify other FINDINGS*.md or FALSIFICATION.md."""
from __future__ import annotations

import sys
from pathlib import Path

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parent))
from md2_common import (  # noqa: E402
    MD2_DIR, MD2_SEED, MD2_BOOT, N_PERM, N_BOOT, N_RANDOM_DIR, MIN_CELLS,
    PREREG_TASK1, PREREG_TASK2, PREREG_TASK3, CLUSTER_SETTINGS, STAR_METHODS_QUOTES,
    PREREG_TASK1_FLAG, PREREG_TASK2_FLAG, PREREG_TASK3_FLAG,
    CLUSTER_SETTINGS_FLAG, DECLARED_BEFORE_SCORES_FLAG, FINDINGS_PATH, FROZEN_RULER,
    STAGE2_VERDICT_FIBRO, FIBRO3_TASK3_AGED_D0D7, FIBRO3_TASK1_SENTENCE,
    load_json, load_manifest, jsonable, progress_snapshot, LOUVAIN_RES,
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


def write_findings():
    man = load_manifest() if (MD2_DIR / "manifest.json").exists() else {}
    gs = load_json(MD2_DIR / "genesets_summary.json") if (MD2_DIR / "genesets_summary.json").exists() else {}
    t1s = load_json(MD2_DIR / "t1_summary.json") if (MD2_DIR / "t1_summary.json").exists() else {}
    t1r = load_json(MD2_DIR / "t1_reading.json") if (MD2_DIR / "t1_reading.json").exists() else {}
    t2s = load_json(MD2_DIR / "t2_summary.json") if (MD2_DIR / "t2_summary.json").exists() else {}
    t2r = load_json(MD2_DIR / "t2_reading.json") if (MD2_DIR / "t2_reading.json").exists() else {}
    t2repro = load_json(MD2_DIR / "t2_reproduction.json") if (MD2_DIR / "t2_reproduction.json").exists() else {}
    t3s = load_json(MD2_DIR / "t3_summary.json") if (MD2_DIR / "t3_summary.json").exists() else {}
    cmeta = load_json(MD2_DIR / "cluster_meta.json") if (MD2_DIR / "cluster_meta.json").exists() else {}

    t1_txt = f"Task 1: `{t1r.get('key', 'not scored')}`."
    t2_txt = f"Task 2: `{t2r.get('key', 'not scored')}`."
    t3_txt = (
        "Task 3: aging signature not run."
        if not t3s.get("aging_signature_built", True)
        else "Task 3: aging signature scored."
    )
    if t3s and t3s.get("aging_signature_built") is False:
        t3_txt = "Task 3: aging signature not run (GSE113957 not DESeq2-rebuildable); MD and frozen ruler tabulated."
    status = (
        f"**Status:** {t1_txt} {t2_txt} {t3_txt} "
        f"Fired: Task 1 `{t1r.get('key', 'NA')}`; Task 2 `{t2r.get('key', 'NA')}`. "
        f"Seed `{MD2_SEED}`. boot `{MD2_BOOT}`. n_perm={N_PERM}. n_boot={N_BOOT}. "
        f"n_random={N_RANDOM_DIR}. Frozen ruler `{FROZEN_RULER.name}` exists={FROZEN_RULER.exists()}."
    )

    n_lv_a = (cmeta.get("louvain_GM00731") or {}).get("n_clusters")
    n_lv_y = (cmeta.get("louvain_GM23815") or {}).get("n_clusters")
    fail_lines = []
    for f in man.get("failures") or []:
        msg = str(f.get("message"))
        resolved = ""
        if "GSE226189" in msg and (MD2_DIR / "t3_instruments.csv").exists():
            t3df = pd.read_csv(MD2_DIR / "t3_instruments.csv")
            hit = t3df[(t3df.instrument.astype(str) == "MD_AddModuleScore") & (t3df.accession.astype(str) == "GSE226189")]
            if len(hit) and pd.notna(hit.rho.iloc[0]):
                resolved = (
                    " Resolved on retry: geneCOUNT Tracking_ID is Ensembl; mapped to symbols via "
                    "GSE297234 10x features + frozen ruler; MD then scored (see Task 3 table)."
                )
        fail_lines.append(f"- **{f.get('step')}:** {f.get('message')}{resolved}")
    if not fail_lines:
        fail_lines = ["- none recorded"]

    cols_used = []
    for k, v in (man.get("columns") or {}).items():
        cols_used.append(f"- **{k}** source=`{v.get('source')}` columns={v.get('columns')}")
    gsets = []
    for g in man.get("gene_sets_read") or []:
        gsets.append(
            f"- **{g.get('key')}** name=`{g.get('name')}` n={g.get('n')} source=`{g.get('source')}`"
        )

    lines = [
        "# FINDINGS_MD2 — Lu et al. cell states, then both instruments on the same cells",
        "",
        status,
        "",
        "Does not modify any existing `FINDINGS_*.md` or `FALSIFICATION.md`. No GSE325735. "
        "Uncalibrated R² is reported and is never a gate. Nothing averaged across donors, clusters, or resolutions. "
        "Refit nothing. d7 is not the primary Stage 2 endpoint.",
        "",
        f"Reproduced by `src/md2_run.py`. Frozen ruler: `{FROZEN_RULER}` exists={FROZEN_RULER.exists()}.",
        "",
        f"Flag: `PREREG_TASK1.flag` exists={PREREG_TASK1_FLAG.exists()}. "
        f"`PREREG_TASK2.flag` exists={PREREG_TASK2_FLAG.exists()}. "
        f"`PREREG_TASK3.flag` exists={PREREG_TASK3_FLAG.exists()}. "
        f"`CLUSTER_SETTINGS.flag` exists={CLUSTER_SETTINGS_FLAG.exists()}. "
        f"`DECLARED_BEFORE_SCORES.flag` exists={DECLARED_BEFORE_SCORES_FLAG.exists()}.",
        "",
        "## Methods as published",
        "",
        "Lu et al., Cell 2025, 188:5895–5911; DOI 10.1016/j.cell.2025.07.031. STAR Methods quotes:",
        "",
        f"> {STAR_METHODS_QUOTES['gene_sets']}",
        "",
        f"> {STAR_METHODS_QUOTES['new_scrna']}",
        "",
        f"> {STAR_METHODS_QUOTES['public_scrna_score']}",
        "",
        f"> {STAR_METHODS_QUOTES['lognormalize_prose']}",
        "",
        f"> {STAR_METHODS_QUOTES['fig3g']}",
        "",
        f"- MSigDB version used: `{gs.get('msigdb_version', 'NA')}`",
        f"- HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION n={gs.get('n_EMT')}",
        f"- TFs added to EMT: {gs.get('MD_TFs_added')} already in EMT: {gs.get('MD_TFs_already_in_EMT')}",
        f"- MD built n={gs.get('n_MD_built')}",
        f"- HALLMARK_TGF_BETA_SIGNALING n={gs.get('n_TGFB')}",
        f"- mmc3 `MD_signatures` column `MD score` n={gs.get('mmc3_n_MD')} equal_to_built={gs.get('mmc3_MD_equal')}",
        f"- mmc3 only_built={gs.get('mmc3_MD_only_built')} only_supplement={gs.get('mmc3_MD_only_supplement')}",
        f"- mmc3 `TGFB score` n={gs.get('mmc3_n_TGFB')} equal_to_built={gs.get('mmc3_TGFB_equal')}",
        f"- mmc3 TGFB only_built={gs.get('mmc3_TGFB_only_built')} only_supplement={gs.get('mmc3_TGFB_only_supplement')}",
        f"- Mismatch is reported, not reconciled. Scoring uses the STAR Methods built set.",
        f"- AddModuleScore: Python implementation of Seurat's control-gene-bin algorithm "
        f"(nbin=24, ctrl=100, quantile bins with 1e-30 Gaussian noise). LogNormalize "
        f"scale.factor=10000 (Seurat NormalizeData default). Paper public-data prose multiplies "
        f"by median total UMIs; that prose is not used.",
        f"- Clustering Louvain: STAR Methods specify Seurat v5 sctransform v2. This environment "
        f"has no R/sctransform. Deviation (not a silent substitute): LogNormalize + quadratic "
        f"HVG n=3000 + ScaleData-style percent.mt residualization + PCA 50 + SNN k.param=20 "
        f"prune.SNN=1/15 + networkx louvain_communities resolution={LOUVAIN_RES} "
        f"(Seurat FindClusters default; STAR Methods omit resolution).",
        f"- Louvain n clusters GM00731={n_lv_a} GM23815={n_lv_y}. Not retuned.",
        f"- Aging signature (Fleischer DESeq2): built={gs.get('aging_signature_built')} "
        f"reason={gs.get('aging_signature_reason')}",
        f"- mmc3 Age up n={gs.get('mmc3_age_up_n')} Age down n={gs.get('mmc3_age_down_n')} "
        f"used_as_substitute={gs.get('mmc3_age_used_as_substitute')}",
        "",
        "## Failures recorded",
        "",
    ]
    lines += fail_lines
    lines += [
        "",
        "## Task 1 pre-registration (verbatim, written before any age score)",
        "",
        CLUSTER_SETTINGS,
        "",
        PREREG_TASK1,
        "",
        "## Task 1 — resolution sweep tables",
        "",
        _md(MD2_DIR / "t1_cluster_table.csv", [
            "resolution", "cell_line", "day", "cluster", "n_cells", "frac_timepoint",
            "below_min_cells", "age_score", "md_score", "tgfb_score",
            "pluri_primary", "pluri_with_OSKM",
        ]),
        "",
        "## Task 1 — frozen-ruler random-direction null (d0→d7 and d0→d10; not averaged)",
        "",
        _md(MD2_DIR / "t1_age_null.csv", [
            "resolution", "cluster", "endpoint", "n_cells_d0", "n_cells_end",
            "score_d0", "score_end", "decline", "p", "n_random_ge_real", "n_random",
            "pass_p", "ok", "reason", "null_kind",
        ]),
        "",
        f"k=3 age scores and nulls reused from `results/fibro3/` (not recomputed; TMM among both donors in FIBRO3). "
        f"k=8, k=15, Louvain TMM among GM00731 cluster×timepoint rows of that resolution only.",
        "",
        "## Task 1 reading",
        "",
        f"Fired key: `{t1r.get('key', 'NA')}`.",
        "",
        t1r.get("text", "_Task 1 not scored._"),
        "",
        "Louvain occupancy: most Louvain clusters are occupied at one timepoint. "
        "Aged GM00731 tests with ≥20 cells at both endpoints (from t1_reading.json): "
        f"k3 n_ok={(t1r.get('by_resolution') or {}).get('k3', {}).get('n_ok')} n_pass={(t1r.get('by_resolution') or {}).get('k3', {}).get('n_pass')}; "
        f"k8 n_ok={(t1r.get('by_resolution') or {}).get('k8', {}).get('n_ok')} n_pass={(t1r.get('by_resolution') or {}).get('k8', {}).get('n_pass')}; "
        f"k15 n_ok={(t1r.get('by_resolution') or {}).get('k15', {}).get('n_ok')} n_pass={(t1r.get('by_resolution') or {}).get('k15', {}).get('n_pass')}; "
        f"louvain n_ok={(t1r.get('by_resolution') or {}).get('louvain', {}).get('n_ok')} n_pass={(t1r.get('by_resolution') or {}).get('louvain', {}).get('n_pass')}.",
        "",
        "d7 is not the Stage 2 endpoint. Rows are not averaged across clusters or resolutions.",
        "",
        "## Task 2 pre-registration (verbatim, written before any MD vs ruler comparison)",
        "",
        PREREG_TASK2,
        "",
        "## Task 2 — Figure 3G reproduction check",
        "",
        f"matches_reported_direction={t2repro.get('matches_reported_direction')}",
        f"n_PartialReprog={t2repro.get('n_PartialReprog')} n_NonReprog={t2repro.get('n_NonReprog')}",
        f"md_mean_PartialReprog={_signed(t2repro.get('md_mean_PartialReprog'))} "
        f"md_mean_NonReprog={_signed(t2repro.get('md_mean_NonReprog'))} "
        f"delta={_signed(t2repro.get('md_Partial_minus_NonReprog'))}",
        f"n_cells_by_label={t2repro.get('n_cells_by_label')}",
        f"criterion={t2repro.get('criterion')}",
        "",
        "Cluster state labels (argmax mean AddModuleScore of mmc3 Reprog_cell_state_signatures):",
        "",
        _md(MD2_DIR / "t2_cluster_labels_GM00731.csv"),
        "",
        _md(MD2_DIR / "t2_cluster_labels_GM23815.csv"),
        "",
        "## Task 2 — both instruments on the same Louvain cluster × timepoint rows",
        "",
        _md(MD2_DIR / "t2_cluster_table.csv", [
            "resolution", "cell_line", "day", "cluster", "label", "n_cells",
            "below_min_cells", "age_score", "md_score", "tgfb_score", "pluri_primary",
        ]),
        "",
        "## Task 2 — cross-correlation (Spearman ρ, age vs MD; permutation shuffles MD)",
        "",
        _md(MD2_DIR / "t2_crosscorr.csv"),
        "",
        "## Task 2 — MD composition test (size-matched random gene-set null, not permuted weights)",
        "",
        _md(MD2_DIR / "t2_md_null.csv", [
            "cell_line", "cluster", "endpoint", "n_cells_d0", "n_cells_end",
            "md_d0", "md_end", "decline", "p", "n_random_ge_real", "n_random",
            "pass_p", "ok", "reason", "null_kind", "n_MD_genes",
        ]),
        "",
        "Null difference: MD is a fixed gene list, so the null is 200 draws of the same n genes "
        "matched on mean-expression bin, scored with AddModuleScore. The frozen ruler's null is "
        "200 permutations of its weights. Those nulls are not interchangeable.",
        "",
        "## Task 2 reading",
        "",
        f"Fired key: `{t2r.get('key', 'NA')}`. comparison_valid={t2r.get('comparison_valid')}.",
        "",
        t2r.get("text", "_Task 2 not scored._"),
        "",
    ]
    mdn = MD2_DIR / "t2_md_null.csv"
    if mdn.exists():
        mddf = pd.read_csv(mdn)
        aged_pass = mddf[(mddf.cell_line.astype(str) == "GM00731") & (mddf.ok == True) & (mddf.pass_p == True)]  # noqa: E712
        labp = MD2_DIR / "t2_cluster_labels_GM00731.csv"
        labmap = {}
        if labp.exists():
            ldf = pd.read_csv(labp)
            labmap = {int(r.cluster): str(r.label) for _, r in ldf.iterrows()}
        bits = []
        for _, r in aged_pass.iterrows():
            bits.append(
                f"GM00731 cluster {int(r.cluster)} label={labmap.get(int(r.cluster), 'NA')} "
                f"{r.endpoint} MD decline={_signed(r.decline)} "
                f"p={_signed(r.p)} n_ge={int(r.n_random_ge_real)}/{int(r.n_random)} "
                f"n_d0={int(r.n_cells_d0)} n_end={int(r.n_cells_end)} "
                f"md_d0={_signed(r.md_d0)} md_end={_signed(r.md_end)}"
            )
        if bits:
            lines += ["Aged-donor MD null passes (size-matched gene-set null): " + "; ".join(bits) + ".", ""]
        else:
            lines += ["Aged-donor MD null passes: none.", ""]
        if labp.exists() and not mddf.empty:
            ok_aged = mddf[(mddf.cell_line.astype(str) == "GM00731") & (mddf.ok == True)]  # noqa: E712
            clusters_ok = sorted({int(c) for c in ok_aged.cluster})
            labels_ok = [f"{c}:{labmap.get(c, 'NA')}" for c in clusters_ok]
            n_partial_ok = sum(1 for c in clusters_ok if labmap.get(c) == "PartialReprog")
            lines += [
                f"Aged-donor Louvain clusters pairable at both endpoints (ok=True): n={len(clusters_ok)} "
                f"({', '.join(labels_ok) if labels_ok else 'none'}). PartialReprog among them: n={n_partial_ok}.",
                "",
            ]
    cc = MD2_DIR / "t2_crosscorr.csv"
    if cc.exists():
        cdf = pd.read_csv(cc)
        bits = []
        for _, r in cdf.iterrows():
            bits.append(
                f"{r.cell_line} n={int(r.n)} ρ={_signed(r.rho)} null={_signed(r.rho_null)} "
                f"p={_signed(r.rho_p)} CI=[{_signed(r.rho_ci_lo)}, {_signed(r.rho_ci_hi)}]"
            )
        if bits:
            lines += ["Cross-correlation frozen age vs MD (cluster×timepoint ≥20 cells): " + "; ".join(bits) + ".", ""]
    lines += [
        "## Task 3 pre-registration (verbatim, written before any external-cohort MD ρ)",
        "",
        PREREG_TASK3,
        "",
        "## Task 3 — instrument validity vs donor age (report-only, no gate)",
        "",
        _md(MD2_DIR / "t3_instruments.csv", [
            "instrument", "cohort", "accession", "n", "n_donors", "rho", "rho_null",
            "rho_p", "rho_ci_lo", "rho_ci_hi", "r", "cal_r2", "r2", "platform",
            "n_MD_mapped", "n_MD_missing", "note", "reason",
        ]),
        "",
        f"Aging signature: not run. {t3s.get('aging_signature_reason', gs.get('aging_signature_reason', ''))}",
        "",
        "Frozen ruler numbers are read from `results/fibro2/`. The ruler is not refit.",
        "",
        "## What this changes in FINDINGS_FIBRO.md and FINDINGS_FIBRO3.md (quote, not edited)",
        "",
        "This file does not edit `FINDINGS_FIBRO.md`, `FINDINGS_FIBRO2.md`, `FINDINGS_FIBRO3.md`, "
        "`FINDINGS_MD.md`, or `FALSIFICATION.md`.",
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
    ]
    if t1r.get("key") == "resolution_artifact":
        lines += [
            "Task 1 fired `resolution_artifact`. That quoted FIBRO3 per-cluster sentence is superseded "
            "by this file. FINDINGS_FIBRO3.md is not edited.",
            "",
        ]
    elif t1r.get("key") == "negative_holds":
        lines += [
            "Task 1 fired `negative_holds`. The quoted FIBRO3 per-cluster sentence stands. "
            "The frozen ruler does not move in these cells at k=3, 8, 15, or Louvain.",
            "",
        ]
    else:
        lines += [
            f"Task 1 fired `{t1r.get('key')}`. See the Task 1 reading. FINDINGS_FIBRO3.md is not edited.",
            "",
        ]
    lines += [
        "FINDINGS_FIBRO3.md Task 1 cell-intrinsic sentence, quoted (d7–d10 composition, a different claim):",
        "",
        f"> {FIBRO3_TASK1_SENTENCE}",
        "",
        "That sentence is about d7→d10 reversal within joint k=3 clusters. It is not the per-cluster "
        "d0→d7/d0→d10 random-direction negative this prompt re-tests. It is not edited.",
        "",
        "FINDINGS_MD.md closed the gate because the MD score's computation was not publicly retrievable. "
        "This file opens that computation from STAR Methods and runs the instrument comparison. "
        "FINDINGS_MD.md is not edited.",
        "",
        "## Limitations",
        "",
        "1. Two donors in GSE297234.",
        "2. Our reproduction of their metric is ours, not theirs. Seurat v5 / sctransform v2 was not run; "
        "LogNormalize + residualized PCA Louvain is a reported deviation. AddModuleScore is a Python "
        "reimplementation of the published algorithm, not Seurat's C++/R object.",
        "3. Clusters are not lineage-tracked. PartialReprog / NonReprog labels are argmax of mmc3 "
        "Reprog_cell_state_signatures AddModuleScore, not Slingshot trajectories from day-0 fibroblasts.",
        "4. Cross-platform shift from GTEx bulk polyA (RNASeQCv2.4.2) to 10x 3' scRNA-seq pseudobulk.",
        "5. Frozen ruler trained on GTEx V10 cultured fibroblasts, public `AGE` 10-year bins.",
        "6. Missing overlap genes at z=0. Nothing fitted on GSE297234 or on the external matrices.",
        "7. GSE113957 is FPKM; the Fleischer DESeq2 aging signature was not rebuilt and was not substituted "
        "from mmc3 Age up/down.",
        "8. MIN_CELLS=20. Cluster×timepoint rows below that threshold are not scored.",
        "9. Seed `20260914` (permutation / random directions / clustering / AddModuleScore sampling). "
        "Bootstrap seed `20260918`. n_perm=200, n_boot=200, n_random=200.",
        "10. GSE325735 not opened. d7 is not the Stage 2 endpoint.",
        "11. Louvain resolution 0.8 is the Seurat default; STAR Methods do not state the resolution used for Figure 3.",
        "12. Louvain cluster occupancy is time-specific: on aged GM00731 only cluster 10 has ≥20 cells at both d0 and d7 "
        "(n_d0=56, n_d7=191). The MD per-cluster null therefore rests on that one pairable cluster.",
        "13. GSE226189 geneCOUNT `Tracking_ID` is Ensembl. Symbols were mapped from the GSE297234 Cell Ranger feature "
        "table plus the frozen ruler; n_MD_mapped=204 n_MD_missing=1 (PRSS2; t3_report.txt).",
        "14. CELLxGENE a19d1667 MD: n_MD_mapped=202 n_MD_missing=3 (GREM1, PRSS2, TWIST2; t3_report.txt).",
        "",
        "## Open questions",
        "",
        "1. Would the authors' sctransform v2 object change which Louvain clusters are PartialReprog vs NonReprog?",
        "2. If the Fleischer integer-count matrix and PCA old/young split were released, how does that aging "
        "signature correlate with the frozen GTEx ruler on CELLxGENE a19d1667 and GSE226189?",
        "3. Slingshot from day-0 fibroblast clusters was not run (optional, report-only).",
        "",
        "## Columns actually read",
        "",
    ]
    lines += (cols_used or ["- none yet"])
    lines += ["", "## Gene-set names actually read", ""]
    lines += (gsets or ["- none yet"])
    lines += [
        "",
        "## Files",
        "",
        "- `src/md2_common.py`, `src/md2_genesets.py`, `src/md2_score.py`, `src/md2_cluster.py`, "
        "`src/md2_task1.py`, `src/md2_task2.py`, `src/md2_task3.py`, `src/md2_findings.py`, `src/md2_run.py`",
        "- `results/md2/`",
        "- `FINDINGS_MD2.md`",
        "- `PROGRESS_MD2.md`",
        "",
    ]
    FINDINGS_PATH.write_text("\n".join(lines) + "\n", encoding="utf-8")
    progress_snapshot(
        "done",
        stop="FINDINGS_WRITTEN",
        extra=(
            "FINDINGS_MD2.md regenerated from disk after Task 1 occupancy / Task 2 MD-null "
            "and cross-corr lines / GSE226189 resolved-failure annotation / limitations 12–14. "
            "Manifest still lists the first GSE226189 0-mapped AddModuleScore failure; "
            "Task 3 later mapped Ensembl→symbol and scored MD (t3_instruments.csv). "
            "Fired: Task 1 negative_holds; Task 2 instruments_disagree. Aging signature not run."
        ),
    )
    return FINDINGS_PATH


if __name__ == "__main__":
    write_findings()
