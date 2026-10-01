"""Write FINDINGS_FIBRO3.md from results/fibro3/. Does not modify other FINDINGS*.md."""
from __future__ import annotations

import sys
from pathlib import Path

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parent))
from fibro3_common import (  # noqa: E402
    FIBRO3_DIR, FIBRO3_SEED, FIBRO3_BOOT, N_PERM, N_BOOT, N_RANDOM_DIR,
    MIN_CELLS, CLUSTER_K, CLUSTER_NPC, CLUSTER_BATCH,
    PREREG_TASK1, PREREG_TASK4, PREREG_TASK2, PREREG_TASK3, CLUSTER_SETTINGS,
    PREREG_TASK1_FLAG, PREREG_TASK4_FLAG, CLUSTER_SETTINGS_FLAG,
    DECLARED_BEFORE_SCORES_FLAG, FINDINGS_PATH, FROZEN_RULER,
    STANDING_CONSTRAINT, STAGE2_VERDICT_SENTENCE, FIBRO2_TA_SENTENCE, FIBRO2_TB_SENTENCE,
    CXG_DATASET_ID, load_json, load_manifest, jsonable, progress_snapshot, fmt_ci,
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


def _status_line(man):
    t4 = man.get("t4_reading") or {}
    t1 = man.get("t1_reading") or {}
    t2_notes = man.get("t2_observations") or []
    t4s = t4.get("key", "Task 4 not scored")
    if t4.get("interpretable") is True:
        t4_txt = f"Task 4: superseded primary passed (`{t4s}`); Stage 2 trajectory interpretable as a statement about the cells."
    elif t4.get("interpretable") is False:
        t4_txt = f"Task 4: superseded primary failed (`{t4s}`); Tasks 1–3 descriptive only."
    else:
        t4_txt = f"Task 4: `{t4s}`."
    t1k = t1.get("key", "Task 1 not scored")
    t1_txt = f"Task 1: `{t1k}`."
    if t2_notes:
        t2_txt = "Task 2: " + t2_notes[0][:180]
    else:
        t2_txt = "Task 2: not scored."
    t3p = FIBRO3_DIR / "t3_summary.json"
    if t3p.exists():
        t3 = load_json(t3p)
        t3_txt = f"Task 3: n_ok={t3.get('n_ok')} / n_rows={t3.get('n_rows')} (d7 not the Stage 2 endpoint)."
    else:
        t3_txt = "Task 3: not scored."
    return (
        f"**Status:** {t4_txt} {t1_txt} {t2_txt} {t3_txt} "
        f"Seed `{FIBRO3_SEED}`. boot `{FIBRO3_BOOT}`. n_perm={N_PERM}. n_boot={N_BOOT}. "
        f"n_random={N_RANDOM_DIR}. Frozen ruler `{FROZEN_RULER.name}` exists={FROZEN_RULER.exists()}."
    )


def write_findings():
    man = load_manifest() if (FIBRO3_DIR / "manifest.json").exists() else {}
    t4s = load_json(FIBRO3_DIR / "t4_summary.json") if (FIBRO3_DIR / "t4_summary.json").exists() else {}
    t1s = load_json(FIBRO3_DIR / "t1_summary.json") if (FIBRO3_DIR / "t1_summary.json").exists() else {}
    t2s = load_json(FIBRO3_DIR / "t2_summary.json") if (FIBRO3_DIR / "t2_summary.json").exists() else {}
    t3s = load_json(FIBRO3_DIR / "t3_summary.json") if (FIBRO3_DIR / "t3_summary.json").exists() else {}
    t4_read = (t4s.get("reading") or man.get("t4_reading") or {})
    t1_read = (t1s.get("reading") or man.get("t1_reading") or {})
    t4_res = t4s.get("result") or {}

    lines = [
        "# FINDINGS_FIBRO3 — what happens between day 7 and day 10?",
        "",
        _status_line(man),
        "",
        "Does not modify any existing `FINDINGS_*.md` or `FALSIFICATION.md`. No GSE325735. "
        "Uncalibrated R² is reported and is never a gate. Nothing averaged across donors, clusters, or regimes. "
        "Refit nothing. d7 is not the primary Stage 2 endpoint.",
        "",
        f"Reproduced by `src/fibro3_run.py`. Frozen ruler: `{FROZEN_RULER}` exists={FROZEN_RULER.exists()}.",
        "",
        f"Flag: `results/fibro3/PREREG_TASK4.flag` exists={PREREG_TASK4_FLAG.exists()}. "
        f"`PREREG_TASK1.flag` exists={PREREG_TASK1_FLAG.exists()}. "
        f"`CLUSTER_SETTINGS.flag` exists={CLUSTER_SETTINGS_FLAG.exists()}. "
        f"`DECLARED_BEFORE_SCORES.flag` exists={DECLARED_BEFORE_SCORES_FLAG.exists()}.",
        "",
        "## Standing constraint",
        "",
        STANDING_CONSTRAINT,
        "",
        "## Task 4 supersession (verbatim, written before any Task 4 ρ)",
        "",
        PREREG_TASK4,
        "",
        "## Task 4 — result",
        "",
    ]
    if t4_res:
        lines += [
            f"choice=`{t4_res.get('choice')}` option_a_run={t4_res.get('option_a_run')} "
            f"accession=`{t4_res.get('accession')}` role=`{t4_res.get('role')}` "
            f"platform_detail=`{t4_res.get('platform_detail')}` used_10x_only={t4_res.get('used_10x_only')}.",
            "",
            f"n_samples={t4_res.get('n_samples')} n_donors={t4_res.get('n_donors')} "
            f"n_excluded_child={t4_res.get('n_excluded_child')} "
            f"n_overlap={t4_res.get('n_overlap')} n_missing_z0={t4_res.get('n_missing_z0')} "
            f"age_source={t4_res.get('age_source')} "
            f"age range {t4_res.get('age_min')}–{t4_res.get('age_max')} "
            f"n_unique_age={t4_res.get('n_unique_age')}.",
            "",
            f"ρ={_signed(t4_res.get('rho'))} null={_signed(t4_res.get('rho_null'))} "
            f"p={_signed(t4_res.get('rho_p'))} CI={fmt_ci(t4_res.get('rho_ci_lo'), t4_res.get('rho_ci_hi'))} "
            f"calR²={_signed(t4_res.get('cal_r2'))} R²={_signed(t4_res.get('r2'))} r={_signed(t4_res.get('r'))} "
            f"pass_rho={t4_res.get('pass_rho')} n_perm={t4_res.get('n_perm')} n_boot={t4_res.get('n_boot')}.",
            "",
            f"Fired key: `{t4_read.get('key')}`. interpretable={t4_read.get('interpretable')}.",
            "",
            t4_read.get("text", ""),
            "",
            _md(FIBRO3_DIR / "t4_result.csv", [
                "accession", "role", "choice", "tag", "n_samples", "n_donors",
                "n_overlap", "n_missing_z0", "rho", "rho_null", "rho_p",
                "rho_ci_lo", "rho_ci_hi", "cal_r2", "r2", "r", "pass_rho",
                "used_10x_only", "platform_detail",
            ]),
            "",
        ]
    else:
        lines += ["_Task 4 not scored._", ""]

    lines += [
        "## Task 1 pre-registration (verbatim, written before any joint-cluster age score)",
        "",
        CLUSTER_SETTINGS,
        "",
        PREREG_TASK1,
        "",
        "## Task 1 — clustering",
        "",
    ]
    meta_p = FIBRO3_DIR / "t1_cluster_meta.json"
    if meta_p.exists():
        meta = load_json(meta_p)
        lines.append(
            f"method=PCA svd_solver=randomized k={CLUSTER_K} n_components={CLUSTER_NPC} "
            f"random_state=`{FIBRO3_SEED}` n_init=10. "
            f"MIN_CELLS={MIN_CELLS}."
        )
        lines.append("")
        for line, rec in meta.items():
            lines.append(
                f"- {line}: n_cells={rec.get('n_cells')} n_genes_kept={rec.get('n_genes_kept')} "
                f"npc={rec.get('npc')} k={rec.get('k')} sizes={rec.get('sizes')} "
                f"days={rec.get('days')} mito_rule=`{rec.get('mito_rule')}` "
                f"n_mito_genes={rec.get('n_mito_genes')}"
            )
        lines.append("")
    else:
        lines += ["_Task 1 clustering not run._", ""]

    lines += [
        "## Task 1 — per donor × cluster × timepoint",
        "",
        _md(FIBRO3_DIR / "t1_cluster_table.csv", [
            "cell_line", "day", "cluster", "gsm", "age_years", "n_cells",
            "frac_timepoint", "median_UMI", "median_genes", "mito_frac_median_cell",
            "below_min_cells", "age_score", "pluri_with_OSKM", "pluri_without_OSKM",
            "pluri_primary",
        ]),
        "",
        "## Task 1 — composition d7 vs d10",
        "",
        _md(FIBRO3_DIR / "t1_composition_d7_d10.csv", [
            "cell_line", "cluster", "n_d7", "frac_d7", "n_d10", "frac_d10",
            "delta_frac", "status", "age_score_d7", "age_score_d10",
            "below_min_d7", "below_min_d10",
        ]),
        "",
        "## Task 1 — reading",
        "",
    ]
    if t1_read:
        lines += [
            f"Fired key: `{t1_read.get('key')}`. donor={t1_read.get('donor')} "
            f"MIN_CELLS={t1_read.get('min_cells')} C_ok={t1_read.get('C_ok')} "
            f"mix_shift={t1_read.get('mix_shift')}.",
            "",
            t1_read.get("text", ""),
            "",
        ]
        if t1_read.get("per_cluster"):
            lines.append("C_ok cluster d7→d10 (aged donor; not averaged):")
            lines.append("")
            lines.append(md_table(pd.DataFrame(t1_read["per_cluster"])))
            lines.append("")
        if t1_read.get("mix_notes"):
            lines.append("mix_notes: " + "; ".join(str(x) for x in t1_read["mix_notes"]))
            lines.append("")
        lines.append(
            "d7 is not the pre-registered Stage 2 endpoint. This reading characterizes "
            "the d7–d10 trajectory. It does not re-decide FINDINGS_FIBRO.md."
        )
        lines.append("")
    else:
        lines += ["_Task 1 not scored._", ""]

    lines += [
        "## Task 2 pre-registration (verbatim, written before any distance)",
        "",
        PREREG_TASK2,
        "",
        "## Task 2 — extrapolation (report-only, no gate)",
        "",
        _md(FIBRO3_DIR / "t2_extrap.csv", [
            "cell_line", "day", "cluster", "n_cells", "below_min_cells",
            "age_score", "nn_euclidean", "mahalanobis_pca", "pca_k",
        ]),
        "",
        "## Task 2 — correlation",
        "",
        _md(FIBRO3_DIR / "t2_correlation.csv", [
            "scope", "cell_line", "distance", "n", "rho", "rho_null", "rho_p",
            "rho_ci_lo", "rho_ci_hi", "r", "r2", "cal_r2", "age_rises_with_distance",
        ]),
        "",
    ]
    if t2s.get("observations"):
        for n in t2s["observations"]:
            lines.append(f"- {n}")
        lines.append("")
    elif not (FIBRO3_DIR / "t2_summary.json").exists():
        lines += ["_Task 2 not scored._", ""]

    lines += [
        "## Task 3 pre-registration (verbatim, written before any per-cluster p)",
        "",
        PREREG_TASK3,
        "",
        "## Task 3 — per-cluster random-direction null",
        "",
        _md(FIBRO3_DIR / "t3_null.csv", [
            "cell_line", "cluster", "endpoint", "primary", "n_cells_d0", "n_cells_end",
            "score_d0", "score_end", "decline", "p", "n_random_ge_real", "n_random",
            "pass_p", "ok", "reason",
        ]),
        "",
        "d7 is not the Stage 2 endpoint. Rows are not averaged across clusters or donors.",
        "",
    ]
    t3csv = FIBRO3_DIR / "t3_null.csv"
    if t3csv.exists():
        t3df = pd.read_csv(t3csv)
        aged7 = t3df[(t3df.cell_line.astype(str) == "GM00731") & (t3df.endpoint.astype(str) == "d0→d7")]
        bits = []
        for _, r in aged7.iterrows():
            if bool(r.ok):
                bits.append(
                    f"cluster {int(r.cluster)} decline={_signed(r.decline)} "
                    f"p={_signed(r.p)} n_ge={int(r.n_random_ge_real)}/{int(r.n_random)} "
                    f"n_d0={int(r.n_cells_d0)} n_d7={int(r.n_cells_end)}"
                )
            else:
                bits.append(f"cluster {int(r.cluster)} skipped ({r.reason})")
        if bits:
            lines += [
                "Aged GM00731 d0→d7 (the all-cell T-B cell had p=+0.005): " + "; ".join(bits) + ".",
                "",
            ]

    lines += [
        "## What this changes in FINDINGS_FIBRO.md and FINDINGS_FIBRO2.md (quote, not edited)",
        "",
        "This file does not edit `FINDINGS_FIBRO.md` or `FINDINGS_FIBRO2.md`.",
        "",
        "FINDINGS_FIBRO.md Stage 2 verdict sentence, quoted:",
        "",
        f"> {STAGE2_VERDICT_SENTENCE}",
        "",
        "That d0→d10 all-cell verdict stands. Nothing here promotes d7 to the primary Stage 2 result.",
        "",
        "FINDINGS_FIBRO2.md T-A sentence, quoted:",
        "",
        f"> {FIBRO2_TA_SENTENCE}",
        "",
    ]
    if t4_read.get("interpretable") is True:
        lines += [
            "Task 4 option (b) supersedes the primary external cell to the CELLxGENE 10x cohort "
            f"`{CXG_DATASET_ID}` and that cell passed (ρ > 0, null ≤ 0.05). "
            "Under the Task 4 reading, the frozen fibroblast ruler is a validated instrument "
            "outside GTEx and the Stage 2 trajectory is interpretable as a statement about the cells. "
            "The quoted T-A sentence is not edited in FINDINGS_FIBRO2.md.",
            "",
        ]
    elif t4_read.get("interpretable") is False:
        lines += [
            "Task 4 option (b) ran and did not pass. Stage 2 remains uninterpretable and every "
            "number in Tasks 1–3 is descriptive only. The quoted T-A sentence still applies.",
            "",
        ]
    else:
        lines += ["Task 4 not scored.", ""]

    lines += [
        "FINDINGS_FIBRO2.md T-B sentence, quoted:",
        "",
        f"> {FIBRO2_TB_SENTENCE}",
        "",
        "Task 1 locates that discrepancy in the joint-cluster tables above. "
        "It does not replace the T-B sentence and does not make d7 the Stage 2 endpoint.",
        "",
        "## Limitations",
        "",
        "1. Two donors in GSE297234.",
        "2. One timepoint spacing (days 0, 3, 7, 10); nothing is measured between d7 and d10.",
        "3. Clustering choice is frozen (PCA svd_solver=randomized n_components=20, KMeans k=3, seed `20260914`) and was not swept. IncrementalPCA was listed first and did not complete; randomized PCA replaced it before any age score.",
        "4. Clusters are not lineage-tracked, so \"the same cells\" is an inference from cluster identity, not from barcodes.",
        "5. Frozen ruler trained on GTEx V10 cultured fibroblasts, public `AGE` 10-year bins.",
        "6. Cross-platform: GTEx bulk polyA RNASeQCv2.4.2 vs 10x 3' scRNA-seq pseudobulk vs CELLxGENE 10x.",
        "7. Missing overlap genes at z=0. Nothing fitted on GSE297234 or on the CELLxGENE matrix.",
        "8. MIN_CELLS=20. Cluster×timepoint rows below that threshold are QC-only.",
        "9. Seed `20260914` (permutation / random directions / k-means). Donor-bootstrap seed `20260918`. n_perm=200, n_boot=200, n_random=200.",
        "10. GSE325735 not opened. Option (a) (GSE113957 FPKM rank-transform) was not run.",
        "11. Task 2 distances have no pass/fail.",
        "12. Uncalibrated R² is reported and is never a gate.",
        "",
        "## Files",
        "",
        "| path | content |",
        "|---|---|",
        "| `src/fibro3_common.py`, `fibro3_task4.py`, `fibro3_task1.py`, `fibro3_task2.py`, `fibro3_task3.py`, `fibro3_findings.py`, `fibro3_run.py` | code |",
        "| `results/fibro3/` | manifest, flags, tables, logs |",
        "| `FINDINGS_FIBRO3.md` | this file |",
        "| `PROGRESS_FIBRO3.md` | resume state |",
        "",
    ]

    fails = man.get("failures") or []
    if fails:
        lines += ["## Failures (manifest)", ""]
        for f in fails:
            lines.append(f"- **{f.get('step')}:** {f.get('message')}")
        lines.append("")

    cols = man.get("columns") or {}
    if cols:
        lines += ["## Columns actually read", ""]
        for k, v in cols.items():
            lines.append(f"- **{k}** source=`{v.get('source')}` columns={v.get('columns')}")
        lines.append("")

    FINDINGS_PATH.write_text("\n".join(lines) + "\n", encoding="utf-8")
    return FINDINGS_PATH
