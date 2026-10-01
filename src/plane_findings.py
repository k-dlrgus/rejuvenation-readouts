"""Write FINDINGS_PLANE.md from results/plane/. Does not modify other FINDINGS*.md or FALSIFICATION.md."""
from __future__ import annotations

import sys
from pathlib import Path

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parent))
from plane_common import (  # noqa: E402
    PLANE_DIR, PLANE_SEED, PLANE_BOOT, N_PERM, N_BOOT, N_RANDOM_DIR, MIN_CELLS,
    PREREG_TASK1, PREREG_TASK2, PREREG_TASK3,
    PREREG_TASK1_FLAG, PREREG_TASK2_FLAG, PREREG_TASK3_FLAG, DECLARED_BEFORE_SCORES_FLAG,
    FINDINGS_PATH, FROZEN_RULER, AMS_NBIN, AMS_CTRL, LOUVAIN_RES, PCA_MAHAL_K,
    AGED_LINE, YOUNG_LINE, TRAJECTORY, OFF_TRAJECTORY, AGE_INSTRUMENTS, INST_LABELS,
    IDENTITY_KEY, PLURI_ENDOGENOUS, FIBRO_IDENTITY, ROOT,
    MD3_MD_VS_PLURI_GM00731, MD3_MD_VS_PLURI_GM23815,
    FIBRO3_EXTRAP_RHO_NN, FIBRO3_EXTRAP_CI, FIBRO3_EXTRAP_N, MD4_EXTRAP_RHO_NN,
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


def _num(x, d=3):
    try:
        v = float(x)
    except (TypeError, ValueError):
        return "NA"
    if not np.isfinite(v):
        return "NA"
    return f"{v:.{d}f}"


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
        inner = []
        for inst in AGE_INSTRUMENTS:
            key = (rec.get(inst) or {}).get("key", "not scored")
            inner.append(f"{INST_LABELS.get(inst, inst)}=`{key}`")
        bits.append(f"{line} " + "; ".join(inner))
    t1_txt = "Task 1: " + ". ".join(bits) + "."
    t2_txt = "Task 2: report-only extrapolation tabulated; no gate."
    t3_txt = "Task 3: one figure per donor."
    return (
        f"**Status:** {t1_txt} {t2_txt} {t3_txt} "
        f"Seed `{PLANE_SEED}`. boot `{PLANE_BOOT}`. n_perm={N_PERM}. n_boot={N_BOOT}. "
        f"n_random={N_RANDOM_DIR}. Frozen ruler `{FROZEN_RULER.name}` exists={FROZEN_RULER.exists()}."
    )


def _md5_ledger_quote():
    path = ROOT / "FINDINGS_MD5.md"
    if not path.exists():
        return "_FINDINGS_MD5.md missing. Ledger not quoted._"
    text = path.read_text(encoding="utf-8")
    start = text.find("## Task 3 — one-line ledger")
    end = text.find("## What this changes")
    if start < 0 or end < 0 or end <= start:
        return "_FINDINGS_MD5.md ledger section not found. Not reconstructing from memory._"
    return text[start:end].strip()


def _reading_block(t1r):
    by = (t1r or {}).get("by_donor") or {}
    lines = []
    for line in (AGED_LINE, YOUNG_LINE):
        rec = by.get(line) or {}
        lines.append(f"### {line}")
        lines.append("")
        for inst in AGE_INSTRUMENTS:
            r = rec.get(inst) or {}
            lines.append(
                f"**{INST_LABELS.get(inst, inst)}** key=`{r.get('key', 'NA')}` "
                f"segment_B_end=`{r.get('segment_B_end')}` "
                f"fallback={r.get('segment_B_end_is_fallback')} "
                f"null=`{r.get('null_kind')}`."
            )
            lines.append("")
            lines.append(
                f"slope_A={_signed(r.get('slope_A'), 4)} slope_B={_signed(r.get('slope_B'), 4)} "
                f"bend={_signed(r.get('bend'), 4)} "
                f"p_more_negative={_num(r.get('p_more_negative'), 4)} "
                f"p_more_positive={_num(r.get('p_more_positive'), 4)}."
            )
            lines.append("")
            lines.append(r.get("text") or "_no text._")
            lines.append("")
            flags = r.get("flags") or {}
            if flags:
                lines.append(
                    "flags: " + "; ".join(f"{k}={v}" for k, v in flags.items())
                )
                lines.append("")
    return lines


def _instrument_comparison(t1r):
    by = (t1r or {}).get("by_donor") or {}
    cmp_ = (t1r or {}).get("comparison") or {}
    lines = [
        "Instruments are not combined. No winner is declared.",
        "",
    ]
    if cmp_.get("note"):
        lines.append(str(cmp_["note"]))
        lines.append("")
    else:
        lines.append(
            f"FINDINGS_MD3.md Task 3: MD vs pluripotency−fibroblast ρ="
            f"{MD3_MD_VS_PLURI_GM00731} (GM00731 n=33) and ρ={MD3_MD_VS_PLURI_GM23815} "
            "(GM23815 n=34). MD is expected to run in lockstep."
        )
        lines.append("")
    for inst in AGE_INSTRUMENTS:
        keys = []
        for line in (AGED_LINE, YOUNG_LINE):
            keys.append(f"{line}=`{((by.get(line) or {}).get(inst) or {}).get('key', 'NA')}`")
        lines.append(f"- {INST_LABELS.get(inst, inst)}: " + "; ".join(keys))
    lines.append("")
    md_bits = []
    for line in (AGED_LINE, YOUNG_LINE):
        r = (by.get(line) or {}).get("md_score") or {}
        flags = r.get("flags") or {}
        md_bits.append(
            f"{line} key=`{r.get('key', 'NA')}` "
            f"age_falls_A={flags.get('age_falls_A')} "
            f"age_falls_B={flags.get('age_falls_B')} "
            f"slopes_overlap={flags.get('slopes_overlap')} "
            f"bend_beats_null={flags.get('bend_beats_null')}"
        )
    lines.append("MD vs the lockstep definition: " + "; ".join(md_bits) + ".")
    md_keys = [
        ((by.get(line) or {}).get("md_score") or {}).get("key")
        for line in (AGED_LINE, YOUNG_LINE)
    ]
    ruler_keys = [
        ((by.get(line) or {}).get("frozen_ruler") or {}).get("key")
        for line in (AGED_LINE, YOUNG_LINE)
    ]
    if all(k == "lockstep" for k in md_keys if k):
        lines.append(
            "MD ran in lockstep on every scored donor, matching the FINDINGS_MD3.md Task 3 "
            "expectation from ρ ≈ −0.94 with the pluripotency−fibroblast score."
        )
    else:
        lines.append(
            "MD did not fire `lockstep` on either donor. FINDINGS_MD3.md Task 3 ρ ≈ −0.94 "
            "predicted lockstep; that is not the key that fired. No winner is declared."
        )
    if any(k == "bend" for k in ruler_keys) and any(k == "lockstep" for k in md_keys):
        lines.append(
            "Observed: the frozen ruler is `bend` on at least one donor and MD is `lockstep` "
            "on at least one donor. That is the separation between the instruments on this "
            "plane, stated as observed. No winner is declared."
        )
    lines.append("")
    return lines


def write_findings():
    man = load_manifest() if (PLANE_DIR / "manifest.json").exists() else {}
    gs = load_json(PLANE_DIR / "genesets_summary.json") if (PLANE_DIR / "genesets_summary.json").exists() else {}
    t1s = load_json(PLANE_DIR / "t1_summary.json") if (PLANE_DIR / "t1_summary.json").exists() else {}
    t1r = load_json(PLANE_DIR / "t1_reading.json") if (PLANE_DIR / "t1_reading.json").exists() else {}
    t1san = load_json(PLANE_DIR / "t1_sanity.json") if (PLANE_DIR / "t1_sanity.json").exists() else {}
    t2s = load_json(PLANE_DIR / "t2_summary.json") if (PLANE_DIR / "t2_summary.json").exists() else {}
    t3f = load_json(PLANE_DIR / "t3_figures.json") if (PLANE_DIR / "t3_figures.json").exists() else {}

    lines = [
        "# FINDINGS_PLANE — does the age–identity curve bend?",
        "",
        _status_line(t1r, t2s),
        "",
        "Does not modify any existing `FINDINGS_*.md` or `FALSIFICATION.md`. No GSE325735. "
        "Uncalibrated R² is reported and is never a gate. Nothing averaged across donors "
        "or instruments. Refit nothing. d7 is not the primary Stage 2 endpoint. "
        "Task 2 is report-only and is not used to discount Task 1.",
        "",
        f"Reproduced by `src/plane_run.py`. Frozen ruler: `{FROZEN_RULER}` exists={FROZEN_RULER.exists()}.",
        "",
        f"Flag: `PREREG_TASK1.flag` exists={PREREG_TASK1_FLAG.exists()}. "
        f"`PREREG_TASK2.flag` exists={PREREG_TASK2_FLAG.exists()}. "
        f"`PREREG_TASK3.flag` exists={PREREG_TASK3_FLAG.exists()}. "
        f"`DECLARED_BEFORE_SCORES.flag` exists={DECLARED_BEFORE_SCORES_FLAG.exists()}.",
        "",
        "## The question",
        "",
        "Place reprogramming cells on the age–identity plane and ask whether the curve "
        "bends — whether age can move before identity is lost — or whether age and identity "
        "move in lockstep. FINDINGS_MD5.md established that partially reprogrammed cells "
        "read younger than their starting fibroblasts on the frozen ruler. The live "
        "alternative is that the ruler reads younger because cells are leaving fibroblast "
        "identity, not because they are rejuvenating.",
        "",
        "- Lockstep (identity loss): the age reading keeps falling as cells move further "
        "from fibroblast, monotonically, all the way to pluripotency.",
        "- Bend (rejuvenation separable from dedifferentiation): the age reading falls from "
        "Fibroblast to PartialReprog, then stops falling or reverses as cells continue "
        "toward EarlyPluripotency and Pluripotency, while identity keeps moving.",
        "",
        "A bend is the result the project was built to find. Lockstep is equally reportable "
        "and is not argued away. Trajectory order, fixed before scoring: "
        f"{' → '.join(TRAJECTORY)}. `{OFF_TRAJECTORY}` is off-trajectory: scored and plotted, "
        "excluded from every trajectory statistic. Unit: cell state per donor; timepoints "
        f"pooled within state. ≥{MIN_CELLS} cells per state per donor or the state is not scored. "
        f"x = identity = `{IDENTITY_KEY}` (FINDINGS_FIBRO.md frozen lists; md4 primary, "
        "drop_oskm=False). y = age. Three instruments, never combined: frozen ruler; MD; "
        "Table S3 age-up minus age-down.",
        "",
        "## Failures recorded",
        "",
        _failures_block(man),
        "## Task 1 pre-registration (verbatim, written before any PLANE state-median, slope, bend, or null p)",
        "",
        PREREG_TASK1,
        "",
        "## Task 1 — state cell counts and timepoint composition (md2 Louvain labels, not re-labelled)",
        "",
        "Assignment basis: argmax of mean AddModuleScore of mmc3 `Reprog_cell_state_signatures` "
        "(Fibroblast, PartialReprog, EarlyPluripotency, Pluripotency, NonReprog). "
        "Labels from `results/md2/t2_cluster_labels_*.csv`. Not re-labelled. "
        "Occupancy matches `results/md4/t1_cell_counts.csv` or this file stops.",
        "",
        "Pooled n cells by state (timepoints pooled within state for scoring; composition reported):",
        "",
        _md(PLANE_DIR / "t1_cell_counts.csv"),
        "",
        "Timepoint composition of every state:",
        "",
        _md(PLANE_DIR / "t1_timepoint_composition.csv", [
            "cell_line", "day", "label", "n_cells", "n_timepoint", "n_state",
            "frac_state", "frac_timepoint", "state_ge30",
        ]),
        "",
        "Qualify (segment A = Fibroblast → PartialReprog; segment B = PartialReprog → Pluripotency, "
        "else EarlyPluripotency and said so):",
        "",
        _md(PLANE_DIR / "t1_qualify.csv"),
        "",
        "## Task 1 — sanity check (FINDINGS_MD2.md pooled MD means; occupancy vs md4)",
        "",
    ]
    if t1san:
        lines += [
            f"- n_PartialReprog={t1san.get('n_PartialReprog')} "
            f"n_NonReprog={t1san.get('n_NonReprog')} "
            f"MD_PR={_signed(t1san.get('md_mean_PartialReprog'), 6)} "
            f"MD_NR={_signed(t1san.get('md_mean_NonReprog'), 6)} "
            f"ok={t1san.get('ok')}",
            f"- md2 json n_PR={t1san.get('md2_json_n_PartialReprog')} "
            f"n_NR={t1san.get('md2_json_n_NonReprog')} "
            f"MD_PR={_signed(t1san.get('md2_json_md_mean_PartialReprog'), 6)} "
            f"MD_NR={_signed(t1san.get('md2_json_md_mean_NonReprog'), 6)}",
            f"- printed 3 d.p. match (+0.150 / +0.439): {t1san.get('printed_3dp_match')}",
            "",
        ]
    else:
        lines += ["_t1_sanity.json missing._", ""]

    lines += [
        "## Task 1 — state medians with bootstrap 95% CIs (B=200, cells within state)",
        "",
        "Identity = pluripotency−fibroblast (`pluri_primary`). Age instruments in separate rows. "
        "NonReprog is scored and plotted; it is not used in slopes or the bend statistic.",
        "",
        _md(PLANE_DIR / "t1_state_medians.csv", [
            "cell_line", "instrument", "label", "n_cells", "scored", "on_trajectory",
            "identity_median", "identity_ci_lo", "identity_ci_hi",
            "age_median", "age_ci_lo", "age_ci_hi",
        ]),
        "",
        "## Task 1 — segment slopes, bend statistic, nulls",
        "",
        "slope = Δage / Δidentity. Bend = slope_A − slope_B. Frozen-ruler null = 200 "
        "permuted-weight random directions, seed 20260914. Gene-list nulls = 200 size-matched, "
        "expression-bin-matched random gene sets. Identity held fixed. "
        "p_more_negative = (n_null ≤ obs + 1) / (n_null + 1). Beats null iff p_more_negative ≤ 0.05.",
        "",
        _md(PLANE_DIR / "t1_slopes.csv", [
            "cell_line", "instrument", "state_B1", "segment_B_end_is_fallback", "ok",
            "delta_identity_A", "delta_age_A", "slope_A", "slope_A_ci_lo", "slope_A_ci_hi",
            "delta_identity_B", "delta_age_B", "slope_B", "slope_B_ci_lo", "slope_B_ci_hi",
            "bend", "bend_ci_lo", "bend_ci_hi",
            "slope_A_young", "slope_B_young", "slope_B_flat", "slope_B_old", "slopes_overlap",
        ]),
        "",
        _md(PLANE_DIR / "t1_null_summary.csv", [
            "cell_line", "instrument", "ok", "null_kind", "n_random", "n_null_finite",
            "obs_bend", "null_bend_median", "null_bend_ci_lo", "null_bend_ci_hi",
            "p_more_negative", "p_more_positive", "segment_B_end", "segment_B_end_is_fallback",
        ]),
        "",
        _md(PLANE_DIR / "t1_bend.csv", [
            "cell_line", "instrument", "key", "ok", "slope_A", "slope_B", "bend",
            "p_more_negative", "p_more_positive", "bend_beats_null",
            "slope_A_young", "slope_B_young", "slope_B_flat", "slope_B_old", "slopes_overlap",
        ]),
        "",
        "## Task 1 — reading (only what fired)",
        "",
    ]
    lines += _reading_block(t1r)
    if t1s:
        lines += [
            f"t1_summary by_donor={t1s.get('by_donor')} sanity_ok={t1s.get('sanity_ok')}",
            "",
        ]
    lines += [
        "## Instrument comparison (report-only)",
        "",
    ]
    lines += _instrument_comparison(t1r)
    lines += [
        "## Task 2 pre-registration (verbatim, written before any PLANE extrapolation distance)",
        "",
        PREREG_TASK2,
        "",
        "## Task 2 — extrapolation to GTEx fibroblast training (report-only, no gate)",
        "",
        f"PCA k={PCA_MAHAL_K} nearest-neighbour Euclidean and Mahalanobis, same measures as "
        "FINDINGS_FIBRO2.md / `src/fibro2_extrap.py`. TMM among that donor's scored-state rows "
        "(not the md4 state×timepoint rows). Frozen μ/σ/w. Not a gate.",
        "",
        _md(PLANE_DIR / "t2_extrap.csv", [
            "cell_line", "label", "n_cells", "scored", "on_trajectory",
            "age_score", "nn_euclidean", "mahalanobis_pca", "pca_k", "reason",
        ]),
        "",
        f"Inherited ρ (not recomputed as a gate): FINDINGS_FIBRO3.md Task 2 nn_euclidean "
        f"ρ={FIBRO3_EXTRAP_RHO_NN} CI=[{FIBRO3_EXTRAP_CI[0]}, {FIBRO3_EXTRAP_CI[1]}] "
        f"n={FIBRO3_EXTRAP_N}. FINDINGS_MD4.md ρ_nn={MD4_EXTRAP_RHO_NN}.",
        "",
        "Direction statement (from the distances in `t2_extrap.csv` and those inherited ρ values):",
        "",
        t2s.get("direction_statement") if t2s else "_t2_summary.json missing._",
        "",
        f"report_only={t2s.get('report_only') if t2s else 'NA'} "
        f"used_to_discount_task1={t2s.get('used_to_discount_task1') if t2s else 'NA'}",
        "",
        "## Task 3 — figure",
        "",
        PREREG_TASK3,
        "",
    ]
    figs = (t3f or {}).get("figures") or {}
    if figs:
        for line in (AGED_LINE, YOUNG_LINE):
            p = figs.get(line)
            lines.append(f"- {line}: `{p}`" if p else f"- {line}: missing")
        lines.append("")
        lines.append(
            "Three panels per donor (frozen ruler, MD, age-up−age-down). x = identity, y = age. "
            "Trajectory states as points with 2D bootstrap CIs joined in order. NonReprog as a "
            "separate marker. Random-direction / random-set envelope shaded. Same identity-axis "
            "scaling within a donor. Y-axis per instrument (instruments are not combined)."
        )
        lines.append("")
    else:
        lines += ["_t3_figures.json missing._", ""]

    lines += [
        "## Updated ledger (quoted from FINDINGS_MD5.md; that file is not edited)",
        "",
        "Append this file's reading. Quote the FINDINGS_MD5.md ledger; edit nothing in it.",
        "",
        _md5_ledger_quote(),
        "",
        "### FINDINGS_PLANE.md Task 1 (this file)",
        "",
        "File: `FINDINGS_PLANE.md`. Not a replacement of any prior ledger entry. Not an edit "
        "of FINDINGS_FIBRO.md Stage 2 or of any FINDINGS_MD5.md ledger line.",
        "",
    ]
    by = (t1r or {}).get("by_donor") or {}
    for line in (AGED_LINE, YOUNG_LINE):
        rec = by.get(line) or {}
        for inst in AGE_INSTRUMENTS:
            r = rec.get(inst) or {}
            lines.append(
                f"- {line} {INST_LABELS.get(inst, inst)}: `{r.get('key', 'NA')}`. "
                f"{r.get('text') or ''}"
            )
            lines.append("")
    lines += [
        "## Limitations",
        "",
        "1. Two donors in GSE297234.",
        "2. States are inferred from expression, not lineage. PartialReprog / Fibroblast / "
        "EarlyPluripotency / Pluripotency / NonReprog labels are argmax of mmc3 "
        "Reprog_cell_state_signatures AddModuleScore, not lineage-tracked cells from day-0 fibroblasts.",
        "3. A state ordering is not a measured trajectory. The order Fibroblast → PartialReprog → "
        "EarlyPluripotency → Pluripotency was fixed before scoring; it is not a Slingshot or "
        "time-resolved path of the same cell.",
        "4. Timepoints are pooled within state. The state, not the day, is the trajectory position. "
        "Timepoint composition of every state is reported and is not a substitute for a day-resolved path.",
        "5. The ruler was trained on GTEx V10 cultured fibroblasts only, public `AGE` 10-year bins. "
        "Per-cell scores are the md4 log2-CPM (TMM nf=1) vectors; not refit. Missing overlap genes at z=0.",
        "6. Louvain is a Python deviation from their R/sctransform pipeline, as recorded in "
        "FINDINGS_MD2.md (LogNormalize + quadratic HVG + percent.mt residualization + PCA 50 + "
        f"SNN + networkx louvain, resolution={LOUVAIN_RES} Seurat default; STAR Methods omit resolution).",
        "7. Our reproduction of their metrics is ours, not theirs. AddModuleScore is a Python "
        f"reimplementation of the published algorithm (nbin={AMS_NBIN} ctrl={AMS_CTRL}), not Seurat's C++/R object.",
        "8. Cross-platform shift from GTEx bulk polyA (RNASeQCv2.4.2) to 10x 3' scRNA-seq.",
        "9. MIN_CELLS=30. The threshold is not lowered. Donors are not pooled. Instruments are not combined.",
        "10. Seed `20260914`. Bootstrap seed `20260918`. n_perm=200, n_boot=200, n_random=200.",
        "11. GSE325735 not opened. d7 is not the Stage 2 endpoint. FINDINGS_FIBRO.md Stage 2 d0→d10 is not changed.",
        "12. Task 2 distances are on TMM state-pooled pseudobulks; Task 1 ruler scores are per-cell log2-CPM. "
        "Those two ruler numbers are not interchangeable. Task 2 is not used to explain away a lockstep result.",
        "13. Frozen pluripotency list "
        f"{list(PLURI_ENDOGENOUS)}; fibroblast identity {list(FIBRO_IDENTITY)}; "
        "identity version md4 pluri_primary drop_oskm=False.",
        "",
        "## Gene-list versions (from disk, never from memory)",
        "",
    ]
    if gs:
        lines += [
            f"- MD n={gs.get('n_MD')} Age up n={gs.get('n_age_up')} Age down n={gs.get('n_age_down')}",
            f"- used_as_deseq2_rebuild={gs.get('used_as_deseq2_rebuild')} "
            f"used_as_published_lists={gs.get('used_as_published_lists')}",
            f"- identity_version={gs.get('identity_version')}",
            f"- pluripotency endogenous: {gs.get('pluri_endogenous')}",
            f"- fibroblast identity: {gs.get('fibro_identity')}",
            "",
        ]
        s3 = gs.get("s3_mapped") or {}
        for line in (AGED_LINE, YOUNG_LINE):
            rec = s3.get(line) or {}
            lines.append(f"- {line} S3 AMS meta: {jsonable(rec)}")
            lines.append("")
    else:
        lines += ["_genesets_summary.json missing._", ""]

    idt = man.get("id_types") or {}
    lines += ["## ID types actually read", ""]
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
    lines.append("")
    cols = man.get("columns") or {}
    lines += ["## Columns actually read", ""]
    if cols:
        for k, v in cols.items():
            lines.append(f"- **{k}** source=`{v.get('source')}` columns={v.get('columns')}")
    else:
        lines.append("- none recorded")
    lines.append("")
    gsets = man.get("gene_sets_read") or []
    lines += ["## Gene-set names actually read", ""]
    if gsets:
        for g in gsets:
            lines.append(
                f"- **{g.get('key')}** name=`{g.get('name')}` n={g.get('n')} "
                f"id_type={g.get('id_type')} source=`{g.get('source')}`"
            )
    else:
        lines.append("- none recorded")
    lines.append("")
    lines += [
        "## Files",
        "",
        "- `src/plane_common.py`, `src/plane_task1.py`, `src/plane_nulls.py`, "
        "`src/plane_task2.py`, `src/plane_task3.py`, `src/plane_findings.py`, `src/plane_run.py`",
        "- `results/plane/`",
        "- `results/plane/figures/plane_GM00731.png`, `results/plane/figures/plane_GM23815.png`",
        "- `FINDINGS_PLANE.md`",
        "- `PROGRESS_PLANE.md`",
        "",
        "Existing `FINDINGS_*.md` and `FALSIFICATION.md` were not modified.",
        "",
    ]
    FINDINGS_PATH.write_text("\n".join(str(x) for x in lines) + "\n", encoding="utf-8")
    progress_snapshot("done", stop="FINDINGS_WRITTEN")
    return FINDINGS_PATH


if __name__ == "__main__":
    write_findings()
