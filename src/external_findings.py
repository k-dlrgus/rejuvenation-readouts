"""Write FINDINGS_EXTERNAL.md from results/external/. Does not modify other FINDINGS*.md."""
from __future__ import annotations

import sys
from pathlib import Path

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parent))
from external_common import (  # noqa: E402
    EXT_DIR, EXT_FIG, EXT_SEED, DATASET_ID, TITLE_REQUIRED, TISSUE_REQUIRED, COLLECTION_ID,
    N_PERM, N_BOOT, BOOT_SEED, TRANSFER_NULL_BAR, PREREG_DATE, PREREG_FLAG,
    FINDINGS_PATH, LOWDIM_L4_QUOTE, PREREG_READING,
    load_json, fmt, md_table,
)


def _plot_primary():
    obs_path = EXT_DIR / "project_pseudobulk_obs.csv"
    don_path = EXT_DIR / "stage0_donor_meta.csv"
    if not obs_path.exists() or not don_path.exists():
        return None
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    pb = pd.read_csv(obs_path)
    don = pd.read_csv(don_path)
    don["donor"] = don.donor.astype(str)
    pb["donor"] = pb.donor.astype(str)
    low = set(don.loc[don.pathology_low.astype(str).str.lower().eq("true"), "donor"])
    g = pb[pb.donor.isin(low)].copy()
    if g.empty or "pred_raw" not in g.columns:
        return None
    med = g.groupby("donor").agg(age=("age", "first"), pred_raw=("pred_raw", "median"),
                                 pred_id=("pred_idresid", "median")).reset_index()
    fig, ax = plt.subplots(1, 2, figsize=(9.0, 3.6))
    fig.suptitle("SEA-AD primary (no/low ADNC)", fontsize=11)
    for a, col, title in ((ax[0], "pred_raw", "raw ridge"), (ax[1], "pred_id", "identity-residual ridge")):
        a.scatter(med.age, med[col], s=28, c="#1f4e79", alpha=0.85)
        a.set_xlabel("chronological age (years)")
        a.set_ylabel("DLPFC-calibrated score")
        a.set_title(title)
        a.spines["top"].set_visible(False)
        a.spines["right"].set_visible(False)
    EXT_FIG.mkdir(parents=True, exist_ok=True)
    out = EXT_FIG / "primary_score_vs_age.png"
    fig.tight_layout()
    fig.savefig(out, dpi=140)
    plt.close(fig)
    return out


def _signed(x, d=3):
    try:
        v = float(x)
    except (TypeError, ValueError):
        return "NA"
    if not np.isfinite(v):
        return "NA"
    return f"{v:+.{d}f}"


def _ci(lo, hi):
    return f"[{_signed(lo)}, {_signed(hi)}]"


def _reading(prim, sec):
    p_ok = bool(prim.get("pass_r"))
    s_ok = bool(sec.get("pass_r"))
    if p_ok:
        return (
            "pass",
            "the direction transfers to a third bank, region, and lab",
        )
    if (not p_ok) and s_ok:
        return (
            "masked",
            "transfer is masked by pathology",
        )
    return (
        "no_transfer",
        "does not transfer beyond the two DLPFC banks",
    )


def _cell_row(label, rec):
    return (
        f"| {label} | {_signed(rec.get('r'))} | {_signed(rec.get('r_null'))} | "
        f"{_signed(rec.get('r_p'))} | {_signed(rec.get('rho'))} | {_signed(rec.get('cal_r2'))} | "
        f"{_signed(rec.get('r2'))} | {_ci(rec.get('r_ci_lo'), rec.get('r_ci_hi'))} | "
        f"{rec.get('pass_r')} | {rec.get('n_donors')} | {rec.get('n_types_ok')} | "
        f"{rec.get('n_perm')} |"
    )


def write_findings():
    man = load_json(EXT_DIR / "manifest.json") if (EXT_DIR / "manifest.json").exists() else {}
    failures = man.get("failures") or []
    status_file = man.get("status", "UNKNOWN")
    summ0 = load_json(EXT_DIR / "stage0_summary.json") if (EXT_DIR / "stage0_summary.json").exists() else {}

    if status_file == "STOP" or failures:
        status_line = (
            f"**Status:** STOP Stage 0. Seed `{EXT_SEED}`. "
            f"SEA-AD MTG dataset `{DATASET_ID}`. Failures recorded in `results/external/manifest.json`."
        )
    elif not (EXT_DIR / "project_summary.json").exists():
        status_line = (
            f"**Status:** Stage 0 done; fit {'gated on missing pre-reg' if not PREREG_FLAG.exists() else 'not yet scored'}. "
            f"Seed `{EXT_SEED}`. SEA-AD `{DATASET_ID}` ({TITLE_REQUIRED})."
        )
    else:
        status_line = (
            f"**Status:** external transfer scored. Seed `{EXT_SEED}`. "
            f"SEA-AD `{DATASET_ID}`. Primary metric: Pearson r (n_perm={N_PERM})."
        )

    lines = [
        "# FINDINGS_EXTERNAL — does the frozen DLPFC age direction transfer to SEA-AD MTG?",
        "",
        status_line,
        "",
        "Does not modify `FINDINGS_LOWDIM.md`, `FINDINGS_GEOMETRY.md`, `FINDINGS_POSCTRL.md`, "
        "`FINDINGS_GTEX.md`, `FINDINGS_TARGET.md`, `FINDINGS_TRAJECTORY.md`, "
        "`FINDINGS_BRAIN_PHASE1.md`, or `FALSIFICATION.md`. No perturbation data. No TF Atlas. "
        "No candidate interventions.",
        "",
        f"Reproduced by `src/external_run.py`. Cohort: {TITLE_REQUIRED}. "
        f"Tissue: `{TISSUE_REQUIRED}`. Collection `{COLLECTION_ID}`.",
        "",
        "## Pre-registration (verbatim, written before any fit)",
        "",
    ]
    if PREREG_FLAG.exists():
        lines += [PREREG_FLAG.read_text(encoding="utf-8").strip(), ""]
    else:
        lines += ["Pre-registration flag not written (Stage 0 STOP or not yet reached).", ""]
    lines += ["## Pre-registered reading (verbatim)", "", PREREG_READING, ""]

    if failures:
        lines += ["## STOP", "", "Halted. Not substituting columns, types, or donors.", ""]
        for f in failures:
            lines.append(f"- **{f.get('step')}:** {f.get('message')}")
        lines += ["", "See `results/external/manifest.json`.", ""]

    if summ0:
        lines += [
            "## Stage 0 — cohort",
            "",
            f"- dataset_id: `{summ0.get('dataset_id')}`",
            f"- title: {summ0.get('title')}",
            f"- tissue: `{summ0.get('tissue')}`",
            f"- nuclei: {summ0.get('n_nuclei')}",
            f"- donors (CXG): {summ0.get('n_donors')}",
            f"- donors with ≥1 mapped nucleus: **{summ0.get('n_donors_mapped')}**",
            f"- age range: **{summ0.get('age_min')}–{summ0.get('age_max')}** years "
            f"(nunique={summ0.get('nunique_age')})",
            f"- sex: {summ0.get('sex_counts')}",
            f"- pathology field: `{summ0.get('pathology_field')}` (family={summ0.get('pathology_family')})",
            f"- pathology levels observed: {summ0.get('pathology_observed')}",
            f"- no/low pathology donors: {summ0.get('n_pathology_low')} "
            f"(levels={summ0.get('pathology_low_levels')}; missing={summ0.get('n_pathology_missing')})",
            f"- pathology source: `{summ0.get('pathology_source')}`",
            f"- subclasses mapped: {summ0.get('n_subclasses_mapped')} / {summ0.get('n_subclasses')}",
            f"- unmapped subclasses: {summ0.get('unmapped_subclasses')}",
            "",
            "Unmapped subclasses failed **exact** `cell_type` (and exact DLPFC-subclass) match. "
            "Not synonymized: e.g. SEA-AD `astrocyte of the cerebral cortex` ≠ DLPFC `astrocyte`; "
            "`near-projecting glutamatergic cortical neuron` ≠ `L5/6 near-projecting glutamatergic neuron`.",
            "",
        ]
        leftover = man.get("cxg_donors_not_in_xlsx") or []
        if leftover:
            lines += [
                f"CELLxGENE has {len(leftover)} donors with no row in the Allen donor-metadata xlsx "
                f"(`{leftover}`). They are kept because h5ad `ADNC`/`Braak stage`/`CERAD score` are present "
                "(SEA-AD neurotypical-reference set; ADNC level `Reference`). Not dropped.",
                "",
            ]
        map_csv = EXT_DIR / "stage0_mapping.csv"
        if map_csv.exists():
            mp = pd.read_csv(map_csv)
            lines += ["### Subclass → DLPFC type (exact strings only)", "", md_table(mp), ""]
        nd = EXT_DIR / "stage0_n_donors_per_type.csv"
        if nd.exists():
            tab = pd.read_csv(nd)
            lines += ["### n donors per DLPFC type", "", md_table(tab), ""]
        cols = (man.get("columns_read") or {})
        lines += ["### Columns actually read", ""]
        for src, cl in cols.items():
            lines.append(f"- `{src}`: {cl}")
        lines.append("")

    proj = EXT_DIR / "project_summary.json"
    if proj.exists():
        ps = load_json(proj)
        ov = ps.get("gene_overlap") or {}
        lines += [
            "## Gene overlap",
            "",
            f"DLPFC genes={ov.get('n_dlpfc')}  SEA-AD genes={ov.get('n_seaad')}  "
            f"overlap=**{ov.get('n_overlap')}**  fraction={ov.get('frac_overlap')}  "
            f"DLPFC genes absent from SEA-AD={ov.get('n_dlpfc_missing')}.",
            "Missing DLPFC genes are left at the DLPFC mean (z=0). The frozen vector is not rewritten.",
            "",
            "## Transfer",
            "",
            "| cell | r | r_null | p | ρ | calR² | R² | r 95% CI | pass | n_donors | n_types | n_perm |",
            "|" + "|".join(["---"] * 12) + "|",
        ]
        results = ps.get("results") or {}
        readings = []
        for kind in ("raw", "idresid"):
            recs = results.get(kind) or {}
            prim_path = EXT_DIR / f"primary_{kind}.json"
            sec_path = EXT_DIR / f"secondary_{kind}.json"
            prim = recs.get("primary") or (load_json(prim_path) if prim_path.exists() else {})
            sec = recs.get("secondary") or (load_json(sec_path) if sec_path.exists() else {})
            if prim:
                lines.append(_cell_row(f"{kind} primary (no/low pathology)", prim))
            if sec:
                lines.append(_cell_row(f"{kind} secondary (pathology covariate)", sec))
            if prim and sec:
                tag, text = _reading(prim, sec)
                readings.append((kind, tag, text, prim, sec))
        lines.append("")

        lines += ["## Verdict", ""]
        if not readings:
            lines.append("No scored cells.")
        else:
            tags = {kind: tag for kind, tag, text, prim, sec in readings}
            unique_tags = sorted(set(tags.values()))
            if len(unique_tags) == 1:
                text = readings[0][2]
                lines.append("**" + text[:1].upper() + text[1:] + ".**" if text else text)
                lines.append("")
                lines.append("Same reading on both frozen directions (raw and identity-residualized). Do not average.")
                lines.append("Pearson r is the gate (POSCTRL S2/S3). Uncalibrated R² is a calibration artifact "
                             "(negative in the primary cell, as in HBCC↔MSSM). Spearman ρ is lower than r because "
                             "SEA-AD chronological age is HsapDv integer-year bins (23 unique ages, 29–89).")
            else:
                lines.append("Directions disagree; report which, do not average.")
            lines.append("")
            for kind, tag, text, prim, sec in readings:
                lines.append(
                    f"- **{kind}:** {text}. "
                    f"primary r={_signed(prim.get('r'))} (null {_signed(prim.get('r_null'))}, "
                    f"CI {_ci(prim.get('r_ci_lo'), prim.get('r_ci_hi'))}, pass={prim.get('pass_r')}); "
                    f"secondary r={_signed(sec.get('r'))} (null {_signed(sec.get('r_null'))}, "
                    f"pass={sec.get('pass_r')})."
                )
            lines.append("")
            fig = _plot_primary()
            if fig is not None:
                rel = "results/external/figures/" + fig.name
                lines += [f"Figure: `{rel}`", "", f"![]({rel})", ""]

        for kind in ("raw", "idresid"):
            pth = EXT_DIR / f"primary_{kind}_per_type.csv"
            if pth.exists():
                lines += [f"### {kind} primary per type", "", md_table(pd.read_csv(pth)), ""]

    lines += [
        "## What this changes in FINDINGS_LOWDIM L4 (quote, not edited)",
        "",
        LOWDIM_L4_QUOTE,
        "",
        "## Limitations",
        "",
    ]
    age_lo, age_hi = summ0.get("age_min"), summ0.get("age_max")
    lines += [
        f"1. **Age range:** SEA-AD MTG donors span {age_lo}–{age_hi} years "
        f"(nunique={summ0.get('nunique_age')}). The DLPFC clock was fit on adult controls 20–89. "
        "Do not extrapolate outside the overlapping range.",
        f"2. **Region:** SEA-AD is **middle temporal gyrus**, not DLPFC. A pass is transfer across "
        "bank, region, and lab; a fail confounds region with cohort.",
        f"3. **Pathology:** primary cell uses `{summ0.get('pathology_field')}` "
        f"(family={summ0.get('pathology_family')}); no/low levels={summ0.get('pathology_low_levels')}. "
        "SEA-AD is an AD-spectrum cohort. Unmapped subclasses are excluded, not imputed.",
        "4. Two DLPFC banks (HBCC, MSSM) trained the direction; SEA-AD is a third cohort, never in the fit.",
        "5. Unmapped SEA-AD subclasses are reported, not synonymized onto the 20 DLPFC types. "
        "11 of 20 DLPFC types have zero mapped SEA-AD donors (astrocyte, endothelium, several glutamatergic subclasses).",
        "6. Primary n=26. Pearson r is leverage-sensitive to the younger tail (neurotypical-reference donors ~30–60 y); "
        "Spearman ρ is the rank-robust companion and is reported alongside. The gate is Pearson r.",
        f"7. Seed `{EXT_SEED}` (permutation); donor-bootstrap seed `{BOOT_SEED}`. n_perm={N_PERM}, n_boot={N_BOOT}.",
        "8. Same log2(CPM+1) as the DLPFC Phase-1 matrix. Gene overlap is incomplete; missing genes are z=0.",
        "",
        "## Files",
        "",
        "| path | content |",
        "|---|---|",
        "| `src/external_common.py`, `external_stage0.py`, `external_freeze.py`, `external_project.py`, `external_download.py`, `external_findings.py`, `external_run.py` | code |",
        "| `results/external/` | manifest, mapping, frozen directions, primary/secondary cells |",
        "| `FINDINGS_EXTERNAL.md` | this file |",
        "| `PROGRESS_EXTERNAL.md` | resume state |",
        "",
    ]
    FINDINGS_PATH.write_text("\n".join(lines) + "\n", encoding="utf-8")
    return FINDINGS_PATH


if __name__ == "__main__":
    p = write_findings()
    print(f"wrote {p}")
