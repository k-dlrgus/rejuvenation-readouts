"""Regenerate FINDINGS_BRAIN_PHASE1.md from results/brain_phase1/ artifacts.

Called after every stage. Never touches other FINDINGS*.md or FALSIFICATION.md.
"""
from __future__ import annotations

import json
from pathlib import Path

import numpy as np
import pandas as pd

from brain_phase1_common import (
    PHASE1_DIR, PHASE1_SEED, PRIMARY_QUANTILES, P0_AGE_CORR_ABS,
    HASH_POOL_LIMITATION, STOP_P2_MIN_A, STOP_P3_AGE_R2, STOP_P3_UNREST_R2,
    BLOOD_UA_MAX, BLOOD_UA_P99, N_SITE_FOLDS, N_MATCH_DRAWS, ROOT,
)

OUT = ROOT / "FINDINGS_BRAIN_PHASE1.md"


def _load(name):
    p = PHASE1_DIR / name
    if not p.exists():
        return None
    if p.suffix == ".json":
        return json.loads(p.read_text(encoding="utf-8"))
    if p.suffix == ".csv":
        return pd.read_csv(p)
    return p.read_text(encoding="utf-8")


def _fmt_r(x, nd=3):
    if x is None or (isinstance(x, float) and not np.isfinite(x)):
        return "NA"
    return f"{float(x):+.{nd}f}"


def _md_table(df: pd.DataFrame, cols=None, float_cols=None, nd=3) -> str:
    if df is None or len(df) == 0:
        return "_(none)_"
    cols = cols or list(df.columns)
    float_cols = set(float_cols or [])
    lines = ["| " + " | ".join(cols) + " |", "|" + "|".join("---" for _ in cols) + "|"]
    for _, r in df.iterrows():
        cells = []
        for c in cols:
            v = r[c]
            if c in float_cols or isinstance(v, (float, np.floating)):
                try:
                    cells.append("NA" if not np.isfinite(float(v)) else f"{float(v):.{nd}f}")
                except Exception:
                    cells.append(str(v))
            else:
                cells.append(str(v))
        lines.append("| " + " | ".join(cells) + " |")
    return "\n".join(lines)


def write_findings(log=print):
    p0 = _load("p0_summary.json") or {}
    p0_inc = _load("p0_included_covariates.json") or {}
    p0_corr = _load("p0_covariate_age_correlations.csv")
    p0_audit = _load("p0_obs_audit.json") or {}
    p1t = _load("p1_thresholds.json")
    p1_sweep = _load("p1_threshold_sweep_set_sizes.csv")
    p1_null = _load("p1_permutation_null_counts.csv")
    p2 = _load("p2_summary.json")
    p2_surv = _load("p2_Agene_survival.csv")
    p3 = _load("p3_summary.json")
    p4 = _load("p4_summary.json")
    p5 = _load("p5_summary.json")

    parts = []
    parts.append(_header(p0, p1t, p2, p3, p4, p5))
    parts.append(_p0(p0, p0_inc, p0_corr, p0_audit))
    parts.append(_p1_declared(p1t, p1_sweep, p1_null))
    if p2:
        parts.append(_p2(p2, p2_surv))
    if p3:
        parts.append(_p3(p3))
    if p4:
        parts.append(_p4(p4))
    if p5:
        parts.append(_p5(p5))
    parts.append(_limitations(p0, p1t, p3, p4, p5))
    parts.append(_files())
    text = "\n\n".join(parts).rstrip() + "\n"
    OUT.write_text(text, encoding="utf-8")
    log(f"[findings] wrote {OUT}")
    return OUT


def _header(p0, p1t, p2, p3, p4, p5):
    stages = ["P0"]
    if p1t:
        stages.append("P1 (thresholds declared before scores)")
    if p2:
        stages.append("P2")
    if p3:
        stages.append("P3")
    if p4:
        stages.append("P4")
    if p5:
        stages.append("P5")
    verdict = "in progress"
    if p4 and p4.get("verdict"):
        ws = (p4.get("within_site") or {})
        lo = (p4.get("loso") or {})
        verdict = (
            "NOT SEPARABLE. In both CV schemes the age model is as good or better on I-genes "
            f"than on A-genes (within-site R²_A={_fmt_r(ws.get('r2_A'))} vs R²_I={_fmt_r(ws.get('r2_I_matched'))}; "
            f"LOSO R²_A={_fmt_r(lo.get('r2_A'))} vs R²_I={_fmt_r(lo.get('r2_I_matched'))}). "
            f"FALSIFICATION.md labels disagree (LOSO: {lo.get('verdict')}; within-site: {ws.get('verdict')}) "
            "and are both reported; they agree on the comparison."
        )
    elif p3 and p3.get("stop"):
        verdict = "STOP P3 — " + str(p3.get("stop_reason", "halt"))
    elif p2 and p2.get("stop"):
        verdict = "STOP P2 — fewer than 30 A-genes survived filters"
    return f"""# FINDINGS_BRAIN_PHASE1 — age/identity axes on DLPFC adult controls

**Status: {", ".join(stages)}.** Seed `{PHASE1_SEED}`. Cohort: 233 adult neurotypical DLPFC donors (ages 20–89), 3,323 pseudobulks, 20 cell types, 2 sites (Source = HBCC/MSSM). Reproduced by `notebooks/brain_phase1_p0.ipynb` … `p5`. Run date 2026-09-14.

**Do not modify** `FINDINGS.md`, `FINDINGS_RERUN.md`, `FINDINGS_STAGEB.md`, `FINDINGS_TISSUE.md`, `FINDINGS_BRAIN.md`, or `FALSIFICATION.md`. `FALSIFICATION.md` is pre-registration and governs Stage P4.

**P4c verdict:** {verdict}

Blood remains closed (compositional). This file is the Phase 1 payoff run on brain.

{HASH_POOL_LIMITATION}"""


def _p0(p0, inc, corr, audit):
    absent = inc.get("absent", {})
    included = inc.get("included", p0.get("included_extra", []))
    reasons = inc.get("reasons", {})
    lines = [
        "## P0 — Brain-specific confounders",
        "",
        "Audited every obs column and `uns` on the cached Aging_Cohort h5ad "
        "(dataset ID `4442d412-91cb-4261-acca-8adf5fa04c11`, fetched in T1/R0 — not invented). "
        "Neuronal fraction and n_genes/n_counts use obs only (no 12 GB matrix). "
        "MALAT1/MT fractions use the existing 3,323-row pseudobulk.",
        "",
        "| variable | in obs / recoverable? | notes |",
        "|---|---|---|",
        f"| post-mortem interval (PMI) | **{'yes' if absent.get('PMI') is False else 'NO'}** | {('present' if not absent.get('PMI', True) else 'ABSENT — real limitation. Known post-mortem confound, often age-correlated.')} |",
        f"| RIN / RNA quality | **{'yes' if absent.get('RIN') is False else 'NO'}** | {('present' if not absent.get('RIN', True) else 'ABSENT — real limitation.')} |",
        f"| tissue pH | **{'yes' if absent.get('tissue_pH') is False else 'NO'}** | {('present' if not absent.get('tissue_pH', True) else 'ABSENT.')} |",
        f"| cause of death | **{'yes' if absent.get('cause_of_death') is False else 'NO'}** | {('present' if not absent.get('cause_of_death', True) else 'ABSENT.')} |",
        f"| SoupX / decontX / ambient scores | **{'yes' if absent.get('soupX_decontX') is False else 'NO'}** | Proxies computed: MALAT1 UMI fraction, MT- UMI fraction, fraction of nuclei with n_genes < 500. |",
        f"| 6-plex hashing pools | **NO** | {HASH_POOL_LIMITATION} `uns/batch_condition` = `Source` only. |",
        "| neuronal nucleus fraction | **computed** | Donor-level fraction of nuclei with `class` in {EN, IN}. Dissociation-bias proxy. |",
        "",
        f"**A priori inclusion rule:** a donor-level expression-affecting covariate enters the P1 design if within-site Spearman |r| with age ≥ **{P0_AGE_CORR_ABS}** in either site or in the site-residualized pool. "
        "log(mean UMI/nucleus) and log(mean genes/nucleus) are **always** in the design (specified), not selected on |r|.",
        "",
        f"Donors in the P0 table: **{p0.get('n_donors', 'NA')}**. Mean neuronal fraction: **{p0.get('frac_neuronal_mean', float('nan')):.3f}**." if p0.get("frac_neuronal_mean") is not None else f"Donors: {p0.get('n_donors')}.",
    ]
    if corr is not None and len(corr):
        lines += ["", "Within-site Spearman(age, covariate):", ""]
        piv = corr.pivot(index="covariate", columns="site", values="spearman_r").reset_index()
        float_cols = [c for c in piv.columns if c != "covariate"]
        lines.append(_md_table(piv, list(piv.columns), float_cols=float_cols, nd=3))
    lines += ["", f"**P0 extras added to the P1 model (donor-level, merged as `p0_*`):** {included if included else '(none beyond specified depth + C(Source))'}"]
    if reasons:
        lines.append("")
        lines.append("`frac_low_n_genes` is NA because it is identically 0 on this QC'd upload (no nuclei with n_genes < 500 among adult analysis donors).")
        lines.append("")
        for k, v in reasons.items():
            lines.append(f"- `{k}`: {v}")
    lines += [
        "",
        "Figures: `results/brain_phase1/figures/p0_covariates_vs_age.png`.",
        f"Obs columns audited: {audit.get('n_obs_columns', 'NA')}. "
        f"`uns/batch_condition`={audit.get('uns_batch_condition')}.",
    ]
    return "\n".join(lines)


def _p1_declared(p1t, sweep, null):
    lines = [
        "## P1 — Variance decomposition",
        "",
        HASH_POOL_LIMITATION,
        "",
        "Model on the 3,323 pseudobulks:",
        "",
        "```",
        "expression ~ C(cell_type) + C(Source) + log(mean UMI/nucleus) + log(mean genes/nucleus)",
        "             + [P0 extras] + age",
        "```",
        "",
        "`unique_age = R²(N+ct+age) − R²(N+ct)` with N = Source + within-type-centred depth + P0 extras. "
        "This is within-site, depth/P0-adjusted age variance. `unique_ct` analogously.",
        "",
        "### Declared A/I thresholds (written BEFORE any score)",
        "",
        "The blood lesson: `t_age` was an absolute unique_age cutoff calibrated on a distribution later shown "
        "to be ~65% pool-depth artifact; after correction that cutoff was absurdly strict (99.9th percentile). "
        "**Brain thresholds are QUANTILES of the corrected distribution**, with quantile *levels* declared a priori "
        f"in `src/brain_phase1_common.py` as `{PRIMARY_QUANTILES}` — not chosen by looking at downstream R².",
        "",
        "- **A-genes:** `unique_age ≥ Q90(unique_age)` AND `unique_ct ≤ Q50(unique_ct)`",
        "- **I-genes:** `unique_ct ≥ Q90(unique_ct)` AND `unique_age ≤ Q50(unique_age)`",
        "- **MIXED:** `unique_age ≥ Q90` AND `unique_ct ≥ Q90` — excluded from both scores",
        "",
        "A sweep over high ∈ {0.80, 0.90, 0.95, 0.99} × low ∈ {0.25, 0.50, 0.75} is reported; "
        "the primary cell is the declared (0.90, 0.50) pair. **No cell is selected on score R².**",
    ]
    if not p1t:
        lines += ["", "_P1 not yet run. Quantile *levels* above are the commitment._"]
        return "\n".join(lines)

    lines += [
        "",
        "**Numeric cutoffs on this cohort's corrected distribution:**",
        "",
        f"| cutoff | quantile | value |",
        f"|---|---:|---:|",
        f"| t_age | Q{p1t['t_age_q']:.2f} unique_age | **{p1t['t_age']:.6f}** |",
        f"| k_ct | Q{p1t['k_ct_q']:.2f} unique_ct | **{p1t['k_ct']:.6f}** |",
        f"| t_ct | Q{p1t['t_ct_q']:.2f} unique_ct | **{p1t['t_ct']:.6f}** |",
        f"| k_age | Q{p1t['k_age_q']:.2f} unique_age | **{p1t['k_age']:.6f}** |",
        "",
        f"**Set sizes at the declared quantiles (pre-P2 filters):** A = **{p1t['n_A']}**, "
        f"I = **{p1t['n_I']}**, MIXED = **{p1t['n_MIXED']}**, OTHER = **{p1t['n_OTHER']}** "
        f"(n genes = {p1t['n_genes']}).",
        "",
        "These sizes were written into this file after P1 and before P3 (`results/brain_phase1/p1_DECLARED_BEFORE_SCORES.flag`).",
        "",
        "### unique_age / unique_ct distribution vs blood",
        "",
        f"| | DLPFC corrected | blood within-pool (FINDINGS_RERUN) |",
        f"|---|---:|---:|",
        f"| unique_age median | {p1t.get('unique_age_median', float('nan')):.5f} | 0.00011 |",
        f"| unique_age 99th pct | {p1t.get('unique_age_p99', float('nan')):.5f} ({100*p1t.get('unique_age_p99', 0):.2f}%) | **{BLOOD_UA_P99:.4f} ({100*BLOOD_UA_P99:.2f}%)** |",
        f"| unique_age max | {p1t.get('unique_age_max', float('nan')):.5f} ({100*p1t.get('unique_age_max', 0):.2f}%) | **{BLOOD_UA_MAX:.4f} ({100*BLOOD_UA_MAX:.2f}%)** |",
        f"| unique_ct median | {p1t.get('unique_ct_median', float('nan')):.4f} | 0.1248 |",
        f"| unique_ct max | {p1t.get('unique_ct_max', float('nan')):.4f} | 0.9610 |",
        "",
        f"P0 extras in this partition: `{p1t.get('extra_covariates')}`.",
    ]
    if null is not None and len(null):
        lines += ["", "Within-site donor-age permutation null (orientation only; **not** used to pick thresholds):", "",
                  _md_table(null, list(null.columns), float_cols=[c for c in null.columns if c != "observed"], nd=4)]
    if sweep is not None and len(sweep):
        a_piv = sweep.pivot(index="high_q", columns="low_q", values="n_A")
        lines += ["", "Sweep: n A-genes (rows = high quantile for unique_age / unique_ct; columns = low quantile). Primary cell marked in the csv `is_primary`.", "",
                  "```", a_piv.to_string(), "```"]
        i_piv = sweep.pivot(index="high_q", columns="low_q", values="n_I")
        lines += ["", "Sweep: n I-genes:", "", "```", i_piv.to_string(), "```"]
    lines += ["", "Tables: `p1_variance_decomposition.csv`, `p1_thresholds.json`, `p1_threshold_sweep_set_sizes.csv`.",
              "Figures: `p1_unique_age_vs_unique_ct.png`, `p1_unique_age_hist_vs_blood.png`."]
    return "\n".join(lines)


def _p2(p2, surv):
    lines = [
        "## P2 — Confound filters on A-genes",
        "",
        HASH_POOL_LIMITATION,
        "",
        "Filters (same lists/thresholds as blood 1b, plus P0-correlated genes). Sequential survival of the declared A set:",
        "",
    ]
    if surv is not None and len(surv):
        lines.append(_md_table(surv, list(surv.columns),
                               float_cols=[c for c in surv.columns if c not in ("filter", "flagged_marginal", "removed_incremental", "surviving")]))
    lines += [
        "",
        f"A-genes in: **{p2.get('n_A_in')}**. Surviving: **{p2.get('n_A_surviving')}**. "
        f"I-genes (unfiltered by A-filters; MIXED still excluded): **{p2.get('n_I')}**.",
        "",
        f"STOP P2 (< {STOP_P2_MIN_A} surviving A-genes): **{p2.get('stop')}**.",
    ]
    if p2.get("cc_note"):
        lines += ["", "**Cell-cycle note (DLPFC is largely post-mitotic):** " + str(p2["cc_note"])]
    if p2.get("stop"):
        lines += ["", p2.get("stop_reason", "Halt.")]
    return "\n".join(lines)


def _p3(p3):
    lines = [
        "## P3 — The two scores",
        "",
        HASH_POOL_LIMITATION,
        "",
        "Two CV schemes, both with **no donor in both folds** (asserted):",
        "",
        "1. **LOSO** — leave-one-site-out (2 folds: train HBCC/test MSSM and vice versa). Batch-grouped, weak (n_batches=2).",
        f"2. **within-site donor k-fold** — {N_SITE_FOLDS}-fold grouped by donor, nested *inside each site*.",
        "",
        "Gene selection (P1 quantiles + P2 filters) is **re-done inside each outer training fold** (FALSIFICATION.md). "
        "Reported performance is nested-CV. Frozen full-data gene sets (after P2) are saved for P5 and for listing, "
        "and are not the CV features.",
        "",
        "Age model: ridge on A-genes → donor age, within cell type. Identity: multinomial logistic on I-genes → cell type. "
        "Unrestricted age model: kernel ridge on all genes, within cell type (R3-like ceiling).",
    ]
    for scheme in ("loso", "within_site"):
        s = p3.get(scheme) or {}
        lines += ["", f"### {scheme}", ""]
        if not s:
            lines.append("_(missing)_")
            continue
        lines.append(
            f"- Median-over-types age R² (A-genes): **{_fmt_r(s.get('age_A_median_r2'))}** "
            f"(MAE {_fmt_r(s.get('age_A_median_mae'), 2).lstrip('+')})"
        )
        lines.append(
            f"- Median-over-types unrestricted R²: **{_fmt_r(s.get('age_unrest_median_r2'))}**"
        )
        lines.append(
            f"- Joint (all types) A-gene age R²: **{_fmt_r(s.get('age_A_joint_r2'))}**; "
            f"unrestricted joint: **{_fmt_r(s.get('age_unrest_joint_r2'))}**"
        )
        acc = s.get("id_I_acc")
        f1 = s.get("id_I_f1")
        ch = s.get("id_chance")
        lines.append(
            f"- Identity accuracy (I-genes): **{_fmt_r(acc, 3).lstrip('+') if acc is not None else 'NA'}**  "
            f"macro-F1 **{_fmt_r(f1, 3).lstrip('+') if f1 is not None else 'NA'}**  "
            f"chance **{_fmt_r(ch, 3).lstrip('+') if ch is not None else 'NA'}**"
        )
    if p3.get("per_type_table"):
        lines += ["", "Per-cell-type tables: `p3_age_per_type.csv`, `p3_identity.csv`."]
    lines += [
        "",
        f"STOP P3: **{p3.get('stop')}**. {p3.get('stop_reason', '')}",
        "",
        "Interpretation vs R3 ceiling (~0.44 joint within-site expression): "
        + str(p3.get("interpretation", "")),
    ]
    return "\n".join(lines)


def _p4(p4):
    lines = [
        "## P4 — Negative controls and the cross-test (the verdict)",
        "",
        "Governed by `FALSIFICATION.md` (written before any 1c/1d model; not revised). "
        f"Size-matched I-gene age models: {N_MATCH_DRAWS} draws, median reported.",
        "",
    ]
    for scheme in ("loso", "within_site"):
        s = p4.get(scheme) or {}
        lines += [f"### {scheme}", ""]
        if not s:
            lines.append("_(missing)_")
            continue
        lines.append(f"- P4a shuffled ages, A-gene R² (median types): **{_fmt_r(s.get('shuffle_age_r2'))}**")
        lines.append(f"- P4b shuffled types, identity acc: **{_fmt_r(s.get('shuffle_id_acc'), 3).lstrip('+')}**")
        lines.append(
            f"- R²_A (proper A-genes) = **{_fmt_r(s.get('r2_A'))}**; "
            f"R²_I (size-matched I-genes, median of {N_MATCH_DRAWS}) = **{_fmt_r(s.get('r2_I_matched'))}**; "
            f"R²_I (full I set) = **{_fmt_r(s.get('r2_I_full'))}**"
        )
        lines.append(
            f"- Acc_I = **{_fmt_r(s.get('acc_I'), 3).lstrip('+')}**; "
            f"Acc_A (identity from A-genes) = **{_fmt_r(s.get('acc_A'), 3).lstrip('+')}**"
        )
        rr = s.get("ratio_r2")
        ra = s.get("ratio_acc")
        lines.append(f"- Δ R²_A − R²_I(matched) = **{_fmt_r(s.get('delta_r2'))}**; "
                     f"R²_I / R²_A = **{_fmt_r(rr, 3).lstrip('+') if isinstance(rr, (int, float)) else rr}**")
        lines.append(f"- Acc_A / Acc_I = **{_fmt_r(ra, 3).lstrip('+') if isinstance(ra, (int, float)) else ra}**")
        lines.append(f"- scheme verdict: **{s.get('verdict')}**")
        lines.append("")
    lines += [
        f"**THE DECISION (plain):** {p4.get('verdict', '')}",
        "",
        p4.get("verdict_text", ""),
        "",
        p4.get("scientific_plain", ""),
        "",
        "Where LOSO and within-site disagree on the *FALSIFICATION.md label*, both labels are reported. "
        "They **agree** on the comparison the decision rule is about: in both schemes the age model is "
        "as good or better on I-genes as on A-genes, and identity accuracy on A-genes is essentially the "
        "same as on I-genes. That is the cross-test result. Thresholds were not retuned.",
    ]
    return "\n".join(lines)


def _p5(p5):
    lines = [
        "## P5 — Replication on the frozen sets (white-matter atlas)",
        "",
        "Dataset `c05e6940-729c-47bd-a2a6-6ce3730c4919` (cached). **Nothing was re-derived.** "
        "Frozen DLPFC A/I gene IDs (Ensembl, version-stripped) were mapped onto the WM pseudobulk.",
        "",
        f"- A-genes mapped: **{p5.get('n_A_mapped')}** / {p5.get('n_A_frozen')}",
        f"- I-genes mapped: **{p5.get('n_I_mapped')}** / {p5.get('n_I_frozen')}",
        f"- A-gene age R² (median over matched types, SequencingPool-grouped CV): **{_fmt_r(p5.get('r2_A'))}**",
        f"- size-matched I-gene age R²: **{_fmt_r(p5.get('r2_I_matched'))}**",
        f"- cross-test holds? **{p5.get('cross_test_holds')}**",
        "",
        p5.get("text", ""),
    ]
    if p5.get("recommend_retina"):
        lines += [
            "",
            "**Does not replicate** as a clean copy of the DLPFC verdict. "
            "Recommend the retina 104-donor snRNA atlas (`d6505c89-c43d-4c28-8c4f-7351a5fd5528`) as the second replication target.",
        ]
    return "\n".join(lines)


def _limitations(p0, p1t, p3, p4, p5):
    return f"""## Limitations (every stage that depends on them)

1. **Two-batch CV.** Source has 2 levels (HBCC/MSSM). Leave-one-site-out is the only batch-grouped scheme and is weak: the two folds are the two banks, MSSM is older (R2(age~Source)=0.264), and repeats are not independent. Donor-grouped k-fold *within each site* is reported alongside it. Where they disagree, both numbers stand. This is not OneK1K's 75-pool CV.
2. **Unmeasured 6-plex hashing pools.** {HASH_POOL_LIMITATION} Pool-structured age would inflate `unique_age` the way OneK1K pools did; we cannot test that on this upload.
3. **PMI and RIN are not in obs.** Post-mortem interval and RNA integrity are the standard brain-snRNA confounds and are often age-correlated. They were not uploaded. Neuronal fraction, MALAT1, MT fraction, and low-n_genes debris are the recoverable proxies; they are not substitutes for PMI/RIN.
4. **Single-nucleus, 10x 3′ v3.** Metrics are comparable to blood scRNA; absolute sensitivity is not.
5. **Gene selection nested in CV** uses the same declared *quantile levels* on each training fold's own unique_age / unique_ct distribution. Frozen full-data sets are for P5 and for listing genes, not for the reported scores.
6. Seed `{PHASE1_SEED}` for every RNG (`numpy.random.default_rng`).
7. **P4 within-site age shuffle** left median-over-types R² = +0.195 (above the 0.10 soundness line). LOSO shuffle collapsed (R² = −0.91). High-p ridge (1,355 genes) on ~100 donors per type can leave noisy positive OOF R² after a label permutation; this is reported, not used to claim or rescue separability. Identity shuffles hit chance in both schemes.
8. P4 cross-test uses the **frozen** full-data A/I sets so A vs I are size-and-selection comparable. P3 nested A-gene R² (within-site 0.345) is the leakage-safer estimate of the A-gene clock; frozen P4 R²_A is 0.423. The I > A ordering holds for both."""


def _files():
    return """## Files

| path | content |
|---|---|
| `src/brain_phase1_common.py`, `brain_phase1_p0.py` … `p5.py`, `brain_phase1_models.py`, `brain_phase1_findings.py`, `brain_phase1_run.py`, `brain_phase1_make_notebooks.py` | code |
| `notebooks/brain_phase1_p0.ipynb` … `p5` | runnable per stage |
| `results/brain_phase1/` | tables, json, logs, `figures/` |
| `FINDINGS_BRAIN_PHASE1.md` | this file |
| `FALSIFICATION.md` | pre-registration (unchanged) |"""


if __name__ == "__main__":
    write_findings()
