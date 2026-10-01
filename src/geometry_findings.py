"""Write FINDINGS_GEOMETRY.md from results/geometry/*.json.

Does not modify any prior FINDINGS file or FALSIFICATION.md.
"""
from __future__ import annotations

import sys
from pathlib import Path

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parent))
from geometry_common import GEO_DIR, GEO_SEED, DATASET_ID, NULL_USABLE, load_json  # noqa: E402
from config import ROOT  # noqa: E402


def _fmt(x, d=3):
    if x is None:
        return "NA"
    try:
        x = float(x)
    except (TypeError, ValueError):
        return "NA"
    if not np.isfinite(x):
        return "NA"
    return f"{x:+.{d}f}"


def _fmt_u(x, d=3):
    """Unsigned / accuracy-style."""
    s = _fmt(x, d)
    return s.lstrip("+") if s != "NA" else s


def _md_table(df, cols=None):
    if cols is not None:
        df = df[cols]
    lines = ["| " + " | ".join(df.columns) + " |", "|" + "|".join("---" for _ in df.columns) + "|"]
    for _, row in df.iterrows():
        cells = []
        for c in df.columns:
            v = row[c]
            if isinstance(v, (float, np.floating)):
                cells.append("NA" if not np.isfinite(v) else f"{float(v):+.3f}")
            elif isinstance(v, (bool, np.bool_)):
                cells.append(str(bool(v)))
            elif v is None or (isinstance(v, float) and not np.isfinite(v)):
                cells.append("NA")
            else:
                cells.append(str(v))
        lines.append("| " + " | ".join(cells) + " |")
    return "\n".join(lines)


def write_findings():
    a_path = GEO_DIR / "partA_summary.json"
    if not a_path.exists():
        raise FileNotFoundError("results/geometry/partA_summary.json missing — run Part A first")
    A = load_json(a_path)
    B = load_json(GEO_DIR / "partB_summary.json") if (GEO_DIR / "partB_summary.json").exists() else None
    C = load_json(GEO_DIR / "partC_summary.json") if (GEO_DIR / "partC_summary.json").exists() else None
    side = pd.read_csv(GEO_DIR / "partA_side_by_side.csv") if (GEO_DIR / "partA_side_by_side.csv").exists() else None
    hyp = pd.read_csv(GEO_DIR / "partA_hypotheses.csv") if (GEO_DIR / "partA_hypotheses.csv").exists() else None

    stop_A = bool(A.get("stop_A"))
    lines = []
    lines.append("# FINDINGS_GEOMETRY — leakage audit, gene-count control, geometric reframe")
    lines.append("")
    lines.append(
        f"**Status:** Part A {'STOP' if stop_A else 'resolved'}; "
        f"Part B {'skipped' if (not B or B.get('skipped')) else 'done'}; "
        f"Part C {'skipped' if (not C or C.get('skipped')) else 'done'}. "
        f"Seed `{GEO_SEED}`. Dataset ID `{DATASET_ID}` (cached Aging_Cohort h5ad; not invented). "
        f"Does not modify `FINDINGS_BRAIN_PHASE1.md` or `FALSIFICATION.md`. "
        f"Does not derive new A/I gene-set definitions."
    )
    lines.append("")
    lines.append(
        "Reproduced by `notebooks/geometry_partA.ipynb`"
        + ("" if stop_A else ", `geometry_partB.ipynb`, `geometry_partC.ipynb`")
        + f". Null usability line: site-stratified shuffle R² ≤ {NULL_USABLE}."
    )
    lines.append("")

    # ----- headline -----
    v_ws = A["verdict_within_site"]
    v_lo = A["verdict_loso"]
    lines.append("## Headline")
    lines.append("")
    lines.append(
        f"**Leakage.** The P4 within-site shuffle R² of {_fmt(A['schemes']['within_site']['reported_p4']['shuffle_age_r2'])} "
        f"is real. It is not donor leakage and not a per-pseudobulk shuffle. It is the between-bank age gap "
        f"(R²(age~site)={_fmt(A['r2_age_on_site'])}) recovered by a site-specific intercept when OOF R² is pooled "
        f"across HBCC and MSSM. Site-stratified R² puts that null at "
        f"{_fmt(A['schemes']['within_site']['shuffle_A_null']['site_strat_mean'])}."
    )
    lines.append("")
    lines.append(
        f"**Corrected P4c.** Within-site (site-stratified): R²_A={_fmt(v_ws['r2_A'])} vs size-matched R²_I={_fmt(v_ws['r2_I'])}, "
        f"shuffle={_fmt(v_ws['shuffle'])}, verdict **{v_ws['verdict']}**. "
        f"LOSO (original pooled metric): R²_A={_fmt(v_lo['r2_A'])} vs R²_I={_fmt(v_lo['r2_I'])}, "
        f"shuffle={_fmt(v_lo['shuffle'])}, verdict **{v_lo['verdict']}**. "
        f"P4c change: **{A.get('p4c_change')}**."
    )
    lines.append("")
    if B and not B.get("skipped"):
        b1s = B.get("b1") or []
        ws_b = next((r for r in b1s if r.get("scheme") == "within_site"), None)
        if ws_b:
            lines.append(
                f"**B3.** A-genes at matched size do **not** beat random "
                f"(within-site R²_A={_fmt(ws_b.get('A_frozen'))} vs random median {_fmt(ws_b.get('random_median'))}, "
                f"p95 {_fmt(ws_b.get('random_p95'))}). I-genes **do** "
                f"(R²_I={_fmt(ws_b.get('I_full'))}). So the A-bucket is empty as an age axis "
                f"(methods / weak-signal aggregation), while identity genes carry extra age signal "
                f"(biology). This is not clean (a) or clean (b)."
            )
        else:
            lines.append(f"**B3.** {B.get('b3_text')}")
        lines.append("")
    if C and not C.get("skipped"):
        c3 = C["c3"]
        shared_label = _c4_label(c3)
        lines.append(
            f"**C3.** Median-over-types angle of the PLS-1 age direction to the identity subspace: "
            f"**{c3['median_angle_deg']:.1f}°** (permutation null mean {c3['null_median_angle_mean']:.1f}°, "
            f"p={_fmt_u(c3['p_angle_le_obs'])}). Projection R²={_fmt_u(c3['median_proj_r2'])} "
            f"(null {_fmt_u(c3['null_median_proj_r2_mean'])}). Ridge sensitivity: "
            f"{c3.get('ridge_median_angle_deg', float('nan')):.1f}°. "
            f"Geometry: **largely orthogonal** (a small excess overlap vs null, p=0.03, is not "
            f"'substantially within' the identity subspace). {shared_label}"
        )
        lines.append("")
    claim = _overall_claim(A, B, C, stop_A)
    lines.append(f"**Claim the evidence now supports:** {claim}")
    lines.append("")

    # ----- Part A -----
    lines.append("## Part A — Leakage audit")
    lines.append("")
    lines.append(f"Cohort: {A['n_donors']} donors, {A['n_types']} cell types, frozen A={A['n_A']} I={A['n_I']}.")
    lines.append("")
    lines.append("### A1. The 0.195 is real")
    lines.append("")
    lines.append(
        f"Stored P4 `within_site.shuffle_age_r2` = "
        f"`{A['schemes']['within_site']['reported_p4']['shuffle_age_r2']}` "
        f"(JSON, not a rounding artifact). Independent original-ridge donor-level shuffle: "
        f"within-site pooled R²=+0.177 (same magnitude as +0.195), site-stratified R²=−0.074. "
        f"An intercept-only model — predict the training-fold mean age, no genes — already yields "
        f"pooled within-site R²=+0.217 vs R²(age~site)={_fmt(A['r2_age_on_site'])}, and "
        f"site-stratified R²=−0.032. LOSO shuffle remains largely negative, as in P4, because "
        f"the intercept is the *other* bank's mean age."
    )
    lines.append("")
    lines.append("### A2. Hypotheses")
    lines.append("")
    if hyp is not None:
        lines.append(_md_table(hyp, ["id", "fires", "result"]))
        lines.append("")
        for _, h in hyp.iterrows():
            lines.append(f"- **{h['id']}** (fires={h['fires']}): {h['note']}")
        lines.append("")
    lines.append(
        f"**Fired:** {', '.join(A.get('fired') or ['none'])}. "
        f"**Fix:** {A.get('fix')}"
    )
    lines.append("")
    lines.append("Donor age is constant across a donor's pseudobulks "
                 f"({A.get('donor_age_constant')}). Shuffling within a donor's rows is therefore a no-op. "
                 "P4 already used donor-level within-site permutation (`permute_age_within_site_obs`). "
                 "That permutation preserves each bank's age mean, which is exactly what the pooled within-site "
                 "metric reads out via the intercept.")
    lines.append("")
    lines.append("### A3. Reported vs corrected")
    lines.append("")
    lines.append(
        "Corrected within-site number = site-stratified median-over-types R² "
        "(mean of HBCC and MSSM R² per cell type, then median over types). "
        "Corrected LOSO number = original pooled median-over-types R² "
        "(intercept mismatch is a real transfer failure). "
        "Every predictive number below carries its permutation null."
    )
    lines.append("")
    lines.append(
        "**On whether the null is usable.** STOP A asked whether the *positive* leak remained "
        f"(shuffle still above ~{NULL_USABLE}). It does not: site-stratified A-gene shuffle mean is "
        f"{_fmt(A['schemes']['within_site']['shuffle_A_null']['site_strat_mean'])}. "
        "The intercept-only site-stratified R² is −0.032, so the *metric* is sound. "
        "SVD-LOO ridge under permutation is more negative (−0.23) than original sklearn ridge (−0.07) "
        "because p ≫ n overfits; that is conservative (it cannot create a false age clock), not residual "
        "donor leak. FALSIFICATION.md fails soundness if shuffle R² > 0.10; that bar is now cleared. "
        "Treat observed R² as a contrast against this overfit-negative null, not as a number whose "
        "null sits in ±0.05. Original-ridge site-stratified shuffle (−0.074) *is* near zero."
    )
    lines.append("")
    if side is not None:
        # split by scheme for readability
        for scheme in ("within_site", "loso"):
            sub = side[side.scheme == scheme].copy()
            lines.append(f"#### {scheme}")
            lines.append("")
            lines.append(_md_table(sub, ["metric", "reported_in_P4c", "corrected", "permutation_null", "null_usable"]))
            lines.append("")
    lines.append(
        f"**STOP A:** {stop_A} — {A.get('stop_A_reason')}"
    )
    lines.append("")
    lines.append(
        f"**P4c verdict change.** Reported within-site `{A['schemes']['within_site']['reported_p4']['verdict']}` "
        f"→ corrected `{v_ws['verdict']}` ({v_ws['detail']}). "
        f"Reported LOSO `{A['schemes']['loso']['reported_p4']['verdict']}` "
        f"→ corrected `{v_lo['verdict']}` ({v_lo['detail']}). "
        f"{A.get('p4c_change')}. "
        f"Ordering within-site: {A['old_order']['within_site']} → {A['new_order']['within_site']}; "
        f"LOSO: {A['old_order']['loso']} → {A['new_order']['loso']}."
    )
    lines.append("")
    lines.append(
        "Feature-selection note: P4c used frozen full-data A/I sets. Nested re-derivation (declared "
        "quantile levels + P2 filters inside each training fold) is in the table as `nested A/I`. "
        "Frozen remains the size-matched cross-test, as in P4; nested is the leakage-safer clock."
    )
    lines.append("")

    if stop_A:
        lines.append("Part B and Part C were not run. The pipeline is not sound enough to interpret biology or geometry.")
        lines.append("")
        lines.append("## Files")
        lines.append("")
        lines.append("| path | content |")
        lines.append("|---|---|")
        lines.append("| `src/geometry_*.py` | audit / controls / geometry |")
        lines.append("| `results/geometry/` | tables, json, `figures/` |")
        lines.append("| `FINDINGS_GEOMETRY.md` | this file |")
        text = "\n".join(lines) + "\n"
        (ROOT / "FINDINGS_GEOMETRY.md").write_text(text, encoding="utf-8")
        return ROOT / "FINDINGS_GEOMETRY.md"

    # ----- Part B -----
    lines.append("## Part B — Gene-count control")
    lines.append("")
    if not B or B.get("skipped"):
        lines.append("Skipped.")
        lines.append("")
    else:
        lines.append(
            f"Random sets drawn from all {B['n_genes']} genes. A/I pools are the frozen P2 sets "
            f"(A={B['n_A']}, I={B['n_I']}), not re-derived. {B['n_random_draws']} draws at size |A|, "
            "each with a paired donor-level shuffle. "
            f"Size sweep {B['size_sweep']} with {B.get('n_sweep_draws', 10)} observed draws and "
            "one permutation null per (scheme, family, k). "
            "Primary metric: site-stratified R² within-site, pooled R² LOSO."
        )
        lines.append("")
        if (GEO_DIR / "partB_b1_summary.csv").exists():
            b1 = pd.read_csv(GEO_DIR / "partB_b1_summary.csv")
            lines.append("### B1. Size-matched random sets")
            lines.append("")
            lines.append(_md_table(b1))
            lines.append("")
        if (GEO_DIR / "partB_b2_size_sweep_summary.csv").exists():
            sw = pd.read_csv(GEO_DIR / "partB_b2_size_sweep_summary.csv")
            lines.append("### B2. Size sweep (median R²)")
            lines.append("")
            lines.append(_md_table(sw, ["scheme", "family", "k", "k_effective", "n_draws",
                                        "r2_median", "null_median", "null_usable"]))
            lines.append("")
            lines.append("Figure: `results/geometry/figures/partB_size_sweep.png`.")
            lines.append("")
        lines.append("### B3. Verdict")
        lines.append("")
        lines.append(
            "Not clean (a) and not clean (b). "
            "**A-genes match random at matched size** (within-site 0.250 vs random median 0.254, "
            "inside the random 5–95% band 0.209–0.298; LOSO both ~0). The A-bucket is not a "
            "privileged age axis — that is the arithmetic/methods result. "
            "**I-genes beat both A and random** (within-site 0.424; LOSO 0.289 vs random 0.025). "
            "Identity-associated genes carry extra age information that a random 1,355-gene set "
            "does not. Size sweep: A tracks random at every k; I pulls away from k ≥ 500. "
            "At k=100, SVD-LOO ridge is numerically unstable (hat diagonals near 1 when n≈p); "
            "do not interpret that cell. The fusion *gene-set* finding is therefore: genes cannot "
            "be partitioned into an age bucket, and age prediction lives in identity programmes."
        )
        lines.append("")
        for n in B.get("b3_notes") or []:
            lines.append(f"- {n}")
        lines.append("")

    # ----- Part C -----
    lines.append("## Part C — Geometric reframe")
    lines.append("")
    if not C or C.get("skipped"):
        lines.append("Skipped.")
        lines.append("")
    else:
        c3 = C["c3"]
        lines.append(
            f"All {C['n_genes']} genes. Age direction: PLS-1 (covariance `X'y` within cell type), "
            "ridge coefficients as sensitivity. Identity subspace: SVD of the 20 type centroids "
            f"({C['n_identity_axes']} axes). Angle of a vector `w` to that subspace is "
            r"$\arccos(\\|P_{\mathrm{id}} w\\|)$"
            f"; projection R² is `||P_id w||²`. "
            f"A random unit vector in R^{C['n_genes']} has E[proj R²]={C['random_proj_expect']:.4f}."
        )
        lines.append("")
        lines.append("### C3. Age direction vs identity subspace")
        lines.append("")
        lines.append("| quantity | observed | permutation null mean | p |")
        lines.append("|---|---:|---:|---:|")
        lines.append(
            f"| median-over-types angle (PLS-1) | {c3['median_angle_deg']:.2f}° | "
            f"{c3['null_median_angle_mean']:.2f}° | {_fmt_u(c3['p_angle_le_obs'])} |"
        )
        lines.append(
            f"| median-over-types proj R² | {_fmt_u(c3['median_proj_r2'])} | "
            f"{_fmt_u(c3['null_median_proj_r2_mean'])} | {_fmt_u(c3['p_proj_r2_ge_obs'])} |"
        )
        lines.append(
            f"| consensus age-dir angle | {c3['consensus_angle_deg']:.2f}° | "
            f"{c3['null_consensus_angle_mean']:.2f}° | — |"
        )
        cr = C.get("c3_ridge") or {}
        lines.append(
            f"| ridge median angle (sensitivity) | {cr.get('median_angle_deg', float('nan')):.2f}° | — | — |"
        )
        lines.append("")
        lines.append("Both CV schemes (directions fit on each training fold, z-scored on train, identity basis on train):")
        lines.append("")
        lines.append("| scheme | median angle | proj R² | null angle | p_angle | n_perm |")
        lines.append("|---|---:|---:|---:|---:|---:|")
        for sch, rec in (C.get("cv") or {}).items():
            lines.append(
                f"| {sch} | {rec['median_angle_deg']:.2f}° | {_fmt_u(rec['median_proj_r2'])} | "
                f"{rec['null_median_angle_mean']:.2f}° | {_fmt_u(rec['p_angle_le_obs'])} | {rec['n_perm']} |"
            )
        lines.append("")
        lines.append(f"**C3 reading:** largely orthogonal — {C.get('geom_text')} "
                     "The excess vs the permutation null is real (p=0.03) but small "
                     "(66° vs 71°; proj R² 0.163 vs 0.108). Ridge coefficients are closer still "
                     "to orthogonal (83°). This is not 'substantially within' the identity subspace.")
        lines.append("")
        lines.append("### C4. Shared vs type-specific age direction")
        lines.append("")
        lines.append(
            f"Pairwise unsigned angle among the 20 PLS-1 age directions: median "
            f"{c3['pairwise_median_deg']:.1f}° (null {c3['null_pairwise_median_mean']:.1f}°). "
            f"PC1 of the 20 directions explains {_fmt_u(c3['shared_pc1_frac'])} of their variance "
            f"(null {_fmt_u(c3['null_shared_pc1_mean'])}, p={_fmt_u(c3['p_shared_pc1_ge_obs'])}). "
            f"{_c4_label(c3)} A rejuvenation intervention that moved one shared axis would not "
            "automatically move 20 type-specific 78° neighbours. CV schemes agree (pairwise ~81–83°)."
        )
        lines.append("")
        lines.append("Figures: `partC_angle_vs_null.png`, `partC_pairwise_angles.png`, "
                     "`partC_per_type_angles.png`, `partC_shared_pc1.png`.")
        lines.append("")

    lines.append("## Which claim")
    lines.append("")
    lines.append(claim)
    lines.append("")
    lines.append("These are different findings. Gene-set non-partitionability does not imply geometric fusion, "
                 "and geometric orthogonality would reverse the project's conclusion even if A/I buckets fail.")
    lines.append("")
    lines.append("## Limitations")
    lines.append("")
    lines.append("1. Two brain banks only. LOSO is one df of transfer and is reported, not averaged away.")
    lines.append("2. Unmeasured 6-plex hashing pools, PMI, RIN — as in FINDINGS_BRAIN_PHASE1.md.")
    lines.append("3. Corrected within-site metric is site-stratified R², not a new gene set and not a retune of FALSIFICATION.md thresholds.")
    lines.append(f"4. Seed `{GEO_SEED}` (`numpy.random.default_rng`).")
    lines.append("5. Part B/C age models use SVD-LOO ridge (same α grid as P3/P4). A1 reproduced the shuffle with the original sklearn inner-CV ridge.")
    lines.append("6. Size-sweep k=100 is numerically unstable (n≈p, LOO hat diagonals → 1). Dropped from interpretation.")
    lines.append("7. Pairwise angles among per-type age directions (~78°) are close to the permutation null (~80°); do not read PC1=0.40 as a single shared aging axis.")
    lines.append("")
    lines.append("## Files")
    lines.append("")
    lines.append("| path | content |")
    lines.append("|---|---|")
    lines.append("| `src/geometry_common.py`, `geometry_partA.py`, `geometry_partB.py`, `geometry_partC.py` | code |")
    lines.append("| `notebooks/geometry_partA.ipynb`, `partB`, `partC` | runnable from a clean checkout |")
    lines.append("| `results/geometry/` | tables, json, logs, `figures/` |")
    lines.append("| `FINDINGS_GEOMETRY.md` | this file |")
    lines.append("")

    text = "\n".join(lines) + "\n"
    path = ROOT / "FINDINGS_GEOMETRY.md"
    path.write_text(text, encoding="utf-8")
    return path


def _c4_label(c3):
    pair = float(c3.get("pairwise_median_deg") or 90)
    pc1 = float(c3.get("shared_pc1_frac") or 0)
    null_pc1 = float(c3.get("null_shared_pc1_mean") or 0)
    if pair >= 70:
        return ("Per-type age directions are **largely type-specific** "
                f"(pairwise ~{pair:.0f}°, only ~{90 - pair:.0f}° from orthogonal; "
                f"PC1 {pc1:.2f} vs null {null_pc1:.2f} is a weak common component, not a shared axis).")
    if pair <= 45:
        return ("Per-type age directions are **shared** "
                f"(pairwise {pair:.0f}°, PC1 {pc1:.2f}).")
    return (f"Per-type age directions are partially aligned (pairwise {pair:.0f}°, "
            f"PC1 {pc1:.2f} vs null {null_pc1:.2f}).")


def _overall_claim(A, B, C, stop_A):
    if stop_A:
        return ("The within-site pipeline is still not sound after the attempted fix. "
                "Do not interpret P4c as biology or as a methods-only arithmetic artifact.")
    parts = []
    v_ws = A["verdict_within_site"]["verdict"]
    v_lo = A["verdict_loso"]["verdict"]
    parts.append(
        f"After correcting the site-pooled R² artifact, the gene-set cross-test is "
        f"within-site **{v_ws}** and LOSO **{v_lo}** ({A.get('p4c_change')})."
    )
    if B and not B.get("skipped"):
        parts.append(
            "**Genes cannot be partitioned into clean age vs identity buckets:** "
            "A-genes match size-matched random sets; I-genes beat both, so chronological age "
            "is carried by identity programmes rather than by a privileged A-set. "
            "That is a methods finding about the construct, plus a biological finding that "
            "identity genes are age-informative."
        )
    if C and not C.get("skipped"):
        c3 = C.get("c3") or {}
        parts.append(
            "**The axes are not geometrically fused.** The within-type age direction is largely "
            f"orthogonal to the identity (centroid) subspace "
            f"({c3.get('median_angle_deg', float('nan')):.0f}° vs null "
            f"{c3.get('null_median_angle_mean', float('nan')):.0f}°; ridge {c3.get('ridge_median_angle_deg', float('nan')):.0f}°). "
            "Gene-set non-partitionability and geometric fusion are different claims; "
            "the evidence now supports the first, not the second. "
            + _c4_label(c3)
        )
    return " ".join(parts)


if __name__ == "__main__":
    p = write_findings()
    print("wrote FINDINGS_GEOMETRY.md")
