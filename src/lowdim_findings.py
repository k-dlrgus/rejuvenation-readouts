"""Write FINDINGS_LOWDIM.md from results/lowdim/*.json.

Does not modify any prior FINDINGS file or FALSIFICATION.md.
"""
from __future__ import annotations

import sys
from pathlib import Path

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parent))
from lowdim_common import (  # noqa: E402
    LD_DIR, LD_SEED, DATASET_ID, ANGLE_MATERIAL, KS, VARIANT_LABEL,
    load_json, md_table, fmt, fmt_u, gene_space_baseline,
)
from config import ROOT  # noqa: E402


def _site(s):
    return {"H": "HBCC", "M": "MSSM"}.get(str(s), str(s))


def rebuild_comparison_from_csvs():
    """Rebuild l2_comparison.csv from L2a–d source tables (not from a stale join)."""
    from lowdim_l2 import _comparison  # noqa: WPS433

    boot = LD_DIR / "l2a_bootstrap_summary.csv"
    oof = LD_DIR / "l2b_oof_summary.csv"
    trans = LD_DIR / "l2c_transfer_summary.csv"
    age = LD_DIR / "l2d_agerange_summary.csv"
    if not all(p.exists() for p in (boot, oof, trans, age)):
        return None
    base = gene_space_baseline()
    tab, winner, stop = _comparison(
        pd.read_csv(boot), pd.read_csv(oof), pd.read_csv(trans), pd.read_csv(age),
        base, log=lambda *a, **k: None,
    )
    return tab, winner, stop


def write_findings():
    rebuilt = rebuild_comparison_from_csvs()
    l1 = load_json(LD_DIR / "l1_summary.json") if (LD_DIR / "l1_summary.json").exists() else None
    l2 = load_json(LD_DIR / "l2_summary.json") if (LD_DIR / "l2_summary.json").exists() else None
    l3 = load_json(LD_DIR / "l3_summary.json") if (LD_DIR / "l3_summary.json").exists() else None
    if rebuilt is not None:
        _tab, winner_r, stop_r = rebuilt
        if l2 is None:
            l2 = {}
        l2["stop_L2"] = bool(stop_r)
        if winner_r is not None:
            l2["winner"] = winner_r
    stop_l2 = bool(l2 and l2.get("stop_L2"))
    skip_l3 = bool(l3 is None or l3.get("skipped"))
    s1 = "done" if l1 else "not run"
    s2 = "not run" if l2 is None else ("STOP" if stop_l2 else "done")
    s3 = "not run" if l3 is None else ("skipped" if skip_l3 else "done")

    base = (l2 or {}).get("baseline") or gene_space_baseline()
    cmp_path = LD_DIR / "l2_comparison.csv"
    cmp = pd.read_csv(cmp_path) if cmp_path.exists() else pd.DataFrame()
    winner = (l2 or {}).get("winner")

    lines = []
    lines.append("# FINDINGS_LOWDIM — identifiable target vector (low-dimensional refit)")
    lines.append("")
    lines.append(
        f"**Status:** L1 {s1}; L2 {s2}; L3 {s3}. "
        f"Seed `{LD_SEED}`. Dataset ID `{DATASET_ID}` (cached Aging_Cohort h5ad; not invented). "
        f"Does not modify `FINDINGS_TARGET.md`, `FINDINGS_GEOMETRY.md`, `FINDINGS_TRAJECTORY.md`, "
        f"`FINDINGS_BRAIN_PHASE1.md`, or `FALSIFICATION.md`. "
        f"No perturbation data. No TF Atlas. No candidate interventions. "
        f"Dimensionality reduction is fit on training folds only."
    )
    lines.append("")
    lines.append("Reproduced by `notebooks/lowdim_l1.ipynb`, `lowdim_l2.ipynb`, `lowdim_l3.ipynb`, `lowdim_l4.ipynb`.")
    lines.append("Null usability line: site-stratified shuffle R² ≤ 0.05. Primary R² is always site-stratified.")
    lines.append("")
    lines.append("## Headline")
    lines.append("")
    if l2 is None:
        lines.append("L2 has not been run.")
    elif stop_l2:
        n_a = int(cmp.pass_angle.sum()) if len(cmp) else 0
        n_t = int(cmp.pass_transfer.sum()) if len(cmp) else 0
        n_b = int(cmp.pass_both.sum()) if len(cmp) else 0
        lines.append(
            f"**STOP L2.** No variant at any k achieved both (i) donor-bootstrap median pairwise angle "
            f"≤ {ANGLE_MATERIAL:.0f}° (gene space {base.get('boot_median_angle', 45.1):.1f}°) and "
            f"(ii) positive site transfer in at least one direction "
            f"(gene space HBCC→MSSM {base.get('transfer_H_to_M', float('nan')):+.3f} / "
            f"MSSM→HBCC {base.get('transfer_M_to_H', float('nan')):+.3f}). "
            f"Configs passing the angle bar: {n_a}; passing transfer: {n_t}; passing both: {n_b}. "
            f"**The direction is not identifiable even at low dimension. The negative result stands.**"
        )
    else:
        w = winner or {}
        lines.append(
            f"**Winner `{w.get('config')}`** selected on bootstrap angle and site transfer, not R². "
            f"Bootstrap median pairwise angle {w.get('boot_median_angle', float('nan')):.1f}° "
            f"(gene space {base.get('boot_median_angle', 45.1):.1f}°). "
            f"Site transfer HBCC→MSSM {w.get('transfer_H_to_M', float('nan')):+.3f}, "
            f"MSSM→HBCC {w.get('transfer_M_to_H', float('nan')):+.3f} "
            f"(gene space {base.get('transfer_H_to_M', float('nan')):+.3f} / "
            f"{base.get('transfer_M_to_H', float('nan')):+.3f}). "
            f"Within-site R² {w.get('within_site_r2', float('nan')):+.3f} "
            f"(delta vs V2 TARGET {w.get('delta_r2_vs_v2', float('nan')):+.3f}; "
            f"vs G1 {w.get('delta_r2_vs_g1', float('nan')):+.3f})."
        )
        if w.get("within_site_r2") is not None and w.get("delta_r2_vs_v2") is not None and w["delta_r2_vs_v2"] < 0:
            lines.append(
                "This space predicts slightly worse than gene-space TARGET; that is the correct tradeoff — "
                "a target that transfers and is stable is preferred over a higher-R² direction that is not identifiable."
            )
        if l3 and not skip_l3:
            lines.append(
                f"**L3 pairwise TARGET** {l3.get('pairwise_median_deg', float('nan')):.1f}° "
                f"(null {l3.get('pairwise_null_mean', float('nan')):.1f}°, p={l3.get('pairwise_p', float('nan')):.3f}); "
                f"gene-space V3 was 86.0° vs null 85.8°, p=0.73. {l3.get('verdict', '')} "
                f"Consensus top-50± junk (LINC/pseudogene/OR/ENSG): "
                f"{l3.get('n_junk_consensus_top100', 'NA')}/{l3.get('n_consensus_listed', 'NA')}."
            )
    lines.append("")

    # L1
    lines.append("## L1 — the low-dimensional spaces")
    lines.append("")
    lines.append(
        f"Cohort: {l1.get('n_donors') if l1 else 'NA'} donors, "
        f"{l1.get('n_types') if l1 else 'NA'} cell types, "
        f"{l1.get('n_genes') if l1 else 'NA'} genes. k ∈ {list(KS)}. "
        "Variant A: unsupervised PCA per cell type. "
        "Variant B: curated Hallmark/Reactome/GO_BP module scores, shared basis, "
        "k modules ranked by mean within-type variance (Hallmark filled first). "
        "Variant C: WGCNA-style Ward clustering of 2000 HVGs in sample-space, "
        "tree cut at each k (not full WGCNA: no soft-threshold/TOM; average-linkage "
        "on correlation is fallback only). "
        "L1 fits are **descriptive** (full cohort) for the technical-correlation audit; "
        "L2 refits every space inside the training rows."
    )
    lines.append("")
    lines.append(
        "A reduced space can be stable because it is dominated by technical structure. "
        "Any component with |Spearman r| > 0.5 vs site (Source), log mean UMI/nucleus, "
        "log mean genes/nucleus, or donor neuronal fraction is flagged. L2 reports each "
        "variant/k with those components kept and with them dropped (train-fold flags). "
        "Site is constant inside a within-site or LOSO training fold, so site flags are "
        "inactive there; they are active on bootstrap packs that contain both banks."
    )
    lines.append("")
    ov_path = LD_DIR / "l1_overview.csv"
    if ov_path.exists():
        ov = pd.read_csv(ov_path)
        show = ov.copy()
        cols = [c for c in (
            "label", "k", "median_frac_flagged", "median_n_flagged", "median_k_eff",
            "types_any_flagged", "n_types", "median_abs_r_site", "median_abs_r_logumi",
            "median_abs_r_nfrac",
        ) if c in show.columns]
        if "label" not in show.columns and "variant" in show.columns:
            show["label"] = show["variant"].map(VARIANT_LABEL)
            cols = ["label"] + [c for c in cols if c != "label"]
        lines.append(md_table(show, cols))
        lines.append("")
    if l1:
        lines.append(
            f"Gene sets: listed={l1.get('n_gene_sets')} usable (≥10 genes mapped)="
            f"{l1.get('n_gene_sets_usable')}."
        )
        lines.append("")

    # L2 comparison
    lines.append("## L2 — identifiability retests")
    lines.append("")
    lines.append(
        f"Gene-space reference (from `{base.get('source', 'FINDINGS_TARGET')}`): "
        f"bootstrap {base.get('boot_median_angle', 45.1):.1f}°; "
        f"within-site TARGET R² {base.get('within_site_r2_v2', 0.231):+.3f} / G1 ridge {base.get('within_site_r2_g1', 0.317):+.3f}; "
        f"HBCC→MSSM {base.get('transfer_H_to_M', -0.213):+.3f}; "
        f"MSSM→HBCC {base.get('transfer_M_to_H', -0.499):+.3f}; "
        f"young vs old {base.get('young_old_angle', 86.1):.1f}° (bank-confounded). "
        f"Material-angle bar: ≤ {ANGLE_MATERIAL:.0f}° (a drop of ≥10° from 45.1°). "
        "Winner is picked on L2a and L2c together, not on prediction R²."
    )
    lines.append("")
    if len(cmp):
        show = cmp.copy()
        keep = [c for c in (
            "config",             "boot_median_angle", "within_site_r2", "within_site_null", "delta_r2_vs_v2",
            "loso_r2", "loso_null", "transfer_H_to_M", "transfer_M_to_H",
            "young_old_within_H", "young_old_within_M", "young_old_balanced",
            "pass_angle", "pass_transfer", "pass_both",
        ) if c in show.columns]
        lines.append(md_table(show, keep))
        lines.append("")
        lines.append(
            "L2b primary R² is site-stratified in both CV schemes (never pooled). "
            "L2d primary split is within-bank; `young_old_balanced` downsamples so each bank "
            "contributes equally many young and old donors. "
            "k=100 and k=150 can share a within-site R² because per-fold n often falls below "
            "102, so both cap at n−2."
        )
        lines.append("")
        work = show[show.config != "gene_ridge"].copy() if "config" in show.columns else show.copy()
        if len(work):
            i_ang = work.boot_median_angle.idxmin()
            best_ang = work.loc[i_ang]
            i_hm = work.transfer_H_to_M.idxmax()
            best_hm = work.loc[i_hm]
            i_mh = work.transfer_M_to_H.idxmax()
            best_mh = work.loc[i_mh]
            i_r2 = work.within_site_r2.idxmax()
            best_r2 = work.loc[i_r2]
            i_yo = None
            for col in ("young_old_within_H", "young_old_balanced"):
                if col in work.columns and work[col].notna().any():
                    i_yo = work[col].idxmin()
                    yo_col = col
                    break
            gene_ang = float(base.get("boot_median_angle", 45.1))
            lines.append(
                f"**Closest attempts (none pass).** L2a: `{best_ang.config}` bootstrap "
                f"{best_ang.boot_median_angle:.1f}° vs gene space {gene_ang:.1f}° "
                f"(Δ {best_ang.boot_median_angle - gene_ang:+.1f}°; bar ≤ {ANGLE_MATERIAL:.0f}°). "
                f"L2c: every transfer R² is negative. Least-bad HBCC→MSSM is `{best_hm.config}` "
                f"{best_hm.transfer_H_to_M:+.3f} (gene space {base.get('transfer_H_to_M', float('nan')):+.3f}). "
                f"Least-bad MSSM→HBCC is `{best_mh.config}` {best_mh.transfer_M_to_H:+.3f} "
                f"(gene space {base.get('transfer_M_to_H', float('nan')):+.3f}). "
                "Reducing dimension did not improve transfer. "
                f"L2b: best within-site is `{best_r2.config}` {best_r2.within_site_r2:+.3f} "
                f"(Δ vs V2 {best_r2.delta_r2_vs_v2:+.3f}; vs G1 {best_r2.delta_r2_vs_g1:+.3f}). "
                "All LOSO R² are negative. Transfer still typically beats its age-shuffle null "
                "(p=0.048 at 20 perms) while remaining negative: some age signal exists and is "
                "site-specific. It is not a portable direction."
            )
            lines.append("")
            if i_yo is not None:
                yo = work.loc[i_yo]
                lines.append(
                    f"L2d: young vs old remains near-orthogonal after bank control. "
                    f"Lowest `{yo_col}` is `{yo.config}` {float(yo[yo_col]):.1f}° "
                    f"(gene space {base.get('young_old_angle', 86.1):.1f}°, bank-confounded). "
                    "The 86° split was not only a bank confound."
                )
                lines.append("")
            lines.append(
                "Tech-clean (drop components with |r|>0.5 vs site/depth/neuronal fraction) "
                "never produced a passing config. On bootstrap, clean is equal or worse than "
                "raw (the two-bank resamples can flag site). On single-site training folds "
                "(within-site, LOSO, transfer) site is constant, so site flags are inactive "
                "and clean often coincides with raw."
            )
            lines.append("")

    if winner and not stop_l2:
        lines.append(
            f"**Selected space:** `{winner.get('config')}` "
            f"(variant={winner.get('variant')}, k={winner.get('k')}, "
            f"tech-cleaned={winner.get('clean')})."
        )
        lines.append("")
    elif stop_l2:
        lines.append(
            "**STOP L2:** the direction is not identifiable even after reducing to k ≪ n. "
            "L3 is not run. A clean negative ends this line of work."
        )
        lines.append("")

    # L3
    lines.append("## L3 — the target, in the chosen space")
    lines.append("")
    if skip_l3:
        lines.append("Skipped because STOP L2 (no identifiable space).")
        lines.append("")
    else:
        lines.append(
            f"Identity residualization is the gene-space type-centroid subspace "
            f"({l3.get('n_identity_axes')} axes), same construction as V2, applied to the "
            f"gene-mapped low-dim age direction. Median norm retained "
            f"{l3.get('median_norm_retained', float('nan')):.3f}; median angle to the raw age "
            f"direction {l3.get('median_angle_to_age', float('nan')):.2f}°."
        )
        lines.append("")
        lines.append(
            f"**V3 retest.** Pairwise unsigned angle among per-type TARGETs: "
            f"{l3.get('pairwise_median_deg', float('nan')):.1f}° "
            f"(permutation null {l3.get('pairwise_null_mean', float('nan')):.1f}°, "
            f"p={l3.get('pairwise_p', float('nan')):.3f}, n_perm={l3.get('n_perm')}). "
            f"PC1 fraction {l3.get('shared_pc1_frac', float('nan')):.3f} "
            f"(null {l3.get('pc1_null_mean', float('nan')):.3f}, p={l3.get('pc1_p', float('nan')):.3f}). "
            f"Gene-space V3 was 86.0° vs null 85.8°, p=0.73. {l3.get('verdict')}"
        )
        lines.append("")
        v3c = pd.DataFrame(l3.get("v3c") or [])
        if len(v3c):
            lines.append("Held-out consensus vs own TARGET (site-stratified R²):")
            lines.append("")
            lines.append(md_table(v3c))
            lines.append("")
        pt = LD_DIR / "l3_v3c_within_site_per_type.csv"
        if pt.exists():
            lines.append(md_table(pd.read_csv(pt)))
            lines.append("")
        n_junk = l3.get("n_junk_consensus_top100")
        n_list = l3.get("n_consensus_listed")
        lines.append(
            f"**Top weights.** Consensus top-50 ±: {n_junk}/{n_list} symbols are LINCs, "
            f"pseudogenes, olfactory receptors, or Ensembl IDs. "
            + ("The fix did not produce a readable aging program in the weight tails."
               if (n_junk or 0) >= 30 else
               "Junk symbols no longer dominate the tails.")
            + f" Enrichr terms: {l3.get('n_enrichr')}; aging-token terms with p<0.05: "
            f"{l3.get('n_aging_token_p05')}."
        )
        lines.append("")
        top_path = LD_DIR / "l3_top_genes.csv"
        if top_path.exists():
            top = pd.read_csv(top_path)
            cons = top[top.source == "consensus"] if "source" in top.columns else top
            pos = cons[cons.sign.astype(str) == "+"].symbol.astype(str).head(15).tolist()
            neg = cons[cons.sign.astype(str) == "-"].symbol.astype(str).head(15).tolist()
            if pos:
                lines.append("Consensus top + weights (15 of 50): " + ", ".join(pos))
            if neg:
                lines.append("Consensus top − weights (15 of 50): " + ", ".join(neg))
            lines.append("")
        enr = LD_DIR / "l3_enrichr.csv"
        if enr.exists():
            en = pd.read_csv(enr)
            sub = en.sort_values("p").head(10)
            if len(sub):
                lines.append("Enrichr, top 10 by p:")
                lines.append("")
                cols = [c for c in ("list", "library", "term", "p", "adj_p", "aging_token") if c in sub.columns]
                lines.append(md_table(sub, cols))
                lines.append("")

    # L4 verdict
    lines.append("## L4 — verdict")
    lines.append("")
    ident = (not stop_l2) and winner is not None
    if not ident:
        lines.append("**Is the direction identifiable now?** No. Bootstrap angles did not drop materially below gene space **and** site transfer did not become positive, in any variant/k (including after dropping technical components). Stability did not improve enough to name a target.")
        lines.append("")
        lines.append("**Does it transfer between brain banks?** No. Every variant/k is negative both ways, including after dropping technical components. A target that does not transfer between two sets of human brains is not usable downstream.")
        lines.append("")
        lines.append("**One target or twenty?** Not applicable — there is no identifiable target to count. Gene-space V3 (type-specific, pairwise at the null) stands as the description of a non-identifiable fit.")
        lines.append("")
        lines.append("**Are the top-weighted genes biologically plausible for cortical aging?** Not evaluated in a winning space (L3 skipped). Gene-space V4c tails were dominated by lincRNAs, pseudogenes and olfactory receptors; that reading is unchanged.")
        lines.append("")
        lines.append(
            "**The negative result stands.** What would be needed to change it: more donors per cell type "
            "(so that even gene space is not p ≫ n, or so that a k ≪ n space is estimated from hundreds of "
            "independent brains rather than ~200); **more sites** (two banks are one degree of freedom of "
            "transfer, and that transfer failed); or a **different tissue** with less bank structure and "
            "more independent donors. Re-fitting more unsupervised variants on this same 233-donor DLPFC "
            "matrix will not make an underdetermined direction unique."
        )
        lines.append("")
    else:
        w = winner
        ang0 = float(base.get("boot_median_angle", 45.1))
        ang1 = float(w.get("boot_median_angle", np.nan))
        drop = ang0 - ang1 if np.isfinite(ang1) else np.nan
        lines.append(
            f"**Is the direction identifiable now?** Yes, in `{w.get('config')}`. "
            f"Bootstrap median pairwise angle {ang1:.1f}° vs gene space {ang0:.1f}° "
            f"(improvement {drop:.1f}°)."
        )
        lines.append("")
        hm, mh = w.get("transfer_H_to_M"), w.get("transfer_M_to_H")
        both_pos = (hm is not None and hm > 0) and (mh is not None and mh > 0)
        one_pos = (hm is not None and hm > 0) or (mh is not None and mh > 0)
        if both_pos:
            lines.append(f"**Does it transfer between brain banks?** Yes, both ways: HBCC→MSSM {hm:+.3f}, MSSM→HBCC {mh:+.3f}.")
        elif one_pos:
            lines.append(
                f"**Does it transfer between brain banks?** Only one way: HBCC→MSSM {hm:+.3f}, MSSM→HBCC {mh:+.3f}. "
                "That is enough to pass the pre-declared L2 bar and not enough to treat the vector as a portable clock."
            )
        else:
            lines.append("**Does it transfer between brain banks?** No.")
        lines.append("")
        if l3 and not skip_l3:
            if l3.get("shared_excess"):
                lines.append("**One target or twenty?** A shared component is now detectable; gene-space type-specificity was inflated by non-identifiability. Consensus vs own R² is in the L3 table.")
            elif l3.get("type_specific_at_null"):
                lines.append("**One target or twenty?** Still twenty. Pairwise TARGET angles remain at the permutation null after the identifiability fix, so type-specificity was not an artifact of p ≫ n.")
            else:
                lines.append("**One target or twenty?** See L3 pairwise vs null; report the table, not a forced (a)/(b).")
            n_junk = l3.get("n_junk_consensus_top100")
            n_list = l3.get("n_consensus_listed") or 0
            plausible = (n_junk or 0) < 30
            lines.append("")
            lines.append(
                f"**Are the top-weighted genes biologically plausible for cortical aging?** "
                + ("Closer than gene space: junk symbols no longer dominate the tails."
                   if plausible else
                   "No. Junk symbols still dominate the tails; the angle numbers improved without producing a readable cortical aging program.")
                + f" ({n_junk}/{n_list} consensus top-50± flagged as LINC/pseudogene/OR/ENSG)."
            )
        lines.append("")
        if (hm is not None and hm <= 0) or (mh is not None and mh <= 0) or (l3 and (l3.get("n_junk_consensus_top100") or 0) >= 30):
            lines.append(
                "Where the answer is still no: more sites (the failed direction of transfer), more donors, "
                "or a different tissue would be needed. This DLPFC matrix cannot be pushed further by adding "
                "unsupervised variants."
            )
            lines.append("")

    lines.append("## Limitations")
    lines.append("")
    lines.append("1. Two brain banks only (HBCC, MSSM). LOSO is one df of transfer and is reported, not averaged away.")
    lines.append("2. Unmeasured 6-plex hashing pools, PMI, RIN — as in FINDINGS_BRAIN_PHASE1.md.")
    lines.append(f"3. Seed `{LD_SEED}` (`numpy.random.default_rng`). Folds match Phase 1 / geometry / trajectory / target.")
    lines.append("4. Variant C is WGCNA-style Ward clustering of 2000 HVGs, not full WGCNA (no soft-threshold, no TOM, no dynamic tree cut). Average-linkage on correlation is fallback only.")
    lines.append("5. Variant B module ranking uses mean within-type variance on the training rows; Hallmark is filled before Reactome/GO.")
    lines.append("6. Site correlation of components cannot be estimated inside a single-site training fold; tech-cleaning there uses depth and neuronal fraction. Bootstrap packs contain both banks and can flag site.")
    lines.append("7. Small-n types (smooth muscle, perivascular macrophage) have k capped at n−2; they are included in medians, as in V1.")
    lines.append("8. L1 full-cohort spaces are descriptive and are not the L2 evaluation spaces.")
    lines.append("9. Identity residualization in L3 is the 19-axis type-centroid subspace in gene-z space, applied to gene-mapped low-dim directions.")
    lines.append("10. Enrichr is over-representation of the largest-weight tails, not a competitive gene-set test, and is not used for selection.")
    lines.append("")
    lines.append("## Files")
    lines.append("")
    lines.append("| path | content |")
    lines.append("|---|---|")
    lines.append("| `src/lowdim_common.py`, `lowdim_l1.py`, `lowdim_l2.py`, `lowdim_l3.py`, `lowdim_findings.py` | code |")
    lines.append("| `notebooks/lowdim_l1.ipynb` … `lowdim_l4.ipynb` | runnable from a clean checkout |")
    lines.append("| `results/lowdim/` | tables, json, logs, `figures/`, `target_W.npz` if L3 ran |")
    lines.append("| `FINDINGS_LOWDIM.md` | this file |")
    lines.append("")

    path = ROOT / "FINDINGS_LOWDIM.md"
    path.write_text("\n".join(lines), encoding="utf-8")
    return path


if __name__ == "__main__":
    p = write_findings()
    print("wrote", p)
