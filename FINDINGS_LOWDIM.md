# FINDINGS_LOWDIM — identifiable target vector (low-dimensional refit)

**Status:** L1 done; L2 STOP; L3 skipped. Seed `20260914`. Dataset ID `4442d412-91cb-4261-acca-8adf5fa04c11` (cached Aging_Cohort h5ad; not invented). Does not modify `FINDINGS_TARGET.md`, `FINDINGS_GEOMETRY.md`, `FINDINGS_TRAJECTORY.md`, `FINDINGS_BRAIN_PHASE1.md`, or `FALSIFICATION.md`. No perturbation data. No TF Atlas. No candidate interventions. Dimensionality reduction is fit on training folds only.

Reproduced by `notebooks/lowdim_l1.ipynb`, `lowdim_l2.ipynb`, `lowdim_l3.ipynb`, `lowdim_l4.ipynb`.
Null usability line: site-stratified shuffle R² ≤ 0.05. Primary R² is always site-stratified.

## Headline

**STOP L2.** No variant at any k achieved both (i) donor-bootstrap median pairwise angle ≤ 35° (gene space 45.1°) and (ii) positive site transfer in at least one direction (gene space HBCC→MSSM -0.213 / MSSM→HBCC -0.499). Configs passing the angle bar: 0; passing transfer: 0; passing both: 0. **The direction is not identifiable even at low dimension. The negative result stands.**

## L1 — the low-dimensional spaces

Cohort: 233 donors, 20 cell types, 25526 genes. k ∈ [20, 50, 100, 150]. Variant A: unsupervised PCA per cell type. Variant B: curated Hallmark/Reactome/GO_BP module scores, shared basis, k modules ranked by mean within-type variance (Hallmark filled first). Variant C: WGCNA-style Ward clustering of 2000 HVGs in sample-space, tree cut at each k (not full WGCNA: no soft-threshold/TOM; average-linkage on correlation is fallback only). L1 fits are **descriptive** (full cohort) for the technical-correlation audit; L2 refits every space inside the training rows.

A reduced space can be stable because it is dominated by technical structure. Any component with |Spearman r| > 0.5 vs site (Source), log mean UMI/nucleus, log mean genes/nucleus, or donor neuronal fraction is flagged. L2 reports each variant/k with those components kept and with them dropped (train-fold flags). Site is constant inside a within-site or LOSO training fold, so site flags are inactive there; they are active on bootstrap packs that contain both banks.

| label | k | median_frac_flagged | median_n_flagged | median_k_eff | types_any_flagged | n_types | median_abs_r_site | median_abs_r_logumi | median_abs_r_nfrac |
|---|---|---|---|---|---|---|---|---|---|
| A_pca | 20 | +0.050 | +1.000 | +20.000 | 13 | 20 | +0.162 | +0.131 | +0.107 |
| A_pca | 50 | +0.020 | +1.000 | +50.000 | 13 | 20 | +0.100 | +0.111 | +0.093 |
| A_pca | 100 | +0.010 | +1.000 | +100.000 | 13 | 20 | +0.067 | +0.075 | +0.071 |
| A_pca | 150 | +0.007 | +1.000 | +150.000 | 13 | 20 | +0.057 | +0.061 | +0.060 |
| B_curated | 20 | +0.000 | +0.000 | +20.000 | 3 | 20 | +0.120 | +0.101 | +0.073 |
| B_curated | 50 | +0.000 | +0.000 | +50.000 | 6 | 20 | +0.104 | +0.083 | +0.096 |
| B_curated | 100 | +0.000 | +0.000 | +100.000 | 7 | 20 | +0.087 | +0.117 | +0.060 |
| B_curated | 150 | +0.000 | +0.000 | +150.000 | 7 | 20 | +0.083 | +0.107 | +0.059 |
| C_wgcna | 20 | +0.050 | +1.000 | +20.000 | 16 | 20 | +0.105 | +0.079 | +0.091 |
| C_wgcna | 50 | +0.020 | +1.000 | +50.000 | 18 | 20 | +0.101 | +0.076 | +0.078 |
| C_wgcna | 100 | +0.020 | +1.500 | +100.000 | 19 | 20 | +0.100 | +0.071 | +0.072 |
| C_wgcna | 150 | +0.013 | +2.000 | +150.000 | 19 | 20 | +0.099 | +0.076 | +0.071 |

Gene sets: listed=7275 usable (≥10 genes mapped)=4686.

## L2 — identifiability retests

Gene-space reference (from `results/target json`): bootstrap 45.1°; within-site TARGET R² +0.231 / G1 ridge +0.317; HBCC→MSSM -0.213; MSSM→HBCC -0.499; young vs old 86.1° (bank-confounded). Material-angle bar: ≤ 35° (a drop of ≥10° from 45.1°). Winner is picked on L2a and L2c together, not on prediction R².

| config | boot_median_angle | within_site_r2 | within_site_null | delta_r2_vs_v2 | loso_r2 | loso_null | transfer_H_to_M | transfer_M_to_H | young_old_within_H | young_old_within_M | young_old_balanced | pass_angle | pass_transfer | pass_both |
|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|
| gene_ridge | +45.074 | +0.231 | NA | +0.000 | NA | NA | -0.213 | -0.499 | NA | NA | +86.111 | False | False | False |
| A_pca_k020 | +56.854 | +0.061 | -0.094 | -0.170 | -0.806 | -1.513 | -0.669 | -0.903 | +86.633 | +82.980 | +82.572 | False | False | False |
| A_pca_k020_clean | +60.781 | +0.057 | -0.087 | -0.174 | -0.800 | -1.495 | -0.669 | -0.903 | +86.633 | +82.980 | +84.345 | False | False | False |
| A_pca_k050 | +55.544 | +0.191 | -0.147 | -0.040 | -0.473 | -1.556 | -0.255 | -0.603 | +86.578 | +85.873 | +85.420 | False | False | False |
| A_pca_k050_clean | +57.783 | +0.187 | -0.127 | -0.044 | -0.527 | -1.529 | -0.255 | -0.603 | +86.578 | +85.873 | +87.845 | False | False | False |
| A_pca_k100 | +49.564 | -0.147 | -0.244 | -0.378 | -0.386 | -1.598 | -0.225 | -0.525 | +85.652 | +84.726 | +85.219 | False | False | False |
| A_pca_k100_clean | +50.674 | +0.132 | -0.173 | -0.099 | -0.411 | -1.569 | -0.225 | -0.525 | +85.652 | +84.726 | +86.861 | False | False | False |
| A_pca_k150 | +44.611 | -0.147 | -0.244 | -0.378 | -0.344 | -1.610 | -0.259 | -0.525 | +85.652 | +84.726 | +85.219 | False | False | False |
| A_pca_k150_clean | +47.814 | +0.132 | -0.173 | -0.099 | -0.384 | -1.564 | -0.259 | -0.525 | +85.652 | +84.726 | +86.861 | False | False | False |
| B_curated_k020 | +52.571 | -0.102 | -0.120 | -0.333 | -0.969 | -1.571 | -0.872 | -1.146 | +73.834 | +78.158 | +79.705 | False | False | False |
| B_curated_k020_clean | +55.080 | -0.120 | -0.118 | -0.351 | -0.996 | -1.569 | -0.872 | -1.146 | +73.834 | +78.158 | +79.127 | False | False | False |
| B_curated_k050 | +48.505 | -0.119 | -0.177 | -0.350 | -0.859 | -1.604 | -0.802 | -0.989 | +71.673 | +81.386 | +73.979 | False | False | False |
| B_curated_k050_clean | +48.753 | -0.122 | -0.164 | -0.353 | -0.863 | -1.601 | -0.802 | -0.989 | +71.673 | +81.386 | +74.257 | False | False | False |
| B_curated_k100 | +55.583 | -0.104 | -0.181 | -0.335 | -0.874 | -1.565 | -0.742 | -0.934 | +74.886 | +83.367 | +73.440 | False | False | False |
| B_curated_k100_clean | +56.460 | -0.131 | -0.168 | -0.362 | -0.875 | -1.569 | -0.742 | -0.934 | +74.886 | +83.367 | +73.064 | False | False | False |
| B_curated_k150 | +63.146 | -0.104 | -0.181 | -0.335 | -0.858 | -1.567 | -0.601 | -0.933 | +74.886 | +83.367 | +73.440 | False | False | False |
| B_curated_k150_clean | +66.129 | -0.131 | -0.168 | -0.362 | -0.858 | -1.569 | -0.601 | -0.933 | +74.886 | +83.367 | +73.064 | False | False | False |
| C_wgcna_k020 | +83.825 | -0.003 | -0.162 | -0.234 | -0.698 | -1.554 | -0.790 | -0.752 | +85.867 | +85.995 | +85.667 | False | False | False |
| C_wgcna_k020_clean | +84.844 | -0.001 | -0.160 | -0.232 | -0.698 | -1.549 | -0.790 | -0.752 | +85.867 | +85.995 | +85.445 | False | False | False |
| C_wgcna_k050 | +82.235 | -0.033 | -0.263 | -0.264 | -0.694 | -1.629 | -0.759 | -0.644 | +87.493 | +87.261 | +87.412 | False | False | False |
| C_wgcna_k050_clean | +82.785 | -0.034 | -0.247 | -0.265 | -0.587 | -1.607 | -0.759 | -0.644 | +87.493 | +87.261 | +87.553 | False | False | False |
| C_wgcna_k100 | +79.135 | -0.465 | -0.830 | -0.696 | -0.500 | -1.699 | -0.531 | -0.619 | +87.012 | +87.252 | +87.412 | False | False | False |
| C_wgcna_k100_clean | +79.082 | -0.314 | -0.597 | -0.545 | -0.441 | -1.683 | -0.531 | -0.619 | +87.012 | +87.252 | +88.015 | False | False | False |
| C_wgcna_k150 | +84.990 | -0.465 | -0.830 | -0.696 | -0.592 | -1.724 | -0.450 | -0.619 | +87.012 | +87.252 | +87.412 | False | False | False |
| C_wgcna_k150_clean | +84.711 | -0.314 | -0.597 | -0.545 | -0.559 | -1.703 | -0.450 | -0.619 | +87.012 | +87.252 | +88.015 | False | False | False |

L2b primary R² is site-stratified in both CV schemes (never pooled). L2d primary split is within-bank; `young_old_balanced` downsamples so each bank contributes equally many young and old donors. k=100 and k=150 can share a within-site R² because per-fold n often falls below 102, so both cap at n−2.

**Closest attempts (none pass).** L2a: `A_pca_k150` bootstrap 44.6° vs gene space 45.1° (Δ -0.5°; bar ≤ 35°). L2c: every transfer R² is negative. Least-bad HBCC→MSSM is `A_pca_k100` -0.225 (gene space -0.213). Least-bad MSSM→HBCC is `A_pca_k100` -0.525 (gene space -0.499). Reducing dimension did not improve transfer. L2b: best within-site is `A_pca_k050` +0.191 (Δ vs V2 -0.040; vs G1 -0.126). All LOSO R² are negative. Transfer still typically beats its age-shuffle null (p=0.048 at 20 perms) while remaining negative: some age signal exists and is site-specific. It is not a portable direction.

L2d: young vs old remains near-orthogonal after bank control. Lowest `young_old_within_H` is `B_curated_k050` 71.7° (gene space 86.1°, bank-confounded). The 86° split was not only a bank confound.

Tech-clean (drop components with |r|>0.5 vs site/depth/neuronal fraction) never produced a passing config. On bootstrap, clean is equal or worse than raw (the two-bank resamples can flag site). On single-site training folds (within-site, LOSO, transfer) site is constant, so site flags are inactive and clean often coincides with raw.

**STOP L2:** the direction is not identifiable even after reducing to k ≪ n. L3 is not run. A clean negative ends this line of work.

## L3 — the target, in the chosen space

Skipped because STOP L2 (no identifiable space).

## L4 — verdict

**Is the direction identifiable now?** No. Bootstrap angles did not drop materially below gene space **and** site transfer did not become positive, in any variant/k (including after dropping technical components). Stability did not improve enough to name a target.

**Does it transfer between brain banks?** No. Every variant/k is negative both ways, including after dropping technical components. A target that does not transfer between two sets of human brains is not usable downstream.

**One target or twenty?** Not applicable — there is no identifiable target to count. Gene-space V3 (type-specific, pairwise at the null) stands as the description of a non-identifiable fit.

**Are the top-weighted genes biologically plausible for cortical aging?** Not evaluated in a winning space (L3 skipped). Gene-space V4c tails were dominated by lincRNAs, pseudogenes and olfactory receptors; that reading is unchanged.

**The negative result stands.** What would be needed to change it: more donors per cell type (so that even gene space is not p ≫ n, or so that a k ≪ n space is estimated from hundreds of independent brains rather than ~200); **more sites** (two banks are one degree of freedom of transfer, and that transfer failed); or a **different tissue** with less bank structure and more independent donors. Re-fitting more unsupervised variants on this same 233-donor DLPFC matrix will not make an underdetermined direction unique.

## Limitations

1. Two brain banks only (HBCC, MSSM). LOSO is one df of transfer and is reported, not averaged away.
2. Unmeasured 6-plex hashing pools, PMI, RIN — as in FINDINGS_BRAIN_PHASE1.md.
3. Seed `20260914` (`numpy.random.default_rng`). Folds match Phase 1 / geometry / trajectory / target.
4. Variant C is WGCNA-style Ward clustering of 2000 HVGs, not full WGCNA (no soft-threshold, no TOM, no dynamic tree cut). Average-linkage on correlation is fallback only.
5. Variant B module ranking uses mean within-type variance on the training rows; Hallmark is filled before Reactome/GO.
6. Site correlation of components cannot be estimated inside a single-site training fold; tech-cleaning there uses depth and neuronal fraction. Bootstrap packs contain both banks and can flag site.
7. Small-n types (smooth muscle, perivascular macrophage) have k capped at n−2; they are included in medians, as in V1.
8. L1 full-cohort spaces are descriptive and are not the L2 evaluation spaces.
9. Identity residualization in L3 is the 19-axis type-centroid subspace in gene-z space, applied to gene-mapped low-dim directions.
10. Enrichr is over-representation of the largest-weight tails, not a competitive gene-set test, and is not used for selection.

## Files

| path | content |
|---|---|
| `src/lowdim_common.py`, `lowdim_l1.py`, `lowdim_l2.py`, `lowdim_l3.py`, `lowdim_findings.py` | code |
| `notebooks/lowdim_l1.ipynb` … `lowdim_l4.ipynb` | runnable from a clean checkout |
| `results/lowdim/` | tables, json, logs, `figures/`, `target_W.npz` if L3 ran |
| `FINDINGS_LOWDIM.md` | this file |
