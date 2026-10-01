# FINDINGS_SENG — do Sengstack et al. (PNAS 2026) "rejuvenating" TF perturbations move fibroblasts TOWARD young cells?

**Status:** DONE.

## Pre-registration (verbatim, written before any SENG statistic)

Written before any SENG statistic. No threshold in this block is re-tuned after numbers exist. The block below is the task text, verbatim.

```
TASK: SENG — do Sengstack et al. (PNAS 2026) "rejuvenating" TF perturbations
move fibroblasts TOWARD young cells, or just away from old ones, and is the
effect explained by proliferation?
Write FINDINGS_SENG.md, PROGRESS_SENG.md, results/seng/*, src/seng_run.py.
Do not modify any existing FINDINGS_*.md, PROGRESS_*.md, FALSIFICATION.md,
or existing src/ files (import from src/toward_run.py and src/same_run.py;
do not change them). Southard/Hs27 is out of scope. GSE325735 out of scope.
Seeds: 20260914 (splits, nulls, permutations), 20260918 (bootstrap).
n_null=200, n_perm=200, n_boot=200. Uncalibrated R2 never a gate.
=== WRITE THE PRE-REGISTRATION FIRST ===
Write this whole block verbatim into FINDINGS_SENG.md and touch
results/seng/PREREG.flag BEFORE computing any statistic. Stage 0 below
(descriptive inventory only) may run after the flag exists. Nothing is
re-tuned after numbers exist.
=== DATA ===
figshare 10.6084/m9.figshare.30898748 (one h5ad, ~1.43 GB, CC BY 4.0).
Download it to data/raw/seng/. Record file size and md5.
=== STAGE 0: INVENTORY (descriptive only; no statistic) ===
From obs, determine and report:
- which column(s) hold perturbation labels (CRA_<TF>, CRI_<TF>, NT_*), guide
  identity, population doubling (PD14/PD26/PD32 or similar), and WT status.
- which NT labels belong to the CRISPRa line and which to the CRISPRi line.
- which PD the perturbed cells are at.
- cell counts for: every perturbation (pooled over its guides, and per guide),
  NT per modality, WT per PD.
- whether X is raw counts. If X is not raw integer counts and no raw layer
  exists, STOP (key `no_raw_counts`).
Branching rules (fixed now):
- If no WT (untransduced) cells exist at >= 2 distinct PDs: STOP, key
  `no_wt_reference`. Report the inventory only.
- If NT cells cannot be assigned to modality: STOP, key `nt_unassignable`.
- EARLY = WT cells at the lowest PD. LATE = WT cells at the PD of the perturbed
  cells; if that PD cannot be determined, LATE = WT at the highest PD, and say so.
- If EARLY or LATE has < 200 cells: STOP, key `wt_too_small`.
=== GENE SPACE AND TRANSFORM (same as TOWARD/SAME, unchanged) ===
Frozen-ruler overlap gene set; pseudobulk -> TMM log2-CPM prior.count=2 across
this task's full pseudobulk panel -> frozen GTEx mu/sd -> z; missing genes z=0.
Report how many ruler genes are present in this file.
=== SPLIT (frozen before any statistic; seed 20260914) ===
Split EARLY and LATE WT cells each 50/50 into halves A and B.
Half A defines the young direction. Half B is used only in positive control P.
=== STATISTICS (per perturbation g, per modality; MIN_CELLS = 30 pooled) ===
v   = z(EARLY_A) - z(LATE_A)          passage-reversal direction ("young")
d_g = z(cells with g) - z(NT same modality)   the TF's own effect (construct
                                               effect cancels in the subtraction)
cos_g   = dot(d_g, v) / (||d_g|| ||v||)
frac_g  = dot(d_g, v) / ||v||^2
delta_g = ||v - d_g|| - ||v||   (effect applied to a LATE cell; negative =
                                 closer to EARLY. Assumes the TF effect is
                                 additive to the construct effect; state this.)
Proliferation axis (defined BEFORE looking at any perturbation):
- Score every NT cell (both modalities) with the Tirosh 2016 S + G2M gene sets
  (Seurat cc.genes; save the list used to data/reference/).
- p = z(NT cells in top quartile of S+G2M score) - z(NT cells in bottom
  quartile), per modality, then normalize.
- Report cos(v, p): how much of "passage reversal" is proliferation.
- For each g also compute the residual versions: project p out of both v and
  d_g, recompute cos_g_resid and delta_g_resid.
- Report each g's mean S+G2M score shift vs NT.
Donor-age axis: u = u_TOWARD from results/toward/anchors.npz (GTEx young-
minus-old donors). Report cos(v, u) with N1 (permute u's entries across genes,
200x). For each g: cos(d_g, u).
Identity check (report-only flag, not a gate): mean log2 expression of
COL1A1, COL1A2, FN1, LUM, PDGFRA, PDGFRB, POSTN, PRRX1, SERPINH1, VIM in g's
cells minus NT. Flag `identity_loss` if the drop is > 0.5.
=== NULLS ===
N_NT (primary, per modality): for each perturbation's cell count n, draw 200
random subsets of n NT cells as fake perturbations (origin = the remaining NT
cells), compute cos, frac, delta, cos_resid, delta_resid. Empirical one-sided
p for each statistic in the direction of "more toward young".
Multiple testing: BH-FDR across all perturbations within each modality.
Bootstrap: resample cells within g and within NT (B=200) for 95% CIs.
Every reported CI must contain its point estimate or be flagged INVALID.
=== POSITIVE CONTROL P (run BEFORE any perturbation vector) ===
Using half-B WT cells only: mixtures of LATE_B with EARLY_B at young fraction
f = 0, 0.25, 0.50 (cell counts; total fixed). d_f = z(mix_f) - z(LATE_B only).
P PASSES if frac increases monotonically with f, frac(0.50) CI excludes 0,
and delta(0.50) < 0. Report f=0 as the noise floor.
If P fails: STOP, key `pos_control_broken`.
=== REPRODUCTION CHECK (their claim, before our critique) ===
Named hits: CRA_E2F3, CRA_EZH2, CRI_STAT3, CRI_ZFX.
If fewer than 3 of 4 have cos_g > 0 with N_NT FDR q <= 0.05: key
`pipeline_disagrees`. Report all numbers, but do not interpret any hit as
toward/away; say our pipeline does not reproduce their direction.
=== PER-HIT CLASSIFICATION (named hits; only if reproduction passed) ===
Each named hit gets exactly one label, checked in this order:
1. `toward_young`: delta_g < 0 with N_NT q <= 0.05 AND delta_g_resid < 0 with
   N_NT p <= 0.05.
2. `proliferation_explained`: cos_g significant (q <= 0.05) but cos_g_resid
   N_NT p > 0.05.
3. `away_from_old_only`: cos_g significant, delta_g >= 0 or not significant.
4. `no_effect`: none of the above.
Also report per hit: cos(d_g, u) with its N1-style p, identity_loss flag.
=== SCREEN-LEVEL READING (all perturbations with n >= 30) ===
- rho_prolif = Spearman(cos_g, S+G2M shift) across all perturbations, each
  modality separately.
- Key `screen_tracks_proliferation` if rho_prolif >= 0.5 in BOTH modalities.
- Report: number of perturbations with cos_g q <= 0.05 ("their-style hits"),
  and how many of those also meet `toward_young`.
- Report cos(v, u) and whether it beats N1. If it does not, add flag
  `passage_axis_not_donor_age`: aging-in-a-dish does not resemble GTEx donor
  aging in this gene space.
=== WHAT DOES NOT COUNT ===
- cos_g alone as evidence of rejuvenation. That is their metric; ours is delta.
- Using NT from the other modality as origin.
- Using half-B cells outside P, or half-A cells inside P.
- Picking a different EARLY/LATE PD, proliferation gene set, or quartile cut
  after seeing results.
- Dropping or re-labelling perturbations because they look odd. Report them.
- Averaging CRA and CRI results.
- Frozen-ruler score changes: context column only.
- Claims about Sengstack's in vivo mouse results. This data cannot test them.
=== OUTPUT ===
FINDINGS_SENG.md: prereg verbatim; Stage 0 inventory; split sizes; P table;
reproduction check; named-hit table (all statistics, CIs, p, q, labels);
screen-level table and rho_prolif; cos(v,p), cos(v,u); fired keys; limitations
(must include: neonatal cells aged by passaging, not donor aging; one cell
line; additivity assumption in delta; pseudobulk per perturbation).
PROGRESS_SENG.md: STOP status, next action.
Record every failure verbatim. Do not substitute columns or repair rows.
```

Flag: `results/seng/PREREG.flag` exists=True.

## Failures (verbatim)

- **pipeline_disagrees:** Fewer than 3 of 4 named hits have cos_g>0 with N_NT FDR q<=0.05 (1/4). Our pipeline does not reproduce their direction. Named hits are not interpreted as toward/away.

## Stage 0 inventory (descriptive; no statistic)

- File: `<repo>\data\raw\seng\exp3_merged_not_normalized.h5ad`
- Size: 1426871504 bytes. md5: `bfccaa94a916c22e4b0793913a04128b`.
- Figshare API size: 1426871504. URL: https://ndownloader.figshare.com/files/60428597.
- Matrix used: X (csr_matrix).
- Raw-count audit: `{'dtype': 'float32', 'integer_dtype': False, 'n_seen': 148239518, 'n_fractional': 0, 'n_negative': 0, 'n_nonfinite': 0, 'min': 1.0, 'max': 2960.0, 'raw_integer_counts': True, 'where': 'X/data'}`.

### Columns

- Perturbation labels: `Perturbation identity is the pooled id {CRA|CRI}_{TF} over guide columns {TF}_{CRA|CRI}{k}. There is no obs column whose values are the strings CRA_<TF>.`.
- Guide identity: `None` (No single guide-label column. Each guide is a UMI count column named {TF}_CRA{k}, {TF}_CRI{k}, or NT_{k}. A cell is that guide when the column is a unique strict maximum and > 0. Ties with a positive maximum: 2680 (not given a winner). TF-guide calls that disagree with the CRA/CRI indicator: 270 (not assigned). Pooled perturbation id is {CRA|CRI}_{TF} from columns {TF}_{CRA|CRI}{k}.).
- Population doubling: `PD14,PD26,PD32`.
- WT status: `WT column is 0/1 and marks the same cells as exactly one of the boolean PD columns PD14,PD26,PD32. n_WT=6809.`.
- NT assignment: NT_1..NT_6 are shared guide-count columns: a unique-maximum call to each NT_* occurs in both CRA and CRI cells, so the NT label is not a line id. Modality of an NT cell is the CRA or CRI indicator (those two columns do not overlap and cover every non-WT cell). Unique-argmax NT cells: CRA=666 CRI=1479.
- Perturbed-cell PD: perturbed PD values=[]; determined=False; blank=24610
- EARLY: EARLY = WT at the lowest PD PD14 (n=581).
- LATE: The PD of the perturbed cells cannot be determined (perturbed PD values=[], blank PD n=24610). LATE = WT at the highest PD, PD32.

### WT cells per PD

| pd | n_cells |
|---|---|
| PD14 | 581 |
| PD26 | 1880 |
| PD32 | 4348 |

### NT cells per modality

| modality | label | n_cells |
|---|---|---|
| CRA | NT_1 | 123 |
| CRA | NT_2 | 47 |
| CRA | NT_3 | 42 |
| CRA | NT_4 | 125 |
| CRA | NT_5 | 326 |
| CRA | NT_6 | 3 |
| CRI | NT_1 | 393 |
| CRI | NT_2 | 70 |
| CRI | NT_3 | 293 |
| CRI | NT_4 | 161 |
| CRI | NT_5 | 561 |
| CRI | NT_6 | 1 |

### Perturbation counts (pooled over guides)

| modality | perturbation | n_cells |
|---|---|---|
| CRA | CRA_ALX1 | 56 |
| CRA | CRA_APOD | 24 |
| CRA | CRA_ARID5A | 68 |
| CRA | CRA_ARNTL2 | 97 |
| CRA | CRA_ARX | 112 |
| CRA | CRA_ATF4 | 37 |
| CRA | CRA_ATF6 | 27 |
| CRA | CRA_ATF7 | 50 |
| CRA | CRA_ATG3 | 94 |
| CRA | CRA_ATG4C | 37 |
| CRA | CRA_ATG4D | 50 |
| CRA | CRA_ATG5 | 63 |
| CRA | CRA_BATF3 | 85 |
| CRA | CRA_BRCA1 | 76 |
| CRA | CRA_BSX | 51 |
| CRA | CRA_CDX1 | 98 |
| CRA | CRA_CDX2 | 49 |
| CRA | CRA_CEBPA | 63 |
| CRA | CRA_CEBPB | 54 |
| CRA | CRA_CEBPZ | 65 |
| CRA | CRA_CPEB1 | 115 |
| CRA | CRA_CREB1 | 60 |
| CRA | CRA_CREB3L1 | 119 |
| CRA | CRA_CTCF | 69 |
| CRA | CRA_CTCFL | 59 |
| CRA | CRA_DLX2 | 115 |
| CRA | CRA_DLX4 | 63 |
| CRA | CRA_DLX5 | 39 |
| CRA | CRA_DLX6 | 89 |
| CRA | CRA_DNMT1 | 52 |
| CRA | CRA_DRGX | 33 |
| CRA | CRA_DUXA | 24 |
| CRA | CRA_E2F1 | 54 |
| CRA | CRA_E2F3 | 128 |
| CRA | CRA_E2F4 | 66 |
| CRA | CRA_E2F5 | 54 |
| CRA | CRA_E2F7 | 44 |
| CRA | CRA_EGR1 | 34 |
| CRA | CRA_EGR2 | 28 |
| CRA | CRA_EN1 | 26 |
| CRA | CRA_EOMES | 30 |
| CRA | CRA_ESR1 | 3 |
| CRA | CRA_ESX1 | 13 |
| CRA | CRA_ETS1 | 29 |
| CRA | CRA_ETS2 | 119 |
| CRA | CRA_ETV5 | 68 |
| CRA | CRA_EZH2 | 107 |
| CRA | CRA_FERD3L | 83 |
| CRA | CRA_FLI1 | 34 |
| CRA | CRA_FOSL1 | 149 |
| CRA | CRA_FOXA2 | 54 |
| CRA | CRA_FOXA3 | 92 |
| CRA | CRA_FOXB2 | 118 |
| CRA | CRA_FOXC1 | 99 |
| CRA | CRA_FOXD2 | 82 |
| CRA | CRA_FOXF2 | 33 |
| CRA | CRA_FOXI1 | 65 |
| CRA | CRA_FOXM1 | 38 |
| CRA | CRA_FOXO3 | 47 |
| CRA | CRA_FOXO4 | 76 |
| CRA | CRA_FOXP1 | 30 |
| CRA | CRA_FOXP2 | 22 |
| CRA | CRA_FOXS1 | 24 |
| CRA | CRA_GABARAPL1 | 88 |
| CRA | CRA_GATA1 | 55 |
| CRA | CRA_GATA2 | 25 |
| CRA | CRA_GATA3 | 85 |
| CRA | CRA_GATA5 | 83 |
| CRA | CRA_GSX2 | 53 |
| CRA | CRA_HINFP | 40 |
| CRA | CRA_HMGB1 | 42 |
| CRA | CRA_HMGB2 | 35 |
| CRA | CRA_HNF1B | 43 |
| CRA | CRA_HOXA1 | 94 |
| CRA | CRA_HOXA2 | 80 |
| CRA | CRA_HOXA6 | 96 |
| CRA | CRA_HOXB13 | 153 |
| CRA | CRA_HOXB4 | 20 |
| CRA | CRA_HOXB7 | 102 |
| CRA | CRA_HOXB9 | 77 |
| CRA | CRA_HOXC10 | 31 |
| CRA | CRA_HOXC5 | 70 |
| CRA | CRA_HOXD11 | 30 |
| CRA | CRA_HOXD3 | 100 |
| CRA | CRA_HSF1 | 49 |
| CRA | CRA_HSF2 | 121 |
| CRA | CRA_HSF4 | 37 |
| CRA | CRA_HSPA8 | 56 |
| CRA | CRA_HSPD1 | 26 |
| CRA | CRA_IRF4 | 89 |
| CRA | CRA_IRF5 | 41 |
| CRA | CRA_IRF7 | 44 |
| CRA | CRA_JUN | 29 |
| CRA | CRA_JUND | 118 |
| CRA | CRA_KLF3 | 61 |
| CRA | CRA_KLF4 | 21 |
| CRA | CRA_KLF6 | 37 |
| CRA | CRA_LCOR | 62 |
| CRA | CRA_LHX1 | 62 |
| CRA | CRA_LHX5 | 115 |
| CRA | CRA_MAZ | 52 |
| CRA | CRA_MEF2A | 98 |
| CRA | CRA_MEF2C | 67 |
| CRA | CRA_MITF | 48 |
| CRA | CRA_MRRF | 53 |
| CRA | CRA_MXD3 | 100 |
| CRA | CRA_MYB | 62 |
| CRA | CRA_MYBL2 | 61 |
| CRA | CRA_MYC | 187 |
| CRA | CRA_NFATC2 | 3 |
| CRA | CRA_NFATC4 | 38 |
| CRA | CRA_NFE2L2 | 20 |
| CRA | CRA_NFIA | 21 |
| CRA | CRA_NFIL3 | 32 |
| CRA | CRA_NFKB1 | 49 |
| CRA | CRA_NFKB2 | 120 |
| CRA | CRA_NKX1-1 | 91 |
| CRA | CRA_NKX2-2 | 84 |
| CRA | CRA_NOBOX | 74 |
| CRA | CRA_NPAS1 | 6 |
| CRA | CRA_NPAS3 | 66 |
| CRA | CRA_NR0B1 | 83 |
| CRA | CRA_NR4A1 | 24 |
| CRA | CRA_NR4A2 | 61 |
| CRA | CRA_NRF1 | 51 |
| CRA | CRA_ONECUT1 | 40 |
| CRA | CRA_PARP1 | 72 |
| CRA | CRA_PAX2 | 71 |
| CRA | CRA_PAX4 | 86 |
| CRA | CRA_PAX8 | 90 |
| CRA | CRA_PBX2 | 31 |
| CRA | CRA_PGR | 106 |
| CRA | CRA_PHOX2A | 3 |
| CRA | CRA_PHOX2B | 30 |
| CRA | CRA_PLAU | 53 |
| CRA | CRA_POU2F1 | 72 |
| CRA | CRA_POU3F2 | 61 |
| CRA | CRA_POU5F1 | 61 |
| CRA | CRA_PPARG | 64 |
| CRA | CRA_PPARGC1A | 40 |
| CRA | CRA_PPARGC1B | 142 |
| CRA | CRA_PRKRIR | 68 |
| CRA | CRA_RUNX2 | 23 |
| CRA | CRA_RXRB | 34 |
| CRA | CRA_SALL4 | 71 |
| CRA | CRA_SIM1 | 106 |
| CRA | CRA_SIM2 | 105 |
| CRA | CRA_SIRT1 | 52 |
| CRA | CRA_SOX10 | 82 |
| CRA | CRA_SOX11 | 56 |
| CRA | CRA_SOX15 | 140 |
| CRA | CRA_SOX17 | 18 |
| CRA | CRA_SOX2 | 18 |
| CRA | CRA_SOX5 | 25 |
| CRA | CRA_SOX7 | 59 |
| CRA | CRA_SOX9 | 66 |
| CRA | CRA_SP1 | 57 |
| CRA | CRA_SREBF1 | 80 |
| CRA | CRA_SREBF2 | 8 |
| CRA | CRA_SRY | 35 |
| CRA | CRA_STAT1 | 49 |
| CRA | CRA_STAT3 | 25 |
| CRA | CRA_STAT4 | 43 |
| CRA | CRA_TBX3 | 30 |
| CRA | CRA_TCF3 | 37 |
| CRA | CRA_TCF7L2 | 12 |
| CRA | CRA_TEAD3 | 119 |
| CRA | CRA_TFAM | 54 |
| CRA | CRA_TFAP2A | 58 |
| CRA | CRA_TFCP2 | 120 |
| CRA | CRA_TFDP1 | 48 |
| CRA | CRA_TFE3 | 82 |
| CRA | CRA_TFEB | 61 |
| CRA | CRA_TFEC | 11 |
| CRA | CRA_THRA | 50 |
| CRA | CRA_TLX2 | 140 |
| CRA | CRA_TWIST1 | 62 |
| CRA | CRA_VAX1 | 79 |
| CRA | CRA_VEZF1 | 93 |
| CRA | CRA_XBP1 | 42 |
| CRA | CRA_XRCC5 | 100 |
| CRA | CRA_ZBTB2 | 66 |
| CRA | CRA_ZBTB7A | 17 |
| CRA | CRA_ZBTB7B | 22 |
| CRA | CRA_ZBTB7C | 35 |
| CRA | CRA_ZFHX3 | 36 |
| CRA | CRA_ZFX | 51 |
| CRA | CRA_ZIC3 | 84 |
| CRA | CRA_ZIC5 | 10 |
| CRA | CRA_ZKSCAN3 | 60 |
| CRA | CRA_ZKSCAN4 | 39 |
| CRA | CRA_ZNF219 | 31 |
| CRA | CRA_ZNF263 | 36 |
| CRA | CRA_ZNF281 | 61 |
| CRA | CRA_ZNF318 | 101 |
| CRA | CRA_ZNF35 | 88 |
| CRA | CRA_ZNF367 | 35 |
| CRA | CRA_ZNF384 | 135 |
| CRA | CRA_ZNF695 | 54 |
| CRI | CRI_ALX1 | 82 |
| CRI | CRI_APOD | 74 |
| CRI | CRI_ARID5A | 72 |
| CRI | CRI_ARNTL2 | 95 |
| CRI | CRI_ARX | 65 |
| CRI | CRI_ATF4 | 73 |
| CRI | CRI_ATF6 | 2 |
| CRI | CRI_ATF7 | 97 |
| CRI | CRI_ATG3 | 98 |
| CRI | CRI_ATG4C | 80 |
| CRI | CRI_ATG4D | 37 |
| CRI | CRI_ATG5 | 65 |
| CRI | CRI_BATF3 | 31 |
| CRI | CRI_BRCA1 | 25 |
| CRI | CRI_BSX | 98 |
| CRI | CRI_CDX1 | 146 |
| CRI | CRI_CDX2 | 46 |
| CRI | CRI_CEBPA | 49 |
| CRI | CRI_CEBPB | 59 |
| CRI | CRI_CEBPZ | 17 |
| CRI | CRI_CPEB1 | 66 |
| CRI | CRI_CREB1 | 68 |
| CRI | CRI_CREB3L1 | 73 |
| CRI | CRI_CTCF | 19 |
| CRI | CRI_CTCFL | 68 |
| CRI | CRI_DLX2 | 27 |
| CRI | CRI_DLX4 | 89 |
| CRI | CRI_DLX5 | 16 |
| CRI | CRI_DLX6 | 87 |
| CRI | CRI_DNMT1 | 58 |
| CRI | CRI_DRGX | 35 |
| CRI | CRI_DUXA | 55 |
| CRI | CRI_E2F1 | 84 |
| CRI | CRI_E2F3 | 74 |
| CRI | CRI_E2F4 | 194 |
| CRI | CRI_E2F5 | 72 |
| CRI | CRI_E2F7 | 99 |
| CRI | CRI_EGR1 | 70 |
| CRI | CRI_EGR2 | 55 |
| CRI | CRI_EN1 | 18 |
| CRI | CRI_EOMES | 155 |
| CRI | CRI_ESR1 | 11 |
| CRI | CRI_ESX1 | 45 |
| CRI | CRI_ETS1 | 31 |
| CRI | CRI_ETS2 | 44 |
| CRI | CRI_ETV5 | 107 |
| CRI | CRI_EZH2 | 86 |
| CRI | CRI_FERD3L | 89 |
| CRI | CRI_FLI1 | 60 |
| CRI | CRI_FOSL1 | 46 |
| CRI | CRI_FOXA2 | 28 |
| CRI | CRI_FOXA3 | 79 |
| CRI | CRI_FOXB2 | 106 |
| CRI | CRI_FOXC1 | 39 |
| CRI | CRI_FOXD2 | 58 |
| CRI | CRI_FOXF2 | 52 |
| CRI | CRI_FOXI1 | 48 |
| CRI | CRI_FOXM1 | 27 |
| CRI | CRI_FOXO3 | 25 |
| CRI | CRI_FOXO4 | 59 |
| CRI | CRI_FOXP1 | 64 |
| CRI | CRI_FOXP2 | 125 |
| CRI | CRI_FOXS1 | 101 |
| CRI | CRI_GABARAPL1 | 58 |
| CRI | CRI_GATA1 | 54 |
| CRI | CRI_GATA2 | 70 |
| CRI | CRI_GATA3 | 56 |
| CRI | CRI_GATA5 | 58 |
| CRI | CRI_GSX2 | 83 |
| CRI | CRI_HINFP | 16 |
| CRI | CRI_HMGB1 | 92 |
| CRI | CRI_HMGB2 | 42 |
| CRI | CRI_HNF1B | 41 |
| CRI | CRI_HOXA1 | 33 |
| CRI | CRI_HOXA13 | 38 |
| CRI | CRI_HOXA2 | 50 |
| CRI | CRI_HOXA6 | 120 |
| CRI | CRI_HOXB13 | 47 |
| CRI | CRI_HOXB4 | 114 |
| CRI | CRI_HOXB7 | 25 |
| CRI | CRI_HOXB9 | 27 |
| CRI | CRI_HOXC10 | 27 |
| CRI | CRI_HOXC5 | 33 |
| CRI | CRI_HOXD11 | 52 |
| CRI | CRI_HOXD3 | 7 |
| CRI | CRI_HSF1 | 89 |
| CRI | CRI_HSF2 | 49 |
| CRI | CRI_HSF4 | 31 |
| CRI | CRI_HSPA8 | 2 |
| CRI | CRI_HSPD1 | 13 |
| CRI | CRI_IRF4 | 46 |
| CRI | CRI_IRF5 | 87 |
| CRI | CRI_IRF7 | 49 |
| CRI | CRI_JUN | 35 |
| CRI | CRI_JUND | 69 |
| CRI | CRI_KLF3 | 75 |
| CRI | CRI_KLF4 | 38 |
| CRI | CRI_KLF6 | 140 |
| CRI | CRI_LCOR | 99 |
| CRI | CRI_LHX1 | 90 |
| CRI | CRI_LHX5 | 94 |
| CRI | CRI_MAZ | 88 |
| CRI | CRI_MEF2A | 35 |
| CRI | CRI_MEF2C | 8 |
| CRI | CRI_MITF | 50 |
| CRI | CRI_MRRF | 75 |
| CRI | CRI_MXD3 | 44 |
| CRI | CRI_MYB | 45 |
| CRI | CRI_MYBL2 | 55 |
| CRI | CRI_MYC | 121 |
| CRI | CRI_NFATC2 | 45 |
| CRI | CRI_NFATC4 | 23 |
| CRI | CRI_NFE2L2 | 91 |
| CRI | CRI_NFIA | 50 |
| CRI | CRI_NFIL3 | 40 |
| CRI | CRI_NFKB1 | 70 |
| CRI | CRI_NFKB2 | 37 |
| CRI | CRI_NKX1-1 | 39 |
| CRI | CRI_NKX2-2 | 46 |
| CRI | CRI_NOBOX | 30 |
| CRI | CRI_NPAS1 | 43 |
| CRI | CRI_NPAS3 | 49 |
| CRI | CRI_NR0B1 | 31 |
| CRI | CRI_NR4A1 | 58 |
| CRI | CRI_NR4A2 | 82 |
| CRI | CRI_NRF1 | 48 |
| CRI | CRI_ONECUT1 | 86 |
| CRI | CRI_PARP1 | 81 |
| CRI | CRI_PAX2 | 66 |
| CRI | CRI_PAX4 | 100 |
| CRI | CRI_PAX8 | 109 |
| CRI | CRI_PBX2 | 73 |
| CRI | CRI_PGR | 89 |
| CRI | CRI_PHOX2A | 58 |
| CRI | CRI_PHOX2B | 60 |
| CRI | CRI_PLAU | 49 |
| CRI | CRI_POU2F1 | 5 |
| CRI | CRI_POU3F2 | 35 |
| CRI | CRI_POU5F1 | 31 |
| CRI | CRI_PPARG | 51 |
| CRI | CRI_PPARGC1A | 15 |
| CRI | CRI_PPARGC1B | 65 |
| CRI | CRI_PRKRIR | 48 |
| CRI | CRI_RUNX2 | 29 |
| CRI | CRI_RXRB | 36 |
| CRI | CRI_SALL4 | 45 |
| CRI | CRI_SIM1 | 87 |
| CRI | CRI_SIM2 | 49 |
| CRI | CRI_SIRT1 | 42 |
| CRI | CRI_SOX10 | 32 |
| CRI | CRI_SOX11 | 94 |
| CRI | CRI_SOX15 | 71 |
| CRI | CRI_SOX17 | 91 |
| CRI | CRI_SOX2 | 125 |
| CRI | CRI_SOX5 | 73 |
| CRI | CRI_SOX7 | 73 |
| CRI | CRI_SOX9 | 113 |
| CRI | CRI_SP1 | 98 |
| CRI | CRI_SREBF1 | 45 |
| CRI | CRI_SREBF2 | 45 |
| CRI | CRI_SRY | 35 |
| CRI | CRI_STAT1 | 80 |
| CRI | CRI_STAT3 | 96 |
| CRI | CRI_STAT4 | 66 |
| CRI | CRI_TBX3 | 78 |
| CRI | CRI_TCF3 | 95 |
| CRI | CRI_TCF7L2 | 13 |
| CRI | CRI_TEAD3 | 33 |
| CRI | CRI_TFAM | 88 |
| CRI | CRI_TFAP2A | 31 |
| CRI | CRI_TFCP2 | 42 |
| CRI | CRI_TFDP1 | 55 |
| CRI | CRI_TFE3 | 37 |
| CRI | CRI_TFEB | 77 |
| CRI | CRI_TFEC | 33 |
| CRI | CRI_THRA | 63 |
| CRI | CRI_TLX2 | 76 |
| CRI | CRI_TWIST1 | 76 |
| CRI | CRI_VAX1 | 110 |
| CRI | CRI_VEZF1 | 97 |
| CRI | CRI_XBP1 | 56 |
| CRI | CRI_XRCC5 | 73 |
| CRI | CRI_ZBTB2 | 110 |
| CRI | CRI_ZBTB7A | 44 |
| CRI | CRI_ZBTB7B | 46 |
| CRI | CRI_ZBTB7C | 65 |
| CRI | CRI_ZFHX3 | 121 |
| CRI | CRI_ZFX | 157 |
| CRI | CRI_ZIC3 | 13 |
| CRI | CRI_ZIC5 | 7 |
| CRI | CRI_ZKSCAN3 | 116 |
| CRI | CRI_ZKSCAN4 | 52 |
| CRI | CRI_ZNF219 | 33 |
| CRI | CRI_ZNF263 | 54 |
| CRI | CRI_ZNF281 | 91 |
| CRI | CRI_ZNF318 | 80 |
| CRI | CRI_ZNF35 | 80 |
| CRI | CRI_ZNF367 | 25 |
| CRI | CRI_ZNF384 | 57 |
| CRI | CRI_ZNF695 | 41 |

### Per-guide counts

Full table: `<repo>\results\seng\inventory_guides.csv`.
Rows: 781.

### Other labels in the perturbation column (not CRA_/CRI_/NT_/WT)

| label | n_cells |
|---|---|
| NA | +2950.000 |

## Gene space

Frozen-ruler overlap genes present in this file: **16781** / 23485.
Missing ruler genes are z=0 after the frozen GTEx mu/sd. prior.count=2.

## Split sizes (frozen before any statistic)

Seed `20260914`. Each of EARLY and LATE uses an independent `Generator(20260914)` (same construction as `split_half`). Half A defines v. Half B is used only in P.

| early_pd | late_pd | n_early | n_early_a | n_early_b | n_late | n_late_a | n_late_b | seed | late_determined |
|---|---|---|---|---|---|---|---|---|---|
| PD14 | PD32 | 581 | 290 | 291 | 4348 | 2174 | 2174 | 20260914 | False |

## Positive control P (half-B only; run before any perturbation vector)

v_B = z(EARLY_B) − z(LATE_B). d_f = z(mix_f) − z(LATE_B). Young fraction f=[0.0, 0.25, 0.5]. Total cells fixed at min(n_EARLY_B, n_LATE_B). Seed `20260914` for the mixture draw. Bootstrap seed `20260918`, B=200, resamples cells inside each mixture; EARLY_B and LATE_B membership stays the full half. The P panel is TMM'd on its own. An INVALID interval does not count as excluding 0.

- P pass=True. monotonic=True. frac(0.50) CI excludes 0 = True [+0.483, +0.643] valid=True. delta(0.50)=-34.053 < 0 = True.
- Noise floor f=0: frac=+0.034 delta=+20.634.
- Delta percentile CIs at f=0, 0.25, and 0.50 do not contain the point estimate. They are **INVALID**. The pass rule uses the point estimate of delta and the frac CI, which does contain its point estimate.

| f | n_late | n_early | n_total | frac | delta | cos | frac_ci_lo | frac_ci_hi | frac_ci_valid | delta_ci_lo | delta_ci_hi | delta_ci_valid | noise_floor |
|---|---|---|---|---|---|---|---|---|---|---|---|---|---|
| +0.000 | 291 | 0 | 291 | +0.034 | +20.634 | +0.050 | -0.033 | +0.160 | True | +32.020 | +51.419 | False | True |
| +0.250 | 218 | 73 | 291 | +0.360 | -9.537 | +0.480 | +0.306 | +0.465 | True | +8.311 | +23.232 | False | False |
| +0.500 | 145 | 146 | 291 | +0.547 | -34.053 | +0.709 | +0.483 | +0.643 | True | -14.724 | +1.363 | False | False |

## Reproduction check

Named hits: CRA_E2F3, CRA_EZH2, CRI_STAT3, CRI_ZFX. A hit counts if cos_g > 0 and N_NT BH q within its modality is <= 0.05.

- 1 / 4 counted. reproduction_pass=False.

Our pipeline does not reproduce their direction. Named hits are not interpreted as toward or away.

## Named hits

delta assumes the TF effect is additive to the construct effect: d_g = z(g) − z(NT same modality) is applied to a LATE cell. Negative delta means that sum is closer to EARLY than LATE was. cos_g alone is not evidence of rejuvenation.

Identity genes in the panel log2-CPM: ['COL1A1', 'COL1A2', 'FN1', 'LUM', 'PDGFRA', 'PDGFRB', 'POSTN', 'PRRX1', 'SERPINH1', 'VIM']. Missing: []. `identity_loss` if the mean drop (NT − g) is > 0.5.

| perturbation | present | n_cells | label | repro_cos | cos | frac | delta | cos_resid | delta_resid | p_cos | q_cos | p_frac | q_frac | p_delta | q_delta | p_cos_resid | q_cos_resid | p_delta_resid | q_delta_resid | cos_ci_lo | cos_ci_hi | cos_ci_valid | frac_ci_lo | frac_ci_hi | frac_ci_valid | delta_ci_lo | delta_ci_hi | delta_ci_valid | cos_resid_ci_lo | cos_resid_ci_hi | cos_resid_ci_valid | delta_resid_ci_lo | delta_resid_ci_hi | delta_resid_ci_valid | prolif_shift | cos_d_u | p_cos_d_u | identity_diff | identity_drop | identity_loss | ruler_score_context |
|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|
| CRA_E2F3 | True | 128 | not_interpreted | False | +0.248 | +0.253 | +21.596 | +0.068 | +35.794 | 0.0050 | 0.0923 | 0.0050 | 0.0831 | 0.0100 | 0.1278 | 0.0050 | 0.0692 | 0.6368 | 1.0000 | +0.119 | +0.268 | True | +0.161 | +0.375 | True | +40.840 | +54.219 | False | +0.020 | +0.083 | True | +55.435 | +68.217 | False | +0.165 | -0.023 | 0.9950 | -0.215 | +0.215 | False | +5.225 |
| CRA_EZH2 | True | 107 | not_interpreted | False | +0.197 | +0.198 | +24.426 | +0.083 | +36.258 | 0.0050 | 0.0923 | 0.0050 | 0.0831 | 0.0149 | 0.1662 | 0.0050 | 0.0692 | 0.0249 | 0.4532 | +0.064 | +0.210 | True | +0.087 | +0.298 | True | +44.544 | +59.962 | False | +0.036 | +0.082 | False | +57.807 | +71.544 | False | +0.148 | -0.034 | 1.0000 | -0.330 | +0.330 | False | +5.552 |
| CRI_STAT3 | True | 96 | not_interpreted | False | +0.155 | +0.157 | +28.007 | -0.002 | +43.241 | 0.0249 | 0.1001 | 0.0249 | 0.0956 | 0.1244 | 0.4391 | 0.5920 | 0.9140 | 0.9403 | 1.0000 | +0.018 | +0.221 | True | +0.024 | +0.328 | True | +48.139 | +67.682 | False | -0.034 | +0.026 | True | +64.137 | +83.166 | False | +0.057 | -0.008 | 0.7960 | -0.148 | +0.148 | False | +5.495 |
| CRI_ZFX | True | 157 | not_interpreted | True | +0.316 | +0.303 | +13.125 | +0.039 | +32.032 | 0.0050 | 0.0344 | 0.0050 | 0.0331 | 0.0149 | 0.1359 | 0.0448 | 0.2671 | 1.0000 | 1.0000 | +0.181 | +0.313 | False | +0.213 | +0.403 | True | +29.310 | +40.616 | False | +0.005 | +0.052 | True | +47.850 | +59.747 | False | +0.129 | -0.054 | 1.0000 | -0.038 | +0.038 | False | +6.939 |

Named-hit intervals with `*_ci_valid=False` are **INVALID** (the 95% percentile interval does not contain the point estimate): 10 column-entries. Across all 340 perturbations in `results/seng/screen.csv`: delta 340/340 INVALID, delta_resid 340/340 INVALID, cos 5/340 INVALID, cos_resid 1/340 INVALID, frac 0/340 INVALID. Gates use the point estimates and the N_NT p/q values, not these intervals.

## Screen-level reading

Modalities are not averaged. Full per-perturbation table: `results/seng/screen.csv`.

cos(v, p_CRA)=+0.567. cos(v, p_CRI)=+0.634.
cos(v, u)=-0.040 N1 p=1.0000 beats N1=False.

Quartile sizes (NT only, cut frozen before P): CRA bottom=167 top=167 n_NT=666, CRI bottom=370 top=370 n_NT=1479.

| modality | n_perturbations | rho_prolif | n_cos_q | n_cos_q_and_positive | n_toward_young | n_toward_among_their_style |
|---|---|---|---|---|---|---|
| CRA | 167 | +0.828 | 0 | 0 | 0 | 0 |
| CRI | 173 | +0.862 | 25 | 25 | 0 | 0 |

Their-style hits in this table are perturbations with cos_g BH q <= 0.05 (the screen bullet). The column n_cos_q_and_positive also requires cos_g > 0, which is the reproduction definition.

## Fired keys

- `screen_tracks_proliferation`
- `passage_axis_not_donor_age`
- `pipeline_disagrees`

key=`pipeline_disagrees`.

`passage_axis_not_donor_age`: aging-in-a-dish does not resemble GTEx donor aging in this gene space.

## Limitations

- Neonatal cells aged by passaging, not donor aging. Passage reversal is not a donor-age effect.
- One cell line.
- Additivity assumption in delta: delta applies the TF-minus-NT vector to a WT LATE cell. That assumes the TF effect is additive to the construct effect (the construct cancels in d_g = z(g) − z(NT)) and that the same additive shift applies on the WT background.
- Pseudobulk per perturbation: each perturbation is one summed count vector, so cell-to-cell heterogeneity inside the guide pool is not the statistic.
- TMM for positive control P uses half-B rows only. The analysis panel (half-A, NT, quartile pseudobulks, qualifying perturbations) is a separate TMM. A joint TMM would put half-A inside P and half-B inside the young direction.
- Null and bootstrap z-scores recompute a TMM factor for the resampled pseudobulk against the analysis-panel reference. v and p stay the panel point estimates.
- Proliferation score is the sum of mean log1p(CP10k) of Seurat cc.genes S genes and of G2M genes present in this file. Quartile cut is the 25th and 75th percentile of NT scores within modality.
- Frozen-ruler score is a context column only.
- This file cannot test Sengstack's in vivo mouse results.
- Southard/Hs27 is out of scope and is not opened.
- GSE325735 is out of scope and is not opened.

## Files

- `src/seng_run.py`
- `src/toward_run.py` and `src/same_run.py` imported, not modified
- `results/seng/`
- `FINDINGS_SENG.md`
- `PROGRESS_SENG.md`
- `data/reference/seurat_cc.genes_tirosh2016.tsv`
