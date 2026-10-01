# FINDINGS_TARGET — the target vector (expression space only)

**Status:** V1 done; V2 done; V3 done; V4 done. Seed `20260914`. Dataset ID `4442d412-91cb-4261-acca-8adf5fa04c11` (cached Aging_Cohort h5ad; not invented). Does not modify `FINDINGS_TRAJECTORY.md`, `FINDINGS_GEOMETRY.md`, `FINDINGS_BRAIN_PHASE1.md`, or `FALSIFICATION.md`. No gene selection: every gene in the Phase-1 matrix enters every model with a weight. No perturbation data. No TF Atlas. Expression space only.

Reproduced by `notebooks/target_v1.ipynb`, `target_v2.ipynb`, `target_v3.ipynb`, `target_v4.ipynb`, `target_v5.ipynb`. Null usability line: site-stratified shuffle R² ≤ 0.05.

## Headline

**V1 ridge bootstrap:** median pairwise angle 45.1° (STOP=False; bar 60°). **V2 ridge TARGET** retained 0.994 of the age direction (∠ 6.3°); within-site R²=+0.231 (null -0.112); types≤null 0/20; STOP V2=False. **V3 (b):** pairwise TARGET 86.0° (null 85.8°); PC1 0.173 (null 0.162); consensus R²=+0.100 vs own +0.231 (loss +0.144). Type-specific targets — consensus loses most of the signal. There is no universal rejuvenation direction in cortex; interventions would have to be cell-type-specific. **V4a site transfer:** HBCC→MSSM R²=-0.213 (null -1.349); MSSM→HBCC R²=-0.499 (null -1.587). **V4b age-range:** young vs old TARGET median angle 86.1°.

## V1 — recover the axes and bootstrap stability

Cohort: 233 donors, 20 cell types, 25526 genes. Identity subspace: 19 axes from type-centroid SVD. Working space: global z-score; ridge coefficients live in that space. Donor bootstrap with replacement (ridge n_boot=20, PLS-1 n_boot=40). Angles are unsigned.

**Ridge median-over-types pairwise bootstrap angle:** 45.1° (to full-data direction 33.7°). **PLS-1:** 38.6°. Types above 60°: 0/20.
Median angle to G1 frozen ridge directions: 23.98° (spaces differ by G1's second per-type z-score).

| celltype | n_boot_ok | median_pairwise_deg | p05_pairwise_deg | p95_pairwise_deg | median_to_full_deg |
|---|---|---|---|---|---|
| GABAergic neuron | 20 | +45.003 | +37.531 | +54.400 | +34.474 |
| L2/3 intratelencephalic projecting glutamatergic neuron | 20 | +43.394 | +37.747 | +49.532 | +33.617 |
| L2/3-6 intratelencephalic projecting glutamatergic neuron | 20 | +41.311 | +35.197 | +47.024 | +31.670 |
| L5/6 near-projecting glutamatergic neuron | 20 | +48.337 | +38.145 | +57.999 | +35.849 |
| L6 corticothalamic-projecting glutamatergic cortical neuron | 20 | +45.144 | +35.918 | +53.337 | +33.001 |
| L6 intratelencephalic projecting glutamatergic neuron | 20 | +45.158 | +37.257 | +52.314 | +33.325 |
| L6b glutamatergic neuron of the primary motor cortex | 20 | +46.334 | +37.246 | +57.255 | +33.783 |
| VIP GABAergic interneuron | 20 | +44.081 | +38.395 | +50.176 | +33.022 |
| astrocyte | 20 | +43.927 | +37.776 | +50.928 | +34.412 |
| endothelial cell | 20 | +49.382 | +40.161 | +56.838 | +35.469 |
| lamp5 GABAergic interneuron | 20 | +46.625 | +39.635 | +53.666 | +34.796 |
| microglial cell | 20 | +43.995 | +38.051 | +50.590 | +32.703 |
| oligodendrocyte | 20 | +42.118 | +35.555 | +47.154 | +31.726 |
| oligodendrocyte precursor cell | 20 | +42.791 | +36.128 | +48.621 | +32.317 |
| pericyte | 20 | +49.726 | +40.627 | +59.806 | +36.860 |
| perivascular macrophage | 20 | +55.455 | +37.447 | +72.947 | +41.433 |
| pvalb GABAergic interneuron | 20 | +44.642 | +37.750 | +51.163 | +32.734 |
| smooth muscle cell | 19 | +52.029 | +25.618 | +80.835 | +40.825 |
| sst GABAergic interneuron | 20 | +42.909 | +36.687 | +48.994 | +31.665 |
| vascular leptomeningeal cell | 20 | +51.931 | +38.288 | +64.338 | +38.616 |

**STOP V1:** False — the ridge age direction is stable enough under donor resampling to define a target.

## V2 — the target vector

TARGET = the age direction with its identity-subspace component removed, then renormalized. Positive TARGET is the older direction; a rejuvenating move would be −TARGET. Held-out scores are 1D OLS of age on the projection, fit on the training fold only.

Full-cohort ridge: median norm retained = 0.994; median ∠(TARGET, age) = 6.32°; median ∠(age, identity) = 83.7°. PLS-1 companion: retained 0.915, ∠(TARGET, age) 23.79°.

| celltype | n | norm_retained | angle_to_age_deg | angle_to_id_deg | proj_r2_id | alpha |
|---|---|---|---|---|---|---|
| GABAergic neuron | 179 | +0.993 | +6.571 | +83.429 | +0.013 | +0.010 |
| L2/3 intratelencephalic projecting glutamatergic neuron | 220 | +0.991 | +7.771 | +82.229 | +0.018 | +0.010 |
| L2/3-6 intratelencephalic projecting glutamatergic neuron | 215 | +0.984 | +10.136 | +79.864 | +0.031 | +0.010 |
| L5/6 near-projecting glutamatergic neuron | 124 | +0.991 | +7.562 | +82.438 | +0.017 | +0.010 |
| L6 corticothalamic-projecting glutamatergic cortical neuron | 153 | +0.992 | +7.102 | +82.898 | +0.015 | +0.010 |
| L6 intratelencephalic projecting glutamatergic neuron | 150 | +0.988 | +8.761 | +81.239 | +0.023 | +0.010 |
| L6b glutamatergic neuron of the primary motor cortex | 109 | +0.990 | +8.208 | +81.792 | +0.020 | +0.010 |
| VIP GABAergic interneuron | 217 | +0.992 | +7.128 | +82.872 | +0.015 | +0.010 |
| astrocyte | 229 | +0.996 | +4.879 | +85.121 | +0.007 | +0.010 |
| endothelial cell | 160 | +0.996 | +4.839 | +85.161 | +0.007 | +0.010 |
| lamp5 GABAergic interneuron | 205 | +0.993 | +6.745 | +83.255 | +0.014 | +0.010 |
| microglial cell | 219 | +0.996 | +4.854 | +85.146 | +0.007 | +0.010 |
| oligodendrocyte | 232 | +0.998 | +4.013 | +85.987 | +0.005 | +0.010 |
| oligodendrocyte precursor cell | 232 | +0.996 | +4.826 | +85.174 | +0.007 | +0.010 |
| pericyte | 136 | +0.998 | +3.966 | +86.034 | +0.005 | +0.010 |
| perivascular macrophage | 31 | +0.997 | +4.156 | +85.844 | +0.005 | +10000.000 |
| pvalb GABAergic interneuron | 209 | +0.994 | +6.079 | +83.921 | +0.011 | +0.010 |
| smooth muscle cell | 16 | +0.997 | +4.470 | +85.530 | +0.006 | +10000.000 |
| sst GABAergic interneuron | 206 | +0.988 | +8.870 | +81.130 | +0.024 | +0.010 |
| vascular leptomeningeal cell | 81 | +0.995 | +6.008 | +83.992 | +0.011 | +0.010 |

| scheme | model | r2 | mae | pooled_r2 | site_strat_r2 | null_mean | null_p95 | p | null_usable | n_perm | types_lose |
|---|---|---|---|---|---|---|---|---|---|---|---|
| within_site | ridge_target | +0.231 | +10.577 | +0.411 | +0.231 | -0.112 | -0.067 | +0.048 | True | 20 | 0/20 |
| within_site | ridge_age | +0.222 | +10.468 | +0.416 | +0.222 | -0.132 | -0.093 | +0.048 | True | 20 | 0/20 |
| within_site | pls1_target | +0.025 | +11.975 | +0.266 | +0.025 | -0.159 | -0.102 | +0.048 | True | 20 | 0/20 |
| loso | ridge_target | -0.077 | +14.970 | -0.077 | -0.410 | -0.844 | -0.606 | +0.048 | True | 20 | 1/20 |
| loso | ridge_age | -0.011 | +14.846 | -0.011 | -0.352 | -0.857 | -0.586 | +0.048 | True | 20 | 1/20 |
| loso | pls1_target | -0.329 | +16.586 | -0.329 | -0.740 | -0.922 | -0.748 | +0.048 | True | 20 | 1/20 |

**Within-site ridge TARGET:** primary R²=+0.231 (MAE 10.58 y; null -0.112, p=0.048). Raw ridge age companion R²=+0.222. Types at or below null: 0/20. G1's within-site ridge R²=+0.317 used a second per-type z-score and the un-normalized ridge predictor; this number is 1D OLS in the global-z space that TARGET and P_id share. TARGET ≈ the age direction (6°), so the gap vs G1 is the metric, not identity projection.
**LOSO ridge TARGET:** primary R²=-0.077 (null -0.844, p=0.048).

**STOP V2:** False — TARGET retains held-out age predictivity in a majority of cell types (consistent with G2c: residual age R² survived identity residualization).

## V3 — one target, or twenty?

**V3a.** Pairwise unsigned angle among per-type ridge TARGETs: median 86.0° (permutation null 85.8°, p=0.725, n_perm=50). Only ~4° from orthogonal. The pairwise angle is **not below the permutation null** — the 20 TARGETs are as orthogonal as chance. (Part C's PLS-1 age directions were 78° / PC1 0.40; ridge TARGET is more type-specific still.)

**V3b.** PC1 of the 20 (ok=20) per-type TARGETs explains 0.173 of their variance (null 0.162, p=0.039). That is a statistically detectable but tiny excess — not a shared axis. Jointly-fit within-type ridge TARGET: median angle to per-type TARGETs 72.0° (norm retained after identity projection 0.997).

**V3c within-site.** Own TARGET R²=+0.231; consensus R²=+0.100 (null -0.332, p=0.048); median per-type loss=+0.144; median R²_cons/R²_own=+0.620.
**V3c LOSO.** Own R²=-0.077; consensus R²=-0.147 (null -1.006, p=0.048); median loss=+0.174.

| celltype | r2_own | r2_consensus | loss | ratio |
|---|---|---|---|---|
| GABAergic neuron | +0.173 | +0.168 | +0.005 | +0.974 |
| L2/3 intratelencephalic projecting glutamatergic neuron | +0.593 | +0.364 | +0.229 | +0.614 |
| L2/3-6 intratelencephalic projecting glutamatergic neuron | +0.582 | +0.373 | +0.209 | +0.641 |
| L5/6 near-projecting glutamatergic neuron | +0.059 | +0.104 | -0.045 | +1.758 |
| L6 corticothalamic-projecting glutamatergic cortical neuron | +0.200 | +0.238 | -0.038 | +1.191 |
| L6 intratelencephalic projecting glutamatergic neuron | +0.271 | +0.169 | +0.101 | +0.625 |
| L6b glutamatergic neuron of the primary motor cortex | +0.043 | +0.105 | -0.062 | +2.435 |
| VIP GABAergic interneuron | +0.294 | +0.168 | +0.126 | +0.571 |
| astrocyte | +0.427 | -0.010 | +0.437 | -0.024 |
| endothelial cell | +0.085 | -0.206 | +0.290 | -2.428 |
| lamp5 GABAergic interneuron | +0.149 | +0.096 | +0.053 | +0.644 |
| microglial cell | +0.277 | -0.195 | +0.472 | -0.702 |
| oligodendrocyte | +0.440 | -0.012 | +0.453 | -0.028 |
| oligodendrocyte precursor cell | +0.339 | -0.095 | +0.434 | -0.281 |
| pericyte | -0.017 | -0.028 | +0.012 | +1.697 |
| perivascular macrophage | -0.064 | -0.002 | -0.062 | +0.024 |
| pvalb GABAergic interneuron | +0.261 | +0.129 | +0.133 | +0.492 |
| smooth muscle cell | -0.107 | -0.273 | +0.165 | +2.539 |
| sst GABAergic interneuron | +0.392 | +0.235 | +0.156 | +0.601 |
| vascular leptomeningeal cell | -0.078 | -0.297 | +0.219 | +3.819 |

Glial types lose the consensus (astrocyte / oligodendrocyte / OPC / microglia own R² med=+0.383 → consensus -0.054). The weak shared component is a neuronal compromise, not a cortex-wide axis.

**V3 verdict: (b)** Type-specific targets — consensus loses most of the signal. There is no universal rejuvenation direction in cortex; interventions would have to be cell-type-specific.

## V4 — how trustworthy is it?

### V4a. Site transfer

TARGET fit on one bank, 1D OLS calibrated on that bank, scored on the other (H=HBCC, M=MSSM). This is the known weak point of the age clock (G1 within-site 0.317 → LOSO 0.101). Median-over-types R² is **negative** in both directions: the vector does not predict chronological age in the other bank. It still beats a severely negative permutation null (p=0.048), which only says it is not pure overfit noise — it does not make it a usable transfer clock. L2/3 IT is the exception (positive R² both ways); glia and most interneurons fail.

| train_site | test_site | median_r2 | null_mean | null_p95 | p | n_perm | n_train_donors | n_test_donors |
|---|---|---|---|---|---|---|---|---|
| HBCC | MSSM | -0.213 | -1.349 | -1.056 | +0.048 | 20 | 120 | 113 |
| MSSM | HBCC | -0.499 | -1.587 | -0.909 | +0.048 | 20 | 113 | 120 |

Types with R²>0 on transfer (all other type×direction cells are in `v4a_per_type.csv`):

| train_site | test_site | celltype | r2 | mae | n_tr | n_te |
|---|---|---|---|---|---|---|
| HBCC | MSSM | L2/3 intratelencephalic projecting glutamatergic neuron | +0.545 | +8.762 | 112 | 108 |
| HBCC | MSSM | L2/3-6 intratelencephalic projecting glutamatergic neuron | +0.312 | +10.951 | 111 | 104 |
| HBCC | MSSM | astrocyte | +0.164 | +12.149 | 117 | 112 |
| HBCC | MSSM | oligodendrocyte | +0.040 | +13.183 | 120 | 112 |
| HBCC | MSSM | oligodendrocyte precursor cell | +0.040 | +13.108 | 119 | 113 |
| MSSM | HBCC | L2/3 intratelencephalic projecting glutamatergic neuron | +0.414 | +10.129 | 108 | 112 |
| MSSM | HBCC | astrocyte | +0.012 | +12.491 | 112 | 117 |

### V4b. Age-range validity

Split: young age<50 (n_donors=102) vs old age≥50 (n_donors=131). Median per-type angle between the two TARGETs: **86.1°**. Consensus-to-consensus angle: 82.1°. Young→old transfer R²=-5.390; old→young R²=-6.967.
This split is **confounded with bank**: young donors are 75 HBCC / 27 MSSM; old donors are 45 HBCC / 86 MSSM. The 86° young–old angle therefore mixes lifespan non-stationarity with the site-transfer failure in V4a. Either way, the direction is not a single line that can be extrapolated across the adult span.
The two half-lifespan directions differ sharply (median > 60°): the target is **not one line across the adult lifespan** and cannot be extrapolated from young to old or from old to young as if it were a single axis.

| celltype | n_young | n_old | angle_deg | young_ok | old_ok |
|---|---|---|---|---|---|
| GABAergic neuron | 83 | 96 | +87.375 | True | True |
| L2/3 intratelencephalic projecting glutamatergic neuron | 95 | 125 | +79.683 | True | True |
| L2/3-6 intratelencephalic projecting glutamatergic neuron | 95 | 120 | +80.352 | True | True |
| L5/6 near-projecting glutamatergic neuron | 59 | 65 | +87.713 | True | True |
| L6 corticothalamic-projecting glutamatergic cortical neuron | 64 | 89 | +86.111 | True | True |
| L6 intratelencephalic projecting glutamatergic neuron | 62 | 88 | +84.959 | True | True |
| L6b glutamatergic neuron of the primary motor cortex | 52 | 57 | +86.035 | True | True |
| VIP GABAergic interneuron | 94 | 123 | +84.480 | True | True |
| astrocyte | 99 | 130 | +82.846 | True | True |
| endothelial cell | 83 | 77 | +88.436 | True | True |
| lamp5 GABAergic interneuron | 91 | 114 | +86.981 | True | True |
| microglial cell | 96 | 123 | +86.275 | True | True |
| oligodendrocyte | 102 | 130 | +86.171 | True | True |
| oligodendrocyte precursor cell | 101 | 131 | +83.442 | True | True |
| pericyte | 70 | 66 | +89.066 | True | True |
| perivascular macrophage | 20 | 11 | +89.054 | True | True |
| pvalb GABAergic interneuron | 92 | 117 | +85.168 | True | True |
| smooth muscle cell | 12 | 4 | NA | True | False |
| sst GABAergic interneuron | 94 | 112 | +83.901 | True | True |
| vascular leptomeningeal cell | 46 | 35 | +89.585 | True | True |

### V4c. Biological sanity (interpretation only — not gene selection)

Technical composition of the top-50 weight lists:
- consensus +: n=50  mito=0  ribo=1  cell-cycle=1
- consensus −: n=50  mito=0  ribo=1  cell-cycle=0
- mean-per-type +: n=50  mito=1  ribo=0  cell-cycle=0
- mean-per-type −: n=50  mito=0  ribo=1  cell-cycle=0

Enrichr terms returned: 950; aging-token terms with p<0.05: 25.

Enrichr, consensus lists (top 10 by p):

| list | library | term | p | adj_p | aging_token |
|---|---|---|---|---|---|
| consensus_pos | Reactome_2022 | Diseases Of DNA Repair R-HSA-9675135 | +0.000 | +0.022 | False |
| consensus_pos | Reactome_2022 | Impaired BRCA2 Binding To PALB2 R-HSA-9709603 | +0.002 | +0.041 | False |
| consensus_pos | Reactome_2022 | Defective HDR Thru Homologous Recombination (HRR) Due To BRCA1 Loss-Of-Function R-HSA-9701192 | +0.002 | +0.041 | False |
| consensus_pos | Reactome_2022 | Resolution Of D-loop Structures Thru Synthesis-Dependent Strand Annealing (SDSA) R-HSA-5693554 | +0.002 | +0.041 | False |
| consensus_pos | Reactome_2022 | Resolution Of D-loop Structures Thru Holliday Junction Intermediates R-HSA-5693568 | +0.003 | +0.043 | False |
| consensus_pos | Reactome_2022 | Resolution Of D-Loop Structures R-HSA-5693537 | +0.003 | +0.043 | False |
| consensus_pos | Reactome_2022 | Presynaptic Phase Of Homologous DNA Pairing And Strand Exchange R-HSA-5693616 | +0.004 | +0.045 | False |
| consensus_pos | KEGG_2021_Human | Homologous recombination | +0.005 | +0.077 | False |
| consensus_pos | Reactome_2022 | Defective Homologous Recombination Repair (HRR) Due To BRCA2 Loss Of Function R-HSA-9701190 | +0.005 | +0.045 | False |
| consensus_pos | GO_Biological_Process_2023 | DNA Recombination (GO:0006310) | +0.005 | +0.237 | False |

Aging-token terms with p<0.05 (top 15; mostly per-type, not consensus):

| list | library | term | p | adj_p | aging_token |
|---|---|---|---|---|---|
| microglial cell_pos | Reactome_2022 | Inflammasomes R-HSA-622312 | +0.001 | +0.095 | True |
| microglial cell_pos | MSigDB_Hallmark_2020 | Inflammatory Response | +0.002 | +0.020 | True |
| microglial cell_pos | GO_Biological_Process_2023 | Positive Regulation Of Inflammatory Response (GO:0050729) | +0.002 | +0.123 | True |
| microglial cell_pos | GO_Biological_Process_2023 | Regulation Of Inflammatory Response (GO:0050727) | +0.003 | +0.128 | True |
| microglial cell_pos | KEGG_2021_Human | Inflammatory bowel disease | +0.012 | +0.215 | True |
| microglial cell_pos | GO_Biological_Process_2023 | Leukocyte Chemotaxis Involved In Inflammatory Response (GO:0002232) | +0.012 | +0.135 | True |
| microglial cell_pos | GO_Biological_Process_2023 | Mitochondrial Protein Catabolic Process (GO:0035694) | +0.012 | +0.135 | True |
| oligodendrocyte_pos | GO_Biological_Process_2023 | Positive Regulation Of ER-associated Ubiquitin-Dependent Protein Catabolic Process (GO:1903071) | +0.015 | +0.177 | True |
| L2/3 intratelencephalic projecting glutamatergic neuron_neg | GO_Biological_Process_2023 | Response To Insulin (GO:0032868) | +0.017 | +0.138 | True |
| mean_per_type_neg | GO_Biological_Process_2023 | Insulin Metabolic Process (GO:1901142) | +0.022 | +0.126 | True |
| consensus_neg | GO_Biological_Process_2023 | Insulin Metabolic Process (GO:1901142) | +0.022 | +0.126 | True |
| mean_per_type_pos | Reactome_2022 | Formation Of Apoptosome R-HSA-111458 | +0.027 | +0.238 | True |
| consensus_pos | Reactome_2022 | Formation Of Apoptosome R-HSA-111458 | +0.027 | +0.121 | True |
| astrocyte_neg | GO_Biological_Process_2023 | Positive Regulation Of Neuroinflammatory Response (GO:0150078) | +0.030 | +0.154 | True |
| microglial cell_pos | Reactome_2022 | Cytokine Signaling In Immune System R-HSA-1280215 | +0.030 | +0.195 | True |

Consensus top + weights (15 of 50): TEKT4P2, ENSG00000223779, HSD17B13, CRX, TMEM62, PRSS55, ENSG00000249588, ENSG00000232855, CLNS1A, ENSG00000197376, LINC01304, APIP, CAPN8, PLA2G10FP, LINC02050
Consensus top − weights (15 of 50): C5orf58, PAX8-AS1, TBC1D3D, PAX8, EPN1, CYP4F29P, ENSG00000267737, ZNF354C, OR2B4P, RNF212, RPL10P9, AGAP12P, FRG1DP, CARS1, TMEM191B

**V4c reading:** The consensus top-weight tails are dominated by sparse/p≫n artifacts (45/100 symbols are LINCs, pseudogenes, olfactory receptors, or Ensembl IDs). That is not a readable aging program. The **microglial** TARGET is the exception: 10 aging-token Enrichr hits (inflammatory response / inflammasome). That is type-specific biology, which is the V3 point. Do not select genes from these tails. A perturbation mapped onto TARGET inherits a correlational linear direction, not a CRISPR menu.

## V5 — verdict

**Is there a usable target vector?** Within one bank, within one cell type, as a correlational linear score: yes (V1 stable, V2 predicts age). **One or many?** **twenty type-specific targets.** It is not a portable rejuvenation axis (V4a negative site transfer; V4b 86° young vs old).

Type-specific targets — consensus loses most of the signal. There is no universal rejuvenation direction in cortex; interventions would have to be cell-type-specific.

**How far can it be trusted?**

- Within-site R²=+0.231 vs LOSO R²=-0.077 (the known site-transfer drop).
- Fit HBCC → test MSSM: R²=-0.213 (null -1.349, p=0.048).
- Fit MSSM → test HBCC: R²=-0.499 (null -1.587, p=0.048).
- Young vs old TARGET median angle 86.1° — not one line across the adult lifespan.

**What would have to be true for a perturbation mapped onto this vector to actually rejuvenate a cell?** A later step inherits all of the following, none of which this analysis tests:

1. **Correlational, not causal.** The direction is fit as a supervised association of expression with chronological age. Moving a cell along −TARGET changes the *score*; it does not by itself mean the cell is younger in any biological sense.
2. **Linear.** TARGET is a single unit vector in log-CPM space. Real interventions are nonlinear, dose-dependent, and constrained to a manifold the linear model never sees.
3. **Slow natural aging, not fast intervention.** Weights are estimated from adult cross-sectional aging (decades), not from a perturbation time course. A two-week OSK pulse need not travel this line.
4. **Identity is the 20-type centroid subspace**, not the rest of what a biologist means by identity (subtype, state, donor, activity). Orthogonal-to-centroids is not orthogonal-to-identity in the full sense.
5. **The right TARGET is the per-type one** if V3 is (b); a consensus move can age one type while missing another. A perturbation applied to mixed cortex would have to be type-resolved, or it is not the vector tested here.
6. **Site and age-range transport.** A vector that fails HBCC↔MSSM or young↔old is a description of this cohort's covariance, not a portable rejuvenation axis. Any mapped perturbation inherits that non-transport.
7. **Held-out prediction is necessary and not sufficient.** V2 only asks whether projection onto TARGET still tracks chronological age. It does not ask whether moving along −TARGET reverses function, pathology, or epigenetic age.
8. **p ≫ n ridge.** Every gene has a weight; many weights are shrinkage artifacts. The top-50 list is interpretation, not a CRISPR menu.

## Limitations

1. Two brain banks only (HBCC, MSSM). LOSO is one df of transfer and is reported, not averaged away.
2. Unmeasured 6-plex hashing pools, PMI, RIN — as in FINDINGS_BRAIN_PHASE1.md.
3. Seed `20260914` (`numpy.random.default_rng`). Folds match Phase 1 / geometry / trajectory.
4. No gene selection. Ridge uses SVD-LOO α on globally z-scored genes (no second per-type z-score, so TARGET and P_id share a metric). PLS-1 is the covariance companion.
5. Donor bootstrap of a p ≫ n direction is a hard test; pairwise unsigned angles near 90° are the random baseline.
6. Smooth-muscle and perivascular-macrophage pseudobulks are small-n; their directions are expected to be noisier and are not dropped.
7. V4c Enrichr is over-representation of the largest-weight tails, not a competitive gene-set test, and is not used for selection.
8. TARGET is aligned with aging (positive = older). An intervention that 'opposes aging' would move against it. Sign is conventional.
9. Identity subspace rank is n_types−1. Residualizing it does not residualize donor, site, or unmeasured batch.
10. The 20–50 vs 50–89 split is confounded with bank (young is HBCC-heavy, old is MSSM-heavy); the 86° young–old angle is not a pure lifespan effect.
11. Within-site TARGET R²=+0.231 is not G1's +0.317: global-z 1D OLS vs per-type-z ridge predictor. TARGET≈age (6°), so identity projection is not the cause of the gap.

## Files

| path | content |
|---|---|
| `src/target_common.py`, `target_v1.py` … `target_v4.py`, `target_findings.py` | code |
| `notebooks/target_v1.ipynb` … `target_v5.ipynb` | runnable from a clean checkout |
| `results/target/` | tables, json, logs, `figures/`, `target_W.npz` |
| `FINDINGS_TARGET.md` | this file |

