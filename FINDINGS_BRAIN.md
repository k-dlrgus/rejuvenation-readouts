# FINDINGS_BRAIN — DLPFC Aging_Cohort at scale (after tissue selection)

**Status: R0–R4 complete as measurement.** Blood (OneK1K) remains closed; Stage C/D were not run; prior `FINDINGS*.md` and `FALSIFICATION.md` were not modified. Seed `20260914`. Adult cutoff **≥ 20**. Reproduced by `notebooks/brain_replication.ipynb` (`src/brain_*.py`). Run date 2026-09-14.

**The number that decides the tissue.** Held-out within-batch age R² from cell-type composition vs from within-cell-type expression, same metric as Stage B / T3.

| tissue (atlas measured) | n donors | batches | composition R² | expression R² | expr − comp | expression > composition in repeats? |
|---|---:|---:|---:|---:|---:|---|
| **blood PBMC (OneK1K, closed)** | 981 | 75 pools | **0.493** | **0.482** | −0.011 | no (means) |
| **brain white matter** (`c05e6940-…`, 20 donors) | 20 | 7 SequencingPool | **−0.290** | **+0.194** | **+0.484** | **10/10** |
| **brain DLPFC Aging_Cohort (adult controls)** | **233** | **2 Source (H/M)** | **+0.053** | **+0.444** | **+0.391** | **3/3** |

**Plain statement.** The white-matter flip **replicates at scale** in DLPFC: within-cell-type expression predicts held-out within-site age (R² = 0.44) and CLR composition does not (R² = 0.05). That is the opposite of blood. The signal is **broad** (18/20 cell types positive; median per-type R² = 0.33), not a single glial type carrying everything.

**R4 verdict: (a).** Brain is the primary tissue. Recommend proceeding to Phase 1 (age/identity gene derivation) on this adult-control DLPFC cohort (233 donors, ages 20–89, 2 brain-bank sites).

This is single-**nucleus** data (10x 3′ v3). Metrics are comparable to blood scRNA; absolute sensitivity is not.

---

## R0 — batch peek (obs only; 12 GB matrix not yet downloaded)

Accession verified against the live CELLxGENE curation API, 2026-09-14, collection `84ce6837-548d-4a1f-919f-0bc0d9a3952f` (*Population-scale cross-disorder atlas of the human prefrontal cortex at single-cell resolution*, DOI **10.1101/2024.10.31.24316513**). Dataset `4442d412-91cb-4261-acca-8adf5fa04c11`, title **Aging_Cohort**. Live H5AD URL matches the fetched candidate table: `https://datasets.cellxgene.cziscience.com/13e8f1dd-962f-47b7-9cf7-c71d2b21b8a5.h5ad` (12,421,671,635 bytes). 1,332,155 nuclei × 34,176 genes, 291 donors, 25 cell types, 10x 3′ v3, nucleus, DLPFC. Nothing invented.

**33 obs columns. No sequencing pool, library, run, 10x lane, or hash-pool column.** The methods describe 6-sample hashing pools; they were not uploaded. Candidate shared-donor factors that do **not** alias donor 1:1:

| column | n levels | donors / level | what it is |
|---|---:|---|---|
| **Source** | 2 | H=178, M=113 | brain bank (HBCC / MSSM) |
| genetic_ancestry | 5 | 2–159 | biology, not a batch |
| sex | 2 | 91 / 200 | biology |
| disease | 6 | mostly `normal` | systemic comorbidity |
| assay, tissue, suspension_type | 1 | 291 | constant; unusable |

STOP R0 does **not** fire: Source is a shared-donor factor, so within-batch R² is defined. It is **not** OneK1K-like (75 multiplexed 10x pools) or white-matter-like (7 SequencingPools). Within-batch on this file means **within brain bank**, leave-one-site-out. Psychiatric/neurodegeneration flags (`Schizophrenia`, `Bipolar_Disorder`, `Parkinson_disease`, `DLBD_status`, `FTD_status`, `Vascular_status`) are constant `No` — this file is already the neurotypical Aging_Cohort subset of the cross-disorder atlas. `is_primary_data` is False for every nucleus (collection duplication into the parent PsychAD file); not applied.

Log: `results/brain/r0_peek.txt`, `r0_brain_aging_batch_columns.csv`.

---

## R1 — cohort construction

**Neurotypical controls.** 0 nuclei dropped for brain disease: every psychiatric/neurodegeneration flag is constant `No`. CXG `disease` is systemic comorbidity (atherosclerosis, T1D, T2D) on 32/291 donors; those donors are the paper’s neurotypical controls and were **kept**. Developmental samples are not aging.

**Adult cutoff ≥ 20.** The source paper treats 12–19 as adolescence (developmental group) and 20–39 as young adulthood; the DLPFC transcriptome “stabilized after age 20.” This is the cutoff, not a value tuned for a flip.

| filter | unit | before → after (removed) |
|---|---|---|
| is_primary_data==True | cells | 1,332,155 → 1,332,155 (0) — **not applied** (all False) |
| brain-disease flag | cells | 1,332,155 → 1,332,155 (0) — already neurotypical |
| numeric donor age | cells | 1,332,155 → 1,330,258 (1,897; 1 donor `unknown`) |
| age ≥ 20 | cells | 1,330,258 → **1,054,845** (275,413 developmental) |
| adult control donors | donors | 291 → **233** (57 aged 0–19; 1 non-numeric) |

**233 adult control donors, ages 20–89 (span 69 y), 1,054,845 nuclei.** STOP R1 does not fire (≥40 donors, span ≥40 y). Source split: H (HBCC) 120 donors, median age 44.5, 20–85; M (MSSM) 113 donors, median age 64, 20–89. MSSM is older — the same kind of age–batch non-randomisation OneK1K taught us to correct.

Pseudobulk (min 20 nuclei/group, min 10 donors/type): 3,323 groups, 20 cell types (B/T/NK/plasma dropped as too rare). Composition CLR: 18 types + `other` (19 categories) from all 1,054,845 nuclei of analysis donors. Mean UMI/nucleus = 13,670.

---

## R2 — batch audit

Batch column used: **Source** (2 brain banks). Not a sequencing pool.

| | DLPFC adult | brain WM | blood (ref.) |
|---|---:|---:|---:|
| donors | 233 | 20 | 981 |
| batch | Source (2) | SequencingPool (7) | pool (75) |
| R²(age ~ batch) | **0.264** | 0.272 | 0.23 |
| R²(log depth ~ batch) | **0.049** | 0.626 | 0.86 |
| corr(age, depth) between / within | n/a (2 batches) / **−0.030** | +0.625 / +0.006 | +0.43 / −0.03 |

**Reading R2.** Age is as confounded with site as it was with pool in blood and white matter (R² ≈ 0.26). Depth is **not**: site explains 5 % of log UMI/nucleus vs 86 % of log UMI/cell in OneK1K. Within-site age–depth correlation is −0.03, the same null as blood’s within-pool value. The blood A-gene collapse was a pool-linked depth artefact; that mechanism is weak here. The residual risk is the **unmeasured 6-plex hashing pools** (not in obs), not Source.

---

## R3 — the turnover test

Batch-grouped kernel ridge, within-batch R², same helpers as `src/tissue_analyze.py` (`pool_center`, `kernel_ridge_oof`, CLR, pool-cluster bootstrap). Outer CV is leave-one-site-out (2 sites). Inner alpha uses donor-grouped CV inside the training site because batch-grouped inner CV is degenerate with 1 training site (the literal T3 path would pick α = 0.01 by accident). Outer split remains batch-grouped; no donor and no site in both folds (asserted). Repeats = 3 (same rule as T3 when n_batches < 5); they are **not independent** — the outer split is determined by the 2 sites. A 2-cluster bootstrap CI is nearly degenerate and is not the evidence.

| | composition (CLR, 19) | expression (20 types jointly) |
|---|---:|---:|
| within-batch R² (mean) | **+0.053** {0.049, 0.062} | **+0.444** {0.444, 0.444} |
| literal T3 path (α stuck at 0.01) | −0.004 | +0.444 |
| naive (between-site) R² | −0.678 | −0.487 |
| expression > composition | **3/3 repeats** | |

Negative naive R² means the between-site age structure does not transfer. The within-site comparison is the claim: **expr − comp = +0.391**. Composition R² = +0.05 is a weak beat of the mean, not “composition carries no age information.” Do not ratio the two (composition is near zero).

**Per-cell-type expression within-batch R²** (not a T-cell / single-glia phenomenon):

| cell type | n donors | R² |
|---|---:|---:|
| L2/3 IT glutamatergic | 220 | **0.710** |
| L2/3-6 IT glutamatergic | 215 | **0.689** |
| oligodendrocyte precursor | 232 | 0.518 |
| sst GABAergic | 206 | 0.511 |
| astrocyte | 229 | 0.487 |
| VIP GABAergic | 217 | 0.434 |
| L6 IT glutamatergic | 150 | 0.410 |
| pvalb GABAergic | 209 | 0.397 |
| oligodendrocyte | 232 | 0.356 |
| microglial cell | 219 | 0.296 |
| endothelial / pericyte / VLMC | 160 / 136 / 81 | ~0.05 |
| smooth muscle (n=16) / PVM (n=31) | | −0.03 / −0.09 |

18/20 types positive; median 0.33. Neurons and the major glia all carry age. Vascular types are near null. Blood’s age clock was a T-cell / composition story; this is not the brain analogue of that.

Figures: `results/brain/figures/r3_comp_vs_expr.png`, `r3_per_type.png`, `r1_age_hist.png`, `r2_age_by_source.png`.

---

## R4 — verdict

**(a) The flip replicates at scale → brain is the primary tissue; proceed to Phase 1 (age/identity gene derivation) on this cohort.**

- **233** adult neurotypical control donors, ages **20–89**, **2** batches (Source = HBCC/MSSM), **20** cell types, **3,323** pseudobulks, **25,526** genes after the log-CPM filter.
- Expression beats composition by **+0.39** R², in **3/3** repeats, and in **18/20** cell types. The 20-donor white-matter result was not an n=20 fluke.
- Retina 104-donor snRNA is no longer required to decide the tissue. It remains the backup if Phase 1 on DLPFC fails for a reason this file cannot see (unmeasured hash pools).

**Will Phase 1 A-genes survive within-batch correction?** In blood, C(pool) + depth reduced 140 pooled A-genes to **2**. That correction was aggressive because 75 pools absorbed age and depth was a pool property (R² = 0.86). Here:

- C(Source) is **1 df**. It removes the between-bank age mean (MSSM older) and leaves a 20–85 / 20–89 span inside each site (~120 and ~113 donors). Age is not swallowed.
- Depth is not a site property (R² = 0.049). The blood mechanism is weak.
- The practical risk is the **opposite** of blood: hashing pools of 6 (described in the methods, absent from obs) could still inflate A-genes if they are age-structured. Flag it now. Phase 1 should treat Source as the recoverable batch and state that pool-level correction is impossible on this upload.

Cohort size is large enough for the A/I derivation that blood could not support once pool was modelled. Do not treat 0.444 as a promise that hundreds of clean A-genes exist; it is a promise that within-type transcriptomes carry age after the correction this file allows.

---

## Files

| path | content |
|---|---|
| `src/brain_r0.py`, `brain_analyze.py`, `brain_common.py`, `brain_make_notebook.py` | R0–R4 |
| `notebooks/brain_replication.ipynb` | runnable end to end |
| `results/brain/r0_*.csv`, `r0_peek.txt` | obs peek (no matrix) |
| `results/brain/r1_cohort.csv`, `r1_donors.csv`, `filter_log.csv` | cohort + every filter |
| `results/brain/r2_batch_audit.csv` | age~batch, depth~batch |
| `results/brain/r3_*.csv`, `r4_comparison.csv` | turnover test |
| `results/brain/figures/` | age histograms, comp vs expr, per-type |
| `data/raw/tissue/brain_aging/` | cached 12.4 GB h5ad |
| `data/processed/tissue_brain_aging_pseudobulk.h5ad` | donor × cell-type sums |

Dataset ID `4442d412-91cb-4261-acca-8adf5fa04c11` and URL `13e8f1dd-962f-47b7-9cf7-c71d2b21b8a5.h5ad` were read from `results/tissue/t1_cxg_primary_candidates.csv` and re-verified against the live curation API. None was assumed.
