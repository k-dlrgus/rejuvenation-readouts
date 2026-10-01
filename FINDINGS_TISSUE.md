# FINDINGS_TISSUE — low-turnover tissue selection (after blood Stage B)

**Status: T1–T4 complete as measurement.** Blood (OneK1K) remains closed; Stage C/D were not run on OneK1K; prior `FINDINGS*.md` and `FALSIFICATION.md` were not modified. Seed `20260914`. Reproduced by `notebooks/tissue_selection.ipynb` (`src/tissue_*.py`). Run date 2026-09-14.

**The number that decides the tissue.** Held-out within-batch age R² from cell-type composition vs from within-cell-type expression, same metric as Stage B (blood: composition **0.493**, expression **0.482**).

| tissue (atlas measured) | composition R² | expression R² | expr − comp | expression beats composition? |
|---|---:|---:|---:|---|
| **blood PBMC (OneK1K, closed)** | **0.493** | **0.482** | −0.011 | no |
| **brain white matter** (`c05e6940-…`, 20 donors) | **−0.290** | **+0.194** | **+0.484** | **yes** |
| skeletal muscle (HCA-harmonized, 27 donors) | +0.032 | +0.009 | −0.023 | no |
| retina scRNA (20-donor file, 16 with numeric age) | −0.209 | −0.004 | +0.205 | yes on sign, **both null** |

**Plain statement.** Compositional dominance is **not general**. In post-mitotic white-matter glia/neurons, 11 within-type transcriptomes predict held-out within-batch age (R² = 0.19) and 14 CLR proportions predict nothing (R² = −0.29). That is the opposite of blood. Muscle and the downloadable retina scRNA file have essentially no age signal of either kind at this n. The 104-donor retina snRNA atlas and the 291-donor DLPFC Aging_Cohort were **not** run through T3 (37 GB / 12 GB); they remain the scale-up resources.

**Recommendation.** **Primary: brain** (measured on the white-matter atlas; scale-up = DLPFC Aging_Cohort). **Backup: retina** (Human Retina Atlas snRNA, 104 donors — T3 not yet run at that scale; Phase 2 mouse/human RPE OSK data exist). Do not use the HCA-harmonized muscle as primary.

---

## T1 — search (nothing invented)

CELLxGENE Discover curation API, fetched 2026-09-14: 2,226 datasets → 1,621 human (`results/tissue/t1_cxg_all_human_datasets.csv`). First-pass regex put kidney `cortex of kidney` in brain; rescored assignment (majority of tissue labels, organoids and cancer atlases dropped) leaves **59** primary candidates (`t1_cxg_primary_candidates.csv`; log `t1_curate.txt`).

| tissue | full atlases (≥10 donors, numeric span ≥40 y, ≥3 types, not organoid/cancer) | best dedicated atlas (from fetched metadata) |
|---|---|---|
| retina | 6 full (+ 10 cell-type subsets of the same two collections) | `d6505c89-c43d-4c28-8c4f-7351a5fd5528` snRNA all cells, **104 donors**, 3.18 M nuclei, ages 3–90 (50 distinct years), 31 types, 10x 3' v3, nucleus, normal. Collection *Single cell atlas of the human retina*, DOI **10.1038/s41588-025-02454-1**. 37.7 GB. |
| brain | 15 | Aging_Cohort DLPFC `4442d412-…`, **291 donors**, 1.33 M nuclei, 0–89, 25 types, 10x 3' v3, nucleus. Collection *Population-scale cross-disorder atlas of the human prefrontal cortex*, DOI **10.1101/2024.10.31.24316513**. 12.4 GB. Also MSSM 1,042 / HBCC 300 (disease-heavy, 36 / 14 GB) and SEA-AD 83 donors (53 GB). |
| skeletal muscle | 1 | `15d374d6-…` **32 donors**, 202 k cells, 19–90, 19 types, mixed 10x + Smart-seq2, cell+nucleus, normal. Collection *Automatic cell-type harmonization…* (HCA integration), DOI **10.1016/j.cell.2023.11.026**. 1.44 GB. **Integrated, not a dedicated aging study.** |
| heart | 2 full | HCA-harmonized Heart `364bd0c7-…` 50 donors, 20–75, 6.4 GB (integrated); myocarditis+controls `fe7aae33-…` 42 donors, 19–70, 5.7 GB (disease-mixed). |
| kidney | 7 full | KPMP v2.0 `7ff0197b-…` 223 donors, 14–76, 1.39 M nuclei, AKI/CKD/normal, 13.8 GB. |
| liver | 7 full | HLiCA `ef3055e1-…` 110 donors, 21–80, normal, integrated, 5.7 GB; Guilliams macrophage niche `e84f2780-…` 19 donors, 28–77, 1.0 GB. High hepatocyte turnover — last resort. |

GEO E-utilities (human aging + sc/sn, then mouse, then OSK/OSKM/partial reprogramming): 297 records, 228 unique accessions (`t1_geo_hits.csv`, `t1_geo_search.txt`). Dedicated human aging sc/sn hits that are actually the tissue (not title collisions): retina AMD atlases GSE137537/GSE137846/GSE137847; brain GSE291605 (aging clocks), GSE243292 (DLPFC); muscle GSE196554 (12 muscle-stem-cell samples), GSE268953 (snMultiome sarcopenia); no clean human heart/kidney aging sc atlas beyond what CXG already has. Every accession above was returned by NCBI; none was assumed.

### Batch recoverability (from metadata, before download)

A dataset without a **shared-donor** technical factor cannot be corrected the way OneK1K was. Library/sample IDs that alias donor (one library per donor) make within-batch age R² undefined.

| atlas | batch recoverable from obs? | notes |
|---|---|---|
| retina snRNA 104-donor (not downloaded) | **unknown until the 37 GB file is opened** | CXG schema has `donor_id`, `assay`, `tissue` (fovea/macula/periphery). Whether a sequencing-pool column exists was not inspected. Flag: not audited. |
| retina scRNA 20-donor (downloaded) | **yes, weak** | `tissue` (3 regions), `assay` (v2/v3), `study_name` (4 studies), `sample_collection_year`. `library_id` aliases donor. 149,852 / 265,767 cells labelled only `adult stage` (no year) — dropped. |
| brain white matter (downloaded) | **yes — OneK1K-like** | `SequencingPool` (7 pools, 2–4 donors each), `10XBatch` (9 levels). |
| muscle HCA (downloaded) | **yes — chemistry** | `assay` (3), `Dataset` (2 source studies), `suspension_type` (2). No sequencing-pool column. |
| DLPFC Aging_Cohort / SEA-AD / kidney / heart | not opened | |

### Reprogramming (Phase 2 coverage)

CELLxGENE human listing: **0** datasets whose title/collection matches OSK/OSKM/Yamanaka/partial reprogramming. GEO SOFT was fetched for every reprogramming-query accession (`t1_geo_reprogramming_meta.csv`, `t1_reprogramming_extra.txt`). GSE145980 was fetched on a guessed Lu-2020 ID and is **zebrafish heart regeneration**, not OSK vision — recorded so it is not reused.

| tissue | aging atlas (human sc/sn) | reprogramming data (fetched) | both? |
|---|---|---|---|
| retina | yes (CXG 104-donor snRNA) | **GSE304044** SuperSeries, *Decoupling oxidative resilience from reprogramming to rejuvenate the RPE and restore vision* (human + mouse; bulk/ATAC subseries). **GSE307031** mouse retina Xenium aging + Müller-glia test. No human retina OSK **scRNA** time course in the returned set. | **yes** (human atlas + mouse/human RPE OSK-adjacent) |
| brain | yes (DLPFC 291; WM 20) | **GSE276656** mouse mPFC engram OSK multiome (PMID 41672073). **GSE224438** mouse SVZ partial-reprogramming scRNA (PMID 38553564). **GSE271794** mouse neocortex transient Yamanaka in development. No human brain OSK scRNA. | **yes** (human atlas + mouse OSK sc) |
| skeletal muscle | yes (HCA 32, integrated) | no tissue OSK scRNA. Closest: GSE176206 murine multi-type partial reprogramming (PMID 35690067). | **no** (atlas only) |
| heart | yes (HCA 50 / myocarditis 42) | no OSK. GSE270268 is PHF7 cardiac reprogramming. | **no** |
| kidney | yes (KPMP 223) | none in the reprogramming queries | **no** |
| liver | yes (HLiCA 110) | mouse in vivo GSE144600, GSE274988, GSE201710 | yes, but high turnover |

A tissue with an aging atlas and no reprogramming data cannot complete Phase 2. That eliminates muscle, heart, and kidney as Phase-2 primaries.

---

## T2 — batch audit (top 3 downloaded)

Trio chosen for downloadability among full, non-cancer, tissue-majority atlases: **brain white matter** (0.43 GB), **skeletal muscle** (1.44 GB), **retina scRNA** (2.54 GB). The scientifically larger retina snRNA (37 GB) and DLPFC Aging_Cohort (12 GB) were listed, not downloaded.

| | brain WM | muscle HCA | retina scRNA | blood (ref.) |
|---|---:|---:|---:|---:|
| donors after adult + numeric-age filters | 20 | 27 (32 in CXG; 5 lost to non-primary / non-numeric age) | 16 (20 in CXG; 4 lost to `adult stage` only) | 981 |
| age span | 34–74 | 19–90 | 42–89 | 19–97 |
| batch column used | **SequencingPool** (7 pools) | **assay** (3 chemistries) | **tissue** (3 regions) | pool (75) |
| R²(age ~ batch) | **0.272** | 0.067 | 0.190 | 0.23 |
| R²(log depth ~ batch) | **0.626** | **0.642** | 0.013 | 0.86 |
| corr(age, depth) between / within batch | **+0.625 / +0.006** | +0.471 / −0.144 | +0.999* / +0.081 | +0.43 / −0.03 |

\*three tissue means — the +0.999 is not interpretable.

**Reading T2.** White matter has the **same confound as OneK1K**: donors are not age-randomised across sequencing pools (R² = 0.27), depth is a pool property (R² = 0.63), and age correlates with depth between pools only. It is still usable — it needs the corrected (within-batch) model from the start, which is what T3 uses. Muscle’s confound is chemistry, not pool, and is weaker on age. Retina scRNA’s “batch” is anatomical region, not a technical pool; library IDs alias donor. A tissue where batch cannot be identified at all was not in this trio.

Filters (logged, seed 20260914):

- brain WM: 45,528 cells, all primary, all numeric adult age; 15 types → 11 with ≥8 donors; 193 pseudobulks. `results/tables/tissue_brain_wm_filter_log.csv`.
- muscle: 201,575 → 195,843 numeric age → 165,109 primary; 16 types → 13 with ≥9 donors; 214 pseudobulks. `results/tables/tissue_muscle_filter_log.csv`.
- retina scRNA: `is_primary_data` would have kept 2.3 % (cells duplicated into subset files of the same collection) — **not applied**; 265,767 → 115,915 with a year-old `development_stage` (149,852 are `adult stage` only); 17 types → 8 with ≥8 donors; 105 pseudobulks.

---

## T3 — the turnover test

Batch-grouped kernel ridge, within-batch R², 10 repeats when ≥5 batches else leave-one-batch-out. Same code path as Stage B (`pool_center`, `kernel_ridge_oof`, pool-cluster bootstrap). Composition = CLR of all cells (not the ≥20-cell pseudobulk groups).

| atlas | n | batches | composition within-batch R² | expression within-batch R² | naive composition / expression | blood-comparable? |
|---|---:|---:|---:|---:|---|---|
| brain WM | 20 | 7 | **−0.290** {−0.69, +0.07} | **+0.194** [0.11, 0.26] {0.11, 0.22} | −0.55 / +0.04 | **yes** |
| muscle | 27 | 3 | +0.032 | +0.009 | −0.16 / −0.14 | yes (LOBO, 3 identical repeats) |
| retina scRNA | 16 | 3 | −0.209 | −0.004 | −0.30 / −0.25 | weak (batch = region) |

White-matter expression is stable across all 10 repeats (0.11–0.22). Composition is negative in 9/10. The two predictors are not interchangeable. **Cells in this tissue are aging in place; the mix is not an age clock.**

Muscle: neither source predicts age (R² ≤ 0.03). An integrated 27-donor, 3-assay atlas is underpowered and/or over-harmonized for this test.

Retina scRNA: both null. Half the cells lack a numeric age; the batch is anatomy. This is **not** evidence that retina is compositional. It is evidence that this file cannot decide. The 104-donor snRNA atlas is the one that can.

Figures: `results/tissue/figures/t3_*_comp_vs_expr.png`, `t4_comp_vs_expr.png`.

---

## T4 — recommendation

| tissue | donors (measured / best available) | age span | batch auditability | composition vs expression (measured) | reprogramming |
|---|---|---|---|---|---|
| **brain** | 20 WM measured; 291 DLPFC available | 34–74 (WM); 0–89 (DLPFC) | **yes** (SequencingPool; same confound as blood) | **expr 0.194 vs comp −0.290** | mouse OSK sc/multiome (GSE276656, GSE224438) |
| **retina** | 16 scRNA measured (null); **104 snRNA not measured** | 42–89 / 3–90 | scRNA: region only; snRNA unaudited | scRNA both null | human+mouse RPE GSE304044; no human OSK scRNA |
| skeletal muscle | 27 measured | 19–90 | assay / Dataset | 0.009 vs 0.032 (both null) | none in tissue |
| heart | not measured (50 / 42) | 20–75 | unknown | — | no OSK |
| kidney | not measured (223) | 14–76 | unknown | — | none |
| liver | not measured (110) | 21–80 | unknown | — | mouse in vivo; high turnover |
| blood | 981 (closed) | 19–97 | 75 pools | 0.482 vs 0.493 | n/a |

**Primary: brain.** The only completed T3 in which expression beats composition, on the same within-batch metric as blood, with a recoverable technical batch and the same age–depth confound OneK1K already taught us to correct. Phase 2 has mouse OSK single-cell in cortex (mPFC) and SVZ. Scale-up inside the same organ: Aging_Cohort DLPFC (`4442d412-…`, 291 donors, 12 GB) — same collection family as the prefrontal-cortex atlas, not yet downloaded.

**Backup: retina.** Priority-1 tissue, post-mitotic neurons that persist for decades, a 104-donor snRNA atlas already on CXG, and a fetched RPE reprogramming SuperSeries (GSE304044). The 20-donor scRNA T3 did **not** fail because retina is compositional; it failed because numeric age is missing for 56 % of cells and there is no sequencing-pool batch. Do not treat that null as a tissue-level result. Next measurement on retina should be the 104-donor snRNA file, after an obs peek for a shared-donor batch column.

**Not recommended as primary.** Muscle (no age signal, integrated, no OSK). Heart/kidney (no OSK scRNA). Liver (turnover). Blood (closed).

**Is compositional dominance general?** No — not on the one low-turnover atlas that had a usable shared-donor batch and a detectable expression age signal. That is a larger finding than “pick a new tissue”: the blood result is consistent with high cell turnover (the mix shifts; individual cells do not get old). White matter is the first counterexample. The project does not need to change again unless the 104-donor retina snRNA T3 comes back compositional.

---

## Files

| path | content |
|---|---|
| `src/tissue_search_cxg.py`, `tissue_curate.py`, `tissue_search_geo.py`, `tissue_search_reprogramming.py`, `tissue_search_reprog_extra.py` | T1 |
| `src/tissue_peek_obs.py`, `tissue_download_top.py`, `tissue_analyze.py`, `tissue_t4.py` | T2–T4 |
| `notebooks/tissue_selection.ipynb` | runnable end to end |
| `results/tissue/t1_*.csv`, `t1_*.txt` | search logs and tables |
| `results/tissue/t2_*.csv`, `t3_*.csv`, `t4_comparison.csv` | audit and turnover numbers |
| `results/tissue/figures/` | age histograms and composition-vs-expression bars |
| `results/tissue/*_analyze.txt` | full filter logs + T2/T3 text |
| `data/raw/tissue/{brain_wm,muscle,retina_sc}/` | cached h5ads |
| `data/processed/tissue_*_pseudobulk.h5ad` | donor × cell-type sums |

Accessions and URLs in this file were read from the CELLxGENE curation API or from GEO SOFT. None were assumed.
