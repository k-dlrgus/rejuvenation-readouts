# FINDINGS — Age–Identity Separability Benchmark, Phases 0–1

**Status: HALTED at STOP CONDITION 1b.** Phase 0 completed; Phase 1a and 1b completed; **1c and 1d were not run.**
Run date: 2026-09-14. Seed: `20260914` (`src/config.py`). Every number below is reproduced by `notebooks/phase0_data_acquisition.ipynb`
and `notebooks/phase1_scores.ipynb` from a clean checkout (downloads cached in `data/raw/`).

---

## TL;DR

1. **The data are adequate.** OneK1K (Yazar et al. 2022): 981 healthy donors, ages 19–97 (78 distinct per-year values), 16 usable PBMC
   cell types, 1.17 M cells → 9,623 donor × cell-type pseudobulks. STOP 0 did not fire.
2. **Per gene, chronological age is a very weak axis relative to identity.** Under the specified model `expression ~ cell_type + age`,
   the *largest* unique age variance of any gene is **3.5 %** (99th percentile 0.7 %); cell type reaches **97 %**. The two axes
   differ in per-gene effect size by roughly 30×.
3. **The A-gene set produced by the specified model is a technical artifact, not an age signature — and the artifact is
   batch-linked sequencing depth, not proliferation.** Of 140 primary A-genes, only 3 are cell-cycle. But 131/140 (94 %) lose more than
   75 % of their age variance once the 10x pool is included in the model (median retained: 3.4 %; a real effect retains ≥ ~77 %, confirmed
   by simulation). Depth covariates alone (mean UMI/cell, mean genes/cell) strip ~77 % of it. Mechanism: donors were **not
   age-randomized across pools** (R²(age ~ pool) = 0.23, p = 4 × 10⁻²¹); per-cell depth is a pool property (R² = 0.86); age correlates
   with depth *between* pools (r = +0.43) and not *within* pools (r = −0.03).
4. **→ STOP CONDITION 1b fired**: the majority of A-genes are removed as technical. Per the operating mode, the run halts here. This is
   the finding: *the standard donor-level age model on the field's largest PBMC cohort selects depth-artifact genes as "age genes".*
5. A genuine within-batch age signal exists — CDKN2A (p16), LRRN3, SOX4, TSC22D3, KLF6, BCL2, CD70, SERPINE2, TXNIP — but it is weak
   (≤ 1.3 % of variance per gene; only 10 genes pass the primary threshold within-pool, vs 283 in the pooled model) and much of it sits
   in genes with high cell-type variance (MIXED, e.g. LRRN3, LGALS1, CR1, BCL2), which is itself a warning for separability.

---

## Phase 0 — Data acquisition

### 0.1 Confirmed accession GSE216481 (verified, downloaded, not analysed)
Metadata fetched from GEO (`results/phase0_GSE216481_metadata.txt`, `_subseries.txt`):
*A Transcription Factor Atlas of Directed Differentiation*, SuperSeries of 15 SubSeries, 281 samples, PubMed 36608654, public 2023-01-05.
The ~670k-cell hESC TF-overexpression atlas is SubSeries **GSE217460** (`SHAREseq_210322_TFAtlas`, 48 samples). Downloaded to
`data/raw/GSE217460/` (segmented, resumable; NCBI throttled to ~0.3–4 MB/s):
- `GSE217460_210322_TFAtlas_raw.h5ad.gz` (4.69 GB) — full raw counts.
- `GSE217460_210322_TFAtlas_differentiated_raw.h5ad.gz` (0.13 GB) — verified loads: **28,825 cells × 37,528 genes**, X stored **dense
  float32** (4.33 GB decompressed), obs columns `TF, batch, louvain, n_counts, n_genes, percent_mito, batchTF, score, temp,
  dpt_pseudotime, m3_pseudotime, difflouvain`.
- Full raw file: the HDF5 superblock (read from the gzip stream) gives an end-of-file address of **172.1 GB** — X is dense float32 for
  **1,145,823 cells × 37,528 genes** (the pre-QC cell count; the paper's ~670k is post-filtering). That does not fit on disk (154 GB free),
  so `src/stream_convert_tf_atlas.py` converts it in a single pass directly from the gzip stream into a sparse CSR h5ad
  (`data/processed/GSE217460_TFAtlas_raw_sparse.h5ad`; nnz = 2,673,739,132, density 6.22 %; ~10.5 min at ~270 MB/s). Shape and obs
  columns of the full file: `results/phase0_GSE217460_full_verify.txt` (reproduced in the Phase 0 addendum at the end of this file).

### 0.2 Search for aging atlases (nothing assumed; everything the servers returned is on disk)
- **CELLxGENE Discover** curation API (`src/search_cxg.py`): 2,226 datasets → 1,621 human → **236** meet (≥10 donors, numeric age span
  ≥40 y, ≥3 cell types, scRNA) → **87** blood/immune. Full tables: `results/tables/phase0_cxg_all_human_datasets.csv`, `phase0_cxg_candidates.csv`.
- **GEO** esearch (`src/search_geo.py`): 34 hits for aging + single-cell + PBMC/immune; 6 for OneK1K/sc-eQTL; 0 for (super)centenarian.
  `results/phase0_geo_search.txt`.

Candidates (all values from API metadata):

| dataset | accession(s) | cells | donors | age range (distinct yrs) | cell types | assay | tissue | age continuous? |
|---|---|---|---|---|---|---|---|---|
| **OneK1K** (Yazar 2022, *Science*) | GEO GSE196830 / GSE196735; CXG `3faad104-2ab8-4434-816d-474d8d2641db` | 1,248,980 | **981** | **19–97** (78) | 29 | 10x 3' v2 | blood | **yes, per year** |
| **AIDA Phase 1 v2** (*Cell* 2025) | CXG `c838aec3-03ef-4398-b882-0e3912abfff0` | 1,265,624 | 625 | 19–77 (55) | 32 | 10x 5' v2 | blood | yes |
| Immunobiology of Aging Cohort (*Nature* 2025) | CXG `ff4235bb-…` | 3,758,514 | 234 | 40–89 (51) | 37 | 10x Flex | blood (CMV+/−) | yes |
| Indonesia PBMC (2026 preprint) | CXG `ca7d95ac-…` | 462,034 | 199 | 18–60 (41) | 13 | 10x 5' v2 | blood | yes (span 42 y, borderline) |
| altra / RA-at-risk (*Sci Transl Med*) | CXG `789ad837-…` | 2,029,864 | 89 | 21–78 (50) | 38 | 10x 3' v3 | blood (RA mixed) | yes |
| Immune aging project (*Nat Immunol* 2025) | CXG `1b350d0a-…` | 1,281,499 | 24 | 20–75 (12) | 31 | 10x | multi-tissue | yes |
| Tabula Sapiens | CXG `53d208b0-…` | 1,136,218 | 24 | 22–74 (19) | 180 | 10x + Smart-seq | multi-tissue | yes |
| Luo et al. frailty (GEO GSE157007) | GSE157007 | ~114k | 17 | cord / young / old / frail | — | 10x | PBMC | **no — binned** (would fire STOP 0) |

**Primary: OneK1K.** **Replication (held out; never opened in Phase 1): AIDA Phase 1 v2** — independent cohort, ancestry, chemistry.
OneK1K profile (`results/phase0_onek1k_profile.txt`, `results/figures/phase0_onek1k_donor_age_distribution.png`): X = integer UMI counts;
donors by decade 9 / 43 / 64 / 75 / 136 / 256 / 269 / 118 / 11 (teens→90s); 565 F / 416 M; median 1,246 cells per donor (333–3,511);
75 multiplexed 10x pools; sex vs age Welch p = 0.59 (not confounded).

**STOP CONDITION 0: not fired** (981 donors ≥ 10; 78 age values, not binned; span 78 y ≥ 40; 29 cell types ≥ 3).

---

## Phase 1 — The two scores

### 1.0 Pseudobulk construction (`src/pseudobulk.py`; log `results/tables/onek1k_filter_log.csv`)

| step | unit | before | after | removed |
|---|---|---|---|---|
| cells with numeric donor age | cells | 1,248,980 | 1,248,980 | 0 |
| cells whose Azimuth L2 label is in the analysis map (dropped: CD4/CD8 Proliferating, Doublet, Eryth, HSPC, Platelet) | cells | 1,248,980 | 1,243,859 | 5,121 |
| (donor, cell type) groups with ≥ 20 cells | groups | 19,754 | 9,681 | 10,073 |
| cells in retained groups | cells | 1,243,859 | 1,170,622 | 73,237 |
| cell types present in ≥ 100 donors (dropped: NK_CD56bright 33, Plasmablast 14, cDC 7, pDC 3, dnT 1) | cell types | 21 | 16 | 5 |
| cells in retained cell types | cells | 1,170,622 | **1,168,986** | 1,636 |
| genes with CPM > 1 in ≥ 10 % of pseudobulks | genes | 35,528 | **13,904** | 21,624 |

Final: **9,623 pseudobulks × 13,904 genes**, 981 donors, 16 cell types (donors per type: CD4_TCM 980, CD4_Naive 974, NK 977, CD8_TEM 966,
B_naive 844, CD4_TEM 698, CD8_Naive 671, B_memory 614, Treg 606, B_intermediate 493, CD14_Mono 493, gdT 329, CD8_TCM 293, CD16_Mono 293,
CD4_CTL 280, MAIT 112). `NK Proliferating` was merged into NK (as the CL ontology does) so that proliferating cells remain in the data
and the cell-cycle confound is *measured*, not pre-emptively removed.

### 1a. Variance decomposition (`src/decomp.py`, `src/run_1a_decomp.py`; `results/phase1a_decomp.txt`)

**Model and why.** Per gene, OLS on log₂-CPM pseudobulks: `expression ~ C(cell_type) + age`, exactly as specified. Variance partition by
commonality analysis: `unique_ct = R²(full) − R²(age)`, `unique_age = R²(full) − R²(ct)`, `shared = R²(ct) + R²(age) − R²(full)`,
`resid = 1 − R²(full)`; plus `interaction = R²(ct*age) − R²(full)` and `partial_age = unique_age / (1 − R²(ct))`.
Pseudobulk rather than per-cell because (i) it removes cell-level Poisson noise and cell-level depth as sources of "variance", (ii) the
unit of inference for age is the donor, so donor × type is the correct replicate, (iii) with 981 donors contributing most types the design
is nearly balanced, so age and cell type are close to orthogonal (**shared ∈ [−0.011, 0.015]**) and the partition is unambiguous.

**Result.** Distribution over 13,904 genes (fraction of total variance):

| component | median | 90 % | 99 % | 99.9 % | max |
|---|---|---|---|---|---|
| unique cell type | 0.126 | 0.516 | 0.854 | 0.921 | **0.969** (FTL, CD74, CD37, NKG7 …) |
| unique age | 0.0002 | 0.0018 | 0.0074 | 0.0215 | **0.035** (ANKRD12, TSC22D3) |
| interaction (type-specific age slopes) | 0.0025 | 0.0041 | 0.0060 | 0.0100 | 0.015 |
| residual | 0.873 | 0.948 | 0.978 | 0.985 | 0.991 |

Mean partition per gene: cell type 20.7 %, age 0.07 %, shared 0.05 %, interaction 0.26 %, residual 78.9 %.
Sensitivity (pseudobulks with ≥ 100 cells only, n = 3,486, 6 types): age max 0.076, 99th pct 0.019 — sampling noise in small
pseudobulks deflates all R² but does not change the picture (Spearman between rankings 0.64).

**Gene sets.** *Transparency note:* the first sweep grid written before seeing data (`t_age ∈ {0.01 … 0.12}`) produced **empty**
A-gene sets above 0.03; it was re-centred on the observed distribution *before any downstream model existed*
(`results/phase1a_decomp_initial.txt`, comment in `src/decomp.py`). Primary thresholds, chosen from set sizes only:
A: `unique_age ≥ 0.005 & unique_ct ≤ 0.20`; I: `unique_ct ≥ 0.50 & unique_age ≤ 0.001`; MIXED: `unique_age ≥ 0.005 & unique_ct ≥ 0.50`.

| class | n (primary) |
|---|---|
| **A-genes** | **140** |
| **I-genes** | **1,164** |
| MIXED (excluded) | 41 (LRRN3, TPT1, RPS10, EEF1A1, SELL, VIM, HLA-E, FOXP1, CR1 … — mostly ribosomal/translation) |
| other | 12,559 |

Sweep of set sizes (`results/tables/phase1a_threshold_sweep_set_sizes.csv`): A-genes = 328/700/872/1096 (t_age 0.002; k_ct 0.1/0.2/0.3/0.5),
65/140/170/242 (0.005), 13/30/39/64 (0.01), 6/10/11/19 (0.02). I-genes 448–3,026 across the grid. *Downstream performance across the sweep
was not evaluated because the run halted at 1b.* Figures: `results/figures/phase1a_joint_variance_distribution.png`, `phase1a_mean_variance_partition.png`.

**Batch structure discovered here:** donor age is associated with 10x pool — `age ~ C(pool)` across 981 donors, 75 pools: R² = 0.234,
adj. R² = 0.172, **p = 3.8 × 10⁻²¹**. Donors were not age-randomized across pools.

### 1b. Confound removal (`src/confounds.py`, `src/run_1b_confounds.py`, `run_1b_batch_verify.py`, `run_1b_depth_mediator.py`)

**Are the confounders themselves age-associated?** (Spearman within cell type, `results/tables/phase1b_covariate_age_correlations.csv`)
- Mean UMI per cell vs age: **positive in all 16 cell types**, r = +0.07 … +0.33 (mean genes/cell likewise).
- Cell-cycle UMI fraction vs age: |r| ≤ 0.17; proliferating-cell fraction vs age: |r| ≤ 0.13. Proliferation is not strongly age-linked in PBMCs.
- Sex vs age: p = 0.59 — not confounded.

**Pre-specified flags** (thresholds fixed in code before running): Tirosh S/G2M + proliferation list (144 expressed genes); data-driven
cycling |r| ≥ 0.3 with cell-cycle covariates (69); ribosomal (157); mitochondrial (13); chrX/chrY via Ensembl (462); within-type sex DE
|d| ≥ 0.5 & p < 10⁻⁶ (25); depth |r| ≥ 0.3 with mean UMI/cell or genes/cell, cell-type-residualized (236).

**Primary A-genes (n = 140), sequential survival:**

| filter | flagged | removed (incremental) | surviving |
|---|---|---|---|
| cell cycle — Tirosh/proliferation list | 1 (LBR) | 1 | 139 |
| cell cycle — data-driven | 2 (HNRNPH1, CNBP) | 2 | 137 |
| ribosomal | 0 | 0 | 137 |
| mitochondrial | 2 (MT-ATP8, MT-ND2) | 2 | 135 |
| sex-linked (chrX/Y) | 7 | 7 | 128 |
| sex DE | 0 | 0 | 128 |
| depth-correlated | 27 | 23 | **105** |

Literal criterion: cell-cycle **or** technical = **29/140 = 20.7 %** → by that count alone STOP 1b would *not* fire. Cell cycle accounts
for 3/140: **proliferation is not the confound here.**

**But the strongest A-genes are all depth genes**, and depth is batch. The top-ranked A-genes — SET, RIPOR2, DCK, HNRNPH1, C6orf62,
C1orf56, SRSF6, PPP1CB, RSL24D1 (all *down* with age) and EIF3F, RBIS, DONSON, GABRE (*up*) — have cell-type-residualized correlations
with mean UMI/cell of |r| = 0.35–0.71. So the batch check was run:

| adjustment added to `ct + age` | median fraction of A-gene age variance retained | A-genes retaining < 25 % | < 50 % |
|---|---|---|---|
| 10x pool (74 df) | **0.034** | **131 / 140** | 136 / 140 |
| depth covariates only (2 df: log mean UMI/cell, log mean genes/cell) | 0.228 | 74 / 140 | 118 / 140 |
| pool + depth | 0.039 | 131 / 140 | 136 / 140 |

Verification that this is not an artefact of the adjustment (`results/phase1b_batch_verify.txt`):
- **Simulation** of a genuine additive age effect with random pool effects: retention 0.86–1.30 (expected ≈ 1 − R²(age~pool) = 0.77). Observed 0.034.
- **Depth is a pool property**: R²(log mean UMI/cell ~ pool) = 0.855; R²(log mean genes/cell ~ pool) = 0.812.
- **Age–depth correlation is entirely between pools**: across 75 pool means r = +0.43; within pools r = −0.03.
- **Direct test**: primary A-genes' cell-type-residualized |r(expression, age)| falls from median 0.090 (pooled) to **0.018 within pool**,
  *below* the 99th percentile of a within-pool donor-level permutation null (0.031). Whatever the A-genes track, it is not age within a batch.
- Genome-wide, genes with `unique_age ≥ 0.002` fall from 1,256 (pooled) to 116 (within pool); ≥ 0.005 from 283 to **10**; ≥ 0.01 from 76 to 2.

Figures: `results/figures/phase1b_Agene_confound_filters.png`, `phase1b_batch_adjustment_age_signal.png`, `phase1b_pool_age_depth_confound.png`.

**STOP CONDITION 1b — FIRED** (`src/run_1b_verdict.py`, `results/phase1b_STOP_verdict.txt`). The condition is "the majority of A-genes
are removed as cell-cycle or technical". Batch-linked sequencing depth is a technical confound in exactly the sense the spec lists ("genes
correlated with sequencing depth or n_genes per cell"); the pre-specified |r| ≥ 0.3 flag was simply too blunt to catch the mechanism,
whereas modelling the mechanism shows that the A-gene set has no within-batch age association: **118/140 (84 %) lose > 50 % of their
age variance to the two depth covariates alone; 136/140 (97 %) lose > 50 % to pool.** Union of all removal criteria: 136/140 (97.1 %);
**4 A-genes survive everything: ETV7, PHLDA3, CNN2, DUXAP8.** I did not tune any threshold to reach this: the retention numbers are an
order of magnitude away from any defensible cutoff (real-effect expectation ≥ 0.77; observed median 0.03).

---

## What it means (assessment)

**1. The field's standard age model, applied to its largest single-cell cohort, selects a technical signature.** `expression ~ cell_type + age`
with donor-level age is the default in single-cell aging papers. On OneK1K it ranks depth-sensitive housekeeping/translation genes
(SET, HNRNPH1, SRSF6, EIF3F, PTMA, NACA, EIF1 …) as the top age genes, because (a) pools were assembled in an age-nonrandom way and
(b) per-cell UMI depth is set at the pool level. A "transcriptomic clock" trained on these genes would validate beautifully under
*donor*-held-out CV (donors from the same pools are in train and test) and be wrong. **Donor-held-out is necessary but not sufficient when
donors nest in batches that are confounded with the outcome; batch-held-out validation is required.** This is a general point, not specific
to OneK1K — any multiplexed cohort recruited over time is at risk.

**2. For Phase 2 this artifact is worse than proliferation would have been.** The spec's stated worry was that reprogramming raises
proliferation and proliferation genes look like age genes. In quiescent PBMCs proliferation contributes almost nothing (3/140). But
partial reprogramming *massively* changes RNA content per cell and hence per-cell depth. An age score that is secretly a depth score would
"detect rejuvenation" in any perturbation that changes transcript abundance — a false positive that no proliferation filter would catch.

**3. Age is a weak axis per gene; identity is a strong one.** Even the batch-inflated ranking tops out at 3.5 % of variance; the honest
within-pool maximum is 1.3 % (RPS24), with CDKN2A, LRRN3, SOX4, TSC22D3, KLF6, BCL2 ≈ 0.4–1.2 %. Cell type explains up to 97 %. The
gene-level "A-gene" concept (high age variance, low identity variance) is close to empty at any threshold that would be called "high"
for identity. Multivariate models can still integrate many weak genes, but that is a claim to be tested in 1c, not assumed.

**4. The real age genes lean towards identity.** Among the top within-pool age genes, LRRN3 (unique_ct 0.61), LGALS1 (0.83), CR1 (0.63),
BCL2 (0.62), KLF6 (0.53), GTSCR1 (0.60) are MIXED, not A. Immune aging in PBMC is substantially a change in *which subsets* exist (naive →
memory, CD28⁻ CTL expansion) and in subset-defining genes — i.e., age acts *through* identity programmes. This is exactly the hypothesis
the 1d cross-test is designed to falsify, and the 1a/1b data already point toward it. It should temper expectations for separability.

**5. Positive controls survived.** CDKN2A/p16 and LRRN3 emerge as within-pool age genes with the expected direction, so the pipeline can
see real aging biology when it is there; the problem is signal size and confounding, not a broken pipeline.

---

## What was *not* done, and what a continuation would need (recommendations, not executed)

Not run: 1c (age regression, identity classifier), 1d (shuffles, cross-test), the downstream-performance half of the threshold sweep,
and anything on the AIDA replication set (still unopened). `FALSIFICATION.md` is nevertheless written now as a pre-registration for a future 1d.

A defensible continuation (a *new* run, not a workaround of this one) would:
1. Redefine the age term as **within-batch**: `expression ~ cell_type + pool + age` (or equivalently pool-centre both age and expression),
   and add log depth covariates; re-derive A-genes on within-pool age variance. Expect **~10–100 A-genes**, not 140–1,000.
2. Use **pool-grouped** CV (donors nest in pools; assert no donor *and* no pool overlap).
3. Treat STOP 1c (R² < 0.3 within cell type) as a live risk: with ≤ 1.3 % variance per gene, R² ≈ 0.3 requires coherent integration of many
   weak genes, and published bulk-blood transcriptomic clocks rarely exceed R² ≈ 0.3–0.5.
4. Test replication of the within-pool age genes in AIDA (which will need the same pool ↔ age audit before use).

---

## Reproducibility
- How to run from a clean checkout (Windows PowerShell shown; adjust paths on POSIX):
  `python -m venv .venv; .\.venv\Scripts\pip install -r requirements.txt; .\.venv\Scripts\python -m ipykernel install --user --name age-identity-venv;`
  then `.\.venv\Scripts\python -m nbconvert --to notebook --execute --inplace --ExecutePreprocessor.kernel_name=age-identity-venv notebooks\phase0_data_acquisition.ipynb`
  and the same for `notebooks\phase1_scores.ipynb`. Set `PYTHONUTF8=1` on non-UTF-8 consoles. Phase 0 downloads ≈ 24 GB (OneK1K 4.4, AIDA 14.2,
  TF Atlas 4.8). The dense 172 GB TF Atlas matrix is never written to disk: it is converted from the gzip stream into a sparse CSR h5ad
  (12.9 GB, `data/processed/`, ~11 min single-threaded; peak RAM < 1 GB); everything is cached and skipped on rerun.
- Environment: Python 3.12.7 venv; scanpy 1.12.4, anndata 0.13.3, scikit-learn 1.9.1, scipy 1.18.1, statsmodels 0.15.0, numpy 2.5.3,
  pandas 3.0.5, h5py 3.16.0 (`requirements.txt`).
- Seed `20260914` set in `src/config.py`; all randomness (simulation, permutation null) uses `numpy.random.default_rng(SEED)`.
- Raw inputs are cached in `data/raw/` (gitignored); processed pseudobulk in `data/processed/`. All tables in `results/tables/`,
  console logs in `results/*.txt`, figures in `results/figures/`.
- Nothing from the replication set (AIDA) was read.

---

## Phase 0 addendum — GSE216481 / GSE217460 full-file verification (`results/phase0_GSE217460_full_verify.txt`)

Verified after the streaming conversion (`src/stream_convert_tf_atlas.py`; HDF5 EOF address 172,105,013,517 bytes, X dense block
2,048 → 172,001,784,224; 0 decompression rewinds in the final pass). Not analysed — reserved for Phase 2.

| property | value |
|---|---|
| shape (cells × genes) | **1,145,823 × 37,528** (pre-QC; the paper's ~670k cells is after the authors' filtering) |
| `obs` columns | `n_genes`, `percent_mito`, `n_counts`, `batch`, `TF` |
| distinct `TF` labels | **3,369** (matches the ~3,500 TF-isoform ORF library minus dropouts) |
| `batch` | 2 levels |
| `var` index | gene symbols (`A1BG`, `A1BG-AS1`, `A1CF`, `A2M`, …); 21 per-sub-batch `n_cells-*` columns |
| X | float32 raw UMI counts (1 M-value sample: min 1, max 113, 100 % integer-valued); nnz 2,673,739,132, density 6.22 % |
| example cells | `R1.01,R2.01,R3.01,P1.22-0-0` (TFORF0867-NR1H2; 2,186 UMI, 1,642 genes), `R1.01,R2.01,R3.02,P1.22-0-0` (TFORF1728-ZNF695) |

Cell barcodes are SHARE-seq round-1/2/3 + plate identifiers; the ORF identity per cell is the `TF` column (`TFORFxxxx-SYMBOL`).
