# FINDINGS_RERUN — Phase 1 re-derived under the corrected within-pool model

**Status: HALTED at STOP CONDITION A.** Stage A (A1–A4) completed; **Stages B, C and D were not run** as part of this halted run.
*Addendum (same day): Stage B was subsequently executed as a stand-alone measurement on the protocol owner's explicit instruction, with gene sets and
thresholds untouched — see `FINDINGS_STAGEB.md`. Stages C and D remain not run.*
Run date: 2026-09-14. Seed: `20260914` (`src/rerun_common.RERUN_SEED`; every random draw uses `numpy.random.default_rng(RERUN_SEED)`).
Every number below is reproduced by `notebooks/rerun_stageA_age_axis.ipynb` from a clean checkout (raw downloads cached in `data/raw/`).
Outputs are namespaced (`results/rerun/`, `results/tables/rerun_*.csv`, `results/figures/rerun_*.png`); nothing from the prior run was
overwritten. `FINDINGS.md` and `FALSIFICATION.md` were not modified.

This run is not a continuation of the prior run. The prior A-gene set (140 genes) was discarded. Reused: the cached raw data, the pseudobulk
construction (`src/pseudobulk.py`, `src/run_pseudobulk_onek1k.py`, filter log `results/tables/onek1k_filter_log.csv`, unchanged), the log-CPM
normalisation, the confound gene lists and flag thresholds (`src/gene_lists.py`, `src/confounds.py`), and the **unchanged** thresholds
`PRIMARY_THRESH = {t_age 0.005, k_ct 0.20, t_ct 0.50, k_age 0.001}` and sweep grid from `src/decomp.py`.

---

## TL;DR

1. **The correction works as intended and changes only the age axis.** Under `expression ~ C(cell_type) + C(pool) + log(mean UMI/cell) +
   log(mean genes/cell) + age`, the identity axis is untouched (Spearman(prior, corrected `unique_ct`) = 0.9997) while per-gene age variance
   collapses: max **1.26 %** (SOX4) vs 3.5 % before; 99th percentile **0.19 %** vs 0.74 %. Summed over all genes, 65 % of the prior model's
   per-gene age variance was batch/depth. The two depth covariates add essentially nothing once pool is in the model (Spearman 0.982 between
   pool-only and pool+depth `unique_age`) — pool is the operative correction; depth is a pool property.
2. **The within-pool age signal is real but diffuse and weak.** A within-pool donor-level permutation null (1,000 draws) shows the observed
   counts are far above chance — 9 genes ≥ 0.005 (null mean 0.12), 118 ≥ 0.002 (null 2.5; empirical FDR 2 %), 666 ≥ 0.001 (FDR 4 %) — and
   **1,180 genes have a detectable within-pool age association at BH q < 0.05**, but with a median unique age variance of 0.10 %. Age is
   spread thinly over many genes; no gene carries a large within-cell-type age effect.
3. **Known positives are recovered with the expected direction** — CDKN2A/p16 up (rank 12 of 13,904, +0.0084 log₂-CPM/yr), LRRN3 down (rank 4,
   −0.0119/yr), both q = 0.02 — **but neither satisfies the A-gene definition**: LRRN3 has 60 % identity variance (MIXED), CDKN2A 29 % (above
   `k_ct = 0.20`). The correction is not over-aggressive; the A-gene construct is what fails.
4. **A-genes at the pre-registered thresholds: 2 (PHLDA3, ETV7).** Both survive every confound filter (0 removed). Across the whole sweep the
   largest surviving set at `t_age = 0.005` is 5 (`k_ct = 0.5`); ≥ 10 survivors appear only at `t_age = 0.002` (16 / 34 / 55 / 82).
   **→ STOP CONDITION A FIRED (2 < 10).** Thresholds were not relaxed. The rerun halts here.
5. **The within-pool age genes lean toward identity.** Of the 9 genes with `unique_age ≥ 0.005`, 2 are A-eligible, 6 sit in the gap
   (0.2 < `unique_ct` < 0.5: SOX4, RPS24, TSC22D3, SERPINE2, EIF2S3, LIMD2) and 1 is MIXED (LRRN3). Among the 1,180 q < 0.05 genes, those with
   `unique_ct ≥ 0.5` are enriched 2.5× relative to the genome (25.5 % vs 10.3 %). Every classic bulk-blood age gene tested (CCR7, LEF1, NELL2,
   GZMH, B3GAT1, KLRG1) is an **I-gene** here (identity variance 0.56–0.79, within-pool age variance ≤ 0.08 %): their bulk age association is
   compositional. This is the Stage B hypothesis showing up inside Stage A, but Stage B itself was not run.
6. **Negative / cautionary observations**: 78 % of the 1,180 significant within-pool age genes *decrease* with age (a possible residual
   donor-level RNA-content effect not captured by pool + mean depth; not interpreted here); the pool + depth nuisance block explains a median
   2.7 % (max 71 %) of a gene's variance, two orders of magnitude more than age (median 0.01 %); cell-type-specific age slopes (interaction) are
   only modestly above their 15-df chance level (median 0.24 % vs ~0.16 % expected).

---

## What was reused, what was discarded, what was changed

| item | status |
|---|---|
| raw OneK1K h5ad (`data/raw/onek1k/`) | reused (cached) |
| pseudobulk construction: 9,623 donor × cell-type pseudobulks, 981 donors, 16 types, 75 pools; log₂-CPM on 13,904 genes | reused unchanged (`results/tables/onek1k_filter_log.csv`) |
| confound gene lists / flags (Tirosh S/G2M + proliferation; data-driven cycling; ribo; mito; chrX/Y via Ensembl; sex-DE; depth \|r\| ≥ 0.3) | reused unchanged; recomputed (`results/tables/rerun_stageA_gene_flags.csv`; genome-wide counts identical to the prior run: 144 / 69 / 157 / 13 / 462 / 25 / 236) |
| thresholds `t_age 0.005, k_ct 0.20, t_ct 0.50, k_age 0.001` and sweep grid | reused **unchanged** (not re-centred on the new distribution) |
| the prior 140-gene A set, the prior I set, all prior 1a/1b numbers | **discarded**; recomputed in the same script only for side-by-side comparison |
| the age model | **changed**: `ct + age` → `ct + pool + log depth (2) + age`; depth covariates are centred within cell type so the full-model fit is identical to raw log depth but cell-type variance stays attributed to cell type |
| inference | **added**: within-pool donor-level permutation null (the nominal F test is anti-conservative because a donor's ≤ 16 pseudobulks share its age) |
| AIDA (replication set) | **not opened** (Stage D was not reached) |

---

## Stage A — re-derive the age axis within batch (`src/rerun_stageA.py`; log `results/rerun/stageA_report.txt`)

### A0. Design audit (re-stated for this run)
981 donors in 75 pools (9–14 donors per pool, median 13); every donor in exactly one pool (asserted). `age ~ C(pool)`: R² = 0.234,
p = 3.8 × 10⁻²¹ → **76.6 % of donor age variance is within pool** and available to the corrected model. Within-type-centred log depth vs pool:
R² = 0.855 (UMI/cell), 0.812 (genes/cell). Corrected model: 16 cell types + 74 pool contrasts + 2 depth + 1 age = 93 parameters on 9,623 rows.

### A1. Per-gene variance partition
Nuisance N = pool + within-type-centred log depth. `unique_age = R²(N+ct+age) − R²(N+ct)` (within-pool, depth-adjusted); `unique_ct = R²(N+ct+age)
− R²(N+age)`; `shared`, `nuisance = R²(N)`, `interaction = R²(N+ct+ct:age) − R²(full)`, `resid`.

Distribution over 13,904 genes (fraction of total variance), corrected model:

| component | median | 90 % | 95 % | 99 % | 99.9 % | max |
|---|---|---|---|---|---|---|
| unique cell type | 0.1248 | 0.5062 | 0.6599 | 0.8336 | 0.9103 | 0.9610 |
| **unique age (within pool)** | **0.0001** | **0.0007** | **0.0010** | **0.0019** | **0.0046** | **0.0126** (SOX4) |
| shared | 0.0000 | 0.0011 | 0.0018 | 0.0036 | 0.0059 | 0.0108 |
| nuisance (pool + depth) | 0.0265 | 0.0733 | 0.1033 | 0.2159 | 0.4599 | 0.7087 |
| interaction (type-specific age slopes) | 0.0024 | 0.0039 | 0.0045 | 0.0058 | 0.0099 | 0.0146 |
| residual | 0.8367 | 0.9271 | 0.9454 | 0.9636 | 0.9746 | 0.9805 |
| partial age (share of within-type-and-pool variance) | 0.0002 | 0.0010 | 0.0015 | 0.0041 | 0.0153 | 0.0462 |

Mean partition per gene: cell type 20.4 %, age 0.03 %, shared 0.03 %, nuisance 3.8 %, interaction 0.25 %, residual 75.7 %.

**`unique_age`: prior pooled model vs pool-only vs corrected** (`results/tables/rerun_stageA_unique_age_quantiles_prior_vs_corrected.csv`):

| model | median | 90 % | 95 % | 99 % | 99.9 % | max | n ≥ 0.001 | n ≥ 0.002 | n ≥ 0.005 | n ≥ 0.01 | n ≥ 0.02 |
|---|---|---|---|---|---|---|---|---|---|---|---|
| pooled (prior run: `ct + age`) | 0.00024 | 0.00183 | 0.00291 | 0.00736 | 0.02154 | 0.03500 | 2,730 | 1,256 | 283 | 76 | 19 |
| pool-only (`ct + pool + age`) | 0.00010 | 0.00062 | 0.00092 | 0.00184 | 0.00454 | 0.01262 | 595 | 116 | 10 | 2 | 0 |
| **corrected (`ct + pool + depth + age`)** | 0.00011 | 0.00066 | 0.00098 | 0.00187 | 0.00457 | 0.01262 | 666 | 118 | **9** | 2 | 0 |

Spearman(pooled, corrected `unique_age`) = 0.329 — the pooled ranking is largely uninformative about within-pool age. Spearman(pool-only,
corrected) = 0.982 — depth covariates add almost nothing once pool is modelled. For the 283 genes with pooled `unique_age ≥ 0.005`, the median
fraction retained under the corrected model is **0.044** (a real effect would retain ≈ 0.77; simulation in the prior run: 0.86–1.30).
Spearman(prior, corrected `unique_ct`) = 0.9997, max |Δ| = 0.097 — the identity axis is unaffected.

Top 15 genes by corrected `unique_age` (full table: `results/tables/rerun_stageA_variance_decomposition.csv`):

| gene | unique_age (corrected) | pool-only | pooled (prior) | unique_ct | slope /yr | class (primary) |
|---|---|---|---|---|---|---|
| SOX4 | 0.0126 | 0.0123 | 0.0215 | 0.447 | −0.0166 | OTHER (ct in gap) |
| RPS24 | 0.0120 | 0.0126 | 0.0006 | 0.378 | −0.0018 | OTHER (ribo-flagged) |
| TSC22D3 | 0.0096 | 0.0098 | 0.0350 | 0.246 | +0.0043 | OTHER (ct in gap) |
| LRRN3 | 0.0080 | 0.0079 | 0.0187 | 0.602 | −0.0119 | **MIXED** |
| PHLDA3 | 0.0055 | 0.0056 | 0.0090 | 0.103 | +0.0040 | **A** |
| SERPINE2 | 0.0053 | 0.0051 | 0.0060 | 0.328 | −0.0063 | OTHER (ct in gap) |
| EIF2S3 | 0.0053 | 0.0054 | 0.0136 | 0.338 | −0.0029 | OTHER (chrX) |
| ETV7 | 0.0052 | 0.0056 | 0.0091 | 0.106 | +0.0064 | **A** |
| LIMD2 | 0.0051 | 0.0043 | 0.0064 | 0.432 | −0.0018 | OTHER (ct in gap) |
| GREM2 | 0.0050 | 0.0050 | 0.0098 | 0.305 | +0.0045 | OTHER |
| KLF6 | 0.0049 | 0.0055 | 0.0068 | 0.510 | +0.0044 | OTHER (t_age missed by 0.0001) |
| CDKN2A | 0.0048 | 0.0050 | 0.0100 | 0.287 | +0.0084 | OTHER |
| PODXL2 | 0.0047 | 0.0046 | 0.0042 | 0.079 | −0.0066 | OTHER |
| DUXAP8 | 0.0046 | 0.0047 | 0.0065 | 0.046 | +0.0049 | OTHER |
| CR1 | 0.0046 | 0.0045 | 0.0054 | 0.621 | −0.0087 | OTHER |

The prior run's top pooled age genes (ANKRD12, EIF1, SET, RIPOR2, UQCRB, DCK, HNRNPH1 …) are absent from this list; 138 of the 140 prior A-genes are
now OTHER, and the 2 that remain A (PHLDA3, ETV7) were 2 of the 4 genes that survived the prior run's union criterion.

### A1-null. Within-pool donor-level permutation (1,000 draws; `results/tables/rerun_stageA_permutation_null_counts.csv`)
Donor ages permuted within pool at the donor level; `unique_age` recomputed for all genes in closed form (verified equal to the QR partition to 1e-9).

| t_age | observed n | null mean | null 95 % | null max | empirical FDR |
|---|---|---|---|---|---|
| 0.001 | 666 | 29.3 | 47 | 82 | 0.044 |
| 0.002 | 118 | 2.5 | 10 | 20 | 0.021 |
| 0.005 | 9 | 0.12 | 1 | 9 | 0.014 |
| 0.010 | 2 | 0.00 | 0 | 0 | 0.000 |

Family-wise: the null maximum over genes has median 0.0027, 95th percentile **0.00525**; **7 genes** exceed it (SOX4, RPS24, TSC22D3, LRRN3, PHLDA3,
SERPINE2 0.00528, EIF2S3 0.00527); ETV7 (0.00521) and LIMD2 (0.00508) fall just short. Per-gene BH on permutation p-values: **1,180 genes at q < 0.05**, 0 at q < 0.01
(the permutation floor p = 1/1001 is hit by 694 genes, so q cannot fall below ≈ 0.02 — a resolution limit, not an absence of signal). The nominal F
test gives 1,902 genes at q < 0.05 and 1,020 at q < 0.01 — anti-conservative, reported for orientation only.

### A2. Gene sets on the corrected partition (thresholds unchanged)

| class | corrected model, primary | prior model, primary (recomputed) |
|---|---|---|
| **A** | **2** — PHLDA3, ETV7 | 140 |
| I | 1,362 | 1,164 (overlap 1,111) |
| MIXED | 1 — LRRN3 | 41 |
| OTHER | 12,539 | 12,559 |

A-gene set sizes across the sweep (rows `t_age`, columns `k_ct`; I thresholds fixed) — prior pooled model in brackets:

| t_age \ k_ct | 0.10 | 0.20 | 0.30 | 0.50 |
|---|---|---|---|---|
| 0.002 | 24 (328) | 42 (700) | 65 (872) | 98 (1,096) |
| **0.005** | 0 (65) | **2 (140)** | 3 (170) | 8 (242) |
| 0.010 | 0 (13) | 0 (30) | 0 (39) | 2 (64) |
| 0.020 | 0 (6) | 0 (10) | 0 (11) | 0 (19) |

MIXED (k_ct 0.2, k_age 0.001): 53 / 20 / 5 at t_age 0.002 for t_ct 0.3 / 0.5 / 0.7; 6 / 1 / 0 at 0.005. I-genes (t_age 0.005, k_ct 0.2): 1,237–1,428
at t_ct 0.5, 502–552 at t_ct 0.7 across k_age.

**Where the within-pool age genes sit on the identity axis** (`results/tables/rerun_stageA_age_genes_by_identity_axis.csv`):

| gene set | n | unique_ct ≤ 0.2 (A-eligible) | 0.2–0.5 (gap) | unique_ct ≥ 0.5 (MIXED-eligible) |
|---|---|---|---|---|
| all expressed genes | 13,904 | 9,423 (67.8 %) | 3,052 (22.0 %) | 1,429 (10.3 %) |
| permutation q < 0.05 | 1,180 | 527 (44.7 %) | 352 (29.8 %) | 301 (**25.5 %**, 2.5× enriched) |
| unique_age ≥ 0.002 | 118 | 42 | 56 | 20 |
| unique_age ≥ 0.005 | 9 | 2 | 6 | 1 |

Sum of `unique_age` over genes by identity bin: ≤ 0.2 → 2.30 (64 %), 0.2–0.5 → 0.94 (26 %), ≥ 0.5 → 0.36 (10 %) of a total 3.59 (prior model total:
10.27). The low-identity bin holds most of the total only because it holds 68 % of the genes; per gene, the strongest within-pool age effects
concentrate in the gap and MIXED bins.

### A3. Known positives (`results/tables/rerun_stageA_known_positives.csv`)

| gene | expected | slope corrected (/yr) | slope pooled | direction | unique_age corrected (rank) | unique_age pooled (rank) | unique_ct | q (perm) | class |
|---|---|---|---|---|---|---|---|---|---|
| **CDKN2A** (gating) | + | **+0.0084** | +0.0106 | ✓ | 0.0048 (**12**) | 0.0100 (77) | 0.287 | 0.020 | OTHER (ct > k_ct) |
| **LRRN3** (gating) | − | **−0.0119** | −0.0159 | ✓ | 0.0080 (**4**) | 0.0187 (21) | 0.602 | 0.020 | MIXED |
| CD248 | − | — | — | — | not expressed (CPM filter) | | | | |
| NELL2 | − | −0.0004 | −0.0003 | ✓ | 0.0000 (12,212) | 0.0000 (13,001) | 0.687 | 0.87 | I |
| LEF1 | − | −0.0022 | −0.0059 | ✓ | 0.0002 (5,560) | 0.0017 (1,527) | 0.702 | 0.14 | I |
| CCR7 | − | −0.0044 | −0.0027 | ✓ | 0.0008 (967) | 0.0004 (5,259) | 0.693 | 0.020 | I |
| GZMH | + | +0.0027 | +0.0065 | ✓ | 0.0002 (5,890) | 0.0012 (2,270) | 0.786 | 0.07 | I |
| B3GAT1 | + | +0.0024 | +0.0037 | ✓ | 0.0007 (1,379) | 0.0021 (1,139) | 0.564 | 0.03 | I |
| KLRG1 | + | −0.0013 | −0.0025 | ✗ | 0.0001 (8,769) | 0.0002 (6,949) | 0.746 | 0.46 | I |

Both gating positives are recovered with the expected direction, at q = 0.02, ranked 4th and 12th genome-wide — they *rise* in rank under the
correction (LRRN3 21 → 4, CDKN2A 77 → 12), i.e. the correction removes artefact around them rather than removing them. **A3 passes; the correction is
not over-aggressive.** Per-cell-type within-pool r (`…_per_celltype.csv`): LRRN3 declines with age most in CD4 Naive (r = −0.42), CD8 TCM (−0.26),
CD8 Naive (−0.26), gdT (−0.23), CD4 TCM (−0.19); CDKN2A rises most in gdT (+0.24), CD8 TCM (+0.24), NK (+0.21), CD8 TEM (+0.16), CD4 TCM (+0.14).
The remaining literature genes (context only) are all I-genes with negligible within-pool age variance: their well-replicated *bulk-blood* age
association is compositional (naive → memory / CTL shift), not a within-cell-type program.

### A4. Prior-run confound filters applied to the new A set (`results/tables/rerun_stageA_Agene_survival_primary.csv`, `…_sweep.csv`)

| filter (same lists, same fixed thresholds) | flagged | removed | surviving |
|---|---|---|---|
| cell cycle — Tirosh/proliferation list | 0 | 0 | 2 |
| cell cycle — data-driven (\|r\| ≥ 0.3) | 0 | 0 | 2 |
| ribosomal | 0 | 0 | 2 |
| mitochondrial | 0 | 0 | 2 |
| sex-linked (chrX/Y) | 0 | 0 | 2 |
| sex DE | 0 | 0 | 2 |
| depth-correlated (\|r\| ≥ 0.3) | 0 | 0 | **2** |

Surviving: **PHLDA3** (+0.0040/yr; unique_age 0.0055; unique_ct 0.103; mean log₂-CPM 0.30) and **ETV7** (+0.0064/yr; 0.0052; 0.106; 0.80). Both are
low-expression genes (PHLDA3 is a p53 target; ETV7 an interferon-inducible ETS factor). Survivors across the sweep (`t_age` × `k_ct`):
16 / 34 / 55 / 82 at 0.002; **0 / 2 / 2 / 5** at 0.005; 0 / 0 / 0 / 1 at 0.01; 0 at 0.02 — essentially all survivors also have permutation q < 0.05
(14 / 32 / 53 / 80 at 0.002). The prior run's mechanism filter (age-variance retention within pool) is moot here: the model is within-pool.

### STOP CONDITION A — VERDICT (`results/rerun/stageA_STOP_verdict.txt`)

```
corrected-model A-genes at primary thresholds: 2
A-genes surviving all confound filters:       2   (STOP A fires if < 10)
of these with within-pool permutation q<0.05:  2
genes with unique_age >= 0.005 regardless of cell-type variance: 9 (permutation null mean 0.12)
genes with ANY within-pool age association (permutation BH q<0.05): 1180
known positives: CDKN2A +0.0084/yr (expected +), rank 12, unique_ct 0.287, OTHER; LRRN3 -0.0119/yr (expected -), rank 4, unique_ct 0.602, MIXED
VERDICT: STOP CONDITION A FIRED.
```

The margin is not a threshold artefact: 2 survivors against a stop line of 10, and the number does not reach 10 anywhere on the `t_age = 0.005` row of
the sweep (max 5 at `k_ct = 0.5`). Per the operating mode the rerun halts here; Stages B, C, D were not run and no threshold was relaxed.

---

## Prior run vs corrected model — side by side

| quantity | prior run (`ct + age`, FINDINGS.md) | corrected (`ct + pool + depth + age`, this run) |
|---|---|---|
| pseudobulks × genes | 9,623 × 13,904 | same (reused) |
| model parameters | 17 | 93 |
| unique age variance: max / 99.9 % / 99 % / 90 % / median | 0.0350 / 0.0215 / 0.0074 / 0.0018 / 0.0002 | **0.0126 / 0.0046 / 0.0019 / 0.0007 / 0.0001** |
| genes with unique_age ≥ 0.002 / 0.005 / 0.01 / 0.02 | 1,256 / 283 / 76 / 19 | **118 / 9 / 2 / 0** |
| sum of unique_age over all genes | 10.27 | **3.59** (−65 %) |
| top age genes | ANKRD12, TSC22D3, EIF1, SET, RIPOR2, UQCRB, DCK, HNRNPH1 … (depth genes) | SOX4, RPS24, TSC22D3, LRRN3, PHLDA3, SERPINE2, EIF2S3, ETV7, LIMD2 |
| unique cell-type variance: max / 99 % / median | 0.969 / 0.854 / 0.126 | 0.961 / 0.834 / 0.125 (Spearman 0.9997) |
| A-genes (primary) | 140 | **2** |
| I-genes (primary) | 1,164 | 1,362 (overlap 1,111) |
| MIXED (primary) | 41 | 1 (LRRN3) |
| A-genes surviving confound filters | 105 by literal flags; 4 by the union criterion (ETV7, PHLDA3, CNN2, DUXAP8) | **2 (PHLDA3, ETV7); 0 removed** |
| CDKN2A: rank / slope / class | 77 / +0.0106 / OTHER | 12 / +0.0084 / OTHER (unique_ct 0.287) |
| LRRN3: rank / slope / class | 21 / −0.0159 / MIXED | 4 / −0.0119 / MIXED |
| inference on the age term | none (set sizes only) | within-pool donor-level permutation: 1,180 genes q < 0.05; FWER-5 % threshold 0.0053 → 7 genes |
| stop condition | STOP 1b fired (A set is a depth artefact) | **STOP A fired (2 < 10 A-genes)** |
| stages not run | 1c, 1d | B, C, D |

The prior run's A-gene inflation was entirely batch: 138/140 prior A-genes are OTHER under the corrected model, and the only two that remain are the
two low-identity genes the prior run's own union criterion had already isolated.

---

## What the halt establishes, and what it does not (assessment)

**1. The A-gene construct — high age variance *and* low identity variance per gene — is nearly empty in PBMC once batch is modelled.** This is the
result the spec anticipated ("EXPECT FAR FEWER A-genes") but at the low end of the projection (2, not 10–100). It is a property of the data at the
pre-registered thresholds, not of the pipeline: the permutation null shows the 9 genes above `t_age` are real (null mean 0.12), and the positives
are recovered and improved in rank. What fails is the conjunction: the genes with the strongest within-pool age effects (SOX4, TSC22D3, LRRN3,
SERPINE2, KLF6, CR1, CDKN2A, SELL, BCL2, LGALS1) carry 25–82 % identity variance.

**2. "Too sparse" means "too weak per gene", not "absent".** 1,180 genes carry a detectable within-pool age association (BH q < 0.05 against a
donor-level within-pool null), but the median such gene explains 0.10 % of variance and the strongest 1.26 %. Age in PBMC within cell type is a
diffuse, many-small-effects signal. A multivariate score over many weak genes may still recover age (bulk-blood clocks reach R² ≈ 0.3–0.5), but a
*gene-set-restricted* score built from genes that are individually age-specific cannot be assembled from 2 genes. That was the question STOP A
gates, and the answer is no.

**3. The age genes lean toward identity — measured, not inferred.** The 2.5× enrichment of q < 0.05 age genes among high-identity genes, the
6-of-9 top genes in the identity gap, LRRN3 as the top MIXED gene, and the fact that every classic bulk-blood age marker tested (CCR7, LEF1, NELL2,
GZMH, B3GAT1, KLRG1) is an I-gene with no within-type age variance all point the same way: much of what bulk transcriptomics calls "aging" in blood
is a change in *which* cells are present and in subset-defining genes. This is exactly the Stage B question. **Stage B was not run** — it is
gated behind STOP A in the protocol — so the compositional share of age (B1/B2) and the MIXED-vs-A share of the within-pool age signal (B3) remain
unmeasured. The observations here are consistent with a largely compositional signal but do not quantify it.

**4. Depth covariates are redundant with pool on this dataset.** Adding log mean UMI/cell and log mean genes/cell to a model that already has pool
changes `unique_age` negligibly (Spearman 0.982; counts 10 → 9 at 0.005). Pool captures depth (R² = 0.86). In a cohort where depth varies *within*
batch this would not hold; the two covariates should stay in the model as a matter of design, but on OneK1K the operative correction is pool.

**5. Cautionary observations, reported not interpreted.** (a) 78 % of the 1,180 significant within-pool age genes have negative slopes. A genuine
global down-regulation with age within cell types is possible, but so is a residual donor-level technical effect (e.g. RNA content or ambient
fraction) that pool and mean depth do not fully absorb. This asymmetry should be checked before any of these genes are called biology.
(b) RPS24 (ribosomal, depth-flagged) is the 2nd-ranked within-pool age gene; it would have been removed by the ribosomal filter had it been
A-eligible, and its presence near the top is a reminder that the within-pool model does not remove every technical axis. (c) PHLDA3 and ETV7 are
expressed at log₂-CPM 0.3 and 0.8 — low-abundance genes in which pseudobulk log-CPM is noisy for small pseudobulks; two such genes are not a score.

**6. The threshold question, stated openly and not acted on.** `t_age = 0.005` was calibrated on the pooled distribution (top ~2.5 %). On the
corrected distribution it is the top 0.06 %; the 99th-percentile-matched value would be ≈ 0.002, where 34 A-genes survive at `k_ct = 0.2` (32 with
q < 0.05). Re-centring the threshold on the new distribution would have avoided STOP A. It was not done, because the rerun spec forbids it and
because the alternative is a moving target: the prior run already re-centred once. If a future run wishes to use quantile-anchored thresholds it
must pre-register them *before* seeing the corrected distribution — which is no longer possible for OneK1K, but is possible for AIDA.

---

## Not run, and what a continuation would need (recommendations, not executed)

Not run: **Stage B** (B1 composition ~ age within pool; B2 composition vs within-type expression under pool-grouped CV; B3 share of within-pool age
signal in MIXED vs A), **Stage C** (A-gene and I-gene scores, unrestricted all-gene age model), **Stage D** (shuffles, cross-test, AIDA audit and
replication). AIDA was not opened. `FALSIFICATION.md` remains un-evaluated and unchanged.

A continuation is a decision for the protocol owner, not for this run. Observations relevant to that decision:
1. **B1 and B2 do not depend on the A-gene set.** Composition ~ age (within pool) and an all-gene within-type age model under pool-grouped CV
   (`rerun_common.pool_grouped_folds` / `assert_disjoint` already enforce HARD RULE 1) can be run without touching any threshold. They would
   quantify what §3 above only indicates. B3 depends on the A set and would be reported as 2 genes vs 1 MIXED gene at the primary thresholds —
   i.e. it is answered trivially and uninformatively by the current sets.
2. If the gene-set construct is kept, the honest route is a **new pre-registration** of quantile-anchored thresholds on a dataset whose corrected
   distribution has not yet been seen (AIDA, after the same pool ↔ age audit), not a re-centring on OneK1K.
3. The unrestricted (all-gene) age model of Stage C is the right diagnostic for "diffuse but real": with 1,180 weak genes, ridge under pool-grouped
   CV would show whether they integrate to a usable age predictor within cell type. STOP C's "unrestricted model at chance" branch is live.
4. Before interpreting any within-pool age gene as biology, check the 78 % negative-slope asymmetry against a donor-level RNA-content covariate.

---

## Reproducibility
- Clean checkout, Windows PowerShell (adjust paths on POSIX):
  `python -m venv .venv; .\.venv\Scripts\pip install -r requirements.txt; .\.venv\Scripts\python -m ipykernel install --user --name age-identity-venv;`
  `$env:PYTHONUTF8=1; .\.venv\Scripts\python src\rerun_make_notebooks.py;`
  `.\.venv\Scripts\python -m nbconvert --to notebook --execute --inplace --ExecutePreprocessor.kernel_name=age-identity-venv --ExecutePreprocessor.timeout=1800 notebooks\rerun_stageA_age_axis.ipynb`
  The notebook rebuilds the pseudobulk from the cached OneK1K h5ad if `data/processed/` is absent (`rerun_common.ensure_logcpm`, via the reused
  `src/run_pseudobulk_onek1k.py`), fetches chromosome annotations from Ensembl if the cache is absent, and re-asserts the STOP A verdict.
  Wall time ≈ 3.5 min for Stage A on the cached pseudobulk (1,000 permutations included).
- Environment as the prior run: Python 3.12 venv; numpy 2.5.3, pandas 3.0.5, scipy 1.18.1, statsmodels 0.15.0, scikit-learn 1.9.1 (`requirements.txt`).
- Seed `20260914` (`rerun_common.RERUN_SEED`); the only randomness in Stage A is the permutation null (`numpy.random.default_rng(RERUN_SEED)`).
- Nothing from AIDA was read.

## Files produced by this run
| path | content |
|---|---|
| `src/rerun_common.py` | corrected model, commonality partition, within-pool permutation, pool-grouped fold builder with donor+pool disjointness assertion |
| `src/rerun_stageA.py` | Stage A runner (A0–A4, permutation null, STOP A verdict, figures) |
| `src/rerun_make_notebooks.py` | generates `notebooks/rerun_stageA_age_axis.ipynb` |
| `notebooks/rerun_stageA_age_axis.ipynb` | executed notebook (all console output embedded) |
| `results/rerun/stageA_report.txt`, `stageA_STOP_verdict.txt`, `stageA_summary.json`, `stageA_null_max_unique_age.npy`, `nb_stageA_exec.log` | logs and verdict |
| `results/tables/rerun_stageA_variance_decomposition.csv` | per-gene partition (corrected + prior + pool-only), slopes, nominal and permutation p/q, ranks, classes |
| `results/tables/rerun_stageA_unique_age_quantiles_prior_vs_corrected.csv`, `…_permutation_null_counts.csv`, `…_threshold_sweep_set_sizes.csv`, `…_age_genes_by_identity_axis.csv` | distributions, null, sweep |
| `results/tables/rerun_stageA_known_positives.csv`, `…_known_positives_per_celltype.csv` | A3 |
| `results/tables/rerun_stageA_gene_flags.csv`, `…_Agene_survival_primary.csv`, `…_Agene_survival_sweep.csv`, `…_A_and_MIXED_genes_primary.csv` | A4 |
| `results/figures/rerun_stageA_pooled_vs_corrected_unique_age.png` | per-gene age variance before vs after correction; rank profiles vs null |
| `results/figures/rerun_stageA_joint_variance_distribution.png` | identity vs within-pool age axis, literature genes circled |
| `results/figures/rerun_stageA_permutation_null.png` | observed counts vs null; family-wise null of the maximum |
| `results/figures/rerun_stageA_Agene_survival.png` | filter survival at primary thresholds; survivors across the sweep |
