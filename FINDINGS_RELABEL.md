# FINDINGS_RELABEL — do the headline numbers survive Lu-style joint clustering?

## Answer

**Yes. Under the joint clustering the pre-registered rule is met: yes. Under Lu et al.'s own published labels it is met: yes.** In both the full and MD gene spaces the aged PartialReprog cells still end farther from the young target, the mixture control still passes, and the asymmetry interval stays far from 0.

- Joint clustering changed the label of 7.8% of aged cells and 8.3% of young cells.
- Same-platform, full space: forward progress 0.6744 → 0.6493; change in distance 130.894 → 133.303; asymmetry 0.7317 → 0.6742 [0.6630, 0.6812].
- Same-platform, MD space: forward progress 0.7901 → 0.7683; change in distance 17.0880 → 17.4599; asymmetry 1.0196 → 0.9761 [0.9492, 0.9995].
- MD score drop ÷ gap: 7.72 → 7.85 [7.4480, 8.3518] (joint); 8.03 [7.5828, 8.4349] (Lu labels).
- GTEx, full space: the change in distance to young GTEx is 9.3697 → 9.4073 (joint), 15.7867 (Lu labels); still > 0.

Under the joint clustering, forward progress, change in distance, final ÷ start, asymmetry and the mixture numbers each moved by less than 8% of their old value. Reverse progress, which is close to 0, moved from -0.0573 to -0.0248 (full) and -0.2295 to -0.2078 (MD). Every sign and every gate verdict is unchanged. Lu's labels move the numbers somewhat more, with the largest shift in the GTEx distances (details below); that doesn't change a sign or a verdict either.

**Sentence the paper could use:** "Re-clustering both donors together, as Lu et al. did, changed the state label of about 8% of cells and left the result unchanged: aged partially reprogrammed cells still ended farther from young fibroblasts than they started (change in distance 133.3 across all genes and 17.5 across MD genes, versus 130.9 and 17.1 with the original labels), and the same held with Lu et al.'s published per-cell state labels."

## Pre-registration (written before the joint clustering and before any statistic)

Frozen copy: results/relabel_check/PREREG.flag. Not rewritten after this point.

### What Step 0 found (done before this block)

Lu et al. deposited two Seurat objects as GEO GSE297234 supplementary files
(GSE297234_HFIB_COMBINED_SEVOSKM.rds, GSE297234_GM00731_SEVOSKM.rds; the third file is the
10x RAW tar we already use). The paper has no separate data-availability paragraph; its key
resources table points to GEO GSE297234. The combined object (both donors, 53,454 cells, active
assay SCT) has a per-cell `cell_state` column with 16 sub-states (Fib_1-3, ParRep_1/2/4 and
ParRepr_3, EarPlu_1-3, Plu_1-2, NonRep_1-4). These are Lu et al.'s published per-cell labels.
The GM00731-only object has cluster letters but no `cell_state`. Cell names are
`<line>_D<day>_<barcode>`, which match our (cell_line, day, barcode).
Disclosure: before writing this block I tabulated Lu's cell_state x sample cell counts (counts
only; no expression, score or distance was computed with any new label).

### Labellings

- OLD: md2 labels (per-donor Louvain, results/md2/t2_cluster_labels_*.csv).
- JOINT (primary; answers the question asked): all cells of both donors in one matrix, md2
  Louvain pipeline unchanged: Gene Expression features; cells with >500 detected genes; genes
  detected in >=5 cells of the merged matrix; LogNormalize (scale 1e4); 3000 HVGs by the md2
  quadratic log-mean/log-variance rule; each HVG regressed on percent-mito and z-scored;
  50 PCs (randomized, random_state 20260914); SNN k=20, prune 1/15; networkx Louvain
  resolution 0.8, seed 20260914. Labels: md2 AddModuleScore state signatures
  (results/md2/genesets.json; nbin 24, ctrl 100, seed 20260914) computed on the joint matrix
  (one set of expression bins over both donors), then src/md2_task2.py::_label_clusters
  (each joint cluster gets the state with the highest mean score over all its cells).
  Row order: aged libraries (day 0, 3, 7, 10) then young.
- LU (second labelling): Lu's `cell_state` collapsed by prefix: Fib_ -> Fibroblast,
  ParRep_/ParRepr_ -> PartialReprog, EarPlu_ -> EarlyPluripotency, Plu_ -> Pluripotency,
  NonRep_ -> NonReprog. Our cells absent from Lu's object are unlabelled and belong to no
  state.

Donors are pooled only to cluster and label. Every statistic stays per donor.

### Memory deviation and its gate

The unchanged code (src/md2_cluster.py::louvain_one_donor) needs several full float64 copies
of the merged matrix and does not fit in this machine's memory for ~56k cells. The joint run
therefore uses a streamed version meant to perform the same arithmetic in the same order
(per-gene means accumulated row by row in the same order scipy uses; per-cell steps done per
donor; the HVG block residualized in column blocks). No cell or gene is subsampled.
Gate: the same streamed code is first run on GM00731 alone and must reproduce md2 exactly
(identical cluster of every cell in results/md2/cluster_labels_louvain_GM00731.csv, identical
state scores in data/processed/md2/louvain_state_GM00731.npz, identical labels in
results/md2/t2_cluster_labels_GM00731.csv). If it does not, the joint run is not done and I
stop and report.

### Statistics (per donor, same functions and seeds)

a. MD score per state (src/newstory.py Part A joint scoring): the per-cell joint MD scores in
   results/newstory/partA_cell_scores.npz do not depend on labels and are reused; state means
   and bootstrap CIs (B=200, Generator(20260918), donors then states in fixed order) are
   recomputed with the new masks. Reported: aged day-0 Fibroblast, aged PartialReprog, young
   day-0 Fibroblast, drop / gap.
b. Same-platform test (src/genespace_run.py sections 2 and bootstrap, FULL and MD spaces):
   O, Y = half A of day-0 Fibroblast cells of the aged and young donor; S = aged PartialReprog;
   S_rev = young PartialReprog. Forward progress, change in distance (delta), final / start
   distance, reverse progress, asymmetry (forward minus reverse progress) with its 95%
   percentile bootstrap interval (B=200, Generator(20260918), O and Y frozen), and the mixture
   control (genespace_run.mixture_gate: progress non-decreasing in f = 0, 0.10, 0.25, 0.50;
   progress CI at f=0.50 excludes 0 and contains the point; delta at f=0.50 < 0).
   Split rule: if the set of day-0 Fibroblast cells of either donor differs from the old one,
   the half-split is redone with same_run.split_half(seed 20260914) on the new set (the
   mixtures are rebuilt from the new halves with the same seed), and this is stated.
c. GTEx test (src/toward_run.py donor panel, genespace_run.gtex_point_block, FULL space,
   aged PartialReprog): cosine, change in distance to young GTEx, change in distance to old.

Harness check before any new-label number: the same code is run with the OLD labels and must
reproduce the stored values (results/newstory/partA_*.csv, results/genespace/same_stats.csv,
same_posctrl.csv, gtex_stats.csv; points and CI bounds) to within 1e-9. If not, stop.

### Decision rule (verbatim)

"The conclusion holds if, in both the full and MD gene spaces, the aged PartialReprog cells end farther from the young target (change in distance > 0), the mixture control passes, and the forward-minus-reverse asymmetry interval excludes 0."

Applied to each labelling separately (JOINT is the answer to the question; LU is reported
alongside, not averaged). "Full" and "MD" are the FULL (23,485 ruler genes) and MD (205
genes) spaces of src/genespace_run.py; "change in distance" is the forward delta; the
interval is the 95% percentile bootstrap interval. If a state needed by the rule is missing
or too small for the code to run, the rule is counted as not met for that labelling.
How far numbers moved is reported descriptively (old, new, difference); no threshold.


## Checks before the new-label numbers

- **Clustering gate (GM00731 alone, streamed code).** Passed: cells and order identical, cluster of every cell identical (19 clusters, 808,714 SNN edges, as in md2), gene means and AMS bins identical, MD and all five state scores identical (max |diff| 0), cluster labels identical. percent-mito max |diff| 7.1e-15 (CSV round-trip of the stored value).
- **Harness (old labels through the new code).** 112/112 stored values reproduced within 1e-9 (max |diff| 5.7e-14): newstory Part A state means, CIs, gap/drop/ratio; genespace same_stats and same_posctrl (FULL, MD; points and CIs); GTEx FULL aged PartialReprog. See results/relabel_check/harness_check.csv.
- **Memory.** The joint run used all 55,990 cells (aged 22,132, young 33,858) and 30,804 genes (detected in ≥5 cells of the merged matrix; 5,797 dropped). Nothing was subsampled.

## Step 1 — joint clustering

24 joint clusters (resolution 0.8; 2,043,044 SNN edges; 50 PCs explain 0.227 of HVG variance). 11 of 24 clusters have 20–80% aged cells; the rest are mostly one donor, so the donors still separate in several clusters even when clustered together. Each cluster's label, its mean state scores, the margin between the top two scores and its donor mix:

| cluster | n_cells | label | top-2 score margin | aged cells | young cells |
|---|---|---|---|---|---|
| 0 | 6348 | Fibroblast | 0.3186 | 11 | 6337 |
| 1 | 4446 | PartialReprog | 0.2026 | 2066 | 2380 |
| 2 | 4406 | EarlyPluripotency | 0.3151 | 191 | 4215 |
| 3 | 4397 | Fibroblast | 0.4728 | 4353 | 44 |
| 4 | 4047 | EarlyPluripotency | 0.4364 | 41 | 4006 |
| 5 | 3462 | EarlyPluripotency | 0.2691 | 62 | 3400 |
| 6 | 3347 | EarlyPluripotency | 0.5022 | 3316 | 31 |
| 7 | 3114 | NonReprog | 0.4148 | 2976 | 138 |
| 8 | 2978 | PartialReprog | 0.1225 | 874 | 2104 |
| 9 | 2441 | PartialReprog | 0.1969 | 1709 | 732 |
| 10 | 2065 | NonReprog | 0.1039 | 2008 | 57 |
| 11 | 2064 | Fibroblast | 0.0147 | 14 | 2050 |
| 12 | 1977 | Pluripotency | 0.05312 | 636 | 1341 |
| 13 | 1845 | EarlyPluripotency | 0.004758 | 493 | 1352 |
| 14 | 1733 | NonReprog | 0.2295 | 30 | 1703 |
| 15 | 1532 | NonReprog | 0.1521 | 905 | 627 |
| 16 | 1425 | PartialReprog | 0.1458 | 625 | 800 |
| 17 | 1361 | Fibroblast | 0.2959 | 15 | 1346 |
| 18 | 841 | PartialReprog | 0.1739 | 449 | 392 |
| 19 | 692 | EarlyPluripotency | 0.2228 | 416 | 276 |
| 20 | 580 | Fibroblast | 0.4806 | 580 | 0 |
| 21 | 414 | EarlyPluripotency | 0.04064 | 227 | 187 |
| 22 | 316 | EarlyPluripotency | 0.5156 | 20 | 296 |
| 23 | 159 | EarlyPluripotency | 0.3875 | 115 | 44 |

Clusters whose top-2 margin is below 0.01 (the label is close to a tie): c13 (EarlyPluripotency, n=1845, margin 0.0048).

## Step 2 — agreement (per donor)

### GM00731 (aged donor)

Cell counts, old md2 label (rows) × joint label (columns). Share of cells with an unchanged label: **92.2%** of 22,132.

|  | Fibroblast | PartialReprog | EarlyPluripotency | Pluripotency | NonReprog |
|---|---|---|---|---|---|
| Fibroblast | 4947 | 3 | 546 | 1 | 491 |
| PartialReprog | 6 | 5632 | 58 | 33 | 103 |
| EarlyPluripotency | 0 | 53 | 4181 | 142 | 38 |
| Pluripotency | 0 | 1 | 13 | 395 | 30 |
| NonReprog | 20 | 34 | 83 | 65 | 5257 |

Cell counts, old md2 label (rows) × Lu label (columns). Share of cells with an unchanged label: **80.0%** (of 21,033 cells that Lu labelled; 1,099 of our cells are not in Lu's object).

|  | Fibroblast | PartialReprog | EarlyPluripotency | Pluripotency | NonReprog | Unlabelled |
|---|---|---|---|---|---|---|
| Fibroblast | 5023 | 25 | 251 | 0 | 656 | 33 |
| PartialReprog | 149 | 4990 | 160 | 12 | 363 | 158 |
| EarlyPluripotency | 183 | 344 | 2402 | 599 | 485 | 401 |
| Pluripotency | 2 | 2 | 185 | 216 | 2 | 32 |
| NonReprog | 550 | 144 | 99 | 1 | 4190 | 475 |

Cells per state (all days, and day 0):

| state | old | joint | lu | old day 0 | joint day 0 | lu day 0 |
|---|---|---|---|---|---|---|
| Fibroblast | 5988 | 4973 | 5907 | 4939 | 4938 | 4937 |
| PartialReprog | 5832 | 5723 | 5505 | 3 | 3 | 1 |
| EarlyPluripotency | 4414 | 4881 | 3097 | 2 | 0 | 2 |
| Pluripotency | 439 | 636 | 828 | 0 | 7 | 1 |
| NonReprog | 5459 | 5919 | 5696 | 65 | 61 | 39 |
| Unlabelled | 0 | 0 | 1099 | 0 | 0 | 29 |

### GM23815 (young donor)

Cell counts, old md2 label (rows) × joint label (columns). Share of cells with an unchanged label: **91.7%** of 33,858.

|  | Fibroblast | PartialReprog | EarlyPluripotency | Pluripotency | NonReprog |
|---|---|---|---|---|---|
| Fibroblast | 9692 | 39 | 1529 | 6 | 89 |
| PartialReprog | 41 | 6328 | 69 | 12 | 31 |
| EarlyPluripotency | 37 | 28 | 12165 | 225 | 506 |
| Pluripotency | 1 | 9 | 15 | 1098 | 126 |
| NonReprog | 6 | 4 | 29 | 0 | 1773 |

Cell counts, old md2 label (rows) × Lu label (columns). Share of cells with an unchanged label: **71.6%** (of 32,419 cells that Lu labelled; 1,439 of our cells are not in Lu's object).

|  | Fibroblast | PartialReprog | EarlyPluripotency | Pluripotency | NonReprog | Unlabelled |
|---|---|---|---|---|---|---|
| Fibroblast | 7711 | 86 | 300 | 1 | 3218 | 39 |
| PartialReprog | 114 | 5875 | 87 | 46 | 230 | 129 |
| EarlyPluripotency | 518 | 276 | 7061 | 3456 | 778 | 872 |
| Pluripotency | 16 | 27 | 40 | 768 | 5 | 393 |
| NonReprog | 6 | 13 | 2 | 0 | 1785 | 6 |

Cells per state (all days, and day 0):

| state | old | joint | lu | old day 0 | joint day 0 | lu day 0 |
|---|---|---|---|---|---|---|
| Fibroblast | 11355 | 9777 | 8365 | 7703 | 7688 | 7741 |
| PartialReprog | 6481 | 6408 | 6277 | 9 | 17 | 0 |
| EarlyPluripotency | 12961 | 13807 | 7490 | 56 | 1 | 0 |
| Pluripotency | 1249 | 1341 | 4271 | 0 | 14 | 1 |
| NonReprog | 1812 | 2525 | 6016 | 0 | 48 | 3 |
| Unlabelled | 0 | 0 | 1439 | 0 | 0 | 23 |

Day-0 Fibroblast membership (the origin cells O/Y of the same-platform test and the GTEx origin):

| set | n_old | n_new | n_both | only_old | only_new | identical |
|---|---|---|---|---|---|---|
| GM00731_joint | 4939 | 4938 | 4936 | 3 | 2 | False |
| GM00731_lu | 4939 | 4937 | 4893 | 46 | 44 | False |
| GM23815_joint | 7703 | 7688 | 7687 | 16 | 1 | False |
| GM23815_lu | 7703 | 7741 | 7697 | 6 | 44 | False |

Membership changed for both donors under both labellings, so the half-split was redone with same_run.split_half(seed 20260914) on each new day-0 Fibroblast set, and the mixtures were rebuilt from the new halves (same seed).

- joint: O=2469, Y=3844, aged half B=2469, young half B=3844 (old: 2469 / 3851 / 2470 / 3852); panel sizes {'O': 2469, 'Y': 3844, 'S': 5723, 'S_rev': 6408, 'Pluri': 636, 'NonReprog': 5919}.
- lu: O=2468, Y=3870, aged half B=2469, young half B=3871 (old: 2469 / 3851 / 2470 / 3852); panel sizes {'O': 2468, 'Y': 3870, 'S': 5505, 'S_rev': 6277, 'Pluri': 828, 'NonReprog': 5696}.

## Step 3 — headline numbers, old vs new

"new" = joint labelling (primary). Lu-label values are in the last columns. diff = new − old; rel = diff / |old|.

### a. MD score per state (joint AddModuleScore)

| quantity | old | new (joint) | diff | rel | Lu labels | Lu diff |
|---|---|---|---|---|---|---|
| aged day-0 Fibroblast MD score | 0.5151 | 0.5151 | -3.47e-05 | -0.0% | 0.5142 | -9.74e-04 |
| aged PartialReprog MD score | 0.1575 | 0.1530 | -4.55e-03 | -2.9% | 0.1458 | -0.0117 |
| young day-0 Fibroblast MD score | 0.4688 | 0.4689 | 1.41e-04 | +0.0% | 0.4683 | -5.39e-04 |
| gap (aged d0 - young d0) | 0.0463 | 0.0462 | -1.76e-04 | -0.4% | 0.0459 | -4.35e-04 |
| drop (aged d0 - aged PartialReprog) | 0.3576 | 0.3621 | 4.51e-03 | +1.3% | 0.3683 | 0.0107 |
| drop / gap | 7.7191 | 7.8462 | 0.1271 | +1.6% | 8.0255 | 0.3064 |
| drop / gap CI low | 7.3319 | 7.4480 | 0.1162 | +1.6% | 7.5828 | 0.2510 |
| drop / gap CI high | 8.2120 | 8.3518 | 0.1398 | +1.7% | 8.4349 | 0.2228 |

### b. Same-platform test, FULL space

| quantity | old | new (joint) | diff | rel | Lu labels | Lu diff |
|---|---|---|---|---|---|---|
| forward progress | 0.6744 | 0.6493 | -0.0251 | -3.7% | 0.6550 | -0.0195 |
| change in distance (delta) | 130.894 | 133.303 | 2.4091 | +1.8% | 135.834 | 4.9407 |
| delta CI low | 130.473 | 133.238 | 2.7646 | +2.1% | 135.896 | 5.4229 |
| delta CI high | 133.918 | 136.390 | 2.4720 | +1.8% | 139.510 | 5.5917 |
| final / start distance | 2.0016 | 2.0176 | 0.0161 | +0.8% | 2.0269 | 0.0254 |
| reverse progress | -0.0573 | -0.0248 | 0.0325 | +56.7% | -0.0341 | 0.0232 |
| asymmetry | 0.7317 | 0.6742 | -0.0576 | -7.9% | 0.6890 | -0.0427 |
| asymmetry CI low | 0.7240 | 0.6630 | -0.0611 | -8.4% | 0.6746 | -0.0494 |
| asymmetry CI high | 0.7404 | 0.6812 | -0.0592 | -8.0% | 0.7209 | -0.0195 |
| mixture progress at f=0.50 | 0.4408 | 0.4354 | -5.40e-03 | -1.2% | 0.4311 | -9.69e-03 |
| mixture delta at f=0.50 | -34.4430 | -34.0295 | 0.4135 | +1.2% | -34.3325 | 0.1105 |
| mixture control passes | pass | pass |  |  | pass |  |

### b. Same-platform test, MD space

| quantity | old | new (joint) | diff | rel | Lu labels | Lu diff |
|---|---|---|---|---|---|---|
| forward progress | 0.7901 | 0.7683 | -0.0218 | -2.8% | 0.7426 | -0.0475 |
| change in distance (delta) | 17.0880 | 17.4599 | 0.3720 | +2.2% | 17.8387 | 0.7507 |
| delta CI low | 16.7333 | 17.1772 | 0.4439 | +2.7% | 17.4529 | 0.7196 |
| delta CI high | 17.4193 | 17.8396 | 0.4203 | +2.4% | 18.5247 | 1.1054 |
| final / start distance | 2.1706 | 2.1884 | 0.0178 | +0.8% | 2.1881 | 0.0175 |
| reverse progress | -0.2295 | -0.2078 | 0.0216 | +9.4% | -0.1906 | 0.0388 |
| asymmetry | 1.0196 | 0.9761 | -0.0434 | -4.3% | 0.9333 | -0.0863 |
| asymmetry CI low | 1.0094 | 0.9492 | -0.0602 | -6.0% | 0.9015 | -0.1079 |
| asymmetry CI high | 1.0547 | 0.9995 | -0.0552 | -5.2% | 1.0266 | -0.0281 |
| mixture progress at f=0.50 | 0.3883 | 0.3793 | -8.98e-03 | -2.3% | 0.4115 | 0.0232 |
| mixture delta at f=0.50 | -5.1430 | -5.1250 | 0.0180 | +0.3% | -5.3603 | -0.2173 |
| mixture control passes | pass | pass |  |  | pass |  |

### c. GTEx test, FULL space, aged PartialReprog

| quantity | old | new (joint) | diff | rel | Lu labels | Lu diff |
|---|---|---|---|---|---|---|
| cosine | 0.0772 | 0.0778 | 5.85e-04 | +0.8% | 0.0815 | 4.26e-03 |
| change in distance to young GTEx | 9.3697 | 9.4073 | 0.0376 | +0.4% | 15.7867 | 6.4170 |
| change in distance to old GTEx | 10.8571 | 10.9147 | 0.0576 | +0.5% | 17.3980 | 6.5410 |

Mixture control detail (progress at f = 0, 0.10, 0.25, 0.50; gate = non-decreasing, CI at 0.50 excludes 0 and contains the point, delta at 0.50 < 0):

| labelling | space | progress | ci50 | delta50 | gate |
|---|---|---|---|---|---|
| old | FULL | 0.1146, 0.1824, 0.2829, 0.4408 | [0.4243, 0.4475] | -34.4430 | pass |
| old | MD | 0.0104, 0.0955, 0.2074, 0.3883 | [0.3671, 0.4072] | -5.1430 | pass |
| joint | FULL | 0.1094, 0.1771, 0.2713, 0.4354 | [0.4168, 0.4442] | -34.0295 | pass |
| joint | MD | -0.0007, 0.0933, 0.2094, 0.3793 | [0.3538, 0.3985] | -5.1250 | pass |
| lu | FULL | 0.1102, 0.1781, 0.2737, 0.4311 | [0.4149, 0.4385] | -34.3325 | pass |
| lu | MD | 0.0338, 0.1033, 0.2096, 0.4115 | [0.3824, 0.4345] | -5.3603 | pass |

MD-score CIs (bootstrap B=200): old: gap [0.0434, 0.0487], drop [0.3530, 0.3628], drop/gap [7.3319, 8.2120]; joint: gap [0.0432, 0.0486], drop [0.3579, 0.3664], drop/gap [7.4480, 8.3518]; lu: gap [0.0436, 0.0484], drop [0.3642, 0.3726], drop/gap [7.5828, 8.4349].

## Decision rule, applied

"The conclusion holds if, in both the full and MD gene spaces, the aged PartialReprog cells end farther from the young target (change in distance > 0), the mixture control passes, and the forward-minus-reverse asymmetry interval excludes 0."

| labelling | space | change in distance | > 0 | mixture control | asymmetry interval | excludes 0 | met |
|---|---|---|---|---|---|---|---|
| old | FULL | 130.894 | True | pass | [0.7240, 0.7404] | True | True |
| old | MD | 17.0880 | True | pass | [1.0094, 1.0547] | True | True |
| joint | FULL | 133.303 | True | pass | [0.6630, 0.6812] | True | True |
| joint | MD | 17.4599 | True | pass | [0.9492, 0.9995] | True | True |
| lu | FULL | 135.834 | True | pass | [0.6746, 0.7209] | True | True |
| lu | MD | 17.8387 | True | pass | [0.9015, 1.0266] | True | True |

Joint labelling (primary): **holds**. Lu labels: **holds**. (Old labels, for reference: holds.)

## Plain summary

- **Does the conclusion hold?** Yes, under both relabellings, in both gene spaces. Clustering the two donors together, as Lu et al. did, does not change what the aged partially reprogrammed cells do: part of their movement runs in the young direction, but they end about twice as far from the young day-0 cells as they started (final ÷ start 2.02 in all genes, 2.19 in MD genes; old 2.00 and 2.17).
- **How far did the numbers move?** With the joint clustering, about 8% of cells changed label. Same-platform forward progress fell by about 0.02–0.03, the change in distance rose by about 2% in both spaces, and the asymmetry fell by about 0.04–0.06, staying near 0.7 (full) and 1.0 (MD) with intervals far from 0. The MD score drop ÷ gap moved from 7.72 to 7.85. The GTEx numbers barely moved (cosine +0.0006, distance changes +0.04 and +0.06).
- **Lu's published labels** agree less with ours (72–80% of cells keep their label) and move the numbers a little further: same-platform deltas rise by about 4%, the MD asymmetry falls by 0.09, and the GTEx distance changes rise from 9.3697 / 10.8571 to 15.7867 / 17.3980 (young / old). That GTEx shift is the largest move anywhere, and it goes in the direction of *more* distance, not less. Lu's aged PartialReprog set is more purely day 3 than ours (cells by day: Lu {0: 1, 3: 5258, 7: 207, 10: 39}, old {0: 3, 3: 5220, 7: 578, 10: 31}, joint {0: 3, 3: 5166, 7: 546, 10: 8}). I did not test whether that difference causes the GTEx shift.
- **One sentence for the paper:** see the Answer section above.

## Notes

- Joint labelling reuses md2 code for the whole pipeline. The only deviation is the streamed memory layout, and it reproduced md2 exactly on GM00731 alone. As in md2, this is LogNormalize plus percent-mito regression, not Lu's sctransform v2, so this answers "clustered together" rather than "clustered with Lu's exact normalization". The Lu-label run uses Lu's own clustering, which was built on sctransform.
- The MD score in part a is the newstory joint per-cell score, which does not depend on labels; only the state memberships changed.
- results/verify/claims.csv (paper sentence P42/P96) describes the analysis as using "the published cell-state assignments". The analysis actually used the md2 re-derived labels, not Lu's published `cell_state`. That wording is worth fixing; this check shows the numbers survive either way. No manuscript file was edited.
- Disclosure (repeated from the pre-registration): Lu's cell_state × sample counts were tabulated before the rule was written.

## Files

- src/relabel_check.py (`step0`, `prereg`, `cluster`, `agree`, `stats`, `findings`)
- results/relabel_check/: PREREG.flag; lu_meta_*.csv, lu_slots_*.json, lu_meta_columns.csv (Step 0); validation_GM00731.json, joint_cells.csv, joint_cluster_labels.csv, joint_ams.npz, joint_meta.json (Step 1); agree_*.csv, agreement_summary.csv, d0_fibroblast_membership.json (Step 2); harness_check.csv, headline_old_new.csv, headline.json, partA_states_*.csv (Step 3); run_log.txt and console logs.
