# FALSIFICATION.md — pre-registered criteria for Phase 1d

**Status: written 2026-09-14 at the Phase 1b halt, BEFORE any 1c or 1d model was fitted. Not yet evaluated.**
This file is a commitment for whichever run reaches 1d. It must not be revised after 1d numbers exist; a continuation run must
re-commit it unchanged (or supersede it in a clearly dated new file *before* running 1d).

## The claim under test
"Chronological age and cell-type identity are separable axes in PBMC transcriptomic space": a gene set with high age variance and low
identity variance (A-genes) carries age information that a gene set with high identity variance and low age variance (I-genes) does not,
and vice-versa.

## Models (fixed)
- **Age model**: ridge regression (α by inner CV) on pseudobulk log₂-CPM of the gene set, predicting donor age; trained within cell type.
  Outer CV: 5 folds grouped by **10x pool** (donors nest in pools → this is donor-held-out *and* batch-held-out; both asserted in code).
  Metrics: R² and MAE on held-out donors, per cell type; plus a pooled cross-cell-type model.
- **Identity model**: multinomial logistic regression on the gene set, predicting cell type, same pool-grouped CV. Metrics: accuracy,
  macro-F1, confusion matrix; chance = 1/16 ≈ 0.0625 (16 types) or the majority-class rate (CD4_TCM ≈ 0.10), whichever is larger.
- Gene selection (1a thresholds, 1b filters) is **re-done inside each outer training fold** to avoid selection leakage; the reported
  performance is nested-CV performance.

## Negative controls (all required)
1. Shuffle donor ages (at donor level, within pool), retrain age model.
2. Shuffle cell-type labels (at pseudobulk level), retrain identity model.
3. Cross test: age model on I-genes; identity model on A-genes. The I-gene set used for the age cross-test is (a) the full I set and
   (b) a size-matched random subset of I-genes (same n as A-genes; 20 draws, median reported) so that set size cannot explain the result.

## Pre-declared thresholds
Let R²_A = held-out age R² using A-genes (median over cell types), R²_I = the same using size-matched I-genes.
Let Acc_I = identity accuracy using I-genes, Acc_A = using A-genes.

| control | expected if the pipeline is sound | falsifies soundness if |
|---|---|---|
| shuffled ages | R² ≤ 0.05 (and |R²| within ±0.05 of 0) | R² > 0.10 |
| shuffled cell types | accuracy ≤ 1.5 × chance | accuracy > 2 × chance |

**Separability is falsified (STOP 1d) if either:**
- **R²_I ≥ 0.8 × R²_A**, or **R²_A − R²_I < 0.10** (age is predicted about as well from identity genes as from age genes), or
- **Acc_A ≥ 0.9 × Acc_I** (identity is classified about as well from age genes as from identity genes).

Separability is *supported* only if R²_I ≤ 0.5 × R²_A **and** R²_A − R²_I ≥ 0.15 **and** Acc_A ≤ 0.7 × Acc_I, with shuffles passing.
Results between the two bands are reported as **inconclusive**, not as support.

## Anticipated failure modes I commit to reporting as failures, not fixing
- Age R² on I-genes is high because many I-genes have small but coherent age slopes (age acting *through* identity programmes,
  cf. LRRN3, BCL2, KLF6 in FINDINGS §4). This would be a real falsification, not a gene-set bug.
- Identity accuracy on A-genes is high because pseudobulk profiles of 100+ genes almost always encode cell type. If so, the honest
  conclusion is that "low identity variance per gene" does not yield "no identity information per set" — and the benchmark's construct is
  flawed at the gene-set level.
- The pooled cross-cell-type age model outperforming within-type models: report as evidence that the age signal is shared, not "generalizing".

## Things that would NOT count as rescuing separability
- Re-tuning t_age / k_ct / t_ct / k_age after seeing 1d.
- Switching CV grouping from pool to donor.
- Removing individual genes from I or A after seeing which ones drive the cross-test.
- Replacing OneK1K with AIDA to get a better number (AIDA is for replication of whatever OneK1K shows, good or bad).
