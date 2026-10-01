# STATE

## current stage
L1 done. L2 STOP. L3 skipped. L4 written (`FINDINGS_LOWDIM.md`).

## L2 STOP
stop_L2=True. winner=NONE. Seed 20260914. Dataset 4442d412-91cb-4261-acca-8adf5fa04c11.
n_pass_angle=0. n_pass_transfer=0. n_pass_both=0.
Bar: bootstrap ≤ 35° AND at least one of HBCC→MSSM / MSSM→HBCC R² > 0.

| config | boot° | within-site R² (null) | H→M | M→H |
|---|---|---|---|---|
| gene_ridge | 45.1 | +0.231 | −0.213 | −0.499 |
| A_pca_k050 | 55.5 | +0.191 (−0.147) | −0.255 | −0.603 |
| A_pca_k150 | 44.6 | −0.147 (−0.244) | −0.259 | −0.525 |
| B_curated_k050 | 48.5 | −0.119 (−0.177) | −0.802 | −0.989 |
| C_wgcna_k150 | 85.0 | −0.465 (−0.830) | −0.450 | −0.619 |

Closest bootstrap: PCA k=150 at 44.6° (Δ −0.5° vs gene, not material).
All transfers negative. Young vs old within-bank still ~72–87°.

## decisions
- STOP L2: direction not identifiable at k ≪ n. Do not add variants.
- L3 not run (no TARGET matrix).
- Variant C: Ward on 2000 HVGs (average-linkage chained to singletons).
- Winner rule: L2a and L2c together, not R².
- Comparison table rebuilt from L2a–d source CSVs.

## remains
None. Negative result stands. Changing it needs more donors, more sites, or a different tissue.
