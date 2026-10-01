# PROGRESS_VERIFY

**Status:** Completed. No stop key fired. Same status as FINDINGS_VERIFY.md.

409 claim rows. match 274; rounding 5; mismatch 0; untraceable 36; not_a_claim 94.

## What must be fixed before posting

1. Replace `[GITHUB REPOSITORY URL]` and `[ZENODO DOI]` in the data-availability paragraph. Both brackets are still in the draft.
2. Change the two printed upper bounds `+0.413` to the value stored in `results/paper_figs/cxg_true_donor_results.csv` (`0.41248891146806693`, which prints as `+0.412` at three decimals). One is the Table 2 cell for the 65-donor MD contrast. The other is the same interval in the results paragraph.
3. The discussion says the partially reprogrammed cells end "twice as far" from the young target as they started. `results/same/stats.csv` stores `delta` `130.89356478446803` and `v_norm` `130.68843692001553` on the forward row, and no field stores a ratio. Decide the wording from those two fields.
4. The abstract says "across four transfers" while Table 1 prints five Spearman rows, and no results file stores a transfer count of 4. The range `0.43 to 0.56` itself matches the lowest and highest stored transfer rhos (`0.43338041757870716` and `0.5611771722875578`).
5. "Two thirds" is used for a stored progress of `0.674440748216413`. The same paragraph also prints `67%` and `0.674`, which match that field. The words do not.
6. GSE226189 is a scored cohort and has no paper in the reference list. The 1,836-factor sentence does not name Southard et al., which is only in the reference list. Per-cohort sample tables are still deferred to the repository.
7. These methods numbers are not in a file under `results/`: the alpha grid other than `0.01`, TMM trims `0.30` and `0.05`, the geometric-mean scale `1`, the standard-deviation floor `1e-12`, the p-value denominator `201`, and the bootstrap percentiles `2.5` and `97.5`. The stored `alpha` is `0.01`, `n_perm` and `n_boot` are `200`, and the gene filter and `prior.count=2` are stored. "Three readings were withdrawn" has no stored count. "Twenty thousand genes" has no stored count of 20000; the frozen ruler has `n_genes` 23485.
