# FINDINGS_VERIFY_V9

Re-run of the VERIFY audit in `results/verify/PREREG.flag` against the current
manuscript. The protocol is unchanged: every quantitative claim in the draft is
traced to the file under `results/` that produced it, and the stored number is
read as it sits on disk. Nothing was recomputed, no analysis was rerun, and
neither the preprint nor `results/verify/claims.csv` was edited.

**Status:** completed. No stop key fired.

## Draft identity

Audited file: `paper/rejuvenation_readouts_preprint_9.docx`

- size: 1,297,833 bytes
- md5: `FC31277DAF8EFABBDB715E4E81EA1BFC`
- modification time: 2026-10-01 09:32:55
- extracted body: 159 paragraphs and 4 tables

The superseded draft audited by `claims.csv` is kept alongside it:

- `paper/rejuvenation_readouts_preprint_3.docx` — 1,100,494 bytes, md5
  `004445ECEA7BD41EC19505B0EAE23D6F`, 119 paragraphs and 3 tables
- `paper/PREPRINT_DRAFT.md` — 25,348 bytes, not audited

## Figure and table numbering changed between v3 and v9

This matters when reading the two audits side by side: a claim filed under
`P51` or `T2R3` in `claims.csv` is not the same sentence as the row with that
location in `claims_v9.csv`. Every figure except Table 3 moved.

| v9 | v3 | subject |
| --- | --- | --- |
| Figure 1 | — (new) | the published MD score falls while the cells move away from young cells |
| Figure 2 | — (new) | the result does not depend on which genes are read |
| Figure 3 | Figure 1 | transfer of the frozen fibroblast ruler |
| Figure 4 | Figure 3 | the frozen ruler against two curated instruments |
| Figure 5 | Figure 2 (panel B) | a better age score makes the same mistake |
| Supplementary Figure S1 | Figure 4 | reprogramming states against a GTEx-defined youth axis |
| Supplementary Table S1 | — (new) | the GTEx-anchored test in three gene spaces |
| — | Figure 5 | the same-platform test, now folded into Figures 1B and 2 |

Tables 1, 2 and 3 keep their numbers, but Table 1 gained the 65-donor atlas row
and Table 2 is now reported on the corrected donor unit throughout.

The paragraph indices also shift: v9 has 40 more paragraphs than v3, so the
`P<n>` locations in the two CSVs are not comparable.

## Stage 1 — extraction

Claim rows in `results/verify/claims_v9.csv`: **546**.

Every numeric token in the body text, tables, figure captions, Methods and
Supplementary Information was enumerated: correlations, intervals, distances,
counts, sample sizes, fractions, cosines, p-values, thresholds, accessions and
software versions.

## Stage 2 — verdicts

| verdict | count |
| --- | --- |
| verified | 346 |
| mismatch | 1 |
| no-source | 55 |
| not_a_claim | 144 |

Provenance for all 346 verified claims was confirmed by opening the stored file
and reading the named field, not by value coincidence.

The bulk of the manuscript's numbers now come from `results/newstory/`
(`partA_md_by_state.csv`, `partA_gap_drop.csv`, `partB_same_distances.csv`,
`partB_mixture.csv`, `partB_gtex_distances.csv`, `partB_random_sets_counts.csv`,
`fig_N5_ruler_same_mistake_data.csv`) and `results/genespace/`
(`same_stats.csv`, `gtex_stats.csv`, `gtex_c3.csv`, `readout_cos.csv`,
`spaces.json`). Tables 1 and 2 come from `results/fibro/stage1_transfer.csv`,
`results/md3/t2_instruments.csv`, `results/md3/t2_delta_rho.csv` and
`results/paper_figs/cxg_true_donor_results.csv`. The relabelling checks come
from `results/relabel_check/headline.json` and `agreement_summary.csv`.

### Mismatch

One, and it is a rounding error in the last printed digit.

**`P21` — `0.514`.** Stored `0.5134646010253517` in
`results/newstory/partA_md_by_state.csv`, field `joint_ci_lo` on row
`GM00731 / Fibroblast_d0`. At three decimals that value prints as **0.513**,
not 0.514. The same number is stored again as `old/a/…joint_ci_lo` in
`results/relabel_check/headline.json`.

> We scored MD per cell with a reimplementation of the published algorithm
> (Methods), scoring both donors against the same control genes so that their
> scores are comparable. The score behaves as reported. In the aged donor it
> falls from 0.515 (95% CI 0.514 to 0.517) in day-0 fibroblasts to 0.158 (0.153
> to 0.162) in partially reprogrammed cells, and continues downward through
> early pluripotency (0.147) to pluripotency (0.022).

Every other bound in that sentence is correct: `0.515` is `0.5151367570932331`,
`0.517` is `0.5169498495629744`, `0.158` is `0.15752579308739909`, `0.153` is
`0.1531354436736409`, `0.162` is `0.1617369394380799`, `0.147` is
`0.1467058603483695` and `0.022` is `0.021948962428081423`.

### No-source

Fifty-five rows, in five groups. None is a wrong number; each is a figure the
repository does not store.

**1. The interval level, `95`, in "95% CI" — 13 rows**
(`P21` ×2, `P22`, `P23`, `P27`, `P44`, `P47`, `P54`, `P56`, `P61`, `P159`,
`T2R1` ×2.) The draft states the level; the result files store `n_boot = 200`
without recording the percentile. Carried over unchanged from the v3 audit.

**2. The alpha grid and normalisation constants in Methods — 14 rows**
`P91`: every grid value except the selected `0.01` — `0.03`, `0.1`, `0.3`, `3`,
`10`, `30`, `100`, `300`, `1000`, `3000`, `10000`. `P92`: the TMM log-ratio trim
`0.30`, the absolute-expression trim `0.05`, and the standard-deviation floor
`1e-12`. `results/fibro/freeze.json` stores the chosen `alpha = 0.01`, the gene
filter and `prior.count = 2`, but not the grid it was chosen from or the trims.
Carried over unchanged from the v3 audit.

**3. Stated thresholds and parameters with no stored counterpart — 12 rows**
`P94` `30` (the minimum-cell threshold for the sensitivity analysis); `P102`
`3,000`, `50` and `0.8` (variable genes, principal components and Louvain
resolution) and `16` (the count of published sub-states in the GSE297234
deposit); `P104` `10,000` (the AddModuleScore scale factor); `P107` `201` (the
p-value denominator); `P108` `2.5` and `97.5` (the bootstrap percentiles);
`P111`, `P148` and `P154` `0.5` (the pre-registered pass bar for the five-fold
control); `P65` `50` (PCA components).

**4. Donor ages — 11 rows**
`96`, `22` and the derived `74` appear in `P7`, `P20`, `P21`, `P27`, `P65`,
`P68` and `P153`. These are GSE297234 sample metadata. No file under `results/`
stores a donor age for either cell line, so the ages themselves cannot be
checked against repository output, although every quantity computed from them
can be and is.

**5. Two counts taken from outside `results/` — 2 rows**
`P30` `3,538`, the size of the Fleischer et al. age list before mapping; the
repository stores the mapped counts (`1,312` of `1,533` age-up and `1,777` of
`2,007` age-down in `results/genespace/spaces.json`) but not the source-list
total. `P73` `1,836`, the transcription-factor count in the Southard et al.
screen; this is recorded in `SURVEY_TF.md` from BioProject PRJNA1108254 and from
the guide HDF5, but no file under `results/` holds it.

## Stage 3 — placeholders

Two bracketed placeholders are still in the draft, both in `P119`, the data and
code availability paragraph:

> Analysis code, the frozen instruments, the pre-registration documents with
> their timestamp flags, and every result table used in this paper are available
> at **[GITHUB REPOSITORY URL]** and archived at **[ZENODO DOI]**.

These are the same two placeholders the v3 audit reported; they remain
unresolved.

`P117` and `P119` both promise that "the repository README lists the entry-point
script that reproduces each table and figure". That README now exists at the
repository root.

## What changed since the v3 audit

- **The `+0.413` rounding error is fixed.** The v3 audit flagged two printed
  upper bounds of `+0.413` against a stored `0.41248891146806693`. v9 prints
  `+0.412` in both places, Table 2 row 2 and `P56`. Resolved.
- **The "two thirds" wording is gone.** v3 described a stored progress of
  `0.674440748216413` as "two thirds" in three places. v9 reports the MD-gene
  figure as `0.790` and `79%` and the all-gene figure as `67%`, each matching
  its stored field. Resolved.
- **The "four transfers" count is gone.** v3's abstract said "across four
  transfers" while Table 1 printed five rows. v9's abstract does not count
  transfers, and Table 1 now has four rows that match the four stored
  transfer statistics. Resolved.
- **"Twice as far" is now sourced.** v3 asserted the cells end "twice as far"
  from the young target with no stored ratio. v9 cites `2.2` times in the MD
  genes and `2.0` in all genes, which are the stored `final_over_start` fields
  in `results/newstory/partB_same_distances.csv`. Resolved.
- **Still open:** the two bracketed placeholders, the alpha grid and TMM trims,
  the bootstrap percentiles and the p-value denominator, and the donor ages.

## Caveats on this audit

- Thirty-one files under `results/` larger than 2 MB were not indexed by value;
  they are per-cell, per-gene and per-draw tables (cell labels, per-gene
  correlations, full perturbation tables) and are not the source of any number
  in the draft. Where a claim's source was one of these, it was opened directly.
- Software versions in `P117` were checked against
  `results/paper_figs/environment_versions.txt`, which is a text file and
  therefore outside the value index; it was read in full.
- `not_a_claim` covers the reference list, the author and back matter,
  cross-references to the paper's own figures and tables, accession and UUID
  fragments, timepoint labels in the GSE297234 design, and constants inside
  stated formulae.

## Note — mismatch corrected in preprint_10

The single mismatch above, **`P21`** (printed `0.514` against stored
`0.5134646010253517`), is corrected in
`paper/rejuvenation_readouts_preprint_10.docx`, which prints the lower bound as
`0.513`: "falls from 0.515 (95% CI 0.513 to 0.517) in day-0 fibroblasts". A
paragraph-by-paragraph text comparison of preprint_9 and preprint_10 finds no
other difference. preprint_10 is 1,297,833 bytes, md5
`4CBE1CF551E20B8BF168D46CC793C127`. The verdict counts above describe
preprint_9 and were not re-run.

## preprint_11 changes

Added 2026-10-02. `paper/rejuvenation_readouts_preprint_11.docx` (1,298,789 bytes, md5
`65F57FAD709672ABEAFB28344E802029`) is preprint_10 with Methods paragraphs **P89** and **P94**
changed and nothing else. The paragraph-by-paragraph, cell-by-cell and image-hash comparison is
in `results/verify/PREPRINT_11_CHANGES.md`. Every quantitative claim in the two changed
paragraphs was re-checked against the stored files, including the sentences that did not
change. Nothing was recomputed. The verdict counts above still describe preprint_9.

**Result: 14 claims checked. 13 were verified against stored outputs and 1 against code (the 30-cell threshold, a constant in `src/paper_figs_cxg_donor.py`). 0 mismatch.**

### P89

| claim | stored source | stored value | verdict |
| --- | --- | --- | --- |
| "Nine readings across the project are marked this way" | `results/verify/WITHDRAWN_READINGS.md`, sections A and B (A1–A2, B1–B7); each marking re-opened at its cited place | 2 withdrawn (A1 `FINDINGS_SOUTH4.md:1`; A2 `FINDINGS_MD4.md:482`) + 7 superseded (B1 `FINDINGS_SOUTH.md:1`, B2 `FINDINGS_SOUTH2.md:1`, B3 `FINDINGS_SOUTH3.md:1`, B4 `FINDINGS_SOUTH5.md:1`, B5 `FINDINGS_MD4.md:366`, B6 `FINDINGS_POSCTRL.md:3`, B7 `FINDINGS_MD3.md:1659`) = 9. Section C (withdrawn in text only) is empty | verified |
| "three of which concern analyses reported here" | `WITHDRAWN_READINGS.md`, table "Which of these bear on claims made in preprint_10" (same Results, tables and figures in preprint_11) | "Yes" rows: A2, B5, B7 = 3 | verified |
| "marked withdrawn or superseded" (for the atlas Δρ, B7) | appended note in `FINDINGS_MD3.md`, `FINDINGS_MD4.md`, `FINDINGS_MD5.md`, `FINDINGS_PLANE.md`, `PROGRESS_MD4.md`, `PROGRESS_MD5.md`, their `paper_package/` copies, and `paper_package/CONTEXT.md` | present in all 13 files | verified |
| "including a transfer that failed its gate" (unchanged) | `FINDINGS_EXTERNAL.md:3`, `:282` | SEA-AD transfer: "≥65-only primary fails" | verified |

P89 has no other number. "pre-registered in the same way, with any deviations recorded
alongside them" (B2) is not quantitative. Its support is the per-analysis PREREG flags, which
the bundle's `RESULTS_SUMMARY.md` section 3 lists, and the recorded deviations, for example
the SOUTH6 full run launched at the user's instruction (`results/south6/manifest.json`).

### P94

| claim | stored source | stored value | verdict |
| --- | --- | --- | --- |
| "giving 65 donors" (unchanged) | `results/paper_figs/cxg_true_donor_results.csv`, `analysis=true_donor` | `n_units` = 65, `n_pseudobulks` = 93 | verified |
| "64 after dropping pseudobulks with fewer than 30 cells" (unchanged) | same file, `analysis=true_donor_pseudobulk_ge30cells` | `n_units` = 64 (81 pseudobulks) | verified |
| "fewer than 30 cells" (unchanged) | `src/paper_figs_cxg_donor.py:26` | `MIN_CELLS = 30` (code constant; claims_v9 recorded it as no-source) | verified in code |
| "which gave the same conclusion" (unchanged) | same file, ge30 rows vs `true_donor` rows | ρ +0.437 (p 0.005) vs +0.433 (p 0.005); ruler − MD +0.177 [−0.105, +0.405] vs +0.171 [−0.118, +0.412], both include zero; ruler − age-up/down +0.382 [+0.080, +0.639] vs +0.372 [+0.094, +0.595], both exclude zero | verified |
| "Seeds and resampling settings were unchanged" (unchanged) | same file, all rows; `results/md3/t2_delta_rho.csv` | `seed_perm` 20260914, `seed_boot` 20260918, `n_perm` 200, `n_boot` 200 in every row; md3 `seed` 20260918, `n_boot` 200 | verified |
| "This grouping was not pre-registered" | `results/paper_figs/` | no `PREREG*` file in the directory or below; results written 2026-09-23, after MD4 quoted the MD3 Task 2 result on 2026-09-19 | verified |
| "each of the 93 pseudobulks as a unit" | `results/fibro2/ta_cxg_a19d1667_result.json`; `results/md3/PREREG_TASK2.flag:6` | `n_samples` = 93, `n_donors` = 93; the flag names the atlas "(10x, n=93)" | verified |
| "gave ρ = +0.455" | `results/fibro2/ta_cxg_a19d1667_result.json`, `rho` | 0.45512666633316906 → +0.455 (three decimals, explicit sign, as in P47 and Table 1) | verified |
| "it gives the lower value" | `cxg_true_donor_results.csv`, `true_donor` / `ruler_spearman` | 0.43338041757870716 → 0.433 < 0.455 | verified |
| "superseded by the per-person value in Table 2" | `cxg_true_donor_results.csv`, `true_donor` / `ruler_minus_MD_AddModuleScore`; preprint_11 Table 2 data row 1 | stored +0.171 [−0.118, +0.412]; printed "+0.171 [−0.118, +0.412]". 93-pseudobulk value (`t2_delta_rho.csv`, `cxg_ruler_minus_MD`): +0.202 [+0.006, +0.440] | verified |

The 93 and +0.455 asked for in the edit instructions matched the stored values and were written
unchanged.
