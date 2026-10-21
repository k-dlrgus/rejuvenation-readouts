# FINDINGS — DONOR

Reading: **DONOR-SPECIFIC**. No STOP fired. The gate passed (M_same = 0.934 >= 0.85; gate specificity 86 of 89 against a chance baseline of 4.45, binomial p = 1.26e-107).
DONOR-SPECIFIC fired because the CI upper bound 0.460 is below 0.5. Its specificity branch did **not** fire: specificity is above chance (22 of 89, baseline 4.45, p = 3.64e-10). Section 6 says what this does and does not mean.

## The pre-registration, verbatim

```
PREREG written 2026-10-21T22:17:04Z (local 2026-10-22 07:17:04 KST).
Written before any DONOR number exists. Reproduced verbatim at the top of
FINDINGS_DONOR.md. Not edited after this timestamp.

Repository state when written: commit 6815e95b850d46d80e37ad70de48b2a6c6b4524c
("SRSFINV: data inventory of GSE325735 (Plesa et al. SRSF1 screen)"), branch main.

TASK: DONOR — does a transgene push different aged donors' fibroblasts in the
same direction?

Written into this flag verbatim, as required:

  "This compares aged donors with each other. No reading licenses a claim
  about transfer between an aged and a young context, nor about the Southard
  screen."

--- 0. WHY THIS TASK EXISTS ---
TRANSFER (FINDINGS_TRANSFER.md, commit 7230271) asked whether a factor's
direction carries between cell contexts and stopped at its gate: at ~60 cells
per factor the Junevity (Sengstack) arm could not recover known agreement.
GSE325735 (inventory: SRSFINV, commit 6815e95) is bulk RNA-seq, so cell depth
is not the limitation, and 89 transgenes were applied to more than one aged
donor line. This task asks the same question where it can be measured.

This is a NARROWER question than TRANSFER's. All six perturbed lines are aged
55-79 (AG12657 55, AG06294A 65, AG09697 67, AG13369 68, AG12587 69,
AG16359 79). It tests agreement between aged donors, not between an aged and
a young context. Under the primary condition fixed in section 3 only three of
the six lines enter (55, 65, 79); see section 3.

What was looked at before this flag, disclosed: the SRSFINV inventory and its
tracked CSVs (results/srsfinv/samples_human_annotated.csv and others) were
cross-tabulated by line, batch, treatment, genotype, title Rep label and
description flag, to settle the metadata below. No expression value, no count
from the GSE325735 matrix and no ruler overlap was computed. The matrix was
not on disk.

Decisions taken by the author in chat before this flag, after the
replicate structure was found (quoted as answered):
  - gate at n = 2 replicates: "Uncorrected, bar 0.85" (section 6);
  - primary scope: "Dox1000 only" (section 3);
  - download of GSE325735_human_raw_counts.txt.gz: "Yes, download".

--- EXECUTION RULES (binding) ---
- src/donor_run.py runs standalone and refuses to compute if this flag is
  missing.
- On a STOP, write FINDINGS_DONOR.md with the STOP key and exit. Never compute
  after a STOP; never resume a stopped run with an adjustment.
- Writes only results/donor/, src/donor_run.py, FINDINGS_DONOR.md, and the
  raw download under data/raw/srsfinv/ (git-ignored). No earlier results
  directory is written. results/transfer/NOT_A_RESULT_gate_failed/ and
  results/reach/run.log are not opened.
- Imported, not modified: src/toward_run.py (tmm_logcpm_quiet, zscore_frozen,
  percentile_ci), src/fibro2_common.py (load_frozen_ruler), src/same_run.py
  (ci_covers_point). Stage 0's SOFT parser is written fresh in
  src/donor_run.py; src/srsfinv_run.py is not imported.
- Seed 20260914 for every random draw (bootstrap; any random split).
- n_boot = 10,000.
- Every failure is printed verbatim. No column is substituted and no row is
  repaired.

--- 1. DATA ---
(1a) Metadata: data/raw/srsfinv/GSE325735_family.soft.gz, md5 must equal the
  SRSFINV value 6257c37dbec6f49a5cff6f110ca641a4.
(1b) Counts: GSE325735_human_raw_counts.txt.gz from
  https://ftp.ncbi.nlm.nih.gov/geo/series/GSE325nnn/GSE325735/suppl/
  saved to data/raw/srsfinv/. Size must equal the 26,007,361 bytes SRSFINV
  recorded by HTTP HEAD; md5 recorded at run time (no stored md5 exists).
  Header: first field `Geneid`, then 807 sample columns whose set must equal
  the 807 human library names (`Library name:` in !Sample_description).
  Columns are matched to GSMs through the library name. Size or header
  mismatch: STOP prior_results_not_reproduced (it reproduces SRSFINV's
  stored file size and header check).

--- 2. METADATA RESOLUTIONS (SRSFINV ambiguities) ---
Every resolution below is a rule applied to the structured record. Where an
ambiguity cannot be settled, the affected samples are excluded. No ambiguity
is resolved by the reading that keeps more data.

A1 Series title. Not depended on.

A2 Donor identity. The unit is the `cell line` value. The six perturbed lines
  have six distinct Coriell-style IDs and six distinct described ages. That
  they are six different people is not verified from the record (A2), and sex
  is not stated. Where this task says "donor" it means "donor line". Not
  settled further; nothing depends on more than "distinct line ID".

A3 NTg. Not defined in the record; title and `genotype` disagree for 22 of
  the 30 NTg-titled GSMs. Every GSM whose `genotype` is `WT` or `NTg`, or
  whose title names an NTg transgene, is EXCLUDED, as transgene and as
  control. Reason: its identity cannot be settled from the record, and a
  dox-matched, consistently labelled control exists for every primary cell
  without it.

A4 What the control is. The control for transgene t in line L is the set of
  `genotype: mTagBFP2` GSMs with the same `cell line`, the same `batch`
  value, the same `treatment` string, and (where the description carries
  one) the same `Batch <n>` sub-label. Reason: mTagBFP2 is a Dox-inducible
  reporter on the same induction, so the displacement isolates which gene is
  induced rather than whether dox and a construct are present. The record
  never calls mTagBFP2 a control; that is this flag's reading and it is
  disclosed. Not used as control: WT/`NT` lines, NTg (A3), rapamycin,
  uninduced (`NT`) transgene GSMs, any GSM from another line or batch.

A5 Schedules. `Dox_1000ng/mL, 3days` and `Dox_1000ng/mL, 3d-ON` are NOT
  pooled. Each exact `treatment` string is its own condition. The primary
  uses `Dox_1000ng/mL, 3days` only (section 3). Every d-ON / d-OFF schedule
  (3d-ON included) is outside this task; those schedules carry only SRSF1,
  OSKM and mTagBFP2.

A6 Rapamycin dose. Not depended on; no rapamycin GSM is used.

A7 Title versus genotype. The transgene identity is the `genotype` field,
  and a GSM is used only if its title token ("Dox-inducible X transgene")
  equals its `genotype`. GSMs where they differ are EXCLUDED: the 2
  OSKM-titled GSMs with `genotype: SOX2` and the 2 TP53 (DN)-titled GSMs with
  `genotype: TP53` (and the NTg-titled WT GSMs, already excluded under A3).

A8 Cell type. Only `cell type: Fibroblast`.

A9 Batch. The batch key is the `batch` field (LIB...), refined by the
  description's `Batch <n>` where present. Neither is defined in the record.
  No primary GSM carries a `Batch <n>` label, so this does not affect the
  primary.

A10 Depth. Counts are used as deposited. Each GSM's library size (sum of its
  ruler-aligned counts) is reported as a diagnostic.

A11 Passage. The 2 GSMs whose title P<n> disagrees with the described
  passage (GSM9611955, GSM9611956; both `redo`, batch LIB053449) are
  EXCLUDED. Passage otherwise differs between lines and, within a cell,
  between transgene and control GSMs in some cases; it is recorded per cell,
  never used to exclude.

A12 Molecule. Only `polyA RNA` GSMs. Total-RNA and polyA samples are never
  compared.

A13 Gene overlap. Counted at run time (section 4).

Further rules the task depends on, fixed now:
R1 Description flags (`High GC`, `<10M reads`). Kept: they are QC notes on
  deposited rows and do not make the sample's identity ambiguous. Every
  per-pair row carries a column saying whether a flagged GSM entered it.
R2 One batch cell per (t, L, condition). If t's GSMs in line L under one
  condition sit in more than one batch cell, the batch cell with the most t
  GSMs is used; a tie excludes (t, L, condition). Controls come from that
  same batch cell only.
R3 Eligible cell. A (t, L, condition) cell is eligible only if, after the
  rules above, it has at least 2 transgene GSMs and at least 2 control GSMs,
  so that disjoint halves exist. Ineligible cells are listed with the reason.
R4 Candidate transgene. A transgene is a candidate if it has eligible cells in
  at least 2 lines under the primary condition. Transgene = any `genotype`
  value other than WT, NTg, mTagBFP2, GFP.

--- 3. PRIMARY CONDITION AND SAMPLE SET ---
Primary condition: `treatment: Dox_1000ng/mL, 3days`. Cross-line pairs are
formed only within this condition. Under it, transgene GSMs exist in three of
the six lines only: AG12657 (55), AG06294A (65), AG16359 (79). The other
three lines (67, 68, 69) carry transgenes only at `Dox_2000ng/mL, 3days`
(8 transgenes in all six lines, OSKM in 67/68/69), which enters as
diagnostic D5 and never as the primary.

Panel for normalisation: every GSM (transgene or control) in an eligible
cell (R3) of the primary condition or of the D5 condition
`Dox_2000ng/mL, 3days`. One joint TMM over that panel. Platform and batch of
every panel GSM are reported; the expectation from the metadata, verified at
run time, is that every eligible primary cell lies in batch LIB052667
(GPL24676, polyA).

--- 4. INSTRUMENTS: READ, NEVER REFIT ---
(4a) Frozen ruler results/fibro/frozen_ruler_ridge_raw.npz via
  fibro2_common.load_frozen_ruler(); keys `ensembl`, `mu`, `sd` (n = 23,485).
  `w` is not used.
(4b) Ruler alignment: a counts-file row maps to ruler gene j iff its Ensembl
  ID with the version suffix stripped (text before the first '.') equals
  ruler `ensembl[j]`. If several file rows map to one ruler gene, the first in
  file order is used (as seng_run.ruler_column_index does); how many such
  duplicates exist, and whether the unused copies carry any nonzero count, is
  reported. Reported: number of file rows, number of ruler genes mapped, number
  unmapped. Unmapped ruler genes get count 0 and therefore z = 0.
(4c) Transform, in this order: counts aligned to the 23,485 ruler columns ->
  joint TMM over the panel (toward_run.tmm_logcpm_quiet: TMM, edgeR log2-CPM,
  prior.count 2) -> frozen z (toward_run.zscore_frozen against the frozen GTEx
  mu/sd; columns all-zero across the panel set to 0) -> per-GSM z vectors.
  z(group) = mean of its GSMs' z vectors. The difference is formed in the
  23,485-column z space, and projection is last.
(4d) Space: results/south7/basis.npz, `L_mid_full[:, :20]`, the PCA loadings
  of the 308 GTEx donors aged 40-59. Orthonormality checked
  (max|L'L - I| <= 1e-8); failure is STOP prior_results_not_reproduced. No
  GSE325735 sample enters the basis, the ruler or mu/sd.

--- 5. THE MEASUREMENT (primary condition) ---
  d(t, L) = z(t's GSMs in its batch cell in L) - z(the cell's mTagBFP2 GSMs),
  projected onto the 20-component basis. All cosines are in that space unless
  a diagnostic says otherwise.

RELIABILITY FLOOR. Replicates allow a disjoint split in every eligible cell
  (R3), so the floor is built from transgene replicate splits, NOT from
  control-versus-control displacements. For a cell with exactly 2 transgene
  and 2 control GSMs: half A = {transgene GSM 1, control GSM x}, half B =
  {transgene GSM 2, control GSM y}, over both pairings (x, y) in {(1,2),(2,1)}
  (GSMs ordered by title Rep label, then accession). c_half(t, L) = mean of
  the two cosine(d_A, d_B). A cell with more GSMs uses 100 random disjoint
  halves (extra GSM to half A), numpy default_rng([20260914, line_index,
  transgene_index, split_index]).
  floor angle(t, L) = arccos(c_half), reported per cell.
  r(t, L) = 2 c_half / (1 + c_half) (Spearman-Brown to the full 2-replicate
  depth). If c_half <= 0, r is UNDEFINED, recorded as undefined, never 0.
  QUALIFYING CELL: r >= 0.3 (undefined fails).

CROSS-LINE AGREEMENT. For each candidate t and each unordered pair of lines
  {A, B} where both cells qualify:
    c_obs(t; A, B) = cosine(d(t, A), d(t, B)); angle = arccos(c_obs).
    c_corr(t; A, B) = c_obs / sqrt(r(t, A) r(t, B)). NOT clipped.
  c_corr(t) = mean of c_corr over t's qualifying pairs (1 or 3 pairs).
  QUALIFYING TRANSGENE: at least one qualifying pair. n_q = their number.

SPECIFICITY. For qualifying t, its directed pairs are (A -> B) and (B -> A)
  for each qualifying pair. s(t) = mean over directed pairs of
  cosine(d(t, A), d(t, B)). Comparison set G_t = transgenes g != t whose cell
  qualifies in every line that is a target B in t's directed pairs.
  m(t, g) = mean over t's directed pairs of cosine(d(t, A), d(g, B)).
  t is specific if s(t) > 95th percentile of {m(t, g) : g in G_t} (strict >,
  numpy linear interpolation). t is testable only if |G_t| >= 19; untestable
  transgenes are reported and left out of the denominator. Report k specific
  of n_test against the chance baseline 0.05 * n_test, with a one-sided exact
  binomial test (p0 = 0.05, H1 greater), LABELLED APPROXIMATE: the
  comparisons share displacement pools and are not independent. The count is
  never reported without the baseline beside it.

PRIMARY STATISTIC. M = median of c_corr(t) over qualifying transgenes, 95%
  bootstrap CI over transgenes (10,000 resamples with replacement, percentile
  2.5/97.5 as toward_run.percentile_ci, seed 20260914). arccos(M) in degrees
  (undefined if |M| > 1, stated as such). n_q reported. Non-finite c_corr
  values are reported and excluded, with their count. INTERVAL CHECK: if the
  CI does not contain M (same_run.ci_covers_point), the interval is INVALID
  and every CI-dependent branch is unreadable.

--- 6. QUALIFICATION GATE (runs before the primary reading) ---
Same-context control. For every candidate transgene t and every line L where
its primary cell is eligible, the two halves of section 5 (one transgene GSM
and one control GSM each, both pairings) are treated as two "donors". True
agreement is 1 by construction.

At 2 replicates per cell the halves are single GSMs, and their reliability
cannot be estimated without splitting one sample. The only estimate the data
allow is the split cosine itself, which would make the corrected value 1 by
construction. Decided by the author before this flag: the gate is
UNCORRECTED.
  c_same(t) = mean over t's eligible cells of c_half(t, L).
  M_same = median of c_same(t) over ALL candidate transgenes. No r >= 0.3
  filter is applied inside the gate, because that filter selects on the same
  split cosine.
  Gate specificity: s_same(t) = c_same(t). G'_t = candidates g != t eligible
  in every line where t is eligible. m_same(t, g) = mean over those lines, both
  pairings and both directions of cosine(d_A(t, L), d_B(g, L)) and
  cosine(d_B(t, L), d_A(g, L)) (halves from the same split, so no GSM is
  shared). Specific if s_same(t) > 95th percentile of {m_same(t, g)}; testable
  if |G'_t| >= 19. Binomial as in section 5, labelled approximate.

GATE PASSES only if M_same >= 0.85 AND gate specificity is above chance at
binomial p < 0.01. It fails otherwise, including an empty candidate set.
If it fails: STOP unmeasurable. Reading = UNMEASURABLE. The primary statistic
is not computed.

This gate is the recovery check: at single-replicate depth, it measures
whether the statistic recovers agreement that is 1 by construction. Its result
is reported beside any null.

--- 7. READINGS: mutually exclusive, wording fixed verbatim ---
Evaluated in this order; the first whose conditions hold is the conclusion.
Placeholders in <> are filled from the run.

(R0) UNMEASURABLE — the gate failed, or n_q < 20.
  "UNMEASURABLE. [The same-context control did not recover agreement at
  single-replicate depth (M_same = <M_same>, bar 0.85; gate specificity <k'>
  of <n'> against a chance baseline of <0.05 n'>, binomial p = <p'>) / Only
  <n_q> transgenes had sufficient within-line reliability, fewer than the 20
  fixed in advance.] No statement is made about whether a transgene's
  direction carries between aged donor lines. This is a statement about the
  measurement, not about the biology. This compares aged donors with each
  other. No reading licenses a claim about transfer between an aged and a
  young context, nor about the Southard screen."

(R1) DONOR-SPECIFIC — CI upper bound of M below 0.5, OR specificity not above
  chance (binomial p >= 0.05). (If the interval is INVALID, only the
  specificity branch can fire.)
  "DONOR-SPECIFIC. A transgene's effect does not carry between aged donors
  even at bulk depth: the median reliability-corrected cosine between the same
  transgene's displacements in different donor lines is <M> (95% CI <lo> to
  <hi>, <arccos M> degrees) over <n_q> qualifying transgenes, and <k> of
  <n_test> transgenes were specific against a chance baseline of
  <0.05 n_test> (binomial p = <p>, approximate). A forward model cannot borrow
  a factor's direction across donors. The combination screen needs effects
  measured per donor and more donors than planned. The same-context control
  recovered M_same = <M_same> at single-replicate depth. Measured in three
  lines aged 55, 65 and 79, under one dox condition, in one library batch.
  This compares aged donors with each other. No reading licenses a claim
  about transfer between an aged and a young context, nor about the Southard
  screen."

(R2) CONSISTENT — CI lower bound of M at or above 0.7 AND specificity above
  chance at binomial p < 0.01.
  "CONSISTENT. For these <n_q> transgenes, the displacement in one aged donor
  line points the same way as in another: median reliability-corrected cosine
  <M> (95% CI <lo> to <hi>, <arccos M> degrees), with <k> of <n_test>
  transgenes specific against a chance baseline of <0.05 n_test> (binomial
  p = <p>, approximate). This holds for these transgenes in these three lines
  aged 55, 65 and 79 (of the six perturbed lines aged 55-79) only, under one
  dox condition in one library batch. It says nothing about transfer to a
  different cell type, to younger donors, or to untested factors. This
  compares aged donors with each other. No reading licenses a claim about
  transfer between an aged and a young context, nor about the Southard
  screen."

(R3) INTERMEDIATE — anything else, including an INVALID interval with
  specificity p < 0.05.
  "INTERMEDIATE. The measurement falls between the pre-registered thresholds:
  median reliability-corrected cosine <M> (95% CI <lo> to <hi>, <arccos M>
  degrees) over <n_q> qualifying transgenes, with <k> of <n_test> specific
  against a chance baseline of <0.05 n_test> (binomial p = <p>, approximate).
  This neither establishes that transgene directions carry between aged donor
  lines nor that they do not. This compares aged donors with each other. No
  reading licenses a claim about transfer between an aged and a young context,
  nor about the Southard screen."

--- 8. DIAGNOSTICS (never gated; each labelled "diagnostic, not a gated
result" where it appears) ---
(D1) Line mean removed. In each line, the mean d over that line's eligible
  primary cells is subtracted from every d(t, L), and each half-displacement
  has the mean of the same half over the line's eligible cells subtracted.
  The whole primary readout (r, qualification, c_corr, M with CI,
  specificity) is rerun. Disclosed where it appears: the line mean includes t
  itself, so this is not independent of the vectors it corrects.
(D2) Batch versus donor. If, as expected (section 3), all eligible primary
  cells share one batch (LIB052667), primary donor differences cannot track
  batch. Batch varies
  only across conditions: under Dox 2000, lines 55/65/79 are in LIB053907 and
  67/68/69 in LIB054839; the same line and transgene appear in LIB052667
  (Dox 1000) and LIB053907 (Dox 2000). For the transgenes present in both,
  c_obs and c_corr are tabulated by category: same line / different batch
  (dose differs too, so batch and dose are confounded); different line / same
  batch; different line / different batch. Medians per category; descriptive.
(D3) Agreement and replicate count. Eligible primary cells have 2 replicates
  by design; if replicate count does not vary, that is stated and the relation
  is reported as undefined. Reported instead: Spearman correlation of c_corr(t)
  with the mean library size of t's GSMs and with t's number of qualifying
  lines, with p-values.
(D4) Full gene space. c_obs for every qualifying pair in the 23,485-column z
  space, unprojected. TRANSFER's expectation was near-orthogonality, because
  noise spread over 23,485 dimensions dilutes shared signal. At bulk depth the
  noise is smaller and the expectation is weaker; the angles are reported as
  they are.
(D5) Dox 2000 six-line set. The section 5 readout (reliability, c_obs, c_corr,
  per-transgene table, median) for the `Dox_2000ng/mL, 3days` transgenes, all
  line pairs. Fewer than 20 transgenes, so it could not be a reading; no
  specificity test (|G_t| < 19).
(D6) Norm bookkeeping: share of each d's squared z-space norm retained by the
  20-component projection.
(D7) Per-transgene and per-pair tables: lines, GSMs, controls, passages,
  flags, c_half, floor angle, r, c_obs, angle, c_corr, qualifying, specific,
  the 95th percentile of the comparison set, |G_t|.

--- 9. STAGE 0: REPRODUCE THE STORED NUMBERS FIRST (gate) ---
Before anything new, with this task's code from the family SOFT:
  839 GSMs; 807 human; 664 human fibroblast GSMs with a genotype other than
  WT/NTg/mTagBFP2/GFP; 96 such genotype values; 89 of them at more than one
  described donor age; and every row and column of
  results/srsfinv/donors.csv (19 cell lines). Also the SOFT md5 (1a), the
  counts-file size and header (1b), and the basis orthonormality (4d).
  Counts exact; strings exact. Any mismatch: STOP
  prior_results_not_reproduced, naming every mismatch with stored and new
  value.

--- 10. WHAT THIS DOES NOT COVER (flag and findings) ---
- Six donor lines, all aged 55-79, and the primary uses three of them (55,
  65, 79). No young donor and no cell type other than dermal fibroblast.
- Passage differs between lines (and between transgene and control GSMs
  within some cells).
- This says nothing about whether any transgene approaches young: no young
  target enters anywhere.
- The authors' own age predictor is not used anywhere.
- One dox condition (1000 ng/mL, 3 days) in one library batch. Agreement
  under other doses, schedules or batches is not measured by the primary.
- Two replicates per cell: the reliability correction rests on one
  single-replicate split per cell (two control pairings), so r is itself a
  noisy estimate and c_corr inherits that noise.
- Within a line every transgene shares the same 2 mTagBFP2 GSMs, so control
  noise is common to all displacements in that line. It cancels from no
  cross-line cosine but adds a shared component within a line (D1 removes the
  line mean).
- The 20-component basis was chosen in SOUTH7/SOUTH9 for its sensitivity on
  the Southard screen; it is applied here unchanged and is not chosen for
  this data.
- The specificity binomial is approximate (dependent comparisons).
- Direction only: nothing about effect magnitude or usefulness.

--- STOP KEYS ---
  prior_results_not_reproduced   Stage 0 mismatch (incl. file size/header,
                                 SOFT md5, basis orthonormality)
  unmeasurable                   the section 6 gate failed
  ambiguity_unresolved           a metadata ambiguity the task depends on
                                 cannot be settled from the record

--- WHAT DOES NOT COUNT ---
- Changing the control rule, the condition, the space, the thresholds 0.3 /
  0.5 / 0.7 / 0.85, the seeds or the readings after seeing a number.
- Clipping c_corr or dropping a transgene because its value is awkward.
- Refitting the ruler, mu/sd or the basis.
- Reading a reading whose conditions did not fire, or reporting a diagnostic
  as a conclusion.
- Pooling 3days with 3d-ON, using NTg or WT GSMs as controls, or comparing
  total-RNA with polyA GSMs.

--- OUTPUT ---
FINDINGS_DONOR.md: this flag verbatim; every deviation and every choice the
flag did not fix; the Stage 0 table; ruler mapping counts; exclusions and
eligibility; the gate table; the per-transgene table; specificity with its
baseline; M with CI and angle; the one reading that fired, verbatim; the
diagnostics, labelled; section 10 verbatim.
results/donor/: this flag, per-stage CSVs, summary JSONs, run log.
```

`results/donor/PREREG.flag` md5 `a5d3b22986b6dab9f5af955c392498f6`.

## 1. Deviations and choices the flag did not fix

- **The gate is uncorrected.** The task text gates on the corrected median. At 2 replicates per cell the gate's halves are single GSMs, and correcting them would divide the split cosine by itself (1 by construction). The author chose "Uncorrected, bar 0.85" in chat before the flag was written. The gate is therefore stricter than a corrected gate on the same halves, and it is computed without the r >= 0.3 filter (flag section 6).
- **The primary is one condition in three lines, not six.** The author chose "Dox1000 only" in chat before the flag. The task's fixed CONSISTENT wording said "these six lines aged 55-79"; the flag's wording says "these three lines aged 55, 65 and 79 (of the six perturbed lines aged 55-79)". Lines 67, 68 and 69 enter only diagnostic D5.
- **Readings carry added sentences.** Each branch keeps the elements the task required and adds the gate result, the line and batch scope, and the verbatim aged-donors sentence (flag section 7).
- **What the task did not specify, fixed in the flag:** the control (mTagBFP2 in the same line, batch, treatment and sub-batch); the batch-cell rule R2; eligibility R3; the reliability threshold r >= 0.3 (TRANSFER's); per-transgene c_corr as the mean over its qualifying line pairs; the specificity construction (comparison set G_t, |G_t| >= 19, strict 95th percentile); z(group) as the mean of per-GSM z vectors rather than a z of pooled counts; one joint TMM over every GSM in an eligible cell of either condition; duplicate Ensembl IDs resolved to the first row in file order; the D1-D7 constructions.
- **STOP keys.** The file-size, header, SOFT-md5 and basis-orthonormality checks are folded into `prior_results_not_reproduced` instead of separate keys. None fired.
- **Before the flag**, the SRSFINV tracked CSVs were cross-tabulated (line x batch x treatment x genotype x Rep x description flag) to settle the metadata and find the replicate structure. No expression value existed then. The flag file was created with a timestamp placeholder and revised three times before it was stamped (panel definition narrowed to eligible cells; a parenthetical claim about platforms removed; D2 made conditional on the batch check). The stamp went in at 2026-10-21T22:17:04Z, before the counts matrix was downloaded. It has not been edited since.
- **Download.** `GSE325735_human_raw_counts.txt.gz` was fetched with curl into `data/raw/srsfinv/` (SRSFINV's raw directory, git-ignored), after the flag. Size 26,007,361 bytes matches SRSFINV's HEAD; md5 8e3ecb46366730a5f6479147bf57895f. `fetch_counts()` (urllib, same URL) was added to the script after the run so the committed code can fetch it. It does not run when the file exists.
- **Imports.** `toward_run` (tmm_logcpm_quiet, zscore_frozen, percentile_ci), `fibro2_common` (load_frozen_ruler), `same_run` (ci_covers_point), none modified. Importing `toward_run` runs its module-level `mkdir` of `results/toward/`, which already exists; nothing is written there.
- **Order of work.** The findings writer (`findings()`, `md()`, `reading_text()`, `duplicate_ids()`, this list) and `fetch_counts()` were written after `read` had produced its numbers. They format stored outputs and compute no statistic. `duplicate_ids()` re-reads the raw matrix only to list which duplicated IDs carry counts. No analysis code was changed after any DONOR number existed, and no stage was rerun.
- **SUMO1 in line 65** has dox GSMs in two batch cells (LIB051945, LIB054839), 2 each. R2's tie rule excludes the cell. Neither batch has a line-65 mTagBFP2 at that treatment, so it would have been ineligible anyway. CYCS and UBE2I in line 65 (LIB053449) have one mTagBFP2 GSM, so they are ineligible under R3.
- **OSKM** enters the primary in lines 55 and 65 only. Its two line-79 Dox 1000 GSMs carry `genotype: SOX2` and are excluded under A7.
- **Files outside results/donor/.** `.gitignore` gained one entry, `results/donor/displacements_full.npy` (26.6 MB, regenerable), with a comment naming its source. `FINDINGS_DONOR.md` was regenerated more than once by the `findings` stage while its wording was revised; the numbers come from the single `read` run.
- **Seeds.** Every eligible cell is 2 + 2, so no random split was drawn. The only random draw is the bootstrap (default_rng(20260914), 10,000 resamples), used for the primary and D1.

## 2. Stage 0 — reproduction of SRSFINV (gate)

210 quantities compared, 0 mismatches. Counts and strings exact. The 190 per-cell comparisons of `results/srsfinv/donors.csv` (19 lines x 10 columns) are in `results/donor/stage0_reproduction.csv`; the other rows are below.

| quantity | stored | new | match |
|---|---|---|---|
| SOFT md5 | 6257c37dbec6f49a5cff6f110ca641a4 | 6257c37dbec6f49a5cff6f110ca641a4 | True |
| n GSMs | 839 | 839 | True |
| n GSMs (inventory_meta n_gsm) | 839 | 839 | True |
| n human GSMs | 807 | 807 | True |
| n human GSMs (inventory_meta n_human) | 807 | 807 | True |
| n human fibroblast GSMs (inventory_meta n_fibroblast) | 798 | 798 | True |
| n transgene GSMs (human fibroblast, genotype not WT/NTg/mTagBFP2/GFP) | 664 | 664 | True |
| n transgene genotype values | 96 | 96 | True |
| n transgenes at more than one described age | 89 | 89 | True |
| transgenes at more than one age, names (inventory_meta) | ACTG1;ACTL6A;AKT1;ATG5;CCNA2;CCNB1;CDC20;CDC34;CDK1;CDK2;CDK4;CDKN2A;CHEK1;CHEK2;CTNNB1;CXCR2;CYCS;DDX3X;EEF2;ESR1;FAU;GART;GSK3B;HDAC1;HDAC6;HES1;HNRNPA2B1;HSP90AB1;HSP90B1;HSPA5;HSPA8;HSPD1;ISG15;KAT7;KDM6A;LRRK2;MAPK14;MAPK7;MCM7;MDM2;NFKB1;NOTCH1;NPM1;OGT;OSKM;P4HB;PARK7;PIK3C2A;PIK3CA;PKM;PLK1;PPIA;PPP1CA;PPP1CB;PPP1CC;PPP2CA;PPP3CA;PPP4C;PSMC6;PSMD4;PTEN;RAD23B;RAD51;RNASEL;RPLP0;RPS27A;RRM1;RUVBL2;SIRT1;SLC2A10;SLC2A13;SOD1;SOD2;SRSF1;SUMO1;TP53 (DN);TXN;UBA52;UBB;UBE2A;UBE2I;UBE2N;UBE2S;UBQLN1;VCP;WDTC1;XPO1;YWHAB;YWHAQ | ACTG1;ACTL6A;AKT1;ATG5;CCNA2;CCNB1;CDC20;CDC34;CDK1;CDK2;CDK4;CDKN2A;CHEK1;CHEK2;CTNNB1;CXCR2;CYCS;DDX3X;EEF2;ESR1;FAU;GART;GSK3B;HDAC1;HDAC6;HES1;HNRNPA2B1;HSP90AB1;HSP90B1;HSPA5;HSPA8;HSPD1;ISG15;KAT7;KDM6A;LRRK2;MAPK14;MAPK7;MCM7;MDM2;NFKB1;NOTCH1;NPM1;OGT;OSKM;P4HB;PARK7;PIK3C2A;PIK3CA;PKM;PLK1;PPIA;PPP1CA;PPP1CB;PPP1CC;PPP2CA;PPP3CA;PPP4C;PSMC6;PSMD4;PTEN;RAD23B;RAD51;RNASEL;RPLP0;RPS27A;RRM1;RUVBL2;SIRT1;SLC2A10;SLC2A13;SOD1;SOD2;SRSF1;SUMO1;TP53 (DN);TXN;UBA52;UBB;UBE2A;UBE2I;UBE2N;UBE2S;UBQLN1;VCP;WDTC1;XPO1;YWHAB;YWHAQ | True |
| donors.csv n rows (cell lines) | 19 | 19 | True |
| donors.csv columns | cell type;n Fibroblast GSMs;n iPSC GSMs;described age(s);title M-code(s);described passage(s);title P token(s);sex;n GSMs;n batches | cell type;n Fibroblast GSMs;n iPSC GSMs;described age(s);title M-code(s);described passage(s);title P token(s);sex;n GSMs;n batches | True |
| counts file bytes | 26007361 | 26007361 | True |
| counts header first field | Geneid | Geneid | True |
| counts header sample columns (inventory_meta counts_header_cols) | 807 | 807 | True |
| counts header columns == human library names (set) | True | True | True |
| counts header duplicate columns | 0 | 0 | True |
| basis max|L'L - I| <= 1e-8 | True | True | True |
| basis rows (ruler genes) | 23485 | 23485 | True |
| basis n_mid | 308 | 308 | True |

Counts file md5 (recorded, no stored value to compare): `8e3ecb46366730a5f6479147bf57895f`, 26,007,361 bytes. Basis max|L'L - I| = 3.55e-15.

## 3. Ruler mapping, sample rules, cells

**Ruler mapping (SRSFINV A13).** The counts file has 58,721 rows (versioned Ensembl IDs). On unversioned IDs, **22,662 of the 23,485 ruler genes map**; 823 are unmapped and set to z = 0. 36,014 file IDs are not ruler genes. Only 5,755 file IDs match a ruler `gene_id` with its version, so matching on versions would have lost most of the ruler. Across the 584 panel GSMs, 954 ruler columns are all-zero and carry z = 0 (823 unmapped + 131 mapped but never counted). Median share of a panel library's counts on mapped ruler genes: 0.9983; median ruler library size 14,466,012. TMM factors 0.616 to 1.416.

Duplicate unversioned IDs: 45 file rows; 21 of them duplicate a ruler gene, and the first row in file order was used (flag 4b). Unused copies with any count in the panel: 2:

| unused_row | used_row | unused_panel_sum | used_panel_sum |
|---|---|---|---|
| ENSG00000214717.11_PAR_Y | ENSG00000214717.11 | 1 | 0 |
| ENSG00000002586.19_PAR_Y | ENSG00000002586.19 | 42 | 145 |

**Sample rules.** Of the 771 GSMs from the six perturbed lines, 597 pass section 2 and 174 are excluded (one row per GSM in `results/donor/excluded_gsms.csv`; a GSM can carry more than one reason). Counts by rule:

| rule | gsm_reasons |
|---|---|
| A11 | 2 |
| A12 | 4 |
| A3 | 34 |
| A5 | 152 |
| A7 | 4 |

GSMs excluded under A7 and A11:

| gsm | title | ch_genotype | exclusion |
|---|---|---|---|
| GSM9611856 | M79 line with Dox-inducible OSKM transgene, Dox_1000ng/mL, P13, Rep1 | SOX2 | A7 title token 'OSKM' != genotype 'SOX2' |
| GSM9611857 | M79 line with Dox-inducible OSKM transgene, Dox_1000ng/mL, P13, Rep2 | SOX2 | A7 title token 'OSKM' != genotype 'SOX2' |
| GSM9611918 | M79 line with Dox-inducible TP53 (DN) transgene, Dox_1000ng/mL, P13, Rep1 | TP53 | A7 title token 'TP53 (DN)' != genotype 'TP53' |
| GSM9611919 | M79 line with Dox-inducible TP53 (DN) transgene, Dox_1000ng/mL, P13, Rep2 | TP53 | A7 title token 'TP53 (DN)' != genotype 'TP53' |
| GSM9611955 | M65 line with Dox-inducible MDM2 transgene, Dox_1000ng/mL, P10, Rep1, redo | MDM2 | A11 title P10 != described passage 9 |
| GSM9611956 | M65 line with Dox-inducible UBA52 transgene, Dox_1000ng/mL, P10, Rep2, redo | UBA52 | A11 title P10 != described passage 9 |

**Cells.** 286 (transgene, line, condition) cells; 283 eligible (2 transgene + 2 mTagBFP2 GSMs in one batch cell; every eligible cell has exactly 2 + 2, so every split is the 2-pairing enumeration). Primary eligible cells by batch: {'LIB052667': 232}, platform {'GPL24676': 232}, all polyA. Ineligible cells:

| transgene | line | cond | bcells_available | n_t | n_c | reason |
|---|---|---|---|---|---|---|
| CYCS | AG06294A | Dox_1000ng/mL, 3days | LIB053449| (2) | 2 | 1 | R3 1 control GSM(s) |
| SUMO1 | AG06294A | Dox_1000ng/mL, 3days | LIB051945| (2); LIB054839| (2) | 0 | 0 | R2 tie between batch cells |
| UBE2I | AG06294A | Dox_1000ng/mL, 3days | LIB053449| (2) | 2 | 1 | R3 1 control GSM(s) |

Primary candidates (eligible in >= 2 lines): 89. Line sets: 55/65/79: 49, 55/79: 37, 55: 3, 55/65: 2, 79: 1, 65/79: 1, 65: 1. Single-line, not candidates: DICER1, MTOR, PIK3C3, SMARCA4, UBXN7.

Every line-65 primary cell uses the same two mTagBFP2 GSMs, one of which (GSM9611762) carries the description flag `High GC`; kept under rule R1. Flagged GSMs per cell are in `cells.csv`.

## 4. Qualification gate (same-context control, uncorrected)

| quantity | value | bar |
|---|---|---|
| candidate transgenes | 89 |  |
| M_same = median c_same (single-replicate halves) | 0.9340 | >= 0.85 |
| angle of M_same (degrees) | 20.9 |  |
| c_same interquartile range | 0.900 to 0.952 |  |
| specific / testable | 86 / 89 |  |
| chance baseline (5%) | 4.45 |  |
| binomial p (one-sided, approximate) | 1.26e-107 | < 0.01 |
| gate | PASS |  |

This is the recovery check this project requires beside any null: at single-replicate depth, inside one line and one batch, the statistic recovers agreement that is 1 by construction at 0.934 (20.9 degrees) and identifies 86 of 89 transgenes. Per-transgene gate table (`results/donor/gate_table.csv`):

| transgene | lines | c_same | angle_same | n_comp | p95_comp | specific |
|---|---|---|---|---|---|---|
| ACTG1 | 55/79 | 0.881 | 28.271 | 85 | 0.697 | True |
| ACTL6A | 55/65/79 | 0.935 | 20.820 | 48 | 0.863 | True |
| AKT1 | 55/79 | 0.914 | 23.980 | 85 | 0.576 | True |
| ATG5 | 55/65/79 | 0.952 | 17.904 | 48 | 0.812 | True |
| CCNA2 | 55/65/79 | 0.946 | 18.885 | 48 | 0.865 | True |
| CCNB1 | 55/65/79 | 0.928 | 21.860 | 48 | 0.863 | True |
| CDC20 | 55/65/79 | 0.904 | 25.306 | 48 | 0.818 | True |
| CDC34 | 55/79 | 0.953 | 17.718 | 85 | 0.801 | True |
| CDK1 | 55/65/79 | 0.943 | 19.507 | 48 | 0.863 | True |
| CDK2 | 55/79 | 0.957 | 16.781 | 85 | 0.795 | True |
| CDK4 | 55/79 | 0.944 | 19.214 | 85 | 0.826 | True |
| CDKN2A | 55/65/79 | 0.886 | 27.636 | 48 | 0.756 | True |
| CHEK1 | 55/79 | 0.702 | 45.395 | 85 | 0.584 | True |
| CHEK2 | 55/65/79 | 0.900 | 25.847 | 48 | 0.794 | True |
| CTNNB1 | 55/79 | 0.946 | 18.951 | 85 | 0.792 | True |
| CXCR2 | 55/65/79 | 0.826 | 34.319 | 48 | 0.796 | True |
| CYCS | 55/79 | 0.956 | 17.001 | 85 | 0.798 | True |
| DDX3X | 55/79 | 0.985 | 10.039 | 85 | 0.862 | True |
| EEF2 | 55/65/79 | 0.947 | 18.766 | 48 | 0.815 | True |
| ESR1 | 55/65/79 | 0.953 | 17.690 | 48 | 0.813 | True |
| FAU | 55/79 | 0.975 | 12.958 | 85 | 0.863 | True |
| GART | 55/79 | 0.957 | 16.833 | 85 | 0.850 | True |
| GSK3B | 55/65/79 | 0.962 | 15.889 | 48 | 0.837 | True |
| HDAC1 | 55/79 | 0.966 | 14.917 | 85 | 0.853 | True |
| HDAC6 | 55/79 | 0.963 | 15.676 | 85 | 0.827 | True |
| HES1 | 55/65/79 | 0.931 | 21.444 | 48 | 0.666 | True |
| HNRNPA2B1 | 55/79 | 0.897 | 26.276 | 85 | 0.699 | True |
| HSP90AB1 | 55/79 | 0.957 | 16.780 | 85 | 0.677 | True |
| HSP90B1 | 55/79 | 0.897 | 26.218 | 85 | 0.718 | True |
| HSPA5 | 55/79 | 0.966 | 15.020 | 85 | 0.809 | True |
| HSPA8 | 55/79 | 0.962 | 15.852 | 85 | 0.839 | True |
| HSPD1 | 55/79 | 0.956 | 17.052 | 85 | 0.749 | True |
| ISG15 | 55/65/79 | 0.786 | 38.227 | 48 | 0.749 | True |
| KAT7 | 55/65/79 | 0.906 | 25.059 | 48 | 0.783 | True |
| KDM6A | 55/79 | 0.958 | 16.626 | 85 | 0.718 | True |
| LRRK2 | 65/79 | 0.939 | 20.038 | 49 | 0.841 | True |
| MAPK14 | 55/79 | 0.967 | 14.702 | 85 | 0.830 | True |
| MAPK7 | 55/65/79 | 0.946 | 18.920 | 48 | 0.721 | True |
| MCM7 | 55/79 | 0.934 | 20.935 | 85 | 0.706 | True |
| MDM2 | 55/65/79 | 0.937 | 20.456 | 48 | 0.815 | True |
| NFKB1 | 55/65/79 | 0.934 | 20.882 | 48 | 0.756 | True |
| NOTCH1 | 55/79 | 0.983 | 10.464 | 85 | 0.690 | True |
| NPM1 | 55/79 | 0.944 | 19.241 | 85 | 0.695 | True |
| OGT | 55/65/79 | 0.946 | 18.869 | 48 | 0.792 | True |
| OSKM | 55/65 | 0.960 | 16.244 | 50 | 0.568 | True |
| P4HB | 55/65/79 | 0.900 | 25.777 | 48 | 0.759 | True |
| PARK7 | 55/65/79 | 0.941 | 19.748 | 48 | 0.827 | True |
| PIK3C2A | 55/65/79 | 0.921 | 22.870 | 48 | 0.797 | True |
| PIK3CA | 55/65/79 | 0.834 | 33.537 | 48 | 0.753 | True |
| PKM | 55/65/79 | 0.932 | 21.176 | 48 | 0.857 | True |
| PLK1 | 55/79 | 0.939 | 20.171 | 85 | 0.740 | True |
| PPIA | 55/65/79 | 0.943 | 19.462 | 48 | 0.879 | True |
| PPP1CA | 55/65/79 | 0.908 | 24.705 | 48 | 0.768 | True |
| PPP1CB | 55/65/79 | 0.763 | 40.291 | 48 | 0.788 | False |
| PPP1CC | 55/65/79 | 0.940 | 19.922 | 48 | 0.836 | True |
| PPP2CA | 55/79 | 0.963 | 15.549 | 85 | 0.854 | True |
| PPP3CA | 55/65/79 | 0.901 | 25.738 | 48 | 0.822 | True |
| PPP4C | 55/65/79 | 0.951 | 18.045 | 48 | 0.780 | True |
| PSMC6 | 55/65/79 | 0.918 | 23.398 | 48 | 0.781 | True |
| PSMD4 | 55/79 | 0.891 | 27.047 | 85 | 0.772 | True |
| PTEN | 55/65/79 | 0.949 | 18.466 | 48 | 0.828 | True |
| RAD23B | 55/65/79 | 0.938 | 20.262 | 48 | 0.852 | True |
| RAD51 | 55/79 | 0.970 | 14.070 | 85 | 0.660 | True |
| RNASEL | 55/79 | 0.906 | 25.100 | 85 | 0.559 | True |
| RPLP0 | 55/65/79 | 0.894 | 26.612 | 48 | 0.813 | True |
| RPS27A | 55/65/79 | 0.917 | 23.447 | 48 | 0.856 | True |
| RRM1 | 55/79 | 0.942 | 19.626 | 85 | 0.806 | True |
| RUVBL2 | 55/79 | 0.852 | 31.576 | 85 | 0.824 | True |
| SIRT1 | 55/65/79 | 0.932 | 21.282 | 48 | 0.802 | True |
| SLC2A10 | 55/65/79 | 0.834 | 33.437 | 48 | 0.734 | True |
| SLC2A13 | 55/65/79 | 0.938 | 20.313 | 48 | 0.672 | True |
| SOD1 | 55/65/79 | 0.928 | 21.914 | 48 | 0.869 | True |
| SOD2 | 55/65/79 | 0.932 | 21.255 | 48 | 0.867 | True |
| SRSF1 | 55/65/79 | 0.963 | 15.709 | 48 | 0.369 | True |
| SUMO1 | 55/79 | 0.933 | 21.031 | 85 | 0.722 | True |
| TP53 (DN) | 55/65 | 0.955 | 17.276 | 50 | 0.823 | True |
| TXN | 55/65/79 | 0.928 | 21.859 | 48 | 0.838 | True |
| UBA52 | 55/65/79 | 0.872 | 29.292 | 48 | 0.775 | True |
| UBB | 55/65/79 | 0.916 | 23.593 | 48 | 0.800 | True |
| UBE2A | 55/79 | 0.930 | 21.527 | 85 | 0.770 | True |
| UBE2I | 55/79 | 0.914 | 23.961 | 85 | 0.753 | True |
| UBE2N | 55/79 | 0.873 | 29.203 | 85 | 0.687 | True |
| UBE2S | 55/65/79 | 0.865 | 30.132 | 48 | 0.752 | True |
| UBQLN1 | 55/65/79 | 0.547 | 56.820 | 48 | 0.604 | False |
| VCP | 55/79 | 0.827 | 34.256 | 85 | 0.387 | True |
| WDTC1 | 55/65/79 | 0.466 | 62.243 | 48 | 0.609 | False |
| XPO1 | 55/79 | 0.835 | 33.420 | 85 | 0.691 | True |
| YWHAB | 55/65/79 | 0.650 | 49.472 | 48 | 0.640 | True |
| YWHAQ | 55/65/79 | 0.724 | 43.610 | 48 | 0.722 | True |

## 5. Primary readout (Dox_1000ng/mL, 3days; lines 55, 65, 79; batch LIB052667)

| quantity | value |
|---|---|
| candidate transgenes | 89 |
| eligible cells / qualifying cells (r >= 0.3) | 232 / 230 |
| qualifying transgenes n_q | 89 |
| non-finite c_corr | 0 |
| M = median c_corr | 0.4120 |
| 95% bootstrap CI (10,000, over transgenes) | 0.3535 to 0.4601 |
| CI contains M | True |
| arccos(M), degrees | 65.7 |
| specific / testable | 22 / 89 |
| chance baseline (5%) | 4.45 |
| binomial p (one-sided, approximate) | 3.64e-10 |
| reading | DONOR-SPECIFIC |

Reliability floor, per line (c_half = cosine between single-replicate halves; floor angle = arccos c_half). Built from transgene replicate splits, not from control-versus-control displacements:

| age | cells | median_c_half | q10_c_half | median_r | qualifying | median_floor_angle |
|---|---|---|---|---|---|---|
| 55 | 91 | 0.923 | 0.814 | 0.960 | 90 | 22.643 |
| 65 | 53 | 0.926 | 0.902 | 0.962 | 53 | 22.156 |
| 79 | 88 | 0.958 | 0.873 | 0.979 | 87 | 16.608 |

Per-transgene table (`results/donor/primary_transgenes.csv`; pairs in `primary_pairs.csv`, cells in `primary_cells.csv`). c_corr is the mean over the transgene's qualifying line pairs; s is the mean uncorrected cosine used for specificity; p95_comp the 95th percentile of its comparison set of size n_comp.

| transgene | lines | n_qual_pairs | c_obs_mean | c_corr | angle_corr | s | n_comp | p95_comp | specific |
|---|---|---|---|---|---|---|---|---|---|
| ACTG1 | 55/79 | 1 | 0.452 | 0.483 | 61.149 | 0.452 | 83 | 0.567 | False |
| ACTL6A | 55/65/79 | 3 | 0.032 | 0.035 | 87.986 | 0.032 | 46 | 0.364 | False |
| AKT1 | 55/79 | 1 | 0.677 | 0.709 | 44.843 | 0.677 | 83 | 0.454 | True |
| ATG5 | 55/65/79 | 3 | 0.255 | 0.263 | 74.762 | 0.255 | 46 | 0.476 | False |
| CCNA2 | 55/65/79 | 3 | 0.327 | 0.335 | 70.409 | 0.327 | 46 | 0.473 | False |
| CCNB1 | 55/65/79 | 3 | 0.020 | 0.019 | 88.930 | 0.020 | 46 | 0.299 | False |
| CDC20 | 55/65/79 | 3 | 0.415 | 0.435 | 64.183 | 0.415 | 46 | 0.494 | False |
| CDC34 | 55/79 | 1 | 0.353 | 0.362 | 68.786 | 0.353 | 83 | 0.486 | False |
| CDK1 | 55/65/79 | 3 | 0.038 | 0.042 | 87.607 | 0.038 | 46 | 0.322 | False |
| CDK2 | 55/79 | 1 | 0.390 | 0.399 | 66.511 | 0.390 | 83 | 0.579 | False |
| CDK4 | 55/79 | 1 | 0.620 | 0.639 | 50.313 | 0.620 | 83 | 0.585 | True |
| CDKN2A | 55/65/79 | 3 | 0.509 | 0.550 | 56.659 | 0.509 | 46 | 0.510 | False |
| CHEK1 | 55/79 | 1 | 0.362 | 0.468 | 62.069 | 0.362 | 83 | 0.408 | False |
| CHEK2 | 55/65/79 | 3 | 0.194 | 0.211 | 77.834 | 0.194 | 46 | 0.385 | False |
| CTNNB1 | 55/79 | 1 | 0.724 | 0.745 | 41.817 | 0.724 | 83 | 0.582 | True |
| CXCR2 | 55/65/79 | 3 | 0.407 | 0.452 | 63.129 | 0.407 | 46 | 0.539 | False |
| CYCS | 55/79 | 1 | 0.770 | 0.788 | 37.996 | 0.770 | 83 | 0.625 | True |
| DDX3X | 55/79 | 1 | 0.330 | 0.332 | 70.590 | 0.330 | 83 | 0.637 | False |
| EEF2 | 55/65/79 | 3 | 0.389 | 0.400 | 66.404 | 0.389 | 46 | 0.521 | False |
| ESR1 | 55/65/79 | 3 | 0.469 | 0.482 | 61.194 | 0.469 | 46 | 0.497 | False |
| FAU | 55/79 | 1 | 0.586 | 0.594 | 53.582 | 0.586 | 83 | 0.609 | False |
| GART | 55/79 | 1 | 0.419 | 0.429 | 64.604 | 0.419 | 83 | 0.676 | False |
| GSK3B | 55/65/79 | 3 | 0.436 | 0.444 | 63.609 | 0.436 | 46 | 0.526 | False |
| HDAC1 | 55/79 | 1 | 0.579 | 0.590 | 53.877 | 0.579 | 83 | 0.587 | False |
| HDAC6 | 55/79 | 1 | 0.229 | 0.234 | 76.497 | 0.229 | 83 | 0.449 | False |
| HES1 | 55/65/79 | 3 | 0.442 | 0.460 | 62.608 | 0.442 | 46 | 0.390 | True |
| HNRNPA2B1 | 55/79 | 1 | 0.333 | 0.353 | 69.299 | 0.333 | 83 | 0.396 | False |
| HSP90AB1 | 55/79 | 1 | 0.213 | 0.218 | 77.409 | 0.213 | 83 | 0.442 | False |
| HSP90B1 | 55/79 | 1 | 0.075 | 0.079 | 85.443 | 0.075 | 83 | 0.418 | False |
| HSPA5 | 55/79 | 1 | 0.603 | 0.614 | 52.132 | 0.603 | 83 | 0.536 | True |
| HSPA8 | 55/79 | 1 | 0.243 | 0.248 | 75.665 | 0.243 | 83 | 0.482 | False |
| HSPD1 | 55/79 | 1 | -0.140 | -0.144 | 98.259 | -0.140 | 83 | 0.335 | False |
| ISG15 | 55/65/79 | 3 | -0.072 | -0.131 | 97.510 | -0.072 | 46 | 0.286 | False |
| KAT7 | 55/65/79 | 3 | 0.320 | 0.335 | 70.452 | 0.320 | 46 | 0.407 | False |
| KDM6A | 55/79 | 1 | 0.467 | 0.477 | 61.483 | 0.467 | 83 | 0.415 | True |
| LRRK2 | 65/79 | 1 | 0.070 | 0.072 | 85.884 | 0.070 | 48 | 0.276 | False |
| MAPK14 | 55/79 | 1 | 0.183 | 0.186 | 79.253 | 0.183 | 83 | 0.451 | False |
| MAPK7 | 55/65/79 | 3 | 0.176 | 0.184 | 79.383 | 0.176 | 46 | 0.392 | False |
| MCM7 | 55/79 | 1 | 0.398 | 0.412 | 65.668 | 0.398 | 83 | 0.452 | False |
| MDM2 | 55/65/79 | 3 | 0.319 | 0.330 | 70.711 | 0.319 | 46 | 0.435 | False |
| NFKB1 | 55/65/79 | 3 | -0.060 | -0.061 | 93.476 | -0.060 | 46 | 0.317 | False |
| NOTCH1 | 55/79 | 1 | 0.811 | 0.818 | 35.145 | 0.811 | 83 | 0.546 | True |
| NPM1 | 55/79 | 1 | 0.413 | 0.425 | 64.855 | 0.413 | 83 | 0.440 | False |
| OGT | 55/65/79 | 3 | 0.180 | 0.186 | 79.286 | 0.180 | 46 | 0.391 | False |
| OSKM | 55/65 | 1 | 0.745 | 0.760 | 40.503 | 0.745 | 49 | 0.480 | True |
| P4HB | 55/65/79 | 3 | 0.289 | 0.307 | 72.140 | 0.289 | 46 | 0.389 | False |
| PARK7 | 55/65/79 | 3 | 0.414 | 0.427 | 64.725 | 0.414 | 46 | 0.495 | False |
| PIK3C2A | 55/65/79 | 3 | 0.627 | 0.653 | 49.202 | 0.627 | 46 | 0.569 | True |
| PIK3CA | 55/65/79 | 3 | 0.223 | 0.251 | 75.459 | 0.223 | 46 | 0.384 | False |
| PKM | 55/65/79 | 3 | 0.497 | 0.516 | 58.939 | 0.497 | 46 | 0.546 | False |
| PLK1 | 55/79 | 1 | 0.482 | 0.498 | 60.120 | 0.482 | 83 | 0.451 | True |
| PPIA | 55/65/79 | 3 | 0.591 | 0.610 | 52.438 | 0.591 | 46 | 0.585 | True |
| PPP1CA | 55/65/79 | 3 | 0.163 | 0.173 | 80.067 | 0.163 | 46 | 0.377 | False |
| PPP1CB | 55/65/79 | 3 | 0.505 | 0.641 | 50.145 | 0.505 | 46 | 0.536 | False |
| PPP1CC | 55/65/79 | 3 | 0.329 | 0.340 | 70.096 | 0.329 | 46 | 0.479 | False |
| PPP2CA | 55/79 | 1 | 0.521 | 0.531 | 57.934 | 0.521 | 83 | 0.585 | False |
| PPP3CA | 55/65/79 | 3 | 0.444 | 0.469 | 62.013 | 0.444 | 46 | 0.498 | False |
| PPP4C | 55/65/79 | 3 | 0.362 | 0.372 | 68.140 | 0.362 | 46 | 0.451 | False |
| PSMC6 | 55/65/79 | 3 | 0.419 | 0.438 | 64.025 | 0.419 | 46 | 0.505 | False |
| PSMD4 | 55/79 | 1 | -0.190 | -0.202 | 101.628 | -0.190 | 83 | 0.433 | False |
| PTEN | 55/65/79 | 3 | 0.376 | 0.385 | 67.327 | 0.376 | 46 | 0.481 | False |
| RAD23B | 55/65/79 | 3 | 0.376 | 0.391 | 66.978 | 0.376 | 46 | 0.460 | False |
| RAD51 | 55/79 | 1 | 0.689 | 0.700 | 45.564 | 0.689 | 83 | 0.448 | True |
| RNASEL | 55/79 | 1 | 0.357 | 0.376 | 67.920 | 0.357 | 83 | 0.339 | True |
| RPLP0 | 55/65/79 | 3 | 0.649 | 0.688 | 46.520 | 0.649 | 46 | 0.636 | True |
| RPS27A | 55/65/79 | 3 | 0.441 | 0.462 | 62.476 | 0.441 | 46 | 0.499 | False |
| RRM1 | 55/79 | 1 | 0.667 | 0.688 | 46.527 | 0.667 | 83 | 0.588 | True |
| RUVBL2 | 55/79 | 1 | 0.767 | 0.840 | 32.842 | 0.767 | 83 | 0.652 | True |
| SIRT1 | 55/65/79 | 3 | 0.560 | 0.579 | 54.593 | 0.560 | 46 | 0.533 | True |
| SLC2A10 | 55/65/79 | 3 | 0.485 | 0.529 | 58.045 | 0.485 | 46 | 0.520 | False |
| SLC2A13 | 55/65/79 | 3 | 0.333 | 0.344 | 69.852 | 0.333 | 46 | 0.464 | False |
| SOD1 | 55/65/79 | 3 | 0.592 | 0.616 | 51.996 | 0.592 | 46 | 0.575 | True |
| SOD2 | 55/65/79 | 3 | 0.566 | 0.587 | 54.062 | 0.566 | 46 | 0.565 | True |
| SRSF1 | 55/65/79 | 3 | 0.771 | 0.785 | 38.235 | 0.771 | 46 | 0.193 | True |
| SUMO1 | 55/79 | 1 | 0.068 | 0.071 | 85.947 | 0.068 | 83 | 0.399 | False |
| TP53 (DN) | 55/65 | 1 | 0.781 | 0.800 | 36.914 | 0.781 | 49 | 0.758 | True |
| TXN | 55/65/79 | 3 | 0.448 | 0.466 | 62.234 | 0.448 | 46 | 0.515 | False |
| UBA52 | 55/65/79 | 3 | 0.042 | 0.047 | 87.324 | 0.042 | 46 | 0.350 | False |
| UBB | 55/65/79 | 3 | 0.249 | 0.259 | 75.004 | 0.249 | 46 | 0.397 | False |
| UBE2A | 55/79 | 1 | 0.558 | 0.579 | 54.625 | 0.558 | 83 | 0.653 | False |
| UBE2I | 55/79 | 1 | 0.156 | 0.163 | 80.596 | 0.156 | 83 | 0.450 | False |
| UBE2N | 55/79 | 1 | 0.067 | 0.072 | 85.868 | 0.067 | 83 | 0.615 | False |
| UBE2S | 55/65/79 | 3 | 0.357 | 0.383 | 67.507 | 0.357 | 46 | 0.543 | False |
| UBQLN1 | 55/65/79 | 1 | -0.033 | -0.034 | 91.946 | -0.033 | 48 | 0.270 | False |
| VCP | 55/79 | 1 | -0.173 | -0.191 | 101.006 | -0.173 | 83 | 0.175 | False |
| WDTC1 | 55/65/79 | 1 | 0.229 | 0.260 | 74.953 | 0.229 | 49 | 0.532 | False |
| XPO1 | 55/79 | 1 | 0.387 | 0.427 | 64.721 | 0.387 | 83 | 0.514 | False |
| YWHAB | 55/65/79 | 3 | 0.549 | 0.768 | 39.855 | 0.549 | 46 | 0.577 | False |
| YWHAQ | 55/65/79 | 3 | -0.096 | -0.091 | 95.226 | -0.096 | 46 | 0.324 | False |

## 6. The reading that fired, verbatim

> DONOR-SPECIFIC. A transgene's effect does not carry between aged donors even at bulk depth: the median reliability-corrected cosine between the same transgene's displacements in different donor lines is 0.412 (95% CI 0.353 to 0.460, 65.7 degrees) over 89 qualifying transgenes, and 22 of 89 transgenes were specific against a chance baseline of 4.45 (binomial p = 3.64e-10, approximate). A forward model cannot borrow a factor's direction across donors. The combination screen needs effects measured per donor and more donors than planned. The same-context control recovered M_same = 0.934 at single-replicate depth. Measured in three lines aged 55, 65 and 79, under one dox condition, in one library batch. This compares aged donors with each other. No reading licenses a claim about transfer between an aged and a young context, nor about the Southard screen.

What the numbers beside it say (descriptive, not gated). The reading fired on its CI branch. The measured agreement is not zero: the median corrected cosine is 0.412 (65.7 degrees), and 22 of 89 transgenes are specific against 4.45 expected by chance. So the directions are partly shared between lines, and transgene-specifically so. What is ruled out, at the pre-registered bar, is agreement close enough to borrow a direction across lines. Within one line the same transgene's replicate halves agree at 0.934, so the shortfall is not measurement noise at this depth. The wording's "does not carry" is stronger than these numbers. It is reproduced because it is the pre-registered text of the branch that fired.

Four things that limit how far the reading reaches:

1. **Line pairs differ a lot.** Median corrected cosine by line pair: 55-65: 0.721 over 50 pairs; 55-79: 0.388 over 84 pairs; 65-79: 0.180 over 49 pairs. Line sets per transgene: 55/65/79 49, 55/79 37, 55/65 2, 65/79 1. So the median over transgenes partly reflects which line pairs each transgene happens to have, and that is an accident of the design.

| pair | n_pairs | median_c_obs | median_c_corr | angle_median_c_corr |
|---|---|---|---|---|
| 55-65 | 50 | 0.685 | 0.721 | 43.852 |
| 55-79 | 84 | 0.357 | 0.388 | 67.197 |
| 65-79 | 49 | 0.170 | 0.180 | 79.620 |

2. **Re-measuring the same line moves the direction too** (D2, diagnostic). For the 9 transgenes measured under both doses, the same line re-measured in another batch at another dose agrees at median corrected 0.257, below two different lines in the same batch (0.504). Batch and dose are confounded there, so this does not separate them. But re-measuring one line in another batch at another dose moved the direction at least as much as changing the line did. "Donor-specific" here cannot be told apart from "specific to the experiment".

3. **What Rep1 and Rep2 are is not stated.** The floor is the agreement between the record's Rep1 and Rep2 GSMs in one line and batch. If those are technical replicates or wells of one culture, the floor is set by technical noise only, and culture-to-culture variation within a donor ends up in the cross-line drift. Item 2 is consistent with that. The gate's 0.934 is therefore a recovery check at replicate depth, not evidence that a second culture of the same donor would agree at 0.934.

4. **The 20-component space holds little of each displacement** (D6): median share of squared norm retained 0.131 (5-95%: 0.039 to 0.274). The angles are measured in a GTEx-derived subspace that holds about an eighth of each displacement. The basis was chosen for the Southard screen, not for this data.

## 7. Diagnostics (each one diagnostic, not a gated result)

**D1 — line mean removed** (diagnostic, not a gated result). n_q 89; M 0.2195 (95% CI 0.1304 to 0.2907, 77.3 degrees); specificity 28 of 89 against 4.45, p = 1.95e-15. Removing each line's mean displacement lowers agreement and raises specificity. Part of the primary agreement is a component that all transgenes in a line share against mTagBFP2. Disclosed: the line mean includes the transgene itself, so this is not independent of the vectors it corrects. Tables `d1_line_mean_removed_*.csv`.

**D2 — batch versus donor** (diagnostic, not a gated result). In the primary all three lines share batch LIB052667, so primary donor differences cannot track batch. Across conditions, for the transgenes present under both doses (HDAC1, OGT, OSKM, PPP1CB, RNASEL, SIRT1, SLC2A13, SRSF1, WDTC1):

| category | n_pairs | median_c_obs | n_corr | median_c_corr | median_angle |
|---|---|---|---|---|---|
| different line / different batch | 188 | 0.162 | 177 | 0.194 | 80.702 |
| different line / same batch | 72 | 0.391 | 68 | 0.504 | 66.975 |
| same line / different batch | 22 | 0.156 | 21 | 0.257 | 81.039 |

**D3 — agreement and replicate count** (diagnostic, not a gated result). Every eligible primary cell has [2] transgene and [2] control GSMs, so replicate count does not vary and its relation to agreement is undefined. Spearman of c_corr(t) with the mean library size of t's GSMs: rho = -0.065, p = 0.5436; with t's number of qualifying lines: rho = -0.085, p = 0.4285 (n = 89). No transgene shortlist is reported, so there is no shortlist to check against cell count.

**D4 — full gene space** (diagnostic, not a gated result). Over 183 qualifying pairs, median cosine in the 23,485-column z space 0.090 (84.8 degrees) against 0.413 (65.6 degrees) in the 20-component space. That is close to orthogonal, as TRANSFER expected: shared signal is diluted by noise spread over 23,485 dimensions. The flag said the expectation was weaker at bulk depth; it held anyway. Table `d4_full_gene_space_pairs.csv`.

**D5 — Dox 2000 six-line set** (diagnostic, not a gated result). 9 qualifying transgenes (fewer than 20, so it could not be a reading); median c_corr 0.385 (67.3 degrees); no specificity test (comparison sets below 19). Two batches: 55/65/79 in LIB053907, 67/68/69 in LIB054839, so most of its pairs cross a batch (see D2).

| transgene | lines | n_qual_pairs | c_obs_mean | c_corr | angle_corr |
|---|---|---|---|---|---|
| HDAC1 | 55/65/67/68/69/79 | 15 | 0.465 | 0.619 | 51.774 |
| OGT | 55/65/67/68/69/79 | 15 | -0.095 | -0.119 | 96.844 |
| OSKM | 67/68/69 | 3 | 0.897 | 1.025 | n/a |
| PPP1CB | 55/65/67/68/69/79 | 15 | 0.232 | 0.385 | 67.350 |
| RNASEL | 55/65/67/68/69/79 | 15 | 0.093 | 0.116 | 83.358 |
| SIRT1 | 55/65/67/68/69/79 | 10 | 0.295 | 0.399 | 66.500 |
| SLC2A13 | 55/65/67/68/69/79 | 15 | 0.270 | 0.360 | 68.901 |
| SRSF1 | 55/65/67/68/69/79 | 15 | 0.923 | 0.941 | 19.715 |
| WDTC1 | 55/65/67/68/69/79 | 15 | -0.032 | -0.096 | 95.502 |

**D6 — norm bookkeeping** (diagnostic, not a gated result). Share of each displacement's squared z-space norm retained by the projection: median 0.131, 5% 0.039, 95% 0.274. Per cell in `d6_norms.csv`.

**D7 — tables.** `cells.csv` (lines, GSMs, controls, passages, flags, platforms, library size), `primary_cells.csv` (c_half, floor angle, r), `primary_pairs.csv` (c_obs, angle, r_A, r_B, c_corr), `primary_transgenes.csv`, `gate_table.csv`, `panel_gsms.csv`.

## 8. What this does not cover (flag section 10, verbatim)

```
- Six donor lines, all aged 55-79, and the primary uses three of them (55,
  65, 79). No young donor and no cell type other than dermal fibroblast.
- Passage differs between lines (and between transgene and control GSMs
  within some cells).
- This says nothing about whether any transgene approaches young: no young
  target enters anywhere.
- The authors' own age predictor is not used anywhere.
- One dox condition (1000 ng/mL, 3 days) in one library batch. Agreement
  under other doses, schedules or batches is not measured by the primary.
- Two replicates per cell: the reliability correction rests on one
  single-replicate split per cell (two control pairings), so r is itself a
  noisy estimate and c_corr inherits that noise.
- Within a line every transgene shares the same 2 mTagBFP2 GSMs, so control
  noise is common to all displacements in that line. It cancels from no
  cross-line cosine but adds a shared component within a line (D1 removes the
  line mean).
- The 20-component basis was chosen in SOUTH7/SOUTH9 for its sensitivity on
  the Southard screen; it is applied here unchanged and is not chosen for
  this data.
- The specificity binomial is approximate (dependent comparisons).
- Direction only: nothing about effect magnitude or usefulness.
```

Added here. None of these lines is young or neonatal, and nothing here transfers to the Southard (Hs27) or Sengstack cells. Of the six lines, 67, 68 and 69 enter only D5.

## 9. Files

`src/donor_run.py`; `results/donor/` (PREREG.flag, run.log, stage0_reproduction.csv, stage0.json, build.json, excluded_gsms.csv, cells.csv, panel_gsms.csv, gate.json, gate_table.csv, primary.json, primary_*.csv, diagnostics.json, d1/d2/d4/d5/d6 CSVs, displacements.npz). `results/donor/displacements_full.npy` (26.6 MB, the full z-space displacements) is git-ignored and regenerated by `python src/donor_run.py build`. Raw inputs, git-ignored under `data/raw/srsfinv/`: the family SOFT (SRSFINV) and `GSE325735_human_raw_counts.txt.gz` (https://ftp.ncbi.nlm.nih.gov/geo/series/GSE325nnn/GSE325735/suppl/GSE325735_human_raw_counts.txt.gz).
