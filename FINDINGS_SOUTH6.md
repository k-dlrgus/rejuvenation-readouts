# FINDINGS_SOUTH6 — Southard ranking with the control set defined by the guide_target label

**Status:** Ran to completion. No STOP key fired.

**toward_young counts (not pooled):** S1 O1_Y1=0; S2_k20 O1_Y1=0; S2_k50 O1_Y1=0; S2_k100 O1_Y1=0; S1 O2_Y1=0; S2_k20 O2_Y1=15; S2_k50 O2_Y1=9; S2_k100 O2_Y1=10; S2_k20 O1_Y2=0; S2_k50 O1_Y2=0; S2_k20 O2_Y2=17; S2_k50 O2_Y2=10.

**Control set:** 78 label-defined non-targeting guides, 14013 cells; 2179 cells in 17 guides moved out of the targeting pool relative to SOUTH5.

**Stage 2B robustness:** the 78-guide and 61-guide control floors differ materially in 12 of 12 settings.

The fenced pre-registration was written to this file and to `results/south6/PREREG.flag` before any statistic. The fence is the flag file, not a paraphrase.

SOUTH6 supersedes SOUTH5. `results/south5/` is kept.

## Pre-registration (verbatim, written before any statistic)

```
TASK: SOUTH6 — SOUTH5 with the control set defined by label, not by flag.

Nature of this task: computational reanalysis of a public human cell dataset.
No lab work, no protocol, no new data. Output is a ranked table.

SOUTH5 stopped with `control_ambiguous`: the file holds 78 non-targeting guides
but `control=True` flags only 61 of them, leaving 17 guides (2,179 cells) in the
targeting pool. SOUTH6 supersedes SOUTH5 and resolves this by defining the
control set from `guide_target`, not from `control`. Reason, stated now before
any statistic: `guide_target` is the per-guide annotation of what each guide is,
`control` is a derived summary that omits 17 of them, and SURVEY_TF.md recorded
78 non-targeting guides before any SOUTH task ran. The two guide sets are
disjoint at guide level, as SOUTH5 verified, so this is a clean relabelling and
not a repair of individual rows.
Add one line at the top of FINDINGS_SOUTH5.md marking it superseded; change
nothing else in it. Keep results/south5/.

EXECUTION RULES, binding (unchanged from SOUTH5):
- On a STOP, write the findings and exit. Never compute after a STOP, never
  resume a stopped run.
- No checkpoint/resume logic unless resumed random draws are verified
  bit-identical to an uninterrupted stream, and that check is reported.
- FINDINGS and PROGRESS must report the same stop status; every counted failure
  printed verbatim.
- src/south6_run.py runs standalone (`python src/south6_run.py`) and refuses to
  compute if PREREG.flag is missing. Write it, smoke-test Stage 0 only, stop.
  Do not run the full job yourself.

Write FINDINGS_SOUTH6.md, PROGRESS_SOUTH6.md, results/south6/*, src/south6_run.py.
Modify no other FINDINGS_*, PROGRESS_*, or src/ file. Import from
src/toward_run.py, src/same_run.py, src/seng_run.py, src/lowdim_common.py
without changing them. Do not refit the frozen ruler.
Seeds: 20260914, 20260918. n_perm=5000, n_boot=200.

WRITE THE PRE-REGISTRATION FIRST: this file verbatim into FINDINGS_SOUTH6.md
and results/south6/PREREG.flag before any statistic.

--- DATA ---
The same file SOUTH5 verified: Southard Hs27 per-cell singlets,
md5 2c926c7cae2eaab1ed2519b06ca21db5, 8,048,055,098 bytes, in data/raw/southard/.
Confirm the md5 and do not re-download. Counts come from layers/counts, which
SOUTH5 confirmed holds raw integers; X is log1p-normalised and is not used.

--- CONTROL SET (the one change) ---
Non-targeting cells are all cells whose `guide_target` is the literal string
'non', across all 78 guides, regardless of `control`.
Report: the 78 guides, their cell counts, and how many are control=True (expect
61) versus control=False (expect 17). STOP `control_count_unexpected` if the
total is not 78 guides, or if any guide has cells under both control states.
Targeting cells are all cells whose `guide_target` is anything other than 'non'.
Report how many cells move from the targeting pool to the control pool relative
to SOUTH5's definition.
Cell filtering, guide and target columns: as SOUTH5 (`guide_identity`,
`guide_target`, `keep`, `single_cell`).

--- ROBUSTNESS (report-only, never a gate) ---
Compute the Stage 2B control floor twice: once on all 78 guides, once on the 61
control=True guides only. Report both and state whether they differ materially.
If they do, say so plainly in the findings; do not change the primary analysis,
which uses 78.

--- STAGE 0 through STAGE 5 ---
Identical to the SOUTH5 pre-registration, reproduced verbatim below this line in
FINDINGS_SOUTH6.md, with the control-set definition above replacing its
"Non-targeting cells are those flagged by `control`" sentence and its
`control_ambiguous` stop, and no other change.

--- STAGE 0: INVENTORY (descriptive only) ---
Report: matrix shape before and after filtering; whether X is raw integer
counts (STOP `no_raw_counts` if not and no raw layer exists); distinct target
genes; guides per target (min/median/max); non-targeting guides and cells;
cells per target.
Ruler coverage three ways (symbol, version-stripped Ensembl, best after
src/md3_idtype.py), each as gene count AND fraction of the ruler's total
absolute weight. Use the highest-weight mapping; name it. Record
`coverage_low` if that fraction is below 0.25. SURVEY_TF.md expects ~0.63;
report the actual number whatever it is.
STOP `too_few_factors` if fewer than 500 targets have 2+ guides.
STOP `no_controls` if fewer than 20 non-targeting guides.

--- STAGE 1 through STAGE 5 ---
Identical to the SOUTH4 pre-registration, which is reproduced verbatim below
this line in FINDINGS_SOUTH5.md, with two changes and no others:
  (a) n_perm = 5000 instead of 20000. Reason stated now, before any number:
      S1 here has ~7x the genes of SOUTH3, so each draw costs ~7x more.
      At 5000 the smallest p is 1/5001 and BH across 1,836 factors reaches
      q <= 0.05 once 8 factors tie at that minimum, against 2 at 20000.
      Report the attainable minimum p and how many factors tie at it.
  (b) The guide columns are those named above.
[Paste the SOUTH4 prereg text from "--- STAGE 1: VECTORS AND SPACES ---"
through "--- OUTPUT ---" here verbatim.]

--- OUTPUT ---
As SOUTH4, with results written to results/south5/ and the candidate table to
results/south5/top20.csv. PROGRESS_SOUTH5.md gives the stop status, the exact
terminal command to run the full job, and the expected runtime.
Record every failure verbatim. Do not substitute columns or repair rows.

================================================================================
SOUTH4 pre-registration, Stages 1–5, reproduced verbatim (governed by the two
changes stated above: n_perm = 5000, and the named guide columns):
================================================================================

--- STAGE 1: VECTORS AND SPACES ---
Pseudobulk per target gene (all its cells) and per non-targeting guide.
Process exactly as TOWARD/SAME: sum counts, TMM, log2 CPM prior.count=2
across this task's own panel, frozen GTEx mu/sd z-score, missing genes z=0.
Effect vector d_g = z(target g) - z(all non-targeting cells pooled).
One gene variant only: all matched genes. No significance masking; SOUTH3
showed masking makes guides orthogonal by construction.
Spaces: S1 gene space; S2 low-dim PCA fit ONLY on the GTEx fibroblast donor
matrix, k in (20, 50, 100), applied unchanged. Report every k.

--- ORIGINS AND TARGETS ---
PRIMARY, declared now and reported first: O1 = GTEx donors aged 60-79,
Y1 = GTEx donors aged 20-39, in S1 and each S2 k.
SECONDARY, reported but never primary: O2 = Lu aged donor GM00731 day-0
pseudobulk; Y2 = young donors as a region (mean, covariance, Mahalanobis,
radius at the 95th percentile of held-out young donors), S2 only, skipped
at any k where the covariance is singular.
Reason for the primary choice, stated now: SOUTH3 found ||O2-Y1|| is 7-15x
||O1-Y1||, so O2 makes any favourable projection lower delta, and the aged
GTEx origin already sits inside the Y2 region, so Y2 membership does not
discriminate at O1. This choice is made before any SOUTH4 number exists.
Report each origin-to-target distance before scoring anything.

--- STAGE 2: DETECTION LIMIT (before any ranking) ---
A. Primary floor: permute each factor's own vector across genes 20,000 times;
   delta = ||(origin + d_perm) - target|| - ||origin - target||; the floor is
   the 5th percentile. Report per-factor and pooled-at-median-norm floors.
B. Control floor: same delta using each non-targeting guide as a pseudo-effect,
   and using random splits of the non-targeting cells. Report both. With 78
   guides this is now a real floor; report it alongside A and say which is
   stricter.
C. Dose-response: add the true young-minus-old direction to the origin at
   f = 0, 0.05, 0.1, 0.25, 0.5, 1.0, scaled to the median factor norm.
   MDA = smallest f clearing the primary floor.
   STOP `no_detection_power` if no f <= 1.0 clears it anywhere.
Never subtract a floor. It is a bar to clear.

--- STAGE 3: RANKING ---
Per factor, per setting: cos, frac, delta, and for Y2 the Mahalanobis change
and membership. Guide agreement = median pairwise cosine between that factor's
own guide-level vectors, computed on unmasked vectors; flag `guides_disagree`
if <= 0. Identity check on COL1A1, COL1A2, FN1, LUM, PDGFRA, PDGFRB, POSTN,
PRRX1, SERPINH1, VIM (report which are present); flag `identity_loss` if the
drop exceeds 0.5. Proliferation-residualized cos and delta as in seng_run.py.
Frozen ruler score change as context only.
Nulls: Stage 2A permutations, one-sided p, BH across factors within setting.
Bootstrap over cells; flag any interval not containing its point estimate as
INVALID and never gate on it.
Label `toward_young` only if delta < 0 AND below the factor's own Stage 2A
floor AND q <= 0.05 AND not identity_loss AND not guides_disagree.
Report counts per setting; never pool. Zero is a valid answer.

--- STAGE 4: DOES DELTA ADD ANYTHING? ---
Repeat SOUTH3's check at this coverage. Report ‖d‖/‖target-origin‖ per factor;
Spearman between delta and cosine rankings overall and by quartile of that
ratio; the residual of delta ≈ -‖d‖cos(theta); and what fraction of survivors
a plain top-n-by-cosine recovers. State plainly whether delta separates from
cosine at this coverage, and in which ratio quartile if anywhere.

--- STAGE 5: COMBINATIONS ---
State at the top of the section that this assumes effects add linearly and
that no public fibroblast screen can test that assumption.
Greedy forward search on summed vectors: add the factor that most reduces
distance; stop when nothing reduces it, or on entering the Y2 region, or at
size 5. Repeat from 50 random starts (seed 20260914); report how often the
same set recurs. Per set: delta, margin over floor, ‖d‖/‖gap‖ ratio, identity
flag, and delta with each member removed.
Report no set unless its delta clears the floor by at least the Stage 2B MDA
margin. Report whether combinations reach a ratio where Stage 4 showed delta
and cosine separating.

--- WHAT DOES NOT COUNT ---
- cosine or ruler-score drop alone as evidence of approach.
- any delta not clearing the measured primary floor.
- subtracting a floor instead of clearing it.
- choosing k, mapping, origin or target after seeing results.
- substituting a smaller file for the pinned one.
- dropping or relabelling factors because they look odd; report them.
- pooling settings into one number.
- claims about the source paper's own conclusions.

--- OUTPUT ---
FINDINGS_SOUTH4.md: prereg verbatim; Stage 0 inventory with all three coverage
numbers; Stage 2 floors and MDA; Stage 3 counts plus a top-20 candidate table
(factor, delta, margin over floor, q, guide agreement, identity flag,
proliferation flag) written to results/south4/top20.csv; Stage 4 verdict;
Stage 5 combinations; fired keys; limitations, which must include: Hs27 is one
neonatal line with no aged cells, so every delta is a counterfactual assuming
the effect transfers; combinations assume additivity, untestable here; unmeasured
genes assumed unmoved.
PROGRESS_SOUTH4.md: stop status, next action.
Record every failure verbatim. Do not substitute columns or repair rows.

--- OUTPUT ---
As SOUTH5, writing to results/south6/ and results/south6/top20.csv.
Limitations must additionally include: the control set was defined by the
`guide_target` label rather than the file's own `control` column, which omits
17 non-targeting guides; both floors are reported.
PROGRESS_SOUTH6.md: stop status, the terminal command, expected runtime.
Record every failure verbatim. Do not substitute columns or repair rows.
```

## Stage 0 — inventory

File `fibroblast_CRISPRa_final_pop_singlets_normalized_log1p.h5ad`, 8048055098 bytes, md5 `2c926c7cae2eaab1ed2519b06ca21db5` (expected `2c926c7cae2eaab1ed2519b06ca21db5`, match). Not re-downloaded.

Matrix source `layers/counts (csr)`. X was audited and is not raw counts; the counts layer is. Shape before filtering [220403, 30395]; after keeping singlets that pass `keep` [220403, 30395].

Distinct target genes (excluding `non`): 1836. Guides per target: min 2, median 6.0, max 6. Targets with >=2 guides: 1836.

Non-targeting guides: 78. Non-targeting cells: 14013. Cells per target gene: min 13, median 106.0, max 341.

Ruler coverage (gene count and fraction of total absolute weight):

| mapping | genes matched | fraction of genes | fraction of \|w\| |
| --- | --- | --- | --- |
| symbol | 16631 / 23485 | 0.7082 | 0.6204 |
| ensembl_version_stripped | 20257 / 23485 | 0.8626 | 0.8177 |
| best_after_md3_idtype | 20265 / 23485 | 0.8629 | 0.8182 |

Mapping used, highest weight fraction: `best_after_md3_idtype` (0.8182). SURVEY_TF.md expected ~0.63. `coverage_low` did not fire.

Identity genes present: COL1A1, COL1A2, FN1, LUM, PDGFRA, PDGFRB, POSTN, PRRX1, SERPINH1, VIM. Absent: none.

## Control set — defined by the `guide_target` label, not by `control`

Definition used: guide_target == 'non' (exact string equality). This is the one change from SOUTH5. `control` is reported below but never partitions the data, and no row was edited.

Matrix shape before filtering: [220403, 30395]. After keeping singlets that pass `keep`: [220403, 30395].

Cells before filtering: 220403. Dropped not-singlet: 0. Dropped singlet-fails-keep: 0. Cells after filtering: 220403.

Distinct values before filtering:

- `control`: {False: 208569, True: 11834}
- `keep`: {True: 220403}
- `single_cell`: {True: 220403}

Non-targeting guides found: **78** (pre-registration expected 78), carrying 14013 kept cells. Of those guides, 61 are `control=True` (expected 61) with 11834 cells, and 17 are `control=False` (expected 17) with 2179 cells.

Guides with cells under both control states: none. Guides carrying both the `non` label and another `guide_target` label: none. Kept cells with `control=True` but `guide_target` != `non`: 0.

### Movement relative to SOUTH5's control-flag definition

SOUTH5 defined the control set as `control=True`: 61 guides, 11834 cells. SOUTH6 defines it by label: 78 guides, 14013 cells. **2179 cells in 17 guides move from the targeting pool to the control pool.** The targeting pool is 206390 cells here against 208569 under SOUTH5's definition.

### The 78 non-targeting guides

| guide | cells | control=True cells | control=False cells |
| --- | --- | --- | --- |
| non_targeting_GCGCCGGAGCACTCGCGTAT | 439 | 439 | 0 |
| non_targeting_GTGCAGCTGGTCGCGCAAAC | 411 | 411 | 0 |
| non_targeting_GAACGGGCCGTGATCGGACC | 400 | 400 | 0 |
| non_targeting_GGCCCTCCCCACCGGCACGA | 325 | 325 | 0 |
| non_targeting_GAAATGGACCGCGCTTACGC | 323 | 323 | 0 |
| non_targeting_GCTCTACACTGGATGATCGT | 278 | 278 | 0 |
| non_targeting_GCCCGTGATAATCGATACGA | 250 | 250 | 0 |
| non_targeting_GCTCTTGGTACGTATTCGAA | 250 | 250 | 0 |
| non_targeting_GGCCGGAGCCGTTCGCTAGC | 238 | 238 | 0 |
| non_targeting_GCTTCGGGGGAATCACCGGT | 236 | 236 | 0 |
| non_targeting_GCCAGACGCGCCCGTAACGG | 235 | 235 | 0 |
| non_targeting_GGCTGTACGCTCGCGTAAAA | 235 | 235 | 0 |
| non_targeting_GCTAGGACCGCACGTCGCGT | 233 | 233 | 0 |
| non_targeting_GAGGAGCGTCATATCACGCG | 231 | 231 | 0 |
| non_targeting_GTGATTAGAGTCCGGTAACG | 226 | 226 | 0 |
| non_targeting_GGGTAGGCCTCGCCAGTCGT | 224 | 224 | 0 |
| non_targeting_GGGACTCCCTCGGCGGTTAT | 219 | 219 | 0 |
| non_targeting_GCAGCCTAAAGCGTACGACG | 216 | 216 | 0 |
| non_targeting_GTTTGATTTAAGCGTAATAC | 210 | 210 | 0 |
| non_targeting_GTTAAGGCACTCGTATGCGA | 206 | 206 | 0 |
| non_targeting_GATAAGTCGTCGTCCAATCG | 205 | 205 | 0 |
| non_targeting_GTTGTGATGTATCACGACTA | 198 | 198 | 0 |
| non_targeting_GAGCTGTCATTGTCGTACGG | 194 | 194 | 0 |
| non_targeting_GAAGATAAAGCGCGTATTGC | 193 | 193 | 0 |
| non_targeting_GACGCGACACGAACTAGCGT | 193 | 193 | 0 |
| non_targeting_GAGCCGATTGACGCATATAG | 192 | 192 | 0 |
| non_targeting_GACCACCGACGCGTTACGCA | 190 | 190 | 0 |
| non_targeting_GCAAAGACCACGTCATCGGT | 190 | 190 | 0 |
| non_targeting_GTGATCTTTGCTCGCATAGA | 188 | 188 | 0 |
| non_targeting_GGTCGTTTACATACTCGCAC | 186 | 186 | 0 |
| non_targeting_GCCGTTCGATACGCGATCTA | 184 | 184 | 0 |
| non_targeting_GCGTTGCAACCGATCGTAAG | 184 | 184 | 0 |
| non_targeting_GTAATGGACGCGTACCGTTC | 180 | 180 | 0 |
| non_targeting_GCTAACCAGACGCCGTCAAT | 179 | 179 | 0 |
| non_targeting_GGACTATGGGCGCCGATCAG | 177 | 177 | 0 |
| non_targeting_GTAGGGAGTATACCTCGCGA | 175 | 175 | 0 |
| non_targeting_GGGTCGTAGAACTAGCGCAA | 173 | 173 | 0 |
| non_targeting_GGCAGGGGAGTCGACTATTC | 168 | 168 | 0 |
| non_targeting_GTCGTGGAGAGCGCTAGCGT | 166 | 166 | 0 |
| non_targeting_GTGCATCAAGTCGACCGGCG | 164 | 164 | 0 |
| non_targeting_GCCCTGTACCACGTACCGTC | 161 | 161 | 0 |
| non_targeting_GGCGGGGAACAATAGGCGTA | 160 | 160 | 0 |
| non_targeting_GACCAGCATTCATACGCCGG | 153 | 153 | 0 |
| non_targeting_GGCAAGGTGACGCCTAGTAC | 151 | 151 | 0 |
| non_targeting_GCATACCAGTCGAACCGGAC | 142 | 142 | 0 |
| non_targeting_GCTCTGTTCGTCGTTCGCTG | 141 | 141 | 0 |
| non_targeting_GGATGTCGCCACGTACGGAT | 141 | 141 | 0 |
| non_targeting_GGGGTGTAGTCGCGTGTAGT | 141 | 141 | 0 |
| non_targeting_GCGGACAAAGCGATATGTTA | 140 | 140 | 0 |
| non_targeting_GTCCACAAGACGTGCTCGCA | 138 | 138 | 0 |
| non_targeting_GTCACCCAGCCGATATATAA | 134 | 134 | 0 |
| non_targeting_GAAACAAACAGACGGACCTA | 130 | 130 | 0 |
| non_targeting_GTGGGTGAATGCGTACGTCG | 125 | 125 | 0 |
| non_targeting_GTTGATCGACCGATTGCGAT | 121 | 121 | 0 |
| non_targeting_GGATGAAAACCCGAACCGAT | 120 | 120 | 0 |
| non_targeting_GCGGGTAGCTTATATAGCGC | 116 | 116 | 0 |
| non_targeting_GTCAATACGACGAGCCGACA | 114 | 114 | 0 |
| non_targeting_GACACTGGACGCGACTTATA | 113 | 113 | 0 |
| non_targeting_GATGCTCGCAGTTCGACCGC | 112 | 112 | 0 |
| non_targeting_GTGTTTTTCGGTCGGCCGAT | 110 | 110 | 0 |
| non_targeting_GTGCTCGCGAATCGTATCGG | 107 | 107 | 0 |
| non_targeting_GAAGCAGCATGAATACGCCG | 376 | 0 | 376 |
| non_targeting_GGGAACAGGGGCGGTCCGTA | 199 | 0 | 199 |
| non_targeting_GTCCATAGGGTCTAGCGCCG | 189 | 0 | 189 |
| non_targeting_GGAATCCAGCTCGACGACCA | 156 | 0 | 156 |
| non_targeting_GACGAATGAAGCGTCGATAA | 153 | 0 | 153 |
| non_targeting_GGCTTCTACACCGCGATGAC | 137 | 0 | 137 |
| non_targeting_GAGATATGAGGCGACGATAT | 127 | 0 | 127 |
| non_targeting_GCCTAAATACTATTCGCGGA | 100 | 0 | 100 |
| non_targeting_GGATTCTGGAAAGTCGGTCG | 98 | 0 | 98 |
| non_targeting_GCTGTTTCGACCCGTCGAAT | 96 | 0 | 96 |
| non_targeting_GTCGTTATCTCGCTATTTCG | 92 | 0 | 92 |
| non_targeting_GGAATCTTCGCGTAACGAGC | 87 | 0 | 87 |
| non_targeting_GTGGATTATCTGACGCGAAT | 82 | 0 | 82 |
| non_targeting_GACTGAAAGCCGATATCGGG | 81 | 0 | 81 |
| non_targeting_GTGGTTATACCCGACTAGAC | 77 | 0 | 77 |
| non_targeting_GTTCCTGTTGGGTCGCGAAT | 77 | 0 | 77 |
| non_targeting_GTAAATTCTCGCGTAACGTT | 52 | 0 | 52 |

The same table is `results/south6/control_set_guides.csv`.

## Origin-to-target distances (before any factor was scored)

Primary is O1->Y1. O2 and Y2 are secondary. Y2 is S2 only and is skipped where the young-donor covariance is singular.

| space | origin | target | distance | Mahalanobis | radius | origin inside Y2 |
| --- | --- | --- | --- | --- | --- | --- |
| S1 | O1 | Y1 | +32.2839 | NA | NA | None |
| S2_k20 | O1 | Y1 | +28.6452 | NA | NA | None |
| S2_k50 | O1 | Y1 | +30.8643 | NA | NA | None |
| S2_k100 | O1 | Y1 | +31.4343 | NA | NA | None |
| S1 | O2 | Y1 | +468.3878 | NA | NA | None |
| S2_k20 | O2 | Y1 | +196.4504 | NA | NA | None |
| S2_k50 | O2 | Y1 | +259.3335 | NA | NA | None |
| S2_k100 | O2 | Y1 | +272.0445 | NA | NA | None |
| S2_k20 | O1 | Y2 | +33.2100 | +1.0554 | +11.8805 | True |
| S2_k50 | O1 | Y2 | +35.1680 | +6.0114 | +59.2987 | True |
| S2_k100 | O1 | Y2 | NA | NA | NA | None |
| S2_k20 | O2 | Y2 | +204.3757 | +14.2055 | +11.8805 | False |
| S2_k50 | O2 | Y2 | +265.1887 | +127.0348 | +59.2987 | False |
| S2_k100 | O2 | Y2 | NA | NA | NA | None |

## Stage 2 — detection limit

The floor is a bar. It is not subtracted from any later number.

### A. Primary permutation floor

Each factor's matched-gene vector was permuted 5000 times (seed 20260914) as one uninterrupted stream, with no checkpoint or resume. The per-factor floor is the 5th percentile. The pooled floor rescales those draws to that setting's median factor norm.

| space | pair | per-factor floor min | median | max | pooled floor |
| --- | --- | --- | --- | --- | --- |
| S1 | O1_Y1 | +42.4155 | +94.9645 | +363.5736 | +94.9656 |
| S2_k20 | O1_Y1 | -0.8636 | -0.3808 | +8.2257 | -0.0750 |
| S2_k50 | O1_Y1 | -0.6454 | +0.4931 | +31.2387 | -0.1110 |
| S2_k100 | O1_Y1 | -0.2893 | +1.9734 | +55.2062 | -0.0747 |
| S1 | O2_Y1 | +3.7869 | +13.3767 | +133.8689 | +13.2262 |
| S2_k20 | O2_Y1 | -1.0420 | -0.7536 | +0.9439 | -0.0799 |
| S2_k50 | O2_Y1 | -0.8929 | -0.3340 | +7.3311 | -0.0789 |
| S2_k100 | O2_Y1 | -1.0906 | -0.8598 | +8.8308 | -0.2129 |
| S2_k20 | O1_Y2 | -1.0660 | -0.8736 | +4.9581 | -0.1190 |
| S2_k50 | O1_Y2 | -0.8156 | -0.2094 | +27.2125 | -0.1953 |
| S2_k20 | O2_Y2 | -1.1783 | -0.8462 | +0.3575 | -0.0878 |
| S2_k50 | O2_Y2 | -0.9273 | -0.4236 | +6.6616 | -0.0909 |

### B. Control floor, computed twice

Primary: all 78 label-defined non-targeting guides. Report-only: the 61 `control=True` guides only. Non-targeting guide pseudo-effects are leave-one-guide-out within each control set. Random splits are 200 disjoint halves of that control set's cells. The report-only floor never enters the primary analysis, Stage 2C, Stage 3 or Stage 5.

| space | pair | control set | guides | cells | NT-guide floor | split floor | stricter |
| --- | --- | --- | --- | --- | --- | --- | --- |
| S1 | O1_Y1 | all78 | 78 | 14013 | +41.8464 | +5.4827 | split |
| S1 | O1_Y1 | control_true61 | 61 | 11834 | +42.0843 | +6.4176 | split |
| S2_k20 | O1_Y1 | all78 | 78 | 14013 | -0.6068 | -0.4372 | nt_guide |
| S2_k20 | O1_Y1 | control_true61 | 61 | 11834 | -0.7303 | -0.3718 | nt_guide |
| S2_k50 | O1_Y1 | all78 | 78 | 14013 | -0.0810 | -0.3263 | split |
| S2_k50 | O1_Y1 | control_true61 | 61 | 11834 | -0.1506 | -0.3144 | split |
| S2_k100 | O1_Y1 | all78 | 78 | 14013 | +0.2455 | -0.3154 | split |
| S2_k100 | O1_Y1 | control_true61 | 61 | 11834 | +0.1283 | -0.3203 | split |
| S1 | O2_Y1 | all78 | 78 | 14013 | +5.0245 | +0.0421 | split |
| S1 | O2_Y1 | control_true61 | 61 | 11834 | +4.8619 | +0.0114 | split |
| S2_k20 | O2_Y1 | all78 | 78 | 14013 | -0.5327 | -0.4390 | nt_guide |
| S2_k20 | O2_Y1 | control_true61 | 61 | 11834 | -0.6598 | -0.6116 | nt_guide |
| S2_k50 | O2_Y1 | all78 | 78 | 14013 | +0.0200 | -0.4483 | split |
| S2_k50 | O2_Y1 | control_true61 | 61 | 11834 | -0.3490 | -0.6380 | split |
| S2_k100 | O2_Y1 | all78 | 78 | 14013 | -0.0805 | -0.4141 | split |
| S2_k100 | O2_Y1 | control_true61 | 61 | 11834 | -0.1548 | -0.5571 | split |
| S2_k20 | O1_Y2 | all78 | 78 | 14013 | +0.2688 | -0.3050 | split |
| S2_k20 | O1_Y2 | control_true61 | 61 | 11834 | +0.0935 | -0.3370 | split |
| S2_k50 | O1_Y2 | all78 | 78 | 14013 | +0.5086 | -0.2790 | split |
| S2_k50 | O1_Y2 | control_true61 | 61 | 11834 | +0.3336 | -0.2924 | split |
| S2_k20 | O2_Y2 | all78 | 78 | 14013 | -0.5042 | -0.4467 | nt_guide |
| S2_k20 | O2_Y2 | control_true61 | 61 | 11834 | -0.5738 | -0.6146 | split |
| S2_k50 | O2_Y2 | all78 | 78 | 14013 | -0.1214 | -0.4623 | split |
| S2_k50 | O2_Y2 | control_true61 | 61 | 11834 | -0.3840 | -0.6252 | split |

Material difference was defined before the numbers: a relative gap |f61 - f78| / |f78| above 0.10 on either floor, or a change in the Stage 2C minimum detectable amplitude computed against the stricter floor.

| space | pair | guide floor 78 | guide floor 61 | rel gap | split floor 78 | split floor 61 | rel gap | MDA 78 | MDA 61 | differs materially |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| S1 | O1_Y1 | +41.8464 | +42.0843 | +0.006 | +5.4827 | +6.4176 | +0.171 | +0.05 | +0.05 | True |
| S2_k20 | O1_Y1 | -0.6068 | -0.7303 | +0.204 | -0.4372 | -0.3718 | +0.149 | +0.10 | +0.10 | True |
| S2_k50 | O1_Y1 | -0.0810 | -0.1506 | +0.860 | -0.3263 | -0.3144 | +0.037 | +0.05 | +0.05 | True |
| S2_k100 | O1_Y1 | +0.2455 | +0.1283 | +0.477 | -0.3154 | -0.3203 | +0.016 | +0.05 | +0.05 | True |
| S1 | O2_Y1 | +5.0245 | +4.8619 | +0.032 | +0.0421 | +0.0114 | +0.728 | +0.05 | +0.05 | True |
| S2_k20 | O2_Y1 | -0.5327 | -0.6598 | +0.239 | -0.4390 | -0.6116 | +0.393 | +0.25 | +0.25 | True |
| S2_k50 | O2_Y1 | +0.0200 | -0.3490 | +18.441 | -0.4483 | -0.6380 | +0.423 | +0.25 | +0.25 | True |
| S2_k100 | O2_Y1 | -0.0805 | -0.1548 | +0.923 | -0.4141 | -0.5571 | +0.345 | +0.10 | +0.25 | True |
| S2_k20 | O1_Y2 | +0.2688 | +0.0935 | +0.652 | -0.3050 | -0.3370 | +0.105 | +0.05 | +0.05 | True |
| S2_k50 | O1_Y2 | +0.5086 | +0.3336 | +0.344 | -0.2790 | -0.2924 | +0.048 | +0.05 | +0.05 | True |
| S2_k20 | O2_Y2 | -0.5042 | -0.5738 | +0.138 | -0.4467 | -0.6146 | +0.376 | +0.25 | +0.25 | True |
| S2_k50 | O2_Y2 | -0.1214 | -0.3840 | +2.163 | -0.4623 | -0.6252 | +0.352 | +0.25 | +0.25 | True |

**The two floors differ materially in 12 of 12 settings: S1 O1_Y1 (guide-floor rel gap +0.006, split-floor rel gap +0.171, MDA +0.05 vs +0.05); S2_k20 O1_Y1 (guide-floor rel gap +0.204, split-floor rel gap +0.149, MDA +0.10 vs +0.10); S2_k50 O1_Y1 (guide-floor rel gap +0.860, split-floor rel gap +0.037, MDA +0.05 vs +0.05); S2_k100 O1_Y1 (guide-floor rel gap +0.477, split-floor rel gap +0.016, MDA +0.05 vs +0.05); S1 O2_Y1 (guide-floor rel gap +0.032, split-floor rel gap +0.728, MDA +0.05 vs +0.05); S2_k20 O2_Y1 (guide-floor rel gap +0.239, split-floor rel gap +0.393, MDA +0.25 vs +0.25); S2_k50 O2_Y1 (guide-floor rel gap +18.441, split-floor rel gap +0.423, MDA +0.25 vs +0.25); S2_k100 O2_Y1 (guide-floor rel gap +0.923, split-floor rel gap +0.345, MDA +0.10 vs +0.25); S2_k20 O1_Y2 (guide-floor rel gap +0.652, split-floor rel gap +0.105, MDA +0.05 vs +0.05); S2_k50 O1_Y2 (guide-floor rel gap +0.344, split-floor rel gap +0.048, MDA +0.05 vs +0.05); S2_k20 O2_Y2 (guide-floor rel gap +0.138, split-floor rel gap +0.376, MDA +0.25 vs +0.25); S2_k50 O2_Y2 (guide-floor rel gap +2.163, split-floor rel gap +0.352, MDA +0.25 vs +0.25).** Stated plainly: which 17 guides are counted as controls changes the measured control floor in those settings. The primary analysis was not changed; it uses all 78 guides, as pre-registered, and the 61-guide floor is reported only.

### C. Dose-response and MDA

The true young-minus-old direction is added to the origin at f = 0.0, 0.05, 0.1, 0.25, 0.5, 1.0, scaled to that setting's median factor norm. A dose clears a floor only when delta is below the floor and negative. The stricter floor here is the primary 78-guide one.

| space | pair | primary | MDA vs pooled | MDA vs stricter | margin | stricter source |
| --- | --- | --- | --- | --- | --- | --- |
| S1 | O1_Y1 | True | +0.05 | +0.05 | +11.6502 | split |
| S2_k20 | O1_Y1 | True | +0.05 | +0.10 | +0.5745 | nt_guide |
| S2_k50 | O1_Y1 | True | +0.05 | +0.05 | +0.5763 | split |
| S2_k100 | O1_Y1 | True | +0.05 | +0.05 | +0.7718 | split |
| S1 | O2_Y1 | False | +0.05 | +0.05 | +0.7832 | split |
| S2_k20 | O2_Y1 | False | +0.05 | +0.25 | +0.3131 | nt_guide |
| S2_k50 | O2_Y1 | False | +0.05 | +0.25 | +0.6022 | split |
| S2_k100 | O2_Y1 | False | +0.05 | +0.10 | +0.0575 | split |
| S2_k20 | O1_Y2 | False | +0.05 | +0.05 | +0.2083 | split |
| S2_k50 | O1_Y2 | False | +0.05 | +0.05 | +0.5113 | split |
| S2_k20 | O2_Y2 | False | +0.05 | +0.25 | +0.3128 | nt_guide |
| S2_k50 | O2_Y2 | False | +0.05 | +0.25 | +0.5655 | split |

## Stage 3 — ranking

Attainable minimum permutation p = 1/(n_perm+1) = 0.000200. Factors tied at that minimum p per setting: S1 O1_Y1=63; S2_k20 O1_Y1=4; S2_k50 O1_Y1=1; S2_k100 O1_Y1=0; S1 O2_Y1=8; S2_k20 O2_Y1=13; S2_k50 O2_Y1=9; S2_k100 O2_Y1=10; S2_k20 O1_Y2=0; S2_k50 O1_Y2=0; S2_k20 O2_Y2=15; S2_k50 O2_Y2=10.

`toward_young` requires delta < 0, delta below that factor's own Stage 2A floor, q <= 0.05, no identity loss, and guide agreement > 0. Counts are not pooled.

| space | pair | primary | delta<0 | below floor | q<=0.05 | identity loss | guides disagree | toward_young | CI invalid |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| S1 | O1_Y1 | True | 0 | 377 | 202 | 0 | 0 | 0 | 1836 |
| S2_k20 | O1_Y1 | True | 108 | 58 | 0 | 0 | 0 | 0 | 1189 |
| S2_k50 | O1_Y1 | True | 18 | 13 | 0 | 0 | 0 | 0 | 1755 |
| S2_k100 | O1_Y1 | True | 5 | 6 | 0 | 0 | 0 | 0 | 1813 |
| S1 | O2_Y1 | False | 0 | 25 | 8 | 0 | 0 | 0 | 1836 |
| S2_k20 | O2_Y1 | False | 96 | 38 | 15 | 0 | 0 | 15 | 333 |
| S2_k50 | O2_Y1 | False | 27 | 18 | 9 | 0 | 0 | 9 | 1360 |
| S2_k100 | O2_Y1 | False | 43 | 18 | 10 | 0 | 0 | 10 | 1329 |
| S2_k20 | O1_Y2 | False | 24 | 3 | 0 | 0 | 0 | 0 | 1546 |
| S2_k50 | O1_Y2 | False | 9 | 3 | 0 | 0 | 0 | 0 | 1802 |
| S2_k20 | O2_Y2 | False | 100 | 43 | 17 | 0 | 0 | 17 | 334 |
| S2_k50 | O2_Y2 | False | 31 | 20 | 10 | 0 | 0 | 10 | 1317 |

The top 20 factors by delta in each setting are in `results/south6/top20.csv`.

## Stage 4 — does delta add anything?

| space | pair | Spearman | median ||d||/||gap|| | median |residual|/||d|| | survivors | recovered by top-n cosine |
| --- | --- | --- | --- | --- | --- | --- |
| S1 | O1_Y1 | -0.169 | +3.8208 | +0.7702 | 0 | NA |
| S2_k20 | O1_Y1 | +0.642 | +0.4124 | +0.1883 | 0 | NA |
| S2_k50 | O1_Y1 | +0.344 | +0.5849 | +0.2672 | 0 | NA |
| S2_k100 | O1_Y1 | +0.217 | +0.6917 | +0.3070 | 0 | NA |
| S1 | O2_Y1 | +0.587 | +0.2634 | +0.1287 | 0 | NA |
| S2_k20 | O2_Y1 | +0.834 | +0.0601 | +0.0254 | 15 | +0.867 |
| S2_k50 | O2_Y1 | +0.755 | +0.0696 | +0.0307 | 9 | +1.000 |
| S2_k100 | O2_Y1 | +0.769 | +0.0799 | +0.0374 | 10 | +0.900 |
| S2_k20 | O1_Y2 | +0.619 | +0.3557 | +0.1555 | 0 | NA |
| S2_k50 | O1_Y2 | +0.484 | +0.5133 | +0.2317 | 0 | NA |
| S2_k20 | O2_Y2 | +0.839 | +0.0578 | +0.0242 | 17 | +0.941 |
| S2_k50 | O2_Y2 | +0.765 | +0.0681 | +0.0299 | 10 | +1.000 |

Across settings the Spearman correlation between the delta ranking and the cosine ranking ranges from -0.169 to +0.839. The median |delta - (-||d|| cos)| / ||d|| ranges from 0.0242 to 0.7702. Where a setting has survivors, top-n cosine recovery ranges from 0.867 to 1.000. Delta and cosine separate (Spearman below 0.8) in: S1 O1_Y1 Q1 (Spearman -0.059, median ratio 3.0970); S1 O1_Y1 Q2 (Spearman +0.110, median ratio 3.6220); S1 O1_Y1 Q3 (Spearman -0.022, median ratio 4.0748); S1 O1_Y1 Q4 (Spearman +0.031, median ratio 4.7617); S2_k20 O1_Y1 Q4 (Spearman +0.779, median ratio 0.5942); S2_k50 O1_Y1 Q4 (Spearman +0.580, median ratio 0.8414); S2_k100 O1_Y1 Q1 (Spearman +0.763, median ratio 0.4753); S2_k100 O1_Y1 Q3 (Spearman +0.749, median ratio 0.8037); S2_k100 O1_Y1 Q4 (Spearman +0.469, median ratio 1.0028); S1 O2_Y1 Q1 (Spearman +0.732, median ratio 0.2135); S1 O2_Y1 Q3 (Spearman +0.786, median ratio 0.2809); S1 O2_Y1 Q4 (Spearman +0.564, median ratio 0.3282); S2_k20 O1_Y2 Q4 (Spearman +0.663, median ratio 0.5125); S2_k50 O1_Y2 Q4 (Spearman +0.502, median ratio 0.7385).

## Stage 5 — combinations

This assumes effects add linearly. No public fibroblast screen can test that assumption.

Unique sets evaluated: 592. Sets reported: 223.

Reported sets:

| space | pair | members | delta | margin over own floor | ratio ||d||/||gap|| | identity flag | hits | leave-one-out delta |
| --- | --- | --- | --- | --- | --- | --- | --- | --- |
| S2_k20 | O1_Y1 | GTF3C2;RELA | -2.7122 | +2.1953 | +0.4854 | False | 1 | GTF3C2:-2.0134;RELA:-1.7703 |
| S2_k20 | O1_Y1 | DMRTA1;GSC2;HEY2 | -2.2402 | +1.7661 | +0.4679 | False | 1 | DMRTA1:-0.6079;GSC2:-1.1577;HEY2:-2.8273 |
| S2_k20 | O1_Y1 | DMRTA1;MAPK8IP1;THAP4 | -1.9209 | +1.8052 | +0.4605 | False | 1 | DMRTA1:-0.6467;MAPK8IP1:-1.4637;THAP4:-2.7255 |
| S2_k20 | O1_Y1 | BNC1;MAPK8IP1;TERF1 | -1.4389 | +1.3123 | +0.5112 | False | 1 | BNC1:+0.4768;MAPK8IP1:-1.2534;TERF1:-2.1914 |
| S2_k20 | O1_Y1 | DMRTA1;MAPK8IP1;VPS72 | -1.0498 | +0.7879 | +0.5254 | False | 1 | DMRTA1:-0.3190;MAPK8IP1:-0.8299;VPS72:-2.7255 |
| S2_k20 | O1_Y1 | DMRTA1;GSC2;MAPK8IP1;ZFP2 | -0.7724 | +0.7677 | +0.5746 | False | 1 | DMRTA1:+0.5445;GSC2:+0.2164;MAPK8IP1:-0.7440;ZFP2:-2.9387 |
| S2_k20 | O1_Y1 | DMRTA1;ERCC8;MAPK8IP1 | -1.1705 | +1.1003 | +0.4430 | False | 1 | DMRTA1:+0.1996;ERCC8:-2.7255;MAPK8IP1:-0.8809 |
| S2_k20 | O1_Y1 | GSC2;IRX4;MAPK8IP1 | -0.8851 | +0.5816 | +0.5025 | False | 1 | GSC2:+0.3377;IRX4:-2.1734;MAPK8IP1:-0.5818 |
| S2_k20 | O1_Y1 | GTF3C2;XBP1 | -0.8174 | +0.7884 | +0.5241 | False | 1 | GTF3C2:-0.4650;XBP1:-1.7703 |
| S2_k50 | O1_Y1 | CRTC1 | -1.1586 | +0.8734 | +0.3153 | False | 1 | CRTC1:+nan |
| S2_k20 | O2_Y1 | HELZ2;JDP2;OSR2;ZEB2;ZNF496 | -10.2371 | +9.4258 | +0.3127 | False | 1 | HELZ2:-8.9144;JDP2:-8.7153;OSR2:-6.9893;ZEB2:-9.5945;ZNF496:-13.1928 |
| S2_k20 | O2_Y1 | HELZ2;JDP2;OSR2;ZEB2;ZNF197 | -10.2667 | +9.7490 | +0.3295 | False | 1 | HELZ2:-9.2928;JDP2:-8.7649;OSR2:-7.1929;ZEB2:-9.8745;ZNF197:-13.1928 |
| S2_k20 | O2_Y1 | FOXF1;HELZ2;JDP2;MEF2A;OSR2 | -6.3070 | +5.4542 | +0.2520 | False | 1 | FOXF1:-4.3531;HELZ2:-4.1031;JDP2:-3.8252;MEF2A:-13.7040;OSR2:-2.3474 |
| S2_k20 | O2_Y1 | FOXF1;HELZ2;JDP2;OSR2;ZKSCAN3 | -6.7007 | +5.7650 | +0.2950 | False | 1 | FOXF1:-5.3957;HELZ2:-4.6598;JDP2:-5.0697;OSR2:-3.2507;ZKSCAN3:-13.7040 |
| S2_k20 | O2_Y1 | HELZ2;ID4;JDP2;OSR2;ZEB2 | -2.5776 | +2.0665 | +0.3071 | False | 1 | HELZ2:-1.3993;ID4:-13.1928;JDP2:-0.5780;OSR2:+0.7797;ZEB2:-1.5871 |
| S2_k20 | O2_Y1 | HELZ2;JDP2;NRF1;OSR2;ZEB2 | -4.3662 | +3.7852 | +0.2995 | False | 1 | HELZ2:-2.8479;JDP2:-2.6189;NRF1:-13.1928;OSR2:-0.9347;ZEB2:-3.2689 |
| S2_k20 | O2_Y1 | HELZ2;HMGB4;JDP2;OSR2;ZEB2 | -7.7937 | +6.9017 | +0.3113 | False | 1 | HELZ2:-6.4926;HMGB4:-13.1928;JDP2:-6.3049;OSR2:-4.4861;ZEB2:-7.0844 |
| S2_k20 | O2_Y1 | HELZ2;JDP2;OSR2;ZEB2;ZNF764 | -11.2499 | +10.3509 | +0.3052 | False | 1 | HELZ2:-9.8359;JDP2:-9.6716;OSR2:-7.8773;ZEB2:-10.4336;ZNF764:-13.1928 |
| S2_k20 | O2_Y1 | HELZ2;JDP2;OSR2;TGIF2;ZEB2 | -8.7466 | +7.9280 | +0.3275 | False | 1 | HELZ2:-7.5617;JDP2:-7.4368;OSR2:-5.7257;TGIF2:-13.1928;ZEB2:-8.3861 |
| S2_k20 | O2_Y1 | HELZ2;JDP2;OSR2;RCOR3;ZEB2 | -7.3832 | +6.5731 | +0.3173 | False | 1 | HELZ2:-6.2392;JDP2:-5.6476;OSR2:-4.2475;RCOR3:-13.1928;ZEB2:-6.7375 |
| S2_k20 | O2_Y1 | FOXF1;HLF;JDP2;OSR2;ZEB2 | -2.9284 | +2.2257 | +0.2888 | False | 1 | FOXF1:-1.8391;HLF:-12.7344;JDP2:-1.5236;OSR2:+0.5365;ZEB2:-1.2843 |
| S2_k20 | O2_Y1 | FOXF1;HELZ2;JDP2;OSR2;SUB1 | -8.1430 | +7.0959 | +0.2411 | False | 1 | FOXF1:-6.4495;HELZ2:-5.5208;JDP2:-5.7610;OSR2:-4.1630;SUB1:-13.7040 |
| S2_k20 | O2_Y1 | FOXF1;JDP2;OSR2;ZEB2;ZNF526 | -7.3509 | +6.5148 | +0.2761 | False | 1 | FOXF1:-6.1120;JDP2:-5.8564;OSR2:-3.6729;ZEB2:-5.3199;ZNF526:-12.7344 |
| S2_k20 | O2_Y1 | GLI4;HELZ2;JDP2;OSR2;ZEB2 | -8.3575 | +7.5195 | +0.2941 | False | 1 | GLI4:-13.1928;HELZ2:-6.9207;JDP2:-6.4528;OSR2:-4.9269;ZEB2:-7.2091 |
| S2_k20 | O2_Y1 | FOXF1;JDP2;OSR2;ZEB2;ZNF683 | -8.3720 | +7.5909 | +0.2723 | False | 1 | FOXF1:-7.1367;JDP2:-6.6210;OSR2:-4.7073;ZEB2:-6.3112;ZNF683:-12.7344 |
| S2_k20 | O2_Y1 | HELZ2;JDP2;OSR2;ZBTB14;ZEB2 | -6.1113 | +5.2228 | +0.3064 | False | 1 | HELZ2:-4.6957;JDP2:-4.5930;OSR2:-2.8165;ZBTB14:-13.1928;ZEB2:-5.2637 |
| S2_k20 | O2_Y1 | HELZ2;JDP2;OSR2;SMARCD1;ZEB2 | -9.4738 | +8.6444 | +0.3022 | False | 1 | HELZ2:-8.1501;JDP2:-7.7200;OSR2:-6.1019;SMARCD1:-13.1928;ZEB2:-8.5590 |
| S2_k20 | O2_Y1 | HELZ2;JDP2;OSR2;ZEB2;ZNF394 | -9.5842 | +8.8904 | +0.3053 | False | 1 | HELZ2:-8.1896;JDP2:-7.9191;OSR2:-6.2489;ZEB2:-8.7583;ZNF394:-13.1928 |
| S2_k20 | O2_Y1 | HELZ2;JDP2;OSR2;POU3F2;ZEB2 | -3.8564 | +3.3802 | +0.3074 | False | 1 | HELZ2:-2.4960;JDP2:-2.1163;OSR2:-0.5276;POU3F2:-13.1928;ZEB2:-2.8677 |
| S2_k20 | O2_Y1 | FOXF1;HELZ2;OSR2;ZEB2;ZSCAN21 | -6.9941 | +6.2774 | +0.2912 | False | 1 | FOXF1:-5.2902;HELZ2:-5.2594;OSR2:-3.5098;ZEB2:-5.6283;ZSCAN21:-13.3970 |
| S2_k20 | O2_Y1 | FOXF1;JDP2;OSR2;TAF6;ZEB2 | -13.1961 | +12.2146 | +0.2800 | False | 1 | FOXF1:-11.9150;JDP2:-11.6869;OSR2:-9.6343;TAF6:-12.7344;ZEB2:-11.3425 |
| S2_k20 | O2_Y1 | FOSB;HELZ2;JDP2;OSR2;ZNF300 | -6.8815 | +6.0844 | +0.2503 | False | 1 | FOSB:-5.0474;HELZ2:-4.3940;JDP2:-4.8339;OSR2:-3.0284;ZNF300:-13.4564 |
| S2_k20 | O2_Y1 | FOXF1;GFI1B;HELZ2;OSR2;ZEB2 | -4.4831 | +3.8813 | +0.2942 | False | 1 | FOXF1:-2.9970;GFI1B:-13.3970;HELZ2:-2.5393;OSR2:-1.1015;ZEB2:-2.9815 |
| S2_k20 | O2_Y1 | FOXF1;HELZ2;OSR2;ZEB2;ZNF3 | -5.2887 | +4.4186 | +0.2713 | False | 1 | FOXF1:-3.6243;HELZ2:-3.1898;OSR2:-1.5871;ZEB2:-3.3380;ZNF3:-13.3970 |
| S2_k20 | O2_Y1 | HELZ2;JDP2;OSR2;ZEB2;ZNF792 | -8.1929 | +7.4942 | +0.2910 | False | 1 | HELZ2:-6.5740;JDP2:-6.2926;OSR2:-4.7854;ZEB2:-6.8585;ZNF792:-13.1928 |
| S2_k20 | O2_Y1 | FOXF1;JDP2;OSR2;ZEB2;ZNF625 | -7.4015 | +6.3328 | +0.2724 | False | 1 | FOXF1:-6.1161;JDP2:-5.9686;OSR2:-3.7234;ZEB2:-5.5045;ZNF625:-12.7344 |
| S2_k20 | O2_Y1 | HELZ2;HNF4A;JDP2;OSR2;ZEB2 | -9.0510 | +8.2922 | +0.3166 | False | 1 | HELZ2:-7.8824;HNF4A:-13.1928;JDP2:-7.4038;OSR2:-5.8660;ZEB2:-8.4503 |
| S2_k20 | O2_Y1 | HELZ2;JDP2;OSR2;POU5F1;ZEB2 | -7.5638 | +6.9154 | +0.3164 | False | 1 | HELZ2:-6.3800;JDP2:-5.9681;OSR2:-4.4056;POU5F1:-13.1928;ZEB2:-6.8679 |
| S2_k20 | O2_Y1 | FOXF1;JDP2;OSR2;ZEB2;ZSCAN5A | -8.5410 | +7.8641 | +0.2809 | False | 1 | FOXF1:-7.2839;JDP2:-6.8622;OSR2:-4.9869;ZEB2:-6.6279;ZSCAN5A:-12.7344 |
| S2_k20 | O2_Y1 | FOXF1;HELZ2;JDP2;OSR2;ZNF397 | -8.2579 | +7.2204 | +0.2428 | False | 1 | FOXF1:-6.4662;HELZ2:-5.8836;JDP2:-5.7771;OSR2:-4.2143;ZNF397:-13.7040 |
| S2_k20 | O2_Y1 | CERS3;HELZ2;JDP2;MAFF;OSR2 | -9.5503 | +8.4912 | +0.2690 | False | 1 | CERS3:-8.1489;HELZ2:-7.6195;JDP2:-7.4683;MAFF:-13.7453;OSR2:-5.7947 |
| S2_k20 | O2_Y1 | EBF4;FOXF1;HELZ2;JDP2;OSR2 | -11.5924 | +10.4736 | +0.2540 | False | 1 | EBF4:-13.7040;FOXF1:-10.2590;HELZ2:-8.5948;JDP2:-9.7667;OSR2:-7.6545 |
| S2_k20 | O2_Y1 | CERS3;HELZ2;JDP2;OSR2;PHOX2A | -9.0876 | +8.0599 | +0.2804 | False | 1 | CERS3:-7.7045;HELZ2:-7.3654;JDP2:-7.1707;OSR2:-5.4748;PHOX2A:-13.7453 |
| S2_k20 | O2_Y1 | HELZ2;JDP2;OSR2;ZEB2;ZNF561 | -9.3914 | +8.8381 | +0.3243 | False | 1 | HELZ2:-8.3394;JDP2:-7.8176;OSR2:-6.2469;ZEB2:-8.9735;ZNF561:-13.1928 |
| S2_k20 | O2_Y1 | FOXF1;HELZ2;JDP2;OSR2;SRF | -4.2543 | +3.1821 | +0.2629 | False | 1 | FOXF1:-2.9854;HELZ2:-1.5752;JDP2:-2.3929;OSR2:-0.4801;SRF:-13.7040 |
| S2_k20 | O2_Y1 | HELZ2;JDP2;KDM5B;OSR2;ZEB2 | -9.4349 | +8.6299 | +0.3315 | False | 1 | HELZ2:-8.4370;JDP2:-8.0641;KDM5B:-13.1928;OSR2:-6.4062;ZEB2:-9.2115 |
| S2_k20 | O2_Y1 | HELZ2;ISL1;JDP2;OSR2;ZEB2 | -9.1038 | +8.3596 | +0.3267 | False | 1 | HELZ2:-7.9025;ISL1:-13.1928;JDP2:-7.7760;OSR2:-6.0498;ZEB2:-8.7440 |
| S2_k20 | O2_Y1 | ELK3;FOXF1;HELZ2;OSR2;ZEB2 | -9.7469 | +8.7276 | +0.2855 | False | 1 | ELK3:-13.3970;FOXF1:-8.0269;HELZ2:-7.9418;OSR2:-6.2333;ZEB2:-8.2784 |
| S2_k20 | O2_Y1 | HELZ2;JDP2;OSR2;ZEB2;ZNF146 | -8.7005 | +7.9017 | +0.3049 | False | 1 | HELZ2:-7.3876;JDP2:-7.0237;OSR2:-5.3603;ZEB2:-7.8737;ZNF146:-13.1928 |
| S2_k20 | O2_Y1 | FOXF1;HELZ2;OSR2;THRA;ZEB2 | -10.4394 | +9.5113 | +0.2890 | False | 1 | FOXF1:-8.3851;HELZ2:-8.9642;OSR2:-6.8959;THRA:-13.3970;ZEB2:-9.0237 |
| S2_k20 | O2_Y1 | FOXF1;HELZ2;MESP1;OSR2;ZEB2 | -11.1578 | +10.1180 | +0.2772 | False | 1 | FOXF1:-9.3740;HELZ2:-9.2105;MESP1:-13.3970;OSR2:-7.5167;ZEB2:-9.4888 |
| S2_k20 | O2_Y1 | ATF1;HELZ2;JDP2;OSR2;ZEB2 | -8.0364 | +7.1447 | +0.3087 | False | 1 | ATF1:-13.1928;HELZ2:-6.8302;JDP2:-6.4118;OSR2:-4.7027;ZEB2:-7.2238 |
| S2_k20 | O2_Y1 | HELZ2;JDP2;OSR2;SMARCA5;ZEB2 | -2.9468 | +2.4042 | +0.3396 | False | 1 | HELZ2:-1.9803;JDP2:-1.6789;OSR2:+0.0335;SMARCA5:-13.1928;ZEB2:-2.8199 |
| S2_k20 | O2_Y1 | FOXA1;HELZ2;JDP2;OSR2;ZEB2 | -2.2978 | +1.6320 | +0.3000 | False | 1 | FOXA1:-13.1928;HELZ2:-0.7141;JDP2:-0.6183;OSR2:+1.1041;ZEB2:-1.1475 |
| S2_k20 | O2_Y1 | CNOT3;FOXF1;JDP2;OSR2;ZEB2 | -9.3580 | +8.4037 | +0.2663 | False | 1 | CNOT3:-12.7344;FOXF1:-8.0713;JDP2:-7.7049;OSR2:-5.6458;ZEB2:-7.1697 |
| S2_k20 | O2_Y1 | HELZ2;JDP2;OSR2;RXRA;ZEB2 | -13.1421 | +12.1494 | +0.2916 | False | 1 | HELZ2:-11.7037;JDP2:-11.1829;OSR2:-9.6027;RXRA:-13.1928;ZEB2:-11.9395 |
| S2_k20 | O2_Y1 | HELZ2;JDP2;OSR2;ZEB2;ZNF492 | -9.5272 | +8.6424 | +0.3002 | False | 1 | HELZ2:-8.2192;JDP2:-7.7454;OSR2:-6.1087;ZEB2:-8.5453;ZNF492:-13.1928 |
| S2_k50 | O2_Y1 | FOXF1;HELZ2;JDP2;OSR2;ZNF84 | -6.0004 | +6.3170 | +0.2168 | False | 1 | FOXF1:-4.2692;HELZ2:-4.2314;JDP2:-4.8349;OSR2:-4.6164;ZNF84:-9.2691 |
| S2_k50 | O2_Y1 | FOXF1;HELZ2;JDP2;NFKB2;OSR2 | -4.1942 | +4.3708 | +0.1994 | False | 1 | FOXF1:-2.2966;HELZ2:-2.3560;JDP2:-2.6429;NFKB2:-9.2691;OSR2:-2.5882 |
| S2_k50 | O2_Y1 | BLZF1;FOSB;FOXF1;HELZ2;OSR2 | -0.6090 | +1.1029 | +0.2034 | False | 1 | BLZF1:-9.3813;FOSB:+0.8229;FOXF1:+1.1319;HELZ2:+1.6170;OSR2:+1.0167 |
| S2_k50 | O2_Y1 | FOXF1;HELZ2;HIC1;KLF5;OSR2 | -7.0890 | +7.2143 | +0.2232 | False | 1 | FOXF1:-4.9569;HELZ2:-5.4810;HIC1:-6.6494;KLF5:-8.9138;OSR2:-5.6980 |
| S2_k50 | O2_Y1 | FOXF1;HELZ2;JDP2;NFKBIL1;OSR2 | -1.7733 | +2.3049 | +0.2231 | False | 1 | FOXF1:-0.1752;HELZ2:+0.0490;JDP2:-0.7230;NFKBIL1:-9.2691;OSR2:-0.4492 |
| S2_k50 | O2_Y1 | CLOCK;FOXF1;HELZ2;HIC1;OSR2 | -7.6873 | +7.3306 | +0.1911 | False | 1 | CLOCK:-8.9138;FOXF1:-5.8037;HELZ2:-5.2990;HIC1:-6.9738;OSR2:-5.9515 |
| S2_k50 | O2_Y1 | FOXF1;HELZ2;JDP2;OSR2;ZNF319 | -2.3308 | +2.5419 | +0.2121 | False | 1 | FOXF1:-0.8131;HELZ2:-0.2097;JDP2:-1.2615;OSR2:-0.8066;ZNF319:-9.2691 |
| S2_k50 | O2_Y1 | FOXF1;HELZ2;JDP2;OSR2;ZNF658 | -1.3648 | +2.0698 | +0.2353 | False | 1 | FOXF1:+0.1344;HELZ2:+0.4662;JDP2:-0.6683;OSR2:-0.1278;ZNF658:-9.2691 |
| S2_k50 | O2_Y1 | CDK9;FOXF1;HELZ2;JDP2;OSR2 | -1.8860 | +2.0978 | +0.2252 | False | 1 | CDK9:-9.2691;FOXF1:-0.3426;HELZ2:-0.1183;JDP2:-0.9328;OSR2:-0.5940 |
| S2_k50 | O2_Y1 | FOXF1;HELZ2;JDP2;OSR2;ZFP1 | -1.3867 | +2.1488 | +0.2339 | False | 1 | FOXF1:+0.2427;HELZ2:+0.1480;JDP2:-0.4010;OSR2:-0.1700;ZFP1:-9.2691 |
| S2_k50 | O2_Y1 | FOXF1;HELZ2;HIC1;OSR2;TCF7 | -4.6587 | +4.6986 | +0.2004 | False | 1 | FOXF1:-2.7511;HELZ2:-2.4857;HIC1:-4.0427;OSR2:-3.0429;TCF7:-8.9138 |
| S2_k50 | O2_Y1 | FOXF1;GATA1;HELZ2;OSR2;ZNF699 | -0.2791 | +1.1121 | +0.2231 | False | 1 | FOXF1:+1.1924;GATA1:+0.8444;HELZ2:+1.8706;OSR2:+1.0555;ZNF699:-9.1793 |
| S2_k50 | O2_Y1 | FOSB;FOXF1;HELZ2;OSR2;TSHZ2 | -1.7312 | +2.2225 | +0.2060 | False | 1 | FOSB:-0.4837;FOXF1:-0.1645;HELZ2:+0.6614;OSR2:-0.1566;TSHZ2:-9.3813 |
| S2_k50 | O2_Y1 | FOSB;FOXF1;HELZ2;OSR2;SOX8 | -4.2726 | +4.5228 | +0.2101 | False | 1 | FOSB:-2.7477;FOXF1:-2.5213;HELZ2:-2.3775;OSR2:-2.8191;SOX8:-9.3813 |
| S2_k50 | O2_Y1 | FOXF1;GATA1;HELZ2;OSR2;TAF10 | -4.2848 | +4.5711 | +0.2214 | False | 1 | FOXF1:-2.7629;GATA1:-3.3322;HELZ2:-2.1465;OSR2:-2.9831;TAF10:-9.1793 |
| S2_k50 | O2_Y1 | FOSB;FOXF1;HELZ2;HEY2;OSR2 | -8.6253 | +8.4259 | +0.2125 | False | 1 | FOSB:-7.6833;FOXF1:-7.2095;HELZ2:-6.2216;HEY2:-9.3813;OSR2:-7.0490 |
| S2_k50 | O2_Y1 | FOSB;FOXF1;HELZ2;OSR2;ZNF799 | -3.6479 | +4.0693 | +0.2015 | False | 1 | FOSB:-2.1599;FOXF1:-1.9906;HELZ2:-1.4493;OSR2:-2.0018;ZNF799:-9.3813 |
| S2_k50 | O2_Y1 | FOXF1;HELZ2;HIC1;OSR2;ZSCAN4 | -6.9455 | +6.8173 | +0.1831 | False | 1 | FOXF1:-5.0099;HELZ2:-4.4260;HIC1:-6.1764;OSR2:-5.1496;ZSCAN4:-8.9138 |
| S2_k50 | O2_Y1 | FOSB;FOXF1;HELZ2;OSR2;ZNF254 | -5.6316 | +5.9970 | +0.2265 | False | 1 | FOSB:-4.6934;FOXF1:-4.3220;HELZ2:-3.4712;OSR2:-4.3671;ZNF254:-9.3813 |
| S2_k50 | O2_Y1 | E4F1;FOXF1;GATA1;HELZ2;OSR2 | -0.7210 | +1.5285 | +0.2317 | False | 1 | E4F1:-9.1793;FOXF1:+0.8069;GATA1:+0.4266;HELZ2:+1.0849;OSR2:+0.5444 |
| S2_k50 | O2_Y1 | FOXF1;GATA1;GMEB2;HELZ2;OSR2 | -1.9754 | +2.5092 | +0.2161 | False | 1 | FOXF1:-0.3985;GATA1:-0.7327;GMEB2:-9.1793;HELZ2:+0.0263;OSR2:-0.5237 |
| S2_k50 | O2_Y1 | FOSB;FOXF1;GSC;HELZ2;OSR2 | -2.9295 | +3.3260 | +0.2116 | False | 1 | FOSB:-1.5001;FOXF1:-1.2980;GSC:-9.3813;HELZ2:-0.9501;OSR2:-1.4203 |
| S2_k50 | O2_Y1 | FOXF1;HELZ2;JDP2;OSR2;ZNF491 | -2.8994 | +3.4283 | +0.2354 | False | 1 | FOXF1:-1.3765;HELZ2:-1.2163;JDP2:-2.0926;OSR2:-1.7913;ZNF491:-9.2691 |
| S2_k50 | O2_Y1 | FOXF1;HELZ2;JDP2;OSR2;ZNF777 | -1.1976 | +1.5616 | +0.2014 | False | 1 | FOXF1:+0.4314;HELZ2:+1.0779;JDP2:+0.1730;OSR2:+0.4032;ZNF777:-9.2691 |
| S2_k50 | O2_Y1 | FOSB;FOXF1;HELZ2;LIN9;OSR2 | -4.9161 | +4.9536 | +0.2045 | False | 1 | FOSB:-3.8152;FOXF1:-3.4408;HELZ2:-2.4249;LIN9:-9.3813;OSR2:-3.3256 |
| S2_k50 | O2_Y1 | FOXF1;HELZ2;HOXA3;JDP2;OSR2 | -0.9157 | +1.5180 | +0.2186 | False | 2 | FOXF1:+0.7770;HELZ2:+0.9041;HOXA3:-9.2691;JDP2:+0.3018;OSR2:+0.5027 |
| S2_k50 | O2_Y1 | DNMT1;FOSB;FOXF1;HELZ2;OSR2 | -8.3532 | +8.2813 | +0.2073 | False | 1 | DNMT1:-9.3813;FOSB:-7.1929;FOXF1:-6.8696;HELZ2:-5.9394;OSR2:-6.8165 |
| S2_k50 | O2_Y1 | FOXF1;HELZ2;JDP2;MLLT10;OSR2 | -4.9880 | +5.0093 | +0.2202 | False | 1 | FOXF1:-3.3915;HELZ2:-3.1595;JDP2:-3.9459;MLLT10:-9.2691;OSR2:-3.6678 |
| S2_k50 | O2_Y1 | FOSB;FOXF1;HELZ2;HOXA5;OSR2 | -0.4310 | +1.1797 | +0.2265 | False | 1 | FOSB:+0.9709;FOXF1:+1.1715;HELZ2:+1.3160;HOXA5:-9.3813;OSR2:+0.9061 |
| S2_k50 | O2_Y1 | FOXF1;HELZ2;JDP2;OSR2;ZNF665 | -4.3702 | +4.4094 | +0.2088 | False | 1 | FOXF1:-2.7116;HELZ2:-2.4724;JDP2:-3.1987;OSR2:-2.8352;ZNF665:-9.2691 |
| S2_k50 | O2_Y1 | FOXF1;HELZ2;JDP2;OSR2;RCOR1 | -5.5349 | +5.3454 | +0.2008 | False | 1 | FOXF1:-3.8933;HELZ2:-3.3679;JDP2:-4.2091;OSR2:-3.9517;RCOR1:-9.2691 |
| S2_k50 | O2_Y1 | FOXF1;HELZ2;JDP2;OSR2;ZNF93 | -5.1264 | +5.1587 | +0.2175 | False | 1 | FOXF1:-3.6135;HELZ2:-2.9423;JDP2:-4.3524;OSR2:-3.6854;ZNF93:-9.2691 |
| S2_k50 | O2_Y1 | BHLHE40;FOXF1;HELZ2;JDP2;OSR2 | -2.4432 | +2.7821 | +0.2239 | False | 1 | BHLHE40:-9.2691;FOXF1:-0.8004;HELZ2:-0.8710;JDP2:-1.1975;OSR2:-1.0814 |
| S2_k50 | O2_Y1 | EPAS1;FOXF1;GATA1;HELZ2;OSR2 | -2.9485 | +2.9149 | +0.1964 | False | 1 | EPAS1:-9.1793;FOXF1:-1.5177;GATA1:-1.7164;HELZ2:-0.2153;OSR2:-1.1857 |
| S2_k50 | O2_Y1 | AHRR;FOSB;FOXF1;HELZ2;OSR2 | -6.0588 | +6.0229 | +0.1925 | False | 1 | AHRR:-9.3813;FOSB:-4.6357;FOXF1:-4.4698;HELZ2:-3.4783;OSR2:-4.3346 |
| S2_k100 | O2_Y1 | AIRE;FOSB;FOXF1;HELZ2;OSR2 | -6.9021 | +6.0749 | +0.2361 | False | 1 | AIRE:-10.5032;FOSB:-6.0763;FOXF1:-4.8571;HELZ2:-4.2990;OSR2:-5.8947 |
| S2_k100 | O2_Y1 | FOXF1;HELZ2;OSR2;REST;SOX10 | -0.4440 | +0.2831 | +0.2431 | False | 1 | FOXF1:+1.6933;HELZ2:+1.7434;OSR2:+0.5429;REST:+0.0233;SOX10:-10.1154 |
| S2_k100 | O2_Y1 | FOXF1;GATA1;HELZ2;IKZF3;OSR2 | -1.3106 | +0.6394 | +0.2331 | False | 1 | FOXF1:+0.6842;GATA1:+0.0170;HELZ2:+0.9927;IKZF3:-10.3939;OSR2:-0.1791 |
| S2_k100 | O2_Y1 | FOXF1;GATA1;HELZ2;OSR2;SETDB2 | -1.2298 | +1.2052 | +0.2641 | False | 1 | FOXF1:+0.4820;GATA1:-0.4106;HELZ2:+0.8678;OSR2:-0.5853;SETDB2:-10.3939 |
| S2_k100 | O2_Y1 | FOXF1;GATA1;HELZ2;OSR2;ZEB1 | -3.7450 | +3.0772 | +0.2270 | False | 1 | FOXF1:-1.6352;GATA1:-2.3349;HELZ2:-1.3727;OSR2:-2.6071;ZEB1:-10.3939 |
| S2_k100 | O2_Y1 | DPRX;FOXF1;GATA1;HELZ2;OSR2 | -6.4710 | +5.6934 | +0.2307 | False | 1 | DPRX:-10.3939;FOXF1:-4.3196;GATA1:-4.9719;HELZ2:-4.3652;OSR2:-5.3701 |
| S2_k100 | O2_Y1 | FOXF1;GATA1;HELZ2;OSR2;ZNF221 | -2.0031 | +1.9475 | +0.2673 | False | 1 | FOXF1:-0.2651;GATA1:-0.9360;HELZ2:-0.3079;OSR2:-1.4105;ZNF221:-10.3939 |
| S2_k100 | O2_Y1 | FOXF1;GATA1;HELZ2;OSR2;ZNF611 | -5.3677 | +5.0379 | +0.2505 | False | 1 | FOXF1:-3.4006;GATA1:-4.4274;HELZ2:-3.2877;OSR2:-4.5938;ZNF611:-10.3939 |
| S2_k100 | O2_Y1 | FOXF1;GATA1;HELZ2;HIF1A;OSR2 | -3.1738 | +2.3143 | +0.2259 | False | 1 | FOXF1:-1.0735;GATA1:-1.8745;HELZ2:-0.7908;HIF1A:-10.3939;OSR2:-1.9483 |
| S2_k100 | O2_Y1 | FOXF1;GATA1;HELZ2;OSR2;ZNF823 | -1.1438 | +1.0252 | +0.2550 | False | 1 | FOXF1:+1.0130;GATA1:+0.2671;HELZ2:+0.5603;OSR2:-0.3865;ZNF823:-10.3939 |
| S2_k100 | O2_Y1 | AP2B1;FOXF1;GATA1;HELZ2;OSR2 | -2.1618 | +1.6606 | +0.2357 | False | 1 | AP2B1:-10.3939;FOXF1:-0.1806;GATA1:-1.0998;HELZ2:+0.3633;OSR2:-1.0332 |
| S2_k100 | O2_Y1 | FOXF1;GATA1;HELZ2;MECOM;OSR2 | -0.7256 | +0.7635 | +0.2690 | False | 1 | FOXF1:+1.0235;GATA1:+0.0808;HELZ2:+1.1612;MECOM:-10.3939;OSR2:-0.1698 |
| S2_k100 | O2_Y1 | AFF3;FOXF1;HELZ2;OSR2;REST | -3.9599 | +3.4730 | +0.2416 | False | 1 | AFF3:-10.1154;FOXF1:-1.6609;HELZ2:-1.9008;OSR2:-3.1342;REST:-3.4708 |
| S2_k100 | O2_Y1 | FOXF1;GATA1;HELZ2;OSR2;TUB | -5.2876 | +5.1780 | +0.2673 | False | 1 | FOXF1:-3.4105;GATA1:-4.6333;HELZ2:-3.3191;OSR2:-4.8063;TUB:-10.3939 |
| S2_k100 | O2_Y1 | FOXF1;GATA1;HELZ2;HOXB6;OSR2 | -0.9833 | +1.0645 | +0.2741 | False | 1 | FOXF1:+0.6050;GATA1:-0.3859;HELZ2:+1.0789;HOXB6:-10.3939;OSR2:-0.3984 |
| S2_k100 | O2_Y1 | CEBPA;FOXF1;HELZ2;JDP2;OSR2 | -0.5827 | +0.8424 | +0.2686 | False | 1 | CEBPA:-10.2792;FOXF1:+1.5719;HELZ2:+1.0019;JDP2:+0.4055;OSR2:-0.0469 |
| S2_k100 | O2_Y1 | FOXF1;GATA1;HELZ2;HSF1;OSR2 | -6.7063 | +6.1837 | +0.2499 | False | 1 | FOXF1:-4.7134;GATA1:-5.7393;HELZ2:-4.6894;HSF1:-10.3939;OSR2:-5.9636 |
| S2_k100 | O2_Y1 | ERCC2;FOSB;FOXF1;HELZ2;OSR2 | -0.6682 | +0.5169 | +0.2493 | False | 1 | ERCC2:-10.5032;FOSB:+0.2829;FOXF1:+1.3083;HELZ2:+1.5641;OSR2:+0.2599 |
| S2_k100 | O2_Y1 | FOXF1;GATA1;HELZ2;OSR2;TEF | -5.7654 | +5.0458 | +0.2396 | False | 1 | FOXF1:-3.8180;GATA1:-4.8015;HELZ2:-3.5364;OSR2:-4.7109;TEF:-10.3939 |
| S2_k100 | O2_Y1 | FOSB;FOXF1;HELZ2;OSR2;ZNF224 | -3.4553 | +3.4945 | +0.2642 | False | 1 | FOSB:-2.4858;FOXF1:-1.5671;HELZ2:-1.6531;OSR2:-2.8251;ZNF224:-10.5032 |
| S2_k100 | O2_Y1 | FOXF1;GATA1;HELZ2;OSR2;SOX13 | -0.9271 | +1.0671 | +0.2701 | False | 1 | FOXF1:+0.7412;GATA1:-0.2859;HELZ2:+1.2389;OSR2:-0.2398;SOX13:-10.3939 |
| S2_k100 | O2_Y1 | FOXF1;GATA1;HELZ2;OSR2;ZNF396 | -5.3870 | +4.9476 | +0.2550 | False | 1 | FOXF1:-3.5630;GATA1:-4.5498;HELZ2:-3.3710;OSR2:-4.6671;ZNF396:-10.3939 |
| S2_k100 | O2_Y1 | FOXF1;GATA1;HELZ2;OSR2;PMS1 | -0.1943 | +0.2087 | +0.2734 | False | 1 | FOXF1:+1.2907;GATA1:+0.4279;HELZ2:+1.8902;OSR2:+0.3566;PMS1:-10.3939 |
| S2_k100 | O2_Y1 | FOXF1;HELZ2;JDP2;OSR2;ZNF852 | -3.6429 | +3.0969 | +0.2487 | False | 1 | FOXF1:-1.6182;HELZ2:-1.5129;JDP2:-3.0654;OSR2:-2.8237;ZNF852:-10.2792 |
| S2_k100 | O2_Y1 | FOXF1;GATA1;HELZ2;OSR2;THAP5 | -3.7471 | +3.2622 | +0.2428 | False | 1 | FOXF1:-1.7356;GATA1:-2.5683;HELZ2:-1.4704;OSR2:-2.9382;THAP5:-10.3939 |
| S2_k100 | O2_Y1 | FOXF1;GATA1;HELZ2;NKX2-3;OSR2 | -3.8636 | +3.2630 | +0.2355 | False | 1 | FOXF1:-2.0218;GATA1:-2.8108;HELZ2:-1.3497;NKX2-3:-10.3939;OSR2:-2.7267 |
| S2_k100 | O2_Y1 | FEZF1;FOXF1;GATA1;HELZ2;OSR2 | -3.5469 | +2.7845 | +0.2306 | False | 1 | FEZF1:-10.3939;FOXF1:-1.4119;GATA1:-2.4354;HELZ2:-1.1181;OSR2:-2.4294 |
| S2_k100 | O2_Y1 | FOXF1;GATA1;HELZ2;OSR2;ZFX | -4.5445 | +4.1144 | +0.2466 | False | 1 | FOXF1:-2.6181;GATA1:-3.5125;HELZ2:-2.3773;OSR2:-3.6114;ZFX:-10.3939 |
| S2_k100 | O2_Y1 | FOXF1;GATA1;HELZ2;OSR2;ZNF138 | -1.0195 | +0.9131 | +0.2613 | False | 1 | FOXF1:+0.9537;GATA1:-0.2389;HELZ2:+0.8733;OSR2:-0.4224;ZNF138:-10.3939 |
| S2_k100 | O2_Y1 | BNC1;FOXF1;GATA1;HELZ2;OSR2 | -5.0821 | +4.3039 | +0.2245 | False | 1 | BNC1:-10.3939;FOXF1:-3.0337;GATA1:-3.5967;HELZ2:-2.5907;OSR2:-3.8772 |
| S2_k100 | O2_Y1 | FOXF1;GATA1;HELZ2;OSR2;ZNF639 | -3.0377 | +2.9242 | +0.2550 | False | 1 | FOXF1:-1.1854;GATA1:-1.9455;HELZ2:-1.0268;OSR2:-2.2713;ZNF639:-10.3939 |
| S2_k100 | O2_Y1 | FOSB;FOXF1;HELZ2;OSR2;ZNF23 | -7.1048 | +6.8044 | +0.2542 | False | 1 | FOSB:-6.3838;FOXF1:-5.1619;HELZ2:-4.9834;OSR2:-6.2705;ZNF23:-10.5032 |
| S2_k100 | O2_Y1 | FOXF1;GATA1;HELZ2;INSM1;OSR2 | -5.2291 | +4.5193 | +0.2330 | False | 1 | FOXF1:-3.1192;GATA1:-3.9097;HELZ2:-2.9386;INSM1:-10.3939;OSR2:-4.2067 |
| S2_k100 | O2_Y1 | CUX1;FOXF1;GATA1;HELZ2;OSR2 | -4.0036 | +2.9919 | +0.2107 | False | 1 | CUX1:-10.3939;FOXF1:-1.8259;GATA1:-2.3697;HELZ2:-1.3973;OSR2:-2.5667 |
| S2_k100 | O2_Y1 | DEAF1;FOXF1;GATA1;HELZ2;OSR2 | -3.8535 | +3.2521 | +0.2530 | False | 1 | DEAF1:-10.3939;FOXF1:-2.0194;GATA1:-3.0961;HELZ2:-1.6771;OSR2:-3.0791 |
| S2_k100 | O2_Y1 | FOXF1;GATA1;HELZ2;OSR2;ZNF208 | -0.9316 | +1.3102 | +0.2848 | False | 1 | FOXF1:+0.9080;GATA1:-0.2717;HELZ2:+0.6799;OSR2:-0.5797;ZNF208:-10.3939 |
| S2_k100 | O2_Y1 | FOXF1;GATA1;HELZ2;ISL1;OSR2 | -4.4415 | +3.9912 | +0.2545 | False | 1 | FOXF1:-2.6106;GATA1:-3.6201;HELZ2:-2.3134;ISL1:-10.3939;OSR2:-3.7780 |
| S2_k100 | O2_Y1 | FOSB;FOXF1;HELZ2;HOXA7;OSR2 | -4.1184 | +3.6769 | +0.2407 | False | 1 | FOSB:-3.2017;FOXF1:-2.1107;HELZ2:-1.7721;HOXA7:-10.5032;OSR2:-3.1468 |
| S2_k100 | O2_Y1 | FOXF1;GATA1;HELZ2;OSR2;ZSCAN25 | -2.8602 | +2.2738 | +0.2438 | False | 1 | FOXF1:-0.8377;GATA1:-1.8488;HELZ2:-0.7519;OSR2:-1.9937;ZSCAN25:-10.3939 |
| S2_k100 | O2_Y1 | ETV3L;FOXF1;HELZ2;OSR2;REST | -10.1355 | +9.0799 | +0.2488 | False | 1 | ETV3L:-10.1154;FOXF1:-8.0288;HELZ2:-8.0741;OSR2:-9.3698;REST:-9.7861 |
| S2_k20 | O1_Y2 | AFF1;MAFA;OSR2 | -1.3405 | +0.3030 | +0.4773 | False | 1 | AFF1:-2.0453;MAFA:-0.8494;OSR2:+0.3282 |
| S2_k20 | O2_Y2 | HELZ2;JDP2;NEUROD1;OSR2;TBXT | -8.6541 | +7.4867 | +0.2740 | False | 1 | HELZ2:-6.9372;JDP2:-6.0402;NEUROD1:-6.6937;OSR2:-4.6631;TBXT:-15.0583 |
| S2_k20 | O2_Y2 | HELZ2;JDP2;OSR2;PHOX2B;ZEB2 | -11.2345 | +10.2383 | +0.2939 | False | 1 | HELZ2:-9.7023;JDP2:-8.8888;OSR2:-7.5367;PHOX2B:-14.9046;ZEB2:-9.6586 |
| S2_k20 | O2_Y2 | FOXF1;HOXA10;JDP2;OSR2;ZEB2 | -11.9791 | +10.6236 | +0.2594 | False | 1 | FOXF1:-10.0895;HOXA10:-14.8646;JDP2:-9.5607;OSR2:-7.9228;ZEB2:-9.3338 |
| S2_k20 | O2_Y2 | FOXF1;JDP2;OSR2;ZEB2;ZNF274 | -12.9946 | +11.6717 | +0.2619 | False | 1 | FOXF1:-11.3146;JDP2:-10.7050;OSR2:-8.9107;ZEB2:-10.2284;ZNF274:-14.8646 |
| S2_k20 | O2_Y2 | FOXF1;JDP2;OSR2;ZEB2;ZNF502 | -11.7571 | +10.5550 | +0.2547 | False | 2 | FOXF1:-9.9860;JDP2:-9.2313;OSR2:-7.6827;ZEB2:-8.8754;ZNF502:-14.8646 |
| S2_k20 | O2_Y2 | HELZ2;JDP2;MYF5;NEUROD1;OSR2 | -5.6990 | +4.5970 | +0.2784 | False | 1 | HELZ2:-4.0589;JDP2:-3.0292;MYF5:-15.0583;NEUROD1:-3.7390;OSR2:-1.8197 |
| S2_k20 | O2_Y2 | HELZ2;HOXB1;JDP2;OSR2;ZEB2 | -13.3191 | +12.1886 | +0.2938 | False | 1 | HELZ2:-11.7284;HOXB1:-14.9046;JDP2:-11.0259;OSR2:-9.6112;ZEB2:-11.7417 |
| S2_k20 | O2_Y2 | HELZ2;JDP2;NEUROD1;OSR2;ZNF667 | -6.9438 | +5.7254 | +0.2827 | False | 1 | HELZ2:-5.1761;JDP2:-4.4770;NEUROD1:-5.0901;OSR2:-3.1443;ZNF667:-15.0583 |
| S2_k20 | O2_Y2 | FOXF1;JDP2;OSR2;SMAD4;ZEB2 | -8.6942 | +7.4677 | +0.2436 | False | 1 | FOXF1:-6.9724;JDP2:-6.0788;OSR2:-4.4962;SMAD4:-14.8646;ZEB2:-5.5897 |
| S2_k20 | O2_Y2 | HELZ2;JDP2;OSR2;ZEB2;ZNF260 | -4.4027 | +3.5140 | +0.2893 | False | 1 | HELZ2:-2.7058;JDP2:-1.9807;OSR2:-0.7004;ZEB2:-2.7391;ZNF260:-14.9046 |
| S2_k20 | O2_Y2 | FOXF1;FUS;JDP2;OSR2;ZEB2 | -2.6043 | +1.6796 | +0.2731 | False | 1 | FOXF1:-0.9870;FUS:-14.8646;JDP2:-0.3798;OSR2:+1.3326;ZEB2:-0.1632 |
| S2_k20 | O2_Y2 | HELZ2;JDP2;NFIL3;OSR2;ZEB2 | -14.3774 | +13.1222 | +0.3167 | False | 1 | HELZ2:-12.9507;JDP2:-12.4029;NFIL3:-14.9046;OSR2:-10.9365;ZEB2:-13.3587 |
| S2_k20 | O2_Y2 | HELZ2;JDP2;NEUROD1;OSR2;ZNF672 | -4.4210 | +3.2912 | +0.2708 | False | 1 | HELZ2:-2.4619;JDP2:-1.9195;NEUROD1:-2.3469;OSR2:-0.4674;ZNF672:-15.0583 |
| S2_k20 | O2_Y2 | FOXF1;JDP2;OSR2;TBX22;ZEB2 | -12.4649 | +11.1212 | +0.2771 | False | 1 | FOXF1:-11.0340;JDP2:-10.5122;OSR2:-8.5549;TBX22:-14.8646;ZEB2:-9.8975 |
| S2_k20 | O2_Y2 | HELZ2;JDP2;NEUROD1;OSR2;ZNF700 | -10.8924 | +9.7596 | +0.2763 | False | 1 | HELZ2:-9.1409;JDP2:-8.3371;NEUROD1:-8.9307;OSR2:-6.9030;ZNF700:-15.0583 |
| S2_k20 | O2_Y2 | FOXF1;JDP2;NFATC3;OSR2;ZEB2 | -10.0650 | +8.8539 | +0.2785 | False | 1 | FOXF1:-8.5396;JDP2:-8.0421;NFATC3:-14.8646;OSR2:-6.2027;ZEB2:-7.7087 |
| S2_k20 | O2_Y2 | BACH1;HELZ2;JDP2;OSR2;ZEB2 | -6.3154 | +5.2111 | +0.2727 | False | 1 | BACH1:-14.9046;HELZ2:-4.5977;JDP2:-3.6147;OSR2:-2.4398;ZEB2:-4.2105 |
| S2_k20 | O2_Y2 | FOXF1;JDP2;OSR2;SMARCD1;ZEB2 | -11.2515 | +10.0415 | +0.2585 | False | 1 | FOXF1:-9.6093;JDP2:-8.9089;OSR2:-7.1968;SMARCD1:-14.8646;ZEB2:-8.4355 |
| S2_k20 | O2_Y2 | ERCC6;HELZ2;JDP2;OSR2;ZEB2 | -6.6726 | +5.6334 | +0.2884 | False | 1 | ERCC6:-14.9046;HELZ2:-4.9691;JDP2:-4.1453;OSR2:-3.0079;ZEB2:-4.8734 |
| S2_k20 | O2_Y2 | FOXF1;HIC1;JDP2;OSR2;ZEB2 | -15.6866 | +14.2002 | +0.2680 | False | 1 | FOXF1:-13.9269;HIC1:-14.8646;JDP2:-13.5760;OSR2:-11.6192;ZEB2:-13.0699 |
| S2_k20 | O2_Y2 | FOXF1;HMGB3;JDP2;OSR2;ZEB2 | -5.6922 | +4.6592 | +0.2623 | False | 1 | FOXF1:-4.2080;HMGB3:-14.8646;JDP2:-3.2973;OSR2:-1.6478;ZEB2:-2.9099 |
| S2_k20 | O2_Y2 | HELZ2;IRX1;JDP2;OSR2;ZEB2 | -12.5363 | +11.3136 | +0.2852 | False | 1 | HELZ2:-10.7443;IRX1:-14.9046;JDP2:-10.1950;OSR2:-8.7689;ZEB2:-10.7543 |
| S2_k20 | O2_Y2 | FOXF1;JDP2;OSR2;TFCP2L1;ZEB2 | -9.4530 | +8.3505 | +0.2900 | False | 1 | FOXF1:-8.0175;JDP2:-7.5556;OSR2:-5.7351;TFCP2L1:-14.8646;ZEB2:-7.3211 |
| S2_k20 | O2_Y2 | FOXF1;JDP2;OSR2;ZBED2;ZEB2 | -8.6690 | +7.5429 | +0.2687 | False | 1 | FOXF1:-7.1629;JDP2:-6.2482;OSR2:-4.7227;ZBED2:-14.8646;ZEB2:-6.0104 |
| S2_k20 | O2_Y2 | HELZ2;JDP2;NEUROD1;OSR2;ZNF705B | -4.6356 | +3.6021 | +0.2763 | False | 1 | HELZ2:-2.8130;JDP2:-2.0706;NEUROD1:-2.5719;OSR2:-0.6492;ZNF705B:-15.0583 |
| S2_k20 | O2_Y2 | FOXF1;JDP2;OSR2;ZEB2;ZNF23 | -11.6484 | +10.3866 | +0.2837 | False | 1 | FOXF1:-10.1388;JDP2:-9.6946;OSR2:-7.7611;ZEB2:-9.4395;ZNF23:-14.8646 |
| S2_k20 | O2_Y2 | FOXF1;JDP2;OSR2;ZEB2;ZNF451 | -9.3877 | +8.0793 | +0.2603 | False | 1 | FOXF1:-7.6403;JDP2:-6.9365;OSR2:-5.4171;ZEB2:-6.7519;ZNF451:-14.8646 |
| S2_k20 | O2_Y2 | HELZ2;JDP2;OSR2;TARBP2;ZEB2 | -13.2275 | +12.1117 | +0.3089 | False | 1 | HELZ2:-11.7106;JDP2:-11.2026;OSR2:-9.6760;TARBP2:-14.9046;ZEB2:-11.9770 |
| S2_k20 | O2_Y2 | FOXF1;JDP2;NKX2-5;OSR2;ZEB2 | -8.8763 | +7.7800 | +0.2604 | False | 1 | FOXF1:-7.0855;JDP2:-6.3938;NKX2-5:-14.8646;OSR2:-4.8850;ZEB2:-6.2428 |
| S2_k20 | O2_Y2 | HELZ2;JDP2;OSR2;ZEB2;ZNF266 | -4.2165 | +3.3916 | +0.3149 | False | 1 | HELZ2:-2.8607;JDP2:-2.0447;OSR2:-0.7728;ZEB2:-3.0364;ZNF266:-14.9046 |
| S2_k20 | O2_Y2 | FOXF1;JDP2;MEF2B;OSR2;ZEB2 | -9.0704 | +7.9387 | +0.2567 | False | 1 | FOXF1:-7.3022;JDP2:-6.5185;MEF2B:-14.8646;OSR2:-5.0286;ZEB2:-6.3053 |
| S2_k20 | O2_Y2 | FOXF1;JDP2;OSR2;ZEB2;ZNF367 | -4.9458 | +3.9409 | +0.2694 | False | 1 | FOXF1:-3.3029;JDP2:-2.7156;OSR2:-1.0293;ZEB2:-2.3735;ZNF367:-14.8646 |
| S2_k20 | O2_Y2 | FOXF1;HOXC4;JDP2;OSR2;ZEB2 | -6.8123 | +5.7680 | +0.2711 | False | 1 | FOXF1:-5.2928;HOXC4:-14.8646;JDP2:-4.5712;OSR2:-2.8707;ZEB2:-4.2874 |
| S2_k20 | O2_Y2 | HELZ2;JDP2;OSR2;ZEB2;ZNF416 | -4.5074 | +3.8026 | +0.3092 | False | 1 | HELZ2:-2.9282;JDP2:-2.2987;OSR2:-1.0109;ZEB2:-3.2460;ZNF416:-14.9046 |
| S2_k20 | O2_Y2 | FOXF1;GRHL1;JDP2;OSR2;ZEB2 | -7.7759 | +6.5911 | +0.2599 | False | 1 | FOXF1:-6.0396;GRHL1:-14.8646;JDP2:-5.4278;OSR2:-3.7837;ZEB2:-5.0829 |
| S2_k20 | O2_Y2 | FOXF1;JDP2;OSR2;ZEB2;ZNF213 | -4.3502 | +3.2682 | +0.2657 | False | 1 | FOXF1:-2.8121;JDP2:-1.9008;OSR2:-0.3265;ZEB2:-1.6612;ZNF213:-14.8646 |
| S2_k20 | O2_Y2 | FOXF1;JDP2;OSR2;ZBTB22;ZEB2 | -3.7914 | +2.6804 | +0.2683 | False | 1 | FOXF1:-2.3559;JDP2:-1.5088;OSR2:+0.1794;ZBTB22:-14.8646;ZEB2:-0.9711 |
| S2_k20 | O2_Y2 | BRCA1;HELZ2;JDP2;NEUROD1;OSR2 | -10.6965 | +9.5781 | +0.2641 | False | 1 | BRCA1:-15.0583;HELZ2:-8.6472;JDP2:-7.9521;NEUROD1:-8.3871;OSR2:-6.7427 |
| S2_k20 | O2_Y2 | DLX2;HELZ2;JDP2;OSR2;ZEB2 | -9.7970 | +8.5805 | +0.2985 | False | 1 | DLX2:-14.9046;HELZ2:-8.2993;JDP2:-7.5698;OSR2:-6.1246;ZEB2:-8.3860 |
| S2_k20 | O2_Y2 | HELZ2;JDP2;NEUROD1;OSR2;PRDM8 | -4.0746 | +3.0198 | +0.2635 | False | 1 | HELZ2:-2.3192;JDP2:-1.2197;NEUROD1:-1.8493;OSR2:+0.0160;PRDM8:-15.0583 |
| S2_k20 | O2_Y2 | ATOH8;HELZ2;JDP2;NEUROD1;OSR2 | -8.4965 | +7.1580 | +0.2808 | False | 1 | ATOH8:-15.0583;HELZ2:-6.7058;JDP2:-6.0394;NEUROD1:-6.6434;OSR2:-4.6364 |
| S2_k20 | O2_Y2 | FOXF1;JDP2;NEUROD1;NKX2-6;OSR2 | -5.9719 | +4.7145 | +0.2562 | False | 1 | FOXF1:-4.4473;JDP2:-3.7573;NEUROD1:-3.8935;NKX2-6:-13.9766;OSR2:-1.7545 |
| S2_k20 | O2_Y2 | HELZ2;HESX1;JDP2;OSR2;ZEB2 | -4.4871 | +3.7087 | +0.2838 | False | 1 | HELZ2:-2.8586;HESX1:-14.9046;JDP2:-1.7123;OSR2:-0.7414;ZEB2:-2.5922 |
| S2_k20 | O2_Y2 | FOXF1;JDP2;OSR2;ZEB2;ZNF844 | -12.3355 | +10.9757 | +0.2983 | False | 1 | FOXF1:-10.9556;JDP2:-10.6457;OSR2:-8.6290;ZEB2:-10.3861;ZNF844:-14.8646 |
| S2_k20 | O2_Y2 | HELZ2;JDP2;MYBL2;OSR2;ZEB2 | -6.1846 | +5.2504 | +0.3070 | False | 1 | HELZ2:-4.7110;JDP2:-4.0151;MYBL2:-14.9046;OSR2:-2.6766;ZEB2:-4.9285 |
| S2_k20 | O2_Y2 | BCL6B;HELZ2;JDP2;OSR2;ZEB2 | -9.6312 | +8.5709 | +0.2937 | False | 1 | BCL6B:-14.9046;HELZ2:-7.9229;JDP2:-7.3462;OSR2:-5.9986;ZEB2:-8.0585 |
| S2_k20 | O2_Y2 | HELZ2;JDP2;OSR2;TBP;ZEB2 | -9.8515 | +8.7842 | +0.2972 | False | 1 | HELZ2:-8.1304;JDP2:-7.6839;OSR2:-6.2106;TBP:-14.9046;ZEB2:-8.3342 |
| S2_k50 | O2_Y2 | FOXF1;HELZ2;JDP2;OSR2;TUB | -5.6385 | +5.9883 | +0.2292 | False | 1 | FOXF1:-3.6740;HELZ2:-3.7940;JDP2:-4.1796;OSR2:-4.1938;TUB:-10.5779 |
| S2_k50 | O2_Y2 | DMRTC1;FOXF1;HELZ2;JDP2;OSR2 | -3.6093 | +3.4348 | +0.1976 | False | 2 | DMRTC1:-10.5779;FOXF1:-1.4224;HELZ2:-1.6332;JDP2:-1.6509;OSR2:-1.7036 |
| S2_k50 | O2_Y2 | FOXF1;GATA1;HELZ2;NR4A1;OSR2 | -1.2465 | +2.0735 | +0.2382 | False | 1 | FOXF1:+0.4697;GATA1:+0.4071;HELZ2:+0.4713;NR4A1:-10.5123;OSR2:+0.1867 |
| S2_k50 | O2_Y2 | FOXF1;GATA1;HELZ2;IRX3;OSR2 | -0.0134 | +0.5674 | +0.2262 | False | 1 | FOXF1:+1.6121;GATA1:+1.6629;HELZ2:+2.0948;IRX3:-10.5123;OSR2:+1.6453 |
| S2_k50 | O2_Y2 | FOXF1;GATA1;HELZ2;NKX1-2;OSR2 | -7.4072 | +7.0929 | +0.1930 | False | 1 | FOXF1:-5.4168;GATA1:-5.4430;HELZ2:-5.0401;NKX1-2:-10.5123;OSR2:-5.4314 |
| S2_k50 | O2_Y2 | FOXF1;GTF2IRD2B;HELZ2;JDP2;OSR2 | -7.3833 | +7.1545 | +0.2204 | False | 1 | FOXF1:-5.4729;GTF2IRD2B:-10.5779;HELZ2:-5.6176;JDP2:-5.8628;OSR2:-5.7470 |
| S2_k50 | O2_Y2 | FOXF1;GATA1;HELZ2;LDB1;OSR2 | -2.6734 | +2.9860 | +0.2166 | False | 1 | FOXF1:-0.5733;GATA1:-0.3124;HELZ2:-1.1683;LDB1:-10.5123;OSR2:-0.9151 |
| S2_k50 | O2_Y2 | FOXF1;GATA1;HELZ2;IRF9;OSR2 | -9.9654 | +9.7325 | +0.2224 | False | 1 | FOXF1:-8.2085;GATA1:-8.2647;HELZ2:-8.1247;IRF9:-10.5123;OSR2:-8.2793 |
| S2_k50 | O2_Y2 | FOXF1;HELZ2;JDP2;OSR2;TCFL5 | -7.0368 | +6.6785 | +0.2030 | False | 1 | FOXF1:-5.0563;HELZ2:-4.7293;JDP2:-5.4690;OSR2:-5.1367;TCFL5:-10.5779 |
| S2_k50 | O2_Y2 | FOXF1;HELZ2;JDP2;OSR2;ZBTB14 | -3.5172 | +3.2454 | +0.2061 | False | 1 | FOXF1:-1.6477;HELZ2:-1.4632;JDP2:-1.7738;OSR2:-1.6668;ZBTB14:-10.5779 |
| S2_k50 | O2_Y2 | FOXF1;HELZ2;JDP2;MIER2;OSR2 | -8.5717 | +8.1794 | +0.1953 | False | 1 | FOXF1:-6.4764;HELZ2:-6.4699;JDP2:-6.6805;MIER2:-10.5779;OSR2:-6.6111 |
| S2_k50 | O2_Y2 | FOXF1;HELZ2;JDP2;OSR2;ZNF569 | -4.7461 | +4.7002 | +0.2052 | False | 1 | FOXF1:-2.8254;HELZ2:-2.4720;JDP2:-3.1181;OSR2:-2.8941;ZNF569:-10.5779 |
| S2_k50 | O2_Y2 | FOXF1;HELZ2;JDP2;MYCBP;OSR2 | -5.7015 | +5.6198 | +0.2129 | False | 1 | FOXF1:-3.6980;HELZ2:-3.8015;JDP2:-4.0393;MYCBP:-10.5779;OSR2:-3.9789 |
| S2_k50 | O2_Y2 | AFF3;FOXF1;HELZ2;JDP2;OSR2 | -4.4435 | +4.4877 | +0.2276 | False | 1 | AFF3:-10.5779;FOXF1:-2.5333;HELZ2:-2.4820;JDP2:-3.1243;OSR2:-2.9499 |
| S2_k50 | O2_Y2 | FOXF1;HELZ2;JDP2;OSR2;TFAP2E | -8.0126 | +7.5754 | +0.1987 | False | 2 | FOXF1:-6.0483;HELZ2:-5.7682;JDP2:-6.2197;OSR2:-6.1049;TFAP2E:-10.5779 |
| S2_k50 | O2_Y2 | FOXF1;GATA1;HELZ2;NAT10;OSR2 | -2.1079 | +2.6460 | +0.2272 | False | 1 | FOXF1:-0.1635;GATA1:-0.3016;HELZ2:-0.4421;NAT10:-10.5123;OSR2:-0.5020 |
| S2_k50 | O2_Y2 | FOXF1;GATA1;HELZ2;NKRF;OSR2 | -4.4027 | +4.3470 | +0.1972 | False | 1 | FOXF1:-2.4713;GATA1:-2.3223;HELZ2:-2.2322;NKRF:-10.5123;OSR2:-2.4204 |
| S2_k50 | O2_Y2 | FOXF1;HELZ2;HOXB4;JDP2;OSR2 | -5.2812 | +5.1680 | +0.2162 | False | 1 | FOXF1:-3.4345;HELZ2:-3.2845;HOXB4:-10.5779;JDP2:-3.7579;OSR2:-3.5843 |
| S2_k50 | O2_Y2 | FOXF1;HELZ2;JDP2;OSR2;ZNF37A | -5.1267 | +5.2908 | +0.2073 | False | 1 | FOXF1:-3.2024;HELZ2:-2.9824;JDP2:-3.3374;OSR2:-3.3242;ZNF37A:-10.5779 |
| S2_k50 | O2_Y2 | FOXF1;HELZ2;JDP2;NFKBIB;OSR2 | -2.2546 | +2.6135 | +0.2368 | False | 1 | FOXF1:-0.4228;HELZ2:-0.6161;JDP2:-0.7950;NFKBIB:-10.5779;OSR2:-0.8043 |
| S2_k50 | O2_Y2 | FOXF1;GATA1;GTF2IRD2;HELZ2;OSR2 | -3.3888 | +4.0578 | +0.2296 | False | 1 | FOXF1:-1.5263;GATA1:-1.4777;GTF2IRD2:-10.5123;HELZ2:-1.7870;OSR2:-1.7846 |
| S2_k50 | O2_Y2 | FOXF1;HELZ2;JDP2;OSR2;ZNF692 | -5.1920 | +5.0283 | +0.2049 | False | 1 | FOXF1:-3.0975;HELZ2:-3.3196;JDP2:-3.2989;OSR2:-3.4193;ZNF692:-10.5779 |
| S2_k50 | O2_Y2 | FOXF1;HELZ2;JDP2;OSR2;ZNF517 | -5.2008 | +5.0723 | +0.2019 | False | 1 | FOXF1:-3.1942;HELZ2:-2.9917;JDP2:-3.4500;OSR2:-3.3668;ZNF517:-10.5779 |
| S2_k50 | O2_Y2 | ATOH7;FOXF1;HELZ2;JDP2;OSR2 | -0.1938 | +0.9440 | +0.2385 | False | 1 | ATOH7:-10.5779;FOXF1:+1.7445;HELZ2:+1.2229;JDP2:+1.5036;OSR2:+1.2296 |
| S2_k50 | O2_Y2 | FOXF1;HELZ2;JDP2;OSR2;PBX3 | -1.5576 | +1.8629 | +0.2185 | False | 1 | FOXF1:+0.4562;HELZ2:+0.0691;JDP2:+0.4551;OSR2:+0.1117;PBX3:-10.5779 |
| S2_k50 | O2_Y2 | FOXF1;HELZ2;JDP2;OSR2;PLAGL1 | -5.0367 | +4.8820 | +0.1912 | False | 1 | FOXF1:-2.7960;HELZ2:-2.8672;JDP2:-2.8812;OSR2:-3.1024;PLAGL1:-10.5779 |
| S2_k50 | O2_Y2 | FOXF1;GATA1;HELZ2;OSR2;SIX3 | -6.4116 | +6.8170 | +0.2323 | False | 1 | FOXF1:-4.5510;GATA1:-4.6708;HELZ2:-4.7619;OSR2:-4.9234;SIX3:-10.5123 |
| S2_k50 | O2_Y2 | FOXF1;HELZ2;JDP2;OSR2;THRB | -8.3096 | +7.8866 | +0.1994 | False | 1 | FOXF1:-6.2747;HELZ2:-6.2626;JDP2:-6.3967;OSR2:-6.4254;THRB:-10.5779 |
| S2_k50 | O2_Y2 | FOXF1;HELZ2;JDP2;OSR2;ZNF32 | -2.9902 | +3.3063 | +0.2338 | False | 1 | FOXF1:-1.2016;HELZ2:-1.0697;JDP2:-1.8264;OSR2:-1.4965;ZNF32:-10.5779 |
| S2_k50 | O2_Y2 | FOSB;FOXF1;HELZ2;JDP2;OSR2 | -11.5869 | +11.1449 | +0.2358 | False | 1 | FOSB:-10.5779;FOXF1:-10.0554;HELZ2:-9.2464;JDP2:-10.6435;OSR2:-10.0276 |
| S2_k50 | O2_Y2 | FOXF1;GLIS2;HELZ2;JDP2;OSR2 | -1.6955 | +1.7477 | +0.2168 | False | 1 | FOXF1:+0.3064;GLIS2:-10.5779;HELZ2:+0.1755;JDP2:-0.1280;OSR2:-0.0620 |
| S2_k50 | O2_Y2 | FOXF1;HELZ2;HMG20B;JDP2;OSR2 | -1.6141 | +2.1009 | +0.2245 | False | 1 | FOXF1:+0.3208;HELZ2:+0.1979;HMG20B:-10.5779;JDP2:+0.0175;OSR2:-0.0881 |
| S2_k50 | O2_Y2 | FOXF1;HELZ2;JDP2;OSR2;RUNX1T1 | -4.5040 | +4.4675 | +0.2157 | False | 1 | FOXF1:-2.5941;HELZ2:-2.5423;JDP2:-2.9181;OSR2:-2.7698;RUNX1T1:-10.5779 |
| S2_k50 | O2_Y2 | FOXF1;GATA3;HELZ2;JDP2;OSR2 | -6.1277 | +5.9086 | +0.2442 | False | 1 | FOXF1:-4.8933;GATA3:-10.5779;HELZ2:-4.0886;JDP2:-5.1099;OSR2:-4.5864 |
| S2_k50 | O2_Y2 | CERS6;FOXF1;HELZ2;JDP2;OSR2 | -4.7387 | +4.7356 | +0.2146 | False | 1 | CERS6:-10.5779;FOXF1:-2.7764;HELZ2:-2.7563;JDP2:-3.1253;OSR2:-3.0225 |
| S2_k50 | O2_Y2 | FOXF1;HELZ2;JDP2;KRBOX5;OSR2 | -4.1368 | +4.0289 | +0.2223 | False | 1 | FOXF1:-2.3183;HELZ2:-2.2289;JDP2:-2.7330;KRBOX5:-10.5779;OSR2:-2.4855 |
| S2_k50 | O2_Y2 | FOXF1;GATA1;HELZ2;OSR2;ZNF541 | -2.4726 | +2.6000 | +0.2061 | False | 1 | FOXF1:-0.5972;GATA1:-0.6732;HELZ2:-0.2042;OSR2:-0.6411;ZNF541:-10.5123 |
| S2_k50 | O2_Y2 | FOXF1;FOXH1;GATA1;HELZ2;OSR2 | -5.5054 | +5.3353 | +0.1984 | False | 1 | FOXF1:-3.5626;FOXH1:-10.5123;GATA1:-3.4844;HELZ2:-3.2234;OSR2:-3.5942 |
| S2_k50 | O2_Y2 | FOXF1;GATA1;HELZ2;MLXIP;OSR2 | -6.8312 | +6.4643 | +0.2010 | False | 1 | FOXF1:-4.9257;GATA1:-4.9015;HELZ2:-4.6187;MLXIP:-10.5123;OSR2:-4.9337 |
| S2_k50 | O2_Y2 | FOXF1;HELZ2;JDP2;OSR2;THAP12 | -6.1440 | +5.7113 | +0.2023 | False | 1 | FOXF1:-4.1246;HELZ2:-4.0542;JDP2:-4.4258;OSR2:-4.2956;THAP12:-10.5779 |
| S2_k50 | O2_Y2 | FOXF1;HELZ2;JDP2;KLF1;OSR2 | -1.7972 | +1.7133 | +0.2104 | False | 1 | FOXF1:-0.0932;HELZ2:+0.3470;JDP2:-0.0338;KLF1:-10.5779;OSR2:+0.0172 |
| S2_k50 | O2_Y2 | FOXF1;HELZ2;JDP2;OSR2;ZNF512B | -2.8227 | +3.0641 | +0.2178 | False | 1 | FOXF1:-0.9197;HELZ2:-0.9418;JDP2:-1.0449;OSR2:-1.1709;ZNF512B:-10.5779 |
| S2_k50 | O2_Y2 | FOXF1;HELZ2;JDP2;NFIX;OSR2 | -3.4183 | +3.7358 | +0.2266 | False | 1 | FOXF1:-1.4277;HELZ2:-1.8197;JDP2:-1.7488;NFIX:-10.5779;OSR2:-1.9014 |

Stage 4 separation quartiles exist; compare each reported set's ratio with `results/south6/stage4_quartiles.csv`.

## Fired keys

No report-only flag fired.

No STOP key fired.

## Limitations

- The control set was defined by the `guide_target` label (the literal string `non`, all 78 guides) rather than the file's own `control` column, which flags only 61 of those guides and omits 17. Both Stage 2B floors are reported: the primary one on all 78 guides and a report-only one on the 61 `control=True` guides.
- Hs27 is one neonatal line with no aged cells, so every delta is a counterfactual assuming the effect transfers.
- Combinations assume additivity, untestable here. No public fibroblast screen can test that assumption.
- Unmeasured genes are assumed unmoved (z = 0) and are not permuted into.
- The primary origin and target are GTEx fibroblast donors, not cells from this line. O2 is one aged donor's day-0 pseudobulk from a different study.
- Non-targeting guides are the control the experiment measured. They are not a proof that the guides are inert.
- The frozen ruler was not refit. Its score is context only and is not evidence of approach.

## Failures (verbatim)

None recorded.

## Flags

None recorded.

## Notes

- launch: Full job launched by the agent at the user's explicit instruction, after the Stage 0 smoke test returned a ruler coverage weight fraction of 0.8182 on the best mapping (symbol-only 0.6204, against the ~0.63 SURVEY_TF.md expected). This departs from the pre-registration line 'Write it, smoke-test Stage 0 only, stop. Do not run the full job yourself.' Nothing else in the pre-registration was changed and no statistic was computed before the pre-registration was written.
- Y2_singular: Y2 at k=100 rank 55 < 100. Skipped.
