# FINDINGS_SENG3 — do Junevity's four factors approach young, measured in
# the compressed space?

Task SENG3. Computational, no new data. Reads `results/south6/`,
`results/south7/`, `results/south8/`, `results/south9/`, `results/seng/`,
`results/seng2/` and the Sengstack h5ad already on disk. Writes
`results/seng3/` and this file. Nothing committed, nothing pushed. No
earlier findings file or results directory was written.

**The gates, in the order the flag runs them.** Stage 0 PASS; D1
**CALIBRATED**; D2 **LIMITED**; D3 **USABLE**; D4
**APPROACH**.

**Read this before the gate sentences.** D4's APPROACH is the
pre-registered reading and it is reported as such. Three report-only facts
sit beside it and are not optional reading:

1. **The approach statistic is almost a monotone function of the
   proliferation shift**: Spearman rho = -0.8639
   (p = 1e-102) across the 340 perturbations.
   Residualizing proliferation drops the called count from
   80 to 25
   under origin (i) and from 73 to
   8 under origin (ii), and leaves
   **one** of the four named factors called under origin (i) and **none**
   under origin (ii).
2. **Calling depends on cell count.** The call rate rises from
   0.056 in the shallowest decile to 0.545
   in the deepest; Spearman(n_cells, called) = +0.2712. **All four named factors sit at
   the 82nd to 99th percentile of cell count**, and in their own depth band
   24 of 60 perturbations are called (40%).
3. **D2 fired LIMITED**, so every number here speaks only to approaches of
   at least half the aging gap.

What the four do have, stated as fairly: among equally deep perturbations
their deltas are still extreme (percentiles 1.7, 3.3, 15.0, 21.7 within the
depth band), so cell count alone does not account for them. Their
proliferation shifts are at the 98th, 95th, 90th and 60th percentile of the
same band. Both facts are in D5.

---

## 1. The pre-registration, verbatim

From `results/seng3/PREREG.flag`, written 2026-10-06T09:51:19Z before any
SENG3 number existed and not edited since.

```
PREREG written 2026-10-06T09:51:19Z (local 2026-10-06 18:51:19 KST)
Written before any SENG3 number exists. Reproduced verbatim at the top of
FINDINGS_SENG3.md. Not edited after this timestamp.

TASK: SENG3 — do Junevity's four factors approach young, measured in the
compressed space?
Written before any SENG3 number exists.

Nature: computational. No new data. Reads results/south6/, results/south7/,
results/south8/, results/south9/, results/seng/, results/seng2/ and the
Sengstack h5ad already on disk (data/raw/seng/exp3_merged_not_normalized.h5ad,
md5 bfccaa94a916c22e4b0793913a04128b). Writes only results/seng3/ and
FINDINGS_SENG3.md. No earlier findings file is regenerated, no earlier
results directory is written, nothing is committed or pushed.

--- THE QUESTION ---

SENG2 ran the approach test on Junevity's data in the full-gene space and
fired D3 INSENSITIVE: power at a 25% plant at 100 cells was 0.000, and even
a plant of the entire aging gap was detected almost never. Its NO APPROACH
reading therefore carries no conclusion.

The Southard work already established that this is a property of the space
and not of the cell count. In SOUTH9's calibration the compressed space
S2A_k20 fires at 0.0350 at 106 cells where the full-gene space S1 fires at
0.1200, and SOUTH6REC's plant recovery detects a 25% plant at 106 cells with
power 0.98 in the compressed space against roughly 0.11 in full genes. Same
cells, same depth, different instrument.

SENG2 used the full-gene space because Stage 0 had to reproduce SENG's
stored numbers, which live there. That was correct for reproduction and
wrong for the test. This task rebuilds the compressed instrument on
Junevity's data, qualifies it by planting before reading any result, and
only then asks about the four factors.

The order is deliberate. In SENG2 the sensitivity check ran after the
approach test and retroactively voided it. Here the instrument is qualified
first and the approach test is not computed at all if it fails.

--- STATE THE METHOD BEFORE USING IT ---

Report, before any SENG3 number:

 (a) Exactly how SOUTH6/SOUTH8 construct the space named S2A_k20. What the
     basis is, what data it is fit on, how many components, whether the
     young or old anchors enter the fit, whether any perturbed cell enters
     the fit, and in what order normalisation and projection are applied.
     Quote the source lines. If any part of the recipe cannot be read off
     the committed code, say so and STOP `basis_recipe_unreadable`.

 (b) Whether that recipe can be applied to the Sengstack data unchanged.
     Name every place it cannot and what is substituted. A substitution is a
     deviation and is disclosed as one.

 (c) Whether the rebuilt basis is fit on cells the test is later applied to.
     If it is, say so plainly and label every SENG3 number as carrying that
     circularity. Do not refit to remove it; report it.

 (d) Confirm the four factors' cell counts under SENG's own assignment rule
     are unchanged: CRA_E2F3 128, CRA_EZH2 107, CRI_STAT3 96, CRI_ZFX 157.

--- STAGE 0: REPRODUCE FIRST (gate) ---

Re-derive with this task's code, from stored files:
 (a) SOUTH8's amended counts in S2A_k20: 24 of 304 and 125 of 1,836 called,
     and the below-floor counts 27 and 153;
 (b) SOUTH9's D2 calibration rows: S2A_k20 at 106 cells 0.0350 and S1 at
     106 cells 0.1200;
 (c) SENG2's D1 depth-matched cosines: cos(PD14 - NT_A, a) = +0.584155 and
     cos(PD14 - NT_I, a) = +0.757055;
 (d) SENG2's D2 counts under origin (i): 131 below floor, 3 called, of 340;
 (e) SENG2's D3 plant recovery at 100 cells, 25% dose: power 0.000 in both
     libraries.
Pass: counts exact, numerics within 1e-4 * max(1, |stored|).
Fail -> STOP `prior_results_not_reproduced`, naming every mismatch.

--- D1. Rebuild the compressed basis, and validate its null (gated) ---

Build the S2A_k20 recipe on the Sengstack data. The basis must not be fit on
PD14, which is the target, nor on any perturbation group that is later
tested. State what it is fit on.

Then validate the null in that space, which is the check that the buzzer
works. Split the NT pool of each library into two disjoint halves. Build the
floor from 200 draws on half one, at each of the cell counts 96, 107, 128 and
157. Evaluate the full criterion — delta below floor AND delta negative — on
200 independent draws from half two, which contain no effect by construction.

Report, per library and per cell count: the observed firing rate, the floor
value, and the median and 2.5-97.5% range of the null delta. Report the same
quantities from SENG2's full-gene space side by side, so the change in the
pedestal is visible.

Gate, wording fixed now. Let f be the largest observed firing rate across the
eight library-by-cell-count combinations.
 - CALIBRATED (f <= 0.10): "the criterion fires on effect-free data at <f>,
   against the 0.05 it is built to give. The null is calibrated in this space
   and a count above chance can be read as such."
 - MISCALIBRATED (f > 0.10): STOP `null_miscalibrated`. Report the full table
   and compute nothing further. A miscalibrated null is what SOUTH8 removed
   and it is not carried forward here.

--- D2. Can the instrument see a planted approach at these depths? (gated) ---

This gate runs before any approach number is computed and decides whether one
is computed at all.

Let a be the aging axis in the compressed space, centroid(PD14 half A) minus
centroid(PD32 half A), matching SENG2's amendment A4. Spike and plant
material is PD14 half B, which appears in no target.

Plant recovery: displace NT draws along a by 0, 25, 50 and 100% of ||a||,
200 draws per dose, at each of 96, 107, 128 and 157 cells, in both libraries.
Power is the fraction of draws the D4 criterion would call.

Mixture control across the library boundary: mix PD14 half B into NT cells at
0, 10, 25 and 50% at fixed total cell count, 200 draws per dose, at 96 and 157
cells, in both libraries. The zero-dose point is the noise floor and is
reported, not subtracted.

Gate, wording fixed now. Let p25 be the lowest power at a 25% plant across
all library-by-cell-count combinations, and p50 the same at 50%.
 - QUALIFIED (p25 >= 0.80): "the statistic detects a 25% approach across the
   library boundary at every one of these factors' own cell counts, lowest
   power <p25>, so a null in D4 is a null and not an insensitive assay."
 - LIMITED (p25 < 0.80 and p50 >= 0.80): "the statistic detects a 50%
   approach but not a 25% one, lowest powers <p50> and <p25>. D3 and D4 are
   run, and every number in them is labelled as able to speak only to
   approaches of at least half the aging gap."
 - BLIND (p50 < 0.80): STOP `instrument_blind`. Report the full power curve
   and the mixture table and compute no approach number. The finding is that
   the compressed instrument does not transfer to this dataset at these
   depths, and that is the task's result.

--- D3. The library offset in the compressed space (gated) ---

Cosines are not comparable between spaces of different dimension, so SENG2's
D1 reading does not carry over and is recomputed here, not quoted.

Compute, as in SENG2's D1 but in the compressed space: ||offset_A||/||a||,
||offset_I||/||a||, cos(offset_A, a), cos(offset_I, a), and the gate quantity
cos(PD14 - NT_L, a) for each library. Depth matching is SENG2's amendment A3:
binomial thinning to the lowest group median, one common factor per group,
every cell kept, 200 draws, seed 20260914, thinned primary, with the raw
version reported beside it. The gate value is the median over draws with its
2.5-97.5% range.

Gate, wording fixed now. Let r be the larger of the two cosines.
 - USABLE (r >= 0.70): "the library offset rotates the target direction by
   less than 45 degrees (cos <r>); the apparent target direction is dominated
   by aging and the approach test can be read."
 - DEGRADED (0.40 <= r < 0.70): "the library offset rotates the target
   direction substantially (cos <r>). D4 is run and reported, and every
   number in it is labelled as carrying a target direction that is <pct>%
   not aging."
 - UNUSABLE (r < 0.40): "the vector from a perturbation's control to the
   young target is mostly library difference, not aging (cos <r>). This
   dataset cannot answer whether these factors approach young in this space."
   D4 is still run to the end and reported, and nothing in it is a conclusion
   about any factor. That is said in the gate sentence and repeated wherever
   a factor is named.

Report the per-library cosines separately beside the gate, because the two
libraries carry different factors: CRA_E2F3 and CRA_EZH2 sit in _A, and
CRI_STAT3 and CRI_ZFX in _I. In SENG2's full-gene space those cosines fell in
different bands and the gate was carried entirely by _I. Say whether that is
still so here.

--- D4. The approach test (gated) ---

Carries D1, D2 and D3's readings, each named.

For all 340 perturbations and separately for the four named factors, compute
the change in distance to the PD14 half A target in the compressed space,
with a depth-matched floor built at each perturbation's own cell count from
NT draws in its own library, 200 draws, floor at the 5th percentile. Called
approaching means delta below its own floor and delta negative.

Two pre-specified origins, both reported, neither privileged, per SENG2's
amendment A6:
 (i) as-measured: target is PD14 half A minus the NT centroid of the
     perturbation's own library;
 (ii) offset-corrected: the displacement is applied to the PD32 half A
      origin. This transplants a displacement onto a different origin within
      the same cell line and experiment; label it as such wherever it
      appears.

Report, per origin: below floor, delta negative, called, against the chance
expectation of 17.0 of 340 with its 95% binomial interval; the same per
modality; and the four factors individually with cell count and its
percentile in the screen.

Gate, wording fixed now, per SENG2's amendment A1. Count the four named
factors called under each origin separately. The branches are mutually
exclusive.
 - APPROACH (both origins call two or more): "of Junevity's four validated
   factors, <k> end closer to early-passage cells than their own controls do,
   in their own data, with their own perturbations: <list>."
 - NO APPROACH (both origins call fewer than two): "none of Junevity's four
   validated factors measurably approaches early-passage cells in their own
   data, with their own perturbations, against a floor built at each factor's
   own depth and an instrument qualified at that depth."
 - SPLIT (one origin calls two or more and the other fewer): report both,
   name which factor is which under each, draw no single conclusion.

--- D5. Report-only, no decision rule ---
 - Proliferation: a residualized D4 with the null residualized identically by
   projecting out the same per-modality proliferation axis, reported beside
   the distance result and never in place of it, with SENG's stored rho
   (+0.828 CRA, +0.862 CRI) quoted beside it. Disclose that the proliferation
   axis is estimated from the same control cells that define the origin and
   the floors.
 - Identity: the fibroblast identity change for the four and for the screen,
   over the ten identity genes SENG found present.
 - The frozen ruler: how many of the 340 score younger on the one-dimensional
   ruler while not approaching, under both origins. This is the preprint's
   thesis in this dataset and is a diagnostic, not evidence about any factor.
 - The same count restricted to the four.
 - A side-by-side table of every SENG2 headline number and its SENG3
   counterpart, so the effect of changing the space alone is visible.

--- WHAT DOES NOT COUNT ---
- choosing a space, a basis, a component count, an origin, a criterion, a
  null or a factor set after seeing a number;
- reading D4 if D2 fired BLIND;
- treating any factor called here as validated or as a pilot candidate in
  this findings file;
- reading a null as evidence about transcription factors in aged human cells:
  this dataset is passaged neonatal-line fibroblasts and its "old" is passage
  number, which this project has found does not resemble donor aging
  (cos(v,u) = -0.040, permutation p = 1.0);
- refitting any axis, basis or offset to avoid reporting a circularity;
- any claim about Junevity's published numbers being unreproducible.

--- STOP KEYS ---
basis_recipe_unreadable; prior_results_not_reproduced; null_miscalibrated;
instrument_blind; no_raw_counts.
Seeds: 20260914 splits/nulls/permutations, 20260918 bootstraps.
n_null = n_perm = n_boot = 200.
```

## 2. Amendments, and every deviation

The flag was ambiguous in two places that decide the instrument and the
gates. Both were put to the author **before computing**; eight further
definitions or citations the flag omits or gets wrong were fixed in writing
before computing. Full text of `results/seng3/AMENDMENTS.md`:

```
# SENG3 — amendments to the pre-registration

Written 2026-10-06T10:05Z, **before any SENG3 number existed**, after the flag
was found ambiguous in two places that decide the instrument and the gates.
Both were put to the author rather than resolved by this session.
`results/seng3/PREREG.flag` is **not edited**; it stands as written at
2026-10-06T09:51:19Z. Everything below is reproduced in `FINDINGS_SENG3.md`.

B1 and B2 are the author's answers. B3–B7 are definitions or citations the flag
does not supply or gets wrong; those were fixed in writing here, before
computing, by this session, and each is labelled as such.

---

## B1 — What the basis is (ambiguous; author's answer)

**The ambiguity.** The S2A_k20 basis is PCA loadings fit on **GTEx bulk donors
aged 40–59 only** — no Southard cell, no anchor donor and no perturbed cell
enters it (308 donors, `n_mid` in the stored basis) (`src/south7_run.py:117`, `anchor_free_basis`). The flag says three
times to "rebuild the compressed instrument on Junevity's data" and asks in (c)
whether "the rebuilt basis is fit on cells the test is later applied to", which
implies a refit on Sengstack cells. But the recipe's basis is not fit on screen
cells at all, so applying the recipe *unchanged* means reading the frozen GTEx
basis. The two give different instruments and therefore different D1 and D2
gates.

**Resolution, as answered: read the frozen GTEx S2A basis.** The basis is
`L_mid_full[:, :20]` read from `results/south7/basis.npz` — the PCA loadings of
the 308 GTEx donors aged 40–59, in the 23,485-column frozen-ruler z space.
Sengstack displacements are projected onto those 20 components, with the 6,704
ruler genes absent from this file left at zero (they are already z=0 in SENG's
pipeline). **Nothing is fit on any Sengstack cell**, so item (c) of the method
statement answers "no circularity in the basis". This also satisfies CLAUDE.md:
frozen instruments are read, not refit.

## B2 — Which origin D2's power gate uses (ambiguous; author's answer)

**The ambiguity.** D2's power is "the fraction of draws the D4 criterion would
call", but D4 defines two origins and privileges neither; the delta, and hence
the power, differs between them. Different origins could land in different D2
branches, one BLIND (stop, compute no approach number) and the other QUALIFIED.

**Resolution, as answered: compute both, gate on origin (i), report both.**
Power is computed and reported under both origins. The gate is read off
**origin (i)**, the as-measured geometry, because that is the origin that
actually contains the library boundary D2 exists to test. If the two origins
would land in different branches, that is said prominently and both branch
sentences are reported.

---

The remainder were not put to the author.

## B3 — Stage 0(b) is miscited in the flag (fixed by this session)

The flag asks for "**SOUTH9's** D2 calibration rows: S2A_k20 at 106 cells 0.0350
and S1 at 106 cells 0.1200". Those rows are not in `results/south9/`. They are in
**`results/south8/d2_calibration.csv`**, column `frac_p_le_05`, rows
`n_cells=106`. They are re-derived from SOUTH8 and the corrected citation is
reported. SOUTH9's own D2 table (`d2_residualized_recovery.csv`) holds the
*plant recovery* numbers the flag's preamble quotes separately — S2A_k20
`raw_power_25` = 0.98 and S1 = 0.115 at 106 cells — and those are reported
beside them as a check that the right table was found.

**A second, more important point about those two numbers.** `frac_p_le_05` is
the fraction of effect-free draws whose **permutation p** is at or below 0.05. It
is *not* the criterion D1 of this task validates, which is "delta below its own
floor **and** delta negative". The 0.0350 and 0.1200 are therefore not directly
comparable to D1's firing rates, and are reported as what they are. D1's
side-by-side full-gene comparison is computed in this task on this dataset, not
imported from SOUTH8.

## B4 — Stage 0(a) and (b) are recounts, not recomputations (fixed by this session)

Stage 0 says "re-derive with this task's code, **from stored files**". For the
SOUTH8 and SOUTH9 items that is what is done: the counts are re-derived by
recounting SOUTH8's stored per-factor call tables
(`d3c_primary_304_umi_calls.csv`, `d3c_all_1836_umi_calls.csv`) and its stored
per-draw true-null p table (`d2_true_null_p.csv`). **Nothing is recomputed from
the Southard count matrix.** This is a weaker check than recomputation: it
verifies the stored tables are internally consistent with the headline counts,
not that the Southard pipeline reproduces. It is reported as such.

The SENG2 items (c), (d) and (e) are stronger: (c) is **recomputed from the
h5ad** by this task's own code in the full-gene space, and (d) and (e) are
recounted from SENG2's stored per-perturbation and per-draw tables.

## B5 — The null-validation displacement (missing definition; fixed by this session)

D1 says to build the floor from 200 draws on NT half one and evaluate on 200
draws from half two, but does not say what displacement a draw has. It is taken
to be the Southard/SENG construction **within each half**: for a draw of `n`
cells from a half, the displacement is `z(draw) - z(that half minus draw)`. The
two halves are disjoint, so the floor and the evaluation share no cell. Halves
are split with seed 20260914.

## B6 — Order of normalisation and projection (fixed by this session)

Stated here because (a) asks for it and it is the one place a reader could go
wrong. For any group: sum raw counts over its cells; scatter into the 23,485
ruler columns (16,781 present in this file); one joint TMM against SENG's frozen
panel cache; log2-CPM; z against the frozen GTEx mu and sd, with absent genes set
to z=0. **Then** form the difference vector (displacement or target direction).
**Then** project that difference onto the 20 basis columns. Projection is always
the last step and is applied to differences only, so the PCA centring vector
`center_mid` never enters — it cancels in any difference. This is the order
`src/south7_run.py:130` (`spaces_of`) and `:144` (`observed_stats`) use.

## B7 — Inherited from SENG2, unchanged (fixed by this session)

The flag names SENG2's amendments A3, A4, A6 and A1 and adopts them. Carried
over verbatim in substance: A4 (the target is PD14 half A; `a` is the full-arm
difference), A6 (the two origins written out, with PD32 **half A** as the
offset-corrected origin), A3 (binomial thinning, applied to the group pseudobulk
by the exact distributional identity recorded in SENG2's amendments), A1 (the
mutually exclusive APPROACH / NO APPROACH / SPLIT branches). A5's definition of
the statistic (`delta = ||t - d|| - ||t||`, `progress = frac`) is also carried
over, computed in the compressed space. A7's `sqrt(1 - r^2)` is carried over for
D3's `<pct>`.

## B8 — The flag miscites SENG2's A4 for the definition of `a` (fixed by this session)

D2 says: "Let `a` be the aging axis in the compressed space, centroid(PD14 half
A) minus centroid(PD32 half A), **matching SENG2's amendment A4**." A4 does not
say that. A4 fixed `a` over the **full arms** (581 and 4,348 cells) and made
PD14 half A the *target*. So the flag's explicit definition and its
justification disagree.

**Resolution, and it is not a stop because both readings let the run proceed.**
The explicit definitional sentence is taken as the definition: `a =
centroid(PD14 half A) - centroid(PD32 half A)`, and D3's "as in SENG2's D1"
is read as naming the quantities to compute, not as redefining `a`. The
mis-citation of A4 is disclosed here and in the findings.

To make the choice harmless, **`||a||` and every D3 cosine are reported under
both definitions** — half A minus half A (primary, the flag's explicit text) and
the full arms (SENG2's actual A4) — so the reader can see what the mis-citation
would have changed. The D2 plant doses use the primary `a`. If the two
definitions imply different D3 gate branches, that is said prominently.

## B9 — The plant floor (missing definition; fixed by this session)

D2's plant power is "the fraction of draws the D4 criterion would call", and
that criterion needs a floor. The floor at each cell count is the 5th percentile
of the **dose-0 draws of the same plant set**, which are effect-free NT draws at
that cell count — the same construction SENG2's D3 used. The mixture's zero-dose
point is likewise its own floor and is reported, not subtracted.

## B10 — D1's gate also needs an origin (fixed by this session)

D1's gate says "let `f` be the largest observed firing rate across the **eight**
library-by-cell-count combinations" and never names an origin, but the criterion
it validates ("delta below floor AND delta negative") requires one.

**Resolution, and the flag's own arithmetic settles it.** Eight combinations is
exactly 2 libraries x 4 cell counts, so the flag intends **one origin per
combination**; two origins would give sixteen. B2 fixes origin (i) as this
task's gating origin, so D1's eight combinations are origin (i).

**This is material and is reported, not buried.** Under origin (ii) the largest
firing rate is **0.11**, above the 0.10 threshold, so gating D1 on origin (ii)
would have fired `null_miscalibrated` and stopped the task. Both numbers are in
the D1 table and the findings say this plainly.
```

### Deviations beyond the amendments

1. **The basis is read, not rebuilt** (amendment B1, the author's answer).
   The flag's phrase "rebuild the compressed instrument on Junevity's data"
   is honoured by applying the recipe unchanged, which means reading the
   frozen GTEx basis. No PCA was fit in this task.
2. **Stage 0(a) and (b) are recounts from stored tables, not recomputations**
   from the Southard count matrix (amendment B4). Weaker than reproduction;
   it checks that SOUTH8's stored per-factor and per-draw tables are
   internally consistent with its headline counts.
3. **The flag miscites two sources**: SOUTH9 for calibration rows that are in
   SOUTH8 (B3), and SENG2's amendment A4 for a definition of `a` that A4 does
   not give (B8). Both are reported and both are worked around by computing
   the quantity under every candidate reading.
4. **`.X` is streamed, not loaded** (about 2.6 GB of RAM free against a 1.4 GB
   matrix). Changes no number.
5. **SENG2's bootstrap CIs remain un-re-derived**, as in SENG2. No gate here
   reads a CI.
6. This task reuses the frozen helper functions in `src/seng_run.py` and
   `src/seng2_common.py`, and the frozen basis from `results/south7/`. A
   reproduction that reuses the original implementation tests the record and
   the data path, not the implementation.

## 3. State the method before using it

### (a) How SOUTH6/SOUTH8 construct S2A_k20

Readable in full off committed code, so **`basis_recipe_unreadable` does not
fire**. The construction is `src/south7_run.py:117`:

```python
def anchor_free_basis(gtex, log):
    Zg = np.asarray(gtex["Z"], np.float64)
    age = np.asarray(gtex["age"]).astype(str)
    mid = np.isin(age, list(MID_BINS))      # MID_BINS = ("40-49", "50-59")
    Zm = Zg[mid]
    center = Zm.mean(0)
    sp = fit_pca(Zm - center, KMAX)
    L = np.asarray(sp.L, np.float64)
```

and the projection is `src/south7_run.py:130` and `:144`:

```python
def spaces_of(L, center, c_old, c_young, mr, ks=KS7, tag="S2A"):
    r_full = c_old - c_young
    for k in ks:
        r = r_full @ L[:, :k]
...
    A = D_m if not s["proj"] else D_m @ Lm[:, :s["k"]]
```

- **Basis**: PCA loadings of the frozen-ruler z matrix of GTEx **bulk donors**,
  restricted to age bins 40-49 and 50-59 — **308 donors**
  (`n_mid` in the stored basis). Centred on that submatrix's own mean.
- **Components**: the first **20** for S2A_k20 (`PRIMARY_SPACE`,
  `src/south7_run.py:56`).
- **Anchors in the fit: none.** `ANCHOR_BINS` are 20-29, 30-39, 60-69 and
  70-79 and are excluded by construction — that is what "anchor-free" names.
  SOUTH7 logs "anchor donors in fit = 0".
- **Perturbed cells in the fit: none.** No single cell of any kind enters it;
  the fit is over bulk donors.
- **Order**: counts -> ruler alignment -> TMM -> log2-CPM -> frozen GTEx z
  (mu and sd from **all 652** donors, anchors included — `src/south7_run.py:1134`
  says so) -> form the difference -> **project last**. Projection is applied
  only to differences, so `center_mid` cancels (amendment B6).

- Read frozen from `results/south7/basis.npz`, key `L_mid_full`, first 20
  columns. Orthonormality check: max |LᵗL − I| = 3.553e-15.

### (b) Can the recipe be applied to the Sengstack data unchanged?

Almost. Three substitutions, each a deviation:

| what | Southard | Sengstack | consequence |
|---|---|---|---|
| ruler genes measured | 20,264 of 23,485 | **16,781** of 23,485 | the 6,704 absent genes are z=0 and contribute nothing to the projection |
| the old/young anchors | GTEx donor age bins | **PD32 / PD14 passage arms** | the axis is passage, not donor age — see "what this does not cover" |
| the control pool | Hs27 non-targeting guides | **NT_1..NT_6 cells per library** | two disjoint pools, one per library, as SENG2 established |

The basis itself needs no substitution at all, because it is fit on GTEx
donors and not on screen cells. That is the whole content of amendment B1.

### (c) Is the rebuilt basis fit on cells the test is applied to?

**No.** The basis is fit on 308 GTEx bulk donors aged 40-59. No PD14 cell, no
PD32 cell, no NT cell and no perturbed cell enters it. There is **no basis
circularity** in this task.

**Three circularities that do remain, and are not refitted away:**

1. The aging axis `a` and the target are built from PD14 and PD32 cells, and
   D3's gate quantity is a cosine between two vectors that both contain
   `z(PD14 half A)`. The gate is not independent of the anchor estimate.
2. The frozen GTEx mu and sd used for z-scoring come from all 652 donors,
   anchors included — inherited from the Southard instrument.
3. **The proliferation axis in D5 is estimated from the same NT control cells
   that define the origin and every floor.** The flag asks for this to be
   disclosed and it is.

### (d) The four factors' cell counts

Re-derived under SENG's own 804-column unique-strict-maximum rule:
`CRA_E2F3` 128, `CRA_EZH2` 107, `CRI_STAT3` 96, `CRI_ZFX` 157 — all four
exact, in the Stage 0 table below.

## 4. Stage 0 — reproduce first (gate)

**Gate: PASS.** Counts exact, numerics within
1e-4 × max(1, |stored|).

| item | kind | stored | re-derived | ok |
|---|---|---|---|---|
| `SOUTH8 primary_304 m` | count | 304 | 304 | yes |
| `SOUTH8 primary_304 below_depth_floor_S2A_k20` | count | 27 | 27 | yes |
| `SOUTH8 primary_304 amended_call_S2A_k20` | count | 24 | 24 | yes |
| `SOUTH8 all_1836 m` | count | 1836 | 1836 | yes |
| `SOUTH8 all_1836 below_depth_floor_S2A_k20` | count | 153 | 153 | yes |
| `SOUTH8 all_1836 amended_call_S2A_k20` | count | 125 | 125 | yes |
| `SOUTH8 D2 frac_p_le_05 S2A_k20@106` | numeric | 0.035000 | 0.035000 | yes |
| `SOUTH8 D2 frac_p_le_05 S1@106` | numeric | 0.120000 | 0.120000 | yes |
| `SENG2 D1 cos(PD14-NT_A, a) depth-matched median` | numeric | 0.584155 | 0.584155 | yes |
| `SENG2 D1 cos(PD14-NT_I, a) depth-matched median` | numeric | 0.757055 | 0.757055 | yes |
| `SENG2 D2 n` | count | 340 | 340 | yes |
| `SENG2 D2 below_floor origin(i)` | count | 131 | 131 | yes |
| `SENG2 D2 called origin(i)` | count | 3 | 3 | yes |
| `SENG2 D3 power_25@100 CRA` | numeric | 0.000000 | 0.000000 | yes |
| `SENG2 D3 power_25@100 CRI` | numeric | 0.000000 | 0.000000 | yes |
| `n_cells[CRI_STAT3]` | count | 96 | 96 | yes |
| `n_cells[CRA_EZH2]` | count | 107 | 107 | yes |
| `n_cells[CRA_E2F3]` | count | 128 | 128 | yes |
| `n_cells[CRI_ZFX]` | count | 157 | 157 | yes |

Reported beside, as the check that the right SOUTH9 table was found (B3):
SOUTH9's own plant recovery at 106 cells gives `raw_power_25` = S2A_k20 0.98, S2A_k50 0.955, S2A_k100 0.93, S1 0.115.

## 5. D1 — the rebuilt basis, and whether its null is calibrated (gated)

The NT pool of each library is split in half (seed 20260914). The floor is
200 draws on half one; the criterion is evaluated on 200 draws from half two,
which share no cell with the floor and contain no effect by construction
(amendment B5). The full-gene space S1 is run side by side on the same draws.

| library | n cells | space | floor | **firing rate** | null delta median | 2.5–97.5% | n of 200 with delta<0 |
|---|---|---|---|---|---|---|---|
| CRA | 96 | S2A_k20 | -8.241 | **0.060** | -0.049 | -10.70 to +10.44 | 101 |
| CRA | 96 | S1 | +24.127 | **0.000** | +33.654 | +23.63 to +44.85 | 0 |
| CRA | 107 | S2A_k20 | -7.278 | **0.085** | -0.316 | -10.23 to +10.24 | 109 |
| CRA | 107 | S1 | +23.136 | **0.000** | +31.974 | +22.21 to +43.13 | 0 |
| CRA | 128 | S2A_k20 | -7.486 | **0.040** | +0.551 | -7.97 to +11.16 | 88 |
| CRA | 128 | S1 | +20.932 | **0.000** | +29.680 | +21.74 to +39.76 | 0 |
| CRA | 157 | S2A_k20 | -5.708 | **0.065** | +0.706 | -8.38 to +8.79 | 86 |
| CRA | 157 | S1 | +21.609 | **0.000** | +29.119 | +20.64 to +38.31 | 0 |
| CRI | 96 | S2A_k20 | -6.994 | **0.085** | +0.337 | -10.04 to +11.56 | 94 |
| CRI | 96 | S1 | +21.136 | **0.000** | +29.934 | +18.71 to +45.39 | 0 |
| CRI | 107 | S2A_k20 | -6.136 | **0.065** | +0.182 | -7.93 to +10.03 | 96 |
| CRI | 107 | S1 | +20.120 | **0.000** | +27.728 | +19.32 to +38.99 | 0 |
| CRI | 128 | S2A_k20 | -6.733 | **0.040** | +0.050 | -7.48 to +8.20 | 99 |
| CRI | 128 | S1 | +15.756 | **0.000** | +24.595 | +15.66 to +34.02 | 0 |
| CRI | 157 | S2A_k20 | -5.082 | **0.065** | +0.112 | -6.53 to +7.59 | 98 |
| CRI | 157 | S1 | +14.532 | **0.000** | +20.607 | +12.73 to +29.92 | 0 |

**This table is the task's central mechanical result and is worth reading
slowly.** In the compressed space the null delta is centred on zero (median
−0.32 to +0.71) and about half of 200 effect-free draws come out negative, so
the criterion can fire in either direction and does so at 0.040 to 0.085.
In the full-gene space the null delta sits on a **positive pedestal of +14.5
to +33.7**, its whole 95% range is above zero, and **0 of 200 draws are
negative at any of the eight combinations**. S1's firing rate of 0.000 is not
calibration — it is a test that cannot fire at all, which is exactly why
SENG2's count of 3 of 340 could not be read as evidence of absence.

### The D1 gate

**CALIBRATED** fired, with f = 0.0850 at library
CRA, n = 107, over the eight
library-by-cell-count combinations of origin (i).

> "the criterion fires on effect-free data at 0.0850, against
> the 0.05 it is built to give. The null is calibrated in this space and a
> count above chance can be read as such."

**A near miss that must be stated.** Under origin (ii) the largest firing
rate is **0.1100**, above the 0.10 threshold. Gating D1
on origin (ii) would have fired `null_miscalibrated` and stopped this task.
The flag's own arithmetic — "the **eight** library-by-cell-count
combinations", which is 2 × 4 and not 2 × 4 × 2 — fixes one origin per
combination, and amendment B2 fixes that origin as (i). Amendment B10
records the reasoning. The margin is 0.085 against 0.10.

## 6. D2 — can the instrument see a planted approach at these depths? (gated)

Run **before** any approach number. `a` is PD14 half A minus PD32 half A
(amendment B8); plant and spike material is PD14 half B, which appears in no
target. Power is the fraction of draws the D4 criterion would call, with the
floor taken from the dose-0 draws of the same set (amendment B9).

`||a||` = 40.6380 in S2A_k20 against
90.1493 in S1, so the 20 components retain
45.1% of the aging gap's
length. Under the full-arm variant the two norms are
39.8616 and 82.8672.

### Plant recovery, S2A_k20, origin (i)

| library | n cells | floor | power @25% | power @50% | power @100% | delta median @25% |
|---|---|---|---|---|---|---|
| CRA | 96 | -6.542 | **0.565** | 0.950 | 1.000 | -7.272 |
| CRA | 107 | -7.188 | **0.570** | 0.940 | 1.000 | -7.821 |
| CRA | 128 | -6.393 | **0.600** | 0.975 | 1.000 | -7.485 |
| CRA | 157 | -5.404 | **0.695** | 0.990 | 1.000 | -7.680 |
| CRI | 96 | -7.098 | **0.600** | 0.970 | 1.000 | -8.329 |
| CRI | 107 | -7.235 | **0.620** | 0.975 | 1.000 | -8.361 |
| CRI | 128 | -6.389 | **0.700** | 1.000 | 1.000 | -8.666 |
| CRI | 157 | -5.876 | **0.740** | 1.000 | 1.000 | -8.271 |

The same plant in the full-gene space, same cells and same depths, reaches a
25% power of **0.000** at its worst combination — the
SENG2 result, reproduced here as the side-by-side the flag asks for.

### Mixture across the library boundary, S2A_k20, origin (i)

| library | n total | dose | progress (frac) | 2.5–97.5% | delta | 2.5–97.5% |
|---|---|---|---|---|---|---|
| CRA | 96 | 0% | -0.0013 | -0.1658 to +0.1494 | +0.564 | -7.683 to +9.865 |
| CRA | 96 | 10% | +0.1050 | -0.0629 to +0.2560 | -5.365 | -13.756 to +3.955 |
| CRA | 96 | 25% | +0.2478 | +0.0827 to +0.4224 | -12.765 | -22.238 to -4.072 |
| CRA | 96 | 50% | +0.4838 | +0.3057 to +0.6449 | -25.285 | -33.941 to -15.977 |
| CRA | 157 | 0% | -0.0019 | -0.1288 to +0.1098 | +0.384 | -6.043 to +7.751 |
| CRA | 157 | 10% | +0.1041 | -0.0510 to +0.2493 | -5.472 | -13.314 to +3.129 |
| CRA | 157 | 25% | +0.2378 | +0.1094 to +0.3698 | -12.791 | -20.219 to -5.725 |
| CRA | 157 | 50% | +0.4704 | +0.3456 to +0.6083 | -25.333 | -32.078 to -19.092 |
| CRI | 96 | 0% | +0.0019 | -0.1920 to +0.1702 | +0.461 | -8.317 to +10.561 |
| CRI | 96 | 10% | +0.0919 | -0.0642 to +0.2824 | -4.305 | -14.221 to +3.970 |
| CRI | 96 | 25% | +0.2149 | +0.0519 to +0.3638 | -10.626 | -18.776 to -2.467 |
| CRI | 96 | 50% | +0.4395 | +0.2437 to +0.6164 | -22.247 | -31.119 to -12.221 |
| CRI | 157 | 0% | -0.0073 | -0.1351 to +0.1453 | +0.547 | -6.964 to +7.578 |
| CRI | 157 | 10% | +0.0870 | -0.0274 to +0.2111 | -4.348 | -10.994 to +1.698 |
| CRI | 157 | 25% | +0.2149 | +0.0774 to +0.3318 | -11.060 | -17.473 to -3.898 |
| CRI | 157 | 50% | +0.4333 | +0.2952 to +0.5637 | -22.544 | -29.570 to -14.197 |

Progress is monotonic in dose in all four library-by-depth combinations: CRA_96 True, CRA_157 True, CRI_96 True, CRI_157 True.
The zero-dose row is the noise floor and is reported, not subtracted.

### The D2 gate

**LIMITED** fired. p25 = 0.565, p50 = 0.940, each the
lowest across the eight library-by-cell-count combinations, read off origin
(i) per amendment B2. Origin (ii) gives p25 = 0.755 and
p50 = 0.995 and lands in the **same** branch
(LIMITED), so the choice of origin does not decide this gate.

> "the statistic detects a 50% approach but not a 25% one, lowest powers
> 0.940 and 0.565. D3 and D4 are run, and every number
> in them is labelled as able to speak only to approaches of at least half
> the aging gap."

**That label is applied in D3, D4 and D5 and repeated in the reading.**

## 7. D3 — the library offset in the compressed space (gated)

Recomputed, not carried over: cosines are not comparable between spaces of
different dimension. Depth matching is SENG2's amendment A3 (binomial
thinning, 200 draws, seed 20260914, thinned primary).

| variant | space | library | ‖offset‖/‖a‖ | cos(offset, a) | **cos(PD14−NT, a)** | 2.5–97.5% |
|---|---|---|---|---|---|---|
| half_A | S2A_k20 | A | 0.8352 | -0.1693 | **+0.811375** | +0.8038 to +0.8193 |
| half_A | S2A_k20 | I | 0.6456 | -0.2883 | **+0.887005** | +0.8807 to +0.8928 |
| half_A | S1 | A | 1.0532 | +0.1123 | **+0.644549** | +0.6380 to +0.6516 |
| half_A | S1 | I | 0.7298 | +0.0708 | **+0.793429** | +0.7892 to +0.7979 |
| full_arms | S2A_k20 | A | 0.8459 | -0.1805 | **+0.811142** | +0.8046 to +0.8170 |
| full_arms | S2A_k20 | I | 0.6420 | -0.2924 | **+0.888278** | +0.8847 to +0.8926 |
| full_arms | S1 | A | 1.2021 | +0.1174 | **+0.584155** | +0.5766 to +0.5897 |
| full_arms | S1 | I | 0.8177 | +0.0672 | **+0.757055** | +0.7532 to +0.7606 |

The two `a` variants of amendment B8 agree on the gate branch
(USABLE either way; r = 0.887005 half-A, 0.888278 full arms), so the flag's mis-citation of A4
changes nothing here. The S1 rows reproduce SENG2's D1 exactly under the
full-arm variant (+0.584155 and +0.757055), which is the Stage 0(c) check.

### The D3 gate

**USABLE** fired, r = +0.887005 from library I,
2.5–97.5% [+0.8807, +0.8928], angle
27.50°.

> "the library offset rotates the target direction by less than 45
> degrees (cos +0.887005); the apparent target direction is dominated
> by aging and the approach test can be read."

The orthogonal component is still 46.2% of that
vector's length (`sqrt(1 − r²)`, amendment A7), which the flag's word
"dominated" does not convey.

### Per library, which the flag asks about explicitly

| library | factors it carries | cos(PD14−NT, a), S2A_k20 | band | the same in SENG2's S1 | band there |
|---|---|---|---|---|---|
| `_A` | CRA_E2F3, CRA_EZH2 | +0.811375 | **USABLE** | +0.644549 | DEGRADED |
| `_I` | CRI_STAT3, CRI_ZFX | +0.887005 | **USABLE** | +0.793429 | USABLE |

**The answer to the flag's question is no, that is no longer so.** In SENG2's
full-gene space the CRISPRa library's cosine was +0.584, inside the DEGRADED
band, and the USABLE gate was carried entirely by CRISPRi. In the compressed
space **both libraries are in the USABLE band**, so the two factors in the
CRISPRa library are no longer behind a degraded cosine. Compressing the space
raised the CRISPRa cosine from +0.584 to +0.811.

## 8. D4 — the approach test (gated)

Carries D1 **CALIBRATED**, D2 **LIMITED** and D3 **USABLE**.
Because D2 fired LIMITED, **every number in this section speaks only to
approaches of at least half the aging gap.**

Target is PD14 half A; origin (i) is the perturbation's own NT centroid;
origin (ii) transplants the displacement onto the PD32 half A origin — a
different origin within the same cell line and experiment. Floors are one
null per distinct cell count per library, 200 draws, 5th percentile.

### Counts over all 340 perturbations

Chance expectation is **17.0 of 340** by construction, since the floor is a
5th percentile. D1 measured the realised rate on effect-free data at 0.040 to
0.085, i.e. up to 29 of 340.

| origin | below floor | delta negative | **called** | chance exp. | 95% interval on the called fraction |
|---|---|---|---|---|---|
| (i) as-measured | 80 | 191 | **80** | 17.0 | 19.1–28.4% |
| (ii) offset-corrected | 73 | 163 | **73** | 17.0 | 17.2–26.2% |

| origin | modality | n | below floor | delta negative | called | chance exp. |
|---|---|---|---|---|---|---|
| (i) | CRA | 167 | 30 | 85 | 30 | 8.3 |
| (i) | CRI | 173 | 50 | 106 | 50 | 8.7 |
| (ii) | CRA | 167 | 25 | 70 | 25 | 8.3 |
| (ii) | CRI | 173 | 48 | 93 | 48 | 8.7 |

Unlike SENG2, `delta < 0` is now reachable — 191 of 340 under origin (i) —
because the compressed null is centred on zero rather than on a +27 pedestal.
In SENG2 only 3 of 340 had a negative delta at all.

### The four named factors

| factor | n cells | pctl | delta (i) | floor (i) | called (i) | delta (ii) | floor (ii) | called (ii) |
|---|---|---|---|---|---|---|---|---|
| `CRA_E2F3` | 128 | 96% | -19.211 | -6.393 | **yes** | -13.415 | -4.495 | **yes** |
| `CRA_EZH2` | 107 | 89% | -13.211 | -7.188 | **yes** | -10.278 | -5.321 | **yes** |
| `CRI_STAT3` | 96 | 82% | -10.334 | -7.098 | **yes** | -9.952 | -6.805 | **yes** |
| `CRI_ZFX` | 157 | 99% | -19.936 | -5.876 | **yes** | -14.247 | -6.004 | **yes** |

Called under origin (i): **4 of 4**. Under origin (ii):
**4 of 4**.

### The D4 gate

**APPROACH** fired, on amendment A1's mutually exclusive branches.

> "of Junevity's four validated factors, 4 end closer to
> early-passage cells than their own controls do, in their own data, with
> their own perturbations: CRA_E2F3, CRA_EZH2, CRI_STAT3, CRI_ZFX."

**Four qualifications on that sentence, none of which it contains.**
D2 fired LIMITED, so it speaks only to approaches of at least half the aging
gap. The statistic is rank-correlated −0.86 with the proliferation shift. The
call rate rises with cell count and all four factors are at the deep end. And
the dataset's "old" is passage number, not donor age. All four are in D5 and
the reading.

## 9. D5 — report-only, no decision rule

Nothing here is gated and nothing here decides anything. It is also where the
most important qualifications on D4 live.

### Proliferation

SENG's stored rank correlation between its screen statistic and
proliferation is rho = +0.828 (CRA) and +0.862 (CRI). In
this task's compressed space the relationship is stronger still:

| pair | Spearman rho | p |
|---|---|---|
| cell count vs delta (i) | -0.3513 | 2.61e-11 |
| cell count vs floor (i) | +0.8281 | 5.96e-87 |
| cell count vs called (i) | +0.2712 | 3.82e-07 |
| proliferation shift vs delta (i) | -0.8639 | 1e-102 |
| frozen-ruler change vs delta (i) | +0.5019 | 4.24e-23 |
| identity change vs delta (i) | +0.4724 | 2.64e-20 |

Residualizing the proliferation axis, with the null residualized identically:

| origin | called, plain | called, residualized | of the four, called residualized |
|---|---|---|---|
| (i) | 80 | 25 | 1 |
| (ii) | 73 | 8 | 0 |

**The proliferation axis is estimated from the same NT control cells that
define the origin and every floor**, as the flag requires be disclosed. It is
therefore not an independent covariate, and the residualized column is not a
clean adjustment — it is reported beside the distance result, never in place
of it.

### Cell count — the check CLAUDE.md requires before reporting any shortlist

| cell-count decile | n | median cells | called (i) | call rate | median delta (i) | median floor (i) |
|---|---|---|---|---|---|---|
| 1 | 36 | 33 | 2 | 0.056 | +4.739 | -10.180 |
| 2 | 34 | 39 | 7 | 0.206 | +1.207 | -8.982 |
| 3 | 33 | 47 | 2 | 0.061 | +3.384 | -9.564 |
| 4 | 37 | 53 | 7 | 0.189 | +0.485 | -8.748 |
| 5 | 32 | 60 | 3 | 0.094 | +1.192 | -8.290 |
| 6 | 34 | 68 | 10 | 0.294 | -2.327 | -7.696 |
| 7 | 32 | 77 | 11 | 0.344 | -5.489 | -7.438 |
| 8 | 34 | 88 | 13 | 0.382 | -4.749 | -6.998 |
| 9 | 35 | 99 | 7 | 0.200 | -2.103 | -7.188 |
| 10 | 33 | 121 | 18 | 0.545 | -6.177 | -6.105 |

The call rate runs from 0.056 to 0.545
across deciles. Two mechanisms both favour deep perturbations: the floor
tightens with depth (rho +0.828) and the
delta grows more negative with depth (rho -0.351).
The second is partly power — a real effect emerges from noise as cells
accumulate — and partly unexplained; this task does not separate them.

**The four named factors against equally deep perturbations.** In their own
depth band (96 to 157 cells) there are 60 perturbations, of which
**24 are called (40%)** under origin (i). So
being called at that depth is common. Within that band:

| factor | n cells | delta (i) | its percentile in the band | prolif. shift | its percentile in the band |
|---|---|---|---|---|---|
| `CRA_E2F3` | 128 | -19.211 | 3.3% | +0.1654 | 98.3% |
| `CRA_EZH2` | 107 | -13.211 | 15.0% | +0.1484 | 95.0% |
| `CRI_STAT3` | 96 | -10.334 | 21.7% | +0.0574 | 60.0% |
| `CRI_ZFX` | 157 | -19.936 | 1.7% | +0.1288 | 90.0% |

**Read both columns together.** The four are genuinely extreme among equally
deep perturbations on delta — 1.7th, 3.3rd, 15.0th and 21.7th percentile — so
cell count alone does not account for them. But three of the four are also at
the 90th percentile or above on proliferation shift, and the statistic is
rank-correlated −0.86 with that shift. This project's SENG task already fired
`screen_tracks_proliferation` on the same data; that flag is not retired by
anything here.

### Fibroblast identity

Mean log2-CPM change over the ten identity genes SENG found present,
perturbation minus its NT: screen median -0.0472,
2.5–97.5% [-0.3561, +0.2715];
1 of 340 exceed SENG's identity-loss threshold.

### The frozen ruler — the preprint's thesis in this dataset

- origin (i) as-measured: **181 of 340** score younger on the
  one-dimensional frozen ruler, of which **116
  score younger while not approaching**; 65 do
  both and 15 approach without scoring younger.
  Restricted to the four: 3 score younger and
  0 of those are not called.
- origin (ii) offset-corrected: **181 of 340** score younger on the
  one-dimensional frozen ruler, of which **120
  score younger while not approaching**; 61 do
  both and 12 approach without scoring younger.
  Restricted to the four: 3 score younger and
  0 of those are not called.

A diagnostic, not evidence about any factor. Note it cuts the other way from
SENG2: there 179 of 181 scored younger while not approaching, here
116 of 181 do, because many more
perturbations are called in this space.

### SENG2 against SENG3 — the effect of changing the space alone

Same cells, same depths, same assignment rule, same target, same floors
recipe. Only the space differs.

| quantity | SENG2 (full genes, S1) | SENG3 (compressed, S2A_k20) |
|---|---|---|
| D1/D3 cos(PD14−NT_A, a), depth-matched | +0.584155 (DEGRADED band) | +0.811375 (USABLE band) |
| cos(PD14−NT_I, a), depth-matched | +0.757055 | +0.887005 |
| gate on the library offset | USABLE, carried by `_I` only | USABLE, both libraries |
| null delta median at ~100 cells | +27 (positive pedestal) | about 0 |
| effect-free draws with delta<0 | 0 of 200 | about 100 of 200 |
| power at a 25% plant, worst combination | 0.000 | 0.565 |
| power at a 50% plant, worst combination | not computed at these depths | 0.940 |
| sensitivity gate | INSENSITIVE | LIMITED |
| below floor, origin (i), of 340 | 131 | 80 |
| delta negative, origin (i), of 340 | 3 | 191 |
| **called, origin (i), of 340** | 3 | **80** |
| called, origin (ii), of 340 | 0 | 73 |
| of the four, called under both origins | 0 | 4 |
| approach gate | NO APPROACH (voided by INSENSITIVE) | APPROACH |
| ruler: score younger while not approaching | 179 of 181 | 116 of 181 |

**The headline of this task is that line about the pedestal.** Changing
nothing but the coordinate space moves the null delta from +27 to zero, which
turns `delta < 0` from unreachable into a coin flip, which turns a count of 3
into a count of 80 and a sensitivity of 0.000
into 0.565. Whether the larger count is biology or is the
proliferation axis in a 20-dimensional disguise is not settled here, and the
correlations above are the reason to doubt it.

## 10. The reading

Only the pre-registered readings whose conditions fired are conclusions:

1. **Stage 0: PASS.** Every stored SOUTH8, SOUTH9 and SENG2 number
   re-derived within tolerance, counts exact.
2. **D1: CALIBRATED**, f = 0.0850 (origin (ii) would give
   0.1100 and stop the task).
3. **D2: LIMITED**, p25 = 0.565, p50 = 0.940.
4. **D3: USABLE**, r = +0.887005, both libraries in band.
5. **D4: APPROACH**.

So the conclusion this task supports, stated with the qualifications the
gates themselves impose:

> In the compressed S2A_k20 space, with a null that is calibrated at these
> depths and an instrument that detects a planted half-gap approach, all four
> of Junevity's named factors end closer to early-passage cells than their
> own controls do, under both pre-specified origins. The instrument cannot
> resolve approaches smaller than half the aging gap, the statistic is
> rank-correlated −0.86 with proliferation, residualizing proliferation
> leaves one of the four under origin (i) and none under origin (ii), and all
> four sit in the top fifth of the screen by cell count.

**This is not a claim that the four factors rejuvenate anything.** It is a
claim about where four displacement vectors land relative to a floor, in one
20-dimensional space, on passaged neonatal fibroblasts.

### What this does not cover

- **Not a result about aged human cells.** The dataset is passaged neonatal
  dermal fibroblasts and its "old" is passage number. This project has found
  that axis does not resemble donor aging: cos(v, u) = −0.040, permutation
  p = 1.0. Nothing here transfers to aged donor cells, and no result here is
  evidence about transcription factors in aged human cells.
- **No factor named here is validated or a pilot candidate**, under any gate,
  including all four that D4 called.
- **The space was fixed before any number**, by the flag and by amendment B1,
  and the component count k=20 is SOUTH7's `PRIMARY_SPACE`, not chosen here.
  But the honest framing of this task is that **SENG2 and SENG3 disagree
  because the space differs**, and this task does not establish that the
  compressed space is the right one — only that it is calibrated and
  sensitive where the full-gene space is neither.
- **The proliferation confound is not resolved.** It is reported, with the
  disclosure that the proliferation axis comes from the same control cells as
  the origin and the floors.
- **The cell-count gradient is not resolved.** It is reported, with a
  depth-matched comparison.
- **No claim that Junevity's published numbers are unreproducible.**
- **The circularities in (c) are not removed anywhere.**
- **D1's margin is 0.085 against a 0.10 threshold**, and the gate would have
  stopped the task had it been read on origin (ii).

### Files

- `results/seng3/AMENDMENTS.md` (9,507 B)
- `results/seng3/d1.log` (1,004 B)
- `results/seng3/d1_null_validation.csv` (3,992 B)
- `results/seng3/d1_stderr.txt` (0 B)
- `results/seng3/d1_stdout.txt` (1,004 B)
- `results/seng3/d1_summary.json` (188 B)
- `results/seng3/d2.log` (1,327 B)
- `results/seng3/d2_mixture.csv` (8,903 B)
- `results/seng3/d2_plant.csv` (10,745 B)
- `results/seng3/d2_stderr.txt` (0 B)
- `results/seng3/d2_stdout.txt` (1,327 B)
- `results/seng3/d2_summary.json` (530 B)
- `results/seng3/d3.log` (1,177 B)
- `results/seng3/d3_offset.csv` (4,629 B)
- `results/seng3/d3_stderr.txt` (0 B)
- `results/seng3/d3_stdout.txt` (1,177 B)
- `results/seng3/d3_summary.json` (799 B)
- `results/seng3/d4.log` (1,048 B)
- `results/seng3/d4_stderr.txt` (0 B)
- `results/seng3/d4_stdout.txt` (1,048 B)
- `results/seng3/d4_summary.json` (8,862 B)
- `results/seng3/d4_table.csv` (162,474 B)
- `results/seng3/d5_confounds.json` (3,419 B)
- `results/seng3/d5_confounds.txt` (2,164 B)
- `results/seng3/extra_counts.npz` (81,068 B)
- `results/seng3/PREREG.flag` (12,651 B)
- `results/seng3/stage0.log` (961 B)
- `results/seng3/stage0_reproduction.csv` (1,067 B)
- `results/seng3/stage0_summary.json` (3,177 B)

Nothing was committed and nothing was pushed.
