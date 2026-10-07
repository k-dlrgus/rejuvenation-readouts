# FINDINGS_SENG2 — do Junevity's four factors approach young, once the
# library structure is accounted for?

Task SENG2. Computational, no new data. Reads `results/seng/`,
`results/seng_why/`, `results/senginv/` and the Sengstack h5ad already on
disk. Writes `results/seng2/` and this file. Nothing committed, nothing
pushed. No earlier findings file or results directory was written.

**Headline, stated once and then supported below.** The gate that governs
every factor-level number in this file is D1, and it fired
**USABLE**. D3 fired **INSENSITIVE**. D2's own gate fired
**NO APPROACH**. What those three together license is set out in
"The reading" at the end, and nowhere else.

---

## 1. The pre-registration, verbatim

Reproduced verbatim from `results/seng2/PREREG.flag`, written
2026-10-06T06:01:09Z before any SENG2 number existed and not edited since.

```
PREREG written 2026-10-06T06:01:09Z (local 2026-10-06 15:01:09 KST)
Written before any SENG2 number exists. Reproduced verbatim at the top of
FINDINGS_SENG2.md. Not edited after this timestamp.

TASK: SENG2 — do Junevity's four factors approach young, once the library
structure is accounted for?
Written before any SENG2 number exists.

Nature: computational. No new data. Reads results/seng/, results/seng_why/,
results/senginv/ and the Sengstack h5ad already on disk
(data/raw/seng/exp3_merged_not_normalized.h5ad, md5
bfccaa94a916c22e4b0793913a04128b). Writes only results/seng2/ and
FINDINGS_SENG2.md. No earlier findings file is regenerated, no earlier
results directory is written, nothing is committed or pushed.

--- THE QUESTION ---

SENG ran the gated test and recorded n_toward_young = 0 across 340
perturbations, but its reading is reproduction_pass = False: only 1 of 4
named hits reproduced Junevity's own direction result, and all four were
labelled not_interpreted. SENG_WHY tested four explanations for that failure
(metric, reference, gene space, file) and fired `unexplained`.

SENGINV established a structure neither task had: the cells sit in three
barcode groups. All 6,809 passage cells (PD14/PD26/PD32, all WT) carry
suffix -1 and zero guide reads; all 14,333 CRISPRa cells carry _A; all 15,372
CRISPRi cells carry _I. No cell population appears in more than one group.
The SENGINV addendum adds two measured signs that the groups were processed
differently: the stored n_counts/n_genes disagree with .X on 7,113 rows, of
which 2,349 are _A and 4,764 are _I and 0 are -1; and median n_counts is
9,140 in PD32 against 19,611 in _A and 21,867 in _I.

The aging axis is therefore fitted entirely inside one group while every
perturbation displacement is measured in another. If the between-group offset
is large relative to the aging gap, the vector from a perturbation's own
control to the young target is mostly library difference, a genuine
anti-aging displacement is near-orthogonal to it, and no approach can
register whatever the biology. That is the same failure the preprint
discloses for its own GTEx-anchored control, in a different dataset. This
task measures the offset, gates on it, and only then reports the four.

--- STATE THE METHOD BEFORE USING IT ---

Report, before any SENG2 number:
 (a) exactly what coordinate space SENG used — what the basis is, what it was
     fit on, whether any anchor cell entered that fit, and whether it is
     refit here or read frozen. If it was fit on cells the test is applied
     to, say so plainly and label every SENG2 number as carrying that
     circularity. Do not refit to remove it; report it.
 (b) exactly how SENG assigned a cell to a perturbation, including the 804
     vs 798 guide-column question and the 270 modality disagreements it
     dropped.
 (c) the documented source for the claim that the NT control cells are
     late-passage cells of the same line as PD32. Cite the paper, by
     section. If no source can be given, STOP `nt_passage_unsourced` — the
     offset estimate in D1 depends on it entirely.

--- STAGE 0: REPRODUCE FIRST (gate) ---

Re-derive with this task's code, from stored files:
 (a) SENG's axis numbers: cos(v, p_CRA) = +0.5667444435453722,
     cos(v, p_CRI) = +0.6337418620832961, cos(v, u) = -0.04027293769593313,
     ||v|| = 90.14934904214671;
 (b) SENG's positive control: frac50 = 0.5469054675795995,
     delta50 = -34.052728101022154, frac0 = 0.033645346601332594,
     delta0 = +20.63359395062885, n_total 291;
 (c) SENG's screen counts: CRA n=167 with n_toward_young 0 and n_cos_q 0;
     CRI n=173 with n_toward_young 0 and n_cos_q 25;
 (d) the four named hits' cell counts under SENG's own rule: CRA_E2F3 128,
     CRA_EZH2 107, CRI_STAT3 96, CRI_ZFX 157.
Pass: counts exact, numerics within 1e-4 * max(1, |stored|).
Fail -> STOP `prior_results_not_reproduced`, naming every mismatch.

--- D1. How much of the target direction is library offset? (gated) ---

Let a = centroid(PD14) - centroid(PD32), the aging axis, the quantity the
whole test is scaled against.
Let offset_A = centroid(NT cells in _A) - centroid(PD32), and offset_I the
same for _I. Under the assumption sourced in (c) above, NT and PD32 are the
same biology in different libraries, so each offset estimates its library's
displacement.

Depth is confounded with group, so compute every centroid twice: raw, and
with cells subsampled to a common median n_counts across the groups entering
that comparison, 200 draws, seed 20260914. Report both. Where they disagree,
the depth-matched version is primary and that is said.

Report:
 - ||offset_A|| / ||a|| and ||offset_I|| / ||a||;
 - cos(offset_A, a) and cos(offset_I, a);
 - the angle between (PD14 - NT_A) and a, and between (PD14 - NT_I) and a.
   This is the gate quantity: it is how far the apparent target direction
   from a perturbation's own origin has been rotated away from aging.

Gate, wording fixed now. Let r be the larger of the two cosines between
(PD14 - NT) and a, over the two libraries, depth-matched.
 - USABLE (r >= 0.70): "the library offset rotates the target direction by
   less than 45 degrees (cos <r>); the apparent target direction is
   dominated by aging and the approach test can be read."
 - DEGRADED (0.40 <= r < 0.70): "the library offset rotates the target
   direction substantially (cos <r>). D2 is run and reported, and every
   number it produces is labelled as carrying a target direction that is
   <pct>% library difference."
 - UNUSABLE (r < 0.40): "the vector from a perturbation's control to the
   young target is mostly library difference, not aging (cos <r>). This
   dataset cannot answer whether these factors approach young, and SENG's
   zero is a statement about its own geometry rather than about the
   factors."
If UNUSABLE fires, D2 and D3 are still run to the end and reported, and
nothing in them is a conclusion about any factor. That is said in the
gate sentence and repeated wherever a factor is named.

--- D2. The approach test, run two ways (gated) ---

For every perturbation, and separately for the four named factors, compute
the change in distance to the PD14 target, with a depth-matched noise floor
built at each perturbation's own cell count from NT draws in its own library
(200 draws, the floor at the 5th percentile), exactly as the Southard work
builds floors. Report the chance baseline beside every count.

Two pre-specified origins, both reported, neither privileged:
 (i) as-measured: origin is the NT centroid of the perturbation's own
     library;
 (ii) offset-corrected: the perturbation-minus-NT displacement is applied to
     the PD32 origin, i.e. the library offset of D1 is subtracted. This
     transplants a displacement onto a different origin within the same cell
     line and experiment; label it as such wherever it appears.

Report for each: count below floor, count with distance decreased, count
called approaching (both), against the 5% chance baseline with its binomial
interval; and the four factors individually, with cell count and its
percentile in the screen.

Gate, wording fixed now, applied to the four named factors under whichever
origin D1's reading licenses:
 - APPROACH (two or more of the four called under both origins): "of
   Junevity's four validated factors, <k> end closer to early-passage cells
   than their own controls do, in their own data, with their own
   perturbations: <list>."
 - NO APPROACH (fewer than two under both origins): "none of Junevity's
   four validated factors measurably approaches early-passage cells in their
   own data, with their own perturbations, against a floor built at each
   factor's own depth."
 - SPLIT (the two origins disagree): report both, name which factor is which
   under each, and draw no single conclusion.

--- D3. Could the statistic have seen an approach across the gap? (gated) ---

SENG's positive control ran entirely inside the -1 group and so never tested
sensitivity across a library boundary. A null is not a null unless the
statistic could have seen the effect.

Split PD14 in half by the stored SENG split (A 290, B 291). Half A is the
target. Half B is the spike material.
 (a) Mixture control across the gap: mix half B into NT_A cells, and
     separately into NT_I cells, at 0/10/25/50%, at fixed total cell count,
     and ask whether progress rises monotonically and distance falls. The
     zero-dose point is the noise floor and is reported, not subtracted.
 (b) Plant recovery: displace NT cells along a by 0/25/50/100% of ||a||,
     200 draws per dose, at 100 cells and at 160 cells, the depths of the
     four factors.

Gate, wording fixed now:
 - SENSITIVE (mixture monotonic with the 50% interval excluding zero, AND
   power at a 25% plant >= 0.80 at 100 cells): "the statistic detects a
   known approach across the library boundary at the depth of these
   factors, power <p>, so a null in D2 is a null and not an insensitive
   assay."
 - INSENSITIVE (either condition fails): "the statistic cannot reliably
   detect a 25% approach across the library boundary at <n> cells, power
   <p>. D2's counts cannot be read as evidence of absence, and that applies
   to SENG's zero as well."
If INSENSITIVE fires, D2's reading is reported but carries no conclusion,
and that is said in the gate sentence.

--- D4. Report-only, no decision rule ---
 - Proliferation: the stored rho_prolif (+0.828 CRA, +0.862 CRI) beside
   every result, and a proliferation-residualized version of D2 computed
   with the null residualized identically. Reported beside the distance
   result, never in place of it.
 - Identity: the fibroblast identity score change for the four and for the
   screen, using the ten identity genes SENG already found present.
 - The frozen ruler applied to the four: how many of the 340 perturbations
   score younger on a one-dimensional ruler while not approaching. This is
   the preprint's thesis in this dataset and is a diagnostic, not evidence
   about any factor.
 - The SENG_WHY finding that recomputed R_rej is roughly half the published
   value for all four, restated with its numbers, labelled post-hoc and
   unexplained. No claim is made about it in this task.

--- WHAT DOES NOT COUNT ---
- choosing a space, an origin, a criterion, a null or a factor set after
  seeing a number;
- treating any factor called here as validated or as a pilot candidate in
  this findings file;
- reading a null as evidence about transcription factors in aged human
  cells: this dataset is passaged neonatal-line fibroblasts and its "old" is
  passage number, which this project has already found does not resemble
  donor aging (cos(v,u) = -0.040, p = 1.0);
- refitting any axis or offset to avoid reporting a circularity;
- any claim about Junevity's published numbers being unreproducible.

--- STOP KEYS ---
prior_results_not_reproduced; nt_passage_unsourced; no_raw_counts.
Seeds: 20260914 splits/nulls/permutations, 20260918 bootstraps.
n_null = n_perm = n_boot = 200.
```

## 2. Amendments, and every deviation

The flag was found ambiguous or self-contradictory in four places. Those
four were put to the author **before computing**, not resolved by the
session; five further definitions the flag does not supply at all were
fixed in writing before computing. All of it is in
`results/seng2/AMENDMENTS.md`, reproduced here in full.

```
# SENG2 — amendments to the pre-registration

Written 2026-10-06T06:14Z, **before any SENG2 number existed**, after the flag
was found ambiguous or self-contradictory in four places and the four were put
to the author rather than resolved by this session. `results/seng2/PREREG.flag`
is **not edited**; it stands as written at 2026-10-06T06:01:09Z. These are
additions to the record, and every one of them is reproduced in
`FINDINGS_SENG2.md`.

Nothing below was chosen by this session. A1–A4 are the author's answers to
four questions asked before computing. A5–A9 are definitions the flag does not
supply at all; those were fixed in writing here, before computing, by this
session, and each is labelled as such.

---

## A1 — D2's gate (self-contradiction; author's answer)

**The contradiction.** D2's preamble applies the gate "under whichever origin
D1's reading licenses", but none of D1's three readings names an origin. The
branches are then written over both origins and overlap: with 3 factors called
under origin (i) and 0 under origin (ii), NO APPROACH ("fewer than two under
both") and SPLIT ("the two origins disagree") fire together with contradictory
instructions.

**Resolution, as answered: SPLIT only on branch disagreement.** The clause
"under whichever origin D1's reading licenses" is dropped, because D1 licenses
none; both origins are co-primary, as D2's own "both reported, neither
privileged" already says. Calls are counted separately per origin. Then:

- **SPLIT** fires if and only if the two origins reach different verdicts (one
  origin calls >= 2 of the four, the other calls < 2).
- **APPROACH** fires if and only if both origins call >= 2.
- **NO APPROACH** fires if and only if both origins call < 2.

The three are now mutually exclusive and exhaustive.

## A2 — What a centroid is (missing definition; author's answer)

D1 writes "centroid(PD14)", "centroid(PD32)", "centroid(NT cells in _A)"
without saying in what space, or whether a centroid is a group pseudobulk or a
mean over per-cell vectors.

**Resolution, as answered: SENG's pseudobulk-z path, read frozen.** A centroid
of a group of cells is: sum the raw counts over the group's cells, scatter the
sum into the 23,485 frozen-ruler columns (16,781 present, the rest left at
zero), run TMM jointly over the analysis panel, take log2-CPM, and z-score
against the frozen GTEx mu and sd. This is the instrument SENG's stored numbers
live in, it is read and not refit, and it makes every SENG2 quantity
commensurable with SENG's `v`.

## A3 — Depth matching (ambiguous; author's amendment, verbatim)

D1 says "with cells subsampled to a common median n_counts across the groups
entering that comparison, 200 draws, seed 20260914", which reads either as
selecting cells by depth or as downsampling UMIs. The author amended the flag.
Quoted verbatim:

> Amendment, disclosed as a deviation from the flag's literal wording:
> depth matching is done by binomial thinning, not by selecting cells.
> For each comparison, take the lowest group median n_counts among the
> groups entering it as the target. Thin every cell in each other group
> by one common factor per group (target median / that group's median),
> so the group median matches the target. Keep every cell; no cell is
> dropped or selected on depth. 200 draws, seed 20260914. The gate value
> r is the median over the 200 draws, with its 2.5-97.5 percentile range
> reported beside it. The raw (unthinned) version is still reported as the
> flag requires, with the thinned version primary. Reason for the
> amendment: in late-passage cells, low read depth is partly biological
> (smaller, non-dividing cells), so selecting cells by depth would
> preferentially remove the most aged cells from PD32 and bias the
> aging axis itself.

**This is a deviation from the flag's literal wording and is reported as one**
wherever a depth-matched number appears.

## A4 — Which PD14 is the target (ambiguous; author's answer)

D1 defines `a` over the full arms; D3 makes PD14 half A the target and half B
the spike material; D2 says only "the PD14 target". Reading D2's target as the
full arm would put D3's spike material inside D2's target.

**Resolution, as answered: `a` over the full arms; PD14 half A is the target in
D2 and D3.** `a = centroid(PD14, all 581) - centroid(PD32, all 4,348)`, as the
flag writes it. The target in D2 and D3 is `centroid(PD14 half A, 290)`, by the
stored SENG split. Half B (291) never appears in any target.

---

The remainder were not put to the author. Each is a definition the flag does
not supply, fixed in writing here before computing, by this session.

## A5 — The statistic (missing definition; fixed by this session)

D2 asks for "the change in distance to the PD14 target" and does not define it.
It is taken to be SENG's own `delta`, read frozen from
`src/seng_run.py:stats_from_d`: for a target vector `t` and a displacement `d`,

    delta = ||t - d|| - ||t||

negative meaning the group ends closer to the target than the origin was.
`progress` is SENG's `frac`, the component of `d` along `t` as a fraction of
`||t||`. Both are the functions SENG's stored `screen.csv` columns were computed
with; neither is re-derived here.

## A6 — The two origins, written out (missing definition; fixed by this session)

For a perturbation `g` in library `L` (A for CRISPRa, I for CRISPRi), with
`z(.)` the centroid of A2:

- displacement, both origins: `d = z(g) - z(NT_L)`
- origin (i), as-measured: `t_i = z(PD14 half A) - z(NT_L)`
- origin (ii), offset-corrected: `t_ii = z(PD14 half A) - z(PD32 half A)`

PD32 **half A** is used as the offset-corrected origin, not the full PD32 arm,
so that both anchors of `t_ii` come from the same stored split as the target.
One consequence is checkable and is reported as a check: `t_ii` is then
identical to SENG's frozen `v` in `results/seng/axes.npz`, so origin (ii) is
literally SENG's geometry and origin (i) is the perturbation's own.

## A7 — D1's undefined percentage (missing definition; fixed by this session)

D1's DEGRADED sentence contains "<pct>% library difference" without a formula.
It is taken to be the fraction of the apparent target direction orthogonal to
aging, `sqrt(1 - r^2)`, reported as a percentage. D1's USABLE sentence says
"less than 45 degrees" while its own threshold `r >= 0.70` corresponds to
45.57 degrees; the exact angle `arccos(r)` is reported beside `r` wherever the
gate is quoted, so the sentence's rounding is visible.

## A8 — D3(a)'s unspecified cell count (missing definition; fixed by this session)

D3(a) mixes at "fixed total cell count" without naming it. The mixture is run at
**160 cells and again at 100 cells**, the two counts D3(b) names as "the depths
of the four factors". Monotonicity and the 50% interval are assessed at 160
cells, with 100 cells reported beside them.

## A9 — D3's "power" (missing definition; fixed by this session)

D3(b) asks for "power at a 25% plant" without defining the event. Power is the
fraction of the 200 draws at that dose that D2 would have called approaching,
i.e. `delta < floor` and `delta < 0`, with the floor built exactly as in D2 at
the same cell count. Intervals on the mixture doses are the 2.5-97.5 percentile
range of the 200 draws, not a separate bootstrap.

---

## Implementation identity used for thinning, and why it is exact

A3 thins every cell of a group by one common factor `p_g`. Thinning is binomial
and independent across cells, and independent binomials with a common `p` sum to
a binomial of the sum, so the thinned **group pseudobulk** is distributed exactly
as `Binomial(pseudobulk, p_g)`. Thinning is therefore applied to the group
pseudobulk rather than cell by cell. This is an exact distributional identity
given a common per-group factor, not an approximation, and nothing downstream of
the pseudobulk sees individual cells. The realised per-cell median after thinning
is checked directly on one draw and reported.
```

### Deviations beyond the amendments

1. **The bootstrap confidence intervals of SENG's screen were not
   re-derived.** Stage 0 checks no CI, and neither `meets_toward` nor
   `n_cos_q` reads one: both depend only on the nulls and the point
   estimate. The `*_ci_*` columns of `results/seng2/stage0_screen.csv` are
   therefore all NaN. This was a deliberate scope reduction taken after the
   first attempt showed the bootstrap would take about four hours and feed
   nothing that is checked. SENG's own CIs stand unexamined by this task.
2. **`.X` is streamed, not loaded.** About 2.6 GB of RAM was free against a
   1.4 GB matrix, so group pseudobulks are accumulated in 2,000-row blocks
   and only the few thousand cells needing per-cell resampling are
   materialised. This changes no number; the panel reproduces SENG's.
3. **Amendment A10, fixed in writing before computing.** Groups that are
   not in SENG's 348-row panel (the full PD14 and PD32 arms, the half-B
   spike material, every mixture and plant draw) enter the frozen space
   through `S.z_from_counts` against the panel's frozen TMM cache. The
   panel is **not** extended and nothing is refit. This is SENG's own
   mechanism for every null and bootstrap draw it computed, and SENG
   verified it agrees with the panel path to under 1e-4.
4. **This session reused SENG's implementation.** Stage 0 calls the frozen
   helper functions in `src/seng_run.py` (`z_from_counts`, `stats_from_d`,
   `bh_q`, `meets_toward`, `parse_wide_obs`, the TMM path) from a SENG2
   driver. A reproduction that reuses the original implementation tests the
   record, the data path and the driver — **not** the implementation. An
   independent reimplementation would be a stronger check and was not done.
5. The positive control was reimplemented in `src/seng2_run.py` rather than
   calling `S.run_posctrl`, because that function writes into
   `results/seng/`, which this task may not do.

## 3. State the method before using it

### (a) The coordinate space, and its circularity

The space is the **23,485-column frozen GTEx fibroblast ruler z space**, the
project's frozen instrument, read and not refit. A group of cells becomes a
point in it by: summing raw counts over the group's cells; scattering the sum
into the 23,485 ruler columns (16,781 matched, by gene symbol only — the file
carries no Ensembl IDs, so 0 matched by Ensembl — and the remaining 6,704 left
at z=0); one joint TMM over the analysis panel; log2-CPM with edgeR's prior;
and z against the frozen GTEx mu and sd.

The basis is not fitted at all — it is the frozen ruler's gene axes. What **is**
fitted on cells is:

- `v = z(EARLY_A) - z(LATE_A)`, SENG's aging axis, from PD14 half A (290
  cells) and PD32 half A (2,174 cells) — `src/seng_run.py:1836`.
- the proliferation axis `p` per modality, from the top and bottom NT
  proliferation quartiles.
- the joint TMM normalisation, fitted over **one panel containing EARLY_A,
  LATE_A, both NT groups, all four proliferation-quartile groups and all 340
  perturbations** — `src/seng_run.py:1795`.

**The circularity, stated plainly and not refitted away.** The anchor cells
that define the target direction also co-determine the normalisation that
every displacement is measured through, because they are rows of the same TMM
panel. Separately, D1's `a` and D2's target are built from PD14 and PD32 cells,
and D1's gate quantity is a cosine between two vectors that both contain
`z(PD14)` — so the gate is not independent of the anchor estimate. **Every
number in D1, D2, D3 and D4 of this file carries that circularity.** SENG
held PD14/PD32 **half B** out of the panel and used it for the positive
control, which is why half B is the spike material in D3 and never appears in
a target (amendment A4). Nothing was refit to remove any of this.

### (b) How SENG assigned a cell to a perturbation

Reproduced exactly by this task's code (`src/seng_run.py:648`,
`parse_wide_obs`), with the counts below re-derived here, not quoted:

- the pool is **804** guide columns: 798 targeting columns
  `<TARGET>_<CRA|CRI><k>` plus `NT_1`..`NT_6`. The six NT columns **are** in
  the argmax pool. SENGINV's 798 excluded them and is a different question;
  SENG's 804 is the rule that produced the stored numbers and is the rule
  used here.
- a non-WT cell is assigned the guide column that is a **unique strict
  maximum and > 0**. Pooled identity is `{CRA|CRI}_{TF}`.
- **2,680 cells have a tied positive maximum** and are given no winner.
- **270 cells** whose winning targeting column disagrees with the cell's own
  `CRA`/`CRI` indicator are **not assigned** — they are dropped to
  `unassigned`, not relabelled and not reassigned to the indicator's
  direction. 2,680 + 270 = 2,950 unassigned cells, 8.1% of the file.
- roles re-derived here: WT 6,809, CRA 12,264, CRI 12,346, NT 2,145,
  unassigned 2,950. Unique-argmax NT cells: CRA 666, CRI 1,479.
- a perturbation enters the panel at `n >= 30` cells: **340** qualify,
  167 CRA and 173 CRI, out of 400 target x direction pairs.

**Recorded as a limitation of the rule, not corrected.** The 2,950 dropped
cells are dropped non-randomly: a tie is likelier at low guide UMI, and the
270 disagreements are exactly the cells where the two encodings of direction
conflict — the conflict SENGINV already recorded as an ambiguity the file
does not resolve.

### (c) The documented source for NT being late-passage, same line as PD32

**STOP key `nt_passage_unsourced` does not fire.** A source can be given.

Sengstack et al., PNAS 2026, Results, section *"Top TF Perturbations
Identified By TRDP Reversed Several Cellular Aging Hallmarks in Passaged
Fibroblasts"*, defines the passage bands verbatim:

> "Based on phenotyping we performed (data in Fig. 2), early passage cells
> were 0 to 20 population doublings (PD)."

> "middle passage cells were 21 to 30 PD, and late passage cells were over
> 31 PDs."

and names the controls in the same section: *"Nontargeting sgRNA were used as
controls and are noted as either CRA NT (for CRISPRa cell line with
nontargeting sgRNA)"* / *"or CRI NT ... with nontargeting sgRNA"*. The same
Results text states the perturbations were applied to late-passage cells
(*"On average, the TF perturbations caused late passage cells to have cellular
phenotypes more like middle passage cells"*). PD32 is over 31 PD, so it is
late passage by the paper's own cut.

**Two gaps in that source, carried into D1 rather than papered over.**

1. The paper gives a **band**, not a number. It never states the PD at which
   the Perturb-seq cells were harvested, only that they were over 31 PD. PD32
   sits at the very bottom of that band. **`offset_A` and `offset_I` are
   therefore upper bounds on library offset**: any real passage difference
   between the screen cells and PD32 is inside them and cannot be separated.
2. The same section says the change *"occurred within about three population
   doublings (duration of the experiment)"*, so the perturbed cells drifted in
   passage during the experiment. The h5ad carries no PD label for any
   perturbed cell (SENG: perturbed PD values empty, 24,610 cells blank), so
   none of this is checkable in the data.

Also: the CRISPRa and CRISPRi lines are engineered derivatives of the
fibroblast line, not the line itself, and the object gives them separate
barcode groups. "Same line" is the paper's framing, adopted here as the
flag instructs, and it is an assumption, not a measurement.

## 4. Stage 0 — reproduce first (gate)

`no_raw_counts` does not fire. Full scan of all 148,239,518 stored
values of `X/data`: 0 fractional, 0 negative,
0 non-finite, min 1, max 2960,
verdict raw integer counts in float32. (This also corrects the SENGINV
inventory, which reported max 2,446 from a 10M-value sample; the full-scan
maximum is 2,960, as SENG already had. That file is not edited.)

**Gate: PASS.** Tolerance as pre-registered: counts exact, numerics within 1e-4 x max(1, |stored|).

| item | kind | stored | re-derived | abs diff | ok |
|---|---|---|---|---|---|
| `cos_v_p_CRA` | numeric | 0.5667444435453722 | 0.5667444435453722 | 0.000e+00 | yes |
| `cos_v_p_CRI` | numeric | 0.6337418620832961 | 0.6337418620832961 | 0.000e+00 | yes |
| `cos_v_u` | numeric | -0.0402729376959331 | -0.0402729376959331 | 0.000e+00 | yes |
| `v_norm` | numeric | 90.14934904214672 | 90.14934904214672 | 0.000e+00 | yes |
| `frac50` | numeric | 0.5469054675795995 | 0.5469054675795995 | 0.000e+00 | yes |
| `delta50` | numeric | -34.05272810102215 | -34.05272810102215 | 0.000e+00 | yes |
| `frac0` | numeric | 0.0336453466013325 | 0.0336453466013325 | 0.000e+00 | yes |
| `delta0` | numeric | 20.63359395062885 | 20.63359395062885 | 0.000e+00 | yes |
| `posctrl_n_total` | count | 291 | 291 | 0 | yes |
| `CRA_n` | count | 167 | 167 | 0 | yes |
| `CRA_n_cos_q` | count | 0 | 0 | 0 | yes |
| `CRA_n_toward_young` | count | 0 | 0 | 0 | yes |
| `CRI_n` | count | 173 | 173 | 0 | yes |
| `CRI_n_cos_q` | count | 25 | 25 | 0 | yes |
| `CRI_n_toward_young` | count | 0 | 0 | 0 | yes |
| `n_cells[CRA_E2F3]` | count | 128 | 128 | 0 | yes |
| `n_cells[CRA_EZH2]` | count | 107 | 107 | 0 | yes |
| `n_cells[CRI_STAT3]` | count | 96 | 96 | 0 | yes |
| `n_cells[CRI_ZFX]` | count | 157 | 157 | 0 | yes |

- `v` re-derived from the data against SENG's stored `axes.npz` `v`:
  max abs difference **0.000e+00** over 23,485
  components. The axis is the same vector, not merely the same summary.
- cos(v, u) permutation p re-derived: 1 (stored 1.0).
- the positive control reproduced to every printed digit of all four stored
  numbers, and the screen counts and all four named-hit cell counts are exact.

## 5. D1 — how much of the target direction is library offset? (gated)

`a = centroid(PD14, 581) - centroid(PD32, 4,348)`; `offset_L = centroid(NT_L)
- centroid(PD32)`; the gate quantity is the cosine between `(PD14 - NT_L)`
and `a`. Depth matching is amendment A3: binomial thinning to the lowest
group median, 200 draws, seed 20260914, thinned primary.

Group median `n_counts`, which set the thinning factors:

| group | n cells | median n_counts |
|---|---|---|
| PD14 (all) | 581 | 19,853.0 |
| PD32 (all) | 4348 | 9,139.5 |
| NT_A | 666 | 19,605.5 |
| NT_I | 1479 | 22,220.0 |
| EARLY_A (PD14 half A) | 290 | 20,330.5 |
| LATE_A (PD32 half A) | 2174 | 8,759.5 |

- comparison A: target median 9,139.5, factors `PD14_FULL`=0.4604, `PD32_FULL`=1.0000, `NT_CRA`=0.4662
- comparison I: target median 9,139.5, factors `PD14_FULL`=0.4604, `PD32_FULL`=1.0000, `NT_CRI`=0.4113

### Results, raw and depth-matched

| quantity | raw | depth-matched median | 2.5–97.5% |
|---|---|---|---|
| \|\|a\|\| | +82.867150 | +89.851446 | +89.141752 to +90.522181 |
| \|\|offset_A\|\| / \|\|a\|\| | +1.245203 | +1.202091 | +1.190028 to +1.215352 |
| \|\|offset_I\|\| / \|\|a\|\| | +0.840694 | +0.817736 | +0.810239 to +0.824738 |
| cos(offset_A, a) | +0.138934 | +0.117405 | +0.108863 to +0.126275 |
| cos(offset_I, a) | +0.079530 | +0.067211 | +0.057857 to +0.076246 |
| **cos(PD14 - NT_A, a)** | +0.556990 | +0.584155 | +0.576564 to +0.589674 |
| **cos(PD14 - NT_I, a)** | +0.744005 | +0.757055 | +0.753166 to +0.760620 |

Angles, the same numbers read as rotations:

| | raw | depth-matched |
|---|---|---|
| angle between (PD14 - NT_A) and a | 56.15° | 54.26° |
| angle between (PD14 - NT_I) and a | 41.93° | 40.79° |

- raw `r` = +0.744005 (gate on raw alone would be **USABLE**)
- **depth-matched `r` = +0.757055**, from library `I`,
  2.5–97.5% [+0.753166, +0.760620], angle 40.79°
- the two agree on the branch, so A3's "where they disagree" clause is not
  invoked.
- cos(a, v) = +0.918322: the full-arm axis `a` against SENG's
  half-A axis `v`. They are not the same vector, and D1 is scaled against `a`
  as the flag defines it.

**Why thinning moves these numbers so little, stated so the reader is not
misled.** log2-CPM already divides out library size, and the prior count
scales with the library, so thinning a pseudobulk by a common factor changes
the z vector mainly by adding sampling noise rather than by removing a depth
effect. A3's thinning is therefore a weak control for the depth confound, not
a strong one, and the raw and thinned numbers being close is not evidence that
depth does not matter. Per-cell median after thinning is in
`results/seng2/d1_thinning_check.json`.

### The two libraries are not in the same condition — read this before the gate

The flag fixes `r` as **the larger** of the two cosines, so the gate is
decided by whichever library is better aligned. The two are not close, and
they fall in different pre-registered bands:

| library | cos(PD14 - NT_L, a), depth-matched | angle | band that cosine sits in | named factors in it |
|---|---|---|---|---|
| `_A` (CRISPRa) | +0.584155 | 54.26° | DEGRADED | CRA_E2F3, CRA_EZH2 |
| `_I` (CRISPRi) | +0.757055 | 40.79° | USABLE | CRI_STAT3, CRI_ZFX |

**This is a diagnostic, not a second gate** — the flag fixed one gate on the
maximum and that is what is reported below. But it bears directly on two of
the four named factors: `CRA_E2F3` and `CRA_EZH2` sit in the library whose
own cosine is +0.584155, inside the DEGRADED band, where the
pre-registered wording would have said the target direction is rotated
substantially. The gate that fired was carried by the CRISPRi library.

The offsets themselves are large and nearly orthogonal to aging, which is how
both of those cosines stay as high as they are:

- `||offset_A||` is **1.20x** the whole aging gap
  `||a||`, and cos(offset_A, a) = +0.117405.
- `||offset_I||` is **0.82x** `||a||`, with
  cos(offset_I, a) = +0.067211.

So the between-library displacement is comparable to or larger than the entire
PD14-to-PD32 aging gap, and it points almost at right angles to it. A large
near-orthogonal offset shortens the apparent target direction's alignment with
aging without destroying it; that is the geometry the gate is measuring.

### The D1 gate

**USABLE** fired, on the pre-registered thresholds, with
`r` = +0.757055.

> "the library offset rotates the target direction by less than 45
> degrees (cos +0.757055); the apparent target direction is dominated by
> aging and the approach test can be read."

The flag's "less than 45 degrees" is the flag's rounding: `r >= 0.70`
is 45.57°, and the measured angle is 40.79° (amendment A7).

**One qualification on the word "dominated", which the flag's wording
does not carry.** A cosine of
+0.757055 means the component of the apparent target direction orthogonal to
aging is still **65.3%** of its length
(`sqrt(1 - r^2)`). Aging is the single largest component and the gate's
threshold is met, but a clear majority of that vector's length is not
aging. Every D2 number carries that.

## 6. D2 — the approach test, run two ways (gated)

**Carries D1: the target direction is dominated by aging.**

Origins, per amendment A6. Target is PD14 half A throughout (A4).
Displacement `d = z(g) - z(NT_L)` under both.

- **origin (i), as-measured**: `t_i = z(PD14 half A) - z(NT_L)`
- **origin (ii), offset-corrected**: `t_ii = z(PD14 half A) - z(PD32 half A)`
  — this **transplants the displacement onto a different origin within the
  same cell line and experiment**, and `t_ii` is identical to SENG's frozen
  `v`, so origin (ii) is literally SENG's geometry.

Floors are the Southard construction, confirmed identical to SENG's null
(`src/south7_run.py:330`): draw `n` control cells from the perturbation's own
library, displacement `z(draw) - z(pool minus draw)`, 200 draws, floor at the
5th percentile. The null depends only on `n`, so one null per distinct cell
count. A perturbation is **called approaching** when its delta is below its
own floor **and** negative.

### Counts over all 340 perturbations, with the chance baseline

The floor is the 5th percentile of each perturbation's own null, so **5% of
true nulls fall below it by construction**: the chance expectation is 17.0 of
340. A count near 17 is what a screen with no real effect produces.

| origin | below floor | distance decreased | **called approaching** | chance expectation | 95% interval on the called fraction |
|---|---|---|---|---|---|
| (i) as-measured | 131 | 3 | **3** | 17.0 | 0.2–2.6% |
| (ii) offset-corrected | 56 | 0 | **0** | 17.0 | 0.0–1.1% |

Per modality:

| origin | modality | n | below floor | decreased | called | chance exp. |
|---|---|---|---|---|---|---|
| (i) | CRA | 167 | 73 | 1 | 1 | 8.3 |
| (i) | CRI | 173 | 58 | 2 | 2 | 8.7 |
| (ii) | CRA | 167 | 24 | 0 | 0 | 8.3 |
| (ii) | CRI | 173 | 32 | 0 | 0 | 8.7 |

**Which half of the criterion is binding, because it decides everything
below.** "Called approaching" is the conjunction of two conditions, and they
are not equally hard to meet here:

- under origin (i), **131 of 340** perturbations fall below
  their own floor, but only **3** have a negative delta.
- under origin (ii), **56** fall below floor and
  **0** have a negative delta.

The floor is easy to clear and `delta < 0` is almost impossible, because at
these cell counts the null delta is **positive**: sampling noise alone moves a
group *away* from the target. D3 measures that pedestal directly — at 100
cells the median null delta is about +27 — and it is why the counts below are
near zero. A count near zero here is a statement about the pedestal at least
as much as about the perturbations, which is exactly what D3 then tests.

### The four named factors

Cell counts are SENG's own rule (b). Percentile is within the 340-perturbation
screen by cell count.

| factor | n cells | pctl | delta (i) | floor (i) | called (i) | delta (ii) | floor (ii) | called (ii) |
|---|---|---|---|---|---|---|---|---|
| `CRA_E2F3` | 128 | 96% | +0.9257 | +14.3142 | no | +21.5961 | +24.3164 | no |
| `CRA_EZH2` | 107 | 89% | +8.3140 | +18.7303 | no | +24.4261 | +28.2714 | no |
| `CRI_STAT3` | 96 | 82% | +16.1092 | +20.0593 | no | +28.0071 | +26.4950 | no |
| `CRI_ZFX` | 157 | 99% | -0.3374 | +12.5131 | **yes** | +13.1248 | +16.6274 | no |

- called under origin (i): **1 of 4**
- called under origin (ii): **0 of 4**

### The D2 gate

**NO APPROACH** fired, on amendment A1's mutually exclusive branches.

> "none of Junevity's four validated factors measurably approaches
> early-passage cells in their own data, with their own perturbations,
> against a floor built at each factor's own depth."

## 7. D3 — could the statistic have seen an approach across the gap? (gated)

SENG's positive control ran entirely inside the `-1` barcode group and so
never crossed a library boundary. Spike material is PD14 **half B** (291
cells), which appears in no target.

`||a||` = 82.8672, the full aging gap the plant doses scale.

### (a) Mixture across the library boundary

Half B mixed into NT cells at 0/10/25/50% at fixed total cell count (160 and
100, amendment A8), 200 draws per dose. The zero-dose point is the noise
floor and is reported, not subtracted.

| n total | library | dose | progress (frac, origin i) | 2.5–97.5% | delta (origin i) | 2.5–97.5% |
|---|---|---|---|---|---|---|
| 100 | CRA | 0% | +0.0099 | -0.0691 to +0.0861 | +27.160 | +19.918 to +35.442 |
| 100 | CRA | 10% | +0.0979 | +0.0050 to +0.1727 | +18.050 | +10.362 to +28.120 |
| 100 | CRA | 25% | +0.2328 | +0.1381 to +0.3331 | +5.378 | -3.079 to +14.318 |
| 100 | CRA | 50% | +0.4364 | +0.3469 to +0.5232 | -11.640 | -17.392 to -3.749 |
| 100 | CRI | 0% | +0.0152 | -0.0874 to +0.1182 | +26.318 | +17.790 to +36.901 |
| 100 | CRI | 10% | +0.0777 | -0.0342 to +0.1987 | +21.616 | +11.230 to +31.736 |
| 100 | CRI | 25% | +0.1973 | +0.0787 to +0.3012 | +11.862 | +3.937 to +21.103 |
| 100 | CRI | 50% | +0.3919 | +0.2772 to +0.4829 | -0.317 | -6.148 to +8.595 |
| 160 | CRA | 0% | -0.0011 | -0.0715 to +0.0788 | +19.896 | +12.274 to +28.046 |
| 160 | CRA | 10% | +0.0932 | +0.0271 to +0.1686 | +9.645 | +1.747 to +17.058 |
| 160 | CRA | 25% | +0.2256 | +0.1575 to +0.2979 | -4.058 | -10.760 to +2.332 |
| 160 | CRA | 50% | +0.4344 | +0.3698 to +0.5086 | -23.189 | -29.126 to -17.513 |
| 160 | CRI | 0% | +0.0084 | -0.0798 to +0.0930 | +17.454 | +10.254 to +26.536 |
| 160 | CRI | 10% | +0.0834 | +0.0038 to +0.1742 | +10.848 | +3.811 to +19.310 |
| 160 | CRI | 25% | +0.1836 | +0.0879 to +0.2904 | +2.120 | -5.941 to +11.729 |
| 160 | CRI | 50% | +0.3866 | +0.2902 to +0.4622 | -12.965 | -18.344 to -4.979 |

Monotonicity of progress, and whether the 50% dose's interval excludes zero:

| key | progress monotonic | delta(50%) interval excludes 0 |
|---|---|---|
| 100_CRA_i | yes | yes |
| 100_CRA_ii | yes | **no** |
| 100_CRI_i | yes | **no** |
| 100_CRI_ii | yes | **no** |
| 160_CRA_i | yes | yes |
| 160_CRA_ii | yes | **no** |
| 160_CRI_i | yes | yes |
| 160_CRI_ii | yes | **no** |

### (b) Plant recovery along `a`

NT draws displaced along `a` by 0/25/50/100% of `||a||`, 200 draws per dose,
at 100 and 160 cells. Power is the fraction of draws D2 would have called
(amendment A9).

| n cells | library | dose | delta median (i) | floor (i) | **power (i)** | power (ii) |
|---|---|---|---|---|---|---|
| 100 | CRA | 0% | +27.160 | +20.801 | **0.000** | 0.000 |
| 100 | CRA | 25% | +19.090 | +20.801 | **0.000** | 0.000 |
| 100 | CRA | 50% | +13.658 | +20.801 | **0.000** | 0.000 |
| 100 | CRA | 100% | +11.440 | +20.801 | **0.000** | 0.005 |
| 100 | CRI | 0% | +26.318 | +19.528 | **0.000** | 0.000 |
| 100 | CRI | 25% | +15.815 | +19.528 | **0.000** | 0.000 |
| 100 | CRI | 50% | +7.876 | +19.528 | **0.010** | 0.005 |
| 100 | CRI | 100% | +2.405 | +19.528 | **0.110** | 0.520 |
| 160 | CRA | 0% | +19.896 | +12.949 | **0.000** | 0.000 |
| 160 | CRA | 25% | +11.312 | +12.949 | **0.000** | 0.000 |
| 160 | CRA | 50% | +5.256 | +12.949 | **0.030** | 0.235 |
| 160 | CRA | 100% | +2.943 | +12.949 | **0.065** | 1.000 |
| 160 | CRI | 0% | +17.454 | +11.606 | **0.000** | 0.000 |
| 160 | CRI | 25% | +6.008 | +11.606 | **0.040** | 0.010 |
| 160 | CRI | 50% | -2.673 | +11.606 | **0.725** | 0.945 |
| 160 | CRI | 100% | -8.490 | +11.606 | **1.000** | 1.000 |

### The D3 gate

**INSENSITIVE** fired. Inputs: mixture monotonic at 160 cells under origin
(i) = True; 50% interval excludes zero = True; minimum power at a 25% plant at 100 cells = 0.000 against the 0.80 requirement.

> "the statistic cannot reliably detect a 25% approach across the library
> boundary at 100 cells, power 0.000. D2's
> counts cannot be read as evidence of absence, and that applies to SENG's
> zero as well."

**D2's reading is reported but carries no conclusion.**

## 8. D4 — report-only, no decision rule

Nothing in this section is gated and nothing in it decides anything.

### Proliferation

SENG's stored rank correlation between the screen statistic and
proliferation is **rho = +0.828 (CRA)** and **+0.862
(CRI)**. It is quoted beside every result in this file, as the flag requires.
A proliferation-residualized D2, with the null residualized identically by
projecting out the same per-modality proliferation axis:

| origin | called, plain delta | called, residualized | called by both |
|---|---|---|---|
| (i) | 3 | 0 | 0 |
| (ii) | 0 | 0 | 0 |

Reported beside the distance result, never in place of it.

### Fibroblast identity

Mean log2-CPM change over the 10 identity genes SENG found
present (COL1A1, COL1A2, FN1, LUM, PDGFRA, PDGFRB, POSTN, PRRX1, SERPINH1, VIM), perturbation minus its NT:
screen median -0.0472, 2.5–97.5%
[-0.3561, +0.2715]; 1 of 340 exceed SENG's identity-loss threshold.

### The frozen ruler — the preprint's thesis in this dataset

- origin (i) as-measured: **181 of 340** perturbations score
  younger on the one-dimensional frozen ruler, of which
  **179 score younger while not approaching**;
  2 both score younger and approach, and
  1 approach without scoring younger.
- origin (ii) offset-corrected: **181 of 340** perturbations score
  younger on the one-dimensional frozen ruler, of which
  **181 score younger while not approaching**;
  0 both score younger and approach, and
  0 approach without scoring younger.

This is a **diagnostic**, not evidence about any factor: it restates in this
dataset the dissociation the preprint is about. It inherits every
circularity in (a) and the D1 gate.

### The four factors, report-only columns

| factor | n cells | identity change | ruler change | younger on ruler | prolif. shift | SENG cos | SENG q_cos |
|---|---|---|---|---|---|---|---|
| `CRA_E2F3` | 128 | -0.2154 | -1.0060 | yes | +0.1654 | +0.2479 | 0.0923 |
| `CRA_EZH2` | 107 | -0.3298 | -0.6785 | yes | +0.1484 | +0.1973 | 0.0923 |
| `CRI_STAT3` | 96 | -0.1477 | -0.9329 | yes | +0.0574 | +0.1546 | 0.1001 |
| `CRI_ZFX` | 157 | -0.0380 | +0.5116 | no | +0.1288 | +0.3165 | 0.0344 |

### R_rej, restated post-hoc and unexplained

SENG_WHY recomputed the authors' own statistic and got roughly half the
published value for all four. Restated with its numbers, **labelled post-hoc
and unexplained. No claim is made about it in this task, and in particular no
claim that Junevity's published numbers are unreproducible.**

| factor | n cells | R_rej recomputed | published Table 1 | ratio |
|---|---|---|---|---|
| `CRA_E2F3` | 128 | -0.2809 | -0.53 | 0.53 |
| `CRA_EZH2` | 107 | -0.2030 | -0.36 | 0.56 |
| `CRI_STAT3` | 96 | -0.1549 | -0.41 | 0.38 |
| `CRI_ZFX` | 157 | -0.3553 | -0.51 | 0.70 |

SENG_WHY's own `assumptions` block records six choices made after numbers
existed (gene selection, Pearson, sign convention, the log2FC recipe, which
PDs count as aging, and the z-scoring of unmapped genes). The ratio above is
not corrected for any of them.

## 9. The reading

Only the pre-registered readings whose conditions fired are reported as
conclusions. They are, in the order the flag gates them:

1. **Stage 0: PASS.** Every stored SENG number re-derived within tolerance,
   counts exact.
2. **D1: USABLE**, `r` = +0.757055.
3. **D2: NO APPROACH**.
4. **D3: INSENSITIVE**.

D1 licenses reading D2, but **D3 fired INSENSITIVE**, so D2's counts
cannot be read as evidence of absence. D2's sentence is reported above
and carries no conclusion.

### What this does not cover

- **Not a result about aged human cells.** This dataset is passaged neonatal
  dermal fibroblasts and its "old" is passage number. This project has
  already found that axis does not resemble donor aging: cos(v, u) =
  -0.040, permutation p = 1.0, re-derived in Stage 0. Nothing here transfers
  to aged donor cells, and no null here is evidence about transcription
  factors in aged human cells.
- **No factor named in this file is validated or a pilot candidate**, under
  any of the gates, including any factor D2 called.
- **No claim that Junevity's published numbers are unreproducible.** The
  R_rej comparison in D4 is post-hoc, unexplained, and explicitly not a claim.
- **The circularity in (a) is not removed anywhere.** The anchors co-determine
  the normalisation; D1's gate quantity shares `z(PD14)` with the target it
  judges. Nothing was refit to hide it.
- **`offset_A` and `offset_I` are upper bounds**, because the paper sources the
  screen cells only to a band (>31 PD) and not to PD32, and the h5ad carries
  no PD label for any perturbed cell.
- **Depth matching is weak.** See the note in D1: log2-CPM already divides out
  library size, so thinning mostly adds noise. The depth confound between the
  barcode groups is reported, not removed.
- **SENG's bootstrap CIs were not re-derived** (deviation 1).
- **Stage 0 reuses SENG's implementation** (deviation 4), so it tests the
  record and the data path, not the implementation.

### Files

- `results/seng2/AMENDMENTS.md` (7,987 B)
- `results/seng2/d1.log` (1,952 B)
- `results/seng2/d1_stderr.txt` (0 B)
- `results/seng2/d1_stdout.txt` (1,952 B)
- `results/seng2/d1_summary.json` (3,123 B)
- `results/seng2/d1_thinned_draws_A.csv` (22,792 B)
- `results/seng2/d1_thinned_draws_I.csv` (22,718 B)
- `results/seng2/d1_thinning_check.json` (878 B)
- `results/seng2/d2.log` (933 B)
- `results/seng2/d2_stderr.txt` (0 B)
- `results/seng2/d2_stdout.txt` (933 B)
- `results/seng2/d2_summary.json` (6,881 B)
- `results/seng2/d2_table.csv` (112,248 B)
- `results/seng2/d3.log` (889 B)
- `results/seng2/d3_mixture_draws.csv` (419,608 B)
- `results/seng2/d3_mixture_summary.csv` (4,161 B)
- `results/seng2/d3_plant_draws.csv` (418,082 B)
- `results/seng2/d3_plant_summary.csv` (1,717 B)
- `results/seng2/d3_stderr.txt` (0 B)
- `results/seng2/d3_stdout.txt` (889 B)
- `results/seng2/d3_summary.json` (910 B)
- `results/seng2/d4.log` (150 B)
- `results/seng2/d4_merged.csv` (276,076 B)
- `results/seng2/d4_summary.json` (4,614 B)
- `results/seng2/extra_counts.npz` (198,437 B)
- `results/seng2/panel.log` (936 B)
- `results/seng2/panel_counts.npz` (9,349,406 B)
- `results/seng2/panel_medians.json` (8,605 B)
- `results/seng2/PREREG.flag` (11,043 B)
- `results/seng2/raw_audit.json` (204 B)
- `results/seng2/stage0.log` (2,078 B)
- `results/seng2/stage0_posctrl.csv` (572 B)
- `results/seng2/stage0_reproduction.csv` (924 B)
- `results/seng2/stage0_screen.csv` (166,242 B)
- `results/seng2/stage0_space.npz` (47,247,445 B)
- `results/seng2/stage0_stderr.txt` (0 B)
- `results/seng2/stage0_stdout.txt` (2,078 B)
- `results/seng2/stage0_summary.json` (3,500 B)

Nothing was committed and nothing was pushed.
