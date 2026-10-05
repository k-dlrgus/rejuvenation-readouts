# FINDINGS_SOUTH9 — is the below-floor excess proliferation?

Status: Stage 0 passed (11 checks). D1, D2 and D3 ran; no STOP fired.

D1 gate: **EXPLAINED**. D2 gate: **RETAINED** (power at a 25% plant after residualization 0.990, against 0.980 raw).

## Pre-registration (verbatim, `results/south9/PREREG.flag`)

```
TASK: SOUTH9 — is the below-floor excess proliferation?
Written before any SOUTH9 number exists.

Nature: computational. No new data. Reads results/south6/, results/south7/,
results/south8/ and the Southard singlets file already on disk. Writes only
results/south9/ and FINDINGS_SOUTH9.md. No earlier findings file is
regenerated and no earlier results directory is written. Nothing is committed
or pushed.

--- THE QUESTION ---

SOUTH8's primary reading, after the permutation condition was removed as
uncalibrated, is 24 of 304 and 125 of 1,836 factors called APPROACHING in
S2A_k20 under the amended criterion, against 15.2 and 91.8 expected by chance
at a 5th-percentile floor. SOUTH8 established that this excess is larger than
chance and is not explained by cell count or by total reads. It did not
establish what it is.

SOUTH6 already computes a proliferation-residualized displacement
(`delta_resid`, `cos_resid`) alongside the raw one. In the SOUTH6 basis, across
all 1,836 factors, the count with a negative raw change in distance is 108 and
the count with a negative residualized change in distance is 51 — below the 5%
a chance baseline would give. That is a diagnostic from a stored table, in the
wrong basis, with no floor attached. This task tests it properly.

A rejuvenation readout that cannot separate approach to young from a change in
proliferative state is the confound the preprint names, and the preprint's own
remedy list asks for a proliferation-residualized version to be reported beside
the headline number. This task runs that remedy at screen scale.

--- STATE THE METHOD BEFORE USING IT ---
Report, before any SOUTH9 number: exactly how SOUTH6 constructs the
proliferation axis and the residual — what the axis is, what data it was fit
on, whether any anchor donor or any Southard cell entered that fit, and whether
the axis is refit here or carried over frozen. If the axis was fit on the same
Southard cells the test is applied to, say so plainly and label every SOUTH9
number as carrying that circularity. Do not refit the axis to remove the
problem; report it.

--- STAGE 0: REPRODUCE FIRST (gate) ---
Re-derive with this task's code, from stored files:
 (a) SOUTH8's amended counts in S2A_k20: 24 of 304, 125 of 1,836, and the
     below-floor counts 27 and 153;
 (b) SOUTH8's D2 calibration row for S2A_k20 at 106 cells (0.0350) and S1 at
     106 cells (0.1200);
 (c) the 24-factor list in primary_304.
Pass: counts exact, numeric items within 1e-4 * max(1, |stored|).
Fail -> STOP `prior_results_not_reproduced`, naming every mismatch.

--- D1. The residualized statistic in the anchor-free space (gated) ---

Build the proliferation-residualized change in distance in S2A_k20, the SOUTH8
primary space, for all 1,836 factors and for the published 304.

The null must be residualized identically. Every depth-matched non-targeting
draw from SOUTH7 Stage A1 is put through the same residualization as a factor,
at the same cell count, and the 5th percentile is recomputed from the
residualized draws. A raw floor applied to a residualized statistic is not a
comparison and must not be reported as one.

Report, per factor set, in S2A_k20 and in the three secondary spaces:
 - below the residualized depth-matched floor;
 - change in residualized distance < 0;
 - called APPROACHING under SOUTH8's amended criterion (both of the above);
 - the same three counts for the raw statistic, side by side;
 - which of SOUTH8's 24, and of its 125, survive.

Gate, wording fixed now. 5% of 304 is 15.2 with a 95% binomial interval of 8 to
23; 5% of 1,836 is 91.8 with an interval of 74 to 110.
 - EXPLAINED (called count falls inside or below the interval, in S2A_k20, on
   both factor sets): "the excess of factors below their own noise floor is
   accounted for by proliferation: residualized, <k> of <m> are called, against
   <e> expected by chance. What the raw statistic scores as approach to young
   is in this screen a change in proliferative state."
 - SURVIVES (above the interval on either set): "the excess survives
   proliferation residualization: <k> of <m> are called against <e> expected.
   Proliferation does not account for it."
 - SPLIT (one set inside, one above): report both counts and draw no single
   conclusion; name which set is which.
Both the raw and residualized counts are reported whichever fires.

--- D2. Does residualization also remove a real approach? (gated) ---
Residualization can remove genuine signal as well as confound. The plant
recovery control answers whether it does.
Repeat SOUTH6REC's plant recovery in S2A_k20 at 106 cells, 200 draws per dose,
doses 0/25/50/100% of the measured gap, with the residualization applied to
plant and null alike.
Gate, wording fixed now:
 - RETAINED (power at a 25% plant >= 0.80): "residualization leaves the
   statistic able to detect a 25% approach, power <p>, so a null result after
   residualization is a null result and not a loss of sensitivity."
 - LOST (power < 0.80): "residualization costs the statistic its sensitivity:
   power at a 25% plant falls from 0.980 to <p>. The residualized counts cannot
   be read as evidence of absence."
If LOST fires, D1's reading is reported but carries no conclusion, and that is
said in the gate sentence.

--- D3. Report-only, no decision rule ---
 - For the factors called under the residualized criterion: cell count and its
   percentile, the high_depth_candidate flag, identity drop, guide agreement,
   and the frozen ruler score change.
 - The raw-versus-residualized change in distance for SOUTH8's 24 named
   factors, individually.
 - Identity: SOUTH8's 24 have a median identity drop of -0.0122 against +0.0408
   for the screen, so they lose less fibroblast identity than average. Report
   the same comparison for whichever set survives here.
 - The count of factors that score younger on the frozen ruler while not being
   called under the residualized criterion, both factor sets. This is the
   preprint's thesis in this screen and is a diagnostic, not evidence about any
   factor.

--- WHAT DOES NOT COUNT ---
- choosing a space, a criterion, a null or a factor set after seeing a number;
- treating any surviving factor as a candidate for the pilot panel in this
  findings file;
- reading a residualized null as evidence that no transcription factor can
  rejuvenate aged cells: the screen uses neonatal foreskin fibroblasts with no
  aged donor and no matched young target;
- any result from an O2 pair;
- refitting the proliferation axis to avoid reporting a circularity.

--- STOP RULES ---
On any STOP: write the findings with the key and exit. Never compute after a
STOP; never resume a stopped run.

--- OUTPUT ---
FINDINGS_SOUTH9.md: this text verbatim; the method statement; the Stage 0
table; D1 and D2 with the gate reading that fired; D3; limitations. Tables in
results/south9/. Code src/south9_run.py. Hand over files, not summaries. Do not
regenerate any earlier findings file. Do not commit or push.
```

## 1. The method, stated before it is used

How SOUTH6 builds the proliferation axis (`src/south6_run.py`, the Stage 1 block that writes `p_unit`):

1. **Gene set.** Seurat's `cc.genes` (Tirosh 2016), parsed from the packaged `.rda` with boundary checks (97 genes: 43 S-phase, 54 G2M; the parser refuses a hand-written substitute). 94 of the 97 are present in the Southard panel.
2. **Score.** For each **non-targeting Southard cell**, the mean of log1p(counts/library x 1e4) over those cc genes.
3. **Contrast.** The non-targeting cells are split at the quartiles of that score: 3,504 cells in the bottom quartile, 3,504 in the top. Each group is summed into a pseudobulk, carried through the same TMM / log2-CPM / frozen-GTEx-z transform as every other pseudobulk, and the axis is z(top) - z(bottom), norm 65.9633, stored normalised to unit length as `p_unit`.
4. **Residual.** SOUTH6's `stats_from_d` removes the component along `p_unit` from both the factor's displacement and the origin-target vector, then recomputes the change in distance.

What entered the fit, and what did not:

- **No GTEx donor entered it.** The axis is built only from Southard cells, so neither anchor (O1, Y1) nor the frozen ruler's training set contributed to its direction.
- **No targeting cell entered it.** Only non-targeting cells were used, so a factor's own cells did not help define the axis it is tested against.
- **The axis is carried over frozen here and is not refit.** It is read from `results/south6/stage1.npz` (`p_unit`), unit norm, 18,849 nonzero coordinates, all of them inside the 20,264 matched genes.

> **The circularity, stated plainly.** The axis was fit on the same 14,013 non-targeting cells that (i) form the reference pseudobulk every factor's displacement is measured against, and (ii) generate every depth-matched null draw used here, including the residualized floors. The axis is therefore not independent of the data it is applied to: it is the leading cell-cycle contrast *within the very pool that defines zero*. Removing that direction removes whatever part of a factor's displacement happens to align with the dominant cell-cycle contrast of the control pool, and that contrast is itself estimated from a finite sample of those cells. Every SOUTH9 number below carries this circularity. The pre-registration forbids refitting the axis to avoid it, and it was not refit.

**One choice the pre-registration left open, fixed here before any number was computed.** SOUTH6's residual is defined in gene space; the SOUTH8 primary space is a 20-component projection, and removing a direction in gene space does not leave the projection free of it. Two readings are therefore possible, and both are reported:

- `in_space` (**primary**): the axis is projected into the space, and its component is removed from the displacement and the origin-target vector *inside that space*. This is the reading that makes the residualized statistic actually free of the proliferation direction as the space sees it.
- `gene_space` (secondary): SOUTH6's own formula - remove in gene space, then project.

In S1 the two coincide. How much of the axis each space retains, and how aligned it is with the young-old direction:

| space | ||proj p|| (axis retention) | cos(origin->target, axis) |
|---|---|---|
| S2A_k20 | 0.5433 | 0.3384 |
| S2A_k50 | 0.5923 | 0.2646 |
| S2A_k100 | 0.6020 | 0.2531 |
| S1 | 1.0000 | 0.1298 |

The cosine matters for reading D2: the more the proliferation axis aligns with the young-old direction, the more of a genuine approach residualization must also remove.

## 2. Stage 0: SOUTH8's numbers re-derived (gate)

11 checks, 11 gating, 0 failures. Full table: `results/south9/stage0_reproduction.csv`.

| block | item | n | mismatches | passed |
|---|---|---|---|---|
| (a) SOUTH8 amended primary_304 S2A_k20 | called | 1 | 0 | yes |
| (a) SOUTH8 amended primary_304 S2A_k20 | below depth floor | 1 | 0 | yes |
| (a) SOUTH8 amended all_1836 S2A_k20 | called | 1 | 0 | yes |
| (a) SOUTH8 amended all_1836 S2A_k20 | below depth floor | 1 | 0 | yes |
| (a) SOUTH8 amended | called, primary_304 | 1 | 0 | yes |
| (a) SOUTH8 amended | called, all_1836 | 1 | 0 | yes |
| (a) SOUTH8 amended | below floor, primary_304 | 1 | 0 | yes |
| (a) SOUTH8 amended | below floor, all_1836 | 1 | 0 | yes |
| (b) SOUTH8 D2 calibration S2A_k20 at 106 cells | fraction p <= 0.05 | 1 | 0 | yes |
| (b) SOUTH8 D2 calibration S1 at 106 cells | fraction p <= 0.05 | 1 | 0 | yes |
| (c) the 24-factor list | primary_304 amended call list | 24 | 0 | yes |

Nothing in `results/south6/`, `results/south7/` or `results/south8/` was written, and no earlier findings file was regenerated.

## 3. D1: the residualized statistic in the anchor-free space (gated)

Every depth-matched non-targeting draw from SOUTH7 Stage A1 is put through the same residualization as a factor, at the same cell count, and the 5th percentile is recomputed from the residualized draws: 212 distinct cell counts, 200 draws each, per space and per variant. No raw floor is applied to a residualized statistic anywhere in this file.

Primary reading (`in_space`), raw and residualized side by side:

| factor set | space | m | raw: below floor | raw: delta<0 | raw: called | residualized: below floor | residualized: delta<0 | residualized: called |
|---|---|---|---|---|---|---|---|---|
| primary_304 | S2A_k20 | 304 | 27 | 35 | 24 | 22 | 24 | 15 |
| primary_304 | S2A_k50 | 304 | 23 | 4 | 4 | 17 | 1 | 1 |
| primary_304 | S2A_k100 | 304 | 20 | 1 | 1 | 14 | 0 | 0 |
| primary_304 | S1 | 304 | 3 | 0 | 0 | 5 | 0 | 0 |
| all_1836 | S2A_k20 | 1836 | 153 | 214 | 125 | 133 | 139 | 95 |
| all_1836 | S2A_k50 | 1836 | 130 | 21 | 19 | 120 | 13 | 13 |
| all_1836 | S2A_k100 | 1836 | 118 | 7 | 6 | 105 | 4 | 4 |
| all_1836 | S1 | 1836 | 48 | 0 | 0 | 49 | 0 | 0 |

Secondary reading (`gene_space`, SOUTH6's own formula):

| factor set | space | m | raw: called | residualized: below floor | residualized: delta<0 | residualized: called |
|---|---|---|---|---|---|---|
| primary_304 | S2A_k20 | 304 | 24 | 26 | 15 | 12 |
| primary_304 | S2A_k50 | 304 | 4 | 20 | 1 | 1 |
| primary_304 | S2A_k100 | 304 | 1 | 17 | 0 | 0 |
| primary_304 | S1 | 304 | 0 | 5 | 0 | 0 |
| all_1836 | S2A_k20 | 1836 | 125 | 133 | 94 | 69 |
| all_1836 | S2A_k50 | 1836 | 19 | 124 | 10 | 10 |
| all_1836 | S2A_k100 | 1836 | 6 | 107 | 2 | 2 |
| all_1836 | S1 | 1836 | 0 | 49 | 0 | 0 |

Which of SOUTH8's calls survive, in the primary space and primary variant:

| factor set | SOUTH8 called (raw) | called residualized | in both | lost | called residualized but NOT by SOUTH8 | in both |
|---|---|---|---|---|---|---|
| primary_304 | 24 | 15 | 8 | 16 | 7 | BNC2;EGR2;EGR3;FOXN3;KLF12;MYRF;PLAGL1;SOX8 |
| all_1836 | 125 | 95 | 54 | 71 | 41 | ATMIN;BAZ2B;BNC2;BTAF1;CRTC1;CTBP1;DBX2;DMRTA1;DPRX;EGR2;EGR3;FOXJ3;FOXN3;FUBP1;GABPB1;GLMP;GRHL2;GSC2;GSX1;GTF3C2;HOXA10;HOXD1;KLF12;MAEL;MAFB;MEF2B;MYRF;NFATC4;OVOL3;PLAGL1;POU5F1;RAD21;RELA;SMAD9;SNAPC2;SOX8;SP9;SPIB;THAP1;TLE1;XRCC4;ZFP90;ZNF184;ZNF222;ZNF224;ZNF268;ZNF275;ZNF286B;ZNF473;ZNF502;ZNF561;ZNF780B;ZNF84;ZSCAN1 |

**The residualized call set is not a subset of the raw one.** In the primary set, 8 of SOUTH8's 24 survive, but 7 of the 15 residualized calls are factors SOUTH8 did not call: BATF3, ESR2, HELZ2, ILF3, NEUROD1, OSR2, ZNF844. Removing a direction moves some factors below their floor as well as lifting others above it, so '15 called' must not be read as '15 of the 24'.

### Gate reading that fired (pre-registered wording)

- primary_304: 15 of 304 called residualized, against 15.2 expected by chance (95% interval 8 to 23) — inside or below the interval.
- all_1836: 95 of 1836 called residualized, against 91.8 expected by chance (95% interval 74 to 110) — inside or below the interval.

**EXPLAINED.** The excess of factors below their own noise floor is accounted for by proliferation: residualized, 15 of 304 and 95 of 1,836 are called, against 15.2 and 91.8 expected by chance. What the raw statistic scores as approach to young is in this screen a change in proliferative state.

## 4. D2: does residualization also remove a real approach? (gated)

SOUTH6REC's plant recovery, repeated in S2A_k20 at 106 cells, 200 draws per dose, doses 0/25/50/100% of the measured gap, with the residualization applied to plant and null alike.

| space | raw power 25% | residualized power 25% | residualized power 50% | residualized power 100% | residualized median at 0% | at 25% | monotone |
|---|---|---|---|---|---|---|---|
| S2A_k20 | 0.9800 | 0.9900 | 1.0000 | 1.0000 | 1.9763 | -3.5908 | yes |
| S2A_k50 | 0.9550 | 0.9750 | 1.0000 | 1.0000 | 4.1685 | -1.2100 | yes |
| S2A_k100 | 0.9300 | 0.9500 | 1.0000 | 1.0000 | 5.5107 | 0.2001 | yes |
| S1 | 0.1150 | 0.1150 | 0.3300 | 0.4650 | 93.7560 | 92.3390 | yes |

### Gate reading that fired (pre-registered wording)

- **RETAINED.** Residualization leaves the statistic able to detect a 25% approach, power 0.990, so a null result after residualization is a null result and not a loss of sensitivity.

The planted gap and the proliferation axis are not orthogonal, so some of a genuine approach lies along the direction being removed; the cosine in section 1 says how much.

## 5. D3: report-only (no decision rule)

| factor set | m | called residualized | of those, high_depth_candidate | identity drop, median (called) | identity drop, median (set) | guide agreement (called) | guide agreement (set) |
|---|---|---|---|---|---|---|---|
| primary_304 | 304 | 15 | 1 | 0.0292 | 0.0518 | 0.1798 | 0.1935 |
| all_1836 | 1836 | 95 | 14 | 0.0284 | 0.0408 | 0.1798 | 0.1993 |

Identity, the comparison the pre-registration asks for: SOUTH8's 24 have a median identity drop of -0.0122 against +0.0408 for the screen, so they lose less fibroblast identity than average. The same comparison for the set that survives here is in the table above (a negative drop means identity genes rose relative to the non-targeting pool).

The frozen ruler against the distance call - the preprint's thesis in this screen, a diagnostic and not evidence about any factor:

| factor set | m | score younger on the ruler | score younger yet NOT called residualized |
|---|---|---|---|
| primary_304 | 304 | 126 | 117 |
| all_1836 | 1836 | 761 | 716 |

SOUTH8's 24, raw against residualized, individually:

| factor | cells | high_depth | raw delta | raw floor | residualized delta | residualized floor | still called | identity drop | guide agreement | ruler score change |
|---|---|---|---|---|---|---|---|---|---|---|
| FOXN3 | 128 | no | -2.1206 | -0.8317 | -1.6828 | -0.2534 | yes | -0.0151 | 0.1635 | -0.4337 |
| EGR3 | 71 | no | -2.3610 | 0.9528 | -1.1886 | 1.2539 | yes | -0.0765 | 0.2597 | -0.1754 |
| SOX8 | 114 | no | -1.8734 | -0.7587 | -1.0162 | -0.1778 | yes | 0.0244 | 0.1999 | -1.8477 |
| KLF12 | 117 | no | -1.5769 | -0.5573 | -0.8074 | -0.1281 | yes | -0.0353 | 0.1884 | 0.3469 |
| AFF1 | 177 | no | -1.8316 | -1.2665 | -0.7721 | -0.8633 | no | -0.0457 | 0.1441 | -0.1530 |
| MYRF | 130 | no | -0.9079 | -0.8448 | -0.4590 | -0.3528 | yes | 0.2207 | 0.1881 | -0.2304 |
| EGR2 | 133 | no | -0.8786 | -0.6818 | -0.3540 | -0.3211 | yes | -0.0865 | 0.1784 | -0.2135 |
| BNC2 | 120 | no | -1.1147 | -0.5986 | -0.2885 | -0.1942 | yes | 0.1896 | 0.1798 | -0.0548 |
| PLAGL1 | 118 | no | -1.6283 | -0.7335 | -0.2842 | -0.2746 | yes | -0.1005 | 0.2029 | -1.2399 |
| KLF15 | 143 | no | -0.9909 | -0.4996 | -0.2341 | -0.3751 | no | -0.0325 | 0.1693 | -0.0941 |
| AFF4 | 199 | yes | -1.3419 | -1.0042 | -0.1896 | -0.5726 | no | -0.0521 | 0.1466 | 0.1591 |
| BNC1 | 152 | no | -1.7147 | -1.1506 | -0.1608 | -0.4554 | no | -0.0320 | 0.1599 | -0.1231 |
| KAT2A | 104 | no | -1.2149 | -0.2712 | -0.0859 | -0.1328 | no | 0.0223 | 0.1869 | 0.6764 |
| ZXDA | 84 | no | -0.9624 | 0.5761 | 0.3842 | 0.6672 | no | -0.0536 | 0.2671 | 0.1704 |
| CDX2 | 80 | no | -1.1307 | -0.0351 | 0.3849 | 0.5988 | no | -0.2897 | 0.2287 | 0.5348 |
| NFATC1 | 142 | no | -0.9403 | -0.7434 | 0.4011 | -0.2930 | no | 0.0181 | 0.1807 | 0.5786 |
| FOXC1 | 120 | no | -0.6239 | -0.5986 | 0.4297 | -0.1942 | no | 0.0112 | 0.2561 | 0.6902 |
| HNF1B | 84 | no | -0.0315 | 0.5761 | 0.4633 | 0.6672 | no | 0.0668 | 0.2537 | -0.3233 |
| NKX3-2 | 94 | no | -0.8888 | 0.0053 | 0.5689 | 0.3421 | no | 0.1360 | 0.1968 | 0.5134 |
| PTF1A | 125 | no | -0.8604 | -0.5991 | 0.6627 | -0.4043 | no | -0.0140 | 0.1588 | 1.3133 |
| CEBPA | 74 | no | -0.3894 | 0.7661 | 0.6686 | 1.3607 | no | 0.0056 | 0.2682 | 0.6765 |
| TSHZ1 | 77 | no | -0.5359 | 0.3335 | 0.7438 | 1.0736 | no | 0.0479 | 0.2479 | -0.2600 |
| NR2F1 | 76 | no | -0.6318 | 0.1834 | 0.7680 | 1.0831 | no | -0.0105 | 0.2276 | 0.5711 |
| HOXA13 | 96 | no | -0.7382 | -0.1069 | 0.9028 | 0.4378 | no | 0.0454 | 0.2263 | -0.1781 |

## 6. Limitations

- **The axis is not independent of the test.** It is the dominant cell-cycle contrast of the same non-targeting pool that defines the reference and every null draw (section 1). A residualized count is therefore a statement about movement orthogonal to one finite-sample estimate of that contrast, not about proliferation in general.
- Residualization removes one direction. A factor whose effect is partly proliferative and partly not keeps only the orthogonal part, so a surviving factor is not proliferation-free and a lost factor is not purely proliferative.
- The cc.genes score is a coarse proxy for proliferative state; it does not distinguish cycle arrest from senescence, quiescence or a change in cycle phase distribution.
- No factor named here is a candidate for a pilot panel, and nothing in this file proposes one.
- A residualized null is not evidence that no transcription factor can rejuvenate aged cells: the screen uses neonatal foreskin fibroblasts with no aged donor and no matched young target.
- O2 pairs were not run.
- The residualized criterion inherits SOUTH8's: a 5th-percentile floor with no multiple-testing control, since the permutation condition was removed as uncalibrated in SOUTH8 D2.
- 3,221 of the ruler's 23,485 genes are outside the panel and are assumed unmoved; the axis is supported only on the 20,264 matched genes.

## 7. Files

`src/south9_run.py`. In `results/south9/`: `PREREG.flag`, `method_statement.json`, `stage0_reproduction.csv`, `d1_residualized_floors.csv`, `d1_counts.csv`, `d1_primary_304_residualized.csv`, `d1_all_1836_residualized.csv`, `d2_residualized_recovery.csv`, `d3_summary.csv`, `d3_south8_24_raw_vs_residualized.csv`, `d3_*_called_detail.csv`, `south9_run.log`.
