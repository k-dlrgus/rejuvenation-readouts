# FINDINGS_AXANG — how far is each aging model from a proliferation axis?

Task AXANG. Computational, descriptive geometry. No new data. Writes
`results/axang/` and this file. Nothing committed, nothing pushed.

**The gates.** Stage 0 PASS; D1 **INTERMEDIATE** on the primary
construction; D3a **QUALIFIED**; D3b **NO APPROACH**.

**The one-scale answer the task was built to give.** Each aging model
against a proliferation axis built from its own units, in S2A_k20, as the
angle from parallel (amendment F10):

| aging model | its own proliferation axis | angle from parallel | % of the aging axis removed by projecting it out |
|---|---|---|---|
| A3, Sengstack passage | C1 Sengstack NT CRA (333) | **19.66°** | 66.4% |
| A1, GTEx donor age | C1 GTEx donors (652, bulk) | **35.83°** | 41.5% |
| A4, Lu donor pair | C1 Lu day-0 (12,803) | **57.63°** | 15.5% |
| A2, frozen ruler | C1 GTEx donors (652, bulk) | **85.51°** | 0.3% |

So the ordering, which is the task's result: **Sengstack passage aging is
the worst case by a wide margin** (19.66°, 66% removed), the GTEx donor-age
axis is intermediate (35.83°, 41.5% removed), and the **frozen ruler and the
preprint's own Lu donor-pair axis are close to orthogonal** to proliferation
(85.51° and 57.63°, 0.3% and 15.5% removed).

**Two findings the flag's recorded expectation did not anticipate, stated
rather than smoothed, as it requires.**

1. **The sign.** `cos(A1, C1_gtex)` = -0.810808, bootstrap
   [-0.939, -0.340] — strongly **negative**.
   GTEx donors with higher cell-cycle scores sit on the *old* side of the
   donor-age axis. The flag's expectation was near-orthogonality, and the
   magnitude is far from it.
2. **The expectation holds for cell-derived axes and fails for the bulk one.**
   A1 against the Southard NT-cell axis is -0.3384
   and against the Sengstack NT-cell axes -0.2891
   / -0.2689 — all near-orthogonal, as
   predicted. Only A1's **own bulk-donor** axis is strongly anti-aligned. The
   inconsistency is therefore between two *constructions* of proliferation,
   not between the aging models, and the bulk construction is the outlier.

---

## 1. The pre-registration, verbatim

From `results/axang/PREREG.flag`, written 2026-10-06T23:34:40Z before any
AXANG number existed and not edited since. Amendments E1-E4 were supplied
with the flag and are part of it.

```
PREREG written 2026-10-06T23:34:40Z (local 2026-10-07 08:34:40 KST)
Written before any AXANG number exists. Reproduced verbatim at the top of
FINDINGS_AXANG.md. Not edited after this timestamp. The amendments E1-E4 at
the end were supplied with the flag and are part of the pre-registration.

TASK: AXANG — how far is each aging model from a proliferation axis?
Written before any AXANG number exists.

Nature: computational, descriptive geometry. No new data. Reads
results/fibro/, results/south6/, results/south7/, results/south9/,
results/seng3/, results/seng4/ and, only if an axis cannot be read from a
stored file, the Southard and Sengstack matrices already on disk. Writes only
results/axang/ and FINDINGS_AXANG.md. Nothing is committed or pushed.

--- THE QUESTION ---

SENG4 fired D1 INSEPARABLE: the Sengstack passage-aging axis sits
19.66 degrees from a proliferation axis built from that dataset's own
non-targeting cells (cos +0.941708, CRA), and 61.68 degrees from a
proliferation axis built from the Tirosh gene list alone (cos +0.474385).
Projecting the first out removed 66% of the aging axis's length and cost the
statistic its sensitivity.

SOUTH9 fired D2 RETAINED on the Southard screen, whose aging axis is the GTEx
donor-age direction: residualizing proliferation there left power at 0.990
against 0.980 raw. So the separation evidently works against donor age and
fails against passage. But that is an inference from sensitivity, not a
measured angle, and no angle exists for the donor-age model.

This task measures the same angle, by the same constructions, for every aging
model this project uses, so the models can be compared directly on one scale.
It decides nothing about any transcription factor and names none as a
candidate.

A geometric expectation is recorded here before any number: this project's
standing finding is cos(passage axis, GTEx donor-age axis) = -0.040273, i.e.
the two aging models are near-orthogonal. Combined with SENG4's finding that
passage is near-parallel to proliferation, the donor-age axis should then be
near-orthogonal to proliferation as well. If the measured angles do not form
a consistent triangle, that inconsistency is itself reportable and must be
stated rather than smoothed.

--- STATE THE METHOD BEFORE USING IT ---

Report, before any AXANG number:
 (a) for every axis used below, the stored file it is read from, or if it
     must be rebuilt, exactly what cells or donors enter it and the recipe,
     quoting the source lines from the committed code;
 (b) which axes are read frozen and which are rebuilt here. Frozen
     instruments are read, not refit;
 (c) the space each comparison is made in, and a statement that cosines are
     only compared between axes living in the same space and dimension. Any
     cross-space number is labelled as not comparable;
 (d) whether any axis is estimated from cells that also define another axis
     in the same comparison. If so, say so plainly and label the affected
     numbers. Do not refit to hide it.

--- STAGE 0: REPRODUCE FIRST (gate) ---

Re-derive with this task's code, from stored files:
 (a) SENG4 D1: cos(P1, a) = +0.941708 (CRA) and +0.936287 (CRI);
     cos(P3, a) = +0.474385;
 (b) SENG4 D1: cos(a, u) = -0.040273 in the full-gene space and
     -0.199293 in S2A_k20;
 (c) SOUTH9's D2 residualized and raw power at 106 cells in S2A_k20:
     0.990 and 0.980;
 (d) SOUTH8's amended counts in S2A_k20: 125 of 1,836 called and 153 below
     floor.
Pass: counts exact, numerics within 1e-4 * max(1, |stored|).
Fail -> STOP `prior_results_not_reproduced`, naming every mismatch.

--- D1. The angles (gated) ---

Three aging axes, each a unit vector:
 A1, GTEx donor age: the young-minus-old centroid direction, young = GTEx
     reference donors aged 20-39 (n = 113), old = 60-79 (n = 231), read from
     the stored anchors, not rebuilt.
 A2, the frozen ruler direction: the unit-length frozen fibroblast ruler
     vector, read from results/fibro/, not refit.
 A3, Sengstack passage: PD14 half A minus PD32 half A, as SENG4 defines it.

Two proliferation axis constructions, each applied wherever the data allows:
 C1, data-derived: score the units of that dataset for cell-cycle activity
     with the Tirosh S and G2M lists at
     results/seng/seurat_cc.genes_tirosh2016.tsv, and take the direction from
     the bottom quartile to the top quartile of that score. The units are
     GTEx donors for A1 and A2, Southard non-targeting cells for the Southard
     comparison, and Sengstack non-targeting cells for A3. State n for each.
 C2, gene-list: the equal-weight unit direction over the same Tirosh genes,
     positive for cycling, built from no cell or donor at all.

Report, in S2A_k20 and again in the full-gene space where both axes exist
there, with every pair labelled by the space it lives in:
 - cos and angle in degrees between each aging axis and each proliferation
   construction available to it;
 - the fraction of the aging axis's length surviving projection of each
   proliferation axis out, sqrt(1 - cos^2), as a percentage;
 - cos and angle between every pair of aging axes;
 - cos and angle between every pair of proliferation axes, including the
   Southard-derived and Sengstack-derived C1 axes against each other;
 - n of cells or donors entering each C1 axis;
 - the three SENG4 values reproduced in Stage 0, in the same table, so the
   new numbers sit beside the known ones.

Also report, as the triangle consistency check the question section names:
the angle between A1 and A3, the angle between A3 and its C1, and the angle
between A1 and its C1, with a statement of whether the three are mutually
consistent to within the bootstrap intervals below.

Intervals: 200-draw percentile bootstraps, resampling the unit that defines
each axis (donors for A1 and A2, cells for C1 in each dataset). Report the
2.5-97.5% range on every cosine. Seed 20260918. Any interval that does not
contain its point estimate is flagged invalid and is not used as a gate.

Gate, wording fixed now. Let d be the angle between A1 and its data-derived
proliferation axis C1, in S2A_k20, and let p = 19.66 degrees be the SENG4
value for A3.
 - SEPARABLE (d >= 45): "donor aging sits <d> degrees from an independently
   measured proliferation axis, against <p> degrees for passage aging in the
   same construction. A rejuvenation signal can be distinguished from a
   change in proliferative state against donor age and cannot be against
   passage number."
 - INTERMEDIATE (25 <= d < 45): "donor aging sits <d> degrees from
   proliferation, against <p> for passage aging. Donor age is the better
   separated of the two models but proliferation still carries <pct>% of it,
   and any donor-age rejuvenation readout must report a
   proliferation-residualized version beside the headline."
 - ENTANGLED (d < 25): "donor aging sits <d> degrees from proliferation,
   close to the <p> of passage aging. Proliferation cannot be separated from
   either aging model here, and this project's own Southard proliferation
   finding requires re-reading."

The gate is read on C1, the data-derived construction, because that is the
construction SENG4's primary axis used and the comparison must be
like-for-like. If C2 would fire a different branch, say so prominently beside
the gate and report both branch sentences, with C1's as the task's reading.
If the ENTANGLED branch fires, state explicitly that SOUTH9's RETAINED
reading and SOUTH8's residualized count are placed in question by it, and
name the files, without revising them here.

--- D2. Report-only, no decision rule ---
 - Whether the Southard-derived and Sengstack-derived C1 axes agree with each
   other, and what that says about whether the recipe measures the same thing
   in two datasets.
 - The fraction of each aging axis's length that survives projecting out C2
   rather than C1, so the dependence on the construction is visible.
 - A single table collecting every angle this project has measured between an
   aging direction and a proliferation direction, across SOUTH9, SENG4 and
   this task.

--- WHAT DOES NOT COUNT ---
- choosing an axis, a space or a construction after seeing an angle;
- naming any transcription factor as a candidate, or reporting any
  per-factor statistic at all; this task is about directions, not factors;
- treating C2 as primary because it gives a preferred answer;
- reading any angle measured in GTEx bulk donors as a statement about single
  cells, or the reverse;
- refitting any frozen instrument;
- any claim about Junevity's published numbers being unreproducible.

--- STOP KEYS ---
prior_results_not_reproduced; axis_not_readable; space_mismatch.
Seeds: 20260914 splits and permutations, 20260918 bootstraps.
n_boot = 200.

--- AMENDMENTS ADDED BEFORE ANY AXANG NUMBER EXISTS ---

These are part of the pre-registration and are written to PREREG.flag with
the rest.

E1. A fourth aging axis, for the preprint's own data.
 A4, Lu donor pair: the direction from the 96-year-old donor's (GM00731)
     day-0 fibroblasts to the 22-year-old donor's (GM23815) day-0
     fibroblasts in GSE297234, as the preprint defines it, read from stored
     results if present. Add A4 to every D1 table, with its C1 built from
     the day-0 fibroblasts of both donors and its C2 as defined. Report the
     stored preprint value cos(A4, GTEx age axis) = 0.039 in Stage 0 and
     reproduce it. A4 is reported, not gated; the D1 gate is unchanged.

E2. Additions to D2, report-only:
 - For each of the three datasets (Lu/GSE297234, Southard, Sengstack), the
   number of the 23,485 reference-panel genes present and absent, and the
   share of the GTEx young-minus-old centroid difference's squared length
   that lies in the genes that dataset does not measure, as a percentage.
   This bounds how much of any approach is unreachable by construction and
   is reported, not corrected.
 - Whether any earlier SENG task computed a GTEx-anchored approach test on
   the Sengstack data. Report what it found if so; state plainly if not.
 - The distance from the Sengstack PD32 pseudobulk to the GTEx young and old
   centroids, against the young-to-old gap, in the same gene space,
   mirroring Supplementary Note 1 of the preprint. A bound, not a test.

E3. D3 — Junevity's four factors against donor aging (gated).
 This section is exempt from the line under WHAT DOES NOT COUNT that forbids
 naming factors or reporting per-factor statistics. That line still applies
 to D1 and D2.

 Method: the Southard GTEx-anchored test exactly as SOUTH8/SOUTH9 run it —
 each perturbation's displacement (perturbed minus its own library's NT
 cells, same machine) applied to the GTEx old centroid, asking whether it
 ends closer to the GTEx young centroid, with a floor built at each
 perturbation's own cell count from NT draws in its own library, in
 S2A_k20 primary. Read the recipe off the committed SOUTH8/SOUTH9 code and
 quote the lines. If it cannot be read, STOP `south_recipe_unreadable`.
 State that this transplants a displacement measured on one platform onto a
 bulk GTEx origin, wherever a number appears.

 D3a, plant gate first, before any factor number. Displace NT draws along
 the GTEx young-minus-old direction by 25, 50 and 100% of the gap, 200 draws
 per dose, at 96, 107, 128 and 157 cells, both libraries. p25 and p50 are the
 lowest powers across the eight combinations.
  - QUALIFIED (p25 >= 0.80): "the GTEx-anchored test detects a 25% approach
    toward donor youth at these factors' depths, lowest power <p25>."
  - LIMITED (p25 < 0.80, p50 >= 0.80): "it detects a 50% approach but not a
    25% one, lowest powers <p50> and <p25>; every D3b number speaks only to
    approaches of at least half the donor-age gap."
  - BLIND (p50 < 0.80): "the GTEx-anchored test cannot see a half-gap
    approach toward donor youth at these depths, lowest power <p50>; D3b is
    reported and carries no conclusion about any factor."

 D3b, for all 340 perturbations and for the four individually:
  - called approaching (below floor and distance decreased), against the
    chance expectation of 17.0 of 340 with its 95% binomial interval, per
    modality;
  - for each of the four: cell count, change in distance, floor, called,
    and the cosine and angle of its displacement with the GTEx age axis,
    with the passage axis A3, and with its library's C1 proliferation axis;
  - SOUTH8's stored result for the same four gene names in the Southard
    screen beside them, labelled as activation-only and, for STAT3 and ZFX,
    the opposite perturbation.

 D3 gate, on the four, applied only if D3a did not fire BLIND:
  - APPROACH (two or more called): "<k> of Junevity's four factors move
    their own cells toward young human donors as measured against GTEx:
    <list>."
  - NO APPROACH (fewer than two called): "none of Junevity's four factors
    measurably moves their cells toward young human donors as measured
    against GTEx, with an instrument qualified at their depth."

 No factor named in D3 is validated or a pilot candidate. Passaged
 neonatal-line cells remain the substrate whatever fires.

E4. Additional STOP key: south_recipe_unreadable.
```

## 2. Amendments and deviations

One self-contradiction in the flag was put to the author before any
conclusion was written (F10); nine further definitions or citations the flag
omits or gets wrong were fixed in writing before computing. Full text of
`results/axang/AMENDMENTS.md`:

```
# AXANG — amendments to the pre-registration

Written 2026-10-06T23:50Z, **before any AXANG number existed**. The flag at
`results/axang/PREREG.flag` (2026-10-06T23:34:40Z) is **not edited**. None of
these needed to go to the author: none is a contradiction that decides a gate.
Each is a definition the flag does not supply, or a citation it gets wrong, and
each was fixed in writing here before computing.

## F1 — The cell-cycle score for bulk GTEx donors (missing definition)

C1 says to "score the units of that dataset for cell-cycle activity with the
Tirosh S and G2M lists", and for cells this project's recipe is
`module_score`: mean `log1p(CP10k)` over the listed genes
(`src/seng_run.py:915`). That function is defined on raw counts and a library
size, and **GTEx units are bulk donors with no counts in the frozen pack** —
`load_gtex_z` exposes log2-CPM, not counts.

**Substituted, and disclosed as a substitution:** for GTEx donors the score is
the **mean log2-CPM over the same Tirosh genes**, the closest available bulk
analogue of a mean expression level over a gene set. Quartiles are then taken
exactly as `quartile_masks` does (`<= q25`, `>= q75`). For every single-cell
dataset (Southard, Sengstack, Lu) the score is this project's own
`module_score` on counts, unchanged.

A consequence that matters and is reported: the GTEx score is a *mean z-free
expression level* while the cell scores are *mean log1p(CP10k)*. The two are
monotone in the same quantity but are not the same statistic, so
`C1_gtex` is not a like-for-like estimate of `C1_sengstack`. Angles involving
`C1_gtex` are labelled as bulk-donor measurements throughout, as the flag's
"what does not count" requires.

## F2 — Which GTEx donors are the units for `C1_gtex` (missing definition)

The flag says "GTEx donors" without naming a set. **Primary: all 652 donors of
the frozen reference panel**, the largest set that is not selected on age.
Reported beside it as a disclosed secondary: the same axis built on only the
344 anchor donors (113 young + 231 old) that define A1.

**Item (d) of the method statement bites here and is not refitted away.** The
donors that define A1 are a subset of the donors that define `C1_gtex` under
either choice, so `cos(A1, C1_gtex)` — the gate quantity — is estimated partly
from the same units on both sides. This is stated again at the gate.

## F3 — The "triangle" has four vectors, not three (flag's own anticipation)

The flag asks for "the angle between A1 and A3, the angle between A3 and its
C1, and the angle between A1 and its C1" as a triangle consistency check. Those
involve **four** distinct vectors — A1, A3, `C1_gtex` and `C1_sengstack` — so
they do not form a triangle and no consistency can be asserted between them.

**Fixed in writing:** two genuine spherical triangles are reported instead,
each on three vectors, with the standard inequality

    | theta(X,Z) - theta(Y,Z) |  <=  theta(X,Y)  <=  theta(X,Z) + theta(Y,Z)

checked explicitly:
- triangle T1: A1, A3, `C1_sengstack`
- triangle T2: A1, A3, `C1_gtex`

The flag's own four-vector version is also tabulated, labelled as not a
triangle. The flag says an inconsistency "is itself reportable and must be
stated rather than smoothed"; the inconsistency here is in the question's
construction, and this is the statement of it.

## F4 — D3a has no zero dose, but its floor needs one (missing definition)

D3a lists doses 25, 50 and 100% and no 0%. The criterion it gates on ("below
floor and distance decreased") requires a floor, and D3's method paragraph says
the floor is built "at each perturbation's own cell count from NT draws in its
own library". **Fixed in writing:** a 0% dose set of 200 effect-free NT draws is
computed at each cell count and library, its 5th percentile is the floor, and
the 0% row is reported in the table as the noise floor rather than subtracted.
This is the same construction SENG3 and SENG4 used.

## F5 — The Southard C1 is stored and is read, not rebuilt

`results/south6/stage1.npz` key `p_unit` is already the Southard NT
bottom-to-top-quartile cell-cycle direction, built from the same Tirosh lists
(`src/south6_run.py:1193-1241`, note "built from NT quartiles n_bottom=3504
n_top=3504"). CLAUDE.md requires frozen instruments be read rather than refit,
so `C1_southard` is read from that file. 14,013 NT cells, 3,504 per quartile.

## F6 — A4 and its stored cosine

A4 is read frozen from `results/same/panel_z.npz` key `v`, which SAME defines
as `z(Y) - z(O)` with `Y = GM23815 d0 Fibroblast, half A` (22-year-old donor)
and `O` the GM00731 (96-year-old) day-0 fibroblast pseudobulk
(`results/same/PREREG.flag` line 16). The stored cosine against the GTEx age
direction is `cos_v_u_toward` in the same file. The flag quotes it as 0.039; the
stored value is 0.038815304035390776, which rounds to it. Reproduced in
Stage 0 against the stored value, with the rounding noted.

**A4 is half-A only.** SAME split the donor's day-0 fibroblasts and used half A
for the target. A4 therefore carries SAME's split, not all day-0 cells.
`C1_lu`, by contrast, is built from **all** day-0 fibroblasts of both donors as
E1 specifies, so A4 and `C1_lu` are not estimated from identical cell sets but
do overlap. Disclosed, not refitted.

## F7 — Spaces, and one place a cosine is not comparable

Every axis in this task lives natively in the 23,485-column frozen-ruler z
space, so `space_mismatch` does not fire and all pairs are computable in both
the full-gene space and S2A_k20. Cosines are reported per space and never
compared across spaces. The S2A_k20 basis is read frozen from
`results/south7/basis.npz` (308 GTEx donors aged 40-59), is fit on no cell of
any of the three single-cell datasets, and adds no circularity.

One exception is flagged where it appears: the Southard axes live natively on
20,264 matched ruler genes and the Sengstack axes on 16,781, both scattered into
the 23,485 frame with the unmeasured coordinates at zero. A cosine between a
Southard-derived and a Sengstack-derived axis is therefore taken over a frame
in which each vector is structurally zero on genes the other dataset measures.
That is reported beside the number, not corrected.

## F8 — E2's "share of squared length" is computed on the GTEx age axis

E2 asks for "the share of the GTEx young-minus-old centroid difference's squared
length that lies in the genes that dataset does not measure". Computed as
`sum(d_i^2 for i not measured) / sum(d_i^2)` where `d = c_young - c_old` read
frozen from `results/toward/anchors.npz`, per dataset's measured-gene set.

## F9 — Bootstrap scope

The flag asks for 200-draw bootstraps "resampling the unit that defines each
axis (donors for A1 and A2, cells for C1 in each dataset)", seed 20260918.
Resampled: A1 (its 113 young and 231 old donors, separately), `C1_gtex` (652 or
344 donors), `C1_sengstack` (the NT half-one cells of each library),
`C1_lu` (day-0 cells of both donors). **Not resampled, and stated as such:**
A2 (the frozen ruler is a fitted instrument with no stored per-donor
resampling pack), A3 (SENG4's definition, whose per-cell resampling is not
part of this task), A4 (read frozen from a stored vector) and `C1_southard`
(read frozen from a stored vector). Those four have point estimates only, and
any gate or claim resting on them is labelled as having no interval. The gate
quantity `cos(A1, C1_gtex)` does have an interval.

## F10 — The D1 gate's threshold variable (self-contradiction; author's answer)

**The contradiction.** The gate is written on "the angle between A1 and its
data-derived proliferation axis C1", SEPARABLE at `d >= 45`. The measured value
is `cos(A1, C1_gtex) = -0.8108`: a **signed** angle of 144.17 degrees, which
fires SEPARABLE and asserts donor-age rejuvenation is distinguishable from
proliferation — while the same row shows projecting that axis out removes
**41.5%** of the donor-age axis's length. Both cannot be reported as true.
Anti-alignment is not separability.

**Resolution, as answered: the magnitude reading.** `d = arccos|cos|`, the angle
from parallel, which is 35.83 degrees. `25 <= d < 45` fires **INTERMEDIATE**.
Three things support it and are recorded rather than assumed: the INTERMEDIATE
branch's own placeholder "proliferation still carries `<pct>`% of it" is exactly
the projection loss and is meaningless under the signed reading; SENG4's
equivalent gate was written on `|cos(P1, a)|`, sign-invariant; and separability
is a property of shared direction, not of sign.

**The sign is not discarded.** That `cos(A1, C1_gtex)` is negative — GTEx donors
with higher cell-cycle scores sitting on the *old* side of the donor-age axis —
is reported as a finding in its own right, and as something the flag's recorded
geometric expectation did not anticipate.

The same magnitude rule is applied to C2 for the comparison the flag asks for.
```

## 3. State the method before using it

### (a) and (b) Every axis, its source, and whether it is read or rebuilt

| axis | what it is | source | read frozen or rebuilt |
|---|---|---|---|
| **A1** | GTEx young-minus-old centroid direction, young 20-39 (n=113), old 60-79 (n=231) | `results/toward/anchors.npz` key `u` | **read frozen** |
| **A2** | unit-length frozen fibroblast ruler weight vector | `results/fibro/frozen_ruler_ridge_raw.npz` key `w` | **read frozen** |
| **A3** | Sengstack PD14 half A minus PD32 half A | SENG3's frozen panel, `results/seng2/panel_counts.npz` | rebuilt by SENG3's own definition |
| **A4** | Lu GM00731 (96y) day-0 to GM23815 (22y) day-0, half A | `results/same/panel_z.npz` key `v` | **read frozen** |
| **C1 GTEx** | bottom-to-top cell-cycle quartile over GTEx donors | built here (F1, F2) | rebuilt |
| **C1 Southard** | bottom-to-top NT cell-cycle quartile, 3,504 per quartile | `results/south6/stage1.npz` key `p_unit` | **read frozen** (F5) |
| **C1 Sengstack** | bottom-to-top NT half-one quartile, per modality | `results/seng4/prep.npz` | rebuilt by SENG4's definition |
| **C1 Lu** | bottom-to-top quartile over day-0 fibroblasts of both donors | built here | rebuilt |
| **C2** | equal-weight Tirosh direction, no cells | built here | rebuilt, deterministic |

The C1 recipe, quoted from committed code — the per-unit score
(`src/seng_run.py:915`) and the quartile rule (`src/seng_run.py:925`):

```python
def module_score(counts_csr, cols, lib):
    sub = counts_csr[:, cols]; scale = np.maximum(lib, 1.0)[:, None]
    return np.log1p(arr / scale * 1e4).mean(axis=1)

q25, q75 = np.quantile(scores[finite], [0.25, 0.75])
bottom = finite & (scores <= q25);  top = finite & (scores >= q75)
```

and the Southard axis, `src/south6_run.py:1240`:

```python
p_unit = pvec / nrm
prolif_note = f"built from NT quartiles n_bottom={int(bottom.sum())} n_top={int(top.sum())} ||p||={nrm:.4f}"
```

Units entering each C1: GTEx **652 donors** (163 per quartile; anchors-only variant 344 donors, 86 per quartile); Southard **14,013 NT cells** (3,504 per quartile); Sengstack **333 CRA / 739 CRI** NT half-one cells; Lu **7,782 + 5,021** day-0 fibroblasts. C2 uses 94 of 97 Tirosh genes present on the ruler (missing ['MLF1IP', 'FAM64A', 'HN1']).

### (c) Spaces

Every axis lives natively in the 23,485-column frozen-ruler z space, so
**`space_mismatch` does not fire** and every pair is computable in both the
full-gene space and S2A_k20. Cosines are reported per space and are **never**
compared across spaces: a 20-dimensional cosine and a 23,485-dimensional one
are different quantities, and `results/axang/d1_angles.csv` carries a `space`
column on every row. The S2A_k20 basis is read frozen from
`results/south7/basis.npz` (308 GTEx donors aged 40-59) and is fit on no cell
of any single-cell dataset here.

One caveat is stated where it appears (F7): Southard axes are structurally
zero on 3,221 ruler genes and Sengstack/Lu axes on 6,704, so a cosine between
a Southard-derived and a Sengstack-derived axis is taken over a frame where
each is zero on genes the other measures.

### (d) Axes estimated from units that also define another axis

**Yes, in two places, and neither is refitted away.**

1. **The gate quantity itself.** The 344 donors defining A1 are a subset of
   the 652 donors defining `C1_gtex`, so `cos(A1, C1_gtex)`
   is estimated partly from the same units on both sides. The anchors-only
   variant makes the overlap total and is reported beside it; the two agree
   (-0.8108 against -0.7646),
   so the overlap is not what produces the value — but the number is labelled.
2. **A4 and `C1_lu`** both come from the Lu day-0 fibroblasts, A4 from SAME's
   half A only and `C1_lu` from all of them (F6). They overlap but are not
   identical sets.

A3 against `C1_seng` is **not** in this list: SENG4 reserved NT half one for
the axis and half two for everything else, and A3 is built from PD14/PD32
cells which enter no C1.

## 4. Stage 0 — reproduce first (gate)

**Gate: PASS.** All 14 checks.

| item | kind | stored | re-derived | ok |
|---|---|---|---|---|
| `SENG4 cos(P1,a) CRA` | numeric | 0.941708 | 0.941708 | yes |
| `SENG4 cos(P1,a) CRI` | numeric | 0.936287 | 0.936287 | yes |
| `SENG4 cos(P3,a)` | numeric | 0.474385 | 0.474385 | yes |
| `SENG4 cos(a,u) full-gene` | numeric | -0.040273 | -0.040273 | yes |
| `SENG4 cos(a,u) S2A_k20` | numeric | -0.199293 | -0.199293 | yes |
| `SOUTH9 res power_25 @106` | numeric | 0.990000 | 0.990000 | yes |
| `SOUTH9 raw power_25 @106` | numeric | 0.980000 | 0.980000 | yes |
| `SOUTH8 all_1836 m` | count | 1836 | 1836 | yes |
| `SOUTH8 all_1836 amended_call_S2A_k20` | count | 125 | 125 | yes |
| `SOUTH8 all_1836 below_depth_floor_S2A_k20` | count | 153 | 153 | yes |
| `E1 stored cos(A4, GTEx age axis)` | numeric | 0.039000 | 0.039000 | yes |
| `E1 cos(A4,u) re-derived vs stored` | numeric | 0.038815 | 0.038815 | yes |
| `anchors n_young` | count | 113 | 113 | yes |
| `anchors n_old` | count | 231 | 231 | yes |

Two notes. The frozen anchors are **exactly** the pack centroids —
max |c_young − mean(Z[young])| = 0.0 and
the same for old — so A1 is the axis this project has always used. And the
E1 stored value is 0.038815304035391, which the flag quotes rounded
to 0.039; both the rounded and the exact forms are checked.

## 5. D1 — the angles (gated)

### Each aging axis against every proliferation construction, S2A_k20

Angle is the signed angle; **angle from parallel** is `arccos|cos|`, the
separability quantity (F10). `% removed` is `1 − sqrt(1 − cos²)`.

| aging axis | proliferation axis | cos | 2.5–97.5% | signed angle | from parallel | % removed |
|---|---|---|---|---|---|---|
| A1 GTEx donor age | C1 GTEx donors (652, bulk) | -0.8108 | -0.939 to -0.340 | 144.17° | **35.83°** | 41.5% |
| A1 GTEx donor age | C1 GTEx anchor donors (344, bulk) | -0.7646 | -0.930 to -0.286 | 139.87° | **40.13°** | 35.6% |
| A1 GTEx donor age | C1 Southard NT cells (14,013) | -0.3384 | -0.633 to +0.119 | 109.78° | **70.22°** | 5.9% |
| A1 GTEx donor age | C1 Sengstack NT CRA (333) | -0.2891 | -0.564 to +0.153 | 106.81° | **73.19°** | 4.3% |
| A1 GTEx donor age | C1 Sengstack NT CRI (739) | -0.2689 | -0.563 to +0.168 | 105.60° | **74.40°** | 3.7% |
| A1 GTEx donor age | C1 Lu day-0 (12,803) | -0.3264 | -0.634 to +0.136 | 109.05° | **70.95°** | 5.5% |
| A1 GTEx donor age | C2 Tirosh gene list (94 genes) | -0.6588 | -0.753 to -0.352 | 131.21° | **48.79°** | 24.8% |
| A2 frozen ruler | C1 GTEx donors (652, bulk) | +0.0782 | +0.054 to +0.097 | 85.51° | **85.51°** | 0.3% |
| A2 frozen ruler | C1 GTEx anchor donors (344, bulk) | +0.0787 | — | 85.49° | **85.49°** | 0.3% |
| A2 frozen ruler | C1 Southard NT cells (14,013) | -0.1653 | — | 99.51° | **80.49°** | 1.4% |
| A2 frozen ruler | C1 Sengstack NT CRA (333) | -0.2345 | -0.267 to -0.202 | 103.56° | **76.44°** | 2.8% |
| A2 frozen ruler | C1 Sengstack NT CRI (739) | -0.2071 | -0.228 to -0.185 | 101.95° | **78.05°** | 2.2% |
| A2 frozen ruler | C1 Lu day-0 (12,803) | -0.1395 | — | 98.02° | **81.98°** | 1.0% |
| A2 frozen ruler | C2 Tirosh gene list (94 genes) | -0.0942 | — | 95.40° | **84.60°** | 0.4% |
| A3 Sengstack passage | C1 GTEx donors (652, bulk) | +0.5681 | +0.520 to +0.615 | 55.38° | **55.38°** | 17.7% |
| A3 Sengstack passage | C1 GTEx anchor donors (344, bulk) | +0.6133 | — | 52.17° | **52.17°** | 21.0% |
| A3 Sengstack passage | C1 Southard NT cells (14,013) | +0.8886 | — | 27.30° | **27.30°** | 54.1% |
| A3 Sengstack passage | C1 Sengstack NT CRA (333) | +0.9417 | +0.916 to +0.953 | 19.66° | **19.66°** | 66.4% |
| A3 Sengstack passage | C1 Sengstack NT CRI (739) | +0.9363 | +0.915 to +0.951 | 20.56° | **20.56°** | 64.9% |
| A3 Sengstack passage | C1 Lu day-0 (12,803) | +0.8421 | — | 32.63° | **32.63°** | 46.1% |
| A3 Sengstack passage | C2 Tirosh gene list (94 genes) | +0.4744 | — | 61.68° | **61.68°** | 12.0% |
| A4 Lu donor pair | C1 GTEx donors (652, bulk) | +0.1101 | +0.082 to +0.138 | 83.68° | **83.68°** | 0.6% |
| A4 Lu donor pair | C1 GTEx anchor donors (344, bulk) | +0.0989 | — | 84.33° | **84.33°** | 0.5% |
| A4 Lu donor pair | C1 Southard NT cells (14,013) | +0.4308 | — | 64.48° | **64.48°** | 9.8% |
| A4 Lu donor pair | C1 Sengstack NT CRA (333) | +0.4163 | +0.353 to +0.478 | 65.40° | **65.40°** | 9.1% |
| A4 Lu donor pair | C1 Sengstack NT CRI (739) | +0.3118 | +0.284 to +0.349 | 71.83° | **71.83°** | 5.0% |
| A4 Lu donor pair | C1 Lu day-0 (12,803) | +0.5355 | — | 57.63° | **57.63°** | 15.5% |
| A4 Lu donor pair | C2 Tirosh gene list (94 genes) | +0.2239 | — | 77.06° | **77.06°** | 2.5% |

Rows with no interval are the four axes F9 could not resample because they
are read from stored vectors with no per-unit pack (A2, A3, A4, `C1_south`).
Every bootstrapped interval contained its point estimate, so none is flagged
invalid.

### Aging axis against aging axis, S2A_k20

| pair | cos | signed angle |
|---|---|---|
| A1 (GTEx donor age) vs A2 (frozen ruler) | -0.2982 | 107.35° |
| A1 (GTEx donor age) vs A3 (Sengstack passage) | -0.1993 | 101.50° |
| A1 (GTEx donor age) vs A4 (Lu donor pair) | -0.0014 | 90.08° |
| A2 (frozen ruler) vs A3 (Sengstack passage) | -0.1872 | 100.79° |
| A2 (frozen ruler) vs A4 (Lu donor pair) | -0.1092 | 96.27° |
| A3 (Sengstack passage) vs A4 (Lu donor pair) | +0.3432 | 69.93° |

The four aging models are mutually near-orthogonal, which is this project's
standing finding extended to four axes: A1 vs A3 is
-0.1993 and A1 vs A4 is -0.0014
— the donor-age axis is essentially orthogonal to both the passage model and
the preprint's own donor-pair model.

### Proliferation axis against proliferation axis, S2A_k20

| pair | cos | signed angle |
|---|---|---|
| C1 GTEx donors (652, bulk) vs C1 GTEx anchor donors (344, bulk) | +0.9944 | 6.07° |
| C1 GTEx donors (652, bulk) vs C1 Southard NT cells (14,013) | +0.6423 | 50.03° |
| C1 GTEx donors (652, bulk) vs C1 Sengstack NT CRA (333) | +0.6008 | 53.07° |
| C1 GTEx donors (652, bulk) vs C1 Sengstack NT CRI (739) | +0.6076 | 52.58° |
| C1 GTEx donors (652, bulk) vs C1 Lu day-0 (12,803) | +0.6346 | 50.61° |
| C1 GTEx donors (652, bulk) vs C2 Tirosh gene list (94 genes) | +0.7198 | 43.96° |
| C1 GTEx anchor donors (344, bulk) vs C1 Southard NT cells (14,013) | +0.6774 | 47.36° |
| C1 GTEx anchor donors (344, bulk) vs C1 Sengstack NT CRA (333) | +0.6368 | 50.45° |
| C1 GTEx anchor donors (344, bulk) vs C1 Sengstack NT CRI (739) | +0.6509 | 49.39° |
| C1 GTEx anchor donors (344, bulk) vs C1 Lu day-0 (12,803) | +0.6614 | 48.59° |
| C1 GTEx anchor donors (344, bulk) vs C2 Tirosh gene list (94 genes) | +0.7193 | 44.00° |
| C1 Southard NT cells (14,013) vs C1 Sengstack NT CRA (333) | +0.9551 | 17.23° |
| C1 Southard NT cells (14,013) vs C1 Sengstack NT CRI (739) | +0.9066 | 24.97° |
| C1 Southard NT cells (14,013) vs C1 Lu day-0 (12,803) | +0.9478 | 18.60° |
| C1 Southard NT cells (14,013) vs C2 Tirosh gene list (94 genes) | +0.6921 | 46.20° |
| C1 Sengstack NT CRA (333) vs C1 Sengstack NT CRI (739) | +0.9784 | 11.94° |
| C1 Sengstack NT CRA (333) vs C1 Lu day-0 (12,803) | +0.8814 | 28.19° |
| C1 Sengstack NT CRA (333) vs C2 Tirosh gene list (94 genes) | +0.6183 | 51.81° |
| C1 Sengstack NT CRI (739) vs C1 Lu day-0 (12,803) | +0.8034 | 36.54° |
| C1 Sengstack NT CRI (739) vs C2 Tirosh gene list (94 genes) | +0.5513 | 56.55° |
| C1 Lu day-0 (12,803) vs C2 Tirosh gene list (94 genes) | +0.6907 | 46.32° |

### The full-gene space, the same aging-against-proliferation pairs

Reported because the flag asks for both spaces. **Not comparable to the table
above**: these are 23,485-dimensional cosines.

| aging axis | proliferation axis | cos | from parallel | % removed |
|---|---|---|---|---|
| A1 | C1 GTEx donors (652, bulk) | -0.7065 | 45.05° | 29.2% |
| A1 | C1 Southard NT cells (14,013) | -0.1298 | 82.54° | 0.8% |
| A1 | C1 Sengstack NT CRA (333) | -0.0760 | 85.64° | 0.3% |
| A1 | C1 Lu day-0 (12,803) | -0.1287 | 82.61° | 0.8% |
| A1 | C2 Tirosh gene list (94 genes) | -0.0873 | 84.99° | 0.4% |
| A2 | C1 GTEx donors (652, bulk) | +0.0073 | 89.58° | 0.0% |
| A2 | C1 Southard NT cells (14,013) | -0.0213 | 88.78° | 0.0% |
| A2 | C1 Sengstack NT CRA (333) | -0.0148 | 89.15° | 0.0% |
| A2 | C1 Lu day-0 (12,803) | -0.0225 | 88.71° | 0.0% |
| A2 | C2 Tirosh gene list (94 genes) | -0.0159 | 89.09° | 0.0% |
| A3 | C1 GTEx donors (652, bulk) | +0.2582 | 75.04° | 3.4% |
| A3 | C1 Southard NT cells (14,013) | +0.4069 | 65.99° | 8.7% |
| A3 | C1 Sengstack NT CRA (333) | +0.5035 | 59.77° | 13.6% |
| A3 | C1 Lu day-0 (12,803) | +0.3585 | 68.99° | 6.6% |
| A3 | C2 Tirosh gene list (94 genes) | +0.0507 | 87.10° | 0.1% |
| A4 | C1 GTEx donors (652, bulk) | +0.0319 | 88.17° | 0.1% |
| A4 | C1 Southard NT cells (14,013) | +0.1711 | 80.15° | 1.5% |
| A4 | C1 Sengstack NT CRA (333) | +0.1571 | 80.96° | 1.2% |
| A4 | C1 Lu day-0 (12,803) | +0.2612 | 74.86° | 3.5% |
| A4 | C2 Tirosh gene list (94 genes) | +0.0275 | 88.43° | 0.0% |

### The triangle consistency check

The flag's version names four vectors, not three, so it is not a triangle
(F3). Two genuine spherical triangles are checked instead, on the standard
inequality:

| triangle | θ(A1,A3) | θ(A1,third) | θ(A3,third) | allowed range for θ(A1,A3) | consistent |
|---|---|---|---|---|---|
| T1_A1_A3_C1seng (third = C1 Sengstack NT CRA (333)) | 101.50° | 106.81° | 19.66° | 87.15° to 126.47° | **yes** |
| T2_A1_A3_C1gtex (third = C1 GTEx donors (652, bulk)) | 101.50° | 144.17° | 55.38° | 88.79° to 199.56° | **yes** |

**Both are consistent**, so there is no geometric contradiction in the
measurements. The flag's expectation failed not because the geometry is
inconsistent but because it assumed one proliferation direction where there
are several, and they differ: the bulk-donor construction is far from the
cell-derived ones (see the proliferation-against-proliferation table).

### The D1 gate

**INTERMEDIATE** fired, on `d` = 35.83° — the angle
from parallel between A1 and `C1_gtex` in S2A_k20, per amendment F10
(`cos` = -0.810808, signed angle 144.17°).

> "donor aging sits 35.83 degrees from
> proliferation, against 19.66 for passage aging. Donor age is
> the better separated of the two models but proliferation still carries
> 41.5% of it, and any donor-age rejuvenation readout must
> report a proliferation-residualized version beside the headline."

**The ENTANGLED branch did not fire**, so the clause about SOUTH9's RETAINED
reading and SOUTH8's residualized count being placed in question does not
apply. Those files are untouched and unquestioned by this task.

**C2 would fire a different branch, and the flag requires that be said
prominently.**

| construction | cos with A1 | from parallel | branch |
|---|---|---|---|
| **C1, data-derived (primary)** | -0.8108 | 35.83° | **INTERMEDIATE** |
| C2, gene-list | -0.6588 | 48.79° | SEPARABLE |

C2's sentence, reported but **not** the task's reading: "donor aging sits
48.79 degrees from an independently measured
proliferation axis, against 19.66 degrees for passage aging in the
same construction. A rejuvenation signal can be distinguished from a change
in proliferative state against donor age and cannot be against passage
number." C1 is primary because the flag fixed it so, and because it is the
construction SENG4's primary axis used.

## 6. D2 — report-only, no decision rule

### Do the Southard- and Sengstack-derived C1 axes agree?

`cos(C1_south, C1_seng_CRA)` = +0.9551 (17.23°)
and `cos(C1_south, C1_seng_CRI)` = +0.9066 (24.97°), in S2A_k20.

The same recipe on two different CRISPR screens in two different fibroblast
systems produces axes that agree substantially but not closely. That is
evidence the recipe measures a shared, real direction rather than noise, and
also that it is not measuring exactly the same thing in both — part of the
difference is the structural-zero issue of F7, since each axis is zero on
genes the other dataset measures.

### Dependence on the construction: C2 against C1

| aging axis | % removed by its C1 | % removed by C2 |
|---|---|---|
| A1 GTEx donor age | 41.5% | 24.8% |
| A2 frozen ruler | 0.3% | 0.4% |
| A3 Sengstack passage | 66.4% | 12.0% |
| A4 Lu donor pair | 15.5% | 2.5% |

The choice of construction changes the answer by tens of percentage points,
which is why the flag fixed C1 as primary before any number.

### Reference-panel coverage, and how much of the donor-age gap is unreachable

| dataset | panel genes present | absent | % of the GTEx age gap's squared length in unmeasured genes |
|---|---|---|---|
| Lu | 16,781 | 6,704 | **21.97%** |
| Southard | 20,264 | 3,221 | **10.00%** |
| Sengstack | 16,781 | 6,704 | **21.97%** |

So between a tenth and a fifth of what defines the GTEx donor-age direction
is in genes these datasets do not measure at all. **This bounds how much of
any approach is unreachable by construction and is reported, not corrected.**

### Did any earlier SENG task run a GTEx-anchored test?

**A direction test, yes; an approach test, no.** `results/seng/screen.csv` carries
columns `cos_d_u, p_cos_d_u`: SENG computed the cosine of each
displacement with the GTEx donor-age direction and a permutation p for all
340 perturbations, and **55 of 340** reached p ≤ 0.05
(median cos -0.0125, max +0.0853). It did
**not** run an approach or distance test anchored on the GTEx centroids. D3
below is the first GTEx-anchored *approach* test on the Sengstack data.

### The Sengstack cells against the GTEx centroids

Mirroring Supplementary Note 1 of the preprint, in the full-gene space:

- distance from the Sengstack PD32 pseudobulk to the GTEx **young** centroid:
  **465.0**
- to the GTEx **old** centroid: **461.1**
- the GTEx young-to-old gap itself: **32.28**
- so PD32 sits **14.4×** the entire donor-age gap away
  from young and **14.3×** away from old — essentially
  equidistant, and far outside the span the gap defines.

**A bound, not a test.** Any displacement of these cells toward GTEx young is
a move of at most a few percent of the distance to it, from a starting point
that is not between the two donor-age anchors in the first place.

## 7. D3 — the four factors against donor aging (gated)

**Every number in this section transplants a displacement measured on 10x
single cells onto a bulk GTEx origin.** The displacement is the perturbation
minus its own library's NT cells, on the same machine; the origin is the GTEx
old centroid and the target the GTEx young centroid, which are bulk
post-mortem tissue. That is the Southard construction
(`src/south6_run.py:506` `block_delta`, with `r = c_old − c_young` from
`src/south7_run.py:134`), and `south_recipe_unreadable` does not fire.

### D3a — the plant gate, run before any factor number

The GTEx gap is 32.28 in the full-gene space and
27.93 in S2A_k20; plant doses are fractions of it. The 0%
row is the noise floor and is reported, not subtracted (amendment F4).

| library | n cells | floor | power @25% | power @50% | power @100% |
|---|---|---|---|---|---|
| CRA | 96 | -1.030 | **0.945** | 1.000 | 1.000 |
| CRA | 107 | -1.223 | **0.955** | 1.000 | 1.000 |
| CRA | 128 | -1.346 | **0.945** | 1.000 | 1.000 |
| CRA | 157 | -1.343 | **0.990** | 1.000 | 1.000 |
| CRI | 96 | -1.024 | **0.945** | 1.000 | 1.000 |
| CRI | 107 | -1.052 | **0.970** | 1.000 | 1.000 |
| CRI | 128 | -1.079 | **0.965** | 1.000 | 1.000 |
| CRI | 157 | -1.202 | **0.980** | 1.000 | 1.000 |

**QUALIFIED** fired: p25 = 0.945, p50 = 1.000.

> "the GTEx-anchored test detects a 25% approach toward donor youth at
> these factors' depths, lowest power 0.945."

This is the opposite of the sensitivity problem SENG2 hit and SENG3 only
partly fixed: against the GTEx anchors the instrument is well powered at
every one of these cell counts.

### D3b — the counts

Chance expectation is 17.0 of 340.

| group | n | below floor | distance decreased | **called** | chance exp. | 95% interval |
|---|---|---|---|---|---|---|
| all | 340 | 5 | 17 | **5** | 17.0 | 0.5–3.4% |
| CRA | 167 | 3 | 10 | 3 | 8.3 | — |
| CRI | 173 | 2 | 7 | 2 | 8.7 | — |

**5 of 340 is below the 17 expected by
chance**, and only 17 of 340 even move closer at all. With an instrument
qualified at p25 = 0.945, that is a positive finding of absence rather than a
failure to detect.

### The four factors

| factor | n cells | pctl | Δdistance to GTEx young | floor | called | cos(d, GTEx age) | cos(d, A3 passage) | cos(d, its C1) |
|---|---|---|---|---|---|---|---|---|
| `CRA_E2F3` | 128 | 96% | **+11.448** | -1.346 | **no** | -0.1647 | +0.7655 | +0.8800 (28.4°) |
| `CRA_EZH2` | 107 | 89% | **+9.022** | -1.223 | **no** | -0.2669 | +0.7224 | +0.7863 (38.2°) |
| `CRI_STAT3` | 96 | 82% | **+5.327** | -1.024 | **no** | -0.1388 | +0.7782 | +0.8402 (32.8°) |
| `CRI_ZFX` | 157 | 99% | **+15.336** | -1.202 | **no** | -0.2051 | +0.7647 | +0.8969 (26.2°) |

**All four move away from young human donors, not toward them** — every delta
is positive against a floor near −1.2, and every `cos(d, GTEx age)` is
negative. The last two columns say where they do point: strongly along the
passage-rejuvenation direction (+0.72 to +0.78) and more strongly still along
their own proliferation axis (+0.79 to +0.90, 26-38°).

### SOUTH8's stored result for the same four gene names

Labelled as the flag requires: the Southard screen is **CRISPR activation
only**, so for `STAT3` and `ZFX` — which Junevity perturbs by
**interference** — this is the opposite perturbation and is not a replication.

| gene | in Southard | n cells | delta, S2A_k20 | below depth floor | amended call |
|---|---|---|---|---|---|
| `E2F3` *(activation, same direction)* | yes | 96 | +1.5621 | no | no |
| `EZH2` *(activation, same direction)* | yes | 140 | +0.4811 | no | no |
| `STAT3` *(opposite perturbation)* | yes | 83 | +2.3197 | no | no |
| `ZFX` *(opposite perturbation)* | yes | 129 | +2.4415 | no | no |

All four are positive-delta and uncalled in Southard too, on a different
cell line with a different screen.

### The D3 gate

**NO APPROACH** fired, on 0 of the four called.

> "none of Junevity's four factors measurably moves their cells toward
> young human donors as measured against GTEx, with an instrument
> qualified at their depth."

**No factor named here is validated or a pilot candidate.** Passaged
neonatal-line cells remain the substrate, the origin is bulk post-mortem
tissue, and between 10% and 22% of the donor-age direction is in genes this
dataset does not measure.

## 8. The reading

Only the pre-registered readings whose conditions fired are conclusions:

1. **Stage 0: PASS.**
2. **D1: INTERMEDIATE** on C1, `d` = 35.83° from
   parallel; C2 would fire SEPARABLE and is reported beside it.
3. **D3a: QUALIFIED**, p25 = 0.945.
4. **D3b: NO APPROACH**, 0 of four.

Put together, on one scale, in S2A_k20, as angle from parallel to a
proliferation axis built from the same units:

| aging model | angle | what that licenses |
|---|---|---|
| Sengstack passage (A3) | **19.66°** | nothing: SENG4 showed residualization here removes the signal |
| GTEx donor age (A1) | **35.83°** | a readout, provided a residualized version is reported beside it |
| Lu donor pair (A4) | **57.63°** | a readout; 15.5% of the axis is proliferation |
| frozen ruler (A2) | **85.51°** | a readout; 0.3% of the axis is proliferation |

And the factor-level result, which D3 is the first test in this project to
produce: with an instrument qualified at these depths, **none of the four
moves toward young human donors, 5 of 340 perturbations do against 17
expected by chance, and all four instead point along passage reversal and
proliferation.**

### What this does not cover

- **No angle measured in GTEx bulk donors is a statement about single cells,
  and no angle measured in cells is a statement about bulk donors.** The
  flag forbids both readings and the bulk/cell split is exactly where the
  expectation broke: A1's bulk-donor proliferation axis behaves unlike the
  three cell-derived ones.
- **What `C1_gtex` actually tracks is not established here.** In bulk
  post-mortem tissue a cell-cycle gene-set score can move with tissue
  quality, ischemic time or sampling as easily as with proliferation, and
  this task ran no diagnostic to separate those. Its strong anti-alignment
  with donor age should be read as a reason to doubt the construction, not
  as a biological claim. That is the main limitation of the D1 gate.
- **No claim about aged human cells from the Sengstack data**: passaged
  neonatal-line fibroblasts whose "old" is passage number, and
  cos(A1, A3) = -0.0403 in the full-gene
  space confirms again that passage is not donor aging.
- **A1, A2, A4 and `C1_south` have no intervals** (F9), so any statement
  resting on them is a point estimate only.
- **The gate quantity shares units on both sides** (method statement (d)).
- **No claim that Junevity's published numbers are unreproducible.**
- **D1 and D2 name no transcription factor**, as the flag requires; D3 is the
  exempted section.

### Files

- `results/axang/AMENDMENTS.md` (8,907 B)
- `results/axang/axes.npz` (5,109,661 B)
- `results/axang/d1.log` (7,507 B)
- `results/axang/d1_angles.csv` (19,596 B)
- `results/axang/d1_stderr.txt` (0 B)
- `results/axang/d1_stdout.txt` (7,507 B)
- `results/axang/d1_summary.json` (2,354 B)
- `results/axang/d2.log` (1,014 B)
- `results/axang/d2_summary.json` (1,117 B)
- `results/axang/d3.log` (1,046 B)
- `results/axang/d3_stderr.txt` (0 B)
- `results/axang/d3_stdout.txt` (1,046 B)
- `results/axang/d3_summary.json` (3,500 B)
- `results/axang/d3a_plant.csv` (1,859 B)
- `results/axang/d3b_table.csv` (70,190 B)
- `results/axang/PREREG.flag` (13,219 B)
- `results/axang/stage0.log` (2,804 B)
- `results/axang/stage0_reproduction.csv` (855 B)
- `results/axang/stage0_summary.json` (3,221 B)

Nothing was committed and nothing was pushed.
