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
