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
