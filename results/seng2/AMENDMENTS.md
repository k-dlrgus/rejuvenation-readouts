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
