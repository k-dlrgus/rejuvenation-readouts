# FINDINGS_SENG4 — does the proliferation adjustment survive an
# independent reference, and can proliferation be separated from passage
# aging at all?

Task SENG4. Computational, no new data. Reads `results/south6/`,
`results/south7/`, `results/south9/`, `results/seng/`, `results/seng2/`,
`results/seng3/` and the Sengstack h5ad already on disk. Writes
`results/seng4/` and this file. Nothing committed, nothing pushed. No
earlier findings file or results directory was written.

**The gates.** Stage 0 PASS; D1 **INSEPARABLE**; D2 **LOST**;
D3 **SPLIT** on the primary axis.

**The answer to the task's own question, stated once.** The prior question
decides it: in this dataset proliferation and passage aging are not
separable directions. |cos(P1, a)| = 0.941708, an angle of
19.66°, measured with an axis built on NT cells that
enter no origin, no target and no floor of the test. Projecting it out
removes 66.4% of the aging axis's length, and D2
confirms the consequence mechanically: power at a planted half-gap approach
falls from 0.940 to 0.010. **No residualization in this dataset can
distinguish a proliferation artefact from passage rejuvenation, so the
question SENG3 left open is not answerable here. That is the result.**

**A second finding, which cuts against this task's own premise.** The
circularity SENG3 disclosed was not what undermined its residualized
column. SENG3's circular axis and the independent P1 are cos
0.9996
apart (CRI; 0.9965 for CRA), and they give nearly the same screen counts —
27 called against SENG3's 25 under origin (i),
9 against 8 under origin (ii). Removing the
circularity changed almost nothing. Inseparability, not circularity, is what
makes the adjustment uninterpretable.

---

## 1. The pre-registration, verbatim

From `results/seng4/PREREG.flag`, written 2026-10-06T11:45:52Z before any
SENG4 number existed and not edited since.

```
PREREG written 2026-10-06T11:45:52Z (local 2026-10-06 20:45:52 KST)
Written before any SENG4 number exists. Reproduced verbatim at the top of
FINDINGS_SENG4.md. Not edited after this timestamp.

TASK: SENG4 — does the proliferation adjustment survive an independent
reference, and can proliferation be separated from passage aging at all?
Written before any SENG4 number exists.

Nature: computational. No new data. Reads results/south6/, results/south7/,
results/south9/, results/seng/, results/seng2/, results/seng3/ and the
Sengstack h5ad already on disk (data/raw/seng/exp3_merged_not_normalized.h5ad,
md5 bfccaa94a916c22e4b0793913a04128b). Writes only results/seng4/ and
FINDINGS_SENG4.md. No earlier findings file is regenerated, no earlier
results directory is written, nothing is committed or pushed.

--- THE QUESTION ---

SENG3 fired D4 APPROACH: all four named factors end closer to early-passage
cells under both origins, in the compressed space S2A_k20, with a calibrated
null and an instrument that detects a planted half-gap approach.

Its strongest qualification is unresolved. The approach statistic is
rank-correlated -0.8639 with the proliferation shift, and residualizing
proliferation drops the screen from 80 called to 25 and the four from 4 to 1
under origin (i) and to 0 under origin (ii). But SENG3's proliferation axis
is estimated from the same NT control cells that define the origin and every
floor. That circularity was disclosed and not removed, so the residualized
column is not a clean adjustment and the claim that rests on it — that what
these screens detect is proliferation — currently cannot be defended.

This task removes the circularity, and asks a prior question SENG3 did not:
whether proliferation and passage aging are separable directions in this
dataset at all. If they are not, no residualization can decide anything, and
that is itself the result.

--- STATE THE METHOD BEFORE USING IT ---

Report, before any SENG4 number:
 (a) exactly how SENG3 constructs its proliferation axis: what it is fit on,
     which cells, which genes, whether per-modality, and the line of code.
     Quote the source. If it cannot be read off committed code, STOP
     `prolif_recipe_unreadable`.
 (b) exactly how SENG3 residualizes: what is projected out of what, and
     whether the null is residualized by the same operation.
 (c) for each of the three axes defined in D1 below, what cells enter its
     fit, and confirm in writing that none of those cells enters the origin,
     the target or any floor of the test the axis is later applied to. If
     any does, say so plainly and label every number downstream as carrying
     that circularity. Do not refit to hide it.
 (d) the four factors' cell counts under SENG's assignment rule are
     unchanged: CRA_E2F3 128, CRA_EZH2 107, CRI_STAT3 96, CRI_ZFX 157.

--- STAGE 0: REPRODUCE FIRST (gate) ---

Re-derive with this task's code, from stored files:
 (a) SENG3 D1 firing rates in S2A_k20, origin (i): the maximum 0.0850 and
     the full eight-row table;
 (b) SENG3 D2 plant recovery, origin (i): p25 worst 0.565, p50 worst 0.940;
 (c) SENG3 D3 gate value r = +0.887005 and the per-library pair
     +0.811375 / +0.887005;
 (d) SENG3 D4 counts: 80 called origin (i), 73 origin (ii), of 340, and all
     four named factors called under both;
 (e) SENG3 D5: rho(proliferation shift, delta origin (i)) = -0.8639, and the
     residualized counts 80 -> 25 and 73 -> 8;
 (f) SOUTH9 D2 plant recovery after residualization at 106 cells in S2A_k20:
     power 0.990 residualized against 0.980 raw.
Pass: counts exact, numerics within 1e-4 * max(1, |stored|).
Fail -> STOP `prior_results_not_reproduced`, naming every mismatch.

--- D1. Three proliferation axes, and whether aging is separable from them
    (gated) ---

All three are built in S2A_k20, the SENG3 primary space, read frozen.

 P1, PRIMARY — disjoint NT halves. Split each library's NT pool in half by
 cell (seed 20260914), disjoint from the halves SENG3 used for its null
 validation if that is possible without reusing a cell; if it is not, say so
 and use a fresh split. Half one is reserved for this axis and for nothing
 else. Score its cells for cell-cycle activity using the Tirosh S and G2M
 lists already present at results/seng/seurat_cc.genes_tirosh2016.tsv, and
 define P1 as the direction from the bottom quartile to the top quartile of
 that score. Half two, and only half two, defines the origin and every floor
 in D2 and D3. Per modality.

 P2, SECONDARY — PD26. The 1,880 middle-passage cells enter no origin, no
 target and no floor anywhere in this project's Sengstack tests. Build the
 axis the same way inside PD26 alone: bottom-quartile to top-quartile
 cell-cycle score. One axis, not per modality, since PD26 has no modality.

 P3, SECONDARY — gene-set direction. The unit vector of the Tirosh S and G2M
 genes with equal weights, positive for cycling, projected into S2A_k20. No
 Sengstack cell enters this at all.

Report:
 - pairwise cosines between P1, P2, P3 and SENG3's circular axis, per
   modality where applicable;
 - the cosine between each axis and the aging axis a (PD14 half A minus PD32
   half A), and the fraction of ||a|| that survives projecting each axis out,
   sqrt(1 - cos^2);
 - the cosine between a and the GTEx donor-age direction u, re-derived, and
   the stored Southard equivalent if one exists in results/south6/ or
   results/south9/, for comparison;
 - how many cells enter each axis.

Gate, wording fixed now, on the PRIMARY axis P1 and its cosine with a, taken
as the larger absolute value across the two modalities. This is a reading,
not a STOP; all branches continue to D2.
 - SEPARABLE (|cos(P1, a)| < 0.40): "the proliferation direction and the
   passage-aging direction are substantially distinct in this dataset
   (cos <c>), so projecting one out leaves <pct>% of the other's length and
   an adjustment is meaningful."
 - ENTANGLED (0.40 <= |cos(P1, a)| < 0.70): "the proliferation direction
   carries a substantial part of the passage-aging direction (cos <c>);
   projecting it out removes <pct>% of the aging axis's length, and every
   residualized number below is reduced by that as well as by any confound."
 - INSEPARABLE (|cos(P1, a)| >= 0.70): "in this dataset passage aging and
   proliferation are substantially the same direction (cos <c>). Projecting
   proliferation out removes most of what defines young, so no
   residualization here can distinguish a proliferation artefact from
   passage rejuvenation. The question is not answerable in this dataset,
   and that is the finding."
If INSEPARABLE fires, D2 and D3 are still run to the end and reported, and
nothing in D3 is a conclusion about any factor. That is said in the gate
sentence and repeated wherever a factor is named.

Also report, whichever fires: the angle between a and each axis in degrees,
and state that this is the quantity the project handoff asked for as a
measure of how much of Junevity's signal is structurally confounded.

--- D2. Does residualization leave the statistic able to see an approach?
    (gated) ---

Runs before any residualized approach number and decides whether one can be
read.

Repeat SENG3's plant recovery with the residualization applied to plant,
null and floor alike, separately for P1, P2 and P3. Displace draws along a
by 0, 25, 50 and 100% of ||a||, 200 draws per dose, at 96, 107, 128 and 157
cells, both libraries, origin (i) per SENG3's amendment B2, with origin (ii)
reported beside it. Report the unresidualized power from SENG3 in the same
table.

Gate, wording fixed now, on P1. p50 is the lowest power at a 50% plant
across the eight library-by-cell-count combinations, p25 the same at 25%.
SENG3 already fired LIMITED, so 50% is the dose that matters.
 - RETAINED (p50 >= 0.80): "residualizing an independent proliferation axis
   leaves the statistic able to detect a planted half-gap approach, lowest
   power <p50>, so a residualized null is a null and not a loss of
   sensitivity."
 - LOST (p50 < 0.80): "residualizing costs the statistic its sensitivity:
   power at a 50% plant falls from 0.940 to <p50>. The residualized counts
   in D3 cannot be read as evidence of absence, and neither can SENG3's."
If LOST fires, D3 is reported but carries no conclusion, and that is said in
the gate sentence.

--- D3. The residualized approach test (gated) ---

Carries D1 and D2's readings, each named, and SENG3's LIMITED.

For all 340 perturbations and separately for the four named factors, recompute
the SENG3 D4 approach test with each of P1, P2 and P3 projected out of the
displacement, the target and the null identically. Floors are rebuilt from the
reserved NT half at each perturbation's own cell count, 200 draws, 5th
percentile. Both origins, neither privileged, gate on origin (i).

Report, per axis and per origin: below floor, delta negative, called, against
the chance expectation of 17.0 of 340 with its 95% binomial interval; per
modality; the four factors individually with cell count and percentile; and
the unresidualized SENG3 counts beside every row.

Also report, for the four factors, their delta percentile among the 60
perturbations in their own 96-157 cell band, residualized, beside SENG3's
unresidualized percentiles (3.3%, 15.0%, 21.7%, 1.7%).

Gate, wording fixed now, on the PRIMARY axis P1, counting the four named
factors called under each origin separately. Branches are mutually exclusive.
 - SURVIVES (both origins call two or more): "<k> of Junevity's four
   factors still end closer to early-passage cells after an independent
   proliferation axis is projected out: <list>. The approach is not
   accounted for by proliferation."
 - EXPLAINED (both origins call fewer than two): "after an independent
   proliferation axis is projected out, <k> of the four are called under
   origin (i) and <k2> under origin (ii), against <e> expected by chance.
   What the distance statistic scores as approach to early-passage cells is
   in this dataset a change in proliferative state."
 - SPLIT (the origins disagree): report both, name which factor is which
   under each, draw no single conclusion.
If the secondary axes P2 or P3 would fire a different branch than P1, say so
prominently beside the gate and report all three branch sentences, with P1's
as the task's reading.

--- D4. Report-only, no decision rule ---
 - The screen-wide count under each of P1, P2, P3 and the circular axis side
   by side, so the effect of the reference choice alone is visible.
 - rho between the residualized delta and the proliferation shift, per axis.
 - The cell-count decile table of SENG3's D5, recomputed residualized.
 - Fibroblast identity for the four and the screen, unchanged recipe.
 - The frozen ruler count: how many of 340 score younger on the
   one-dimensional ruler while not approaching, residualized, both origins.
   A diagnostic, not evidence about any factor.
 - A single table of every SENG2, SENG3 and SENG4 headline number.

--- WHAT DOES NOT COUNT ---
- choosing an axis, a space, an origin, a criterion, a null or a factor set
  after seeing a number;
- treating the secondary axes as primary because they give a preferred
  answer;
- treating any factor called here as validated or as a pilot candidate in
  this findings file;
- reading a null as evidence about transcription factors in aged human
  cells: this dataset is passaged neonatal-line fibroblasts and its "old" is
  passage number, which this project has found does not resemble donor aging
  (cos(v,u) = -0.040, permutation p = 1.0);
- refitting any axis to avoid reporting a circularity;
- any claim about Junevity's published numbers being unreproducible.

--- STOP KEYS ---
prolif_recipe_unreadable; prior_results_not_reproduced; no_raw_counts.
Seeds: 20260914 splits/nulls/permutations, 20260918 bootstraps.
n_null = n_perm = n_boot = 200.
```

## 2. Deviations and decisions

The flag was not ambiguous in a way that needed asking. Nine points of
record, each fixed in writing before computing:

1. **The NT half split could not be made disjoint from SENG3's halves.**
   The flag anticipates this ("if it is not, say so and use a fresh split").
   SENG3's two halves exhaust each NT pool, so no half-split of the same pool
   can avoid reusing cells. A fresh split at seed 20260914 is therefore used,
   and it is **identical to SENG3's split** because the seed and the pool are
   the same. Half one (333 CRA, 739 CRI) builds P1 and nothing else; half two
   (333 CRA, 740 CRI) is the origin and every floor. SENG3 used these same
   halves for its D1 null validation, which is a different test from the one
   an axis is applied to here.
2. **The origin changes, so an unresidualized control column was added.**
   Because half two alone may define the origin, every D3 number uses the
   half-two NT centroid, not SENG3's full NT pool. That alone would move the
   counts, so an unresidualized column (`NONE`) is computed with the same
   half-two origin. **This is an addition beyond the flag**, and it is what
   makes the comparison isolate residualization: `NONE` gives
   82 called against SENG3's 80, so the origin
   change costs about 2 of 340.
3. **P3 is built on measured ruler genes only.** 94 of the 97 Tirosh genes
   are present (`MLF1IP`, `FAM64A`, `HN1` absent from this file — the same
   three SENG recorded). Genes absent from the file are z=0 and carry no
   information, so including them would put weight on coordinates the data
   cannot move. Disclosed as a choice.
4. **The Southard comparison the flag asks for is degenerate and is reported
   as such.** Southard's own old-minus-young axis is built from GTEx donor age
   bins, so its cosine with the donor-age direction u is -1.0
   — the same vector up to sign. It is a definitional identity, not an
   independent comparison, and nothing is inferred from it.
5. **`cos(a, u)` is reported in two spaces** because the flag does not say
   which: full-gene (comparable to SENG's stored −0.040) and S2A_k20.
6. **The axis is a pseudobulk difference, not a per-cell regression**, in all
   four cases, because that is SENG's recipe (item (a) below).
7. **`.X` is streamed, not loaded.** Library sizes are recomputed over all
   36,601 file genes rather than taken from the stored `n_counts` column,
   which SENGINV showed disagrees with `.X` on 7,113 rows; the maximum
   disagreement is 22 counts and it is confirmed again here.
8. **SENG2's bootstrap CIs remain un-re-derived.** No gate reads a CI.
9. This task reuses the frozen helpers in `src/seng_run.py`,
   `src/seng2_common.py` and `src/seng3_run.py`, and the frozen basis from
   `results/south7/`. A reproduction that reuses the original implementation
   tests the record and the data path, not the implementation.

## 3. State the method before using it

### (a) How SENG3 constructs its proliferation axis

**`prolif_recipe_unreadable` does not fire.** Readable end to end off
committed code, in three steps.

Per-cell cell-cycle score, `src/seng_run.py:915` and `:1706`:

```python
def module_score(counts_csr, cols, lib):
    """Mean log1p(CP10k) over the listed columns. counts_csr is cells × genes."""
    sub = counts_csr[:, cols]
    scale = np.maximum(lib, 1.0)[:, None]
    return np.log1p(arr / scale * 1e4).mean(axis=1)
...
score = module_score(csr, s_sub, lib_total) + module_score(csr, g_sub, lib_total)
```

Quartile sets within each modality's NT cells, `src/seng_run.py:925` and `:1716`:

```python
q25, q75 = np.quantile(scores[finite], [0.25, 0.75])
bottom = finite & (scores <= q25);  top = finite & (scores >= q75)
...
nt_idx = np.flatnonzero(nt_mask & (modality == mod))
bottom, top, info = quartile_masks(score[nt_idx])
q_idx[mod] = dict(bottom=nt_idx[bottom], top=nt_idx[top])
```

The axis itself, `src/seng3_run.py` (`space3`) and `src/seng_run.py:1844`:

```python
pv = Z[ni[f"TOP_{m}"]] - Z[ni[f"BOT_{m}"]]
p_unit[m] = pv / float(np.linalg.norm(pv))
```

- **fit on**: the NT control cells of that modality, **all of them** — 666
  CRA and 1,479 CRI.
- **which cells**: the top and bottom quartiles by score, 167/167 for CRA and
  370/370 for CRI (`results/seng/prolif_quartiles.json`).
- **which genes**: the score uses the Tirosh S and G2M lists (43 + 54, of
  which 42 and 52 are in this file); the **axis** is a difference of full
  pseudobulk z vectors over all 23,485 ruler coordinates, not restricted to
  cell-cycle genes.
- **per modality**: yes, one axis each for CRA and CRI.
- in SENG3 it is then projected into S2A_k20 and renormalised.

**The circularity this task exists to remove**: those 666 and 1,479 NT cells
are exactly the cells whose centroid is SENG3's origin (i) and from which
every SENG3 floor is drawn.

### (b) How SENG3 residualizes

From `src/seng3_run.py`, inside `d4`'s `stats`:

```python
dr = d20 - np.outer(d20 @ phat, phat)      # the displacement
rr = r20 - float(r20 @ phat) * phat        # the origin-to-target vector
dlr = ||rr + dr|| - ||rr||
```

So the unit proliferation direction is projected out of **both** the
displacement and the origin-to-target vector, and the delta is recomputed in
the residual geometry. **The null is residualized by the same operation**:
the null draws are passed through the identical `stats` call, so the floor is
a residualized floor. That part of SENG3 is sound; this task's finding is not
that the operation was wrong.

### (c) What cells enter each axis, and whether they enter the test

| axis | fit on | n cells | enters origin? | enters target? | enters any floor? |
|---|---|---|---|---|---|
| **P1** (primary) | top vs bottom quartile of NT **half one**, per modality | 168 CRA + 370 CRI (of 333 / 739) | **no** | **no** | **no** |
| **P2** | top vs bottom quartile of PD26 | 940 (of 1880) | **no** | **no** | **no** |
| **P3** | the Tirosh gene set, equal weights | **0 cells** | **no** | **no** | **no** |
| CIRC (SENG3's) | top vs bottom quartile of the **whole** NT pool | 334 CRA + 740 CRI | **yes** | no | **yes** |

**Confirmed in writing: for P1, P2 and P3, none of the cells entering the
axis enters the origin, the target or any floor of the test it is applied
to.** Half one builds P1 and nothing else; half two is the origin and every
floor. PD26 appears nowhere else in this project's Sengstack tests. P3 uses
no cell at all.

**Circularities that remain, not refitted away:**

1. The aging axis `a` and the target are both built from PD14 and PD32 cells,
   as in SENG2 and SENG3. D1's gate quantity is a cosine involving `a`.
2. The frozen GTEx mu and sd used for z-scoring come from all 652 donors,
   anchors included — inherited from the Southard instrument.
3. The S2A_k20 basis is read frozen from 308 GTEx donors aged 40–59 and is
   fit on no Sengstack cell, so it adds no circularity.
4. **The CIRC column carries SENG3's circularity by construction** and is
   reported only as the comparison the task is about.

### (d) The four factors' cell counts

Re-derived under SENG's 804-column unique-strict-maximum rule: 128, 107, 96,
157 — all exact, in the Stage 0 table.

## 4. Stage 0 — reproduce first (gate)

**Gate: PASS.** All 22 checks;
counts exact, numerics within 1e-4 × max(1, |stored|).

| item | kind | stored | re-derived | ok |
|---|---|---|---|---|
| `SENG3 D1 eight rows` | count | 8 | 8 | yes |
| `SENG3 D1 max firing rate` | numeric | 0.085000 | 0.085000 | yes |
| `SENG3 D2 p25 worst` | numeric | 0.565000 | 0.565000 | yes |
| `SENG3 D2 p50 worst` | numeric | 0.940000 | 0.940000 | yes |
| `SENG3 D3 r` | numeric | 0.887005 | 0.887005 | yes |
| `SENG3 D3 cos per library A` | numeric | 0.811375 | 0.811375 | yes |
| `SENG3 D3 cos per library I` | numeric | 0.887005 | 0.887005 | yes |
| `SENG3 D4 n` | count | 340 | 340 | yes |
| `SENG3 D4 called origin(i)` | count | 80 | 80 | yes |
| `SENG3 D4 called origin(ii)` | count | 73 | 73 | yes |
| `SENG3 D4 named called origin(i)` | count | 4 | 4 | yes |
| `SENG3 D4 named called origin(ii)` | count | 4 | 4 | yes |
| `SENG3 D5 rho(prolif_shift, delta_i)` | numeric | -0.863900 | -0.863943 | yes |
| `SENG3 D5 residualized called origin(i)` | count | 25 | 25 | yes |
| `SENG3 D5 residualized called origin(ii)` | count | 8 | 8 | yes |
| `SOUTH9 res power_25 @106 S2A_k20` | numeric | 0.990000 | 0.990000 | yes |
| `SOUTH9 raw power_25 @106 S2A_k20` | numeric | 0.980000 | 0.980000 | yes |
| `n_cells[CRA_E2F3]` | count | 128 | 128 | yes |
| `n_cells[CRA_EZH2]` | count | 107 | 107 | yes |
| `n_cells[CRI_STAT3]` | count | 96 | 96 | yes |
| `n_cells[CRI_ZFX]` | count | 157 | 157 | yes |
| `raw integer counts (no_raw_counts)` | count | 0 | 0 | yes |

The eight-row SENG3 D1 table is reproduced in full at
`results/seng4/stage0_seng3_d1_eight_rows.csv`. `no_raw_counts` does not fire.

## 5. D1 — three proliferation axes, and whether aging is separable (gated)

All four axes are unit vectors in S2A_k20, read frozen.

### Each axis against the aging axis `a`

| axis | modality | cos with `a` | angle | fraction of ‖a‖ surviving | % of ‖a‖ removed | cells in the axis |
|---|---|---|---|---|---|---|
| P1 | CRA | +0.941708 | 19.66° | 0.3364 | 66.4% | 333 |
| P1 | CRI | +0.936287 | 20.56° | 0.3512 | 64.9% | 739 |
| P2 | CRA | +0.948360 | 18.49° | 0.3172 | 68.3% | 1880 |
| P2 | CRI | +0.948360 | 18.49° | 0.3172 | 68.3% | 1880 |
| P3 | CRA | +0.474385 | 61.68° | 0.8803 | 12.0% | 0 |
| P3 | CRI | +0.474385 | 61.68° | 0.8803 | 12.0% | 0 |
| CIRC | CRA | +0.937209 | 20.41° | 0.3488 | 65.1% | 666 |
| CIRC | CRI | +0.938874 | 20.14° | 0.3443 | 65.6% | 1479 |

**This is the quantity the project handoff asked for** as a measure of how
much of Junevity's signal is structurally confounded. The angle between the
passage-aging direction and an independently estimated proliferation
direction is **19.66°** (P1, CRA; 20.56°
for CRI). For the axis that touches no cell of this dataset at all, P3, it is
61.68°.

### Pairwise cosines between the axes

**CRA**

| | P1 | P2 | P3 | CIRC |
|---|---|---|---|---|
| **P1** | +1.0000 | +0.9641 | +0.6183 | +0.9965 |
| **P2** | +0.9641 | +1.0000 | +0.5968 | +0.9586 |
| **P3** | +0.6183 | +0.5968 | +1.0000 | +0.6049 |
| **CIRC** | +0.9965 | +0.9586 | +0.6049 | +1.0000 |

**CRI**

| | P1 | P2 | P3 | CIRC |
|---|---|---|---|---|
| **P1** | +1.0000 | +0.9385 | +0.5513 | +0.9996 |
| **P2** | +0.9385 | +1.0000 | +0.5968 | +0.9412 |
| **P3** | +0.5513 | +0.5968 | +1.0000 | +0.5495 |
| **CIRC** | +0.9996 | +0.9412 | +0.5495 | +1.0000 |

**The line that matters most here.** P1 and SENG3's circular CIRC are cos
0.9965 (CRA) and 0.9996 (CRI) apart — essentially the same direction.
Reserving half the NT pool to estimate the axis changed the axis almost not
at all. P2, built on PD26 cells that appear nowhere else, is also close
(0.964 / 0.939). Only P3, the pure gene-set direction, is distinct from the
three cell-based axes (0.55–0.62).

### `a` against the GTEx donor-age direction

- full-gene space: cos(a, u) = -0.040273, which reproduces
  SENG's stored −0.040 and is the basis of this project's standing finding
  that passage is not donor aging.
- S2A_k20: cos(a, u) = -0.199293.
- the Southard equivalent: -1.0, which
  is a definitional identity rather than a comparison — Southard's
  old-minus-young axis **is** the donor-age direction up to sign. Reported
  because the flag asks; nothing is inferred from it.

So the aging axis in this dataset is near-orthogonal to donor age and
near-parallel to proliferation. Both facts are about the dataset, not about
any factor.

### The D1 gate

**INSEPARABLE** fired, on |cos(P1, a)| = 0.941708, the larger
absolute value across modalities (CRA).

> "in this dataset passage aging and proliferation are substantially the
> same direction (cos +0.941708). Projecting proliferation out removes most
> of what defines young, so no residualization here can distinguish a
> proliferation artefact from passage rejuvenation. The question is not
> answerable in this dataset, and that is the finding."

**D2 and D3 are still run to the end and reported, and nothing in D3 is a
conclusion about any factor.** That is repeated wherever a factor is named.

## 6. D2 — does residualization leave the statistic able to see an approach? (gated)

Run **before** any residualized approach number. Plant, null and floor are
all residualized by the same projection. Draws come from the reserved NT half
two; origin (i) is its centroid.

### Lowest power across the eight library-by-cell-count combinations

| axis projected out | origin | p25 | **p50** | p100 |
|---|---|---|---|---|
| none (control) | (i) | 0.450 | **0.895** | 1.000 |
| none (control) | (ii) | 0.635 | **0.985** | 1.000 |
| P1 | (i) | 0.030 | **0.010** | 0.000 |
| P1 | (ii) | 0.385 | **0.760** | 0.945 |
| P2 | (i) | 0.060 | **0.040** | 0.000 |
| P2 | (ii) | 0.310 | **0.715** | 0.880 |
| P3 | (i) | 0.440 | **0.885** | 1.000 |
| P3 | (ii) | 0.655 | **0.990** | 1.000 |
| none, SENG3's full-NT origin | (i) | 0.565 | **0.940** | 1.000 |

Read the control rows first. Changing the origin to the reserved half costs
almost nothing: p50 goes from 0.940 to
0.895. Residualizing P3, which sits
62° from `a`, also costs almost
nothing (0.885). Residualizing P1 or P2,
which sit about 19–21° from `a`, destroys the statistic.

Full dose curves are in `results/seng4/d2_plant_residualized.csv`. Per-depth
detail for P1 under origin (i):

| library | n cells | floor | power @25% | power @50% | power @100% | SENG3 unresid. @50% |
|---|---|---|---|---|---|---|
| CRA | 96 | -2.693 | 0.050 | **0.020** | 0.000 | 0.950 |
| CRA | 107 | -2.486 | 0.035 | **0.030** | 0.000 | 0.940 |
| CRA | 128 | -2.858 | 0.030 | **0.020** | 0.000 | 0.975 |
| CRA | 157 | -2.860 | 0.030 | **0.010** | 0.000 | 0.990 |
| CRI | 96 | -2.481 | 0.060 | **0.035** | 0.000 | 0.970 |
| CRI | 107 | -2.239 | 0.065 | **0.055** | 0.000 | 0.975 |
| CRI | 128 | -1.992 | 0.085 | **0.070** | 0.005 | 1.000 |
| CRI | 157 | -2.424 | 0.070 | **0.045** | 0.010 | 1.000 |

### The D2 gate

**LOST** fired, on P1, p50 = 0.010 against the 0.80
requirement.

> "residualizing costs the statistic its sensitivity: power at a 50% plant
> falls from 0.940 to 0.010. The residualized counts in D3 cannot
> be read as evidence of absence, and neither can SENG3's."

**D3 is reported but carries no conclusion.** That sentence also retires
the only quantitative support SENG3 had for the claim that its approach
signal is proliferation: SENG3's 80 → 25 and 73 → 8 are drops produced by
an operation that removes the signal along with the confound.

## 7. D3 — the residualized approach test (gated)

Carries D1 **INSEPARABLE**, D2 **LOST** and SENG3's **LIMITED**.
**Because D2 fired LOST and D1 fired INSEPARABLE, nothing in this section is
a conclusion about any factor.**

Chance expectation is 17.0 of 340 by construction.

| axis projected out | origin | below floor | delta negative | **called** | CRA | CRI | 95% interval |
|---|---|---|---|---|---|---|---|
| none (control) | (i) | 82 | 208 | **82** | 36 | 46 | 19.7–29.0% |
| none (control) | (ii) | 78 | 182 | **78** | 32 | 46 | 18.6–27.8% |
| P1 | (i) | 27 | 97 | **27** | 24 | 3 | 5.3–11.3% |
| P1 | (ii) | 13 | 20 | **9** | 4 | 5 | 1.2–5.0% |
| P2 | (i) | 45 | 141 | **45** | 26 | 19 | 9.8–17.3% |
| P2 | (ii) | 19 | 13 | **10** | 1 | 9 | 1.4–5.3% |
| P3 | (i) | 75 | 192 | **75** | 34 | 41 | 17.8–26.8% |
| P3 | (ii) | 72 | 160 | **72** | 27 | 45 | 17.0–25.9% |
| SENG3, circular axis | (i) | — | — | **25** | — | — | — |
| SENG3, circular axis | (ii) | — | — | **8** | — | — | — |
| SENG3, unresidualized | (i) | 80 | 191 | **80** | 30 | 50 | — |

**The comparison the task was built to make.** The independent P1 gives
27 called against SENG3's circular 25 under
origin (i), and 9 against 8 under origin (ii).
Removing the circularity moved the count by about 2 and 1. The circularity
was real and was correctly disclosed by SENG3, but it was not what made the
adjustment questionable.

### The four named factors

| factor | n cells | pctl | delta, none | delta, P1 | floor, P1 | called P1 (i) | called P1 (ii) | called P2 (i) | called P3 (i) |
|---|---|---|---|---|---|---|---|---|---|
| `CRA_E2F3` | 128 | 96% | -21.757 | -4.162 | -2.858 | **yes** | no | yes | yes |
| `CRA_EZH2` | 107 | 89% | -15.756 | -2.975 | -2.486 | **yes** | no | yes | yes |
| `CRI_STAT3` | 96 | 82% | -10.021 | +1.574 | -2.481 | no | no | no | yes |
| `CRI_ZFX` | 157 | 99% | -19.623 | +2.362 | -2.424 | no | no | yes | yes |

### Their delta percentile among the 60 perturbations in their own depth band

| factor | SENG3, unresidualized | SENG4, P1-residualized |
|---|---|---|
| `CRA_E2F3` | 3.3% | 1.7% |
| `CRA_EZH2` | 15.0% | 10.0% |
| `CRI_STAT3` | 21.7% | 75.0% |
| `CRI_ZFX` | 1.7% | 78.3% |

In the band, 7 of 60 are called after P1 is
projected out, against 24 of 60 unresidualized in SENG3. `CRI_ZFX` moves from
the most extreme of the band (1.7%) to the 78th percentile; `CRI_STAT3` from
21.7% to 75.0%. `CRA_E2F3` stays extreme (1.7%).

### The D3 gate

**SPLIT** fired on the primary axis P1: origin (i) calls 2 of
the four, origin (ii) calls 0.

> SPLIT. Under origin (i) 2 are called: CRA_E2F3, CRA_EZH2.
> Under origin (ii) 0 are called: none. Both
> are reported and no single conclusion is drawn.

**The secondary axes do not all agree, and the flag requires that be said
prominently.**

| axis | branch it would fire | origin (i) / (ii) of the four |
|---|---|---|
| P1 (primary) | **SPLIT** | 2 / 0 |
| P2 | **SPLIT** | 3 / 0 |
| P3 | **SURVIVES** | 4 / 4 |
| none (control) | **SURVIVES** | 4 / 4 |

P2 fires the same branch as P1. **P3 fires SURVIVES**, a different branch,
and it is the axis furthest from `a`. That is exactly what D1 and D2 predict:
an axis at 62° from aging can be projected out without removing the aging
signal, so the factors survive it; the axes that are near-parallel to aging
cannot. P1's branch is this task's reading, and under D2's LOST it is not a
conclusion about any factor.

## 8. D4 — report-only, no decision rule

### Screen-wide counts by reference choice, and rho with proliferation

| axis projected out | called, origin (i) | called, origin (ii) | rho(prolif. shift, delta) |
|---|---|---|---|
| none (control) | 82 | 78 | -0.8583 |
| P1 | 27 | 9 | -0.0667 |
| P2 | 45 | 10 | -0.3126 |
| P3 | 75 | 72 | -0.8123 |

The rho column is the cleanest summary of what each projection does: P1
removes the proliferation correlation almost exactly (-0.0667
from -0.8583), P2 partly
(-0.3126), P3 barely
(-0.8123). P1 does the job it was built to do.
It also removes the aging signal, which is D2.

### Cell-count deciles, residualized

| decile | n | median cells | call rate, none | call rate, P1 |
|---|---|---|---|---|
| 1 | 36 | 33 | 0.056 | 0.083 |
| 2 | 34 | 39 | 0.206 | 0.059 |
| 3 | 33 | 47 | 0.061 | 0.061 |
| 4 | 37 | 53 | 0.189 | 0.081 |
| 5 | 32 | 60 | 0.125 | 0.156 |
| 6 | 34 | 68 | 0.265 | 0.029 |
| 7 | 32 | 77 | 0.312 | 0.062 |
| 8 | 34 | 88 | 0.382 | 0.029 |
| 9 | 35 | 99 | 0.286 | 0.029 |
| 10 | 33 | 121 | 0.545 | 0.212 |

The cell-count gradient SENG3 reported (0.056 to 0.545) is still present
unresidualized and is largely flattened after P1 is projected out. That is
consistent with the gradient being partly a proliferation gradient, and is
not evidence for it, because D2 shows the same projection removes the signal.

### Fibroblast identity

Unchanged recipe, ten identity genes: screen median
-0.0472; 1 of
340 exceed SENG's identity-loss threshold. For the four:

| factor | identity change | younger on the frozen ruler |
|---|---|---|
| `CRA_E2F3` | -0.2154 | yes |
| `CRA_EZH2` | -0.3298 | yes |
| `CRI_STAT3` | -0.1477 | yes |
| `CRI_ZFX` | -0.0380 | no |

### The frozen ruler — the preprint's thesis in this dataset

- unresidualized, origin (i): **119 of 181** perturbations
  score younger on the one-dimensional ruler while not approaching
  (82 called in total).
- unresidualized, origin (ii): **121 of 181** perturbations
  score younger on the one-dimensional ruler while not approaching
  (78 called in total).
- P1-residualized, origin (i): **168 of 181** perturbations
  score younger on the one-dimensional ruler while not approaching
  (27 called in total).
- P1-residualized, origin (ii): **176 of 181** perturbations
  score younger on the one-dimensional ruler while not approaching
  (9 called in total).

A diagnostic, not evidence about any factor. The dissociation widens under
residualization only because fewer things are called at all.

### Every headline number, SENG2 → SENG3 → SENG4

| quantity | SENG2 (full genes) | SENG3 (S2A_k20) | SENG4 (S2A_k20, independent axis) |
|---|---|---|---|
| cos(PD14−NT_A, a) | +0.584155 | +0.811375 | unchanged (read frozen) |
| cos(PD14−NT_I, a) | +0.757055 | +0.887005 | unchanged |
| offset gate | USABLE (by `_I` only) | USABLE (both) | — |
| null firing rate on effect-free data | 0.000 (degenerate) | 0.040–0.085 | — |
| power, 25% plant, worst | 0.000 | 0.565 | 0.030 (P1) |
| power, 50% plant, worst | not computed | 0.940 | **0.010** (P1) |
| sensitivity gate | INSENSITIVE | LIMITED | **LOST** |
| called of 340, origin (i) | 3 | 80 | 27 (P1), 82 (none) |
| called of 340, origin (ii) | 0 | 73 | 9 (P1) |
| of the four, called (i) / (ii) | 1 / 0 | 4 / 4 | 2 / 0 (P1) |
| approach gate | NO APPROACH (voided) | APPROACH | **SPLIT** (carries no conclusion) |
| rho(prolif, delta) | — | −0.8639 | -0.0667 after P1 |
| |cos(prolif axis, a)| | — | not measured | **0.9417** |
| separability gate | — | — | **INSEPARABLE** |

## 9. The reading

Only the pre-registered readings whose conditions fired are conclusions:

1. **Stage 0: PASS.** Every stored SENG3 and SOUTH9 number re-derived.
2. **D1: INSEPARABLE**, |cos(P1, a)| = 0.941708,
   19.66°.
3. **D2: LOST**, p50 = 0.010 against 0.940.
4. **D3: SPLIT**, and under D2's LOST it carries no conclusion.

So the conclusion of this task is D1's, and it is a conclusion about the
dataset rather than about any factor:

> In this dataset passage aging and proliferation are substantially the same
> direction — 19.7° apart, measured with a proliferation axis
> estimated on control cells that enter no origin, no target and no floor of
> the test. Projecting proliferation out removes 66% of the
> aging axis's length and costs the statistic its ability to detect a planted
> half-gap approach (0.010 against 0.940). No residualization here can
> distinguish a proliferation artefact from passage rejuvenation. The
> question is not answerable in this dataset.

**What this does and does not do to SENG3.** It does not restore SENG3's
APPROACH as a positive result: SENG3's own LIMITED still stands, the
cell-count gradient still stands, and the near-parallelism measured here
means its approach statistic and a proliferation score are not
distinguishable in this dataset. What it does remove is the one piece of
evidence that was being used *against* SENG3 — the residualized 80 → 25 and
73 → 8. Those drops are now shown to be what residualizing a near-parallel
axis does to any signal, confounded or not, because the same operation
destroys a planted approach that is known to be real.

### What this does not cover

- **Not a result about aged human cells.** Passaged neonatal dermal
  fibroblasts, whose "old" is passage number. cos(a, u) =
  -0.040 in full genes, reproduced here; this project
  has already found passage does not resemble donor aging (permutation
  p = 1.0). Nothing here transfers to aged donor cells.
- **No factor named here is validated or a pilot candidate**, under any gate,
  including the two that P1 called under origin (i).
- **The secondary axes are not promoted.** P3 fires SURVIVES and is reported
  prominently, but P1 was fixed as primary before any number and remains so;
  P3's branch is explained by its 62° angle, not preferred for its answer.
- **Inseparability is a property of this dataset**, not a general claim about
  proliferation and aging. It is measured here on one passage series in one
  cell line in one 20-dimensional space.
- **The circularities in (c) are not removed**: `a` and the target are built
  from the same PD14/PD32 cells, and the z-scoring mu and sd come from all
  652 GTEx donors including the anchors.
- **The origin changed from SENG3's** (reserved half rather than the full NT
  pool), which is why an unresidualized control column is reported beside
  every residualized one.
- **No claim that Junevity's published numbers are unreproducible.**

### Files

- `results/seng4/d1.log` (656 B)
- `results/seng4/d1_pairwise.csv` (573 B)
- `results/seng4/d1_summary.json` (1,988 B)
- `results/seng4/d1_vs_aging.csv` (933 B)
- `results/seng4/d2.log` (882 B)
- `results/seng4/d2_plant_residualized.csv` (17,312 B)
- `results/seng4/d2_summary.json` (838 B)
- `results/seng4/d3.log` (1,219 B)
- `results/seng4/d3_deciles.csv` (775 B)
- `results/seng4/d3_stderr.txt` (0 B)
- `results/seng4/d3_stdout.txt` (1,219 B)
- `results/seng4/d3_summary.json` (11,036 B)
- `results/seng4/d3_table.csv` (215,953 B)
- `results/seng4/nt2_rows_CRA.npy` (2,792 B)
- `results/seng4/nt2_rows_CRI.npy` (6,048 B)
- `results/seng4/prep.log` (1,171 B)
- `results/seng4/prep.npz` (445,072 B)
- `results/seng4/prep_quartiles.json` (439 B)
- `results/seng4/PREREG.flag` (11,917 B)
- `results/seng4/stage0.log` (532 B)
- `results/seng4/stage0_reproduction.csv` (1,168 B)
- `results/seng4/stage0_seng3_d1_eight_rows.csv` (1,129 B)
- `results/seng4/stage0_summary.json` (3,243 B)

Nothing was committed and nothing was pushed.
