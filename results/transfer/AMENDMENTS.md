# TRANSFER — amendments to the pre-registration

Written **2026-10-13T19:13Z (local 2026-10-14 04:13 KST), before any TRANSFER
number existed.** `results/transfer/PREREG.flag` is **not edited**; it stands as
written at 2026-10-13T18:44:34Z. Everything below is reproduced verbatim in
`FINDINGS_TRANSFER.md`.

These are the author's answers to the hand-off questions this session raised
after writing the flag, plus the few definitions the answers required. The flag
is the primary document; where this file narrows or adds to it, that is marked
as a deviation and is disclosed again where the affected number appears.

---

## A1 — Reading order (author's answer; flag section 5 unchanged)

The ordered evaluation stands: UNMEASURABLE, then LARGE DRIFT, then SMALL
DRIFT, then INTERMEDIATE. The kill is evaluated first.

Recorded for the reader: **when the bootstrap interval is valid the branches are
already mutually exclusive**, so the order changes nothing. SMALL DRIFT requires
specificity at p < 0.01 while LARGE DRIFT's specificity arm requires p >= 0.05;
LARGE DRIFT's interval arm requires the CI upper bound below 0.5 while SMALL
DRIFT requires the lower bound at or above 0.7. The two cannot both hold. The
order matters only for an INVALID interval, where the CI arms are unreadable and
the kill may still fire on specificity alone.

## A2 — Stage 0(b) is narrowed (author's answer; **deviation from the flag**)

The flag's Stage 0(b) required the SENG3 floors and the floor-dependent boolean
columns to be recomputed. **They are not.** Stage 0(b) reproduces:

- the number of qualifying CRA perturbations (stored: 167);
- every CRA perturbation's `n_cells`;
- `delta_i`, `frac_i`, `delta_ii`, `frac_ii` for the CRA rows of
  `results/seng3/d4_table.csv`.

It does **not** recompute `floor_delta_i`, `floor_delta_ii` or the
`below_floor_*`, `decreased_*` and `called_*` columns.

Reason, as answered: this task never uses the floor machinery. It measures
direction agreement between two displacements and never compares a distance to a
null floor. Reproduction should test the path this task actually uses — the
assignment rule, the panel transform, and the projection into the 20-component
space — and the four reproduced numerics exercise exactly that path. The
narrowing is a **deviation from the flag** and is disclosed here, in the Stage 0
section of the findings, and beside the Stage 0 table.

Consequence to state plainly: Stage 0 no longer re-verifies SENG3's floors.
Nothing in TRANSFER depends on them, and no TRANSFER number is a check on them.

## A3 — The matched-gene count was already resolved (author's answer)

The apparent conflict this session raised between "20,265 matched genes"
(`results/coordspace/PREREG.flag`) and "20,264 genes" (the check labels in
`src/south6rec_run.py`) is **already resolved in `FINDINGS_COORDSPACE.md`**,
under the deviation headed "matched coordinates": 20,265 ruler genes match the
Southard file, two of them map to the same source column, SOUTH6 collapses
duplicate source columns (`src/south6_run.py`, the `seen` dict in
`stage0_and_1`), and the panel therefore has 20,264 columns. Both texts are
correct about different quantities. Nothing in TRANSFER changes.

`|G|` is still computed from the arrays at run time exactly as the flag requires
and is not quoted from any earlier findings file. That note is cited where `|G|`
is reported.

## A4 — No depth-matched control variant (author's answer)

The flag's `d_J` compares a perturbation of tens to hundreds of cells against
the pooled NT cells of the whole CRA arm, and `d_S` likewise against the pooled
78-guide non-targeting set. **No depth-matched variant is added, not even as a
diagnostic.** Unequal depth between factor and control is precisely what the
reliability correction in flag section 3 exists to absorb: the control side
contributes almost no sampling noise, the factor side contributes nearly all of
it, and the split-half reliability measures that noise at the factor's own
depth.

Stated as a limitation, not corrected: the reliability correction assumes the
noise it measures by splitting is the same noise that limits the full-depth
displacement. Control-side noise, being shared between both halves, is not
measured by the split and is not corrected for.

## A5 — No separate feasibility probe (author's answer)

No cheap pre-flight probe of whether `|F_q| >= 20` is pre-registered. With A2's
narrowing, Stage 0 is cheap, and **the Stage 4 same-context gate is itself the
feasibility check**. If the readout cannot recover zero drift within a single
screen, or if fewer than 20 factors reach the reliability floor, the result is
UNMEASURABLE and that is a finding about the instrument, reported as such.

## A6 — The 0.3 reliability floor, and why that number (author's answer)

`r >= 0.3` in both screens stands, and the rationale is recorded now rather than
left unstated: the correction divides by `sqrt(r_S * r_J)`, so a floor of 0.3 in
both screens **caps the inflation factor at about 3.3x**
(`1 / sqrt(0.3 * 0.3) = 3.33`). A factor whose split halves barely agree in
either screen therefore cannot produce a large corrected cosine out of a small
observed one. The floor was fixed before any data existed and is not tuned to
any result.

---

## Carried over unchanged from the hand-off

Accepted as this session fixed them in the flag, with no change: the Southard
factor universe (all 1,836 `keep` labels, SOUTH7 primary-set membership as a
descriptive column only); the exact-symbol matching rule with no alias table;
the half-split displacement convention (factor half against an independent
control half, which differs from SENG3 amendment B5's within-set paired
construction and is disclosed as a deviation from it); the pooled 78-guide NT
control on the Southard side; the split mechanics, the extra cell to half A, the
shared per-split control halves and the per-factor seed streams; seed 20260914
for every draw including the bootstrap, departing from the project's usual
20260918 bootstrap seed; the treatment of non-finite `c_corr`; the specificity
test details and its disclosure as approximate; the quartering of both factor
and control cells inside the gate; the real re-accumulation of the Southard
panel in Stage 0(a); the CRI cell-count checks in Stage 0(c) as a structural
check on an arm that is not tested; no PROGRESS file; reading the Sengstack
panel counts from `results/seng2/panel_counts.npz` while re-deriving cell
assignment and counts from the h5ad; and no activity filter on either screen's
side, so a shared gene is in `F` whether or not Southard's guides for it were in
the published active set.

## Implementation definitions fixed here, before computing

These are mechanical and are listed so the record is complete.

1. **Control split stream.** The flag fixes the per-factor stream as
   `default_rng([20260914, screen_id, factor_index, split_index])` but does not
   give the form for the control halves, which are shared across factors within
   a split index. They are drawn from `default_rng([20260914, screen_id,
   split_index])` — a three-element key, so it cannot collide with any
   four-element per-factor key.
2. **The two "contexts" in the Stage 4 gate.** The gate's two context vectors
   come from split index 0 of a dedicated gate stream
   (`default_rng([20260914, screen_id, 7, factor_index])` for factor cells and
   `default_rng([20260914, screen_id, 7])` for control cells); the quarter-split
   reliabilities of each context are averaged over the same 100 splits the
   primary readout uses, with quarter streams keyed
   `[20260914, screen_id, 8, factor_index, split_index]`.
3. **Per-stage caching.** Each stage writes its own outputs to
   `results/transfer/` and a later stage may read them. This is not
   checkpoint/resume of a random stream: every draw stream inside a stage runs
   uninterrupted from its seed, and re-running a stage re-draws it identically.
4. **Source-file verification.** The Sengstack h5ad is opened, so its size and
   md5 are verified against `src/seng2_common.py` `MD5_EXPECT`
   (`bfccaa94a916c22e4b0793913a04128b`). The **Southard h5ad is never opened by
   this task** — counts come from the SOUTH6 matched CSR cache — so its md5 is
   not recomputed; Stage 0(a)'s re-accumulation of `C_f`, `C_nt` and `C_pool`
   from that cache, compared against the stored panel at tolerance 0, is the
   check that stands in its place. Disclosed as a deviation from the flag's
   `file_mismatch` stop key, which can therefore fire only on the Sengstack
   file.
