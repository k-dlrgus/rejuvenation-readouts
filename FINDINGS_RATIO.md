# FINDINGS_RATIO — is the typical perturbation short enough, against its own experiment's gap, for arrival to register?

**Task:** RATIO. **Status:** complete for what was workable; no pre-registered STOP fired.
GSE92872 was excluded by the user before the flag (U1). GSE278572 was analysed in full. GSE314342
was **not analysed**: under the flag's D13 its 171.2 GB deposit projects to ~46.7 h of download
alone against a 24 h limit, so the task stopped after GSE278572, and no subset of it was analysed.

**Readings that fired (pre-registered wording; one space only, see below):**

```
RATIO-OPEN: in GSE278572 Treg, 19 of 19 perturbations have ratio < 2 in full measured-gene space; arrival was geometrically attainable for them. [within_floor = 9]
RATIO-OPEN: in GSE278572 Teff, 17 of 17 perturbations have ratio < 2 in full measured-gene space; arrival was geometrically attainable for them. [within_floor = 11]
```

RATIO-BLOCKED did not fire in either sub-screen: no perturbation of either reached ratio 2, and none
came within a factor 2.6 of it. RATIO-UNMEASURABLE did not fire in either: 19 of 28 (Treg) and 17 of
28 (Teff) perturbations exceed their size-matched control floor, against a chance expectation of 1.4
(one-sided binomial p = 8.5e-19 and 9.6e-16).

**Descriptive lines (not gated):** 28 of 28 eligible perturbations have ratio < 2 and therefore could
in principle have registered an arrival, in each sub-screen; 9 (Treg) and 11 (Teff) of those are
within_floor.

**The part that undercuts this, stated up front.** In this screen, ratio < 2 is close to guaranteed
by depth alone. Random control-guide cells, which carry no effect, already reach ratio 0.45 (Treg) and
0.43 (Teff) at the median perturbation size (~1,000 cells). The perturbations sit only just above that:
median ratio / own-size floor P95 = 1.03 and 1.02, maximum 1.73 and 1.48. At 30 cells, the pre-registered
minimum, control cells alone reach 1.74 / 1.52 (P50) and 1.82 / 1.58 (P95), so a 30-cell screen of this
experiment would sit near the ratio-2 line on noise alone. RATIO-OPEN here is a true statement about
geometry. It is **not** evidence that perturbations are short relative to the gap: most of each
measured length is sampling noise. A rough noise-corrected length (diagnostic, not pre-registered,
below) has medians of 0.15 and 0.13 gap lengths.

**The task's premise, checked against the stored record (R1b).** "Every perturbation's displacement
exceeded twice the old-to-young gap" holds for SOUTH6 in the full-gene space S1 (1,836 of 1,836,
minimum ratio 2.097). It does **not** hold in SOUTH6's compressed spaces: 1, 2 and 2 of 1,836 reach 2
in S2_k20/k50/k100 (minima 0.13–0.24).

**What this does not answer.** Whether short-enough perturbations are "a general property of
perturbation screens". One experiment with two sub-screens was measured. Its gap is a T-cell
restimulation gap, ~145–152 log2-CPM units in its own space; it is not an aging gap, and it is
measured in a different coordinate system from SOUTH6's. The SOUTH6 S1 ratio was never set beside a
size-matched control floor of this kind, so whether its ≥2 was itself depth-driven is not known from
this task.

Files: `results/ratio/PREREG.flag`; `readings.txt`; `summary.csv`; per sub-screen
`gse278572_{treg,teff}_{perturbations,floor,posctrl,excluded,diagnostic}.csv` and
`_floor_draws.npz`; `gse278572_assignment.json`, `gse278572_cells_per_target.csv`;
`gene_overlap.csv`; `reproduction.csv`; `diagnostic_summary.csv`; `gse314342_workability.json`
(plus the invalid first attempt, `gse314342_workability_attempt1_badurl.json`) and
`gse314342_filelist.txt`; logs `ratio_run_*.log`, `gse278572_stdout.txt` and the crashed first
attempt `gse278572_stdout_attempt1_crashed.txt`. Code: `src/ratio_run.py`.

## Pre-registration (verbatim copy of `results/ratio/PREREG.flag`)

```
RATIO pre-registration flag
Written 2026-10-25T19:24:21+09:00 (2026-10-25T10:24:21Z), before any statistic of this task
exists. No ratio, norm, pseudobulk, normalisation or floor has been computed. The expression
matrix of no dataset has been downloaded or opened.

What had been read or seen before this file was written (disclosed):
 - Repository files: CLAUDE.md, SCOPING_DISEASE.md (whole), FINDINGS_SOUTH6.md lines 1-511
   (this includes its stored Stage 4 table, so the stored SOUTH6 median ||d||/||gap|| per
   setting had been seen: S1 O1_Y1 3.8208; S2_k20/50/100 O1_Y1 0.4124/0.5849/0.6917),
   FINDINGS_SENG3.md lines 312-345 (what the S2A_k20 basis is), results/scopedis/datasets.csv
   rows for GSE92872/GSE278572/GSE314342, results/scopedis/opened/gse92872/arm_counts.txt, the
   headers of results/south6/stage3_full_table.csv, stage4_overall.csv, stage4_quartiles.csv,
   origin_target_distances.csv, stage1_notes.json, the key list of results/south6/stage1.npz,
   and code in src/south6_run.py, src/south7_run.py, src/seng_run.py, src/toward_run.py,
   src/fibro2_common.py, .gitignore.
 - NCBI GEO FTP directory listings (file names and sizes) for GSE92872, GSE278572, GSE314342,
   and GSE314342's filelist.txt.
 - Three small GSE278572 metadata files were downloaded into data/raw/ratio/gse278572/
   (barcodes.tsv.gz 1.1 MB, features.tsv.gz 0.3 MB, protospacer_calls_per_cell.csv.gz 1.9 MB) and
   only their structure was read: feature types and names (36,601 Gene Expression, 141 Antibody
   Capture including HTO_1 Resting_Treg / HTO_2 Stimulated_Treg / HTO_3 Resting_Teff /
   HTO_4 Stimulated_Teff, 65 CRISPR Guide Capture = 28 targets x 2 guides + 9 Non-Targeting),
   barcode suffixes -1..-8 with per-suffix barcode totals, and the calls file's header and first
   3 rows. No per-cell guide or HTO assignment was tabulated.
 - Four questions were put to the user in chat; their rulings are reproduced verbatim below (U1-U4).

== Task (verbatim, as given in chat) ==

Task: `RATIO`. Deliverables: `results/ratio/PREREG.flag`, `FINDINGS_RATIO.md`, `PROGRESS_RATIO.md`. One local commit, no push. Modify no existing FINDINGS_*.md and no SCOPING_*.md.

## The question

In the Southard screen, every perturbation's displacement exceeded twice the old-to-young gap, so no approach could register whatever its direction. This task asks whether that is a general property of perturbation screens: is the typical perturbation effect short enough, relative to the experiment's own diseased-to-unaffected gap, that arrival could register at all?

This question needs only two lengths, both measurable inside one experiment: the displacement of each perturbation from its control arm, and the gap between the diseased/activated state and the unaffected state. It does NOT need the unaffected state to serve as a rejuvenation target, and it does NOT depend on perturbation timing. A screen whose design cannot support a rejuvenation claim can still answer this.

## Datasets

The immune-activation screens SCOPEDIS recorded as `where_control_state_lives = in_run`, in this order, smallest raw data first:
1. GSE92872 (CROP-seq, 32 genes, stimulated and naive arms of one pool)
2. GSE278572 (28 regulators, TotalSeq hashtags encoding Rest/Stim, same wells)
3. GSE314342 (genome-scale CRISPRi, Rest/Stim8hr/Stim48hr, same 10x Flex lane)

Start with 1. Proceed to the next only when the previous is complete. GSE249595 is excluded: its readout is a 374-gene panel, flagged `insufficient_for_distance` in SCOPEDIS.

Raw-data size is the binding constraint. State each dataset's download size before downloading. If a dataset exceeds what is workable, say so, complete what is workable, and stop — do not substitute a different dataset.

## Pre-register before any statistic

Write the flag before computing anything, including the statistic, the spaces, the controls and the exact wording of every possible conclusion.

Statistic, per perturbation: ‖d‖ / ‖gap‖, where d is the perturbation pseudobulk minus its own screen's control-guide pseudobulk, and gap is the activated-state control pseudobulk minus the resting-state control pseudobulk, both from control-guide cells in the same run. Report the distribution: median, quartiles, min, max, and the fraction with ratio >= 2.

Two coordinate spaces, both reported, neither primary over the other without saying so: the full measured-gene space, and a compressed space built the way paper 2 builds it. If the compressed space cannot be built the same way for a non-fibroblast cell type, say so plainly and report the full space only — do NOT invent a substitute basis and do NOT reuse the fibroblast basis on T cells.

Fix now, before any number:
- A ratio >= 2 means no displacement of that length can end closer to the target, whatever its direction. State this as geometry, not as a claim about the perturbation's biology.
- Report, as a separate descriptive line, how many perturbations have ratio < 2 and therefore could in principle have registered an arrival.
- Minimum cells per perturbation: 30, per D7 carried from the SCOPEFIB flag. Perturbations below it are excluded and counted.
- A noise floor: compute the same ratio for random samples of control-guide cells at the same group sizes, which carry no effect. Report it beside the distribution. A perturbation ratio indistinguishable from the control floor means the screen cannot measure displacement length at that depth, and that reading must be available as a conclusion.

Possible conclusions, worded now:
  "RATIO-BLOCKED: in <dataset>, <n> of <N> perturbations have ratio >= 2 in <space>; arrival could not register for them whatever their direction."
  "RATIO-OPEN: in <dataset>, <n> of <N> perturbations have ratio < 2 in <space>; arrival was geometrically attainable for them."
  "RATIO-UNMEASURABLE: in <dataset>, the perturbation ratio distribution is not distinguishable from the control floor in <space>; this screen cannot measure displacement length at this depth."
Report every reading whose conditions fired, per dataset and per space.

## Hard rules

No conclusion about any named gene or factor. This is a property of experiments, not of biology; do not name hits.
State the gene overlap with the frozen 23,485-gene panel for each dataset, and say plainly where coordinates are not comparable to the fibroblast work.
Disclose every deviation, including anything decided that the flag did not specify.

== Rulings by the user in chat, before this file was written ==

The session raised four points before writing this flag. Questions are paraphrased; the text
after "Ruling:" is the user's answer, verbatim.

[U1] GSE92872. The task names GSE92872's CROP-seq arm as an in_run screen, but SCOPEDIS
recorded that arm `separate_experiment` (one Drop-seq run per condition) and only its 20-gene
arrayed bulk arm `in_run`. Which arm?
Ruling: Skip GSE92872.
Consequence: GSE92872 is not analysed in this task, under either arm. The run starts with
GSE278572 and then GSE314342.

[U2] Firing of RATIO-BLOCKED / RATIO-OPEN and relation to RATIO-UNMEASURABLE.
Ruling: Any count; UNMEASURABLE supersedes — as in the recommended option: BLOCKED fires if
n>=1 perturbations have ratio>=2; OPEN fires if n>=1 have ratio<2; both can fire. Where
UNMEASURABLE fires for a dataset x space, BLOCKED and OPEN are not reported as conclusions there
and their counts stay as descriptive lines.
Addition, fixed now: at the perturbation level, a perturbation whose ratio does not exceed the
95th percentile of the control floor at its own group size is counted as `within_floor` and
excluded from both n and N for BLOCKED and OPEN. Report the within_floor count beside each
reading. Reason: a ratio below 2 that is no larger than random control cells' is not a measured
displacement, and counting it as OPEN would inflate OPEN with perturbations that did nothing
measurable.

[U3] Deciding "indistinguishable from the control floor".
Ruling: Per-pert P95, binomial (Recommended) — the option as offered: for each perturbation,
1,000 random draws of n control-guide cells from the activated arm (n = its cell count),
referenced against the remaining control cells. It 'exceeds' if its ratio > the P95 of its
draws. UNMEASURABLE fires if the exceed count is not above the 5% chance rate (one-sided exact
binomial p>=0.05). The 5% baseline is reported beside the count.

[U4] Scaling in the full measured-gene space.
Ruling: log2-CPM, no z (Recommended) — the option as offered: sum counts per pseudobulk, TMM,
log2 CPM with prior.count=2 (as TOWARD/SAME), over every gene the dataset measures; no per-gene
z-score, since the frozen GTEx fibroblast mu/sd would reuse a fibroblast instrument on T cells.
Coordinates are stated as not comparable to the fibroblast work.

== The compressed space: not built (decided now, stated to the user in chat before U1-U4) ==

Paper 2's compressed space (S2A_k20; FINDINGS_SENG3.md B1, src/south7_run.py anchor_free_basis)
is the first 20 PCA loadings of the 308 GTEx cultured-fibroblast donors aged 40-59, in the
frozen-ruler 23,485-gene GTEx-fibroblast z space. Building it "the same way" for T cells would
need a GTEx donor matrix of the same cell type; GTEx has no T-cell tissue (whole blood is a
mixture; EBV-transformed lymphocytes are B cells), and the z space itself is fibroblast-scaled.
It therefore cannot be built the same way. Per the task: the full space only is reported; no
substitute basis is invented and the fibroblast basis is not applied to T cells. <space> in
every reading is "full measured-gene space".

== Definitions fixed by the session now, because the task and rulings do not specify them
   (disclosed; not shown to the user beforehand) ==

[D1] Sub-screens. A sub-screen is one cell type x one activated condition, with that cell
     type's resting condition as the unaffected state. GSE278572 has two: Treg (HTO_2
     Stimulated_Treg vs HTO_1 Resting_Treg) and Teff (HTO_4 Stimulated_Teff vs HTO_3
     Resting_Teff). GSE314342 has two: Stim8hr vs Rest and Stim48hr vs Rest (the Rest arm is
     shared). Each sub-screen has its own gap, distribution, floor, positive control and
     readings, and nothing is pooled across sub-screens. In the reading wording, <dataset> is
     filled as "<accession> <sub-screen>", e.g. "GSE278572 Teff".
[D2] Which cells give d. d is computed in the activated arm only: the pseudobulk of
     activated-arm cells carrying guides for one target minus the pseudobulk of activated-arm
     control-guide cells. Perturbed cells of the resting arm are not used (they are not
     perturbations of the activated state; SCOPEFIB D6). The gap is activated-arm control-guide
     pseudobulk minus resting-arm control-guide pseudobulk, all control cells of each arm.
[D3] Cell assignment.
     Guides: the deposit's own per-cell guide calls. A cell is assigned only if exactly one guide
     is called (GSE278572: num_features == 1). Cells with zero or >=2 calls are excluded and
     counted. Guides of the same target are pooled (carried D7). Control guides are the guides
     the deposit labels non-targeting (GSE278572: target label "Non-Targeting"; GSE314342: the
     deposit's non-targeting label, read from its guide feature names and named in the findings).
     If a deposit gives guide UMI counts but no per-cell calls, a cell is assigned the guide with
     the highest UMI count when that count is >= 5 and >= 3x the second-highest; otherwise it is
     excluded.
     Conditions: GSE278572, from the four HTO counts in the matrix: the HTO with the highest count
     is assigned when its count is >= 10 and >= 3x the second-highest HTO count (a second-highest
     of 0 satisfies the 3x rule); otherwise the cell is unassigned and excluded, and counted.
     GSE314342: the condition in the deposited per-sample file name (Cell Ranger probe-barcode
     demultiplexing, one h5 per sample).
     No cell QC is added beyond the deposit's filtered barcodes. Only Gene Expression features
     form the space; antibody, HTO and guide features are excluded.
[D4] Pseudobulk and coordinates. Raw counts are summed over the cells of each group. The panel
     of a sub-screen is: every eligible perturbation pseudobulk, the activated control
     pseudobulk and the resting control pseudobulk. The panel is normalised by
     toward_run.tmm_logcpm_quiet (TMM, log2 CPM, prior.count = 2), imported unchanged. Any
     further pseudobulk (floor draws, floor references, positive-control draws) is normalised
     against that panel's cached TMM reference by seng_run.tmm_cache_from_panel and
     seng_run.z_from_counts with mu = 0, sd = 1 and no missing mask, using its log2-CPM output
     only; imported unchanged. ratio = ||x_pert - x_act_ctl||_2 / ||x_act_ctl - x_rest_ctl||_2
     over all Gene Expression features of the deposit (U4).
[D5] Eligibility. A perturbation is eligible in a sub-screen when it has >= 30 assigned cells in
     that sub-screen's activated arm. Below 30: excluded and counted. The distribution, the
     descriptive line and every count are over eligible perturbations.
[D6] Floor (implements U3). For each distinct group size n among eligible perturbations, 1,000
     draws of n cells without replacement from the sub-screen's activated-arm control cells;
     each draw's reference is the pseudobulk of the full activated control pool minus the drawn
     cells; floor ratio = ||x_draw - x_ref|| / ||gap|| with the gap fixed (full pools). P95 is
     the 95th percentile (numpy default linear interpolation) of the 1,000 floor ratios at that
     n. Exact n is used when a sub-screen has <= 200 distinct n; otherwise a grid of 60
     log-spaced integer sizes from 30 to the largest eligible n (unique after rounding), and each
     perturbation is compared with the P95 at the largest grid size <= its n (a smaller group
     has a higher floor, so this direction makes `exceeds` harder). If a perturbation's n
     exceeds half the activated control pool, it is `floor_unavailable`: excluded from n, N and
     the binomial, and counted. A draw pool may be a seeded uniform subsample of at most 60,000
     control cells per arm if the full pool does not fit in memory (GSE314342 only); the
     reference is still the full pool minus the draw. Seed 20261025.
     exceeds: ratio > own-n P95. within_floor (U2): ratio <= own-n P95.
     UNMEASURABLE test (U3): k = eligible, floor-available perturbations that exceed; N_f = eligible,
     floor-available perturbations; p = P(X >= k), X ~ Binomial(N_f, 0.05). UNMEASURABLE fires when
     p >= 0.05. Reported: k, N_f, the chance expectation 0.05*N_f, p.
[D7] Positive control (recovery power; descriptive, gates nothing). A mixing series per
     sub-screen: draws of n cells of which round(f*n) come from resting-arm control cells and the
     rest from activated-arm control cells, f in {0, 0.1, 0.25, 0.5, 1.0}; reference = full
     activated control pool minus the activated cells drawn. 200 draws per (f, n), at n = 30 and
     at n = the sub-screen's median eligible perturbation size (rounded down; and at the grid
     size used for it, where D6's grid applies). Power(f, n) = fraction of draws with ratio > the
     D6 P95 at that n. Also reported: the median ratio per (f, n). f = 0 is a floor calibration
     check (power should be near 0.05). f = 1.0 is a displacement of roughly gap length. Seed
     20261026. Reason: CLAUDE.md requires recovery power beside any null, and a reading of
     UNMEASURABLE or a within_floor count is that null here.
[D8] Descriptive line (task). Over all eligible perturbations of a sub-screen, including
     within_floor ones: "<m> of <M> eligible perturbations have ratio < 2 and therefore could in
     principle have registered an arrival", followed by how many of those <m> are within_floor.
     Labelled descriptive. The distribution (median, Q1, Q3, min, max, fraction >= 2) is over
     all eligible perturbations, and beside it the floor: the median over perturbations of
     their own-n P95, and the floor ratios' median and fraction >= 2 at the sub-screen's median n.
[D9] Readings, per dataset x sub-screen x space, implementing U2:
     RATIO-UNMEASURABLE fires when D6's binomial p >= 0.05.
     Where it does not fire: among eligible, floor-available perturbations that exceed (not
     within_floor), RATIO-BLOCKED fires if n >= 1 have ratio >= 2, with <N> = that set's size;
     RATIO-OPEN fires if n >= 1 have ratio < 2, same <N>. Both may fire. Each fired reading is
     followed by the within_floor count (and floor_unavailable count, if any).
     Where UNMEASURABLE fires, BLOCKED and OPEN are not reported as conclusions there; the counts
     they would have had are reported as descriptive lines, labelled as such.
     If a sub-screen has zero eligible perturbations, no reading fires there and that is said.
[D10] Donors and wells. GSE278572 carries no per-cell donor label in its barcodes or calls
     files as seen; cells are pooled across both donors and all 8 wells. A perturbation's donor
     mix can differ from the control pool's; that difference enters d and not the floor. This is
     disclosed, not corrected. GSE314342 is pooled across donors and lanes under the same rule,
     for consistency, even though its file names carry donors.
[D11] Gene overlap (task hard rule). Count of the deposit's Gene Expression features also in the
     frozen 23,485-gene panel, by version-stripped Ensembl id with symbol fallback (as SCOPEFIB
     D4). A count only. Coordinates here are log2-CPM in each deposit's own gene set, not
     GTEx-fibroblast z in the ruler's gene set, so no number here is on the fibroblast work's
     scale, and a ratio here is not comparable to a SOUTH ratio except as a dimensionless
     ratio of two lengths measured in the same space.
[D12] Cell-count check (CLAUDE.md, shortlist vs cell count). Spearman correlation between ratio
     and n over eligible perturbations, per sub-screen, reported beside the distribution;
     descriptive. Noise lengthens d at small n, so a ratio distribution partly reflects depth.
[D13] Workability (task). GSE278572: 3.6 GB matrix.mtx.gz plus three small files; workable.
     GSE314342: 171.2 GB in 283 per-sample filtered_feature_bc_matrix.h5 files (GSE314342_RAW.tar
     is not downloaded; files are taken individually from the same deposit); 52 GB free disk.
     It is processed by streaming, one file at a time (download, accumulate, delete). It is
     workable only if all three hold: (a) peak disk use <= 40 GB; (b) the wall-clock estimate from
     the measured throughput over its first three files is <= 24 h for the whole deposit; (c)
     the accumulators fit in 15.6 GB RAM. If any fails, the findings say so and which, GSE314342 is
     not analysed, and the task stops after GSE278572. No subset of GSE314342's files is analysed
     in its place (a subset of lanes would change the depth the question is asked at).
[D14] Proceeding. GSE314342 is started only after GSE278572 is complete (all its readings written
     to results/ratio/). A STOP in GSE278572 ends the task.

== STOP conditions, fixed now ==

 [S1] `reproduction_failed` — any gating check in the reproduction plan below fails. Every
      mismatch named.
 [S2] `flag_ambiguous` — the flag is ambiguous or self-contradictory in a way U1-U4 and D1-D14
      do not settle. Write the findings with the STOP key, name the ambiguity, ask.
 [S3] `control_pool_empty` — a sub-screen has zero control-guide cells in its activated or its
      resting arm, so the gap is undefined. Ends the task (D14).
A STOP is never resumed and never restarted with an adjustment.

== Reproduction of stored numbers before anything new (CLAUDE.md) ==

Run first, by src/ratio_run.py reproduce, into results/ratio/reproduction.csv.
 R1 (gating) SOUTH6 S1 O1_Y1 ratio, recomputed with this task's ratio function from stored
    vectors: ||d_f|| from results/south6/stage1.npz `d_f` (1,836 rows) and ||gap|| from the frozen
    anchors results/toward/anchors.npz (c_old - c_young), read not refit. Compared with
    results/south6/stage3_full_table.csv `d_over_v` (S1, O1_Y1) per factor, and its median with
    results/south6/stage4_overall.csv `ratio_median`. Count exact (1,836); numerics within
    relative 1e-9. The S2 settings are not reproduced: their PCA basis is fit at run time in
    SOUTH6 and not stored, so reproducing them would refit an instrument.
 R1b (descriptive, not gating) the task's premise. From the stored stage3_full_table.csv (read,
    not recomputed), per SOUTH6 O1_Y1 setting: min ratio and the count with ratio >= 2.
 R2 (gating) GSE278572 distinct targeting genes in the deposit's guide features = 28 (SCOPEDIS
    n_pert_d1), exact.
 R3 (descriptive, not gating) GSE278572's article states every target exceeded 150 cells in each
    condition (SCOPEDIS units_per_pert). Reported: per target, cells per condition under D3. Not
    gating, because D3's HTO and guide rules are this session's, not the authors'.
 R4 (gating, run only when GSE314342 is reached) distinct targeting genes = 12,748 and guides =
    26,504 in the deposit (SCOPEDIS n_pert_d1), exact.

== Honesty obligations (CLAUDE.md), applied ==

- Chance baseline beside every count: the exceed count carries 0.05*N_f; the fraction >= 2 is
  shown beside the floor's own fraction >= 2 at the median n.
- Recovery power beside any null: D7.
- Circularity, disclosed where results appear: the activated control pool defines zero for d,
  one end of the gap, and the pool every floor draw comes from. The gap's activated end and the
  floor references are therefore not independent of the floor draws. Nothing is refit to avoid it.
- What a null does not cover: these are CD4+ T cells from young or unstated-age donors (GSE314342
  22-34 y; GSE278572 2 donors, age unstated), the resting arm is `returned`, not naive (SCOPEDIS
  U4), and perturbations precede the stimulus (timing `before`). Nothing transfers to aged cells or
  to fibroblasts.
- No conclusion about any named gene; no hit is named. Per-perturbation tables carry the target
  label because they are data, but no target is discussed.

== Process ==

Code: src/ratio_run.py (modes: reproduce, gse278572, gse314342). Imports, unchanged:
toward_run.tmm_logcpm_quiet, toward_run.ANCHORS_PATH, seng_run.tmm_cache_from_panel,
seng_run.z_from_counts, fibro2_common.load_frozen_ruler. Results: results/ratio/. Raw downloads:
data/raw/ratio/ (git-ignored). Deliverables: this file; FINDINGS_RATIO.md with this file
reproduced verbatim at the top; PROGRESS_RATIO.md. One local commit, no push. Not modified: any
FINDINGS_*.md, any SCOPING_*.md, any previous results directory.
```

## Reproduction of stored numbers, before anything new

Run as `python -I src/ratio_run.py reproduce` → `results/ratio/reproduction.csv`, with this task's ratio
function. **6 gating checks, 0 mismatches.** No STOP fired.

| check | expected | observed | match | gating |
|---|---|---|---|---|
| R1 SOUTH6 S1 O1_Y1 factor count | 1836 | 1836 | yes | yes |
| R1 factors matched by name | 1836 | 1836 | yes | yes |
| R1 max relative error, per-factor ratio ‖d_f‖/‖gap‖ vs stored `d_over_v` | ≤1e-9 | 2.6e-15 | yes | yes |
| R1 ‖gap‖ from frozen anchors (`results/toward/anchors.npz`, read not refit) vs stored v_norm | 32.283949049860 | 32.283949049860 | yes | yes |
| R1 median ratio vs `stage4_overall.csv` | 3.820798669451 | 3.820798669451 | yes | yes |
| R2 GSE278572 distinct targeting genes (SCOPEDIS n_pert_d1) | 28 | 28 | yes | yes |
| R1b premise, S1 O1_Y1: min ratio; n ≥ 2 of N | descriptive | min 2.0966; 1836 of 1836 | — | no |
| R1b premise, S2_k20 O1_Y1 | descriptive | min 0.1312; 1 of 1836 | — | no |
| R1b premise, S2_k50 O1_Y1 | descriptive | min 0.1940; 2 of 1836 | — | no |
| R1b premise, S2_k100 O1_Y1 | descriptive | min 0.2384; 2 of 1836 | — | no |
| R3 GSE278572 targets with >150 cells per condition (article claims all) | descriptive | Rest Treg 27, Stim Treg 28, Rest Teff 28, Stim Teff 28 of 28 | — | no |
| R4 GSE314342 genes/guides | — | not run: GSE314342 not reached (D13) | — | — |

R1 reproduces only SOUTH6's S1. Its S2 bases are fit at run time and not stored, so reproducing them
would refit an instrument. R1b reads the S2 minima from the stored table without recomputing them.
R3's one shortfall (a single target with 134 resting-Treg cells) is under this session's D3 hashtag and
guide rules, not the authors'; it is not gating and does not touch the activated arms the statistic uses.

## Datasets

| dataset | download size | status | reason |
|---|---|---|---|
| GSE92872 | 17 MB (CROP-seq digital expression) | **not analysed** | U1: the user ruled "Skip GSE92872". SCOPEDIS recorded its CROP-seq arm `separate_experiment`; only its 20-gene bulk arm was `in_run`. |
| GSE278572 | 3.86 GB (`matrix.mtx.gz` 3,857,480,145 bytes) + 3.4 MB of metadata | **analysed** | gzip integrity checked; 1,027,798,128 entries read, equal to the header nnz |
| GSE314342 | 171.15 GB in 284 per-sample h5 files | **not analysed: not workable (D13 b)** | 3 files (1.92 GB) downloaded together in 1,885 s = 1.02 MB/s → ~46.7 h for the deposit's download alone, against 24 h. (a) disk passed; (c) RAM not assessed, since (b) decides. Timed files deleted unopened. |

No substitute dataset was used, and no subset of GSE314342 was analysed.

## GSE278572 inventory

Deposit: 249,799 filtered barcodes over 8 wells; 36,601 Gene Expression features, 141 antibody
features including four hashtags, and 65 guides (28 targets × 2, plus 9 non-targeting).

| step (flag D3) | cells |
|---|---|
| exactly one guide called | 140,641 |
| zero guide calls (excluded) | 58,179 |
| ≥2 guide calls (excluded) | 50,979 |
| hashtag assigned (top ≥ 10 and ≥ 3× the second) | 214,773 |
| hashtag unassigned (excluded) | 35,026 |
| assigned: Resting_Treg / Stimulated_Treg / Resting_Teff / Stimulated_Teff | 40,526 / 77,502 / 43,176 / 53,569 |
| control-guide cells, single guide and hashtag-assigned, all four arms | 19,001 |

Cells per target per arm: `results/ratio/gse278572_cells_per_target.csv`. No per-cell donor label
exists in the deposit as read, so the two donors are pooled (D10). A perturbation's donor mix can differ
from the control pool's; that difference enters d and not the floor.

**Gene overlap with the frozen 23,485-gene panel:** 20,459 of the deposit's 36,601 Gene Expression
features (20,443 by version-stripped Ensembl id, 16 by symbol fallback). **These coordinates are not
comparable to the fibroblast work.** Here the space is log2-CPM over the deposit's own 36,601 genes,
with no per-gene scaling (U4). SOUTH6's S1 is GTEx-fibroblast z over the ruler's genes. Only the ratio,
a dimensionless quotient of two lengths each measured within one space, has the same meaning in both.

**Compressed space: not built.** Paper 2's compressed space is a PCA on 308 GTEx cultured-fibroblast
donors aged 40–59 in the fibroblast z space, and GTEx has no T-cell tissue to repeat that on. Per the
task, no substitute basis was invented and the fibroblast basis was not applied to T cells. Every number
below is in the full measured-gene space.

## GSE278572 — the ratio, per sub-screen (full measured-gene space)

d is computed in the activated arm only (D2). The gap is the activated-arm control pseudobulk minus the
resting-arm control pseudobulk. The panel is normalised by TMM and log2-CPM with prior 2. Every extra
pseudobulk (floor and positive-control draws) is normalised against the panel's cached reference. A
self-check confirmed the cached and panel normalisations agree to 0.0 on the control pseudobulk.

| | Treg (Stim vs Rest) | Teff (Stim vs Rest) |
|---|---|---|
| perturbations / eligible (≥30 cells) / excluded | 28 / 28 / 0 | 28 / 28 / 0 |
| activated / resting control-guide cells | 5,805 / 4,330 | 4,489 / 4,377 |
| perturbation cells, activated arm: min / median / max | 472 / 1,145 / 2,041 | 369 / 895 / 1,507 |
| ‖gap‖ (log2-CPM units) | 144.56 | 151.83 |
| **ratio median [Q1, Q3]** | **0.481 [0.439, 0.559]** | **0.451 [0.415, 0.497]** |
| ratio min / max | 0.397 / 0.763 | 0.387 / 0.628 |
| **fraction with ratio ≥ 2** | **0.000** | **0.000** |
| *floor beside it:* P50 / P95 at the median size | 0.446 / 0.452 | 0.431 / 0.437 |
| floor fraction ≥ 2 at the median size | 0.000 | 0.000 |
| floor P95, median over perturbations (own size) | 0.452 | 0.438 |
| floor P50 / P95 at 30 cells | 1.736 / 1.815 | 1.517 / 1.578 |
| exceed floor P95 (k) of N_f | **19 of 28** | **17 of 28** |
| chance expectation (5%) | 1.40 | 1.40 |
| one-sided binomial p | 8.5e-19 | 9.6e-16 |
| within_floor / floor_unavailable | 9 / 0 | 11 / 0 |
| exceeders with ratio ≥ 2 / < 2 | 0 / 19 | 0 / 17 |
| Spearman(ratio, cells) | −0.56 | −0.45 |

All floors used exact sizes: 30 sizes per sub-screen (28 perturbation sizes plus 30 cells and the
median), 1,000 draws each. No perturbation exceeded half the control pool, so floor_unavailable is 0.
Per-size floors: `gse278572_{treg,teff}_floor.csv`; every draw: `_floor_draws.npz`.

**Geometry, not biology.** A ratio ≥ 2 would mean no displacement of that length could end closer to
the resting state, whatever its direction. None here is ≥ 2. A ratio < 2 means a displacement of that
length *could* end closer, if its direction were right. This task measured no direction, so it says
nothing about whether any perturbation moved toward rest.

**Borderline case, reported.** One Treg perturbation (472 cells) exceeds its floor by less than 0.001
(ratio 0.6411 against P95 0.6409). It is counted as an exceeder under the rule as written.

### Readings, per sub-screen (pre-registered wording; `results/ratio/readings.txt`)

```
== GSE278572 Treg / full measured-gene space ==
CONCLUSION  RATIO-OPEN: in GSE278572 Treg, 19 of 19 perturbations have ratio < 2 in full measured-gene space; arrival was geometrically attainable for them. [within_floor = 9]
            (descriptive) 28 of 28 eligible perturbations have ratio < 2 and therefore could in principle have registered an arrival; 9 of those are within_floor.
            exceed k = 19 of N_f = 28; chance expectation 1.40; one-sided binomial p = 8.5e-19
== GSE278572 Teff / full measured-gene space ==
CONCLUSION  RATIO-OPEN: in GSE278572 Teff, 17 of 17 perturbations have ratio < 2 in full measured-gene space; arrival was geometrically attainable for them. [within_floor = 11]
            (descriptive) 28 of 28 eligible perturbations have ratio < 2 and therefore could in principle have registered an arrival; 11 of those are within_floor.
            exceed k = 17 of N_f = 28; chance expectation 1.40; one-sided binomial p = 9.6e-16
```

RATIO-BLOCKED: did not fire in either sub-screen (0 exceeders with ratio ≥ 2).
RATIO-UNMEASURABLE: did not fire in either sub-screen (binomial p < 0.05).
`python -I src/ratio_run.py verify` recomputes k, N_f, p and the ≥2/<2 split from the tracked
per-perturbation CSVs alone: 0 mismatches.

### Recovery power (positive control, D7; descriptive)

Draws of n cells, a fraction f of them from resting-arm control cells and the rest from activated-arm
control cells. Power is the share of 200 draws above the floor's P95 at that n.

| sub-screen | n | f = 0 | 0.1 | 0.25 | 0.5 | 1.0 |
|---|---|---|---|---|---|---|
| Treg | 30 | 0.065 (1.74) | 0.215 (1.77) | 0.70 (1.84) | 1.00 (1.95) | 1.00 (2.59) |
| Treg | 1,145 | 0.030 (0.45) | 0.955 (0.46) | 1.00 (0.49) | 1.00 (0.58) | 1.00 (1.15) |
| Teff | 30 | 0.070 (1.52) | 0.245 (1.55) | 0.76 (1.61) | 1.00 (1.73) | 1.00 (2.35) |
| Teff | 895 | 0.055 (0.43) | 0.965 (0.45) | 1.00 (0.48) | 1.00 (0.57) | 1.00 (1.15) |

Cells show power, with the median ratio in parentheses. At f = 0 the power sits near the nominal 0.05,
so the floor is calibrated. At the median depth the statistic separates a 10% mixture of resting cells
from noise in ≥ 95% of draws. Two readings of this table matter:

- At 30 cells, a displacement of the full gap (f = 1) has a median measured ratio of 2.59 (Treg) and
  2.35 (Teff), **above 2**. Noise lengthens d enough that a gap-length move toward rest would itself
  read as "blocked" at this depth. The ratio-2 line is therefore depth-sensitive, not a fixed property
  of the perturbation.
- At the median depth, f = 1 reads 1.15, not 1.0, for the same reason.

The within_floor counts (9 and 11) sit next to a power of ≥ 0.955 at f = 0.1. Those perturbations are
therefore within the noise *at a depth where a displacement of about a tenth of the gap would usually
have been seen*. This statement leans on the mixing control's direction (toward rest) and length; it does
not establish power for displacements in other directions.

### Cell-count check (CLAUDE.md)

Spearman(ratio, cells) is −0.56 (Treg) and −0.45 (Teff): smaller groups have longer measured d, as noise
predicts. Exceeding is not just "the smallest groups": the median size is 1,039 for exceeders against
1,300 for within_floor perturbations (Treg), and 886 against 1,020 (Teff), so exceeders are somewhat
smaller. Against a size-matched floor that is not by construction an advantage, but it is reported.

## Diagnostic, NOT pre-registered (gates nothing; `python -I src/ratio_run.py diagnostic`)

Added after the readings were written, because the floor sat so close to the perturbation ratios. If
noise adds roughly orthogonally to an effect, then ‖d_measured‖² ≈ ‖d_effect‖² + ‖noise‖². So
sqrt(max(ratio² − floor_P50(n)², 0)) is a rough noise-corrected length, in gap units. That
orthogonality is an assumption and is not tested.

| sub-screen | ratio / own P95: median, max | noise-corrected length: median [Q1, Q3], max |
|---|---|---|
| Treg | 1.03, 1.73 | 0.146 [0.000, 0.302], 0.626 |
| Teff | 1.02, 1.48 | 0.128 [0.000, 0.210], 0.414 |

Per perturbation: `gse278572_{treg,teff}_diagnostic.csv`; summary: `diagnostic_summary.csv`.

## Deviations and decisions not in the flag (disclosed)

1. **Four points were put to the user before the flag was written**, and their rulings are reproduced
   verbatim in it (U1–U4). U1 removed GSE92872 from the task. That is the user's ruling, not a
   substitution.
2. **Pre-flag reads, beyond repository files.** Before the flag existed, three small GSE278572 metadata
   files were downloaded and their *structure* read (feature types and names, barcode suffixes and
   per-suffix totals, the calls file's header and three rows). FTP listings and file sizes were also
   fetched. SOUTH6's stored Stage 4 medians had been seen in `FINDINGS_SOUTH6.md`. All of this is
   listed at the top of the flag. No per-cell assignment was tabulated before the flag.
3. **The first GSE278572 analysis attempt crashed** (`gse278572_stdout_attempt1_crashed.txt`). Windows
   multiprocessing could not pickle the ~500 MB control pools into the worker processes (OSError 22).
   This was a software failure, not a pre-registered STOP. Logged before the crash: the streaming
   counts, the Treg sub-screen's cell counts and its ‖gap‖ (144.5565). No per-perturbation ratio had
   been written or seen. The fix only changed how pools reach the workers (via a temporary file), and
   the worker count went from 12 to 8 to fit RAM. Draw indices are generated in the main process from
   the fixed seed, so the worker count cannot change a number. The rerun read the cached stream output
   (`data/raw/ratio/gse278572/_streamed.npz`, git-ignored). It reproduced the logged ‖gap‖ exactly, and
   the synthetic smoke test gave identical output before and after the fix. Nothing in the statistic,
   seeds, rules or thresholds was changed. I judged this a crash to rerun, not a STOP to honour.
   CLAUDE.md's "never resume" rule is about STOPs, and nothing had been computed that could have
   steered the rerun. Disclosed so a reader can judge it.
4. **The first GSE314342 workability measurement was invalid**
   (`gse314342_workability_attempt1_badurl.json`, kept, not used). The per-sample files are served
   from each GSM's own directory; the series `suppl/` directory holds them only inside `RAW.tar`. The
   first attempt fetched 2,970 bytes of error pages and reported a meaningless ~13,850 h. The code was
   fixed to use the GSM paths and to refuse any conclusion unless every byte arrives. The second
   measurement received all 1,917,820,392 bytes and is the one D13 is read from.
5. **File count in the flag was wrong by one.** D13 says 283 per-sample files; GSE314342's
   `filelist.txt` lists 284 (my miscount when writing the flag). The total size, 171.15 GB, is as
   stated, and the decision does not depend on the count.
6. **Modes added beyond the flag's list**: `verify` (re-derives the gated quantities from tracked CSVs),
   `diagnostic` (above, not pre-registered) and `gse314342_workability` (implements D13).
7. **Floor sizes for the positive control.** D7 compares mixing draws with "the D6 P95 at that n". For
   30 cells and the median size, neither of which is necessarily a perturbation's size, the floor was
   computed at those sizes as well. In Treg the median size is 1,145, between two perturbations' sizes.
8. **One guard added in code**: if a mixing draw needed more resting cells than the pool holds, the code
   would STOP `flag_ambiguous`. It did not fire.
9. **Imports, unchanged:** `toward_run.tmm_logcpm_quiet`, `toward_run.ANCHORS_PATH`,
   `seng_run.tmm_cache_from_panel`, `seng_run.z_from_counts` (with mu = 0, sd = 1, using its log2-CPM
   output only), `fibro2_common.load_frozen_ruler`. The panel's TMM reference sample is chosen by the
   imported upper-quartile rule.
10. **Where things were written.** Raw downloads and the 994 MB stream cache are in `data/raw/ratio/`
    (git-ignored). Temporary pool files were written there during the floor and deleted. The three
    timed GSE314342 files were deleted unopened. Results are in `results/ratio/` only; no previous
    results directory, `FINDINGS_*.md` or `SCOPING_*.md` was modified.
11. **Untrusted inputs.** Downloaded files sit in their own directories under `data/raw/ratio/`. All
    Python was run with `-I` from the repository root, from scripts outside those directories.

## Honesty obligations

- **Chance baseline beside counts:** the exceed counts (19, 17) sit beside 1.4 expected by chance. The
  fraction ≥ 2 (0.000) sits beside the floor's own fraction ≥ 2 at the median size (also 0.000).
- **Recovery power:** in the positive-control table. Power is ≥ 0.955 at f = 0.1 at the median depth and
  0.215–0.245 at f = 0.1 at 30 cells.
- **Circularity (disclosed, not refitted):** the activated control pool defines zero for every d, one end
  of the gap, and the pool every floor draw comes from. Floor references are that pool minus each draw.
  The gap's activated end and the floor are therefore not independent.
- **What a null would not cover / what this does not transfer to:** these are primary CD4+ Tregs and Teffs
  from two donors of unstated age. The resting arm is `returned` (activated at isolation, rested), not
  naive (SCOPEDIS U4). Perturbations precede restimulation (timing `before`). The gap is a restimulation
  gap, not an aging or disease gap. Nothing here transfers to aged cells, to fibroblasts, or to the
  Southard screen's coordinates.
- **No gene is named or discussed.** Per-perturbation CSVs carry target labels because they are data.

## What this does not cover

- **Generality.** One experiment, two sub-screens of 28 perturbations each. The task's question ("a
  general property of perturbation screens") is not answered by it. GSE314342, the genome-scale screen,
  was not workable here; GSE92872 was excluded by ruling.
- **Depth.** All 56 perturbation groups have 369–2,041 cells, far above the 30-cell minimum. The floor
  shows the ratio is strongly depth-dependent: at 30 cells, noise alone reads 1.5–1.8 and a full-gap
  move reads 2.4–2.6. A screen at Southard-like depth (median 106 cells per target there) would sit much
  higher on this scale. Whether SOUTH6's S1 ratio ≥ 2 was itself depth-driven is **not** known: SOUTH6
  never set its ratio beside a size-matched control floor, and this task did not compute one for it.
- **Gap size.** The denominator here is a large restimulation contrast. A smaller gap, such as an aging
  gap, would raise every ratio in proportion. "Short enough" is relative to the gap, and the gap differs
  between the fibroblast-aging work and this screen by an unknown factor in incomparable units.
- **Direction.** Nothing here measures whether any perturbation points toward the resting state.

## Next

None within this task (one task per session). Two follow-ups would settle what this leaves open,
neither chosen here:
(1) the same size-matched control floor for SOUTH6's S1 ratio, from its stored non-targeting cells, to
learn whether "every factor ≥ 2" there was depth;
(2) GSE314342 on a connection or host where its 171 GB is workable.
