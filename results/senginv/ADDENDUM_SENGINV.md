# SENGINV addendum - counting and listing only

Counts and lists. No distance, no fitted axis, no statistic that needs a
null, no pre-registration, nothing committed. Item 1 is in the companion
section below this one. Ambiguities are recorded, not resolved.

Inputs: the same untracked `data/raw/seng/exp3_merged_not_normalized.h5ad`,
md5 `bfccaa94a916c22e4b0793913a04128b`, verified in the SENGINV run. `.X` was
streamed in 2000-row blocks through its CSR `indptr`; the matrix was never
loaded. Script: `src/senginv_addendum.py`.

Definitions used throughout, stated once:

- **barcode group**: the suffix on the `.obs` index string, one of `-1`,
  `_A`, `_I`. The file never says these are libraries; the task calls them
  libraries and that reading is adopted as the task's, not the file's.
- **NT-positive**: at least one of `NT_1`..`NT_6` is > 0.
- **targeting-guide reads**: the sum over the 798 `<TARGET>_<CRA|CRI><k>`
  columns. `NT_*` columns are excluded from that sum by construction, which
  is the definitional choice item 3 turns on.
- **WT cell**: `WT` > 0.
- **gene detected in a group**: at least one cell of that group has a
  nonzero entry for it in `.X`.

Group sizes, fixed here and referred to below:

| group | n cells |
|---|---|
| barcode `-1` | 6809 |
| barcode `_A` | 14333 |
| barcode `_I` | 15372 |
| `PD14` | 581 |
| `PD26` | 1880 |
| `PD32` | 4348 |
| all | 36514 |

## 2. Control cells by library (barcode group)

| barcode group | n cells | NT-positive | WT>0 | NT-positive and WT>0 | neither |
|---|---|---|---|---|---|
| `-1` | 6809 | 0 | 6809 | 0 | 0 |
| `_A` | 14333 | 1156 | 0 | 0 | 13177 |
| `_I` | 15372 | 3261 | 0 | 0 | 12111 |
| **total** | 36514 | 4417 | 6809 | 0 | 25288 |

Per `NT_*` column, by barcode group:

| column | `-1` | `_A` | `_I` | total |
|---|---|---|---|---|
| `NT_1` | 0 | 231 | 910 | 1141 |
| `NT_2` | 0 | 184 | 679 | 863 |
| `NT_3` | 0 | 188 | 828 | 1016 |
| `NT_4` | 0 | 256 | 943 | 1199 |
| `NT_5` | 0 | 491 | 1299 | 1790 |
| `NT_6` | 0 | 60 | 124 | 184 |
| **any `NT_*`** | 0 | 1156 | 3261 | 4417 |

NT UMI totals, by barcode group (sum over the six `NT_*` columns):

| barcode group | total NT UMI | NT UMI in NT-positive cells: min / median / max |
|---|---|---|
| `-1` | 0 | no NT-positive cells |
| `_A` | 13754 | 1 / 4 / 212 |
| `_I` | 40023 | 1 / 5 / 130 |

WT and the PD arms against barcode group (WT is already known to coincide
with the union of the PD arms; repeated here because item 2 asks for it):

| barcode group | `WT`>0 | `PD14` | `PD26` | `PD32` |
|---|---|---|---|---|
| `-1` | 6809 | 581 | 1880 | 4348 |
| `_A` | 0 | 0 | 0 | 0 |
| `_I` | 0 | 0 | 0 | 0 |

## 3. NT-positive cells with zero targeting-guide reads

- NT-positive cells: 4417
- of those, **1182 have zero reads in all 798 targeting-guide columns**
  and **3235 have at least one targeting-guide read**.

| barcode group | NT-positive | zero targeting-guide reads | at least one | % zero |
|---|---|---|---|---|
| `-1` | 0 | 0 | 0 | n/a |
| `_A` | 1156 | 439 | 717 | 38.0% |
| `_I` | 3261 | 743 | 2518 | 22.8% |
| **total** | 4417 | 1182 | 3235 | 26.8% |

For the NT-positive cells that do carry targeting-guide reads, how many
targeting-guide columns and how much targeting UMI they carry:

| barcode group | n such cells | n targeting cols: min/median/max | targeting UMI: min/median/max | NT UMI: median |
|---|---|---|---|---|
| `-1` | 0 | - | - | - |
| `_A` | 717 | 1/2/8 | 1/9/279 | 2 |
| `_I` | 2518 | 1/3/46 | 1/19/195 | 3 |

**Recorded as an ambiguity.** Whether an NT-positive cell carrying
targeting-guide reads is an NT control, a perturbed cell, a doublet or
ambient contamination is not something the object states, and no threshold
on NT UMI or targeting UMI is recorded in it. The split above uses >0 on
both sides because that is the only cut the file supports without a choice
being imported. A different cut would move these numbers.

## 4. Read depth by group

From the stored `.obs` QC columns `n_counts` and `n_genes`, not recomputed.
Percentiles are linear-interpolation `numpy.percentile`.

### `n_counts`

| group | n | mean | p0 | p1 | p5 | p10 | p25 | p50 | p75 | p90 | p95 | p99 | p100 |
|---|---|---|---|---|---|---|---|---|---|---|---|---|---|
| barcode `-1` | 6809 | 15422 | 500 | 744 | 1510 | 2214 | 5471 | 11872 | 22520 | 32712 | 40150 | 59528 | 89974 |
| barcode `_A` | 14333 | 19978 | 500 | 1145 | 4589 | 7166 | 13702 | 19611 | 25634 | 32061 | 36311 | 47334 | 90340 |
| barcode `_I` | 15372 | 21577 | 502 | 792 | 3489 | 7229 | 15630 | 21867 | 27690 | 33529 | 37522 | 47607 | 79652 |
| `PD14` | 581 | 20721 | 547 | 725 | 2411 | 5287 | 10335 | 19853 | 28368 | 37143 | 43339 | 60571 | 88657 |
| `PD26` | 1880 | 17439 | 590 | 1510 | 2909 | 4336 | 8390 | 15533 | 23305 | 31924 | 40161 | 58481 | 81481 |
| `PD32` | 4348 | 13841 | 500 | 693 | 1307 | 1799 | 3760 | 9140 | 20494 | 32130 | 39733 | 59377 | 89974 |
| all cells | 36514 | 19801 | 500 | 843 | 2676 | 5193 | 12335 | 19777 | 26324 | 32870 | 37430 | 50209 | 90340 |

### `n_genes`

| group | n | mean | p0 | p1 | p5 | p10 | p25 | p50 | p75 | p90 | p95 | p99 | p100 |
|---|---|---|---|---|---|---|---|---|---|---|---|---|---|
| barcode `-1` | 6809 | 3317 | 33 | 423 | 901 | 1191 | 1999 | 3134 | 4652 | 5584 | 6103 | 7022 | 8430 |
| barcode `_A` | 14333 | 4162 | 84 | 605 | 1592 | 2306 | 3572 | 4342 | 4973 | 5560 | 5950 | 6712 | 8602 |
| barcode `_I` | 15372 | 4294 | 99 | 453 | 1430 | 2312 | 3845 | 4533 | 5101 | 5617 | 5920 | 6559 | 7994 |
| `PD14` | 581 | 4151 | 41 | 330 | 1242 | 2018 | 2955 | 4325 | 5290 | 5975 | 6357 | 7354 | 8049 |
| `PD26` | 1880 | 3756 | 67 | 730 | 1424 | 1844 | 2600 | 3766 | 4810 | 5665 | 6196 | 7072 | 8274 |
| `PD32` | 4348 | 3016 | 33 | 400 | 811 | 1024 | 1683 | 2628 | 4317 | 5476 | 6007 | 6898 | 8430 |
| all cells | 36514 | 4060 | 33 | 483 | 1262 | 1854 | 3342 | 4336 | 5006 | 5592 | 5960 | 6713 | 8602 |

The two QC columns are checked against `.X` itself, exactly and for every
row, in the same streaming pass as section 5. The result is reported after
section 5 because it needs that pass.

## 5. Genes detected per group, and the three-way barcode overlap

- streamed 19 row blocks of at most 2000 rows; peak slice read was the
  `indices` and `data` range of one block, never the whole 148,239,518-element
  arrays.

### Do the stored QC columns describe this `.X`?

Exact, every row, from the stream above.

| comparison | n rows equal | n rows differing | min diff | max diff |
|---|---|---|---|---|
| row sum of `.X` vs stored `n_counts` | 29401 | 7113 | -22 | 0 |
| row nnz of `.X` vs stored `n_genes` | 29401 | 7113 | -5 | 0 |

| barcode group | n rows where `n_counts` differs | where `n_genes` differs |
|---|---|---|
| `-1` | 0 | 0 |
| `_A` | 2349 | 2349 |
| `_I` | 4764 | 4764 |

The shortfall against candidate sums the file does hold, as a count of
exact matches. This is arithmetic, not an explanation:

| candidate | rows where shortfall equals it (of 7113 mismatched) |
|---|---|
| NT UMI total | 282 |
| targeting-guide UMI total | 347 |
| NT + targeting UMI | 303 |

| candidate | rows where the `n_genes` shortfall equals it |
|---|---|
| n nonzero `NT_*` columns | 648 |
| n nonzero targeting-guide columns | 1067 |
| n nonzero guide columns of either kind | 914 |

| shortfall in `n_counts` | n rows |   | shortfall in `n_genes` | n rows |
|---|---|---|---|---|
| 1 | 5277 |   | 1 | 6202 |
| 2 | 1316 |   | 2 | 841 |
| 3 | 352 |   | 3 | 64 |
| 4 | 96 |   | 4 | 5 |
| 5 | 30 |   | 5 | 1 |
| 6 | 15 |   |  |  |
| 7 | 7 |   |  |  |
| 8 | 8 |   |  |  |
| 9 | 2 |   |  |  |
| 11 | 2 |   |  |  |
| 12 | 1 |   |  |  |
| 13 | 1 |   |  |  |
| 16 | 1 |   |  |  |
| 17 | 2 |   |  |  |
| 18 | 1 |   |  |  |
| 21 | 1 |   |  |  |
| 22 | 1 |   |  |  |

**Recorded as an ambiguity.** Where the two disagree the stored value
exceeds the matrix on 7113 rows and falls short on 0, by at most 22 counts
and 5 genes. The object does not record what matrix `n_counts` and
`n_genes` were computed on, so whether the difference is the guide UMI
columns, a feature set that was later trimmed, or something else cannot be
read off the file. Section 4 uses the stored columns as the task asked; this
check is reported so that choice is visible.

Genes detected, with the cell count beside every number because group size
sets the ceiling on how many genes can be detected at all:

| group | n cells | genes in >=1 cell | >=3 cells | >=10 cells | >=1% of cells | total UMI-bearing entries |
|---|---|---|---|---|---|---|
| barcode -1 | 6809 | 22353 | 19033 | 16270 | 12814 | 22585389 |
| barcode _A | 14333 | 24635 | 21289 | 18235 | 13058 | 59658027 |
| barcode _I | 15372 | 25489 | 22126 | 18940 | 13373 | 65996102 |
| PD14 | 581 | 17221 | 14349 | 12283 | 13075 | 2411574 |
| PD26 | 1880 | 19526 | 16431 | 14010 | 12965 | 7060392 |
| PD32 | 4348 | 21224 | 17934 | 15250 | 12692 | 13113423 |
| barcode -1 matched n=6809 | 6809 | 22353 | 19033 | 16270 | 12814 | 22585389 |
| barcode _A matched n=6809 | 6809 | 22979 | 19570 | 16620 | 13075 | 28345082 |
| barcode _I matched n=6809 | 6809 | 23631 | 20189 | 17112 | 13377 | 29258448 |
| PD14 matched n=581 | 581 | 17221 | 14349 | 12283 | 13075 | 2411574 |
| PD26 matched n=581 | 581 | 17290 | 14339 | 12239 | 13073 | 2180580 |
| PD32 matched n=581 | 581 | 17020 | 14080 | 11924 | 12794 | 1770110 |

- written: `results/senginv/addendum_genes_detected.csv`
- the four `matched n=` rows are a size-matched diagnostic, not a requested
  result: each is one subsample without replacement at seed 20260914, drawn to
  the smallest group of its family (6809 cells for the barcode groups, 581 for
  the PD arms). One draw is not a null distribution and no interval is given.

### Three-way overlap between the barcode groups

Detection is 'nonzero in at least one cell of the group'.

| region | n genes |
|---|---|
| in all three (`-1` and `_A` and `_I`) | 21260 |
| `-1` and `_A` only | 219 |
| `-1` and `_I` only | 569 |
| `_A` and `_I` only | 2110 |
| `-1` only | 305 |
| `_A` only | 1046 |
| `_I` only | 1550 |
| in none of the three | 9542 |
| **all `.var` genes** | 36601 |

| | `-1` | `_A` | `_I` |
|---|---|---|---|
| `-1` | 22353 | 21479 | 21829 |
| `_A` | 21479 | 24635 | 23370 |
| `_I` | 21829 | 23370 | 25489 |

Jaccard between barcode groups, on the detected-gene sets:

| pair | intersection | union | Jaccard |
|---|---|---|---|
| `-1` vs `_A` | 21479 | 25509 | 0.8420 |
| `-1` vs `_I` | 21829 | 26013 | 0.8392 |
| `_A` vs `_I` | 23370 | 26754 | 0.8735 |

Same overlap on the size-matched subsamples (diagnostic, one draw, seed 20260914):

| region | n genes |
|---|---|
| in all three | 20420 |
| `-1` only | 663 |
| `_A` only | 985 |
| `_I` only | 1247 |
| in none | 10882 |

### Genes detected per PD arm, and their overlap

| region | n genes |
|---|---|
| in all three PD arms | 16468 |
| `PD14` only | 247 |
| `PD26` only | 770 |
| `PD32` only | 2186 |
| in no PD arm | 14248 |

| | `PD14` | `PD26` | `PD32` |
|---|---|---|---|
| `PD14` | 17221 | 16580 | 16862 |
| `PD26` | 16580 | 19526 | 18644 |
| `PD32` | 16862 | 18644 | 21224 |

**Recorded as an ambiguity, and it is the one that governs this whole
section.** The three barcode groups hold 6809, 14333 and 15372 cells and the three
PD arms hold 581, 1880 and 4348. Detection at a >=1-cell threshold is a
saturating function of cell count and of depth per cell, so a larger group
detects more genes whether or not anything biological differs. The matched
rows above hold cell count fixed but not depth per cell, and section 4 shows
the depth distributions differ between groups. No comparison here is
corrected for either, and none of these counts is a result about the groups.

Genes detected in no cell at all: 9542 of 36601 `.var` columns.


---

# SENGINV addendum, item 1 — audit of the earlier SENG work

Listing only. Nothing was rerun. Every number below is quoted from a file
already in the tree, with the file named, and is **not** reproduced against
this task's code — that is the chief limitation of this section and the
reason none of it is a check on anything.

## 1.1 Git entries

| commit | what it carried |
|---|---|
| `d79ecf3` "Initial public release: code and results for preprint v12" | **every** SENG artefact below. `git log` on all SENG paths returns this one commit. |

There is no per-task SENG commit. The SENG, SENG_WHY, SURVEY and PAPER_FIGS
work was squashed into the initial public-release commit, so the repository's
own history cannot order those tasks or show what each one added. Dates below
are file mtimes, which survived the squash, not commit dates.

**Recorded as an ambiguity.** CLAUDE.md requires one commit per task, each
carrying its findings file, `src/` entry and results directory. For the SENG
family that requirement is not met in the committed history, and this session
is not fixing it (no commits in this task).

## 1.2 Tasks, in file-mtime order

### SURVEY — remote peek, before any download (2025-09-22 13:53 to 14:05)

| file | kind | what it is |
|---|---|---|
| `src/survey_coverage.py` | script | the survey's screen-dataset triage (not SENG-specific) |
| `results/survey/opened/esearch_seng.json` | result | NCBI esearch for the paper |
| `results/survey/opened/sengstack_head.bin` (4,096 B) | result | HDF5 superblock, read by HTTP range |
| `results/survey/opened/sengstack_obs.bin` (64,512 B) | result | `obs` group bytes |
| `results/survey/opened/sengstack_obs_heap.bin` (22,529 B) | result | `obs` local heap — where the 200 target names were read from |
| `results/survey/opened/sengstack_obs_leaf.bin` (16,192 B) | result | `obs` B-tree leaf |
| `results/survey/opened/sengstack_var*.bin` (4 files) | result | `var` group bytes; gene symbols were **not** extracted |
| `results/survey/opened/sengstack_xcont.bin` (2,048 B) | result | `X` continuation block, where the shape attribute was read |
| `SURVEY_TF.md` | findings | the triage table and its reasoning |

**Pre-registration: none.** `SURVEY_TF.md` is a scoping document and there is
no `PREREG` flag under `results/survey/`.

**Numbers computed on this h5ad** (all by HTTP range request, without
downloading it):

- `X` shape attribute 36,514 × 36,601
- 200 target symbols in the `obs` name heap, each with a CRA and a CRI label,
  plus `NT_1` through `NT_6`
- 179 of the 200 symbols on Lambert v1.01; the other 21 include ATG3, BRCA1,
  EZH2, HMGB1, PARP1, SIRT1
- `PD14`, `PD26`, `PD32` seen among the obs labels
- ruler coverage: **not computed** (gene symbols not extracted)

### SENG — the gated task (2025-09-22 15:09 to 17:25)

| file | kind | what it is |
|---|---|---|
| `results/seng/PREREG.flag` (7,578 B, 15:09) | flag | the pre-registration, written before any number existed |
| `src/seng_run.py` (2,464 lines) | script | download, Stage 0 inventory, split, axes, positive control, proliferation quartiles, the screen |
| `results/seng/download.json` | result | url, size, md5 of the h5ad |
| `results/seng/h5_tree.txt` (50,666 B) | result | full HDF5 tree dump |
| `results/seng/inventory.json` + `inventory_{guides,nt,other,perturbations,wt_pd}.csv` | result | Stage 0 descriptive inventory |
| `results/seng/split.json`, `split_table.csv`, `split_idx.npz` | result | early/late cohort split and its A/B halves |
| `results/seng/axes.json`, `axes.npz` | result | the fitted passage axis `v` and its cosines |
| `results/seng/posctrl.csv`, `posctrl_summary.json` | result | positive control on the split |
| `results/seng/prolif_quartiles.json`, `prolif_quartile_idx.npz` | result | NT proliferation quartiles and the Tirosh gene lists found |
| `results/seng/seurat_cc.genes_tirosh2016.tsv` (mtime 09-27) | result | the cell-cycle gene list; **mtime is five days after the rest of the directory** |
| `results/seng/screen.csv` (228,652 B, 46 columns), `screen_summary.json` | result | the per-perturbation screen, the task's primary table |
| `results/seng/named_hits.csv` | result | the four named hits |
| `results/seng/manifest.json`, `reading.json`, `run_report.txt` | result | failure record, fired reading, run log |
| `FINDINGS_SENG.md` (694 lines) | findings | the pre-registration reproduced verbatim plus the readings |
| `PROGRESS_SENG.md` (22 lines) | progress | status DONE, 1 failure recorded, elapsed 4,973.8 s |

**Pre-registration: yes.** `results/seng/PREREG.flag`, timestamped 15:09,
before the first number. It fixes the seeds (20260914 splits/nulls/permutations,
20260918 bootstrap), `n_null=n_perm=n_boot=200`, the STOP keys
(`no_raw_counts`, `no_wt_reference`), that uncalibrated R² is never a gate, and
that Southard/Hs27 and GSE325735 are out of scope.

**Numbers computed on this h5ad:**

*Stage 0 inventory* — file 1,426,871,504 B, md5 `bfccaa94a916c22e4b0793913a04128b`;
`X` is CSR float32, all 148,239,518 stored values scanned, 0 fractional, 0
negative, 0 non-finite, min 1.0, **max 2960.0**, verdict raw integer counts
true; `n_WT`=6,809; PD14 581 / PD26 1,880 / PD32 4,348; perturbed-cell PD
**not determinable** (perturbed PD values empty, blank PD n=24,610); ties with a
positive maximum 2,680 (no winner given); TF-guide calls disagreeing with the
CRA/CRI indicator 270 (not assigned); unique-argmax NT cells CRA 666, CRI 1,479;
NT per modality in `inventory_nt.csv` (CRA 123/47/42/125/326/3, CRI
393/70/293/161/561/1 for NT_1..NT_6); unassigned 2,950; ruler overlap
23,485 panel → 16,781 present, 6,704 missing, 0 by Ensembl, 16,781 by symbol.

*Split* — early PD14 n=581 (A 290, B 291), late PD32 n=4,348 (A 2,174,
B 2,174), `late_determined=False`.

*Axes* — `cos(v, p_CRA)` = +0.5667444435453722, `cos(v, p_CRI)` =
+0.6337418620832961, `cos(v, u)` = −0.04027293769593313, N1 p = 1.0,
`beats_n1=False`, `||v||` = 90.14934904214671.

*Positive control* — passed; monotonic; `frac50` = 0.5469054675795995, CI
[0.48253102547211657, 0.6431535824253682]; `delta50` = −34.052728101022154;
`frac0` = 0.033645346601332594; `delta0` = +20.63359395062885; n_total 291.

*Proliferation quartiles* — library = total UMI in the raw matrix; CRA q25
0.15076169486925256, q75 0.9263986645926925, bottom 167, top 167, n_NT 666;
CRI q25 0.15641638362658056, q75 0.9633279093128779, bottom 370, top 370,
n_NT 1,479. Tirosh S genes found 42, missing `MLF1IP`; G2M found 52, missing
`FAM64A`, `HN1`.

*Screen* — CRA n=167 perturbations, `rho_prolif` = +0.8280530161707566
(p = 2.56e−43), `n_cos_q` 0, `n_toward_young` 0; CRI n=173,
`rho_prolif` = +0.8620017521821158 (p = 2.46e−52), `n_cos_q` 25, all 25 with
positive cosine, `n_toward_young` 0. Identity genes present: COL1A1, COL1A2,
FN1, LUM, PDGFRA, PDGFRB, POSTN, PRRX1, SERPINH1, VIM (0 missing).

*Reading* — `reproduction_pass = False`, 1 of 4 named hits reproduced; fired
flags `screen_tracks_proliferation`, `passage_axis_not_donor_age`,
`pipeline_disagrees`; all four named hits labelled `not_interpreted`.
Recorded failure: "Fewer than 3 of 4 named hits have cos_g>0 with N_NT FDR
q<=0.05 (1/4)."

### SENG_WHY — diagnostic, no gate (2025-09-22 19:02 to 19:31)

| file | kind | what it is |
|---|---|---|
| `src/seng_why.py` (697 lines) | script | four hypotheses for why the four named hits missed the gate |
| `results/seng_why/figshare_{30898748,30899018}.json` + `_files.csv` | result | both figshare records, re-fetched |
| `results/seng_why/obs_columns.csv` (127,926 B), `obs_low_cardinality.csv` | result | a second, independent `.obs` inventory |
| `results/seng_why/sample_name_by_group.csv` | result | `SampleName` against role |
| `results/seng_why/axis_by_late.csv` | result | the axis refit with PD26 as "late" |
| `results/seng_why/cos_by_late.csv`, `cos_all_genes.csv`, `rrej_vs_cos.csv` | result | per-perturbation cosines and the authors' R_rej recomputed |
| `results/seng_why/named_*.csv` | result | the same three tables restricted to the four named hits |
| `results/seng_why/run_report.txt`, `summary.json` | result | run log and fired keys |
| `DEBUG_SENG.md` (164 lines) | findings | the diagnostic write-up |

**Pre-registration: none, and that is stated in the artefacts themselves.**
There is no flag file under `results/seng_why/`; the run log opens with
"diagnostic only; FINDINGS_SENG.md is not written", and `DEBUG_SENG.md` says
the fired key `pipeline_disagrees` stands and no gate was rerun. The
`assumptions` block in `summary.json` lists six choices made **after** numbers
existed (R_rej gene selection, Pearson, sign convention, the log2FC recipe,
which PDs count as aging, and how unmapped genes were z-scored). Those are
post-hoc by construction. They are disclosed there, and are disclosed again
here.

**Numbers computed on this h5ad:**

- `X` enc csr_matrix, shape (36,514, 36,601), nnz 148,239,518, obs_cols 814,
  var_cols 0
- `n_WT` 6,809, `n_CRA` 14,333, `n_CRI` 15,372, **guide_cols 804**
- roles: WT 6,809, perturbed 24,610, NT 2,145, unassigned 2,950; ties 2,680;
  modality disagreements 270
- `SampleName` by role: 0/WT/PD14 581, 1/WT/PD26 1,880, 2/WT/PD32 4,348,
  3/CRA 12,264, 3/CRI 12,346, 3/NT 2,145, 3/unassigned 2,950
- Spearman of the authors' recomputed R_rej against this project's cosine:
  pooled −0.9838760054407149 (p = 9.7e−255); CRA −0.9825873210549642 (n=167);
  CRI −0.9864159532009773 (n=173)
- axis with PD26 as late: `cos(v,u)` = −0.005987805105684607, p = 0.736318407960199,
  `rho_prolif` CRA +0.81907871004978, CRI +0.8588612670408982, `||v||` 86.15849167004308
- PD32 cosine reproduced against `results/seng/screen.csv` to max abs
  difference 9.71445146547012e−17
- gene mapping: 23,485 panel, 16,781 present, 36,601 file genes, 16,771 mapped,
  19,830 unmapped, 10 duplicate maps
- named hits, recomputed R_rej vs the paper's Table 1: CRA_E2F3 n=128
  R_rej −0.2809089293905073 (published −0.53); CRA_EZH2 n=107 −0.20297679394214646
  (−0.36); CRI_STAT3 n=96 −0.1548712170415503 (−0.41); CRI_ZFX n=157
  −0.35530223525825205 (−0.51)
- fired key: `unexplained`; `metric_explains`, `reference_explains`,
  `gene_space_explains`, `file_explains` all False; elapsed 1,512.1 s

### PAPER_FIGS — verification of the stored SENG axis (undated in-directory)

| file | kind | what it is |
|---|---|---|
| `src/paper_figs_seng_null.py` (152 lines) | script | recomputes the SENG null for the figure |
| `results/paper_figs/seng_cos_vu_null.csv` | result | 200 permutation draws per late-PD choice |
| `results/paper_figs/seng_cos_vu_null_report.txt` | result | stored-vs-recomputed comparison |

**Pre-registration: none** (a figure script; `results/paper_figs/` has no flag).

**Numbers computed on this h5ad:** PD32 `max |v − axes.npz v|` = 0.00e+00;
`cos(v,u)` stored −0.0402729376959 vs recomputed identical to 3.47e−17,
p = 1 both, `||v||` differing by 1.42e−14, null mean −0.00111, sd 0.00635,
`n_null_ge_obs` 200/200. PD26 `cos(v,u)` −0.00598780510568 reproduced to
6.94e−18, p 0.73631840796, null mean −0.00084, sd 0.00690,
`n_null_ge_obs` 147/200.

### SOUTH6 (2025-09-27) — one file only

`results/seng/seurat_cc.genes_tirosh2016.tsv` has mtime 2025-09-27 16:34, five
days after the rest of `results/seng/`, inside the SOUTH5–SOUTH6 window.
**Recorded as an ambiguity:** a later task wrote into a previous task's results
directory, which CLAUDE.md forbids. Whether the SENG numbers above were
computed with this file or with an earlier copy of it cannot be determined from
the tree, because the file was overwritten rather than added beside the
original. Nothing in this session touches it.

## 1.3 This task's own overlap with the earlier work

`results/senginv/INVENTORY_SENGINV.md` (the SENGINV run, this session, earlier
today) is the **third** inventory of the same `.obs`, after
`results/seng/inventory.json` and `results/seng_why/obs_columns.csv`. Where the
three can be compared they agree, except where a definition differs:

| quantity | SENG (`inventory.json`) | SENG_WHY (`run_report.txt`) | SENGINV (this session) | same? |
|---|---|---|---|---|
| shape | — | 36,514 × 36,601 | 36,514 × 36,601 | yes |
| nnz | 148,239,518 | 148,239,518 | 148,239,518 | yes |
| `.obs` columns | — | 814 | 814 | yes |
| guide columns | 804 (NT included) | 804 | **798** (NT excluded) | **definition differs** |
| n_WT | 6,809 | 6,809 | 6,809 | yes |
| PD14 / PD26 / PD32 | 581 / 1,880 / 4,348 | 581 / 1,880 / 4,348 | 581 / 1,880 / 4,348 | yes |
| ties at a positive max | 2,680 | 2,680 | **2,624** | **definition differs** |
| max value in `X` | **2,960.0** (full scan) | — | 2,446 (10M-value sample) | **sample vs full scan** |
| panel genes present / missing | 16,781 / 6,704 | 16,781 / 6,704 | 16,781 / 6,704 | yes |
| distinct mapped symbols | — | 16,771 | 16,771 | yes |
| unmapped file genes | — | 19,830 | 19,830 | yes |

The three rows that differ, and why, stated rather than reconciled:

- **798 vs 804 guide columns.** SENGINV counted only columns matching
  `<TARGET>_<CRA|CRI><k>`; SENG and SENG_WHY also counted `NT_1`..`NT_6` as
  guide columns. 798 + 6 = 804. Neither is wrong; they answer different
  questions, and the task wording ("distinct target TFs") is what led SENGINV
  to exclude NT.
- **2,624 vs 2,680 ties.** Same cause: a tie over 798 columns is not a tie over
  804. Counted directly in this session (the only computation item 1 performed,
  on the guide columns only): **199** cells are tied over 804 columns but have a
  unique maximum over 798; **45** are tied over 804 and have no positive value at
  all over 798; and **188** run the other way — tied over 798 but unique over
  804, because an `NT_*` column exceeded the tied targeting maximum and became
  the unique winner. 199 + 45 − 188 = 56, the difference between the two
  published tie counts. Neither pool is the file's own; both are task choices.
- **max 2,446 vs 2,960.** SENGINV sampled 10,000,000 of the 148,239,518 stored
  values; SENG scanned all of them. **The SENG number is the correct maximum
  and the SENGINV figure should be read as a sample maximum only.** SENGINV's
  integrality claim is likewise a sample claim; SENG's `n_fractional=0` over
  the full nnz is the stronger statement. This is recorded here as a
  limitation of the SENGINV file, which is not edited.

The per-target cell counts also differ from SENG's, by construction. SENGINV's
"convention B" argmax over 798 columns gives CRA_E2F3 128, CRA_EZH2 109,
CRI_STAT3 101, CRI_ZFX 163. SENG's assignment gives 128, 107, 96, 157 — it used
804 columns and additionally dropped the 270 cells whose guide call disagreed
with the `CRA`/`CRI` indicator. E2F3 agrees exactly; the other three do not.
**Recorded as an ambiguity, not resolved:** the object supports neither rule
over the other, as section 4 of the SENGINV file already says.

## 1.4 What this audit does not cover

- Nothing was rerun, so no number above is verified by this session. The
  agreements in 1.3 are agreements between *stored* numbers, which is weaker
  than reproduction and is not the reproduction check CLAUDE.md requires before
  computing something new. This task computed nothing that depends on them.
- `results/seng/screen.csv` has 46 columns and 340 data rows of per-perturbation
  statistics that are not summarised here beyond the two `screen_summary.json`
  aggregates. Anyone needing those should read the CSV.
- The SENG reading is `reproduction_pass=False` with three fired flags. That
  reading stands as written in `FINDINGS_SENG.md`; nothing here supersedes,
  edits or re-opens it.

