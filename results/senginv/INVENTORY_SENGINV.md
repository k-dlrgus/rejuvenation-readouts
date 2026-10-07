# SENGINV - data inventory, Sengstack et al. Perturb-seq h5ad

Inventory only. Nothing here is a pre-registered result and nothing is
interpreted. Ambiguities are recorded as ambiguities, not resolved.

## 0. Provenance

- file: `data/raw/seng/exp3_merged_not_normalized.h5ad`
  (untracked; `.gitignore` rule `data/raw/`)
- bytes on disk: 1,426,871,504
- figshare article 30898748 v2: Systematic identification of single transcription factor perturbations that drive cellular and tissue rejuvenation.
  - doi `10.6084/m9.figshare.30898748.v2`, article size 1,426,871,504 B, CC BY 4.0, published 2025-12-25T06:54:05Z
  - file `exp3_merged_not_normalized.h5ad` 1,426,871,504 B md5 `bfccaa94a916c22e4b0793913a04128b`
- figshare article 30899018 v2: Systematic identification of single transcription factor perturbations that drive cellular and tissue rejuvenation: Mouse liver RNAseq data
  - doi `10.6084/m9.figshare.30899018.v2`, article size 1,669,672 B, CC BY 4.0, published 2025-12-25T07:02:58Z
  - file `mouse_rnaseq_sample_name.csv` 233 B md5 `7f43663064a09f8e0e3aadf5918b1f9a`
  - file `mouse_liver_raw_count_figshare.csv` 1,669,439 B md5 `818a9e405cfe609f7c7f56f0702d6400`
- md5 of the local copy: `bfccaa94a916c22e4b0793913a04128b`

Download status, stated plainly:

- **article 30898748 was not downloaded in this session.** A byte-identical
  copy was already on disk at `data/raw/seng/` from the earlier SENG task.
  Its md5 was recomputed here and equals the md5 figshare reports for the
  v2 file, so the local copy is the v2 file named in the task.
- **article 30899018 was not downloaded.** Checked by API only, per the
  task instruction: its title and description say mouse liver RNA-seq raw
  counts (two CSVs, 1,669,672 B in all), so it was left alone. Title and size
  are reported above and nothing else about it is reported.

## 1. Shape and structure

- `.shape` = (36514, 36601)  ->  n_cells = 36,514, n_genes = 36,601
- top-level HDF5 keys: ['X', 'obs', 'var']
- root attrs: {}  (no `encoding-type` / anndata-version written)
- `.X` encoding: csr_matrix v0.1.0, nnz = 148,239,518, dtype = float32
- density = 0.1109
- `.obsm` keys: **NONE** (no `obsm` group in the file)
- `.layers` keys: **NONE** (no `layers` group in the file)
- `.uns` keys: **NONE** (no `uns` group in the file)
- `.varm`, `.obsp`, `.varp`, `.raw`: **NONE**
- `.var` columns: [] (index only, index name `_index`)
- `.obs` index name: `_index`
- number of `.obs` columns: 814
- `.obs` also holds a legacy `__categories` group, which is not a column:
  the levels of the categorical columns live there, h5ad's pre-0.8 layout.
  Its members: ['SampleName']

### Are raw integer counts present, and where?

- `.X` is the only matrix in the file: CSR, float32.
- sampled 10,000,000 stored values (first 5M and last 5M of `X/data`):
  - every sampled value equals its own rounding: True
  - min 1, max 2446, n distinct 965
  - first 20 stored values: [2.0, 2.0, 2.0, 5.0, 1.0, 3.0, 2.0, 1.0, 1.0, 8.0, 1.0, 1.0, 7.0, 9.0, 9.0, 1.0, 1.0, 1.0, 1.0, 1.0]
- so `.X` holds integer-valued counts stored as float32, consistent with
  the file name `not_normalized`.
- there is no second copy of the counts: no `.raw`, no `.layers`.

## 1b. Full `.obs` column list with dtype

`.obs` has 814 columns. 798 match `<NAME>_(CRA|CRI)<n>` and are per-guide
UMI count columns (section 4); the other 16 are listed below in full.

### Non-guide columns, all 16 of them, in stored order

| # | column | stored dtype | n nonzero | min | max | n distinct |
|---|---|---|---|---|---|---|
| 1 | `PD14` | bool | 581 | 0 | 1 | 2 |
| 2 | `PD26` | bool | 1880 | 0 | 1 | 2 |
| 3 | `PD32` | bool | 4348 | 0 | 1 | 2 |
| 4 | `CRI` | int64 | 15372 | 0 | 1 | 2 |
| 5 | `CRA` | int64 | 14333 | 0 | 1 | 2 |
| 6 | `n_counts` | float32 | 36514 | 500 | 90340 | 23366 |
| 7 | `log_counts` | float32 | 36514 | 6.21461 | 11.4113 | 23366 |
| 8 | `n_genes` | int64 | 36514 | 33 | 8602 | 6402 |
| 9 | `NT_5` | float64 | 1790 | 0 | 105 | 76 |
| 10 | `NT_4` | float64 | 1199 | 0 | 212 | 65 |
| 11 | `NT_1` | float64 | 1141 | 0 | 67 | 49 |
| 12 | `NT_2` | float64 | 863 | 0 | 33 | 30 |
| 13 | `NT_3` | float64 | 1016 | 0 | 76 | 63 |
| 14 | `NT_6` | float64 | 184 | 0 | 4 | 5 |
| 15 | `WT` | int64 | 6809 | 0 | 1 | 2 |
| 16 | `SampleName` | categorical codes [int8] n_cat=4 (legacy obs/__categories) | - | - | - | 4 |

### Distinct values and counts for the non-guide columns

One `.obs` column is categorical: `SampleName`, stored in the legacy h5ad
layout as int8 codes whose levels are in `obs/__categories/SampleName`.
Its distinct values and counts are given below by level name. The other
non-guide columns are plain bool or numeric; value counts are given for
every one with at most 50 distinct values, percentiles for the rest.

**`PD14`** (stored bool)

| value | n cells |
|---|---|
| `False` | 35933 |
| `True` | 581 |

**`PD26`** (stored bool)

| value | n cells |
|---|---|
| `False` | 34634 |
| `True` | 1880 |

**`PD32`** (stored bool)

| value | n cells |
|---|---|
| `False` | 32166 |
| `True` | 4348 |

**`CRI`** (stored int64)

| value | n cells |
|---|---|
| `0` | 21142 |
| `1` | 15372 |

**`CRA`** (stored int64)

| value | n cells |
|---|---|
| `0` | 22181 |
| `1` | 14333 |

**`n_counts`** (stored float32)

- 23366 distinct values; percentiles 0/1/25/50/75/99/100 = 500 / 843.13 / 12335 / 19777 / 26324 / 50208.9 / 90340
- n zero = 0, mean = 19801.4

**`log_counts`** (stored float32)

- 23366 distinct values; percentiles 0/1/25/50/75/99/100 = 6.21461 / 6.73712 / 9.4202 / 9.89227 / 10.1782 / 10.8239 / 11.4113
- n zero = 0, mean = 9.66572

**`n_genes`** (stored int64)

- 6402 distinct values; percentiles 0/1/25/50/75/99/100 = 33 / 483.13 / 3342 / 4335.5 / 5006 / 6713.09 / 8602
- n zero = 0, mean = 4060.02

**`NT_5`** (stored float64)

- 76 distinct values; percentiles 0/1/25/50/75/99/100 = 0 / 0 / 0 / 0 / 0 / 19 / 105
- n zero = 34724, mean = 0.507504

**`NT_4`** (stored float64)

- 65 distinct values; percentiles 0/1/25/50/75/99/100 = 0 / 0 / 0 / 0 / 0 / 10 / 212
- n zero = 35315, mean = 0.301282

**`NT_1`** (stored float64)

| value | n cells |
|---|---|
| `0.0` | 35373 |
| `1.0` | 384 |
| `2.0` | 104 |
| `3.0` | 45 |
| `4.0` | 42 |
| `5.0` | 37 |
| `6.0` | 38 |
| `7.0` | 47 |
| `8.0` | 33 |
| `9.0` | 36 |
| `10.0` | 38 |
| `11.0` | 33 |
| `12.0` | 27 |
| `13.0` | 30 |
| `14.0` | 22 |
| `15.0` | 24 |
| `16.0` | 25 |
| `17.0` | 22 |
| `18.0` | 17 |
| `19.0` | 14 |
| `20.0` | 13 |
| `21.0` | 15 |
| `22.0` | 17 |
| `23.0` | 6 |
| `24.0` | 6 |
| `25.0` | 7 |
| `26.0` | 4 |
| `27.0` | 7 |
| `28.0` | 4 |
| `29.0` | 8 |
| `30.0` | 4 |
| `31.0` | 2 |
| `32.0` | 5 |
| `33.0` | 3 |
| `34.0` | 3 |
| `36.0` | 2 |
| `37.0` | 1 |
| `39.0` | 2 |
| `41.0` | 2 |
| `42.0` | 1 |
| `43.0` | 1 |
| `44.0` | 2 |
| `47.0` | 1 |
| `49.0` | 1 |
| `50.0` | 2 |
| `52.0` | 1 |
| `54.0` | 1 |
| `58.0` | 1 |
| `67.0` | 1 |

**`NT_2`** (stored float64)

| value | n cells |
|---|---|
| `0.0` | 35651 |
| `1.0` | 272 |
| `2.0` | 95 |
| `3.0` | 71 |
| `4.0` | 68 |
| `5.0` | 62 |
| `6.0` | 48 |
| `7.0` | 34 |
| `8.0` | 41 |
| `9.0` | 31 |
| `10.0` | 26 |
| `11.0` | 27 |
| `12.0` | 10 |
| `13.0` | 9 |
| `14.0` | 10 |
| `15.0` | 11 |
| `16.0` | 6 |
| `17.0` | 7 |
| `18.0` | 8 |
| `19.0` | 2 |
| `20.0` | 5 |
| `21.0` | 2 |
| `22.0` | 2 |
| `23.0` | 4 |
| `24.0` | 5 |
| `26.0` | 2 |
| `28.0` | 1 |
| `30.0` | 2 |
| `31.0` | 1 |
| `33.0` | 1 |

**`NT_3`** (stored float64)

- 63 distinct values; percentiles 0/1/25/50/75/99/100 = 0 / 0 / 0 / 0 / 0 / 10 / 76
- n zero = 35498, mean = 0.284603

**`NT_6`** (stored float64)

| value | n cells |
|---|---|
| `0.0` | 36330 |
| `1.0` | 160 |
| `2.0` | 20 |
| `3.0` | 2 |
| `4.0` | 2 |

**`WT`** (stored int64)

| value | n cells |
|---|---|
| `0` | 29705 |
| `1` | 6809 |

**`SampleName`** (stored categorical codes [int8] n_cat=4 (legacy obs/__categories))

| value | n cells |
|---|---|
| `UNKNOWN` | 29705 |
| `PD32` | 4348 |
| `PD26` | 1880 |
| `PD14` | 581 |

## 2. Are early-passage (young) cells present in this object?

- `.obs` columns whose name is `PD<number>`: ['PD14', 'PD26', 'PD32']
- the numbers in those names are [14, 26, 32].
- **Recorded as an ambiguity.** PD conventionally means population
  doublings, so `PD14` would be the earliest-passage arm of the three, but
  the object itself carries no definition: `.uns` is absent and there is no
  description field anywhere in `.obs` or `.var`. That `PD14` is the
  early-passage (young) arm rests on the column name alone. The three
  columns are stored as bool, so they are indicators and not counts, but
  what the indicator means - doublings, passage number, culture age - is
  not stated in the object.

### Per-PD-column occupancy

| column | n cells nonzero | n distinct values | value counts (value:n cells) |
|---|---|---|---|
| `PD14` | 581 | 2 | 0:35933, 1:581 |
| `PD26` | 1880 | 2 | 0:34634, 1:1880 |
| `PD32` | 4348 | 2 | 0:32166, 1:4348 |

### Do the PD columns partition the cells?

| n PD columns nonzero for a cell | n cells |
|---|---|
| 0 | 29705 |
| 1 | 6809 |

- cells with exactly one nonzero PD column: 6809
- cells with no nonzero PD column: 29705
- cells with more than one nonzero PD column: 0

### Lane / batch / library / sample / 10x-channel columns

- the only such column in the object is `SampleName` (stored categorical codes [int8] n_cat=4 (legacy obs/__categories)).
- there is no column named lane, batch, library, channel, 10x, gem, run,
  donor, replicate or experiment. The 16 non-guide columns listed in 1b
  are the whole of what exists.

| SampleName | n cells |
|---|---|
| `UNKNOWN` | 29705 |
| `PD32` | 4348 |
| `PD26` | 1880 |
| `PD14` | 581 |

### The `.obs` index itself (10x barcode suffixes)

The `.obs` index is 10x cell barcodes carrying a suffix. That suffix is the
only other per-cell grouping in the file. It is part of the index string and
not a column, so it is reported here as what it is.

- first 3 index entries: ['AAACCTGAGCAATCTC-1', 'AAACCTGAGCTGGAAC-1', 'AAACCTGAGGCCATAG-1']
- duplicated index entries: 0
- distinct bare barcodes with the suffix stripped: 36,202 of 36,514 rows

| barcode suffix | n cells |
|---|---|
| `_I` | 15372 |
| `_A` | 14333 |
| `-1` | 6809 |

| SampleName | `_I` | `_A` | `-1` |
|---|---|---|---|
| `UNKNOWN` | 15372 | 14333 | 0 |
| `PD32` | 0 | 0 | 4348 |
| `PD26` | 0 | 0 | 1880 |
| `PD14` | 0 | 0 | 581 |

| barcode suffix | `PD14`>0 | `PD26`>0 | `PD32`>0 | `CRA`>0 | `CRI`>0 | `WT`>0 |
|---|---|---|---|---|---|---|
| `_I` | 0 | 0 | 0 | 0 | 15372 | 0 |
| `_A` | 0 | 0 | 0 | 14333 | 0 | 0 |
| `-1` | 581 | 1880 | 4348 | 0 | 0 | 6809 |

**Recorded as an ambiguity.** A barcode suffix is what `cellranger aggr` or
an anndata concatenation writes to keep barcodes unique across inputs, so
the suffix groups are candidates for separate libraries. The file does not
say that, and does not say what `-1`, `_A` or `_I` stand for. Nothing here
is a claim that the suffix is a lane, a library or a sequencing run.

### Cross-tabulation: PD columns against SampleName

| SampleName | `PD14`>0 | `PD26`>0 | `PD32`>0 | no PD>0 | total |
|---|---|---|---|---|---|
| `UNKNOWN` | 0 | 0 | 0 | 29705 | 29705 |
| `PD32` | 0 | 0 | 4348 | 0 | 4348 |
| `PD26` | 0 | 1880 | 0 | 0 | 1880 |
| `PD14` | 581 | 0 | 0 | 0 | 581 |

### Answer to question 2 as asked

- **early-passage cells, reading `PD14` as the earliest arm: present,**
  **n = 581.**
- identified by the bool `.obs` column `PD14`, True for those 581 cells;
  equivalently by `SampleName == 'PD14'` (581 cells; the two labels agree
  on 36514 of 36514 cells).
- the three PD arms together are 6809 cells. That is also exactly the number
  with `WT`>0 (6809) and exactly the number whose barcode suffix is `-1`.
- columns that could show a sequencing run: only `SampleName`, whose four
  levels are the three PD arms and `UNKNOWN`, with 29705 cells in `UNKNOWN`.
- measured, not assumed: 28523 cells have a nonzero value in at least one of
  the 798 guide columns. Of those, 28523 are in `SampleName == 'UNKNOWN'` and
  0 are in a PD arm. So `SampleName` separates the PD cells from the
  cells carrying guides rather than placing them in one group. The barcode
  suffixes split the same way: `-1` for all PD cells, `_A` and `_I` for all
  `UNKNOWN` cells. There is no lane, batch, library, channel or run column.
- **Recorded as an ambiguity, not resolved.** Whether the PD cells were
  sequenced in the same run as the perturbed cells cannot be read off this
  object. A separate `SampleName` level and a separate barcode suffix are
  consistent with a separate run and equally with a separate `cellranger`
  input within one run. The file records no run, flowcell, chemistry or
  date identifier that would decide it.

## 3. Non-targeting control cells

- `.obs` columns matching `NT_<n>`: ['NT_1', 'NT_2', 'NT_3', 'NT_4', 'NT_5', 'NT_6']
- there is also a column `WT` (stored int64).
- **Recorded as an ambiguity.** The object does not say that `NT_*` means
  non-targeting guide, nor what `WT` is, nor whether `WT` is a separate
  unperturbed arm or an overlapping label. The tables below report the
  numbers without choosing a reading.

| column | n cells nonzero | min nonzero | median nonzero | max |
|---|---|---|---|---|
| `NT_1` | 1141 | 1 | 4 | 67 |
| `NT_2` | 863 | 1 | 3 | 33 |
| `NT_3` | 1016 | 1 | 3 | 76 |
| `NT_4` | 1199 | 1 | 2 | 212 |
| `NT_5` | 1790 | 1 | 3 | 105 |
| `NT_6` | 184 | 1 | 1 | 4 |
| `WT` | 6809 | 1 | 1 | 1 |

- cells with at least one nonzero `NT_*` column: 4417
- cells with nonzero `WT`: 6809
- cells nonzero in both an `NT_*` and `WT`: 0
- number of nonzero `NT_*` columns per cell:

| n NT_* columns nonzero | n cells |
|---|---|
| 0 | 32097 |
| 1 | 2893 |
| 2 | 1320 |
| 3 | 160 |
| 4 | 40 |
| 5 | 4 |

### `NT_*` and `WT` against the PD columns

| group | `PD14`>0 | `PD26`>0 | `PD32`>0 | no PD>0 | total |
|---|---|---|---|---|---|
| any `NT_*`>0 | 0 | 0 | 0 | 4417 | 4417 |
| `WT`>0 | 581 | 1880 | 4348 | 0 | 6809 |
| `NT_1`>0 | 0 | 0 | 0 | 1141 | 1141 |
| `NT_2`>0 | 0 | 0 | 0 | 863 | 863 |
| `NT_3`>0 | 0 | 0 | 0 | 1016 | 1016 |
| `NT_4`>0 | 0 | 0 | 0 | 1199 | 1199 |
| `NT_5`>0 | 0 | 0 | 0 | 1790 | 1790 |
| `NT_6`>0 | 0 | 0 | 0 | 184 | 184 |

### `NT_*` and `WT` against SampleName

| group | `UNKNOWN` | `PD32` | `PD26` | `PD14` |
|---|---|---|---|---|
| any `NT_*`>0 | 4417 | 0 | 0 | 0 |
| `WT`>0 | 0 | 4348 | 1880 | 581 |
| `NT_1`>0 | 1141 | 0 | 0 | 0 |
| `NT_2`>0 | 863 | 0 | 0 | 0 |
| `NT_3`>0 | 1016 | 0 | 0 | 0 |
| `NT_4`>0 | 1199 | 0 | 0 | 0 |
| `NT_5`>0 | 1790 | 0 | 0 | 0 |
| `NT_6`>0 | 184 | 0 | 0 | 0 |

## 4. Guide / perturbation assignment

- the per-guide UMI counts are `.obs` columns, one column per guide, named
  `<TARGET>_<CRA|CRI><1|2>`. There is no `obsm['guide_counts']`, no guide
  `layers` entry and no separate guide-by-cell matrix: `.obs` is where they
  live. This matches the figshare description, which says the two guides per
  target were counted separately.
- number of guide columns: 798

- per-guide-column table: `results/senginv/guide_columns.csv` (798 rows)
- distinct stored dtypes among guide columns: ['float64']
- replicate suffixes present: [np.int64(1), np.int64(2)]
- all guide UMI values integral: True
- distinct target names across guide columns: **200**
- distinct target names appearing in CRA columns: 200
- distinct target names appearing in CRI columns: 200
- guide columns per target x direction:
  - 1 column(s): 2 target x direction pairs
  - 2 column(s): 398 target x direction pairs
- 200 x 2 directions x 2 guides would be 800 columns; there are 798.
- pairs not holding exactly 2 guide columns: [('DLX5', 'CRI', 1), ('PHOX2A', 'CRA', 1)]
- targets with both CRA and CRI guide columns: 200
- targets with CRA columns only: []
- targets with CRI columns only: []
- **Recorded as an ambiguity.** 'Target TF' cannot be read off the object:
  the target names are just the prefixes of the guide column names, and the
  object says nothing about which of them are transcription factors. Section
  9 reports whether each prefix is also a symbol in `.var`, which is not the
  same question.

### How a cell gets assigned to a target

**Nothing in the object assigns a cell to a target.** There is no
`guide_identity`, `perturbation`, `target`, `gene_target`, `condition` or
similar single-value column - the 16 non-guide columns are listed in full
in 1b and none of them is an assignment. An assignment would have to be
derived from the 798 UMI columns by a rule that is not recorded in the file.
**Recorded as an ambiguity and not resolved here.** The counts below are
reported under two explicitly stated counting conventions chosen only so
that a number can be printed. Neither is an assignment rule and neither is
endorsed:

- convention A (any UMI): a cell is counted toward target T if any guide
  column of T has UMI > 0. A cell may be counted toward many targets, so
  the column does not sum to n_cells.
- convention B (argmax): a cell is counted toward the single target holding
  the strictly largest guide UMI over all 798 guide columns. Cells with a
  tied maximum, and cells whose guide columns are all zero, are dropped.

- cells with every guide column zero: 7991
- cells with a tied maximum across guide columns: 2624
- cells with a unique argmax guide: 25899

- full sorted table also written to `results/senginv/cells_per_target.csv`
  (400 target x direction rows)

- convention A (any UMI): min 2, median 147, max 686, sum over rows 68078
- convention B (argmax): min 0, median 59.5, max 212, sum over rows 25899

Aggregated over the two directions, per target name:

- convention A (any UMI): min 56, median 336.5, max 791
- convention B (argmax): min 15, median 128.5, max 323

### Full sorted table: cells per target x direction

| rank | target | direction | n guide cols | cells any UMI>0 | cells argmax | total guide UMI |
|---|---|---|---|---|---|---|
| 1 | `E2F4` | CRI | 2 | 686 | 212 | 5910 |
| 2 | `ZFX` | CRI | 2 | 522 | 163 | 4647 |
| 3 | `KLF6` | CRI | 2 | 502 | 154 | 5410 |
| 4 | `EOMES` | CRI | 2 | 480 | 162 | 5838 |
| 5 | `SOX2` | CRI | 2 | 480 | 131 | 5718 |
| 6 | `ETV5` | CRI | 2 | 459 | 114 | 3604 |
| 7 | `HOXB4` | CRI | 2 | 455 | 120 | 4297 |
| 8 | `ZFHX3` | CRI | 2 | 452 | 130 | 4392 |
| 9 | `ARNTL2` | CRI | 2 | 434 | 107 | 3372 |
| 10 | `BSX` | CRI | 2 | 425 | 105 | 3019 |
| 11 | `VAX1` | CRI | 2 | 413 | 116 | 2828 |
| 12 | `SP1` | CRI | 2 | 408 | 107 | 3463 |
| 13 | `ATG3` | CRI | 2 | 406 | 106 | 3828 |
| 14 | `ZNF281` | CRI | 2 | 406 | 114 | 2949 |
| 15 | `HOXA6` | CRI | 2 | 405 | 127 | 2919 |
| 16 | `PAX8` | CRI | 2 | 401 | 118 | 3241 |
| 17 | `CDX1` | CRI | 2 | 397 | 158 | 3763 |
| 18 | `MYC` | CRI | 2 | 394 | 129 | 3852 |
| 19 | `SIM1` | CRI | 2 | 387 | 91 | 2497 |
| 20 | `SOX17` | CRI | 2 | 375 | 97 | 2408 |
| 21 | `ZKSCAN3` | CRI | 2 | 368 | 126 | 3975 |
| 22 | `ALX1` | CRI | 2 | 365 | 86 | 2560 |
| 23 | `SOX11` | CRI | 2 | 361 | 97 | 1858 |
| 24 | `TFAM` | CRI | 2 | 359 | 97 | 3343 |
| 25 | `FOXP2` | CRI | 2 | 358 | 130 | 4603 |
| 26 | `PARP1` | CRI | 2 | 354 | 89 | 2112 |
| 27 | `SOX9` | CRI | 2 | 354 | 131 | 2342 |
| 28 | `ZNF318` | CRI | 2 | 342 | 90 | 2654 |
| 29 | `E2F7` | CRI | 2 | 337 | 101 | 3498 |
| 30 | `ATG4C` | CRI | 2 | 334 | 90 | 2669 |
| 31 | `TBX3` | CRI | 2 | 334 | 90 | 2601 |
| 32 | `MAZ` | CRI | 2 | 332 | 95 | 2277 |
| 33 | `LHX1` | CRI | 2 | 330 | 94 | 2640 |
| 34 | `DLX6` | CRI | 2 | 329 | 94 | 3287 |
| 35 | `ATF7` | CRI | 2 | 328 | 100 | 3443 |
| 36 | `PAX4` | CRI | 2 | 324 | 106 | 4251 |
| 37 | `TFEB` | CRI | 2 | 323 | 80 | 2678 |
| 38 | `ZNF35` | CRI | 2 | 322 | 83 | 2414 |
| 39 | `FOXB2` | CRI | 2 | 318 | 113 | 3417 |
| 40 | `ZBTB2` | CRI | 2 | 318 | 116 | 2925 |
| 41 | `LHX5` | CRI | 2 | 315 | 99 | 1910 |
| 42 | `NFKB1` | CRI | 2 | 313 | 75 | 3066 |
| 43 | `E2F3` | CRI | 2 | 309 | 79 | 2229 |
| 44 | `HMGB1` | CRI | 2 | 306 | 101 | 2546 |
| 45 | `NR4A2` | CRI | 2 | 305 | 87 | 3064 |
| 46 | `DNMT1` | CRI | 2 | 303 | 63 | 1930 |
| 47 | `MYC` | CRA | 2 | 297 | 194 | 4498 |
| 48 | `PGR` | CRI | 2 | 297 | 96 | 2596 |
| 49 | `HSF1` | CRI | 2 | 296 | 95 | 3692 |
| 50 | `TWIST1` | CRI | 2 | 295 | 81 | 2910 |
| 51 | `SOX7` | CRI | 2 | 293 | 78 | 1468 |
| 52 | `TCF3` | CRI | 2 | 293 | 102 | 3228 |
| 53 | `VEZF1` | CRI | 2 | 292 | 103 | 2740 |
| 54 | `PBX2` | CRI | 2 | 291 | 75 | 2591 |
| 55 | `CREB3L1` | CRI | 2 | 289 | 74 | 1879 |
| 56 | `ZNF384` | CRA | 2 | 289 | 136 | 3395 |
| 57 | `FOXS1` | CRI | 2 | 287 | 113 | 2909 |
| 58 | `ONECUT1` | CRI | 2 | 284 | 88 | 2489 |
| 59 | `FOSL1` | CRA | 2 | 283 | 156 | 4183 |
| 60 | `HOXB13` | CRA | 2 | 282 | 157 | 3276 |
| 61 | `EZH2` | CRI | 2 | 281 | 93 | 3294 |
| 62 | `FOXO4` | CRI | 2 | 280 | 66 | 2024 |
| 63 | `PHOX2A` | CRI | 2 | 277 | 62 | 1548 |
| 64 | `STAT1` | CRI | 2 | 277 | 86 | 2450 |
| 65 | `TLX2` | CRA | 2 | 273 | 143 | 5012 |
| 66 | `NFE2L2` | CRI | 2 | 272 | 93 | 2148 |
| 67 | `SOX15` | CRI | 2 | 270 | 78 | 2481 |
| 68 | `IRF5` | CRI | 2 | 269 | 90 | 2007 |
| 69 | `E2F1` | CRI | 2 | 268 | 96 | 2454 |
| 70 | `GATA2` | CRI | 2 | 262 | 72 | 1878 |
| 71 | `FERD3L` | CRI | 2 | 256 | 96 | 1939 |
| 72 | `E2F3` | CRA | 2 | 254 | 128 | 2722 |
| 73 | `FLI1` | CRI | 2 | 254 | 62 | 2105 |
| 74 | `GSX2` | CRI | 2 | 253 | 91 | 2211 |
| 75 | `JUND` | CRI | 2 | 252 | 72 | 2324 |
| 76 | `KLF3` | CRI | 2 | 251 | 79 | 2125 |
| 77 | `LCOR` | CRI | 2 | 250 | 103 | 2368 |
| 78 | `STAT4` | CRI | 2 | 250 | 68 | 3602 |
| 79 | `MRRF` | CRI | 2 | 249 | 82 | 2144 |
| 80 | `XBP1` | CRI | 2 | 245 | 60 | 1931 |
| 81 | `SOX15` | CRA | 2 | 244 | 143 | 2720 |
| 82 | `PPARGC1B` | CRA | 2 | 243 | 143 | 3924 |
| 83 | `LHX5` | CRA | 2 | 241 | 120 | 3731 |
| 84 | `SOX5` | CRI | 2 | 241 | 80 | 2361 |
| 85 | `PPARGC1B` | CRI | 2 | 238 | 73 | 1702 |
| 86 | `ARID5A` | CRI | 2 | 234 | 87 | 1447 |
| 87 | `TEAD3` | CRA | 2 | 234 | 119 | 1681 |
| 88 | `CREB1` | CRI | 2 | 231 | 70 | 2235 |
| 89 | `NRF1` | CRI | 2 | 231 | 49 | 1561 |
| 90 | `ATF4` | CRI | 2 | 227 | 75 | 1893 |
| 91 | `CPEB1` | CRA | 2 | 227 | 116 | 3293 |
| 92 | `ARX` | CRI | 2 | 225 | 68 | 2150 |
| 93 | `CEBPB` | CRI | 2 | 225 | 62 | 2189 |
| 94 | `DLX4` | CRI | 2 | 224 | 92 | 2452 |
| 95 | `XRCC5` | CRI | 2 | 224 | 80 | 1557 |
| 96 | `CREB3L1` | CRA | 2 | 223 | 125 | 3601 |
| 97 | `PLAU` | CRI | 2 | 223 | 55 | 1340 |
| 98 | `NFKB2` | CRA | 2 | 222 | 122 | 3870 |
| 99 | `PAX4` | CRA | 2 | 222 | 87 | 1479 |
| 100 | `TFCP2` | CRA | 2 | 222 | 123 | 2835 |
| 101 | `EZH2` | CRA | 2 | 221 | 109 | 4683 |
| 102 | `IRF7` | CRI | 2 | 221 | 56 | 1997 |
| 103 | `GATA3` | CRA | 2 | 220 | 86 | 1304 |
| 104 | `KLF4` | CRI | 2 | 219 | 39 | 1205 |
| 105 | `NPAS3` | CRI | 2 | 218 | 53 | 2053 |
| 106 | `STAT3` | CRI | 2 | 218 | 101 | 2260 |
| 107 | `ARNTL2` | CRA | 2 | 216 | 103 | 3257 |
| 108 | `FOXA3` | CRI | 2 | 216 | 82 | 1868 |
| 109 | `TLX2` | CRI | 2 | 212 | 90 | 629 |
| 110 | `PAX2` | CRI | 2 | 210 | 73 | 2288 |
| 111 | `ARX` | CRA | 2 | 209 | 115 | 2658 |
| 112 | `EGR1` | CRI | 2 | 209 | 72 | 1656 |
| 113 | `MEF2A` | CRI | 2 | 206 | 37 | 1390 |
| 114 | `GATA3` | CRI | 2 | 205 | 66 | 2298 |
| 115 | `DLX2` | CRA | 2 | 204 | 118 | 2841 |
| 116 | `FOXA3` | CRA | 2 | 204 | 97 | 1164 |
| 117 | `FOXF2` | CRI | 2 | 203 | 55 | 1849 |
| 118 | `CDX1` | CRA | 2 | 202 | 99 | 2982 |
| 119 | `HOXA2` | CRI | 2 | 202 | 52 | 1141 |
| 120 | `THRA` | CRI | 2 | 200 | 66 | 1221 |
| 121 | `MEF2A` | CRA | 2 | 199 | 99 | 2817 |
| 122 | `MYBL2` | CRI | 2 | 199 | 59 | 1245 |
| 123 | `GABARAPL1` | CRI | 2 | 198 | 60 | 1434 |
| 124 | `TFE3` | CRI | 2 | 198 | 47 | 1324 |
| 125 | `FOXB2` | CRA | 2 | 196 | 120 | 3367 |
| 126 | `POU3F2` | CRI | 2 | 196 | 50 | 764 |
| 127 | `SREBF1` | CRI | 2 | 196 | 46 | 1499 |
| 128 | `E2F5` | CRI | 2 | 195 | 77 | 1247 |
| 129 | `PGR` | CRA | 2 | 194 | 108 | 2440 |
| 130 | `SIM2` | CRA | 2 | 194 | 105 | 3727 |
| 131 | `NKX2-2` | CRA | 2 | 193 | 89 | 1922 |
| 132 | `POU5F1` | CRI | 2 | 193 | 36 | 977 |
| 133 | `GATA5` | CRI | 2 | 189 | 59 | 1599 |
| 134 | `MXD3` | CRA | 2 | 189 | 102 | 3169 |
| 135 | `PHOX2B` | CRI | 2 | 189 | 61 | 1822 |
| 136 | `APOD` | CRI | 2 | 188 | 75 | 2270 |
| 137 | `JUND` | CRA | 2 | 188 | 121 | 2410 |
| 138 | `FOXP1` | CRI | 2 | 187 | 69 | 1832 |
| 139 | `HSF2` | CRA | 2 | 187 | 124 | 2614 |
| 140 | `PRKRIR` | CRI | 2 | 187 | 54 | 1424 |
| 141 | `DLX6` | CRA | 2 | 186 | 91 | 2755 |
| 142 | `HOXA6` | CRA | 2 | 186 | 100 | 3225 |
| 143 | `ZBTB7C` | CRI | 2 | 186 | 66 | 1201 |
| 144 | `HOXD3` | CRA | 2 | 185 | 104 | 3203 |
| 145 | `ETS2` | CRA | 2 | 184 | 121 | 2104 |
| 146 | `GATA5` | CRA | 2 | 184 | 87 | 1741 |
| 147 | `HOXB7` | CRA | 2 | 184 | 104 | 1830 |
| 148 | `CTCFL` | CRI | 2 | 183 | 75 | 838 |
| 149 | `HOXC5` | CRI | 2 | 183 | 33 | 1482 |
| 150 | `CPEB1` | CRI | 2 | 182 | 68 | 2286 |
| 151 | `NKX1-1` | CRA | 2 | 181 | 94 | 1773 |
| 152 | `SREBF1` | CRA | 2 | 180 | 95 | 2320 |
| 153 | `VEZF1` | CRA | 2 | 180 | 97 | 1602 |
| 154 | `ARID5A` | CRA | 2 | 179 | 70 | 1157 |
| 155 | `SIM1` | CRA | 2 | 179 | 107 | 1710 |
| 156 | `FOXD2` | CRA | 2 | 178 | 84 | 1723 |
| 157 | `ZBTB7A` | CRI | 2 | 177 | 47 | 1316 |
| 158 | `GABARAPL1` | CRA | 2 | 176 | 91 | 2441 |
| 159 | `ATG3` | CRA | 2 | 175 | 98 | 1638 |
| 160 | `DUXA` | CRI | 2 | 175 | 62 | 1280 |
| 161 | `SREBF2` | CRI | 2 | 175 | 46 | 701 |
| 162 | `ZBTB7B` | CRI | 2 | 175 | 50 | 1685 |
| 163 | `XRCC5` | CRA | 2 | 172 | 101 | 2316 |
| 164 | `FOXC1` | CRA | 2 | 171 | 100 | 1759 |
| 165 | `DRGX` | CRI | 2 | 170 | 36 | 1408 |
| 166 | `TFDP1` | CRI | 2 | 170 | 59 | 1363 |
| 167 | `NFIL3` | CRI | 2 | 168 | 44 | 1501 |
| 168 | `NKX2-2` | CRI | 2 | 168 | 48 | 1977 |
| 169 | `PPARG` | CRI | 2 | 167 | 54 | 1417 |
| 170 | `BATF3` | CRA | 2 | 166 | 96 | 1238 |
| 171 | `FOXO3` | CRI | 2 | 166 | 28 | 664 |
| 172 | `TFE3` | CRA | 2 | 162 | 82 | 3447 |
| 173 | `HOXA1` | CRA | 2 | 161 | 96 | 2383 |
| 174 | `ZIC3` | CRA | 2 | 160 | 85 | 1951 |
| 175 | `HSF2` | CRI | 2 | 159 | 50 | 1389 |
| 176 | `HOXB13` | CRI | 2 | 158 | 54 | 1204 |
| 177 | `ZKSCAN4` | CRI | 2 | 158 | 53 | 1445 |
| 178 | `ZNF384` | CRI | 2 | 158 | 59 | 1388 |
| 179 | `MXD3` | CRI | 2 | 157 | 50 | 738 |
| 180 | `NR0B1` | CRA | 2 | 157 | 84 | 3809 |
| 181 | `POU2F1` | CRA | 2 | 157 | 80 | 1339 |
| 182 | `EGR2` | CRI | 2 | 156 | 57 | 1643 |
| 183 | `ZNF318` | CRA | 2 | 156 | 105 | 1588 |
| 184 | `CEBPA` | CRI | 2 | 153 | 50 | 784 |
| 185 | `FOXI1` | CRI | 2 | 153 | 50 | 1239 |
| 186 | `FOXC1` | CRI | 2 | 152 | 42 | 1159 |
| 187 | `HOXD11` | CRI | 2 | 152 | 54 | 1544 |
| 188 | `JUN` | CRI | 2 | 152 | 39 | 935 |
| 189 | `SRY` | CRI | 2 | 152 | 42 | 823 |
| 190 | `HOXB9` | CRA | 2 | 151 | 80 | 2310 |
| 191 | `SIM2` | CRI | 2 | 151 | 52 | 1214 |
| 192 | `ATG4D` | CRI | 2 | 150 | 40 | 897 |
| 193 | `FERD3L` | CRA | 2 | 150 | 85 | 2658 |
| 194 | `HOXC5` | CRA | 2 | 150 | 70 | 990 |
| 195 | `NR4A1` | CRI | 2 | 150 | 64 | 1549 |
| 196 | `PAX8` | CRA | 2 | 149 | 93 | 1427 |
| 197 | `POU5F1` | CRA | 2 | 149 | 62 | 1497 |
| 198 | `ZBTB2` | CRA | 2 | 149 | 67 | 1588 |
| 199 | `NFATC2` | CRI | 2 | 148 | 48 | 770 |
| 200 | `CEBPZ` | CRA | 2 | 147 | 71 | 1539 |
| 201 | `ETS2` | CRI | 2 | 147 | 49 | 1354 |
| 202 | `HOXA13` | CRI | 2 | 147 | 38 | 909 |
| 203 | `MEF2C` | CRA | 2 | 147 | 71 | 1953 |
| 204 | `ETV5` | CRA | 2 | 144 | 71 | 980 |
| 205 | `SOX10` | CRA | 2 | 144 | 83 | 1972 |
| 206 | `CEBPB` | CRA | 2 | 143 | 57 | 415 |
| 207 | `MITF` | CRI | 2 | 143 | 53 | 870 |
| 208 | `MYB` | CRI | 2 | 142 | 45 | 928 |
| 209 | `TFCP2` | CRI | 2 | 142 | 43 | 1109 |
| 210 | `ZNF263` | CRI | 2 | 141 | 57 | 587 |
| 211 | `ZNF35` | CRA | 2 | 141 | 89 | 2126 |
| 212 | `PRKRIR` | CRA | 2 | 140 | 70 | 1327 |
| 213 | `FOXD2` | CRI | 2 | 138 | 64 | 1929 |
| 214 | `TWIST1` | CRA | 2 | 138 | 64 | 1368 |
| 215 | `DLX2` | CRI | 2 | 137 | 29 | 609 |
| 216 | `HMGB2` | CRI | 2 | 137 | 43 | 1215 |
| 217 | `FOSL1` | CRI | 2 | 136 | 55 | 1544 |
| 218 | `GATA1` | CRI | 2 | 136 | 62 | 1137 |
| 219 | `NFKB2` | CRI | 2 | 135 | 42 | 1022 |
| 220 | `ZNF219` | CRI | 2 | 135 | 36 | 564 |
| 221 | `ATG4D` | CRA | 2 | 133 | 55 | 1245 |
| 222 | `LCOR` | CRA | 2 | 133 | 64 | 1231 |
| 223 | `GSX2` | CRA | 2 | 132 | 53 | 1089 |
| 224 | `NOBOX` | CRA | 2 | 132 | 74 | 1600 |
| 225 | `SALL4` | CRI | 2 | 132 | 48 | 1044 |
| 226 | `SOX9` | CRA | 2 | 132 | 66 | 1779 |
| 227 | `TFAP2A` | CRI | 2 | 132 | 35 | 916 |
| 228 | `ZKSCAN3` | CRA | 2 | 132 | 61 | 1948 |
| 229 | `HOXA1` | CRI | 2 | 131 | 34 | 728 |
| 230 | `NPAS1` | CRI | 2 | 131 | 46 | 968 |
| 231 | `ESX1` | CRI | 2 | 130 | 48 | 1149 |
| 232 | `PAX2` | CRA | 2 | 130 | 71 | 1258 |
| 233 | `SALL4` | CRA | 2 | 130 | 72 | 2344 |
| 234 | `CTCF` | CRA | 2 | 129 | 73 | 1318 |
| 235 | `HINFP` | CRI | 2 | 129 | 17 | 719 |
| 236 | `HOXA2` | CRA | 2 | 129 | 81 | 2426 |
| 237 | `FOXO4` | CRA | 2 | 128 | 76 | 2818 |
| 238 | `ZNF695` | CRA | 2 | 128 | 57 | 1379 |
| 239 | `SIRT1` | CRI | 2 | 127 | 48 | 1008 |
| 240 | `ATG5` | CRI | 2 | 126 | 67 | 1150 |
| 241 | `BRCA1` | CRA | 2 | 125 | 77 | 1464 |
| 242 | `FOXA2` | CRA | 2 | 125 | 58 | 722 |
| 243 | `ZNF367` | CRI | 2 | 125 | 26 | 597 |
| 244 | `NR4A2` | CRA | 2 | 124 | 63 | 1358 |
| 245 | `SOX10` | CRI | 2 | 124 | 32 | 731 |
| 246 | `CREB1` | CRA | 2 | 123 | 62 | 1384 |
| 247 | `IRF4` | CRA | 2 | 123 | 90 | 1498 |
| 248 | `DLX4` | CRA | 2 | 122 | 64 | 1641 |
| 249 | `MYB` | CRA | 2 | 122 | 63 | 1156 |
| 250 | `PLAU` | CRA | 2 | 122 | 54 | 572 |
| 251 | `SOX7` | CRA | 2 | 122 | 60 | 577 |
| 252 | `NOBOX` | CRI | 2 | 121 | 34 | 1174 |
| 253 | `GATA1` | CRA | 2 | 120 | 58 | 1538 |
| 254 | `HOXB9` | CRI | 2 | 120 | 32 | 763 |
| 255 | `LHX1` | CRA | 2 | 120 | 62 | 1081 |
| 256 | `VAX1` | CRA | 2 | 120 | 80 | 1102 |
| 257 | `CEBPA` | CRA | 2 | 119 | 66 | 1099 |
| 258 | `HSF1` | CRA | 2 | 119 | 49 | 1153 |
| 259 | `HSF4` | CRI | 2 | 118 | 34 | 974 |
| 260 | `MYBL2` | CRA | 2 | 118 | 68 | 1176 |
| 261 | `FOXO3` | CRA | 2 | 117 | 48 | 1005 |
| 262 | `NPAS3` | CRA | 2 | 117 | 68 | 1590 |
| 263 | `TFAM` | CRA | 2 | 117 | 58 | 1244 |
| 264 | `CDX2` | CRI | 2 | 116 | 47 | 868 |
| 265 | `IRF4` | CRI | 2 | 116 | 47 | 1011 |
| 266 | `MITF` | CRA | 2 | 114 | 48 | 974 |
| 267 | `MRRF` | CRA | 2 | 113 | 63 | 1014 |
| 268 | `PARP1` | CRA | 2 | 113 | 73 | 1335 |
| 269 | `TEAD3` | CRI | 2 | 113 | 35 | 645 |
| 270 | `ATG5` | CRA | 2 | 112 | 63 | 1001 |
| 271 | `KLF3` | CRA | 2 | 112 | 62 | 1052 |
| 272 | `NFKB1` | CRA | 2 | 112 | 50 | 2292 |
| 273 | `ZNF695` | CRI | 2 | 112 | 42 | 1657 |
| 274 | `CTCFL` | CRA | 2 | 111 | 59 | 1081 |
| 275 | `FOXM1` | CRI | 2 | 111 | 32 | 603 |
| 276 | `NKX1-1` | CRI | 2 | 111 | 41 | 341 |
| 277 | `TFEC` | CRI | 2 | 111 | 34 | 698 |
| 278 | `FOXI1` | CRA | 2 | 110 | 73 | 1130 |
| 279 | `NFATC4` | CRA | 2 | 110 | 38 | 494 |
| 280 | `PPARG` | CRA | 2 | 109 | 64 | 905 |
| 281 | `PPARGC1A` | CRI | 2 | 109 | 22 | 674 |
| 282 | `CDX2` | CRA | 2 | 108 | 50 | 1184 |
| 283 | `ZKSCAN4` | CRA | 2 | 108 | 41 | 470 |
| 284 | `HOXB7` | CRI | 2 | 107 | 27 | 592 |
| 285 | `SIRT1` | CRA | 2 | 107 | 54 | 970 |
| 286 | `XBP1` | CRA | 2 | 107 | 44 | 718 |
| 287 | `E2F1` | CRA | 2 | 105 | 54 | 879 |
| 288 | `E2F4` | CRA | 2 | 105 | 67 | 1408 |
| 289 | `MAZ` | CRA | 2 | 105 | 53 | 1071 |
| 290 | `POU3F2` | CRA | 2 | 105 | 63 | 996 |
| 291 | `ZFX` | CRA | 2 | 104 | 51 | 822 |
| 292 | `NRF1` | CRA | 2 | 103 | 52 | 1249 |
| 293 | `TFAP2A` | CRA | 2 | 103 | 59 | 1568 |
| 294 | `RXRB` | CRI | 2 | 102 | 39 | 452 |
| 295 | `STAT1` | CRA | 2 | 102 | 51 | 1155 |
| 296 | `CEBPZ` | CRI | 2 | 101 | 20 | 538 |
| 297 | `HSPA8` | CRA | 2 | 101 | 56 | 552 |
| 298 | `BATF3` | CRI | 2 | 100 | 32 | 703 |
| 299 | `ZNF281` | CRA | 2 | 100 | 61 | 1338 |
| 300 | `E2F7` | CRA | 2 | 99 | 46 | 1266 |
| 301 | `TFEB` | CRA | 2 | 99 | 62 | 981 |
| 302 | `SOX11` | CRA | 2 | 98 | 57 | 663 |
| 303 | `ALX1` | CRA | 2 | 97 | 58 | 1415 |
| 304 | `ZIC5` | CRI | 2 | 97 | 8 | 661 |
| 305 | `ATF4` | CRA | 2 | 96 | 45 | 461 |
| 306 | `TFDP1` | CRA | 2 | 96 | 50 | 1089 |
| 307 | `DNMT1` | CRA | 2 | 95 | 53 | 1204 |
| 308 | `NFIA` | CRI | 2 | 95 | 51 | 760 |
| 309 | `PPARGC1A` | CRA | 2 | 95 | 44 | 704 |
| 310 | `BSX` | CRA | 2 | 93 | 52 | 1399 |
| 311 | `EN1` | CRI | 2 | 93 | 19 | 636 |
| 312 | `IRF5` | CRA | 2 | 93 | 43 | 716 |
| 313 | `IRF7` | CRA | 2 | 92 | 46 | 806 |
| 314 | `DLX5` | CRA | 2 | 91 | 41 | 609 |
| 315 | `HNF1B` | CRI | 2 | 91 | 44 | 755 |
| 316 | `ZBTB7B` | CRA | 2 | 90 | 22 | 811 |
| 317 | `THRA` | CRA | 2 | 89 | 51 | 482 |
| 318 | `RXRB` | CRA | 2 | 88 | 36 | 937 |
| 319 | `E2F5` | CRA | 2 | 86 | 54 | 1244 |
| 320 | `ONECUT1` | CRA | 2 | 86 | 43 | 2031 |
| 321 | `BRCA1` | CRI | 2 | 84 | 27 | 458 |
| 322 | `FLI1` | CRA | 2 | 84 | 36 | 327 |
| 323 | `ATF7` | CRA | 2 | 83 | 50 | 786 |
| 324 | `HMGB1` | CRA | 2 | 83 | 43 | 1377 |
| 325 | `RUNX2` | CRI | 2 | 83 | 29 | 505 |
| 326 | `EGR1` | CRA | 2 | 82 | 35 | 939 |
| 327 | `FOXM1` | CRA | 2 | 82 | 38 | 264 |
| 328 | `ZIC3` | CRI | 2 | 82 | 14 | 362 |
| 329 | `FOXF2` | CRA | 2 | 79 | 35 | 763 |
| 330 | `EOMES` | CRA | 2 | 78 | 30 | 866 |
| 331 | `HOXC10` | CRI | 2 | 78 | 28 | 627 |
| 332 | `MEF2C` | CRI | 2 | 78 | 8 | 419 |
| 333 | `NFATC4` | CRI | 2 | 78 | 25 | 316 |
| 334 | `ETS1` | CRI | 2 | 77 | 33 | 516 |
| 335 | `NR0B1` | CRI | 2 | 77 | 32 | 812 |
| 336 | `SP1` | CRA | 2 | 77 | 58 | 335 |
| 337 | `ZNF367` | CRA | 2 | 77 | 35 | 551 |
| 338 | `FOXA2` | CRI | 2 | 76 | 29 | 761 |
| 339 | `HSPD1` | CRI | 2 | 75 | 16 | 576 |
| 340 | `HNF1B` | CRA | 2 | 73 | 44 | 764 |
| 341 | `RUNX2` | CRA | 2 | 73 | 27 | 922 |
| 342 | `KLF6` | CRA | 2 | 71 | 38 | 1242 |
| 343 | `ZFHX3` | CRA | 2 | 71 | 37 | 894 |
| 344 | `STAT4` | CRA | 2 | 70 | 43 | 465 |
| 345 | `TCF3` | CRA | 2 | 70 | 40 | 326 |
| 346 | `ATG4C` | CRA | 2 | 69 | 38 | 1087 |
| 347 | `EGR2` | CRA | 2 | 69 | 29 | 448 |
| 348 | `APOD` | CRA | 2 | 68 | 24 | 615 |
| 349 | `ZNF219` | CRA | 2 | 68 | 34 | 530 |
| 350 | `HOXC10` | CRA | 2 | 67 | 33 | 299 |
| 351 | `JUN` | CRA | 2 | 67 | 29 | 702 |
| 352 | `PBX2` | CRA | 2 | 67 | 32 | 893 |
| 353 | `ZBTB7C` | CRA | 2 | 67 | 36 | 697 |
| 354 | `EN1` | CRA | 2 | 66 | 30 | 296 |
| 355 | `ETS1` | CRA | 2 | 66 | 31 | 1189 |
| 356 | `SRY` | CRA | 2 | 66 | 37 | 886 |
| 357 | `HSF4` | CRA | 2 | 65 | 39 | 240 |
| 358 | `NFIL3` | CRA | 2 | 65 | 36 | 371 |
| 359 | `HINFP` | CRA | 2 | 64 | 40 | 754 |
| 360 | `TCF7L2` | CRI | 2 | 64 | 15 | 475 |
| 361 | `HOXD11` | CRA | 2 | 62 | 30 | 499 |
| 362 | `HOXD3` | CRI | 2 | 61 | 8 | 188 |
| 363 | `HSPD1` | CRA | 2 | 61 | 27 | 634 |
| 364 | `NR4A1` | CRA | 2 | 61 | 24 | 611 |
| 365 | `ZNF263` | CRA | 2 | 61 | 37 | 874 |
| 366 | `FOXS1` | CRA | 2 | 60 | 25 | 355 |
| 367 | `HOXB4` | CRA | 2 | 60 | 27 | 299 |
| 368 | `HMGB2` | CRA | 2 | 59 | 36 | 842 |
| 369 | `DRGX` | CRA | 2 | 58 | 33 | 282 |
| 370 | `STAT3` | CRA | 2 | 58 | 25 | 255 |
| 371 | `GATA2` | CRA | 2 | 56 | 25 | 350 |
| 372 | `SOX5` | CRA | 2 | 56 | 25 | 348 |
| 373 | `KLF4` | CRA | 2 | 54 | 21 | 333 |
| 374 | `FOXP1` | CRA | 2 | 52 | 30 | 468 |
| 375 | `NFE2L2` | CRA | 2 | 52 | 21 | 439 |
| 376 | `TBX3` | CRA | 2 | 51 | 30 | 170 |
| 377 | `PHOX2B` | CRA | 2 | 50 | 34 | 666 |
| 378 | `ATF6` | CRA | 2 | 48 | 27 | 173 |
| 379 | `ESR1` | CRI | 2 | 48 | 11 | 372 |
| 380 | `SOX17` | CRA | 2 | 48 | 19 | 156 |
| 381 | `ZBTB7A` | CRA | 2 | 46 | 18 | 283 |
| 382 | `FOXP2` | CRA | 2 | 45 | 23 | 574 |
| 383 | `SOX2` | CRA | 2 | 45 | 20 | 740 |
| 384 | `CTCF` | CRI | 2 | 42 | 20 | 96 |
| 385 | `NFIA` | CRA | 2 | 40 | 21 | 372 |
| 386 | `DUXA` | CRA | 2 | 39 | 26 | 338 |
| 387 | `DLX5` | CRI | 1 | 37 | 19 | 377 |
| 388 | `ESX1` | CRA | 2 | 36 | 13 | 502 |
| 389 | `TFEC` | CRA | 2 | 33 | 11 | 103 |
| 390 | `TCF7L2` | CRA | 2 | 27 | 12 | 184 |
| 391 | `POU2F1` | CRI | 2 | 25 | 5 | 73 |
| 392 | `ZIC5` | CRA | 2 | 23 | 10 | 224 |
| 393 | `ESR1` | CRA | 2 | 20 | 4 | 52 |
| 394 | `NFATC2` | CRA | 2 | 16 | 3 | 16 |
| 395 | `NPAS1` | CRA | 2 | 14 | 6 | 108 |
| 396 | `SREBF2` | CRA | 2 | 11 | 8 | 99 |
| 397 | `ATF6` | CRI | 2 | 8 | 2 | 9 |
| 398 | `PHOX2A` | CRA | 1 | 5 | 3 | 23 |
| 399 | `HOXA13` | CRA | 2 | 2 | 0 | 6 |
| 400 | `HSPA8` | CRI | 2 | 2 | 2 | 3 |

### Guide-UMI columns per cell

- number of guide columns with UMI>0 per cell: percentiles 0/25/50/75/90/99/100 = 0 / 1 / 2 / 3 / 4 / 21 / 46
- total guide UMI per cell: percentiles 0/25/50/75/99/100 = 0 / 1 / 8 / 25 / 119.87 / 486

## 5. CRISPRa vs CRISPRi

- **both are in this one object.**
- direction is encoded in two places, which pick out different cell sets
  (tables below):
  1. in the guide column **name**: the `CRA` / `CRI` infix.
  2. as two non-guide `.obs` columns, `CRA` (stored int64) and `CRI` (int64).

| column | n cells nonzero | n distinct | min nonzero | median nonzero | max |
|---|---|---|---|---|---|
| `CRA` | 14333 | 2 | 1 | 1 | 1 |
| `CRI` | 15372 | 2 | 1 | 1 | 1 |

| `CRA`>0 | `CRI`>0 | n cells |
|---|---|---|
| False | False | 6809 |
| False | True | 15372 |
| True | False | 14333 |
| True | True | 0 |

### Counts per direction taken from the guide columns instead

- CRA: 399 guide columns, 200 distinct target names, 14187 cells with any UMI>0, total UMI 278254
- CRI: 399 guide columns, 200 distinct target names, 15060 cells with any UMI>0, total UMI 362038

| any CRA guide UMI>0 | any CRI guide UMI>0 | n cells |
|---|---|---|
| False | False | 7991 |
| False | True | 14336 |
| True | False | 13463 |
| True | True | 724 |

- cells assigned to a CRA guide by convention B (argmax): 12646
- cells assigned to a CRI guide by convention B (argmax): 13253

**Recorded as an ambiguity.** The `CRA`/`CRI` `.obs` columns and the
`_CRA`/`_CRI` guide-column names give different cell sets (tables above),
and the file carries no statement of what the two `.obs` columns count or
which encoding is authoritative.

### Direction against SampleName and the PD columns

| group | `UNKNOWN` | `PD32` | `PD26` | `PD14` | `PD14`>0 | `PD26`>0 | `PD32`>0 | total |
|---|---|---|---|---|---|---|---|---|
| `CRA`>0 | 14333 | 0 | 0 | 0 | 0 | 0 | 0 | 14333 |
| `CRI`>0 | 15372 | 0 | 0 | 0 | 0 | 0 | 0 | 15372 |
| any CRA guide UMI>0 | 14187 | 0 | 0 | 0 | 0 | 0 | 0 | 14187 |
| any CRI guide UMI>0 | 15060 | 0 | 0 | 0 | 0 | 0 | 0 | 15060 |

## 6. E2F3, EZH2, STAT3, ZFX

| factor | guide columns present | direction(s) | cells any UMI>0 | cells argmax | total guide UMI |
|---|---|---|---|---|---|
| `E2F3` | `E2F3_CRA1`, `E2F3_CRA2`, `E2F3_CRI1`, `E2F3_CRI2` | CRA, CRI | 563 | 207 | 4951 |
| `EZH2` | `EZH2_CRA1`, `EZH2_CRA2`, `EZH2_CRI1`, `EZH2_CRI2` | CRA, CRI | 501 | 202 | 7977 |
| `STAT3` | `STAT3_CRA1`, `STAT3_CRA2`, `STAT3_CRI1`, `STAT3_CRI2` | CRA, CRI | 276 | 126 | 2515 |
| `ZFX` | `ZFX_CRA1`, `ZFX_CRA2`, `ZFX_CRI1`, `ZFX_CRI2` | CRA, CRI | 626 | 214 | 5469 |

Per direction:

| factor | direction | guide columns | cells any UMI>0 | cells argmax |
|---|---|---|---|---|
| `E2F3` | CRA | `E2F3_CRA1`, `E2F3_CRA2` | 254 | 128 |
| `E2F3` | CRI | `E2F3_CRI1`, `E2F3_CRI2` | 309 | 79 |
| `EZH2` | CRA | `EZH2_CRA1`, `EZH2_CRA2` | 221 | 109 |
| `EZH2` | CRI | `EZH2_CRI2`, `EZH2_CRI1` | 281 | 93 |
| `STAT3` | CRA | `STAT3_CRA1`, `STAT3_CRA2` | 58 | 25 |
| `STAT3` | CRI | `STAT3_CRI2`, `STAT3_CRI1` | 218 | 101 |
| `ZFX` | CRA | `ZFX_CRA2`, `ZFX_CRA1` | 104 | 51 |
| `ZFX` | CRI | `ZFX_CRI2`, `ZFX_CRI1` | 522 | 163 |

Both counts are under the stated conventions of section 4, not an
assignment the object makes. Section 9 reports these four as `.var` symbols.

## 7. Cell-cycle phase annotation

- `.obs` columns matching /phase|cycle|S_score|G2M|G1S|cc_/i: NONE
- no `phase`, `S_score` or `G2M_score` column exists; `.uns` is absent, so
  no cell-cycle gene list is stored either.
- **absent.** There is no cell-cycle phase annotation in this object.

## 8. Author-computed scores already in the object

- `.obs` columns matching /R_rej|rej|score|age|clock|senes|pseudo|umap|pca|harmony|leiden|louvain|cluster/i: NONE
- `.obsm` is absent, so no embedding (PCA, UMAP, Harmony) is stored.
- `.uns` is absent, so no author parameters, rankings or score tables are
  stored.
- the only derived per-cell quantities in `.obs` are the three QC columns
  `n_counts`, `log_counts`, `n_genes`, summarised in 1b.
- **no `R_rej` and no comparable author-computed score is present.**

## 9. Gene identifier format, and overlap with the frozen GTEx panel

- `.var` index length 36,601; first 10 entries: ['MIR1302-2HG', 'FAM138A', 'OR4F5', 'AL627309.1', 'AL627309.3', 'AL627309.2', 'AL627309.5', 'AL627309.4', 'AP006222.2', 'AL732372.1']
- entries matching `^ENSG\d{11}`: 0
- format: **gene symbols**, not Ensembl IDs.
- duplicated `.var` index entries: 0 (n unique = 36,601)
- symbols appearing more than once: 0

- reference panel **read, not refit**, from
  `results/fibro/frozen_ruler_ridge_raw.npz` (the frozen GTEx fibroblast
  ruler of this project), keys `symbol`, `ensembl`, `gene_id`; n = 23,485.

| panel key | panel rows matched in `.var` | panel rows missing | distinct ids in common |
|---|---|---|---|
| `symbol` | 16,781 | 6,704 | 16,771 |
| `ensembl` | 0 | 23,485 | 0 |
| `gene_id` | 0 | 23,485 | 0 |

- **answer to the question as asked** (panel `symbol` against the `.var`
  index, plain string match): of the 23,485 panel genes, **16,781 are present**
  in this object and **6,704 are missing**.
- the panel's own `symbol` column has 18 duplicated values; distinct panel
  symbols = 23,467, of which 16,771 are present.
- `.var` symbols not on the panel: 19,830 of 36,601 distinct.
- written: `panel_symbols_missing.csv`, `panel_symbols_present.csv`,
  `var_index.csv` (all under `results/senginv/`).
- first 30 missing panel symbols: ['AARS1', 'AARSD1P1', 'ABCA11P', 'ABCA17P', 'ABCC6P1', 'ABCC6P2', 'ABHD11-AS1', 'ABHD17AP1', 'ACP3', 'ACTBP11', 'ACTBP2', 'ACTBP8', 'ACTBP9', 'ACTG1P1', 'ACTG1P10', 'ACTG1P14', 'ACTG1P15', 'ACTG1P17', 'ACTG1P20', 'ACTG1P24', 'ACTG1P25', 'ACTG1P3', 'ACTG1P9', 'ACTN1-DT', 'ACTN4P1', 'ACTR3-AS1', 'ADAM1B', 'ADAM20P1', 'ADAM21P1', 'ADAMTS7P1']

**Recorded as an ambiguity.** The match above is a plain string match on
gene symbols. No alias resolution and no Ensembl-mediated remapping was
attempted. The panel carries Ensembl IDs and this file does not, so the
`ensembl` and `gene_id` rows above are zero by construction rather than by
disagreement, and the symbol row is the only one that can be read. The
`.var` index has 0 duplicated entries and the panel `symbol` column has 18,
so 16781 panel rows can map onto fewer distinct `.var` columns. Nothing was
deduplicated and no decision about the duplicates was made.

### The guide-column target prefixes as `.var` symbols

- of the 200 distinct guide target prefixes, 199 appear as a symbol in `.var` and 1 do not.
- prefixes absent from `.var`: ['PRKRIR']

### The four named factors as genes in `.var`

| factor | in `.var` | n `.var` columns with that symbol | on the frozen panel |
|---|---|---|---|
| `E2F3` | True | 1 | True |
| `EZH2` | True | 1 | True |
| `STAT3` | True | 1 | True |
| `ZFX` | True | 1 | True |

## Files written by this run

- `results/senginv/cells_per_target.csv` (10,579 B)
- `results/senginv/figshare_30898748.json` (1,117 B)
- `results/senginv/figshare_30899018.json` (1,156 B)
- `results/senginv/guide_columns.csv` (43,331 B)
- `results/senginv/INVENTORY_SENGINV.md` (44,077 B)
- `results/senginv/md5_local.txt` (33 B)
- `results/senginv/panel_symbols_missing.csv` (95,502 B)
- `results/senginv/panel_symbols_present.csv` (128,333 B)
- `results/senginv/run_stderr.txt` (0 B)
- `results/senginv/run_stdout.txt` (41,003 B)
- `results/senginv/var_index.csv` (344,572 B)

Nothing was committed and nothing was pushed. No pre-registration was
written, no statistic was computed on `.X`, and no analysis plan is proposed.
