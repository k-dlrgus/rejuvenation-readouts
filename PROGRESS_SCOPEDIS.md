# PROGRESS_SCOPEDIS

**STOP status:** complete. Neither `bar_ambiguous` (S1) nor `reproduction_failed` (S2) fired as a
task STOP. Two carried-bar ambiguities did surface and were **put to the user mid-task rather than
resolved by the session** (U4, U5 in the flag's addendum); the user ruled, and the run continued on
their rulings. Area 1 reached `cap_reached` (U3), which closes one area and is not a STOP.

**Outcome:** the fibrosis gap is general in four areas of five and fails in one.

- `GO` fired three times, all in area 1 (immune-cell activation): **GSE314342**, **GSE278572**,
  **GSE92872**. `INCOMPLETE` fired for **GSE249595** (bar C unverified).
- `NO-GO` fired separately in areas 2 (genetic disease), 3 (senescence), 4 (cancer) and 5
  (neurodegeneration). It does **not** fire overall.
- `PATTERN` (descriptive): of 8 screens with ≥20 perturbations of a diseased/activated state across
  5 disease areas, 4 have the control state `in_run`, 1 `separate_experiment`, 3 `absent`, 0
  `unstated`. Human-only: 4 / 1 / 1 / 0. All 8 are in area 1, the capped area; areas 2–5 contributed
  none.
- Why area 1 differs is mechanical, not biological: activation is induced and reversible, so the
  resting arm is the same cells held back, and the field multiplexes conditions into one lane (10x
  Flex probe barcodes, TotalSeq hashtags, cell hashing). Once conditions share a lane the control
  state is in the run by construction and C3 `split` cannot arise.
- The price, carried in the GO readings themselves per U4: in both primary-cell GO datasets the
  in-run target is `returned` (activated once, then rested), not `naive`.

**Single next action:** none in this task. Whether a `returned` resting arm is an acceptable target
for a distance measurement, and what to do about both GO datasets having `before` timing (C1 makes
that inventory-only for the follow-up), are decisions outside a scoping session.

## Record of what was decided when (times +09:00, 2026-10-24)

| time | event |
|---|---|
| before 11:28 | Read CLAUDE.md, `results/scopefib/PREREG.flag`, `SCOPING_FIBROSIS.md`, `SCOPING_FIBROSIS2.md`, `PROGRESS_SCOPEFIB2.md`, part of `src/scopefib2_run.py`, `.gitignore`, `git log -1`. Repository files only, no query. |
| ~11:25 | Three ambiguities put to the user before the flag was written; rulings U1 (area 4 target), U2 (knock-in pairs), U3 (the 8-record cap) recorded verbatim in the flag. |
| 11:29:12 | `results/scopedis/PREREG.flag` written, with lines 10–87 of the SCOPEFIB flag spliced in byte-for-byte by script. Not edited afterwards; U4/U5 were later **appended** as a timestamped addendum. |
| ~11:31 | Pre-registered reproduction run: 12 checks, 0 mismatches → `reproduction.csv`. |
| ~11:32 | First queries of the session (GEO DataSets, area 1). |
| ~11:33 | `candidates_1.csv` frozen, then records opened in stated-count order. |
| ~11:40 | GSE314342 resolved: resting and stimulated multiplexed in the same 10x Flex lane. First `in_run`. |
| ~11:45 | **Stopped and asked.** GSE314342's resting arm had been activated once on day 0 (D5 ambiguous) and GSE249595 reads out a 374-gene panel (bar B's "transcriptomic readout", a question SCOPEFIB left open). Both could decide a reading, so neither was resolved by the session. |
| ~11:48 | Rulings U4 and U5 received and appended to the flag as a disclosed addendum, before any further record was opened. Mandatory columns `target_is_naive` and `readout_scale` added, plus the `insufficient_for_distance` flag at U5's 2,000-gene floor. |
| ~11:50–12:05 | Area 1 completed to its 8-record cap. GSE92872's library size and per-arm median cells per gene resolved by parsing the header rows of two of its own supplementary files (1.55 MB and 17.3 MB, both under the limit). |
| ~12:05–12:20 | Areas 2–5 searched and their records opened. `candidates_2.csv` frozen; areas 3–5 candidate files were **not** written before opening (disclosed as a deviation from D9). |
| ~12:25 | `results/scopedis/_tables.py` → `datasets.csv` (30 rows); `SCOPING_DISEASE.md` written with the flag embedded verbatim; this file. |

## Queries run (do not repeat)

All stored with term, hit count and summaries in `results/scopedis/searches.json`; opened records in
`geo_records.json` and `url_records.json`.

- **GEO DataSets (`db=gds`)**: 13 queries for area 1, 13 for area 2, 9 for area 3, 8 for area 4, 9
  for area 5, each capped at 30–40 esummary records. Search listings do not count as opened records
  (D8).
- **Europe PMC**: 3 queries for area 1, 2 for area 2, plus single-record lookups for
  EXT_ID 40011559, 39663454, 42062284, the CROP-seq 2017 title, and DOI 10.64898/2025.12.23.696273.
- **Full texts read**: PMC11906366 (PS/GSE249595), PMC11754113 (IL2RA/GSE278572), PMC5334791
  (CROP-seq 2017/GSE92872), bioRxiv `10.64898/2025.12.23.696273v1.full` (GSE314342). The journal
  version of the last (doi 10.1016/j.cell.2026.08.002) was not fetched; the preprint was.
- **Failed by design, logged**: `europepmc .../MED/42664972/fullTextXML` 404,
  `.../PPR/PPR1230028/fullTextXML` 404.
- **Supplementary files opened**: `GSE92872_CROP-seq_Jurkat_TCR.count_matrix.csv.gz` (1,621,695
  bytes) and `GSE92872_CROP-seq_Jurkat_TCR.digital_expression.csv.gz` (18,181,357 bytes), header
  rows only, saved under `results/scopedis/opened/gse92872/`.
- Nothing over 50 MB was downloaded. No controlled-access request was submitted. No statistic was
  computed; the only arithmetic was counting, one median per arm, and two means labelled estimates.

## Open items (not resolved in this task)

- **Area 1's cap.** Four candidates on its frozen list were never opened: GSE203593, GSE119450,
  GSE190604, GSE137554. Area 1 is the area that decides the overall verdict, so the GO is a
  statement about records opened, not about the area.
- **Bar C for GSE249595.** Achieved cells per gene are nowhere public; the per-cell guide and hash
  assignment in `GSE249595_RAW.tar` would settle it. The post-QC total divides to ~28 cells per gene
  in the activated arm, which is below D7's 30 — so this could resolve either way.
- **Distinct-agent counts not in four records**: GSE221578, GSE305037 (186 sgRNAs, genes unstated),
  GSE213511, GSE272093. Each fails bar A independently, so they cannot change a reading, but they
  would change the PATTERN denominator.
- **D9 was not followed for areas 3–5**: their candidate lists were not frozen to disk before
  opening began. The ordering is reconstructible from `searches.json` but was not pre-committed.
- **Two definitional questions SCOPEFIB left open, one now closed.** U5 closes the targeted-panel
  question and this file applies it retrospectively to GSE227391 without editing
  `SCOPING_FIBROSIS.md`. Whether CRISPRi/a of regulatory **elements** counts as a D1 agent is still
  open; it bears on GSE309511 here and on GSE243756 there, and neither can pass bar A anyway.
- **Repositories not searched**: ArrayExpress/BioStudies, GSA, EGA, CELLxGENE, Single Cell Portal,
  Zenodo, figshare. Only GEO and Europe PMC were queried in this task.
- **Aged cells.** No GO dataset has an aged donor: GSE314342's donors are 22–34, GSE278572's two
  donors are of unstated age, and GSE92872 / GSE249595 are a leukaemia line. Nothing here transfers
  to aged cells.
