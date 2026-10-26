# WITHDRAWN_READINGS — every withdrawn or superseded reading in the repository

Written 2026-10-02. Read-only search; nothing was re-run and no source file was edited. Searched every project file outside `.venv/`, `bundle/`, `data/raw/` and `data/processed/` for `withdraw*`, `supersed*`, `retract*`, `retire*`, `replaced by`, `no longer`, `correction`, `stale`, `obsolete` and `overturn`. Copies under `paper_package/` repeat the root files and are not listed separately. Third-party text (saved GEO pages, gene-set catalogs) was dropped.

preprint_10 paragraph P89 (1-indexed top-level body paragraphs, the numbering used by `results/verify/claims_v9.csv`) says:

> Where a later analysis superseded an earlier reading, the earlier text was retained and marked withdrawn rather than edited; three readings were withdrawn in this way.

`README.md` lines 145–146 repeat the same count.

> **Updated 2026-10-02 (preprint_11 fixes).** The atlas ruler − MD Δρ on 93 pseudobulks, listed as C1 in the first version of this file, is now marked superseded in the repository. A dated note was appended to every ledger or FINDINGS file that said MD3 Task 2 "stands", and to the file holding the reading. It is now B7, and section C is empty. The counts below were updated to match. preprint_11 P89 now gives these counts (see "Counts used by preprint_11 P89" at the end). The rest of the file is as first written, except where a sentence says it was updated.

> **Updated 2026-10-26 (CORR1).** A statement in paper 1 v1.2, Supplementary Note 3, was corrected: the transcription-factor screen's length statement in space C. The correction notice is `CORRECTIONS_PAPER1.md`; the check behind it is `FINDINGS_FLOORCHK.md`. The statement is descriptive, not a gated reading, so it is listed in section D and the counts below (nine, three) do not change. Nothing else in this file was changed.

## Result in one paragraph

The repository does not record three. It labels **two** readings *withdrawn* (SOUTH4's stop and MD3 Task 1) and, since the 2026-10-02 update, **seven** *superseded* (SOUTH, SOUTH2, SOUTH3, SOUTH5, MD2 Task 2, POSCTRL P1 and the atlas ruler − MD Δρ on 93 pseudobulks). That makes **nine** readings marked withdrawn or superseded. Before the update, the ninth, the atlas Δρ, was withdrawn only in draft manuscript text, and the repository's ledgers listed its parent reading as standing. No file stores the count three (also found by `FINDINGS_VERIFY.md` item 36 and `PROGRESS_VERIFY.md` line 15). **Three** of the nine concern analyses reported in the manuscript (MD3 Task 1, MD2 Task 2 and the atlas Δρ), but they are not three *withdrawn-and-marked* readings. preprint_10's P89 wording is therefore not accurate. A suggested replacement is at the end of this file. preprint_11 uses a different replacement, which states nine and three.

The "Where marked" column gives `file:line`. Line numbers are 1-indexed.

## A. Readings labelled withdrawn

| # | Reading (file it lives in) | What replaced it | Where marked, and the wording |
| --- | --- | --- | --- |
| A1 | SOUTH4 TF-screen run: STOP `guide_assignment` (`FINDINGS_SOUTH4.md`; `results/south4/` deleted) | SOUTH5 (columns named explicitly) | `FINDINGS_SOUTH4.md:1` "**SUPERSEDED AND WITHDRAWN by SOUTH5 (FINDINGS_SOUTH5.md).** … `results/south4/` has been deleted as untrustworthy. Nothing below is changed." Also `results/south5/PREREG.flag:9-12` "Delete results/south4/ entirely; its run continued past its own recorded STOP and its outputs are not trustworthy. Add one line at the top of FINDINGS_SOUTH4.md marking it superseded and withdrawn". |
| A2 | MD3 Task 1 `instruments_disagree_on_claim_population`: "Ruler does not decline in the pooled PartialReprog state while MD does … No winner is declared." (`FINDINGS_MD3.md:3`, `:1466`; not edited there) | MD4 Task 1 `both_lower` (PartialReprog vs NonReprog at timepoints with ≥ 50 cells each) | `FINDINGS_MD4.md:482` "… rested on n₀=3 PartialReprog cells at day 0 and is withdrawn as unsupported, superseded by MD4 Task 1." Also `FINDINGS_MD4.md:493`, `results/md4/PREREG_TASK3.flag:1`, `results/md4/t1_reading.json:45` (key `md3_status_withdrawn`), the ledgers `FINDINGS_MD5.md:433,446-448` and `FINDINGS_PLANE.md:327,340-342` ("FINDINGS_MD3.md Task 1 (withdrawn in md4)"), and `paper_package/CONTEXT.md:9,25` ("MD3's disagreement claim is withdrawn and replaced by MD4"). |

## B. Readings labelled superseded (not withdrawn)

| # | Reading (file it lives in) | Replaced by | Where marked, and the wording |
| --- | --- | --- | --- |
| B1 | SOUTH: STOP `no_raw_counts` (`FINDINGS_SOUTH.md:5`) | SOUTH2 | `FINDINGS_SOUTH.md:1` "**SUPERSEDED BY SOUTH2 (see FINDINGS_SOUTH2.md)** … Nothing below is changed." |
| B2 | SOUTH2: ran to completion; flags `scale_unknown`, `coverage_low` (`FINDINGS_SOUTH2.md:5`) | SOUTH3 | `FINDINGS_SOUTH2.md:1` "**SUPERSEDED BY SOUTH3 (see FINDINGS_SOUTH3.md)** … Nothing below is changed." |
| B3 | SOUTH3: ran to completion on the 4,914-gene subset (`FINDINGS_SOUTH3.md:5`) | SOUTH4 | `FINDINGS_SOUTH3.md:1` "**SUPERSEDED BY SOUTH4 (see FINDINGS_SOUTH4.md)** … Nothing below is changed." |
| B4 | SOUTH5: STOP `control_ambiguous` (`FINDINGS_SOUTH5.md:5`) | SOUTH6 | `FINDINGS_SOUTH5.md:1` "**SUPERSEDED by SOUTH6 (FINDINGS_SOUTH6.md)** … nothing else in this file is changed and `results/south5/` is kept." |
| B5 | MD2 Task 2 `instruments_disagree`: "MD moves per cluster, the frozen ruler does not. The instruments disagree on the same cells. …" (`FINDINGS_MD2.md:3`, `:620`) | MD4 Task 1 `both_lower` | `FINDINGS_MD4.md:366-370` and `:380-384` "Sentences superseded (quoted, not edited in the source files):", followed by this sentence and the A2 sentence. `FINDINGS_MD4.md:515-517` adds that it "was also a within-cluster trajectory test". Not in the MD5/PLANE ledgers. |
| B6 | POSCTRL P1 "pipeline broken" (C2a bootstrap angle 41.3° > 35°) (`FINDINGS_POSCTRL.md:10`) | C1/C3 under new bars | `FINDINGS_POSCTRL.md:3` "P1 fired under a miscalibrated angle bar; superseded 2026-09-16"; `:76-82` "## Supersession 2026-09-16 … The ≤35° angle bar is **retired**." Also `PROGRESS_POSCTRL.md:3`, `results/posctrl/c1_report.txt:11`, `c3_report.txt:11`. |
| B7 (was C1) | Atlas Δρ (ruler − MD) on 93 pseudobulks: +0.202, CI +0.006 to +0.440, `excludes_zero=True` (`results/md3/t2_delta_rho.csv` row `cxg_ruler_minus_MD`). Part of MD3 Task 2 `ruler_stronger_age_in_ge2`, "Δρ CI excludes zero in favour of the ruler in 2 of 3 named cohorts" (`FINDINGS_MD3.md:1521`). | The per-person re-analysis, with the atlas regrouped into 65 donors: +0.171, CI −0.118 to +0.412, `excludes_zero=False` (`results/paper_figs/cxg_true_donor_results.csv` line 6, written 2026-09-23). **Not pre-registered** (no PREREG flag in `results/paper_figs/`; run after the 93-pseudobulk result was known). | Marked 2026-10-02 by an appended section "Supersession note, appended 2026-10-02: atlas ruler − MD Δρ on 93 pseudobulks" at the end of `FINDINGS_MD3.md`, `FINDINGS_MD4.md`, `FINDINGS_MD5.md`, `FINDINGS_PLANE.md`, `PROGRESS_MD4.md` and `PROGRESS_MD5.md`, and their byte-identical copies under `paper_package/`. Wording: "The ruler − MD Δρ computed with each of the 93 CELLxGENE atlas pseudobulks (a19d1667) as a unit … is **superseded** by the per-person re-analysis … **The per-person re-analysis was not pre-registered**". A dated paragraph was also appended at the end of `paper_package/CONTEXT.md` ("'The paired Δρ result still stands' no longer holds for the fibroblast atlas"). Nothing above those notes was changed. The PREREG flags that say Task 2 stands (`results/md4/PREREG_TASK3.flag:1`, `results/md5/PREREG_TASK3.flag:6`) were not edited. Before 2026-10-02 it was withdrawn only in the preprint_3 draft text quoted in `FINDINGS_VERIFY.md:40` and `:106`. preprint_11 also marks it in the manuscript (P94: "The ruler-minus-MD difference computed on the 93 pseudobulks is superseded by the per-person value in Table 2."). With the atlas row regrouped, only one of the three named cohorts (GSE226189) still has a ruler − MD interval excluding zero. |

## C. Withdrawn in manuscript text only, not marked in the repository

None since 2026-10-02. The one reading first listed here (C1, the atlas ruler − MD Δρ on 93 pseudobulks) is now marked superseded in the repository and appears as B7 above.

## D. Recorded as superseded or corrected, but not readings (listed so they are not counted)

- `FINDINGS_GTEX.md:124`, `:141`: the original Prompt B bars were "superseded before any GTEx Stage 1 number". This is a decision rule changed before any number existed.
- `FINDINGS_FIBRO3.md:19`: "The primary external cell is superseded from GSE113957 to CELLxGENE dataset `a19d1667-…`", written before any Task 4 ρ. GSE113957 was never scored (FPKM deposit), so no reading was replaced. This design change is why the atlas is preprint_10's external single-cell cohort (Table 1 row 4, Table 3).
- `FINDINGS_SOUTH.md:228-232`: "Corrected here: 1,836 transcription factors and 17 control rows, so `no_controls` does **not** fire." This corrects which keys fired within the same file. It "does not change the task's outcome".
- The atlas regrouping (93 pseudobulks into 65 donors) of the ruler ρ, 0.455 to 0.433. preprint_3 described it as "One later correction is reported alongside the original result rather than in place of it" (`FINDINGS_VERIFY.md:274`). The transfer still passes, so it is a correction, not a withdrawal. preprint_10 Table 1 now shows only the 65-donor row, and P89 no longer has the sentence disclosing that the regrouping came after the pre-registered result was known. (Updated 2026-10-02: preprint_11 P94 states that the grouping was not pre-registered and gives the pre-registered 93-pseudobulk ρ = +0.455. The ρ correction is still not counted as a reading.)
- `paper_package/CONTEXT.md:25`: FIBRO's "no decline", PLANE's bend-versus-lockstep reading and the SEA-AD transfer "were not carried forward". That is a choice about what to write up, not a withdrawal. The ledgers list FIBRO Stage 2 as standing.
- Conditional supersession clauses that never fired, so nothing was superseded:
  - MD2 Task 1 would have superseded the FIBRO3 per-cluster reading if a cluster declined (`FINDINGS_MD2.md:77`). It fired `negative_holds` instead (`FINDINGS_MD2.md:333`).
  - MD3 Task 1 would have superseded MD2's `negative_holds` if the ruler declined (`FINDINGS_MD3.md:56`). It fired `instruments_disagree_on_claim_population` instead.
- Document-level supersessions: `paper/PREPRINT_DRAFT.md:1` ("superseded by the docx"), preprint_3 and preprint_9 (`results/verify/FINDINGS_VERIFY_V9.md:20`).
- (Added 2026-10-26, CORR1; **paragraph 142 added 2026-10-26, CORR3**.) Paper 1 v1.2, Supplementary Note 3 (PDF page 13; docx body paragraphs 144 and 142):
  - **What it said:** (a, paragraph 144) "every one of the 1,836 perturbations moved more than twice the distance between the two centroids (smallest ratio 2.10 …) … it is the same limitation the main text discloses for the same-platform test." (b, paragraph 142) "In space C both perturbations move away from the young centroid (+106.6 and +113.8, q = 1.0), and in space B both end farther from it (+32.9 and +35.1)."
  - **What FLOORCHK found** (`FINDINGS_FLOORCHK.md`, reading `CLAIM1-PARTIAL: 18 of 1,836`): non-targeting cells with no effect read more than twice the gap at every group size in the screen. After the noise correction, 18 of 1,836 are estimated to move more than twice the gap, against about 5 expected by chance.
  - **What was wrong with (b)** (`FINDINGS_CORR2CHK.md` step 4, from stored COORDSPACE tables): neither clause separated those two perturbations from the screen. In space C no perturbation's change in distance could have been negative (`n_delta_negative` = 0 of 1,836), and in space B none is negative either while both of these two sit *inside* their own cell-count-matched floors, as 1,088 of 1,836 do. The quoted `q = 1.0` came from a permutation null SOUTH8 found anti-conservative in that space.
  - **Where it is corrected:** `CORRECTIONS_PAPER1.md`, sections (a) and (b), and the revised source v1.3 (diff in `results/corr1/paper1_diff.txt`). v1.3 carries no in-document correction note and inserts no paragraph, so the bullet is paragraph 142 in both v1.2 and v1.3.
  - **Why it is not counted:** both are descriptive statements, not gated readings. The gated SOUTH6 result, that no perturbation is `toward_young` in S1 O1→Y1, stands, and so does the decision not to treat the two space A labels as candidates — it rests on the paragraph's other three diagnostics.
  - **What was not edited:** the same raw statement in `FINDINGS_RATIO.md:35-36`, as a reproduction check (R1b), and as a premise at `FINDINGS_RATIO.md:89` and `results/ratio/PREREG.flag:33`. For (b): `SUPPLEMENTARY_NOTE_COORDSPACE.md`, which carries the bullet verbatim in a file headed "text for the preprint". Its wording is unedited; a dated one-line marker above the bullet names v1.3 paragraph 142 and points to `CORRECTIONS_PAPER1.md` passage (b), and a one-line note under the file's title states that the manuscript is authoritative where the two differ. Row 11 of that notice's Appendix A records it.

- (Added 2026-10-26, CORR3.) Paper 2 v1 (https://doi.org/10.5281/zenodo.23591623), paragraphs 31 and 30. The corrected text is the paper's v2:
  - **What it said:** (a, paragraph 31) "the four genes are also uncalled and move away from young donors in the Southard screen". (b, paragraph 30) "None of the 80 perturbations called closer to early-passage cells is called closer to young donors, and all 80 move away from them."
  - **What CORR2CHK found** (`FINDINGS_CORR2CHK.md`, readings `PAPER2-DEPTH: A1b` and `PAPER2-DEPTH: C5`): both clauses read the sign of a change in distance against zero, when the size-matched no-effect baseline at those depths is itself positive. For (a), 0 of 4 factors exceed the 95th percentile of no-effect samples at their own cell count (0.2 expected by chance) and 3 fall below the median (2.0 expected) — the clause has no support. For (b), 56 of 80 exceed their own 95th percentile (4.0 expected) and 1 falls below its median (40.0 expected) — the direction holds for the set and only the universal "all" fails, for the 24 smallest groups.
  - **Where it is corrected:** `CORRECTIONS_PAPER2.md` and the paper's v2, whose source is `paper/paper2_passage_draft_v0.8.docx` (diff in `results/corr2/paper2_diff.txt`).
  - **Why it is not counted:** both are descriptive statements, not gated readings. The paper's conclusion rests on the calls — 80 of 340 against 5 of 340, 17 expected by chance each — which are gated by size-matched thresholds and are unaffected. The other 17 inventoried statements are `PAPER2-NOT-AFFECTED`, and A1a ("uncalled") stands at power 0.925 to 1.000.
  - **How strong the pre-registration was:** weaker than FLOORCHK's. The flag (`results/corr2chk/PREREG.flag`, 2026-10-26T11:21:38+09:00) was written after an inventory that already implied both answers, because the inventory was asked for first; its Part 0 discloses this. The four conclusion wordings were fixed before the inventory ran, and both clauses are universal claims, so the finding does not turn on where a threshold sits. `CORRECTIONS_PAPER2.md` states this in full.
  - **What was not edited:** `results/corr2chk/readings.txt` and `FINDINGS_CORR2CHK.md`, which keep "1 move away less far" where the v2 source now reads "1 moves away less far"; the figure-provenance records in `results/paper2_figs_fix/`, which define themselves as the sign of `delta_gtex_young`; and `FINDINGS_AXANG.md:869`, whose directional half CORR2CHK supersedes. Also `FINDINGS_CORR2.md`, `PROGRESS_CORR2.md` and `src/corr2_run.py`, which state that paper 2 is not public — the premise CORR2 and CORR3 worked under, corrected here on 2026-10-26 when the author gave the DOI. Those files are the record of what those tasks believed and are marked superseded here rather than edited.

## Which of these bear on claims made in preprint_10

preprint_11 has the same Results, tables and figures as preprint_10, so this table holds for it too. Its paragraph numbers are the same. (Updated 2026-10-02: C1 is now B7.)

| # | Bears on preprint_10? | Where | How |
| --- | --- | --- | --- |
| A2 | **Yes** | P20–P23, P60, P65, P68, P84, P89 | Same population: the aged donor's partially reprogrammed cells in GSE297234. preprint_10 has the ruler scoring those cells younger than day-0 (+0.65 vs +3.22, P60), consistent with the superseding MD4 Task 1 and contrary to the withdrawn "ruler does not decline". |
| B5 | **Yes** | same as A2 | Superseded by the same MD4 Task 1 contrast. Its "instruments disagree" no longer stands. |
| B7 (was C1) | **Yes** | Table 2 row 2, P47, P54, P56 (preprint_11 also P94) | preprint_10 reports the regrouped interval (includes zero). The superseded 93-sample interval is the pre-registered Task 2 result for this cohort. |
| A1, B1–B4 | Only P89's statement about the archive | P89 ("a published transcription-factor screen … not reported in this paper") | The TF screen is not reported in preprint_10. |
| B6 | No | — | Brain positive control; brain analyses are not reported (P89). |

## Suggested corrected wording for P89

Replace the sentence "Where a later analysis superseded an earlier reading, the earlier text was retained and marked withdrawn rather than edited; three readings were withdrawn in this way." with:

> Where a later analysis superseded an earlier reading, the earlier text was kept unedited and the later record quotes it and states the change. Three earlier readings that bear on this paper were superseded. Two concern the reprogramming dataset: a finding that the ruler and MD disagree on partially reprogrammed cells, withdrawn because it rested on three day-0 cells, and an earlier per-cluster disagreement between the two instruments; both were replaced by a direct comparison of partially reprogrammed with non-reprogrammed cells at the same timepoints. The third is the fibroblast-atlas comparison of the ruler with MD: as first run, on 93 samples, its interval excluded zero, but that counted repeated samples from the same person as independent; after the atlas was regrouped into 65 donors, which was done once the pre-registered result was known, the interval includes zero (Table 2). Superseded runs of the transcription-factor screen are recorded in the archive.

The rest of P89 can stay. If this wording is used:

- (a) `README.md` lines 145–146 should match it.
- (b) The C1 withdrawal should also be recorded in the repository, for example as a ledger line next to `FINDINGS_PLANE.md` Task 3. Today no repository file marks it, and the ledgers say the opposite.
- (c) A shorter alternative that avoids a count: "… the earlier text was kept unedited and the superseding record quotes it; the superseded readings are listed in results/verify/WITHDRAWN_READINGS.md."

## Counts used by preprint_11 P89 (added 2026-10-02)

| | Count | Readings |
| --- | --- | --- |
| Marked withdrawn | 2 | A1, A2 |
| Marked superseded | 7 | B1, B2, B3, B4, B5, B6, B7 |
| **Marked withdrawn or superseded, across the project** | **9** ("Nine") | all of the above |
| Of these, concern analyses reported in the manuscript | **3** ("three") | A2, B5, B7 (the "Yes" rows of the table above) |
| Withdrawn in manuscript text only | 0 | section C is empty |

preprint_11 P89 now reads: "Where a later analysis superseded an earlier reading, the earlier text was retained and marked withdrawn or superseded rather than edited. Nine readings across the project are marked this way, three of which concern analyses reported here; all are listed in the repository." The suggested wording above was not used. Follow-ups (a) and (b) are done: `README.md` matches P89, and B7 is marked in the ledgers. (c) was not used, because preprint_11 states the counts.
