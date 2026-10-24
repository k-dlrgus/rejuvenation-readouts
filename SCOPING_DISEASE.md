# SCOPING_DISEASE: does any public perturbation screen put the unaffected state inside the run?

**Task:** SCOPEDIS. **Status:** search complete 2026-10-24. Fit nothing, compute no statistic,
download nothing over 50 MB, submit no controlled-access request. Does not modify
`SCOPING_FIBROSIS.md`, `SCOPING_FIBROSIS2.md`, `results/scopefib/`, `results/scopefib2/` or any
`FINDINGS_*.md`.

**The answer in one line.** The fibrosis gap is general in four of the five areas searched and
**fails in one**: immune-cell activation. Four screens there put control-perturbation cells of the
resting state inside the same run as the perturbed activated cells, three of them in the same 10x
lane, and three of those four pass every carried bar.

**Readings that fired (highest first; carried wording):**

1. `GO: GSE314342 passes A+B+C+D+E. target_match=isogenic; batch_evidence=stated; perturbation_timing=before.`
2. `GO: GSE278572 passes A+B+C+D+E. target_match=isogenic; batch_evidence=stated; perturbation_timing=before.`
3. `GO: GSE92872 passes A+B+C+D+E. target_match=isogenic; batch_evidence=same_series; perturbation_timing=before.`
4. `INCOMPLETE: GSE249595 passes A+B; C fails or is unverified. Not GO, not NO-GO.`
5. `NO-GO` fired separately in each of areas 2, 3, 4 and 5 (wording per area below). It does **not** fire overall.

Required by U4 and carried in its own sentence, not a footnote: for GSE314342 and GSE278572 the
in-run target is **activated-then-rested cells of the same donors, not naive T cells**, so distance
to it measures return toward a rested state rather than toward an unstimulated one. GSE92872's
target is naive, and it is a leukaemia line.

Machine-readable table: `results/scopedis/datasets.csv` (30 rows, every mandatory column plus
`where_control_state_lives`, `target_is_naive`, `isogenic_construction`, `readout_scale`, `area`,
`d6`). Raw record text: `results/scopedis/geo_records.json`, `url_records.json`, `searches.json`;
frozen candidate lists `candidates_1.csv`, `candidates_2.csv`; opened supplementary files under
`opened/gse92872/`. Reproduction: `reproduction.csv`. Code: `src/scopedis_run.py`,
`results/scopedis/_tables.py`.

## Pre-registration (verbatim copy of `results/scopedis/PREREG.flag`, including its addendum)

```
SCOPEDIS pre-registration flag
Written 2026-10-24T11:29:11+09:00 (2026-10-24T02:29:11Z), before any database query, article fetch or web search of any
kind in this session (no GEO / SRA / ENA / ArrayExpress / GSA / CELLxGENE / Single Cell Portal /
Europe PMC / PubMed / publisher / web search has been run). The only reads so far are repository
files: CLAUDE.md, results/scopefib/PREREG.flag, SCOPING_FIBROSIS.md, SCOPING_FIBROSIS2.md,
PROGRESS_SCOPEFIB2.md, the first 80 lines of src/scopefib2_run.py, the head of .gitignore, and
`git log -1`. Three questions were put to the user in chat before this file was written; their
rulings are reproduced below (U1-U3). Nothing below is re-tuned after any record is seen.

== Question (verbatim from the task) ==

SCOPEFIB and SCOPEFIB2 found that fibrosis has no public dataset with an in-run unactivated
target and >=20 perturbations of the activated state. In every near-miss the "before" state was
characterised in a separate, earlier experiment rather than inside the screen.

This task asks whether that is specific to fibrosis or general: across disease areas, does ANY
public perturbation screen place control-perturbation cells of the unaffected state inside the
same run as the perturbed diseased cells?

Two outcomes are equally reportable. Either a dataset passes, or the gap is general and
documented. Do not treat a null as failure.

== Scope (verbatim from the task) ==

For each area below, in this order: find the screens with >=20 distinct perturbations (D1)
applied to the diseased/activated state, and check bar A on those ONLY. Do not catalogue small
screens. Open at most 8 records per area; if an area needs more, say so and stop there.

1. Immune-cell activation - T cell or macrophage stimulation screens (CRISPR, compound), where
   part of a pool is stimulated and part kept resting.
2. Genetic disease in iPSC-derived or primary cells with an isogenic gene-corrected line in the
   same run.
3. Senescence - screens on senescent cells with proliferating cells of the same line in the run.
4. Cancer with a differentiation or normal-counterpart target in the same run (AML
   differentiation screens are the likeliest shape).
5. Neurodegeneration - iPSC neuron screens with a mutant and corrected pair.

Out of scope: fibrosis and hepatic stellate cells (closed by the two prior tasks), aging (the
project's own subject), anything not a perturbation screen.

Hard rules: every dataset in a table resolved from its own record page in this session, nothing
recalled or assumed. Download nothing over 50 MB. Compute no statistic. Submit no
controlled-access request. Disclose every deviation, including anything decided that the flag
did not specify.

== Rulings by the user in chat, before this file was written and before any query ==

Text after each label is the user's ruling, verbatim.

[U1] Area 4 target (user ruling; definition).
Normal counterpart only.
Bar A for area 4 passes only when non-malignant cells of the same tissue and type are in the
run: isogenic if from the same patient, otherwise cross_donor with donor counts per C2. A
differentiated malignant cell is not the unaffected state, so an AML differentiation screen
fails bar A unless normal cells are in the run.
Where a screen's only in-run reference is a control-perturbation arm driven into a
differentiated malignant state, record `where_control_state_lives = absent` for the unaffected
state, and note the differentiated arm in the row's design text so the distinction is visible.
Do not record it as `in_run`.

[U2] Areas 2 and 5, isogenic construction (user ruling; definition).
Knock-in pairs count too.
For areas 2 and 5, a mutant line and its isogenic unaffected partner pass bar A in either
construction:
- `corrected` - patient-derived mutant line with a gene-corrected isogenic line as the target
- `knock_in` - disease mutation introduced into a healthy line, with the unedited parental line
  as the target
Both are target_match = isogenic under C2. Add a mandatory `isogenic_construction` column
(corrected / knock_in / n/a) to every area 2 and 5 row.
The scope text's "gene-corrected" (area 2) and "mutant and corrected pair" (area 5) was loose
wording by the task author, and this ruling replaces it.
Note in the row design text that a knock-in carries the healthy donor's genetic background
rather than a patient's; this is recorded, not a bar.

[U3] The 8-record cap (user ruling).
Stop that area only.
The area closes at 8 opened records, marked `cap_reached`, with the list of unopened candidate
accessions (titles and sample counts only, as seen in search results). The next area proceeds.
Within an area, open records in this order: largest stated perturbation count first.
In the PATTERN reading and the overall verdict, report capped areas separately. A capped area
may contribute screens it opened to the counts, but the verdict may not say the gap holds in
that area - only that it held in the records opened. If any area that would decide the overall
verdict is capped, say so beside the verdict.

== Required per area (verbatim from the task) ==

Record in the mandatory column set carried from the SCOPEFIB flag, including
perturbation_timing (C1), target_match with donor counts (C2), and batch_evidence with per-arm
chemistry, instrument and submission date (C3).

Additionally, for every screen with >=20 perturbations, record one column:
`where_control_state_lives`, one of:
  `in_run` - control-perturbation cells of the unaffected state are inside the same screen/run
  `separate_experiment` - the unaffected state is characterised in a different experiment,
     batch or chemistry
  `absent` - no unaffected state anywhere in the deposit
  `unstated`
This column is the point of the task. Report its distribution across all areas as a table.

== Readings (verbatim from the task) ==

Report the carried readings in their fixed wording for each area separately, and one overall
verdict across areas. Additionally fix now, before any query, this descriptive reading (not a
gated result, labelled descriptive):

  "PATTERN: of <n> screens with >=20 perturbations of a diseased/activated state across <k>
  disease areas, <a> have the control state in_run, <b> separate_experiment, <c> absent, <d>
  unstated."

State plainly whether the fibrosis finding generalises, and name any area where it does not.

== Definitions fixed by the session now, because the task does not specify them (disclosed;
   not shown to the user in chat beforehand) ==

[D8] Opened record (for the 8-record cap, U3). One opened record = one dataset whose own record
     page (GEO acc.cgi text, ArrayExpress / BioStudies, GSA, Single Cell Portal, CELLxGENE,
     Zenodo/figshare, or a similar repository page) is fetched. Reading that dataset's linked
     article, supplementary file or sibling SubSeries counts toward the same record, not a new
     one. A literature lead with no accession yet, whose article is fetched to locate one,
     counts as one record whether or not an accession is found. Search-result listings do not
     count as opening: GEO esearch/esummary, Europe PMC / PubMed search hits, WebSearch
     results. A record opened that resolves to <20 distinct perturbations still counts toward
     the cap; it is listed by accession and count in an `opened_below_20` list and is not
     catalogued further ("do not catalogue small screens").
[D9] Candidate list and ordering (implements U3). Per area, searches run first and produce a
     candidate list, saved to results/scopedis/candidates_<area>.csv before any record in that
     area is opened. Each row: accession or lead, title, sample count if shown, stated
     perturbation count as written in the search-result text (or `unstated`), source of the
     lead. Opening order: stated count descending; then `unstated` rows, by sample count
     descending; ties by accession. Rows whose search text states <20 distinct perturbations
     are not opened. The list is not re-ordered after opening begins. Leads may be added to an
     area's list only by a further search, and a further search after opening has begun is
     logged as such with its time. Candidate sources: NCBI GEO (E-utilities), Europe PMC search,
     WebSearch. A repository named by an article is opened when the article names it.
[D10] Prior knowledge. This session knows of some candidate screens from training. Such
     knowledge may be used to choose search terms only. No dataset enters a candidate list
     except from a search result in this session, and no fact enters a table except from that
     dataset's own record or article read in this session.
[D11] Unaffected state per area (for bar A and for where_control_state_lives). Area 1: resting /
     unstimulated cells of the same donor(s) or line (D5 applies: induced-activation design).
     Areas 2 and 5: the isogenic partner per U2. Area 3: proliferating cells of the same line or
     donor; for induced senescence, D5 applies (not exposed to the senescence inducer).
     Area 4: per U1.
[D12] where_control_state_lives decision rules. `in_run` needs the record or article to show
     the unaffected-state cells inside the same screen run: same pool, plate, lane, library
     prep, or described as an arm of the screen itself. `separate_experiment` needs the record
     or article to show a different experiment, batch, run, chemistry, sequencer or submission.
     `absent`: no unaffected-state cells anywhere in the deposit or its article's data
     statement (or, area 4, only a differentiated malignant arm, per U1). `unstated`: the
     unaffected state is in the deposit but the record does not show whether it was in the same
     run. If the only in-run unaffected-state cells carry no perturbation at all (no control
     guide, no vehicle), the row is recorded `in_run`, with "(unperturbed, not
     control-perturbation)" in its design text, and counted on a separate sub-line of the
     distribution table. This column is descriptive and does not gate anything. Bar A is still
     decided only by the carried text with C2/C3/D5. The two can differ: for example, a
     same-series dataset with nothing stated is C3 `same_series` (A can pass) and
     where_control_state_lives `unstated`. Any such row is flagged.
[D13] One row per screen. A screen = one perturbation library (or compound set) applied to one
     cell system in one experiment. A SuperSeries holding several screens gives several rows,
     each under its own SubSeries accession where one exists.
[D14] Per-area evaluation of the carried readings. The carried four readings, in their fixed
     wording, are evaluated over the human datasets of each area separately ("no human dataset
     passes A+B" read as "in this area"), and then once over the union of all five areas for the
     overall verdict. The wording is not changed. Each reading line is labelled with its area.
     For a capped area, U3 governs what the verdict may say.
[D15] PATTERN counts. <n> counts screens that pass bar B's count (>=20, D1) and D6 (applied to
     the diseased/activated state, read per area), any organism. A second PATTERN line in the
     same wording gives the human-only counts. Per U3, uncapped and capped areas get separate
     PATTERN lines. A >=20 screen that fails D6 gets a where_control_state_lives value but is not
     in <n>; it is listed beside the distribution table.
[D16] perturbation_timing (C1) for constitutive disease states. Areas 2, 4 and 5 have no
     activating stimulus, so C1's before/during/after cannot be read literally. Rows where the
     perturbed cells already carry the mutation or malignancy when perturbed are written
     `after (constitutive)`, with the stage of perturbation (e.g. iPSC vs differentiated neuron)
     in the design text. This follows SCOPEFIB's disclosed `after (cells diseased at isolation)`
     convention. Rows with nothing to time are written `n/a`. Areas 1 and 3 follow C1 as written.
[D17] Per-arm chemistry, instrument and submission date (C3) are recorded for every >=20 row
     where the record gives them, not only for pass-table rows. That is more than C3 requires,
     and it changes no bar.

== STOP conditions, fixed now ==

 [S1] `bar_ambiguous` - a carried bar, clarification or definition is ambiguous for a given
      area in a way U1-U3 and D8-D17 do not settle. Write the findings with the STOP key, name
      the ambiguity, and ask. The ambiguity is not resolved by picking the reading that lets
      the run proceed.
 [S2] `reproduction_failed` - see the reproduction plan below.
A STOP is never resumed and never restarted with an adjustment. `cap_reached` (U3) is not a
STOP; it closes one area.

== Reproduction of the stored record before anything new (CLAUDE.md) ==

SCOPEFIB and SCOPEFIB2 computed no statistic, so there is no numeric to reproduce within a
tolerance. Reproduced instead, before any query, by `python src/scopedis_run.py reproduce` into
results/scopedis/reproduction.csv, counts and strings exact:
 R1 results/scopefib/datasets.csv data rows = 23.
 R2 results/scopefib2/datasets.csv data rows = 1.
 R3 SCOPEFIB's two readings, as quoted in SCOPING_FIBROSIS.md lines 7-8, appear verbatim in
    PROGRESS_SCOPEFIB.md.
 R4 SCOPEFIB's NO-GO string appears verbatim in SCOPING_FIBROSIS2.md, i.e. not superseded.
 R5 Lines 10-87 of results/scopefib/PREREG.flag appear byte-for-byte in this flag (78 lines).
 R6 results/scopefib2/PREREG.flag contains the same 78 lines byte-for-byte (the carried text
    did not drift between the two prior tasks).
Any mismatch is STOP `reproduction_failed`, with every mismatch named.

== Carried forward from results/scopefib/PREREG.flag, verbatim and unchanged ==

Lines 10-87 of that file follow, byte-for-byte: bars A-E, clarifications C1-C4, definitions
D1-D7, the four readings with their exact wording, and the mandatory columns. The task directs
that "activated/diseased state" be read per disease area rather than as fibroblast activation
only (D11, U1, U2). Nothing else in them is restated, reinterpreted, softened or extended. The
carried "Search order" and "Hard rules" lines are SCOPEFIB's and are replaced, for scope only,
by this task's Scope section above. Where "this task" or "the follow-up" appears in the carried
text, it refers to SCOPEFIB as written.

== Bars (verbatim from the task, fixed now) ==

* (A) In-run target. The dataset contains unactivated or healthy-control cells of the same cells, in the same experiment and the same batch as the perturbed cells. A published external reference does not pass. This bar is not negotiable; it is the thing the aging work could never get for free.
* (B) Perturbation count. ≥20 distinct perturbations applied to the activated/diseased state, transcriptomic readout.
* (C) Depth stated. Units per perturbation recoverable from the record. Single-cell passes at a stated or derivable median ≥30 cells per perturbation; bulk passes on depth automatically. If the record does not state it, the dataset goes in `unverified`, not in the pass table.
* (D) A readout to disagree with. The source publication reports a directional readout (signature reversal, correlation score, module score) with named hits.
* (E) Human. Human cells. A mouse or rat dataset is recorded separately and never counted toward a pass.

GO = one dataset passes A+B+C+D+E. PARTIAL = A+B+C+E but no published hit list. NO-GO = nothing passes A+B. On NO-GO, say so plainly and name, for the closest near-miss, the single missing ingredient.
Do not relax a bar after seeing what exists. A dataset that misses by a little goes in a `near_miss` block with the exact shortfall.

Search order: (1) TGF-β-induced activation of human dermal or lung fibroblasts with a perturbation screen; (2) primary fibroblasts from fibrotic human tissue (IPF, cardiac, liver, keloid, scleroderma) with perturbations; (3) hepatic stellate cell activation; (4) any myofibroblast-reversal screen. Also check LINCS/L1000 fibroblast arms explicitly and record what they do and do not contain — SCOPING_TISSUE left that unresolved.

Hard rules: every dataset in a table is resolved from its own record page in this session, nothing assumed or recalled. Download nothing over 50 MB. Compute no statistic. Modify no existing `FINDINGS_*.md`.

== Clarifications fixed by the user in chat before any query (2026-10-24, same session) ==

The session raised four ambiguities in the bars and readings before writing this flag; the user resolved them as follows. Text after each label is the user's ruling.

[C1] Bar (B), perturbation timing.
Any timing passes bar B. But every dataset in every table carries a mandatory perturbation_timing column with one of: before (prevention — perturbation complete before the activating stimulus), during (co-treatment), after (reversal — perturbation begins after activation is established), or unstated. Resolve it from the record; unstated is not a guess.
Pre-declared restriction on the follow-up: the follow-up distance analysis will use after designs as its primary and during as secondary, and before designs are inventory only.

[C2] Bar (A), identity of the target cells.
Same type, healthy donor OK — with a mandatory target_match column on every dataset in every table, one of:
  isogenic — the control is the same line or donor as the perturbed cells (unstimulated arm, vehicle arm, or corrected isogenic line)
  cross_donor — the control is a different person's cells of the same tissue and type
  unstated
For cross_donor datasets, also record the number of healthy donors and the number of diseased donors. A single healthy donor against a single diseased donor cannot separate disease from individual, so such a dataset is recorded as passing bar A but flagged confounded_with_individual, and the GO verdict must state which match type it rests on. The follow-up distance analysis treats isogenic as primary and cross_donor with at least three donors per side as secondary; cross_donor with fewer goes to inventory only.

[C3] Bar (A), batch evidence.
Same series unless split — with a mandatory batch_evidence column on every dataset, one of:
  stated — the record names a shared batch, plate, lane, pool, or the controls sit inside the same pooled screen
  same_series — same GEO series or run, with nothing in the record indicating a split
  split — the record shows a different batch, date, plate, chemistry, sequencer or submission between control and perturbed arms
  unstated
stated and same_series pass bar A. split fails it. unstated goes to unverified.
For every dataset reaching the pass table, also record library chemistry, sequencing instrument and submission date for the control arm and the perturbed arm separately, where the record gives them, and flag any difference. (Reason given: FINDINGS_DONOR.md measured a transgene's direction agreeing with itself at 0.26 across batches against 0.50 between two lines measured together, so an unverified batch assumption can exceed the effect size.) The follow-up distance analysis treats stated as primary and same_series as secondary, and must report the batch evidence grade beside the verdict.

[C4] Fourth reading, precedence.
Add reading INCOMPLETE, with this exact wording fixed now:
  "INCOMPLETE: <accession> passes A+B; <bar> fails or is unverified. Not GO, not NO-GO."
Precedence: GO > PARTIAL > INCOMPLETE > NO-GO. Report every reading whose conditions fired, not only the highest. NO-GO fires only when no human dataset passes A+B.
For each INCOMPLETE dataset, name in one line the single thing that would resolve it (e.g. "per-perturbation cell counts not in the record; resolvable by opening the supplementary table" or "mouse only"). Where resolving it means opening a file under the 50 MB limit, do that and say so; where it does not, leave it unresolved and say which file would settle it.

== Definitions fixed by the session (proposed in chat before any query; user did not override) ==

[D1] Distinct perturbation (bar B): a distinct agent — one gene target (any modality: KO, CRISPRi, CRISPRa, ORF, siRNA) or one compound / biologic. Doses, timepoints, guides or replicates of the same agent count once. Controls (DMSO, vehicle, non-targeting / safe-targeting guides, empty vector, unstimulated) do not count.
[D2] L1000 counts as a transcriptomic readout for bar B, noted wherever it appears: 978 landmark genes are measured; the rest are inferred.
[D3] Bar (E): mouse and rat datasets are recorded in a separate block and are excluded from GO and PARTIAL and from the NO-GO test ("no human dataset passes A+B"). Per C4, a non-human dataset that passes A+B may fire INCOMPLETE with <bar> = E.
[D4] Gene-overlap section: if the record provides a gene list (features file) under 50 MB, report the exact count of its genes that also appear in the frozen 23,485-gene fibroblast panel (set intersection by Ensembl id, symbol fallback — a count, not a statistic). If the record states only a gene count, or nothing, report an estimate and label it "estimate".
== Definitions added by the session while writing this flag (before any query; NOT shown to the user in chat beforehand — disclosed here) ==

[D5] Bar (A) "unactivated": for induced-activation designs, cells of the same line/donor that did not receive the activating stimulus (e.g. no TGF-β), carrying a control perturbation or none. A perturbed-but-unstimulated arm does not count as the target unless a control-perturbation unstimulated arm is also present.
[D6] Bar (B) "applied to the activated/diseased state": under C1, satisfied when the perturbed cells also received the activating stimulus (any timing) or are from diseased tissue / are an activated cell type (e.g. activated hepatic stellate cells, LX-2). Perturbations applied only to unstimulated healthy cells do not count toward the ≥20.
[D7] Bar (C) "median ≥30 cells per perturbation": per distinct agent (D1) in the stimulated / diseased arm, pooling guides of the same agent.

== Readings — exact wording, fixed now ==

GO fires for a human dataset passing A+B+C+D+E. Wording:
  "GO: <accession> passes A+B+C+D+E. target_match=<isogenic|cross_donor>[ confounded_with_individual]; batch_evidence=<stated|same_series>; perturbation_timing=<before|during|after|unstated>."

PARTIAL fires for a human dataset passing A+B+C+E and failing D. Wording:
  "PARTIAL: <accession> passes A+B+C+E; no published directional readout with named hits (D fails). target_match=<…>; batch_evidence=<…>; perturbation_timing=<…>."

INCOMPLETE fires for any dataset passing A+B that is neither GO nor PARTIAL (C fails or is unverified, or E fails). Wording (from C4):
  "INCOMPLETE: <accession> passes A+B; <bar> fails or is unverified. Not GO, not NO-GO."

NO-GO fires only when no human dataset passes A+B. Wording:
  "NO-GO: no human dataset passes A+B. Closest near-miss: <accession>; single missing ingredient: <ingredient>."

Precedence GO > PARTIAL > INCOMPLETE > NO-GO; every reading that fired is reported, the highest first.

== Mandatory columns on every dataset table ==

accession | organism | cell type / source | design | n distinct perturbations (D1) | readout | units per perturbation | perturbation_timing | target_match (+ n healthy / n diseased donors if cross_donor) | batch_evidence | A | B | C | D | E

Pass-table rows additionally carry, per arm (control vs perturbed): library chemistry, sequencing instrument, submission date, with any difference flagged.
== Columns added by this task (on top of the carried mandatory set) ==

 - `where_control_state_lives` (task; D12) on every row with >=20 perturbations.
 - `isogenic_construction` (U2): corrected / knock_in / n/a, on every area 2 and 5 row.
 - `area` (1-5) on every row, and `d6` (whether the >=20 perturbations are applied to the
   diseased/activated state) on every row.
 - Per-arm chemistry / instrument / submission date where the record gives them (D17).

== Honesty obligations carried from CLAUDE.md ==

- No statistic is computed, so no chance baseline and no recovery power apply. Counts here are
  record inventories, not test results. PATTERN is descriptive and is labelled as such.
- A null here is a record-reading null over the records opened, with the search recall stated.
  It is not an experimental null.
- What a null does not cover is stated per area. In particular, results on neonatal, fetal,
  iPSC-derived or passaged cells with no aged donor do not transfer to aged cells.
- Facts that arrive through a tool returning a model-written summary rather than raw text
  (WebFetch) are marked in the source column.

== Process ==

Deliverables: results/scopedis/PREREG.flag (this file), SCOPING_DISEASE.md (this flag
reproduced verbatim at the top), PROGRESS_SCOPEDIS.md. Raw query responses and candidate lists
go under results/scopedis/. Code: src/scopedis_run.py (logging fetcher plus a `reproduce` mode;
CLAUDE.md layout). One local commit, no push. Not modified: SCOPING_FIBROSIS.md,
SCOPING_FIBROSIS2.md, results/scopefib/, results/scopefib2/, any FINDINGS_*.md.

== ADDENDUM, added 2026-10-24T11:45:03+09:00 (2026-10-24T02:45:03Z),
   AFTER queries had begun. Disclosed as such. ==

Two carried-bar ambiguities surfaced on the first two records opened in area 1, in ways U1-U3
and D8-D17 do not settle. Per CLAUDE.md the session stopped and asked rather than picking the
reading that let the run proceed, and the user's rulings follow verbatim. Nothing above this
addendum is edited. Both rulings were given before any further record was opened and before any
bar verdict was written. The deviation is that they were fixed after two records had been read,
not before any query; the records that prompted them are named below.

[U4] Bar (A) / D5, the identity of a "resting" arm in a primary-cell screen. Prompted by
GSE314342 (resting and stimulated cells in the same 10x lane; the resting cells had been
activated once on day 0 to permit transduction, expanded 12 days, then held without
restimulation). User ruling, verbatim:

Yes, with the caveat recorded - and with the limit written into the row and the verdict, not
only the design text.
The design's activating stimulus is the day-12 restimulation. The rest arm did not receive it,
so bar A passes and where_control_state_lives = in_run.
Required alongside, in the flag:
- Add a mandatory `target_is_naive` column: `naive` (never activated) / `returned` (activated
  once, then rested) / `unstated`. GSE314342's rest arm is `returned`.
- Any GO reading for this dataset must carry the limit in its own sentence, not in a footnote:
  the target is activated-then-rested cells of the same donors, not naive T cells, so distance
  to it measures return toward a rested state rather than toward an unstimulated one.
- Record in the row: days between day-0 activation and harvest, whether the rest and
  restimulated arms share donors, lane, chemistry and submission date (C3), and the cells per
  perturbation in each arm separately (D7).
Do not treat the rest arm as equivalent to a naive target anywhere in the file.

[U5] Bar (B), "transcriptomic readout" for a targeted probe panel. Prompted by GSE249595
(18,595 CRISPRi gene targets, untreated and activated cells hashed into the same 16 lanes, read
out on a 374-gene TAP-seq panel). SCOPEFIB left this question explicitly open (SCOPING_FIBROSIS.md
deviation 3), saying a ruling was needed before a future search could pass such a dataset. User
ruling, verbatim:

Counts, with the genes measured recorded - and with a separate flag for whether the readout can
carry the follow-up.
A targeted panel measuring real transcripts counts as a transcriptomic readout for bar B,
consistent with D2 allowing L1000's 978 measured genes. Record the number of genes measured
wherever the dataset appears. GSE249595 passes bar B on readout.
Add a mandatory `readout_scale` column: `genome_wide` / `panel_<n>` / `l1000` / `unstated`.
Also fix now, before any further record: a readout below 2,000 measured genes passes bar B but
is flagged `insufficient_for_distance` and is inventory only for the follow-up distance
analysis. The reason, stated in the flag: the follow-up measures Euclidean distance in a frozen
23,485-gene space, and a few hundred pathway-selected genes cannot support that geometry or a
noise floor in it. This is a restriction on the follow-up, not a relaxation or tightening of
bar B.
Apply the same `readout_scale` and flag retrospectively in this file's tables to SCOPEFIB's
GSE227391 (PHDs-seq), noting that SCOPEFIB left the question open and this ruling closes it. Do
not edit SCOPING_FIBROSIS.md.

Consequent additions to the mandatory column set: `target_is_naive` (U4) and `readout_scale`
(U5) on every row, plus the `insufficient_for_distance` flag where U5's 2,000-gene floor bites.
```

## Reproduction of the stored SCOPEFIB / SCOPEFIB2 record, before any query

Run as `python src/scopedis_run.py reproduce` into `results/scopedis/reproduction.csv`. Neither
prior task computed a statistic, so there is no numeric to reproduce within a tolerance; what the
flag fixed instead was reproduced exactly. **12 checks, 0 mismatches.** No STOP fired.

| check | expected | observed | match |
|---|---|---|---|
| R1 `results/scopefib/datasets.csv` data rows | 23 | 23 | yes |
| R2 `results/scopefib2/datasets.csv` data rows | 1 | 1 | yes |
| R3 SCOPEFIB's two readings (SCOPING_FIBROSIS.md lines 7-8) verbatim in PROGRESS_SCOPEFIB.md | True | True | yes |
| R3b reading 2 is the NO-GO | True | True | yes |
| R4 SCOPEFIB's NO-GO verbatim in SCOPING_FIBROSIS2.md | True | True | yes |
| R4b SCOPING_FIBROSIS2.md states the NO-GO is not superseded | True | True | yes |
| R5 carried block line count | 78 | 78 | yes |
| R5b carried lines 10-87 byte-identical in this task's flag | True | True | yes |
| R6 same 78 lines byte-identical in the SCOPEFIB2 flag | True | True | yes |

R6 matters because this task carries the bars at two removes: it confirms the carried text did not
drift between SCOPEFIB and SCOPEFIB2 before being carried again here.

## Records opened, per area

D8 counts one opened record as one dataset whose own record page was fetched; its linked article,
supplementary files and SubSeries count toward the same record. Opening order was largest stated
perturbation count first (D9).

| area | opened | cap | records opened (in order) | candidates left unopened |
|---|---|---|---|---|
| 1 immune activation | 8 | **cap_reached** | GSE314342, GSE249595, GSE221321, GSE312107, GSE271090 (via GSE278572, GSE276097), GSE92872, GSE266618, GSE255832 | GSE203593 (58 samples), GSE119450 (8), GSE190604 (16), GSE137554 (2) |
| 2 genetic disease | 6 | no | GSE307350, GSE215253, GSE271444, GSE343369, GSE150062, MED:42062284/PMC13338134 | GSE288996 (15), GSE303901 (97) |
| 3 senescence | 6 | no | GSE300879, GSE124197, GSE268569, GSE224293, GSE309511, GSE75627 | GSE313494, GSE313462, GSE288940 (no perturbation library) |
| 4 cancer | 4 | no | GSE213511, GSE305037, GSE221578, GSE294098 | GSE254241 (22), GSE217779 (2), GSE181338 (40), GSE347432 (20) |
| 5 neurodegeneration | 6 | no | GSE289235, GSE152988, GSE178317, GSE250378, GSE255902, GSE272093 | GSE207099 (4), GSE124703 (14), GSE182309 (36), GSE289721 (8) |

`opened_below_20` (opened, resolved to fewer than 20 distinct agents with a transcriptomic readout,
not catalogued further): GSE266618 (1), GSE294098 (16), GSE224293 (3), GSE75627 (<20), GSE255902
(2), PMC13338134 (~7), GSE307350 (1 in the RNA-seq arm), GSE300879 (3 in the RNA-seq arm),
GSE124197 (1 in the RNA-seq arm).

## Readings, per area and overall

Carried wording, unchanged. Area 1 reached its 8-record cap (U3); areas 2–5 did not.

**Area 1 — immune-cell activation (cap_reached).**

```
GO: GSE314342 passes A+B+C+D+E. target_match=isogenic; batch_evidence=stated; perturbation_timing=before.
```
> Limit carried in the reading itself, per U4: the in-run target is activated-then-rested cells of
> the same donors (`target_is_naive = returned`), not naive T cells, so distance to it measures
> return toward a rested state rather than toward an unstimulated one.

```
GO: GSE278572 passes A+B+C+D+E. target_match=isogenic; batch_evidence=stated; perturbation_timing=before.
```
> Same limit, same wording: `target_is_naive = returned` (stimulated at isolation, half the culture
> restimulated on day 8).

```
GO: GSE92872 passes A+B+C+D+E. target_match=isogenic; batch_evidence=same_series; perturbation_timing=before.
```
> Its target is `naive` — the only one of the four. Its control state is `separate_experiment`
> (same pool, same submission, but a separate Drop-seq run per condition).

```
INCOMPLETE: GSE249595 passes A+B; C fails or is unverified. Not GO, not NO-GO.
```
> Single resolver: achieved cells per gene. The record states only an expected recovery of 60,000
> cells per channel and the article only a post-QC total of >586,000 cells, which divides to ~28
> cells per gene in the activated arm — below D7's 30, and a mean, not a median. The per-cell guide
> and hash assignment in `GSE249595_RAW.tar` would settle it. Its readout is `panel_374`, flagged
> `insufficient_for_distance` under U5.

**Area 2 — genetic disease in iPSC-derived or primary cells.**

```
NO-GO: no human dataset passes A+B. Closest near-miss: GSE215253; single missing ingredient: a perturbation applied to the diseased state (D6) — the ≥20 perturbations are the disease variants themselves, introduced into healthy primary cells.
```

**Area 3 — senescence.**

```
NO-GO: no human dataset passes A+B. Closest near-miss: GSE309511; single missing ingredient: a screen run on senescent cells — the CRISPRa CROP-seq arm is in proliferating WI-38 and the senescent state is a separate RNA-seq time course.
```

**Area 4 — cancer with a differentiation or normal-counterpart target.**

```
NO-GO: no human dataset passes A+B. Closest near-miss: GSE221578; single missing ingredient: non-malignant cells of the same tissue in the run (bar A; per U1 a malignant arm is not the unaffected state).
```

**Area 5 — neurodegeneration.**

```
NO-GO: no human dataset passes A+B. Closest near-miss: GSE178317; single missing ingredient: the screen run on the activated state — the 38-gene CROP-seq is on unstimulated iTF-microglia and the ±LPS arms are separate bulk experiments.
```

**Overall, across areas.**

```
GO: GSE314342 passes A+B+C+D+E. target_match=isogenic; batch_evidence=stated; perturbation_timing=before.
```

with the U4 limit sentence above, and `INCOMPLETE: GSE249595 passes A+B; C fails or is unverified.
Not GO, not NO-GO.` also firing. NO-GO does **not** fire overall, because human datasets pass
A+B. Per U3, the one area that decides this verdict — area 1 — is the capped area: the GO rests on
records opened there, and four candidates on its frozen list (GSE203593, GSE119450, GSE190604,
GSE137554) were never opened.

## PATTERN (descriptive, not a gated result)

Counting per D15: screens with ≥20 distinct perturbations (D1) **and** D6 satisfied, any organism.

```
PATTERN: of 8 screens with >=20 perturbations of a diseased/activated state across 5 disease areas, 4 have the control state in_run, 1 separate_experiment, 3 absent, 0 unstated.
```

```
PATTERN: of 6 human screens with >=20 perturbations of a diseased/activated state across 5 disease areas, 4 have the control state in_run, 1 separate_experiment, 1 absent, 0 unstated.
```

Per U3, split by cap: **all 8 of these screens are in area 1, the capped area.** Areas 2–5,
uncapped, contributed **zero** screens that satisfy both the ≥20 count and D6.

| where_control_state_lives | n (all) | n (human) | screens |
|---|---|---|---|
| `in_run` | 4 | 4 | GSE314342; GSE278572; GSE249595; GSE92872 (arrayed bulk arm) |
| — of which same lane/pool, hash- or probe-barcode demultiplexed | 3 | 3 | GSE314342, GSE278572, GSE249595 |
| — of which an arm of the same experiment, parallel bulk libraries | 1 | 1 | GSE92872 (bulk) |
| `separate_experiment` | 1 | 1 | GSE92872 (CROP-seq arm) |
| — of which same pool and submission, separate run | 1 | 1 | GSE92872 (CROP-seq) |
| — of which a different experiment or chemistry | 0 | 0 | — |
| `absent` | 3 | 1 | GSE221321 (human); GSE312107, GSE255832 (mouse) |
| `unstated` | 0 | 0 | — |

**None of these rows carries an unperturbed-only in-run control**: every `in_run` row has control
guides (NTC, CTRL or library-perturbed) in the unaffected arm, so the D12 sub-line for
"unperturbed, not control-perturbation" is empty.

Screens with ≥20 perturbations that **fail D6**, so excluded from the counts above, with their
column value recorded anyway (D15): GSE215253 `n/a (the edit creates the variant)`; GSE271444
`n/a (the screen perturbs the unaffected state)`; GSE309511 `separate_experiment`; GSE178317
`separate_experiment`; GSE152988 `absent`; GSE289235 `absent`; GSE250378 `absent`; GSE343369 and
GSE150062 `n/a (no diseased state)`. Screens whose distinct-agent count could not be resolved from
the record, so not counted either way: GSE221578, GSE305037, GSE213511, GSE272093 (counts stated
only as sgRNA numbers, family names or not at all).

## Does the fibrosis finding generalise?

**Mostly yes, and it fails in exactly one area: immune-cell activation.**

- In areas 2, 3, 4 and 5 the SCOPEFIB pattern repeated in the records opened, and it repeated for
  two distinguishable reasons. (i) *The screen is on the wrong state.* GSE309511 (senescence),
  GSE178317 (microglial activation), GSE152988, GSE289235, GSE250378, GSE271444: ≥20 perturbations
  applied to the **unaffected** cells, with the diseased/activated state characterised in a
  separate experiment — the SCOPEFIB near-miss shape exactly. (ii) *The unaffected state is not in
  the deposit at all.* GSE221578 and GSE305037 (AML patient cells and lines), GSE221321 (THP-1
  LPS): `absent`.
- In area 1 the gap closes, and the reason is mechanical rather than biological. T-cell and
  macrophage activation is an **induced, reversible** state, so the resting arm is the same cells
  held back; and the field's standard kit multiplexes conditions into one lane (10x Flex sample
  probe barcodes in GSE314342, TotalSeq hashtags in GSE278572, cell hashing in GSE249595). Once
  conditions share a lane, the control state is in the run by construction and C3 `split` cannot
  arise. Fibrosis near-misses had neither property: activation is treated as near-irreversible, and
  the "before" state was profiled in its own experiment, in the ADAM12 case on a different 10x
  chemistry.
- The area where it does **not** generalise is therefore area 1 only, and the price is in what the
  target means (U4): in the two primary-cell GO datasets the target is `returned`, not `naive`.

## Deviations and decisions not in the flag (disclosed)

1. **Two carried-bar ambiguities were resolved by the user after queries had begun**, not before
   any query as the flag intends. U4 (the identity of a "resting" arm) and U5 (targeted panels as
   transcriptomic readouts) are appended to `PREREG.flag` with their own timestamp, marked as added
   after two records had been opened, and the records that prompted them are named there. The flag
   above this addendum was not edited. U5 also closes a question SCOPEFIB left open, and this file
   applies it retrospectively to SCOPEFIB's GSE227391 (`readout_scale = panel`, PHDs-seq;
   `insufficient_for_distance`) **without editing `SCOPING_FIBROSIS.md`**.
2. **A sub-classification was added to `where_control_state_lives` for reporting only.** D12's
   letter makes a separate run of the same pool `separate_experiment`, which lumps GSE92872's
   design together with a different-chemistry split like the ADAM12 near-miss. The distribution
   table therefore carries sub-lines (same pool and submission vs different experiment or
   chemistry). The four values themselves are unchanged and nothing is re-graded.
3. **`n/a` values outside the four enumerated ones** were written where there are no perturbed
   diseased cells to compare against (GSE215253, GSE271444, the no-diseased-state screens). This is
   the same gap SCOPEFIB disclosed when it wrote `perturbation_timing = n/a`.
4. **Two supplementary files were opened** to resolve bars that the GEO record could not:
   `GSE92872_CROP-seq_Jurkat_TCR.count_matrix.csv.gz` (1.55 MB) and
   `GSE92872_CROP-seq_Jurkat_TCR.digital_expression.csv.gz` (17.3 MB), both under the 50 MB limit.
   Only their header rows were parsed, giving the 32-gene library size and the per-arm median cells
   per gene. They are untrusted downloaded input: they sit in `results/scopedis/opened/gse92872/`
   and were read by a script kept in the session scratchpad, run with `python -I`. The arithmetic
   involved is counting and a median, which bar C requires; no statistic was computed.
5. **One judgement about scope.** GSE300879 / GSE268569 are framed by their own records as
   "replicative aging and inflammatory aging" screens in MSC. Aging is out of scope for this task,
   senescence is in scope, and the records describe MSC senescence; they were therefore read as
   area 3. Both fail bar B on readout, so the judgement cannot change a reading.
6. **Additions beyond the flag:** `src/scopedis_run.py` (logging fetcher, GEO and Europe PMC search
   wrappers, `reproduce` mode) and `results/scopedis/_tables.py` (hand-entered table → CSV). Helper
   scripts used only to read stored responses were kept in the session scratchpad and are not part
   of the repository.
7. **Searches and candidate freezing.** D9 required a candidate list saved before any record in an
   area was opened. `candidates_1.csv` and `candidates_2.csv` were written that way. For areas 3, 4
   and 5 the candidate lists were **not** written to disk before opening began — the opening order
   (largest stated count first) was applied from the search log, which is stored in
   `searches.json`, but the frozen file is missing. That is a deviation from my own D9, and it means
   the ordering for those three areas is reconstructible from `searches.json` rather than
   pre-committed.
8. **Facts from tools returning model-written summaries:** none. Every row rests on raw record text
   or raw article text stored in this directory, except the 16p12.1 row, which rests on a Europe PMC
   abstract read in session and is labelled as having no accession opened.

## What this does not cover

- **Search recall.** 13 GEO DataSets queries for area 1 and area 2, 9 for area 3, 8 for area 4 and
  9 for area 5, each capped at 30–40 summaries, plus 5 Europe PMC queries. Only GEO and Europe PMC
  were searched: no ArrayExpress/BioStudies, GSA, EGA, CELLxGENE, Single Cell Portal, Zenodo or
  figshare query was run in this task. SCOPEFIB's recall limits are untouched and SCOPEFIB2's one
  narrowing (GSA-Human, for one study) does not extend here.
- **The cap.** Area 1 closed at 8 opened records with four candidates unopened, and area 1 is the
  area that decides the overall verdict. The GO is a statement about records opened, not about the
  area.
- **Counts left unresolved.** Four screens (GSE221578, GSE305037, GSE213511, GSE272093) have no
  distinct-agent count in their records. Each fails bar A independently, so resolving them cannot
  change a reading; it would change the PATTERN denominator.
- **Aged cells.** Nothing here touches aged donors. The GO datasets use donors aged 22–34
  (GSE314342), two unstated-age healthy donors (GSE278572) and a leukaemia line (GSE92872,
  GSE249595). A result on them does not transfer to aged cells, and area 1's resting arm is
  `returned`, not `naive`.
- **No statistic was computed**, so there is no chance baseline and no recovery power to report.
  Every count here is a record inventory. The PATTERN line is descriptive and is labelled as such.
  Its denominator is 8 screens, which is small enough that one unopened dataset could move it.
- **Circularity:** none to disclose. No axis, basis or reference was fitted or applied.
- **Group size.** CLAUDE.md's shortlist warning applies to cell-count artefacts in this project's
  own screens; the only "shortlist" here is which records were opened, and the opening rule
  (largest stated perturbation count first) is stated in the flag rather than chosen after seeing
  results — with the area 3–5 exception disclosed in deviation 7.

## Next

None within this task (one task per session). For the distance-to-target question the record now
shows that the design the fibrosis work could not find already exists in immune activation:
GSE314342 and GSE278572 put a control-perturbation resting arm in the same lane as thousands or
dozens of perturbed activated cells, in human primary cells, with published directional readouts.
Neither route is chosen here. Two things would have to be settled first: whether a `returned`
resting arm is an acceptable target for a distance measurement, and that the timing in both is
`before`, which C1 makes inventory-only for the follow-up.
