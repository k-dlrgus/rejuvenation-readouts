# SCOPING_FIBROSIS2: the one `unverified` human dataset SCOPEFIB left

**Task:** SCOPEFIB2. **Status:** resolved 2026-10-24. Fit nothing, compute no statistic beyond
the cell counts bar C requires, download nothing over 50 MB. Modifies no `FINDINGS_*.md`, and
does not touch `SCOPING_FIBROSIS.md`, `PROGRESS_SCOPEFIB.md` or `results/scopefib/`.

**Scope:** the ADAM12 study's human fibroblast Perturb-seq and sibling human accessions of that
study. Nothing else was searched.

## The answer in one block

The human deposit exists and was located: **HRA014194**, in the GSA-Human archive (Genome
Sequence Archive for Human), under BioProject **PRJCA048014**. It is **controlled access**.

Bars, read from that record and from the article's Open Access full text:

| bar | verdict | why |
|---|---|---|
| **A** in-run target | **fail** | Every screened cell received TGF-β. The STAR Methods give one screen condition: transduce at MOI 0.5, add TGF-β 10 ng/mL two days later, culture three more days, harvest. There is no unstimulated control-perturbation arm inside the screen, which is what D5 requires. Taking instead the unperturbed ±TGF-β arms as the target, those are different experiments on a different library chemistry (10x 5′ kit vs the screens' 10x 3′ HT v3.1) and different samples, which is C3 `split`, and `split` fails bar A. |
| **B** ≥20 perturbations | **pass** | 590 CRISPRi and 265 CRISPRa gene targets by design (667 distinct genes, 188 in both libraries); 432 and 210 recovered after the screen's own cell filter (532 distinct). All perturbed cells received the activating stimulus, so D6 is satisfied. |
| **C** depth stated | **unverified** | The only depth statement is the retention rule "targets detected in more than 5 cells". No median is given. Mean cells per recovered target is 80.8 (CRISPRi) and 187.5 (CRISPRa) — an **estimate**, not the median D7 asks for. |
| **D** readout to disagree with | **pass** | Directional readout: log2FoldChange of each perturbation against the non-targeting controls, averaged over gene-module genes; plus a TGF-β response module score. GSVA (v1.38.2) is used for descriptive module scoring. Named hits include ADAM12, FAP, NTM, HAS2, TCF4, ATF3, ANXA1, NOX4. |
| **E** human | **pass** | Human. 8 PDF donors, ages 35–79; the screen itself ran on cells from two donors, aged 35 and 54. |

**Reading. No carried reading changes.** Bar A fails, so the dataset is neither GO, PARTIAL nor
INCOMPLETE — all three require A. The line fixed in this task's flag for that case:

> UNRESOLVED (SCOPEFIB2 scope statement, not a reading): the ADAM12 human Perturb-seq stays in
> `unverified`; <bars> could not be read because <reason>. SCOPEFIB's NO-GO stands unchanged.

does **not** apply either, because the bars *were* read. What holds instead is the carried
NO-GO, unchanged and now resting on one more resolved dataset rather than on an unread one:

```
NO-GO: no human dataset passes A+B. Closest near-miss: GSE305412; single missing ingredient: perturbation count (5 distinct knockdown agents applied across TGF-β activation; bar B needs ≥20).
```

**SCOPEFIB's NO-GO is not superseded.** The ADAM12 dataset moves out of `unverified` and into
`near_miss`, where it is now the strongest near-miss on the human side: it is the only human
dataset found in either task with ≥20 perturbations applied to an activated state (bar B) on
adult patient-derived cells. It fails on the bar SCOPEFIB called non-negotiable.

Exact shortfall, one line: **the screen has no unstimulated arm; the unactivated target was
profiled as a separate experiment on a different 10x chemistry.**

**Access.** The reading above used only the public record and the Open Access article. The
FASTQ data behind bar C are controlled access and were not requested; see §Access.

Machine-readable: `results/scopefib2/datasets.csv` (1 row, all mandatory columns),
`hra014194_experiments.csv` (39 rows, the arm-level record), `library_counts.csv` (13 counts),
`reproduction.csv`, `browser_pane_urls.csv`. Raw record text: `url_records.json` (54 requests
with URL and HTTP status), `geo_records.json` (6 GEO series), `article_excerpts.json`,
`opened/HRA014194.xlsx`. Code: `src/scopefib2_run.py`, `results/scopefib2/_tables.py`.

## Pre-registration (verbatim copy of `results/scopefib2/PREREG.flag`)

```
SCOPEFIB2 pre-registration flag
Written 2026-10-24T10:47:10+09:00 (2026-10-24T01:47:10Z), before any database query, article
fetch or web search of any kind in this session (no GSA / GEO / Europe PMC / PMC / publisher /
web search has been run; the only reads so far are repository files: CLAUDE.md,
results/scopefib/PREREG.flag, SCOPING_FIBROSIS.md, PROGRESS_SCOPEFIB.md, src/scopefib_run.py).
Nothing below is re-tuned after any record is seen.

== What this task is ==

SCOPEFIB (commit f93b244) returned, in precedence order:
  INCOMPLETE: GSE261783 passes A+B; E fails or is unverified. Not GO, not NO-GO.
  NO-GO: no human dataset passes A+B. Closest near-miss: GSE305412; single missing ingredient:
  perturbation count (5 distinct knockdown agents applied across TGF-beta activation; bar B
  needs >=20).

It left one human dataset in `unverified`: the ADAM12 study's patient-derived fibroblast
CRISPRi/CRISPRa Perturb-seq (Li et al., Cancer Cell 2026, doi 10.1016/j.ccell.2025.12.018,
PMID 41544628). SCOPEFIB could not read bars A-C: GEO SuperSeries GSE278990 holds mouse data
only, and the Cancer Cell article page returned HTTP 403. The article's data-availability
statement names two repositories, NCBI GEO and the GSA (Genome Sequence Archive, China
National Center for Bioinformation, ngdc.cncb.ac.cn/gsa). SCOPEFIB searched only GEO.

SCOPEFIB2 resolves that one entry and nothing else.

== Scope, fixed now ==

In scope: the ADAM12 study's human fibroblast Perturb-seq, and any sibling human accession
deposited by the same study.

Out of scope and not to be searched: new disease areas, new tissues, new screens, any dataset
not deposited by this study. SCOPEFIB's search recall limits (SCOPING_FIBROSIS.md, "What this
NO-GO does not cover") stand as written; this task does not extend them. A promising unrelated
dataset seen incidentally is recorded as an accession in a `leads` block and not opened.

Resolution order, fixed now:
 1. The GSA search interface, for the study, the project, or the authors.
 2. The article's data-availability and key-resources sections, via Europe PMC full text, PMC,
    or a publisher mirror (SCOPEFIB found cell.com returned content to a search index even
    though a direct fetch returned 403).
 3. The paper's supplementary key-resources table, if under the 50 MB limit.
Every URL tried and its HTTP status goes in results/scopefib2/url_records.json.

Bar A is read from the record and the methods, to this precise question: do the perturbed cells
and an unperturbed-or-control-perturbation arm that did NOT receive the activating stimulus
both appear in the same experiment and batch? D5 governs. The record is also to state whether
the screen was run on cells already activated by TGF-beta, on cells activated after
perturbation, or on unstimulated fibroblasts with activation characterised separately. A screen
run entirely on unstimulated cells fails D6 even at 54 targets.

Bar B: distinct agents per D1, counted across the CRISPRi and CRISPRa arms separately and
together. Bar C: median cells per agent per D7, computed only if a per-cell guide assignment or
a per-guide table under 50 MB is available; otherwise `unverified`. Bar D: the named hits and
the directional readout used (the article uses GSVA on gene modules; that is recorded as the
readout to disagree with).

If the deposit is controlled-access or request-gated: say so plainly, state what the request
requires, and stop. No access request is submitted in this session.

Gene-overlap section per D4, only if a reading passes.

== Reproduction of the stored record before anything new (CLAUDE.md) ==

SCOPEFIB produced no statistic, so there is no numeric to reproduce within tolerance. What is
reproduced instead, before any query, and reported as a table: the row count of
results/scopefib/datasets.csv against the 23 rows SCOPING_FIBROSIS.md claims; the two readings
recorded in PROGRESS_SCOPEFIB.md against the two in SCOPING_FIBROSIS.md, character-for-
character; the presence of the ADAM12 entry in the `unverified` block; and that lines 10-87 of
results/scopefib/PREREG.flag, spliced verbatim into this flag below, match byte-for-byte.
Counts exact, strings exact. Any mismatch is a STOP, keyed `reproduction_failed`, and every
mismatch is named.

== STOP conditions, fixed now ==

 [S1] `bar_ambiguous` - a carried bar, clarification or definition turns out to be ambiguous
      for this dataset. Write the findings with the STOP key and ask; do not resolve it by
      picking the reading that lets the run proceed.
 [S2] `reproduction_failed` - see above.
 [S3] `access_gated` - the human deposit exists but is controlled-access or request-gated.
      Record what the request requires and stop. This is a stop on bar reading, not a STOP on
      the task: the dataset stays in `unverified` and the carried NO-GO stands.
A STOP is never resumed and never restarted with an adjustment.

== Carried forward from results/scopefib/PREREG.flag, verbatim and unchanged ==

Lines 10-87 of that file follow, byte-for-byte: bars A-E, clarifications C1-C4, definitions
D1-D7, the four readings with their exact wording, and the mandatory columns. They are not
restated, reinterpreted, softened or extended anywhere in this task. Where "this task" or
"the follow-up" appears in the carried text, it refers to SCOPEFIB as written; SCOPEFIB2
changes none of it.

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
== Outcome wording for the no-change case (added by this task; disclosed) ==

The four carried readings above are the only conclusions. This task reports which of them fires
for this dataset. The carried set has no wording for "the unverified entry could not be
resolved", and the task requires that case to be stated, so one line is fixed now for it. It is
not a fifth reading and it changes none of the carried four:

  "UNRESOLVED (SCOPEFIB2 scope statement, not a reading): the ADAM12 human Perturb-seq stays in
  `unverified`; <bars> could not be read because <reason>. SCOPEFIB's NO-GO stands unchanged."

If GO, PARTIAL or INCOMPLETE fires for this dataset, SCOPEFIB's NO-GO is superseded. That is
said in SCOPING_FIBROSIS2.md and marked there. SCOPING_FIBROSIS.md, PROGRESS_SCOPEFIB.md,
results/scopefib/ and every FINDINGS_*.md are not modified, per CLAUDE.md ("A reading that is
superseded is marked superseded, never edited") and the task.

== Honesty obligations carried from CLAUDE.md ==

- No statistic is computed beyond the cell counts bar C requires, so no chance baseline and no
  recovery power apply. Counts reported here are record inventories, not test results.
- A null here is a record-reading null, not an experimental one. What this task does not cover
  is stated explicitly: donor age of the fibroblasts, and whether any result on them transfers
  to aged cells.
- Every fact is resolved from its own record page or article text in this session. Nothing is
  recalled. Where a fact arrives through a tool that returns a model-written summary rather
  than raw text, that is marked in the source column.
- Nothing over 50 MB is downloaded.

== Process ==

Deliverables: results/scopefib2/PREREG.flag (this file), SCOPING_FIBROSIS2.md (this flag
reproduced verbatim at the top), PROGRESS_SCOPEFIB2.md. Raw query responses under
results/scopefib2/. Code at src/scopefib2_run.py. One local commit, no push.
```

## Reproduction of the stored SCOPEFIB record, before any query

Run as `python src/scopefib2_run.py reproduce`; written to `results/scopefib2/reproduction.csv`.
SCOPEFIB computed no statistic, so there is no numeric to reproduce within a tolerance. What
the flag fixed instead was reproduced exactly. 10 checks, 0 mismatches.

| check | expected | observed | match |
|---|---|---|---|
| R1 `results/scopefib/datasets.csv` data rows | 23 | 23 | yes |
| R2a readings in `PROGRESS_SCOPEFIB.md` | 2 | 2 | yes |
| R2b reading 1 also in `SCOPING_FIBROSIS.md` verbatim | True | True | yes |
| R2b reading 2 also in `SCOPING_FIBROSIS.md` verbatim | True | True | yes |
| R3 ADAM12 entry present in the `unverified` block | True | True | yes |
| R3b that block names GSE278990 as mouse only | True | True | yes |
| R4 carried lines 10–87 byte-identical in this task's flag | True | True | yes |
| R4b carried block line count | 78 | 78 | yes |

The two reading strings reproduced are the ones quoted above. No STOP fired.

I also re-resolved, in this session, the SCOPEFIB fact this task depends on: **GSE278990 holds
no human data.** Its own record gives `Series_sample_organism = Mus musculus`, `taxid 10090`,
and five SubSeries — GSE278974, GSE278975, GSE278980, GSE278983, GSE313761 — each of which I
fetched and each of which is `Mus musculus` only. So SCOPEFIB's reason for leaving this entry
unverified was correct, and the human data were never in GEO.

## How the deposit was located

The flag fixed the order: GSA first, then the article, then supplementary. Both worked. Every
request with its HTTP status is in `url_records.json` (54 requests: 42×200, 8×404, 3×403, 1×400)
and `browser_pane_urls.csv` (4 requests made through the browser pane). The route that found it:

1. `https://ngdc.cncb.ac.cn/gsa-human/search?searchTerm=…` → **404**. GSA-Human has no search of
   its own; its home page links out to NGDC BIG Search (`/search/specific?db=hra&q=`).
2. `https://ngdc.cncb.ac.cn/gsa/search?searchTerm=ADAM12` → 200, but GSA proper (the non-human
   archive) returns only two INSDC-mirrored mouse *Adam12* KO T-cell runs from an unrelated 2020
   study. This is why a GSA search alone would have missed it.
3. The BIG Search page is server-rendered around one JSON endpoint,
   `https://ngdc.cncb.ac.cn/search/api/specific?db=hra&q=ADAM12&sort=desc` → 200, one hit:
   **HRA014194**, BioProject **PRJCA048014**, *Homo sapiens*, title identical to the paper.
4. `https://ngdc.cncb.ac.cn/gsa-human/browse/HRA014194` → 200. Record read.
5. The article's own statement, read in the browser pane, names both repositories and says which
   data went where (verbatim in `article_excerpts.json`): "Data generated for this study can be
   downloaded from NCBI GEO (accession no. GSE278990) and GSA (accession no. PRJCA048014). …
   The raw FASTQ files of human scRNA-seq data in this study are released on Genome Sequence
   Archive in compliance with Chinese regulations, under accession no. PRJCA048014."

Sibling human accessions: I queried `db=hra`, `db=omix` and `db=gsa` for `PRJCA048014`,
`fibroblast checkpoint` and `Jianan Li`. **HRA014194 is the only one.** There is no OMIX
(processed-data) deposit and no CRA deposit for this BioProject, so there is no public processed
human matrix anywhere.

The article page itself: a direct urllib fetch of `https://www.cell.com/cancer-cell/fulltext/…`
returned **403**, and `https://www.sciencedirect.com/science/article/pii/S1535610825005525`
returned **400**, exactly as SCOPEFIB found. The same URL loads fine in the Claude desktop app's
built-in **browser pane**, where the article is Open Access. That is how the full text and the
supplementary workbooks were read; see §Deviations.

## The record, arm by arm

`results/scopefib2/hra014194_experiments.csv` has all 39 experiments. 101 runs, 202 FASTQ files.
Every experiment lists platform `Illumina NovaSeq X Plus` and release date 2025-12-16.

**The screen** — 6 experiments, all on one sample, `HRS2102648`, individual `HRI1628169`,
described in the record as "Mixed sample of adjacent-tumor sample from patient P1 and
adjacent-tumor and tumor-core sample from patient P4":

| experiment | title | strategy | runs | lanes |
|---|---|---|---|---|
| HRX2334923 | CRISPR_1_RNA | RNA-Seq | 4 | P0926-1, P0926-2 |
| HRX2334924 | CRISPR_2_RNA | RNA-Seq | 12 | P0926-3 … P0926-8 |
| HRX2334925 | CRISPR_1_CellPlex | Targeted-Capture | 2 | Gu0926-1,2 |
| HRX2334926 | CRISPR_2_CellPlex | Targeted-Capture | 6 | Gu0926-3 … 8 |
| HRX2334927 | CRISPR_1_sgRNA | Targeted-Capture | 2 | Ce0926-1,2 |
| HRX2334928 | CRISPR_2_sgRNA | Targeted-Capture | 6 | Ce0926-3 … 8 |

So: 8 10x lanes, each with a gene-expression, a CellPlex (CMO) and an sgRNA library. The record
does **not** say what the CMO tags encode. The methods say CMOs were assigned per cell and cells
with none or more than one were dropped; they do not list the multiplexed conditions. Given the
sample is a two-patient, three-region mix, patient and region are the obvious candidates, but
that is an inference and nothing in the record states it.

**The unperturbed state-definition arms** — 12 experiments, `HRX2334929`–`HRX2334940`, on six
different samples from three different individuals: P1 (`HRS2102642/3`), P2 (`HRS2102644/5`),
P3 (`HRS2102646/7`), each as `{NT, T} × {−, TGFb}`, e.g. `P1_NT` and `P1_NT_TGFb`. Their
sequencing-run names are all `20210825-LJN-*`, against `P0926-*` for the screen.

**Downstream validation** — 21 experiments, `HRX2363326`–`HRX2363346`, on P5 (`HRS2127694`:
DMSO / TGFb1 / IFNb1, three replicates each) and P6 (`HRS2127695`: EV, AAV, KO, OE, three each).
These are bulk arms, not part of the screen.

## Bar A, in detail — this is where it fails

The question D5 poses is precise: are the perturbed cells and a control-perturbation arm that
did **not** get the activating stimulus in the same experiment and the same batch?

**What the screen was run on: cells activated *after* perturbation, all of them.** The STAR
Methods section "Ex vivo CRISPRi and CRISPRa screenings" describes exactly one condition:

> For screening, PDFs were transduced with CRISPRi or CRISPRa lentivirus and their respective
> gRNA libraries at an optimal multiplicity of infection (MOI) of 0.5. The culture medium was
> refreshed 16 h later, and TGFβ was added at a final concentration of 10 ng/mL two days
> post-transduction. PDFs were cultured for another three days and then harvested.

There is no parallel no-TGF-β screen arm in that description, and `Table S1A` lists the screens
as exactly two rows — "CRISPRi screen" (34,890 cells) and "CRISPRa screen" (39,367 cells) — with
no condition split, while it lists the unperturbed cells as four rows that *do* split on TGF-β.
The Perturb-seq analysis section likewise scores cells "under TGFβ treatment".

**The one sentence that points the other way, disclosed rather than resolved away.** In the
Results, the design is described as: "CRISPRi/a-expressing cells were cultured under the same ex
vivo conditions as the aforementioned unperturbed cells and then subjected to scRNA-seq
harvesting." The unperturbed cells were "cultured and exposed to TGF-β or alone". Read
literally, "the same ex vivo conditions" could mean the perturbed cells were also split ±TGF-β.
I did not take that reading, for three reasons: the STAR Methods state a single condition with a
concentration and a schedule; Table S1A gives the screens one row each and the unperturbed cells
four; and the analysis section describes the scored cells as TGF-β-treated. Both texts are in
`article_excerpts.json` so the next reader can disagree. **If that sentence is right and the
screen did include an unstimulated control-perturbation arm, bar A could pass and this verdict
is wrong** — and the thing that would settle it is the CMO-to-condition mapping, which is in
neither the record nor the paper.

**And even on the most generous reading, C3 bites.** Suppose the unperturbed ±TGF-β arms are
accepted as the target. Then control and perturbed arms differ in library chemistry — 10x
Chromium **5′ kit** for all four unperturbed rows, 10x Chromium **3′ HT v3.1** for both screen
rows (Table S1A, and the STAR Methods name the 3′ HT kit for the screen) — and also in sample,
in individual (P1/P2/P3 vs the P1+P4 pool) and in sequencing-run name. C3 defines exactly this
as `split`, and `split` fails bar A. The reason C3 exists is in the carried flag: FINDINGS_DONOR
measured a transgene's direction agreeing with itself at 0.26 across batches against 0.50 for
two lines measured together, so a chemistry difference between target and perturbed arms can
exceed the effect being measured. A 5′/3′ difference is larger than a plate difference.

`perturbation_timing` = **before**. Transduction precedes TGF-β by two days, so this is
prevention, not reversal. Under C1 that makes it inventory only for the follow-up distance
analysis even if every other bar passed.

## Bar B, in detail

Counts in `results/scopefib2/library_counts.csv`, all resolved this session from Table S3
(`mmc4.xlsx`) and Table S4 (`mmc5.xlsx`), parsed in the browser pane:

| set | n |
|---|---|
| CRISPRi library, distinct gene targets | 590 |
| CRISPRa library, distinct gene targets | 265 |
| union of the two libraries (D1: one gene = one agent, any modality) | **667** |
| genes present in both libraries | 188 |
| targeting sgRNA rows, CRISPRi / CRISPRa | 2,450 / 1,150 |
| targets recovered after the screen's >5-cell filter, CRISPRi / CRISPRa | 432 / 210 |
| union of recovered targets | **532** |
| recovered in both arms | 110 |

The gene-list sheet and the sgRNA sheet agree on 590 / 265 / 667 / 188. Control guides
(AAVS1_1–47, NonTargetingControl_1–43) are excluded, as D1 requires. A few gene labels are spelt
differently between the two sheets (`ATFs3` vs `ATF3`), so the union is exact only up to those
label differences. 2,450/590 = 4.15 and 1,150/265 = 4.34 sgRNAs per gene, against "4 sgRNAs per
gene" in the text; the surplus rows are the control guides counted above.

Bar B passes on count by a wide margin, and D6 is satisfied because every perturbed cell
received TGF-β.

One correction to the task brief, which named "54 targets": the 54 is the **in silico**
prognostic shortlist (1,358 genes discriminating fibroblasts → 211 low in normal tissue → 54
with hazard ratio > 1 and p < 0.05), not the Perturb-seq library size. The brief's note that
"a screen run entirely on unstimulated cells fails D6 even at 54 targets" is therefore moot
twice over: the library is 667 genes, and the cells were stimulated.

## Bar C, in detail — unverified, and not resolvable from public files

The only per-target depth statement in the paper is a floor, not a median:

> For the inference of perturbation effects, we only used cells detected with sgRNAs targeting
> one single gene and retained those targets detected in more than 5 cells.

A floor of 6 cells is far below D7's median of 30, and a floor says nothing about the median.
What can be computed from stated totals is a mean, and only as an estimate:

- CRISPRi: 34,890 cells ÷ 432 recovered targets = **80.8 cells/target (estimate)**
- CRISPRa: 39,367 cells ÷ 210 recovered targets = **187.5 cells/target (estimate)**

These are upper-leaning estimates, not D7's statistic: the numerators are all retained
high-quality cells, including cells with no sgRNA, more than one sgRNA, or a control guide, none
of which belong to a target. The real median could be well below 30 while the mean is 80, since
a >5-cell floor with a long tail is exactly the shape that pulls a mean above a median.

**What would settle it, and why it cannot be opened.** A per-cell guide assignment. That exists
only inside the controlled-access FASTQ data at HRA014194. I checked the alternatives and none
carries a per-target cell count: Table S4D (156,651 rows) gives log2FoldChange per perturbation
per gene-module gene; Table S4E (3,213 rows) gives the mean per programme; the study's GitHub
repository `Qile-Guo/Fibroblast` holds seven figure scripts and a README and no data; and there
is no OMIX or CRA deposit. So bar C stays `unverified` and stays that way for any reader without
DAC approval.

Because bar A fails, resolving C could not change any reading. I say that explicitly so the
`unverified` is not mistaken for a gap that a download would close.

## Access — controlled, and nothing was requested

- **Accessibility**: the HRA014194 record states `Data Accessibility: Controlled access`, with a
  "Request Data" action.
- **Committee**: DAC no. **HDAC007040**, "ADAM12 DAC". Contact Li Jianan, lijianan@pku.edu.cn.
  Its stated remit: "This Data Access Committee (DAC), composed of members from the
  collaborating institutions that generated the data, governs all access to the patient-derived
  fibroblast (PDF) samples and the associated single-cell data."
- **Submitter**: Li Jianan, Peking University. Submission date 2025-10-27; release 2025-12-16.
  The record showed `Requests: 1` at the time it was read.
- **What a request requires**, from the GSA-Human policy page: a registered GSA-Human account;
  an application through the GSA-Human system naming the research scope; approval by that DAC,
  which has "the right to Agree or Disagree"; and acceptance of the data-user responsibilities —
  use confined to the promised research scope and to the applicant's own group and
  collaborators, no redistribution, no attempt to identify individuals, accession cited in any
  publication, responsibility extending to the head of the applicant's organisation, and
  NGDC publishing the dataset ID, applicant name, affiliation, usage period and research
  programme as statistics.
- **No access request was submitted in this session**, and none should be read into this file.
  Deciding whether to apply is outside a scoping task.

## Deviations and decisions not in the flag (disclosed)

1. **The article was read through the Claude desktop app's built-in browser pane**, not with a
   plain HTTP fetch, because `www.cell.com` returns 403 to the fetcher and ScienceDirect returns
   400. The flag anticipated a "publisher mirror that SCOPEFIB found reachable"; the browser
   pane is not a mirror, it is the same URL from a different client. The text is the published
   Open Access article. Consequence: the article text is **not** stored verbatim in this
   directory the way GEO and GSA text is. What is stored is `article_excerpts.json` — the
   passages each bar rests on, quoted verbatim with their character offsets into the article's
   `innerText` (total length 129,332), plus the URL and how it was read. A reader who wants the
   whole text must open the URL. Browser-pane requests are listed in `browser_pane_urls.csv`
   rather than `url_records.json`, because they did not go through the logging fetcher.
2. **Three supplementary workbooks were parsed inside the browser page**, not downloaded, for
   the same 403 reason: `mmc2.xlsx` (Table S1, 18,108 bytes), `mmc4.xlsx` (Table S3, 125,281
   bytes), `mmc5.xlsx` (Table S4, 5,449,018 bytes). All well under 50 MB. I unzipped and parsed
   them with JavaScript in the page and brought back counts and the specific rows quoted in
   `article_excerpts.json`, not the whole sheets. Nothing was written to disk from them. The
   counts are therefore reproducible only by repeating that parse; the parsing logic is
   described in `article_excerpts.json` and the row-level evidence is quoted there.
3. **One file was downloaded**: `results/scopefib2/opened/HRA014194.xlsx` (12,956 bytes), the
   study's own public metadata workbook, from the GSA-Human record page. Its page warns that
   "The detailed metadata can only be downloaded by authorized users, while other users can
   access partial information", and indeed it contains the same two tables the web record shows
   (individuals & samples, files) and no library protocol. It added nothing beyond the AJAX
   endpoints, and is kept as provenance. It is untrusted downloaded input: it sits in its own
   directory and was parsed by a script kept elsewhere.
4. **Additions beyond the flag**: `src/scopefib2_run.py` (logging fetcher with a `reproduce`
   mode, a `file` mode for binary record files, and UTF-8 console reconfiguration — the Windows
   console here is cp949 and SCOPEFIB lost a query batch to exactly that), and
   `results/scopefib2/_tables.py` (hand-entered table → CSVs). Two WebSearch queries were run
   for the article's data-availability statement before the browser pane was tried; both
   returned nothing usable and nothing in this file rests on a search snippet.
5. **`target_match` is recorded as `n/a` with an explanation**, not as one of C2's three values.
   C2 offers `isogenic`, `cross_donor`, `unstated`. None fits a screen with no in-run
   unstimulated arm at all: there is no control arm whose donor identity could be compared. The
   full string in `datasets.csv` names what the nearest unstimulated cells are instead (P1
   isogenic to part of the screen pool, P2 and P3 cross-donor) and the donor counts (0 healthy,
   8 diseased PDF donors). This is the same kind of gap SCOPEFIB disclosed when it wrote
   `perturbation_timing = n/a` for datasets with nothing to time.
6. **An outcome line for the no-change case was fixed in this task's flag** before any query,
   because the four carried readings have no wording for "the entry could not be resolved" and
   the task required that case to be stated. In the event it was not needed: the bars were read.
   It changes none of the carried four and is labelled a scope statement, not a reading.
7. **Two record inconsistencies noticed and not resolved**, because neither can move a bar:
   - Table S1A's accession column attributes the six **human** ex vivo fibroblast rows to
     "GEO: GSE278990". GSE278990 and all five of its SubSeries are `Mus musculus` only, as I
     re-checked this session. The data-availability statement is the one that matches the
     archives: human raw FASTQ are at GSA under PRJCA048014. A reader following Table S1A alone
     would look for human data in GEO and not find it — which is very likely what happened to
     SCOPEFIB.
   - The GSA record's description says "We first defined fibroblast states by scRNA-seq of
     unperturbed PDFs under TGFb stimulation. We then introduced optimized CRISPR perturbations
     into this model", which reads as the screen being run in the TGF-β-stimulated model and is
     consistent with the STAR Methods, but is a submitter summary rather than a protocol.
8. **The CMO-to-condition mapping was not found** in the record, the article, the supplementary
   tables I opened, or the GitHub repository. It is the single thing that would turn the bar A
   ambiguity in §Bar A into a fact. I did not open Tables S5 and S6 (`mmc6.xlsx` 7.41 MB,
   `mmc7.xlsx` 5.49 MB), Data S1 (`mmc8.zip`, 49.73 KB), or the two PDFs (`mmc1.pdf` 5.20 MB,
   `mmc9.pdf` 36.47 MB — the latter is within the 50 MB limit). Since bar A fails on C3 `split`
   even under the generous reading, none of these can change the reading; but `mmc1.pdf`
   (supplementary figures, which include Figure S1's screen schematic) is the file I would open
   first if someone wanted the CMO mapping settled.

## Honesty notes the flag and CLAUDE.md require

- **No statistic was computed.** The only arithmetic here is division for the two mean-cells
  estimates and set sizes for the gene counts. There is therefore no chance baseline and no
  recovery power to report: nothing here is a test result, and no count is being compared
  against a null.
- **This is a record-reading failure, not an experimental null.** Bar A failing means the
  experiment as deposited cannot pose the distance-to-target question. It says nothing about
  whether perturbations in this system would move cells toward an unactivated state.
- **A superseded reading would have been marked, not edited.** Nothing was superseded.
  `SCOPING_FIBROSIS.md`, `PROGRESS_SCOPEFIB.md`, `results/scopefib/` and every `FINDINGS_*.md`
  are untouched; `git status` on the commit for this task shows only new files plus this one.
- **Circularity**: none to disclose. No axis, basis or reference was fitted or applied here.
- **The chance-baseline and cell-count warnings in CLAUDE.md** apply to shortlists; the one
  shortlist in play, the paper's 54 prognostic genes, is the paper's, is not used for any bar,
  and is reported only to correct the task brief's reading of it.

## near_miss (this task's one dataset)

| accession | aimed | exact shortfall |
|---|---|---|
| HRA014194 (PRJCA048014) | A+B+C+D+E | **A fails.** 667 library genes / 532 recovered targets applied to TGF-β-activated patient-derived fibroblasts (B passes, D6 satisfied), published directional readout with named hits (D passes), human adult donors (E passes) — but no unstimulated control-perturbation arm inside the screen, and the unperturbed ±TGF-β arms are a separate experiment on 10x 5′ chemistry against the screen's 3′ HT v3.1, which is C3 `split`. C is additionally `unverified` (floor of >5 cells per target; no median; per-cell data controlled access). Timing `before`. |

## leads (seen incidentally, not opened, not in scope)

Recorded only as accessions, per the flag. These came back from `ADAM12` queries against NGDC
databases and are other groups' studies:

- GSA / INSDC mirror: `SRX8431572`, `SRX8431573`, `SRR11883857`, `SRR11883858` (mouse T-cell
  *Adam12* KO, SRA1081565 / SRP265377, released 2020-07-31).
- BioProject: `PRJNA852211`, `PRJNA852197` (ADAM12⁺ MSCs in B16-OVA melanoma, Institut Pasteur),
  `PRJNA565443` (ADAM12 KO mouse heart failure), `PRJNA327112` (ADAM12 knockdown in SUM159PT
  breast cancer cells), `PRJNA437473` (Adam12 / lnc015192 ceRNA in breast cancer).
- BioSample: `SAMEA115610216`, `SAMN15065521`, `SAMN15065522`.
- The study's own mouse siblings, in scope but mouse and so excluded from any pass by D3:
  `PRJNA1169870` (= GSE278990), `PRJNA1169876`, `PRJNA1381552` (= GSE313761).

None was opened beyond the one-line search result above.

## What this does not cover

- **Search recall is SCOPEFIB's, unchanged.** This task extended nothing. The limits listed
  under "What this NO-GO does not cover" in `SCOPING_FIBROSIS.md` stand as written: GEO
  `esearch` capped at 100 hits per query, one ArrayExpress query, and no systematic search of
  EGA, Synapse, CELLxGENE, Zenodo/figshare or the Broad Single Cell Portal. One recall limit is
  now narrower in one specific way: **GSA-Human was not searched by SCOPEFIB and now has been,
  for this study only.** Other Chinese-deposited human fibroblast screens would still be missed
  by both tasks, and `db=hra` was queried with three study-specific terms, not with
  disease-area or assay terms.
- **Donor age, and what the adult donors do and do not buy.** This is the first human dataset in
  either task whose fibroblasts are not fetal, neonatal or unstated: the 8 PDF donors are aged
  35, 79, 74, 54, 59, 68, 50 and 61. But the **screen** ran on a pool of two donors, aged 35 (F)
  and 54 (M). So a result from this screen would be a result on middle-aged colorectal-cancer
  fibroblasts from two people, not on aged cells, and not separable from those two individuals.
  The older donors (P2 79, P3 74, P9 86, P11 77) appear in the unperturbed and validation arms,
  not in the screen.
- **Disease context.** These are cancer-associated fibroblasts from colorectal cancer, activated
  by TGF-β in culture. Nothing here transfers to idiopathic pulmonary fibrosis, cardiac or liver
  fibrosis, or to fibroblast activation in aged tissue without its own evidence.
- **Bar A rests on a reading of two passages that do not fully agree**, as §Bar A says plainly.
  The CMO mapping would settle it and was not found.
- **Bar C cannot be closed by anyone without DAC approval**, so a future task cannot resolve it
  by opening a file.

## Next

None within this task (one task per session). The record now shows, for the fibrosis
distance-to-target question:

1. The mouse GSE261783 route, unchanged from SCOPEFIB: needs an ortholog map, a new
   pre-registration, and accepts `before` timing.
2. HRA014194 is the closest human screen that exists, and it is one design change away: the same
   screen with a no-TGF-β control-perturbation arm in the same lanes would pass A and B
   together. It does not have one.
3. If someone wants to pursue HRA014194 anyway, the two things to settle first are the CMO
   mapping (`mmc1.pdf`, Figure S1) and DAC approval for the per-cell data. Both are decisions
   outside a scoping session, and neither is taken here.
