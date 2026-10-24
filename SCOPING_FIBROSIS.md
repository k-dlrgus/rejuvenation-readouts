# SCOPING_FIBROSIS: is there a public fibrosis dataset where the distance-to-target question can be asked?

**Task:** SCOPEFIB. **Status:** search complete 2026-10-24. Fit nothing, compute no statistic, download nothing over 50 MB. Does not modify any `FINDINGS_*.md`.

**Readings that fired (in precedence order; wording as fixed in the flag):**

1. `INCOMPLETE: GSE261783 passes A+B; E fails or is unverified. Not GO, not NO-GO.`
2. `NO-GO: no human dataset passes A+B. Closest near-miss: GSE305412; single missing ingredient: perturbation count (5 distinct knockdown agents applied across TGF-β activation; bar B needs ≥20).`

GO and PARTIAL did not fire. The only dataset that passes A+B is mouse. In the human records resolved here, the in-run unactivated target (bar A) and ≥20 perturbations of the activated state (bar B) never occur in the same experiment. Datasets with an isogenic in-run target have 0–5 perturbations. Datasets with ≥20 perturbations either have no activated state (LINCS fibroblast lines, Ginkgo GDPx2) or have no healthy or unactivated target in the run (GSE227391, keloid fibroblasts).

Machine-readable table: `results/scopefib/datasets.csv` (23 rows, all mandatory columns). Raw record text: `results/scopefib/geo_records.json`, `url_records.json`, `ncbi_searches.json`, `lincs_*.csv`, `sigcom_fibroblast_signatures.csv`. Fetcher: `src/scopefib_run.py`.

## Pre-registration (verbatim copy of `results/scopefib/PREREG.flag`)

```
SCOPEFIB pre-registration flag
Written 2026-10-24T10:06:33+09:00 (2026-10-24T01:06:33Z), before any database query of any kind
(no GEO / SRA / ENA / ArrayExpress / CELLxGENE / LINCS / clue.io / PubMed / web search has been run in this task).
Nothing below is re-tuned after any record is seen.

== Question (verbatim from the task) ==

Is there a public dataset in which fibroblast activation / fibrosis can be asked the distance-to-target question — does a perturbation move cells closer to the unactivated state, with the target measured in the same experiment?

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

== Process ==

Deliverables: results/scopefib/PREREG.flag (this file), SCOPING_FIBROSIS.md (this flag reproduced verbatim at the top), PROGRESS_SCOPEFIB.md. Raw query responses saved under results/scopefib/. One local commit, no push.
```

## Deviations and decisions not in the flag (disclosed)

1. **The flag was edited once after it was first written, before any query.** I added the header line that separates D5–D7 from D1–D4. It says D5–D7 were written by the session and were not shown to the user in chat beforehand. The original header would have implied the user saw them. No bar, reading or definition text changed. The timestamp is the time of the first write (10:06:33 +09:00). The edit left the file mtime at 10:07:07.99 +09:00. The first database query was logged at 01:07:46Z = 10:07:46 +09:00 (`ncbi_searches.json`, first `fetched_utc`).
2. **`perturbation_timing = n/a`.** The flag allows `before / during / after / unstated`. Where a dataset applies **no** perturbation to an activated state (activation atlases, LINCS, GDPx2), there is nothing to time, so I wrote `n/a` rather than `unstated`. For primary cells from diseased tissue I wrote `after (cells diseased at isolation)`. That follows from the design (the perturbation is applied to cells already diseased), not from a sentence in the record.
3. **Two definitional questions were left open because they cannot change any reading:**
   - whether CRISPRi of regulatory elements counts as "gene targets" under D1 (GSE243756);
   - whether a targeted probe-panel RNA readout (PHDs-seq, GSE227391) is "transcriptomic" (D2 covers only L1000).

   Both datasets fail bar A on other grounds, so the answers would not move them. These questions would need a ruling before a future search could pass either dataset.
4. **LINCS 2020 signature counts come from the SigCom LINCS metadata API** (Ma'ayan lab mirror, `signatures/find`), not from clue.io's `siginfo_beta.txt`. That file is 465 MB, over the 50 MB limit. `cellinfo_beta.txt` (38 KB) was read directly.
5. **Some record facts came through WebFetch**, which returns a model-written summary of the page, not the raw text. These are: the Ginkgo GDPx2 portal page, the GDPx2 bioRxiv full text, and the Nature Communications article for GSE261742 (target count, timing, data-availability text). They are marked in the `source` column of `datasets.csv`. GEO, LINCS and SigCom facts are from raw text stored in this directory.
6. **Additions beyond the flag:**
   - `src/scopefib_run.py`: a metadata fetcher, plus count-only LINCS/SigCom inventories.
   - `results/scopefib/_datasets_table.py`: hand-entered table → `datasets.csv`.
   - WebSearch queries for literature leads. Those results are leads, not records; nothing in a table rests on a search snippet.
7. **Untracked raw inputs.** The GEO LINCS metadata files (`GSE92742_Broad_LINCS_cell_info.txt.gz`, `…sig_info.txt.gz` 11.1 MB, `GSE70138_…cell_info_2017-04-28.txt.gz`, `…sig_info_2017-03-06.txt.gz` 1.9 MB) were downloaded to the session scratchpad, not the repo. `python src/scopefib_run.py lincs <dir>` re-downloads them from `ftp.ncbi.nlm.nih.gov/geo/series/GSE92nnn/GSE92742/suppl/` and `…/GSE70nnn/GSE70138/suppl/`. The repo already holds copies under `results/survey/opened/lincs/` from the SURVEY task; those were not read or reused.
8. `url_records.json` is 3.2 MB of record text (mostly GSM brief dumps). It is text, not a regenerable binary, and is tracked as provenance.
9. **One query batch crashed partway** on a console encoding error. Query 3 had already been logged; queries 4–6 were re-run. No query appears twice in `ncbi_searches.json`.
10. **Zenodo record 10.5281/zenodo.14794723** (GSE261742 data bundle) was requested once. The API call returned non-JSON and was not retried. It would only bear on bar C for the mouse dataset, and that cannot change the reading (see §Readings).

## Readings, with the conditions that fired them

| reading | fired | dataset(s) | why |
|---|---|---|---|
| GO | no | — | no human dataset passes A+B |
| PARTIAL | no | — | no human dataset passes A+B |
| INCOMPLETE | **yes** | GSE261783 (mouse) | A pass (isogenic, same_series), B pass (30 chromatin-factor KOs, timing `before`); C unverified; D pass; **E fails** |
| NO-GO | **yes** | — | no human dataset passes A+B |

**INCOMPLETE, single resolver:** `mouse only`. No file can resolve E. Bar C is also unverified for this dataset. The design aimed at 200 cells per sgRNA; the achieved count is not in the GEO record or the article text read. The per-cell guide assignments in `GSE261783_RAW.tar` or Zenodo 14794723 would settle C. They were not opened because settling C cannot change the reading.

**NO-GO, closest near-miss: GSE305412.** I chose it because it misses only bar B, and misses it on a count, while it passes A, the non-negotiable bar:
- **Design:** IMR-90; shSCR and five knockdown agents; each in an Unactivated arm and an Activated arm (TGF-β 2 ng/mL plus serum starvation, 24 h).
- **Bar A:** passes. The shSCR Unactivated target and every sh-perturbed arm were sequenced on the same instrument.
- **What disqualifies it anyway:**
  - It has 5 agents, and bar B needs 20.
  - Its timing is `before`, so under C1 even a ≥20-agent version would be inventory only for the follow-up.
  - IMR-90 is a fetal line.

## Tables (mandatory columns; full text in `datasets.csv`)

Abbreviations: `iso` = isogenic, `xd` = cross_donor (healthy/diseased donors), `ss` = same_series, `st` = stated.

### Pass table (human, A+B+C+E)

None. No human dataset passes A+B, so there are no per-arm chemistry / instrument / submission-date rows to report here. The per-arm record for the one dataset that passes A+B, GSE261783 (mouse), is in the non-human block.

### (1) TGF-β-induced activation of human dermal or lung fibroblasts

| accession | organism | cell type / source | design | n pert (D1) | readout | units / pert | timing | target_match | batch | A | B | C | D | E |
|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|
| GSE305412 | human | IMR-90 fetal lung line | WT, shSCR + 5 KD agents × Unactivated / Activated (TGF-β 2 ng/mL, 24 h) | 5 | bulk RNA-seq | 2–5 reps | before | iso | ss (WT arms split: NextSeq 1000 vs NovaSeq X Plus; shSCR target not split) | pass | **fail (5)** | pass | n.a. | pass |
| GSE342214 | human | primary dermal fibroblasts | activation by LPS / TGF-β / CRISPRa IL6, ACTA2 / PBMC / CM; matched controls | 0 on activated state | bulk RNA-seq | n/a | n/a | iso | ss | pass | fail (0) | pass | n.a. | pass |
| GSE330910 | human | skin fibroblasts (immortalized + primary) | shRUNX2 vs shControl ± TGF-β; TGF-β ± F0565-0303 | 2 | bulk RNA-seq | bulk | before / during | iso | ss | pass | fail (2) | pass | n.a. | pass |
| GSE300230 | human | dermal fibroblasts, 10 y and 75 y | 3D gel, tension ± TGF-β | 0 | bulk RNA-seq | n/a | n/a | iso | ss | pass | fail (0) | pass | n.a. | pass |
| GSE243756 | human | HFF, dCas9-KRAB | CRISPRi of 87 mechano-enhancers on rigid TCP | 87 elements (D1 open) | scRNA-seq | mean 159 cells/gRNA (median not stated) | unstated | unstated | ss | **fail** (no soft/unactivated arm) | not resolved | unverified | n.a. | pass |
| GSE306853 | human | PDAC CAF lines | Perturb-seq of 5 genes, mono vs co-culture | 5 | scRNA-seq | not stated | unstated | unstated | ss | fail | fail (5) | unverified | n.a. | pass |
| GSE308912 | human | engineered intestinal tissue ± fibroblasts | 4 inhibitors | 4 | bulk RNA-seq | bulk | unstated | unstated | ss | fail | fail (4) | pass | n.a. | pass |
| GSE294260 | human | lung fibroblasts + SAEC co-culture | injury models with TGF-β | 0 | bulk RNA-seq | n/a | n/a | unstated | ss | n.r. | fail (0) | pass | n.a. | pass |
| GSE152250 | human + rodent | primary human cells ± TGF-β; rodent UUO/CCl4/BDL | TGF-β vs control only | 0 | bulk RNA-seq | n/a | n/a | iso | ss | pass | fail (0) | pass | n.a. | pass (human arm) |
| GDPx2 (Ginkgo) | human | primary dermal fibroblasts | 85–89 compounds × 6 conc. on **unstimulated** cells; DMSO every plate | 85–89, none on activated state | DRUG-seq | 4 reps | n/a | n/a | st | **fail** (no activated state) | fail (D6) | pass | n.a. | pass |

### (2) Primary fibroblasts from fibrotic human tissue

| accession | organism | cell type / source | design | n pert (D1) | readout | units / pert | timing | target_match | batch | A | B | C | D | E |
|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|
| GSE328392 | human | lung fibroblasts, healthy + IPF | irradiation time course; serum starvation; bleomycin | 3 | bulk RNA-seq | bulk | after | xd 7/9 | ss | pass | fail (3) | pass | n.a. | pass |
| GSE289407 | human | SSc + healthy dermal fibroblasts, 3D rafts | nintedanib vs vehicle | 1 | bulk RNA-seq | bulk | after | xd (counts not in record) | ss | pass | fail (1) | pass | n.a. | pass |
| GSE141905 | human | IPF lung fibroblasts | TGF-β1, cmp4, both | 1 | microarray | bulk | unstated | unstated | ss | fail | fail (1) | pass | n.a. | pass |
| GSE227391 | human | keloid fibroblasts only | 128 compounds in adipogenic MDI medium | 128 | targeted probe RNA-seq (PHDs-seq) | 1 well | after | unstated (no healthy cells) | ss | **fail** (no target) | count ≥20; readout question open | pass | hit ABT869 (toward adipocyte) | pass |
| GSE280438 | human | cardiac fibroblasts, 2 patients | Kat5 inhibitor vs vehicle, then TGF-β | 1 | bulk RNA-seq | bulk | before | unstated | ss | n.r. | fail (1) | pass | n.a. | pass |
| GSE276445 | human | 6 primary fibroblast lines | vehicle / TGF-β / TGF-β + artesunate | 1 | bulk RNA-seq | bulk | unstated | iso | ss | pass | fail (1) | pass | n.a. | pass |

### (3) Hepatic stellate cell activation

| accession | organism | cell type / source | design | n pert (D1) | readout | units / pert | timing | target_match | batch | A | B | C | D | E |
|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|
| GSE266939 | human | primary HSC, 1 donor | RNP KO of 11 TFs + COL1A1 vs NT, 12 d | 12 | bulk RNA-seq | 3 reps | unstated | unstated (no quiescent arm) | ss | fail | fail (12) | pass | n.a. | pass |
| GSE78853 | human | HSC myofibroblasts | nortriptyline ± TGF-β; ceramide | 2 | bulk RNA-seq | 2 reps | unstated | unstated (no quiescent arm) | **split** (ceramide arms 2017, HiSeq 2000) | fail | fail (2) | pass | n.a. | pass |
| GSE256398 | human + mouse | liver snRNA-seq, patients vs controls | no perturbation | 0 | snRNA-seq | n/a | n/a | xd 6/20 | ss | pass | fail (0) | n/a | n.a. | pass |

### (4) Myofibroblast-reversal screens

No reversal screen with a transcriptomic readout of ≥20 agents was found:
- The GSE78853 study's upstream screen (1,600 compounds) and the GSE141905 study's upstream screen (FDA/EMA library) used imaging readouts. Their GEO records hold 1–2 follow-up compounds.
- The GSE152250 compound screen is described as phenotypic.
- GSE227391 is the only ≥20-compound transcriptomic screen on human diseased fibroblasts. Its direction is fibroblast → adipocyte, and it has no healthy target in the run.

### Non-human block (never counted toward a pass)

| accession | organism | cell type / source | design | n pert (D1) | readout | units / pert | timing | target_match | batch | A | B | C | D | E |
|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|
| GSE261783 (SuperSeries GSE261742) | mouse | primary Cas9 cardiac fibroblasts | pooled KO of chromatin factors; Quiescent vs Activated (OP1L); resting / IL-1β / TGF-β (OP2) | 30 (article text) | 10x 3' v3.1 scRNA-seq | designed 200 cells/sgRNA; achieved not stated | before | iso (NTC in every condition) | ss | pass | pass | **unverified** | pass | **fail** |
| GSE278990 (ADAM12, mouse SubSeries) | mouse | tumour models | in vivo KO / tumour | 0 screens | sc/bulk RNA-seq | n/a | n/a | n/a | n/a | fail | fail | n/a | n/a | fail |

GSE261783 per arm: all 10 GSMs list the Chromium Next GEM 3' v3.1 Dual Index kit, NextSeq 2000, and submission date Mar 17 2024. There is no difference between the control (Quiescent / resting) and perturbed (Activated / TGF-β) arms.

Caveats for GSE261783:
- The target count (30) comes from the article, not from GEO.
- Whether experiments OP1L and OP2 used the same library is not stated.
- The perturbation is a knockout made before stimulation. That is prevention, not reversal, so under C1 the follow-up would treat it as inventory only even if it were human.

## near_miss

| accession | aimed | exact shortfall |
|---|---|---|
| GSE305412 | A+B | 5 distinct knockdown agents; needs ≥20. A passes (isogenic, same instrument for target and perturbed arms). Timing `before`. |
| GSE227391 | A+B | 128 compounds on keloid fibroblasts, but no healthy or unactivated fibroblasts in the run (A fails). Readout is a targeted probe panel; whether that counts as transcriptomic was not pre-defined. |
| LINCS 2020 LUNG001 (IMR90 / WI38 / HFL1) | A+B | 319 distinct compounds per line, vehicle on the plate, but no activated state: no TGF-β or other stimulus in any signature. Perturbations are on unstimulated cells (D6). |
| Ginkgo GDPx2 | A+B | 85–89 compounds × 6 concentrations on primary dermal fibroblasts, but on unstimulated cells; no activated state. |
| GSE266939 | A+B | 12 KOs in primary HSC; needs ≥20, and there is no quiescent HSC arm. |
| GSE328392 | B | cross_donor 7 healthy / 9 IPF (≥3 per side), but 3 perturbations, all senescence inducers. |

## unverified

- **ADAM12 human patient-derived fibroblast CRISPRi/CRISPRa Perturb-seq** (Cancer Cell 2026, doi 10.1016/j.ccell.2025.12.018; Europe PMC id 41544628). The abstract describes an IFN-I program antagonising TGF-β-driven myofibroblast activation. The GEO SuperSeries GSE278990 contains mouse data only. No accession for the human screen was located, and the article page returned HTTP 403. Bars A–C cannot be read. **What would settle it:** the article's data-availability / key-resources table.
- **GSE261783, bar C** (see above).

Seen in search results but not opened (titles and sample counts only, so not assessed): GSE317322, GSE293962, GSE298860, GSE299164, GSE272972, GSE325500, GSE299851, GSE299852, GSE290015, GSE287614, GSE269412, E-MEXP-2750, E-MTAB-15932. Their titles name one agent each. That is an inference from titles, not a resolved record.

## LINCS / L1000 fibroblast arms (SCOPING_TISSUE left this unresolved)

What they contain:
- **Phase I / II GEO (GSE92742, GSE70138).** 98 cell ids in each `cell_info`. The only fibroblast line is **MCH58** ("skin fibroblast"). It has 5 perturbagen signatures on plate T2D001 at 144 h: tunicamycin plus shRNA against EIF2AK3, LRPPRC, TCF7L2 and WFS1. This agrees with SURVEY_TF.md. The fibroblast regex also matched `NPC`, `NPC.CAS9` and `NPC.TAK` ("normal stem fibroblast-derived iPScs"). Those are neural progenitors and are excluded.
- **LINCS 2020 (clue.io beta).** `cellinfo_beta.txt` adds **IMR90, WI38, HFL1** (lung fibroblast lines; IMR90 donor age "16 weeks", WI38 "3 months") and **1HAE** ("normal fibroblast sample"). SigCom lists 347 level-5 signatures for each lung line and 348 for 1HAE, all on plate **LUNG001**, all chemical, 24 h, mostly 10 µM. That is 319 distinct compounds per lung line and 320 for 1HAE.
- **TGF-β in LINCS.** TGF-β ligand (`trt_lig`: TGFB1/2/3) and TGF-β-pathway shRNA/ORF/CRISPR signatures exist only in non-fibroblast lines (A375, A549, HA1E, HCC515, HEPG2, HT29, MCF7, PC3, VCAP and their Phase II variants); see `lincs_tgfb_signatures.csv`.

What they do not contain:
- any TGF-β-stimulated (activated) fibroblast arm;
- any compound applied on top of an activating stimulus;
- any diseased-donor fibroblast.

So no LINCS fibroblast arm can pose the distance-to-target question in the flag's sense. Two further limits apply:
- L1000 measures 978 genes and infers the rest.
- The lung lines are fetal.

Not checked: `siginfo_beta.txt` itself (465 MB), and therefore whether LUNG001 has extra non-signature (level 3) wells.

## Frozen fibroblast coordinate system (23,485-gene GTEx ruler)

**Which passing datasets could be scored without a new reference panel: none.** No dataset reaches GO or PARTIAL, so there are no passing datasets to score, and no gene overlap is reported or estimated (D4 applies only to passing datasets).

The one INCOMPLETE dataset, GSE261783, is mouse. Scoring it in the human GTEx-derived ruler would first need a mouse-to-human ortholog map. That would be a new instrument decision, not a reading of the frozen one.

## What this NO-GO does not cover

- **Search recall.**
  - GEO `esearch` returned the first 100 hits per query; several queries had more (191, 214, 263, 511, 599). Hits were screened by title and sample count before opening.
  - ArrayExpress was queried once (8 hits). EGA, Synapse, CELLxGENE, Zenodo/figshare and the Broad Single Cell Portal were not searched systematically.
  - Paywalled articles (Cancer Cell 403) and supplementary tables were mostly not read.
  - A dataset that exists but was not returned by these queries is outside this file.
- **Bar B "applied to the activated state" (D6) is the binding constraint for the large compound sets** (LINCS LUNG001, GDPx2). If the follow-up instead wanted "does a compound move unstimulated fibroblasts toward or away from an activated reference", those sets would be candidates. That is a different question from the flag's, and pursuing it here would relax a bar after seeing what exists.
- **Aged cells.** Every human near-miss uses fetal (IMR-90, WI38, HFL1), neonatal (HFF) or unstated-age cells, except the GSE328392 donors (29–75). A null or a hit on these does not transfer to aged fibroblasts.
- **No statistic was computed**, so no recovery-power or chance-baseline statement applies. The counts here are record inventories, not test results.

## Next

None within this task (one task per session). If the distance-to-target question is to be run on fibrosis data, the record shows two routes. Neither is chosen here:
1. Use the mouse GSE261783 with an ortholog map. This needs a new pre-registration and accepts `before` timing.
2. Wait for or generate a human ≥20-agent screen applied **after** TGF-β activation, with an unstimulated control-perturbation arm in the same run.
