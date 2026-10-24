# PROGRESS_SCOPEFIB

**STOP status:** complete. No STOP key fired.

**Readings fired:**
- `INCOMPLETE: GSE261783 passes A+B; E fails or is unverified. Not GO, not NO-GO.`
- `NO-GO: no human dataset passes A+B. Closest near-miss: GSE305412; single missing ingredient: perturbation count (5 distinct knockdown agents applied across TGF-β activation; bar B needs ≥20).`

**Single next action:** none in this task. Picking a route (mouse GSE261783 with an ortholog map, or a new human after-activation screen) is a decision outside this session and needs its own pre-registration.

## Record of what was decided when

| time (+09:00, 2026-10-24) | event |
|---|---|
| before 10:06 | Read CLAUDE.md, SCOPING_TISSUE.md, PROGRESS_SCOPING.md, STATE.md; grep for the 23,485-gene panel. No database query. |
| before 10:06 | Four ambiguities raised in chat (B timing, A identity, A batch evidence, outcome gap). The user ruled on each (C1–C4). D1–D4 proposed in chat; not overridden. |
| 10:06:33 | `results/scopefib/PREREG.flag` written. |
| 10:07:07 | Flag edited once: header added stating D5–D7 were not shown in chat. No bar, reading or definition text changed. |
| 10:07:46 | First database query (NCBI GDS esearch). |
| 10:07–10:20 | GEO searches in search order 1→4, LINCS GEO + LINCS 2020 + SigCom, ArrayExpress, Europe PMC, WebSearch/WebFetch for literature and portals. |
| ~10:21 | `datasets.csv` written; SCOPING_FIBROSIS.md written with the flag embedded verbatim (checked by string match). |

## Searches run (do not repeat)

1. `ncbi_searches.json`: 17 GDS esearch + esummary queries, each with its `fetched_utc`.
2. `geo_records.json`: GEO series text for 32 accessions.
3. `url_records.json`: GSM brief dumps (GSE261783, GSE266939, GSE305412, GSE227391, GSE78853, GSE152250), GPL20795, LINCS 2020 `cellinfo_beta.txt`, ArrayExpress search, Europe PMC search.
4. `lincs_fibroblast_cells.csv`, `lincs_fibroblast_arms.csv`, `lincs_tgfb_signatures.csv`: from GEO GSE92742 / GSE70138 `cell_info` + `sig_info` (raw files in the session scratchpad, not tracked; re-download with `python src/scopefib_run.py lincs <dir>`).
5. `sigcom_fibroblast_signatures.csv`: SigCom LINCS metadata for IMR90, WI38, HFL1, 1HAE, MCH58.
6. WebFetch (model summaries, not stored verbatim): Ginkgo GDPx2 portal page; GDPx2 bioRxiv full text; Nature Communications s41467-025-66597-9. The Cancer Cell ADAM12 page returned 403. The Zenodo 14794723 API returned non-JSON.

## Open items (not resolved in this task)

- GSE261783 achieved cells per sgRNA (bar C). Cannot change the reading, because E fails.
- The ADAM12 human Perturb-seq accession.
- Whether targeted probe panels count as transcriptomic, and whether regulatory elements count as D1 agents. Both moot here, because A fails for those datasets.
- LINCS 2020 `siginfo_beta.txt` (465 MB) was not opened.
