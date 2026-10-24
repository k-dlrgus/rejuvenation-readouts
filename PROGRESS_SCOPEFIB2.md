# PROGRESS_SCOPEFIB2

**STOP status:** complete. No STOP key fired. Neither `bar_ambiguous` (S1),
`reproduction_failed` (S2) nor `access_gated` (S3) was triggered as a task STOP — the deposit is
controlled access, but the bars were readable from the public record and the Open Access article,
so S3's "stop on bar reading" applies only to bar C.

**Outcome:** the one `unverified` entry SCOPEFIB left is resolved.

- Human deposit located: **HRA014194**, GSA-Human, BioProject **PRJCA048014**, controlled access.
- Bars: **A fail, B pass, C unverified, D pass, E pass.**
- **No carried reading changes.** GO, PARTIAL and INCOMPLETE all require bar A. SCOPEFIB's
  NO-GO stands verbatim and is **not** superseded:
  `NO-GO: no human dataset passes A+B. Closest near-miss: GSE305412; single missing ingredient: perturbation count (5 distinct knockdown agents applied across TGF-β activation; bar B needs ≥20).`
- The dataset moves from `unverified` to `near_miss`. Exact shortfall: the screen has no
  unstimulated arm, and the unactivated target was profiled as a separate experiment on a
  different 10x chemistry (5′ vs 3′ HT v3.1), which is C3 `split`.

**Single next action:** none in this task. Settling the CMO-to-condition mapping (`mmc1.pdf`,
Figure S1) and deciding whether to apply to the ADAM12 DAC are decisions outside a scoping
session. No access request was submitted.

## Record of what was decided when (times +09:00, 2026-10-24)

| time | event |
|---|---|
| before 10:47 | Read CLAUDE.md, `results/scopefib/PREREG.flag`, `SCOPING_FIBROSIS.md`, `PROGRESS_SCOPEFIB.md`, `src/scopefib_run.py`. Repository files only. No query of any kind. |
| 10:47:11 | `results/scopefib2/PREREG.flag` written: scope, resolution order, reproduction plan, STOP keys S1–S3, the no-change outcome line, and lines 10–87 of the SCOPEFIB flag spliced in byte-for-byte (bars A–E, C1–C4, D1–D7, the four readings, the mandatory columns). Written by script so the carried block cannot drift. Not edited afterwards. |
| 10:48:09 | Pre-registered reproduction check run: 10 checks, 0 mismatches → `results/scopefib2/reproduction.csv`. |
| 10:48:20 | First query of the session (Europe PMC, DOI lookup). |
| 10:48–10:52 | Europe PMC record (not Open Access, no PMC full text). NGDC search endpoints probed; GSA-Human has no own search (404); GSA proper returns only unrelated mouse *Adam12* runs. BIG Search JSON endpoint found. |
| ~10:52 | **HRA014194 found** (`db=hra&q=ADAM12`). Record page, BioProject page, DAC page, individuals and runs listings read; public metadata workbook downloaded. |
| 10:53–10:56 | GEO: GSE278990 and all five SubSeries re-resolved — `Mus musculus` only, confirming SCOPEFIB. Sibling-accession queries (`hra`/`omix`/`gsa` × PRJCA048014, "fibroblast checkpoint", "Jianan Li"): HRA014194 is the only one. |
| ~10:57 | Direct fetches of the article 403 / 400, as SCOPEFIB found. Two WebSearch queries returned nothing usable. |
| 10:57–11:00 | Article read in the built-in browser pane (Open Access). Data-availability statement, screen design, STAR Methods, Perturb-seq analysis, GSVA section and patient section extracted verbatim. Tables S1, S3 and S4 parsed in the page: bar B counts and the 5′/3′ chemistry split. GitHub repo listed (figure scripts only). |
| 11:01:57 | Last query (GSA-Human policy page, for what an access request requires). |
| ~11:07 | `results/scopefib2/_tables.py` → `datasets.csv`, `hra014194_experiments.csv`, `library_counts.csv`, `browser_pane_urls.csv`; `article_excerpts.json` written; `SCOPING_FIBROSIS2.md` written with the flag embedded verbatim (checked by string match, as was the carried SCOPEFIB block). |

## Queries run (do not repeat)

54 requests through the logging fetcher: 42×200, 8×404, 3×403, 1×400. Every URL and status is in
`results/scopefib2/url_records.json`, with the response text stored for the ones that returned it.

1. Europe PMC search by DOI (one hit, PMID 41544628, `isOpenAccess: N`, no PMC id).
2. NGDC: `/search/api/specific?db={hra,omix,gsa,bioproject,biosample,node}` for `ADAM12`, and
   `db={hra,omix,gsa}` for `PRJCA048014`, `fibroblast checkpoint`, `Jianan Li`;
   `/search/databases/ngdc`; `/gsa/search?searchTerm=ADAM12`; `/search/specific?db=hra&q=ADAM12`.
3. GSA-Human: `/browse/HRA014194`, `/browse/dac/HDAC007040`, `/policy`, `/`,
   `/ajaxb/runinstudy?accession=HRA014194&pageSize=2000` (202 run-file rows),
   `/ajaxb/indinstudy?...` (9 individuals), `/browse/experimentDetail/HRX…` for 9 experiments,
   `/file/exportExcelFile?...HRA014194.xlsx`. The endpoints `expinstudy` and `sampleinstudy` do
   not exist (404). BioProject `/bioproject/browse/PRJCA048014`.
4. GEO series text: GSE278990, GSE278974, GSE278975, GSE278980, GSE278983, GSE313761 →
   `results/scopefib2/geo_records.json`.
5. GitHub API: `repos/Qile-Guo/Fibroblast` and its recursive tree (8 blobs, no data files).
6. Failed by design, logged: `www.cell.com/cancer-cell/fulltext/…` 403,
   `sciencedirect.com/science/article/pii/S1535610825005525` 400, `mmc2.xlsx` 403, `mmc4.xlsx`
   403, `gsa-human/search?searchTerm=…` 404 ×4, two NGDC search API guesses 404.
7. Not through the fetcher, so listed separately in `browser_pane_urls.csv`: the article full
   text and `mmc2/mmc4/mmc5.xlsx`, read in the built-in browser pane.
8. Two WebSearch queries for the data-availability statement. Nothing in this task rests on a
   search snippet.

## Open items (not resolved in this task)

- **The CMO-to-condition mapping for the screen's 8 lanes.** Not in the GSA record, the article,
  Tables S1/S3/S4, or the GitHub repository. It is the one thing that would turn the bar A
  reading from an inference into a fact. `mmc1.pdf` (5.20 MB, supplementary figures including
  Figure S1's screen schematic) is the file to open first. Not opened here, because bar A fails
  on C3 `split` even under the reading most favourable to the dataset.
- **Bar C.** Floor of >5 cells per target; no median anywhere public. Needs the per-cell guide
  assignment inside the controlled-access FASTQ data. Not resolvable by any reader without DAC
  approval, and it cannot change a reading while A fails.
- **Table S1A attributes the human ex vivo rows to GEO: GSE278990**, which holds mouse data only.
  Recorded as a record inconsistency; the data-availability statement is the accurate one.
- **Supplementary files not opened**: `mmc1.pdf` (5.20 MB), `mmc6.xlsx` (7.41 MB), `mmc7.xlsx`
  (5.49 MB), `mmc8.zip` (49.73 KB, Data S1), `mmc9.pdf` (36.47 MB). All within the 50 MB limit.
- **SCOPEFIB's own open items are untouched** and remain open: GSE261783's achieved cells per
  sgRNA, whether targeted probe panels count as transcriptomic, whether regulatory elements
  count as D1 agents, and LINCS 2020 `siginfo_beta.txt`.
- **GSA-Human was searched for this study only**, with three study-specific terms. Other
  Chinese-deposited human fibroblast screens would still be missed by SCOPEFIB and by this task.
