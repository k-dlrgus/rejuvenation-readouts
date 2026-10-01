"""Write FINDINGS_MD.md from results/md/. Does not modify other FINDINGS*.md."""
from __future__ import annotations

import json
import sys
from pathlib import Path

import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parent))
from md_common import (  # noqa: E402
    MD_DIR, FINDINGS_PATH, FROZEN_RULER, RETRIEVAL_DATE, PMID, DOI, GSE_SC, GSE_BULK,
    TITLE, PREREG_STAGE2, PREREG_STAGE2_FLAG, GATE_PATH, SOURCES_CSV,
    MD_SEED, MD_BOOT, N_PERM, N_BOOT, N_RANDOM_DIR,
    STAGE2_VERDICT_FIBRO, FIBRO3_TASK1_SENTENCE, FIBRO3_TASK3_AGED_D0D7,
    load_json, dump_json, jsonable, load_manifest, progress_snapshot,
)
from target_common import md_table  # noqa: E402


def _clip(s, n=500):
    s = "" if s is None else str(s)
    s = s.replace("\n", " ").strip()
    return s if len(s) <= n else s[: n - 1] + "…"


def _load_series(acc):
    p = MD_DIR / f"geo_{acc}_series.json"
    if not p.exists():
        return {}
    return load_json(p)


def _first(d, key, default="not determinable from public sources"):
    v = d.get(key)
    if not v:
        return default
    if isinstance(v, list):
        return v[0] if len(v) == 1 else v
    return v


def _join(d, key, sep=" "):
    v = d.get(key) or []
    if isinstance(v, str):
        return _repair_soft_encoding(v)
    return _repair_soft_encoding(sep.join(str(x) for x in v))


def _repair_soft_encoding(s):
    """GEO SOFT was stored as UTF-8 after latin-1 decode of UTF-8 curly quotes.

    On disk the summary contains U+00E2 U+0080 U+009C (UTF-8 bytes of U+201C
    misread as latin-1). Repair those sequences; leave the rest unchanged.
    """
    if not s or "\xe2" not in s:
        return s
    try:
        return s.encode("latin-1").decode("utf-8")
    except (UnicodeDecodeError, UnicodeEncodeError):
        return s


def _sources_md():
    if not SOURCES_CSV.exists():
        return "_sources.csv missing._"
    df = pd.read_csv(SOURCES_CSV)
    keep = [c for c in ("source", "url", "retrieval_date", "found", "http_status", "error") if c in df.columns]
    df = df[keep].copy()
    df["found"] = df["found"].map(lambda x: _clip(x, 280))
    df["url"] = df["url"].map(lambda x: _clip(x, 180))
    if "http_status" in df.columns:
        def _st(x):
            try:
                if x is None or (isinstance(x, float) and not (x == x)):
                    return "NA"
                return str(int(float(x)))
            except (TypeError, ValueError):
                return "NA"
        df["http_status"] = df["http_status"].map(_st)
    return md_table(df)


def _what_authors_measured(geo_sc, pubmed, gate):
    geo_url = "https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE297234"
    pubmed_url = f"https://pubmed.ncbi.nlm.nih.gov/{PMID}/"
    cell_url = "https://www.cell.com/cell/fulltext/S0092-8674(25)00853-0"
    doi_url = f"https://doi.org/{DOI}"
    epmc_url = f"https://europepmc.org/article/med/{PMID}"
    mmc3_url = "https://ars.els-cdn.com/content/image/1-s2.0-S0092867425008530-mmc3.xlsx"
    mmc3 = load_json(MD_DIR / "mmc3_MD_signatures.json") if (MD_DIR / "mmc3_MD_signatures.json").exists() else {}
    summary = _join(geo_sc, "Series_summary")
    design = _join(geo_sc, "Series_overall_design")
    summary_note = (
        "; SOFT on disk had latin-1-misdecoded UTF-8 quotation marks, repaired to U+201C/U+201D for this quote. "
        "PubMed abstract below is the un-repaired public abstract"
        if "\u201c" in summary or "\u201d" in summary
        else ""
    )
    abstract = (pubmed or {}).get("abstract") or ""
    gsm_path = MD_RAW_GSM = MD_DIR / "raw" / "geo_GSE297234_gsm.soft.txt"
    data_proc = ""
    if gsm_path.exists():
        txt = gsm_path.read_text(encoding="utf-8", errors="replace")
        procs = []
        for line in txt.splitlines():
            if line.startswith("!Sample_data_processing = ") and "Cell Ranger" in line:
                procs.append(line.split(" = ", 1)[1])
                break
        data_proc = procs[0] if procs else ""
    lines = []
    lines.append("Every sentence below is tied to a URL retrieved this session, or is marked not determinable.")
    lines.append("")
    lines.append("### What is the mesenchymal-drift (MD) score?")
    lines.append("")
    n_md = mmc3.get("n_MD_score")
    genes = mmc3.get("MD_score_genes") or []
    if gate.get("gene_list_recovered") and genes:
        lines.append(
            f"- **Genes:** n={n_md} symbols in open supplement mmc3.xlsx sheet `MD_signatures` column `MD score` "
            f"([mmc3.xlsx]({mmc3_url})). "
            "The xlsx does not state weights or a formula. Direction of the gene set is the column name "
            "`MD score`; the abstract describes MD as upregulation of mesenchymal genes "
            f"([PubMed {PMID}]({pubmed_url}); [GEO GSE297234]({geo_url}))."
        )
        lines.append(f"- MD score genes, verbatim from [mmc3.xlsx]({mmc3_url}): " + ", ".join(genes) + ".")
        n_tg = mmc3.get("n_TGFB_score")
        tg = mmc3.get("TGFB_score_genes") or []
        lines.append(
            f"- Same sheet, separate column `TGFB score`: n={n_tg} genes. Not the MD list. "
            + (", ".join(tg) + "." if tg else "")
        )
        lines.append(
            f"- Other sheets in the same file ([mmc3.xlsx]({mmc3_url})): `{mmc3.get('sheets')}`. "
            "Those are labelled gene lists (aging / fibroblast subtype / reprogramming cell state). "
            "They are not a scoring formula."
        )
    else:
        lines.append("- Gene list: **not determinable from public sources.**")
    lines.append("")
    lines.append(
        "- **Computation (how those genes become a number, per cell vs per sample):** "
        "**not determinable from public sources.** "
        f"mmc3.xlsx ([url]({mmc3_url})) is a gene list only. "
        "GEO `Sample_data_processing` (same text on each GSE297234 GSM) is: "
        f"`{data_proc}` "
        f"([GEO GSE297234 SOFT GSM]({geo_url}&targ=gsm&form=text&view=brief)). "
        "That describes Cell Ranger v7.1.0 and Seurat v5.0.0 sctransform v2. It does not mention MD, "
        "AddModuleScore, ssGSEA, GSVA, or a module score. "
        f"Cell landing / STAR Methods: HTTP 403 ([Cell]({cell_url})). "
        "PMC deposit: none "
        f"([PMC esearch](https://www.ncbi.nlm.nih.gov/pmc/?term=40816266); Europe PMC isOpenAccess=N, pmcid=None, "
        f"[Europe PMC]({epmc_url})). "
        "GitHub repositories named for this dataset or 'mesenchymal drift': total_count=0. "
        "Approximating AddModuleScore or Hallmark EMT is disallowed. "
        "The source that would supply the computation is the paywalled STAR Methods of "
        f"Lu et al., *Cell* {doi_url}."
    )
    lines.append("")
    lines.append("### What did they use as the rejuvenation readout?")
    lines.append("")
    lines.append(
        "- **not determinable from public sources** as a named transcriptomic clock vs MD itself vs the "
        "`Aging_signatures` sheet in mmc3.xlsx. "
        "The public abstract says suppression of \"key MD transcription factors leads to epigenetic rejuvenation\" "
        "and that partial reprogramming can \"markedly reduce MD before dedifferentiation and gain of pluripotency, "
        "rejuvenating the aging transcriptome at the cellular and tissue levels\" "
        f"([PubMed {PMID}]({pubmed_url}); [GEO GSE297234]({geo_url})). "
        f"mmc3.xlsx has a sheet `Aging_signatures` with columns `Age up` / `Age down` ([mmc3.xlsx]({mmc3_url})). "
        "Whether that sheet, the MD list, a clock, or something else carried the published rejuvenation claim "
        "is not stated in GEO, PubMed, or the xlsx headers. "
        "Which clock, fit on what, is **not determinable from public sources**."
    )
    lines.append("")
    lines.append("### Was the rejuvenation claim made on individual cells, cell-type-stratified data, or whole-sample averages?")
    lines.append("")
    lines.append(
        "- **not determinable from public sources** for the reported claim. "
        f"The GSE297234 record is 10x Genomics scRNA-seq ([GEO overall design]({geo_url}): "
        f"`{design}`), so the authors had per-cell data. GEO processing built Seurat objects "
        "(see `Sample_data_processing` above). "
        "The abstract claims effects \"at the cellular and tissue levels\" "
        f"([PubMed {PMID}]({pubmed_url})). "
        f"mmc3.xlsx sheet `Reprog_cell_state_signatures` has columns "
        "`Fibroblast`, `PartialReprog`, `EarlyPluripotency`, `Pluripotency`, `NonReprog` "
        f"([mmc3.xlsx]({mmc3_url})), which is a cell-state gene-list table, not a statement that the MD "
        "reduction was reported per state rather than as a whole-sample average. "
        "Sibling bulk series GSE297233 is a different experiment "
        "([GEO GSE297233](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE297233); "
        "overall design: RNA-seq of dox OSK and OCT4YR-SK in GM00731 at 4 days). "
        "GEO GDS esummary `relations` is empty for both series; there is no GEO SuperSeries. "
        "They share PMID 40816266."
    )
    lines.append("")
    lines.append("### Which timepoints did they compare, and did they pre-specify them?")
    lines.append("")
    gsm_json = MD_DIR / f"geo_{GSE_SC}_gsm.json"
    gsm_note = "GSM table not retrieved"
    if gsm_json.exists():
        gsm = load_json(gsm_json)
        titles = []
        for s in gsm:
            t = s.get("Sample_title") or []
            titles.append(t[0] if t else s.get("gsm"))
        gsm_note = "; ".join(str(x) for x in titles)
    lines.append(
        f"- Sample titles retrieved from GEO SOFT GSM ([GSE297234]({geo_url})): `{gsm_note}`. "
        "Days present in those titles: 0, 3, 7, 10. "
        "Overall design says \"treated with Sendai virus OSKM for up to 10 days\" "
        f"([GEO]({geo_url})). Treatment protocol (same on each GSM): CytoTune-iPS 2.0 Sendai, "
        "MOI 5:5:3 for SeV-KLF4-OCT4-SOX2 / SeV-MYC / SeV-KLF4; day 7 replating onto vitronectin; "
        "day 8 Essential 8 Medium. "
        "Whether those four days were pre-specified as the MD contrast, and which pair is the reported "
        "comparison, is **not determinable from public sources**."
    )
    lines.append("")
    lines.append("### What did they claim about identity / dedifferentiation, and how was it measured?")
    lines.append("")
    lines.append(
        "- Claim in the public abstract: Yamanaka-factor partial reprogramming can "
        "\"markedly reduce MD **before** dedifferentiation and gain of pluripotency\" "
        f"([PubMed {PMID}]({pubmed_url}); [GEO GSE297234]({geo_url}); Cell landing HTTP 403). "
        f"mmc3.xlsx sheet `Reprog_cell_state_signatures` lists 200 genes each for "
        "`Fibroblast`, `PartialReprog`, `EarlyPluripotency`, `Pluripotency`, `NonReprog` "
        f"([mmc3.xlsx]({mmc3_url})). "
        "How those lists were scored, and whether they are the published identity/dedifferentiation readout, "
        "is **not determinable from public sources**."
    )
    lines.append("")
    if summary:
        lines.append(f"GEO `Series_summary` (retrieved this session{summary_note}):")
        lines.append("")
        lines.append(f"> {summary}")
        lines.append("")
        lines.append(f"Source: {geo_url}")
        lines.append("")
    if abstract:
        lines.append("PubMed abstract (verbatim, retrieved this session):")
        lines.append("")
        lines.append(f"> {abstract}")
        lines.append("")
        lines.append(f"Source: {pubmed_url}")
        lines.append("")
    return "\n".join(lines)


def write_findings():
    man = load_manifest()
    gate = load_json(GATE_PATH) if GATE_PATH.exists() else {}
    pubmed = load_json(MD_DIR / "pubmed_parsed.json") if (MD_DIR / "pubmed_parsed.json").exists() else {}
    geo_sc = _load_series(GSE_SC)
    geo_bulk = _load_series(GSE_BULK)
    s1 = load_json(MD_DIR / "s1_summary.json") if (MD_DIR / "s1_summary.json").exists() else {}
    if man.get("s1_summary"):
        s1 = {**s1, **man["s1_summary"]}
    if SOURCES_CSV.exists():
        s1["n_sources"] = int(len(pd.read_csv(SOURCES_CSV)))
    dump_json(MD_DIR / "s1_summary.json", jsonable(s1))
    s2 = load_json(MD_DIR / "s2_summary.json") if (MD_DIR / "s2_summary.json").exists() else {}
    stage2_ran = bool(s2.get("ran"))
    if gate.get("open") and stage2_ran:
        status = f"Stage 1 + Stage 2. Gate open. Reading: {s2.get('reading_key', 'not written')}."
    elif gate.get("open") and not stage2_ran:
        status = "Stage 1 complete, gate open, Stage 2 not yet run."
    elif gate.get("gene_list_recovered") and not gate.get("computation_recovered"):
        status = (
            f"Stage 1 only. Gene list recovered (n={gate.get('n_MD_genes')} from open mmc3.xlsx). "
            "Computation not recovered. Gate closed. Stage 2 did not run. No reading fired."
        )
    else:
        status = (
            "Stage 1 only. Gate closed: MD score gene list and computation not recoverable "
            "from public sources. Stage 2 did not run. No reading fired."
        )
    lines = [
        "# FINDINGS_MD — what the GSE297234 authors measured, then (if recoverable) our instrument",
        "",
        f"**Status:** {status} Seed `{MD_SEED}`. boot `{MD_BOOT}`. n_perm={N_PERM}. n_boot={N_BOOT}. n_random={N_RANDOM_DIR}. "
        f"Retrieval date `{RETRIEVAL_DATE}`. Frozen ruler `{FROZEN_RULER}` exists={FROZEN_RULER.exists()}.",
        "",
        "Does not modify any existing `FINDINGS_*.md` or `FALSIFICATION.md`. No GSE325735. "
        "Public sources only. Refit nothing. Nothing averaged across donors, clusters, or regimes.",
        "",
        f"Flag: `results/md/PREREG_STAGE2.flag` exists={PREREG_STAGE2_FLAG.exists()}. "
        f"Gate open={gate.get('open')}.",
        "",
        "## Sources table",
        "",
        f"Retrieved {RETRIEVAL_DATE}. Each row is one HTTP/API call this session.",
        "",
        _sources_md(),
        "",
        "## What the authors measured, as far as public sources show",
        "",
        _what_authors_measured(geo_sc, pubmed, gate),
        "",
        "## Stage 1 gate outcome",
        "",
        f"- open={gate.get('open')}",
        f"- gene_list_recovered={gate.get('gene_list_recovered')}",
        f"- computation_recovered={gate.get('computation_recovered')}",
        f"- gene_list_source={gate.get('gene_list_source')}",
        f"- computation_source={gate.get('computation_source')}",
        f"- reason: {gate.get('reason')}",
        "",
    ]
    if gate.get("missing"):
        lines.append("Missing (verbatim from gate.json):")
        lines.append("")
        for m in gate["missing"]:
            lines.append(f"- {m}")
        lines.append("")
    notes = (gate.get("notes") or {})
    if notes:
        lines.append(
            "Gate notes (not used as a substitute MD score): "
            f"do_not_approximate_AddModuleScore={notes.get('do_not_approximate_addmodulescore')} "
            f"GEO_data_processing_mentions_MD_score={notes.get('geo_data_processing_mentions_md_score')} "
            f"mmc3_sheets={notes.get('mmc3_sheets')}."
        )
        lines.append("")
        lines.append(
            "A later paper using Seurat `AddModuleScore` on an MD set, or the phrase "
            "\"Hallmark EMT plus EMT transcription factors\", is not Lu et al.'s method "
            "and is not implemented."
        )
        lines.append("")
    lines += [
        "Paper identity as retrieved:",
        "",
        f"- Title: {pubmed.get('title') or _first(geo_sc, 'Series_title', TITLE)}",
        f"- PMID: {PMID} — https://pubmed.ncbi.nlm.nih.gov/{PMID}/",
        f"- DOI: {DOI} — https://doi.org/{DOI}",
        f"- Journal: {pubmed.get('journal')} {pubmed.get('year')};"
        f"{pubmed.get('volume')}({pubmed.get('issue')}):{pubmed.get('pages')}",
        f"- PMCID: {pubmed.get('pmc')!r} (None = no PMC id in the PubMed XML this session)",
        f"- Authors (PubMed XML): {pubmed.get('authors')}",
        f"- GEO scRNA-seq: {GSE_SC} — https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc={GSE_SC}",
        f"- GEO bulk sibling: {GSE_BULK} — https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc={GSE_BULK}",
        f"- GEO pubmed_id field: {_first(geo_sc, 'Series_pubmed_id', 'not in SOFT')}",
        f"- GEO supplementary files: {geo_sc.get('Series_supplementary_file')}",
        f"- GEO relation: {geo_sc.get('Series_relation')}",
        f"- GSE297233 title: {_first(geo_bulk, 'Series_title')}",
        f"- GSE297233 design: {_join(geo_bulk, 'Series_overall_design')}",
        "- GEO SuperSeries: none. GDS esummary `relations`=[] for GSE297234 and GSE297233; title search returned exactly those two series (uids 200297234, 200297233).",
        "- BioProject scRNA: PRJNA1263211 (GSE297234). BioProject bulk: PRJNA1263275 (GSE297233).",
        "- SRA runs for PRJNA1263211: count=8 (one per GSM).",
        "- GDS esummary: GSE297234 n_samples=8 suppfile=H5, RDS pdat=2025/08/14 bioproject=PRJNA1263211 relations=[]. GSE297233 n_samples=8 suppfile=CSV pdat=2025/08/14 bioproject=PRJNA1263275 relations=[].",
        "- Unpaywall: is_oa=False, oa_status=closed, has_repository_copy=False. OpenAlex oa_status=closed. Crossref relation={}.",
        "- Preprint: Europe PMC TITLE search hitCount=0; HAS_PREPRINT AND PMID hitCount=0; bioRxiv HTML search HTTP 403; bioRxiv pubs API 'Server not recognized'; Research Square `/browse` page returned the generic latest-preprint listing (473,989), not a title hit.",
        "- Code: GitHub repo search total_count=0 for GSE297234 and for 'mesenchymal drift'. figshare n=0. Zenodo hit is an unrelated record.",
        "- Open Excel supplements: mmc1.xlsx (dataset accessions + cell composition), mmc2.xlsx (plasma-protein tables labelled as reproduction of other PMIDs), mmc3.xlsx (MD / aging / fibroblast / reprogramming gene lists). mmc PDF/docx 404 or Cell HTTP 403.",
        "",
        "## Stage 2 pre-registration (verbatim, written before any MD score)",
        "",
        PREREG_STAGE2,
        "",
        "## Stage 2 tables",
        "",
    ]
    if not gate.get("open"):
        lines.append("Stage 2 did not run. Gate closed.")
        lines.append("")
        if s2:
            lines.append(f"s2_summary.json: `{json.dumps(s2)[:1500]}`")
            lines.append("")
    else:
        for name, path in (
            ("MD all-cell trajectory", MD_DIR / "s2_md_allcell.csv"),
            ("MD cluster trajectory", MD_DIR / "s2_md_cluster.csv"),
            ("Four-cell comparison", MD_DIR / "s2_four_cell.csv"),
            ("MD vs frozen age ρ", MD_DIR / "s2_md_vs_age_rho.csv"),
            ("Composition", MD_DIR / "s2_composition.csv"),
            ("Random-direction null", MD_DIR / "s2_null.csv"),
        ):
            lines.append(f"### {name}")
            lines.append("")
            if path.exists():
                lines.append(md_table(pd.read_csv(path)))
            else:
                lines.append("_file missing._")
            lines.append("")
    lines += [
        "## Reading",
        "",
    ]
    if not gate.get("open"):
        lines.append(
            "No Stage 2 reading fired. The pre-registered instrument/composition/limitation bullets "
            "are not evaluated because the authors' MD score could not be implemented exactly."
        )
        lines.append("")
    else:
        lines.append(s2.get("reading_text") or "Stage 2 reading not written.")
        lines.append("")
    lines += [
        "## What this changes in FINDINGS_FIBRO.md / FINDINGS_FIBRO3.md (quote, do not edit)",
        "",
        "This file does not edit `FINDINGS_FIBRO.md` or `FINDINGS_FIBRO3.md`.",
        "",
        "FINDINGS_FIBRO.md Stage 2 verdict sentence, quoted:",
        "",
        f"> {STAGE2_VERDICT_FIBRO}",
        "",
        "FINDINGS_FIBRO3.md Task 1 reading sentence, quoted:",
        "",
        f"> {FIBRO3_TASK1_SENTENCE}",
        "",
        "FINDINGS_FIBRO3.md Task 3 aged-donor d0→d7 line, quoted:",
        "",
        f"> {FIBRO3_TASK3_AGED_D0D7}",
        "",
    ]
    if not gate.get("open"):
        lines.append(
            "Because the authors' MD computation is not recoverable from public sources, this file does not "
            "adjudicate whether FINDINGS_FIBRO.md / FINDINGS_FIBRO3.md disagree with their *measurement* "
            "or only with their *words*. The gene list is public (mmc3.xlsx, n=205); the scoring formula is not. "
            "That adjudication was the point of Stage 2 and did not run."
        )
        lines.append("")
    lines += [
        "## Limitations",
        "",
        "1. Two donors in GSE297234.",
        "2. Our reproduction of their metric, if Stage 2 had run, would be ours, not theirs.",
        "3. Clusters in FINDINGS_FIBRO3.md are not lineage-tracked.",
        "4. The Cell article is paywalled on the journal landing page retrieved this session; "
        "we did not circumvent that paywall.",
        "5. A third-party PDF of the journal article appeared in web search (rapamycin.news) and was not used.",
        "6. Citing papers (Frontiers review; SENOMORPHIC bioRxiv) that paraphrase MD were not treated as the authors' methods.",
        "7. Frozen ruler trained on GTEx V10 cultured fibroblasts, public AGE 10-year bins. Refit nothing.",
        "8. Seed `20260914`. boot `20260918`. n_perm=200, n_boot=200, n_random=200.",
        "9. GSE325735 not opened.",
        "10. Where a fact could not be retrieved, this file writes \"not determinable from public sources\" rather than an inference.",
        "",
        "## Open questions",
        "",
        "- How is the mmc3.xlsx `MD score` gene list turned into a number (AddModuleScore / ssGSEA / mean / other; per cell vs pseudobulk)? STAR Methods of the paywalled Cell article is the source that would supply it.",
        "- Was the scRNA MD reduction reported per cell state (`Reprog_cell_state_signatures`) or on mixed-population averages?",
        "- Which timepoint pair is the pre-specified MD contrast (d0→d3 vs d0→d7 vs other)?",
        "- Which rejuvenation readout (MD score vs `Aging_signatures` vs a clock vs both) carried the \"aging transcriptome\" claim?",
        "- A legal OA copy (PMC author manuscript or preprint) that includes STAR Methods would reopen Stage 2 without approximation.",
        "",
        "## Files",
        "",
        "| path | content |",
        "|---|---|",
        "| `src/md_common.py`, `md_stage1.py`, `md_stage2.py`, `md_findings.py`, `md_run.py`, `md_finalize_stage1.py` | code |",
        "| `results/md/` | sources, gate, mmc3 gene list, GDS two-series, extra probes, raw retrievals |",
        "| `FINDINGS_MD.md` | this file |",
        "| `PROGRESS_MD.md` | resume state |",
        "",
        f"s1_summary: `{json.dumps(s1)[:1500]}`",
        "",
    ]
    FINDINGS_PATH.write_text("\n".join(lines) + "\n", encoding="utf-8")
    next_action = (
        "STOP. Gate closed: gene list public (mmc3 n=205), computation not in public sources. "
        "Stage 2 does not run until STAR Methods or author code is public. Do not approximate. Do not re-run the searches in PROGRESS_MD.md."
    )
    progress_snapshot(next_action, stop=man.get("status", "FINDINGS_WRITTEN"))
    return FINDINGS_PATH


if __name__ == "__main__":
    print(write_findings())
