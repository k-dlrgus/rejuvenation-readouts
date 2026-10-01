"""Stage 1 — retrieve public sources for the GSE297234 / Lu et al. MD claim. No fitting.

Public sources only: GEO, PubMed, PMC, preprint servers, code repositories, open supplements.
Fail loudly. Do not paraphrase from memory. Do not circumvent a paywall.
"""
from __future__ import annotations

import json
import re
import sys
import time
import zipfile
from pathlib import Path
from urllib.parse import urlparse
from xml.etree import ElementTree as ET

import pandas as pd
import requests

sys.path.insert(0, str(Path(__file__).resolve().parent))
from geo_meta import fetch_soft, parse_series  # noqa: E402
from md_common import (  # noqa: E402
    MD_DIR, MD_RAW, RETRIEVAL_DATE, PMID, DOI, PII, PII_SD, GSE_SC, GSE_BULK,
    TITLE, GATE_PATH, SOURCES_CSV, SOURCES_JSON, DISALLOWED_HOST_SUBSTR,
    PREREG_STAGE2, write_prereg_stage2,
    StopStep, Logger, dump_json, load_json, jsonable, load_manifest, save_manifest,
    md_log_banner, progress_snapshot, record_failure,
)

EUTILS = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils"
UA = (
    "AgeIdentitySeparabilityBenchmark/0.1 "
    "(Prompt H public-metadata retrieval; +https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE297234)"
)
HEADERS = {"User-Agent": UA, "Accept": "*/*"}
NCBI_SLEEP = 0.4
MAX_SAVE = 8_000_000
MAX_EXCERPT = 4000


def _disallowed(url: str) -> bool:
    host = (urlparse(url).hostname or "").lower()
    return any(x in host or x in url.lower() for x in DISALLOWED_HOST_SUBSTR)


def _safe_name(url: str, prefix: str) -> str:
    path = urlparse(url).path.rstrip("/") or "index"
    leaf = path.split("/")[-1] or "index"
    leaf = re.sub(r"[^A-Za-z0-9._-]+", "_", leaf)[:120]
    host = re.sub(r"[^A-Za-z0-9._-]+", "_", (urlparse(url).hostname or "host"))[:40]
    q = urlparse(url).query
    if q:
        leaf = leaf + "_" + re.sub(r"[^A-Za-z0-9._-]+", "_", q)[:80]
    name = f"{prefix}__{host}__{leaf}"
    return re.sub(r"[^A-Za-z0-9._-]+", "_", name)[:180]


def http_get(url, dest_name, log, timeout=90, max_save=MAX_SAVE, allow_binary=True):
    """GET a URL. Record status. Save body when it is a public document, not a paywall interstitial."""
    rec = dict(
        url=url,
        retrieval_date=RETRIEVAL_DATE,
        status=None,
        ok=False,
        content_type=None,
        nbytes=0,
        path=None,
        excerpt="",
        error=None,
        note="",
    )
    if _disallowed(url):
        rec["error"] = "disallowed host (paywall circumvention / third-party full-text mirror)"
        rec["note"] = "not retrieved"
        log(f"[skip-disallowed] {url}")
        return rec
    dest_name = re.sub(r"[^A-Za-z0-9._-]+", "_", str(dest_name))[:180]
    dest = MD_RAW / dest_name
    try:
        r = requests.get(url, headers=HEADERS, timeout=timeout, allow_redirects=True)
        rec["status"] = int(r.status_code)
        rec["content_type"] = r.headers.get("Content-Type", "")
        rec["final_url"] = r.url
        rec["nbytes"] = len(r.content or b"")
        rec["ok"] = r.status_code == 200
        ctype = (rec["content_type"] or "").lower()
        body = r.content or b""
        is_text = any(x in ctype for x in (
            "text/", "json", "xml", "javascript", "csv", "html",
        )) or dest_name.endswith((".txt", ".json", ".xml", ".html", ".csv"))
        magic_pdf = body[:5] == b"%PDF-"
        magic_zip = body[:2] == b"PK"
        too_big = rec["nbytes"] > max_save
        if rec["ok"] and (is_text or magic_pdf or magic_zip or allow_binary) and not too_big:
            try:
                dest.write_bytes(body)
                rec["path"] = str(dest)
            except OSError as e:
                rec["error"] = f"write {type(e).__name__}: {e}"
                rec["ok"] = False
        elif rec["ok"] and too_big:
            try:
                dest.write_bytes(body[:max_save])
                rec["path"] = str(dest)
                rec["note"] = f"truncated to {max_save} bytes"
            except OSError as e:
                rec["error"] = f"write {type(e).__name__}: {e}"
                rec["ok"] = False
        text = ""
        if is_text or (not magic_pdf and not magic_zip):
            try:
                text = body.decode("utf-8", errors="replace")
            except Exception:
                text = ""
        rec["excerpt"] = text[:MAX_EXCERPT].replace("\r", "")
        if magic_pdf:
            rec["note"] = (rec.get("note") or "") + " PDF magic"
            rec["excerpt"] = _pdf_strings(body)[:MAX_EXCERPT]
        if magic_zip:
            rec["note"] = (rec.get("note") or "") + " ZIP/XLSX magic"
            if rec.get("path"):
                rec["excerpt"] = _xlsx_preview(Path(rec["path"]))[:MAX_EXCERPT]
        log(f"[get] {r.status_code} {rec['nbytes']}B {ctype[:40]} {url}")
    except requests.RequestException as e:
        rec["error"] = f"{type(e).__name__}: {e}"
        rec["ok"] = False
        log(f"[get-fail] {url} {rec['error']}")
    return rec


def _pdf_strings(body: bytes, min_len=6) -> str:
    """Crude printable-string scrape. Not a PDF parser. Used only to see if a gene list is in an open PDF."""
    out = []
    buf = []
    for b in body:
        if 32 <= b < 127:
            buf.append(chr(b))
        else:
            if len(buf) >= min_len:
                out.append("".join(buf))
            buf = []
    if len(buf) >= min_len:
        out.append("".join(buf))
    text = "\n".join(out)
    # Drop obvious PDF operator noise.
    keep = []
    for line in text.splitlines():
        if re.fullmatch(r"[\s0-9.\[\]/]+", line):
            continue
        keep.append(line)
    return "\n".join(keep)[: 50_000]


def _xlsx_preview(path: Path) -> str:
    try:
        with zipfile.ZipFile(path) as z:
            names = z.namelist()
            bits = ["xlsx members: " + ", ".join(names[:40])]
            if "xl/sharedStrings.xml" in names:
                raw = z.read("xl/sharedStrings.xml")
                # Strip tags crudely.
                text = re.sub(rb"<[^>]+>", b" ", raw).decode("utf-8", errors="replace")
                text = re.sub(r"\s+", " ", text)
                bits.append("sharedStrings: " + text[:3000])
            sheets = [n for n in names if n.startswith("xl/worksheets/") and n.endswith(".xml")]
            for s in sheets[:4]:
                raw = z.read(s)
                text = re.sub(rb"<[^>]+>", b" ", raw).decode("utf-8", errors="replace")
                text = re.sub(r"\s+", " ", text)
                bits.append(f"{s}: " + text[:1500])
            return "\n".join(bits)
    except Exception as e:
        return f"xlsx preview failed: {type(e).__name__}: {e}"


def _parse_series_full(text: str) -> dict:
    return parse_series(text)


def _parse_gsm(text: str) -> list[dict]:
    recs, cur = [], {}
    for line in text.splitlines():
        if line.startswith("^SAMPLE"):
            if cur:
                recs.append(cur)
            cur = dict(gsm=line.split("=")[-1].strip())
            continue
        if not line.startswith("!Sample_"):
            continue
        k, _, v = line[1:].partition(" = ")
        cur.setdefault(k.strip(), []).append(v.strip())
    if cur:
        recs.append(cur)
    return recs


def record_source(rows, source, url, found, rec=None, extra=None):
    r = dict(
        source=source,
        url=url,
        retrieval_date=RETRIEVAL_DATE,
        found=found,
        http_status=(rec or {}).get("status"),
        nbytes=(rec or {}).get("nbytes"),
        path=(rec or {}).get("path"),
        error=(rec or {}).get("error"),
        note=(rec or {}).get("note") or "",
    )
    if extra:
        r.update(extra)
    rows.append(r)
    return r


def eutils_get(path, params, log, dest_name):
    params = dict(params)
    params.setdefault("tool", "age_identity_separability")
    params.setdefault("email", "none@example.invalid")
    url = f"{EUTILS}/{path}?" + "&".join(f"{k}={requests.utils.quote(str(v))}" for k, v in params.items())
    # requests.get with params is cleaner
    rec = dict(url=f"{EUTILS}/{path}", retrieval_date=RETRIEVAL_DATE, status=None, ok=False,
               content_type=None, nbytes=0, path=None, excerpt="", error=None, note="eutils")
    try:
        time.sleep(NCBI_SLEEP)
        r = requests.get(f"{EUTILS}/{path}", params=params, headers=HEADERS, timeout=90)
        rec["status"] = int(r.status_code)
        rec["content_type"] = r.headers.get("Content-Type", "")
        rec["final_url"] = r.url
        rec["nbytes"] = len(r.content or b"")
        rec["ok"] = r.status_code == 200
        rec["url"] = r.url
        dest_name = re.sub(r"[^A-Za-z0-9._-]+", "_", str(dest_name))[:180]
        dest = MD_RAW / dest_name
        dest.write_bytes(r.content or b"")
        rec["path"] = str(dest)
        rec["excerpt"] = (r.text or "")[:MAX_EXCERPT]
        log(f"[eutils] {r.status_code} {path} {rec['nbytes']}B")
    except requests.RequestException as e:
        rec["error"] = f"{type(e).__name__}: {e}"
        log(f"[eutils-fail] {path} {rec['error']}")
    return rec


# ---------------------------------------------------------------------------
# GEO
# ---------------------------------------------------------------------------

def fetch_geo(log, rows, searches):
    out = {}
    for acc in (GSE_SC, GSE_BULK):
        log(f"[geo] SOFT series {acc}")
        try:
            text = fetch_soft(acc, targ_gsm=False)
        except Exception as e:
            record_source(rows, f"GEO {acc} SOFT series",
                          f"https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc={acc}&targ=self&form=text&view=brief",
                          f"NOT RETRIEVED: {type(e).__name__}: {e}")
            continue
        p = MD_RAW / f"geo_{acc}_series.soft.txt"
        p.write_text(text, encoding="utf-8", errors="replace")
        series = _parse_series_full(text)
        dump_json(MD_DIR / f"geo_{acc}_series.json", jsonable(series))
        found_bits = []
        for k in (
            "Series_title", "Series_summary", "Series_overall_design",
            "Series_pubmed_id", "Series_type", "Series_platform_id",
            "Series_sample_organism", "Series_supplementary_file",
            "Series_sample_id", "Series_relation", "Series_contributor",
            "Series_web_link", "Series_contact_institute", "Series_status",
            "Series_submission_date", "Series_last_update_date",
        ):
            if k in series:
                found_bits.append(f"{k}={series[k] if k != 'Series_summary' else ['<see json>']}")
        record_source(
            rows, f"GEO {acc} SOFT series",
            f"https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc={acc}&targ=self&form=text&view=brief",
            "retrieved. keys=" + ",".join(sorted(series.keys())),
            extra=dict(path=str(p), n_samples=len(series.get("Series_sample_id", []))),
        )
        out[acc] = series
        rec_full = http_get(
            f"https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc={acc}&targ=self&form=text&view=full",
            f"geo_{acc}_series_full.soft.txt", log, timeout=120,
        )
        if rec_full.get("ok") and rec_full.get("path"):
            series_full = _parse_series_full(Path(rec_full["path"]).read_text(encoding="utf-8", errors="replace"))
            dump_json(MD_DIR / f"geo_{acc}_series_full.json", jsonable(series_full))
            for k in ("Series_relation", "Series_supplementary_file", "Series_web_link",
                      "Series_contributor", "Series_contact_institute"):
                if k in series_full and k not in series:
                    series[k] = series_full[k]
            out[acc] = series
            record_source(
                rows, f"GEO {acc} SOFT series full", rec_full["url"],
                "retrieved. keys=" + ",".join(sorted(series_full.keys()))
                + f" relation={series_full.get('Series_relation')} "
                + f"web_link={series_full.get('Series_web_link')}",
                rec_full,
            )
        else:
            record_source(rows, f"GEO {acc} SOFT series full", rec_full["url"],
                          f"NOT RETRIEVED status={rec_full.get('status')} error={rec_full.get('error')}", rec_full)
        log(f"[geo] fetching GSM SOFT {acc}")
        try:
            gsm_text = fetch_soft(acc, targ_gsm=True)
        except Exception as e:
            record_source(rows, f"GEO {acc} SOFT GSM",
                          f"https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc={acc}&targ=gsm&form=text&view=brief",
                          f"NOT RETRIEVED: {type(e).__name__}: {e}")
            continue
        gp = MD_RAW / f"geo_{acc}_gsm.soft.txt"
        gp.write_text(gsm_text, encoding="utf-8", errors="replace")
        samples = _parse_gsm(gsm_text)
        dump_json(MD_DIR / f"geo_{acc}_gsm.json", jsonable(samples))
        char_keys = set()
        suppl = []
        for s in samples:
            for c in s.get("Sample_characteristics_ch1", []):
                if ":" in c:
                    char_keys.add(c.split(":", 1)[0].strip().lower())
            for k, v in s.items():
                if "supplementary" in k.lower():
                    suppl.extend(v if isinstance(v, list) else [v])
        record_source(
            rows, f"GEO {acc} SOFT GSM",
            f"https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc={acc}&targ=gsm&form=text&view=brief",
            f"retrieved n_gsm={len(samples)} characteristic_keys={sorted(char_keys)} "
            f"sample_supplementary_n={len(suppl)}",
            extra=dict(path=str(gp)),
        )
        # Full HTML landing (may timeout; record either way).
        rec = http_get(
            f"https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc={acc}",
            f"geo_{acc}_html.txt", log, timeout=120,
        )
        record_source(rows, f"GEO {acc} HTML", rec["url"],
                      "retrieved" if rec["ok"] else f"NOT RETRIEVED status={rec.get('status')} error={rec.get('error')}",
                      rec)

    # SuperSeries / siblings via PubMed-linked GEO
    searches.append(dict(
        query=f"{PMID}[PMID] in db=gds", via="eutils esearch gds",
        result="queued",
    ))
    rec = eutils_get("esearch.fcgi", dict(db="gds", term=f"{PMID}[PMID]", retmode="json", retmax=50),
                     log, "eutils_gds_pmid.json")
    record_source(rows, "NCBI GDS esearch by PMID", rec["url"],
                  "retrieved" if rec["ok"] else f"NOT RETRIEVED {rec.get('error')}", rec)
    gds_ids = []
    if rec.get("path"):
        try:
            js = json.loads(Path(rec["path"]).read_text(encoding="utf-8"))
            gds_ids = js.get("esearchresult", {}).get("idlist", [])
            searches[-1]["result"] = f"n={js.get('esearchresult', {}).get('count')} ids={gds_ids}"
        except Exception as e:
            searches[-1]["result"] = f"parse fail {e}"
    if gds_ids:
        rec2 = eutils_get("esummary.fcgi",
                          dict(db="gds", id=",".join(gds_ids), retmode="json"),
                          log, "eutils_gds_pmid_summary.json")
        accs = []
        if rec2.get("path"):
            try:
                js = json.loads(Path(rec2["path"]).read_text(encoding="utf-8"))
                res = js.get("result", {})
                for uid in res.get("uids", []):
                    accs.append(dict(
                        uid=uid,
                        accession=res.get(uid, {}).get("accession"),
                        title=res.get(uid, {}).get("title"),
                        n_samples=res.get(uid, {}).get("n_samples"),
                        gdstype=res.get(uid, {}).get("gdstype"),
                        summary=(res.get(uid, {}).get("summary") or "")[:500],
                    ))
            except Exception as e:
                accs = [dict(error=str(e))]
        dump_json(MD_DIR / "geo_gds_by_pmid.json", jsonable(accs))
        record_source(rows, "NCBI GDS esummary by PMID", rec2["url"],
                      f"retrieved accessions={[a.get('accession') for a in accs]}", rec2)
        # Any SuperSeries not already fetched
        for a in accs:
            acc = a.get("accession") or ""
            if acc.startswith("GSE") and acc not in (GSE_SC, GSE_BULK) and acc not in out:
                log(f"[geo] additional series from PMID link: {acc}")
                try:
                    text = fetch_soft(acc, targ_gsm=False)
                    p = MD_RAW / f"geo_{acc}_series.soft.txt"
                    p.write_text(text, encoding="utf-8", errors="replace")
                    series = _parse_series_full(text)
                    dump_json(MD_DIR / f"geo_{acc}_series.json", jsonable(series))
                    rel = series.get("Series_relation", [])
                    record_source(
                        rows, f"GEO {acc} SOFT series (PMID-linked sibling/super)",
                        f"https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc={acc}&targ=self&form=text&view=brief",
                        f"retrieved title={series.get('Series_title')} relation={rel} "
                        f"suppl={series.get('Series_supplementary_file')}",
                        extra=dict(path=str(p)),
                    )
                    out[acc] = series
                except Exception as e:
                    record_source(
                        rows, f"GEO {acc} SOFT series (PMID-linked)",
                        f"https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc={acc}",
                        f"NOT RETRIEVED: {type(e).__name__}: {e}",
                    )

    # BioProject / SRA from GSE297234 relation fields
    series = out.get(GSE_SC) or {}
    rels = series.get("Series_relation") or []
    dump_json(MD_DIR / "geo_GSE297234_relation.json", jsonable(rels))
    for rel in rels:
        m = re.search(r"(PRJ[A-Z]{2}\d+|SRP\d+)", rel, re.I)
        if m:
            acc_id = m.group(1)
            rec = http_get(
                f"https://www.ncbi.nlm.nih.gov/bioproject/?term={acc_id}",
                f"bioproject_{acc_id}.html", log, timeout=60,
            )
            record_source(rows, f"BioProject lookup {acc_id}", rec["url"],
                          "retrieved" if rec["ok"] else f"NOT RETRIEVED {rec.get('error')}", rec)
    rec = eutils_get("esearch.fcgi", dict(db="bioproject", term="GSE297234", retmode="json", retmax=20),
                     log, "eutils_bioproject_gse297234.json")
    record_source(rows, "NCBI BioProject esearch GSE297234", rec["url"],
                  "retrieved" if rec["ok"] else f"NOT RETRIEVED {rec.get('error')}", rec)
    rec = eutils_get("esearch.fcgi", dict(db="sra", term="GSE297234", retmode="json", retmax=20),
                     log, "eutils_sra_gse297234.json")
    record_source(rows, "NCBI SRA esearch GSE297234", rec["url"],
                  "retrieved" if rec["ok"] else f"NOT RETRIEVED {rec.get('error')}", rec)
    return out


# ---------------------------------------------------------------------------
# PubMed / journal
# ---------------------------------------------------------------------------

def fetch_pubmed_journal(log, rows, searches):
    rec = eutils_get(
        "efetch.fcgi",
        dict(db="pubmed", id=PMID, rettype="xml", retmode="xml"),
        log, "pubmed_40816266.xml",
    )
    record_source(rows, "PubMed XML", rec["url"],
                  "retrieved" if rec["ok"] else f"NOT RETRIEVED {rec.get('error')}", rec)
    pubmed_parsed = {}
    if rec.get("path") and rec.get("ok"):
        try:
            root = ET.parse(rec["path"]).getroot()
            art = root.find(".//PubmedArticle")
            pubmed_parsed["pmid"] = PMID
            pubmed_parsed["doi"] = None
            for aid in root.findall(".//ArticleId"):
                if aid.get("IdType") == "doi":
                    pubmed_parsed["doi"] = (aid.text or "").strip()
            pubmed_parsed["title"] = "".join(root.findtext(".//ArticleTitle") or "")
            abs_bits = ["".join(a.itertext()) for a in root.findall(".//AbstractText")]
            pubmed_parsed["abstract"] = "\n".join(abs_bits)
            pubmed_parsed["journal"] = root.findtext(".//Journal/Title")
            pubmed_parsed["volume"] = root.findtext(".//JournalIssue/Volume")
            pubmed_parsed["issue"] = root.findtext(".//JournalIssue/Issue")
            pubmed_parsed["pages"] = root.findtext(".//Pagination/MedlinePgn")
            pubmed_parsed["year"] = root.findtext(".//PubDate/Year")
            pubmed_parsed["authors"] = [
                f"{(a.findtext('LastName') or '')} {(a.findtext('ForeName') or '')}".strip()
                for a in root.findall(".//Author")
            ]
            pubmed_parsed["pmc"] = None
            for aid in root.findall(".//ArticleId"):
                if aid.get("IdType") == "pmc":
                    pubmed_parsed["pmc"] = (aid.text or "").strip()
            pubmed_parsed["status_pubstatus"] = root.findtext(".//PublicationStatus")
            pubmed_parsed["grants"] = [
                dict(id=g.findtext("GrantID"), agency=g.findtext("Agency"), country=g.findtext("Country"))
                for g in root.findall(".//Grant")
            ]
        except Exception as e:
            pubmed_parsed["parse_error"] = f"{type(e).__name__}: {e}"
        dump_json(MD_DIR / "pubmed_parsed.json", jsonable(pubmed_parsed))
        record_source(
            rows, "PubMed parsed fields",
            f"https://pubmed.ncbi.nlm.nih.gov/{PMID}/",
            "abstract=" + ("yes" if pubmed_parsed.get("abstract") else "NO")
            + f" doi={pubmed_parsed.get('doi')} pmc={pubmed_parsed.get('pmc')!r} "
            + f"n_authors={len(pubmed_parsed.get('authors') or [])} "
            + f"journal={pubmed_parsed.get('journal')} "
            + f"{pubmed_parsed.get('year')};{pubmed_parsed.get('volume')}({pubmed_parsed.get('issue')}):{pubmed_parsed.get('pages')}",
        )

    rec = http_get(f"https://pubmed.ncbi.nlm.nih.gov/{PMID}/", "pubmed_html.html", log)
    record_source(rows, "PubMed HTML", rec["url"],
                  "retrieved" if rec["ok"] else f"NOT RETRIEVED {rec.get('error')}", rec)

    rec = http_get(f"https://www.cell.com/cell/fulltext/{PII}", "cell_landing.html", log, timeout=120)
    paywall = False
    if rec.get("excerpt"):
        paywall = bool(re.search(r"Get full text access|Log in, subscribe or purchase", rec["excerpt"], re.I))
    record_source(
        rows, "Cell journal landing page", rec["url"],
        ("retrieved; landing page has highlights + summary; "
         + ("paywall interstitial present" if paywall else "paywall interstitial not seen in excerpt")
         + "; STAR Methods listed as a site nav item, not article methods text"),
        rec,
        extra=dict(paywall_interstitial=paywall),
    )
    rec = http_get(f"https://doi.org/{DOI}", "doi_resolve.html", log, timeout=60)
    record_source(rows, "DOI resolver", rec["url"],
                  f"retrieved status={rec.get('status')} final={rec.get('final_url')}"
                  if rec.get("status") else f"NOT RETRIEVED {rec.get('error')}", rec)
    rec = http_get(
        f"https://www.sciencedirect.com/science/article/pii/{PII_SD}",
        "sciencedirect_landing.html", log, timeout=120,
    )
    record_source(rows, "ScienceDirect landing", rec["url"],
                  "retrieved" if rec["ok"] else f"NOT RETRIEVED {rec.get('error')}", rec)
    return pubmed_parsed


# ---------------------------------------------------------------------------
# PMC / Europe PMC / Unpaywall / OpenAlex / Crossref
# ---------------------------------------------------------------------------

def fetch_pmc_oa(log, rows, searches, pubmed_parsed):
    searches.append(dict(query=f"{PMID}[pmid] in db=pmc", via="eutils esearch pmc", result="queued"))
    rec = eutils_get("esearch.fcgi", dict(db="pmc", term=f"{PMID}[pmid]", retmode="json", retmax=20),
                     log, "eutils_pmc_pmid.json")
    pmc_ids = []
    if rec.get("path"):
        try:
            js = json.loads(Path(rec["path"]).read_text(encoding="utf-8"))
            pmc_ids = js.get("esearchresult", {}).get("idlist", [])
            searches[-1]["result"] = f"n={js.get('esearchresult', {}).get('count')} ids={pmc_ids}"
        except Exception as e:
            searches[-1]["result"] = f"parse fail {e}"
    record_source(
        rows, "PMC esearch by PMID", rec["url"],
        f"retrieved pmc_ids={pmc_ids} (empty list = no PMC deposit found by this query)"
        if rec["ok"] else f"NOT RETRIEVED {rec.get('error')}",
        rec,
    )
    rec = eutils_get("elink.fcgi",
                     dict(dbfrom="pubmed", db="pmc", id=PMID, retmode="json"),
                     log, "eutils_elink_pubmed_pmc.json")
    record_source(rows, "elink PubMed→PMC", rec["url"],
                  "retrieved" if rec["ok"] else f"NOT RETRIEVED {rec.get('error')}", rec)

    rec = http_get(
        f"https://www.ebi.ac.uk/europepmc/webservices/rest/search?query=EXT_ID:{PMID}&resultType=core&format=json",
        "europepmc_core.json", log,
    )
    epmc = {}
    if rec.get("path") and rec.get("ok"):
        try:
            epmc = json.loads(Path(rec["path"]).read_text(encoding="utf-8"))
        except Exception as e:
            epmc = dict(parse_error=str(e))
        dump_json(MD_DIR / "europepmc_core_parsed.json", jsonable(_trim_epmc(epmc)))
    rec_hit = ((epmc.get("resultList") or {}).get("result") or [{}])
    hit = rec_hit[0] if rec_hit else {}
    record_source(
        rows, "Europe PMC core JSON", rec["url"],
        "retrieved "
        + f"isOpenAccess={hit.get('isOpenAccess')!r} pmcid={hit.get('pmcid')!r} "
        + f"hasTextMinedTerms={hit.get('hasTextMinedTerms')!r} "
        + f"fullTextUrlList_n={len(((hit.get('fullTextUrlList') or {}).get('fullTextUrl') or []))} "
        + f"license={hit.get('license')!r}",
        rec,
    )
    # Follow OA full-text URLs only (not publisher subscription PDFs that are not marked OA).
    for u in ((hit.get("fullTextUrlList") or {}).get("fullTextUrl") or []):
        avail = (u.get("availability") or "").lower()
        docu = u.get("documentStyle")
        url = u.get("url")
        if not url:
            continue
        # Europe PMC sometimes lists the publisher PDF as 'Subscription'. Skip those.
        if avail in ("subscription", "restricted"):
            record_source(rows, f"Europe PMC fullTextUrl ({avail} {docu})", url,
                          f"listed but not fetched (availability={avail})")
            continue
        recu = http_get(url, _safe_name(url, "epmc_ft"), log, timeout=120)
        record_source(rows, f"Europe PMC fullTextUrl ({avail} {docu})", url,
                      "retrieved" if recu["ok"] else f"NOT RETRIEVED {recu.get('error')}", recu)

    rec = http_get(
        f"https://www.ebi.ac.uk/europepmc/webservices/rest/search?query=HAS_PREPRINT:Y%20AND%20EXT_ID:{PMID}&format=json",
        "europepmc_has_preprint.json", log,
    )
    record_source(rows, "Europe PMC HAS_PREPRINT AND this PMID", rec["url"],
                  "retrieved" if rec["ok"] else f"NOT RETRIEVED {rec.get('error')}", rec)

    rec = http_get(
        f"https://api.unpaywall.org/v2/{DOI}?email=none@example.invalid",
        "unpaywall.json", log,
    )
    oa = {}
    if rec.get("path") and rec.get("ok"):
        try:
            oa = json.loads(Path(rec["path"]).read_text(encoding="utf-8"))
        except Exception as e:
            oa = dict(parse_error=str(e))
        dump_json(MD_DIR / "unpaywall_parsed.json", jsonable(oa))
    record_source(
        rows, "Unpaywall", rec["url"],
        "retrieved "
        + f"oa_status={oa.get('oa_status')!r} is_oa={oa.get('is_oa')!r} "
        + f"has_repository_copy={oa.get('has_repository_copy')!r} "
        + f"best_oa={((oa.get('best_oa_location') or {}) or {}).get('url_for_pdf')!r}"
        if rec["ok"] else f"NOT RETRIEVED {rec.get('error')}",
        rec,
    )
    best = oa.get("best_oa_location") or {}
    if oa.get("is_oa") and best.get("url_for_pdf"):
        recp = http_get(best["url_for_pdf"], "unpaywall_best_oa.pdf", log, timeout=120)
        record_source(rows, "Unpaywall best OA PDF", best["url_for_pdf"],
                      "retrieved" if recp["ok"] else f"NOT RETRIEVED {recp.get('error')}", recp)

    rec = http_get(f"https://api.openalex.org/works/https://doi.org/{DOI}", "openalex.json", log)
    ox = {}
    if rec.get("path") and rec.get("ok"):
        try:
            ox = json.loads(Path(rec["path"]).read_text(encoding="utf-8"))
        except Exception as e:
            ox = dict(parse_error=str(e))
        dump_json(MD_DIR / "openalex_parsed.json", jsonable({
            k: ox.get(k) for k in (
                "id", "doi", "title", "publication_year", "open_access",
                "primary_location", "locations", "ids", "authorships",
            )
        }))
    oa2 = ox.get("open_access") or {}
    record_source(
        rows, "OpenAlex", rec["url"],
        "retrieved "
        + f"is_oa={oa2.get('is_oa')!r} oa_status={oa2.get('oa_status')!r} "
        + f"oa_url={oa2.get('oa_url')!r}"
        if rec["ok"] else f"NOT RETRIEVED {rec.get('error')}",
        rec,
    )

    rec = http_get(f"https://api.crossref.org/works/{DOI}", "crossref.json", log)
    cr = {}
    if rec.get("path") and rec.get("ok"):
        try:
            cr = json.loads(Path(rec["path"]).read_text(encoding="utf-8")).get("message") or {}
        except Exception as e:
            cr = dict(parse_error=str(e))
        dump_json(MD_DIR / "crossref_parsed.json", jsonable({
            "DOI": cr.get("DOI"),
            "title": cr.get("title"),
            "license": cr.get("license"),
            "link": cr.get("link"),
            "relation": cr.get("relation"),
            "URL": cr.get("URL"),
            "abstract": cr.get("abstract"),
            "funder": cr.get("funder"),
            "container-title": cr.get("container-title"),
            "published": cr.get("published"),
        }))
    rel = cr.get("relation") or {}
    record_source(
        rows, "Crossref", rec["url"],
        "retrieved "
        + f"n_license={len(cr.get('license') or [])} "
        + f"relation_keys={list(rel.keys())} "
        + f"has_abstract={bool(cr.get('abstract'))}"
        if rec["ok"] else f"NOT RETRIEVED {rec.get('error')}",
        rec,
    )
    # Crossref preprint relation
    for kind, items in rel.items():
        if "preprint" in kind.lower() or kind in ("is-preprint-of", "has-preprint", "is-version-of"):
            record_source(rows, f"Crossref relation {kind}", rec["url"],
                          f"items={json.dumps(items)[:800]}")
    return dict(epmc=hit, unpaywall=oa, openalex_oa=oa2, crossref_rel=rel, pmc_ids=pmc_ids)


def _trim_epmc(epmc):
    hits = ((epmc.get("resultList") or {}).get("result") or [])
    out = []
    keep = (
        "id", "pmid", "pmcid", "doi", "title", "authorString", "journalTitle",
        "pubYear", "isOpenAccess", "inEPMC", "inPMC", "hasPDF", "hasBook",
        "hasSuppl", "hasData", "hasPreprint", "hasTextMinedTerms", "license",
        "fullTextUrlList", "commentCorrectionList",
    )
    for h in hits[:5]:
        out.append({k: h.get(k) for k in keep})
    return dict(hitCount=epmc.get("hitCount"), result=out)


# ---------------------------------------------------------------------------
# Preprints
# ---------------------------------------------------------------------------

def fetch_preprints(log, rows, searches):
    queries = [
        ("biorxiv_title",
         "https://www.biorxiv.org/search/"
         "prevalent%20mesenchymal%20drift%20in%20aging%20jcode%3Abiorxiv%20numresults%3A25%20sort%3Arelevance-rank"),
        ("biorxiv_author_lu",
         "https://www.biorxiv.org/search/"
         "author%3ALu%20mesenchymal%20drift%20jcode%3Abiorxiv%20numresults%3A25"),
        ("biorxiv_phrase",
         "https://www.biorxiv.org/search/"
         "%22mesenchymal%20drift%22%20jcode%3Abiorxiv%20numresults%3A25%20sort%3Arelevance-rank"),
        ("medrxiv_phrase",
         "https://www.medrxiv.org/search/"
         "%22mesenchymal%20drift%22%20jcode%3Amedrxiv%20numresults%3A25"),
        ("researchsquare_title",
         "https://www.researchsquare.com/browse?q=Prevalent%20mesenchymal%20drift%20in%20aging"),
        ("osf_preprints",
         "https://osf.io/preprints/discover?q=Prevalent%20mesenchymal%20drift"),
        ("ssrn",
         "https://www.ssrn.com/index.cfm/en/search/?q=Prevalent%20mesenchymal%20drift"),
    ]
    for name, url in queries:
        searches.append(dict(query=name, via=url, result="queued"))
        rec = http_get(url, f"preprint_{name}.html", log, timeout=90)
        n_hits = None
        excerpt = rec.get("excerpt") or ""
        m = re.search(r"(\d+)\s+Results", excerpt, re.I)
        if m:
            n_hits = m.group(1)
        # bioRxiv empty result markers
        empty = bool(re.search(r"We were unable to find results|No results|0 Results", excerpt, re.I))
        found = (
            f"retrieved status={rec.get('status')} nbytes={rec.get('nbytes')} "
            f"results_token={n_hits!r} empty_marker={empty}"
        )
        if rec.get("error"):
            found = f"NOT RETRIEVED {rec.get('error')}"
        searches[-1]["result"] = found
        record_source(rows, f"preprint search {name}", url, found, rec)

    # Europe PMC preprint search by title
    rec = http_get(
        "https://www.ebi.ac.uk/europepmc/webservices/rest/search?"
        "query=TITLE:%22Prevalent%20mesenchymal%20drift%22%20AND%20(SRC:PPR%20OR%20PUB_TYPE:preprint)&format=json",
        "europepmc_preprint_title.json", log,
    )
    n = None
    if rec.get("path") and rec.get("ok"):
        try:
            js = json.loads(Path(rec["path"]).read_text(encoding="utf-8"))
            n = js.get("hitCount")
            dump_json(MD_DIR / "europepmc_preprint_title.json", jsonable(js) if n and n <= 20 else dict(hitCount=n))
        except Exception as e:
            n = f"parse fail {e}"
    record_source(rows, "Europe PMC preprint title search", rec["url"],
                  f"retrieved hitCount={n}", rec)
    searches.append(dict(query='TITLE:"Prevalent mesenchymal drift" preprint', via="europepmc", result=f"hitCount={n}"))

    rec = http_get(
        "https://www.ebi.ac.uk/europepmc/webservices/rest/search?"
        "query=%22mesenchymal%20drift%22%20AND%20AUTHOR:Lu%20AND%20(SRC:PPR%20OR%20PUB_TYPE:preprint)&format=json",
        "europepmc_preprint_lu.json", log,
    )
    n = None
    if rec.get("path") and rec.get("ok"):
        try:
            js = json.loads(Path(rec["path"]).read_text(encoding="utf-8"))
            n = js.get("hitCount")
        except Exception as e:
            n = f"parse fail {e}"
    record_source(rows, "Europe PMC preprint author Lu + mesenchymal drift", rec["url"],
                  f"retrieved hitCount={n}", rec)


# ---------------------------------------------------------------------------
# Code / data repositories
# ---------------------------------------------------------------------------

def fetch_code(log, rows, searches):
    gh_q = [
        ("GSE297234", "https://api.github.com/search/repositories?q=GSE297234"),
        ("mesenchymal_drift", "https://api.github.com/search/repositories?q=mesenchymal+drift"),
        ("lu_mesenchymal", "https://api.github.com/search/repositories?q=mesenchymal+drift+Lu"),
        ("altos_md", "https://api.github.com/search/repositories?q=mesenchymal+drift+altos"),
    ]
    for name, url in gh_q:
        time.sleep(1.0)
        searches.append(dict(query=name, via="github repos", result="queued"))
        rec = http_get(url, f"github_repos_{name}.json", log)
        n = None
        items = []
        if rec.get("path") and rec.get("ok"):
            try:
                js = json.loads(Path(rec["path"]).read_text(encoding="utf-8"))
                n = js.get("total_count")
                items = [
                    dict(full_name=i.get("full_name"), html_url=i.get("html_url"),
                         description=i.get("description"), stars=i.get("stargazers_count"))
                    for i in (js.get("items") or [])[:10]
                ]
            except Exception as e:
                n = f"parse fail {e}"
        dump_json(MD_DIR / f"github_repos_{name}.json", jsonable(dict(total_count=n, items=items)))
        searches[-1]["result"] = f"total_count={n} items={[i.get('full_name') for i in items]}"
        record_source(rows, f"GitHub repo search {name}", url,
                      f"retrieved total_count={n} names={[i.get('full_name') for i in items]}"
                      if rec["ok"] else f"NOT RETRIEVED {rec.get('error')} status={rec.get('status')}",
                      rec)

    rec = http_get(
        "https://api.github.com/search/code?q=GSE297234",
        "github_code_GSE297234.json", log,
    )
    record_source(
        rows, "GitHub code search GSE297234", rec["url"],
        "retrieved" if rec["ok"] else
        f"NOT RETRIEVED status={rec.get('status')} (code search often requires auth; not substituting)",
        rec,
    )

    rec = http_get(
        "https://zenodo.org/api/records?q=%22mesenchymal%20drift%22%20AND%20(Lu%20OR%20GSE297234)&size=20",
        "zenodo.json", log,
    )
    n = None
    hits = []
    if rec.get("path") and rec.get("ok"):
        try:
            js = json.loads(Path(rec["path"]).read_text(encoding="utf-8"))
            n = js.get("hits", {}).get("total")
            for h in (js.get("hits", {}).get("hits") or [])[:10]:
                m = h.get("metadata") or {}
                hits.append(dict(id=h.get("id"), title=m.get("title"), doi=m.get("doi"),
                                 url=(h.get("links") or {}).get("html")))
        except Exception as e:
            n = f"parse fail {e}"
    dump_json(MD_DIR / "zenodo_parsed.json", jsonable(dict(total=n, hits=hits)))
    record_source(rows, "Zenodo search", rec["url"],
                  f"retrieved total={n} titles={[h.get('title') for h in hits]}"
                  if rec["ok"] else f"NOT RETRIEVED {rec.get('error')}", rec)

    # figshare
    fig_url = "https://api.figshare.com/v2/articles/search"
    rec = dict(url=fig_url, retrieval_date=RETRIEVAL_DATE, status=None, ok=False,
               content_type=None, nbytes=0, path=None, excerpt="", error=None, note="POST")
    try:
        r = requests.post(
            fig_url,
            json={"search_for": "mesenchymal drift Lu GSE297234"},
            headers={**HEADERS, "Content-Type": "application/json"},
            timeout=60,
        )
        rec["status"] = int(r.status_code)
        rec["ok"] = r.status_code == 200
        rec["nbytes"] = len(r.content or b"")
        rec["content_type"] = r.headers.get("Content-Type", "")
        dest = MD_RAW / "figshare.json"
        dest.write_bytes(r.content or b"")
        rec["path"] = str(dest)
        rec["excerpt"] = (r.text or "")[:MAX_EXCERPT]
        titles = []
        if rec["ok"]:
            try:
                arr = r.json()
                titles = [a.get("title") for a in (arr if isinstance(arr, list) else [])[:10]]
            except Exception:
                titles = []
        record_source(rows, "figshare search POST", fig_url,
                      f"retrieved n={len(titles) if rec['ok'] else 'NA'} titles={titles}", rec)
        log(f"[post] {r.status_code} figshare")
    except requests.RequestException as e:
        rec["error"] = f"{type(e).__name__}: {e}"
        record_source(rows, "figshare search POST", fig_url, f"NOT RETRIEVED {rec['error']}", rec)

    rec = http_get(
        "https://codeocean.com/explore?query=mesenchymal%20drift",
        "codeocean.html", log, timeout=60,
    )
    record_source(rows, "Code Ocean explore", rec["url"],
                  "retrieved" if rec["ok"] else f"NOT RETRIEVED {rec.get('error')}", rec)


# ---------------------------------------------------------------------------
# Open supplements (Elsevier mmc attachments, Cell supplemental page)
# ---------------------------------------------------------------------------

def fetch_open_supplements(log, rows):
    urls = [
        f"https://www.cell.com/cell/supplemental/{PII}",
        f"https://www.cell.com/cms/10.1016/j.cell.2025.07.031/attachment/mmc1.pdf",
        f"https://www.cell.com/cell/fulltext/{PII}#supplementaryMaterial",
        f"https://ars.els-cdn.com/content/image/1-s2.0-{PII_SD}-mmc1.pdf",
        f"https://ars.els-cdn.com/content/image/1-s2.0-{PII_SD}-mmc1.xlsx",
        f"https://ars.els-cdn.com/content/image/1-s2.0-{PII_SD}-mmc2.xlsx",
        f"https://ars.els-cdn.com/content/image/1-s2.0-{PII_SD}-mmc3.xlsx",
        f"https://ars.els-cdn.com/content/image/1-s2.0-{PII_SD}-mmc4.xlsx",
        f"https://ars.els-cdn.com/content/image/1-s2.0-{PII_SD}-mmc5.xlsx",
        f"https://ars.els-cdn.com/content/image/1-s2.0-{PII_SD}-mmc2.pdf",
        f"https://ars.els-cdn.com/content/image/1-s2.0-{PII_SD}-mmc3.pdf",
        f"https://ars.els-cdn.com/content/image/1-s2.0-{PII_SD}-mmc6.xlsx",
        f"https://ars.els-cdn.com/content/image/1-s2.0-{PII_SD}-mmc7.xlsx",
        f"https://www.cell.com/cms/attachment/stateless/pii/{PII}/mmc1.pdf",
    ]
    for url in urls:
        rec = http_get(url, _safe_name(url, "suppl"), log, timeout=60)
        kind = "unknown"
        if rec.get("path"):
            b = Path(rec["path"]).read_bytes()[:8]
            if b.startswith(b"%PDF"):
                kind = "pdf"
            elif b.startswith(b"PK"):
                kind = "xlsx/zip"
            elif b.lower().startswith(b"<!do") or b.lower().startswith(b"<html"):
                kind = "html (likely paywall or 404 page)"
        found = (
            f"status={rec.get('status')} nbytes={rec.get('nbytes')} kind={kind} "
            f"ctype={rec.get('content_type')!r}"
        )
        if rec.get("error"):
            found = f"NOT RETRIEVED {rec.get('error')}"
        # A tiny HTML 404 is not an open supplement.
        if rec.get("ok") and kind.startswith("html"):
            found += " — not an open supplement file"
        record_source(rows, "open-supplement probe", url, found, rec)


# ---------------------------------------------------------------------------
# Gate
# ---------------------------------------------------------------------------

def _text_blob():
    """Concatenate retrieved public text for gate inspection. Skip disallowed hosts."""
    parts = []
    for p in sorted(MD_RAW.glob("*")):
        if p.stat().st_size > 9_000_000:
            continue
        try:
            raw = p.read_bytes()
        except Exception:
            continue
        if raw[:5] == b"%PDF-":
            parts.append(f"\n\n===== {p.name} =====\n" + _pdf_strings(raw)[:80_000])
        elif raw[:2] == b"PK":
            parts.append(f"\n\n===== {p.name} =====\n" + _xlsx_preview(p)[:80_000])
        else:
            parts.append(f"\n\n===== {p.name} =====\n" + raw.decode("utf-8", errors="replace")[:80_000])
    for p in MD_DIR.glob("*.json"):
        try:
            parts.append(f"\n\n===== {p.name} =====\n" + p.read_text(encoding="utf-8", errors="replace")[:80_000])
        except Exception:
            pass
    return "\n".join(parts)


def _xlsx_columns(path: Path) -> dict:
    """Return {sheet_name: {header: [values...]}} from an xlsx using stdlib zip/xml."""
    ns = {"m": "http://schemas.openxmlformats.org/spreadsheetml/2006/main"}

    def col_row(cell_ref):
        m = re.match(r"([A-Z]+)(\d+)", cell_ref or "")
        if not m:
            return None, None
        col, row = m.group(1), int(m.group(2))
        n = 0
        for ch in col:
            n = n * 26 + (ord(ch) - 64)
        return n - 1, row - 1

    out = {}
    with zipfile.ZipFile(path) as z:
        strings = []
        if "xl/sharedStrings.xml" in z.namelist():
            root = ET.fromstring(z.read("xl/sharedStrings.xml"))
            for si in root.findall("m:si", ns):
                strings.append("".join(t.text or "" for t in si.findall(".//m:t", ns)))
        wb = ET.fromstring(z.read("xl/workbook.xml"))
        names = [sh.attrib.get("name") for sh in wb.findall("m:sheets/m:sheet", ns)]
        sheet_files = sorted(
            n for n in z.namelist()
            if n.startswith("xl/worksheets/sheet") and n.endswith(".xml")
        )
        for i, sf in enumerate(sheet_files):
            root = ET.fromstring(z.read(sf))
            grid = {}
            max_r, max_c = 0, 0
            for c in root.findall(".//m:c", ns):
                col, row = col_row(c.attrib.get("r"))
                if col is None:
                    continue
                v = c.find("m:v", ns)
                val = None
                if c.attrib.get("t") == "s" and v is not None and v.text is not None:
                    idx = int(v.text)
                    val = strings[idx] if idx < len(strings) else v.text
                elif v is not None:
                    val = v.text
                if val is not None:
                    grid[(row, col)] = str(val).strip()
                    max_r = max(max_r, row)
                    max_c = max(max_c, col)
            cols = {}
            for c in range(max_c + 1):
                header = grid.get((0, c), f"col{c}")
                vals = [grid[(r, c)] for r in range(1, max_r + 1) if grid.get((r, c))]
                cols[header] = vals
            out[names[i] if i < len(names) else sf] = cols
    return out


def extract_mmc3_md_list(log=None):
    """Parse the open Elsevier mmc3.xlsx. Write gene list. Return dict or None."""
    path = MD_RAW / "suppl__ars.els-cdn.com__1-s2.0-S0092867425008530-mmc3.xlsx"
    url = "https://ars.els-cdn.com/content/image/1-s2.0-S0092867425008530-mmc3.xlsx"
    if not path.exists():
        return None
    sheets = _xlsx_columns(path)
    md_sheet = sheets.get("MD_signatures") or {}
    md_genes = [g for g in (md_sheet.get("MD score") or []) if g]
    tgfb_genes = [g for g in (md_sheet.get("TGFB score") or []) if g]
    rec = dict(
        url=url,
        file=path.name,
        sheets={k: {hk: len(hv) for hk, hv in v.items()} for k, v in sheets.items()},
        n_MD_score=len(md_genes),
        n_TGFB_score=len(tgfb_genes),
        MD_score_genes=md_genes,
        TGFB_score_genes=tgfb_genes,
        computation_in_xlsx=False,
        note=(
            "xlsx contains gene symbols under column headers 'MD score' and 'TGFB score'. "
            "No formula, no AddModuleScore, no ssGSEA, no per-cell vs per-sample instruction in the file."
        ),
    )
    dump_json(MD_DIR / "mmc3_MD_signatures.json", jsonable(rec))
    rows = []
    for g in md_genes:
        rows.append(dict(gene=g, set="MD score", source_url=url, sheet="MD_signatures", column="MD score"))
    for g in tgfb_genes:
        rows.append(dict(gene=g, set="TGFB score", source_url=url, sheet="MD_signatures", column="TGFB score"))
    pd.DataFrame(rows).to_csv(MD_DIR / "mmc3_MD_score_genes.csv", index=False)
    if log:
        log(f"[mmc3] MD score n={len(md_genes)} TGFB n={len(tgfb_genes)} sheets={list(sheets)}")
    return rec


def evaluate_gate(log, rows):
    """Open only if gene list AND computation are recovered from public sources in this session."""
    mmc3 = extract_mmc3_md_list(log)
    gene_list_recovered = bool(mmc3 and mmc3.get("n_MD_score", 0) >= 10)
    gene_list_source = (
        mmc3.get("url") if gene_list_recovered else None
    )
    computation_recovered = False
    computation_source = None
    missing = []

    # Computation: require an explicit method in GEO/PubMed/open-suppl/code, not a citing paper.
    # GEO Sample_data_processing describes Cell Ranger + Seurat sctransform, not MD scoring.
    geo_gsm = (MD_RAW / "geo_GSE297234_gsm.soft.txt").read_text(encoding="utf-8", errors="replace") if (MD_RAW / "geo_GSE297234_gsm.soft.txt").exists() else ""
    methods_ok_files = []
    for p in MD_RAW.glob("*"):
        name = p.name.lower()
        if name.startswith(("geo_", "pubmed_", "suppl_", "epmc_ft", "unpaywall_best", "github_", "zenodo_")):
            methods_ok_files.append(p)
    author_blob = geo_gsm.lower()
    for p in methods_ok_files:
        try:
            raw = p.read_bytes()
        except Exception:
            continue
        if raw[:5] == b"%PDF-":
            author_blob += "\n" + _pdf_strings(raw).lower()
        elif raw[:2] == b"PK":
            # gene lists are not a computation
            continue
        else:
            author_blob += "\n" + raw.decode("utf-8", errors="replace").lower()[:80_000]

    methods_phrases = (
        r"addmodulescore",
        r"aucell",
        r"ssgsea",
        r"\bgsva\b",
        r"module score",
        r"mean z",
        r"z-score of the mean",
        r"geometric mean",
        r"singscore",
        r"\bucell\b",
    )
    for pat in methods_phrases:
        for m in re.finditer(pat, author_blob):
            window = author_blob[max(0, m.start() - 250): m.end() + 250]
            if re.search(r"\bmd\b|mesenchymal drift|md score", window):
                computation_recovered = True
                computation_source = pat
                break
        if computation_recovered:
            break

    if not gene_list_recovered:
        missing.append(
            "MD score gene list (which genes). Public abstract/GEO summary names the phenomenon "
            "but does not list the genes. Would be in a labelled supplementary table."
        )
    if not computation_recovered:
        missing.append(
            "MD score computation (per cell vs per sample; mean / ssGSEA / AddModuleScore / other). "
            "Open supplement mmc3.xlsx lists genes under column 'MD score' but contains no formula. "
            "GEO Sample_data_processing is Cell Ranger v7.1.0 + Seurat v5.0.0 sctransform v2, "
            "and does not mention an MD score. "
            "Would be in STAR Methods of the paywalled article "
            "(https://www.cell.com/cell/fulltext/S0092-8674(25)00853-0 , HTTP 403 this session) "
            "or in an author code repository (GitHub repo search total_count=0)."
        )

    open_gate = bool(gene_list_recovered and computation_recovered)
    if open_gate:
        reason = "MD gene list and computation recovered from public sources"
    elif gene_list_recovered and not computation_recovered:
        reason = (
            f"MD gene list recovered from open supplement mmc3.xlsx "
            f"(n={mmc3.get('n_MD_score')} genes, sheet MD_signatures, column 'MD score'; {mmc3.get('url')}). "
            "Computation of the score from those genes is not in any public source retrieved this session. "
            "STOP after Stage 1."
        )
    else:
        reason = "MD score gene list and computation cannot be recovered from public sources. STOP after Stage 1."

    gate = dict(
        open=open_gate,
        gene_list_recovered=gene_list_recovered,
        computation_recovered=computation_recovered,
        gene_list_source=gene_list_source,
        computation_source=computation_source,
        n_MD_genes=(mmc3 or {}).get("n_MD_score"),
        n_TGFB_genes=(mmc3 or {}).get("n_TGFB_score"),
        mmc3_url=(mmc3 or {}).get("url"),
        missing=missing,
        reason=reason,
        notes=dict(
            citing_methods_are_not_the_authors_methods=True,
            do_not_approximate_addmodulescore=True,
            geo_data_processing_mentions_md_score=bool(re.search(r"md score|mesenchymal drift score", geo_gsm, re.I)),
            mmc3_sheets=(mmc3 or {}).get("sheets"),
        ),
        retrieval_date=RETRIEVAL_DATE,
    )
    dump_json(GATE_PATH, jsonable(gate))
    log(f"[gate] open={open_gate} gene_list={gene_list_recovered} n_MD={gate.get('n_MD_genes')} computation={computation_recovered}")
    log(f"[gate] {reason}")
    return gate


def run_stage1(log=None):
    close_log = False
    if log is None:
        log = Logger(MD_DIR / "s1_report.txt")
        close_log = True
    md_log_banner(log, "STAGE1")
    write_prereg_stage2()
    log("PREREG_STAGE2.flag written before any MD score.")
    log(PREREG_STAGE2[:500] + " …")
    man = load_manifest()
    man["status"] = "STAGE1_RUNNING"
    man["searches_run"] = man.get("searches_run") or []
    save_manifest(man)
    rows = []
    searches = man["searches_run"]
    try:
        geo = fetch_geo(log, rows, searches)
        pubmed_parsed = fetch_pubmed_journal(log, rows, searches)
        oa = fetch_pmc_oa(log, rows, searches, pubmed_parsed)
        fetch_preprints(log, rows, searches)
        fetch_code(log, rows, searches)
        fetch_open_supplements(log, rows)
        gate = evaluate_gate(log, rows)
        df = pd.DataFrame(rows)
        df.to_csv(SOURCES_CSV, index=False)
        dump_json(SOURCES_JSON, jsonable(rows))
        summary = dict(
            n_sources=len(rows),
            n_ok=int(sum(1 for r in rows if str(r.get("found", "")).startswith("retrieved"))),
            n_fail=int(sum(1 for r in rows if "NOT RETRIEVED" in str(r.get("found", "")))),
            geo_accessions=sorted(geo.keys()),
            pubmed_pmc=pubmed_parsed.get("pmc") if isinstance(pubmed_parsed, dict) else None,
            gate_open=gate.get("open"),
            gate_reason=gate.get("reason"),
        )
        dump_json(MD_DIR / "s1_summary.json", jsonable(summary))
        man = load_manifest()
        man["status"] = "STAGE1_DONE_GATE_OPEN" if gate.get("open") else "STAGE1_DONE_GATE_CLOSED"
        man["searches_run"] = searches
        man["s1_summary"] = summary
        save_manifest(man)
        next_a = (
            "Stage 2: implement the recovered MD score exactly and compare to the frozen ruler"
            if gate.get("open") else
            "STOP. Write FINDINGS_MD.md Stage 1 only. Do not run Stage 2."
        )
        progress_snapshot(next_a, stop=man["status"])
        log(f"[stage1] sources={len(rows)} gate_open={gate.get('open')}")
        return gate
    except Exception as e:
        record_failure("stage1", f"{type(e).__name__}: {e}")
        progress_snapshot("diagnose Stage 1 retrieval failure", stop="STOP")
        raise
    finally:
        if close_log and hasattr(log, "close"):
            log.close()


if __name__ == "__main__":
    run_stage1()
