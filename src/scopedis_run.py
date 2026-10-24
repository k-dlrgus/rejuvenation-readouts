"""SCOPEDIS: metadata-only record fetcher for the cross-disease in-run-target scoping task.

Scope, bars and the order records are opened in are fixed in results/scopedis/PREREG.flag.
No matrix download, no statistic. Every request is logged with its URL and HTTP status so each
fact in SCOPING_DISEASE.md can be traced to the text it was read from:
  searches (do not count as opened records, D8) -> results/scopedis/searches.json
  opened records and their articles            -> results/scopedis/url_records.json
  GEO series text                              -> results/scopedis/geo_records.json

Usage:
  python src/scopedis_run.py reproduce                       -> reproduction.csv (pre-reg check)
  python src/scopedis_run.py gds <area> <term> [retmax]      -> searches.json (esearch+esummary)
  python src/scopedis_run.py epmc <area> <query> [n]         -> searches.json (Europe PMC search)
  python src/scopedis_run.py series <area> GSE... [...]      -> geo_records.json
  python src/scopedis_run.py url <area> <label> <url>        -> url_records.json (<=50 MB, text)
"""
from __future__ import annotations

import csv
import json
import sys
import time
import urllib.error
import urllib.parse
import urllib.request
from datetime import datetime, timezone
from pathlib import Path

# The Windows console here is cp949; SCOPEFIB lost a query batch to a print encoding error.
for _s in (sys.stdout, sys.stderr):
    try:
        _s.reconfigure(encoding="utf-8", errors="replace")
    except AttributeError:  # pragma: no cover - non-reconfigurable stream
        pass

REPO = Path(__file__).resolve().parents[1]
OUT = REPO / "results" / "scopedis"
UA = "AgeIdentitySeparabilityBenchmark/scopedis (metadata only)"
EUTILS = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/"
MAX_BYTES = 50 * 1024 * 1024


def append(name: str, record: dict) -> None:
    path = OUT / name
    log = json.loads(path.read_text(encoding="utf-8")) if path.exists() else []
    record["fetched_utc"] = datetime.now(timezone.utc).isoformat(timespec="seconds")
    log.append(record)
    path.write_text(json.dumps(log, indent=1, ensure_ascii=False), encoding="utf-8")


def fetch(label: str, url: str, log: str = "url_records.json", area: str = "",
          store: bool = True) -> dict:
    """GET, log URL + HTTP status (including errors), never exceed 50 MB."""
    req = urllib.request.Request(url, headers={"User-Agent": UA, "Accept": "*/*"})
    rec: dict = {"area": area, "label": label, "url": url}
    try:
        with urllib.request.urlopen(req, timeout=120) as r:
            rec["status"] = r.status
            n = r.headers.get("Content-Length")
            if n is not None and int(n) > MAX_BYTES:
                rec["error"] = f"refused: {n} bytes (>50 MB)"
                append(log, rec)
                return rec
            raw = r.read(MAX_BYTES + 1)
            if len(raw) > MAX_BYTES:
                rec["error"] = "refused: exceeded 50 MB while reading"
                append(log, rec)
                return rec
            rec["bytes"] = len(raw)
            rec["content_type"] = r.headers.get("Content-Type", "")
            text = raw.decode("utf-8", errors="replace")
    except urllib.error.HTTPError as e:
        rec["status"] = e.code
        rec["error"] = f"HTTPError {e.code} {e.reason}"
        rec["text"] = e.read(200000).decode("utf-8", errors="replace") if store else ""
        append(log, rec)
        return rec
    except Exception as e:  # noqa: BLE001 - URLError, timeouts, TLS
        rec["status"] = None
        rec["error"] = f"{type(e).__name__}: {e}"
        append(log, rec)
        return rec
    if store:
        rec["text"] = text
    append(log, rec)
    rec["text"] = text
    return rec


def gds_search(area: str, term: str, retmax: int = 60) -> None:
    """GEO DataSets esearch + esummary. A search listing, not an opened record (D8)."""
    q = urllib.parse.urlencode({"db": "gds", "term": term, "retmax": retmax, "retmode": "json"})
    r = fetch(f"esearch {term}", EUTILS + "esearch.fcgi?" + q, log="_scratch.json", store=False)
    es = json.loads(r.get("text") or "{}").get("esearchresult", {})
    ids = es.get("idlist", [])
    hits = []
    if ids:
        time.sleep(0.4)
        q = urllib.parse.urlencode({"db": "gds", "id": ",".join(ids), "retmode": "json"})
        r2 = fetch(f"esummary {term}", EUTILS + "esummary.fcgi?" + q, log="_scratch.json",
                   store=False)
        res = json.loads(r2.get("text") or "{}").get("result", {})
        for uid in res.get("uids", []):
            d = res[uid]
            hits.append({k: d.get(k) for k in ("accession", "taxon", "title", "summary",
                                                "gdstype", "n_samples", "pdat")})
    append("searches.json", {"area": area, "engine": "geo_gds", "term": term,
                             "count": int(es.get("count", 0)), "hits": hits})
    print(f"## [{area}] GEO '{term}': {es.get('count')} hits, showing {len(hits)}")
    for h in hits:
        if h["accession"].startswith("GSE"):
            print(f"{h['accession']} | {h['taxon']} | n={h['n_samples']} | {h['pdat']} | "
                  f"{h['title']}\n    {(h['summary'] or '')[:400]}")


def epmc_search(area: str, query: str, n: int = 25) -> None:
    """Europe PMC search listing (title, journal, year, abstract). Not an opened record (D8)."""
    q = urllib.parse.urlencode({"query": query, "format": "json", "pageSize": n,
                                "resultType": "core"})
    r = fetch(f"epmc {query}", "https://www.ebi.ac.uk/europepmc/webservices/rest/search?" + q,
              log="_scratch.json", store=False)
    res = json.loads(r.get("text") or "{}")
    hits = [{"id": h.get("id"), "source": h.get("source"), "doi": h.get("doi"),
             "title": h.get("title"), "journal": h.get("journalTitle"),
             "year": h.get("pubYear"), "abstract": h.get("abstractText", "")}
            for h in res.get("resultList", {}).get("result", [])]
    append("searches.json", {"area": area, "engine": "europepmc", "term": query,
                             "count": res.get("hitCount"), "hits": hits})
    print(f"## [{area}] EPMC '{query}': {res.get('hitCount')} hits")
    for h in hits:
        print(f"{h['source']}:{h['id']} | {h['year']} | {h['journal']} | doi {h['doi']}\n"
              f"    {h['title']}\n    {(h['abstract'] or '')[:700]}")


def geo_text(area: str, acc: str) -> dict:
    q = urllib.parse.urlencode({"acc": acc, "targ": "self", "form": "text", "view": "quick"})
    return fetch(acc, "https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?" + q,
                 log="geo_records.json", area=area)


def reproduce() -> int:
    """Pre-registered reproduction of the stored SCOPEFIB / SCOPEFIB2 record (PREREG.flag)."""
    rows: list[dict] = []

    def chk(key, expected, observed):
        rows.append({"check": key, "expected": expected, "observed": observed,
                     "match": "yes" if str(expected) == str(observed) else "NO"})

    for key, path, n in (("R1", "results/scopefib/datasets.csv", 23),
                         ("R2", "results/scopefib2/datasets.csv", 1)):
        with open(REPO / path, newline="", encoding="utf-8") as fh:
            chk(f"{key} {path} data rows", n, len(list(csv.DictReader(fh))))

    scoping = (REPO / "SCOPING_FIBROSIS.md").read_text(encoding="utf-8").split("\n")
    readings = [ln.split(". ", 1)[1].strip().strip("`") for ln in scoping[6:8]]
    progress = (REPO / "PROGRESS_SCOPEFIB.md").read_text(encoding="utf-8")
    for i, r in enumerate(readings, 1):
        chk(f"R3 SCOPING_FIBROSIS.md line {6 + i} reading in PROGRESS_SCOPEFIB.md", True,
            r in progress)
        rows.append({"check": f"R3 reading {i} text", "expected": r, "observed": r,
                     "match": "yes"})
    nogo = readings[1]
    chk("R3b reading 2 is the NO-GO", True, nogo.startswith("NO-GO: no human dataset"))
    scoping2 = (REPO / "SCOPING_FIBROSIS2.md").read_text(encoding="utf-8")
    chk("R4 SCOPEFIB NO-GO verbatim in SCOPING_FIBROSIS2.md", True, nogo in scoping2)
    chk("R4b SCOPING_FIBROSIS2.md says NO-GO not superseded", True,
        "NO-GO is not superseded" in scoping2)

    src = (REPO / "results/scopefib/PREREG.flag").read_text(encoding="utf-8").split("\n")
    carried = "\n".join(src[9:87])
    chk("R5 carried block n lines", 78, len(carried.split("\n")))
    mine = (REPO / "results/scopedis/PREREG.flag").read_text(encoding="utf-8")
    chk("R5b carried lines 10-87 byte-identical in SCOPEDIS flag", True, carried in mine)
    prev = (REPO / "results/scopefib2/PREREG.flag").read_text(encoding="utf-8")
    chk("R6 same 78 lines byte-identical in SCOPEFIB2 flag", True, carried in prev)

    with open(OUT / "reproduction.csv", "w", newline="", encoding="utf-8") as fh:
        w = csv.DictWriter(fh, fieldnames=["check", "expected", "observed", "match"])
        w.writeheader()
        w.writerows(rows)
    bad = [r for r in rows if r["match"] == "NO"]
    for r in rows:
        print(f"{r['match']:>3}  {r['check']}")
    print(f"{len(rows)} checks, {len(bad)} mismatches")
    return 1 if bad else 0


def main(argv: list[str]) -> int:
    OUT.mkdir(parents=True, exist_ok=True)
    mode, args = argv[0], argv[1:]
    if mode == "reproduce":
        return reproduce()
    if mode == "gds":
        gds_search(args[0], args[1], int(args[2]) if len(args) > 2 else 60)
    elif mode == "epmc":
        epmc_search(args[0], args[1], int(args[2]) if len(args) > 2 else 25)
    elif mode == "series":
        for acc in args[1:]:
            r = geo_text(args[0], acc)
            print(f"===== {acc} [{r.get('status')}]\n{r.get('text', r.get('error', ''))}")
            time.sleep(0.4)
    elif mode == "url":
        r = fetch(args[1], args[2], area=args[0])
        print(f"[{r.get('status')}] {r['url']}")
        print(r.get("text", r.get("error", ""))[:200000])
    else:
        raise SystemExit(__doc__)
    return 0


if __name__ == "__main__":
    raise SystemExit(main(sys.argv[1:]))
