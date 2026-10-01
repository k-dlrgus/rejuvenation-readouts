"""Fetch remaining high-priority GEO records and dump sample titles."""
from __future__ import annotations

import json
import time
import urllib.parse
import urllib.request
from collections import defaultdict
from pathlib import Path

OUT = Path(__file__).resolve().parent
UA = "AgeIdentitySeparabilityBenchmark/scoping (metadata only)"
accessions = [
    "GSE325735",
    "GSE326022",
    "GSE297984",
    "GSE201710",
    "GSE247179",
    "GSE247178",
    "GSE300230",
    "GSE307377",
    "GSE304044",
    "GSE190665",
    "GSE237269",
    "GSE254389",
    "GSE261783",
    "GSE280438",
    "GSE276291",
    "GSE117720",
    "GSE102983",
    "GSE184546",
    "GSE165450",
]


def get(url: str, timeout: int = 120) -> bytes:
    req = urllib.request.Request(url, headers={"User-Agent": UA})
    with urllib.request.urlopen(req, timeout=timeout) as r:
        return r.read()


def esearch(term: str, db: str = "gds") -> dict:
    url = (
        "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi?"
        + urllib.parse.urlencode({"db": db, "term": term, "retmax": 5, "retmode": "json"})
    )
    return json.loads(get(url).decode())


def esummary_gds(acc: str) -> dict:
    d = esearch(f"{acc}[ACCN]")
    ids = d.get("esearchresult", {}).get("idlist", [])
    if not ids:
        return {"accession": acc, "error": "no_gds_id"}
    url = (
        "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esummary.fcgi?"
        + urllib.parse.urlencode({"db": "gds", "id": ids[0], "retmode": "json"})
    )
    sm = json.loads(get(url).decode())
    rec = sm.get("result", {}).get(ids[0], {})
    samples = rec.get("samples") or []
    return {
        "accession": rec.get("accession") or acc,
        "gds_id": ids[0],
        "title": rec.get("title"),
        "summary": rec.get("summary"),
        "n_samples": rec.get("n_samples"),
        "gpl": rec.get("gpl"),
        "gdstype": rec.get("gdstype"),
        "taxon": rec.get("taxon"),
        "pdat": rec.get("pdat"),
        "pubmedids": rec.get("pubmedids"),
        "n_sample_titles_returned": len(samples),
        "sample_titles": [{"accession": s.get("accession"), "title": s.get("title")} for s in samples],
        "ftplink": rec.get("ftplink"),
        "bioproject": rec.get("bioproject"),
        "subsetinfo": rec.get("subsetinfo"),
        "ssinfo": rec.get("ssinfo"),
    }


def geo_text(acc: str, view: str = "quick") -> str:
    url = (
        "https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?"
        + urllib.parse.urlencode({"acc": acc, "targ": "self", "form": "text", "view": view})
    )
    return get(url).decode("utf-8", errors="replace")


def parse_series(txt: str) -> dict:
    d = defaultdict(list)
    for line in txt.splitlines():
        if line.startswith("!") and " = " in line:
            k, v = line.split(" = ", 1)
            d[k].append(v)
    return dict(d)


out = {}
for acc in accessions:
    print("RESOLVE", acc, flush=True)
    rec = {"accession": acc}
    try:
        rec["esummary"] = esummary_gds(acc)
        print(" ", rec["esummary"].get("title"), "n", rec["esummary"].get("n_samples"),
              "titles_returned", rec["esummary"].get("n_sample_titles_returned"), flush=True)
        time.sleep(0.34)
    except Exception as e:
        rec["esummary_error"] = f"{type(e).__name__}: {e}"
        print("  esummary ERR", e, flush=True)
        time.sleep(1)
    try:
        txt = geo_text(acc, "quick")
        s = parse_series(txt)
        rec["title"] = (s.get("!Series_title") or [None])[0]
        rec["overall_design"] = (s.get("!Series_overall_design") or [None])[0]
        rec["summary"] = (s.get("!Series_summary") or [None])[0]
        rec["type"] = s.get("!Series_type")
        rec["relation"] = s.get("!Series_relation")
        rec["n_sample_ids"] = len(s.get("!Series_sample_id") or [])
        rec["sample_ids_head"] = (s.get("!Series_sample_id") or [])[:20]
        rec["pubmed"] = s.get("!Series_pubmed_id")
        rec["subseries"] = [x for x in (s.get("!Series_relation") or []) if "SubSeries" in x or "SuperSeries" in x or "GSE" in x]
        print("  design:", (rec["overall_design"] or "")[:300], flush=True)
        time.sleep(0.34)
    except Exception as e:
        rec["geo_text_error"] = f"{type(e).__name__}: {e}"
        print("  text ERR", e, flush=True)
        time.sleep(1)
    out[acc] = rec
    (OUT / "geo_records_round2.json").write_text(json.dumps(out, indent=2), encoding="utf-8")

print("done", len(out))
