"""Fetch GEO series text + SRA runinfo for a short accession list. Metadata only."""
from __future__ import annotations

import json
import time
import urllib.parse
import urllib.request
from pathlib import Path

OUT = Path(__file__).resolve().parent
UA = "AgeIdentitySeparabilityBenchmark/scoping (metadata only)"
accessions = [
    "GSE297234",
    "GSE297233",
    "GSE165180",
    "GSE165177",
    "GSE165176",
    "GSE142439",
    "GSE216481",
    "GSE217460",
    "GSE113957",
    "GSE226189",
    "GSE307377",
    "GSE261783",
    "GSE280438",
    "GSE47489",
    "GSE28688",
    "GSE254389",
    "GSE223748",
    "GSE190665",
    "GSE237269",
    "GSE247179",
    "GSE325735",
    "GSE326022",
    "GSE300477",
    "GSE299018",
    "GSE338065",
    "GSE310818",
    "GSE287206",
    "GSE107654",
    "GSE110544",
]


def get(url: str, timeout: int = 120) -> bytes:
    req = urllib.request.Request(url, headers={"User-Agent": UA, "Accept": "text/plain, text/html, */*"})
    with urllib.request.urlopen(req, timeout=timeout) as r:
        return r.read()


def geo_text(acc: str, view: str = "quick") -> str:
    url = (
        "https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?"
        + urllib.parse.urlencode({"acc": acc, "targ": "self", "form": "text", "view": view})
    )
    return get(url).decode("utf-8", errors="replace")


def esearch(term: str, db: str = "gds") -> dict:
    url = (
        "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi?"
        + urllib.parse.urlencode({"db": db, "term": term, "retmax": 20, "retmode": "json"})
    )
    return json.loads(get(url).decode())


def esummary_gds(acc: str) -> dict:
    d = esearch(f"{acc}[ACCN]")
    ids = d.get("esearchresult", {}).get("idlist", [])
    if not ids:
        return {"accession": acc, "error": "no_gds_id", "esearch": d}
    url = (
        "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esummary.fcgi?"
        + urllib.parse.urlencode({"db": "gds", "id": ids[0], "retmode": "json"})
    )
    sm = json.loads(get(url).decode())
    rec = sm.get("result", {}).get(ids[0], {})
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
        "suppfile": rec.get("suppfile"),
        "gse": rec.get("gse"),
        "ssinfo": rec.get("ssinfo"),
        "extrelations": rec.get("extrelations"),
        "samples": rec.get("samples"),
        "subsetinfo": rec.get("subsetinfo"),
        "relations": rec.get("relations"),
        "projects": rec.get("projects"),
        "all_keys": sorted(rec.keys()),
    }


out = {}
for acc in accessions:
    print("RESOLVE", acc, flush=True)
    rec = {"accession": acc}
    try:
        rec["esummary"] = esummary_gds(acc)
        time.sleep(0.34)
    except Exception as e:
        rec["esummary_error"] = f"{type(e).__name__}: {e}"
        print("  esummary ERR", e, flush=True)
        time.sleep(1)
    try:
        txt = geo_text(acc, "quick")
        rec["geo_text_quick"] = txt
        rec["geo_text_nchars"] = len(txt)
        print("  text", len(txt), "samples_line", [ln for ln in txt.splitlines() if "Sample" in ln or "GSM" in ln][:8], flush=True)
        time.sleep(0.34)
    except Exception as e:
        rec["geo_text_error"] = f"{type(e).__name__}: {e}"
        print("  text ERR", e, flush=True)
        time.sleep(1)
    out[acc] = rec
    (OUT / "geo_records.json").write_text(json.dumps(out, indent=2), encoding="utf-8")

print("done", len(out))
