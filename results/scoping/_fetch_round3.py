"""Resolve Joung subseries designs, GTEx public URLs, Sarkar SRA, CXG skin metadata."""
from __future__ import annotations

import json
import time
import urllib.parse
import urllib.request
from collections import defaultdict
from pathlib import Path

OUT = Path(__file__).resolve().parent
UA = "AgeIdentitySeparabilityBenchmark/scoping (metadata only)"


def get_bytes(url: str, timeout: int = 90) -> bytes:
    req = urllib.request.Request(url, headers={"User-Agent": UA})
    with urllib.request.urlopen(req, timeout=timeout) as r:
        return r.read()


def geo_text(acc: str, view: str = "quick") -> str:
    url = (
        "https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?"
        + urllib.parse.urlencode({"acc": acc, "targ": "self", "form": "text", "view": view})
    )
    return get_bytes(url).decode("utf-8", errors="replace")


def parse_kv(txt: str) -> dict:
    d = defaultdict(list)
    for line in txt.splitlines():
        if line.startswith("!") and " = " in line:
            k, v = line.split(" = ", 1)
            d[k].append(v)
    return dict(d)


out = {"joung_subseries": {}, "gtex": {}, "sra": {}, "cxg_skin": []}

joung = [
    "GSE217460",
    "GSE216479",
    "GSE216595",
    "GSE219000",
    "GSE219058",
    "GSE216463",
    "GSE216601",
    "GSE216602",
]
for acc in joung:
    print("JOUNG", acc, flush=True)
    try:
        s = parse_kv(geo_text(acc, "quick"))
        out["joung_subseries"][acc] = {
            "title": (s.get("!Series_title") or [None])[0],
            "overall_design": (s.get("!Series_overall_design") or [None])[0],
            "summary": (s.get("!Series_summary") or [None])[0],
            "n_sample_ids": len(s.get("!Series_sample_id") or []),
            "type": s.get("!Series_type"),
        }
        print(" ", out["joung_subseries"][acc]["title"], flush=True)
        print("  design:", (out["joung_subseries"][acc]["overall_design"] or "")[:240], flush=True)
        time.sleep(0.35)
    except Exception as e:
        out["joung_subseries"][acc] = {"error": f"{type(e).__name__}: {e}"}
        print("  ERR", e, flush=True)
        time.sleep(1)

# GTEx public annotation URLs (tiny phenotype file)
gtex_urls = {
    "subject_phenotypes": "https://storage.googleapis.com/adult-gtex/annotations/v10/metadata-files/GTEx_Analysis_v10_Annotations_SubjectPhenotypesDS.txt",
    "sample_attributes": "https://storage.googleapis.com/adult-gtex/annotations/v10/metadata-files/GTEx_Analysis_v10_Annotations_SampleAttributesDS.txt",
    "portal_adult_metadata": "https://gtexportal.org/home/downloads/adult-gtex/metadata",
    "dbgap": "https://www.ncbi.nlm.nih.gov/projects/gap/cgi-bin/study.cgi?study_id=phs000424.v10.p1",
}
for name, url in gtex_urls.items():
    print("GTEX", name, flush=True)
    try:
        req = urllib.request.Request(url, headers={"User-Agent": UA}, method="GET")
        with urllib.request.urlopen(req, timeout=60) as r:
            # only first 800 bytes
            chunk = r.read(800)
            out["gtex"][name] = {
                "url": url,
                "final_url": r.geturl(),
                "status": getattr(r, "status", None),
                "content_type": r.headers.get("Content-Type"),
                "content_length": r.headers.get("Content-Length"),
                "head": chunk.decode("utf-8", errors="replace"),
            }
        print("  ok", out["gtex"][name]["content_type"], "len", out["gtex"][name]["content_length"], flush=True)
        time.sleep(0.2)
    except Exception as e:
        out["gtex"][name] = {"url": url, "error": f"{type(e).__name__}: {e}"}
        print("  ERR", e, flush=True)

# Sarkar RNA SRA
for acc in ("PRJNA598923", "PRJNA1108254", "PRJNA1263211"):
    print("SRA", acc, flush=True)
    url = (
        "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi?"
        + urllib.parse.urlencode({"db": "bioproject", "term": acc, "retmode": "json"})
    )
    try:
        d = json.loads(get_bytes(url).decode())
        ids = d.get("esearchresult", {}).get("idlist", [])
        rec = {"accession": acc, "bioproject_esearch": d.get("esearchresult")}
        if ids:
            sm_url = (
                "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esummary.fcgi?"
                + urllib.parse.urlencode({"db": "bioproject", "id": ids[0], "retmode": "json"})
            )
            rec["esummary"] = json.loads(get_bytes(sm_url).decode())
        # SRA experiment count
        sra_url = (
            "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi?"
            + urllib.parse.urlencode({"db": "sra", "term": acc, "retmode": "json", "retmax": 5})
        )
        rec["sra_esearch"] = json.loads(get_bytes(sra_url).decode()).get("esearchresult")
        out["sra"][acc] = rec
        print("  sra count", rec.get("sra_esearch", {}).get("count"), flush=True)
        time.sleep(0.35)
    except Exception as e:
        out["sra"][acc] = {"error": f"{type(e).__name__}: {e}"}
        print("  ERR", e, flush=True)

# CXG datasets API: search skin in titles via collections is heavy; use datasets endpoint page
print("CXG", flush=True)
try:
    url = "https://api.cellxgene.cziscience.com/curation/v1/datasets"
    raw = get_bytes(url, timeout=120)
    datasets = json.loads(raw.decode())
    hits = []
    for ds in datasets:
        title = (ds.get("title") or "") + " " + (ds.get("name") or "")
        tissues = []
        for t in ds.get("tissue") or []:
            if isinstance(t, dict):
                tissues.append(t.get("label") or "")
            else:
                tissues.append(str(t))
        tissue_s = " ".join(tissues).lower()
        blob = (title + " " + tissue_s).lower()
        if any(k in blob for k in ("skin", "dermis", "fibroblast", "keratinocyte")):
            n_donors = None
            for x in ds.get("donor_id") or []:
                pass
            # schema varies; keep key fields
            hits.append(
                {
                    "dataset_id": ds.get("dataset_id") or ds.get("id"),
                    "title": ds.get("title") or ds.get("name"),
                    "n_cells": ds.get("cell_count") or ds.get("n_cells"),
                    "collection_id": ds.get("collection_id"),
                    "tissues": tissues[:12],
                    "organism": ds.get("organism"),
                    "assay": ds.get("assay"),
                    "keys": sorted(ds.keys())[:40],
                }
            )
    out["cxg_skin"] = hits[:40]
    out["cxg_n_datasets_total"] = len(datasets)
    out["cxg_n_skinish"] = len(hits)
    print("  cxg datasets", len(datasets), "skinish", len(hits), flush=True)
except Exception as e:
    out["cxg_error"] = f"{type(e).__name__}: {e}"
    print("  ERR", e, flush=True)

(OUT / "round3.json").write_text(json.dumps(out, indent=2), encoding="utf-8")
print("wrote round3.json")
