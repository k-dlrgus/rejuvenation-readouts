"""Resolve GEO/SRA records for SCOPING_TISSUE. No large downloads."""
from __future__ import annotations

import json
import time
import urllib.parse
import urllib.request
from pathlib import Path

OUT = Path(__file__).resolve().parent
OUT.mkdir(parents=True, exist_ok=True)

UA = "AgeIdentitySeparabilityBenchmark/scoping (local survey; metadata only)"


def get(url: str, timeout: int = 90) -> bytes:
    req = urllib.request.Request(url, headers={"User-Agent": UA})
    with urllib.request.urlopen(req, timeout=timeout) as r:
        return r.read()


def geo_text(acc: str) -> str:
    url = (
        "https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?"
        + urllib.parse.urlencode({"acc": acc, "targ": "self", "form": "text", "view": "full"})
    )
    return get(url).decode("utf-8", errors="replace")


def esearch(term: str, db: str = "gds", retmax: int = 80) -> dict:
    url = (
        "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi?"
        + urllib.parse.urlencode(
            {"db": db, "term": term, "retmax": retmax, "retmode": "json"}
        )
    )
    return json.loads(get(url).decode())


def esummary(ids: list[str], db: str = "gds") -> dict:
    if not ids:
        return {}
    url = (
        "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esummary.fcgi?"
        + urllib.parse.urlencode(
            {"db": db, "id": ",".join(ids), "retmode": "json"}
        )
    )
    return json.loads(get(url).decode())


queries = {
    "partial_reprogramming_human": '"partial reprogramming" AND Homo sapiens[Organism]',
    "transient_reprogramming_human": '"transient reprogramming" AND Homo sapiens[Organism]',
    "OSKM_fibroblast_human": "OSKM AND fibroblast AND Homo sapiens[Organism]",
    "OSK_fibroblast_human": "OSK AND fibroblast AND Homo sapiens[Organism]",
    "Yamanaka_fibroblast_human": "Yamanaka AND fibroblast AND Homo sapiens[Organism] AND (time course OR time-course OR day)",
    "rejuvenation_fibroblast_rnaseq": "rejuvenation AND fibroblast AND Homo sapiens[Organism] AND rna seq",
    "CRISPRa_fibroblast_human": "CRISPRa AND fibroblast AND Homo sapiens[Organism]",
    "CRISPR_screen_fibroblast_human": "CRISPR AND screen AND fibroblast AND Homo sapiens[Organism] AND rna",
    "TF_atlas_joung": "transcription factor atlas directed differentiation Joung",
    "GSE216481": "GSE216481[ACCN]",
    "GSE165180": "GSE165180[ACCN]",
    "GSE297234": "GSE297234[ACCN]",
    "GSE142439": "GSE142439[ACCN]",
    "muscle_OSKM_human": "(OSKM OR OSK OR \"partial reprogramming\") AND (muscle OR myoblast OR satellite) AND Homo sapiens[Organism]",
    "adipose_OSKM_human": "(OSKM OR OSK OR \"partial reprogramming\") AND (adipose OR adipocyte) AND Homo sapiens[Organism]",
    "lung_OSKM_human": "(OSKM OR OSK OR \"partial reprogramming\") AND (lung OR fibroblast) AND Homo sapiens[Organism] AND (aged OR aging OR ageing)",
    "PBMC_OSKM_human": "(OSKM OR OSK OR \"partial reprogramming\") AND (PBMC OR blood) AND Homo sapiens[Organism]",
    "ipsc_neuron_aged": "(iPSC OR iPS) AND neuron AND (aged OR aging) AND Homo sapiens[Organism] AND (perturbation OR CRISPR OR overexpression)",
    "LINCS_fibroblast": "LINCS AND fibroblast AND Homo sapiens[Organism]",
}

search_out = {}
for name, term in queries.items():
    print("SEARCH", name, flush=True)
    try:
        d = esearch(term)
        ids = d.get("esearchresult", {}).get("idlist", [])
        count = d.get("esearchresult", {}).get("count")
        search_out[name] = {"term": term, "count": count, "ids": ids, "error": None}
        print(" ", count, "ids", ids[:15], flush=True)
        time.sleep(0.35)
        if ids:
            sm = esummary(ids[:40])
            recs = []
            result = sm.get("result", {})
            uids = result.get("uids", ids[:40])
            for uid in uids:
                rec = result.get(uid, {})
                recs.append(
                    {
                        "uid": uid,
                        "accession": rec.get("accession"),
                        "title": rec.get("title"),
                        "n_samples": rec.get("n_samples"),
                        "gpl": rec.get("gpl"),
                        "gds_type": rec.get("gdstype"),
                        "taxon": rec.get("taxon"),
                        "summary": (rec.get("summary") or "")[:400],
                        "pdat": rec.get("pdat"),
                        "gse": rec.get("gse"),
                        "suppfile": rec.get("suppfile"),
                    }
                )
            search_out[name]["summaries"] = recs
            time.sleep(0.35)
    except Exception as e:
        search_out[name] = {"term": term, "error": f"{type(e).__name__}: {e}"}
        print(" ERR", e, flush=True)
        time.sleep(1)

(OUT / "ncbi_searches.json").write_text(json.dumps(search_out, indent=2), encoding="utf-8")
print("wrote ncbi_searches.json")
