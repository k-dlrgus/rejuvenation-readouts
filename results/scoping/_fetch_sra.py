"""SRA runinfo + BioProject summaries for Sarkar RNA and Southard Perturb-seq."""
from __future__ import annotations

import csv
import io
import json
import urllib.parse
import urllib.request
from pathlib import Path

OUT = Path(__file__).resolve().parent
UA = "AgeIdentitySeparabilityBenchmark/scoping (metadata only)"


def get_bytes(url: str, timeout: int = 90) -> bytes:
    req = urllib.request.Request(url, headers={"User-Agent": UA})
    with urllib.request.urlopen(req, timeout=timeout) as r:
        return r.read()


def esummary(db: str, term_acc: str) -> dict:
    url = (
        "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi?"
        + urllib.parse.urlencode({"db": db, "term": term_acc, "retmode": "json"})
    )
    d = json.loads(get_bytes(url).decode())
    ids = d.get("esearchresult", {}).get("idlist", [])
    if not ids:
        return {"error": "no id", "esearch": d}
    url = (
        "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esummary.fcgi?"
        + urllib.parse.urlencode({"db": db, "id": ids[0], "retmode": "json"})
    )
    return json.loads(get_bytes(url).decode())


out = {}
for acc in ("PRJNA598923", "PRJNA1108254", "PRJNA1263211"):
    print("BIOSUM", acc, flush=True)
    rec = {"accession": acc}
    try:
        rec["bioproject"] = esummary("bioproject", acc)
    except Exception as e:
        rec["bioproject_error"] = f"{type(e).__name__}: {e}"
    # runinfo
    url = f"https://trace.ncbi.nlm.nih.gov/Traces/sra-db-be/sra/sra.cgi?save=efetch&db=sra&rettype=runinfo&term={acc}"
    # NCBI documented:
    url2 = (
        "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi?"
        + urllib.parse.urlencode({"db": "sra", "term": acc, "retmode": "json", "retmax": 30})
    )
    try:
        sra = json.loads(get_bytes(url2).decode())
        ids = sra.get("esearchresult", {}).get("idlist", [])
        rec["sra_count"] = sra.get("esearchresult", {}).get("count")
        rec["sra_ids_head"] = ids
        if ids:
            sm_url = (
                "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esummary.fcgi?"
                + urllib.parse.urlencode({"db": "sra", "id": ",".join(ids[:20]), "retmode": "json"})
            )
            sm = json.loads(get_bytes(sm_url).decode())
            rec["sra_summaries_keys"] = list((sm.get("result") or {}).keys())[:5]
            uids = (sm.get("result") or {}).get("uids", ids[:20])
            rows = []
            for uid in uids:
                r = (sm.get("result") or {}).get(uid, {})
                expxml = r.get("expxml") or ""
                runs = r.get("runs") or ""
                rows.append(
                    {
                        "uid": uid,
                        "title": r.get("title"),
                        "expxml_head": expxml[:800],
                        "runs_head": runs[:400],
                        "keys": sorted(r.keys()),
                    }
                )
            rec["sra_rows"] = rows
            for row in rows[:8]:
                print(" ", row.get("title"), flush=True)
                print("   xml", (row.get("expxml_head") or "")[:200], flush=True)
    except Exception as e:
        rec["sra_error"] = f"{type(e).__name__}: {e}"
        print(" ERR", e, flush=True)
    out[acc] = rec

(OUT / "sra_meta.json").write_text(json.dumps(out, indent=2), encoding="utf-8")
print("wrote sra_meta.json")
