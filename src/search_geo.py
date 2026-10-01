"""Search NCBI GEO (gds database) for human single-cell PBMC/immune aging datasets.
Prints what GEO actually returns; makes no assumptions about accessions.
"""
import time
import requests
import xml.etree.ElementTree as ET

EUTILS = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils"

QUERIES = [
    '(aging[Title] OR ageing[Title] OR "age-related"[Title] OR "healthy aging"[All Fields]) AND ("single-cell"[All Fields] OR "single cell"[All Fields] OR scRNA-seq[All Fields]) AND (PBMC[All Fields] OR "peripheral blood"[All Fields] OR immune[All Fields]) AND "Homo sapiens"[Organism] AND gse[Entry Type]',
    '("single-cell eQTL"[All Fields] OR OneK1K[All Fields]) AND "Homo sapiens"[Organism] AND gse[Entry Type]',
    '(supercentenarian[All Fields] OR centenarian[All Fields]) AND ("single-cell"[All Fields] OR scRNA-seq[All Fields]) AND "Homo sapiens"[Organism] AND gse[Entry Type]',
]


def esearch(q, retmax=40):
    r = requests.get(f"{EUTILS}/esearch.fcgi", params=dict(db="gds", term=q, retmax=retmax, retmode="json"), timeout=60)
    r.raise_for_status()
    return r.json()["esearchresult"]


def esummary(ids):
    if not ids:
        return []
    r = requests.get(f"{EUTILS}/esummary.fcgi", params=dict(db="gds", id=",".join(ids), retmode="json"), timeout=60)
    r.raise_for_status()
    res = r.json()["result"]
    return [res[i] for i in res["uids"]]


def main():
    seen = set()
    for q in QUERIES:
        print("=" * 100)
        print("QUERY:", q)
        es = esearch(q)
        print(f"  GEO returned {es['count']} records (showing up to {len(es['idlist'])})")
        time.sleep(0.4)
        for s in esummary(es["idlist"]):
            acc = s.get("accession")
            if acc in seen:
                continue
            seen.add(acc)
            print(f"  {acc:<11} n_samples={s.get('n_samples'):>5}  {s.get('pdat')}  {s.get('gdstype','')[:40]:<40} | {s.get('title','')[:110]}")
        time.sleep(0.4)


if __name__ == "__main__":
    main()
