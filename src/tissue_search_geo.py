"""T1 GEO search: human (and separately mouse) single-cell / single-nucleus aging atlases in the
priority low-turnover tissues, plus a dedicated OSK/OSKM / partial-reprogramming query.

Every accession printed was returned by NCBI E-utilities. Nothing is assumed. Rate-limited.
Outputs: results/tissue/t1_geo_search.txt, results/tissue/t1_geo_hits.csv
"""
import sys
import time
from pathlib import Path

import pandas as pd
import requests

sys.path.insert(0, str(Path(__file__).resolve().parent))
from tissue_common import TISSUE_DIR, Logger  # noqa: E402

EUTILS = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils"

TISSUE_TERMS = {
    "retina": 'retina OR macular OR "retinal pigment" OR photoreceptor',
    "brain": 'brain OR cortex OR cortical OR hippocampus OR "prefrontal cortex" OR DLPFC',
    "skeletal muscle": '"skeletal muscle" OR myofiber OR vastus OR quadriceps',
    "heart": 'heart OR cardiac OR cardiomyocyte OR myocardium',
    "kidney": 'kidney OR renal OR nephron',
    "liver": 'liver OR hepatic OR hepatocyte',
}

AGING = '(aging[Title] OR ageing[Title] OR "age-related"[Title] OR "healthy aging" OR "human aging" OR "ageing atlas")'
SC = '("single-cell" OR "single cell" OR "single-nucleus" OR "single nucleus" OR scRNA-seq OR snRNA-seq OR "snRNA-seq")'
REPROG = '(OSKM OR OSK OR "Oct4 Sox2 Klf4" OR "partial reprogramming" OR "Yamanaka factor" OR "in vivo reprogramming" OR "OSK treatment")'


def esearch(q, retmax=40):
    r = requests.get(f"{EUTILS}/esearch.fcgi",
                     params=dict(db="gds", term=q, retmax=retmax, retmode="json"), timeout=90)
    r.raise_for_status()
    return r.json()["esearchresult"]


def esummary(ids):
    if not ids:
        return []
    r = requests.get(f"{EUTILS}/esummary.fcgi",
                     params=dict(db="gds", id=",".join(ids), retmode="json"), timeout=90)
    r.raise_for_status()
    res = r.json()["result"]
    return [res[i] for i in res["uids"] if i in res]


def run_query(log, rows, seen, query_name, tissue, organism, q, retmax=40):
    log("=" * 100)
    log("QUERY:", query_name)
    log("TERM: ", q)
    try:
        es = esearch(q, retmax=retmax)
    except requests.RequestException as e:
        log(f"  NCBI error: {e}")
        return
    n = int(es.get("count", 0))
    ids = es.get("idlist", [])
    log(f"  GEO returned {n} records (showing {len(ids)})")
    time.sleep(0.4)
    try:
        summaries = esummary(ids)
    except requests.RequestException as e:
        log(f"  esummary error: {e}")
        return
    for s in summaries:
        acc = s.get("accession")
        title = s.get("title", "")
        rec = dict(
            query=query_name, tissue=tissue, organism=organism, accession=acc,
            uid=s.get("uid"), n_samples=s.get("n_samples"), pdat=s.get("pdat"),
            gdstype=s.get("gdstype", ""), taxon=s.get("taxon", ""),
            title=title, summary=(s.get("summary") or "")[:400],
        )
        rows.append(rec)
        flag = "  [dup]" if acc in seen else ""
        seen.add(acc)
        def _safe(x):
            return str(x).encode("ascii", "replace").decode("ascii")
        log(f"  {acc:<12} n={str(s.get('n_samples')):>5}  {s.get('pdat')}  "
            f"{_safe((s.get('gdstype') or '')[:36]):<36} | {_safe(title[:110])}{flag}")
    time.sleep(0.4)


def main():
    log = Logger(TISSUE_DIR / "t1_geo_search.txt")
    rows, seen = [], set()

    for tissue, tterm in TISSUE_TERMS.items():
        q = f"{AGING} AND {SC} AND ({tterm}) AND \"Homo sapiens\"[Organism] AND gse[Entry Type]"
        run_query(log, rows, seen, f"human aging sc {tissue}", tissue, "human", q)

    for tissue, tterm in TISSUE_TERMS.items():
        q = f"{AGING} AND {SC} AND ({tterm}) AND \"Mus musculus\"[Organism] AND gse[Entry Type]"
        run_query(log, rows, seen, f"mouse aging sc {tissue}", tissue, "mouse", q, retmax=20)

    qh = f"{REPROG} AND {SC} AND \"Homo sapiens\"[Organism] AND gse[Entry Type]"
    run_query(log, rows, seen, "human OSK/OSKM/reprogramming sc", "any", "human", qh, retmax=50)
    qm = f"{REPROG} AND {SC} AND \"Mus musculus\"[Organism] AND gse[Entry Type]"
    run_query(log, rows, seen, "mouse OSK/OSKM/reprogramming sc", "any", "mouse", qm, retmax=50)
    # tissue-restricted reprogramming (human + mouse together so we do not miss a series)
    for tissue, tterm in TISSUE_TERMS.items():
        q = f"{REPROG} AND ({tterm}) AND {SC} AND gse[Entry Type]"
        run_query(log, rows, seen, f"reprogramming sc {tissue}", tissue, "any", q, retmax=20)

    df = pd.DataFrame(rows)
    df.to_csv(TISSUE_DIR / "t1_geo_hits.csv", index=False)
    log(f"\n{len(df)} GEO records written ({df.accession.nunique()} unique accessions)")
    log.close()
    return df


if __name__ == "__main__":
    main()
