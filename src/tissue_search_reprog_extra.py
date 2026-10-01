"""Extra GEO queries for OSK-in-tissue series that the scRNA-restricted search can miss
(e.g. Lu et al. 2020 OSK vision restoration is often bulk). Every accession is from E-utilities.
"""
import sys
import time
from pathlib import Path

import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parent))
from tissue_common import TISSUE_DIR, Logger  # noqa: E402
from tissue_search_geo import esearch, esummary  # noqa: E402
from geo_meta import fetch_soft, parse_series  # noqa: E402

QUERIES = [
    ("OSK/OSKM + retina/vision (any assay)",
     '(OSK OR OSKM OR Yamanaka) AND (retina OR photoreceptor OR vision OR RPE) AND gse[Entry Type]'),
    ("partial reprogramming [Title]",
     '"partial reprogramming"[Title] AND gse[Entry Type]'),
    ("in vivo OSK + aging",
     '(OSK OR OSKM) AND (aging OR ageing) AND "in vivo" AND gse[Entry Type]'),
]


def main():
    log = Logger(TISSUE_DIR / "t1_reprogramming_extra.txt")
    rows = []
    for name, q in QUERIES:
        log("=" * 90)
        log(name)
        log(q)
        es = esearch(q, retmax=20)
        log(f"  GEO returned {es.get('count')} (showing {len(es.get('idlist', []))})")
        time.sleep(0.4)
        for s in esummary(es.get("idlist", [])):
            acc = s.get("accession")
            title = s.get("title", "")
            log(f"  {acc:<12} n={s.get('n_samples')}  {s.get('taxon','')[:24]:<24} {title[:110]}")
            rows.append(dict(query=name, accession=acc, n_samples=s.get("n_samples"),
                             taxon=s.get("taxon"), title=title, pdat=s.get("pdat")))
        time.sleep(0.4)
    pd.DataFrame(rows).to_csv(TISSUE_DIR / "t1_geo_reprogramming_extra.csv", index=False)
    # SOFT for a shortlist of titles that actually mention a priority tissue + OSK/partial
    want = []
    for r in rows:
        t = (r["title"] or "").lower()
        if any(k in t for k in ("retina", "photoreceptor", "vision", "brain", "cortex", "neuron",
                                "muscle", "heart", "cardiac", "kidney", "liver", "aging", "ageing",
                                "partial reprogramming")):
            want.append(r["accession"])
    want = sorted(set(want))
    log(f"\nSOFT fetch for {len(want)} extra accessions: {want}")
    meta = []
    for acc in want:
        try:
            text = fetch_soft(acc)
        except Exception as e:
            log(f"[{acc}] fetch failed: {e}")
            continue
        if "!Series_title" not in text:
            log(f"[{acc}] NOT FOUND")
            continue
        m = parse_series(text)
        title = (m.get("Series_title") or [""])[0]
        org = "; ".join(m.get("Series_sample_organism", ["?"]))
        n = len(m.get("Series_sample_id", []))
        typ = "; ".join(m.get("Series_type", []))
        pubmed = "; ".join(m.get("Series_pubmed_id", ["none"]))
        summ = " ".join(m.get("Series_summary", []))[:400]
        log(f"\n{acc} | {org} | n={n} | {typ} | PMID {pubmed}\n  {title}\n  {summ}")
        meta.append(dict(accession=acc, title=title, organism=org, n_samples=n,
                         series_type=typ, pubmed=pubmed, summary=summ))
        time.sleep(0.35)
    pd.DataFrame(meta).to_csv(TISSUE_DIR / "t1_geo_reprogramming_extra_meta.csv", index=False)
    log.close()


if __name__ == "__main__":
    main()
