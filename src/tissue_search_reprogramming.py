"""T1 reprogramming: scan the already-fetched CXG human listing for OSK/OSKM / Yamanaka / partial
reprogramming / time-course titles, and fetch GEO SOFT metadata for every unique accession the GEO
search returned under a reprogramming query. Nothing is assumed — if GEO has no Series record we say so.

Outputs: results/tissue/t1_reprogramming.txt, t1_cxg_reprogramming.csv, t1_geo_reprogramming_meta.csv
"""
import re
import sys
import time
from pathlib import Path

import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parent))
from tissue_common import TISSUE_DIR, Logger  # noqa: E402
from geo_meta import summarize, parse_series, fetch_soft  # noqa: E402

RX = re.compile(
    r"OSKM|\bOSK\b|Yamanaka|partial reprogramming|in vivo reprogramming|"
    r"reprogramming time|OSK treatment|Oct4.*Sox2.*Klf4",
    re.I,
)
TISSUE_RX = {
    "retina": re.compile(r"retina|macula|photoreceptor|RPE|choroid", re.I),
    "brain": re.compile(r"brain|cortex|neuron|hippocamp|DLPFC|glia", re.I),
    "skeletal muscle": re.compile(r"muscle|myofiber|myocyte|vastus", re.I),
    "heart": re.compile(r"heart|cardiac|cardiomyocyte|myocardium", re.I),
    "kidney": re.compile(r"kidney|renal|nephron", re.I),
    "liver": re.compile(r"liver|hepatic|hepatocyte", re.I),
}


def tissues_mentioned(text: str):
    return [g for g, rx in TISSUE_RX.items() if rx.search(text or "")]


def main():
    log = Logger(TISSUE_DIR / "t1_reprogramming.txt")
    cxg = pd.read_csv(TISSUE_DIR / "t1_cxg_all_human_datasets.csv")
    blob = (cxg.title.fillna("") + " | " + cxg.collection_name.fillna("") + " | " + cxg.tissues.fillna(""))
    hit = cxg[blob.str.contains(RX)]
    hit = hit.copy()
    hit["tissues_mentioned"] = [
        ";".join(tissues_mentioned(b)) for b in blob[hit.index]
    ]
    hit.to_csv(TISSUE_DIR / "t1_cxg_reprogramming.csv", index=False)
    log(f"CXG human listing: {len(cxg)} datasets; title/collection match OSK/OSKM/reprogramming: {len(hit)}")
    if len(hit):
        log(hit[["dataset_id", "title", "collection_name", "n_donors", "cell_count", "tissues", "tissues_mentioned"]].to_string(index=False))
    else:
        log("  (none — CELLxGENE Discover has no human OSK/OSKM time-course under these title terms)")

    geo = pd.read_csv(TISSUE_DIR / "t1_geo_hits.csv") if (TISSUE_DIR / "t1_geo_hits.csv").exists() else pd.DataFrame()
    repro_q = geo[geo["query"].str.contains("reprogramming|OSK", case=False, na=False)] if len(geo) else geo
    accs = sorted(repro_q.accession.dropna().unique()) if len(repro_q) else []
    log(f"\nGEO reprogramming-query accessions to fetch SOFT for: {len(accs)}")
    meta_rows = []
    for acc in accs:
        log("\n")
        try:
            text = fetch_soft(acc)
        except Exception as e:
            log(f"[{acc}] fetch failed: {e}")
            continue
        if "!Series_title" not in text:
            log(f"[{acc}] NOT FOUND in GEO (no Series record). Head:\n{text[:300]}")
            meta_rows.append(dict(accession=acc, found=False))
            time.sleep(0.35)
            continue
        m = parse_series(text)
        title = (m.get("Series_title") or [""])[0]
        organism = "; ".join(m.get("Series_sample_organism", m.get("Series_platform_organism", ["?"])))
        n_samples = len(m.get("Series_sample_id", []))
        pubmed = "; ".join(m.get("Series_pubmed_id", ["none listed"]))
        summary = " ".join(m.get("Series_summary", []))
        design = " ".join(m.get("Series_overall_design", []))
        supp = " | ".join(m.get("Series_supplementary_file", [])[:8])
        mentioned = tissues_mentioned(title + " " + summary + " " + design)
        is_sc = bool(re.search(r"single.cell|single.nucleus|scRNA|snRNA", title + summary + " ".join(m.get("Series_type", [])), re.I))
        log(f"ACCESSION: {acc}")
        log(f"  Title:     {title}")
        log(f"  Organism:  {organism}")
        log(f"  Type:      {'; '.join(m.get('Series_type', ['?']))}")
        log(f"  N samples: {n_samples}   PubMed: {pubmed}")
        log(f"  Tissues mentioned (priority list): {mentioned or ['(none of retina/brain/muscle/heart/kidney/liver)']}")
        log(f"  single-cell/nucleus in title/summary/type: {is_sc}")
        log(f"  Summary:   {summary[:500]}")
        log(f"  Design:    {design[:400]}")
        log(f"  Suppl:     {supp[:400]}")
        meta_rows.append(dict(
            accession=acc, found=True, title=title, organism=organism, n_samples=n_samples,
            pubmed=pubmed, tissues_mentioned=";".join(mentioned), is_scrna=is_sc,
            series_type="; ".join(m.get("Series_type", [])),
            summary=summary[:800], design=design[:500],
        ))
        time.sleep(0.35)
    pd.DataFrame(meta_rows).to_csv(TISSUE_DIR / "t1_geo_reprogramming_meta.csv", index=False)

    # coverage table: which priority tissues have an aging atlas (from curated CXG) AND reprogramming data
    prim = TISSUE_DIR / "t1_cxg_primary_candidates.csv"
    aging = pd.read_csv(prim) if prim.exists() else pd.DataFrame()
    log("\n=== tissue coverage: aging atlas (CXG primary) vs reprogramming (GEO SOFT, any organism) ===")
    for g in ["retina", "brain", "skeletal muscle", "heart", "kidney", "liver"]:
        n_aging = int((aging.tissue_group == g).sum()) if len(aging) else 0
        n_full = int(((aging.tissue_group == g) & (~aging.is_subset)).sum()) if len(aging) else 0
        geo_hits = [r for r in meta_rows if r.get("found") and g in (r.get("tissues_mentioned") or "")]
        log(f"  {g:<18} CXG primary datasets={n_aging} (full={n_full});  "
            f"GEO reprogramming series mentioning tissue={len(geo_hits)} "
            f"{[r['accession'] for r in geo_hits]}")
    log.close()


if __name__ == "__main__":
    main()
