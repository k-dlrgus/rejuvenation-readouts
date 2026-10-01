"""G3 search — verify the three previously returned reprogramming accessions and
look for human-preferred alternatives plus an independent mouse aging brain
dataset for cross-species transfer validation.

Nothing is assumed. Accessions are taken from the fetched GEO / CELLxGENE
records. Does not invent an accession.

Outputs: results/trajectory/g3_search_*.csv / g3_search.txt

Usage: python src/trajectory_g3_search.py
"""
from __future__ import annotations

import re
import sys
import time
from pathlib import Path

import pandas as pd
import requests

sys.path.insert(0, str(Path(__file__).resolve().parent))
from trajectory_common import TRAJ_DIR, TRAJ_SEED, Logger, dump_json, traj_log_banner  # noqa: E402
from geo_meta import fetch_soft, parse_series  # noqa: E402
from search_cxg import API, parse_age_years  # noqa: E402
from search_geo import esearch, esummary  # noqa: E402

# Previously returned by the tissue-stage GEO reprogramming search (FINDINGS_TISSUE.md).
# Re-fetched here; if GEO has no Series record we say so.
CANDIDATE_REPROG = ("GSE276656", "GSE224438", "GSE271794")

GEO_QUERIES = [
    ('human_brain_OSK_sc',
     '(OSK[All Fields] OR OSKM[All Fields] OR Yamanaka[All Fields] OR "partial reprogramming"[All Fields]) '
     'AND (brain[All Fields] OR cortex[All Fields] OR DLPFC[All Fields] OR neuron[All Fields] OR glia[All Fields]) '
     'AND ("single-cell"[All Fields] OR scRNA[All Fields] OR snRNA[All Fields] OR multiome[All Fields]) '
     'AND "Homo sapiens"[Organism] AND gse[Entry Type]'),
    ('human_OSK_sc_any_tissue',
     '(OSK[All Fields] OR OSKM[All Fields] OR "partial reprogramming"[All Fields]) '
     'AND ("single-cell"[All Fields] OR scRNA[All Fields] OR snRNA[All Fields]) '
     'AND "Homo sapiens"[Organism] AND gse[Entry Type]'),
    ('mouse_brain_aging_sc',
     '(aging[Title] OR ageing[Title] OR "aged brain"[All Fields] OR "aging brain"[All Fields]) '
     'AND (brain[All Fields] OR cortex[All Fields] OR prefrontal[All Fields] OR hippocampus[All Fields]) '
     'AND ("single-cell"[All Fields] OR scRNA[All Fields] OR snRNA[All Fields]) '
     'AND "Mus musculus"[Organism] AND gse[Entry Type]'),
    ('mouse_brain_OSK_sc',
     '(OSK[All Fields] OR OSKM[All Fields] OR Yamanaka[All Fields] OR "partial reprogramming"[All Fields]) '
     'AND (brain[All Fields] OR cortex[All Fields] OR SVZ[All Fields] OR neuron[All Fields]) '
     'AND ("single-cell"[All Fields] OR scRNA[All Fields] OR snRNA[All Fields] OR multiome[All Fields]) '
     'AND "Mus musculus"[Organism] AND gse[Entry Type]'),
]

BRAIN_RX = re.compile(
    r"brain|cortex|cortical|hippocamp|prefrontal|DLPFC|mPFC|SVZ|neocortex|"
    r"neuron|glia|astrocyte|oligodendrocyt",
    re.I,
)
OSK_RX = re.compile(r"OSKM|\bOSK\b|Yamanaka|partial reprogramming|Oct4.*Sox2.*Klf4", re.I)
RNA_RX = re.compile(r"10x|Seq-Well|Drop-seq|Smart-seq|sci-RNA|inDrop|CEL-seq|snRNA|scRNA|multiome", re.I)


def _soft_row(acc: str, log) -> dict:
    log(f"\n[GEO] fetching SOFT {acc}")
    try:
        text = fetch_soft(acc)
    except Exception as e:
        log(f"  fetch failed: {e}")
        return dict(accession=acc, found=False, error=str(e))
    if "!Series_title" not in text:
        log(f"  NOT FOUND in GEO (no Series record). Head:\n{text[:240]}")
        return dict(accession=acc, found=False)
    m = parse_series(text)
    title = (m.get("Series_title") or [""])[0]
    organism = "; ".join(m.get("Series_sample_organism", m.get("Series_platform_organism", ["?"])))
    n_samples = len(m.get("Series_sample_id", []))
    pubmed = "; ".join(m.get("Series_pubmed_id", ["none listed"]))
    summary = " ".join(m.get("Series_summary", []))
    design = " ".join(m.get("Series_overall_design", []))
    stype = "; ".join(m.get("Series_type", ["?"]))
    supp = list(m.get("Series_supplementary_file", []) or [])
    is_sc = bool(re.search(r"single.cell|single.nucleus|scRNA|snRNA|multiome",
                           title + summary + stype, re.I))
    brain = bool(BRAIN_RX.search(title + " " + summary + " " + design))
    osk = bool(OSK_RX.search(title + " " + summary + " " + design))
    log(f"  Title:    {title}")
    log(f"  Organism: {organism}")
    log(f"  Type:     {stype}")
    log(f"  N samples:{n_samples}  PubMed: {pubmed}")
    log(f"  brain={brain}  OSK/Yamanaka={osk}  sc/sn={is_sc}")
    log(f"  Summary:  {summary[:400]}")
    log(f"  Design:   {design[:400]}")
    log(f"  Suppl ({len(supp)}):")
    for f in supp[:12]:
        log(f"    {f}")
    return dict(
        accession=acc, found=True, title=title, organism=organism, n_samples=n_samples,
        pubmed=pubmed, series_type=stype, is_scrna=is_sc, brain=brain, osk=osk,
        summary=summary[:1200], design=design[:800],
        n_suppl=len(supp), suppl="; ".join(supp[:20]),
    )


def _geo_search(qid, query, log, retmax=40):
    log(f"\n[GEO search] {qid}")
    log(f"  QUERY: {query}")
    try:
        es = esearch(query, retmax=retmax)
    except Exception as e:
        log(f"  esearch failed: {e}")
        return []
    log(f"  GEO returned {es.get('count')} records (showing up to {len(es.get('idlist', []))})")
    time.sleep(0.35)
    rows = []
    try:
        summaries = esummary(es.get("idlist", []))
    except Exception as e:
        log(f"  esummary failed: {e}")
        return []
    for s in summaries:
        acc = s.get("accession")
        title = s.get("title", "")
        rows.append(dict(
            query_id=qid, accession=acc, n_samples=s.get("n_samples"),
            pdat=s.get("pdat"), gdstype=s.get("gdstype", ""),
            taxon=s.get("taxon", ""), title=title,
            brain=bool(BRAIN_RX.search(title or "")),
            osk=bool(OSK_RX.search(title or "")),
        ))
        log(f"  {acc:<11} n={s.get('n_samples')}  {s.get('pdat')}  | {title[:110]}")
    time.sleep(0.35)
    return rows


def _parse_mouse_age(stage: str):
    if stage is None:
        return None
    a = parse_age_years(stage)
    if a is not None:
        return a
    s = str(stage).lower()
    m = re.match(r"^(\d+(?:\.\d+)?)\s*-?\s*month", s)
    if m:
        return float(m.group(1)) / 12.0 * 12.0  # keep months as months? parse_age_years already /12 to years
    # parse_age_years already handles month-old -> years. Also:
    m = re.match(r"^(\d+)\s*-?\s*week", s)
    if m:
        return float(m.group(1)) / 52.0
    return a


def _cxg_mouse_brain(log):
    log("\n[CXG] fetching curation API /datasets (all organisms; filter Mus musculus + brain)")
    r = requests.get(f"{API}/datasets", timeout=600)
    r.raise_for_status()
    ds = r.json()
    log(f"  {len(ds)} datasets returned")
    rows = []
    for d in ds:
        orgs = [o.get("label") for o in d.get("organism", [])]
        if "Mus musculus" not in orgs:
            continue
        tissues = sorted({t.get("label") for t in d.get("tissue", [])})
        tjoin = "; ".join(tissues)
        if not BRAIN_RX.search(tjoin):
            continue
        assays = sorted({a.get("label") for a in d.get("assay", [])})
        if not RNA_RX.search("; ".join(assays)):
            continue
        diseases = sorted({x.get("label") for x in d.get("disease", [])})
        cts = sorted({c.get("label") for c in d.get("cell_type", [])})
        stages = [s.get("label") for s in d.get("development_stage", [])]
        donors = d.get("donor_id", []) or []
        ages = sorted(a for a in (_parse_mouse_age(s) for s in stages) if a is not None)
        title = d.get("title") or ""
        rows.append(dict(
            dataset_id=d.get("dataset_id"),
            collection_id=d.get("collection_id"),
            title=title,
            collection_name=d.get("collection_name"),
            cell_count=d.get("cell_count"),
            n_donors=len(donors),
            n_cell_types=len(cts),
            assays="; ".join(assays),
            tissues=tjoin[:300],
            diseases="; ".join(diseases)[:200],
            n_dev_stages=len(stages),
            n_stages_numeric=len(ages),
            age_min=ages[0] if ages else None,
            age_max=ages[-1] if ages else None,
            age_span=(ages[-1] - ages[0]) if ages else None,
            ages_numeric=";".join(f"{a:g}" for a in ages[:20]),
            dev_stages_sample="; ".join(stages[:15]),
            cell_types_sample="; ".join(cts[:20]),
            h5ad_url=next((a.get("url") for a in d.get("assets", []) if a.get("filetype") == "H5AD"), None),
            h5ad_bytes=next((a.get("filesize") for a in d.get("assets", []) if a.get("filetype") == "H5AD"), None),
            osk=bool(OSK_RX.search(title + " " + (d.get("collection_name") or ""))),
            collection_doi=d.get("collection_doi"),
        ))
    df = pd.DataFrame(rows)
    log(f"  Mus musculus + brain + RNA assay: {len(df)}")
    return df


def run():
    log = Logger(TRAJ_DIR / "g3_search.txt")
    try:
        return _run(log)
    finally:
        log.close()


def _run(log):
    traj_log_banner(log, "G3 SEARCH — reprogramming + mouse aging (no downloads yet)")
    log(f"seed={TRAJ_SEED}  (search only; no model is fit)")

    cand_rows = []
    for acc in CANDIDATE_REPROG:
        cand_rows.append(_soft_row(acc, log))
        time.sleep(0.35)
    cand = pd.DataFrame(cand_rows)
    cand.to_csv(TRAJ_DIR / "g3_search_candidates_soft.csv", index=False)

    geo_rows = []
    for qid, q in GEO_QUERIES:
        geo_rows.extend(_geo_search(qid, q, log))
    geo = pd.DataFrame(geo_rows)
    geo.to_csv(TRAJ_DIR / "g3_search_geo_hits.csv", index=False)

    # Fetch SOFT for new human-brain OSK hits not already in the candidate list
    extra_acc = []
    if len(geo):
        human_brain = geo[(geo.query_id == "human_brain_OSK_sc") | (
            (geo.query_id == "human_OSK_sc_any_tissue") & geo.brain
        )]
        known = set(CANDIDATE_REPROG)
        extra_acc = [a for a in human_brain.accession.dropna().unique() if a not in known]
        extra_acc = extra_acc[:12]
    extra_rows = []
    if extra_acc:
        log(f"\n[GEO] SOFT for additional human-brain-ish OSK hits ({len(extra_acc)})")
        for acc in extra_acc:
            extra_rows.append(_soft_row(acc, log))
            time.sleep(0.35)
    extra = pd.DataFrame(extra_rows)
    extra.to_csv(TRAJ_DIR / "g3_search_human_extra_soft.csv", index=False)

    cxg = _cxg_mouse_brain(log)
    cxg.to_csv(TRAJ_DIR / "g3_search_cxg_mouse_brain.csv", index=False)

    # Rank mouse aging candidates: numeric age span, n_donors, not OSK (independent aging)
    aging = cxg.copy()
    if len(aging):
        aging = aging[~aging.osk.astype(bool)].copy()
        aging["usable_span"] = aging.age_span.fillna(0)
        aging["usable_donors"] = aging.n_donors.fillna(0)
        aging = aging.sort_values(["usable_span", "usable_donors", "cell_count"],
                                  ascending=False)
        aging.head(25).to_csv(TRAJ_DIR / "g3_search_cxg_mouse_aging_ranked.csv", index=False)
        log("\n[CXG] top mouse brain RNA datasets by numeric age span (OSK titles excluded):")
        cols = [c for c in ("dataset_id", "title", "n_donors", "cell_count",
                            "age_min", "age_max", "age_span", "h5ad_bytes") if c in aging.columns]
        log(aging[cols].head(12).to_string(index=False))

    # Human CXG OSK (should be empty per T1; re-check)
    log("\n[CXG] human OSK/Yamanaka title scan (re-check)")
    r = requests.get(f"{API}/datasets", timeout=600)
    r.raise_for_status()
    ds = r.json()
    human_osk = []
    for d in ds:
        orgs = [o.get("label") for o in d.get("organism", [])]
        if "Homo sapiens" not in orgs:
            continue
        blob = (d.get("title") or "") + " " + (d.get("collection_name") or "")
        if OSK_RX.search(blob):
            tissues = "; ".join(sorted({t.get("label") for t in d.get("tissue", [])}))
            human_osk.append(dict(
                dataset_id=d.get("dataset_id"), title=d.get("title"),
                collection_name=d.get("collection_name"),
                tissues=tissues[:200], cell_count=d.get("cell_count"),
                n_donors=len(d.get("donor_id") or []),
                h5ad_url=next((a.get("url") for a in d.get("assets", []) if a.get("filetype") == "H5AD"), None),
                h5ad_bytes=next((a.get("filesize") for a in d.get("assets", []) if a.get("filetype") == "H5AD"), None),
                brain=bool(BRAIN_RX.search(blob + " " + tissues)),
            ))
    hdf = pd.DataFrame(human_osk)
    hdf.to_csv(TRAJ_DIR / "g3_search_cxg_human_osk.csv", index=False)
    log(f"  human OSK/Yamanaka title matches: {len(hdf)}")
    if len(hdf):
        log(hdf.to_string(index=False)[:3000])

    dump_json(TRAJ_DIR / "g3_search_summary.json", dict(
        seed=TRAJ_SEED,
        candidates=CANDIDATE_REPROG,
        n_geo_hits=int(len(geo)),
        n_human_extra_soft=int(len(extra)),
        n_cxg_mouse_brain=int(len(cxg)),
        n_cxg_human_osk=int(len(hdf)),
        candidate_found={r["accession"]: bool(r.get("found")) for r in cand_rows},
    ))
    log("\n[G3 search] wrote results/trajectory/g3_search_*")
    log("[G3 search] done.")
    return dict(candidates=cand, geo=geo, cxg=cxg)


if __name__ == "__main__":
    run()
