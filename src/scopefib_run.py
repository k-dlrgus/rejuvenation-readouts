"""SCOPEFIB: metadata-only record fetcher for the fibrosis scoping search.

No matrix download, no statistic. Every response is appended verbatim (truncated
only where noted) to a JSON log under results/scopefib/ so each table row can be
traced to the record text it was resolved from.

Usage:
  python src/scopefib_run.py search "<GDS term>" [...]      -> ncbi_searches.json
  python src/scopefib_run.py series GSE... [...]            -> geo_records.json
  python src/scopefib_run.py gsm GSM... [...]               -> gsm_records.json
  python src/scopefib_run.py url <label> <url> [...]         -> url_records.json (≤50 MB, text only)
"""
from __future__ import annotations

import json
import sys
import time
import urllib.parse
import urllib.request
from datetime import datetime, timezone
from pathlib import Path

OUT = Path(__file__).resolve().parents[1] / "results" / "scopefib"
UA = "AgeIdentitySeparabilityBenchmark/scopefib (metadata only)"
MAX_BYTES = 50 * 1024 * 1024
EUTILS = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/"


def get(url: str, timeout: int = 120) -> bytes:
    req = urllib.request.Request(url, headers={"User-Agent": UA})
    with urllib.request.urlopen(req, timeout=timeout) as r:
        n = r.headers.get("Content-Length")
        if n is not None and int(n) > MAX_BYTES:
            raise RuntimeError(f"refused: {url} is {n} bytes (>50 MB)")
        data = r.read(MAX_BYTES + 1)
        if len(data) > MAX_BYTES:
            raise RuntimeError(f"refused: {url} exceeded 50 MB while reading")
        return data


def append(name: str, record: dict) -> None:
    path = OUT / name
    log = json.loads(path.read_text(encoding="utf-8")) if path.exists() else []
    record["fetched_utc"] = datetime.now(timezone.utc).isoformat(timespec="seconds")
    log.append(record)
    path.write_text(json.dumps(log, indent=1, ensure_ascii=False), encoding="utf-8")


def esearch_esummary(term: str, retmax: int = 100) -> dict:
    q = urllib.parse.urlencode({"db": "gds", "term": term, "retmax": retmax, "retmode": "json"})
    es = json.loads(get(EUTILS + "esearch.fcgi?" + q))["esearchresult"]
    ids = es.get("idlist", [])
    hits = []
    if ids:
        time.sleep(0.4)
        q = urllib.parse.urlencode({"db": "gds", "id": ",".join(ids), "retmode": "json"})
        res = json.loads(get(EUTILS + "esummary.fcgi?" + q))["result"]
        for uid in res.get("uids", []):
            d = res[uid]
            hits.append({k: d.get(k) for k in ("accession", "entrytype", "taxon", "title", "gdstype", "n_samples", "pdat")})
    return {"term": term, "count": int(es.get("count", 0)), "hits": hits}


def geo_text(acc: str) -> str:
    q = urllib.parse.urlencode({"acc": acc, "targ": "self", "form": "text", "view": "quick"})
    return get("https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?" + q).decode("utf-8", errors="replace")


def main(argv: list[str]) -> None:
    OUT.mkdir(parents=True, exist_ok=True)
    mode, args = argv[0], argv[1:]
    if mode == "search":
        for term in args:
            r = esearch_esummary(term)
            append("ncbi_searches.json", r)
            print(f"## {term}  (count={r['count']})")
            for h in r["hits"]:
                if h["entrytype"] == "GSE":
                    print(f"  {h['accession']}\t{h['taxon']}\tn={h['n_samples']}\t{h['pdat']}\t{h['title']}")
            time.sleep(0.4)
    elif mode in ("series", "gsm"):
        for acc in args:
            txt = geo_text(acc)
            append("geo_records.json" if mode == "series" else "gsm_records.json", {"accession": acc, "text": txt})
            print(f"===== {acc}\n{txt}")
            time.sleep(0.4)
    elif mode == "url":
        label, url = args[0], args[1]
        txt = get(url).decode("utf-8", errors="replace")
        append("url_records.json", {"label": label, "url": url, "text": txt})
        print(txt)
    elif mode == "lincs":
        lincs_inventory(Path(args[0]))
    elif mode == "sigcom":
        sigcom_inventory(args)
    else:
        raise SystemExit(__doc__)


# LINCS L1000 metadata tables (GEO supplementary, all < 12 MB). Raw files go to a
# scratch directory given on the command line and are not tracked; they are
# deterministic downloads from the GEO FTP paths below.
LINCS_FILES = {
    "GSE92742": ("GSE92nnn", ["GSE92742_Broad_LINCS_cell_info.txt.gz", "GSE92742_Broad_LINCS_sig_info.txt.gz"]),
    "GSE70138": ("GSE70nnn", ["GSE70138_Broad_LINCS_cell_info_2017-04-28.txt.gz", "GSE70138_Broad_LINCS_sig_info_2017-03-06.txt.gz"]),
}
FIBRO_PAT = r"fibro|stellate|myofib|IMR|WI-?38|BJ|MRC|HFL|NHDF|HDF|LX-?2"


def lincs_inventory(raw: Path) -> None:
    """Inventory of LINCS L1000 fibroblast arms: per cell line, signatures by pert_type,
    distinct perturbagens, and whether any TGF-beta ligand arm exists. Counts only."""
    import csv
    import gzip
    import re

    raw.mkdir(parents=True, exist_ok=True)
    rows_cells, rows_arms, rows_tgf = [], [], []
    for gse, (stub, files) in LINCS_FILES.items():
        tabs = {}
        for f in files:
            p = raw / f
            if not p.exists():
                p.write_bytes(get(f"https://ftp.ncbi.nlm.nih.gov/geo/series/{stub}/{gse}/suppl/{f}"))
            with gzip.open(p, "rt", encoding="utf-8", errors="replace") as fh:
                tabs["cell" if "cell_info" in f else "sig"] = list(csv.DictReader(fh, delimiter="\t"))
        cells = tabs["cell"]
        fib = {}
        for c in cells:
            text = " ".join(str(v) for v in c.values())
            if re.search(FIBRO_PAT, text, re.I):
                fib[c.get("cell_id") or c.get("cell_iname")] = c
                rows_cells.append({"gse": gse, **{k: c.get(k, "") for k in ("cell_id", "cell_type", "cell_lineage", "cell_histology", "cell_source_id", "primary_site", "subtype", "original_growth_pattern", "donor_age", "donor_sex", "donor_ethnicity", "provider_catalog_id")}})
        per = {}
        for s in tabs["sig"]:
            cid = s.get("cell_id")
            name = s.get("pert_iname", "")
            if re.search(r"TGFB|TGF-?beta|TGFbeta", name, re.I):
                rows_tgf.append({"gse": gse, "cell_id": cid, "pert_type": s.get("pert_type"), "pert_iname": name, "pert_dose": s.get("pert_idose", s.get("pert_dose", "")), "pert_time": s.get("pert_itime", s.get("pert_time", "")), "sig_id": s.get("sig_id")})
            if cid not in fib:
                continue
            d = per.setdefault((cid, s.get("pert_type")), {"n_sig": 0, "perts": set()})
            d["n_sig"] += 1
            d["perts"].add(name)
        for (cid, pt), d in sorted(per.items()):
            rows_arms.append({"gse": gse, "cell_id": cid, "pert_type": pt, "n_signatures": d["n_sig"], "n_distinct_pert_iname": len(d["perts"])})
    for name, rows in (("lincs_fibroblast_cells.csv", rows_cells), ("lincs_fibroblast_arms.csv", rows_arms), ("lincs_tgfb_signatures.csv", rows_tgf)):
        with open(OUT / name, "w", newline="", encoding="utf-8") as fh:
            if rows:
                w = csv.DictWriter(fh, fieldnames=list(rows[0].keys()))
                w.writeheader()
                w.writerows(rows)
        print(f"wrote {name}: {len(rows)} rows")


def sigcom_inventory(cells: list[str]) -> None:
    """LINCS 2020 signatures per fibroblast line via the SigCom LINCS metadata API
    (siginfo_beta.txt is 465 MB, over the 50 MB limit). Writes one row per signature's
    metadata (no expression values) and prints per-line counts."""
    import csv

    url = "https://maayanlab.cloud/sigcom-lincs/metadata-api/signatures/find"
    keys = ("cell_line", "tissue", "cmap_id", "local_id", "pert_name", "pert_type", "pert_dose", "pert_time", "data_level")
    rows = []
    for c in cells:
        skip = 0
        while True:
            body = json.dumps({"filter": {"where": {"meta.cell_line": c}, "limit": 1000, "skip": skip}}).encode()
            req = urllib.request.Request(url, data=body, headers={"User-Agent": UA, "Content-Type": "application/json"})
            with urllib.request.urlopen(req, timeout=180) as r:
                page = json.loads(r.read())
            rows += [{k: s["meta"].get(k, "") for k in keys} for s in page]
            if len(page) < 1000:
                break
            skip += 1000
            time.sleep(0.5)
    with open(OUT / "sigcom_fibroblast_signatures.csv", "w", newline="", encoding="utf-8") as fh:
        w = csv.DictWriter(fh, fieldnames=list(keys))
        w.writeheader()
        w.writerows(rows)
    for c in cells:
        sub = [r for r in rows if r["cell_line"] == c]
        plates = sorted({r["cmap_id"].split("_")[0] for r in sub})
        print(c, "n_sig", len(sub), "distinct pert_name", len({r["pert_name"] for r in sub}),
              "pert_types", sorted({r["pert_type"] for r in sub}), "plates", plates,
              "times", sorted({r["pert_time"] for r in sub}),
              "TGF-matching", sorted({r["pert_name"] for r in sub if "tgf" in r["pert_name"].lower()}))


if __name__ == "__main__":
    main(sys.argv[1:])
