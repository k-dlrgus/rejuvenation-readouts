"""Fetch and print GEO series metadata (SOFT format) without downloading data.

Usage: python src/geo_meta.py GSE216481 [GSE...]
Prints title, summary, organism, platform, sample count, and supplementary files.
Never guesses accessions: if GEO returns nothing, it says so.
"""
import sys
import re
import requests

GEO_SOFT = "https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc={acc}&targ=self&form=text&view=brief"
GEO_FULL = "https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc={acc}&targ=gsm&form=text&view=brief"


def fetch_soft(acc: str, targ_gsm: bool = False) -> str:
    url = (GEO_FULL if targ_gsm else GEO_SOFT).format(acc=acc)
    r = requests.get(url, timeout=60)
    r.raise_for_status()
    return r.text


def parse_series(text: str) -> dict:
    out = {}
    for line in text.splitlines():
        if not line.startswith("!Series_"):
            continue
        k, _, v = line[1:].partition(" = ")
        out.setdefault(k, []).append(v.strip())
    return out


def summarize(acc: str) -> dict:
    text = fetch_soft(acc)
    if "!Series_title" not in text:
        print(f"[{acc}] NOT FOUND in GEO (no Series record). Raw response head:\n{text[:300]}")
        return {}
    meta = parse_series(text)
    print("=" * 80)
    print(f"ACCESSION: {acc}")
    print("Title:      ", meta.get("Series_title", [""])[0])
    print("Organism:   ", "; ".join(meta.get("Series_sample_organism", meta.get("Series_platform_organism", ["?"]))))
    print("Type:       ", "; ".join(meta.get("Series_type", ["?"])))
    print("Platforms:  ", "; ".join(meta.get("Series_platform_id", ["?"])))
    print("Status:     ", meta.get("Series_status", ["?"])[0])
    print("Submission: ", meta.get("Series_submission_date", ["?"])[0])
    print("PubMed:     ", "; ".join(meta.get("Series_pubmed_id", ["none listed"])))
    print("N samples:  ", len(meta.get("Series_sample_id", [])))
    print("Contributors:", ", ".join(meta.get("Series_contributor", [])[:6]), "...")
    summ = " ".join(meta.get("Series_summary", []))
    print("Summary:    ", summ[:1200] + ("..." if len(summ) > 1200 else ""))
    design = " ".join(meta.get("Series_overall_design", []))
    print("Design:     ", design[:800] + ("..." if len(design) > 800 else ""))
    print("Supplementary files:")
    for f in meta.get("Series_supplementary_file", []):
        print("   ", f)
    print("=" * 80)
    return meta


if __name__ == "__main__":
    accs = sys.argv[1:] or ["GSE216481"]
    for a in accs:
        summarize(a)
