"""Parse already-fetched GEO series text; fetch GSM characteristics for key series."""
from __future__ import annotations

import json
import re
import time
import urllib.parse
import urllib.request
from collections import defaultdict
from pathlib import Path

OUT = Path(__file__).resolve().parent
UA = "AgeIdentitySeparabilityBenchmark/scoping (metadata only)"


def get(url: str, timeout: int = 90) -> bytes:
    req = urllib.request.Request(url, headers={"User-Agent": UA})
    with urllib.request.urlopen(req, timeout=timeout) as r:
        return r.read()


def geo_text(acc: str, view: str = "full") -> str:
    url = (
        "https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?"
        + urllib.parse.urlencode({"acc": acc, "targ": "self", "form": "text", "view": view})
    )
    return get(url).decode("utf-8", errors="replace")


def parse_series(txt: str) -> dict:
    d = defaultdict(list)
    for line in txt.splitlines():
        if line.startswith("!") and " = " in line:
            k, v = line.split(" = ", 1)
            d[k].append(v)
    return dict(d)


def parse_gsm(txt: str) -> dict:
    fields = {}
    chars = []
    title = None
    source = None
    organism = None
    for line in txt.splitlines():
        if line.startswith("!Sample_title = "):
            title = line.split(" = ", 1)[1]
        elif line.startswith("!Sample_source_name_ch1 = "):
            source = line.split(" = ", 1)[1]
        elif line.startswith("!Sample_organism_ch1 = "):
            organism = line.split(" = ", 1)[1]
        elif line.startswith("!Sample_characteristics_ch1 = "):
            chars.append(line.split(" = ", 1)[1])
        elif line.startswith("!Sample_data_processing") or line.startswith("!Sample_description"):
            pass
    return {
        "title": title,
        "source": source,
        "organism": organism,
        "characteristics": chars,
    }


# Load whatever records exist
rec_path = OUT / "geo_records.json"
records = json.loads(rec_path.read_text(encoding="utf-8")) if rec_path.exists() else {}

parsed = {}
for acc, rec in records.items():
    txt = rec.get("geo_text_quick") or ""
    if not txt:
        continue
    s = parse_series(txt)
    parsed[acc] = {
        "title": (s.get("!Series_title") or [None])[0],
        "summary": (s.get("!Series_summary") or [None])[0],
        "overall_design": (s.get("!Series_overall_design") or [None])[0],
        "type": s.get("!Series_type"),
        "sample_organism": s.get("!Series_sample_organism"),
        "n_sample_ids": len(s.get("!Series_sample_id") or []),
        "sample_ids": s.get("!Series_sample_id") or [],
        "platform_id": s.get("!Series_platform_id"),
        "relation": s.get("!Series_relation"),
        "supplementary_file": s.get("!Series_supplementary_file"),
        "n_samples_field": (s.get("!Series_sample_id") and len(s["!Series_sample_id"])),
        "pubmed_id": s.get("!Series_pubmed_id"),
        "web_link": s.get("!Series_web_link"),
    }

(OUT / "series_parsed.json").write_text(json.dumps(parsed, indent=2), encoding="utf-8")
print("parsed series", list(parsed))

# GSM fetch for small series that could be (B)
gsm_series = {
    "GSE297234": parsed.get("GSE297234", {}).get("sample_ids") or [
        f"GSM{n}" for n in range(8986586, 8986594)
    ],
    "GSE297233": parsed.get("GSE297233", {}).get("sample_ids") or [
        f"GSM{n}" for n in range(8986578, 8986586)
    ],
}

# Gill RNA-seq series: fetch ALL GSMs if n is manageable
for acc in ("GSE165176", "GSE165177", "GSE165178", "GSE165179"):
    ids = parsed.get(acc, {}).get("sample_ids") or []
    if ids and len(ids) <= 80:
        gsm_series[acc] = ids

gsm_out = {}
for series, ids in gsm_series.items():
    print(f"GSM {series} n={len(ids)}", flush=True)
    rows = []
    for i, gsm in enumerate(ids):
        try:
            txt = geo_text(gsm, "full")
            rows.append({"accession": gsm, **parse_gsm(txt)})
            print(" ", gsm, rows[-1]["title"], rows[-1]["characteristics"][:6], flush=True)
        except Exception as e:
            rows.append({"accession": gsm, "error": f"{type(e).__name__}: {e}"})
            print("  ERR", gsm, e, flush=True)
            time.sleep(1)
        time.sleep(0.34)
    gsm_out[series] = rows
    (OUT / "gsm_tables.json").write_text(json.dumps(gsm_out, indent=2), encoding="utf-8")

print("wrote gsm_tables.json")
