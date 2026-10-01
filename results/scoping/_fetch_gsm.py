"""Fetch GSM characteristics for GSE325735 (a few) plus GSE297234 all 8 and Gill/Joung."""
from __future__ import annotations

import json
import time
import urllib.parse
import urllib.request
from collections import defaultdict
from pathlib import Path

OUT = Path(__file__).resolve().parent
UA = "AgeIdentitySeparabilityBenchmark/scoping (metadata only)"

round2 = json.loads((OUT / "geo_records_round2.json").read_text(encoding="utf-8"))
samples_325 = round2["GSE325735"]["esummary"]["sample_titles"]

# pick representative GSMs: one M55, M65, M79, M94, one mouse, one iPSC if present
picked = []
want = {"M55": None, "M65": None, "M79": None, "M94": None, "M68": None, "Mouse": None, "iPSC": None}
for s in samples_325:
    t = s.get("title") or ""
    acc = s.get("accession")
    for key in list(want):
        if want[key] is None and key in t:
            want[key] = acc
picked = [v for v in want.values() if v]

# GSE297234
picked += [
    "GSM8986586", "GSM8986587", "GSM8986588", "GSM8986589",
    "GSM8986590", "GSM8986591", "GSM8986592", "GSM8986593",
]
# Gill
picked += ["GSM5391793", "GSM5391824", "GSM5391847"]
# GSE142439
picked += ["GSM4227385"]
# GSE307377 all 9 if we have them
if "GSE307377" in round2:
    picked += [s["accession"] for s in round2["GSE307377"]["esummary"].get("sample_titles") or []]
# GSE226189 one
picked += ["GSM7067633"]

# unique preserve order
seen = set()
uniq = []
for x in picked:
    if x and x not in seen:
        seen.add(x)
        uniq.append(x)


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


def parse_kv(txt: str) -> dict:
    d = defaultdict(list)
    for line in txt.splitlines():
        if line.startswith("!") and " = " in line:
            k, v = line.split(" = ", 1)
            d[k].append(v)
    return dict(d)


print("n GSM to fetch", len(uniq), uniq, flush=True)
out = {}
for gsm in uniq:
    print("GSM", gsm, flush=True)
    try:
        txt = geo_text(gsm, "full")
        s = parse_kv(txt)
        out[gsm] = {
            "title": (s.get("!Sample_title") or [None])[0],
            "source": (s.get("!Sample_source_name_ch1") or [None])[0],
            "organism": (s.get("!Sample_organism_ch1") or [None])[0],
            "characteristics": s.get("!Sample_characteristics_ch1") or [],
            "description": s.get("!Sample_description") or [],
            "molecule": s.get("!Sample_molecule_ch1"),
            "extract_protocol": (s.get("!Sample_extract_protocol_ch1") or [None])[0],
            "growth_protocol": (s.get("!Sample_growth_protocol_ch1") or [None])[0],
            "series_id": s.get("!Sample_series_id"),
        }
        print(" ", out[gsm]["title"], out[gsm]["characteristics"][:10], flush=True)
        time.sleep(0.35)
    except Exception as e:
        out[gsm] = {"error": f"{type(e).__name__}: {e}"}
        print("  ERR", e, flush=True)
        time.sleep(1)

(OUT / "gsm_key.json").write_text(json.dumps(out, indent=2), encoding="utf-8")
print("wrote gsm_key.json")

# Joung relations
print("JOUNG", flush=True)
txt = geo_text("GSE216481", "quick")
s = parse_kv(txt)
rels = s.get("!Series_relation") or []
(OUT / "joung_relations.json").write_text(
    json.dumps(
        {
            "title": (s.get("!Series_title") or [None])[0],
            "summary": (s.get("!Series_summary") or [None])[0],
            "overall_design": (s.get("!Series_overall_design") or [None])[0],
            "relations": rels,
            "pubmed": s.get("!Series_pubmed_id"),
        },
        indent=2,
    ),
    encoding="utf-8",
)
print("n relations", len(rels))
for r in rels:
    print(" ", r)
