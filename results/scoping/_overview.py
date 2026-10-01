import json
from pathlib import Path
from collections import defaultdict

raw = json.loads(Path("results/scoping/geo_records.json").read_text(encoding="utf-8"))

def parse_series(txt: str):
    d = defaultdict(list)
    for line in txt.splitlines():
        if line.startswith("!") and " = " in line:
            k, v = line.split(" = ", 1)
            d[k].append(v)
    return dict(d)

out = {}
for acc, rec in raw.items():
    es = rec.get("esummary") or {}
    txt = rec.get("geo_text_quick") or ""
    s = parse_series(txt) if txt else {}
    samples = es.get("samples") or []
    out[acc] = {
        "title": es.get("title") or (s.get("!Series_title") or [None])[0],
        "n_samples": es.get("n_samples"),
        "gdstype": es.get("gdstype"),
        "taxon": es.get("taxon"),
        "pdat": es.get("pdat"),
        "summary": (es.get("summary") or "")[:500],
        "overall_design": (s.get("!Series_overall_design") or [None])[0],
        "type": s.get("!Series_type"),
        "relation": s.get("!Series_relation"),
        "pubmed": s.get("!Series_pubmed_id"),
        "sample_titles": [{"accession": x.get("accession"), "title": x.get("title")} for x in samples],
        "n_sample_titles": len(samples),
        "esummary_error": rec.get("esummary_error"),
        "geo_text_error": rec.get("geo_text_error"),
        "geo_text_nchars": rec.get("geo_text_nchars"),
        "text_sample_ids": s.get("!Series_sample_id"),
    }
    print("=" * 80)
    print(acc, out[acc]["title"])
    print(" n", out[acc]["n_samples"], out[acc]["gdstype"], out[acc]["taxon"])
    print(" design:", (out[acc]["overall_design"] or "")[:400])
    print(" n_sample_titles", out[acc]["n_sample_titles"])
    for st in out[acc]["sample_titles"][:25]:
        print("  ", st["accession"], st["title"])
    if out[acc]["n_sample_titles"] > 25:
        print("  ...")

Path("results/scoping/series_overview.json").write_text(json.dumps(out, indent=2), encoding="utf-8")
print("\nwrote series_overview.json", len(out))
