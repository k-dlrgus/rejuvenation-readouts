import json
import re
import sys
from collections import Counter, defaultdict
from pathlib import Path

sys.stdout.reconfigure(encoding="utf-8", errors="replace")

d = json.loads(Path("results/scoping/geo_records_round2.json").read_text(encoding="utf-8"))
titles = [s.get("title") or "" for s in (d["GSE325735"]["esummary"].get("sample_titles") or [])]
print("n titles", len(titles))

human = []
mouse = []
other = []
for t in titles:
    if re.search(r"\bMouse#", t, re.I) or t.startswith("Brain from") or t.startswith("Liver from"):
        mouse.append(t)
    elif re.search(r"\bM\d{2}\b", t) or "iPSC" in t or "NHDF" in t or "fibroblast" in t.lower() or "WT line" in t:
        human.append(t)
    else:
        other.append(t)

print("human-like", len(human), "mouse-like", len(mouse), "other", len(other))
print("other examples:")
for t in other[:30]:
    print(" ", t)

donors = Counter()
transgenes = Counter()
wt = 0
ipsc = 0
for t in human:
    ms = re.findall(r"\bM(\d{2})\b", t)
    for m in ms:
        donors[int(m)] += 1
    if "iPSC" in t:
        ipsc += 1
    if re.search(r"\bWT line\b", t):
        wt += 1
        transgenes["WT"] += 1
        continue
    m = re.search(r"Dox-inducible ([A-Za-z0-9]+) transgene", t)
    if m:
        transgenes[m.group(1)] += 1
    else:
        m2 = re.search(r"overexpressing ([A-Za-z0-9]+)", t, re.I)
        if m2:
            transgenes["OE:" + m2.group(1)] += 1
        else:
            transgenes["UNPARSED:" + t[:80]] += 1

print("\nDonors Mxx (interpreted as age codes from titles):")
for k, v in sorted(donors.items()):
    print(f"  M{k} n={v}  ge60={k>=60}")
print("n distinct M-codes", len(donors), "any>=60", any(k >= 60 for k in donors))
print("WT samples", wt, "iPSC-mentioned samples", ipsc)
print("n distinct transgene tokens", len(transgenes))
print("transgenes:")
for k, v in sorted(transgenes.items(), key=lambda kv: (-kv[1], kv[0])):
    print(f"  {v:4d} {k}")

# GSE297984 ages and timepoints
print("\n##### GSE297984 #####")
t2 = [s.get("title") or "" for s in (d["GSE297984"]["esummary"].get("sample_titles") or [])]
print("n", len(t2))
ages = Counter()
days = Counter()
rx = Counter()
for t in t2:
    print(" ", t)
    m = re.search(r"_D(\d+)_(\d+)y_([A-Za-z0-9]+)", t)
    if m:
        days[m.group(1)] += 1
        ages[m.group(2)] += 1
        rx[m.group(3)] += 1
print("days", dict(days), "ages", dict(ages), "rx", dict(rx))

# remaining series short
for acc in d:
    if acc in ("GSE325735", "GSE297984"):
        continue
    rec = d[acc]
    es = rec.get("esummary") or {}
    titles_a = [s.get("title") or "" for s in (es.get("sample_titles") or [])]
    print(f"\n==== {acc} n={es.get('n_samples')} taxon={es.get('taxon')}")
    print(" title:", es.get("title") or rec.get("title"))
    print(" design:", (rec.get("overall_design") or "")[:400])
    print(" n_titles", len(titles_a))
    for t in titles_a[:12]:
        print("  ", t)

out = {
    "GSE325735": {
        "n_titles": len(titles),
        "n_human_like": len(human),
        "n_mouse_like": len(mouse),
        "n_other": len(other),
        "donor_M_codes": {str(k): v for k, v in sorted(donors.items())},
        "n_distinct_transgenes_incl_WT_unparsed": len(transgenes),
        "transgenes": dict(sorted(transgenes.items(), key=lambda kv: (-kv[1], kv[0]))),
        "other_titles": other,
        "mouse_n": len(mouse),
    }
}
Path("results/scoping/gse325735_parse.json").write_text(json.dumps(out, indent=2), encoding="utf-8")
print("\nwrote gse325735_parse.json")
