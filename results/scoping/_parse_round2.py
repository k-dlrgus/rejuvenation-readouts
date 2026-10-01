import json
import re
from collections import Counter, defaultdict
from pathlib import Path

d = json.loads(Path("results/scoping/geo_records_round2.json").read_text(encoding="utf-8"))

for acc, rec in d.items():
    es = rec.get("esummary") or {}
    print("=" * 80)
    print(acc)
    print(" title:", es.get("title") or rec.get("title"))
    print(" n:", es.get("n_samples"), "type:", es.get("gdstype"), "taxon:", es.get("taxon"))
    print(" design:", (rec.get("overall_design") or "")[:500])
    print(" summary:", ((es.get("summary") or rec.get("summary") or "")[:350]))
    titles = [s.get("title") or "" for s in (es.get("sample_titles") or [])]
    print(" n titles returned", len(titles))
    # print unique-ish prefixes
    for t in titles[:15]:
        print("  ", t)
    if len(titles) > 15:
        print("  ...", len(titles) - 15, "more")

# Special parse GSE325735
print("\n\n##### GSE325735 TITLE PARSE #####")
titles = [s.get("title") or "" for s in (d["GSE325735"]["esummary"].get("sample_titles") or [])]
print("n titles", len(titles))
# dump all unique titles sorted for inspection
Path("results/scoping/gse325735_titles.txt").write_text("\n".join(sorted(titles)), encoding="utf-8")

# guess tokens
c = Counter(titles)
print("unique titles", len(c), "dup max", c.most_common(3))

# try split on common separators
patterns = defaultdict(int)
for t in titles:
    # collapse numbers of replicates
    t2 = re.sub(r"rep(licate)?\s*\d+", "repN", t, flags=re.I)
    t2 = re.sub(r"_r\d+\b", "_rN", t2, flags=re.I)
    patterns[t2] += 1
print("collapsed unique", len(patterns))
for k, v in list(sorted(patterns.items(), key=lambda kv: -kv[1]))[:40]:
    print(f"  {v:4d} {k}")

# donor-like tokens
donors = Counter()
genes = Counter()
ages = Counter()
for t in titles:
    for m in re.findall(r"\bM\d{2}\b", t):
        donors[m] += 1
    for m in re.findall(r"\b(\d{2,3})\s*(y|yr|year)", t, flags=re.I):
        ages[m[0]] += 1
    for m in re.findall(r"\b(NHDF|HDF|iPSC|IMR90|BJ|GM\d+)\b", t, flags=re.I):
        donors[m] += 1
print("donor-like", donors.most_common(30))
print("ages-like", ages.most_common(20))

# token after OE/overexpression/cDNA
oe = Counter()
for t in titles:
    m = re.search(r"(?:OE_|overexpression_|cDNA_|tg_|transgene[:= ]+)([A-Za-z0-9-]+)", t, flags=re.I)
    if m:
        oe[m.group(1)] += 1
print("OE-like", len(oe), oe.most_common(20))
