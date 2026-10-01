"""Extract donor-age fields from GEO sample SOFT already saved."""
import gzip
import re
from collections import Counter
from pathlib import Path

text = Path("results/survey/opened/GSE165177_gsm.txt").read_text(encoding="utf-8", errors="replace")
# split samples
blocks = text.split("^SAMPLE")
print("BLOCKS", len(blocks) - 1)
ages = []
titles = []
for b in blocks[1:]:
    title = ""
    m = re.search(r"!Sample_title = (.+)", b)
    if m:
        title = m.group(1).strip()
        titles.append(title)
    chars = re.findall(r"!Sample_characteristics_ch1 = (.+)", b)
    age = [c for c in chars if "age" in c.lower() or "donor" in c.lower()]
    if age:
        ages.append((title, age))
print("TITLES", len(titles))
print("WITH_AGE_OR_DONOR", len(ages))
for t, a in ages[:8]:
    print("EX", t, "|", a)
# unique age values
vals = []
for t, a in ages:
    for item in a:
        if "age" in item.lower():
            vals.append(item)
print("AGE_FIELD_COUNTS")
for k, v in Counter(vals).most_common():
    print(v, k)
pairs = []
for t, a in ages:
    age = next((x for x in a if "age" in x.lower()), "")
    prefix = t.split("_")[0]
    pairs.append((prefix, age))
print("PREFIX_AGE")
for k, v in Counter(pairs).most_common():
    print(v, k)

lam = set(
    x.strip()
    for x in Path("data/reference/TF_names_v_1.01.txt").read_text(encoding="utf-8").splitlines()
    if x.strip()
)
pheno = Path("results/survey/opened/GSE107185_up-tf-stem_pheno_dict.csv.gz")
lines = gzip.open(pheno, "rt", encoding="utf-8", errors="replace").read().splitlines()
keys = []
for ln in lines:
    if not ln.strip():
        continue
    keys.append(ln.split(",")[0].strip().strip('"'))
print("PHENO_ROWS", len(keys))
print("PHENO_KEYS", keys)
print("PHENO_LAMBERT", len(set(keys) & lam))
