import json
from pathlib import Path
d = json.loads(Path("results/scoping/gse325735_parse.json").read_text(encoding="utf-8"))
tg = d["GSE325735"]["transgenes"]
skip = {"WT", "mTagBFP2", "NTg"}
keep = []
unp = []
for k, v in tg.items():
    if k.startswith("UNPARSED"):
        unp.append((k, v))
    elif k not in skip:
        keep.append((k, v))
print("n keep", len(keep), "n skip", 3, "n unparsed keys", len(unp))
print("keep names", sorted(k for k,_ in keep))
print("unparsed", unp)
print("OSKM" in dict(keep))
