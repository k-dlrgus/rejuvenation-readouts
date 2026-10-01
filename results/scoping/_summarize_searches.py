import json
from pathlib import Path
p = Path("results/scoping/ncbi_searches.json")
d = json.loads(p.read_text(encoding="utf-8"))
print("KEYS", list(d))
for k, v in d.items():
    print("\n====", k, "count", v.get("count"), "err", v.get("error"))
    for s in (v.get("summaries") or [])[:30]:
        acc = s.get("accession") or ""
        if not str(acc).startswith("GSE"):
            continue
        title = (s.get("title") or "")[:120]
        print(f"  {acc:12} n={s.get('n_samples')}  {s.get('gds_type')}  {title}")
