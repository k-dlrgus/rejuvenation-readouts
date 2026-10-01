"""Print titles and a few series fields from pages already saved this session."""
import json
from pathlib import Path

root = Path("results/survey/opened")
for name in [
    "bioproject_PRJNA1108254.json",
    "bioproject_PRJNA831566.json",
    "figshare_30898748.json",
]:
    data = json.loads((root / name).read_text(encoding="utf-8"))
    print("FILE", name)
    if isinstance(data, dict):
        print(" keys", list(data)[:12])
        for k in ("title", "description", "doi", "url", "figshare_url"):
            if k in data and not isinstance(data[k], (dict, list)):
                text = str(data[k]).replace("\n", " ")
                print(k, text[:500])
    print("---")

for gse in [
    "GSE92742",
    "GSE70138",
    "GSE216481",
    "GSE216595",
    "GSE159786",
    "GSE159779",
    "GSE159778",
    "GSE133344",
    "GSE107185",
    "GSE80676",
    "GSE33816",
    "GSE297234",
    "GSE165180",
    "GSE165177",
    "GSE297984",
    "GSE226189",
]:
    p = root / f"{gse}_quick.txt"
    if not p.exists():
        print("MISSING", gse)
        continue
    text = p.read_text(encoding="utf-8", errors="replace")
    want = ("!Series_title", "!Series_pubmed_id", "!Series_summary", "!Series_overall_design", "!Series_type")
    print("FILE", gse)
    for line in text.splitlines():
        if line.startswith(want):
            print(line[:400])
    print("---")
