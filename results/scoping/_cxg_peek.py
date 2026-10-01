import json
import sys
from pathlib import Path
sys.stdout.reconfigure(encoding="utf-8", errors="replace")
d = json.loads(Path("results/scoping/round3.json").read_text(encoding="utf-8"))
print("cxg total", d.get("cxg_n_datasets_total"), "skinish", d.get("cxg_n_skinish"))
hits = d.get("cxg_skin") or []
print("stored hits", len(hits), "keys example", hits[0].get("keys") if hits else None)
# sort by n_cells
hits_sorted = sorted(hits, key=lambda x: (x.get("n_cells") or 0), reverse=True)
for h in hits_sorted[:15]:
    print(h.get("n_cells"), (h.get("title") or "")[:110], h.get("dataset_id"))
    print("  tissues", h.get("tissues")[:8] if h.get("tissues") else None)

# Sarkar sra titles from sra_meta
s = json.loads(Path("results/scoping/sra_meta.json").read_text(encoding="utf-8"))
print("\nSarkar PRJNA598923 count", s["PRJNA598923"].get("sra_count"))
for row in s["PRJNA598923"].get("sra_rows") or []:
    xml = row.get("expxml_head") or ""
    # extract Title
    i = xml.find("<Title>")
    j = xml.find("</Title>")
    title = xml[i+7:j] if i>=0 and j>i else None
    print(" ", title)
print("Southard error", s.get("PRJNA1108254", {}).get("sra_error") or s.get("PRJNA1108254", {}).get("bioproject_error"))
