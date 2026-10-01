"""Donor ages from GEO sample pages already saved."""
from pathlib import Path

for name in ["GSE165177_gsm.txt", "GSE297984_gsm.txt", "GSE297234_quick.txt"]:
    text = Path("results/survey/opened") .joinpath(name).read_text(encoding="utf-8", errors="replace")
    print("FILE", name, "chars", len(text))
    hits = []
    for line in text.splitlines():
        low = line.lower()
        if "age" in low or "donor" in low or "characteristics" in low or "source_name" in low or "title =" in low:
            if line.startswith("!") or line.startswith("^"):
                hits.append(line[:220])
    print("nhits", len(hits))
    for h in hits[:40]:
        print(h)
    print("---")
