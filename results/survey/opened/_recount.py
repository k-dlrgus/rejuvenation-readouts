"""Recount Lambert overlaps from small files already opened."""
import csv
import gzip
from pathlib import Path

lambert = set(Path("data/reference/TF_names_v_1.01.txt").read_text(encoding="utf-8").split())
out = []

def nlines(path):
    opener = gzip.open if str(path).endswith(".gz") else open
    with opener(path, "rt", encoding="utf-8", errors="replace") as fh:
        return sum(1 for _ in fh)

for p in [
    Path("results/survey/gene_lists/GSE216595_genes.tsv.gz"),
    Path("results/survey/gene_lists/GSE133344_filtered_genes.tsv.gz"),
]:
    if p.exists():
        out.append(f"{p.name} lines {nlines(p)}")
    else:
        out.append(f"MISSING {p}")

# Joung TFmap
p = Path("results/survey/opened/joung/GSM6681047_180124_TFmap.csv.gz")
if p.exists():
    genes = set()
    with gzip.open(p, "rt", encoding="utf-8", errors="replace") as fh:
        for i, line in enumerate(fh):
            if i == 0:
                out.append("tfmap header " + line.strip()[:200])
                continue
            parts = line.strip().split(",")
            if len(parts) >= 2:
                genes.add(parts[1].strip())
    named = {g for g in genes if g and g not in {"NA", "AMB", "nan", "None"}}
    hit = sorted(g for g in named if g in lambert)
    miss = sorted(g for g in named if g not in lambert)
    out.append(f"joung unique col2 {len(genes)} named {len(named)} lambert {len(hit)} not {len(miss)}")
    out.append("joung not " + ",".join(miss))

p = Path("results/survey/opened/GSE133344_filtered_cell_identities.csv.gz")
if p.exists():
    syms = set()
    with gzip.open(p, "rt", encoding="utf-8", errors="replace") as fh:
        header = fh.readline()
        out.append("norman header " + header.strip()[:200])
        for line in fh:
            ident = line.split(",")[1] if "," in line else ""
            for half in ident.split("__"):
                tok = half.split("_")[0]
                if tok and tok not in {"NegCtrl0", "NegCtrl1", "NegCtrl"}:
                    syms.add(tok)
    hit = {s for s in syms if s in lambert}
    out.append(f"norman symbols {len(syms)} lambert {len(hit)} not {sorted(syms-hit)}")

p = Path("results/survey/opened/GSE107185_up-tf-stem_pheno_dict.csv.gz")
if p.exists():
    keys = []
    with gzip.open(p, "rt", encoding="utf-8", errors="replace") as fh:
        for line in fh:
            k = line.split(",")[0].strip().strip('"')
            if k and k.lower() != "index":
                keys.append(k)
    hit = [k for k in keys if k in lambert]
    out.append(f"parekh keys {len(keys)} lambert {len(hit)} not {[k for k in keys if k not in lambert]}")

Path("results/survey/opened/_recount.txt").write_text("\n".join(out), encoding="utf-8")
print("wrote", len(out))
