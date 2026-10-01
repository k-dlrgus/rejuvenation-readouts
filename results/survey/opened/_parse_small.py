"""One-off parsers for survey metadata already downloaded. Not part of the analysis."""
import gzip
from collections import Counter
from pathlib import Path

lam = set(
    x.strip()
    for x in Path("data/reference/TF_names_v_1.01.txt").read_text(encoding="utf-8").splitlines()
    if x.strip()
)

print("LAMBERT", len(lam))
for name in ["SPI1", "MNDA", "CEBPA", "IRF8", "PRRX1", "OSR1", "LHX9", "TWIST2", "SRF", "EGR1", "ATF3", "TCF7L2"]:
    print("GENE", name, name in lam)

p = Path("results/survey/opened/GSE133344_filtered_cell_identities.csv.gz")
genes = set()
n = 0
with gzip.open(p, "rt", encoding="utf-8", errors="replace") as f:
    header = f.readline()
    for line in f:
        n += 1
        ident = line.split(",")[1]
        # "GENE_guide__GENE_guide" or single
        for half in ident.split("__"):
            sym = half.split("_")[0]
            if sym and sym not in ("NegCtrl", "neg", "NT"):
                genes.add(sym)
print("NORMAN_ROWS", n, "SYMBOLS", len(genes), "LAMBERT", len(genes & lam))
print("NON_LAMBERT", sorted(genes - lam)[:30], "n", len(genes - lam))
