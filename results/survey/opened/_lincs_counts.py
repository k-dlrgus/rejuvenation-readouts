"""Counts for the LINCS rows, from the metadata files opened this session."""
import csv
import gzip
from collections import Counter, defaultdict
from pathlib import Path

root = Path("results/survey/opened/lincs")
lambert = set(Path("data/reference/TF_names_v_1.01.txt").read_text(encoding="utf-8").split())

lm = 0
inferred = 0
with gzip.open(root / "GSE92742_Broad_LINCS_gene_info.txt.gz", "rt", encoding="utf-8", errors="replace") as fh:
    for rec in csv.DictReader(fh, delimiter="\t"):
        if rec["pr_is_lm"] == "1":
            lm += 1
        else:
            inferred += 1
print("genes", lm + inferred, "landmark", lm, "not_landmark", inferred)

def cells(path):
    with gzip.open(path, "rt", encoding="utf-8", errors="replace") as fh:
        return list(csv.DictReader(fh, delimiter="\t"))

c1 = cells(root / "GSE92742_Broad_LINCS_cell_info.txt.gz")
c2 = cells(root / "GSE70138_Broad_LINCS_cell_info_2017-04-28.txt.gz")
print("n_cells_phase1", len(c1), "phase2", len(c2))
print("same_ids", {r["cell_id"] for r in c1} == {r["cell_id"] for r in c2})
for rec in c1:
    blob = " ".join(rec.values()).lower()
    if "fibroblast" in blob or rec["primary_site"].lower() == "skin":
        print("CELL", rec["cell_id"], rec["sample_type"], rec["primary_site"], rec["subtype"], rec["cell_type"], "age", rec["donor_age"])

oe_types = {"trt_oe", "trt_oe.mut"}
kd_types = {"trt_sh", "trt_sh.cgs", "trt_sh.css"}
oe_by_cell = defaultdict(set)
kd_by_cell = defaultdict(set)
oe_all = set()
kd_all = set()
ctrl_cells = set()
n_sig = 0
pert_types = Counter()
with gzip.open(root / "GSE92742_Broad_LINCS_sig_info.txt.gz", "rt", encoding="utf-8", errors="replace") as fh:
    for rec in csv.DictReader(fh, delimiter="\t"):
        n_sig += 1
        pt = rec["pert_type"]
        pert_types[pt] += 1
        name = rec["pert_iname"]
        cell = rec["cell_id"]
        if pt.startswith("ctl"):
            ctrl_cells.add(cell)
        if pt in oe_types and name in lambert:
            oe_by_cell[cell].add(name)
            oe_all.add(name)
        if pt in kd_types and name in lambert:
            kd_by_cell[cell].add(name)
            kd_all.add(name)
print("n_sig", n_sig)
print("pert_types", pert_types.most_common(12))
print("oe_lambert", len(oe_all), "kd_lambert", len(kd_all))
print("OE_BY_CELL")
for cell, genes in sorted(oe_by_cell.items(), key=lambda kv: -len(kv[1]))[:12]:
    print(cell, len(genes), "ctrl" if cell in ctrl_cells else "noctrl")
print("KD_TOP")
for cell, genes in sorted(kd_by_cell.items(), key=lambda kv: -len(kv[1]))[:8]:
    print(cell, len(genes))
for cid in ("MCH58", "FIBRNPC", "HA1E", "HUES3", "NPC"):
    print("FOCUS", cid, "oe", len(oe_by_cell.get(cid, ())), "kd", len(kd_by_cell.get(cid, ())), "ctrl", cid in ctrl_cells)
