"""Table 1 of paper2_passage_draft_v0.6, rebuilt from stored angles. Recomputes nothing.

usage (from the repository root): python results/paper2_figs_v06/table1.py
writes results/paper2_figs_v06/table1.csv

Angle from parallel = min(angle_deg, 180 - angle_deg). Lu day-0 values come from the
corrected file (results/axang/CORRECTION.md). Each cell is "20-component (full gene)".
"""
import csv
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
AX = ROOT / "results" / "axang"
ROWS = [("C1_south", "Southard control cells"),
        ("C1_seng_CRA", "Sengstack activation-screen control cells (half)"),
        ("C1_seng_CRI", "Sengstack interference-screen control cells (half)"),
        ("C1_lu", "Lu day-0 fibroblasts"),
        ("C2", "Cell-cycle gene list, equal weights"),
        ("C1_gtex", "GTEx bulk tissue samples")]
COLS = [("A3", "Passage (Sengstack)"), ("A1", "Donor age (GTEx)"),
        ("A4", "Donor pair (Lu, 96 vs 22 y)"), ("A2", "Age ruler (trained on GTEx)")]


def parallel(deg):
    return min(deg, 180.0 - deg)


ang = {}
with open(AX / "d1_angles.csv", newline="", encoding="utf-8") as fh:
    for r in csv.DictReader(fh):
        if r["kind"] == "aging-prolif":
            ang[(r["space"], r["axis_a"], r["axis_b"])] = parallel(float(r["angle_deg"]))
with open(AX / "d1_angles_corrected_C1lu.csv", newline="", encoding="utf-8") as fh:
    for r in csv.DictReader(fh):
        ang[(r["space"], r["axis_a"], r["axis_b"])] = float(r["angle_from_parallel_corrected"])

out = Path(__file__).with_name("table1.csv")
with open(out, "w", newline="", encoding="utf-8") as fh:
    w = csv.writer(fh)
    w.writerow(["cell_division_direction"] + [lab for _, lab in COLS])
    for key, lab in ROWS:
        w.writerow([lab] + [f"{ang[('S2A_k20', a, key)]:.1f} ({ang[('full_gene', a, key)]:.1f})" for a, _ in COLS])
print(out.read_text(encoding="utf-8"))
