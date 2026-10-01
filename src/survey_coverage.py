"""Ruler-gene coverage for the TF-perturbation survey.

Reads the frozen fibroblast ruler at results/fibro/frozen_ruler_ridge_raw.npz
(gene symbols and weights only). Does not refit, rescale, or rescore.

Measured genes are one-symbol-per-line text files listed in
results/survey/gene_lists/manifest.csv. Those lists must come from a
feature or metadata file that was actually opened. Inferred genes
(for example L1000 genes that are not landmarks) must not be written
into those lists.

Writes results/survey/coverage.csv.
"""
from __future__ import annotations

import csv
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
RULER = ROOT / "results" / "fibro" / "frozen_ruler_ridge_raw.npz"
MANIFEST = ROOT / "results" / "survey" / "gene_lists" / "manifest.csv"
OUT = ROOT / "results" / "survey" / "coverage.csv"


def load_ruler():
    with np.load(RULER, allow_pickle=False) as z:
        symbols = [str(s) for s in z["symbol"]]
        weights = np.abs(np.asarray(z["w"], dtype=np.float64))
    if len(symbols) != weights.shape[0]:
        raise SystemExit("ruler symbol and weight lengths differ")
    # Collapse duplicate symbols by summing |weight|. The first occurrence
    # keeps the symbol; ENSG-placeholder symbols stay as stored.
    order = {}
    for sym, w in zip(symbols, weights):
        key = sym.strip()
        order[key] = order.get(key, 0.0) + float(w)
    keys = list(order)
    w = np.array([order[k] for k in keys], dtype=np.float64)
    return keys, w


def load_measured(path: Path) -> set[str]:
    genes = set()
    for line in path.read_text(encoding="utf-8", errors="replace").splitlines():
        g = line.strip()
        if not g or g.startswith("#"):
            continue
        genes.add(g)
    return genes


def main() -> None:
    symbols, weights = load_ruler()
    total_w = float(weights.sum())
    n_ruler = len(symbols)
    if total_w <= 0 or n_ruler == 0:
        raise SystemExit("ruler weights are empty")
    rows = []
    with MANIFEST.open(encoding="utf-8", newline="") as fh:
        for rec in csv.DictReader(fh):
            path = ROOT / rec["gene_list"]
            measured = load_measured(path)
            # Match on the symbol string as stored. Case-fold only the
            # comparison key so ENSG placeholders are unchanged in the report.
            measured_key = {g.upper() for g in measured}
            hit = np.array([s.upper() in measured_key for s in symbols], dtype=bool)
            n_hit = int(hit.sum())
            w_hit = float(weights[hit].sum())
            rows.append(
                {
                    "dataset": rec["dataset"],
                    "tier": rec["tier"],
                    "gene_list": rec["gene_list"],
                    "n_measured_symbols": len(measured),
                    "n_ruler_genes": n_ruler,
                    "n_ruler_genes_measured": n_hit,
                    "fraction_ruler_genes": n_hit / n_ruler,
                    "fraction_abs_weight": w_hit / total_w,
                    "note": rec.get("note", ""),
                }
            )
    OUT.parent.mkdir(parents=True, exist_ok=True)
    fields = [
        "dataset",
        "tier",
        "gene_list",
        "n_measured_symbols",
        "n_ruler_genes",
        "n_ruler_genes_measured",
        "fraction_ruler_genes",
        "fraction_abs_weight",
        "note",
    ]
    with OUT.open("w", encoding="utf-8", newline="") as fh:
        writer = csv.DictWriter(fh, fieldnames=fields)
        writer.writeheader()
        writer.writerows(rows)
    print(f"wrote coverage.csv n={len(rows)} ruler_genes={n_ruler}")


if __name__ == "__main__":
    main()
