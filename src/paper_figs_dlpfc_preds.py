"""Donor-level (age, predicted age) for POSCTRL S3 gene_ridge HBCC?묺SSM and MSSM?묱BCC.

Reads frozen predictions (results/posctrl/s3_gene_ridge_preds.npz); no refit. The stored transfer ?
(s3_gene_ridge.json rho_H_to_M / rho_M_to_H) is the median over cell types of per-type Spearman ?,
where each type's rows are test-bank donors. That is reproduced to 1e-6 before anything is written.
"""
from __future__ import annotations

import json
import sys
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import spearmanr

sys.path.insert(0, str(Path(__file__).resolve().parent))
from external_common import median_over_types  # noqa: E402

ROOT = Path(__file__).resolve().parents[1]
PC = ROOT / "results" / "posctrl"
OUT = ROOT / "results" / "paper_figs"
TOL = 1e-6
DIRS = {"H_to_M": "HBCC?묺SSM", "M_to_H": "MSSM?묱BCC"}


def spearman(a, b):
    """No minimum n: the stored per-type table has ? for types with n_te < 8 (e.g. smooth muscle, n=4)."""
    a, b = np.asarray(a, float), np.asarray(b, float)
    m = np.isfinite(a) & np.isfinite(b)
    if int(m.sum()) < 3 or np.std(a[m]) < 1e-12 or np.std(b[m]) < 1e-12:
        return np.nan
    return float(spearmanr(a[m], b[m])[0])


def main():
    OUT.mkdir(parents=True, exist_ok=True)
    z = np.load(PC / "s3_gene_ridge_preds.npz", allow_pickle=True)
    stored = json.loads((PC / "s3_gene_ridge.json").read_text(encoding="utf-8"))
    per_stored = pd.read_csv(PC / "s3_gene_ridge_transfer_per_type.csv")
    lines, longs, dons = [], [], []
    for key, label in DIRS.items():
        df = pd.DataFrame(dict(direction=label, donor=z[f"{key}_donor"].astype(str),
                               celltype=z[f"{key}_celltype"].astype(str),
                               age=z[f"{key}_y"].astype(float), pred=z[f"{key}_pred"].astype(float)))
        if df.duplicated(["donor", "celltype"]).any():
            raise SystemExit(f"STOP: {key} has duplicate donor횞type rows")
        rho_t = df.groupby("celltype").apply(lambda g: spearman(g.age, g.pred), include_groups=False)
        rho = median_over_types(rho_t.to_numpy())
        tr, te = key.split("_to_")
        ps = per_stored[(per_stored.train_site == tr) & (per_stored.test_site == te)].set_index("celltype")
        d_type = max(abs(float(ps.loc[t, "rho"]) - float(v)) for t, v in rho_t.items() if np.isfinite(v))
        d = abs(rho - float(stored[f"rho_{key}"]))
        lines.append(f"[{label}] median-over-types rho stored={stored[f'rho_{key}']:.12g} recomputed={rho:.12g} "
                     f"|diff|={d:.2e}; per-type max |diff| vs s3_gene_ridge_transfer_per_type.csv={d_type:.2e}")
        if d > TOL or d_type > TOL:
            print("\n".join(lines))
            raise SystemExit(f"STOP: {label} stored rho not reproduced")
        g = df.groupby("donor").agg(age=("age", "first"), n_types=("celltype", "nunique"),
                                    pred_median_over_types=("pred", "median")).reset_index()
        g.insert(0, "direction", label)
        rho_don = spearman(g.age, g.pred_median_over_types)
        lines.append(f"[{label}] n_donors={len(g)} n_rows={len(df)}; rho(age, donor median pred)={rho_don:.6f} "
                     "(descriptive only; not the stored statistic)")
        longs.append(df)
        dons.append(g)
    pd.concat(longs, ignore_index=True).to_csv(OUT / "dlpfc_gene_ridge_donor_type_preds.csv", index=False)
    pd.concat(dons, ignore_index=True).to_csv(OUT / "dlpfc_gene_ridge_donor_preds.csv", index=False)
    (OUT / "dlpfc_gene_ridge_report.txt").write_text("\n".join(lines) + "\n", encoding="utf-8")
    print("\n".join(lines))


if __name__ == "__main__":
    main()
