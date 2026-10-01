"""Compile the T4 comparison table and the cross-tissue figure from the per-atlas T2/T3 CSVs."""
import sys
from pathlib import Path

import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

sys.path.insert(0, str(Path(__file__).resolve().parent))
from tissue_common import TISSUE_DIR, TISSUE_FIG, BLOOD_COMP_R2, BLOOD_EXPR_R2  # noqa: E402

REPROG = {
    "retina": "human+mouse RPE (GSE304044 SuperSeries); mouse retina aging multiome (GSE307031). No human retina OSK scRNA found.",
    "brain": "mouse mPFC OSK multiome (GSE276656); mouse SVZ partial reprogramming scRNA (GSE224438). No human brain OSK scRNA found.",
    "skeletal muscle": "no tissue OSK scRNA. Closest: C2C12 / multi-type murine cells (GSE176206). Phase 2 incomplete.",
    "heart": "no OSK scRNA. Closest: PHF7 cardiac reprogramming (GSE270268), not OSK.",
    "kidney": "none found in GEO reprogramming queries.",
    "liver": "mouse in vivo (GSE144600, GSE274988, GSE201710). High-turnover tissue.",
    "blood": "closed — OneK1K Stage B. Not a candidate.",
}


def main():
    t2s, t3s = [], []
    for p in sorted(TISSUE_DIR.glob("t2_*_batch_audit.csv")):
        if p.name == "t2_batch_audit.csv":
            continue
        t2s.append(pd.read_csv(p))
    for p in sorted(TISSUE_DIR.glob("t3_*_summary.csv")):
        if "_cv_" in p.name:
            continue
        t3s.append(pd.read_csv(p))
    t2 = pd.concat(t2s, ignore_index=True) if t2s else pd.DataFrame()
    t3 = pd.concat(t3s, ignore_index=True) if t3s else pd.DataFrame()
    t2.to_csv(TISSUE_DIR / "t2_batch_audit.csv", index=False)
    t3.to_csv(TISSUE_DIR / "t3_turnover.csv", index=False)

    blood = dict(name="blood_OneK1K", tissue="blood (closed)", n_donors=981, age_span="19–97",
                 batch_auditability="yes — 75 10x pools", batch_column="pool",
                 r2_age_batch=0.23, r2_depth_batch=0.86,
                 composition_r2=BLOOD_COMP_R2, expression_r2=BLOOD_EXPR_R2,
                 delta=BLOOD_EXPR_R2 - BLOOD_COMP_R2,
                 expression_beats_composition=False,
                 reprogramming=REPROG["blood"])
    rows = [blood]
    for _, r in t3.iterrows():
        t2r = t2[t2.name == r["name"]]
        t2r = t2r.iloc[0] if len(t2r) else None
        span = f"{t2r.age_min:.0f}–{t2r.age_max:.0f}" if t2r is not None else ""
        aud = ("yes — " + str(r.batch_column)) if r.batch_usable else "NO shared-donor batch"
        rows.append(dict(
            name=r["name"], tissue=r.tissue, n_donors=int(r.n_donors), age_span=span,
            batch_auditability=aud, batch_column=r.batch_column,
            r2_age_batch=None if t2r is None else t2r.r2_age_batch,
            r2_depth_batch=None if t2r is None else t2r.r2_depth_batch,
            composition_r2=float(r.composition_r2), expression_r2=float(r.expression_r2),
            delta=float(r.expression_r2) - float(r.composition_r2),
            expression_beats_composition=bool(r.expression_beats_composition),
            reprogramming=REPROG.get(r.tissue, ""),
        ))
    tab = pd.DataFrame(rows).sort_values("delta", ascending=False)
    tab.to_csv(TISSUE_DIR / "t4_comparison.csv", index=False)

    fig, ax = plt.subplots(figsize=(8.2, 4.4))
    names = [f"{r.tissue}\n({r['name']})" for _, r in tab.iterrows()]
    x = np.arange(len(tab))
    w = 0.38
    ax.bar(x - w/2, tab.composition_r2, w, color="tab:orange", label="composition (CLR)")
    ax.bar(x + w/2, tab.expression_r2, w, color="tab:blue", label="within-type expression")
    ax.axhline(0, c="k", lw=0.6)
    ax.set_xticks(x)
    ax.set_xticklabels(names, fontsize=8)
    ax.set_ylabel("held-out within-batch R²")
    ax.set_title("T3 turnover test vs blood Stage B (composition 0.493, expression 0.482)")
    ax.legend(fontsize=8)
    fig.tight_layout()
    fig.savefig(TISSUE_FIG / "t4_comp_vs_expr.png", dpi=140)
    plt.close(fig)
    print(tab.to_string(index=False))
    print("wrote", TISSUE_DIR / "t4_comparison.csv")


if __name__ == "__main__":
    main()
