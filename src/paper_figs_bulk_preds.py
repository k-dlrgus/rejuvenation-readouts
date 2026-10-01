"""Donor-level (age, predicted age) for GTEx T3 BA9↔Cortex and GTEx fibroblast B1↔C1, raw ridge.

The original runs stored only summary statistics. This re-executes the original observed-label code
path unchanged (gtex_stage1.run_transfer_pair lines before the permutation loop; fibro_stage1._transfer_one
likewise): same packs, same z-scoring, same closed-form ridge and alpha grid, no seeds. Each direction must
reproduce stored ρ, r, R², calibrated R² and n to 1e-6 or the run stops. Permutation nulls are not re-run.
"""
from __future__ import annotations

import sys
from pathlib import Path

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parent))
from gtex_common import PRIMARY_TISSUE, PAIRED_TISSUE, TISSUE_SLUG, GTEX_DIR, pred_scores  # noqa: E402
from gtex_stage1 import load_tissue_pack, _transfer_align, transfer_prepare, transfer_score  # noqa: E402
import fibro_stage1 as fs1  # noqa: E402
from fibro_common import FIBRO_DIR, BOOT_SEED, N_BOOT, bootstrap_rho_ci  # noqa: E402

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "results" / "paper_figs"
TOL = 1e-6
CHECK_KEYS = ("rho", "r", "r2", "cal_r2")


def _compare(label, stored_row, got, extra, lines):
    worst = 0.0
    for k in CHECK_KEYS:
        d = abs(float(stored_row[k]) - float(got[k]))
        worst = max(worst, d)
        lines.append(f"[{label}] {k}: stored={float(stored_row[k]):.12g} recomputed={float(got[k]):.12g} |diff|={d:.2e}")
    for k, v in extra.items():
        d = abs(float(stored_row[k]) - float(v))
        worst = max(worst, d)
        lines.append(f"[{label}] {k}: stored={float(stored_row[k]):.12g} recomputed={float(v):.12g} |diff|={d:.2e}")
    if worst > TOL:
        print("\n".join(lines))
        raise SystemExit(f"STOP: {label} stored statistic not reproduced (max |diff|={worst:.3e})")


def gtex_t3(lines):
    stored = pd.read_csv(GTEX_DIR / "t3.csv")
    ba9 = load_tissue_pack(TISSUE_SLUG[PRIMARY_TISSUE])
    ctx = load_tissue_pack(TISSUE_SLUG[PAIRED_TISSUE])
    out = []
    for src, dst in ((ba9, ctx), (ctx, ba9)):
        label = f"{src['tissue']} -> {dst['tissue']}"
        align = _transfer_align(src, dst)
        prep = transfer_prepare(src, dst, "raw", "ridge", align)
        sc = transfer_score(prep, "ridge")
        row = stored[(stored.label == label) & (stored.regime == "raw") & (stored.method == "ridge")].iloc[0]
        _compare(f"GTEx {label}", row, sc,
                 dict(n_train=prep["n_train"], n_test=prep["n_test"], n_genes=align["n_genes"]), lines)
        te_obs = dst["obs"].iloc[align["te"]].reset_index(drop=True)
        df = pd.DataFrame(dict(
            direction=f"{TISSUE_SLUG[src['tissue']]}_to_{TISSUE_SLUG[dst['tissue']]}",
            train_tissue=src["tissue"], test_tissue=dst["tissue"],
            donor=te_obs.donor.astype(str), sampid=te_obs.SAMPID.astype(str),
            age_bracket=te_obs.AGE.astype(str) if "AGE" in te_obs else "",
            age_mid=np.asarray(prep["yte"], float), pred=np.asarray(sc["pred"], float),
        ))
        if not np.allclose(df.age_mid.to_numpy(), te_obs.age_mid.to_numpy(float)):
            raise SystemExit(f"STOP: {label} test age order mismatch")
        out.append(df)
    return pd.concat(out, ignore_index=True)


def fibro_b1c1(lines):
    stored = pd.read_csv(FIBRO_DIR / "stage1_transfer.csv")
    pack = fs1.load_pack()
    obs, X = pack["obs"], pack["X"]
    cvec = obs.SMCENTER.astype(str).to_numpy()
    out = []
    for tr_s, te_s in (("B1", "C1"), ("C1", "B1")):
        tr = np.flatnonzero(cvec == tr_s)
        te = np.flatnonzero(cvec == te_s)
        Xtr, Xte, ytr, yte, _ = fs1.prepare_fold_X(X, obs, tr, te, "raw", "ridge")
        svd = fs1.ridge_prestd(Xtr, ytr)[3]
        pred, _, _, _ = fs1.predict_fold(Xtr, ytr, Xte, "ridge", svd=svd)
        sc = pred_scores(yte, pred)
        ci = bootstrap_rho_ci(yte, pred, np.random.default_rng(BOOT_SEED), n_boot=N_BOOT)
        row = stored[(stored.direction == f"{tr_s}→{te_s}") & (stored.regime == "raw") & (stored.method == "ridge")].iloc[0]
        _compare(f"fibro {tr_s}→{te_s}", row, sc,
                 dict(n_train=len(tr), n_test=len(te), rho_ci_lo=ci["p025"], rho_ci_hi=ci["p975"]), lines)
        te_obs = obs.iloc[te].reset_index(drop=True)
        out.append(pd.DataFrame(dict(
            direction=f"{tr_s}_to_{te_s}", train_center=tr_s, test_center=te_s,
            donor=te_obs.donor.astype(str), sampid=te_obs.SAMPID.astype(str),
            age_bracket=te_obs.AGE.astype(str) if "AGE" in te_obs else "",
            age_mid=np.asarray(yte, float), pred=np.asarray(pred, float),
        )))
    return pd.concat(out, ignore_index=True)


def main():
    OUT.mkdir(parents=True, exist_ok=True)
    lines = []
    g = gtex_t3(lines)
    f = fibro_b1c1(lines)
    for df, name in ((g, "gtex_cortex_ridge_donor_preds.csv"), (f, "gtex_fibroblast_ridge_donor_preds.csv")):
        if df.groupby("direction").donor.apply(lambda s: s.duplicated().any()).any():
            raise SystemExit(f"STOP: {name} has more than one sample per donor within a direction")
        df.to_csv(OUT / name, index=False)
    lines.append(f"[out] GTEx T3 rows={len(g)} fibro rows={len(f)} (one sample per donor per direction)")
    (OUT / "bulk_ridge_preds_report.txt").write_text("\n".join(lines) + "\n", encoding="utf-8")
    print("\n".join(lines))


if __name__ == "__main__":
    main()
