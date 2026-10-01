"""SEA-AD Supersession 1 (FINDINGS_EXTERNAL.md, 2026-09-17): ρ-gated primary and ≥65-only primary.

Spearman ρ = median over mapped cell types (n ≥ MIN_N_R donors) of per-type ρ on frozen DLPFC
predictions (results/external/project_pseudobulk_obs.csv). Donor-level age permutation null
(n_perm=200, fresh Generator on 20260914 per cell) and donor bootstrap percentile CI (B=200, 20260918).
No DLPFC or SEA-AD refit. Reads only; writes to results/paper_figs/.
"""
from __future__ import annotations

import json
import sys
from pathlib import Path

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parent))
from external_common import EXT_DIR, EXT_SEED, BOOT_SEED, N_PERM, N_BOOT, MIN_N_R  # noqa: E402
from external_project import (  # noqa: E402
    _type_vectors, per_type_metrics, permute_median_r, bootstrap_median_r,
)
from trajectory_common import permutation_p, summarize_null_col  # noqa: E402

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "results" / "paper_figs"
TOL = 1e-6
AGE_MIN_S1 = 65
PB_OBS = EXT_DIR / "project_pseudobulk_obs.csv"
DONOR_META = EXT_DIR / "stage0_donor_meta.csv"

# FINDINGS_EXTERNAL.md "S1 ≥65-only primary" and "S1 ρ — original primary" (3-decimal prose only).
S1_PROSE = {
    ("primary", "raw"): dict(rho=0.478, rho_null=0.002, rho_p=0.005, rho_ci_lo=0.026, rho_ci_hi=0.742, n_donors=26),
    ("primary", "idresid"): dict(rho=0.491, rho_null=0.002, rho_p=0.005, rho_ci_lo=0.074, rho_ci_hi=0.741, n_donors=26),
    ("ge65", "raw"): dict(rho=-0.009, rho_null=0.017, rho_p=0.557, rho_ci_lo=-0.348, rho_ci_hi=0.337, n_donors=21),
    ("ge65", "idresid"): dict(rho=-0.014, rho_null=0.015, rho_p=0.567, rho_ci_lo=-0.301, rho_ci_hi=0.333, n_donors=21),
}
PROSE_TOL = 5e-4 + 1e-12


def _log(fh, msg):
    print(msg, flush=True)
    fh.write(msg + "\n")
    fh.flush()


def permute_median_rho(type_dfs, rng, n_perm):
    """external_project.permute_median_r with med['rho'] in place of med['r']."""
    donors = sorted(set().union(*[set(d.donor.astype(str)) for d in type_dfs.values()]))
    age_map0 = {}
    for d in type_dfs.values():
        for _, row in d.iterrows():
            age_map0[str(row.donor)] = float(row.age)
    ages0 = np.array([age_map0[d] for d in donors], float)
    null = []
    for _ in range(n_perm):
        amap = dict(zip(donors, rng.permutation(ages0)))
        tdfs = {}
        for t, d in type_dfs.items():
            dd = d.copy()
            dd["age"] = dd.donor.astype(str).map(amap).astype(float)
            tdfs[t] = dd
        _, med = per_type_metrics(tdfs)
        null.append(med["rho"])
    return np.asarray(null, float)


def bootstrap_median_rho(type_dfs, rng, n_boot):
    """external_project.bootstrap_median_r with med['rho'] in place of med['r']."""
    donors = np.array(sorted(set().union(*[set(d.donor.astype(str)) for d in type_dfs.values()])), dtype=object)
    by = {t: d.set_index(d.donor.astype(str)) for t, d in type_dfs.items()}
    vals = []
    for _ in range(n_boot):
        draw = rng.choice(donors, size=len(donors), replace=True)
        tdfs = {}
        for t, d in by.items():
            rows = []
            for k, don in enumerate(draw):
                if don in d.index:
                    rec = d.loc[don]
                    if isinstance(rec, pd.DataFrame):
                        rec = rec.iloc[0]
                    rec = rec.copy()
                    rec["donor"] = f"{don}#{k}"
                    rows.append(rec)
            tdfs[t] = pd.DataFrame(rows).reset_index(drop=True) if rows else pd.DataFrame(columns=d.reset_index().columns)
        _, med = per_type_metrics(tdfs)
        vals.append(med["rho"])
    v = np.asarray(vals, float)
    v = v[np.isfinite(v)]
    if v.size < 2:
        return dict(p025=np.nan, p975=np.nan, n=int(v.size))
    return dict(p025=float(np.percentile(v, 2.5)), p975=float(np.percentile(v, 97.5)), n=int(v.size))


def load_inputs():
    pb = pd.read_csv(PB_OBS)
    pb["donor"] = pb.donor.astype(str)
    pb["celltype"] = pb.celltype.astype(str)
    donors = pd.read_csv(DONOR_META)
    donors["donor"] = donors.donor.astype(str)
    donors["age"] = pd.to_numeric(donors["age"], errors="coerce")
    pl = donors["pathology_low"]
    if pl.dtype != bool:
        donors["pathology_low"] = pl.map(lambda x: str(x).strip().lower() in ("true", "1", "yes"))
    return pb, donors


def check_stored_r(pb, prim_don, fh):
    """Frozen inputs must reproduce primary_{kind}.json (full precision) before any ρ cell."""
    worst = 0.0
    for kind in ("raw", "idresid"):
        stored = json.loads((EXT_DIR / f"primary_{kind}.json").read_text(encoding="utf-8"))
        pred = pb[f"pred_{kind}"].to_numpy(float)
        tdf = _type_vectors(pb.assign(pred=pred), pred, donors_keep=set(prim_don))
        _, med = per_type_metrics(tdf)
        null = permute_median_r(tdf, np.random.default_rng(EXT_SEED), N_PERM)
        ci = bootstrap_median_r(tdf, np.random.default_rng(BOOT_SEED), N_BOOT)
        got = dict(r=med["r"], rho=med["rho"], cal_r2=med["cal_r2"], r2=med["r2"],
                   r_null=summarize_null_col(null)["mean"], r_p=permutation_p(med["r"], null, greater=True),
                   r_ci_lo=ci["p025"], r_ci_hi=ci["p975"])
        for k, v in got.items():
            d = abs(float(stored[k]) - float(v))
            worst = max(worst, d)
            _log(fh, f"[check primary_{kind}.json] {k}: stored={float(stored[k]):.12g} recomputed={float(v):.12g} |diff|={d:.2e}")
    if worst > TOL:
        raise SystemExit(f"STOP: frozen SEA-AD primary cell not reproduced (max |diff|={worst:.3e} > {TOL}).")
    return worst


def rho_cell(pb, kind, donors_keep):
    pred = pb[f"pred_{kind}"].to_numpy(float)
    tdf = _type_vectors(pb.assign(pred=pred), pred, donors_keep=set(donors_keep))
    tab, med = per_type_metrics(tdf)
    null = permute_median_rho(tdf, np.random.default_rng(EXT_SEED), N_PERM)
    ns = summarize_null_col(null)
    p = permutation_p(med["rho"], null, greater=True)
    ci = bootstrap_median_rho(tdf, np.random.default_rng(BOOT_SEED), N_BOOT)
    n_don = int(len(set().union(*[set(d.donor.astype(str)) for d in tdf.values()])))
    rec = dict(rho=med["rho"], rho_null=ns["mean"], rho_null_median=ns["median"], rho_p=p,
               rho_ci_lo=ci["p025"], rho_ci_hi=ci["p975"], r=med["r"], cal_r2=med["cal_r2"],
               n_donors=n_don, n_types=med["n_types_ok"], n_perm=N_PERM, n_boot=N_BOOT,
               pass_rho=bool(np.isfinite(med["rho"]) and med["rho"] > 0 and p <= 0.05))
    return rec, tab


def main():
    OUT.mkdir(parents=True, exist_ok=True)
    fh = open(OUT / "seaad_s1_report.txt", "w", encoding="utf-8")
    try:
        pb, donors = load_inputs()
        prim = donors.loc[donors.pathology_low.astype(bool) & donors.age.notna(), ["donor", "age"]]
        ge65 = prim[prim.age >= AGE_MIN_S1]
        _log(fh, f"[donors] primary no/low ADNC={prim.donor.nunique()}  >= {AGE_MIN_S1}={ge65.donor.nunique()} "
                 f"unique ages >= {AGE_MIN_S1}={ge65.age.nunique()}")
        if ge65.donor.nunique() < MIN_N_R:
            raise SystemExit(f"STOP: >= {AGE_MIN_S1} cell has {ge65.donor.nunique()} < {MIN_N_R} donors")
        worst = check_stored_r(pb, prim.donor, fh)
        _log(fh, f"[check] PASS primary r/rho/null/CI max |diff|={worst:.2e}")

        rows, tabs = [], []
        for cell, keep in (("primary", prim.donor), ("ge65", ge65.donor)):
            for kind in ("raw", "idresid"):
                rec, tab = rho_cell(pb, kind, keep)
                rows.append(dict(cell=cell, direction=kind, **rec))
                tabs.append(tab.assign(cell=cell, direction=kind))
                _log(fh, f"[S1 {cell} {kind}] rho={rec['rho']:+.6f} null={rec['rho_null']:+.4f} p={rec['rho_p']:.4f} "
                         f"CI=[{rec['rho_ci_lo']:+.4f}, {rec['rho_ci_hi']:+.4f}] r={rec['r']:+.4f} "
                         f"n={rec['n_donors']} types={rec['n_types']} pass={rec['pass_rho']}")
        res = pd.DataFrame(rows)

        # ρ primary raw is stored at full precision in primary_raw.json; the rest exists only as 3-dp prose.
        stored_rho = json.loads((EXT_DIR / "primary_raw.json").read_text(encoding="utf-8"))["rho"]
        got = float(res.loc[(res.cell == "primary") & (res.direction == "raw"), "rho"].iloc[0])
        if abs(got - stored_rho) > TOL:
            raise SystemExit(f"STOP: primary raw rho {got} != stored {stored_rho}")
        _log(fh, f"[check] primary raw rho stored={stored_rho:.12g} recomputed={got:.12g} |diff|={abs(got - stored_rho):.2e}")
        prose_fail = []
        for (cell, kind), exp in S1_PROSE.items():
            r = res[(res.cell == cell) & (res.direction == kind)].iloc[0]
            for k, v in exp.items():
                d = abs(float(r[k]) - float(v))
                ok = d <= (0 if k == "n_donors" else PROSE_TOL)
                _log(fh, f"[check prose {cell} {kind}] {k}: FINDINGS={v} recomputed={float(r[k]):.6f} ok={ok}")
                if not ok:
                    prose_fail.append((cell, kind, k, v, float(r[k])))
        if prose_fail:
            res.to_csv(OUT / "seaad_s1_rho_UNMATCHED.csv", index=False)
            raise SystemExit(f"STOP: S1 prose values not reproduced: {prose_fail}")
        _log(fh, "[check] PASS all S1 values match FINDINGS_EXTERNAL.md at the printed 3 decimals")
        res.to_csv(OUT / "seaad_s1_rho.csv", index=False)
        pd.concat(tabs, ignore_index=True).to_csv(OUT / "seaad_s1_rho_per_type.csv", index=False)

        long = pb.merge(donors[["donor", "pathology_low"]], on="donor", how="left", suffixes=("", "_meta"))
        long["in_primary"] = long.donor.isin(set(prim.donor))
        long["in_ge65"] = long.donor.isin(set(ge65.donor))
        long[["donor", "celltype", "n_cells", "age", "sex", "pathology_low", "in_primary", "in_ge65",
              "pred_raw", "score_raw", "pred_idresid", "score_idresid"]].to_csv(
            OUT / "seaad_donor_type_scores.csv", index=False)
        don = long.groupby("donor").agg(
            age=("age", "first"), sex=("sex", "first"), pathology_low=("pathology_low", "first"),
            in_primary=("in_primary", "first"), in_ge65=("in_ge65", "first"), n_types=("celltype", "nunique"),
            raw_pred_median_over_types=("pred_raw", "median"), raw_score_median_over_types=("score_raw", "median"),
            idresid_pred_median_over_types=("pred_idresid", "median"),
            idresid_score_median_over_types=("score_idresid", "median"),
        ).reset_index()
        don.to_csv(OUT / "seaad_donor_scores.csv", index=False)
        _log(fh, f"[out] seaad_donor_scores.csv donors={len(don)}  seaad_donor_type_scores.csv rows={len(long)}")
    finally:
        fh.close()


if __name__ == "__main__":
    main()
