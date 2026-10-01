"""CELLxGENE a19d1667 donor-identity audit and true-donor rerun of the ruler ρ and MD3 paired Δρ.

Reads frozen per-pseudobulk scores (results/md3/t2_cxg_a19d1667_S3_scores.csv) and the h5ad obs.
No model is refit. Stored statistics are reproduced first (tolerance 1e-6); any mismatch stops the run.
"""
from __future__ import annotations

import sys
from pathlib import Path

import h5py
import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parent))
from md3_common import MD3_DIR, MD3_SEED, MD3_BOOT, N_PERM, N_BOOT, ADULT_MIN_AGE  # noqa: E402
from md3_task2 import paired_delta_rho  # noqa: E402
from md2_common import CXG_H5AD  # noqa: E402
from fibro2_common import FIBRO2_DIR, rho_null_donor, strict_age_years  # noqa: E402
from fibro_common import bootstrap_rho_ci  # noqa: E402
from tissue_peek_obs import read_obs_column  # noqa: E402

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "results" / "paper_figs"
TOL = 1e-6
MIN_CELLS = 30

SCORES = MD3_DIR / "t2_cxg_a19d1667_S3_scores.csv"
RULER_RESULT = FIBRO2_DIR / "ta_cxg_a19d1667_result.json"
DELTA_CSV = MD3_DIR / "t2_delta_rho.csv"

OBS_FIELDS = {
    "Author": "source_study",
    "Donor identifier": "orig_donor_identifier",
    "Sample identifier": "orig_sample_identifier",
    "Internal sample identifier": "orig_internal_sample_identifier",
    "Accession (Sample)": "orig_accession_sample",
    "Age": "orig_author_age",
    "development_stage": "development_stage",
    "sex": "sex",
    "Sample location": "sample_location",
    "tissue": "tissue",
    "Condition": "condition",
    "disease": "disease",
    "assay": "assay",
}


def _log(fh, msg):
    print(msg, flush=True)
    fh.write(msg + "\n")
    fh.flush()


def _join(x):
    return " | ".join(sorted(set(map(str, x))))


def read_obs():
    with h5py.File(CXG_H5AD, "r") as h:
        cols = ["donor_id"] + list(OBS_FIELDS)
        obs = pd.DataFrame({c: np.asarray(read_obs_column(h, c)).astype(str) for c in cols})
    return obs


def build_audit(obs, scores):
    """One row per donor_id. n_cells_used applies the fibro2 filter: stated age, 10x assay, age ≥ 18."""
    age = obs["development_stage"].map(strict_age_years)
    age_num = pd.to_numeric(age, errors="coerce")
    is_10x = obs["assay"].str.lower().str.contains("10x")
    used = age_num.notna() & is_10x & (age_num >= ADULT_MIN_AGE)
    obs = obs.assign(_age=age_num, _used=used)
    agg = {c: (c, _join) for c in OBS_FIELDS}
    aud = obs.groupby("donor_id", sort=True).agg(
        **agg,
        age_years=("_age", lambda x: float(x.iloc[0]) if x.notna().all() and x.nunique() == 1 else np.nan),
        n_cells_total=("_age", "size"),
        n_cells_used=("_used", "sum"),
    ).rename(columns=OBS_FIELDS).reset_index()
    scored = set(scores["sample"].astype(str))
    aud["in_scored_93"] = aud.donor_id.isin(scored)
    aud["true_donor"] = aud.source_study + "::" + aud.orig_donor_identifier
    n_ids = aud[aud.in_scored_93].groupby("true_donor").donor_id.nunique()
    aud["n_scored_ids_same_true_donor"] = aud.true_donor.map(n_ids).fillna(0).astype(int)
    single_year = aud.orig_author_age.str.fullmatch(r"\d{1,3} y")
    aud["author_age_is_single_year"] = single_year
    aud["age_is_placeholder"] = (~single_year) & aud.age_years.notna()
    notes = []
    for _, r in aud.iterrows():
        n = []
        if r.n_scored_ids_same_true_donor > 1:
            n.append(f"{r.n_scored_ids_same_true_donor} scored donor_ids share Author+Donor identifier")
        if r.age_is_placeholder:
            n.append(f"author Age '{r.orig_author_age}' is not a single year; CxG stage '{r.development_stage}'")
        if " | " in r.orig_donor_identifier:
            n.append("multiple Donor identifier values within donor_id")
        notes.append("; ".join(n))
    aud["notes"] = notes
    return aud


def ruler_pack(y, pred, donor):
    sc, rho_null, p, _ = rho_null_donor(y, pred, donor, n_perm=N_PERM, seed=MD3_SEED)
    ci = bootstrap_rho_ci(y, pred, np.random.default_rng(MD3_BOOT), n_boot=N_BOOT)
    return dict(rho=sc["rho"], rho_null=rho_null, rho_p=p, rho_ci_lo=ci.get("p025"), rho_ci_hi=ci.get("p975"))


def _quiet(_msg):
    return None


def run_unit(df, unit_col, label):
    """Average pseudobulk scores within unit, then ruler ρ and both paired Δρ. Row order = first appearance."""
    g = df.groupby(unit_col, sort=False)
    ages = g.age_years.nunique()
    if (ages > 1).any():
        raise SystemExit(f"[{label}] age not constant within {unit_col}: {list(ages[ages > 1].index)}")
    d = g.agg(age_years=("age_years", "first"), age_score=("age_score", "mean"),
              md_score=("md_score", "mean"), s3=("age_up_minus_age_down", "mean"),
              n_pseudobulks=("sample", "size"), n_cells=("n_cells", "sum")).reset_index()
    y = d.age_years.to_numpy(float)
    unit = d[unit_col].astype(str).to_numpy()
    rows = []
    rp = ruler_pack(y, d.age_score.to_numpy(float), unit)
    rows.append(dict(analysis=label, contrast="ruler_spearman", n_units=len(d),
                     n_pseudobulks=int(d.n_pseudobulks.sum()), **rp))
    for other, col in (("MD_AddModuleScore", "md_score"), ("age_up_minus_age_down", "s3")):
        rec = paired_delta_rho(y, d.age_score.to_numpy(float), d[col].to_numpy(float), unit, _quiet, label)
        rows.append(dict(analysis=label, contrast=f"ruler_minus_{other}", n_units=rec["n_donors"],
                         n_pseudobulks=int(d.n_pseudobulks.sum()),
                         rho_ruler=rec["rho_ruler"], rho_other=rec["rho_other"], delta_rho=rec["delta_rho"],
                         delta_ci_lo=rec["delta_ci_lo"], delta_ci_hi=rec["delta_ci_hi"],
                         excludes_zero=rec["excludes_zero"], favours=rec["favours"]))
    return rows, d


def check_stored(rows, fh):
    """Pseudobulk-as-donor rows must equal fibro2 ruler JSON and md3 Δρ CSV to TOL."""
    import json
    stored_r = json.loads(RULER_RESULT.read_text(encoding="utf-8"))
    stored_d = pd.read_csv(DELTA_CSV)
    by = {r["contrast"]: r for r in rows}
    checks = []
    for k in ("rho", "rho_null", "rho_p", "rho_ci_lo", "rho_ci_hi"):
        checks.append((f"ruler.{k}", stored_r[k], by["ruler_spearman"][k], str(RULER_RESULT)))
    for tag, contrast in (("cxg_ruler_minus_MD", "ruler_minus_MD_AddModuleScore"),
                          ("cxg_ruler_minus_S3", "ruler_minus_age_up_minus_age_down")):
        s = stored_d[stored_d.tag == tag].iloc[0]
        for k in ("rho_ruler", "rho_other", "delta_rho", "delta_ci_lo", "delta_ci_hi"):
            checks.append((f"{tag}.{k}", float(s[k]), by[contrast][k], str(DELTA_CSV)))
        checks.append((f"{tag}.n_donors", float(s["n_donors"]), float(by[contrast]["n_units"]), str(DELTA_CSV)))
    worst = 0.0
    for name, a, b, src in checks:
        diff = abs(float(a) - float(b))
        worst = max(worst, diff)
        _log(fh, f"[check] {name}: stored={float(a):.12g} recomputed={float(b):.12g} |diff|={diff:.2e} ({Path(src).name})")
    if worst > TOL:
        raise SystemExit(f"STOP: stored CXG statistic not reproduced (max |diff|={worst:.3e} > {TOL}).")
    return worst


def main():
    OUT.mkdir(parents=True, exist_ok=True)
    fh = open(OUT / "cxg_true_donor_report.txt", "w", encoding="utf-8")
    try:
        scores = pd.read_csv(SCORES)
        if len(scores) != 93:
            raise SystemExit(f"STOP: {SCORES.name} has {len(scores)} rows, expected 93")
        obs = read_obs()
        _log(fh, f"[obs] {CXG_H5AD.name} n_cells={len(obs):,} donor_ids={obs.donor_id.nunique()}")
        aud = build_audit(obs, scores)

        m = scores.merge(aud[["donor_id", "n_cells_used", "true_donor", "age_years", "age_is_placeholder"]]
                         .rename(columns={"age_years": "age_audit"}),
                         left_on="sample", right_on="donor_id", how="left")
        bad_n = m[m.n_cells != m.n_cells_used]
        bad_a = m[m.age_years != m.age_audit]
        if len(bad_n) or len(bad_a):
            raise SystemExit(f"STOP: audit disagrees with scored pseudobulks n_cells={len(bad_n)} age={len(bad_a)}")
        _log(fh, "[audit] n_cells_used and age match all 93 scored pseudobulks")
        aud.to_csv(OUT / "cxg_donor_audit.csv", index=False)

        rey = aud[aud.in_scored_93 & (aud.source_study == "Reynolds_2021")]
        _log(fh, f"[Reynolds_2021] scored donor_ids={len(rey)} ages={sorted(rey.age_years.unique())} "
                 f"orig Donor identifier={sorted(rey.orig_donor_identifier.unique())} "
                 f"author Age={sorted(rey.orig_author_age.unique())} "
                 f"accessions unique={rey.orig_accession_sample.nunique()}")
        for td, sub in rey.groupby("true_donor"):
            _log(fh, f"   {td}: {len(sub)} donor_ids, n_cells_used={int(sub.n_cells_used.sum())}")

        all_rows = []
        rows0, _ = run_unit(m, "sample", "stored_reproduction_pseudobulk_as_donor")
        worst = check_stored(rows0, fh)
        _log(fh, f"[check] PASS max |diff|={worst:.2e}")
        all_rows += rows0

        rows1, d1 = run_unit(m, "true_donor", "true_donor")
        all_rows += rows1
        m30 = m[m.n_cells >= MIN_CELLS]
        rows2, d2 = run_unit(m30, "true_donor", f"true_donor_pseudobulk_ge{MIN_CELLS}cells")
        all_rows += rows2
        mnp = m[~m.age_is_placeholder.astype(bool)]
        rows3, d3 = run_unit(mnp, "true_donor", "supplement_true_donor_excl_placeholder_age")
        all_rows += rows3

        res = pd.DataFrame(all_rows)
        res.insert(1, "unit", res.analysis.map(lambda a: "donor_id (pseudobulk)" if a.startswith("stored") else "Author::Donor identifier"))
        res["seed_perm"] = MD3_SEED
        res["seed_boot"] = MD3_BOOT
        res["n_perm"] = N_PERM
        res["n_boot"] = N_BOOT
        res.to_csv(OUT / "cxg_true_donor_results.csv", index=False)
        d1.assign(analysis="true_donor").to_csv(OUT / "cxg_true_donor_scores.csv", index=False)
        for r in all_rows:
            _log(fh, "[result] " + ", ".join(f"{k}={v:.4f}" if isinstance(v, float) else f"{k}={v}" for k, v in r.items()))
        _log(fh, f"[n] pseudobulks=93 true_donors={len(d1)} true_donors_ge{MIN_CELLS}={len(d2)} "
                 f"pseudobulks_ge{MIN_CELLS}={len(m30)} excl_placeholder_true_donors={len(d3)}")
    finally:
        fh.close()


if __name__ == "__main__":
    main()
