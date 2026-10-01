"""Donor-level profile of a CELLxGENE h5ad: donors, ages, sex, cell types, counts check.
Evaluates STOP CONDITION 0 (>=10 donors, continuous age, >=40y age range, >=3 cell types).
"""
import sys
import numpy as np
import pandas as pd
import h5py
import anndata as ad

sys.path.insert(0, str(__import__("pathlib").Path(__file__).resolve().parent))
from config import TAB, FIG  # noqa
from pseudobulk import parse_age  # noqa


def main(path, name, ct_cols=("cell_type", "predicted.celltype.l2")):
    a = ad.read_h5ad(path, backed="r")
    obs = a.obs
    print(f"=== {name}: {a.shape[0]:,} cells x {a.shape[1]:,} genes ===")
    with h5py.File(path, "r") as f:
        d = f["X/data"][:2_000_000]
        print(f"X sample: dtype={d.dtype} min={d.min()} max={d.max()} integer-fraction={np.mean(np.abs(d-np.round(d))<1e-6):.5f}")
    age = obs["development_stage"].map(parse_age) if "age" not in obs else pd.to_numeric(obs["age"], errors="coerce")
    donors = obs.assign(age_num=age).groupby("donor_id", observed=True).agg(
        n_cells=("age_num", "size"), age=("age_num", "first"), age_nunique=("age_num", "nunique"),
        sex=("sex", "first"), pool=("pool_number", "first") if "pool_number" in obs else ("sex", "first"))
    assert (donors.age_nunique == 1).all(), "a donor has >1 age value"
    n_don = len(donors)
    print(f"donors: {n_don}")
    print(f"age: min={donors.age.min():.0f} max={donors.age.max():.0f} mean={donors.age.mean():.1f} median={donors.age.median():.0f} "
          f"distinct values={donors.age.nunique()}  -> {'CONTINUOUS (per-year)' if donors.age.nunique() > 10 else 'BINNED'}")
    print("age decade distribution (n donors):")
    print(pd.cut(donors.age, bins=range(10, 101, 10)).value_counts().sort_index().to_string())
    print("sex:", donors.sex.value_counts().to_dict())
    print(f"cells per donor: median={donors.n_cells.median():.0f} min={donors.n_cells.min()} max={donors.n_cells.max()}")
    if "pool" in donors:
        print(f"pools: {donors.pool.nunique()}")
    for c in ct_cols:
        if c in obs:
            vc = obs[c].value_counts()
            print(f"\ncell types [{c}]: {len(vc)}")
            print(vc.to_string())
    # sex vs age confound (checked again in 1b)
    m = donors[donors.sex == "male"].age
    f_ = donors[donors.sex == "female"].age
    from scipy import stats
    t = stats.ttest_ind(m, f_, equal_var=False)
    print(f"\nsex-age check: male mean age {m.mean():.1f} (n={len(m)}), female {f_.mean():.1f} (n={len(f_)}); Welch t p={t.pvalue:.3g}")
    donors.to_csv(TAB / f"{name}_donors.csv")
    # STOP 0
    span = donors.age.max() - donors.age.min()
    stop = []
    if n_don < 10:
        stop.append(f"fewer than 10 donors ({n_don})")
    if donors.age.nunique() <= 3:
        stop.append(f"age binned into {donors.age.nunique()} groups")
    if span < 40:
        stop.append(f"age range {span:.0f}y < 40y")
    n_ct = obs[ct_cols[0]].nunique()
    if n_ct < 3:
        stop.append(f"only {n_ct} cell types")
    print("\nSTOP CONDITION 0:", "FIRED -> " + "; ".join(stop) if stop else "NOT fired (>=10 donors, continuous age, span>=40y, >=3 cell types)")
    # age histogram
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    fig, ax = plt.subplots(1, 2, figsize=(10, 3.5))
    ax[0].hist(donors.age, bins=np.arange(15, 101, 2), color="steelblue")
    ax[0].set_xlabel("donor age (years)"); ax[0].set_ylabel("n donors"); ax[0].set_title(f"{name}: {n_don} donors")
    for s, col in (("male", "tab:blue"), ("female", "tab:red")):
        ax[1].hist(donors[donors.sex == s].age, bins=np.arange(15, 101, 4), alpha=0.5, label=s, color=col)
    ax[1].legend(); ax[1].set_xlabel("donor age (years)"); ax[1].set_title("age by sex")
    fig.tight_layout(); fig.savefig(FIG / f"phase0_{name}_donor_age_distribution.png", dpi=130)
    return donors


if __name__ == "__main__":
    main(sys.argv[1], sys.argv[2])
