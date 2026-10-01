"""STAGE P5 — apply FROZEN DLPFC gene sets to the 20-donor white-matter atlas.

Do not re-derive A/I. Dataset ID c05e6940-729c-47bd-a2a6-6ce3730c4919 from the T1 registry.

Usage: python src/brain_phase1_p5.py
"""
from __future__ import annotations

import sys
from pathlib import Path

import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

sys.path.insert(0, str(Path(__file__).resolve().parent))
from brain_phase1_common import (  # noqa: E402
    PHASE1_DIR, PHASE1_FIG, PHASE1_SEED, Logger, HASH_POOL_LIMITATION,
    N_MATCH_DRAWS, N_INNER, ALPHAS, DATA_PROC, dump_json, load_json,
    strip_ensembl, r2_mae, bootstrap_r2_mae,
)
from brain_phase1_models import ridge_predict, size_matched_I_age, age_oof_for_genes, identity_oof_for_genes  # noqa: E402
from brain_phase1_findings import write_findings  # noqa: E402
from decomp import normalize_logcpm  # noqa: E402
from tissue_analyze import make_pool_folds  # noqa: E402
from config import DATA_PROC as _DP  # noqa: E402

WM_PB = DATA_PROC / "tissue_brain_wm_pseudobulk.h5ad"
WM_DATASET_ID = "c05e6940-729c-47bd-a2a6-6ce3730c4919"


def wm_batch_folds(obs, rng, n_outer=5):
    """SequencingPool-grouped folds; assert no donor in both."""
    if "SequencingPool" in obs.columns:
        pool = obs.SequencingPool.astype(str)
    elif "batch" in obs.columns:
        pool = obs.batch.astype(str)
    else:
        raise RuntimeError("WM obs has no SequencingPool")
    pools = pool.unique()
    n_f = int(min(n_outer, max(2, len(pools))))
    fold_of = make_pool_folds(pools, n_f, rng)
    fold = pool.map(fold_of).to_numpy()
    idx = np.arange(len(obs))
    folds = []
    for k in range(n_f):
        te = idx[fold == k]
        tr = idx[fold != k]
        ov = set(obs.iloc[tr].donor.astype(str)) & set(obs.iloc[te].donor.astype(str))
        assert not ov, f"WM donor leak {sorted(ov)[:5]}"
        p_ov = set(pool.iloc[tr]) & set(pool.iloc[te])
        assert not p_ov, f"WM pool leak {p_ov}"
        folds.append((tr, te, k))
    return folds


def map_genes(frozen_ids, wm_genes):
    frozen0 = {strip_ensembl(g) for g in frozen_ids}
    gid = wm_genes.gene_id.astype(str)
    gid0 = gid.map(strip_ensembl)
    hit = np.flatnonzero(gid0.isin(frozen0).to_numpy())
    return hit, int(len(frozen0)), int(len(hit))


def run():
    rng = np.random.default_rng(PHASE1_SEED)
    log = Logger(PHASE1_DIR / "p5_report.txt")
    try:
        return _run(log, rng)
    finally:
        log.close()


def _run(log, rng):
    import anndata as ad
    log("=" * 100)
    log("BRAIN PHASE 1  STAGE P5 — WM replication of FROZEN sets.  seed =", PHASE1_SEED)
    log("=" * 100)
    log(HASH_POOL_LIMITATION)
    log(f"WM dataset_id={WM_DATASET_ID} (T1 registry; not invented). Nothing re-derived.")
    frozen = load_json(PHASE1_DIR / "p2_frozen_gene_sets.json")
    if not WM_PB.exists():
        raise FileNotFoundError(WM_PB)
    pb = ad.read_h5ad(WM_PB)
    log(f"[P5] WM pseudobulk {pb.shape} donors={pb.obs.donor.nunique()} types={pb.obs.celltype.nunique()}")
    Y, gid, sym = normalize_logcpm(pb, log=log)
    genes = pd.DataFrame({"gene_id": gid, "symbol": np.asarray(sym).astype(str)})
    obs = pb.obs.copy().reset_index(drop=True)
    obs["donor"] = obs.donor.astype(str)
    obs["celltype"] = obs.celltype.astype(str)
    obs["age"] = obs.age.astype(float)
    if "SequencingPool" not in obs.columns:
        raise RuntimeError("SequencingPool missing on WM pb")
    obs["Source"] = obs.SequencingPool.astype(str)  # models permute-within-site uses Source
    obs["batch"] = obs.SequencingPool.astype(str)

    A_idx, nA_f, nA_m = map_genes(frozen["A_genes"], genes)
    I_idx, nI_f, nI_m = map_genes(frozen["I_genes"], genes)
    log(f"[P5] mapped A {nA_m}/{nA_f}  I {nI_m}/{nI_f}")

    folds = wm_batch_folds(obs, rng)
    log(f"[P5] SequencingPool-grouped folds={len(folds)}  pools={obs.SequencingPool.nunique()}  "
        f"donors={obs.donor.nunique()}")

    age_A = age_oof_for_genes(Y, obs, A_idx, folds, rng, method="ridge")
    age_I = age_oof_for_genes(Y, obs, I_idx, folds, rng, method="ridge")
    matched = size_matched_I_age(Y, obs, I_idx, n_A=max(len(A_idx), 1), folds=folds,
                                 rng=rng, n_draws=min(N_MATCH_DRAWS, 20), method="ridge")
    # identity on WM types with frozen I vs A (retrained; types differ from DLPFC)
    id_I = identity_oof_for_genes(Y, obs, I_idx, folds, rng)
    id_A = identity_oof_for_genes(Y, obs, A_idx, folds, rng)

    r2_A, r2_I = age_A["median_r2"], matched["median_r2"]
    delta = (r2_A - r2_I) if (np.isfinite(r2_A) and np.isfinite(r2_I)) else np.nan
    # Does the DLPFC cross-test direction hold? A better than I for age.
    holds_dir = bool(np.isfinite(r2_A) and np.isfinite(r2_I) and r2_A > r2_I + 0.05)
    # FALSIFICATION support band is stricter; n=20 is small — report plainly.
    p4 = load_json(PHASE1_DIR / "p4_summary.json") if (PHASE1_DIR / "p4_summary.json").exists() else {}
    dlpfc_v = p4.get("verdict", "")
    if nA_m < 10:
        text = f"Too few A-genes mapped ({nA_m}) to test replication."
        holds = False
        recommend_retina = True
    elif not np.isfinite(r2_A) or r2_A < 0.05:
        text = (
            f"Frozen A-genes do **not** retain a clear age association on WM "
            f"(median-over-types R²={r2_A:+.3f}, n_donors=20). Cross-test cannot confirm DLPFC."
        )
        holds = False
        recommend_retina = True
    elif holds_dir:
        text = (
            f"A-genes retain age association on WM (R²_A={r2_A:+.3f}) and beat size-matched I-genes "
            f"(R²_I={r2_I:+.3f}, Δ={delta:+.3f}). Direction of the DLPFC cross-test holds. "
            f"n=20 / 7 pools: CIs are wide; this is not a powered copy of the 233-donor result. "
            f"DLPFC verdict was: {dlpfc_v}."
        )
        holds = True
        recommend_retina = False
    else:
        text = (
            f"WM R²_A={r2_A:+.3f} vs R²_I(matched)={r2_I:+.3f} (Δ={delta:+.3f}). "
            f"The DLPFC cross-test verdict does **not** cleanly replicate at n=20. "
            f"DLPFC verdict was: {dlpfc_v}."
        )
        holds = False
        recommend_retina = True
    log("[P5] " + text)
    age_A["per_type"].assign(which="A").to_csv(PHASE1_DIR / "p5_wm_age_per_type.csv", index=False)
    dump_json(PHASE1_DIR / "p5_summary.json", dict(
        dataset_id=WM_DATASET_ID, n_donors=int(obs.donor.nunique()),
        n_types=int(obs.celltype.nunique()), n_pb=int(len(obs)),
        n_A_frozen=nA_f, n_A_mapped=nA_m, n_I_frozen=nI_f, n_I_mapped=nI_m,
        r2_A=r2_A, r2_I_matched=r2_I, r2_I_full=age_I["median_r2"],
        acc_I=id_I["acc"], acc_A=id_A["acc"],
        joint_r2_A=age_A["joint_r2"],
        cross_test_holds=holds, recommend_retina=recommend_retina,
        text=text, seed=PHASE1_SEED, dlpfc_verdict=dlpfc_v,
    ))

    fig, ax = plt.subplots(figsize=(5.6, 3.6))
    ax.bar(["A-genes", "I matched", "I full"], [r2_A, r2_I, age_I["median_r2"]],
           color=["tab:red", "tab:blue", "tab:cyan"])
    ax.axhline(0, c="k", lw=0.5)
    ax.set_ylabel("WM median-over-types age R²")
    ax.set_title("P5 frozen DLPFC genes on white matter (n=20)")
    fig.tight_layout()
    fig.savefig(PHASE1_FIG / "p5_wm_cross_test.png", dpi=140)
    plt.close(fig)

    write_findings(log=log)
    log("[P5] done.")
    return dict(holds=holds)


if __name__ == "__main__":
    run()
