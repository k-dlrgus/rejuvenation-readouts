"""Freeze DLPFC raw ridge and identity-residualized ridge on all 233 donors.

Must not run unless results/external/PREREG_20260917.flag exists (written after
Stage 0, before this file computes a coefficient).
"""
from __future__ import annotations

import sys
import time
from pathlib import Path

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parent))
from external_common import (  # noqa: E402
    EXT_DIR, EXT_SEED, PREREG_FLAG, StopStep, Logger, dump_json, jsonable,
    ext_log_banner, load_dlpfc_types, refresh_progress,
)
from target_common import per_type_directions, load_phase1_matrix  # noqa: E402
from trajectory_common import zscore_train, calibrate_1d, residualize_identity  # noqa: E402
from geometry_partC import identity_basis  # noqa: E402


def assert_prereg():
    if not PREREG_FLAG.exists():
        raise StopStep(
            "prereg",
            f"{PREREG_FLAG.name} missing — write the pre-registration flag after Stage 0 "
            "and before freezing any direction.",
            dict(path=str(PREREG_FLAG)),
        )


def freeze():
    assert_prereg()
    log = Logger(EXT_DIR / "freeze_report.txt")
    try:
        return _freeze(log)
    finally:
        log.close()


def _freeze(log):
    ext_log_banner(log, "freeze DLPFC directions on all 233 donors")
    types, _ = load_dlpfc_types(log)
    data = load_phase1_matrix(log)
    Y, obs = data["Y"], data["obs"]
    genes = data["genes"]
    ct = obs.celltype.astype(str).to_numpy()
    y = obs.age.to_numpy(float)
    types_obs = sorted(pd.unique(ct))
    if types_obs != types:
        raise StopStep(
            "freeze",
            f"Phase-1 matrix types {types_obs} != obs-csv types {types}. Not reordering silently.",
            dict(matrix=types_obs, csv=types),
        )
    t0 = time.time()
    Xs, mu, sd = zscore_train(Y)
    log(f"[freeze] global z-score  rows={Xs.shape[0]} genes={Xs.shape[1]}")
    W_raw, meta_raw, _ = per_type_directions(Xs, y, ct, types, method="ridge")
    log("[freeze] raw ridge per type:")
    log(meta_raw.to_string(index=False))
    basis, S_id, used = identity_basis(Xs, ct, types)
    log(f"[freeze] identity basis rank={basis.shape[1]}  singular={np.asarray(S_id)[:5].tolist()}  used_types={len(used)}")
    R = residualize_identity(Xs, basis)
    W_id, meta_id, _ = per_type_directions(R, y, ct, types, method="ridge")
    log("[freeze] identity-residual ridge per type:")
    log(meta_id.to_string(index=False))

    cal_raw, cal_id = [], []
    for j, t in enumerate(types):
        m = ct == t
        rec_r = dict(celltype=t, n=int(m.sum()), intercept=np.nan, slope=np.nan, ok=False)
        rec_i = dict(celltype=t, n=int(m.sum()), intercept=np.nan, slope=np.nan, ok=False)
        w = W_raw[:, j]
        if np.isfinite(w).all() and float(np.linalg.norm(w)) > 1e-12 and int(m.sum()) >= 8:
            _, a, b = calibrate_1d(Xs[m] @ w, y[m], Xs[m] @ w)
            rec_r.update(intercept=float(a), slope=float(b), ok=True)
        w = W_id[:, j]
        if np.isfinite(w).all() and float(np.linalg.norm(w)) > 1e-12 and int(m.sum()) >= 8:
            _, a, b = calibrate_1d(R[m] @ w, y[m], R[m] @ w)
            rec_i.update(intercept=float(a), slope=float(b), ok=True)
        cal_raw.append(rec_r)
        cal_id.append(rec_i)
    cal_raw = pd.DataFrame(cal_raw)
    cal_id = pd.DataFrame(cal_id)
    cal_raw.to_csv(EXT_DIR / "freeze_raw_calibration.csv", index=False)
    cal_id.to_csv(EXT_DIR / "freeze_idresid_calibration.csv", index=False)
    meta_raw.to_csv(EXT_DIR / "freeze_raw_meta.csv", index=False)
    meta_id.to_csv(EXT_DIR / "freeze_idresid_meta.csv", index=False)

    out = EXT_DIR / "frozen_directions.npz"
    np.savez_compressed(
        out,
        gene_id=genes.gene_id.astype(str).to_numpy(),
        symbol=genes.symbol.astype(str).to_numpy(),
        types=np.array(types, dtype=object),
        mu=np.asarray(mu, np.float64),
        sd=np.asarray(sd, np.float64),
        W_raw=np.asarray(W_raw, np.float64),
        W_idresid=np.asarray(W_id, np.float64),
        identity_basis=np.asarray(basis, np.float64),
        identity_S=np.asarray(S_id, np.float64),
        cal_raw_intercept=cal_raw.intercept.to_numpy(float),
        cal_raw_slope=cal_raw.slope.to_numpy(float),
        cal_id_intercept=cal_id.intercept.to_numpy(float),
        cal_id_slope=cal_id.slope.to_numpy(float),
        seed=np.array(EXT_SEED),
        n_donors=np.array(int(obs.donor.nunique())),
        n_id_axes=np.array(int(basis.shape[1])),
    )
    rec = dict(
        path=str(out),
        n_donors=int(obs.donor.nunique()),
        n_genes=int(Y.shape[1]),
        n_rows=int(Y.shape[0]),
        n_types=len(types),
        n_id_axes=int(basis.shape[1]),
        n_raw_ok=int(cal_raw.ok.sum()),
        n_id_ok=int(cal_id.ok.sum()),
        elapsed_s=float(time.time() - t0),
        seed=EXT_SEED,
        note="trained on all 233 DLPFC donors; not fit on SEA-AD",
    )
    dump_json(EXT_DIR / "freeze.json", jsonable(rec))
    log(f"[freeze] wrote {out} in {rec['elapsed_s']:.0f}s  raw_ok={rec['n_raw_ok']} id_ok={rec['n_id_ok']}")
    refresh_progress(next_action="project SEA-AD pseudobulks; primary + secondary cells")
    return rec


if __name__ == "__main__":
    freeze()
