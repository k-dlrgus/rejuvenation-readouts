"""Task 2 — extrapolation of state-pooled pseudobulks to GTEx fibroblast training. Report-only."""
from __future__ import annotations

import sys
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import sparse

sys.path.insert(0, str(Path(__file__).resolve().parent))
from plane_common import (  # noqa: E402
    PLANE_DIR, MD2_DIR, MD4_DIR, MD4_PROC, PLANE_SEED,
    AGED_LINE, YOUNG_LINE, STATE_NAMES, TRAJECTORY, OFF_TRAJECTORY,
    MIN_CELLS, PCA_MAHAL_K, FROZEN_RULER, IDENTITY_KEY, AGE_INSTRUMENTS,
    FIBRO3_EXTRAP_RHO_NN, FIBRO3_EXTRAP_CI, FIBRO3_EXTRAP_N, MD4_EXTRAP_RHO_NN,
    PREREG_TASK2, PREREG_TASK2_FLAG,
    StopStep, Logger, dump_json, jsonable, load_json, load_frozen_ruler,
    plane_log_banner, load_manifest, save_manifest, log_columns,
    progress_snapshot, flag_true,
)
from fibro2_extrap import _pca_mahal_fit, _nn_and_mahal, _gtex_z  # noqa: E402
from fibro_stage1 import load_pack  # noqa: E402
from fibro2_common import tmm_logcpm_rows, score_frozen  # noqa: E402


def _load_Y(line):
    path = MD4_PROC / f"louvain_rulerY_{line}.npz"
    if not path.exists():
        raise StopStep("task2", f"missing {path}. Not substituting.")
    z = np.load(path)
    return sparse.csr_matrix(
        (z["data"], z["indices"], z["indptr"]),
        shape=tuple(int(x) for x in z["shape"]),
    )


def _direction_from_disk(df, t1r):
    """State the direction extrapolation pushes, from Task 2 distances + repo ρ. Not a gate."""
    by = (t1r or {}).get("by_donor") or {}
    bits = []
    bits.append(
        f"FINDINGS_FIBRO3.md Task 2 (same PCA k={PCA_MAHAL_K} NN / Mahalanobis measures as "
        f"FINDINGS_FIBRO2.md) ρ(nn_euclidean, frozen ruler)={FIBRO3_EXTRAP_RHO_NN} "
        f"CI=[{FIBRO3_EXTRAP_CI[0]}, {FIBRO3_EXTRAP_CI[1]}] n={FIBRO3_EXTRAP_N}. "
        f"FINDINGS_MD4.md ρ_nn={MD4_EXTRAP_RHO_NN} on state×timepoint rows. "
        "Negative ρ: the ruler score falls as distance from GTEx fibroblast training rises."
    )
    for line in (AGED_LINE, YOUNG_LINE):
        g = df[(df.cell_line == line) & (df.scored.astype(str).str.lower() == "true")]
        if not len(g):
            g = df[(df.cell_line == line) & (df["scored"] == True)]  # noqa: E712
        scored = g[g.nn_euclidean.notna()] if len(g) else g
        if not len(scored):
            bits.append(f"{line}: no scored-state distances.")
            continue
        order = scored.sort_values("nn_euclidean")
        seq = "; ".join(
            f"{r.label} nn={float(r.nn_euclidean):+.3f} mahal={float(r.mahalanobis_pca):+.3f}"
            for _, r in order.iterrows()
        )
        bits.append(f"{line} scored states by nn_euclidean (near→far): {seq}.")
        plu = scored[scored.label == "Pluripotency"]
        fib = scored[scored.label == "Fibroblast"]
        if len(plu) and len(fib):
            nn_p = float(plu.nn_euclidean.iloc[0])
            nn_f = float(fib.nn_euclidean.iloc[0])
            far = bool(nn_p > nn_f)
            farthest = str(order.iloc[-1].label)
            bits.append(
                f"{line}: Pluripotency nn={nn_p:+.3f} Fibroblast nn={nn_f:+.3f}; "
                f"Pluripotency farther than Fibroblast={far}; farthest scored state={farthest}."
            )
        ruler_key = ((by.get(line) or {}).get("frozen_ruler") or {}).get("key")
        if ruler_key == "bend":
            bits.append(
                f"{line} frozen ruler reading=`bend`. Extrapolation that pushes age readings "
                "down with distance manufactures lockstep, not a bend, so a bend on the frozen "
                "ruler is conservative with respect to extrapolation."
            )
        elif ruler_key == "lockstep":
            bits.append(
                f"{line} frozen ruler reading=`lockstep`. Extrapolation that pushes age readings "
                "down with distance manufactures lockstep. Lockstep on this instrument is "
                "uninterpretable with respect to extrapolation. This is not used to rescue "
                "the lockstep result."
            )
        else:
            bits.append(
                f"{line} frozen ruler reading=`{ruler_key}`. Extrapolation that pushes age "
                "readings down with distance manufactures lockstep, not a bend. Not used to "
                "reclassify Task 1."
            )
    bits.append(
        "Task 2 is report-only. It is not a gate and is not used to explain away a lockstep result."
    )
    return " ".join(bits)


def run_task2(log=None):
    if not PREREG_TASK2_FLAG.exists():
        raise StopStep("prereg", "PREREG_TASK2.flag missing")
    if not FROZEN_RULER.exists():
        raise StopStep("frozen_ruler", f"missing {FROZEN_RULER}")
    close_log = False
    if log is None:
        log = Logger(PLANE_DIR / "t2_report.txt")
        close_log = True
    plane_log_banner(log, "TASK2")
    log(PREREG_TASK2)
    counts_p = PLANE_DIR / "t1_cell_counts.csv"
    if not counts_p.exists():
        raise StopStep("task2", f"missing {counts_p} — run Task 1")
    counts = pd.read_csv(counts_p)
    log_columns("t1_cell_counts", list(counts.columns), str(counts_p))

    frozen = load_frozen_ruler()
    pack = load_pack()
    log(f"[t2] GTEx pack X shape={np.asarray(pack['X']).shape}")
    Xz = _gtex_z(pack, frozen)
    log(f"[t2] GTEx train z n={Xz.shape[0]} p={Xz.shape[1]}")
    pca = _pca_mahal_fit(Xz, PCA_MAHAL_K, PLANE_SEED)
    log(f"[t2] PCA k={pca['k']}")

    rows = []
    Z_blocks = []
    for line in (AGED_LINE, YOUNG_LINE):
        obs = pd.read_csv(MD2_DIR / f"louvain_obs_{line}.csv")
        lab = pd.read_csv(MD2_DIR / f"t2_cluster_labels_{line}.csv")
        log_columns(f"t2_louvain_obs_{line}", list(obs.columns), str(MD2_DIR / f"louvain_obs_{line}.csv"))
        log_columns(f"t2_labels_{line}", list(lab.columns), str(MD2_DIR / f"t2_cluster_labels_{line}.csv"))
        lab_map = {int(r.cluster): str(r.label) for _, r in lab.iterrows()}
        obs = obs.reset_index(drop=True)
        obs["label"] = obs["cluster"].map(lambda c: lab_map.get(int(c)))
        if obs["label"].isna().any():
            raise StopStep("task2", f"{line}: unlabeled clusters")
        Y = _load_Y(line)
        if Y.shape[0] != len(obs):
            raise StopStep("task2", f"{line} Y n={Y.shape[0]} obs n={len(obs)}")
        mats, meta = [], []
        for st in STATE_NAMES:
            hit = counts[(counts.cell_line == line) & (counts.label == st)]
            if not len(hit):
                raise StopStep("task2", f"t1_cell_counts missing {line} {st}")
            n_all = int(hit.n_all.iloc[0])
            n_obs = int((obs.label == st).sum())
            if n_all != n_obs:
                raise StopStep(
                    "task2",
                    f"occupancy drift {line} {st}: t1 n_all={n_all} obs={n_obs}",
                )
            scored = bool(n_all >= MIN_CELLS)
            if not scored:
                rows.append(dict(
                    cell_line=line, label=st, n_cells=n_all, scored=False,
                    reason=f"n<{MIN_CELLS}", on_trajectory=bool(st in TRAJECTORY),
                    nn_euclidean=np.nan, mahalanobis_pca=np.nan, age_score=np.nan,
                    pca_k=int(pca["k"]),
                ))
                log(f"[t2] {line} {st} n={n_all} not scored (n<{MIN_CELLS})")
                continue
            idx = np.flatnonzero((obs.label == st).to_numpy())
            pb = np.asarray(Y[idx].sum(axis=0), dtype=np.float64).ravel()
            mats.append(pb)
            meta.append(dict(
                cell_line=line, label=st, n_cells=n_all, scored=True, reason=None,
                on_trajectory=bool(st in TRAJECTORY),
            ))
        if not mats:
            log(f"[t2] {line}: no scored-state rows")
            continue
        C = np.vstack(mats)
        logcpm, nf = tmm_logcpm_rows(C, log, tag=f"t2_state_{line}")
        age, Z, missing = score_frozen(logcpm, frozen, counts=C)
        nn, nn_idx, mahal = _nn_and_mahal(Z, Xz, pca)
        for i, rec in enumerate(meta):
            rec.update(
                age_score=float(age[i]),
                nn_euclidean=float(nn[i]),
                nn_train_index=int(nn_idx[i]),
                mahalanobis_pca=float(mahal[i]),
                pca_k=int(pca["k"]),
                age_source="TMM_log2CPM_frozen_ruler_state_pooled_per_donor",
                n_missing_z0=int(np.asarray(missing).sum()),
                tmm_n_rows=int(C.shape[0]),
            )
            rows.append(rec)
            log(f"[t2] {line} {rec['label']} n={rec['n_cells']} "
                f"age={rec['age_score']:+.3f} nn={rec['nn_euclidean']:+.3f} "
                f"mahal={rec['mahalanobis_pca']:+.3f}")
        Z_blocks.append(Z)

    if not rows:
        raise StopStep("task2", "no state rows")
    df = pd.DataFrame(rows)
    df.to_csv(PLANE_DIR / "t2_extrap.csv", index=False)
    if Z_blocks:
        np.savez_compressed(PLANE_DIR / "t2_state_z.npz", Z=np.vstack(Z_blocks))

    t1r = load_json(PLANE_DIR / "t1_reading.json") if (PLANE_DIR / "t1_reading.json").exists() else {}
    statement = _direction_from_disk(df, t1r)
    log(f"[t2 direction] {statement}")

    summary = dict(
        pca_k=int(pca["k"]),
        n_rows=int(len(df)),
        n_scored=int(df.scored.astype(bool).sum()) if "scored" in df.columns else int(len(df)),
        fibro3_rho_nn=FIBRO3_EXTRAP_RHO_NN,
        fibro3_rho_ci=list(FIBRO3_EXTRAP_CI),
        fibro3_n=FIBRO3_EXTRAP_N,
        md4_rho_nn=MD4_EXTRAP_RHO_NN,
        report_only=True,
        used_to_discount_task1=False,
        tmm="among that donor's scored-state rows (not md4 state×timepoint)",
        direction_statement=statement,
        identity=IDENTITY_KEY,
        instruments=list(AGE_INSTRUMENTS),
        NonReprog_scored_but_off_trajectory=True,
        NonReprog_in_trajectory_statistics=False,
    )
    dump_json(PLANE_DIR / "t2_summary.json", jsonable(summary))
    dump_json(PLANE_DIR / "t2_direction.json", jsonable(dict(statement=statement)))
    man = load_manifest()
    man["status"] = "TASK2_DONE"
    save_manifest(man)
    progress_snapshot("Task 3 figures", stop="TASK2_DONE")
    if close_log:
        log.close()
    return summary


if __name__ == "__main__":
    try:
        run_task2()
    except StopStep as e:
        from plane_common import record_failure  # noqa: E402
        record_failure(e.step, e.message, e.details)
        progress_snapshot(f"fix STOP [{e.step}]", stop=f"STOP {e.step}: {e.message}")
        raise
