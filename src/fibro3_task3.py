"""Task 3 — per-cluster random-direction null on d0→d7 and d0→d10. d7 is not the Stage 2 endpoint."""
from __future__ import annotations

import sys
from pathlib import Path

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parent))
from fibro3_common import (  # noqa: E402
    FIBRO3_DIR, FIBRO3_SEED, N_RANDOM_DIR, MIN_CELLS, AGED_LINE, YOUNG_LINE,
    PREREG_TASK3, FROZEN_RULER,
    StopStep, Logger, dump_json, jsonable, fibro3_log_banner,
    load_manifest, save_manifest, record_failure, log_columns,
    load_frozen_ruler, progress_snapshot,
)
from fibro_stage2 import random_directions  # noqa: E402
from trajectory_common import permutation_p  # noqa: E402


def _pick_row(meta, line, cl, day):
    m = (
        (meta.cell_line.astype(str).str.upper() == str(line).upper())
        & (meta.cluster.astype(int) == int(cl))
        & (meta.day.astype(int) == int(day))
    )
    idx = np.flatnonzero(m.to_numpy())
    if len(idx) == 0:
        return None
    if len(idx) != 1:
        raise StopStep("task3", f"{line} cluster={cl} day={day} n_rows={len(idx)} (want 0 or 1)")
    return int(idx[0])


def run_task3(log=None):
    if not (FIBRO3_DIR / "PREREG_TASK3.flag").exists():
        raise StopStep("prereg", "PREREG_TASK3.flag missing")
    close_log = False
    if log is None:
        log = Logger(FIBRO3_DIR / "t3_report.txt")
        close_log = True
    fibro3_log_banner(log, "TASK3")
    log(PREREG_TASK3)
    zpath = FIBRO3_DIR / "t1_cluster_z.npz"
    tpath = FIBRO3_DIR / "t1_cluster_table.csv"
    if not zpath.exists():
        raise StopStep("task3", f"missing {zpath} — run Task 1")
    if not tpath.exists():
        raise StopStep("task3", f"missing {tpath}")
    if not FROZEN_RULER.exists():
        raise StopStep("frozen_ruler", f"missing {FROZEN_RULER}")
    frozen = load_frozen_ruler()
    z = np.load(zpath, allow_pickle=True)
    log(f"[task3] t1_cluster_z keys actually read: {list(z.files)}")
    log_columns("t3_t1_z", list(z.files), str(zpath))
    need = {"Z", "age"}
    missing = sorted(need - set(z.files))
    if missing:
        raise StopStep("task3", f"{zpath.name} missing {missing}")
    Z = np.asarray(z["Z"], np.float64)
    age = np.asarray(z["age"], np.float64)
    meta = pd.read_csv(tpath)
    log(f"[task3] t1_cluster_table columns actually read: {list(meta.columns)}")
    log_columns("t3_t1_table", list(meta.columns), str(tpath))
    if len(meta) != len(Z) or len(meta) != len(age):
        raise StopStep("task3", f"length mismatch table={len(meta)} Z={len(Z)} age={len(age)}")
    if "n_cells" not in meta.columns:
        raise StopStep("task3", "t1_cluster_table missing n_cells")

    null = random_directions(frozen["w"], Z, seed=FIBRO3_SEED, n=N_RANDOM_DIR)
    log(f"[task3] random directions n={N_RANDOM_DIR} seed={FIBRO3_SEED} Z={Z.shape}")
    np.savez_compressed(FIBRO3_DIR / "t3_random_dir_scores.npz", null=null)

    rows = []
    for line in (AGED_LINE, YOUNG_LINE):
        sub = meta[meta.cell_line.astype(str).str.upper() == line]
        clusters = sorted(sub.cluster.astype(int).unique().tolist())
        for cl in clusters:
            for dend in (7, 10):
                i0 = _pick_row(meta, line, cl, 0)
                i1 = _pick_row(meta, line, cl, dend)
                rec = dict(
                    cell_line=line, cluster=int(cl), endpoint=f"d0→d{dend}", d_end=int(dend),
                    primary=bool(line == AGED_LINE), n_random=N_RANDOM_DIR, seed=FIBRO3_SEED,
                    min_cells=MIN_CELLS,
                )
                if i0 is None or i1 is None:
                    rec.update(ok=False, reason="missing cluster row at an endpoint",
                               n_cells_d0=np.nan if i0 is None else int(meta.n_cells.iloc[i0]),
                               n_cells_end=np.nan if i1 is None else int(meta.n_cells.iloc[i1]))
                    rows.append(rec)
                    log(f"[task3] skip {line} c{cl} d0→d{dend}: missing row")
                    continue
                n0 = int(meta.n_cells.iloc[i0])
                n1 = int(meta.n_cells.iloc[i1])
                rec["n_cells_d0"] = n0
                rec["n_cells_end"] = n1
                if n0 < MIN_CELLS or n1 < MIN_CELLS:
                    rec.update(ok=False, reason=f"n_cells<{MIN_CELLS} at an endpoint")
                    rows.append(rec)
                    log(f"[task3] skip {line} c{cl} d0→d{dend}: n0={n0} n1={n1} < {MIN_CELLS}")
                    continue
                real = float(age[i0] - age[i1])
                null_d = np.asarray(null)[i0] - np.asarray(null)[i1]
                p = permutation_p(real, null_d, greater=True)
                n_ge = int(np.sum(null_d >= real))
                rec.update(
                    ok=True, score_d0=float(age[i0]), score_end=float(age[i1]),
                    decline=real, p=p, n_random_ge_real=n_ge,
                    pass_p=bool(np.isfinite(p) and p <= 0.05),
                )
                rows.append(rec)
                log(f"[task3] {line} c{cl} d0→d{dend} decline={real:+.4f} "
                    f"p={p:+.3f} n_ge={n_ge}/{N_RANDOM_DIR} n0={n0} n1={n1}")

    df = pd.DataFrame(rows)
    df.to_csv(FIBRO3_DIR / "t3_null.csv", index=False)
    aged_ok = df[(df.cell_line == AGED_LINE) & (df.ok == True)]
    summary = dict(
        n_rows=int(len(df)), n_ok=int(df.ok.sum()) if "ok" in df.columns else 0,
        min_cells=MIN_CELLS, n_random=N_RANDOM_DIR, seed=FIBRO3_SEED,
        note="joint cluster IDs pair across days; d7 is not the Stage 2 endpoint; not averaged",
        aged_ok=aged_ok.to_dict(orient="records") if len(aged_ok) else [],
    )
    dump_json(FIBRO3_DIR / "t3_summary.json", jsonable(summary))
    man = load_manifest()
    man["status"] = "TASK3_DONE"
    save_manifest(man)
    progress_snapshot("write FINDINGS_FIBRO3.md from disk", stop="TASK3_DONE")
    if close_log:
        log.close()
    return summary


if __name__ == "__main__":
    try:
        run_task3()
    except StopStep as e:
        record_failure(e.step, e.message, e.details)
        progress_snapshot(f"fix STOP [{e.step}]", stop=f"STOP {e.step}: {e.message}")
        raise
