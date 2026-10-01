"""Task 2 — is NonReprog the odd population? Report-only. No gate."""
from __future__ import annotations

import sys
from pathlib import Path

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parent))
from md5_common import (  # noqa: E402
    MD5_DIR, MD2_DIR, MD2_PROC, MD4_DIR, MD4_PROC,
    AGED_LINE, YOUNG_LINE, DAYS, INSTRUMENTS,
    TASK2_PAIR, TASK2_STATES, QC_COLS, CELL_CYCLE_COLS, MIN_CELLS_A,
    PREREG_TASK2, PREREG_TASK2_FLAG, FROZEN_RULER, N_PERM, N_BOOT,
    PCA_MAHAL_K, PAPER_NONREPROG_QUOTE,
    perm_seed_t2, boot_seed_t2,
    StopStep, Logger, dump_json, jsonable, load_json, md5_log_banner,
    load_manifest, save_manifest, record_failure, record_unusable, log_columns,
    progress_snapshot,
)
from md4_task1 import _one_contrast  # noqa: E402
from md5_task1 import _load_labels, _load_obs_labelled, load_md4_cell_scores  # noqa: E402


def _cell_cycle_from_md2(obs, ams_files, log, line):
    obs_hits = [c for c in CELL_CYCLE_COLS if c in obs.columns]
    ams_hits = [c for c in CELL_CYCLE_COLS if c in ams_files]
    rec = dict(
        cell_line=line,
        obs_columns=list(obs.columns),
        ams_keys=list(ams_files),
        obs_cell_cycle_hits=obs_hits,
        ams_cell_cycle_hits=ams_hits,
        computable=bool(obs_hits or ams_hits),
    )
    if not rec["computable"]:
        rec["reason"] = (
            "md2 louvain_obs / louvain_ams have no phase / S.Score / G2M.Score columns. "
            "Cell-cycle phase fractions are not computable from the md2 object. "
            "Not substituting Seurat CellCycleScoring or a new Tirosh score."
        )
        log(f"[t2 cellcycle] {line} not computable. {rec['reason']}")
        record_unusable("cell_cycle", rec["reason"], rec)
    else:
        rec["reason"] = None
        log(f"[t2 cellcycle] {line} hits obs={obs_hits} ams={ams_hits}")
    return rec


def run_task2(log=None):
    if not PREREG_TASK2_FLAG.exists():
        raise StopStep("prereg", "PREREG_TASK2.flag missing")
    if not FROZEN_RULER.exists():
        raise StopStep("frozen_ruler", f"missing {FROZEN_RULER}")
    close_log = False
    if log is None:
        log = Logger(MD5_DIR / "t2_report.txt")
        close_log = True
    md5_log_banner(log, "TASK2")
    log(PREREG_TASK2)

    extrap_p = MD4_DIR / "t2_extrap.csv"
    if not extrap_p.exists():
        raise StopStep("task2", f"missing {extrap_p}. Not recomputing distances.")
    extrap_all = pd.read_csv(extrap_p)
    log_columns("md4_t2_extrap", list(extrap_all.columns), str(extrap_p))
    need_ex = {"cell_line", "day", "label", "n_cells", "nn_euclidean", "mahalanobis_pca", "pca_k"}
    missing_ex = sorted(need_ex - set(extrap_all.columns))
    if missing_ex:
        raise StopStep("task2", f"md4 t2_extrap.csv missing {missing_ex}. Not substituting.")

    qc_rows = []
    extrap_rows = []
    contrast_rows = []
    cc_meta = {}
    scores_by, obs_by = {}, {}

    for line in (AGED_LINE, YOUNG_LINE):
        lab = _load_labels(line, log)
        obs = _load_obs_labelled(line, lab, log)
        obs_by[line] = obs
        ams_path = MD2_PROC / f"louvain_ams_{line}.npz"
        ams = np.load(ams_path, allow_pickle=True)
        log_columns(f"t2_louvain_ams_{line}", list(ams.files), str(ams_path))
        cc_meta[line] = _cell_cycle_from_md2(obs, list(ams.files), log, line)
        miss_qc = [c for c in QC_COLS if c not in obs.columns]
        if miss_qc:
            raise StopStep("task2", f"louvain_obs_{line} missing QC columns {miss_qc}. Not substituting.")
        scores_by[line] = load_md4_cell_scores(line, obs, log)

        for d in DAYS:
            for st in TASK2_STATES:
                m = (obs.day.astype(int) == int(d)) & (obs.label == st)
                n = int(m.sum())
                rec = dict(cell_line=line, day=int(d), label=st, n_cells=n, min_cells=MIN_CELLS_A)
                if n < MIN_CELLS_A:
                    rec.update(
                        ok=False, reason=f"n<{MIN_CELLS_A}",
                        median_umi=np.nan, median_genes=np.nan, median_mito_frac=np.nan,
                        cell_cycle_computable=bool(cc_meta[line]["computable"]),
                        cell_cycle_note=cc_meta[line]["reason"],
                    )
                    qc_rows.append(rec)
                    continue
                sub = obs.loc[m]
                rec.update(
                    ok=True, reason=None,
                    median_umi=float(np.median(sub["umi"].to_numpy(float))),
                    median_genes=float(np.median(sub["n_genes"].to_numpy(float))),
                    median_mito_frac=float(np.median(sub["mito_frac"].to_numpy(float))),
                    cell_cycle_computable=bool(cc_meta[line]["computable"]),
                    cell_cycle_note=cc_meta[line]["reason"],
                )
                if cc_meta[line]["computable"]:
                    phase_col = None
                    for c in ("phase", "Phase"):
                        if c in sub.columns:
                            phase_col = c
                            break
                    if phase_col is not None:
                        vc = sub[phase_col].astype(str).value_counts(dropna=False)
                        rec["phase_fractions"] = {str(k): float(v / n) for k, v in vc.items()}
                qc_rows.append(rec)
                log(f"[t2 qc] {line} d{d} {st} n={n} median_UMI={rec['median_umi']:.1f} "
                    f"median_genes={rec['median_genes']:.1f} mito={rec['median_mito_frac']:.4f}")

                hit = extrap_all[
                    (extrap_all.cell_line == line)
                    & (extrap_all.day.astype(int) == int(d))
                    & (extrap_all.label == st)
                ]
                if not len(hit):
                    raise StopStep(
                        "task2",
                        f"md4 t2_extrap.csv missing {line} d{d} {st}. Not recomputing.",
                    )
                er = hit.iloc[0]
                n_ex = int(er.n_cells)
                if n_ex != n:
                    raise StopStep(
                        "task2",
                        f"extrap n_cells drift {line} d{d} {st}: md4={n_ex} md5={n}",
                    )
                pca_k = int(er.pca_k)
                if pca_k != int(PCA_MAHAL_K):
                    raise StopStep("task2", f"{line} d{d} {st} pca_k={pca_k} != {PCA_MAHAL_K}")
                extrap_rows.append(dict(
                    cell_line=line, day=int(d), label=st, n_cells=n,
                    age_score=float(er.age_score) if "age_score" in er.index else np.nan,
                    nn_euclidean=float(er.nn_euclidean),
                    mahalanobis_pca=float(er.mahalanobis_pca),
                    pca_k=pca_k,
                    age_source=str(er.age_source) if "age_source" in er.index else None,
                    source=str(extrap_p),
                    recomputed=False,
                ))

        for d in DAYS:
            day_m = obs.day.astype(int) == int(d)
            a, b = TASK2_PAIR
            n_a = int((day_m & (obs.label == a)).sum())
            n_b = int((day_m & (obs.label == b)).sum())
            qualifies = bool(n_a >= MIN_CELLS_A and n_b >= MIN_CELLS_A)
            base = dict(
                cell_line=line, day=int(d), state_a=a, state_b=b,
                n_a=n_a, n_b=n_b, qualifies=qualifies, min_cells=MIN_CELLS_A,
                permutation="md4_task1._one_contrast_unstratified_within_timepoint",
            )
            if not qualifies:
                for inst in INSTRUMENTS:
                    contrast_rows.append(dict(
                        **base, instrument=inst, ok=False,
                        reason=f"n<{MIN_CELLS_A} in at least one state",
                    ))
                continue
            sub = day_m & obs.label.isin([a, b])
            idx = np.flatnonzero(sub.to_numpy())
            labels_sub = obs.loc[sub, "label"].to_numpy()
            mask_a = labels_sub == a
            sc = scores_by[line]
            for inst in INSTRUMENTS:
                x = np.asarray(sc[inst], float)[idx]
                rec = _one_contrast(
                    x, mask_a, N_PERM, N_BOOT,
                    perm_seed_t2(line, d), boot_seed_t2(line, d),
                )
                rec.update(base)
                rec["instrument"] = inst
                rec["nr_older"] = bool(
                    rec.get("ok") and rec.get("delta_mean") is not None
                    and float(rec["delta_mean"]) > 0
                    and rec.get("p_mean_two_sided") is not None
                    and float(rec["p_mean_two_sided"]) <= 0.05
                )
                contrast_rows.append(rec)
                log(f"[t2 NR-Fib] {line} d{d} {inst} Δmean={rec.get('delta_mean')} "
                    f"p_two={rec.get('p_mean_two_sided')} nr_older={rec.get('nr_older')}")

    qc_df = pd.DataFrame(qc_rows)
    qc_df.to_csv(MD5_DIR / "t2_qc.csv", index=False)
    ex_df = pd.DataFrame(extrap_rows)
    ex_df.to_csv(MD5_DIR / "t2_extrap.csv", index=False)
    cdf = pd.DataFrame(contrast_rows)
    cdf.to_csv(MD5_DIR / "t2_contrasts_nr_fib.csv", index=False)

    dump_json(MD5_DIR / "t2_cellcycle.json", jsonable(cc_meta))

    q = cdf[(cdf.qualifies == True) & (cdf.ok == True)]  # noqa: E712
    obs_rows = []
    consistently = {}
    for inst in INSTRUMENTS:
        sub = q[q.instrument == inst]
        n_q = int(len(sub))
        n_older = int(sub["nr_older"].sum()) if n_q and "nr_older" in sub.columns else 0
        consistently[inst] = bool(n_q > 0 and n_older == n_q)
        obs_rows.append(dict(
            instrument=inst, n_qualifying=n_q, n_nr_older=n_older,
            consistently_older=consistently[inst],
        ))
        log(f"[t2 obs] {inst} n_qualifying={n_q} n_nr_older={n_older} consistently={consistently[inst]}")
    obs_df = pd.DataFrame(obs_rows)
    obs_df.to_csv(MD5_DIR / "t2_nr_older_summary.csv", index=False)

    ruler_cons = bool(consistently.get("frozen_ruler"))
    md_cons = bool(consistently.get("md_score"))
    note = None
    if ruler_cons:
        note = (
            "NonReprog reads consistently older than Fibroblast on the frozen ruler "
            f"at every qualifying donor × timepoint (n={int(obs_df.loc[obs_df.instrument=='frozen_ruler','n_qualifying'].iloc[0])}). "
            "Lu et al. (Cell 2025, 188:5895–5911) on the non-reprogrammed state: "
            f"\"{PAPER_NONREPROG_QUOTE}\". "
            "Observation only. Not used to discount Task 1."
        )
    else:
        note = (
            "NonReprog does not read consistently older than Fibroblast on the frozen ruler "
            "across all qualifying donor × timepoints. Observation only. Not used to discount Task 1. "
            f"Paper remark (not used as a gate): \"{PAPER_NONREPROG_QUOTE}\"."
        )
    if md_cons:
        note += " MD likewise has NonReprog older than Fibroblast at every qualifying timepoint."

    summary = dict(
        n_qc_rows=int(len(qc_df)),
        n_qc_ok=int(qc_df.ok.sum()) if len(qc_df) else 0,
        n_extrap_rows=int(len(ex_df)),
        n_contrast_rows=int(len(cdf)),
        pca_k=int(PCA_MAHAL_K),
        cell_cycle_computable={ln: bool(cc_meta[ln]["computable"]) for ln in cc_meta},
        consistently_older=consistently,
        paper_nonreprog_quote=PAPER_NONREPROG_QUOTE,
        note=note,
        report_only=True,
        used_to_discount_task1=False,
    )
    dump_json(MD5_DIR / "t2_summary.json", jsonable(summary))
    (MD5_DIR / "t2_note.txt").write_text(note, encoding="utf-8")
    log(f"[t2] {note}")
    man = load_manifest()
    man["status"] = "TASK2_DONE"
    save_manifest(man)
    progress_snapshot("write FINDINGS_MD5.md", stop="TASK2_DONE")
    if close_log:
        log.close()
    return summary


if __name__ == "__main__":
    try:
        run_task2()
    except StopStep as e:
        record_failure(e.step, e.message, e.details)
        progress_snapshot(f"fix STOP [{e.step}]", stop=f"STOP {e.step}: {e.message}")
        raise
