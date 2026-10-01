"""Task 1 — occupancy, sanity, state medians, segment slopes, bend statistic."""
from __future__ import annotations

import sys
from pathlib import Path

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parent))
from plane_common import (  # noqa: E402
    PLANE_DIR, MD2_DIR, MD2_PROC, MD3_DIR, MD4_DIR, MD4_PROC,
    PLANE_SEED, PLANE_BOOT, AGED_LINE, YOUNG_LINE, DAYS, STATE_NAMES,
    TRAJECTORY, OFF_TRAJECTORY, SEGMENT_A, SEGMENT_B_PREFERRED_END,
    SEGMENT_B_FALLBACK_END, IDENTITY_KEY, AGE_INSTRUMENTS, MIN_CELLS,
    N_BOOT, N_RANDOM_DIR, FROZEN_RULER, PLURI_ENDOGENOUS, FIBRO_IDENTITY,
    MD3_MD_VS_PLURI_GM00731, MD3_MD_VS_PLURI_GM23815,
    PREREG_TASK1, PREREG_TASK1_FLAG, DECLARED_BEFORE_SCORES_FLAG,
    MD2_POOLED_MD_PR, MD2_POOLED_MD_NR, MD2_N_PR, MD2_N_NR,
    MMC3_XLSX, StopStep, Logger, dump_json, jsonable, load_json,
    plane_log_banner, load_manifest, save_manifest, log_columns, log_geneset,
    progress_snapshot, detect_id_type, require_mappable, boot_seed,
    ci_includes_0, cis_overlap, slope_young_side, ci_excludes_0_on_side,
    _slope, segment_B_end, flag_true,
)
from md3_idtype import csr_from_npz  # noqa: E402


def _require_prereg():
    if not PREREG_TASK1_FLAG.exists():
        raise StopStep("prereg", "PREREG_TASK1.flag missing")
    if not FROZEN_RULER.exists():
        raise StopStep("frozen_ruler", f"missing {FROZEN_RULER}")
    for line in (AGED_LINE, YOUNG_LINE):
        for p in (
            MD2_DIR / f"t2_cluster_labels_{line}.csv",
            MD2_DIR / f"louvain_obs_{line}.csv",
            MD2_PROC / f"louvain_counts_{line}.npz",
            MD2_PROC / f"louvain_ams_{line}.npz",
            MD4_PROC / f"cell_scores_{line}.npz",
            MD4_PROC / f"louvain_rulerY_{line}.npz",
        ):
            if not p.exists():
                raise StopStep("task1", f"missing artifact {p}. Not substituting.")
    for p in (
        MD2_DIR / "t2_reproduction.json",
        MD2_DIR / "genesets.json",
        MD3_DIR / "genesets.json",
        MD4_DIR / "t1_cell_counts.csv",
        MD4_DIR / "genesets.json",
        MD4_DIR / "genesets_summary.json",
    ):
        if not p.exists():
            raise StopStep("task1", f"missing artifact {p}. Not substituting.")


def load_genesets(log):
    md2_gs = load_json(MD2_DIR / "genesets.json")
    md3_gs = load_json(MD3_DIR / "genesets.json")
    md4_gs = load_json(MD4_DIR / "genesets.json")
    log_columns("md2_genesets_json", list(md2_gs.keys()), str(MD2_DIR / "genesets.json"))
    log_columns("md3_genesets_json", list(md3_gs.keys()), str(MD3_DIR / "genesets.json"))
    log_columns("md4_genesets_json", list(md4_gs.keys()), str(MD4_DIR / "genesets.json"))
    md = [str(g).upper() for g in md2_gs["MD"]]
    age_up = [str(g).upper() for g in md4_gs["age_up"]]
    age_down = [str(g).upper() for g in md4_gs["age_down"]]
    md3_md = [str(g).upper() for g in md3_gs["MD"]]
    md3_up = [str(g).upper() for g in md3_gs["age_up"]]
    md3_dn = [str(g).upper() for g in md3_gs["age_down"]]
    if set(md) != set(md3_md):
        raise StopStep("genesets", "md2 vs md3 MD sets differ. Not reconciling.")
    if set(age_up) != set(md3_up) or set(age_down) != set(md3_dn):
        raise StopStep("genesets", "md4 vs md3 Age up/down differ. Not reconciling.")
    log_geneset("MD_built_from_md2", "MD (STAR Methods built; md2)", md,
                str(MD2_DIR / "genesets.json"), id_type="symbol")
    log_geneset("mmc3_Aging_signatures_Age up", "Aging_signatures:Age up", age_up,
                str(MMC3_XLSX) + " via md4/genesets.json", id_type="symbol")
    log_geneset("mmc3_Aging_signatures_Age down", "Aging_signatures:Age down", age_down,
                str(MMC3_XLSX) + " via md4/genesets.json", id_type="symbol")
    log_geneset("PLURI_ENDOGENOUS", "FINDINGS_FIBRO.md PLURI_ENDOGENOUS",
                list(PLURI_ENDOGENOUS), "fibro_common.PLURI_ENDOGENOUS", id_type="symbol")
    log_geneset("FIBRO_IDENTITY", "FINDINGS_FIBRO.md FIBRO_IDENTITY",
                list(FIBRO_IDENTITY), "fibro_common.FIBRO_IDENTITY", id_type="symbol")
    md4_sum = load_json(MD4_DIR / "genesets_summary.json")
    summary = dict(
        n_MD=len(md), n_age_up=len(age_up), n_age_down=len(age_down),
        md4_summary=md4_sum,
        used_as_deseq2_rebuild=False, used_as_published_lists=True,
        seed=PLANE_SEED,
        source_md2=str(MD2_DIR / "genesets.json"),
        source_md3=str(MD3_DIR / "genesets.json"),
        source_md4=str(MD4_DIR / "genesets.json"),
        pluri_endogenous=list(PLURI_ENDOGENOUS),
        fibro_identity=list(FIBRO_IDENTITY),
        identity_version="md4_pluri_primary_drop_oskm_False",
    )
    s3_meta = {}
    for line in (AGED_LINE, YOUNG_LINE):
        p = MD4_DIR / f"t1_s3_ams_meta_{line}.json"
        if not p.exists():
            raise StopStep("genesets", f"missing {p}. Not substituting mapped counts.")
        rec = load_json(p)
        s3_meta[line] = rec
        log_columns(f"t1_s3_ams_meta_{line}", list(rec.keys()), str(p))
    summary["s3_mapped"] = s3_meta
    dump_json(PLANE_DIR / "genesets_summary.json", jsonable(summary))
    dump_json(PLANE_DIR / "genesets.json", jsonable(dict(
        MD=md, age_up=age_up, age_down=age_down,
    )))
    log(f"[genesets] MD n={len(md)} Age up n={len(age_up)} Age down n={len(age_down)}")
    return dict(MD=md, age_up=age_up, age_down=age_down, s3_meta=s3_meta)


def _load_labels(line, log):
    path = MD2_DIR / f"t2_cluster_labels_{line}.csv"
    lab = pd.read_csv(path)
    log(f"[t1 labels] {line} columns actually read: {list(lab.columns)}")
    log_columns(f"t2_cluster_labels_{line}", list(lab.columns), str(path))
    need = {"cluster", "label", "n_cells"}
    missing = sorted(need - set(lab.columns))
    if missing:
        raise StopStep("task1", f"t2_cluster_labels_{line}.csv missing {missing}")
    return lab


def _load_obs_labelled(line, lab, log):
    path = MD2_DIR / f"louvain_obs_{line}.csv"
    obs = pd.read_csv(path)
    log(f"[t1 obs] {line} columns actually read: {list(obs.columns)}")
    log_columns(f"louvain_obs_{line}", list(obs.columns), str(path))
    if "cluster" not in obs.columns or "day" not in obs.columns:
        raise StopStep("task1", f"louvain_obs_{line} missing cluster/day")
    lab_map = {int(r.cluster): str(r.label) for _, r in lab.iterrows()}
    obs = obs.reset_index(drop=True)
    obs["label"] = obs["cluster"].map(lambda c: lab_map.get(int(c)))
    if obs["label"].isna().any():
        raise StopStep(
            "task1",
            f"{line}: clusters without md2 labels: "
            f"{sorted(obs.loc[obs.label.isna(), 'cluster'].unique())}",
        )
    return obs


def _load_ams_md(line, obs, log):
    path = MD2_PROC / f"louvain_ams_{line}.npz"
    z = np.load(path, allow_pickle=True)
    log(f"[t1 ams] {line} keys actually read: {list(z.files)}")
    log_columns(f"louvain_ams_{line}", list(z.files), str(path))
    for k in ("md", "bins", "cluster", "day"):
        if k not in z.files:
            raise StopStep("task1", f"{path.name} missing field {k!r}. Not substituting.")
    md = np.asarray(z["md"], float)
    cl = np.asarray(z["cluster"], int)
    day = np.asarray(z["day"], int)
    bins = np.asarray(z["bins"])
    if len(md) != len(obs):
        raise StopStep("task1", f"{line} louvain_ams n_md={len(md)} obs n={len(obs)}")
    if not np.array_equal(cl, obs["cluster"].to_numpy(int)):
        raise StopStep("task1", f"{line} louvain_ams cluster does not match louvain_obs")
    if not np.array_equal(day, obs["day"].to_numpy(int)):
        raise StopStep("task1", f"{line} louvain_ams day does not match louvain_obs")
    gene_mean = np.asarray(z["gene_mean"]) if "gene_mean" in z.files else None
    return dict(md=md, bins=bins, gene_mean=gene_mean)


def load_md4_cell_scores(line, obs, log):
    path = MD4_PROC / f"cell_scores_{line}.npz"
    z = np.load(path, allow_pickle=True)
    log(f"[t1 cell_scores] {line} keys actually read: {list(z.files)}")
    log_columns(f"md4_cell_scores_{line}", list(z.files), str(path))
    need = {"frozen_ruler", "md_score", "age_up_minus_age_down", "pluri_primary",
            "cluster", "day", "label"}
    missing = sorted(need - set(z.files))
    if missing:
        raise StopStep("task1", f"{path.name} missing fields {missing}. Not substituting.")
    n = len(obs)
    cl = np.asarray(z["cluster"], int)
    day = np.asarray(z["day"], int)
    lab = np.asarray(z["label"]).astype(str)
    if len(cl) != n:
        raise StopStep("task1", f"{line} cell_scores n={len(cl)} obs n={n}")
    if not np.array_equal(cl, obs["cluster"].to_numpy(int)):
        raise StopStep("task1", f"{line} cell_scores cluster does not match louvain_obs")
    if not np.array_equal(day, obs["day"].to_numpy(int)):
        raise StopStep("task1", f"{line} cell_scores day does not match louvain_obs")
    if not np.array_equal(lab, obs["label"].to_numpy().astype(str)):
        raise StopStep("task1", f"{line} cell_scores label does not match md2 labels")
    out = {}
    for k in list(AGE_INSTRUMENTS) + [IDENTITY_KEY, "age_up", "age_down", "tgfb_score"]:
        if k not in z.files:
            continue
        arr = np.asarray(z[k], float)
        if len(arr) != n:
            raise StopStep("task1", f"{line} cell_scores {k} n={len(arr)} obs n={n}")
        out[k] = arr
    return out


def occupancy_and_qualify(log):
    count_rows = []
    comp_rows = []
    pooled = []
    for line in (AGED_LINE, YOUNG_LINE):
        lab = _load_labels(line, log)
        obs = _load_obs_labelled(line, lab, log)
        for st in STATE_NAMES:
            rec = dict(cell_line=line, label=st, n_clusters=int((lab.label.astype(str) == st).sum()))
            n_all = 0
            for d in DAYS:
                n = int(((obs.label == st) & (obs.day.astype(int) == int(d))).sum())
                rec[f"n_d{d}"] = n
                n_all += n
            rec["n_all"] = n_all
            rec["ge30"] = bool(n_all >= MIN_CELLS)
            rec["scored"] = bool(n_all >= MIN_CELLS)
            rec["reason"] = None if rec["scored"] else f"n<{MIN_CELLS}"
            pooled.append(rec)
            log(f"[t1 counts] {line} {st} n0={rec['n_d0']} n3={rec['n_d3']} "
                f"n7={rec['n_d7']} n10={rec['n_d10']} n_all={n_all} scored={rec['scored']}")
        for d in DAYS:
            n_tp = int((obs.day.astype(int) == int(d)).sum())
            for st in STATE_NAMES:
                n = int(((obs.label == st) & (obs.day.astype(int) == int(d))).sum())
                n_state = int((obs.label == st).sum())
                comp_rows.append(dict(
                    cell_line=line, day=int(d), label=st, n_cells=n,
                    n_timepoint=n_tp, n_state=n_state,
                    frac_state=float(n / n_state) if n_state else np.nan,
                    frac_timepoint=float(n / n_tp) if n_tp else np.nan,
                    state_ge30=bool(n_state >= MIN_CELLS),
                ))
                count_rows.append(dict(
                    cell_line=line, day=int(d), label=st, n_cells=n,
                    n_timepoint=n_tp, n_state=n_state,
                ))
    counts = pd.DataFrame(count_rows)
    pooled_df = pd.DataFrame(pooled)
    comp = pd.DataFrame(comp_rows)
    counts.to_csv(PLANE_DIR / "t1_cell_counts_by_day.csv", index=False)
    pooled_df.to_csv(PLANE_DIR / "t1_cell_counts.csv", index=False)
    comp.to_csv(PLANE_DIR / "t1_timepoint_composition.csv", index=False)

    md4 = pd.read_csv(MD4_DIR / "t1_cell_counts.csv")
    log_columns("md4_t1_cell_counts", list(md4.columns), str(MD4_DIR / "t1_cell_counts.csv"))
    for _, r in pooled_df.iterrows():
        hit = md4[(md4.cell_line == r.cell_line) & (md4.label == r.label)]
        if not len(hit):
            raise StopStep("task1", f"md4 occupancy missing {r.cell_line} {r.label}")
        md4_sum = int(hit.n_cells.sum())
        if int(r.n_all) != md4_sum:
            raise StopStep(
                "task1",
                f"occupancy drift {r.cell_line} {r.label}: plane n_all={int(r.n_all)} "
                f"md4 sum={md4_sum}. Labels have drifted.",
            )
        for d in DAYS:
            a = int(r[f"n_d{d}"])
            b = int(hit.loc[hit.day.astype(int) == int(d), "n_cells"].iloc[0])
            if a != b:
                raise StopStep(
                    "task1",
                    f"occupancy drift {r.cell_line} {r.label} d{d}: plane={a} md4={b}",
                )
    log("[t1] occupancy matches md4 t1_cell_counts.csv")

    qrows = []
    for line in (AGED_LINE, YOUNG_LINE):
        g = pooled_df[pooled_df.cell_line == line]
        scored = {str(r.label): bool(r.scored) for _, r in g.iterrows()}
        nmap = {str(r.label): int(r.n_all) for _, r in g.iterrows()}
        end = SEGMENT_B_PREFERRED_END if scored.get(SEGMENT_B_PREFERRED_END) else (
            SEGMENT_B_FALLBACK_END if scored.get(SEGMENT_B_FALLBACK_END) else None
        )
        end_is_fallback = bool(end == SEGMENT_B_FALLBACK_END and not scored.get(SEGMENT_B_PREFERRED_END))
        a0, a1 = SEGMENT_A
        ok_a = bool(scored.get(a0) and scored.get(a1))
        ok_b = bool(scored.get(a1) and end is not None)
        qrows.append(dict(
            cell_line=line,
            n_Fibroblast=nmap.get("Fibroblast", 0),
            n_PartialReprog=nmap.get("PartialReprog", 0),
            n_EarlyPluripotency=nmap.get("EarlyPluripotency", 0),
            n_Pluripotency=nmap.get("Pluripotency", 0),
            n_NonReprog=nmap.get("NonReprog", 0),
            segment_A_ok=ok_a,
            segment_B_ok=ok_b,
            segment_B_end=end if end is not None else "NA",
            segment_B_end_is_fallback=end_is_fallback,
            NonReprog_scored=bool(scored.get(OFF_TRAJECTORY)),
            min_cells=MIN_CELLS,
            reason=None if (ok_a and ok_b) else "a trajectory state n<30",
        ))
        log(f"[t1 qualify] {line} A_ok={ok_a} B_ok={ok_b} B_end={end} "
            f"fallback={end_is_fallback}")
    qdf = pd.DataFrame(qrows)
    qdf.to_csv(PLANE_DIR / "t1_qualify.csv", index=False)
    return counts, pooled_df, comp, qdf


def sanity_md2(log, obs_by, ams_by):
    repro = load_json(MD2_DIR / "t2_reproduction.json")
    log_columns("t2_reproduction", list(repro.keys()), str(MD2_DIR / "t2_reproduction.json"))
    obs = obs_by[AGED_LINE]
    md = ams_by[AGED_LINE]["md"]
    n_pr = int((obs.label == "PartialReprog").sum())
    n_nr = int((obs.label == "NonReprog").sum())
    md_pr = float(np.mean(md[obs.label.eq("PartialReprog").to_numpy()])) if n_pr else np.nan
    md_nr = float(np.mean(md[obs.label.eq("NonReprog").to_numpy()])) if n_nr else np.nan
    rec = dict(
        n_PartialReprog=n_pr, n_NonReprog=n_nr,
        md_mean_PartialReprog=md_pr, md_mean_NonReprog=md_nr,
        md2_json_n_PartialReprog=int(repro["n_PartialReprog"]),
        md2_json_n_NonReprog=int(repro["n_NonReprog"]),
        md2_json_md_mean_PartialReprog=float(repro["md_mean_PartialReprog"]),
        md2_json_md_mean_NonReprog=float(repro["md_mean_NonReprog"]),
        findings_printed_PartialReprog="+0.150",
        findings_printed_NonReprog="+0.439",
        n_match=bool(n_pr == int(repro["n_PartialReprog"]) and n_nr == int(repro["n_NonReprog"])),
        md_match_1e12=bool(
            abs(md_pr - float(repro["md_mean_PartialReprog"])) < 1e-12
            and abs(md_nr - float(repro["md_mean_NonReprog"])) < 1e-12
        ),
        printed_3dp_match=bool(f"{md_pr:+.3f}" == "+0.150" and f"{md_nr:+.3f}" == "+0.439"),
        expected_constants=dict(n_pr=MD2_N_PR, n_nr=MD2_N_NR, md_pr=MD2_POOLED_MD_PR, md_nr=MD2_POOLED_MD_NR),
    )
    rec["ok"] = bool(rec["n_match"] and rec["md_match_1e12"] and rec["printed_3dp_match"])
    dump_json(PLANE_DIR / "t1_sanity.json", jsonable(rec))
    log(f"[t1 sanity] n_PR={n_pr} n_NR={n_nr} MD_PR={md_pr:+.6f} MD_NR={md_nr:+.6f} ok={rec['ok']}")
    if not rec["ok"]:
        raise StopStep(
            "sanity",
            "FINDINGS_MD2.md pooled MD means did not reproduce. "
            "Cell-state labels have drifted and nothing downstream is valid. "
            f"got n_PR={n_pr} n_NR={n_nr} MD_PR={md_pr} MD_NR={md_nr}",
            details=rec,
        )
    return rec


def _finite_pair(x, y):
    x = np.asarray(x, float)
    y = np.asarray(y, float)
    m = np.isfinite(x) & np.isfinite(y)
    return x[m], y[m]


def _median_ci(vals, boot_meds):
    vals = np.asarray(vals, float)
    rec = dict(
        n=int(len(vals)),
        median=float(np.median(vals)) if len(vals) else np.nan,
    )
    boot_meds = np.asarray(boot_meds, float)
    boot_meds = boot_meds[np.isfinite(boot_meds)]
    if len(boot_meds):
        rec["ci_lo"] = float(np.percentile(boot_meds, 2.5))
        rec["ci_hi"] = float(np.percentile(boot_meds, 97.5))
        rec["n_boot_finite"] = int(len(boot_meds))
    else:
        rec["ci_lo"] = np.nan
        rec["ci_hi"] = np.nan
        rec["n_boot_finite"] = 0
    return rec


def bootstrap_states(obs, x, y, n_boot, seed, states):
    """Resample cells within each state. Returns dict state -> (x_meds, y_meds)."""
    rng = np.random.default_rng(int(seed))
    out = {}
    packed = {}
    for st in states:
        m = obs.label.eq(st).to_numpy()
        xx, yy = _finite_pair(x[m], y[m])
        packed[st] = (xx, yy)
        n = int(len(xx))
        xm = np.full(int(n_boot), np.nan)
        ym = np.full(int(n_boot), np.nan)
        if n >= 1:
            for b in range(int(n_boot)):
                idx = rng.choice(n, size=n, replace=True)
                xm[b] = float(np.median(xx[idx]))
                ym[b] = float(np.median(yy[idx]))
        out[st] = (xm, ym)
    return out, packed


def classify_one(srec):
    flags = dict(
        slope_A_young=flag_true(srec.get("slope_A_young")),
        slope_B_young=flag_true(srec.get("slope_B_young")),
        slope_B_flat=flag_true(srec.get("slope_B_flat")),
        slope_B_old=flag_true(srec.get("slope_B_old")),
        slopes_overlap=flag_true(srec.get("slopes_overlap")),
        bend_beats_null=flag_true(srec.get("bend_beats_null")),
        age_falls_A=flag_true(srec.get("age_falls_A")),
        age_falls_B=flag_true(srec.get("age_falls_B")),
        ok=flag_true(srec.get("ok")),
    )
    if not flags["ok"]:
        return "not_scored", flags, f"not scored: {srec.get('reason')}"
    young_A = flags["slope_A_young"]
    B_flat = flags["slope_B_flat"]
    B_old = flags["slope_B_old"]
    B_young = flags["slope_B_young"]
    overlap = flags["slopes_overlap"]
    beats = flags["bend_beats_null"]
    if young_A and (B_flat or B_old) and beats:
        return "bend", flags, (
            "age falls in segment A (slope_A CI excludes 0 in the young direction) "
            "and segment B is flat or reversed and the bend statistic beats its null. "
            "Age moves before identity is lost; rejuvenation is separable from "
            "dedifferentiation on this instrument."
        )
    if young_A and B_young and overlap and (not beats):
        return "lockstep", flags, (
            "age falls in both segments with overlapping slopes and the bend statistic "
            "does not beat its null. This instrument cannot separate rejuvenation from "
            "identity loss."
        )
    bits = []
    bits.append(f"slope_A_young={young_A}")
    bits.append(f"slope_B_young={B_young}")
    bits.append(f"slope_B_flat={B_flat}")
    bits.append(f"slope_B_old={B_old}")
    bits.append(f"slopes_overlap={overlap}")
    bits.append(f"bend_beats_null={beats}")
    bits.append(f"age_falls_A={srec.get('age_falls_A')}")
    bits.append(f"age_falls_B={srec.get('age_falls_B')}")
    bits.append(f"p_more_negative={srec.get('p_more_negative')}")
    return "mixed", flags, "mixed / other: " + "; ".join(bits)


def run_task1(log=None):
    close_log = False
    if log is None:
        log = Logger(PLANE_DIR / "t1_report.txt")
        close_log = True
    plane_log_banner(log, "TASK1")
    log(PREREG_TASK1)
    _require_prereg()
    sets = load_genesets(log)
    progress_snapshot("Task 1 occupancy from md2 labels (no new scores)", stop="GENESETS_DONE")

    obs_by, ams_by, lab_by = {}, {}, {}
    for line in (AGED_LINE, YOUNG_LINE):
        lab_by[line] = _load_labels(line, log)
        obs_by[line] = _load_obs_labelled(line, lab_by[line], log)
        ams_by[line] = _load_ams_md(line, obs_by[line], log)
        Xl, symbols, gene_id, _ = csr_from_npz(MD2_PROC / f"louvain_counts_{line}.npz")
        id_sym = detect_id_type(symbols, log, f"louvain_{line}_symbols",
                                str(MD2_PROC / f"louvain_counts_{line}.npz") + ":symbols")
        require_mappable(id_sym, {"symbol"}, log, "task1")
        if gene_id is not None and len(gene_id):
            detect_id_type(gene_id, log, f"louvain_{line}_gene_id",
                           str(MD2_PROC / f"louvain_counts_{line}.npz") + ":gene_id")
        if Xl.shape[0] != len(obs_by[line]):
            raise StopStep("task1", f"{line} louvain counts n={Xl.shape[0]} obs n={len(obs_by[line])}")

    counts, pooled, comp, qdf = occupancy_and_qualify(log)
    sanity = sanity_md2(log, obs_by, ams_by)
    if not DECLARED_BEFORE_SCORES_FLAG.exists():
        DECLARED_BEFORE_SCORES_FLAG.write_text(
            PREREG_TASK1
            + "\n\nThis flag was written after cell-count tables, qualify tables, and the "
              "md2 MD-mean sanity check, and before any state median, slope, bend statistic, "
              "or null p-value.\n"
            + f"\nsanity_ok={sanity['ok']}\n",
            encoding="utf-8",
        )
        log(f"[prereg] wrote {DECLARED_BEFORE_SCORES_FLAG}")
    else:
        log(f"[prereg] {DECLARED_BEFORE_SCORES_FLAG.name} already exists; not rewriting")
    progress_snapshot("state medians, slopes, then nulls (observed scores from md4/md2 disk)",
                      stop="DECLARED_BEFORE_SCORES")

    scores_by = {}
    for line in (AGED_LINE, YOUNG_LINE):
        scores_by[line] = load_md4_cell_scores(line, obs_by[line], log)
        md_ams = ams_by[line]["md"]
        md_sc = scores_by[line]["md_score"]
        if not np.allclose(md_ams, md_sc, equal_nan=True, atol=0, rtol=0):
            maxd = float(np.nanmax(np.abs(md_ams - md_sc)))
            raise StopStep(
                "task1",
                f"{line}: md4 cell_scores md_score does not equal md2 louvain_ams md. "
                f"max_abs_delta={maxd}. Not substituting.",
            )
        log(f"[t1] {line} md_score matches louvain_ams exactly")

    median_rows = []
    boot_pack = {}
    for line in (AGED_LINE, YOUNG_LINE):
        obs = obs_by[line]
        sc = scores_by[line]
        x = np.asarray(sc[IDENTITY_KEY], float)
        qrow = qdf[qdf.cell_line == line].iloc[0].to_dict()
        boot_pack[line] = {}
        for inst in AGE_INSTRUMENTS:
            y = np.asarray(sc[inst], float)
            boots, packed = bootstrap_states(
                obs, x, y, N_BOOT, boot_seed(line, inst), STATE_NAMES,
            )
            boot_pack[line][inst] = dict(boots=boots, packed=packed)
            for st in STATE_NAMES:
                xx, yy = packed[st]
                n_state = int((obs.label == st).sum())
                scored = bool(n_state >= MIN_CELLS)
                rec = dict(
                    cell_line=line, instrument=inst, label=st,
                    n_cells=n_state, n_finite=int(len(xx)),
                    scored=scored, min_cells=MIN_CELLS,
                    on_trajectory=bool(st in TRAJECTORY),
                    reason=None if scored else f"n<{MIN_CELLS}",
                )
                if not scored:
                    rec.update(
                        identity_median=np.nan, identity_ci_lo=np.nan, identity_ci_hi=np.nan,
                        age_median=np.nan, age_ci_lo=np.nan, age_ci_hi=np.nan,
                    )
                    median_rows.append(rec)
                    continue
                xm, ym = boots[st]
                ix = _median_ci(xx, xm)
                iy = _median_ci(yy, ym)
                rec.update(
                    identity_median=ix["median"], identity_ci_lo=ix["ci_lo"],
                    identity_ci_hi=ix["ci_hi"],
                    age_median=iy["median"], age_ci_lo=iy["ci_lo"],
                    age_ci_hi=iy["ci_hi"],
                    n_boot=int(N_BOOT), boot_seed=int(boot_seed(line, inst)),
                )
                median_rows.append(rec)
                log(f"[t1 median] {line} {inst} {st} n={n_state} "
                    f"x={ix['median']:+.4f} [{ix['ci_lo']:+.4f},{ix['ci_hi']:+.4f}] "
                    f"y={iy['median']:+.4f} [{iy['ci_lo']:+.4f},{iy['ci_hi']:+.4f}]")
    mdf = pd.DataFrame(median_rows)
    mdf.to_csv(PLANE_DIR / "t1_state_medians.csv", index=False)

    slope_rows = []
    for line in (AGED_LINE, YOUNG_LINE):
        qrow = qdf[qdf.cell_line == line].iloc[0].to_dict()
        end, fallback = segment_B_end(qrow)
        for inst in AGE_INSTRUMENTS:
            packed = boot_pack[line][inst]["packed"]
            boots = boot_pack[line][inst]["boots"]
            rec = dict(
                cell_line=line, instrument=inst,
                state_A0=SEGMENT_A[0], state_A1=SEGMENT_A[1],
                state_B0=SEGMENT_A[1], state_B1=end if end else "NA",
                segment_B_end_is_fallback=fallback,
                n_boot=int(N_BOOT), boot_seed=int(boot_seed(line, inst)),
                min_cells=MIN_CELLS,
            )
            a0, a1 = SEGMENT_A
            ok_a = bool(qrow.get("segment_A_ok"))
            ok_b = bool(qrow.get("segment_B_ok")) and end is not None
            if not (ok_a and ok_b):
                rec.update(ok=False, reason="trajectory state n<30")
                slope_rows.append(rec)
                continue
            med = {}
            for st in (a0, a1, end):
                xx, yy = packed[st]
                if len(xx) < MIN_CELLS:
                    rec.update(ok=False, reason=f"{st} n_finite<{MIN_CELLS}")
                    slope_rows.append(rec)
                    ok_a = False
                    break
                med[st] = dict(x=float(np.median(xx)), y=float(np.median(yy)), n=int(len(xx)))
            if not ok_a:
                continue
            dx_A = med[a1]["x"] - med[a0]["x"]
            dy_A = med[a1]["y"] - med[a0]["y"]
            dx_B = med[end]["x"] - med[a1]["x"]
            dy_B = med[end]["y"] - med[a1]["y"]
            sA = _slope(dy_A, dx_A)
            sB = _slope(dy_B, dx_B)
            bend = (sA - sB) if (np.isfinite(sA) and np.isfinite(sB)) else np.nan
            rec.update(
                ok=bool(np.isfinite(sA) and np.isfinite(sB)),
                n_A0=med[a0]["n"], n_A1=med[a1]["n"], n_B1=med[end]["n"],
                x_A0=med[a0]["x"], x_A1=med[a1]["x"], x_B1=med[end]["x"],
                y_A0=med[a0]["y"], y_A1=med[a1]["y"], y_B1=med[end]["y"],
                delta_identity_A=dx_A, delta_age_A=dy_A, slope_A=sA,
                delta_identity_B=dx_B, delta_age_B=dy_B, slope_B=sB,
                bend=bend,
                reason=None if np.isfinite(bend) else "undefined slope (|Δidentity|<eps)",
            )
            xmA0, ymA0 = boots[a0]
            xmA1, ymA1 = boots[a1]
            xmB1, ymB1 = boots[end]
            boot_sA = np.full(N_BOOT, np.nan)
            boot_sB = np.full(N_BOOT, np.nan)
            boot_bend = np.full(N_BOOT, np.nan)
            boot_dyA = ymA1 - ymA0
            boot_dxA = xmA1 - xmA0
            boot_dyB = ymB1 - ymA1
            boot_dxB = xmB1 - xmA1
            for b in range(N_BOOT):
                boot_sA[b] = _slope(boot_dyA[b], boot_dxA[b])
                boot_sB[b] = _slope(boot_dyB[b], boot_dxB[b])
                if np.isfinite(boot_sA[b]) and np.isfinite(boot_sB[b]):
                    boot_bend[b] = boot_sA[b] - boot_sB[b]
            def _ci(v):
                v = np.asarray(v, float)
                v = v[np.isfinite(v)]
                if len(v) < 100:
                    return np.nan, np.nan, int(len(v))
                return float(np.percentile(v, 2.5)), float(np.percentile(v, 97.5)), int(len(v))
            sa_lo, sa_hi, sa_n = _ci(boot_sA)
            sb_lo, sb_hi, sb_n = _ci(boot_sB)
            be_lo, be_hi, be_n = _ci(boot_bend)
            dya_lo, dya_hi, _ = _ci(boot_dyA)
            dyb_lo, dyb_hi, _ = _ci(boot_dyB)
            rec.update(
                slope_A_ci_lo=sa_lo, slope_A_ci_hi=sa_hi, n_boot_slope_A=sa_n,
                slope_B_ci_lo=sb_lo, slope_B_ci_hi=sb_hi, n_boot_slope_B=sb_n,
                bend_ci_lo=be_lo, bend_ci_hi=be_hi, n_boot_bend=be_n,
                delta_age_A_ci_lo=dya_lo, delta_age_A_ci_hi=dya_hi,
                delta_age_B_ci_lo=dyb_lo, delta_age_B_ci_hi=dyb_hi,
            )
            side_A = slope_young_side(dx_A)
            side_B = slope_young_side(dx_B)
            rec["slope_A_young"] = ci_excludes_0_on_side(sa_lo, sa_hi, side_A)
            rec["slope_B_young"] = ci_excludes_0_on_side(sb_lo, sb_hi, side_B)
            rec["slope_B_flat"] = ci_includes_0(sb_lo, sb_hi)
            rec["slope_B_old"] = ci_excludes_0_on_side(sb_lo, sb_hi, None if side_B is None else -side_B)
            rec["slopes_overlap"] = cis_overlap(sa_lo, sa_hi, sb_lo, sb_hi)
            rec["age_falls_A"] = bool(np.isfinite(dya_hi) and dya_hi < 0.0)
            rec["age_falls_B"] = bool(np.isfinite(dyb_hi) and dyb_hi < 0.0)
            rec["age_rises_B"] = bool(np.isfinite(dyb_lo) and dyb_lo > 0.0)
            if not rec["ok"]:
                rec["reason"] = rec.get("reason") or "undefined observed slope"
            slope_rows.append(rec)
            log(f"[t1 slope] {line} {inst} A={sA} B={sB} bend={bend} "
                f"B_end={end} fallback={fallback}")
    sdf = pd.DataFrame(slope_rows)
    sdf.to_csv(PLANE_DIR / "t1_slopes.csv", index=False)

    progress_snapshot("Task 1 nulls (random directions / random gene sets)",
                      stop="MEDIANS_DONE")
    from plane_nulls import run_nulls  # noqa: E402
    null_sum, null_state, null_bend = run_nulls(
        log, obs_by, ams_by, scores_by, sets, qdf, sdf,
    )

    bend_rows = []
    reading_by = {}
    for line in (AGED_LINE, YOUNG_LINE):
        reading_by[line] = {}
        for inst in AGE_INSTRUMENTS:
            srow = sdf[(sdf.cell_line == line) & (sdf.instrument == inst)]
            if not len(srow):
                raise StopStep("task1", f"missing slope row {line} {inst}")
            srec = srow.iloc[0].to_dict()
            nrow = null_sum[(null_sum.cell_line == line) & (null_sum.instrument == inst)]
            rec = dict(srec)
            if len(nrow):
                nrec = nrow.iloc[0].to_dict()
                rec.update(
                    n_random=nrec.get("n_random"),
                    n_null_finite=nrec.get("n_null_finite"),
                    null_kind=nrec.get("null_kind"),
                    null_bend_median=nrec.get("null_bend_median"),
                    p_more_negative=nrec.get("p_more_negative"),
                    p_more_positive=nrec.get("p_more_positive"),
                    null_seed=nrec.get("null_seed"),
                )
            else:
                rec.update(
                    n_random=N_RANDOM_DIR, n_null_finite=0, null_kind="missing",
                    p_more_negative=np.nan, p_more_positive=np.nan,
                )
            p_neg = rec.get("p_more_negative")
            rec["bend_beats_null"] = bool(
                rec.get("ok") and np.isfinite(p_neg) and float(p_neg) <= 0.05
            )
            key, flags, text = classify_one(rec)
            rec["key"] = key
            rec["text"] = text
            rec.update({f"flag_{k}": v for k, v in flags.items()})
            bend_rows.append(rec)
            reading_by[line][inst] = dict(
                key=key, text=text, flags=flags,
                slope_A=rec.get("slope_A"), slope_B=rec.get("slope_B"),
                bend=rec.get("bend"),
                p_more_negative=rec.get("p_more_negative"),
                p_more_positive=rec.get("p_more_positive"),
                segment_B_end=rec.get("state_B1"),
                segment_B_end_is_fallback=rec.get("segment_B_end_is_fallback"),
                null_kind=rec.get("null_kind"),
            )
            log(f"[t1 reading] {line} {inst} key={key} p_neg={rec.get('p_more_negative')} "
                f"bend={rec.get('bend')}")
    bdf = pd.DataFrame(bend_rows)
    bdf.to_csv(PLANE_DIR / "t1_bend.csv", index=False)

    inst_cmp = []
    for inst in AGE_INSTRUMENTS:
        keys = {line: (reading_by[line][inst]["key"]) for line in (AGED_LINE, YOUNG_LINE)}
        inst_cmp.append(dict(instrument=inst, **{f"key_{k}": v for k, v in keys.items()}))
    comparison = dict(
        by_instrument=inst_cmp,
        md_expected_lockstep_from_md3=True,
        md3_rho={
            AGED_LINE: MD3_MD_VS_PLURI_GM00731,
            YOUNG_LINE: MD3_MD_VS_PLURI_GM23815,
        },
        md3_n={AGED_LINE: 33, YOUNG_LINE: 34},
        note=(
            "FINDINGS_MD3.md Task 3: MD vs pluri_primary "
            f"ρ={MD3_MD_VS_PLURI_GM00731} (GM00731 n=33) and "
            f"ρ={MD3_MD_VS_PLURI_GM23815} (GM23815 n=34). "
            "MD is expected to run in lockstep. No winner is declared."
        ),
    )

    reading = dict(
        by_donor=reading_by,
        instruments=list(AGE_INSTRUMENTS),
        trajectory=list(TRAJECTORY),
        min_cells=MIN_CELLS,
        n_boot=N_BOOT,
        n_random=N_RANDOM_DIR,
        seed=PLANE_SEED,
        boot_seed=PLANE_BOOT,
        comparison=comparison,
        identity=IDENTITY_KEY,
        NonReprog_excluded_from_trajectory_statistics=True,
    )
    dump_json(PLANE_DIR / "t1_reading.json", jsonable(reading))
    dump_json(PLANE_DIR / "t1_summary.json", jsonable(dict(
        n_median_rows=int(len(mdf)),
        n_slope_rows=int(len(sdf)),
        n_bend_rows=int(len(bdf)),
        sanity_ok=bool(sanity["ok"]),
        by_donor={ln: {inst: reading_by[ln][inst]["key"] for inst in AGE_INSTRUMENTS}
                  for ln in (AGED_LINE, YOUNG_LINE)},
    )))
    man = load_manifest()
    man["status"] = "TASK1_DONE"
    save_manifest(man)
    progress_snapshot("Task 2 extrapolation (report-only)", stop="TASK1_DONE")
    if close_log:
        log.close()
    return reading


if __name__ == "__main__":
    try:
        run_task1()
    except StopStep as e:
        from plane_common import record_failure  # noqa: E402
        record_failure(e.step, e.message, e.details)
        progress_snapshot(f"fix STOP [{e.step}]", stop=f"STOP {e.step}: {e.message}")
        raise
