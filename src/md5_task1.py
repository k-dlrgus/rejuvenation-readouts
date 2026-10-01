"""Task 1 — PartialReprog vs Fibroblast under scheme (a) and scheme (b)."""
from __future__ import annotations

import sys
from pathlib import Path

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parent))
from md5_common import (  # noqa: E402
    MD5_DIR, MD2_DIR, MD2_PROC, MD3_DIR, MD4_DIR, MD4_PROC,
    MD5_SEED, AGED_LINE, YOUNG_LINE, DAYS, STATE_NAMES, INSTRUMENTS,
    PRIMARY_PAIR, MIN_CELLS_A, SCHEME_B_FIB_DAYS, SCHEME_B_PR_DAYS,
    PREREG_TASK1, PREREG_TASK1_FLAG, DECLARED_BEFORE_SCORES_FLAG,
    FROZEN_RULER, N_PERM, N_BOOT, PAIR_INDEX,
    MD2_POOLED_MD_PR, MD2_POOLED_MD_NR, MD2_N_PR, MD2_N_NR,
    MD4_TASK1_YOUNGER_CLAIM, WHAT_WOULD_SETTLE,
    md4_perm_seed, md4_boot_seed, perm_seed_b, boot_seed_b,
    StopStep, Logger, dump_json, jsonable, load_json, md5_log_banner,
    load_manifest, save_manifest, record_failure, log_columns, log_geneset,
    progress_snapshot, detect_id_type, require_mappable,
)
from md4_task1 import _one_contrast, cohens_d  # noqa: E402
from trajectory_common import permutation_p  # noqa: E402


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
        ):
            if not p.exists():
                raise StopStep("task1", f"missing artifact {p}. Not substituting.")
    for p in (
        MD2_DIR / "t2_reproduction.json",
        MD2_DIR / "genesets.json",
        MD3_DIR / "genesets.json",
        MD4_DIR / "t1_cell_counts.csv",
        MD4_DIR / "t1_contrasts.csv",
        MD4_DIR / "genesets.json",
        MD4_DIR / "genesets_summary.json",
    ):
        if not p.exists():
            raise StopStep("task1", f"missing artifact {p}. Not substituting.")
    if PRIMARY_PAIR not in PAIR_INDEX:
        raise StopStep("task1", f"{PRIMARY_PAIR} missing from md4 PAIR_INDEX. Not substituting.")


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
    for k in ("md", "cluster", "day"):
        if k not in z.files:
            raise StopStep("task1", f"{path.name} missing field {k!r}. Not substituting.")
    md = np.asarray(z["md"], float)
    cl = np.asarray(z["cluster"], int)
    day = np.asarray(z["day"], int)
    if len(md) != len(obs):
        raise StopStep("task1", f"{line} louvain_ams n_md={len(md)} obs n={len(obs)}")
    if not np.array_equal(cl, obs["cluster"].to_numpy(int)) or not np.array_equal(day, obs["day"].to_numpy(int)):
        raise StopStep(
            "task1",
            f"{line} louvain_ams cluster/day arrays do not match louvain_obs. "
            "Labels have drifted. Not scoring.",
        )
    return md


def inherit_genesets(log):
    md2_gs = load_json(MD2_DIR / "genesets.json")
    md3_gs = load_json(MD3_DIR / "genesets.json")
    md4_gs = load_json(MD4_DIR / "genesets.json")
    md4_sum = load_json(MD4_DIR / "genesets_summary.json")
    log_columns("md2_genesets_json", list(md2_gs.keys()), str(MD2_DIR / "genesets.json"))
    log_columns("md3_genesets_json", list(md3_gs.keys()), str(MD3_DIR / "genesets.json"))
    log_columns("md4_genesets_json", list(md4_gs.keys()), str(MD4_DIR / "genesets.json"))
    md = [str(g).upper() for g in md2_gs["MD"]]
    tgfb = [str(g).upper() for g in md2_gs["TGFB"]]
    age_up = [str(g).upper() for g in md3_gs["age_up"]]
    age_down = [str(g).upper() for g in md3_gs["age_down"]]
    md4_md = [str(g).upper() for g in md4_gs["MD"]]
    md4_tgfb = [str(g).upper() for g in md4_gs["TGFB"]]
    md4_up = [str(g).upper() for g in md4_gs["age_up"]]
    md4_dn = [str(g).upper() for g in md4_gs["age_down"]]
    if set(md) != set(md4_md) or set(tgfb) != set(md4_tgfb) or set(age_up) != set(md4_up) or set(age_down) != set(md4_dn):
        raise StopStep("genesets", "md2/md3 gene lists do not equal md4/genesets.json. Not reconciling.")
    log_geneset("MD_built_from_md2", "MD (STAR Methods built; md2)", md,
                str(MD2_DIR / "genesets.json"), id_type="symbol")
    log_geneset("TGFB_built_from_md2", "TGFB (STAR Methods built; md2)", tgfb,
                str(MD2_DIR / "genesets.json"), id_type="symbol")
    log_geneset("mmc3_Aging_signatures_Age up", "Aging_signatures:Age up", age_up,
                str(MD4_DIR / "genesets.json"), id_type="symbol")
    log_geneset("mmc3_Aging_signatures_Age down", "Aging_signatures:Age down", age_down,
                str(MD4_DIR / "genesets.json"), id_type="symbol")
    s3_meta = {}
    for line in (AGED_LINE, YOUNG_LINE):
        p = MD4_DIR / f"t1_s3_ams_meta_{line}.json"
        if not p.exists():
            raise StopStep("genesets", f"missing {p}. Not reconstructing mapped counts from memory.")
        rec = load_json(p)
        log_columns(f"t1_s3_ams_meta_{line}", list(rec.keys()), str(p))
        s3_meta[line] = {
            k: {kk: rec[k][kk] for kk in rec[k] if kk != "missing"}
            for k in rec
        }
    summary = dict(md4_genesets_summary=md4_sum, s3_mapped=s3_meta, inherited_from_md4=True)
    dump_json(MD5_DIR / "genesets_summary.json", jsonable(summary))
    dump_json(MD5_DIR / "genesets.json", jsonable(dict(
        MD=md, TGFB=tgfb, age_up=age_up, age_down=age_down,
    )))
    log(f"[genesets] inherited md4 MD n={len(md)} TGFB n={len(tgfb)} "
        f"Age up n={len(age_up)} Age down n={len(age_down)}")
    return dict(MD=md, TGFB=tgfb, age_up=age_up, age_down=age_down, s3_meta=s3_meta)


def log_id_types_louvain(log):
    for line in (AGED_LINE, YOUNG_LINE):
        path = MD2_PROC / f"louvain_counts_{line}.npz"
        z = np.load(path, allow_pickle=True)
        log(f"[id_type] {line} louvain_counts keys={list(z.files)}")
        if "symbols" not in z.files:
            raise StopStep("task1", f"{path} missing field 'symbols'. Not substituting.")
        id_sym = detect_id_type(
            z["symbols"], log, f"louvain_{line}_symbols", str(path) + ":symbols",
        )
        require_mappable(id_sym, {"symbol"}, log, "task1")
        if "gene_id" in z.files:
            detect_id_type(
                z["gene_id"], log, f"louvain_{line}_gene_id", str(path) + ":gene_id",
            )


def occupancy(log, obs_by):
    rows = []
    pooled = []
    fib_pr = []
    for line in (AGED_LINE, YOUNG_LINE):
        obs = obs_by[line]
        lab = pd.read_csv(MD2_DIR / f"t2_cluster_labels_{line}.csv")
        for st in STATE_NAMES:
            rec = dict(cell_line=line, label=st, n_clusters=int((lab.label.astype(str) == st).sum()))
            for d in DAYS:
                rec[f"n_d{d}"] = int(((obs.label == st) & (obs.day.astype(int) == int(d))).sum())
            rec["n_all"] = int((obs.label == st).sum())
            pooled.append(rec)
        for d in DAYS:
            n_tp = int((obs.day.astype(int) == int(d)).sum())
            for st in STATE_NAMES:
                n = int(((obs.label == st) & (obs.day.astype(int) == int(d))).sum())
                rows.append(dict(
                    cell_line=line, day=int(d), label=st, n_cells=n,
                    n_timepoint=n_tp,
                    frac_timepoint=float(n / n_tp) if n_tp else np.nan,
                    ge30=bool(n >= MIN_CELLS_A),
                    ge50=bool(n >= 50),
                ))
            n_fib = int(((obs.label == "Fibroblast") & (obs.day.astype(int) == int(d))).sum())
            n_pr = int(((obs.label == "PartialReprog") & (obs.day.astype(int) == int(d))).sum())
            fib_pr.append(dict(
                cell_line=line, day=int(d), n_Fibroblast=n_fib, n_PartialReprog=n_pr,
                n_timepoint=n_tp, both_ge30=bool(n_fib >= MIN_CELLS_A and n_pr >= MIN_CELLS_A),
            ))
            log(f"[t1 counts] {line} d{d} Fibroblast n={n_fib} PartialReprog n={n_pr} "
                f"both_ge30={n_fib >= MIN_CELLS_A and n_pr >= MIN_CELLS_A}")
    counts = pd.DataFrame(rows)
    pooled_df = pd.DataFrame(pooled)
    fib_pr_df = pd.DataFrame(fib_pr)
    counts.to_csv(MD5_DIR / "t1_cell_counts.csv", index=False)
    pooled_df.to_csv(MD5_DIR / "t1_cell_counts_pooled.csv", index=False)
    fib_pr_df.to_csv(MD5_DIR / "t1_fib_pr_counts.csv", index=False)

    md4 = pd.read_csv(MD4_DIR / "t1_cell_counts.csv")
    log_columns("md4_t1_cell_counts", list(md4.columns), str(MD4_DIR / "t1_cell_counts.csv"))
    need = {"cell_line", "day", "label", "n_cells"}
    if not need <= set(md4.columns):
        raise StopStep("task1", f"md4 t1_cell_counts.csv missing {sorted(need - set(md4.columns))}")
    for _, r in counts.iterrows():
        hit = md4[
            (md4.cell_line == r.cell_line) & (md4.day.astype(int) == int(r.day)) & (md4.label == r.label)
        ]
        if not len(hit):
            raise StopStep("task1", f"md4 occupancy missing {r.cell_line} d{r.day} {r.label}")
        a, b = int(r.n_cells), int(hit.iloc[0].n_cells)
        if a != b:
            raise StopStep(
                "task1",
                f"occupancy drift {r.cell_line} {r.label} d{r.day}: md5={a} md4={b}",
            )
    log("[t1] occupancy matches md4 t1_cell_counts.csv")

    q_a = []
    for line in (AGED_LINE, YOUNG_LINE):
        for d in DAYS:
            n_pr = int(fib_pr_df.loc[
                (fib_pr_df.cell_line == line) & (fib_pr_df.day == int(d)), "n_PartialReprog",
            ].iloc[0])
            n_fib = int(fib_pr_df.loc[
                (fib_pr_df.cell_line == line) & (fib_pr_df.day == int(d)), "n_Fibroblast",
            ].iloc[0])
            ok = bool(n_pr >= MIN_CELLS_A and n_fib >= MIN_CELLS_A)
            q_a.append(dict(
                scheme="a", cell_line=line, day=int(d),
                state_a="PartialReprog", state_b="Fibroblast",
                n_a=n_pr, n_b=n_fib, min_cells=MIN_CELLS_A,
                qualifies=ok, primary=True,
                reason=None if ok else f"n<{MIN_CELLS_A} in at least one state",
            ))
            log(f"[t1 qualify a] {line} d{d} PR n={n_pr} Fib n={n_fib} qualifies={ok}")
    q_a_df = pd.DataFrame(q_a)
    q_a_df.to_csv(MD5_DIR / "t1_qualify_a.csv", index=False)

    comp_b = []
    q_b = []
    for line in (AGED_LINE, YOUNG_LINE):
        obs = obs_by[line]
        n_fib_days = {}
        n_pr_days = {}
        for d in DAYS:
            n_fib_days[int(d)] = int(((obs.label == "Fibroblast") & (obs.day.astype(int) == int(d))).sum())
            n_pr_days[int(d)] = int(((obs.label == "PartialReprog") & (obs.day.astype(int) == int(d))).sum())
        n_fib = int(sum(n_fib_days[d] for d in SCHEME_B_FIB_DAYS))
        n_pr = int(sum(n_pr_days[d] for d in SCHEME_B_PR_DAYS))
        for d in SCHEME_B_FIB_DAYS:
            n = n_fib_days[int(d)]
            comp_b.append(dict(
                cell_line=line, arm="Fibroblast", day=int(d), n_cells=n,
                n_arm=n_fib, frac_arm=float(n / n_fib) if n_fib else np.nan,
            ))
        for d in SCHEME_B_PR_DAYS:
            n = n_pr_days[int(d)]
            comp_b.append(dict(
                cell_line=line, arm="PartialReprog", day=int(d), n_cells=n,
                n_arm=n_pr, frac_arm=float(n / n_pr) if n_pr else np.nan,
            ))
        ok = bool(n_fib >= MIN_CELLS_A and n_pr >= MIN_CELLS_A)
        a_qualifies_any = bool(q_a_df[(q_a_df.cell_line == line) & (q_a_df.qualifies)].shape[0] > 0)
        q_b.append(dict(
            scheme="b", cell_line=line,
            state_a="PartialReprog", state_b="Fibroblast",
            n_a=n_pr, n_b=n_fib, min_cells=MIN_CELLS_A,
            qualifies=ok,
            primary=bool(ok and not a_qualifies_any),
            weaker_design=True,
            scheme_a_qualifies_on_this_donor=a_qualifies_any,
            fib_days=list(SCHEME_B_FIB_DAYS), pr_days=list(SCHEME_B_PR_DAYS),
            n_Fibroblast_d0=n_fib_days[0], n_Fibroblast_d3=n_fib_days[3],
            n_PartialReprog_d3=n_pr_days[3], n_PartialReprog_d7=n_pr_days[7],
            reason=None if ok else f"n<{MIN_CELLS_A} in at least one pooled arm",
            note=("weaker design; permutation shuffles state label within timepoint strata; "
                  "not primary where scheme (a) qualified"),
        ))
        log(f"[t1 qualify b] {line} Fib d0+d3 n={n_fib} PR d3+d7 n={n_pr} qualifies={ok} "
            f"primary={ok and not a_qualifies_any}")
    comp_b_df = pd.DataFrame(comp_b)
    q_b_df = pd.DataFrame(q_b)
    comp_b_df.to_csv(MD5_DIR / "t1_scheme_b_composition.csv", index=False)
    q_b_df.to_csv(MD5_DIR / "t1_qualify_b.csv", index=False)
    return counts, pooled_df, fib_pr_df, q_a_df, q_b_df, comp_b_df


def sanity_md2(log, obs_by, md_by):
    repro = load_json(MD2_DIR / "t2_reproduction.json")
    log_columns("t2_reproduction", list(repro.keys()), str(MD2_DIR / "t2_reproduction.json"))
    obs = obs_by[AGED_LINE]
    md = md_by[AGED_LINE]
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
    dump_json(MD5_DIR / "t1_sanity.json", jsonable(rec))
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


def load_md4_cell_scores(line, obs, log):
    path = MD4_PROC / f"cell_scores_{line}.npz"
    z = np.load(path, allow_pickle=True)
    log(f"[t1 cell_scores] {line} keys actually read: {list(z.files)}")
    log_columns(f"md4_cell_scores_{line}", list(z.files), str(path))
    need = set(INSTRUMENTS) | {"cluster", "day", "label"}
    missing = sorted(need - set(z.files))
    if missing:
        raise StopStep("task1", f"{path.name} missing fields {missing}. Not rescoring.")
    n = len(obs)
    cl = np.asarray(z["cluster"], int)
    day = np.asarray(z["day"], int)
    lab = np.asarray(z["label"]).astype(str)
    if len(cl) != n or len(day) != n or len(lab) != n:
        raise StopStep("task1", f"{line} cell_scores n={len(cl)} obs n={n}")
    if not np.array_equal(cl, obs["cluster"].to_numpy(int)):
        raise StopStep("task1", f"{line} cell_scores cluster does not match louvain_obs")
    if not np.array_equal(day, obs["day"].to_numpy(int)):
        raise StopStep("task1", f"{line} cell_scores day does not match louvain_obs")
    if not np.array_equal(lab, obs["label"].astype(str).to_numpy()):
        raise StopStep("task1", f"{line} cell_scores label does not match md2 labels")
    out = {inst: np.asarray(z[inst], float) for inst in INSTRUMENTS}
    for inst, arr in out.items():
        if len(arr) != n:
            raise StopStep("task1", f"{line} cell_scores {inst} n={len(arr)} obs n={n}")
    return out


def _one_contrast_stratified(x, mask_a, strata, n_perm, n_boot, seed_p, seed_b):
    """Same statistics as md4_task1._one_contrast; permutation shuffles within strata."""
    x = np.asarray(x, float)
    mask_a = np.asarray(mask_a, bool)
    strata = np.asarray(strata)
    ok = np.isfinite(x)
    x = x[ok]
    mask_a = mask_a[ok]
    strata = strata[ok]
    a = x[mask_a]
    b = x[~mask_a]
    n_a, n_b = int(len(a)), int(len(b))
    rec = dict(
        n_a=n_a, n_b=n_b, n_perm=int(n_perm), n_boot=int(n_boot),
        perm_seed=int(seed_p), boot_seed=int(seed_b),
        effect_size="cohens_d_pooled_nminus1",
        permutation="state_label_shuffled_within_timepoint_strata",
    )
    if n_a < 2 or n_b < 2:
        rec.update(ok=False, reason="n<2 after dropping non-finite scores")
        return rec
    dmean = float(a.mean() - b.mean())
    dmed = float(np.median(a) - np.median(b))
    rec.update(
        ok=True,
        mean_a=float(a.mean()), mean_b=float(b.mean()), delta_mean=dmean,
        median_a=float(np.median(a)), median_b=float(np.median(b)), delta_median=dmed,
        cohens_d=cohens_d(a, b),
    )
    rng = np.random.default_rng(int(seed_p))
    null_mean = np.empty(int(n_perm), dtype=np.float64)
    null_med = np.empty(int(n_perm), dtype=np.float64)
    lab0 = mask_a.copy()
    strata_ids = np.unique(strata)
    rec["n_strata"] = int(len(strata_ids))
    rec["n_strata_mixed"] = int(sum(
        (lab0[strata == s].any() and (~lab0[strata == s]).any()) for s in strata_ids
    ))
    for i in range(int(n_perm)):
        lab = lab0.copy()
        for s in strata_ids:
            idx = np.flatnonzero(strata == s)
            if len(idx) > 1:
                lab[idx] = rng.permutation(lab[idx])
        aa, bb = x[lab], x[~lab]
        null_mean[i] = float(aa.mean() - bb.mean())
        null_med[i] = float(np.median(aa) - np.median(bb))
    rec["p_mean_lower"] = permutation_p(dmean, null_mean, greater=False)
    rec["p_mean_two_sided"] = float(
        (np.sum(np.abs(null_mean) >= abs(dmean)) + 1) / (len(null_mean) + 1)
    )
    rec["n_perm_mean_le_obs"] = int(np.sum(null_mean <= dmean))
    rec["p_median_lower"] = permutation_p(dmed, null_med, greater=False)
    rec["p_median_two_sided"] = float(
        (np.sum(np.abs(null_med) >= abs(dmed)) + 1) / (len(null_med) + 1)
    )
    rec["null_mean_median"] = float(np.median(null_mean))
    rec["beats_null_lower"] = bool(dmean < 0 and rec["p_mean_lower"] <= 0.05)

    rng_b = np.random.default_rng(int(seed_b))
    boot_m = np.empty(int(n_boot), dtype=np.float64)
    boot_d = np.empty(int(n_boot), dtype=np.float64)
    for i in range(int(n_boot)):
        aa = rng_b.choice(a, size=n_a, replace=True)
        bb = rng_b.choice(b, size=n_b, replace=True)
        boot_m[i] = float(aa.mean() - bb.mean())
        boot_d[i] = float(np.median(aa) - np.median(bb))
    rec["delta_mean_ci_lo"] = float(np.percentile(boot_m, 2.5))
    rec["delta_mean_ci_hi"] = float(np.percentile(boot_m, 97.5))
    rec["delta_median_ci_lo"] = float(np.percentile(boot_d, 2.5))
    rec["delta_median_ci_hi"] = float(np.percentile(boot_d, 97.5))
    return rec


def _classify(md_lower, ruler_lower):
    if md_lower and ruler_lower:
        return "pr_below_fibroblast"
    if md_lower != ruler_lower:
        return "instruments_disagree"
    return "pr_not_below_fibroblast"


def _reading_text(key, line, day, md_row, ru_row, scheme):
    md_d = md_row.get("delta_mean")
    ru_d = ru_row.get("delta_mean")
    md_p = md_row.get("p_mean_lower")
    ru_p = ru_row.get("p_mean_lower")
    n_pr = md_row.get("n_a")
    n_fib = md_row.get("n_b")
    bits = (
        f"{line} scheme={scheme} d{day} PartialReprog vs Fibroblast n_PR={n_pr} n_Fib={n_fib} "
        f"MD Δmean={md_d} p_lower={md_p} Cohen_d={md_row.get('cohens_d')} "
        f"ruler Δmean={ru_d} p_lower={ru_p} Cohen_d={ru_row.get('cohens_d')}."
    )
    if key == "pr_below_fibroblast":
        return (
            "PartialReprog below Fibroblast on the frozen ruler, beating its null. "
            "Responding cells read younger than their own starting population. "
            "The rejuvenation reading holds and is not an artifact of NonReprog reading old. "
            + bits
        )
    if key == "instruments_disagree":
        return (
            "Instruments disagree (MD below Fibroblast, ruler not, or the reverse). "
            "No winner is declared. What would settle it: " + WHAT_WOULD_SETTLE + " "
            + bits
        )
    return (
        "PartialReprog not below Fibroblast while FINDINGS_MD4.md Task 1 "
        "PartialReprog < NonReprog stands. The md4 result is driven by NonReprog "
        "reading old, not by PartialReprog reading young. The rejuvenation reading "
        "is not supported by this data, and md4's reading must be restated in those terms. "
        "Quoted from FINDINGS_MD4.md: "
        f"\"{MD4_TASK1_YOUNGER_CLAIM}\" "
        + bits
    )


def _blank_inst_row(base, inst, reason):
    return dict(**base, instrument=inst, ok=False, reason=reason)


def run_scheme_a(log, obs_by, scores_by, q_a_df):
    rows = []
    keys = []
    for line in (AGED_LINE, YOUNG_LINE):
        obs = obs_by[line]
        sc = scores_by[line]
        for d in DAYS:
            q = q_a_df[(q_a_df.cell_line == line) & (q_a_df.day == int(d))].iloc[0]
            day_m = obs.day.astype(int) == int(d)
            a, b = PRIMARY_PAIR
            n_a, n_b = int(q.n_a), int(q.n_b)
            qualifies = bool(q.qualifies)
            base = dict(
                scheme="a", cell_line=line, day=int(d), state_a=a, state_b=b,
                n_a=n_a, n_b=n_b, qualifies=qualifies, primary=True,
                min_cells=MIN_CELLS_A, weaker_design=False,
                permutation="md4_task1._one_contrast_unstratified_within_timepoint",
            )
            if not qualifies:
                for inst in INSTRUMENTS:
                    rows.append(_blank_inst_row(base, inst, str(q.reason)))
                continue
            sub = day_m & obs.label.isin([a, b])
            idx = np.flatnonzero(sub.to_numpy())
            labels_sub = obs.loc[sub, "label"].to_numpy()
            mask_a = labels_sub == a
            pair_inst = {}
            for inst in INSTRUMENTS:
                x = np.asarray(sc[inst], float)[idx]
                rec = _one_contrast(
                    x, mask_a, N_PERM, N_BOOT,
                    md4_perm_seed(line, d, PRIMARY_PAIR),
                    md4_boot_seed(line, d, PRIMARY_PAIR),
                )
                rec.update(base)
                rec["instrument"] = inst
                rows.append(rec)
                pair_inst[inst] = rec
                log(f"[t1 a] {line} d{d} {inst} Δmean={rec.get('delta_mean')} "
                    f"p_lower={rec.get('p_mean_lower')} d={rec.get('cohens_d')} "
                    f"beats={rec.get('beats_null_lower')}")
            md_row = pair_inst.get("md_score") or {}
            ru_row = pair_inst.get("frozen_ruler") or {}
            key = _classify(bool(md_row.get("beats_null_lower")), bool(ru_row.get("beats_null_lower")))
            keys.append(dict(
                scheme="a", cell_line=line, day=int(d), key=key,
                md={k: md_row.get(k) for k in (
                    "n_a", "n_b", "delta_mean", "p_mean_lower", "p_mean_two_sided",
                    "delta_mean_ci_lo", "delta_mean_ci_hi", "cohens_d", "beats_null_lower",
                    "mean_a", "mean_b",
                )},
                ruler={k: ru_row.get(k) for k in (
                    "n_a", "n_b", "delta_mean", "p_mean_lower", "p_mean_two_sided",
                    "delta_mean_ci_lo", "delta_mean_ci_hi", "cohens_d", "beats_null_lower",
                    "mean_a", "mean_b",
                )},
                text=_reading_text(key, line, d, md_row, ru_row, "a"),
            ))
    df = pd.DataFrame(rows)
    df.to_csv(MD5_DIR / "t1_contrasts_a.csv", index=False)
    return df, keys


def run_scheme_b(log, obs_by, scores_by, q_b_df):
    rows = []
    keys = []
    a, b = PRIMARY_PAIR
    for line in (AGED_LINE, YOUNG_LINE):
        obs = obs_by[line]
        sc = scores_by[line]
        q = q_b_df[q_b_df.cell_line == line].iloc[0]
        fib_m = (obs.label == "Fibroblast") & obs.day.astype(int).isin(list(SCHEME_B_FIB_DAYS))
        pr_m = (obs.label == "PartialReprog") & obs.day.astype(int).isin(list(SCHEME_B_PR_DAYS))
        sub = fib_m | pr_m
        n_a = int(pr_m.sum())
        n_b = int(fib_m.sum())
        qualifies = bool(q.qualifies)
        base = dict(
            scheme="b", cell_line=line, day="d0+d3_Fib_vs_d3+d7_PR",
            state_a=a, state_b=b, n_a=n_a, n_b=n_b,
            qualifies=qualifies, primary=bool(q.primary),
            min_cells=MIN_CELLS_A, weaker_design=True,
            permutation="state_label_shuffled_within_timepoint_strata",
            n_Fibroblast_d0=int(q.n_Fibroblast_d0), n_Fibroblast_d3=int(q.n_Fibroblast_d3),
            n_PartialReprog_d3=int(q.n_PartialReprog_d3), n_PartialReprog_d7=int(q.n_PartialReprog_d7),
        )
        if not qualifies:
            for inst in INSTRUMENTS:
                rows.append(_blank_inst_row(base, inst, str(q.reason)))
            continue
        idx = np.flatnonzero(sub.to_numpy())
        labels_sub = obs.loc[sub, "label"].to_numpy()
        days_sub = obs.loc[sub, "day"].to_numpy(int)
        mask_a = labels_sub == a
        pair_inst = {}
        for inst in INSTRUMENTS:
            x = np.asarray(sc[inst], float)[idx]
            rec = _one_contrast_stratified(
                x, mask_a, days_sub, N_PERM, N_BOOT,
                perm_seed_b(line), boot_seed_b(line),
            )
            rec.update(base)
            rec["instrument"] = inst
            rows.append(rec)
            pair_inst[inst] = rec
            log(f"[t1 b] {line} {inst} Δmean={rec.get('delta_mean')} "
                f"p_lower={rec.get('p_mean_lower')} d={rec.get('cohens_d')} "
                f"beats={rec.get('beats_null_lower')} mixed_strata={rec.get('n_strata_mixed')}")
        md_row = pair_inst.get("md_score") or {}
        ru_row = pair_inst.get("frozen_ruler") or {}
        key = _classify(bool(md_row.get("beats_null_lower")), bool(ru_row.get("beats_null_lower")))
        keys.append(dict(
            scheme="b", cell_line=line, day="pooled_d0d3_Fib_vs_d3d7_PR", key=key,
            primary=bool(q.primary),
            md={k: md_row.get(k) for k in (
                "n_a", "n_b", "delta_mean", "p_mean_lower", "p_mean_two_sided",
                "delta_mean_ci_lo", "delta_mean_ci_hi", "cohens_d", "beats_null_lower",
                "mean_a", "mean_b", "n_strata_mixed",
            )},
            ruler={k: ru_row.get(k) for k in (
                "n_a", "n_b", "delta_mean", "p_mean_lower", "p_mean_two_sided",
                "delta_mean_ci_lo", "delta_mean_ci_hi", "cohens_d", "beats_null_lower",
                "mean_a", "mean_b", "n_strata_mixed",
            )},
            text=_reading_text(key, line, "pooled", md_row, ru_row, "b"),
        ))
    df = pd.DataFrame(rows)
    df.to_csv(MD5_DIR / "t1_contrasts_b.csv", index=False)
    return df, keys


def _verify_md4_overlap(log, cdf_a):
    md4 = pd.read_csv(MD4_DIR / "t1_contrasts.csv")
    log_columns("md4_t1_contrasts", list(md4.columns), str(MD4_DIR / "t1_contrasts.csv"))
    hit = md4[
        (md4.state_a == "PartialReprog") & (md4.state_b == "Fibroblast")
        & (md4.qualifies == True)  # noqa: E712
    ]
    mismatches = []
    for _, r in hit.iterrows():
        ours = cdf_a[
            (cdf_a.cell_line.astype(str) == str(r.cell_line))
            & (cdf_a.day.astype(int) == int(r.day))
            & (cdf_a.instrument.astype(str) == str(r.instrument))
            & (cdf_a.qualifies == True)  # noqa: E712
        ]
        if not len(ours):
            mismatches.append(dict(
                cell_line=str(r.cell_line), day=int(r.day), instrument=str(r.instrument),
                col="missing_scheme_a_row", md5=None, md4=float(r["delta_mean"]) if "delta_mean" in r.index else None,
            ))
            continue
        o = ours.iloc[0]
        for col in ("delta_mean", "p_mean_lower", "cohens_d"):
            if col not in r.index or col not in o.index:
                continue
            a, b = float(o[col]), float(r[col])
            if np.isfinite(a) and np.isfinite(b) and abs(a - b) > 1e-12:
                mismatches.append(dict(
                    cell_line=str(r.cell_line), day=int(r.day), instrument=str(r.instrument),
                    col=col, md5=a, md4=b,
                ))
    rec = dict(n_md4_pr_fib_qualifying=int(len(hit)), n_mismatch=len(mismatches), mismatches=mismatches)
    dump_json(MD5_DIR / "t1_md4_overlap_check.json", jsonable(rec))
    log(f"[t1 overlap] md4 PR vs Fib qualifying rows={len(hit)} mismatches={len(mismatches)}")
    if mismatches:
        raise StopStep(
            "sanity",
            "scheme (a) did not reproduce md4 PartialReprog vs Fibroblast context rows. "
            "Not using drifted scores.",
            details=rec,
        )
    return rec


def _donor_reading(line, keys_a, keys_b, q_a_df, q_b_df):
    a_hits = [h for h in keys_a if h["cell_line"] == line]
    b_hits = [h for h in keys_b if h["cell_line"] == line]
    a_ok = bool(q_a_df[(q_a_df.cell_line == line) & (q_a_df.qualifies)].shape[0] > 0)
    b_ok = bool(q_b_df[(q_b_df.cell_line == line) & (q_b_df.qualifies)].shape[0] > 0)
    if a_ok:
        used = a_hits
        scheme = "a"
        primary = True
    elif b_ok:
        used = b_hits
        scheme = "b"
        primary = True
    else:
        used = []
        scheme = None
        primary = False
    if line == AGED_LINE and not a_ok and not b_ok:
        return dict(
            key="neither_scheme_qualifies",
            scheme_used=None,
            n_qualifying_a=0, n_qualifying_b=0,
            timepoints=[],
            text=(
                f"{line}: neither scheme qualifies. The decisive test cannot be run in this "
                "dataset. The threshold is not lowered further. The young donor does not "
                "stand in for the aged one."
            ),
        )
    if not used:
        return dict(
            key="no_qualifying_timepoint",
            scheme_used=scheme,
            n_qualifying_a=int(q_a_df[(q_a_df.cell_line == line) & (q_a_df.qualifies)].shape[0]),
            n_qualifying_b=int(b_ok),
            timepoints=[],
            text=f"{line}: no qualifying PartialReprog vs Fibroblast contrast.",
        )
    keys = sorted({h["key"] for h in used})
    rec = dict(
        scheme_used=scheme,
        scheme_a_is_primary=bool(a_ok),
        n_qualifying_a=len(a_hits),
        n_qualifying_b=len(b_hits),
        timepoints=[h["day"] for h in used],
        per_timepoint=[{k: h[k] for k in ("day", "key", "text", "scheme")} for h in used],
        scheme_b_tabulated_not_primary=bool(a_ok and b_ok),
    )
    if len(keys) == 1:
        rec["key"] = keys[0]
        rec["text"] = " ".join(h["text"] for h in used)
    else:
        rec["key"] = "timepoints_disagree_do_not_pick"
        rec["keys"] = keys
        rec["text"] = (
            "Qualifying scheme-(a) timepoints fire different bullets; reported separately, "
            "not picked, not averaged. " + " ".join(h["text"] for h in used)
        )
    rec["primary"] = primary
    return rec


def run_task1(log=None):
    close_log = False
    if log is None:
        log = Logger(MD5_DIR / "t1_report.txt")
        close_log = True
    md5_log_banner(log, "TASK1")
    log(PREREG_TASK1)
    _require_prereg()
    inherit_genesets(log)
    log_id_types_louvain(log)
    progress_snapshot("Task 1 occupancy from md2 labels (no new scores)", stop="GENESETS_DONE")

    obs_by, md_by, lab_by = {}, {}, {}
    for line in (AGED_LINE, YOUNG_LINE):
        lab_by[line] = _load_labels(line, log)
        obs_by[line] = _load_obs_labelled(line, lab_by[line], log)
        md_by[line] = _load_ams_md(line, obs_by[line], log)

    counts, pooled, fib_pr, q_a, q_b, comp_b = occupancy(log, obs_by)
    sanity = sanity_md2(log, obs_by, md_by)
    if not DECLARED_BEFORE_SCORES_FLAG.exists():
        DECLARED_BEFORE_SCORES_FLAG.write_text(
            PREREG_TASK1
            + "\n\nThis flag was written after cell-count tables, scheme (a)/(b) qualify tables, "
              "and the md2 MD-mean sanity check, and before any MD5 contrast p-value or "
              "loading of md4 cell_scores for PartialReprog vs Fibroblast.\n"
            + f"\nsanity_ok={sanity['ok']}\n",
            encoding="utf-8",
        )
        log(f"[prereg] wrote {DECLARED_BEFORE_SCORES_FLAG}")
    else:
        log(f"[prereg] {DECLARED_BEFORE_SCORES_FLAG.name} already exists; not rewriting")
    progress_snapshot(
        "load md4 cell_scores; scheme (a) then scheme (b) contrasts",
        stop="DECLARED_BEFORE_SCORES",
    )

    scores_by = {}
    for line in (AGED_LINE, YOUNG_LINE):
        scores_by[line] = load_md4_cell_scores(line, obs_by[line], log)

    cdf_a, keys_a = run_scheme_a(log, obs_by, scores_by, q_a)
    _verify_md4_overlap(log, cdf_a)
    cdf_b, keys_b = run_scheme_b(log, obs_by, scores_by, q_b)

    by_donor = {
        line: _donor_reading(line, keys_a, keys_b, q_a, q_b)
        for line in (AGED_LINE, YOUNG_LINE)
    }
    reading = dict(
        by_donor=by_donor,
        primary_pair=list(PRIMARY_PAIR),
        min_cells_a=MIN_CELLS_A,
        scheme_b_fib_days=list(SCHEME_B_FIB_DAYS),
        scheme_b_pr_days=list(SCHEME_B_PR_DAYS),
        md4_sentence_at_issue=MD4_TASK1_YOUNGER_CLAIM,
        what_would_settle=WHAT_WOULD_SETTLE,
        effect_size="cohens_d_pooled_nminus1",
        delta_definition="mean(PartialReprog) - mean(Fibroblast)",
        lower_rule="delta_mean < 0 and p_mean_lower <= 0.05",
        scheme_b_is_weaker=True,
        cells_rescored=False,
        contrast_function="md4_task1._one_contrast for scheme (a); stratified variant for scheme (b)",
    )
    dump_json(MD5_DIR / "t1_reading.json", jsonable(reading))
    dump_json(MD5_DIR / "t1_summary.json", jsonable(dict(
        n_contrast_rows_a=int(len(cdf_a)),
        n_contrast_rows_b=int(len(cdf_b)),
        n_qualifying_a=int(q_a.qualifies.sum()) if len(q_a) else 0,
        n_qualifying_b=int(q_b.qualifies.sum()) if len(q_b) else 0,
        sanity_ok=bool(sanity["ok"]),
        by_donor={k: v.get("key") for k, v in by_donor.items()},
        scheme_used={k: v.get("scheme_used") for k, v in by_donor.items()},
    )))
    log("[t1 reading] " + " ".join(f"{k}={v.get('key')} scheme={v.get('scheme_used')}"
                                   for k, v in by_donor.items()))
    man = load_manifest()
    man["status"] = "TASK1_DONE"
    save_manifest(man)
    progress_snapshot("Task 2 QC / extrap / NR − Fib (report-only)", stop="TASK1_DONE")
    if close_log:
        log.close()
    return reading


if __name__ == "__main__":
    try:
        run_task1()
    except StopStep as e:
        record_failure(e.step, e.message, e.details)
        progress_snapshot(f"fix STOP [{e.step}]", stop=f"STOP {e.step}: {e.message}")
        raise
