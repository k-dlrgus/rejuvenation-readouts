"""NEWSTORY — numbers and two figures for the restructured paper.

Part A: the published MD score (Seurat AddModuleScore, src/md2_score.py) per md2 cell state in
GSE297234, scored jointly over both donors so aged and young cells share one set of expression
bins and control genes.
Part B: extra point-estimate geometry from the GENESPACE panels (functions from
src/genespace_run.py; nothing stored there is rewritten).
Part C: Figures N1, N2 and N5.

Writes only results/newstory/ and NUMBERS_NEWSTORY.md. paper/ and paper_package/ are not read.

Usage: python src/newstory.py [a|b|figs|numbers|all]
"""
from __future__ import annotations

import gc
import json
import sys
import time
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import sparse

sys.path.insert(0, str(Path(__file__).resolve().parent))
from config import RESULTS, ROOT, DATA_PROC  # noqa: E402
from md2_common import AMS_NBIN, AMS_CTRL, MD2_SEED  # noqa: E402
from md2_score import lognormalize_csr, add_module_score, symbol_index, _cut_number  # noqa: E402
from md3_idtype import csr_from_npz  # noqa: E402
from fibro3_common import AGED_LINE, YOUNG_LINE  # noqa: E402
from fibro2_common import load_frozen_ruler  # noqa: E402
import genespace_run as gs  # noqa: E402
from toward_run import load_donor, load_gtex_z, donor_panel, zscore_frozen, sum_rows  # noqa: E402
from same_run import (  # noqa: E402
    check_gene_space, split_half, d0_fibroblast_idx, state_idx, mix_indices, mix_sum, panel_z,
)

OUT = RESULTS / "newstory"
NUMBERS_PATH = ROOT / "NUMBERS_NEWSTORY.md"
MD2_DIR = RESULTS / "md2"
MD2_PROC = DATA_PROC / "md2"
MD3_GENESETS = RESULTS / "md3" / "genesets.json"
MD4_MEANS = RESULTS / "md4" / "t1_state_timepoint_means.csv"
GS_DIR = RESULTS / "genespace"
TOWARD_STATS = RESULTS / "toward" / "stats.csv"
FIBRO3_EXTRAP = RESULTS / "fibro3" / "t2_extrap.csv"
FIBRO3_CORR = RESULTS / "fibro3" / "t2_correlation.csv"

BOOT_SEED = 20260918
N_BOOT = 200
LINES = (AGED_LINE, YOUNG_LINE)
DONOR_WORDS = {AGED_LINE: "aged", YOUNG_LINE: "young"}
STATES = ("Fibroblast_d0", "PartialReprog", "EarlyPluripotency", "Pluripotency", "NonReprog")
F_VALUES = (0.0, 0.10, 0.25, 0.50)
SPACES = ("FULL", "MD", "AGE")
RAND_SPACES = ("MD", "AGE")
RAND_MODE = "A1_mu_zero"

OI = dict(blue="#0072B2", orange="#D55E00", pink="#CC79A7", grey="#999999",
          lightgrey="#C8C8C8", darkgrey="#555555", black="#000000")
SPACE_COLOR = dict(MD=OI["blue"], AGE=OI["orange"], FULL=OI["pink"])
SPACE_WORDS = dict(MD="MD genes (205)", AGE="Age genes (3,089)", FULL="All genes (23,485)")
STATE_WORDS = dict(Fibroblast_d0="Day 0", PartialReprog="Partial\nreprogramming",
                   EarlyPluripotency="Early\npluripotency", Pluripotency="Pluripotency",
                   NonReprog="Non-\nreprogrammed")


class Log:
    def __init__(self, path):
        self.f = open(path, "a", encoding="utf-8")
        self("=" * 90)
        self(f"newstory run {time.strftime('%Y-%m-%d %H:%M:%S')}")

    def __call__(self, msg):
        print(msg, flush=True)
        self.f.write(str(msg) + "\n")
        self.f.flush()

    def close(self):
        self.f.close()


# --------------------------------------------------------------------------- integrity
def snapshot():
    """(size, mtime) of every file outside results/newstory that this script could touch."""
    roots = [RESULTS, ROOT / "src", DATA_PROC / "md2", DATA_PROC / "md4"]
    snap = {}
    for r in roots:
        for p in r.rglob("*"):
            if OUT in p.parents or p == OUT or "__pycache__" in p.parts or p.name == "newstory.py":
                continue
            st = p.stat()
            snap[str(p.relative_to(ROOT))] = (p.is_dir(), st.st_size if p.is_file() else 0,
                                              st.st_mtime_ns if p.is_file() else 0)
    for p in ROOT.glob("*.md"):
        if p != NUMBERS_PATH:
            st = p.stat()
            snap[str(p.relative_to(ROOT))] = (False, st.st_size, st.st_mtime_ns)
    return snap


def snapshot_diff(a, b):
    return sorted(k for k in set(a) | set(b) if a.get(k) != b.get(k))


# --------------------------------------------------------------------------- Part A
def load_labelled_obs(line):
    obs = pd.read_csv(MD2_DIR / f"louvain_obs_{line}.csv").reset_index(drop=True)
    lab = pd.read_csv(MD2_DIR / f"t2_cluster_labels_{line}.csv")
    lab_map = {int(c): str(s) for c, s in zip(lab["cluster"], lab["label"])}
    obs["label"] = obs["cluster"].map(lambda c: lab_map.get(int(c)))
    if obs["label"].isna().any():
        raise RuntimeError(f"{line}: clusters without md2 labels")
    return obs


def state_masks(obs):
    lab = obs["label"].astype(str).to_numpy()
    day = obs["day"].astype(int).to_numpy()
    out = {"Fibroblast_d0": (lab == "Fibroblast") & (day == 0)}
    for st in STATES[1:]:
        out[st] = lab == st
    return out


def load_lognorm(line):
    X, symbols, gene_id, _ = csr_from_npz(MD2_PROC / f"louvain_counts_{line}.npz")
    logX = lognormalize_csr(X)
    del X
    gc.collect()
    return logX, np.asarray(symbols).astype(str), np.asarray(gene_id).astype(str)


def boot_mean_draws(x, rng):
    x = np.asarray(x, float)
    return np.array([rng.choice(x, size=x.size, replace=True).mean() for _ in range(N_BOOT)])


def md4_weighted_means():
    df = pd.read_csv(MD4_MEANS)
    rows = []
    for line in LINES:
        d = df[df.cell_line == line]
        for st in STATES:
            if st == "Fibroblast_d0":
                s = d[(d.label == "Fibroblast") & (d.day == 0)]
            else:
                s = d[(d.label == st) & (d.n_cells > 0)]
            n = int(s.n_cells.sum())
            m = float((s.n_cells * s.md_score_mean).sum() / n) if n else np.nan
            rows.append(dict(cell_line=line, state=st, n_cells_md4=n, md4_nweighted_mean=m))
    return pd.DataFrame(rows)


def part_a(log):
    log("[A] PART A — MD score per state, joint AddModuleScore over both donors")
    md_list = [str(g).upper() for g in json.loads(MD3_GENESETS.read_text(encoding="utf-8"))["MD"]]
    log(f"[A] MD genes from {MD3_GENESETS.relative_to(ROOT)}: n={len(md_list)} unique={len(set(md_list))}")

    obs = {line: load_labelled_obs(line) for line in LINES}
    masks = {line: state_masks(obs[line]) for line in LINES}
    for line in LINES:
        log(f"[A] {line} n_cells={len(obs[line])} " +
            " ".join(f"{st}={int(m.sum())}" for st, m in masks[line].items()))

    # Pass 1: per-donor rescoring (reproduces the stored md2 scores) + column sums for the joint mean.
    per_donor, colsum, genes, n_cells, repro = {}, {}, {}, {}, {}
    for line in LINES:
        t0 = time.time()
        logX, sym, gid = load_lognorm(line)
        if logX.shape[0] != len(obs[line]):
            raise RuntimeError(f"{line}: counts rows {logX.shape[0]} != obs rows {len(obs[line])}")
        ams, meta = add_module_score(logX, sym, {"MD": md_list}, seed=MD2_SEED, tag=f"{line}_perdonor")
        stored = np.load(MD2_PROC / f"louvain_ams_{line}.npz", allow_pickle=True)
        diff = float(np.max(np.abs(ams["MD"] - np.asarray(stored["md"], float))))
        bins_eq = bool(np.array_equal(np.asarray(meta["bins"]), np.asarray(stored["bins"])))
        repro[line] = dict(max_abs_diff_vs_stored_md=diff, bins_equal_stored=bins_eq,
                           n_ctrl=meta["per_set"]["MD"]["n_ctrl"], n_genes=int(logX.shape[1]))
        log(f"[A] {line} per-donor rescoring vs data/processed/md2/louvain_ams_{line}.npz: "
            f"max|diff|={diff:.3e} bins_equal={bins_eq} n_genes={logX.shape[1]} "
            f"n_ctrl={repro[line]['n_ctrl']} ({time.time() - t0:.0f}s)")
        per_donor[line] = ams["MD"]
        colsum[line] = np.asarray(logX.sum(axis=0), np.float64).ravel()
        genes[line] = (sym, gid)
        n_cells[line] = int(logX.shape[0])
        del logX, ams, stored
        gc.collect()

    # Union gene axis keyed on Ensembl id: aged order, then young-only genes appended.
    sym_a, gid_a = genes[AGED_LINE]
    sym_y, gid_y = genes[YOUNG_LINE]
    pos_a = {g: i for i, g in enumerate(gid_a)}
    extra = [j for j, g in enumerate(gid_y) if g not in pos_a]
    u_gid = np.concatenate([gid_a, gid_y[extra]])
    u_sym = np.concatenate([sym_a, sym_y[extra]])
    u_pos = {g: i for i, g in enumerate(u_gid)}
    colmap = {AGED_LINE: np.arange(gid_a.size), YOUNG_LINE: np.array([u_pos[g] for g in gid_y])}
    shared_j = [j for j, g in enumerate(gid_y) if g in pos_a]
    shared = [gid_y[j] for j in shared_j]
    sym_mismatch = int(np.sum(np.char.upper(sym_a[[pos_a[g] for g in shared]]) != np.char.upper(sym_y[shared_j])))
    log(f"[A] union genes={u_gid.size} (aged {gid_a.size}, young {gid_y.size}, shared {len(shared)}, "
        f"aged-only {gid_a.size - len(shared)}, young-only {len(extra)}); shared ids with different symbol={sym_mismatch}")

    n_tot = n_cells[AGED_LINE] + n_cells[YOUNG_LINE]
    tot = np.zeros(u_gid.size)
    for line in LINES:
        tot[colmap[line]] += colsum[line]
    joint_mean = tot / n_tot
    joint_bins = _cut_number(joint_mean, AMS_NBIN, np.random.default_rng(int(MD2_SEED)))

    # MD feature columns on the union axis must be the same Ensembl ids each donor used.
    u_index = symbol_index(u_sym)
    feat_u = [u_index[g] for g in md_list if g in u_index]
    feat_check = {}
    for line, (sym, gid) in genes.items():
        d_index = symbol_index(sym)
        same = sum(1 for g in md_list if g in d_index and g in u_index and gid[d_index[g]] == u_gid[u_index[g]])
        feat_check[line] = dict(n_md_in_donor=sum(1 for g in md_list if g in d_index), n_same_ensembl=same,
                                md_absent_from_donor_matrix=[g for g in md_list if g not in d_index])
    log(f"[A] MD features on union axis n={len(feat_u)}; same Ensembl id as per-donor mapping: {feat_check}")

    # Control genes, drawn once exactly as add_module_score draws them (Generator(seed + 1)).
    rng_ctrl = np.random.default_rng(int(MD2_SEED) + 1)
    ctrl = []
    for j in feat_u:
        pool = np.flatnonzero(joint_bins == joint_bins[j])
        ctrl.extend(rng_ctrl.choice(pool, size=int(AMS_CTRL), replace=False).tolist())
    ctrl_u = np.unique(ctrl)

    # Pass 2: score each donor's rows on the union axis with the joint means / bins.
    joint = {}
    for line in LINES:
        t0 = time.time()
        logX, sym, gid = load_lognorm(line)
        cm = colmap[line].astype(np.int32)
        Xu = sparse.csr_matrix((logX.data, cm[logX.indices], logX.indptr), shape=(logX.shape[0], u_gid.size))
        del logX
        gc.collect()
        ams, meta = add_module_score(Xu, u_sym, {"MD": md_list}, seed=MD2_SEED, tag=f"{line}_joint",
                                     bins=joint_bins, gene_mean=joint_mean)
        direct = (np.asarray(Xu[:, feat_u].mean(axis=1)).ravel() - np.asarray(Xu[:, ctrl_u].mean(axis=1)).ravel())
        dchk = float(np.max(np.abs(direct - ams["MD"])))
        log(f"[A] {line} joint score: n_ctrl={meta['per_set']['MD']['n_ctrl']} (shared set n={ctrl_u.size}); "
            f"max|AMS - direct(shared ctrl)|={dchk:.3e} ({time.time() - t0:.0f}s)")
        if meta["per_set"]["MD"]["n_ctrl"] != ctrl_u.size or dchk > 1e-12:
            raise RuntimeError(f"{line}: joint control set is not the shared set")
        joint[line] = ams["MD"]
        del Xu, ams
        gc.collect()

    np.savez_compressed(OUT / "partA_cell_scores.npz",
                        **{f"joint__{l}": joint[l] for l in LINES},
                        **{f"perdonor__{l}": per_donor[l] for l in LINES},
                        **{f"label__{l}": obs[l]["label"].to_numpy().astype(str) for l in LINES},
                        **{f"day__{l}": obs[l]["day"].to_numpy(int) for l in LINES},
                        ctrl_union_cols=ctrl_u, feat_union_cols=np.asarray(feat_u),
                        union_gene_id=u_gid, union_symbol=u_sym)

    # Per-state statistics. One Generator(BOOT_SEED); donors then states in fixed order.
    rng = np.random.default_rng(BOOT_SEED)
    rows, draws = [], {}
    md4 = md4_weighted_means()
    for line in LINES:
        for st in STATES:
            m = masks[line][st]
            x = joint[line][m]
            b = boot_mean_draws(x, rng)
            draws[(line, st)] = b
            r4 = md4[(md4.cell_line == line) & (md4.state == st)].iloc[0]
            rows.append(dict(
                cell_line=line, donor=DONOR_WORDS[line], state=st, n_cells=int(m.sum()),
                joint_mean=float(x.mean()), joint_ci_lo=float(np.percentile(b, 2.5)),
                joint_ci_hi=float(np.percentile(b, 97.5)),
                perdonor_mean=float(per_donor[line][m].mean()),
                n_cells_md4=int(r4.n_cells_md4), md4_nweighted_mean=float(r4.md4_nweighted_mean),
            ))
            log(f"[A] {line} {st:17s} n={int(m.sum()):6d} joint={x.mean():+.4f} "
                f"CI=[{np.percentile(b, 2.5):+.4f}, {np.percentile(b, 97.5):+.4f}] "
                f"per-donor={per_donor[line][m].mean():+.4f} md4={r4.md4_nweighted_mean:+.4f}")
    st_df = pd.DataFrame(rows)
    st_df.to_csv(OUT / "partA_md_by_state.csv", index=False)

    def get(col, line, st):
        return float(st_df[(st_df.cell_line == line) & (st_df.state == st)][col].iloc[0])

    grows = []
    for scoring, col in (("joint", "joint_mean"), ("per-donor (md2 stored)", "perdonor_mean"),
                         ("md4 n-weighted", "md4_nweighted_mean")):
        a0, y0 = get(col, AGED_LINE, "Fibroblast_d0"), get(col, YOUNG_LINE, "Fibroblast_d0")
        apr = get(col, AGED_LINE, "PartialReprog")
        gap, drop = a0 - y0, a0 - apr
        rec = dict(scoring=scoring, aged_d0=a0, young_d0=y0, aged_PartialReprog=apr,
                   gap=gap, drop=drop, drop_over_gap=drop / gap if gap != 0 else np.nan)
        if scoring == "joint":
            bg = draws[(AGED_LINE, "Fibroblast_d0")] - draws[(YOUNG_LINE, "Fibroblast_d0")]
            bd = draws[(AGED_LINE, "Fibroblast_d0")] - draws[(AGED_LINE, "PartialReprog")]
            br = bd / bg
            rec.update(gap_ci_lo=np.percentile(bg, 2.5), gap_ci_hi=np.percentile(bg, 97.5),
                       drop_ci_lo=np.percentile(bd, 2.5), drop_ci_hi=np.percentile(bd, 97.5),
                       ratio_ci_lo=np.percentile(br, 2.5), ratio_ci_hi=np.percentile(br, 97.5))
        grows.append(rec)
        log(f"[A] {scoring}: gap={gap:+.4f} drop={drop:+.4f} drop/gap={rec['drop_over_gap']:.2f}")
    pd.DataFrame(grows).to_csv(OUT / "partA_gap_drop.csv", index=False)

    meta = dict(
        md_genes=len(md_list), nbin=AMS_NBIN, ctrl=AMS_CTRL, seed=int(MD2_SEED), boot_seed=BOOT_SEED,
        n_boot=N_BOOT, union_genes=int(u_gid.size), shared_genes=len(shared), young_only=len(extra),
        aged_only=int(gid_a.size - len(shared)), n_ctrl_joint=int(ctrl_u.size),
        n_cells=n_cells, perdonor_repro=repro, feature_check=feat_check,
        md4_scoring=("md_score_mean in results/md4/t1_state_timepoint_means.csv is the md2 "
                     "louvain_ams vector, computed by add_module_score on each donor's own Louvain "
                     "matrix separately (src/md2_task1.py for GM00731, src/md2_task2.py for GM23815); "
                     "gene means, bins and control genes differ between donors."),
    )
    (OUT / "partA_meta.json").write_text(json.dumps(meta, indent=2, default=str), encoding="utf-8")
    return meta


# --------------------------------------------------------------------------- Part B
def build_geometry_ctx(log):
    frozen_npz = load_frozen_ruler()
    frozen = {k: np.asarray(frozen_npz[k]) for k in ("mu", "sd", "w", "symbol", "ensembl")}
    genesets = json.loads(gs.GENESETS_PATH.read_text(encoding="utf-8"))
    gm = gs.map_gene_spaces(frozen, genesets, log)
    spaces = gm["spaces"]
    obs_a, Y_a = load_donor(AGED_LINE, log)
    obs_y, Y_y = load_donor(YOUNG_LINE, log)
    check_gene_space(Y_a, frozen, AGED_LINE)
    check_gene_space(Y_y, frozen, YOUNG_LINE)
    sp_file = np.load(gs.SAME_DIR / "split_idx.npz")
    o_idx, y_idx = np.asarray(sp_file["O"], int), np.asarray(sp_file["Y"], int)
    aged_b, young_b = np.asarray(sp_file["aged_B"], int), np.asarray(sp_file["young_B"], int)
    ra, rb = split_half(d0_fibroblast_idx(obs_a), gs.SEED)
    ya, yb = split_half(d0_fibroblast_idx(obs_y), gs.SEED)
    if not (np.array_equal(ra, o_idx) and np.array_equal(rb, aged_b)
            and np.array_equal(ya, y_idx) and np.array_equal(yb, young_b)):
        raise RuntimeError("split_half does not equal results/same/split_idx.npz")
    main_spec = [
        ("O", Y_a, o_idx), ("Y", Y_y, y_idx), ("S", Y_a, state_idx(obs_a, "PartialReprog")),
        ("S_rev", Y_y, state_idx(obs_y, "PartialReprog")), ("Pluri", Y_a, state_idx(obs_a, "Pluripotency")),
        ("NonReprog", Y_a, state_idx(obs_a, "NonReprog")),
    ]
    n_cells = {nm: int(ix.size) for nm, _, ix in main_spec}
    main_counts = [sum_rows(Ym, ix) for _, Ym, ix in main_spec]
    zero = (np.vstack(main_counts) == 0).all(0)
    bins, _ = gs.mu_bins(frozen["mu"])
    mixes, _ = mix_indices(aged_b, young_b, gs.SEED)
    mix_counts = [sum_rows(Y_a, o_idx), sum_rows(Y_y, y_idx)]
    for m in mixes:
        mix_counts.append(mix_sum(Y_a, Y_y, m["aged_idx"], m["young_idx"]))
    Zmain, _, _, _ = panel_z(main_counts, frozen)
    Zmix, _, _, _ = panel_z(mix_counts, frozen)
    gtex = load_gtex_z(log)
    age = gtex["age"]
    anchors = dict(c_young=gtex["Z"][np.isin(age, gs.YOUNG_BINS)].mean(0),
                   c_old=gtex["Z"][np.isin(age, gs.OLD_BINS)].mean(0))
    anchor_n = dict(n_young=int(np.isin(age, gs.YOUNG_BINS).sum()), n_old=int(np.isin(age, gs.OLD_BINS).sum()))
    donor_z = {}
    for line, (obs, Ym) in {AGED_LINE: (obs_a, Y_a), YOUNG_LINE: (obs_y, Y_y)}.items():
        panel = donor_panel(obs, Ym, log, line)
        Zd, miss_d = zscore_frozen(panel["logcpm"], frozen, panel["C"])
        donor_z[line] = dict(Z=Zd, names=list(panel["names"]), idx=panel["idx"])
    del gtex
    gc.collect()
    return dict(frozen=frozen, spaces=spaces, bins=bins, zero=zero, Zmain=Zmain, Zmix=Zmix,
                mixes=mixes, n_cells=n_cells, anchors=anchors, anchor_n=anchor_n, donor_z=donor_z)


def part_b(log):
    log("[B] PART B — geometry point estimates (genespace_run functions)")
    ctx = build_geometry_ctx(log)
    spaces = ctx["spaces"]
    ref_same = pd.read_csv(GS_DIR / "same_stats.csv")
    ref_mix = pd.read_csv(GS_DIR / "same_posctrl.csv")
    ref_gtex = pd.read_csv(GS_DIR / "gtex_stats.csv")
    checks = []

    same_rows, mix_rows = [], []
    for sp in SPACES:
        main, mix = gs.same_point_block(ctx["Zmain"], ctx["Zmix"], spaces[sp], ctx["mixes"], sp)
        for test in ("forward", "reverse", "Pluri", "NonReprog"):
            st = main[test]
            same_rows.append(dict(space=sp, n_genes=int(spaces[sp].size), test=test,
                                  n_cells=ctx["n_cells"][gs.MAIN_NAMES[gs.MAIN_TESTS[test][2]]],
                                  dist_start=st["dist_start"], dist_final=st["dist_final"],
                                  final_over_start=st["dist_final"] / st["dist_start"],
                                  progress=st["progress"], perp=st["perp"], delta=st["delta"],
                                  rel_delta=st["rel_delta"], axis_length=st["v_norm"]))
            r = ref_same[(ref_same.space == sp) & (ref_same.test == test)].iloc[0]
            for k in ("dist_start", "dist_final", "progress", "perp", "delta"):
                checks.append(("same_stats.csv", f"{sp} {test} {k}", st[k], float(r[k])))
        for m, st in zip(ctx["mixes"], mix):
            mix_rows.append(dict(space=sp, f=m["f"], n_aged=m["n_aged"], n_young=m["n_young"],
                                 dist_start=st["dist_start"], dist_final=st["dist_final"],
                                 progress=st["progress"], perp=st["perp"], delta=st["delta"],
                                 rel_delta=st["rel_delta"]))
            r = ref_mix[(ref_mix.space == sp) & np.isclose(ref_mix.f, m["f"])].iloc[0]
            for k in ("progress", "perp", "delta", "rel_delta"):
                checks.append(("same_posctrl.csv", f"{sp} f={m['f']:.2f} {k}", st[k], float(r[k])))
    pd.DataFrame(same_rows).to_csv(OUT / "partB_same_distances.csv", index=False)
    pd.DataFrame(mix_rows).to_csv(OUT / "partB_mixture.csv", index=False)

    gtex_rows = []
    for sp in SPACES:
        g = gs.gtex_point_block(ctx["donor_z"], ctx["anchors"], spaces[sp])
        for line in LINES:
            for nm, rec in g[line].items():
                gtex_rows.append(dict(space=sp, n_genes=int(spaces[sp].size), cell_line=line, state=nm,
                                      n_cells=int(ctx["donor_z"][line]["idx"][nm].size),
                                      dist_d0_to_gtex_young=rec["dist_young_origin"],
                                      dist_d0_to_gtex_old=rec["dist_old_origin"],
                                      gtex_young_to_old=rec["w_norm"],
                                      young_to_old_over_d0_to_young=rec["w_norm"] / rec["dist_young_origin"],
                                      dist_state_to_gtex_young=rec["dist_young"],
                                      dist_state_to_gtex_old=rec["dist_old"],
                                      delta_young=rec["delta_young"], delta_old=rec["delta_old"],
                                      rel_delta_young=rec["rel_delta_young"], rel_delta_old=rec["rel_delta_old"]))
                r = ref_gtex[(ref_gtex.space == sp) & (ref_gtex.cell_line == line) & (ref_gtex.state == nm)].iloc[0]
                for k in ("dist_young_origin", "dist_old_origin", "w_norm", "delta_young", "delta_old"):
                    checks.append(("gtex_stats.csv", f"{sp} {line} {nm} {k}", rec[k], float(r[k])))
    pd.DataFrame(gtex_rows).to_csv(OUT / "partB_gtex_distances.csv", index=False)

    # Random gene sets (Amendment 1): stored draws, per-set point statistics.
    sets_npz = np.load(GS_DIR / "size_baseline_sets.npz")
    ref_base = pd.read_csv(GS_DIR / "size_baseline.csv")
    rand_rows, pct_rows = [], []
    for sp in RAND_SPACES:
        sets = np.asarray(sets_npz[f"{RAND_MODE}__{sp}"], int)
        redraw, _ = gs.draw_baseline_sets(spaces[sp], ctx["bins"], ctx["zero"], RAND_MODE)
        sets_equal = bool(np.array_equal(redraw, sets))
        log(f"[B] {RAND_MODE} {sp}: stored sets {sets.shape}; redrawn with draw_baseline_sets equal: {sets_equal}")
        checks.append(("size_baseline_sets.npz", f"{RAND_MODE}__{sp} redraw identical", float(sets_equal), 1.0))
        for i, s in enumerate(sets):
            main, mix = gs.same_point_block(ctx["Zmain"], ctx["Zmix"], np.sort(s), ctx["mixes"], f"{sp} set{i}")
            f50 = next(st for m, st in zip(ctx["mixes"], mix) if abs(m["f"] - 0.50) < 1e-12)
            rand_rows.append(dict(
                space_size=sp, n_genes=int(s.size), set_index=i,
                overlap_with_real=int(np.isin(s, spaces[sp]).sum()),
                fwd_progress=main["forward"]["progress"], fwd_perp=main["forward"]["perp"],
                fwd_delta=main["forward"]["delta"], fwd_rel_delta=main["forward"]["rel_delta"],
                fwd_final_over_start=main["forward"]["dist_final"] / main["forward"]["dist_start"],
                rev_progress=main["reverse"]["progress"], rev_perp=main["reverse"]["perp"],
                rev_delta=main["reverse"]["delta"], rev_rel_delta=main["reverse"]["rel_delta"],
                mix50_delta=f50["delta"], mix50_rel_delta=f50["rel_delta"], mix50_progress=f50["progress"],
            ))
        rd = pd.DataFrame([r for r in rand_rows if r["space_size"] == sp])
        for test, state, stat, col in (("same", "forward", "progress", "fwd_progress"),
                                       ("same", "forward", "perp", "fwd_perp"),
                                       ("same", "forward", "delta", "fwd_delta"),
                                       ("same", "reverse", "progress", "rev_progress"),
                                       ("same", "reverse", "perp", "rev_perp"),
                                       ("same", "reverse", "delta", "rev_delta"),
                                       ("mixture", "f=0.50", "delta", "mix50_delta")):
            r = ref_base[(ref_base.baseline == RAND_MODE) & (ref_base.space == sp) & (ref_base.test == test)
                         & (ref_base.state == state) & (ref_base.stat == stat)].iloc[0]
            v = rd[col].to_numpy(float)
            v = v[np.isfinite(v)]
            real = float(r["real"])
            pos = (float(np.sum(v < real)) + 0.5 * float(np.sum(v == real))) / gs.N_BASELINE
            rec = dict(space_size=sp, test=test, state=state, stat=stat, n_finite=int(v.size), real=real,
                       p2_5=np.percentile(v, 2.5), p50=np.percentile(v, 50), p97_5=np.percentile(v, 97.5),
                       position=pos, stored_p2_5=float(r.rand_p2_5), stored_p50=float(r.rand_p50),
                       stored_p97_5=float(r.rand_p97_5), stored_position=float(r.position))
            rec["max_abs_diff"] = max(abs(rec["p2_5"] - rec["stored_p2_5"]), abs(rec["p50"] - rec["stored_p50"]),
                                      abs(rec["p97_5"] - rec["stored_p97_5"]), abs(pos - rec["stored_position"]))
            pct_rows.append(rec)
            log(f"[B] {sp} {test} {state} {stat}: p2.5/50/97.5 = {rec['p2_5']:.6f}/{rec['p50']:.6f}/{rec['p97_5']:.6f} "
                f"max|diff vs size_baseline.csv|={rec['max_abs_diff']:.2e}")
    rand_df = pd.DataFrame(rand_rows)
    rand_df.to_csv(OUT / "partB_random_sets_per_set.csv", index=False)
    pd.DataFrame(pct_rows).to_csv(OUT / "partB_random_sets_repro.csv", index=False)

    same_df = pd.DataFrame(same_rows)
    cnt_rows = []
    for sp in RAND_SPACES:
        rd = rand_df[rand_df.space_size == sp]
        real_f = same_df[(same_df.space == sp) & (same_df.test == "forward")].iloc[0]
        real_r = same_df[(same_df.space == sp) & (same_df.test == "reverse")].iloc[0]
        cnt_rows.append(dict(
            space_size=sp, n_sets=int(len(rd)),
            fwd_closer=int((rd.fwd_delta < 0).sum()), fwd_farther=int((rd.fwd_delta > 0).sum()),
            rev_closer=int((rd.rev_delta < 0).sum()), mix50_closer=int((rd.mix50_delta < 0).sum()),
            fwd_min_final_over_start=float(rd.fwd_final_over_start.min()),
            real_fwd_final_over_start=float(real_f.final_over_start),
            real_fwd_progress=float(real_f.progress), real_fwd_perp=float(real_f.perp),
            real_rev_progress=float(real_r.progress), real_rev_perp=float(real_r.perp),
        ))
        log(f"[B] {sp}-size random sets: forward closer to young target in {cnt_rows[-1]['fwd_closer']}/{len(rd)}")
    pd.DataFrame(cnt_rows).to_csv(OUT / "partB_random_sets_counts.csv", index=False)

    ck = pd.DataFrame(checks, columns=["reference_file", "quantity", "recomputed", "stored"])
    ck["abs_diff"] = (ck.recomputed.astype(float) - ck.stored.astype(float)).abs()
    ck.to_csv(OUT / "partB_repro_checks.csv", index=False)
    log(f"[B] reproduction checks vs results/genespace: n={len(ck)} max|diff|={ck.abs_diff.max():.3e}")
    meta = dict(anchor_n=ctx["anchor_n"], n_cells=ctx["n_cells"],
                spaces={sp: int(spaces[sp].size) for sp in SPACES},
                mixes=[dict(f=m["f"], n_aged=m["n_aged"], n_young=m["n_young"]) for m in ctx["mixes"]],
                max_abs_diff_point=float(ck[ck.reference_file != "size_baseline_sets.npz"].abs_diff.max()),
                max_abs_diff_percentiles=float(pd.DataFrame(pct_rows).max_abs_diff.max()))
    (OUT / "partB_meta.json").write_text(json.dumps(meta, indent=2, default=str), encoding="utf-8")
    return meta


# --------------------------------------------------------------------------- Part C
def fig_style():
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    plt.rcParams.update({
        "font.family": "sans-serif", "font.sans-serif": ["Arial", "Helvetica", "DejaVu Sans"],
        "font.size": 7, "axes.labelsize": 7, "xtick.labelsize": 6.5, "ytick.labelsize": 6.5,
        "legend.fontsize": 6.5, "axes.spines.top": False, "axes.spines.right": False,
        "axes.linewidth": 0.6, "xtick.major.width": 0.6, "ytick.major.width": 0.6,
        "xtick.major.size": 2.5, "ytick.major.size": 2.5, "lines.linewidth": 0.9,
        "pdf.fonttype": 42, "ps.fonttype": 42, "savefig.dpi": 300,
    })
    return plt


def panel_letter(ax, s, x=-0.2):
    ax.text(x, 1.04, s, transform=ax.transAxes, fontsize=9, fontweight="bold", va="bottom", ha="left")


def aligned_panel_letters(fig, letters, x=-0.3):
    # Equal-aspect axes shrink inside their gridspec slot, so letters are placed at one figure height.
    from matplotlib.transforms import blended_transform_factory
    for ax, _ in letters:
        ax.apply_aspect()
    pos = [ax.get_position() for ax, _ in letters]
    y = max(p.y1 + 0.04 * p.height for p in pos)
    for ax, s in letters:
        fig.text(x, y, s, transform=blended_transform_factory(ax.transAxes, fig.transFigure), fontsize=9,
                 fontweight="bold", va="bottom", ha="left")


def draw_plane_frame(ax, plt, xlim, ylim, circle_label_xy=(1.3, 0.5)):
    from matplotlib.patches import Circle
    ax.add_patch(Circle((1.0, 0.0), 1.0, fill=False, ls=(0, (4, 2.5)), lw=0.8, ec=OI["grey"], zorder=1))
    if circle_label_xy:
        ax.text(*circle_label_xy, "closer than\nat start", ha="center", va="center", fontsize=6, color=OI["darkgrey"])
    ax.axhline(0, color=OI["black"], lw=0.5, zorder=0)
    ax.plot([0], [0], marker="s", ms=4.5, mfc="white", mec=OI["black"], mew=0.8, ls="none", zorder=5, clip_on=False)
    ax.plot([1], [0], marker="*", ms=7.5, mfc=OI["black"], mec=OI["black"], ls="none", zorder=5, clip_on=False)
    ax.annotate("start", (0, 0), xytext=(-4, 3), textcoords="offset points", fontsize=5.8, ha="right", va="bottom")
    ax.annotate("target", (1, 0), xytext=(4.5, 3), textcoords="offset points", fontsize=5.8, ha="left", va="bottom")
    ax.set_xlim(*xlim)
    ax.set_ylim(*ylim)
    ax.set_aspect("equal", adjustable="box")
    ax.set_xlabel("Progress along start-to-target axis\n(axis lengths)")
    ax.set_ylabel("Sideways distance (axis lengths)")


def save_fig(fig, stem):
    fig.savefig(OUT / f"{stem}.png", dpi=300)
    fig.savefig(OUT / f"{stem}.pdf")


def figure_n1(log):
    plt = fig_style()
    st = pd.read_csv(OUT / "partA_md_by_state.csv")
    same = pd.read_csv(OUT / "partB_same_distances.csv")
    mix = pd.read_csv(OUT / "partB_mixture.csv")
    data = []

    fig_w, left_in, ax_in, gap_in = 4.8, 0.525, 1.5925, 0.79625
    fig = plt.figure(figsize=(fig_w, 3.0))
    gsp = fig.add_gridspec(1, 2, wspace=gap_in / ax_in, left=left_in / fig_w,
                           right=(left_in + 2 * ax_in + gap_in) / fig_w, bottom=0.27, top=0.92)
    # A
    ax = fig.add_subplot(gsp[0])
    order = ("Fibroblast_d0", "PartialReprog", "EarlyPluripotency", "Pluripotency")
    aged = st[st.cell_line == AGED_LINE].set_index("state").loc[list(order)]
    young0 = st[(st.cell_line == YOUNG_LINE) & (st.state == "Fibroblast_d0")].iloc[0]
    x = np.arange(len(order))
    ax.axhline(young0.joint_mean, color=OI["darkgrey"], ls=(0, (4, 2.5)), lw=0.9, zorder=1)
    ax.text(len(order) - 0.55, young0.joint_mean + 0.012, "22-year-old donor, day 0", ha="right", va="bottom",
            fontsize=6, color=OI["darkgrey"])
    ax.errorbar(x, aged.joint_mean, yerr=[aged.joint_mean - aged.joint_ci_lo, aged.joint_ci_hi - aged.joint_mean],
                fmt="o", ms=4, color=OI["blue"], ecolor=OI["blue"], elinewidth=0.9, capsize=2, zorder=3)
    ax.set_xticks(x)
    ax.set_xticklabels([STATE_WORDS[s].replace("\n", " ") for s in order], rotation=30, ha="right",
                       rotation_mode="anchor")
    ax.set_xlim(-0.5, len(order) - 0.5)
    ax.set_ylim(-0.05, 0.6)
    ax.set_ylabel("MD score (mean, 95% CI)")
    ax.set_xlabel("Aged donor's cells")
    panel_letter(ax, "A", x=-0.3)
    for s, r in aged.iterrows():
        data.append(dict(panel="A", series="aged donor", label=s, x=float(order.index(s)), y=r.joint_mean,
                         y_lo=r.joint_ci_lo, y_hi=r.joint_ci_hi, n_cells=int(r.n_cells)))
    data.append(dict(panel="A", series="young donor day 0 (dashed line)", label="Fibroblast_d0", x=np.nan,
                     y=young0.joint_mean, y_lo=young0.joint_ci_lo, y_hi=young0.joint_ci_hi,
                     n_cells=int(young0.n_cells)))

    # B
    ax = fig.add_subplot(gsp[1])
    draw_plane_frame(ax, plt, (-0.75, 2.05), (0, 4.25), circle_label_xy=(1.0, 0.57))
    md = same[same.space == "MD"].set_index("test")
    mm = mix[mix.space == "MD"].sort_values("f")
    for _, r in mm.iterrows():
        ax.plot(r.progress, r.perp, marker="o", ms=3.6, ls="none", mec=OI["grey"], mew=0.8,
                mfc="white" if r.f == 0 else OI["grey"], zorder=3)
        data.append(dict(panel="B", series="mixture", label=f"f={r.f:.2f}", x=r.progress, y=r.perp))
    ax.annotate("young-cell\nmixtures\n(0–50%)", xy=(mm.progress.iloc[1], mm.perp.iloc[1]), xytext=(-0.22, 0.78),
                fontsize=5.8, color=OI["darkgrey"], ha="center", va="bottom",
                arrowprops=dict(arrowstyle="-", color=OI["grey"], lw=0.5, shrinkA=1, shrinkB=3))
    f = md.loc["forward"]
    ax.plot(f.progress, f.perp, marker="o", ms=5.5, color=OI["blue"], ls="none", zorder=4)
    ax.text(f.progress + 0.12, f.perp, "Partial\nreprogramming", fontsize=6, color=OI["blue"], va="center")
    for test, word in (("Pluri", "Pluripotency"), ("NonReprog", "Non-reprogrammed")):
        r = md.loc[test]
        ax.plot(r.progress, r.perp, marker="o", ms=4.5, mfc="white", mec=OI["black"], mew=0.8, ls="none", zorder=4)
        ax.text(r.progress + (0.1 if test == "NonReprog" else -0.1), r.perp + (-0.02 if test == "NonReprog" else 0.0),
                word, fontsize=6, ha="left" if test == "NonReprog" else "right", va="center")
    r = md.loc["reverse"]
    ax.plot(r.progress, r.perp, marker="o", ms=5.0, color=OI["orange"], ls="none", zorder=4)
    ax.text(r.progress, r.perp + 0.2, "reverse", fontsize=6, color=OI["orange"], ha="center", va="bottom")
    ax.set_xlabel("Progress toward young donor\n(axis lengths)")
    panel_letter(ax, "B", x=-0.32)
    for test in ("forward", "reverse", "Pluri", "NonReprog"):
        r = md.loc[test]
        data.append(dict(panel="B", series=test, label=test, x=r.progress, y=r.perp, n_cells=int(r.n_cells),
                         note=f"final/start={r.final_over_start:.4f}"))
    save_fig(fig, "fig_N1_md_own_genes")
    plt.close(fig)
    pd.DataFrame(data).to_csv(OUT / "fig_N1_md_own_genes_data.csv", index=False)
    log("[C] wrote fig_N1_md_own_genes.{png,pdf,_data.csv}")


def figure_n2(log):
    plt = fig_style()
    from matplotlib.lines import Line2D
    same = pd.read_csv(OUT / "partB_same_distances.csv")
    mix = pd.read_csv(OUT / "partB_mixture.csv")
    rnd = pd.read_csv(OUT / "partB_random_sets_per_set.csv")
    data = []
    rand_style = dict(MD=dict(marker="o", color=OI["lightgrey"], ms=2.2, word="random 205-gene sets (200)"),
                      AGE=dict(marker="^", color=OI["darkgrey"], ms=2.2, word="random 3,089-gene sets (200)"))

    fig = plt.figure(figsize=(7.0, 3.1))
    gsp = fig.add_gridspec(1, 3, width_ratios=[1.0, 1.0, 0.95], wspace=0.5, left=0.07, right=0.985,
                           bottom=0.3, top=0.92)
    xlim, ylim = (-1.15, 2.05), (0, 3.05)
    letters = []
    for k, (test, pre, letter) in enumerate((("forward", "fwd", "A"), ("reverse", "rev", "B"))):
        ax = fig.add_subplot(gsp[k])
        draw_plane_frame(ax, plt, xlim, ylim)
        for sp in RAND_SPACES:
            rs = rnd[rnd.space_size == sp]
            sty = rand_style[sp]
            ax.plot(rs[f"{pre}_progress"], rs[f"{pre}_perp"], ls="none", marker=sty["marker"], ms=sty["ms"],
                    mfc=sty["color"], mec="none", alpha=0.9, zorder=2)
            for _, r in rs.iterrows():
                data.append(dict(panel=letter, series=f"random {sp}-size", label=f"set {int(r.set_index)}",
                                 x=r[f"{pre}_progress"], y=r[f"{pre}_perp"], note=f"delta={r[f'{pre}_delta']:.4f}"))
        for sp in ("FULL", "AGE", "MD"):
            r = same[(same.space == sp) & (same.test == test)].iloc[0]
            ax.plot(r.progress, r.perp, marker="o", ms=5.5, color=SPACE_COLOR[sp], mec="white", mew=0.5,
                    ls="none", zorder=4)
            data.append(dict(panel=letter, series=f"real {sp}", label=SPACE_WORDS[sp], x=r.progress, y=r.perp,
                             n_cells=int(r.n_cells), note=f"final/start={r.final_over_start:.4f}"))
        if test == "forward":
            ax.set_xlabel("Progress toward young donor\n(axis lengths)")
        else:
            ax.set_xlabel("Progress toward aged donor\n(axis lengths)")
            ax.text(0.97, 0.97, "reverse", transform=ax.transAxes, fontsize=6.5, color=OI["orange"],
                    ha="right", va="top")
        letters.append((ax, letter))

    ax = fig.add_subplot(gsp[2])
    for sp in ("MD", "AGE", "FULL"):
        mm = mix[mix.space == sp].sort_values("f")
        ax.plot(mm.f, mm.rel_delta, color=SPACE_COLOR[sp], lw=0.9, zorder=2)
        for _, r in mm.iterrows():
            ax.plot(r.f, r.rel_delta, marker="o", ms=4, ls="none", mec=SPACE_COLOR[sp], mew=0.8,
                    mfc="white" if r.f == 0 else SPACE_COLOR[sp], zorder=3)
            data.append(dict(panel="C", series=f"mixture {sp}", label=f"f={r.f:.2f}", x=r.f, y=r.rel_delta,
                             note=f"delta={r.delta:.4f}"))
    ax.axhline(0, color=OI["black"], lw=0.5)
    ax.set_xticks(F_VALUES)
    ax.set_xticklabels(["0", "0.10", "0.25", "0.50"])
    ax.set_xlim(-0.03, 0.53)
    ax.set_xlabel("Fraction of young cells in mixture")
    ax.set_ylabel("Change in distance ÷ starting distance")
    letters.append((ax, "C"))
    aligned_panel_letters(fig, letters, x=-0.3)

    handles = [Line2D([], [], marker="o", ls="none", ms=5.5, color=SPACE_COLOR[sp], label=SPACE_WORDS[sp])
               for sp in ("MD", "AGE", "FULL")]
    handles += [Line2D([], [], marker=rand_style[sp]["marker"], ls="none", ms=3.5, mfc=rand_style[sp]["color"],
                       mec="none", label=rand_style[sp]["word"]) for sp in RAND_SPACES]
    handles += [Line2D([], [], marker="s", ls="none", ms=4.5, mfc="white", mec=OI["black"], label="start"),
                Line2D([], [], marker="*", ls="none", ms=7, color=OI["black"], label="target (young donor, day 0; "
                       "reverse: aged donor)"),
                Line2D([], [], ls=(0, (4, 2.5)), color=OI["grey"], label="closer than at start")]
    fig.legend(handles=handles, loc="lower center", ncol=4, frameon=False, bbox_to_anchor=(0.5, 0.0),
               columnspacing=1.2, handletextpad=0.4)
    save_fig(fig, "fig_N2_any_genes")
    plt.close(fig)
    pd.DataFrame(data).to_csv(OUT / "fig_N2_any_genes_data.csv", index=False)
    log("[C] wrote fig_N2_any_genes.{png,pdf,_data.csv}")


def figure_n5(log):
    plt = fig_style()
    from matplotlib.lines import Line2D
    from scipy.stats import spearmanr
    tw = pd.read_csv(TOWARD_STATS, float_precision="round_trip")
    ex = pd.read_csv(FIBRO3_EXTRAP, float_precision="round_trip")
    corr = pd.read_csv(FIBRO3_CORR, float_precision="round_trip")
    data = []

    fig = plt.figure(figsize=(7.0, 3.0))
    gsp = fig.add_gridspec(1, 2, width_ratios=[1.0, 1.45], wspace=0.35, left=0.075, right=0.985,
                           bottom=0.27, top=0.92)
    # A
    ax = fig.add_subplot(gsp[0])
    order = ("Fibroblast", "PartialReprog", "EarlyPluripotency", "Pluripotency")
    aged = tw[tw.cell_line == AGED_LINE].set_index("state").loc[list(order)]
    young0 = tw[(tw.cell_line == YOUNG_LINE) & (tw.state == "Fibroblast")].iloc[0]
    x = np.arange(len(order))
    ax.axhline(young0.ruler_score_context, color=OI["darkgrey"], ls=(0, (4, 2.5)), lw=0.9, zorder=1)
    ax.text(len(order) - 0.55, young0.ruler_score_context + 0.35, "22-year-old donor, day 0", ha="right",
            va="bottom", fontsize=6, color=OI["darkgrey"])
    ax.plot(x, aged.ruler_score_context, marker="o", ms=4, ls="none", color=OI["blue"], zorder=3)
    ax.set_xticks(x)
    ax.set_xticklabels([STATE_WORDS["Fibroblast_d0" if s == "Fibroblast" else s].replace("\n", " ") for s in order],
                       rotation=30, ha="right", rotation_mode="anchor")
    ax.set_xlim(-0.5, len(order) - 0.5)
    ax.set_ylim(-12, 5)
    ax.set_ylabel("Frozen ruler score (lower = scored younger)")
    ax.set_xlabel("Aged donor's cells")
    letter_in = -0.48
    panel_letter(ax, "A", x=letter_in / (ax.get_position().width * fig.get_figwidth()))
    for s, r in aged.iterrows():
        data.append(dict(panel="A", series="aged donor", label=s, x=float(order.index(s)), y=r.ruler_score_context,
                         n_cells=int(r.n_cells), note="Fibroblast = day-0 cells only" if s == "Fibroblast" else ""))
    data.append(dict(panel="A", series="young donor day 0 (dashed line)", label="Fibroblast", x=np.nan,
                     y=young0.ruler_score_context, n_cells=int(young0.n_cells), note="Fibroblast = day-0 cells only"))

    # B
    ax = fig.add_subplot(gsp[1])
    sc = ex[~ex.below_min_cells.astype(bool) & np.isfinite(ex.age_score.to_numpy(float))]
    ref = corr[(corr.distance == "nn_euclidean") & (corr.scope == "all_donor_cluster_timepoint_n_ge_20")].iloc[0]
    rho_chk = float(spearmanr(sc.nn_euclidean, sc.age_score).statistic)
    if len(sc) != int(ref.n) or abs(rho_chk - float(ref.rho)) > 1e-12:
        raise RuntimeError(f"t2_extrap rows (n={len(sc)}, rho={rho_chk}) do not match t2_correlation.csv")
    donor_color = {AGED_LINE: OI["orange"], YOUNG_LINE: OI["blue"]}
    donor_word = {AGED_LINE: "96-y donor", YOUNG_LINE: "22-y donor"}
    day_marker = {0: "o", 3: "s", 7: "^", 10: "D"}
    day_ms = {0: 4.2, 3: 4.0, 7: 4.4, 10: 3.8}
    for _, r in sc.iterrows():
        d = int(r.day)
        ax.plot(r.nn_euclidean, r.age_score, marker=day_marker[d], ms=day_ms[d], ls="none", mec="none",
                color=donor_color[r.cell_line], alpha=0.85, zorder=3)
        data.append(dict(panel="B", series=donor_word[r.cell_line], label=f"{r.cell_line} day {d} cluster {int(r.cluster)}",
                         x=r.nn_euclidean, y=r.age_score, n_cells=int(r.n_cells), note=f"day={d}"))
    ax.text(0.98, 0.98, f"ρ = {ref.rho:.3f}\n95% CI [{ref.rho_ci_lo:.3f}, {ref.rho_ci_hi:.3f}]\n"
            f"n = {int(ref.n)} pseudobulks".replace("-", "−"), transform=ax.transAxes, ha="right", va="top",
            fontsize=6.5)
    ax.set_ylim(-20, 10)
    ax.set_xlabel("Nearest-neighbour distance to GTEx training donors\n(PCA, k = 50)")
    ax.set_ylabel("Frozen-ruler score (a.u.)")
    panel_letter(ax, "B", x=letter_in / (ax.get_position().width * fig.get_figwidth()))
    donor_h = [Line2D([], [], marker="o", ls="none", ms=4.2, mec="none", color=donor_color[l], label=donor_word[l])
               for l in LINES]
    day_h = [Line2D([], [], marker=day_marker[d], ls="none", ms=day_ms[d], mec="none", color=OI["grey"],
                    label=f"day {d}") for d in day_marker]
    leg = ax.legend(handles=donor_h, loc="lower left", frameon=False, handletextpad=0.3, borderaxespad=0.3)
    ax.add_artist(leg)
    ax.legend(handles=day_h, loc="lower left", bbox_to_anchor=(0.36, 0.0), ncol=2, frameon=False,
              handletextpad=0.3, columnspacing=1.0, borderaxespad=0.3)
    data.append(dict(panel="B", series="annotation", label="spearman nn_euclidean", x=np.nan, y=float(ref.rho),
                     y_lo=float(ref.rho_ci_lo), y_hi=float(ref.rho_ci_hi), n_cells=np.nan,
                     note=f"n={int(ref.n)}; from {FIBRO3_CORR.relative_to(ROOT).as_posix()}"))
    save_fig(fig, "fig_N5_ruler_same_mistake")
    plt.close(fig)
    pd.DataFrame(data).to_csv(OUT / "fig_N5_ruler_same_mistake_data.csv", index=False)
    log("[C] wrote fig_N5_ruler_same_mistake.{png,pdf,_data.csv}")


# --------------------------------------------------------------------------- numbers file
def md_tbl(df, fmt=None):
    fmt = fmt or {}
    cols = list(df.columns)
    out = ["| " + " | ".join(cols) + " |", "|" + "|".join("---" for _ in cols) + "|"]
    for _, r in df.iterrows():
        cells = []
        for c in cols:
            v = r[c]
            if c in fmt and isinstance(v, (int, float, np.integer, np.floating)) and np.isfinite(v):
                cells.append(format(v, fmt[c]))
            elif isinstance(v, (float, np.floating)):
                cells.append("NA" if not np.isfinite(v) else f"{v:.4g}")
            else:
                cells.append(str(v))
        out.append("| " + " | ".join(cells) + " |")
    return "\n".join(out)


def write_numbers(log, integrity=None):
    ma = json.loads((OUT / "partA_meta.json").read_text(encoding="utf-8"))
    mb = json.loads((OUT / "partB_meta.json").read_text(encoding="utf-8"))
    st = pd.read_csv(OUT / "partA_md_by_state.csv")
    gd = pd.read_csv(OUT / "partA_gap_drop.csv")
    same = pd.read_csv(OUT / "partB_same_distances.csv")
    mix = pd.read_csv(OUT / "partB_mixture.csv")
    gt = pd.read_csv(OUT / "partB_gtex_distances.csv")
    pct = pd.read_csv(OUT / "partB_random_sets_repro.csv")
    cnt = pd.read_csv(OUT / "partB_random_sets_counts.csv")
    rnd = pd.read_csv(OUT / "partB_random_sets_per_set.csv")
    f4 = ".4f"

    L = ["# NUMBERS_NEWSTORY", "",
         "Numbers for the restructured paper. Produced by `python src/newstory.py all`. Writes only "
         "`results/newstory/` and this file. paper/ and paper_package/ were not read.", ""]

    # ---- A
    rp = ma["perdonor_repro"]
    L += ["## Part A. The MD score per state (GSE297234)", "",
          f"MD = `results/md3/genesets.json` \"MD\" ({ma['md_genes']} genes, all mapped). Seurat AddModuleScore "
          f"from `src/md2_score.py` (LogNormalize scale 10,000; nbin {ma['nbin']}; ctrl {ma['ctrl']}; seed {ma['seed']}). "
          "Labels: md2 Louvain clusters (`results/md2/louvain_obs_*.csv` + `t2_cluster_labels_*.csv`). "
          "Day-0 Fibroblast = label Fibroblast at day 0; other states pooled over days 0/3/7/10. Donors never pooled.", "",
          "**Joint scoring.** Both donors' Louvain count matrices were put on one gene axis (union by Ensembl id: "
          f"{ma['union_genes']:,} genes; {ma['shared_genes']:,} shared, {ma['aged_only']} aged-only, "
          f"{ma['young_only']} young-only, zero-filled where a donor lacks the gene, which leaves every cell's "
          "library size and LogNormalize values unchanged). Gene means were taken over all "
          f"{sum(ma['n_cells'].values()):,} cells of both donors; the 24 bins and the control genes "
          f"({ma['n_ctrl_joint']:,} unique) were drawn once from those means and used for every cell. Because "
          "the score of a cell depends only on its own row, this equals scoring one stacked matrix; it was run "
          "donor by donor for memory (3.8 GB free). MD genes absent from a donor's own Louvain matrix (fewer "
          "than 5 cells; zero in the joint matrix for that donor): "
          + "; ".join(f"{k}: {', '.join(v['md_absent_from_donor_matrix']) or 'none'}"
                      for k, v in ma["feature_check"].items()) + ".", "",
          "**How md4 was scored: per donor, not jointly.** " + ma["md4_scoring"] + " Rerunning that per-donor "
          f"scoring reproduces the stored vectors: max |diff| {rp[AGED_LINE]['max_abs_diff_vs_stored_md']:.1e} "
          f"(GM00731), {rp[YOUNG_LINE]['max_abs_diff_vs_stored_md']:.1e} (GM23815); bins identical "
          f"({rp[AGED_LINE]['bins_equal_stored']}, {rp[YOUNG_LINE]['bins_equal_stored']}); per-donor control sets "
          f"{rp[AGED_LINE]['n_ctrl']:,} and {rp[YOUNG_LINE]['n_ctrl']:,} genes.", "",
          "### A1. MD score per state", "",
          f"Mean over cells; 95% percentile bootstrap CI over cells (B = {ma['n_boot']}, one Generator({ma['boot_seed']}) "
          "used in order aged then young, states in table order). `per-donor` = md2 stored score "
          "(each donor scored alone). `md4` = n-weighted mean of `md_score_mean` over days in "
          "`results/md4/t1_state_timepoint_means.csv`.", "",
          "Source: `results/newstory/partA_md_by_state.csv`", ""]
    t = st[["donor", "cell_line", "state", "n_cells", "joint_mean", "joint_ci_lo", "joint_ci_hi",
            "perdonor_mean", "md4_nweighted_mean", "n_cells_md4"]].copy()
    L += [md_tbl(t, {c: f4 for c in ("joint_mean", "joint_ci_lo", "joint_ci_hi", "perdonor_mean",
                                     "md4_nweighted_mean")}), ""]
    L += ["### A2. Old–young gap, drop, drop ÷ gap", "",
          "gap = aged day-0 − young day-0; drop = aged day-0 − aged PartialReprog. CIs (joint only) from the same "
          "bootstrap draws, states resampled independently.", "",
          "Source: `results/newstory/partA_gap_drop.csv`", ""]
    g = gd.copy()
    for a, d in (("gap", 4), ("drop", 4), ("ratio", 2)):
        lo, hi = f"{a}_ci_lo", f"{a}_ci_hi"
        if lo in g:
            g[f"{a}_CI"] = [f"[{x:.{d}f}, {y:.{d}f}]" if np.isfinite(x) else "" for x, y in zip(g[lo], g[hi])]
    L += [md_tbl(g[["scoring", "aged_d0", "young_d0", "aged_PartialReprog", "gap", "gap_CI", "drop", "drop_CI",
                    "drop_over_gap", "ratio_CI"]],
                 {c: f4 for c in ("aged_d0", "young_d0", "aged_PartialReprog", "gap", "drop")} | {"drop_over_gap": ".2f"}),
          ""]

    # ---- B
    L += ["## Part B. Extra geometry numbers (point estimates)", "",
          "Recomputed with `src/genespace_run.py` functions (`map_gene_spaces`, `same_point_block`, "
          "`gtex_point_block`, `draw_baseline_sets`) on the same panels. Every recomputed quantity that is also "
          f"stored in `results/genespace/` matches it: max |diff| {mb['max_abs_diff_point']:.1e} "
          "(`results/newstory/partB_repro_checks.csv`). Distances are Euclidean in frozen-ruler z units within "
          "each gene space.", "",
          "### B1. Same-platform test: starting and final distance", "",
          "Forward: start = aged day-0 half A (O), target = young day-0 half A (Y), cells = aged PartialReprog. "
          "Reverse: start = Y, target = O, cells = young PartialReprog. Starting distance = |O − Y| (same for both "
          "tests). Pluripotency and NonReprog (aged, forward frame) are context rows.", "",
          "Source: `results/newstory/partB_same_distances.csv`", ""]
    L += [md_tbl(same[["space", "n_genes", "test", "n_cells", "dist_start", "dist_final", "final_over_start",
                       "progress", "perp", "delta"]],
                 {"dist_start": ".3f", "dist_final": ".3f", "final_over_start": f4, "progress": f4, "perp": f4,
                  "delta": ".3f"}), ""]
    L += ["Mixture control (aged half B + young half B cells, f = young fraction; start O, target Y):", "",
          "Source: `results/newstory/partB_mixture.csv`", "",
          md_tbl(mix[["space", "f", "n_aged", "n_young", "dist_start", "dist_final", "delta", "rel_delta",
                      "progress", "perp"]],
                 {"f": ".2f", "dist_start": ".3f", "dist_final": ".3f", "delta": ".3f", "rel_delta": f4,
                  "progress": f4, "perp": f4}), ""]
    L += ["### B2. GTEx test: distances to the GTEx centroids", "",
          f"GTEx young centroid = mean z of GTEx fibroblast donors aged 20–39 (n = {mb['anchor_n']['n_young']}); "
          f"old = 60–79 (n = {mb['anchor_n']['n_old']}). Day-0 pseudobulk = all day-0 Fibroblast cells of the "
          "donor, TMM with that donor's state rows (the GTEx-test panel, not the same-platform half A).", "",
          "Source: `results/newstory/partB_gtex_distances.csv`", ""]
    g0 = gt[gt.state == "PartialReprog"][["space", "n_genes", "cell_line", "dist_d0_to_gtex_young",
                                           "dist_d0_to_gtex_old", "gtex_young_to_old",
                                           "young_to_old_over_d0_to_young"]]
    L += [md_tbl(g0, {"dist_d0_to_gtex_young": ".3f", "dist_d0_to_gtex_old": ".3f", "gtex_young_to_old": ".3f",
                      "young_to_old_over_d0_to_young": f4}), "",
          "Change in distance from day 0, aged donor (positive = farther; Figure N1-C uses the MD rows):", "",
          md_tbl(gt[gt.cell_line == AGED_LINE][["space", "state", "n_cells", "delta_young", "rel_delta_young",
                                                "delta_old", "rel_delta_old"]],
                 {"delta_young": ".3f", "rel_delta_young": f4, "delta_old": ".3f", "rel_delta_old": f4}), ""]
    L += ["### B3. Random gene sets (Amendment 1: matched on 24 GTEx-mean bins × zero/nonzero; 200 per size)", "",
          "Per-set values were not stored, so they were recomputed from `results/genespace/size_baseline_sets.npz` "
          "(`A1_mu_zero__MD`, `A1_mu_zero__AGE`; redrawing with `draw_baseline_sets` gives identical sets). "
          "They reproduce the stored percentiles and positions in `results/genespace/size_baseline.csv`: "
          f"max |diff| {mb['max_abs_diff_percentiles']:.1e}.", "",
          "Per-set values: `results/newstory/partB_random_sets_per_set.csv`. Percentile check: "
          "`results/newstory/partB_random_sets_repro.csv`.", "",
          md_tbl(pct[["space_size", "test", "state", "stat", "real", "p2_5", "p50", "p97_5", "position",
                      "max_abs_diff"]],
                 {"real": f4, "p2_5": f4, "p50": f4, "p97_5": f4, "position": ".3f", "max_abs_diff": ".1e"}), "",
          "How many random sets end closer to the target (change in distance < 0):", "",
          "Source: `results/newstory/partB_random_sets_counts.csv`", "",
          md_tbl(cnt[["space_size", "n_sets", "fwd_closer", "rev_closer", "mix50_closer", "fwd_min_final_over_start",
                      "real_fwd_final_over_start"]],
                 {"fwd_min_final_over_start": f4, "real_fwd_final_over_start": f4}), ""]
    rr = rnd.groupby("space_size")[["fwd_progress", "fwd_perp", "rev_progress", "rev_perp", "mix50_delta",
                                    "mix50_rel_delta"]].agg(["min", "median", "max"])
    rr.columns = [f"{a}_{b}" for a, b in rr.columns]
    L += ["Range of per-set values (min / median / max):", "",
          md_tbl(rr.reset_index(), {c: f4 for c in rr.columns}), ""]

    L += ["## Figures", "",
          "All in `results/newstory/` (PNG 300 dpi, 7 in wide; PDF; plotting data CSV).", "",
          "- Figure N1, \"MD's own genes\": `fig_N1_md_own_genes.png`, `fig_N1_md_own_genes.pdf`, "
          "`fig_N1_md_own_genes_data.csv`",
          "  - A: MD score (joint scoring) for aged day-0, PartialReprog, EarlyPluripotency, Pluripotency; "
          "dashed line = young day-0.",
          "  - B: same-platform plane in the MD genes (x = progress, y = sideways, axis lengths).",
          "  - C: change in distance to GTEx young / old in the MD genes, aged donor.",
          "- Figure N2, \"Any genes\": `fig_N2_any_genes.png`, `fig_N2_any_genes.pdf`, `fig_N2_any_genes_data.csv`",
          "  - A / B: forward / reverse plane; real MD, age, all genes; 200 random sets per size.",
          "  - C: mixture control, change in distance ÷ starting distance vs f, three spaces.", "",
          "Other outputs: `partA_cell_scores.npz` (per-cell joint and per-donor MD scores, labels, shared control "
          "columns), `partA_meta.json`, `partB_meta.json`, `run_log.txt`.", ""]
    if integrity is not None:
        L += ["## Integrity", "",
              f"Files outside `results/newstory/` (results/, src/, data/processed/md2, data/processed/md4, root *.md) "
              f"compared by size and mtime before and after the run: changed = {integrity or 'none'}.", ""]
    NUMBERS_PATH.write_text("\n".join(L) + "\n", encoding="utf-8")
    log(f"[numbers] wrote {NUMBERS_PATH.name}")


def main():
    step = sys.argv[1] if len(sys.argv) > 1 else "all"
    OUT.mkdir(parents=True, exist_ok=True)
    log = Log(OUT / "run_log.txt")
    before = snapshot()
    if step in ("a", "all"):
        part_a(log)
    if step in ("b", "all"):
        part_b(log)
    if step in ("figs", "all"):
        figure_n1(log)
        figure_n2(log)
        figure_n5(log)
    changed = snapshot_diff(before, snapshot())
    log(f"[integrity] files outside results/newstory changed during this step: {changed or 'none'}")
    if step in ("numbers", "all", "figs"):
        write_numbers(log, integrity=changed)
    log.close()


if __name__ == "__main__":
    main()
