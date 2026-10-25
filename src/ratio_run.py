"""RATIO: is the typical perturbation short enough, against its own experiment's
activated-to-resting gap, that arrival could register at all?

Pre-registration: results/ratio/PREREG.flag (refuses to run without it).

Modes
  python -I src/ratio_run.py reproduce    R1/R1b (SOUTH6), R2 (GSE278572 target count)
  python -I src/ratio_run.py gse278572    stream the deposit, then every sub-screen
  python -I src/ratio_run.py verify       recompute every reading from tracked CSVs only

Raw inputs (git-ignored, data/raw/ratio/):
  gse278572/  https://ftp.ncbi.nlm.nih.gov/geo/series/GSE278nnn/GSE278572/suppl/
              GSE278572_{matrix.mtx.gz, features.tsv.gz, barcodes.tsv.gz,
              protospacer_calls_per_cell.csv.gz}

Statistic, per perturbation of one sub-screen (flag D1-D4):
  ratio = ||x_pert - x_act_ctl|| / ||x_act_ctl - x_rest_ctl||
in log2-CPM (TMM, prior.count 2) over every Gene Expression feature (U4).
Floor, positive control and readings: flag D6, D7, D9.
"""
from __future__ import annotations

import gzip
import json
import multiprocessing as mp
import sys
import time
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import sparse, stats as sps

sys.path.insert(0, str(Path(__file__).resolve().parent))

from toward_run import ANCHORS_PATH, tmm_logcpm_quiet  # noqa: E402
from seng_run import tmm_cache_from_panel, z_from_counts  # noqa: E402
from fibro2_common import load_frozen_ruler  # noqa: E402

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "results" / "ratio"
FLAG = OUT / "PREREG.flag"
RAW = ROOT / "data" / "raw" / "ratio"
SOUTH6 = ROOT / "results" / "south6"

MIN_CELLS = 30          # D5
N_FLOOR = 1000          # D6 / U3
P_FLOOR = 95.0
EXACT_N_MAX = 200       # D6: exact sizes up to this many distinct n
GRID_N = 60
POOL_CAP = 60000        # D6, GSE314342 only
SEED_FLOOR = 20261025
SEED_PC = 20261026
PC_F = (0.0, 0.1, 0.25, 0.5, 1.0)
N_PC = 200
CHANCE = 0.05
HTO_MIN, HTO_FOLD = 10, 3.0    # D3
N_WORKERS = 8

GSE278572_HTO = {"HTO_1": "Resting_Treg", "HTO_2": "Stimulated_Treg",
                 "HTO_3": "Resting_Teff", "HTO_4": "Stimulated_Teff"}
GSE278572_SUBS = {"Treg": ("Stimulated_Treg", "Resting_Treg"),
                  "Teff": ("Stimulated_Teff", "Resting_Teff")}
GSE278572_NT = "Non-Targeting"
EXPECTED_TARGETS = {"GSE278572": 28}

READ_BLOCKED = ("RATIO-BLOCKED: in {ds}, {n} of {N} perturbations have ratio >= 2 in {sp}; "
                "arrival could not register for them whatever their direction.")
READ_OPEN = ("RATIO-OPEN: in {ds}, {n} of {N} perturbations have ratio < 2 in {sp}; "
             "arrival was geometrically attainable for them.")
READ_UNMEAS = ("RATIO-UNMEASURABLE: in {ds}, the perturbation ratio distribution is not "
               "distinguishable from the control floor in {sp}; this screen cannot measure "
               "displacement length at this depth.")
SPACE = "full measured-gene space"


class Stop(Exception):
    def __init__(self, key, msg):
        super().__init__(f"STOP {key}: {msg}")
        self.key = key


LOG_FH = None


def log(msg):
    line = f"[{time.strftime('%H:%M:%S')}] {msg}"
    print(line, flush=True)
    if LOG_FH is not None:
        LOG_FH.write(line + "\n")
        LOG_FH.flush()


def dump(obj, path):
    Path(path).write_text(json.dumps(obj, indent=1, default=lambda o: o.item()
                                     if hasattr(o, "item") else str(o)), encoding="utf-8")


def require_flag():
    if not FLAG.exists():
        raise SystemExit(f"{FLAG} missing; refusing to compute.")
    log(f"[prereg] {FLAG.relative_to(ROOT)} bytes={FLAG.stat().st_size}")


# ------------------------------------------------------------------ normalisation

class Space:
    """One sub-screen's panel normalisation (flag D4)."""

    def __init__(self, C):
        C = np.asarray(C, np.float64)
        self.logcpm, nf = tmm_logcpm_quiet(C)
        self.cache = tmm_cache_from_panel(C, nf)
        g = C.shape[1]
        self.mu0, self.sd1, self.nomiss = np.zeros(g), np.ones(g), np.zeros(g, bool)

    def lc(self, counts):
        return z_from_counts(np.asarray(counts, np.float64).ravel(), self.cache,
                             self.mu0, self.sd1, self.nomiss)[1]


_W = {}


def _winit(cache, gap_norm, act_path, rest_path):
    # Pools travel by file: pickling ~500 MB of initargs to spawned workers fails on Windows
    # (OSError 22, first GSE278572 attempt).
    act_pool, rest_pool = sparse.load_npz(act_path).tocsr(), sparse.load_npz(rest_path).tocsr()
    g = act_pool.shape[1]
    _W.update(cache=cache, mu0=np.zeros(g), sd1=np.ones(g), nomiss=np.zeros(g, bool),
              gap_norm=gap_norm, A=act_pool, R=rest_pool,
              A_total=np.asarray(act_pool.sum(0), np.float64).ravel())


def _wratio(spec):
    """spec = (activated-control rows drawn, resting-control rows drawn).
    draw = their sum; reference = activated pool minus the activated rows drawn (D6, D7)."""
    ia, ir = spec
    g = _W["A_total"].size
    a = np.asarray(_W["A"][ia].sum(0), np.float64).ravel() if ia.size else np.zeros(g)
    r = np.asarray(_W["R"][ir].sum(0), np.float64).ravel() if ir.size else np.zeros(g)
    la = z_from_counts(a + r, _W["cache"], _W["mu0"], _W["sd1"], _W["nomiss"])[1]
    lb = z_from_counts(_W["A_total"] - a, _W["cache"], _W["mu0"], _W["sd1"], _W["nomiss"])[1]
    return float(np.linalg.norm(la - lb) / _W["gap_norm"])


def ratios_parallel(space, gap_norm, act_pool, rest_pool, specs):
    """Ratio for each (ia, ir) draw spec, order preserved."""
    out = np.empty(len(specs))
    act_path, rest_path = RAW / "_pool_act.npz", RAW / "_pool_rest.npz"
    sparse.save_npz(act_path, act_pool, compressed=False)
    sparse.save_npz(rest_path, rest_pool, compressed=False)
    with mp.Pool(N_WORKERS, initializer=_winit,
                 initargs=(space.cache, gap_norm, str(act_path), str(rest_path))) as pool:
        for i, r in enumerate(pool.imap(_wratio, specs, chunksize=32)):
            out[i] = r
    act_path.unlink()
    rest_path.unlink()
    return out


# ------------------------------------------------------------------ reproduce

def reproduce():
    rows = []

    def add(check, expected, observed, match, gating):
        rows.append(dict(check=check, expected=expected, observed=observed,
                         match=bool(match), gating=gating))
        log(f"[reproduce] {check}: expected={expected} observed={observed} match={match}")

    s1 = np.load(SOUTH6 / "stage1.npz", allow_pickle=True)
    d_f = np.asarray(s1["d_f"], np.float64)
    names = [str(x) for x in s1["keep"]]
    anc = np.load(ANCHORS_PATH, allow_pickle=True)
    gap = np.asarray(anc["c_old"], float).ravel() - np.asarray(anc["c_young"], float).ravel()
    if gap.size != d_f.shape[1]:
        raise Stop("reproduction_failed", f"anchor length {gap.size} != d_f width {d_f.shape[1]}")
    ratio = np.linalg.norm(d_f, axis=1) / np.linalg.norm(gap)
    st = pd.read_csv(SOUTH6 / "stage3_full_table.csv",
                     usecols=["factor", "space", "pair", "d_over_v", "v_norm"])
    ref = st[(st.space == "S1") & (st.pair == "O1_Y1")].set_index("factor")
    add("R1 SOUTH6 S1 O1_Y1 factor count", len(ref), len(names), len(ref) == len(names) == 1836, True)
    common = [n for n in names if n in ref.index]
    add("R1 factors matched by name", len(names), len(common), len(common) == len(names), True)
    obs = pd.Series(ratio, index=names).loc[common].to_numpy()
    exp = ref.loc[common, "d_over_v"].to_numpy(float)
    rel = np.abs(obs - exp) / np.abs(exp)
    add("R1 max relative error, per-factor ratio", "<=1e-9", f"{rel.max():.3e}", rel.max() <= 1e-9, True)
    gn_ref = float(ref["v_norm"].iloc[0])
    add("R1 ||gap|| from frozen anchors vs stored v_norm", f"{gn_ref:.12f}",
        f"{np.linalg.norm(gap):.12f}", abs(np.linalg.norm(gap) - gn_ref) / gn_ref <= 1e-9, True)
    ov = pd.read_csv(SOUTH6 / "stage4_overall.csv")
    med_ref = float(ov[(ov.space == "S1") & (ov.pair == "O1_Y1")].ratio_median.iloc[0])
    med = float(np.median(obs))
    add("R1 median ratio vs stage4_overall", f"{med_ref:.12f}", f"{med:.12f}",
        abs(med - med_ref) / med_ref <= 1e-9, True)

    full = pd.read_csv(SOUTH6 / "stage3_full_table.csv", usecols=["space", "pair", "d_over_v"])
    for (sp, pr), g in full[full.pair == "O1_Y1"].groupby(["space", "pair"], sort=False):
        r = g.d_over_v.to_numpy(float)
        add(f"R1b premise {sp} {pr}: min ratio; n>=2 of N", "descriptive",
            f"min={r.min():.4f}; {int((r >= 2).sum())} of {r.size}", True, False)

    feats = pd.read_csv(RAW / "gse278572" / "GSE278572_features.tsv.gz", sep="\t", header=None,
                        names=["id", "name", "type"])
    g = feats[feats.type == "CRISPR Guide Capture"].name.str.replace(r"_\d+$", "", regex=True)
    targets = sorted(set(g) - {GSE278572_NT})
    add("R2 GSE278572 distinct targeting genes", 28, len(targets), len(targets) == 28, True)

    df = pd.DataFrame(rows)
    df.to_csv(OUT / "reproduction.csv", index=False)
    bad = df[df.gating & ~df.match]
    if len(bad):
        raise Stop("reproduction_failed", "; ".join(f"{r.check}: expected {r.expected}, "
                                                    f"observed {r.observed}" for r in bad.itertuples()))
    log(f"[reproduce] {int(df.gating.sum())} gating checks, 0 mismatches")
    return df


# ------------------------------------------------------------------ GSE278572 I/O

def mtx_chunks(path, chunksize=40_000_000):
    fh = gzip.open(path, "rt")
    while True:
        pos_line = fh.readline()
        if not pos_line.startswith("%"):
            break
    n_feat, n_bc, nnz = (int(x) for x in pos_line.split())
    reader = pd.read_csv(fh, sep=" ", header=None, names=["f", "c", "v"],
                         dtype={"f": np.int32, "c": np.int32, "v": np.int32},
                         chunksize=chunksize, engine="c")
    return (n_feat, n_bc, nnz), reader, fh


def gse278572_load():
    d = RAW / "gse278572"
    feats = pd.read_csv(d / "GSE278572_features.tsv.gz", sep="\t", header=None,
                        names=["id", "name", "type"])
    bcs = pd.read_csv(d / "GSE278572_barcodes.tsv.gz", header=None, names=["bc"]).bc.to_numpy()
    calls = pd.read_csv(d / "GSE278572_protospacer_calls_per_cell.csv.gz")
    log(f"[278572] features={len(feats)} barcodes={len(bcs)} call rows={len(calls)}")

    ge = np.flatnonzero(feats.type.to_numpy() == "Gene Expression")
    hto_names = list(GSE278572_HTO)
    hto_idx = [int(np.flatnonzero(feats.id.to_numpy() == h)[0]) for h in hto_names]
    feat_to_ge = np.full(len(feats), -1, np.int64)
    feat_to_ge[ge] = np.arange(ge.size)
    feat_to_hto = np.full(len(feats), -1, np.int64)
    feat_to_hto[hto_idx] = np.arange(4)

    # guide calls (D3): exactly one called guide
    bc_index = pd.Series(np.arange(len(bcs)), index=bcs)
    calls = calls[calls.cell_barcode.isin(bc_index.index)]
    one = calls[calls.num_features == 1]
    target_of_guide = lambda s: s.replace("_CRISPRi", "").rsplit("_", 1)[0]  # noqa: E731
    cell_target = np.full(len(bcs), "", dtype=object)
    cell_target[bc_index[one.cell_barcode].to_numpy()] = [target_of_guide(s) for s in one.feature_call]
    n_calls = np.zeros(len(bcs), np.int64)
    n_calls[bc_index[calls.cell_barcode].to_numpy()] = calls.num_features.to_numpy()
    guide_stats = dict(barcodes=len(bcs), zero_calls=int((n_calls == 0).sum()),
                       one_call=int((n_calls == 1).sum()), multi_calls=int((n_calls >= 2).sum()))
    log(f"[278572] guide calls {guide_stats}")

    # pass 1: HTO counts
    t0 = time.time()
    hto = np.zeros((len(bcs), 4), np.int64)
    (nf, nb, nnz), reader, fh = mtx_chunks(d / "GSE278572_matrix.mtx.gz")
    if nf != len(feats) or nb != len(bcs):
        raise Stop("flag_ambiguous", f"mtx dims {nf}x{nb} do not match features/barcodes")
    seen = 0
    for ch in reader:
        f = ch.f.to_numpy() - 1
        h = feat_to_hto[f]
        m = h >= 0
        np.add.at(hto, (ch.c.to_numpy()[m] - 1, h[m]), ch.v.to_numpy()[m])
        seen += len(ch)
    fh.close()
    log(f"[278572] pass 1 (HTO) entries={seen} of nnz={nnz} in {time.time() - t0:.0f}s")
    if seen != nnz:
        raise Stop("flag_ambiguous", f"mtx entries read {seen} != header nnz {nnz}")

    # D3 HTO rule
    srt = np.sort(hto, axis=1)
    top, second = srt[:, -1], srt[:, -2]
    arg = hto.argmax(1)
    ok = (top >= HTO_MIN) & ((second == 0) | (top >= HTO_FOLD * second))
    cond = np.where(ok, np.array([GSE278572_HTO[h] for h in hto_names], dtype=object)[arg], "")
    hto_stats = dict(assigned=int(ok.sum()), unassigned=int((~ok).sum()),
                     **{c: int((cond == c).sum()) for c in GSE278572_HTO.values()})
    log(f"[278572] HTO {hto_stats}")

    # groups: (condition, target) for single-guide, HTO-assigned cells
    use = ok & (cell_target != "")
    keys = pd.Series([f"{c}|{t}" for c, t in zip(cond[use], cell_target[use])])
    gnames = sorted(set(keys))
    gid_of = {k: i for i, k in enumerate(gnames)}
    cell_gid = np.full(len(bcs), -1, np.int64)
    cell_gid[np.flatnonzero(use)] = keys.map(gid_of).to_numpy()
    # control cells kept at cell level, all four conditions
    ctl_cells = np.flatnonzero(use & (cell_target == GSE278572_NT))
    ctl_row = np.full(len(bcs), -1, np.int64)
    ctl_row[ctl_cells] = np.arange(ctl_cells.size)

    # pass 2: group sums and control-cell rows
    t0 = time.time()
    G = ge.size
    sums = np.zeros(len(gnames) * G, np.float64)
    rr, cc, vv = [], [], []
    (_, _, _), reader, fh = mtx_chunks(d / "GSE278572_matrix.mtx.gz")
    for ch in reader:
        f = ch.f.to_numpy() - 1
        c = ch.c.to_numpy() - 1
        v = ch.v.to_numpy()
        gcol = feat_to_ge[f]
        gid = cell_gid[c]
        m = (gcol >= 0) & (gid >= 0)
        sums += np.bincount(gid[m] * G + gcol[m], weights=v[m], minlength=sums.size)
        cr = ctl_row[c]
        m2 = (gcol >= 0) & (cr >= 0)
        rr.append(cr[m2]); cc.append(gcol[m2]); vv.append(v[m2])
    fh.close()
    sums = sums.reshape(len(gnames), G)
    ctl = sparse.csr_matrix((np.concatenate(vv).astype(np.float64),
                             (np.concatenate(rr), np.concatenate(cc))),
                            shape=(ctl_cells.size, G))
    log(f"[278572] pass 2 (sums) groups={len(gnames)} control cells={ctl_cells.size} "
        f"in {time.time() - t0:.0f}s")

    n_cells = pd.Series(keys).value_counts()
    counts = pd.DataFrame([dict(condition=k.split("|")[0], target=k.split("|")[1],
                                cells=int(n_cells[k])) for k in gnames])
    return dict(feats=feats, ge=ge, gnames=gnames, sums=sums, counts=counts,
                ctl=ctl, ctl_cond=cond[ctl_cells], guide_stats=guide_stats,
                hto_stats=hto_stats)


# ------------------------------------------------------------------ one sub-screen

def floor_sizes(ns):
    ns = np.asarray(sorted(set(int(n) for n in ns)))
    if ns.size <= EXACT_N_MAX:
        return ns, "exact"
    grid = np.unique(np.round(np.geomspace(MIN_CELLS, ns.max(), GRID_N)).astype(int))
    return grid, "grid"


def size_used(n, sizes, mode):
    if mode == "exact":
        return int(n)
    return int(sizes[sizes <= n].max())


def analyse_subscreen(ds, sub, pert_counts, pert_sums, act_pool, rest_pool, act_sum, rest_sum):
    """pert_counts: Series target->cells (activated arm). Pools are CSR cell x gene."""
    tag = f"{ds} {sub}"
    if act_pool.shape[0] == 0 or rest_pool.shape[0] == 0:
        raise Stop("control_pool_empty", f"{tag}: activated control cells={act_pool.shape[0]}, "
                                         f"resting control cells={rest_pool.shape[0]}")
    elig = pert_counts[pert_counts >= MIN_CELLS]
    excluded = pert_counts[pert_counts < MIN_CELLS]
    log(f"[{tag}] perturbations={len(pert_counts)} eligible={len(elig)} "
        f"excluded(<{MIN_CELLS})={len(excluded)} act_ctl={act_pool.shape[0]} "
        f"rest_ctl={rest_pool.shape[0]}")
    names = list(elig.index)
    C = np.vstack([pert_sums[names], act_sum[None, :], rest_sum[None, :]])
    sp = Space(C)
    x = sp.logcpm
    xa, xr = x[-2], x[-1]
    gap_norm = float(np.linalg.norm(xa - xr))
    selfcheck = float(np.max(np.abs(sp.lc(act_sum) - xa)))
    log(f"[{tag}] ||gap||={gap_norm:.4f}; cache-vs-panel logCPM max abs diff={selfcheck:.2e}")
    if selfcheck > 1e-8:
        raise Stop("flag_ambiguous", f"{tag}: cached normalisation differs from panel by {selfcheck}")
    ratio = np.linalg.norm(x[:-2] - xa, axis=1) / gap_norm

    # floor (D6)
    P = act_pool.shape[0]
    med_n = int(np.floor(np.median(elig.to_numpy()))) if len(elig) else MIN_CELLS
    sizes, mode = floor_sizes(list(elig.to_numpy()) + [MIN_CELLS, med_n])
    sizes = sizes[sizes <= P // 2]
    rng = np.random.default_rng(SEED_FLOOR)
    empty = np.zeros(0, np.int64)
    idx_sets = [(rng.choice(P, int(n), replace=False), empty) for n in sizes for _ in range(N_FLOOR)]
    t0 = time.time()
    fr = ratios_parallel(sp, gap_norm, act_pool, rest_pool, idx_sets)
    log(f"[{tag}] floor: {len(sizes)} sizes ({mode}) x {N_FLOOR} draws in {time.time() - t0:.0f}s")
    fr = fr.reshape(len(sizes), N_FLOOR)
    p95 = {int(n): float(np.percentile(fr[i], P_FLOOR)) for i, n in enumerate(sizes)}
    floor_tab = pd.DataFrame(dict(n=sizes, p50=np.median(fr, 1), p95=[p95[int(n)] for n in sizes],
                                  frac_ge2=(fr >= 2).mean(1), max=fr.max(1)))

    rows = []
    for t, n, r in zip(names, elig.to_numpy(), ratio):
        avail = n <= P // 2
        nu = size_used(n, sizes, mode) if avail else None
        f95 = p95[nu] if avail else np.nan
        rows.append(dict(dataset=ds, sub_screen=sub, target=t, cells=int(n), ratio=float(r),
                         ratio_ge2=bool(r >= 2), floor_n_used=nu, floor_p95=f95,
                         floor_unavailable=not avail,
                         exceeds=bool(avail and r > f95), within_floor=bool(avail and r <= f95)))
    pt = pd.DataFrame(rows)

    # positive control (D7)
    rng2 = np.random.default_rng(SEED_PC)
    R = rest_pool.shape[0]
    pc_sizes = sorted(set(n for n in (MIN_CELLS, med_n) if n <= P // 2))
    pc_specs = []
    for n in pc_sizes:
        nu = size_used(n, sizes, mode)
        for f in PC_F:
            m = int(round(f * nu))
            if m > R:
                raise Stop("flag_ambiguous", f"{tag}: positive control needs {m} resting cells, pool has {R}")
            for _ in range(N_PC):
                ia = rng2.choice(P, nu - m, replace=False)
                ir = rng2.choice(R, m, replace=False)
                pc_specs.append((n, nu, f, ia, ir))

    pcr = ratios_parallel(sp, gap_norm, act_pool, rest_pool, [(s[3], s[4]) for s in pc_specs])
    pcd = pd.DataFrame([dict(n=s[0], n_used=s[1], f=s[2]) for s in pc_specs])
    pcd["ratio"] = pcr
    pc_tab = (pcd.groupby(["n", "n_used", "f"])
              .agg(median_ratio=("ratio", "median"),
                   power=("ratio", lambda v: float(np.mean(v > p95[int(pcd.loc[v.index[0], 'n_used'])]))),
                   draws=("ratio", "size")).reset_index())

    summ = summarise(ds, sub, pt, len(pert_counts), len(excluded), gap_norm, floor_tab, med_n,
                     act_pool.shape[0], rest_pool.shape[0], mode)
    stem = f"{ds.lower()}_{sub.lower()}"
    pt.to_csv(OUT / f"{stem}_perturbations.csv", index=False)
    floor_tab.to_csv(OUT / f"{stem}_floor.csv", index=False)
    pc_tab.to_csv(OUT / f"{stem}_posctrl.csv", index=False)
    pd.DataFrame(dict(target=list(excluded.index), cells=excluded.to_numpy())).to_csv(
        OUT / f"{stem}_excluded.csv", index=False)
    np.savez_compressed(OUT / f"{stem}_floor_draws.npz", sizes=sizes, ratios=fr)
    return summ


def summarise(ds, sub, pt, n_all, n_excl, gap_norm, floor_tab, med_n, n_act, n_rest, mode):
    r = pt.ratio.to_numpy()
    fa = pt[~pt.floor_unavailable]
    k, Nf = int(fa.exceeds.sum()), len(fa)
    p = float(sps.binomtest(k, Nf, CHANCE, alternative="greater").pvalue) if Nf else np.nan
    ex = fa[fa.exceeds]
    fl_med = floor_tab[floor_tab.n == med_n] if (floor_tab.n == med_n).any() else floor_tab.iloc[[0]]
    rho = float(sps.spearmanr(pt.ratio, pt.cells).statistic) if len(pt) > 2 else np.nan
    return dict(
        dataset=ds, sub_screen=sub, space=SPACE, perturbations=n_all, excluded_lt30=n_excl,
        eligible=len(pt), act_ctl_cells=n_act, rest_ctl_cells=n_rest, gap_norm=gap_norm,
        ratio_median=float(np.median(r)) if r.size else np.nan,
        ratio_q1=float(np.percentile(r, 25)) if r.size else np.nan,
        ratio_q3=float(np.percentile(r, 75)) if r.size else np.nan,
        ratio_min=float(r.min()) if r.size else np.nan, ratio_max=float(r.max()) if r.size else np.nan,
        frac_ge2=float((r >= 2).mean()) if r.size else np.nan,
        n_lt2_all=int((r < 2).sum()), n_lt2_within_floor=int(((r < 2) & pt.within_floor).sum()),
        floor_mode=mode, floor_p95_median_over_perts=float(fa.floor_p95.median()) if Nf else np.nan,
        floor_median_n=int(fl_med.n.iloc[0]), floor_p50_at_median_n=float(fl_med.p50.iloc[0]),
        floor_p95_at_median_n=float(fl_med.p95.iloc[0]),
        floor_frac_ge2_at_median_n=float(fl_med.frac_ge2.iloc[0]),
        spearman_ratio_vs_cells=rho,
        exceed_k=k, floor_available_N=Nf, chance_expected=CHANCE * Nf, binom_p=p,
        within_floor=int(fa.within_floor.sum()), floor_unavailable=int(pt.floor_unavailable.sum()),
        exceed_ge2=int((ex.ratio >= 2).sum()), exceed_lt2=int((ex.ratio < 2).sum()),
    )


def readings_of(s):
    """Flag D9 / U2. Returns (conclusions, descriptive)."""
    ds = f"{s['dataset']} {s['sub_screen']}"
    concl, desc = [], []
    tail = f" [within_floor = {s['within_floor']}" + (
        f"; floor_unavailable = {s['floor_unavailable']}" if s["floor_unavailable"] else "") + "]"
    if s["eligible"] == 0:
        desc.append(f"No reading fires in {ds}: zero eligible perturbations.")
        return concl, desc
    unmeas = not (s["binom_p"] < 0.05)
    N = s["exceed_ge2"] + s["exceed_lt2"]
    b = READ_BLOCKED.format(ds=ds, n=s["exceed_ge2"], N=N, sp=SPACE)
    o = READ_OPEN.format(ds=ds, n=s["exceed_lt2"], N=N, sp=SPACE)
    if unmeas:
        concl.append(READ_UNMEAS.format(ds=ds, sp=SPACE))
        desc.append(f"(descriptive, superseded by UNMEASURABLE) BLOCKED count {s['exceed_ge2']} of {N}; "
                    f"OPEN count {s['exceed_lt2']} of {N}" + tail)
    else:
        if s["exceed_ge2"] >= 1:
            concl.append(b + tail)
        if s["exceed_lt2"] >= 1:
            concl.append(o + tail)
        if not concl:
            desc.append(f"Neither BLOCKED nor OPEN fires in {ds}: no perturbation exceeds its floor" + tail)
    desc.append(f"(descriptive) {s['n_lt2_all']} of {s['eligible']} eligible perturbations have ratio < 2 "
                f"and therefore could in principle have registered an arrival; "
                f"{s['n_lt2_within_floor']} of those are within_floor.")
    return concl, desc


# ------------------------------------------------------------------ gene overlap

def gene_overlap(ds, feats_ge):
    fr = load_frozen_ruler()
    ens = set(str(e).split(".")[0] for e in fr["ensembl"])
    sym = set(str(s) for s in fr["symbol"])
    ids = feats_ge.id.astype(str).str.split(".").str[0]
    by_ens = ids.isin(ens)
    by_sym = ~by_ens & feats_ge.name.astype(str).isin(sym)
    return dict(dataset=ds, features=len(feats_ge), in_panel_by_ensembl=int(by_ens.sum()),
                in_panel_by_symbol_fallback=int(by_sym.sum()),
                in_panel_total=int(by_ens.sum() + by_sym.sum()), panel_size=int(len(fr["ensembl"])))


# ------------------------------------------------------------------ modes

def gse278572_cached():
    """Two passes over a 3.6 GB mtx are slow; their output is cached beside the raw data."""
    cp = RAW / "gse278572" / "_streamed.npz"
    if cp.exists():
        z = np.load(cp, allow_pickle=True)
        L = {k: z[k] for k in ("ge", "sums", "ctl_cond")}
        L["gnames"] = list(z["gnames"])
        L["ctl"] = sparse.csr_matrix((z["ctl_data"], z["ctl_indices"], z["ctl_indptr"]),
                                     shape=tuple(z["ctl_shape"]))
        L["counts"] = pd.DataFrame(json.loads(str(z["counts"])))
        L["guide_stats"], L["hto_stats"] = json.loads(str(z["guide_stats"])), json.loads(str(z["hto_stats"]))
        L["feats"] = pd.read_csv(RAW / "gse278572" / "GSE278572_features.tsv.gz", sep="\t",
                                 header=None, names=["id", "name", "type"])
        log(f"[278572] read cached stream output {cp.relative_to(ROOT)}")
        return L
    L = gse278572_load()
    np.savez(cp, ge=L["ge"], sums=L["sums"], ctl_cond=L["ctl_cond"], gnames=np.array(L["gnames"]),
             ctl_data=L["ctl"].data, ctl_indices=L["ctl"].indices, ctl_indptr=L["ctl"].indptr,
             ctl_shape=np.array(L["ctl"].shape), counts=json.dumps(L["counts"].to_dict("list")),
             guide_stats=json.dumps(L["guide_stats"]), hto_stats=json.dumps(L["hto_stats"]))
    return L


def run_gse278572():
    ds = "GSE278572"
    L = gse278572_cached()
    counts = L["counts"]
    piv = counts.pivot_table(index="target", columns="condition", values="cells", fill_value=0)
    piv.to_csv(OUT / "gse278572_cells_per_target.csv")
    tg = piv.drop(index=GSE278572_NT, errors="ignore")
    r3 = {c: int((tg[c] > 150).sum()) for c in tg.columns}
    dump(dict(guide_calls=L["guide_stats"], hto=L["hto_stats"],
              hto_rule=f"top >= {HTO_MIN} and >= {HTO_FOLD}x second",
              R3_targets_over_150_cells_per_condition=r3, n_targets=int(len(tg))),
         OUT / "gse278572_assignment.json")
    log(f"[278572] R3 (descriptive) targets with >150 cells per condition: {r3} of {len(tg)}")
    ov = gene_overlap(ds, L["feats"].iloc[L["ge"]])
    log(f"[278572] gene overlap {ov}")

    gidx = {k: i for i, k in enumerate(L["gnames"])}
    summaries = []
    for sub, (act, rest) in GSE278572_SUBS.items():
        tgt = counts[(counts.condition == act) & (counts.target != GSE278572_NT)]
        pert_counts = pd.Series(tgt.cells.to_numpy(), index=tgt.target.to_numpy())
        pert_sums = pd.Series({t: L["sums"][gidx[f"{act}|{t}"]] for t in pert_counts.index})
        pert_sums = np.vstack(pert_sums.to_numpy()) if len(pert_sums) else np.zeros((0, L["ge"].size))
        ps = {t: pert_sums[i] for i, t in enumerate(pert_counts.index)}
        act_pool = L["ctl"][np.flatnonzero(L["ctl_cond"] == act)]
        rest_pool = L["ctl"][np.flatnonzero(L["ctl_cond"] == rest)]
        act_sum = L["sums"][gidx[f"{act}|{GSE278572_NT}"]] if f"{act}|{GSE278572_NT}" in gidx else np.zeros(L["ge"].size)
        rest_sum = L["sums"][gidx[f"{rest}|{GSE278572_NT}"]] if f"{rest}|{GSE278572_NT}" in gidx else np.zeros(L["ge"].size)
        chk = max(np.abs(np.asarray(act_pool.sum(0)).ravel() - act_sum).max(),
                  np.abs(np.asarray(rest_pool.sum(0)).ravel() - rest_sum).max()) if act_pool.shape[0] and rest_pool.shape[0] else 0.0
        if chk > 0:
            raise Stop("flag_ambiguous", f"{sub}: control pool rows do not sum to the control pseudobulk ({chk})")
        summaries.append(analyse_subscreen(ds, sub, pert_counts, _Rows(ps), act_pool, rest_pool,
                                           act_sum, rest_sum))
    write_summary(summaries, [ov])
    return summaries


class _Rows:
    def __init__(self, d):
        self.d = d

    def __getitem__(self, names):
        return np.vstack([self.d[n] for n in names]) if names else np.zeros((0, 0))


def write_summary(summaries, overlaps):
    s_path, o_path = OUT / "summary.csv", OUT / "gene_overlap.csv"
    old = pd.read_csv(s_path) if s_path.exists() else pd.DataFrame()
    new = pd.DataFrame(summaries)
    if len(old):
        old = old[~old.dataset.isin(new.dataset.unique())]
    pd.concat([old, new], ignore_index=True).to_csv(s_path, index=False)
    oo = pd.read_csv(o_path) if o_path.exists() else pd.DataFrame()
    on = pd.DataFrame(overlaps)
    if len(oo):
        oo = oo[~oo.dataset.isin(on.dataset)]
    pd.concat([oo, on], ignore_index=True).to_csv(o_path, index=False)
    write_readings()


def write_readings():
    summ = pd.read_csv(OUT / "summary.csv")
    lines = []
    for s in summ.to_dict("records"):
        concl, desc = readings_of(s)
        lines.append(f"== {s['dataset']} {s['sub_screen']} / {s['space']} ==")
        lines += [f"CONCLUSION  {c}" for c in concl] + [f"            {d}" for d in desc]
        lines.append(f"            exceed k = {s['exceed_k']} of N_f = {s['floor_available_N']}; "
                     f"chance expectation {s['chance_expected']:.2f}; one-sided binomial p = {s['binom_p']:.3g}")
    (OUT / "readings.txt").write_text("\n".join(lines) + "\n", encoding="utf-8")
    log("\n".join(lines))


def verify():
    """Recompute every gated quantity from tracked per-perturbation CSVs alone."""
    summ = pd.read_csv(OUT / "summary.csv")
    bad = []
    for s in summ.to_dict("records"):
        pt = pd.read_csv(OUT / f"{s['dataset'].lower()}_{s['sub_screen'].lower()}_perturbations.csv")
        fa = pt[~pt.floor_unavailable]
        exc = fa.ratio > fa.floor_p95
        k, Nf = int(exc.sum()), len(fa)
        p = float(sps.binomtest(k, Nf, CHANCE, alternative="greater").pvalue) if Nf else np.nan
        ge2 = int((fa.ratio[exc] >= 2).sum())
        lt2 = int((fa.ratio[exc] < 2).sum())
        for name, a, b in (("exceed_k", k, s["exceed_k"]), ("N_f", Nf, s["floor_available_N"]),
                           ("exceed_ge2", ge2, s["exceed_ge2"]), ("exceed_lt2", lt2, s["exceed_lt2"])):
            if a != b:
                bad.append(f"{s['dataset']} {s['sub_screen']} {name}: csv {a} vs summary {b}")
        if not np.isclose(p, s["binom_p"], rtol=1e-12, atol=0):
            bad.append(f"{s['dataset']} {s['sub_screen']} binom_p {p} vs {s['binom_p']}")
        log(f"[verify] {s['dataset']} {s['sub_screen']}: k={k} N_f={Nf} p={p:.3g} ge2={ge2} lt2={lt2}")
    if bad:
        raise SystemExit("verify FAILED: " + "; ".join(bad))
    log("[verify] all readings re-derived from tracked CSVs; 0 mismatches")


def diagnostic():
    """NOT pre-registered; gates nothing. From tracked CSVs only. If sampling noise adds
    roughly orthogonally to an effect, ||d_measured||^2 ~ ||d_effect||^2 + ||noise||^2, so
    sqrt(max(ratio^2 - floor_p50(n)^2, 0)) is a rough noise-corrected length, in gap units."""
    rows = []
    for s in pd.read_csv(OUT / "summary.csv").to_dict("records"):
        stem = f"{s['dataset'].lower()}_{s['sub_screen'].lower()}"
        p = pd.read_csv(OUT / f"{stem}_perturbations.csv")
        f = pd.read_csv(OUT / f"{stem}_floor.csv").set_index("n")
        p["floor_p50"] = f.loc[p.floor_n_used.astype(int), "p50"].to_numpy()
        p["ratio_over_p95"] = p.ratio / p.floor_p95
        p["noise_corrected"] = np.sqrt(np.clip(p.ratio ** 2 - p.floor_p50 ** 2, 0, None))
        p[["dataset", "sub_screen", "target", "cells", "ratio", "floor_p50", "floor_p95",
           "ratio_over_p95", "noise_corrected"]].to_csv(OUT / f"{stem}_diagnostic.csv", index=False)
        nc = p.noise_corrected
        rows.append(dict(dataset=s["dataset"], sub_screen=s["sub_screen"],
                         ratio_over_p95_median=float(p.ratio_over_p95.median()),
                         ratio_over_p95_max=float(p.ratio_over_p95.max()),
                         noise_corrected_median=float(nc.median()),
                         noise_corrected_q1=float(nc.quantile(.25)),
                         noise_corrected_q3=float(nc.quantile(.75)),
                         noise_corrected_max=float(nc.max()),
                         floor_p50_at_30=float(f.loc[MIN_CELLS, "p50"]),
                         floor_p95_at_30=float(f.loc[MIN_CELLS, "p95"])))
    pd.DataFrame(rows).to_csv(OUT / "diagnostic_summary.csv", index=False)
    log(pd.DataFrame(rows).to_string(index=False))


GSE314342_FTP ="https://ftp.ncbi.nlm.nih.gov/geo/series/GSE314nnn/GSE314342/suppl/"
WORK_DISK_GB, WORK_HOURS, WORK_RAM_GB = 40.0, 24.0, 15.6   # D13


def gse314342_workability():
    """D13: time the first three per-sample files (downloaded together, as streaming would),
    extrapolate to the whole deposit. The files are not opened and are deleted afterwards."""
    import shutil
    import subprocess
    import urllib.request
    d = RAW / "gse314342"
    d.mkdir(parents=True, exist_ok=True)
    fl = urllib.request.urlopen(GSE314342_FTP + "filelist.txt", timeout=60).read().decode()
    (OUT / "gse314342_filelist.txt").write_text(fl, encoding="utf-8")
    rows = [r.split("\t") for r in fl.splitlines()[1:] if r.startswith("File\t")]
    files = [(r[1], int(r[3])) for r in rows]
    total = sum(s for _, s in files)
    first = files[:3]
    log(f"[314342] {len(files)} per-sample files, {total / 1e9:.1f} GB; timing the first 3: "
        + ", ".join(f"{n} ({s / 1e9:.2f} GB)" for n, s in first))
    # Per-sample files live in each GSM's own directory, not the series suppl/ (which holds
    # them only inside RAW.tar); the first attempt used the series path and got error pages.
    def url(n):
        gsm = n.split("_")[0]
        return f"https://ftp.ncbi.nlm.nih.gov/geo/samples/{gsm[:-3]}nnn/{gsm}/suppl/{n}"

    t0 = time.time()
    procs = [subprocess.Popen(["curl", "-s", "-S", "-f", "-o", str(d / n), url(n)]) for n, _ in first]
    rc = [p.wait() for p in procs]
    secs = time.time() - t0
    got = sum((d / n).stat().st_size for n, _ in first if (d / n).exists())
    want = sum(s for _, s in first)
    if got != want or any(rc):
        for n, _ in first:
            (d / n).unlink(missing_ok=True)
        raise SystemExit(f"timing download incomplete ({got} of {want} bytes, curl exit {rc}); "
                         "no workability conclusion drawn")
    rate = got / secs
    est_h = total / rate / 3600
    free_gb = shutil.disk_usage(d).free / 1e9
    per_file_gb = max(s for _, s in files) / 1e9
    for n, _ in first:
        (d / n).unlink(missing_ok=True)
    res = dict(files=len(files), total_bytes=total, timed_files=[n for n, _ in first], curl_exit=rc,
               bytes_expected=want, bytes_received=got, seconds=secs, bytes_per_s=rate,
               est_hours_download_only=est_h, limit_hours=WORK_HOURS,
               peak_disk_gb_streaming=per_file_gb * 3, limit_disk_gb=WORK_DISK_GB, free_disk_gb=free_gb,
               a_disk_ok=per_file_gb * 3 <= WORK_DISK_GB, b_time_ok=est_h <= WORK_HOURS,
               c_ram="not assessed: (b) decides" if est_h > WORK_HOURS else "to assess")
    res["workable"] = bool(res["a_disk_ok"] and res["b_time_ok"] and got == want)
    dump(res, OUT / "gse314342_workability.json")
    log(f"[314342] received {got} of {want} bytes in {secs:.0f}s = {rate / 1e6:.2f} MB/s; "
        f"whole deposit (download only) ~{est_h:.1f} h against the {WORK_HOURS:.0f} h limit; "
        f"workable={res['workable']}")
    return res


def main():
    global LOG_FH
    OUT.mkdir(parents=True, exist_ok=True)
    mode = sys.argv[1] if len(sys.argv) > 1 else ""
    LOG_FH = open(OUT / f"ratio_run_{mode}.log", "a", encoding="utf-8")
    require_flag()
    try:
        if mode == "reproduce":
            reproduce()
        elif mode == "gse278572":
            if not (OUT / "reproduction.csv").exists():
                raise SystemExit("run `reproduce` first")
            run_gse278572()
        elif mode == "verify":
            verify()
        elif mode == "diagnostic":
            diagnostic()
        elif mode == "gse314342_workability":
            if not (OUT / "summary.csv").exists():
                raise SystemExit("GSE278572 must be complete first (flag D14)")
            gse314342_workability()
        else:
            raise SystemExit(__doc__)
    except Stop as e:
        log(str(e))
        (OUT / "STOP.txt").write_text(f"{e.key}\n{e}\n", encoding="utf-8")
        raise SystemExit(1)


if __name__ == "__main__":
    main()
