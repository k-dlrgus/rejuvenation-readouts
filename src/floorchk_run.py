"""FLOORCHK — do paper 1's two length claims survive a size-matched noise floor?

Pre-registration: results/floorchk/PREREG.flag, written before any statistic. This script
refuses to run without it and never rewrites it. It imports, unchanged, the normalisation
and coordinate functions of the pipelines that produced the stored numbers (toward_run,
seng_run, south6_run, same_run, lowdim_common, fibro2_common). It writes only under
results/floorchk/.

    python -I src/floorchk_run.py run      reproduce, floors, gates, corrected ratios, readings
    python -I src/floorchk_run.py verify   re-derive every reading from the tracked CSVs only

Order inside `run` (flag D5): reproduction (STOP on any gating miss) -> claim-1 floor ->
claim-1 gate -> claim-2 floor and 1/n fit -> claim-2 gate -> claim-2 gap spread -> corrected
ratios of the real units -> readings. Gate results are on disk before any real unit's
corrected ratio is computed.
"""
from __future__ import annotations

import json
import sys
import time
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import spearmanr

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from toward_run import ANCHORS_PATH, load_donor, load_gtex_z, sum_rows, tmm_logcpm_quiet  # noqa: E402
from seng_run import tmm_cache_from_panel, z_from_counts  # noqa: E402
from south6_run import KMAX, S2_KS, gather_sum, open_matched  # noqa: E402
from lowdim_common import fit_pca  # noqa: E402
from same_run import d0_fibroblast_idx, panel_z, state_idx  # noqa: E402
from fibro2_common import load_frozen_ruler  # noqa: E402
from fibro3_common import AGED_LINE, YOUNG_LINE  # noqa: E402

OUT = ROOT / "results" / "floorchk"
PREREG = OUT / "PREREG.flag"
SOUTH6 = ROOT / "results" / "south6"
SAME = ROOT / "results" / "same"
GENESPACE = ROOT / "results" / "genespace"

SEED_FLOOR = 20261027
SEED_GATE = 20261028
SEED_GAP = 20261029
N_DRAWS = 200
LENGTHS = (0.5, 1.0, 2.0, 3.0)
GATE_TOL = 0.20
GATE_PCTS = (5, 25, 50, 95)
GATE_SIZES_C1 = (59, 84, 106, 189)
LADDER = (154, 309, 617, 1235, 2470)
N_UNIT_C2 = 5832
FIT_TOL = 0.10

C1_SPACES = ["S1"] + [f"S2_k{k}" for k in S2_KS]
C2_SPACES = ["FULL", "MD", "AGE"]
C2_LABEL = {"FULL": "all genes (23,485)", "MD": "MD genes (205)", "AGE": "age genes (3,089)"}

LOG_PATH = OUT / "floorchk_run.log"


class Stop(Exception):
    def __init__(self, key, message):
        super().__init__(message)
        self.key = key
        self.message = message


def log(msg):
    line = f"[{time.strftime('%H:%M:%S')}] {msg}"
    print(line, flush=True)
    with open(LOG_PATH, "a", encoding="utf-8") as f:
        f.write(line + "\n")


def rel_err(a, b):
    a = np.asarray(a, float)
    b = np.asarray(b, float)
    return np.abs(a - b) / np.maximum(np.abs(b), 1e-300)


# ----------------------------------------------------------------------------------------
# Reproduction (flag R1-R4)
# ----------------------------------------------------------------------------------------

class Repro:
    def __init__(self):
        self.rows = []

    def add(self, check, expected, observed, ok, gating=True):
        self.rows.append(dict(check=check, expected=str(expected), observed=str(observed),
                              match=("yes" if ok else "NO") if gating else "-",
                              gating="yes" if gating else "no"))
        log(f"[repro] {check}: expected {expected} observed {observed} -> {'ok' if ok else 'MISMATCH'}")

    def failures(self):
        return [r for r in self.rows if r["gating"] == "yes" and r["match"] != "yes"]

    def write(self):
        pd.DataFrame(self.rows).to_csv(OUT / "reproduction.csv", index=False)


def claim1_load(rp):
    s1 = np.load(SOUTH6 / "stage1.npz", allow_pickle=True)
    names = [str(x) for x in s1["keep"]]
    d_f = np.asarray(s1["d_f"], np.float64)
    anc = np.load(ANCHORS_PATH, allow_pickle=True)
    c_old = np.asarray(anc["c_old"], float).ravel()
    c_young = np.asarray(anc["c_young"], float).ravel()
    gap_vec = c_old - c_young
    st = pd.read_csv(SOUTH6 / "stage3_full_table.csv",
                     usecols=["factor", "space", "pair", "d_over_v", "v_norm", "d_norm"])
    st = st[st.pair == "O1_Y1"]

    # R1
    ref = st[st.space == "S1"].set_index("factor").reindex(names)
    gap_s1 = float(np.linalg.norm(gap_vec))
    ratio_s1 = np.linalg.norm(d_f, axis=1) / gap_s1
    rp.add("R1 claim 1 factor count (stage1.npz vs stage3 S1 O1_Y1)", 1836,
           f"{len(names)} / {int(ref['d_over_v'].notna().sum())}",
           len(names) == 1836 and int(ref["d_over_v"].notna().sum()) == 1836)
    e = float(rel_err(ratio_s1, ref["d_over_v"].values).max())
    rp.add("R1 S1 per-factor ratio vs stored d_over_v, max rel err", "<=1e-9", f"{e:.2e}", e <= 1e-9)
    e = float(rel_err(gap_s1, ref["v_norm"].values[0]))
    rp.add("R1 S1 ||gap|| (anchors.npz, read not refit) vs stored v_norm", f"{ref['v_norm'].values[0]:.10f}",
           f"{gap_s1:.10f}", e <= 1e-9)
    mn = float(ratio_s1.min())
    rp.add("R1 S1 smallest ratio (published 2.10)", "2.0966 (rounds to 2.10)", f"{mn:.4f}",
           round(mn, 4) == 2.0966 and round(mn, 2) == 2.10)
    rp.add("R1 S1 factors with ratio >= 2 (published: all 1,836)", 1836, int((ratio_s1 >= 2).sum()),
           int((ratio_s1 >= 2).sum()) == 1836)

    # R2: panel, normalisation and per-cell sums
    C_f = np.asarray(s1["C_f"], np.float64)
    C_nt = np.asarray(s1["C_nt"], np.float64)
    C_pool = np.asarray(s1["C_pool"], np.float64)
    matched_ruler = np.asarray(s1["matched_ruler"], int)
    n_ruler = int(np.asarray(s1["n_ruler"]).ravel()[0])
    frozen = load_frozen_ruler()
    mu = np.asarray(frozen["mu"], float)[matched_ruler]
    sd = np.asarray(frozen["sd"], float)[matched_ruler]
    C_panel = np.vstack([C_f, C_nt, C_pool.reshape(1, -1)])
    logcpm, nf = tmm_logcpm_quiet(C_panel)
    sd_safe = np.where(sd < 1e-12, 1.0, sd)
    missing = (C_panel == 0).all(0)
    Z_m = (logcpm - mu) / sd_safe
    Z_m[:, missing] = 0.0
    Z = np.zeros((C_panel.shape[0], n_ruler))
    Z[:, matched_ruler] = Z_m
    n_f = C_f.shape[0]
    rp.add("R2 missing mask vs stored", "equal", f"{int(missing.sum())} vs {int(np.asarray(s1['missing']).sum())}",
           bool(np.array_equal(missing, np.asarray(s1["missing"], bool))))
    e1 = float(np.max(np.abs(Z[:n_f] - np.asarray(s1["z_f"]))))
    e2 = float(np.max(np.abs(Z[-1] - np.asarray(s1["z_pool"]))))
    e3 = float(np.max(np.abs((Z[:n_f] - Z[-1]) - d_f)))
    rp.add("R2 panel z_f rebuilt vs stored, max abs", "<=1e-8", f"{e1:.2e}", e1 <= 1e-8)
    rp.add("R2 panel z_pool rebuilt vs stored, max abs", "<=1e-8", f"{e2:.2e}", e2 <= 1e-8)
    rp.add("R2 d_f rebuilt vs stored, max abs", "<=1e-8", f"{e3:.2e}", e3 <= 1e-8)
    cache = tmm_cache_from_panel(C_panel, nf)
    zc, _, _, _ = z_from_counts(C_pool, cache, mu, sd, missing)
    e4 = float(np.max(np.abs(zc - Z_m[-1])))
    rp.add("R2 cache self-check z_from_counts(C_pool) vs panel, max abs", "<=1e-6", f"{e4:.2e}", e4 <= 1e-6)

    nnz = int(np.asarray(s1["nnz"]).ravel()[0])
    data_m, ind_m, ip_m = open_matched(nnz)
    n_matched = int(C_pool.size)
    nt_cells = np.asarray(s1["nt_cells"], np.int64)
    pool_sum = gather_sum(data_m, ind_m, ip_m, nt_cells, n_matched)
    rp.add("R2 sum over the 14,013 non-targeting cells == C_pool (exact)", "14013 cells, equal",
           f"{nt_cells.size} cells, max abs diff {float(np.max(np.abs(pool_sum - C_pool))):.1f}",
           nt_cells.size == 14013 and bool(np.array_equal(pool_sum, C_pool)))
    flat = np.asarray(s1["cell_flat"], np.int64)
    indptr = np.asarray(s1["cell_indptr"], np.int64)
    n_cells = np.diff(indptr)
    cpt = pd.read_csv(SOUTH6 / "stage0_cells_per_target.csv").set_index("target")
    cpt_n = cpt.reindex(names)["n_cells"].values
    rp.add("R2 cells per factor == stage0_cells_per_target (exact)", "1836 equal",
           f"{int((n_cells == cpt_n).sum())} equal", bool(np.array_equal(n_cells, cpt_n)))
    log("[repro] summing every factor's cells to check C_f (206,390 cells)")
    bad = 0
    for i in range(n_f):
        s = gather_sum(data_m, ind_m, ip_m, flat[indptr[i]:indptr[i + 1]], n_matched)
        if not np.array_equal(s, C_f[i]):
            bad += 1
    rp.add("R2 per-factor cell sums == C_f rows (exact)", "0 mismatches", f"{bad} mismatches", bad == 0)
    pct = np.percentile(n_cells, GATE_PCTS)
    rp.add("R2 gate sizes = P5/P25/P50/P95 of cells per factor", list(GATE_SIZES_C1), [float(x) for x in pct],
           bool(np.allclose(pct, GATE_SIZES_C1, rtol=0, atol=1e-9)))

    # R3: rebuilt S2 basis (flag D9)
    gtex = load_gtex_z(log)
    Zg = np.asarray(gtex["Z"], np.float64)
    center = Zg.mean(0)
    sp = fit_pca(Zg - center, KMAX)
    L = np.asarray(sp.L, np.float64)
    log(f"[S2] PCA rebuilt on GTEx donors only n={Zg.shape[0]} p={Zg.shape[1]} k={L.shape[1]}")
    del Zg, gtex
    DL = d_f @ L
    gL = gap_vec @ L
    gaps = {"S1": gap_s1}
    raw = {"S1": ratio_s1}
    for k in S2_KS:
        sname = f"S2_k{k}"
        g = float(np.linalg.norm(gL[:k]))
        r = np.linalg.norm(DL[:, :k], axis=1) / g
        refk = st[st.space == sname].set_index("factor").reindex(names)
        e = float(rel_err(r, refk["d_over_v"].values).max())
        rp.add(f"R3 {sname} per-factor ratio (rebuilt basis) vs stored, max rel err", "<=1e-6", f"{e:.2e}", e <= 1e-6)
        e = float(rel_err(g, refk["v_norm"].values[0]))
        rp.add(f"R3 {sname} ||gap|| vs stored v_norm", f"{refk['v_norm'].values[0]:.6f}", f"{g:.6f}", e <= 1e-6)
        rp.add(f"R3 {sname} stored min ratio; n >= 2 (descriptive)", "-",
               f"min {r.min():.4f}; {int((r >= 2).sum())} of 1836", True, gating=False)
        gaps[sname] = g
        raw[sname] = r
    return dict(names=names, d_f=d_f, gap_vec=gap_vec, gaps=gaps, raw=raw, L=L,
                C_pool=C_pool, cache=cache, mu=mu, sd=sd, missing=missing, matched_ruler=matched_ruler,
                n_ruler=n_ruler, n_matched=n_matched, nt_cells=nt_cells, n_cells=n_cells,
                mm=(data_m, ind_m, ip_m))


def claim2_load(rp):
    frozen = load_frozen_ruler()
    obs_a, Y_a = load_donor(AGED_LINE, log)
    obs_y, Y_y = load_donor(YOUNG_LINE, log)
    sp = np.load(SAME / "split_idx.npz")
    o_idx, y_idx, aged_b = sp["O"], sp["Y"], sp["aged_B"]
    s_idx = state_idx(obs_a, "PartialReprog")
    srev_idx = state_idx(obs_y, "PartialReprog")
    pl_idx = state_idx(obs_a, "Pluripotency")
    nr_idx = state_idx(obs_a, "NonReprog")
    d0_a = d0_fibroblast_idx(obs_a)
    d0_y = d0_fibroblast_idx(obs_y)
    got = (o_idx.size, y_idx.size, aged_b.size, s_idx.size, d0_a.size, d0_y.size)
    rp.add("R4 cells O / Y / aged half B / S / aged d0 / young d0", "(2469, 3851, 2470, 5832, 4939, 7703)",
           str(got), got == (2469, 3851, 2470, 5832, 4939, 7703))
    rp.add("R4 halves are disjoint and inside d0", "yes",
           str(bool(not (set(o_idx) & set(aged_b)) and set(o_idx) | set(aged_b) == set(d0_a))),
           bool(not (set(o_idx) & set(aged_b)) and set(o_idx) | set(aged_b) == set(d0_a)))
    counts = [sum_rows(Y_a, o_idx), sum_rows(Y_y, y_idx), sum_rows(Y_a, s_idx),
              sum_rows(Y_y, srev_idx), sum_rows(Y_a, pl_idx), sum_rows(Y_a, nr_idx)]
    Z, missing, nf, C = panel_z(counts, frozen)
    stored = np.load(SAME / "panel_z.npz", allow_pickle=True)
    e = float(np.max(np.abs(Z - np.asarray(stored["Z"]))))
    rp.add("R4 SAME main panel z rebuilt vs stored panel_z.npz, max abs", "<=1e-8", f"{e:.2e}", e <= 1e-8)
    gm = pd.read_csv(GENESPACE / "gene_map.csv")
    gm = gm[gm.record == "ruler_column"].sort_values("col")
    cols = {"FULL": np.arange(Z.shape[1]),
            "MD": gm.loc[gm.in_MD.astype(bool), "col"].astype(int).values,
            "AGE": gm.loc[gm.in_AGE.astype(bool), "col"].astype(int).values}
    rp.add("R4 space sizes FULL / MD / AGE", "(23485, 205, 3089)",
           str(tuple(len(cols[s]) for s in C2_SPACES)),
           tuple(len(cols[s]) for s in C2_SPACES) == (23485, 205, 3089))
    ss = pd.read_csv(GENESPACE / "same_stats.csv")
    ss = ss[ss.test == "forward"].set_index("space")
    d = Z[2] - Z[0]
    gap = Z[1] - Z[0]
    raw, dn, gn = {}, {}, {}
    for s in C2_SPACES:
        dn[s] = float(np.linalg.norm(d[cols[s]]))
        gn[s] = float(np.linalg.norm(gap[cols[s]]))
        raw[s] = dn[s] / gn[s]
        e = float(max(rel_err(dn[s], ss.loc[s, "d_norm"]), rel_err(gn[s], ss.loc[s, "v_norm"])))
        rp.add(f"R4 {s} ||d||, ||gap|| vs genespace same_stats, max rel err", "<=1e-9",
               f"{e:.2e} (||d|| {dn[s]:.4f}, ||gap|| {gn[s]:.4f}, ratio {raw[s]:.4f})", e <= 1e-9)
    rp.add("R4 published ratios: MD 2.3 (main text); AGE, FULL > 2 (Limitations)", "2.30 / >2 / >2",
           f"{raw['MD']:.3f} / {raw['AGE']:.3f} / {raw['FULL']:.3f}",
           round(raw["MD"], 1) == 2.3 and raw["AGE"] > 2 and raw["FULL"] > 2)
    cache = tmm_cache_from_panel(C, nf)
    mu = np.asarray(frozen["mu"], float)
    sd = np.asarray(frozen["sd"], float)
    for row, nm in ((0, "O"), (2, "S")):
        zc, _, _, _ = z_from_counts(C[row], cache, mu, sd, missing)
        e = float(np.max(np.abs(zc - Z[row])))
        rp.add(f"R4 cache self-check z_from_counts({nm}) vs panel, max abs", "<=1e-6", f"{e:.2e}", e <= 1e-6)
    return dict(Y_a=Y_a, Y_y=Y_y, aged_b=np.asarray(aged_b, int), d0_a=d0_a, d0_y=d0_y,
                z_O=Z[0].copy(), d=d, gap=gap, cols=cols, raw=raw, dn=dn, gn=gn, cache=cache,
                mu=mu, sd=sd, missing=missing)


# ----------------------------------------------------------------------------------------
# Claim 1: floor, gate, corrected ratios
# ----------------------------------------------------------------------------------------

def c1_z(c1, counts):
    zm, _, _, _ = z_from_counts(counts, c1["cache"], c1["mu"], c1["sd"], c1["missing"])
    full = np.zeros(c1["n_ruler"])
    full[c1["matched_ruler"]] = zm
    return full


def c1_d0(c1, cells):
    data_m, ind_m, ip_m = c1["mm"]
    c_draw = gather_sum(data_m, ind_m, ip_m, cells, c1["n_matched"])
    c_ref = c1["C_pool"] - c_draw
    return c1_z(c1, c_draw) - c1_z(c1, c_ref)


def c1_sq_norms(c1, d0):
    """Squared norm of d0 in S1 and in each S2_k; also the projection for the gate."""
    p = d0 @ c1["L"]
    out = {"S1": float(d0 @ d0)}
    for k in S2_KS:
        out[f"S2_k{k}"] = float(p[:k] @ p[:k])
    return out, p


def claim1_floor(c1):
    sizes = sorted(set(int(x) for x in np.unique(c1["n_cells"])) | set(GATE_SIZES_C1))
    log(f"[c1 floor] {len(sizes)} sizes x {N_DRAWS} draws (seed {SEED_FLOOR}, stream [seed, 1, n])")
    sq = {s: np.zeros((len(sizes), N_DRAWS)) for s in C1_SPACES}
    for si, n in enumerate(sizes):
        rng = np.random.default_rng([SEED_FLOOR, 1, n])
        for j in range(N_DRAWS):
            cells = rng.choice(c1["nt_cells"], size=n, replace=False)
            q, _ = c1_sq_norms(c1, c1_d0(c1, cells))
            for s in C1_SPACES:
                sq[s][si, j] = q[s]
        if si % 20 == 0 or si == len(sizes) - 1:
            log(f"[c1 floor] size {n} ({si + 1}/{len(sizes)}) S1 median ratio "
                f"{np.median(np.sqrt(sq['S1'][si])) / c1['gaps']['S1']:.4f}")
    np.savez_compressed(OUT / "claim1_floor_draws.npz", sizes=np.array(sizes),
                        **{f"sq_{s}": sq[s] for s in C1_SPACES})
    factor_sizes = set(int(x) for x in c1["n_cells"])
    rows = []
    for s in C1_SPACES:
        g = c1["gaps"][s]
        for si, n in enumerate(sizes):
            r = np.sqrt(sq[s][si]) / g
            rows.append(dict(space=s, n=n, is_factor_size=n in factor_sizes,
                             n_factors=int((c1["n_cells"] == n).sum()), gap=g,
                             floor_median=float(np.median(r)), floor_p95=float(np.percentile(r, 95)),
                             m_d0sq=float(np.median(sq[s][si])), frac_ge2=float((r >= 2).mean())))
    df = pd.DataFrame(rows)
    df.to_csv(OUT / "claim1_floor.csv", index=False)
    return df


def claim1_gate(c1, floor):
    S1_coords = c1["matched_ruler"][~c1["missing"]]
    rows = []
    for n in GATE_SIZES_C1:
        rng = np.random.default_rng([SEED_GATE, 1, n])
        m = {s: float(floor[(floor.space == s) & (floor.n == n)]["m_d0sq"].iloc[0]) for s in C1_SPACES}
        cr = {(s, L): [] for s in C1_SPACES for L in LENGTHS}
        for j in range(N_DRAWS):
            cells = rng.choice(c1["nt_cells"], size=n, replace=False)
            d0 = c1_d0(c1, cells)
            p = d0 @ c1["L"]
            for L in LENGTHS:
                for s in C1_SPACES:
                    g = c1["gaps"][s]
                    if s == "S1":
                        u = rng.standard_normal(S1_coords.size)
                        u /= np.linalg.norm(u)
                        x = d0.copy()
                        x[S1_coords] += L * g * u
                    else:
                        k = int(s.split("k")[1])
                        u = rng.standard_normal(k)
                        u /= np.linalg.norm(u)
                        x = p[:k] + L * g * u
                    cr[(s, L)].append(np.sqrt(max(0.0, float(x @ x) - m[s])) / g)
        for s in C1_SPACES:
            for L in LENGTHS:
                v = np.asarray(cr[(s, L)])
                med = float(np.median(v))
                rows.append(dict(space=s, n=n, planted=L, median_cr=med, q25=float(np.percentile(v, 25)),
                                 q75=float(np.percentile(v, 75)), rel_err=(med - L) / L,
                                 gating=L >= 1, within_20pct=abs(med - L) <= GATE_TOL * L))
        log(f"[c1 gate] size {n} done: S1 medians " +
            ", ".join(f"L{L}={np.median(cr[('S1', L)]):.3f}" for L in LENGTHS))
    df = pd.DataFrame(rows)
    df.to_csv(OUT / "claim1_gate.csv", index=False)
    summ = gate_summary(df, C1_SPACES)
    summ.to_csv(OUT / "claim1_gate_summary.csv", index=False)
    return df, summ


def gate_summary(df, spaces):
    out = []
    for s in spaces:
        g = df[(df.space == s) & (df.gating)]
        fails = g[~g.within_20pct.astype(bool)]
        out.append(dict(space=s, gate_pass=bool(len(fails) == 0), n_gating_cells=len(g),
                        n_failing=len(fails),
                        failing=";".join(f"n={int(r.n)} L={r.planted:g} med={r.median_cr:.3f}"
                                         for r in fails.itertuples())))
    return pd.DataFrame(out)


def claim1_units(c1, floor):
    rows = []
    for s in C1_SPACES:
        f = floor[floor.space == s].set_index("n")
        g = c1["gaps"][s]
        dn = c1["raw"][s] * g
        m = f.loc[c1["n_cells"], "m_d0sq"].values
        cr = np.sqrt(np.maximum(0.0, dn ** 2 - m)) / g
        for i, nm in enumerate(c1["names"]):
            rows.append(dict(factor=nm, space=s, n_cells=int(c1["n_cells"][i]), d_norm=dn[i], gap=g,
                             raw_ratio=c1["raw"][s][i], m_d0sq=m[i],
                             floor_median=f.loc[c1["n_cells"][i], "floor_median"],
                             floor_p95=f.loc[c1["n_cells"][i], "floor_p95"], corrected_ratio=cr[i]))
    df = pd.DataFrame(rows)
    df.to_csv(OUT / "claim1_units.csv", index=False)
    return df


# ----------------------------------------------------------------------------------------
# Claim 2: floor ladder, fit, gate, gap spread, corrected ratio
# ----------------------------------------------------------------------------------------

def c2_z(c2, counts):
    z, _, _, _ = z_from_counts(counts, c2["cache"], c2["mu"], c2["sd"], c2["missing"])
    return z


def claim2_floor(c2):
    rows, sq_all = [], {}
    for n in LADDER:
        rng = np.random.default_rng([SEED_FLOOR, 2, n])
        sq = {s: np.zeros(N_DRAWS) for s in C2_SPACES}
        sets = set()
        for j in range(N_DRAWS):
            cells = np.sort(rng.choice(c2["aged_b"], size=n, replace=False))
            sets.add(cells.tobytes())
            d0 = c2_z(c2, sum_rows(c2["Y_a"], cells)) - c2["z_O"]
            for s in C2_SPACES:
                v = d0[c2["cols"][s]]
                sq[s][j] = float(v @ v)
        for s in C2_SPACES:
            r = np.sqrt(sq[s]) / c2["gn"][s]
            rows.append(dict(space=s, n=n, n_distinct_draws=len(sets), gap=c2["gn"][s],
                             floor_median=float(np.median(r)), floor_p95=float(np.percentile(r, 95)),
                             m_d0sq=float(np.median(sq[s])), frac_ge2=float((r >= 2).mean())))
            sq_all[f"sq_{s}_{n}"] = sq[s]
        log(f"[c2 floor] n={n} distinct draws {len(sets)}; FULL median ratio "
            f"{rows[-3]['floor_median']:.4f}")
    np.savez_compressed(OUT / "claim2_floor_draws.npz", **sq_all)
    df = pd.DataFrame(rows)
    df.to_csv(OUT / "claim2_floor.csv", index=False)
    return df


def claim2_fit(floor):
    rows = []
    for s in C2_SPACES:
        f = floor[floor.space == s].sort_values("n")
        n = f.n.values.astype(float)
        m = f.m_d0sq.values
        A = np.column_stack([np.ones_like(n), 1.0 / n])
        (a, b), *_ = np.linalg.lstsq(A, m, rcond=None)
        fit = a + b / n
        resid = (m - fit) / fit
        ok = bool(np.all(np.abs(resid) <= FIT_TOL))
        rec = dict(space=s, a=float(a), b=float(b), fit_pass=ok, F_ext=float(a + b / N_UNIT_C2),
                   F_obs=float(m[-1]), max_abs_rel_resid=float(np.max(np.abs(resid))),
                   extrapolation_factor=N_UNIT_C2 / float(n.max()))
        for nn, rr, mm, ff in zip(f.n.values, resid, m, fit):
            rec[f"m_{int(nn)}"] = float(mm)
            rec[f"fit_{int(nn)}"] = float(ff)
            rec[f"resid_{int(nn)}"] = float(rr)
        rows.append(rec)
        log(f"[c2 fit] {s} a={a:.4f} b={b:.2f} max|rel resid|={np.max(np.abs(resid)):.4f} pass={ok} "
            f"F_ext={a + b / N_UNIT_C2:.4f} F_obs={m[-1]:.4f}")
    df = pd.DataFrame(rows)
    df.to_csv(OUT / "claim2_fit.csv", index=False)
    return df


def claim2_gate(c2, floor):
    coords = {s: np.asarray([c for c in c2["cols"][s] if not c2["missing"][c]], int) for s in C2_SPACES}
    rows = []
    for n in LADDER:
        rng = np.random.default_rng([SEED_GATE, 2, n])
        m = {s: float(floor[(floor.space == s) & (floor.n == n)]["m_d0sq"].iloc[0]) for s in C2_SPACES}
        cr = {(s, L): [] for s in C2_SPACES for L in LENGTHS}
        for j in range(N_DRAWS):
            cells = np.sort(rng.choice(c2["aged_b"], size=n, replace=False))
            d0 = c2_z(c2, sum_rows(c2["Y_a"], cells)) - c2["z_O"]
            for L in LENGTHS:
                for s in C2_SPACES:
                    g = c2["gn"][s]
                    u = rng.standard_normal(coords[s].size)
                    u /= np.linalg.norm(u)
                    x = d0.copy()
                    x[coords[s]] += L * g * u
                    v = x[c2["cols"][s]]
                    cr[(s, L)].append(np.sqrt(max(0.0, float(v @ v) - m[s])) / g)
        for s in C2_SPACES:
            for L in LENGTHS:
                v = np.asarray(cr[(s, L)])
                med = float(np.median(v))
                rows.append(dict(space=s, n=n, planted=L, median_cr=med, q25=float(np.percentile(v, 25)),
                                 q75=float(np.percentile(v, 75)), rel_err=(med - L) / L,
                                 gating=L >= 1, within_20pct=abs(med - L) <= GATE_TOL * L))
        log(f"[c2 gate] n={n} FULL medians " + ", ".join(f"L{L}={np.median(cr[('FULL', L)]):.3f}" for L in LENGTHS))
    df = pd.DataFrame(rows)
    df.to_csv(OUT / "claim2_gate.csv", index=False)
    summ = gate_summary(df, C2_SPACES)
    summ.to_csv(OUT / "claim2_gate_summary.csv", index=False)
    return df, summ


def claim2_gap_spread(c2):
    rng = np.random.default_rng(SEED_GAP)
    na, ny = c2["d0_a"].size // 2, c2["d0_y"].size // 2
    vals = {s: [] for s in C2_SPACES}
    for j in range(N_DRAWS):
        pa = rng.permutation(c2["d0_a"].size)
        py = rng.permutation(c2["d0_y"].size)
        o = np.sort(c2["d0_a"][pa[:na]])
        y = np.sort(c2["d0_y"][py[:ny]])
        g = c2_z(c2, sum_rows(c2["Y_y"], y)) - c2_z(c2, sum_rows(c2["Y_a"], o))
        for s in C2_SPACES:
            vals[s].append(float(np.linalg.norm(g[c2["cols"][s]])))
    rows = []
    for s in C2_SPACES:
        v = np.asarray(vals[s])
        rows.append(dict(space=s, observed_gap=c2["gn"][s], n_aged_half=na, n_young_half=ny,
                         median=float(np.median(v)), p5=float(np.percentile(v, 5)),
                         p95=float(np.percentile(v, 95)), min=float(v.min()), max=float(v.max()),
                         p95_minus_p5_over_observed=float((np.percentile(v, 95) - np.percentile(v, 5)) / c2["gn"][s])))
    df = pd.DataFrame(rows)
    df.to_csv(OUT / "claim2_gap_spread.csv", index=False)
    log("[c2 gap] " + df[["space", "observed_gap", "median", "p5", "p95"]].to_string(index=False).replace("\n", " | "))
    return df


def claim2_unit(c2, fit, floor):
    rows = []
    for s in C2_SPACES:
        fr = fit.set_index("space").loc[s]
        g, dn = c2["gn"][s], c2["dn"][s]
        cr_obs = np.sqrt(max(0.0, dn ** 2 - fr.F_obs)) / g
        cr_ext = np.sqrt(max(0.0, dn ** 2 - fr.F_ext)) / g if bool(fr.fit_pass) else np.nan
        f2470 = floor[(floor.space == s) & (floor.n == 2470)].iloc[0]
        rows.append(dict(space=s, unit="GM00731 PartialReprog", n_cells=N_UNIT_C2, d_norm=dn, gap=g,
                         raw_ratio=c2["raw"][s], fit_pass=bool(fr.fit_pass), F_ext=fr.F_ext, F_obs=fr.F_obs,
                         floor_ratio_ext=np.sqrt(max(fr.F_ext, 0.0)) / g, floor_ratio_obs_2470=f2470.floor_median,
                         corrected_ratio_ext=cr_ext, corrected_ratio_obs=cr_obs))
    df = pd.DataFrame(rows)
    df.to_csv(OUT / "claim2_unit.csv", index=False)
    return df


# ----------------------------------------------------------------------------------------
# Readings (flag D7) — derived from the CSVs only, so `verify` reproduces them exactly
# ----------------------------------------------------------------------------------------

def readings_from_tables(t):
    out = []
    c1f, c1g, c1u = t["claim1_floor"], t["claim1_gate_summary"], t["claim1_units"]
    c2f, c2g, c2u = t["claim2_floor"], t["claim2_gate_summary"], t["claim2_unit"]

    # Claim 1
    for s in C1_SPACES:
        out.append(f"== claim 1 / {s}{' (primary)' if s == 'S1' else ' (secondary)'} ==")
        f = c1f[(c1f.space == s) & (c1f.is_factor_size.astype(bool))]
        n2 = int((f.floor_p95 >= 2).sum())
        if n2 >= 1:
            out.append(f"CONCLUSION  FLOOR-REACHES-2: in claim 1, {s}, the control floor's P95 is >= 2 at {n2} "
                       f"of {len(f)} group sizes; at those depths a ratio of 2 is reached by sampling noise alone.")
        else:
            out.append(f"            (descriptive) floor P95 >= 2 at 0 of {len(f)} group sizes; FLOOR-REACHES-2 did not fire.")
        gp = bool(c1g.set_index("space").loc[s, "gate_pass"])
        u = c1u[c1u.space == s]
        if not gp:
            out.append(f"CONCLUSION  LENGTH-UNRESOLVED: in claim 1, {s}, the correction failed its gate; only raw "
                       f"ratio and floor are reported and no statement about the published claim is made.")
            if s == "S1":
                out.append("            CLAIM1-* not fired: the S1 gate failed.")
            continue
        n_st = int((u.corrected_ratio >= 2).sum())
        n_dp = int(((u.raw_ratio >= 2) & (u.corrected_ratio < 2)).sum())
        if n_st >= 1:
            out.append(f"CONCLUSION  LENGTH-STANDS: in claim 1, {s}, {n_st} of 1836 units have corrected ratio >= 2; "
                       f"the published 'more than twice' holds for them after removing the noise contribution.")
        else:
            out.append("            (descriptive) 0 of 1836 units have corrected ratio >= 2; LENGTH-STANDS did not fire.")
        if n_dp >= 1:
            out.append(f"CONCLUSION  LENGTH-DEPTH: in claim 1, {s}, {n_dp} of 1836 units have raw ratio >= 2 but "
                       f"corrected ratio < 2; for them the published 'more than twice' was reached partly by sampling noise.")
        else:
            out.append(f"            (descriptive) 0 of 1836 units have raw >= 2 and corrected < 2; LENGTH-DEPTH did not fire.")
        if s == "S1":
            if n_st == 1836:
                out.append("CONCLUSION  CLAIM1-HOLDS")
            elif n_st >= 1:
                out.append(f"CONCLUSION  CLAIM1-PARTIAL: {n_st} of 1,836")
            else:
                out.append("CONCLUSION  CLAIM1-FAILS")

    # Claim 2
    for s in C2_SPACES:
        lab = C2_LABEL[s]
        out.append(f"== claim 2 / {lab}{' (primary)' if s == 'FULL' else ' (secondary)'} ==")
        f = c2f[c2f.space == s]
        n2 = int((f.floor_p95 >= 2).sum())
        if n2 >= 1:
            out.append(f"CONCLUSION  FLOOR-REACHES-2: in claim 2, {lab}, the control floor's P95 is >= 2 at {n2} "
                       f"of {len(f)} group sizes; at those depths a ratio of 2 is reached by sampling noise alone.")
        else:
            out.append(f"            (descriptive) floor P95 >= 2 at 0 of {len(f)} group sizes; FLOOR-REACHES-2 did not fire.")
        gp = bool(c2g.set_index("space").loc[s, "gate_pass"])
        u = c2u.set_index("space").loc[s]
        if not gp:
            out.append(f"CONCLUSION  LENGTH-UNRESOLVED: in claim 2, {lab}, the correction failed its gate; only raw "
                       f"ratio and floor are reported and no statement about the published claim is made.")
            if s == "FULL":
                out.append("            CLAIM2-* not fired: the all-genes gate failed.")
            continue
        fit_ok = bool(u.fit_pass)
        if fit_ok:
            ge = [u.corrected_ratio_ext >= 2, u.corrected_ratio_obs >= 2]
            basis = (f"corrected ratio {u.corrected_ratio_ext:.4f} (floor extrapolated to 5,832) and "
                     f"{u.corrected_ratio_obs:.4f} (observed 2,470 floor)")
        else:
            ge = [u.corrected_ratio_obs >= 2]
            basis = (f"corrected ratio {u.corrected_ratio_obs:.4f} (observed 2,470 floor alone: the 1/n fit "
                     f"failed its 10% residual check, so no extrapolation)")
        if all(ge):
            verdict = "stands"
        elif not any(ge):
            verdict = "below"
        else:
            verdict = "split"
        out.append(f"            {basis}; raw ratio {u.raw_ratio:.4f}")
        if verdict == "stands":
            out.append(f"CONCLUSION  LENGTH-STANDS: in claim 2, {lab}, 1 of 1 units have corrected ratio >= 2; "
                       f"the published 'more than twice' holds for them after removing the noise contribution.")
        elif verdict == "below":
            if u.raw_ratio >= 2:
                out.append(f"CONCLUSION  LENGTH-DEPTH: in claim 2, {lab}, 1 of 1 units have raw ratio >= 2 but "
                           f"corrected ratio < 2; for them the published 'more than twice' was reached partly by sampling noise.")
            else:
                out.append("            (descriptive) raw ratio < 2 and corrected < 2; neither LENGTH-STANDS nor LENGTH-DEPTH fires.")
        else:
            out.append(f"CONCLUSION  LENGTH-UNRESOLVED: in claim 2, {lab}, the corrected ratio is "
                       f"{u.corrected_ratio_ext:.3f} with the floor extrapolated to 5,832 cells and "
                       f"{u.corrected_ratio_obs:.3f} with the observed 2,470-cell floor, on opposite sides of 2; "
                       f"the verdict depends on the extrapolation and no statement about the published claim is made.")
        if s == "FULL":
            if verdict == "stands":
                out.append("CONCLUSION  CLAIM2-HOLDS")
            elif verdict == "below":
                out.append("CONCLUSION  CLAIM2-FAILS")
            else:
                out.append("            CLAIM2-* not fired: the two floors disagree (LENGTH-UNRESOLVED above).")
    return out


def load_tables():
    names = ["claim1_floor", "claim1_gate_summary", "claim1_units", "claim2_floor", "claim2_gate_summary", "claim2_unit"]
    return {n: pd.read_csv(OUT / f"{n}.csv") for n in names}


def write_readings(lines):
    (OUT / "readings.txt").write_text("\n".join(lines) + "\n", encoding="utf-8")
    for ln in lines:
        log("[reading] " + ln)


def descriptive(c1u):
    rows = []
    for s in C1_SPACES:
        u = c1u[c1u.space == s]
        rho_raw = spearmanr(u.raw_ratio, u.n_cells).correlation
        rho_cr = spearmanr(u.corrected_ratio, u.n_cells).correlation
        rows.append(dict(claim="claim 1", space=s, n_units=len(u), raw_min=u.raw_ratio.min(),
                         raw_median=u.raw_ratio.median(), raw_n_ge2=int((u.raw_ratio >= 2).sum()),
                         cr_min=u.corrected_ratio.min(), cr_median=u.corrected_ratio.median(),
                         cr_max=u.corrected_ratio.max(), cr_n_ge2=int((u.corrected_ratio >= 2).sum()),
                         floor_median_at_own_n_median=u.floor_median.median(),
                         raw_over_floor_p95_median=(u.raw_ratio / u.floor_p95).median(),
                         n_raw_within_floor_p95=int((u.raw_ratio <= u.floor_p95).sum()),
                         spearman_raw_vs_n=rho_raw, spearman_cr_vs_n=rho_cr))
    df = pd.DataFrame(rows)
    df.to_csv(OUT / "claim1_summary.csv", index=False)
    return df


def run():
    if not PREREG.exists():
        raise SystemExit("results/floorchk/PREREG.flag is missing. Refusing to compute.")
    OUT.mkdir(parents=True, exist_ok=True)
    log("=" * 90)
    log(f"FLOORCHK run. PREREG.flag bytes={PREREG.stat().st_size}")
    rp = Repro()
    status = dict(stop=None)
    try:
        c1 = claim1_load(rp)
        c2 = claim2_load(rp)
        rp.write()
        fails = rp.failures()
        if fails:
            raise Stop("reproduction_failed", "; ".join(f"{r['check']}: expected {r['expected']} observed {r['observed']}"
                                                         for r in fails))
        log("[repro] all gating checks passed")
        c1_floor = claim1_floor(c1)
        _, c1_gs = claim1_gate(c1, c1_floor)
        log("[c1 gate] " + c1_gs.to_string(index=False).replace("\n", " | "))
        c2_floor = claim2_floor(c2)
        c2_fit = claim2_fit(c2_floor)
        _, c2_gs = claim2_gate(c2, c2_floor)
        log("[c2 gate] " + c2_gs.to_string(index=False).replace("\n", " | "))
        claim2_gap_spread(c2)
        log("[units] gates are on disk; computing corrected ratios of the real units")
        c1u = claim1_units(c1, c1_floor)
        claim2_unit(c2, c2_fit, c2_floor)
        descriptive(c1u)
        write_readings(readings_from_tables(load_tables()))
        status["stop"] = None
    except Stop as e:
        rp.write()
        status["stop"] = e.key
        status["message"] = e.message
        (OUT / "readings.txt").write_text(f"STOP {e.key}: {e.message}\n", encoding="utf-8")
        log(f"STOP {e.key}: {e.message}")
    (OUT / "summary.json").write_text(json.dumps(status, indent=1), encoding="utf-8")


def verify():
    t = load_tables()
    want = readings_from_tables(t)
    have = (OUT / "readings.txt").read_text(encoding="utf-8").rstrip("\n").split("\n")
    mism = [(a, b) for a, b in zip(want, have) if a != b] + ([("length", f"{len(want)} vs {len(have)}")]
                                                            if len(want) != len(have) else [])
    # gate summaries re-derived from the per-cell gate tables
    for c, spaces in (("claim1", C1_SPACES), ("claim2", C2_SPACES)):
        g = pd.read_csv(OUT / f"{c}_gate.csv")
        g["within_20pct"] = (g.median_cr - g.planted).abs() <= GATE_TOL * g.planted
        s1 = gate_summary(g, spaces).set_index("space").gate_pass
        s0 = t[f"{c}_gate_summary"].set_index("space").gate_pass.astype(bool)
        if not s1.equals(s0):
            mism.append((f"{c} gate", f"{s1.to_dict()} vs {s0.to_dict()}"))
    print(f"verify: {len(want)} reading lines, {len(mism)} mismatches")
    for m in mism:
        print("  MISMATCH", m)
    return len(mism)


def diagnostic():
    """NOT pre-registered; gates nothing. Chance baseline for the claim-1 LENGTH-STANDS count:
    treat each floor draw as a unit of that size with no effect and apply the same correction.
    Expected count = sum over the 1,836 factors of P_null(CR >= 2 | own n). Reads the floor draws
    (results/floorchk/claim1_floor_draws.npz, regenerable) and claim1_units.csv."""
    fd = np.load(OUT / "claim1_floor_draws.npz")
    sizes = fd["sizes"]
    u = pd.read_csv(OUT / "claim1_units.csv")
    fl = pd.read_csv(OUT / "claim1_floor.csv")
    rows = []
    for s in C1_SPACES:
        sq = fd[f"sq_{s}"]
        g = float(fl[fl.space == s].gap.iloc[0])
        m = np.median(sq, axis=1, keepdims=True)
        cr0 = np.sqrt(np.maximum(0.0, sq - m)) / g
        p_ge2 = dict(zip(sizes.tolist(), (cr0 >= 2).mean(1)))
        p95 = dict(zip(sizes.tolist(), np.percentile(cr0, 95, axis=1)))
        x = u[u.space == s]
        exp = float(sum(p_ge2[int(n)] for n in x.n_cells))
        obs = int((x.corrected_ratio >= 2).sum())
        st = x[x.corrected_ratio >= 2]
        rows.append(dict(space=s, observed_cr_ge2=obs, expected_cr_ge2_if_all_null=exp,
                         null_cr_p95_at_n59=p95[59], null_cr_p95_at_n106=p95[106], null_cr_p95_at_n189=p95[189],
                         sizes_where_null_cr_p95_ge2=int(sum(v >= 2 for v in p95.values())),
                         standing_median_n=float(st.n_cells.median()) if len(st) else np.nan,
                         all_median_n=float(x.n_cells.median()),
                         standing_exceeding_own_null_cr_p95=int(sum(r.corrected_ratio > p95[int(r.n_cells)]
                                                                    for r in st.itertuples()))))
    df = pd.DataFrame(rows)
    df.to_csv(OUT / "claim1_diagnostic.csv", index=False)
    print(df.T.to_string())


if __name__ == "__main__":
    mode = sys.argv[1] if len(sys.argv) > 1 else "run"
    if mode == "run":
        run()
    elif mode == "diagnostic":
        diagnostic()
    elif mode == "verify":
        sys.exit(1 if verify() else 0)
    else:
        raise SystemExit(f"unknown mode {mode!r}")
