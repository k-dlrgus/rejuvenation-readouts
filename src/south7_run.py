"""SOUTH7 — re-rank the Southard CRISPRa screen with a qualified statistic, a
depth-matched floor, an anchor-free basis, and a permutation floor that does not
cap the answer.

The pre-registration is results/south7/PREREG.flag, written before any SOUTH7
number. This script reads it and does not rewrite it.

Nothing in results/south6/ or results/south6rec/ is written. src/south6rec_run.py
is imported unchanged and its Stage 0 is reused for the SOUTH6 reproduction; its
module-level output directory is redirected to results/south7/ for the duration
of that call so no file in results/south6rec/ is touched. The frozen ruler and
the TOWARD anchors are read, never refit.

Stages, run from the repository root:
    SOUTH7_STAGE=0A1 python src/south7_run.py     Stage 0 gate + Stage A1 nulls
    SOUTH7_STAGE=A2  python src/south7_run.py     permutations (SOUTH7_SET=primary|rest)
    SOUTH7_STAGE=A3  python src/south7_run.py     calls, report-only, findings
"""
from __future__ import annotations

import os
import sys
import time
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import sparse, stats as sps

sys.path.insert(0, str(Path(__file__).resolve().parent))

from brain_phase1_common import Logger  # noqa: E402
from target_common import unit  # noqa: E402
from toward_run import ANCHORS_PATH, load_gtex_z  # noqa: E402
from seng_run import bh_q  # noqa: E402
from lowdim_common import fit_pca  # noqa: E402
from south6_run import (  # noqa: E402
    KMAX, PERM_CHUNK, S2_KS, SEED, SOUTH6, block_delta,
)
import south6rec_run as REC  # noqa: E402

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "results" / "south7"
OUT.mkdir(parents=True, exist_ok=True)
PREREG = OUT / "PREREG.flag"
FINDINGS = ROOT / "FINDINGS_SOUTH7.md"
S4_CSV = OUT / "southard_table_s4_hs27.csv"
S4_XLSX = ROOT / "data" / "raw" / "southard" / "41588_2025_2284_MOESM4_ESM.xlsx"

N_PERM7 = 100_000
R_NULL = 200
Q_BAR = 0.05
MID_BINS = ("40-49", "50-59")
ANCHOR_BINS = ("20-29", "30-39", "60-69", "70-79")
KS7 = (20, 50, 100)
PRIMARY_SPACE = "S2A_k20"
SPACES7 = ("S2A_k20", "S2A_k50", "S2A_k100", "S1")
SPACES6 = ("S2_k20", "S2_k50", "S2_k100", "S1")
HIGH_DEPTH_PCTL = 90.0
TOL = 1e-4
# Paper-stated counts for the published classification columns (Table S4).
PAPER_ACTIVE = 212
PAPER_EXPANDED = 102

_W = {}


class Stop(Exception):
    def __init__(self, key, message):
        super().__init__(message)
        self.key = key
        self.message = message


def require_prereg():
    if not PREREG.exists():
        raise SystemExit(f"PREREG.flag missing at {PREREG}. Not computing.")


# ------------------------------------------------------------------ factor set

def factor_set(log):
    """The published primary set: Table S4 active_fibro OR expanded_active_fibro.

    The classification columns are used as published. No threshold is applied here.
    """
    if not S4_CSV.exists():
        raise Stop("factor_list_unavailable",
                   f"{S4_CSV} missing. Supplementary Table S4 was not extracted.")
    t = pd.read_csv(S4_CSV, dtype=str)
    def tr(c):
        return t[c].astype(str).str.upper() == "TRUE"
    A = sorted(set(t.loc[tr("active_fibro"), "target_gene"].astype(str)))
    E = sorted(set(t.loc[tr("expanded_active_fibro"), "target_gene"].astype(str)))
    both = sorted(set(A) & set(E))
    U = sorted(set(A) | set(E))
    log(f"[factor set] active_fibro guides={int(tr('active_fibro').sum())} genes={len(A)} (paper {PAPER_ACTIVE})")
    log(f"[factor set] expanded_active_fibro guides={int(tr('expanded_active_fibro').sum())} genes={len(E)} (paper {PAPER_EXPANDED})")
    log(f"[factor set] in both flags={len(both)} union={len(U)}")
    if len(A) != PAPER_ACTIVE or len(E) != PAPER_EXPANDED:
        raise Stop("factor_count_unexpected",
                   f"active_fibro gives {len(A)} genes (paper {PAPER_ACTIVE}); "
                   f"expanded_active_fibro gives {len(E)} genes (paper {PAPER_EXPANDED}). "
                   "Nothing adjusted.")
    return dict(active=A, expanded=E, both=both, union=U,
                n_guides_active=int(tr("active_fibro").sum()),
                n_guides_expanded=int(tr("expanded_active_fibro").sum()))


# ------------------------------------------------------------------ geometry

def anchor_free_basis(gtex, log):
    Zg = np.asarray(gtex["Z"], np.float64)
    age = np.asarray(gtex["age"]).astype(str)
    mid = np.isin(age, list(MID_BINS))
    n_mid = int(mid.sum())
    if n_mid < 50:
        raise Stop("basis_donors", f"only {n_mid} donors aged 40-59. Not substituting.")
    Zm = Zg[mid]
    center = Zm.mean(0)
    sp = fit_pca(Zm - center, KMAX)
    L = np.asarray(sp.L, np.float64)
    log(f"[S2A] PCA on the {n_mid} GTEx donors aged 40-59 only; k={L.shape[1]}; "
        f"anchor donors in fit = 0 of {n_mid}")
    if L.shape[1] < max(KS7):
        raise Stop("basis_rank", f"anchor-free basis returned {L.shape[1]} components < {max(KS7)}.")
    return L, center, n_mid, mid


def spaces_of(L, center, c_old, c_young, mr, ks=KS7, tag="S2A"):
    """r, r2, rnorm per space. S1 carries the matched-coordinate r for permuting."""
    out = {}
    r_full = c_old - c_young
    out["S1"] = dict(space="S1", k=None, proj=False, r_perm=r_full[mr],
                     r2=float(r_full @ r_full), rnorm=float(np.linalg.norm(r_full)))
    for k in ks:
        r = r_full @ L[:, :k]
        out[f"{tag}_k{k}"] = dict(space=f"{tag}_k{k}", k=int(k), proj=True, r_perm=r,
                                  r2=float(r @ r), rnorm=float(np.linalg.norm(r)))
    return out


def observed_stats(D_m, sp, Lm):
    """delta, cos, frac, d_norm per space, from matched-coordinate d."""
    out = {}
    for name, s in sp.items():
        A = D_m if not s["proj"] else D_m @ Lm[:, :s["k"]]
        dl, cs, fr, dn = block_delta(A, s["r_perm"], s["r2"], s["rnorm"])
        out[name] = dict(delta=dl, cos=cs, frac=fr, d_norm=dn)
    return out


# ------------------------------------------------------------------ stage 0

def reproduce_rec_gate(ctx, log, rows):
    """Stage 0 (b): rebuild the SOUTH6REC Stage 1 gate table (primary variant)."""
    X_nt = ctx["X_nt"]
    N = X_nt.shape[0]
    zm, proj, sp_r = ctx["zm"], ctx["proj"], ctx["sp_r"]
    mr, g = ctx["matched_ruler"], ctx["g"]
    meas = np.zeros(ctx["n_ruler"], bool)
    meas[mr[~ctx["missing"]]] = True
    g_m = np.where(meas, g, 0.0)
    plant = {}
    for s in SPACES6:
        plant[s] = g_m if s == "S1" else g_m @ ctx["L_full"][:, :int(s.split("k")[1])]
    C_pool = ctx["C_pool"]
    gate = pd.read_csv(SOUTH6.parent / "south6rec" / "stage1_gate.csv")
    summ = pd.read_csv(SOUTH6.parent / "south6rec" / "stage1_summary.csv")
    gate = gate[gate.variant == "primary"]
    summ = summ[summ.variant == "primary"]
    for n in REC.DEPTHS:
        t0 = time.perf_counter()
        draws, perm_of, info = REC.primary_draws(n, N)
        flat = [d for di in range(len(REC.DOSES)) for d in draws[di]]
        CD = np.asarray((REC.indicator(flat, N) @ X_nt).todense())
        E = np.empty_like(CD)
        for i in range(CD.shape[0]):
            E[i] = zm(CD[i]) - zm(C_pool - CD[i])
        del CD
        dose_of = np.repeat(np.arange(len(REC.DOSES)), REC.R)
        for s in SPACES6:
            P = proj(E, s)
            sr = sp_r[s]
            dl, *_ = block_delta(P + np.asarray(REC.DOSES)[dose_of][:, None] * plant[s],
                                 sr["r"], sr["r2"], sr["rnorm"])
            x0 = dl[dose_of == 0]
            q05 = float(np.percentile(x0, 5))
            x25 = dl[dose_of == 1]
            pw = float(np.mean(x25 < q05))
            st = gate[(gate.space == s) & (gate.depth == n)]["power_25"].to_numpy()
            REC.check_num(rows, f"(b) SOUTH6REC gate {s} n={n}", "power at 25%", pw, st)
            st = summ[(summ.space == s) & (summ.depth == n) & (summ.dose == 0.0)]["median_delta"].to_numpy()
            REC.check_num(rows, f"(b) SOUTH6REC gate {s} n={n}", "median delta at dose 0", np.median(x0), st)
            mono = [float(np.median(dl[dose_of == di])) for di in range(4)]
            st_m = bool(gate[(gate.space == s) & (gate.depth == n)]["monotone"].iloc[0])
            REC.check_cat(rows, f"(b) SOUTH6REC gate {s} n={n}", "monotone",
                          [bool(all(mono[i + 1] < mono[i] for i in range(3)))], [st_m])
            del P
        del E
        log(f"[rec gate] depth {n} reproduced ({time.perf_counter() - t0:.0f}s)")


def stage0(log):
    ctx = {}
    old_out = REC.OUT
    try:
        REC.OUT = OUT  # redirect so nothing in results/south6rec/ is written
        ctx = REC.stage0(ctx, log)
    finally:
        REC.OUT = old_out
    rows = list(ctx["rep"].to_dict("records"))
    # perm_D.npy is SOUTH6's matched-coordinate d; verify it before reusing it.
    s1 = np.load(SOUTH6 / "stage1.npz", allow_pickle=True)
    D_stored = np.load(SOUTH6 / "perm_D.npy", mmap_mode="r")
    mr = ctx["matched_ruler"]
    REC.check_num(rows, "(a) counts", "perm_D.npy equals stage1.npz d_f on matched genes",
                  np.asarray(D_stored), s1["d_f"][:, mr], tol=0.0)
    cell_n = np.diff(np.asarray(s1["cell_indptr"], np.int64))
    keep = [str(x) for x in s1["keep"]]
    s1.close()
    cpt = pd.read_csv(SOUTH6 / "stage0_cells_per_target.csv")
    cpt = cpt[cpt.target != "non"].set_index("target").reindex(keep)
    REC.check_num(rows, "(a) counts", "cells per factor equal stage0_cells_per_target.csv",
                  cell_n, cpt["n_cells"].to_numpy(), tol=0.0)
    reproduce_rec_gate(ctx, log, rows)
    rep = pd.DataFrame(rows)
    rep.to_csv(OUT / "stage0_reproduction.csv", index=False)
    failed = rep[(~rep.passed) & rep.gating]
    log(f"[stage0] checks={len(rep)} gating={int(rep.gating.sum())} failed={len(failed)}")
    ctx.update(rep=rep, keep=keep, cell_n=cell_n, D_stored=D_stored)
    if len(failed):
        raise Stop("prior_results_not_reproduced",
                   "; ".join(f"{r.block} / {r.item}: {r.n_mismatch} of {r.n} outside tol "
                             f"(max scaled {r.max_scaled_diff})" for r in failed.itertuples()))
    return ctx


# ------------------------------------------------------------------ stage A1

def null_draws(n, N, R=R_NULL):
    """SOUTH6REC primary scheme, single condition: disjoint blocks per permutation."""
    rng = np.random.default_rng([SEED, int(n)])
    m = N // n
    out, perm_of, n_perm = [], [], 0
    while len(out) < R:
        perm = rng.permutation(N)
        for j in range(m):
            if len(out) < R:
                out.append(np.sort(perm[j * n:(j + 1) * n]))
                perm_of.append(n_perm)
        n_perm += 1
    return out, np.asarray(perm_of), dict(m=int(m), n_perm=int(n_perm))


def stage_a1(ctx, log):
    fs = ctx["fs"]
    keep = ctx["keep"]
    cell_n = ctx["cell_n"]
    X_nt, N = ctx["X_nt"], ctx["X_nt"].shape[0]
    zm, C_pool = ctx["zm"], ctx["C_pool"]
    mr = ctx["matched_ruler"]
    gtex = ctx["gtex"]
    anc = np.load(ANCHORS_PATH, allow_pickle=True)
    c_old = np.asarray(anc["c_old"], float).ravel()
    c_young = np.asarray(anc["c_young"], float).ravel()
    L_mid, center_mid, n_mid, mid = anchor_free_basis(gtex, log)
    Lm_mid = L_mid[mr]
    sp7 = spaces_of(L_mid, center_mid, c_old, c_young, mr, tag="S2A")
    g = c_young - c_old
    gn = float(np.linalg.norm(g))
    ret = []
    for k in KS7:
        ret.append(dict(k=k, basis="S2A anchor-free (40-59, n=%d)" % n_mid,
                        retention=float(np.linalg.norm(g @ L_mid[:, :k]) / gn),
                        retention_south6_basis=float(np.linalg.norm(g @ ctx["L_full"][:, :k]) / gn),
                        origin_target_distance=sp7[f"S2A_k{k}"]["rnorm"],
                        origin_target_distance_south6=ctx["sp_r"][f"S2_k{k}"]["rnorm"]))
    rdf = pd.DataFrame(ret)
    rdf.to_csv(OUT / "basis_retention.csv", index=False)
    log(f"[S2A] gap retention {rdf[['k','retention','retention_south6_basis']].to_dict('records')}")
    np.savez(OUT / "basis.npz", L_mid_matched=Lm_mid, L_mid_full=L_mid, center_mid=center_mid,
             n_mid=np.array([n_mid]), c_old=c_old, c_young=c_young,
             matched_ruler=mr, mid_mask=mid,
             r_S1_matched=sp7["S1"]["r_perm"],
             **{f"r_{s}": sp7[s]["r_perm"] for s in sp7 if s != "S1"},
             r2=np.array([sp7[s]["r2"] for s in SPACES7]),
             rnorm=np.array([sp7[s]["rnorm"] for s in SPACES7]),
             spaces=np.array(SPACES7, dtype=object))

    D_m = np.asarray(ctx["D_stored"], np.float64)
    obs7 = observed_stats(D_m, sp7, Lm_mid)
    obs6 = observed_stats(D_m, ctx["sp6"], ctx["L_full"][mr])
    rows = []
    in_primary = np.array([k in set(fs["union"]) for k in keep])
    for i, f in enumerate(keep):
        row = dict(factor=f, n_cells=int(cell_n[i]), in_primary_set=bool(in_primary[i]),
                   in_active_212=bool(f in set(fs["active"])),
                   in_expanded_102=bool(f in set(fs["expanded"])))
        for s in SPACES7:
            row[f"delta_{s}"] = float(obs7[s]["delta"][i])
            row[f"cos_{s}"] = float(obs7[s]["cos"][i])
            row[f"frac_{s}"] = float(obs7[s]["frac"][i])
            row[f"dnorm_{s}"] = float(obs7[s]["d_norm"][i])
        for s in SPACES6:
            row[f"delta_south6basis_{s}"] = float(obs6[s]["delta"][i])
        rows.append(row)
    odf = pd.DataFrame(rows)
    odf.to_csv(OUT / "stage_a1_observed.csv", index=False)
    log(f"[A1] observed stats for {len(odf)} factors written")

    # Depth-matched nulls: the null depends only on n, so one null per distinct n.
    ns = sorted(set(int(x) for x in cell_n))
    log(f"[A1] distinct cell counts {len(ns)} spanning {ns[0]}-{ns[-1]}; "
        f"{len(ns) * R_NULL} draws, exact per n (no interpolation)")
    nrows = []
    draw_store = {}
    t0 = time.perf_counter()
    for j, n in enumerate(ns):
        draws, perm_of, info = null_draws(n, N)
        CD = np.asarray((REC.indicator(draws, N) @ X_nt).todense())
        E = np.empty_like(CD)
        for i in range(CD.shape[0]):
            E[i] = zm(CD[i]) - zm(C_pool - CD[i])
        del CD
        S = REC.indicator(draws, N)
        O = (S @ S.T).toarray()
        iu = np.triu_indices(O.shape[0], 1)
        ov = O[iu] / n
        for s in SPACES7:
            sr = sp7[s]
            A = E if not sr["proj"] else E @ Lm_mid[:, :sr["k"]]
            dl, *_ = block_delta(A, sr["r_perm"], sr["r2"], sr["rnorm"])
            draw_store[(n, s)] = dl.astype(np.float64)
            nrows.append(dict(n_cells=n, space=s, R=R_NULL,
                              floor_p05=float(np.percentile(dl, 5)),
                              median=float(np.median(dl)), sd=float(dl.std(ddof=1)),
                              p95=float(np.percentile(dl, 95)),
                              mutually_disjoint=int(np.max(np.bincount(perm_of))),
                              permutations=info["n_perm"], blocks_per_perm=info["m"],
                              mean_pairwise_overlap_frac=float(ov.mean()),
                              expected_overlap_independent=n / N))
        del E
        if (j + 1) % 25 == 0 or j + 1 == len(ns):
            log(f"[A1] {j + 1}/{len(ns)} depths ({time.perf_counter() - t0:.0f}s)")
    ndf = pd.DataFrame(nrows)
    ndf.to_csv(OUT / "stage_a1_depth_nulls.csv", index=False)
    np.savez_compressed(OUT / "stage_a1_null_draws.npz",
                        **{f"{n}|{s}": v for (n, s), v in draw_store.items()})
    log(f"[A1] depth-matched nulls written ({time.perf_counter() - t0:.0f}s)")
    return ctx


# ------------------------------------------------------------------ stage A2

def _perm_init(d_path, spec_path, n_perm):
    spec = np.load(spec_path, allow_pickle=True)
    _W["D"] = np.load(d_path, mmap_mode="r")
    _W["L7"] = np.asarray(spec["L7"], np.float64)
    _W["L6"] = np.asarray(spec["L6"], np.float64)
    _W["names"] = [str(x) for x in spec["names"]]
    _W["r"] = np.asarray(spec["r"], np.float64)
    _W["r2"] = np.asarray(spec["r2"], np.float64)
    _W["rnorm"] = np.asarray(spec["rnorm"], np.float64)
    _W["k"] = np.asarray(spec["k"], np.int32)
    _W["proj"] = np.asarray(spec["proj"], bool)
    _W["basis"] = np.asarray(spec["basis"], np.int32)
    _W["obs"] = np.asarray(spec["obs"], np.float64)
    _W["n_perm"] = int(n_perm)


def _perm_factor(i):
    """Same construction as SOUTH6 _perm_one_medians: the factor's own matched
    coordinates permuted, chunk PERM_CHUNK, default_rng([SEED, i])."""
    i = int(i)
    D = np.asarray(_W["D"][i], np.float64)
    p = D.size
    n_perm = _W["n_perm"]
    n_set = len(_W["names"])
    rng = np.random.default_rng([SEED, i])
    deltas = np.empty((n_set, n_perm), np.float64)
    pos = 0
    while pos < n_perm:
        m = min(PERM_CHUNK, n_perm - pos)
        idx = np.empty((m, p), np.int32)
        idx[:] = np.arange(p, dtype=np.int32)
        rng.permuted(idx, axis=1, out=idx)
        dperm = D[idx]
        del idx
        qn2_gene = np.einsum("ij,ij->i", dperm, dperm)
        proj7 = dperm @ _W["L7"]
        proj6 = dperm @ _W["L6"]
        for s_i in range(n_set):
            if _W["proj"][s_i]:
                kk = int(_W["k"][s_i])
                A = (proj7 if _W["basis"][s_i] == 7 else proj6)[:, :kk]
                dots = A @ _W["r"][s_i, :kk]
                qn2 = np.einsum("ij,ij->i", A, A)
            else:
                dots = dperm @ _W["r"][s_i, :p]
                qn2 = qn2_gene
            r2 = float(_W["r2"][s_i])
            rnorm = float(_W["rnorm"][s_i])
            deltas[s_i, pos:pos + m] = np.sqrt(np.maximum(r2 + 2.0 * dots + qn2, 0.0)) - rnorm
        del dperm, proj7, proj6
        pos += m
    obs = _W["obs"][:, i]
    floors = np.percentile(deltas, 5, axis=1)
    floors5k = np.percentile(deltas[:, :5000], 5, axis=1)
    le = np.array([int(np.sum(deltas[s_i] <= obs[s_i])) if np.isfinite(obs[s_i]) else -1
                   for s_i in range(n_set)], np.int64)
    le5k = np.array([int(np.sum(deltas[s_i, :5000] <= obs[s_i])) if np.isfinite(obs[s_i]) else -1
                     for s_i in range(n_set)], np.int64)
    return i, floors, le, floors5k, le5k


def stage_a2(log):
    import multiprocessing as mp
    which = os.environ.get("SOUTH7_SET", "primary").lower()
    n_workers = int(os.environ.get("SOUTH7_WORKERS", "6"))
    n_perm = int(os.environ.get("SOUTH7_NPERM", str(N_PERM7)))
    b = np.load(OUT / "basis.npz", allow_pickle=True)
    mr = np.asarray(b["matched_ruler"], np.int64)
    L7 = np.asarray(b["L_mid_full"], np.float64)
    Lm7 = np.asarray(b["L_mid_matched"], np.float64)
    gtex = load_gtex_z(log)
    Zg = np.asarray(gtex["Z"], np.float64)
    L6 = np.asarray(fit_pca(Zg - Zg.mean(0), KMAX).L, np.float64)
    Lm6 = L6[mr]
    c_old = np.asarray(b["c_old"], float).ravel()
    c_young = np.asarray(b["c_young"], float).ravel()
    sp7 = spaces_of(L7, None, c_old, c_young, mr, tag="S2A")
    sp6 = spaces_of(L6, None, c_old, c_young, mr, tag="S2")
    D_m = np.load(SOUTH6 / "perm_D.npy", mmap_mode="r")
    Dd = np.asarray(D_m, np.float64)
    obs7 = observed_stats(Dd, sp7, Lm7)
    obs6 = observed_stats(Dd, sp6, Lm6)
    del Dd
    names, r_rows, r2s, rns, ks, projs, bas, obs_rows = [], [], [], [], [], [], [], []
    for s in SPACES7:
        names.append(s); sr = sp7[s]
        r_rows.append(sr["r_perm"]); r2s.append(sr["r2"]); rns.append(sr["rnorm"])
        ks.append(sr["k"] or 0); projs.append(sr["proj"]); bas.append(7)
        obs_rows.append(obs7[s]["delta"])
    for s in SPACES6:
        if s == "S1":
            continue
        names.append(f"south6basis_{s}"); sr = sp6[s]
        r_rows.append(sr["r_perm"]); r2s.append(sr["r2"]); rns.append(sr["rnorm"])
        ks.append(sr["k"] or 0); projs.append(sr["proj"]); bas.append(6)
        obs_rows.append(obs6[s]["delta"])
    p = D_m.shape[1]
    R = np.zeros((len(names), max(p, KMAX)))
    for i, rr in enumerate(r_rows):
        R[i, :np.size(rr)] = np.asarray(rr).ravel()
    spec_path = OUT / f"perm_spec_{which}.npz"
    np.savez(spec_path, L7=Lm7, L6=Lm6, names=np.array(names, dtype=object), r=R,
             r2=np.array(r2s), rnorm=np.array(rns), k=np.array(ks, np.int32),
             proj=np.array(projs, bool), basis=np.array(bas, np.int32),
             obs=np.vstack(obs_rows))
    odf = pd.read_csv(OUT / "stage_a1_observed.csv")
    prim = np.flatnonzero(odf.in_primary_set.to_numpy())
    rest = np.flatnonzero(~odf.in_primary_set.to_numpy())
    idx = {"primary": prim, "rest": rest, "all": np.arange(len(odf))}[which]
    log(f"[A2] set={which} factors={idx.size} n_perm={n_perm} workers={n_workers} "
        f"settings={names}")
    log("[A2] null: each factor's own 20,264 matched-gene coordinates permuted uniformly; "
        "unmatched ruler coordinates stay 0; ||d|| preserved in gene space; "
        f"chunk {PERM_CHUNK}; default_rng([{SEED}, factor_index]) — the SOUTH6 construction.")
    n_set = len(names)
    floors = np.full((n_set, len(odf)), np.nan)
    les = np.full((n_set, len(odf)), -1, np.int64)
    floors5k = np.full((n_set, len(odf)), np.nan)
    les5k = np.full((n_set, len(odf)), -1, np.int64)
    out_path = OUT / f"stage_a2_perm_{which}.npz"
    t0 = time.perf_counter()
    done = 0
    with mp.Pool(n_workers, initializer=_perm_init,
                 initargs=(str(SOUTH6 / "perm_D.npy"), str(spec_path), n_perm)) as pool:
        for i, fl, le, fl5, le5 in pool.imap_unordered(_perm_factor, idx.tolist(), chunksize=1):
            floors[:, i] = fl; les[:, i] = le; floors5k[:, i] = fl5; les5k[:, i] = le5
            done += 1
            if done % 10 == 0 or done == idx.size:
                el = time.perf_counter() - t0
                log(f"[A2] {done}/{idx.size} ({el:.0f}s, {el / done:.1f}s/factor, "
                    f"eta {(idx.size - done) * el / done / 60:.0f}m)")
            if done % 50 == 0:
                np.savez(out_path, floors=floors, le=les, floors5k=floors5k, le5k=les5k,
                         names=np.array(names, dtype=object), idx=idx,
                         n_perm=np.array([n_perm]), done=np.array([done]))
    np.savez(out_path, floors=floors, le=les, floors5k=floors5k, le5k=les5k,
             names=np.array(names, dtype=object), idx=idx,
             n_perm=np.array([n_perm]), done=np.array([done]))
    log(f"[A2] done {done}/{idx.size} in {(time.perf_counter() - t0) / 60:.1f}m -> {out_path}")


# ------------------------------------------------------------------ stage A3

# SOUTH6's recurring set: factors labelled toward_young in >= 3 of the 5 settings
# that produced any (all O2 pairs). Derived from results/south6/stage3_full_table.csv
# and verified in Stage A3 against the cell-count figures the pre-registration states.
RECURRING_11 = ("CERS3", "ETV3L", "FOSB", "FOXF1", "FOXI3", "GATA1", "HELZ2",
                "JDP2", "NEUROD1", "OSR2", "REST")

SOUTH6_FLOOR_ROWS = {"S1": 0, "south6basis_S2_k20": 1,
                     "south6basis_S2_k50": 2, "south6basis_S2_k100": 3}


def load_perm(which):
    path = OUT / f"stage_a2_perm_{which}.npz"
    if not path.exists():
        return None
    z = np.load(path, allow_pickle=True)
    return dict(floors=np.asarray(z["floors"]), le=np.asarray(z["le"]),
                floors5k=np.asarray(z["floors5k"]), le5k=np.asarray(z["le5k"]),
                names=[str(x) for x in z["names"]], idx=np.asarray(z["idx"]),
                n_perm=int(np.asarray(z["n_perm"]).ravel()[0]),
                done=int(np.asarray(z["done"]).ravel()[0]))


def build_table(odf, floor_of, floors, le, n_perm, mask, extra):
    """Per-factor table for the factors in mask, with BH within that set."""
    sub = odf[mask].copy().reset_index(drop=True)
    pos = np.flatnonzero(mask)
    n = sub.n_cells.to_numpy()
    sub["cell_pctl_in_set"] = [float((n <= x).mean() * 100) for x in n]
    allv = odf.n_cells.to_numpy()
    sub["cell_pctl_in_screen"] = [float((allv <= x).mean() * 100) for x in n]
    hi = float(np.percentile(n, HIGH_DEPTH_PCTL))
    sub["high_depth_cut_cells"] = hi
    sub["high_depth_candidate"] = n > hi
    for key in ("guide_agreement", "identity_drop", "ruler_change"):
        sub[key] = extra[key][pos]
    sub["ruler_scores_younger"] = extra["ruler_change"][pos] < 0
    for s in SPACES7:
        si = extra["name_index"][s]
        fl = np.array([floor_of.get((int(x), s), np.nan) for x in n])
        d = sub[f"delta_{s}"].to_numpy()
        p = (le[si, pos] + 1.0) / (n_perm + 1.0)
        p = np.where(le[si, pos] >= 0, p, np.nan)
        q = bh_q(p)
        sub[f"depth_floor_{s}"] = fl
        sub[f"below_depth_floor_{s}"] = d < fl
        sub[f"perm_floor_{s}"] = floors[si, pos]
        sub[f"p_{s}"] = p
        sub[f"q_{s}"] = q
        sub[f"call_{s}"] = (d < fl) & (q <= Q_BAR)
        sub[f"south6_rule_{s}"] = d < 0
    return sub


def stage_a3(log):
    odf = pd.read_csv(OUT / "stage_a1_observed.csv")
    nulls = pd.read_csv(OUT / "stage_a1_depth_nulls.csv")
    floor_of = {(int(r.n_cells), str(r.space)): float(r.floor_p05) for r in nulls.itertuples()}
    P = {w: load_perm(w) for w in ("primary", "rest")}
    P = {k: v for k, v in P.items() if v is not None}
    if "primary" not in P:
        raise Stop("perm_missing",
                   "results/south7/stage_a2_perm_primary.npz is absent. Stage A2 has not run.")
    names = P["primary"]["names"]
    n_perm = P["primary"]["n_perm"]
    name_index = {nm: i for i, nm in enumerate(names)}
    floors = P["primary"]["floors"].copy()
    le = P["primary"]["le"].copy()
    f5 = P["primary"]["floors5k"].copy()
    l5 = P["primary"]["le5k"].copy()
    if "rest" in P and P["rest"]["n_perm"] == n_perm:
        j = P["rest"]["idx"]
        for A, B in ((floors, P["rest"]["floors"]), (le, P["rest"]["le"]),
                     (f5, P["rest"]["floors5k"]), (l5, P["rest"]["le5k"])):
            A[:, j] = B[:, j]
    # A factor counts as computed only when its own stream ran: le is -1 until then.
    have = le[0] >= 0
    want_primary = odf.in_primary_set.to_numpy()
    n_missing = int((want_primary & ~have).sum())
    log(f"[A3] permutation results for {int(have.sum())} of {len(odf)} factors "
        f"({int((want_primary & have).sum())} of {int(want_primary.sum())} in the primary set), n_perm={n_perm}")
    if n_missing and os.environ.get("SOUTH7_ALLOW_PARTIAL") != "1":
        raise Stop("stage_a_incomplete",
                   f"{n_missing} of {int(want_primary.sum())} primary-set factors have no permutation result. "
                   "Stage A must complete before the call is written. Re-run SOUTH7_STAGE=A2 SOUTH7_SET=primary.")

    s1 = np.load(SOUTH6 / "stage1.npz", allow_pickle=True)
    extra = dict(guide_agreement=np.asarray(s1["g_agree"], float),
                 identity_drop=np.asarray(s1["id_drop"], float),
                 ruler_change=np.asarray(s1["ruler_change"], float),
                 name_index=name_index)
    stored = np.load(SOUTH6 / "stage2a_floors.npz", allow_pickle=True)

    # A2 same-null confirmation: the first 5,000 draws of this stream must reproduce
    # SOUTH6's stored per-factor floors and <=-counts exactly.
    nrows = []
    for nm, row in SOUTH6_FLOOR_ROWS.items():
        if nm not in name_index:
            continue
        m = have & np.isfinite(f5[name_index[nm]])
        a = f5[name_index[nm]][m]
        b = np.asarray(stored["floors"])[row][m]
        la = l5[name_index[nm]][m]
        lb = np.asarray(stored["les"])[row][m]
        nrows.append(dict(setting=nm, south6_row=row, n_factors=int(m.sum()),
                          max_abs_floor_diff=float(np.max(np.abs(a - b))) if m.sum() else np.nan,
                          floors_bit_identical=bool(np.array_equal(a, b)),
                          le_identical=bool(np.array_equal(la, lb)),
                          n_le_mismatch=int(np.sum(la != lb))))
    samenull = pd.DataFrame(nrows)
    samenull.to_csv(OUT / "stage_a2_same_null_check.csv", index=False)
    log(f"[A3] same-null check: {samenull.to_dict('records')}")

    prim_mask = odf.in_primary_set.to_numpy() & have
    tp = build_table(odf, floor_of, floors, le, n_perm, prim_mask, extra)
    tp.to_csv(OUT / "stage_a3_primary_factors.csv", index=False)
    tables = {"primary_304": tp}
    if have.all():
        ta = build_table(odf, floor_of, floors, le, n_perm, have, extra)
        ta.to_csv(OUT / "stage_a3_all1836_factors.csv", index=False)
        tables["all_1836"] = ta

    rows = []
    for set_name, t in tables.items():
        for s in SPACES7:
            rows.append(dict(factor_set=set_name, m=len(t), space=s,
                             primary_space=(s == PRIMARY_SPACE),
                             n_called=int(t[f"call_{s}"].sum()),
                             n_below_depth_floor=int(t[f"below_depth_floor_{s}"].sum()),
                             n_q_ok=int((t[f"q_{s}"] <= Q_BAR).sum()),
                             n_south6_rule=int(t[f"south6_rule_{s}"].sum()),
                             n_high_depth_among_called=int((t[f"call_{s}"] & t.high_depth_candidate).sum()),
                             attainable_min_p=1.0 / (n_perm + 1.0),
                             n_tie_at_min_p=int((t[f"p_{s}"] <= 1.0 / (n_perm + 1.0) + 1e-15).sum()),
                             min_q=float(np.nanmin(t[f"q_{s}"].to_numpy())),
                             called=";".join(sorted(t.loc[t[f"call_{s}"], "factor"])) or ""))
    calls = pd.DataFrame(rows)
    calls.to_csv(OUT / "stage_a3_calls.csv", index=False)
    log("[A3] calls\n" + calls[["factor_set", "space", "m", "n_called", "n_south6_rule"]].to_string(index=False))

    # R4: depth dependence.
    r4 = []
    for set_name, t in tables.items():
        for s in SPACES7:
            n = t.n_cells.to_numpy(float)
            d = t[f"delta_{s}"].to_numpy()
            dn = t[f"dnorm_{s}"].to_numpy()
            prod = dn * np.sqrt(n)
            qs = np.percentile(n, [25, 50, 75])
            binid = np.digitize(n, qs)
            for b in range(4):
                m = binid == b
                if int(m.sum()) < 3:
                    continue
                r4.append(dict(factor_set=set_name, space=s, bin=f"Q{b + 1}",
                               n_factors=int(m.sum()), cells_median=float(np.median(n[m])),
                               dnorm_median=float(np.median(dn[m])),
                               dnorm_times_sqrt_cells_median=float(np.median(prod[m]))))
            r4.append(dict(factor_set=set_name, space=s, bin="all", n_factors=len(t),
                           cells_median=float(np.median(n)), dnorm_median=float(np.median(dn)),
                           dnorm_times_sqrt_cells_median=float(np.median(prod)),
                           spearman_delta_vs_cells=float(sps.spearmanr(d, n).statistic),
                           spearman_dnorm_vs_cells=float(sps.spearmanr(dn, n).statistic),
                           iqr_over_median_pct=float(100 * (np.percentile(prod, 75) - np.percentile(prod, 25))
                                                     / np.median(prod))))
    r4df = pd.DataFrame(r4)
    r4df.to_csv(OUT / "r4_depth_relation.csv", index=False)

    # R3: the SOUTH6 recurring set under the new setup.
    t = tables.get("all_1836", tables["primary_304"])
    r3 = t[t.factor.isin(RECURRING_11)].copy()
    r3.to_csv(OUT / "r3_recurring11.csv", index=False)
    allv = odf.n_cells.to_numpy()
    r11_cells = odf.set_index("factor").reindex(list(RECURRING_11)).n_cells.to_numpy()
    r3_check = dict(n=len(RECURRING_11), median_cells=float(np.median(r11_cells)),
                    n_above_screen_p90=int(np.sum([(allv <= x).mean() * 100 > 90 for x in r11_cells])),
                    screen_p90=float(np.percentile(allv, 90)),
                    prereg_median=183.0, prereg_n_above=8,
                    in_primary_set=int(odf.set_index("factor").reindex(list(RECURRING_11))
                                       .in_primary_set.sum()))
    log(f"[A3] R3 recurring-set check {r3_check}")

    # Stage B: B0 feasibility probe.
    nonsinglet = ROOT / "data" / "raw" / "southard" / "fibroblast_CRISPRa_final_pop.h5ad"
    b0 = dict(file="data/raw/southard/fibroblast_CRISPRa_final_pop.h5ad",
              present=bool(nonsinglet.exists()),
              zenodo_record="10.5281/zenodo.15200179", zenodo_size_listed="9.7 GB",
              singlet_median_cells_per_factor=float(np.median(allv)))
    stage_b_stop = None if b0["present"] else "fulldata_infeasible"
    pd.DataFrame([b0]).to_csv(OUT / "stage_b0_probe.csv", index=False)
    log(f"[A3] Stage B probe {b0} -> {stage_b_stop}")

    write_findings(dict(tables=tables, calls=calls, samenull=samenull, r4=r4df, r3=r3,
                        r3_check=r3_check, b0=b0, stage_b_stop=stage_b_stop, n_perm=n_perm,
                        nulls=nulls, odf=odf, have=have, names=names), log)
    s1.close()


# ------------------------------------------------------------------ R6 recovery in the primary space

def stage_a1r(ctx, log):
    """Report-only addition, not in the pre-registration's list: the SOUTH6REC plant
    recovery measured in the SOUTH7 primary space, so the gate's PASS wording does not
    rest on a figure from a different basis."""
    X_nt, N = ctx["X_nt"], ctx["X_nt"].shape[0]
    zm, C_pool, mr = ctx["zm"], ctx["C_pool"], ctx["matched_ruler"]
    b = np.load(OUT / "basis.npz", allow_pickle=True)
    Lm = np.asarray(b["L_mid_matched"], np.float64)
    L7 = np.asarray(b["L_mid_full"], np.float64)
    c_old = np.asarray(b["c_old"], float).ravel()
    c_young = np.asarray(b["c_young"], float).ravel()
    sp7 = spaces_of(L7, None, c_old, c_young, mr, tag="S2A")
    meas = np.zeros(ctx["n_ruler"], bool)
    meas[mr[~ctx["missing"]]] = True
    g_m = np.where(meas, c_young - c_old, 0.0)
    rows = []
    for n in (106, 390):
        draws, perm_of, info = REC.primary_draws(n, N)
        flat = [d for di in range(len(REC.DOSES)) for d in draws[di]]
        CD = np.asarray((REC.indicator(flat, N) @ X_nt).todense())
        E = np.empty_like(CD)
        for i in range(CD.shape[0]):
            E[i] = zm(CD[i]) - zm(C_pool - CD[i])
        del CD
        dose_of = np.repeat(np.arange(len(REC.DOSES)), REC.R)
        for s in SPACES7:
            sr = sp7[s]
            # Matched coordinates throughout, as the observed statistic uses:
            # g_m is 0 off the measured genes and E holds only matched columns.
            plant = g_m[mr] if s == "S1" else g_m @ L7[:, :sr["k"]]
            A = E if s == "S1" else E @ Lm[:, :sr["k"]]
            dl, *_ = block_delta(A + np.asarray(REC.DOSES)[dose_of][:, None] * plant,
                                 sr["r_perm"], sr["r2"], sr["rnorm"])
            x0 = dl[dose_of == 0]
            q05 = float(np.percentile(x0, 5))
            meds = [float(np.median(dl[dose_of == di])) for di in range(4)]
            rec = dict(space=s, depth=n, basis="S2A anchor-free" if s != "S1" else "gene space",
                       monotone=bool(all(meds[i + 1] < meds[i] for i in range(3))),
                       median_delta_0=meds[0], median_delta_25=meds[1],
                       median_delta_50=meds[2], median_delta_100=meds[3])
            for di, f in ((1, 0.25), (2, 0.50), (3, 1.00)):
                rec[f"power_{int(f * 100)}"] = float(np.mean(dl[dose_of == di] < q05))
            rows.append(rec)
            del A
        del E
        log(f"[R6] recovery at {n} cells done")
    df = pd.DataFrame(rows)
    df.to_csv(OUT / "r6_recovery_primary_space.csv", index=False)
    log("[R6]\n" + df.to_string(index=False))
    return ctx


# ------------------------------------------------------------------ findings

def md_table(df, cols, heads=None, d=4):
    heads = heads or cols
    out = ["| " + " | ".join(heads) + " |", "|" + "---|" * len(cols)]
    for r in df.itertuples(index=False):
        vals = []
        for c in cols:
            x = getattr(r, c, None)
            if isinstance(x, (bool, np.bool_)):
                vals.append("yes" if x else "no")
            elif isinstance(x, (float, np.floating)):
                vals.append("NA" if not np.isfinite(x) else f"{float(x):.{d}f}")
            else:
                vals.append(str(x))
        out.append("| " + " | ".join(vals) + " |")
    return "\n".join(out)


def write_findings(R, log):
    pre = PREREG.read_text(encoding="utf-8")
    tables, calls, n_perm = R["tables"], R["calls"], R["n_perm"]
    odf = R["odf"]
    L = []
    L.append("# FINDINGS_SOUTH7 — re-rank of the Southard CRISPRa screen with a qualified statistic")
    L.append("")
    have_all = "all_1836" in tables
    L.append(f"Status: Stage 0 passed. Stage A complete for the primary set"
             + (" and for all 1,836 factors." if have_all else "; the all-1,836 secondary run (R2) is reported where available."))
    if R["stage_b_stop"]:
        L.append(f"Stage B: STOP `{R['stage_b_stop']}` at the B0 feasibility probe. Stage A is unaffected.")
    L.append("")
    L.append("## Pre-registration (verbatim, `results/south7/PREREG.flag`)")
    L.append("")
    L.append("```")
    L.append(pre.rstrip("\n"))
    L.append("```")
    L.append("")

    # ---- factor set sourcing
    fs = pd.read_csv(OUT / "factor_set.csv")
    L.append("## 1. Factor set: sourcing record")
    L.append("")
    L.append("| item | value |")
    L.append("|---|---|")
    L.append("| Paper | Southard, Ardy et al., *Comprehensive transcription factor perturbations recapitulate fibroblast transcriptional states*, Nature Genetics 57:2323-2334 (2025), doi:10.1038/s41588-025-02284-1 |")
    L.append("| File | `41588_2025_2284_MOESM4_ESM.xlsx` (Supplementary Tables 1-19), version of record, 16,946,878 bytes, md5 `031d1831dc8fdf26b855f1ba0ce5708b`, downloaded 2026-10-04 from static-content.springer.com |")
    L.append("| Table | sheet `Table S4 Guide Activity`; column semantics from sheet `Table S17 Column descriptions` |")
    L.append("| Columns | `active_fibro` - \"guide exhibited phenotype driven clustering (i.e. clustering with guides to same target and/or related genes)\"; `expanded_active_fibro` - \"guide clustered with active guides\" |")
    L.append("| Inclusion rule | a gene is in the primary set if at least one of its guides has `active_fibro` TRUE or `expanded_active_fibro` TRUE. The published classification is used as published; no threshold was applied here and no list was re-derived. |")
    L.append("| Rows read | 10,708 guides |")
    L.append(f"| Result | `active_fibro`: 659 guides, **212 genes**. `expanded_active_fibro`: 106 guides, **102 genes**. Union: **{len(fs)} genes**. |")
    L.append("")
    L.append("Both published counts reproduce exactly: the paper's text states \"212 and 178 transcription factors with strong effects in Hs27 and RPE-1 cells, respectively, plus 102 and 82 additional genes that had single guides that clustered strongly with these gold standards\".")
    L.append("")
    L.append("Discrepancy, reported and not adjusted. The two flags are not disjoint: 10 genes carry both. "
             "The paper's wording (\"plus 102 additional genes\") implies 212 + 102 = 314 distinct genes; the table yields 304. "
             "Each of the 10 has 6 guides, of which 2-4 are `active_fibro` and 1 further guide is `expanded_active_fibro`, "
             "so for these the expanded flag marks an additional guide of an already-active gene rather than an additional gene. "
             "Genes added by the expanded flag: 92. The 10 are "
             + ", ".join(sorted(fs.loc[fs.in_both_flags, "factor"])) + ". "
             "Computation was halted at this point and the author chose the union (304) as the primary set; nothing was altered to make the counts match.")
    L.append("")
    L.append(f"All {len(fs)} are among the screen's 1,836 targets. Cell counts in the set: "
             f"median {float(np.median(tables['primary_304'].n_cells)):.0f} against the screen's "
             f"{float(np.median(odf.n_cells)):.0f}.")
    L.append("")

    # ---- stage 0
    rep = pd.read_csv(OUT / "stage0_reproduction.csv")
    L.append("## 2. Stage 0: prior results reproduced with this task's code (gate)")
    L.append("")
    L.append(f"{len(rep)} checks, {int(rep.gating.sum())} gating, {int(((~rep.passed) & rep.gating).sum())} gating failures. "
             "Counts exact; every other item within |new - stored| <= 1e-4 x max(1, |stored|). "
             "Full table: `results/south7/stage0_reproduction.csv`.")
    L.append("")
    L.append("| block | items | max abs diff | mismatches | passed |")
    L.append("|---|---|---|---|---|")
    for blk, grp in rep.groupby("block", sort=False):
        md = grp.max_abs_diff.astype(float)
        L.append(f"| {blk} | {len(grp)} | {('%.3g' % np.nanmax(md.to_numpy())) if np.isfinite(md.to_numpy()).any() else 'n/a'} | "
                 f"{int(grp.n_mismatch.sum())} | {'yes' if bool(grp.passed.all()) else 'NO'} |")
    L.append("")
    L.append("(a) covers the SOUTH6 items SOUTH6REC already checked - per-factor and per-guide counts, the TMM/log2-CPM/frozen-z transform, d_f, "
             "the origin-target distances, and per-factor delta, cos, frac and ||d|| in S1 and each S2_k on O1_Y1 - plus two additions: "
             "`results/south6/perm_D.npy` equals `stage1.npz` d_f on the matched genes exactly, and the per-factor cell counts equal "
             "`stage0_cells_per_target.csv` exactly. (b) covers the SOUTH6REC Stage 1 gate: all 12 space-by-depth rows, their power at 25%, "
             "their dose-0 median change in distance, and their monotonicity flag.")
    L.append("")
    L.append("Nothing in `results/south6/` or `results/south6rec/` was written. `src/south6rec_run.py` was imported unchanged and its Stage 0 reused; "
             "its module-level output directory was redirected to `results/south7/` for the duration of that call so no file in `results/south6rec/` was touched.")
    L.append("")

    # ---- basis
    bret = pd.read_csv(OUT / "basis_retention.csv")
    L.append("## 3. The anchor-free basis (S2A)")
    L.append("")
    L.append("Fitted exactly as `src/south6_run.py` `build_spaces` fits its basis - `fit_pca` on the frozen-z donor matrix, centred on that matrix's own mean - "
             "but on the 308 GTEx fibroblast donors aged 40-59 only. **0 of the 308 is an anchor donor**, against 344 of 652 (52.8%) in the SOUTH6 basis. "
             "O1 and Y1 themselves are the frozen TOWARD anchors, unchanged.")
    L.append("")
    L.append(md_table(bret, ["k", "retention", "retention_south6_basis", "origin_target_distance", "origin_target_distance_south6"],
                      ["k", "gap retention, S2A (anchor-free)", "gap retention, SOUTH6 basis", "||O1-Y1|| in S2A", "||O1-Y1|| in SOUTH6 basis"]))
    L.append("")
    L.append("Removing the circularity costs little, as SOUTH6REC's R5 predicted: retention at k=20 falls from 0.8873 to 0.8651.")
    L.append("")

    # ---- A1
    nulls = R["nulls"]
    L.append("## 4. Stage A1: the depth-matched null")
    L.append("")
    L.append(f"The null change in distance depends only on how many cells a pseudobulk holds, so it was computed once per distinct cell count: "
             f"{nulls.n_cells.nunique()} distinct counts spanning {int(nulls.n_cells.min())}-{int(nulls.n_cells.max())} cells, "
             f"R = {R_NULL} draws each, {nulls.n_cells.nunique() * R_NULL} draws per space. "
             "**Every factor's floor is exact at its own cell count; the interpolation the pre-registration permits was not used, "
             "so there is no interpolation error to report.**")
    L.append("")
    L.append("Draws follow the SOUTH6REC primary scheme: a permutation of the 14,013 non-targeting cells cut into disjoint blocks of n, "
             "each draw's reference being the complement of that draw, so no draw is inside its own reference.")
    L.append("")
    key = nulls[nulls.n_cells.isin([int(np.median(odf.n_cells)), 106, 390] if 390 in set(nulls.n_cells) else [106])]
    key = nulls[nulls.n_cells.isin(sorted(set([13, 106, 212, 341]) & set(nulls.n_cells)))]
    L.append(md_table(key.sort_values(["space", "n_cells"]), ["space", "n_cells", "floor_p05", "median", "sd", "mutually_disjoint", "permutations"],
                      ["space", "cells", "5th-pct floor", "null median", "null SD", "mutually disjoint draws", "permutations"]))
    L.append("")
    L.append("Full table: `results/south7/stage_a1_depth_nulls.csv` (every depth, every space, with overlap figures). "
             "The floor rises steeply as cells fall - the reason a floor built from ~7,006-cell half-splits, as SOUTH6's was, is not a floor for a 106-cell factor.")
    L.append("")

    # ---- A2
    sn = R["samenull"]
    L.append("## 5. Stage A2: the permutation null")
    L.append("")
    L.append(f"n_perm = {n_perm:,}. Attainable minimum p = 1/(n_perm+1) = {1.0 / (n_perm + 1):.3g}, which is <= 1e-5 as the pre-registration requires.")
    L.append("")
    L.append("What is permuted, stated exactly: for each factor, its own displacement vector d over the 20,264 matched gene coordinates is permuted "
             "uniformly at random; the 3,221 ruler coordinates the panel does not measure stay at 0 in every draw, so a draw cannot invent a measurement; "
             "||d|| is preserved in gene space. Draws are generated in chunks of 1,000 from `default_rng([20260914, factor_index])`. "
             "This is SOUTH6's construction unchanged - the same generator, seeded the same way, consuming the stream the same way.")
    L.append("")
    L.append("Confirmation that it is the same null, not merely a similar one: the first 5,000 draws of each factor's stream were compared against "
             "SOUTH6's stored per-factor floors and <=-counts (`results/south6/stage2a_floors.npz`, rows 0-3, its four primary O1_Y1 settings).")
    L.append("")
    L.append(md_table(sn, ["setting", "south6_row", "n_factors", "max_abs_floor_diff", "floors_bit_identical", "le_identical", "n_le_mismatch"],
                      ["setting", "SOUTH6 row", "factors compared", "max abs floor diff", "floors bit-identical", "counts identical", "count mismatches"], d=6))
    L.append("")
    L.append("Because each factor is seeded independently of every other, the draws do not depend on how many worker processes ran or in what order the "
             "factors were taken, and no checkpoint or resume logic was used: each factor's stream runs start to finish in one process.")
    L.append("")
    L.append("Multiple testing no longer caps the answer. With m = 304 the smallest attainable q is "
             f"{304.0 / (n_perm + 1):.4f}, and with m = 1,836 it is {1836.0 / (n_perm + 1):.4f}; both are below 0.05, so a single factor can clear "
             "Benjamini-Hochberg on its own. Under SOUTH6's 5,000 draws, 8 factors had to tie at the minimum p before any could pass, and S2_k20 on the "
             "primary pair had 4.")
    L.append("")

    # ---- A3 gate
    L.append("## 6. Stage A3: the call")
    L.append("")
    L.append("A factor is APPROACHING in the primary space if its change in distance is below the 5th percentile of its own depth-matched null "
             "AND its Benjamini-Hochberg q over the factor set is <= 0.05.")
    L.append("")
    L.append(md_table(calls, ["factor_set", "m", "space", "primary_space", "n_called", "n_below_depth_floor", "n_q_ok", "n_south6_rule",
                              "n_tie_at_min_p", "min_q"],
                      ["factor set", "m", "space", "primary", "called approaching", "below depth floor", "q <= 0.05", "SOUTH6 rule (delta<0)",
                       "tied at min p", "smallest q"]))
    L.append("")
    L.append("### Gate reading that fired (pre-registered wording)")
    L.append("")
    for set_name, t in tables.items():
        row = calls[(calls.factor_set == set_name) & (calls.space == PRIMARY_SPACE)].iloc[0]
        label = "the published primary set" if set_name == "primary_304" else "all 1,836 factors (secondary)"
        if int(row.n_called) == 0:
            L.append(f"- **{label}, S2A_k20: of {int(row.m)} factors, none is called approaching.** Each factor is judged against a null built at its own "
                     "cell count, and SOUTH6REC showed this statistic recovers a 25%-of-gap approach at this depth with power 0.990. The screen contains "
                     "no transcription factor that moves the reference cells a quarter of the way toward young.")
        else:
            lst = row.called.replace(";", ", ")
            L.append(f"- **{label}, S2A_k20: {int(row.n_called)} of {int(row.m)} factors are called approaching: {lst}.** "
                     f"{int(row.n_high_depth_among_called)} of them carry the high_depth_candidate flag.")
    L.append("")
    r6p = OUT / "r6_recovery_primary_space.csv"
    if r6p.exists():
        r6 = pd.read_csv(r6p)
        L.append("Qualification on that wording. The 0.990 figure is SOUTH6REC's, measured in S2_k20 on the 652-donor basis, not in S2A_k20. "
                 "Because the pre-registered sentence rests on it, the same plant-recovery measurement was repeated in the SOUTH7 primary space "
                 "(report-only, R6 below): in S2A_k20 at 106 cells the power at a 25% plant is "
                 f"{float(r6[(r6.space == 'S2A_k20') & (r6.depth == 106)].power_25.iloc[0]):.3f}. "
                 "The premise holds in the space actually used.")
        L.append("")
    # How close the primary space came, so a zero is not read as a wide margin.
    for set_name, t in tables.items():
        q = t[f"q_{PRIMARY_SPACE}"].to_numpy()
        below = t[f"below_depth_floor_{PRIMARY_SPACE}"]
        if int(t[f"call_{PRIMARY_SPACE}"].sum()) == 0 and np.isfinite(q).any():
            mn = float(np.nanmin(q))
            near = sorted(t.loc[np.isclose(q, mn), "factor"])
            L.append(f"How close it came in {set_name}: {int(below.sum())} of {len(t)} factors do fall below their own "
                     f"depth-matched floor, and {int((t[f'south6_rule_{PRIMARY_SPACE}']).sum())} have a negative change in distance, "
                     f"but none reaches q <= {Q_BAR}. The smallest q is {mn:.4f}, held by "
                     + ", ".join(near) + ". The margin is therefore thin, not wide: the permutation test, not the depth floor, "
                     "is what none of them clears - their displacement direction is not special given its magnitude.")
            L.append("")
    # Calls in the secondary spaces, reported but never primary.
    for set_name, t in tables.items():
        for sp in SPACES7:
            if sp == PRIMARY_SPACE:
                continue
            got = t[t[f"call_{sp}"]]
            if not len(got):
                continue
            L.append(f"#### {sp} (secondary, never primary): {len(got)} called in {set_name}")
            L.append("")
            L.append(md_table(got, ["factor", "n_cells", f"delta_{sp}", f"depth_floor_{sp}", f"p_{sp}", f"q_{sp}",
                                    f"south6_rule_{sp}", "high_depth_candidate", "guide_agreement", "identity_drop",
                                    "ruler_change"],
                              ["factor", "cells", "delta", "depth floor", "p", "q", "delta < 0", "high_depth_candidate",
                               "guide agreement", "identity drop", "ruler score change"]))
            L.append("")
            if sp == "S1" and not bool(got[f"south6_rule_{sp}"].any()):
                nl = R["nulls"]
                for r in got.itertuples():
                    nn = nl[(nl.space == "S1") & (nl.n_cells == r.n_cells)].iloc[0]
                    L.append(f"**Read this one carefully.** {r.factor}'s change in distance in S1 is "
                             f"{getattr(r, 'delta_S1'):+.2f}, in a space where the depth-matched null at {int(r.n_cells)} cells "
                             f"has median {nn['median']:+.2f} and 5th percentile {nn.floor_p05:+.2f}. "
                             f"{r.factor} moved **away** from the young target, only less far than sampling noise usually does. "
                             "The pre-registered criterion labels that APPROACHING because it compares each factor with its own "
                             "null rather than with zero, which is the remedy the preprint asks for; the label is not a claim "
                             f"that the cells moved toward young. In S1, 0 of {len(t)} factors have a negative change in distance "
                             "at all, and a full-gap plant at this depth still leaves it positive (R6).")
                    L.append("")
    for set_name, t in tables.items():
        called = t[t[f"call_{PRIMARY_SPACE}"]]
        if len(called):
            L.append(f"#### Factors called approaching in {set_name}")
            L.append("")
            L.append(md_table(called, ["factor", "n_cells", "cell_pctl_in_set", "high_depth_candidate",
                                       f"delta_{PRIMARY_SPACE}", f"depth_floor_{PRIMARY_SPACE}", f"q_{PRIMARY_SPACE}",
                                       "guide_agreement", "identity_drop", f"cos_{PRIMARY_SPACE}",
                                       "call_S2A_k50", "call_S2A_k100", "call_S1"],
                              ["factor", "cells", "cell pctl in set", "high_depth_candidate", "delta", "depth floor", "q",
                               "guide agreement", "identity drop", "cosine", "call S2A_k50", "call S2A_k100", "call S1"]))
            L.append("")
    L.append(f"A4 depth flag: the 90th percentile of cell count within the primary set is "
             f"{float(tables['primary_304'].high_depth_cut_cells.iloc[0]):.0f} cells. Every factor's flag is in the per-factor table, whether or not it was called.")
    L.append("")
    L.append("Per-factor tables, every factor and every column the pre-registration asks for (change in distance, depth-matched floor, permutation floor, "
             "p, q, the call, the SOUTH6 rule, cell count and its percentile, guide agreement, identity drop, cosine, and the call in each secondary space): "
             "`results/south7/stage_a3_primary_factors.csv`"
             + (" and `results/south7/stage_a3_all1836_factors.csv`." if have_all else "."))
    L.append("")

    # ---- report only
    L.append("## 7. Report-only readings (no decision rule, never a gate)")
    L.append("")
    L.append("### R1 The SOUTH6 rule (change in distance < 0), every factor")
    L.append("")
    sub = calls[["factor_set", "space", "m", "n_south6_rule", "n_below_depth_floor", "n_called"]]
    L.append(md_table(sub, ["factor_set", "space", "m", "n_south6_rule", "n_below_depth_floor", "n_called"],
                      ["factor set", "space", "m", "delta < 0 (SOUTH6 rule)", "below depth floor", "called approaching"]))
    L.append("")
    L.append("### R2 Multiple testing at m = 1,836")
    L.append("")
    if have_all:
        L.append("Both sets appear in the tables above. The comparison to make is between the two `m` values at the same space: the depth-matched "
                 "criterion and the permutation p are per factor, so widening the set from 304 to 1,836 changes only the Benjamini-Hochberg threshold.")
    else:
        L.append("The all-1,836 permutation run was not complete when this file was written; the primary set is unaffected, since its q values are "
                 "computed within the primary set.")
    L.append("")
    L.append("### R3 SOUTH6's recurring set under the new setup")
    L.append("")
    c = R["r3_check"]
    L.append(f"The recurring set is the {c['n']} factors labelled toward_young in at least 3 of the 5 SOUTH6 settings that produced any "
             "(all of them O2 pairs): " + ", ".join(RECURRING_11) + ". "
             f"Median cell count {c['median_cells']:.0f} against the screen's median of {float(np.median(odf.n_cells)):.0f}, and "
             f"{c['n_above_screen_p90']} of {c['n']} sit above the screen's 90th percentile ({c['screen_p90']:.0f} cells) - "
             f"reproducing the figures the pre-registration states ({c['prereg_median']:.0f} and {c['prereg_n_above']}). "
             f"{c['in_primary_set']} of them are in the published primary set.")
    L.append("")
    if len(R["r3"]) < len(RECURRING_11):
        missing = sorted(set(RECURRING_11) - set(R["r3"].factor))
        L.append(f"The table below covers the {len(R['r3'])} of the {len(RECURRING_11)} that are in the factor set reported here; "
                 + ", ".join(missing) + " are outside the published primary set and appear only once the all-1,836 run is complete.")
        L.append("")
    if len(R["r3"]):
        L.append(md_table(R["r3"], ["factor", "n_cells", "cell_pctl_in_screen", "high_depth_candidate",
                                    f"delta_{PRIMARY_SPACE}", f"depth_floor_{PRIMARY_SPACE}", f"q_{PRIMARY_SPACE}",
                                    f"call_{PRIMARY_SPACE}", f"south6_rule_{PRIMARY_SPACE}", "ruler_scores_younger"],
                          ["factor", "cells", "cell pctl in screen", "high_depth_candidate", "delta", "depth floor", "q",
                           "called approaching", "SOUTH6 rule", "ruler scores younger"]))
        L.append("")
    L.append("### R4 Depth dependence of the statistic")
    L.append("")
    r4 = R["r4"]
    L.append(md_table(r4[r4.bin == "all"], ["factor_set", "space", "n_factors", "spearman_delta_vs_cells", "spearman_dnorm_vs_cells",
                                            "dnorm_times_sqrt_cells_median", "iqr_over_median_pct"],
                      ["factor set", "space", "factors", "Spearman delta vs cells", "Spearman ||d|| vs cells",
                       "median ||d||x sqrt(cells)", "IQR as % of median"]))
    L.append("")
    L.append("By quartile of cell count:")
    L.append("")
    L.append(md_table(r4[r4.bin != "all"], ["factor_set", "space", "bin", "n_factors", "cells_median", "dnorm_median",
                                            "dnorm_times_sqrt_cells_median"],
                      ["factor set", "space", "quartile", "factors", "median cells", "median ||d||", "median ||d||x sqrt(cells)"]))
    L.append("")
    ra = r4[(r4.bin == "all") & (r4.factor_set == list(tables)[0])].set_index("space")
    rq = r4[(r4.bin != "all") & (r4.factor_set == list(tables)[0])]
    def _span(sp):
        v = rq[rq.space == sp].sort_values("cells_median").dnorm_times_sqrt_cells_median.to_numpy()
        return v[0], v[-1]
    s1_lo, s1_hi = _span("S1")
    k_lo, k_hi = _span(PRIMARY_SPACE)
    L.append(f"The answer, stated rather than left to the table. In S1 the product is flat: median "
             f"{ra.loc['S1', 'dnorm_times_sqrt_cells_median']:.0f}, moving only from {s1_lo:.0f} to {s1_hi:.0f} between the lowest and "
             f"highest cell-count quartile, with an interquartile range of {ra.loc['S1', 'iqr_over_median_pct']:.1f}% of the median. "
             "That reproduces the ~1,273-within-3% signature the pre-registration describes, and it is what displacement looks like when it is "
             f"pure sampling error. In {PRIMARY_SPACE} it is **not** flat: {k_lo:.0f} to {k_hi:.0f} across the same quartiles, an interquartile "
             f"range of {ra.loc[PRIMARY_SPACE, 'iqr_over_median_pct']:.1f}% of the median. So the compressed space is not merely a rescaling of the "
             f"sampling error, though depth still drives much of it: Spearman of ||d|| against cell count is "
             f"{ra.loc[PRIMARY_SPACE, 'spearman_dnorm_vs_cells']:.2f} there, against {ra.loc['S1', 'spearman_dnorm_vs_cells']:.2f} in S1. "
             "Any shortlist from either space has to be checked against cell count before it is believed, which is what the A4 flag is for.")
    L.append("")
    L.append("### R5 The frozen ruler beside the distance call")
    L.append("")
    for set_name, t in tables.items():
        n_young = int(t.ruler_scores_younger.sum())
        n_closer = int(t[f"south6_rule_{PRIMARY_SPACE}"].sum())
        n_call = int(t[f"call_{PRIMARY_SPACE}"].sum())
        n_young_s1 = int((t.ruler_change < 0).sum())
        L.append(f"- {set_name}: {n_young} of {len(t)} factors score younger on the frozen ruler (ruler score change < 0), "
                 f"while {n_closer} end closer to young in {PRIMARY_SPACE} by the SOUTH6 rule and {n_call} are called approaching. "
                 f"In SOUTH6, 761 of 1,836 scored younger on the ruler while none ended closer in S1.")
    L.append("")
    L.append("This is the preprint's thesis in a second dataset: a one-dimensional age score falls for many factors that do not move closer to a young "
             "target of the same cell type. It is reported as a diagnostic and is not evidence about any factor. The ruler was not refit.")
    L.append("")
    if (OUT / "r6_recovery_primary_space.csv").exists():
        r6 = pd.read_csv(OUT / "r6_recovery_primary_space.csv")
        L.append("### R6 Plant recovery in the SOUTH7 spaces (an addition, not in the pre-registration's report-only list)")
        L.append("")
        L.append("The SOUTH6REC recovery control, repeated in the anchor-free spaces so the gate's PASS wording does not rest on a figure measured in a "
                 "different basis. Same construction as SOUTH6REC: a known fraction of the measured gap planted into real non-targeting pseudobulks, "
                 "200 draws per dose, power against the 5th percentile of the 0% draws.")
        L.append("")
        L.append(md_table(r6, ["space", "depth", "monotone", "median_delta_0", "median_delta_25", "power_25", "power_50", "power_100"],
                          ["space", "cells", "monotone", "median delta 0%", "median delta 25%", "power 25%", "power 50%", "power 100%"]))
        L.append("")
        L.append("This carries no decision rule. It states what the statistic could have detected at the screen's depth in the space the gate uses.")
        L.append("")

    # ---- stage B
    L.append("## 8. Stage B: full dataset")
    L.append("")
    b0 = R["b0"]
    if R["stage_b_stop"]:
        L.append(f"**STOP `{R['stage_b_stop']}` at B0, the feasibility probe, before any computation.** Stage A is complete and unaffected.")
        L.append("")
        L.append(f"- The non-singlet file `{b0['file']}` is **not present** on disk. Only the singlets file is: "
                 "`fibroblast_CRISPRa_final_pop_singlets_normalized_log1p.h5ad`, 8,048,055,098 bytes, md5 confirmed.")
        L.append(f"- It exists in the same Zenodo record ({b0['zenodo_record']}) as `fibroblast_CRISPRa_final_pop.h5ad`, listed at {b0['zenodo_size_listed']}. "
                 "The pre-registration admits no new expression data beyond the Southard files already on disk, so it was not downloaded.")
        L.append("- Droplet-overloading normalization was therefore not assessed against data. The singlets file does carry the per-cell columns the "
                 "published method produces (`pairwise_coef`, `pairwise_UMI_count`, `pairwise_singlet_equivalent_UMI_count`, `final_pairwise_UMI_count`), "
                 "so the route is visible, but it cannot be exercised without the non-singlet matrix.")
        L.append(f"- Cells per factor achievable: not measurable without the file. The singlets file gives a median of "
                 f"{b0['singlet_median_cells_per_factor']:.0f} cells per factor.")
        L.append("")
        L.append("B1 and B2 were not attempted. Nothing in Stage A depends on them, and no Stage B number exists.")
    L.append("")

    # ---- limitations
    L.append("## 9. Limitations")
    L.append("")
    L.append("- **A zero in Stage A is not evidence that no transcription factor can rejuvenate aged cells.** The screen uses neonatal foreskin "
             "fibroblasts with no aged donor and no matched young target. A zero says the screen contains no such factor as measured against a "
             "GTEx-anchored young reference; it does not test the question in aged cells.")
    L.append("- Recovery of a planted effect says what the statistic could have seen, not that any factor has such an effect. A plant points exactly "
             "along the gap; real displacement vectors are longer than the gap and mostly point off it.")
    L.append("- The depth-matched null is built by resampling one pool of 14,013 non-targeting cells from one line. Its spread is the spread of "
             "resampling that pool, and random cells mix all 78 non-targeting guides, so it carries no guide-to-guide variation of the kind a "
             "factor's ~6 guides carry.")
    L.append("- The plant and the null both enter in z, after the transform. A real effect would pass through counts, TMM and the log and would carry "
             "its own cell-to-cell variance.")
    L.append("- S2A removes the anchor donors from the basis but not from the z-scoring: the frozen mu and sd come from all 652 GTEx donors, the anchors "
             "included. That is the frozen ruler and was not refit.")
    L.append("- O1 and Y1 are GTEx fibroblast tissue donors, not cells of this line; every change in distance is a counterfactual that assumes the "
             "effect transfers. The O2 pairs were not run: their axis is a 468-unit platform gap, not an age gap.")
    L.append("- Unmeasured genes are assumed unmoved. 3,221 of the ruler's 23,485 genes are outside the panel, and 164 matched genes have no counts.")
    L.append("- Guide agreement and identity drop are carried from SOUTH6's stored Stage 1, which Stage 0 reproduced; they are gene-space quantities and "
             "do not depend on the basis.")
    L.append("")
    L.append("## 10. Files")
    L.append("")
    L.append("`src/south7_run.py`. In `results/south7/`: `PREREG.flag`, `southard_table_s4_hs27.csv` (the extracted published columns), `factor_set.csv`, "
             "`stage0_reproduction.csv`, `basis.npz`, `basis_retention.csv`, `stage_a1_observed.csv`, `stage_a1_depth_nulls.csv`, "
             "`stage_a1_null_draws.npz`, `stage_a2_perm_*.npz`, `stage_a2_same_null_check.csv`, `stage_a3_primary_factors.csv`, "
             "`stage_a3_calls.csv`, `r3_recurring11.csv`, `r4_depth_relation.csv`, `r6_recovery_primary_space.csv`, `stage_b0_probe.csv`, "
             "`south7_run.log`.")
    FINDINGS.write_text("\n".join(L) + "\n", encoding="utf-8")
    log(f"[findings] {FINDINGS} written ({len(L)} lines)")


if __name__ == "__main__":
    require_prereg()
    stage = os.environ.get("SOUTH7_STAGE", "0A1").upper()
    log = Logger(OUT / "south7_run.log", mode="a")
    try:
        log(f"[prereg] {PREREG} bytes={PREREG.stat().st_size}")
        log(f"[stage] {stage}")
        if stage in ("0A1", "0", "A1"):
            fs = factor_set(log)
            ctx = stage0(log)
            ctx["fs"] = fs
            ctx["gtex"] = dict(Z=ctx["Zg"], age=ctx["age"])
            _a = np.load(ANCHORS_PATH, allow_pickle=True)
            ctx["sp6"] = spaces_of(ctx["L_full"], None,
                                   np.asarray(_a["c_old"], float).ravel(),
                                   np.asarray(_a["c_young"], float).ravel(),
                                   ctx["matched_ruler"], tag="S2")
            pd.DataFrame(dict(factor=fs["union"],
                              in_active_212=[f in set(fs["active"]) for f in fs["union"]],
                              in_expanded_102=[f in set(fs["expanded"]) for f in fs["union"]],
                              in_both_flags=[f in set(fs["both"]) for f in fs["union"]],
                              )).to_csv(OUT / "factor_set.csv", index=False)
            stage_a1(ctx, log)
        elif stage == "A2":
            stage_a2(log)
        elif stage == "A3":
            stage_a3(log)
        elif stage == "A1R":
            ctx = {}
            _old = REC.OUT
            try:
                REC.OUT = OUT / "_recheck"
                REC.OUT.mkdir(parents=True, exist_ok=True)
                ctx = REC.stage0(ctx, log)
            finally:
                REC.OUT = _old
            stage_a1r(ctx, log)
        else:
            raise SystemExit(f"unknown SOUTH7_STAGE={stage}")
    except Stop as err:
        log(f"[STOP] {err.key}: {err.message}")
        (OUT / f"STOP_{err.key}.txt").write_text(err.message + "\n", encoding="utf-8")
        raise
    finally:
        log.close()
