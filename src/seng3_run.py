"""SENG3 driver. Stages: stage0, d1, d2, d3, d4.

The compressed instrument is the frozen S2A_k20 basis, READ from
results/south7/basis.npz (PCA loadings of the 308 GTEx donors aged 40-59) and
never refit. Projection is the last step and is applied only to differences,
so the PCA centring vector cancels (amendment B6).

Reads results/south*/, results/seng/, results/seng2/ and the h5ad.
Writes only results/seng3/.
"""
from __future__ import annotations

import json
import sys
from pathlib import Path

import numpy as np
import pandas as pd

SRC = Path(__file__).resolve().parent
if str(SRC) not in sys.path:
    sys.path.insert(0, str(SRC))

import seng2_common as K  # noqa: E402
import seng_run as S  # noqa: E402
from same_run import monotonic_increasing  # noqa: E402

ROOT = SRC.parent
OUT = ROOT / "results" / "seng3"
OUT.mkdir(parents=True, exist_ok=True)
SENG2 = ROOT / "results" / "seng2"
SOUTH7 = ROOT / "results" / "south7"
SOUTH8 = ROOT / "results" / "south8"
SOUTH9 = ROOT / "results" / "south9"
KDIM = 20
DEPTHS = (96, 107, 128, 157)
NAMED_N = {"CRI_STAT3": 96, "CRA_EZH2": 107, "CRA_E2F3": 128, "CRI_ZFX": 157}


def _js(o):
    if isinstance(o, dict):
        return {str(k): _js(v) for k, v in o.items()}
    if isinstance(o, (list, tuple)):
        return [_js(v) for v in o]
    if isinstance(o, np.integer):
        return int(o)
    if isinstance(o, np.floating):
        return None if not np.isfinite(o) else float(o)
    if isinstance(o, np.bool_):
        return bool(o)
    if isinstance(o, np.ndarray):
        return _js(o.tolist())
    if isinstance(o, float) and not np.isfinite(o):
        return None
    return o


def dump(name, obj):
    (OUT / name).write_text(json.dumps(_js(obj), indent=2), encoding="utf-8")


# --------------------------------------------------------------------------
def load_basis(log):
    z = np.load(SOUTH7 / "basis.npz", allow_pickle=True)
    L = np.asarray(z["L_mid_full"], np.float64)[:, :KDIM]
    err = float(np.max(np.abs(L.T @ L - np.eye(KDIM))))
    log(f"[basis] S2A_k{KDIM} read frozen from results/south7/basis.npz "
        f"L={L.shape} n_mid={int(np.asarray(z['n_mid']).ravel()[0])} "
        f"max|LtL-I|={err:.3e}")
    if err > 1e-8:
        raise RuntimeError("basis columns are not orthonormal")
    return L, int(np.asarray(z["n_mid"]).ravel()[0]), err


def stats_block(d, r, proj):
    """SOUTH's block_delta semantics. r = origin - target. delta<0 is closer."""
    d = np.atleast_2d(np.asarray(d, np.float64))
    r = np.asarray(r, np.float64).ravel()
    if proj is not None:
        d = d @ proj
        r = r @ proj
    r2 = float(r @ r)
    rn = float(np.sqrt(r2))
    dots = d @ r
    qn2 = np.einsum("ij,ij->i", d, d)
    delta = np.sqrt(np.maximum(r2 + 2.0 * dots + qn2, 0.0)) - rn
    with np.errstate(divide="ignore", invalid="ignore"):
        frac = np.where(r2 > 1e-12, -dots / r2, np.nan)
    return delta, frac


def space3(log, need_extra=True):
    z = np.load(SENG2 / "panel_counts.npz", allow_pickle=True)
    names = [str(x) for x in z["names"]]
    ni = {n: i for i, n in enumerate(names)}
    C = z["C"]
    ctx = K.load_context(log)
    frozen = ctx["frozen"]
    Z, missing, nf, Cc = K.panel_z_from_counts(C, frozen)
    cache = S.tmm_cache_from_panel(Cc, nf)
    mu = np.asarray(frozen["mu"], float)
    sd = np.asarray(frozen["sd"], float)
    miss = np.asarray(missing, bool)
    p_unit = {}
    for m in ("CRA", "CRI"):
        pv = Z[ni[f"TOP_{m}"]] - Z[ni[f"BOT_{m}"]]
        p_unit[m] = pv / float(np.linalg.norm(pv))
    zf = lambda c: S.z_from_counts(c, cache, mu, sd, miss)[0]
    L, n_mid, berr = load_basis(log)
    ex = {}
    if need_extra:
        ex = extra3(ctx, log)
    return dict(ctx=ctx, names=names, ni=ni, C=C, Z=Z, cache=cache, mu=mu, sd=sd,
                miss=miss, zf=zf, p_unit=p_unit, L=L, n_mid=n_mid, berr=berr,
                dst=z["dst"], n_ruler=ctx["n_ruler"], frozen=frozen, ex=ex)


def extra3(ctx, log):
    path = OUT / "extra_counts.npz"
    if path.exists():
        d = np.load(path, allow_pickle=True)
        return {k: d[k] for k in d.files}
    is_wt, pdl = ctx["is_wt"], ctx["pd_labels"]
    idx = {f"{l}_FULL": np.flatnonzero(is_wt & (pdl == l)) for l in ("PD14", "PD32")}
    keys = list(idx)
    loc = K.stream_pseudobulk(ctx["f"], [idx[k] for k in keys], ctx["exp"], log)
    out = {k: K.to_ruler(loc[i], ctx["dst"], ctx["n_ruler"]) for i, k in enumerate(keys)}
    np.savez_compressed(path, **out)
    log(f"[extra] {keys} sizes={[idx[k].size for k in keys]}")
    return out


def nt_blocks(ctx, log, with_eb=True):
    out = {}
    for m in ("CRA", "CRI"):
        out[m] = np.asarray(K.submatrix(ctx["f"], ctx["nt_mod"][m], ctx["exp"]).todense(),
                            np.float64)
        log(f"[nt dense] {m} {out[m].shape}")
    if with_eb:
        eb = np.asarray(ctx["sidx"]["early_b"], int)
        out["EB"] = np.asarray(K.submatrix(ctx["f"], eb, ctx["exp"]).todense(), np.float64)
        log(f"[eb dense] {out['EB'].shape}")
    return out


# ==========================================================================
def stage0():
    log = K.Log(OUT / "stage0.log")
    checks = []

    def num(item, stored, got):
        tol = 1e-4 * max(1.0, abs(stored))
        checks.append(dict(item=item, kind="numeric", stored=stored, got=got,
                           ok=bool(np.isfinite(got) and abs(got - stored) <= tol)))

    def cnt(item, stored, got):
        checks.append(dict(item=item, kind="count", stored=stored, got=got,
                           ok=bool(got == stored)))

    # (a) SOUTH8 amended counts, recounted from the stored per-factor tables (B4)
    for fs, m_exp, below_exp, called_exp in (("primary_304", 304, 27, 24),
                                             ("all_1836", 1836, 153, 125)):
        f = pd.read_csv(SOUTH8 / f"d3c_{fs}_umi_calls.csv")
        cnt(f"SOUTH8 {fs} m", m_exp, int(len(f)))
        cnt(f"SOUTH8 {fs} below_depth_floor_S2A_k20", below_exp,
            int(f["below_depth_floor_S2A_k20"].sum()))
        cnt(f"SOUTH8 {fs} amended_call_S2A_k20", called_exp,
            int(f["amended_call_S2A_k20"].sum()))
    # (b) SOUTH8 D2 calibration, recounted from the stored per-draw p table (B3, B4)
    tn = pd.read_csv(SOUTH8 / "d2_true_null_p.csv")
    log(f"[s0b] d2_true_null_p.csv cols={list(tn.columns)[:8]} rows={len(tn)}")
    for sp, exp in (("S2A_k20", 0.035), ("S1", 0.12)):
        sub = tn[(tn.space == sp) & (tn.n_cells == 106)]
        pc = "p" if "p" in sub.columns else [c for c in sub.columns if c.startswith("p")][0]
        num(f"SOUTH8 D2 frac_p_le_05 {sp}@106", exp, float((sub[pc] <= 0.05).mean()))
    # SOUTH9's own plant numbers, reported beside as the citation check (B3)
    s9 = pd.read_csv(SOUTH9 / "d2_residualized_recovery.csv")
    s9 = s9[(s9.depth == 106) & (s9.variant == "in_space")]
    aux = {str(r.space): float(r.raw_power_25) for r in s9.itertuples()}
    log(f"[s0b] SOUTH9 raw_power_25@106 {aux}")

    # (c) SENG2 D1 cosines, RECOMPUTED from the h5ad in the full-gene space
    sp = space3(log)
    ctx, zf, ex = sp["ctx"], sp["zf"], sp["ex"]
    C, ni = sp["C"], sp["ni"]
    nc = ctx["n_counts"]
    is_wt, pdl = ctx["is_wt"], ctx["pd_labels"]
    med = {"PD14_FULL": float(np.median(nc[np.flatnonzero(is_wt & (pdl == "PD14"))])),
           "PD32_FULL": float(np.median(nc[np.flatnonzero(is_wt & (pdl == "PD32"))])),
           "NT_CRA": float(np.median(nc[ctx["nt_mod"]["CRA"]])),
           "NT_CRI": float(np.median(nc[ctx["nt_mod"]["CRI"]]))}
    cnt_v = {"PD14_FULL": ex["PD14_FULL"], "PD32_FULL": ex["PD32_FULL"],
             "NT_CRA": C[ni["NT_CRA"]], "NT_CRI": C[ni["NT_CRI"]]}
    for Lk, g, stored in (("A", "NT_CRA", 0.584155), ("I", "NT_CRI", 0.757055)):
        gl = ["PD14_FULL", "PD32_FULL", g]
        tgt, fac = K.thin_factors(med, gl)
        rng = np.random.default_rng(K.SEED)
        vals = []
        for _ in range(K.N_NULL):
            zz = {k: zf(K.thin_counts(cnt_v[k], fac[k], rng)) for k in gl}
            a = zz["PD14_FULL"] - zz["PD32_FULL"]
            t = zz["PD14_FULL"] - zz[g]
            vals.append(float(np.dot(t, a) / (np.linalg.norm(t) * np.linalg.norm(a))))
        num(f"SENG2 D1 cos(PD14-NT_{Lk}, a) depth-matched median", stored,
            float(np.median(vals)))
    # (d) SENG2 D2 counts, recounted from its stored table
    d2 = pd.read_csv(SENG2 / "d2_table.csv")
    cnt("SENG2 D2 n", 340, int(len(d2)))
    cnt("SENG2 D2 below_floor origin(i)", 131, int(d2.below_floor_delta_i.sum()))
    cnt("SENG2 D2 called origin(i)", 3, int(d2.called_delta_i.sum()))
    # (e) SENG2 D3 plant power at 100 cells, 25% dose
    ps = pd.read_csv(SENG2 / "d3_plant_summary.csv")
    for m in ("CRA", "CRI"):
        r = ps[(ps.n_cells == 100) & (ps.modality == m) & (ps.dose == 0.25)].iloc[0]
        num(f"SENG2 D3 power_25@100 {m}", 0.0, float(r.power_i))
    # (d of the method statement) the four factors' cell counts
    sizes = {lab: int(idx.size) for m, lab, idx in ctx["qual"]}
    for lab, n in NAMED_N.items():
        cnt(f"n_cells[{lab}]", n, sizes.get(lab))

    passed = all(c["ok"] for c in checks)
    pd.DataFrame(checks).to_csv(OUT / "stage0_reproduction.csv", index=False)
    dump("stage0_summary.json", dict(passed=bool(passed), checks=checks,
                                     south9_raw_power_25_at_106=aux,
                                     basis_orthonormality_err=sp["berr"],
                                     basis_n_mid=sp["n_mid"], medians=med))
    log(f"[stage0] PASS={passed}")
    if not passed:
        bad = [c["item"] for c in checks if not c["ok"]]
        log(f"[stage0] STOP prior_results_not_reproduced {bad}")
        dump("STOP.json", dict(key="prior_results_not_reproduced", mismatches=bad))
    ctx["f"].close()
    log.close()


# ==========================================================================
def _origins(sp, m):
    """r = origin - target, per SOUTH's convention, in full ruler space."""
    Z, ni = sp["Z"], sp["ni"]
    z_target = Z[ni["EARLY_A"]]
    return {"i": Z[ni[f"NT_{m}"]] - z_target, "ii": Z[ni["LATE_A"]] - z_target}


def _draw_deltas(blk, rows, n, reps, zf, dst, n_ruler, rs, L, rng, plant=None):
    """Within-set paired displacement, both origins, both spaces."""
    pool = K.to_ruler(blk[rows].sum(axis=0), dst, n_ruler)
    out = {}
    D = []
    for _ in range(reps):
        take = rng.choice(rows.size, size=n, replace=False)
        sub = K.to_ruler(blk[rows[take]].sum(axis=0), dst, n_ruler)
        d = zf(sub) - zf(pool - sub)
        if plant is not None:
            d = d + plant
        D.append(d)
    D = np.vstack(D)
    for tag in ("i", "ii"):
        for spc, proj in (("S2A_k20", L), ("S1", None)):
            dl, fr = stats_block(D, rs[tag], proj)
            out[(tag, spc)] = dict(delta=dl, frac=fr)
    return out


def d1():
    log = K.Log(OUT / "d1.log")
    sp = space3(log)
    ctx, zf, L = sp["ctx"], sp["zf"], sp["L"]
    dst, n_ruler = sp["dst"], sp["n_ruler"]
    blk = nt_blocks(ctx, log, with_eb=False)
    rows = []
    for m in ("CRA", "CRI"):
        rs = _origins(sp, m)
        nt_n = blk[m].shape[0]
        rng_split = np.random.default_rng(K.SEED)
        perm = rng_split.permutation(nt_n)
        h1, h2 = np.sort(perm[: nt_n // 2]), np.sort(perm[nt_n // 2:])
        log(f"[d1] {m} NT={nt_n} half1={h1.size} half2={h2.size}")
        for n in DEPTHS:
            rngf = np.random.default_rng(K.SEED)
            F = _draw_deltas(blk[m], h1, n, K.N_NULL, zf, dst, n_ruler, rs, L, rngf)
            rnge = np.random.default_rng(K.SEED + 1)
            E = _draw_deltas(blk[m], h2, n, K.N_NULL, zf, dst, n_ruler, rs, L, rnge)
            for tag in ("i", "ii"):
                for spc in ("S2A_k20", "S1"):
                    fl = float(np.percentile(F[(tag, spc)]["delta"], 5))
                    ed = E[(tag, spc)]["delta"]
                    rows.append(dict(
                        library=m, n_cells=n, origin=tag, space=spc, floor=fl,
                        firing_rate=float(np.mean((ed < fl) & (ed < 0))),
                        eval_median=float(np.median(ed)),
                        eval_lo=float(np.percentile(ed, 2.5)),
                        eval_hi=float(np.percentile(ed, 97.5)),
                        floor_median=float(np.median(F[(tag, spc)]["delta"])),
                        n_eval_negative=int((ed < 0).sum()),
                        n_eval_below_floor=int((ed < fl).sum())))
            log(f"[d1] {m} n={n} done")
    df = pd.DataFrame(rows)
    df.to_csv(OUT / "d1_null_validation.csv", index=False)
    prim = df[(df.space == "S2A_k20") & (df.origin == "i")]
    f_max = float(prim.firing_rate.max())
    gate = "CALIBRATED" if f_max <= 0.10 else "MISCALIBRATED"
    worst = prim.loc[prim.firing_rate.idxmax()]
    dump("d1_summary.json", dict(
        gate=gate, f_max=f_max,
        f_max_at=dict(library=str(worst.library), n_cells=int(worst.n_cells)),
        f_max_origin_ii=float(df[(df.space == "S2A_k20") & (df.origin == "ii")].firing_rate.max()),
        f_max_S1_origin_i=float(df[(df.space == "S1") & (df.origin == "i")].firing_rate.max()),
        n_combos=int(len(prim))))
    log(f"[d1] gate={gate} f_max={f_max:.4f} at {worst.library} n={int(worst.n_cells)}")
    if gate == "MISCALIBRATED":
        dump("STOP.json", dict(key="null_miscalibrated", f_max=f_max))
        log("[d1] STOP null_miscalibrated")
    ctx["f"].close()
    log.close()


# ==========================================================================
def d2():
    log = K.Log(OUT / "d2.log")
    sp = space3(log)
    ctx, zf, L, ex = sp["ctx"], sp["zf"], sp["L"], sp["ex"]
    dst, n_ruler, Z, ni = sp["dst"], sp["n_ruler"], sp["Z"], sp["ni"]
    blk = nt_blocks(ctx, log)
    # a: half A minus half A (amendment B8, primary), and the full-arm variant
    a_half = Z[ni["EARLY_A"]] - Z[ni["LATE_A"]]
    a_full = zf(ex["PD14_FULL"]) - zf(ex["PD32_FULL"])
    norms = {}
    for nm, av in (("half_A", a_half), ("full_arms", a_full)):
        norms[nm] = dict(S2A_k20=float(np.linalg.norm(av @ L)),
                         S1=float(np.linalg.norm(av)))
    a = a_half
    a_hat = {"S2A_k20": None, "S1": None}
    log(f"[d2] ||a|| {json.dumps(norms)}")

    # --- plant recovery
    prow = []
    for m in ("CRA", "CRI"):
        rs = _origins(sp, m)
        nt_n = blk[m].shape[0]
        allrows = np.arange(nt_n)
        for n in DEPTHS:
            rng = np.random.default_rng(K.SEED)
            base = _draw_deltas(blk[m], allrows, n, K.N_NULL, zf, dst, n_ruler, rs, L, rng)
            # recover the raw displacement draws once, then plant on them
            rng2 = np.random.default_rng(K.SEED)
            pool = K.to_ruler(blk[m].sum(axis=0), dst, n_ruler)
            D0 = []
            for _ in range(K.N_NULL):
                take = rng2.choice(nt_n, size=n, replace=False)
                sub = K.to_ruler(blk[m][take].sum(axis=0), dst, n_ruler)
                D0.append(zf(sub) - zf(pool - sub))
            D0 = np.vstack(D0)
            for spc, proj in (("S2A_k20", L), ("S1", None)):
                an = norms["half_A"][spc]
                ah = (a @ proj) if proj is not None else a
                ah = ah / float(np.linalg.norm(ah))
                for tag in ("i", "ii"):
                    d0p = D0 @ proj if proj is not None else D0
                    r = rs[tag] @ proj if proj is not None else rs[tag]
                    r2 = float(r @ r)
                    rn = float(np.sqrt(r2))
                    fl = None
                    for dose in (0.0, 0.25, 0.50, 1.0):
                        Dp = d0p + dose * an * ah
                        dots = Dp @ r
                        qn2 = np.einsum("ij,ij->i", Dp, Dp)
                        dl = np.sqrt(np.maximum(r2 + 2.0 * dots + qn2, 0.0)) - rn
                        if dose == 0.0:
                            fl = float(np.percentile(dl, 5))
                        prow.append(dict(library=m, n_cells=n, space=spc, origin=tag,
                                         dose=dose, floor=fl,
                                         delta_median=float(np.median(dl)),
                                         power=float(np.mean((dl < fl) & (dl < 0))),
                                         frac_median=float(np.median(-dots / r2))))
            log(f"[d2 plant] {m} n={n} done")
    pdf = pd.DataFrame(prow)
    pdf.to_csv(OUT / "d2_plant.csv", index=False)

    # --- mixture across the library boundary
    mrow = []
    for m in ("CRA", "CRI"):
        rs = _origins(sp, m)
        nt_n = blk[m].shape[0]
        eb_n = blk["EB"].shape[0]
        pool = K.to_ruler(blk[m].sum(axis=0), dst, n_ruler)
        for n_tot in (96, 157):
            for dose in (0.0, 0.10, 0.25, 0.50):
                n_e = int(round(dose * n_tot))
                n_c = n_tot - n_e
                rng = np.random.default_rng(K.SEED)
                D = []
                for _ in range(K.N_NULL):
                    tc = rng.choice(nt_n, size=n_c, replace=False)
                    cm = blk[m][tc].sum(axis=0)
                    if n_e:
                        te = rng.choice(eb_n, size=n_e, replace=False)
                        cm = cm + blk["EB"][te].sum(axis=0)
                    cmr = K.to_ruler(cm, dst, n_ruler)
                    rest = pool - K.to_ruler(blk[m][tc].sum(axis=0), dst, n_ruler)
                    D.append(zf(cmr) - zf(rest))
                D = np.vstack(D)
                for spc, proj in (("S2A_k20", L), ("S1", None)):
                    for tag in ("i", "ii"):
                        dl, fr = stats_block(D, rs[tag], proj)
                        mrow.append(dict(library=m, n_total=n_tot, space=spc, origin=tag,
                                         dose=dose,
                                         delta_median=float(np.median(dl)),
                                         delta_lo=float(np.percentile(dl, 2.5)),
                                         delta_hi=float(np.percentile(dl, 97.5)),
                                         frac_median=float(np.median(fr)),
                                         frac_lo=float(np.percentile(fr, 2.5)),
                                         frac_hi=float(np.percentile(fr, 97.5))))
            log(f"[d2 mix] {m} n_tot={n_tot} done")
    mdf = pd.DataFrame(mrow)
    mdf.to_csv(OUT / "d2_mixture.csv", index=False)

    prim = pdf[(pdf.space == "S2A_k20") & (pdf.origin == "i")]
    p25 = float(prim[prim.dose == 0.25].power.min())
    p50 = float(prim[prim.dose == 0.50].power.min())
    alt = pdf[(pdf.space == "S2A_k20") & (pdf.origin == "ii")]
    p25b = float(alt[alt.dose == 0.25].power.min())
    p50b = float(alt[alt.dose == 0.50].power.min())
    branch = lambda a25, a50: ("QUALIFIED" if a25 >= 0.80
                               else "LIMITED" if a50 >= 0.80 else "BLIND")
    gate = branch(p25, p50)
    gate_ii = branch(p25b, p50b)
    mono = {}
    for m in ("CRA", "CRI"):
        for n_tot in (96, 157):
            sub = mdf[(mdf.library == m) & (mdf.n_total == n_tot) &
                      (mdf.space == "S2A_k20") & (mdf.origin == "i")].sort_values("dose")
            mono[f"{m}_{n_tot}"] = bool(monotonic_increasing(sub.frac_median.tolist()))
    dump("d2_summary.json", dict(gate=gate, gate_origin_ii=gate_ii,
                                 origins_agree=bool(gate == gate_ii),
                                 p25=p25, p50=p50, p25_origin_ii=p25b,
                                 p50_origin_ii=p50b, a_norms=norms,
                                 mixture_frac_monotonic=mono,
                                 p25_S1_origin_i=float(
                                     pdf[(pdf.space == "S1") & (pdf.origin == "i") &
                                         (pdf.dose == 0.25)].power.min())))
    log(f"[d2] gate={gate} (origin ii: {gate_ii}) p25={p25:.3f} p50={p50:.3f}")
    if gate == "BLIND":
        dump("STOP.json", dict(key="instrument_blind", p25=p25, p50=p50))
        log("[d2] STOP instrument_blind")
    ctx["f"].close()
    log.close()


# ==========================================================================
def _proj(w, L):
    w = np.asarray(w, np.float64)
    return w if L is None else w @ L


def d3():
    log = K.Log(OUT / "d3.log")
    sp = space3(log)
    ctx, zf, L, ex = sp["ctx"], sp["zf"], sp["L"], sp["ex"]
    C, ni, Z = sp["C"], sp["ni"], sp["Z"]
    nc = ctx["n_counts"]
    is_wt, pdl = ctx["is_wt"], ctx["pd_labels"]
    sidx = ctx["sidx"]
    med = {
        "PD14_HA": float(np.median(nc[np.asarray(sidx["early_a"], int)])),
        "PD32_HA": float(np.median(nc[np.asarray(sidx["late_a"], int)])),
        "PD14_FULL": float(np.median(nc[np.flatnonzero(is_wt & (pdl == "PD14"))])),
        "PD32_FULL": float(np.median(nc[np.flatnonzero(is_wt & (pdl == "PD32"))])),
        "NT_CRA": float(np.median(nc[ctx["nt_mod"]["CRA"]])),
        "NT_CRI": float(np.median(nc[ctx["nt_mod"]["CRI"]])),
    }
    cnt_v = {"PD14_HA": C[ni["EARLY_A"]], "PD32_HA": C[ni["LATE_A"]],
             "PD14_FULL": ex["PD14_FULL"], "PD32_FULL": ex["PD32_FULL"],
             "NT_CRA": C[ni["NT_CRA"]], "NT_CRI": C[ni["NT_CRI"]]}
    variants = {"half_A": ("PD14_HA", "PD32_HA"), "full_arms": ("PD14_FULL", "PD32_FULL")}

    def quants(zz, y, o, g, L):
        a = _proj(zz[y] - zz[o], L)
        off = _proj(zz[g] - zz[o], L)
        t = _proj(zz[y] - zz[g], L)
        na = float(np.linalg.norm(a))
        cosf = lambda u, w: float(np.dot(u, w) / (np.linalg.norm(u) * np.linalg.norm(w)))
        return dict(a_norm=na, offset_norm=float(np.linalg.norm(off)),
                    offset_over_a=float(np.linalg.norm(off) / na),
                    cos_offset_a=cosf(off, a), cos_t_a=cosf(t, a),
                    t_norm=float(np.linalg.norm(t)))

    rows = []
    for vname, (y, o) in variants.items():
        for spc, Lp in (("S2A_k20", L), ("S1", None)):
            for Lk, g in (("A", "NT_CRA"), ("I", "NT_CRI")):
                gl = [y, o, g]
                raw = quants({k: zf(cnt_v[k]) for k in gl}, y, o, g, Lp)
                tgt, fac = K.thin_factors(med, gl)
                rng = np.random.default_rng(K.SEED)
                draws = []
                for _ in range(K.N_NULL):
                    zz = {k: zf(K.thin_counts(cnt_v[k], fac[k], rng)) for k in gl}
                    draws.append(quants(zz, y, o, g, Lp))
                dd = pd.DataFrame(draws)
                rec = dict(variant=vname, space=spc, library=Lk,
                           thin_target_median=tgt,
                           thin_factors=";".join(f"{k}={v:.4f}" for k, v in fac.items()))
                for c in dd.columns:
                    rec[f"raw_{c}"] = raw[c]
                    rec[f"thin_{c}"] = float(np.median(dd[c]))
                    rec[f"thin_{c}_lo"] = float(np.percentile(dd[c], 2.5))
                    rec[f"thin_{c}_hi"] = float(np.percentile(dd[c], 97.5))
                rows.append(rec)
                log(f"[d3] {vname} {spc} {Lk} cos_t_a raw={raw['cos_t_a']:.6f} "
                    f"thin={rec['thin_cos_t_a']:.6f}")
    df = pd.DataFrame(rows)
    df.to_csv(OUT / "d3_offset.csv", index=False)

    def gate_of(r):
        return "USABLE" if r >= 0.70 else "DEGRADED" if r >= 0.40 else "UNUSABLE"

    prim = df[(df.variant == "half_A") & (df.space == "S2A_k20")]
    per_lib = {str(r.library): float(r.thin_cos_t_a) for r in prim.itertuples()}
    r = max(per_lib.values())
    r_lib = max(per_lib, key=per_lib.get)
    row = prim[prim.library == r_lib].iloc[0]
    alt = df[(df.variant == "full_arms") & (df.space == "S2A_k20")]
    r_alt = float(alt.thin_cos_t_a.max())
    s1 = df[(df.variant == "half_A") & (df.space == "S1")]
    dump("d3_summary.json", dict(
        gate=gate_of(r), r=r, r_library=r_lib,
        r_ci=[float(row.thin_cos_t_a_lo), float(row.thin_cos_t_a_hi)],
        r_angle_deg=float(np.degrees(np.arccos(np.clip(r, -1, 1)))),
        pct_not_aging=100.0 * float(np.sqrt(max(0.0, 1.0 - r ** 2))),
        per_library=per_lib,
        per_library_gate={k: gate_of(v) for k, v in per_lib.items()},
        r_raw=float(prim.raw_cos_t_a.max()),
        gate_raw=gate_of(float(prim.raw_cos_t_a.max())),
        r_full_arms_variant=r_alt, gate_full_arms_variant=gate_of(r_alt),
        variants_agree=bool(gate_of(r) == gate_of(r_alt)),
        S1_per_library={str(x.library): float(x.thin_cos_t_a) for x in s1.itertuples()},
        medians=med))
    log(f"[d3] gate={gate_of(r)} r={r:.6f} ({r_lib}) per_lib={per_lib} "
        f"full_arms_r={r_alt:.6f}")
    ctx["f"].close()
    log.close()


# ==========================================================================
def d4():
    log = K.Log(OUT / "d4.log")
    sp = space3(log, need_extra=False)
    ctx, zf, L = sp["ctx"], sp["zf"], sp["L"]
    C, ni, Z = sp["C"], sp["ni"], sp["Z"]
    dst, n_ruler = sp["dst"], sp["n_ruler"]
    blk = nt_blocks(ctx, log, with_eb=False)
    pool = {m: C[ni[f"NT_{m}"]] for m in ("CRA", "CRI")}
    p20 = {m: _proj(sp["p_unit"][m], L) for m in ("CRA", "CRI")}
    p20 = {m: v / float(np.linalg.norm(v)) for m, v in p20.items()}

    def stats(d20, r20, phat):
        r2 = float(r20 @ r20)
        rn = float(np.sqrt(r2))
        dots = d20 @ r20
        qn2 = np.einsum("ij,ij->i", d20, d20)
        dl = np.sqrt(np.maximum(r2 + 2.0 * dots + qn2, 0.0)) - rn
        dr = d20 - np.outer(d20 @ phat, phat)
        rr = r20 - float(r20 @ phat) * phat
        r2r = float(rr @ rr)
        rnr = float(np.sqrt(r2r))
        dotsr = dr @ rr
        qn2r = np.einsum("ij,ij->i", dr, dr)
        dlr = np.sqrt(np.maximum(r2r + 2.0 * dotsr + qn2r, 0.0)) - rnr
        return dl, dlr, (-dots / r2 if r2 > 0 else np.nan * dots)

    by_mod = {"CRA": [], "CRI": []}
    for m, lab, idx in ctx["qual"]:
        by_mod[m].append((lab, int(idx.size)))
    floors = {}
    for m in ("CRA", "CRI"):
        rs = _origins(sp, m)
        r20 = {t: _proj(rs[t], L) for t in ("i", "ii")}
        ns = sorted({n for _, n in by_mod[m]})
        log(f"[d4] {m}: {len(ns)} distinct cell counts {ns[0]}-{ns[-1]}")
        nt_n = blk[m].shape[0]
        for j, n in enumerate(ns):
            rng = np.random.default_rng(K.SEED)
            D = []
            for _ in range(K.N_NULL):
                take = rng.choice(nt_n, size=n, replace=False)
                sub = K.to_ruler(blk[m][take].sum(axis=0), dst, n_ruler)
                D.append(zf(sub) - zf(pool[m] - sub))
            D20 = _proj(np.vstack(D), L)
            for t in ("i", "ii"):
                dl, dlr, _ = stats(D20, r20[t], p20[m])
                floors[(m, n, t)] = dict(delta=float(np.percentile(dl, 5)),
                                         delta_resid=float(np.percentile(dlr, 5)),
                                         null_median=float(np.median(dl)))
            if (j + 1) % 20 == 0:
                log(f"[d4 floors] {m} {j + 1}/{len(ns)}")
    rows = []
    for m, lab, idx in ctx["qual"]:
        n = int(idx.size)
        rs = _origins(sp, m)
        d20 = _proj((Z[ni[lab]] - Z[ni[f"NT_{m}"]])[None, :], L)
        rec = dict(modality=m, perturbation=lab, n_cells=n)
        for t in ("i", "ii"):
            r20 = _proj(rs[t], L)
            dl, dlr, fr = stats(d20, r20, p20[m])
            fl = floors[(m, n, t)]
            rec[f"delta_{t}"] = float(dl[0])
            rec[f"frac_{t}"] = float(fr[0])
            rec[f"delta_resid_{t}"] = float(dlr[0])
            for stat, val in (("delta", float(dl[0])), ("delta_resid", float(dlr[0]))):
                k = f"{stat}_{t}"
                rec[f"floor_{k}"] = fl[stat]
                rec[f"below_floor_{k}"] = bool(val < fl[stat])
                rec[f"decreased_{k}"] = bool(val < 0)
                rec[f"called_{k}"] = bool(val < fl[stat] and val < 0)
            rec[f"null_median_{t}"] = fl["null_median"]
        rows.append(rec)
    df = pd.DataFrame(rows)
    allv = df.n_cells.to_numpy()
    df["cell_pctl_in_screen"] = [float((allv <= x).mean() * 100) for x in allv]
    # D5 report-only joins
    scr = pd.read_csv(SENG2 / "stage0_screen.csv")
    df = df.merge(scr[["modality", "perturbation", "identity_diff", "identity_loss",
                       "prolif_shift", "ruler_score_context",
                       "ruler_score_nt_context", "cos", "q_cos"]],
                  on=["modality", "perturbation"], how="left")
    df["ruler_change"] = df.ruler_score_context - df.ruler_score_nt_context
    df["ruler_scores_younger"] = df.ruler_change < 0
    df.to_csv(OUT / "d4_table.csv", index=False)

    summ = {}
    for t in ("i", "ii"):
        for stat in ("delta", "delta_resid"):
            k = f"{stat}_{t}"
            kc = int(df[f"called_{k}"].sum())
            lo, hi = K.binom_ci(kc, len(df))
            summ[k] = dict(n=len(df), below_floor=int(df[f"below_floor_{k}"].sum()),
                           decreased=int(df[f"decreased_{k}"].sum()), called=kc,
                           called_ci=[lo, hi], chance_expected=0.05 * len(df))
            for m in ("CRA", "CRI"):
                s = df[df.modality == m]
                summ[k][m] = dict(n=len(s),
                                  below_floor=int(s[f"below_floor_{k}"].sum()),
                                  decreased=int(s[f"decreased_{k}"].sum()),
                                  called=int(s[f"called_{k}"].sum()),
                                  chance_expected=0.05 * len(s))
    named = {}
    for lab in K.NAMED:
        r = df[df.perturbation == lab]
        named[lab] = {} if r.empty else {c: _js(r.iloc[0][c]) for c in df.columns
                                         if c not in ("modality", "perturbation")}
    k_i = sum(1 for l in K.NAMED if named[l] and named[l]["called_delta_i"])
    k_ii = sum(1 for l in K.NAMED if named[l] and named[l]["called_delta_ii"])
    gate = ("SPLIT" if (k_i >= 2) != (k_ii >= 2)
            else "APPROACH" if (k_i >= 2 and k_ii >= 2) else "NO APPROACH")
    ruler = {}
    for t in ("i", "ii"):
        y = df.ruler_scores_younger
        c = df[f"called_delta_{t}"]
        ruler[t] = dict(n_younger=int(y.sum()), n_younger_not_called=int((y & ~c).sum()),
                        n_younger_and_called=int((y & c).sum()), n_called=int(c.sum()),
                        n_called_not_younger=int((~y & c).sum()))
        nf = df[df.perturbation.isin(K.NAMED)]
        ruler[f"{t}_named"] = dict(
            n_younger=int(nf.ruler_scores_younger.sum()),
            n_younger_not_called=int((nf.ruler_scores_younger & ~nf[f"called_delta_{t}"]).sum()))
    ident = dict(screen_median=float(df.identity_diff.median()),
                 screen_lo=float(np.percentile(df.identity_diff, 2.5)),
                 screen_hi=float(np.percentile(df.identity_diff, 97.5)),
                 n_identity_loss=int(df.identity_loss.sum()))
    dump("d4_summary.json", dict(summary=summ, named=named, k_origin_i=k_i,
                                 k_origin_ii=k_ii, gate=gate, ruler=ruler,
                                 identity=ident,
                                 rho_prolif_stored=dict(CRA=0.8280530161707566,
                                                        CRI=0.8620017521821158)))
    log(f"[d4] gate={gate} k_i={k_i} k_ii={k_ii} called_i={summ['delta_i']['called']} "
        f"called_ii={summ['delta_ii']['called']} of {len(df)}")
    ctx["f"].close()
    log.close()


STAGES = dict(stage0=stage0, d1=d1, d2=d2, d3=d3, d4=d4)

if __name__ == "__main__":
    STAGES[sys.argv[1]]()
