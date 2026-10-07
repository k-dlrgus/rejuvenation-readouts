"""SENG4 driver. Stages: prep, stage0, d1, d2, d3.

Removes SENG3's proliferation-axis circularity by reserving one NT half for
the axis and the other for the origin and every floor, and adds two axes that
touch no cell of the test at all. Reads frozen instruments; writes only
results/seng4/.
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
import seng3_run as R3  # noqa: E402
import seng_run as S  # noqa: E402

ROOT = SRC.parent
OUT = ROOT / "results" / "seng4"
OUT.mkdir(parents=True, exist_ok=True)
SENG = ROOT / "results" / "seng"
SENG2 = ROOT / "results" / "seng2"
SENG3 = ROOT / "results" / "seng3"
SOUTH9 = ROOT / "results" / "south9"
DEPTHS = (96, 107, 128, 157)
AXES = ("P1", "P2", "P3", "CIRC")


def dump(name, obj):
    (OUT / name).write_text(json.dumps(R3._js(obj), indent=2), encoding="utf-8")


# --------------------------------------------------------------------------
def cc_lists():
    t = pd.read_csv(SENG / "seurat_cc.genes_tirosh2016.tsv", sep="\t")
    return (t.loc[t.phase == "S", "symbol"].tolist(),
            t.loc[t.phase == "G2M", "symbol"].tolist())


def prep():
    """One streaming pass for library sizes; then the cells the axes need."""
    log = K.Log(OUT / "prep.log")
    ctx = K.load_context(log)
    f, exp = ctx["f"], ctx["exp"]
    X = f["X"]
    n_cells = int(X.attrs["shape"][0])
    indptr = X["indptr"][:]
    dat = X["data"]
    lib = np.zeros(n_cells, np.float64)
    for a in range(0, n_cells, K.BLOCK):
        b = min(a + K.BLOCK, n_cells)
        lo, hi = int(indptr[a]), int(indptr[b])
        if hi == lo:
            continue
        rn = np.diff(indptr[a:b + 1]).astype(np.int64)
        offs = (indptr[a:b] - lo).astype(np.int64)
        ne = rn > 0
        rs = np.zeros(b - a, np.float64)
        if ne.any():
            rs[ne] = np.add.reduceat(dat[lo:hi].astype(np.float64), offs[ne])
        lib[a:b] = rs
    log(f"[prep] lib_total over all {ctx['var'].shape[0]} file genes; "
        f"median={np.median(lib):.0f}; vs stored n_counts max abs diff="
        f"{float(np.max(np.abs(lib - ctx['n_counts']))):.0f}")

    # cells the axes and the reserved origin need
    nt = {m: np.asarray(ctx["nt_mod"][m], int) for m in ("CRA", "CRI")}
    pd26 = np.flatnonzero(ctx["is_wt"] & (ctx["pd_labels"] == "PD26"))
    rows = np.concatenate([nt["CRA"], nt["CRI"], pd26])
    sub = K.submatrix(f, rows, exp)
    loc_of = {int(r): i for i, r in enumerate(rows)}
    log(f"[prep] submatrix {sub.shape} (NT {nt['CRA'].size}+{nt['CRI'].size}, PD26 {pd26.size})")

    # cell-cycle score, SENG's recipe exactly: mean log1p(CP10k) over S genes
    # plus the same over G2M genes, CP10k using total UMI over all file genes.
    s_genes, g_genes = cc_lists()
    sym = np.array([str(x).upper() for x in np.asarray(ctx["var"].index, dtype=object)])
    first = {}
    for i, s in enumerate(sym):
        first.setdefault(s, i)
    src = np.asarray(ctx["src"], np.int64)
    file_to_local = {int(c): i for i, c in enumerate(src)}

    def cols(gl):
        out, miss_file, miss_ruler = [], [], []
        for g in gl:
            fc = first.get(g.upper())
            if fc is None:
                miss_file.append(g)
            elif int(fc) in file_to_local:
                out.append(file_to_local[int(fc)])
            else:
                miss_ruler.append(g)
        return np.array(out, int), miss_file, miss_ruler

    sc, sm_f, sm_r = cols(s_genes)
    gc, gm_f, gm_r = cols(g_genes)
    log(f"[prep] S genes used={sc.size}/{len(s_genes)} missing_file={sm_f} not_in_ruler={sm_r}")
    log(f"[prep] G2M genes used={gc.size}/{len(g_genes)} missing_file={gm_f} not_in_ruler={gm_r}")
    libsub = lib[rows]
    score = (S.module_score(sub, sc.tolist(), libsub)
             + S.module_score(sub, gc.tolist(), libsub))

    out = dict(lib=lib, rows=rows, score=score, pd26=pd26,
               nt_CRA=nt["CRA"], nt_CRI=nt["CRI"],
               s_cols=sc, g_cols=gc, dst=np.asarray(ctx["dst"], int))
    # half split of each NT pool, seed 20260914 -- identical to SENG3's split
    for m in ("CRA", "CRI"):
        rng = np.random.default_rng(K.SEED)
        perm = rng.permutation(nt[m].size)
        h1 = np.sort(perm[: nt[m].size // 2])
        h2 = np.sort(perm[nt[m].size // 2:])
        out[f"h1_{m}"] = nt[m][h1]
        out[f"h2_{m}"] = nt[m][h2]
        log(f"[prep] NT_{m} half1={h1.size} half2={h2.size} "
            f"(identical to SENG3's seed-{K.SEED} split)")

    # pseudobulks: quartile groups for P1 (half one only) and P2 (PD26), and
    # the reserved half-two pool.
    nr = ctx["n_ruler"]
    pb = {}
    qinfo = {}

    def pseudo(global_rows):
        loc = [loc_of[int(r)] for r in global_rows]
        return K.to_ruler(np.asarray(sub[loc].sum(axis=0)).ravel(), ctx["dst"], nr)

    pos_in_rows = {int(r): i for i, r in enumerate(rows)}
    for m in ("CRA", "CRI"):
        h1 = out[f"h1_{m}"]
        sc1 = score[[pos_in_rows[int(r)] for r in h1]]
        bot, top, info = S.quartile_masks(sc1)
        qinfo[f"P1_{m}"] = dict(info, n_half1=int(h1.size))
        pb[f"P1_{m}_TOP"] = pseudo(h1[top])
        pb[f"P1_{m}_BOT"] = pseudo(h1[bot])
        pb[f"NT2_{m}"] = pseudo(out[f"h2_{m}"])
        log(f"[prep] P1_{m} half1={h1.size} bottom={info['n_bottom']} top={info['n_top']} "
            f"q25={info['q25']:.4f} q75={info['q75']:.4f}")
    sc26 = score[[pos_in_rows[int(r)] for r in pd26]]
    bot, top, info = S.quartile_masks(sc26)
    qinfo["P2"] = dict(info, n_cells=int(pd26.size))
    pb["P2_TOP"] = pseudo(pd26[top])
    pb["P2_BOT"] = pseudo(pd26[bot])
    log(f"[prep] P2 PD26 n={pd26.size} bottom={info['n_bottom']} top={info['n_top']} "
        f"q25={info['q25']:.4f} q75={info['q75']:.4f}")

    np.savez_compressed(OUT / "prep.npz", **out,
                        **{f"pb_{k}": v for k, v in pb.items()})
    dump("prep_quartiles.json", qinfo)
    # cache the half-two cell blocks for drawing
    np.save(OUT / "nt2_rows_CRA.npy", out["h2_CRA"])
    np.save(OUT / "nt2_rows_CRI.npy", out["h2_CRI"])
    log("[prep] wrote prep.npz")
    ctx["f"].close()
    log.close()


# --------------------------------------------------------------------------
def stage0():
    log = K.Log(OUT / "stage0.log")
    checks = []

    def num(item, stored, got):
        tol = 1e-4 * max(1.0, abs(stored))
        checks.append(dict(item=item, kind="numeric", stored=stored, got=got,
                           ok=bool(got is not None and np.isfinite(got)
                                   and abs(got - stored) <= tol)))

    def cnt(item, stored, got):
        checks.append(dict(item=item, kind="count", stored=stored, got=got,
                           ok=bool(got == stored)))

    v1 = pd.read_csv(SENG3 / "d1_null_validation.csv")
    p = v1[(v1.space == "S2A_k20") & (v1.origin == "i")]
    cnt("SENG3 D1 eight rows", 8, int(len(p)))
    num("SENG3 D1 max firing rate", 0.0850, float(p.firing_rate.max()))
    pl = pd.read_csv(SENG3 / "d2_plant.csv")
    pp = pl[(pl.space == "S2A_k20") & (pl.origin == "i")]
    num("SENG3 D2 p25 worst", 0.565, float(pp[pp.dose == 0.25].power.min()))
    num("SENG3 D2 p50 worst", 0.940, float(pp[pp.dose == 0.50].power.min()))
    off = pd.read_csv(SENG3 / "d3_offset.csv")
    o = off[(off.variant == "half_A") & (off.space == "S2A_k20")]
    num("SENG3 D3 r", 0.887005, float(o.thin_cos_t_a.max()))
    for lib, st in (("A", 0.811375), ("I", 0.887005)):
        num(f"SENG3 D3 cos per library {lib}", st,
            float(o[o.library == lib].thin_cos_t_a.iloc[0]))
    d4 = pd.read_csv(SENG3 / "d4_table.csv")
    cnt("SENG3 D4 n", 340, int(len(d4)))
    cnt("SENG3 D4 called origin(i)", 80, int(d4.called_delta_i.sum()))
    cnt("SENG3 D4 called origin(ii)", 73, int(d4.called_delta_ii.sum()))
    nm = d4[d4.perturbation.isin(K.NAMED)]
    cnt("SENG3 D4 named called origin(i)", 4, int(nm.called_delta_i.sum()))
    cnt("SENG3 D4 named called origin(ii)", 4, int(nm.called_delta_ii.sum()))
    from scipy.stats import spearmanr
    num("SENG3 D5 rho(prolif_shift, delta_i)", -0.8639,
        float(spearmanr(d4.prolif_shift, d4.delta_i)[0]))
    cnt("SENG3 D5 residualized called origin(i)", 25, int(d4.called_delta_resid_i.sum()))
    cnt("SENG3 D5 residualized called origin(ii)", 8, int(d4.called_delta_resid_ii.sum()))
    s9 = pd.read_csv(SOUTH9 / "d2_residualized_recovery.csv")
    r9 = s9[(s9.space == "S2A_k20") & (s9.depth == 106) & (s9.variant == "in_space")].iloc[0]
    num("SOUTH9 res power_25 @106 S2A_k20", 0.990, float(r9.res_power_25))
    num("SOUTH9 raw power_25 @106 S2A_k20", 0.980, float(r9.raw_power_25))
    # (d) the four cell counts, re-derived from the h5ad
    sp = R3.space3(log, need_extra=False)
    sizes = {lab: int(idx.size) for m, lab, idx in sp["ctx"]["qual"]}
    for lab, n in (("CRA_E2F3", 128), ("CRA_EZH2", 107),
                   ("CRI_STAT3", 96), ("CRI_ZFX", 157)):
        cnt(f"n_cells[{lab}]", n, sizes.get(lab))
    audit = json.loads((SENG2 / "raw_audit.json").read_text())
    cnt("raw integer counts (no_raw_counts)", 0, int(audit["n_fractional"]))

    passed = all(c["ok"] for c in checks)
    pd.DataFrame(checks).to_csv(OUT / "stage0_reproduction.csv", index=False)
    p.to_csv(OUT / "stage0_seng3_d1_eight_rows.csv", index=False)
    dump("stage0_summary.json", dict(passed=bool(passed), checks=checks))
    log(f"[stage0] PASS={passed}")
    if not passed:
        bad = [c["item"] for c in checks if not c["ok"]]
        dump("STOP.json", dict(key="prior_results_not_reproduced", mismatches=bad))
        log(f"[stage0] STOP {bad}")
    sp["ctx"]["f"].close()
    log.close()


# --------------------------------------------------------------------------
def _axes(sp, pr, log):
    """The four axes, unit vectors in S2A_k20."""
    L = sp["L"]
    zf, ni, Z = sp["zf"], sp["ni"], sp["Z"]
    ax = {}
    for m in ("CRA", "CRI"):
        pv = zf(pr[f"pb_P1_{m}_TOP"]) - zf(pr[f"pb_P1_{m}_BOT"])
        ax[("P1", m)] = pv
        ax[("CIRC", m)] = Z[ni[f"TOP_{m}"]] - Z[ni[f"BOT_{m}"]]
    p2 = zf(pr["pb_P2_TOP"]) - zf(pr["pb_P2_BOT"])
    e = np.zeros(sp["n_ruler"], np.float64)
    dst = np.asarray(pr["dst"], int)
    cc_local = np.concatenate([np.asarray(pr["s_cols"], int), np.asarray(pr["g_cols"], int)])
    e[dst[cc_local]] = 1.0
    e /= float(np.linalg.norm(e))
    for m in ("CRA", "CRI"):
        ax[("P2", m)] = p2
        ax[("P3", m)] = e
    out = {}
    for k, v in ax.items():
        p20 = v @ L
        n = float(np.linalg.norm(p20))
        out[k] = p20 / n if n > 1e-12 else p20
    log(f"[axes] P3 built on {cc_local.size} measured ruler cc genes")
    return out, cc_local.size


def d1():
    log = K.Log(OUT / "d1.log")
    sp = R3.space3(log)
    pr = np.load(OUT / "prep.npz", allow_pickle=True)
    L, zf, ni, Z = sp["L"], sp["zf"], sp["ni"], sp["Z"]
    axes, n_cc = _axes(sp, pr, log)
    a20 = (Z[ni["EARLY_A"]] - Z[ni["LATE_A"]]) @ L
    a_full = Z[ni["EARLY_A"]] - Z[ni["LATE_A"]]
    cosf = lambda x, y: float(np.dot(x, y) / (np.linalg.norm(x) * np.linalg.norm(y)))

    rows = []
    for m in ("CRA", "CRI"):
        for i, A in enumerate(AXES):
            for B in AXES[i:]:
                rows.append(dict(modality=m, axis_a=A, axis_b=B,
                                 cos=cosf(axes[(A, m)], axes[(B, m)])))
    pd.DataFrame(rows).to_csv(OUT / "d1_pairwise.csv", index=False)

    n_in = {"P1_CRA": int(pr["h1_CRA"].size), "P1_CRI": int(pr["h1_CRI"].size),
            "P2": int(pr["pd26"].size), "P3": 0, "CIRC_CRA": 666, "CIRC_CRI": 1479}
    arows = []
    for A in AXES:
        for m in ("CRA", "CRI"):
            c = cosf(axes[(A, m)], a20)
            arows.append(dict(axis=A, modality=m, cos_with_a=c, abs_cos=abs(c),
                              angle_deg=float(np.degrees(np.arccos(np.clip(c, -1, 1)))),
                              frac_a_surviving=float(np.sqrt(max(0.0, 1 - c ** 2))),
                              pct_a_removed=100.0 * (1 - float(np.sqrt(max(0.0, 1 - c ** 2)))),
                              n_cells_in_axis=n_in.get(f"{A}_{m}", n_in.get(A, 0))))
    adf = pd.DataFrame(arows)
    adf.to_csv(OUT / "d1_vs_aging.csv", index=False)

    # a against the GTEx donor-age direction u
    u = np.asarray(np.load(S.ANCHORS_PATH, allow_pickle=True)["u"], float).ravel()
    cos_au_full = cosf(a_full, u)
    cos_au_k20 = cosf(a20, u @ L)
    south = {}
    for p, col in (("south9/d3_summary.csv", None),):
        pass
    stored_u = None
    try:
        b = np.load(ROOT / "results" / "south7" / "basis.npz", allow_pickle=True)
        r_old_young = np.asarray(b["c_old"], float) - np.asarray(b["c_young"], float)
        south = dict(southard_axis_vs_u_full=cosf(r_old_young, u),
                     southard_axis_vs_u_k20=cosf(r_old_young @ L, u @ L),
                     note="Southard's own old-minus-young axis is built from GTEx "
                          "donor age bins, so this is a consistency check, not an "
                          "independent comparison")
    except Exception as e:  # pragma: no cover
        south = dict(error=str(e))

    p1abs = float(adf[(adf.axis == "P1")].abs_cos.max())
    p1row = adf[(adf.axis == "P1")].loc[adf[(adf.axis == "P1")].abs_cos.idxmax()]
    gate = ("SEPARABLE" if p1abs < 0.40 else
            "ENTANGLED" if p1abs < 0.70 else "INSEPARABLE")
    dump("d1_summary.json", dict(
        gate=gate, p1_abs_cos=p1abs, p1_modality=str(p1row.modality),
        p1_cos=float(p1row.cos_with_a), p1_angle_deg=float(p1row.angle_deg),
        p1_frac_a_surviving=float(p1row.frac_a_surviving),
        p1_pct_a_removed=float(p1row.pct_a_removed),
        per_axis={f"{r.axis}_{r.modality}": dict(
            cos=float(r.cos_with_a), angle_deg=float(r.angle_deg),
            frac_surviving=float(r.frac_a_surviving)) for r in adf.itertuples()},
        cos_a_u_full_gene=cos_au_full, cos_a_u_S2A_k20=cos_au_k20,
        southard=south, n_cc_genes_in_P3=n_cc, n_cells_in_axis=n_in,
        a_norm_k20=float(np.linalg.norm(a20))))
    log(f"[d1] gate={gate} |cos(P1,a)|={p1abs:.6f} ({p1row.modality}) "
        f"angle={float(p1row.angle_deg):.2f} cos(a,u)full={cos_au_full:.6f}")
    sp["ctx"]["f"].close()
    log.close()


# --------------------------------------------------------------------------
def _resid_stats(d20, r20, phat):
    """delta under a residualization. phat=None means no residualization.

    The displacement, the origin-to-target vector and therefore the null are
    all put through the same projection, exactly as SENG3 does.
    """
    if phat is not None:
        d20 = d20 - np.outer(d20 @ phat, phat)
        r20 = r20 - float(r20 @ phat) * phat
    r2 = float(r20 @ r20)
    rn = float(np.sqrt(r2))
    dots = d20 @ r20
    qn2 = np.einsum("ij,ij->i", d20, d20)
    return np.sqrt(np.maximum(r2 + 2.0 * dots + qn2, 0.0)) - rn


def _setup(log):
    sp = R3.space3(log)
    pr = np.load(OUT / "prep.npz", allow_pickle=True)
    axes, n_cc = _axes(sp, pr, log)
    ctx, zf, L, ni, Z = sp["ctx"], sp["zf"], sp["L"], sp["ni"], sp["Z"]
    dst, nr = np.asarray(pr["dst"], int), sp["n_ruler"]
    # reserved half two: cell block, pool pseudobulk, origin
    h2 = {m: np.asarray(pr[f"h2_{m}"], int) for m in ("CRA", "CRI")}
    blk = {m: np.asarray(K.submatrix(ctx["f"], h2[m], ctx["exp"]).todense(), np.float64)
           for m in ("CRA", "CRI")}
    pool = {m: K.to_ruler(blk[m].sum(axis=0), dst, nr) for m in ("CRA", "CRI")}
    z_nt2 = {m: zf(pool[m]) for m in ("CRA", "CRI")}
    for m in ("CRA", "CRI"):
        log(f"[setup] reserved NT half two {m}: {h2[m].size} cells")
    z_target, z_late = Z[ni["EARLY_A"]], Z[ni["LATE_A"]]
    r20 = {}
    for m in ("CRA", "CRI"):
        r20[(m, "i")] = (z_nt2[m] - z_target) @ L
        r20[(m, "ii")] = (z_late - z_target) @ L
    a20 = (z_target - z_late) @ L
    return dict(sp=sp, pr=pr, axes=axes, blk=blk, pool=pool, dst=dst, nr=nr,
                r20=r20, a20=a20, zf=zf, L=L, h2=h2, n_cc=n_cc)


def _null_draws(E, n, reps, seed):
    """Displacements from the reserved half: z(draw) - z(half minus draw)."""
    blk, pool, zf, dst, nr = E["blk"], E["pool"], E["zf"], E["dst"], E["nr"]
    out = {}
    for m in ("CRA", "CRI"):
        rng = np.random.default_rng(seed)
        D = []
        for _ in range(reps):
            take = rng.choice(blk[m].shape[0], size=n, replace=False)
            sub = K.to_ruler(blk[m][take].sum(axis=0), dst, nr)
            D.append(zf(sub) - zf(pool[m] - sub))
        out[m] = np.vstack(D) @ E["L"]
    return out


def d2():
    log = K.Log(OUT / "d2.log")
    E = _setup(log)
    a20 = E["a20"]
    an = float(np.linalg.norm(a20))
    ahat = a20 / an
    seng3 = pd.read_csv(SENG3 / "d2_plant.csv")
    rows = []
    for n in DEPTHS:
        D = _null_draws(E, n, K.N_NULL, K.SEED)
        for m in ("CRA", "CRI"):
            for ax in ("NONE", "P1", "P2", "P3"):
                ph = None if ax == "NONE" else E["axes"][(ax, m)]
                for tag in ("i", "ii"):
                    r = E["r20"][(m, tag)]
                    fl = None
                    for dose in (0.0, 0.25, 0.50, 1.0):
                        dl = _resid_stats(D[m] + dose * an * ahat, r, ph)
                        if dose == 0.0:
                            fl = float(np.percentile(dl, 5))
                        s3 = seng3[(seng3.space == "S2A_k20") & (seng3.origin == tag)
                                   & (seng3.library == m) & (seng3.n_cells == n)
                                   & (seng3.dose == dose)]
                        rows.append(dict(
                            axis=ax, library=m, n_cells=n, origin=tag, dose=dose,
                            floor=fl, delta_median=float(np.median(dl)),
                            power=float(np.mean((dl < fl) & (dl < 0))),
                            seng3_power_unresidualized=(float(s3.power.iloc[0])
                                                        if len(s3) else np.nan)))
        log(f"[d2] n={n} done")
    df = pd.DataFrame(rows)
    df.to_csv(OUT / "d2_plant_residualized.csv", index=False)
    g = {}
    for ax in ("NONE", "P1", "P2", "P3"):
        for tag in ("i", "ii"):
            sub = df[(df.axis == ax) & (df.origin == tag)]
            g[f"{ax}_{tag}"] = dict(
                p25=float(sub[sub.dose == 0.25].power.min()),
                p50=float(sub[sub.dose == 0.50].power.min()),
                p100=float(sub[sub.dose == 1.0].power.min()))
    p50 = g["P1_i"]["p50"]
    p25 = g["P1_i"]["p25"]
    gate = "RETAINED" if p50 >= 0.80 else "LOST"
    dump("d2_summary.json", dict(gate=gate, p50=p50, p25=p25, by_axis_origin=g,
                                 seng3_p50=0.940, seng3_p25=0.565, a_norm_k20=an))
    log(f"[d2] gate={gate} P1 p50={p50:.3f} p25={p25:.3f} | "
        f"NONE p50={g['NONE_i']['p50']:.3f} P2 p50={g['P2_i']['p50']:.3f} "
        f"P3 p50={g['P3_i']['p50']:.3f}")
    E["sp"]["ctx"]["f"].close()
    log.close()


def d3():
    log = K.Log(OUT / "d3.log")
    E = _setup(log)
    sp, L = E["sp"], E["L"]
    ctx, zf, ni, Z = sp["ctx"], sp["zf"], sp["ni"], sp["Z"]
    by_mod = {"CRA": [], "CRI": []}
    for m, lab, idx in ctx["qual"]:
        by_mod[m].append((lab, int(idx.size)))
    floors = {}
    for m in ("CRA", "CRI"):
        ns = sorted({n for _, n in by_mod[m]})
        log(f"[d3] {m}: {len(ns)} distinct cell counts")
        for j, n in enumerate(ns):
            D = _null_draws(E, n, K.N_NULL, K.SEED)[m]
            for ax in ("NONE", "P1", "P2", "P3"):
                ph = None if ax == "NONE" else E["axes"][(ax, m)]
                for tag in ("i", "ii"):
                    dl = _resid_stats(D, E["r20"][(m, tag)], ph)
                    floors[(m, n, ax, tag)] = float(np.percentile(dl, 5))
            if (j + 1) % 20 == 0:
                log(f"[d3 floors] {m} {j + 1}/{len(ns)}")
    rows = []
    for m, lab, idx in ctx["qual"]:
        n = int(idx.size)
        d_full = Z[ni[lab]] - zf(E["pool"][m])
        d20 = (d_full @ L)[None, :]
        rec = dict(modality=m, perturbation=lab, n_cells=n)
        for ax in ("NONE", "P1", "P2", "P3"):
            ph = None if ax == "NONE" else E["axes"][(ax, m)]
            for tag in ("i", "ii"):
                dl = float(_resid_stats(d20, E["r20"][(m, tag)], ph)[0])
                fl = floors[(m, n, ax, tag)]
                rec[f"delta_{ax}_{tag}"] = dl
                rec[f"floor_{ax}_{tag}"] = fl
                rec[f"below_{ax}_{tag}"] = bool(dl < fl)
                rec[f"neg_{ax}_{tag}"] = bool(dl < 0)
                rec[f"called_{ax}_{tag}"] = bool(dl < fl and dl < 0)
        rows.append(rec)
    df = pd.DataFrame(rows)
    allv = df.n_cells.to_numpy()
    df["cell_pctl_in_screen"] = [float((allv <= x).mean() * 100) for x in allv]
    scr = pd.read_csv(SENG2 / "stage0_screen.csv")
    df = df.merge(scr[["modality", "perturbation", "identity_diff", "identity_loss",
                       "prolif_shift", "ruler_score_context",
                       "ruler_score_nt_context"]],
                  on=["modality", "perturbation"], how="left")
    df["ruler_change"] = df.ruler_score_context - df.ruler_score_nt_context
    df["ruler_scores_younger"] = df.ruler_change < 0
    s3 = pd.read_csv(SENG3 / "d4_table.csv")[
        ["modality", "perturbation", "delta_i", "delta_ii",
         "called_delta_i", "called_delta_ii"]].rename(
        columns={"delta_i": "seng3_delta_i", "delta_ii": "seng3_delta_ii",
                 "called_delta_i": "seng3_called_i",
                 "called_delta_ii": "seng3_called_ii"})
    df = df.merge(s3, on=["modality", "perturbation"], how="left")
    df.to_csv(OUT / "d3_table.csv", index=False)

    from scipy.stats import spearmanr
    summ, rho = {}, {}
    for ax in ("NONE", "P1", "P2", "P3"):
        for tag in ("i", "ii"):
            kc = int(df[f"called_{ax}_{tag}"].sum())
            lo, hi = K.binom_ci(kc, len(df))
            summ[f"{ax}_{tag}"] = dict(
                n=len(df), below=int(df[f"below_{ax}_{tag}"].sum()),
                neg=int(df[f"neg_{ax}_{tag}"].sum()), called=kc,
                called_ci=[lo, hi], chance_expected=0.05 * len(df),
                CRA=int(df[(df.modality == "CRA")][f"called_{ax}_{tag}"].sum()),
                CRI=int(df[(df.modality == "CRI")][f"called_{ax}_{tag}"].sum()))
        rho[ax] = float(spearmanr(df.prolif_shift, df[f"delta_{ax}_i"])[0])
    named = {}
    band = df[(df.n_cells >= 96) & (df.n_cells <= 157)]
    for lab in K.NAMED:
        r = df[df.perturbation == lab]
        if r.empty:
            continue
        r = r.iloc[0]
        e = dict(n_cells=int(r.n_cells), pctl=float(r.cell_pctl_in_screen),
                 prolif_shift=float(r.prolif_shift),
                 identity_diff=float(r.identity_diff),
                 ruler_scores_younger=bool(r.ruler_scores_younger),
                 seng3_delta_i=float(r.seng3_delta_i),
                 seng3_called_i=bool(r.seng3_called_i),
                 seng3_called_ii=bool(r.seng3_called_ii))
        for ax in ("NONE", "P1", "P2", "P3"):
            for tag in ("i", "ii"):
                e[f"delta_{ax}_{tag}"] = float(r[f"delta_{ax}_{tag}"])
                e[f"floor_{ax}_{tag}"] = float(r[f"floor_{ax}_{tag}"])
                e[f"called_{ax}_{tag}"] = bool(r[f"called_{ax}_{tag}"])
            e[f"band_pctl_{ax}_i"] = float(
                (band[f"delta_{ax}_i"] <= r[f"delta_{ax}_i"]).mean() * 100)
        named[lab] = e
    ks = {}
    for ax in ("NONE", "P1", "P2", "P3"):
        ks[ax] = {t: sum(1 for l in K.NAMED if named[l][f"called_{ax}_{t}"])
                  for t in ("i", "ii")}
    br = lambda k: ("SPLIT" if (k["i"] >= 2) != (k["ii"] >= 2)
                    else "SURVIVES" if k["i"] >= 2 else "EXPLAINED")
    gate = br(ks["P1"])
    # decile table, residualized on P1
    df["dec"] = pd.qcut(df.n_cells, 10, labels=False, duplicates="drop")
    dec = df.groupby("dec").agg(
        n=("n_cells", "size"), median_cells=("n_cells", "median"),
        called_P1_i=("called_P1_i", "sum"), called_NONE_i=("called_NONE_i", "sum"),
        median_delta_P1_i=("delta_P1_i", "median")).reset_index()
    dec["rate_P1_i"] = dec.called_P1_i / dec.n
    dec["rate_NONE_i"] = dec.called_NONE_i / dec.n
    dec.to_csv(OUT / "d3_deciles.csv", index=False)
    ruler = {}
    for ax in ("NONE", "P1"):
        for tag in ("i", "ii"):
            y, c = df.ruler_scores_younger, df[f"called_{ax}_{tag}"]
            ruler[f"{ax}_{tag}"] = dict(n_younger=int(y.sum()),
                                        n_younger_not_called=int((y & ~c).sum()),
                                        n_called=int(c.sum()))
    dump("d3_summary.json", dict(
        gate=gate, gate_by_axis={a: br(ks[a]) for a in ("P1", "P2", "P3", "NONE")},
        k_named=ks, summary=summ, rho_prolif_vs_delta=rho, named=named,
        band_n=int(len(band)),
        band_called_P1_i=int(band.called_P1_i.sum()),
        deciles=dec.round(4).to_dict("records"), ruler=ruler,
        identity=dict(screen_median=float(df.identity_diff.median()),
                      n_identity_loss=int(df.identity_loss.sum()))))
    log(f"[d3] gate={gate} k_named={ks} called P1_i={summ['P1_i']['called']} "
        f"NONE_i={summ['NONE_i']['called']} of {len(df)}")
    E["sp"]["ctx"]["f"].close()
    log.close()


STAGES = dict(prep=prep, stage0=stage0, d1=d1, d2=d2, d3=d3)

if __name__ == "__main__":
    STAGES[sys.argv[1]]()
