"""AXANG driver. Stages: stage0, d1, d2, d3.

Descriptive geometry: the angle between every aging axis this project uses and
a proliferation axis, by two constructions, in two spaces. Frozen instruments
are read, not refit. Writes only results/axang/.
"""
from __future__ import annotations

import json
import sys
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import sparse

SRC = Path(__file__).resolve().parent
if str(SRC) not in sys.path:
    sys.path.insert(0, str(SRC))

import seng2_common as K  # noqa: E402
import seng3_run as R3  # noqa: E402
import seng_run as S  # noqa: E402
from fibro2_common import load_frozen_ruler  # noqa: E402
from target_common import unit  # noqa: E402
from toward_run import load_gtex_z  # noqa: E402

ROOT = SRC.parent
OUT = ROOT / "results" / "axang"
OUT.mkdir(parents=True, exist_ok=True)
FIBRO, SOUTH6, SOUTH7 = (ROOT / "results" / x for x in ("fibro", "south6", "south7"))
SOUTH8, SOUTH9 = ROOT / "results" / "south8", ROOT / "results" / "south9"
SAME, TOWARD = ROOT / "results" / "same", ROOT / "results" / "toward"
SENG, SENG2, SENG3, SENG4 = (ROOT / "results" / x for x in
                             ("seng", "seng2", "seng3", "seng4"))
MD2P = ROOT / "data" / "processed" / "md2"
MD2R = ROOT / "results" / "md2"
KDIM = 20
DEPTHS = (96, 107, 128, 157)
BOOT_SEED = 20260918
N_BOOT = 200


def dump(n, o):
    (OUT / n).write_text(json.dumps(R3._js(o), indent=2), encoding="utf-8")


def cosang(x, y):
    x, y = np.asarray(x, float).ravel(), np.asarray(y, float).ravel()
    nx, ny = np.linalg.norm(x), np.linalg.norm(y)
    if nx < 1e-12 or ny < 1e-12:
        return np.nan, np.nan
    c = float(np.dot(x, y) / (nx * ny))
    return c, float(np.degrees(np.arccos(np.clip(c, -1, 1))))


def cc_lists():
    t = pd.read_csv(SENG / "seurat_cc.genes_tirosh2016.tsv", sep="\t")
    return (t.loc[t.phase == "S", "symbol"].tolist(),
            t.loc[t.phase == "G2M", "symbol"].tolist())


def ruler_cc_index(frozen):
    """Ruler rows of the Tirosh genes, first occurrence by upper-case symbol."""
    sym = np.array([str(s).upper() for s in np.asarray(frozen["symbol"])], dtype=object)
    pos = {}
    for i, s in enumerate(sym):
        pos.setdefault(s, i)
    s_g, g_g = cc_lists()
    idx, missing = [], []
    for g in list(s_g) + list(g_g):
        p = pos.get(g.upper())
        (idx.append(p) if p is not None else missing.append(g))
    return np.array(idx, int), missing, len(s_g) + len(g_g)


def quartile_axis(score, pseudo_fn):
    """bottom-to-top-quartile direction, SENG's quartile_masks rule."""
    bot, top, info = S.quartile_masks(np.asarray(score, float))
    return pseudo_fn(top) - pseudo_fn(bot), info


# --------------------------------------------------------------------------
def build(log):
    frozen = load_frozen_ruler()
    n_ruler = int(np.asarray(frozen["mu"]).shape[0])
    L = np.asarray(np.load(SOUTH7 / "basis.npz", allow_pickle=True)["L_mid_full"],
                   np.float64)[:, :KDIM]
    anch = np.load(TOWARD / "anchors.npz", allow_pickle=True)
    cc_idx, cc_missing, cc_total = ruler_cc_index(frozen)
    log(f"[cc] {cc_idx.size}/{cc_total} Tirosh genes on the ruler; missing={cc_missing}")

    ax = {}
    # ---- aging axes
    ax["A1"] = np.asarray(anch["u"], float).ravel()
    ax["A2"] = unit(np.asarray(frozen["w"], float))[0]
    sp3 = R3.space3(log, need_extra=False)
    Z3, ni3 = sp3["Z"], sp3["ni"]
    ax["A3"] = Z3[ni3["EARLY_A"]] - Z3[ni3["LATE_A"]]
    sm = np.load(SAME / "panel_z.npz", allow_pickle=True)
    ax["A4"] = np.asarray(sm["v"], float).ravel()
    stored_a4_cos = float(np.asarray(sm["cos_v_u_toward"]).ravel()[0])

    # ---- C2, two variants: ruler-wide (this task) and SENG4's P3 (Sengstack-measured)
    e = np.zeros(n_ruler)
    e[cc_idx] = 1.0
    ax["C2"] = e / np.linalg.norm(e)
    pr4 = np.load(SENG4 / "prep.npz", allow_pickle=True)
    dst4 = np.asarray(pr4["dst"], int)
    loc4 = np.concatenate([np.asarray(pr4["s_cols"], int), np.asarray(pr4["g_cols"], int)])
    e4 = np.zeros(n_ruler)
    e4[dst4[loc4]] = 1.0
    ax["C2_seng4_P3"] = e4 / np.linalg.norm(e4)

    # ---- C1 Southard, read frozen
    s1 = np.load(SOUTH6 / "stage1.npz", allow_pickle=True)
    ax["C1_south"] = np.asarray(s1["p_unit"], float).ravel()

    # ---- C1 Sengstack, from SENG4's reserved-half quartile pseudobulks
    zf = sp3["zf"]
    for m in ("CRA", "CRI"):
        ax[f"C1_seng_{m}"] = zf(pr4[f"pb_P1_{m}_TOP"]) - zf(pr4[f"pb_P1_{m}_BOT"])

    # ---- C1 GTEx, built here (amendments F1, F2)
    g = load_gtex_z(log)
    Xg = np.asarray(g["pack"]["X"], np.float64)          # log2-CPM, donors x ruler
    age = np.asarray(g["age"]).astype(str)
    Zg = np.asarray(g["Z"], np.float64)
    sc_g = Xg[:, cc_idx].mean(axis=1)
    pseudo_g = lambda mask: Zg[mask].mean(axis=0)
    ax["C1_gtex"], qi_g = quartile_axis(sc_g, pseudo_g)
    anchor_m = np.isin(age, ["20-29", "30-39", "60-69", "70-79"])
    sc_ga = sc_g[anchor_m]
    Zga = Zg[anchor_m]
    ax["C1_gtex_anchors"], qi_ga = quartile_axis(sc_ga, lambda mk: Zga[mk].mean(axis=0))
    log(f"[C1_gtex] n_donors={Xg.shape[0]} bottom={qi_g['n_bottom']} top={qi_g['n_top']}; "
        f"anchors-only n={int(anchor_m.sum())} bottom={qi_ga['n_bottom']}")

    # ---- C1 Lu, built here from day-0 fibroblasts of both donors
    lu = lu_c1(frozen, cc_idx, log)
    ax["C1_lu"] = lu["axis"]

    meta = dict(cc_on_ruler=int(cc_idx.size), cc_total=cc_total,
                cc_missing_from_ruler=cc_missing,
                n_gtex_donors=int(Xg.shape[0]),
                n_gtex_anchor_donors=int(anchor_m.sum()),
                gtex_quartiles=qi_g, gtex_anchor_quartiles=qi_ga,
                n_south_nt_cells=int(np.asarray(s1["nt_cells"]).size),
                n_seng_half1=dict(CRA=int(pr4["h1_CRA"].size), CRI=int(pr4["h1_CRI"].size)),
                lu=lu["meta"], stored_a4_cos=stored_a4_cos,
                n_ruler=n_ruler)
    return ax, L, frozen, anch, meta, sp3


def lu_c1(frozen, cc_idx, log):
    """C1 for the Lu donor pair: day-0 fibroblasts of GM23815 and GM00731."""
    sym_r = np.array([str(s).upper() for s in np.asarray(frozen["symbol"])], dtype=object)
    pos_r = {}
    for i, s in enumerate(sym_r):
        pos_r.setdefault(s, i)
    mu = np.asarray(frozen["mu"], float)
    sd = np.asarray(frozen["sd"], float)
    tops, bots, meta = [], [], {}
    for donor in ("GM23815", "GM00731"):
        obs = pd.read_csv(MD2R / f"allcell_obs_{donor}.csv")
        z = np.load(MD2P / f"allcell_counts_{donor}.npz", allow_pickle=True)
        shp = tuple(int(x) for x in z["shape"])
        X = sparse.csr_matrix((z["data"], z["indices"], z["indptr"]), shape=shp)
        d0 = np.flatnonzero(obs.day.to_numpy() == 0)
        X0 = X[d0].astype(np.float64)
        del X
        sym_f = np.array([str(s).upper() for s in np.asarray(z["symbols"])], dtype=object)
        first = {}
        for i, s in enumerate(sym_f):
            first.setdefault(s, i)
        cc_file = np.array([first[sym_r[j]] for j in cc_idx if sym_r[j] in first], int)
        lib = obs.umi.to_numpy(np.float64)[d0]
        sc = S.module_score(X0, cc_file.tolist(), lib)
        bot, top, info = S.quartile_masks(sc)
        ruler_file = np.array([first.get(s, -1) for s in sym_r], int)
        have = ruler_file >= 0

        def pb(mask):
            v = np.zeros(mu.shape[0])
            s = np.asarray(X0[np.flatnonzero(mask)].sum(axis=0)).ravel()
            v[have] = s[ruler_file[have]]
            return v
        tops.append(pb(top))
        bots.append(pb(bot))
        meta[donor] = dict(n_day0=int(d0.size), n_bottom=info["n_bottom"],
                           n_top=info["n_top"], cc_genes_found=int(cc_file.size),
                           ruler_genes_present=int(have.sum()))
        log(f"[C1_lu] {donor} day0={d0.size} bottom={info['n_bottom']} top={info['n_top']} "
            f"ruler_present={int(have.sum())}")
        del X0
    ct, cb = sum(tops), sum(bots)
    miss = ~(ct + cb > 0)
    from toward_run import tmm_logcpm_quiet, zscore_frozen
    C = np.vstack([ct, cb])
    lc, nf = tmm_logcpm_quiet(C)
    Zl, missing = zscore_frozen(lc, frozen, C)
    return dict(axis=Zl[0] - Zl[1], meta=meta)


# --------------------------------------------------------------------------
def stage0():
    log = K.Log(OUT / "stage0.log")
    ax, L, frozen, anch, meta, sp3 = build(log)
    ck = []

    def num(item, stored, got):
        ck.append(dict(item=item, kind="numeric", stored=stored, got=got,
                       ok=bool(got is not None and np.isfinite(got)
                               and abs(got - stored) <= 1e-4 * max(1.0, abs(stored)))))

    def cnt(item, stored, got):
        ck.append(dict(item=item, kind="count", stored=stored, got=got,
                       ok=bool(got == stored)))

    a3_20 = ax["A3"] @ L
    for m, st in (("CRA", 0.941708), ("CRI", 0.936287)):
        num(f"SENG4 cos(P1,a) {m}", st, cosang(ax[f"C1_seng_{m}"] @ L, a3_20)[0])
    num("SENG4 cos(P3,a)", 0.474385, cosang(ax["C2_seng4_P3"] @ L, a3_20)[0])
    num("SENG4 cos(a,u) full-gene", -0.040273, cosang(ax["A3"], ax["A1"])[0])
    num("SENG4 cos(a,u) S2A_k20", -0.199293, cosang(a3_20, ax["A1"] @ L)[0])
    s9 = pd.read_csv(SOUTH9 / "d2_residualized_recovery.csv")
    r9 = s9[(s9.space == "S2A_k20") & (s9.depth == 106) & (s9.variant == "in_space")].iloc[0]
    num("SOUTH9 res power_25 @106", 0.990, float(r9.res_power_25))
    num("SOUTH9 raw power_25 @106", 0.980, float(r9.raw_power_25))
    f8 = pd.read_csv(SOUTH8 / "d3c_all_1836_umi_calls.csv")
    cnt("SOUTH8 all_1836 m", 1836, int(len(f8)))
    cnt("SOUTH8 all_1836 amended_call_S2A_k20", 125, int(f8["amended_call_S2A_k20"].sum()))
    cnt("SOUTH8 all_1836 below_depth_floor_S2A_k20", 153,
        int(f8["below_depth_floor_S2A_k20"].sum()))
    num("E1 stored cos(A4, GTEx age axis)", 0.039, round(meta["stored_a4_cos"], 3))
    num("E1 cos(A4,u) re-derived vs stored", meta["stored_a4_cos"],
        cosang(ax["A4"], ax["A1"])[0])
    cnt("anchors n_young", 113, int(np.asarray(anch["n_young"]).ravel()[0]))
    cnt("anchors n_old", 231, int(np.asarray(anch["n_old"]).ravel()[0]))

    passed = all(c["ok"] for c in ck)
    pd.DataFrame(ck).to_csv(OUT / "stage0_reproduction.csv", index=False)
    dump("stage0_summary.json", dict(passed=bool(passed), checks=ck, meta=meta))
    np.savez_compressed(OUT / "axes.npz", **{k: v for k, v in ax.items()}, L=L)
    log(f"[stage0] PASS={passed}")
    if not passed:
        bad = [c["item"] for c in ck if not c["ok"]]
        dump("STOP.json", dict(key="prior_results_not_reproduced", mismatches=bad))
        log(f"[stage0] STOP {bad}")
    sp3["ctx"]["f"].close()
    log.close()


# --------------------------------------------------------------------------
AG = [("A1", "GTEx donor age"), ("A2", "frozen ruler"),
      ("A3", "Sengstack passage"), ("A4", "Lu donor pair")]
PR = [("C1_gtex", "C1 GTEx donors"), ("C1_gtex_anchors", "C1 GTEx anchor donors"),
      ("C1_south", "C1 Southard NT cells"), ("C1_seng_CRA", "C1 Sengstack NT CRA"),
      ("C1_seng_CRI", "C1 Sengstack NT CRI"), ("C1_lu", "C1 Lu day-0"),
      ("C2", "C2 Tirosh gene list"), ("C2_seng4_P3", "C2 as SENG4's P3")]


def _boots(ax, L, frozen, log):
    """200-draw bootstraps on the axes whose defining units are available (F9)."""
    rng = np.random.default_rng(BOOT_SEED)
    g = load_gtex_z(log)
    Zg = np.asarray(g["Z"], np.float64)
    Xg = np.asarray(g["pack"]["X"], np.float64)
    age = np.asarray(g["age"]).astype(str)
    cc_idx, _, _ = ruler_cc_index(frozen)
    sc_g = Xg[:, cc_idx].mean(axis=1)
    yi = np.flatnonzero(np.isin(age, ["20-29", "30-39"]))
    oi = np.flatnonzero(np.isin(age, ["60-69", "70-79"]))
    n_d = Zg.shape[0]
    B = {k: [] for k in ("A1", "C1_gtex")}
    for _ in range(N_BOOT):
        B["A1"].append(Zg[rng.choice(yi, yi.size, True)].mean(0)
                       - Zg[rng.choice(oi, oi.size, True)].mean(0))
        d = rng.choice(n_d, n_d, True)
        v, _ = quartile_axis(sc_g[d], lambda mk: Zg[d][mk].mean(axis=0))
        B["C1_gtex"].append(v)
    # Sengstack NT half-one cells
    sp3 = R3.space3(log, need_extra=False)
    ctx, zf = sp3["ctx"], sp3["zf"]
    pr4 = np.load(SENG4 / "prep.npz", allow_pickle=True)
    dst4 = np.asarray(pr4["dst"], int)
    nr = sp3["n_ruler"]
    cc_loc = np.concatenate([np.asarray(pr4["s_cols"], int),
                             np.asarray(pr4["g_cols"], int)])
    for m in ("CRA", "CRI"):
        rows = np.asarray(pr4[f"h1_{m}"], int)
        blk = np.asarray(K.submatrix(ctx["f"], rows, ctx["exp"]).todense(), np.float64)
        lib = np.asarray(pr4["lib"], float)[rows]
        sc = S.module_score(sparse.csr_matrix(blk), cc_loc.tolist(), lib)
        B[f"C1_seng_{m}"] = []
        for _ in range(N_BOOT):
            d = rng.choice(blk.shape[0], blk.shape[0], True)
            bot, top, _i = S.quartile_masks(sc[d])
            pb = lambda mk: K.to_ruler(blk[d][np.flatnonzero(mk)].sum(axis=0), dst4, nr)
            B[f"C1_seng_{m}"].append(zf(pb(top)) - zf(pb(bot)))
        log(f"[boot] C1_seng_{m} done")
        del blk
    ctx["f"].close()
    return B


def d1():
    log = K.Log(OUT / "d1.log")
    ax, L, frozen, anch, meta, sp3 = build(log)
    sp3["ctx"]["f"].close()
    # consistency: the frozen anchors equal the pack centroids
    gg = load_gtex_z(log)
    Zg, age = np.asarray(gg["Z"], np.float64), np.asarray(gg["age"]).astype(str)
    cy = Zg[np.isin(age, ["20-29", "30-39"])].mean(0)
    co = Zg[np.isin(age, ["60-69", "70-79"])].mean(0)
    chk = dict(max_abs_c_young=float(np.max(np.abs(cy - np.asarray(anch["c_young"], float)))),
               max_abs_c_old=float(np.max(np.abs(co - np.asarray(anch["c_old"], float)))))
    log(f"[check] anchors vs pack centroids {chk}")
    B = _boots(ax, L, frozen, log)

    rows = []
    for space, P in (("full_gene", None), ("S2A_k20", L)):
        pj = lambda v: (v if P is None else v @ P)
        names = [k for k, _ in AG] + [k for k, _ in PR]
        for i, a in enumerate(names):
            for b in names[i:]:
                if a == b:
                    continue
                c, d = cosang(pj(ax[a]), pj(ax[b]))
                kind = ("aging-aging" if a in dict(AG) and b in dict(AG) else
                        "prolif-prolif" if a not in dict(AG) and b not in dict(AG)
                        else "aging-prolif")
                lo = hi = np.nan
                if a in B and b in B:
                    v = [cosang(pj(x), pj(y))[0] for x, y in zip(B[a], B[b])]
                elif a in B:
                    v = [cosang(pj(x), pj(ax[b]))[0] for x in B[a]]
                elif b in B:
                    v = [cosang(pj(ax[a]), pj(y))[0] for y in B[b]]
                else:
                    v = None
                if v:
                    lo, hi = float(np.percentile(v, 2.5)), float(np.percentile(v, 97.5))
                rows.append(dict(space=space, kind=kind, axis_a=a, axis_b=b, cos=c,
                                 angle_deg=d,
                                 frac_a_surviving=float(np.sqrt(max(0.0, 1 - c ** 2))),
                                 pct_removed=100.0 * (1 - float(np.sqrt(max(0.0, 1 - c ** 2)))),
                                 ci_lo=lo, ci_hi=hi,
                                 ci_valid=bool(np.isnan(lo) or (lo <= c <= hi)),
                                 bootstrapped=bool(v is not None)))
        log(f"[d1] {space} done")
    df = pd.DataFrame(rows)
    df.to_csv(OUT / "d1_angles.csv", index=False)

    get = lambda sp_, a, b: df[(df.space == sp_) &
                               (((df.axis_a == a) & (df.axis_b == b)) |
                                ((df.axis_a == b) & (df.axis_b == a)))].iloc[0]
    g1 = get("S2A_k20", "A1", "C1_gtex")
    d_ang = float(g1.angle_deg)
    p_ang = 19.66
    gate = ("SEPARABLE" if d_ang >= 45 else
            "INTERMEDIATE" if d_ang >= 25 else "ENTANGLED")
    g2 = get("S2A_k20", "A1", "C2")
    d_c2 = float(g2.angle_deg)
    gate_c2 = ("SEPARABLE" if d_c2 >= 45 else
               "INTERMEDIATE" if d_c2 >= 25 else "ENTANGLED")
    # triangles (F3)
    tri = {}
    for nm, third in (("T1_A1_A3_C1seng", "C1_seng_CRA"), ("T2_A1_A3_C1gtex", "C1_gtex")):
        t_ab = float(get("S2A_k20", "A1", "A3").angle_deg)
        t_az = float(get("S2A_k20", "A1", third).angle_deg)
        t_bz = float(get("S2A_k20", "A3", third).angle_deg)
        tri[nm] = dict(third=third, angle_A1_A3=t_ab, angle_A1_third=t_az,
                       angle_A3_third=t_bz, lower=abs(t_az - t_bz),
                       upper=t_az + t_bz,
                       consistent=bool(abs(t_az - t_bz) - 1e-9 <= t_ab <= t_az + t_bz + 1e-9))
    dump("d1_summary.json", dict(
        gate=gate, d_angle=d_ang, d_cos=float(g1.cos),
        d_ci=[float(g1.ci_lo), float(g1.ci_hi)], p_angle=p_ang,
        pct_removed=float(g1.pct_removed), gate_c2=gate_c2, d_angle_c2=d_c2,
        gates_agree=bool(gate == gate_c2), triangles=tri,
        anchors_vs_pack=chk, meta=meta))
    log(f"[d1] gate={gate} d={d_ang:.2f} (C2 {gate_c2} {d_c2:.2f}) "
        f"tri={[(k, v['consistent']) for k, v in tri.items()]}")
    log.close()


# --------------------------------------------------------------------------
def d2():
    log = K.Log(OUT / "d2.log")
    ax = np.load(OUT / "axes.npz", allow_pickle=True)
    L = np.asarray(ax["L"], float)
    frozen = load_frozen_ruler()
    anch = np.load(TOWARD / "anchors.npz", allow_pickle=True)
    d_age = np.asarray(anch["c_young"], float) - np.asarray(anch["c_old"], float)
    tot = float(np.dot(d_age, d_age))
    sp3 = R3.space3(log, need_extra=False)
    nr = sp3["n_ruler"]
    out = {}
    # measured-gene sets per dataset
    meas = {}
    meas["Sengstack"] = np.zeros(nr, bool)
    meas["Sengstack"][np.asarray(np.load(SENG2 / "panel_counts.npz")["dst"], int)] = True
    s1 = np.load(SOUTH6 / "stage1.npz", allow_pickle=True)
    meas["Southard"] = np.zeros(nr, bool)
    meas["Southard"][np.asarray(s1["matched_ruler"], int)] = True
    sym_r = np.array([str(x).upper() for x in np.asarray(frozen["symbol"])], dtype=object)
    zl = np.load(MD2P / "allcell_counts_GM23815.npz", allow_pickle=True)
    luset = set(str(x).upper() for x in np.asarray(zl["symbols"]))
    meas["Lu"] = np.array([s in luset for s in sym_r], bool)
    for k, m in meas.items():
        out[k] = dict(present=int(m.sum()), absent=int((~m).sum()),
                      pct_gap_sq_in_unmeasured=100.0 * float(np.dot(d_age[~m], d_age[~m]) / tot))
        log(f"[d2] {k} present={int(m.sum())} absent={int((~m).sum())} "
            f"pct_gap_sq_unmeasured={out[k]['pct_gap_sq_in_unmeasured']:.3f}")
    # PD32 against the GTEx centroids (E2, Supplementary Note 1 mirror)
    Z3, ni3 = sp3["Z"], sp3["ni"]
    cy, co = np.asarray(anch["c_young"], float), np.asarray(anch["c_old"], float)
    z32 = Z3[ni3["LATE_A"]]
    gap = float(np.linalg.norm(cy - co))
    out["PD32_vs_GTEx"] = dict(
        dist_to_young=float(np.linalg.norm(z32 - cy)),
        dist_to_old=float(np.linalg.norm(z32 - co)),
        young_old_gap=gap,
        ratio_to_gap_young=float(np.linalg.norm(z32 - cy) / gap),
        ratio_to_gap_old=float(np.linalg.norm(z32 - co) / gap))
    log(f"[d2] PD32 {out['PD32_vs_GTEx']}")
    # whether an earlier SENG task ran a GTEx-anchored test
    scr = pd.read_csv(SENG / "screen.csv")
    out["earlier_gtex_anchored"] = dict(
        file="results/seng/screen.csv",
        columns=["cos_d_u", "p_cos_d_u"],
        n=int(len(scr)),
        n_p_le_05=int((scr.p_cos_d_u <= 0.05).sum()),
        median_cos_d_u=float(scr.cos_d_u.median()),
        max_cos_d_u=float(scr.cos_d_u.max()),
        note=("SENG computed the cosine of each displacement with the GTEx "
              "donor-age direction u and a permutation p, which is a DIRECTION "
              "test against donor age. It did not run an approach/distance test "
              "anchored on the GTEx centroids. D3 of this task is the first."))
    log(f"[d2] earlier: {out['earlier_gtex_anchored']['n_p_le_05']} of "
        f"{out['earlier_gtex_anchored']['n']} with p<=0.05 on cos(d,u)")
    dump("d2_summary.json", out)
    sp3["ctx"]["f"].close()
    log.close()


# --------------------------------------------------------------------------
def d3():
    log = K.Log(OUT / "d3.log")
    axz = np.load(OUT / "axes.npz", allow_pickle=True)
    L = np.asarray(axz["L"], float)
    anch = np.load(TOWARD / "anchors.npz", allow_pickle=True)
    cy, co = np.asarray(anch["c_young"], float), np.asarray(anch["c_old"], float)
    r_full = co - cy                      # origin minus target, SOUTH's convention
    r20 = r_full @ L
    gap20 = float(np.linalg.norm((cy - co) @ L))
    ahat20 = ((cy - co) @ L) / gap20
    sp3 = R3.space3(log, need_extra=False)
    ctx, zf, Z3, ni3 = sp3["ctx"], sp3["zf"], sp3["Z"], sp3["ni"]
    dst, nr = np.asarray(np.load(SENG2 / "panel_counts.npz")["dst"], int), sp3["n_ruler"]
    C = np.load(SENG2 / "panel_counts.npz")["C"]
    blk, pool = {}, {}
    for m in ("CRA", "CRI"):
        blk[m] = np.asarray(K.submatrix(ctx["f"], ctx["nt_mod"][m], ctx["exp"]).todense(),
                            np.float64)
        pool[m] = C[ni3[f"NT_{m}"]]
        log(f"[d3] NT {m} {blk[m].shape}")

    def draws(m, n, reps=200, seed=K.SEED):
        rng = np.random.default_rng(seed)
        D = []
        for _ in range(reps):
            t = rng.choice(blk[m].shape[0], size=n, replace=False)
            sub = K.to_ruler(blk[m][t].sum(axis=0), dst, nr)
            D.append(zf(sub) - zf(pool[m] - sub))
        return np.vstack(D) @ L

    def delta(d20):
        r2 = float(r20 @ r20)
        rn = float(np.sqrt(r2))
        dots = d20 @ r20
        qn2 = np.einsum("ij,ij->i", d20, d20)
        return np.sqrt(np.maximum(r2 + 2.0 * dots + qn2, 0.0)) - rn

    # D3a plant gate (0% added per amendment F4)
    prow = []
    for n in DEPTHS:
        for m in ("CRA", "CRI"):
            D = draws(m, n)
            fl = None
            for dose in (0.0, 0.25, 0.50, 1.0):
                dl = delta(D + dose * gap20 * ahat20)
                if dose == 0.0:
                    fl = float(np.percentile(dl, 5))
                prow.append(dict(library=m, n_cells=n, dose=dose, floor=fl,
                                 delta_median=float(np.median(dl)),
                                 power=float(np.mean((dl < fl) & (dl < 0)))))
        log(f"[d3a] n={n} done")
    pdf = pd.DataFrame(prow)
    pdf.to_csv(OUT / "d3a_plant.csv", index=False)
    p25 = float(pdf[pdf.dose == 0.25].power.min())
    p50 = float(pdf[pdf.dose == 0.50].power.min())
    gate_a = ("QUALIFIED" if p25 >= 0.80 else
              "LIMITED" if p50 >= 0.80 else "BLIND")
    log(f"[d3a] gate={gate_a} p25={p25:.3f} p50={p50:.3f}")

    # D3b
    floors = {}
    by = {"CRA": [], "CRI": []}
    for m, lab, idx in ctx["qual"]:
        by[m].append((lab, int(idx.size)))
    for m in ("CRA", "CRI"):
        ns = sorted({n for _, n in by[m]})
        for j, n in enumerate(ns):
            floors[(m, n)] = float(np.percentile(delta(draws(m, n)), 5))
            if (j + 1) % 25 == 0:
                log(f"[d3b floors] {m} {j + 1}/{len(ns)}")
    rows = []
    for m, lab, idx in ctx["qual"]:
        n = int(idx.size)
        d_full = Z3[ni3[lab]] - Z3[ni3[f"NT_{m}"]]
        d20 = (d_full @ L)[None, :]
        dl = float(delta(d20)[0])
        fl = floors[(m, n)]
        rec = dict(modality=m, perturbation=lab, n_cells=n, delta=dl, floor=fl,
                   below_floor=bool(dl < fl), decreased=bool(dl < 0),
                   called=bool(dl < fl and dl < 0))
        for nm, v in (("u", axz["A1"]), ("A3", axz["A3"]),
                      (f"C1", axz[f"C1_seng_{m}"])):
            c, a = cosang(d_full @ L, np.asarray(v, float) @ L)
            rec[f"cos_d_{nm}"], rec[f"angle_d_{nm}"] = c, a
        rows.append(rec)
    df = pd.DataFrame(rows)
    allv = df.n_cells.to_numpy()
    df["cell_pctl_in_screen"] = [float((allv <= x).mean() * 100) for x in allv]
    df.to_csv(OUT / "d3b_table.csv", index=False)
    kc = int(df.called.sum())
    lo, hi = K.binom_ci(kc, len(df))
    summ = dict(n=len(df), below_floor=int(df.below_floor.sum()),
                decreased=int(df.decreased.sum()), called=kc,
                called_ci=[lo, hi], chance_expected=0.05 * len(df))
    for m in ("CRA", "CRI"):
        sub = df[df.modality == m]
        summ[m] = dict(n=len(sub), below_floor=int(sub.below_floor.sum()),
                       decreased=int(sub.decreased.sum()), called=int(sub.called.sum()),
                       chance_expected=0.05 * len(sub))
    named = {}
    for lab in K.NAMED:
        r = df[df.perturbation == lab]
        if len(r):
            named[lab] = {c: R3._js(r.iloc[0][c]) for c in df.columns
                          if c != "perturbation"}
    k_named = sum(1 for l in named if named[l]["called"])
    gate_b = "APPROACH" if k_named >= 2 else "NO APPROACH"
    # SOUTH8's stored result for the same gene names
    s8 = pd.read_csv(SOUTH8 / "d3c_all_1836_umi_calls.csv")
    south_named = {}
    for gname in ("E2F3", "EZH2", "STAT3", "ZFX"):
        r = s8[s8.factor == gname]
        south_named[gname] = (dict(in_southard=False) if r.empty else
                              dict(in_southard=True, n_cells=int(r.n_cells.iloc[0]),
                                   delta_S2A_k20=float(r.delta_S2A_k20.iloc[0]),
                                   below_depth_floor=bool(r.below_depth_floor_S2A_k20.iloc[0]),
                                   amended_call=bool(r.amended_call_S2A_k20.iloc[0])))
    dump("d3_summary.json", dict(gate_a=gate_a, p25=p25, p50=p50, gate_b=gate_b,
                                 k_named=k_named, summary=summ, named=named,
                                 southard_same_names=south_named,
                                 gtex_gap_full=float(np.linalg.norm(cy - co)),
                                 gtex_gap_k20=gap20))
    log(f"[d3] gate_a={gate_a} gate_b={gate_b} k_named={k_named} called={kc} of {len(df)}")
    ctx["f"].close()
    log.close()


STAGES = dict(stage0=stage0, d1=d1, d2=d2, d3=d3)

if __name__ == "__main__":
    STAGES[sys.argv[1]]()
