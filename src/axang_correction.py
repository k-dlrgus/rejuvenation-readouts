"""AXANG correction: rebuild C1_lu with the project-standard Ensembl-first
ruler alignment and recompute only the affected quantities.

Writes results/axang/CORRECTION.md and results/axang/d1_angles_corrected_C1lu.csv.
The original d1_angles.csv, d1_summary.json and d2_summary.json are NOT edited;
their values are marked superseded in CORRECTION.md.
"""
from __future__ import annotations
import json, sys
from pathlib import Path
import numpy as np, pandas as pd
from scipy import sparse

SRC = Path(__file__).resolve().parent
sys.path.insert(0, str(SRC))
ROOT = SRC.parent
OUT = ROOT / "results" / "axang"
MD2P, MD2R = ROOT / "data" / "processed" / "md2", ROOT / "results" / "md2"

import seng_run as S
from fibro2_common import load_frozen_ruler
from brain_phase1_common import strip_ensembl
from toward_run import tmm_logcpm_quiet, zscore_frozen


def cosang(x, y):
    x, y = np.asarray(x, float).ravel(), np.asarray(y, float).ravel()
    c = float(np.dot(x, y) / (np.linalg.norm(x) * np.linalg.norm(y)))
    return c, float(np.degrees(np.arccos(np.clip(c, -1, 1))))


def main():
    fr = load_frozen_ruler()
    rs = np.array([str(s).upper() for s in np.asarray(fr["symbol"])], dtype=object)
    re_ = np.asarray(fr["ensembl"]).astype(str)
    nr = rs.size
    cc = pd.read_csv(ROOT / "results" / "seng" / "seurat_cc.genes_tirosh2016.tsv", sep="\t")
    cc_sym = [s.upper() for s in cc.symbol.tolist()]
    rpos = {}
    for i, s in enumerate(rs):
        rpos.setdefault(s, i)
    cc_ruler = np.array([rpos[g] for g in cc_sym if g in rpos], int)

    tops, bots, meta = [], [], {}
    align = {}
    for donor in ("GM23815", "GM00731"):
        obs = pd.read_csv(MD2R / f"allcell_obs_{donor}.csv")
        z = np.load(MD2P / f"allcell_counts_{donor}.npz", allow_pickle=True)
        shp = tuple(int(x) for x in z["shape"])
        X = sparse.csr_matrix((z["data"], z["indices"], z["indptr"]), shape=shp)
        d0 = np.flatnonzero(obs.day.to_numpy() == 0)
        X0 = X[d0].astype(np.float64)
        del X
        fsym = np.array([str(x).upper() for x in np.asarray(z["symbols"])], dtype=object)
        fens = np.array([strip_ensembl(str(x)) for x in np.asarray(z["gene_id"])])
        pos_e, pos_s = {}, {}
        for i, e in enumerate(fens):
            if e.startswith("ENS"):
                pos_e.setdefault(e, i)
        for i, s in enumerate(fsym):
            pos_s.setdefault(s, i)
        idx = np.full(nr, -1, int)
        n_e = n_s = 0
        for j in range(nr):
            i = pos_e.get(re_[j])
            if i is not None:
                idx[j] = i; n_e += 1; continue
            i = pos_s.get(rs[j])
            if i is not None:
                idx[j] = i; n_s += 1
        have = idx >= 0
        align[donor] = dict(present=int(have.sum()), absent=int((~have).sum()),
                            by_ensembl=n_e, by_symbol=n_s)
        cc_file = np.array([idx[j] for j in cc_ruler if idx[j] >= 0], int)
        lib = obs.umi.to_numpy(np.float64)[d0]
        sc = S.module_score(X0, cc_file.tolist(), lib)
        bot, top, info = S.quartile_masks(sc)

        def pb(mask):
            v = np.zeros(nr)
            ssum = np.asarray(X0[np.flatnonzero(mask)].sum(axis=0)).ravel()
            v[have] = ssum[idx[have]]
            return v
        tops.append(pb(top)); bots.append(pb(bot))
        meta[donor] = dict(n_day0=int(d0.size), n_bottom=info["n_bottom"],
                           n_top=info["n_top"], cc_genes_found=int(cc_file.size),
                           **align[donor])
        print(f"[C1_lu*] {donor} day0={d0.size} ruler_present={int(have.sum())} "
              f"(ensembl {n_e} + symbol {n_s}) cc_genes={cc_file.size}")
        del X0
    C = np.vstack([sum(tops), sum(bots)])
    lc, nf = tmm_logcpm_quiet(C)
    Zl, _ = zscore_frozen(lc, fr, C)
    c1_lu_new = Zl[0] - Zl[1]

    axz = np.load(OUT / "axes.npz", allow_pickle=True)
    L = np.asarray(axz["L"], float)
    ax = {k: np.asarray(axz[k], float) for k in axz.files if k != "L"}
    c1_lu_old = ax["C1_lu"]
    ax["C1_lu"] = c1_lu_new

    rows = []
    others = [k for k in ax if k != "C1_lu"]
    for space, P in (("full_gene", None), ("S2A_k20", L)):
        pj = lambda v: (v if P is None else v @ P)
        for k in others:
            c_new, a_new = cosang(pj(c1_lu_new), pj(ax[k]))
            c_old, a_old = cosang(pj(c1_lu_old), pj(ax[k]))
            rows.append(dict(space=space, axis_a=k, axis_b="C1_lu",
                             cos_superseded=c_old, angle_superseded=a_old,
                             cos_corrected=c_new, angle_corrected=a_new,
                             angle_from_parallel_corrected=float(
                                 np.degrees(np.arccos(min(1.0, abs(c_new))))),
                             frac_surviving_corrected=float(np.sqrt(max(0.0, 1 - c_new ** 2))),
                             pct_removed_corrected=100.0 * (1 - float(np.sqrt(max(0.0, 1 - c_new ** 2))))))
    df = pd.DataFrame(rows)
    df.to_csv(OUT / "d1_angles_corrected_C1lu.csv", index=False)
    c_self, a_self = cosang(c1_lu_old, c1_lu_new)

    anch = np.load(ROOT / "results" / "toward" / "anchors.npz", allow_pickle=True)
    dage = np.asarray(anch["c_young"], float) - np.asarray(anch["c_old"], float)
    tot = float(np.dot(dage, dage))
    have_r = np.zeros(nr, bool)
    # a ruler gene is measured in Lu if either donor file carries it
    for donor in ("GM23815", "GM00731"):
        pass
    # recompute the union mask from GM23815 (both files share the 10x reference)
    z = np.load(MD2P / "allcell_counts_GM23815.npz", allow_pickle=True)
    fens = np.array([strip_ensembl(str(x)) for x in np.asarray(z["gene_id"])])
    fsym = np.array([str(x).upper() for x in np.asarray(z["symbols"])], dtype=object)
    se, ss = set(e for e in fens if e.startswith("ENS")), set(fsym.tolist())
    have_r = np.array([(re_[j] in se) or (rs[j] in ss) for j in range(nr)], bool)
    pct_new = 100.0 * float(np.dot(dage[~have_r], dage[~have_r]) / tot)
    d2 = json.loads((OUT / "d2_summary.json").read_text())
    out = dict(lu_present=int(have_r.sum()), lu_absent=int((~have_r).sum()),
               pct_gap_sq_unmeasured_corrected=pct_new,
               pct_gap_sq_unmeasured_superseded=d2["Lu"]["pct_gap_sq_in_unmeasured"],
               cos_old_new_C1lu=c_self, angle_old_new_C1lu=a_self,
               per_donor=meta)
    (OUT / "correction_values.json").write_text(json.dumps(out, indent=2))
    print(f"[coverage*] present={int(have_r.sum())} absent={int((~have_r).sum())} "
          f"pct_gap_unmeasured {d2['Lu']['pct_gap_sq_in_unmeasured']:.3f} -> {pct_new:.3f}")
    print(f"[C1_lu*] cos(old, new) = {c_self:.6f} ({a_self:.2f} deg)")


if __name__ == "__main__":
    main()
