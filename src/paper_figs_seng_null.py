"""SENG_WHY N1 null: 200 gene-permuted cos(v, u) for LATE = PD32 and PD26.

v = z(EARLY_A) − z(LATE_A) on the seng_why panel (same groups, TMM, frozen z; src/seng_why.py lines
for run_late). u = frozen TOWARD axis. Null = cos(v, perm_i(u)), perm_i from Generator(20260914).
No fit. PD32 v must equal results/seng/axes.npz v; cos and p must equal results/seng_why/axis_by_late.csv.
"""
from __future__ import annotations

import sys
from pathlib import Path

import h5py
import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parent))
from fibro2_common import load_frozen_ruler  # noqa: E402
from seng_run import (  # noqa: E402
    ANCHORS_PATH, H5_PATH, MIN_CELLS, N_PERM, SEED,
    cos_against_perms, parse_wide_obs, quartile_masks, read_obs_var, read_sparse_group,
    ruler_column_index, split_half, sum_cells,
)
from seng_why import load_cc_from_tsv, module_score  # noqa: E402
from toward_run import tmm_logcpm_quiet, zscore_frozen  # noqa: E402

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "results" / "paper_figs"
AXIS_CSV = ROOT / "results" / "seng_why" / "axis_by_late.csv"
AXES_NPZ = ROOT / "results" / "seng" / "axes.npz"
TOL = 1e-6


def _quiet(_msg):
    return None


def build_panel_inputs():
    frozen = load_frozen_ruler()
    n_ruler = int(np.asarray(frozen["mu"]).shape[0])
    with h5py.File(H5_PATH, "r") as f:
        obs = read_obs_var(f["obs"])
        var = read_obs_var(f["var"])
        mat, _ = read_sparse_group(f["X"])
    parsed = parse_wide_obs(obs, _quiet)
    labels = parsed["labels"].to_numpy()
    kind = parsed["kind"].to_numpy()
    pd_labels = parsed["pd_labels"].to_numpy()
    is_wt = np.asarray(parsed["is_wt"], bool)
    modality = np.asarray(parsed["modality"], dtype=object)

    symbols = var.index.astype(str).to_numpy()
    ruler_idx = ruler_column_index(symbols, symbols, frozen, _quiet)
    present = np.flatnonzero(ruler_idx >= 0)
    src = ruler_idx[present]
    lib_total = np.asarray(mat.sum(axis=1), np.float64).ravel()
    csr = mat[:, src].tocsr()
    dst = present

    early_idx = np.flatnonzero(is_wt & (pd_labels == "PD14"))
    late = {p: np.flatnonzero(is_wt & (pd_labels == p)) for p in ("PD32", "PD26")}
    early_a, _ = split_half(early_idx, SEED)
    late_a = {p: split_half(idx, SEED)[0] for p, idx in late.items()}

    s_genes, g2m_genes = load_cc_from_tsv()
    pos = {}
    for i, s in enumerate(np.array([str(s).upper() for s in symbols])):
        pos.setdefault(s, i)
    file_to_sub = {int(src[i]): i for i in range(len(src))}
    s_sub = [file_to_sub[pos[g]] for g in s_genes if g in pos and pos[g] in file_to_sub]
    g_sub = [file_to_sub[pos[g]] for g in g2m_genes if g in pos and pos[g] in file_to_sub]
    score = module_score(csr, s_sub, lib_total) + module_score(csr, g_sub, lib_total)

    nt_idx, q_idx = {}, {}
    for mod in ("CRA", "CRI"):
        nt = np.flatnonzero((kind == "NT") & (modality == mod))
        nt_idx[mod] = nt
        bottom, top, info = quartile_masks(score[nt])
        if bottom is None:
            raise SystemExit(f"quartile failed {mod} {info}")
        q_idx[mod] = dict(bottom=nt[bottom], top=nt[top])

    qual = []
    for lab, sub in pd.Series(np.arange(len(labels))).groupby(labels, sort=True).groups.items():
        lab = str(lab)
        if not (lab.startswith("CRA_") or lab.startswith("CRI_")):
            continue
        idx = np.asarray(sub, int)
        mod = "CRA" if lab.startswith("CRA_") else "CRI"
        if not np.all(modality[idx] == mod) or idx.size < MIN_CELLS:
            continue
        qual.append((mod, lab, idx))
    qual.sort(key=lambda t: (t[0], t[1]))
    return dict(frozen=frozen, n_ruler=n_ruler, csr=csr, dst=dst, early_a=early_a, late_a=late_a,
                nt_idx=nt_idx, q_idx=q_idx, qual=qual)


def v_for_late(P, late_name):
    groups = [("EARLY_A", P["early_a"]), ("LATE_A", P["late_a"][late_name])]
    for mod in ("CRA", "CRI"):
        groups.append((f"NT_{mod}", P["nt_idx"][mod]))
        groups.append((f"TOP_{mod}", P["q_idx"][mod]["top"]))
        groups.append((f"BOT_{mod}", P["q_idx"][mod]["bottom"]))
    for _, lab, idx in P["qual"]:
        groups.append((lab, idx))
    C = np.vstack([sum_cells(P["csr"], idx, P["dst"], P["n_ruler"]) for _, idx in groups])
    logcpm, _ = tmm_logcpm_quiet(C)
    Z, _ = zscore_frozen(logcpm, P["frozen"], C)
    return Z[0] - Z[1]


def main():
    OUT.mkdir(parents=True, exist_ok=True)
    stored = pd.read_csv(AXIS_CSV).set_index("late")
    axes = np.load(AXES_NPZ, allow_pickle=True)
    u = np.asarray(np.load(ANCHORS_PATH, allow_pickle=True)["u"], float).ravel()
    rng_u = np.random.default_rng(SEED)
    Uperm = np.empty((N_PERM, u.shape[0]), np.float64)
    for i in range(N_PERM):
        Uperm[i] = rng_u.permutation(u)

    P = build_panel_inputs()
    lines, cols = [], {}
    for late in ("PD32", "PD26"):
        v = v_for_late(P, late)
        if late == "PD32":
            dv = float(np.max(np.abs(v - np.asarray(axes["v"], float))))
            lines.append(f"[PD32] max |v - results/seng/axes.npz v|={dv:.2e}")
            if dv > TOL:
                raise SystemExit(f"STOP: PD32 v does not match axes.npz (max {dv:.3e})")
        cos, p, null = cos_against_perms(v, u, Uperm)
        dc = abs(cos - float(stored.loc[late, "cos_v_u"]))
        dp = abs(p - float(stored.loc[late, "p_n1"]))
        dn = abs(float(np.linalg.norm(v)) - float(stored.loc[late, "v_norm"]))
        lines.append(f"[{late}] cos(v,u) stored={stored.loc[late, 'cos_v_u']:.12g} recomputed={cos:.12g} |diff|={dc:.2e}; "
                     f"p stored={stored.loc[late, 'p_n1']:.12g} recomputed={p:.12g} |diff|={dp:.2e}; "
                     f"||v|| |diff|={dn:.2e}; null mean={np.mean(null):+.5f} sd={np.std(null):.5f} "
                     f"n_null_ge_obs={int(np.sum(null >= cos))}")
        if max(dc, dp) > TOL or dn > 1e-6 * max(1.0, float(stored.loc[late, "v_norm"])):
            print("\n".join(lines))
            raise SystemExit(f"STOP: {late} cos(v,u) not reproduced")
        cols[f"null_cos_v_u_{late}"] = null
        cols[f"obs_cos_v_u_{late}"] = np.full(N_PERM, cos)
        cols[f"p_n1_{late}"] = np.full(N_PERM, p)
    df = pd.DataFrame(dict(perm_index=np.arange(N_PERM), **cols))
    df["seed"] = SEED
    df.to_csv(OUT / "seng_cos_vu_null.csv", index=False)
    (OUT / "seng_cos_vu_null_report.txt").write_text("\n".join(lines) + "\n", encoding="utf-8")
    print("\n".join(lines))


if __name__ == "__main__":
    main()
