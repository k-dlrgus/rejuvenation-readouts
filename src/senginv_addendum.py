"""SENGINV addendum: counting and listing only.

Items 2-5 of the addendum task: control cells by library, NT cells with no
targeting-guide reads, read depth per group, and genes detected per group with
the three-way overlap. `.X` is streamed in row blocks and never materialised.

No distance, no statistic on a fitted axis, no pre-registration, no commit.
Item 1 (the earlier-SENG audit) is assembled by `senginv_addendum_item1.py`
and the two outputs are concatenated into ADDENDUM_SENGINV.md.
"""
import json
import re
from pathlib import Path

import h5py
import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
H5 = ROOT / "data" / "raw" / "seng" / "exp3_merged_not_normalized.h5ad"
OUT = ROOT / "results" / "senginv"
SEED = 20260914  # the seed the SENG task fixed; reused for the matched subsamples
BLOCK = 2000
PCTS = [0, 1, 5, 10, 25, 50, 75, 90, 95, 99, 100]

lines = []


def say(s=""):
    print(s)
    lines.append(s)


def dec(x):
    return x.decode() if isinstance(x, bytes) else x


def main():
    f = h5py.File(H5, "r")
    obs = f["obs"]
    cols = [dec(c) for c in obs.attrs["column-order"]]
    X = f["X"]
    n_cells, n_genes = (int(v) for v in X.attrs["shape"])
    indptr = X["indptr"][:]
    idx_ds, dat_ds = X["indices"], X["data"]

    oidx = np.array([dec(x) for x in obs[dec(obs.attrs["_index"])][:]])
    sfx = np.array([re.sub(r"^[ACGT]+", "", s) or "<none>" for s in oidx])
    bgroups = ["-1", "_A", "_I"]
    bmask = {g: sfx == g for g in bgroups}

    pd_cols = ["PD14", "PD26", "PD32"]
    pdm = {c: obs[c][:].astype(bool) for c in pd_cols}
    wt = obs["WT"][:].astype(float) > 0
    nt_cols = ["NT_1", "NT_2", "NT_3", "NT_4", "NT_5", "NT_6"]
    ntv = {c: obs[c][:].astype(float) for c in nt_cols}
    nt_any = np.zeros(n_cells, dtype=bool)
    for c in nt_cols:
        nt_any |= ntv[c] > 0
    gp = re.compile(r"^(.+)_(CRA|CRI)(\d+)$")
    guide_cols = [c for c in cols if gp.match(c)]
    tg_sum = np.zeros(n_cells, dtype=np.float64)
    tg_any = np.zeros(n_cells, dtype=bool)
    tg_ncol = np.zeros(n_cells, dtype=np.int32)
    for c in guide_cols:
        v = obs[c][:].astype(np.float64)
        tg_sum += v
        tg_any |= v > 0
        tg_ncol += (v > 0).astype(np.int32)
    n_counts = obs["n_counts"][:].astype(np.float64)
    n_genes_obs = obs["n_genes"][:].astype(np.float64)

    say("# SENGINV addendum - counting and listing only")
    say()
    say("Counts and lists. No distance, no fitted axis, no statistic that needs a")
    say("null, no pre-registration, nothing committed. Item 1 is in the companion")
    say("section below this one. Ambiguities are recorded, not resolved.")
    say()
    say("Inputs: the same untracked `data/raw/seng/exp3_merged_not_normalized.h5ad`,")
    say("md5 `bfccaa94a916c22e4b0793913a04128b`, verified in the SENGINV run. `.X` was")
    say("streamed in {}-row blocks through its CSR `indptr`; the matrix was never".format(BLOCK))
    say("loaded. Script: `src/senginv_addendum.py`.")
    say()
    say("Definitions used throughout, stated once:")
    say()
    say("- **barcode group**: the suffix on the `.obs` index string, one of `-1`,")
    say("  `_A`, `_I`. The file never says these are libraries; the task calls them")
    say("  libraries and that reading is adopted as the task's, not the file's.")
    say("- **NT-positive**: at least one of `NT_1`..`NT_6` is > 0.")
    say("- **targeting-guide reads**: the sum over the {} `<TARGET>_<CRA|CRI><k>`".format(
        len(guide_cols)))
    say("  columns. `NT_*` columns are excluded from that sum by construction, which")
    say("  is the definitional choice item 3 turns on.")
    say("- **WT cell**: `WT` > 0.")
    say("- **gene detected in a group**: at least one cell of that group has a")
    say("  nonzero entry for it in `.X`.")
    say()
    say("Group sizes, fixed here and referred to below:")
    say()
    say("| group | n cells |")
    say("|---|---|")
    for g in bgroups:
        say("| barcode `{}` | {} |".format(g, int(bmask[g].sum())))
    for c in pd_cols:
        say("| `{}` | {} |".format(c, int(pdm[c].sum())))
    say("| all | {} |".format(n_cells))
    say()

    # ---- 2 ------------------------------------------------------------------
    say("## 2. Control cells by library (barcode group)")
    say()
    say("| barcode group | n cells | NT-positive | WT>0 | NT-positive and WT>0 | "
        "neither |")
    say("|---|---|---|---|---|---|")
    for g in bgroups:
        m = bmask[g]
        say("| `{}` | {} | {} | {} | {} | {} |".format(
            g, int(m.sum()), int((m & nt_any).sum()), int((m & wt).sum()),
            int((m & nt_any & wt).sum()), int((m & ~nt_any & ~wt).sum())))
    say("| **total** | {} | {} | {} | {} | {} |".format(
        n_cells, int(nt_any.sum()), int(wt.sum()), int((nt_any & wt).sum()),
        int((~nt_any & ~wt).sum())))
    say()
    say("Per `NT_*` column, by barcode group:")
    say()
    say("| column | " + " | ".join("`{}`".format(g) for g in bgroups) + " | total |")
    say("|---|---|---|---|---|")
    for c in nt_cols:
        m = ntv[c] > 0
        say("| `{}` | ".format(c)
            + " | ".join(str(int((m & bmask[g]).sum())) for g in bgroups)
            + " | {} |".format(int(m.sum())) + "")
    say("| **any `NT_*`** | "
        + " | ".join(str(int((nt_any & bmask[g]).sum())) for g in bgroups)
        + " | {} |".format(int(nt_any.sum())))
    say()
    say("NT UMI totals, by barcode group (sum over the six `NT_*` columns):")
    say()
    ntsum = np.zeros(n_cells, dtype=np.float64)
    for c in nt_cols:
        ntsum += ntv[c]
    say("| barcode group | total NT UMI | NT UMI in NT-positive cells: "
        "min / median / max |")
    say("|---|---|---|")
    for g in bgroups:
        m = bmask[g] & nt_any
        if m.sum():
            say("| `{}` | {:.0f} | {:g} / {:g} / {:g} |".format(
                g, ntsum[bmask[g]].sum(), ntsum[m].min(), np.median(ntsum[m]),
                ntsum[m].max()))
        else:
            say("| `{}` | {:.0f} | no NT-positive cells |".format(
                g, ntsum[bmask[g]].sum()))
    say()
    say("WT and the PD arms against barcode group (WT is already known to coincide")
    say("with the union of the PD arms; repeated here because item 2 asks for it):")
    say()
    say("| barcode group | `WT`>0 | `PD14` | `PD26` | `PD32` |")
    say("|---|---|---|---|---|")
    for g in bgroups:
        m = bmask[g]
        say("| `{}` | {} | {} | {} | {} |".format(
            g, int((m & wt).sum()), *[int((m & pdm[c]).sum()) for c in pd_cols]))
    say()

    # ---- 3 ------------------------------------------------------------------
    say("## 3. NT-positive cells with zero targeting-guide reads")
    say()
    nt_clean = nt_any & ~tg_any
    say("- NT-positive cells: {}".format(int(nt_any.sum())))
    say("- of those, **{} have zero reads in all {} targeting-guide columns**".format(
        int(nt_clean.sum()), len(guide_cols)))
    say("  and **{} have at least one targeting-guide read**.".format(
        int((nt_any & tg_any).sum())))
    say()
    say("| barcode group | NT-positive | zero targeting-guide reads | "
        "at least one | % zero |")
    say("|---|---|---|---|---|")
    for g in bgroups:
        m = bmask[g] & nt_any
        z = int((m & ~tg_any).sum())
        pct = "{:.1f}%".format(100.0 * z / m.sum()) if m.sum() else "n/a"
        say("| `{}` | {} | {} | {} | {} |".format(
            g, int(m.sum()), z, int((m & tg_any).sum()), pct))
    say("| **total** | {} | {} | {} | {:.1f}% |".format(
        int(nt_any.sum()), int(nt_clean.sum()), int((nt_any & tg_any).sum()),
        100.0 * nt_clean.sum() / nt_any.sum()))
    say()
    say("For the NT-positive cells that do carry targeting-guide reads, how many")
    say("targeting-guide columns and how much targeting UMI they carry:")
    say()
    say("| barcode group | n such cells | n targeting cols: min/median/max | "
        "targeting UMI: min/median/max | NT UMI: median |")
    say("|---|---|---|---|---|")
    for g in bgroups:
        m = bmask[g] & nt_any & tg_any
        if not m.sum():
            say("| `{}` | 0 | - | - | - |".format(g))
            continue
        say("| `{}` | {} | {}/{:g}/{} | {:g}/{:g}/{:g} | {:g} |".format(
            g, int(m.sum()), tg_ncol[m].min(), np.median(tg_ncol[m]), tg_ncol[m].max(),
            tg_sum[m].min(), np.median(tg_sum[m]), tg_sum[m].max(),
            np.median(ntsum[m])))
    say()
    say("**Recorded as an ambiguity.** Whether an NT-positive cell carrying")
    say("targeting-guide reads is an NT control, a perturbed cell, a doublet or")
    say("ambient contamination is not something the object states, and no threshold")
    say("on NT UMI or targeting UMI is recorded in it. The split above uses >0 on")
    say("both sides because that is the only cut the file supports without a choice")
    say("being imported. A different cut would move these numbers.")
    say()

    # ---- 4 ------------------------------------------------------------------
    say("## 4. Read depth by group")
    say()
    say("From the stored `.obs` QC columns `n_counts` and `n_genes`, not recomputed.")
    say("Percentiles are linear-interpolation `numpy.percentile`.")
    say()
    groups = [("barcode `{}`".format(g), bmask[g]) for g in bgroups] \
        + [("`{}`".format(c), pdm[c]) for c in pd_cols] \
        + [("all cells", np.ones(n_cells, dtype=bool))]
    for label, col in [("n_counts", n_counts), ("n_genes", n_genes_obs)]:
        say("### `{}`".format(label))
        say()
        say("| group | n | mean | " + " | ".join("p{}".format(p) for p in PCTS) + " |")
        say("|---" * (len(PCTS) + 3) + "|")
        for gl, m in groups:
            v = col[m]
            q = np.percentile(v, PCTS)
            say("| {} | {} | {:.0f} | ".format(gl, int(m.sum()), v.mean())
                + " | ".join("{:.0f}".format(x) for x in q) + " |")
        say()
    say("The two QC columns are checked against `.X` itself, exactly and for every")
    say("row, in the same streaming pass as section 5. The result is reported after")
    say("section 5 because it needs that pass.")
    say()

    # ---- 5: stream .X -------------------------------------------------------
    say("## 5. Genes detected per group, and the three-way barcode overlap")
    say()
    rng = np.random.default_rng(SEED)
    n_match_b = min(int(bmask[g].sum()) for g in bgroups)
    n_match_p = min(int(pdm[c].sum()) for c in pd_cols)
    sub_b, sub_p = {}, {}
    for g in bgroups:
        ii = np.flatnonzero(bmask[g])
        pick = rng.choice(ii, size=n_match_b, replace=False)
        m = np.zeros(n_cells, dtype=bool)
        m[pick] = True
        sub_b[g] = m
    for c in pd_cols:
        ii = np.flatnonzero(pdm[c])
        pick = rng.choice(ii, size=n_match_p, replace=False)
        m = np.zeros(n_cells, dtype=bool)
        m[pick] = True
        sub_p[c] = m

    acc_names, acc_masks = [], []
    for g in bgroups:
        acc_names.append("barcode {}".format(g))
        acc_masks.append(bmask[g])
    for c in pd_cols:
        acc_names.append(c)
        acc_masks.append(pdm[c])
    for g in bgroups:
        acc_names.append("barcode {} matched n={}".format(g, n_match_b))
        acc_masks.append(sub_b[g])
    for c in pd_cols:
        acc_names.append("{} matched n={}".format(c, n_match_p))
        acc_masks.append(sub_p[c])
    acc = np.zeros((len(acc_masks), n_genes), dtype=np.int32)
    xrowsum = np.zeros(n_cells, dtype=np.float64)
    xrownnz = np.zeros(n_cells, dtype=np.int64)

    nblocks = 0
    for a in range(0, n_cells, BLOCK):
        b = min(a + BLOCK, n_cells)
        lo, hi = int(indptr[a]), int(indptr[b])
        if hi == lo:
            continue
        ind = idx_ds[lo:hi]
        rnnz = np.diff(indptr[a:b + 1]).astype(np.int64)
        dvals = dat_ds[lo:hi]
        offs = (indptr[a:b] - lo).astype(np.int64)
        nonempty = rnnz > 0
        rs = np.zeros(b - a, dtype=np.float64)
        if nonempty.any():
            rs[nonempty] = np.add.reduceat(dvals.astype(np.float64),
                                           offs[nonempty])
        xrowsum[a:b] = rs
        xrownnz[a:b] = rnnz
        for k, m in enumerate(acc_masks):
            sel = m[a:b]
            if not sel.any():
                continue
            keep = np.repeat(sel, rnnz)
            if keep.any():
                acc[k] += np.bincount(ind[keep], minlength=n_genes).astype(np.int32)
        nblocks += 1
    say("- streamed {} row blocks of at most {} rows; peak slice read was the".format(
        nblocks, BLOCK))
    say("  `indices` and `data` range of one block, never the whole {:,}-element".format(
        int(indptr[-1])))
    say("  arrays.")
    say()
    say("### Do the stored QC columns describe this `.X`?")
    say()
    say("Exact, every row, from the stream above.")
    say()
    dcnt = xrowsum - n_counts
    dgen = xrownnz - n_genes_obs
    say("| comparison | n rows equal | n rows differing | min diff | max diff |")
    say("|---|---|---|---|---|")
    say("| row sum of `.X` vs stored `n_counts` | {} | {} | {:g} | {:g} |".format(
        int((dcnt == 0).sum()), int((dcnt != 0).sum()), dcnt.min(), dcnt.max()))
    say("| row nnz of `.X` vs stored `n_genes` | {} | {} | {:g} | {:g} |".format(
        int((dgen == 0).sum()), int((dgen != 0).sum()), dgen.min(), dgen.max()))
    say()
    say("| barcode group | n rows where `n_counts` differs | where `n_genes` differs |")
    say("|---|---|---|")
    for g in bgroups:
        m = bmask[g]
        say("| `{}` | {} | {} |".format(g, int((m & (dcnt != 0)).sum()),
                                        int((m & (dgen != 0)).sum())))
    say()
    say("The shortfall against candidate sums the file does hold, as a count of")
    say("exact matches. This is arithmetic, not an explanation:")
    say()
    defc = n_counts - xrowsum
    defg = n_genes_obs - xrownnz
    mm = defc != 0
    cand = [("NT UMI total", ntsum), ("targeting-guide UMI total", tg_sum),
            ("NT + targeting UMI", ntsum + tg_sum)]
    say("| candidate | rows where shortfall equals it (of {} mismatched) |".format(
        int(mm.sum())))
    say("|---|---|")
    for nm, v in cand:
        say("| {} | {} |".format(nm, int((mm & (defc == v)).sum())))
    say()
    say("| candidate | rows where the `n_genes` shortfall equals it |")
    say("|---|---|")
    say("| n nonzero `NT_*` columns | {} |".format(
        int((mm & (defg == sum((ntv[c] > 0).astype(np.int64) for c in nt_cols))).sum())))
    say("| n nonzero targeting-guide columns | {} |".format(
        int((mm & (defg == tg_ncol)).sum())))
    say("| n nonzero guide columns of either kind | {} |".format(
        int((mm & (defg == tg_ncol + sum((ntv[c] > 0).astype(np.int64)
                                         for c in nt_cols))).sum())))
    say()
    say("| shortfall in `n_counts` | n rows |   | shortfall in `n_genes` | n rows |")
    say("|---|---|---|---|---|")
    vc1 = pd.Series(defc[mm]).value_counts().sort_index()
    vc2 = pd.Series(defg[mm]).value_counts().sort_index()
    k1, k2 = list(vc1.items()), list(vc2.items())
    for i in range(max(len(k1), len(k2))):
        lft = "{:g} | {}".format(*k1[i]) if i < len(k1) else " | "
        rgt = "{:g} | {}".format(*k2[i]) if i < len(k2) else " | "
        say("| {} |   | {} |".format(lft, rgt))
    say()
    say("**Recorded as an ambiguity.** Where the two disagree the stored value")
    say("exceeds the matrix on {} rows and falls short on {}, by at most {:g} counts".format(
        int((dcnt < 0).sum()), int((dcnt > 0).sum()), max(abs(dcnt.min()), abs(dcnt.max()))))
    say("and {:g} genes. The object does not record what matrix `n_counts` and".format(
        max(abs(dgen.min()), abs(dgen.max()))))
    say("`n_genes` were computed on, so whether the difference is the guide UMI")
    say("columns, a feature set that was later trimmed, or something else cannot be")
    say("read off the file. Section 4 uses the stored columns as the task asked; this")
    say("check is reported so that choice is visible.")
    say()
    say("Genes detected, with the cell count beside every number because group size")
    say("sets the ceiling on how many genes can be detected at all:")
    say()
    say("| group | n cells | genes in >=1 cell | >=3 cells | >=10 cells | "
        ">=1% of cells | total UMI-bearing entries |")
    say("|---|---|---|---|---|---|---|")
    rows = []
    for k, nm in enumerate(acc_names):
        m = acc_masks[k]
        n = int(m.sum())
        c1 = int((acc[k] >= 1).sum())
        rows.append(dict(group=nm, n_cells=n, genes_ge1=c1,
                         genes_ge3=int((acc[k] >= 3).sum()),
                         genes_ge10=int((acc[k] >= 10).sum()),
                         genes_ge1pct=int((acc[k] >= max(1, round(0.01 * n))).sum()),
                         nnz_entries=int(acc[k].sum())))
        say("| {} | {} | {} | {} | {} | {} | {} |".format(
            nm, n, c1, rows[-1]["genes_ge3"], rows[-1]["genes_ge10"],
            rows[-1]["genes_ge1pct"], rows[-1]["nnz_entries"]))
    pd.DataFrame(rows).to_csv(OUT / "addendum_genes_detected.csv", index=False)
    say()
    say("- written: `results/senginv/addendum_genes_detected.csv`")
    say("- the four `matched n=` rows are a size-matched diagnostic, not a requested")
    say("  result: each is one subsample without replacement at seed {}, drawn to".format(SEED))
    say("  the smallest group of its family ({} cells for the barcode groups, {} for".format(
        n_match_b, n_match_p))
    say("  the PD arms). One draw is not a null distribution and no interval is given.")
    say()
    say("### Three-way overlap between the barcode groups")
    say()
    det = {g: acc[acc_names.index("barcode {}".format(g))] >= 1 for g in bgroups}
    A, B, C = det["-1"], det["_A"], det["_I"]
    say("Detection is 'nonzero in at least one cell of the group'.")
    say()
    say("| region | n genes |")
    say("|---|---|")
    cells = [
        ("in all three (`-1` and `_A` and `_I`)", A & B & C),
        ("`-1` and `_A` only", A & B & ~C),
        ("`-1` and `_I` only", A & ~B & C),
        ("`_A` and `_I` only", ~A & B & C),
        ("`-1` only", A & ~B & ~C),
        ("`_A` only", ~A & B & ~C),
        ("`_I` only", ~A & ~B & C),
        ("in none of the three", ~A & ~B & ~C),
    ]
    for lbl, m in cells:
        say("| {} | {} |".format(lbl, int(m.sum())))
    say("| **all `.var` genes** | {} |".format(n_genes))
    say()
    say("| | `-1` | `_A` | `_I` |")
    say("|---|---|---|---|")
    for g1 in bgroups:
        say("| `{}` | ".format(g1)
            + " | ".join(str(int((det[g1] & det[g2]).sum())) for g2 in bgroups) + " |")
    say()
    say("Jaccard between barcode groups, on the detected-gene sets:")
    say()
    say("| pair | intersection | union | Jaccard |")
    say("|---|---|---|---|")
    for i, g1 in enumerate(bgroups):
        for g2 in bgroups[i + 1:]:
            inter = int((det[g1] & det[g2]).sum())
            uni = int((det[g1] | det[g2]).sum())
            say("| `{}` vs `{}` | {} | {} | {:.4f} |".format(g1, g2, inter, uni,
                                                             inter / uni))
    say()
    say("Same overlap on the size-matched subsamples (diagnostic, one draw, seed "
        "{}):".format(SEED))
    say()
    detm = {g: acc[acc_names.index("barcode {} matched n={}".format(g, n_match_b))] >= 1
            for g in bgroups}
    Am, Bm, Cm = detm["-1"], detm["_A"], detm["_I"]
    say("| region | n genes |")
    say("|---|---|")
    for lbl, m in [("in all three", Am & Bm & Cm),
                   ("`-1` only", Am & ~Bm & ~Cm),
                   ("`_A` only", ~Am & Bm & ~Cm),
                   ("`_I` only", ~Am & ~Bm & Cm),
                   ("in none", ~Am & ~Bm & ~Cm)]:
        say("| {} | {} |".format(lbl, int(m.sum())))
    say()
    say("### Genes detected per PD arm, and their overlap")
    say()
    detp = {c: acc[acc_names.index(c)] >= 1 for c in pd_cols}
    say("| region | n genes |")
    say("|---|---|")
    P1, P2, P3 = detp["PD14"], detp["PD26"], detp["PD32"]
    for lbl, m in [("in all three PD arms", P1 & P2 & P3),
                   ("`PD14` only", P1 & ~P2 & ~P3),
                   ("`PD26` only", ~P1 & P2 & ~P3),
                   ("`PD32` only", ~P1 & ~P2 & P3),
                   ("in no PD arm", ~P1 & ~P2 & ~P3)]:
        say("| {} | {} |".format(lbl, int(m.sum())))
    say()
    say("| | `PD14` | `PD26` | `PD32` |")
    say("|---|---|---|---|")
    for c1 in pd_cols:
        say("| `{}` | ".format(c1)
            + " | ".join(str(int((detp[c1] & detp[c2]).sum())) for c2 in pd_cols) + " |")
    say()
    say("**Recorded as an ambiguity, and it is the one that governs this whole")
    say("section.** The three barcode groups hold {}, {} and {} cells and the three".format(
        *[int(bmask[g].sum()) for g in bgroups]))
    say("PD arms hold {}, {} and {}. Detection at a >=1-cell threshold is a".format(
        *[int(pdm[c].sum()) for c in pd_cols]))
    say("saturating function of cell count and of depth per cell, so a larger group")
    say("detects more genes whether or not anything biological differs. The matched")
    say("rows above hold cell count fixed but not depth per cell, and section 4 shows")
    say("the depth distributions differ between groups. No comparison here is")
    say("corrected for either, and none of these counts is a result about the groups.")
    say()
    say("Genes detected in no cell at all: {} of {} `.var` columns.".format(
        int((acc[:3].sum(axis=0) == 0).sum()), n_genes))
    say()

    meta = dict(seed=SEED, block=BLOCK, n_cells=n_cells, n_genes=n_genes,
                nnz=int(indptr[-1]), n_guide_cols=len(guide_cols),
                n_nt_positive=int(nt_any.sum()),
                n_nt_positive_no_targeting=int(nt_clean.sum()),
                n_match_barcode=n_match_b, n_match_pd=n_match_p,
                n_rows_qc_mismatch=int((xrowsum != n_counts).sum()),
                qc_count_deficit_max=float((n_counts - xrowsum).max()),
                qc_gene_deficit_max=float((n_genes_obs - xrownnz).max()))
    (OUT / "addendum_meta.json").write_text(json.dumps(meta, indent=2))
    (OUT / "_addendum_body.md").write_text("\n".join(lines) + "\n", encoding="utf-8")
    f.close()
    print("\nwrote results/senginv/_addendum_body.md")


if __name__ == "__main__":
    main()
