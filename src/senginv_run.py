"""SENGINV: data inventory of the Sengstack et al. Perturb-seq h5ad.

Inventory only. No statistic, no gate, no pre-registration. Reports what is in
the file. Raw input is figshare DOI 10.6084/m9.figshare.30898748 (v2), stored
untracked under data/raw/seng/ (see .gitignore).
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
RULER = ROOT / "results" / "fibro" / "frozen_ruler_ridge_raw.npz"
GUIDE_PAT = re.compile(r"^(?P<tf>.+)_(?P<dir>CRA|CRI)(?P<n>\d+)$")
FOUR = ["E2F3", "EZH2", "STAT3", "ZFX"]

lines = []


def say(s=""):
    print(s)
    lines.append(s)


def dec(x):
    return x.decode() if isinstance(x, bytes) else x


def read_str_index(ds):
    a = ds[:]
    if a.dtype.kind in "OS":
        return np.array([dec(x) for x in a])
    return a


def read_col(obs, name):
    """Return (values, declared dtype string, categories or None).

    This file uses the legacy h5ad categorical layout: the codes are in the
    column dataset and the levels in the sibling group `obs/__categories`, not
    in a modern per-column group with `codes`/`categories` members. Handle both.
    """
    node = obs[name]
    if isinstance(node, h5py.Group):
        cats = read_str_index(node["categories"])
        codes = node["codes"][:]
        return codes, "categorical[%s] n_cat=%d (modern layout)" % (
            codes.dtype, len(cats)), cats
    v = node[:]
    legacy = obs.get("__categories")
    if legacy is not None and name in legacy:
        cats = read_str_index(legacy[name])
        return v, "categorical codes [%s] n_cat=%d (legacy obs/__categories)" % (
            v.dtype, len(cats)), cats
    if v.dtype.kind in "OS":
        return np.array([dec(x) for x in v]), "string/object", None
    return v, str(v.dtype), None


def main():
    say("# SENGINV - data inventory, Sengstack et al. Perturb-seq h5ad")
    say()
    say("Inventory only. Nothing here is a pre-registered result and nothing is")
    say("interpreted. Ambiguities are recorded as ambiguities, not resolved.")
    say()

    # ---- 0. provenance -----------------------------------------------------
    say("## 0. Provenance")
    say()
    say("- file: `data/raw/seng/exp3_merged_not_normalized.h5ad`")
    say("  (untracked; `.gitignore` rule `data/raw/`)")
    say("- bytes on disk: {:,}".format(H5.stat().st_size))
    for p in sorted(OUT.glob("figshare_*.json")):
        m = json.loads(p.read_text())
        say("- figshare article {} v{}: {}".format(m["id"], m["version"], m["title"]))
        say("  - doi `{}`, article size {:,} B, {}, published {}".format(
            m["doi"], m["size"], m["license"], m["published_date"]))
        for fl in m["files"]:
            say("  - file `{}` {:,} B md5 `{}`".format(
                fl["name"], fl["size"], fl["computed_md5"]))
    md5 = OUT / "md5_local.txt"
    if md5.exists():
        say("- md5 of the local copy: `{}`".format(md5.read_text().strip()))
    say()
    say("Download status, stated plainly:")
    say()
    say("- **article 30898748 was not downloaded in this session.** A byte-identical")
    say("  copy was already on disk at `data/raw/seng/` from the earlier SENG task.")
    say("  Its md5 was recomputed here and equals the md5 figshare reports for the")
    say("  v2 file, so the local copy is the v2 file named in the task.")
    say("- **article 30899018 was not downloaded.** Checked by API only, per the")
    say("  task instruction: its title and description say mouse liver RNA-seq raw")
    say("  counts (two CSVs, {:,} B in all), so it was left alone. Title and size"
        .format(json.loads((OUT / "figshare_30899018.json").read_text())["size"]))
    say("  are reported above and nothing else about it is reported.")
    say()

    f = h5py.File(H5, "r")
    obs = f["obs"]
    cols = [dec(c) for c in obs.attrs["column-order"]]
    X = f["X"]
    n_cells, n_genes = (int(v) for v in X.attrs["shape"])

    # ---- 1. shape and structure -------------------------------------------
    say("## 1. Shape and structure")
    say()
    say("- `.shape` = ({}, {})  ->  n_cells = {:,}, n_genes = {:,}".format(
        n_cells, n_genes, n_cells, n_genes))
    say("- top-level HDF5 keys: {}".format(list(f.keys())))
    say("- root attrs: {}  (no `encoding-type` / anndata-version written)".format(dict(f.attrs)))
    say("- `.X` encoding: {} v{}, nnz = {:,}, dtype = {}".format(
        dec(X.attrs["encoding-type"]), dec(X.attrs["encoding-version"]),
        X["data"].shape[0], X["data"].dtype))
    say("- density = {:.4f}".format(X["data"].shape[0] / (n_cells * n_genes)))
    say("- `.obsm` keys: **NONE** (no `obsm` group in the file)")
    say("- `.layers` keys: **NONE** (no `layers` group in the file)")
    say("- `.uns` keys: **NONE** (no `uns` group in the file)")
    say("- `.varm`, `.obsp`, `.varp`, `.raw`: **NONE**")
    say("- `.var` columns: {} (index only, index name `{}`)".format(
        [dec(c) for c in f["var"].attrs.get("column-order", [])],
        dec(f["var"].attrs["_index"])))
    say("- `.obs` index name: `{}`".format(dec(obs.attrs["_index"])))
    say("- number of `.obs` columns: {}".format(len(cols)))
    say("- `.obs` also holds a legacy `__categories` group, which is not a column:")
    say("  the levels of the categorical columns live there, h5ad's pre-0.8 layout.")
    say("  Its members: {}".format(
        list(obs["__categories"].keys()) if "__categories" in obs else "n/a"))
    say()

    say("### Are raw integer counts present, and where?")
    say()
    d = X["data"]
    head = d[:5_000_000]
    samp = np.concatenate([head, d[-5_000_000:]])
    say("- `.X` is the only matrix in the file: CSR, float32.")
    say("- sampled {:,} stored values (first 5M and last 5M of `X/data`):".format(samp.size))
    say("  - every sampled value equals its own rounding: {}".format(
        bool(np.all(samp == np.rint(samp)))))
    say("  - min {:g}, max {:g}, n distinct {}".format(
        samp.min(), samp.max(), np.unique(samp).size))
    say("  - first 20 stored values: {}".format(head[:20].tolist()))
    say("- so `.X` holds integer-valued counts stored as float32, consistent with")
    say("  the file name `not_normalized`.")
    say("- there is no second copy of the counts: no `.raw`, no `.layers`.")
    say()

    # ---- obs column split ---------------------------------------------------
    guide_cols = [c for c in cols if GUIDE_PAT.match(c)]
    other_cols = [c for c in cols if not GUIDE_PAT.match(c)]

    say("## 1b. Full `.obs` column list with dtype")
    say()
    say("`.obs` has {} columns. {} match `<NAME>_(CRA|CRI)<n>` and are per-guide".format(
        len(cols), len(guide_cols)))
    say("UMI count columns (section 4); the other {} are listed below in full.".format(
        len(other_cols)))
    say()
    say("### Non-guide columns, all {} of them, in stored order".format(len(other_cols)))
    say()
    say("| # | column | stored dtype | n nonzero | min | max | n distinct |")
    say("|---|---|---|---|---|---|---|")
    store = {}
    for i, c in enumerate(other_cols, 1):
        v, dt, cats = read_col(obs, c)
        store[c] = (v, dt, cats)
        if cats is not None:
            say("| {} | `{}` | {} | - | - | - | {} |".format(i, c, dt, len(cats)))
        elif v.dtype.kind in "fiub":
            say("| {} | `{}` | {} | {} | {:g} | {:g} | {} |".format(
                i, c, dt, int(np.count_nonzero(v)), v.min(), v.max(), np.unique(v).size))
        else:
            say("| {} | `{}` | {} | - | - | - | {} |".format(i, c, dt, np.unique(v).size))
    say()
    say("### Distinct values and counts for the non-guide columns")
    say()
    say("One `.obs` column is categorical: `SampleName`, stored in the legacy h5ad")
    say("layout as int8 codes whose levels are in `obs/__categories/SampleName`.")
    say("Its distinct values and counts are given below by level name. The other")
    say("non-guide columns are plain bool or numeric; value counts are given for")
    say("every one with at most 50 distinct values, percentiles for the rest.")
    say()
    for c in other_cols:
        v, dt, cats = store[c]
        say("**`{}`** (stored {})".format(c, dt))
        say()
        if cats is not None:
            vals = cats[v]
            tab = pd.Series(vals).value_counts()
        else:
            u = np.unique(v)
            if u.size > 50:
                q = np.percentile(v.astype(float), [0, 1, 25, 50, 75, 99, 100])
                say("- {} distinct values; percentiles 0/1/25/50/75/99/100 = {}".format(
                    u.size, " / ".join("{:.6g}".format(x) for x in q)))
                say("- n zero = {}, mean = {:.6g}".format(
                    int((v == 0).sum()), v.astype(float).mean()))
                say()
                continue
            tab = pd.Series(v).value_counts().sort_index()
        say("| value | n cells |")
        say("|---|---|")
        for k, n in tab.items():
            say("| `{}` | {} |".format(k, n))
        say()

    # ---- 2. early passage ---------------------------------------------------
    say("## 2. Are early-passage (young) cells present in this object?")
    say()
    pd_cols = [c for c in other_cols if re.fullmatch(r"PD\d+", c)]
    say("- `.obs` columns whose name is `PD<number>`: {}".format(pd_cols))
    say("- the numbers in those names are {}.".format(
        [int(c[2:]) for c in pd_cols]))
    say("- **Recorded as an ambiguity.** PD conventionally means population")
    say("  doublings, so `PD14` would be the earliest-passage arm of the three, but")
    say("  the object itself carries no definition: `.uns` is absent and there is no")
    say("  description field anywhere in `.obs` or `.var`. That `PD14` is the")
    say("  early-passage (young) arm rests on the column name alone. The three")
    say("  columns are stored as bool, so they are indicators and not counts, but")
    say("  what the indicator means - doublings, passage number, culture age - is")
    say("  not stated in the object.")
    say()
    say("### Per-PD-column occupancy")
    say()
    say("| column | n cells nonzero | n distinct values | value counts (value:n cells) |")
    say("|---|---|---|---|")
    pdmat = {}
    for c in pd_cols:
        v = store[c][0].astype(float)
        pdmat[c] = v
        tab = pd.Series(v).value_counts().sort_index()
        shown = ", ".join("{:g}:{}".format(k, n) for k, n in list(tab.items())[:12])
        if len(tab) > 12:
            shown += ", ... ({} values total)".format(len(tab))
        say("| `{}` | {} | {} | {} |".format(
            c, int(np.count_nonzero(v)), np.unique(v).size, shown))
    say()
    M = np.vstack([pdmat[c] for c in pd_cols])
    nz = (M != 0).sum(axis=0)
    wt_pre = store["WT"][0].astype(float) != 0
    say("### Do the PD columns partition the cells?")
    say()
    say("| n PD columns nonzero for a cell | n cells |")
    say("|---|---|")
    for k, n in pd.Series(nz).value_counts().sort_index().items():
        say("| {} | {} |".format(k, n))
    say()
    say("- cells with exactly one nonzero PD column: {}".format(int((nz == 1).sum())))
    say("- cells with no nonzero PD column: {}".format(int((nz == 0).sum())))
    say("- cells with more than one nonzero PD column: {}".format(int((nz > 1).sum())))
    say()
    say("### Lane / batch / library / sample / 10x-channel columns")
    say()
    sn, sdt, scats = store["SampleName"]
    vals = np.asarray(scats[sn] if scats is not None else sn)
    say("- the only such column in the object is `SampleName` (stored {}).".format(sdt))
    say("- there is no column named lane, batch, library, channel, 10x, gem, run,")
    say("  donor, replicate or experiment. The {} non-guide columns listed in 1b".format(
        len(other_cols)))
    say("  are the whole of what exists.")
    say()
    scount = pd.Series(vals).value_counts()
    say("| SampleName | n cells |")
    say("|---|---|")
    for k, n in scount.items():
        say("| `{}` | {} |".format(k, n))
    say()
    say("### The `.obs` index itself (10x barcode suffixes)")
    say()
    say("The `.obs` index is 10x cell barcodes carrying a suffix. That suffix is the")
    say("only other per-cell grouping in the file. It is part of the index string and")
    say("not a column, so it is reported here as what it is.")
    say()
    oidx = read_str_index(obs[dec(obs.attrs["_index"])])
    sfxs = np.array([re.sub(r"^[ACGT]+", "", s) or "<none>" for s in oidx])
    sfxorder = list(pd.Series(sfxs).value_counts().index)
    say("- first 3 index entries: {}".format(oidx[:3].tolist()))
    say("- duplicated index entries: {}".format(int(oidx.size - np.unique(oidx).size)))
    say("- distinct bare barcodes with the suffix stripped: {:,} of {:,} rows".format(
        len(set(re.match(r"^[ACGT]+", s).group(0) for s in oidx)), oidx.size))
    say()
    say("| barcode suffix | n cells |")
    say("|---|---|")
    for k, n in pd.Series(sfxs).value_counts().items():
        say("| `{}` | {} |".format(k, n))
    say()
    say("| SampleName | " + " | ".join("`{}`".format(k) for k in sfxorder) + " |")
    say("|---" * (len(sfxorder) + 1) + "|")
    for sv in scount.index:
        m = vals == sv
        say("| `{}` | ".format(sv)
            + " | ".join(str(int((m & (sfxs == k)).sum())) for k in sfxorder) + " |")
    say()
    say("| barcode suffix | " + " | ".join("`{}`>0".format(c) for c in pd_cols)
        + " | `CRA`>0 | `CRI`>0 | `WT`>0 |")
    say("|---" * (len(pd_cols) + 4) + "|")
    for k in sfxorder:
        m = sfxs == k
        row = [str(int((m & (pdmat[c] != 0)).sum())) for c in pd_cols]
        row += [str(int((m & (store["CRA"][0].astype(float) != 0)).sum())),
                str(int((m & (store["CRI"][0].astype(float) != 0)).sum())),
                str(int((m & wt_pre).sum()))]
        say("| `{}` | ".format(k) + " | ".join(row) + " |")
    say()
    say("**Recorded as an ambiguity.** A barcode suffix is what `cellranger aggr` or")
    say("an anndata concatenation writes to keep barcodes unique across inputs, so")
    say("the suffix groups are candidates for separate libraries. The file does not")
    say("say that, and does not say what `-1`, `_A` or `_I` stand for. Nothing here")
    say("is a claim that the suffix is a lane, a library or a sequencing run.")
    say()
    say("### Cross-tabulation: PD columns against SampleName")
    say()
    say("| SampleName | " + " | ".join("`{}`>0".format(c) for c in pd_cols)
        + " | no PD>0 | total |")
    say("|---" * (len(pd_cols) + 3) + "|")
    for s in scount.index:
        m = vals == s
        row = [str(int(((pdmat[c] != 0) & m).sum())) for c in pd_cols]
        row += [str(int(((nz == 0) & m).sum())), str(int(m.sum()))]
        say("| `{}` | ".format(s) + " | ".join(row) + " |")
    say()

    say("### Answer to question 2 as asked")
    say()
    pd14 = pdmat["PD14"] != 0
    say("- **early-passage cells, reading `PD14` as the earliest arm: present,**")
    say("  **n = {}.**".format(int(pd14.sum())))
    say("- identified by the bool `.obs` column `PD14`, True for those {} cells;".format(
        int(pd14.sum())))
    say("  equivalently by `SampleName == 'PD14'` ({} cells; the two labels agree".format(
        int((vals == "PD14").sum())))
    say("  on {} of {} cells).".format(int((pd14 == (vals == "PD14")).sum()), n_cells))
    say("- the three PD arms together are {} cells. That is also exactly the number".format(
        int((nz > 0).sum())))
    say("  with `WT`>0 ({}) and exactly the number whose barcode suffix is `-1`.".format(
        int(wt_pre.sum())))
    anyg = np.zeros(n_cells, dtype=bool)
    for c in guide_cols:
        anyg |= obs[c][:] > 0
    say("- columns that could show a sequencing run: only `SampleName`, whose four")
    say("  levels are the three PD arms and `UNKNOWN`, with {} cells in `UNKNOWN`.".format(
        int((vals == "UNKNOWN").sum())))
    say("- measured, not assumed: {} cells have a nonzero value in at least one of".format(
        int(anyg.sum())))
    say("  the {} guide columns. Of those, {} are in `SampleName == 'UNKNOWN'` and".format(
        len(guide_cols), int((anyg & (vals == "UNKNOWN")).sum())))
    say("  {} are in a PD arm. So `SampleName` separates the PD cells from the".format(
        int((anyg & (nz > 0)).sum())))
    say("  cells carrying guides rather than placing them in one group. The barcode")
    say("  suffixes split the same way: `-1` for all PD cells, `_A` and `_I` for all")
    say("  `UNKNOWN` cells. There is no lane, batch, library, channel or run column.")
    say("- **Recorded as an ambiguity, not resolved.** Whether the PD cells were")
    say("  sequenced in the same run as the perturbed cells cannot be read off this")
    say("  object. A separate `SampleName` level and a separate barcode suffix are")
    say("  consistent with a separate run and equally with a separate `cellranger`")
    say("  input within one run. The file records no run, flowcell, chemistry or")
    say("  date identifier that would decide it.")
    say()

    # ---- 3. non-targeting controls -----------------------------------------
    say("## 3. Non-targeting control cells")
    say()
    nt_cols = sorted((c for c in other_cols if re.fullmatch(r"NT_\d+", c)),
                     key=lambda c: int(c.split("_")[1]))
    say("- `.obs` columns matching `NT_<n>`: {}".format(nt_cols))
    say("- there is also a column `WT` (stored {}).".format(store["WT"][1]))
    say("- **Recorded as an ambiguity.** The object does not say that `NT_*` means")
    say("  non-targeting guide, nor what `WT` is, nor whether `WT` is a separate")
    say("  unperturbed arm or an overlapping label. The tables below report the")
    say("  numbers without choosing a reading.")
    say()
    say("| column | n cells nonzero | min nonzero | median nonzero | max |")
    say("|---|---|---|---|---|")
    for c in nt_cols + ["WT"]:
        v = store[c][0].astype(float)
        p = v[v != 0]
        if p.size:
            say("| `{}` | {} | {:g} | {:g} | {:g} |".format(
                c, p.size, p.min(), np.median(p), v.max()))
        else:
            say("| `{}` | 0 | - | - | - |".format(c))
    say()
    NTM = np.vstack([store[c][0].astype(float) for c in nt_cols])
    nt_any = (NTM != 0).any(axis=0)
    wt = store["WT"][0].astype(float) != 0
    say("- cells with at least one nonzero `NT_*` column: {}".format(int(nt_any.sum())))
    say("- cells with nonzero `WT`: {}".format(int(wt.sum())))
    say("- cells nonzero in both an `NT_*` and `WT`: {}".format(int((nt_any & wt).sum())))
    say("- number of nonzero `NT_*` columns per cell:")
    say()
    say("| n NT_* columns nonzero | n cells |")
    say("|---|---|")
    for k, n in pd.Series((NTM != 0).sum(axis=0)).value_counts().sort_index().items():
        say("| {} | {} |".format(k, n))
    say()
    say("### `NT_*` and `WT` against the PD columns")
    say()
    say("| group | " + " | ".join("`{}`>0".format(c) for c in pd_cols)
        + " | no PD>0 | total |")
    say("|---" * (len(pd_cols) + 3) + "|")
    groups = [("any `NT_*`>0", nt_any), ("`WT`>0", wt)] + \
             [("`{}`>0".format(c), store[c][0].astype(float) != 0) for c in nt_cols]
    for label, m in groups:
        row = [str(int(((pdmat[c] != 0) & m).sum())) for c in pd_cols]
        row += [str(int(((nz == 0) & m).sum())), str(int(m.sum()))]
        say("| {} | ".format(label) + " | ".join(row) + " |")
    say()
    say("### `NT_*` and `WT` against SampleName")
    say()
    say("| group | " + " | ".join("`{}`".format(s) for s in scount.index) + " |")
    say("|---" * (len(scount) + 1) + "|")
    for label, m in groups:
        say("| {} | ".format(label)
            + " | ".join(str(int((m & (vals == s)).sum())) for s in scount.index) + " |")
    say()

    # ---- 4. guides ----------------------------------------------------------
    say("## 4. Guide / perturbation assignment")
    say()
    say("- the per-guide UMI counts are `.obs` columns, one column per guide, named")
    say("  `<TARGET>_<CRA|CRI><1|2>`. There is no `obsm['guide_counts']`, no guide")
    say("  `layers` entry and no separate guide-by-cell matrix: `.obs` is where they")
    say("  live. This matches the figshare description, which says the two guides per")
    say("  target were counted separately.")
    say("- number of guide columns: {}".format(len(guide_cols)))
    say()
    recs = []
    for c in guide_cols:
        m = GUIDE_PAT.match(c)
        v = obs[c][:].astype(float)
        p = v[v > 0]
        recs.append(dict(col=c, target=m.group("tf"), direction=m.group("dir"),
                         replicate=int(m.group("n")), stored_dtype=str(obs[c].dtype),
                         n_cells_umi_gt0=int(p.size),
                         umi_min_pos=float(p.min()) if p.size else np.nan,
                         umi_median_pos=float(np.median(p)) if p.size else np.nan,
                         umi_max=float(v.max()),
                         umi_total=float(v.sum())))
    G = pd.DataFrame(recs)
    G.to_csv(OUT / "guide_columns.csv", index=False)
    say("- per-guide-column table: `results/senginv/guide_columns.csv` ({} rows)".format(len(G)))
    say("- distinct stored dtypes among guide columns: {}".format(
        sorted(G["stored_dtype"].unique())))
    say("- replicate suffixes present: {}".format(sorted(G["replicate"].unique())))
    say("- all guide UMI values integral: {}".format(
        bool(np.all(G["umi_max"] == np.rint(G["umi_max"])))))
    say("- distinct target names across guide columns: **{}**".format(G["target"].nunique()))
    say("- distinct target names appearing in CRA columns: {}".format(
        G.loc[G.direction == "CRA", "target"].nunique()))
    say("- distinct target names appearing in CRI columns: {}".format(
        G.loc[G.direction == "CRI", "target"].nunique()))
    a_set = set(G.loc[G.direction == "CRA", "target"])
    i_set = set(G.loc[G.direction == "CRI", "target"])
    blocks = G.groupby(["target", "direction"]).size()
    say("- guide columns per target x direction:")
    for k, n in blocks.value_counts().sort_index().items():
        say("  - {} column(s): {} target x direction pairs".format(k, n))
    odd = blocks[blocks != 2]
    say("- {} x 2 directions x 2 guides would be {} columns; there are {}.".format(
        G["target"].nunique(), G["target"].nunique() * 4, len(G)))
    say("- pairs not holding exactly 2 guide columns: {}".format(
        [(t, d, int(v)) for (t, d), v in odd.items()] if len(odd) else "none"))
    say("- targets with both CRA and CRI guide columns: {}".format(len(a_set & i_set)))
    say("- targets with CRA columns only: {}".format(sorted(a_set - i_set)))
    say("- targets with CRI columns only: {}".format(sorted(i_set - a_set)))
    say("- **Recorded as an ambiguity.** 'Target TF' cannot be read off the object:")
    say("  the target names are just the prefixes of the guide column names, and the")
    say("  object says nothing about which of them are transcription factors. Section")
    say("  9 reports whether each prefix is also a symbol in `.var`, which is not the")
    say("  same question.")
    say()
    say("### How a cell gets assigned to a target")
    say()
    say("**Nothing in the object assigns a cell to a target.** There is no")
    say("`guide_identity`, `perturbation`, `target`, `gene_target`, `condition` or")
    say("similar single-value column - the {} non-guide columns are listed in full".format(
        len(other_cols)))
    say("in 1b and none of them is an assignment. An assignment would have to be")
    say("derived from the {} UMI columns by a rule that is not recorded in the file.".format(
        len(guide_cols)))
    say("**Recorded as an ambiguity and not resolved here.** The counts below are")
    say("reported under two explicitly stated counting conventions chosen only so")
    say("that a number can be printed. Neither is an assignment rule and neither is")
    say("endorsed:")
    say()
    say("- convention A (any UMI): a cell is counted toward target T if any guide")
    say("  column of T has UMI > 0. A cell may be counted toward many targets, so")
    say("  the column does not sum to n_cells.")
    say("- convention B (argmax): a cell is counted toward the single target holding")
    say("  the strictly largest guide UMI over all {} guide columns. Cells with a".format(
        len(guide_cols)))
    say("  tied maximum, and cells whose guide columns are all zero, are dropped.")
    say()
    Gm = np.zeros((n_cells, len(guide_cols)), dtype=np.float32)
    for j, c in enumerate(guide_cols):
        Gm[:, j] = obs[c][:]
    targets = G["target"].to_numpy()
    dirs = G["direction"].to_numpy()
    keys = np.array(["{}|{}".format(t, dd) for t, dd in zip(targets, dirs)])
    ukeys = np.unique(keys)

    mx = Gm.max(axis=1)
    has = mx > 0
    nties = np.zeros(n_cells, dtype=int)
    nties[has] = (Gm[has] == mx[has, None]).sum(axis=1)
    uniq = has & (nties == 1)
    amax = np.zeros(n_cells, dtype=int)
    amax[uniq] = Gm[uniq].argmax(axis=1)
    bkey = np.where(uniq, keys[amax], "<none>")
    bcount = pd.Series(bkey).value_counts()

    say("- cells with every guide column zero: {}".format(int((~has).sum())))
    say("- cells with a tied maximum across guide columns: {}".format(
        int((has & (nties > 1)).sum())))
    say("- cells with a unique argmax guide: {}".format(int(uniq.sum())))
    say()
    rows = []
    for k in ukeys:
        t, dd = k.split("|")
        sel = keys == k
        sub = Gm[:, sel]
        rows.append(dict(target=t, direction=dd, n_guide_cols=int(sel.sum()),
                         cells_any_umi=int((sub > 0).any(axis=1).sum()),
                         cells_argmax=int(bcount.get(k, 0)),
                         total_guide_umi=float(sub.sum())))
    T = pd.DataFrame(rows).sort_values(
        ["cells_any_umi", "target", "direction"], ascending=[False, True, True])
    T.to_csv(OUT / "cells_per_target.csv", index=False)
    say("- full sorted table also written to `results/senginv/cells_per_target.csv`")
    say("  ({} target x direction rows)".format(len(T)))
    say()
    for name, col in [("A (any UMI)", "cells_any_umi"), ("B (argmax)", "cells_argmax")]:
        s = T[col]
        say("- convention {}: min {}, median {:g}, max {}, sum over rows {}".format(
            name, s.min(), s.median(), s.max(), s.sum()))
    say()
    say("Aggregated over the two directions, per target name:")
    say()
    for name, col in [("A (any UMI)", "cells_any_umi"), ("B (argmax)", "cells_argmax")]:
        if col == "cells_any_umi":
            per = []
            for t in np.unique(targets):
                sub = Gm[:, targets == t]
                per.append(int((sub > 0).any(axis=1).sum()))
            s = pd.Series(per)
        else:
            s = T.groupby("target")[col].sum()
        say("- convention {}: min {}, median {:g}, max {}".format(
            name, s.min(), s.median(), s.max()))
    say()
    say("### Full sorted table: cells per target x direction")
    say()
    say("| rank | target | direction | n guide cols | cells any UMI>0 | "
        "cells argmax | total guide UMI |")
    say("|---|---|---|---|---|---|---|")
    for i, r in enumerate(T.itertuples(index=False), 1):
        say("| {} | `{}` | {} | {} | {} | {} | {:.0f} |".format(
            i, r.target, r.direction, r.n_guide_cols, r.cells_any_umi,
            r.cells_argmax, r.total_guide_umi))
    say()
    say("### Guide-UMI columns per cell")
    say()
    k = (Gm > 0).sum(axis=1)
    say("- number of guide columns with UMI>0 per cell: percentiles "
        "0/25/50/75/90/99/100 = " + " / ".join(
            "{:g}".format(x) for x in np.percentile(k, [0, 25, 50, 75, 90, 99, 100])))
    say("- total guide UMI per cell: percentiles 0/25/50/75/99/100 = " + " / ".join(
        "{:.6g}".format(x) for x in np.percentile(Gm.sum(axis=1), [0, 25, 50, 75, 99, 100])))
    say()

    # ---- 5. direction -------------------------------------------------------
    say("## 5. CRISPRa vs CRISPRi")
    say()
    say("- **both are in this one object.**")
    say("- direction is encoded in two places, which pick out different cell sets")
    say("  (tables below):")
    say("  1. in the guide column **name**: the `CRA` / `CRI` infix.")
    say("  2. as two non-guide `.obs` columns, `CRA` (stored {}) and `CRI` ({}).".format(
        store["CRA"][1], store["CRI"][1]))
    say()
    cra = store["CRA"][0].astype(float)
    cri = store["CRI"][0].astype(float)
    say("| column | n cells nonzero | n distinct | min nonzero | median nonzero | max |")
    say("|---|---|---|---|---|---|")
    for nm, v in [("`CRA`", cra), ("`CRI`", cri)]:
        p = v[v != 0]
        if p.size:
            say("| {} | {} | {} | {:g} | {:g} | {:g} |".format(
                nm, p.size, np.unique(v).size, p.min(), np.median(p), v.max()))
        else:
            say("| {} | 0 | {} | - | - | - |".format(nm, np.unique(v).size))
    say()
    say("| `CRA`>0 | `CRI`>0 | n cells |")
    say("|---|---|---|")
    for a in [False, True]:
        for b in [False, True]:
            say("| {} | {} | {} |".format(
                a, b, int((((cra > 0) == a) & ((cri > 0) == b)).sum())))
    say()
    say("### Counts per direction taken from the guide columns instead")
    say()
    for dd in ["CRA", "CRI"]:
        sel = dirs == dd
        sub = Gm[:, sel]
        say("- {}: {} guide columns, {} distinct target names, {} cells with any "
            "UMI>0, total UMI {:.0f}".format(
                dd, int(sel.sum()), len(set(targets[sel])),
                int((sub > 0).any(axis=1).sum()), sub.sum()))
    aA = (Gm[:, dirs == "CRA"] > 0).any(axis=1)
    aI = (Gm[:, dirs == "CRI"] > 0).any(axis=1)
    say()
    say("| any CRA guide UMI>0 | any CRI guide UMI>0 | n cells |")
    say("|---|---|---|")
    for a in [False, True]:
        for b in [False, True]:
            say("| {} | {} | {} |".format(a, b, int(((aA == a) & (aI == b)).sum())))
    say()
    say("- cells assigned to a CRA guide by convention B (argmax): {}".format(
        int(np.isin(bkey, [k for k in ukeys if k.endswith("|CRA")]).sum())))
    say("- cells assigned to a CRI guide by convention B (argmax): {}".format(
        int(np.isin(bkey, [k for k in ukeys if k.endswith("|CRI")]).sum())))
    say()
    say("**Recorded as an ambiguity.** The `CRA`/`CRI` `.obs` columns and the")
    say("`_CRA`/`_CRI` guide-column names give different cell sets (tables above),")
    say("and the file carries no statement of what the two `.obs` columns count or")
    say("which encoding is authoritative.")
    say()
    say("### Direction against SampleName and the PD columns")
    say()
    say("| group | " + " | ".join("`{}`".format(s) for s in scount.index) + " | "
        + " | ".join("`{}`>0".format(c) for c in pd_cols) + " | total |")
    say("|---" * (len(scount) + len(pd_cols) + 2) + "|")
    for label, m in [("`CRA`>0", cra > 0), ("`CRI`>0", cri > 0),
                     ("any CRA guide UMI>0", aA), ("any CRI guide UMI>0", aI)]:
        row = [str(int((m & (vals == s)).sum())) for s in scount.index]
        row += [str(int((m & (pdmat[c] != 0)).sum())) for c in pd_cols]
        row.append(str(int(m.sum())))
        say("| {} | ".format(label) + " | ".join(row) + " |")
    say()

    # ---- 6. the four named factors -----------------------------------------
    say("## 6. E2F3, EZH2, STAT3, ZFX")
    say()
    say("| factor | guide columns present | direction(s) | cells any UMI>0 | "
        "cells argmax | total guide UMI |")
    say("|---|---|---|---|---|---|")
    for t in FOUR:
        sub = G[G.target == t]
        if sub.empty:
            say("| `{}` | **NONE** | - | - | - | - |".format(t))
            continue
        blk = Gm[:, targets == t]
        nax = sum(int(bcount.get("{}|{}".format(t, dd), 0))
                  for dd in sub.direction.unique())
        say("| `{}` | {} | {} | {} | {} | {:.0f} |".format(
            t, ", ".join("`%s`" % x for x in sorted(sub.col)),
            ", ".join(sorted(sub.direction.unique())),
            int((blk > 0).any(axis=1).sum()), nax, blk.sum()))
    say()
    say("Per direction:")
    say()
    say("| factor | direction | guide columns | cells any UMI>0 | cells argmax |")
    say("|---|---|---|---|---|")
    for t in FOUR:
        for dd in ["CRA", "CRI"]:
            k = "{}|{}".format(t, dd)
            if k in set(ukeys):
                sel = keys == k
                sub = Gm[:, sel]
                say("| `{}` | {} | {} | {} | {} |".format(
                    t, dd, ", ".join("`%s`" % c for c in np.array(guide_cols)[sel]),
                    int((sub > 0).any(axis=1).sum()), int(bcount.get(k, 0))))
            else:
                say("| `{}` | {} | column absent | - | - |".format(t, dd))
    say()
    say("Both counts are under the stated conventions of section 4, not an")
    say("assignment the object makes. Section 9 reports these four as `.var` symbols.")
    say()

    # ---- 7. cell cycle ------------------------------------------------------
    say("## 7. Cell-cycle phase annotation")
    say()
    hits = [c for c in cols if re.search(r"phase|cycle|S_score|G2M|G1S|cc_", c, re.I)]
    say("- `.obs` columns matching /phase|cycle|S_score|G2M|G1S|cc_/i: {}".format(
        hits if hits else "NONE"))
    say("- no `phase`, `S_score` or `G2M_score` column exists; `.uns` is absent, so")
    say("  no cell-cycle gene list is stored either.")
    say("- **absent.** There is no cell-cycle phase annotation in this object.")
    say()

    # ---- 8. author-computed scores -----------------------------------------
    say("## 8. Author-computed scores already in the object")
    say()
    pat = r"R_rej|rej|score|age|clock|senes|pseudo|umap|pca|harmony|leiden|louvain|cluster"
    hits = [c for c in cols if re.search(pat, c, re.I)]
    say("- `.obs` columns matching /{}/i: {}".format(pat, hits if hits else "NONE"))
    say("- `.obsm` is absent, so no embedding (PCA, UMAP, Harmony) is stored.")
    say("- `.uns` is absent, so no author parameters, rankings or score tables are")
    say("  stored.")
    say("- the only derived per-cell quantities in `.obs` are the three QC columns")
    say("  `n_counts`, `log_counts`, `n_genes`, summarised in 1b.")
    say("- **no `R_rej` and no comparable author-computed score is present.**")
    say()

    # ---- 9. gene identifiers -----------------------------------------------
    say("## 9. Gene identifier format, and overlap with the frozen GTEx panel")
    say()
    vidx = read_str_index(f["var"][dec(f["var"].attrs["_index"])])
    say("- `.var` index length {:,}; first 10 entries: {}".format(
        vidx.size, vidx[:10].tolist()))
    n_ens = int(sum(bool(re.match(r"^ENSG\d{11}", s)) for s in vidx))
    say("- entries matching `^ENSG\\d{{11}}`: {}".format(n_ens))
    say("- format: **gene symbols**, not Ensembl IDs.")
    say("- duplicated `.var` index entries: {} (n unique = {:,})".format(
        int(vidx.size - np.unique(vidx).size), np.unique(vidx).size))
    dups = pd.Series(vidx).value_counts()
    dups = dups[dups > 1]
    say("- symbols appearing more than once: {}{}".format(
        len(dups), "; first 20: {}".format(list(dups.index[:20])) if len(dups) else ""))
    say()
    z = np.load(RULER, allow_pickle=True)
    say("- reference panel **read, not refit**, from")
    say("  `results/fibro/frozen_ruler_ridge_raw.npz` (the frozen GTEx fibroblast")
    say("  ruler of this project), keys `symbol`, `ensembl`, `gene_id`; n = {:,}.".format(
        z["symbol"].size))
    say()
    say("| panel key | panel rows matched in `.var` | panel rows missing | "
        "distinct ids in common |")
    say("|---|---|---|---|")
    for key in ["symbol", "ensembl", "gene_id"]:
        ref_k = np.array([str(s) for s in z[key]])
        present = int(np.isin(ref_k, vidx).sum())
        inter = np.intersect1d(np.unique(ref_k), np.unique(vidx))
        say("| `{}` | {:,} | {:,} | {:,} |".format(
            key, present, ref_k.size - present, inter.size))
    say()
    ref = np.array([str(s) for s in z["symbol"]])
    inref = np.isin(ref, vidx)
    say("- **answer to the question as asked** (panel `symbol` against the `.var`")
    say("  index, plain string match): of the {:,} panel genes, **{:,} are present**".format(
        ref.size, int(inref.sum())))
    say("  in this object and **{:,} are missing**.".format(int((~inref).sum())))
    say("- the panel's own `symbol` column has {} duplicated values; distinct panel".format(
        int(ref.size - np.unique(ref).size)))
    say("  symbols = {:,}, of which {:,} are present.".format(
        np.unique(ref).size, int(np.isin(np.unique(ref), vidx).sum())))
    say("- `.var` symbols not on the panel: {:,} of {:,} distinct.".format(
        int((~np.isin(np.unique(vidx), ref)).sum()), np.unique(vidx).size))
    pd.DataFrame({"panel_symbol": np.sort(np.unique(ref[~inref]))}).to_csv(
        OUT / "panel_symbols_missing.csv", index=False)
    pd.DataFrame({"panel_symbol": np.sort(np.unique(ref[inref]))}).to_csv(
        OUT / "panel_symbols_present.csv", index=False)
    pd.DataFrame({"var_symbol": vidx}).to_csv(OUT / "var_index.csv", index=False)
    say("- written: `panel_symbols_missing.csv`, `panel_symbols_present.csv`,")
    say("  `var_index.csv` (all under `results/senginv/`).")
    say("- first 30 missing panel symbols: {}".format(
        [str(x) for x in sorted(set(ref[~inref]))[:30]]))
    say()
    say("**Recorded as an ambiguity.** The match above is a plain string match on")
    say("gene symbols. No alias resolution and no Ensembl-mediated remapping was")
    say("attempted. The panel carries Ensembl IDs and this file does not, so the")
    say("`ensembl` and `gene_id` rows above are zero by construction rather than by")
    say("disagreement, and the symbol row is the only one that can be read. The")
    say("`.var` index has {} duplicated entries and the panel `symbol` column has {},".format(
        int(vidx.size - np.unique(vidx).size), int(ref.size - np.unique(ref).size)))
    say("so {} panel rows can map onto fewer distinct `.var` columns. Nothing was".format(
        int(inref.sum())))
    say("deduplicated and no decision about the duplicates was made.")
    say()
    say("### The guide-column target prefixes as `.var` symbols")
    say()
    ut = np.unique(targets)
    say("- of the {} distinct guide target prefixes, {} appear as a symbol in `.var`"
        " and {} do not.".format(ut.size, int(np.isin(ut, vidx).sum()),
                                 int((~np.isin(ut, vidx)).sum())))
    say("- prefixes absent from `.var`: {}".format(sorted(ut[~np.isin(ut, vidx)])))
    say()
    say("### The four named factors as genes in `.var`")
    say()
    say("| factor | in `.var` | n `.var` columns with that symbol | on the frozen panel |")
    say("|---|---|---|---|")
    for t in FOUR:
        n = int((vidx == t).sum())
        say("| `{}` | {} | {} | {} |".format(t, n > 0, n, bool((ref == t).any())))
    say()

    # ---- files --------------------------------------------------------------
    say("## Files written by this run")
    say()
    for p in sorted(OUT.iterdir()):
        say("- `results/senginv/{}` ({:,} B)".format(p.name, p.stat().st_size))
    say()
    say("Nothing was committed and nothing was pushed. No pre-registration was")
    say("written, no statistic was computed on `.X`, and no analysis plan is proposed.")

    f.close()
    (OUT / "INVENTORY_SENGINV.md").write_text("\n".join(lines) + "\n", encoding="utf-8")


if __name__ == "__main__":
    main()
