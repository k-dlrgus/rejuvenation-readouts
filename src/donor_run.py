"""DONOR — does a transgene push different aged donors' fibroblasts in the
same direction?

The pre-registration is results/donor/PREREG.flag. This script reads it and
never rewrites it. Frozen instruments are read, never refit: the frozen ruler
and GTEx mu/sd (results/fibro/frozen_ruler_ridge_raw.npz) and the
20-component basis results/south7/basis.npz L_mid_full[:, :20].

Data: GEO GSE325735, family SOFT and the human raw counts matrix, both under
data/raw/srsfinv/ (git-ignored; URLs below, fetched by SRSFINV and by this
task respectively).

Stages, run from the repository root:
    python src/donor_run.py stage0     reproduce SRSFINV's stored counts
    python src/donor_run.py build      sample rules, ruler alignment, z, displacements
    python src/donor_run.py read       gate, primary readout, diagnostics
    python src/donor_run.py findings   FINDINGS_DONOR.md
    python src/donor_run.py all        the four in order
"""
from __future__ import annotations

import collections
import gzip
import hashlib
import json
import re
import sys
import time
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import stats as sps

SRC = Path(__file__).resolve().parent
if str(SRC) not in sys.path:
    sys.path.insert(0, str(SRC))

from fibro2_common import load_frozen_ruler  # noqa: E402
from same_run import ci_covers_point      # noqa: E402
from toward_run import percentile_ci, tmm_logcpm_quiet, zscore_frozen  # noqa: E402

ROOT = SRC.parent
OUT = ROOT / "results" / "donor"
PREREG = OUT / "PREREG.flag"
FINDINGS = ROOT / "FINDINGS_DONOR.md"
RAW = ROOT / "data" / "raw" / "srsfinv"
SRSFINV = ROOT / "results" / "srsfinv"
SOUTH7 = ROOT / "results" / "south7"

SOFT = RAW / "GSE325735_family.soft.gz"
COUNTS = RAW / "GSE325735_human_raw_counts.txt.gz"
COUNTS_URL = ("https://ftp.ncbi.nlm.nih.gov/geo/series/GSE325nnn/GSE325735/suppl/"
              "GSE325735_human_raw_counts.txt.gz")
SOFT_MD5 = "6257c37dbec6f49a5cff6f110ca641a4"
COUNTS_BYTES = 26_007_361

SEED = 20260914
N_BOOT = 10_000
N_RANDOM_SPLIT = 100
KDIM = 20
R_FLOOR = 0.3
M_SAME_BAR = 0.85
GATE_P = 0.01
MIN_NQ = 20
MIN_COMP = 19
UPPER_BAR = 0.5
LOWER_BAR = 0.7
SPEC_P_KILL = 0.05
SPEC_P_PASS = 0.01

LINES = {"AG12657": 55, "AG06294A": 65, "AG09697": 67, "AG13369": 68, "AG12587": 69, "AG16359": 79}
LINE_ORDER = list(LINES)
PRIMARY = "Dox_1000ng/mL, 3days"
D5_COND = "Dox_2000ng/mL, 3days"
CONTROL = "mTagBFP2"
NOT_TRANSGENE = {"WT", "NTg", "mTagBFP2", "GFP"}
FLAGS = ("High GC", "<10M reads")

STAGES = ("stage0", "build", "read", "findings")


# --------------------------------------------------------------------------
# bookkeeping
# --------------------------------------------------------------------------
LOG_LINES = []


def log(s=""):
    print(s, flush=True)
    LOG_LINES.append(s)
    with open(OUT / "run.log", "a", encoding="utf-8") as fh:
        fh.write(s + "\n")


def require_prereg():
    if not PREREG.exists():
        raise SystemExit(f"{PREREG} missing. Not computing.")


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


def load(name):
    return json.loads((OUT / name).read_text(encoding="utf-8"))


def md5(p):
    h = hashlib.md5()
    with open(p, "rb") as fh:
        for chunk in iter(lambda: fh.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


class Stop(Exception):
    def __init__(self, key, message, extra=None):
        super().__init__(message)
        self.key, self.message, self.extra = key, message, extra or {}


def stopped():
    p = OUT / "STOP.json"
    return json.loads(p.read_text(encoding="utf-8")) if p.exists() else None


def fmt(x, d=3):
    if x is None:
        return "n/a"
    try:
        x = float(x)
    except (TypeError, ValueError):
        return str(x)
    if not np.isfinite(x):
        return "n/a"
    return f"{x:.{d}f}"


def fmt_p(p):
    if p is None or not np.isfinite(p):
        return "n/a"
    return f"{p:.2e}" if p < 1e-3 else f"{p:.4f}"


def angle(c):
    if c is None or not np.isfinite(c) or abs(c) > 1:
        return np.nan
    return float(np.degrees(np.arccos(c)))


def cos(a, b):
    a = np.asarray(a, np.float64).ravel()
    b = np.asarray(b, np.float64).ravel()
    na, nb = np.linalg.norm(a), np.linalg.norm(b)
    if na < 1e-12 or nb < 1e-12:
        return np.nan
    return float(a @ b / (na * nb))


# --------------------------------------------------------------------------
# SOFT parsing (written fresh for this task; src/srsfinv_run.py not imported)
# --------------------------------------------------------------------------
def parse_soft(path):
    samples, cur = [], None
    with gzip.open(path, "rt", encoding="utf-8") as fh:
        for line in fh:
            line = line.rstrip("\r\n")
            if line.startswith("^"):
                kind = line[1:].split(" = ")[0]
                cur = collections.defaultdict(list) if kind == "SAMPLE" else None
                if cur is not None:
                    samples.append(cur)
                continue
            if cur is not None and line.startswith("!") and " = " in line:
                k, _, v = line[1:].partition(" = ")
                cur[k].append(v)
    rows = []
    for s in samples:
        r = dict(gsm=s["Sample_geo_accession"][0], title=s["Sample_title"][0],
                 organism=s["Sample_organism_ch1"][0], platform=s["Sample_platform_id"][0],
                 molecule=s["Sample_molecule_ch1"][0])
        for c in s["Sample_characteristics_ch1"]:
            k, _, v = c.partition(": ")
            r["ch_" + k.strip().replace(" ", "_")] = v.strip()
        desc = s["Sample_description"]
        m = re.match(r"Library name: (\S+)$", desc[0]) if desc else None
        r["library"] = m.group(1) if m else ""
        r["desc_rest"] = " || ".join(desc[1:])
        rows.append(r)
    df = pd.DataFrame(rows)
    rest = df["desc_rest"]
    df["age"] = pd.to_numeric(rest.str.extract(r"Donor age (\d+)")[0])
    df["passage"] = pd.to_numeric(rest.str.extract(r"Passage (\d+)")[0])
    df["title_M"] = pd.to_numeric(df["title"].str.extract(r"\bM(\d+)\b")[0])
    df["title_P"] = df["title"].str.extract(r", (P\+?\d+)\b")[0]
    df["rep"] = pd.to_numeric(df["title"].str.extract(r"\bRep(\d+)\b")[0])
    df["flags"] = (rest.str.replace(r"^Donor age \d+(, Passage \d+)?", "", regex=True)
                   .str.lstrip(", ").str.strip())
    df["sub_batch"] = df["flags"].where(df["flags"].str.match(r"^Batch \d+$"), "")
    df["title_token"] = df["title"].str.extract(r"Dox-inducible (.+?) transgene")[0]
    return df


def human_fibroblast(df):
    hum = df[df.organism == "Homo sapiens"]
    return hum, hum[hum.get("ch_cell_type") != "iPSC"]


# --------------------------------------------------------------------------
# Stage 0
# --------------------------------------------------------------------------
def donors_table(hum):
    def jn(x):
        return "/".join(str(int(v)) for v in sorted(set(x.dropna())))
    g = hum.groupby("ch_cell_line")
    don = pd.DataFrame({
        "cell type": g["ch_cell_type"].agg(lambda x: "/".join(sorted(set(x)))),
        "n Fibroblast GSMs": g["ch_cell_type"].agg(lambda x: int((x == "Fibroblast").sum())),
        "n iPSC GSMs": g["ch_cell_type"].agg(lambda x: int((x == "iPSC").sum())),
        "described age(s)": g["age"].agg(jn),
        "title M-code(s)": g["title_M"].agg(lambda x: jn(x) or "-"),
        "described passage(s)": g["passage"].agg(lambda x: jn(x) or "-"),
        "title P token(s)": g["title_P"].agg(lambda x: "/".join(sorted(set(x.dropna()))) or "-"),
        "sex": "not stated",
        "n GSMs": g.size(),
        "n batches": g["ch_batch"].nunique(),
    })
    don.index.name = "cell line"
    return don.astype(str)


def counts_header():
    with gzip.open(COUNTS, "rt", encoding="utf-8", newline="") as fh:
        head = fh.readline().rstrip("\r\n").split("\t")
    return head


def fetch_counts():
    """The matrix is git-ignored. Fetched once from GEO if absent (curl in the original run)."""
    if COUNTS.exists():
        return False
    import urllib.request
    RAW.mkdir(parents=True, exist_ok=True)
    urllib.request.urlretrieve(COUNTS_URL, COUNTS)
    return True


def stage0():
    if fetch_counts():
        log(f"[stage0] fetched {COUNTS_URL}")
    rows = []

    def chk(q, stored, new):
        rows.append(dict(quantity=q, stored=str(stored), new=str(new), match=str(stored) == str(new)))

    chk("SOFT md5", SOFT_MD5, md5(SOFT))
    df = parse_soft(SOFT)
    hum, hf = human_fibroblast(df)
    meta = json.loads((SRSFINV / "inventory_meta.json").read_text(encoding="utf-8"))
    chk("n GSMs", 839, len(df))
    chk("n GSMs (inventory_meta n_gsm)", meta["n_gsm"], len(df))
    chk("n human GSMs", 807, len(hum))
    chk("n human GSMs (inventory_meta n_human)", meta["n_human"], len(hum))
    chk("n human fibroblast GSMs (inventory_meta n_fibroblast)", meta["n_fibroblast"], len(hf))
    tg = hf[~hf.ch_genotype.isin(NOT_TRANSGENE)]
    chk("n transgene GSMs (human fibroblast, genotype not WT/NTg/mTagBFP2/GFP)", 664, len(tg))
    chk("n transgene genotype values", 96, tg.ch_genotype.nunique())
    n_age = tg.groupby("ch_genotype")["age"].nunique()
    multi = sorted(n_age[n_age > 1].index)
    chk("n transgenes at more than one described age", 89, len(multi))
    stored_multi = sorted(set(meta["genotypes_multi_age_fibroblast"]) - NOT_TRANSGENE)
    chk("transgenes at more than one age, names (inventory_meta)", ";".join(stored_multi), ";".join(multi))

    stored = pd.read_csv(SRSFINV / "donors.csv", dtype=str).set_index("cell line")
    new = donors_table(hum)
    chk("donors.csv n rows (cell lines)", len(stored), len(new))
    chk("donors.csv columns", ";".join(stored.columns), ";".join(new.columns))
    for line in sorted(set(stored.index) | set(new.index)):
        for col in stored.columns:
            s = stored.loc[line, col] if line in stored.index else "<absent>"
            n = new.loc[line, col] if (line in new.index and col in new.columns) else "<absent>"
            chk(f"donors.csv {line} / {col}", s, n)

    chk("counts file bytes", COUNTS_BYTES, COUNTS.stat().st_size)
    head = counts_header()
    chk("counts header first field", "Geneid", head[0])
    chk("counts header sample columns (inventory_meta counts_header_cols)", meta["counts_header_cols"], len(head) - 1)
    libs = set(hum["library"])
    chk("counts header columns == human library names (set)", True, set(head[1:]) == libs)
    chk("counts header duplicate columns", 0, len(head) - 1 - len(set(head[1:])))

    z = np.load(SOUTH7 / "basis.npz", allow_pickle=True)
    L = np.asarray(z["L_mid_full"], np.float64)[:, :KDIM]
    err = float(np.max(np.abs(L.T @ L - np.eye(KDIM))))
    chk("basis max|L'L - I| <= 1e-8", True, err <= 1e-8)
    chk("basis rows (ruler genes)", 23485, L.shape[0])
    chk("basis n_mid", 308, int(np.asarray(z["n_mid"]).ravel()[0]))

    tab = pd.DataFrame(rows)
    tab.to_csv(OUT / "stage0_reproduction.csv", index=False)
    bad = tab[~tab["match"]]
    n_donor = int(tab.quantity.str.startswith("donors.csv").sum())
    log(f"[stage0] compared {len(tab)} quantities ({n_donor} donors.csv cells); mismatches {len(bad)}")
    log(f"[stage0] counts md5 {md5(COUNTS)} (recorded; no stored md5); basis err {err:.3e}")
    dump("stage0.json", dict(n_compared=len(tab), n_mismatch=len(bad), counts_md5=md5(COUNTS),
                             counts_bytes=COUNTS.stat().st_size, basis_err=err,
                             mismatches=bad.to_dict("records")))
    if len(bad):
        for r in bad.itertuples():
            log(f"  MISMATCH {r.quantity}: stored={r.stored!r} new={r.new!r}")
        raise Stop("prior_results_not_reproduced",
                   f"{len(bad)} of {len(tab)} Stage 0 quantities did not reproduce.",
                   dict(mismatches=bad.to_dict("records")))


# --------------------------------------------------------------------------
# build: sample rules, cells, ruler alignment, z, displacements
# --------------------------------------------------------------------------
def select_samples(df):
    """Apply flag section 2 to the six lines' GSMs. Returns (kept, excluded)."""
    _, hf = human_fibroblast(df)
    six = hf[hf.ch_cell_line.isin(LINES)].copy()
    reasons = []
    for r in six.itertuples():
        why = []
        if r.ch_cell_type != "Fibroblast":
            why.append("A8 cell type")
        if r.molecule != "polyA RNA":
            why.append("A12 not polyA")
        if r.ch_genotype in ("WT", "NTg") or (isinstance(r.title_token, str) and r.title_token == "NTg"):
            why.append("A3 NTg/WT")
        elif not (isinstance(r.title_token, str) and r.title_token == r.ch_genotype):
            why.append(f"A7 title token {r.title_token!r} != genotype {r.ch_genotype!r}")
        if isinstance(r.title_P, str) and np.isfinite(r.passage) and r.title_P != f"P{int(r.passage)}":
            why.append(f"A11 title {r.title_P} != described passage {int(r.passage)}")
        if r.ch_treatment not in (PRIMARY, D5_COND):
            why.append(f"A5 treatment {r.ch_treatment!r} outside the task")
        reasons.append("; ".join(why))
    six["exclusion"] = reasons
    return six[six.exclusion == ""].copy(), six[six.exclusion != ""].copy()


def build_cells(kept):
    """Flag R2-R4: one batch cell per (t, L, condition), mTagBFP2 controls, eligibility."""
    kept = kept.assign(bcell=kept.ch_batch + "|" + kept.sub_batch)
    tg = kept[kept.ch_genotype != CONTROL]
    ctl = kept[kept.ch_genotype == CONTROL]
    order = lambda d: d.sort_values(["rep", "gsm"], na_position="last").gsm.tolist()  # noqa: E731
    cells = []
    for (t, line, cond), g in tg.groupby(["ch_genotype", "ch_cell_line", "ch_treatment"]):
        vc = g.bcell.value_counts()
        rec = dict(transgene=t, line=line, age=LINES[line], cond=cond,
                   bcells_available="; ".join(f"{k} ({v})" for k, v in vc.items()))
        if len(vc) > 1 and vc.iloc[0] == vc.iloc[1]:
            rec.update(bcell="", t_gsms=[], c_gsms=[], eligible=False, reason="R2 tie between batch cells")
            cells.append(rec)
            continue
        bc = vc.index[0]
        tt = g[g.bcell == bc]
        cc = ctl[(ctl.ch_cell_line == line) & (ctl.ch_treatment == cond) & (ctl.bcell == bc)]
        rec.update(bcell=bc, batch=bc.split("|")[0], t_gsms=order(tt), c_gsms=order(cc))
        why = []
        if len(tt) < 2:
            why.append(f"R3 {len(tt)} transgene GSM(s)")
        if len(cc) < 2:
            why.append(f"R3 {len(cc)} control GSM(s)")
        rec.update(eligible=not why, reason="; ".join(why))
        cells.append(rec)
    return pd.DataFrame(cells)


def align_counts(gsm_lib, panel_gsms, frozen):
    """Read the counts matrix, align rows to the ruler on unversioned Ensembl IDs."""
    t0 = time.time()
    raw = pd.read_csv(COUNTS, sep="\t", index_col=0, compression="gzip")
    raw.columns = [c.strip() for c in raw.columns]
    log(f"[counts] read {raw.shape[0]} rows x {raw.shape[1]} columns in {time.time() - t0:.1f}s; "
        f"dtypes {sorted(set(map(str, raw.dtypes)))}")
    ids = np.asarray(raw.index).astype(str)
    stripped = np.array([i.split(".")[0] for i in ids])
    ruler_ens = np.asarray(frozen["ensembl"]).astype(str)
    pos = {}
    dup_rows = []
    for i, e in enumerate(stripped):
        if e in pos:
            dup_rows.append(i)
        else:
            pos[e] = i
    src = np.array([pos.get(e, -1) for e in ruler_ens])
    mapped = src >= 0
    ruler_set = set(ruler_ens)
    dup_on_ruler = [i for i in dup_rows if stripped[i] in ruler_set]
    cols = [gsm_lib[g] for g in panel_gsms]
    M = raw[cols].to_numpy(np.float64)
    dup_nonzero = int(sum(M[i].sum() > 0 for i in dup_on_ruler))
    C = np.zeros((len(panel_gsms), ruler_ens.size), np.float64)
    C[:, mapped] = M[src[mapped]].T
    lib_all = M.sum(0)
    info = dict(file_rows=int(ids.size), file_rows_unique_unversioned=int(len(pos)),
                file_rows_duplicate_unversioned=len(dup_rows),
                duplicates_on_ruler=len(dup_on_ruler), duplicates_on_ruler_nonzero_in_panel=dup_nonzero,
                ruler_genes=int(ruler_ens.size), ruler_mapped=int(mapped.sum()),
                ruler_unmapped=int((~mapped).sum()),
                file_rows_not_on_ruler=int(sum(1 for e in pos if e not in ruler_set)),
                versioned_exact_matches=int(np.isin(ids, np.asarray(frozen["gene_id"]).astype(str)).sum()),
                panel_library_size_all_rows_median=float(np.median(lib_all)),
                panel_library_size_ruler_median=float(np.median(C.sum(1))),
                panel_share_of_counts_on_ruler_median=float(np.median(C.sum(1) / lib_all)))
    return C, mapped, info


def build():
    df = parse_soft(SOFT)
    kept, excl = select_samples(df)
    excl[["gsm", "title", "ch_cell_line", "ch_genotype", "ch_treatment", "ch_batch", "exclusion"]] \
        .to_csv(OUT / "excluded_gsms.csv", index=False)
    cells = build_cells(kept)
    if cells.empty:
        raise Stop("ambiguity_unresolved", "no transgene cell survives the metadata rules")
    cells["n_t"] = cells.t_gsms.map(len)
    cells["n_c"] = cells.c_gsms.map(len)
    elig = cells[cells.eligible]
    panel = sorted(set(g for lst in elig.t_gsms for g in lst) | set(g for lst in elig.c_gsms for g in lst))
    meta = df.set_index("gsm")
    gsm_lib = meta["library"].to_dict()
    log(f"[rules] six-line GSMs {len(kept) + len(excl)}; kept {len(kept)}; excluded {len(excl)}")
    for k, v in collections.Counter(r.split(" ")[0] for e in excl.exclusion for r in e.split("; ")).items():
        log(f"  exclusion rule {k}: {v} GSM-reasons")
    log(f"[cells] {len(cells)} (t, L, condition) cells; eligible {int(cells.eligible.sum())}; panel GSMs {len(panel)}")

    frozen = load_frozen_ruler()
    C, mapped, info = align_counts(gsm_lib, panel, frozen)
    log(f"[align] file rows {info['file_rows']}; ruler genes mapped {info['ruler_mapped']} of "
        f"{info['ruler_genes']}; unmapped {info['ruler_unmapped']} (z = 0); duplicate unversioned IDs "
        f"{info['file_rows_duplicate_unversioned']} (on ruler {info['duplicates_on_ruler']}, nonzero "
        f"in panel {info['duplicates_on_ruler_nonzero_in_panel']})")
    logcpm, nf = tmm_logcpm_quiet(C)
    Z, missing = zscore_frozen(logcpm, frozen, C)
    info.update(panel_gsms=len(panel), z_zero_columns=int(missing.sum()),
                z_zero_columns_mapped_but_all_zero=int((missing & mapped).sum()),
                tmm_factor_min=float(nf.min()), tmm_factor_max=float(nf.max()))
    log(f"[z] all-zero ruler columns across the panel {int(missing.sum())} "
        f"({int((missing & mapped).sum())} of them mapped but never counted)")

    zb = np.load(SOUTH7 / "basis.npz", allow_pickle=True)
    L20 = np.asarray(zb["L_mid_full"], np.float64)[:, :KDIM]
    row = {g: i for i, g in enumerate(panel)}
    lib = dict(zip(panel, C.sum(1)))
    pan = meta.loc[panel, ["title", "ch_cell_line", "ch_genotype", "ch_treatment", "ch_batch",
                          "sub_batch", "platform", "molecule", "passage", "rep", "flags"]].copy()
    pan["library_size_ruler"] = [lib[g] for g in panel]
    pan["tmm_factor"] = nf
    pan.to_csv(OUT / "panel_gsms.csv")

    # displacements: full z-space and projected, plus split halves
    full, proj, halves, split_mode = {}, {}, {}, {}
    for r in elig.itertuples():
        key = f"{r.transgene}|{r.line}|{r.cond}"
        zt = Z[[row[g] for g in r.t_gsms]]
        zc = Z[[row[g] for g in r.c_gsms]]
        d = zt.mean(0) - zc.mean(0)
        full[key], proj[key] = d, d @ L20
        hs = []
        if len(r.t_gsms) == 2 and len(r.c_gsms) == 2:
            split_mode[key] = "enumerated (2 control pairings)"
            for x, y in ((0, 1), (1, 0)):
                dA = zt[0] - zc[x]
                dB = zt[1] - zc[y]
                hs.append((dA @ L20, dB @ L20))
        else:
            split_mode[key] = f"{N_RANDOM_SPLIT} random disjoint splits"
            li = LINE_ORDER.index(r.line)
            ti = sorted(cells.transgene.unique()).index(r.transgene)
            for s in range(N_RANDOM_SPLIT):
                rng = np.random.default_rng([SEED, li, ti, s])
                pt, pc = rng.permutation(len(zt)), rng.permutation(len(zc))
                ht, hc = (len(zt) + 1) // 2, (len(zc) + 1) // 2
                dA = zt[pt[:ht]].mean(0) - zc[pc[:hc]].mean(0)
                dB = zt[pt[ht:]].mean(0) - zc[pc[hc:]].mean(0)
                hs.append((dA @ L20, dB @ L20))
        halves[key] = np.array(hs)
    keys = sorted(proj)
    np.savez_compressed(OUT / "displacements.npz", keys=np.array(keys),
                        proj=np.array([proj[k] for k in keys]),
                        halves=np.array([halves[k] for k in keys], dtype=object),
                        norm_full=np.array([np.linalg.norm(full[k]) for k in keys]))
    np.save(OUT / "displacements_full.npy", np.array([full[k] for k in keys]).astype(np.float32))

    out = cells.copy()
    out["t_gsms"] = out.t_gsms.map(";".join)
    out["c_gsms"] = out.c_gsms.map(";".join)
    out["split_mode"] = [split_mode.get(f"{r.transgene}|{r.line}|{r.cond}", "") for r in cells.itertuples()]
    out["passages_t"] = [";".join(str(int(meta.loc[g, "passage"])) for g in lst) for lst in cells.t_gsms]
    out["passages_c"] = [";".join(str(int(meta.loc[g, "passage"])) for g in lst) for lst in cells.c_gsms]
    out["flagged_gsms"] = [";".join(g for g in a + b if meta.loc[g, "flags"] in FLAGS)
                           for a, b in zip(cells.t_gsms, cells.c_gsms)]
    out["platforms"] = [";".join(sorted(set(meta.loc[a + b, "platform"]))) if (a + b) else ""
                        for a, b in zip(cells.t_gsms, cells.c_gsms)]
    out["lib_t_mean"] = [float(np.mean([lib[g] for g in a])) if (a and all(g in lib for g in a)) else np.nan
                         for a in cells.t_gsms]
    out.to_csv(OUT / "cells.csv", index=False)
    info["split_modes"] = dict(collections.Counter(split_mode.values()))
    info["eligible_primary_batches"] = dict(collections.Counter(elig[elig.cond == PRIMARY].batch))
    info["eligible_primary_platforms"] = dict(collections.Counter(
        out[out.eligible & (out.cond == PRIMARY)].platforms))
    dump("build.json", info)
    log(f"[build] split modes {info['split_modes']}; primary eligible batches "
        f"{info['eligible_primary_batches']}; platforms {info['eligible_primary_platforms']}")


# --------------------------------------------------------------------------
# readout
# --------------------------------------------------------------------------
def load_vectors():
    z = np.load(OUT / "displacements.npz", allow_pickle=True)
    keys = [str(k) for k in z["keys"]]
    proj = {k: v for k, v in zip(keys, z["proj"])}
    halves = {k: np.asarray(h, np.float64) for k, h in zip(keys, z["halves"])}
    return proj, halves


def sb(c):
    return 2 * c / (1 + c) if (np.isfinite(c) and c > 0) else np.nan


def binom_p(k, n):
    return float(sps.binomtest(int(k), int(n), 0.05, alternative="greater").pvalue) if n > 0 else np.nan


def p95(v):
    return float(np.percentile(np.asarray(v, float), 95))


def boot_median(vals):
    v = np.asarray(vals, float)
    rng = np.random.default_rng(SEED)
    idx = rng.integers(0, v.size, size=(N_BOOT, v.size))
    meds = np.median(v[idx], axis=1)
    lo, hi = percentile_ci(meds)
    return lo, hi


def by_line(cond, proj):
    """{transgene: {line: key}} for one condition."""
    out = collections.defaultdict(dict)
    for k in proj:
        t, line, c = k.split("|")
        if c == cond:
            out[t][line] = k
    return out


def readout(cond, proj, halves, specificity=True, ci=True):
    """Flag section 5 on one condition. Vectors are 20-component."""
    tl = by_line(cond, proj)
    cell_rows = {}
    for t, lines in tl.items():
        for line, k in lines.items():
            h = halves[k]
            ch = float(np.nanmean([cos(a, b) for a, b in h]))
            r = sb(ch)
            cell_rows[k] = dict(transgene=t, line=line, age=LINES[line], c_half=ch, floor_angle=angle(ch),
                                r=r, qualifying=bool(np.isfinite(r) and r >= R_FLOOR))
    cand = sorted(t for t, lines in tl.items() if len(lines) >= 2)
    qcell = {k for k, v in cell_rows.items() if v["qualifying"]}
    pair_rows, tg_rows = [], []
    for t in cand:
        ls = [l for l in LINE_ORDER if l in tl[t]]
        pairs = []
        for i in range(len(ls)):
            for j in range(i + 1, len(ls)):
                kA, kB = tl[t][ls[i]], tl[t][ls[j]]
                co = cos(proj[kA], proj[kB])
                q = kA in qcell and kB in qcell
                cc = co / np.sqrt(cell_rows[kA]["r"] * cell_rows[kB]["r"]) if q else np.nan
                pr = dict(transgene=t, line_A=ls[i], line_B=ls[j], age_A=LINES[ls[i]], age_B=LINES[ls[j]],
                          c_obs=co, angle=angle(co), r_A=cell_rows[kA]["r"], r_B=cell_rows[kB]["r"],
                          qualifying_pair=q, c_corr=cc)
                pair_rows.append(pr)
                if q:
                    pairs.append(pr)
        rec = dict(transgene=t, n_lines=len(ls), lines="/".join(str(LINES[l]) for l in ls),
                   n_qual_pairs=len(pairs), qualifying=bool(pairs),
                   c_corr=float(np.mean([p["c_corr"] for p in pairs])) if pairs else np.nan,
                   c_obs_mean=float(np.mean([p["c_obs"] for p in pairs])) if pairs else np.nan)
        rec["n_qual_lines"] = len({p["line_A"] for p in pairs} | {p["line_B"] for p in pairs})
        tg_rows.append(rec)
    tgdf = pd.DataFrame(tg_rows)
    res = dict(cond=cond, n_candidates=len(cand), n_cells=len(cell_rows), n_cells_qualifying=len(qcell))
    if tgdf.empty:
        res.update(n_q=0)
        return res, pd.DataFrame(cell_rows.values()), pd.DataFrame(pair_rows), tgdf
    q = tgdf[tgdf.qualifying]
    fin = q[np.isfinite(q.c_corr)]
    res.update(n_q=int(len(q)), n_nonfinite=int(len(q) - len(fin)))
    M = float(np.median(fin.c_corr)) if len(fin) else np.nan
    res.update(M=M, angle_M=angle(M))
    if ci and len(fin) >= 2:
        lo, hi = boot_median(fin.c_corr.to_numpy())
        res.update(ci_lo=lo, ci_hi=hi, ci_valid=ci_covers_point(M, lo, hi))
    if specificity:
        spec = {}
        for t in q.transgene:
            dirs = []
            for p in pair_rows:
                if p["transgene"] == t and p["qualifying_pair"]:
                    dirs += [(p["line_A"], p["line_B"]), (p["line_B"], p["line_A"])]
            targets = {b for _, b in dirs}
            G = [g for g in tl if g != t and all(b in tl[g] and tl[g][b] in qcell for b in targets)]
            s = float(np.mean([cos(proj[tl[t][a]], proj[tl[t][b]]) for a, b in dirs]))
            m = [float(np.mean([cos(proj[tl[t][a]], proj[tl[g][b]]) for a, b in dirs])) for g in G]
            testable = len(G) >= MIN_COMP
            thr = p95(m) if m else np.nan
            spec[t] = dict(s=s, n_comp=len(G), p95=thr, testable=testable,
                           specific=bool(testable and s > thr))
        tgdf["s"] = tgdf.transgene.map(lambda t: spec.get(t, {}).get("s", np.nan))
        tgdf["n_comp"] = tgdf.transgene.map(lambda t: spec.get(t, {}).get("n_comp", np.nan))
        tgdf["p95_comp"] = tgdf.transgene.map(lambda t: spec.get(t, {}).get("p95", np.nan))
        tgdf["spec_testable"] = tgdf.transgene.map(lambda t: spec.get(t, {}).get("testable", False))
        tgdf["specific"] = tgdf.transgene.map(lambda t: spec.get(t, {}).get("specific", False))
        n_test = int(sum(v["testable"] for v in spec.values()))
        k = int(sum(v["specific"] for v in spec.values()))
        res.update(spec_k=k, spec_n=n_test, spec_baseline=0.05 * n_test, spec_p=binom_p(k, n_test),
                   spec_untestable=[t for t, v in spec.items() if not v["testable"]])
    return res, pd.DataFrame(cell_rows.values()), pd.DataFrame(pair_rows), tgdf


def gate(proj, halves):
    """Flag section 6: same-context control, uncorrected."""
    tl = by_line(PRIMARY, proj)
    cand = sorted(t for t, lines in tl.items() if len(lines) >= 2)
    rows = []
    for t in cand:
        ls = [l for l in LINE_ORDER if l in tl[t]]
        ch = [float(np.mean([cos(a, b) for a, b in halves[tl[t][l]]])) for l in ls]
        G = [g for g in cand if g != t and all(l in tl[g] for l in ls)]
        m = []
        for g in G:
            v = []
            for l in ls:
                ht, hg = halves[tl[t][l]], halves[tl[g][l]]
                for s in range(min(len(ht), len(hg))):
                    v.append(cos(ht[s][0], hg[s][1]))
                    v.append(cos(ht[s][1], hg[s][0]))
            m.append(float(np.mean(v)))
        c_same = float(np.mean(ch))
        thr = p95(m) if m else np.nan
        rows.append(dict(transgene=t, lines="/".join(str(LINES[l]) for l in ls), n_lines=len(ls),
                         c_same=c_same, angle_same=angle(c_same), n_comp=len(G), p95_comp=thr,
                         testable=len(G) >= MIN_COMP, specific=bool(len(G) >= MIN_COMP and c_same > thr)))
    gd = pd.DataFrame(rows)
    res = dict(n_candidates=len(cand))
    if gd.empty:
        res.update(passed=False, reason="empty candidate set")
        return res, gd
    M_same = float(np.median(gd.c_same))
    n_test, k = int(gd.testable.sum()), int(gd.specific.sum())
    p = binom_p(k, n_test)
    passed = bool(M_same >= M_SAME_BAR and np.isfinite(p) and p < GATE_P)
    q = gd.c_same.quantile([0.25, 0.75])
    res.update(M_same=M_same, angle_M_same=angle(M_same), q25=float(q.iloc[0]), q75=float(q.iloc[1]),
               spec_k=k, spec_n=n_test, spec_baseline=0.05 * n_test, spec_p=p, passed=passed,
               fail_reasons=[s for s, bad in (
                   (f"M_same {M_same:.4f} < {M_SAME_BAR}", M_same < M_SAME_BAR),
                   (f"gate specificity p {p:.3g} not < {GATE_P}", not (np.isfinite(p) and p < GATE_P)))
                   if bad])
    return res, gd


def reading(g, pr):
    if not g["passed"]:
        return "UNMEASURABLE", "gate"
    if pr["n_q"] < MIN_NQ:
        return "UNMEASURABLE", "n_q"
    valid = pr.get("ci_valid", False)
    if (valid and pr["ci_hi"] < UPPER_BAR) or not (pr["spec_p"] < SPEC_P_KILL):
        return "DONOR-SPECIFIC", ""
    if valid and pr["ci_lo"] >= LOWER_BAR and pr["spec_p"] < SPEC_P_PASS:
        return "CONSISTENT", ""
    return "INTERMEDIATE", ""


def remove_line_mean(cond, proj, halves):
    """D1: subtract each line's mean displacement (and each half's mean, per split)."""
    tl = by_line(cond, proj)
    p2, h2 = dict(proj), dict(halves)
    for line in LINE_ORDER:
        ks = [lines[line] for lines in tl.values() if line in lines]
        if not ks:
            continue
        mu = np.mean([proj[k] for k in ks], axis=0)
        ns = min(len(halves[k]) for k in ks)
        hmu = np.mean([halves[k][:ns] for k in ks], axis=0)
        for k in ks:
            p2[k] = proj[k] - mu
            h2[k] = halves[k][:ns] - hmu
    return p2, h2


def read():
    proj, halves = load_vectors()
    g, gd = gate(proj, halves)
    gd.to_csv(OUT / "gate_table.csv", index=False)
    dump("gate.json", g)
    log(f"[gate] candidates {g['n_candidates']}; M_same {fmt(g.get('M_same'), 4)} (bar {M_SAME_BAR}); "
        f"specificity {g.get('spec_k')} of {g.get('spec_n')} vs baseline {fmt(g.get('spec_baseline'), 2)}, "
        f"p {fmt_p(g.get('spec_p'))}; passed {g['passed']}")
    if not g["passed"]:
        raise Stop("unmeasurable", "the same-context gate failed: " + "; ".join(g.get("fail_reasons", [])), g)

    pr, cdf, pdf, tdf = readout(PRIMARY, proj, halves)
    cdf.to_csv(OUT / "primary_cells.csv", index=False)
    pdf.to_csv(OUT / "primary_pairs.csv", index=False)
    tdf.to_csv(OUT / "primary_transgenes.csv", index=False)
    rd, why = reading(g, pr)
    pr.update(reading=rd, reading_detail=why)
    dump("primary.json", pr)
    log(f"[primary] n_q {pr['n_q']}; M {fmt(pr.get('M'), 4)} CI {fmt(pr.get('ci_lo'), 4)} to "
        f"{fmt(pr.get('ci_hi'), 4)} valid {pr.get('ci_valid')}; angle {fmt(pr.get('angle_M'), 1)}; "
        f"specificity {pr.get('spec_k')} of {pr.get('spec_n')} vs {fmt(pr.get('spec_baseline'), 2)}, "
        f"p {fmt_p(pr.get('spec_p'))}; reading {rd}")
    diagnostics(proj, halves, pdf, tdf)


def diagnostics(proj, halves, pdf, tdf):
    diag = {}
    # D1
    p1, h1 = remove_line_mean(PRIMARY, proj, halves)
    r1, c1, pp1, t1 = readout(PRIMARY, p1, h1)
    t1.to_csv(OUT / "d1_line_mean_removed_transgenes.csv", index=False)
    pp1.to_csv(OUT / "d1_line_mean_removed_pairs.csv", index=False)
    diag["D1"] = r1
    # D5 and D2
    r5, c5, pp5, t5 = readout(D5_COND, proj, halves, specificity=False, ci=False)
    c5.to_csv(OUT / "d5_dox2000_cells.csv", index=False)
    pp5.to_csv(OUT / "d5_dox2000_pairs.csv", index=False)
    t5.to_csv(OUT / "d5_dox2000_transgenes.csv", index=False)
    diag["D5"] = r5
    cells = pd.read_csv(OUT / "cells.csv")
    cells = cells[cells.eligible]
    batch = {f"{r.transgene}|{r.line}|{r.cond}": r.batch for r in cells.itertuples()}
    rr = {}
    for cond in (PRIMARY, D5_COND):
        _, cc, _, _ = readout(cond, proj, halves, specificity=False, ci=False)
        for r in cc.itertuples():
            rr[f"{r.transgene}|{r.line}|{cond}"] = r.r
    both = sorted(set(by_line(PRIMARY, proj)) & set(by_line(D5_COND, proj)))
    rows = []
    for t in both:
        ks = sorted(k for k in proj if k.split("|")[0] == t)
        for i in range(len(ks)):
            for j in range(i + 1, len(ks)):
                a, b = ks[i], ks[j]
                la, lb = a.split("|")[1], b.split("|")[1]
                same_line, same_batch = la == lb, batch[a] == batch[b]
                cat = ("same line / different batch" if same_line else
                       "different line / same batch" if same_batch else "different line / different batch")
                co = cos(proj[a], proj[b])
                ra, rb = rr[a], rr[b]
                q = np.isfinite(ra) and np.isfinite(rb) and ra >= R_FLOOR and rb >= R_FLOOR
                rows.append(dict(transgene=t, cell_A=a, cell_B=b, batch_A=batch[a], batch_B=batch[b],
                                 category=cat, c_obs=co, angle=angle(co),
                                 c_corr=co / np.sqrt(ra * rb) if q else np.nan))
    d2 = pd.DataFrame(rows)
    d2.to_csv(OUT / "d2_batch_vs_donor_pairs.csv", index=False)
    if not d2.empty:
        summ = d2.groupby("category").agg(n_pairs=("c_obs", "size"), median_c_obs=("c_obs", "median"),
                                          n_corr=("c_corr", lambda x: int(np.isfinite(x).sum())),
                                          median_c_corr=("c_corr", "median")).reset_index()
        summ["median_angle"] = summ.median_c_obs.map(angle)
        summ.to_csv(OUT / "d2_batch_vs_donor_summary.csv", index=False)
        diag["D2"] = dict(transgenes=both, summary=summ.to_dict("records"))
    # D3
    q = tdf[tdf.qualifying & np.isfinite(tdf.c_corr)].copy()
    prim = cells[cells.cond == PRIMARY]
    nrep = sorted(set(prim.n_t)) if len(prim) else []
    libmean = prim.groupby("transgene").lib_t_mean.mean()
    q["lib_mean"] = q.transgene.map(libmean)
    d3 = dict(replicate_counts_in_primary_cells=nrep, control_counts=sorted(set(prim.n_c)))
    if len(q) >= 3:
        s1 = sps.spearmanr(q.c_corr, q.lib_mean)
        s2 = sps.spearmanr(q.c_corr, q.n_qual_lines)
        d3.update(rho_lib=float(s1.statistic), p_lib=float(s1.pvalue),
                  rho_nlines=float(s2.statistic), p_nlines=float(s2.pvalue), n=int(len(q)))
    diag["D3"] = d3
    # D4 full gene space
    zf = np.load(OUT / "displacements_full.npy").astype(np.float64)
    keys = [str(k) for k in np.load(OUT / "displacements.npz", allow_pickle=True)["keys"]]
    F = dict(zip(keys, zf))
    qp = pdf[pdf.qualifying_pair].copy()
    qp["c_obs_full"] = [cos(F[f"{r.transgene}|{r.line_A}|{PRIMARY}"], F[f"{r.transgene}|{r.line_B}|{PRIMARY}"])
                        for r in qp.itertuples()]
    qp["angle_full"] = qp.c_obs_full.map(angle)
    qp.to_csv(OUT / "d4_full_gene_space_pairs.csv", index=False)
    diag["D4"] = dict(n_pairs=int(len(qp)), median_c_obs_full=float(qp.c_obs_full.median()) if len(qp) else None,
                      median_angle_full=angle(float(qp.c_obs_full.median())) if len(qp) else None,
                      median_c_obs_20=float(qp.c_obs.median()) if len(qp) else None,
                      median_angle_20=angle(float(qp.c_obs.median())) if len(qp) else None)
    # D6 norm bookkeeping
    nf = np.load(OUT / "displacements.npz", allow_pickle=True)["norm_full"]
    share = np.array([np.linalg.norm(proj[k]) ** 2 for k in keys]) / np.asarray(nf) ** 2
    pd.DataFrame(dict(cell=keys, norm_full=nf, share_retained_20=share)).to_csv(OUT / "d6_norms.csv", index=False)
    diag["D6"] = dict(median_share=float(np.median(share)), q05=float(np.quantile(share, 0.05)),
                      q95=float(np.quantile(share, 0.95)))
    dump("diagnostics.json", diag)
    log(f"[D1 line mean removed] n_q {r1['n_q']}; M {fmt(r1.get('M'), 4)} CI {fmt(r1.get('ci_lo'), 4)} to "
        f"{fmt(r1.get('ci_hi'), 4)}; spec {r1.get('spec_k')} of {r1.get('spec_n')} p {fmt_p(r1.get('spec_p'))}")
    log(f"[D5 dox2000] n_q {r5['n_q']}; M {fmt(r5.get('M'), 4)}")
    log(f"[D4 full gene space] median angle {fmt(diag['D4']['median_angle_full'], 1)} vs 20-comp "
        f"{fmt(diag['D4']['median_angle_20'], 1)}")


# --------------------------------------------------------------------------
# findings
# --------------------------------------------------------------------------
def md(df, cols=None, d=3):
    df = df if cols is None else df[cols]
    out = ["| " + " | ".join(map(str, df.columns)) + " |", "|" + "---|" * len(df.columns)]
    for row in df.itertuples(index=False):
        cells = []
        for v in row:
            if isinstance(v, (bool, np.bool_)):
                cells.append(str(bool(v)))
            elif isinstance(v, (float, np.floating)):
                cells.append(fmt(v, d) if np.isfinite(v) else "n/a")
            else:
                cells.append("" if (v is None or (isinstance(v, float) and np.isnan(v))) else str(v))
        out.append("| " + " | ".join(cells) + " |")
    return "\n".join(out)


def duplicate_ids():
    """The duplicated unversioned IDs that sit on the ruler, with their panel count sums."""
    raw = pd.read_csv(COUNTS, sep="\t", index_col=0, compression="gzip")
    raw.columns = [c.strip() for c in raw.columns]
    panel = pd.read_csv(OUT / "panel_gsms.csv", index_col=0)
    lib = parse_soft(SOFT).set_index("gsm")["library"]
    sub = raw[[lib[g] for g in panel.index]]
    ids = pd.Series(raw.index.astype(str))
    st = ids.str.split(".").str[0]
    ruler = set(np.asarray(load_frozen_ruler()["ensembl"]).astype(str))
    dup = st.duplicated(keep="first") & st.isin(ruler)
    rows = []
    for i in np.flatnonzero(dup.to_numpy()):
        first = int(np.flatnonzero((st == st[i]).to_numpy())[0])
        rows.append(dict(unused_row=ids[i], used_row=ids[first],
                         unused_panel_sum=int(sub.iloc[i].sum()), used_panel_sum=int(sub.iloc[first].sum())))
    return pd.DataFrame(rows)


def reading_text(rd, why, g, pr):
    base = ("This compares aged donors with each other. No reading licenses a claim about transfer "
            "between an aged and a young context, nor about the Southard screen.")
    stats = (f"median reliability-corrected cosine {fmt(pr.get('M'))} (95% CI {fmt(pr.get('ci_lo'))} to "
             f"{fmt(pr.get('ci_hi'))}, {fmt(pr.get('angle_M'), 1)} degrees)")
    spec = (f"{pr.get('spec_k')} of {pr.get('spec_n')} transgenes were specific against a chance baseline of "
            f"{fmt(pr.get('spec_baseline'), 2)} (binomial p = {fmt_p(pr.get('spec_p'))}, approximate)")
    if rd == "UNMEASURABLE":
        if why == "gate":
            inner = (f"The same-context control did not recover agreement at single-replicate depth "
                     f"(M_same = {fmt(g.get('M_same'))}, bar 0.85; gate specificity {g.get('spec_k')} of "
                     f"{g.get('spec_n')} against a chance baseline of {fmt(g.get('spec_baseline'), 2)}, "
                     f"binomial p = {fmt_p(g.get('spec_p'))}).")
        else:
            inner = (f"Only {pr['n_q']} transgenes had sufficient within-line reliability, fewer than the 20 "
                     "fixed in advance.")
        return (f"UNMEASURABLE. {inner} No statement is made about whether a transgene's direction carries "
                "between aged donor lines. This is a statement about the measurement, not about the biology. "
                + base)
    if rd == "DONOR-SPECIFIC":
        return ("DONOR-SPECIFIC. A transgene's effect does not carry between aged donors even at bulk depth: "
                f"the median reliability-corrected cosine between the same transgene's displacements in "
                f"different donor lines is {fmt(pr['M'])} (95% CI {fmt(pr['ci_lo'])} to {fmt(pr['ci_hi'])}, "
                f"{fmt(pr['angle_M'], 1)} degrees) over {pr['n_q']} qualifying transgenes, and {spec}. "
                "A forward model cannot borrow a factor's direction across donors. The combination screen "
                "needs effects measured per donor and more donors than planned. The same-context control "
                f"recovered M_same = {fmt(g['M_same'])} at single-replicate depth. Measured in three lines "
                "aged 55, 65 and 79, under one dox condition, in one library batch. " + base)
    if rd == "CONSISTENT":
        return (f"CONSISTENT. For these {pr['n_q']} transgenes, the displacement in one aged donor line points "
                f"the same way as in another: {stats}, with {spec}. This holds for these transgenes in these "
                "three lines aged 55, 65 and 79 (of the six perturbed lines aged 55-79) only, under one dox "
                "condition in one library batch. It says nothing about transfer to a different cell type, to "
                "younger donors, or to untested factors. " + base)
    return (f"INTERMEDIATE. The measurement falls between the pre-registered thresholds: {stats} over "
            f"{pr['n_q']} qualifying transgenes, with {spec}. This neither establishes that transgene "
            "directions carry between aged donor lines nor that they do not. " + base)


def findings():
    flag = PREREG.read_text(encoding="utf-8")
    stop = stopped()
    s0 = pd.read_csv(OUT / "stage0_reproduction.csv") if (OUT / "stage0_reproduction.csv").exists() else None
    W = []
    w = W.append
    w("# FINDINGS — DONOR")
    w("")
    if stop:
        w(f"Reading: **UNMEASURABLE** — STOP key `{stop['key']}` at stage `{stop['stage']}`: {stop['message']}")
        w("")
        w("The run stopped. Nothing was computed after the STOP and the run was not resumed.")
        w("")
    g = load("gate.json") if (OUT / "gate.json").exists() else {}
    pr = load("primary.json") if (OUT / "primary.json").exists() else {}
    if pr:
        rd = pr["reading"]
        w(f"Reading: **{rd}**. No STOP fired. The gate passed (M_same = {fmt(g['M_same'])} >= 0.85; gate "
          f"specificity {g['spec_k']} of {g['spec_n']} against a chance baseline of {fmt(g['spec_baseline'], 2)}, "
          f"binomial p = {fmt_p(g['spec_p'])}).")
        branch = []
        if pr.get("ci_valid") and pr["ci_hi"] < UPPER_BAR:
            branch.append(f"the CI upper bound {fmt(pr['ci_hi'])} is below 0.5")
        if not pr["spec_p"] < SPEC_P_KILL:
            branch.append(f"specificity p = {fmt_p(pr['spec_p'])} >= 0.05")
        if rd == "DONOR-SPECIFIC":
            w(f"DONOR-SPECIFIC fired because {' and '.join(branch)}. Its specificity branch did **not** fire: "
              f"specificity is above chance ({pr['spec_k']} of {pr['spec_n']}, baseline "
              f"{fmt(pr['spec_baseline'], 2)}, p = {fmt_p(pr['spec_p'])}). Section 6 says what this does and "
              "does not mean.")
        w("")
    w("## The pre-registration, verbatim")
    w("")
    w("```")
    w(flag.rstrip("\n"))
    w("```")
    w("")
    w(f"`results/donor/PREREG.flag` md5 `{md5(PREREG)}`.")
    w("")

    # ---- deviations ----
    w("## 1. Deviations and choices the flag did not fix")
    w("")
    for item in DEVIATIONS:
        w(f"- {item}")
    w("")

    # ---- stage 0 ----
    w("## 2. Stage 0 — reproduction of SRSFINV (gate)")
    w("")
    if s0 is not None:
        s0j = load("stage0.json")
        w(f"{s0j['n_compared']} quantities compared, {s0j['n_mismatch']} mismatches. Counts and strings exact. "
          f"The {int(s0.quantity.str.startswith('donors.csv ').sum() - 2)} per-cell comparisons of "
          "`results/srsfinv/donors.csv` (19 lines x 10 columns) are in "
          "`results/donor/stage0_reproduction.csv`; the other rows are below.")
        w("")
        w(md(s0[~s0.quantity.str.match(r"donors\.csv (?!n rows|columns)")]))
        w("")
        w(f"Counts file md5 (recorded, no stored value to compare): `{s0j['counts_md5']}`, "
          f"{s0j['counts_bytes']:,} bytes. Basis max|L'L - I| = {s0j['basis_err']:.2e}.")
        w("")
    if stop and stop["stage"] == "stage0":
        FINDINGS.write_text("\n".join(W) + "\n", encoding="utf-8")
        return

    # ---- build ----
    b = load("build.json")
    w("## 3. Ruler mapping, sample rules, cells")
    w("")
    w(f"**Ruler mapping (SRSFINV A13).** The counts file has {b['file_rows']:,} rows (versioned Ensembl IDs). "
      f"On unversioned IDs, **{b['ruler_mapped']:,} of the {b['ruler_genes']:,} ruler genes map**; "
      f"{b['ruler_unmapped']} are unmapped and set to z = 0. {b['file_rows_not_on_ruler']:,} file IDs are not "
      f"ruler genes. Only {b['versioned_exact_matches']:,} file IDs match a ruler `gene_id` with its version, "
      "so matching on versions would have lost most of the ruler. Across the "
      f"{b['panel_gsms']} panel GSMs, {b['z_zero_columns']} ruler columns are all-zero and carry z = 0 "
      f"({b['ruler_unmapped']} unmapped + {b['z_zero_columns_mapped_but_all_zero']} mapped but never counted). "
      f"Median share of a panel library's counts on mapped ruler genes: "
      f"{b['panel_share_of_counts_on_ruler_median']:.4f}; median ruler library size "
      f"{b['panel_library_size_ruler_median']:,.0f}. TMM factors {b['tmm_factor_min']:.3f} to "
      f"{b['tmm_factor_max']:.3f}.")
    w("")
    w(f"Duplicate unversioned IDs: {b['file_rows_duplicate_unversioned']} file rows; "
      f"{b['duplicates_on_ruler']} of them duplicate a ruler gene, and the first row in file order was used "
      f"(flag 4b). Unused copies with any count in the panel: {b['duplicates_on_ruler_nonzero_in_panel']}:")
    w("")
    dd = duplicate_ids()
    w(md(dd[dd.unused_panel_sum > 0]) if (dd.unused_panel_sum > 0).any() else "(none)")
    w("")
    ex = pd.read_csv(OUT / "excluded_gsms.csv")
    w(f"**Sample rules.** Of the 771 GSMs from the six perturbed lines, {771 - len(ex)} pass section 2 and "
      f"{len(ex)} are excluded (one row per GSM in `results/donor/excluded_gsms.csv`; a GSM can carry more than "
      "one reason). Counts by rule:")
    w("")
    cnt = collections.Counter(r.split(" ")[0] for e in ex.exclusion for r in e.split("; "))
    w(md(pd.DataFrame([dict(rule=k, gsm_reasons=v) for k, v in sorted(cnt.items())])))
    w("")
    a7 = ex[ex.exclusion.str.contains("A7|A11")]
    w("GSMs excluded under A7 and A11:")
    w("")
    w(md(a7[["gsm", "title", "ch_genotype", "exclusion"]]))
    w("")
    cells = pd.read_csv(OUT / "cells.csv")
    w(f"**Cells.** {len(cells)} (transgene, line, condition) cells; {int(cells.eligible.sum())} eligible "
      f"(2 transgene + 2 mTagBFP2 GSMs in one batch cell; every eligible cell has exactly 2 + 2, so every split "
      f"is the 2-pairing enumeration). Primary eligible cells by batch: {b['eligible_primary_batches']}, "
      f"platform {b['eligible_primary_platforms']}, all polyA. Ineligible cells:")
    w("")
    w(md(cells[~cells.eligible][["transgene", "line", "cond", "bcells_available", "n_t", "n_c", "reason"]]))
    w("")
    prim = cells[cells.eligible & (cells.cond == PRIMARY)]
    nl = prim.groupby("transgene").line.apply(lambda x: "/".join(str(LINES[l]) for l in LINE_ORDER if l in set(x)))
    w(f"Primary candidates (eligible in >= 2 lines): {int((prim.groupby('transgene').size() >= 2).sum())}. "
      "Line sets: " + ", ".join(f"{k}: {v}" for k, v in nl.value_counts().items()) + ". Single-line, not "
      "candidates: " + ", ".join(sorted(nl[~nl.str.contains('/')].index)) + ".")
    w("")
    w("Every line-65 primary cell uses the same two mTagBFP2 GSMs, one of which (GSM9611762) carries the "
      "description flag `High GC`; kept under rule R1. Flagged GSMs per cell are in `cells.csv`.")
    w("")
    if stop and stop["stage"] == "read" and not pr:
        w("## 4. Gate")
        w("")
        gd = pd.read_csv(OUT / "gate_table.csv")
        w(md(gd))
        FINDINGS.write_text("\n".join(W) + "\n", encoding="utf-8")
        return

    # ---- gate ----
    gd = pd.read_csv(OUT / "gate_table.csv")
    w("## 4. Qualification gate (same-context control, uncorrected)")
    w("")
    w(md(pd.DataFrame([
        dict(quantity="candidate transgenes", value=g["n_candidates"], bar=""),
        dict(quantity="M_same = median c_same (single-replicate halves)", value=fmt(g["M_same"], 4), bar=">= 0.85"),
        dict(quantity="angle of M_same (degrees)", value=fmt(g["angle_M_same"], 1), bar=""),
        dict(quantity="c_same interquartile range", value=f"{fmt(g['q25'])} to {fmt(g['q75'])}", bar=""),
        dict(quantity="specific / testable", value=f"{g['spec_k']} / {g['spec_n']}", bar=""),
        dict(quantity="chance baseline (5%)", value=fmt(g["spec_baseline"], 2), bar=""),
        dict(quantity="binomial p (one-sided, approximate)", value=fmt_p(g["spec_p"]), bar="< 0.01"),
        dict(quantity="gate", value="PASS" if g["passed"] else "FAIL", bar=""),
    ])))
    w("")
    w("This is the recovery check this project requires beside any null: at single-replicate depth, inside "
      "one line and one batch, the statistic recovers agreement that is 1 by construction at "
      f"{fmt(g['M_same'])} ({fmt(g['angle_M_same'], 1)} degrees) and identifies {g['spec_k']} of "
      f"{g['spec_n']} transgenes. Per-transgene gate table (`results/donor/gate_table.csv`):")
    w("")
    w(md(gd, ["transgene", "lines", "c_same", "angle_same", "n_comp", "p95_comp", "specific"]))
    w("")

    # ---- primary ----
    t = pd.read_csv(OUT / "primary_transgenes.csv")
    pp = pd.read_csv(OUT / "primary_pairs.csv")
    pc = pd.read_csv(OUT / "primary_cells.csv")
    w("## 5. Primary readout (Dox_1000ng/mL, 3days; lines 55, 65, 79; batch LIB052667)")
    w("")
    w(md(pd.DataFrame([
        dict(quantity="candidate transgenes", value=pr["n_candidates"]),
        dict(quantity="eligible cells / qualifying cells (r >= 0.3)", value=f"{pr['n_cells']} / {pr['n_cells_qualifying']}"),
        dict(quantity="qualifying transgenes n_q", value=pr["n_q"]),
        dict(quantity="non-finite c_corr", value=pr["n_nonfinite"]),
        dict(quantity="M = median c_corr", value=fmt(pr["M"], 4)),
        dict(quantity="95% bootstrap CI (10,000, over transgenes)", value=f"{fmt(pr['ci_lo'], 4)} to {fmt(pr['ci_hi'], 4)}"),
        dict(quantity="CI contains M", value=pr["ci_valid"]),
        dict(quantity="arccos(M), degrees", value=fmt(pr["angle_M"], 1)),
        dict(quantity="specific / testable", value=f"{pr['spec_k']} / {pr['spec_n']}"),
        dict(quantity="chance baseline (5%)", value=fmt(pr["spec_baseline"], 2)),
        dict(quantity="binomial p (one-sided, approximate)", value=fmt_p(pr["spec_p"])),
        dict(quantity="reading", value=pr["reading"]),
    ])))
    w("")
    w(f"Reliability floor, per line (c_half = cosine between single-replicate halves; floor angle = "
      f"arccos c_half). Built from transgene replicate splits, not from control-versus-control displacements:")
    w("")
    fl = pc.groupby("age").agg(cells=("c_half", "size"), median_c_half=("c_half", "median"),
                               q10_c_half=("c_half", lambda x: x.quantile(0.1)),
                               median_r=("r", "median"), qualifying=("qualifying", "sum")).reset_index()
    fl["median_floor_angle"] = fl.median_c_half.map(angle)
    w(md(fl))
    w("")
    w("Per-transgene table (`results/donor/primary_transgenes.csv`; pairs in `primary_pairs.csv`, cells in "
      "`primary_cells.csv`). c_corr is the mean over the transgene's qualifying line pairs; s is the mean "
      "uncorrected cosine used for specificity; p95_comp the 95th percentile of its comparison set of size n_comp.")
    w("")
    t["angle_corr"] = t.c_corr.map(angle)
    w(md(t.sort_values("transgene"), ["transgene", "lines", "n_qual_pairs", "c_obs_mean", "c_corr", "angle_corr",
                                      "s", "n_comp", "p95_comp", "specific"]))
    w("")

    # ---- reading ----
    w("## 6. The reading that fired, verbatim")
    w("")
    w("> " + reading_text(pr["reading"], pr.get("reading_detail", ""), g, pr))
    w("")
    w("What the numbers beside it say (descriptive, not gated). The reading fired on its CI branch. The "
      "measured agreement is not zero: the median corrected cosine is "
      f"{fmt(pr['M'])} ({fmt(pr['angle_M'], 1)} degrees), and {pr['spec_k']} of {pr['spec_n']} transgenes "
      f"are specific against {fmt(pr['spec_baseline'], 2)} expected by chance. So the directions are partly "
      "shared between lines, and transgene-specifically so. What is ruled out, at the pre-registered bar, is "
      "agreement close enough to borrow a direction across lines. Within one line the same transgene's "
      f"replicate halves agree at {fmt(g['M_same'])}, so the shortfall is not measurement noise at this "
      "depth. The wording's \"does not carry\" is stronger than these numbers. It is reproduced because it is "
      "the pre-registered text of the branch that fired.")
    w("")
    w("Four things that limit how far the reading reaches:")
    w("")
    lp = pp[pp.qualifying_pair].assign(pair=lambda d: d.age_A.astype(str) + "-" + d.age_B.astype(str)) \
        .groupby("pair").agg(n_pairs=("c_corr", "size"), median_c_obs=("c_obs", "median"),
                             median_c_corr=("c_corr", "median")).reset_index()
    lp["angle_median_c_corr"] = lp.median_c_corr.map(angle)
    d = load("diagnostics.json")
    w(f"1. **Line pairs differ a lot.** Median corrected cosine by line pair: " +
      "; ".join(f"{r.pair}: {fmt(r.median_c_corr)} over {r.n_pairs} pairs" for r in lp.itertuples()) +
      ". Line sets per transgene: " + ", ".join(f"{k} {v}" for k, v in t.lines.value_counts().items()) +
      ". So the median over transgenes partly reflects which line pairs each transgene happens to have, "
      "and that is an accident of the design.")
    w("")
    w(md(lp))
    w("")
    d2 = pd.DataFrame(d["D2"]["summary"])
    w("2. **Re-measuring the same line moves the direction too** (D2, diagnostic). For the 9 transgenes measured "
      "under both doses, the same line re-measured in another batch at another dose agrees at median corrected "
      + fmt(d2.set_index("category").loc["same line / different batch", "median_c_corr"]) +
      ", below two different lines in the same batch ("
      + fmt(d2.set_index("category").loc["different line / same batch", "median_c_corr"]) +
      "). Batch and dose are confounded there, so this does not separate them. But re-measuring one line "
      "in another batch at another dose moved the direction at least as much as changing the line did. "
      "\"Donor-specific\" here cannot be told apart from \"specific to the experiment\".")
    w("")
    w("3. **What Rep1 and Rep2 are is not stated.** The floor is the agreement between the record's Rep1 and "
      "Rep2 GSMs in one line and batch. If those are technical replicates or wells of one culture, the floor "
      "is set by technical noise only, and culture-to-culture variation within a donor ends up in the "
      "cross-line drift. Item 2 is consistent with that. The gate's 0.934 is therefore a recovery check at "
      "replicate depth, not evidence that a second culture of the same donor would agree at 0.934.")
    w("")
    w(f"4. **The 20-component space holds little of each displacement** (D6): median share of squared norm "
      f"retained {fmt(d['D6']['median_share'])} (5-95%: {fmt(d['D6']['q05'])} to {fmt(d['D6']['q95'])}). "
      "The angles are measured in a GTEx-derived subspace that holds about an eighth of each displacement. "
      "The basis was chosen for the Southard screen, not for this data.")
    w("")

    # ---- diagnostics ----
    w("## 7. Diagnostics (each one diagnostic, not a gated result)")
    w("")
    r1 = d["D1"]
    w(f"**D1 — line mean removed** (diagnostic, not a gated result). n_q {r1['n_q']}; M {fmt(r1['M'], 4)} "
      f"(95% CI {fmt(r1['ci_lo'], 4)} to {fmt(r1['ci_hi'], 4)}, {fmt(r1['angle_M'], 1)} degrees); specificity "
      f"{r1['spec_k']} of {r1['spec_n']} against {fmt(r1['spec_baseline'], 2)}, p = {fmt_p(r1['spec_p'])}. "
      "Removing each line's mean displacement lowers agreement and raises specificity. Part of the primary "
      "agreement is a component that all transgenes in a line share against mTagBFP2. Disclosed: the line "
      "mean includes the transgene itself, so this is not independent of the vectors it corrects. Tables "
      "`d1_line_mean_removed_*.csv`.")
    w("")
    w("**D2 — batch versus donor** (diagnostic, not a gated result). In the primary all three lines share batch "
      "LIB052667, so primary donor differences cannot track batch. Across conditions, for the transgenes "
      f"present under both doses ({', '.join(d['D2']['transgenes'])}):")
    w("")
    w(md(d2))
    w("")
    d3 = d["D3"]
    w(f"**D3 — agreement and replicate count** (diagnostic, not a gated result). Every eligible primary cell has "
      f"{d3['replicate_counts_in_primary_cells']} transgene and {d3['control_counts']} control GSMs, so "
      "replicate count does not vary and its relation to agreement is undefined. Spearman of c_corr(t) with "
      f"the mean library size of t's GSMs: rho = {fmt(d3['rho_lib'])}, p = {fmt_p(d3['p_lib'])}; with t's "
      f"number of qualifying lines: rho = {fmt(d3['rho_nlines'])}, p = {fmt_p(d3['p_nlines'])} (n = {d3['n']}). "
      "No transgene shortlist is reported, so there is no shortlist to check against cell count.")
    w("")
    d4 = d["D4"]
    w(f"**D4 — full gene space** (diagnostic, not a gated result). Over {d4['n_pairs']} qualifying pairs, median "
      f"cosine in the 23,485-column z space {fmt(d4['median_c_obs_full'])} ({fmt(d4['median_angle_full'], 1)} "
      f"degrees) against {fmt(d4['median_c_obs_20'])} ({fmt(d4['median_angle_20'], 1)} degrees) in the "
      "20-component space. That is close to orthogonal, as TRANSFER expected: shared signal is diluted by "
      "noise spread over 23,485 dimensions. The flag said the expectation was weaker at bulk depth; it held "
      "anyway. Table `d4_full_gene_space_pairs.csv`.")
    w("")
    t5 = pd.read_csv(OUT / "d5_dox2000_transgenes.csv")
    t5["angle_corr"] = t5.c_corr.map(angle)
    r5 = d["D5"]
    w(f"**D5 — Dox 2000 six-line set** (diagnostic, not a gated result). {r5['n_q']} qualifying transgenes "
      f"(fewer than 20, so it could not be a reading); median c_corr {fmt(r5['M'])} "
      f"({fmt(r5['angle_M'], 1)} degrees); no specificity test (comparison sets below 19). Two batches: "
      "55/65/79 in LIB053907, 67/68/69 in LIB054839, so most of its pairs cross a batch (see D2).")
    w("")
    w(md(t5, ["transgene", "lines", "n_qual_pairs", "c_obs_mean", "c_corr", "angle_corr"]))
    w("")
    w(f"**D6 — norm bookkeeping** (diagnostic, not a gated result). Share of each displacement's squared "
      f"z-space norm retained by the projection: median {fmt(d['D6']['median_share'])}, 5% "
      f"{fmt(d['D6']['q05'])}, 95% {fmt(d['D6']['q95'])}. Per cell in `d6_norms.csv`.")
    w("")
    w("**D7 — tables.** `cells.csv` (lines, GSMs, controls, passages, flags, platforms, library size), "
      "`primary_cells.csv` (c_half, floor angle, r), `primary_pairs.csv` (c_obs, angle, r_A, r_B, c_corr), "
      "`primary_transgenes.csv`, `gate_table.csv`, `panel_gsms.csv`.")
    w("")

    # ---- section 10 verbatim ----
    w("## 8. What this does not cover (flag section 10, verbatim)")
    w("")
    m = re.search(r"--- 10\. WHAT THIS DOES NOT COVER.*?---\n(.*?)\n--- STOP KEYS", flag, re.S)
    w("```")
    w(m.group(1).rstrip() if m else "(section 10 not found in the flag)")
    w("```")
    w("")
    w("Added here. None of these lines is young or neonatal, and nothing here transfers to the Southard "
      "(Hs27) or Sengstack cells. Of the six lines, 67, 68 and 69 enter only D5.")
    w("")
    w("## 9. Files")
    w("")
    w("`src/donor_run.py`; `results/donor/` (PREREG.flag, run.log, stage0_reproduction.csv, stage0.json, "
      "build.json, excluded_gsms.csv, cells.csv, panel_gsms.csv, gate.json, gate_table.csv, primary.json, "
      "primary_*.csv, diagnostics.json, d1/d2/d4/d5/d6 CSVs, displacements.npz). "
      "`results/donor/displacements_full.npy` (26.6 MB, the full z-space displacements) is git-ignored and "
      "regenerated by `python src/donor_run.py build`. Raw inputs, git-ignored under `data/raw/srsfinv/`: "
      f"the family SOFT (SRSFINV) and `{COUNTS.name}` ({COUNTS_URL}).")
    FINDINGS.write_text("\n".join(W) + "\n", encoding="utf-8")
    log(f"[findings] wrote {FINDINGS}")


DEVIATIONS = [
    "**The gate is uncorrected.** The task text gates on the corrected median. At 2 replicates per cell the "
    "gate's halves are single GSMs, and correcting them would divide the split cosine by itself (1 by "
    "construction). The author chose \"Uncorrected, bar 0.85\" in chat before the flag was written. The "
    "gate is therefore stricter than a corrected gate on the same halves, and it is computed without the "
    "r >= 0.3 filter (flag section 6).",
    "**The primary is one condition in three lines, not six.** The author chose \"Dox1000 only\" in chat before "
    "the flag. The task's fixed CONSISTENT wording said \"these six lines aged 55-79\"; the flag's wording says "
    "\"these three lines aged 55, 65 and 79 (of the six perturbed lines aged 55-79)\". Lines 67, 68 and 69 enter "
    "only diagnostic D5.",
    "**Readings carry added sentences.** Each branch keeps the elements the task required and adds the gate "
    "result, the line and batch scope, and the verbatim aged-donors sentence (flag section 7).",
    "**What the task did not specify, fixed in the flag:** the control (mTagBFP2 in the same line, batch, "
    "treatment and sub-batch); the batch-cell rule R2; eligibility R3; the reliability threshold r >= 0.3 "
    "(TRANSFER's); per-transgene c_corr as the mean over its qualifying line pairs; the specificity "
    "construction (comparison set G_t, |G_t| >= 19, strict 95th percentile); z(group) as the mean of per-GSM z "
    "vectors rather than a z of pooled counts; one joint TMM over every GSM in an eligible cell of either "
    "condition; duplicate Ensembl IDs resolved to the first row in file order; the D1-D7 constructions.",
    "**STOP keys.** The file-size, header, SOFT-md5 and basis-orthonormality checks are folded into "
    "`prior_results_not_reproduced` instead of separate keys. None fired.",
    "**Before the flag**, the SRSFINV tracked CSVs were cross-tabulated (line x batch x treatment x genotype x "
    "Rep x description flag) to settle the metadata and find the replicate structure. No expression value "
    "existed then. The flag file was created with a timestamp placeholder and revised three times before it "
    "was stamped (panel definition narrowed to eligible cells; a parenthetical claim about platforms removed; "
    "D2 made conditional on the batch check). The stamp went in at 2026-10-21T22:17:04Z, before the counts "
    "matrix was downloaded. It has not been edited since.",
    "**Download.** `GSE325735_human_raw_counts.txt.gz` was fetched with curl into `data/raw/srsfinv/` "
    "(SRSFINV's raw directory, git-ignored), after the flag. Size 26,007,361 bytes matches SRSFINV's HEAD; "
    "md5 8e3ecb46366730a5f6479147bf57895f. `fetch_counts()` (urllib, same URL) was added to the script after "
    "the run so the committed code can fetch it. It does not run when the file exists.",
    "**Imports.** `toward_run` (tmm_logcpm_quiet, zscore_frozen, percentile_ci), `fibro2_common` "
    "(load_frozen_ruler), `same_run` (ci_covers_point), none modified. Importing `toward_run` runs its "
    "module-level `mkdir` of `results/toward/`, which already exists; nothing is written there.",
    "**Order of work.** The findings writer (`findings()`, `md()`, `reading_text()`, `duplicate_ids()`, this "
    "list) and `fetch_counts()` were written after `read` had produced its numbers. They format stored outputs "
    "and compute no statistic. `duplicate_ids()` re-reads the raw matrix only to list which duplicated IDs "
    "carry counts. No analysis code was changed after any DONOR number existed, and no stage was rerun.",
    "**SUMO1 in line 65** has dox GSMs in two batch cells (LIB051945, LIB054839), 2 each. R2's tie rule "
    "excludes the cell. Neither batch has a line-65 mTagBFP2 at that treatment, so it would have been "
    "ineligible anyway. CYCS and UBE2I in line 65 (LIB053449) have one mTagBFP2 GSM, so they are ineligible "
    "under R3.",
    "**OSKM** enters the primary in lines 55 and 65 only. Its two line-79 Dox 1000 GSMs carry `genotype: SOX2` "
    "and are excluded under A7.",
    "**Files outside results/donor/.** `.gitignore` gained one entry, `results/donor/displacements_full.npy` "
    "(26.6 MB, regenerable), with a comment naming its source. `FINDINGS_DONOR.md` was regenerated more than once by "
    "the `findings` stage while its wording was revised; the numbers come from the single `read` run.",
    "**Seeds.** Every eligible cell is 2 + 2, so no random split was drawn. The only random draw is the "
    "bootstrap (default_rng(20260914), 10,000 resamples), used for the primary and D1.",
]


def main():
    require_prereg()
    OUT.mkdir(parents=True, exist_ok=True)
    args = sys.argv[1:] or ["all"]
    todo = STAGES if args == ["all"] else args
    for st in todo:
        if st not in STAGES:
            raise SystemExit(f"unknown stage {st}")
        s = stopped()
        if s and st != "findings":
            log(f"[{st}] not run: STOP {s['key']} already recorded. Never resumed.")
            continue
        log(f"===== {st} {time.strftime('%Y-%m-%dT%H:%M:%SZ', time.gmtime())}")
        try:
            {"stage0": stage0, "build": build, "read": read, "findings": findings}[st]()
        except Stop as e:
            log(f"STOP {e.key}: {e.message}")
            dump("STOP.json", dict(key=e.key, message=e.message, stage=st, extra=e.extra))


if __name__ == "__main__":
    main()
