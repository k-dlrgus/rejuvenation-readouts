"""Stage 2 — project frozen GTEx fibroblast ridge onto GSE297234.

Gated on Stage 1 pass. Writes PREREG_STAGE2.flag before the first projection.
Nothing is fitted on GSE297234. Does not open GSE325735.
"""
from __future__ import annotations

import gzip
import io
import re
import sys
import tarfile
import time
from pathlib import Path

import h5py
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from scipy import sparse
from scipy.io import mmread
from sklearn.cluster import KMeans
from sklearn.decomposition import PCA

sys.path.insert(0, str(Path(__file__).resolve().parent))
from fibro_common import (  # noqa: E402
    FIBRO_DIR, FIBRO_FIG, FIBRO_RAW, FIBRO_PROC, FIBRO_SEED, BOOT_SEED,
    N_PERM, N_RANDOM_DIR, SMTSD_FIBRO, FROZEN_RULER, PREREG_STAGE2_FLAG,
    PREREG_STAGE2, PLURI_ENDOGENOUS, FIBRO_IDENTITY, OSKM_FAMILY,
    GSE_ACCESSION, GSE_EXPECTED_DONORS, GSE_EXPECTED_DAYS, GSE_SUPPL,
    StopStep, Logger, dump_json, load_json, strip_ensembl, spearman_safe,
    tmm_norm_factors, log2_cpm_edger, unit, jsonable, load_manifest,
    save_manifest, record_failure, write_prereg_stage2, progress_snapshot,
    fibro_log_banner,
)
from fibro_stage1 import freeze_ruler  # noqa: E402
from download import download as http_download  # noqa: E402
from geo_meta import fetch_soft, parse_series  # noqa: E402
from trajectory_common import permutation_p as _perm_p  # noqa: E402
from gtex_common import EDGE_R_PRIOR  # noqa: E402
from scipy import stats as _stats


def spearman_traj(a, b):
    """Spearman on a short time course. gtex spearman_safe requires n≥8 (donor CV)."""
    a = np.asarray(a, float)
    b = np.asarray(b, float)
    m = np.isfinite(a) & np.isfinite(b)
    if int(m.sum()) < 3:
        return np.nan
    a, b = a[m], b[m]
    if float(np.std(a)) < 1e-12 or float(np.std(b)) < 1e-12:
        return np.nan
    r, _ = _stats.spearmanr(a, b)
    return float(r) if np.isfinite(r) else np.nan


SEV_RE = re.compile(
    r"sendai|\bsev\b|\bsev[-_]|kozak|transgene|cytotune|cyto.?tune|"
    r"vector|orfeome|\boskm\b|sevoskm",
    re.I,
)
CLUSTER_MIN_CELLS = 50
CLUSTER_K = 3
CLUSTER_NPC = 20
# RAW.tar holds the 10x matrices. The two .rds files are Seurat objects; this
# pipeline does not load them (no R / pyreadr) and does not need them.
GSE_SUPPL_REQUIRED = (
    "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE297nnn/GSE297234/suppl/GSE297234_RAW.tar",
)
AGED_LINE = "GM00731"
YOUNG_LINE = "GM23815"


def _gate_pass():
    p = FIBRO_DIR / "stage1_gate.json"
    if not p.exists():
        raise StopStep("stage2_gate", "stage1_gate.json missing — run Stage 1")
    rec = load_json(p)
    if not rec.get("pass"):
        raise StopStep(
            "stage2_gate",
            "Stage 1 failed. Refusing to open GSE297234. "
            + str(rec.get("reason", "")),
        )
    return rec


def fetch_gse_sample_table(log):
    log(f"[geo] fetching SOFT for {GSE_ACCESSION}")
    text = fetch_soft(GSE_ACCESSION, targ_gsm=False)
    if "!Series_title" not in text:
        raise StopStep("geo", f"{GSE_ACCESSION} SOFT has no Series_title. head={text[:300]!r}")
    series = parse_series(text)
    log(f"[geo] Series_title: {series.get('Series_title', [''])[0]}")
    log(f"[geo] Series_overall_design: {(' '.join(series.get('Series_overall_design', [])))[:1500]}")
    log(f"[geo] Series_supplementary_file: {series.get('Series_supplementary_file', [])}")
    log(f"[geo] Series_sample_id: {series.get('Series_sample_id', [])}")
    log(f"[geo] series fields actually read: {sorted(series.keys())}")
    gsm_text = fetch_soft(GSE_ACCESSION, targ_gsm=True)
    samples = _parse_gsm_block(gsm_text, log)
    dump_json(FIBRO_DIR / "stage2_geo_series.json", jsonable({
        k: series.get(k) for k in (
            "Series_title", "Series_overall_design", "Series_summary",
            "Series_supplementary_file", "Series_sample_id", "Series_pubmed_id",
            "Series_type", "Series_platform_id", "Series_sample_organism",
        )
    }))
    samples.to_csv(FIBRO_DIR / "stage2_geo_samples.csv", index=False)
    (FIBRO_DIR / "stage2_columns_read_geo.txt").write_text(
        "series_fields:\n" + "\n".join(sorted(series.keys())) + "\n\n"
        + "gsm_fields_union:\n" + "\n".join(sorted(samples.columns.astype(str))) + "\n",
        encoding="utf-8",
    )
    return series, samples


def _parse_gsm_block(text: str, log) -> pd.DataFrame:
    recs, cur = [], {}
    for line in text.splitlines():
        if line.startswith("^SAMPLE"):
            if cur:
                recs.append(cur)
            cur = dict(gsm=line.split("=")[-1].strip())
            continue
        if not line.startswith("!Sample_"):
            continue
        k, _, v = line[1:].partition(" = ")
        k, v = k.strip(), v.strip()
        cur.setdefault(k, [])
        cur[k].append(v)
    if cur:
        recs.append(cur)
    if not recs:
        raise StopStep("geo", "no GSM sample blocks in SOFT")
    rows = []
    for rec in recs:
        chars = rec.get("Sample_characteristics_ch1", [])
        char_map = {}
        for c in chars:
            if ":" in c:
                a, b = c.split(":", 1)
                char_map[a.strip().lower()] = b.strip()
            else:
                char_map.setdefault("_unparsed", []).append(c)
        growth = " ".join(rec.get("Sample_growth_protocol", []))
        title = " ".join(rec.get("Sample_title", []))
        treatment = char_map.get("treatment", "")
        day = _parse_day(treatment) or _parse_day(title)
        cell_line = char_map.get("cell line") or char_map.get("cell_line") or _parse_line(title)
        age_from_char = char_map.get("age") or char_map.get("donor age") or char_map.get("age (years)")
        age_from_growth = _parse_age_years(growth)
        rows.append(dict(
            gsm=rec.get("gsm"),
            title=title,
            source="; ".join(rec.get("Sample_source_name_ch1", [])),
            organism="; ".join(rec.get("Sample_organism_ch1", [])),
            molecule="; ".join(rec.get("Sample_molecule_ch1", [])),
            description="; ".join(rec.get("Sample_description", [])),
            characteristics="; ".join(chars),
            growth_protocol=growth,
            cell_line=cell_line,
            cell_type=char_map.get("cell type"),
            tissue=char_map.get("tissue"),
            treatment=treatment,
            day=day,
            age_from_characteristics=age_from_char,
            age_from_growth_protocol=age_from_growth,
            library_name=_parse_library(rec.get("Sample_description", [])),
            fields_read=sorted(rec.keys()),
        ))
    df = pd.DataFrame(rows)
    log(f"[geo] parsed {len(df)} GSM rows. columns actually read: {list(df.columns)}")
    log(df[["gsm", "title", "cell_line", "day", "age_from_characteristics",
            "age_from_growth_protocol"]].to_string(index=False))
    return df


def _parse_day(text):
    if text is None or (isinstance(text, float) and not np.isfinite(text)):
        return np.nan
    s = str(text)
    m = re.search(r"\bday\s*(\d+)\b", s, re.I)
    if m:
        return int(m.group(1))
    m = re.search(r"\bD(\d+)\b", s)
    if m:
        return int(m.group(1))
    return np.nan


def _parse_line(text):
    m = re.search(r"(GM\d+)", str(text), re.I)
    return m.group(1).upper() if m else None


def _parse_age_years(text):
    s = str(text)
    m = re.search(r"(\d+)\s*[-\s]?year", s, re.I)
    if m:
        return int(m.group(1))
    m = re.search(r"\((\d+)\s*yr\)", s, re.I)
    if m:
        return int(m.group(1))
    return np.nan


def _parse_library(desc_list):
    for d in desc_list:
        m = re.search(r"Library name:\s*(\S+)", str(d), re.I)
        if m:
            return m.group(1)
    return None


def verify_sample_table(series, samples, log):
    discrepancies = []
    design = " ".join(series.get("Series_overall_design", []))
    n = int(len(samples))
    if n != 8:
        discrepancies.append(f"n_samples={n} (expected 8)")
    days = sorted({int(d) for d in samples.day.dropna().astype(int).tolist()})
    if tuple(days) != GSE_EXPECTED_DAYS:
        discrepancies.append(f"days observed={days} expected={list(GSE_EXPECTED_DAYS)}")
    lines = sorted(set(samples.cell_line.dropna().astype(str).str.upper()))
    expected_lines = sorted(d["cell_line"] for d in GSE_EXPECTED_DONORS)
    if lines != expected_lines:
        discrepancies.append(f"cell_line observed={lines} expected={expected_lines}")
    # ages: characteristics may be empty; use the record field that actually has years
    age_map = {}
    for _, r in samples.iterrows():
        cl = str(r.cell_line).upper() if pd.notna(r.cell_line) else None
        if not cl:
            continue
        age = r.age_from_characteristics
        src = "Sample_characteristics_ch1"
        if age is None or (isinstance(age, float) and not np.isfinite(age)) or str(age).strip() in ("", "nan"):
            age = r.age_from_growth_protocol
            src = "Sample_growth_protocol"
        if age is None or (isinstance(age, float) and not np.isfinite(age)):
            m = re.search(re.escape(cl) + r"[^0-9]{0,40}(\d+)\s*yr", design, re.I)
            if m:
                age = int(m.group(1))
                src = "Series_overall_design"
        if cl not in age_map and age is not None and str(age) not in ("", "nan"):
            try:
                age_map[cl] = dict(age_years=int(float(age)), source=src)
            except (TypeError, ValueError):
                discrepancies.append(f"unparsable age for {cl}: {age!r} from {src}")
    expected_age = {d["cell_line"]: d["age_years"] for d in GSE_EXPECTED_DONORS}
    for cl, exp in expected_age.items():
        got = age_map.get(cl, {}).get("age_years")
        if got != exp:
            discrepancies.append(f"{cl} age in record={got} (source={age_map.get(cl)}) expected={exp}")
    # attach ages from the record onto samples — never inferred beyond parsed fields
    samples = samples.copy()
    samples["age_years"] = samples.cell_line.astype(str).str.upper().map(
        lambda c: age_map.get(c, {}).get("age_years")
    )
    samples["age_source"] = samples.cell_line.astype(str).str.upper().map(
        lambda c: age_map.get(c, {}).get("source")
    )
    if samples["age_years"].isna().any():
        bad = samples.loc[samples.age_years.isna(), "gsm"].tolist()
        raise StopStep(
            "age_record",
            f"donor age not in the record for GSM {bad}. Not inferring.",
            dict(age_map=age_map, discrepancies=discrepancies),
        )
    rec = dict(
        n_samples=n, days=days, cell_lines=lines, age_map=age_map,
        discrepancies=discrepancies, used_record=True,
        overall_design=design,
    )
    dump_json(FIBRO_DIR / "stage2_sample_verify.json", jsonable(rec))
    if discrepancies:
        log(f"[geo] DISCREPANCY vs scoping description: {discrepancies}. Using the record.")
    else:
        log("[geo] sample table matches scoping description (2 donors, days 0/3/7/10, ages 22 and 96).")
    samples.to_csv(FIBRO_DIR / "stage2_geo_samples.csv", index=False)
    return samples, rec


def download_suppl(log):
    FIBRO_RAW.mkdir(parents=True, exist_ok=True)
    raw_tar = FIBRO_RAW / "GSE297234_RAW.tar"
    dest = FIBRO_RAW / "RAW"
    h5_existing = list(dest.glob("*filtered_feature_bc_matrix.h5")) if dest.exists() else []
    if raw_tar.exists() and len(h5_existing) >= 8:
        log(f"[download] RAW.tar already on disk size={raw_tar.stat().st_size}; "
            f"{len(h5_existing)} h5 already extracted. Not re-fetching FTP.")
        got = {
            GSE_SUPPL_REQUIRED[0]: dict(
                path=str(raw_tar), size=int(raw_tar.stat().st_size), cached=True
            )
        }
        for url in GSE_SUPPL:
            if url not in got:
                got[url] = dict(skipped=True, reason="Seurat RDS not required; 10x matrices are in RAW.tar")
                log(f"[download] skip (Seurat RDS not loaded): {url}")
        dump_json(FIBRO_DIR / "stage2_downloads.json", jsonable(got))
        return raw_tar, got
    got = {}
    skip = [u for u in GSE_SUPPL if u not in GSE_SUPPL_REQUIRED]
    for url in skip:
        log(f"[download] skip (Seurat RDS not loaded; counts are in RAW.tar): {url}")
        got[url] = dict(skipped=True, reason="Seurat RDS not required; 10x matrices are in RAW.tar")
    for url in GSE_SUPPL_REQUIRED:
        log(f"[download] {url}")
        try:
            path = Path(http_download(url, str(FIBRO_RAW)))
        except Exception as e:
            log(f"[download] FAIL {url}: {type(e).__name__}: {e}")
            record_failure("download", f"{url}: {type(e).__name__}: {e}")
            continue
        got[url] = dict(path=str(path), size=int(path.stat().st_size) if path.exists() else 0)
        log(f"[download] {path.name} size={got[url]['size']}")
    dump_json(FIBRO_DIR / "stage2_downloads.json", jsonable(got))
    if not raw_tar.exists():
        raise StopStep("download", f"GSE297234_RAW.tar missing after download. got={list(got)}")
    return raw_tar, got


def extract_tar(tar_path: Path, log) -> Path:
    dest = FIBRO_RAW / "RAW"
    dest.mkdir(parents=True, exist_ok=True)
    marker = dest / "_extracted.ok"
    if marker.exists() and any(dest.iterdir()):
        log(f"[tar] already extracted in {dest}")
        return dest
    log(f"[tar] extracting {tar_path}")
    with tarfile.open(tar_path, "r") as tf:
        names = tf.getnames()
        log(f"[tar] n_members={len(names)} example={names[:20]}")
        tf.extractall(dest)
    marker.write_text("ok\n", encoding="utf-8")
    (FIBRO_DIR / "stage2_tar_members.txt").write_text("\n".join(names) + "\n", encoding="utf-8")
    return dest


def _gsm_from_name(name: str):
    m = re.search(r"(GSM\d+)", name)
    return m.group(1) if m else None


def find_mtx_triplets(extract: Path, log):
    mats = list(extract.rglob("*matrix.mtx.gz")) + list(extract.rglob("*matrix.mtx"))
    log(f"[10x] matrix.mtx files: {len(mats)}")
    trips = []
    for mat in mats:
        stem = mat.name.replace("_matrix.mtx.gz", "").replace("_matrix.mtx", "")
        feat = mat.with_name(stem + "_features.tsv.gz")
        if not feat.exists():
            feat = mat.with_name(stem + "_features.tsv")
        if not feat.exists():
            feat = mat.with_name(stem + "_genes.tsv.gz")
        bc = mat.with_name(stem + "_barcodes.tsv.gz")
        if not bc.exists():
            bc = mat.with_name(stem + "_barcodes.tsv")
        # also try sibling filtered_feature_bc_matrix
        if (not feat.exists() or not bc.exists()) and mat.parent.name.lower().startswith("filtered"):
            feat = mat.parent / "features.tsv.gz"
            bc = mat.parent / "barcodes.tsv.gz"
        if not feat.exists() or not bc.exists():
            log(f"   skip {mat.name}: missing features/barcodes "
                f"(feat={feat.exists()} bc={bc.exists()})")
            continue
        trips.append(dict(matrix=mat, features=feat, barcodes=bc, stem=stem,
                          gsm=_gsm_from_name(mat.name) or _gsm_from_name(stem)))
    if not trips:
        # directory-style 10x
        for feat in list(extract.rglob("features.tsv.gz")) + list(extract.rglob("genes.tsv.gz")):
            parent = feat.parent
            mat = parent / "matrix.mtx.gz"
            if not mat.exists():
                mat = parent / "matrix.mtx"
            bc = parent / "barcodes.tsv.gz"
            if not bc.exists():
                bc = parent / "barcodes.tsv"
            if mat.exists() and bc.exists():
                trips.append(dict(matrix=mat, features=feat, barcodes=bc, stem=parent.name,
                                  gsm=_gsm_from_name(parent.name) or _gsm_from_name(str(parent))))
    log(f"[10x] usable triplets: {len(trips)}")
    dump_json(FIBRO_DIR / "stage2_mtx_triplets.json", jsonable(
        [{k: str(v) if isinstance(v, Path) else v for k, v in t.items()} for t in trips]
    ))
    return trips


def _h5_decode(arr) -> np.ndarray:
    out = []
    for x in arr:
        if isinstance(x, (bytes, np.bytes_)):
            out.append(x.decode("utf-8"))
        else:
            out.append(str(x))
    return np.asarray(out, dtype=object)


def find_10x_h5(extract: Path, log):
    files = sorted(extract.rglob("*filtered_feature_bc_matrix.h5")) + sorted(
        extract.rglob("*feature_bc_matrix.h5")
    )
    seen = set()
    out = []
    for p in files:
        if p.resolve() in seen:
            continue
        seen.add(p.resolve())
        gsm = _gsm_from_name(p.name)
        out.append(dict(path=p, gsm=gsm))
        log(f"   h5 {p.name} gsm={gsm} size={p.stat().st_size}")
    log(f"[10x] h5 matrices: {len(out)}")
    dump_json(FIBRO_DIR / "stage2_h5_files.json", jsonable(
        [dict(path=str(x["path"]), gsm=x["gsm"]) for x in out]
    ))
    return out


def load_one_h5(rec, log):
    path = rec["path"]
    log(f"   reading {path.name}")
    with h5py.File(path, "r") as f:
        if "matrix" not in f:
            raise StopStep("10x", f"{path.name} has no /matrix group")
        g = f["matrix"]
        shape = tuple(int(x) for x in g["shape"][()])
        data = g["data"][()]
        indices = g["indices"][()]
        indptr = g["indptr"][()]
        barcodes = _h5_decode(g["barcodes"][()])
        feat = g["features"]
        gene_id = _h5_decode(feat["id"][()])
        symbol = _h5_decode(feat["name"][()]) if "name" in feat else gene_id.copy()
        if "feature_type" in feat:
            ftype = _h5_decode(feat["feature_type"][()])
        else:
            ftype = np.array(["Gene Expression"] * len(gene_id), dtype=object)
        genome = _h5_decode(feat["genome"][()]) if "genome" in feat else np.array([""] * len(gene_id), dtype=object)
        n_genes, n_cells = shape
        # Cell Ranger: CSC genes × cells
        X = sparse.csc_matrix((data, indices, indptr), shape=(n_genes, n_cells)).T.tocsr()
        lib_ids = None
        if "library_ids" in f.attrs:
            lib_ids = _h5_decode(np.atleast_1d(f.attrs["library_ids"]))
        chemistry = f.attrs.get("chemistry_description", None)
        if isinstance(chemistry, (bytes, np.bytes_)):
            chemistry = chemistry.decode("utf-8")
        software = f.attrs.get("software_version", None)
        if isinstance(software, (bytes, np.bytes_)):
            software = software.decode("utf-8")
    obs = pd.DataFrame({
        "barcode": barcodes.astype(str),
        "gsm": rec.get("gsm"),
        "file": path.name,
    })
    var = pd.DataFrame({
        "gene_id": gene_id.astype(str),
        "symbol": symbol.astype(str),
        "feature_type": ftype.astype(str),
        "genome": genome.astype(str),
    })
    log(f"   {path.name}: {X.shape[0]} cells × {X.shape[1]} genes gsm={rec.get('gsm')} "
        f"chemistry={chemistry} software={software} library_ids={None if lib_ids is None else list(lib_ids)}")
    log(f"   columns actually read: matrix/barcodes, matrix/data, matrix/indices, matrix/indptr, "
        f"matrix/shape, matrix/features/id, matrix/features/name, matrix/features/feature_type, "
        f"matrix/features/genome, attrs chemistry_description/software_version/library_ids")
    return dict(X=X, obs=obs, var=var, gsm=rec.get("gsm"),
                chemistry=chemistry, software=software, library_ids=lib_ids)


def load_one_mtx(trip, log):
    mat, feat, bc = trip["matrix"], trip["features"], trip["barcodes"]
    log(f"   reading {mat.name}")
    opener = gzip.open if str(mat).endswith(".gz") else open
    with opener(mat, "rb") as f:
        X = mmread(f).T.tocsr()  # cells × genes
    genes = pd.read_csv(feat, sep="\t", header=None)
    barcodes = pd.read_csv(bc, sep="\t", header=None)
    if genes.shape[1] >= 2:
        gene_id = genes[0].astype(str).to_numpy()
        symbol = genes[1].astype(str).to_numpy()
        ftype = genes[2].astype(str).to_numpy() if genes.shape[1] >= 3 else np.array(["Gene Expression"] * len(genes))
    else:
        gene_id = genes[0].astype(str).to_numpy()
        symbol = gene_id.copy()
        ftype = np.array(["unknown"] * len(genes))
    obs = pd.DataFrame({
        "barcode": barcodes[0].astype(str).to_numpy(),
        "gsm": trip.get("gsm"),
        "file": mat.name,
    })
    var = pd.DataFrame({"gene_id": gene_id, "symbol": symbol, "feature_type": ftype})
    log(f"   {mat.name}: {X.shape[0]} cells × {X.shape[1]} genes gsm={trip.get('gsm')}")
    return dict(X=X, obs=obs, var=var, gsm=trip.get("gsm"))


def sendai_check(var_frames, log):
    """Report Sendai/vector features. Does not drop them silently."""
    hits = []
    all_names = []
    genomes = set()
    ftypes = set()
    columns_scanned = []
    n_ensg = n_non_ensg = n_feat = None
    for i, v in enumerate(var_frames):
        cols = [c for c in ("gene_id", "symbol", "feature_type", "genome") if c in v.columns]
        if i == 0:
            columns_scanned = cols
            n_feat = int(len(v))
            if "gene_id" in v.columns:
                gid = v["gene_id"].astype(str)
                n_ensg = int(gid.str.startswith("ENS").sum())
                n_non_ensg = int((~gid.str.startswith("ENS")).sum())
        if "genome" in v.columns:
            genomes.update(v["genome"].astype(str).unique().tolist())
        if "feature_type" in v.columns:
            ftypes.update(v["feature_type"].astype(str).unique().tolist())
        for col in cols:
            for name in v[col].astype(str):
                all_names.append(name)
                if SEV_RE.search(name):
                    hits.append(dict(source_var=i, column=col, name=name))
    uniq = sorted({h["name"] for h in hits})
    rec = dict(
        n_hits=len(uniq), names=uniq, n_features_scanned=len(set(all_names)),
        n_features_first_library=n_feat,
        n_ensg=n_ensg, n_non_ensg=n_non_ensg, n_libraries=len(var_frames),
        genomes=sorted(genomes), feature_types=sorted(ftypes),
        columns_actually_read=columns_scanned,
        pattern=SEV_RE.pattern,
        note="Sendai/vector features reported; not used in the frozen GTEx ruler. "
             "Pluripotency is scored with and without OSKM-family endogenous genes. "
             "Zero hits means the Cell Ranger matrix has no labelled Sendai/vector features; "
             "vector reads can still map onto endogenous OSKM gene models.",
    )
    dump_json(FIBRO_DIR / "stage2_sendai_check.json", jsonable(rec))
    log(f"[sendai] hits={len(uniq)} names={uniq[:40]} n_ensg={n_ensg} n_non_ensg={n_non_ensg} "
        f"genomes={sorted(genomes)} feature_types={sorted(ftypes)}")
    return rec


def attach_meta(libs, samples, log):
    by_gsm = samples.set_index(samples.gsm.astype(str))
    out = []
    for lib in libs:
        gsm = str(lib["gsm"]) if lib.get("gsm") else None
        if gsm is None or gsm not in by_gsm.index:
            raise StopStep("10x", f"library gsm={gsm!r} not in GEO sample table. Not guessing.")
        row = by_gsm.loc[gsm]
        if isinstance(row, pd.DataFrame):
            row = row.iloc[0]
        obs = lib["obs"].copy()
        obs["gsm"] = gsm
        obs["cell_line"] = str(row.cell_line).upper()
        obs["day"] = int(row.day)
        obs["age_years"] = int(row.age_years)
        obs["age_source"] = str(row.age_source)
        lib = dict(lib)
        lib["obs"] = obs
        lib["cell_line"] = str(row.cell_line).upper()
        lib["day"] = int(row.day)
        lib["age_years"] = int(row.age_years)
        out.append(lib)
        log(f"[meta] {gsm} {lib['cell_line']} day={lib['day']} age={lib['age_years']} "
            f"n_cells={len(obs)} age_source={row.age_source}")
    return out


def _var_ensembl_and_symbol(var: pd.DataFrame):
    gid = var["gene_id"].astype(str).to_numpy()
    sym = var["symbol"].astype(str).to_numpy()
    ens = np.array([strip_ensembl(g) for g in gid])
    is_ens = np.array([s.startswith("ENS") for s in ens])
    return ens, np.array([str(s).upper() for s in sym], dtype=object), is_ens


def align_counts_to_ruler(X, var, frozen, log):
    """Return cells × n_ruler counts on overlap genes; missing genes stay 0."""
    ruler_ens = np.asarray(frozen["ensembl"]).astype(str)
    ruler_sym = np.array([str(s).upper() for s in np.asarray(frozen["symbol"])], dtype=object)
    ens, sym_u, is_ens = _var_ensembl_and_symbol(var)
    pos_e, pos_s = {}, {}
    for i, e in enumerate(ens):
        if is_ens[i]:
            pos_e.setdefault(e, i)
        pos_s.setdefault(sym_u[i], i)
    n_r = len(ruler_ens)
    idx = np.full(n_r, -1, dtype=int)
    n_ens = n_sym = 0
    for j, e in enumerate(ruler_ens):
        i = pos_e.get(e)
        if i is not None:
            idx[j] = i
            n_ens += 1
            continue
        i = pos_s.get(ruler_sym[j])
        if i is not None:
            idx[j] = i
            n_sym += 1
    n_ov = int((idx >= 0).sum())
    log(f"[align] ruler_genes={n_r} overlap={n_ov} (ensembl={n_ens} symbol_fallback={n_sym}) missing={n_r-n_ov}")
    X = X.tocsr()
    take_r = np.flatnonzero(idx >= 0)
    take_x = idx[take_r]
    n_cells = int(X.shape[0])
    if take_r.size:
        sub = X[:, take_x].tocoo()
        Y = sparse.csr_matrix(
            (sub.data.astype(np.float64), (sub.row, take_r[sub.col])),
            shape=(n_cells, n_r),
        )
    else:
        Y = sparse.csr_matrix((n_cells, n_r), dtype=np.float64)
    rec = dict(n_ruler=n_r, n_overlap=n_ov, n_ensembl=n_ens, n_symbol_fallback=n_sym,
               n_missing=int(n_r - n_ov))
    return Y, rec, idx


def _row_sum(mat, mask=None):
    if mask is not None:
        mat = mat[mask]
    s = mat.sum(axis=0)
    return np.asarray(s, dtype=np.float64).ravel()


def pseudobulk_allcell(libs, frozen, log):
    rows, mats, aligns = [], [], []
    for lib in libs:
        Y, rec, idx = align_counts_to_ruler(lib["X"], lib["var"], frozen, log)
        pb = _row_sum(Y)
        rows.append(dict(
            cell_line=lib["cell_line"], day=lib["day"], age_years=lib["age_years"],
            gsm=lib["obs"]["gsm"].iloc[0], n_cells=int(Y.shape[0]),
            variant="all_cell", cluster=-1,
            age_source=lib["obs"]["age_source"].iloc[0],
        ))
        mats.append(pb)
        aligns.append(rec)
    C = np.vstack(mats)  # samples × genes
    # TMM wants genes × samples
    nf = tmm_norm_factors(C.T)
    logcpm = log2_cpm_edger(C.T, nf).T
    obs = pd.DataFrame(rows)
    dump_json(FIBRO_DIR / "stage2_allcell_align.json", jsonable(aligns))
    return obs, logcpm, C


def _log1p_cp10k(Y):
    lib = Y.sum(1, keepdims=True)
    lib = np.where(lib <= 0, 1.0, lib)
    return np.log1p(Y / lib * 1e4)


def cluster_and_pseudobulk(libs, frozen, log):
    """Unsupervised k-means within each donor × timepoint. Never mixed into (a)."""
    rng = np.random.default_rng(FIBRO_SEED)
    rows, mats, skipped = [], [], []
    for lib in libs:
        Y, rec, idx = align_counts_to_ruler(lib["X"], lib["var"], frozen, log)
        n = int(Y.shape[0])
        if n < CLUSTER_MIN_CELLS:
            skipped.append(dict(cell_line=lib["cell_line"], day=lib["day"], n_cells=n,
                                reason=f"n_cells<{CLUSTER_MIN_CELLS}"))
            log(f"[cluster] skip {lib['cell_line']} d{lib['day']}: n={n}")
            continue
        colsum = _row_sum(Y)
        keep = colsum > 0
        Ykeep = Y[:, keep]
        if sparse.issparse(Ykeep):
            Ykeep = Ykeep.astype(np.float32).toarray()
        Z = _log1p_cp10k(np.asarray(Ykeep, np.float32))
        npc = int(min(CLUSTER_NPC, n - 1, Z.shape[1]))
        pca = PCA(n_components=npc, random_state=int(FIBRO_SEED))
        P = pca.fit_transform(Z)
        km = KMeans(n_clusters=CLUSTER_K, random_state=int(FIBRO_SEED), n_init=10)
        lab = km.fit_predict(P)
        for k in range(CLUSTER_K):
            m = lab == k
            if int(m.sum()) == 0:
                continue
            pb = _row_sum(Y, m)
            rows.append(dict(
                cell_line=lib["cell_line"], day=lib["day"], age_years=lib["age_years"],
                gsm=lib["obs"]["gsm"].iloc[0], n_cells=int(m.sum()),
                variant="cluster", cluster=int(k),
                age_source=lib["obs"]["age_source"].iloc[0],
            ))
            mats.append(pb)
        log(f"[cluster] {lib['cell_line']} d{lib['day']} n={n} sizes={np.bincount(lab, minlength=CLUSTER_K).tolist()}")
    if not mats:
        n_g = int(np.asarray(frozen["w"]).shape[0])
        return pd.DataFrame(), np.zeros((0, n_g)), skipped, np.zeros((0, n_g))
    C = np.vstack(mats)
    nf = tmm_norm_factors(C.T)
    logcpm = log2_cpm_edger(C.T, nf).T
    dump_json(FIBRO_DIR / "stage2_cluster_skipped.json", jsonable(skipped))
    return pd.DataFrame(rows), logcpm, skipped, C


def score_age(logcpm, frozen, counts=None):
    """Missing genes: columns that are all-zero in the aligned count matrix → z=0."""
    mu = np.asarray(frozen["mu"], float)
    sd = np.asarray(frozen["sd"], float)
    sd = np.where(sd < 1e-12, 1.0, sd)
    w = np.asarray(frozen["w"], float)
    w, _ = unit(w)
    if counts is not None:
        missing = (np.asarray(counts) == 0).all(0)
    else:
        missing = (np.asarray(logcpm) == 0).all(0)
    Z = (np.asarray(logcpm, float) - mu) / sd
    Z[:, missing] = 0.0
    return Z @ w, Z, missing


def gene_index_by_symbol(frozen):
    sym = np.array([str(s).upper() for s in np.asarray(frozen["symbol"])], dtype=object)
    ens = np.asarray(frozen["ensembl"]).astype(str)
    pos = {}
    for i, s in enumerate(sym):
        pos.setdefault(s, i)
    return pos, ens, sym


def pluri_score(Z, frozen, gene_list, log, tag):
    pos, ens, sym = gene_index_by_symbol(frozen)
    idx, missing = [], []
    for g in gene_list:
        j = pos.get(str(g).upper())
        if j is None:
            missing.append(g)
        else:
            idx.append(j)
    log(f"[pluri {tag}] requested={list(gene_list)} mapped={len(idx)} missing={missing}")
    if not idx:
        return np.full(Z.shape[0], np.nan), dict(mapped=[], missing=list(gene_list), n_mapped=0)
    s = Z[:, np.array(idx, dtype=int)].mean(1)
    mapped = [g for g in gene_list if str(g).upper() in pos]
    return s, dict(mapped=mapped, missing=missing, n_mapped=len(idx), indices=idx)


def differentiation_score(Z, frozen, drop_oskm, log):
    pluri_list = list(PLURI_ENDOGENOUS)
    if drop_oskm:
        oskm_u = {x.upper() for x in OSKM_FAMILY}
        pluri_list = [g for g in pluri_list if g.upper() not in oskm_u]
    p, pmeta = pluri_score(Z, frozen, pluri_list, log, "endogenous" + ("_noOSKM" if drop_oskm else ""))
    f, fmeta = pluri_score(Z, frozen, list(FIBRO_IDENTITY), log, "fibro_identity")
    return p - f, dict(pluri=pmeta, fibro=fmeta, drop_oskm=drop_oskm, pluri_list=pluri_list)


def random_directions(w, Z, seed=FIBRO_SEED, n=N_RANDOM_DIR):
    rng = np.random.default_rng(int(seed))
    w = np.asarray(w, float).ravel()
    scores = np.zeros((Z.shape[0], n), dtype=np.float64)
    for i in range(n):
        wp = rng.permutation(w)
        wp, _ = unit(wp)
        scores[:, i] = Z @ wp
    return scores


def _pick(obs, logcpm_score, cell_line, day):
    m = (obs.cell_line.astype(str).str.upper() == cell_line) & (obs.day.astype(int) == int(day))
    if int(m.sum()) != 1:
        raise StopStep("score", f"{cell_line} day={day} n={int(m.sum())} (want 1) in all-cell table")
    return float(np.asarray(logcpm_score)[m.to_numpy()][0])


def first_change_day(days, values, kind):
    """kind='drop' -> first day < day0; kind='rise' -> first day > day0."""
    days = list(days)
    values = list(values)
    if 0 not in days:
        return None
    v0 = values[days.index(0)]
    later = [d for d in days if d > 0]
    for d in sorted(later):
        v = values[days.index(d)]
        if kind == "drop" and v < v0:
            return int(d)
        if kind == "rise" and v > v0:
            return int(d)
    return None


def read_outcome(age_days, age_s, pluri_days, pluri_s, p_end, sanity_ok):
    """Only the outcome that fired. Aged donor, all-cell (a)."""
    d0 = age_s[age_days.index(0)] if 0 in age_days else np.nan
    d10 = age_s[age_days.index(10)] if 10 in age_days else np.nan
    decline = float(d0 - d10) if np.isfinite(d0) and np.isfinite(d10) else np.nan
    mono = all(age_s[age_days.index(a)] >= age_s[age_days.index(b)]
               for a, b in zip(sorted(age_days), sorted(age_days)[1:])) if len(age_days) >= 2 else False
    drop_t = first_change_day(age_days, age_s, "drop")
    rise_t = first_change_day(pluri_days, pluri_s, "rise")
    beats = bool(np.isfinite(p_end) and p_end <= 0.05)
    declined = bool(np.isfinite(decline) and decline > 0)
    if not declined:
        key = "no_decline"
        text = (
            "age score does not decline → the ruler reads a static donor property, not a "
            "modifiable state. Step 3 is not supported by this data."
        )
    elif declined and not beats:
        key = "uninformative"
        text = (
            "age score declines but does not beat the random-direction null → uninformative; "
            "wholesale transcriptome change moves every direction, including this one. Do not report a bend."
        )
    elif declined and beats and rise_t is not None and drop_t is not None and rise_t > drop_t:
        key = "bend"
        text = (
            "age score declines, beats the random-direction null, and the pluripotency score rises "
            "later than the age drop → the curve bends; there is a measurable window where age "
            "moves before identity does."
        )
    elif declined and beats and rise_t is not None and drop_t is not None and rise_t == drop_t:
        key = "lockstep"
        text = (
            "age score declines and pluripotency rises in lockstep → no window; the field's core "
            "assumption is challenged. Report it as such, not as a failure."
        )
    elif declined and beats:
        key = "decline_other"
        text = (
            f"age score declines and beats the random-direction null; pluripotency rise_time={rise_t} "
            f"age drop_time={drop_t}. Not a listed bend/lockstep/uninformative/no-decline bullet. "
            "Not averaged into a headline."
        )
    else:
        key = "unclassified"
        text = "pattern did not match a pre-registered bullet. See tables. Do not average."
    if not sanity_ok:
        text = ("DAY-0 SANITY FAILED (aged donor not higher than young on the age ruler). "
                "Every downstream number is suspect. " + text)
    return dict(
        key=key, text=text, decline=decline, monotonic=bool(mono),
        drop_time=drop_t, rise_time=rise_t, beats_null=beats, declined=declined,
        p_endpoint=p_end, sanity_ok=bool(sanity_ok),
    )


def make_plot(obs, age, pluri, null_age_2d, path):
    fig, ax = plt.subplots(figsize=(6.2, 5.0))
    colors = {AGED_LINE: "#8c2d04", YOUNG_LINE: "#2171b5"}
    labels = {AGED_LINE: "GM00731 (96 y)", YOUNG_LINE: "GM23815 (22 y)"}
    null_age_2d = np.asarray(null_age_2d, float)
    for line in (YOUNG_LINE, AGED_LINE):
        m = obs.cell_line.astype(str).str.upper() == line
        if not m.any():
            continue
        idx = np.flatnonzero(m.to_numpy())
        order = np.argsort(obs.iloc[idx].day.astype(int).to_numpy())
        idx = idx[order]
        xs = np.asarray(age, float)[idx]
        ys = np.asarray(pluri, float)[idx]
        days = obs.iloc[idx].day.astype(int).to_numpy()
        lo = np.percentile(null_age_2d[idx], 5, axis=1)
        hi = np.percentile(null_age_2d[idx], 95, axis=1)
        ax.fill_betweenx(ys, lo, hi, color=colors[line], alpha=0.15, linewidth=0)
        ax.plot(xs, ys, "-o", color=colors[line], label=labels[line], ms=6)
        for x, y, day in zip(xs, ys, days):
            ax.annotate(f"d{int(day)}", (x, y), textcoords="offset points",
                        xytext=(5, 4), fontsize=7, color=colors[line])
    ax.set_xlabel("age score (frozen GTEx fibroblast ridge)")
    ax.set_ylabel("pluripotency − fibroblast identity")
    ax.legend(frameon=False, fontsize=8)
    ax.set_title("GSE297234 on frozen fibroblast age ruler")
    fig.tight_layout()
    fig.savefig(path, dpi=140)
    plt.close(fig)
    return path


def project(log):
    if not FROZEN_RULER.exists():
        raise StopStep("project", f"missing {FROZEN_RULER}")
    frozen = np.load(FROZEN_RULER, allow_pickle=True)
    series, samples_raw = fetch_gse_sample_table(log)
    samples, verify = verify_sample_table(series, samples_raw, log)
    tar, got = download_suppl(log)
    extract = extract_tar(tar, log)
    h5s = find_10x_h5(extract, log)
    trips = find_mtx_triplets(extract, log)
    if h5s:
        log(f"[10x] using Cell Ranger HDF5 n={len(h5s)} (not MTX; MTX triplets={len(trips)})")
        libs = [load_one_h5(h, log) for h in h5s]
        man = load_manifest()
        man["stage2_matrix_format"] = "10x_h5_filtered_feature_bc_matrix"
        save_manifest(man)
    elif trips:
        log(f"[10x] using MTX triplets n={len(trips)}")
        libs = [load_one_mtx(t, log) for t in trips]
        man = load_manifest()
        man["stage2_matrix_format"] = "10x_mtx"
        save_manifest(man)
    else:
        raise StopStep(
            "10x",
            "no 10x HDF5 (filtered_feature_bc_matrix.h5) and no matrix/features/barcodes "
            "triplets in GSE297234_RAW.tar. Not substituting a fitted Seurat object.",
            dict(extract=str(extract), downloads=got),
        )
    sendai = sendai_check([lib["var"] for lib in libs], log)
    libs = attach_meta(libs, samples, log)
    # require all 8 donor×day
    pairs = sorted({(l["cell_line"], l["day"]) for l in libs})
    expect = sorted((d["cell_line"], day) for d in GSE_EXPECTED_DONORS for day in GSE_EXPECTED_DAYS)
    # use the record's cell lines / days, not the expected if they differed
    rec_pairs = sorted({(str(r.cell_line).upper(), int(r.day)) for _, r in samples.iterrows()})
    if pairs != rec_pairs:
        raise StopStep("10x", f"loaded libraries {pairs} != GEO sample table {rec_pairs}")

    obs_a, logcpm_a, counts_a = pseudobulk_allcell(libs, frozen, log)
    age_a, Z_a, missing_a = score_age(logcpm_a, frozen, counts=counts_a)
    log(f"[score] all-cell missing-gene columns set to z=0: {int(missing_a.sum())}")
    pluri_with, meta_w = differentiation_score(Z_a, frozen, drop_oskm=False, log=log)
    pluri_wo, meta_wo = differentiation_score(Z_a, frozen, drop_oskm=True, log=log)
    # disagreement on aged donor days
    aged_m = obs_a.cell_line.astype(str).str.upper() == AGED_LINE
    days_a = obs_a.loc[aged_m, "day"].astype(int).tolist()
    w_s = np.asarray(pluri_with)[aged_m.to_numpy()]
    o_s = np.asarray(pluri_wo)[aged_m.to_numpy()]
    sign_w = np.sign(w_s[days_a.index(10)] - w_s[days_a.index(0)]) if (0 in days_a and 10 in days_a) else np.nan
    sign_o = np.sign(o_s[days_a.index(10)] - o_s[days_a.index(0)]) if (0 in days_a and 10 in days_a) else np.nan
    rise_w = first_change_day(days_a, list(w_s), "rise")
    rise_o = first_change_day(days_a, list(o_s), "rise")
    disagree = bool((sign_w != sign_o) or (rise_w != rise_o))
    pluri_a = pluri_wo if disagree else pluri_with
    pluri_primary = "without_OSKM" if disagree else "with_OSKM"
    log(f"[pluri] with/without OSKM disagree={disagree} → primary={pluri_primary}")

    null_a = random_directions(frozen["w"], Z_a, seed=FIBRO_SEED, n=N_RANDOM_DIR)

    def series_for(line, score):
        m = obs_a.cell_line.astype(str).str.upper() == line
        d = obs_a.loc[m, "day"].astype(int).to_numpy()
        s = np.asarray(score)[m.to_numpy()]
        order = np.argsort(d)
        return [int(x) for x in d[order]], [float(x) for x in s[order]]

    age_days_aged, age_s_aged = series_for(AGED_LINE, age_a)
    age_days_young, age_s_young = series_for(YOUNG_LINE, age_a)
    pluri_days_aged, pluri_s_aged = series_for(AGED_LINE, pluri_a)
    pluri_days_young, pluri_s_young = series_for(YOUNG_LINE, pluri_a)

    # day-0 sanity
    a0 = _pick(obs_a, age_a, AGED_LINE, 0)
    y0 = _pick(obs_a, age_a, YOUNG_LINE, 0)
    sanity_ok = bool(a0 > y0)
    log(f"[sanity] day0 aged={a0:+.4f} young={y0:+.4f} aged>young={sanity_ok}")

    # endpoint null on aged donor
    def _idx(line, day):
        m = (obs_a.cell_line.astype(str).str.upper() == line) & (obs_a.day.astype(int) == int(day))
        return int(np.flatnonzero(m.to_numpy())[0])

    i0, i10 = _idx(AGED_LINE, 0), _idx(AGED_LINE, 10)
    real_decl = float(age_a[i0] - age_a[i10])
    null_decl = null_a[i0] - null_a[i10]
    p_end = _perm_p(real_decl, null_decl, greater=True)
    # monotonicity: Spearman of score vs day on aged (a)
    real_mono = spearman_traj(np.asarray(age_days_aged, float), np.asarray(age_s_aged, float))
    # null mono: each random direction, spearman vs day on the same four points
    aged_idx = [int(x) for x in np.flatnonzero((obs_a.cell_line.astype(str).str.upper() == AGED_LINE).to_numpy())]
    days_ord = obs_a.iloc[aged_idx].day.astype(int).to_numpy()
    null_mono = []
    for j in range(null_a.shape[1]):
        null_mono.append(spearman_traj(days_ord.astype(float), null_a[aged_idx, j]))
    p_mono = _perm_p(real_mono, np.asarray(null_mono, float), greater=False)

    obs_a = obs_a.copy()
    obs_a["age_score"] = age_a
    obs_a["pluri_with_OSKM"] = pluri_with
    obs_a["pluri_without_OSKM"] = pluri_wo
    obs_a["pluri_primary"] = pluri_a
    obs_a.to_csv(FIBRO_DIR / "stage2_allcell_scores.csv", index=False)
    np.savez_compressed(FIBRO_DIR / "stage2_allcell_null.npz", null_age=null_a,
                        days=obs_a.day.to_numpy(), cell_line=obs_a.cell_line.astype(str).to_numpy())

    # cluster variant (b)
    obs_b, logcpm_b, skipped_b, counts_b = cluster_and_pseudobulk(libs, frozen, log)
    cluster_tbl = None
    if len(obs_b):
        age_b, Z_b, missing_b = score_age(logcpm_b, frozen, counts=counts_b)
        pluri_b, _ = differentiation_score(Z_b, frozen, drop_oskm=(pluri_primary == "without_OSKM"), log=log)
        obs_b = obs_b.copy()
        obs_b["age_score"] = age_b
        obs_b["pluri_primary"] = pluri_b
        obs_b.to_csv(FIBRO_DIR / "stage2_cluster_scores.csv", index=False)
        cluster_tbl = obs_b

    outcome = read_outcome(age_days_aged, age_s_aged, pluri_days_aged, pluri_s_aged, p_end, sanity_ok)
    fig_path = FIBRO_FIG / "plane.png"
    make_plot(obs_a, age_a, pluri_a, null_a, fig_path)

    young_decl = float(_pick(obs_a, age_a, YOUNG_LINE, 0) - _pick(obs_a, age_a, YOUNG_LINE, 10))
    summary = dict(
        sanity=dict(day0_aged=a0, day0_young=y0, aged_gt_young=sanity_ok),
        aged=dict(
            cell_line=AGED_LINE, days=age_days_aged, age_score=age_s_aged,
            pluri=pluri_s_aged, decline_0_to_10=real_decl, spearman_vs_day=real_mono,
            p_endpoint=float(p_end), p_monotonicity=float(p_mono),
            n_random=N_RANDOM_DIR, n_null_ge_real=int(np.sum(null_decl >= real_decl)),
        ),
        young=dict(
            cell_line=YOUNG_LINE, days=age_days_young, age_score=age_s_young,
            pluri=pluri_s_young, decline_0_to_10=young_decl,
            spearman_vs_day=spearman_traj(np.asarray(age_days_young, float), np.asarray(age_s_young, float)),
        ),
        pluri_primary=pluri_primary, pluri_disagree=disagree,
        sendai=sendai, verify=verify, outcome=outcome,
        figure=str(fig_path),
        matrix_format="10x_h5_filtered_feature_bc_matrix",
        n_allcell=int(len(obs_a)),
        n_cluster_rows=int(len(obs_b)) if cluster_tbl is not None else 0,
        cluster_skipped=skipped_b,
        missing_genes_z0=int(missing_a.sum()),
        n_ruler_genes=int(len(frozen["w"])),
        seed=FIBRO_SEED, n_perm_random=N_RANDOM_DIR,
        note="nothing fitted on GSE297234",
    )
    dump_json(FIBRO_DIR / "stage2_summary.json", jsonable(summary))
    dump_json(FIBRO_DIR / "stage2_pluri_meta.json", jsonable(dict(with_OSKM=meta_w, without_OSKM=meta_wo)))
    (FIBRO_DIR / "stage2_verdict.txt").write_text(outcome["text"] + "\n", encoding="utf-8")
    man = load_manifest()
    man["stage2_opened"] = True
    man["stage2_outcome"] = jsonable(outcome)
    man["status"] = "STAGE2_DONE"
    save_manifest(man)
    log(f"[verdict] {outcome['key']}: {outcome['text']}")
    return summary


def main():
    log = Logger(FIBRO_DIR / "stage2_report.txt")
    fibro_log_banner(log, "STAGE 2")
    try:
        gate = _gate_pass()
        log(f"[gate] Stage 1 passed. {gate.get('reason')}")
        write_prereg_stage2()
        log(f"[prereg] {PREREG_STAGE2_FLAG} exists BEFORE freeze/projection")
        progress_snapshot("freeze ruler on all GTEx fibroblast donors", stop="STAGE1_PASS Stage 2 prereg written")
        if FROZEN_RULER.exists() and (FIBRO_DIR / "freeze.json").exists():
            log(f"[freeze] already written {FROZEN_RULER}; not refitting")
        else:
            freeze_ruler(log)
        progress_snapshot("verify GSE297234 sample table; download; project", stop="STAGE1_PASS ruler frozen")
        project(log)
        progress_snapshot("write FINDINGS_FIBRO.md", stop="STAGE2_DONE")
    except StopStep as e:
        log(f"STOP [{e.step}] {e.message}")
        record_failure(e.step, e.message, e.details)
        dump_json(FIBRO_DIR / "stage2_STOP.json", dict(step=e.step, message=e.message, details=e.details))
        progress_snapshot(f"fix STOP [{e.step}]", stop=f"STOP {e.step}: {e.message}")
        raise
    except Exception as e:
        log(f"FAIL {type(e).__name__}: {e}")
        record_failure("stage2_fail", f"{type(e).__name__}: {e}")
        progress_snapshot(f"fix {type(e).__name__}", stop=f"FAIL {type(e).__name__}: {e}")
        raise
    finally:
        log.close()


if __name__ == "__main__":
    main()
