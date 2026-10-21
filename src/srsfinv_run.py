"""SRSFINV: data inventory of GEO series GSE325735 (Plesa et al., SRSF1 screen).

Inventory only. No statistic, no gate, no pre-registration, no distance, no
fitted axis. Reports what the GEO record and the sample metadata say.

Raw inputs, untracked under data/raw/srsfinv/ (.gitignore rule `data/raw/`),
fetched by this script when absent:
  - GSE325735_family.soft.gz  GEO family SOFT, every GSM's metadata
    https://ftp.ncbi.nlm.nih.gov/geo/series/GSE325nnn/GSE325735/soft/
  - sra_runinfo.csv           SRA run table for BioProject PRJNA1441554
    https://trace.ncbi.nlm.nih.gov/Traces/sra-db-be/runinfo?acc=PRJNA1441554
  - human_counts_head.gz.part first 262,144 bytes of
    GSE325735_human_raw_counts.txt.gz (HTTP range request; header and
    first rows only, the matrix itself is not downloaded)

Run from the repository root: python -I src/srsfinv_run.py
"""
import collections
import gzip
import hashlib
import json
import re
import zlib
from pathlib import Path

import numpy as np
import pandas as pd
import requests

ROOT = Path(__file__).resolve().parents[1]
RAW = ROOT / "data" / "raw" / "srsfinv"
OUT = ROOT / "results" / "srsfinv"
RULER = ROOT / "results" / "fibro" / "frozen_ruler_ridge_raw.npz"
SCOPING = ROOT / "results" / "scoping" / "gse325735_parse.json"

ACC = "GSE325735"
FTP = "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE325nnn/GSE325735"
SOFT_URL = FTP + "/soft/GSE325735_family.soft.gz"
SRA_URL = "https://trace.ncbi.nlm.nih.gov/Traces/sra-db-be/runinfo?acc=PRJNA1441554"
HUMAN_URL = FTP + "/suppl/GSE325735_human_raw_counts.txt.gz"
MOUSE_URL = FTP + "/suppl/GSE325735_mouse_raw_counts.txt.gz"
HEAD_BYTES = 262144
TASK_TITLE = ("Transcriptomic reprogramming screen identifies SRSF1 as a central "
              "regulator of aging and lifespan")

lines = []


def say(s=""):
    print(s)
    lines.append(s)


def md5(p):
    return hashlib.md5(p.read_bytes()).hexdigest()


def fetch(url, dest, byte_range=None):
    if dest.exists():
        return False
    headers = {"Range": "bytes=0-%d" % (byte_range - 1)} if byte_range else {}
    r = requests.get(url, headers=headers, timeout=120)
    r.raise_for_status()
    dest.write_bytes(r.content)
    return True


def head_size(url):
    r = requests.head(url, timeout=60, allow_redirects=True)
    r.raise_for_status()
    return int(r.headers["Content-Length"]), r.headers.get("Last-Modified", "")


def table(df, index=True):
    """Markdown table from a DataFrame."""
    cols = ([df.index.name or ""] if index else []) + [str(c) for c in df.columns]
    say("| " + " | ".join(cols) + " |")
    say("|" + "---|" * len(cols))
    for idx, row in df.iterrows():
        cells = ([str(idx)] if index else []) + [
            "" if (isinstance(v, float) and np.isnan(v)) else str(v) for v in row]
        say("| " + " | ".join(cells) + " |")


def parse_soft(path):
    series = collections.defaultdict(list)
    platforms = {}
    samples = []
    cur, kind = None, None
    with gzip.open(path, "rt", encoding="utf-8") as fh:
        for line in fh:
            line = line.rstrip("\n")
            if line.startswith("^"):
                kind, _, acc = line[1:].partition(" = ")
                if kind == "SAMPLE":
                    cur = collections.defaultdict(list)
                    samples.append(cur)
                elif kind == "PLATFORM":
                    cur = collections.defaultdict(list)
                    platforms[acc] = cur
                else:
                    cur = series
                continue
            if line.startswith("!") and " = " in line:
                k, _, v = line[1:].partition(" = ")
                cur[k].append(v)
    return series, platforms, samples


def sample_rows(samples):
    rows = []
    for s in samples:
        r = {
            "gsm": s["Sample_geo_accession"][0],
            "title": s["Sample_title"][0],
            "source_name": s["Sample_source_name_ch1"][0],
            "organism": s["Sample_organism_ch1"][0],
            "platform": s["Sample_platform_id"][0],
            "instrument": s["Sample_instrument_model"][0],
            "molecule": s["Sample_molecule_ch1"][0],
            "library_strategy": s["Sample_library_strategy"][0],
            "library_source": s["Sample_library_source"][0],
            "library_selection": s["Sample_library_selection"][0],
            "sample_type": s["Sample_type"][0],
            "n_char_lines": len(s["Sample_characteristics_ch1"]),
        }
        for c in s["Sample_characteristics_ch1"]:
            k, _, v = c.partition(": ")
            key = "ch_" + k.strip().replace(" ", "_")
            if key in r:
                raise SystemExit("duplicate characteristic %s in %s" % (key, r["gsm"]))
            r[key] = v.strip()
        desc = s["Sample_description"]
        r["description"] = " || ".join(desc)
        m = re.match(r"Library name: (\S+)$", desc[0])
        r["library_name"] = m.group(1) if m else ""
        r["description_rest"] = " || ".join(desc[1:])
        rows.append(r)
    df = pd.DataFrame(rows)
    rest = df["description_rest"]
    df["age_description"] = pd.to_numeric(rest.str.extract(r"Donor age (\d+)")[0])
    df["passage_description"] = pd.to_numeric(rest.str.extract(r"Passage (\d+)")[0])
    df["title_Mcode"] = pd.to_numeric(df["title"].str.extract(r"\bM(\d+)\b")[0])
    df["title_P"] = df["title"].str.extract(r", (P\+?\d+)\b")[0]
    df["title_Rep"] = df["title"].str.extract(r"\b(Rep\d+)\b")[0]
    df["description_flags"] = (rest.str.replace(r"^Donor age \d+(, Passage \d+)?", "", regex=True)
                               .str.lstrip(", ").str.strip())
    df["mouse_age_days"] = pd.to_numeric(rest.str.extract(r"Age at death (\d+) days")[0])
    df["mouse_rin"] = pd.to_numeric(rest.str.extract(r"RIN-([\d.]+)")[0])
    df["mouse_id"] = df["title"].str.extract(r"(Mouse#\d+)")[0]
    return df


def main():
    RAW.mkdir(parents=True, exist_ok=True)
    OUT.mkdir(parents=True, exist_ok=True)
    soft_p = RAW / "GSE325735_family.soft.gz"
    sra_p = RAW / "sra_runinfo.csv"
    head_p = RAW / "human_counts_head.gz.part"
    fetched = {p.name: fetch(u, p, b) for u, p, b in [
        (SOFT_URL, soft_p, None), (SRA_URL, sra_p, None), (HUMAN_URL, head_p, HEAD_BYTES)]}
    human_len, human_mod = head_size(HUMAN_URL)
    mouse_len, mouse_mod = head_size(MOUSE_URL)

    series, platforms, samples = parse_soft(soft_p)
    df = sample_rows(samples)
    df.to_csv(OUT / "samples.csv", index=False)
    hum = df[df.organism == "Homo sapiens"].copy()
    mou = df[df.organism == "Mus musculus"].copy()

    say("# SRSFINV - data inventory, GEO series GSE325735")
    say()
    say("Inventory only. Nothing here is a pre-registered result and nothing is")
    say("interpreted. No statistic, distance or fitted axis was computed.")
    say("Ambiguities are recorded as ambiguities, not resolved. Every count below")
    say("is a count of GEO sample records (GSMs) unless it says otherwise.")
    say()

    # ---- 0. provenance -----------------------------------------------------
    say("## 0. Provenance")
    say()
    say("Files read, all untracked under `data/raw/srsfinv/` (`.gitignore` rule `data/raw/`):")
    say()
    say("| file | source | bytes | md5 | fetched this run |")
    say("|---|---|---|---|---|")
    for p, u in [(soft_p, SOFT_URL), (sra_p, SRA_URL), (head_p, HUMAN_URL + " (bytes 0-%d)" % (HEAD_BYTES - 1))]:
        say("| `%s` | %s | %s | `%s` | %s |" % (p.name, u, format(p.stat().st_size, ","), md5(p), fetched[p.name]))
    say()
    say("Supplementary files on the GEO FTP (sizes from HTTP HEAD, not downloaded):")
    say()
    say("- `GSE325735_human_raw_counts.txt.gz` %s B, Last-Modified %s" % (format(human_len, ","), human_mod))
    say("- `GSE325735_mouse_raw_counts.txt.gz` %s B, Last-Modified %s" % (format(mouse_len, ","), mouse_mod))
    say()
    say("Download status, stated plainly: **no expression matrix was downloaded.**")
    say("The first %s bytes of the human counts file were read with an HTTP" % format(HEAD_BYTES, ","))
    say("range request (%.1f%% of the file) to see its header and identifier" % (100 * HEAD_BYTES / human_len))
    say("type. The mouse counts file was not opened. The SRA run table is")
    say("metadata, fetched for per-sample read depth, which GEO does not carry.")
    say()
    say("Series record fields:")
    say()
    for k in ["Series_title", "Series_status", "Series_submission_date", "Series_last_update_date",
              "Series_web_link", "Series_type", "Series_overall_design", "Series_relation"]:
        for v in series.get(k, []):
            say("- `%s`: %s" % (k, v))
    say("- `Series_sample_id`: %d GSMs; family SOFT holds %d `^SAMPLE` blocks"
        % (len(series["Series_sample_id"]), len(samples)))
    say("- `Series_supplementary_file`: " + "; ".join(series["Series_supplementary_file"]))
    say()
    gtitle = series["Series_title"][0]
    say("The GEO title today is \"%s\". The task names the series" % gtitle)
    say("\"%s\". They differ; the GEO record" % TASK_TITLE)
    say("was last updated %s, after `SCOPING_TISSUE.md` (2026-09-17)." % series["Series_last_update_date"][0])
    say("This file does not know which title was current on which date.")
    say()

    # ---- 0b. reproduce the stored scoping parse ------------------------------
    say("## 0b. Comparison with the stored scoping parse")
    say()
    say("`results/scoping/gse325735_parse.json` was built from GSM titles by")
    say("`results/scoping/_parse_325735.py`. The same title rules, re-applied here to")
    say("the titles in the family SOFT:")
    say()
    sc = json.loads(SCOPING.read_text(encoding="utf-8"))[ACC]
    titles = df["title"].tolist()
    t_h, t_m, t_o = [], [], []
    for t in titles:
        if re.search(r"\bMouse#", t, re.I) or t.startswith("Brain from") or t.startswith("Liver from"):
            t_m.append(t)
        elif re.search(r"\bM\d{2}\b", t) or "iPSC" in t or "NHDF" in t or "fibroblast" in t.lower() or "WT line" in t:
            t_h.append(t)
        else:
            t_o.append(t)
    mc = collections.Counter()
    tg = collections.Counter()
    for t in t_h:
        for m in re.findall(r"\bM(\d{2})\b", t):
            mc[str(int(m))] += 1
        if re.search(r"\bWT line\b", t):
            tg["WT"] += 1
            continue
        m = re.search(r"Dox-inducible ([A-Za-z0-9]+) transgene", t)
        if m:
            tg[m.group(1)] += 1
        else:
            m2 = re.search(r"overexpressing ([A-Za-z0-9]+)", t, re.I)
            tg[("OE:" + m2.group(1)) if m2 else ("UNPARSED:" + t[:80])] += 1
    comp = [("n_titles", sc["n_titles"], len(titles)),
            ("n_human_like", sc["n_human_like"], len(t_h)),
            ("n_mouse_like", sc["n_mouse_like"], len(t_m)),
            ("n_other", sc["n_other"], len(t_o)),
            ("n distinct title transgene tokens", sc["n_distinct_transgenes_incl_WT_unparsed"], len(tg))]
    for k, v in sc["donor_M_codes"].items():
        comp.append(("title M%s count" % k, v, mc.get(k, 0)))
    for k in sorted(set(sc["transgenes"]) | set(tg), key=lambda k: (-sc["transgenes"].get(k, 0), k)):
        comp.append(("title token %s" % k, sc["transgenes"].get(k, 0), tg.get(k, 0)))
    cdf = pd.DataFrame(comp, columns=["quantity", "stored", "this run"])
    cdf["match"] = cdf["stored"] == cdf["this run"]
    cdf.to_csv(OUT / "scoping_parse_comparison.csv", index=False)
    say("- quantities compared: %d; exact matches: %d; mismatches: %d"
        % (len(cdf), int(cdf["match"].sum()), int((~cdf["match"]).sum())))
    if (~cdf["match"]).any():
        table(cdf[~cdf["match"]], index=False)
    say("- full comparison: `results/srsfinv/scoping_parse_comparison.csv`")
    say()
    say("The scoping parse reads titles. Sections 2-4 below read the structured")
    say("`!Sample_characteristics_ch1` fields, which give different counts for")
    say("some labels (WT, NTg, OSKM, SRSF1); section 3c tabulates the difference.")
    say()

    # ---- 1. assay ------------------------------------------------------------
    say("## 1. Assay per sample")
    say()
    say("Fields that say what kind of assay each GSM is, all 839 GSMs:")
    say()
    for col in ["sample_type", "library_strategy", "library_source", "library_selection"]:
        vc = df[col].value_counts()
        say("- `%s`: %s" % (col, ", ".join("%s (%d)" % (k, v) for k, v in vc.items())))
    plat = (df.groupby(["organism", "platform", "instrument", "molecule"]).size()
            .rename("n GSMs").reset_index())
    plat["platform title"] = [platforms[p]["Platform_title"][0] for p in plat["platform"]]
    say()
    table(plat, index=False)
    say()
    tot = df[df.molecule == "total RNA"]
    say("`total RNA` GSMs (%d): batches %s; groups %s; ages %s." % (
        len(tot), ", ".join("%s %d" % kv for kv in tot.ch_batch.value_counts().items()),
        ", ".join("%s %d" % kv for kv in tot.organism.value_counts().items()),
        "/".join(str(int(a)) for a in sorted(tot.age_description.dropna().unique()))))
    say("The extract protocol below names a polyA mRNA workflow for the NEBNext libraries")
    say("and does not mention total RNA (ambiguity A12).")
    say()
    soft_text = gzip.open(soft_p, "rt", encoding="utf-8").read()
    sc_hits = sorted(set(m.lower() for m in re.findall(
        r"\b(?:single[- ]cell|scRNA|10x|Chromium|droplet|snRNA|cell barcode|UMIs?)\b", soft_text, re.I)))
    say("Whole-word search of the family SOFT text for single-cell terms (`single-cell`,")
    say("`single cell`, `scRNA`, `10x`, `Chromium`, `droplet`, `snRNA`, `cell barcode`, `UMI`):")
    say("%s." % ("no hits" if not sc_hits else "hits: " + ", ".join(sc_hits)))
    say("`Series_type` is \"%s\". Every GSM is `RNA-Seq`, `cDNA`, one" % series["Series_type"][0])
    say("library per GSM, and the processing text below describes STAR + per-gene")
    say("read counting of bulk libraries. **No single-cell assay is recorded.**")
    say("Counts below are therefore not split by assay; they are split by organism.")
    say()
    proc = samples[0]["Sample_data_processing"]
    same_proc = all(s["Sample_data_processing"] == proc for s in samples)
    say("`!Sample_data_processing`, identical in all 839 GSMs: %s" % same_proc)
    say()
    for p in proc:
        say("> %s" % p)
    say()
    ext = samples[0]["Sample_extract_protocol_ch1"]
    say("`!Sample_extract_protocol_ch1` (identical in all GSMs: %s):"
        % all(s["Sample_extract_protocol_ch1"] == ext for s in samples))
    say()
    for p in ext:
        say("> %s" % p)
    say()

    say("### 1b. Read depth")
    say()
    say("GEO carries no per-sample read-count field (`!Sample_data_row_count` is")
    say("%s in every GSM). Depth below is from the SRA run table." %
        ", ".join(sorted(set(s["Sample_data_row_count"][0] for s in samples))))
    sra = pd.read_csv(sra_p, dtype=str)
    sra = sra[sra["Run"].str.startswith("SRR", na=False)].copy()
    sra["spots"] = pd.to_numeric(sra["spots"])
    sra["bases"] = pd.to_numeric(sra["bases"])
    say("- SRA runs: %d; distinct `SampleName` (= GSM): %d; GSMs with no run: %d; runs whose `SampleName` is not a GSM of this series: %d"
        % (len(sra), sra["SampleName"].nunique(),
           len(set(df.gsm) - set(sra["SampleName"])), int((~sra["SampleName"].isin(df.gsm)).sum())))
    say("- `LibraryLayout`: %s" % ", ".join("%s (%d runs)" % kv for kv in sra["LibraryLayout"].value_counts().items()))
    say("- `Sex` column in the SRA table: %s" % ", ".join(
        "%r (%d)" % (k, v) for k, v in sra["Sex"].fillna("<empty>").value_counts().items()))
    per = sra.groupby("SampleName").agg(n_runs=("Run", "size"), spots=("spots", "sum"),
                                        bases=("bases", "sum"),
                                        layout=("LibraryLayout", lambda x: "/".join(sorted(set(x)))))
    per.index.name = "gsm"
    per = per.join(df.set_index("gsm")[["organism", "platform", "ch_batch"]])
    per.to_csv(OUT / "sra_depth_per_gsm.csv")
    say("- per-GSM depth table: `results/srsfinv/sra_depth_per_gsm.csv`. `spots` are")
    say("  read pairs for PAIRED runs, summed over a GSM's runs. These are depths")
    say("  **before** the downsampling the processing text describes.")
    say()
    q = (per.groupby(["organism", "platform"])["spots"]
         .describe(percentiles=[.05, .5, .95])[["count", "min", "5%", "50%", "95%", "max"]]
         .astype(int))
    q.index = [" ".join(i) for i in q.index]
    q.index.name = "organism platform"
    q = q.apply(lambda c: c.map("{:,}".format))
    table(q)
    say()
    lowd = per[per.spots < 20e6]
    say("GSMs with fewer than 20,000,000 spots in SRA: %d (%s)." % (
        len(lowd), ", ".join("%s %d" % kv for kv in lowd.organism.value_counts().items())))
    say()
    ds_libs = ["LIB058979", "LIB092645", "LIB083102"]
    say("GSMs in the three libraries the processing text exempts from 20 M downsampling")
    say("(`ch_batch` value): " + ", ".join("%s %d" % (b, int((df.ch_batch == b).sum())) for b in ds_libs)
        + "; all other GSMs %d." % int((~df.ch_batch.isin(ds_libs)).sum()))
    say()
    say("`!Sample_description` free-text flags beyond donor age and passage (human GSMs):")
    say()
    fl = hum["description_flags"].replace("", "<none>").value_counts()
    fl = fl.rename("n GSMs").to_frame()
    fl["batches"] = [", ".join(sorted(hum.loc[hum.description_flags.replace("", "<none>") == k, "ch_batch"].unique()))
                     for k in fl.index]
    table(fl)
    say()
    for k in ["<10M reads", "High GC"]:
        sub = hum[hum.description_flags == k]
        say("- `%s`: %s" % (k, ", ".join("%s (%s, %s spots)" % (r.gsm, r.title, format(int(per.loc[r.gsm, "spots"]), ","))
                                         for r in sub.itertuples())))
    redo = hum[hum.title.str.contains(r"\bredo\b")]
    say("- titles containing `redo`: %s" % ", ".join("%s (%s)" % (r.gsm, r.title) for r in redo.itertuples()))
    say()

    # ---- 2. donors -----------------------------------------------------------
    say("## 2. Donor table (human)")
    say()
    say("Where each field comes from:")
    say()
    say("- donor line: `!Sample_characteristics_ch1` `cell line:` (equals")
    say("  `!Sample_source_name_ch1` in %d of %d human GSMs)" %
        (int((hum["ch_cell_line"] == hum["source_name"]).sum()), len(hum)))
    say("- donor age: `!Sample_description` text `Donor age N`. Not in")
    say("  `!Sample_characteristics_ch1` (no `age` key exists there).")
    say("- title M-code: the `M<n>` token in `!Sample_title`")
    say("- passage: `!Sample_description` `Passage N`; titles also carry a `P<n>` token")
    say("- sex: **not stated** anywhere in the GEO record (no `sex`/`gender` key in")
    say("  characteristics; see 1b for the SRA `Sex` column). Not looked up elsewhere.")
    say()
    chkeys = collections.Counter(k for s in samples for c in s["Sample_characteristics_ch1"]
                                 for k in [c.partition(": ")[0]])
    say("Characteristic keys present, with GSM counts: " +
        ", ".join("`%s` %d" % kv for kv in chkeys.most_common()))
    sexhits = sorted(set(re.findall(r"\b(sex|gender|male|female)\b", soft_text, re.I)))
    say("Whole-SOFT search for `sex|gender|male|female`: %s." % (sexhits if sexhits else "no hits"))
    say()
    say("Age provenance per human GSM:")
    say()
    has_age = hum["age_description"].notna()
    has_m = hum["title_Mcode"].notna()
    agree = has_age & has_m & (hum["age_description"] == hum["title_Mcode"])
    say("- age stated in `!Sample_description`: %d of %d" % (int(has_age.sum()), len(hum)))
    say("- title M-code present: %d; of these M-code equals described age: %d; disagrees: %d"
        % (int(has_m.sum()), int(agree.sum()), int((has_m & has_age & ~agree).sum())))
    say("- age stated in description but no M-code in title: %d (lines %s)"
        % (int((has_age & ~has_m).sum()),
           ", ".join(sorted(hum.loc[has_age & ~has_m, "ch_cell_line"].unique()))))
    say("- age only inferable from title (no described age): %d" % int((~has_age & has_m).sum()))
    say("- no age from either source: %d" % int((~has_age & ~has_m).sum()))
    say()
    say("So no human GSM's age is inferred from the title alone within this record.")
    say("What the record does not say is the meaning of the `M` prefix; it equals")
    say("the described age wherever both exist (see ambiguity A2).")
    say()
    g = hum.groupby("ch_cell_line")
    don = pd.DataFrame({
        "cell type": g["ch_cell_type"].agg(lambda x: "/".join(sorted(set(x)))),
        "n Fibroblast GSMs": g["ch_cell_type"].agg(lambda x: int((x == "Fibroblast").sum())),
        "n iPSC GSMs": g["ch_cell_type"].agg(lambda x: int((x == "iPSC").sum())),
        "described age(s)": g["age_description"].agg(lambda x: "/".join(str(int(v)) for v in sorted(set(x.dropna())))),
        "title M-code(s)": g["title_Mcode"].agg(lambda x: "/".join(str(int(v)) for v in sorted(set(x.dropna()))) or "-"),
        "described passage(s)": g["passage_description"].agg(lambda x: "/".join(str(int(v)) for v in sorted(set(x.dropna()))) or "-"),
        "title P token(s)": g["title_P"].agg(lambda x: "/".join(sorted(set(x.dropna()))) or "-"),
        "sex": "not stated",
        "n GSMs": g.size(),
        "n batches": g["ch_batch"].nunique(),
    })
    don.index.name = "cell line"
    don = don.sort_values(["described age(s)", "n GSMs"])
    don.to_csv(OUT / "donors.csv")
    n_fib_lines = int((don["n Fibroblast GSMs"] > 0).sum())
    n_ips_lines = int((don["n iPSC GSMs"] > 0).sum())
    say("Distinct `cell line` values among human GSMs: %d. With at least one `Fibroblast`"
        % len(don))
    say("GSM: %d. With at least one `iPSC` GSM: %d (%s; %s has both)." % (
        n_fib_lines, n_ips_lines, ", ".join(don.index[don["n iPSC GSMs"] > 0]),
        ", ".join(don.index[(don["n iPSC GSMs"] > 0) & (don["n Fibroblast GSMs"] > 0)]) or "none"))
    say("Distinct described ages among human fibroblast GSMs: %d. (`SCOPING_TISSUE.md`"
        % hum.loc[hum.ch_cell_type == "Fibroblast", "age_description"].nunique())
    say("wrote \"16 M-coded NHDF lines\"; 16 is the number of distinct M-codes, i.e. ages,")
    say("not of lines.)")
    say()
    table(don)
    say()
    ip = hum[hum.ch_cell_type == "iPSC"]
    say("The %d iPSC GSMs, verbatim titles: %s." % (len(ip), "; ".join(
        "%s `%s` (line %s, described age %d)" % (r.gsm, r.title, r.ch_cell_line, r.age_description)
        for r in ip.itertuples())))
    iso = hum[hum.title.str.contains("isogenic")]
    say("Titles containing `isogenic` (%d): %s." % (len(iso), "; ".join(
        "%s `%s` (library %s)" % (r.gsm, r.title, r.library_name) for r in iso.itertuples())))
    say()
    lines_per_age = hum.groupby("age_description")["ch_cell_line"].agg(lambda x: sorted(set(x)))
    multi = {int(a): l for a, l in lines_per_age.items() if len(l) > 1}
    say("Ages carried by more than one `cell line` value: %s." %
        ("; ".join("%d: %s" % (a, ", ".join(l)) for a, l in multi.items()) if multi else "none"))
    age_per_line = hum.groupby("ch_cell_line")["age_description"].nunique()
    say("Cell lines with more than one described age: %s." %
        (", ".join(age_per_line[age_per_line > 1].index) or "none"))
    say("Whether two `cell line` values are the same human donor is not stated; the")
    say("table counts `cell line` values, not verified donors.")
    say()

    # ---- 3. perturbations ----------------------------------------------------
    say("## 3. Perturbation table")
    say()
    say("Source: `!Sample_characteristics_ch1` `genotype:`. Every value is listed as")
    say("written. Which ones are controls is section 4.")
    say()
    df["group"] = np.where(df.organism == "Mus musculus", "mouse",
                           np.where(df.get("ch_cell_type") == "iPSC", "human iPSC", "human fibroblast"))
    hum["group"] = df.loc[hum.index, "group"]
    gt = pd.crosstab(df["ch_genotype"], df["group"])
    gt["total"] = gt.sum(axis=1)
    hf = hum[hum.group == "human fibroblast"]
    gt["n fibroblast donor ages"] = hf.groupby("ch_genotype")["age_description"].nunique().reindex(gt.index).fillna(0).astype(int)
    gt["n fibroblast cell lines"] = hf.groupby("ch_genotype")["ch_cell_line"].nunique().reindex(gt.index).fillna(0).astype(int)
    gt = gt.sort_values(["total"], ascending=False)
    gt.index.name = "genotype"
    gt.to_csv(OUT / "genotype_counts.csv")
    say("- distinct `genotype` values, all GSMs: %d" % len(gt))
    say("- distinct `genotype` values in human fibroblasts: %d; in human iPSC: %d; in mouse: %d"
        % (hf.ch_genotype.nunique(), hum[hum.group == "human iPSC"].ch_genotype.nunique(), mou.ch_genotype.nunique()))
    nonctrl = sorted(set(hf.ch_genotype) - {"WT", "NTg", "mTagBFP2", "GFP"})
    say("- human-fibroblast `genotype` values other than `WT`, `NTg`, `mTagBFP2`, `GFP`: %d" % len(nonctrl))
    say()
    table(gt)
    say()
    say("### 3b. Genotype x described donor age, human fibroblasts")
    say()
    xt = pd.crosstab(hf["ch_genotype"], hf["age_description"].astype(int))
    xt.columns = ["age %d" % c for c in xt.columns]
    xt = xt.loc[gt.index.intersection(xt.index, sort=False)]
    xt.to_csv(OUT / "genotype_by_age_fibroblast.csv")
    table(xt)
    say()
    multi_age = gt[gt["n fibroblast donor ages"] > 1].index.tolist()
    say("**Genotypes applied to more than one donor age in human fibroblasts (%d):**" % len(multi_age))
    say()
    for k in multi_age:
        ages = sorted(int(a) for a in hf.loc[hf.ch_genotype == k, "age_description"].unique())
        say("- `%s`: ages %s" % (k, ", ".join(map(str, ages))))
    say()
    tg_multi = [k for k in multi_age if k not in ("WT", "NTg", "mTagBFP2", "GFP")]
    say("Of these, %d are not `WT`, `NTg` or a reporter (`mTagBFP2`, `GFP`)." % len(tg_multi))
    tgf = hf[~hf.ch_genotype.isin(["WT", "NTg", "mTagBFP2", "GFP"])]
    say("Human fibroblast GSMs carrying a genotype other than `WT`/`NTg`/reporter: %d," % len(tgf))
    say("from cell lines %s, at described ages %s. **No such GSM comes from a donor"
        % (", ".join(sorted(tgf.ch_cell_line.unique())),
           ", ".join(str(int(a)) for a in sorted(tgf.age_description.unique()))))
    say("younger than %d.**" % int(tgf.age_description.min()))
    say()
    one_age = gt[(gt["n fibroblast donor ages"] == 1)].index.tolist()
    oa = collections.Counter(int(hf.loc[hf.ch_genotype == k, "age_description"].iloc[0]) for k in one_age)
    say("Genotypes seen at exactly one fibroblast donor age: %d, by that age: %s."
        % (len(one_age), ", ".join("%d (%d)" % kv for kv in sorted(oa.items()))))
    say("Genotypes with no human fibroblast GSM: %s."
        % (", ".join(gt[gt["human fibroblast"] == 0].index) if "human fibroblast" in gt else "-"))
    say()
    say("### 3c. Title token versus `genotype` field")
    say()
    say("Title token = the scoping rule (`WT line` -> WT, else `Dox-inducible X transgene` -> X).")
    say("Human GSMs where it differs from `genotype`:")
    say()
    def title_token(t):
        if re.search(r"\bWT line\b", t) or re.search(r"\bWT iPSC line\b", t):
            return "WT"
        m = re.search(r"Dox-inducible (.+?) transgene", t)
        return m.group(1) if m else "<none>"
    hum["title_token"] = hum["title"].map(title_token)
    diff = hum[hum["title_token"] != hum["ch_genotype"]]
    dt = diff.groupby(["title_token", "ch_genotype"]).size().rename("n GSMs").reset_index()
    if len(dt):
        table(dt, index=False)
    else:
        say("none")
    say()
    say("Examples (first GSM of each pair):")
    say()
    for (a, b), sub in diff.groupby(["title_token", "ch_genotype"]):
        r = sub.iloc[0]
        say("- %s title `%s`; genotype `%s`; treatment `%s`" % (r.gsm, r.title, r.ch_genotype, r.ch_treatment))
    say()

    # ---- 4. controls ---------------------------------------------------------
    say("## 4. Controls")
    say()
    say("The record never uses the words control, untreated or vector. The")
    say("candidate control conditions below are the stated field values; which of")
    say("them counts as \"untreated\" is not stated (ambiguity A4).")
    say()
    say("Values of `treatment:` among human GSMs: " + ", ".join(
        "`%s` %d" % kv for kv in hum["ch_treatment"].value_counts().items()))
    say()
    cat = []
    for _, r in hum.iterrows():
        if r.ch_genotype == "WT" and r.ch_treatment == "NT" and "NTg" in r.title:
            c = "C1b genotype WT, treatment NT, title says NTg transgene"
        elif r.ch_genotype == "WT" and r.ch_treatment == "NT":
            c = "C1a genotype WT, treatment NT, title says WT/plain line"
        elif r.ch_genotype == "WT" and r.ch_treatment.startswith("Rapa"):
            c = "C2 genotype WT, rapamycin"
        elif r.ch_genotype == "WT":
            c = "C3 genotype WT, dox (title says NTg transgene)" if "NTg" in r.title else "C3' genotype WT, dox"
        elif r.ch_genotype == "NTg":
            c = "C4 genotype NTg, dox" if r.ch_treatment.startswith("Dox") else "C4' genotype NTg, " + r.ch_treatment
        elif r.ch_genotype in ("mTagBFP2", "GFP") and r.ch_treatment == "NT":
            c = "C5' reporter, treatment NT"
        elif r.ch_genotype in ("mTagBFP2", "GFP"):
            c = "C5 reporter (%s), dox" % r.ch_genotype
        elif r.ch_treatment == "NT":
            c = "C6 other transgene, treatment NT (not induced)"
        else:
            c = "-"
        cat.append(c)
    hum["control_candidate"] = cat
    cc = hum[hum.control_candidate != "-"]
    ct = pd.crosstab([cc.control_candidate, cc.group], cc.age_description.astype(int))
    ct.columns = ["age %d" % c for c in ct.columns]
    ct["total"] = ct.sum(axis=1)
    ct.index = [" / ".join(i) for i in ct.index]
    ct.index.name = "candidate / group"
    ct.to_csv(OUT / "control_candidates_by_age.csv")
    table(ct)
    say()
    say("Labels are this file's shorthand for field combinations, not the authors'.")
    say("Primed labels (C3', C4', C5') would be rows with zero GSMs if absent above.")
    for lab in sorted(cc.control_candidate.unique()):
        sub = cc[cc.control_candidate == lab]
        say("- %s: treatments %s; batches %s" % (lab, ", ".join("`%s` %d" % kv for kv in sub.ch_treatment.value_counts().items()),
                                                 ", ".join("%s %d" % kv for kv in sub.ch_batch.value_counts().items())))
    say()
    say("Human GSMs with `genotype: WT` and a dox treatment: %d; of these, titles say"
        % int((cc.control_candidate.str.startswith("C3")).sum()))
    say("\"Dox-inducible NTg transgene\": %d (ambiguity A3)."
        % int((cc.control_candidate == "C3 genotype WT, dox (title says NTg transgene)").sum()))
    nt = hum[hum.ch_treatment == "NT"]
    nt0 = nt[nt.title.str.contains(r"Dox_0ng/mL", regex=False)]
    say("Of the %d human `treatment: NT` GSMs, %d have `Dox_0ng/mL` in the title (all C1b);"
        % (len(nt), len(nt0)))
    say("the rest have no dox token in the title. That is the only place the record")
    say("hints at what `NT` means.")
    say()
    mol = hum.groupby("control_candidate")["molecule"].agg(lambda x: ", ".join(
        "%s %d" % kv for kv in x.value_counts().items()))
    say("`molecule` by candidate: " + "; ".join("%s: %s" % kv for kv in mol.items() if kv[0] != "-")
        + "; transgene GSMs: " + ", ".join("%s %d" % kv for kv in tgf.molecule.value_counts().items()) + ".")
    say()
    say("### 4b. Every human GSM from a donor under 40")
    say()
    u40 = hum[hum.age_description < 40].sort_values(["age_description", "ch_cell_line", "gsm"])
    u40[["gsm", "library_name", "title", "ch_cell_line", "ch_cell_type", "age_description",
         "passage_description", "ch_genotype", "ch_treatment", "ch_batch", "description_flags",
         "control_candidate"]].to_csv(OUT / "under40_samples.csv", index=False)
    say("Human GSMs with described age < 40: %d (fibroblast %d, iPSC %d)."
        % (len(u40), int((u40.group == "human fibroblast").sum()), int((u40.group == "human iPSC").sum())))
    say()
    show = u40[["gsm", "library_name", "title", "ch_cell_line", "ch_cell_type", "age_description",
                "passage_description", "ch_genotype", "ch_treatment", "ch_batch", "control_candidate"]].copy()
    show["age_description"] = show["age_description"].astype("Int64")
    show["passage_description"] = show["passage_description"].astype("Int64").astype(str).replace("<NA>", "-")
    table(show.rename(columns=lambda c: c.replace("ch_", "")), index=False)
    say()
    u40f = u40[u40.group == "human fibroblast"]
    say("**Untreated (`genotype WT`, `treatment NT`) human fibroblast GSMs from donors under 40: %d.**"
        % int(((u40f.ch_genotype == "WT") & (u40f.ch_treatment == "NT")).sum()))
    say("Their batches: %s." % (", ".join("%s %d" % kv for kv in
        u40f[(u40f.ch_genotype == "WT") & (u40f.ch_treatment == "NT")].ch_batch.value_counts().items()) or "-"))
    say()
    say("### 4c. Batch co-measurement of the under-40 fibroblast GSMs")
    say()
    say("For each batch holding an under-40 human fibroblast GSM, all human GSMs in")
    say("that batch by group and age (descriptive; the record does not say batch is a")
    say("sequencing run, a plate or a culture round):")
    say()
    bs = sorted(u40f.ch_batch.unique())
    bt = pd.crosstab([hum.ch_batch, hum.group], hum.age_description.astype(int)).loc[bs]
    bt = bt.loc[:, (bt.sum() > 0)]
    bt.columns = ["age %d" % c for c in bt.columns]
    bt.index = [" / ".join(i) for i in bt.index]
    bt.index.name = "batch / group"
    table(bt)
    say()
    bt_all = pd.crosstab(df.ch_batch, df.group)
    bt_all["described ages"] = hum.groupby("ch_batch")["age_description"].agg(
        lambda x: "/".join(str(int(v)) for v in sorted(set(x.dropna())))).reindex(bt_all.index).fillna("-")
    bt_all["genotypes"] = df.groupby("ch_batch")["ch_genotype"].agg(
        lambda x: ", ".join(sorted(set(x))) if x.nunique() <= 5 else "%d values" % x.nunique())
    bt_all["molecule"] = df.groupby("ch_batch")["molecule"].agg(lambda x: "/".join(sorted(set(x))))
    bt_all["platform"] = df.groupby("ch_batch")["platform"].agg(lambda x: "/".join(sorted(set(x))))
    bt_all["fibroblast transgene GSMs"] = tgf.groupby("ch_batch").size().reindex(bt_all.index).fillna(0).astype(int)
    wtnt = hf[(hf.ch_genotype == "WT") & (hf.ch_treatment == "NT")]
    bt_all["fibroblast WT+NT GSMs"] = wtnt.groupby("ch_batch").size().reindex(bt_all.index).fillna(0).astype(int)
    bt_all["fibroblast WT+NT ages"] = wtnt.groupby("ch_batch")["age_description"].agg(
        lambda x: "/".join(str(int(v)) for v in sorted(set(x)))).reindex(bt_all.index).fillna("-")
    bt_all.to_csv(OUT / "batches.csv")
    say("All batches (\"transgene\" = genotype other than `WT`/`NTg`/reporter):")
    say()
    table(bt_all)
    say()
    shared = sorted(set(u40f.ch_batch) & set(tgf.ch_batch))
    say("Batches holding both an under-40 human fibroblast GSM and a fibroblast")
    say("transgene GSM: %s." % (", ".join(shared) if shared else "**none**"))
    say()

    # ---- 5. timepoints, doses ------------------------------------------------
    say("## 5. Timepoints, doses, induction conditions")
    say()
    say("Source: `!Sample_characteristics_ch1` `treatment:`. Treatment protocol text")
    say("(identical in all GSMs: %s):" % all(
        s["Sample_treatment_protocol_ch1"] == samples[0]["Sample_treatment_protocol_ch1"] for s in samples))
    say()
    for p in samples[0]["Sample_treatment_protocol_ch1"]:
        say("> %s" % p)
    say()
    say("The protocol text says rapamycin at 100 nM; the field values say")
    say("`Rapa_100ng/mL`. Recorded, not resolved (ambiguity A6).")
    say()
    tr = pd.crosstab(df["ch_treatment"], df["group"])
    tr["total"] = tr.sum(axis=1)
    tr["n genotypes"] = df.groupby("ch_treatment")["ch_genotype"].nunique()
    tr["genotypes (if <= 6)"] = df.groupby("ch_treatment")["ch_genotype"].agg(
        lambda x: ", ".join(sorted(set(x))) if x.nunique() <= 6 else "")
    tr["fibroblast ages"] = hf.groupby("ch_treatment")["age_description"].agg(
        lambda x: "/".join(str(int(v)) for v in sorted(set(x)))).reindex(tr.index).fillna("-")
    tr = tr.sort_values("total", ascending=False)
    tr.index.name = "treatment"
    tr.to_csv(OUT / "treatments.csv")
    table(tr)
    say()
    say("Reading of the codes, as written and not interpreted further: `Dox_<c>` is")
    say("doxycycline concentration; `<n>d-ON/<m>d-OFF` reads as n days on dox then m")
    say("days off; `3days` and `3d-ON` are both written and may or may not be the same")
    say("schedule (ambiguity A5); `NT` is not defined in the record.")
    say()
    say("Title tokens: `P<n>` / `P+<n>` (%d human GSMs) and `Rep<n>` (%d human GSMs)."
        % (int(hum.title_P.notna().sum()), int(hum.title_Rep.notna().sum())))
    pp = hum[hum.title_P.notna() & hum.passage_description.notna()]
    same_p = (pp.title_P.str.replace("P", "", regex=False) == pp.passage_description.astype(int).astype(str))
    say("Where both a title `P<n>` and a described `Passage N` exist (%d GSMs), they are equal in %d."
        % (len(pp), int(same_p.sum())))
    for r in pp[~same_p].itertuples():
        say("- unequal: %s title `%s`; description `%s`" % (r.gsm, r.title, r.description_rest))
    say("Described passage present for %d of %d human GSMs." % (int(hum.passage_description.notna().sum()), len(hum)))
    say()

    # ---- 6. human vs mouse ---------------------------------------------------
    say("## 6. Human versus mouse")
    say()
    say("| organism | n GSMs | platforms | sample kind |")
    say("|---|---|---|---|")
    say("| Homo sapiens | %d | %s | `cell line` + `cell type` (%s) |" % (
        len(hum), ", ".join(sorted(hum.platform.unique())),
        ", ".join("%s %d" % kv for kv in hum.ch_cell_type.value_counts().items())))
    say("| Mus musculus | %d | %s | `tissue` (%s) |" % (
        len(mou), ", ".join(sorted(mou.platform.unique())),
        ", ".join("%s %d" % kv for kv in mou.ch_tissue.value_counts().items())))
    say()
    say("Human and mouse are separated by `!Sample_organism_ch1`; no GSM lacks it.")
    say("The human counts file is a separate supplementary file from the mouse one;")
    say("its header is compared with the human GSMs in section 7.")
    say()
    say("Mouse GSMs (all AAV8, dox %s; `treatment` %s):" % (
        "200 mg/kg", ", ".join("`%s`" % v for v in mou.ch_treatment.unique())))
    say()
    mt = mou.groupby(["ch_genotype", "mouse_id"]).agg(
        tissues=("ch_tissue", lambda x: "/".join(sorted(x))),
        age_at_death_days=("mouse_age_days", lambda x: "/".join(str(int(v)) for v in sorted(set(x)))),
        RIN=("mouse_rin", lambda x: "/".join(str(v) for v in x)),
        batch=("ch_batch", lambda x: "/".join(sorted(set(x)))))
    mt.index = [" / ".join(i) for i in mt.index]
    mt.index.name = "genotype / mouse"
    table(mt)
    say()
    say("Mouse: %d animals (%s), 4 tissues each. The mouse ages are \"Age at death\"," % (
        mou.mouse_id.nunique(), ", ".join("%s %d" % kv for kv in
                                         mou.groupby("ch_genotype").mouse_id.nunique().items())))
    say("in days, from `!Sample_description`. No young mouse group is recorded.")
    say()

    # ---- 7. gene identifiers -------------------------------------------------
    say("## 7. Gene identifiers")
    say()
    raw = head_p.read_bytes()
    d = zlib.decompressobj(16 + zlib.MAX_WBITS)
    txt = d.decompress(raw).decode("utf-8", errors="replace")
    crlf = "\r\n" in txt
    rows = txt.replace("\r\n", "\n").split("\n")
    header = rows[0].split("\t")
    complete = [r.split("\t") for r in rows[1:-1]]
    complete = [r for r in complete if len(r) == len(header)]
    ids = [r[0] for r in complete]
    say("From the %s-byte range read of `GSE325735_human_raw_counts.txt.gz`:" % format(HEAD_BYTES, ","))
    say()
    say("- first header field: `%s`; sample columns: %d; line endings CRLF: %s"
        % (header[0], len(header) - 1, crlf))
    say("- complete data rows decoded: %d (the file's total row count is not known" % len(ids))
    say("  without reading it all)")
    pat = collections.Counter(re.sub(r"\d", "#", re.sub(r"\.\d+$", ".v", i)) for i in ids)
    say("- identifier pattern of those rows: %s" % ", ".join("`%s` %d" % kv for kv in pat.items()))
    say("- first 5 identifiers: %s" % ", ".join("`%s`" % i for i in ids[:5]))
    say("- identifier type: versioned Ensembl gene IDs (`ENSG...`.version); no")
    say("  symbol column in the rows read")
    hl = set(hum.library_name)
    hc = set(header[1:])
    say("- header columns vs human GSM library names (`Library name:` in")
    say("  `!Sample_description`): equal sets: %s; in header only: %d; in GEO only: %d"
        % (hc == hl, len(hc - hl), len(hl - hc)))
    if hc != hl:
        say("  - header only: %s" % sorted(hc - hl)[:20])
        say("  - GEO only: %s" % sorted(hl - hc)[:20])
    say()
    z = np.load(RULER, allow_pickle=True)
    rg_v, rg_u = set(z["gene_id"].tolist()), set(z["ensembl"].tolist())
    n_v = sum(i in rg_v for i in ids)
    n_u = sum(i.split(".")[0] in rg_u for i in ids)
    say("**Frozen-ruler overlap: not counted.** The ruler (`results/fibro/frozen_ruler_ridge_raw.npz`,")
    say("n = %d genes, keys `gene_id` versioned and `ensembl` unversioned) can only be" % len(z["gene_id"]))
    say("matched against the full identifier column, which needs the whole")
    say("%s-byte matrix. That is a matrix download, so it was not done." % format(human_len, ","))
    say()
    say("Key-format check only, not an overlap estimate. Among the %d identifiers read," % len(ids))
    say("%d match a ruler `ensembl` once the version suffix is dropped; %d of those also" % (n_u, n_v))
    say("match the ruler's versioned `gene_id` exactly. So the Ensembl versions in this")
    say("file and in the ruler's annotation differ for some genes, and any later match")
    say("would have to be made on unversioned IDs (a choice for that task, not made")
    say("here). These rows are the first in file order, not a sample; the numbers say")
    say("nothing about how many of the 23,485 ruler genes are present.")
    say()
    say("The mouse counts file (%s B) was not opened." % format(mouse_len, ","))
    say()

    # ---- 8. ambiguities ------------------------------------------------------
    say("## 8. Ambiguities, recorded and not resolved")
    say()
    n_ntg_wt = int(((hum.title_token == "NTg") & (hum.ch_genotype == "WT")).sum())
    n_ntg = int((hum.ch_genotype == "NTg").sum())
    amb = [
        ("A1", "Series title. GEO now reads \"%s\"; the task names a different title. "
               "The record was updated %s." % (gtitle, series["Series_last_update_date"][0])),
        ("A2", "Donor identity and age. Age is free text in `!Sample_description`, not a "
               "structured characteristic. The title `M<n>` equals it wherever both exist, but "
               "the record never says what `M` means. Sex is not stated. Whether each `cell line` "
               "value is a distinct donor is not stated; ages 36 and 55 each belong to two lines. "
               "Coriell catalog entries for the AG IDs were not consulted."),
        ("A3", "`NTg`. %d human GSMs have `genotype: NTg`; %d more have `genotype: WT` but a "
               "title reading \"line with Dox-inducible NTg transgene\" (%d with dox, %d with "
               "`NT`/`Dox_0ng/mL`). The record does not define NTg, nor say whether these lines "
               "carry a construct. So `genotype: WT` covers both plain \"WT line\" GSMs and "
               "NTg-titled GSMs." % (n_ntg, n_ntg_wt,
                                     int(((hum.title_token == "NTg") & (hum.ch_genotype == "WT") & (hum.ch_treatment != "NT")).sum()),
                                     int(((hum.title_token == "NTg") & (hum.ch_genotype == "WT") & (hum.ch_treatment == "NT")).sum()))),
        ("A4", "What is \"untreated\". The record does not use the word, and `NT` is undefined "
               "(NTg-titled `NT` rows say `Dox_0ng/mL`). Candidates: C1a plain WT line, no dox; "
               "C1b NTg-titled, 0 dox; C3/C4 NTg with dox; C5 `mTagBFP2` reporter with dox. "
               "Below age %d only C1a (and rapamycin-treated C2) exist." % int(tgf.age_description.min())),
        ("A5", "Schedules. `Dox_1000ng/mL, 3days` and `Dox_1000ng/mL, 3d-ON` are both used; "
               "whether they are the same exposure, and what precedes harvest after a `d-OFF` "
               "period, is not stated beyond the code."),
        ("A6", "Rapamycin dose: protocol text says 100 nM, field values say 100 ng/mL."),
        ("A7", "Label disagreements between title and `genotype`: " + "; ".join(
            "%d GSMs titled `%s` carry `genotype: %s`" % (n, a, b) for (a, b), n in
            diff[diff.title_token != "NTg"].groupby(["title_token", "ch_genotype"]).size().items() if a != "<none>")
            + ". (`TP53 (DN)` is the genotype of %d other GSMs.)" % int((hum.ch_genotype == "TP53 (DN)").sum())),
        ("A8", "Cell type. The design mentions \"iPSCs and their derivatives\"; `cell type` holds "
               "only `Fibroblast` and `iPSC`, so no derivative is labelled. The iPSC GSMs carry a "
               "`Donor age`, presumably of the source donor; the record does not say so. "
               "AG07141 contributes one Fibroblast and one iPSC GSM."),
        ("A9", "Batch. The `batch` field holds `LIB...` identifiers; the description adds "
               "`Batch 1`-`Batch 12` for 96 GSMs. Neither is defined as culture, plate or "
               "sequencing batch."),
        ("A10", "Depth. SRA spots are pre-downsampling. The processing text says libraries were "
                "downsampled to 20 million reads before counting (three LIB IDs exempt) and does "
                "not say reads or read pairs. GEO counts were not opened, so post-downsampling "
                "depth is unknown here."),
        ("A11", "Passage. Two GSMs have a title `P<n>` that disagrees with the described "
                "passage (listed in section 5). B6 titles use `P0`, `P+2`, `P+10`."),
        ("A12", "Molecule. %d GSMs are `total RNA`, %d `polyA RNA`; the extract protocol "
                "describes a polyA workflow for NEBNext libraries and never mentions total-RNA "
                "libraries. Every under-40 fibroblast GSM is `total RNA`." % (
                    len(tot), int((df.molecule == "polyA RNA").sum()))),
        ("A13", "Gene overlap. Not counted (section 7). Ensembl versions differ between the "
                "file and the ruler annotation for some genes."),
    ]
    for k, v in amb:
        say("- **%s.** %s" % (k, v))
    say()

    # ---- 9. plain statement --------------------------------------------------
    say("## 9. Plain statement (not a result)")
    say()
    yu = u40f[(u40f.ch_genotype == "WT") & (u40f.ch_treatment == "NT")]
    say("What a distance-to-young readout would need from this series:")
    say()
    say("**(a) Perturbations in donor-aged human fibroblasts: present.** %d human"
        % len(tgf))
    say("fibroblast GSMs carry %d distinct genotype values other than `WT`/`NTg`/reporter,"
        % tgf.ch_genotype.nunique())
    say("from %d lines at stated donor ages %s. %d of those values occur at more than one"
        % (tgf.ch_cell_line.nunique(), ", ".join(str(int(a)) for a in sorted(tgf.age_description.unique())),
           len(tg_multi)))
    say("donor age. The ages are stated in `!Sample_description`, not inferred.")
    say()
    say("**(b) Untreated young human fibroblasts of the same cell type, same series:")
    say("present but thin, and ambiguous.** Fibroblast GSMs with `genotype: WT`,")
    say("`treatment: NT` from donors under 40: %d (%s)." % (len(yu), "; ".join(
        "%s age %d, line %s, passage %s, batch %s, %s %s, %s SRA spots" % (
            r.gsm, r.age_description, r.ch_cell_line,
            "-" if np.isnan(r.passage_description) else int(r.passage_description),
            r.ch_batch, r.platform, r.molecule, format(int(per.loc[r.gsm, "spots"]), ","))
        for r in yu.itertuples())))
    say("They are labelled `cell type: Fibroblast`, as the perturbed GSMs are. What")
    say("makes them ambiguous as a reference, stated without resolving it:")
    say()
    tg_pass = tgf.passage_description.dropna()
    say("- no under-40 fibroblast GSM shares a `batch` value with any transgene GSM (section 4c);")
    say("- none is a reporter (`mTagBFP2`) or dox-exposed control: those exist only at ages %s;"
        % ", ".join(str(int(a)) for a in sorted(hf[hf.ch_genotype == "mTagBFP2"].age_description.unique())))
    say("- their described passages are %s; transgene GSMs' passages span %d-%d;" % (
        ", ".join(str(int(p)) for p in yu.passage_description.dropna()), int(tg_pass.min()), int(tg_pass.max())))
    say("- all %d are `total RNA` on %s; the transgene GSMs are %s;" % (
        len(yu), "/".join(sorted(yu.platform.unique())) if (yu.molecule == "total RNA").all() else "(mixed)",
        ", ".join("%s %s %d" % (a, b, n) for (a, b), n in tgf.groupby(["platform", "molecule"]).size().items())))
    older = hf[(hf.ch_genotype == "WT") & (hf.ch_treatment == "NT") & (hf.age_description >= 40)
               & (hf.age_description < tgf.age_description.min())]
    say("- untreated fibroblasts aged 40 to %d: %s; the youngest transgene line is %d." % (
        int(tgf.age_description.min()) - 1,
        ", ".join("age %d %d GSM(s)" % kv for kv in older.age_description.astype(int).value_counts().sort_index().items()) or "none",
        int(tgf.age_description.min())))
    yu_u = sorted(set(u40f.ch_batch))
    say("- the under-40 fibroblast batches (%s) hold only plain \"WT line\" GSMs, iPSC GSMs" % ", ".join(yu_u))
    say("  and rapamycin GSMs (section 4c).")
    hfc = hum[hum.group == "human fibroblast"]
    arm = hfc[hfc.control_candidate.isin(["C1a genotype WT, treatment NT, title says WT/plain line",
                                        "C2 genotype WT, rapamycin"])]
    say("- read together: the plain \"WT line\" GSMs (C1a, C2; %d fibroblast GSMs, %d lines," % (
        len(arm), arm.ch_cell_line.nunique()))
    say("  ages %d-%d) are all `%s`, in batches %s. The transgene GSMs come from" % (
        int(arm.age_description.min()), int(arm.age_description.max()),
        "/".join(sorted(arm.molecule.unique())), ", ".join(sorted(arm.ch_batch.unique()))))
    say("  ages %s; the untreated rows inside the screen's batches (C1b, polyA) come" % ", ".join(
        str(int(a)) for a in sorted(tgf.age_description.unique())))
    say("  only from ages %s. The record" % ", ".join(
        str(int(a)) for a in sorted(hfc[hfc.control_candidate.str.startswith("C1b")].age_description.unique())))
    say("  does not say whether the WT-line set and the screen are one experiment.")
    say()
    say("The %d iPSC GSMs from donors under 40 are a different cell type and are not" %
        int((u40.group == "human iPSC").sum()))
    say("counted under (b).")
    say()

    # ---- 10. what this inventory did beyond the brief ------------------------
    say("## 10. Disclosures")
    say()
    say("- No flag file: the task is an inventory with no pre-registered statistic.")
    say("  No count here is a gated result, so the chance-baseline, cell-count and")
    say("  recovery-power obligations of `CLAUDE.md` have nothing to attach to.")
    say("- Raw metadata stored under `data/raw/srsfinv/` (git-ignored); the script")
    say("  refetches a missing file. md5s are in section 0 and `inventory_meta.json`.")
    say("- Beyond the GEO record: the SRA run table (for depth; GEO has none) and a")
    say("  262,144-byte HTTP range read of the human counts file (for identifier type).")
    say("  Neither is an expression matrix download.")
    say("- The frozen ruler npz was opened for its `gene_id`/`ensembl` arrays only; no")
    say("  weight, mean or sd was read or used.")
    say("- Added beyond the brief: section 0b (scoping-parse comparison), section 3c")
    say("  (title vs genotype), section 4c (batch co-measurement), the C1-C6 shorthand.")
    say("- Decided here because the brief did not specify: \"donor\" is counted as a")
    say("  distinct `cell line` value; \"transgene\" excludes `WT`, `NTg`, `mTagBFP2`, `GFP`")
    say("  only where a sentence says so, and the full genotype table lists all values;")
    say("  \"under 40\" means described age < 40.")
    say("- Not covered: anything about donors younger than 55 receiving a transgene")
    say("  (there are none); the mouse arm beyond its sample table; expression values.")
    say()
    say("Machine-readable tables in `results/srsfinv/`: `samples.csv`,")
    say("`samples_human_annotated.csv`, `donors.csv`, `genotype_counts.csv`,")
    say("`genotype_by_age_fibroblast.csv`, `control_candidates_by_age.csv`,")
    say("`under40_samples.csv`, `batches.csv`, `treatments.csv`, `sra_depth_per_gsm.csv`,")
    say("`scoping_parse_comparison.csv`, `inventory_meta.json`.")

    meta = {"n_gsm": len(df), "n_human": len(hum), "n_mouse": len(mou),
            "n_fibroblast": int((hum.group == "human fibroblast").sum()),
            "n_ipsc": int((hum.group == "human iPSC").sum()),
            "n_genotypes": int(df.ch_genotype.nunique()),
            "n_cell_lines": int(hum.ch_cell_line.nunique()),
            "genotypes_multi_age_fibroblast": multi_age,
            "n_u40_fibroblast": int(len(u40f)),
            "n_u40_fibroblast_WT_NT": int(((u40f.ch_genotype == "WT") & (u40f.ch_treatment == "NT")).sum()),
            "counts_header_cols": len(header) - 1,
            "ruler_overlap_counted": False,
            "md5": {p.name: md5(p) for p in [soft_p, sra_p, head_p]}}
    (OUT / "inventory_meta.json").write_text(json.dumps(meta, indent=2), encoding="utf-8")
    hum.to_csv(OUT / "samples_human_annotated.csv", index=False)
    (OUT / "INVENTORY_SRSFINV.md").write_text("\n".join(lines) + "\n", encoding="utf-8")


if __name__ == "__main__":
    main()
