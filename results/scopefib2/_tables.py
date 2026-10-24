"""SCOPEFIB2: build the result CSVs from records already logged in this directory.

Inputs (all written by src/scopefib2_run.py in this session):
  url_records.json   - every HTTP request with URL + status
  geo_records.json   - GEO series text for GSE278990 and its five SubSeries

Outputs:
  hra014194_experiments.csv - one row per HRA experiment (arm-level record for bar A / C3)
  library_counts.csv        - bar B counts, from Tables S3A/S3B/S4E parsed in the browser pane
  datasets.csv              - the one dataset row, all mandatory columns from the flag
  browser_pane_urls.csv     - requests made through the built-in browser pane, not urllib

Hand-entered values are only those read from the article text or a supplementary table in
this session; each carries its source in the `source` column. No statistic is computed here.
"""
import csv
import json
from collections import OrderedDict
from pathlib import Path

HERE = Path(__file__).resolve().parent

ART = ("https://www.cell.com/cancer-cell/fulltext/S1535-6108(25)00552-5"
       " (read in the built-in browser pane; direct urllib fetch returned 403)")


def write(name, fieldnames, rows):
    with open(HERE / name, "w", newline="", encoding="utf-8") as fh:
        w = csv.DictWriter(fh, fieldnames=fieldnames)
        w.writeheader()
        w.writerows(rows)
    print(f"wrote {name}: {len(rows)} rows")


# ---------------------------------------------------------------- experiments
url_log = json.loads((HERE / "url_records.json").read_text(encoding="utf-8"))
runs = json.loads([x for x in url_log
                   if x["label"] == "hra_runinstudy_all"][0]["text"])["runViews"]
inds = json.loads([x for x in url_log
                   if x["label"] == "hra_indinstudy_all"][0]["text"])["individualViews"]
sample_ind = {i["sampleAcc"]: (i["accession"], i["name"], i["gender"], i["sampleDesc"])
              for i in inds}

agg = OrderedDict()
for r in runs:
    k = r["expName"]
    a = agg.setdefault(k, {"experiment": k, "title": r["expTitle"], "platform": r["platform"],
                           "strategy": r["strategyName"], "sample": r["sampleAcc"],
                           "runs": set(), "files": 0, "gb": 0.0, "run_titles": set()})
    a["runs"].add(r["runAcc"])
    a["files"] += 1
    a["gb"] += float(r["runFileSize"])
    a["run_titles"].add(r["runTitle"])

rows = []
for a in agg.values():
    ind = sample_ind.get(a["sample"], ("", "", "", ""))
    rows.append({
        "experiment": a["experiment"], "title": a["title"], "arm":
            "screen (perturbed, TGFb)" if a["title"].startswith("CRISPR_")
            else ("unperturbed state definition" if a["title"][:2] in ("P1", "P2", "P3")
                  else "downstream validation"),
        "platform": a["platform"], "strategy": a["strategy"],
        "sample": a["sample"], "individual": ind[0], "individual_name": ind[1],
        "gender": ind[2], "sample_description": ind[3],
        "n_runs": len(a["runs"]), "n_files": a["files"], "total_GB": round(a["gb"], 1),
        "run_title_prefixes": ";".join(sorted({t.split("_")[0] for t in a["run_titles"]})),
        "release_date": "2025-12-16",
    })
write("hra014194_experiments.csv", list(rows[0].keys()), rows)

# ---------------------------------------------------------------- bar B counts
# Parsed inside the browser pane from the published supplementary workbooks:
#   mmc4.xlsx = Table S3 (gene lists + sgRNA sequences), mmc5.xlsx = Table S4 (screen analysis).
# Control guides (AAVS1_1-47, NonTargetingControl_1-43) are excluded per D1.
bar_b = [
    {"set": "CRISPRi library, distinct gene targets", "n": 590,
     "source": "Table S3A gene list and Table S3B sgRNA table agree (mmc4.xlsx)"},
    {"set": "CRISPRa library, distinct gene targets", "n": 265,
     "source": "Table S3A gene list and Table S3B sgRNA table agree (mmc4.xlsx)"},
    {"set": "union of the two libraries (D1: one gene = one agent, any modality)", "n": 667,
     "source": "Table S3A/S3B (mmc4.xlsx); 188 genes appear in both libraries"},
    {"set": "genes in both libraries", "n": 188, "source": "Table S3A/S3B (mmc4.xlsx)"},
    {"set": "CRISPRi sgRNA rows (targeting, controls excluded)", "n": 2450,
     "source": "Table S3B (mmc4.xlsx)"},
    {"set": "CRISPRa sgRNA rows (targeting, controls excluded)", "n": 1150,
     "source": "Table S3B (mmc4.xlsx)"},
    {"set": "CRISPRi targets recovered after the >5-cell filter", "n": 432,
     "source": "article Results and Table S4E (mmc5.xlsx) agree"},
    {"set": "CRISPRa targets recovered after the >5-cell filter", "n": 210,
     "source": "article Results and Table S4E (mmc5.xlsx) agree"},
    {"set": "union of recovered targets", "n": 532,
     "source": "Table S4E (mmc5.xlsx); 110 recovered in both arms"},
    {"set": "recovered in both arms", "n": 110, "source": "Table S4E (mmc5.xlsx)"},
    {"set": "CRISPRi cells retained after QC", "n": 34890,
     "source": "article Results; Table S1A (mmc2.xlsx) agrees"},
    {"set": "CRISPRa cells retained after QC", "n": 39367,
     "source": "article Results; Table S1A (mmc2.xlsx) agrees"},
    {"set": "unperturbed cells retained after QC", "n": 67268,
     "source": "article Results; Table S1A rows sum to 26826+18722+13127+8593=67268"},
]
write("library_counts.csv", ["set", "n", "source"], bar_b)

# ---------------------------------------------------------------- dataset row
cols = ["accession", "organism", "cell type / source", "design",
        "n distinct perturbations (D1)", "readout", "units per perturbation",
        "perturbation_timing", "target_match", "batch_evidence",
        "A", "B", "C", "D", "E", "source"]
row = {
    "accession": "HRA014194 (BioProject PRJCA048014; GSA-Human, controlled access)",
    "organism": "human",
    "cell type / source": ("patient-derived fibroblasts (PDFs) from adjacent-tumour (NT) and "
                           "tumour-core (T) colorectal cancer tissue; the screen ran on one "
                           "pooled sample HRS2102648 = P1 NT + P4 NT/T (donors aged 35 F and "
                           "54 M, Table S1D); 8 PDF donors in the study, ages 35-79, 0 healthy "
                           "donors"),
    "design": ("direct-capture CRISPRi (dCas9-KRAB) and CRISPRa (dCas9-SunTag-VP64) Perturb-seq "
               "at MOI 0.5; TGFb 10 ng/mL added two days post-transduction, cells cultured three "
               "more days, FACS-sorted, 10x 3' HT v3.1 + Feature Barcode, CellPlex CMO "
               "multiplexing, 8 lanes (CRISPR_1 2 lanes, CRISPR_2 6 lanes), NovaSeq X Plus. "
               "Unperturbed PDFs +/- TGFb were profiled separately on the 10x 5' kit"),
    "n distinct perturbations (D1)": ("library 590 CRISPRi / 265 CRISPRa / 667 union; recovered "
                                      "432 / 210 / 532 after the >5-cell filter"),
    "readout": "scRNA-seq (10x 3' HT v3.1) with direct sgRNA capture",
    "units per perturbation": ("not stated; retention threshold is >5 cells per target. Mean "
                              "cells/recovered target = 34890/432 = 80.8 (CRISPRi) and "
                              "39367/210 = 187.5 (CRISPRa) - estimate, not the median D7 asks for"),
    "perturbation_timing": "before (transduction two days before TGFb)",
    "target_match": ("n/a - no unstimulated arm exists inside the screen. The nearest "
                     "unstimulated cells are unperturbed PDFs from P1/P2/P3: isogenic for P1, "
                     "cross_donor for P2 and P3; 0 healthy / 8 diseased PDF donors"),
    "batch_evidence": ("split - the unperturbed +/- TGFb arms are 10x 5' chemistry (Table S1A) "
                       "while both screen arms are 10x 3' HT v3.1 (Table S1A and STAR Methods); "
                       "different samples and individuals, and different sequencing-run names "
                       "(20210825-LJN-* vs P0926-*). Platform (NovaSeq X Plus), study and "
                       "release date (2025-12-16) are shared"),
    "A": "fail", "B": "pass", "C": "unverified", "D": "pass", "E": "pass",
    "source": ART,
}
write("datasets.csv", cols, [row])

# ------------------------------------------------- browser-pane request record
pane = [
    {"url": "https://www.cell.com/cancer-cell/fulltext/S1535-6108(25)00552-5",
     "status": "200", "what": "article full text (129,332 characters of article innerText)"},
    {"url": "https://www.cell.com/cms/10.1016/j.ccell.2025.12.018/attachment/"
            "0d1d4c8e-4950-4aae-a111-77bfb5f2ac76/mmc2.xlsx",
     "status": "200", "what": "Table S1 (18,108 bytes), parsed in the page"},
    {"url": "https://www.cell.com/cms/10.1016/j.ccell.2025.12.018/attachment/"
            "c0385583-38d9-4607-9990-d9f3fb8a4ba5/mmc4.xlsx",
     "status": "ok (status not captured)",
     "what": "Table S3 (125,281 bytes), parsed in the page"},
    {"url": "https://www.cell.com/cms/10.1016/j.ccell.2025.12.018/attachment/"
            "8d4319d8-73c5-4bd8-bbfc-f8fdc06b430c/mmc5.xlsx",
     "status": "ok (status not captured)",
     "what": "Table S4 (5,449,018 bytes), parsed in the page"},
]
write("browser_pane_urls.csv", ["url", "status", "what"], pane)
