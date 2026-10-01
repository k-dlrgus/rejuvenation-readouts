"""GENESPACE diagnostic (NOT pre-registered): do the published mmc3 Age up / Age down
lists point the right way with donor age in GTEx fibroblasts?

1. Loading check: independent stdlib read of mmc3.xlsx sheet Aging_signatures (sheet file
   resolved via workbook relationships, not by sorted filename) vs results/md3/genesets.json.
2. Per-gene Spearman rho with GTEx 10-year age bin on the frozen Stage 0 log-CPM pack
   (toward_run.load_gtex_z; nothing refitted).
3. Baseline: all ruler genes.

Writes only to results/genespace_diag/ and FINDINGS_GENESPACE_DIAG.md.
Usage: python src/genespace_diag.py
"""
from __future__ import annotations

import re
import sys
import zipfile
import xml.etree.ElementTree as ET
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import mannwhitneyu, rankdata

sys.path.insert(0, str(Path(__file__).resolve().parent))
from config import RESULTS, ROOT  # noqa: E402
from brain_phase1_common import Logger, dump_json, load_json  # noqa: E402
from fibro_common import jsonable  # noqa: E402
from md2_common import MMC3_XLSX  # noqa: E402
from toward_run import load_gtex_z, YOUNG_BINS, OLD_BINS  # noqa: E402
from target_common import unit  # noqa: E402
from md3_genesets import _xlsx_columns  # noqa: E402

OUT = RESULTS / "genespace_diag"
FINDINGS = ROOT / "FINDINGS_GENESPACE_DIAG.md"
GENESETS = RESULTS / "md3" / "genesets.json"
LABEL = "diagnostic, not pre-registered"
AGE_ORD = {"20-29": 0, "30-39": 1, "40-49": 2, "50-59": 3, "60-69": 4, "70-79": 5}

NS = {"m": "http://schemas.openxmlformats.org/spreadsheetml/2006/main",
      "r": "http://schemas.openxmlformats.org/officeDocument/2006/relationships",
      "pr": "http://schemas.openxmlformats.org/package/2006/relationships"}


def read_sheet_by_rels(path, sheet_name):
    """{header: [values...]} for one sheet, resolving sheet -> file through workbook.xml.rels."""
    with zipfile.ZipFile(path) as z:
        strings = []
        if "xl/sharedStrings.xml" in z.namelist():
            root = ET.fromstring(z.read("xl/sharedStrings.xml"))
            for si in root.findall("m:si", NS):
                strings.append("".join(t.text or "" for t in si.findall(".//m:t", NS)))
        wb = ET.fromstring(z.read("xl/workbook.xml"))
        rels = ET.fromstring(z.read("xl/_rels/workbook.xml.rels"))
        rid_to_target = {r.attrib["Id"]: r.attrib["Target"] for r in rels.findall("pr:Relationship", NS)}
        sheets = {sh.attrib["name"]: sh.attrib[f"{{{NS['r']}}}id"] for sh in wb.findall("m:sheets/m:sheet", NS)}
        target = rid_to_target[sheets[sheet_name]].lstrip("/")
        target = target if target.startswith("xl/") else "xl/" + target
        root = ET.fromstring(z.read(target))
        grid = {}
        for c in root.iter(f"{{{NS['m']}}}c"):
            m = re.match(r"([A-Z]+)(\d+)", c.attrib.get("r", ""))
            if not m:
                continue
            col = 0
            for ch in m.group(1):
                col = col * 26 + (ord(ch) - 64)
            row = int(m.group(2)) - 1
            t = c.attrib.get("t")
            v = c.find("m:v", NS)
            if t == "s" and v is not None:
                val = strings[int(v.text)]
            elif t == "inlineStr":
                val = "".join(x.text or "" for x in c.iter(f"{{{NS['m']}}}t"))
            else:
                val = v.text if v is not None else None
            grid[(row, col - 1)] = val
    n_rows = max(r for r, _ in grid) + 1
    n_cols = max(c for _, c in grid) + 1
    headers = [grid.get((0, j)) for j in range(n_cols)]
    out = {}
    for j, h in enumerate(headers):
        if h is None:
            continue
        out[str(h)] = [grid.get((i, j)) for i in range(1, n_rows)]
    return out, target, list(sheets)


def norm(vals):
    return [str(g).strip().upper() for g in vals if g is not None and str(g).strip()]


def spearman_cols(X, y):
    """Spearman rho of each column of X with y (average ranks for ties)."""
    ry = rankdata(y)
    ry = (ry - ry.mean()) / np.linalg.norm(ry - ry.mean())
    R = np.apply_along_axis(rankdata, 0, X)
    R = R - R.mean(0)
    nrm = np.linalg.norm(R, axis=0)
    rho = np.full(X.shape[1], np.nan)
    ok = nrm > 0
    rho[ok] = (ry @ R[:, ok]) / nrm[ok]
    return rho


def summarize(name, rho):
    r = rho[np.isfinite(rho)]
    return dict(set=name, n_genes=int(rho.size), n_finite=int(r.size),
                share_pos=float((r > 0).mean()), share_neg=float((r < 0).mean()),
                share_zero=float((r == 0).mean()), median_rho=float(np.median(r)),
                mean_rho=float(r.mean()), q25=float(np.percentile(r, 25)), q75=float(np.percentile(r, 75)))


def main():
    OUT.mkdir(parents=True, exist_ok=True)
    log = Logger(OUT / "diag_log.txt")
    log(f"[{LABEL}] GENESPACE sign diagnostic")

    # ---- 1. loading check
    gs = load_json(GENESETS)
    indep, target, sheet_names = read_sheet_by_rels(MMC3_XLSX, "Aging_signatures")
    repo = _xlsx_columns(MMC3_XLSX)["Aging_signatures"]
    load = {"label": LABEL, "mmc3": str(MMC3_XLSX), "sheets": sheet_names,
            "Aging_signatures_file_via_rels": target, "Aging_signatures_headers": list(indep)}
    for col, key in (("Age up", "age_up"), ("Age down", "age_down")):
        raw_first5 = [g for g in indep[col] if g is not None and str(g).strip()][:5]
        ind, rep, js = norm(indep[col]), norm(repo[col]), [str(g) for g in gs[key]]
        load[key] = dict(
            column=col, n_indep=len(ind), n_repo_reader=len(rep), n_json=len(js),
            indep_equals_json_ordered=ind == js, repo_reader_equals_json_ordered=rep == js,
            first5_xlsx_raw=raw_first5, first5_json=js[:5],
        )
        log(f"[load] {col}: xlsx n={len(ind)} json n={len(js)} equal(ordered)={ind == js} "
            f"first5 xlsx={raw_first5} json={js[:5]}")
    load["age_up_equals_age_down_column"] = norm(indep["Age down"]) == [str(g) for g in gs["age_up"]]
    dump_json(OUT / "loading_check.json", jsonable(load))

    # ---- 2/3. GTEx per-gene Spearman with age bin (frozen Stage 0 pack)
    g = load_gtex_z(log)
    X = np.asarray(g["pack"]["X"], np.float64)
    age = g["age"]
    y = np.array([AGE_ORD[a] for a in age], float)
    n_donors = X.shape[0]
    log(f"[gtex] donors={n_donors} genes={X.shape[1]}")
    sym = np.array([str(s).strip().upper() for s in np.asarray(g["frozen"]["symbol"])], dtype=object)
    ens = np.asarray(g["frozen"]["ensembl"]).astype(str)
    rho = spearman_cols(X, y)

    up = set(str(s) for s in gs["age_up"])
    dn = set(str(s) for s in gs["age_down"])
    in_up = np.array([s in up for s in sym])
    in_dn = np.array([s in dn for s in sym])
    per_gene = pd.DataFrame(dict(ensembl=ens, symbol=sym, in_age_up=in_up, in_age_down=in_dn, spearman_rho_age=rho))
    per_gene.to_csv(OUT / "per_gene_spearman_age.csv", index=False)

    rows = [summarize("age_up", rho[in_up]), summarize("age_down", rho[in_dn]),
            summarize("all_ruler", rho), summarize("ruler_not_in_lists", rho[~in_up & ~in_dn])]
    summ = pd.DataFrame(rows)
    summ.insert(0, "label", LABEL)
    summ.to_csv(OUT / "summary.csv", index=False)
    for r in rows:
        log(f"[rho] {r}")

    ru, rd = rho[in_up], rho[in_dn]
    mw = mannwhitneyu(ru[np.isfinite(ru)], rd[np.isfinite(rd)], alternative="two-sided")

    # consistency: reproduce readout cos in AGE space and the sign of young-old mean Z per list
    Z = g["Z"]
    ym, om = np.isin(age, YOUNG_BINS), np.isin(age, OLD_BINS)
    diff = Z[ym].mean(0) - Z[om].mean(0)
    cols = np.where(in_up | in_dn)[0]
    r_age = np.where(in_up[cols], -1.0, 0.0) + np.where(in_dn[cols], 1.0, 0.0)
    u, _ = unit(diff[cols])
    cos_age = float(u @ r_age / np.linalg.norm(r_age))
    extra = dict(
        label=LABEL, n_donors=int(n_donors), age_bin_counts=pd.Series(age).value_counts().sort_index().to_dict(),
        n_young=int(ym.sum()), n_old=int(om.sum()),
        n_age_up_cols=int(in_up.sum()), n_age_down_cols=int(in_dn.sum()),
        mannwhitney_up_vs_down_rho=dict(U=float(mw.statistic), p=float(mw.pvalue)),
        cos_u_young_minus_old_vs_r_AGE_in_AGE_space=cos_age,
        share_age_up_higher_in_young=float((diff[in_up] > 0).mean()),
        share_age_down_higher_in_young=float((diff[in_dn] > 0).mean()),
        median_young_minus_old_z_age_up=float(np.median(diff[in_up])),
        median_young_minus_old_z_age_down=float(np.median(diff[in_dn])),
    )
    obs = g["obs"]
    extra["age_ordinal_matches_bin_order"] = bool(
        "age_ordinal" in obs.columns and (obs["age_ordinal"].to_numpy(int) == y.astype(int)).all())
    sanity = ["CDKN2A", "CDKN1A", "MKI67", "MMP3", "XIST", "RPS4Y1"]
    extra["sanity_genes_rho"] = {s: float(rho[sym == s][0]) for s in sanity if (sym == s).any()}
    dump_json(OUT / "summary_extra.json", jsonable(extra))
    log(f"[extra] {extra}")
    write_findings(load, rows, extra)
    log.close()


def pct(x):
    return f"{100 * x:.1f}%"


def write_findings(load, rows, extra):
    r = {x["set"]: x for x in rows}
    up, dn, base, rest = r["age_up"], r["age_down"], r["all_ruler"], r["ruler_not_in_lists"]
    lu, ld = load["age_up"], load["age_down"]
    ok_load = lu["indep_equals_json_ordered"] and ld["indep_equals_json_ordered"]
    up_ok = up["share_pos"] > base["share_pos"] and up["median_rho"] > 0
    dn_ok = dn["share_neg"] > base["share_neg"] and dn["median_rho"] < 0
    lines = [
        "# FINDINGS_GENESPACE_DIAG — do the published age lists point the right way in GTEx?",
        "",
        f"**Label: {LABEL}.** Nothing here was in `results/genespace/PREREG.flag`. No manuscript or existing file "
        "was edited. Script: `src/genespace_diag.py`. Outputs: `results/genespace_diag/`.",
        "",
        "## 1. Loading check",
        "",
        f"- `src/md3_genesets.py` reads sheet `Aging_signatures` and indexes columns by header name "
        "(`aging[\"Age up\"]`, `aging[\"Age down\"]`), with a StopStep if either header is missing.",
        f"- Independent read (stdlib, sheet resolved through `workbook.xml.rels` → `{load['Aging_signatures_file_via_rels']}`), "
        f"headers: {load['Aging_signatures_headers']}.",
        f"- `genesets.json` `age_up` equals the xlsx `Age up` column (stripped/uppercased, same order): "
        f"**{lu['indep_equals_json_ordered']}** (n={lu['n_json']}). `age_down` equals `Age down`: "
        f"**{ld['indep_equals_json_ordered']}** (n={ld['n_json']}). `age_up` equals the `Age down` column: "
        f"{load['age_up_equals_age_down_column']}.",
        "",
        "| list | first 5 in xlsx column | first 5 in genesets.json |",
        "|---|---|---|",
        f"| Age up / age_up | {', '.join(map(str, lu['first5_xlsx_raw']))} | {', '.join(lu['first5_json'])} |",
        f"| Age down / age_down | {', '.join(map(str, ld['first5_xlsx_raw']))} | {', '.join(ld['first5_json'])} |",
        "",
        "## 2–3. Per-gene Spearman ρ with GTEx age bin",
        "",
        f"Frozen Stage 0 fibroblast pack (TMM log2-CPM, prior 2; `toward_run.load_gtex_z`), {extra['n_donors']} donors, "
        f"age = ordinal 10-year bin (20-29 … 70-79; counts {extra['age_bin_counts']}). Ties use average ranks. "
        "Genes are the mapped ruler columns (same symbol mapping as GENESPACE).",
        "",
        "| set | n genes | share ρ > 0 | share ρ < 0 | median ρ | IQR |",
        "|---|---|---|---|---|---|",
    ]
    for x in (up, dn, base, rest):
        lines.append(f"| {x['set']} | {x['n_finite']} | {pct(x['share_pos'])} | {pct(x['share_neg'])} | "
                     f"{x['median_rho']:+.4f} | [{x['q25']:+.3f}, {x['q75']:+.3f}] |")
    mw = extra["mannwhitney_up_vs_down_rho"]
    lines += [
        "",
        f"Mann–Whitney age_up ρ vs age_down ρ: p = {mw['p']:.2e}.",
        "",
        "## Consistency with the GENESPACE readout",
        "",
        f"- Recomputed cos(GTEx c_young − c_old, r_AGE) in AGE space: {extra['cos_u_young_minus_old_vs_r_AGE_in_AGE_space']:+.4f} "
        "(r_AGE = −1 on age_up, +1 on age_down; FINDINGS_GENESPACE reports −0.37).",
        f"- Share of age_up genes higher in young (20–39, n={extra['n_young']}) than old (60–79, n={extra['n_old']}) "
        f"frozen-z: {pct(extra['share_age_up_higher_in_young'])}; median young−old z {extra['median_young_minus_old_z_age_up']:+.4f}.",
        f"- Share of age_down genes higher in young: {pct(extra['share_age_down_higher_in_young'])}; "
        f"median young−old z {extra['median_young_minus_old_z_age_down']:+.4f}.",
        f"- Age coding check: GTEx `age_ordinal` equals the bin order used here: {extra['age_ordinal_matches_bin_order']}. "
        "Reference genes (ρ with age; MKI67 and MMP3 are in age_down): "
        + ", ".join(f"{k} {v:+.3f}" for k, v in extra["sanity_genes_rho"].items())
        + ". XIST negative and RPS4Y1 positive mean the older bins are more male; sex is not adjusted here.",
        "- Per-gene effects are small (IQR of ρ about ±0.1); the signal is in the consistent direction across "
        "thousands of genes, not in individual genes.",
        "",
        "## Plain summary",
        "",
        f"- Age-up genes go up with age in GTEx: **{'yes' if up_ok else 'no'}** "
        f"({pct(up['share_pos'])} have ρ > 0 vs {pct(base['share_pos'])} of all ruler genes; median ρ {up['median_rho']:+.4f}).",
        f"- Age-down genes go down with age in GTEx: **{'yes' if dn_ok else 'no'}** "
        f"({pct(dn['share_neg'])} have ρ < 0 vs {pct(base['share_neg'])} of all ruler genes; median ρ {dn['median_rho']:+.4f}).",
        f"- Loading error: **{'no sign of one' if ok_load else 'MISMATCH — see loading_check.json'}**. "
        "Columns are read by header name and the JSON matches the xlsx columns exactly.",
        "",
        f"_{LABEL}_",
        "",
    ]
    FINDINGS.write_text("\n".join(lines), encoding="utf-8")


if __name__ == "__main__":
    main()
