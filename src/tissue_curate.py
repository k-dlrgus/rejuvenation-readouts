"""T1 curate: re-assign tissue groups (the first-match regex put kidney 'cortex of kidney' in brain),
flag organoids / cancer atlases / cell-type subsets / multi-tissue incidentals, and write the
primary-atlas table used for shortlisting. All numbers come from the already-fetched CXG listing
(results/tissue/t1_cxg_all_human_datasets.csv); this script does not invent accessions.
"""
import re
import sys
from pathlib import Path

import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parent))
from tissue_common import TISSUE_DIR, Logger, PRIORITY_TISSUES  # noqa: E402

# More specific patterns; scored independently (not first-match).
TISSUE_PATTERNS = {
    "retina": r"\bretina|macula|fovea|retinal|choroid|chorioretinal|pigmented layer of retina|optic choroid",
    "brain": (
        r"\bbrain\b|cerebral|dorsolateral|prefrontal|hippocamp|entorhinal|cerebell|"
        r"striatum|substantia nigra|thalam|hypothalam|midbrain|white matter|"
        r"pons\b|medulla oblongata|amygdala|putamen|caudate nucleus|dentate gyrus|"
        r"basal ganglia|brodmann|spinal cord|cortical plate|neocortex|"
        r"temporal lobe|frontal lobe|parietal|occipital|DLPFC|pfc\b"
    ),
    "skeletal muscle": r"skeletal muscle|vastus|quadriceps|gastrocnemius|rectus abdominis|muscle of abdomen|muscle of pelvic|muscle tissue",
    "heart": r"\bheart\b|cardiac|ventricle|atrium|atrial|myocardium|septum|apex of heart",
    "kidney": r"kidney|renal medulla|renal papilla|renal pelvis|nephron|cortex of kidney",
    "liver": r"\bliver\b|hepatic|hepatocyte|caudate lobe|right lobe of liver|left lobe of liver",
}

# Tokens that claim a tissue we do NOT want as the primary label.
EXCLUDE_PRIMARY = re.compile(
    r"blood|pbmc|bone marrow|spleen|lymph node|thymus|tonsil|cell culture|organoid|"
    r"lung|colon|breast|skin|prostate|pancreas|intestine|stomach",
    re.I,
)
ORGANOID = re.compile(r"organoid|HNOCA|in vitro|iPSC-derived|directed differentiation", re.I)
CANCER = re.compile(
    r"cancer|carcinoma|glioblastoma|GBmap|melanoma|tumor|tumour|CRC|metastasis|"
    r"LuCA|lung cancer|checkpoint blockade",
    re.I,
)
SUBSET = re.compile(
    r"subset|bipolar cell|amacrine|ganglion cell|immune cells|T cells|myeloid|"
    r"lymphocyte|endothelial cells|epithelial cells|stromal cells|Schwann|"
    r"vascular associated|cholangiocyte cells|hepatocyte cells|mesenchyme cells|"
    r"non PT|non-tumor",
    re.I,
)
INTEGRATED = re.compile(r"harmonization|integrated|integration across|cell landscape|Tabula Sapiens", re.I)


def score_tissues(tissues: str) -> dict:
    labels = [x.strip() for x in str(tissues).split(";") if x.strip()]
    scores = {g: 0 for g in TISSUE_PATTERNS}
    for lab in labels:
        for g, pat in TISSUE_PATTERNS.items():
            if re.search(pat, lab, flags=re.I):
                scores[g] += 1
    return scores, labels


def assign_group(tissues: str):
    scores, labels = score_tissues(tissues)
    ranked = sorted(scores.items(), key=lambda kv: -kv[1])
    top_g, top_n = ranked[0]
    second_n = ranked[1][1]
    n_lab = max(len(labels), 1)
    if top_n == 0:
        return None, 0.0, scores, False
    frac = top_n / n_lab
    # kidney 'cortex of kidney' must beat brain 'cortex'
    return top_g, frac, scores, (top_n > second_n)


def main():
    log = Logger(TISSUE_DIR / "t1_curate.txt")
    src = TISSUE_DIR / "t1_cxg_all_human_datasets.csv"
    if not src.exists():
        raise SystemExit("run src/tissue_search_cxg.py first")
    df = pd.read_csv(src)
    log(f"loaded {len(df)} human CXG datasets from {src.name} (fetched, not assumed)")

    rows = []
    for r in df.itertuples(index=False):
        g, frac, scores, unique = assign_group(r.tissues)
        title = str(r.title or "")
        coll = str(r.collection_name or "")
        blob = f"{title} | {coll} | {r.tissues}"
        rows.append(dict(
            dataset_id=r.dataset_id, collection_id=r.collection_id, title=title,
            collection_name=coll, cell_count=r.cell_count, n_donors=r.n_donors,
            n_cell_types=r.n_cell_types, age_min=r.age_min, age_max=r.age_max,
            age_span=r.age_span, n_stages_numeric_years=r.n_stages_numeric_years,
            assays=r.assays, suspension_type=r.suspension_type, tissues=r.tissues,
            diseases=r.diseases, h5ad_url=r.h5ad_url, h5ad_bytes=r.h5ad_bytes,
            collection_doi=r.collection_doi, tissue_group=g, tissue_frac=frac,
            tissue_unique=unique, score_retina=scores["retina"], score_brain=scores["brain"],
            score_muscle=scores["skeletal muscle"], score_heart=scores["heart"],
            score_kidney=scores["kidney"], score_liver=scores["liver"],
            is_organoid=bool(ORGANOID.search(blob)),
            is_cancer=bool(CANCER.search(blob)),
            is_subset=bool(SUBSET.search(title)),
            is_integrated=bool(INTEGRATED.search(blob)),
            has_normal=("normal" in str(r.diseases).lower()),
        ))
    out = pd.DataFrame(rows)

    is_rna = out.assays.fillna("").str.contains(
        "10x|Seq-Well|Drop-seq|Smart-seq|sci-RNA|inDrop|CEL-seq|BD Rhapsody|snRNA|scRNA|microwell",
        case=False, regex=True,
    )
    cand = out[(out.tissue_group.isin(PRIORITY_TISSUES)) & is_rna].copy()
    log(f"[filter] human RNA in a priority tissue (rescored)           {len(out):,} -> {len(cand):,}")
    n = len(cand); cand = cand[cand.n_donors >= 10]
    log(f"[filter] >= 10 donors                                        {n:,} -> {len(cand):,}")
    n = len(cand); cand = cand[cand.age_span.fillna(0) >= 40]
    log(f"[filter] numeric age span >= 40 y                            {n:,} -> {len(cand):,}")
    n = len(cand); cand = cand[cand.n_cell_types >= 3]
    log(f"[filter] >= 3 cell types                                     {n:,} -> {len(cand):,}")
    n = len(cand); cand = cand[cand.tissue_frac >= 0.5]
    log(f"[filter] >= 50% of tissue labels belong to the assigned group {n:,} -> {len(cand):,}")
    n = len(cand); cand = cand[~cand.is_organoid]
    log(f"[filter] drop organoid / in-vitro                            {n:,} -> {len(cand):,}")
    n = len(cand); cand = cand[~cand.is_cancer]
    log(f"[filter] drop cancer / tumor atlases                         {n:,} -> {len(cand):,}")

    order = {g: i for i, g in enumerate(PRIORITY_TISSUES)}
    cand["priority"] = cand.tissue_group.map(order)
    cand = cand.sort_values(["priority", "is_subset", "is_integrated", "n_donors", "cell_count"],
                            ascending=[True, True, True, False, False])
    cand.to_csv(TISSUE_DIR / "t1_cxg_primary_candidates.csv", index=False)

    # keep a full rescored table too (including subsets / cancer — for the audit trail)
    out.to_csv(TISSUE_DIR / "t1_cxg_rescored.csv", index=False)

    pd.set_option("display.width", 280)
    pd.set_option("display.max_colwidth", 70)
    pd.set_option("display.max_rows", 400)
    cols = ["tissue_group", "dataset_id", "title", "collection_name", "cell_count", "n_donors",
            "n_cell_types", "age_min", "age_max", "assays", "suspension_type", "is_subset",
            "is_integrated", "has_normal", "h5ad_bytes"]
    log(f"\n=== primary (non-organoid, non-cancer, tissue-majority) candidates n={len(cand)} ===")
    for g in PRIORITY_TISSUES:
        sub = cand[cand.tissue_group == g]
        log(f"\n--- {g}: {len(sub)} datasets ({int((~sub.is_subset).sum())} full, {int(sub.is_subset.sum())} labelled subset) ---")
        if len(sub):
            log(sub[cols].to_string(index=False))
    log.close()
    return cand


if __name__ == "__main__":
    main()
