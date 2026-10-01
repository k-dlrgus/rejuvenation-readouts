# PROGRESS_MD

**STOP status:** STAGE1_DONE_GATE_CLOSED
**Next action:** STOP. Gate closed: gene list public (mmc3 n=205), computation not in public sources. Stage 2 does not run until STAR Methods or author code is public. Do not approximate. Do not re-run the searches in PROGRESS_MD.md.

## Seeds / gates

- seed `20260914`  boot `20260918`  n_perm=200  n_boot=200  n_random=200
- primary metric: Spearman ρ with permutation null; uncalibrated R² is never a gate
- frozen ruler: `<repo>\results\fibro\frozen_ruler_ridge_raw.npz` exists=True
- PREREG_STAGE2.flag exists=True
- gate.json exists=True
- GSE325735: out of scope
- Refit nothing. Do not circumvent a paywall.

## Gate status

- open=False  reason=MD gene list recovered from open supplement mmc3.xlsx (n=205 genes, sheet MD_signatures, column 'MD score'; https://ars.els-cdn.com/content/image/1-s2.0-S0092867425008530-mmc3.xlsx). Computation of the score from those genes is not in any public source retrieved this session. STOP after Stage 1.
- gene_list_recovered=True
- computation_recovered=False
- missing=["MD score computation (per cell vs per sample; mean / ssGSEA / AddModuleScore / other). Open supplement mmc3.xlsx lists genes under column 'MD score' but contains no formula. GEO Sample_data_processing is Cell Ranger v7.1.0 + Seurat v5.0.0 sctransform v2, and does not mention an MD score. Would be in STAR Methods of the paywalled article (https://www.cell.com/cell/fulltext/S0092-8674(25)00853-0 , HTTP 403 this session) or in an author code repository (GitHub repo search total_count=0)."]

## Searches already run (do not repeat)

- `40816266[PMID] in db=gds` via eutils esearch gds → n=0 ids=[]
- `40816266[pmid] in db=pmc` via eutils esearch pmc → n=0 ids=[]
- `biorxiv_title` via https://www.biorxiv.org/search/prevalent%20mesenchymal%20drift%20in%20aging%20jcode%3Abiorxiv%20numresults%3A25%20sort%3Arelevance-rank → retrieved status=403 nbytes=6187 results_token=None empty_marker=False
- `biorxiv_author_lu` via https://www.biorxiv.org/search/author%3ALu%20mesenchymal%20drift%20jcode%3Abiorxiv%20numresults%3A25 → retrieved status=403 nbytes=6039 results_token=None empty_marker=False
- `biorxiv_phrase` via https://www.biorxiv.org/search/%22mesenchymal%20drift%22%20jcode%3Abiorxiv%20numresults%3A25%20sort%3Arelevance-rank → retrieved status=403 nbytes=6108 results_token=None empty_marker=False
- `medrxiv_phrase` via https://www.medrxiv.org/search/%22mesenchymal%20drift%22%20jcode%3Amedrxiv%20numresults%3A25 → retrieved status=403 nbytes=5994 results_token=None empty_marker=False
- `researchsquare_title` via https://www.researchsquare.com/browse?q=Prevalent%20mesenchymal%20drift%20in%20aging → retrieved status=200 nbytes=234130 results_token=None empty_marker=False
- `osf_preprints` via https://osf.io/preprints/discover?q=Prevalent%20mesenchymal%20drift → retrieved status=200 nbytes=4207 results_token=None empty_marker=False
- `ssrn` via https://www.ssrn.com/index.cfm/en/search/?q=Prevalent%20mesenchymal%20drift → retrieved status=403 nbytes=5989 results_token=None empty_marker=False
- `TITLE:"Prevalent mesenchymal drift" preprint` via europepmc → hitCount=0
- `GSE297234` via github repos → total_count=0 items=[]
- `mesenchymal_drift` via github repos → total_count=0 items=[]
- `lu_mesenchymal` via github repos → total_count=0 items=[]
- `altos_md` via github repos → total_count=0 items=[]

## Finished cells (from disk)

### sources (`sources.csv`)
- source=GEO GSE297234 SOFT series url=https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE297234&targ=self&form=text&view=brief retrieval_date=2026-09-19 found=retrieved. keys=Series_contact_address,Series_contact_city,Series_contact_country,Series_contact_institute,Series_contact_name,Series_contact_state,Series_contact_zip/postal_code,Series_contributor,Series_geo_accession,Series_last_update_date,Series_overall_design,Series_platform_id,Series_platform_organism,Series_platform_taxid,Series_pubmed_id,Series_relation,Series_sample_id,Series_sample_organism,Series_sample_taxid,Series_status,Series_submission_date,Series_summary,Series_supplementary_file,Series_title,Series_type,Series_web_link http_status=nan nbytes=nan path=<repo>\results\md\raw\geo_GSE297234_series.soft.txt error=nan note=nan n_samples=8.0 paywall_interstitial=nan
- source=GEO GSE297234 SOFT series full url=https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE297234&targ=self&form=text&view=full retrieval_date=2026-09-19 found=retrieved. keys=Series_contact_address,Series_contact_city,Series_contact_country,Series_contact_institute,Series_contact_name,Series_contact_state,Series_contact_zip/postal_code,Series_contributor,Series_geo_accession,Series_last_update_date,Series_overall_design,Series_platform_id,Series_platform_organism,Series_platform_taxid,Series_pubmed_id,Series_relation,Series_sample_id,Series_sample_organism,Series_sample_taxid,Series_status,Series_submission_date,Series_summary,Series_supplementary_file,Series_title,Series_type,Series_web_link relation=['BioProject: https://www.ncbi.nlm.nih.gov/bioproject/PRJNA1263211'] web_link=['https://doi.org/10.1016/j.cell.2025.07.031'] http_status=200.0 nbytes=2960.0 path=<repo>\results\md\raw\geo_GSE297234_series_full.soft.txt error=nan note=nan n_samples=nan paywall_interstitial=nan
- source=GEO GSE297234 SOFT GSM url=https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE297234&targ=gsm&form=text&view=brief retrieval_date=2026-09-19 found=retrieved n_gsm=8 characteristic_keys=['batch', 'cell line', 'cell type', 'genotype', 'tissue', 'treatment'] sample_supplementary_n=8 http_status=nan nbytes=nan path=<repo>\results\md\raw\geo_GSE297234_gsm.soft.txt error=nan note=nan n_samples=nan paywall_interstitial=nan
- source=GEO GSE297234 HTML url=https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE297234 retrieval_date=2026-09-19 found=retrieved http_status=200.0 nbytes=21383.0 path=<repo>\results\md\raw\geo_GSE297234_html.txt error=nan note=nan n_samples=nan paywall_interstitial=nan
- source=GEO GSE297233 SOFT series url=https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE297233&targ=self&form=text&view=brief retrieval_date=2026-09-19 found=retrieved. keys=Series_contact_address,Series_contact_city,Series_contact_country,Series_contact_institute,Series_contact_name,Series_contact_state,Series_contact_zip/postal_code,Series_contributor,Series_geo_accession,Series_last_update_date,Series_overall_design,Series_platform_id,Series_platform_organism,Series_platform_taxid,Series_pubmed_id,Series_relation,Series_sample_id,Series_sample_organism,Series_sample_taxid,Series_status,Series_submission_date,Series_summary,Series_supplementary_file,Series_title,Series_type,Series_web_link http_status=nan nbytes=nan path=<repo>\results\md\raw\geo_GSE297233_series.soft.txt error=nan note=nan n_samples=8.0 paywall_interstitial=nan
- source=GEO GSE297233 SOFT series full url=https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE297233&targ=self&form=text&view=full retrieval_date=2026-09-19 found=retrieved. keys=Series_contact_address,Series_contact_city,Series_contact_country,Series_contact_institute,Series_contact_name,Series_contact_state,Series_contact_zip/postal_code,Series_contributor,Series_geo_accession,Series_last_update_date,Series_overall_design,Series_platform_id,Series_platform_organism,Series_platform_taxid,Series_pubmed_id,Series_relation,Series_sample_id,Series_sample_organism,Series_sample_taxid,Series_status,Series_submission_date,Series_summary,Series_supplementary_file,Series_title,Series_type,Series_web_link relation=['BioProject: https://www.ncbi.nlm.nih.gov/bioproject/PRJNA1263275'] web_link=['https://doi.org/10.1016/j.cell.2025.07.031'] http_status=200.0 nbytes=2688.0 path=<repo>\results\md\raw\geo_GSE297233_series_full.soft.txt error=nan note=nan n_samples=nan paywall_interstitial=nan
- source=GEO GSE297233 SOFT GSM url=https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE297233&targ=gsm&form=text&view=brief retrieval_date=2026-09-19 found=retrieved n_gsm=8 characteristic_keys=['batch', 'cell line', 'cell type', 'genotype', 'tissue', 'treatment'] sample_supplementary_n=8 http_status=nan nbytes=nan path=<repo>\results\md\raw\geo_GSE297233_gsm.soft.txt error=nan note=nan n_samples=nan paywall_interstitial=nan
- source=GEO GSE297233 HTML url=https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE297233 retrieval_date=2026-09-19 found=retrieved http_status=200.0 nbytes=21383.0 path=<repo>\results\md\raw\geo_GSE297233_html.txt error=nan note=nan n_samples=nan paywall_interstitial=nan
- source=NCBI GDS esearch by PMID url=https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi?db=gds&term=40816266%5BPMID%5D&retmode=json&retmax=50&tool=age_identity_separability&email=none%40example.invalid retrieval_date=2026-09-19 found=retrieved http_status=200.0 nbytes=357.0 path=<repo>\results\md\raw\eutils_gds_pmid.json error=nan note=eutils n_samples=nan paywall_interstitial=nan
- source=BioProject lookup PRJNA1263211 url=https://www.ncbi.nlm.nih.gov/bioproject/?term=PRJNA1263211 retrieval_date=2026-09-19 found=retrieved http_status=200.0 nbytes=60372.0 path=<repo>\results\md\raw\bioproject_PRJNA1263211.html error=nan note=nan n_samples=nan paywall_interstitial=nan
- source=NCBI BioProject esearch GSE297234 url=https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi?db=bioproject&term=GSE297234&retmode=json&retmax=20&tool=age_identity_separability&email=none%40example.invalid retrieval_date=2026-09-19 found=retrieved http_status=200.0 nbytes=297.0 path=<repo>\results\md\raw\eutils_bioproject_gse297234.json error=nan note=eutils n_samples=nan paywall_interstitial=nan
- source=NCBI SRA esearch GSE297234 url=https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi?db=sra&term=GSE297234&retmode=json&retmax=20&tool=age_identity_separability&email=none%40example.invalid retrieval_date=2026-09-19 found=retrieved http_status=200.0 nbytes=347.0 path=<repo>\results\md\raw\eutils_sra_gse297234.json error=nan note=eutils n_samples=nan paywall_interstitial=nan
- source=PubMed XML url=https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=pubmed&id=40816266&rettype=xml&retmode=xml&tool=age_identity_separability&email=none%40example.invalid retrieval_date=2026-09-19 found=retrieved http_status=200.0 nbytes=8321.0 path=<repo>\results\md\raw\pubmed_40816266.xml error=nan note=eutils n_samples=nan paywall_interstitial=nan
- source=PubMed parsed fields url=https://pubmed.ncbi.nlm.nih.gov/40816266/ retrieval_date=2026-09-19 found=abstract=yes doi=10.1016/j.cell.2025.07.031 pmc=None n_authors=9 journal=Cell 2025;188(21):5895-5911.e17 http_status=nan nbytes=nan path=nan error=nan note=nan n_samples=nan paywall_interstitial=nan
- source=PubMed HTML url=https://pubmed.ncbi.nlm.nih.gov/40816266/ retrieval_date=2026-09-19 found=NOT RETRIEVED None http_status=203.0 nbytes=5565.0 path=nan error=nan note=nan n_samples=nan paywall_interstitial=nan
- source=Cell journal landing page url=https://www.cell.com/cell/fulltext/S0092-8674(25)00853-0 retrieval_date=2026-09-19 found=retrieved; landing page has highlights + summary; paywall interstitial not seen in excerpt; STAR Methods listed as a site nav item, not article methods text http_status=403.0 nbytes=5871.0 path=nan error=nan note=nan n_samples=nan paywall_interstitial=False
- source=DOI resolver url=https://doi.org/10.1016/j.cell.2025.07.031 retrieval_date=2026-09-19 found=retrieved status=200 final=https://linkinghub.elsevier.com/retrieve/pii/S0092867425008530 http_status=200.0 nbytes=2973.0 path=<repo>\results\md\raw\doi_resolve.html error=nan note=nan n_samples=nan paywall_interstitial=nan
- source=ScienceDirect landing url=https://www.sciencedirect.com/science/article/pii/S0092867425008530 retrieval_date=2026-09-19 found=NOT RETRIEVED None http_status=403.0 nbytes=832805.0 path=nan error=nan note=nan n_samples=nan paywall_interstitial=nan
- source=PMC esearch by PMID url=https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi?db=pmc&term=40816266%5Bpmid%5D&retmode=json&retmax=20&tool=age_identity_separability&email=none%40example.invalid retrieval_date=2026-09-19 found=retrieved pmc_ids=[] (empty list = no PMC deposit found by this query) http_status=200.0 nbytes=276.0 path=<repo>\results\md\raw\eutils_pmc_pmid.json error=nan note=eutils n_samples=nan paywall_interstitial=nan
- source=elink PubMed→PMC url=https://eutils.ncbi.nlm.nih.gov/entrez/eutils/elink.fcgi?dbfrom=pubmed&db=pmc&id=40816266&retmode=json&tool=age_identity_separability&email=none%40example.invalid retrieval_date=2026-09-19 found=retrieved http_status=200.0 nbytes=429.0 path=<repo>\results\md\raw\eutils_elink_pubmed_pmc.json error=nan note=eutils n_samples=nan paywall_interstitial=nan
- source=Europe PMC core JSON url=https://www.ebi.ac.uk/europepmc/webservices/rest/search?query=EXT_ID:40816266&resultType=core&format=json retrieval_date=2026-09-19 found=retrieved isOpenAccess='N' pmcid=None hasTextMinedTerms='Y' fullTextUrlList_n=1 license=None http_status=200.0 nbytes=6929.0 path=<repo>\results\md\raw\europepmc_core.json error=nan note=nan n_samples=nan paywall_interstitial=nan
- source=Europe PMC fullTextUrl (subscription required doi) url=https://doi.org/10.1016/j.cell.2025.07.031 retrieval_date=2026-09-19 found=retrieved http_status=200.0 nbytes=2973.0 path=<repo>\results\md\raw\epmc_ft__doi.org__j.cell.2025.07.031 error=nan note=nan n_samples=nan paywall_interstitial=nan
- source=Europe PMC HAS_PREPRINT AND this PMID url=https://www.ebi.ac.uk/europepmc/webservices/rest/search?query=HAS_PREPRINT:Y%20AND%20EXT_ID:40816266&format=json retrieval_date=2026-09-19 found=retrieved http_status=200.0 nbytes=197.0 path=<repo>\results\md\raw\europepmc_has_preprint.json error=nan note=nan n_samples=nan paywall_interstitial=nan
- source=Unpaywall url=https://api.unpaywall.org/v2/10.1016/j.cell.2025.07.031?email=none@example.invalid retrieval_date=2026-09-19 found=retrieved oa_status='closed' is_oa=False has_repository_copy=False best_oa=None http_status=200.0 nbytes=2931.0 path=<repo>\results\md\raw\unpaywall.json error=nan note=nan n_samples=nan paywall_interstitial=nan
- source=OpenAlex url=https://api.openalex.org/works/https://doi.org/10.1016/j.cell.2025.07.031 retrieval_date=2026-09-19 found=retrieved is_oa=False oa_status='closed' oa_url=None http_status=200.0 nbytes=22498.0 path=<repo>\results\md\raw\openalex.json error=nan note=nan n_samples=nan paywall_interstitial=nan
- source=Crossref url=https://api.crossref.org/works/10.1016/j.cell.2025.07.031 retrieval_date=2026-09-19 found=retrieved n_license=7 relation_keys=[] has_abstract=False http_status=200.0 nbytes=38660.0 path=<repo>\results\md\raw\crossref.json error=nan note=nan n_samples=nan paywall_interstitial=nan
- source=preprint search biorxiv_title url=https://www.biorxiv.org/search/prevalent%20mesenchymal%20drift%20in%20aging%20jcode%3Abiorxiv%20numresults%3A25%20sort%3Arelevance-rank retrieval_date=2026-09-19 found=retrieved status=403 nbytes=6187 results_token=None empty_marker=False http_status=403.0 nbytes=6187.0 path=nan error=nan note=nan n_samples=nan paywall_interstitial=nan
- source=preprint search biorxiv_author_lu url=https://www.biorxiv.org/search/author%3ALu%20mesenchymal%20drift%20jcode%3Abiorxiv%20numresults%3A25 retrieval_date=2026-09-19 found=retrieved status=403 nbytes=6039 results_token=None empty_marker=False http_status=403.0 nbytes=6039.0 path=nan error=nan note=nan n_samples=nan paywall_interstitial=nan
- source=preprint search biorxiv_phrase url=https://www.biorxiv.org/search/%22mesenchymal%20drift%22%20jcode%3Abiorxiv%20numresults%3A25%20sort%3Arelevance-rank retrieval_date=2026-09-19 found=retrieved status=403 nbytes=6108 results_token=None empty_marker=False http_status=403.0 nbytes=6108.0 path=nan error=nan note=nan n_samples=nan paywall_interstitial=nan
- source=preprint search medrxiv_phrase url=https://www.medrxiv.org/search/%22mesenchymal%20drift%22%20jcode%3Amedrxiv%20numresults%3A25 retrieval_date=2026-09-19 found=retrieved status=403 nbytes=5994 results_token=None empty_marker=False http_status=403.0 nbytes=5994.0 path=nan error=nan note=nan n_samples=nan paywall_interstitial=nan
- source=preprint search researchsquare_title url=https://www.researchsquare.com/browse?q=Prevalent%20mesenchymal%20drift%20in%20aging retrieval_date=2026-09-19 found=retrieved status=200 nbytes=234130 results_token=None empty_marker=False http_status=200.0 nbytes=234130.0 path=<repo>\results\md\raw\preprint_researchsquare_title.html error=nan note=nan n_samples=nan paywall_interstitial=nan
- source=preprint search osf_preprints url=https://osf.io/preprints/discover?q=Prevalent%20mesenchymal%20drift retrieval_date=2026-09-19 found=retrieved status=200 nbytes=4207 results_token=None empty_marker=False http_status=200.0 nbytes=4207.0 path=<repo>\results\md\raw\preprint_osf_preprints.html error=nan note=nan n_samples=nan paywall_interstitial=nan
- source=preprint search ssrn url=https://www.ssrn.com/index.cfm/en/search/?q=Prevalent%20mesenchymal%20drift retrieval_date=2026-09-19 found=retrieved status=403 nbytes=5989 results_token=None empty_marker=False http_status=403.0 nbytes=5989.0 path=nan error=nan note=nan n_samples=nan paywall_interstitial=nan
- source=Europe PMC preprint title search url=https://www.ebi.ac.uk/europepmc/webservices/rest/search?query=TITLE:%22Prevalent%20mesenchymal%20drift%22%20AND%20(SRC:PPR%20OR%20PUB_TYPE:preprint)&format=json retrieval_date=2026-09-19 found=retrieved hitCount=0 http_status=200.0 nbytes=235.0 path=<repo>\results\md\raw\europepmc_preprint_title.json error=nan note=nan n_samples=nan paywall_interstitial=nan
- source=Europe PMC preprint author Lu + mesenchymal drift url=https://www.ebi.ac.uk/europepmc/webservices/rest/search?query=%22mesenchymal%20drift%22%20AND%20AUTHOR:Lu%20AND%20(SRC:PPR%20OR%20PUB_TYPE:preprint)&format=json retrieval_date=2026-09-19 found=retrieved hitCount=0 http_status=200.0 nbytes=233.0 path=<repo>\results\md\raw\europepmc_preprint_lu.json error=nan note=nan n_samples=nan paywall_interstitial=nan
- source=GitHub repo search GSE297234 url=https://api.github.com/search/repositories?q=GSE297234 retrieval_date=2026-09-19 found=retrieved total_count=0 names=[] http_status=200.0 nbytes=55.0 path=<repo>\results\md\raw\github_repos_GSE297234.json error=nan note=nan n_samples=nan paywall_interstitial=nan
- source=GitHub repo search mesenchymal_drift url=https://api.github.com/search/repositories?q=mesenchymal+drift retrieval_date=2026-09-19 found=retrieved total_count=0 names=[] http_status=200.0 nbytes=55.0 path=<repo>\results\md\raw\github_repos_mesenchymal_drift.json error=nan note=nan n_samples=nan paywall_interstitial=nan
- source=GitHub repo search lu_mesenchymal url=https://api.github.com/search/repositories?q=mesenchymal+drift+Lu retrieval_date=2026-09-19 found=retrieved total_count=0 names=[] http_status=200.0 nbytes=55.0 path=<repo>\results\md\raw\github_repos_lu_mesenchymal.json error=nan note=nan n_samples=nan paywall_interstitial=nan
- source=GitHub repo search altos_md url=https://api.github.com/search/repositories?q=mesenchymal+drift+altos retrieval_date=2026-09-19 found=retrieved total_count=0 names=[] http_status=200.0 nbytes=55.0 path=<repo>\results\md\raw\github_repos_altos_md.json error=nan note=nan n_samples=nan paywall_interstitial=nan
- source=GitHub code search GSE297234 url=https://api.github.com/search/code?q=GSE297234 retrieval_date=2026-09-19 found=NOT RETRIEVED status=401 (code search often requires auth; not substituting) http_status=401.0 nbytes=120.0 path=nan error=nan note=nan n_samples=nan paywall_interstitial=nan
- … 27 more rows

### Stage 2 MD trajectory all-cell (`s2_md_allcell.csv`)
- missing or empty

### Stage 2 MD trajectory cluster (`s2_md_cluster.csv`)
- missing or empty

### Stage 2 four-cell comparison (`s2_four_cell.csv`)
- missing or empty

### Stage 2 MD vs age ρ (`s2_md_vs_age_rho.csv`)
- missing or empty

### Stage 2 composition (`s2_composition.csv`)
- missing or empty

### Stage 2 random-direction null (`s2_null.csv`)
- missing or empty

### gate.json
- {"open": false, "gene_list_recovered": true, "computation_recovered": false, "gene_list_source": "https://ars.els-cdn.com/content/image/1-s2.0-S0092867425008530-mmc3.xlsx", "computation_source": null, "n_MD_genes": 205, "n_TGFB_genes": 54, "mmc3_url": "https://ars.els-cdn.com/content/image/1-s2.0-S0092867425008530-mmc3.xlsx", "missing": ["MD score computation (per cell vs per sample; mean / ssGSEA / AddModuleScore / other). Open supplement mmc3.xlsx lists genes under column 'MD score' but contains no formula. GEO Sample_data_processing is Cell Ranger v7.1.0 + Seurat v5.0.0 sctransform v2, and does not mention an MD score. Would be in STAR Methods of the paywalled article (https://www.cell.com/cell/fulltext/S0092-8674(25)00853-0 , HTTP 403 this session) or in an author code repository (GitHub repo search total_count=0)."], "reason": "MD gene list recovered from open supplement mmc3.xlsx (n=205 genes, sheet MD_signatures, column 'MD score'; https://ars.els-cdn.com/content/image/1-s2.0-S0092867425008530-mmc3.xlsx). Computation of the score from those genes is not in any public source retrieved this session. STOP after Stage 1.", "notes": {"citing_methods_are_not_the_authors_methods": true, "do_not_approximate_addmodulescore": true, "geo_data_processing_mentions_md_score": false, "mmc3_sheets": {"MD_signatures": {"MD score": 205, "TGFB score": 54}, "Aging_signatures": {"Age up": 1533, "Age down": 2007}, "Fibroblast_subtype_signatures": {"PI16_univ": 778, "LRRC15_myo": 472, "COL3A1_myo": 278}, "Reprog_cell_state_signatures": {"Fibroblast": 200, "PartialReprog": 200, "EarlyPluripotency": 200, "Pluripotency": 200, "NonReprog": 200}}}, "retrieval_date": "2026-09-19"}

### s1_summary.json
- {"n_sources": 67, "n_ok": 38, "n_fail": 4, "geo_accessions": ["GSE297233", "GSE297234"], "pubmed_pmc": null, "gate_open": false, "gate_reason": "MD gene list recovered from open supplement mmc3.xlsx (n=205 genes, sheet MD_signatures, column 'MD score'; https://ars.els-cdn.com/content/image/1-s2.0-S0092867425008530-mmc3.xlsx). Computation of the score from those genes is not in any public source retrieved this session. STOP after Stage 1.", "gene_list_recovered": true, "n_MD_genes": 205, "computation_recovered": false}

### s2_summary.json
- {"ran": false, "reason": "MD score gene list and computation cannot be recovered from public sources. STOP after Stage 1.", "missing": ["MD score computation (per cell vs per sample; mean / ssGSEA / AddModuleScore / other). Not in GEO SOFT, PubMed abstract, or Cell landing summary. Would be in STAR Methods of the paywalled article or in an author code repository."]}

### mmc3_MD_signatures.json
- {"url": "https://ars.els-cdn.com/content/image/1-s2.0-S0092867425008530-mmc3.xlsx", "file": "suppl__ars.els-cdn.com__1-s2.0-S0092867425008530-mmc3.xlsx", "sheets": {"MD_signatures": {"MD score": 205, "TGFB score": 54}, "Aging_signatures": {"Age up": 1533, "Age down": 2007}, "Fibroblast_subtype_signatures": {"PI16_univ": 778, "LRRC15_myo": 472, "COL3A1_myo": 278}, "Reprog_cell_state_signatures": {"Fibroblast": 200, "PartialReprog": 200, "EarlyPluripotency": 200, "Pluripotency": 200, "NonReprog": 200}}, "n_MD_score": 205, "n_TGFB_score": 54, "gene_lists_omitted_from_progress": true}

### geo_gds_two_series.json
- [{"accession": "GSE297234", "title": "Prevalent mesenchymal drift in aging and disease is reversed by partial reprogramming [scRNA-seq]", "n_samples": 8, "pubmedids": ["40816266"], "relations": [], "extrelations": [], "bioproject": "PRJNA1263211", "suppfile": "H5, RDS", "pdat": "2025/08/14"}, {"accession": "GSE297233", "title": "Prevalent mesenchymal drift in aging and disease is reversed by partial reprogramming [RNA-seq]", "n_samples": 8, "pubmedids": ["40816266"], "relations": [], "extrelations": [], "bioproject": "PRJNA1263275", "suppfile": "CSV", "pdat": "2025/08/14"}]

## Failures (manifest)

- **stage1:** OSError: [Errno 22] Invalid argument: '<repo>\\results\\md\\raw\\bioproject_BioProject: https:\\www.ncbi.nlm.nih.gov\\bioproject\\PRJNA1263211.html'

## Files

- `src/md_common.py`, `md_stage1.py`, `md_stage2.py`, `md_findings.py`, `md_run.py`, `md_finalize_stage1.py`
- `results/md/`
- `FINDINGS_MD.md`
- `PROGRESS_MD.md`

