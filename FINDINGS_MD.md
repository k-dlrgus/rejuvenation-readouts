# FINDINGS_MD — what the GSE297234 authors measured, then (if recoverable) our instrument

**Status:** Stage 1 only. Gene list recovered (n=205 from open mmc3.xlsx). Computation not recovered. Gate closed. Stage 2 did not run. No reading fired. Seed `20260914`. boot `20260918`. n_perm=200. n_boot=200. n_random=200. Retrieval date `2026-09-19`. Frozen ruler `<repo>\results\fibro\frozen_ruler_ridge_raw.npz` exists=True.

Does not modify any existing `FINDINGS_*.md` or `FALSIFICATION.md`. No GSE325735. Public sources only. Refit nothing. Nothing averaged across donors, clusters, or regimes.

Flag: `results/md/PREREG_STAGE2.flag` exists=True. Gate open=False.

## Sources table

Retrieved 2026-09-19. Each row is one HTTP/API call this session.

| source | url | retrieval_date | found | http_status | error |
|---|---|---|---|---|---|
| GEO GSE297234 SOFT series | https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE297234&targ=self&form=text&view=brief | 2026-09-19 | retrieved. keys=Series_contact_address,Series_contact_city,Series_contact_country,Series_contact_institute,Series_contact_name,Series_contact_state,Series_contact_zip/postal_code,Series_contributor,Series_geo_accession,Series_last_update_date,Series_overall_design,Series_platfor… | NA | NA |
| GEO GSE297234 SOFT series full | https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE297234&targ=self&form=text&view=full | 2026-09-19 | retrieved. keys=Series_contact_address,Series_contact_city,Series_contact_country,Series_contact_institute,Series_contact_name,Series_contact_state,Series_contact_zip/postal_code,Series_contributor,Series_geo_accession,Series_last_update_date,Series_overall_design,Series_platfor… | 200 | NA |
| GEO GSE297234 SOFT GSM | https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE297234&targ=gsm&form=text&view=brief | 2026-09-19 | retrieved n_gsm=8 characteristic_keys=['batch', 'cell line', 'cell type', 'genotype', 'tissue', 'treatment'] sample_supplementary_n=8 | NA | NA |
| GEO GSE297234 HTML | https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE297234 | 2026-09-19 | retrieved | 200 | NA |
| GEO GSE297233 SOFT series | https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE297233&targ=self&form=text&view=brief | 2026-09-19 | retrieved. keys=Series_contact_address,Series_contact_city,Series_contact_country,Series_contact_institute,Series_contact_name,Series_contact_state,Series_contact_zip/postal_code,Series_contributor,Series_geo_accession,Series_last_update_date,Series_overall_design,Series_platfor… | NA | NA |
| GEO GSE297233 SOFT series full | https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE297233&targ=self&form=text&view=full | 2026-09-19 | retrieved. keys=Series_contact_address,Series_contact_city,Series_contact_country,Series_contact_institute,Series_contact_name,Series_contact_state,Series_contact_zip/postal_code,Series_contributor,Series_geo_accession,Series_last_update_date,Series_overall_design,Series_platfor… | 200 | NA |
| GEO GSE297233 SOFT GSM | https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE297233&targ=gsm&form=text&view=brief | 2026-09-19 | retrieved n_gsm=8 characteristic_keys=['batch', 'cell line', 'cell type', 'genotype', 'tissue', 'treatment'] sample_supplementary_n=8 | NA | NA |
| GEO GSE297233 HTML | https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE297233 | 2026-09-19 | retrieved | 200 | NA |
| NCBI GDS esearch by PMID | https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi?db=gds&term=40816266%5BPMID%5D&retmode=json&retmax=50&tool=age_identity_separability&email=none%40example.invalid | 2026-09-19 | retrieved | 200 | NA |
| BioProject lookup PRJNA1263211 | https://www.ncbi.nlm.nih.gov/bioproject/?term=PRJNA1263211 | 2026-09-19 | retrieved | 200 | NA |
| NCBI BioProject esearch GSE297234 | https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi?db=bioproject&term=GSE297234&retmode=json&retmax=20&tool=age_identity_separability&email=none%40example.invalid | 2026-09-19 | retrieved | 200 | NA |
| NCBI SRA esearch GSE297234 | https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi?db=sra&term=GSE297234&retmode=json&retmax=20&tool=age_identity_separability&email=none%40example.invalid | 2026-09-19 | retrieved | 200 | NA |
| PubMed XML | https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=pubmed&id=40816266&rettype=xml&retmode=xml&tool=age_identity_separability&email=none%40example.invalid | 2026-09-19 | retrieved | 200 | NA |
| PubMed parsed fields | https://pubmed.ncbi.nlm.nih.gov/40816266/ | 2026-09-19 | abstract=yes doi=10.1016/j.cell.2025.07.031 pmc=None n_authors=9 journal=Cell 2025;188(21):5895-5911.e17 | NA | NA |
| PubMed HTML | https://pubmed.ncbi.nlm.nih.gov/40816266/ | 2026-09-19 | NOT RETRIEVED None | 203 | NA |
| Cell journal landing page | https://www.cell.com/cell/fulltext/S0092-8674(25)00853-0 | 2026-09-19 | retrieved; landing page has highlights + summary; paywall interstitial not seen in excerpt; STAR Methods listed as a site nav item, not article methods text | 403 | NA |
| DOI resolver | https://doi.org/10.1016/j.cell.2025.07.031 | 2026-09-19 | retrieved status=200 final=https://linkinghub.elsevier.com/retrieve/pii/S0092867425008530 | 200 | NA |
| ScienceDirect landing | https://www.sciencedirect.com/science/article/pii/S0092867425008530 | 2026-09-19 | NOT RETRIEVED None | 403 | NA |
| PMC esearch by PMID | https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi?db=pmc&term=40816266%5Bpmid%5D&retmode=json&retmax=20&tool=age_identity_separability&email=none%40example.invalid | 2026-09-19 | retrieved pmc_ids=[] (empty list = no PMC deposit found by this query) | 200 | NA |
| elink PubMed→PMC | https://eutils.ncbi.nlm.nih.gov/entrez/eutils/elink.fcgi?dbfrom=pubmed&db=pmc&id=40816266&retmode=json&tool=age_identity_separability&email=none%40example.invalid | 2026-09-19 | retrieved | 200 | NA |
| Europe PMC core JSON | https://www.ebi.ac.uk/europepmc/webservices/rest/search?query=EXT_ID:40816266&resultType=core&format=json | 2026-09-19 | retrieved isOpenAccess='N' pmcid=None hasTextMinedTerms='Y' fullTextUrlList_n=1 license=None | 200 | NA |
| Europe PMC fullTextUrl (subscription required doi) | https://doi.org/10.1016/j.cell.2025.07.031 | 2026-09-19 | retrieved | 200 | NA |
| Europe PMC HAS_PREPRINT AND this PMID | https://www.ebi.ac.uk/europepmc/webservices/rest/search?query=HAS_PREPRINT:Y%20AND%20EXT_ID:40816266&format=json | 2026-09-19 | retrieved | 200 | NA |
| Unpaywall | https://api.unpaywall.org/v2/10.1016/j.cell.2025.07.031?email=none@example.invalid | 2026-09-19 | retrieved oa_status='closed' is_oa=False has_repository_copy=False best_oa=None | 200 | NA |
| OpenAlex | https://api.openalex.org/works/https://doi.org/10.1016/j.cell.2025.07.031 | 2026-09-19 | retrieved is_oa=False oa_status='closed' oa_url=None | 200 | NA |
| Crossref | https://api.crossref.org/works/10.1016/j.cell.2025.07.031 | 2026-09-19 | retrieved n_license=7 relation_keys=[] has_abstract=False | 200 | NA |
| preprint search biorxiv_title | https://www.biorxiv.org/search/prevalent%20mesenchymal%20drift%20in%20aging%20jcode%3Abiorxiv%20numresults%3A25%20sort%3Arelevance-rank | 2026-09-19 | retrieved status=403 nbytes=6187 results_token=None empty_marker=False | 403 | NA |
| preprint search biorxiv_author_lu | https://www.biorxiv.org/search/author%3ALu%20mesenchymal%20drift%20jcode%3Abiorxiv%20numresults%3A25 | 2026-09-19 | retrieved status=403 nbytes=6039 results_token=None empty_marker=False | 403 | NA |
| preprint search biorxiv_phrase | https://www.biorxiv.org/search/%22mesenchymal%20drift%22%20jcode%3Abiorxiv%20numresults%3A25%20sort%3Arelevance-rank | 2026-09-19 | retrieved status=403 nbytes=6108 results_token=None empty_marker=False | 403 | NA |
| preprint search medrxiv_phrase | https://www.medrxiv.org/search/%22mesenchymal%20drift%22%20jcode%3Amedrxiv%20numresults%3A25 | 2026-09-19 | retrieved status=403 nbytes=5994 results_token=None empty_marker=False | 403 | NA |
| preprint search researchsquare_title | https://www.researchsquare.com/browse?q=Prevalent%20mesenchymal%20drift%20in%20aging | 2026-09-19 | retrieved status=200 nbytes=234130 results_token=None empty_marker=False | 200 | NA |
| preprint search osf_preprints | https://osf.io/preprints/discover?q=Prevalent%20mesenchymal%20drift | 2026-09-19 | retrieved status=200 nbytes=4207 results_token=None empty_marker=False | 200 | NA |
| preprint search ssrn | https://www.ssrn.com/index.cfm/en/search/?q=Prevalent%20mesenchymal%20drift | 2026-09-19 | retrieved status=403 nbytes=5989 results_token=None empty_marker=False | 403 | NA |
| Europe PMC preprint title search | https://www.ebi.ac.uk/europepmc/webservices/rest/search?query=TITLE:%22Prevalent%20mesenchymal%20drift%22%20AND%20(SRC:PPR%20OR%20PUB_TYPE:preprint)&format=json | 2026-09-19 | retrieved hitCount=0 | 200 | NA |
| Europe PMC preprint author Lu + mesenchymal drift | https://www.ebi.ac.uk/europepmc/webservices/rest/search?query=%22mesenchymal%20drift%22%20AND%20AUTHOR:Lu%20AND%20(SRC:PPR%20OR%20PUB_TYPE:preprint)&format=json | 2026-09-19 | retrieved hitCount=0 | 200 | NA |
| GitHub repo search GSE297234 | https://api.github.com/search/repositories?q=GSE297234 | 2026-09-19 | retrieved total_count=0 names=[] | 200 | NA |
| GitHub repo search mesenchymal_drift | https://api.github.com/search/repositories?q=mesenchymal+drift | 2026-09-19 | retrieved total_count=0 names=[] | 200 | NA |
| GitHub repo search lu_mesenchymal | https://api.github.com/search/repositories?q=mesenchymal+drift+Lu | 2026-09-19 | retrieved total_count=0 names=[] | 200 | NA |
| GitHub repo search altos_md | https://api.github.com/search/repositories?q=mesenchymal+drift+altos | 2026-09-19 | retrieved total_count=0 names=[] | 200 | NA |
| GitHub code search GSE297234 | https://api.github.com/search/code?q=GSE297234 | 2026-09-19 | NOT RETRIEVED status=401 (code search often requires auth; not substituting) | 401 | NA |
| Zenodo search | https://zenodo.org/api/records?q=%22mesenchymal%20drift%22%20AND%20(Lu%20OR%20GSE297234)&size=20 | 2026-09-19 | retrieved total=1 titles=['The Broken Biological Arrow: Disease, Cellular Reprogramming, and the Ontology of Dynamic State-Space Medicine'] | 200 | NA |
| figshare search POST | https://api.figshare.com/v2/articles/search | 2026-09-19 | retrieved n=0 titles=[] | 200 | NA |
| Code Ocean explore | https://codeocean.com/explore?query=mesenchymal%20drift | 2026-09-19 | NOT RETRIEVED None | 403 | NA |
| open-supplement probe | https://www.cell.com/cell/supplemental/S0092-8674(25)00853-0 | 2026-09-19 | status=403 nbytes=5883 kind=unknown ctype='text/html; charset=UTF-8' | 403 | NA |
| open-supplement probe | https://www.cell.com/cms/10.1016/j.cell.2025.07.031/attachment/mmc1.pdf | 2026-09-19 | status=403 nbytes=5937 kind=unknown ctype='text/html; charset=UTF-8' | 403 | NA |
| open-supplement probe | https://www.cell.com/cell/fulltext/S0092-8674(25)00853-0#supplementaryMaterial | 2026-09-19 | status=403 nbytes=5871 kind=unknown ctype='text/html; charset=UTF-8' | 403 | NA |
| open-supplement probe | https://ars.els-cdn.com/content/image/1-s2.0-S0092867425008530-mmc1.pdf | 2026-09-19 | status=404 nbytes=189 kind=unknown ctype='text/xml' | 404 | NA |
| open-supplement probe | https://ars.els-cdn.com/content/image/1-s2.0-S0092867425008530-mmc1.xlsx | 2026-09-19 | status=200 nbytes=20033 kind=xlsx/zip ctype='application/excel' | 200 | NA |
| open-supplement probe | https://ars.els-cdn.com/content/image/1-s2.0-S0092867425008530-mmc2.xlsx | 2026-09-19 | status=200 nbytes=166098 kind=xlsx/zip ctype='application/excel' | 200 | NA |
| open-supplement probe | https://ars.els-cdn.com/content/image/1-s2.0-S0092867425008530-mmc3.xlsx | 2026-09-19 | status=200 nbytes=84749 kind=xlsx/zip ctype='application/excel' | 200 | NA |
| open-supplement probe | https://ars.els-cdn.com/content/image/1-s2.0-S0092867425008530-mmc4.xlsx | 2026-09-19 | status=404 nbytes=189 kind=unknown ctype='text/xml' | 404 | NA |
| open-supplement probe | https://ars.els-cdn.com/content/image/1-s2.0-S0092867425008530-mmc5.xlsx | 2026-09-19 | status=404 nbytes=189 kind=unknown ctype='text/xml' | 404 | NA |
| open-supplement probe | https://ars.els-cdn.com/content/image/1-s2.0-S0092867425008530-mmc2.pdf | 2026-09-19 | status=404 nbytes=189 kind=unknown ctype='text/xml' | 404 | NA |
| open-supplement probe | https://ars.els-cdn.com/content/image/1-s2.0-S0092867425008530-mmc3.pdf | 2026-09-19 | status=404 nbytes=189 kind=unknown ctype='text/xml' | 404 | NA |
| open-supplement probe | https://ars.els-cdn.com/content/image/1-s2.0-S0092867425008530-mmc6.xlsx | 2026-09-19 | status=404 nbytes=189 kind=unknown ctype='text/xml' | 404 | NA |
| open-supplement probe | https://ars.els-cdn.com/content/image/1-s2.0-S0092867425008530-mmc7.xlsx | 2026-09-19 | status=404 nbytes=189 kind=unknown ctype='text/xml' | 404 | NA |
| open-supplement probe | https://www.cell.com/cms/attachment/stateless/pii/S0092-8674(25)00853-0/mmc1.pdf | 2026-09-19 | status=403 nbytes=5985 kind=unknown ctype='text/html; charset=UTF-8' | 403 | NA |
| extra probe | https://api.biorxiv.org/pubs/10.1016/j.cell.2025.07.031 | 2026-09-19 | status=200 nbytes=74 excerpt={"messages":[{"status":"Server not recognized"}], "collection":[]} | 200 | NA |
| extra probe | https://api.crossref.org/works?query.bibliographic=Prevalent%20mesenchymal%20drift%20Lu&filter=type:posted-content&rows=20 | 2026-09-19 | status=200 nbytes=83305 excerpt={"status":"ok","message-type":"work-list","message-version":"1.0.0","message":{"facets":{},"total-results":65390,"items":[{"indexed":{"date-parts":[[2023,8,19]],"date-time":"2023-0 | 200 | NA |
| extra probe | https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi?db=gds&term=GSE297234[ACCN]&retmode=json | 2026-09-19 | status=200 nbytes=392 excerpt={"header":{"type":"esearch","version":"0.3"},"esearchresult":{"count":"10","retmax":"10","retstart":"0","idlist":["200297234","100024676","308986593","308986592","308986591","30898 | 200 | NA |
| extra probe | https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi?db=gds&term=Prevalent%20mesenchymal%20drift&retmode=json&retmax=20 | 2026-09-19 | status=200 nbytes=540 excerpt={"header":{"type":"esearch","version":"0.3"},"esearchresult":{"count":"2","retmax":"2","retstart":"0","idlist":["200297234","200297233"],"translationset":[],"translationstack":[{"t | 200 | NA |
| extra probe | https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi?db=bioproject&term=PRJNA1263211&retmode=json | 2026-09-19 | status=200 nbytes=303 excerpt={"header":{"type":"esearch","version":"0.3"},"esearchresult":{"count":"1","retmax":"1","retstart":"0","idlist":["1263211"],"translationset":[],"translationstack":[{"term":"PRJNA126 | 200 | NA |
| extra probe | https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi?db=sra&term=PRJNA1263211&retmode=json&retmax=5 | 2026-09-19 | status=200 nbytes=348 excerpt={"header":{"type":"esearch","version":"0.3"},"esearchresult":{"count":"8","retmax":"5","retstart":"0","idlist":["38661437","38661436","38661435","38661434","38661433"],"translation | 200 | NA |
| extra probe | https://ars.els-cdn.com/content/image/1-s2.0-S0092867425008530-mmc4.docx | 2026-09-19 | status=404 nbytes=189 excerpt=<?xml version="1.0" encoding="UTF-8" standalone="yes"?><service-error><status><statusCode>NOT FOUND</statusCode><statusText>Attachment/Metadata missing</statusText></status></servi | 404 | NA |
| extra probe | https://ars.els-cdn.com/content/image/1-s2.0-S0092867425008530-mmc1.docx | 2026-09-19 | status=404 nbytes=189 excerpt=<?xml version="1.0" encoding="UTF-8" standalone="yes"?><service-error><status><statusCode>NOT FOUND</statusCode><statusText>Attachment/Metadata missing</statusText></status></servi | 404 | NA |
| extra probe | https://ars.els-cdn.com/content/image/1-s2.0-S0092867425008530-si1.pdf | 2026-09-19 | status=404 nbytes=189 excerpt=<?xml version="1.0" encoding="UTF-8" standalone="yes"?><service-error><status><statusCode>NOT FOUND</statusCode><statusText>Attachment/Metadata missing</statusText></status></servi | 404 | NA |
| extra probe | https://www.biorxiv.org/search/text_abstract_title%3A%22mesenchymal%20drift%22%20text_abstract_title_flags%3Amatch-phrase%20jcode%3Abiorxiv%20numresults%3A10%20sort%3Arelevance-ra… | 2026-09-19 | status=403 nbytes=6238 excerpt=<!DOCTYPE html><html lang="en-US"><head><title>Just a moment...</title><meta http-equiv="Content-Type" content="text/html; charset=UTF-8"><meta http-equiv="X-UA-Compatible" content | 403 | NA |

## What the authors measured, as far as public sources show

Every sentence below is tied to a URL retrieved this session, or is marked not determinable.

### What is the mesenchymal-drift (MD) score?

- **Genes:** n=205 symbols in open supplement mmc3.xlsx sheet `MD_signatures` column `MD score` ([mmc3.xlsx](https://ars.els-cdn.com/content/image/1-s2.0-S0092867425008530-mmc3.xlsx)). The xlsx does not state weights or a formula. Direction of the gene set is the column name `MD score`; the abstract describes MD as upregulation of mesenchymal genes ([PubMed 40816266](https://pubmed.ncbi.nlm.nih.gov/40816266/); [GEO GSE297234](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE297234)).
- MD score genes, verbatim from [mmc3.xlsx](https://ars.els-cdn.com/content/image/1-s2.0-S0092867425008530-mmc3.xlsx): ABI3BP, ACTA2, ADAM12, ANPEP, APLP1, AREG, BASP1, BDNF, BGN, BMP1, CADM1, CALD1, CALU, CAP2, CAPG, CCN1, CCN2, CD44, CD59, CDH11, CDH2, CDH6, COL11A1, COL12A1, COL16A1, COL1A1, COL1A2, COL3A1, COL4A1, COL4A2, COL5A1, COL5A2, COL5A3, COL6A2, COL6A3, COL7A1, COL8A2, COLGALT1, COMP, COPA, CRLF1, CTHRC1, CXCL1, CXCL12, CXCL6, CXCL8, DAB2, DCN, DKK1, DPYSL3, DST, ECM1, ECM2, EDIL3, EFEMP2, ELN, EMP3, ENO2, FAP, FAS, FBLN1, FBLN2, FBLN5, FBN1, FBN2, FERMT2, FGF2, FLNA, FMOD, FN1, FOXC2, FSTL1, FSTL3, FUCA1, FZD8, GADD45A, GADD45B, GAS1, GEM, GJA1, GLIPR1, GPC1, GPX7, GREM1, HTRA1, ID2, IGFBP2, IGFBP3, IGFBP4, IL15, IL32, IL6, INHBA, ITGA2, ITGA5, ITGAV, ITGB1, ITGB3, ITGB5, JUN, LAMA1, LAMA2, LAMA3, LAMC1, LAMC2, LGALS1, LOX, LOXL1, LOXL2, LRP1, LRRC15, LUM, MAGEE1, MATN2, MATN3, MCM7, MEST, MFAP5, MGP, MMP1, MMP14, MMP2, MMP3, MSX1, MXRA5, MYL9, MYLK, NID2, NNMT, NOTCH2, NT5E, NTM, OXTR, P3H1, PCOLCE, PCOLCE2, PDGFRB, PDLIM4, PFN2, PLAUR, PLOD1, PLOD2, PLOD3, PMEPA1, PMP22, POSTN, PPIB, PRRX1, PRSS2, PTHLH, PTX3, PVR, QSOX1, RGS4, RHOB, SAT1, SCG2, SDC1, SDC4, SERPINE1, SERPINE2, SERPINH1, SFRP1, SFRP4, SGCB, SGCD, SGCG, SLC6A8, SLIT2, SLIT3, SNAI1, SNAI2, SNTB1, SPARC, SPOCK1, SPP1, TAGLN, TFPI2, TGFB1, TGFBI, TGFBR3, TGM2, THBS1, THBS2, THY1, TIMP1, TIMP3, TNC, TNFAIP3, TNFRSF11B, TNFRSF12A, TPM1, TPM2, TPM4, TWIST1, TWIST2, VCAM1, VCAN, VEGFA, VEGFC, VIM, WIPF1, WNT5A, ZEB1, ZEB2.
- Same sheet, separate column `TGFB score`: n=54 genes. Not the MD list. ACVR1, APC, ARID4B, BCAR3, BMP2, BMPR1A, BMPR2, CDH1, CDK9, CDKN1C, CTNNB1, ENG, FKBP1A, FNTA, FURIN, HDAC1, HIPK2, ID1, ID2, ID3, IFNGR2, JUNB, KLF10, LEFTY2, LTBP2, MAP3K7, NCOR2, NOG, PMEPA1, PPM1A, PPP1CA, PPP1R15A, RAB31, RHOA, SERPINE1, SKI, SKIL, SLC20A1, SMAD1, SMAD3, SMAD6, SMAD7, SMURF1, SMURF2, SPTBN1, TGFB1, TGFBR1, TGIF1, THBS1, TJP1, TRIM33, UBE2D3, WWTR1, XIAP.
- Other sheets in the same file ([mmc3.xlsx](https://ars.els-cdn.com/content/image/1-s2.0-S0092867425008530-mmc3.xlsx)): `{'MD_signatures': {'MD score': 205, 'TGFB score': 54}, 'Aging_signatures': {'Age up': 1533, 'Age down': 2007}, 'Fibroblast_subtype_signatures': {'PI16_univ': 778, 'LRRC15_myo': 472, 'COL3A1_myo': 278}, 'Reprog_cell_state_signatures': {'Fibroblast': 200, 'PartialReprog': 200, 'EarlyPluripotency': 200, 'Pluripotency': 200, 'NonReprog': 200}}`. Those are labelled gene lists (aging / fibroblast subtype / reprogramming cell state). They are not a scoring formula.

- **Computation (how those genes become a number, per cell vs per sample):** **not determinable from public sources.** mmc3.xlsx ([url](https://ars.els-cdn.com/content/image/1-s2.0-S0092867425008530-mmc3.xlsx)) is a gene list only. GEO `Sample_data_processing` (same text on each GSE297234 GSM) is: `Samples were demultiplexed and fastq files were generated from the mkfastq pipeline using Cell Ranger v7.1.0. Fastq files were aligned to the human genome GRCh38 using the count pipeline with intron mode included of Cell Ranger v7.1.0. Count matrices were loaded to create Seurat objects using Seurat v5.0.0. Cells with more than 500 genes and genes detected in at least 5 cells were kept for further processing. Separate samples were merged and normalized using sctransform v2, with regression for percent mitochondrial gene.` ([GEO GSE297234 SOFT GSM](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE297234&targ=gsm&form=text&view=brief)). That describes Cell Ranger v7.1.0 and Seurat v5.0.0 sctransform v2. It does not mention MD, AddModuleScore, ssGSEA, GSVA, or a module score. Cell landing / STAR Methods: HTTP 403 ([Cell](https://www.cell.com/cell/fulltext/S0092-8674(25)00853-0)). PMC deposit: none ([PMC esearch](https://www.ncbi.nlm.nih.gov/pmc/?term=40816266); Europe PMC isOpenAccess=N, pmcid=None, [Europe PMC](https://europepmc.org/article/med/40816266)). GitHub repositories named for this dataset or 'mesenchymal drift': total_count=0. Approximating AddModuleScore or Hallmark EMT is disallowed. The source that would supply the computation is the paywalled STAR Methods of Lu et al., *Cell* https://doi.org/10.1016/j.cell.2025.07.031.

### What did they use as the rejuvenation readout?

- **not determinable from public sources** as a named transcriptomic clock vs MD itself vs the `Aging_signatures` sheet in mmc3.xlsx. The public abstract says suppression of "key MD transcription factors leads to epigenetic rejuvenation" and that partial reprogramming can "markedly reduce MD before dedifferentiation and gain of pluripotency, rejuvenating the aging transcriptome at the cellular and tissue levels" ([PubMed 40816266](https://pubmed.ncbi.nlm.nih.gov/40816266/); [GEO GSE297234](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE297234)). mmc3.xlsx has a sheet `Aging_signatures` with columns `Age up` / `Age down` ([mmc3.xlsx](https://ars.els-cdn.com/content/image/1-s2.0-S0092867425008530-mmc3.xlsx)). Whether that sheet, the MD list, a clock, or something else carried the published rejuvenation claim is not stated in GEO, PubMed, or the xlsx headers. Which clock, fit on what, is **not determinable from public sources**.

### Was the rejuvenation claim made on individual cells, cell-type-stratified data, or whole-sample averages?

- **not determinable from public sources** for the reported claim. The GSE297234 record is 10x Genomics scRNA-seq ([GEO overall design](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE297234): `10x Genomics scRNA-seq of human fibroblasts from a young (GM23815, 22yr) and aged donor (GM00731, 96yr) treated with Sendai virus OSKM for up to 10 days.`), so the authors had per-cell data. GEO processing built Seurat objects (see `Sample_data_processing` above). The abstract claims effects "at the cellular and tissue levels" ([PubMed 40816266](https://pubmed.ncbi.nlm.nih.gov/40816266/)). mmc3.xlsx sheet `Reprog_cell_state_signatures` has columns `Fibroblast`, `PartialReprog`, `EarlyPluripotency`, `Pluripotency`, `NonReprog` ([mmc3.xlsx](https://ars.els-cdn.com/content/image/1-s2.0-S0092867425008530-mmc3.xlsx)), which is a cell-state gene-list table, not a statement that the MD reduction was reported per state rather than as a whole-sample average. Sibling bulk series GSE297233 is a different experiment ([GEO GSE297233](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE297233); overall design: RNA-seq of dox OSK and OCT4YR-SK in GM00731 at 4 days). GEO GDS esummary `relations` is empty for both series; there is no GEO SuperSeries. They share PMID 40816266.

### Which timepoints did they compare, and did they pre-specify them?

- Sample titles retrieved from GEO SOFT GSM ([GSE297234](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE297234)): `GM00731, Sendai OSKM, day 0; GM00731, Sendai OSKM, day 3; GM00731, Sendai OSKM, day 7; GM00731, Sendai OSKM, day 10; GM23815, Sendai OSKM, day 0; GM23815, Sendai OSKM, day 3; GM23815, Sendai OSKM, day 7; GM23815, Sendai OSKM, day 10`. Days present in those titles: 0, 3, 7, 10. Overall design says "treated with Sendai virus OSKM for up to 10 days" ([GEO](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE297234)). Treatment protocol (same on each GSM): CytoTune-iPS 2.0 Sendai, MOI 5:5:3 for SeV-KLF4-OCT4-SOX2 / SeV-MYC / SeV-KLF4; day 7 replating onto vitronectin; day 8 Essential 8 Medium. Whether those four days were pre-specified as the MD contrast, and which pair is the reported comparison, is **not determinable from public sources**.

### What did they claim about identity / dedifferentiation, and how was it measured?

- Claim in the public abstract: Yamanaka-factor partial reprogramming can "markedly reduce MD **before** dedifferentiation and gain of pluripotency" ([PubMed 40816266](https://pubmed.ncbi.nlm.nih.gov/40816266/); [GEO GSE297234](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE297234); Cell landing HTTP 403). mmc3.xlsx sheet `Reprog_cell_state_signatures` lists 200 genes each for `Fibroblast`, `PartialReprog`, `EarlyPluripotency`, `Pluripotency`, `NonReprog` ([mmc3.xlsx](https://ars.els-cdn.com/content/image/1-s2.0-S0092867425008530-mmc3.xlsx)). How those lists were scored, and whether they are the published identity/dedifferentiation readout, is **not determinable from public sources**.

GEO `Series_summary` (retrieved this session; SOFT on disk had latin-1-misdecoded UTF-8 quotation marks, repaired to U+201C/U+201D for this quote. PubMed abstract below is the un-repaired public abstract):

> The loss of cellular and tissue identity is a hallmark of aging and numerous diseases, but the underlying mechanisms are not well understood. Our analysis of gene expression data from over 40 human tissues and 20 diseases reveals a pervasive upregulation of mesenchymal genes across multiple cell types, along with an altered composition of stromal cell populations, denoting a “mesenchymal drift” (MD). Increased MD correlates with disease progression, reduced patient survival, and an elevated mortality risk, whereas suppression of key MD transcription factors leads to epigenetic rejuvenation. Notably, Yamanaka factor-induced partial reprogramming can markedly reduce MD before dedifferentiation and gain of pluripotency, rejuvenating the aging transcriptome at the cellular and tissue levels. These findings provide mechanistic insight into the underlying beneficial effects of partial reprogramming and offer a framework for developing interventions to reverse age-related diseases using the partial reprogramming approach.

Source: https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE297234

PubMed abstract (verbatim, retrieved this session):

> The loss of cellular and tissue identity is a hallmark of aging and numerous diseases, but the underlying mechanisms are not well understood. Our analysis of gene expression data from over 40 human tissues and 20 diseases reveals a pervasive upregulation of mesenchymal genes across multiple cell types, along with an altered composition of stromal cell populations, denoting a "mesenchymal drift" (MD). Increased MD correlates with disease progression, reduced patient survival, and an elevated mortality risk, whereas suppression of key MD transcription factors leads to epigenetic rejuvenation. Notably, Yamanaka factor-induced partial reprogramming can markedly reduce MD before dedifferentiation and gain of pluripotency, rejuvenating the aging transcriptome at the cellular and tissue levels. These findings provide mechanistic insight into the underlying beneficial effects of partial reprogramming and offer a framework for developing interventions to reverse age-related diseases using the partial reprogramming approach.

Source: https://pubmed.ncbi.nlm.nih.gov/40816266/


## Stage 1 gate outcome

- open=False
- gene_list_recovered=True
- computation_recovered=False
- gene_list_source=https://ars.els-cdn.com/content/image/1-s2.0-S0092867425008530-mmc3.xlsx
- computation_source=None
- reason: MD gene list recovered from open supplement mmc3.xlsx (n=205 genes, sheet MD_signatures, column 'MD score'; https://ars.els-cdn.com/content/image/1-s2.0-S0092867425008530-mmc3.xlsx). Computation of the score from those genes is not in any public source retrieved this session. STOP after Stage 1.

Missing (verbatim from gate.json):

- MD score computation (per cell vs per sample; mean / ssGSEA / AddModuleScore / other). Open supplement mmc3.xlsx lists genes under column 'MD score' but contains no formula. GEO Sample_data_processing is Cell Ranger v7.1.0 + Seurat v5.0.0 sctransform v2, and does not mention an MD score. Would be in STAR Methods of the paywalled article (https://www.cell.com/cell/fulltext/S0092-8674(25)00853-0 , HTTP 403 this session) or in an author code repository (GitHub repo search total_count=0).

Gate notes (not used as a substitute MD score): do_not_approximate_AddModuleScore=True GEO_data_processing_mentions_MD_score=False mmc3_sheets={'MD_signatures': {'MD score': 205, 'TGFB score': 54}, 'Aging_signatures': {'Age up': 1533, 'Age down': 2007}, 'Fibroblast_subtype_signatures': {'PI16_univ': 778, 'LRRC15_myo': 472, 'COL3A1_myo': 278}, 'Reprog_cell_state_signatures': {'Fibroblast': 200, 'PartialReprog': 200, 'EarlyPluripotency': 200, 'Pluripotency': 200, 'NonReprog': 200}}.

A later paper using Seurat `AddModuleScore` on an MD set, or the phrase "Hallmark EMT plus EMT transcription factors", is not Lu et al.'s method and is not implemented.

Paper identity as retrieved:

- Title: Prevalent mesenchymal drift in aging and disease is reversed by partial reprogramming.
- PMID: 40816266 — https://pubmed.ncbi.nlm.nih.gov/40816266/
- DOI: 10.1016/j.cell.2025.07.031 — https://doi.org/10.1016/j.cell.2025.07.031
- Journal: Cell 2025;188(21):5895-5911.e17
- PMCID: None (None = no PMC id in the PubMed XML this session)
- Authors (PubMed XML): ['Lu Jinlong Y', 'Tu William B', 'Li Ronghui', 'Weng Mingxi', 'Sanketi Bhargav D', 'Yuan Baolei', 'Reddy Pradeep', 'Rodriguez Esteban Concepcion', 'Izpisua Belmonte Juan Carlos']
- GEO scRNA-seq: GSE297234 — https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE297234
- GEO bulk sibling: GSE297233 — https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE297233
- GEO pubmed_id field: 40816266
- GEO supplementary files: ['ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE297nnn/GSE297234/suppl/GSE297234_GM00731_SEVOSKM.rds', 'ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE297nnn/GSE297234/suppl/GSE297234_HFIB_COMBINED_SEVOSKM.rds', 'ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE297nnn/GSE297234/suppl/GSE297234_RAW.tar']
- GEO relation: ['BioProject: https://www.ncbi.nlm.nih.gov/bioproject/PRJNA1263211']
- GSE297233 title: Prevalent mesenchymal drift in aging and disease is reversed by partial reprogramming [RNA-seq]
- GSE297233 design: RNA-seq of dox induced OSK and OCT4YR mutant+SK expression in human fibroblasts from an aged donor (GM00731, 96yr male).
- GEO SuperSeries: none. GDS esummary `relations`=[] for GSE297234 and GSE297233; title search returned exactly those two series (uids 200297234, 200297233).
- BioProject scRNA: PRJNA1263211 (GSE297234). BioProject bulk: PRJNA1263275 (GSE297233).
- SRA runs for PRJNA1263211: count=8 (one per GSM).
- GDS esummary: GSE297234 n_samples=8 suppfile=H5, RDS pdat=2025/08/14 bioproject=PRJNA1263211 relations=[]. GSE297233 n_samples=8 suppfile=CSV pdat=2025/08/14 bioproject=PRJNA1263275 relations=[].
- Unpaywall: is_oa=False, oa_status=closed, has_repository_copy=False. OpenAlex oa_status=closed. Crossref relation={}.
- Preprint: Europe PMC TITLE search hitCount=0; HAS_PREPRINT AND PMID hitCount=0; bioRxiv HTML search HTTP 403; bioRxiv pubs API 'Server not recognized'; Research Square `/browse` page returned the generic latest-preprint listing (473,989), not a title hit.
- Code: GitHub repo search total_count=0 for GSE297234 and for 'mesenchymal drift'. figshare n=0. Zenodo hit is an unrelated record.
- Open Excel supplements: mmc1.xlsx (dataset accessions + cell composition), mmc2.xlsx (plasma-protein tables labelled as reproduction of other PMIDs), mmc3.xlsx (MD / aging / fibroblast / reprogramming gene lists). mmc PDF/docx 404 or Cell HTTP 403.

## Stage 2 pre-registration (verbatim, written before any MD score)

Written **before any MD score or frozen-ruler comparison**, 2026-09-19. No threshold in this block is re-tuned after numbers exist. Stage 2 runs only if Stage 1 recovered the authors' MD gene list and computation from public sources (GEO, PubMed, PMC, preprint servers, code repositories, open supplements). Exact or not at all. Refit nothing. Recalibrate nothing. Rescale nothing. Frozen ruler `results/fibro/frozen_ruler_ridge_raw.npz` used as-is. Clustering and cell filters already frozen in `results/fibro3/` used as-is. GSE325735 is out of scope.

On the GSE297234 matrices already used by FINDINGS_FIBRO / FINDINGS_FIBRO3 (Cell Ranger `filtered_feature_bc_matrix.h5`), using the joint-cluster labels and MIN_CELLS=20 already frozen in `results/fibro3/`:

1. **Reproduce their metric.** Compute the MD score as published, per cell and per (donor × cluster × timepoint). Report its trajectory. State plainly whether our reproduction reproduces their reported direction and magnitude. If it does not, say so and stop — a failed reproduction of their metric means Stage 2's comparison is not valid.
2. **Run their comparison both ways.** For the exact timepoint comparison they made, report side by side: their MD score, and our frozen age ruler. Both on whole-sample averages, and both per cluster. Four cells. Never averaged together.
3. **Cross-check the metrics against each other.** Across all donor × cluster × timepoint rows with ≥20 cells, report Spearman ρ between MD score and frozen age score, with permutation null and bootstrap CI. Two age-related metrics on the same cells should agree if they measure the same thing.
4. **Composition test on their metric.** Apply the FINDINGS_FIBRO3.md per-cluster decomposition to *their* score: does the MD change survive when populations are analysed separately, or is it a mix shift? Report with the same 200-direction null where applicable.

Seed `20260914`. boot `20260918`. n_perm=200. n_boot=200. n_random=200. Spearman ρ with permutation null is the reported statistic. Nothing averaged across donors, clusters, or regimes.

**Pre-registered reading (only the outcome that fired):**
- MD moves per cluster, our ruler does not → the two instruments disagree on the same cells. Report as an instrument disagreement, describe both, and do not declare a winner. Name the test that would settle it.
- Neither moves per cluster, both move on whole-sample averages → the published claim rests on composition, as ours did before decomposition. This is a methodological finding about how reprogramming time courses are read, and it must be stated carefully and without overreach: one dataset, two donors, our reproduction of their metric.
- Both move per cluster → their claim holds and our ruler is insensitive to what OSKM does in these cells. Report it as a limitation of our ruler. This is a real possible outcome and must not be argued away.
- MD and our ruler correlate strongly across rows but disagree on the trajectory → report the discrepancy as unexplained.


## Stage 2 tables

Stage 2 did not run. Gate closed.

s2_summary.json: `{"ran": false, "reason": "MD score gene list and computation cannot be recovered from public sources. STOP after Stage 1.", "missing": ["MD score computation (per cell vs per sample; mean / ssGSEA / AddModuleScore / other). Not in GEO SOFT, PubMed abstract, or Cell landing summary. Would be in STAR Methods of the paywalled article or in an author code repository."]}`

## Reading

No Stage 2 reading fired. The pre-registered instrument/composition/limitation bullets are not evaluated because the authors' MD score could not be implemented exactly.

## What this changes in FINDINGS_FIBRO.md / FINDINGS_FIBRO3.md (quote, do not edit)

This file does not edit `FINDINGS_FIBRO.md` or `FINDINGS_FIBRO3.md`.

FINDINGS_FIBRO.md Stage 2 verdict sentence, quoted:

> age score does not decline → the ruler reads a static donor property, not a modifiable state. Step 3 is not supported by this data.

FINDINGS_FIBRO3.md Task 1 reading sentence, quoted:

> Cell-intrinsic. The same cells' clusters individually reverse from d7 to d10. Then reprogramming genuinely pushes cells back along the age axis after d7, and that is a biological claim worth its own test. clusters=[0, 1, 2].

FINDINGS_FIBRO3.md Task 3 aged-donor d0→d7 line, quoted:

> Aged GM00731 d0→d7 (the all-cell T-B cell had p=+0.005): cluster 0 decline=-0.582 p=+0.652 n_ge=130/200 n_d0=4941 n_d7=2719; cluster 1 decline=+2.215 p=+0.144 n_ge=28/200 n_d0=76 n_d7=850; cluster 2 skipped (n_cells<20 at an endpoint).

Because the authors' MD computation is not recoverable from public sources, this file does not adjudicate whether FINDINGS_FIBRO.md / FINDINGS_FIBRO3.md disagree with their *measurement* or only with their *words*. The gene list is public (mmc3.xlsx, n=205); the scoring formula is not. That adjudication was the point of Stage 2 and did not run.

## Limitations

1. Two donors in GSE297234.
2. Our reproduction of their metric, if Stage 2 had run, would be ours, not theirs.
3. Clusters in FINDINGS_FIBRO3.md are not lineage-tracked.
4. The Cell article is paywalled on the journal landing page retrieved this session; we did not circumvent that paywall.
5. A third-party PDF of the journal article appeared in web search (rapamycin.news) and was not used.
6. Citing papers (Frontiers review; SENOMORPHIC bioRxiv) that paraphrase MD were not treated as the authors' methods.
7. Frozen ruler trained on GTEx V10 cultured fibroblasts, public AGE 10-year bins. Refit nothing.
8. Seed `20260914`. boot `20260918`. n_perm=200, n_boot=200, n_random=200.
9. GSE325735 not opened.
10. Where a fact could not be retrieved, this file writes "not determinable from public sources" rather than an inference.

## Open questions

- How is the mmc3.xlsx `MD score` gene list turned into a number (AddModuleScore / ssGSEA / mean / other; per cell vs pseudobulk)? STAR Methods of the paywalled Cell article is the source that would supply it.
- Was the scRNA MD reduction reported per cell state (`Reprog_cell_state_signatures`) or on mixed-population averages?
- Which timepoint pair is the pre-specified MD contrast (d0→d3 vs d0→d7 vs other)?
- Which rejuvenation readout (MD score vs `Aging_signatures` vs a clock vs both) carried the "aging transcriptome" claim?
- A legal OA copy (PMC author manuscript or preprint) that includes STAR Methods would reopen Stage 2 without approximation.

## Files

| path | content |
|---|---|
| `src/md_common.py`, `md_stage1.py`, `md_stage2.py`, `md_findings.py`, `md_run.py`, `md_finalize_stage1.py` | code |
| `results/md/` | sources, gate, mmc3 gene list, GDS two-series, extra probes, raw retrievals |
| `FINDINGS_MD.md` | this file |
| `PROGRESS_MD.md` | resume state |

s1_summary: `{"n_sources": 67, "n_ok": 38, "n_fail": 4, "geo_accessions": ["GSE297233", "GSE297234"], "pubmed_pmc": null, "gate_open": false, "gate_reason": "MD gene list recovered from open supplement mmc3.xlsx (n=205 genes, sheet MD_signatures, column 'MD score'; https://ars.els-cdn.com/content/image/1-s2.0-S0092867425008530-mmc3.xlsx). Computation of the score from those genes is not in any public source retrieved this session. STOP after Stage 1.", "gene_list_recovered": true, "n_MD_genes": 205, "computation_recovered": false}`

