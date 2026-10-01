"""Curated gene lists for confound flagging (Phase 1b).

Cell-cycle lists are the Tirosh et al. 2016 (Science) S and G2/M signatures as
distributed with Seurat/scanpy (`cc.genes.updated.2019` symbols). Aliases for
renamed symbols are included so matching works on current HGNC names.
"""
import re

S_GENES = [
    "MCM5", "PCNA", "TYMS", "FEN1", "MCM2", "MCM4", "RRM1", "UNG", "GINS2", "MCM6", "CDCA7", "DTL",
    "PRIM1", "UHRF1", "MLF1IP", "CENPU", "HELLS", "RFC2", "RPA2", "NASP", "RAD51AP1", "GMNN", "WDR76",
    "SLBP", "CCNE2", "UBR7", "POLD3", "MSH2", "ATAD2", "RAD51", "RRM2", "CDC45", "CDC6", "EXO1", "TIPIN",
    "DSCC1", "BLM", "CASP8AP2", "USP1", "CLSPN", "POLA1", "CHAF1B", "BRIP1", "E2F8",
]
G2M_GENES = [
    "HMGB2", "CDK1", "NUSAP1", "UBE2C", "BIRC5", "TPX2", "TOP2A", "NDC80", "CKS2", "NUF2", "CKS1B", "MKI67",
    "TMPO", "CENPF", "TACC3", "FAM64A", "PIMREG", "SMC4", "CCNB2", "CKAP2L", "CKAP2", "AURKB", "BUB1", "KIF11",
    "ANP32E", "TUBB4B", "GTSE1", "KIF20B", "HJURP", "CDCA3", "HN1", "JPT1", "CDC20", "TTK", "CDC25C", "KIF2C",
    "RANGAP1", "NCAPD2", "DLGAP5", "CDCA2", "CDCA8", "ECT2", "KIF23", "HMMR", "AURKA", "PSRC1", "ANLN", "LBR",
    "CKAP5", "CENPE", "CTCF", "NEK2", "G2E3", "GAS2L3", "CBX5", "CENPA",
]
# Broader proliferation / cell-cycle machinery not in the Tirosh lists (conservative superset used for flagging)
PROLIF_EXTRA = [
    "STMN1", "TUBA1B", "TUBB", "TUBA1C", "H2AFZ", "H2AZ1", "HMGB1", "HMGN2", "HMGB3", "CCNA2", "CCNB1", "CCND1",
    "CCND2", "CCND3", "CCNE1", "CDK2", "CDK4", "CDK6", "CDKN3", "E2F1", "E2F2", "MYBL2", "FOXM1", "PLK1", "BUB1B",
    "KIF4A", "KIF14", "KIF15", "KIF18B", "KIF22", "KIFC1", "CENPM", "CENPW", "CENPN", "CENPK", "CENPH", "MCM3",
    "MCM7", "MCM10", "ORC1", "ORC6", "CDT1", "GINS1", "PBK", "SPC24", "SPC25", "ZWINT", "MAD2L1", "SGO1", "SGO2",
    "ESPL1", "PTTG1", "TK1", "DHFR", "RFC4", "RFC5", "POLE2", "POLA2", "DUT", "H1-0", "H1-2", "H1-3", "H1-4", "H1-5",
    "H2BC12", "H4C3", "HIST1H4C", "HIST1H1B", "HIST1H1E", "HIST1H2BK", "MKI67IP", "NUSAP1", "ASPM", "DEPDC1", "DEPDC1B",
    "SHCBP1", "PRC1", "ARHGAP11A", "CIT", "KPNA2", "LMNB1", "TOP2B", "SMC2", "NCAPG", "NCAPG2", "NCAPH", "PARPBP",
    "RACGAP1", "FANCI", "FANCD2", "CHEK1", "WEE1", "CDC25A", "CDC25B", "TFDP1", "TFDP2", "RB1", "MYC",
]
CELL_CYCLE_ALL = sorted(set(S_GENES) | set(G2M_GENES) | set(PROLIF_EXTRA))

RIBO_RE = re.compile(r"^(RPL\d+[A-Z]?\d*|RPS\d+[A-Z]?\d*|RPLP\d|RPSA|MRPL\d+|MRPS\d+)$")
MITO_RE = re.compile(r"^MT-")

# Y-linked genes commonly expressed in PBMC and X-inactivation transcripts; chromosome-based
# flagging (Ensembl) is the primary sex filter, this list is a fallback / sanity check.
SEX_CORE = [
    "XIST", "TSIX", "JPX", "FTX", "RPS4Y1", "DDX3Y", "UTY", "KDM5D", "EIF1AY", "USP9Y", "ZFY", "TXLNGY", "NLGN4Y",
    "PRKY", "TMSB4Y", "LINC00278", "TTTY14", "TTTY15", "TTTY10", "ZFX", "KDM6A", "KDM5C", "EIF1AX", "DDX3X",
    "RPS4X", "USP9X", "PRKX", "NLGN4X", "TXLNG", "TMSB4X", "SMC1A", "PUDP", "STS", "EIF2S3", "ZRSR2", "CA5B",
]


def is_ribo(sym: str) -> bool:
    return bool(RIBO_RE.match(sym))


def is_mito(sym: str) -> bool:
    return bool(MITO_RE.match(sym))
