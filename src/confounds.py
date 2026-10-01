"""Phase 1b: confound flags for genes (cell cycle, ribo/mito, sex, depth, batch).

All flags are computed on the same pseudobulk log-CPM matrix used in 1a. Data-driven flags residualize
expression on cell type first (so a gene is not flagged merely because monocytes have more UMIs than T cells).
"""
import sys
import numpy as np
import pandas as pd
from scipy import stats

sys.path.insert(0, str(__import__("pathlib").Path(__file__).resolve().parent))
from gene_lists import CELL_CYCLE_ALL, S_GENES, G2M_GENES, is_ribo, is_mito, SEX_CORE  # noqa
from decomp import design_matrices  # noqa


def residualize_on_celltype(Y, ct):
    D = pd.get_dummies(pd.Categorical(ct), dtype=float).to_numpy()
    Q, _ = np.linalg.qr(D)
    return Y - Q @ (Q.T @ Y)


def corr_with_covariate(Yr, cov_r):
    """Pearson r between each column of Yr and a (residualized) covariate vector."""
    c = cov_r - cov_r.mean()
    Yc = Yr - Yr.mean(0, keepdims=True)
    num = Yc.T @ c
    den = np.sqrt((Yc ** 2).sum(0) * (c ** 2).sum())
    den[den == 0] = np.nan
    return num / den


def within_pool_age_r2(Y, ct, pool, age):
    """unique age variance after adjusting for cell type AND pool: R2(ct+pool+age) - R2(ct+pool)."""
    Y = np.asarray(Y, float)
    Yc = Y - Y.mean(0, keepdims=True)
    ss_tot = (Yc ** 2).sum(0)
    D_ct = pd.get_dummies(pd.Categorical(ct), dtype=float).to_numpy()
    D_pool = pd.get_dummies(pd.Categorical(pool), dtype=float).to_numpy()[:, 1:]  # drop one level (intercept in D_ct)
    age_c = (np.asarray(age, float) - np.mean(age))[:, None]
    X0 = np.hstack([D_ct, D_pool])
    X1 = np.hstack([X0, age_c])
    Q0, _ = np.linalg.qr(X0)
    Q1, _ = np.linalg.qr(X1)
    r2_0 = ((Q0.T @ Yc) ** 2).sum(0) / ss_tot
    r2_1 = ((Q1.T @ Yc) ** 2).sum(0) / ss_tot
    return r2_1 - r2_0


def sex_de(Y, ct, sex):
    """Within-cell-type sex effect: Cohen's d (male - female) on cell-type-residualized expression, and Welch t p."""
    Yr = residualize_on_celltype(Y, ct)
    m = (np.asarray(sex) == "male")
    f = (np.asarray(sex) == "female")
    mu_m, mu_f = Yr[m].mean(0), Yr[f].mean(0)
    sd = np.sqrt((Yr[m].var(0, ddof=1) + Yr[f].var(0, ddof=1)) / 2)
    sd[sd == 0] = np.nan
    d = (mu_m - mu_f) / sd
    t, p = stats.ttest_ind(Yr[m], Yr[f], equal_var=False)
    return d, p


def build_flags(Y, obs, gene_ids, symbols, chrom: pd.Series, log=print):
    """Return DataFrame indexed by gene_id with boolean flag columns and the underlying statistics."""
    ct = obs.celltype.to_numpy()
    F = pd.DataFrame(index=gene_ids)
    F["symbol"] = symbols
    F["chrom"] = chrom.reindex(gene_ids).fillna("NA").to_numpy()

    # --- curated lists ---
    F["flag_cc_list"] = np.isin(symbols, CELL_CYCLE_ALL)
    F["cc_list_source"] = np.where(np.isin(symbols, S_GENES), "Tirosh_S", np.where(np.isin(symbols, G2M_GENES), "Tirosh_G2M",
                                   np.where(F.flag_cc_list, "prolif_extra", "")))
    F["flag_ribo"] = [is_ribo(s) for s in symbols]
    F["flag_mito"] = [is_mito(s) for s in symbols] 
    F["flag_mito"] |= (F.chrom == "MT")
    F["flag_sex_chrom"] = F.chrom.isin(["X", "Y"])
    F["flag_sex_core_list"] = np.isin(symbols, SEX_CORE)

    # --- data-driven, residualized on cell type ---
    Yr = residualize_on_celltype(Y, ct)
    def rcov(name):
        v = obs[name].to_numpy(float)
        return residualize_on_celltype(v[:, None], ct)[:, 0]
    F["r_cc_umi_frac"] = corr_with_covariate(Yr, rcov("cc_umi_frac"))
    F["r_prolif_cell_frac"] = corr_with_covariate(Yr, rcov("prolif_marker_cell_frac"))
    F["r_mean_counts"] = corr_with_covariate(Yr, rcov("mean_counts_per_cell"))
    F["r_mean_genes"] = corr_with_covariate(Yr, rcov("mean_genes_per_cell"))
    F["r_log_n_cells"] = corr_with_covariate(Yr, residualize_on_celltype(np.log(obs.n_cells.to_numpy(float))[:, None], ct)[:, 0])
    d, p = sex_de(Y, ct, obs.sex.to_numpy())
    F["sex_cohen_d"] = d
    F["sex_p"] = p
    # thresholds (fixed a priori; see FINDINGS)
    F["flag_cc_data"] = (F.r_cc_umi_frac.abs() >= 0.30) | (F.r_prolif_cell_frac.abs() >= 0.30)
    F["flag_depth"] = (F.r_mean_counts.abs() >= 0.30) | (F.r_mean_genes.abs() >= 0.30)
    F["flag_sex_data"] = (F.sex_cohen_d.abs() >= 0.50) & (F.sex_p < 1e-6)

    # --- batch: does the age signal survive within-pool? ---
    if "pool_number" in obs:
        F["unique_age_within_pool"] = within_pool_age_r2(Y, ct, obs.pool_number.to_numpy(), obs.age.to_numpy())
    F["flag_cell_cycle"] = F.flag_cc_list | F.flag_cc_data
    F["flag_sex"] = F.flag_sex_chrom | F.flag_sex_core_list | F.flag_sex_data
    F["flag_technical"] = F.flag_ribo | F.flag_mito | F.flag_depth
    return F


def covariate_age_report(obs, log=print):
    """How do the confound covariates themselves relate to age (within cell type)?"""
    rows = []
    for c in ("cc_umi_frac", "prolif_marker_cell_frac", "mean_counts_per_cell", "mean_genes_per_cell", "n_cells"):
        for ct, g in obs.groupby("celltype", observed=True):
            r, p = stats.spearmanr(g.age, g[c])
            rows.append(dict(covariate=c, celltype=ct, spearman_r=r, p=p, n=len(g)))
    rep = pd.DataFrame(rows)
    piv = rep.pivot(index="celltype", columns="covariate", values="spearman_r")
    log("Spearman(age, covariate) within cell type:")
    log(piv.to_string(float_format=lambda x: f"{x:+.3f}"))
    return rep
