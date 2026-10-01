"""Ridge / kernel-ridge / logistic helpers and nested gene selection for P3–P4."""
from __future__ import annotations

import numpy as np
import pandas as pd
import warnings
from sklearn.exceptions import ConvergenceWarning

from sklearn.linear_model import Ridge, LogisticRegression
from sklearn.metrics import accuracy_score, f1_score
from sklearn.preprocessing import LabelEncoder

from brain_phase1_common import (
    ALPHAS, LOGREG_CS, N_INNER, PRIMARY_QUANTILES, PHASE1_DIR,
    decompose_brain, quantiles_to_thresholds, apply_thresholds,
    donor_grouped_inner_splits, r2_mae,
)
from brain_phase1_p2 import FILTER_STEPS, survival_table, add_p0_flags
from confounds import build_flags
from gene_lists import CELL_CYCLE_ALL  # noqa: F401


def _std_apply(Xtr, Xte):
    Xtr = np.asarray(Xtr, float)
    Xte = np.asarray(Xte, float)
    mu = Xtr.mean(0)
    sd = Xtr.std(0)
    sd[sd < 1e-12] = 1.0
    return (Xtr - mu) / sd, (Xte - mu) / sd


def select_ridge_alpha(X, y, groups, rng, alphas=ALPHAS, n_inner=N_INNER):
    X = np.asarray(X, float)
    y = np.asarray(y, float)
    if len(y) < 8 or X.shape[1] == 0:
        return 1.0
    splits = donor_grouped_inner_splits(groups, rng, n_inner=n_inner)
    if not splits:
        return 1.0
    scores = np.zeros(len(alphas))
    n_ok = 0
    for tr, te in splits:
        Xtr, Xte = _std_apply(X[tr], X[te])
        ytr, yte = y[tr], y[te]
        for i, a in enumerate(alphas):
            pred = Ridge(alpha=float(a), fit_intercept=True).fit(Xtr, ytr).predict(Xte)
            scores[i] += -((yte - pred) ** 2).sum()
        n_ok += 1
    if n_ok == 0:
        return 1.0
    return float(alphas[int(np.argmax(scores))])


def ridge_predict(Xtr, ytr, Xte, groups_tr, rng, alphas=ALPHAS):
    if Xtr.shape[1] == 0 or len(ytr) < 5:
        return np.full(len(Xte), np.nan), np.nan
    a = select_ridge_alpha(Xtr, ytr, groups_tr, rng, alphas=alphas)
    Xtr_s, Xte_s = _std_apply(Xtr, Xte)
    pred = Ridge(alpha=a, fit_intercept=True).fit(Xtr_s, ytr).predict(Xte_s)
    return pred, a


def kernel_ridge_predict(Xtr, ytr, Xte, groups_tr, rng, alphas=ALPHAS):
    """Kernel ridge (linear kernel / n_features) — unrestricted high-p ceiling."""
    Xtr = np.asarray(Xtr, float)
    Xte = np.asarray(Xte, float)
    ytr = np.asarray(ytr, float)
    if len(ytr) < 5:
        return np.full(len(Xte), np.nan), np.nan
    Xtr_s, Xte_s = _std_apply(Xtr, Xte)
    p = max(Xtr_s.shape[1], 1)
    Ktr = Xtr_s @ Xtr_s.T / p
    Kte = Xte_s @ Xtr_s.T / p
    splits = donor_grouped_inner_splits(groups_tr, rng, n_inner=N_INNER)
    scores = np.zeros(len(alphas))
    if splits:
        for tr, te in splits:
            mu = ytr[tr].mean()
            S, U = np.linalg.eigh(Ktr[np.ix_(tr, tr)])
            Uy = U.T @ (ytr[tr] - mu)
            K_te = Ktr[np.ix_(te, tr)]
            for i, a in enumerate(alphas):
                dual = U @ (Uy / (S + a))
                pred = mu + K_te @ dual
                scores[i] += -((ytr[te] - pred) ** 2).sum()
        a_best = float(alphas[int(np.argmax(scores))])
    else:
        a_best = 1.0
    mu = ytr.mean()
    S, U = np.linalg.eigh(Ktr)
    dual = U @ ((U.T @ (ytr - mu)) / (S + a_best))
    pred = mu + Kte @ dual
    return pred, a_best


def select_logreg_C(X, y, groups, rng, Cs=LOGREG_CS):
    """Fixed C=1.0. Inner C search is ~12× multinomial lbfgs fits per outer fold and
    was not a FALSIFICATION.md hyperparameter; C is not tuned on the cross-test."""
    return 1.0


def _logreg(C):
    # saga is used because lbfgs on ~2k I-genes × 20 types does not converge in 200
    # iters and dominated runtime; C=1.0 is fixed (not tuned on the cross-test).
    return LogisticRegression(
        C=float(C), solver="saga", max_iter=80, class_weight=None, tol=1e-3,
        random_state=0,
    )


def logreg_predict(Xtr, ytr, Xte, groups_tr, rng):
    if Xtr.shape[1] == 0 or len(np.unique(ytr)) < 2:
        return np.array(["__none__"] * len(Xte)), np.full((len(Xte), 1), np.nan), None, np.nan
    C = select_logreg_C(Xtr, ytr, groups_tr, rng)
    Xtr_s, Xte_s = _std_apply(Xtr, Xte)
    with warnings.catch_warnings():
        warnings.simplefilter("ignore", ConvergenceWarning)
        clf = _logreg(C).fit(Xtr_s, ytr)
    pred = clf.predict(Xte_s)
    proba = clf.predict_proba(Xte_s)
    return pred, proba, clf, C


def load_chrom_series(gene_ids):
    p = PHASE1_DIR / "p2_gene_flags.csv"
    if p.exists():
        F = pd.read_csv(p, index_col=0)
        s = F["chrom"].astype(str)
        return s.reindex(gene_ids).fillna("NA")
    from ensembl_chrom import fetch_chromosomes
    chrom = fetch_chromosomes(gene_ids, verbose=False)
    m = pd.Series(chrom.chrom.astype(str).to_numpy(),
                  index=[str(i).split(".")[0] for i in chrom.gene_id])
    return m.reindex([str(g).split(".")[0] for g in gene_ids]).fillna("NA")


def select_genes_on_train(Y, obs, genes, extra_cols, train_idx, chrom_s=None):
    """Re-derive A/I on training pseudobulks only (declared quantile LEVELS + P2 filters)."""
    obs_tr = obs.iloc[train_idx].copy()
    Ytr = Y[train_idx]
    df = decompose_brain(Ytr, obs_tr, extra_cols=extra_cols,
                         gene_ids=genes.gene_id.to_numpy(), symbols=genes.symbol.to_numpy())
    thresh = quantiles_to_thresholds(df, PRIMARY_QUANTILES)
    lab = apply_thresholds(df, thresh)
    df["gene_class"] = lab
    gid = genes.gene_id.to_numpy()
    sym = genes.symbol.to_numpy()
    if chrom_s is None:
        chrom_s = load_chrom_series(gid)
    chrom_s = pd.Series(chrom_s.to_numpy(), index=gid)
    F = build_flags(Ytr, obs_tr, gid, sym, chrom_s, log=lambda *a, **k: None)
    F = add_p0_flags(F, Ytr, obs_tr, extra_cols, log=lambda *a, **k: None)
    A = df.index[lab == "A"]
    _, removed = survival_table(A, F)
    A_keep = [g for g in A if not bool(removed.loc[g])]
    I_keep = list(df.index[lab == "I"])
    # column indices into Y
    id_to_j = {g: i for i, g in enumerate(gid)}
    A_idx = np.array([id_to_j[g] for g in A_keep], dtype=int)
    I_idx = np.array([id_to_j[g] for g in I_keep], dtype=int)
    return dict(A_idx=A_idx, I_idx=I_idx, n_A=len(A_idx), n_I=len(I_idx),
                thresh=thresh, A_ids=A_keep, I_ids=I_keep)


def gene_index(genes, ids):
    m = {g: i for i, g in enumerate(genes.gene_id.astype(str))}
    return np.array([m[str(g)] for g in ids if str(g) in m], dtype=int)


def age_oof_for_genes(Y, obs, gene_idx, folds, rng, method="ridge"):
    """Within-cell-type OOF age predictions using a (possibly fold-specific) gene set.

    folds: list of (tr, te, *rest). gene_idx is an int array OR a list of arrays aligned to folds
    (nested selection).
    Returns dict with row-level pred, per-type metrics, donor-mean joint pred.
    """
    pred = np.full(len(obs), np.nan)
    alphas = []
    types = sorted(obs.celltype.astype(str).unique())
    ct = obs.celltype.astype(str).to_numpy()
    y = obs.age.to_numpy(float)
    donors = obs.donor.astype(str).to_numpy()
    for fi, fold in enumerate(folds):
        tr, te = fold[0], fold[1]
        gi = gene_idx[fi] if isinstance(gene_idx, (list, tuple)) else gene_idx
        gi = np.asarray(gi, dtype=int)
        for t in types:
            tr_t = tr[ct[tr] == t]
            te_t = te[ct[te] == t]
            if len(tr_t) < 8 or len(te_t) < 1:
                continue
            Xtr = Y[np.ix_(tr_t, gi)] if gi.size else np.zeros((len(tr_t), 0))
            Xte = Y[np.ix_(te_t, gi)] if gi.size else np.zeros((len(te_t), 0))
            fn = ridge_predict if method == "ridge" else kernel_ridge_predict
            p, a = fn(Xtr, y[tr_t], Xte, donors[tr_t], rng)
            pred[te_t] = p
            alphas.append(a)
    per_type = []
    for t in types:
        m = ct == t
        mm = m & np.isfinite(pred)
        met = r2_mae(y[mm], pred[mm])
        per_type.append(dict(celltype=t, n=int(mm.sum()), r2=met["r2"], mae=met["mae"]))
    pt = pd.DataFrame(per_type)
    # donor-level: mean of available type predictions
    tmp = pd.DataFrame({"donor": donors, "y": y, "pred": pred, "Source": obs.Source.astype(str).to_numpy()})
    tmp = tmp[np.isfinite(tmp.pred)]
    don = tmp.groupby("donor").agg(y=("y", "first"), pred=("pred", "mean"), Source=("Source", "first"))
    joint = r2_mae(don.y, don.pred)
    median_r2 = float(pt.r2.median()) if len(pt) else np.nan
    median_mae = float(pt.mae.median()) if len(pt) else np.nan
    return dict(pred=pred, per_type=pt, joint_r2=joint["r2"], joint_mae=joint["mae"],
                median_r2=median_r2, median_mae=median_mae,
                donor_pred=don, alpha_mean=float(np.nanmean(alphas) if alphas else np.nan))


def nested_age_oof(Y, obs, genes, extra_cols, folds, rng, chrom_s=None, method="ridge"):
    gene_idx_list = []
    infos = []
    for fold in folds:
        tr = fold[0]
        info = select_genes_on_train(Y, obs, genes, extra_cols, tr, chrom_s=chrom_s)
        gene_idx_list.append(info["A_idx"])
        infos.append(info)
    out = age_oof_for_genes(Y, obs, gene_idx_list, folds, rng, method=method)
    out["fold_info"] = infos
    out["median_n_A"] = float(np.median([i["n_A"] for i in infos]))
    return out


def identity_oof_for_genes(Y, obs, gene_idx, folds, rng):
    y = obs.celltype.astype(str).to_numpy()
    donors = obs.donor.astype(str).to_numpy()
    pred = np.array([None] * len(obs), dtype=object)
    proba_true = np.full(len(obs), np.nan)
    dist_own = np.full(len(obs), np.nan)
    classes_ref = None
    for fi, fold in enumerate(folds):
        tr, te = fold[0], fold[1]
        gi = gene_idx[fi] if isinstance(gene_idx, (list, tuple)) else gene_idx
        gi = np.asarray(gi, dtype=int)
        Xtr = Y[np.ix_(tr, gi)] if gi.size else np.zeros((len(tr), 0))
        Xte = Y[np.ix_(te, gi)] if gi.size else np.zeros((len(te), 0))
        yhat, proba, clf, C = logreg_predict(Xtr, y[tr], Xte, donors[tr], rng)
        pred[te] = yhat
        if clf is not None:
            cls = list(clf.classes_)
            classes_ref = cls
            # P(true type)
            yte = y[te]
            col = {c: j for j, c in enumerate(cls)}
            for i, lab in enumerate(yte):
                j = col.get(lab)
                proba_true[te[i]] = proba[i, j] if j is not None else np.nan
            # centroid distance on standardized features
            Xtr_s, Xte_s = _std_apply(Xtr, Xte)
            cents = {}
            for c in cls:
                cents[c] = Xtr_s[y[tr] == c].mean(0)
            for i, lab in enumerate(yte):
                if lab not in cents:
                    continue
                dist_own[te[i]] = float(np.linalg.norm(Xte_s[i] - cents[lab]))
    mask = np.array([p is not None for p in pred])
    y_m, p_m = y[mask], pred[mask].astype(str)
    acc = float(accuracy_score(y_m, p_m)) if mask.any() else np.nan
    f1 = float(f1_score(y_m, p_m, average="macro")) if mask.any() else np.nan
    vc = pd.Series(y).value_counts(normalize=True)
    chance = float(max(1.0 / max(obs.celltype.nunique(), 1), vc.max()))
    return dict(pred=pred, proba_true=proba_true, dist_own=dist_own,
                acc=acc, macro_f1=f1, chance=chance,
                mean_proba_true=float(np.nanmean(proba_true)),
                mean_dist_own=float(np.nanmean(dist_own)), n=int(mask.sum()))


def nested_identity_oof(Y, obs, genes, extra_cols, folds, rng, chrom_s=None, which="I"):
    gene_idx_list = []
    for fold in folds:
        info = select_genes_on_train(Y, obs, genes, extra_cols, fold[0], chrom_s=chrom_s)
        gene_idx_list.append(info["I_idx"] if which == "I" else info["A_idx"])
    out = identity_oof_for_genes(Y, obs, gene_idx_list, folds, rng)
    out["median_n_genes"] = float(np.median([len(g) for g in gene_idx_list]))
    return out


def size_matched_I_age(Y, obs, I_idx, n_A, folds, rng, n_draws, method="ridge"):
    """Age model on random subsets of I-genes with size n_A; median R² over draws."""
    I_idx = np.asarray(I_idx, dtype=int)
    if I_idx.size == 0 or n_A <= 0:
        return dict(median_r2=np.nan, median_mae=np.nan, draws=[])
    k = int(min(n_A, I_idx.size))
    rows = []
    for d in range(n_draws):
        sub = rng.choice(I_idx, size=k, replace=False)
        out = age_oof_for_genes(Y, obs, sub, folds, rng, method=method)
        rows.append(dict(draw=d, median_r2=out["median_r2"], joint_r2=out["joint_r2"],
                         median_mae=out["median_mae"]))
    tab = pd.DataFrame(rows)
    return dict(median_r2=float(tab.median_r2.median()),
                joint_r2=float(tab.joint_r2.median()),
                median_mae=float(tab.median_mae.median()),
                draws=tab)


def permute_age_within_site_obs(obs, rng):
    don = obs.groupby("donor").agg(age=("age", "first"), site=("Source", "first"))
    perm = don.age.copy()
    for _, idx in don.groupby("site").groups.items():
        perm.loc[idx] = rng.permutation(don.age.loc[idx].to_numpy())
    out = obs.copy()
    out["age"] = out.donor.map(perm).to_numpy(float)
    return out


def permute_celltype(obs, rng):
    out = obs.copy()
    out["celltype"] = rng.permutation(out.celltype.astype(str).to_numpy())
    return out
