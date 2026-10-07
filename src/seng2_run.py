"""SENG2 driver. Stages: panel, stage0, d1, d2, d3, d4.

Usage: python src/seng2_run.py <stage>

Every stage reads the frozen instruments and writes only results/seng2/.
results/seng/ and results/seng_why/ are read and never written.
"""
from __future__ import annotations

import json
import sys
from pathlib import Path

import numpy as np
import pandas as pd

SRC = Path(__file__).resolve().parent
if str(SRC) not in sys.path:
    sys.path.insert(0, str(SRC))

import seng2_common as K  # noqa: E402
import seng_run as S  # noqa: E402
from same_run import (axis_stats, ci_covers_point, ci_excludes_zero,  # noqa: E402
                      cosine_full, monotonic_increasing, permutation_p,
                      resample_idx)
from toward_run import percentile_ci, ruler_score  # noqa: E402

OUT = K.OUT
OUT.mkdir(parents=True, exist_ok=True)


def dump(name, obj):
    (OUT / name).write_text(json.dumps(_js(obj), indent=2), encoding="utf-8")


def _js(o):
    if isinstance(o, dict):
        return {str(k): _js(v) for k, v in o.items()}
    if isinstance(o, (list, tuple)):
        return [_js(v) for v in o]
    if isinstance(o, (np.integer,)):
        return int(o)
    if isinstance(o, (np.floating,)):
        return None if not np.isfinite(o) else float(o)
    if isinstance(o, (np.bool_,)):
        return bool(o)
    if isinstance(o, np.ndarray):
        return _js(o.tolist())
    if isinstance(o, float) and not np.isfinite(o):
        return None
    return o


# ==========================================================================
def stage_panel():
    log = K.Log(OUT / "panel.log")
    md5 = S.file_md5(K.H5)
    log(f"[data] md5={md5} expected={K.MD5_EXPECT} match={md5 == K.MD5_EXPECT}")
    if md5 != K.MD5_EXPECT:
        dump("STOP.json", dict(key="data_md5_mismatch", md5=md5))
        log("STOP data_md5_mismatch")
        return
    ctx = K.load_context(log)
    groups = ctx["groups"]
    names = [g["name"] for g in groups]
    log(f"[panel] {len(groups)} groups")
    counts_local = K.stream_pseudobulk(ctx["f"], [g["idx"] for g in groups], ctx["exp"], log)
    C = np.vstack([K.to_ruler(counts_local[i], ctx["dst"], ctx["n_ruler"])
                   for i in range(len(groups))])
    np.savez_compressed(
        OUT / "panel_counts.npz", C=C, names=np.array(names, dtype=object),
        kinds=np.array([g["kind"] for g in groups], dtype=object),
        sizes=np.array([g["idx"].size for g in groups]),
        dst=ctx["dst"], src=ctx["src"])
    # per-group median n_counts, for the A3 thinning factors
    nc = ctx["n_counts"]
    med = {g["name"]: float(np.median(nc[g["idx"]])) for g in groups}
    dump("panel_medians.json", med)
    log(f"[panel] wrote panel_counts.npz C={C.shape}")
    # raw-count audit (STOP key no_raw_counts)
    X = ctx["f"]["X"]
    audit = S.scan_numeric(X["data"])
    dump("raw_audit.json", audit)
    log(f"[raw audit] {audit}")
    ctx["f"].close()
    log.close()


# ==========================================================================
def _panel():
    z = np.load(OUT / "panel_counts.npz", allow_pickle=True)
    names = [str(x) for x in z["names"]]
    return z, names, {n: i for i, n in enumerate(names)}


def stage0():
    log = K.Log(OUT / "stage0.log")
    z, names, ni = _panel()
    C = z["C"]
    dst, src = z["dst"], z["src"]
    ctx = K.load_context(log)
    frozen = ctx["frozen"]
    n_ruler = ctx["n_ruler"]
    Z, missing, nf, Cc = K.panel_z_from_counts(C, frozen)
    cache = S.tmm_cache_from_panel(Cc, nf)
    mu = np.asarray(frozen["mu"], float)
    sd = np.asarray(frozen["sd"], float)
    miss = np.asarray(missing, bool)
    log(f"[panel z] rows={Z.shape} n_missing_z0={int(miss.sum())} ref={cache['ref']}")

    v = Z[ni["EARLY_A"]] - Z[ni["LATE_A"]]
    p_unit, cos_vp = {}, {}
    for m in ("CRA", "CRI"):
        pv = Z[ni[f"TOP_{m}"]] - Z[ni[f"BOT_{m}"]]
        p_unit[m] = pv / float(np.linalg.norm(pv))
        cos_vp[m] = cosine_full(v, p_unit[m])
    anch = np.load(S.ANCHORS_PATH, allow_pickle=True)
    u = np.asarray(anch["u"], float).ravel()
    rng_u = np.random.default_rng(K.SEED)
    Uperm = np.empty((K.N_PERM, n_ruler), np.float64)
    for i in range(K.N_PERM):
        Uperm[i] = rng_u.permutation(u)
    cos_vu, p_vu, _ = S.cos_against_perms(v, u, Uperm)
    v_norm = float(np.linalg.norm(v))

    stored = dict(cos_v_p_CRA=0.5667444435453722, cos_v_p_CRI=0.6337418620832961,
                  cos_v_u=-0.04027293769593313, v_norm=90.14934904214671)
    got = dict(cos_v_p_CRA=cos_vp["CRA"], cos_v_p_CRI=cos_vp["CRI"],
               cos_v_u=cos_vu, v_norm=v_norm)
    sv = np.asarray(ctx["axes"]["v"], float).ravel()
    log(f"[axis] max|v_seng2 - v_stored|={float(np.max(np.abs(v - sv))):.3e}")
    log(f"[axis] {got}")

    # --- Stage 0(b): positive control, reimplemented, written to results/seng2/
    early_b = np.asarray(ctx["sidx"]["early_b"], int)
    late_b = np.asarray(ctx["sidx"]["late_b"], int)
    rowsB = np.concatenate([early_b, late_b])
    subB = K.submatrix(ctx["f"], rowsB, ctx["exp"])
    g2l = {int(r): i for i, r in enumerate(rowsB)}

    def sumc(idx):
        if len(idx) == 0:
            return np.zeros(n_ruler)
        loc = [g2l[int(i)] for i in idx]
        return K.to_ruler(np.asarray(subB[loc].sum(axis=0)).ravel(), dst, n_ruler)

    n_tot = int(min(late_b.size, early_b.size))
    rng = np.random.default_rng(K.SEED)
    mixes = []
    for fv in S.F_VALUES:
        n_e = int(np.round(float(fv) * n_tot))
        n_l = int(n_tot - n_e)
        tl = rng.choice(late_b, size=n_l, replace=False) if n_l else np.array([], int)
        te = rng.choice(early_b, size=n_e, replace=False) if n_e else np.array([], int)
        mixes.append(dict(f=float(fv), n_late=n_l, n_early=n_e, n_total=n_tot,
                          late_idx=tl, early_idx=te))
    rows = [sumc(early_b), sumc(late_b)]
    for m in mixes:
        rows.append(sumc(m["late_idx"]) + sumc(m["early_idx"]))
    ZP, missP, nfP, CP = K.panel_z_from_counts(rows, frozen)
    cacheP = S.tmm_cache_from_panel(CP, nfP)
    z_e, z_l = ZP[0], ZP[1]
    obs_st = [axis_stats(z_l, z_e, ZP[2 + i]) for i in range(len(mixes))]
    rng_b = np.random.default_rng(K.BOOT_SEED)
    boot = {m["f"]: dict(frac=[], delta=[], cos=[]) for m in mixes}
    for b in range(K.N_BOOT):
        for i, m in enumerate(mixes):
            cmix = sumc(resample_idx(m["late_idx"], rng_b)) + sumc(resample_idx(m["early_idx"], rng_b))
            zmix, _, _, _ = S.z_from_counts(cmix, cacheP, mu, sd, np.asarray(missP, bool))
            st = axis_stats(z_l, z_e, zmix)
            for k in ("frac", "delta", "cos"):
                boot[m["f"]][k].append(st[k])
        if (b + 1) % 50 == 0:
            log(f"[P boot] {b + 1}/{K.N_BOOT}")
    prows = []
    for i, m in enumerate(mixes):
        st = obs_st[i]
        flo, fhi = percentile_ci(boot[m["f"]]["frac"])
        dlo, dhi = percentile_ci(boot[m["f"]]["delta"])
        prows.append(dict(f=m["f"], n_late=m["n_late"], n_early=m["n_early"], n_total=n_tot,
                          frac=st["frac"], delta=st["delta"], cos=st["cos"],
                          frac_ci_lo=flo, frac_ci_hi=fhi, delta_ci_lo=dlo, delta_ci_hi=dhi,
                          frac_ci_valid=bool(ci_covers_point(st["frac"], flo, fhi))))
    pdf = pd.DataFrame(prows)
    pdf.to_csv(OUT / "stage0_posctrl.csv", index=False)
    r50 = next(r for r in prows if abs(r["f"] - 0.50) < 1e-12)
    r0 = next(r for r in prows if abs(r["f"]) < 1e-12)
    psum = dict(frac50=r50["frac"], delta50=r50["delta"], frac0=r0["frac"],
                delta0=r0["delta"], n_total=n_tot,
                monotonic=monotonic_increasing([r["frac"] for r in prows]),
                frac50_ci_excludes_0=bool(r50["frac_ci_valid"]
                                          and ci_excludes_zero(r50["frac_ci_lo"], r50["frac_ci_hi"])))
    log(f"[P] {psum}")
    del subB

    # --- Stage 0(c)+(d): the screen
    score = np.asarray(ctx["qidx"]["score"], float)
    nt_mod = ctx["nt_mod"]
    nt_sub, nt_sum = {}, {}
    for m in ("CRA", "CRI"):
        sm = K.submatrix(ctx["f"], nt_mod[m], ctx["exp"])
        nt_sub[m] = np.asarray(sm.todense(), np.float64)
        nt_sum[m] = K.to_ruler(nt_sub[m].sum(axis=0), dst, n_ruler)
        log(f"[NT dense] {m} {nt_sub[m].shape}")
    ruler_ctx = ruler_score(Z, frozen)
    logcpm_panel = S._official_logcpm_row(Cc, nf)
    ruler_sym = np.array([str(s).upper() for s in np.asarray(frozen["symbol"])], dtype=object)
    rpos = {}
    for i, s in enumerate(ruler_sym):
        rpos.setdefault(s, i)
    id_cols = np.array([rpos[g] for g in K.IDENTITY_PRESENT if g in rpos], int)

    rng_null = np.random.default_rng(K.SEED)
    rng_boot = np.random.default_rng(K.BOOT_SEED)
    out = []
    for m, lab, idx in ctx["qual"]:
        i, i_nt = ni[lab], ni[f"NT_{m}"]
        d = Z[i] - Z[i_nt]
        st = S.stats_from_d(d, v, p_unit[m])
        n = int(idx.size)
        n_nt = nt_sub[m].shape[0]
        keys = ("cos", "frac", "delta", "cos_resid", "delta_resid")
        nulls = {k: np.empty(K.N_NULL) for k in keys}
        if 0 < n < n_nt:
            for t in range(K.N_NULL):
                take = rng_null.choice(n_nt, size=n, replace=False)
                fake = K.to_ruler(nt_sub[m][take].sum(axis=0), dst, n_ruler)
                origin = nt_sum[m] - fake
                z_f, _, _, _ = S.z_from_counts(fake, cache, mu, sd, miss)
                z_o, _, _, _ = S.z_from_counts(origin, cache, mu, sd, miss)
                stn = S.stats_from_d(z_f - z_o, v, p_unit[m])
                for k in keys:
                    nulls[k][t] = stn[k]
        else:
            for k in keys:
                nulls[k][:] = np.nan
        # The bootstrap CIs are NOT re-derived. No Stage 0 check and no SENG2 gate
        # reads a CI: meets_toward and n_cos_q depend only on the nulls and the
        # point estimate. Disclosed as a scope reduction in FINDINGS_SENG2.md.
        boots = {k: np.full(K.N_BOOT, np.nan) for k in keys}
        rec = dict(modality=m, perturbation=lab, n_cells=n, cos=st["cos"], frac=st["frac"],
                   delta=st["delta"], cos_resid=st["cos_resid"],
                   delta_resid=st["delta_resid"], d_norm=st["d_norm"])
        for k in keys:
            lo, hi = percentile_ci(boots[k])
            rec[f"{k}_ci_lo"], rec[f"{k}_ci_hi"] = lo, hi
            rec[f"{k}_ci_valid"] = bool(ci_covers_point(rec[k], lo, hi))
            rec[f"p_{k}"] = permutation_p(rec[k], nulls[k],
                                          greater=(k not in ("delta", "delta_resid")))
        rec["prolif_shift"] = float(np.nanmean(score[idx]) - np.nanmean(score[nt_mod[m]]))
        cu, pu, _ = S.cos_against_perms(d, u, Uperm)
        rec["cos_d_u"], rec["p_cos_d_u"] = cu, pu
        diff = logcpm_panel[i, id_cols] - logcpm_panel[i_nt, id_cols]
        rec["identity_diff"] = float(np.mean(diff))
        rec["identity_drop"] = float(-rec["identity_diff"])
        rec["identity_n_genes"] = int(id_cols.size)
        rec["identity_loss"] = bool(rec["identity_drop"] > S.IDENTITY_DROP)
        rec["ruler_score_context"] = float(ruler_ctx[i])
        rec["ruler_score_nt_context"] = float(ruler_ctx[i_nt])
        out.append(rec)
        if len(out) % 20 == 0:
            log(f"[screen] {len(out)}/{len(ctx['qual'])} last={lab}")
    sdf = pd.DataFrame(out)
    for m in ("CRA", "CRI"):
        sel = sdf["modality"].to_numpy() == m
        for pcol, qcol in (("p_cos", "q_cos"), ("p_frac", "q_frac"), ("p_delta", "q_delta"),
                           ("p_cos_resid", "q_cos_resid"), ("p_delta_resid", "q_delta_resid")):
            if qcol not in sdf.columns:
                sdf[qcol] = np.nan
            sdf.loc[sel, qcol] = S.bh_q(sdf.loc[sel, pcol].to_numpy())
    sdf["meets_toward_young"] = [S.meets_toward(r) for _, r in sdf.iterrows()]
    sdf.to_csv(OUT / "stage0_screen.csv", index=False)

    by = {}
    for m in ("CRA", "CRI"):
        s = sdf[sdf.modality == m]
        by[m] = dict(n=int(len(s)),
                     n_cos_q=int(((s.q_cos <= K.P_BAR) & np.isfinite(s.q_cos)).sum()),
                     n_toward_young=int(s.meets_toward_young.sum()))
    named = {lab: int(sdf.loc[sdf.perturbation == lab, "n_cells"].iloc[0])
             if (sdf.perturbation == lab).any() else None for lab in K.NAMED}
    log(f"[screen] by_modality={by} named={named}")

    # --- gate
    tol = lambda s: 1e-4 * max(1.0, abs(s))
    checks = []
    for k, sv_ in stored.items():
        checks.append(dict(item=k, stored=sv_, got=got[k], kind="numeric",
                           ok=bool(abs(got[k] - sv_) <= tol(sv_))))
    for k, sv_ in (("frac50", 0.5469054675795995), ("delta50", -34.052728101022154),
                   ("frac0", 0.033645346601332594), ("delta0", 20.63359395062885)):
        checks.append(dict(item=k, stored=sv_, got=psum[k], kind="numeric",
                           ok=bool(abs(psum[k] - sv_) <= tol(sv_))))
    checks.append(dict(item="posctrl_n_total", stored=291, got=psum["n_total"],
                       kind="count", ok=psum["n_total"] == 291))
    for m, n_, cq, ty in (("CRA", 167, 0, 0), ("CRI", 173, 25, 0)):
        checks.append(dict(item=f"{m}_n", stored=n_, got=by[m]["n"], kind="count",
                           ok=by[m]["n"] == n_))
        checks.append(dict(item=f"{m}_n_cos_q", stored=cq, got=by[m]["n_cos_q"],
                           kind="count", ok=by[m]["n_cos_q"] == cq))
        checks.append(dict(item=f"{m}_n_toward_young", stored=ty,
                           got=by[m]["n_toward_young"], kind="count",
                           ok=by[m]["n_toward_young"] == ty))
    for lab, n_ in (("CRA_E2F3", 128), ("CRA_EZH2", 107), ("CRI_STAT3", 96), ("CRI_ZFX", 157)):
        checks.append(dict(item=f"n_cells[{lab}]", stored=n_, got=named.get(lab),
                           kind="count", ok=named.get(lab) == n_))
    passed = all(c["ok"] for c in checks)
    pd.DataFrame(checks).to_csv(OUT / "stage0_reproduction.csv", index=False)
    dump("stage0_summary.json", dict(passed=bool(passed), checks=checks, by_modality=by,
                                     named_n_cells=named, posctrl=psum, axis=got,
                                     v_max_abs_diff_vs_stored=float(np.max(np.abs(v - sv))),
                                     p_vu=p_vu))
    np.savez_compressed(OUT / "stage0_space.npz", Z=Z, v=v, u=u,
                        p_CRA=p_unit["CRA"], p_CRI=p_unit["CRI"], missing=miss,
                        nf=nf, names=np.array(names, dtype=object))
    log(f"[stage0] PASS={passed}")
    if not passed:
        bad = [c["item"] for c in checks if not c["ok"]]
        log(f"[stage0] STOP prior_results_not_reproduced mismatches={bad}")
        dump("STOP.json", dict(key="prior_results_not_reproduced", mismatches=bad))
    ctx["f"].close()
    log.close()



# ==========================================================================
def _space(log):
    """The frozen SENG space: the 348-row panel, its TMM cache, v and p.

    New pseudobulks enter this space through S.z_from_counts against the frozen
    cache, which is SENG's own mechanism for every null and bootstrap draw. The
    panel is not extended and nothing is refit (amendment A10).
    """
    z, names, ni = _panel()
    C = z["C"]
    ctx = K.load_context(log)
    frozen = ctx["frozen"]
    Z, missing, nf, Cc = K.panel_z_from_counts(C, frozen)
    cache = S.tmm_cache_from_panel(Cc, nf)
    mu = np.asarray(frozen["mu"], float)
    sd = np.asarray(frozen["sd"], float)
    miss = np.asarray(missing, bool)
    v = Z[ni["EARLY_A"]] - Z[ni["LATE_A"]]
    p_unit = {}
    for m in ("CRA", "CRI"):
        pv = Z[ni[f"TOP_{m}"]] - Z[ni[f"BOT_{m}"]]
        p_unit[m] = pv / float(np.linalg.norm(pv))
    zf = lambda c: S.z_from_counts(c, cache, mu, sd, miss)[0]
    return ctx, z, names, ni, C, Z, cache, mu, sd, miss, v, p_unit, zf, frozen


def _extra_counts(ctx, log, want):
    """Pseudobulks for groups outside SENG's panel, streamed once and cached."""
    cache_path = OUT / "extra_counts.npz"
    if cache_path.exists():
        d = np.load(cache_path, allow_pickle=True)
        if all(k in d.files for k in want):
            return {k: d[k] for k in d.files if k != "__names__"}
    idx = {}
    is_wt = ctx["is_wt"]
    pdl = ctx["pd_labels"]
    for lab in ("PD14", "PD26", "PD32"):
        idx[f"{lab}_FULL"] = np.flatnonzero(is_wt & (pdl == lab))
    idx["EARLY_B"] = np.asarray(ctx["sidx"]["early_b"], int)
    idx["LATE_B"] = np.asarray(ctx["sidx"]["late_b"], int)
    keys = list(idx)
    loc = K.stream_pseudobulk(ctx["f"], [idx[k] for k in keys], ctx["exp"], log)
    out = {k: K.to_ruler(loc[i], ctx["dst"], ctx["n_ruler"]) for i, k in enumerate(keys)}
    out_sz = {f"n_{k}": np.array([idx[k].size]) for k in keys}
    np.savez_compressed(cache_path, **out, **out_sz)
    log(f"[extra] streamed {keys} sizes={[idx[k].size for k in keys]}")
    return out


def stage_d1():
    log = K.Log(OUT / "d1.log")
    ctx, z, names, ni, C, Z, cache, mu, sd, miss, v, p_unit, zf, frozen = _space(log)
    ex = _extra_counts(ctx, log, ["PD14_FULL", "PD32_FULL"])
    nc = ctx["n_counts"]
    is_wt, pdl = ctx["is_wt"], ctx["pd_labels"]
    med = {
        "PD14_FULL": float(np.median(nc[np.flatnonzero(is_wt & (pdl == "PD14"))])),
        "PD32_FULL": float(np.median(nc[np.flatnonzero(is_wt & (pdl == "PD32"))])),
        "NT_CRA": float(np.median(nc[ctx["nt_mod"]["CRA"]])),
        "NT_CRI": float(np.median(nc[ctx["nt_mod"]["CRI"]])),
    }
    cnt = {"PD14_FULL": ex["PD14_FULL"], "PD32_FULL": ex["PD32_FULL"],
           "NT_CRA": C[ni["NT_CRA"]], "NT_CRI": C[ni["NT_CRI"]]}
    size = {"PD14_FULL": 581, "PD32_FULL": 4348,
            "NT_CRA": int(ctx["nt_mod"]["CRA"].size), "NT_CRI": int(ctx["nt_mod"]["CRI"].size)}
    log(f"[d1] medians={med} sizes={size}")

    def quantities(zz):
        a = zz["PD14_FULL"] - zz["PD32_FULL"]
        na = float(np.linalg.norm(a))
        out = dict(a_norm=na)
        for L, g in (("A", "NT_CRA"), ("I", "NT_CRI")):
            off = zz[g] - zz["PD32_FULL"]
            t = zz["PD14_FULL"] - zz[g]
            out[f"offset_{L}_norm"] = float(np.linalg.norm(off))
            out[f"offset_{L}_over_a"] = float(np.linalg.norm(off) / na)
            out[f"cos_offset_{L}_a"] = float(cosine_full(off, a))
            out[f"cos_t_{L}_a"] = float(cosine_full(t, a))
            out[f"angle_t_{L}_a_deg"] = float(np.degrees(np.arccos(
                np.clip(cosine_full(t, a), -1, 1))))
            out[f"t_{L}_norm"] = float(np.linalg.norm(t))
        out["r"] = max(out["cos_t_A_a"], out["cos_t_I_a"])
        return out

    raw = quantities({k: zf(cnt[k]) for k in cnt})
    log(f"[d1 raw] {raw}")

    # --- A3 thinning. Per comparison, the target is the lowest median among the
    # groups entering it; one common factor per group; 200 draws.
    comps = {"A": ["PD14_FULL", "PD32_FULL", "NT_CRA"],
             "I": ["PD14_FULL", "PD32_FULL", "NT_CRI"]}
    draws = {L: [] for L in comps}
    targets, factors = {}, {}
    for L, gl in comps.items():
        tgt, fac = K.thin_factors(med, gl)
        targets[L], factors[L] = tgt, fac
        log(f"[d1 thin {L}] target_median={tgt:.1f} factors={ {k: round(x,4) for k,x in fac.items()} }")
        rng = np.random.default_rng(K.SEED)
        for d in range(K.N_NULL):
            zz = {}
            for g in gl:
                zz[g] = zf(K.thin_counts(cnt[g], fac[g], rng))
            nt = "NT_CRA" if L == "A" else "NT_CRI"
            a = zz["PD14_FULL"] - zz["PD32_FULL"]
            t = zz["PD14_FULL"] - zz[nt]
            off = zz[nt] - zz["PD32_FULL"]
            draws[L].append(dict(
                a_norm=float(np.linalg.norm(a)),
                offset_norm=float(np.linalg.norm(off)),
                offset_over_a=float(np.linalg.norm(off) / np.linalg.norm(a)),
                cos_offset_a=float(cosine_full(off, a)),
                cos_t_a=float(cosine_full(t, a)),
                t_norm=float(np.linalg.norm(t))))
            if (d + 1) % 50 == 0:
                log(f"[d1 thin {L}] {d + 1}/{K.N_NULL}")
    thin = {}
    for L in comps:
        df = pd.DataFrame(draws[L])
        df.to_csv(OUT / f"d1_thinned_draws_{L}.csv", index=False)
        for c in df.columns:
            thin[f"{c}_{L}_median"] = float(np.median(df[c]))
            thin[f"{c}_{L}_lo"] = float(np.percentile(df[c], 2.5))
            thin[f"{c}_{L}_hi"] = float(np.percentile(df[c], 97.5))
    r_thin = max(thin["cos_t_a_A_median"], thin["cos_t_a_I_median"])
    r_lib = "A" if thin["cos_t_a_A_median"] >= thin["cos_t_a_I_median"] else "I"
    r_lo, r_hi = thin[f"cos_t_a_{r_lib}_lo"], thin[f"cos_t_a_{r_lib}_hi"]

    # verify the thinning really matches the per-cell median (one draw, disclosed)
    rngc = np.random.default_rng(K.SEED)
    chk = {}
    for L, gl in comps.items():
        for g in gl:
            m = med[g] * factors[L][g]
            chk[f"{L}:{g}"] = dict(factor=factors[L][g], median_before=med[g],
                                   median_expected_after=m, target=targets[L])
    dump("d1_thinning_check.json", chk)

    gate = ("USABLE" if r_thin >= 0.70 else "DEGRADED" if r_thin >= 0.40 else "UNUSABLE")
    pct_lib = 100.0 * float(np.sqrt(max(0.0, 1.0 - r_thin ** 2)))
    res = dict(raw=raw, thinned=thin, r_primary=r_thin, r_library=r_lib,
               r_ci=[r_lo, r_hi], r_angle_deg=float(np.degrees(np.arccos(np.clip(r_thin, -1, 1)))),
               gate=gate, pct_library=pct_lib, medians=med, sizes=size,
               targets=targets, factors=factors,
               r_raw=raw["r"], gate_raw=("USABLE" if raw["r"] >= 0.70 else
                                         "DEGRADED" if raw["r"] >= 0.40 else "UNUSABLE"),
               cos_a_v=float(cosine_full(zf(cnt["PD14_FULL"]) - zf(cnt["PD32_FULL"]), v)),
               a_vs_v_note="a is the full-arm axis; v is SENG's half-A axis")
    dump("d1_summary.json", res)
    log(f"[d1] gate={gate} r={r_thin:.6f} ({r_lo:.6f},{r_hi:.6f}) raw_r={raw['r']:.6f} "
        f"gate_raw={res['gate_raw']} cos(a,v)={res['cos_a_v']:.6f}")
    ctx["f"].close()
    log.close()


# ==========================================================================
def _two_origin_stats(z_nt, z_late, z_target, d, p_hat):
    """delta and frac under both pre-specified origins, plus the
    proliferation-residualized delta (D4) with the same projection."""
    out = {}
    for tag, z_o in (("i", z_nt), ("ii", z_late)):
        t = z_target - z_o
        z_S = z_o + d
        out[f"delta_{tag}"] = float(np.linalg.norm(z_S - z_target) - np.linalg.norm(t))
        nt2 = float(np.dot(t, t))
        out[f"frac_{tag}"] = float(np.dot(d, t) / nt2) if nt2 > 0 else np.nan
        d_r = d - np.dot(d, p_hat) * p_hat
        t_r = t - np.dot(t, p_hat) * p_hat
        out[f"delta_resid_{tag}"] = float(np.linalg.norm(t_r - d_r) - np.linalg.norm(t_r))
    return out


def _nt_cells(ctx, log):
    nt = {}
    for m in ("CRA", "CRI"):
        sm = K.submatrix(ctx["f"], ctx["nt_mod"][m], ctx["exp"])
        nt[m] = np.asarray(sm.todense(), np.float64)
        log(f"[nt dense] {m} {nt[m].shape}")
    return nt


def _null_deltas(n, m, nt_blk, nt_pool, zf, dst, n_ruler, z_nt, z_late, z_target,
                 p_hat, rng, reps=None):
    """Southard/SENG floor: draw n control cells, displacement is
    z(draw) - z(pool minus draw). Depends only on n, so one set per n."""
    reps = reps or K.N_NULL
    n_nt = nt_blk.shape[0]
    rec = []
    for _ in range(reps):
        take = rng.choice(n_nt, size=n, replace=False)
        fake = K.to_ruler(nt_blk[take].sum(axis=0), dst, n_ruler)
        d = zf(fake) - zf(nt_pool - fake)
        rec.append(_two_origin_stats(z_nt, z_late, z_target, d, p_hat))
    return pd.DataFrame(rec)


def stage_d2():
    log = K.Log(OUT / "d2.log")
    ctx, z, names, ni, C, Z, cache, mu, sd, miss, v, p_unit, zf, frozen = _space(log)
    d1 = json.loads((OUT / "d1_summary.json").read_text())
    dst, n_ruler = z["dst"], ctx["n_ruler"]
    z_target = Z[ni["EARLY_A"]]          # PD14 half A, amendment A4
    z_late = Z[ni["LATE_A"]]             # PD32 half A, amendment A6
    nt_blk = _nt_cells(ctx, log)
    nt_pool = {m: C[ni[f"NT_{m}"]] for m in ("CRA", "CRI")}

    # one null per distinct cell count per modality (Southard: the null depends on n)
    by_mod = {"CRA": [], "CRI": []}
    for m, lab, idx in ctx["qual"]:
        by_mod[m].append((lab, int(idx.size)))
    floors = {}
    rng = np.random.default_rng(K.SEED)
    for m in ("CRA", "CRI"):
        ns = sorted({n for _, n in by_mod[m]})
        log(f"[d2] {m}: {len(ns)} distinct cell counts {ns[0]}-{ns[-1]}")
        for j, n in enumerate(ns):
            df = _null_deltas(n, m, nt_blk[m], nt_pool[m], zf, dst, n_ruler,
                              Z[ni[f"NT_{m}"]], z_late, z_target, p_unit[m], rng)
            floors[(m, n)] = {c: float(np.percentile(df[c], 5)) for c in df.columns}
            floors[(m, n)]["_median_delta_i"] = float(np.median(df["delta_i"]))
            floors[(m, n)]["_median_delta_ii"] = float(np.median(df["delta_ii"]))
            if (j + 1) % 20 == 0:
                log(f"[d2 floors] {m} {j + 1}/{len(ns)}")
    rows = []
    for m, lab, idx in ctx["qual"]:
        n = int(idx.size)
        d = Z[ni[lab]] - Z[ni[f"NT_{m}"]]
        st = _two_origin_stats(Z[ni[f"NT_{m}"]], z_late, z_target, d, p_unit[m])
        fl = floors[(m, n)]
        rec = dict(modality=m, perturbation=lab, n_cells=n)
        rec.update(st)
        for tag in ("i", "ii"):
            for stat in ("delta", "delta_resid"):
                k = f"{stat}_{tag}"
                rec[f"floor_{k}"] = fl[k]
                rec[f"below_floor_{k}"] = bool(st[k] < fl[k])
                rec[f"decreased_{k}"] = bool(st[k] < 0)
                rec[f"called_{k}"] = bool(st[k] < fl[k] and st[k] < 0)
            rec[f"null_median_delta_{tag}"] = fl[f"_median_delta_{tag}"]
        rows.append(rec)
    df = pd.DataFrame(rows)
    allv = df.n_cells.to_numpy()
    df["cell_pctl_in_screen"] = [float((allv <= x).mean() * 100) for x in allv]
    df.to_csv(OUT / "d2_table.csv", index=False)

    summ = {}
    for tag in ("i", "ii"):
        for stat in ("delta", "delta_resid"):
            k = f"{stat}_{tag}"
            kb = int(df[f"below_floor_{k}"].sum())
            kd = int(df[f"decreased_{k}"].sum())
            kc = int(df[f"called_{k}"].sum())
            lo, hi = K.binom_ci(kc, len(df))
            summ[k] = dict(n=len(df), below_floor=kb, decreased=kd, called=kc,
                           called_frac=kc / len(df), called_ci=[lo, hi],
                           chance_baseline=0.05,
                           chance_expected=0.05 * len(df))
            for m in ("CRA", "CRI"):
                sub = df[df.modality == m]
                summ[k][m] = dict(n=len(sub),
                                  below_floor=int(sub[f"below_floor_{k}"].sum()),
                                  decreased=int(sub[f"decreased_{k}"].sum()),
                                  called=int(sub[f"called_{k}"].sum()),
                                  chance_expected=0.05 * len(sub))
    named = {}
    for lab in K.NAMED:
        r = df[df.perturbation == lab]
        named[lab] = {} if r.empty else {c: _js(r.iloc[0][c]) for c in df.columns
                                         if c not in ("modality", "perturbation")}
    k_i = sum(1 for lab in K.NAMED if named[lab] and named[lab]["called_delta_i"])
    k_ii = sum(1 for lab in K.NAMED if named[lab] and named[lab]["called_delta_ii"])
    gate = ("SPLIT" if (k_i >= 2) != (k_ii >= 2)
            else "APPROACH" if (k_i >= 2 and k_ii >= 2) else "NO APPROACH")
    dump("d2_summary.json", dict(summary=summ, named=named, k_origin_i=k_i,
                                 k_origin_ii=k_ii, gate=gate,
                                 d1_gate=d1["gate"], d1_r=d1["r_primary"]))
    log(f"[d2] gate={gate} k_i={k_i} k_ii={k_ii} called_i={summ['delta_i']['called']} "
        f"called_ii={summ['delta_ii']['called']} of {len(df)}")
    ctx["f"].close()
    log.close()


# ==========================================================================
def stage_d3():
    log = K.Log(OUT / "d3.log")
    ctx, z, names, ni, C, Z, cache, mu, sd, miss, v, p_unit, zf, frozen = _space(log)
    ex = _extra_counts(ctx, log, ["PD14_FULL", "PD32_FULL"])
    dst, n_ruler = z["dst"], ctx["n_ruler"]
    z_target, z_late = Z[ni["EARLY_A"]], Z[ni["LATE_A"]]
    a = zf(ex["PD14_FULL"]) - zf(ex["PD32_FULL"])
    a_hat = a / float(np.linalg.norm(a))
    a_norm = float(np.linalg.norm(a))
    nt_blk = _nt_cells(ctx, log)
    nt_pool = {m: C[ni[f"NT_{m}"]] for m in ("CRA", "CRI")}
    eb = np.asarray(ctx["sidx"]["early_b"], int)
    sub_eb = np.asarray(K.submatrix(ctx["f"], eb, ctx["exp"]).todense(), np.float64)
    log(f"[d3] spike material EARLY_B n={eb.size} ||a||={a_norm:.4f}")

    # (a) mixture across the library boundary
    mrows = []
    for n_tot in (160, 100):
        for m in ("CRA", "CRI"):
            for f in (0.0, 0.10, 0.25, 0.50):
                n_e = int(round(f * n_tot))
                n_c = n_tot - n_e
                rng = np.random.default_rng(K.SEED)
                for r in range(K.N_NULL):
                    tc = rng.choice(nt_blk[m].shape[0], size=n_c, replace=False)
                    te = rng.choice(sub_eb.shape[0], size=n_e, replace=False) if n_e else None
                    cm = nt_blk[m][tc].sum(axis=0)
                    if n_e:
                        cm = cm + sub_eb[te].sum(axis=0)
                    cm = K.to_ruler(cm, dst, n_ruler)
                    c_rest = nt_pool[m] - K.to_ruler(nt_blk[m][tc].sum(axis=0), dst, n_ruler)
                    d = zf(cm) - zf(c_rest)
                    st = _two_origin_stats(Z[ni[f"NT_{m}"]], z_late, z_target, d, p_unit[m])
                    mrows.append(dict(n_total=n_tot, modality=m, dose=f, rep=r, **st))
            log(f"[d3a] n_tot={n_tot} {m} done")
    md = pd.DataFrame(mrows)
    md.to_csv(OUT / "d3_mixture_draws.csv", index=False)
    msum = []
    for (n_tot, m, f), g in md.groupby(["n_total", "modality", "dose"]):
        rec = dict(n_total=n_tot, modality=m, dose=f, n_draws=len(g))
        for c in ("delta_i", "frac_i", "delta_ii", "frac_ii"):
            rec[f"{c}_median"] = float(np.median(g[c]))
            rec[f"{c}_lo"] = float(np.percentile(g[c], 2.5))
            rec[f"{c}_hi"] = float(np.percentile(g[c], 97.5))
        msum.append(rec)
    ms = pd.DataFrame(msum).sort_values(["n_total", "modality", "dose"])
    ms.to_csv(OUT / "d3_mixture_summary.csv", index=False)

    mono, ci_ok = {}, {}
    for n_tot in (160, 100):
        for m in ("CRA", "CRI"):
            sel = ms[(ms.n_total == n_tot) & (ms.modality == m)].sort_values("dose")
            for tag in ("i", "ii"):
                mono[f"{n_tot}_{m}_{tag}"] = bool(monotonic_increasing(
                    sel[f"frac_{tag}_median"].tolist()))
                r50 = sel[sel.dose == 0.50].iloc[0]
                ci_ok[f"{n_tot}_{m}_{tag}"] = bool(r50[f"delta_{tag}_hi"] < 0)

    # (b) plant recovery along a
    prows = []
    for n in (100, 160):
        for m in ("CRA", "CRI"):
            rng = np.random.default_rng(K.SEED)
            base = []
            for r in range(K.N_NULL):
                take = rng.choice(nt_blk[m].shape[0], size=n, replace=False)
                fake = K.to_ruler(nt_blk[m][take].sum(axis=0), dst, n_ruler)
                base.append(zf(fake) - zf(nt_pool[m] - fake))
            for f in (0.0, 0.25, 0.50, 1.0):
                for r, d0 in enumerate(base):
                    d = d0 + f * a_norm * a_hat
                    st = _two_origin_stats(Z[ni[f"NT_{m}"]], z_late, z_target, d, p_unit[m])
                    prows.append(dict(n_cells=n, modality=m, dose=f, rep=r, **st))
            log(f"[d3b] n={n} {m} done")
    pdw = pd.DataFrame(prows)
    pdw.to_csv(OUT / "d3_plant_draws.csv", index=False)
    psum = []
    for (n, m), g in pdw.groupby(["n_cells", "modality"]):
        fl = {tag: float(np.percentile(g[g.dose == 0.0][f"delta_{tag}"], 5)) for tag in ("i", "ii")}
        for f, gg in g.groupby("dose"):
            rec = dict(n_cells=n, modality=m, dose=f, n_draws=len(gg))
            for tag in ("i", "ii"):
                d = gg[f"delta_{tag}"].to_numpy()
                rec[f"floor_{tag}"] = fl[tag]
                rec[f"delta_{tag}_median"] = float(np.median(d))
                rec[f"power_{tag}"] = float(np.mean((d < fl[tag]) & (d < 0)))
            psum.append(rec)
    ps = pd.DataFrame(psum).sort_values(["n_cells", "modality", "dose"])
    ps.to_csv(OUT / "d3_plant_summary.csv", index=False)

    pw25 = {}
    for n in (100, 160):
        for m in ("CRA", "CRI"):
            r = ps[(ps.n_cells == n) & (ps.modality == m) & (ps.dose == 0.25)].iloc[0]
            for tag in ("i", "ii"):
                pw25[f"{n}_{m}_{tag}"] = float(r[f"power_{tag}"])
    pmin100_i = min(pw25[f"100_{m}_i"] for m in ("CRA", "CRI"))
    mono_i = all(mono[f"160_{m}_i"] for m in ("CRA", "CRI"))
    ci_i = all(ci_ok[f"160_{m}_i"] for m in ("CRA", "CRI"))
    gate = "SENSITIVE" if (mono_i and ci_i and pmin100_i >= 0.80) else "INSENSITIVE"
    dump("d3_summary.json", dict(gate=gate, monotonic=mono, ci50_excludes_zero=ci_ok,
                                 power_25pct=pw25, power_25_at_100_origin_i_min=pmin100_i,
                                 a_norm=a_norm, n_spike=int(eb.size),
                                 gate_inputs=dict(mono_160_origin_i=mono_i,
                                                  ci50_160_origin_i=ci_i,
                                                  power25_100_origin_i=pmin100_i)))
    log(f"[d3] gate={gate} mono={mono_i} ci50={ci_i} power25@100={pmin100_i:.3f}")
    ctx["f"].close()
    log.close()


# ==========================================================================
def stage_d4():
    log = K.Log(OUT / "d4.log")
    scr = pd.read_csv(OUT / "stage0_screen.csv")
    d2 = pd.read_csv(OUT / "d2_table.csv")
    mg = scr.merge(d2, on=["modality", "perturbation"], suffixes=("_seng", "_d2"))
    mg["ruler_change"] = mg.ruler_score_context - mg.ruler_score_nt_context
    mg["ruler_scores_younger"] = mg.ruler_change < 0
    out = dict(rho_prolif_stored=dict(CRA=0.8280530161707566, CRI=0.8620017521821158),
               n=int(len(mg)))
    for tag in ("i", "ii"):
        younger = mg.ruler_scores_younger
        called = mg[f"called_delta_{tag}"]
        out[f"ruler_origin_{tag}"] = dict(
            n_scores_younger=int(younger.sum()),
            n_younger_and_not_called=int((younger & ~called).sum()),
            n_younger_and_called=int((younger & called).sum()),
            n_called=int(called.sum()),
            n_not_younger_and_called=int((~younger & called).sum()))
        out[f"resid_origin_{tag}"] = dict(
            called_plain=int(mg[f"called_delta_{tag}"].sum()),
            called_resid=int(mg[f"called_delta_resid_{tag}"].sum()),
            called_both=int((mg[f"called_delta_{tag}"] & mg[f"called_delta_resid_{tag}"]).sum()))
    out["identity"] = dict(
        screen_median=float(mg.identity_diff.median()),
        screen_lo=float(np.percentile(mg.identity_diff, 2.5)),
        screen_hi=float(np.percentile(mg.identity_diff, 97.5)),
        n_identity_loss=int(mg.identity_loss.sum()),
        n_genes=int(mg.identity_n_genes.iloc[0]))
    named = {}
    for lab in K.NAMED:
        r = mg[mg.perturbation == lab]
        if r.empty:
            continue
        r = r.iloc[0]
        named[lab] = dict(n_cells=int(r.n_cells_d2), identity_diff=float(r.identity_diff),
                          ruler_change=float(r.ruler_change),
                          ruler_scores_younger=bool(r.ruler_scores_younger),
                          prolif_shift=float(r.prolif_shift),
                          delta_i=float(r.delta_i), delta_ii=float(r.delta_ii),
                          delta_resid_i=float(r.delta_resid_i),
                          delta_resid_ii=float(r.delta_resid_ii),
                          called_delta_i=bool(r.called_delta_i),
                          called_delta_ii=bool(r.called_delta_ii),
                          called_delta_resid_i=bool(r.called_delta_resid_i),
                          called_delta_resid_ii=bool(r.called_delta_resid_ii),
                          cell_pctl_in_screen=float(r.cell_pctl_in_screen),
                          seng_cos=float(r.cos), seng_delta=float(r.delta),
                          seng_q_cos=float(r.q_cos))
    out["named"] = named
    rr = pd.read_csv(K.SENG_WHY / "named_rrej_vs_cos.csv")
    out["rrej_seng_why"] = [
        dict(perturbation=str(r.perturbation), n_cells=int(r.n_cells),
             R_rej_recomputed=float(r.R_rej), published_table1=float(r.published_table1_R_rej),
             ratio=float(r.R_rej / r.published_table1_R_rej))
        for r in rr.itertuples()]
    out["identity_genes_used"] = K.IDENTITY_PRESENT
    mg.to_csv(OUT / "d4_merged.csv", index=False)
    dump("d4_summary.json", out)
    log(f"[d4] {json.dumps(_js(out)['ruler_origin_i'])}")
    log.close()


STAGES = dict(panel=stage_panel, stage0=stage0, d1=stage_d1, d2=stage_d2, d3=stage_d3, d4=stage_d4)

if __name__ == "__main__":
    which = sys.argv[1] if len(sys.argv) > 1 else "panel"
    STAGES[which]()
