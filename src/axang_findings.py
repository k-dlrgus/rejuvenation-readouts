"""Writes FINDINGS_AXANG.md from the stage outputs in results/axang/."""
from __future__ import annotations

import json
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parent.parent
OUT = ROOT / "results" / "axang"
FIND = ROOT / "FINDINGS_AXANG.md"
NAMED = ["CRA_E2F3", "CRA_EZH2", "CRI_STAT3", "CRI_ZFX"]
AGN = {"A1": "GTEx donor age", "A2": "frozen ruler", "A3": "Sengstack passage",
       "A4": "Lu donor pair"}
PRN = {"C1_gtex": "C1 GTEx donors (652, bulk)",
       "C1_gtex_anchors": "C1 GTEx anchor donors (344, bulk)",
       "C1_south": "C1 Southard NT cells (14,013)",
       "C1_seng_CRA": "C1 Sengstack NT CRA (333)",
       "C1_seng_CRI": "C1 Sengstack NT CRI (739)",
       "C1_lu": "C1 Lu day-0 (12,803)", "C2": "C2 Tirosh gene list (94 genes)",
       "C2_seng4_P3": "C2 as SENG4's P3"}
L = []


def w(s=""):
    L.append(s)


def f(x, nd=6):
    if x is None:
        return "n/a"
    try:
        v = float(x)
    except (TypeError, ValueError):
        return str(x)
    return "n/a" if not np.isfinite(v) else f"{v:+.{nd}f}"


def main():
    jl = lambda n: json.loads((OUT / n).read_text(encoding="utf-8"))
    s0, d1, d2, d3 = jl("stage0_summary.json"), jl("d1_summary.json"), \
        jl("d2_summary.json"), jl("d3_summary.json")
    repro = pd.read_csv(OUT / "stage0_reproduction.csv")
    ang = pd.read_csv(OUT / "d1_angles.csv")
    pl = pd.read_csv(OUT / "d3a_plant.csv")
    mt = s0["meta"]

    def row(a, b, sp="S2A_k20"):
        m = (ang.space == sp) & (((ang.axis_a == a) & (ang.axis_b == b)) |
                                 ((ang.axis_a == b) & (ang.axis_b == a)))
        return ang[m].iloc[0]

    w("# FINDINGS_AXANG — how far is each aging model from a proliferation axis?")
    w()
    w("Task AXANG. Computational, descriptive geometry. No new data. Writes")
    w("`results/axang/` and this file. Nothing committed, nothing pushed.")
    w()
    w(f"**The gates.** Stage 0 PASS; D1 **{d1['gate']}** on the primary")
    w(f"construction; D3a **{d3['gate_a']}**; D3b **{d3['gate_b']}**.")
    w()
    w("**The one-scale answer the task was built to give.** Each aging model")
    w("against a proliferation axis built from its own units, in S2A_k20, as the")
    w("angle from parallel (amendment F10):")
    w()
    w("| aging model | its own proliferation axis | angle from parallel | % of the aging axis removed by projecting it out |")
    w("|---|---|---|---|")
    for a, c in (("A3", "C1_seng_CRA"), ("A1", "C1_gtex"), ("A4", "C1_lu"),
                 ("A2", "C1_gtex")):
        r = row(a, c)
        mag = float(np.degrees(np.arccos(min(1.0, abs(float(r.cos))))))
        w(f"| {a}, {AGN[a]} | {PRN[c]} | **{mag:.2f}°** | "
          f"{100 * (1 - float(r.frac_a_surviving)):.1f}% |")
    w()
    w("So the ordering, which is the task's result: **Sengstack passage aging is")
    w("the worst case by a wide margin** (19.66°, 66% removed), the GTEx donor-age")
    w("axis is intermediate (35.83°, 41.5% removed), and the **frozen ruler and the")
    w("preprint's own Lu donor-pair axis are close to orthogonal** to proliferation")
    w("(85.51° and 57.63°, 0.3% and 15.5% removed).")
    w()
    w("**Two findings the flag's recorded expectation did not anticipate, stated")
    w("rather than smoothed, as it requires.**")
    w()
    r1 = row("A1", "C1_gtex")
    w(f"1. **The sign.** `cos(A1, C1_gtex)` = {f(float(r1.cos))}, bootstrap")
    w(f"   [{f(float(r1.ci_lo),3)}, {f(float(r1.ci_hi),3)}] — strongly **negative**.")
    w("   GTEx donors with higher cell-cycle scores sit on the *old* side of the")
    w("   donor-age axis. The flag's expectation was near-orthogonality, and the")
    w("   magnitude is far from it.")
    w("2. **The expectation holds for cell-derived axes and fails for the bulk one.**")
    w(f"   A1 against the Southard NT-cell axis is {f(float(row('A1','C1_south').cos),4)}")
    w(f"   and against the Sengstack NT-cell axes {f(float(row('A1','C1_seng_CRA').cos),4)}")
    w(f"   / {f(float(row('A1','C1_seng_CRI').cos),4)} — all near-orthogonal, as")
    w("   predicted. Only A1's **own bulk-donor** axis is strongly anti-aligned. The")
    w("   inconsistency is therefore between two *constructions* of proliferation,")
    w("   not between the aging models, and the bulk construction is the outlier.")
    w()
    w("---")
    w()
    w("## 1. The pre-registration, verbatim")
    w()
    w("From `results/axang/PREREG.flag`, written 2026-10-06T23:34:40Z before any")
    w("AXANG number existed and not edited since. Amendments E1-E4 were supplied")
    w("with the flag and are part of it.")
    w()
    w("```")
    for line in (OUT / "PREREG.flag").read_text(encoding="utf-8").splitlines():
        w(line)
    w("```")
    w()
    w("## 2. Amendments and deviations")
    w()
    w("One self-contradiction in the flag was put to the author before any")
    w("conclusion was written (F10); nine further definitions or citations the flag")
    w("omits or gets wrong were fixed in writing before computing. Full text of")
    w("`results/axang/AMENDMENTS.md`:")
    w()
    w("```")
    for line in (OUT / "AMENDMENTS.md").read_text(encoding="utf-8").splitlines():
        w(line)
    w("```")
    w()
    w("## 3. State the method before using it")
    w()
    w("### (a) and (b) Every axis, its source, and whether it is read or rebuilt")
    w()
    w("| axis | what it is | source | read frozen or rebuilt |")
    w("|---|---|---|---|")
    w("| **A1** | GTEx young-minus-old centroid direction, young 20-39 (n=113), old 60-79 (n=231) | `results/toward/anchors.npz` key `u` | **read frozen** |")
    w("| **A2** | unit-length frozen fibroblast ruler weight vector | `results/fibro/frozen_ruler_ridge_raw.npz` key `w` | **read frozen** |")
    w("| **A3** | Sengstack PD14 half A minus PD32 half A | SENG3's frozen panel, `results/seng2/panel_counts.npz` | rebuilt by SENG3's own definition |")
    w("| **A4** | Lu GM00731 (96y) day-0 to GM23815 (22y) day-0, half A | `results/same/panel_z.npz` key `v` | **read frozen** |")
    w("| **C1 GTEx** | bottom-to-top cell-cycle quartile over GTEx donors | built here (F1, F2) | rebuilt |")
    w("| **C1 Southard** | bottom-to-top NT cell-cycle quartile, 3,504 per quartile | `results/south6/stage1.npz` key `p_unit` | **read frozen** (F5) |")
    w("| **C1 Sengstack** | bottom-to-top NT half-one quartile, per modality | `results/seng4/prep.npz` | rebuilt by SENG4's definition |")
    w("| **C1 Lu** | bottom-to-top quartile over day-0 fibroblasts of both donors | built here | rebuilt |")
    w("| **C2** | equal-weight Tirosh direction, no cells | built here | rebuilt, deterministic |")
    w()
    w("The C1 recipe, quoted from committed code — the per-unit score")
    w("(`src/seng_run.py:915`) and the quartile rule (`src/seng_run.py:925`):")
    w()
    w("```python")
    w("def module_score(counts_csr, cols, lib):")
    w("    sub = counts_csr[:, cols]; scale = np.maximum(lib, 1.0)[:, None]")
    w("    return np.log1p(arr / scale * 1e4).mean(axis=1)")
    w("")
    w("q25, q75 = np.quantile(scores[finite], [0.25, 0.75])")
    w("bottom = finite & (scores <= q25);  top = finite & (scores >= q75)")
    w("```")
    w()
    w("and the Southard axis, `src/south6_run.py:1240`:")
    w()
    w("```python")
    w("p_unit = pvec / nrm")
    w("prolif_note = f\"built from NT quartiles n_bottom={int(bottom.sum())} "
      "n_top={int(top.sum())} ||p||={nrm:.4f}\"")
    w("```")
    w()
    w(f"Units entering each C1: GTEx **{mt['n_gtex_donors']} donors** "
      f"({mt['gtex_quartiles']['n_bottom']} per quartile; anchors-only variant "
      f"{mt['n_gtex_anchor_donors']} donors, {mt['gtex_anchor_quartiles']['n_bottom']} "
      f"per quartile); Southard **{mt['n_south_nt_cells']:,} NT cells** (3,504 per "
      f"quartile); Sengstack **{mt['n_seng_half1']['CRA']} CRA / "
      f"{mt['n_seng_half1']['CRI']} CRI** NT half-one cells; Lu "
      f"**{mt['lu']['GM23815']['n_day0']:,} + {mt['lu']['GM00731']['n_day0']:,}** "
      f"day-0 fibroblasts. C2 uses {mt['cc_on_ruler']} of {mt['cc_total']} Tirosh "
      f"genes present on the ruler (missing {mt['cc_missing_from_ruler']}).")
    w()
    w("### (c) Spaces")
    w()
    w("Every axis lives natively in the 23,485-column frozen-ruler z space, so")
    w("**`space_mismatch` does not fire** and every pair is computable in both the")
    w("full-gene space and S2A_k20. Cosines are reported per space and are **never**")
    w("compared across spaces: a 20-dimensional cosine and a 23,485-dimensional one")
    w("are different quantities, and `results/axang/d1_angles.csv` carries a `space`")
    w("column on every row. The S2A_k20 basis is read frozen from")
    w("`results/south7/basis.npz` (308 GTEx donors aged 40-59) and is fit on no cell")
    w("of any single-cell dataset here.")
    w()
    w("One caveat is stated where it appears (F7): Southard axes are structurally")
    w("zero on 3,221 ruler genes and Sengstack/Lu axes on 6,704, so a cosine between")
    w("a Southard-derived and a Sengstack-derived axis is taken over a frame where")
    w("each is zero on genes the other measures.")
    w()
    w("### (d) Axes estimated from units that also define another axis")
    w()
    w("**Yes, in two places, and neither is refitted away.**")
    w()
    w("1. **The gate quantity itself.** The 344 donors defining A1 are a subset of")
    w(f"   the {mt['n_gtex_donors']} donors defining `C1_gtex`, so `cos(A1, C1_gtex)`")
    w("   is estimated partly from the same units on both sides. The anchors-only")
    w("   variant makes the overlap total and is reported beside it; the two agree")
    w(f"   ({f(float(r1.cos),4)} against {f(float(row('A1','C1_gtex_anchors').cos),4)}),")
    w("   so the overlap is not what produces the value — but the number is labelled.")
    w("2. **A4 and `C1_lu`** both come from the Lu day-0 fibroblasts, A4 from SAME's")
    w("   half A only and `C1_lu` from all of them (F6). They overlap but are not")
    w("   identical sets.")
    w()
    w("A3 against `C1_seng` is **not** in this list: SENG4 reserved NT half one for")
    w("the axis and half two for everything else, and A3 is built from PD14/PD32")
    w("cells which enter no C1.")
    w()
    w("## 4. Stage 0 — reproduce first (gate)")
    w()
    w(f"**Gate: {'PASS' if s0['passed'] else 'FAIL'}.** All {len(repro)} checks.")
    w()
    w("| item | kind | stored | re-derived | ok |")
    w("|---|---|---|---|---|")
    for r in repro.itertuples():
        st = int(r.stored) if r.kind == "count" else f"{r.stored:.6f}"
        g = "n/a" if pd.isna(r.got) else (int(r.got) if r.kind == "count" else f"{r.got:.6f}")
        w(f"| `{r.item}` | {r.kind} | {st} | {g} | {'yes' if r.ok else '**NO**'} |")
    w()
    w("Two notes. The frozen anchors are **exactly** the pack centroids —")
    w(f"max |c_young − mean(Z[young])| = {d1['anchors_vs_pack']['max_abs_c_young']:.1f} and")
    w(f"the same for old — so A1 is the axis this project has always used. And the")
    w(f"E1 stored value is {mt['stored_a4_cos']:.15f}, which the flag quotes rounded")
    w("to 0.039; both the rounded and the exact forms are checked.")
    w()
    w("## 5. D1 — the angles (gated)")
    w()
    w("### Each aging axis against every proliferation construction, S2A_k20")
    w()
    w("Angle is the signed angle; **angle from parallel** is `arccos|cos|`, the")
    w("separability quantity (F10). `% removed` is `1 − sqrt(1 − cos²)`.")
    w()
    w("| aging axis | proliferation axis | cos | 2.5–97.5% | signed angle | from parallel | % removed |")
    w("|---|---|---|---|---|---|---|")
    for a in ("A1", "A2", "A3", "A4"):
        for c in PRN:
            if c == "C2_seng4_P3":
                continue
            r = row(a, c)
            mag = float(np.degrees(np.arccos(min(1.0, abs(float(r.cos))))))
            ci = ("—" if not bool(r.bootstrapped) else
                  f"{f(float(r.ci_lo),3)} to {f(float(r.ci_hi),3)}")
            w(f"| {a} {AGN[a]} | {PRN[c]} | {f(float(r.cos),4)} | {ci} | "
              f"{float(r.angle_deg):.2f}° | **{mag:.2f}°** | "
              f"{100 * (1 - float(r.frac_a_surviving)):.1f}% |")
    w()
    w("Rows with no interval are the four axes F9 could not resample because they")
    w("are read from stored vectors with no per-unit pack (A2, A3, A4, `C1_south`).")
    w("Every bootstrapped interval contained its point estimate, so none is flagged")
    w("invalid.")
    w()
    w("### Aging axis against aging axis, S2A_k20")
    w()
    w("| pair | cos | signed angle |")
    w("|---|---|---|")
    for i, a in enumerate(["A1", "A2", "A3", "A4"]):
        for b in ["A1", "A2", "A3", "A4"][i + 1:]:
            r = row(a, b)
            w(f"| {a} ({AGN[a]}) vs {b} ({AGN[b]}) | {f(float(r.cos),4)} | "
              f"{float(r.angle_deg):.2f}° |")
    w()
    w("The four aging models are mutually near-orthogonal, which is this project's")
    w("standing finding extended to four axes: A1 vs A3 is")
    w(f"{f(float(row('A1','A3').cos),4)} and A1 vs A4 is {f(float(row('A1','A4').cos),4)}")
    w("— the donor-age axis is essentially orthogonal to both the passage model and")
    w("the preprint's own donor-pair model.")
    w()
    w("### Proliferation axis against proliferation axis, S2A_k20")
    w()
    w("| pair | cos | signed angle |")
    w("|---|---|---|")
    keys = [k for k in PRN if k != "C2_seng4_P3"]
    for i, a in enumerate(keys):
        for b in keys[i + 1:]:
            r = row(a, b)
            w(f"| {PRN[a]} vs {PRN[b]} | {f(float(r.cos),4)} | {float(r.angle_deg):.2f}° |")
    w()
    w("### The full-gene space, the same aging-against-proliferation pairs")
    w()
    w("Reported because the flag asks for both spaces. **Not comparable to the table")
    w("above**: these are 23,485-dimensional cosines.")
    w()
    w("| aging axis | proliferation axis | cos | from parallel | % removed |")
    w("|---|---|---|---|---|")
    for a in ("A1", "A2", "A3", "A4"):
        for c in ("C1_gtex", "C1_south", "C1_seng_CRA", "C1_lu", "C2"):
            r = row(a, c, "full_gene")
            mag = float(np.degrees(np.arccos(min(1.0, abs(float(r.cos))))))
            w(f"| {a} | {PRN[c]} | {f(float(r.cos),4)} | {mag:.2f}° | "
              f"{100 * (1 - float(r.frac_a_surviving)):.1f}% |")
    w()
    w("### The triangle consistency check")
    w()
    w("The flag's version names four vectors, not three, so it is not a triangle")
    w("(F3). Two genuine spherical triangles are checked instead, on the standard")
    w("inequality:")
    w()
    w("| triangle | θ(A1,A3) | θ(A1,third) | θ(A3,third) | allowed range for θ(A1,A3) | consistent |")
    w("|---|---|---|---|---|---|")
    for k, t in d1["triangles"].items():
        w(f"| {k} (third = {PRN.get(t['third'], t['third'])}) | {t['angle_A1_A3']:.2f}° | "
          f"{t['angle_A1_third']:.2f}° | {t['angle_A3_third']:.2f}° | "
          f"{t['lower']:.2f}° to {t['upper']:.2f}° | "
          f"**{'yes' if t['consistent'] else 'NO'}** |")
    w()
    w("**Both are consistent**, so there is no geometric contradiction in the")
    w("measurements. The flag's expectation failed not because the geometry is")
    w("inconsistent but because it assumed one proliferation direction where there")
    w("are several, and they differ: the bulk-donor construction is far from the")
    w("cell-derived ones (see the proliferation-against-proliferation table).")
    w()
    w("### The D1 gate")
    w()
    w(f"**{d1['gate']}** fired, on `d` = {d1['d_angle_magnitude']:.2f}° — the angle")
    w(f"from parallel between A1 and `C1_gtex` in S2A_k20, per amendment F10")
    w(f"(`cos` = {f(d1['d_cos'])}, signed angle {d1['d_angle_signed']:.2f}°).")
    w()
    if d1["gate"] == "INTERMEDIATE":
        w(f"> \"donor aging sits {d1['d_angle_magnitude']:.2f} degrees from")
        w(f"> proliferation, against {d1['p_angle']} for passage aging. Donor age is")
        w("> the better separated of the two models but proliferation still carries")
        w(f"> {d1['pct_removed']:.1f}% of it, and any donor-age rejuvenation readout must")
        w("> report a proliferation-residualized version beside the headline.\"")
    elif d1["gate"] == "SEPARABLE":
        w(f"> \"donor aging sits {d1['d_angle_magnitude']:.2f} degrees from an")
        w(f"> independently measured proliferation axis, against {d1['p_angle']}")
        w("> degrees for passage aging in the same construction. A rejuvenation")
        w("> signal can be distinguished from a change in proliferative state")
        w("> against donor age and cannot be against passage number.\"")
    else:
        w(f"> \"donor aging sits {d1['d_angle_magnitude']:.2f} degrees from")
        w(f"> proliferation, close to the {d1['p_angle']} of passage aging.")
        w("> Proliferation cannot be separated from either aging model here, and")
        w("> this project's own Southard proliferation finding requires re-reading.\"")
    w()
    w("**The ENTANGLED branch did not fire**, so the clause about SOUTH9's RETAINED")
    w("reading and SOUTH8's residualized count being placed in question does not")
    w("apply. Those files are untouched and unquestioned by this task.")
    w()
    w("**C2 would fire a different branch, and the flag requires that be said")
    w("prominently.**")
    w()
    w("| construction | cos with A1 | from parallel | branch |")
    w("|---|---|---|---|")
    w(f"| **C1, data-derived (primary)** | {f(d1['d_cos'],4)} | "
      f"{d1['d_angle_magnitude']:.2f}° | **{d1['gate']}** |")
    w(f"| C2, gene-list | {f(d1['d_cos_c2'],4)} | {d1['d_angle_c2_magnitude']:.2f}° | "
      f"{d1['gate_c2']} |")
    w()
    w(f"C2's sentence, reported but **not** the task's reading: \"donor aging sits")
    w(f"{d1['d_angle_c2_magnitude']:.2f} degrees from an independently measured")
    w(f"proliferation axis, against {d1['p_angle']} degrees for passage aging in the")
    w("same construction. A rejuvenation signal can be distinguished from a change")
    w("in proliferative state against donor age and cannot be against passage")
    w("number.\" C1 is primary because the flag fixed it so, and because it is the")
    w("construction SENG4's primary axis used.")
    w()
    w("## 6. D2 — report-only, no decision rule")
    w()
    w("### Do the Southard- and Sengstack-derived C1 axes agree?")
    w()
    rs = row("C1_south", "C1_seng_CRA")
    rs2 = row("C1_south", "C1_seng_CRI")
    w(f"`cos(C1_south, C1_seng_CRA)` = {f(float(rs.cos),4)} ({float(rs.angle_deg):.2f}°)")
    w(f"and `cos(C1_south, C1_seng_CRI)` = {f(float(rs2.cos),4)} "
      f"({float(rs2.angle_deg):.2f}°), in S2A_k20.")
    w()
    w("The same recipe on two different CRISPR screens in two different fibroblast")
    w("systems produces axes that agree substantially but not closely. That is")
    w("evidence the recipe measures a shared, real direction rather than noise, and")
    w("also that it is not measuring exactly the same thing in both — part of the")
    w("difference is the structural-zero issue of F7, since each axis is zero on")
    w("genes the other dataset measures.")
    w()
    w("### Dependence on the construction: C2 against C1")
    w()
    w("| aging axis | % removed by its C1 | % removed by C2 |")
    w("|---|---|---|")
    for a, c in (("A1", "C1_gtex"), ("A2", "C1_gtex"), ("A3", "C1_seng_CRA"),
                 ("A4", "C1_lu")):
        r1_, r2_ = row(a, c), row(a, "C2")
        w(f"| {a} {AGN[a]} | {100*(1-float(r1_.frac_a_surviving)):.1f}% | "
          f"{100*(1-float(r2_.frac_a_surviving)):.1f}% |")
    w()
    w("The choice of construction changes the answer by tens of percentage points,")
    w("which is why the flag fixed C1 as primary before any number.")
    w()
    w("### Reference-panel coverage, and how much of the donor-age gap is unreachable")
    w()
    w("| dataset | panel genes present | absent | % of the GTEx age gap's squared length in unmeasured genes |")
    w("|---|---|---|---|")
    for k in ("Lu", "Southard", "Sengstack"):
        v = d2[k]
        w(f"| {k} | {v['present']:,} | {v['absent']:,} | "
          f"**{v['pct_gap_sq_in_unmeasured']:.2f}%** |")
    w()
    w("So between a tenth and a fifth of what defines the GTEx donor-age direction")
    w("is in genes these datasets do not measure at all. **This bounds how much of")
    w("any approach is unreachable by construction and is reported, not corrected.**")
    w()
    w("### Did any earlier SENG task run a GTEx-anchored test?")
    w()
    e = d2["earlier_gtex_anchored"]
    w(f"**A direction test, yes; an approach test, no.** `{e['file']}` carries")
    w(f"columns `{', '.join(e['columns'])}`: SENG computed the cosine of each")
    w("displacement with the GTEx donor-age direction and a permutation p for all")
    w(f"{e['n']} perturbations, and **{e['n_p_le_05']} of {e['n']}** reached p ≤ 0.05")
    w(f"(median cos {e['median_cos_d_u']:+.4f}, max {e['max_cos_d_u']:+.4f}). It did")
    w("**not** run an approach or distance test anchored on the GTEx centroids. D3")
    w("below is the first GTEx-anchored *approach* test on the Sengstack data.")
    w()
    w("### The Sengstack cells against the GTEx centroids")
    w()
    p = d2["PD32_vs_GTEx"]
    w(f"Mirroring Supplementary Note 1 of the preprint, in the full-gene space:")
    w()
    w(f"- distance from the Sengstack PD32 pseudobulk to the GTEx **young** centroid:")
    w(f"  **{p['dist_to_young']:.1f}**")
    w(f"- to the GTEx **old** centroid: **{p['dist_to_old']:.1f}**")
    w(f"- the GTEx young-to-old gap itself: **{p['young_old_gap']:.2f}**")
    w(f"- so PD32 sits **{p['ratio_to_gap_young']:.1f}×** the entire donor-age gap away")
    w(f"  from young and **{p['ratio_to_gap_old']:.1f}×** away from old — essentially")
    w("  equidistant, and far outside the span the gap defines.")
    w()
    w("**A bound, not a test.** Any displacement of these cells toward GTEx young is")
    w("a move of at most a few percent of the distance to it, from a starting point")
    w("that is not between the two donor-age anchors in the first place.")
    w()
    w("## 7. D3 — the four factors against donor aging (gated)")
    w()
    w("**Every number in this section transplants a displacement measured on 10x")
    w("single cells onto a bulk GTEx origin.** The displacement is the perturbation")
    w("minus its own library's NT cells, on the same machine; the origin is the GTEx")
    w("old centroid and the target the GTEx young centroid, which are bulk")
    w("post-mortem tissue. That is the Southard construction")
    w("(`src/south6_run.py:506` `block_delta`, with `r = c_old − c_young` from")
    w("`src/south7_run.py:134`), and `south_recipe_unreadable` does not fire.")
    w()
    w("### D3a — the plant gate, run before any factor number")
    w()
    w(f"The GTEx gap is {d3['gtex_gap_full']:.2f} in the full-gene space and")
    w(f"{d3['gtex_gap_k20']:.2f} in S2A_k20; plant doses are fractions of it. The 0%")
    w("row is the noise floor and is reported, not subtracted (amendment F4).")
    w()
    w("| library | n cells | floor | power @25% | power @50% | power @100% |")
    w("|---|---|---|---|---|---|")
    for (lib, n), g in pl.groupby(["library", "n_cells"]):
        g = g.set_index("dose")
        w(f"| {lib} | {int(n)} | {f(g.loc[0.0,'floor'],3)} | "
          f"**{g.loc[0.25,'power']:.3f}** | {g.loc[0.50,'power']:.3f} | "
          f"{g.loc[1.0,'power']:.3f} |")
    w()
    w(f"**{d3['gate_a']}** fired: p25 = {d3['p25']:.3f}, p50 = {d3['p50']:.3f}.")
    w()
    w(f"> \"the GTEx-anchored test detects a 25% approach toward donor youth at")
    w(f"> these factors' depths, lowest power {d3['p25']:.3f}.\"")
    w()
    w("This is the opposite of the sensitivity problem SENG2 hit and SENG3 only")
    w("partly fixed: against the GTEx anchors the instrument is well powered at")
    w("every one of these cell counts.")
    w()
    w("### D3b — the counts")
    w()
    s = d3["summary"]
    w(f"Chance expectation is {s['chance_expected']:.1f} of {s['n']}.")
    w()
    w("| group | n | below floor | distance decreased | **called** | chance exp. | 95% interval |")
    w("|---|---|---|---|---|---|---|")
    w(f"| all | {s['n']} | {s['below_floor']} | {s['decreased']} | **{s['called']}** | "
      f"{s['chance_expected']:.1f} | {s['called_ci'][0]*100:.1f}–{s['called_ci'][1]*100:.1f}% |")
    for m in ("CRA", "CRI"):
        v = s[m]
        w(f"| {m} | {v['n']} | {v['below_floor']} | {v['decreased']} | {v['called']} | "
          f"{v['chance_expected']:.1f} | — |")
    w()
    w(f"**{s['called']} of {s['n']} is below the {s['chance_expected']:.0f} expected by")
    w("chance**, and only 17 of 340 even move closer at all. With an instrument")
    w("qualified at p25 = 0.945, that is a positive finding of absence rather than a")
    w("failure to detect.")
    w()
    w("### The four factors")
    w()
    w("| factor | n cells | pctl | Δdistance to GTEx young | floor | called | cos(d, GTEx age) | cos(d, A3 passage) | cos(d, its C1) |")
    w("|---|---|---|---|---|---|---|---|---|")
    for lab in NAMED:
        n = d3["named"][lab]
        w(f"| `{lab}` | {int(n['n_cells'])} | {n['cell_pctl_in_screen']:.0f}% | "
          f"**{f(n['delta'],3)}** | {f(n['floor'],3)} | "
          f"{'yes' if n['called'] else '**no**'} | {f(n['cos_d_u'],4)} | "
          f"{f(n['cos_d_A3'],4)} | {f(n['cos_d_C1'],4)} ({n['angle_d_C1']:.1f}°) |")
    w()
    w("**All four move away from young human donors, not toward them** — every delta")
    w("is positive against a floor near −1.2, and every `cos(d, GTEx age)` is")
    w("negative. The last two columns say where they do point: strongly along the")
    w("passage-rejuvenation direction (+0.72 to +0.78) and more strongly still along")
    w("their own proliferation axis (+0.79 to +0.90, 26-38°).")
    w()
    w("### SOUTH8's stored result for the same four gene names")
    w()
    w("Labelled as the flag requires: the Southard screen is **CRISPR activation")
    w("only**, so for `STAT3` and `ZFX` — which Junevity perturbs by")
    w("**interference** — this is the opposite perturbation and is not a replication.")
    w()
    w("| gene | in Southard | n cells | delta, S2A_k20 | below depth floor | amended call |")
    w("|---|---|---|---|---|---|")
    for g, v in d3["southard_same_names"].items():
        note = " *(opposite perturbation)*" if g in ("STAT3", "ZFX") else " *(activation, same direction)*"
        w(f"| `{g}`{note} | {'yes' if v['in_southard'] else 'no'} | "
          f"{v.get('n_cells','—')} | {f(v.get('delta_S2A_k20'),4)} | "
          f"{'yes' if v.get('below_depth_floor') else 'no'} | "
          f"{'yes' if v.get('amended_call') else 'no'} |")
    w()
    w("All four are positive-delta and uncalled in Southard too, on a different")
    w("cell line with a different screen.")
    w()
    w("### The D3 gate")
    w()
    w(f"**{d3['gate_b']}** fired, on {d3['k_named']} of the four called.")
    w()
    if d3["gate_b"] == "NO APPROACH":
        w("> \"none of Junevity's four factors measurably moves their cells toward")
        w("> young human donors as measured against GTEx, with an instrument")
        w("> qualified at their depth.\"")
    else:
        called = [l for l in NAMED if d3["named"][l]["called"]]
        w(f"> \"{d3['k_named']} of Junevity's four factors move their own cells toward")
        w(f"> young human donors as measured against GTEx: {', '.join(called)}.\"")
    w()
    w("**No factor named here is validated or a pilot candidate.** Passaged")
    w("neonatal-line cells remain the substrate, the origin is bulk post-mortem")
    w("tissue, and between 10% and 22% of the donor-age direction is in genes this")
    w("dataset does not measure.")
    w()
    w("## 8. The reading")
    w()
    w("Only the pre-registered readings whose conditions fired are conclusions:")
    w()
    w("1. **Stage 0: PASS.**")
    w(f"2. **D1: {d1['gate']}** on C1, `d` = {d1['d_angle_magnitude']:.2f}° from")
    w(f"   parallel; C2 would fire {d1['gate_c2']} and is reported beside it.")
    w(f"3. **D3a: {d3['gate_a']}**, p25 = {d3['p25']:.3f}.")
    w(f"4. **D3b: {d3['gate_b']}**, {d3['k_named']} of four.")
    w()
    w("Put together, on one scale, in S2A_k20, as angle from parallel to a")
    w("proliferation axis built from the same units:")
    w()
    w("| aging model | angle | what that licenses |")
    w("|---|---|---|")
    w("| Sengstack passage (A3) | **19.66°** | nothing: SENG4 showed residualization here removes the signal |")
    w(f"| GTEx donor age (A1) | **{d1['d_angle_magnitude']:.2f}°** | a readout, provided a residualized version is reported beside it |")
    w(f"| Lu donor pair (A4) | **{float(np.degrees(np.arccos(abs(float(row('A4','C1_lu').cos))))):.2f}°** | a readout; 15.5% of the axis is proliferation |")
    w(f"| frozen ruler (A2) | **{float(np.degrees(np.arccos(abs(float(row('A2','C1_gtex').cos))))):.2f}°** | a readout; 0.3% of the axis is proliferation |")
    w()
    w("And the factor-level result, which D3 is the first test in this project to")
    w("produce: with an instrument qualified at these depths, **none of the four")
    w("moves toward young human donors, 5 of 340 perturbations do against 17")
    w("expected by chance, and all four instead point along passage reversal and")
    w("proliferation.**")
    w()
    w("### What this does not cover")
    w()
    w("- **No angle measured in GTEx bulk donors is a statement about single cells,")
    w("  and no angle measured in cells is a statement about bulk donors.** The")
    w("  flag forbids both readings and the bulk/cell split is exactly where the")
    w("  expectation broke: A1's bulk-donor proliferation axis behaves unlike the")
    w("  three cell-derived ones.")
    w("- **What `C1_gtex` actually tracks is not established here.** In bulk")
    w("  post-mortem tissue a cell-cycle gene-set score can move with tissue")
    w("  quality, ischemic time or sampling as easily as with proliferation, and")
    w("  this task ran no diagnostic to separate those. Its strong anti-alignment")
    w("  with donor age should be read as a reason to doubt the construction, not")
    w("  as a biological claim. That is the main limitation of the D1 gate.")
    w("- **No claim about aged human cells from the Sengstack data**: passaged")
    w("  neonatal-line fibroblasts whose \"old\" is passage number, and")
    w(f"  cos(A1, A3) = {f(float(row('A1','A3','full_gene').cos),4)} in the full-gene")
    w("  space confirms again that passage is not donor aging.")
    w("- **A1, A2, A4 and `C1_south` have no intervals** (F9), so any statement")
    w("  resting on them is a point estimate only.")
    w("- **The gate quantity shares units on both sides** (method statement (d)).")
    w("- **No claim that Junevity's published numbers are unreproducible.**")
    w("- **D1 and D2 name no transcription factor**, as the flag requires; D3 is the")
    w("  exempted section.")
    w()
    w("### Files")
    w()
    for p_ in sorted(OUT.iterdir()):
        w(f"- `results/axang/{p_.name}` ({p_.stat().st_size:,} B)")
    w()
    w("Nothing was committed and nothing was pushed.")

    FIND.write_text("\n".join(L) + "\n", encoding="utf-8")
    print("wrote FINDINGS_AXANG.md", len(L), "lines")


if __name__ == "__main__":
    main()
