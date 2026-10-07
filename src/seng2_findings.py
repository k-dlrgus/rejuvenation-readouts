"""Writes FINDINGS_SENG2.md from the stage outputs in results/seng2/.

Only the pre-registered reading whose conditions fired is reported as a
conclusion. Everything else is labelled diagnostic or report-only.
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
ROOT = SRC.parent
OUT = ROOT / "results" / "seng2"
FIND = ROOT / "FINDINGS_SENG2.md"
NAMED = ["CRA_E2F3", "CRA_EZH2", "CRI_STAT3", "CRI_ZFX"]
L = []


def w(s=""):
    L.append(s)


def jl(name):
    return json.loads((OUT / name).read_text(encoding="utf-8"))


def f(x, nd=6):
    if x is None:
        return "n/a"
    try:
        xf = float(x)
    except (TypeError, ValueError):
        return str(x)
    if not np.isfinite(xf):
        return "n/a"
    return f"{xf:+.{nd}f}" if nd and abs(xf) < 1e6 else f"{xf:.{nd}f}"


def main():
    s0 = jl("stage0_summary.json")
    d1 = jl("d1_summary.json")
    d2 = jl("d2_summary.json")
    d3 = jl("d3_summary.json")
    d4 = jl("d4_summary.json")
    audit = jl("raw_audit.json")
    med = jl("panel_medians.json")
    repro = pd.read_csv(OUT / "stage0_reproduction.csv")
    d2t = pd.read_csv(OUT / "d2_table.csv")
    ms = pd.read_csv(OUT / "d3_mixture_summary.csv")
    ps = pd.read_csv(OUT / "d3_plant_summary.csv")

    w("# FINDINGS_SENG2 — do Junevity's four factors approach young, once the")
    w("# library structure is accounted for?")
    w()
    w("Task SENG2. Computational, no new data. Reads `results/seng/`,")
    w("`results/seng_why/`, `results/senginv/` and the Sengstack h5ad already on")
    w("disk. Writes `results/seng2/` and this file. Nothing committed, nothing")
    w("pushed. No earlier findings file or results directory was written.")
    w()
    w("**Headline, stated once and then supported below.** The gate that governs")
    w("every factor-level number in this file is D1, and it fired")
    w(f"**{d1['gate']}**. D3 fired **{d3['gate']}**. D2's own gate fired")
    w(f"**{d2['gate']}**. What those three together license is set out in")
    w("\"The reading\" at the end, and nowhere else.")
    w()

    # ---------------- pre-registration, verbatim -------------------------
    w("---")
    w()
    w("## 1. The pre-registration, verbatim")
    w()
    w("Reproduced verbatim from `results/seng2/PREREG.flag`, written")
    w("2026-10-06T06:01:09Z before any SENG2 number existed and not edited since.")
    w()
    w("```")
    for line in (OUT / "PREREG.flag").read_text(encoding="utf-8").splitlines():
        w(line)
    w("```")
    w()

    # ---------------- amendments -----------------------------------------
    w("## 2. Amendments, and every deviation")
    w()
    w("The flag was found ambiguous or self-contradictory in four places. Those")
    w("four were put to the author **before computing**, not resolved by the")
    w("session; five further definitions the flag does not supply at all were")
    w("fixed in writing before computing. All of it is in")
    w("`results/seng2/AMENDMENTS.md`, reproduced here in full.")
    w()
    w("```")
    for line in (OUT / "AMENDMENTS.md").read_text(encoding="utf-8").splitlines():
        w(line)
    w("```")
    w()
    w("### Deviations beyond the amendments")
    w()
    w("1. **The bootstrap confidence intervals of SENG's screen were not")
    w("   re-derived.** Stage 0 checks no CI, and neither `meets_toward` nor")
    w("   `n_cos_q` reads one: both depend only on the nulls and the point")
    w("   estimate. The `*_ci_*` columns of `results/seng2/stage0_screen.csv` are")
    w("   therefore all NaN. This was a deliberate scope reduction taken after the")
    w("   first attempt showed the bootstrap would take about four hours and feed")
    w("   nothing that is checked. SENG's own CIs stand unexamined by this task.")
    w("2. **`.X` is streamed, not loaded.** About 2.6 GB of RAM was free against a")
    w("   1.4 GB matrix, so group pseudobulks are accumulated in 2,000-row blocks")
    w("   and only the few thousand cells needing per-cell resampling are")
    w("   materialised. This changes no number; the panel reproduces SENG's.")
    w("3. **Amendment A10, fixed in writing before computing.** Groups that are")
    w("   not in SENG's 348-row panel (the full PD14 and PD32 arms, the half-B")
    w("   spike material, every mixture and plant draw) enter the frozen space")
    w("   through `S.z_from_counts` against the panel's frozen TMM cache. The")
    w("   panel is **not** extended and nothing is refit. This is SENG's own")
    w("   mechanism for every null and bootstrap draw it computed, and SENG")
    w("   verified it agrees with the panel path to under 1e-4.")
    w("4. **This session reused SENG's implementation.** Stage 0 calls the frozen")
    w("   helper functions in `src/seng_run.py` (`z_from_counts`, `stats_from_d`,")
    w("   `bh_q`, `meets_toward`, `parse_wide_obs`, the TMM path) from a SENG2")
    w("   driver. A reproduction that reuses the original implementation tests the")
    w("   record, the data path and the driver — **not** the implementation. An")
    w("   independent reimplementation would be a stronger check and was not done.")
    w("5. The positive control was reimplemented in `src/seng2_run.py` rather than")
    w("   calling `S.run_posctrl`, because that function writes into")
    w("   `results/seng/`, which this task may not do.")
    w()

    # ---------------- method before use ----------------------------------
    w("## 3. State the method before using it")
    w()
    w("### (a) The coordinate space, and its circularity")
    w()
    w("The space is the **23,485-column frozen GTEx fibroblast ruler z space**, the")
    w("project's frozen instrument, read and not refit. A group of cells becomes a")
    w("point in it by: summing raw counts over the group's cells; scattering the sum")
    w("into the 23,485 ruler columns (16,781 matched, by gene symbol only — the file")
    w("carries no Ensembl IDs, so 0 matched by Ensembl — and the remaining 6,704 left")
    w("at z=0); one joint TMM over the analysis panel; log2-CPM with edgeR's prior;")
    w("and z against the frozen GTEx mu and sd.")
    w()
    w("The basis is not fitted at all — it is the frozen ruler's gene axes. What **is**")
    w("fitted on cells is:")
    w()
    w("- `v = z(EARLY_A) - z(LATE_A)`, SENG's aging axis, from PD14 half A (290")
    w("  cells) and PD32 half A (2,174 cells) — `src/seng_run.py:1836`.")
    w("- the proliferation axis `p` per modality, from the top and bottom NT")
    w("  proliferation quartiles.")
    w("- the joint TMM normalisation, fitted over **one panel containing EARLY_A,")
    w("  LATE_A, both NT groups, all four proliferation-quartile groups and all 340")
    w("  perturbations** — `src/seng_run.py:1795`.")
    w()
    w("**The circularity, stated plainly and not refitted away.** The anchor cells")
    w("that define the target direction also co-determine the normalisation that")
    w("every displacement is measured through, because they are rows of the same TMM")
    w("panel. Separately, D1's `a` and D2's target are built from PD14 and PD32 cells,")
    w("and D1's gate quantity is a cosine between two vectors that both contain")
    w("`z(PD14)` — so the gate is not independent of the anchor estimate. **Every")
    w("number in D1, D2, D3 and D4 of this file carries that circularity.** SENG")
    w("held PD14/PD32 **half B** out of the panel and used it for the positive")
    w("control, which is why half B is the spike material in D3 and never appears in")
    w("a target (amendment A4). Nothing was refit to remove any of this.")
    w()
    w("### (b) How SENG assigned a cell to a perturbation")
    w()
    w("Reproduced exactly by this task's code (`src/seng_run.py:648`,")
    w("`parse_wide_obs`), with the counts below re-derived here, not quoted:")
    w()
    w("- the pool is **804** guide columns: 798 targeting columns")
    w("  `<TARGET>_<CRA|CRI><k>` plus `NT_1`..`NT_6`. The six NT columns **are** in")
    w("  the argmax pool. SENGINV's 798 excluded them and is a different question;")
    w("  SENG's 804 is the rule that produced the stored numbers and is the rule")
    w("  used here.")
    w("- a non-WT cell is assigned the guide column that is a **unique strict")
    w("  maximum and > 0**. Pooled identity is `{CRA|CRI}_{TF}`.")
    w("- **2,680 cells have a tied positive maximum** and are given no winner.")
    w("- **270 cells** whose winning targeting column disagrees with the cell's own")
    w("  `CRA`/`CRI` indicator are **not assigned** — they are dropped to")
    w("  `unassigned`, not relabelled and not reassigned to the indicator's")
    w("  direction. 2,680 + 270 = 2,950 unassigned cells, 8.1% of the file.")
    w("- roles re-derived here: WT 6,809, CRA 12,264, CRI 12,346, NT 2,145,")
    w("  unassigned 2,950. Unique-argmax NT cells: CRA 666, CRI 1,479.")
    w("- a perturbation enters the panel at `n >= 30` cells: **340** qualify,")
    w("  167 CRA and 173 CRI, out of 400 target x direction pairs.")
    w()
    w("**Recorded as a limitation of the rule, not corrected.** The 2,950 dropped")
    w("cells are dropped non-randomly: a tie is likelier at low guide UMI, and the")
    w("270 disagreements are exactly the cells where the two encodings of direction")
    w("conflict — the conflict SENGINV already recorded as an ambiguity the file")
    w("does not resolve.")
    w()
    w("### (c) The documented source for NT being late-passage, same line as PD32")
    w()
    w("**STOP key `nt_passage_unsourced` does not fire.** A source can be given.")
    w()
    w("Sengstack et al., PNAS 2026, Results, section *\"Top TF Perturbations")
    w("Identified By TRDP Reversed Several Cellular Aging Hallmarks in Passaged")
    w("Fibroblasts\"*, defines the passage bands verbatim:")
    w()
    w("> \"Based on phenotyping we performed (data in Fig. 2), early passage cells")
    w("> were 0 to 20 population doublings (PD).\"")
    w()
    w("> \"middle passage cells were 21 to 30 PD, and late passage cells were over")
    w("> 31 PDs.\"")
    w()
    w("and names the controls in the same section: *\"Nontargeting sgRNA were used as")
    w("controls and are noted as either CRA NT (for CRISPRa cell line with")
    w("nontargeting sgRNA)\"* / *\"or CRI NT ... with nontargeting sgRNA\"*. The same")
    w("Results text states the perturbations were applied to late-passage cells")
    w("(*\"On average, the TF perturbations caused late passage cells to have cellular")
    w("phenotypes more like middle passage cells\"*). PD32 is over 31 PD, so it is")
    w("late passage by the paper's own cut.")
    w()
    w("**Two gaps in that source, carried into D1 rather than papered over.**")
    w()
    w("1. The paper gives a **band**, not a number. It never states the PD at which")
    w("   the Perturb-seq cells were harvested, only that they were over 31 PD. PD32")
    w("   sits at the very bottom of that band. **`offset_A` and `offset_I` are")
    w("   therefore upper bounds on library offset**: any real passage difference")
    w("   between the screen cells and PD32 is inside them and cannot be separated.")
    w("2. The same section says the change *\"occurred within about three population")
    w("   doublings (duration of the experiment)\"*, so the perturbed cells drifted in")
    w("   passage during the experiment. The h5ad carries no PD label for any")
    w("   perturbed cell (SENG: perturbed PD values empty, 24,610 cells blank), so")
    w("   none of this is checkable in the data.")
    w()
    w("Also: the CRISPRa and CRISPRi lines are engineered derivatives of the")
    w("fibroblast line, not the line itself, and the object gives them separate")
    w("barcode groups. \"Same line\" is the paper's framing, adopted here as the")
    w("flag instructs, and it is an assumption, not a measurement.")
    w()

    # ---------------- stage 0 ---------------------------------------------
    w("## 4. Stage 0 — reproduce first (gate)")
    w()
    w(f"`no_raw_counts` does not fire. Full scan of all {audit['n_seen']:,} stored")
    w(f"values of `X/data`: {audit['n_fractional']} fractional, {audit['n_negative']} negative,")
    w(f"{audit['n_nonfinite']} non-finite, min {audit['min']:g}, max {audit['max']:g},")
    w("verdict raw integer counts in float32. (This also corrects the SENGINV")
    w("inventory, which reported max 2,446 from a 10M-value sample; the full-scan")
    w("maximum is 2,960, as SENG already had. That file is not edited.)")
    w()
    w(f"**Gate: {'PASS' if s0['passed'] else 'FAIL'}.** "
      f"Tolerance as pre-registered: counts exact, numerics within 1e-4 x max(1, |stored|).")
    w()
    w("| item | kind | stored | re-derived | abs diff | ok |")
    w("|---|---|---|---|---|---|")
    for r in repro.itertuples():
        if r.kind == "count":
            w(f"| `{r.item}` | count | {int(r.stored)} | "
              f"{'n/a' if pd.isna(r.got) else int(r.got)} | "
              f"{'n/a' if pd.isna(r.got) else abs(int(r.got) - int(r.stored))} | "
              f"{'yes' if r.ok else '**NO**'} |")
        else:
            w(f"| `{r.item}` | numeric | {r.stored:.16g} | {r.got:.16g} | "
              f"{abs(r.got - r.stored):.3e} | {'yes' if r.ok else '**NO**'} |")
    w()
    w(f"- `v` re-derived from the data against SENG's stored `axes.npz` `v`:")
    w(f"  max abs difference **{s0['v_max_abs_diff_vs_stored']:.3e}** over 23,485")
    w("  components. The axis is the same vector, not merely the same summary.")
    w(f"- cos(v, u) permutation p re-derived: {s0['p_vu']:.6g} (stored 1.0).")
    w("- the positive control reproduced to every printed digit of all four stored")
    w("  numbers, and the screen counts and all four named-hit cell counts are exact.")
    w()
    if not s0["passed"]:
        bad = [c["item"] for c in s0["checks"] if not c["ok"]]
        w(f"**STOP `prior_results_not_reproduced`.** Mismatches: {bad}.")
        w()
        FIND.write_text("\n".join(L) + "\n", encoding="utf-8")
        return

    # ---------------- D1 --------------------------------------------------
    w("## 5. D1 — how much of the target direction is library offset? (gated)")
    w()
    w("`a = centroid(PD14, 581) - centroid(PD32, 4,348)`; `offset_L = centroid(NT_L)")
    w("- centroid(PD32)`; the gate quantity is the cosine between `(PD14 - NT_L)`")
    w("and `a`. Depth matching is amendment A3: binomial thinning to the lowest")
    w("group median, 200 draws, seed 20260914, thinned primary.")
    w()
    w("Group median `n_counts`, which set the thinning factors:")
    w()
    w("| group | n cells | median n_counts |")
    w("|---|---|---|")
    for g, lbl in (("PD14_FULL", "PD14 (all)"), ("PD32_FULL", "PD32 (all)"),
                   ("NT_CRA", "NT_A"), ("NT_CRI", "NT_I")):
        w(f"| {lbl} | {d1['sizes'][g]} | {d1['medians'][g]:,.1f} |")
    w(f"| EARLY_A (PD14 half A) | 290 | {med['EARLY_A']:,.1f} |")
    w(f"| LATE_A (PD32 half A) | 2174 | {med['LATE_A']:,.1f} |")
    w()
    for Lk in ("A", "I"):
        w(f"- comparison {Lk}: target median {d1['targets'][Lk]:,.1f}, factors "
          + ", ".join(f"`{k}`={v:.4f}" for k, v in d1["factors"][Lk].items()))
    w()
    w("### Results, raw and depth-matched")
    w()
    w("| quantity | raw | depth-matched median | 2.5–97.5% |")
    w("|---|---|---|---|")
    t = d1["thinned"]
    rows = [
        ("\\|\\|a\\|\\|", d1["raw"]["a_norm"], t["a_norm_A_median"], (t["a_norm_A_lo"], t["a_norm_A_hi"])),
        ("\\|\\|offset_A\\|\\| / \\|\\|a\\|\\|", d1["raw"]["offset_A_over_a"],
         t["offset_over_a_A_median"], (t["offset_over_a_A_lo"], t["offset_over_a_A_hi"])),
        ("\\|\\|offset_I\\|\\| / \\|\\|a\\|\\|", d1["raw"]["offset_I_over_a"],
         t["offset_over_a_I_median"], (t["offset_over_a_I_lo"], t["offset_over_a_I_hi"])),
        ("cos(offset_A, a)", d1["raw"]["cos_offset_A_a"],
         t["cos_offset_a_A_median"], (t["cos_offset_a_A_lo"], t["cos_offset_a_A_hi"])),
        ("cos(offset_I, a)", d1["raw"]["cos_offset_I_a"],
         t["cos_offset_a_I_median"], (t["cos_offset_a_I_lo"], t["cos_offset_a_I_hi"])),
        ("**cos(PD14 - NT_A, a)**", d1["raw"]["cos_t_A_a"],
         t["cos_t_a_A_median"], (t["cos_t_a_A_lo"], t["cos_t_a_A_hi"])),
        ("**cos(PD14 - NT_I, a)**", d1["raw"]["cos_t_I_a"],
         t["cos_t_a_I_median"], (t["cos_t_a_I_lo"], t["cos_t_a_I_hi"])),
    ]
    for nm, rv, tv, (lo, hi) in rows:
        w(f"| {nm} | {f(rv)} | {f(tv)} | {f(lo)} to {f(hi)} |")
    w()
    w("Angles, the same numbers read as rotations:")
    w()
    w("| | raw | depth-matched |")
    w("|---|---|---|")
    for Lk in ("A", "I"):
        ra = np.degrees(np.arccos(np.clip(d1["raw"][f"cos_t_{Lk}_a"], -1, 1)))
        ta = np.degrees(np.arccos(np.clip(t[f"cos_t_a_{Lk}_median"], -1, 1)))
        w(f"| angle between (PD14 - NT_{Lk}) and a | {ra:.2f}° | {ta:.2f}° |")
    w()
    w(f"- raw `r` = {f(d1['r_raw'])} (gate on raw alone would be **{d1['gate_raw']}**)")
    w(f"- **depth-matched `r` = {f(d1['r_primary'])}**, from library `{d1['r_library']}`,")
    w(f"  2.5–97.5% [{f(d1['r_ci'][0])}, {f(d1['r_ci'][1])}], angle "
      f"{d1['r_angle_deg']:.2f}°")
    w(f"- the two agree on the branch, so A3's \"where they disagree\" clause is not")
    w("  invoked." if d1["gate"] == d1["gate_raw"] else
      f"- **the raw and depth-matched readings disagree on the branch.** Per A3 the")
    w(f"- cos(a, v) = {f(d1['cos_a_v'])}: the full-arm axis `a` against SENG's")
    w("  half-A axis `v`. They are not the same vector, and D1 is scaled against `a`")
    w("  as the flag defines it.")
    w()
    w("**Why thinning moves these numbers so little, stated so the reader is not")
    w("misled.** log2-CPM already divides out library size, and the prior count")
    w("scales with the library, so thinning a pseudobulk by a common factor changes")
    w("the z vector mainly by adding sampling noise rather than by removing a depth")
    w("effect. A3's thinning is therefore a weak control for the depth confound, not")
    w("a strong one, and the raw and thinned numbers being close is not evidence that")
    w("depth does not matter. Per-cell median after thinning is in")
    w("`results/seng2/d1_thinning_check.json`.")
    w()
    w("### The two libraries are not in the same condition — read this before the gate")
    w()
    w("The flag fixes `r` as **the larger** of the two cosines, so the gate is")
    w("decided by whichever library is better aligned. The two are not close, and")
    w("they fall in different pre-registered bands:")
    w()
    w("| library | cos(PD14 - NT_L, a), depth-matched | angle | band that cosine sits in | "
      "named factors in it |")
    w("|---|---|---|---|---|")
    for Lk, nm, facs in (("A", "`_A` (CRISPRa)", "CRA_E2F3, CRA_EZH2"),
                         ("I", "`_I` (CRISPRi)", "CRI_STAT3, CRI_ZFX")):
        cv = t[f"cos_t_a_{Lk}_median"]
        band = "USABLE" if cv >= 0.70 else "DEGRADED" if cv >= 0.40 else "UNUSABLE"
        ang = np.degrees(np.arccos(np.clip(cv, -1, 1)))
        w(f"| {nm} | {f(cv)} | {ang:.2f}° | {band} | {facs} |")
    w()
    w("**This is a diagnostic, not a second gate** — the flag fixed one gate on the")
    w("maximum and that is what is reported below. But it bears directly on two of")
    w("the four named factors: `CRA_E2F3` and `CRA_EZH2` sit in the library whose")
    w(f"own cosine is {f(t['cos_t_a_A_median'])}, inside the DEGRADED band, where the")
    w("pre-registered wording would have said the target direction is rotated")
    w("substantially. The gate that fired was carried by the CRISPRi library.")
    w()
    w("The offsets themselves are large and nearly orthogonal to aging, which is how")
    w("both of those cosines stay as high as they are:")
    w()
    w(f"- `||offset_A||` is **{t['offset_over_a_A_median']:.2f}x** the whole aging gap")
    w(f"  `||a||`, and cos(offset_A, a) = {f(t['cos_offset_a_A_median'])}.")
    w(f"- `||offset_I||` is **{t['offset_over_a_I_median']:.2f}x** `||a||`, with")
    w(f"  cos(offset_I, a) = {f(t['cos_offset_a_I_median'])}.")
    w()
    w("So the between-library displacement is comparable to or larger than the entire")
    w("PD14-to-PD32 aging gap, and it points almost at right angles to it. A large")
    w("near-orthogonal offset shortens the apparent target direction's alignment with")
    w("aging without destroying it; that is the geometry the gate is measuring.")
    w()
    w("### The D1 gate")
    w()
    r = d1["r_primary"]
    w(f"**{d1['gate']}** fired, on the pre-registered thresholds, with")
    w(f"`r` = {f(r)}.")
    w()
    if d1["gate"] == "USABLE":
        w(f"> \"the library offset rotates the target direction by less than 45")
        w(f"> degrees (cos {f(r)}); the apparent target direction is dominated by")
        w("> aging and the approach test can be read.\"")
        w()
        w(f"The flag's \"less than 45 degrees\" is the flag's rounding: `r >= 0.70`")
        w(f"is 45.57°, and the measured angle is {d1['r_angle_deg']:.2f}°"
          f" (amendment A7).")
        w()
        w("**One qualification on the word \"dominated\", which the flag's wording")
        w("does not carry.** A cosine of")
        w(f"{f(r)} means the component of the apparent target direction orthogonal to")
        w(f"aging is still **{d1['pct_library']:.1f}%** of its length")
        w("(`sqrt(1 - r^2)`). Aging is the single largest component and the gate's")
        w("threshold is met, but a clear majority of that vector's length is not")
        w("aging. Every D2 number carries that.")
    elif d1["gate"] == "DEGRADED":
        w(f"> \"the library offset rotates the target direction substantially (cos")
        w(f"> {f(r)}). D2 is run and reported, and every number it produces is")
        w(f"> labelled as carrying a target direction that is")
        w(f"> {d1['pct_library']:.1f}% library difference.\"")
        w()
        w(f"The percentage is amendment A7: `sqrt(1 - r^2)`, the fraction of the")
        w("apparent target direction orthogonal to aging.")
    else:
        w(f"> \"the vector from a perturbation's control to the young target is")
        w(f"> mostly library difference, not aging (cos {f(r)}). This dataset")
        w("> cannot answer whether these factors approach young, and SENG's zero is")
        w("> a statement about its own geometry rather than about the factors.\"")
        w()
        w("**D2 and D3 are still run to the end and reported, and nothing in them is")
        w("a conclusion about any factor.** That is said again wherever a factor is")
        w("named.")
    w()

    # ---------------- D2 --------------------------------------------------
    unus = d1["gate"] == "UNUSABLE"
    carry = ("**Carries D1: " + (
        "the target direction is mostly library difference, so nothing below is a "
        "conclusion about any factor.**" if unus else
        f"the target direction is {d1['pct_library']:.1f}% library difference.**"
        if d1["gate"] == "DEGRADED" else "the target direction is dominated by aging.**"))
    w("## 6. D2 — the approach test, run two ways (gated)")
    w()
    w(carry)
    w()
    w("Origins, per amendment A6. Target is PD14 half A throughout (A4).")
    w("Displacement `d = z(g) - z(NT_L)` under both.")
    w()
    w("- **origin (i), as-measured**: `t_i = z(PD14 half A) - z(NT_L)`")
    w("- **origin (ii), offset-corrected**: `t_ii = z(PD14 half A) - z(PD32 half A)`")
    w("  — this **transplants the displacement onto a different origin within the")
    w("  same cell line and experiment**, and `t_ii` is identical to SENG's frozen")
    w("  `v`, so origin (ii) is literally SENG's geometry.")
    w()
    w("Floors are the Southard construction, confirmed identical to SENG's null")
    w("(`src/south7_run.py:330`): draw `n` control cells from the perturbation's own")
    w("library, displacement `z(draw) - z(pool minus draw)`, 200 draws, floor at the")
    w("5th percentile. The null depends only on `n`, so one null per distinct cell")
    w("count. A perturbation is **called approaching** when its delta is below its")
    w("own floor **and** negative.")
    w()
    w("### Counts over all 340 perturbations, with the chance baseline")
    w()
    w("The floor is the 5th percentile of each perturbation's own null, so **5% of")
    w("true nulls fall below it by construction**: the chance expectation is 17.0 of")
    w("340. A count near 17 is what a screen with no real effect produces.")
    w()
    w("| origin | below floor | distance decreased | **called approaching** | "
      "chance expectation | 95% interval on the called fraction |")
    w("|---|---|---|---|---|---|")
    for tag, nm in (("i", "(i) as-measured"), ("ii", "(ii) offset-corrected")):
        s = d2["summary"][f"delta_{tag}"]
        w(f"| {nm} | {s['below_floor']} | {s['decreased']} | **{s['called']}** | "
          f"{s['chance_expected']:.1f} | "
          f"{s['called_ci'][0] * 100:.1f}–{s['called_ci'][1] * 100:.1f}% |")
    w()
    w("Per modality:")
    w()
    w("| origin | modality | n | below floor | decreased | called | chance exp. |")
    w("|---|---|---|---|---|---|---|")
    for tag, nm in (("i", "(i)"), ("ii", "(ii)")):
        for m in ("CRA", "CRI"):
            s = d2["summary"][f"delta_{tag}"][m]
            w(f"| {nm} | {m} | {s['n']} | {s['below_floor']} | {s['decreased']} | "
              f"{s['called']} | {s['chance_expected']:.1f} |")
    w()
    w("**Which half of the criterion is binding, because it decides everything")
    w("below.** \"Called approaching\" is the conjunction of two conditions, and they")
    w("are not equally hard to meet here:")
    w()
    si = d2["summary"]["delta_i"]
    sii = d2["summary"]["delta_ii"]
    w(f"- under origin (i), **{si['below_floor']} of 340** perturbations fall below")
    w(f"  their own floor, but only **{si['decreased']}** have a negative delta.")
    w(f"- under origin (ii), **{sii['below_floor']}** fall below floor and")
    w(f"  **{sii['decreased']}** have a negative delta.")
    w()
    w("The floor is easy to clear and `delta < 0` is almost impossible, because at")
    w("these cell counts the null delta is **positive**: sampling noise alone moves a")
    w("group *away* from the target. D3 measures that pedestal directly — at 100")
    w("cells the median null delta is about +27 — and it is why the counts below are")
    w("near zero. A count near zero here is a statement about the pedestal at least")
    w("as much as about the perturbations, which is exactly what D3 then tests.")
    w()
    w("### The four named factors")
    w()
    w("Cell counts are SENG's own rule (b). Percentile is within the 340-perturbation")
    w("screen by cell count.")
    w()
    w("| factor | n cells | pctl | delta (i) | floor (i) | called (i) | "
      "delta (ii) | floor (ii) | called (ii) |")
    w("|---|---|---|---|---|---|---|---|---|")
    for lab in NAMED:
        n = d2["named"][lab]
        w(f"| `{lab}` | {int(n['n_cells'])} | {n['cell_pctl_in_screen']:.0f}% | "
          f"{f(n['delta_i'], 4)} | {f(n['floor_delta_i'], 4)} | "
          f"{'**yes**' if n['called_delta_i'] else 'no'} | "
          f"{f(n['delta_ii'], 4)} | {f(n['floor_delta_ii'], 4)} | "
          f"{'**yes**' if n['called_delta_ii'] else 'no'} |")
    w()
    w(f"- called under origin (i): **{d2['k_origin_i']} of 4**")
    w(f"- called under origin (ii): **{d2['k_origin_ii']} of 4**")
    w()
    w("### The D2 gate")
    w()
    w(f"**{d2['gate']}** fired, on amendment A1's mutually exclusive branches.")
    w()
    k_i, k_ii = d2["k_origin_i"], d2["k_origin_ii"]
    called_i = [l for l in NAMED if d2["named"][l]["called_delta_i"]]
    called_ii = [l for l in NAMED if d2["named"][l]["called_delta_ii"]]
    if d2["gate"] == "APPROACH":
        w(f"> \"of Junevity's four validated factors, {min(k_i, k_ii)} end closer to")
        w("> early-passage cells than their own controls do, in their own data, with")
        w(f"> their own perturbations: {', '.join(sorted(set(called_i) & set(called_ii)))}.\"")
    elif d2["gate"] == "NO APPROACH":
        w("> \"none of Junevity's four validated factors measurably approaches")
        w("> early-passage cells in their own data, with their own perturbations,")
        w("> against a floor built at each factor's own depth.\"")
    else:
        w("> SPLIT. Both origins are reported and no single conclusion is drawn.")
        w(f"> Under origin (i) {k_i} of four are called"
          f"{': ' + ', '.join(called_i) if called_i else ''}.")
        w(f"> Under origin (ii) {k_ii} of four are called"
          f"{': ' + ', '.join(called_ii) if called_ii else ''}.")
    w()
    if unus:
        w("**This sentence is not a conclusion about any factor.** D1 fired UNUSABLE,")
        w("so the target direction these deltas are measured against is mostly")
        w("library difference. The gate wording is reported because it was")
        w("pre-registered, not because it licenses a claim.")
        w()

    # ---------------- D3 --------------------------------------------------
    w("## 7. D3 — could the statistic have seen an approach across the gap? (gated)")
    w()
    w("SENG's positive control ran entirely inside the `-1` barcode group and so")
    w("never crossed a library boundary. Spike material is PD14 **half B** (291")
    w("cells), which appears in no target.")
    w()
    w(f"`||a||` = {d3['a_norm']:.4f}, the full aging gap the plant doses scale.")
    w()
    w("### (a) Mixture across the library boundary")
    w()
    w("Half B mixed into NT cells at 0/10/25/50% at fixed total cell count (160 and")
    w("100, amendment A8), 200 draws per dose. The zero-dose point is the noise")
    w("floor and is reported, not subtracted.")
    w()
    w("| n total | library | dose | progress (frac, origin i) | 2.5–97.5% | "
      "delta (origin i) | 2.5–97.5% |")
    w("|---|---|---|---|---|---|---|")
    for r in ms.itertuples():
        w(f"| {int(r.n_total)} | {r.modality} | {r.dose:.0%} | {f(r.frac_i_median, 4)} | "
          f"{f(r.frac_i_lo, 4)} to {f(r.frac_i_hi, 4)} | {f(r.delta_i_median, 3)} | "
          f"{f(r.delta_i_lo, 3)} to {f(r.delta_i_hi, 3)} |")
    w()
    w("Monotonicity of progress, and whether the 50% dose's interval excludes zero:")
    w()
    w("| key | progress monotonic | delta(50%) interval excludes 0 |")
    w("|---|---|---|")
    for k in sorted(d3["monotonic"]):
        w(f"| {k} | {'yes' if d3['monotonic'][k] else '**no**'} | "
          f"{'yes' if d3['ci50_excludes_zero'][k] else '**no**'} |")
    w()
    w("### (b) Plant recovery along `a`")
    w()
    w("NT draws displaced along `a` by 0/25/50/100% of `||a||`, 200 draws per dose,")
    w("at 100 and 160 cells. Power is the fraction of draws D2 would have called")
    w("(amendment A9).")
    w()
    w("| n cells | library | dose | delta median (i) | floor (i) | **power (i)** | "
      "power (ii) |")
    w("|---|---|---|---|---|---|---|")
    for r in ps.itertuples():
        w(f"| {int(r.n_cells)} | {r.modality} | {r.dose:.0%} | {f(r.delta_i_median, 3)} | "
          f"{f(r.floor_i, 3)} | **{r.power_i:.3f}** | {r.power_ii:.3f} |")
    w()
    w("### The D3 gate")
    w()
    gi = d3["gate_inputs"]
    w(f"**{d3['gate']}** fired. Inputs: mixture monotonic at 160 cells under origin")
    w(f"(i) = {gi['mono_160_origin_i']}; 50% interval excludes zero = "
      f"{gi['ci50_160_origin_i']}; minimum power at a 25% plant at 100 cells = "
      f"{gi['power25_100_origin_i']:.3f} against the 0.80 requirement.")
    w()
    if d3["gate"] == "SENSITIVE":
        w("> \"the statistic detects a known approach across the library boundary at")
        w(f"> the depth of these factors, power {gi['power25_100_origin_i']:.3f}, so a")
        w("> null in D2 is a null and not an insensitive assay.\"")
    else:
        w("> \"the statistic cannot reliably detect a 25% approach across the library")
        w(f"> boundary at 100 cells, power {gi['power25_100_origin_i']:.3f}. D2's")
        w("> counts cannot be read as evidence of absence, and that applies to SENG's")
        w("> zero as well.\"")
        w()
        w("**D2's reading is reported but carries no conclusion.**")
    w()

    # ---------------- D4 --------------------------------------------------
    w("## 8. D4 — report-only, no decision rule")
    w()
    w("Nothing in this section is gated and nothing in it decides anything.")
    w()
    w("### Proliferation")
    w()
    rp = d4["rho_prolif_stored"]
    w(f"SENG's stored rank correlation between the screen statistic and")
    w(f"proliferation is **rho = +{rp['CRA']:.3f} (CRA)** and **+{rp['CRI']:.3f}")
    w("(CRI)**. It is quoted beside every result in this file, as the flag requires.")
    w("A proliferation-residualized D2, with the null residualized identically by")
    w("projecting out the same per-modality proliferation axis:")
    w()
    w("| origin | called, plain delta | called, residualized | called by both |")
    w("|---|---|---|---|")
    for tag in ("i", "ii"):
        r = d4[f"resid_origin_{tag}"]
        w(f"| ({tag}) | {r['called_plain']} | {r['called_resid']} | {r['called_both']} |")
    w()
    w("Reported beside the distance result, never in place of it.")
    w()
    w("### Fibroblast identity")
    w()
    idn = d4["identity"]
    w(f"Mean log2-CPM change over the {idn['n_genes']} identity genes SENG found")
    w(f"present ({', '.join(d4['identity_genes_used'])}), perturbation minus its NT:")
    w(f"screen median {f(idn['screen_median'], 4)}, 2.5–97.5%")
    w(f"[{f(idn['screen_lo'], 4)}, {f(idn['screen_hi'], 4)}]; "
      f"{idn['n_identity_loss']} of {d4['n']} exceed SENG's identity-loss threshold.")
    w()
    w("### The frozen ruler — the preprint's thesis in this dataset")
    w()
    for tag, nm in (("i", "(i) as-measured"), ("ii", "(ii) offset-corrected")):
        r = d4[f"ruler_origin_{tag}"]
        w(f"- origin {nm}: **{r['n_scores_younger']} of {d4['n']}** perturbations score")
        w(f"  younger on the one-dimensional frozen ruler, of which")
        w(f"  **{r['n_younger_and_not_called']} score younger while not approaching**;")
        w(f"  {r['n_younger_and_called']} both score younger and approach, and")
        w(f"  {r['n_not_younger_and_called']} approach without scoring younger.")
    w()
    w("This is a **diagnostic**, not evidence about any factor: it restates in this")
    w("dataset the dissociation the preprint is about. It inherits every")
    w("circularity in (a) and the D1 gate.")
    w()
    w("### The four factors, report-only columns")
    w()
    w("| factor | n cells | identity change | ruler change | younger on ruler | "
      "prolif. shift | SENG cos | SENG q_cos |")
    w("|---|---|---|---|---|---|---|---|")
    for lab in NAMED:
        n = d4["named"][lab]
        w(f"| `{lab}` | {n['n_cells']} | {f(n['identity_diff'], 4)} | "
          f"{f(n['ruler_change'], 4)} | {'yes' if n['ruler_scores_younger'] else 'no'} | "
          f"{f(n['prolif_shift'], 4)} | {f(n['seng_cos'], 4)} | {n['seng_q_cos']:.4f} |")
    w()
    w("### R_rej, restated post-hoc and unexplained")
    w()
    w("SENG_WHY recomputed the authors' own statistic and got roughly half the")
    w("published value for all four. Restated with its numbers, **labelled post-hoc")
    w("and unexplained. No claim is made about it in this task, and in particular no")
    w("claim that Junevity's published numbers are unreproducible.**")
    w()
    w("| factor | n cells | R_rej recomputed | published Table 1 | ratio |")
    w("|---|---|---|---|---|")
    for r in d4["rrej_seng_why"]:
        w(f"| `{r['perturbation']}` | {r['n_cells']} | {r['R_rej_recomputed']:.4f} | "
          f"{r['published_table1']:.2f} | {r['ratio']:.2f} |")
    w()
    w("SENG_WHY's own `assumptions` block records six choices made after numbers")
    w("existed (gene selection, Pearson, sign convention, the log2FC recipe, which")
    w("PDs count as aging, and the z-scoring of unmapped genes). The ratio above is")
    w("not corrected for any of them.")
    w()

    # ---------------- the reading ----------------------------------------
    w("## 9. The reading")
    w()
    w("Only the pre-registered readings whose conditions fired are reported as")
    w("conclusions. They are, in the order the flag gates them:")
    w()
    w(f"1. **Stage 0: PASS.** Every stored SENG number re-derived within tolerance,")
    w("   counts exact.")
    w(f"2. **D1: {d1['gate']}**, `r` = {f(d1['r_primary'])}.")
    w(f"3. **D2: {d2['gate']}**.")
    w(f"4. **D3: {d3['gate']}**.")
    w()
    if unus:
        w("Because D1 fired UNUSABLE, **the only conclusion this task supports is")
        w("D1's**: in this dataset the vector from a perturbation's own control to")
        w("the young target is mostly library difference rather than aging, so the")
        w("dataset cannot answer whether these factors approach young, and SENG's")
        w("zero is a statement about its own geometry rather than about the factors.")
        w("D2's and D3's sentences are reported above as pre-registered, and neither")
        w("is a conclusion about any transcription factor.")
    elif d3["gate"] == "INSENSITIVE":
        w("D1 licenses reading D2, but **D3 fired INSENSITIVE**, so D2's counts")
        w("cannot be read as evidence of absence. D2's sentence is reported above")
        w("and carries no conclusion.")
    else:
        w("D1 licenses reading D2 and D3 confirms the statistic could have seen the")
        w("effect, so D2's sentence above is the conclusion of this task.")
    w()
    w("### What this does not cover")
    w()
    w("- **Not a result about aged human cells.** This dataset is passaged neonatal")
    w("  dermal fibroblasts and its \"old\" is passage number. This project has")
    w("  already found that axis does not resemble donor aging: cos(v, u) =")
    w("  -0.040, permutation p = 1.0, re-derived in Stage 0. Nothing here transfers")
    w("  to aged donor cells, and no null here is evidence about transcription")
    w("  factors in aged human cells.")
    w("- **No factor named in this file is validated or a pilot candidate**, under")
    w("  any of the gates, including any factor D2 called.")
    w("- **No claim that Junevity's published numbers are unreproducible.** The")
    w("  R_rej comparison in D4 is post-hoc, unexplained, and explicitly not a claim.")
    w("- **The circularity in (a) is not removed anywhere.** The anchors co-determine")
    w("  the normalisation; D1's gate quantity shares `z(PD14)` with the target it")
    w("  judges. Nothing was refit to hide it.")
    w("- **`offset_A` and `offset_I` are upper bounds**, because the paper sources the")
    w("  screen cells only to a band (>31 PD) and not to PD32, and the h5ad carries")
    w("  no PD label for any perturbed cell.")
    w("- **Depth matching is weak.** See the note in D1: log2-CPM already divides out")
    w("  library size, so thinning mostly adds noise. The depth confound between the")
    w("  barcode groups is reported, not removed.")
    w("- **SENG's bootstrap CIs were not re-derived** (deviation 1).")
    w("- **Stage 0 reuses SENG's implementation** (deviation 4), so it tests the")
    w("  record and the data path, not the implementation.")
    w()
    w("### Files")
    w()
    for p in sorted(OUT.iterdir()):
        w(f"- `results/seng2/{p.name}` ({p.stat().st_size:,} B)")
    w()
    w("Nothing was committed and nothing was pushed.")

    FIND.write_text("\n".join(L) + "\n", encoding="utf-8")
    print("wrote FINDINGS_SENG2.md", len(L), "lines")


if __name__ == "__main__":
    main()
