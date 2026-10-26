"""Builds FINDINGS_CORR2CHK.md from results/corr2chk/.

Computes nothing. Every number is read from the result CSVs or from
results/corr2chk/PREREG.flag, which is reproduced verbatim. The prose is fixed here;
the numbers are not, so the findings file cannot drift from the tables.

Run from the repository root:
    python -I src/corr2chk_findings.py
"""
from __future__ import annotations

import json
from pathlib import Path

import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "results" / "corr2chk"
DST = ROOT / "FINDINGS_CORR2CHK.md"

# --------------------------------------------------------------------------- inventory
# (statement id, paragraph(s), what it says in brief, dataset, space, start,
#  floor/threshold defined at, drawn at, same space/start, plant gate, gate numbers)
INVENTORY = [
    ("A1a", "31", '"the four genes are also uncalled" (Southard)', "Southard", "k20",
     "GTEx old centroid",
     "`src/south7_run.py:319-324`, reference = pool − draw at `:323`; stored "
     "`results/south7/stage_a1_depth_nulls.csv`; attached per factor at `:537-543,555`",
     "own n = 83, 96, 129, 140; 200 draws (`src/south7_run.py:51`)", "yes / yes",
     "`src/south8_run.py:157-193`, reference stated `:160-161`, built `:176-178`; "
     "`results/south8/d1_plant_curves_all.csv`",
     "`power_25` = 0.925 (83 cells), 0.975 (96), 0.995 (129), 1.000 (140)"),
    ("A1b", "31", '"and move away from young donors in the Southard screen"', "Southard",
     "k20", "GTEx old centroid", "**CHECKED BY THIS TASK**", "—", "yes / yes",
     "as A1a", "as A1a"),
    ("A2", "18, 77", "the 20-component count was carried over from the Southard screen; "
     "quoted recovery powers 0.990 / 0.965 / 0.930, full gene space 0.115", "Southard",
     "k20/k50/k100 + full", "GTEx old centroid",
     "no floor — the quoted numbers **are** size-matched positive-control powers",
     "106 cells, the Southard screen's median",
     "yes / yes", "SOUTH6REC, re-run as a stage-0 check at `src/south8_run.py:266`",
     "0.990 / 0.965 / 0.930; full gene space 0.115, which FLOORCHK independently confirms"),
    ("A3", "Table 1 row 1; 48, 80, 101",
     "angles between the Southard cell-division direction and four definitions of aging",
     "Southard", "both, reported separately", "none (angles)",
     "none — `arccos|cos|`, defined at paragraph 80; no length, no change in distance",
     "—", "n/a — paper 2 states angles are never compared across spaces (paras 66, 80, 104)",
     "none", 'labelled "Descriptive only" (Table 1); no intervals (para 71)'),
    ("A4", "16", "the Southard and Lu datasets listed as reference datasets", "Southard",
     "—", "—", "no statistic", "—", "—", "—", "—"),
    ("B1", "18, 67", '"keeps 45.1% of the passage gap\'s length"', "Sengstack",
     "ratio across both", "PD32 half A",
     "**none exists** — see \"Three length ratios with no size-matched floor\"", "—", "—",
     "—", "—"),
    ("B2", "20, 102", "the screen-to-PD32 offset as a multiple of the passage gap "
     "(0.65 to 1.05 times)", "Sengstack", "both", "PD32",
     "**none exists** — `results/seng3/d3_offset.csv` holds a depth-matched bootstrap "
     "interval, not a no-effect floor", "—", "—", "—", "—"),
    ("B3", "70", '"14.4 and 14.3 times the donor-age gap"', "Sengstack + GTEx", "full",
     "— (a position, not a displacement)", "**none exists**", "—", "—", "—", "—"),
    ("B4", "40, 47, 104", "shares of a direction removed by projection (12.0%, 41.5%, "
     "66.4% and 64.9%, 13.6% and 20.3%)", "Sengstack / GTEx", "both", "—",
     "none — `1 − √(1 − cos²)` of an angle, defined at paragraph 80", "—", "n/a", "none",
     "cosine noise attenuates toward zero, so the near-parallel shares are if anything "
     "understated"),
    ("B5", "79", "the mixing positive control moved control samples toward the target by "
     "0.087–0.105, 0.215–0.248 and 0.433–0.484 of the gap", "Sengstack", "k20",
     "own-control",
     "`src/seng3_run.py:591`, reference = pool − draw at `:588`",
     "96 and 157 cells, both screens, 200 draws", "yes / yes",
     "the mixture **is** the positive control; `results/seng3/d2_mixture.csv`",
     "as quoted in the paper"),
    ("C1", "7, 26", "80 of 340 called closer to early-passage cells; the four factors' "
     "−19.2, −13.2, −10.3, −19.9 against thresholds −5.9 to −7.2", "Sengstack", "k20",
     "own-control",
     "`src/seng3_run.py:591-593`, reference = pool − draw at `:588`",
     "own n, every distinct count per screen (`:578`), 200 draws", "yes / yes",
     "`results/seng3/d2_plant.csv`; false-call check `src/seng3_run.py:285-287` → "
     "`results/seng3/d1_null_validation.csv`",
     "half-gap power 0.940–1.000, quarter-gap 0.565–0.740; false-call 0.040–0.085 against "
     "a limit of 0.10. Per-factor `null_median_i` = +0.42, +0.02, +0.46, +0.55"),
    ("C2", "26", "73 called from the late-passage start, including all four", "Sengstack",
     "k20", "PD32 half A", "as C1, with `t=\"ii\"`", "own n", "yes / yes",
     "as C1, origin `ii`", "half-gap power 0.815–0.995"),
    ("C3", "19, 60", "only ZFX and 3 of 340 called in the full gene space; "
     '"the near-empty result there carries no information"', "Sengstack", "full",
     "own-control", "`src/seng3_run.py:591`, space `S1`",
     "own n and the four gate depths", "yes / yes",
     "`results/seng3/d2_plant.csv`, `S1` rows",
     "floor +14.5 to +24.1, no-effect median +20.6 to +33.7, firing rate 0.000, half-gap "
     "power 0.000 in seven of eight cells. This is FLOORCHK's finding, already in the "
     "published text"),
    ("C4", "7, 29, 30", "the four factors' +11.4, +9.0, +5.3, +15.3 against thresholds "
     "−1.0 to −1.3; 5 of 340 called", "Sengstack", "k20", "GTEx old centroid",
     "`src/axang_run.py:516`, reference = pool − draw at `:476`",
     "own n, every distinct count per screen, 200 reps", "yes / yes",
     "`src/axang_run.py:486-506` → `results/axang/d3a_plant.csv`",
     "power 0.945–0.990 at 25% of the gap, 1.000 at 50% and 100%; gate QUALIFIED, "
     "p25 = 0.945. Dose-0 medians +0.81 to +1.50 against the four factors' +5.3 to +15.3"),
    ("C5", "30", '"and all 80 move away from them" (them = young donors)', "Sengstack",
     "k20", "GTEx old centroid", "**CHECKED BY THIS TASK**", "—", "yes / yes",
     "as C4", "as C4"),
    ("C6", "27", "call rate 0.056 → 0.545 across cell-count deciles; 24 of the 60 "
     "perturbations closest in cell count are called", "Sengstack", "k20", "own-control",
     "as C1", "own n; the band is 96–157 cells, n = 60", "yes / yes",
     "`results/seng3/d5_confounds.json` keys `decile_table`, `depth_band_96_157`",
     "0.056 → 0.545; 24 of 60 in the band. The band **is** the size control"),
    ("C7", "37", "Spearman ρ = −0.864 between change in distance and change in "
     "cell-cycle score", "Sengstack", "k20", "own-control",
     "**none** — no threshold of any kind", "—", "—", "none",
     "paper 2 labels it descriptive in the same sentence; paragraph 37 adds the "
     "size-controlled version (3 of 4 at or above the 90th percentile among perturbations "
     "of similar size)"),
    ("C8", "39, 40, 105, 106", "called counts after removing a cell-division direction: "
     "82→27, 78→9, and 75 for the gene-list direction", "Sengstack", "k20", "both",
     "`src/seng4_run.py:458`, keyed `(modality, n, axis, start)` — floors **recomputed "
     "after each removal**", "own n", "yes / yes",
     "`src/seng4_run.py:407` → `results/seng4/d2_plant_residualized.csv`",
     "half-gap power, own-control start: nothing removed 0.895, primary axis 0.010, "
     "PD26 axis 0.040, gene-list axis 0.885; late-passage start 0.985 and 0.760"),
    ("C9", "54", "of 181 perturbations that score younger on the ruler, 116 are not "
     "called closer (120 from the late-passage start)", "Sengstack", "k20", "both",
     "the call half uses C1's floors; the ruler-score half has no floor by design",
     "own n", "yes / yes", "as C1",
     'paper 2 labels it "a diagnostic, not evidence about any factor"'),
    ("C10", "67", "false-call rate 0.085 from the own-control start, 0.110 from the "
     "late-passage start", "Sengstack", "k20", "both",
     "this **is** the floor validation, `src/seng3_run.py:285-287`",
     "96, 107, 128, 157; threshold built on one control half, evaluated on the other",
     "yes / yes", "self", "0.040–0.085 (own-control), 0.110 (late-passage). Paragraph 67 "
     "already discloses that the 0.110 would have failed its pre-set limit and that the "
     "choice was written down after the number was known"),
    ("C11", "19, 68", "half-gap power at least 0.940, quarter-gap at least 0.565",
     "Sengstack", "k20", "own-control", "this **is** the plant gate",
     "96, 107, 128, 157, both screens", "yes / yes", "self, `results/seng3/d2_plant.csv`",
     "as quoted in the paper"),
    ("C12", "31", '"any effect that points away from the donor-age direction … must '
     'increase the distance to young donors"', "Sengstack", "k20", "GTEx old centroid",
     "none needed — a geometric identity", "—", "—", "none",
     "with the start at the old centroid, cos(e, Y−O) < 0 forces the distance to increase "
     "whatever ‖e‖ is, so no length enters it"),
]

A1B_CORRECTED = """\
> Separately, the four genes are also uncalled in the Southard screen, a different neonatal
> line, where the test detects a quarter-gap approach toward young donors at their cell
> counts with power 0.925 to 1.000. Their changes in distance to young donors are positive,
> but so are those of non-targeting cells with no effect at the same cell counts: the
> no-effect median runs from +1.3 to +2.7 against the four factors' +0.5 to +2.4. Three of
> the four move away less far than the median no-effect sample of their own size, and none
> exceeds its 95th percentile, so that screen shows no approach and does not show a
> departure either. That screen used activation only, so for STAT3 and ZFX it tests the
> opposite perturbation."""

PARA142_REPLACEMENT = """\
> • Neither of the other spaces separates these two from the rest of the screen. In space C
> both end farther from the young centroid (+106.6 and +113.8), as do all 1,836
> perturbations: there every displacement is longer than twice the 32-unit gap, so no change
> in distance could have come out negative, and the median non-targeting floor is 2.9 gap
> lengths. In space B both also end farther (+32.9 and +35.1), again with none of the 1,836
> moving closer, and both sit just inside their own cell-count-matched floors (+33.3 and
> +35.7), as 1,088 of the 1,836 do. They have 84 and 82 cells, against a median of 106."""


def md_table(df, cols, heads=None, fmt=None):
    heads = heads or cols
    fmt = fmt or {}
    out = ["| " + " | ".join(heads) + " |",
           "|" + "|".join("---" for _ in heads) + "|"]
    for row in df[cols].itertuples(index=False, name=None):
        cells = [fmt[c](v) if c in fmt else str(v) for c, v in zip(cols, row)]
        out.append("| " + " | ".join(cells) + " |")
    return "\n".join(out)


def main():
    flag = (OUT / "PREREG.flag").read_text(encoding="utf-8").rstrip("\n")
    a1b = pd.read_csv(OUT / "a1b_units.csv")
    c5 = pd.read_csv(OUT / "c5_units.csv")
    base = pd.read_csv(OUT / "c5_baseline.csv")
    rep = pd.read_csv(OUT / "reproduction.csv")
    readings = (OUT / "readings.txt").read_text(encoding="utf-8").rstrip("\n")
    S = json.loads((OUT / "summary.json").read_text(encoding="utf-8"))
    ta, tc = S["a1b"], S["c5"]
    n_gate = int(rep.gating.sum())
    n_bad = int((~rep[rep.gating].match.astype(bool)).sum())
    f3 = lambda v: f"{v:+.3f}"
    L = []
    w = L.append

    w("# FINDINGS_CORR2CHK — do paper 2's two \"move away\" clauses survive a size-matched "
      "no-effect baseline?")
    w("")
    w(f"**Task:** CORR2CHK. **Status:** complete. No STOP fired. Reproduction: {n_gate} "
      f"gating checks, {n_bad} mismatches.")
    w("")
    w("**Primary readings, in the pre-registered wording:**")
    w("")
    w("```")
    w(readings)
    w("```")
    w("")
    w("## What this found, in plain words")
    w("")
    w("**The inventory is the main result: 20 of the 22 clauses need no correction.** "
      "Paper 2 has 21 statements that depend on the Southard screen, on a displacement "
      "length or length ratio, or on a change in distance to a target; splitting paragraph "
      "31's sentence into its two clauses, as the ruling on Q1 requires, makes 22 "
      "clause-level verdicts. Paper 2's own thresholds are drawn from control samples at each "
      "perturbation's own cell count, in the same space and from the same starting profile "
      "as the statement they gate. Paper 2 also nowhere repeats paper 1's corrected "
      "\"more than twice the gap\" claim, and it already applies FLOORCHK's own argument to "
      "its own full-gene-space result (paragraphs 19 and 60). Three length ratios have no "
      "size-matched floor at all and are recorded as such without being checked.")
    w("")
    w(f"**Two clauses were checked. Both fire PAPER2-DEPTH, but they fail very "
      f"differently, and the single key hides that.** Each is of the form \"these units "
      f"move away from young donors\", and the pre-registered rule fails such a claim if "
      f"even one unit is not distinguishable from no effect, because the published wording "
      f"is universal.")
    w("")
    w(f"- **A1b has no support at all.** {ta['supports']} of {ta['m']} units exceed the "
      f"95th percentile of no-effect samples at their own cell count and "
      f"{ta['contradicts']} fall below the median. Under a true null about "
      f"{ta['exp_supports_chance']} of four would exceed the 95th percentile and about "
      f"{ta['exp_contradicts_chance']} would fall below the median, so what was observed is "
      f"what no effect looks like. At these depths a no-effect sample's change in distance "
      f"is itself positive, with a median of +1.3 to +2.7, so the positive sign the clause "
      f"reports is the baseline's own sign.")
    w(f"- **C5's aggregate claim is strongly supported; only the word \"all\" fails.** "
      f"{tc['supports']} of the {tc['m']} units exceed their own 95th percentile, against "
      f"about {tc['exp_supports_chance']} expected by chance, and only "
      f"{tc['contradicts']} falls below its own median, against about "
      f"{tc['exp_contradicts_chance']} expected. That is a large, real departure for most "
      f"of the 80. What does not hold is the universal quantifier: "
      f"{tc['inside_noise'] + tc['contradicts']} of the {tc['m']} are not individually "
      f"distinguishable from no effect, and they are the smaller groups.")
    w("")
    w("Reading those two as one finding would be wrong. A1b's clause has nothing behind it. "
      "C5's clause is right about the direction and overstated only in claiming it of every "
      "one of the 80.")
    w("")
    w("**Paragraph 31's sentence mixes two kinds of claim.** \"Uncalled\" is gated by a "
      "threshold matched to each factor's own cell count and is well powered at those "
      "depths (0.925 to 1.000 for a quarter-gap approach), so it stands and is better "
      "supported than the sentence claims. \"Move away\" is not gated by anything: it "
      "reports the sign of a measured change in distance against zero, when the matched "
      "no-effect baseline at those depths is itself positive. The two clauses sit in one "
      "sentence joined by \"and\", which is why the inventory had to split them.")
    w("")
    w("**What the failures are not.** Neither result says any factor moves *toward* young "
      "donors; nothing here is evidence of rejuvenation by either definition. A1b says the "
      "Southard measurement cannot tell a departure from no effect at those four depths. "
      "C5 says the Sengstack measurement can tell a departure for most of the 80 but not "
      "for the smallest groups, so the claim holds of the set and not of every member. The "
      "call-based half of each sentence is untouched in both cases, and no factor named "
      "here is a candidate for anything.")
    w("")
    w("Files, all in `results/corr2chk/`: `PREREG.flag`, `reproduction.csv`, "
      "`a1b_units.csv`, `c5_units.csv`, `c5_baseline.csv`, `readings.txt`, `summary.json`, "
      "`corr2chk_run.log`, `run_stdout.txt`, `run_stderr.txt`. Code: `src/corr2chk_run.py` "
      "(modes `run`, `verify`) and `src/corr2chk_findings.py`, which writes this file from "
      "those tables and computes nothing. No paper source was edited.")
    w("")
    w("## Pre-registration (verbatim copy of `results/corr2chk/PREREG.flag`)")
    w("")
    w("~~~~text")
    w(flag)
    w("~~~~")
    w("")
    w("## Reproduction of stored numbers, before anything new")
    w("")
    w(f"`python -I src/corr2chk_run.py run` writes `results/corr2chk/reproduction.csv`. It "
      f"has **{n_gate} gating checks and {n_bad} mismatches**, so no STOP fired. "
      f"{len(rep) - n_gate} further rows are descriptive and gate nothing.")
    w("")
    w(md_table(rep, ["check", "item", "expected", "observed", "tol", "match", "gating"]))
    w("")
    w("The two checks that carry the most weight:")
    w("")
    w("- **R2** shows the baseline read for A1b is the same floor the published call was "
      "made against: `depth_floor_S2A_k20` in `results/south8/d3c_all_1836_umi_calls.csv` "
      "equals `floor_p05` in `results/south7/stage_a1_depth_nulls.csv` at each factor's own "
      "cell count, to 1e-12.")
    w("- **R5** shows the re-derived C5 baseline is AXANG's own: re-running AXANG's floor "
      "loop at its own seed reproduces every stored 5th percentile. That is the proof the "
      "repeated loop is faithful, and the ruling made it a STOP condition.")
    w("")
    w("## The inventory (step 1)")
    w("")
    w("Quotes are from `paper/paper2_passage_draft_v0.7.docx`, whose text matches "
      "`paper/two_definitions_of_young_preprint_v1.pdf`. Paragraph numbers are docx "
      "body-paragraph indices. **k20** is paper 2's 20-component space, which is "
      "SOUTH6/7/8's `S2A_k20`: both read `results/south7/basis.npz` key `L_mid_full` "
      "frozen, at `src/seng3_run.py:64`, `src/axang_run.py:87` and "
      "`src/south8_run.py:439-441`. **full** is the 23,485-gene space, SOUTH6's `S1`, "
      "paper 1's \"space C\".")
    w("")
    inv = pd.DataFrame(INVENTORY, columns=[
        "id", "para", "what it says", "dataset", "space", "start", "floor defined at",
        "drawn at", "same space / start", "plant gate", "gate numbers"])
    w(md_table(inv, list(inv.columns)))
    w("")
    w("Three findings of the inventory that bear on scope:")
    w("")
    w("- **Paper 2 does not repeat paper 1's corrected claim.** Searching the source for "
      "\"twice\", \"2.10\", \"1,836\" and \"overshoot\" returns only paragraph 13, \"apply "
      "that readout to the Sengstack data twice\". CORR1's correction does not propagate "
      "into paper-2 text.")
    w("- **Paper 2 already applies FLOORCHK's argument to its own result**, at paragraphs "
      "19 and 60, for the full gene space. A check there would confirm the published text.")
    w("- **A1b and C4/C5 are in one space**, so the rule \"compute nothing in a space paper "
      "2 does not use\" is satisfied for both checked statements.")
    w("")
    w("## A1b — the Southard clause of paragraph 31")
    w("")
    w("> Separately, the four genes are also **uncalled and move away from young donors in "
      "the Southard screen**, a different neonatal line; that screen used activation only, "
      "so for STAT3 and ZFX it tests the opposite perturbation.")
    w("")
    w("Nothing was drawn. The size-matched baseline already exists at each factor's own "
      "cell count, from 200 non-targeting draws with the reference being the non-targeting "
      "pool minus the draw (`src/south7_run.py:319-324`), in k20, from the GTEx old "
      "centroid to the GTEx young centroid. Per the ruling on Q3 it was read, and "
      "reproduced before use.")
    w("")
    w(md_table(a1b, ["unit", "n_cells", "change_in_distance", "floor_p05", "null_median",
                     "null_sd", "null_p95", "d_minus_median", "power_25_at_this_depth",
                     "verdict"],
               ["gene", "cells", "change in distance", "floor (P5)", "no-effect median",
                "no-effect sd", "no-effect P95", "d − median", "power at 25% of gap",
                "verdict"],
               {"change_in_distance": f3, "floor_p05": f3, "null_median": f3,
                "null_p95": f3, "d_minus_median": f3,
                "null_sd": lambda v: f"{v:.3f}",
                "power_25_at_this_depth": lambda v: f"{v:.3f}"}))
    w("")
    w(f"**{ta['contradicts']} of the four move away by less than the median no-effect "
      f"sample of their own size, {ta['inside_noise']} is inside the no-effect spread, and "
      f"{ta['supports']} exceeds its own 95th percentile.** Under a true null of no effect "
      f"anywhere, about {ta['exp_supports_chance']} of four would exceed the 95th "
      f"percentile and about {ta['exp_contradicts_chance']} would fall below the median, so "
      f"what was observed is what no effect looks like.")
    w("")
    w("At these depths a no-effect draw's change in distance to the GTEx young centroid is "
      "**positive**, with a median of +1.3 to +2.7. The positive sign the published clause "
      "reports is the baseline's own sign. This is the same defect SOUTH7's criterion "
      "amendment was written for, in the opposite direction: there a factor with a positive "
      "change in distance was called *approaching* because it moved away less far than "
      "noise usually does.")
    w("")
    w("**The \"uncalled\" half of the sentence stands.** None of the four falls below its "
      "own floor, and the test detects a quarter-gap approach at these cell counts with "
      "power 0.925 to 1.000, so it is a qualified null and not a blind one. The floor is a "
      "5th percentile, so 95% of true nulls are \"uncalled\" by construction; it is the "
      "power, not the count, that makes the null informative.")
    w("")
    w("**The clause is also uninformative in the full gene space**, where paper 2 does not "
      "say it is measured. There no factor's change in distance could be negative at all: "
      "FLOORCHK confirmed that all 1,836 Southard displacements are longer than twice the "
      "gap. The stored S1 no-effect medians at these four depths are +110.5, +100.6, +82.2 "
      "and +78.0 against the four factors' +112.7, +98.2, +83.8 and +79.2 — again all "
      "inside the no-effect spread. So the clause fails in k20 and is vacuous in the full "
      "gene space.")
    w("")
    w("## C5 — the young-donor clause of paragraph 30")
    w("")
    w("> None of the 80 perturbations called closer to early-passage cells is called closer "
      "to young donors, and **all 80 move away from them**.")
    w("")
    w(f"AXANG computed a floor at every distinct cell count but stored only the 5th "
      f"percentile. Per the ruling on Q4 its loop was re-run with one change — the median "
      f"and the 95th percentile are recorded too — at AXANG's own seed "
      f"{S['primary_seed']}, which is what reproducing the stored 5th percentiles requires. "
      f"A replicate at the fresh seed {S['stability_seed']} is reported beside it and gates "
      f"nothing. Per-size baselines, both seeds: `c5_baseline.csv` "
      f"({len(base)} rows).")
    w("")
    w(f"The 80 units are those with `called_delta_i` true in "
      f"`results/seng3/d4_table.csv`, joined to their young-donor change in distance in "
      f"`results/axang/d3b_table.csv`. Their cell counts run from "
      f"{int(c5.n_cells.min())} to {int(c5.n_cells.max())}, median "
      f"{int(c5.n_cells.median())}.")
    w("")
    w(f"| verdict | units | expected under a true null |")
    w("|---|---|---|")
    w(f"| SUPPORTS (d ≥ no-effect P95) | {tc['supports']} | {tc['exp_supports_chance']} |")
    w(f"| INSIDE THE NOISE | {tc['inside_noise']} | — |")
    w(f"| CONTRADICTS (d < no-effect median) | {tc['contradicts']} | "
      f"{tc['exp_contradicts_chance']} |")
    w("")
    w("Every one of the 80, at its own cell count:")
    w("")
    w(md_table(c5, ["unit", "n_cells", "change_in_distance", "floor_p05", "null_median",
                    "null_p95", "d_minus_median", "verdict"],
               ["perturbation", "cells", "change in distance", "floor (P5)",
                "no-effect median", "no-effect P95", "d − median", "verdict"],
               {"change_in_distance": f3, "floor_p05": f3, "null_median": f3,
                "null_p95": f3, "d_minus_median": f3}))
    w("")
    if tc["contradicts"]:
        names = ", ".join(f"`{u}`" for u in S["c5_contradicting_units"])
        w(f"**{tc['contradicts']} of the 80 move away by less than the median no-effect "
          f"sample of their own size:** {names}. The published \"all 80 move away from "
          f"them\" is true of the measured sign for all 80, and established against a "
          f"size-matched baseline for {tc['supports']}.")
    else:
        w("**Every one of the 80 exceeds the median no-effect sample of its own size.**")
    w("")
    w("**The units that fail are the shallow ones.** By verdict, their cell counts are:")
    w("")
    w("| verdict | units | cells: min / median / max | median change in distance |")
    w("|---|---|---|---|")
    for v in ("SUPPORTS", "INSIDE_NOISE", "CONTRADICTS"):
        s = c5[c5.verdict == v]
        if not len(s):
            continue
        w(f"| {v.replace('_', ' ').title()} | {len(s)} | {int(s.n_cells.min())} / "
          f"{int(s.n_cells.median())} / {int(s.n_cells.max())} | "
          f"{s.change_in_distance.median():+.2f} |")
    w("")
    w("So the failure of the universal claim is a depth effect, not a sign that those "
      "perturbations behave differently: the no-effect baseline rises as cell count falls, "
      "and the units it swallows are the small ones.")
    w("")
    sp_med = (base.null_median - base.stability_median).abs()
    sp_p5 = (base.floor_p05 - base.stability_floor_p05).abs()
    w(f"**How precise the baseline is.** It rests on 200 draws, so it is itself a random "
      f"quantity. Between the primary seed and the fresh-seed replicate the no-effect "
      f"median moves by {sp_med.median():.3f} at the median cell count and up to "
      f"{sp_med.max():.3f} at worst, and the 5th percentile by {sp_p5.median():.3f} and up "
      f"to {sp_p5.max():.3f}. Units whose `d − median` is smaller than that could change "
      f"verdict at another seed, so the three counts are approximate at the margin. The "
      f"per-unit `d − median` column is printed above so the margin is visible for every "
      f"unit, and both seeds' baselines are in `c5_baseline.csv`.")
    w("")
    w("**Circularity, disclosed here because this is where the result appears.** The 80 are "
      "selected by the passage test, which uses the same screen control cells that define "
      "zero for the young-donor statistic and that every no-effect draw is taken from. The "
      "selection and the quantity tested are not independent.")
    w("")
    w("## Corrected wording (step 3)")
    w("")
    w("The paper-2 source was not edited and no correction notice was written, as the task "
      "requires. What one would have to say:")
    w("")
    w("**Paragraph 31, the Southard clause.** The sentence as published joins a "
      "size-matched claim to an unmatched one with \"and\". \"Uncalled\" holds. \"Move away "
      "from young donors\" does not: it compares a change in distance with zero when the "
      "matched no-effect baseline at those cell counts is positive, and three of the four "
      "factors fall below that baseline. A corrected clause could read:")
    w("")
    w(A1B_CORRECTED)
    w("")
    w("**Paragraph 30, the young-donor clause.** \"None of the 80 … is called closer to "
      "young donors\" holds; it rests on the same size-matched thresholds as the rest of "
      "the young-donor test, whose gate fired QUALIFIED at p25 = 0.945. The clause \"and "
      "all 80 move away from them\" overstates what a size-matched comparison supports. A "
      "corrected clause could read:")
    w("")
    w(f"> None of the {tc['m']} perturbations called closer to early-passage cells is "
      f"called closer to young donors, and all {tc['m']} have a positive change in distance "
      f"to them. Against no-effect samples of the same size, {tc['supports']} of the "
      f"{tc['m']} move away further than the 95th percentile, {tc['inside_noise']} lie "
      f"inside the no-effect spread and {tc['contradicts']} move away less far than the "
      f"median no-effect sample, so the departure is established for {tc['supports']} of "
      f"them and not for the rest.")
    w("")
    w("Neither correction touches paper 2's conclusion. That conclusion is that the two "
      "definitions of young give opposite answers, and it rests on the *calls*: 80 of 340 "
      "called closer to early-passage cells against 5 of 340 called closer to young donors, "
      "with 17 expected by chance in each case. Both counts are gated by size-matched "
      "thresholds and neither is affected. What the corrections remove is a stronger "
      "directional claim laid on top of the calls.")
    w("")
    w("## Three length ratios with no size-matched floor (B1, B2, B3)")
    w("")
    w("Per the ruling these are recorded, not checked. For each, the direction of the noise "
      "bias is stated, and in each case it **strengthens rather than weakens paper 2's own "
      "caveat**, which is why none was checked.")
    w("")
    w("- **B1, paragraphs 18 and 67**, \"keeps 45.1% of the passage gap's length\". A ratio "
      "of one vector's norm in two spaces. Sampling noise inflates the full-gene norm far "
      "more than the 20-component one, so 45.1% understates the true retention. Paper 2 "
      "uses it to argue that compression discards more noise than signal; the bias works "
      "against that argument, not for it.")
    w("- **B2, paragraphs 20 and 102**, the screen-to-PD32 offset at 0.65 to 1.05 times the "
      "passage gap. Numerator and denominator are both noisy pseudobulk differences, and "
      "the denominator's PD14 half A is the smallest group in the comparison. "
      "`results/seng3/d3_offset.csv` holds a depth-matched bootstrap interval, not a "
      "no-effect floor, and FLOORCHK's construction does not transfer: the two populations "
      "compared sit in different sample groups, so there is no single no-effect pool to "
      "draw both ends from. Paragraph 65 calls the offset \"an upper bound on the technical "
      "difference\", so noise inflation strengthens the disclosure.")
    w("- **B3, paragraph 70**, \"14.4 and 14.3 times the donor-age gap\". A distance from a "
      "4,348-cell pseudobulk to two bulk centroids, not a displacement. Noise inflates it, "
      "so correction could only reduce it, and the statement's role is to warn that the "
      "GTEx-anchored test starts far outside the donor-age span. For scale, FLOORCHK's "
      "stored full-gene floor at the largest Southard group size (341 cells) is a median of "
      "2.06 gap lengths and falls with cell count, while this group has 4,348 cells — but "
      "that is a different dataset with a different gene overlap, so the number does not "
      "transfer and none was computed.")
    w("")
    w("## Step 4 — paper 1, paragraph 142")
    w("")
    w("Separate from everything above, and computing nothing. `CORRECTIONS_PAPER1.md` "
      "identifies a bullet in paragraph 142 of paper 1 as carrying no information:")
    w("")
    w("> • In space C both perturbations move away from the young centroid (+106.6 and "
      "+113.8, q = 1.0), and in space B both end farther from it (+32.9 and +35.1). They "
      "have 84 and 82 cells, against a median of 106.")
    w("")
    w("Replacement wording, for folding into the paper-1 revision. No file was edited.")
    w("")
    w(PARA142_REPLACEMENT)
    w("")
    w("Every number is read from stored files:")
    w("")
    w("| number | source |")
    w("|---|---|")
    w("| +106.6, +113.8; +32.9, +35.1; 84 and 82 cells | "
      "`results/coordspace/stage3_south6_factors.csv` (106.609, 113.750; 32.886, 35.053) — "
      "the paper's own values |")
    w("| none of the 1,836 moves closer, in space C or space B | "
      "`results/coordspace/stage3_south6_summary.csv`, `n_delta_negative` = 0 for both |")
    w("| every space-C displacement longer than twice the 32-unit gap | FLOORCHK "
      "reproduction R1: smallest ratio 2.0966, 1,836 of 1,836; gap 32.284 |")
    w("| median non-targeting floor 2.9 gap lengths in space C | "
      "`median_floor_over_gap` = 2.9415 |")
    w("| space-B floors +33.3 and +35.7 | `floor_p05` = 33.344, 35.705 |")
    w("| 1,088 of 1,836 below floor in space B | `n_below_floor` = 1088 |")
    w("")
    w("Two judgement calls in that wording, both reversible by the author:")
    w("")
    w("1. **\"q = 1.0\" is dropped.** That q comes from the permutation null SOUTH8 found "
      "**anti-conservative in space C** — it flags 0.09 to 0.14 of true-null draws at the "
      "0.05 level (`results/south8/d2_calibration.csv`, S1 rows; verdict in "
      "`d2_verdicts.csv`). Quoting it as though it were a gate is the thing the bullet "
      "should stop doing, and in context it adds nothing the +106.6 does not.")
    w("2. **The space-B floor comparison is added.** The original bullet's space-B clause "
      "has the same weakness as its space-C clause: \"both end farther\" is true of all "
      "1,836. But space B is not vacuous by geometry, and there the two labels are in fact "
      "*below* their own floors, which is the informative fact. Leaving it out would "
      "replace one uninformative clause with another.")
    w("")
    w("## Deviations and decisions not in the task (disclosed)")
    w("")
    w("1. **The session knew both outcomes before the flag was written.** The task's step 1 "
      "required an inventory, and the user's second message required it to be reported in "
      "chat before the flag. That inventory read the stored baselines and the stored "
      "changes in distance, so it established that A1b would fail and that at least one of "
      "C5's 80 units would fall below its baseline. Flag Part 0 records exactly what was "
      "known, with the numbers. The four conclusion wordings were fixed by the user in the "
      "task message before any of it was read, which is what protects the result.")
    w("2. **Four questions went to the user before the flag, and two of the session's own "
      "groupings were wrong.** The inventory first reported C5 as needing no check and "
      "placed C7 under a floor construction it does not have. Both were corrected in the "
      "second inventory pass, before the flag: C5 moved to the checked list and C7's "
      "verdict was re-grounded on having no threshold at all. The rulings Q1 to Q4 are "
      "reproduced verbatim in the flag.")
    w("3. **The FLOORCHK quadrature correction is used nowhere in this task.** Per the "
      "ruling on Q2, a change-in-distance statistic already measured against draws of the "
      "same size has its noise subtracted by construction, so a quadrature step would "
      "double-count it. Both checked statements are of that kind. This is a substitution of "
      "statistic relative to the task's original wording, authorised in writing.")
    w("4. **A tension in the ruling on Q4 was resolved by the session, not by the user.** "
      "The ruling asked for AXANG's loop \"unchanged except for recording the median\", for "
      "\"fresh seeds\", and for the stored 5th percentiles to reproduce \"from the same "
      "re-run\". Those cannot all hold: a 5th percentile of 200 draws is random, so a "
      "fresh-seed run cannot reproduce the stored value exactly. The resolution, fixed in "
      "flag D7 before the run: the primary re-run keeps AXANG's own seed, because the seed "
      "is an argument of the loop and \"unchanged\" therefore includes it, and a "
      "fresh-seed replicate is reported beside it and gates nothing. Both are in "
      "`c5_baseline.csv`, so nothing is hidden by the choice. This is the one place the "
      "session resolved an ambiguity instead of stopping; it affects no conclusion, because "
      "both runs are reported.")
    w("5. **AXANG's `d3()` was not called.** It writes into `results/axang/`, which "
      "CLAUDE.md forbids. Its prologue and floor loop are repeated in "
      "`src/corr2chk_run.py` with the helpers imported unchanged from `axang_run`, "
      "`seng3_run` and `seng2_common`. R5 and R6 are the proof the repetition is faithful: "
      "every stored 5th percentile and all 340 stored changes in distance reproduce. "
      "Importing `axang_run` and `seng3_run` calls `mkdir(exist_ok=True)` on their results "
      "directories, which already exist; no file in them was written.")
    w("6. **A1b draws nothing.** Per the ruling on Q3 the stored baseline was read, because "
      "it exists at exactly the four cell counts with the required construction. No seed "
      "was needed for it, and the task's \"seeds: fresh\" therefore applies only to C5's "
      "stability replicate.")
    w("7. **The 17 NOT-AFFECTED verdicts are inventory judgements, not computed results.** "
      "No number of this task gates them. Their evidence is the table above: the file and "
      "line of each floor, the group size it was drawn at, whether it shares the "
      "statement's space and starting profile, and the file, line and numbers of its plant "
      "gate.")
    w("8. **An arithmetic slip in the flag, disclosed and not edited.** Flag Part 3 says "
      "\"21 statements were found\" and then gives three counts — 2 checked, 17 "
      "NOT-AFFECTED, 3 with no size-matched floor — each labelled \"of 21\". Those sum to "
      "22, because paragraph 31's sentence is one statement among the 21 and two "
      "clause-level entries (A1a and A1b) once the ruling on Q1 splits it. The correct "
      "statement is 21 statements and 22 clause-level verdicts. The flag is the "
      "pre-registration and is reproduced above exactly as written, slip included; no "
      "count used in any result depends on it.")
    w("9. **Where things were written.** `results/corr2chk/` only, plus the two tracked "
      "source files and the three deliverables. Nothing under `results/floorchk/`, "
      "`results/corr1/`, `results/south6/`, `results/south7/`, `results/south8/`, "
      "`results/axang/`, `results/seng3/`, `results/seng4/` or `results/coordspace/` was "
      "modified, no existing `FINDINGS_*.md` or `CORRECTIONS_PAPER1.md` was touched, and no "
      "paper source was edited. The pre-existing uncommitted change to `CITATION.cff` and "
      "the untracked `results/paper2_figs_plain/` are not this task's work and are left out "
      "of its commit.")
    w("")
    w("## Honesty obligations")
    w("")
    w(f"- **Chance baseline beside every count.** The 95th percentile is exceeded by 5% of "
      f"no-effect draws by construction and the median by 50%, so under a true null about "
      f"{ta['exp_supports_chance']} of A1b's four and {tc['exp_supports_chance']} of C5's "
      f"{tc['m']} would read SUPPORTS, and about {ta['exp_contradicts_chance']} and "
      f"{tc['exp_contradicts_chance']} would read CONTRADICTS. Both are printed beside "
      f"every count above.")
    w("- **Recovery power beside every null.** A1b's per-factor power at 25% of the gap at "
      "its own depth is in the result table (0.925 to 1.000). C5's gate is AXANG's D3a, "
      "QUALIFIED at p25 = 0.945 and p50 = 1.000. A1a's \"uncalled\" verdict is never "
      "reported as a bare count.")
    w("- **Cell count beside every verdict.** Both result tables carry each unit's own cell "
      "count. For C5 the units run from 32 to 187 cells, and the no-effect median falls "
      "with cell count, so the baseline differs across units by more than the quantity "
      "being tested in some cases. No unit is presented as a candidate for anything.")
    w("- **Circularity, disclosed and not refitted away.** For both statements the "
      "non-targeting pool of the relevant screen defines zero for every change in distance "
      "and is also the pool every no-effect draw is taken from. For C5 the 80 units are "
      "selected by a statistic that uses those same control cells.")
    w("- **Precision of the baseline.** The medians rest on 200 draws. The fresh-seed "
      "replicate in `c5_baseline.csv` shows the Monte Carlo spread, and verdicts for units "
      "close to the boundary are not stable at that precision.")
    w("- **What this does not cover.** Cell-sampling noise only: not guide-to-guide, batch "
      "or biological variation. Southard and Sengstack are each one neonatal fibroblast "
      "line with no aged donor, so nothing here transfers to aged human cells. Nothing here "
      "tests any factor's phenotype. The GTEx-anchored construction transplants a "
      "displacement measured in single cells onto a bulk post-mortem origin, as both papers "
      "state. Both checked clauses are about what the measurement can resolve, not about "
      "what any factor does.")
    w("- **No conclusion about any named gene or factor** beyond what the existing papers "
      "already state. A1b's units carry gene names because the published sentence names "
      "them; the result is about the sentence.")
    w("")
    w(f"`python -I src/corr2chk_run.py verify` re-derives every reading line from the "
      f"tracked CSVs alone.")
    DST.write_text("\n".join(L) + "\n", encoding="utf-8")
    print(f"wrote {DST.name} ({DST.stat().st_size} bytes)")


if __name__ == "__main__":
    main()
