Task PAPER2FIGS-fix4 (plan: PLAN.flag, written before drawing). Cosmetic only.
Command: py -X utf8 src/paper2_figs_fix4_run.py  (Fig 3 from run5 settings; Fig S2 from run3 settings)

Strings changed (old -> new), src/paper2_figs.py:
 1. Fig 3 panel C x-label: "Cell-division shift (perturbation − NT)" -> "Cell-cycle score change (perturbation − NT)"
 2. Fig S2 panel D x-label: "Ruler genes (grey present, white absent)" -> "Ruler genes (gray present, white absent)"

US-spelling search: every string literal in src/paper2_figs.py was scanned (grey, colour, centre, labelled,
-ise/-isation, -our, -yse, ageing, etc.). Only "grey" in string 2 was found; "grey" in code (ff.GREY) is an identifier, not drawn text.
Only Fig 3 and Fig S2 text changed, so only they are redrawn.

Hard rule: Fig3_data.csv LF-equal to run5; FigS2_data.csv LF-equal to run3 (both PASS, no lone CR).
Code beyond the flag: src/paper2_figs_fix4_run.py. Fig 3 visual check: labels/leaders unchanged from run5.
