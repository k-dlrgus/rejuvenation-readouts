PLAN written 2026-10-10T07:40:29Z
Task PAPER2FIGS-fix3. Cosmetic only. Written before Fig 3 is redrawn. Not edited after.

TASK: redraw Fig 3 into results/paper2_figs_fix/run5/ from run4's settings (FIX on, LABEL_FALLBACK on,
CRA_E2F3 manual at (0.22, -3.6)); one change in panel C: swap label positions of CRI_ZFX and CRA_EZH2.

PLACEMENT (drawing choice, not a result; not a STOP). Text anchors in data coords, ha=left, va=center,
straight leaders, read off run4/Fig3.png:
  CRA_EZH2: text at (0.16, 1.4)   [where CRI_ZFX's label is in run4], leader to point (0.148, -13.2).
  CRI_ZFX:  text at (0.008, -20.4) [where CRA_EZH2's label is in run4], leader to point (0.129, -19.9).
  CRA_E2F3 (0.22, -3.6) and CRI_STAT3 (search, as run 4) unchanged.
  Leaders must not pass through another factor's circled point; checked visually and reported.
  Any residual overlap is reported with the image, not adjusted.

HARD RULE: run5/Fig3_data.csv equals run4/Fig3_data.csv after CRLF->LF (no lone CR).
Code beyond the flag: new src/paper2_figs_fix3_run.py (copy of fix2 runner with new manual dict).

=== RESULT (written after drawing) ===
Command: py src/paper2_figs_fix3_run.py (Fig 3 only). Hard rule PASS: Fig3_data.csv LF-equal to run4, no lone CR.
Placement as planned. Label check (reported, not gated): 4 labels, 0 points / 0 labels / 0 legend-rho overlaps, all inside axes.
Visual: CRA_EZH2 leader runs up-right from its point (0.148,-13.2) to the upper label; CRI_ZFX leader is a short
horizontal line to its point (0.129,-19.9). No leader passes through another factor's circled point.
CRI_STAT3 stayed at search offset and CRA_E2F3 at (0.22,-3.6): unchanged from run4. The CRA_E2F3 leader and the
CRA_EZH2 leader do not cross. Cosmetic only; no number changed.
