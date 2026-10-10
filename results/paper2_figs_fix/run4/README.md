PLAN written 2026-10-10T07:36:53Z
Task PAPER2FIGS-fix2. Cosmetic only. Written before Fig 3 is redrawn. Not edited after.

TASK: redraw Fig 3 only into results/paper2_figs_fix/run4/ with src/paper2_figs.py, FIX on,
LABEL_FALLBACK on (as run 3). Move the CRA_E2F3 label in panel C so its leader goes directly to
CRA_E2F3's point (0.1654, -19.21) and not through or beside CRI_ZFX's point (0.1288, -19.94).

PLACEMENT (drawing choice, not a result; not a STOP):
  CRA_E2F3: text anchor at data (0.205, -21.2), ha=left, va=center, straight leader to the point,
  below and to the right of the point, away from CRI_ZFX (which lies to the left).
  The other three labels (CRA_EZH2, CRI_STAT3, CRI_ZFX) keep run 3's deterministic search;
  CRA_E2F3 is excluded from it and its box is added to the occupied list afterwards.
  If this position overlaps a point or label it is reported with the image, not adjusted.

HARD RULE (only STOP besides CLAUDE.md's): run4/Fig3_data.csv, after CRLF->LF, equals
results/paper2_figs_fix/run3/Fig3_data.csv byte for byte (also no lone CR). Fail -> STOP data_changed.

Code change beyond the flag: place_labels gets an optional 'manual' argument (name -> data xy);
new src/paper2_figs_fix2_run.py draws only Fig 3 and compares. Run 3 files are not written to.

--- REVISION 2026-10-10T07:37:30Z (appended; text above is superseded for placement, not edited) ---
Attempt 1 (CRA_E2F3 text at (0.205, -21.2)) was drawn and reported: its text overlapped 2 points and
crossed the axes bottom, and putting it last in the order let the search move the other three labels
(CRA_EZH2's label then sat beside CRI_ZFX's point). Attempt 1 is superseded; its output is overwritten.
Attempt 2, drawn after this line: CRA_E2F3 text anchor at data (0.22, -3.6), ha=left, va=center,
in empty space above-right of the cluster; leader may cross points but not pass beside CRI_ZFX.
Code change: manual labels keep their original position in the search order (so the other three
are searched as in run 3). Attempt 2 is final; any remaining overlap is reported, not adjusted.

=== RESULT (written after drawing) ===
Command: python src/paper2_figs_fix2_run.py (anaconda python; only Fig 3 drawn).
Hard rule: run4/Fig3_data.csv equals run3/Fig3_data.csv after CRLF->LF, no lone CR: PASS.
Image: results/paper2_figs_fix/run4/Fig3.png (and .pdf).
Attempt 2 (final) placement: CRA_E2F3 text at (0.22, -3.6) data. Its leader runs from the text down-left
to CRA_E2F3's own point (0.165, -19.2); CRI_ZFX's point (0.129, -19.9) lies to the left of the leader
and is not touched by it. The leader crosses some orange points (allowed by the flag).
Label check (Fig3C, reported not gated): all four labels overlap 0 points, 0 labels, 0 legend/rho text.
Reported failure: CRA_E2F3's text box extends past the right edge of the axes (inside_axes = False;
it ends in the right margin, nothing else is there). Not adjusted, per the flag.
Observation, not asked: CRA_EZH2's label sits close to CRI_ZFX's point (no overlap by the check).
Attempt 1 (superseded in the flag): text overlapped 2 points and left the axes; its output was overwritten,
and the other three labels moved because of search order. Disclosed deviations: place_labels got a 'manual'
argument in src/paper2_figs.py (FIG3_MANUAL = None by default, so earlier runs are unchanged); new
src/paper2_figs_fix2_run.py; Fig 3 only was redrawn, so run4 has no other figures, no reproduction stage
and no paper2_numbers rows. Cosmetic only; no number changed.
