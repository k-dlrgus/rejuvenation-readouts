# Paper 2 figures (working copy for paper2_passage_draft_v0.6)

This folder is updated in place each time paper/paper2_passage_draft_v0.6 is
edited. Plotting only: every value is read from stored results; nothing is
recomputed. Earlier figure folders are unchanged.

| figure | content |
|---|---|
| Fig 1 | study design |
| Fig 2 | both targets, threshold-relative axes |
| Fig 3 | cell division (3A: every direction in Table 1; 3B: removal of each screen's own control-cell direction) |
| Table 1 | every cell-division direction (rows) against each definition of aging (columns); built by table1.py |
| S1 | test validation, full gene vs 20-component space |
| S2 | structure of the data |
| S3 | everything from the late-passage start |
| S4 | directional age score vs distance |

Table 1 values: results/axang/d1_angles.csv and d1_angles_corrected_C1lu.csv.
PD26 and gene-list removal results (FINDINGS_SENG4.md, D1-D3) are given in the text only.

    python fig12.py results/paper2_figs_fix/run3 <out>                 # Fig1, FigS1
    python fig23.py results/paper2_figs_fix/run3 results/paper2_figs_fix/run6 \
        results/seng3/d4_table.csv results/axang/d3b_table.csv <out>   # Fig2, Fig3
    python figS.py  results/paper2_figs_fix/run3 results/paper2_figs_fix/run6 <out>  # S2-S4
    python results/paper2_figs_v06/table1.py                          # Table 1 -> table1.csv

