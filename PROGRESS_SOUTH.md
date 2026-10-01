# PROGRESS_SOUTH

**STOP status:** STOP at Stage 0; keys fired: `no_raw_counts`.

**Next action:** Halt as pre-registered on `no_raw_counts`. This file holds per-guide 'mean population' regression output (X is float64, 57% of its entries negative, min -9.16) with no raw layer and no counts layer, so the frozen z-pipeline cannot start and neither the floor nor any delta can be formed. The 17 `off-target` control rows do exist, so `no_controls` does not fire and a floor would be constructible from a counts table. Reported rather than acted on: `L4A_matrix.h5`, in the same folder, does hold raw integer gene counts (36,601 genes x 18,197 cells), but it is not the table this pre-registration pins under DATA, and it is one lane holding 2.4% of the screen's cells at roughly 1.7 cells per guide. Swapping it in now would be choosing the data after seeing the pinned table fail. Stages 1-4 would need a new pre-registration naming a count-level table and the remaining lanes; that is a new-data decision outside this task.

