# PROGRESS_SOUTH3

**STOP status:** no STOP key fired; the task ran to completion.

**Report-only flags:** `scale_unknown`, `coverage_low`.

**Next action:** Nothing further is required by this pre-registration: Stages 0-4 and the three SOUTH3 patches all ran and are reported in FINDINGS_SOUTH3.md. SOUTH3's own three changes are spent — the null now resolves q, the guide rule is single-valued, and the displacement geometry is reported. What would change the answer rather than decorate it is unchanged from SOUTH2: (a) an input whose units are stated, since `scale_unknown` fired and the MDA is not scale-invariant, and (b) a perturbation dataset measured in aged cells, since every delta here is a counterfactual built from effects measured in one neonatal line. Both need a new pre-registration.

**Result in one line:** `toward_young` is 0 in 2 of the 28 settings and non-zero in 26, ranging 40-277 of 1836 factors. Non-zero settings are V_all, V_sig. The counts are not pooled and no single setting is primary. **No.** Across the 26 settings that have survivors, the Spearman correlation between the delta ranking and the cosine ranking over all 1,836 factors is at least +0.871, and within every quartile of ||d|| / ||target - origin|| in those settings it is at least +0.878. Taking the top-n factors by cosine alone, n being that setting's number of survivors, recovers 2,329 of the 3,326 survivors (70.0%) — and taking the top-n by delta itself, the statistic the label is built on, recovers 2,334 (70.2%). The two agree to within a percentage point, so the ~30% that neither recovers is the conjunction of the floor and q gates, not information the cosine lacks. ||d|| alone recovers only 542 (16.3%), so the label is not merely picking the biggest displacements. The mechanism is the identity: in every setting that has a survivor the median residual of delta = -||d||cos(theta) is at most 5.6% of ||d||, so delta there is the cosine re-weighted by the displacement's own length and carries nothing beyond it. No survivor in any setting is distinguishable on direction from what a plain cosine ranking would give. What the delta machinery does buy is not a different ordering but a bar: it says which of those cosines are larger than the factor's own reshuffled null, and a cosine ranking on its own never says that. `coverage_low` fired, so every candidate is labelled exploratory.

**Attainable minimum permutation p:** 4.999750e-05 = 1/(20,000+1).

**Supersession:** SOUTH3 supersedes SOUTH2. One line was added at the top of `FINDINGS_SOUTH2.md` marking it superseded; nothing else in that file was changed.

