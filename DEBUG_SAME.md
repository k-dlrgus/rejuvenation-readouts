# DEBUG_SAME

Task: SAME_DEBUG. Pre-registration (`results/same/PREREG.flag`) and `FINDINGS_SAME.md` were not edited. `fire_reading` was not changed. This file is the debug report. `src/same_run.py` contains the one bootstrap fix; the run was repeated; new numbers live in `results/same/`.

## 1. Why the bootstrap CIs did not contain the point estimates

They missed because **O and Y were rebuilt from resampled cells in the bootstrap, and from the frozen split (each half-A cell once) in the point estimate.** That is not a TMM-panel-membership difference.

### TMM panel: same rows, re-run on different counts

Not a different pseudobulk panel.

Shared transform in both paths is `panel_z`:

```248:253:src/same_run.py
def panel_z(count_rows, frozen):
    C = np.vstack([np.asarray(r, np.float64).ravel() for r in count_rows])
    logcpm, nf = tmm_logcpm_quiet(C)
    Z, missing = zscore_frozen(logcpm, frozen, C)
    return Z, missing, nf, C
```

`tmm_logcpm_quiet` is TMM then edgeR log2-CPM `prior.count=2` on the rows it is given. It does not change which states are in the panel.

**Point estimate (main), unchanged:**

```964:976:src/same_run.py
    counts = [sum_rows(Y, ix) for _, Y, ix, _, _ in panel_spec]
    names = [p[0] for p in panel_spec]
    Z, missing, nf, C = panel_z(counts, frozen)
    ...
    z_O, z_Y = Z[pos["O"]], Z[pos["Y"]]
    z_S, z_Srev = Z[pos["S"]], Z[pos["S_rev"]]
    ...
    fwd = axis_stats(z_O, z_Y, z_S)
```

`panel_spec` is always six rows: `O, Y, S, S_rev, Pluri, NonReprog`. TMM reference on the original panel is **S** (upper-quartile / library closest to the mean). Library sizes: O 4.07e7, Y 3.70e7, S 7.33e7, S_rev 8.02e7, Pluri 2.56e6, NonReprog 5.95e7. nf = `[0.777, 0.878, 1.041, 1.076, 1.488, 0.880]`.

**Bootstrap (buggy, as run for FINDINGS_SAME.md):**

```python
b_counts = [sum_rows(Y, resample_idx(ix, rng_b)) for _, Y, ix, _, _ in panel_spec]
Zb, _, _, _ = panel_z(b_counts, frozen)
```

Same six names. TMM is re-run on resampled counts of those six, not on a different set of states. In 40/40 diagnostic draws the TMM reference stayed S. Applying the *original* TMM nf to bootstrap counts (no re-TMM) left the bias unchanged (`frac_S` bias −0.030 vs −0.030; `frac_rev` bias +0.102 vs +0.102). Re-TMM is not the shift.

Sparse CSR fancy-index with replacement is not the shift either (scipy 1.18.0 duplicates rows; unique-cell fraction in an O resample is 0.633 ≈ 1−1/e; mean bootstrap library size / original ≈ 1).

### O/Y: rebuilt in bootstrap, not in the point estimate

Yes. That is the bug.

| | O, Y | S, S_rev, Pluri, NonReprog | TMM |
|---|---|---|---|
| Point estimate | `sum_rows(Y, frozen_halfA)` — each cell once | `sum_rows(Y, frozen_state_idx)` | `panel_z` on those 6 rows |
| Bootstrap (old) | `sum_rows(Y, resample_idx(halfA, replace=True))` | `sum_rows(Y, resample_idx(state, replace=True))` | `panel_z` on those 6 rows |

`frac = dot(d, v) / ||v||²` with `v = z(Y)−z(O)` and `d = z(S)−z(O)`. Resampling O injects a high-dimensional origin error `ε_O = z(O) − z(O*)` into **both** `d*` and `v*`:

```
d* = d + ε_O + ε_S
v* = v + ε_O + ε_Y
E[dot(d*, v*)] ≈ dot(d, v) + E[||ε_O||²]
```

The extra `E[||ε_O||²]` term (p = 23,485) shifts the bootstrap. It is large enough to **sign-flip** `frac_rev`. Empirical confirmation on the original `boot.npz` (B=200): **0/200** draws of `frac_S` were ≥ the point; **200/200** draws of `frac_rev` were ≥ the point. That is a location shift, not a wide percentile interval.

Diagnostic B=40, same seed, holding O/Y **counts** at the point-estimate vectors and still re-TMM'ing after resampling S / S_rev / Pluri / NonReprog:

| | boot mean | point | point in CI? | bias |
|---|---|---|---|---|
| re-TMM, O/Y resampled | `frac_S` +0.644 | +0.674 | no | −0.031 |
| re-TMM, O/Y resampled | `frac_rev` +0.045 | −0.057 | no | +0.102 |
| frozen nf, O/Y resampled | `frac_S` +0.645 | +0.674 | no | −0.030 |
| **re-TMM, O/Y held** | `frac_S` **+0.674** | +0.674 | **yes** | **−0.0003** |
| **re-TMM, O/Y held** | `frac_rev` **−0.058** | −0.057 | **yes** | **−0.0006** |

The pre-registration text said to resample O, Y, S, S_rev. Implementing that literally rebuilds the frozen-split anchors. The point estimate never does. The CIs were therefore not intervals around the reported statistic.

### Exact code path after the fix (main)

Point estimate: still `counts = [sum_rows(Y, ix) for ... panel_spec]` then `panel_z(counts)`.

Bootstrap: same `counts[i]` for O and Y; destinations still resampled; still `panel_z` (re-TMM on the same six names):

```1011:1028:src/same_run.py
    rng_b = np.random.default_rng(BOOT_SEED)
    ...
    for b in range(N_BOOT):
        b_counts = []
        for i, (_, Y, ix, _, _) in enumerate(panel_spec):
            if names[i] in ("O", "Y"):
                b_counts.append(counts[i])
            else:
                b_counts.append(sum_rows(Y, resample_idx(ix, rng_b)))
        Zb, _, _, _ = panel_z(b_counts, frozen)
```

Nothing else was retuned (TMM, prior.count, gene filter, split, f values, `fire_reading`).

## 2. Why positive-control f=0 gives frac=+0.115, not ~0

Not library size. Not (in any material way) extra mixture rows in the TMM panel.

`f=0` is all 2470 aged half-B d0 Fibroblast cells, scored as S against half-A anchors. Random 50/50 split, seed 20260914.

### Mixture rows do change the TMM panel — but that is not why frac(f=0) is 0.115

P point estimate TMM panel is six rows: `O, Y, mix_f0.00, mix_f0.10, mix_f0.25, mix_f0.50` (`run_posctrl` → `panel_z`). The main-task panel is a different six (`O, Y, S, S_rev, Pluri, NonReprog`). P is run before any PartialReprog vector, so it cannot TMM on the main panel.

TMM reference on the full P panel is mix_f0.25. nf = `[0.965, 1.091, 0.969, 0.976, 0.991, 1.014]`.

Scoring mix_f0 on a **3-row** panel `[O, Y, mix0]` (other mixtures absent) gives frac=**+0.117**, delta=+0.053. Full 6-row P panel: frac=**+0.115**, delta=+0.068. Dropping the other mixes does not send the noise floor to 0.

P bootstrap (old) also used that same 6-row membership, with O, Y, and every mix resampled.

### Half-B vs half-A aged cells do not differ in library size

On the frozen-ruler gene space (sum of `louvain_rulerY` counts):

| | n | total UMI | mean UMI/cell | median |
|---|---|---|---|---|
| aged half A (O) | 2469 | 4.071e7 | 16490 | 14411 |
| aged half B (mix f=0) | 2470 | 4.086e7 | 16541 | 14371 |
| young half A (Y) | 3851 | 3.702e7 | 9612 | 8118 |
| young half B | 3852 | 3.716e7 | 9646 | 8121 |

Aged half-B / half-A library ratio = **1.0035**. Mean UMI/cell difference z = **−0.23**. Young halves similarly (ratio 1.0037, z = −0.28). Forcing TMM nf of mix0 equal to nf of O leaves frac(f=0) at **+0.114**.

Raw count correlation `corr(c_O, c_mix0) = 0.999949`. Relative L1 difference = 1.76%.

### What actually produces +0.115

The two halves are almost the same counts. The **transform** stretches that residual onto `v = z(Y)−z(O)`:

| transform | frac(mix0 vs O toward Y) | \|\|mix0−O\|\| |
|---|---|---|
| log2(CPM+2), no TMM | **+0.018** | 11.05 (log-CPM) |
| TMM + edgeR log2-CPM prior.count=2, then frozen GTEx z | **+0.115** | 49.4 (TMM log-CPM); 59.6 (z) |

`~0` is what you get before TMM+z. After the pre-registered transform, a 1.8% L1 residual between independent halves of the same d0 Fibroblast pool becomes `||d||=60` in z-space and frac=+0.115, because edgeR prior.count=2 puts zeros and singleton detections far apart in log-CPM, and TMM is estimated in a panel that includes Y (compositionally different from O/mix0). That is the noise floor of this statistic. It was not tuned.

The old P bootstrap then resampled O and Y as well, so the f=0 **CI** `[+0.168, +0.195]` sat above the point +0.115 for the same shared-origin reason as §1.

## 3. Fix, rerun, old vs new

**Fix (only this):** stop rebuilding O and Y inside the bootstrap. Use the same half-A count vectors as the point estimate. Still resample destinations (S, S_rev, Pluri, NonReprog) and P mixture cells. Still re-TMM the same panel membership. Seeds, n_boot, split, gene sets, `fire_reading` unchanged.

Rerun: `python src/same_run.py` (exit 0). Point estimates identical to the original run (same split, same sums, same TMM on the original panel). Only CIs moved.

### Primary statistics

| stat | point (old = new) | old 95% CI | new 95% CI | old covers point? | new covers point? |
|---|---|---|---|---|---|
| frac_S | +0.674 | [+0.624, +0.666] | [+0.667, +0.680] | no | **OK** |
| frac_rev | −0.057 | [+0.026, +0.063] | [−0.063, −0.054] | no | **OK** |
| Asymmetry | +0.732 | [+0.562, +0.640] | [+0.724, +0.740] | no | **OK** |
| delta_S | +130.894 | [+119.669, +125.224] | [+130.473, +133.918] | no | **OK** |
| delta_rev | +163.342 | [+152.488, +156.719] | [+163.312, +166.330] | no | **OK** |
| cos_S | +0.323 | (boot mean was +0.336; 196/200 draws > point) | not a gate | no | not re-gated |
| p_N1 forward | 0.0050 | — | — | — | unchanged |
| p_N1 reverse | 1.0000 | — | — | — | unchanged |
| n_cells_S | 5832 | — | — | — | unchanged |

`frac_rev` CI is now entirely negative, matching the negative point. That is the shared-origin sign flip going away, not a change in the reverse test.

### Positive control P

P still passes: monotonic frac, f=0.50 frac CI excludes 0, delta(f=0.50) point < 0.

| f | frac point | old frac CI | new frac CI | delta point | old delta CI | new delta CI |
|---|---|---|---|---|---|---|
| 0.00 | +0.115 | [+0.168, +0.195] **INVALID** | [+0.099, +0.121] **OK** | +0.068 | [−2.060, +1.972] OK | [+5.807, +8.334] **INVALID** |
| 0.10 | +0.182 | [+0.225, +0.252] **INVALID** | [+0.167, +0.191] **OK** | −7.870 | [−9.312, −5.248] OK | [−1.900, +1.007] **INVALID** |
| 0.25 | +0.283 | [+0.309, +0.335] **INVALID** | [+0.267, +0.291] **OK** | −18.913 | [−18.807, −14.793] INVALID | [−11.793, −9.372] **INVALID** |
| 0.50 | +0.441 | [+0.438, +0.463] OK | [+0.424, +0.447] **OK** | −34.443 | [−31.045, −26.909] INVALID | [−25.403, −22.414] **INVALID** |

Holding O/Y fixes P **frac** CIs (the P gate that uses a CI). P **delta** CIs remain off the point; see §4. Frozen original TMM nf on P bootstrap draws does not fix delta (B=40: f=0 delta bias +6.85 re-TMM vs +6.89 frozen nf). That leftover is high-p Euclidean inflation `E[||z(S*)−z(Y)||] > ||z(S)−z(Y)||` when S is resampled, not a second implementation mismatch. It was not tuned.

### Context (not the reading)

| | point | old CI | new CI | new sanity |
|---|---|---|---|---|
| Pluri frac | +1.403 | [+1.165, +1.294] | [+1.365, +1.442] | **OK** |
| Pluri delta | +382.106 | [+379.979, +406.190] | [+391.791, +418.643] | **INVALID** |
| NonReprog frac | +0.555 | [+0.532, +0.561] | [+0.543, +0.561] | **OK** |
| NonReprog delta | +60.331 | [+52.577, +56.671] | [+62.804, +64.575] | **INVALID** |

Pluri n=439. Distance inflation is largest there. Not tuned.

### Pre-registered reading, re-applied unchanged

Same function `fire_reading`, same thresholds, same order.

Conditions on the **fixed** numbers:

- P pass = True (monotonic, frac(f=0.50) CI [+0.424, +0.447] excludes 0, delta(f=0.50)=−34.443 < 0)
- frac_S CI excludes 0 and is positive: True `[+0.667, +0.680]`
- delta_S < 0: **False** (delta_S = +130.894)
- N1 p ≤ 0.05: True (p = 0.0050)
- Asymmetry CI excludes 0 and is positive: True `[+0.724, +0.740]`
- frac_S CI includes 0: False
- delta_S ≥ 0: True
- Pluripotency overshoot (delta < 0 and frac CI > 0): False (delta_Pluri = +382.106)

Rule 2 needs delta_S < 0; it does not fire. Rule 4: `delta_S >= 0 OR frac_S CI includes 0` → **fires**.

**Fired key is unchanged: `no_approach_same_platform`.** Flags still `[]`. Text still: the TOWARD negative holds with the platform gap removed; GM00731 PartialReprog n_cells=5832, frac_S=+0.674, delta_S=+130.894 (not closer to the young donor).

## 4. CI sanity (every reported interval)

Rule: the point estimate must lie inside the 95% percentile CI, else the row is **INVALID**. No interval was widened or recentered to pass this.

| row | stat | point | CI | |
|---|---|---|---|---|
| P f=0.00 | frac | +0.115 | [+0.099, +0.121] | OK |
| P f=0.00 | delta | +0.068 | [+5.807, +8.334] | **INVALID** |
| P f=0.10 | frac | +0.182 | [+0.167, +0.191] | OK |
| P f=0.10 | delta | −7.870 | [−1.900, +1.007] | **INVALID** |
| P f=0.25 | frac | +0.283 | [+0.267, +0.291] | OK |
| P f=0.25 | delta | −18.913 | [−11.793, −9.372] | **INVALID** |
| P f=0.50 | frac | +0.441 | [+0.424, +0.447] | OK |
| P f=0.50 | delta | −34.443 | [−25.403, −22.414] | **INVALID** |
| forward | frac_S | +0.674 | [+0.667, +0.680] | OK |
| forward | delta_S | +130.894 | [+130.473, +133.918] | OK |
| reverse | frac_rev | −0.057 | [−0.063, −0.054] | OK |
| reverse | delta_rev | +163.342 | [+163.312, +166.330] | OK |
| Asymmetry | frac_S − frac_rev | +0.732 | [+0.724, +0.740] | OK |
| Pluri | frac | +1.403 | [+1.365, +1.442] | OK |
| Pluri | delta | +382.106 | [+391.791, +418.643] | **INVALID** |
| NonReprog | frac | +0.555 | [+0.543, +0.561] | OK |
| NonReprog | delta | +60.331 | [+62.804, +64.575] | **INVALID** |

Gates that use a CI (P frac at f=0.50; frac_S; Asymmetry) are OK. INVALID leftover rows are all **delta** (Euclidean distances of a resampled high-p vector). The reading does not use those CIs as gates (`delta_S < 0` is the point; P's `delta(f=0.50) < 0` is the point).

## What was not done

- PREREG.flag, `FINDINGS_SAME.md` body, and `fire_reading` were not edited.
- `src/toward_run.py` was not modified.
- TMM, prior.count, panel membership, split, f grid, and n_boot were not retuned to hide INVALID deltas.
- GSE325735 was not opened.
