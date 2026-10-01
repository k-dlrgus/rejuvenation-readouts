# PROGRESS_POSCTRL

**STOP status:** P1 fired under a miscalibrated angle bar; superseded 2026-09-16; C1/C3 below. S2 2026-09-17: direction transfers, scale does not. S3 2026-09-17: the age axis transfers after removing the identity subspace.
**Next action:** done

## Pre-registered block (verbatim)

- **P1 pipeline broken:** sex control (C2a, with sex-chromosome genes) fails bootstrap ≤35° or transfer >0 → STOP, report, do not proceed to C1/C3. Everything downstream of the target pipeline is suspect.
- **P2 negative is informative:** planted dense-shared direction at R²≈0.30 (C1) reaches bootstrap ≤35° AND transfer >0 in both directions → the pipeline can find a real direction at this n; the age negative in FINDINGS_LOWDIM stands as a statement about the data.
- **P3 negative is a power statement:** planted dense-shared at R²≈0.30 gives bootstrap >35° OR transfer ≤0 in either direction → the LOWDIM bars were unattainable at n=233; "not identifiable" must be rewritten as "not identifiable at n=233"; the power curve (C3) becomes the primary result.
- **Mixed** (passes angle, fails transfer, or vice versa) → report which, do not average. Transfer failing on a planted shared direction means bank covariance differences alone defeat ridge transfer; say so.

## Supersession 2026-09-16 (verbatim, written before any C1/C3 fit)

The ≤35° angle bar is retired. New pre-registered bars for C1 and C3:
- (a) transfer R² > 0 in both directions with permutation null ≤ 0.05 — primary
- (b) bootstrap median pairwise angle reported against the sex reference 41.3°, not against 35°. Angle is not a pass/fail gate.

Restated P2 / P3 (these bars; C1 dense_shared at R²≈0.30):
- **P2** = planted dense_shared at R²≈0.30 transfers both ways
- **P3** = it does not

Sex is sparse and huge, so it is not a fair reference for a planted dense continuous direction at R²≈0.30. C1 is.

Verdict rule (only the outcome that fired):
- if planted transfers and real age does not → the age negative is real and it is specifically a transfer failure, not an angle failure
- If planted does not transfer either → bank covariance differences defeat transfer for any direction at this n; the age negative is uninformative and the power curve is the result
- If mixed across regimes, report per regime, do not average.

No threshold in this block is re-tuned after numbers exist.

## Supersession 2 — 2026-09-17 (verbatim, written before any S2 fit)

Transfer is now also reported as held-out-bank Pearson r and Spearman ρ between predicted and true target, plus R² after refitting intercept and slope on the test bank (calibrated R²), each with the same donor-level permutation null. Pass = r > 0 with null ≤ 0.05 in both directions.

Rerun only the transfer cells (no bootstrap, no C3 refit) for: gene_ridge real age, C2a sex, all four C1 regimes at 0.30, dense_shared at 0.15 and 0.50. Also reload the C3 transfer cells from results/posctrl/*.json if predictions were saved; if not, rerun C3 transfer only (10 draws per n, B/2) for planted 0.30 and real age.

Verdict rule (only the outcome that fired):
- (i) planted r>0 and real age r>0 both ways → direction transfers, scale does not; the LOWDIM negative is a calibration failure
- (ii) planted r>0, real age r≤0 → the age negative is real
- (iii) planted r≤0 → uninformative, power curve stands

No threshold in this block is re-tuned after numbers exist.

## Supersession 3 — 2026-09-17 (verbatim, written before any S3 fit)

S2 used n_perm=20. Three transfer-only tasks; no bootstrap-angle; no C3 refit. Same L2c / V4a path; age direction, not TARGET. Aggregation: median over cell types. Permutation seed `20260914`. Donor-bootstrap seed `20260918`. Draw order: HBCC→MSSM test bank first, then MSSM→HBCC.

(1) Harden: rerun transfer r, ρ, calibrated R² with n_perm=200 for gene_ridge real age, C2a sex, dense_shared @0.30, both directions. Donor-bootstrap 95% CI on r (B=200, resample the test bank). Pass remains S2: r > 0 with null ≤ 0.05 both ways.

(2) Bank-specific share: real-age r / planted-@0.30 r per direction, with the same paired donor-bootstrap 95% CI. "fraction of transferable age signal". If planted r ≤ 0, share is undefined.

(3) Identity control: identity-residualized age fit from FINDINGS_GEOMETRY Part C (the 0.317→0.300 fit) — train-only type-centroid identity basis, residualize expression, refit ridge age per type. Transfer r, ρ, calibrated R² both directions, n_perm=200, null. Pass: r > 0 with null ≤ 0.05 both ways.

Verdict rule (only the outcome that fired):
- if (3) passes → the age axis transfers after removing the identity subspace
- if it fails → transferable age signal lives in the identity subspace

No threshold in this block is re-tuned after numbers exist.

## Seeds

- main / battery: `20260914`
- planted w: `20260915`
- β calibration: `20260916`
- C3 subsample: `20260917`
- S3 donor-bootstrap: `20260918`
- dataset: `4442d412-91cb-4261-acca-8adf5fa04c11`

## Frozen β per regime

- `dense_shared_r030`: β=0.24683784844272705
- `dense_typespec_r030`: β=0.24479061664766935
- `sparse_shared_r030`: β=0.3226823122214063
- `sparse_typespec_r030`: β=0.2667642179173281
- `dense_shared_r015`: β=0.18666913590500295
- `dense_shared_r050`: β=0.34711538643104906

## Finished cells

### c2a
- boot_median_angle: 41.336132228431666
- achieved_within_site_r2: 0.3595875112208249
- within_site_null: -0.1368177906451416
- loso_r2: 0.4097116614066296
- transfer_H_to_M: 0.5117665952957187
- transfer_M_to_H: 0.35822098443172434
- transfer_H_to_M_null: -0.1193174174500941
- transfer_M_to_H_null: -0.13971036646664056
- young_old_balanced: 75.64187843667933
- pass_angle: False
- pass_transfer: True
- pass_both: False
- regime: c2a
- n_boot: 20
- n_perm: 20
- site_strat_auc: 0.9391268216259451
- path: `results/posctrl/c2a.json`

### c2b
- boot_median_angle: 49.19344218893373
- achieved_within_site_r2: 0.007017012086602048
- within_site_null: -0.13863231961763728
- loso_r2: 0.03361296488359594
- transfer_H_to_M: 0.05322837514680673
- transfer_M_to_H: 0.027465053219301927
- transfer_H_to_M_null: -0.13791466931193258
- transfer_M_to_H_null: -0.14524753175084928
- young_old_balanced: 85.8956847377361
- pass_angle: False
- pass_transfer: True
- pass_both: False
- regime: c2b
- n_boot: 20
- n_perm: 20
- site_strat_auc: 0.6025712337413494
- path: `results/posctrl/c2b.json`

### c1_dense_shared_r015
- boot_median_angle: 44.171859850003344
- angle_vs_sex_ref: 2.8718598500033465
- sex_angle_ref: 41.3
- recovery_angle: 64.66385756099184
- achieved_within_site_r2: 0.14174260474379086
- within_site_null: -0.13824817897477262
- loso_r2: -0.5550483713501075
- transfer_H_to_M: -0.44981857253997515
- transfer_M_to_H: -0.5359331653787816
- transfer_H_to_M_null: -1.5311938653944783
- transfer_M_to_H_null: -1.629410121176419
- young_old_balanced: 86.485322377569
- pass_angle: False
- pass_transfer: False
- pass_transfer_primary: False
- pass_both: False
- beta: 0.18666913590500295
- target_r2: 0.15
- regime: dense_shared
- n_boot: 10
- n_perm: 10
- path: `results/posctrl/c1_dense_shared_r015.json`

### c1_dense_shared_r030
- boot_median_angle: 40.64229332178904
- angle_vs_sex_ref: -0.6577066782109569
- sex_angle_ref: 41.3
- recovery_angle: 57.71408518361808
- achieved_within_site_r2: 0.29079956034761606
- within_site_null: -0.15925201156702481
- loso_r2: -0.2629843642918844
- transfer_H_to_M: -0.08741143846402355
- transfer_M_to_H: -0.4148615968605822
- transfer_H_to_M_null: -1.6050643361600918
- transfer_M_to_H_null: -1.4459493875383989
- young_old_balanced: 84.74852621025812
- pass_angle: False
- pass_transfer: False
- pass_transfer_primary: False
- pass_both: False
- beta: 0.24683784844272705
- target_r2: 0.3
- regime: dense_shared
- n_boot: 20
- n_perm: 20
- path: `results/posctrl/c1_dense_shared_r030.json`

### c1_dense_shared_r050
- boot_median_angle: 36.21824273532739
- angle_vs_sex_ref: -5.081757264672611
- sex_angle_ref: 41.3
- recovery_angle: 49.831452157934436
- achieved_within_site_r2: 0.5097057772331364
- within_site_null: -0.14138283487623649
- loso_r2: 0.09952068643337378
- transfer_H_to_M: 0.3087454711771438
- transfer_M_to_H: -0.004080027381451634
- transfer_H_to_M_null: -1.474069884908141
- transfer_M_to_H_null: -1.6725931649177501
- young_old_balanced: 80.55488349564683
- pass_angle: False
- pass_transfer: False
- pass_transfer_primary: False
- pass_both: False
- beta: 0.34711538643104906
- target_r2: 0.5
- regime: dense_shared
- n_boot: 10
- n_perm: 10
- path: `results/posctrl/c1_dense_shared_r050.json`

### c1_dense_typespec_r030
- boot_median_angle: 40.680893100514446
- angle_vs_sex_ref: -0.6191068994855513
- sex_angle_ref: 41.3
- recovery_angle: 57.851277447277994
- achieved_within_site_r2: 0.29083609089792695
- within_site_null: -0.15886675604880246
- loso_r2: -0.2851969487447955
- transfer_H_to_M: -0.12149089243929823
- transfer_M_to_H: -0.4343552882051325
- transfer_H_to_M_null: -1.60567518808809
- transfer_M_to_H_null: -1.4482479230528438
- young_old_balanced: 84.81439731631446
- pass_angle: False
- pass_transfer: False
- pass_transfer_primary: False
- pass_both: False
- beta: 0.24479061664766935
- target_r2: 0.3
- regime: dense_typespec
- n_boot: 20
- n_perm: 20
- path: `results/posctrl/c1_dense_typespec_r030.json`

### c1_sparse_shared_r030
- boot_median_angle: 41.00626640175669
- angle_vs_sex_ref: -0.2937335982433069
- sex_angle_ref: 41.3
- recovery_angle: 57.425901780368065
- achieved_within_site_r2: 0.30441941650197835
- within_site_null: -0.15906562227227922
- loso_r2: -0.2493854344626356
- transfer_H_to_M: -0.06492845093976074
- transfer_M_to_H: -0.4136181609414218
- transfer_H_to_M_null: -1.606520591983644
- transfer_M_to_H_null: -1.444969105499693
- young_old_balanced: 83.38258723352521
- pass_angle: False
- pass_transfer: False
- pass_transfer_primary: False
- pass_both: False
- beta: 0.3226823122214063
- target_r2: 0.3
- regime: sparse_shared
- n_boot: 20
- n_perm: 20
- path: `results/posctrl/c1_sparse_shared_r030.json`

### c1_sparse_typespec_r030
- boot_median_angle: 40.45825888932583
- angle_vs_sex_ref: -0.8417411106741639
- sex_angle_ref: 41.3
- recovery_angle: 58.04140995758933
- achieved_within_site_r2: 0.29578407961886966
- within_site_null: -0.15888221038306988
- loso_r2: -0.2600393753250141
- transfer_H_to_M: -0.08446023949278847
- transfer_M_to_H: -0.4150243081474245
- transfer_H_to_M_null: -1.6054839537027359
- transfer_M_to_H_null: -1.448166325515419
- young_old_balanced: 84.66473772766307
- pass_angle: False
- pass_transfer: False
- pass_transfer_primary: False
- pass_both: False
- beta: 0.2667642179173281
- target_r2: 0.3
- regime: sparse_typespec
- n_boot: 20
- n_perm: 20
- path: `results/posctrl/c1_sparse_typespec_r030.json`

### c3_planted_n060_d00
- boot_median_angle: 50.0904365769124
- angle_vs_sex_ref: 8.790436576912406
- sex_angle_ref: 41.3
- recovery_angle: 71.59138963835548
- transfer_H_to_M: -1.475740319227603
- transfer_M_to_H: -1.2187470244988758
- transfer_H_to_M_null: -2.6039176890826203
- transfer_M_to_H_null: -1.6249415157898546
- pass_transfer: False
- pass_transfer_primary: False
- kind: planted
- n_requested: 60
- draw: 0
- n_boot: 10
- n_perm: 10
- path: `results/posctrl/c3_planted_n060_d00.json`

### c3_planted_n060_d01
- boot_median_angle: 45.98375145451136
- angle_vs_sex_ref: 4.683751454511366
- sex_angle_ref: 41.3
- recovery_angle: 71.80833646662064
- transfer_H_to_M: -0.8762386031643874
- transfer_M_to_H: -1.1173779279823193
- transfer_H_to_M_null: -1.1176604486401909
- transfer_M_to_H_null: -1.6248092066304423
- pass_transfer: False
- pass_transfer_primary: False
- kind: planted
- n_requested: 60
- draw: 1
- n_boot: 10
- n_perm: 10
- path: `results/posctrl/c3_planted_n060_d01.json`

### c3_planted_n060_d02
- boot_median_angle: 44.343053005907365
- angle_vs_sex_ref: 3.0430530059073675
- sex_angle_ref: 41.3
- recovery_angle: 70.45499540500208
- transfer_H_to_M: -1.1823380047477894
- transfer_M_to_H: -0.9618679430924912
- transfer_H_to_M_null: -1.722246472556399
- transfer_M_to_H_null: -1.8806490638109388
- pass_transfer: False
- pass_transfer_primary: False
- kind: planted
- n_requested: 60
- draw: 2
- n_boot: 10
- n_perm: 10
- path: `results/posctrl/c3_planted_n060_d02.json`

### c3_planted_n060_d03
- boot_median_angle: 48.019873900757965
- angle_vs_sex_ref: 6.7198739007579675
- sex_angle_ref: 41.3
- recovery_angle: 72.23376489475011
- transfer_H_to_M: -1.6745736863782819
- transfer_M_to_H: -2.235234444288181
- transfer_H_to_M_null: -1.892389395985703
- transfer_M_to_H_null: -2.844295468315394
- pass_transfer: False
- pass_transfer_primary: False
- kind: planted
- n_requested: 60
- draw: 3
- n_boot: 10
- n_perm: 10
- path: `results/posctrl/c3_planted_n060_d03.json`

### c3_planted_n060_d04
- boot_median_angle: 48.2782322095262
- angle_vs_sex_ref: 6.978232209526205
- sex_angle_ref: 41.3
- recovery_angle: 71.47000845042005
- transfer_H_to_M: -0.08167483179465451
- transfer_M_to_H: -0.1279065082931934
- transfer_H_to_M_null: -0.946049163091169
- transfer_M_to_H_null: -0.6051464510960519
- pass_transfer: False
- pass_transfer_primary: False
- kind: planted
- n_requested: 60
- draw: 4
- n_boot: 10
- n_perm: 10
- path: `results/posctrl/c3_planted_n060_d04.json`

### c3_planted_n060_d05
- boot_median_angle: 50.74208997007976
- angle_vs_sex_ref: 9.442089970079763
- sex_angle_ref: 41.3
- recovery_angle: 72.9819725789902
- transfer_H_to_M: -1.0350090529808276
- transfer_M_to_H: -0.7379024001513672
- transfer_H_to_M_null: -1.6960994584207305
- transfer_M_to_H_null: -1.153241102248338
- pass_transfer: False
- pass_transfer_primary: False
- kind: planted
- n_requested: 60
- draw: 5
- n_boot: 10
- n_perm: 10
- path: `results/posctrl/c3_planted_n060_d05.json`

### c3_planted_n060_d06
- boot_median_angle: 52.19731926768233
- angle_vs_sex_ref: 10.897319267682334
- sex_angle_ref: 41.3
- recovery_angle: 72.3996708052366
- transfer_H_to_M: -0.5324384073432573
- transfer_M_to_H: -0.8557194243326407
- transfer_H_to_M_null: -1.3952468900526607
- transfer_M_to_H_null: -1.0053468062248798
- pass_transfer: False
- pass_transfer_primary: False
- kind: planted
- n_requested: 60
- draw: 6
- n_boot: 10
- n_perm: 10
- path: `results/posctrl/c3_planted_n060_d06.json`

### c3_planted_n060_d07
- boot_median_angle: 48.553114074502474
- angle_vs_sex_ref: 7.253114074502477
- sex_angle_ref: 41.3
- recovery_angle: 71.38596699564737
- transfer_H_to_M: -1.2878500954598628
- transfer_M_to_H: -1.5397038631844215
- transfer_H_to_M_null: -2.2324409460884156
- transfer_M_to_H_null: -1.6472849827874008
- pass_transfer: False
- pass_transfer_primary: False
- kind: planted
- n_requested: 60
- draw: 7
- n_boot: 10
- n_perm: 10
- path: `results/posctrl/c3_planted_n060_d07.json`

### c3_planted_n060_d08
- boot_median_angle: 51.45953170686112
- angle_vs_sex_ref: 10.159531706861124
- sex_angle_ref: 41.3
- recovery_angle: 72.48906053660083
- transfer_H_to_M: -0.9127640416927367
- transfer_M_to_H: -1.013596713977367
- transfer_H_to_M_null: -1.6656528276813962
- transfer_M_to_H_null: -1.278793696048361
- pass_transfer: False
- pass_transfer_primary: False
- kind: planted
- n_requested: 60
- draw: 8
- n_boot: 10
- n_perm: 10
- path: `results/posctrl/c3_planted_n060_d08.json`

### c3_planted_n060_d09
- boot_median_angle: 50.00875564397133
- angle_vs_sex_ref: 8.70875564397133
- sex_angle_ref: 41.3
- recovery_angle: 74.16711594335388
- transfer_H_to_M: -0.6935609948882024
- transfer_M_to_H: -1.053505935402093
- transfer_H_to_M_null: -1.4751906608335947
- transfer_M_to_H_null: -1.9524833443171132
- pass_transfer: False
- pass_transfer_primary: False
- kind: planted
- n_requested: 60
- draw: 9
- n_boot: 10
- n_perm: 10
- path: `results/posctrl/c3_planted_n060_d09.json`

### c3_planted_n100_d00
- boot_median_angle: 45.920735653992736
- angle_vs_sex_ref: 4.620735653992739
- sex_angle_ref: 41.3
- recovery_angle: 65.42821975823944
- transfer_H_to_M: -0.21371159422138053
- transfer_M_to_H: -0.2005613970033261
- transfer_H_to_M_null: -1.2083130155882666
- transfer_M_to_H_null: -1.4441870569664244
- pass_transfer: False
- pass_transfer_primary: False
- kind: planted
- n_requested: 100
- draw: 0
- n_boot: 10
- n_perm: 10
- path: `results/posctrl/c3_planted_n100_d00.json`

### c3_planted_n100_d01
- boot_median_angle: 43.895221731786386
- angle_vs_sex_ref: 2.5952217317863884
- sex_angle_ref: 41.3
- recovery_angle: 66.83000787669695
- transfer_H_to_M: -1.0836502381879876
- transfer_M_to_H: -0.8376354427533399
- transfer_H_to_M_null: -1.7635642718426658
- transfer_M_to_H_null: -2.6154088680188496
- pass_transfer: False
- pass_transfer_primary: False
- kind: planted
- n_requested: 100
- draw: 1
- n_boot: 10
- n_perm: 10
- path: `results/posctrl/c3_planted_n100_d01.json`

### c3_planted_n100_d02
- boot_median_angle: 43.39119405088488
- angle_vs_sex_ref: 2.09119405088488
- sex_angle_ref: 41.3
- recovery_angle: 66.81538314566218
- transfer_H_to_M: -1.1837409692791074
- transfer_M_to_H: -1.0066460472398284
- transfer_H_to_M_null: -2.3075712474128656
- transfer_M_to_H_null: -1.59914883345945
- pass_transfer: False
- pass_transfer_primary: False
- kind: planted
- n_requested: 100
- draw: 2
- n_boot: 10
- n_perm: 10
- path: `results/posctrl/c3_planted_n100_d02.json`

### c3_planted_n100_d03
- boot_median_angle: 42.646961384548966
- angle_vs_sex_ref: 1.346961384548969
- sex_angle_ref: 41.3
- recovery_angle: 68.776885093211
- transfer_H_to_M: -0.8074468035510445
- transfer_M_to_H: -1.149528690182914
- transfer_H_to_M_null: -1.27231543169088
- transfer_M_to_H_null: -1.2751435178854094
- pass_transfer: False
- pass_transfer_primary: False
- kind: planted
- n_requested: 100
- draw: 3
- n_boot: 10
- n_perm: 10
- path: `results/posctrl/c3_planted_n100_d03.json`

### c3_planted_n100_d04
- boot_median_angle: 46.155713683319796
- angle_vs_sex_ref: 4.855713683319799
- sex_angle_ref: 41.3
- recovery_angle: 66.76566937004858
- transfer_H_to_M: -1.2847690682429755
- transfer_M_to_H: -1.3092828755215489
- transfer_H_to_M_null: -2.50426311501087
- transfer_M_to_H_null: -2.0904825882075695
- pass_transfer: False
- pass_transfer_primary: False
- kind: planted
- n_requested: 100
- draw: 4
- n_boot: 10
- n_perm: 10
- path: `results/posctrl/c3_planted_n100_d04.json`

### c3_planted_n100_d05
- boot_median_angle: 45.596467592125016
- angle_vs_sex_ref: 4.296467592125019
- sex_angle_ref: 41.3
- recovery_angle: 67.46127893540574
- transfer_H_to_M: -0.9015699651062803
- transfer_M_to_H: -0.8094039625590561
- transfer_H_to_M_null: -1.5484116629474514
- transfer_M_to_H_null: -2.7413732519075795
- pass_transfer: False
- pass_transfer_primary: False
- kind: planted
- n_requested: 100
- draw: 5
- n_boot: 10
- n_perm: 10
- path: `results/posctrl/c3_planted_n100_d05.json`

### c3_planted_n100_d06
- boot_median_angle: 47.42362394356448
- angle_vs_sex_ref: 6.123623943564482
- sex_angle_ref: 41.3
- recovery_angle: 70.8230907745145
- transfer_H_to_M: -0.9372621352438163
- transfer_M_to_H: -1.1497328451000637
- transfer_H_to_M_null: -1.6400909397208217
- transfer_M_to_H_null: -2.105720835525114
- pass_transfer: False
- pass_transfer_primary: False
- kind: planted
- n_requested: 100
- draw: 6
- n_boot: 10
- n_perm: 10
- path: `results/posctrl/c3_planted_n100_d06.json`

### c3_planted_n100_d07
- boot_median_angle: 46.672195941294866
- angle_vs_sex_ref: 5.372195941294869
- sex_angle_ref: 41.3
- recovery_angle: 69.43232372309114
- transfer_H_to_M: -0.7628345729892891
- transfer_M_to_H: -0.6807944142444431
- transfer_H_to_M_null: -1.8422048400183406
- transfer_M_to_H_null: -1.8035520503153055
- pass_transfer: False
- pass_transfer_primary: False
- kind: planted
- n_requested: 100
- draw: 7
- n_boot: 10
- n_perm: 10
- path: `results/posctrl/c3_planted_n100_d07.json`

### c3_planted_n100_d08
- boot_median_angle: 47.51037116751299
- angle_vs_sex_ref: 6.210371167512996
- sex_angle_ref: 41.3
- recovery_angle: 69.13550970835362
- transfer_H_to_M: -1.099836518857101
- transfer_M_to_H: -1.3128641821910214
- transfer_H_to_M_null: -2.845018676430722
- transfer_M_to_H_null: -1.6614371584705112
- pass_transfer: False
- pass_transfer_primary: False
- kind: planted
- n_requested: 100
- draw: 8
- n_boot: 10
- n_perm: 10
- path: `results/posctrl/c3_planted_n100_d08.json`

### c3_planted_n100_d09
- boot_median_angle: 43.70974979916155
- angle_vs_sex_ref: 2.4097497991615526
- sex_angle_ref: 41.3
- recovery_angle: 66.99544347810861
- transfer_H_to_M: -0.35440999376891336
- transfer_M_to_H: -0.8264905341212101
- transfer_H_to_M_null: -1.2892128083360346
- transfer_M_to_H_null: -1.6020761418853002
- pass_transfer: False
- pass_transfer_primary: False
- kind: planted
- n_requested: 100
- draw: 9
- n_boot: 10
- n_perm: 10
- path: `results/posctrl/c3_planted_n100_d09.json`

### c3_planted_n150_d00
- boot_median_angle: 44.74087913432575
- angle_vs_sex_ref: 3.4408791343257548
- sex_angle_ref: 41.3
- recovery_angle: 64.67007086957605
- transfer_H_to_M: -0.059240249666237976
- transfer_M_to_H: -0.7270154836173054
- transfer_H_to_M_null: -1.187443251800601
- transfer_M_to_H_null: -1.18440499140658
- pass_transfer: False
- pass_transfer_primary: False
- kind: planted
- n_requested: 150
- draw: 0
- n_boot: 10
- n_perm: 10
- path: `results/posctrl/c3_planted_n150_d00.json`

### c3_planted_n150_d01
- boot_median_angle: 39.31377023322781
- angle_vs_sex_ref: -1.9862297667721904
- sex_angle_ref: 41.3
- recovery_angle: 62.2313013483502
- transfer_H_to_M: -0.4731115101518364
- transfer_M_to_H: -0.4523089627971696
- transfer_H_to_M_null: -1.9419008148328445
- transfer_M_to_H_null: -1.404861106733088
- pass_transfer: False
- pass_transfer_primary: False
- kind: planted
- n_requested: 150
- draw: 1
- n_boot: 10
- n_perm: 10
- path: `results/posctrl/c3_planted_n150_d01.json`

### c3_planted_n150_d02
- boot_median_angle: 40.65193309709262
- angle_vs_sex_ref: -0.6480669029073738
- sex_angle_ref: 41.3
- recovery_angle: 63.75139319644588
- transfer_H_to_M: -0.16933321550204927
- transfer_M_to_H: -0.364997530601203
- transfer_H_to_M_null: -1.3139530479472943
- transfer_M_to_H_null: -1.4655279404009645
- pass_transfer: False
- pass_transfer_primary: False
- kind: planted
- n_requested: 150
- draw: 2
- n_boot: 10
- n_perm: 10
- path: `results/posctrl/c3_planted_n150_d02.json`

### c3_planted_n150_d03
- boot_median_angle: 43.675151789677884
- angle_vs_sex_ref: 2.3751517896778864
- sex_angle_ref: 41.3
- recovery_angle: 63.95936873999103
- transfer_H_to_M: -0.21556185570991637
- transfer_M_to_H: -0.5156783948539748
- transfer_H_to_M_null: -1.3075779972310941
- transfer_M_to_H_null: -1.9498020585024556
- pass_transfer: False
- pass_transfer_primary: False
- kind: planted
- n_requested: 150
- draw: 3
- n_boot: 10
- n_perm: 10
- path: `results/posctrl/c3_planted_n150_d03.json`

### c3_planted_n150_d04
- boot_median_angle: 43.83894759965693
- angle_vs_sex_ref: 2.538947599656936
- sex_angle_ref: 41.3
- recovery_angle: 64.9518439630448
- transfer_H_to_M: -1.0918414885925514
- transfer_M_to_H: -1.520047658772837
- transfer_H_to_M_null: -2.0242702239584855
- transfer_M_to_H_null: -2.349002614477244
- pass_transfer: False
- pass_transfer_primary: False
- kind: planted
- n_requested: 150
- draw: 4
- n_boot: 10
- n_perm: 10
- path: `results/posctrl/c3_planted_n150_d04.json`

### c3_planted_n150_d05
- boot_median_angle: 42.22615753423671
- angle_vs_sex_ref: 0.9261575342367152
- sex_angle_ref: 41.3
- recovery_angle: 63.14221620964994
- transfer_H_to_M: -0.5080643280931159
- transfer_M_to_H: -0.4896082484450708
- transfer_H_to_M_null: -1.7715145512194848
- transfer_M_to_H_null: -1.7962870347218853
- pass_transfer: False
- pass_transfer_primary: False
- kind: planted
- n_requested: 150
- draw: 5
- n_boot: 10
- n_perm: 10
- path: `results/posctrl/c3_planted_n150_d05.json`

### c3_planted_n150_d06
- boot_median_angle: 44.705904969751586
- angle_vs_sex_ref: 3.405904969751589
- sex_angle_ref: 41.3
- recovery_angle: 64.92566384226852
- transfer_H_to_M: -0.9321427478810143
- transfer_M_to_H: -1.155136668940893
- transfer_H_to_M_null: -1.5116383976195515
- transfer_M_to_H_null: -2.2265858206968487
- pass_transfer: False
- pass_transfer_primary: False
- kind: planted
- n_requested: 150
- draw: 6
- n_boot: 10
- n_perm: 10
- path: `results/posctrl/c3_planted_n150_d06.json`

### c3_planted_n150_d07
- boot_median_angle: 39.53585555437519
- angle_vs_sex_ref: -1.764144445624808
- sex_angle_ref: 41.3
- recovery_angle: 62.22998995731462
- transfer_H_to_M: -0.8275535078103327
- transfer_M_to_H: -0.8850534573931146
- transfer_H_to_M_null: -2.278710108435954
- transfer_M_to_H_null: -1.9746301228106937
- pass_transfer: False
- pass_transfer_primary: False
- kind: planted
- n_requested: 150
- draw: 7
- n_boot: 10
- n_perm: 10
- path: `results/posctrl/c3_planted_n150_d07.json`

### c3_planted_n150_d08
- boot_median_angle: 44.050938587797944
- angle_vs_sex_ref: 2.750938587797947
- sex_angle_ref: 41.3
- recovery_angle: 62.44637451376187
- transfer_H_to_M: -0.17418615618408584
- transfer_M_to_H: -0.13546817793004728
- transfer_H_to_M_null: -1.5140999027949205
- transfer_M_to_H_null: -1.2557680040427461
- pass_transfer: False
- pass_transfer_primary: False
- kind: planted
- n_requested: 150
- draw: 8
- n_boot: 10
- n_perm: 10
- path: `results/posctrl/c3_planted_n150_d08.json`

### c3_planted_n150_d09
- boot_median_angle: 42.20058188262027
- angle_vs_sex_ref: 0.9005818826202727
- sex_angle_ref: 41.3
- recovery_angle: 61.31178370526221
- transfer_H_to_M: -0.7987685308422507
- transfer_M_to_H: -0.8092328865166322
- transfer_H_to_M_null: -2.585706640175746
- transfer_M_to_H_null: -1.8323245861072786
- pass_transfer: False
- pass_transfer_primary: False
- kind: planted
- n_requested: 150
- draw: 9
- n_boot: 10
- n_perm: 10
- path: `results/posctrl/c3_planted_n150_d09.json`

### c3_planted_n233_d00
- boot_median_angle: 41.1913944240051
- angle_vs_sex_ref: -0.1086055759949005
- sex_angle_ref: 41.3
- recovery_angle: 58.25824021790141
- transfer_H_to_M: -0.08741143846402355
- transfer_M_to_H: -0.4148615968605822
- transfer_H_to_M_null: -1.7374077865949669
- transfer_M_to_H_null: -1.6132943150035703
- pass_transfer: False
- pass_transfer_primary: False
- kind: planted
- n_requested: 233
- draw: 0
- n_boot: 10
- n_perm: 10
- path: `results/posctrl/c3_planted_n233_d00.json`

### c3_real_n060_d00
- boot_median_angle: 49.63798359799575
- angle_vs_sex_ref: 8.337983597995752
- sex_angle_ref: 41.3
- transfer_H_to_M: -0.677137418173989
- transfer_M_to_H: -1.768424670926283
- transfer_H_to_M_null: -1.382966599705205
- transfer_M_to_H_null: -1.659036006692509
- pass_transfer: False
- pass_transfer_primary: False
- kind: real
- n_requested: 60
- draw: 0
- n_boot: 20
- n_perm: 20
- path: `results/posctrl/c3_real_n060_d00.json`

### c3_real_n060_d01
- boot_median_angle: 49.3921192930696
- angle_vs_sex_ref: 8.0921192930696
- sex_angle_ref: 41.3
- transfer_H_to_M: -1.0681579264096883
- transfer_M_to_H: -1.1277033138897936
- transfer_H_to_M_null: -2.229204078951778
- transfer_M_to_H_null: -1.9513808780687085
- pass_transfer: False
- pass_transfer_primary: False
- kind: real
- n_requested: 60
- draw: 1
- n_boot: 20
- n_perm: 20
- path: `results/posctrl/c3_real_n060_d01.json`

### c3_real_n060_d02
- boot_median_angle: 49.828011032905515
- angle_vs_sex_ref: 8.528011032905518
- sex_angle_ref: 41.3
- transfer_H_to_M: -2.05440078231991
- transfer_M_to_H: -3.2690503401148883
- transfer_H_to_M_null: -3.101023283965581
- transfer_M_to_H_null: -3.9491639668444813
- pass_transfer: False
- pass_transfer_primary: False
- kind: real
- n_requested: 60
- draw: 2
- n_boot: 20
- n_perm: 20
- path: `results/posctrl/c3_real_n060_d02.json`

### c3_real_n060_d03
- boot_median_angle: 49.141869202888785
- angle_vs_sex_ref: 7.8418692028887875
- sex_angle_ref: 41.3
- transfer_H_to_M: -1.7728506641336192
- transfer_M_to_H: -2.7878314702103513
- transfer_H_to_M_null: -2.444728488383487
- transfer_M_to_H_null: -3.5511908278024285
- pass_transfer: False
- pass_transfer_primary: False
- kind: real
- n_requested: 60
- draw: 3
- n_boot: 20
- n_perm: 20
- path: `results/posctrl/c3_real_n060_d03.json`

### c3_real_n060_d04
- boot_median_angle: 47.43627904036954
- angle_vs_sex_ref: 6.136279040369544
- sex_angle_ref: 41.3
- transfer_H_to_M: -0.5226932649417043
- transfer_M_to_H: -0.5207668017807467
- transfer_H_to_M_null: -1.3616069245069045
- transfer_M_to_H_null: -1.1545066102966337
- pass_transfer: False
- pass_transfer_primary: False
- kind: real
- n_requested: 60
- draw: 4
- n_boot: 20
- n_perm: 20
- path: `results/posctrl/c3_real_n060_d04.json`

### c3_real_n060_d05
- boot_median_angle: 48.40774934267068
- angle_vs_sex_ref: 7.1077493426706795
- sex_angle_ref: 41.3
- transfer_H_to_M: -0.5204100675392063
- transfer_M_to_H: -1.0283087894128433
- transfer_H_to_M_null: -0.7843869446424404
- transfer_M_to_H_null: -1.4839709122488711
- pass_transfer: False
- pass_transfer_primary: False
- kind: real
- n_requested: 60
- draw: 5
- n_boot: 20
- n_perm: 20
- path: `results/posctrl/c3_real_n060_d05.json`

### c3_real_n060_d06
- boot_median_angle: 50.08374008061295
- angle_vs_sex_ref: 8.78374008061295
- sex_angle_ref: 41.3
- transfer_H_to_M: -0.8043805854292813
- transfer_M_to_H: -1.0159400936046712
- transfer_H_to_M_null: -0.9712297463843449
- transfer_M_to_H_null: -1.8257139029466534
- pass_transfer: False
- pass_transfer_primary: False
- kind: real
- n_requested: 60
- draw: 6
- n_boot: 20
- n_perm: 20
- path: `results/posctrl/c3_real_n060_d06.json`

### c3_real_n060_d07
- boot_median_angle: 46.69695576760575
- angle_vs_sex_ref: 5.396955767605753
- sex_angle_ref: 41.3
- transfer_H_to_M: -0.4355861986373868
- transfer_M_to_H: -0.6959889227762192
- transfer_H_to_M_null: -1.6342928621828121
- transfer_M_to_H_null: -0.9682279893444861
- pass_transfer: False
- pass_transfer_primary: False
- kind: real
- n_requested: 60
- draw: 7
- n_boot: 20
- n_perm: 20
- path: `results/posctrl/c3_real_n060_d07.json`

### c3_real_n060_d08
- boot_median_angle: 50.38210100566871
- angle_vs_sex_ref: 9.082101005668711
- sex_angle_ref: 41.3
- transfer_H_to_M: -1.0041897877192674
- transfer_M_to_H: -1.5355974235582712
- transfer_H_to_M_null: -2.188626839873534
- transfer_M_to_H_null: -1.8213465455047604
- pass_transfer: False
- pass_transfer_primary: False
- kind: real
- n_requested: 60
- draw: 8
- n_boot: 20
- n_perm: 20
- path: `results/posctrl/c3_real_n060_d08.json`

### c3_real_n060_d09
- boot_median_angle: 49.36047086254707
- angle_vs_sex_ref: 8.06047086254707
- sex_angle_ref: 41.3
- transfer_H_to_M: -1.32883463268904
- transfer_M_to_H: -2.3708764914092675
- transfer_H_to_M_null: -2.0796230663210644
- transfer_M_to_H_null: -2.3661335652558675
- pass_transfer: False
- pass_transfer_primary: False
- kind: real
- n_requested: 60
- draw: 9
- n_boot: 20
- n_perm: 20
- path: `results/posctrl/c3_real_n060_d09.json`

### c3_real_n100_d00
- boot_median_angle: 46.642780848985865
- angle_vs_sex_ref: 5.342780848985868
- sex_angle_ref: 41.3
- transfer_H_to_M: -1.8648356031192694
- transfer_M_to_H: -1.8679461979616743
- transfer_H_to_M_null: -3.1582891285651375
- transfer_M_to_H_null: -2.833051251835775
- pass_transfer: False
- pass_transfer_primary: False
- kind: real
- n_requested: 100
- draw: 0
- n_boot: 20
- n_perm: 20
- path: `results/posctrl/c3_real_n100_d00.json`

### c3_real_n100_d01
- boot_median_angle: 46.291203188185364
- angle_vs_sex_ref: 4.9912031881853665
- sex_angle_ref: 41.3
- transfer_H_to_M: -0.915180260958548
- transfer_M_to_H: -0.9911111690388451
- transfer_H_to_M_null: -1.527139209435799
- transfer_M_to_H_null: -1.8572996593487292
- pass_transfer: False
- pass_transfer_primary: False
- kind: real
- n_requested: 100
- draw: 1
- n_boot: 20
- n_perm: 20
- path: `results/posctrl/c3_real_n100_d01.json`

### c3_real_n100_d02
- boot_median_angle: 48.8405935306679
- angle_vs_sex_ref: 7.5405935306679055
- sex_angle_ref: 41.3
- transfer_H_to_M: -0.6708158034532319
- transfer_M_to_H: -0.7819778850389255
- transfer_H_to_M_null: -1.32913909177533
- transfer_M_to_H_null: -1.584971358776616
- pass_transfer: False
- pass_transfer_primary: False
- kind: real
- n_requested: 100
- draw: 2
- n_boot: 20
- n_perm: 20
- path: `results/posctrl/c3_real_n100_d02.json`

### c3_real_n100_d03
- boot_median_angle: 46.15658332672952
- angle_vs_sex_ref: 4.85658332672952
- sex_angle_ref: 41.3
- transfer_H_to_M: -0.40888348052207146
- transfer_M_to_H: -0.867338936756895
- transfer_H_to_M_null: -1.0502786585578385
- transfer_M_to_H_null: -1.5941074468690717
- pass_transfer: False
- pass_transfer_primary: False
- kind: real
- n_requested: 100
- draw: 3
- n_boot: 20
- n_perm: 20
- path: `results/posctrl/c3_real_n100_d03.json`

### c3_real_n100_d04
- boot_median_angle: 47.1613829995311
- angle_vs_sex_ref: 5.861382999531102
- sex_angle_ref: 41.3
- transfer_H_to_M: -0.1795328178883724
- transfer_M_to_H: -0.17549389870621845
- transfer_H_to_M_null: -1.1149532403522187
- transfer_M_to_H_null: -1.0171871372174437
- pass_transfer: False
- pass_transfer_primary: False
- kind: real
- n_requested: 100
- draw: 4
- n_boot: 20
- n_perm: 20
- path: `results/posctrl/c3_real_n100_d04.json`

### c3_real_n100_d05
- boot_median_angle: 47.74246645704453
- angle_vs_sex_ref: 6.442466457044532
- sex_angle_ref: 41.3
- transfer_H_to_M: -0.1473223715603379
- transfer_M_to_H: -0.23837278084935853
- transfer_H_to_M_null: -0.9590253526809798
- transfer_M_to_H_null: -0.9237048287044918
- pass_transfer: False
- pass_transfer_primary: False
- kind: real
- n_requested: 100
- draw: 5
- n_boot: 20
- n_perm: 20
- path: `results/posctrl/c3_real_n100_d05.json`

### c3_real_n100_d06
- boot_median_angle: 47.109970545304414
- angle_vs_sex_ref: 5.809970545304417
- sex_angle_ref: 41.3
- transfer_H_to_M: -0.6911257663116419
- transfer_M_to_H: -0.3636270494831101
- transfer_H_to_M_null: -2.0012237966259
- transfer_M_to_H_null: -1.8242885795349664
- pass_transfer: False
- pass_transfer_primary: False
- kind: real
- n_requested: 100
- draw: 6
- n_boot: 20
- n_perm: 20
- path: `results/posctrl/c3_real_n100_d06.json`

### c3_real_n100_d07
- boot_median_angle: 47.62638379438317
- angle_vs_sex_ref: 6.326383794383176
- sex_angle_ref: 41.3
- transfer_H_to_M: -0.9316317808033894
- transfer_M_to_H: -1.238905939569402
- transfer_H_to_M_null: -2.178131668071973
- transfer_M_to_H_null: -1.8647686654184714
- pass_transfer: False
- pass_transfer_primary: False
- kind: real
- n_requested: 100
- draw: 7
- n_boot: 20
- n_perm: 20
- path: `results/posctrl/c3_real_n100_d07.json`

### c3_real_n100_d08
- boot_median_angle: 45.342563986078034
- angle_vs_sex_ref: 4.0425639860780365
- sex_angle_ref: 41.3
- transfer_H_to_M: -0.6666134214785698
- transfer_M_to_H: -0.7102383175608039
- transfer_H_to_M_null: -1.7814209556321376
- transfer_M_to_H_null: -2.373214027306814
- pass_transfer: False
- pass_transfer_primary: False
- kind: real
- n_requested: 100
- draw: 8
- n_boot: 20
- n_perm: 20
- path: `results/posctrl/c3_real_n100_d08.json`

### c3_real_n100_d09
- boot_median_angle: 44.520127884846595
- angle_vs_sex_ref: 3.2201278848465975
- sex_angle_ref: 41.3
- transfer_H_to_M: -2.1125678669703847
- transfer_M_to_H: -1.3676606241250773
- transfer_H_to_M_null: -3.4943030909912274
- transfer_M_to_H_null: -2.3948152883181386
- pass_transfer: False
- pass_transfer_primary: False
- kind: real
- n_requested: 100
- draw: 9
- n_boot: 20
- n_perm: 20
- path: `results/posctrl/c3_real_n100_d09.json`

### c3_real_n150_d00
- boot_median_angle: 45.72520272064433
- angle_vs_sex_ref: 4.425202720644336
- sex_angle_ref: 41.3
- transfer_H_to_M: -0.5583191751914665
- transfer_M_to_H: -0.840345377189119
- transfer_H_to_M_null: -1.6870201639173998
- transfer_M_to_H_null: -1.9176173138821646
- pass_transfer: False
- pass_transfer_primary: False
- kind: real
- n_requested: 150
- draw: 0
- n_boot: 20
- n_perm: 20
- path: `results/posctrl/c3_real_n150_d00.json`

### c3_real_n150_d01
- boot_median_angle: 44.634830603372144
- angle_vs_sex_ref: 3.334830603372147
- sex_angle_ref: 41.3
- transfer_H_to_M: -0.3947368942042324
- transfer_M_to_H: -0.6443208011647114
- transfer_H_to_M_null: -1.4957963063431539
- transfer_M_to_H_null: -1.5159833869081032
- pass_transfer: False
- pass_transfer_primary: False
- kind: real
- n_requested: 150
- draw: 1
- n_boot: 20
- n_perm: 20
- path: `results/posctrl/c3_real_n150_d01.json`

### c3_real_n150_d02
- boot_median_angle: 45.332070953473576
- angle_vs_sex_ref: 4.032070953473578
- sex_angle_ref: 41.3
- transfer_H_to_M: -0.8539605629221725
- transfer_M_to_H: -1.1210001588971443
- transfer_H_to_M_null: -1.9012947511939475
- transfer_M_to_H_null: -1.8489436916244402
- pass_transfer: False
- pass_transfer_primary: False
- kind: real
- n_requested: 150
- draw: 2
- n_boot: 20
- n_perm: 20
- path: `results/posctrl/c3_real_n150_d02.json`

### c3_real_n150_d03
- boot_median_angle: 45.50975130894723
- angle_vs_sex_ref: 4.209751308947233
- sex_angle_ref: 41.3
- transfer_H_to_M: -0.48518171478914196
- transfer_M_to_H: -0.486372967056402
- transfer_H_to_M_null: -1.4257992329652978
- transfer_M_to_H_null: -1.5638037370427764
- pass_transfer: False
- pass_transfer_primary: False
- kind: real
- n_requested: 150
- draw: 3
- n_boot: 20
- n_perm: 20
- path: `results/posctrl/c3_real_n150_d03.json`

### c3_real_n150_d04
- boot_median_angle: 44.756426574302935
- angle_vs_sex_ref: 3.456426574302938
- sex_angle_ref: 41.3
- transfer_H_to_M: -0.19200159990432253
- transfer_M_to_H: -0.3346648114850548
- transfer_H_to_M_null: -1.4346465787680187
- transfer_M_to_H_null: -1.472856150920567
- pass_transfer: False
- pass_transfer_primary: False
- kind: real
- n_requested: 150
- draw: 4
- n_boot: 20
- n_perm: 20
- path: `results/posctrl/c3_real_n150_d04.json`

### c3_real_n150_d05
- boot_median_angle: 46.35166574213108
- angle_vs_sex_ref: 5.051665742131085
- sex_angle_ref: 41.3
- transfer_H_to_M: -0.06523525676410924
- transfer_M_to_H: -0.3379572108125384
- transfer_H_to_M_null: -1.2187560727386044
- transfer_M_to_H_null: -1.0955639813646605
- pass_transfer: False
- pass_transfer_primary: False
- kind: real
- n_requested: 150
- draw: 5
- n_boot: 20
- n_perm: 20
- path: `results/posctrl/c3_real_n150_d05.json`

### c3_real_n150_d06
- boot_median_angle: 46.27087885228402
- angle_vs_sex_ref: 4.970878852284024
- sex_angle_ref: 41.3
- transfer_H_to_M: -0.2606420185929563
- transfer_M_to_H: -0.48146086630479445
- transfer_H_to_M_null: -1.4761591538837084
- transfer_M_to_H_null: -1.499466188169135
- pass_transfer: False
- pass_transfer_primary: False
- kind: real
- n_requested: 150
- draw: 6
- n_boot: 20
- n_perm: 20
- path: `results/posctrl/c3_real_n150_d06.json`

### c3_real_n150_d07
- boot_median_angle: 45.48764883767355
- angle_vs_sex_ref: 4.1876488376735495
- sex_angle_ref: 41.3
- transfer_H_to_M: -0.15586853937160194
- transfer_M_to_H: -0.6112350594258662
- transfer_H_to_M_null: -1.1845499716033125
- transfer_M_to_H_null: -1.6354078455950656
- pass_transfer: False
- pass_transfer_primary: False
- kind: real
- n_requested: 150
- draw: 7
- n_boot: 20
- n_perm: 20
- path: `results/posctrl/c3_real_n150_d07.json`

### c3_real_n150_d08
- boot_median_angle: 47.43016752196102
- angle_vs_sex_ref: 6.13016752196102
- sex_angle_ref: 41.3
- transfer_H_to_M: -0.1472925996275909
- transfer_M_to_H: -0.21472518174091526
- transfer_H_to_M_null: -1.0863757099284652
- transfer_M_to_H_null: -1.177114329029662
- pass_transfer: False
- pass_transfer_primary: False
- kind: real
- n_requested: 150
- draw: 8
- n_boot: 20
- n_perm: 20
- path: `results/posctrl/c3_real_n150_d08.json`

### c3_real_n150_d09
- boot_median_angle: 45.28434514962676
- angle_vs_sex_ref: 3.984345149626762
- sex_angle_ref: 41.3
- transfer_H_to_M: -0.46939392271550795
- transfer_M_to_H: -0.9691178640805016
- transfer_H_to_M_null: -1.683415645927828
- transfer_M_to_H_null: -2.2174195803386474
- pass_transfer: False
- pass_transfer_primary: False
- kind: real
- n_requested: 150
- draw: 9
- n_boot: 20
- n_perm: 20
- path: `results/posctrl/c3_real_n150_d09.json`

### c3_real_n233_d00
- boot_median_angle: 45.07377681100883
- angle_vs_sex_ref: 3.7737768110088297
- sex_angle_ref: 41.3
- transfer_H_to_M: -0.12401031967221865
- transfer_M_to_H: -0.463305482028161
- transfer_H_to_M_null: -1.4101336187131237
- transfer_M_to_H_null: -1.6022160208464342
- pass_transfer: False
- pass_transfer_primary: False
- kind: real
- n_requested: 233
- draw: 0
- n_boot: 20
- n_perm: 20
- path: `results/posctrl/c3_real_n233_d00.json`

## Paths written

- `results/posctrl/c1_dense_shared_r015.json`
- `results/posctrl/c1_dense_shared_r030.json`
- `results/posctrl/c1_dense_shared_r050.json`
- `results/posctrl/c1_dense_typespec_r030.json`
- `results/posctrl/c1_sparse_shared_r030.json`
- `results/posctrl/c1_sparse_typespec_r030.json`
- `results/posctrl/c2a.json`
- `results/posctrl/c2b.json`
- `results/posctrl/c3_planted_n060_d00.json`
- `results/posctrl/c3_planted_n060_d01.json`
- `results/posctrl/c3_planted_n060_d02.json`
- `results/posctrl/c3_planted_n060_d03.json`
- `results/posctrl/c3_planted_n060_d04.json`
- `results/posctrl/c3_planted_n060_d05.json`
- `results/posctrl/c3_planted_n060_d06.json`
- `results/posctrl/c3_planted_n060_d07.json`
- `results/posctrl/c3_planted_n060_d08.json`
- `results/posctrl/c3_planted_n060_d09.json`
- `results/posctrl/c3_planted_n100_d00.json`
- `results/posctrl/c3_planted_n100_d01.json`
- `results/posctrl/c3_planted_n100_d02.json`
- `results/posctrl/c3_planted_n100_d03.json`
- `results/posctrl/c3_planted_n100_d04.json`
- `results/posctrl/c3_planted_n100_d05.json`
- `results/posctrl/c3_planted_n100_d06.json`
- `results/posctrl/c3_planted_n100_d07.json`
- `results/posctrl/c3_planted_n100_d08.json`
- `results/posctrl/c3_planted_n100_d09.json`
- `results/posctrl/c3_planted_n150_d00.json`
- `results/posctrl/c3_planted_n150_d01.json`
- `results/posctrl/c3_planted_n150_d02.json`
- `results/posctrl/c3_planted_n150_d03.json`
- `results/posctrl/c3_planted_n150_d04.json`
- `results/posctrl/c3_planted_n150_d05.json`
- `results/posctrl/c3_planted_n150_d06.json`
- `results/posctrl/c3_planted_n150_d07.json`
- `results/posctrl/c3_planted_n150_d08.json`
- `results/posctrl/c3_planted_n150_d09.json`
- `results/posctrl/c3_planted_n233_d00.json`
- `results/posctrl/c3_real_n060_d00.json`
- `results/posctrl/c3_real_n060_d01.json`
- `results/posctrl/c3_real_n060_d02.json`
- `results/posctrl/c3_real_n060_d03.json`
- `results/posctrl/c3_real_n060_d04.json`
- `results/posctrl/c3_real_n060_d05.json`
- `results/posctrl/c3_real_n060_d06.json`
- `results/posctrl/c3_real_n060_d07.json`
- `results/posctrl/c3_real_n060_d08.json`
- `results/posctrl/c3_real_n060_d09.json`
- `results/posctrl/c3_real_n100_d00.json`
- `results/posctrl/c3_real_n100_d01.json`
- `results/posctrl/c3_real_n100_d02.json`
- `results/posctrl/c3_real_n100_d03.json`
- `results/posctrl/c3_real_n100_d04.json`
- `results/posctrl/c3_real_n100_d05.json`
- `results/posctrl/c3_real_n100_d06.json`
- `results/posctrl/c3_real_n100_d07.json`
- `results/posctrl/c3_real_n100_d08.json`
- `results/posctrl/c3_real_n100_d09.json`
- `results/posctrl/c3_real_n150_d00.json`
- `results/posctrl/c3_real_n150_d01.json`
- `results/posctrl/c3_real_n150_d02.json`
- `results/posctrl/c3_real_n150_d03.json`
- `results/posctrl/c3_real_n150_d04.json`
- `results/posctrl/c3_real_n150_d05.json`
- `results/posctrl/c3_real_n150_d06.json`
- `results/posctrl/c3_real_n150_d07.json`
- `results/posctrl/c3_real_n150_d08.json`
- `results/posctrl/c3_real_n150_d09.json`
- `results/posctrl/c3_real_n233_d00.json`
- `results/posctrl/SUPERSESSION_20260916.flag`
- `results/posctrl/supersession_20260916.json`
- `results/posctrl/STOP_P1.txt`
- `FINDINGS_POSCTRL.md`
- `PROGRESS_POSCTRL.md`

## Bars (frozen)

### Original (P1 record; 35° retired for C1/C3)
- angle: bootstrap median pairwise ≤ 35°
- planted transfer (P2): HBCC→MSSM > 0 AND MSSM→HBCC > 0
- C2a transfer (P1 / L2): at least one direction > 0
- calibration: within-site site-stratified ridge R² = 0.30 ± 0.03

### Supersession 2026-09-16 (C1/C3)
- (a) primary: transfer R² > 0 in both directions with permutation null ≤ 0.05
- (b) bootstrap median pairwise angle reported against sex reference 41.3°, not 35°
- P2 = planted dense_shared at R²≈0.30 transfers both ways
- P3 = it does not
- calibration unchanged: within-site site-stratified ridge R² = 0.30 ± 0.03

Fired (supersession): planted does not transfer either. Bank covariance differences defeat transfer for any direction at this n; the age negative is uninformative and the power curve is the result.

### Supersession 2 — 2026-09-17 (C1/C3 transfer reporting)
- transfer also reported as held-out-bank Pearson r, Spearman ρ, and calibrated R² (intercept+slope refit on the test bank)
- each with the same donor-level permutation null
- Pass = r > 0 with null ≤ 0.05 in both directions
- planted = C1 dense_shared at R²≈0.30; real age = gene_ridge chronological age
- no bootstrap; no C3 refit; C3 transfer-only at B/2 if predictions not saved

Fired (supersession 2): planted r>0 and real age r>0 both ways. Direction transfers, scale does not; the LOWDIM negative is a calibration failure.

### Supersession 3 — 2026-09-17
- (1) n_perm=200 on gene_ridge, C2a, dense_shared @0.30; donor-bootstrap 95% CI on r, B=200, resample the test bank
- (2) share = real-age r / planted-@0.30 r per direction, paired bootstrap CI
- (3) identity-residualized age fit transfer; pass = r > 0 with null ≤ 0.05 both ways
- permutation seed 20260914; bootstrap seed 20260918
- no bootstrap-angle; no C3 refit; S2 bar not retuned

Fired (supersession 3): (3) passes. The age axis transfers after removing the identity subspace.

## S2 finished cells (2026-09-17)

Transfer only. No bootstrap. C3 predictions were not saved; C3 transfer rerun at B/2.

### s2_gene_ridge
- r2_H_to_M: -0.12401031967221865
- r2_M_to_H: -0.463305482028161
- cal_r2_H_to_M: 0.3416253228027104
- cal_r2_M_to_H: 0.334092904334764
- r_H_to_M: 0.5681984972828398
- r_M_to_H: 0.5780077026604088
- rho_H_to_M: 0.5734828038569333
- rho_M_to_H: 0.5667128038080262
- r_H_to_M_null: -0.005143838892140646
- r_M_to_H_null: -0.02167705162433014
- pass_r: True
- n_perm: 20
- path: `results/posctrl/s2_gene_ridge.json`

### s2_c2a
- r2_H_to_M: 0.5117665952957187
- r2_M_to_H: 0.35822098443172434
- cal_r2_H_to_M: 0.7157130261395845
- cal_r2_M_to_H: 0.6122296854210548
- r_H_to_M: 0.8459897392663815
- r_M_to_H: 0.7824510754168936
- pass_r: True
- n_perm: 20
- path: `results/posctrl/s2_c2a.json`

### s2_c1_dense_shared_r030
- r2_H_to_M: -0.08741143846402355
- r2_M_to_H: -0.4148615968605822
- cal_r2_H_to_M: 0.7479063157375722
- cal_r2_M_to_H: 0.5896237818883039
- r_H_to_M: 0.8648156017122377
- r_M_to_H: 0.7678696386030015
- r_H_to_M_null: 0.015358459749279787
- r_M_to_H_null: 0.013189263548918983
- pass_r: True
- n_perm: 20
- path: `results/posctrl/s2_c1_dense_shared_r030.json`

### s2_c1 remaining @ 0.30 / 0.15 / 0.50
- dense_typespec_r030 pass_r=True r H→M=+0.853 M→H=+0.750
- sparse_shared_r030 pass_r=True r H→M=+0.865 M→H=+0.768
- sparse_typespec_r030 pass_r=True r H→M=+0.871 M→H=+0.725
- dense_shared_r015 pass_r=True r H→M=+0.740 M→H=+0.584 n_perm=10
- dense_shared_r050 pass_r=True r H→M=+0.952 M→H=+0.910 n_perm=10

### s2_c3 (median over draws; pass = n_pass/n_draws; B/2)
- planted n=60: r +0.375 / +0.334  6/10
- planted n=100: r +0.572 / +0.462  7/10
- planted n=150: r +0.753 / +0.627  10/10
- planted n=233: r +0.865 / +0.768  1/1
- real n=60: r +0.199 / +0.255  8/10
- real n=100: r +0.489 / +0.362  7/10
- real n=150: r +0.479 / +0.512  8/10
- real n=233: r +0.568 / +0.578  0/1 (r_null M→H +0.068 > 0.05 at B/2; gene_ridge at n_perm=20 is the real-age cell)

## S2 paths written
- `results/posctrl/s2_gene_ridge.json`
- `results/posctrl/s2_c2a.json`
- `results/posctrl/s2_c1_dense_shared_r015.json`
- `results/posctrl/s2_c1_dense_shared_r030.json`
- `results/posctrl/s2_c1_dense_shared_r050.json`
- `results/posctrl/s2_c1_dense_typespec_r030.json`
- `results/posctrl/s2_c1_sparse_shared_r030.json`
- `results/posctrl/s2_c1_sparse_typespec_r030.json`
- `results/posctrl/s2_c3_*.json` (62 draw cells)
- `results/posctrl/s2_report.txt`
- `FINDINGS_POSCTRL.md`
- `PROGRESS_POSCTRL.md`

## S3 finished cells (2026-09-17)

Transfer only. n_perm=200. Donor-bootstrap B=200, seed `20260918`. No bootstrap-angle. No C3 refit.

### s3_gene_ridge
- r2_H_to_M: -0.12401031967221865
- r2_M_to_H: -0.463305482028161
- cal_r2_H_to_M: 0.3416253228027104
- cal_r2_M_to_H: 0.334092904334764
- r_H_to_M: 0.5681984972828398
- r_M_to_H: 0.5780077026604088
- rho_H_to_M: 0.5734828038569333
- rho_M_to_H: 0.5667128038080262
- r_H_to_M_null: -0.009656105344413166
- r_M_to_H_null: 0.0037189950233844755
- r_H_to_M_boot_lo: 0.42536547354737786
- r_H_to_M_boot_hi: 0.6731965395301741
- r_M_to_H_boot_lo: 0.47789303329237365
- r_M_to_H_boot_hi: 0.6699703420212719
- pass_r: True
- n_perm: 200
- n_boot_r: 200
- path: `results/posctrl/s3_gene_ridge.json`

### s3_c2a
- r2_H_to_M: 0.5117665952957187
- r2_M_to_H: 0.35822098443172434
- cal_r2_H_to_M: 0.7157130261395845
- cal_r2_M_to_H: 0.6122296854210548
- r_H_to_M: 0.8459897392663815
- r_M_to_H: 0.7824510754168936
- r_H_to_M_null: 0.004587423365229684
- r_M_to_H_null: -0.0021821481457122683
- r_H_to_M_boot_lo: 0.7814594003578668
- r_H_to_M_boot_hi: 0.8979655752921772
- r_M_to_H_boot_lo: 0.7213132451743297
- r_M_to_H_boot_hi: 0.8505780624619199
- pass_r: True
- n_perm: 200
- n_boot_r: 200
- path: `results/posctrl/s3_c2a.json`

### s3_c1_dense_shared_r030
- r2_H_to_M: -0.08741143846402355
- r2_M_to_H: -0.4148615968605822
- cal_r2_H_to_M: 0.7479063157375722
- cal_r2_M_to_H: 0.5896237818883039
- r_H_to_M: 0.8648156017122377
- r_M_to_H: 0.7678696386030015
- r_H_to_M_null: 0.005221340038184346
- r_M_to_H_null: 0.0062082998337348
- r_H_to_M_boot_lo: 0.815998182092321
- r_H_to_M_boot_hi: 0.9088001673551334
- r_M_to_H_boot_lo: 0.6747511345745536
- r_M_to_H_boot_hi: 0.812949829547625
- pass_r: True
- n_perm: 200
- n_boot_r: 200
- path: `results/posctrl/s3_c1_dense_shared_r030.json`

### s3_share (fraction of transferable age signal)
- share_H_to_M: 0.6570169365097839
- share_H_to_M_boot_lo: 0.4990788564519759
- share_H_to_M_boot_hi: 0.7659189615159108
- share_M_to_H: 0.7527419676495977
- share_M_to_H_boot_lo: 0.6418621179358367
- share_M_to_H_boot_hi: 0.9072215694359107
- path: `results/posctrl/s3_share.json`

### s3_idresid
- r2_H_to_M: -0.14912718296696337
- r2_M_to_H: -0.4806828482152581
- cal_r2_H_to_M: 0.32875971002783555
- cal_r2_M_to_H: 0.33175170268216336
- r_H_to_M: 0.5733756975762597
- r_M_to_H: 0.5759789081921002
- rho_H_to_M: 0.574325951946812
- rho_M_to_H: 0.5587521507290429
- r_H_to_M_null: -0.00904514903227908
- r_M_to_H_null: 0.007228155059873649
- pass_r: True
- n_perm: 200
- path: `results/posctrl/s3_idresid.json`

## S3 paths written
- `results/posctrl/s3_gene_ridge.json`
- `results/posctrl/s3_c2a.json`
- `results/posctrl/s3_c1_dense_shared_r030.json`
- `results/posctrl/s3_idresid.json`
- `results/posctrl/s3_share.json`
- `results/posctrl/s3_summary.json`
- `results/posctrl/s3_report.txt`
- `results/posctrl/SUPERSESSION3_20260917.flag`
- `FINDINGS_POSCTRL.md`
- `PROGRESS_POSCTRL.md`

