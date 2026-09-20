# PI Gate C — final bounded electrostatic-probe gate

## Purpose

Gate C is the **last planned scientific addition before paper closure**.

The paper premise remains:

> repulsive WCA baseline -> attraction creates CO2 preference -> framework-specific electrostatic differentiation.

Gates A and B already completed the PI's existing-data requests:
- continuous mechanism vector;
- endpoint-matched / mechanism-distant pairs;
- mechanism-distance distribution;
- expanded rank consequences;
- framework bootstrap;
- predictability gap;
- leverage controls;
- geometry-conditioned shuffle as a retained non-supporting/bounding diagnostic.

Gate C addresses the one remaining positive mechanistic question:

> What physically defined electrostatic property of the accessible pore explains the framework-to-framework variation in Delta ln S_elec?

This is **not** a new open-ended phase.

## Frozen calculation

Population: the same 238-MOF common cohort.

Temperature and force field: exactly the production 298 K rigid-framework UFF + ARC-MOF REPEAT-charge model.

Electrostatic convention: identical full-mode Ewald convention used in the production Henry calculations.

For each framework and each gas (CO2 and N2):

1. Generate uniform random center-of-mass positions in the replicated production cell.
2. At each position evaluate 16 pre-generated deterministic molecular orientations.
3. Evaluate the same pose twice:
   - LJ-only;
   - full.
4. Define the pose-level electrostatic stabilization coordinate

   z_elec = -beta (U_full - U_LJ)

   so positive values mean electrostatic stabilization.
5. A center-of-mass position is accepted as pore-accessible when its orientation-averaged LJ Boltzmann weight is finite and strictly positive. Continue until 6,000 accepted positions are collected.
6. Position-level statistical weighting uses the orientation-averaged LJ Boltzmann weight. No arbitrary pore-energy threshold is introduced.

Random seed: 20260920.

## Frozen implementation clarification

Before any Gate C result is generated, the orientation-level definition is fixed as follows.

For position r with the same 16 molecular axes in LJ-only and full modes,

z_pose(r,Omega) = log w_full(r,Omega) - log w_LJ(r,Omega).

The primary position-level electrostatic coordinate is the orientation-integrated local increment

z_eff(r) = log mean_Omega[w_full(r,Omega)] - log mean_Omega[w_LJ(r,Omega)].

This is the same full-minus-LJ log-weight construction used in the retained spatial analysis, now sampled sparsely across all 238 frameworks. Position weights are proportional to mean_Omega[w_LJ(r,Omega)]. The orientation-sensitivity descriptor is the within-position IQR of z_pose, summarized across positions with the same LJ position weights.

No result was inspected before fixing this implementation.

## Frozen probe descriptors

For each gas:
- weighted median z_elec;
- weighted 90th percentile z_elec (favorable tail);
- weighted IQR z_elec;
- weighted fraction z_elec > 1 (electrostatic stabilization greater than kBT);
- weighted median within-position orientation IQR.

Cross-gas coordinates:
- CO2 minus N2 weighted-median stabilization;
- CO2 minus N2 favorable-tail stabilization.

Total probe feature count: 12.

No additional charge statistics, RACs, GNNs, topology labels, or automatic feature generation are allowed in Gate C.

## Numerical preflight

Before all 238 frameworks, the exact implementation must reproduce the direction of the already-retained M029/M085 electrostatic contrast and pass a fixed 8-framework convergence panel.

Convergence check:
- compare 3,000 vs 6,000 accepted positions using the same nested sample stream;
- every primary descriptor must have absolute standardized change <= 0.10 after scaling across the 8-framework panel;
- favorable-tail and median CO2-N2 gaps must preserve framework ordering to Spearman rho >= 0.90.

If this fails, increase only to 10,000 accepted positions once. No further tuning.

## Frozen model comparison

Target: Delta ln S_elec only.

Use exactly the same five-fold shuffled split rule and ExtraTrees settings as the existing geometry control:
- 5 folds;
- shuffle=True;
- random_state=20260918;
- 500 trees;
- min_samples_leaf=3;
- max_features=0.7.

Compare only:
1. geometry-only (14 frozen descriptors);
2. probe-only (12 frozen descriptors);
3. geometry + probe.

Metrics:
- out-of-fold R2;
- out-of-fold MAE;
- median-predictor MAE.

Paired framework bootstrap: 5,000 resamples of the frozen OOF predictions.

## Pre-specified decision rule

**MAIN-TEXT PASS** only if all conditions hold:

1. combined R2 >= geometry R2 + 0.10;
2. combined MAE <= 0.90 * geometry MAE;
3. bootstrap 95% CI for (MAE_geometry - MAE_combined) is entirely > 0;
4. probe-only R2 >= geometry R2 + 0.05;
5. probe-only MAE <= 0.95 * geometry MAE.

If all five pass:
- Figure 3 becomes a two-coordinate physical design figure;
- the manuscript may state that an accessible-pore electrostatic probe recovers information missing from conventional geometry.

If any condition fails:
- Gate C is reported in SI as non-supporting/partial;
- the existing manuscript claim remains unchanged;
- **the project moves to paper closure.**

## Hard stop

Gate C does not automatically open:
- a second force-field campaign;
- CO2/CH4 transfer;
- GNNs;
- larger descriptor searches;
- finite-pressure simulations;
- new mechanism classes.

Gates D and E remain closed unless the PI or a reviewer explicitly requests them.

This rule exists to finish the PI premise rather than turn it into a new research program.

## Closure outcome

Gate C is **CLOSED AS NON-SUPPORTING AT NUMERICAL PREFLIGHT**.

The frozen 3,000 -> 6,000 accepted-position check failed the all-descriptor convergence criterion, so the one permitted retry at 10,000 positions was run. The 6,000 -> 10,000 comparison still failed the pre-specified requirement that every primary descriptor have absolute standardized change <= 0.10. Only the CO2 favorable-tail descriptor passed that threshold.

Two directional checks remained stable and are recorded only as diagnostics:
- CO2-N2 median-gap ordering: Spearman rho = 0.976 between 6,000 and 10,000 positions;
- CO2-N2 favorable-tail-gap ordering: Spearman rho = 0.976;
- M085 > M029 for both median-gap and favorable-tail-gap coordinates.

Per the frozen hard-stop rule, no 238-framework production probe campaign, downstream probe-vs-geometry model comparison, descriptor tuning, larger sample count, or replacement probe definition is allowed. The existing manuscript claim remains unchanged and the project moves to paper closure.
