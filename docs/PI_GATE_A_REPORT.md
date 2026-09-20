# PI Gate A — existing-data results

Status: **PASS**

Branch: `pi-additions`

No new molecular simulation was used. Results below were computed from the frozen 238-framework authoritative table and the existing geometry out-of-fold predictions.

## 1. Continuous mechanism fingerprint

For each framework:
`m = (Delta ln S_attr, Delta ln S_elec)`

Population summaries:
- median mechanism magnitude = **2.586**
- median mechanism angle = **0.128 rad**
- median electrostatic share = **0.114**

Interpretation:
the attractive increment dominates the magnitude for most frameworks, consistent with the sequential population story, while the electrostatic fraction varies enough to create important framework-level differentiation.

## 2. Performance twins / mechanism strangers

Pre-specified selection rule:
1. compute all **28,203** framework pairs;
2. retain the closest **5%** in final endpoint difference `|Delta ln S_full|`;
3. rank only that frozen subset by Euclidean mechanism distance in `(Delta ln S_attr, Delta ln S_elec)`.

This retains **1,411** endpoint-matched pairs. The endpoint-difference cutoff is **0.05384** in `ln S_full`.

Mechanism distance among these endpoint-matched pairs:
- median = **0.189**
- 90th percentile = **0.518**
- 95th percentile = **0.750**
- 99th percentile = **1.534**

The strongest pair has:
- `|Delta ln S_full| = 0.0160`
- mechanism distance = **2.585**
- attractive increments = **2.211 vs 4.238**
- electrostatic increments = **2.001 vs 0.396**

This is strong evidence that nearly identical final selectivity can be assembled from very different balances of attraction and electrostatics.

Important framing:
this is a **secondary mechanistic result**. It does not replace the main population-level sequence:
repulsive WCA baseline -> attractive creation of CO2 preference -> framework-specific electrostatic differentiation.

## 3. Ranking consequences — upgraded population statistics

### WCA -> LJ-only
- Spearman rho = **-0.751**
- Kendall tau = **-0.560**
- median absolute rank shift = **102**
- >=10% of rank range (24 positions): **213/238 = 89.5%**
- >=25% (60 positions): **173/238 = 72.7%**
- >=50% (119 positions): **98/238 = 41.2%**
- top-10 overlap = **0/10**
- top-25 overlap = **0/25**
- top-10% overlap = **0/24**

This is unusually strong population-level evidence that restoring attraction does not merely perturb the repulsive ordering; it reconstructs the candidate landscape.

### LJ-only -> full
- Spearman rho = **0.905**
- Kendall tau = **0.801**
- median absolute rank shift = **9**
- >=10% of rank range: **41/238 = 17.2%**
- >=25%: **9/238 = 3.8%**
- >=50%: **6/238 = 2.5%**
- top-10 overlap = **6/10**
- top-25 overlap = **18/25**
- top-10% overlap = **17/24**

Interpretation:
electrostatics mostly preserves the attraction-established global ordering but creates a smaller set of large framework-specific exceptions.

## 4. Bootstrap uncertainty

A deterministic **5,000-resample framework bootstrap** was applied to the existing framework-level quantities. For geometry predictability, this is a bootstrap of the already frozen out-of-fold predictions; the model is not refit inside each bootstrap.

### Central physical quantities
- median Delta ln S_attr = **2.537**
  - 95% CI: **[2.437, 2.624]**
- median Delta ln S_elec = **0.309**
  - 95% CI: **[0.294, 0.323]**
- Spearman rho(attr, elec) = **-0.089**
  - 95% CI: **[-0.224, 0.044]**
- fraction with positive attractive increment = **1.000**
  - 95% CI: **[1.000, 1.000]**
- fraction with positive electrostatic increment = **1.000**
  - 95% CI: **[1.000, 1.000]**

The correlation interval includes zero, which reinforces the intended wording: **weakly associated**, not independent.

### Geometry predictability gap
Existing OOF geometry-only results:
- R2_attr = **0.874**
- R2_elec = **0.188**
- gap G = R2_attr - R2_elec = **0.686**

Fixed-OOF bootstrap:
- R2_attr 95% CI: **[0.836, 0.905]**
- R2_elec 95% CI: **[-0.009, 0.298]**
- predictability-gap 95% CI: **[0.563, 0.893]**

This is a strong statistical control supporting the qualitative asymmetry between attraction and electrostatic differentiation.

## 5. Additional leverage controls

Removing the top 5% of the attractive increment:
- rho(attr, elec) = **-0.092**

Removing the top 5% of the electrostatic increment:
- rho = **-0.084**

Removing the top 5% of final selectivity:
- rho = **-0.179**

Removing the ten largest attraction-driven rank movers:
- rho = **-0.105**

Removing the ten largest electrostatic rank movers:
- rho = **-0.083**

The weak attraction/electrostatic association is therefore not created by a few spectacular outliers.

## Decision after Gate A

The PI roadmap passes its first decision gate.

The existing 238-framework data already support two additions worth keeping:

1. **endpoint-mechanism non-uniqueness is real and population-level enough to quantify**, not just a hand-picked pair;
2. the geometry-predictability asymmetry has a large bootstrap-supported gap.

These strengthen the paper without changing its premise.

### Recommended manuscript effect

- Figure 1: unchanged.
- Figure 2: add endpoint-twin/mechanism-stranger information, but preserve the continuous mechanism map.
- Figure 3: retain geometry as a mechanistic control and add the predictability-gap CI.
- Figure 4: keep M029/M085 as the geometry-matched electrostatic counterexample.
- Figure 5: add the top-set-overlap / population rank-movement statistics.

## Next gate

Do **not** start CO2/CH4 or a new force-field campaign yet.

Next:
1. regenerate the Gate A source-data files with `analysis/final/39_pi_existing_data_additions.py`;
2. inspect whether the endpoint-twin result can be presented without confusing the sequential population mechanism;
3. then decide whether the PI's proposed accessible-pore electrostatic descriptor across all 238 frameworks is worth opening as the one genuinely new calculation.
