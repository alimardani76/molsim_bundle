# Scientific freeze — pre-PI publication version

## Central thesis

Across the 238-MOF common analysis set, CO2/N2 Henry selectivity emerges through a non-equivalent sequence of interaction stages:

1. Repulsive WCA accessibility establishes an N2-favoring baseline in every framework.
2. Restoring Lennard-Jones attraction reverses that baseline toward CO2 throughout the panel.
3. Electrostatics supplies a further positive but strongly framework-dependent enhancement.
4. Conventional pore geometry predicts the repulsive and attractive coordinates strongly and predicts the final endpoint well, yet provides essentially no improvement in absolute prediction error for the electrostatic increment relative to a median predictor.
5. M029/M085 provides the strongest local counterexample: the two structures are mutual nearest neighbors in the nine-descriptor geometry space but differ by 1.860 in electrostatic selectivity increment; the difference is overwhelmingly CO2-driven and is accompanied by a markedly different probability-supported spatial response.
6. Attraction and electrostatics alter screening order differently: attraction reorganizes the ranking globally, whereas electrostatics mostly preserves the global LJ-only order while creating important framework-specific outliers.

## Frozen quantitative anchors

- n = 238
- median ln S_WCA = -0.203749485237
- median S_WCA = 0.81566668214
- WCA-positive fraction = 0
- median Delta ln S_attr = 2.53737562365
- median attraction multiplier = 12.6464383687
- attraction-positive fraction = 1
- median Delta ln S_elec = 0.309492036602
- median electrostatic multiplier = 1.36273271995
- electrostatic-positive fraction = 1
- Spearman rho(Delta ln S_attr, Delta ln S_elec) = -0.0885675386621
- maximum telescoping residual = 8.881784197e-16

Geometry-only models:
- ln S_WCA: R2 = 0.835027; MAE = 0.032935; baseline MAE = 0.092870
- Delta ln S_attr: R2 = 0.874487; MAE = 0.177359; baseline MAE = 0.549456
- Delta ln S_elec: R2 = 0.188040; MAE = 0.146024; baseline MAE = 0.146993
- ln S_full: R2 = 0.845241; MAE = 0.188768; baseline MAE = 0.501407

Electrostatic matched pair:
- M029 ln S_WCA = -0.2013; Delta ln S_elec = 0.1407
- M085 ln S_WCA = -0.1982; Delta ln S_elec = 2.0007
- electrostatic contrast = 1.8599934028
- multiplicative contrast = 6.42369439287
- M029/M085 are mutual nearest neighbors in the nine-descriptor standardized geometry space
- geometry distance = 0.257529

Per-gas global electrostatic increments:
- M029 CO2: 0.162; N2: 0.021
- M085 CO2: 2.032; N2: 0.032
- between-framework CO2 contrast: 1.871 log units
- between-framework N2 contrast: 0.011 log units

Spatial summaries:
- M029 positive probability-weight fraction = 0.576
- M029 weighted median local delta_elec = 0.042
- M029 weighted IQR = 0.614
- M085 positive probability-weight fraction = 0.791
- M085 weighted median local delta_elec = 0.534
- M085 weighted IQR = 2.019

Population rank transitions:
- WCA -> LJ rank Spearman rho = -0.7508597634
- median |Delta rank| after attraction = 102
- max attraction promotion/demotion = 234 / 226
- LJ -> full rank Spearman rho = 0.9050550573
- median |Delta rank| after electrostatics = 9
- max electrostatic promotion/demotion = 175 / 55

Bounded WCA force-field bridge:
- CO2 MAE shift = 0.0227 kcal/mol; max = 0.0563; rho = 0.996
- N2 MAE shift = 0.0178 kcal/mol; max = 0.0507; rho = 1.000

## Publication boundaries

Do not claim:
- attraction and electrostatics are statistically independent;
- geometry contains no electrostatic information;
- the nested differences are unique experimental free-energy components;
- a specific chemical binding motif has been identified from the spatial maps;
- universal force-field independence;
- finite-pressure mixture performance;
- transport selectivity;
- framework-flexibility behavior;
- transfer to other gas pairs;
- universal behavior across all MOFs.

Use:
- "in this 238-MOF analysis set"
- "under the selected rigid-framework molecular models"
- "Henry-regime"
- "weakly associated"
- "little practical improvement in absolute electrostatic prediction error"
- "overwhelmingly CO2-driven" for the M029/M085 electrostatic contrast

## Main-figure architecture

1. Figure 1 — nested decomposition and population distributions.
2. Figure 2 — population-wide attraction/electrostatic relationship; matched pairs highlighted, no categorical mechanism regions.
3. Figure 3 — geometry information boundary.
4. Figure 4 — M029/M085: trajectory, gas-specific global increments, and local electrostatic maps.
5. Figure 5 — all-238 rank transitions.

Supporting Information:
- precision/exclusion robustness;
- geometry-model details;
- M089/M192 scalar attraction pair;
- per-gas M029/M085 spatial maps;
- population rank statistics;
- bounded WCA force-field bridge;
- source-data/reproducibility inventory.

## Status

The scientific analysis is closed for the current claim set.
No new simulation is required before PI additions.
Local LaTeX compilation is explicitly not a scientific gate; final document compilation will be performed in Overleaf.
