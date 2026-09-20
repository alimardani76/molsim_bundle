# Frozen publication claim hierarchy

This document is the claim spine for the final manuscript. Changes require an explicit evidence-based decision.

## Central claim

CO2/N2 Henry-regime selectivity in the frozen 238-MOF panel is not a single structural property. Repulsive accessibility establishes an N2-favoring baseline, attractive Lennard-Jones interactions create CO2 preference across the panel, and electrostatics provides a distinct framework-dependent enhancement that is poorly encoded by conventional pore geometry.

## Primary claims

1. Repulsive baseline
   - None of the 238 MOFs is CO2-selective under the WCA reference.
   - Median S_WCA = 0.8157.
   - Interpretation: repulsive accessibility alone is not the successful CO2-selectivity route in this population.

2. Attraction creates CO2 preference
   - Delta ln S_attr > 0 for every MOF.
   - Median Delta ln S_attr = 2.5374 (12.65x multiplier).
   - Attraction is therefore the stage that reverses the N2-favoring repulsive baseline.

3. Electrostatics differentiates frameworks
   - Delta ln S_elec > 0 for every MOF.
   - Median Delta ln S_elec = 0.3095 (1.36x multiplier).
   - Attraction and electrostatics are weakly associated (Spearman rho = -0.0886).
   - Do not call them statistically or causally independent.

4. Geometry explains the first two coordinates much better than electrostatics
   - Geometry-only R2: WCA 0.835; attraction 0.874; electrostatics 0.188.
   - Electrostatic MAE = 0.146 versus baseline MAE = 0.147.
   - Allowed conclusion: conventional pore geometry adds little practically useful predictive information for the electrostatic increment.

5. Matched frameworks expose different electrostatic mechanisms
   - M029/M085 have nearly identical WCA log-selectivities (-0.2013 vs -0.1982).
   - Delta ln S_elec = 0.1407 vs 2.0007; difference = 1.8600 (6.42x multiplicative contrast).
   - The contrast is overwhelmingly CO2-driven from the authoritative per-gas Henry coefficients.
   - Probability-supported full-minus-LJ fields provide spatial corroboration.
   - Do not infer a specific atom type or binding motif from the maps alone.

## Supporting claims

6. Interaction stages substantially reorder candidate priorities.
   - Maximum promotion after attraction: 234 ranks.
   - Maximum promotion after electrostatics: 175 ranks.

7. The WCA coordinate is robust to the tested DREIDING-linker/UFF-Co perturbation.
   - This is bounded evidence on a preselected 15-MOF bridge.
   - It is not evidence of universal force-field independence.

## Explicitly excluded claims / analyses

- No family leave-one-out claim: the retained result file was empty and the authoritative table lacks the required family columns.
- No circular permutation result.
- No invalid automated chemistry/RAC hierarchy.
- No quantitative full/LJ/electrostatic force-field robustness from the failed shared-output RASPA bridge.
- No finite-pressure mixture selectivity.
- No transport selectivity.
- No framework-flexibility claim.
- No CO2/CH4 transfer claim.
- No universal gas-pair or universal MOF law.
