# Publication story — frozen cleanup version

## One-sentence thesis

CO2/N2 selectivity in this 238-MOF panel is built in stages: repulsion sets an N2-favoring baseline, attraction creates CO2 preference, and electrostatics supplies an additional framework-dependent recognition channel that conventional pore geometry largely fails to encode.

## Claim hierarchy

### Claim 1 — baseline
All 238 WCA references favor N2. Repulsion is therefore a baseline, not a successful CO2-selectivity route in this panel.

### Claim 2 — creation of selectivity
Restoring Lennard-Jones attraction makes every framework CO2-selective. This is the stage that creates CO2 preference.

### Claim 3 — distinct electrostatic information
Electrostatics increases selectivity in every framework, but the increment varies strongly and is weakly associated with the attractive increment.

### Claim 4 — representation boundary
Pore geometry predicts the WCA and attractive coordinates strongly, but its electrostatic MAE is essentially identical to a median-predictor baseline. A representation can therefore predict the final endpoint while missing an important mechanism that contributes to it.

### Claim 5 — microscopic counterexample
M029 and M085 have nearly identical WCA selectivities but electrostatic increments of 0.141 and 2.001. The difference is overwhelmingly CO2-driven, and probability-supported full-minus-LJ maps show a much stronger and more heterogeneous local response in M085.

### Claim 6 — screening consequence
Restoring attraction or electrostatics can reorder a material by most of the 238-candidate ranking. Endpoint rank and physical mechanism answer different questions.

### Bounded robustness claim
The WCA coordinate is stable to the tested DREIDING-linker/UFF-Co perturbation on a frozen 15-MOF panel. This is not a force-field-independence claim for attraction or electrostatics.

## Figure logic

1. Figure 1 — What is being decomposed and what does each population distribution look like?
2. Figure 2 — Are attraction and electrostatics the same energetic axis? No.
3. Figure 3 — Which coordinates are already encoded by pore geometry? WCA/attraction yes; electrostatics practically no.
4. Figure 4 — Can similar geometry hide different mechanisms? Yes; M029/M085 provides the strongest spatially resolved counterexample.
5. Figure 5 — Does preserving mechanism matter for screening? Yes; rankings change drastically.
6. SI — convergence/precision, per-gas spatial attribution, WCA force-field bridge, full source tables.

## What is intentionally not part of the paper

- circular permutation analysis;
- invalid automated chemistry/RAC hierarchy;
- family leave-one-out claim;
- M089/M192 LJ map as spatial evidence of attraction;
- quantitative full/LJ alternative-force-field robustness from the failed shared-output bridge;
- CO2/CH4 transfer;
- finite-pressure, transport, flexibility, or universal gas-pair claims;
- narrative about debugging failures.

These remain in project/Git provenance where useful, but they do not belong in publication prose.
