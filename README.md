# MOF CO2/N2 selectivity mechanisms — publication workspace

This repository contains the curated publication state, not the full exploratory simulation history.

## Scientific story

Repulsion establishes an N2-favoring baseline across the 238-MOF population. Restoring Lennard-Jones attraction creates CO2 selectivity in every framework. Electrostatics then provides a further framework-dependent enhancement that is weakly associated with attraction and poorly encoded by conventional pore geometry. Geometry-matched M029/M085 provides the strongest microscopic counterexample, with a 1.860-unit electrostatic-selectivity contrast that is overwhelmingly CO2-driven and accompanied by a markedly different probability-supported spatial response.

See `docs/STORY.md` for the claim hierarchy and figure logic.

## Repository structure

- `manuscript/` — clean main text, SI, bibliography
- `figures/main/` — five main publication figures
- `figures/si/` — supporting figures
- `data/authoritative/` — frozen 238-MOF table
- `data/source_data/` — figure/statistical source data
- `data/spatial_fields/` — retained final M029/M085 spatial fields only
- `analysis/final/` — scripts used by the final evidence path
- `reproducibility/` — protocols, manifests, hashes
- `docs/` — claims, decisions, exclusion ledger, PI additions, changelog

## Scope

No new simulation is required by the current claim set. New analyses should be added only if they resolve a concrete PI/reviewer request that cannot be answered from the frozen evidence.
