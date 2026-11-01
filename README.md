# CO₂/N₂ selectivity in MOFs — reproducible analysis and JCTC v11 figure redraws

This repository contains the curated **final** numerical evidence for a
238-framework rigid-MOF Henry-regime selectivity study, plus Python code to
recompute headline checks and generate **ten JCTC v11 statistical redraws**.

**Version information**
- Reference manuscript: corrected **JCTC v11** (2026-11-01), supplied as
  \`JCTC_Selectivity_Revised_v11.zip\` with six main and four SI PDFs.
- Plotting implementation: **new data-equivalent redraws**, not the unarchived
  original v11 plotting code; figures are not byte-identical.
- RASPA production protocol: **RASPA3 3.0.29**, conda build
  \`h82e5c02_0\`; full binary SHA256 in
  \`reproducibility/production_protocol.json\`.
- Repo manuscript and original PDFs in \`manuscript/\` and \`figures/\`
  are from an **older revision**; they are *not* the JCTC v11 source.

## Rebuild the final figures

\`\`\`bash
git clone https://github.com/alimardani76/molsim_bundle.git
cd molsim_bundle
git checkout jctc-v11-figure-redraw
python -m pip install -r requirements-figures.txt
python scripts/verify_publication_data.py
python scripts/redraw_jctc_v11.py
\`\`\`

Ten editable **SVG** and publication **PDF** redraws appear in
\`generated_figures/jctc_v11/\`. Run a subset using \`--only main\` or
\`--only si\`. Alternatively download the PDFs and SVGs from the successful
[GitHub Actions build](https://github.com/alimardani76/molsim_bundle/actions/workflows/jctc-v11-plots.yml).

## Scientific outputs and sources

| Stage | Final evidence |
| --- | --- |
| 238-framework Henry coefficient / selectivity analysis | \`data/authoritative/path_b_final_238_authoritative.csv\` |
| Mechanism, geometry, matched pair, ranks, shortlist | \`data/source_data/fig*.csv\`, \`matched_pair_geometry_context.csv\` |
| Descriptor model predictions and metrics | \`data/source_data/descriptor_hierarchy_predictions.csv\`, \`descriptor_hierarchy_metrics.csv\` |
| WCA precision and force-field bridge | final WCA fields in 238-row table; \`data/source_data/wca_*csv\` |
| M029 / M085 spatial diagnostics | retained \`data/spatial_fields/\` arrays and \`spatial_*summary.csv\` |
| Source plot implementation | \`scripts/redraw_jctc_v11.py\` |
| Original analysis scripts | \`analysis/final/\` (some legacy paths still need repair) |

The key reported medians: WCA log S −0.203749485237, attractive increment
2.537375623648, electrostatic increment 0.309492036602. The verification
script checks all 238 records and their telescoping thermodynamic identity.

## Simulation settings (archived production protocol)

Widom insertion with RASPA3 3.0.29; rigid frameworks, 298 K, shifted LJ
cutoff 12 Å, Lorentz–Berthelot mixing, no tail corrections, UFF framework
LJ parameters and source ARC-MOF partial charges; Ewald in the full model,
charges disabled in LJ-only. Full: 50,000 cycles/seed; LJ-only:
10,000 cycles/seed; nominal seeds 58101–58103. WCA is an adaptive
repulsive-reference calculation. The recorded protocol is marked
\`frozen_candidate_pending_cross_machine_regression\`.

**Reproducibility limitation:** the analysis-data checks and figure redraw
were executed by GitHub Actions. The original CIF inputs, complete RASPA
seed-level jobs, WCA ledger, and portable from-structure job generator
are not included, so raw simulations are **not yet independently rerunnable**.
The v11 corrected manuscript's LaTeX and its exact ten original PDFs also
have **not** been copied into the GitHub repository. Read
\`docs/JCTC_V11_REDRAW.md\` and \`docs/JCTC_V11_ALIGNMENT.md\` before citing
this as an exact paper reproduction.
