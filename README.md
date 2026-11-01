# CO₂/N₂ selectivity in MOFs — reproducibility guide

This repository contains the **final archived numerical evidence** for the 238-framework Henry-regime study, together with portable scripts to calculate paper-level summaries and redraw the figures. Start here to locate the simulation setup, the retained outputs, the numbers reported in the paper, and the figures.

**Reference version:** corrected JCTC manuscript **v11 (2026-11-01)**. The older `manuscript/` and `figures/` folders in this repository are **not** that v11 submission. The updated figure code generates *data-based redraws*, not exact copies of the original v11 artwork.

## Quick start — results and figures (working, no RASPA needed)

Python 3.11 is the tested version.

```bash
git clone https://github.com/alimardani76/molsim_bundle.git
cd molsim_bundle
python -m pip install -r requirements-figures.txt
python scripts/verify_publication_data.py
python scripts/rebuild_paper_results.py
python scripts/redraw_jctc_v11.py
```

The commands produce `generated_results/` (CSV) and `generated_figures/jctc_v11/` (ten PDFs and ten editable SVGs). These output folders are generated locally rather than tracked in Git.

## 1. Simulation — what was run, using which software?

The study used **RASPA3 3.0.29** (conda build `h82e5c02_0`), Widom test-particle insertion for **CO₂ and N₂**, with rigid MOFs at **298 K**. The reference models are nested: **WCA repulsive reference → LJ-only → full LJ + electrostatics**. The reference protocol specifies a 12 Å shifted LJ cutoff, Lorentz–Berthelot mixing, UFF framework LJ parameters, ARC-MOF source charges, and Ewald electrostatics in the full model. Standard runs used three seeds (58101–58103); full-mode 50,000 cycles per seed, LJ-only 10,000; WCA sampling was adaptive.

**Retained simulation settings and inputs:**
- [RASPA version, executable hash, seeds and protocol](reproducibility/production_protocol.json)
- [CO₂ geometry/model reference](reproducibility/CO2.json) and [N₂ geometry/model reference](reproducibility/N2.json)
- [Full force-field reference](reproducibility/full_force_field_reference.json), [LJ-only reference](reproducibility/lj_only_force_field_reference.json), [UFF registry](reproducibility/uff_framework_registry.csv)

**Current reproducibility boundary:** This repository does **not** contain the full original CIF set, complete per-seed RASPA input/output directories, the WCA job ledger, or a verified portable job-generation and execution script. Consequently, **a from-CIF simulation rerun cannot currently be performed from this checkout alone**. These archived JSONs document the configuration; they are not a complete `raspa3` executable workflow. The production protocol is labelled `frozen_candidate_pending_cross_machine_regression`. Do not claim full simulation replay until these files have been restored and tested.

## 2. Simulation results — where are the computed Henry coefficients?

**Primary retained result:** [`data/authoritative/path_b_final_238_authoritative.csv`](data/authoritative/path_b_final_238_authoritative.csv) — **238 rows**, one per included MOF, with the aggregated guest-specific log Henry coefficients, per-mode log selectivities, WCA and sampling uncertainties, and geometry descriptors. It is the final analysis input rather than the raw individual RASPA output files.

The relevant column families are `lnkh_mean_*_co2`, `lnkh_mean_*_n2`, `lnS_wca`, `lnS_lj_only`, `lnS_full`, `delta_lnS_attr`, and `delta_lnS_elec`.

For the specific M029/M085 spatial analysis, final retained spatial arrays are in [`data/spatial_fields/`](data/spatial_fields/). The WCA 15-MOF force-field test is in [`data/source_data/wca_mof_gas_comparison_30.csv`](data/source_data/wca_mof_gas_comparison_30.csv).

## 3. Paper numbers — how are reported results obtained?

Run:

```bash
python scripts/rebuild_paper_results.py
```

This reads the final 238-row table and the saved out-of-fold ML prediction file. It **does not rerun simulations or refit ML models**. It checks the log-selectivity decomposition and generates:

| Output | Contents |
| --- | --- |
| `generated_results/final_selectivity_238.csv` | Clean per-MOF selectivities, interaction increments, Henry values, uncertainty columns and stage-wise ranks |
| `generated_results/paper_numbers.csv` | Calculated population medians, enhancement factors, rank shifts, top-24 overlap, M029/M085 contrast, and recalculated ML (R²)/MAEs |

Archived original inputs to this calculation are [`data/source_data/descriptor_hierarchy_predictions.csv`](data/source_data/descriptor_hierarchy_predictions.csv) (1,904 out-of-fold predictions), [`descriptor_hierarchy_metrics.csv`](data/source_data/descriptor_hierarchy_metrics.csv), and [`reproducibility/manuscript_numbers.csv`](reproducibility/manuscript_numbers.csv), the original frozen broader numerical ledger (its source paths are historical Windows paths). Additional final source-data tables are under [`data/source_data/`](data/source_data/).

The exact identity checked is `lnS_full = lnS_wca + delta_lnS_attr + delta_lnS_elec`. The matched-pair **geometry** calculation uses **nine features**, whereas the descriptor-based ML control uses **14**.

## 4. Figure sources and plotting code

Run:

```bash
python scripts/redraw_jctc_v11.py
```

- **Code:** [`scripts/redraw_jctc_v11.py`](scripts/redraw_jctc_v11.py).
- **Numerical source:** [238-MOF table](data/authoritative/path_b_final_238_authoritative.csv) plus the supporting [figure source CSVs](data/source_data/).
- **Output:** `generated_figures/jctc_v11/` — six main figures and four SI figures, each as PDF and editable SVG.
- **Download prebuilt files:** [successful GitHub Actions figure build](https://github.com/alimardani76/molsim_bundle/actions/runs/41309112396) → *Artifacts* → `jctc-v11-statistical-redraws`. A new build can be launched from [Actions](https://github.com/alimardani76/molsim_bundle/actions/workflows/jctc-v11-plots.yml).

| Figure | Direct numerical inputs |
| --- | --- |
| Main 1, 2, 6 | `data/authoritative/path_b_final_238_authoritative.csv` |
| Main 3 | `data/source_data/fig3_geometry_information_source.csv` |
| Main 4 | `data/source_data/fig4_matched_pair_scalar_source.csv` and `fig4_per_gas_global_increment_source.csv` |
| Main 5 | `data/source_data/fig5_rank_transition_source.csv` |
| SI S1, S2, S3, S4 | `data/authoritative/path_b_final_238_authoritative.csv` |

**Figure fidelity note:** The generator is a *new implementation of the final statistical plots*. It is verified to execute and produce all ten files, but visual identity with the original JCTC v11 PDFs has **not** been established. The original corrected v11 LaTeX files and PDF figures are in the separately supplied `JCTC_Selectivity_Revised_v11.zip`, not the older repository `manuscript/` or `figures/` folders.

## Scope and outstanding items

The current repository supports **(2) final aggregated numerical results → (3) paper statistics → (4) data-based figure redraws**. It does not yet provide a complete **(1) raw simulation replay**, nor a pixel-faithful rebuild of the submitted v11 figure artwork. The 18-job WCA accounting difference noted in [the alignment ledger](docs/JCTC_V11_ALIGNMENT.md) remains unresolved and should be reconciled from the primary run inventory.

See [JCTC v11 figure documentation](docs/JCTC_V11_REDRAW.md) for plot provenance and limitations.
