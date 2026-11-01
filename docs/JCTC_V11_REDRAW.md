# JCTC v11 figure reproduction

**Source of truth:** the corrected JCTC v11 package dated 2026-11-01 is separate from the older repository manuscript. This repository provides a tested *new statistical redraw*, not the original plotting implementation.

## One-command steps

```bash
python -m pip install -r requirements-figures.txt
python scripts/verify_publication_data.py
python scripts/rebuild_paper_results.py
python scripts/redraw_jctc_v11.py
```

Source data: `data/authoritative/path_b_final_238_authoritative.csv` plus final CSV files in `data/source_data/`.

Created files:
- `generated_results/final_selectivity_238.csv`
- `generated_results/paper_numbers.csv`
- `generated_figures/jctc_v11/*.pdf` (ten PDFs)
- `generated_figures/jctc_v11/*.svg` (ten editable SVGs)

A GitHub Actions workflow regenerates and uploads these files. This output documents **numerical regeneration and visualization**, not a rerun of RASPA or exact original figure artistry.

## Simulation provenance

RASPA3 3.0.29 (conda build `h82e5c02_0`) is recorded in `reproducibility/production_protocol.json`. Original CIF/seed-level run directories and the full portable job runner are not archived here. This is explicitly a **limitation of end-to-end reproducibility**.

## Source-data mapping

Main 1/2/6 and SI S1–S4 derive directly from the final 238-MOF table. Main 3 uses `fig3_geometry_information_source.csv`, main 4 uses `fig4_matched_pair_scalar_source.csv` and `fig4_per_gas_global_increment_source.csv`, and main 5 uses `fig5_rank_transition_source.csv`.

**Pending for exact JCTC v11 release:** original v11 figures/TeX in GitHub, exact panel styling verification, CIF + original job inputs + raw outputs, WCA job-count reconciliation.
