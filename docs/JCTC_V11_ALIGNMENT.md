# JCTC v11 / repository alignment ledger

## What was actually inspected
- Public repository default branch: `mechanism-first`, SHA `82ad362725c673527e38fc81c27ab60673d0303c` (2026-09-20); this has **five** main figures and **two** SI figures.
- The repository also contains a `main` branch with some different files. It should not be assumed that `main` is the publication release.
- User-provided 2026-11-01 `v10_P01_Selectivity_Origin_Map_PCCP_LaTeX_Submission.zip`: six main figure PDFs and four SI figure PDFs.
- Corrected `JCTC_Selectivity_Revised_v11.zip`: same ten figure PDFs (v10 figure bytes preserved) and corrected ACS/JCTC LaTeX main/SI. This ZIP does **not** contain the original simulation outputs, original CIFs, or figure-generation scripts.

## Positive matches with archived numerical report
The following medians in `reproducibility/manuscript_numbers_report.txt` agree with the v11 reported rounded quantities: WCA lnS -0.204; attraction increment 2.537; electrostatic increment 0.309; matched M029/M085 electrostatic increment difference 1.860; geometry-only CV R2 for attraction ~0.874 and electrostatics ~0.188.

## Known mismatches / omissions
1. **Figure-set mismatch:** six main and four SI PDFs in v11 versus five and two on `mechanism-first`. Never present the older PDFs as exact v11 plots.
2. **Data/code status:** `analysis/final/generate_path_b_publication_figures.py` reads `path_b_analysis/paper_figures/figure_case_values.csv` and `figure_pair_values.csv`; these are not in this repository tree. Its source figures are not an exact v11 replacement.
3. **Descriptor data:** `data/source_data/descriptor_hierarchy_predictions.csv` and `descriptor_hierarchy_metrics.csv` exist, contrary to the earlier statement that no prediction data were available anywhere. Their actual model/split correspondence to v11 requires row-level audit.
4. **Simulation replay:** protocol gives RASPA3 3.0.29 and seed/cycle parameters but lacks the complete one-command runner, source structure bundle and full per-job outputs. The protocol itself says cross-machine regression is pending.
5. **WCA jobs:** manuscript count of 2052 versus 2034 implied by listed retained population seed counts leaves 18 jobs to account for. Do not quietly allocate them.
6. **Structural claims:** no chemistry/topology attribution for M029/M085 is allowed without primary CIF evidence.
7. **Release lock:** no manuscript version tag or SHA256 lock for the complete v11 numerical and figure provenance has been established.

## Required for an actual paper-ready release
- Locate original JCTC v11 panel-generation scripts; port to repository-relative paths, ensure deterministic figure regeneration, and compare outputs/values against all ten retained PDFs.
- Recompute all tables, intervals, rank movements, matched pairs and ML metrics from archived 238-row table / prediction rows, checking fold provenance and source-data hashes.
- Resolve the 18-job WCA discrepancy from primary records, and confirm the actual final simulation binaries, charge/guest pseudotypes, CIF inputs and seed outputs.
- Add portable simulation input generation, run, parse, summarize commands with a pinned software environment.
- Only after tests pass, assemble a *separate* small release containing final numerical products and plot source files, without exploratory intermediates; tag it with a manuscript-specific version.

**Status:** alignment work started; not a completed publication-grade reproduction claim.
