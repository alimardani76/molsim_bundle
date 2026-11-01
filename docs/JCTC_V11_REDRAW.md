# JCTC v11 plot reproduction: exact provenance and status

This branch contains a NEW portable Python redraw script for ten statistical
plots, sourced from the frozen final 238-framework table and selected final
source-data CSVs already committed in the repository.

## Versions
- Revised manuscript: **JCTC v11**, 2026-11-01.
- Original revised package: \`JCTC_Selectivity_Revised_v11.zip\` (user-provided);
  ten PDFs: figs 1–6, S1–S4, corrected main.tex/SI.
- Original repository paper: older \`mechanism-first\` main.tex and 5 main/2 SI plots.
- Simulation: **RASPA3 3.0.29**; exact build/hash and model details in
  \`reproducibility/production_protocol.json\`; protocol status remains
  \`frozen_candidate_pending_cross_machine_regression\`.

## Regenerate the ten statistical figure redraws

\`\`\`bash
git clone https://github.com/alimardani76/molsim_bundle.git
cd molsim_bundle
git checkout jctc-v11-figure-redraw
python -m pip install -r requirements-figures.txt
python scripts/verify_publication_data.py
python scripts/redraw_jctc_v11.py
\`\`\`

Output: \`generated_figures/jctc_v11/\` with **10 SVG and 10 PDF files**.
SVGs contain editable vectors and text. Output files are deliberately untracked.
GitHub Actions builds/upload the redraws as artifacts.

## Scope and integrity

- This script redraws reported **statistical content**, using *final*
  selectivity, geometry-model metrics, matched contrast, rank transitions,
  WCA SEs and per-gas Henry increment fields.
- The full 238-row file and final summary files in \`data/source_data/\` are
  inputs. No exploratory runs or raw RASPA trajectory dumps are necessary
  for redrawing these statistical panels.
- **Important:** the script is a new implementation, **not** the original
  unarchived v11 figure-source code. Layout, jitter, fonts, choice of
  illustrative cases and several panel styles can differ.
- The provided JCTC v11 *reference PDFs*, and corrected \`main.tex\` and SI,
  are **not yet present in this GitHub branch**. Do not misrepresent the
  older \`manuscript/\` or \`figures/\` directories as v11 artifacts.
- Exact rerun from CIF/seed/force-field input files remains pending.
- Count discrepancy in WCA jobs (2052 vs 2034) remains unresolved in primary records.

For manuscript evidence comparison and limitations, see
\`docs/JCTC_V11_ALIGNMENT.md\`.
