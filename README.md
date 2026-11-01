# CO2/N2 Henry-selectivity mechanisms in metal–organic frameworks

**Publication state:** JCTC v11 reconciliation branch (draft; not a validated end-to-end reproducibility release).  
**Manuscript provenance:** user-supplied `v10_P01_Selectivity_Origin_Map_PCCP_LaTeX_Submission.zip` and corrected `JCTC_Selectivity_Revised_v11.zip`, both dated 2026-11-01. The submitted v10 contains six main and four SI figures, unlike the older five-figure repository manuscript.

## Key result
On the 238-framework common panel, the WCA reference favors N2 (median ln S = -0.203749485237). Restoring LJ attraction raises log-selectivity (median increment 2.53737562365) and electrostatics raises it further (median increment 0.309492036602). All values refer to simulated Henry-regime CO2/N2 selectivity, **not** finite-pressure mixture performance.

## Repository contents
- `data/authoritative/path_b_final_238_authoritative.csv`: 238-framework calculation-level result table.
- `data/source_data/`: archived final descriptor CV predictions, metrics, matched-pair and rank-transition sources, WCA robustness and bridge sources. File provenance should be checked against each v11 panel before publication.
- `data/spatial_fields/`: retained final M029/M085 fields for the older mechanism-first visualization, not a complete simulation trajectory archive.
- `analysis/final/`: legacy final-analysis scripts. **Some contain machine-specific absolute/relative paths and do not run directly in this checkout.**
- `reproducibility/production_protocol.json`: simulation configuration; `reproducibility/*force_field_reference.json`, `CO2.json`, `N2.json`, and `uff_framework_registry.csv`: retained force-field reference artifacts.
- `manuscript/`, `figures/`: older five-main-figure branch artifacts. **Do not cite as the 2026-11-01 JCTC v11 manuscript.**
- `scripts/verify_publication_data.py`: lightweight integrity/thermodynamic arithmetic test requiring pandas and numpy.
- `docs/JCTC_V11_ALIGNMENT.md`: precise version mapping and missing reproduction prerequisites.

## Run the verified available checks
```bash
git clone https://github.com/alimardani76/molsim_bundle.git
cd molsim_bundle
python -m pip install numpy pandas
python scripts/verify_publication_data.py
```
This checks the authoritative table, the decomposition identity and the headline medians. It **does not** rerun RASPA, reconstruct cross-validation folds, or redraw all ten JCTC panels.

## Simulation used
According to the archived `production_protocol.json`: RASPA3 3.0.29 (conda build `h82e5c02_0`; archived Windows x86-64 executable SHA256 `37ab1e47c0fddaf70a28e5c27c838cc9e4bee3a9f12a87505793f5059e32df91`). Rigid frameworks at 298 K; 12 Å shifted LJ cutoff; Lorentz–Berthelot mixing; no tail corrections; framework LJ from UFF and framework charges from the ARC-MOF source; charged full model uses Ewald, LJ-only charges disabled. Full Widom sampling: 50,000 cycles per seed; LJ-only: 10,000 cycles per seed; default seeds 58101–58103. WCA sampling was adaptive. The archived protocol status is `frozen_candidate_pending_cross_machine_regression`, **not a proven exact cross-machine reproduction**.

## Reproduction levels
1. **Recompute article-level summaries:** run `scripts/verify_publication_data.py`; consult the frozen source-data CSVs.
2. **Inspect existing plots:** `figures/` contains the older mechanism-first publication figures. Source scripts `analysis/final/` need path normalization and missing external inputs for a clean rerun.
3. **Redraw exact JCTC v11 figures:** **not yet verified**; the v11 figure PDFs are available in the separate manuscript bundle, but the authoring scripts for each of its six main/four SI plots have not been traced and tested against final sources.
4. **Regenerate raw simulations:** **not yet available end-to-end**: the original CIF set, exact input/job-builder scripts, WCA job registry, and full seed-level outputs must be linked with hashes, and RASPA rerun must be validated.

The repository is a documented evidence archive, not yet a complete reproducibility capsule. See `docs/JCTC_V11_ALIGNMENT.md` for blockers.
