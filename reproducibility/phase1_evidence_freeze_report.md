# Phase 1 Evidence Freeze Report

Project root: `C:\Users\KaraPardazesh\Desktop\bundle\v1\mof_selectivity_mechanisms`
Generated: `2026-09-18T12:41:21`

## Evidence files

- authoritative_238: **FOUND** — `C:\Users\KaraPardazesh\Desktop\bundle\v1\mof_selectivity_mechanisms\path_b_analysis\path_b_final_238_authoritative.csv`
- negative_claim_matrix: **FOUND** — `C:\Users\KaraPardazesh\Desktop\bundle\v1\mof_selectivity_mechanisms\path_b_analysis\negative_controls\negative_control_claim_matrix.csv`
- exclusion_sensitivity: **FOUND** — `C:\Users\KaraPardazesh\Desktop\bundle\v1\mof_selectivity_mechanisms\path_b_analysis\negative_controls\exclusion_sensitivity_results.csv`
- family_leave_one_out: **FOUND / CSV PARSE FAILED** — `C:\Users\KaraPardazesh\Desktop\bundle\v1\mof_selectivity_mechanisms\path_b_analysis\negative_controls\family_leave_one_out_results.csv`
- descriptor_metrics: **FOUND** — `C:\Users\KaraPardazesh\Desktop\bundle\v1\mof_selectivity_mechanisms\path_b_analysis\leakage_free_descriptors\descriptor_hierarchy_metrics.csv`
- descriptor_predictions: **FOUND** — `C:\Users\KaraPardazesh\Desktop\bundle\v1\mof_selectivity_mechanisms\path_b_analysis\leakage_free_descriptors\descriptor_hierarchy_predictions.csv`
- descriptor_manifest: **FOUND** — `C:\Users\KaraPardazesh\Desktop\bundle\v1\mof_selectivity_mechanisms\path_b_analysis\leakage_free_descriptors\descriptor_feature_manifest.csv`
- excluded_leakage_columns: **FOUND** — `C:\Users\KaraPardazesh\Desktop\bundle\v1\mof_selectivity_mechanisms\path_b_analysis\leakage_free_descriptors\excluded_leakage_columns.csv`
- wca_paired_seed: **FOUND** — `C:\Users\KaraPardazesh\Desktop\bundle\v1\mof_selectivity_mechanisms\path_b_analysis\forcefield_bridge\final_analysis\wca_paired_seed_comparison_90.csv`
- wca_mof_gas: **FOUND** — `C:\Users\KaraPardazesh\Desktop\bundle\v1\mof_selectivity_mechanisms\path_b_analysis\forcefield_bridge\final_analysis\wca_mof_gas_comparison_30.csv`
- wca_robustness: **FOUND** — `C:\Users\KaraPardazesh\Desktop\bundle\v1\mof_selectivity_mechanisms\path_b_analysis\forcefield_bridge\final_analysis\wca_robustness_statistics.csv`
- wca_largest_shifts: **FOUND** — `C:\Users\KaraPardazesh\Desktop\bundle\v1\mof_selectivity_mechanisms\path_b_analysis\forcefield_bridge\final_analysis\wca_largest_forcefield_shifts.csv`
- bridge_panel: **FOUND** — `C:\Users\KaraPardazesh\Desktop\bundle\v1\mof_selectivity_mechanisms\path_b_analysis\forcefield_bridge\forcefield_bridge_panel_compact.csv`
- bridge_coverage: **FOUND** — `C:\Users\KaraPardazesh\Desktop\bundle\v1\mof_selectivity_mechanisms\path_b_analysis\forcefield_bridge\forcefield_coverage_by_mof.csv`
- bridge_fallback: **FOUND** — `C:\Users\KaraPardazesh\Desktop\bundle\v1\mof_selectivity_mechanisms\path_b_analysis\forcefield_bridge\forcefield_fallback_manifest.csv`

## Summary

- Expected files: 15
- Found: 15
- Missing: 0
- CSV parse failures: 1
- Stale/prohibited TeX hits: 2

## CSV parse failures

### family_leave_one_out

`C:\Users\KaraPardazesh\Desktop\bundle\v1\mof_selectivity_mechanisms\path_b_analysis\negative_controls\family_leave_one_out_results.csv`

Error: `default: EmptyDataError: No columns to parse from file | utf8sig: EmptyDataError: No columns to parse from file | latin1: EmptyDataError: No columns to parse from file | python_auto_sep: ValueError: bad delimiter value`

Raw preview saved as `raw_preview__family_leave_one_out.txt`.


## Gate

**PASS WITH PARSE WARNINGS — all expected files exist and are hashed; one or more CSVs need format diagnosis before using their contents.**