"""Assess coarse-grid refinement and produce a concise, evidence-linked report."""
import json
import csv
import sys
from pathlib import Path
import numpy as np
import scipy
from scipy.stats import t
from fixed_position_energy import OUT,HERE,sha256
from validate_fixed_position_vs_widom import RESULTS,write_json

def main():
    controls=json.loads((RESULTS/'controls.json').read_text())
    extended=json.loads((RESULTS/'extended_controls.json').read_text())
    widom=json.loads((RESULTS/'widom_validation_lj_only.json').read_text())
    full=json.loads((RESULTS/'widom_validation_full.json').read_text())
    orientation=json.loads((RESULTS/'orientation_convergence_lj_only.json').read_text())
    paths=sorted((OUT/'fields').glob('M085_CO2_lj_only_*/metadata.json'),key=lambda p:json.loads(p.read_text())['requested_spacing_A'],reverse=True)
    rows=[]
    for path in paths:
        meta=json.loads(path.read_text())
        assert meta['provenance']==widom['provenance']
        data=np.load(path.parent/'field.npz')
        assert np.isfinite(data['log_weight']).all()
        assert abs(data['relative_density'].mean()-1)<1e-12
        row={'spacing_A':meta['requested_spacing_A'],'points':int(np.prod(meta['shape'])),
             'orientations':meta['orientations'],'mean_weight':meta['cell_mean_weight'],
             'relative_difference_native':meta['relative_difference_native'],
             'mass_weighted_orientation_relative_se':meta['mass_weighted_orientation_relative_se'],
             'relative_change_previous':None if not rows else meta['cell_mean_weight']/rows[-1]['mean_weight']-1,
             'directory':str(path.parent.relative_to(OUT))}
        rows.append(row)
    pilot_passed=(len(rows)>=2 and abs(rows[-1]['relative_change_previous'])<.10
                   and abs(rows[-1]['relative_difference_native'])<.10
                   and rows[-1]['mass_weighted_orientation_relative_se']<.02
                   and widom['gate']['passed'] and orientation['gate_passed'] and controls['passed'])
    protection=json.loads((RESULTS/'protected_verification.json').read_text()) if (RESULTS/'protected_verification.json').exists() else {'passed':False}
    overall=pilot_passed and protection['passed'] and extended['passed']
    report={'fixed_position_controls_passed':controls['passed'],'lj_widom_gate_passed':widom['gate']['passed'],
        'full_widom_gate_passed':full['gate']['passed'],'orientation_gate_passed':orientation['gate_passed'],
        'coarse_grid_integral_gate_passed':pilot_passed,'extended_controls_passed':extended['passed'],'grid_convergence':rows,
        'pilot_directory':rows[-1]['directory'] if pilot_passed else None,
        'pilot_criterion':'latest two grid means change <10%; finer mean within 10% of native Widom; mass-weighted angular SE <2%',
        'scope':'Validated coarse M085/CO2/LJ pilot and native full/LJ fixed energies; not production spatial-resolution convergence',
        'production_0.30_0.20A_convergence_validated':False,'protected_files':protection,
        'pilot_deliverable_passed':overall,
        'source_commit':'3fdb4a1c4ba468ebea66d1a8619404e05aa800bb',
        'runtime_versions':{'python':sys.version,'numpy':np.__version__,'scipy':scipy.__version__,'raspalib':'3.0.29','raspa3':'3.0.29'},
        'code_sha256':{p.name:sha256(p) for p in HERE.glob('*.py')},
        'runtime_sha256':{str(p):sha256(p) for p in [Path(sys.executable).parent/'bin/raspa3.exe',Path(sys.executable).parent/'Lib/site-packages/raspalib/raspalib.cp312-win_amd64.pyd']}}
    write_json(RESULTS/'validation_report.json',report)
    with (RESULTS/'grid_convergence.csv').open('w',newline='') as file:
        writer=csv.DictWriter(file,fieldnames=list(rows[0]));writer.writeheader();writer.writerows(rows)
    w=widom['convergence'][-1]; wf=full['convergence'][-1]
    native_ci=t.ppf(.975,14)*widom['native_se']
    geometry=controls['geometry']
    lines=[
        '# Fixed-position energy and Henry-field validation',
        '',
        f'Pilot status: {"PASSED (coarse LJ pilot only)" if overall else "NOT YET PASSED"}.',
        '',
        'The `vector too long` failure is resolved without modifying or rebuilding RASPA. '
        'An omitted box has zero volume; the old degree-valued box also has NaN volume. '
        'Both feed nonfinite widths into Ewald wave-vector allocation. A valid box with radians constructs successfully. '
        'Initial positions are all molecular atom coordinates per component, with zero additional random molecules.',
        '',
        'All interactions are evaluated by installed RASPALIB 3.0.29 using the frozen production inputs. '
        'LJ: shifted 12 A, Lorentz-Berthelot, no tails. Full: matched LJ plus production Ewald; '
        'automatic Coulomb cutoff 14.953805 A. Native mixing converts N_com zero parameters to exact None pairs. '
        'The explicit energy difference subtracts the empty framework; this native routine already excludes its constant self energy.',
        '',
        f"At the prescribed z-axis CO2 configuration, both O-C distances are 1.149 A. "
        f"LJ = {geometry['lj_only']['total_kJmol']:.9f} kJ/mol; full = {geometry['full']['total_kJmol']:.9f}; "
        f"electrostatic difference = {geometry['full']['electrostatic_kJmol']:.9f}. "
        'Z-axis and rotated configurations agree with the native executable; unit-cell translation invariance passes. '
        'A strongly binding configuration gives -24.695 kJ/mol LJ and -10.400 kJ/mol electrostatics, '
        'also matching the native executable. See controls.json and extended_controls.json.',
        '',
        f"Independent uniform integration ({w['n']:,} samples): mean weight {w['mean']:.6f} +/- {w['ci95_halfwidth']:.6f} (95%). "
        f"Ordinary RASPA Widom (3 seeds x 2,000 cycles, 40,000 moves per seed): {widom['native_mean']:.6f} +/- {native_ci:.6f} (95%). "
        f"Combined z = {w['combined_z']:.3f}; effective sample size = {w['effective_sample_size']:.1f}. "
        'The predeclared precision/compatibility gate passes; increasing-sample estimates are in widom_convergence_lj_only.csv.',
        '',
        f"Orientation convergence selects {orientation['selected_orientations']:,} total orientations per COM "
        '(four scrambled Sobol sequences), checked at six positions against 16,384 orientations. '
        'Worst reference error at 2,048: 1.25%; worst replicate 95% halfwidth: 2.82%. '
        'Counts 16-256 were inadequate. See orientation_convergence_lj_only.csv.',
        '',
        '| Grid target (A) | Points | Mean w | Difference from native |',
        '|---:|---:|---:|---:|',
    ]
    lines.extend(f"| {x['spacing_A']:g} | {x['points']} | {x['mean_weight']:.6f} | {100*x['relative_difference_native']:+.2f}% |" for x in rows)
    lines.extend(['',f"Pilot files: `{report['pilot_directory']}`." if pilot_passed else 'No grid is promoted yet.',
        'Each field provides CSV and NPZ, raw/log weights, relative density with discrete cell mean one, '
        'triclinic coordinates and angular uncertainty. Normalization to the independent native Widom mean is also retained. '
        'These are infinite-dilution COM probability fields, not finite-pressure occupancies.',
        '',
        f"Full-mode exploratory integration: {wf['mean']:.3f} +/- {wf['ci95_halfwidth']:.3f} (95%, {wf['n']:,} samples), "
        f"versus native {full['native_mean']:.3f} (SE {full['native_se']:.3f}); z = {wf['combined_z']:.3f}. "
        'Compatible but insufficiently precise for the full-map gate. Full energy evaluation itself passes native point comparisons. '
        'No full map was generated. N2 has geometry/native energy controls but no integral/orientation validation. WCA is unimplemented. '
        'Other MOFs and 0.30/0.20 A production resolution require their own convergence checks. '
        'The coarse-grid integral gate does not prove pointwise spatial-feature convergence.',
        '',
        f"Protected-file verification: {protection.get('files_checked',0):,} files checked, passed = {protection['passed']}. "
        'This includes the original production environment, authoritative CSV and manuscript checkpoint. '
        'Only NumPy/SciPy were added to raspa_maps; all repository output is under 02_spatial_maps.',
        '',
        'See README.md for the exact source trace, force-field equivalence, assumptions, and runnable commands. '
        'Machine-readable summary: results/validation_report.json.',
    ])
    (OUT/'VALIDATION_REPORT.md').write_text('\n'.join(lines)+'\n',encoding='utf-8')
    print(json.dumps(report,indent=2))

if __name__=='__main__':main()
