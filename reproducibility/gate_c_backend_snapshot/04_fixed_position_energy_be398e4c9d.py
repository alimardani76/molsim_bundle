"""RASPA 3.0.29 native fixed-site energies. Run with raspa_maps/python.exe.

No Python interaction potential is used. Only input translation, geometry and
unit conversion occur here. Supported input subset is deliberately strict.
"""
from __future__ import annotations
import hashlib
import json
import math
import shlex
from pathlib import Path
import numpy as np
import raspalib as r

HERE = Path(__file__).resolve().parent
OUT = HERE.parents[1]
ROOT = OUT.parents[2]
TEMPLATES = ROOT / "distribution/ARCMOF_MECH_V1_PC2/templates"
SHARED = ROOT / "distribution/ARCMOF_MECH_V1_PC2/shared"
# Exact constants from RASPA 3.0.29 units.ixx / units.cpp (not modern replacements).
AMU = 1.6605402e-27
NA = 6.0221419947e23
KB_SI = 1.380650324e-23
R = 8.314464919
ENERGY_SI = AMU * 1e-10**2 / 1e-12**2
KB = KB_SI / ENERGY_SI
TO_KJMOL = ENERGY_SI * NA / 1000
TO_K = ENERGY_SI * NA / R

def sha256(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()

def read_p1(path):
    """Read explicit P1 CIF; fail closed for symmetry/occupancy we don't handle."""
    lines = Path(path).read_text().splitlines()
    scalars, loops = {}, []
    i = 0
    while i < len(lines):
        s = lines[i].strip()
        i += 1
        if not s or s.startswith('#'):
            continue
        if s == 'loop_':
            heads, rows = [], []
            while i < len(lines) and lines[i].strip().startswith('_'):
                heads.append(lines[i].strip()); i += 1
            while i < len(lines):
                s = lines[i].strip()
                if s.startswith('_') or s == 'loop_' or s.startswith('data_'):
                    break
                i += 1
                if s and not s.startswith('#'):
                    row = shlex.split(s)
                    if len(row) != len(heads):
                        raise ValueError('Unsupported multiline CIF row')
                    rows.append(dict(zip(heads, row)))
            loops.extend(rows)
        elif s.startswith('_'):
            k, v = shlex.split(s); scalars[k] = v
    assert scalars['_symmetry_Int_Tables_number'] == '1', 'Only explicit P1 supported'
    for row in loops:
        if '_symmetry_equiv_pos_as_xyz' in row:
            assert row['_symmetry_equiv_pos_as_xyz'].replace(' ', '') == 'x,y,z'
        assert float(row.get('_atom_site_occupancy', 1)) == 1
    lengths = [float(scalars['_cell_length_' + a]) for a in 'abc']
    angles = np.radians([float(scalars['_cell_angle_' + a]) for a in ['alpha','beta','gamma']])
    a,b,c = lengths; alpha,beta,gamma = angles
    temp = (math.cos(alpha)-math.cos(gamma)*math.cos(beta))/math.sin(gamma)
    cell = np.array([[a,0,0], [b*math.cos(gamma),b*math.sin(gamma),0],
                     [c*math.cos(beta),c*temp,c*math.sqrt(1-math.cos(beta)**2-temp**2)]])
    rows = [x for x in loops if '_atom_site_fract_x' in x]
    return lengths, angles, cell, rows

class FixedPositionEnergy:
    def __init__(self, mof='M085', gas='CO2', mode='full'):
        import sys
        records = list((Path(sys.executable).parent/'conda-meta').glob('raspalib-*.json'))
        if len(records)!=1 or json.loads(records[0].read_text())['version']!='3.0.29':
            raise RuntimeError('This evaluator is validated only with conda RASPALIB 3.0.29')
        if gas not in ('CO2','N2') or mode not in ('full','lj_only'):
            raise ValueError('Supported: CO2/N2, full/lj_only. WCA is not implemented.')
        self.mof,self.gas,self.mode = mof,gas,mode
        template = TEMPLATES / mof / mode
        self.ff_path = template / 'force_field.json'
        self.sim_path = template / 'simulation_base.json'
        self.gas_path = SHARED / (gas+'.json')
        cifs = list(template.glob('*.cif')); assert len(cifs) == 1
        self.cif_path = cifs[0]
        raw = json.loads(self.ff_path.read_text())
        sim = json.loads(self.sim_path.read_text())['Systems'][0]
        guest = json.loads(self.gas_path.read_text())
        assert set(raw) <= {'MixingRule','TruncationMethod','TailCorrections','CutOffVDW','CutOffCoulomb','PseudoAtoms','SelfInteractions'}
        assert raw['MixingRule'] == 'Lorentz-Berthelot'
        assert raw['TruncationMethod'] == 'shifted' and raw['TailCorrections'] is False
        assert sim['UseChargesFrom'] == 'PseudoAtoms'
        assert sim['ChargeMethod'] == ('Ewald' if mode == 'full' else 'None')
        assert guest['Type'] == 'rigid'
        self.temperature = sim['ExternalTemperature']
        self.replication = np.array(sim['NumberOfUnitCells'])
        lengths,angles,self.unit_cell,rows = read_p1(self.cif_path)
        self.cell = self.unit_cell * self.replication[:,None]
        self.volume = abs(np.linalg.det(self.cell))
        widths = self.volume / np.array([np.linalg.norm(np.cross(self.cell[1],self.cell[2])),
                    np.linalg.norm(np.cross(self.cell[0],self.cell[2])),np.linalg.norm(np.cross(self.cell[0],self.cell[1]))])
        self.cutoff_vdw = raw['CutOffVDW']
        assert widths.min() >= 2*self.cutoff_vdw
        # JSON defaults to automatic Coulomb cutoff. Compute BEFORE constructing
        # System so its cached framework structure factors use the final Ewald alpha.
        self.cutoff_coulomb = raw.get('CutOffCoulomb', 0.5*float(widths.min())-np.finfo(float).eps)
        pseudo = raw['PseudoAtoms']; self.names = {p['name']:i for i,p in enumerate(pseudo)}
        self.pseudo = pseudo
        elements = {'H':1,'C':6,'N':7,'O':8,'Co':27,'Zn':30,'Cu':29,'Zr':40}
        pa = [r.PseudoAtom(p['name'],p['framework'],p['mass'],p['charge'],0.,elements[p['element']],p['print_to_output'],p.get('source','')) for p in pseudo]
        interactions = {p['name']:p for p in raw['SelfInteractions']}
        params = []
        for p in pseudo:
            v = interactions[p['name']]
            if v['type'].lower() == 'none':
                assert p['name'] == 'N_com'
                params.append(r.VDWParameters(0.,0.))
            else:
                assert v['type'].lower() == 'lennard-jones'
                params.append(r.VDWParameters(*v['parameters']))
        if mode == 'lj_only':
            assert all(p['charge'] == 0 for p in pseudo)
        self.ff = r.ForceField(pa,params,cutoff_framework_vdw=self.cutoff_vdw,
            cutoff_molecule_vdw=self.cutoff_vdw,cutoff_coulomb=self.cutoff_coulomb,
            shifted=True,tail_corrections=False,use_charge=mode=='full')
        def atom(name,position,framework=False):
            idx = self.names[name]
            return r.Atom(list(map(float,position)),pseudo[idx]['charge'],1.,0,idx,0,0,int(framework))
        fw_atoms = [atom(x['_atom_site_label'],[float(x['_atom_site_fract_'+a]) for a in 'xyz'],True) for x in rows]
        self.framework = r.Framework(self.ff,self.cif_path.stem,r.SimulationBox(*lengths,*map(float,angles)),1,fw_atoms,self.replication.tolist())
        guest_atoms = [atom(name,pos) for name,pos in guest['pseudoAtoms']]
        self.relative = np.array([pos for name,pos in guest['pseudoAtoms']],dtype=float)
        self.masses = np.array([pseudo[self.names[name]]['mass'] for name,pos in guest['pseudoAtoms']])
        assert np.max(np.abs(np.average(self.relative,axis=0,weights=self.masses))) < 1e-12
        self.component = r.Component(self.ff,gas,guest['CriticalTemperature'],guest['CriticalPressure'],guest['AcentricFactor'],
            defined_atoms=guest_atoms,particle_probabilities=r.MCMoveProbabilities(widom_probability=1.))
        common = dict(force_field=self.ff,framework_components=self.framework,components=[self.component],
                      initial_number_of_molecules=[0],external_temperature=self.temperature,external_pressure=sim['ExternalPressure'])
        self.empty = r.System(**common)
        self.baseline = self._parts(self.empty.compute_total_energies())
        initial = self.relative + np.array([5.,5.,5.])
        self.system = r.System(**common,initial_positions=[initial.tolist()])
        assert self.system.number_of_molecules_per_component == [1]
        assert np.isclose(self.system.simulation_box.volume,self.volume,rtol=1e-12)
        self.atoms = self.system.atom_data
        assert np.max(np.abs(np.array([a.position for a in self.atoms[-len(initial):]])-initial)) < 1e-12
        self.guest_atoms = self.atoms[-len(initial):]
        # In this installed pybind build the list is new, but Atom wrappers refer
        # to native vector elements. Prove write-through before relying on it.
        probe=initial[-1]+[.125,0,0]
        self.guest_atoms[-1].position=probe.tolist()
        assert np.allclose(self.system.atom_data[-1].position,probe,atol=1e-13,rtol=0)
        self.guest_atoms[-1].position=initial[-1].tolist()
        self.nframework = len(self.atoms)-len(self.relative)
        assert self.nframework == len(rows)*int(np.prod(self.replication))
        self.provenance = {'mof':mof,'gas':gas,'mode':mode,'temperature_K':self.temperature,
            'replication':self.replication.tolist(),'cell_rows_A':self.cell.tolist(),'volume_A3':self.volume,
            'cutoff_vdw_A':self.cutoff_vdw,'cutoff_coulomb_A':self.cutoff_coulomb,
            'framework_atoms':self.nframework,'files':{str(p.relative_to(ROOT)):sha256(p) for p in [self.ff_path,self.sim_path,self.gas_path,self.cif_path]}}

    @staticmethod
    def _parts(e):
        return np.array([e.potential_energy(),e.van_der_waals_energy(),e.coulomb_energy()])

    def positions(self, com, orientation=(0.,0.,1.)):
        com = np.asarray(com,dtype=float); axis = np.asarray(orientation,dtype=float)
        if com.shape != (3,) or axis.shape != (3,) or not np.isfinite(com).all() or not np.isfinite(axis).all():
            raise ValueError('COM and linear-molecule axis must be finite length-three vectors')
        if np.linalg.norm(axis) == 0: raise ValueError('Zero molecular axis')
        axis = axis/np.linalg.norm(axis)
        assert np.max(np.abs(self.relative[:,:2])) == 0
        return com + self.relative[:,2,None]*axis

    def energy(self, com, orientation=(0.,0.,1.), verify=False):
        """Return E(host+guest)-E(host); axis is a unit-vector direction, not Euler angles.

        This system is for energies only: atom_data updates do not update the
        separate rigid-body dynamics state. Never use it to run MC/MD moves.
        """
        pos = self.positions(com,orientation)
        for a,p in zip(self.guest_atoms,pos): a.position = p.tolist()
        values = self._parts(self.system.compute_total_energies()) - self.baseline
        if not np.isfinite(values).all(): raise FloatingPointError('Nonfinite native energy')
        result = dict(zip(['total_internal','vdw_internal','electrostatic_internal'],values.tolist()))
        result.update(dict(zip(['total_kJmol','vdw_kJmol','electrostatic_kJmol'],(values*TO_KJMOL).tolist())))
        result['log_weight'] = -values[0]/(KB*self.temperature)
        if verify:
            actual = np.array([a.position for a in self.system.atom_data[-len(pos):]])
            assert np.max(np.abs(actual-pos)) < 1e-12
            assert np.max(np.abs(np.average(actual,axis=0,weights=self.masses)-com)) < 1e-12
            result['actual_positions_A'] = actual.tolist()
            result['COM_A'] = np.average(actual,axis=0,weights=self.masses).tolist()
            result['baseline_internal'] = self.baseline.tolist()
        return result

    def log_weights(self, com, axes):
        """Batch interface, native energy for every orientation; no interpolation.

        Do not resize or replace system.atom_data while this evaluator exists:
        the retained guest Atom wrappers reference its native elements.
        """
        com=np.asarray(com,dtype=float);axes=np.asarray(axes,dtype=float)
        assert com.shape==(3,) and axes.ndim==2 and axes.shape[1]==3
        assert np.isfinite(com).all() and np.isfinite(axes).all()
        norms=np.linalg.norm(axes,axis=1)
        assert np.all(norms>0)
        axes=axes/norms[:,None]
        positions=com+axes[:,None,:]*self.relative[None,:,2,None]
        result=np.empty(len(axes))
        beta=1/(KB*self.temperature)
        for i,pos in enumerate(positions):
            for atom,p in zip(self.guest_atoms,pos):atom.position=p.tolist()
            result[i]=-beta*(self.system.compute_total_energies().potential_energy()-self.baseline[0])
        if not np.isfinite(result).all():raise FloatingPointError('Nonfinite native log weight')
        return result

def uniform_axes(rng,n):
    z = rng.uniform(-1,1,n); phi = rng.uniform(0,2*np.pi,n)
    s = np.sqrt(1-z*z)
    return np.column_stack([s*np.cos(phi),s*np.sin(phi),z])

if __name__ == '__main__':
    import time
    start = time.perf_counter()
    for mode in ['lj_only','full']:
        e = FixedPositionEnergy(mode=mode)
        print(json.dumps(e.provenance,indent=2))
        print(json.dumps(e.energy([9.833953244182423,12.196323102710696,4.845113206703169],verify=True),indent=2))
    print('elapsed',time.perf_counter()-start)
