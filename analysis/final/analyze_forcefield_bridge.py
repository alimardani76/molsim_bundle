#!/usr/bin/env python3
from pathlib import Path
import json,re,math
import numpy as np
import pandas as pd
from scipy.stats import spearmanr, pearsonr
import matplotlib.pyplot as plt

R=Path.cwd(); O=R/'path_b_analysis/forcefield_bridge/final_analysis'; O.mkdir(parents=True,exist_ok=True)
OLD=R/'distribution/ARCMOF_WCA_LAMMPS_V1'; NEW=R/'distribution/ARCMOF_DREIDING_UFF_WCA_BRIDGE_V1'
pat=re.compile(r'FINAL\s+job=(\S+)\s+mof=(\S+)\s+gas=(\S+)\s+seed=(\d+)\s+insertions=(\d+)\s+mu_ex_kcal_mol=([-+0-9.eE]+)\s+mean_dU_kcal_mol=([-+0-9.eE]+)\s+insertion_volume_A3=([-+0-9.eE]+)')
def collect(root,label):
 rows=[]
 # Original UFF package preserves one result file per job, while log.lammps is
 # overwritten within each MOF directory. New isolated jobs each have one log.
 files=sorted(root.rglob('*.result.txt'))
 if not files: files=sorted(root.rglob('log.lammps'))
 seen=set()
 for fp in files:
  ms=pat.findall(fp.read_text(errors='replace'))
  if len(ms)!=1: continue
  job,mof,gas,seed,n,mu,du,vol=ms[0]
  key=(mof,gas,int(seed))
  if key in seen: raise SystemExit(f'Duplicate WCA key {key} from {fp}')
  seen.add(key)
  rows.append({'ff':label,'mof_id':mof,'gas':gas,'seed':int(seed),'insertions':int(n),'mu_ex_kcal_mol':float(mu),'mean_dU_kcal_mol':float(du),'volume_A3':float(vol),'source':str(fp.relative_to(root))})
 return pd.DataFrame(rows)
a_all=collect(OLD,'UFF'); b=collect(NEW,'DREIDING_UFF')
# The original UFF campaign contains 240 MOF-gas-seed rows. Restrict it to the
# exact 15-MOF bridge panel represented by the new 90-row campaign.
bridge_mofs=set(b['mof_id'].astype(str))
a=a_all[a_all['mof_id'].astype(str).isin(bridge_mofs)].copy()
keys=['mof_id','gas','seed']
expected=set(map(tuple,b[keys].astype({'mof_id':str,'gas':str,'seed':int}).to_records(index=False)))
observed=set(map(tuple,a[keys].astype({'mof_id':str,'gas':str,'seed':int}).to_records(index=False)))
missing=expected-observed; extra=observed-expected
if len(b)!=90 or len(a)!=90 or missing or extra:
 raise SystemExit(f'Bridge key mismatch: old_all={len(a_all)}, old_panel={len(a)}, new={len(b)}, missing={len(missing)}, extra={len(extra)}')
p=a.merge(b,on=keys,suffixes=('_uff','_dreiding_uff'),validate='one_to_one')
p['delta_mu_ex_kcal_mol']=p.mu_ex_kcal_mol_dreiding_uff-p.mu_ex_kcal_mol_uff
p['abs_delta_mu_ex_kcal_mol']=p.delta_mu_ex_kcal_mol.abs()
p['delta_mean_dU_kcal_mol']=p.mean_dU_kcal_mol_dreiding_uff-p.mean_dU_kcal_mol_uff
p.to_csv(O/'wca_paired_seed_comparison_90.csv',index=False)

g=p.groupby(keys[:2]).agg(mu_uff=('mu_ex_kcal_mol_uff','mean'),mu_dreiding_uff=('mu_ex_kcal_mol_dreiding_uff','mean'),delta_mu=('delta_mu_ex_kcal_mol','mean'),sd_delta_mu=('delta_mu_ex_kcal_mol','std'),max_abs_seed_delta=('abs_delta_mu_ex_kcal_mol','max'),n_seeds=('seed','nunique')).reset_index()
g['rank_uff']=g.groupby('gas').mu_uff.rank(method='average');g['rank_dreiding_uff']=g.groupby('gas').mu_dreiding_uff.rank(method='average');g['rank_shift']=g.rank_dreiding_uff-g.rank_uff
g.to_csv(O/'wca_mof_gas_comparison_30.csv',index=False)

stats=[]
for gas,x in g.groupby('gas'):
 rho,rhop=spearmanr(x.mu_uff,x.mu_dreiding_uff); r,rp=pearsonr(x.mu_uff,x.mu_dreiding_uff)
 stats.append({'scope':gas,'n':len(x),'mean_delta_mu':x.delta_mu.mean(),'median_delta_mu':x.delta_mu.median(),'mae_delta_mu':x.delta_mu.abs().mean(),'max_abs_delta_mu':x.delta_mu.abs().max(),'pearson_r':r,'pearson_p':rp,'spearman_rho':rho,'spearman_p':rhop,'median_abs_rank_shift':x.rank_shift.abs().median(),'max_abs_rank_shift':x.rank_shift.abs().max()})
allx=g;rho,rhop=spearmanr(allx.mu_uff,allx.mu_dreiding_uff);r,rp=pearsonr(allx.mu_uff,allx.mu_dreiding_uff)
stats.append({'scope':'ALL','n':len(allx),'mean_delta_mu':allx.delta_mu.mean(),'median_delta_mu':allx.delta_mu.median(),'mae_delta_mu':allx.delta_mu.abs().mean(),'max_abs_delta_mu':allx.delta_mu.abs().max(),'pearson_r':r,'pearson_p':rp,'spearman_rho':rho,'spearman_p':rhop,'median_abs_rank_shift':allx.rank_shift.abs().median(),'max_abs_rank_shift':allx.rank_shift.abs().max()})
sdf=pd.DataFrame(stats);sdf.to_csv(O/'wca_robustness_statistics.csv',index=False)
g.sort_values('delta_mu',key=lambda z:z.abs(),ascending=False).to_csv(O/'wca_largest_forcefield_shifts.csv',index=False)

# Figures
fig,ax=plt.subplots(figsize=(6,5))
for gas,x in g.groupby('gas'): ax.scatter(x.mu_uff,x.mu_dreiding_uff,label=gas)
lo=min(g.mu_uff.min(),g.mu_dreiding_uff.min());hi=max(g.mu_uff.max(),g.mu_dreiding_uff.max());ax.plot([lo,hi],[lo,hi],'k--',lw=1);ax.set_xlabel('UFF WCA excess chemical potential (kcal/mol)');ax.set_ylabel('DREIDING/UFF WCA excess chemical potential (kcal/mol)');ax.legend();fig.tight_layout();fig.savefig(O/'fig_wca_uff_vs_dreiding_uff.png',dpi=300);plt.close(fig)
fig,ax=plt.subplots(figsize=(7,5)); order=g.assign(a=g.delta_mu.abs()).sort_values('a').mof_id.drop_duplicates().tolist(); pos={m:i for i,m in enumerate(order)}
for gas,x in g.groupby('gas'): ax.scatter([pos[m] for m in x.mof_id],x.delta_mu,label=gas)
ax.axhline(0,color='k',lw=1);ax.set_xlabel('MOFs ordered by absolute force-field shift');ax.set_ylabel('Δμex: DREIDING/UFF − UFF (kcal/mol)');ax.legend();fig.tight_layout();fig.savefig(O/'fig_wca_forcefield_shift.png',dpi=300);plt.close(fig)

# RASPA integrity statement
rr=R/'distribution/ARCMOF_DREIDING_UFF_RASPA_BRIDGE_V1'; raspa={'complete_markers':len(list(rr.rglob('*.complete'))),'output_json':len([p for p in rr.rglob('output_*.s0.json') if 'restart' not in p.name.lower()]),'scientific_use':'directory-level sensitivity only','limitation':'shared output directories do not preserve three independent seed outputs'}
summary={'status':'FORCEFIELD_BRIDGE_ANALYSIS_COMPLETE','wca_old_all_rows':len(a_all),'wca_old_panel_rows':len(a),'wca_new_rows':len(b),'paired_seed_rows':len(p),'mof_gas_rows':len(g),'wca_statistics':stats,'raspa':raspa,'interpretation_policy':'WCA supports seed-level paired robustness analysis. RASPA is not treated as three-seed evidence because only 30 shared final outputs remain.','outputs':['wca_paired_seed_comparison_90.csv','wca_mof_gas_comparison_30.csv','wca_robustness_statistics.csv','wca_largest_forcefield_shifts.csv','fig_wca_uff_vs_dreiding_uff.png','fig_wca_forcefield_shift.png'],'complete':True}
(O/'forcefield_bridge_analysis_summary.json').write_text(json.dumps(summary,indent=2));print(json.dumps(summary,indent=2))
