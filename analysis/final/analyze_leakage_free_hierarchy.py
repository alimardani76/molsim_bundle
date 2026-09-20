#!/usr/bin/env python3
"""Leakage-free hierarchical descriptor analysis for electrostatic selectivity.

Uses only pre-simulation descriptor families already present in the canonical
master-derived dataset. It excludes process outputs, RASPA outputs, convergence,
replication, ranking, IDs, filenames, and all target-derived quantities.
"""
from pathlib import Path
import json, re
import numpy as np
import pandas as pd
from scipy.stats import spearmanr
from sklearn.compose import ColumnTransformer
from sklearn.impute import SimpleImputer
from sklearn.linear_model import RidgeCV
from sklearn.model_selection import RepeatedKFold, cross_validate
from sklearn.pipeline import Pipeline
from sklearn.preprocessing import OneHotEncoder, StandardScaler

ROOT=Path('.').resolve()
SRC=ROOT/'analysis_results/mechanism_dataset_primary_converged.csv'
OUT=ROOT/'analysis_results/leakage_free_descriptor_analysis'
OUT.mkdir(parents=True,exist_ok=True)
d=pd.read_csv(SRC)
y=d['delta_selectivity_elec'].to_numpy(float)

# Absolute deny-list. Any match is forbidden regardless of other naming rules.
DENY_TOKENS=(
 'process_','kh','selectivity','delta','electrostatic','factor','full','lj_only',
 'uncert','warning','converg','status','inclusion','rank','seed','cycle','elapsed',
 'widom','rep_','replicated','job','output','protocol','panel_id','filename','cif',
 'production_eligible','exclusion_reason'
)
DENY_EXACT={'mof_id','panel_index','mechanism_class','native_warning_rows','unresolved_mode_gas_rows','all_mode_gas_converged'}

def allowed(c):
 s=c.lower()
 return c not in DENY_EXACT and not any(t in s for t in DENY_TOKENS)

numeric=[c for c in d.columns if allowed(c) and pd.api.types.is_numeric_dtype(d[c]) and d[c].notna().mean()>=.80 and d[c].nunique(dropna=True)>1]
categorical=[c for c in d.columns if allowed(c) and not pd.api.types.is_numeric_dtype(d[c]) and d[c].notna().mean()>=.80 and 2<=d[c].nunique(dropna=True)<=40]

# Strict, inspectable descriptor families.
basic_geometry=[c for c in ['geom_df','geom_di','geom_density'] if c in numeric]
expanded_geometry=[c for c in numeric if c.lower().startswith('geom_')]
chem_prefix=('rac_','mc_','lc_','func_','metal_','linker_','element_','chem_','comp_')
chemistry=[c for c in numeric if c.lower().startswith(chem_prefix) and not c.lower().startswith('geom_')]
structure_cat=[c for c in categorical if any(k in c.lower() for k in ('topolog','cluster','dimensional','net','metal','linker','functional'))]

# Guard against accidental leakage surviving the family rules.
selected=set(expanded_geometry+chemistry+structure_cat)
leaks=[c for c in selected if not allowed(c)]
if leaks: raise SystemExit(f'Leakage guard failed: {leaks}')

inventory=[]
for c in d.columns:
 if c in basic_geometry:fam='basic_geometry'
 elif c in expanded_geometry:fam='expanded_geometry'
 elif c in chemistry:fam='chemistry_or_RAC'
 elif c in structure_cat:fam='topology_cluster_categorical'
 elif allowed(c):fam='allowed_but_not_whitelisted'
 else:fam='excluded_leakage_or_identifier'
 inventory.append({'column':c,'dtype':str(d[c].dtype),'nonmissing_fraction':float(d[c].notna().mean()),'unique_values':int(d[c].nunique(dropna=True)),'disposition':fam})
pd.DataFrame(inventory).to_csv(OUT/'feature_inventory.csv',index=False)

# FDR correlations only for strict pre-simulation numerical descriptors.
rows=[]
for family,features in [('basic_geometry',basic_geometry),('expanded_geometry',expanded_geometry),('chemistry_or_RAC',chemistry)]:
 for c in features:
  x=pd.to_numeric(d[c],errors='coerce');m=x.notna()&np.isfinite(x)&np.isfinite(y)
  rho,p=spearmanr(x[m],y[m]);rows.append({'family':family,'descriptor':c,'n':int(m.sum()),'spearman_rho':float(rho),'p_value':float(p)})
corr=pd.DataFrame(rows)
if len(corr):
 corr=corr.sort_values('p_value').reset_index(drop=True);rank=np.arange(1,len(corr)+1)
 corr['bh_q_value']=np.minimum.accumulate((corr.p_value*len(corr)/rank)[::-1])[::-1].clip(upper=1)
 corr['abs_rho']=corr.spearman_rho.abs();corr=corr.sort_values(['bh_q_value','abs_rho'],ascending=[True,False])
corr.to_csv(OUT/'strict_numeric_spearman_fdr.csv',index=False)

cv=RepeatedKFold(n_splits=5,n_repeats=10,random_state=20260915)
def evaluate(name,num,cat):
 if not num and not cat:return {'model':name,'status':'SKIPPED_NO_FEATURES','numeric_features':0,'categorical_features':0}
 parts=[]
 if num:parts.append(('num',Pipeline([('impute',SimpleImputer(strategy='median')),('scale',StandardScaler())]),num))
 if cat:parts.append(('cat',Pipeline([('impute',SimpleImputer(strategy='most_frequent')),('onehot',OneHotEncoder(handle_unknown='ignore',min_frequency=3))]),cat))
 model=Pipeline([('prep',ColumnTransformer(parts)),('ridge',RidgeCV(alphas=np.logspace(-4,4,33)))])
 scores=cross_validate(model,d[num+cat],y,cv=cv,scoring={'r2':'r2','mae':'neg_mean_absolute_error'},n_jobs=-1,error_score='raise')
 r2=scores['test_r2'];mae=-scores['test_mae']
 return {'model':name,'status':'COMPLETE','numeric_features':len(num),'categorical_features':len(cat),'mean_test_r2':float(r2.mean()),'sd_test_r2':float(r2.std(ddof=1)),'mean_test_mae':float(mae.mean()),'sd_test_mae':float(mae.std(ddof=1))}

specs=[
 ('basic_geometry',basic_geometry,[]),
 ('expanded_geometry',expanded_geometry,[]),
 ('chemistry_RAC_only',chemistry,[]),
 ('topology_clusters_only',[],structure_cat),
 ('basic_geometry_plus_chemistry_RAC',list(dict.fromkeys(basic_geometry+chemistry)),[]),
 ('all_strict_pre_simulation',list(dict.fromkeys(expanded_geometry+chemistry)),structure_cat)
]
models=[evaluate(*s) for s in specs]
pd.DataFrame(models).to_csv(OUT/'hierarchical_repeated_cv.csv',index=False)

summary={'status':'LEAKAGE_FREE_HIERARCHICAL_ANALYSIS_COMPLETE','mofs':len(d),'target':'delta_selectivity_elec','feature_counts':{'basic_geometry':len(basic_geometry),'expanded_geometry':len(expanded_geometry),'chemistry_RAC':len(chemistry),'topology_cluster_categorical':len(structure_cat)},'features':{'basic_geometry':basic_geometry,'expanded_geometry':expanded_geometry,'chemistry_RAC':chemistry,'topology_cluster_categorical':structure_cat},'models':models,'top_strict_associations':corr.head(20).to_dict('records') if len(corr) else [],'interpretation':'Held-out predictive associations only. No process outputs, RASPA-derived quantities, uncertainty/convergence fields, replication artifacts, ranks, IDs, or filenames are used as predictors.'}
(OUT/'leakage_free_analysis_summary.json').write_text(json.dumps(summary,indent=2),encoding='utf-8')
print(json.dumps(summary,indent=2));print('LEAKAGE-FREE DESCRIPTOR ANALYSIS COMPLETE')
