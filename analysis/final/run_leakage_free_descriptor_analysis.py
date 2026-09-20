#!/usr/bin/env python3
from pathlib import Path
import json, re, warnings
import numpy as np
import pandas as pd
from sklearn.compose import ColumnTransformer
from sklearn.impute import SimpleImputer
from sklearn.pipeline import Pipeline
from sklearn.preprocessing import StandardScaler, OneHotEncoder
from sklearn.ensemble import ExtraTreesRegressor
from sklearn.model_selection import GroupKFold, KFold, cross_val_predict
from sklearn.metrics import r2_score, mean_absolute_error
warnings.filterwarnings('ignore')

R=Path.cwd(); IN=R/'path_b_analysis/path_b_final_238_authoritative.csv'; O=R/'path_b_analysis/leakage_free_descriptors'; O.mkdir(parents=True,exist_ok=True)
if not IN.exists(): raise SystemExit(f'Missing {IN}')
df=pd.read_csv(IN)
idcol=next((c for c in ['mof_id','MOF','name','structure'] if c in df.columns),None)
if not idcol: raise SystemExit('No MOF identifier column found')

def find_target(tokens):
 for c in df.columns:
  s=c.lower()
  if all(t in s for t in tokens) and pd.api.types.is_numeric_dtype(df[c]): return c
 return None

targets={
 'lnS_wca':find_target(['lns','wca']) or find_target(['wca','select']),
 'delta_attr':find_target(['delta','attr']) or find_target(['attr','increment']),
 'delta_elec':find_target(['delta','elec']) or find_target(['elec','increment']),
 'lnS_full':find_target(['lns','full']) or find_target(['full','select'])
}
targets={k:v for k,v in targets.items() if v}
if len(targets)<3:
 pd.DataFrame({'column':df.columns}).to_csv(O/'available_columns.csv',index=False)
 raise SystemExit('Could not identify at least three mechanism targets; see available_columns.csv')

leak=re.compile(r'(henry|widom|loading|select|lns|delta|increment|mu_|chemical|energy|wca|lj.only|full|rank|prediction|target|seed|cycle|uncert|stderr|precision|converg|runtime|volume_a3)',re.I)
geom=[c for c in df.columns if c.startswith('geom_') and not leak.search(c)]
rac=[c for c in df.columns if re.search(r'(^rac_|_rac_|mc-|lc-|f-|func|chem)',c,re.I) and not leak.search(c)]
topo=[c for c in df.columns if re.search(r'(topology|topo_|metal_cluster|functional_cluster|geometry_cluster|dimensional)',c,re.I) and not leak.search(c)]
# remove identifiers and constants
sets={'geometry':geom,'chemistry_rac':sorted(set(rac)),'topology_clusters':sorted(set(topo))}
sets['combined']=sorted(set(geom+rac+topo))
sets={k:[c for c in v if c!=idcol and df[c].nunique(dropna=True)>1] for k,v in sets.items()}

# Group on strongest available structural grouping; fallback to KFold if unavailable.
group_col=next((c for c in ['topology','topology_family','metal_cluster','geometry_cluster'] if c in df.columns and df[c].nunique()>4),None)
rows=[]; preds=[]
for tname,tcol in targets.items():
 mask=df[tcol].notna(); y=df.loc[mask,tcol].astype(float)
 for sname,cols in sets.items():
  if not cols: continue
  X=df.loc[mask,cols].copy(); nums=[c for c in cols if pd.api.types.is_numeric_dtype(X[c])]; cats=[c for c in cols if c not in nums]
  prep=ColumnTransformer([('num',Pipeline([('imp',SimpleImputer(strategy='median')),('sc',StandardScaler())]),nums),('cat',Pipeline([('imp',SimpleImputer(strategy='most_frequent')),('oh',OneHotEncoder(handle_unknown='ignore'))]),cats)])
  model=ExtraTreesRegressor(n_estimators=500,min_samples_leaf=3,max_features=0.7,random_state=20260918,n_jobs=-1)
  pipe=Pipeline([('prep',prep),('model',model)])
  if group_col:
   g=df.loc[mask,group_col].fillna('MISSING').astype(str); cv=GroupKFold(n_splits=min(5,g.nunique())); yp=cross_val_predict(pipe,X,y,cv=cv,groups=g,n_jobs=1)
   cvname=f'GroupKFold({group_col})'
  else:
   cv=KFold(n_splits=5,shuffle=True,random_state=20260918); yp=cross_val_predict(pipe,X,y,cv=cv,n_jobs=1); cvname='KFold5'
  rows.append({'target':tname,'target_column':tcol,'feature_set':sname,'n':len(y),'n_features':len(cols),'cv':cvname,'r2':r2_score(y,yp),'mae':mean_absolute_error(y,yp),'baseline_mae':mean_absolute_error(y,np.repeat(y.median(),len(y)))})
  preds += [{'mof_id':m,'target':tname,'feature_set':sname,'observed':o,'predicted':p} for m,o,p in zip(df.loc[mask,idcol],y,yp)]

pd.DataFrame(rows).to_csv(O/'descriptor_hierarchy_metrics.csv',index=False)
pd.DataFrame(preds).to_csv(O/'descriptor_hierarchy_predictions.csv',index=False)
pd.DataFrame([(k,c) for k,v in sets.items() for c in v],columns=['feature_set','column']).to_csv(O/'descriptor_feature_manifest.csv',index=False)
pd.DataFrame({'excluded_leakage_column':[c for c in df.columns if leak.search(c)]}).to_csv(O/'excluded_leakage_columns.csv',index=False)
summary={'status':'LEAKAGE_FREE_DESCRIPTOR_ANALYSIS_COMPLETE','rows':len(df),'targets':targets,'group_column':group_col,'feature_counts':{k:len(v) for k,v in sets.items()},'models':len(rows),'outputs':['descriptor_hierarchy_metrics.csv','descriptor_hierarchy_predictions.csv','descriptor_feature_manifest.csv','excluded_leakage_columns.csv']}
(O/'descriptor_analysis_summary.json').write_text(json.dumps(summary,indent=2)); print(json.dumps(summary,indent=2))
