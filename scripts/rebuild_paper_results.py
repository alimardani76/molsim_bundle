#!/usr/bin/env python3
"""Rebuild the compact paper results from FINAL archived MOF-level simulation
aggregates and saved cross-validation predictions (no simulation or ML refitting).

Run: python scripts/rebuild_paper_results.py
Outputs:
  generated_results/final_selectivity_238.csv
  generated_results/paper_numbers.csv
"""
from pathlib import Path
import argparse
import math
import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
AUTHORITATIVE = ROOT / "data/authoritative/path_b_final_238_authoritative.csv"
SRC = ROOT / "data/source_data"
EXPECTED_MEDIANS = {
    "lnS_wca": -0.2037494852370207,
    "delta_lnS_attr": 2.537375623648152,
    "delta_lnS_elec": 0.3094920366015232,
}
PAIR_IDS = (
    "DB0-m15_o100_o150_f0_fsc",  # M029
    "DB0-m15_o11_o30_jea",        # M085
)

def check(condition, message):
    if not condition:
        raise ValueError(message)

def record(rows, section, quantity, value, source, units=""):
    rows.append(dict(section=section, quantity=quantity,
                     value=float(value), units=units, source_file=source))

def rebuild(output):
    d = pd.read_csv(AUTHORITATIVE)
    required = ["mof_id","lnS_wca","lnS_lj_only","lnS_full","S_wca",
                "delta_lnS_attr","delta_lnS_elec",
                "lnkh_mean_full_co2","lnkh_mean_full_n2",
                "lnkh_mean_lj_only_co2","lnkh_mean_lj_only_n2"]
    check(set(required).issubset(d.columns), "Missing authoritative fields")
    check(len(d) == 238 and d.mof_id.is_unique, "Cohort is not 238 unique MOFs")
    check(np.isfinite(d[required[1:]].to_numpy(dtype=float)).all(), "Invalid numeric result")
    err1=np.max(np.abs(d.lnS_wca+d.delta_lnS_attr-d.lnS_lj_only))
    err2=np.max(np.abs(d.lnS_lj_only+d.delta_lnS_elec-d.lnS_full))
    check(max(err1,err2) < 1e-9, "Thermodynamic telescope does not reconcile")
    for column,want in EXPECTED_MEDIANS.items():
        check(abs(d[column].median()-want)<2e-9, f"Frozen {column} median disagrees")
    src=str(AUTHORITATIVE.relative_to(ROOT))
    rows=[]
    record(rows,"cohort","n_complete",len(d),src,"MOFs")
    for key,column in (("median_lnS_wca","lnS_wca"),
                       ("median_delta_lnS_attr","delta_lnS_attr"),
                       ("median_delta_lnS_elec","delta_lnS_elec")):
        record(rows,"mechanism",key,d[column].median(),src,"log selectivity")
    record(rows,"mechanism","median_S_wca",d.S_wca.median(),src,"dimensionless")
    for step,col in (("attraction","delta_lnS_attr"),("electrostatic","delta_lnS_elec")):
        record(rows,"mechanism",f"median_{step}_multiplier",
               float(np.median(np.exp(d[col]))),src,"fold")
        record(rows,"mechanism",f"fraction_{step}_positive",
               float((d[col]>0).mean()),src,"fraction")
    record(rows,"mechanism","fraction_WCA_positive",(d.lnS_wca>0).mean(),src,"fraction")
    record(rows,"mechanism","spearman_attr_vs_elec",
           d.delta_lnS_attr.corr(d.delta_lnS_elec,method="spearman"),src)
    record(rows,"mechanism","max_telescoping_error",max(err1,err2),src)

    # Ranks use 1 = highest selectivity; the frozen manuscript uses 'min' ties.
    d=d.copy()
    for key,col in (("wca","lnS_wca"),("lj","lnS_lj_only"),("full","lnS_full")):
        d[f"rank_{key}"]=d[col].rank(method="min",ascending=False).astype(int)
    d["attraction_promotion"]=d.rank_wca-d.rank_lj
    d["electrostatic_promotion"]=d.rank_lj-d.rank_full
    record(rows,"ranking","max_attraction_promotion",d.attraction_promotion.max(),src,"places")
    record(rows,"ranking","max_electrostatic_promotion",d.electrostatic_promotion.max(),src,"places")
    record(rows,"ranking","max_attraction_demotion",(-d.attraction_promotion).max(),src,"places")
    record(rows,"ranking","max_electrostatic_demotion",(-d.electrostatic_promotion).max(),src,"places")
    top_lj=set(d.nlargest(24,"lnS_lj_only").mof_id)
    top_full=set(d.nlargest(24,"lnS_full").mof_id)
    record(rows,"ranking","top24_LJ_full_overlap",len(top_lj & top_full),src,"MOFs")
    pair=d.set_index("mof_id").loc[list(PAIR_IDS)]
    de=abs(float(pair.iloc[0].delta_lnS_elec-pair.iloc[1].delta_lnS_elec))
    record(rows,"matched_pair","M029_M085_electrostatic_contrast",de,src,"log selectivity")
    record(rows,"matched_pair","M029_M085_electrostatic_factor",math.exp(de),src,"fold")
    record(rows,"matched_pair","M029_M085_abs_WCA_difference",
           abs(float(pair.iloc[0].lnS_wca-pair.iloc[1].lnS_wca)),src,"log selectivity")

    # Verify reported CV scores from *saved out-of-fold predictions*.
    pred_file=SRC/"descriptor_hierarchy_predictions.csv"
    metrics_file=SRC/"descriptor_hierarchy_metrics.csv"
    preds=pd.read_csv(pred_file)
    metrics=pd.read_csv(metrics_file)
    check(len(preds)==238*4*2,"Unexpected number of archived OOF prediction rows")
    metrics=metrics.set_index(["target","feature_set"])
    for (target,feat),group in preds.groupby(["target","feature_set"],sort=True):
        check(len(group)==238 and group.mof_id.is_unique,"Incomplete CV panel: "+target+"/"+feat)
        yy=group.observed.to_numpy(dtype=float)
        pp=group.predicted.to_numpy(dtype=float)
        check(np.isfinite(yy).all() and np.isfinite(pp).all(),"Non-finite OOF predictions")
        r2=1-np.square(yy-pp).sum()/np.square(yy-yy.mean()).sum()
        mae=np.abs(yy-pp).mean()
        m=metrics.loc[(target,feat)]
        check(abs(r2-float(m.r2)) < 1e-10 and abs(mae-float(m.mae)) < 1e-10,
              f"OOF metric mismatch for {target}/{feat}")
        stem=f"{target}_{feat}"
        record(rows,"cv",f"r2_{stem}",r2,str(pred_file.relative_to(ROOT)),"R2")
        record(rows,"cv",f"mae_{stem}",mae,str(pred_file.relative_to(ROOT)),"MAE")
        # This baseline is a *reported archive metric*, not recomputed by OOF prediction file.
        record(rows,"cv_archived",f"baseline_mae_{stem}",m.baseline_mae,
               str(metrics_file.relative_to(ROOT)),"MAE")

    desired=["mof_id","lnS_wca","lnS_lj_only","lnS_full","S_wca",
             "delta_lnS_attr","delta_lnS_elec",
             "lnS_wca_se","lnS_lj_only_se","lnS_full_se",
             "lnkh_mean_full_co2","lnkh_mean_full_n2",
             "lnkh_mean_lj_only_co2","lnkh_mean_lj_only_n2",
             "rank_wca","rank_lj","rank_full",
             "attraction_promotion","electrostatic_promotion"]
    output.mkdir(parents=True,exist_ok=True)
    d[[c for c in desired if c in d.columns]].to_csv(
        output/"final_selectivity_238.csv",index=False)
    pd.DataFrame(rows).to_csv(output/"paper_numbers.csv",index=False)
    print(f"PASS: {len(d)} frameworks, {len(rows)} numerical results, "
          "8 OOF predictor panels, no refits.")
    print("Outputs:",output.resolve())

def main():
    ap=argparse.ArgumentParser()
    ap.add_argument("--output",type=Path,default=ROOT/"generated_results")
    opt=ap.parse_args()
    rebuild(opt.output)

if __name__=="__main__":
    main()
