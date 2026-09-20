#!/usr/bin/env python3
"""Generate the five frozen Path B manuscript figures and source-data files."""
from pathlib import Path
import json
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.patches import FancyBboxPatch, FancyArrowPatch

ROOT = Path.cwd()
DATA = ROOT / "path_b_analysis" / "path_b_final_238_authoritative.csv"
CASES = ROOT / "path_b_analysis" / "paper_figures" / "figure_case_values.csv"
PAIRS = ROOT / "path_b_analysis" / "paper_figures" / "figure_pair_values.csv"
DECISIVE = ROOT / "path_b_analysis" / "decisive_results"
ROBUST = ROOT / "path_b_analysis" / "robustness_results"
OUT = ROOT / "path_b_analysis" / "publication_figures"
SRC = OUT / "source_data"
OUT.mkdir(parents=True, exist_ok=True)
SRC.mkdir(parents=True, exist_ok=True)

for p in [DATA, CASES, PAIRS]:
    if not p.is_file():
        raise SystemExit(f"Missing required input: {p}")

d = pd.read_csv(DATA)
cases = pd.read_csv(CASES)
pairs = pd.read_csv(PAIRS)

required = ["mof_id", "lnS_wca", "delta_lnS_attr", "delta_lnS_elec", "lnS_lj_only", "lnS_full"]
missing = [c for c in required if c not in d.columns]
if missing:
    raise SystemExit(f"Missing authoritative columns: {missing}")

plt.rcParams.update({
    "font.size": 9,
    "axes.titlesize": 10,
    "axes.labelsize": 9,
    "legend.fontsize": 8,
    "xtick.labelsize": 8,
    "ytick.labelsize": 8,
    "pdf.fonttype": 42,
    "ps.fonttype": 42,
})
COL = {"wca": "#6B7280", "attr": "#2563EB", "elec": "#DC2626", "full": "#111827"}

def save(fig, stem):
    fig.savefig(OUT / f"{stem}.pdf", bbox_inches="tight")
    fig.savefig(OUT / f"{stem}.png", dpi=300, bbox_inches="tight")
    plt.close(fig)

# Figure 1: schematic plus distributions.
fig = plt.figure(figsize=(12, 6.6))
gs = fig.add_gridspec(2, 3, height_ratios=[0.72, 1.28], hspace=0.35, wspace=0.28)
ax = fig.add_subplot(gs[0, :]); ax.axis("off")
boxes = [(0.08, "WCA", "Repulsive reference", COL["wca"]), (0.40, "LJ-only", "Restore attraction", COL["attr"]), (0.72, "Full", "Restore electrostatics", COL["elec"])]
for x, title, sub, color in boxes:
    patch = FancyBboxPatch((x, 0.28), 0.20, 0.46, boxstyle="round,pad=0.02", linewidth=1.5, edgecolor=color, facecolor="white")
    ax.add_patch(patch); ax.text(x+0.10, 0.57, title, ha="center", va="center", weight="bold", color=color, fontsize=12); ax.text(x+0.10, 0.40, sub, ha="center", va="center", fontsize=9)
for x1, x2, lab in [(0.28,0.40,r"$\Delta\ln S_{attr}$"),(0.60,0.72,r"$\Delta\ln S_{elec}$")]:
    ax.add_patch(FancyArrowPatch((x1,0.51),(x2,0.51),arrowstyle="-|>",mutation_scale=14,linewidth=1.4,color="#374151")); ax.text((x1+x2)/2,0.62,lab,ha="center")
for j,(c,title,color) in enumerate([("lnS_wca", "Repulsive reference", COL["wca"]), ("delta_lnS_attr", "Attractive-LJ increment", COL["attr"]), ("delta_lnS_elec", "Electrostatic increment", COL["elec"])]):
    a=fig.add_subplot(gs[1,j]); a.hist(d[c], bins=24, color=color, alpha=.84, edgecolor="white"); med=d[c].median(); a.axvline(med,color="black",ls="--",lw=1.2,label=f"Median = {med:.2f}"); a.axvline(0,color="#9CA3AF",lw=.8); a.set_title(title); a.set_xlabel("Log-selectivity contribution"); a.set_ylabel("Number of MOFs" if j==0 else ""); a.legend(frameon=False)
fig.suptitle("Nested reference-potential decomposition of CO$_2$/N$_2$ selectivity", y=.99, fontsize=13, weight="bold")
save(fig, "fig1_decomposition")
d[["mof_id","lnS_wca","delta_lnS_attr","delta_lnS_elec"]].to_csv(SRC/"figure_1_source.csv",index=False)

# Figure 2: mechanism map.
fig,ax=plt.subplots(figsize=(7.2,5.7)); sc=ax.scatter(d.delta_lnS_attr,d.delta_lnS_elec,c=d.lnS_wca,cmap="viridis",s=34,alpha=.82,linewidth=.25,edgecolor="white")
labels={"attraction_dominated":"Attraction-dominated","electrostatic_enhanced":"Electrostatically enhanced","jointly_reinforced":"Jointly reinforced","weak_on_both_energetic_axes":"Weak on both"}
for _,r in cases[cases.case_label.isin(labels)].iterrows():
    ax.scatter(r.delta_lnS_attr,r.delta_lnS_elec,s=85,facecolor="none",edgecolor="black",linewidth=1.2); ax.annotate(labels[r.case_label],(r.delta_lnS_attr,r.delta_lnS_elec),xytext=(5,6),textcoords="offset points",fontsize=8)
ax.set_xlabel(r"Attractive-LJ increment, $\Delta\ln S_{attr}$");ax.set_ylabel(r"Electrostatic increment, $\Delta\ln S_{elec}$");ax.set_title("Mechanism space across 238 MOFs");cb=fig.colorbar(sc,ax=ax);cb.set_label(r"Repulsive reference, $\ln S_{WCA}$");fig.tight_layout();save(fig,"fig2_mechanism_map")
d[["mof_id","lnS_wca","delta_lnS_attr","delta_lnS_elec","lnS_full"]].to_csv(SRC/"figure_2_source.csv",index=False)

# Figure 3: random and grouped geometry-only R2.
met=pd.read_csv(ROBUST/"group_holdout_metrics.csv")
target_order=["lnS_wca","delta_lnS_attr","delta_lnS_elec"]
target_labels=["WCA reference","Attractive-LJ increment","Electrostatic increment"]
random=[met[(met.target==t)&(met.group_column.fillna("")=="")].r2.iloc[0] for t in target_order]
grp=[]
for t in target_order:
    x=met[(met.target==t)&(met.group_column.fillna("")!="")].r2
    grp.append((x.min(),x.max()))
fig,ax=plt.subplots(figsize=(7.2,4.8));x=np.arange(3);ax.bar(x,random,width=.50,color=[COL["wca"],COL["attr"],COL["elec"]],alpha=.86,label="Random 5-fold")
for i,(lo,hi) in enumerate(grp): ax.vlines(i,lo,hi,color="black",lw=2);ax.plot(i,lo,"v",color="black",ms=5);ax.plot(i,hi,"^",color="black",ms=5)
ax.axhline(0,color="#6B7280",lw=.9);ax.set_xticks(x,target_labels);ax.set_ylabel(r"Geometry-only cross-validated $R^2$");ax.set_title("Geometry captures accessibility and confinement, but not electrostatics");ax.text(.98,.03,"Black ranges: grouped holdouts",transform=ax.transAxes,ha="right",va="bottom",fontsize=8);fig.tight_layout();save(fig,"fig3_geometry_mechanism")
met.to_csv(SRC/"figure_3_source.csv",index=False)

# Figure 4: matched pairs.
fig,axs=plt.subplots(1,2,figsize=(11,4.8)); stage_labels=["WCA","+ Attraction","+ Electrostatics"]
for ax,(_,p) in zip(axs,pairs.iterrows()):
    A=[p.lnS_wca_a,p.lnS_wca_a+p.delta_lnS_attr_a,p.lnS_full_a];B=[p.lnS_wca_b,p.lnS_wca_b+p.delta_lnS_attr_b,p.lnS_full_b]
    ax.plot(stage_labels,A,"o-",lw=2,label=p.mof_a);ax.plot(stage_labels,B,"s-",lw=2,label=p.mof_b);ax.axhline(0,color="#9CA3AF",lw=.8);ax.set_ylabel(r"$\ln S_{CO_2/N_2}$");ax.set_title("Attractive contrast" if "attraction" in p.case_label else "Electrostatic contrast");ax.legend(frameon=False,fontsize=7)
fig.suptitle("Geometrically similar MOFs follow different energetic routes",fontsize=12,weight="bold");fig.tight_layout();save(fig,"fig4_matched_pairs")
pairs.to_csv(SRC/"figure_4_source.csv",index=False)

# Figure 5: selected rank trajectories.
sel=cases[cases.case_label.isin(["largest_attraction_rank_promotion","largest_electrostatic_rank_promotion","weak_on_both_energetic_axes","jointly_reinforced"])].copy()
fig,ax=plt.subplots(figsize=(8,5.2));stages=np.arange(3)
for _,r in sel.iterrows():
    ax.plot(stages,[r.rank_wca,r.rank_lj,r.rank_full],"o-",lw=2,label=r.case_label.replace("_"," "))
ax.set_xticks(stages,["WCA rank","LJ-only rank","Full rank"]);ax.set_ylabel("Rank among 238 MOFs");ax.invert_yaxis();ax.set_title("Interaction stages reorder candidate priorities");ax.legend(frameon=False,fontsize=7,loc="best");fig.tight_layout();save(fig,"fig5_rank_trajectories")
sel[["case_label","mof_id","rank_wca","rank_lj","rank_full","S_wca","S_lj","S_full"]].to_csv(SRC/"figure_5_source.csv",index=False)

summary={"status":"PATH_B_PUBLICATION_FIGURES_COMPLETE","figures":[f"fig{i}_{n}.pdf" for i,n in [(1,"decomposition"),(2,"mechanism_map"),(3,"geometry_mechanism"),(4,"matched_pairs"),(5,"rank_trajectories")]],"source_data_files":5,"output_directory":str(OUT.relative_to(ROOT))}
(OUT/"figure_generation_summary.json").write_text(json.dumps(summary,indent=2),encoding="utf-8");print(json.dumps(summary,indent=2))
