#!/usr/bin/env python3
"""Redraw the ten statistical panels in JCTC v11 from archived final numeric tables.

Produces editable SVG and PDF. Panels are DATA-EQUIVALENT REDRAWS, not claimed
to be byte-for-byte copies of the unarchived original plotting implementation.
No model fitting or RASPA simulation happens in this script.
"""
from pathlib import Path
import argparse
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

ROOT = Path(__file__).resolve().parents[1]
DATA = ROOT / "data" / "authoritative" / "path_b_final_238_authoritative.csv"
SOURCE = ROOT / "data" / "source_data"
CORAL, TEAL, BLUE, PURPLE = "#D87359", "#0F8D92", "#4C7FAF", "#80659D"
PALETTE = (BLUE, CORAL, TEAL, PURPLE)

def read(name):
    path = SOURCE / name
    if not path.exists():
        raise FileNotFoundError(path)
    return pd.read_csv(path)

def save(fig, name, output):
    fig.savefig(output / (name + ".pdf"), bbox_inches="tight")
    fig.savefig(output / (name + ".svg"), bbox_inches="tight")
    plt.close(fig)
    print(name)

def theme():
    plt.rcParams.update({"font.size": 9, "axes.spines.top": False,
        "axes.spines.right": False, "pdf.fonttype": 42, "svg.fonttype": "none"})

def annotate(ax, letter, title):
    ax.set_title(title, loc="left", fontsize=10, weight="bold", pad=12)
    ax.text(-.06, 1.11, letter, transform=ax.transAxes, color=TEAL,
            weight="bold", fontsize=12)

def figure1(d, o):
    fig, axes = plt.subplots(2, 2, figsize=(10.3, 7))
    a,b,c,z=axes.flat
    for i,(col,lab) in enumerate([("lnS_wca","WCA repulsion"),
                                    ("lnS_lj_only","LJ-only attraction"),
                                    ("lnS_full","Full electrostatics")]):
        vals=d[col].to_numpy()
        a.scatter(np.full(len(vals),i)+np.random.default_rng(11+i).uniform(-.16,.16,len(vals)),
                  vals, s=6, alpha=.28, color=PALETTE[i])
        a.plot(i,np.median(vals),"_",ms=20,color="#222222",mew=2)
    a.set_xticks(range(3),["WCA\nrepulsive","LJ-only\n+ attraction","Full\n+ charges"])
    a.set_ylabel("Log CO2/N2 Henry selectivity")
    a.axhline(0,color="#888888",lw=.8,ls=":")
    annotate(a,"a","Interaction-state distributions (238 MOFs)")
    x=np.array([0,1,2])
    Y=d[["lnS_wca","lnS_lj_only","lnS_full"]].to_numpy()
    for y in Y: b.plot(x,y,color=BLUE,alpha=.055,lw=.55)
    b.plot(x,np.median(Y,axis=0),color=BLUE,lw=2.5,marker="o")
    b.set_xticks(x,["WCA","LJ-only","Full"])
    b.set_ylabel("Log Henry selectivity")
    annotate(b,"b","Trajectory through nested potentials")
    items=[d["delta_lnS_attr"].to_numpy(),d["delta_lnS_elec"].to_numpy()]
    v=c.violinplot(items,positions=[0,1],showextrema=False,showmedians=True)
    for body,color in zip(v["bodies"],[CORAL,TEAL]):
        body.set_facecolor(color);body.set_edgecolor(color);body.set_alpha(.40)
    c.set_xticks([0,1],["Attraction","Electrostatics"])
    c.set_ylabel("Increment to log selectivity")
    annotate(c,"c","Both increments positive, different scales")
    z.axis("off")
    summary=[("Median WCA preference",f"{np.exp(d.lnS_wca.median()):.3f} ×"),
             ("Median attraction factor",f"{np.exp(d.delta_lnS_attr.median()):.2f} ×"),
             ("Median electrostatic factor",f"{np.exp(d.delta_lnS_elec.median()):.2f} ×"),
             ("CO2-selective in full model",f"{(d.lnS_full>0).sum()} / {len(d)}")]
    for i,(label,value) in enumerate(summary):
        z.text(.04,.85-i*.23,label,fontsize=9,transform=z.transAxes,color="#555555")
        z.text(.58,.85-i*.23,value,fontsize=14,weight="bold",transform=z.transAxes,color="#28435C")
    annotate(z,"d","Population-level summary")
    fig.tight_layout(pad=2)
    save(fig,"fig1_interaction_pathway",o)

def figure2(d,o):
    fig,(a,b)=plt.subplots(1,2,figsize=(10,4))
    sc=a.scatter(d.delta_lnS_attr,d.delta_lnS_elec,c=d.lnS_full,
                 s=20,cmap="viridis",alpha=.82,edgecolors="none")
    fig.colorbar(sc,ax=a,label="Full log selectivity")
    rho=d[["delta_lnS_attr","delta_lnS_elec"]].corr(method="spearman").iloc[0,1]
    annotate(a,"a",f"Attraction and charge are weakly coupled (ρ={rho:.3f})")
    a.set_xlabel("Attractive increment, Δln S")
    a.set_ylabel("Electrostatic increment, Δln S")
    b.hist(d.delta_lnS_elec,bins=25,color=TEAL,edgecolor="white",lw=.6)
    b.axvline(d.delta_lnS_elec.median(),color=CORAL,lw=2,label="Median")
    b.legend(frameon=False)
    b.set_xlabel("Electrostatic selectivity increment")
    b.set_ylabel("MOF count")
    annotate(b,"b","A small median share with a long tail")
    fig.tight_layout(pad=2)
    save(fig,"fig2_mechanism_diversity",o)

def figure3(d,o):
    m=read("fig3_geometry_information_source.csv")
    order=["lnS_wca","delta_attr","delta_elec","lnS_full"]
    labels=["WCA","Attraction","Electrostatics","Final"]
    m=m.set_index("target").loc[order]
    y=np.arange(4)
    fig,(a,b)=plt.subplots(1,2,figsize=(10,3.7))
    a.barh(y,m.r2,color=PALETTE)
    a.set_yticks(y,labels);a.invert_yaxis();a.set_xlim(0,1)
    a.set_xlabel("Cross-validated R²")
    for i,v in enumerate(m.r2):a.text(v+.01,i,f"{v:.3f}",va="center",fontsize=8)
    annotate(a,"a","What pore descriptors predict")
    vals=m.mae/m.baseline_mae
    b.barh(y,vals,color=PALETTE)
    b.set_yticks(y,labels);b.invert_yaxis()
    b.axvline(1,ls="--",color="#888888",lw=1)
    b.set_xlim(0,1.2);b.set_xlabel("MAE / median-baseline MAE (lower better)")
    for i,v in enumerate(vals):b.text(v+.013,i,f"{v:.2f}",va="center",fontsize=8)
    annotate(b,"b","Electrostatics gains almost nothing over baseline")
    fig.tight_layout(pad=2)
    save(fig,"fig3_geometry_information",o)

def figure4(d,o):
    p=read("fig4_matched_pair_scalar_source.csv").set_index("mof_alias")
    g=read("fig4_per_gas_global_increment_source.csv").set_index("mof")
    pair=["M029","M085"];f,axs=plt.subplots(1,3,figsize=(11,3.7))
    a,b,c=axs
    for name,color in zip(pair,[BLUE,CORAL]):
        a.plot([0,1,2],p.loc[name,["wca","lj","full"]].astype(float),
               label=name,color=color,lw=2,marker="o")
    a.set_xticks([0,1,2],["WCA","LJ-only","Full"])
    a.set_ylabel("Log Henry selectivity");a.legend(frameon=False)
    annotate(a,"a","Two geometry-nearest neighbours")
    xx=np.arange(2);w=.34
    for j,(name,color) in enumerate(zip(pair,[BLUE,CORAL])):
        b.bar(xx+(j-.5)*w,g.loc[name,["CO2","N2"]],width=w,color=color,label=name)
    b.set_xticks(xx,["CO₂","N₂"]);b.set_ylabel("Full – LJ-only log Henry coefficient")
    b.legend(frameon=False);annotate(b,"b","The difference is CO₂-driven")
    c.axis("off")
    lines=[("Matched-space distance","0.258"),
           ("WCA log S difference",f"{abs(p.loc['M029','wca']-p.loc['M085','wca']):.3f}"),
           ("Electrostatic contrast",f"{abs(p.loc['M029','elec']-p.loc['M085','elec']):.3f}"),
           ("Multiplier contrast",f"{np.exp(abs(p.loc['M029','elec']-p.loc['M085','elec'])):.2f} ×")]
    for i,(label,value) in enumerate(lines):
        c.text(.03,.88-i*.23,label,color="#66727F",transform=c.transAxes,fontsize=9)
        c.text(.08,.74-i*.23,value,color="#263D54",weight="bold",fontsize=14,transform=c.transAxes)
    annotate(c,"c","Matched-pair audit")
    f.tight_layout(pad=2);save(f,"fig4_matched_contrast",o)

def figure5(d,o):
    r=read("fig5_rank_transition_source.csv")
    fig,axes=plt.subplots(2,2,figsize=(10,7))
    a,b,c,z=axes.flat
    a.scatter(r.rank_wca,r.rank_lj,s=8,alpha=.45,color=CORAL)
    a.plot([1,238],[1,238],ls=":",color="#888888")
    a.set(xlabel="WCA rank (1 best)",ylabel="LJ-only rank (1 best)")
    annotate(a,"a","Attraction reorganizes the order")
    b.scatter(r.rank_lj,r.rank_full,s=8,alpha=.45,color=TEAL)
    b.plot([1,238],[1,238],ls=":",color="#888888")
    b.set(xlabel="LJ-only rank (1 best)",ylabel="Full rank (1 best)")
    annotate(b,"b","Electrostatics mostly preserves order")
    selected=pd.concat([r.nlargest(2,"attraction_promotion"),
                        r.nlargest(2,"electrostatic_promotion")]).drop_duplicates("mof_id")
    for _,row in selected.iterrows():
        c.plot([0,1,2],[row.rank_wca,row.rank_lj,row.rank_full],lw=1.5,marker="o",
               label=row.mof_id)
    c.set_xticks([0,1,2],["WCA","LJ-only","Full"])
    c.invert_yaxis();c.set_ylabel("Rank (1 best)")
    c.legend(fontsize=6,frameon=False,loc="upper right")
    annotate(c,"c","Illustrative candidate trajectories")
    z.axis("off");z.text(.06,.75,
        "Ranking shifts from final source table\n"
        f"Max attraction promotion: {int(r.attraction_promotion.max())} places\n"
        f"Max electrostatic promotion: {int(r.electrostatic_promotion.max())} places",
        fontsize=11,transform=z.transAxes)
    fig.tight_layout(pad=2);save(fig,"fig5_ranking_consequences",o)

def figure6(d,o):
    f,axes=plt.subplots(1,3,figsize=(11,3.7))
    d=d.copy()
    d["family"]=d.mof_id.str.extract(r"^(DB0-m\d+)")[0]
    m=d.groupby("family")[["delta_lnS_attr","delta_lnS_elec"]].median()
    n=d.family.value_counts().sort_index()
    a,b,c=axes
    a.barh(n.index,n.to_numpy(),color=[BLUE,TEAL])
    for i,v in enumerate(n):a.text(v+2,i,f"{v}/238",va="center",fontsize=8)
    a.set_xlim(0,260);a.set_xlabel("MOF count")
    annotate(a,"a","Unequal identifier-family sizes")
    x=np.arange(len(m));w=.35
    b.bar(x-w/2,m.delta_lnS_attr,width=w,color=CORAL,label="Attraction")
    b.bar(x+w/2,m.delta_lnS_elec,width=w,color=TEAL,label="Electrostatics")
    b.set_xticks(x,m.index);b.set_ylabel("Median log-selectivity increment")
    b.legend(frameon=False,fontsize=7);annotate(b,"b","Within-family responses differ")
    n_lj=set(d.nlargest(24,"lnS_lj_only").mof_id)
    n_full=set(d.nlargest(24,"lnS_full").mof_id)
    c.barh(["LJ-only to full"],[len(n_lj&n_full)],color=TEAL,label="Retained")
    c.barh(["LJ-only to full"],[len(n_full-n_lj)],left=[len(n_lj&n_full)],color=CORAL,label="Changed")
    c.set_xlim(0,24);c.set_xlabel("Members of the full top-24")
    c.legend(frameon=False,fontsize=7);annotate(c,"c","Top-list sensitivity to each stage")
    f.tight_layout(pad=2);save(f,"fig6_family_and_shortlist",o)

def figureS1(d,o):
    fig,(a,b)=plt.subplots(1,2,figsize=(10,3.7))
    a.hist(d.lnS_wca_se,bins=25,color=BLUE,edgecolor="white")
    a.axvline(.10,color=CORAL,ls="--",label="0.10 precision threshold")
    a.legend(frameon=False,fontsize=8);a.set_xlabel("SE of WCA log selectivity")
    a.set_ylabel("MOF count");annotate(a,"a","Repulsive-reference precision")
    b.scatter(d.lnS_wca,d.lnS_wca_se,c=np.where(d.lnS_wca+1.96*d.lnS_wca_se>=0,CORAL,TEAL),
              s=14,alpha=.7)
    b.set(xlabel="WCA log-selectivity estimate",ylabel="Estimated standard error")
    annotate(b,"b","Point estimates and approximate upper bounds")
    fig.tight_layout(pad=2);save(fig,"figS1_precision_diagnostics",o)

def figureS2(d,o):
    co2=d.lnkh_mean_full_co2-d.lnkh_mean_lj_only_co2
    n2=d.lnkh_mean_full_n2-d.lnkh_mean_lj_only_n2
    fig,(a,b)=plt.subplots(1,2,figsize=(10,3.8))
    a.scatter(co2,n2,c=TEAL,s=13,alpha=.65)
    a.set(xlabel="CO₂: full – LJ-only ln Kₕ",ylabel="N₂: full – LJ-only ln Kₕ")
    a.axhline(0,ls="--",color="#999999");a.axvline(0,ls="--",color="#999999")
    annotate(a,"a","Guest-specific electrostatic response (n = 238)")
    b.hist([co2,n2],bins=28,color=[CORAL,BLUE],label=["CO₂","N₂"],alpha=.72)
    b.set(xlabel="Full – LJ-only log Henry coefficient",ylabel="MOF count")
    b.legend(frameon=False);annotate(b,"b","Electrostatic response distribution by gas")
    fig.tight_layout(pad=2);save(fig,"figS2_gas_specific_increment",o)

def figureS3(d,o):
    d=d.copy();d["family"]=d.mof_id.str.extract(r"^(DB0-m\d+)")[0]
    families=sorted(d.family.unique())
    fig,axs=plt.subplots(1,3,figsize=(11,3.7))
    for a,field,letter,title in zip(axs,["lnS_wca","delta_lnS_attr","delta_lnS_elec"],
      "abc",["WCA reference by identifier family","Attraction by identifier family",
             "Electrostatics by identifier family"]):
        groups=[d.loc[d.family==family,field].to_numpy() for family in families]
        v=a.violinplot(groups,positions=range(1,len(families)+1),showmedians=True,showextrema=False)
        for i,body in enumerate(v["bodies"]):
            body.set_facecolor([BLUE,TEAL][i]);body.set_alpha(.34)
        a.set_xticks(range(1,len(families)+1),families)
        a.set_ylabel("Log-selectivity coordinate")
        annotate(a,letter,title)
    fig.tight_layout(pad=2);save(fig,"figS3_family_details",o)

def figureS4(d,o):
    sets=[set(d.nlargest(24,c).mof_id) for c in ["lnS_wca","lnS_lj_only","lnS_full"]]
    ids=sorted(set().union(*sets),key=lambda x:(x not in sets[2],x))
    matrix=np.asarray([[int(id in subset) for subset in sets] for id in ids])
    fig,ax=plt.subplots(figsize=(7,max(6,.19*len(ids)+1)))
    from matplotlib.colors import ListedColormap
    ax.imshow(matrix,aspect="auto",interpolation="nearest",cmap=ListedColormap(["#E7ECF0",TEAL]),vmin=0,vmax=1)
    ax.set_xticks(range(3),["WCA","LJ-only","Full"]);ax.set_yticks(range(len(ids)),ids,fontsize=5)
    ax.xaxis.tick_top();ax.set_xlabel("Dark: in stage top-24",labelpad=10)
    ax.xaxis.set_label_position("top")
    ax.set_title("Membership of the union of three top-24 lists",loc="left",pad=40)
    fig.tight_layout(pad=1);save(fig,"figS4_shortlist_membership",o)

def main():
    ap=argparse.ArgumentParser()
    ap.add_argument("--output",type=Path,default=ROOT/"generated_figures"/"jctc_v11")
    ap.add_argument("--only",choices=["all","main","si"],default="all")
    opt=ap.parse_args();opt.output.mkdir(parents=True,exist_ok=True)
    d=pd.read_csv(DATA)
    assert len(d)==238 and d.mof_id.is_unique
    assert np.allclose(d.lnS_wca+d.delta_lnS_attr,d.lnS_lj_only,atol=1e-8)
    assert np.allclose(d.lnS_lj_only+d.delta_lnS_elec,d.lnS_full,atol=1e-8)
    theme()
    if opt.only in ("all","main"):
        for fn in [figure1,figure2,figure3,figure4,figure5,figure6]:fn(d,opt.output)
    if opt.only in ("all","si"):
        for fn in [figureS1,figureS2,figureS3,figureS4]:fn(d,opt.output)
    print("Saved final statistical redraws to",opt.output)

if __name__=="__main__":main()
