#!/usr/bin/env python3
"""
Integrate validated PI Gate A results into the mechanism-first manuscript.

This script:
  * regenerates Figure 2 with objective endpoint-twin evidence,
  * regenerates Figure 3 with the bootstrap predictability gap,
  * regenerates Figure 5 with population-wide rank-transition summaries,
  * patches main.tex and si.tex with only validated Gate A/B conclusions,
  * keeps the non-supporting geometry-conditioned shuffle as a bounded SI diagnostic.

It does NOT start any new molecular simulation and does NOT change the central
paper premise:
repulsive WCA baseline -> attraction creates CO2 preference ->
framework-specific electrostatic differentiation.
"""

from __future__ import annotations

import json
from pathlib import Path

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

ROOT = Path(__file__).resolve().parents[2]
AUTH = ROOT / "data" / "authoritative" / "path_b_final_238_authoritative.csv"
SRC = ROOT / "data" / "source_data"
FIG = ROOT / "figures" / "main"
MAIN = ROOT / "manuscript" / "main.tex"
SI = ROOT / "manuscript" / "si.tex"

TWIN = SRC / "pi_performance_twin_pairs.csv"
BOOT = SRC / "pi_bootstrap_summary.json"
RANK = SRC / "pi_rank_extension_summary.json"
GEO = SRC / "fig3_geometry_information_source.csv"
PERM = SRC / "pi_geometry_controlled_permutation_summary.csv"

FIG.mkdir(parents=True, exist_ok=True)

required = [AUTH, TWIN, BOOT, RANK, GEO, PERM, MAIN, SI]
missing = [str(p) for p in required if not p.exists()]
if missing:
    raise SystemExit("Missing required inputs:\n" + "\n".join(missing))

d = pd.read_csv(AUTH)
d = d.loc[d["three_axis_common"].astype(bool)].copy()
assert len(d) == 238

if "delta_lnS_elec_final" in d.columns:
    d["delta_elec"] = d["delta_lnS_elec_final"].astype(float)
else:
    d["delta_elec"] = d["delta_lnS_elec"].astype(float)
d["delta_attr"] = d["delta_lnS_attr"].astype(float)

tw = pd.read_csv(TWIN)
with open(BOOT, "r", encoding="utf-8") as f:
    boot = json.load(f)
with open(RANK, "r", encoding="utf-8") as f:
    rank = json.load(f)
geo = pd.read_csv(GEO)
perm = pd.read_csv(PERM)

best = tw.iloc[0]
q95 = float(tw["mechanism_distance"].quantile(0.95))
q99 = float(tw["mechanism_distance"].quantile(0.99))
med_d = float(tw["mechanism_distance"].median())
endpoint_cutoff = float(tw["endpoint_diff_abs_lnS_full"].max())

gap = float(boot["point"]["predictability_gap"])
gap_lo, gap_hi = map(float, boot["ci95"]["predictability_gap"])
rho_lo, rho_hi = map(float, boot["ci95"]["rho_attr_elec"])

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

def save(fig, name):
    fig.savefig(FIG / f"{name}.pdf", bbox_inches="tight")
    fig.savefig(FIG / f"{name}.png", dpi=300, bbox_inches="tight")
    plt.close(fig)

# -------------------------------------------------------------------------
# Figure 2 — continuous mechanism map + endpoint-matched mechanism distance
# -------------------------------------------------------------------------
fig, axes = plt.subplots(
    1, 2, figsize=(10.6, 4.55),
    gridspec_kw={"width_ratios": [1.18, 1.0]}
)
ax = axes[0]
sc = ax.scatter(
    d["delta_attr"], d["delta_elec"],
    c=d["lnS_full"], s=30, alpha=0.78,
    edgecolors="none", cmap="viridis"
)
a_xy = (float(best["delta_attr_i"]), float(best["delta_elec_i"]))
b_xy = (float(best["delta_attr_j"]), float(best["delta_elec_j"]))
ax.plot([a_xy[0], b_xy[0]], [a_xy[1], b_xy[1]], "--", lw=1.2)
ax.scatter([a_xy[0], b_xy[0]], [a_xy[1], b_xy[1]], s=95, facecolors="none", linewidths=1.5)
ax.annotate("endpoint twin A", a_xy, xytext=(6, 7), textcoords="offset points")
ax.annotate("endpoint twin B", b_xy, xytext=(6, -13), textcoords="offset points")
ax.set_xlabel(r"Attractive increment, $\Delta\ln S_{\rm attr}$")
ax.set_ylabel(r"Electrostatic increment, $\Delta\ln S_{\rm elec}$")
ax.set_title("(a) Continuous mechanism space")
cb = fig.colorbar(sc, ax=ax, pad=0.02)
cb.set_label(r"Final $\ln S_{\rm full}$")

ax = axes[1]
ax.hist(tw["mechanism_distance"], bins=30, alpha=0.82)
ax.axvline(med_d, ls="--", lw=1.2, label=f"median = {med_d:.3f}")
ax.axvline(q95, ls=":", lw=1.3, label=f"95th pct = {q95:.3f}")
ax.axvline(float(best["mechanism_distance"]), ls="-.", lw=1.3,
           label=f"strongest = {float(best['mechanism_distance']):.3f}")
ax.set_xlabel("Mechanism distance among endpoint-matched pairs")
ax.set_ylabel("Number of pairs")
ax.set_title("(b) Same endpoint, different interaction balance")
ax.legend(frameon=False)
fig.tight_layout()
save(fig, "fig2_mechanism_map")

pd.DataFrame({
    "metric": [
        "endpoint_match_fraction_of_all_pairs",
        "endpoint_match_count",
        "endpoint_difference_cutoff",
        "mechanism_distance_median",
        "mechanism_distance_q95",
        "mechanism_distance_q99",
        "strongest_endpoint_difference",
        "strongest_mechanism_distance",
        "strongest_attr_i",
        "strongest_elec_i",
        "strongest_attr_j",
        "strongest_elec_j",
    ],
    "value": [
        0.05, len(tw), endpoint_cutoff, med_d, q95, q99,
        float(best["endpoint_diff_abs_lnS_full"]),
        float(best["mechanism_distance"]),
        float(best["delta_attr_i"]), float(best["delta_elec_i"]),
        float(best["delta_attr_j"]), float(best["delta_elec_j"]),
    ]
}).to_csv(SRC / "fig2_pi_endpoint_twin_summary.csv", index=False)

# -------------------------------------------------------------------------
# Figure 3 — geometry information + bootstrap predictability gap
# -------------------------------------------------------------------------
fig, axes = plt.subplots(1, 2, figsize=(10.4, 4.4))

order = ["lnS_wca", "delta_attr", "delta_elec", "lnS_full"]
labels = ["WCA\nreference", "Attractive\nincrement", "Electrostatic\nincrement", "Final\nselectivity"]
gmap = geo.set_index("target")
r2s = [float(gmap.loc[k, "r2"]) for k in order]
axes[0].bar(np.arange(4), r2s, alpha=0.84)
axes[0].axhline(0, lw=0.8)
axes[0].set_xticks(np.arange(4), labels)
axes[0].set_ylabel(r"Geometry-only OOF $R^2$")
axes[0].set_title("(a) Geometry information by interaction stage")

mae_ratio_attr = float(gmap.loc["delta_attr", "mae_ratio"])
mae_ratio_elec = float(gmap.loc["delta_elec", "mae_ratio"])
axes[1].bar([0, 1], [mae_ratio_attr, mae_ratio_elec], alpha=0.84)
axes[1].axhline(1.0, ls="--", lw=1.0)
axes[1].set_xticks([0, 1], ["Attractive\nincrement", "Electrostatic\nincrement"])
axes[1].set_ylabel("MAE / median-baseline MAE")
axes[1].set_title("(b) Practical predictive gain")
axes[1].text(
    0.04, 0.96,
    rf"$G=R^2_{{attr}}-R^2_{{elec}}={gap:.3f}$" + "\n" +
    rf"95% bootstrap CI [{gap_lo:.3f}, {gap_hi:.3f}]",
    transform=axes[1].transAxes, va="top"
)
fig.tight_layout()
save(fig, "fig3_geometry_information")

pd.DataFrame({
    "metric": [
        "predictability_gap",
        "predictability_gap_ci95_low",
        "predictability_gap_ci95_high",
        "attr_r2",
        "elec_r2",
        "attr_mae_ratio",
        "elec_mae_ratio",
    ],
    "value": [
        gap, gap_lo, gap_hi,
        float(gmap.loc["delta_attr", "r2"]),
        float(gmap.loc["delta_elec", "r2"]),
        mae_ratio_attr, mae_ratio_elec,
    ]
}).to_csv(SRC / "fig3_pi_predictability_gap_summary.csv", index=False)

# -------------------------------------------------------------------------
# Figure 5 — population-wide rank transitions
# -------------------------------------------------------------------------
def descending_rank(x):
    return pd.Series(x).rank(method="average", ascending=False).to_numpy()

rw = descending_rank(d["lnS_wca"])
rl = descending_rank(d["lnS_lj_only"])
rf = descending_rank(d["lnS_full"])

fig, axes = plt.subplots(1, 2, figsize=(10.5, 4.7), sharex=True, sharey=True)
transitions = [
    (
        axes[0], rw, rl, "WCA rank", "LJ-only rank",
        rank["wca_to_lj"], "(a) Restore attraction"
    ),
    (
        axes[1], rl, rf, "LJ-only rank", "Full rank",
        rank["lj_to_full"], "(b) Restore electrostatics"
    ),
]
for ax, x, y, xl, yl, stats, title in transitions:
    ax.scatter(x, y, s=20, alpha=0.65, edgecolors="none")
    ax.plot([1, 238], [1, 238], "--", lw=1)
    ax.set_xlim(242, -3)
    ax.set_ylim(242, -3)
    ax.set_xlabel(xl)
    ax.set_ylabel(yl)
    ax.set_title(title)
    mc = stats["movement_counts"]
    ov = stats["top_overlap"]
    txt = (
        rf"$\rho={stats['spearman']:.3f}$, $\tau={stats['kendall_tau']:.3f}$" + "\n"
        + f"median |Δrank| = {stats['median_abs_rank_shift']:.0f}\n"
        + f"move ≥10%: {100*mc['10pct']['fraction']:.1f}%\n"
        + f"move ≥25%: {100*mc['25pct']['fraction']:.1f}%\n"
        + f"move ≥50%: {100*mc['50pct']['fraction']:.1f}%\n"
        + f"top-10 overlap: {ov['top10']}/10"
    )
    ax.text(0.04, 0.96, txt, transform=ax.transAxes, va="top")
fig.tight_layout()
save(fig, "fig5_rank_trajectories")

rows = []
for name, stats in [("WCA_to_LJ", rank["wca_to_lj"]), ("LJ_to_full", rank["lj_to_full"])]:
    rows.append({
        "transition": name,
        "spearman": stats["spearman"],
        "kendall_tau": stats["kendall_tau"],
        "median_abs_rank_shift": stats["median_abs_rank_shift"],
        "move_ge_10pct_fraction": stats["movement_counts"]["10pct"]["fraction"],
        "move_ge_25pct_fraction": stats["movement_counts"]["25pct"]["fraction"],
        "move_ge_50pct_fraction": stats["movement_counts"]["50pct"]["fraction"],
        "top10_overlap": stats["top_overlap"]["top10"],
        "top25_overlap": stats["top_overlap"]["top25"],
        "top10pct_overlap": stats["top_overlap"]["top10pct"],
    })
pd.DataFrame(rows).to_csv(SRC / "fig5_pi_rank_summary.csv", index=False)

# -------------------------------------------------------------------------
# Manuscript patching
# -------------------------------------------------------------------------
def replace_once(text, old, new, label):
    n = text.count(old)
    if n != 1:
        raise RuntimeError(f"{label}: expected exactly one match, found {n}")
    return text.replace(old, new, 1)

main = MAIN.read_text(encoding="utf-8")

old = """The attractive and electrostatic increments are only weakly associated across the population (Spearman $\\rho=-0.089$; Figure~\\ref{fig:map}). A framework that receives a large attractive enhancement does not necessarily receive a large electrostatic enhancement. This does not imply statistical independence; rather, once attractive interactions have established CO$_2$ preference, electrostatics introduces an additional dimension of framework-specific recognition."""
new = """The attractive and electrostatic increments are only weakly associated across the population (Spearman $\\rho=-0.089$; framework-bootstrap 95\\% interval $-0.222$ to $0.050$; Figure~\\ref{fig:map}). A framework that receives a large attractive enhancement does not necessarily receive a large electrostatic enhancement. This does not imply statistical independence; rather, once attractive interactions have established CO$_2$ preference, electrostatics introduces an additional dimension of framework-specific recognition.

The final selectivity is also mechanistically non-unique. We froze an endpoint-matching rule before inspecting pair identities: the closest 5\\% of all 28,203 framework pairs in $|\\Delta\\ln S_{\\rm full}|$ were retained (1,411 pairs; cutoff 0.0538) and then ranked by Euclidean distance in $(\\Delta\\ln S_{\\rm attr},\\Delta\\ln S_{\\rm elec})$. Within this endpoint-matched set, the median mechanism distance is 0.189 and the 95th percentile is 0.750. The strongest pair differs by only 0.016 in final $\\ln S$ but has a mechanism distance of 2.585: its attractive/electrostatic increments are 2.211/2.001 versus 4.238/0.396. Thus nearly identical endpoint performance can be assembled from sharply different balances of attraction and electrostatics without altering the population-level sequential mechanism."""
main = replace_once(main, old, new, "Results 2.2")

old = """\\caption{\\textbf{Attractive and electrostatic enhancements are weakly associated across the MOF population.} Each point is one framework. The preselected attraction- and electrostatic-matched pairs are highlighted as case studies. The absence of imposed mechanism classes is deliberate: the central result is a continuous nested interaction pathway rather than categorical ``steric'', ``dispersion'', and ``electrostatic'' MOF types.}"""
new = """\\caption{\\textbf{Final selectivity can be similar despite different attraction/electrostatic balances.} (a) Continuous mechanism space for 238 MOFs, colored by final $\\ln S_{\\rm full}$. The connected open circles mark the strongest pair selected by the frozen endpoint-matching rule. (b) Distribution of mechanism distance among the 1,411 pairs in the closest 5\\% of final-selectivity differences. The analysis remains continuous; no categorical mechanism classes are imposed.}"""
main = replace_once(main, old, new, "Figure 2 caption")

old = """The electrostatic increment behaves differently. Geometry gives $R^2=0.188$, but its MAE of 0.146 is essentially identical to the median-predictor MAE of 0.147. Thus some population-level association remains, but geometry provides almost no practical improvement in absolute error for $\\Delta\\ln S_{\\rm elec}$. Importantly, the same geometry descriptors still predict the final $\\ln S_{\\rm full}$ well ($R^2=0.845$; MAE 0.189 versus a 0.501 baseline). The role of this analysis is therefore mechanistic rather than predictive: the framework-to-framework variation in electrostatic recognition is not adequately reduced to the same pore-geometry information that captures the preceding stages."""
new = """The electrostatic increment behaves differently. Geometry gives $R^2=0.188$, but its MAE of 0.146 is essentially identical to the median-predictor MAE of 0.147. Thus some population-level association remains, but geometry provides almost no practical improvement in absolute error for $\\Delta\\ln S_{\\rm elec}$. The attraction--electrostatics predictability gap, $G=R^2_{\\rm attr}-R^2_{\\rm elec}$, is 0.686; a 5,000-resample framework bootstrap of the frozen out-of-fold predictions gives a 95\\% interval of 0.560--0.905. Importantly, the same geometry descriptors still predict the final $\\ln S_{\\rm full}$ well ($R^2=0.845$; MAE 0.189 versus a 0.501 baseline). The role of this analysis is therefore mechanistic rather than predictive: the framework-to-framework variation in electrostatic recognition is not adequately reduced to the same pore-geometry information that captures the preceding stages. We do not interpret this as geometry being irrelevant to electrostatics."""
main = replace_once(main, old, new, "Results 2.3")

old = """\\caption{\\textbf{Pore geometry captures the repulsive and attractive stages much better than electrostatic variation.} Left: five-fold cross-validated $R^2$ of geometry-only models for each interaction coordinate and the final endpoint. Right: model MAE normalized by the corresponding median-predictor MAE; a value of one indicates no practical improvement over the baseline. This is used as a mechanistic control, not as the primary objective of the study.}"""
new = """\\caption{\\textbf{Pore geometry captures the attractive stage much more strongly than electrostatic variation.} (a) Five-fold out-of-fold $R^2$ of geometry-only models for the nested coordinates and final endpoint. (b) MAE normalized by the corresponding median-predictor MAE for the attractive and electrostatic increments; a value of one indicates no practical improvement over the baseline. The attraction--electrostatics $R^2$ gap is 0.686 with a framework-bootstrap 95\\% interval of 0.560--0.905. This analysis is a mechanistic control rather than an ML benchmark.}"""
main = replace_once(main, old, new, "Figure 3 caption")

old = """The physical distinction between attraction and electrostatics has a direct consequence for candidate ordering. Restoring attraction reorganizes the population globally: the Spearman correlation between WCA and LJ-only ranks is $-0.751$, the median absolute displacement is 102 positions, and 183 of 238 frameworks move by at least 50 positions. The largest attraction-driven promotion and demotion are 234 and 226 positions.

Electrostatics acts differently. LJ-only and full ranks remain strongly correlated ($\\rho=0.905$), with a median absolute displacement of only 9 positions, and only 12 of 238 frameworks move by at least 50 positions. Nevertheless, selected frameworks respond strongly, with promotions up to 175 positions. The ranking result is therefore a consequence of the mechanism identified above: attraction creates the global CO$_2$-selective landscape, whereas electrostatics mainly differentiates particular frameworks within it."""
new = """The physical distinction between attraction and electrostatics has a direct consequence for candidate ordering. Restoring attraction reorganizes the population globally: WCA and LJ-only ranks have Spearman $\\rho=-0.751$ and Kendall $\\tau=-0.560$, with a median absolute displacement of 102 positions. Fully 89.5\\% of frameworks move by at least 10\\% of the 238-position range, 72.7\\% move by at least 25\\%, and 41.2\\% move by at least 50\\%. The top-10, top-25, and top-10\\% sets have zero overlap between WCA and LJ-only rankings. The largest attraction-driven promotion and demotion are 234 and 226 positions.

Electrostatics acts differently. LJ-only and full ranks remain strongly correlated ($\\rho=0.905$, $\\tau=0.801$), with a median absolute displacement of only 9 positions. Only 17.2\\%, 3.8\\%, and 2.5\\% of frameworks move by at least 10\\%, 25\\%, and 50\\% of the ranking range, respectively; the top-10 overlap is 6/10 and the top-25 overlap is 18/25. Nevertheless, selected frameworks respond strongly, with promotions up to 175 positions. The ranking result is therefore a consequence of the mechanism identified above: attraction creates the global CO$_2$-selective landscape, whereas electrostatics mainly differentiates particular frameworks within it."""
main = replace_once(main, old, new, "Results 2.5")

old = """\\caption{\\textbf{The attractive and electrostatic stages affect candidate ordering differently.} Every point is one MOF. Restoring attraction broadly reorganizes the WCA ranking, consistent with its role in creating CO$_2$ preference. Restoring electrostatics preserves the global LJ-only ordering much more strongly while producing large framework-specific outliers. Dashed lines denote unchanged rank.}"""
new = """\\caption{\\textbf{The attractive and electrostatic stages affect candidate ordering differently.} Each point is one MOF and the dashed line denotes unchanged rank. Restoring attraction produces broad rank reconstruction, including zero overlap of the top-10, top-25, and top-10\\% sets with the WCA ordering. Restoring electrostatics preserves the global LJ-only order much more strongly while producing selected large rank changes. Insets report correlation, displacement, and top-set statistics.}"""
main = replace_once(main, old, new, "Figure 5 caption")

anchor = """\\subsection{Matched-pair analysis}

Matched pairs were selected from the closest 2\\% of all framework pairs in a standardized nine-descriptor geometry space and then ranked by contrast in the target interaction coordinate. The retained attraction pair is M089/M192 and the retained electrostatic pair is M029/M085. M029 and M085 are mutual nearest neighbors in the matching geometry space. Pair selection preceded the spatial interpretation."""
insert = """\\subsection{Mechanism fingerprint, endpoint matching, and uncertainty}

For each framework we define the continuous mechanism vector
\\begin{equation}
\\mathbf{m}_i=(\\Delta\\ln S_{{\\rm attr},i},\\Delta\\ln S_{{\\rm elec},i}).
\\end{equation}
Endpoint-mechanism non-uniqueness was evaluated before inspecting pair identities. All 28,203 framework pairs were ranked by $|\\ln S_{{\\rm full},i}-\\ln S_{{\\rm full},j}|$; the closest 5\\% were retained, and only within this frozen subset were pairs ranked by Euclidean distance $||\\mathbf{m}_i-\\mathbf{m}_j||_2$.

Framework-level uncertainty intervals for the median increments, the attractive--electrostatic Spearman correlation, and the geometry predictability gap were obtained from 5,000 bootstrap resamples with random seed 20260920. For the predictability gap, the already-frozen out-of-fold predictions were resampled; the model was not refitted within each bootstrap.

\\subsection{Matched-pair analysis}

Matched pairs were selected from the closest 2\\% of all framework pairs in a standardized nine-descriptor geometry space and then ranked by contrast in the target interaction coordinate. The retained attraction pair is M089/M192 and the retained electrostatic pair is M029/M085. M029 and M085 are mutual nearest neighbors in the matching geometry space. Pair selection preceded the spatial interpretation."""
main = replace_once(main, anchor, insert, "Methods mechanism fingerprint")

MAIN.write_text(main, encoding="utf-8")

# -------------------------------------------------------------------------
# SI patching
# -------------------------------------------------------------------------
si = SI.read_text(encoding="utf-8")

anchor = """\\section{Population-wide rank changes}"""
addition = f"""\\section{{PI-directed existing-data robustness additions}}

\\subsection{{Framework bootstrap}}

Five thousand framework-level bootstrap resamples were used to quantify uncertainty in the existing-data additions. The median attractive increment is 2.537 (95\\% interval 2.438--2.624), the median electrostatic increment is 0.309 (0.293--0.323), and the attractive--electrostatic Spearman correlation is $-0.089$ with a 95\\% interval of $-0.222$ to 0.050. The attractive and electrostatic increments remain positive in every bootstrap resample.

For the frozen geometry-only out-of-fold predictions, $R^2_{{\\rm attr}}=0.874$ and $R^2_{{\\rm elec}}=0.188$, giving $G=0.686$. Resampling those frozen out-of-fold predictions gives a 95\\% interval of {gap_lo:.3f}--{gap_hi:.3f} for $G$. This interval quantifies uncertainty in the realized out-of-fold prediction set; the model was not refitted inside each bootstrap.

\\subsection{{Endpoint-matched mechanism non-uniqueness}}

Among all 28,203 framework pairs, the closest 5\\% in final $\\ln S$ were retained before mechanism distance was inspected. This gives 1,411 endpoint-matched pairs with a maximum endpoint difference of 0.0538. Their median mechanism distance is {med_d:.3f}, the 95th percentile is {q95:.3f}, and the 99th percentile is {q99:.3f}. The strongest selected pair differs by {float(best["endpoint_diff_abs_lnS_full"]):.3f} in final $\\ln S$ but by {float(best["mechanism_distance"]):.3f} in the two-dimensional mechanism vector.

\\subsection{{Additional leverage sensitivity}}

Removing the upper 5\\% of the attractive increment gives $\\rho_{{\\rm attr,elec}}=-0.092$; removing the upper 5\\% of the electrostatic increment gives $-0.084$; removing the upper 5\\% of final selectivity gives $-0.179$. Removing the ten largest attraction-driven rank movers gives $-0.105$, and removing the ten largest electrostatic rank movers gives $-0.083$. The weak population association is therefore not created by a few extreme frameworks.

\\subsection{{Geometry-conditioned shuffle diagnostic}}

A secondary geometry-conditioned shuffle diagnostic was also run using deterministic balanced strata in the first two principal components of the 14 frozen geometry descriptors. The first two components explain 92.87\\% of the descriptor variance. Eight-, ten-, and twelve-stratum definitions were tested, with 10,000 within-stratum shuffles for each target while the original out-of-fold prediction was held fixed.

This diagnostic did \\emph{{not}} support the stronger hypothesis that the attraction coordinate should show a larger excess-over-geometry-null association than the electrostatic coordinate. For the primary ten-stratum analysis, observed prediction--target Spearman correlations are 0.924 for attraction and 0.293 for electrostatics, whereas the corresponding shuffle-null 97.5th percentiles are 0.894 and 0.154. The same qualitative result occurs for eight and twelve strata. We therefore do not use this diagnostic as evidence that electrostatics is independent of geometry. Its role is to bound the claim: conventional geometry predicts attraction much more strongly in absolute out-of-fold performance, but some finer-grained geometry association with the electrostatic increment remains detectable.

\\section{{Population-wide rank changes}}"""
si = replace_once(si, anchor, addition, "SI PI robustness section")

SI.write_text(si, encoding="utf-8")

# -------------------------------------------------------------------------
# Static integrity checks
# -------------------------------------------------------------------------
checks = {
    "fig2_pdf": (FIG / "fig2_mechanism_map.pdf").exists(),
    "fig3_pdf": (FIG / "fig3_geometry_information.pdf").exists(),
    "fig5_pdf": (FIG / "fig5_rank_trajectories.pdf").exists(),
    "main_endpoint_twins": "The final selectivity is also mechanistically non-unique." in MAIN.read_text(encoding="utf-8"),
    "main_predictability_gap": "predictability gap" in MAIN.read_text(encoding="utf-8"),
    "main_rank_10pct": "89.5\\%" in MAIN.read_text(encoding="utf-8"),
    "si_negative_diagnostic_retained": "did \\emph{not} support" in SI.read_text(encoding="utf-8"),
}
if not all(checks.values()):
    raise RuntimeError("Static integration checks failed: " + json.dumps(checks, indent=2))

report = {
    "status": "PASS",
    "central_claim_changed": False,
    "new_simulation_used": False,
    "figures_regenerated": [2, 3, 5],
    "main_manuscript_updated": True,
    "si_updated": True,
    "gate_b_non_supporting_result_retained": True,
    "predictability_gap": gap,
    "predictability_gap_ci95": [gap_lo, gap_hi],
    "endpoint_twin_count": len(tw),
    "strongest_endpoint_difference": float(best["endpoint_diff_abs_lnS_full"]),
    "strongest_mechanism_distance": float(best["mechanism_distance"]),
    "checks": checks,
}
(ROOT / "docs" / "PI_INTEGRATION_REPORT.json").write_text(
    json.dumps(report, indent=2) + "\n", encoding="utf-8"
)

print(json.dumps(report, indent=2))
