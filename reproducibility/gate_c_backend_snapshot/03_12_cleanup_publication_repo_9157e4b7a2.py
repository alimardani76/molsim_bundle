from __future__ import annotations

import csv
import json
import math
import shutil
import subprocess
import sys
from pathlib import Path

import numpy as np
import pandas as pd

try:
    import matplotlib.pyplot as plt
except ImportError:
    raise SystemExit("ERROR: matplotlib is required. Install with: python -m pip install matplotlib")

PROJECT_ROOT = Path(
    r"C:\Users\KaraPardazesh\Desktop\bundle\v1\mof_selectivity_mechanisms"
)
REPO = PROJECT_ROOT / "mof_selectivity_publication"

# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

def run(cmd, cwd=REPO, check=True):
    p = subprocess.run(cmd, cwd=cwd, text=True, capture_output=True)
    if check and p.returncode != 0:
        raise RuntimeError(
            f"Command failed ({p.returncode}): {' '.join(map(str, cmd))}\n"
            f"STDOUT:\n{p.stdout}\nSTDERR:\n{p.stderr}"
        )
    return p

def require(path: Path):
    if not path.exists():
        raise FileNotFoundError(path)
    return path

def git(*args, check=True):
    return run(["git", *args], check=check)

def safe_unlink(path: Path):
    if path.exists() and path.is_file():
        path.unlink()
        print(f"[removed redundant] {path.relative_to(REPO)}")

def safe_rmtree(path: Path):
    if path.exists() and path.is_dir():
        shutil.rmtree(path)
        print(f"[removed redundant] {path.relative_to(REPO)}")

def weighted_quantile(x, w, qs):
    x = np.asarray(x, float)
    w = np.asarray(w, float)
    ok = np.isfinite(x) & np.isfinite(w) & (w > 0)
    x, w = x[ok], w[ok]
    idx = np.argsort(x)
    x, w = x[idx], w[idx]
    c = (np.cumsum(w) - 0.5 * w) / w.sum()
    return np.interp(np.asarray(qs), c, x)

# ---------------------------------------------------------------------------
# Spatial-field helpers
# ---------------------------------------------------------------------------

FIELDS = REPO / "data" / "spatial_fields"

def field_candidates(mof: str, gas: str, mode: str):
    base = FIELDS / mof
    if mode == "lj_only":
        return [base / f"{mof}_{gas}_lj_only_0.75A_2048orient"]
    return [
        base / f"{mof}_{gas}_full_0.75A_adaptive",
        base / f"{mof}_{gas}_full_0.75A_32768orient",
    ]

def load_field(mof: str, gas: str, mode: str):
    folder = next(
        (p for p in field_candidates(mof, gas, mode) if (p / "field.npz").exists()),
        None,
    )
    if folder is None:
        raise FileNotFoundError(f"No final field found for {mof}/{gas}/{mode}")

    d = np.load(folder / "field.npz")
    shape = tuple(int(x) for x in d["shape"])
    logw = np.asarray(d["log_weight"], float).reshape(shape)
    rel = np.asarray(d["relative_density"], float).reshape(shape)
    cell = np.asarray(d["unit_cell_vectors_A"], float)
    return {"folder": folder, "shape": shape, "logw": logw, "rel": rel, "cell": cell}

def assert_same_grid(*fs):
    for f in fs[1:]:
        if f["shape"] != fs[0]["shape"]:
            raise RuntimeError("Spatial-field grid mismatch")
        if not np.allclose(f["cell"], fs[0]["cell"], rtol=0, atol=1e-12):
            raise RuntimeError("Spatial-field unit-cell mismatch")

def support_mask(rel, coverage=0.995):
    x = np.asarray(rel, float).ravel()
    total = float(x.sum())
    if not np.isfinite(total) or total <= 0:
        raise RuntimeError("Invalid relative-density field")
    order = np.argsort(x)[::-1]
    cum = np.cumsum(x[order]) / total
    k = int(np.searchsorted(cum, coverage, side="left"))
    mask = np.zeros(x.size, dtype=bool)
    mask[order[: k + 1]] = True
    return mask.reshape(rel.shape)

def support_union(*rels):
    m = np.zeros_like(rels[0], dtype=bool)
    for r in rels:
        m |= support_mask(r)
    return m

def weighted_z_projection(values, weights, support):
    v = np.asarray(values, float)
    w = np.where(support, np.asarray(weights, float), 0.0)
    num = np.sum(np.where(support, v, 0.0) * w, axis=2)
    den = np.sum(w, axis=2)
    out = np.full_like(num, np.nan, dtype=float)
    ok = den > 0
    out[ok] = num[ok] / den[ok]
    return out

def extent_xy(cell):
    lengths = np.linalg.norm(cell, axis=1)
    return [0.0, float(lengths[0]), 0.0, float(lengths[1])]

# ---------------------------------------------------------------------------
# Derived final source data
# ---------------------------------------------------------------------------

def build_spatial_summary(auth: pd.DataFrame):
    rows = []
    gas_rows = []

    id_map = {
        "M029": "DB0-m15_o100_o150_f0_fsc",
        "M085": "DB0-m15_o11_o30_jea",
    }

    for mof in ["M029", "M085"]:
        clj = load_field(mof, "CO2", "lj_only")
        nlj = load_field(mof, "N2", "lj_only")
        cfu = load_field(mof, "CO2", "full")
        nfu = load_field(mof, "N2", "full")
        assert_same_grid(clj, nlj, cfu, nfu)

        support = support_union(clj["rel"], nlj["rel"], cfu["rel"], nfu["rel"])
        weights = clj["rel"] + nlj["rel"] + cfu["rel"] + nfu["rel"]
        delta = (cfu["logw"] - nfu["logw"]) - (clj["logw"] - nlj["logw"])

        x = delta[support]
        w = weights[support]
        w = w / w.sum()
        q25, med, q75 = weighted_quantile(x, w, [0.25, 0.5, 0.75])

        ar = auth.loc[auth["mof_id"].eq(id_map[mof])].iloc[0]
        dco2 = float(ar["lnkh_mean_full_co2"] - ar["lnkh_mean_lj_only_co2"])
        dn2 = float(ar["lnkh_mean_full_n2"] - ar["lnkh_mean_lj_only_n2"])

        rows.append({
            "mof": mof,
            "positive_probability_weighted_fraction": float(w[x > 0].sum()),
            "weighted_mean_local_delta_elec": float(np.sum(w * x)),
            "weighted_q25_local_delta_elec": float(q25),
            "weighted_median_local_delta_elec": float(med),
            "weighted_q75_local_delta_elec": float(q75),
            "weighted_iqr_local_delta_elec": float(q75 - q25),
            "production_delta_lnKH_CO2_full_minus_LJ": dco2,
            "production_delta_lnKH_N2_full_minus_LJ": dn2,
            "production_delta_lnS_elec": float(ar["delta_lnS_elec"]),
        })

        for gas in ["CO2", "N2"]:
            lj = load_field(mof, gas, "lj_only")
            fu = load_field(mof, gas, "full")
            assert_same_grid(lj, fu)
            sg = support_union(lj["rel"], fu["rel"])
            wg = (lj["rel"] + fu["rel"])[sg]
            dg = (fu["logw"] - lj["logw"])[sg]
            wg = wg / wg.sum()
            a, b, c = weighted_quantile(dg, wg, [0.25, 0.5, 0.75])
            gas_rows.append({
                "mof": mof,
                "gas": gas,
                "positive_probability_weighted_fraction": float(wg[dg > 0].sum()),
                "weighted_mean_local_increment": float(np.sum(wg * dg)),
                "weighted_q25_local_increment": float(a),
                "weighted_median_local_increment": float(b),
                "weighted_q75_local_increment": float(c),
                "weighted_iqr_local_increment": float(c - a),
            })

    summary = pd.DataFrame(rows)
    per_gas = pd.DataFrame(gas_rows)

    out = REPO / "data" / "source_data"
    summary.to_csv(out / "spatial_quantitative_summary.csv", index=False)
    per_gas.to_csv(out / "spatial_per_gas_local_summary.csv", index=False)
    return summary, per_gas

# ---------------------------------------------------------------------------
# Publication figures
# ---------------------------------------------------------------------------

def make_fig3():
    metrics_path = require(REPO / "data" / "source_data" / "descriptor_hierarchy_metrics.csv")
    df = pd.read_csv(metrics_path)
    g = df.loc[df["feature_set"].eq("geometry")].copy()

    order = ["lnS_wca", "delta_attr", "delta_elec", "lnS_full"]
    labels = [
        r"Repulsive" + "\n" + r"$\ln S_{\mathrm{WCA}}$",
        r"Attractive" + "\n" + r"$\Delta\ln S_{\mathrm{attr}}$",
        r"Electrostatic" + "\n" + r"$\Delta\ln S_{\mathrm{elec}}$",
        r"Full" + "\n" + r"$\ln S_{\mathrm{full}}$",
    ]
    g = g.set_index("target").loc[order].reset_index()
    g["mae_ratio"] = g["mae"] / g["baseline_mae"]

    source = g[["target", "r2", "mae", "baseline_mae", "mae_ratio"]].copy()
    source.to_csv(REPO / "data" / "source_data" / "fig3_geometry_information_source.csv", index=False)

    fig, axs = plt.subplots(1, 2, figsize=(9.2, 3.8), constrained_layout=True)

    x = np.arange(len(g))
    axs[0].bar(x, g["r2"].to_numpy())
    axs[0].set_xticks(x, labels)
    axs[0].set_ylim(0, 1.0)
    axs[0].set_ylabel(r"Five-fold cross-validated $R^2$")
    axs[0].set_title("Predictive information in pore geometry")
    for i, v in enumerate(g["r2"]):
        axs[0].text(i, min(v + 0.035, 0.97), f"{v:.2f}", ha="center", va="bottom", fontsize=9)

    axs[1].bar(x, g["mae_ratio"].to_numpy())
    axs[1].axhline(1.0, linestyle="--", linewidth=1)
    axs[1].set_xticks(x, labels)
    axs[1].set_ylim(0, 1.12)
    axs[1].set_ylabel("MAE / baseline MAE")
    axs[1].set_title("Error relative to a median predictor")
    for i, v in enumerate(g["mae_ratio"]):
        axs[1].text(i, min(v + 0.035, 1.08), f"{v:.2f}", ha="center", va="bottom", fontsize=9)

    fig.suptitle(
        "Geometry captures repulsive accessibility and attraction, but adds little for electrostatics",
        fontsize=12,
    )

    out = REPO / "figures" / "main" / "fig3_geometry_information.pdf"
    fig.savefig(out, bbox_inches="tight")
    plt.close(fig)
    return out

def make_fig4(auth: pd.DataFrame, spatial_summary: pd.DataFrame):
    ids = {
        "M089": "DB0-m15_o11_o75_f0_fsc",
        "M192": "DB0-m15_o14_o77_f0_fsc",
        "M029": "DB0-m15_o100_o150_f0_fsc",
        "M085": "DB0-m15_o11_o30_jea",
    }

    vals = {}
    for alias, mid in ids.items():
        r = auth.loc[auth["mof_id"].eq(mid)].iloc[0]
        vals[alias] = {
            "wca": float(r["lnS_wca"]),
            "lj": float(r["lnS_wca"] + r["delta_lnS_attr"]),
            "full": float(r["lnS_full"]),
            "attr": float(r["delta_lnS_attr"]),
            "elec": float(r["delta_lnS_elec"]),
        }

    spatial = {}
    plotted_values = []
    for mof in ["M029", "M085"]:
        clj = load_field(mof, "CO2", "lj_only")
        nlj = load_field(mof, "N2", "lj_only")
        cfu = load_field(mof, "CO2", "full")
        nfu = load_field(mof, "N2", "full")
        assert_same_grid(clj, nlj, cfu, nfu)

        support = support_union(clj["rel"], nlj["rel"], cfu["rel"], nfu["rel"])
        weights = clj["rel"] + nlj["rel"] + cfu["rel"] + nfu["rel"]
        delta = (cfu["logw"] - nfu["logw"]) - (clj["logw"] - nlj["logw"])
        proj = weighted_z_projection(delta, weights, support)
        spatial[mof] = (proj, extent_xy(clj["cell"]))
        plotted_values.append(proj[np.isfinite(proj)])

    allx = np.concatenate(plotted_values)
    lim = float(np.quantile(np.abs(allx), 0.99))
    if lim <= 0:
        lim = 1.0

    fig = plt.figure(figsize=(10.8, 7.6), constrained_layout=True)
    gs = fig.add_gridspec(2, 2, height_ratios=[0.78, 1.15])
    a0 = fig.add_subplot(gs[0, 0])
    a1 = fig.add_subplot(gs[0, 1])
    a2 = fig.add_subplot(gs[1, 0])
    a3 = fig.add_subplot(gs[1, 1])

    stages = np.arange(3)
    stage_labels = ["WCA", "+ attraction", "+ electrostatics"]

    for alias in ["M089", "M192"]:
        v = vals[alias]
        a0.plot(stages, [v["wca"], v["lj"], v["full"]], marker="o", label=alias)
    a0.axhline(0, linewidth=0.8)
    a0.set_xticks(stages, stage_labels)
    a0.set_ylabel(r"$\ln S_{\mathrm{CO_2/N_2}}$")
    a0.set_title(r"Attraction-selected pair: $\Delta\ln S_{\rm attr}$ differs by 0.658")
    a0.legend(frameon=False)

    for alias in ["M029", "M085"]:
        v = vals[alias]
        a1.plot(stages, [v["wca"], v["lj"], v["full"]], marker="o", label=alias)
    a1.axhline(0, linewidth=0.8)
    a1.set_xticks(stages, stage_labels)
    a1.set_ylabel(r"$\ln S_{\mathrm{CO_2/N_2}}$")
    a1.set_title(r"Electrostatic pair: $\Delta\ln S_{\rm elec}$ differs by 1.860")
    a1.legend(frameon=False)

    ims = []
    for ax, mof in [(a2, "M029"), (a3, "M085")]:
        arr, extent = spatial[mof]
        im = ax.imshow(
            np.ma.masked_invalid(arr).T,
            origin="lower",
            aspect="auto",
            extent=extent,
            vmin=-lim,
            vmax=lim,
        )
        ims.append(im)
        row = spatial_summary.set_index("mof").loc[mof]
        ax.set_xlabel(r"$x$ (\AA)")
        ax.set_ylabel(r"$y$ (\AA)")
        ax.set_title(
            f"{mof}: local electrostatic selectivity increment\n"
            f"global ΔlnS = {row['production_delta_lnS_elec']:.3f}; "
            f"positive weight = {100*row['positive_probability_weighted_fraction']:.0f}%"
        )

    cbar = fig.colorbar(ims[-1], ax=[a2, a3], shrink=0.86)
    cbar.set_label(r"Probability-weighted local $\delta_{\rm elec}$")

    fig.suptitle(
        "Geometry-matched MOFs can follow different energetic routes",
        fontsize=12,
    )

    out = REPO / "figures" / "main" / "fig4_matched_electrostatic.pdf"
    fig.savefig(out, bbox_inches="tight")
    plt.close(fig)

    rows = []
    for alias in ["M089", "M192", "M029", "M085"]:
        rows.append({"mof_alias": alias, **vals[alias]})
    pd.DataFrame(rows).to_csv(
        REPO / "data" / "source_data" / "fig4_matched_pair_scalar_source.csv",
        index=False,
    )
    return out

def make_forcefield_si():
    path = require(REPO / "data" / "source_data" / "wca_mof_gas_comparison_30.csv")
    d = pd.read_csv(path)

    fig, ax = plt.subplots(figsize=(5.6, 5.0), constrained_layout=True)
    for gas, g in d.groupby("gas"):
        ax.scatter(g["mu_uff"], g["mu_dreiding_uff"], label=gas)

    lo = float(min(d["mu_uff"].min(), d["mu_dreiding_uff"].min()))
    hi = float(max(d["mu_uff"].max(), d["mu_dreiding_uff"].max()))
    pad = 0.03 * (hi - lo if hi > lo else 1.0)
    ax.plot([lo - pad, hi + pad], [lo - pad, hi + pad], linestyle="--", linewidth=1)
    ax.set_xlim(lo - pad, hi + pad)
    ax.set_ylim(lo - pad, hi + pad)
    ax.set_xlabel(r"UFF WCA $\mu^{\rm ex}$ (kcal mol$^{-1}$)")
    ax.set_ylabel(r"DREIDING-linker/UFF-Co WCA $\mu^{\rm ex}$ (kcal mol$^{-1}$)")
    ax.set_title("WCA response to the tested linker-side force-field perturbation")
    ax.legend(frameon=False)

    out = REPO / "figures" / "si" / "figS2_wca_forcefield_bridge.pdf"
    fig.savefig(out, bbox_inches="tight")
    plt.close(fig)
    return out

# ---------------------------------------------------------------------------
# Clean manuscript
# ---------------------------------------------------------------------------

MAIN_TEX = r"""
\documentclass[11pt,a4paper]{article}
\usepackage[utf8]{inputenc}
\usepackage[T1]{fontenc}
\usepackage[margin=2.05cm]{geometry}
\usepackage{microtype,booktabs,tabularx,array,amsmath,amssymb,mathtools,graphicx,float,caption,subcaption,hyperref,siunitx,bm}
\graphicspath{{../figures/main/}{../figures/si/}}
\hypersetup{colorlinks=true,linkcolor=blue,citecolor=teal,urlcolor=teal}

\title{What Makes a MOF Selective?\\From Repulsive Baselines to Attractive and Electrostatic CO$_2$/N$_2$ Recognition}
\author{Author list to be finalized}
\date{}

\begin{document}
\maketitle

\begin{abstract}
High-throughput screening ranks metal--organic frameworks by adsorption performance, but a final selectivity does not reveal which intermolecular interactions produced it. We construct a nested reference-potential decomposition of Henry-regime CO$_2$/N$_2$ selectivity across 238 MOFs, separating a repulsive Weeks--Chandler--Andersen (WCA) reference from the increments obtained by restoring Lennard--Jones attraction and then electrostatics. Repulsive accessibility alone favored N$_2$ in every framework (median $S_{\rm WCA}=0.816$). Restoring attraction reversed this preference in every case, with a median $\Delta\ln S_{\rm attr}=2.537$ (12.65-fold), whereas electrostatics supplied a further positive but strongly framework-dependent increment (median $\Delta\ln S_{\rm elec}=0.309$, 1.36-fold). The attractive and electrostatic increments were only weakly associated ($\rho=-0.089$). Conventional pore geometry contained strong predictive information for the WCA and attractive coordinates ($R^2=0.835$ and 0.874), but added almost no practical predictive value for the electrostatic increment (MAE 0.146 versus a 0.147 median-predictor baseline). A geometry-matched pair makes this distinction explicit: M029 and M085 have nearly identical repulsive selectivities but electrostatic increments of 0.141 and 2.001. Probability-supported full-minus-LJ fields show substantially stronger and more heterogeneous electrostatic reweighting in M085, while the production Henry coefficients identify the pair contrast as overwhelmingly CO$_2$-driven. Restoring attraction and electrostatics also promoted individual materials by up to 234 and 175 ranking positions. These results show that total selectivity is an endpoint rather than a complete mechanistic description and provide a controlled route for separating accessibility, attraction, and electrostatic recognition in MOF screening.
\end{abstract}

\section{Introduction}

Metal--organic frameworks (MOFs) offer broad control over pore dimensions, topology, and internal chemistry, making them attractive candidates for adsorption-based separations. Large computational databases such as ARC--MOF have made it possible to screen chemically diverse structures at scale \cite{burner2023arcmof}. Screening, however, usually compresses each material to an endpoint such as uptake, Henry coefficient, or selectivity. Two frameworks can therefore receive similar performance scores even when the underlying molecular recognition arises from different physical interactions.

For CO$_2$/N$_2$ adsorption, steric accessibility, dispersive confinement, and electrostatic interactions are all physically plausible contributors. Their effects are coupled in a conventional full force field, so a high selectivity alone cannot establish which interaction class created the preference. This is especially important when converting screening results into design rules: a pore-size trend and a local electrostatic recognition motif are not interchangeable explanations even if they produce the same final selectivity.

Here we separate these effects through a nested sequence of reference potentials. The Weeks--Chandler--Andersen construction retains the repulsive branch of the Lennard--Jones interaction \cite{weeks1971wca}; restoring the full Lennard--Jones potential introduces attraction, and restoring framework--guest charges subsequently introduces electrostatics. For each MOF, these states define
\begin{equation}
\ln S_{\rm full}
=
\ln S_{\rm WCA}
+
\Delta\ln S_{\rm attr}
+
\Delta\ln S_{\rm elec},
\end{equation}
with
$\Delta\ln S_{\rm attr}=\ln S_{\rm LJ}-\ln S_{\rm WCA}$
and
$\Delta\ln S_{\rm elec}=\ln S_{\rm full}-\ln S_{\rm LJ}$.
These are controlled finite differences between reference potentials, not uniquely separable experimental free-energy components.

We use this framework to ask four questions across a frozen 238-MOF panel. Which interaction stage first creates CO$_2$ preference? Are attractive and electrostatic enhancements redundant or distinct across materials? How much of each coordinate is already encoded by conventional pore geometry? Finally, can geometry-matched frameworks exhibit different microscopic electrostatic responses? The resulting picture is not one of three equivalent routes: repulsion establishes an N$_2$-favoring baseline, attraction creates CO$_2$ selectivity, and electrostatics subsequently differentiates frameworks in a way that pore geometry alone does not adequately capture.

\section{Results}

\subsection{Repulsion sets the baseline and attraction creates CO$_2$ preference}

The authoritative common population contains 238 MOFs evaluated in WCA, LJ-only, and full interaction states. The nested decomposition closes numerically to a maximum residual of $8.88\times10^{-16}$.

Under the WCA reference, all 238 frameworks remained N$_2$-favoring. The median $\ln S_{\rm WCA}$ was $-0.204$, corresponding to a median $S_{\rm WCA}=0.816$. Repulsive accessibility therefore does not constitute a successful CO$_2$-selectivity route in this population.

Restoring attractive Lennard--Jones interactions changed the result qualitatively. $\Delta\ln S_{\rm attr}$ was positive for every framework, with a median of 2.537, corresponding to a 12.65-fold multiplicative contribution. The LJ-only selectivity was CO$_2$-favoring in every framework. Thus attraction is the interaction stage that reverses the repulsive N$_2$ preference and establishes CO$_2$ selectivity across the panel.

\begin{figure}[H]
\centering
\includegraphics[width=\linewidth]{fig1_decomposition.pdf}
\caption{\textbf{Nested reference-potential decomposition of CO$_2$/N$_2$ selectivity.} The WCA state retains repulsive accessibility, the LJ-only state restores attractive Lennard--Jones interactions, and the full state additionally restores electrostatics. Population distributions show an N$_2$-favoring repulsive baseline, a universally positive attractive increment, and a positive but more heterogeneous electrostatic increment across the 238-MOF panel.}
\label{fig:decomp}
\end{figure}

\subsection{Electrostatics forms a distinct recognition axis}

Electrostatics provided a further positive enhancement in every framework. The median $\Delta\ln S_{\rm elec}$ was 0.309, equivalent to a 1.36-fold multiplicative contribution. Its magnitude nevertheless varied much more strongly among structures than the sign alone suggests, reaching approximately 2.0 in the strongest case.

The attractive and electrostatic increments were only weakly associated across the population (Spearman $\rho=-0.089$). A framework's gain from attractive confinement therefore provides little information about how strongly it will benefit from electrostatic interactions. The final selectivity can consequently conceal materially different energetic routes.

\begin{figure}[H]
\centering
\includegraphics[width=.82\linewidth]{fig2_mechanism_map.pdf}
\caption{\textbf{Attractive and electrostatic enhancements define weakly associated coordinates.} Each point is one MOF. Position along the horizontal and vertical axes gives the attractive and electrostatic log-selectivity increments, respectively; color denotes the repulsive WCA reference. Selected structures illustrate that similar final performance can be assembled from different combinations of interaction contributions.}
\label{fig:map}
\end{figure}

\subsection{Pore geometry captures accessibility and attraction but not electrostatic recognition}

We next asked whether conventional presimulation pore descriptors already contain the information represented by each coordinate. Leakage-free ExtraTrees models using 14 geometry descriptors gave five-fold cross-validated $R^2=0.835$ for $\ln S_{\rm WCA}$ and $R^2=0.874$ for $\Delta\ln S_{\rm attr}$. Their mean absolute errors were 0.033 and 0.177, compared with median-predictor baselines of 0.093 and 0.549, respectively.

The electrostatic increment behaved differently. Geometry produced $R^2=0.188$, but its MAE of 0.146 was essentially identical to the baseline MAE of 0.147. Thus geometry contains some statistical association with $\Delta\ln S_{\rm elec}$, but almost no practically useful improvement in absolute prediction error. Notably, geometry still predicted the final $\ln S_{\rm full}$ well ($R^2=0.845$). A model can therefore reproduce the endpoint while failing to encode an important mechanistic contribution to that endpoint.

\begin{figure}[H]
\centering
\includegraphics[width=\linewidth]{fig3_geometry_information.pdf}
\caption{\textbf{Conventional pore geometry contains strongly unequal information about the interaction coordinates.} Left: five-fold cross-validated $R^2$ of leakage-free geometry-only ExtraTrees models. Right: model MAE normalized by the corresponding median-predictor MAE; a value of one indicates no practical improvement over the baseline. Geometry strongly captures the WCA, attractive, and full-selectivity responses, whereas the electrostatic MAE remains essentially at baseline.}
\label{fig:geometry}
\end{figure}

\subsection{Geometry-matched frameworks expose different energetic mechanisms}

Matched frameworks provide direct counterexamples to the assumption that similar pore geometry implies similar recognition physics. In the attraction-selected pair, M089 and M192 have similar WCA selectivities but attractive increments of 2.626 and 3.283, a difference of 0.658 (1.93-fold in the corresponding multiplicative contribution). This scalar comparison establishes an attraction contrast, but the available LJ spatial maps are not used to localize that increment because corresponding WCA spatial fields were not computed.

The electrostatic pair provides a stronger mechanistic case. M029 and M085 have nearly identical repulsive log-selectivities ($-0.201$ and $-0.198$), yet their electrostatic increments are 0.141 and 2.001, a difference of 1.860 (6.42-fold in the multiplicative electrostatic contribution). Probability-supported full-minus-LJ maps show that M029 remains close to locally neutral over much of the insertion-supported pore region, whereas M085 develops a substantially more positive and heterogeneous electrostatic selectivity response. The probability-weighted median local increment is 0.042 for M029 and 0.534 for M085; 57.6\% and 79.1\% of the corresponding probability weight lies in regions with positive local electrostatic selectivity increments.

The production Henry coefficients identify which gas drives the contrast. Restoring electrostatics changes $\ln K_H$ for CO$_2$ by 0.162 in M029 and 2.032 in M085, whereas the N$_2$ changes are only 0.021 and 0.032. The between-framework difference is therefore 1.871 log units for CO$_2$ but only 0.011 for N$_2$; the much larger electrostatic selectivity of M085 is overwhelmingly CO$_2$-driven.

\begin{figure}[H]
\centering
\includegraphics[width=\linewidth]{fig4_matched_electrostatic.pdf}
\caption{\textbf{Geometry-matched MOFs can follow different energetic routes.} Top: trajectories through the WCA, LJ-only, and full interaction states for the attraction-selected and electrostatic matched pairs. Bottom: probability-supported full-minus-LJ local electrostatic selectivity increments for M029 and M085. The electrostatic pair begins from nearly identical repulsive selectivities, but M085 develops a much larger and more spatially heterogeneous electrostatic enhancement. Spatial fields are used to localize the perturbation; authoritative global increments are taken from the production Henry calculations.}
\label{fig:matched}
\end{figure}

\subsection{Interaction stages reorder screening priorities}

Mechanistic differences are large enough to alter which materials appear promising. Restoring attraction promoted individual frameworks by as many as 234 positions among the 238 candidates, while restoring electrostatics produced promotions of up to 175 positions. The largest corresponding demotions were 226 and 55 positions.

These shifts show why endpoint ranking and mechanism are not interchangeable. A framework can appear weak under a repulsive reference and become competitive after attraction is restored, or remain relatively ordinary in the LJ-only state and rise sharply once electrostatic recognition is included.

\begin{figure}[H]
\centering
\includegraphics[width=.80\linewidth]{fig5_rank_trajectories.pdf}
\caption{\textbf{Successive interaction stages reorder candidate priorities.} Representative rank trajectories illustrate how restoring attraction and electrostatics changes the ordering of the same 238 frameworks. The trajectories are consequences of the interaction-dependent selectivities, not an independent microscopic mechanism.}
\label{fig:ranks}
\end{figure}

\subsection{The repulsive baseline is robust to the tested linker-side force-field perturbation}

A preselected 15-MOF bridge panel tested whether the WCA conclusion was fragile to the framework Lennard--Jones parameterization. Supported linker and main-group framework atoms were changed from UFF to DREIDING while cobalt retained UFF; guest models, structures, charges, mixing rules, and the WCA construction were otherwise unchanged.

The mean absolute change in WCA excess chemical potential was 0.0227~kcal\,mol$^{-1}$ for CO$_2$ and 0.0178~kcal\,mol$^{-1}$ for N$_2$. Spearman rank correlations were 0.996 and 1.000, and the maximum rank displacement was one position. The repulsive-accessibility conclusion is therefore stable to this specific linker-side perturbation. This test does not establish force-field independence of the attractive or electrostatic coordinates, nor does it test polarization, framework flexibility, or all MOF chemistries.

\section{Discussion}

The nested reference potentials reveal a simple population-level sequence. Repulsive accessibility establishes an N$_2$-favoring baseline; attractive Lennard--Jones interactions are sufficient to reverse that baseline and create CO$_2$ preference throughout the panel; electrostatics then supplies an additional enhancement whose magnitude varies strongly among structures. In this sense, the decomposition refines the common idea of multiple equivalent ``routes'' to selectivity: for the present CO$_2$/N$_2$ panel, the stages play different roles rather than acting as interchangeable alternatives.

The contrast between geometry and electrostatics is particularly consequential. Conventional pore descriptors encode most of the information required to reproduce the repulsive and attractive coordinates, and they also predict the final selectivity well. Yet the same descriptors add essentially no improvement over a median predictor for the absolute electrostatic error. Good endpoint prediction therefore does not imply that a representation captures every mechanism contributing to the endpoint. For screening and design, this distinction matters because optimizing geometric confinement alone cannot specify which structures will obtain unusually large additional electrostatic gains.

The matched electrostatic pair converts that population-level statement into a microscopic example. M029 and M085 begin from nearly identical repulsive selectivities but diverge sharply when electrostatics are restored. The global contrast is almost entirely associated with CO$_2$, and the spatial fields show that the larger M085 response is distributed heterogeneously through the probability-supported pore region. These maps do not identify a specific chemical motif; doing so would require a validated chemistry-aware representation or a separate local-structure analysis. Their role here is narrower: they demonstrate that the scalar electrostatic increment corresponds to a genuinely different spatial response rather than merely a change in a final number.

The ranking analysis gives the mechanism decomposition a practical screening consequence. Interaction stages can move individual candidates across most of the 238-material ranking. A final rank therefore answers ``which material performs best under this model,'' whereas the mechanism coordinates answer ``why did it become competitive.'' Preserving both questions provides more transferable information than endpoint screening alone.

The scope is deliberately bounded. The results describe Henry-regime thermodynamic recognition at 298~K in rigid frameworks under the selected molecular models. The WCA coordinate was insensitive to one preselected DREIDING-linker/UFF-Co perturbation, but equivalent force-field sensitivity has not been established for the LJ-only or electrostatic increments. The study does not address finite-pressure mixtures, transport, framework flexibility, or transfer to other gas pairs.

\section{Methods}

\subsection{Framework panel and molecular models}

The starting panel contained 240 ARC--MOF structures \cite{burner2023arcmof}; 238 formed the authoritative common population after the all-mode/all-gas requirements used for the nested comparison. Frameworks were treated as rigid at 298~K. Full and LJ-only Henry calculations used RASPA3 \cite{ran2024raspa3,dubbeldam2016raspa}. The baseline framework Lennard--Jones model used element-level UFF parameters \cite{rappe1992uff} with Lorentz--Berthelot mixing. Full calculations included framework--guest electrostatics; LJ-only calculations set framework and guest pseudoatom charges to zero while preserving the Lennard--Jones definitions.

\subsection{Nested reference potentials}

For gas $g$ and interaction state $x$, the Henry-regime selectivity is
\begin{equation}
\ln S_x
=
\ln K_{H,\mathrm{CO_2}}^x
-
\ln K_{H,\mathrm{N_2}}^x,
\qquad
x\in\{\mathrm{WCA},\mathrm{LJ},\mathrm{full}\}.
\end{equation}
The attractive and electrostatic coordinates are finite differences between successive states:
\begin{align}
\Delta\ln S_{\rm attr} &= \ln S_{\rm LJ}-\ln S_{\rm WCA},\\
\Delta\ln S_{\rm elec} &= \ln S_{\rm full}-\ln S_{\rm LJ}.
\end{align}

For the WCA state, each Lennard--Jones pair potential was truncated at
$r_{c,ij}=2^{1/6}\sigma_{ij}$ and shifted to zero at the cutoff. Molecular Widom insertion \cite{widom1963} was used to evaluate the repulsive reference. Independent equal-length WCA seeds were pooled in Boltzmann-factor space before conversion back to the excess chemical potential. Detailed seed counts, precision diagnostics, and exclusion tests are provided in the Supporting Information.

\subsection{Geometry information analysis}

The leakage-free geometry analysis used 14 conventional presimulation geometry descriptors. Target-derived adsorption quantities, identifiers, filenames, cluster labels, and other leakage-prone variables were excluded before model fitting. ExtraTrees regressors used 500 trees, a minimum leaf size of three, $0.7$ maximum-feature fraction, and five-fold shuffled cross-validation with random seed 20260918. Prediction was used as an information-content diagnostic rather than as the primary objective. Performance was evaluated with $R^2$, MAE, and the MAE of a median predictor.

\subsection{Matched-pair analysis}

Matched frameworks were selected to be close in conventional geometry space while displaying a large contrast along a target energetic coordinate. The frozen attraction pair was M089/M192 and the frozen electrostatic pair was M029/M085. Pair selection preceded the spatial interpretation.

\subsection{Spatial insertion fields}

Spatial fields were computed on a fixed 0.75~\AA{} grid using the validated fixed-position RASPA/RASPALIB energy evaluator and the same rigid guest geometries, framework structures, Lennard--Jones parameters, charge assignments, mixing rules, and temperature as the production calculations. At each position $\mathbf r$,
\begin{equation}
w(\mathbf r)=
\left\langle
\exp[-\beta U(\mathbf r,\Omega)]
\right\rangle_{\Omega}.
\end{equation}
LJ-only fields used 2,048 orientations per voxel. Full fields used 16,384 or 32,768 orientations; targeted nested 16k$\rightarrow$32k$\rightarrow$64k checks at probability-important voxels changed the tested values by at most approximately 0.04\%.

For the electrostatic pair, the local selectivity perturbation was
\begin{equation}
\delta_{\rm elec}(\mathbf r)
=
\left[
\ln w_{\rm full}^{\mathrm{CO_2}}(\mathbf r)
-\ln w_{\rm full}^{\mathrm{N_2}}(\mathbf r)
\right]
-
\left[
\ln w_{\rm LJ}^{\mathrm{CO_2}}(\mathbf r)
-\ln w_{\rm LJ}^{\mathrm{N_2}}(\mathbf r)
\right].
\end{equation}
Visualization and local summary statistics were restricted to the union of the minimal voxel sets carrying 99.5\% of the insertion probability of the relevant fields. Local values were probability weighted; inaccessible framework-overlap voxels were not averaged into logarithmic maps. The 0.75~\AA{} fields are used for localization rather than as replacement quadrature estimators for the production Henry coefficients.

\subsection{Statistical robustness and force-field bridge}

The principal signs and the weak attractive--electrostatic association were recomputed after excluding the four lower-precision WCA cases and after removing the upper 1\% of each response separately. A preselected 15-MOF bridge additionally replaced supported linker/main-group UFF parameters with DREIDING while retaining UFF for cobalt; all other WCA definitions were preserved. This bridge was used only to test sensitivity of the repulsive-accessibility coordinate to the specified linker-side Lennard--Jones perturbation.

\section{Conclusions}

Across 238 MOFs, repulsive accessibility alone favored N$_2$, attractive interactions universally reversed this baseline toward CO$_2$, and electrostatics supplied a further framework-dependent enhancement. Conventional pore geometry captured the repulsive and attractive coordinates strongly but provided almost no practical reduction in electrostatic prediction error relative to a baseline predictor.

The M029/M085 matched pair shows why retaining the individual coordinates matters: nearly identical repulsive selectivities evolve into a 1.860-unit difference in electrostatic selectivity, driven overwhelmingly by CO$_2$ and accompanied by a markedly different spatial response. Successive interaction stages also reorder screening priorities by up to hundreds of positions.

Total selectivity is therefore a useful endpoint but an incomplete mechanistic descriptor. Nested reference potentials provide a controlled way to preserve the physical route by which a MOF becomes selective, separating accessibility, attraction, and electrostatic recognition without requiring those contributions to be treated as uniquely separable experimental free energies.

\bibliographystyle{unsrt}
\bibliography{references}
\end{document}
""".strip() + "\n"

SI_TEX = r"""
\documentclass[11pt,a4paper]{article}
\usepackage[utf8]{inputenc}
\usepackage[T1]{fontenc}
\usepackage[margin=2cm]{geometry}
\usepackage{microtype,booktabs,longtable,tabularx,array,amsmath,graphicx,float,caption,hyperref,siunitx}
\graphicspath{{../figures/si/}{../figures/main/}}
\hypersetup{colorlinks=true,linkcolor=blue,citecolor=teal,urlcolor=teal}

\title{Supporting Information\\What Makes a MOF Selective? From Repulsive Baselines to Attractive and Electrostatic CO$_2$/N$_2$ Recognition}
\author{Author list to be finalized}
\date{}

\begin{document}
\maketitle
\tableofcontents

\section{Panel definition and evidence scope}

The starting ARC--MOF panel contained 240 structures. The authoritative common population contains 238 MOFs with the required WCA, LJ-only, and full CO$_2$/N$_2$ quantities. All population-level mechanism statistics in the manuscript refer to these 238 structures.

The final evidence set includes the nested mechanism table, valid precision and upper-tail exclusion tests, leakage-free geometry analysis, frozen matched pairs, probability-supported electrostatic spatial fields, ranking trajectories, and the bounded 15-MOF WCA force-field bridge. Exploratory analyses that were invalid, unrecoverable, or outside the final claim scope are not used as scientific evidence and are listed in the repository exclusion ledger rather than reproduced here.

\section{Simulation definitions}

Full and LJ-only Henry calculations used RASPA3 at 298~K with rigid frameworks. The van der Waals cutoff was 12~\AA{}, Lorentz--Berthelot mixing was used, and tail corrections were disabled. Full calculations used Ewald electrostatics; LJ-only calculations set framework and guest pseudoatom charges to zero while preserving all Lennard--Jones definitions.

The WCA reference retained the repulsive branch of each framework--guest Lennard--Jones pair,
\begin{equation}
r_{c,ij}=2^{1/6}\sigma_{ij},
\end{equation}
with the potential shifted to zero at the cutoff.

\section{WCA sampling, pooling, and precision}

The final WCA archive contains 2,052 jobs. Independent equal-length seeds were pooled through
\begin{equation}
B_k=\exp[-\mu_k^{\rm ex}/(RT)],
\end{equation}
followed by averaging in Boltzmann-factor space and conversion through $\ln\bar B$. The final population contains 178 MOFs with three seeds per gas, 38 with six, three with nine, and 19 with twelve.

Four frameworks retained $\mathrm{SE}(\ln S_{\rm WCA})>0.10$. Excluding them did not change any principal qualitative conclusion.

\section{Exclusion robustness}

\begin{table}[H]
\centering
\caption{Sensitivity of the principal mechanism summaries to predefined exclusions.}
\begin{tabular}{lrrrrr}
\toprule
Population & $n$ & med. $\ln S_{\rm WCA}$ & med. $\Delta\ln S_{\rm attr}$ & med. $\Delta\ln S_{\rm elec}$ & $\rho_{\rm attr,elec}$\\
\midrule
All & 238 & -0.204 & 2.537 & 0.309 & -0.089\\
Precision pass & 234 & -0.203 & 2.510 & 0.309 & -0.094\\
Exclude top 1\% WCA & 235 & -0.207 & 2.551 & 0.309 & -0.082\\
Exclude top 1\% attraction & 235 & -0.203 & 2.510 & 0.309 & -0.091\\
Exclude top 1\% electrostatics & 235 & -0.204 & 2.542 & 0.308 & -0.087\\
Exclude top 1\% full selectivity & 235 & -0.203 & 2.510 & 0.308 & -0.117\\
\bottomrule
\end{tabular}
\end{table}

In every listed population the WCA-positive fraction remained zero, while the attractive- and electrostatic-positive fractions remained one.

\section{Leakage-free geometry models}

The final information-content analysis used 14 presimulation geometry descriptors. ExtraTrees regressors contained 500 trees, minimum leaf size three, and maximum-feature fraction 0.7. Five-fold shuffled cross-validation used random seed 20260918.

\begin{table}[H]
\centering
\caption{Five-fold cross-validated geometry-only prediction. Baseline MAE is the error of a median predictor.}
\begin{tabular}{lrrr}
\toprule
Target & $R^2$ & MAE & Baseline MAE\\
\midrule
$\ln S_{\rm WCA}$ & 0.835 & 0.033 & 0.093\\
$\Delta\ln S_{\rm attr}$ & 0.874 & 0.177 & 0.549\\
$\Delta\ln S_{\rm elec}$ & 0.188 & 0.146 & 0.147\\
$\ln S_{\rm full}$ & 0.845 & 0.189 & 0.501\\
\bottomrule
\end{tabular}
\end{table}

Target-derived adsorption quantities, identifiers, filenames, cluster labels, and other leakage-prone columns were excluded before fitting. The complete feature manifest and prediction table are provided as source data.

\section{Matched-pair values}

\begin{table}[H]
\centering
\caption{Frozen geometry-matched case studies.}
\begin{tabular}{lrrrr}
\toprule
MOF & $\ln S_{\rm WCA}$ & $\Delta\ln S_{\rm attr}$ & $\Delta\ln S_{\rm elec}$ & $\ln S_{\rm full}$\\
\midrule
M089 & -0.321 & 2.626 & 0.466 & 2.770\\
M192 & -0.248 & 3.283 & 0.140 & 3.175\\
M029 & -0.201 & 2.839 & 0.141 & 2.779\\
M085 & -0.198 & 2.211 & 2.001 & 4.014\\
\bottomrule
\end{tabular}
\end{table}

The attraction-selected pair differs by 0.658 in $\Delta\ln S_{\rm attr}$. The electrostatic pair differs by 1.860 in $\Delta\ln S_{\rm elec}$. LJ spatial maps for M089/M192 are not used to localize attraction because corresponding WCA spatial fields were not computed.

\section{Spatial mapping and orientation checks}

Spatial insertion fields used a fixed 0.75~\AA{} grid and the validated fixed-position RASPA/RASPALIB evaluator. LJ-only fields used 2,048 orientations per voxel. Full electrostatic fields used 16,384 or 32,768 orientations. Targeted nested 16k$\rightarrow$32k$\rightarrow$64k checks at probability-important voxels gave maximum observed 16k$\rightarrow$32k changes of approximately 0.04\% for M029/CO$_2$, 0.01\% for M029/N$_2$, and 0.002\% for M085/N$_2$; subsequent 32k$\rightarrow$64k changes were smaller.

Probability support was defined as the minimal voxel set containing 99.5\% of insertion probability mass. Local log-selectivity quantities were summarized only within the union of the relevant supports and were weighted by normalized insertion probability. This avoids allowing framework-overlap voxels with effectively zero insertion probability to dominate logarithmic averages.

\begin{table}[H]
\centering
\caption{Probability-supported summaries for the electrostatic matched pair. These are descriptive spatial statistics, not replacements for the production Henry coefficients.}
\begin{tabular}{lrr}
\toprule
Quantity & M029 & M085\\
\midrule
Positive probability-weight fraction & 0.576 & 0.791\\
Weighted mean local $\delta_{\rm elec}$ & 0.093 & 1.025\\
Weighted median local $\delta_{\rm elec}$ & 0.042 & 0.534\\
Weighted IQR local $\delta_{\rm elec}$ & 0.614 & 2.019\\
\bottomrule
\end{tabular}
\end{table}

\subsection{Per-gas electrostatic response}

\begin{figure}[H]
\centering
\includegraphics[width=.95\linewidth]{figS1_electrostatic_per_gas.pdf}
\caption{\textbf{Per-gas probability-supported electrostatic increments for the M029/M085 pair.} Full-minus-LJ fields are shown separately for CO$_2$ and N$_2$. These maps support the scalar production result that the large between-framework electrostatic-selectivity contrast is dominated by CO$_2$.}
\label{fig:si_spatial_gas}
\end{figure}

From the production Henry coefficients, $\Delta\ln K_H^{\rm full-LJ}$ for CO$_2$ is 0.162 in M029 and 2.032 in M085. The corresponding N$_2$ values are 0.021 and 0.032. Thus the M085--M029 difference is 1.871 log units for CO$_2$ but only 0.011 for N$_2$.

\section{WCA force-field bridge}

The 15-MOF bridge panel was frozen before the alternative-force-field results were inspected. Supported linker and main-group framework atoms were reassigned from UFF to DREIDING while cobalt retained UFF. Guest models, charges, structures, and the WCA construction were unchanged.

\begin{table}[H]
\centering
\caption{WCA sensitivity to the tested DREIDING-linker/UFF-Co perturbation.}
\begin{tabular}{lrrrr}
\toprule
Gas & $n$ & MAE $|\Delta\mu^{\rm ex}|$ & Max. $|\Delta\mu^{\rm ex}|$ & Spearman $\rho$\\
 & & (kcal mol$^{-1}$) & (kcal mol$^{-1}$) & \\
\midrule
CO$_2$ & 15 & 0.0227 & 0.0563 & 0.996\\
N$_2$ & 15 & 0.0178 & 0.0507 & 1.000\\
\bottomrule
\end{tabular}
\end{table}

\begin{figure}[H]
\centering
\includegraphics[width=.62\linewidth]{figS2_wca_forcefield_bridge.pdf}
\caption{\textbf{WCA excess chemical potentials under the baseline and tested alternative framework Lennard--Jones assignments.} The dashed line denotes equality. The bridge supports robustness of the repulsive-accessibility coordinate to this bounded linker-side perturbation only.}
\label{fig:si_ff_bridge}
\end{figure}

\section{Ranking changes}

Across the 238-MOF population, the maximum promotion after restoring attraction was 234 positions and the maximum demotion was 226. Restoring electrostatics produced a maximum promotion of 175 positions and a maximum demotion of 55. Complete per-MOF values are provided with the source data.

\section{Source data and reproducibility}

The publication repository contains the authoritative 238-MOF mechanism table, figure source data, descriptor metrics and predictions, spatial-field files for the retained electrostatic pair, WCA bridge tables, final analysis scripts, protocol records, and SHA-256 manifests. The 0.75~\AA{} spatial fields are retained for localization; production Henry coefficients remain authoritative for integrated thermodynamic increments.

\bibliographystyle{unsrt}
\bibliography{references}
\end{document}
""".strip() + "\n"

STORY_MD = r"""
# Publication story — frozen cleanup version

## One-sentence thesis

CO2/N2 selectivity in this 238-MOF panel is built in stages: repulsion sets an N2-favoring baseline, attraction creates CO2 preference, and electrostatics supplies an additional framework-dependent recognition channel that conventional pore geometry largely fails to encode.

## Claim hierarchy

### Claim 1 — baseline
All 238 WCA references favor N2. Repulsion is therefore a baseline, not a successful CO2-selectivity route in this panel.

### Claim 2 — creation of selectivity
Restoring Lennard-Jones attraction makes every framework CO2-selective. This is the stage that creates CO2 preference.

### Claim 3 — distinct electrostatic information
Electrostatics increases selectivity in every framework, but the increment varies strongly and is weakly associated with the attractive increment.

### Claim 4 — representation boundary
Pore geometry predicts the WCA and attractive coordinates strongly, but its electrostatic MAE is essentially identical to a median-predictor baseline. A representation can therefore predict the final endpoint while missing an important mechanism that contributes to it.

### Claim 5 — microscopic counterexample
M029 and M085 have nearly identical WCA selectivities but electrostatic increments of 0.141 and 2.001. The difference is overwhelmingly CO2-driven, and probability-supported full-minus-LJ maps show a much stronger and more heterogeneous local response in M085.

### Claim 6 — screening consequence
Restoring attraction or electrostatics can reorder a material by most of the 238-candidate ranking. Endpoint rank and physical mechanism answer different questions.

### Bounded robustness claim
The WCA coordinate is stable to the tested DREIDING-linker/UFF-Co perturbation on a frozen 15-MOF panel. This is not a force-field-independence claim for attraction or electrostatics.

## Figure logic

1. Figure 1 — What is being decomposed and what does each population distribution look like?
2. Figure 2 — Are attraction and electrostatics the same energetic axis? No.
3. Figure 3 — Which coordinates are already encoded by pore geometry? WCA/attraction yes; electrostatics practically no.
4. Figure 4 — Can similar geometry hide different mechanisms? Yes; M029/M085 provides the strongest spatially resolved counterexample.
5. Figure 5 — Does preserving mechanism matter for screening? Yes; rankings change drastically.
6. SI — convergence/precision, per-gas spatial attribution, WCA force-field bridge, full source tables.

## What is intentionally not part of the paper

- circular permutation analysis;
- invalid automated chemistry/RAC hierarchy;
- family leave-one-out claim;
- M089/M192 LJ map as spatial evidence of attraction;
- quantitative full/LJ alternative-force-field robustness from the failed shared-output bridge;
- CO2/CH4 transfer;
- finite-pressure, transport, flexibility, or universal gas-pair claims;
- narrative about debugging failures.

These remain in project/Git provenance where useful, but they do not belong in publication prose.
""".strip() + "\n"

EXCLUDED_MD = r"""
# Excluded analyses and claims

This ledger prevents exploratory or invalid results from re-entering the publication.

| Item | Publication status | Reason |
|---|---|---|
| Circular permutation analysis | Excluded | Algebraically breaks the nested identity and is not a valid negative control. |
| Automated chemistry/RAC hierarchy | Excluded | Feature set contained target-derived/leakage-prone quantities and identifiers; reported high scores are not interpretable as chemistry information. |
| Family leave-one-out | Not completed; excluded | Retained CSV is empty and the authoritative table lacks the family columns required by the original script. |
| Alternative-force-field RASPA full/LJ bridge | Not evidence | Shared output directories did not preserve recoverable seed-resolved Henry observables. |
| M089/M192 LJ spatial maps as attraction localization | Excluded from mechanistic claim | Without WCA spatial maps, LJ maps cannot isolate the spatial origin of Delta ln S_attr. |
| Grid-integrated spatial fields as replacement Henry coefficients | Excluded | The 0.75 A fields are localization estimators and do not exactly reproduce production scalar increments. |
| CO2/CH4 transfer | Not tested | Outside current CO2/N2 claim scope. |
| Finite-pressure mixture/transport/flexibility claims | Prohibited | Not calculated in this study. |
| Universal force-field independence | Prohibited | Only a bounded WCA linker-side perturbation was validated. |
""".strip() + "\n"

README_MD = r"""
# MOF CO2/N2 selectivity mechanisms — publication workspace

This repository contains the curated publication state, not the full exploratory simulation history.

## Scientific story

Repulsion establishes an N2-favoring baseline across the 238-MOF population. Restoring Lennard-Jones attraction creates CO2 selectivity in every framework. Electrostatics then provides a further framework-dependent enhancement that is weakly associated with attraction and poorly encoded by conventional pore geometry. Geometry-matched M029/M085 provides the strongest microscopic counterexample, with a 1.860-unit electrostatic-selectivity contrast that is overwhelmingly CO2-driven and accompanied by a markedly different probability-supported spatial response.

See `docs/STORY.md` for the claim hierarchy and figure logic.

## Repository structure

- `manuscript/` — clean main text, SI, bibliography
- `figures/main/` — five main publication figures
- `figures/si/` — supporting figures
- `data/authoritative/` — frozen 238-MOF table
- `data/source_data/` — figure/statistical source data
- `data/spatial_fields/` — retained final M029/M085 spatial fields only
- `analysis/final/` — scripts used by the final evidence path
- `reproducibility/` — protocols, manifests, hashes
- `docs/` — claims, decisions, exclusion ledger, PI additions, changelog

## Scope

No new simulation is required by the current claim set. New analyses should be added only if they resolve a concrete PI/reviewer request that cannot be answered from the frozen evidence.
""".strip() + "\n"

def write_publication_files():
    (REPO / "manuscript" / "main.tex").write_text(MAIN_TEX, encoding="utf-8")
    (REPO / "manuscript" / "si.tex").write_text(SI_TEX, encoding="utf-8")
    (REPO / "docs" / "STORY.md").write_text(STORY_MD, encoding="utf-8")
    (REPO / "docs" / "EXCLUDED_ANALYSES.md").write_text(EXCLUDED_MD, encoding="utf-8")
    (REPO / "README.md").write_text(README_MD, encoding="utf-8")

# ---------------------------------------------------------------------------
# Repository cleanup
# ---------------------------------------------------------------------------

def clean_redundant_artifacts():
    # PDF is the publication artifact; PNG duplicates are not versioned here.
    for p in (REPO / "figures").rglob("*.png"):
        safe_unlink(p)

    # Remove stale main figures that are superseded by clean rebuilt versions.
    safe_unlink(REPO / "figures" / "main" / "fig3_geometry_mechanism.pdf")
    safe_unlink(REPO / "figures" / "main" / "fig4_matched_pairs.pdf")
    safe_unlink(REPO / "figures" / "main" / "figure_electrostatic_pair_supported_xy.pdf")

    # Rename the retained per-gas SI figure to its final SI identity.
    old = REPO / "figures" / "si" / "figure_electrostatic_gas_increment_supported_xy.pdf"
    new = REPO / "figures" / "si" / "figS1_electrostatic_per_gas.pdf"
    if old.exists():
        if new.exists():
            new.unlink()
        old.rename(new)
        print(f"[renamed] {old.relative_to(REPO)} -> {new.relative_to(REPO)}")

    # Checkpoints are restart artifacts, not final source data. Original project
    # remains untouched and the accepted final field.npz files stay versioned.
    for p in FIELDS.rglob("checkpoint.npz"):
        safe_unlink(p)
    for p in FIELDS.rglob("field.csv"):
        safe_unlink(p)

    # Drop field directories that no longer contain an accepted final field.
    for mof_dir in [FIELDS / "M029", FIELDS / "M085"]:
        if not mof_dir.exists():
            continue
        for d in list(mof_dir.iterdir()):
            if d.is_dir() and not (d / "field.npz").exists():
                safe_rmtree(d)

def update_changelog():
    path = REPO / "docs" / "CHANGELOG.md"
    old = path.read_text(encoding="utf-8") if path.exists() else "# Changelog\n"
    entry = """
## Publication cleanup v1
- Replaced checkpoint prose with a single claim-driven manuscript narrative.
- Rebuilt Figure 3 from the authoritative leakage-free descriptor metrics.
- Rebuilt Figure 4 as the matched-pair trajectory plus retained electrostatic spatial evidence.
- Moved per-gas electrostatic maps and WCA force-field bridge to SI.
- Removed stale grouped-holdout Figure 3, redundant PNGs, spatial restart checkpoints, and non-evidentiary internal failure narratives from publication files.
- Added `STORY.md` and `EXCLUDED_ANALYSES.md`.
"""
    if "## Publication cleanup v1" not in old:
        path.write_text(old.rstrip() + "\n\n" + entry.strip() + "\n", encoding="utf-8")

def try_compile():
    report = []
    latexmk = shutil.which("latexmk")
    pdflatex = shutil.which("pdflatex")

    if latexmk:
        for tex in ["main.tex", "si.tex"]:
            p = subprocess.run(
                [latexmk, "-pdf", "-interaction=nonstopmode", tex],
                cwd=REPO / "manuscript",
                text=True,
                encoding="utf-8",
                errors="replace",
                capture_output=True,
            )
            pdf = (REPO / "manuscript" / tex.replace(".tex", ".pdf"))
            if p.returncode == 0:
                report.append(f"{tex}: PASS")
            elif pdf.exists():
                report.append(
                    f"{tex}: PDF produced, but latexmk returned {p.returncode}; "
                    "inspect bibliography/build warnings (Overleaf can perform the final bibliography pass)."
                )
            else:
                report.append(f"{tex}: FAIL, latexmk return code {p.returncode}")
    elif pdflatex:
        for tex in ["main.tex", "si.tex"]:
            p = subprocess.run(
                [pdflatex, "-interaction=nonstopmode", tex],
                cwd=REPO / "manuscript",
                text=True,
                encoding="utf-8",
                errors="replace",
                capture_output=True,
            )
            pdf = (REPO / "manuscript" / tex.replace(".tex", ".pdf"))
            if p.returncode == 0:
                report.append(f"{tex}: PDF produced by single pdflatex pass")
            elif pdf.exists():
                report.append(f"{tex}: PDF produced with warnings in single pdflatex pass")
            else:
                report.append(f"{tex}: FAIL, pdflatex return code {p.returncode}")
    else:
        report.append("No LaTeX executable found on PATH; compile check skipped.")

    (REPO / "reproducibility" / "cleanup_compile_report.txt").write_text(
        "\n".join(report) + "\n", encoding="utf-8"
    )
    return report

def main():
    if not REPO.exists():
        raise SystemExit(f"Publication repository not found:\n{REPO}")

    # Require a baseline commit. This is the safety boundary.
    head = git("rev-parse", "--verify", "HEAD", check=False)
    if head.returncode != 0:
        raise SystemExit(
            "\nSTOP: the scientific baseline has not been committed yet.\n\n"
            "Configure your real Git identity and commit the already staged baseline first:\n\n"
            f'  cd "{REPO}"\n'
            '  git config user.name "YOUR REAL NAME"\n'
            '  git config user.email "YOUR REAL EMAIL"\n'
            '  git commit -m "Freeze scientific checkpoint v4 before publication revisions"\n\n'
            "Then rerun this cleanup script."
        )

    dirty = git("status", "--porcelain").stdout.strip()
    if dirty:
        raise SystemExit(
            "\nSTOP: publication repository has uncommitted changes.\n"
            "Commit/stash them before running cleanup.\n\n"
            + dirty
        )

    # Tag the immutable scientific baseline if not already tagged.
    tags = set(git("tag").stdout.split())
    if "scientific-checkpoint-v4" not in tags:
        git("tag", "scientific-checkpoint-v4")
        print("[git] tagged baseline: scientific-checkpoint-v4")

    branches = set(
        line.strip().lstrip("* ").strip()
        for line in git("branch", "--format=%(refname:short)").stdout.splitlines()
    )
    current = git("branch", "--show-current").stdout.strip()

    if current != "publication-cleanup":
        if "publication-cleanup" in branches:
            git("switch", "publication-cleanup")
        else:
            git("switch", "-c", "publication-cleanup")
        print("[git] switched to publication-cleanup")

    auth = pd.read_csv(require(REPO / "data" / "authoritative" / "path_b_final_238_authoritative.csv"))

    print("\n[1/6] Cleaning redundant repository artifacts...")
    clean_redundant_artifacts()

    print("[2/6] Recomputing final spatial summary from retained fields...")
    spatial_summary, _ = build_spatial_summary(auth)

    print("[3/6] Rebuilding publication figures...")
    make_fig3()
    make_fig4(auth, spatial_summary)
    make_forcefield_si()

    print("[4/6] Replacing manuscript/SI with clean claim-driven versions...")
    write_publication_files()
    update_changelog()

    print("[5/6] Attempting LaTeX build check...")
    compile_report = try_compile()
    for x in compile_report:
        print("  " + x)

    print("[6/6] Staging cleanup...")
    git("add", "-A")

    name = git("config", "--get", "user.name", check=False).stdout.strip()
    email = git("config", "--get", "user.email", check=False).stdout.strip()

    if name and email:
        c = git(
            "commit",
            "-m",
            "Clean publication narrative and rebuild core figures",
            check=False,
        )
        if c.returncode == 0:
            print("[git] cleanup commit created")
        else:
            print("[git] changes staged but commit failed:")
            print(c.stderr)
    else:
        print("[git] changes staged; Git identity missing, so cleanup was not committed.")

    print("\n" + "=" * 88)
    print("PUBLICATION CLEANUP COMPLETE")
    print("=" * 88)
    print(f"Repository: {REPO}")
    print(f"Branch:     {git('branch', '--show-current').stdout.strip()}")
    print()
    print("Core story:")
    print("  repulsive N2-favoring baseline")
    print("        -> attraction creates CO2 preference")
    print("        -> electrostatics differentiates frameworks")
    print("        -> geometry misses most practically useful electrostatic information")
    print("        -> M029/M085 provides the spatially resolved counterexample")
    print("        -> interaction stages reorder screening priorities")
    print()
    print("Inspect next:")
    print("  manuscript/main.tex")
    print("  manuscript/si.tex")
    print("  docs/STORY.md")
    print("  docs/EXCLUDED_ANALYSES.md")
    print("  figures/main/fig3_geometry_information.pdf")
    print("  figures/main/fig4_matched_electrostatic.pdf")
    print()
    print("Next scientific step after inspection: literature/novelty pass, then PI additions.")
    print()
    print("Git status:")
    print(git("status", "--short").stdout or "  clean")
    print("Recent log:")
    print(git("log", "--oneline", "-4").stdout)

if __name__ == "__main__":
    main()
