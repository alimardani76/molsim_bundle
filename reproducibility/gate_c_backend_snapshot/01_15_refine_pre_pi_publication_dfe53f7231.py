from __future__ import annotations

import json
import math
import re
import shutil
import subprocess
import sys
from pathlib import Path

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

PROJECT_ROOT = Path(
    r"C:\Users\KaraPardazesh\Desktop\bundle\v1\mof_selectivity_mechanisms"
)
REPO = PROJECT_ROOT / "mof_selectivity_publication"
MANUSCRIPT = REPO / "manuscript"
FIG_MAIN = REPO / "figures" / "main"
FIG_SI = REPO / "figures" / "si"
SRC = REPO / "data" / "source_data"
AUTH = REPO / "data" / "authoritative" / "path_b_final_238_authoritative.csv"
FIELDS = REPO / "data" / "spatial_fields"
REPRO = REPO / "reproducibility"

TARGET_BRANCH = "pre-pi-publication"
COMMIT_MESSAGE = "Refine novelty framing and population-wide publication figures"

# The 9 descriptors used in the original nearest-neighbor pair-selection analysis.
PAIR_GEOMETRY_FEATURES = [
    "geom_df",
    "geom_di",
    "geom_dif",
    "geom_density",
    "geom_asa",
    "geom_avaf",
    "geom_navaf",
    "geom_poavaf",
    "geom_npoavaf",
]

IDS = {
    "M089": "DB0-m15_o11_o75_f0_fsc",
    "M192": "DB0-m15_o14_o77_f0_fsc",
    "M029": "DB0-m15_o100_o150_f0_fsc",
    "M085": "DB0-m15_o11_o30_jea",
}


def run(cmd, cwd=REPO, check=False):
    p = subprocess.run(
        [str(x) for x in cmd],
        cwd=cwd,
        text=True,
        encoding="utf-8",
        errors="replace",
        capture_output=True,
    )
    if check and p.returncode != 0:
        raise RuntimeError(
            f"Command failed ({p.returncode}): {' '.join(map(str, cmd))}\n"
            f"STDOUT:\n{p.stdout}\nSTDERR:\n{p.stderr}"
        )
    return p


def git(*args, check=False):
    return run(["git", *args], cwd=REPO, check=check)


def require(path: Path) -> Path:
    if not path.exists():
        raise FileNotFoundError(path)
    return path


def copy_if_exists(src: Path, dst: Path):
    if src.exists():
        dst.parent.mkdir(parents=True, exist_ok=True)
        shutil.copy2(src, dst)
        print(f"[copied reproducibility input] {src}")
        return True
    print(f"[optional reproducibility input not found] {src}")
    return False


# ============================================================================
# Spatial helpers
# ============================================================================

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
        raise FileNotFoundError(f"Accepted field missing for {mof}/{gas}/{mode}")

    d = np.load(folder / "field.npz")
    shape = tuple(int(x) for x in d["shape"])
    return {
        "folder": folder,
        "shape": shape,
        "logw": np.asarray(d["log_weight"], float).reshape(shape),
        "rel": np.asarray(d["relative_density"], float).reshape(shape),
        "cell": np.asarray(d["unit_cell_vectors_A"], float),
    }


def assert_same_grid(*fs):
    for f in fs[1:]:
        if f["shape"] != fs[0]["shape"]:
            raise RuntimeError("Spatial field grid-shape mismatch")
        if not np.allclose(f["cell"], fs[0]["cell"], atol=1e-12, rtol=0):
            raise RuntimeError("Spatial field cell mismatch")


def support_mask(rel, coverage=0.995):
    x = np.asarray(rel, float).ravel()
    order = np.argsort(x)[::-1]
    total = x.sum()
    if not np.isfinite(total) or total <= 0:
        raise RuntimeError("Invalid relative-density field")
    cum = np.cumsum(x[order]) / total
    k = int(np.searchsorted(cum, coverage, side="left"))
    m = np.zeros(x.size, dtype=bool)
    m[order[: k + 1]] = True
    return m.reshape(rel.shape)


def support_union(*rels):
    m = np.zeros_like(rels[0], dtype=bool)
    for rel in rels:
        m |= support_mask(rel)
    return m


def weighted_z_projection(values, weights, support):
    values = np.asarray(values, float)
    weights = np.asarray(weights, float)
    w = np.where(support, weights, 0.0)
    num = np.sum(np.where(support, values, 0.0) * w, axis=2)
    den = np.sum(w, axis=2)
    out = np.full_like(num, np.nan, dtype=float)
    ok = den > 0
    out[ok] = num[ok] / den[ok]
    return out


def xy_extent(cell):
    lengths = np.linalg.norm(cell, axis=1)
    return [0.0, float(lengths[0]), 0.0, float(lengths[1])]


# ============================================================================
# Derived source data
# ============================================================================

def build_pair_geometry_context(d: pd.DataFrame) -> pd.DataFrame:
    missing = [c for c in PAIR_GEOMETRY_FEATURES if c not in d.columns]
    if missing:
        raise RuntimeError(f"Pair geometry columns missing: {missing}")

    X = d[PAIR_GEOMETRY_FEATURES].astype(float).copy()
    X = X.fillna(X.median())
    # Population standardization, matching sklearn StandardScaler(ddof=0).
    Z = (X - X.mean(axis=0)) / X.std(axis=0, ddof=0)
    arr = Z.to_numpy()

    rows = []
    for a, b, target in [
        ("M089", "M192", "attraction"),
        ("M029", "M085", "electrostatic"),
    ]:
        ia = int(d.index[d["mof_id"].eq(IDS[a])][0])
        ib = int(d.index[d["mof_id"].eq(IDS[b])][0])

        dist = float(np.linalg.norm(arr[ia] - arr[ib]))

        da = np.linalg.norm(arr - arr[ia], axis=1)
        da[ia] = np.inf
        order_a = np.argsort(da)
        rank_a_to_b = int(np.where(order_a == ib)[0][0]) + 1

        db = np.linalg.norm(arr - arr[ib], axis=1)
        db[ib] = np.inf
        order_b = np.argsort(db)
        rank_b_to_a = int(np.where(order_b == ia)[0][0]) + 1

        rows.append({
            "pair_type": target,
            "mof_a": a,
            "mof_b": b,
            "mof_id_a": IDS[a],
            "mof_id_b": IDS[b],
            "geometry_distance_standardized_9d": dist,
            "neighbor_rank_b_from_a": rank_a_to_b,
            "neighbor_rank_a_from_b": rank_b_to_a,
            "geometry_features": ";".join(PAIR_GEOMETRY_FEATURES),
        })

    out = pd.DataFrame(rows)
    out.to_csv(SRC / "matched_pair_geometry_context.csv", index=False)
    return out


def build_rank_source(d: pd.DataFrame) -> tuple[pd.DataFrame, dict]:
    x = d[["mof_id", "lnS_wca", "lnS_lj_only", "lnS_full"]].copy()
    x["rank_wca"] = x["lnS_wca"].rank(method="min", ascending=False).astype(int)
    x["rank_lj"] = x["lnS_lj_only"].rank(method="min", ascending=False).astype(int)
    x["rank_full"] = x["lnS_full"].rank(method="min", ascending=False).astype(int)
    x["attraction_promotion"] = x["rank_wca"] - x["rank_lj"]
    x["electrostatic_promotion"] = x["rank_lj"] - x["rank_full"]

    stats = {
        "rho_wca_lj": float(x["rank_wca"].corr(x["rank_lj"], method="spearman")),
        "rho_lj_full": float(x["rank_lj"].corr(x["rank_full"], method="spearman")),
        "median_abs_attr_shift": float(x["attraction_promotion"].abs().median()),
        "median_abs_elec_shift": float(x["electrostatic_promotion"].abs().median()),
        "max_attr_promotion": int(x["attraction_promotion"].max()),
        "max_attr_demotion": int((-x["attraction_promotion"]).max()),
        "max_elec_promotion": int(x["electrostatic_promotion"].max()),
        "max_elec_demotion": int((-x["electrostatic_promotion"]).max()),
        "n_abs_attr_shift_ge_50": int((x["attraction_promotion"].abs() >= 50).sum()),
        "n_abs_elec_shift_ge_50": int((x["electrostatic_promotion"].abs() >= 50).sum()),
    }
    x.to_csv(SRC / "fig5_rank_transition_source.csv", index=False)
    (SRC / "fig5_rank_transition_stats.json").write_text(
        json.dumps(stats, indent=2), encoding="utf-8"
    )
    return x, stats


# ============================================================================
# Figure rebuilding
# ============================================================================

def make_figure2(d: pd.DataFrame):
    fig, ax = plt.subplots(figsize=(7.2, 5.6), constrained_layout=True)

    ax.scatter(
        d["delta_lnS_attr"],
        d["delta_lnS_elec"],
        s=26,
        alpha=0.42,
    )

    pair_specs = [
        ("M089", "M192", "attraction-matched pair"),
        ("M029", "M085", "electrostatic-matched pair"),
    ]
    for a, b, label in pair_specs:
        ra = d.loc[d["mof_id"].eq(IDS[a])].iloc[0]
        rb = d.loc[d["mof_id"].eq(IDS[b])].iloc[0]
        ax.plot(
            [ra["delta_lnS_attr"], rb["delta_lnS_attr"]],
            [ra["delta_lnS_elec"], rb["delta_lnS_elec"]],
            marker="o",
            linewidth=1.5,
            label=label,
        )
        ax.annotate(
            a,
            (ra["delta_lnS_attr"], ra["delta_lnS_elec"]),
            xytext=(5, 5),
            textcoords="offset points",
            fontsize=8,
        )
        ax.annotate(
            b,
            (rb["delta_lnS_attr"], rb["delta_lnS_elec"]),
            xytext=(5, 5),
            textcoords="offset points",
            fontsize=8,
        )

    rho = float(
        d["delta_lnS_attr"].corr(d["delta_lnS_elec"], method="spearman")
    )
    ax.text(
        0.03,
        0.96,
        rf"Spearman $\rho={rho:.3f}$",
        transform=ax.transAxes,
        va="top",
    )
    ax.set_xlabel(r"Attractive increment, $\Delta\ln S_{\rm attr}$")
    ax.set_ylabel(r"Electrostatic increment, $\Delta\ln S_{\rm elec}$")
    ax.set_title("Attractive and electrostatic enhancements are weakly associated")
    ax.legend(frameon=False)

    fig.savefig(FIG_MAIN / "fig2_mechanism_map.pdf", bbox_inches="tight")
    plt.close(fig)

    d[
        ["mof_id", "lnS_wca", "delta_lnS_attr", "delta_lnS_elec", "lnS_full"]
    ].to_csv(SRC / "fig2_mechanism_map_source.csv", index=False)


def make_figure4(d: pd.DataFrame, pair_context: pd.DataFrame):
    r29 = d.loc[d["mof_id"].eq(IDS["M029"])].iloc[0]
    r85 = d.loc[d["mof_id"].eq(IDS["M085"])].iloc[0]

    # Spatial fields.
    c29_lj = load_field("M029", "CO2", "lj_only")
    n29_lj = load_field("M029", "N2", "lj_only")
    c29_fu = load_field("M029", "CO2", "full")
    n29_fu = load_field("M029", "N2", "full")
    c85_lj = load_field("M085", "CO2", "lj_only")
    n85_lj = load_field("M085", "N2", "lj_only")
    c85_fu = load_field("M085", "CO2", "full")
    n85_fu = load_field("M085", "N2", "full")

    assert_same_grid(c29_lj, n29_lj, c29_fu, n29_fu)
    assert_same_grid(c85_lj, n85_lj, c85_fu, n85_fu)

    spatial = {}
    all_proj = []
    for mof, clj, nlj, cfu, nfu in [
        ("M029", c29_lj, n29_lj, c29_fu, n29_fu),
        ("M085", c85_lj, n85_lj, c85_fu, n85_fu),
    ]:
        support = support_union(clj["rel"], nlj["rel"], cfu["rel"], nfu["rel"])
        weights = clj["rel"] + nlj["rel"] + cfu["rel"] + nfu["rel"]
        delta = (cfu["logw"] - nfu["logw"]) - (clj["logw"] - nlj["logw"])
        proj = weighted_z_projection(delta, weights, support)
        spatial[mof] = (proj, xy_extent(clj["cell"]))
        all_proj.append(proj[np.isfinite(proj)])

    lim = float(np.quantile(np.abs(np.concatenate(all_proj)), 0.99))
    if not np.isfinite(lim) or lim <= 0:
        lim = 1.0

    # Global per-gas electrostatic increments.
    pergas = pd.DataFrame([
        {
            "mof": "M029",
            "CO2": float(r29["lnkh_mean_full_co2"] - r29["lnkh_mean_lj_only_co2"]),
            "N2": float(r29["lnkh_mean_full_n2"] - r29["lnkh_mean_lj_only_n2"]),
        },
        {
            "mof": "M085",
            "CO2": float(r85["lnkh_mean_full_co2"] - r85["lnkh_mean_lj_only_co2"]),
            "N2": float(r85["lnkh_mean_full_n2"] - r85["lnkh_mean_lj_only_n2"]),
        },
    ])
    pergas.to_csv(SRC / "fig4_per_gas_global_increment_source.csv", index=False)

    fig = plt.figure(figsize=(10.5, 7.5), constrained_layout=True)
    gs = fig.add_gridspec(2, 2, height_ratios=[0.78, 1.18])
    ax0 = fig.add_subplot(gs[0, 0])
    ax1 = fig.add_subplot(gs[0, 1])
    ax2 = fig.add_subplot(gs[1, 0])
    ax3 = fig.add_subplot(gs[1, 1])

    stages = np.arange(3)
    labels = ["WCA", "+ attraction", "+ electrostatics"]
    for alias, r in [("M029", r29), ("M085", r85)]:
        ax0.plot(
            stages,
            [r["lnS_wca"], r["lnS_lj_only"], r["lnS_full"]],
            marker="o",
            label=alias,
        )
    ax0.axhline(0, linewidth=0.8)
    ax0.set_xticks(stages, labels)
    ax0.set_ylabel(r"$\ln S_{\mathrm{CO_2/N_2}}$")
    ax0.set_title("Matched frameworks diverge when electrostatics is restored")
    ax0.legend(frameon=False)

    x = np.arange(2)
    width = 0.34
    ax1.bar(
        x - width / 2,
        pergas["CO2"].to_numpy(),
        width=width,
        label=r"CO$_2$",
    )
    ax1.bar(
        x + width / 2,
        pergas["N2"].to_numpy(),
        width=width,
        label=r"N$_2$",
    )
    ax1.set_xticks(x, pergas["mof"])
    ax1.set_ylabel(r"$\Delta\ln K_H$ (full $-$ LJ)")
    ax1.set_title("The between-framework contrast is CO$_2$-driven")
    ax1.legend(frameon=False)

    images = []
    for ax, mof in [(ax2, "M029"), (ax3, "M085")]:
        arr, extent = spatial[mof]
        im = ax.imshow(
            np.ma.masked_invalid(arr).T,
            origin="lower",
            aspect="auto",
            extent=extent,
            vmin=-lim,
            vmax=lim,
        )
        images.append(im)
        rr = r29 if mof == "M029" else r85
        ax.set_xlabel("x (Å)")
        ax.set_ylabel("y (Å)")
        ax.set_title(
            rf"{mof}: local full$-$LJ selectivity response; "
            rf"$\Delta\ln S_{{elec}}={rr['delta_lnS_elec']:.3f}$"
        )

    cb = fig.colorbar(images[-1], ax=[ax2, ax3], shrink=0.84)
    cb.set_label(r"Projected local $\delta_{\rm elec}$")

    ctx = pair_context.loc[pair_context["pair_type"].eq("electrostatic")].iloc[0]
    fig.suptitle(
        "A geometry-matched pair exposes framework-specific electrostatic recognition\n"
        f"M029 and M085 are mutual nearest neighbors in the 9-descriptor geometry space "
        f"(distance {ctx['geometry_distance_standardized_9d']:.3f})",
        fontsize=12,
    )

    fig.savefig(FIG_MAIN / "fig4_matched_electrostatic.pdf", bbox_inches="tight")
    plt.close(fig)


def make_figure5(rank_df: pd.DataFrame, stats: dict):
    fig, axes = plt.subplots(1, 2, figsize=(10.2, 4.5), constrained_layout=True)

    specs = [
        (
            axes[0],
            "rank_wca",
            "rank_lj",
            "WCA rank",
            "LJ-only rank",
            stats["rho_wca_lj"],
            stats["median_abs_attr_shift"],
            stats["max_attr_promotion"],
            stats["max_attr_demotion"],
            "Attraction globally reorganizes the ranking",
        ),
        (
            axes[1],
            "rank_lj",
            "rank_full",
            "LJ-only rank",
            "Full rank",
            stats["rho_lj_full"],
            stats["median_abs_elec_shift"],
            stats["max_elec_promotion"],
            stats["max_elec_demotion"],
            "Electrostatics preserves global order but creates outliers",
        ),
    ]

    n = len(rank_df)
    for ax, xc, yc, xl, yl, rho, med, pmax, dmax, title in specs:
        ax.scatter(rank_df[xc], rank_df[yc], s=22, alpha=0.55)
        ax.plot([1, n], [1, n], linestyle="--", linewidth=1)
        ax.set_xlim(0, n + 1)
        ax.set_ylim(0, n + 1)
        ax.set_xlabel(xl)
        ax.set_ylabel(yl)
        ax.set_aspect("equal", adjustable="box")
        ax.set_title(title, fontsize=10)
        ax.text(
            0.03,
            0.97,
            rf"$\rho={rho:.3f}$" + "\n"
            + rf"median $|\Delta rank|={med:.0f}$" + "\n"
            + f"max promotion/demotion = {pmax}/{dmax}",
            transform=ax.transAxes,
            va="top",
            fontsize=9,
        )

    fig.suptitle("Successive interaction stages have different effects on candidate ordering")
    fig.savefig(FIG_MAIN / "fig5_rank_trajectories.pdf", bbox_inches="tight")
    plt.close(fig)


# ============================================================================
# Publication text
# ============================================================================

MAIN_TEX = r"""
\documentclass[11pt,a4paper]{article}
\usepackage[utf8]{inputenc}
\usepackage[T1]{fontenc}
\usepackage[margin=2.05cm]{geometry}
\usepackage{microtype,booktabs,tabularx,array,amsmath,amssymb,mathtools,graphicx,float,caption,subcaption,hyperref,siunitx,bm}
\graphicspath{{../figures/main/}{../figures/si/}}
\hypersetup{colorlinks=true,linkcolor=blue,citecolor=teal,urlcolor=teal}

\title{How CO$_2$/N$_2$ Selectivity Emerges in Metal--Organic Frameworks:\\Repulsive Baselines, Attractive Reversal, and Electrostatic Recognition}
\author{Author list to be finalized}
\date{}

\begin{document}
\maketitle

\begin{abstract}
High-throughput screening ranks metal--organic frameworks (MOFs) by adsorption performance, but a final selectivity does not reveal which intermolecular interactions produced it. We construct a nested reference-potential decomposition of Henry-regime CO$_2$/N$_2$ selectivity across 238 MOFs, separating a repulsive Weeks--Chandler--Andersen (WCA) reference from the increments obtained by restoring Lennard--Jones attraction and then electrostatics. Repulsive accessibility alone favored N$_2$ in every framework (median $S_{\rm WCA}=0.816$). Restoring attraction reversed this preference throughout the panel, with a median $\Delta\ln S_{\rm attr}=2.537$ (12.65-fold), whereas electrostatics supplied a further positive but strongly framework-dependent increment (median $\Delta\ln S_{\rm elec}=0.309$, 1.36-fold). The attractive and electrostatic increments were only weakly associated ($\rho=-0.089$). Conventional pore geometry predicted the WCA and attractive coordinates strongly ($R^2=0.835$ and 0.874) and also predicted the final selectivity ($R^2=0.845$), yet provided essentially no reduction in mean absolute error for the electrostatic increment relative to a median predictor (0.146 versus 0.147). A geometry-matched pair makes this representation boundary explicit: M029 and M085 are mutual nearest neighbors in the geometry space but have electrostatic increments of 0.141 and 2.001. Global Henry coefficients and probability-supported spatial fields show that this contrast is overwhelmingly CO$_2$-driven and is accompanied by a markedly different local electrostatic response. The results show that accurate endpoint prediction can conceal missing mechanistic information and provide a controlled route for separating accessibility, attraction, and electrostatic recognition in MOF screening.
\end{abstract}

\section{Introduction}

Metal--organic frameworks (MOFs) provide a large design space for adsorption-based separations, and high-throughput molecular simulation has become a standard way to screen that space \cite{haldoupis2012charges,burner2023arcmof}. For post-combustion CO$_2$/N$_2$ separation, screening studies commonly reduce each framework to quantities such as uptake, Henry coefficient, working capacity, or selectivity. These endpoints are useful for ranking materials, but they do not state which intermolecular interactions made a framework selective.

The importance of electrostatics in CO$_2$/N$_2$ adsorption is already well established. Haldoupis, Nair, and Sholl showed that reliable charge assignment is central to large-scale screening of quadrupolar CO$_2$ and N$_2$ \cite{haldoupis2012charges}. In a later screening of 3,806 MOFs, Daglar and Keskin explicitly repeated simulations with framework electrostatics removed and found substantial changes in CO$_2$/N$_2$ adsorption performance and in the identity of top candidates \cite{daglar2018screening}. Such comparisons establish that electrostatics matters, but a full-versus-uncharged comparison still combines steric exclusion and dispersive attraction in the same non-electrostatic reference.

A parallel literature asks how much adsorption performance can be inferred from structure alone. Geometry-based descriptors can provide effective low-cost filters for CO$_2$ and N$_2$ adsorption \cite{fernandez2016geometry}, whereas large-scale machine-learning studies have shown that chemistry-aware descriptors can substantially improve low-pressure CO$_2$ working-capacity and CO$_2$/N$_2$ selectivity prediction beyond geometry-only models \cite{burner2020deeplearning}. These results establish a representation hierarchy for endpoint prediction. They do not, however, identify which physical contribution is captured or missed by a representation.

Here we connect these two questions. We construct a nested sequence of reference potentials that separates repulsive accessibility, Lennard--Jones attraction, and electrostatic recognition for the same framework population. The WCA construction provides the repulsive Lennard--Jones reference \cite{weeks1971wca}; related WCA splitting has previously been used to isolate repulsive contributions in MOF adsorption models \cite{campbell2018polar}. Restoring the full Lennard--Jones potential then introduces attraction, and restoring framework--guest electrostatics gives the complete model:
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

Across 238 MOFs, we ask which interaction stage first creates CO$_2$ preference, whether attractive and electrostatic enhancements track one another, and which of these coordinates is encoded by conventional pore geometry. We then use a geometry-matched pair and local insertion fields to test whether a large electrostatic increment corresponds to a distinct spatial response. The resulting picture is sequential rather than symmetric: repulsion sets an N$_2$-favoring baseline, attraction creates CO$_2$ selectivity throughout this panel, and electrostatics provides an additional framework-specific recognition channel that conventional pore geometry does not adequately resolve.

\section{Results}

\subsection{Repulsion sets the baseline and attraction creates CO$_2$ preference}

The common analysis set contains 238 MOFs evaluated in WCA, LJ-only, and full interaction states. The nested identity closes numerically to a maximum residual of $8.88\times10^{-16}$.

Under the WCA reference, all 238 frameworks are N$_2$-favoring. The median $\ln S_{\rm WCA}$ is $-0.204$, corresponding to a median $S_{\rm WCA}=0.816$. Repulsive accessibility therefore does not constitute a successful CO$_2$-selectivity route in this population.

Restoring attractive Lennard--Jones interactions changes the result qualitatively. $\Delta\ln S_{\rm attr}$ is positive for every framework, with a median of 2.537, corresponding to a 12.65-fold multiplicative contribution. The LJ-only selectivity is CO$_2$-favoring throughout the panel. Attraction is therefore the interaction stage that reverses the repulsive N$_2$ preference and establishes CO$_2$ selectivity under these models.

\begin{figure}[H]
\centering
\includegraphics[width=\linewidth]{fig1_decomposition.pdf}
\caption{\textbf{Nested reference-potential decomposition of CO$_2$/N$_2$ selectivity.} The WCA state retains repulsive accessibility, the LJ-only state restores attractive Lennard--Jones interactions, and the full state additionally restores electrostatics. The population distributions reveal an N$_2$-favoring repulsive baseline, a positive attractive increment throughout the panel, and a positive but more heterogeneous electrostatic increment.}
\label{fig:decomp}
\end{figure}

\subsection{Attraction and electrostatics are weakly associated across the population}

Restoring electrostatics provides a further positive enhancement in every framework. The median $\Delta\ln S_{\rm elec}$ is 0.309, equivalent to a 1.36-fold multiplicative contribution, while the largest increments approach 2 log units.

The attractive and electrostatic increments are only weakly associated across the population (Spearman $\rho=-0.089$; Figure~\ref{fig:map}). A framework's gain from attractive confinement therefore provides little information about the additional gain it receives from electrostatic interactions. This does not imply statistical independence; it shows that the two finite-difference coordinates encode substantially different framework-to-framework variation.

\begin{figure}[H]
\centering
\includegraphics[width=.78\linewidth]{fig2_mechanism_map.pdf}
\caption{\textbf{Attractive and electrostatic enhancements are weakly associated across 238 MOFs.} Each point is one framework. The frozen attraction- and electrostatic-matched pairs are highlighted only as case studies; no categorical mechanism boundaries are imposed on the population.}
\label{fig:map}
\end{figure}

\subsection{Geometry predicts the endpoint while missing the electrostatic increment}

We next ask how much information about each mechanistic coordinate is already available from conventional pre-simulation pore geometry. ExtraTrees models using 14 geometry descriptors give five-fold cross-validated $R^2=0.835$ for $\ln S_{\rm WCA}$ and $R^2=0.874$ for $\Delta\ln S_{\rm attr}$. Their mean absolute errors are 0.033 and 0.177, compared with median-predictor baselines of 0.093 and 0.549.

The electrostatic increment behaves differently. Geometry gives $R^2=0.188$, but its MAE of 0.146 is essentially identical to the baseline MAE of 0.147. Thus the descriptors retain some population-level association with $\Delta\ln S_{\rm elec}$, yet provide almost no practical gain in absolute prediction error. In contrast, the same geometry set predicts the final $\ln S_{\rm full}$ well ($R^2=0.845$; MAE 0.189 versus a 0.501 baseline). Accurate prediction of the endpoint therefore does not imply that the representation resolves every interaction contribution used to construct that endpoint.

\begin{figure}[H]
\centering
\includegraphics[width=\linewidth]{fig3_geometry_information.pdf}
\caption{\textbf{Geometry predicts the final endpoint while providing little practical gain for the electrostatic increment.} Left: five-fold cross-validated $R^2$ of geometry-only ExtraTrees models. Right: model MAE normalized by the corresponding median-predictor MAE; a value of one indicates no improvement over that baseline.}
\label{fig:geometry}
\end{figure}

\subsection{A geometry-matched pair exposes framework-specific electrostatic recognition}

M029 and M085 form the strongest local counterexample. In the standardized nine-descriptor geometry space used for matched-pair selection, they are mutual nearest neighbors (distance 0.258). Their repulsive log-selectivities are correspondingly close, $-0.201$ and $-0.198$. Yet their electrostatic increments are 0.141 and 2.001, a difference of 1.860, corresponding to a 6.42-fold contrast in the multiplicative electrostatic contribution.

The global Henry coefficients identify which gas drives this divergence. Restoring electrostatics changes $\ln K_H$ for CO$_2$ by 0.162 in M029 and 2.032 in M085, whereas the N$_2$ changes are only 0.021 and 0.032. The between-framework difference is therefore 1.871 log units for CO$_2$ but only 0.011 for N$_2$: the much larger electrostatic selectivity of M085 is overwhelmingly CO$_2$-driven.

Local insertion fields provide a spatial counterpart to this scalar result. Within the 99.5\%-probability-supported pore region, M029 remains close to locally neutral over much of the accessible volume, whereas M085 shows a substantially more positive and heterogeneous full-minus-LJ selectivity response. The probability-weighted median local increment is 0.042 for M029 and 0.534 for M085; 57.6\% and 79.1\% of the corresponding probability weight lies in regions with positive local electrostatic selectivity increments. These maps localize the electrostatic perturbation but do not identify a specific binding motif.

\begin{figure}[H]
\centering
\includegraphics[width=\linewidth]{fig4_matched_electrostatic.pdf}
\caption{\textbf{A geometry-matched pair exposes framework-specific electrostatic recognition.} M029 and M085 are mutual nearest neighbors in the geometry space used for pair selection. Top left: selectivity trajectories through WCA, LJ-only, and full interaction states. Top right: gas-specific global full-minus-LJ changes in $\ln K_H$, showing that the between-framework contrast is dominated by CO$_2$. Bottom: probability-supported local full-minus-LJ selectivity increments. Global increments are taken from the Henry calculations; spatial fields are used only to localize the response.}
\label{fig:matched}
\end{figure}

\subsection{Attraction and electrostatics reorder candidates in different ways}

Restoring attraction changes the candidate ordering globally. The Spearman correlation between WCA and LJ-only ranks is $-0.751$, and the median absolute rank displacement is 102 positions among 238 frameworks. The largest attraction-driven promotion is 234 positions and the largest demotion is 226.

Electrostatics has a different population-level effect. LJ-only and full ranks remain strongly correlated ($\rho=0.905$), with a median absolute displacement of only 9 positions. Nevertheless, specific frameworks are strongly affected: the largest electrostatic promotion is 175 positions and the largest demotion is 55. Attraction therefore reorganizes the ranking broadly, whereas electrostatics mostly preserves the global ordering while creating important framework-specific outliers.

\begin{figure}[H]
\centering
\includegraphics[width=\linewidth]{fig5_rank_trajectories.pdf}
\caption{\textbf{Successive interaction stages have different effects on candidate ordering.} Every point is one MOF. Restoring attraction broadly reorganizes the WCA ranking, whereas restoring electrostatics leaves the global LJ-only ordering much more intact but produces large outlier promotions. Dashed lines denote unchanged rank.}
\label{fig:ranks}
\end{figure}

\section{Discussion}

Previous high-throughput studies established that electrostatics is important for CO$_2$/N$_2$ adsorption and that turning framework charges off can change both Henry selectivities and the identity of top candidates \cite{haldoupis2012charges,daglar2018screening}. The present decomposition refines that conclusion. In this panel, electrostatics is not the interaction stage that first creates CO$_2$ preference: the WCA reference is N$_2$-favoring throughout, and restoring Lennard--Jones attraction alone reverses that baseline throughout. Electrostatics then acts as a second, strongly framework-dependent recognition coordinate.

The representation result provides a complementary distinction. Geometry-based screening has long been useful for low-pressure adsorption \cite{fernandez2016geometry}, while chemistry-aware descriptors improve CO$_2$ working-capacity and selectivity prediction in large MOF libraries \cite{burner2020deeplearning}. Our results locate one reason for that hierarchy. The same conventional pore geometry that predicts $\ln S_{\rm full}$ accurately captures the repulsive and attractive coordinates but barely improves absolute prediction error for $\Delta\ln S_{\rm elec}$. Endpoint accuracy can therefore hide a representation boundary: a model may predict the final number without resolving one of the physical increments that constructs it.

The M029/M085 pair turns that population result into a controlled local comparison. The two structures are mutual nearest neighbors in the geometry space used for pair selection and have almost identical WCA selectivities, yet their electrostatic responses differ by 1.860 log units. The contrast is almost entirely associated with CO$_2$, and the insertion fields show a correspondingly stronger and more heterogeneous local response in M085. The maps do not establish which atom type or chemical motif is responsible; that requires a validated chemistry-aware local representation beyond the scope of the present study.

The ranking analysis shows why retaining the interaction coordinates matters for screening. Attraction changes the ordering of the panel globally, consistent with its role in reversing the repulsive baseline. Electrostatics leaves the overall LJ-only order more recognizable but can still rescue individual frameworks by more than 100 positions. A final rank therefore answers which material performs well under the complete model, whereas the nested coordinates reveal how it became competitive.

The conclusions are deliberately bounded to Henry-regime thermodynamic recognition at 298~K in rigid frameworks under the selected molecular models. A preselected 15-MOF DREIDING-linker/UFF-Co bridge produced only small changes in the WCA excess chemical potentials and essentially unchanged WCA rankings (Supporting Information), supporting the repulsive-baseline interpretation within that perturbation. Equivalent force-field robustness has not been established for the LJ-only or electrostatic increments. The study does not address finite-pressure mixtures, transport, framework flexibility, or transfer to other gas pairs.

\section{Methods}

\subsection{Framework panel and molecular models}

The starting panel contained 240 ARC--MOF structures \cite{burner2023arcmof}; 238 had the complete WCA, LJ-only, and full CO$_2$/N$_2$ quantities required for the common analysis. Frameworks were treated as rigid at 298~K. Full and LJ-only Henry calculations used RASPA3 version 3.0.29 \cite{ran2024raspa3,dubbeldam2016raspa}. The framework Lennard--Jones model used element-level UFF parameters \cite{rappe1992uff} with Lorentz--Berthelot mixing, a 12~\AA{} cutoff, shifted truncation, and no tail correction. Framework charges were taken from the ARC--MOF structure records. Full calculations used Ewald electrostatics; LJ-only calculations set framework and guest pseudoatom charges to zero while preserving the Lennard--Jones definitions.

The full calculation used 50,000 Widom cycles per seed and the LJ-only calculation used 10,000 cycles per seed, with three matched default seeds and adaptive extension of unresolved cases. The exact protocol, guest geometries, and retained RASPA force-field definitions are archived with the publication source data.

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
The two increments are
\begin{align}
\Delta\ln S_{\rm attr} &= \ln S_{\rm LJ}-\ln S_{\rm WCA},\\
\Delta\ln S_{\rm elec} &= \ln S_{\rm full}-\ln S_{\rm LJ}.
\end{align}

For the WCA state, each framework--guest Lennard--Jones pair potential was truncated at
$r_{c,ij}=2^{1/6}\sigma_{ij}$ and shifted to zero at the cutoff. Molecular Widom insertion \cite{widom1963} was used to evaluate the repulsive reference. Independent equal-length WCA seeds were pooled in Boltzmann-factor space before conversion back to the excess chemical potential. Seed counts, precision diagnostics, and exclusion sensitivity are reported in the Supporting Information.

\subsection{Geometry information analysis}

The geometry-only analysis used 14 conventional pre-simulation descriptors. Target-derived adsorption quantities, identifiers, filenames, cluster labels, and other variables that could encode the response were excluded before fitting. ExtraTrees regressors used 500 trees, minimum leaf size three, a maximum-feature fraction of 0.7, and five-fold shuffled cross-validation with random seed 20260918. Prediction was used as an information-content diagnostic rather than as the primary objective. Performance was evaluated with $R^2$, MAE, and the MAE of a median predictor.

\subsection{Matched-pair analysis}

Matched pairs were selected from the closest 2\% of all framework pairs in a standardized nine-descriptor geometry space and then ranked by contrast in the target mechanistic coordinate. The retained attraction pair is M089/M192 and the retained electrostatic pair is M029/M085. For the electrostatic pair, M029 and M085 are mutual nearest neighbors in that geometry space. Pair selection preceded the spatial interpretation.

\subsection{Spatial insertion fields}

Spatial fields were computed on a fixed 0.75~\AA{} grid using the validated fixed-position RASPA/RASPALIB evaluator and the same rigid guest geometries, framework structures, Lennard--Jones parameters, charge assignments, mixing rules, and temperature as the corresponding Henry calculations. At each position $\mathbf r$,
\begin{equation}
w(\mathbf r)=
\left\langle
\exp[-\beta U(\mathbf r,\Omega)]
\right\rangle_{\Omega}.
\end{equation}
LJ-only fields used 2,048 orientations per voxel. Full fields used 16,384 or 32,768 orientations; targeted nested 16k$\rightarrow$32k$\rightarrow$64k checks at probability-important voxels changed the tested values by at most approximately 0.04\%.

For the electrostatic pair,
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
Visualization and local summaries were restricted to the union of the minimal voxel sets carrying 99.5\% of the insertion probability of the relevant fields. Local quantities were weighted by normalized insertion probability, preventing inaccessible framework-overlap voxels from dominating logarithmic averages. The spatial fields are localization diagnostics; the global Henry calculations remain the estimator for integrated selectivity increments.

\subsection{Robustness analyses}

The principal signs and the weak attractive--electrostatic association were recomputed after excluding the four WCA cases with $\mathrm{SE}(\ln S_{\rm WCA})>0.10$ and after removing the upper 1\% of each response separately. A preselected 15-MOF bridge additionally replaced supported linker/main-group UFF parameters with DREIDING while retaining UFF for cobalt; all other WCA definitions were preserved. This bridge was used only to test sensitivity of the repulsive-accessibility coordinate to that specified linker-side Lennard--Jones perturbation.

\section{Conclusions}

Across the 238-MOF analysis set, repulsive accessibility alone favors N$_2$, Lennard--Jones attraction reverses this baseline toward CO$_2$ throughout the panel, and electrostatics supplies a further framework-dependent enhancement. The attractive and electrostatic increments are only weakly associated.

Conventional pore geometry captures the repulsive and attractive coordinates strongly and predicts the final selectivity well, yet provides essentially no improvement in absolute electrostatic prediction error over a median baseline. The M029/M085 pair shows the consequence directly: two mutual nearest neighbors in geometry space develop a 1.860-unit difference in electrostatic selectivity, driven overwhelmingly by CO$_2$ and accompanied by markedly different local spatial responses.

Total selectivity is therefore a useful endpoint but an incomplete mechanistic descriptor. The nested reference potentials preserve the route by which a framework becomes selective and reveal when an apparently adequate structural representation predicts the endpoint without resolving the underlying interaction coordinate.

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

\title{Supporting Information\\How CO$_2$/N$_2$ Selectivity Emerges in Metal--Organic Frameworks}
\author{Author list to be finalized}
\date{}

\begin{document}
\maketitle
\tableofcontents

\section{Analysis population and simulation definitions}

The starting ARC--MOF panel contained 240 structures \cite{burner2023arcmof}. The common analysis population contains 238 MOFs with the complete WCA, LJ-only, and full CO$_2$/N$_2$ quantities required by the nested comparison.

Full and LJ-only Henry calculations used RASPA3 version 3.0.29 \cite{ran2024raspa3,dubbeldam2016raspa} at 298~K with rigid frameworks. The framework Lennard--Jones parameters were taken from UFF \cite{rappe1992uff}; the cutoff was 12~\AA{}, Lorentz--Berthelot mixing was used, interactions were shifted at the cutoff, and tail corrections were disabled. Full calculations used Ewald electrostatics; LJ-only calculations set framework and guest pseudoatom charges to zero while retaining the same Lennard--Jones definitions.

The full calculation used 50,000 Widom cycles per seed and LJ-only used 10,000 cycles per seed. Three matched seeds were used by default, with adaptive extension for unresolved cases. Exact protocol and force-field files are retained in the publication repository.

The WCA reference retained the repulsive branch of each framework--guest Lennard--Jones pair \cite{weeks1971wca},
\begin{equation}
r_{c,ij}=2^{1/6}\sigma_{ij},
\end{equation}
with the potential shifted to zero at the cutoff.

\section{WCA sampling, pooling, and precision}

The final WCA archive contains 2,052 jobs. Independent equal-length seeds were pooled through
\begin{equation}
B_k=\exp[-\mu_k^{\rm ex}/(RT)],
\end{equation}
followed by averaging in Boltzmann-factor space and conversion through $\ln\bar B$. The population contains 178 MOFs with three seeds per gas, 38 with six, three with nine, and 19 with twelve.

Four frameworks retained $\mathrm{SE}(\ln S_{\rm WCA})>0.10$. Excluding them does not change the sign or interpretation of the population-level conclusions.

\section{Exclusion sensitivity}

\begin{table}[H]
\centering
\caption{Sensitivity of the mechanism summaries to predefined exclusions.}
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

In every listed population, the WCA-positive fraction remains zero while the attractive- and electrostatic-positive fractions remain one.

\section{Geometry-only models}

The final information-content analysis used 14 pre-simulation geometry descriptors. ExtraTrees regressors contained 500 trees, minimum leaf size three, and maximum-feature fraction 0.7. Five-fold shuffled cross-validation used random seed 20260918.

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

Target-derived adsorption quantities, identifiers, filenames, cluster labels, and other variables capable of encoding the target were excluded before fitting. The complete descriptor manifest and out-of-fold prediction table are provided as source data.

\section{Geometry-matched pairs}

Pairs were selected from the closest 2\% of all pairs in the standardized nine-descriptor geometry space and ranked by contrast in the target mechanism coordinate.

\begin{table}[H]
\centering
\caption{Geometry-matched case studies.}
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

M089/M192 have standardized geometry distance 0.199; M192 is the nearest geometry neighbor of M089 and M089 is the third-nearest neighbor of M192. Their attractive increments differ by 0.658.

M029/M085 have standardized geometry distance 0.258 and are mutual nearest neighbors. Their electrostatic increments differ by 1.860. The attraction pair is retained as a scalar case study only; its available LJ spatial maps cannot isolate the attractive increment because WCA spatial fields were not calculated.

\section{Spatial mapping and orientation checks}

Spatial insertion fields used a fixed 0.75~\AA{} grid and the validated fixed-position RASPA/RASPALIB evaluator. LJ-only fields used 2,048 orientations per voxel. Full electrostatic fields used 16,384 or 32,768 orientations. Targeted nested 16k$\rightarrow$32k$\rightarrow$64k checks at probability-important voxels gave maximum observed 16k$\rightarrow$32k changes of approximately 0.04\% for M029/CO$_2$, 0.01\% for M029/N$_2$, and 0.002\% for M085/N$_2$; subsequent 32k$\rightarrow$64k changes were smaller.

Probability support was defined as the minimal voxel set containing 99.5\% of insertion probability mass. Local log-selectivity quantities were summarized only within the union of the relevant supports and were weighted by normalized insertion probability.

\begin{table}[H]
\centering
\caption{Probability-supported local electrostatic summaries. These are descriptive spatial statistics, not replacements for global Henry coefficients.}
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

\subsection{Gas-specific electrostatic response}

\begin{figure}[H]
\centering
\includegraphics[width=.95\linewidth]{figS1_electrostatic_per_gas.pdf}
\caption{\textbf{Gas-specific probability-supported full-minus-LJ increments for M029 and M085.} The local fields support the global result that the large between-framework electrostatic-selectivity contrast is dominated by CO$_2$.}
\label{fig:si_spatial_gas}
\end{figure}

From the global Henry calculations, $\Delta\ln K_H^{\rm full-LJ}$ for CO$_2$ is 0.162 in M029 and 2.032 in M085. The corresponding N$_2$ values are 0.021 and 0.032. Thus the M085--M029 difference is 1.871 log units for CO$_2$ but only 0.011 for N$_2$.

\section{Population-wide rank changes}

\begin{table}[H]
\centering
\caption{Changes in candidate ordering after each interaction stage. Rank 1 denotes the highest selectivity.}
\begin{tabular}{lrrrr}
\toprule
Transition & Spearman $\rho$ & Median $|\Delta rank|$ & Max promotion & Max demotion\\
\midrule
WCA $\rightarrow$ LJ-only & -0.751 & 102 & 234 & 226\\
LJ-only $\rightarrow$ full & 0.905 & 9 & 175 & 55\\
\bottomrule
\end{tabular}
\end{table}

The attraction step therefore reorganizes the ranking throughout the population, whereas electrostatics mostly preserves the global ordering but creates a smaller number of large rank changes.

\section{WCA force-field bridge}

The 15-MOF bridge panel was selected before the alternative-force-field results were inspected. Supported linker and main-group framework atoms were reassigned from UFF to DREIDING while cobalt retained UFF. Guest models, charges, structures, and the WCA construction were unchanged.

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
\caption{\textbf{WCA excess chemical potentials under the baseline and tested alternative framework Lennard--Jones assignments.} The dashed line denotes equality. This supports robustness only to the bounded linker-side perturbation tested here.}
\label{fig:si_ff_bridge}
\end{figure}

\section{Source data and reproducibility}

The publication repository contains the 238-MOF mechanism table, source data for all final figures, the geometry-only prediction tables, retained M029/M085 spatial fields, WCA bridge results, simulation protocol, guest definitions, force-field definitions, final analysis scripts, and file hashes.

\bibliographystyle{unsrt}
\bibliography{references}
\end{document}
""".strip() + "\n"


NEW_BIB = r"""
@article{haldoupis2012charges,
author={Haldoupis, Emmanuel and Nair, Sankar and Sholl, David S.},
title={Finding MOFs for Highly Selective {CO2/N2} Adsorption Using Materials Screening Based on Efficient Assignment of Atomic Point Charges},
journal={Journal of the American Chemical Society},
year={2012},
volume={134},
number={9},
pages={4313--4323},
doi={10.1021/ja2108239}
}

@article{fernandez2016geometry,
author={Fernandez, Michael and Barnard, Amanda S.},
title={Geometrical Properties Can Predict {CO2} and {N2} Adsorption Performance of Metal--Organic Frameworks ({MOFs}) at Low Pressure},
journal={ACS Combinatorial Science},
year={2016},
volume={18},
number={5},
pages={243--252},
doi={10.1021/acscombsci.5b00188}
}

@article{daglar2018screening,
author={Daglar, Hilal and Keskin, Seda},
title={Computational Screening of Metal--Organic Frameworks for Membrane-Based {CO2/N2/H2O} Separations: Best Materials for Flue Gas Separation},
journal={The Journal of Physical Chemistry C},
year={2018},
volume={122},
number={30},
pages={17347--17357},
doi={10.1021/acs.jpcc.8b05416}
}

@article{campbell2018polar,
author={Campbell, Christopher and Gomes, Jos{\'e} R. B. and Fischer, Michael and Jorge, Miguel},
title={New Model for Predicting Adsorption of Polar Molecules in Metal--Organic Frameworks with Unsaturated Metal Sites},
journal={The Journal of Physical Chemistry Letters},
year={2018},
volume={9},
number={12},
pages={3544--3553},
doi={10.1021/acs.jpclett.8b00967}
}

@article{burner2020deeplearning,
author={Burner, Jake and Schwiedrzik, Ludwig and Krykunov, Mykhaylo and Luo, Jun and Boyd, Peter G. and Woo, Tom K.},
title={High-Performing Deep Learning Regression Models for Predicting Low-Pressure {CO2} Adsorption Properties of Metal--Organic Frameworks},
journal={The Journal of Physical Chemistry C},
year={2020},
volume={124},
number={51},
pages={27996--28005},
doi={10.1021/acs.jpcc.0c06334}
}
""".strip() + "\n"


NOVELTY_MD = r"""
# Novelty framing for the pre-PI manuscript

## Established prior art

1. Electrostatics and charge assignment matter for quadrupolar CO2/N2 adsorption in MOFs.
   - Haldoupis, Nair, Sholl, JACS 2012.
2. Large-scale with/without-charge screening already shows electrostatics can change CO2/N2 performance and candidate identity.
   - Daglar, Keskin, JPCC 2018.
3. Geometry can be useful for low-pressure adsorption screening.
   - Fernandez, Barnard, ACS Combinatorial Science 2016.
4. Chemistry-aware descriptors can outperform geometry-only descriptors for low-pressure CO2/selectivity prediction.
   - Burner et al., JPCC 2020.
5. WCA splitting itself is established statistical mechanics and has precedent in MOF adsorption model construction.
   - Weeks, Chandler, Andersen 1971; Campbell et al., JPCL 2018.

## Manuscript contribution

The manuscript should not claim novelty for any item above. Its contribution is the combination of:

- a population-wide nested WCA -> LJ-only -> full mechanism coordinate;
- a panel-wide repulsive-to-attractive sign reversal;
- weak association between attractive and electrostatic increments;
- a representation-boundary result showing that geometry predicts the endpoint while barely improving electrostatic absolute error;
- a mutual-nearest-neighbor M029/M085 counterexample with gas-specific and spatial attribution;
- population-wide rank-transition behavior showing attraction and electrostatics affect screening order differently.

## Strongest discussion sentence

"Endpoint accuracy can therefore hide a representation boundary: a model may predict the final number without resolving one of the physical increments that constructs it."

## Boundary

Use "in this panel", "under these models", and "Henry regime" consistently. Do not claim universal interaction hierarchy, universal force-field robustness, or finite-pressure separation performance.
""".strip() + "\n"


def append_bib_entries():
    path = MANUSCRIPT / "references.bib"
    old = path.read_text(encoding="utf-8")
    additions = []
    for key_block in NEW_BIB.strip().split("\n\n@"):
        block = key_block if key_block.startswith("@") else "@" + key_block
        m = re.search(r"@\w+\{([^,]+),", block)
        if m and ("{" + m.group(1) + ",") not in old:
            additions.append(block)
    if additions:
        path.write_text(old.rstrip() + "\n\n" + "\n\n".join(additions) + "\n", encoding="utf-8")


def copy_reproducibility_inputs():
    # Actual retained RASPA force field definitions from the completed spatial/protocol freeze.
    base = PROJECT_ROOT / "path_b_analysis" / "paper_closure" / "02_spatial_maps" / "protocol_freeze"
    copy_if_exists(
        base / "M029" / "full_force_field.json",
        REPRO / "full_force_field_reference.json",
    )
    copy_if_exists(
        base / "M029" / "lj_only_force_field.json",
        REPRO / "lj_only_force_field_reference.json",
    )

    copy_if_exists(
        PROJECT_ROOT / "config" / "uff_framework_registry.csv",
        REPRO / "uff_framework_registry.csv",
    )


def manual_latex_build(texname: str) -> tuple[bool, str]:
    stem = Path(texname).stem
    cmds = [
        ["pdflatex", "-interaction=nonstopmode", "-halt-on-error", texname],
        ["bibtex", stem],
        ["pdflatex", "-interaction=nonstopmode", "-halt-on-error", texname],
        ["pdflatex", "-interaction=nonstopmode", "-halt-on-error", texname],
    ]

    logs = []
    ok = True
    for cmd in cmds:
        p = run(cmd, cwd=MANUSCRIPT)
        logs.append(
            f"$ {' '.join(cmd)}\nRETURN={p.returncode}\n"
            f"STDOUT:\n{p.stdout[-12000:]}\nSTDERR:\n{p.stderr[-4000:]}\n"
        )
        if p.returncode != 0:
            ok = False
            break

    pdf = MANUSCRIPT / f"{stem}.pdf"
    if not pdf.exists():
        ok = False

    return ok, "\n".join(logs)


def main():
    require(REPO / ".git")
    require(AUTH)

    status = git("status", "--porcelain").stdout.strip()
    if status:
        raise SystemExit(
            "STOP: publication repository is not clean. Commit/stash current changes first:\n\n"
            + status
        )

    current = git("branch", "--show-current").stdout.strip()
    branches = set(git("branch", "--format=%(refname:short)").stdout.split())

    if current != TARGET_BRANCH:
        if TARGET_BRANCH in branches:
            p = git("switch", TARGET_BRANCH)
        else:
            p = git("switch", "-c", TARGET_BRANCH)
        if p.returncode != 0:
            raise SystemExit(p.stderr)
        print(f"[git] switched to {TARGET_BRANCH}")

    d = pd.read_csv(AUTH)

    print("[1/9] Building objective matched-pair geometry context...")
    pair_context = build_pair_geometry_context(d)
    print(pair_context.to_string(index=False))

    print("\n[2/9] Building population-wide rank-transition source data...")
    rank_df, rank_stats = build_rank_source(d)
    print(json.dumps(rank_stats, indent=2))

    print("\n[3/9] Rebuilding Figure 2...")
    make_figure2(d)

    print("[4/9] Rebuilding Figure 4 around M029/M085...")
    make_figure4(d, pair_context)

    print("[5/9] Rebuilding Figure 5 from all 238 frameworks...")
    make_figure5(rank_df, rank_stats)

    print("[6/9] Rewriting main manuscript and SI...")
    (MANUSCRIPT / "main.tex").write_text(MAIN_TEX, encoding="utf-8")
    (MANUSCRIPT / "si.tex").write_text(SI_TEX, encoding="utf-8")
    append_bib_entries()
    (REPO / "docs" / "NOVELTY_AUDIT.md").write_text(NOVELTY_MD, encoding="utf-8")

    print("[7/9] Adding missing molecular-model provenance...")
    copy_reproducibility_inputs()

    print("[8/9] Running explicit LaTeX/BibTeX builds...")
    all_ok = True
    build_report = []
    for tex in ["main.tex", "si.tex"]:
        ok, log = manual_latex_build(tex)
        all_ok &= ok
        build_report.append(f"\n{'='*30} {tex} {'='*30}\n{log}")
        print(f"  {tex}: {'PASS' if ok else 'FAIL'}")

    report_path = REPRO / "pre_pi_manual_build_report.txt"
    report_path.write_text("\n".join(build_report), encoding="utf-8")

    # Build products are ignored by .gitignore. Keep PDFs locally for inspection.
    if not all_ok:
        print("\nBUILD DID NOT FULLY PASS.")
        print(f"Inspect: {report_path}")
        print("No commit will be created. Fix the reported TeX/BibTeX error and rerun.")
        sys.exit(2)

    print("[9/9] Staging and committing pre-PI publication state...")
    git("add", "-A", check=True)
    commit = git("commit", "-m", COMMIT_MESSAGE)
    if commit.returncode != 0:
        print(commit.stdout)
        print(commit.stderr)
        raise SystemExit("Changes are staged, but Git commit failed.")

    # Tag this reviewable pre-PI state.
    tags = set(git("tag").stdout.split())
    if "pre-pi-publication-v1" not in tags:
        tag = git("tag", "pre-pi-publication-v1")
        if tag.returncode != 0:
            print("Warning: could not create pre-pi-publication-v1 tag:")
            print(tag.stderr)

    print("\n" + "=" * 88)
    print("PRE-PI PUBLICATION REFINEMENT COMPLETE")
    print("=" * 88)
    print(f"Branch: {git('branch','--show-current').stdout.strip()}")
    print("Build:  main.tex PASS; si.tex PASS")
    print("Git status:")
    print(git("status", "--short").stdout or "  clean")
    print("\nRecent commits:")
    print(git("log", "--oneline", "-4").stdout)
    print("\nKey new source data:")
    print("  data/source_data/matched_pair_geometry_context.csv")
    print("  data/source_data/fig5_rank_transition_source.csv")
    print("  data/source_data/fig5_rank_transition_stats.json")
    print("\nKey revised figures:")
    print("  figures/main/fig2_mechanism_map.pdf")
    print("  figures/main/fig4_matched_electrostatic.pdf")
    print("  figures/main/fig5_rank_trajectories.pdf")
    print("\nNext step: package this commit and perform the final pre-PI line-by-line paper audit.")


if __name__ == "__main__":
    main()
