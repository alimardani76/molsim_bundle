from __future__ import annotations

import csv
import json
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np

HERE = Path(__file__).resolve()
ROOT = HERE.parents[1]

SPATIAL = ROOT / "path_b_analysis" / "paper_closure" / "02_spatial_maps"
FIELDS = SPATIAL / "fields"
OUT = SPATIAL / "results" / "final_spatial_figures_supported"
OUT.mkdir(parents=True, exist_ok=True)

MASS_COVERAGE = 0.995  # retain voxels carrying 99.5% of insertion probability mass

def candidates(mof, gas, mode):
    if mode == "lj_only":
        return [FIELDS / f"{mof}_{gas}_lj_only_0.75A_2048orient"]
    return [
        FIELDS / f"{mof}_{gas}_full_0.75A_adaptive",
        FIELDS / f"{mof}_{gas}_full_0.75A_32768orient",
    ]

def load_field(mof, gas, mode):
    folder = next((p for p in candidates(mof, gas, mode)
                   if (p/"field.npz").exists()), None)
    if folder is None:
        raise FileNotFoundError(f"Missing field {mof}/{gas}/{mode}")

    d = np.load(folder/"field.npz")
    shape = tuple(int(x) for x in d["shape"])
    logw = np.asarray(d["log_weight"], float).reshape(shape)
    rel = np.asarray(d["relative_density"], float).reshape(shape)
    cell = np.asarray(d["unit_cell_vectors_A"], float)

    if not np.isfinite(logw).all():
        raise RuntimeError(f"Non-finite log_weight in {folder}")
    if not np.isfinite(rel).all():
        raise RuntimeError(f"Non-finite relative_density in {folder}")

    return {
        "mof": mof, "gas": gas, "mode": mode,
        "folder": folder, "shape": shape,
        "logw": logw, "rel": rel, "cell": cell
    }

def assert_same_grid(*fs):
    s = fs[0]["shape"]
    cell = fs[0]["cell"]
    for f in fs[1:]:
        if f["shape"] != s:
            raise RuntimeError("Grid-shape mismatch")
        if not np.allclose(f["cell"], cell, rtol=0, atol=1e-12):
            raise RuntimeError("Unit-cell mismatch")

def support_mask(rel, coverage=MASS_COVERAGE):
    """
    Minimal voxel set carrying `coverage` of total insertion probability mass.
    Zero/underflow voxels are never selected unless all mass is zero.
    """
    x = np.asarray(rel, float).ravel()
    total = float(x.sum())
    if not np.isfinite(total) or total <= 0:
        raise RuntimeError("Relative-density field has no positive finite mass")

    order = np.argsort(x)[::-1]
    cum = np.cumsum(x[order]) / total
    k = int(np.searchsorted(cum, coverage, side="left"))
    keep = order[:k+1]

    mask = np.zeros_like(x, dtype=bool)
    mask[keep] = True
    return mask.reshape(rel.shape)

def support_union(*rels):
    m = np.zeros_like(rels[0], dtype=bool)
    for r in rels:
        m |= support_mask(r)
    return m

def weighted_z_projection(values, weights, support):
    """
    Probability-weighted mean along z.
    Regions with no supported probability become NaN rather than huge artifacts.
    """
    v = np.asarray(values, float)
    w = np.where(support, np.asarray(weights, float), 0.0)

    num = np.sum(np.where(support, v, 0.0) * w, axis=2)
    den = np.sum(w, axis=2)

    out = np.full_like(num, np.nan, dtype=float)
    ok = den > 0
    out[ok] = num[ok] / den[ok]
    return out

def projected_log_density(rel):
    """
    Mean relative density along z, then log10 for dynamic-range compression.
    This is a probability projection, not an average of log weights.
    """
    p = np.mean(rel, axis=2)
    out = np.full_like(p, np.nan, dtype=float)
    ok = p > 0
    out[ok] = np.log10(p[ok])
    return out

def extent_xy(cell):
    lengths = np.linalg.norm(cell, axis=1)
    return [0.0, float(lengths[0]), 0.0, float(lengths[1])]

def finite_concat(arrays):
    xs = []
    for a in arrays:
        x = np.asarray(a, float)
        x = x[np.isfinite(x)]
        if len(x):
            xs.append(x)
    if not xs:
        raise RuntimeError("No finite plotted values")
    return np.concatenate(xs)

def robust_limits(arrays, symmetric=False):
    x = finite_concat(arrays)
    if symmetric:
        q = float(np.quantile(np.abs(x), 0.99))
        if q <= 0:
            q = 1.0
        return -q, q
    lo, hi = np.quantile(x, [0.01, 0.99])
    if lo == hi:
        hi = lo + 1.0
    return float(lo), float(hi)

def draw(ax, arr, extent, title, limits):
    m = np.ma.masked_invalid(arr)
    im = ax.imshow(
        m.T, origin="lower", aspect="auto", extent=extent,
        vmin=limits[0], vmax=limits[1]
    )
    ax.set_title(title, fontsize=9)
    ax.set_xlabel("x (Å)")
    ax.set_ylabel("y (Å)")
    return im

def save(fig, stem):
    fig.savefig(OUT/f"{stem}.png", dpi=300, bbox_inches="tight")
    fig.savefig(OUT/f"{stem}.pdf", bbox_inches="tight")
    plt.close(fig)

def support_stats(mask):
    return {
        "voxels": int(mask.sum()),
        "fraction_of_grid": float(mask.mean()),
    }

def main():
    print("="*105)
    print("PROBABILITY-SUPPORTED SPATIAL FIGURES")
    print("No simulation. Inaccessible overlap voxels excluded by 99.5% insertion-mass support.")
    print("="*105)

    stats = []

    # ------------------------------------------------------------------
    # FIGURE A: attraction-selected matched pair
    # Probability projections + supported local LJ-stage CO2/N2 contrast.
    # This does NOT isolate attraction spatially (no WCA spatial maps).
    # ------------------------------------------------------------------
    attraction = {}
    for mof in ["M089", "M192"]:
        c = load_field(mof, "CO2", "lj_only")
        n = load_field(mof, "N2", "lj_only")
        assert_same_grid(c, n)

        support = support_union(c["rel"], n["rel"])
        weights = c["rel"] + n["rel"]

        co2_proj = projected_log_density(c["rel"])
        n2_proj = projected_log_density(n["rel"])

        local_ratio = c["logw"] - n["logw"]
        ratio_proj = weighted_z_projection(local_ratio, weights, support)

        attraction[(mof, "co2")] = (co2_proj, extent_xy(c["cell"]))
        attraction[(mof, "n2")] = (n2_proj, extent_xy(c["cell"]))
        attraction[(mof, "ratio")] = (ratio_proj, extent_xy(c["cell"]))

        s = support_stats(support)
        stats.append({
            "figure_group": "attraction_selected",
            "mof": mof,
            "support_voxels": s["voxels"],
            "support_fraction_of_grid": s["fraction_of_grid"],
        })

    lim_c = robust_limits([attraction[(m,"co2")][0] for m in ["M089","M192"]])
    lim_n = robust_limits([attraction[(m,"n2")][0] for m in ["M089","M192"]])
    lim_r = robust_limits([attraction[(m,"ratio")][0] for m in ["M089","M192"]], symmetric=True)

    fig, axs = plt.subplots(2, 3, figsize=(12, 7), constrained_layout=True)
    for r, mof in enumerate(["M089","M192"]):
        ims = [
            draw(axs[r,0], *attraction[(mof,"co2")],
                 f"{mof}: CO$_2$ LJ projected log$_{{10}}$ density", lim_c),
            draw(axs[r,1], *attraction[(mof,"n2")],
                 f"{mof}: N$_2$ LJ projected log$_{{10}}$ density", lim_n),
            draw(axs[r,2], *attraction[(mof,"ratio")],
                 f"{mof}: supported LJ local CO$_2$/N$_2$ contrast", lim_r),
        ]
        for cidx, im in enumerate(ims):
            fig.colorbar(im, ax=axs[r,cidx], shrink=0.82)

    fig.suptitle(
        "Attraction-selected matched pair: probability-supported LJ spatial organization",
        fontsize=12
    )
    save(fig, "figure_attraction_selected_pair_supported_xy")

    # ------------------------------------------------------------------
    # FIGURE B: electrostatic pair
    # LJ ratio, full ratio, and direct electrostatic selectivity increment.
    # All projected only where any of the four fields carries probability.
    # ------------------------------------------------------------------
    electro = {}
    for mof in ["M029","M085"]:
        clj = load_field(mof, "CO2", "lj_only")
        nlj = load_field(mof, "N2",  "lj_only")
        cfu = load_field(mof, "CO2", "full")
        nfu = load_field(mof, "N2",  "full")
        assert_same_grid(clj, nlj, cfu, nfu)

        support = support_union(clj["rel"], nlj["rel"], cfu["rel"], nfu["rel"])
        weights = clj["rel"] + nlj["rel"] + cfu["rel"] + nfu["rel"]

        lj_ratio = clj["logw"] - nlj["logw"]
        fu_ratio = cfu["logw"] - nfu["logw"]
        delta = fu_ratio - lj_ratio

        electro[(mof,"lj")] = (
            weighted_z_projection(lj_ratio, weights, support),
            extent_xy(clj["cell"])
        )
        electro[(mof,"full")] = (
            weighted_z_projection(fu_ratio, weights, support),
            extent_xy(clj["cell"])
        )
        electro[(mof,"delta")] = (
            weighted_z_projection(delta, weights, support),
            extent_xy(clj["cell"])
        )

        s = support_stats(support)
        stats.append({
            "figure_group": "electrostatic_pair",
            "mof": mof,
            "support_voxels": s["voxels"],
            "support_fraction_of_grid": s["fraction_of_grid"],
        })

    lim_abs = robust_limits(
        [electro[(m,k)][0] for m in ["M029","M085"] for k in ["lj","full"]],
        symmetric=True
    )
    lim_delta = robust_limits(
        [electro[(m,"delta")][0] for m in ["M029","M085"]],
        symmetric=True
    )

    fig, axs = plt.subplots(2, 3, figsize=(12, 7), constrained_layout=True)
    for r, mof in enumerate(["M029","M085"]):
        ims = [
            draw(axs[r,0], *electro[(mof,"lj")],
                 f"{mof}: supported LJ local CO$_2$/N$_2$ contrast", lim_abs),
            draw(axs[r,1], *electro[(mof,"full")],
                 f"{mof}: supported full local CO$_2$/N$_2$ contrast", lim_abs),
            draw(axs[r,2], *electro[(mof,"delta")],
                 f"{mof}: supported electrostatic selectivity increment", lim_delta),
        ]
        for cidx, im in enumerate(ims):
            fig.colorbar(im, ax=axs[r,cidx], shrink=0.82)

    fig.suptitle(
        "Electrostatic matched pair: probability-supported full − LJ selectivity increment",
        fontsize=12
    )
    save(fig, "figure_electrostatic_pair_supported_xy")

    # ------------------------------------------------------------------
    # FIGURE C: per-gas electrostatic increment
    # ------------------------------------------------------------------
    gas_inc = {}
    for mof in ["M029","M085"]:
        for gas in ["CO2","N2"]:
            lj = load_field(mof, gas, "lj_only")
            fu = load_field(mof, gas, "full")
            assert_same_grid(lj, fu)

            support = support_union(lj["rel"], fu["rel"])
            weights = lj["rel"] + fu["rel"]
            delta = fu["logw"] - lj["logw"]

            gas_inc[(mof,gas)] = (
                weighted_z_projection(delta, weights, support),
                extent_xy(lj["cell"])
            )

    lim_g = robust_limits(
        [gas_inc[(m,g)][0] for m in ["M029","M085"] for g in ["CO2","N2"]],
        symmetric=True
    )

    fig, axs = plt.subplots(2, 2, figsize=(9, 7), constrained_layout=True)
    for r, mof in enumerate(["M029","M085"]):
        im0 = draw(axs[r,0], *gas_inc[(mof,"CO2")],
                   f"{mof}: supported CO$_2$ electrostatic increment", lim_g)
        im1 = draw(axs[r,1], *gas_inc[(mof,"N2")],
                   f"{mof}: supported N$_2$ electrostatic increment", lim_g)
        fig.colorbar(im0, ax=axs[r,0], shrink=0.82)
        fig.colorbar(im1, ax=axs[r,1], shrink=0.82)

    fig.suptitle(
        "Per-gas electrostatic spatial increment on probability-supported regions",
        fontsize=12
    )
    save(fig, "figure_electrostatic_gas_increment_supported_xy")

    with (OUT/"support_statistics.csv").open("w", newline="", encoding="utf-8") as fh:
        w = csv.DictWriter(fh, fieldnames=list(stats[0]))
        w.writeheader()
        w.writerows(stats)

    report = {
        "passed": True,
        "mass_coverage_per_field": MASS_COVERAGE,
        "figures": [
            "figure_attraction_selected_pair_supported_xy.png",
            "figure_attraction_selected_pair_supported_xy.pdf",
            "figure_electrostatic_pair_supported_xy.png",
            "figure_electrostatic_pair_supported_xy.pdf",
            "figure_electrostatic_gas_increment_supported_xy.png",
            "figure_electrostatic_gas_increment_supported_xy.pdf",
        ],
        "method": [
            "Probability support is defined separately for each field as the minimal voxel set carrying 99.5% of insertion probability mass.",
            "Pairwise support is the union of the relevant field supports.",
            "Density panels project linear relative density along z and display log10 of that projected density.",
            "Local selectivity and full-minus-LJ increments are probability-weighted along z rather than averaged uniformly over inaccessible space.",
            "Display limits use robust 1st-99th percentile limits; underlying data are not clipped.",
        ],
        "scope": [
            "M089/M192 illustrate LJ-stage spatial organization for the attraction-selected matched pair; they do not spatially isolate attraction because no WCA field was computed.",
            "M029/M085 full-minus-LJ directly isolates the electrostatic spatial increment.",
        ],
    }
    (OUT/"SUPPORTED_FIGURE_REPORT.json").write_text(
        json.dumps(report, indent=2) + "\n", encoding="utf-8"
    )

    print("="*105)
    print("PASS — probability-supported figures created.")
    for name in report["figures"]:
        print(" ", OUT/name)
    print(" ", OUT/"support_statistics.csv")
    print(" ", OUT/"SUPPORTED_FIGURE_REPORT.json")
    print("="*105)

if __name__ == "__main__":
    main()
