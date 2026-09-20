from __future__ import annotations

import csv
import json
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve()
ROOT = HERE.parents[1]

SPATIAL = ROOT / "path_b_analysis" / "paper_closure" / "02_spatial_maps"
FIELDS = SPATIAL / "fields"
OUT = SPATIAL / "results" / "final_spatial_postprocess"
MAPS = OUT / "maps"

EXPECTED = [
    ("M089", "CO2", "lj_only", "attraction_pair"),
    ("M089", "N2",  "lj_only", "attraction_pair"),
    ("M192", "CO2", "lj_only", "attraction_pair"),
    ("M192", "N2",  "lj_only", "attraction_pair"),
    ("M029", "CO2", "lj_only", "electrostatic_pair"),
    ("M029", "N2",  "lj_only", "electrostatic_pair"),
    ("M029", "CO2", "full",    "electrostatic_pair"),
    ("M029", "N2",  "full",    "electrostatic_pair"),
    ("M085", "CO2", "lj_only", "electrostatic_pair"),
    ("M085", "N2",  "lj_only", "electrostatic_pair"),
    ("M085", "CO2", "full",    "electrostatic_pair"),
    ("M085", "N2",  "full",    "electrostatic_pair"),
]

def dump_json(path: Path, obj):
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(obj, indent=2, allow_nan=False) + "\n", encoding="utf-8")

def candidates(mof, gas, mode):
    if mode == "lj_only":
        return [FIELDS / f"{mof}_{gas}_lj_only_0.75A_2048orient"]
    return [
        FIELDS / f"{mof}_{gas}_full_0.75A_adaptive",
        FIELDS / f"{mof}_{gas}_full_0.75A_32768orient",
    ]

def load_field(mof, gas, mode, role):
    folder = next((p for p in candidates(mof, gas, mode)
                   if (p/"field.npz").exists() and (p/"metadata.json").exists()), None)
    if folder is None:
        raise FileNotFoundError(f"Missing completed field: {mof}/{gas}/{mode}")

    d = np.load(folder/"field.npz")
    meta = json.loads((folder/"metadata.json").read_text(encoding="utf-8"))

    shape = tuple(int(x) for x in d["shape"])
    n = int(np.prod(shape))
    frac = np.asarray(d["fractional_coordinates"], float)
    cell = np.asarray(d["unit_cell_vectors_A"], float)
    logw = np.asarray(d["log_weight"], float).reshape(shape)

    if not np.isfinite(logw).all():
        raise RuntimeError(f"{mof}/{gas}/{mode}: log_weight contains non-finite values")

    log_mean = meta.get("log_cell_mean_weight")
    if log_mean is None:
        # Stable log(mean(exp(logw))) without exponentiating huge/small values.
        flat = logw.ravel()
        m = float(np.max(flat))
        log_mean = float(m + np.log(np.mean(np.exp(flat - m))))
    else:
        log_mean = float(log_mean)

    # IMPORTANT: never use np.log(relative_density), because relative_density may
    # underflow to zero in strongly repulsive voxels.
    logrel = logw - log_mean

    if not np.isfinite(logrel).all():
        raise RuntimeError(f"{mof}/{gas}/{mode}: stable log-relative-density failed")

    return {
        "mof": mof, "gas": gas, "mode": mode, "role": role,
        "folder": str(folder.relative_to(SPATIAL)),
        "shape": shape, "frac": frac, "cell": cell,
        "logw": logw, "logrel": logrel, "meta": meta,
    }

def same_grid(a, b, label):
    if a["shape"] != b["shape"]:
        raise RuntimeError(f"{label}: shape mismatch")
    if not np.allclose(a["frac"], b["frac"], rtol=0, atol=1e-12):
        raise RuntimeError(f"{label}: fractional grid mismatch")
    if not np.allclose(a["cell"], b["cell"], rtol=0, atol=1e-12):
        raise RuntimeError(f"{label}: unit-cell mismatch")

def projections(a):
    return {
        "xy_mean": np.mean(a, axis=2),
        "xz_mean": np.mean(a, axis=1),
        "yz_mean": np.mean(a, axis=0),
        "xy_max": np.max(a, axis=2),
        "xz_max": np.max(a, axis=1),
        "yz_max": np.max(a, axis=0),
    }

def extents(cell):
    L = np.linalg.norm(cell, axis=1)
    return {
        "xy": np.array([0.0, L[0], 0.0, L[1]], float),
        "xz": np.array([0.0, L[0], 0.0, L[2]], float),
        "yz": np.array([0.0, L[1], 0.0, L[2]], float),
        "lengths": L,
    }

def summarize(a):
    x = np.asarray(a, float).ravel()
    if not np.isfinite(x).all():
        raise RuntimeError("Attempted to summarize non-finite array")
    q = np.quantile(x, [0.01, 0.05, 0.50, 0.95, 0.99])
    return {
        "min": float(x.min()), "q01": float(q[0]), "q05": float(q[1]),
        "median": float(q[2]), "q95": float(q[3]), "q99": float(q[4]),
        "max": float(x.max()), "mean": float(x.mean()), "std": float(x.std()),
    }

def save_map(name, array3d, cell, kind):
    if not np.isfinite(array3d).all():
        raise RuntimeError(f"{name}: non-finite derived map")
    e = extents(cell)
    p = projections(array3d)
    path = MAPS / f"{name}.npz"
    np.savez_compressed(
        path,
        array_3d=array3d,
        **p,
        extent_xy=e["xy"], extent_xz=e["xz"], extent_yz=e["yz"],
        lengths_A=e["lengths"],
    )
    return {
        "name": name,
        "kind": kind,
        "path": str(path.relative_to(OUT)),
        "summary": summarize(array3d),
    }

def main():
    OUT.mkdir(parents=True, exist_ok=True)
    MAPS.mkdir(parents=True, exist_ok=True)

    print("="*104)
    print("POSTPROCESS SPATIAL FIELDS V2 — STABLE LOG-SPACE")
    print("No simulation. No log(0).")
    print("="*104)

    fields = {}
    inventory = []
    derived = []

    for mof, gas, mode, role in EXPECTED:
        f = load_field(mof, gas, mode, role)
        fields[(mof, gas, mode)] = f
        inventory.append({
            "role": role, "mof": mof, "gas": gas, "mode": mode,
            "folder": f["folder"],
            "shape": "x".join(map(str, f["shape"])),
            "log_cell_mean_weight": f["meta"].get("log_cell_mean_weight"),
            "cell_mean_weight": f["meta"].get("cell_mean_weight"),
        })
        derived.append(save_map(
            f"{mof}_{gas}_{mode}_log_relative_density",
            f["logrel"], f["cell"], "log_relative_density"
        ))
        print(f"Loaded {mof}/{gas}/{mode}", flush=True)

    with (OUT/"field_inventory.csv").open("w", newline="", encoding="utf-8") as fh:
        w = csv.DictWriter(fh, fieldnames=list(inventory[0]))
        w.writeheader(); w.writerows(inventory)

    rows = []

    # Local CO2/N2 maps.
    # Two forms are saved:
    # 1) log_weight_ratio: absolute local Boltzmann-weight contrast.
    # 2) localization_contrast: normalized redistribution only.
    for mof in ["M089", "M192", "M029", "M085"]:
        modes = ["lj_only"] if mof in ("M089","M192") else ["lj_only","full"]
        for mode in modes:
            c = fields[(mof,"CO2",mode)]
            n = fields[(mof,"N2",mode)]
            same_grid(c, n, f"{mof}/{mode}/CO2-N2")

            log_ratio = c["logw"] - n["logw"]
            loc_contrast = c["logrel"] - n["logrel"]

            for suffix, arr, kind in [
                ("ln_local_weight_ratio_CO2_over_N2", log_ratio,
                 "ln_local_weight_ratio_CO2_over_N2"),
                ("ln_localization_contrast_CO2_over_N2", loc_contrast,
                 "ln_localization_contrast_CO2_over_N2"),
            ]:
                rec = save_map(f"{mof}_{mode}_{suffix}", arr, c["cell"], kind)
                derived.append(rec)
                rows.append({
                    "mof": mof, "gas": "CO2_over_N2", "mode": mode,
                    "kind": kind, **rec["summary"]
                })

    # Electrostatic increment = FULL - LJ on the exact same spatial grid.
    for mof in ["M029","M085"]:
        for gas in ["CO2","N2"]:
            lj = fields[(mof,gas,"lj_only")]
            full = fields[(mof,gas,"full")]
            same_grid(lj, full, f"{mof}/{gas}/full-LJ")

            delta_logw = full["logw"] - lj["logw"]
            delta_logrel = full["logrel"] - lj["logrel"]

            for suffix, arr, kind in [
                ("delta_log_weight_full_minus_lj", delta_logw,
                 "electrostatic_increment_log_weight"),
                ("delta_log_relative_density_full_minus_lj", delta_logrel,
                 "electrostatic_redistribution_log_relative_density"),
            ]:
                rec = save_map(f"{mof}_{gas}_{suffix}", arr, full["cell"], kind)
                derived.append(rec)
                rows.append({
                    "mof": mof, "gas": gas, "mode": "full_minus_lj",
                    "kind": kind, **rec["summary"]
                })

        # Most direct spatial analogue of electrostatic selectivity increment:
        # [(CO2-N2)_full - (CO2-N2)_LJ]
        c_lj = fields[(mof,"CO2","lj_only")]
        n_lj = fields[(mof,"N2","lj_only")]
        c_fu = fields[(mof,"CO2","full")]
        n_fu = fields[(mof,"N2","full")]

        delta_sel = (c_fu["logw"] - n_fu["logw"]) - (c_lj["logw"] - n_lj["logw"])
        rec = save_map(
            f"{mof}_delta_ln_local_weight_ratio_full_minus_lj",
            delta_sel, c_fu["cell"],
            "electrostatic_increment_ln_local_CO2_N2_weight_ratio"
        )
        derived.append(rec)
        rows.append({
            "mof": mof, "gas": "CO2_over_N2", "mode": "full_minus_lj",
            "kind": rec["kind"], **rec["summary"]
        })

    with (OUT/"derived_map_summary.csv").open("w", newline="", encoding="utf-8") as fh:
        w = csv.DictWriter(fh, fieldnames=list(rows[0]))
        w.writeheader(); w.writerows(rows)

    report = {
        "postprocess_passed": True,
        "stable_log_space": True,
        "fields_loaded": len(fields),
        "expected_fields": len(EXPECTED),
        "derived_maps": derived,
        "definitions": {
            "log_relative_density": "log_weight - log(cell mean weight); computed directly in log-space.",
            "ln_local_weight_ratio_CO2_over_N2": "log(w_CO2) - log(w_N2) at the same voxel.",
            "ln_localization_contrast_CO2_over_N2": "log(rho_CO2/<rho_CO2>) - log(rho_N2/<rho_N2>); localization only.",
            "electrostatic_increment_log_weight": "log(w_full) - log(w_LJ) at the same voxel.",
            "electrostatic_increment_ln_local_CO2_N2_weight_ratio": "[(CO2-N2)_full - (CO2-N2)_LJ] in log weight space.",
        },
        "important_scope": [
            "M089/M192 LJ maps illustrate spatial organization of the attraction-selected matched pair; without WCA spatial maps they do not isolate attraction itself.",
            "M029/M085 full-minus-LJ maps directly isolate the electrostatic spatial increment at fixed LJ definitions.",
        ],
        "next_step": "Run 37_make_matched_pair_spatial_figure_v2.py",
    }
    dump_json(OUT/"POSTPROCESS_REPORT.json", report)

    print("="*104)
    print("PASS — stable spatial post-processing complete.")
    print("Next: run script 37 v2.")
    print("="*104)

if __name__ == "__main__":
    main()
