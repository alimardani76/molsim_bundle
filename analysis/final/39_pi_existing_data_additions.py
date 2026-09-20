#!/usr/bin/env python3
"""
PI roadmap, existing-data additions only.

Reads the frozen 238-framework mechanism table and existing geometry OOF
predictions. Generates:
  - continuous mechanism fingerprint
  - endpoint-matched / mechanism-distant pair table
  - upgraded rank-transition statistics
  - fixed-OOF bootstrap uncertainty summary
  - top-5% / leverage sensitivity summary

No molecular simulation is run here.

Selection rule for endpoint twins is frozen before inspection:
take the closest 5% of all pairwise |Delta ln S_full| values, then rank
those pairs by Euclidean distance in (Delta ln S_attr, Delta ln S_elec).
"""

from __future__ import annotations

import json
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import kendalltau, spearmanr

ROOT = Path(__file__).resolve().parents[2]
AUTH = ROOT / "data" / "authoritative" / "path_b_final_238_authoritative.csv"
PRED = ROOT / "data" / "source_data" / "descriptor_hierarchy_predictions.csv"
OUT = ROOT / "data" / "source_data"

SEED = 20260920
N_BOOT = 5000
TWIN_FRACTION = 0.05


def median(x):
    return float(np.median(np.asarray(x, dtype=float)))


def r2(y, p):
    y = np.asarray(y, dtype=float)
    p = np.asarray(p, dtype=float)
    return float(1.0 - np.sum((y - p) ** 2) / np.sum((y - y.mean()) ** 2))


def ci95(x):
    a = np.asarray(x, dtype=float)
    return [float(np.quantile(a, 0.025)), float(np.quantile(a, 0.975))]


def descending_rank(x):
    return pd.Series(x).rank(method="average", ascending=False).to_numpy()


def main():
    OUT.mkdir(parents=True, exist_ok=True)

    df = pd.read_csv(AUTH)
    df = df.loc[df["three_axis_common"].astype(bool)].copy()
    assert len(df) == 238

    elec_col = "delta_lnS_elec_final" if "delta_lnS_elec_final" in df.columns else "delta_lnS_elec"
    df["delta_attr"] = df["delta_lnS_attr"].astype(float)
    df["delta_elec"] = df[elec_col].astype(float)
    df["mechanism_magnitude"] = np.hypot(df["delta_attr"], df["delta_elec"])
    df["mechanism_angle_rad"] = np.arctan2(df["delta_elec"], df["delta_attr"])
    df["electrostatic_share"] = (
        df["delta_elec"].abs()
        / (df["delta_attr"].abs() + df["delta_elec"].abs())
    )

    fp_cols = [
        "mof_id", "lnS_wca", "lnS_lj_only", "lnS_full",
        "delta_attr", "delta_elec",
        "mechanism_magnitude", "mechanism_angle_rad", "electrostatic_share",
    ]
    df[fp_cols].to_csv(OUT / "pi_mechanism_fingerprint.csv", index=False)

    # Endpoint twins: freeze closest 5% of all endpoint differences, then rank by mechanism distance.
    vals = df[["mof_id", "lnS_full", "delta_attr", "delta_elec"]].reset_index(drop=True)
    recs = []
    for i in range(len(vals)):
        for j in range(i + 1, len(vals)):
            endpoint_diff = abs(vals.loc[i, "lnS_full"] - vals.loc[j, "lnS_full"])
            da = vals.loc[i, "delta_attr"] - vals.loc[j, "delta_attr"]
            de = vals.loc[i, "delta_elec"] - vals.loc[j, "delta_elec"]
            recs.append({
                "mof_i": vals.loc[i, "mof_id"],
                "mof_j": vals.loc[j, "mof_id"],
                "endpoint_diff_abs_lnS_full": endpoint_diff,
                "mechanism_distance": float(np.hypot(da, de)),
                "abs_delta_attr_difference": abs(da),
                "abs_delta_elec_difference": abs(de),
                "lnS_full_i": vals.loc[i, "lnS_full"],
                "lnS_full_j": vals.loc[j, "lnS_full"],
                "delta_attr_i": vals.loc[i, "delta_attr"],
                "delta_attr_j": vals.loc[j, "delta_attr"],
                "delta_elec_i": vals.loc[i, "delta_elec"],
                "delta_elec_j": vals.loc[j, "delta_elec"],
            })
    pairs = pd.DataFrame(recs).sort_values("endpoint_diff_abs_lnS_full", kind="stable")
    n_keep = int(np.ceil(TWIN_FRACTION * len(pairs)))
    twins = pairs.iloc[:n_keep].copy()
    endpoint_cutoff = float(twins["endpoint_diff_abs_lnS_full"].max())
    twins = twins.sort_values("mechanism_distance", ascending=False, kind="stable")
    twins.to_csv(OUT / "pi_performance_twin_pairs.csv", index=False)

    # Ranking extension.
    rw = descending_rank(df["lnS_wca"].to_numpy())
    rlj = descending_rank(df["lnS_lj_only"].to_numpy())
    rf = descending_rank(df["lnS_full"].to_numpy())
    sh_attr = rw - rlj
    sh_elec = rlj - rf
    abs_attr = np.abs(sh_attr)
    abs_elec = np.abs(sh_elec)
    n = len(df)

    def top_overlap(a, b, k):
        ia = set(np.argsort(-np.asarray(a))[:k])
        ib = set(np.argsort(-np.asarray(b))[:k])
        return len(ia & ib)

    thresholds = {
        "10pct": int(np.ceil(0.10 * n)),
        "25pct": int(np.ceil(0.25 * n)),
        "50pct": int(np.ceil(0.50 * n)),
    }

    rank_summary = {
        "n": n,
        "wca_to_lj": {
            "spearman": float(spearmanr(df["lnS_wca"], df["lnS_lj_only"]).statistic),
            "kendall_tau": float(kendalltau(df["lnS_wca"], df["lnS_lj_only"]).statistic),
            "median_abs_rank_shift": median(abs_attr),
            "max_promotion": float(np.max(sh_attr)),
            "max_demotion": float(-np.min(sh_attr)),
            "movement_counts": {
                name: {
                    "threshold_positions": t,
                    "count": int(np.sum(abs_attr >= t)),
                    "fraction": float(np.mean(abs_attr >= t)),
                }
                for name, t in thresholds.items()
            },
            "top_overlap": {
                "top10": top_overlap(df["lnS_wca"], df["lnS_lj_only"], 10),
                "top25": top_overlap(df["lnS_wca"], df["lnS_lj_only"], 25),
                "top10pct": top_overlap(df["lnS_wca"], df["lnS_lj_only"], int(np.ceil(0.10*n))),
            },
        },
        "lj_to_full": {
            "spearman": float(spearmanr(df["lnS_lj_only"], df["lnS_full"]).statistic),
            "kendall_tau": float(kendalltau(df["lnS_lj_only"], df["lnS_full"]).statistic),
            "median_abs_rank_shift": median(abs_elec),
            "max_promotion": float(np.max(sh_elec)),
            "max_demotion": float(-np.min(sh_elec)),
            "movement_counts": {
                name: {
                    "threshold_positions": t,
                    "count": int(np.sum(abs_elec >= t)),
                    "fraction": float(np.mean(abs_elec >= t)),
                }
                for name, t in thresholds.items()
            },
            "top_overlap": {
                "top10": top_overlap(df["lnS_lj_only"], df["lnS_full"], 10),
                "top25": top_overlap(df["lnS_lj_only"], df["lnS_full"], 25),
                "top10pct": top_overlap(df["lnS_lj_only"], df["lnS_full"], int(np.ceil(0.10*n))),
            },
        },
    }
    (OUT / "pi_rank_extension_summary.json").write_text(
        json.dumps(rank_summary, indent=2) + "\n", encoding="utf-8"
    )

    # Existing geometry OOF predictions. Bootstrap the realized OOF metrics without refitting.
    pred = pd.read_csv(PRED)
    pred = pred.loc[pred["feature_set"].eq("geometry")].copy()
    wide_obs = pred.pivot(index="mof_id", columns="target", values="observed")
    wide_pre = pred.pivot(index="mof_id", columns="target", values="predicted")
    ids = df["mof_id"].tolist()
    ao = wide_obs.loc[ids, "delta_attr"].to_numpy()
    ap = wide_pre.loc[ids, "delta_attr"].to_numpy()
    eo = wide_obs.loc[ids, "delta_elec"].to_numpy()
    ep = wide_pre.loc[ids, "delta_elec"].to_numpy()

    rng = np.random.default_rng(SEED)
    boot = {
        "median_attr": [], "median_elec": [], "rho_attr_elec": [],
        "fraction_attr_positive": [], "fraction_elec_positive": [],
        "r2_attr": [], "r2_elec": [], "predictability_gap": [],
    }
    a = df["delta_attr"].to_numpy()
    e = df["delta_elec"].to_numpy()
    for _ in range(N_BOOT):
        ix = rng.integers(0, n, n)
        ba = a[ix]
        be = e[ix]
        ra = r2(ao[ix], ap[ix])
        re = r2(eo[ix], ep[ix])
        boot["median_attr"].append(median(ba))
        boot["median_elec"].append(median(be))
        boot["rho_attr_elec"].append(float(spearmanr(ba, be).statistic))
        boot["fraction_attr_positive"].append(float(np.mean(ba > 0)))
        boot["fraction_elec_positive"].append(float(np.mean(be > 0)))
        boot["r2_attr"].append(ra)
        boot["r2_elec"].append(re)
        boot["predictability_gap"].append(ra - re)

    point_r2_attr = r2(ao, ap)
    point_r2_elec = r2(eo, ep)
    bootstrap_summary = {
        "method": "5000 framework-level bootstrap resamples of fixed existing OOF predictions; models are not refit inside bootstrap",
        "seed": SEED,
        "n": n,
        "point": {
            "median_attr": median(a),
            "median_elec": median(e),
            "rho_attr_elec": float(spearmanr(a, e).statistic),
            "fraction_attr_positive": float(np.mean(a > 0)),
            "fraction_elec_positive": float(np.mean(e > 0)),
            "r2_attr": point_r2_attr,
            "r2_elec": point_r2_elec,
            "predictability_gap": point_r2_attr - point_r2_elec,
        },
        "ci95": {k: ci95(v) for k, v in boot.items()},
    }
    (OUT / "pi_bootstrap_summary.json").write_text(
        json.dumps(bootstrap_summary, indent=2) + "\n", encoding="utf-8"
    )

    # Minimal extra outlier/leverage sensitivity.
    def trim_top(column, frac):
        k = int(np.ceil(n * frac))
        keep = df.drop(df.nlargest(k, column).index)
        return {
            "population": f"exclude_top_{int(frac*100)}pct_{column}",
            "n": len(keep),
            "median_delta_attr": median(keep["delta_attr"]),
            "median_delta_elec": median(keep["delta_elec"]),
            "attr_elec_rho": float(spearmanr(keep["delta_attr"], keep["delta_elec"]).statistic),
        }

    sens = [
        trim_top("delta_attr", 0.05),
        trim_top("delta_elec", 0.05),
        trim_top("lnS_full", 0.05),
    ]
    for name, shifts in [
        ("exclude_10_largest_attr_rank_movers", abs_attr),
        ("exclude_10_largest_elec_rank_movers", abs_elec),
    ]:
        drop = np.argsort(-shifts)[:10]
        keep = df.drop(df.index[drop])
        sens.append({
            "population": name,
            "n": len(keep),
            "median_delta_attr": median(keep["delta_attr"]),
            "median_delta_elec": median(keep["delta_elec"]),
            "attr_elec_rho": float(spearmanr(keep["delta_attr"], keep["delta_elec"]).statistic),
        })
    pd.DataFrame(sens).to_csv(OUT / "pi_outlier_sensitivity.csv", index=False)

    report = {
        "endpoint_twin_selection": {
            "fraction_of_all_pairs": TWIN_FRACTION,
            "number_of_all_pairs": int(len(pairs)),
            "number_selected": int(n_keep),
            "endpoint_difference_cutoff": endpoint_cutoff,
            "mechanism_distance_median": float(twins["mechanism_distance"].median()),
            "mechanism_distance_q90": float(twins["mechanism_distance"].quantile(0.90)),
            "mechanism_distance_q95": float(twins["mechanism_distance"].quantile(0.95)),
            "mechanism_distance_q99": float(twins["mechanism_distance"].quantile(0.99)),
        },
        "mechanism_fingerprint": {
            "median_magnitude": median(df["mechanism_magnitude"]),
            "median_angle_rad": median(df["mechanism_angle_rad"]),
            "median_electrostatic_share": median(df["electrostatic_share"]),
        },
    }
    (OUT / "pi_existing_data_report.json").write_text(
        json.dumps(report, indent=2) + "\n", encoding="utf-8"
    )

    print("PASS")
    print(json.dumps(report, indent=2))
    print(json.dumps(rank_summary, indent=2))
    print(json.dumps(bootstrap_summary, indent=2))


if __name__ == "__main__":
    main()
