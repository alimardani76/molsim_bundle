#!/usr/bin/env python3
"""
Evaluate the frozen PI Gate C electrostatic-probe descriptor table.

Input expected:
data/source_data/pi_electrostatic_probe_descriptors.csv

The generator must follow docs/PI_GATE_C_FREEZE.md and
reproducibility/pi_gate_c_frozen_protocol.json.

This script performs the complete pre-specified decision test and closes Gate C.
"""

from __future__ import annotations

import json
from pathlib import Path

import numpy as np
import pandas as pd
from sklearn.compose import ColumnTransformer
from sklearn.ensemble import ExtraTreesRegressor
from sklearn.impute import SimpleImputer
from sklearn.metrics import mean_absolute_error, r2_score
from sklearn.model_selection import KFold, cross_val_predict
from sklearn.pipeline import Pipeline
from sklearn.preprocessing import StandardScaler

ROOT = Path(__file__).resolve().parents[2]
AUTH = ROOT / "data" / "authoritative" / "path_b_final_238_authoritative.csv"
PROBE = ROOT / "data" / "source_data" / "pi_electrostatic_probe_descriptors.csv"
MANIFEST = ROOT / "data" / "source_data" / "descriptor_feature_manifest.csv"
OUT = ROOT / "data" / "source_data"
DOCS = ROOT / "docs"

SEED = 20260918
BOOT_SEED = 20260920
N_BOOT = 5000

def make_model():
    return ExtraTreesRegressor(
        n_estimators=500,
        min_samples_leaf=3,
        max_features=0.7,
        random_state=SEED,
        n_jobs=-1,
    )

def cv_predict(X: pd.DataFrame, y: np.ndarray):
    nums = list(X.columns)
    prep = ColumnTransformer([
        ("num", Pipeline([
            ("imp", SimpleImputer(strategy="median")),
            ("sc", StandardScaler()),
        ]), nums)
    ])
    pipe = Pipeline([("prep", prep), ("model", make_model())])
    cv = KFold(n_splits=5, shuffle=True, random_state=SEED)
    return cross_val_predict(pipe, X, y, cv=cv, n_jobs=1)

def metrics(y, yp):
    baseline = np.repeat(np.median(y), len(y))
    return {
        "r2": float(r2_score(y, yp)),
        "mae": float(mean_absolute_error(y, yp)),
        "baseline_mae": float(mean_absolute_error(y, baseline)),
    }

def ci95(x):
    return [float(np.quantile(x, 0.025)), float(np.quantile(x, 0.975))]

def main():
    for p in [AUTH, PROBE, MANIFEST]:
        if not p.exists():
            raise FileNotFoundError(p)

    df = pd.read_csv(AUTH)
    df = df.loc[df["three_axis_common"].astype(bool)].copy()
    assert len(df) == 238

    probe = pd.read_csv(PROBE)
    if probe["mof_id"].duplicated().any():
        raise RuntimeError("Probe table contains duplicate mof_id rows.")

    probe_cols = [c for c in probe.columns if c.startswith("probe_")]
    if len(probe_cols) != 12:
        raise RuntimeError(f"Frozen Gate C requires exactly 12 probe features; found {len(probe_cols)}")

    manifest = pd.read_csv(MANIFEST)
    geom_cols = manifest.loc[manifest["feature_set"].eq("geometry"), "column"].tolist()
    if len(geom_cols) != 14:
        raise RuntimeError(f"Expected 14 frozen geometry descriptors, found {len(geom_cols)}")

    z = df.merge(probe, on="mof_id", how="left", validate="one_to_one")
    if z[probe_cols].isna().all(axis=1).any():
        miss = z.loc[z[probe_cols].isna().all(axis=1), "mof_id"].tolist()
        raise RuntimeError(f"Missing probe descriptors for {len(miss)} MOFs: {miss[:10]}")

    target_col = "delta_lnS_elec_final" if "delta_lnS_elec_final" in z.columns else "delta_lnS_elec"
    y = z[target_col].to_numpy(float)

    feature_sets = {
        "geometry_only": geom_cols,
        "probe_only": probe_cols,
        "geometry_plus_probe": geom_cols + probe_cols,
    }

    preds = {}
    rows = []
    for name, cols in feature_sets.items():
        yp = cv_predict(z[cols], y)
        preds[name] = yp
        m = metrics(y, yp)
        rows.append({"feature_set": name, "n_features": len(cols), **m})

    met = pd.DataFrame(rows).set_index("feature_set")
    geom = met.loc["geometry_only"]
    probe_m = met.loc["probe_only"]
    comb = met.loc["geometry_plus_probe"]

    rng = np.random.default_rng(BOOT_SEED)
    mae_gain_comb = np.empty(N_BOOT)
    mae_gain_probe = np.empty(N_BOOT)
    r2_gain_comb = np.empty(N_BOOT)

    n = len(y)
    for b in range(N_BOOT):
        ix = rng.integers(0, n, size=n)
        yy = y[ix]
        pg = preds["geometry_only"][ix]
        pp = preds["probe_only"][ix]
        pc = preds["geometry_plus_probe"][ix]

        mae_gain_comb[b] = mean_absolute_error(yy, pg) - mean_absolute_error(yy, pc)
        mae_gain_probe[b] = mean_absolute_error(yy, pg) - mean_absolute_error(yy, pp)
        r2_gain_comb[b] = r2_score(yy, pc) - r2_score(yy, pg)

    combined_r2_gain = float(comb.r2 - geom.r2)
    combined_mae_ratio = float(comb.mae / geom.mae)
    probe_r2_gain = float(probe_m.r2 - geom.r2)
    probe_mae_ratio = float(probe_m.mae / geom.mae)
    comb_mae_gain_ci = ci95(mae_gain_comb)

    conditions = {
        "combined_r2_gain_ge_0p10": combined_r2_gain >= 0.10,
        "combined_mae_ratio_le_0p90": combined_mae_ratio <= 0.90,
        "combined_mae_gain_ci95_low_gt_0": comb_mae_gain_ci[0] > 0.0,
        "probe_only_r2_gain_ge_0p05": probe_r2_gain >= 0.05,
        "probe_only_mae_ratio_le_0p95": probe_mae_ratio <= 0.95,
    }
    passed = all(conditions.values())

    met.reset_index().to_csv(OUT / "pi_gate_c_probe_model_metrics.csv", index=False)

    pred_out = pd.DataFrame({
        "mof_id": z["mof_id"],
        "observed_delta_elec": y,
        "pred_geometry": preds["geometry_only"],
        "pred_probe": preds["probe_only"],
        "pred_geometry_plus_probe": preds["geometry_plus_probe"],
    })
    pred_out.to_csv(OUT / "pi_gate_c_probe_oof_predictions.csv", index=False)

    result = {
        "status": "MAIN_TEXT_PASS" if passed else "SI_ONLY_CLOSE_GATE_C",
        "n": n,
        "probe_feature_count": len(probe_cols),
        "geometry_feature_count": len(geom_cols),
        "metrics": met.to_dict(orient="index"),
        "decision_statistics": {
            "combined_r2_gain": combined_r2_gain,
            "combined_mae_ratio": combined_mae_ratio,
            "probe_only_r2_gain": probe_r2_gain,
            "probe_only_mae_ratio": probe_mae_ratio,
            "combined_mae_gain_bootstrap_ci95": comb_mae_gain_ci,
            "probe_mae_gain_bootstrap_ci95": ci95(mae_gain_probe),
            "combined_r2_gain_bootstrap_ci95": ci95(r2_gain_comb),
        },
        "conditions": conditions,
        "hard_stop_after_gate_c": True,
        "next_action": (
            "Integrate probe coordinate into Figure 3 and main text, then close paper."
            if passed else
            "Record Gate C in SI without strengthening the main claim, then close paper."
        ),
        "do_not_auto_open": [
            "alternative force-field campaign",
            "CO2/CH4 transfer",
            "GNNs",
            "larger descriptor search",
            "finite-pressure simulation",
        ],
    }

    (DOCS / "PI_GATE_C_RESULT.json").write_text(
        json.dumps(result, indent=2) + "\n", encoding="utf-8"
    )

    print("=" * 88)
    print("PI GATE C — FROZEN DECISION TEST")
    print("=" * 88)
    print(pd.DataFrame(rows).to_string(index=False))
    print()
    print(json.dumps(result, indent=2))
    print()
    print("GATE C CLOSED:", result["status"])

if __name__ == "__main__":
    main()
