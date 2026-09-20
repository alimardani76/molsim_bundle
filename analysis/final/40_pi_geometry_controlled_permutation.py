#!/usr/bin/env python3
"""
PI Gate B v2: balanced geometry-conditioned shuffle diagnostic.

WHY V2
------
K-means produced singleton geometry clusters in this dataset. Singleton strata
cannot provide a meaningful within-stratum shuffle and arbitrary post-hoc
merging would make the result depend on an ad hoc repair.

This version therefore uses geometry-only PCA followed by deterministic
equal-frequency stratification in PC1 and PC2. The strata are built entirely
from the 14 frozen pre-simulation geometry descriptors and contain roughly
equal numbers of frameworks.

Primary stratification: 5 x 2 = 10 strata.
Sensitivity:
    4 x 2 = 8 strata
    5 x 2 = 10 strata
    4 x 3 = 12 strata

Within each stratum, the observed target is shuffled 10,000 times while the
existing geometry-only OOF prediction is held fixed.

IMPORTANT INTERPRETATION
------------------------
This is a secondary geometry-conditioned SHUFFLE DIAGNOSTIC for the frozen OOF
geometry score. It is not presented as an exact conditional-randomization
p-value because the OOF scores were learned from the original dataset.

The main paper claim remains based on:
  * the observed OOF predictability gap,
  * its framework bootstrap,
  * outlier sensitivity,
  * geometry-matched counterexamples.

This diagnostic is retained only if its qualitative conclusion is stable for
8, 10, and 12 balanced strata.
"""

from __future__ import annotations

from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import spearmanr
from sklearn.decomposition import PCA
from sklearn.impute import SimpleImputer
from sklearn.metrics import mean_absolute_error, r2_score
from sklearn.preprocessing import StandardScaler

ROOT = Path(__file__).resolve().parents[2]
AUTH = ROOT / "data" / "authoritative" / "path_b_final_238_authoritative.csv"
PRED = ROOT / "data" / "source_data" / "descriptor_hierarchy_predictions.csv"
FEATURES = ROOT / "data" / "source_data" / "descriptor_feature_manifest.csv"
OUT = ROOT / "data" / "source_data"
DOCS = ROOT / "docs"

SEED = 20260920
N_PERM = 10_000
TARGETS = ["delta_attr", "delta_elec"]

STRATA_DESIGNS = {
    8: (4, 2),
    10: (5, 2),
    12: (4, 3),
}
PRIMARY_K = 10


def balanced_bins(values: np.ndarray, n_bins: int) -> np.ndarray:
    """Deterministic equal-frequency bins, robust to tied values."""
    order = np.argsort(values, kind="mergesort")
    labels = np.empty(len(values), dtype=int)
    for rank, idx in enumerate(order):
        labels[idx] = min(n_bins - 1, (rank * n_bins) // len(values))
    return labels


def make_balanced_2d_strata(pc1, pc2, n1, n2):
    b1 = balanced_bins(pc1, n1)
    strata = np.empty(len(pc1), dtype=int)

    for first_bin in range(n1):
        ix = np.flatnonzero(b1 == first_bin)
        local_b2 = balanced_bins(pc2[ix], n2)
        strata[ix] = first_bin * n2 + local_b2

    return strata


def within_stratum_permutation(values, strata, rng):
    out = np.asarray(values, dtype=float).copy()
    for s in np.unique(strata):
        ix = np.flatnonzero(strata == s)
        out[ix] = rng.permutation(out[ix])
    return out


def q(a, p):
    return float(np.quantile(np.asarray(a, dtype=float), p))


def upper_fraction(null, obs):
    null = np.asarray(null, dtype=float)
    return float((1 + np.sum(null >= obs)) / (len(null) + 1))


def lower_fraction(null, obs):
    null = np.asarray(null, dtype=float)
    return float((1 + np.sum(null <= obs)) / (len(null) + 1))


def main():
    OUT.mkdir(parents=True, exist_ok=True)
    DOCS.mkdir(parents=True, exist_ok=True)

    if not AUTH.exists():
        raise FileNotFoundError(f"Missing authoritative table: {AUTH}")
    if not PRED.exists():
        raise FileNotFoundError(f"Missing OOF predictions: {PRED}")
    if not FEATURES.exists():
        raise FileNotFoundError(f"Missing feature manifest: {FEATURES}")

    df = pd.read_csv(AUTH)
    df = df.loc[df["three_axis_common"].astype(bool)].copy()
    assert len(df) == 238, f"Expected 238 common frameworks, found {len(df)}"

    manifest = pd.read_csv(FEATURES)
    geom_cols = manifest.loc[
        manifest["feature_set"].eq("geometry"), "column"
    ].tolist()
    assert len(geom_cols) == 14, f"Expected 14 geometry descriptors, found {geom_cols}"

    X = df[geom_cols].copy()
    X = SimpleImputer(strategy="median").fit_transform(X)
    X = StandardScaler().fit_transform(X)

    pca = PCA(n_components=2, random_state=SEED)
    Z = pca.fit_transform(X)
    pc1, pc2 = Z[:, 0], Z[:, 1]

    pred = pd.read_csv(PRED)
    pred = pred.loc[
        pred["feature_set"].eq("geometry") &
        pred["target"].isin(TARGETS)
    ].copy()

    ids = df["mof_id"].tolist()
    observed = {}
    predicted = {}

    for target in TARGETS:
        p = pred.loc[pred["target"].eq(target)].set_index("mof_id")
        missing = sorted(set(ids) - set(p.index))
        if missing:
            raise RuntimeError(f"Missing OOF rows for {target}: {missing[:10]}")
        p = p.loc[ids]
        observed[target] = p["observed"].to_numpy(float)
        predicted[target] = p["predicted"].to_numpy(float)

    assignments = pd.DataFrame({
        "mof_id": ids,
        "geometry_PC1": pc1,
        "geometry_PC2": pc2,
    })

    summary_rows = []
    quantile_rows = []

    rng_master = np.random.default_rng(SEED)

    for nominal_k, (n1, n2) in STRATA_DESIGNS.items():
        strata = make_balanced_2d_strata(pc1, pc2, n1, n2)
        assignments[f"balanced_geometry_stratum_k{nominal_k}"] = strata

        sizes = pd.Series(strata).value_counts().sort_index()
        actual_k = int(sizes.size)

        if actual_k != nominal_k:
            raise RuntimeError(
                f"Expected {nominal_k} strata but obtained {actual_k}"
            )

        print(
            f"K={nominal_k}: min={int(sizes.min())}, "
            f"max={int(sizes.max())}, sizes={sizes.tolist()}"
        )

        for target in TARGETS:
            y = observed[target]
            yp = predicted[target]

            obs_rho = float(spearmanr(y, yp).statistic)
            obs_r2 = float(r2_score(y, yp))
            obs_mae = float(mean_absolute_error(y, yp))

            null_rho = np.empty(N_PERM, dtype=float)
            null_r2 = np.empty(N_PERM, dtype=float)
            null_mae = np.empty(N_PERM, dtype=float)

            local_seed = int(rng_master.integers(0, 2**32 - 1))
            rng = np.random.default_rng(local_seed)

            for b in range(N_PERM):
                yperm = within_stratum_permutation(y, strata, rng)
                null_rho[b] = float(spearmanr(yperm, yp).statistic)
                null_r2[b] = float(r2_score(yperm, yp))
                null_mae[b] = float(mean_absolute_error(yperm, yp))

            summary_rows.append({
                "nominal_k": nominal_k,
                "pc1_bins": n1,
                "pc2_bins_per_pc1": n2,
                "n": len(y),
                "n_permutations": N_PERM,
                "min_stratum_size": int(sizes.min()),
                "max_stratum_size": int(sizes.max()),
                "target": target,
                "observed_spearman_pred_vs_target": obs_rho,
                "null_spearman_median": float(np.median(null_rho)),
                "null_spearman_q025": q(null_rho, 0.025),
                "null_spearman_q975": q(null_rho, 0.975),
                "upper_exceedance_fraction_spearman": upper_fraction(null_rho, obs_rho),
                "observed_r2": obs_r2,
                "null_r2_median": float(np.median(null_r2)),
                "null_r2_q025": q(null_r2, 0.025),
                "null_r2_q975": q(null_r2, 0.975),
                "upper_exceedance_fraction_r2": upper_fraction(null_r2, obs_r2),
                "observed_mae": obs_mae,
                "null_mae_median": float(np.median(null_mae)),
                "null_mae_q025": q(null_mae, 0.025),
                "null_mae_q975": q(null_mae, 0.975),
                "lower_exceedance_fraction_mae": lower_fraction(null_mae, obs_mae),
            })

            for stat_name, arr in [
                ("spearman", null_rho),
                ("r2", null_r2),
                ("mae", null_mae),
            ]:
                quantile_rows.append({
                    "nominal_k": nominal_k,
                    "target": target,
                    "statistic": stat_name,
                    "q001": q(arr, 0.001),
                    "q01": q(arr, 0.01),
                    "q025": q(arr, 0.025),
                    "q50": q(arr, 0.50),
                    "q975": q(arr, 0.975),
                    "q99": q(arr, 0.99),
                    "q999": q(arr, 0.999),
                })

    summary = pd.DataFrame(summary_rows)
    summary_path = OUT / "pi_geometry_controlled_permutation_summary.csv"
    quantile_path = OUT / "pi_geometry_controlled_permutation_null_quantiles.csv"
    assignments_path = OUT / "pi_geometry_strata_assignments.csv"

    summary.to_csv(summary_path, index=False)
    pd.DataFrame(quantile_rows).to_csv(quantile_path, index=False)
    assignments.to_csv(assignments_path, index=False)

    primary = summary.loc[
        summary["nominal_k"].eq(PRIMARY_K)
    ].set_index("target")

    stable = True
    stability_lines = []
    for k in STRATA_DESIGNS:
        sub = summary.loc[summary["nominal_k"].eq(k)].set_index("target")
        attr_gap = (
            sub.loc["delta_attr", "observed_spearman_pred_vs_target"]
            - sub.loc["delta_attr", "null_spearman_q975"]
        )
        elec_gap = (
            sub.loc["delta_elec", "observed_spearman_pred_vs_target"]
            - sub.loc["delta_elec", "null_spearman_q975"]
        )
        contrast = attr_gap - elec_gap
        stability_lines.append(
            f"- K={k}: attraction excess-over-null={attr_gap:.3f}; "
            f"electrostatic excess-over-null={elec_gap:.3f}; "
            f"contrast={contrast:.3f}"
        )
        if contrast <= 0:
            stable = False

    report = f"""# PI Gate B — balanced geometry-conditioned shuffle diagnostic

Status: {"PASS" if stable else "REVIEW"}

## Why the original K-means version was rejected

The first implementation produced singleton geometry clusters. Rather than
merge those outliers after seeing the result, this version uses deterministic
equal-frequency strata in the first two principal components of the same
14 frozen geometry descriptors.

This keeps the diagnostic geometry-conditioned while ensuring every stratum
contains enough frameworks to shuffle.

## Geometry stratification

PCA explained-variance ratios:
- PC1: {pca.explained_variance_ratio_[0]:.4f}
- PC2: {pca.explained_variance_ratio_[1]:.4f}
- PC1+PC2: {pca.explained_variance_ratio_[:2].sum():.4f}

Sensitivity designs:
- K=8: 4 PC1 bins x 2 within-bin PC2 bins
- K=10: 5 x 2
- K=12: 4 x 3

Each target is shuffled 10,000 times within these geometry-only strata.

## Primary K=10

### Attractive increment

Observed Spearman(prediction, target):
**{primary.loc["delta_attr", "observed_spearman_pred_vs_target"]:.4f}**

Conditional-shuffle 95% interval:
**[{primary.loc["delta_attr", "null_spearman_q025"]:.4f},
{primary.loc["delta_attr", "null_spearman_q975"]:.4f}]**

Observed R2:
**{primary.loc["delta_attr", "observed_r2"]:.4f}**

Shuffle-null R2 95% interval:
**[{primary.loc["delta_attr", "null_r2_q025"]:.4f},
{primary.loc["delta_attr", "null_r2_q975"]:.4f}]**

### Electrostatic increment

Observed Spearman(prediction, target):
**{primary.loc["delta_elec", "observed_spearman_pred_vs_target"]:.4f}**

Conditional-shuffle 95% interval:
**[{primary.loc["delta_elec", "null_spearman_q025"]:.4f},
{primary.loc["delta_elec", "null_spearman_q975"]:.4f}]**

Observed R2:
**{primary.loc["delta_elec", "observed_r2"]:.4f}**

Shuffle-null R2 95% interval:
**[{primary.loc["delta_elec", "null_r2_q025"]:.4f},
{primary.loc["delta_elec", "null_r2_q975"]:.4f}]**

## Sensitivity contrast

{chr(10).join(stability_lines)}

## Interpretation boundary

This is a secondary shuffle diagnostic, not an exact conditional-randomization
p-value. The OOF geometry score was learned on the original dataset and is held
fixed during the shuffles.

Use it only if the same qualitative attraction-vs-electrostatics asymmetry
appears for K=8, K=10, and K=12.

The main claim remains supported independently by:
1. the observed OOF predictability gap;
2. the framework bootstrap of that gap;
3. outlier/leverage sensitivity;
4. geometry-matched mechanistic counterexamples.

Do not use this diagnostic to claim that electrostatics is independent of
geometry or that chemistry has been uniquely isolated.
"""

    report_path = DOCS / "PI_GATE_B_PERMUTATION_REPORT.md"
    report_path.write_text(report, encoding="utf-8")

    print("=" * 88)
    print("PI GATE B V2 — BALANCED GEOMETRY-CONDITIONED SHUFFLE")
    print("=" * 88)
    print()
    print(summary.to_string(index=False))
    print()
    print("PCA explained variance:")
    print(
        f"  PC1={pca.explained_variance_ratio_[0]:.4f}, "
        f"PC2={pca.explained_variance_ratio_[1]:.4f}, "
        f"sum={pca.explained_variance_ratio_[:2].sum():.4f}"
    )
    print()
    print("Qualitative sensitivity:")
    for line in stability_lines:
        print(line)
    print()
    print("FINAL GATE:", "PASS" if stable else "REVIEW")
    print()
    print("Saved:")
    print(" ", summary_path)
    print(" ", quantile_path)
    print(" ", assignments_path)
    print(" ", report_path)


if __name__ == "__main__":
    main()
