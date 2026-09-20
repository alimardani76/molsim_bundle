from __future__ import annotations

import json
import math
import subprocess
import sys
import time
from pathlib import Path

HERE = Path(__file__).resolve()
ROOT = HERE.parents[1]
CONDA = Path.home() / "miniforge3" / "condabin" / "conda.bat"

if "--inside-raspa-maps" not in sys.argv:
    cmd = [
        "cmd.exe", "/d", "/c", str(CONDA),
        "run", "-n", "raspa_maps", "--no-capture-output",
        "python", "-u", str(HERE), "--inside-raspa-maps",
        *sys.argv[1:],
    ]
    raise SystemExit(subprocess.run(cmd).returncode)

import numpy as np
from scipy.special import logsumexp, xlogy
from scipy.stats import spearmanr

SPATIAL = ROOT / "path_b_analysis" / "paper_closure" / "02_spatial_maps"
SCRIPTS = SPATIAL / "scripts" / "spatial_maps"
sys.path.insert(0, str(SCRIPTS))

from fixed_position_energy import FixedPositionEnergy
from build_spatial_field import orientation_axes, average_orientation

MOF = "M089"
GAS = "CO2"
MODE = "lj_only"
SPACING = 0.75
TEST_ORIENT = 16384
REF_ORIENT = 32768

REF_DIR = SPATIAL / "fields" / f"{MOF}_{GAS}_{MODE}_0.75A_{REF_ORIENT}orient"
TEST_DIR = SPATIAL / "fields" / f"{MOF}_{GAS}_{MODE}_0.75A_{TEST_ORIENT}orient"
PANEL = SPATIAL / "results" / "panel_expansion" / f"{MOF}_{GAS}_{MODE}"

# Deliberately strict because this is an orientation-count equivalence test,
# not a grid-convergence test.
MAX_MEAN_DIFF = 0.01
MAX_JS_BITS = 0.005
MAX_TV = 0.03
MIN_SUPPORT_SPEARMAN = 0.995
MAX_WEIGHTED_LOG_MAE = 0.05
SUPPORT_MASS = 0.999


def dump(path, obj):
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(obj, indent=2, allow_nan=False) + "\n", encoding="utf-8")


def normalize(x):
    x = np.asarray(x, dtype=float)
    s = float(x.sum())
    if not np.isfinite(x).all() or np.any(x < 0) or not s > 0:
        raise RuntimeError("Invalid probability weights.")
    return x / s


def js_bits(p, q):
    p = normalize(p)
    q = normalize(q)
    m = 0.5 * (p + q)
    val = 0.5 * np.sum(xlogy(p, p) - xlogy(p, m))
    val += 0.5 * np.sum(xlogy(q, q) - xlogy(q, m))
    return float(val / math.log(2.0))


def support_mask(p, q):
    masks = []
    for x in (p, q):
        order = np.argsort(x)[::-1]
        cs = np.cumsum(x[order])
        k = int(np.searchsorted(cs, SUPPORT_MASS, side="left")) + 1
        m = np.zeros(len(x), dtype=bool)
        m[order[:k]] = True
        masks.append(m)
    return masks[0] | masks[1]


def main():
    if not (REF_DIR / "field.npz").exists():
        raise SystemExit(f"Missing completed 32k reference field:\n{REF_DIR}")

    ref_meta = json.loads((REF_DIR / "metadata.json").read_text(encoding="utf-8"))
    ref = np.load(REF_DIR / "field.npz")

    e = FixedPositionEnergy(mof=MOF, gas=GAS, mode=MODE)
    dims = np.asarray(ref_meta["shape"], dtype=int)
    n = int(np.prod(dims))

    frac = np.asarray(ref["fractional_coordinates"], dtype=float)
    if frac.shape != (n, 3):
        raise RuntimeError("Reference fractional grid shape mismatch.")

    # Exact same spatial grid as the 32k reference.
    com = frac @ e.unit_cell
    axes = orientation_axes(TEST_ORIENT)

    TEST_DIR.mkdir(parents=True, exist_ok=True)
    checkpoint = TEST_DIR / "checkpoint.npz"

    log_w = np.full(n, np.nan)
    rel_se = np.full(n, np.nan)

    if checkpoint.exists():
        saved = np.load(checkpoint)
        if saved["log_weight"].shape != (n,):
            raise SystemExit("Checkpoint shape mismatch.")
        log_w[:] = saved["log_weight"]
        rel_se[:] = saved["orientation_relative_se"]
        print(f"Resuming {np.isfinite(log_w).sum()}/{n} completed points.")

    print("=" * 92)
    print("M089 0.75 A FIELD-LEVEL ORIENTATION EQUIVALENCE")
    print(f"Test: {TEST_ORIENT:,} orientations")
    print(f"Reference: {REF_ORIENT:,} orientations")
    print(f"Same grid points: {n:,}")
    print("=" * 92)

    start = time.perf_counter()
    for i, point in enumerate(com):
        if np.isfinite(log_w[i]):
            continue
        log_w[i], rel_se[i], _ = average_orientation(e, point, axes)
        if (i + 1) % 25 == 0 or i + 1 == n:
            np.savez_compressed(
                checkpoint,
                log_weight=log_w,
                orientation_relative_se=rel_se,
            )
            print(
                f"{np.isfinite(log_w).sum():,}/{n:,} "
                f"elapsed={time.perf_counter()-start:.1f}s",
                flush=True,
            )

    ref_log = np.asarray(ref["log_weight"], dtype=float)
    if not np.isfinite(log_w).all() or not np.isfinite(ref_log).all():
        raise RuntimeError("Non-finite log weights.")

    test_logmean = float(logsumexp(log_w) - math.log(n))
    ref_logmean = float(logsumexp(ref_log) - math.log(n))
    test_mean = float(np.exp(test_logmean))
    ref_mean = float(np.exp(ref_logmean))
    mean_diff = test_mean / ref_mean - 1.0

    p = normalize(np.exp(log_w - logsumexp(log_w)))
    q = normalize(np.exp(ref_log - logsumexp(ref_log)))

    js = js_bits(p, q)
    tv = float(0.5 * np.abs(p - q).sum())

    support = support_mask(p, q)
    tiny = np.finfo(float).tiny
    lp = np.log(np.maximum(p[support], tiny))
    lq = np.log(np.maximum(q[support], tiny))
    spear = float(spearmanr(lp, lq).statistic)

    w = 0.5 * (p[support] + q[support])
    w /= w.sum()
    weighted_log_mae = float(np.sum(w * np.abs(lp - lq)))

    same_peak = int(np.argmax(p)) == int(np.argmax(q))
    mass_weighted_se = float(
        np.sum((p / p.sum()) * rel_se)
    )

    passed = (
        abs(mean_diff) < MAX_MEAN_DIFF
        and js < MAX_JS_BITS
        and tv < MAX_TV
        and spear > MIN_SUPPORT_SPEARMAN
        and weighted_log_mae < MAX_WEIGHTED_LOG_MAE
    )

    relative = np.exp(log_w - test_logmean)
    raw = np.exp(log_w)

    np.savez_compressed(
        TEST_DIR / "field.npz",
        fractional_coordinates=frac,
        cartesian_coordinates_A=com,
        unit_cell_vectors_A=e.unit_cell,
        shape=dims,
        log_weight=log_w,
        weight=raw,
        relative_density=relative,
        orientation_relative_se=rel_se,
    )

    metadata = {
        "mof": MOF,
        "gas": GAS,
        "mode": MODE,
        "orientations": TEST_ORIENT,
        "shape": dims.tolist(),
        "grid_points": n,
        "requested_spacing_A": SPACING,
        "cell_mean_weight": test_mean,
        "reference_32768_mean_weight": ref_mean,
        "relative_difference_reference": mean_diff,
        "mass_weighted_orientation_relative_se": mass_weighted_se,
        "purpose": "field-level orientation equivalence against completed 32768-orientation field",
    }
    dump(TEST_DIR / "metadata.json", metadata)

    report = {
        "target": {"mof": MOF, "gas": GAS, "mode": MODE},
        "test_orientations": TEST_ORIENT,
        "reference_orientations": REF_ORIENT,
        "same_spatial_grid": True,
        "metrics": {
            "mean_weight_relative_difference": mean_diff,
            "js_divergence_bits": js,
            "total_variation_distance": tv,
            "support_log_density_spearman": spear,
            "probability_weighted_abs_log_error": weighted_log_mae,
            "same_peak_voxel": same_peak,
            "test_mass_weighted_orientation_relative_se": mass_weighted_se,
        },
        "gate": {
            "passed": passed,
            "requirements": {
                "abs_mean_weight_difference_lt": MAX_MEAN_DIFF,
                "js_bits_lt": MAX_JS_BITS,
                "total_variation_lt": MAX_TV,
                "support_spearman_gt": MIN_SUPPORT_SPEARMAN,
                "weighted_abs_log_error_lt": MAX_WEIGHTED_LOG_MAE,
            },
        },
        "interpretation": (
            "If passed, 16384 orientations are empirically equivalent to 32768 "
            "for the complete M089/CO2/LJ 0.75 A probability field and may be "
            "used for the next 0.5 A grid-refinement test. If failed, retain 32768."
        ),
    }

    out = PANEL / "M089_16k_vs_32k_FIELD_EQUIVALENCE.json"
    dump(out, report)

    print("")
    print(f"Mean-weight difference:       {mean_diff:+.4%}")
    print(f"JS divergence:                {js:.8f} bits")
    print(f"Total variation:              {tv:.8f}")
    print(f"Support Spearman:             {spear:.8f}")
    print(f"Weighted abs log error:       {weighted_log_mae:.8f}")
    print(f"Same peak voxel:              {same_peak}")
    print(f"Mass-weighted orient. SE:     {mass_weighted_se:.4%}")
    print("")
    print(
        "RESULT:",
        "PASS — 16,384 orientations validated for field refinement"
        if passed
        else "FAIL — retain 32,768 orientations",
    )
    print("Report:", out)
    print("=" * 92)


if __name__ == "__main__":
    main()
