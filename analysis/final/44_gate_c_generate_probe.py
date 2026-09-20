#!/usr/bin/env python3
"""
PI Gate C — final bounded electrostatic-probe generator.

This is the only production generator for Gate C.

Run with the validated RASPALIB 3.0.29 environment, for example:
    conda run -n raspa_maps python .\analysis\final\44_gate_c_generate_probe.py --pilot

After the pilot passes:
    conda run -n raspa_maps python .\analysis\final\44_gate_c_generate_probe.py --all

The script is resumable. Per-framework checkpoints are stored under
work/gate_c_probe/<accepted_positions>/ and are not part of the manuscript.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import math
import os
import sys
from concurrent.futures import ProcessPoolExecutor, as_completed
from pathlib import Path

import csv
import numpy as np

ROOT = Path(__file__).resolve().parents[2]
PROJECT = ROOT.parent
BACKEND_DIR = PROJECT / "path_b_analysis" / "paper_closure" / "02_spatial_maps" / "scripts" / "spatial_maps"
BACKEND_FILE = BACKEND_DIR / "fixed_position_energy.py"

MANIFEST = ROOT / "data" / "source_data" / "pi_gate_c_template_manifest.csv"
PREFLIGHT = ROOT / "docs" / "PI_GATE_C_PREFLIGHT_REPORT.json"
PROTOCOL = ROOT / "reproducibility" / "pi_gate_c_frozen_protocol.json"
OUT = ROOT / "data" / "source_data"
DOCS = ROOT / "docs"
WORK = ROOT / "work" / "gate_c_probe"

BASE_SEED = 20260920
DEFAULT_ACCEPTED = 6000
FALLBACK_ACCEPTED = 10000
N_ORIENT = 16
PILOT_ALIASES = ["M001", "M029", "M050", "M085", "M120", "M160", "M200", "M238"]

EXPECTED_PROBE_COLUMNS = [
    "probe_co2_median_z",
    "probe_co2_q90_z",
    "probe_co2_iqr_z",
    "probe_co2_fraction_z_gt_1",
    "probe_co2_median_orientation_iqr_z",
    "probe_n2_median_z",
    "probe_n2_q90_z",
    "probe_n2_iqr_z",
    "probe_n2_fraction_z_gt_1",
    "probe_n2_median_orientation_iqr_z",
    "probe_co2_minus_n2_median_z",
    "probe_co2_minus_n2_q90_z",
]

def logsumexp(a, axis=1):
    a = np.asarray(a, dtype=float)
    m = np.max(a, axis=axis, keepdims=True)
    out = m + np.log(np.sum(np.exp(a - m), axis=axis, keepdims=True))
    return np.squeeze(out, axis=axis)

def _rankdata(x):
    x = np.asarray(x, dtype=float)
    order = np.argsort(x, kind="mergesort")
    ranks = np.empty(len(x), dtype=float)
    i = 0
    while i < len(x):
        j = i + 1
        while j < len(x) and x[order[j]] == x[order[i]]:
            j += 1
        rank = 0.5 * ((i + 1) + j)
        ranks[order[i:j]] = rank
        i = j
    return ranks

def spearman_rho(x, y):
    rx = _rankdata(x)
    ry = _rankdata(y)
    if np.std(rx) == 0 or np.std(ry) == 0:
        return float("nan")
    return float(np.corrcoef(rx, ry)[0, 1])

def sha256(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as fh:
        for chunk in iter(lambda: fh.read(1024 * 1024), b""):
            h.update(chunk)
    return h.hexdigest()

def fibonacci_axes(n=N_ORIENT):
    i = np.arange(n, dtype=float)
    z = 1.0 - 2.0 * (i + 0.5) / n
    phi = (math.pi * (3.0 - math.sqrt(5.0))) * i
    r = np.sqrt(np.maximum(0.0, 1.0 - z * z))
    a = np.column_stack([r * np.cos(phi), r * np.sin(phi), z])
    a /= np.linalg.norm(a, axis=1)[:, None]
    return a

AXES = fibonacci_axes()

def weighted_quantile(x, w, q):
    x = np.asarray(x, float)
    w = np.asarray(w, float)
    ok = np.isfinite(x) & np.isfinite(w) & (w > 0)
    x = x[ok]
    w = w[ok]
    if len(x) == 0:
        raise RuntimeError("No finite positive-weight observations.")
    order = np.argsort(x)
    x = x[order]
    w = w[order]
    c = np.cumsum(w)
    c /= c[-1]
    return float(np.interp(q, c, x))

def stable_position_stats(log_lj, log_full):
    log_lj = np.asarray(log_lj, float)
    log_full = np.asarray(log_full, float)
    if log_lj.shape != log_full.shape or log_lj.ndim != 2 or log_lj.shape[1] != N_ORIENT:
        raise RuntimeError("Unexpected orientation matrix shape.")
    if not np.isfinite(log_lj).all() or not np.isfinite(log_full).all():
        raise RuntimeError("Non-finite pose log weights.")

    log_pos_lj = logsumexp(log_lj, axis=1) - math.log(N_ORIENT)
    log_pos_full = logsumexp(log_full, axis=1) - math.log(N_ORIENT)
    z_eff = log_pos_full - log_pos_lj
    z_pose = log_full - log_lj
    ori_iqr = np.quantile(z_pose, 0.75, axis=1) - np.quantile(z_pose, 0.25, axis=1)

    m = float(np.max(log_pos_lj))
    w = np.exp(log_pos_lj - m)
    if not np.isfinite(w).all() or not np.any(w > 0):
        raise RuntimeError("Invalid LJ statistical weights.")

    return z_eff, ori_iqr, w, log_pos_lj

def summarize_gas(log_lj, log_full):
    z, oi, w, log_pos_lj = stable_position_stats(log_lj, log_full)
    q25 = weighted_quantile(z, w, 0.25)
    q50 = weighted_quantile(z, w, 0.50)
    q75 = weighted_quantile(z, w, 0.75)
    q90 = weighted_quantile(z, w, 0.90)
    frac = float(np.sum(w[z > 1.0]) / np.sum(w))
    med_oi = weighted_quantile(oi, w, 0.50)
    ess = float(np.sum(w) ** 2 / np.sum(w ** 2))
    return {
        "median_z": q50,
        "q90_z": q90,
        "iqr_z": q75 - q25,
        "fraction_z_gt_1": frac,
        "median_orientation_iqr_z": med_oi,
        "effective_position_sample_size": ess,
        "max_log_lj_position_weight": float(np.max(log_pos_lj)),
    }

def descriptors_from_matrices(mats, n):
    out = {}
    gas_summ = {}
    for gas in ("CO2", "N2"):
        s = summarize_gas(mats[gas]["lj"][:n], mats[gas]["full"][:n])
        gas_summ[gas] = s
        p = gas.lower()
        out[f"probe_{p}_median_z"] = s["median_z"]
        out[f"probe_{p}_q90_z"] = s["q90_z"]
        out[f"probe_{p}_iqr_z"] = s["iqr_z"]
        out[f"probe_{p}_fraction_z_gt_1"] = s["fraction_z_gt_1"]
        out[f"probe_{p}_median_orientation_iqr_z"] = s["median_orientation_iqr_z"]

    out["probe_co2_minus_n2_median_z"] = gas_summ["CO2"]["median_z"] - gas_summ["N2"]["median_z"]
    out["probe_co2_minus_n2_q90_z"] = gas_summ["CO2"]["q90_z"] - gas_summ["N2"]["q90_z"]
    assert set(EXPECTED_PROBE_COLUMNS) == set(out)
    return out, gas_summ

def load_backend():
    """
    Load the frozen validated fixed-position backend without modifying it.

    The retained backend originally hard-coded only the elements present in the
    M029/M085 validation cases. Gate C spans all 238 MOFs, so this loader makes
    one input-translation-only compatibility extension: the local element-to-
    atomic-number table is expanded to the full periodic table before execution.
    No force-field parameters, charges, geometry, Ewald settings, energies, or
    RASPALIB calls are changed.
    """
    if not BACKEND_FILE.exists():
        raise FileNotFoundError(f"Validated backend not found: {BACKEND_FILE}")

    source = BACKEND_FILE.read_text(encoding="utf-8")

    old = """        elements = {'H':1,'C':6,'N':7,'O':8,'Co':27,'Zn':30,'Cu':29,'Zr':40}
"""
    symbols = [
        "H","He","Li","Be","B","C","N","O","F","Ne",
        "Na","Mg","Al","Si","P","S","Cl","Ar","K","Ca",
        "Sc","Ti","V","Cr","Mn","Fe","Co","Ni","Cu","Zn",
        "Ga","Ge","As","Se","Br","Kr","Rb","Sr","Y","Zr",
        "Nb","Mo","Tc","Ru","Rh","Pd","Ag","Cd","In","Sn",
        "Sb","Te","I","Xe","Cs","Ba","La","Ce","Pr","Nd",
        "Pm","Sm","Eu","Gd","Tb","Dy","Ho","Er","Tm","Yb",
        "Lu","Hf","Ta","W","Re","Os","Ir","Pt","Au","Hg",
        "Tl","Pb","Bi","Po","At","Rn","Fr","Ra","Ac","Th",
        "Pa","U","Np","Pu","Am","Cm","Bk","Cf","Es","Fm",
        "Md","No","Lr","Rf","Db","Sg","Bh","Hs","Mt","Ds",
        "Rg","Cn","Nh","Fl","Mc","Lv","Ts","Og",
    ]
    mapping = {sym: i + 1 for i, sym in enumerate(symbols)}
    new = "        elements = " + repr(mapping) + "\n"

    if old not in source:
        raise RuntimeError(
            "Frozen backend element-map line changed unexpectedly; refusing to patch."
        )
    source = source.replace(old, new, 1)

    ns = {
        "__name__": "_gate_c_fixed_position_energy_compat",
        "__file__": str(BACKEND_FILE),
    }
    exec(compile(source, str(BACKEND_FILE), "exec"), ns, ns)
    return ns["FixedPositionEnergy"]

def expected_backend_sha():
    rep = json.loads(PREFLIGHT.read_text(encoding="utf-8"))
    for rec in rep["backend_candidates_snapshotted"]:
        if "fixed_position_energy" in rec["source"].lower():
            return rec["sha256"]
    raise RuntimeError("Preflight report does not contain fixed_position_energy backend.")

def seed_for_alias(alias):
    try:
        n = int(alias[1:])
    except Exception as e:
        raise ValueError(f"Unexpected template alias: {alias}") from e
    return BASE_SEED + n * 1009

def run_framework(row, accepted):
    FixedPositionEnergy = load_backend()

    alias = row["template_alias"]
    mof_id = row["mof_id"]
    rng = np.random.default_rng(seed_for_alias(alias))
    frac = rng.random((accepted, 3))

    mats = {}
    provenance = {}
    cell_ref = None
    temp_ref = None

    for gas in ("CO2", "N2"):
        lj = FixedPositionEnergy(mof=alias, gas=gas, mode="lj_only")
        fu = FixedPositionEnergy(mof=alias, gas=gas, mode="full")

        if not np.allclose(lj.cell, fu.cell, atol=1e-12, rtol=0):
            raise RuntimeError(f"{alias}/{gas}: full/LJ cell mismatch")
        if not np.isclose(lj.temperature, fu.temperature, atol=0, rtol=0):
            raise RuntimeError(f"{alias}/{gas}: full/LJ temperature mismatch")
        if cell_ref is None:
            cell_ref = np.array(lj.cell, float)
            temp_ref = float(lj.temperature)
        else:
            if not np.allclose(cell_ref, lj.cell, atol=1e-12, rtol=0):
                raise RuntimeError(f"{alias}: gas cell mismatch")
            if temp_ref != float(lj.temperature):
                raise RuntimeError(f"{alias}: gas temperature mismatch")

        coms = frac @ np.asarray(lj.cell, float)
        log_lj = np.empty((accepted, N_ORIENT), dtype=np.float64)
        log_fu = np.empty((accepted, N_ORIENT), dtype=np.float64)

        for i, com in enumerate(coms):
            log_lj[i] = lj.log_weights(com, AXES)
            log_fu[i] = fu.log_weights(com, AXES)

        mats[gas] = {"lj": log_lj, "full": log_fu}
        provenance[gas] = {"lj_only": lj.provenance, "full": fu.provenance}

    lower = 3000 if accepted == DEFAULT_ACCEPTED else DEFAULT_ACCEPTED
    upper_desc, upper_aux = descriptors_from_matrices(mats, accepted)
    lower_desc, lower_aux = descriptors_from_matrices(mats, lower)

    return {
        "mof_id": mof_id,
        "template_alias": alias,
        "accepted_positions": accepted,
        "nested_lower_positions": lower,
        "orientations": N_ORIENT,
        "seed": seed_for_alias(alias),
        "descriptors_upper": upper_desc,
        "descriptors_lower": lower_desc,
        "aux_upper": upper_aux,
        "aux_lower": lower_aux,
        "provenance": provenance,
    }

def checkpoint_dir(accepted):
    p = WORK / str(accepted)
    p.mkdir(parents=True, exist_ok=True)
    return p

def checkpoint_path(alias, accepted):
    return checkpoint_dir(accepted) / f"{alias}.json"

def load_or_run(rows, accepted, workers):
    done = {}
    todo = []
    for row in rows:
        cp = checkpoint_path(row["template_alias"], accepted)
        if cp.exists():
            rec = json.loads(cp.read_text(encoding="utf-8"))
            if rec.get("accepted_positions") == accepted and rec.get("orientations") == N_ORIENT:
                done[row["template_alias"]] = rec
                print(f"RESUME {row['template_alias']}", flush=True)
                continue
        todo.append(row)

    if todo:
        if workers <= 1:
            for k, row in enumerate(todo, 1):
                print(f"RUN {row['template_alias']} ({k}/{len(todo)})", flush=True)
                rec = run_framework(row, accepted)
                checkpoint_path(row["template_alias"], accepted).write_text(
                    json.dumps(rec, indent=2) + "\n", encoding="utf-8"
                )
                done[row["template_alias"]] = rec
        else:
            with ProcessPoolExecutor(max_workers=workers) as ex:
                futs = {ex.submit(run_framework, row, accepted): row for row in todo}
                for k, fut in enumerate(as_completed(futs), 1):
                    row = futs[fut]
                    rec = fut.result()
                    checkpoint_path(row["template_alias"], accepted).write_text(
                        json.dumps(rec, indent=2) + "\n", encoding="utf-8"
                    )
                    done[row["template_alias"]] = rec
                    print(f"DONE {row['template_alias']} ({k}/{len(todo)})", flush=True)

    return [done[row["template_alias"]] for row in rows]

def pilot_check(records, accepted):
    if [r["template_alias"] for r in records] != PILOT_ALIASES:
        raise RuntimeError("Pilot aliases/order changed.")

    upper = {r["template_alias"]: r["descriptors_upper"] for r in records}
    lower = {r["template_alias"]: r["descriptors_lower"] for r in records}

    standardized = {}
    all_desc_pass = True
    for c in EXPECTED_PROBE_COLUMNS:
        up = np.array([upper[a][c] for a in PILOT_ALIASES], dtype=float)
        lo = np.array([lower[a][c] for a in PILOT_ALIASES], dtype=float)
        scale = float(np.std(up, ddof=0))
        dif = np.abs(up - lo)
        if scale < 1e-12:
            vals = dif
            pass_c = bool(np.max(dif) <= 1e-12)
        else:
            vals = dif / scale
            pass_c = bool(np.max(vals) <= 0.10)
        standardized[c] = {
            "scale_upper": scale,
            "max_abs_standardized_change": float(np.max(vals)),
            "passed": pass_c,
        }
        all_desc_pass &= pass_c

    rho_med = spearman_rho(
        [lower[a]["probe_co2_minus_n2_median_z"] for a in PILOT_ALIASES],
        [upper[a]["probe_co2_minus_n2_median_z"] for a in PILOT_ALIASES],
    )
    rho_q90 = spearman_rho(
        [lower[a]["probe_co2_minus_n2_q90_z"] for a in PILOT_ALIASES],
        [upper[a]["probe_co2_minus_n2_q90_z"] for a in PILOT_ALIASES],
    )

    direction = {
        "median_gap_M085_gt_M029": bool(
            upper["M085"]["probe_co2_minus_n2_median_z"]
            > upper["M029"]["probe_co2_minus_n2_median_z"]
        ),
        "q90_gap_M085_gt_M029": bool(
            upper["M085"]["probe_co2_minus_n2_q90_z"]
            > upper["M029"]["probe_co2_minus_n2_q90_z"]
        ),
    }

    ordering_pass = bool(
        np.isfinite(rho_med) and np.isfinite(rho_q90)
        and rho_med >= 0.90 and rho_q90 >= 0.90
    )
    direction_pass = all(direction.values())
    passed = bool(all_desc_pass and ordering_pass and direction_pass)

    report = {
        "status": "PASS" if passed else "REVIEW",
        "accepted_positions_upper": accepted,
        "nested_lower_positions": 3000 if accepted == DEFAULT_ACCEPTED else DEFAULT_ACCEPTED,
        "pilot_aliases": PILOT_ALIASES,
        "descriptor_convergence": standardized,
        "cross_gas_ordering": {
            "median_gap_spearman": rho_med,
            "q90_gap_spearman": rho_q90,
            "passed": ordering_pass,
        },
        "M029_M085_direction": {**direction, "passed": direction_pass},
        "allowed_next_action": (
            "RUN_ALL"
            if passed else
            ("RERUN_PILOT_AT_10000_ONCE"
             if accepted == DEFAULT_ACCEPTED
             else "CLOSE_GATE_C_AS_NON_SUPPORTING_NUMERICAL_PREFLIGHT")
        ),
    }

    p = DOCS / f"PI_GATE_C_PILOT_{accepted}.json"
    p.write_text(json.dumps(report, indent=2) + "\n", encoding="utf-8")

    csv_path = OUT / f"pi_gate_c_pilot_descriptors_{accepted}.csv"
    with csv_path.open("w", newline="", encoding="utf-8") as fh:
        fields = ["template_alias"] + EXPECTED_PROBE_COLUMNS
        w = csv.DictWriter(fh, fieldnames=fields)
        w.writeheader()
        for a in PILOT_ALIASES:
            w.writerow({"template_alias": a, **upper[a]})

    print(json.dumps(report, indent=2))
    return passed, report

def build_final_table(records, accepted):
    rows = []
    for r in records:
        row = {
            "mof_id": r["mof_id"],
            "template_alias": r["template_alias"],
            "accepted_positions": accepted,
            "orientations": N_ORIENT,
            **r["descriptors_upper"],
            "diagnostic_co2_effective_position_sample_size":
                r["aux_upper"]["CO2"]["effective_position_sample_size"],
            "diagnostic_n2_effective_position_sample_size":
                r["aux_upper"]["N2"]["effective_position_sample_size"],
        }
        rows.append(row)

    if len(rows) != 238:
        raise RuntimeError("Final probe table does not contain 238 frameworks.")
    ids = [r["mof_id"] for r in rows]
    if len(set(ids)) != len(ids):
        raise RuntimeError("Final probe table contains duplicate mof_id values.")
    if len([c for c in rows[0] if c.startswith("probe_")]) != 12:
        raise RuntimeError("Final probe feature count changed.")

    path = OUT / "pi_electrostatic_probe_descriptors.csv"
    with path.open("w", newline="", encoding="utf-8") as fh:
        w = csv.DictWriter(fh, fieldnames=list(rows[0]))
        w.writeheader()
        w.writerows(rows)

    rep = {
        "status": "PASS",
        "n": 238,
        "accepted_positions": accepted,
        "orientations": N_ORIENT,
        "feature_count": 12,
        "backend_file": str(BACKEND_FILE),
        "backend_sha256": sha256(BACKEND_FILE),
        "seed": BASE_SEED,
        "output": str(path.relative_to(ROOT)),
        "min_effective_position_sample_size_co2": float(min(
            r["diagnostic_co2_effective_position_sample_size"] for r in rows
        )),
        "min_effective_position_sample_size_n2": float(min(
            r["diagnostic_n2_effective_position_sample_size"] for r in rows
        )),
    }
    (DOCS / "PI_GATE_C_PROBE_GENERATION_REPORT.json").write_text(
        json.dumps(rep, indent=2) + "\n", encoding="utf-8"
    )
    print(json.dumps(rep, indent=2))

def main():
    ap = argparse.ArgumentParser()
    mode = ap.add_mutually_exclusive_group(required=True)
    mode.add_argument("--pilot", action="store_true")
    mode.add_argument("--all", action="store_true")
    ap.add_argument("--accepted", type=int, choices=[DEFAULT_ACCEPTED, FALLBACK_ACCEPTED], default=DEFAULT_ACCEPTED)
    ap.add_argument("--workers", type=int, default=2)
    args = ap.parse_args()

    os.environ.setdefault("OMP_NUM_THREADS", "1")

    for p in [MANIFEST, PREFLIGHT, PROTOCOL, BACKEND_FILE]:
        if not p.exists():
            raise FileNotFoundError(p)

    expected_sha = expected_backend_sha()
    actual_sha = sha256(BACKEND_FILE)
    if actual_sha != expected_sha:
        raise RuntimeError(
            "Validated fixed-position backend changed since preflight:\n"
            f"expected {expected_sha}\nactual   {actual_sha}"
        )

    try:
        import raspalib
    except Exception as e:
        raise RuntimeError(
            "RASPALIB is not importable in this Python. Run with the validated "
            "raspa_maps environment, for example: conda run -n raspa_maps python ..."
        ) from e

    with MANIFEST.open("r", newline="", encoding="utf-8") as fh:
        manifest_rows = list(csv.DictReader(fh))
    if len(manifest_rows) != 238 or not all(r["status"] == "MATCH" for r in manifest_rows):
        raise RuntimeError("Gate C template manifest is not the frozen 238/238 match.")

    records = [
        {"mof_id": r["mof_id"], "template_alias": r["template_alias"]}
        for r in manifest_rows
    ]

    if args.pilot:
        index = {r["template_alias"]: r for r in records}
        rows = [index[a] for a in PILOT_ALIASES]
        recs = load_or_run(rows, args.accepted, args.workers)
        passed, report = pilot_check(recs, args.accepted)
        if not passed:
            if args.accepted == DEFAULT_ACCEPTED:
                raise SystemExit(
                    "PILOT REVIEW — one allowed retry only:\n"
                    "conda run -n raspa_maps python .\\analysis\\final\\44_gate_c_generate_probe.py --pilot --accepted 10000"
                )
            raise SystemExit("PILOT REVIEW at 10,000 — Gate C must close as non-supporting; do not tune further.")
        print("PASS — Gate C numerical pilot passed. Production run is authorized.")
        return

    pilot_path = DOCS / f"PI_GATE_C_PILOT_{args.accepted}.json"
    if not pilot_path.exists():
        raise RuntimeError(f"Missing matching pilot report: {pilot_path}")
    pilot = json.loads(pilot_path.read_text(encoding="utf-8"))
    if pilot.get("status") != "PASS":
        raise RuntimeError("Matching Gate C pilot did not pass.")

    recs = load_or_run(records, args.accepted, args.workers)
    build_final_table(recs, args.accepted)
    print("PASS — 238-framework Gate C probe descriptor table complete.")
    print("Next: python .\\analysis\\final\\43_gate_c_evaluate_probe.py")

if __name__ == "__main__":
    main()
