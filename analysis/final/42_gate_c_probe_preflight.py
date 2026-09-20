#!/usr/bin/env python3
"""
PI Gate C preflight.

No simulation is run here.

This script:
1. verifies the 238-MOF publication cohort and production template package;
2. maps publication mof_id values to production template directories/CIFs;
3. searches the local parent project for the previously validated RASPALIB
   fixed-position energy evaluator used for the retained spatial maps;
4. snapshots the most likely backend scripts into the publication repository
   for audit/wiring;
5. writes a frozen Gate C protocol JSON and manifest.

Run from the publication repository with:
    python .\analysis\final\42_gate_c_probe_preflight.py
"""

from __future__ import annotations

import hashlib
import importlib.util
import json
import os
import re
import shutil
from pathlib import Path

import pandas as pd

HERE = Path(__file__).resolve()
REPO = HERE.parents[2]
PROJECT = REPO.parent
AUTH = REPO / "data" / "authoritative" / "path_b_final_238_authoritative.csv"
DIST = PROJECT / "distribution" / "ARCMOF_MECH_V1_PC2"
TEMPLATES = DIST / "templates"
OUT = REPO / "data" / "source_data"
SNAP = REPO / "reproducibility" / "gate_c_backend_snapshot"
DOCS = REPO / "docs"

SEED = 20260920
ACCEPTED_POSITIONS = 6000
ORIENTATIONS = 16

SKIP_PARTS = {
    ".git", "__pycache__", "node_modules", ".venv", "venv",
    "site-packages", "mof_selectivity_publication"
}
TOKENS = {
    "raspalib": 7,
    "field.npz": 5,
    "log_weight": 5,
    "relative_density": 3,
    "orientation": 3,
    "simulation_base.json": 3,
    "force_field.json": 3,
    "fixed": 2,
    "ewald": 2,
    "molecule": 1,
}

def sha256(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as fh:
        for chunk in iter(lambda: fh.read(1024 * 1024), b""):
            h.update(chunk)
    return h.hexdigest()

def safe_read(path: Path, max_bytes=1_500_000) -> str:
    if path.stat().st_size > max_bytes:
        return ""
    try:
        return path.read_text(encoding="utf-8", errors="ignore")
    except Exception:
        return ""

def score_backend(text: str) -> tuple[int, list[str]]:
    low = text.lower()
    hits = []
    score = 0
    for token, weight in TOKENS.items():
        if token in low:
            hits.append(token)
            score += weight
    return score, hits

def find_template_rows(ids: list[str]) -> pd.DataFrame:
    if not TEMPLATES.exists():
        raise FileNotFoundError(f"Missing production templates: {TEMPLATES}")

    rows = []
    dirs = sorted([p for p in TEMPLATES.iterdir() if p.is_dir()])
    by_id = {}

    for td in dirs:
        full = td / "full"
        cifs = sorted(full.glob("*_repeat.cif")) if full.exists() else []
        if not cifs:
            cifs = sorted(td.rglob("*_repeat.cif"))
        framework_id = None
        cif = None
        if cifs:
            cif = cifs[0]
            framework_id = re.sub(r"_repeat$", "", cif.stem)
            by_id.setdefault(framework_id, []).append((td, cif))

    for mid in ids:
        cand = by_id.get(mid, [])
        if len(cand) == 1:
            td, cif = cand[0]
            status = "MATCH"
        elif len(cand) == 0:
            td, cif, status = None, None, "MISSING"
        else:
            td, cif, status = cand[0][0], cand[0][1], "AMBIGUOUS"
        rows.append({
            "mof_id": mid,
            "template_alias": td.name if td else "",
            "template_dir": str(td) if td else "",
            "full_cif": str(cif) if cif else "",
            "status": status,
        })
    return pd.DataFrame(rows)

def find_backend_candidates() -> pd.DataFrame:
    rows = []
    for p in PROJECT.rglob("*.py"):
        rel_parts = set(p.relative_to(PROJECT).parts)
        if rel_parts & SKIP_PARTS:
            continue
        text = safe_read(p)
        if not text:
            continue
        score, hits = score_backend(text)
        if score >= 10:
            rows.append({
                "score": score,
                "path": str(p),
                "relative_path": str(p.relative_to(PROJECT)),
                "size_bytes": p.stat().st_size,
                "sha256": sha256(p),
                "hits": ";".join(hits),
            })
    if not rows:
        return pd.DataFrame(columns=[
            "score","path","relative_path","size_bytes","sha256","hits"
        ])
    return pd.DataFrame(rows).sort_values(
        ["score", "size_bytes"], ascending=[False, True]
    ).reset_index(drop=True)

def snapshot_candidates(cands: pd.DataFrame, n=5):
    # Windows-safe short snapshot names. The source path is preserved in the
    # report, so the copied filename does not need to encode the full tree.
    if SNAP.exists():
        shutil.rmtree(SNAP)
    SNAP.mkdir(parents=True, exist_ok=True)
    copied = []
    seen_sha = set()
    rank = 0
    for _, r in cands.iterrows():
        if len(copied) >= n:
            break
        if r["sha256"] in seen_sha:
            continue
        # Handoff-package copies are archival duplicates, not preferred backends.
        if "handoff_package" in str(r["relative_path"]).lower():
            continue
        seen_sha.add(r["sha256"])
        rank += 1
        src = Path(r["path"])
        stem = re.sub(r"[^A-Za-z0-9_.-]+", "_", src.stem)[:48]
        dst = SNAP / f"{rank:02d}_{stem}_{str(r['sha256'])[:10]}.py"
        shutil.copy2(src, dst)
        copied.append({
            "source": r["relative_path"],
            "snapshot": str(dst.relative_to(REPO)),
            "sha256": r["sha256"],
            "score": int(r["score"]),
        })
    return copied

def main():
    OUT.mkdir(parents=True, exist_ok=True)
    DOCS.mkdir(parents=True, exist_ok=True)

    if not AUTH.exists():
        raise FileNotFoundError(AUTH)

    df = pd.read_csv(AUTH)
    df = df.loc[df["three_axis_common"].astype(bool)].copy()
    if len(df) != 238:
        raise RuntimeError(f"Expected 238 common MOFs, found {len(df)}")

    ids = df["mof_id"].astype(str).tolist()
    manifest = find_template_rows(ids)
    manifest_path = OUT / "pi_gate_c_template_manifest.csv"
    manifest.to_csv(manifest_path, index=False)

    n_match = int((manifest["status"] == "MATCH").sum())
    n_missing = int((manifest["status"] == "MISSING").sum())
    n_amb = int((manifest["status"] == "AMBIGUOUS").sum())

    cands = find_backend_candidates()
    cand_path = OUT / "pi_gate_c_backend_candidates.csv"
    cands.to_csv(cand_path, index=False)
    copied = snapshot_candidates(cands)

    raspalib_present = importlib.util.find_spec("raspalib") is not None

    protocol = {
        "status": "GATE_C_PREFLIGHT",
        "central_claim_changed": False,
        "new_scientific_scope": "electrostatic probe only",
        "population_n": 238,
        "seed": SEED,
        "accepted_positions_per_framework_per_gas": ACCEPTED_POSITIONS,
        "orientations_per_position": ORIENTATIONS,
        "gases": ["CO2", "N2"],
        "pose_evaluation": [
            "LJ-only and full evaluated at the identical COM position and orientation",
            "z_elec = -beta * (U_full - U_LJ)",
            "position accepted when orientation-averaged LJ Boltzmann weight is finite and >0",
            "position weights use orientation-averaged LJ Boltzmann weight",
        ],
        "probe_features_per_gas": [
            "weighted_median_z_elec",
            "weighted_q90_z_elec",
            "weighted_iqr_z_elec",
            "weighted_fraction_z_elec_gt_1",
            "weighted_median_orientation_iqr_z_elec",
        ],
        "cross_gas_features": [
            "CO2_minus_N2_weighted_median_z_elec",
            "CO2_minus_N2_weighted_q90_z_elec",
        ],
        "feature_count": 12,
        "model_target": "delta_lnS_elec_final",
        "model_comparisons": ["geometry_only", "probe_only", "geometry_plus_probe"],
        "cv": {
            "type": "KFold",
            "n_splits": 5,
            "shuffle": True,
            "random_state": 20260918,
        },
        "extra_trees": {
            "n_estimators": 500,
            "min_samples_leaf": 3,
            "max_features": 0.7,
            "random_state": 20260918,
        },
        "main_text_pass_rule": {
            "combined_r2_gain_min": 0.10,
            "combined_mae_ratio_max": 0.90,
            "combined_mae_gain_bootstrap_ci95_low_gt": 0.0,
            "probe_only_r2_gain_min": 0.05,
            "probe_only_mae_ratio_max": 0.95,
        },
        "hard_stop_after_gate_c": True,
        "do_not_auto_open": [
            "alternative force-field campaign",
            "CO2/CH4 transfer",
            "GNNs",
            "feature expansion",
            "finite-pressure simulations",
        ],
    }
    protocol_path = REPO / "reproducibility" / "pi_gate_c_frozen_protocol.json"
    protocol_path.parent.mkdir(parents=True, exist_ok=True)
    protocol_path.write_text(json.dumps(protocol, indent=2) + "\n", encoding="utf-8")

    report = {
        "status": "PASS" if (n_match == 238 and len(cands) > 0) else "REVIEW",
        "publication_population": len(df),
        "template_matches": n_match,
        "template_missing": n_missing,
        "template_ambiguous": n_amb,
        "distribution_exists": DIST.exists(),
        "raspalib_importable_in_current_python": raspalib_present,
        "backend_candidate_count": int(len(cands)),
        "backend_candidates_snapshotted": copied,
        "next_action": (
            "Wire the highest-scoring validated fixed-position backend to the "
            "frozen Gate C protocol; do not change probe definitions."
        ),
    }
    report_path = DOCS / "PI_GATE_C_PREFLIGHT_REPORT.json"
    report_path.write_text(json.dumps(report, indent=2) + "\n", encoding="utf-8")

    print("=" * 92)
    print("PI GATE C — FINAL BOUNDED PROBE PREFLIGHT")
    print("=" * 92)
    print(json.dumps(report, indent=2))
    print()
    if len(cands):
        print("Top backend candidates:")
        print(cands.head(10).to_string(index=False))
    print()
    print("Saved:")
    print(" ", manifest_path)
    print(" ", cand_path)
    print(" ", protocol_path)
    print(" ", report_path)
    print()
    if n_match != 238:
        raise SystemExit("REVIEW: production template mapping is incomplete.")
    if len(cands) == 0:
        raise SystemExit("REVIEW: no local validated fixed-position backend was discovered.")
    print("PASS — Gate C scope and local backend evidence are frozen.")

if __name__ == "__main__":
    main()
