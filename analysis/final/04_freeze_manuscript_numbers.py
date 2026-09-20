from __future__ import annotations

import json
from pathlib import Path
import numpy as np
import pandas as pd

ROOT = Path(r"C:\Users\KaraPardazesh\Desktop\bundle\v1\mof_selectivity_mechanisms")
CLOSURE = ROOT / "path_b_analysis" / "paper_closure" / "01_evidence_freeze"

AUTH = ROOT / "path_b_analysis" / "path_b_final_238_authoritative.csv"
DESC = ROOT / "path_b_analysis" / "leakage_free_descriptors" / "descriptor_hierarchy_metrics.csv"
WCA = ROOT / "path_b_analysis" / "forcefield_bridge" / "final_analysis" / "wca_robustness_statistics.csv"

OUTCSV = CLOSURE / "manuscript_numbers.csv"
REPORT = CLOSURE / "manuscript_numbers_report.txt"
SNAPSHOT = CLOSURE / "authoritative_238_columns.txt"

ATTR_PAIR = [
    "DB0-m15_o11_o75_f0_fsc",
    "DB0-m15_o14_o77_f0_fsc",
]

ELEC_PAIR = [
    "DB0-m15_o100_o150_f0_fsc",
    "DB0-m15_o11_o30_jea",
]

REQ = ["lnS_wca", "delta_lnS_attr", "delta_lnS_elec", "lnS_full"]


def add(rows, section, quantity, value, units="", source="", note=""):
    rows.append({
        "section": section,
        "quantity": quantity,
        "value": value,
        "units": units,
        "source": source,
        "note": note,
    })


def find_id_column(df: pd.DataFrame, ids: list[str]) -> str | None:
    for c in df.columns:
        if df[c].dtype != object:
            continue
        vals = set(df[c].dropna().astype(str))
        if all(x in vals for x in ids):
            return c
    return None


def rank_desc(s: pd.Series) -> pd.Series:
    return s.rank(method="min", ascending=False).astype(int)


def main():
    CLOSURE.mkdir(parents=True, exist_ok=True)

    if not AUTH.exists():
        raise SystemExit(f"Missing authoritative table: {AUTH}")

    d = pd.read_csv(AUTH)
    missing = [c for c in REQ if c not in d.columns]
    if missing:
        raise SystemExit("Missing required columns:\n" + "\n".join(missing))

    SNAPSHOT.write_text("\n".join(map(str, d.columns)), encoding="utf-8")

    rows = []

    # ------------------------------------------------------------
    # Core population
    # ------------------------------------------------------------
    add(rows, "population", "authoritative_n", int(len(d)), "MOFs", str(AUTH))

    # ------------------------------------------------------------
    # Core mechanism coordinates
    # ------------------------------------------------------------
    med_wca = float(d["lnS_wca"].median())
    med_attr = float(d["delta_lnS_attr"].median())
    med_elec = float(d["delta_lnS_elec"].median())
    rho = float(d[["delta_lnS_attr","delta_lnS_elec"]].corr(method="spearman").iloc[0,1])

    add(rows, "mechanism", "median_lnS_wca", med_wca, "", str(AUTH))
    add(rows, "mechanism", "median_S_wca", float(np.exp(med_wca)), "dimensionless", str(AUTH))
    add(rows, "mechanism", "fraction_WCA_positive", float((d["lnS_wca"] > 0).mean()), "fraction", str(AUTH))
    add(rows, "mechanism", "median_delta_lnS_attr", med_attr, "", str(AUTH))
    add(rows, "mechanism", "median_attr_multiplier", float(np.exp(med_attr)), "fold", str(AUTH))
    add(rows, "mechanism", "fraction_attr_positive", float((d["delta_lnS_attr"] > 0).mean()), "fraction", str(AUTH))
    add(rows, "mechanism", "median_delta_lnS_elec", med_elec, "", str(AUTH))
    add(rows, "mechanism", "median_elec_multiplier", float(np.exp(med_elec)), "fold", str(AUTH))
    add(rows, "mechanism", "fraction_elec_positive", float((d["delta_lnS_elec"] > 0).mean()), "fraction", str(AUTH))
    add(rows, "mechanism", "spearman_attr_elec", rho, "", str(AUTH))

    # telescoping residual
    lnS_lj = d["lnS_wca"] + d["delta_lnS_attr"]
    reconstructed_full = lnS_lj + d["delta_lnS_elec"]
    residual = reconstructed_full - d["lnS_full"]
    add(rows, "mechanism", "max_abs_telescoping_residual", float(residual.abs().max()), "", str(AUTH))

    # ------------------------------------------------------------
    # Ranking changes
    # ------------------------------------------------------------
    r_wca = rank_desc(d["lnS_wca"])
    r_lj = rank_desc(lnS_lj)
    r_full = rank_desc(d["lnS_full"])

    attr_promotion = r_wca - r_lj
    elec_promotion = r_lj - r_full

    add(rows, "ranking", "max_attraction_rank_promotion", int(attr_promotion.max()), "positions", str(AUTH))
    add(rows, "ranking", "max_electrostatic_rank_promotion", int(elec_promotion.max()), "positions", str(AUTH))
    add(rows, "ranking", "max_attraction_rank_demotion", int((-attr_promotion).max()), "positions", str(AUTH))
    add(rows, "ranking", "max_electrostatic_rank_demotion", int((-elec_promotion).max()), "positions", str(AUTH))

    # ------------------------------------------------------------
    # Matched pairs
    # ------------------------------------------------------------
    id_col = find_id_column(d, ATTR_PAIR + ELEC_PAIR)

    if id_col is None:
        add(rows, "matched_pairs", "id_column_status", "NOT_FOUND", "", str(AUTH),
            "Could not identify a string column containing all four frozen matched-pair IDs.")
    else:
        add(rows, "matched_pairs", "id_column", id_col, "", str(AUTH))

        x = d.set_index(id_col)

        a1, a2 = ATTR_PAIR
        e1, e2 = ELEC_PAIR

        attr_diff = abs(float(x.loc[a1, "delta_lnS_attr"]) - float(x.loc[a2, "delta_lnS_attr"]))
        elec_diff = abs(float(x.loc[e1, "delta_lnS_elec"]) - float(x.loc[e2, "delta_lnS_elec"]))

        attr_wca_diff = abs(float(x.loc[a1, "lnS_wca"]) - float(x.loc[a2, "lnS_wca"]))
        elec_wca_diff = abs(float(x.loc[e1, "lnS_wca"]) - float(x.loc[e2, "lnS_wca"]))

        add(rows, "matched_pairs", "attraction_pair_delta_lnS_attr_difference", attr_diff, "", str(AUTH),
            f"{a1} vs {a2}")
        add(rows, "matched_pairs", "attraction_pair_multiplier_contrast", float(np.exp(attr_diff)), "fold", str(AUTH),
            f"{a1} vs {a2}")
        add(rows, "matched_pairs", "attraction_pair_abs_lnS_wca_difference", attr_wca_diff, "", str(AUTH),
            f"{a1} vs {a2}")

        add(rows, "matched_pairs", "electrostatic_pair_delta_lnS_elec_difference", elec_diff, "", str(AUTH),
            f"{e1} vs {e2}")
        add(rows, "matched_pairs", "electrostatic_pair_multiplier_contrast", float(np.exp(elec_diff)), "fold", str(AUTH),
            f"{e1} vs {e2}")
        add(rows, "matched_pairs", "electrostatic_pair_abs_lnS_wca_difference", elec_wca_diff, "", str(AUTH),
            f"{e1} vs {e2}")

    # ------------------------------------------------------------
    # Descriptor metrics: save exact table + parse common fields
    # ------------------------------------------------------------
    desc_df = None
    if DESC.exists():
        desc_df = pd.read_csv(DESC)
        desc_df.to_csv(CLOSURE / "descriptor_hierarchy_metrics_frozen_copy.csv", index=False)

        # Record every row exactly in a machine-readable long-ish form
        for i, rec in desc_df.iterrows():
            for c in desc_df.columns:
                add(rows, "descriptor_metrics_raw", f"row{i}__{c}", rec[c], "", str(DESC))

    # ------------------------------------------------------------
    # WCA robustness: save exact table + raw rows
    # ------------------------------------------------------------
    wca_df = None
    if WCA.exists():
        wca_df = pd.read_csv(WCA)
        wca_df.to_csv(CLOSURE / "wca_robustness_statistics_frozen_copy.csv", index=False)

        for i, rec in wca_df.iterrows():
            for c in wca_df.columns:
                add(rows, "wca_bridge_raw", f"row{i}__{c}", rec[c], "", str(WCA))

    # ------------------------------------------------------------
    # Family leave-one-out status
    # ------------------------------------------------------------
    add(
        rows,
        "robustness",
        "family_leave_one_out_status",
        "NOT_COMPLETED",
        "",
        str(ROOT / "path_b_analysis" / "negative_controls" / "family_leave_one_out_results.csv"),
        "Original output is a 2-byte empty file; authoritative 238-MOF table contains no family columns used by the original script."
    )

    out = pd.DataFrame(rows)
    out.to_csv(OUTCSV, index=False)

    # ------------------------------------------------------------
    # Report
    # ------------------------------------------------------------
    lines = []
    lines.append("=" * 88)
    lines.append("MANUSCRIPT NUMBER FREEZE")
    lines.append("=" * 88)
    lines.append("")
    lines.append(f"Authoritative population: {len(d)}")
    lines.append("")
    lines.append("CORE RESULTS")
    lines.append(f"median lnS_wca              = {med_wca:.12g}")
    lines.append(f"median S_wca                = {np.exp(med_wca):.12g}")
    lines.append(f"fraction WCA positive       = {(d['lnS_wca'] > 0).mean():.12g}")
    lines.append(f"median delta_lnS_attr       = {med_attr:.12g}")
    lines.append(f"median attr multiplier      = {np.exp(med_attr):.12g}")
    lines.append(f"fraction attr positive      = {(d['delta_lnS_attr'] > 0).mean():.12g}")
    lines.append(f"median delta_lnS_elec       = {med_elec:.12g}")
    lines.append(f"median elec multiplier      = {np.exp(med_elec):.12g}")
    lines.append(f"fraction elec positive      = {(d['delta_lnS_elec'] > 0).mean():.12g}")
    lines.append(f"rho(attr,elec)              = {rho:.12g}")
    lines.append(f"max telescoping residual    = {residual.abs().max():.12g}")
    lines.append("")
    lines.append("RANKING")
    lines.append(f"max attraction promotion    = {int(attr_promotion.max())}")
    lines.append(f"max electrostatic promotion = {int(elec_promotion.max())}")
    lines.append(f"max attraction demotion     = {int((-attr_promotion).max())}")
    lines.append(f"max electrostatic demotion  = {int((-elec_promotion).max())}")
    lines.append("")
    lines.append("MATCHED PAIRS")
    if id_col is None:
        lines.append("ID column: NOT FOUND")
    else:
        lines.append(f"ID column: {id_col}")
        lines.append(f"attraction pair delta difference = {attr_diff:.12g}")
        lines.append(f"attraction pair multiplier       = {np.exp(attr_diff):.12g}")
        lines.append(f"attraction pair |delta lnS_wca|  = {attr_wca_diff:.12g}")
        lines.append(f"electrostatic pair delta difference = {elec_diff:.12g}")
        lines.append(f"electrostatic pair multiplier       = {np.exp(elec_diff):.12g}")
        lines.append(f"electrostatic pair |delta lnS_wca|  = {elec_wca_diff:.12g}")
    lines.append("")
    lines.append("DESCRIPTOR METRICS TABLE")
    if desc_df is None:
        lines.append("MISSING")
    else:
        lines.append(desc_df.to_string(index=False))
    lines.append("")
    lines.append("WCA ROBUSTNESS TABLE")
    if wca_df is None:
        lines.append("MISSING")
    else:
        lines.append(wca_df.to_string(index=False))
    lines.append("")
    lines.append("FAMILY LEAVE-ONE-OUT")
    lines.append("NOT COMPLETED — remove from final evidence set.")
    lines.append("")
    lines.append(f"Saved CSV: {OUTCSV}")

    REPORT.write_text("\n".join(lines), encoding="utf-8")
    print("\n".join(lines))
    print("")
    print(f"Saved: {REPORT}")


if __name__ == "__main__":
    main()
