#!/usr/bin/env python3
"""Check *archived* paper headline statistics, not raw-simulation reproducibility."""
from pathlib import Path
import numpy as np
import pandas as pd
ROOT = Path(__file__).resolve().parents[1]
p = ROOT / "data/authoritative/path_b_final_238_authoritative.csv"
df = pd.read_csv(p)
required = ["lnS_wca", "delta_lnS_attr", "delta_lnS_elec", "lnS_full", "lnS_lj_only"]
missing = sorted(set(required)-set(df.columns))
assert not missing, f"Missing columns: {missing}"
assert len(df) == 238, f"Expected 238 rows; got {len(df)}"
assert df[required].notna().all().all(), "Non-finite/missing numerical entries"
assert np.isfinite(df[required].to_numpy(dtype=float)).all()
assert np.allclose(df["lnS_wca"] + df["delta_lnS_attr"], df["lnS_lj_only"], atol=1e-8)
assert np.allclose(df["lnS_lj_only"] + df["delta_lnS_elec"], df["lnS_full"], atol=1e-8)
expected = {"lnS_wca": -0.203749485237, "delta_lnS_attr": 2.53737562365, "delta_lnS_elec": 0.309492036602}
for name, want in expected.items():
    got = df[name].median()
    assert abs(got-want) < 2e-6, (name, got, want)
    print(f"median {name}: {got:.12f}")
print("PASS: 238 records; finite results; telescoping identity; frozen headline medians.")
