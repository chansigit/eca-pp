"""obs-profile tests (identify-columns spec §3, §10)."""

from __future__ import annotations

import anndata as ad
import numpy as np
import pandas as pd
import scipy.sparse as sp

from eca_pp.identify_columns import obsprofile

N = 120


def make_adata():
    rng = np.random.default_rng(0)
    A = ad.AnnData(X=sp.csr_matrix(rng.poisson(1.0, size=(N, 30)).astype(np.float32)))
    donor = np.repeat(["d1", "d2", "d3"], N // 3)
    lane = np.array([f"{d}_lane{i % 2}" for i, d in enumerate(donor)])
    A.obs_names = pd.Index([f"{d}-CELL{i:04d}" for i, d in enumerate(donor)],
                           dtype=object)
    A.obs["donor"] = donor
    A.obs["lane"] = lane
    A.obs["donor_name"] = np.char.upper(donor.astype(str))  # equivalent to donor
    A.obs["day"] = np.tile(["day1", "day2"], N // 2)        # orthogonal to donor
    A.obs["const"] = "x"
    A.obs["cell_id"] = [f"c{i}" for i in range(N)]
    A.obs["pct_junk"] = rng.random(N)
    return A


def _col(profile, name):
    return next(e for e in profile["columns"] if e["column"] == name)


def test_column_entries_and_flags():
    p = obsprofile.profile_obs(make_adata())
    assert p["n_obs"] == N
    assert _col(p, "const")["is_constant"]
    assert _col(p, "cell_id")["is_per_cell_unique"]
    assert _col(p, "pct_junk")["dtype"] == "float"
    donor = _col(p, "donor")
    assert donor["n_unique"] == 3
    assert donor["group_sizes"]["n_groups"] == 3
    assert donor["group_sizes"]["min"] == N // 3
    assert donor["group_sizes"]["n_tiny"] == 0
    assert 0.99 <= donor["entropy"] <= 1.0  # perfectly balanced
    assert donor["examples"]["d1"] == N // 3


def test_relations_nested_and_equivalent():
    p = obsprofile.profile_obs(make_adata())
    rel = {(r["finer"], r["coarser"], r["kind"]) for r in p["relations"]}
    assert ("lane", "donor", "nested") in rel
    assert any(k == "equivalent" and {f, c} == {"donor", "donor_name"}
               for f, c, k in rel)


def test_barcode_and_composite_candidates():
    p = obsprofile.profile_obs(make_adata())
    labels = {d["label"]: d for d in p["derived"]}
    assert labels["barcode:prefix:-"]["n_groups"] == 3
    assert set(labels["barcode:prefix:-"]["equivalent_with"]) >= {
        "donor", "donor_name"
    }
    # composite only appears for orthogonal pairs (donor x day); a pair where
    # one column already refines the other adds nothing and must NOT appear
    assert "composite:day+donor" in labels
    assert labels["composite:day+donor"]["n_groups"] == 6
    assert "composite:donor+lane" not in labels


def test_derive_values():
    A = make_adata()
    v = obsprofile.derive_values(A, "barcode:prefix:-")
    assert list(v.unique()) == ["d1", "d2", "d3"]
    c = obsprofile.derive_values(A, "composite:donor+lane")
    assert c.str.contains(r"\|").all()


def test_tiny_groups_flagged():
    A = make_adata()
    tiny = np.array([f"g{i}" for i in range(N // 4)] * 4)[:N]  # 30 groups of 4
    A.obs["micro"] = tiny
    p = obsprofile.profile_obs(A)
    gs = _col(p, "micro")["group_sizes"]
    assert gs["n_tiny"] == gs["n_groups"]  # every group < 25 cells
    assert gs["tiny_cell_frac"] == 1.0


def test_blank_and_common_missing_strings_are_not_groups():
    A = make_adata()
    values = np.array(["b1"] * 40 + ["b2"] * 40
                      + ["", "  ", "NA", "missing"] * 10, dtype=object)
    A.obs["batch_with_blanks"] = values
    p = obsprofile.profile_obs(A)
    col = _col(p, "batch_with_blanks")
    assert col["n_unique"] == 2
    assert col["missing_frac"] == 0.3333
    assert col["group_sizes"]["n_groups"] == 2
    assert "" not in col["examples"] and col["examples"]["<NA>"] == 40


def test_derived_candidates_carry_group_health():
    p = obsprofile.profile_obs(make_adata())
    barcode = next(d for d in p["derived"]
                   if d["label"] == "barcode:prefix:-")
    assert barcode["group_sizes"]["n_groups"] == 3
    assert 0 <= barcode["entropy"] <= 1
