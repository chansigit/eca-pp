"""Storage dtype must not hide discrete batch IDs or admit QC measurements."""

import anndata as ad
import numpy as np
import pandas as pd

from eca_pp.identify_columns.cli import build_candidates
from eca_pp.identify_columns.obsprofile import profile_obs


def test_float_group_ids_are_candidates_but_qc_is_not():
    adata = ad.AnnData(np.zeros((100, 1)))
    adata.obs["batch"] = [1.0] * 50 + [2.0] * 49 + [np.nan]
    adata.obs["total_counts"] = [100.0, 200.0] * 50
    adata.obs["continuous"] = np.linspace(0.1, 0.9, 100)
    candidates = build_candidates(profile_obs(adata))
    existing = {c["label"] for c in candidates["batch"] if c["kind"] == "existing"}
    assert "batch" in existing
    assert "total_counts" not in existing
    assert "continuous" not in existing
    assert all("total_counts" not in c["label"] for c in candidates["batch"])
    # Profiling must not coerce the author's stored values.
    assert pd.api.types.is_float_dtype(adata.obs["batch"])
