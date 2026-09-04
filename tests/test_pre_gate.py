"""Alternate counts sources must survive the provisional X gate."""

import anndata as ad
import numpy as np
import pytest

from eca_pp.standardize.cli import _pre_gate_matrix


@pytest.mark.parametrize("alternate", ["layer", "raw"])
def test_filtered_x_cannot_reject_alternate_counts(alternate):
    adata = ad.AnnData(np.ones((3, 4)))
    if alternate == "raw":
        adata.raw = adata.copy()
    else:
        adata.layers["counts"] = adata.X.copy()
    adata.X[:] = 0
    assert _pre_gate_matrix(adata, None) == (None, False, False)


def test_explicit_counts_still_supports_early_gate():
    adata = ad.AnnData(np.zeros((3, 4)))
    adata.layers["counts"] = np.ones((3, 4))
    matrix, trusted, exact = _pre_gate_matrix(adata, "counts")
    assert trusted and exact
    assert matrix is adata.layers["counts"]
