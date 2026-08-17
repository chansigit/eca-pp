"""Synthetic datasets for the integration-probe / identify-columns tests.

Two cell types (disjoint marker genes) so cLISI is meaningful, plus an
optional multiplicative batch effect on a third gene block so the true batch
column separates in PCA space while a shuffled label does not.
"""

from __future__ import annotations

import anndata as ad
import numpy as np
import pandas as pd
import scipy.sparse as sp


def make_integration_h5ad(path, *, n_per_batch=300, n_genes=400, n_batches=2,
                          effect=4.0, seed=0, obs_extra=None, barcode_batch=False):
    rng = np.random.default_rng(seed)
    n = n_per_batch * n_batches
    counts = rng.poisson(1.0, size=(n, n_genes)).astype(np.float32)

    batch = np.repeat([f"b{i}" for i in range(n_batches)], n_per_batch)
    cell_type = np.tile(
        np.array(["typeA", "typeB"]).repeat(n_per_batch // 2), n_batches)[:n]
    counts[cell_type == "typeA", 0:40] += rng.poisson(6.0, size=(int((cell_type == "typeA").sum()), 40))
    counts[cell_type == "typeB", 40:80] += rng.poisson(6.0, size=(int((cell_type == "typeB").sum()), 80 - 40))
    if effect and effect != 1.0:
        for i in range(1, n_batches, 2):  # odd batches carry the effect
            rows = batch == f"b{i}"
            counts[np.ix_(rows, np.arange(100, 160))] *= effect

    C = counts.astype(np.float64)
    s = C.sum(axis=1, keepdims=True)
    s[s == 0] = 1.0
    X = sp.csr_matrix(np.log1p(C / s * 1e4).astype(np.float32))

    A = ad.AnnData(X=X)
    A.var_names = pd.Index([f"gene{i}" for i in range(n_genes)], dtype=object)
    A.layers["counts"] = sp.csr_matrix(counts)
    if barcode_batch:
        A.obs_names = pd.Index([f"{b}-CELL{i:05d}" for i, b in enumerate(batch)],
                               dtype=object)
    else:
        A.obs_names = pd.Index([f"CELL{i:05d}" for i in range(n)], dtype=object)
        A.obs["batch"] = batch
    A.obs["cell_type"] = cell_type
    for k, v in (obs_extra or {}).items():
        A.obs[k] = v
    A.write_h5ad(path)
    return path
