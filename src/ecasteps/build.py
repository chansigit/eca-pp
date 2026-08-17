"""F7 — build the standard form and write it atomically (spec §4 ⑩, contract I1/I2).

``attach_counts`` runs right after F2 so the counts layer is subset together with
the rest of the AnnData by F4's gene dropping; ``finalize`` + ``write_h5ad`` run
last. The lognorm X is computed here directly (scipy), so ecasteps does not
depend on scanpy.
"""

from __future__ import annotations

import os

import numpy as np
import scipy.sparse as sp

from ecasteps.atomic_io import atomic_write

H5AD_FILENAME = "standardized.h5ad"
TARGET_SUM = 1e4


def lognorm(counts, target_sum: float = TARGET_SUM):
    """``log1p(normalize_total(counts, target_sum))`` as float32; sparse in →
    CSR out, dense in → ndarray out. All-zero cells stay all-zero."""
    if sp.issparse(counts):
        X = counts.tocsr().astype(np.float64)
        rs = np.asarray(X.sum(axis=1)).ravel()
        scale = np.divide(target_sum, rs, out=np.zeros_like(rs), where=rs > 0)
        X = sp.diags(scale) @ X
        X.data = np.log1p(X.data)
        return X.astype(np.float32).tocsr()
    C = np.asarray(counts, dtype=np.float64)
    rs = C.sum(axis=1, keepdims=True)
    scale = np.divide(target_sum, rs, out=np.zeros_like(rs), where=rs > 0)
    return np.log1p(C * scale).astype(np.float32)


def attach_counts(adata, counts, source: str) -> None:
    """``layers["counts"] = counts``; a differently-named source layer is renamed
    (i.e. deleted) so counts is never stored twice. Other layers are untouched."""
    adata.layers["counts"] = counts
    if source.startswith("layer:"):
        name = source.split(":", 1)[1]
        if name != "counts" and name in adata.layers:
            del adata.layers[name]


def finalize(adata, provenance: dict) -> bool:
    """X ← lognorm(counts); stash provenance in uns; drop ``.raw`` (its full
    gene space would be misaligned after F4's dropping — counts live in layers).
    Returns whether a raw was dropped."""
    adata.X = lognorm(adata.layers["counts"])
    raw_dropped = adata.raw is not None
    if raw_dropped:
        adata.raw = None
    adata.uns["ecasteps_standardize"] = provenance
    return raw_dropped


def write_h5ad(adata, outdir: str) -> str:
    """Atomically write ``outdir/standardized.h5ad``; returns the path."""
    path = os.path.join(outdir, H5AD_FILENAME)
    with atomic_write(path) as tmp:
        adata.write_h5ad(tmp)
    return path
