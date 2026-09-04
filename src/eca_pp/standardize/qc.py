"""Matrix checks for the gates + F5's authoritative per-cell QC (spec §5.4).

``count_n_genes_detected`` matches the v2 pipeline's gate metric exactly. The sampled
checks (``is_integer_matrix``, ``has_negative``) mirror stancounts' sampling approach
so the two never disagree on what "integer" means. ``apply_qc`` is F5: the four
conventional QC obs columns, computed on the final gene space's true counts —
same-named columns brought by the data are preserved under ``__original`` and
overwritten (contract I4: this step's values are authoritative).
"""

from __future__ import annotations

import numpy as np
import scipy.sparse as sp

QC_COLS = ("pct_counts_mt", "pct_counts_hb", "total_counts", "n_genes_by_counts")


def count_n_genes_detected(counts) -> int:
    """Number of genes with counts>0 in >=1 cell (name-independent; cheap gate)."""
    if sp.issparse(counts):
        # scipy's getnnz counts stored entries, including explicitly stored
        # zeros.  Remove those representation-only entries before applying a
        # biological "counts > 0" gate.
        C = counts.tocsr(copy=True)
        C.eliminate_zeros()
        return int((np.asarray(C.getnnz(axis=0)).ravel() > 0).sum())
    C = np.asarray(counts)
    return int(((C != 0).sum(axis=0) > 0).sum())


def _sampled_nonzeros(M, n_sample: int = 200, seed: int = 0) -> np.ndarray:
    """Nonzero values from up to ``n_sample`` sampled rows, as float64."""
    n = M.shape[0]
    rng = np.random.RandomState(seed)
    idx = np.arange(n) if n <= n_sample else rng.choice(n, n_sample, replace=False)
    if sp.issparse(M):
        csr = M.tocsr()
        parts = [csr.data[csr.indptr[i]:csr.indptr[i + 1]] for i in idx]
        data = np.concatenate(parts) if parts else np.array([], dtype=float)
        data = data[data != 0]
    else:
        sub = np.asarray(M[idx])
        data = sub[sub != 0].ravel()
    return data.astype(np.float64)


def is_integer_matrix(M, *, n_sample: int = 200, seed: int = 0) -> bool:
    """True if sampled nonzero values are finite, non-negative, near-integer."""
    if M is None:
        return False
    data = _sampled_nonzeros(M, n_sample=n_sample, seed=seed)
    if data.size == 0:
        return False
    if not np.all(np.isfinite(data)) or np.any(data < 0):
        return False
    return bool(np.allclose(data, np.round(data), rtol=0, atol=1e-6))


def has_negative(M, *, n_sample: int = 200, seed: int = 0) -> bool:
    """True if sampled values contain negatives (scaled/z-scored data: the zero
    structure is not trustworthy for the provisional genes gate)."""
    if M is None:
        return False
    data = _sampled_nonzeros(M, n_sample=n_sample, seed=seed)
    return bool(data.size and np.any(data < 0))


def _masked_sum_per_cell(C, mask, is_sparse: bool, n_cells: int) -> np.ndarray:
    if not mask.any():
        return np.zeros(n_cells)
    if is_sparse:
        return np.asarray(C[:, mask].sum(axis=1)).ravel()
    return C[:, mask].sum(axis=1)


def apply_qc(adata, species: str) -> dict:
    """F5: write the four conventional QC obs columns from ``layers["counts"]``.

    mt/hb numerators come from stangene's exact per-species gene sets looked up
    on the harmonized ``var_names`` — no guessing. Pre-existing same-named obs
    columns are renamed to ``<name>__original`` before the authoritative values
    are written. Returns the result.json ``qc`` metrics block.
    """
    import stangene

    counts = adata.layers["counts"]
    is_sparse = sp.issparse(counts)
    C = counts.tocsr(copy=True) if is_sparse else np.asarray(counts)
    if is_sparse:
        # Keep sparse and dense QC semantics identical when an input CSR/CSC
        # happens to contain explicitly stored zeros.
        C.eliminate_zeros()
    n_cells = int(C.shape[0])
    names = list(adata.var_names)

    mt_mask = stangene.mito_mask(names, species)
    hb_mask = stangene.hb_mask(names, species)

    if is_sparse:
        counts_per_cell = np.asarray(C.sum(axis=1)).ravel().astype(np.float64)
        genes_per_cell = C.getnnz(axis=1).astype(np.int64)
    else:
        counts_per_cell = C.sum(axis=1).astype(np.float64)
        genes_per_cell = (C != 0).sum(axis=1).astype(np.int64)
    mt_per_cell = _masked_sum_per_cell(C, mt_mask, is_sparse, n_cells)
    hb_per_cell = _masked_sum_per_cell(C, hb_mask, is_sparse, n_cells)
    with np.errstate(divide="ignore", invalid="ignore"):
        pct_mt = np.where(counts_per_cell > 0, 100.0 * mt_per_cell / counts_per_cell, 0.0)
        pct_hb = np.where(counts_per_cell > 0, 100.0 * hb_per_cell / counts_per_cell, 0.0)

    overwritten = []
    for col, vals in zip(QC_COLS, (pct_mt, pct_hb, counts_per_cell, genes_per_cell)):
        if col in adata.obs.columns:
            adata.obs[f"{col}__original"] = adata.obs[col]
            overwritten.append(col)
        adata.obs[col] = vals

    med = (lambda a: float(np.median(a)) if n_cells else 0.0)
    return {
        "n_mt_genes": int(mt_mask.sum()),
        "n_hb_genes": int(hb_mask.sum()),
        "overwritten_obs_cols": overwritten,
        "median_pct_counts_mt": med(pct_mt),
        "median_pct_counts_hb": med(pct_hb),
        "median_counts_per_cell": med(counts_per_cell),
        "median_genes_per_cell": med(genes_per_cell),
    }
