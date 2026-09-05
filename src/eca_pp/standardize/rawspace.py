"""F1b — prefer ``.raw`` when it carries more genes than X (spec §5.0).

Scanpy workflows routinely leave the full gene space only in ``adata.raw``
(log-normalized) and reduce ``X`` / ``layers`` to a few thousand highly
variable genes. Standardizing the HVG subset would throw away most of the
transcriptome and usually fails the gene gate. When ``raw`` is a strict
superset of ``var_names`` the AnnData is therefore rebuilt on raw's gene
space before counts location: ``X = raw.X`` (counts are then found or
recovered from it by F2 exactly as for any other X), obs/obsm/uns are kept,
and the HVG-space layers are dropped — but an integer counts matrix found in
the HVG space is retained as a *reference* and compared value-by-value with
the recovered counts on the shared genes, so a silent mismatch (raw made
from different counts, wrong normalization target) is flagged for review.
"""

from __future__ import annotations

from dataclasses import dataclass, field

import anndata as ad
import numpy as np
import pandas as pd
import scipy.sparse as sp

from eca_pp.core.layers import layer_names

MATCH_FRAC_REVIEW = (
    0.99  # below this the recovered counts disagree with the HVG reference
)


@dataclass
class Expansion:
    applied: bool = False
    n_vars_x: int | None = None
    n_vars_raw: int | None = None
    reason: str = ""
    dropped_layers: list = field(default_factory=list)
    reference_source: str | None = None
    reference: object = None  # (matrix, var_names) in the HVG space
    check: dict | None = None

    def as_dict(self) -> dict:
        out = {
            "applied": self.applied,
            "n_vars_x": self.n_vars_x,
            "n_vars_raw": self.n_vars_raw,
            "reason": self.reason,
            "dropped_layers": self.dropped_layers,
            "reference_source": self.reference_source,
        }
        if self.check is not None:
            out["counts_check"] = self.check
        return out


def _reference_counts(adata):
    """An integer counts matrix in the current (subset) gene space, if any."""
    import stancounts

    try:
        got = stancounts.get_counts(adata, allow_recovery=False)
    except stancounts.CountsUnavailable:
        return None
    return got["counts"], got["source"]


def maybe_expand(adata, *, enabled: bool = True) -> tuple[ad.AnnData, Expansion]:
    """Return ``(adata_on_raw_gene_space, expansion)``; the input when not applied."""
    exp = Expansion()
    raw = adata.raw
    if raw is None:
        exp.reason = "no .raw"
        return adata, exp
    exp.n_vars_x, exp.n_vars_raw = int(adata.n_vars), int(raw.n_vars)
    if not enabled:
        exp.reason = "disabled by --no-raw-expand"
        return adata, exp
    if raw.n_vars <= adata.n_vars:
        exp.reason = "raw has no more genes than X"
        return adata, exp
    if raw.X is None:
        exp.reason = "raw has no matrix"
        return adata, exp
    raw_names = pd.Index(raw.var_names)
    missing = adata.var_names.difference(raw_names)
    if len(missing):
        exp.reason = (
            f"raw lacks {len(missing)} of X's {adata.n_vars} genes "
            f"(first: {missing[0]!r}); gene spaces do not align"
        )
        return adata, exp

    ref = _reference_counts(adata)
    if ref is not None:
        exp.reference_source = ref[1]
        exp.reference = (ref[0], pd.Index(adata.var_names))
    X = raw.X
    if sp.issparse(X) and not sp.isspmatrix_csr(X):
        X = X.tocsr()
    full = ad.AnnData(
        X=X,
        obs=adata.obs.copy(),
        var=raw.var.copy(),
        obsm=dict(adata.obsm) if adata.obsm is not None else None,
        uns=dict(adata.uns),
    )
    full.var_names = raw_names
    exp.applied = True
    exp.dropped_layers = layer_names(adata)
    exp.reason = (
        f"raw carries {raw.n_vars} genes vs {adata.n_vars} in X; "
        f"rebuilt on raw's gene space"
    )
    return full, exp


def verify_against_reference(
    counts, var_names, exp: Expansion, *, n_rows: int = 200, seed: int = 0
) -> dict | None:
    """Compare ``counts`` (full gene space) with the retained HVG-space
    reference on their shared genes; fills ``exp.check`` and returns it."""
    if exp.reference is None:
        return None
    ref, ref_names = exp.reference
    names = pd.Index(var_names)
    shared = ref_names.intersection(names)
    if not len(shared):
        exp.check = {
            "reference": exp.reference_source,
            "n_shared_genes": 0,
            "match_frac": None,
        }
        return exp.check
    full_idx = names.get_indexer(shared)
    ref_idx = ref_names.get_indexer(shared)
    n = counts.shape[0]
    rng = np.random.default_rng(seed)
    rows = np.arange(n) if n <= n_rows else rng.choice(n, n_rows, replace=False)
    C = counts.tocsr() if sp.issparse(counts) else np.asarray(counts)
    R = ref.tocsr() if sp.issparse(ref) else np.asarray(ref)
    matched = compared = 0
    for i in rows:
        a = C[i]
        b = R[i]
        a = np.asarray(a.todense()).ravel() if sp.issparse(a) else np.asarray(a).ravel()
        b = np.asarray(b.todense()).ravel() if sp.issparse(b) else np.asarray(b).ravel()
        a, b = np.round(a[full_idx]), np.round(b[ref_idx])
        mask = (a != 0) | (b != 0)
        compared += int(mask.sum())
        matched += int((a[mask] == b[mask]).sum())
    frac = (matched / compared) if compared else 1.0
    exp.check = {
        "reference": exp.reference_source,
        "n_shared_genes": int(len(shared)),
        "n_cells_sampled": int(len(rows)),
        "n_values_compared": compared,
        "match_frac": round(float(frac), 4),
    }
    return exp.check
