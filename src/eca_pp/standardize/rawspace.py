"""F1b — prefer ``.raw`` when it carries more genes than X (spec §5.0).

Scanpy workflows routinely leave the full gene space only in ``adata.raw``
(log-normalized) and reduce ``X`` / ``layers`` to a few thousand highly
variable genes. Standardizing the HVG subset would throw away most of the
transcriptome and usually fails the gene gate. When ``raw`` is a strict
superset of ``var_names`` the AnnData is therefore rebuilt on raw's gene
space before counts location: ``X = raw.X`` (counts are then found or
recovered from it by F2 exactly as for any other X), obs/obsm/uns are kept,
and the HVG-space layers are dropped. raw is trusted as the fuller record of
the same experiment; no cross-check against the HVG-space counts is made
(project decision, 2026-09-04). ``--no-raw-expand`` restores the old
behaviour.
"""

from __future__ import annotations

from dataclasses import dataclass, field

import anndata as ad
import pandas as pd
import scipy.sparse as sp

from eca_pp.core.layers import layer_names


@dataclass
class Expansion:
    applied: bool = False
    n_vars_x: int | None = None
    n_vars_raw: int | None = None
    reason: str = ""
    dropped_layers: list = field(default_factory=list)

    def as_dict(self) -> dict:
        return {
            "applied": self.applied,
            "n_vars_x": self.n_vars_x,
            "n_vars_raw": self.n_vars_raw,
            "reason": self.reason,
            "dropped_layers": self.dropped_layers,
        }


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
