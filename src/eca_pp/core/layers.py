"""Layer enumeration that is safe across anndata versions.

GOTCHA (anndata >= 0.13): ``adata.layers`` exposes ``X`` under the key
``None``. Consequences that have surprised several people/agents already:

* ``len(adata.layers)`` is **1, not 0**, for a plain AnnData with only X;
* ``for name in adata.layers`` yields ``None`` first;
* ``list(adata.layers.keys())`` contains ``None``.

Any "does this object have other layers?" test therefore MUST go through
:func:`layer_names`, never ``len(adata.layers)`` / ``bool(adata.layers)``.
Verified on anndata 0.13.2 (the dl2025 venv on Sherlock); anndata 0.10/0.11
do not show the ``None`` key, so code written against them silently breaks
here (this is how the provisional genes gate in standardize stopped firing).
"""

from __future__ import annotations


def layer_names(adata) -> list[str]:
    """Names of the real layers of ``adata``, excluding the ``None`` alias of X."""
    return [name for name in adata.layers if name is not None]


def has_layers(adata) -> bool:
    """True when ``adata`` carries at least one real layer (X does not count)."""
    return bool(layer_names(adata))
