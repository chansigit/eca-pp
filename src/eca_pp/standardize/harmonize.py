"""F4 — gene-name harmonization: rename to canonical symbols, drop unmappable
features by default (spec §5.3).

stangene's five-tier cascade annotates every feature; this module applies the
policy: canonical symbols become ``var_names`` (make_unique on collisions, no
automatic merging), the original names and the full mapping provenance go into
``var``, and features with ``mapping_status`` in ``DROP_STATUSES`` are removed
from the matrix and every layer — unless ``keep_unmapped`` is set, in which
case they stay under their original names. An unmappable fraction above
``DROP_FRAC_REVIEW`` flags ``needs_review`` (usually a wrong species or a
non-gene feature table).
"""

from __future__ import annotations

import numpy as np
import pandas as pd

import stangene

MAPPING_COLS = [
    "gene_id_harmonized", "gene_symbol_harmonized", "mapping_status",
    "mapping_confidence", "mapping_source", "mapping_notes",
]
DROP_STATUSES = ("unmapped", "ambiguous", "non_gene_feature")
DROP_FRAC_REVIEW = 0.30


def _canonical_name(symbol, status, original: str) -> str:
    """Canonical symbol, or the original name when unmappable/blank."""
    if status in DROP_STATUSES:
        return original
    if symbol is None or (isinstance(symbol, float) and pd.isna(symbol)):
        return original
    s = str(symbol).strip()
    return s if s else original


def harmonize_genes(adata, species: str, *, keep_unmapped: bool = False):
    """Harmonize ``adata``'s features for ``species``.

    Returns ``(adata_out, stats, review_reasons)``. ``adata_out`` is a new
    object when features were dropped (X and all layers subset together), else
    the mutated input. ``stats`` matches the result.json ``harmonization``
    block; ``review_reasons`` is non-empty when the driver should look.
    """
    harmon = stangene.harmonize_anndata(adata, species)
    mt = harmon.mapping_table
    if len(mt) != adata.n_vars:
        raise ValueError(
            f"harmonization rows ({len(mt)}) != adata.n_vars ({adata.n_vars})")

    # Provenance first (row-aligned), so it travels through the subset below.
    adata.var["original_feature_name"] = list(adata.var_names)
    for col in MAPPING_COLS:
        if col in mt.columns:
            vals = pd.Series(list(mt[col]), index=adata.var.index)
            if vals.dtype == object:
                # stangene leaves None/NaN in text columns (e.g. mapping_notes);
                # h5py cannot write a mixed None/str column, so blank them.
                vals = vals.where(vals.notna(), "").map(str)
            adata.var[col] = vals

    status = np.array(
        [s if isinstance(s, str) and s else "unmapped"
         for s in mt["mapping_status"]], dtype=object)
    n = len(status)
    unmappable = {s: int((status == s).sum()) for s in DROP_STATUSES}
    n_unmappable = sum(unmappable.values())
    frac = (n_unmappable / n) if n else 0.0

    if keep_unmapped or n_unmappable == 0:
        out = adata
        dropped = {s: 0 for s in DROP_STATUSES}
        n_dropped = 0
    else:
        keep = ~np.isin(status, DROP_STATUSES)
        out = adata[:, keep].copy()
        dropped, n_dropped = unmappable, n_unmappable

    # Keys carry the unit (genes, never cells) so result.json is unambiguous.
    stats = {
        "genes_kept": int(out.n_vars),
        "genes_dropped": dropped,
        "genes_dropped_frac": round((n_dropped / n) if n else 0.0, 4),
        "genes_unmappable": unmappable,
        "keep_unmapped": bool(keep_unmapped),
    }
    reasons = []
    if frac > DROP_FRAC_REVIEW:
        fate = "kept by --keep-unmapped" if keep_unmapped else "dropped"
        reasons.append(
            f"{frac:.0%} of features are unmappable for species {species!r} "
            f"({fate}) — check the species call and the feature table")

    new_names = [
        _canonical_name(sym, st, orig)
        for sym, st, orig in zip(out.var["gene_symbol_harmonized"],
                                 out.var["mapping_status"],
                                 out.var["original_feature_name"])
    ]
    out.var_names = pd.Index(new_names, dtype=object)
    out.var_names_make_unique()
    return out, stats, reasons
