"""Column SPEC resolution — the `--batch-col` / `--cell-type-col` convention
shared by every consumer (identify-columns spec §5/§8): a SPEC is either an
obs column name or the path of a TSV value file (`cell_id<TAB>value`, header
row optional). Derived columns are exchanged as TSV files; input h5ads are
never modified.
"""

from __future__ import annotations

import os

import pandas as pd


class ColumnSpecError(ValueError):
    """SPEC cannot be resolved against this AnnData (a permanent input problem)."""


def resolve_spec(adata, spec: str) -> tuple[pd.Series, str]:
    """Resolve ``spec`` to per-cell values aligned to ``adata.obs_names``.

    Returns ``(values, kind)`` with kind ``"column"`` or ``"tsv"``. Raises
    :class:`ColumnSpecError` when the column is absent or the TSV does not
    cover every cell.
    """
    if spec in adata.obs.columns:
        return adata.obs[spec].astype("string").fillna("<NA>"), "column"
    if os.path.isfile(spec):
        mapping = read_values_tsv(spec)
        missing = [c for c in map(str, adata.obs_names) if c not in mapping.index]
        if missing:
            raise ColumnSpecError(
                f"value file {spec!r} lacks {len(missing)} of {adata.n_obs} "
                f"cells (first missing: {missing[0]!r})")
        vals = mapping.reindex([str(c) for c in adata.obs_names])
        vals.index = adata.obs_names
        return vals.astype("string"), "tsv"
    raise ColumnSpecError(
        f"{spec!r} is neither an obs column ({list(adata.obs.columns)[:10]}...) "
        f"nor an existing value file")


def read_values_tsv(path: str) -> pd.Series:
    df = pd.read_csv(path, sep="\t", header=None, dtype=str,
                     names=["cell_id", "value"])
    if len(df) and df.iloc[0, 0] == "cell_id":  # optional header row
        df = df.iloc[1:]
    if df["cell_id"].duplicated().any():
        raise ColumnSpecError(f"value file {path!r} has duplicated cell_ids")
    return df.set_index("cell_id")["value"]


def write_values_tsv(path: str, obs_names, values) -> str:
    from ecasteps.core.atomic_io import atomic_write

    df = pd.DataFrame({"cell_id": [str(x) for x in obs_names],
                       "value": [str(v) for v in values]})
    with atomic_write(path) as tmp:
        df.to_csv(tmp, sep="\t", index=False)
    return path
