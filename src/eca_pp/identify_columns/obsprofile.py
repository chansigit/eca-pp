"""obs 画像 — the three-layer deterministic evidence base (identify-columns
spec §3): per-column stats with sampled values, group-size health for grouping
columns, and the nesting/equivalence graph among them, plus derived-candidate
enumeration (barcode prefix/suffix, two-column composites). Pure pandas/numpy;
agent and human reviewers read the same JSON.
"""

from __future__ import annotations

import itertools
import re

import numpy as np
import pandas as pd

TINY_GROUP_CELLS = 25     # groups below this are "tiny" (spec §6)
MAX_GROUPING_CARD = 1000  # above this a column is not a plausible grouping
MAX_EXAMPLES = 20
BARCODE_DELIMS = ("-", "_", ".", ":")
MAX_COMPOSITE_BASE = 50   # composite parts must each have <= this many groups
MAX_COMPOSITES = 10


def _entropy_norm(sizes: np.ndarray) -> float:
    p = sizes / sizes.sum()
    if len(p) < 2:
        return 0.0
    return float(-(p * np.log(p)).sum() / np.log(len(p)))


def _dtype_of(s: pd.Series) -> str:
    if isinstance(s.dtype, pd.CategoricalDtype):
        return "categorical"
    if pd.api.types.is_bool_dtype(s):
        return "bool"
    if pd.api.types.is_float_dtype(s):
        return "float"
    if pd.api.types.is_integer_dtype(s):
        return "int"
    return "string"


def _group_sizes_block(sizes: pd.Series, n_obs: int) -> dict:
    tiny = sizes[sizes < TINY_GROUP_CELLS]
    return {
        "n_groups": int(len(sizes)),
        "min": int(sizes.min()),
        "median": float(sizes.median()),
        "max": int(sizes.max()),
        "n_tiny": int(len(tiny)),
        "tiny_group_frac": round(len(tiny) / len(sizes), 4),
        "tiny_cell_frac": round(float(tiny.sum()) / n_obs, 4),
    }


def is_grouping_candidate(entry: dict) -> bool:
    """A column that could partition cells into usable groups."""
    return (entry["dtype"] in ("categorical", "string", "int", "bool")
            and 2 <= entry["n_unique"] <= MAX_GROUPING_CARD
            and not entry["is_per_cell_unique"])


def profile_obs(adata) -> dict:
    n = int(adata.n_obs)
    columns, grouping = [], {}
    for name in adata.obs.columns:
        s = adata.obs[name]
        vc = s.astype("string").fillna("<NA>").value_counts()
        entry = {
            "column": str(name),
            "dtype": _dtype_of(s),
            "n_unique": int(s.nunique(dropna=True)),
            "missing_frac": round(float(s.isna().mean()), 4),
            "is_constant": int(s.nunique(dropna=True)) <= 1,
            "is_per_cell_unique": int(s.nunique(dropna=True)) >= n,
            "examples": {str(k): int(v) for k, v in vc.head(MAX_EXAMPLES).items()},
        }
        if is_grouping_candidate(entry):
            entry["group_sizes"] = _group_sizes_block(vc, n)
            entry["entropy"] = round(_entropy_norm(vc.to_numpy(dtype=float)), 4)
            grouping[entry["column"]] = s.astype("string").fillna("<NA>")
        columns.append(entry)

    return {
        "n_obs": n,
        "columns": columns,
        "relations": _relations(grouping),
        "derived": (barcode_candidates(list(adata.obs_names))
                    + composite_candidates(grouping)),
    }


def _refines(fine: pd.Series, coarse: pd.Series) -> bool:
    """Every fine-group maps into exactly one coarse-group."""
    return bool((pd.DataFrame({"f": fine, "c": coarse})
                 .groupby("f", observed=True)["c"].nunique() <= 1).all())


def _relations(grouping: dict) -> list[dict]:
    out = []
    for a, b in itertools.combinations(sorted(grouping), 2):
        sa, sb = grouping[a], grouping[b]
        a_ref_b, b_ref_a = _refines(sa, sb), _refines(sb, sa)
        if a_ref_b and b_ref_a:
            out.append({"finer": a, "coarser": b, "kind": "equivalent"})
        elif a_ref_b and sa.nunique() > sb.nunique():
            out.append({"finer": a, "coarser": b, "kind": "nested"})
        elif b_ref_a and sb.nunique() > sa.nunique():
            out.append({"finer": b, "coarser": a, "kind": "nested"})
    return out


def barcode_candidates(obs_names: list) -> list[dict]:
    s = pd.Series([str(x) for x in obs_names])
    n = len(s)
    out = []
    for d in BARCODE_DELIMS:
        if not s.str.contains(re.escape(d), regex=True).all():
            continue
        for pos in ("prefix", "suffix"):
            vals = _split_barcode(s, pos, d)
            nun = int(vals.nunique())
            if 2 <= nun < n and nun <= MAX_GROUPING_CARD:
                out.append({"label": f"barcode:{pos}:{d}", "kind": "barcode",
                            "n_groups": nun})
    return out


def composite_candidates(grouping: dict) -> list[dict]:
    small = {k: v for k, v in grouping.items() if v.nunique() <= MAX_COMPOSITE_BASE}
    out = []
    for a, b in itertools.combinations(sorted(small), 2):
        combo = small[a].str.cat(small[b], sep="|")
        nun = int(combo.nunique())
        if (nun <= MAX_GROUPING_CARD
                and nun > max(small[a].nunique(), small[b].nunique())):
            out.append({"label": f"composite:{a}+{b}", "kind": "composite",
                        "n_groups": nun})
        if len(out) >= MAX_COMPOSITES:
            break
    return out


def _split_barcode(s: pd.Series, pos: str, delim: str) -> pd.Series:
    if pos == "prefix":
        return s.str.split(delim, n=1, regex=False).str[0]
    return s.str.rsplit(delim, n=1).str[-1]  # rsplit is always literal


def derive_values(adata, label: str) -> pd.Series:
    """Materialize a derived candidate's per-cell values (deterministic)."""
    kind, _, rest = label.partition(":")
    if kind == "barcode":
        pos, _, delim = rest.partition(":")
        s = pd.Series([str(x) for x in adata.obs_names], index=adata.obs_names)
        return _split_barcode(s, pos, delim)
    if kind == "composite":
        a, _, b = rest.partition("+")
        return (adata.obs[a].astype("string").fillna("<NA>")
                .str.cat(adata.obs[b].astype("string").fillna("<NA>"), sep="|"))
    raise ValueError(f"unknown derived-candidate label: {label!r}")
