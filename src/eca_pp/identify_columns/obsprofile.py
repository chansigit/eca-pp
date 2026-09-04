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

from eca_pp.core.values import normalize_missing

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
        original = adata.obs[name]
        s = normalize_missing(original)
        vc = s.astype("string").value_counts(dropna=True)
        missing_n = int(s.isna().sum())
        examples = {str(k): int(v) for k, v in vc.head(MAX_EXAMPLES).items()}
        if missing_n and len(examples) < MAX_EXAMPLES:
            examples["<NA>"] = missing_n
        entry = {
            "column": str(name),
            "dtype": _dtype_of(original),
            "n_unique": int(s.nunique(dropna=True)),
            "missing_frac": round(float(s.isna().mean()), 4),
            "is_constant": int(s.nunique(dropna=True)) <= 1,
            "is_per_cell_unique": int(s.nunique(dropna=True)) >= n,
            "examples": examples,
        }
        if is_grouping_candidate(entry):
            entry["group_sizes"] = _group_sizes_block(vc, n)
            entry["entropy"] = round(_entropy_norm(vc.to_numpy(dtype=float)), 4)
            grouping[entry["column"]] = s.astype("string")
        columns.append(entry)

    return {
        "n_obs": n,
        "columns": columns,
        "relations": _relations(grouping),
        "derived": (barcode_candidates(list(adata.obs_names), grouping)
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


def _equivalent_with(values: pd.Series, grouping: dict) -> list[str]:
    """Existing columns that encode exactly the same cell partition."""
    left = pd.factorize(values.astype("string").fillna("<NA>"), sort=False)[0]
    matches = []
    for name, other in grouping.items():
        if other.nunique(dropna=False) != values.nunique(dropna=False):
            continue
        right = pd.factorize(
            other.astype("string").fillna("<NA>"), sort=False
        )[0]
        if np.array_equal(left, right):
            matches.append(name)
    return sorted(matches)


def barcode_candidates(obs_names: list, grouping: dict | None = None) -> list[dict]:
    s = pd.Series([str(x) for x in obs_names])
    n = len(s)
    grouping = grouping or {}
    out = []
    for d in BARCODE_DELIMS:
        if not s.str.contains(re.escape(d), regex=True).all():
            continue
        for pos in ("prefix", "suffix"):
            vals = _split_barcode(s, pos, d)
            nun = int(vals.nunique())
            if 2 <= nun < n and nun <= MAX_GROUPING_CARD:
                sizes = vals.value_counts()
                out.append({"label": f"barcode:{pos}:{d}", "kind": "barcode",
                            "n_groups": nun,
                            "missing_frac": 0.0,
                            "equivalent_with": _equivalent_with(vals, grouping),
                            "group_sizes": _group_sizes_block(sizes, n),
                            "entropy": round(
                                _entropy_norm(sizes.to_numpy(dtype=float)), 4)})
    return out


def composite_candidates(grouping: dict) -> list[dict]:
    small = {k: v for k, v in grouping.items() if v.nunique() <= MAX_COMPOSITE_BASE}
    out = []
    for a, b in itertools.combinations(sorted(small), 2):
        combo = small[a].str.cat(small[b], sep="|")
        nun = int(combo.nunique())
        if (nun <= MAX_GROUPING_CARD
                and nun > max(small[a].nunique(), small[b].nunique())):
            sizes = combo.value_counts(dropna=True)
            out.append({"label": f"composite:{a}+{b}", "kind": "composite",
                        "n_groups": nun,
                        "missing_frac": round(float(combo.isna().mean()), 4),
                        "equivalent_with": _equivalent_with(combo, grouping),
                        "group_sizes": _group_sizes_block(sizes, len(combo)),
                        "entropy": round(
                            _entropy_norm(sizes.to_numpy(dtype=float)), 4)})
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
        return normalize_missing(adata.obs[a]).astype("string").str.cat(
            normalize_missing(adata.obs[b]).astype("string"), sep="|")
    raise ValueError(f"unknown derived-candidate label: {label!r}")
