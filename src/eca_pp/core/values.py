"""Shared normalization for categorical metadata values."""

from __future__ import annotations

import pandas as pd


# Conservative spellings that unambiguously mean "not recorded" in public
# single-cell metadata.  Deliberately exclude labels such as "unknown" and
# "unassigned": those can be meaningful author annotation categories.
MISSING_STRINGS = frozenset({"", "na", "n/a", "nan", "none", "null", "<na>",
                             "missing"})


def normalize_missing(values: pd.Series) -> pd.Series:
    """Return values with blank/common missing strings represented as ``NA``.

    Numeric and boolean series retain their dtype.  String-like/categorical
    series become pandas' nullable string dtype so callers can consistently
    exclude missing values from group counts and model fitting.
    """
    if not (isinstance(values.dtype, pd.CategoricalDtype)
            or pd.api.types.is_object_dtype(values)
            or pd.api.types.is_string_dtype(values)):
        return values.copy()
    out = values.astype("string")
    normalized = out.str.strip().str.lower()
    return out.mask(normalized.isin(MISSING_STRINGS), pd.NA)
