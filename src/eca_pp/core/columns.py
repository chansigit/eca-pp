"""Collision-safe preservation of ECAPP-managed metadata columns."""

from __future__ import annotations

import re

import pandas as pd


def preserve_column(frame: pd.DataFrame, name: str) -> str | None:
    """Preserve ``frame[name]`` in its ``__original`` backup family.

    Existing backups with identical values are reduced to the shortest column
    name. If one of them already matches the current column, it is reused.
    Otherwise the current values are copied to the first free numbered backup.
    Unrelated columns are never considered for deduplication.

    Returns the backup column name, or ``None`` when ``name`` is absent.
    """
    if name not in frame.columns:
        return None

    base = f"{name}__original"
    family_pattern = re.compile(rf"^{re.escape(base)}(?:_[1-9][0-9]*)?$")
    family = sorted(
        (col for col in frame.columns if family_pattern.fullmatch(str(col))),
        key=lambda col: (len(str(col)), str(col)),
    )

    # Keep only the shortest name for each exact backup value. Series.equals
    # deliberately requires matching indexes, dtypes, values, and missingness.
    kept: list[str] = []
    for column in family:
        if any(frame[column].equals(frame[other]) for other in kept):
            del frame[column]
        else:
            kept.append(column)

    for column in kept:
        if frame[name].equals(frame[column]):
            return column

    destination = base
    suffix = 2
    while destination in frame.columns:
        destination = f"{base}_{suffix}"
        suffix += 1
    frame[destination] = frame[name].copy()
    return destination
