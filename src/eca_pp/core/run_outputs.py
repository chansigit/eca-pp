"""Retain previous step outputs before reusing an output directory."""

import logging
import os
import re
import tempfile
from pathlib import Path

log = logging.getLogger(__name__)


def archive_outputs(outdir: str, step: str, src: str) -> None:
    """Move only known step outputs into a recoverable history directory.

    An input inside the previous outputs requires a separate destination.
    Validate every target before moving anything; never archive the input.
    """
    root = Path(outdir)
    if not root.is_dir():
        return
    names = {"result.json"}
    if step == "standardize":
        names.add("standardized.h5ad")
    elif step == "identify_columns":
        names.update({"batch.tsv", "candidates"})
    targets = [p for p in root.iterdir() if p.name in names or (
        step == "identify_columns" and re.fullmatch(r"trial_[0-9]+", p.name))]
    source = Path(src).resolve()
    for target in targets:
        resolved = target.resolve()
        if source == resolved or resolved in source.parents or (
            source.exists() and target.is_file() and os.path.samefile(source, target)
        ):
            raise ValueError("input overlaps previous outputs; use a different output directory")
    if not targets:
        return
    history = root / ".history"
    history.mkdir(exist_ok=True)
    destination = Path(tempfile.mkdtemp(prefix=f"{step}-", dir=history))
    for target in targets:
        target.rename(destination / target.name)
    log.info("previous %s outputs archived at %s", step, destination)
