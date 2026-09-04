"""Reused output directories must not expose earlier successful results."""

import pytest

from eca_pp.core.run_outputs import archive_outputs


def test_archives_only_step_outputs(tmp_path):
    for name in ("result.json", "standardized.h5ad", "notes.txt"):
        (tmp_path / name).write_text(name)
    archive_outputs(str(tmp_path), "standardize", str(tmp_path / "input.h5ad"))
    assert not (tmp_path / "standardized.h5ad").exists()
    assert not (tmp_path / "result.json").exists()
    assert (tmp_path / "notes.txt").read_text() == "notes.txt"
    saved = list((tmp_path / ".history").glob("*/standardized.h5ad"))
    assert len(saved) == 1
    assert saved[0].read_text() == "standardized.h5ad"


def test_input_overlap_is_checked_before_archiving(tmp_path):
    source = tmp_path / "standardized.h5ad"
    source.write_text("input")
    (tmp_path / "result.json").write_text("previous")
    with pytest.raises(ValueError, match="input overlaps"):
        archive_outputs(str(tmp_path), "standardize", str(source))
    assert source.read_text() == "input"
    assert (tmp_path / "result.json").read_text() == "previous"
