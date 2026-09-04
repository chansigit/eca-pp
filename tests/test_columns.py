"""Tests for collision-safe metadata preservation."""

from types import SimpleNamespace

import anndata as ad
import numpy as np
import pandas as pd
import scipy.sparse as sp

from eca_pp.core.columns import preserve_column
from eca_pp.standardize.harmonize import harmonize_genes
from eca_pp.standardize.qc import apply_qc


def test_identical_backups_reduce_to_shortest_name():
    frame = pd.DataFrame({
        "metric": [1, 2],
        "metric__original": [1, 2],
        "metric__original_2": [1, 2],
    })

    destination = preserve_column(frame, "metric")

    assert destination == "metric__original"
    assert list(frame.columns) == ["metric", "metric__original"]


def test_distinct_values_use_next_free_backup():
    frame = pd.DataFrame({
        "metric": [3, 4],
        "metric__original": [1, 2],
    })

    destination = preserve_column(frame, "metric")

    assert destination == "metric__original_2"
    assert frame["metric__original"].tolist() == [1, 2]
    assert frame["metric__original_2"].tolist() == [3, 4]


def test_unrelated_identical_column_is_not_deduplicated():
    frame = pd.DataFrame({
        "donor": ["a", "b"],
        "sample_id": ["a", "b"],
    })

    preserve_column(frame, "donor")

    assert "sample_id" in frame.columns
    assert "donor__original" in frame.columns


def test_missing_source_is_a_noop():
    frame = pd.DataFrame({"other": [1, 2]})

    assert preserve_column(frame, "metric") is None
    assert list(frame.columns) == ["other"]


def test_qc_preserves_distinct_existing_backup():
    adata = ad.AnnData(X=sp.csr_matrix([[1, 2], [3, 0]]))
    adata.var_names = ["MT-ND1", "GAPDH"]
    adata.layers["counts"] = adata.X.copy()
    adata.obs["pct_counts_mt"] = [77.0, 77.0]
    adata.obs["pct_counts_mt__original"] = [66.0, 66.0]

    apply_qc(adata, "human")

    assert adata.obs["pct_counts_mt__original"].tolist() == [66.0, 66.0]
    assert adata.obs["pct_counts_mt__original_2"].tolist() == [77.0, 77.0]
    np.testing.assert_allclose(adata.obs["pct_counts_mt"], [100 / 3, 100])


def test_harmonization_preserves_distinct_provenance(monkeypatch):
    import eca_pp.standardize.harmonize as harmonize_module

    adata = ad.AnnData(X=sp.csr_matrix([[1, 2]]))
    adata.var_names = ["G1", "G2"]
    adata.var["original_feature_name"] = ["legacy-a", "legacy-b"]
    adata.var["original_feature_name__original"] = ["older-a", "older-b"]
    adata.var["mapping_status"] = ["legacy", "legacy"]
    adata.var["mapping_status__original"] = ["older", "older"]
    mapping = pd.DataFrame({
        "gene_symbol_harmonized": ["NEW1", "NEW2"],
        "mapping_status": ["approved_symbol", "approved_symbol"],
    })
    monkeypatch.setattr(
        harmonize_module.stangene,
        "harmonize_anndata",
        lambda *_args, **_kwargs: SimpleNamespace(mapping_table=mapping),
    )

    out, _stats, _reasons = harmonize_genes(adata, "human")

    assert out.var["original_feature_name__original"].tolist() == [
        "older-a", "older-b"]
    assert out.var["original_feature_name__original_2"].tolist() == [
        "legacy-a", "legacy-b"]
    assert out.var["mapping_status__original"].tolist() == ["older", "older"]
    assert out.var["mapping_status__original_2"].tolist() == ["legacy", "legacy"]
    assert out.var["mapping_status"].tolist() == [
        "approved_symbol", "approved_symbol"]
