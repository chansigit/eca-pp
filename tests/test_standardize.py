"""Acceptance tests, v0.1 set (spec §9): F1 → F3 → F2, result.json, exit codes.

Since v0.2 these run through the FULL pipeline (species → harmonize → QC →
write), so the synthetic datasets use real human reference genes (see dsets.py);
every counts-location assertion is unchanged.

Run on a compute node:  bash run.sh test tests -q
"""

from __future__ import annotations

import anndata as ad
import numpy as np
import scipy.sparse as sp
from dsets import RNG, G, N, lognorm, make_counts, run_cli, write_h5ad

from eca_pp.core.result import EXIT_BLOCKED, EXIT_OK, EXIT_REJECTED
from eca_pp.standardize import countsloc
from eca_pp.standardize.qc import count_n_genes_detected, is_integer_matrix

# ---------------------------------------------------------------- happy paths

def test_integer_counts_in_X(tmp_path):
    src = write_h5ad(tmp_path / "s.h5ad", sp.csr_matrix(make_counts()))
    code, res = run_cli(tmp_path, src)
    assert code == EXIT_OK and res["status"] == "ok"
    assert res["metrics"]["counts_source"] == "X"
    assert res["metrics"]["n_cells"] == N
    assert res["metrics"]["n_genes_detected"] >= 5000

def test_counts_in_whitelist_layer(tmp_path):
    c = make_counts()
    src = write_h5ad(tmp_path / "s.h5ad", lognorm(c), {"counts": sp.csr_matrix(c)})
    code, res = run_cli(tmp_path, src)
    assert code == EXIT_OK and res["status"] == "ok"
    assert res["metrics"]["counts_source"] == "layer:counts"

def test_odd_named_layer_adopted_by_consistency(tmp_path):
    c = make_counts()
    src = write_h5ad(tmp_path / "s.h5ad", lognorm(c), {"RNA_raw": sp.csr_matrix(c)})
    code, res = run_cli(tmp_path, src)
    assert code == EXIT_OK and res["status"] == "ok"
    assert res["metrics"]["counts_source"] == "layer:RNA_raw"
    assert res["metrics"]["counts_adopted_by"] == "consistency_check"
    assert res["metrics"]["counts_name_recognized"] is False
    (entry,) = [e for e in res["layers"] if e["name"] == "RNA_raw"]
    assert entry["consistent_with_X"] is True

def test_lognorm_only_recovers(tmp_path):
    src = write_h5ad(tmp_path / "s.h5ad", lognorm(make_counts()))
    code, res = run_cli(tmp_path, src)
    assert code == EXIT_OK and res["status"] == "ok"
    assert res["metrics"]["counts_source"].startswith("recovered")
    assert res["metrics"]["x_normalization"]["is_log1p"] is True

def test_scaled_X_with_counts_layer(tmp_path):
    c = make_counts()
    Xs = np.asarray(lognorm(c).todense())
    Xs = ((Xs - Xs.mean(0)) / (Xs.std(0) + 1e-6)).astype(np.float32)  # z-scored
    src = write_h5ad(tmp_path / "s.h5ad", Xs, {"counts": sp.csr_matrix(c)})
    code, res = run_cli(tmp_path, src)
    assert code == EXIT_OK and res["status"] == "ok"
    assert res["metrics"]["counts_source"] == "layer:counts"

def test_counts_layer_override(tmp_path):
    c = make_counts()
    src = write_h5ad(tmp_path / "s.h5ad", lognorm(c), {"mystery": sp.csr_matrix(c)})
    code, res = run_cli(tmp_path, src, "--counts-layer", "mystery")
    assert code == EXIT_OK and res["status"] == "ok"
    assert res["metrics"]["counts_source"] == "layer:mystery"
    assert res["metrics"]["counts_adopted_by"] == "override"

def test_no_gate_lets_tiny_sample_through(tmp_path):
    src = write_h5ad(tmp_path / "s.h5ad", sp.csr_matrix(make_counts(n=50)))
    code, res = run_cli(tmp_path, src, "--no-gate")
    assert code == EXIT_OK and res["status"] == "ok"
    assert res["metrics"]["n_cells"] == 50


def test_explicit_sparse_zeros_are_not_detected_genes_or_counts():
    stored_zeros = sp.csr_matrix(
        (np.array([4.0, 0.0, 0.0]), np.array([0, 1, 2]), np.array([0, 3])),
        shape=(1, 4),
    )
    assert count_n_genes_detected(stored_zeros) == 1

    all_zero = sp.csr_matrix(
        (np.array([0.0, 0.0]), np.array([0, 1]), np.array([0, 2])),
        shape=(1, 3),
    )
    assert is_integer_matrix(all_zero) is False


def test_counts_discovery_ignores_all_zero_sparse_layer():
    A = ad.AnnData(X=sp.csr_matrix([[1.0, 2.0, 0.0]]))
    A.layers["counts"] = sp.csr_matrix(
        (np.array([0.0, 0.0]), np.array([0, 1]), np.array([0, 2])),
        shape=(1, 3),
    )
    resolution = countsloc.resolve(A)
    assert resolution.source == "X"


# ---------------------------------------------------------- review & blocking

def test_recovered_with_inconsistent_layer_needs_review(tmp_path):
    c = make_counts()
    shuffled = c[:, RNG.permutation(G)]  # integer, but not X's counts
    src = write_h5ad(tmp_path / "s.h5ad", lognorm(c), {"weird": sp.csr_matrix(shuffled)})
    code, res = run_cli(tmp_path, src)
    assert code == EXIT_OK and res["status"] == "needs_review"
    assert res["metrics"]["counts_source"].startswith("recovered")
    assert any("weird" in r for r in res["reasons"])
    (entry,) = [e for e in res["layers"] if e["name"] == "weird"]
    assert entry["consistent_with_X"] is False

def test_ambiguous_candidates_block(tmp_path):
    c = make_counts()
    src = write_h5ad(tmp_path / "s.h5ad", lognorm(c),
                     {"a": sp.csr_matrix(c), "b": sp.csr_matrix(c.copy())})
    code, res = run_cli(tmp_path, src)
    assert code == EXIT_BLOCKED and res["status"] == "needs_review"
    assert res["rejected_at"] is None
    assert any("--counts-layer" in r for r in res["reasons"])

def test_missing_designated_layer_blocks(tmp_path):
    src = write_h5ad(tmp_path / "s.h5ad", sp.csr_matrix(make_counts()))
    code, res = run_cli(tmp_path, src, "--counts-layer", "nope")
    assert code == EXIT_BLOCKED and res["status"] == "needs_review"


# ----------------------------------------------------------------- rejections

def test_too_few_cells_rejects_before_load(tmp_path, monkeypatch):
    src = write_h5ad(tmp_path / "s.h5ad", sp.csr_matrix(make_counts(n=50)))

    def boom(*a, **k):  # the pre-gate must fire from HDF5 metadata alone
        raise AssertionError("read_h5ad must not be called for a pre-gate reject")

    monkeypatch.setattr(ad, "read_h5ad", boom)
    code, res = run_cli(tmp_path, src)
    assert code == EXIT_REJECTED and res["status"] == "rejected"
    assert res["rejected_at"] == "pre_gate"
    assert res["metrics"]["n_cells"] == 50

def test_too_few_genes_rejects(tmp_path):
    c = np.zeros((N, G), dtype=np.float32)
    c[:, :100] = RNG.poisson(1.0, size=(N, 100))
    src = write_h5ad(tmp_path / "s.h5ad", sp.csr_matrix(c))
    code, res = run_cli(tmp_path, src)
    assert code == EXIT_REJECTED and res["status"] == "rejected"
    # integer X with no other counts source: the provisional gate must fire
    # before counts discovery (it silently stopped doing so on anndata 0.13)
    assert res["rejected_at"] == "pre_gate"
    assert "f2_counts" not in res["metrics"]["timings"]


def test_misnamed_counts_layer_is_ignored_and_dropped(tmp_path):
    """A float layer called "counts" is treated as absent: counts come from
    elsewhere, the layer is dropped from the output, and result.json says so."""
    c = make_counts()
    src = write_h5ad(tmp_path / "s.h5ad", lognorm(c), {"counts": lognorm(c)})
    code, res = run_cli(tmp_path, src)
    assert code == EXIT_OK and res["status"] == "needs_review"
    assert res["metrics"]["counts_source"].startswith("recovered")
    assert res["metrics"]["ignored_counts_layers"] == ["counts"]
    assert any("named like counts" in r for r in res["reasons"])
    out = ad.read_h5ad(res["output"])
    L = out.layers["counts"]
    assert is_integer_matrix(L)
    assert out.uns["eca_pp_standardize"]["ignored_counts_layers"] == "counts"


def test_review_notes_survive_a_later_block(tmp_path, monkeypatch):
    """Counts-stage doubts must reach result.json even when the species
    stage blocks afterwards (they used to live in a local list)."""
    import eca_pp.standardize.species as species_ladder

    c = make_counts()
    shuffled = c[:, RNG.permutation(G)]
    src = write_h5ad(tmp_path / "s.h5ad", lognorm(c), {"weird": sp.csr_matrix(shuffled)})
    monkeypatch.setattr(species_ladder, "_t1_infer", lambda adata: {
        "species": None, "confidence": 0.0, "evidence": {}})
    code, res = run_cli(tmp_path, src)
    assert code == EXIT_BLOCKED and res["status"] == "needs_review"
    assert "--species" in res["reasons"][0]          # the stop reason leads
    assert any("weird" in r for r in res["reasons"])  # the earlier note follows


def test_unknown_species_code_blocks_before_loading(tmp_path, monkeypatch):
    src = write_h5ad(tmp_path / "s.h5ad", sp.csr_matrix(make_counts()))

    def boom(*a, **k):
        raise AssertionError("a bad --species must be rejected before read_h5ad")

    monkeypatch.setattr(ad, "read_h5ad", boom)
    code, res = run_cli(tmp_path, src, "--species", "martian")
    assert code == EXIT_BLOCKED and res["status"] == "needs_review"
    assert any("martian" in r for r in res["reasons"])

def test_not_h5ad_rejects(tmp_path):
    junk = tmp_path / "junk.h5ad"
    junk.write_text("this is not hdf5")
    code, res = run_cli(tmp_path, junk)
    assert code == EXIT_REJECTED and res["status"] == "rejected"
    assert res["rejected_at"] == "input"
    assert res["reasons"]
