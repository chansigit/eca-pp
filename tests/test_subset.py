"""subset tests (subset spec): the subset keeps the standard form (I1-I4), says why each
dropped cell went, re-checks the source's gates, and a product whose X holds raw counts
is rejected by the probe and by identify-columns before any model call."""

from __future__ import annotations

import json

import anndata as ad
import numpy as np
import pandas as pd
import scipy.sparse as sp
from dsets import lognorm, make_counts, run_cli, write_h5ad
from intdata import make_integration_h5ad

from eca_pp.core.result import EXIT_OK, EXIT_REJECTED
from eca_pp.core.run_outputs import archive_outputs
from eca_pp.subset.cli import main


def product(tmp_path, n=600, **obs):
    """A standardize product (its OUTDIR) of n synthetic cells with the given obs columns."""
    c = make_counts(n=n)
    src = write_h5ad(tmp_path / "s.h5ad", lognorm(c), {"counts": sp.csr_matrix(c)}, obs=obs)
    code, _ = run_cli(tmp_path, src)
    assert code == EXIT_OK
    return tmp_path / "out"


def subset(tmp_path, src, *extra, name="sub"):
    out = tmp_path / name
    code = main([str(src), "-o", str(out), *map(str, extra)])
    res = json.loads((out / "result.json").read_text())
    assert res["exit_code"] == code
    return code, res, out


def test_subset_keeps_the_standard_form_and_says_why_each_cell_went(tmp_path):
    method = np.array(["10X", "smartseq"] * 300)
    src = product(tmp_path, method=method)
    full = ad.read_h5ad(src / "standardized.h5ad")
    C = sp.csr_matrix(full.layers["counts"])
    genes, total = np.diff(C.indptr), np.asarray(C.sum(axis=1)).ravel()
    # thresholds equal to real cells' values: those cells must go (strictly greater keeps)
    g_cell, c_cell = np.argsort(genes)[60], np.argsort(total)[60]
    min_genes, min_counts = int(genes[g_cell]), float(total[c_cell])
    code, res, out = subset(tmp_path, src, "--obs", "method=10X", "--min-genes", min_genes,
                            "--min-counts", min_counts)
    assert code == EXIT_OK and res["step"] == "subset"
    keep = (method == "10X") & (genes > min_genes) & (total > min_counts)
    sub = ad.read_h5ad(out / "standardized.h5ad")
    assert list(sub.obs_names) == list(full.obs_names[keep])
    # I1: one integer counts layer, the source's counts of the kept cells
    assert {k for k in sub.layers.keys() if k is not None} == {"counts"}
    L = sub.layers["counts"].toarray()
    np.testing.assert_array_equal(L, C[keep].toarray())
    # I2: X = log1p(normalize_total(counts, 1e4)), float32, on the subset
    assert sub.X.dtype == np.float32
    np.testing.assert_allclose(sub.X[0].toarray().ravel(), np.log1p(L[0] / L[0].sum() * 1e4), rtol=1e-5, atol=1e-5)
    # I4: QC recomputed, the source's values not backed up as __original
    np.testing.assert_allclose(sub.obs["total_counts"], L.sum(axis=1))
    assert not any(c.endswith("__original") for c in sub.obs.columns)
    assert sub.uns["eca_pp_standardize"]["species"] == "human" and "eca_pp_subset" in sub.uns
    d = res["derived_from"]
    assert d["n_source_cells"] == 600 and d["n_dropped"] == int((~keep).sum()) and len(d["src_sha256"]) == 64
    assert res["metrics"]["n_cells"] == int(keep.sum()) and res["metrics"]["x_normalization"]["is_log1p"]
    dropped = pd.read_csv(out / "dropped_cells.csv.gz", dtype=str).set_index("cell_id").reason
    assert set(dropped.index) == set(full.obs_names[~keep])
    names = list(full.obs_names)
    assert f"n_genes<={min_genes}" in dropped[names[g_cell]].split("+")
    assert f"n_counts<={min_counts:g}" in dropped[names[c_cell]].split("+")
    smart = names[1]
    assert dropped[smart].split("+")[0] == "obs method!=10X"


def test_graphs_always_go_and_counts_only_drops_the_rest(tmp_path):
    src = product(tmp_path)
    path = src / "standardized.h5ad"
    A = ad.read_h5ad(path)
    A.obsm["X_umap"] = np.zeros((A.n_obs, 2))
    A.varm["loadings"] = np.zeros((A.n_vars, 2))
    A.obsp["connectivities"] = sp.identity(A.n_obs, format="csr")
    A.layers["scale_data"] = A.X.copy()
    A.write_h5ad(path)
    code, _, out = subset(tmp_path, src)
    sub = ad.read_h5ad(out / "standardized.h5ad")
    assert code == EXIT_OK and not len(sub.obsp) and "X_umap" in sub.obsm and "scale_data" in sub.layers
    code, _, out = subset(tmp_path, src, "--counts-only", name="lean")
    lean = ad.read_h5ad(out / "standardized.h5ad")
    assert code == EXIT_OK and not len(lean.obsm) and not len(lean.varm)
    assert {k for k in lean.layers.keys() if k is not None} == {"counts"}
    # read without the other matrices, the same cells, counts, X, QC and provenance
    assert (lean.layers["counts"] != sub.layers["counts"]).nnz == 0 and np.allclose(lean.X.toarray(), sub.X.toarray())
    pd.testing.assert_frame_equal(lean.obs, sub.obs)
    assert lean.uns["eca_pp_standardize"] == sub.uns["eca_pp_standardize"]


def test_a_subset_below_the_source_gates_is_rejected(tmp_path):
    src = product(tmp_path, n=300, method=np.array(["rare"] * 50 + ["10X"] * 250))
    code, res, out = subset(tmp_path, src, "--obs", "method=rare")
    assert code == EXIT_REJECTED and res["rejected_at"] == "final_gate"
    assert "n_cells 50 < min_cells 100" in res["reasons"]
    assert not (out / "standardized.h5ad").exists()


def test_a_source_that_is_not_ready_is_rejected(tmp_path):
    src = product(tmp_path, n=300)
    rj = src / "result.json"
    rj.write_text(json.dumps({**json.loads(rj.read_text()), "status": "rejected"}))
    code, res, _ = subset(tmp_path, src)
    assert code == EXIT_REJECTED and res["rejected_at"] == "input"


def test_rerun_archives_the_previous_subset(tmp_path):
    out = tmp_path / "o"
    out.mkdir()
    for name in ("result.json", "standardized.h5ad", "dropped_cells.csv.gz"):
        (out / name).write_text(name)
    archive_outputs(str(out), "subset", str(tmp_path / "in.h5ad"))
    assert len(list((out / ".history").glob("*/dropped_cells.csv.gz"))) == 1


class _NoModel:
    def classify(self, state):
        raise AssertionError("the model must not be called on a product with raw counts in X")


def test_raw_counts_in_x_are_rejected_by_probe_and_identify_columns(tmp_path):
    from eca_pp.identify_columns.cli import main as identify
    from eca_pp.probe.cli import main as probe

    src = make_integration_h5ad(tmp_path / "ok.h5ad")
    A = ad.read_h5ad(src)
    A.X = A.layers["counts"].copy()  # what the 2026-09-07 hand-made subset wrote
    raw = tmp_path / "raw.h5ad"
    A.write_h5ad(raw)
    code = probe([str(raw), "-o", str(tmp_path / "p"), "--batch-col", "batch"])
    res = json.loads((tmp_path / "p" / "result.json").read_text())
    assert code == EXIT_REJECTED and "raw counts" in res["reasons"][0]
    code = identify([str(raw), "-o", str(tmp_path / "i")], classifier=_NoModel())
    res = json.loads((tmp_path / "i" / "result.json").read_text())
    assert code == EXIT_REJECTED and res["status"] == "rejected" and "raw counts" in res["reasons"][0]
