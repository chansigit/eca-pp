"""integration-probe tests (identify-columns spec §10)."""

from __future__ import annotations

import json

import numpy as np

from ecasteps.probe.cli import main
from intdata import make_integration_h5ad

RNG = np.random.default_rng(7)


def run(tmp_path, src, *extra):
    out = tmp_path / "probe"
    code = main([str(src), "-o", str(out), *map(str, extra)])
    res = json.loads((out / "result.json").read_text())
    assert res["exit_code"] == code
    return code, res


def test_true_batch_column_shows_effect_and_gain(tmp_path):
    src = make_integration_h5ad(tmp_path / "s.h5ad", effect=4.0)
    code, res = run(tmp_path, src, "--batch-col", "batch",
                    "--cell-type-col", "cell_type", "--n-cells", 600)
    assert code == 0 and res["status"] == "ok"
    m = res["metrics"]
    assert m["n_batches"] == 2
    assert m["clisi_labels"] == "annotated"
    assert m["ilisi_norm_pre"] < 0.6          # real batch effect separates
    assert m["harmony_converged"] is True
    assert m["ilisi_norm_post"] - m["ilisi_norm_pre"] >= 0.1
    assert m["clisi_norm_post"] >= m["clisi_norm_pre"] - 0.05
    assert (tmp_path / "probe" / "umap.png").exists()


def test_shuffled_labels_show_no_effect(tmp_path):
    src = make_integration_h5ad(tmp_path / "s.h5ad", effect=4.0, obs_extra={
        "shuffled": RNG.permutation(np.repeat(["b0", "b1"], 300))})
    code, res = run(tmp_path, src, "--batch-col", "shuffled", "--n-cells", 600)
    assert code == 0
    m = res["metrics"]
    assert m["ilisi_norm_pre"] > 0.8          # random labels are pre-mixed
    assert m["pc_regression_r2"] < 0.05


def test_no_effect_data_reports_premixed(tmp_path):
    src = make_integration_h5ad(tmp_path / "s.h5ad", effect=1.0)
    code, res = run(tmp_path, src, "--batch-col", "batch", "--n-cells", 600)
    assert code == 0
    m = res["metrics"]
    assert m["ilisi_norm_pre"] > 0.8
    assert m["pc_regression_r2"] < 0.05


def test_pathological_micro_batches_recorded_not_crashed(tmp_path):
    n = 600
    src = make_integration_h5ad(tmp_path / "s.h5ad", effect=1.0, obs_extra={
        "micro": np.array([f"g{i % 150}" for i in range(n)])})
    code, res = run(tmp_path, src, "--batch-col", "micro", "--n-cells", 600)
    assert code == 0 and res["status"] == "ok"
    assert isinstance(res["metrics"]["harmony_converged"], bool)


def test_single_batch_rejected(tmp_path):
    src = make_integration_h5ad(tmp_path / "s.h5ad", obs_extra={
        "mono": np.repeat("only", 600)})
    code, res = run(tmp_path, src, "--batch-col", "mono")
    assert code == 2 and res["status"] == "rejected"


def test_missing_column_rejected(tmp_path):
    src = make_integration_h5ad(tmp_path / "s.h5ad")
    code, res = run(tmp_path, src, "--batch-col", "nope")
    assert code == 2 and any("nope" in r for r in res["reasons"])


def test_tsv_spec_and_pseudo_labels(tmp_path):
    src = make_integration_h5ad(tmp_path / "s.h5ad", effect=4.0)
    import anndata as ad
    A = ad.read_h5ad(src)
    tsv = tmp_path / "batch.tsv"
    lines = ["cell_id\tvalue"] + [f"{c}\t{b}" for c, b in
                                  zip(A.obs_names, A.obs["batch"])]
    tsv.write_text("\n".join(lines) + "\n")
    code, res = run(tmp_path, src, "--batch-col", tsv, "--n-cells", 600)
    assert code == 0
    assert res["metrics"]["n_batches"] == 2
    assert res["metrics"]["clisi_labels"] == "pseudo"  # no --cell-type-col
    assert res["metrics"]["n_cell_types"] >= 2


def test_reproducible(tmp_path):
    src = make_integration_h5ad(tmp_path / "s.h5ad", effect=4.0)
    _, r1 = run(tmp_path / "a", src, "--batch-col", "batch", "--n-cells", 500)
    _, r2 = run(tmp_path / "b", src, "--batch-col", "batch", "--n-cells", 500)
    assert r1["metrics"]["ilisi_norm_pre"] == r2["metrics"]["ilisi_norm_pre"]
    assert r1["metrics"]["ilisi_norm_post"] == r2["metrics"]["ilisi_norm_post"]
