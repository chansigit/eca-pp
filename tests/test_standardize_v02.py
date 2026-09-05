"""Acceptance tests, v0.2 set (spec §9): F4a species ladder, F4 harmonize+drop,
F5 authoritative QC, F7 standardized.h5ad.

Run on a compute node:  bash run.sh test tests -q
"""

from __future__ import annotations

import anndata as ad
import numpy as np
import pytest
import scipy.sparse as sp
from dsets import G, N, lognorm, make_counts, ref_genes, run_cli, write_h5ad

import eca_pp.standardize.species as species_ladder
from eca_pp.core.result import EXIT_BLOCKED, EXIT_OK, EXIT_REJECTED

_T1_NONE = {"species": None, "confidence": 0.0, "evidence": {"rule": "insufficient"}}


def _junk(n, tag="FAKEGENE"):
    return [f"{tag}{i:05d}" for i in range(n)]


# ------------------------------------------------------------- F4a · species

def test_species_inferred_human_from_ids(tmp_path):
    src = write_h5ad(tmp_path / "s.h5ad", sp.csr_matrix(make_counts()))
    code, res = run_cli(tmp_path, src)
    assert code == EXIT_OK
    assert res["species"]["resolved"] == "human"
    assert res["species"]["code"] == "hs"
    assert res["species"]["source"] == "inferred"
    assert res["species"]["confidence"] >= 0.99

def test_species_inferred_mouse_end_to_end(tmp_path):
    src = write_h5ad(tmp_path / "s.h5ad", sp.csr_matrix(make_counts()),
                     species="mouse")
    code, res = run_cli(tmp_path, src)
    assert code == EXIT_OK
    assert res["species"]["resolved"] == "mouse"
    assert res["metrics"]["harmonization"]["genes_dropped_frac"] < 0.05

def test_species_inferred_from_symbols_only(tmp_path):
    src = write_h5ad(tmp_path / "s.h5ad", sp.csr_matrix(make_counts()),
                     with_ids=False)
    code, res = run_cli(tmp_path, src)
    assert code == EXIT_OK
    assert res["species"]["resolved"] == "human"
    assert res["species"]["source"] == "inferred"

def test_species_conflict_blocks_with_evidence(tmp_path):
    names = [f"ENSG{i:011d}" for i in range(G // 2)] + \
            [f"ENSMUSG{i:011d}" for i in range(G - G // 2)]
    src = write_h5ad(tmp_path / "s.h5ad", sp.csr_matrix(make_counts()),
                     var_names=names, with_ids=False)
    code, res = run_cli(tmp_path, src)
    assert code == EXIT_BLOCKED and res["status"] == "needs_review"
    assert res["species"]["resolved"] is None
    hits = res["species"]["evidence"]["ensembl_prefix_hits"]
    assert set(hits) == {"human", "mouse"}
    assert any("--species" in r for r in res["reasons"])

def test_species_cli_override(tmp_path):
    src = write_h5ad(tmp_path / "s.h5ad", sp.csr_matrix(make_counts()))
    code, res = run_cli(tmp_path, src, "--species", "hs")
    assert code == EXIT_OK
    assert res["species"]["source"] == "cli"
    assert res["species"]["confidence"] == 1.0

def test_llm_tier_used_only_with_flag(tmp_path, monkeypatch):
    src = write_h5ad(tmp_path / "s.h5ad", sp.csr_matrix(make_counts()))
    monkeypatch.setattr(species_ladder, "_t1_infer", lambda adata: dict(_T1_NONE))
    monkeypatch.setattr(species_ladder, "_llm_infer",
                        lambda symbols, evidence: ("human", 0.9))
    code, res = run_cli(tmp_path, src, "--llm")
    assert code == EXIT_OK
    assert res["species"]["source"] == "llm"
    assert res["species"]["confidence"] == 0.9

def test_t1_failure_without_llm_blocks(tmp_path, monkeypatch):
    src = write_h5ad(tmp_path / "s.h5ad", sp.csr_matrix(make_counts()))
    monkeypatch.setattr(species_ladder, "_t1_infer", lambda adata: dict(_T1_NONE))
    code, res = run_cli(tmp_path, src)
    assert code == EXIT_BLOCKED
    assert res["species"]["resolved"] is None

def test_llm_failure_falls_through_to_block(tmp_path, monkeypatch):
    src = write_h5ad(tmp_path / "s.h5ad", sp.csr_matrix(make_counts()))
    monkeypatch.setattr(species_ladder, "_t1_infer", lambda adata: dict(_T1_NONE))
    monkeypatch.setattr(species_ladder, "_llm_infer", lambda symbols, evidence: None)
    code, res = run_cli(tmp_path, src, "--llm")
    assert code == EXIT_BLOCKED


def test_llm_species_rejects_out_of_range_confidence(monkeypatch):
    from eca_pp import agent

    monkeypatch.setattr(agent, "check_available", lambda: None)

    def fake_ask_json(**kwargs):
        validate = kwargs["validate"]
        for invalid in (80, -1, float("nan"), True):
            with pytest.raises((TypeError, ValueError), match="confidence"):
                validate({"species": "human", "confidence": invalid,
                          "reason": "test"})
        accepted = validate({"species": "human", "confidence": 0.8,
                             "reason": "human gene symbols"})
        return accepted, None, [], {}

    monkeypatch.setattr(agent, "ask_json", fake_ask_json)
    resolved = species_ladder._llm_infer(["GAPDH", "ACTB"], {})
    assert resolved is not None
    assert resolved[:2] == ("human", 0.8)


# --------------------------------------------------------- F4 · harmonize/drop

def test_alias_gene_renamed_with_provenance(tmp_path):
    syms, ids = ref_genes("human", G)
    names, gids = list(syms), list(ids)
    names[0], gids[0] = "SEPT1", ""  # previous symbol of SEPTIN1; no ID shortcut
    src = write_h5ad(tmp_path / "s.h5ad", sp.csr_matrix(make_counts()),
                     var_names=names, gene_ids=gids)
    code, res = run_cli(tmp_path, src)
    assert code == EXIT_OK
    out = ad.read_h5ad(res["output"])
    (orig_idx,) = np.where(out.var["original_feature_name"] == "SEPT1")[0]
    assert out.var_names[orig_idx].startswith("SEPTIN1")
    assert out.var["mapping_status"].iloc[orig_idx] in ("previous_symbol", "alias_symbol")

def test_unmapped_dropped_by_default(tmp_path):
    syms, ids = ref_genes("human", G)
    names = list(syms) + _junk(50)
    gids = list(ids) + [""] * 50
    src = write_h5ad(tmp_path / "s.h5ad", sp.csr_matrix(make_counts(g=G + 50)),
                     var_names=names, gene_ids=gids)
    code, res = run_cli(tmp_path, src)
    assert code == EXIT_OK and res["status"] == "ok"
    h = res["metrics"]["harmonization"]
    assert h["genes_dropped"]["unmapped"] == 50 and h["genes_kept"] == G
    out = ad.read_h5ad(res["output"])
    assert out.n_vars == G
    assert not any(n.startswith("FAKEGENE") for n in out.var_names)

def test_keep_unmapped_flag(tmp_path):
    syms, ids = ref_genes("human", G)
    names = list(syms) + _junk(50)
    gids = list(ids) + [""] * 50
    src = write_h5ad(tmp_path / "s.h5ad", sp.csr_matrix(make_counts(g=G + 50)),
                     var_names=names, gene_ids=gids)
    code, res = run_cli(tmp_path, src, "--keep-unmapped")
    assert code == EXIT_OK
    h = res["metrics"]["harmonization"]
    assert h["genes_kept"] == G + 50 and sum(h["genes_dropped"].values()) == 0
    assert h["genes_unmappable"]["unmapped"] == 50
    out = ad.read_h5ad(res["output"])
    assert sum(n.startswith("FAKEGENE") for n in out.var_names) == 50

def test_heavy_drop_flags_needs_review(tmp_path):
    n_real, n_junk = 5500, 4000  # 42% unmappable, but 5500 kept still passes the gate
    syms, ids = ref_genes("human", n_real)
    names = list(syms) + _junk(n_junk)
    gids = list(ids) + [""] * n_junk
    src = write_h5ad(tmp_path / "s.h5ad",
                     sp.csr_matrix(make_counts(g=n_real + n_junk)),
                     var_names=names, gene_ids=gids)
    code, res = run_cli(tmp_path, src)
    assert code == EXIT_OK and res["status"] == "needs_review"
    assert any("unmappable" in r for r in res["reasons"])

def test_drop_below_gate_rejects(tmp_path):
    n_real, n_junk = 4000, 3000  # passes the pre-gate (7000), fails after dropping
    syms, ids = ref_genes("human", n_real)
    names = list(syms) + _junk(n_junk)
    gids = list(ids) + [""] * n_junk
    src = write_h5ad(tmp_path / "s.h5ad",
                     sp.csr_matrix(make_counts(g=n_real + n_junk)),
                     var_names=names, gene_ids=gids)
    code, res = run_cli(tmp_path, src)
    assert code == EXIT_REJECTED and res["rejected_at"] == "final_gate"
    assert any("dropping" in r for r in res["reasons"])


# ----------------------------------------------------------------- F5 · QC

def test_qc_values_exact(tmp_path):
    names = ["MT-ND1", "HBB", "GAPDH", "ACTB"]
    counts = np.array([[10, 5, 85, 0],
                       [2, 0, 8, 90]], dtype=np.float32)
    src = write_h5ad(tmp_path / "s.h5ad", sp.csr_matrix(counts),
                     var_names=names, with_ids=False)
    code, res = run_cli(tmp_path, src, "--no-gate", "--species", "hs")
    assert code == EXIT_OK
    q = res["metrics"]["qc"]
    assert q["n_mt_genes"] == 1 and q["n_hb_genes"] == 1
    out = ad.read_h5ad(res["output"])
    np.testing.assert_allclose(out.obs["pct_counts_mt"], [10.0, 2.0])
    np.testing.assert_allclose(out.obs["pct_counts_hb"], [5.0, 0.0])
    np.testing.assert_allclose(out.obs["total_counts"], [100.0, 100.0])
    np.testing.assert_allclose(out.obs["n_genes_by_counts"], [3, 3])

def test_existing_qc_column_preserved_as_original(tmp_path):
    src = write_h5ad(tmp_path / "s.h5ad", sp.csr_matrix(make_counts()),
                     obs={"pct_counts_mt": np.full(N, 77.0)})
    code, res = run_cli(tmp_path, src)
    assert code == EXIT_OK
    assert "pct_counts_mt" in res["metrics"]["qc"]["overwritten_obs_cols"]
    out = ad.read_h5ad(res["output"])
    assert np.all(out.obs["pct_counts_mt__original"] == 77.0)
    assert not np.all(out.obs["pct_counts_mt"] == 77.0)  # authoritative recompute

def test_zero_mt_hb_hits_are_normal(tmp_path):
    # A matrix with no mito/hb genes at all (pre-filtered upstream, like Tabula
    # Muris) is NORMAL: recorded in metrics.qc, but never flagged for review.
    import stangene

    syms, ids = ref_genes("human", G + 60)
    keep = ~(stangene.mito_mask(syms, "human") | stangene.hb_mask(syms, "human"))
    names = [s for s, k in zip(syms, keep) if k][:G]
    gids = [i for i, k in zip(ids, keep) if k][:G]
    src = write_h5ad(tmp_path / "s.h5ad", sp.csr_matrix(make_counts()),
                     var_names=names, gene_ids=gids)
    code, res = run_cli(tmp_path, src)
    assert code == EXIT_OK and res["status"] == "ok"
    q = res["metrics"]["qc"]
    assert q["n_mt_genes"] == 0 and q["n_hb_genes"] == 0
    assert not res["reasons"]


# ------------------------------------------------------------- F7 · standard form

def test_standard_form_contract(tmp_path):
    c = make_counts()
    src = write_h5ad(tmp_path / "s.h5ad", lognorm(c), {"RNA_raw": sp.csr_matrix(c)})
    code, res = run_cli(tmp_path, src)
    assert code == EXIT_OK
    out = ad.read_h5ad(res["output"])
    # I1: single counts layer (the odd-named source layer was renamed away);
    # anndata >=0.13 additionally exposes X as layers[None] — not a real layer.
    assert {k for k in out.layers.keys() if k is not None} == {"counts"}
    L = out.layers["counts"]
    Ld = np.asarray(L.todense()) if sp.issparse(L) else np.asarray(L)
    assert np.allclose(Ld, np.round(Ld))
    # I2: X = log1p(normalize_total(counts, 1e4)), float32
    assert out.X.dtype == np.float32
    row = np.asarray(out.X[0].todense()).ravel() if sp.issparse(out.X) \
        else np.asarray(out.X[0]).ravel()
    crow = Ld[0]
    expect = np.log1p(crow / crow.sum() * 1e4)
    np.testing.assert_allclose(row, expect, rtol=1e-5, atol=1e-5)
    # I3/I8: provenance in var + uns
    for col in ("original_feature_name", "mapping_status", "gene_symbol_harmonized"):
        assert col in out.var.columns
    prov = out.uns["eca_pp_standardize"]
    assert prov["counts_source"] == "layer:RNA_raw"
    assert prov["species"] == "human"
    # I4: the four QC columns
    for col in ("pct_counts_mt", "pct_counts_hb", "total_counts", "n_genes_by_counts"):
        assert col in out.obs.columns

def test_rerun_is_reproducible(tmp_path):
    src = write_h5ad(tmp_path / "s.h5ad", sp.csr_matrix(make_counts()))
    _, res1 = run_cli(tmp_path / "a", src)
    _, res2 = run_cli(tmp_path / "b", src)
    for r in (res1, res2):
        r.pop("finished_at")
        r.pop("output")  # differs by outdir only
        r["metrics"].pop("timings")  # wall-clock, varies by nature
    assert res1 == res2

def test_no_stale_tmp_on_success(tmp_path):
    src = write_h5ad(tmp_path / "s.h5ad", sp.csr_matrix(make_counts()))
    code, res = run_cli(tmp_path, src)
    assert code == EXIT_OK
    leftovers = list((tmp_path / "out").glob("*.tmp"))
    assert leftovers == []


# ------------------------------------------------------ F1b · raw expansion

def _hvg_object(tmp_path, *, raw_from_same_counts=True, raw_subset=False):
    """scanpy-style object: X = scaled 2000-HVG subset, layers[counts] = HVG
    counts, raw = log-normalized FULL gene space (6000 genes)."""
    c = make_counts()  # N x G integer
    syms, ids = ref_genes("human", G)
    full = ad.AnnData(X=lognorm(c if raw_from_same_counts else make_counts()))
    full.var_names = np.array(syms, dtype=object)
    full.var["gene_ids"] = ids
    hvg = np.sort(np.random.default_rng(1).choice(G, 2000, replace=False))
    sub = full[:, hvg].copy()
    Xs = np.asarray(sub.X.todense())
    sub.X = ((Xs - Xs.mean(0)) / (Xs.std(0) + 1e-6)).astype(np.float32)
    sub.layers["counts"] = sp.csr_matrix(c[:, hvg])
    sub.raw = full[:, :G - 100].copy() if raw_subset else full
    path = tmp_path / "hvg.h5ad"
    sub.write_h5ad(path)
    return path


def test_raw_with_more_genes_is_preferred(tmp_path):
    src = _hvg_object(tmp_path)
    code, res = run_cli(tmp_path, src)
    assert code == EXIT_OK and res["status"] == "ok", res["reasons"]
    exp = res["metrics"]["raw_expansion"]
    assert exp["applied"] and exp["n_vars_x"] == 2000 and exp["n_vars_raw"] == G
    assert exp["dropped_layers"] == ["counts"]
    assert res["metrics"]["counts_source"].startswith("recovered")
    out = ad.read_h5ad(res["output"])
    assert out.n_vars >= G - 10 and out.raw is None
    assert {k for k in out.layers.keys() if k is not None} == {"counts"}
    assert out.uns["eca_pp_standardize"]["raw_expanded"] == "true"


def test_raw_is_trusted_over_the_hvg_counts_layer(tmp_path):
    """raw built from other counts is still taken as-is: no cross-check."""
    src = _hvg_object(tmp_path, raw_from_same_counts=False)
    code, res = run_cli(tmp_path, src)
    assert code == EXIT_OK and res["status"] == "ok"
    assert res["metrics"]["raw_expansion"]["applied"]
    assert "counts_check" not in res["metrics"]["raw_expansion"]


def test_raw_expansion_can_be_disabled(tmp_path):
    src = _hvg_object(tmp_path)
    code, res = run_cli(tmp_path, src, "--no-raw-expand", "--min-genes", "1000")
    assert code == EXIT_OK
    exp = res["metrics"]["raw_expansion"]
    assert not exp["applied"] and "disabled" in exp["reason"]
    assert res["metrics"]["counts_source"] == "layer:counts"
    assert res["metrics"]["n_vars"] == 2000


def test_raw_that_does_not_cover_x_is_left_alone(tmp_path):
    src = _hvg_object(tmp_path, raw_subset=True)
    code, res = run_cli(tmp_path, src, "--min-genes", "1000")
    exp = res["metrics"]["raw_expansion"]
    assert not exp["applied"] and "do not align" in exp["reason"]
    assert res["metrics"]["counts_source"] == "layer:counts"
