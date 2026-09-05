"""identify-columns tests (identify-columns spec §10) — the agent is always a
test double here; the deterministic spine (profile → pre-checks → trials →
verdict) is what's under test.
"""

from __future__ import annotations

import json
from pathlib import Path

import numpy as np
import pytest
from intdata import make_integration_h5ad

from eca_pp.core.atomic_io import copyfile_atomic
from eca_pp.identify_columns.cli import build_candidates, classify_column, main
from eca_pp.identify_columns.policies import HeuristicPolicy


class ScriptedPolicy:
    """Replays a fixed decision list (the mock 'agent')."""

    def __init__(self, decisions):
        self.decisions = list(decisions)
        self.seen_states = []

    def decide(self, state):
        self.seen_states.append(state)
        return self.decisions.pop(0)


def run(tmp_path, src, policy, *extra):
    out = tmp_path / "idc"
    code = main([str(src), "-o", str(out), *map(str, extra)], policy=policy)
    res = json.loads((out / "result.json").read_text())
    assert res["exit_code"] == code
    return code, res, out


def test_no_probe_profiles_and_succeeds_with_null_batch(tmp_path):
    src = make_integration_h5ad(tmp_path / "s.h5ad")
    code, res, _ = run(tmp_path, src, None, "--no-probe")
    assert code == 0 and res["status"] == "ok"
    assert res["profile"]["columns"]
    assert any(c["label"] == "batch" for c in res["candidates"]["batch"])
    assert res["columns"]["batch"] is None
    assert any(w["code"] == "probe_disabled" for w in res["warnings"])


def test_probe_internal_error_is_not_a_scientific_rejection(monkeypatch, tmp_path):
    from eca_pp.probe import cli as probe_cli

    def fail_probe(argv):
        outdir = Path(argv[argv.index("-o") + 1])
        outdir.mkdir(parents=True)
        (outdir / "result.json").write_text(json.dumps({
            "status": "error",
            "reasons": ["leiden unavailable"],
            "metrics": {},
        }))
        return 1

    monkeypatch.setattr(probe_cli, "main", fail_probe)
    src = make_integration_h5ad(tmp_path / "s.h5ad")
    policy = ScriptedPolicy([{
        "action": "probe", "candidate": "batch", "cell_type": "cell_type",
        "reason": "probe the technical batch",
    }])
    code, res, _ = run(tmp_path, src, policy, "--n-cells", 600)
    assert code == 1 and res["status"] == "error"
    assert any("leiden unavailable" in reason for reason in res["reasons"])
    assert res["trials"] == []


def test_scripted_adopt_flow(tmp_path):
    src = make_integration_h5ad(tmp_path / "s.h5ad", effect=4.0)
    policy = ScriptedPolicy([
        {"action": "probe", "candidate": "batch", "reason": "finest technical"},
        {"action": "adopt", "candidate": "batch", "cell_type": "cell_type",
         "reason": "converged with iLISI gain"},
    ])
    code, res, out = run(tmp_path, src, policy, "--n-cells", 600)
    assert code == 0 and res["status"] == "ok"
    assert res["columns"]["batch"]["value"] == "batch"
    assert res["columns"]["batch"]["correction"] == "recommended"
    assert res["columns"]["cell_type"]["value"] == "cell_type"
    (trial,) = res["trials"]
    assert trial["verdict"] == "adopted"
    assert trial["metrics"]["harmony_converged"] is True
    assert "umap" not in trial
    assert not list(out.glob("*.png"))
    # the policy saw candidates with the annotation column excluded from batch
    batch_labels = {c["label"] for c in res["candidates"]["batch"]}
    assert "cell_type" not in batch_labels


def test_policy_cannot_ignore_a_qualifying_trial(tmp_path):
    src = make_integration_h5ad(tmp_path / "s.h5ad", effect=4.0)
    policy = ScriptedPolicy([
        {"action": "probe", "candidate": "batch", "cell_type": "cell_type",
         "reason": "try technical batch"},
        {"action": "give_up", "candidate": None, "cell_type": "cell_type",
         "reason": "incorrectly ignored the successful trial"},
    ])
    code, res, _ = run(tmp_path, src, policy, "--n-cells", 600)
    assert code == 0
    assert res["columns"]["batch"]["value"] == "batch"
    assert any(w["code"] == "invalid_policy_decision"
               for w in res["warnings"])


def test_heuristic_policy_adopts_true_batch(tmp_path):
    src = make_integration_h5ad(tmp_path / "s.h5ad", effect=4.0)
    code, res, _ = run(tmp_path, src, HeuristicPolicy(), "--n-cells", 600)
    assert code == 0
    assert res["columns"]["batch"]["value"] == "batch"
    assert res["columns"]["batch"]["correction"] == "recommended"


def test_heuristic_policy_concludes_unnecessary_without_effect(tmp_path):
    src = make_integration_h5ad(tmp_path / "s.h5ad", effect=1.0)
    code, res, _ = run(tmp_path, src, HeuristicPolicy(), "--n-cells", 600)
    assert code == 0
    assert res["columns"]["batch"]["correction"] == "unnecessary"


def test_agent_metric_fast_path_skips_second_model_round(monkeypatch, tmp_path):
    from eca_pp import agent
    from eca_pp.identify_columns.policies import AgentPolicy

    monkeypatch.setattr(agent, "check_available", lambda: None)
    policy = AgentPolicy(str(tmp_path))
    calls = 0

    def choose_once(state):
        nonlocal calls
        calls += 1
        if calls > 1:
            raise AssertionError("clear probe should not trigger a second agent round")
        return {"action": "probe", "candidate": "batch",
                "cell_type": "cell_type", "reason": "technical batch"}

    monkeypatch.setattr(policy, "decide", choose_once)
    src = make_integration_h5ad(tmp_path / "s.h5ad", effect=4.0)
    code, res, _ = run(tmp_path, src, policy, "--n-cells", 600)
    assert code == 0 and calls == 1
    assert res["columns"]["batch"]["correction"] == "recommended"
    assert res["metrics"]["metric_fast_path"] is True
    assert [d["source"] for d in res["decisions"]] == [
        "agent", "metric_fast_path"]
    assert "decision_1" in res["metrics"]["timings"]
    assert "decision_2" not in res["metrics"]["timings"]


def test_metric_fast_path_keeps_ambiguous_cases_for_agent_review():
    from eca_pp.identify_columns.cli import clear_metric_decision

    trial = {
        "batch_col": "batch", "cell_type_col": "cell_type",
        "verdict": "adopted",
        "metrics": {
            "ilisi_norm_pre": 0.1, "ilisi_norm_post": 0.16,
            "clisi_norm_pre": 0.9, "clisi_norm_post": 0.9,
            "clisi_labels": "annotated", "harmony_converged": True,
        },
    }
    primary = {"tier": "primary", "missing_frac": 0.0}
    assert clear_metric_decision(trial, primary) is None  # gain is borderline
    trial["metrics"]["ilisi_norm_post"] = 0.3
    assert clear_metric_decision(trial, {"tier": "fallback"}) is None
    assert clear_metric_decision(
        trial, {"tier": "primary", "missing_frac": 0.1}) is None


def test_derived_barcode_batch_materialized_as_tsv(tmp_path):
    src = make_integration_h5ad(tmp_path / "s.h5ad", effect=4.0,
                                barcode_batch=True)
    policy = ScriptedPolicy([
        {"action": "probe", "candidate": "barcode:prefix:-",
         "reason": "batch lives in the barcode prefix"},
        {"action": "adopt", "candidate": "barcode:prefix:-",
         "cell_type": "cell_type", "reason": "gain confirmed"},
    ])
    code, res, out = run(tmp_path, src, policy, "--n-cells", 600)
    assert code == 0
    assert res["columns"]["batch"]["kind"] == "derived"
    assert res["columns"]["batch"]["value"].endswith("batch.tsv")
    assert (out / "batch.tsv").exists()
    header, first = (out / "batch.tsv").read_text().splitlines()[:2]
    assert header == "cell_id\tvalue"
    assert first.split("\t")[1] in ("b0", "b1")


def test_atomic_copy_preserves_existing_batch_on_interrupted_copy(
    monkeypatch, tmp_path
):
    from eca_pp.core import atomic_io

    src = tmp_path / "candidate.tsv"
    dst = tmp_path / "batch.tsv"
    src.write_text("new complete content")
    dst.write_text("previous complete content")

    def interrupted_copy(_src, temporary_dst):
        Path(temporary_dst).write_text("partial")
        raise OSError("simulated interrupted copy")

    monkeypatch.setattr(atomic_io.shutil, "copyfile", interrupted_copy)
    with pytest.raises(OSError, match="interrupted"):
        copyfile_atomic(src, dst)

    assert dst.read_text() == "previous complete content"
    assert not list(tmp_path.glob("*.tmp"))


def test_no_grouping_columns_concludes_no_batch(tmp_path):
    src = make_integration_h5ad(tmp_path / "s.h5ad", effect=1.0)
    import anndata as ad
    A = ad.read_h5ad(src)
    del A.obs["batch"]
    del A.obs["cell_type"]
    A.obs_names = [f"C{i:06d}" for i in range(A.n_obs)]  # no barcode structure
    src2 = tmp_path / "bare.h5ad"
    A.write_h5ad(src2)
    code, res, _ = run(tmp_path, src2, HeuristicPolicy())
    assert code == 0
    assert res["columns"]["batch"] is None


def test_agent_unavailable_continues_deterministically(tmp_path):
    src = make_integration_h5ad(tmp_path / "s.h5ad", effect=1.0)
    import anndata as ad
    A = ad.read_h5ad(src)
    del A.obs["batch"]
    del A.obs["cell_type"]
    A.obs_names = [f"C{i:06d}" for i in range(A.n_obs)]
    A.write_h5ad(src)
    code, res, _ = run(tmp_path, src, None)
    assert code == 0 and res["columns"]["batch"] is None
    assert any(w["code"] == "agent_unavailable" for w in res["warnings"])


def test_midrun_agent_failure_falls_back_without_blocking(tmp_path):
    from eca_pp.identify_columns.policies import PolicyUnavailable

    class FailingPolicy:
        def decide(self, state):
            raise PolicyUnavailable("temporary agent failure")

    src = make_integration_h5ad(tmp_path / "s.h5ad", effect=1.0)
    import anndata as ad
    A = ad.read_h5ad(src)
    del A.obs["batch"]
    A.obs_names = [f"C{i:06d}" for i in range(A.n_obs)]
    A.write_h5ad(src)
    code, res, _ = run(tmp_path, src, FailingPolicy())
    assert code == 0 and res["columns"]["batch"] is None
    assert any(w["code"] == "agent_failed" for w in res["warnings"])


def test_agent_timeout_is_counted_as_failed_llm_attempt(monkeypatch, tmp_path):
    from eca_pp import agent
    from eca_pp.identify_columns.policies import AgentPolicy, PolicyUnavailable

    monkeypatch.setattr(agent, "check_available", lambda: None)
    policy = AgentPolicy(str(tmp_path), model="test-model")
    monkeypatch.setattr(
        policy, "decide",
        lambda state: (_ for _ in ()).throw(
            PolicyUnavailable(
                "[identify columns] agent run exceeded 2 min (AGENT_WALL_MIN)",
                kind="timeout",
            )
        ),
    )
    src = make_integration_h5ad(tmp_path / "s.h5ad", effect=4.0)
    code, res, _ = run(tmp_path, src, policy, "--n-cells", 600,
                       "--model", "test-model")
    assert code == 0
    llm = res["metrics"]["llm"]
    assert llm["calls"] == 1
    assert llm["successful_calls"] == 0
    assert llm["failed_calls"] == 1
    assert llm["timeout_calls"] == 1
    assert llm["failures"][0]["model"] == "test-model"


def test_pathological_column_excluded_before_trials(tmp_path):
    n = 600
    src = make_integration_h5ad(tmp_path / "s.h5ad", effect=4.0, obs_extra={
        "micro_id": np.array([f"g{i % 200}" for i in range(n)])})
    code, res, _ = run(tmp_path, src, None, "--no-probe")
    micro = next(c for c in res["candidates"]["batch"]
                 if c["label"] == "micro_id")
    assert micro["excluded"] and "pathological" in micro["note"]


def test_probe_budget_exhaustion_succeeds_with_null_batch(tmp_path):
    src = make_integration_h5ad(tmp_path / "s.h5ad", effect=4.0)
    policy = ScriptedPolicy([
        {"action": "probe", "candidate": "batch", "reason": "try"},
    ])
    code, res, _ = run(tmp_path, src, policy, "--max-probes", 0)
    assert code == 0 and res["status"] == "ok"
    assert res["columns"]["batch"] is None
    assert res["decisions"] == []
    assert not any(w["code"] == "invalid_policy_decision"
                   for w in res["warnings"])
    assert any(w["code"] == "batch_evidence_insufficient"
               for w in res["warnings"])
    assert not res["trials"]


def test_equivalent_candidates_are_marked_and_skipped(tmp_path):
    import anndata as ad

    from eca_pp.identify_columns import obsprofile

    src = make_integration_h5ad(tmp_path / "s.h5ad", effect=4.0)
    A = ad.read_h5ad(src)
    A.obs["patient"] = A.obs["batch"].astype(str).to_numpy()
    candidates = build_candidates(obsprofile.profile_obs(A))["batch"]
    batch = next(c for c in candidates if c["label"] == "batch")
    patient = next(c for c in candidates if c["label"] == "patient")
    assert "equivalent_to" not in batch
    assert patient["equivalent_to"] == "batch"


# ---------------------------------------------------- cell-type column choice

def _entry(column, examples, dtype="string"):
    return {"column": column, "dtype": dtype, "n_unique": len(examples),
            "is_constant": len(examples) <= 1, "is_per_cell_unique": False,
            "examples": examples}


def test_classify_cluster_columns_separately_from_annotation():
    named = {"T cell": 10, "B cell": 5}
    ids = {"0": 10, "1": 5}
    assert classify_column(_entry("cell_type", named)) == "annotation"
    assert classify_column(_entry("cell_ontology_class", named)) == "annotation"
    assert classify_column(_entry("free_annotation", named)) == "annotation"
    assert classify_column(_entry("seurat_clusters", ids, "int")) == "cluster"
    assert classify_column(_entry("leiden", ids)) == "cluster"
    assert classify_column(_entry("louvain_res1", ids)) == "cluster"
    # a cluster-level annotation is still an annotation
    assert classify_column(_entry("cluster_annotation", named)) == "annotation"


def test_cell_type_ranking_prefers_annotation_over_clusters(tmp_path):
    """Seurat-style obs: seurat_clusters precedes the manual annotation."""
    n = 600
    src = make_integration_h5ad(tmp_path / "s.h5ad", effect=1.0, obs_extra={
        "seurat_clusters": np.arange(n) % 6})
    import anndata as ad
    A = ad.read_h5ad(src)
    A.obs = A.obs[["batch", "seurat_clusters", "cell_type"]]  # cluster first
    A.write_h5ad(src)
    code, res, _ = run(tmp_path, src, None, "--no-probe")
    labels = [c["label"] for c in res["candidates"]["cell_type"]]
    assert labels == ["cell_type", "seurat_clusters"]
    clus = res["candidates"]["cell_type"][1]
    assert clus["class"] == "cluster" and clus["numeric_labels"]
    assert res["columns"]["cell_type"]["value"] == "cell_type"
    assert res["columns"]["cell_type"]["confidence"] == 0.7
    # cluster columns are neither batch candidates nor composite parts
    batch_labels = [c["label"] for c in res["candidates"]["batch"]]
    assert "seurat_clusters" not in batch_labels
    assert not any("seurat_clusters" in lab for lab in batch_labels)


def test_cluster_column_is_probe_support_not_reported_as_cell_type(tmp_path):
    n = 600
    src = make_integration_h5ad(tmp_path / "s.h5ad", effect=1.0, obs_extra={
        "leiden": np.array([str(i % 4) for i in range(n)])})
    import anndata as ad
    A = ad.read_h5ad(src)
    del A.obs["cell_type"]
    A.write_h5ad(src)
    code, res, _ = run(tmp_path, src, None, "--no-probe")
    assert res["columns"]["cell_type"] is None
    cluster = next(c for c in res["candidates"]["cell_type"]
                   if c["label"] == "leiden")
    assert not cluster["output_eligible"] and cluster["usable_for_clisi"]
    warning = next(w for w in res["warnings"]
                   if w["code"] == "cell_type_not_found")
    assert warning["details"]["cluster_columns"] == ["leiden"]


def test_policy_cell_type_choice_drives_trials(tmp_path):
    """The agent's cell_type on the PROBE round is what the trial uses for
    cLISI, and it is kept through adoption (null does not reset it)."""
    n = 600
    src = make_integration_h5ad(tmp_path / "s.h5ad", effect=4.0, obs_extra={
        "free_annotation": np.array(
            ["alpha", "beta", "gamma"] * (n // 3))})
    policy = ScriptedPolicy([
        {"action": "probe", "candidate": "batch",
         "cell_type": "free_annotation", "reason": "names look like types"},
        {"action": "adopt", "candidate": "batch", "cell_type": None,
         "reason": "converged with iLISI gain"},
    ])
    code, res, _ = run(tmp_path, src, policy, "--n-cells", 600)
    assert code == 0
    assert policy.seen_states[0]["best_cell_type"] == "cell_type"
    assert policy.seen_states[1]["best_cell_type"] == "free_annotation"
    (trial,) = res["trials"]
    assert trial["cell_type_col"] == "free_annotation"
    assert res["columns"]["cell_type"]["value"] == "free_annotation"
    assert "cLISI" in res["columns"]["cell_type"]["evidence"]


def test_cell_type_evidence_does_not_claim_pseudo_clisi_used_annotation():
    from eca_pp.identify_columns.cli import _ct_block

    candidates = {"cell_type": [{
        "label": "cell_type", "class": "annotation",
        "numeric_labels": False, "usable_for_clisi": True,
    }]}
    trials = [{
        "cell_type_col": "cell_type",
        "metrics": {"clisi_labels": "pseudo"},
    }]
    block = _ct_block("cell_type", candidates, trials)
    assert block is not None
    assert "used as cLISI labels" not in block["evidence"]


def test_policy_unknown_cell_type_is_ignored(tmp_path):
    src = make_integration_h5ad(tmp_path / "s.h5ad", effect=4.0)
    policy = ScriptedPolicy([
        {"action": "probe", "candidate": "batch",
         "cell_type": "no_such_column", "reason": "hallucinated"},
        {"action": "adopt", "candidate": "batch", "cell_type": "cell_type",
         "reason": "gain"},
    ])
    code, res, _ = run(tmp_path, src, policy, "--n-cells", 600)
    assert code == 0
    assert res["trials"][0]["cell_type_col"] == "cell_type"
    assert res["columns"]["cell_type"]["value"] == "cell_type"
    assert any(w["code"] == "invalid_cell_type_choice"
               for w in res["warnings"])


def test_classify_ann_affix_is_annotation_but_not_channel():
    named = {"T cell": 10, "B cell": 5}
    assert classify_column(_entry("ann_major_v260516", named)) == "annotation"
    assert classify_column(_entry("cell.ann", named)) == "annotation"
    assert classify_column(_entry("channel", {"c1": 10, "c2": 5})) == "technical"


def test_short_ct_alias_and_constant_annotation_are_valid_output(tmp_path):
    src = make_integration_h5ad(tmp_path / "s.h5ad", effect=1.0)
    import anndata as ad
    A = ad.read_h5ad(src)
    del A.obs["cell_type"]
    A.obs["ct"] = "Embryonic stem cell"
    A.write_h5ad(src)
    code, res, _ = run(tmp_path, src, None, "--no-probe")
    assert code == 0
    ct = res["columns"]["cell_type"]
    assert ct["value"] == "ct"
    candidate = next(c for c in res["candidates"]["cell_type"]
                     if c["label"] == "ct")
    assert candidate["output_eligible"]
    assert not candidate["usable_for_clisi"]
    assert "not used for cLISI" in ct["evidence"]


def test_batch_candidates_have_strict_primary_then_fallback_tiers():
    profile = {"columns": [
        {**_entry("library", {"l1": 50, "l2": 50}),
         "missing_frac": 0.0,
         "group_sizes": {"n_groups": 2, "n_tiny": 0,
                         "tiny_group_frac": 0.0, "tiny_cell_frac": 0.0}},
        {**_entry("stage", {"E1": 50, "E2": 50}),
         "missing_frac": 0.0,
         "group_sizes": {"n_groups": 2, "n_tiny": 0,
                         "tiny_group_frac": 0.0, "tiny_cell_frac": 0.0}},
    ], "derived": []}
    candidates = build_candidates(profile)
    assert [(c["label"], c["tier"]) for c in candidates["batch"]] == [
        ("library", "primary"), ("stage", "fallback")]
    from eca_pp.identify_columns.cli import _active_batch_tier
    assert _active_batch_tier(candidates, []) == "primary"
    assert _active_batch_tier(candidates, [
        {"batch_col": "library", "verdict": "rejected"}]) == "fallback"
    assert _active_batch_tier(candidates, [
        {"batch_col": "library", "verdict": "adopted"}]) is None


def test_fallback_batch_is_allowed_but_warned(tmp_path):
    src = make_integration_h5ad(tmp_path / "s.h5ad", effect=4.0)
    import anndata as ad
    A = ad.read_h5ad(src)
    A.obs["stage"] = A.obs.pop("batch")
    A.write_h5ad(src)
    policy = ScriptedPolicy([
        {"action": "probe", "candidate": "stage", "cell_type": "cell_type",
         "reason": "no technical candidate exists"},
        {"action": "adopt", "candidate": "stage", "cell_type": "cell_type",
         "reason": "fallback preserved cell structure"},
    ])
    code, res, _ = run(tmp_path, src, policy, "--n-cells", 600)
    assert code == 0 and res["columns"]["batch"]["value"] == "stage"
    assert any(w["code"] == "biological_batch_fallback"
               for w in res["warnings"])


def test_pseudo_clisi_uses_weaker_but_still_conservative_veto():
    from eca_pp.identify_columns.cli import qualifies
    base = {"harmony_converged": True, "ilisi_norm_pre": 0.1,
            "ilisi_norm_post": 0.3, "clisi_norm_pre": 0.9,
            "clisi_norm_post": 0.8}
    assert not qualifies({**base, "clisi_labels": "annotated"})
    assert qualifies({**base, "clisi_labels": "pseudo"})
    assert not qualifies({**base, "clisi_labels": "pseudo",
                           "clisi_norm_post": 0.7})


def test_dataset_below_probe_minimum_skips_trials(tmp_path):
    """standardize accepts >=100 cells but the probe needs >=300: say so once
    instead of probing (and re-reading the h5ad for) every candidate."""
    from eca_pp.probe.cli import MIN_CELLS

    src = make_integration_h5ad(tmp_path / "s.h5ad", n_per_batch=100)
    code, res, _ = run(tmp_path, src, HeuristicPolicy())
    assert code == 0 and res["status"] == "ok"
    assert res["columns"]["batch"] is None
    assert res["trials"] == []
    warning = next(w for w in res["warnings"]
                   if w["code"] == "dataset_too_small_to_probe")
    assert str(MIN_CELLS) in warning["message"]
    assert res["columns"]["cell_type"]["value"] == "cell_type"


def test_agent_can_promote_unnamed_text_column_to_cell_type(tmp_path):
    """A column the name heuristic cannot place (class "other") but whose
    values are cell-type names may be chosen by the policy; it is then used
    as the cLISI label column and reported as the cell type."""
    n = 600
    src = make_integration_h5ad(tmp_path / "s.h5ad", effect=4.0, obs_extra={
        "ann0608": np.array(["proB", "CDP", "ILC2P.4"] * (n // 3))})
    import anndata as ad
    A = ad.read_h5ad(src)
    del A.obs["cell_type"]  # no name-recognized annotation at all
    A.write_h5ad(src)
    policy = ScriptedPolicy([
        {"action": "probe", "candidate": "batch", "cell_type": "ann0608",
         "reason": "values proB, CDP, ILC2P.4 are lineage names"},
        {"action": "adopt", "candidate": "batch", "cell_type": "ann0608",
         "reason": "gain confirmed"},
    ])
    code, res, _ = run(tmp_path, src, policy, "--n-cells", 600)
    assert code == 0
    listed = next(c for c in policy.seen_states[0]["candidates"]["cell_type"]
                  if c["label"] == "ann0608")
    assert listed["class"] == "other" and not listed["output_eligible"]
    assert policy.seen_states[0]["best_cell_type"] is None
    assert policy.seen_states[1]["best_cell_type"] == "ann0608"
    (trial,) = res["trials"]
    assert trial["cell_type_col"] == "ann0608"
    assert trial["metrics"]["clisi_labels"] == "annotated"
    assert res["columns"]["cell_type"]["value"] == "ann0608"
    assert any(w["code"] == "cell_type_identified_from_values"
               for w in res["warnings"])
    assert not any(w["code"] == "cell_type_not_found" for w in res["warnings"])


def test_batch_candidates_expose_nesting_parents(tmp_path):
    src = make_integration_h5ad(tmp_path / "s.h5ad", effect=4.0)
    import anndata as ad
    A = ad.read_h5ad(src)
    A.obs["batch-ann"] = (A.obs["batch"].astype(str) + "-"
                          + A.obs["cell_type"].astype(str)).to_numpy()
    candidates = build_candidates(obsprofile_profile(A))["batch"]
    composite = next(c for c in candidates if c["label"] == "batch-ann")
    parents = {p["column"]: p["class"] for p in composite["nested_within"]}
    assert parents == {"batch": "technical", "cell_type": "annotation"}
    assert "nested_within" not in next(c for c in candidates if c["label"] == "batch")


def obsprofile_profile(adata):
    from eca_pp.identify_columns import obsprofile
    return obsprofile.profile_obs(adata)
