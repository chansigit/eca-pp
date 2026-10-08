"""identify-columns tests (identify-columns spec §10) — the model is always a
test double here; the deterministic spine (profile → evidence → one
classification → ordered probes → verdict) is what's under test.
"""

from __future__ import annotations

import json
from pathlib import Path

import numpy as np
import pytest
from intdata import make_integration_h5ad

from eca_pp.core.atomic_io import copyfile_atomic
from eca_pp.identify_columns import obsprofile
from eca_pp.identify_columns.cli import (
    build_candidates,
    build_evidence,
    classify_column,
    main,
    qualifies,
    trial_verdict,
)
from eca_pp.identify_columns.policies import HeuristicClassifier, PolicyUnavailable


class ScriptedClassifier:
    """Returns one fixed classification (the mock 'model')."""

    def __init__(self, batch_ranked=(), cell_type=None, reason="scripted"):
        self.answer = {
            "batch_ranked": [{"column": c, "class": "technical",
                              "reason": f"scripted: {c}"} for c in batch_ranked],
            "cell_type": cell_type, "cell_type_reason": reason,
            "columns": {}, "notes": "scripted classifier"}
        self.seen_states = []

    def classify(self, state):
        self.seen_states.append(state)
        return dict(self.answer)


def run(tmp_path, src, classifier, *extra):
    out = tmp_path / "idc"
    code = main([str(src), "-o", str(out), *map(str, extra)], classifier=classifier)
    res = json.loads((out / "result.json").read_text())
    assert res["exit_code"] == code
    return code, res, out


def _bare(tmp_path, src):
    """Strip batch / cell_type / barcode structure from a synthetic h5ad."""
    import anndata as ad
    A = ad.read_h5ad(src)
    for col in ("batch", "cell_type"):
        if col in A.obs:
            del A.obs[col]
    A.obs_names = [f"C{i:06d}" for i in range(A.n_obs)]
    dst = tmp_path / "bare.h5ad"
    A.write_h5ad(dst)
    return dst


# ------------------------------------------------------------- basic flows

def test_no_probe_profiles_and_succeeds_with_null_batch(tmp_path):
    src = make_integration_h5ad(tmp_path / "s.h5ad")
    code, res, _ = run(tmp_path, src, None, "--no-probe")
    assert code == 0 and res["status"] == "ok"
    assert res["profile"]["columns"]
    assert any(c["label"] == "batch" for c in res["candidates"]["batch"])
    assert res["columns"]["batch"] is None
    assert res["columns"]["cell_type"]["value"] == "cell_type"
    assert res["classification"]["source"] == "deterministic"
    assert {w["code"] for w in res["warnings"]} >= {"probe_disabled", "agent_unavailable"}


def test_probe_internal_error_is_not_a_scientific_rejection(monkeypatch, tmp_path):
    from eca_pp.probe import cli as probe_cli

    def fail_probe(argv):
        outdir = Path(argv[argv.index("-o") + 1])
        outdir.mkdir(parents=True)
        (outdir / "result.json").write_text(json.dumps({
            "status": "error", "reasons": ["leiden unavailable"], "metrics": {}}))
        return 1

    monkeypatch.setattr(probe_cli, "main", fail_probe)
    src = make_integration_h5ad(tmp_path / "s.h5ad")
    code, res, _ = run(tmp_path, src, ScriptedClassifier(["batch"], "cell_type"),
                       "--n-cells", 600)
    assert code == 1 and res["status"] == "error"
    assert any("leiden unavailable" in reason for reason in res["reasons"])
    assert res["trials"] == []


def test_scripted_adopt_flow(tmp_path):
    src = make_integration_h5ad(tmp_path / "s.h5ad", effect=4.0)
    clf = ScriptedClassifier(["batch"], "cell_type", "values typeA/typeB")
    code, res, out = run(tmp_path, src, clf, "--n-cells", 600)
    assert code == 0 and res["status"] == "ok"
    assert res["columns"]["batch"]["value"] == "batch"
    assert res["columns"]["batch"]["correction"] == "recommended"
    assert "probe:" in res["columns"]["batch"]["evidence"]
    assert res["columns"]["cell_type"]["value"] == "cell_type"
    (trial,) = res["trials"]
    assert trial["verdict"] == "adopted" and trial["metrics"]["harmony_converged"] is True
    assert trial["metrics"]["clisi_labels"] == "annotated"
    assert res["classification"]["batch_ranked"][0]["column"] == "batch"
    assert len(res["decisions"]) == 1 and res["decisions"][0]["action"] == "classify"
    assert "classification" in res["metrics"]["timings"]
    assert not list(out.glob("*.png"))
    # the classifier saw the evidence table with value counts
    evidence = clf.seen_states[0]["evidence"]
    batch_row = next(r for r in evidence["columns"] if r["column"] == "batch")
    assert batch_row["value_counts"] == {"b0": 300, "b1": 300}
    assert batch_row["probeable_as_batch"] is True
    assert "batch" in evidence["probeable_batch_columns"]
    assert "cell_type" not in evidence["probeable_batch_columns"]


def test_heuristic_classifier_adopts_true_batch(tmp_path):
    src = make_integration_h5ad(tmp_path / "s.h5ad", effect=4.0)
    code, res, _ = run(tmp_path, src, HeuristicClassifier(), "--n-cells", 600)
    assert code == 0
    assert res["columns"]["batch"]["value"] == "batch"
    assert res["columns"]["batch"]["correction"] == "recommended"
    assert res["columns"]["cell_type"]["confidence"] == 0.6


def test_heuristic_classifier_concludes_unnecessary_without_effect(tmp_path):
    src = make_integration_h5ad(tmp_path / "s.h5ad", effect=1.0)
    code, res, _ = run(tmp_path, src, HeuristicClassifier(), "--n-cells", 600)
    assert code == 0
    assert res["columns"]["batch"]["correction"] == "unnecessary"
    assert "already mixed" in res["columns"]["batch"]["evidence"]


def test_ranked_candidates_are_probed_in_order_until_one_qualifies(tmp_path):
    n = 600
    src = make_integration_h5ad(tmp_path / "s.h5ad", effect=4.0, obs_extra={
        "lane_micro": np.array([f"g{i % 12}" for i in range(n)])})  # 12 groups of 50
    # "lane_micro" is a random partition on the same rung as "batch": neither a
    # gain nor pre-mixed enough with 12 groups -> rejected; then "batch" is adopted.
    clf = ScriptedClassifier(["lane_micro", "batch"], "cell_type")
    code, res, _ = run(tmp_path, src, clf, "--n-cells", 600)
    assert code == 0
    assert [t["batch_col"] for t in res["trials"]] == ["lane_micro", "batch"]
    assert res["trials"][0]["verdict"] != "adopted"
    assert res["columns"]["batch"]["value"] == "batch"


def test_probe_budget_limits_verification(tmp_path):
    src = make_integration_h5ad(tmp_path / "s.h5ad", effect=4.0)
    code, res, _ = run(tmp_path, src, ScriptedClassifier(["batch"], "cell_type"),
                       "--max-probes", 0)
    assert code == 0 and res["status"] == "ok"
    assert res["columns"]["batch"] is None and not res["trials"]
    warning = next(w for w in res["warnings"] if w["code"] == "batch_evidence_insufficient")
    assert "batch" in warning["message"]


def test_unprobeable_batch_choice_is_skipped_with_warning(tmp_path):
    src = make_integration_h5ad(tmp_path / "s.h5ad", effect=4.0)
    clf = ScriptedClassifier(["cell_type", "batch"], "cell_type")  # annotation first
    code, res, _ = run(tmp_path, src, clf, "--n-cells", 600)
    assert code == 0
    assert [t["batch_col"] for t in res["trials"]] == ["batch"]
    warning = next(w for w in res["warnings"] if w["code"] == "invalid_batch_choice")
    assert warning["details"]["candidates"] == ["cell_type"]


def test_empty_ranking_still_probes_every_rung(tmp_path):
    """2026-10-03 ladder: the probes, not the classifier, decide that a dataset
    has no batch; an empty ranking leaves each rung to the heuristic's pick."""
    src = make_integration_h5ad(tmp_path / "s.h5ad", effect=4.0)
    code, res, _ = run(tmp_path, src, ScriptedClassifier([], "cell_type"), "--n-cells", 600)
    assert code == 0 and res["columns"]["batch"]["value"] == "batch"
    assert res["ladder"] == [{"rung": 1, "label": "batch", "verdict": "adopted"}]


def test_ladder_probes_a_higher_rung_before_the_classifiers_first_pick(tmp_path):
    n = 600
    src = make_integration_h5ad(tmp_path / "s.h5ad", effect=4.0, obs_extra={
        "Gender": np.array(["Male", "Female"] * (n // 2))})
    clf = ScriptedClassifier(["Gender", "batch"], "cell_type")
    code, res, _ = run(tmp_path, src, clf, "--n-cells", 600)
    assert code == 0
    assert [t["batch_col"] for t in res["trials"]] == ["batch"]
    assert res["columns"]["batch"]["rung"] == 1
    assert [(s["rung"], s["label"], s["verdict"]) for s in res["ladder"]] == [
        (1, "batch", "adopted"), (3, "Gender", "not probed")]


def test_ladder_falls_through_to_the_cell_names(tmp_path):
    """No usable obs batch column: rung 2 finds the batch hidden in the names."""
    n = 600
    src = make_integration_h5ad(tmp_path / "s.h5ad", effect=4.0, barcode_batch=True,
                                obs_extra={"lane_micro": np.array([f"g{i % 12}" for i in range(n)])})
    clf = ScriptedClassifier(["lane_micro"], "cell_type")
    code, res, _ = run(tmp_path, src, clf, "--n-cells", 600)
    assert code == 0
    assert [t["batch_col"] for t in res["trials"]] == ["lane_micro", "barcode:prefix:-"]
    assert res["columns"]["batch"]["kind"] == "derived"
    assert res["columns"]["batch"]["rung"] == 2
    assert res["columns"]["library"] is None  # 600 cells fit one library


def test_platform_detection():
    from eca_pp.identify_columns.cli import detect_platform

    def prof(**cols):
        return {"columns": [_entry(name, values) for name, values in cols.items()]}

    assert detect_platform(prof(bc1_well={"A1": 5}, sublibrary={"s1": 5}))["value"] == "split-pool"
    assert detect_platform(prof(reagent={"10X-V3": 85, "10X-V2": 75}))["value"] == "droplet"
    assert detect_platform(prof(technology={"MARS-seq": 300, "": 20}))["value"] == "plate"
    assert detect_platform(prof(technology={"modified inDrop platform": 9}))["value"] == "droplet"
    assert detect_platform(prof(sample={"s1": 5}))["value"] == "unknown"
    flagged = detect_platform(prof(reagent={"10X-V3": 1}), "split-pool")
    assert flagged == {"value": "split-pool", "source": "flag", "evidence": "--platform split-pool"}


def test_split_pool_never_uses_wells_or_barcode_segments(tmp_path):
    import anndata as ad

    n = 600
    src = make_integration_h5ad(tmp_path / "s.h5ad", effect=4.0, barcode_batch=True, obs_extra={
        "bc1_well": np.array([f"A{i % 4}" for i in range(n)]),
        "bc2_wind": np.array([f"{i % 5}" for i in range(n)]),  # Parse well index, not "technical" by name
        "sublibrary": np.array([f"s{i % 3}" for i in range(n)]),
        "Gender": np.array(["Male", "Female"] * (n // 2))})
    profile = obsprofile.profile_obs(ad.read_h5ad(src))
    batch = {c["label"]: c for c in build_candidates(profile, "split-pool")["batch"]}
    for label in ("bc1_well", "bc2_wind", "sublibrary", "barcode:prefix:-"):
        assert batch[label]["excluded"] and "split-pool" in batch[label]["note"]
    assert not batch["Gender"]["excluded"] and batch["Gender"]["rung"] == 3
    droplet = {c["label"]: c for c in build_candidates(profile, "droplet")["batch"]}
    assert not droplet["barcode:prefix:-"]["excluded"]
    code, res, _ = run(tmp_path, src, None, "--no-probe", "--platform", "auto")
    assert res["platform"]["value"] == "split-pool"


def test_library_splits_a_batch_too_big_for_one_library():
    """Hua Heart: two donors of 75-85k cells on 10x; the 23 libraries live in the
    cell names ("Donor1.M1-1.<barcode>")."""
    import anndata as ad
    import pandas as pd
    import scipy.sparse as sp
    from eca_pp.identify_columns.cli import candidate_values, find_library

    n = 600
    donor = np.repeat(["Donor1", "Donor2"], n // 2)
    lib = np.array([f"M{(i % 300) // 100}" for i in range(n)])  # 3 libraries per donor
    A = ad.AnnData(X=sp.csr_matrix((n, 5), dtype=np.float32))
    A.obs_names = pd.Index([f"{d}.{m}.CELL{i:05d}" for i, (d, m) in enumerate(zip(donor, lib))],
                           dtype=object)
    A.obs["donor"] = donor
    A.obs["lane"] = np.array([f"g{i % 6}" for i in range(n)])  # crosses donors: not nested
    candidates = build_candidates(obsprofile.profile_obs(A), "droplet")
    by_label = {c["label"]: c for c in candidates["batch"]}
    donors = candidate_values(A, by_label["donor"])
    cand, values = find_library(A, candidates, donors, cap=200, floor=50)
    assert cand["label"] == "barcode:head:." and values.nunique() == 6
    # finer than the batch even when the batch would fit (a patient batch keeps its samples)
    assert find_library(A, candidates, donors, cap=300, floor=50)[0]["label"] == "barcode:head:."
    libraries = candidate_values(A, by_label["barcode:head:."])
    assert find_library(A, candidates, libraries, cap=300, floor=50) == (None, None)  # nothing finer
    assert find_library(A, candidates, donors, cap=200, floor=150) == (None, None)  # groups too small
    cand, _ = find_library(A, candidates, None, cap=200, floor=50)  # no batch: anything that fits
    assert cand["label"] == "lane"


def test_derived_barcode_batch_materialized_as_tsv(tmp_path):
    src = make_integration_h5ad(tmp_path / "s.h5ad", effect=4.0,
                                barcode_batch=True)
    clf = ScriptedClassifier(["barcode:prefix:-"], "cell_type")
    code, res, out = run(tmp_path, src, clf, "--n-cells", 600)
    assert code == 0
    assert res["columns"]["batch"]["kind"] == "derived"
    assert res["columns"]["batch"]["value"].endswith("batch.tsv")
    assert (out / "batch.tsv").exists()
    header, first = (out / "batch.tsv").read_text().splitlines()[:2]
    assert header == "cell_id\tvalue"
    assert first.split("\t")[1] in ("b0", "b1")
    derived = clf.seen_states[0]["evidence"]["derived_candidates"]
    assert any(d["label"] == "barcode:prefix:-" for d in derived)


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


# ------------------------------------------------------- degraded operation

def test_no_grouping_columns_concludes_no_batch(tmp_path):
    src = _bare(tmp_path, make_integration_h5ad(tmp_path / "s.h5ad", effect=1.0))
    code, res, _ = run(tmp_path, src, HeuristicClassifier())
    assert code == 0
    assert res["columns"]["batch"] is None
    assert any(w["code"] == "no_batch_candidate" for w in res["warnings"])
    assert res["sample_unit"]["value"] == "whole" and res["n_obs"] <= 30000  # fits one library


def test_sample_unit_is_the_library_then_the_batch_then_the_whole_dataset_or_stop():
    from eca_pp.identify_columns.cli import sample_unit
    library, batch = {"evidence": "23 libraries"}, {"label": "donor"}
    assert sample_unit(batch, library, "droplet", 160_000)["value"] == "library"
    sample = {"label": "sample_id", "evidence": "rung-1 donor column"}
    assert sample_unit(batch, library, "split-pool", 160_000, sample)["value"] == "library"
    assert sample_unit(batch, None, "split-pool", 160_000, sample)["value"] == "sample"
    assert sample_unit(None, None, "split-pool", 160_000, sample)["value"] == "sample"
    assert sample_unit(batch, None, "droplet", 160_000)["value"] == "batch"
    assert sample_unit(None, None, "split-pool", 465_000)["value"] == "whole"
    assert sample_unit(None, None, "plate", 200_000)["value"] == "whole"
    assert sample_unit(None, None, "droplet", 8_000)["value"] == "whole"
    stop = sample_unit(None, None, "unknown", 200_000)
    assert stop["value"] == "stop" and "--platform" in stop["reason"]


def test_agent_unavailable_continues_deterministically(tmp_path):
    src = _bare(tmp_path, make_integration_h5ad(tmp_path / "s.h5ad", effect=1.0))
    code, res, _ = run(tmp_path, src, None)
    assert code == 0 and res["columns"]["batch"] is None
    assert any(w["code"] == "agent_unavailable" for w in res["warnings"])


def test_classifier_failure_falls_back_to_heuristics(tmp_path):
    class FailingClassifier:
        def classify(self, state):
            raise PolicyUnavailable("temporary agent failure")

    src = make_integration_h5ad(tmp_path / "s.h5ad", effect=4.0)
    code, res, _ = run(tmp_path, src, FailingClassifier(), "--n-cells", 600)
    assert code == 0
    assert res["classification"]["source"] == "deterministic_fallback"
    assert res["columns"]["batch"]["value"] == "batch"  # heuristics still rank it
    assert any(w["code"] == "agent_failed" for w in res["warnings"])


def test_agent_timeout_is_counted_as_failed_llm_attempt(monkeypatch, tmp_path):
    from eca_pp import agent
    from eca_pp.identify_columns.policies import AgentClassifier

    monkeypatch.setattr(agent, "check_available", lambda: None)
    clf = AgentClassifier(str(tmp_path), model="test-model")
    monkeypatch.setattr(
        clf, "classify",
        lambda state: (_ for _ in ()).throw(PolicyUnavailable(
            "[identify columns] agent run exceeded 6 min (AGENT_WALL_MIN)",
            kind="timeout")))
    src = make_integration_h5ad(tmp_path / "s.h5ad", effect=4.0)
    code, res, _ = run(tmp_path, src, clf, "--n-cells", 600, "--model", "test-model")
    assert code == 0
    llm = res["metrics"]["llm"]
    assert (llm["calls"], llm["successful_calls"], llm["failed_calls"],
            llm["timeout_calls"]) == (1, 0, 1, 1)
    assert llm["failures"][0]["model"] == "test-model"
    assert res["classification"]["source"] == "deterministic_fallback"


def test_dataset_below_probe_minimum_skips_trials(tmp_path):
    """standardize accepts >=100 cells but the probe needs >=300: say so once
    instead of probing (and re-reading the h5ad for) every candidate."""
    from eca_pp.probe.cli import MIN_CELLS

    src = make_integration_h5ad(tmp_path / "s.h5ad", n_per_batch=100)
    code, res, _ = run(tmp_path, src, HeuristicClassifier())
    assert code == 0 and res["status"] == "ok"
    assert res["columns"]["batch"] is None and res["trials"] == []
    warning = next(w for w in res["warnings"] if w["code"] == "dataset_too_small_to_probe")
    assert str(MIN_CELLS) in warning["message"]
    assert res["columns"]["cell_type"]["value"] == "cell_type"


# ------------------------------------------------------ candidate building

def test_pathological_column_excluded_before_trials(tmp_path):
    n = 600
    src = make_integration_h5ad(tmp_path / "s.h5ad", effect=4.0, obs_extra={
        "micro_id": np.array([f"g{i % 200}" for i in range(n)])})
    code, res, _ = run(tmp_path, src, None, "--no-probe")
    micro = next(c for c in res["candidates"]["batch"] if c["label"] == "micro_id")
    assert micro["excluded"] and "pathological" in micro["note"]
    assert "micro_id" not in next(
        c for c in res["classification"]["batch_ranked"])["column"]


def test_equivalent_candidates_are_marked(tmp_path):
    import anndata as ad

    src = make_integration_h5ad(tmp_path / "s.h5ad", effect=4.0)
    A = ad.read_h5ad(src)
    A.obs["patient"] = A.obs["batch"].astype(str).to_numpy()
    candidates = build_candidates(obsprofile.profile_obs(A))["batch"]
    batch = next(c for c in candidates if c["label"] == "batch")
    patient = next(c for c in candidates if c["label"] == "patient")
    assert "equivalent_to" not in batch
    assert patient["equivalent_to"] == "batch"
    # the heuristic fallback ranks only one of an equivalent pair
    ranked = HeuristicClassifier().classify(
        {"candidates": {"batch": candidates, "cell_type": []}})["batch_ranked"]
    assert [b["column"] for b in ranked if b["column"] in ("batch", "patient")] == ["batch"]


def test_batch_candidates_expose_nesting_parents(tmp_path):
    import anndata as ad

    src = make_integration_h5ad(tmp_path / "s.h5ad", effect=4.0)
    A = ad.read_h5ad(src)
    A.obs["batch-lbl"] = (A.obs["batch"].astype(str) + "-"
                          + A.obs["cell_type"].astype(str)).to_numpy()
    profile = obsprofile.profile_obs(A)
    candidates = build_candidates(profile)
    composite = next(c for c in candidates["batch"] if c["label"] == "batch-lbl")
    parents = {p["column"]: p["class"] for p in composite["nested_within"]}
    assert parents == {"batch": "technical", "cell_type": "annotation"}
    assert "nested_within" not in next(c for c in candidates["batch"] if c["label"] == "batch")
    row = next(r for r in build_evidence(profile, candidates)["columns"]
               if r["column"] == "batch-lbl")
    assert row["nested_within"] == composite["nested_within"]
    assert row["heuristic_class"] == "technical"


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
    assert classify_column(_entry("cluster_annotation", named)) == "annotation"


def test_classify_ann_affix_is_annotation_but_not_channel():
    named = {"T cell": 10, "B cell": 5}
    assert classify_column(_entry("ann_major_v260516", named)) == "annotation"
    assert classify_column(_entry("cell.ann", named)) == "annotation"
    assert classify_column(_entry("ann0608", named)) == "annotation"  # date-suffixed
    assert classify_column(_entry("ann_v2", named)) == "annotation"
    assert classify_column(_entry("channel", {"c1": 10, "c2": 5})) == "technical"


def test_cell_state_columns_are_never_probeable(tmp_path):
    n = 600
    src = make_integration_h5ad(tmp_path / "s.h5ad", effect=4.0, obs_extra={
        "cell_cycle_phase": np.array(["G1", "S", "G2M"] * (n // 3))})
    code, res, _ = run(tmp_path, src, None, "--no-probe")
    assert classify_column(_entry("cell_cycle_phase", {"G1": 1, "S": 1})) == "state"
    assert "cell_cycle_phase" not in [c["label"] for c in res["candidates"]["batch"]]
    assert "cell_cycle_phase" not in [c["label"] for c in res["candidates"]["cell_type"]]


def test_sex_and_gender_are_last_rung_batch_candidates(tmp_path):
    """issue #8 kept sex out of the donor class ('Gender' classified donor was
    adopted on mouse-pansci/BAT_WT). Since 2026-10-03 (owner) sex is a rung-3
    batch candidate: an atlas aligns one cell type across sexes, and sex is
    the fallback when no technical or cell-name batch exists. It is still
    never a cell type."""
    n = 600
    src = make_integration_h5ad(tmp_path / "s.h5ad", effect=4.0, obs_extra={
        "sex": np.array(["M", "F"] * (n // 2)),
        "Gender": np.array(["Male", "Female"] * (n // 2))})
    code, res, _ = run(tmp_path, src, None, "--no-probe")
    assert classify_column(_entry("sex", {"M": 1, "F": 1})) == "sex"
    assert classify_column(_entry("Gender", {"Male": 1, "Female": 1})) == "sex"
    batch = {c["label"]: c for c in res["candidates"]["batch"]}
    assert batch["batch"]["rung"] == 1
    for col in ("sex", "Gender"):
        assert batch[col]["class"] == "sex" and batch[col]["rung"] == 3
        assert col not in [c["label"] for c in res["candidates"]["cell_type"]]


def test_float_score_is_never_annotation_even_with_a_matching_name(tmp_path):
    """A float dtype check must win over the name heuristic: a continuous
    score named like an annotation (e.g. a per-cell classifier confidence
    score called *_class or *_celltype_score) is qc_numeric, not annotation
    -- it must never become a cell_type candidate just because of its name."""
    assert classify_column(_entry("cell_type_score", {"0.5": 1, "0.9": 1}, dtype="float")) == "qc_numeric"
    assert classify_column(_entry("prediction.class", {"0.1": 1, "0.8": 1}, dtype="float")) == "qc_numeric"
    # a genuine annotation column (string dtype) is unaffected
    assert classify_column(_entry("cell_type", {"T cell": 10, "B cell": 5})) == "annotation"


def test_per_sample_cell_count_metadata_is_never_probeable(tmp_path):
    """Real corpus finding: n_cells is an author-supplied per-sample total
    copied onto every cell of that sample -- a disguised sample id, not a
    technical/donor factor or a cell type. Recurs across 59 3CA datasets
    classified 'other' with no consistent verdict."""
    n = 600
    src = make_integration_h5ad(tmp_path / "s.h5ad", effect=4.0, obs_extra={
        "n_cells": np.array([100, 200, 300] * (n // 3))})
    code, res, _ = run(tmp_path, src, None, "--no-probe")
    assert classify_column(_entry("n_cells", {"100": 1, "200": 1}, dtype="int")) == "qc_numeric"
    assert "n_cells" not in [c["label"] for c in res["candidates"]["batch"]]
    assert "n_cells" not in [c["label"] for c in res["candidates"]["cell_type"]]


def test_organ_role_is_deterministic_from_a_name_and_value_whitelist(tmp_path):
    """issue #3 part 1: organ/tissue is resolved by whitelist alone, no agent
    call -- name match (tissue/organ) and/or content match (known organ
    names) each count as evidence; a column that hits neither stays null."""
    from eca_pp.identify_columns.cli import identify_organ
    n = 600
    src = make_integration_h5ad(tmp_path / "s.h5ad", effect=4.0, obs_extra={
        "tissue": np.array(["Kidney"] * n),
        "region": np.array(["adultheart", "adultheart", "unrecognized_value"] * (n // 3)),
        "notes": np.array(["free text a", "free text b"] * (n // 2))})
    code, res, _ = run(tmp_path, src, None, "--no-probe")
    assert res["columns"]["organ"]["label"] == "tissue"
    assert res["columns"]["organ"]["confidence"] == 0.9   # name + value both hit
    assert identify_organ({"columns": [
        {"column": "region", "dtype": "string", "n_unique": 2, "is_per_cell_unique": False,
         "examples": {"adultheart": 1, "unrecognized_value": 1}}]})["label"] == "region"
    assert identify_organ({"columns": [
        {"column": "notes", "dtype": "string", "n_unique": 2, "is_per_cell_unique": False,
         "examples": {"free text a": 1, "free text b": 1}}]}) is None


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
    batch_labels = [c["label"] for c in res["candidates"]["batch"]]
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
    cluster = next(c for c in res["candidates"]["cell_type"] if c["label"] == "leiden")
    assert cluster["class"] == "cluster" and cluster["usable_for_clisi"]
    warning = next(w for w in res["warnings"] if w["code"] == "cell_type_not_found")
    assert warning["details"]["cluster_columns"] == ["leiden"]


def test_classifier_cell_type_choice_drives_trials(tmp_path):
    n = 600
    src = make_integration_h5ad(tmp_path / "s.h5ad", effect=4.0, obs_extra={
        "free_annotation": np.array(["alpha", "beta", "gamma"] * (n // 3))})
    clf = ScriptedClassifier(["batch"], "free_annotation", "alpha/beta/gamma look like types")
    code, res, _ = run(tmp_path, src, clf, "--n-cells", 600)
    assert code == 0
    (trial,) = res["trials"]
    assert trial["cell_type_col"] == "free_annotation"
    assert res["columns"]["cell_type"]["value"] == "free_annotation"
    assert "cLISI" in res["columns"]["cell_type"]["evidence"]
    assert "cell_type" in res["columns"]["cell_type"]["evidence"]  # listed as other candidate


def test_cell_type_evidence_does_not_claim_pseudo_clisi_used_annotation():
    from eca_pp.identify_columns.cli import _ct_block

    ct = {"label": "cell_type", "n_groups": 2, "usable_for_clisi": True,
          "numeric_labels": False, "heuristic_class": "annotation"}
    trials = [{"cell_type_col": "cell_type", "metrics": {"clisi_labels": "pseudo"}}]
    block = _ct_block(ct, {"cell_type_reason": "r"}, "agent",
                      {"cell_type": []}, trials)
    assert block is not None and block["confidence"] == 0.8
    assert "used as cLISI labels" not in block["evidence"]


def test_classifier_unknown_cell_type_is_ignored(tmp_path):
    src = make_integration_h5ad(tmp_path / "s.h5ad", effect=4.0)
    clf = ScriptedClassifier(["batch"], "no_such_column", "hallucinated")
    code, res, _ = run(tmp_path, src, clf, "--n-cells", 600)
    assert code == 0
    assert res["trials"][0]["cell_type_col"] is None  # pseudo-labels instead
    assert res["columns"]["cell_type"] is None
    assert {w["code"] for w in res["warnings"]} >= {
        "invalid_cell_type_choice", "cell_type_not_found"}


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
    assert "not used for cLISI" in ct["evidence"]


def test_unnamed_text_column_chosen_from_values_is_accepted(tmp_path):
    """A column the name heuristic cannot place but whose values are
    cell-type names may be named by the classifier; it is then used as the
    cLISI label column and reported, with a warning noting the promotion."""
    n = 600
    src = make_integration_h5ad(tmp_path / "s.h5ad", effect=4.0, obs_extra={
        "labels_v2": np.array(["proB", "CDP", "ILC2P.4"] * (n // 3))})
    import anndata as ad
    A = ad.read_h5ad(src)
    del A.obs["cell_type"]
    A.write_h5ad(src)
    clf = ScriptedClassifier(["batch"], "labels_v2", "proB, CDP, ILC2P.4 are lineages")
    code, res, _ = run(tmp_path, src, clf, "--n-cells", 600)
    assert code == 0
    listed = next(c for c in clf.seen_states[0]["candidates"]["cell_type"]
                  if c["label"] == "labels_v2")
    assert listed["class"] == "other"
    (trial,) = res["trials"]
    assert trial["cell_type_col"] == "labels_v2"
    assert trial["metrics"]["clisi_labels"] == "annotated"
    assert res["columns"]["cell_type"]["value"] == "labels_v2"
    assert any(w["code"] == "cell_type_identified_from_values" for w in res["warnings"])
    assert not any(w["code"] == "cell_type_not_found" for w in res["warnings"])


def test_condition_batch_is_adopted_with_a_warning(tmp_path):
    """2026-09-10: #1's class-based guard was reverted by the repo owner --
    cross-sample cell-type alignment matters more here than preserving
    biological variance, so a condition-classified column that genuinely
    mixes the data is adopted like any other, just flagged."""
    src = make_integration_h5ad(tmp_path / "s.h5ad", effect=4.0)
    import anndata as ad
    A = ad.read_h5ad(src)
    A.obs["stage"] = A.obs.pop("batch")
    A.write_h5ad(src)
    clf = ScriptedClassifier(["stage"], "cell_type")
    clf.answer["batch_ranked"][0]["class"] = "condition"
    code, res, _ = run(tmp_path, src, clf, "--n-cells", 600)
    assert code == 0 and res["columns"]["batch"]["value"] == "stage"
    (trial,) = res["trials"]
    assert trial["class"] == "condition" and trial["verdict"] == "adopted"
    assert any(w["code"] == "biological_batch_fallback" for w in res["warnings"])


def test_trial_verdict_is_class_blind():
    """2026-09-10: the class a candidate was given (technical/donor/condition/
    other) no longer changes its verdict -- only the metrics do."""
    good = {"harmony_converged": True, "ilisi_norm_pre": 0.1, "ilisi_norm_post": 0.4,
            "clisi_norm_pre": 0.9, "clisi_norm_post": 0.9, "clisi_labels": "annotated",
            "pc_regression_r2": 0.3}
    mixed = {**good, "ilisi_norm_pre": 0.9, "ilisi_norm_post": 0.9, "pc_regression_r2": 0.01}
    bad = {**good, "ilisi_norm_post": 0.12}
    for cls in ("technical", "donor", "condition", "other", None):
        assert trial_verdict(good, cls) == "adopted"
        assert trial_verdict(mixed, cls) == "correction_unnecessary"
        assert trial_verdict(bad, cls) == "rejected"


def test_pseudo_clisi_uses_weaker_but_still_conservative_veto():
    from eca_pp.identify_columns.cli import qualifies
    base = {"harmony_converged": True, "ilisi_norm_pre": 0.1,
            "ilisi_norm_post": 0.3, "clisi_norm_pre": 0.9,
            "clisi_norm_post": 0.8}
    assert not qualifies({**base, "clisi_labels": "annotated"})
    assert qualifies({**base, "clisi_labels": "pseudo"})
    assert not qualifies({**base, "clisi_labels": "pseudo",
                           "clisi_norm_post": 0.7})


def test_the_sample_is_the_rung1_column_when_harmony_takes_a_condition(tmp_path):
    """Owner 2026-10-08 (eca-rsi #56), PanSci heart_Prkdc: the probe adopted Age_group and rejected sample_id
    (age x sex); the batch stays Age_group, the per-sample QC unit is sample_id. Without a rung-1 column,
    or with a rung-1 batch, nothing changes."""
    import anndata as ad
    n = 600
    src = make_integration_h5ad(tmp_path / "s.h5ad", effect=4.0)
    A = ad.read_h5ad(src)
    A.obs["Age_group"] = A.obs.pop("batch")
    A.obs["sample_id"] = np.array([f"g{i % 12}" for i in range(n)])  # random: rejected
    A.write_h5ad(src)
    code, res, _ = run(tmp_path, src, ScriptedClassifier(["sample_id", "Age_group"], "cell_type"),
                       "--n-cells", 600)
    assert code == 0
    assert [(s["label"], s["verdict"]) for s in res["ladder"]] == [("sample_id", "rejected"), ("Age_group", "adopted")]
    assert res["columns"]["batch"]["label"] == "Age_group"
    assert res["columns"]["sample"]["value"] == "sample_id" and "rejected" in res["columns"]["sample"]["evidence"]
    assert res["sample_unit"]["value"] == "sample"
    _, res, _ = run(tmp_path, make_integration_h5ad(tmp_path / "t.h5ad", effect=4.0),
                    ScriptedClassifier(["batch"], "cell_type"), "--n-cells", 600)
    assert res["columns"]["sample"] is None and res["sample_unit"]["value"] == "batch"
