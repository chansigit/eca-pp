"""identify-columns tests (identify-columns spec §10) — the agent is always a
test double here; the deterministic spine (profile → pre-checks → trials →
verdict) is what's under test.
"""

from __future__ import annotations

import json

import numpy as np

from eca_pp.identify_columns.cli import main
from eca_pp.identify_columns.policies import HeuristicPolicy
from intdata import make_integration_h5ad


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


def test_degraded_mode_profiles_and_blocks(tmp_path):
    src = make_integration_h5ad(tmp_path / "s.h5ad")
    code, res, _ = run(tmp_path, src, None)
    assert code == 3 and res["status"] == "needs_review"
    assert res["profile"]["columns"]
    assert any(c["label"] == "batch" for c in res["candidates"]["batch"])
    assert res["columns"]["batch"] is None


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
    assert (out / "trial_1_umap.png").exists()
    # the policy saw candidates with the annotation column excluded from batch
    batch_labels = {c["label"] for c in res["candidates"]["batch"]}
    assert "cell_type" not in batch_labels


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


def test_pathological_column_excluded_before_trials(tmp_path):
    n = 600
    src = make_integration_h5ad(tmp_path / "s.h5ad", effect=4.0, obs_extra={
        "micro_id": np.array([f"g{i % 200}" for i in range(n)])})
    code, res, _ = run(tmp_path, src, None)  # degraded: just inspect candidates
    micro = next(c for c in res["candidates"]["batch"]
                 if c["label"] == "micro_id")
    assert micro["excluded"] and "pathological" in micro["note"]


def test_probe_budget_exhaustion_blocks(tmp_path):
    src = make_integration_h5ad(tmp_path / "s.h5ad", effect=4.0)
    policy = ScriptedPolicy([
        {"action": "probe", "candidate": "batch", "reason": "try"},
        {"action": "give_up", "reason": "not satisfied"},
    ])
    code, res, _ = run(tmp_path, src, policy, "--n-cells", 600)
    assert code == 3 and res["status"] == "needs_review"
    assert res["trials"]
