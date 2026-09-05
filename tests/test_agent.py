"""Offline tests for the backend-neutral submit-tool interface."""

from __future__ import annotations

import sys
from types import SimpleNamespace

import pytest

from eca_pp import agent, harness


def test_backend_specific_model_defaults(monkeypatch):
    monkeypatch.delenv("ECA_PP_AGENT_MODEL", raising=False)
    monkeypatch.delenv("MODEL", raising=False)
    monkeypatch.setenv("HARNESS", "deepseek")
    assert agent.model_name() == "doubao-seed-2-1-turbo-260628"
    monkeypatch.setenv("HARNESS", "claude")
    assert agent.model_name() == "claude-sonnet-5"
    monkeypatch.setenv("HARNESS", "openai")
    assert agent.model_name() == "doubao-seed-2-1-turbo-260628"
    monkeypatch.setenv("ECA_PP_AGENT_MODEL", "chosen-model")
    assert agent.model_name() == "chosen-model"


def test_default_backend_is_openai_turbo(monkeypatch):
    monkeypatch.delenv("HARNESS", raising=False)
    monkeypatch.delenv("ECA_PP_AGENT_MODEL", raising=False)
    monkeypatch.delenv("MODEL", raising=False)
    assert agent.backend_name() == "openai"
    assert agent.model_name() == "doubao-seed-2-1-turbo-260628"


def test_unknown_backend_fails_explicitly(monkeypatch):
    monkeypatch.setenv("HARNESS", "mystery")
    with pytest.raises(harness.AgentUnavailable, match="unknown HARNESS"):
        harness.check_available()


def test_deepseek_patch_is_read_only_and_disables_builtins(monkeypatch):
    import yaml

    from eca_pp._harness_deepseek import _render_patch

    monkeypatch.setenv("ARK_API_KEY", "test-only")
    rows = yaml.safe_load(_render_patch(
        "http://127.0.0.1:1234/mcp", (), "doubao", "doubao-test"
    ))
    inserted = rows[0]["insert"]
    assert inserted[0]["config"]["transport"] == "streamable-http"
    assert inserted[0]["config"]["failOnStartupError"] is True
    assert inserted[1]["config"]["providers"]["doubao"]["api"] == \
        "openai-completions"
    assert rows[1] == {"id": "sandbox-policy", "config": {"mode": "read-only"}}
    assert {row["id"] for row in rows[2:]} >= {
        "persistent-bash", "str-replace-editor",
    }


def test_deepseek_runtime_is_reused_for_matching_config(monkeypatch, tmp_path):
    import eca_pp._harness_deepseek as backend

    instances = []

    class FakeRuntime:
        def __init__(self, **kwargs):
            self.kwargs = kwargs
            self.started = 0
            self.closed = 0
            instances.append(self)

        def start(self):
            self.started += 1

        def close(self):
            self.closed += 1

    monkeypatch.setitem(
        sys.modules, "deepseek_harness",
        SimpleNamespace(DeepSeekHarness=FakeRuntime),
    )
    backend._close_dsh_runtime()
    kwargs = {
        "dsh_bin": "/tmp/dsh", "cwd": "/tmp/work", "root": str(tmp_path),
        "provider": "doubao", "model": "model", "effort": None,
        "system_prompt": "system", "patch_text": "patch",
    }
    try:
        first, first_home, reused, _ = backend._ensure_dsh_runtime(**kwargs)
        second, second_home, reused_again, _ = backend._ensure_dsh_runtime(**kwargs)
        assert first is second and first_home == second_home
        assert reused is False and reused_again is True
        assert first.started == 1 and first.closed == 0

        changed, _, changed_reused, _ = backend._ensure_dsh_runtime(
            **{**kwargs, "model": "other-model"})
        assert changed is not first and changed_reused is False
        assert first.closed == 1
    finally:
        backend._close_dsh_runtime()


def test_transient_backend_failure_is_retried(monkeypatch):
    import anyio

    attempts = 0

    async def flaky():
        nonlocal attempts
        attempts += 1
        if attempts == 1:
            raise RuntimeError("broken pipe during initialize")
        return 7

    monkeypatch.setattr(harness, "TRANSIENT_BACKOFF_SECONDS", 0)
    assert anyio.run(harness.retry_transient, flaky, "test") == 7
    assert attempts == 2


def test_unavailable_backend_is_not_retried():
    import anyio

    attempts = 0

    async def unavailable():
        nonlocal attempts
        attempts += 1
        raise harness.AgentUnavailable("credential missing")

    with pytest.raises(harness.AgentUnavailable, match="credential missing"):
        anyio.run(harness.retry_transient, unavailable, "test")
    assert attempts == 1


def test_timeout_is_not_retried():
    import anyio

    attempts = 0

    async def timed_out():
        nonlocal attempts
        attempts += 1
        raise harness.AgentTimeout("model request stalled")

    with pytest.raises(harness.AgentTimeout, match="stalled"):
        anyio.run(harness.retry_transient, timed_out, "test")
    assert attempts == 1


def test_ask_json_validates_inside_submit_tool(monkeypatch, tmp_path):
    seen = {}

    async def fake_run_agent(**kwargs):
        seen.update(kwargs)
        submit = kwargs["tools"][0].handler
        bad = await submit({"response_json": '{"answer": 0}'})
        assert bad["is_error"] is True
        good = await submit({"response_json": '{"answer": 7}'})
        assert good["is_error"] is False
        return harness.AgentRunResult(
            submitted=good["_submitted"],
            transcript_text="done",
            tools_used=[{"tool": kwargs["submit_tool"], "target": "..."}],
            usage={"backend": "deepseek", "model": kwargs["model"]},
        )

    monkeypatch.setattr(harness, "run_agent", fake_run_agent)
    payload, transcript, tools, usage = agent.ask_json(
        system_prompt="system",
        message="question",
        cwd=str(tmp_path),
        submit_tool="submit_answer",
        schema={"type": "object", "required": ["answer"]},
        validate=lambda value: value if value["answer"] > 0 else (
            (_ for _ in ()).throw(ValueError("answer must be positive"))
        ),
        model="test-model",
    )
    assert payload == {"answer": 7}
    assert transcript == "done"
    assert tools[0]["tool"] == "submit_answer"
    assert usage == {"backend": "deepseek", "model": "test-model"}
    assert "Finish by calling submit_answer" in seen["prompt"]


def test_ask_json_rejects_non_object(monkeypatch, tmp_path):
    async def fake_run_agent(**kwargs):
        result = await kwargs["tools"][0].handler({"response_json": "[]"})
        assert result["is_error"] is True
        assert "JSON object" in result["content"][0]["text"]
        raise harness.AgentIncompleteError("no successful submit")

    monkeypatch.setattr(harness, "run_agent", fake_run_agent)
    with pytest.raises(agent.AgentUnavailable, match="no successful submit"):
        agent.ask_json(
            system_prompt="system",
            message="question",
            cwd=str(tmp_path),
            submit_tool="submit_answer",
            schema={"type": "object"},
            validate=lambda value: value,
        )


def _state(**overrides):
    profile = {"columns": [
        {"column": "library", "dtype": "string", "n_unique": 2, "missing_frac": 0.0,
         "is_constant": False, "is_per_cell_unique": False,
         "examples": {"l1": 50, "l2": 50},
         "group_sizes": {"n_groups": 2, "min": 50, "median": 50.0, "max": 50,
                         "n_tiny": 0, "tiny_group_frac": 0.0, "tiny_cell_frac": 0.0}},
        {"column": "labels_v2", "dtype": "string", "n_unique": 3, "missing_frac": 0.0,
         "is_constant": False, "is_per_cell_unique": False,
         "examples": {"proB": 40, "CDP": 30, "ILC2P": 30}},
        {"column": "leiden", "dtype": "string", "n_unique": 3, "missing_frac": 0.0,
         "is_constant": False, "is_per_cell_unique": False,
         "examples": {"0": 40, "1": 30, "2": 30}},
        {"column": "cell_id", "dtype": "string", "n_unique": 100, "missing_frac": 0.0,
         "is_constant": False, "is_per_cell_unique": True, "examples": {}},
    ], "relations": [], "derived": [], "n_obs": 100}
    candidates = {
        "batch": [{"label": "library", "kind": "existing", "class": "technical",
                   "n_groups": 2, "excluded": False, "note": ""},
                  {"label": "micro", "kind": "existing", "class": "other",
                   "n_groups": 90, "excluded": True,
                   "note": "pathological: 88/90 groups are tiny (<25 cells)"}],
        "cell_type": [{"label": "labels_v2", "class": "other"},
                      {"label": "leiden", "class": "cluster"}]}
    state = {"profile": profile, "candidates": candidates, "evidence": {},
             "heuristic_class": {"library": "technical", "labels_v2": "other",
                                 "leiden": "cluster", "cell_id": "identifier"}}
    state.update(overrides)
    return state


def test_agent_classifier_preserves_timeout_kind(monkeypatch, tmp_path):
    from eca_pp.identify_columns.policies import AgentClassifier, PolicyUnavailable

    monkeypatch.setattr(agent, "check_available", lambda: None)

    def timeout(**kwargs):
        del kwargs
        try:
            raise harness.AgentTimeout("request timed out")
        except harness.AgentTimeout as exc:
            raise agent.AgentUnavailable(str(exc)) from exc

    monkeypatch.setattr(agent, "ask_json", timeout)
    clf = AgentClassifier(str(tmp_path))
    with pytest.raises(PolicyUnavailable) as caught:
        clf.classify(_state())
    assert caught.value.kind == "timeout"


def test_classifier_uses_validated_submit_tool(monkeypatch, tmp_path):
    from eca_pp.identify_columns.policies import AgentClassifier

    monkeypatch.setattr(agent, "check_available", lambda: None)
    captured = {}

    def fake_ask_json(**kwargs):
        captured.update(kwargs)
        answer = {"batch_ranked": [{"column": "library", "reason": "l1/l2 are lanes"}],
                  "cell_type": "labels_v2", "cell_type_reason": "proB, CDP, ILC2P",
                  "columns": {"library": "technical", "bogus": "technical",
                              "labels_v2": "annotation"}}
        return kwargs["validate"](answer), None, [], {"backend": "openai", "model": "doubao"}

    monkeypatch.setattr(agent, "ask_json", fake_ask_json)
    answer = AgentClassifier(str(tmp_path)).classify(_state())
    assert answer["batch_ranked"][0]["class"] == "technical"  # filled from candidates
    assert answer["cell_type"] == "labels_v2"
    assert answer["columns"] == {"library": "technical", "labels_v2": "annotation"}
    assert answer["usage"]["backend"] == "openai"
    assert captured["submit_tool"] == "submit_column_classification"
    assert captured["allowed_builtin"] == ()
    assert "value_counts" in captured["system_prompt"] or "value" in captured["system_prompt"]


def test_classifier_validation_rejects_bad_columns(monkeypatch, tmp_path):
    from eca_pp.identify_columns.policies import AgentClassifier

    monkeypatch.setattr(agent, "check_available", lambda: None)
    captured = {}

    def fake_ask_json(**kwargs):
        captured["validate"] = kwargs["validate"]
        return {"batch_ranked": [], "cell_type": None, "cell_type_reason": ""}, None, [], {}

    monkeypatch.setattr(agent, "ask_json", fake_ask_json)
    AgentClassifier(str(tmp_path)).classify(_state())
    validate = captured["validate"]
    base = {"cell_type": None, "cell_type_reason": ""}
    with pytest.raises(ValueError, match="excluded before probing"):
        validate({**base, "batch_ranked": [{"column": "micro", "reason": "x"}]})
    with pytest.raises(ValueError, match="not a probeable"):
        validate({**base, "batch_ranked": [{"column": "labels_v2", "reason": "x"}]})
    with pytest.raises(ValueError, match="repeats"):
        validate({**base, "batch_ranked": [{"column": "library", "reason": "x"},
                                           {"column": "library", "reason": "y"}]})
    with pytest.raises(ValueError, match="cluster IDs"):
        validate({"batch_ranked": [], "cell_type": "leiden", "cell_type_reason": "0/1/2"})
    with pytest.raises(ValueError, match="per-cell identifier"):
        validate({"batch_ranked": [], "cell_type": "cell_id", "cell_type_reason": "ids"})
    with pytest.raises(ValueError, match="not an obs column"):
        validate({"batch_ranked": [], "cell_type": "nope", "cell_type_reason": "x"})
    ok = validate({"batch_ranked": [{"column": "library", "reason": "lanes"}],
                   "cell_type": "labels_v2", "cell_type_reason": "proB, CDP"})
    assert ok["cell_type"] == "labels_v2"
