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
    monkeypatch.setenv("ECA_PP_AGENT_MODEL", "chosen-model")
    assert agent.model_name() == "chosen-model"


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


def test_policy_uses_validated_submit_tool(monkeypatch, tmp_path):
    from eca_pp.identify_columns.policies import AgentPolicy

    monkeypatch.setattr(agent, "check_available", lambda: None)
    captured = {}

    def fake_ask_json(**kwargs):
        captured.update(kwargs)
        decision = {
            "action": "probe", "candidate": "library", "cell_type": "cell_type",
            "reason": "technical library is the next primary candidate",
        }
        return kwargs["validate"](decision), None, [], {
            "backend": "deepseek", "model": "doubao",
        }

    monkeypatch.setattr(agent, "ask_json", fake_ask_json)
    policy = AgentPolicy(str(tmp_path))
    state = {
        "profile": {},
        "candidates": {
            "batch": [],
            "cell_type": [{"label": "cell_type", "class": "annotation"}],
        },
        "trials": [],
        "thresholds": {},
        "best_cell_type": "cell_type",
        "active_batch_tier": "primary",
        "eligible_batch_candidates": ["library"],
        "probes_left": 2,
    }
    decision = policy.decide(state)
    assert decision["action"] == "probe"
    assert decision["usage"]["backend"] == "deepseek"
    assert captured["submit_tool"] == "submit_column_decision"
    assert captured["allowed_builtin"] == ("read", "glob", "grep")


def test_policy_submit_validation_rejects_out_of_tier(monkeypatch, tmp_path):
    from eca_pp.identify_columns.policies import AgentPolicy

    monkeypatch.setattr(agent, "check_available", lambda: None)

    def fake_ask_json(**kwargs):
        kwargs["validate"]({
            "action": "probe", "candidate": "condition", "cell_type": None,
            "reason": "skip the technical tier",
        })
        raise AssertionError("unreachable")

    monkeypatch.setattr(agent, "ask_json", fake_ask_json)
    policy = AgentPolicy(str(tmp_path))
    state = {
        "profile": {}, "trials": [], "thresholds": {}, "best_cell_type": None,
        "candidates": {"batch": [], "cell_type": []},
        "active_batch_tier": "primary",
        "eligible_batch_candidates": ["library"], "probes_left": 2,
    }
    with pytest.raises(ValueError, match="currently eligible"):
        policy.decide(state)
