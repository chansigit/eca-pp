"""Offline contract tests for the OpenAI Agents SDK backend."""

from __future__ import annotations

import json
from types import SimpleNamespace

import pytest

from eca_pp import _harness_openai as backend
from eca_pp import harness


def test_openai_params_schema_is_strict():
    async def unused(args):
        return args

    spec = harness.ToolSpec("submit", "submit", {"text": str, "count": int}, unused)
    assert backend._params_schema(spec) == {
        "type": "object",
        "properties": {
            "text": {"type": "string"},
            "count": {"type": "integer"},
        },
        "required": ["text", "count"],
        "additionalProperties": False,
    }


def test_openai_backend_requires_ark_key(monkeypatch):
    monkeypatch.delenv("ARK_API_KEY", raising=False)
    with pytest.raises(harness.AgentUnavailable, match="ARK_API_KEY"):
        backend.check_available()


def test_openai_runner_validates_then_stops_on_valid_submit(monkeypatch, tmp_path):
    agents = pytest.importorskip("agents")
    openai = pytest.importorskip("openai")
    monkeypatch.setenv("ARK_API_KEY", "test-only")
    monkeypatch.setattr(backend, "check_available", lambda: None)
    seen = {}

    class FakeClient:
        def __init__(self, **kwargs):
            seen["client"] = kwargs

        async def close(self):
            seen["closed"] = True

    async def fake_run(starting_agent, prompt, *, max_turns, run_config):
        seen.update(
            {
                "prompt": prompt,
                "max_turns": max_turns,
                "run_config": run_config,
                "agent": starting_agent,
            }
        )
        tool = starting_agent.tools[0]
        bad = await tool.on_invoke_tool(
            None, json.dumps({"response_json": '{"answer": 0}'})
        )
        assert "invalid response" in bad
        assert "value" not in seen
        good = await tool.on_invoke_tool(
            None, json.dumps({"response_json": '{"answer": 7}'})
        )
        result = SimpleNamespace(tool=tool, output=good)
        stop = await starting_agent.tool_use_behavior(None, [result])
        assert stop.is_final_output is True
        usage = SimpleNamespace(
            requests=2,
            input_tokens=101,
            output_tokens=17,
            output_tokens_details=SimpleNamespace(reasoning_tokens=13),
            input_tokens_details=SimpleNamespace(
                cached_tokens=11, cache_write_tokens=3
            ),
        )
        return SimpleNamespace(
            final_output=stop.final_output,
            new_items=[],
            context_wrapper=SimpleNamespace(usage=usage),
        )

    monkeypatch.setattr(openai, "AsyncOpenAI", FakeClient)
    monkeypatch.setattr(agents.Runner, "run", fake_run)

    async def submit(args):
        payload = json.loads(args["response_json"])
        if payload["answer"] <= 0:
            return {
                "content": [{"type": "text", "text": "invalid response"}],
                "is_error": True,
            }
        return {
            "content": [{"type": "text", "text": "response accepted"}],
            "is_error": False,
            "_submitted": payload,
        }

    tool = harness.ToolSpec(
        "submit_answer", "submit answer", {"response_json": str}, submit
    )
    import anyio

    async def exercise():
        return await backend.run_agent(
            tools=[tool],
            submit_tool="submit_answer",
            prompt="question",
            system_prompt="system",
            cwd=str(tmp_path),
            model="doubao-test",
            effort=None,
            max_turns=6,
            allowed_builtin=("read", "glob", "grep"),
            label="test",
            max_buffer_size=None,
            wall_seconds=30,
        )

    value = anyio.run(exercise)
    assert value.submitted == {"answer": 7}
    assert value.tools_used == [
        {"tool": "submit_answer", "target": '{"answer": 0}'},
        {"tool": "submit_answer", "target": '{"answer": 7}'},
    ]
    assert value.usage["backend"] == "openai"
    assert value.usage["input_tokens"] == 101
    assert value.usage["reasoning_tokens"] == 13
    assert value.usage["cache_read_tokens"] == 11
    assert seen["client"]["api_key"] == "test-only"
    assert seen["run_config"].tracing_disabled is True
    assert seen["closed"] is True


def test_openai_runner_enforces_whole_run_timeout(monkeypatch, tmp_path):
    agents = pytest.importorskip("agents")
    openai = pytest.importorskip("openai")
    monkeypatch.setenv("ARK_API_KEY", "test-only")
    monkeypatch.setattr(backend, "check_available", lambda: None)
    seen = {}

    class FakeClient:
        def __init__(self, **_kwargs):
            pass

        async def close(self):
            seen["closed"] = True

    async def slow_run(*_args, **_kwargs):
        import asyncio

        await asyncio.sleep(1)

    monkeypatch.setattr(openai, "AsyncOpenAI", FakeClient)
    monkeypatch.setattr(agents.Runner, "run", slow_run)

    async def submit(args):
        return {"content": [], "is_error": False, "_submitted": args}

    tool = harness.ToolSpec("submit", "submit", {"value": str}, submit)
    import anyio

    async def exercise():
        await backend.run_agent(
            tools=[tool],
            submit_tool="submit",
            prompt="question",
            system_prompt="system",
            cwd=str(tmp_path),
            model="doubao-test",
            effort=None,
            max_turns=6,
            allowed_builtin=(),
            label="test",
            max_buffer_size=None,
            wall_seconds=0.001,
        )

    with pytest.raises(harness.AgentTimeout, match="wall-clock limit"):
        anyio.run(exercise)
    assert seen["closed"] is True
