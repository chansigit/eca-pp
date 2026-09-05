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

    async def fake_run(starting_agent, prompt, *, max_turns, run_config, **kwargs):
        seen.update(
            {
                "prompt": prompt,
                "max_turns": max_turns,
                "run_config": run_config,
                "runner_kwargs": kwargs,
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
    assert value.usage["reasoning_effort"] == "medium"
    assert value.usage["server_state"] is True
    assert value.usage["parallel_tool_calls"] is False
    assert value.usage["nudges"] == 0
    assert value.usage["input_tokens"] == 101
    assert value.usage["reasoning_tokens"] == 13
    assert value.usage["cache_read_tokens"] == 11
    assert seen["client"]["api_key"] == "test-only"
    assert seen["client"]["max_retries"] == 0
    assert seen["run_config"].tracing_disabled is True
    assert seen["runner_kwargs"]["auto_previous_response_id"] is True
    assert seen["agent"].model_settings.parallel_tool_calls is False
    assert seen["agent"].model_settings.reasoning.effort == "medium"
    assert seen["agent"].model_settings.store is True
    assert seen["closed"] is True


def test_openai_runner_nudges_same_response_chain(monkeypatch, tmp_path):
    agents = pytest.importorskip("agents")
    openai = pytest.importorskip("openai")
    monkeypatch.setenv("ARK_API_KEY", "test-only")
    monkeypatch.setenv("OPENAI_AGENTS_REASONING_EFFORT", "medium")
    monkeypatch.setattr(backend, "check_available", lambda: None)
    calls = []

    class FakeClient:
        def __init__(self, **_kwargs):
            pass

        async def close(self):
            pass

    class FakeResult:
        def __init__(self, final_output, response_id, usage):
            self.final_output = final_output
            self.last_response_id = response_id
            self.new_items = []
            self.context_wrapper = SimpleNamespace(usage=usage)

        def to_input_list(self):
            return [{"role": "user", "content": "full history"}]

    def usage(requests, input_tokens, output_tokens, reasoning_tokens):
        return SimpleNamespace(
            requests=requests,
            input_tokens=input_tokens,
            output_tokens=output_tokens,
            output_tokens_details=SimpleNamespace(
                reasoning_tokens=reasoning_tokens
            ),
            input_tokens_details=SimpleNamespace(
                cached_tokens=0, cache_write_tokens=0
            ),
        )

    async def fake_run(starting_agent, run_input, **kwargs):
        calls.append((run_input, kwargs, starting_agent))
        if len(calls) == 1:
            return FakeResult("I should submit next", "resp-1", usage(1, 20, 6, 4))
        tool = starting_agent.tools[0]
        await tool.on_invoke_tool(
            None, json.dumps({"response_json": '{"answer": 7}'})
        )
        return FakeResult("done", "resp-2", usage(1, 8, 3, 1))

    monkeypatch.setattr(openai, "AsyncOpenAI", FakeClient)
    monkeypatch.setattr(agents.Runner, "run", fake_run)

    async def submit(args):
        return {
            "content": [{"type": "text", "text": "accepted"}],
            "is_error": False,
            "_submitted": json.loads(args["response_json"]),
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
            allowed_builtin=(),
            label="test",
            max_buffer_size=None,
            wall_seconds=30,
        )

    value = anyio.run(exercise)
    assert value.submitted == {"answer": 7}
    assert value.usage["num_turns"] == 2
    assert value.usage["input_tokens"] == 28
    assert value.usage["output_tokens"] == 9
    assert value.usage["reasoning_tokens"] == 5
    assert value.usage["nudges"] == 1
    assert calls[0][1]["auto_previous_response_id"] is True
    assert calls[1][1]["previous_response_id"] == "resp-1"
    assert "previous turn ended" in calls[1][0][0]["content"]
    assert calls[0][2].model_settings.reasoning.effort == "medium"


def test_openai_server_state_can_be_disabled(monkeypatch, tmp_path):
    agents = pytest.importorskip("agents")
    openai = pytest.importorskip("openai")
    monkeypatch.setenv("ARK_API_KEY", "test-only")
    monkeypatch.setenv("OPENAI_AGENTS_SERVER_STATE", "0")
    monkeypatch.setenv("OPENAI_AGENTS_REASONING_EFFORT", "off")
    monkeypatch.setattr(backend, "check_available", lambda: None)
    calls = []

    class FakeClient:
        def __init__(self, **_kwargs):
            pass

        async def close(self):
            pass

    class FakeResult:
        final_output = "paused"
        last_response_id = "resp-1"
        new_items = []
        context_wrapper = SimpleNamespace(usage=SimpleNamespace(
            requests=1,
            input_tokens=10,
            output_tokens=2,
            output_tokens_details=SimpleNamespace(reasoning_tokens=0),
            input_tokens_details=SimpleNamespace(
                cached_tokens=0, cache_write_tokens=0
            ),
        ))

        def to_input_list(self):
            return [{"role": "user", "content": "full history"}]

    async def fake_run(starting_agent, run_input, **kwargs):
        calls.append((run_input, kwargs, starting_agent))
        if len(calls) == 2:
            await starting_agent.tools[0].on_invoke_tool(
                None, json.dumps({"response_json": '{"answer": 7}'})
            )
        return FakeResult()

    monkeypatch.setattr(openai, "AsyncOpenAI", FakeClient)
    monkeypatch.setattr(agents.Runner, "run", fake_run)

    async def submit(args):
        return {
            "content": [{"type": "text", "text": "accepted"}],
            "is_error": False,
            "_submitted": json.loads(args["response_json"]),
        }

    tool = harness.ToolSpec(
        "submit_answer", "submit answer", {"response_json": str}, submit
    )
    import anyio

    async def exercise():
        return await backend.run_agent(
            tools=[tool], submit_tool="submit_answer", prompt="question",
            system_prompt="system", cwd=str(tmp_path), model="doubao-test",
            effort=None, max_turns=6, allowed_builtin=(), label="test",
            max_buffer_size=None, wall_seconds=30,
        )

    value = anyio.run(exercise)
    assert value.submitted == {"answer": 7}
    assert "auto_previous_response_id" not in calls[0][1]
    assert calls[1][0][0]["content"] == "full history"
    assert calls[0][2].model_settings.reasoning is None
    assert calls[0][2].model_settings.store is False


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
