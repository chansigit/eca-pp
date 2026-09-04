"""``HARNESS=openai`` backend using the OpenAI Agents SDK with Doubao.

Unlike DSH, the SDK can invoke eca-pp's Python submit handlers directly.  No
CLI subprocess, MCP bridge, listening port, or persistent session is needed.
The runner stops immediately after a *valid* submit call; an invalid call is
returned to the model so it can correct and resubmit in the same run.
"""

from __future__ import annotations

import asyncio
import json
import os
import time
from importlib.metadata import PackageNotFoundError, version

from .harness import (
    AgentIncompleteError,
    AgentRunResult,
    AgentTimeout,
    AgentUnavailable,
    ToolSpec,
)

MIN_OPENAI_AGENTS_SDK = (0, 22, 0)
DOUBAO_BASE_URL_DEFAULT = "https://ark.cn-beijing.volces.com/api/v3"


def _version_tuple(value: str) -> tuple[int, ...]:
    result = []
    for part in value.split("."):
        digits = "".join(ch for ch in part if ch.isdigit())
        if not digits:
            break
        result.append(int(digits))
    return tuple(result)


def check_available() -> None:
    if not os.environ.get("ARK_API_KEY"):
        raise AgentUnavailable("HARNESS=openai needs ARK_API_KEY")
    try:
        installed = version("openai-agents")
    except PackageNotFoundError:
        raise AgentUnavailable(
            "HARNESS=openai needs openai-agents>=0.22,<1; install eca-pp[openai]"
        ) from None
    if _version_tuple(installed) < MIN_OPENAI_AGENTS_SDK:
        floor = ".".join(map(str, MIN_OPENAI_AGENTS_SDK))
        raise AgentUnavailable(
            f"openai-agents {installed} is too old; install >={floor},<1"
        )
    try:
        import agents  # noqa: F401
        import openai  # noqa: F401
    except ImportError as exc:
        raise AgentUnavailable(
            f"HARNESS=openai dependencies are missing: {exc}"
        ) from None


_JSON_TYPE = {
    str: "string",
    int: "integer",
    float: "number",
    bool: "boolean",
    dict: "object",
    list: "array",
}


def _params_schema(spec: ToolSpec) -> dict:
    properties = {}
    for name, kind in spec.input_schema.items():
        json_type = _JSON_TYPE.get(kind)
        if json_type is None:
            raise TypeError(
                f"HARNESS=openai cannot expose tool parameter {name!r} "
                f"with Python type {kind!r}"
            )
        properties[name] = {"type": json_type}
    return {
        "type": "object",
        "properties": properties,
        "required": list(properties),
        "additionalProperties": False,
    }


def _usage_dict(
    result, model: str | None, runtime_init: float, agent_run: float, total: float
) -> dict:
    raw = getattr(getattr(result, "context_wrapper", None), "usage", None)
    input_details = getattr(raw, "input_tokens_details", None)
    output_details = getattr(raw, "output_tokens_details", None)
    return {
        "backend": "openai",
        "model": model,
        "cost_usd": None,
        "input_tokens": getattr(raw, "input_tokens", None),
        "output_tokens": getattr(raw, "output_tokens", None),
        "reasoning_tokens": getattr(output_details, "reasoning_tokens", None),
        "cache_creation_tokens": getattr(input_details, "cache_write_tokens", None),
        "cache_read_tokens": getattr(input_details, "cached_tokens", None),
        "num_turns": getattr(raw, "requests", None),
        "timings": {
            "runtime_init": round(runtime_init, 3),
            "agent_run": round(agent_run, 3),
            "total": round(total, 3),
        },
    }


async def run_agent(
    *,
    tools: list[ToolSpec],
    submit_tool: str,
    prompt: str,
    system_prompt: str | None,
    cwd: str,
    model: str | None,
    effort: str | None,
    max_turns: int,
    allowed_builtin: tuple[str, ...],
    label: str,
    max_buffer_size: int | None,
    wall_seconds: float | None = None,
) -> AgentRunResult:
    check_available()
    runtime_started = time.perf_counter()
    from agents import (
        Agent,
        FunctionTool,
        ItemHelpers,
        ModelSettings,
        RunConfig,
        Runner,
        ToolsToFinalOutputResult,
    )
    from agents.exceptions import MaxTurnsExceeded
    from agents.models.openai_responses import OpenAIResponsesModel
    from openai import APITimeoutError, AsyncOpenAI

    os.makedirs(cwd, exist_ok=True)
    submitted: dict = {}
    tools_used: list[dict] = []

    def wrap(spec: ToolSpec):
        is_submit = spec.name == submit_tool

        async def invoke(_context, arguments_json: str):
            try:
                arguments = json.loads(arguments_json)
                if not isinstance(arguments, dict):
                    raise TypeError("tool arguments must be a JSON object")
            except (json.JSONDecodeError, TypeError, ValueError) as exc:
                return f"invalid tool arguments, fix and retry: {exc}"

            target = str(next(iter(arguments.values()), ""))[:200]
            tools_used.append({"tool": spec.name, "target": target})
            print(f"== [{label}] agent: {spec.name}({target[:80]})", flush=True)
            result = await spec.handler(arguments)
            content = " ".join(
                str(item.get("text", "")) for item in result.get("content", [])
            ).strip()
            if result.get("is_error"):
                print(
                    f"== [{label}] tool error in {spec.name}: {content[:200]!r}",
                    flush=True,
                )
            elif is_submit:
                submitted["value"] = result.get("_submitted", arguments)
            return content or ("tool failed" if result.get("is_error") else "ok")

        return FunctionTool(
            name=spec.name,
            description=spec.description,
            params_json_schema=_params_schema(spec),
            on_invoke_tool=invoke,
            strict_json_schema=True,
        )

    sdk_tools = [wrap(spec) for spec in tools]

    async def stop_after_valid_submit(_context, results):
        if "value" in submitted and any(
            item.tool.name == submit_tool for item in results
        ):
            return ToolsToFinalOutputResult(
                is_final_output=True,
                final_output=json.dumps(submitted["value"], ensure_ascii=False),
            )
        return ToolsToFinalOutputResult(is_final_output=False)

    client = AsyncOpenAI(
        api_key=os.environ["ARK_API_KEY"],
        base_url=os.environ.get("DOUBAO_BASE_URL", DOUBAO_BASE_URL_DEFAULT),
    )
    settings = ModelSettings(timeout=wall_seconds) if wall_seconds else ModelSettings()
    sdk_agent = Agent(
        name="eca-pp",
        instructions=system_prompt,
        model=OpenAIResponsesModel(model=model, openai_client=client),
        model_settings=settings,
        tools=sdk_tools,
        tool_use_behavior=stop_after_valid_submit,
    )
    run_config = RunConfig(
        tracing_disabled=True,
        workflow_name=f"eca-pp: {label}",
    )
    _ = effort, allowed_builtin, max_buffer_size
    agent_started = time.perf_counter()
    runtime_init = agent_started - runtime_started
    try:
        run = Runner.run(
            sdk_agent,
            prompt,
            max_turns=max_turns,
            run_config=run_config,
        )
        result = (
            await asyncio.wait_for(run, timeout=wall_seconds)
            if wall_seconds
            else await run
        )
    except (APITimeoutError, asyncio.TimeoutError) as exc:
        message = (
            f"[{label}] agent run exceeded wall-clock limit ({wall_seconds:g}s)"
            if wall_seconds
            else f"[{label}] agent request timed out"
        )
        raise AgentTimeout(message) from exc
    except MaxTurnsExceeded as exc:
        raise AgentIncompleteError(
            f"[{label}] HARNESS=openai exceeded max_turns={max_turns} "
            "without a successful submit call"
        ) from exc
    finally:
        await client.close()
    agent_run = time.perf_counter() - agent_started
    total = time.perf_counter() - runtime_started

    if "value" not in submitted:
        raise AgentIncompleteError(
            f"[{label}] agent finished without a successful {submit_tool} call. "
            f"Final reply:\n{result.final_output}"
        )
    transcript = ItemHelpers.text_message_outputs(result.new_items).strip()
    if not transcript:
        transcript = str(result.final_output)
    print(
        f"== [{label}] HARNESS=openai provider=doubao model={model} "
        f"runtime_init={runtime_init:.3f}s agent_run={agent_run:.3f}s",
        flush=True,
    )
    return AgentRunResult(
        submitted["value"],
        transcript,
        tools_used,
        _usage_dict(result, model, runtime_init, agent_run, total),
    )
