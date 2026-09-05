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
from typing import Any

from .harness import (
    AgentIncompleteError,
    AgentRateLimited,
    AgentRunResult,
    AgentTimeout,
    AgentTransient,
    AgentUnavailable,
    ToolSpec,
)

MIN_OPENAI_AGENTS_SDK = (0, 22, 0)
DOUBAO_BASE_URL_DEFAULT = "https://ark.cn-beijing.volces.com/api/v3"
DEFAULT_MAX_NUDGES = 2
# ECAPP calls are short, schema-constrained decisions over evidence already
# computed by the host.  Doubao's ``low`` setting still spent most output
# tokens on hidden reasoning in real runs; ``minimal`` produced the same valid
# tool submission without that overhead.  Callers can raise it for harder work.
DEFAULT_REASONING_EFFORT = "minimal"


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
    totals: dict[str, int],
    model: str | None,
    reasoning_effort: str | None,
    server_state: bool,
    nudges: int,
    runtime_init: float,
    agent_run: float,
    total: float,
) -> dict:
    return {
        "backend": "openai",
        "model": model,
        "reasoning_effort": reasoning_effort,
        "server_state": server_state,
        "parallel_tool_calls": False,
        "nudges": nudges,
        "cost_usd": None,
        "input_tokens": totals["input_tokens"],
        "output_tokens": totals["output_tokens"],
        "reasoning_tokens": totals["reasoning_tokens"],
        "cache_creation_tokens": totals["cache_creation_tokens"],
        "cache_read_tokens": totals["cache_read_tokens"],
        "num_turns": totals["requests"],
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
    from openai import (
        APIConnectionError,
        APIStatusError,
        APITimeoutError,
        AsyncOpenAI,
        RateLimitError,
    )

    os.makedirs(cwd, exist_ok=True)
    max_nudges = int(os.environ.get("OPENAI_AGENTS_MAX_NUDGES", DEFAULT_MAX_NUDGES))
    if max_nudges < 0:
        raise ValueError("OPENAI_AGENTS_MAX_NUDGES must be >= 0")
    server_state = os.environ.get("OPENAI_AGENTS_SERVER_STATE", "1").strip().lower() \
        not in {"0", "false", "no", "off"}
    configured_effort = os.environ.get(
        "OPENAI_AGENTS_REASONING_EFFORT", DEFAULT_REASONING_EFFORT
    ).strip().lower()
    effective_effort = effort or configured_effort
    if effective_effort in {"", "none", "off"}:
        effective_effort = None
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
        max_retries=0,
    )
    settings = ModelSettings(
        timeout=wall_seconds,
        # ECAPP exposes one terminal submit tool. Parallel calls provide no
        # useful concurrency and could race two independently valid decisions.
        parallel_tool_calls=False,
        reasoning={"effort": effective_effort} if effective_effort else None,
        # Ark implements Responses server-side continuation.  Let the SDK send
        # only the incremental turn instead of replaying prior tool results.
        store=server_state,
    )
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
    _ = allowed_builtin, max_buffer_size
    agent_started = time.perf_counter()
    runtime_init = agent_started - runtime_started

    async def run_loop():
        run_input: str | list = prompt
        previous_response_id: str | None = None
        turns_used = 0
        nudges = 0
        transcript_parts: list[str] = []
        usage_totals = {
            "requests": 0,
            "input_tokens": 0,
            "output_tokens": 0,
            "reasoning_tokens": 0,
            "cache_creation_tokens": 0,
            "cache_read_tokens": 0,
        }

        while True:
            remaining = max_turns - turns_used
            if remaining <= 0:
                raise AgentIncompleteError(
                    f"[{label}] HARNESS=openai exhausted max_turns={max_turns} "
                    f"without a successful {submit_tool} call"
                )
            runner_kwargs: dict[str, Any] = {}
            if server_state:
                runner_kwargs["auto_previous_response_id"] = True
                if previous_response_id is not None:
                    runner_kwargs["previous_response_id"] = previous_response_id
            try:
                result = await Runner.run(
                    sdk_agent,
                    run_input,
                    max_turns=remaining,
                    run_config=run_config,
                    **runner_kwargs,
                )
            except MaxTurnsExceeded as exc:
                raise AgentIncompleteError(
                    f"[{label}] HARNESS=openai exceeded max_turns={max_turns} "
                    "without a successful submit call"
                ) from exc

            raw = getattr(getattr(result, "context_wrapper", None), "usage", None)
            input_details = getattr(raw, "input_tokens_details", None)
            output_details = getattr(raw, "output_tokens_details", None)
            requests = int(getattr(raw, "requests", 0) or 0)
            turns_used += max(1, requests)
            usage_totals["requests"] += requests
            usage_totals["input_tokens"] += int(getattr(raw, "input_tokens", 0) or 0)
            usage_totals["output_tokens"] += int(getattr(raw, "output_tokens", 0) or 0)
            usage_totals["reasoning_tokens"] += int(
                getattr(output_details, "reasoning_tokens", 0) or 0
            )
            usage_totals["cache_creation_tokens"] += int(
                getattr(input_details, "cache_write_tokens", 0) or 0
            )
            usage_totals["cache_read_tokens"] += int(
                getattr(input_details, "cached_tokens", 0) or 0
            )
            text = ItemHelpers.text_message_outputs(result.new_items).strip()
            if text:
                transcript_parts.append(text)

            if "value" in submitted:
                return result, transcript_parts, usage_totals, nudges

            if nudges >= max_nudges or turns_used >= max_turns:
                final_text = text or str(result.final_output or "")
                raise AgentIncompleteError(
                    f"[{label}] agent finished without a successful {submit_tool} "
                    f"call after {turns_used} model request(s) and {nudges} "
                    f"nudge(s). Final reply:\n{final_text}"
                )

            nudges += 1
            print(
                f"== [{label}] turn ended without {submit_tool} after "
                f"{turns_used} model request(s) - nudging the same session "
                f"({nudges}/{max_nudges})",
                flush=True,
            )
            nudge_input = {
                "role": "user",
                "content": (
                    f"Your previous turn ended without calling {submit_tool}. "
                    f"Continue exactly where you left off and finish by calling "
                    f"{submit_tool}."
                ),
            }
            if server_state and getattr(result, "last_response_id", None):
                previous_response_id = result.last_response_id
                run_input = [nudge_input]
            else:
                run_input = result.to_input_list()
                run_input.append(nudge_input)

    try:
        outcome = (
            await asyncio.wait_for(run_loop(), timeout=wall_seconds)
            if wall_seconds
            else await run_loop()
        )
    except (APITimeoutError, asyncio.TimeoutError) as exc:
        message = (
            f"[{label}] agent run exceeded wall-clock limit ({wall_seconds:g}s)"
            if wall_seconds
            else f"[{label}] agent request timed out"
        )
        raise AgentTimeout(message) from exc
    except RateLimitError as exc:  # HTTP 429 — typed, so the harness waits
        raise AgentRateLimited(f"[{label}] provider rate limit: {exc}") from exc
    except APIConnectionError as exc:  # "Connection error." — retry with backoff
        raise AgentTransient(f"[{label}] connection failure: {exc}") from exc
    except APIStatusError as exc:
        status = getattr(exc, "status_code", None) or 0
        if status >= 500:  # 500/502/503/504 from Ark — retry with backoff
            raise AgentTransient(f"[{label}] provider {status}: {exc}") from exc
        raise
    finally:
        await client.close()
    result, transcript_parts, usage_totals, nudges = outcome
    agent_run = time.perf_counter() - agent_started
    total = time.perf_counter() - runtime_started

    transcript = "\n\n".join(transcript_parts).strip()
    if not transcript:
        transcript = str(result.final_output)
    print(
        f"== [{label}] HARNESS=openai provider=doubao model={model} "
        f"effort={effective_effort or 'provider-default'} "
        f"server_state={'on' if server_state else 'off'} nudges={nudges} "
        f"runtime_init={runtime_init:.3f}s agent_run={agent_run:.3f}s",
        flush=True,
    )
    return AgentRunResult(
        submitted["value"],
        transcript,
        tools_used,
        _usage_dict(
            usage_totals,
            model,
            effective_effort,
            server_state,
            nudges,
            runtime_init,
            agent_run,
            total,
        ),
    )
