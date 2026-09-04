"""``HARNESS=claude`` backend for :mod:`eca_pp.harness`."""

from __future__ import annotations

import asyncio
import os
import shutil

from .harness import (
    AgentIncompleteError,
    AgentRunResult,
    AgentTimeout,
    AgentUnavailable,
    ToolSpec,
)

MIN_CLAUDE_AGENT_SDK = (0, 2, 152)
CLI_ENV = "ECA_PP_CLAUDE_CLI"
_BUILTIN = {"read": ["Read"], "glob": ["Glob"], "grep": ["Grep"]}


def _version_tuple(value: str) -> tuple[int, ...]:
    result = []
    for part in value.split("."):
        digits = "".join(ch for ch in part if ch.isdigit())
        if not digits:
            break
        result.append(int(digits))
    return tuple(result)


def check_available() -> None:
    from importlib.metadata import PackageNotFoundError, version

    if not os.environ.get("ANTHROPIC_API_KEY") and not os.path.isfile(
        os.path.expanduser("~/.claude/.credentials.json")
    ):
        raise AgentUnavailable("HARNESS=claude has no API key or Claude CLI credentials")
    try:
        installed = version("claude-agent-sdk")
    except PackageNotFoundError:
        raise AgentUnavailable(
            "HARNESS=claude needs claude-agent-sdk>=0.2.152"
        ) from None
    if _version_tuple(installed) < MIN_CLAUDE_AGENT_SDK:
        floor = ".".join(map(str, MIN_CLAUDE_AGENT_SDK))
        raise AgentUnavailable(
            f"claude-agent-sdk {installed} is too old; install >={floor}"
        )


def _cli_path() -> str | None:
    return os.environ.get(CLI_ENV) or shutil.which("claude")


async def _bounded(stream, deadline, label: str, wall_seconds: float | None):
    iterator = stream.__aiter__()
    try:
        while True:
            timeout = None if deadline is None else max(
                0.0, deadline - asyncio.get_running_loop().time()
            )
            try:
                yield await asyncio.wait_for(iterator.__anext__(), timeout=timeout)
            except StopAsyncIteration:
                return
            except asyncio.TimeoutError:
                raise AgentTimeout(
                    f"[{label}] agent run exceeded {wall_seconds / 60:g} min "
                    "(AGENT_WALL_MIN)"
                ) from None
    finally:
        await iterator.aclose()


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
    from claude_agent_sdk import (
        AssistantMessage,
        ClaudeAgentOptions,
        ResultMessage,
        ToolResultBlock,
        ToolUseBlock,
        UserMessage,
        create_sdk_mcp_server,
        query,
        tool,
    )

    os.makedirs(cwd, exist_ok=True)
    server_name = "eca_pp_tools"
    submitted: dict = {}

    def wrap(spec: ToolSpec):
        is_submit = spec.name == submit_tool

        @tool(spec.name, spec.description, spec.input_schema)
        async def handler(args):
            result = await spec.handler(args)
            if is_submit and not result.get("is_error"):
                submitted["value"] = result.get("_submitted", args)
            return {key: value for key, value in result.items() if key != "_submitted"}

        return handler

    server = create_sdk_mcp_server(
        name=server_name, version="1.0.0", tools=[wrap(spec) for spec in tools]
    )
    allowed_tools = [
        name for builtin in allowed_builtin for name in _BUILTIN.get(builtin, [])
    ] + [f"mcp__{server_name}__{spec.name}" for spec in tools]
    options = ClaudeAgentOptions(
        mcp_servers={server_name: server},
        allowed_tools=allowed_tools,
        disallowed_tools=[
            "Bash", "Write", "Edit", "MultiEdit", "NotebookEdit", "WebFetch",
            "WebSearch", "Agent", "Task",
        ],
        permission_mode="bypassPermissions",
        cwd=cwd,
        max_turns=max_turns,
        system_prompt=system_prompt,
        model=model,
        effort=effort,
        cli_path=_cli_path(),
        setting_sources=[],
        strict_mcp_config=True,
        max_buffer_size=max_buffer_size or 32 * 1024 * 1024,
    )

    transcript = None
    tools_used: list[dict] = []
    usage = {
        "backend": "claude", "model": model, "cost_usd": None,
        "input_tokens": None, "output_tokens": None,
        "cache_creation_tokens": None, "cache_read_tokens": None,
        "num_turns": None,
    }
    pending: dict[str, str] = {}
    deadline = None if wall_seconds is None else \
        asyncio.get_running_loop().time() + wall_seconds
    async for message in _bounded(query(prompt=prompt, options=options), deadline, label, wall_seconds):
        if isinstance(message, AssistantMessage):
            if getattr(message, "model", None):
                usage["model"] = message.model
            for block in message.content:
                if isinstance(block, ToolUseBlock):
                    pending[block.id] = block.name
                    target = str(next(iter((block.input or {}).values()), ""))[:200]
                    tools_used.append({"tool": block.name, "target": target})
                    print(f"== [{label}] agent: {block.name}({target[:80]})", flush=True)
        elif isinstance(message, UserMessage) and isinstance(message.content, list):
            for block in message.content:
                if isinstance(block, ToolResultBlock) and block.is_error:
                    detail = block.content if isinstance(block.content, str) else str(block.content)
                    print(
                        f"== [{label}] tool error in {pending.get(block.tool_use_id, '?')}: "
                        f"{detail[:200]!r}", flush=True,
                    )
        elif isinstance(message, ResultMessage):
            transcript = message.result
            if message.is_error or message.subtype != "success":
                raise RuntimeError(
                    f"[{label}] Claude run ended with {message.subtype}: {message.result}"
                )
            raw = message.usage or {}
            get = raw.get if isinstance(raw, dict) else lambda key, default=None: getattr(raw, key, default)
            usage.update({
                "cost_usd": message.total_cost_usd,
                "input_tokens": get("input_tokens"),
                "output_tokens": get("output_tokens"),
                "cache_creation_tokens": get("cache_creation_input_tokens"),
                "cache_read_tokens": get("cache_read_input_tokens"),
                "num_turns": getattr(message, "num_turns", None),
            })

    if "value" not in submitted:
        raise AgentIncompleteError(
            f"[{label}] agent finished without a successful {submit_tool} call. "
            f"Final reply:\n{transcript}"
        )
    return AgentRunResult(submitted["value"], transcript, tools_used, usage)
