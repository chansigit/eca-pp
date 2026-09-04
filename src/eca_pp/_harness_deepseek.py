"""``HARNESS=deepseek`` backend for :mod:`eca_pp.harness`.

DSH can only connect to external MCP servers, while eca-pp's submit handlers
are Python closures.  A short-lived localhost FastMCP server bridges those
handlers into one DSH session.  The DSH sandbox is read-only; tool-less calls
also disable its shell/editor tools.
"""

from __future__ import annotations

import asyncio
import inspect
import logging
import os
import socket
import tempfile
import threading

import yaml

from .harness import (
    AgentIncompleteError,
    AgentRunResult,
    AgentTimeout,
    AgentUnavailable,
    ToolSpec,
)

for _logger_name in (
    "mcp", "mcp.server", "mcp.server.streamable_http",
    "mcp.server.streamable_http_manager", "mcp.server.lowlevel.server",
):
    logging.getLogger(_logger_name).setLevel(logging.WARNING)

_BUILTIN_DISABLE_IDS = (
    "persistent-bash", "terminal-bash", "persistent-pwsh", "terminal-pwsh",
    "str-replace-editor",
)
DOUBAO_BASE_URL_DEFAULT = "https://ark.cn-beijing.volces.com/api/v3"


def _default_dsh_bin() -> str | None:
    scratch = os.environ.get("SCRATCH")
    if scratch:
        candidate = os.path.join(
            scratch, "tools", "deepseek-harness-src", "apps", "cli", "lib", "bin.js"
        )
        if os.path.isfile(candidate):
            return candidate
    return None


def _dsh_bin() -> str | None:
    return os.environ.get("DSH_BIN") or _default_dsh_bin()


def check_available() -> None:
    dsh_bin = _dsh_bin()
    if not dsh_bin or not os.path.isfile(dsh_bin):
        raise AgentUnavailable(
            "HARNESS=deepseek needs DSH_BIN pointing at a built dsh CLI "
            "entrypoint (apps/cli/lib/bin.js)"
        )
    try:
        import deepseek_harness  # noqa: F401
        import mcp  # noqa: F401
        import uvicorn  # noqa: F401
    except ImportError as exc:
        raise AgentUnavailable(
            f"HARNESS=deepseek dependencies are missing: {exc}"
        ) from None
    provider = os.environ.get("DSH_PROVIDER", "doubao")
    credential = "DEEPSEEK_API_KEY" if provider == "deepseek-official" else "ARK_API_KEY"
    if not os.environ.get(credential):
        raise AgentUnavailable(
            f"HARNESS=deepseek provider {provider!r} needs {credential}"
        )


def _free_port() -> int:
    with socket.socket(socket.AF_INET, socket.SOCK_STREAM) as sock:
        sock.bind(("127.0.0.1", 0))
        return sock.getsockname()[1]


def _tool_fn(
    spec: ToolSpec,
    submitted: dict,
    tools_used: list[dict],
    is_submit: bool,
    label: str,
):
    import mcp.types as types

    parameters = [
        inspect.Parameter(name, inspect.Parameter.POSITIONAL_OR_KEYWORD, annotation=kind)
        for name, kind in spec.input_schema.items()
    ]

    async def function(**kwargs):
        target = str(next(iter(kwargs.values()), ""))[:200]
        tools_used.append({"tool": spec.name, "target": target})
        print(f"== [{label}] agent: {spec.name}({target[:80]})", flush=True)
        result = await spec.handler(kwargs)
        if result.get("is_error"):
            detail = " ".join(str(item.get("text", "")) for item in result.get("content", []))
            print(f"== [{label}] tool error in {spec.name}: {detail[:200]!r}", flush=True)
        if is_submit and not result.get("is_error"):
            submitted["value"] = result.get("_submitted", kwargs)
        content: list[types.ContentBlock] = [
            types.TextContent(**item) for item in result["content"]
        ]
        return types.CallToolResult(
            content=content, isError=bool(result.get("is_error", False))
        )

    function.__signature__ = inspect.Signature(parameters)  # type: ignore[attr-defined]
    function.__name__ = spec.name
    return function


def _render_patch(
    mcp_url: str, allowed_builtin: tuple[str, ...], provider: str, model: str | None
) -> str:
    insert: list[dict] = [{
        "id": "eca-pp-tools",
        "name": "@deepseek-ai/dsh-mcp-client",
        "config": {
            "serverName": "eca-pp", "transport": "streamable-http", "url": mcp_url,
        },
    }]
    rows: list[dict] = [{"id": "sandbox-policy", "config": {"mode": "read-only"}}]
    if provider != "deepseek-official":
        if not model:
            raise ValueError(
                f"HARNESS=deepseek with DSH_PROVIDER={provider!r} needs a model id"
            )
        insert.append({
            "id": "eca-pp-llm-provider",
            "name": "@deepseek-ai/dsh-llm-pi-ai",
            "config": {"providers": {provider: {
                "apiKeyEnv": "ARK_API_KEY",
                "api": "openai-completions",
                "baseURL": os.environ.get("DOUBAO_BASE_URL", DOUBAO_BASE_URL_DEFAULT),
                "models": [{"id": model}],
            }}},
        })
    if not allowed_builtin:
        rows.extend({"id": plugin_id, "disabled": True} for plugin_id in _BUILTIN_DISABLE_IDS)
    return yaml.safe_dump([{"insert": insert}, *rows], sort_keys=False)


class _TurnsExceeded(RuntimeError):
    pass


def _run_sync(
    *,
    dsh_bin: str,
    cwd: str,
    dsh_home: str,
    provider: str,
    model: str | None,
    effort: str | None,
    system_prompt: str | None,
    prompt: str,
    session_id: str,
    patch_path: str,
    label: str,
    max_turns: int,
    wall_seconds: float | None,
):
    from deepseek_harness import DeepSeekHarness

    trace = os.environ.get("DSH_TRACE_EVENTS", "") not in ("", "0")
    turns = 0
    timed_out = threading.Event()

    def on_notification(notification):
        nonlocal turns
        if notification.method != "session.event":
            return
        event = notification.payload.get("event") or {}
        kind = event.get("type")
        if trace:
            print(f"== [{label}] dsh event: {kind}", flush=True)
        if kind == "assistant/message":
            turns += 1
            if turns > max_turns:
                raise _TurnsExceeded(
                    f"[{label}] HARNESS=deepseek exceeded max_turns={max_turns}"
                )

    with DeepSeekHarness(
        provider=provider,
        model=model,
        reasoning_effort=effort,
        cwd=cwd,
        dsh_home=dsh_home,
        profile="sdk-minimal",
        patches=(patch_path,),
        dsh_bin=dsh_bin,
        env={"DSH_SYSTEM_PROMPT": system_prompt} if system_prompt else {},
        initialize_timeout_seconds=90.0,
    ) as harness:
        watchdog = None
        if wall_seconds is not None:
            def kill():
                timed_out.set()
                print(
                    f"== [{label}] wall-clock budget of {wall_seconds / 60:g} min hit; "
                    "closing dsh runtime", flush=True,
                )
                harness.close()

            watchdog = threading.Timer(wall_seconds, kill)
            watchdog.daemon = True
            watchdog.start()
        try:
            result = harness.run(
                prompt, session_id=session_id, on_notification=on_notification
            )
        except Exception as exc:
            if timed_out.is_set():
                raise AgentTimeout(
                    f"[{label}] agent run exceeded {wall_seconds / 60:g} min "
                    "(AGENT_WALL_MIN)"
                ) from None
            if isinstance(exc, _TurnsExceeded):
                raise AgentIncompleteError(f"{exc} without a successful submit call") from None
            raise
        finally:
            if watchdog is not None:
                watchdog.cancel()
        print(f"== [{label}] dsh run: {turns} model turn(s)", flush=True)
        return result, turns


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
    from mcp.server.fastmcp import FastMCP
    import uvicorn

    check_available()
    dsh_bin = _dsh_bin()
    assert dsh_bin is not None
    provider = os.environ.get("DSH_PROVIDER", "doubao")
    submitted: dict = {}
    tools_used: list[dict] = []
    port = _free_port()
    mcp_server = FastMCP(
        name=f"eca-pp-{label}", host="127.0.0.1", port=port, stateless_http=True
    )
    for spec in tools:
        mcp_server.add_tool(
            _tool_fn(spec, submitted, tools_used, spec.name == submit_tool, label),
            name=spec.name,
            description=spec.description,
        )
    mcp_url = f"http://127.0.0.1:{port}{mcp_server.settings.streamable_http_path}"

    server = uvicorn.Server(uvicorn.Config(
        mcp_server.streamable_http_app(), host="127.0.0.1", port=port,
        log_level="warning", lifespan="on",
    ))
    server_task = asyncio.create_task(server.serve())
    try:
        for _ in range(50):
            await asyncio.sleep(0.05)
            with socket.socket(socket.AF_INET, socket.SOCK_STREAM) as probe:
                if probe.connect_ex(("127.0.0.1", port)) == 0:
                    break

        root = os.environ.get("DSH_HOME_ROOT") or os.environ.get("SCRATCH") or tempfile.gettempdir()
        with tempfile.TemporaryDirectory(prefix="dsh-home-", dir=root) as dsh_home:
            patch_text = _render_patch(mcp_url, allowed_builtin, provider, model)
            with tempfile.NamedTemporaryFile(
                "w", suffix=".patch.yml", dir=dsh_home, delete=False
            ) as patch_file:
                patch_file.write(patch_text)
                patch_path = patch_file.name
            print(
                f"== [{label}] HARNESS=deepseek provider={provider} model={model} "
                f"dsh_home={dsh_home} mcp={mcp_url}", flush=True,
            )
            _ = max_buffer_size
            result, turns = await asyncio.to_thread(
                _run_sync,
                dsh_bin=dsh_bin,
                cwd=cwd,
                dsh_home=dsh_home,
                provider=provider,
                model=model,
                effort=effort,
                system_prompt=system_prompt,
                prompt=prompt,
                session_id=f"{label}-{os.getpid()}",
                patch_path=patch_path,
                label=label,
                max_turns=max_turns,
                wall_seconds=wall_seconds,
            )
    finally:
        server.should_exit = True
        try:
            await asyncio.wait_for(server_task, timeout=10)
        except (asyncio.TimeoutError, asyncio.CancelledError, Exception):
            server_task.cancel()

    if result.finish_reason == "error":
        errors = [
            event for event in result.events
            if event.get("type") == "turn/end"
            and event.get("data", {}).get("reason", {}).get("kind") == "error"
        ]
        detail = errors[-1]["data"]["reason"]["error"] if errors else result.final_response
        raise RuntimeError(f"[{label}] HARNESS=deepseek run ended in error: {detail}")
    if "value" not in submitted:
        raise AgentIncompleteError(
            f"[{label}] agent finished without a successful {submit_tool} call. "
            f"Final reply:\n{result.final_response}"
        )
    usage = {
        "backend": "deepseek", "model": model, "cost_usd": None,
        "input_tokens": None, "output_tokens": None,
        "cache_creation_tokens": None, "cache_read_tokens": None,
        "num_turns": turns,
    }
    return AgentRunResult(
        submitted["value"], result.final_response, tools_used, usage
    )
