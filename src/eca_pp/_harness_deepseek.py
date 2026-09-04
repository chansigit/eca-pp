"""``HARNESS=deepseek`` backend for :mod:`eca_pp.harness`.

DSH can only connect to external MCP servers, while eca-pp's submit handlers
are Python closures. One process-local FastMCP server bridges those handlers
for all sequential DSH sessions. Keeping it alive matters: repeatedly tearing
down its long-lived SSE connection can poison the next FastMCP lifespan and
make DSH's initial tools/list request wait for 60 seconds. The DSH sandbox is
read-only; calls with no allowed builtins also disable its shell/editor tools.
"""

from __future__ import annotations

import asyncio
import atexit
import glob
import inspect
import logging
import os
import shutil
import socket
import tempfile
import threading
import time
import uuid

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


def _keep_session_log(dsh_home: str, cwd: str, label: str) -> None:
    """Keep DSH's transcript when diagnosis is needed after its temp home goes away."""
    safe_label = label.replace(" ", "_").replace("/", "_")
    sources = glob.glob(os.path.join(
        dsh_home, "sessions", "*", "*", "session.jsonl"))
    if not sources:
        return
    source = max(sources, key=os.path.getmtime)
    # Each decision round has the same label; retain every session separately.
    destination = os.path.join(cwd, f"dsh_session_{safe_label}_{uuid.uuid4().hex}.jsonl")
    try:
        shutil.copy(source, destination)
        print(f"== [{label}] dsh session transcript kept at {destination}", flush=True)
    except OSError as exc:
        print(f"== [{label}] could not keep dsh transcript: {exc}", flush=True)


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


_SERVER_INIT_LOCK = threading.Lock()
_CALL_LOCK = threading.Lock()
_MCP_SERVER = None
_MCP_URL: str | None = None
_REGISTERED_TOOLS: set[str] = set()
_ACTIVE_LISTED: threading.Event | None = None
_ACTIVE_HTTP_TRACE: dict[str, list[int]] | None = None
_DSH_RUNTIME = None
_DSH_RUNTIME_KEY: tuple | None = None
_DSH_RUNTIME_HOME: str | None = None
_DSH_RUNTIME_MCP_READY = False


def _ensure_mcp_server():
    """Start one daemonized FastMCP server and reuse it until process exit."""
    global _MCP_SERVER, _MCP_URL
    if _MCP_SERVER is not None:
        return _MCP_SERVER, _MCP_URL
    with _SERVER_INIT_LOCK:
        if _MCP_SERVER is not None:
            return _MCP_SERVER, _MCP_URL

        from mcp.server.fastmcp import FastMCP
        import uvicorn

        port = _free_port()
        mcp_server = FastMCP(
            name="eca-pp", host="127.0.0.1", port=port, stateless_http=True
        )
        original_list_tools = mcp_server._tool_manager.list_tools

        def list_tools_hook(*args, **kwargs):
            if _ACTIVE_LISTED is not None:
                _ACTIVE_LISTED.set()
            return original_list_tools(*args, **kwargs)

        mcp_server._tool_manager.list_tools = list_tools_hook  # type: ignore[method-assign]
        inner_app = mcp_server.streamable_http_app()

        async def app(scope, receive, send):
            trace = _ACTIVE_HTTP_TRACE
            record = None
            if scope["type"] == "http" and trace is not None:
                key = f"{scope['method']} {scope['path']}"
                record = trace.setdefault(key, [0, 0])
                record[0] += 1

            async def traced_send(message):
                if record is not None and message["type"] == "http.response.start":
                    record[1] += 1
                await send(message)

            return await inner_app(scope, receive, traced_send)

        server = uvicorn.Server(uvicorn.Config(
            app, host="127.0.0.1", port=port, log_level="warning", lifespan="on"
        ))

        def serve():
            asyncio.run(server.serve())

        thread = threading.Thread(target=serve, name="eca-pp-fastmcp", daemon=True)
        thread.start()
        for _ in range(200):
            if server.started:
                break
            if not thread.is_alive():
                raise RuntimeError("eca-pp FastMCP server stopped during startup")
            threading.Event().wait(0.05)
        else:
            raise RuntimeError("eca-pp FastMCP server did not start within 10 seconds")

        _MCP_SERVER = mcp_server
        _MCP_URL = f"http://127.0.0.1:{port}{mcp_server.settings.streamable_http_path}"
        return _MCP_SERVER, _MCP_URL


def _close_dsh_runtime() -> None:
    """Close and forget the one reusable DSH subprocess, if any."""
    global _DSH_RUNTIME, _DSH_RUNTIME_KEY, _DSH_RUNTIME_HOME
    global _DSH_RUNTIME_MCP_READY
    runtime, home = _DSH_RUNTIME, _DSH_RUNTIME_HOME
    _DSH_RUNTIME = None
    _DSH_RUNTIME_KEY = None
    _DSH_RUNTIME_HOME = None
    _DSH_RUNTIME_MCP_READY = False
    if runtime is not None:
        try:
            runtime.close()
        except Exception as exc:  # noqa: BLE001  # pragma: no cover - teardown
            logging.getLogger(__name__).debug("DSH teardown failed: %s", exc)
    if home is not None:
        shutil.rmtree(home, ignore_errors=True)


atexit.register(_close_dsh_runtime)


def _ensure_dsh_runtime(
    *, dsh_bin: str, cwd: str, root: str, provider: str, model: str | None,
    effort: str | None, system_prompt: str | None, patch_text: str,
):
    """Return a started runtime, reusing it for an identical route/config."""
    global _DSH_RUNTIME, _DSH_RUNTIME_KEY, _DSH_RUNTIME_HOME
    global _DSH_RUNTIME_MCP_READY
    from deepseek_harness import DeepSeekHarness

    key = (dsh_bin, cwd, provider, model, effort, system_prompt, patch_text)
    if _DSH_RUNTIME is not None and _DSH_RUNTIME_KEY == key:
        return _DSH_RUNTIME, _DSH_RUNTIME_HOME, True, _DSH_RUNTIME_MCP_READY

    _close_dsh_runtime()
    dsh_home = tempfile.mkdtemp(prefix="dsh-home-", dir=root)
    patch_path = os.path.join(dsh_home, "eca-pp.patch.yml")
    try:
        with open(patch_path, "w") as patch_file:
            patch_file.write(patch_text)
        runtime = DeepSeekHarness(
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
        )
        runtime.start()
    except BaseException:
        shutil.rmtree(dsh_home, ignore_errors=True)
        raise
    _DSH_RUNTIME = runtime
    _DSH_RUNTIME_KEY = key
    _DSH_RUNTIME_HOME = dsh_home
    _DSH_RUNTIME_MCP_READY = bool(
        _ACTIVE_LISTED is not None and _ACTIVE_LISTED.is_set()
    )
    return runtime, dsh_home, False, _DSH_RUNTIME_MCP_READY


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
            "failOnStartupError": True,
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


MCP_LIST_GRACE_SECONDS = 90.0


def _run_sync(
    *,
    harness,
    prompt: str,
    session_id: str,
    label: str,
    max_turns: int,
    wall_seconds: float | None,
    listed: threading.Event,
    http_trace: dict,
):
    trace = os.environ.get("DSH_TRACE_EVENTS", "") not in ("", "0")
    turns = 0
    timed_out = threading.Event()
    no_tools = threading.Event()

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

    timers = []
    if wall_seconds is not None:
        def kill():
            timed_out.set()
            print(
                f"== [{label}] wall-clock budget of {wall_seconds / 60:g} min hit; "
                "closing dsh runtime", flush=True,
            )
            harness.close()

        timers.append(threading.Timer(wall_seconds, kill))

    def check_listed():
        if not listed.is_set():
            no_tools.set()
            print(
                f"== [{label}] dsh did not request tools/list within "
                f"{MCP_LIST_GRACE_SECONDS:g}s; closing runtime",
                flush=True,
            )
            harness.close()

    timers.append(threading.Timer(MCP_LIST_GRACE_SECONDS, check_listed))
    for timer in timers:
        timer.daemon = True
        timer.start()

    def stderr_tail(count=25):
        lines = list(getattr(harness.client, "_stderr_lines", []))[-count:]
        return "\n".join(lines)

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
        if no_tools.is_set():
            print(
                f"== [{label}] MCP HTTP trace: {dict(http_trace)}\n"
                f"== [{label}] dsh stderr tail:\n{stderr_tail()}",
                flush=True,
            )
            raise RuntimeError(
                f"[{label}] mcp tools never listed by dsh; mcp-client failed to attach"
            ) from None
        if isinstance(exc, _TurnsExceeded):
            raise AgentIncompleteError(f"{exc} without a successful submit call") from None
        raise
    finally:
        for timer in timers:
            timer.cancel()
    print(
        f"== [{label}] dsh run: {turns} model turn(s); MCP HTTP trace "
        f"(requests, responses started): {dict(http_trace)}",
        flush=True,
    )
    if not listed.is_set():
        print(f"== [{label}] dsh stderr tail:\n{stderr_tail()}", flush=True)
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
    global _ACTIVE_LISTED, _ACTIVE_HTTP_TRACE, _DSH_RUNTIME_MCP_READY
    check_available()
    dsh_bin = _dsh_bin()
    assert dsh_bin is not None
    provider = os.environ.get("DSH_PROVIDER", "doubao")
    submitted: dict = {}
    tools_used: list[dict] = []
    listed = threading.Event()
    http_trace: dict[str, list[int]] = {}
    runtime = None
    dsh_home = None
    runtime_reused = False
    runtime_init_seconds = 0.0
    agent_run_seconds = 0.0
    if not _CALL_LOCK.acquire(blocking=False):
        raise RuntimeError("concurrent DSH calls in one eca-pp process are not supported")
    try:
        mcp_server, mcp_url = await asyncio.to_thread(_ensure_mcp_server)
        assert mcp_url is not None
        for name in list(_REGISTERED_TOOLS):
            mcp_server.remove_tool(name)
        _REGISTERED_TOOLS.clear()
        for spec in tools:
            mcp_server.add_tool(
                _tool_fn(spec, submitted, tools_used, spec.name == submit_tool, label),
                name=spec.name,
                description=spec.description,
            )
            _REGISTERED_TOOLS.add(spec.name)
        _ACTIVE_LISTED = listed
        _ACTIVE_HTTP_TRACE = http_trace

        root = os.environ.get("DSH_HOME_ROOT") or os.environ.get("SCRATCH") or tempfile.gettempdir()
        patch_text = _render_patch(mcp_url, allowed_builtin, provider, model)
        runtime_started = time.perf_counter()
        runtime, dsh_home, runtime_reused, mcp_ready = await asyncio.to_thread(
            _ensure_dsh_runtime,
            dsh_bin=dsh_bin, cwd=cwd, root=root, provider=provider, model=model,
            effort=effort, system_prompt=system_prompt, patch_text=patch_text,
        )
        runtime_init_seconds = time.perf_counter() - runtime_started
        if mcp_ready:
            listed.set()
        print(
            f"== [{label}] HARNESS=deepseek provider={provider} model={model} "
            f"dsh_home={dsh_home} mcp={mcp_url} "
            f"runtime={'reused' if runtime_reused else 'started'}",
            flush=True,
        )
        _ = max_buffer_size
        keep = os.environ.get("DSH_KEEP_SESSIONS", "") not in ("", "0")
        session_id = f"{label.replace(' ', '_')}-{os.getpid()}-{uuid.uuid4().hex}"
        agent_started = time.perf_counter()
        try:
            result, turns = await asyncio.to_thread(
                _run_sync,
                harness=runtime,
                prompt=prompt,
                session_id=session_id,
                label=label,
                max_turns=max_turns,
                wall_seconds=wall_seconds,
                listed=listed,
                http_trace=http_trace,
            )
            agent_run_seconds = time.perf_counter() - agent_started
            _DSH_RUNTIME_MCP_READY = _DSH_RUNTIME_MCP_READY or listed.is_set()
            keep = keep or "value" not in submitted
        except BaseException:
            agent_run_seconds = time.perf_counter() - agent_started
            keep = True
            raise
        finally:
            if keep and dsh_home is not None:
                _keep_session_log(dsh_home, cwd, label)
    except BaseException:
        # A failed turn may poison either the process or its MCP client. Give
        # the existing retry policy a genuinely fresh runtime.
        if runtime is not None:
            await asyncio.to_thread(_close_dsh_runtime)
        raise
    finally:
        _ACTIVE_LISTED = None
        _ACTIVE_HTTP_TRACE = None
        _CALL_LOCK.release()

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
        "timings": {
            "runtime_init": round(runtime_init_seconds, 3),
            "agent_run": round(agent_run_seconds, 3),
            "total": round(runtime_init_seconds + agent_run_seconds, 3),
        },
        "runtime_reused": runtime_reused,
    }
    return AgentRunResult(
        submitted["value"], result.final_response, tools_used, usage
    )
