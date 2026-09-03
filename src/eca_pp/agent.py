"""The ONE place eca-pp talks to Claude.

Every LLM use in eca-pp goes through the Claude Agent SDK (which spawns the
``claude`` CLI and uses its stored credentials or ``ANTHROPIC_API_KEY``) —
never the raw messages API, never ``claude -p``. One model for everything
(:data:`DEFAULT_MODEL`), overridable per run via ``--model`` /
``ECA_PP_AGENT_MODEL``. Callers get one self-contained exchange per
:func:`ask`; transient transport failures (CLI cold-start handshake timeouts,
broken pipes) are retried a couple of times, anything else surfaces as
:class:`AgentUnavailable` so the caller can degrade deterministically.
"""

from __future__ import annotations

import logging
import os
import shutil
import time

log = logging.getLogger(__name__)

DEFAULT_MODEL = "claude-sonnet-5"
MODEL_ENV = "ECA_PP_AGENT_MODEL"
CLI_ENV = "ECA_PP_CLAUDE_CLI"
RETRY_DELAYS = (5.0, 20.0)  # seconds before attempt 2 and 3
_TRANSIENT_MARKERS = ("timeout", "timed out", "broken pipe", "connection",
                      "eof", "exited with", "stream closed")


class AgentUnavailable(Exception):
    """No credentials / no SDK / no CLI, or an exchange failed after retries."""


def model_name(model: str | None = None) -> str:
    return model or os.environ.get(MODEL_ENV) or DEFAULT_MODEL


def cli_path() -> str | None:
    """Prefer an external ``claude`` (env override, then PATH): the SDK's
    bundled native binary needs a newer glibc than old cluster OSes ship,
    while the npm-installed JS CLI runs anywhere node runs."""
    return os.environ.get(CLI_ENV) or shutil.which("claude")


def check_available() -> None:
    """Raise :class:`AgentUnavailable` unless the SDK is importable and some
    credential exists (explicit API key or the Claude CLI's stored login)."""
    cli_creds = os.path.expanduser("~/.claude/.credentials.json")
    if not os.environ.get("ANTHROPIC_API_KEY") and not os.path.isfile(cli_creds):
        raise AgentUnavailable("no ANTHROPIC_API_KEY and no Claude CLI credentials")
    try:
        import claude_agent_sdk  # noqa: F401
    except ImportError as exc:
        raise AgentUnavailable(f"claude-agent-sdk not installed: {exc}")


def make_options(*, system_prompt: str, cwd: str, allowed_tools=(), max_turns: int = 1,
                 model: str | None = None):
    """Build ``ClaudeAgentOptions`` with eca-pp's conventions."""
    from claude_agent_sdk import ClaudeAgentOptions
    os.makedirs(cwd, exist_ok=True)  # the SDK needs an existing cwd
    # max_buffer_size: a wide obs profile can exceed the SDK's 1 MiB default
    # decode buffer.
    # Lean session: eca-pp brings its own prompt and needs at most the Read
    # tool, so skip the user's filesystem settings (plugins, hooks, skills)
    # and every configured MCP server. Measured on Sherlock: initialize
    # drops from ~3.4 s to ~0.7 s and no longer touches plugin dirs on NFS
    # or api.anthropic.com/v1/mcp_servers — the handshake cannot be stalled
    # by an unrelated plugin or a network hiccup.
    return ClaudeAgentOptions(
        system_prompt=system_prompt, allowed_tools=list(allowed_tools),
        max_turns=max_turns, cwd=cwd, permission_mode="default",
        cli_path=cli_path(), model=model_name(model),
        setting_sources=[], strict_mcp_config=True,
        max_buffer_size=32 * 1024 * 1024)


def _is_transient(exc: Exception) -> bool:
    try:
        from claude_agent_sdk import CLIConnectionError, CLINotFoundError, ProcessError
        if isinstance(exc, CLINotFoundError):
            return False
        if isinstance(exc, (CLIConnectionError, ProcessError)):
            return True
    except ImportError:
        pass
    if isinstance(exc, (TimeoutError, ConnectionError, BrokenPipeError)):
        return True
    msg = str(exc).lower()
    return any(m in msg for m in _TRANSIENT_MARKERS)


def ask(options, message: str, *, retries: int = len(RETRY_DELAYS)) -> tuple[str, list, dict]:
    """One self-contained session: send ``message``, collect the reply.
    Returns ``(text, tools_used, usage)``; ``usage`` always has the model
    actually used (from the reply stream) and, when the CLI reports them,
    cost and token counts. Retries transient transport failures with
    backoff; raises :class:`AgentUnavailable` when it gives up."""
    import anyio
    from claude_agent_sdk import ClaudeSDKClient

    async def go():
        async with ClaudeSDKClient(options=options) as client:
            await client.query(message)
            chunks, tools = [], []
            usage = {"model": None, "cost_usd": None, "input_tokens": None,
                     "output_tokens": None, "cache_creation_tokens": None,
                     "cache_read_tokens": None, "num_turns": None}
            async for msg in client.receive_response():
                model = getattr(msg, "model", None)  # AssistantMessage
                if model:
                    usage["model"] = model
                for block in getattr(msg, "content", []) or []:
                    text = getattr(block, "text", None)
                    if text:
                        chunks.append(text)
                    name = getattr(block, "name", None)
                    if name:  # ToolUseBlock — record what the agent looked at
                        inp = getattr(block, "input", {}) or {}
                        tools.append({"tool": name, "target": inp.get("file_path", "")})
                if hasattr(msg, "total_cost_usd"):  # final ResultMessage
                    u = getattr(msg, "usage", None) or {}
                    get = (u.get if isinstance(u, dict)
                           else lambda k, d=None: getattr(u, k, d))
                    usage.update({
                        "cost_usd": getattr(msg, "total_cost_usd", None),
                        "input_tokens": get("input_tokens"),
                        "output_tokens": get("output_tokens"),
                        "cache_creation_tokens": get("cache_creation_input_tokens"),
                        "cache_read_tokens": get("cache_read_input_tokens"),
                        "num_turns": getattr(msg, "num_turns", None)})
            return "".join(chunks), tools, usage

    last: Exception | None = None
    for attempt in range(retries + 1):
        try:
            return anyio.run(go)
        except Exception as exc:  # noqa: BLE001
            last = exc
            if attempt >= retries or not _is_transient(exc):
                break
            delay = RETRY_DELAYS[min(attempt, len(RETRY_DELAYS) - 1)]
            log.warning("agent exchange failed (%s: %s) — retry %d/%d in %.0fs",
                        type(exc).__name__, exc, attempt + 1, retries, delay)
            time.sleep(delay)
    raise AgentUnavailable(f"agent exchange failed: {last}")


def extract_json(reply: str) -> str:
    """The text inside the first fenced ```json block of ``reply``."""
    return reply.split("```json", 1)[1].split("```", 1)[0]
