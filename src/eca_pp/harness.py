"""Pluggable execution harness for every model-facing call in eca-pp.

Callers expose one validated submit tool and never depend on an agent's
free-text final response.  ``HARNESS=deepseek`` uses DeepSeek Harness (dsh)
and is the default; ``HARNESS=claude`` keeps the Claude Agent SDK backend as
an explicit fallback.
"""

from __future__ import annotations

import asyncio
import os
import re
import time
from dataclasses import dataclass, field
from typing import Awaitable, Callable, TypeVar

ToolHandler = Callable[[dict], Awaitable[dict]]
T = TypeVar("T")

LIMIT_PATTERN = re.compile(
    r"usage limit|rate[ _-]?limit|limit will reset|resets at|too many requests|"
    r"overloaded|quota|429|capacity|out of extra usage|spend limit",
    re.IGNORECASE,
)
TRANSIENT_PATTERN = re.compile(
    r"control request timeout|broken pipe|connection reset|econnreset|epipe|"
    r"process exited unexpectedly|failed to start|connection closed|stdout closed|"
    r"transportclosed|initialize timed out|timed out waiting|529|"
    r"mcp tools never listed|initial connection or tool synchronization failed",
    re.IGNORECASE,
)
MAX_TRANSIENT_ATTEMPTS = 5
TRANSIENT_BACKOFF_SECONDS = 20
# A timed-out model request is unlikely to benefit from immediately consuming
# another full wall-clock budget. Callers already degrade to a deterministic
# policy, which keeps unattended dataset processing moving.
MAX_TIMEOUT_ATTEMPTS = 1
DEFAULT_WALL_MINUTES = 2.0


class AgentError(RuntimeError):
    """Base class for harness failures."""


class AgentUnavailable(AgentError):
    """The selected backend, executable, or credentials are unavailable."""


class AgentIncompleteError(AgentError):
    """The run ended without a successful submit-tool call."""


class AgentTimeout(AgentError):
    """The run exceeded its wall-clock budget."""


class AgentLimitExhausted(AgentError):
    """A provider limit outlasted the configured wait budget."""


def wall_seconds() -> float | None:
    raw = os.environ.get("AGENT_WALL_MIN", "")
    try:
        minutes = float(raw) if raw.strip() else DEFAULT_WALL_MINUTES
    except ValueError:
        return None
    return minutes * 60 if minutes > 0 else None


async def retry_transient(coro_fn: Callable[[], Awaitable[T]], label: str) -> T:
    """Retry startup/transient failures and wait through bounded rate limits."""
    wait_min = float(os.environ.get("AGENT_LIMIT_WAIT_MIN", "10"))
    max_h = float(os.environ.get("AGENT_LIMIT_WAIT_MAX_H", "12"))
    waited = 0.0
    transient_attempts = 0
    timeout_attempts = 0
    limit_attempt = 0
    while True:
        try:
            return await coro_fn()
        except AgentUnavailable:
            raise
        except AgentTimeout as exc:
            timeout_attempts += 1
            if timeout_attempts >= MAX_TIMEOUT_ATTEMPTS:
                raise
            print(f"== [{label}] {exc} - one fresh attempt", flush=True)
        except Exception as exc:
            message = str(exc)
            if TRANSIENT_PATTERN.search(message):
                transient_attempts += 1
                if transient_attempts >= MAX_TRANSIENT_ATTEMPTS:
                    raise AgentError(
                        f"[{label}] transient agent failure persisted after "
                        f"{transient_attempts} attempts: {message}"
                    ) from None
                delay = TRANSIENT_BACKOFF_SECONDS * transient_attempts
                print(
                    f"== [{label}] transient agent failure "
                    f"({transient_attempts}/{MAX_TRANSIENT_ATTEMPTS}): "
                    f"{message[:160]!r} - retrying in {delay}s",
                    flush=True,
                )
                await asyncio.sleep(delay)
                continue
            if LIMIT_PATTERN.search(message):
                limit_attempt += 1
                if waited / 3600 >= max_h:
                    raise AgentLimitExhausted(
                        f"[{label}] usage limit remained after {waited / 3600:.1f} h: {message}"
                    ) from None
                print(
                    f"== [{label}] usage/rate limit (attempt {limit_attempt}): "
                    f"{message[:160]!r} - waiting {wait_min:g} min",
                    flush=True,
                )
                started = time.time()
                await asyncio.sleep(wait_min * 60)
                waited += time.time() - started
                continue
            raise


@dataclass(frozen=True)
class ToolSpec:
    name: str
    description: str
    input_schema: dict[str, type]
    handler: ToolHandler


@dataclass
class AgentRunResult:
    submitted: dict | None
    transcript_text: str | None
    tools_used: list[dict] = field(default_factory=list)
    usage: dict = field(default_factory=dict)


DEFAULT_BACKEND = "deepseek"
MODEL_ENV = "ECA_PP_AGENT_MODEL"
_DEFAULT_MODEL = {
    "claude": "claude-sonnet-5",
    "deepseek": "doubao-seed-2-1-turbo-260628",
}


def backend_name() -> str:
    return os.environ.get("HARNESS", DEFAULT_BACKEND).strip().lower()


def default_model(model: str | None = None) -> str:
    backend = backend_name()
    return model or os.environ.get(MODEL_ENV) or os.environ.get("MODEL") or \
        _DEFAULT_MODEL.get(backend, _DEFAULT_MODEL["claude"])


def check_available() -> None:
    backend = backend_name()
    try:
        if backend == "claude":
            from ._harness_claude import check_available as check
        elif backend == "deepseek":
            from ._harness_deepseek import check_available as check
        else:
            raise AgentUnavailable(
                f"unknown HARNESS backend {backend!r} (expected 'claude' or 'deepseek')"
            )
    except ImportError as exc:
        raise AgentUnavailable(
            f"HARNESS={backend} dependencies are missing: {exc}"
        ) from None
    check()


async def run_agent(
    *,
    tools: list[ToolSpec],
    submit_tool: str,
    prompt: str,
    system_prompt: str | None,
    cwd: str,
    model: str | None = None,
    effort: str | None = None,
    max_turns: int = 30,
    allowed_builtin: tuple[str, ...] = ("read", "glob", "grep"),
    label: str = "agent",
    max_buffer_size: int | None = None,
) -> AgentRunResult:
    """Run one bounded session and require ``submit_tool`` to succeed."""
    backend = backend_name()
    try:
        if backend == "claude":
            from ._harness_claude import run_agent as run
        elif backend == "deepseek":
            from ._harness_deepseek import run_agent as run
        else:
            raise AgentUnavailable(
                f"unknown HARNESS backend {backend!r} (expected 'claude' or 'deepseek')"
            )
    except ImportError as exc:
        raise AgentUnavailable(
            f"HARNESS={backend} dependencies are missing: {exc}"
        ) from None

    wall = wall_seconds()

    async def attempt() -> AgentRunResult:
        coro = run(
            tools=tools,
            submit_tool=submit_tool,
            prompt=prompt,
            system_prompt=system_prompt,
            cwd=cwd,
            model=default_model(model),
            effort=effort,
            max_turns=max_turns,
            allowed_builtin=allowed_builtin,
            label=label,
            max_buffer_size=max_buffer_size,
            wall_seconds=wall,
        )
        if wall is None:
            return await coro
        try:
            return await asyncio.wait_for(coro, timeout=wall + 120)
        except asyncio.TimeoutError:
            raise AgentTimeout(
                f"[{label}] agent run exceeded {wall / 60:g} min (AGENT_WALL_MIN)"
            ) from None

    return await retry_transient(attempt, label)
