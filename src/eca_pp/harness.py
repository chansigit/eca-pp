"""Pluggable execution harness for every model-facing call in eca-pp.

Callers expose one validated submit tool and never depend on an agent's
free-text final response.  ``HARNESS=openai`` is the default and drives Doubao
through the OpenAI Agents SDK; ``HARNESS=deepseek`` keeps the DeepSeek Harness
(dsh) path as a fallback; ``HARNESS=claude`` keeps the Claude Agent SDK backend
as an explicit fallback.
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

# Message-based classification is a FALLBACK for exceptions whose only
# information is a string (SDK/CLI errors, DSH stderr). It must never be applied
# to text that may contain model output: an AgentIncompleteError embeds the
# model's final reply, and words like "capacity", "quota", or a number
# containing "429" in that reply used to be mistaken for a provider limit
# (10 min sleeps, up to 12 h). Backends therefore raise the typed
# AgentRateLimited / AgentTransient themselves whenever they can.
LIMIT_PATTERN = re.compile(
    r"usage limit|rate[ _-]?limit|limit will reset|resets at|too many requests|"
    r"overloaded|quota|\b429\b|capacity|out of extra usage|spend limit",
    re.IGNORECASE,
)
TRANSIENT_PATTERN = re.compile(
    r"control request timeout|broken pipe|connection reset|econnreset|epipe|"
    r"process exited unexpectedly|failed to start|connection closed|stdout closed|"
    r"transportclosed|initialize timed out|timed out waiting|\b5(?:00|02|03|04|29)\b|"
    r"connection error|bad gateway|service unavailable|gateway time-?out|"
    r"internal server error|"
    r"mcp tools never listed|initial connection or tool synchronization failed",
    re.IGNORECASE,
)
MAX_TRANSIENT_ATTEMPTS = 5
TRANSIENT_BACKOFF_SECONDS = 20
# A timed-out model request is unlikely to benefit from immediately consuming
# another full wall-clock budget. Callers already degrade to a deterministic
# policy, which keeps unattended dataset processing moving.
MAX_TIMEOUT_ATTEMPTS = 1
# Wall clock per agent run (AGENT_WALL_MIN). 2 min was tuned for ``minimal``
# reasoning; with the ``medium`` default a Doubao decision over a full obs
# profile exceeded it (abm-ilcp, 2026-09-04) and silently degraded the run.
DEFAULT_WALL_MINUTES = 6.0
# Provider limits: wait AGENT_LIMIT_WAIT_MIN minutes per attempt, give up after
# AGENT_LIMIT_WAIT_MAX_H hours. The default ceiling is 1 h — inside a Slurm job
# a longer sleep just burns the allocation; callers fall back deterministically.
DEFAULT_LIMIT_WAIT_MIN = 10.0
DEFAULT_LIMIT_WAIT_MAX_H = 1.0


class AgentError(RuntimeError):
    """Base class for harness failures. Subclasses not listed as retryable
    below are terminal for the current call."""


class AgentUnavailable(AgentError):
    """The selected backend, executable, or credentials are unavailable."""


class AgentIncompleteError(AgentError):
    """The run ended without a successful submit-tool call. Its message may
    contain the model's reply — never pattern-match on it."""


class AgentTimeout(AgentError):
    """The run exceeded its wall-clock budget."""


class AgentLimitExhausted(AgentError):
    """A provider limit outlasted the configured wait budget."""


class AgentTransient(AgentError):
    """Retryable infrastructure failure (connection reset, 5xx, CLI died)."""


class AgentRateLimited(AgentError):
    """Provider usage/rate limit; wait, then retry."""


def classify_error_message(message: str) -> str | None:
    """``"limit"`` / ``"transient"`` / ``None`` for an ERROR string (never a
    model reply). Shared by backends that only get strings from their SDK."""
    if LIMIT_PATTERN.search(message):
        return "limit"
    if TRANSIENT_PATTERN.search(message):
        return "transient"
    return None


def wall_seconds() -> float | None:
    raw = os.environ.get("AGENT_WALL_MIN", "")
    try:
        minutes = float(raw) if raw.strip() else DEFAULT_WALL_MINUTES
    except ValueError:
        return None
    return minutes * 60 if minutes > 0 else None


async def retry_transient(coro_fn: Callable[[], Awaitable[T]], label: str) -> T:
    """Retry transient failures with backoff and wait through bounded rate limits.

    Classification order: typed harness exceptions first (what the backends
    raise when they can see the SDK exception type), then the message-based
    fallback for foreign exceptions only. Every ``AgentError`` subclass other
    than :class:`AgentTransient` / :class:`AgentRateLimited` is terminal —
    in particular :class:`AgentIncompleteError`, whose text is model output.
    """
    wait_min = float(os.environ.get("AGENT_LIMIT_WAIT_MIN", "") or DEFAULT_LIMIT_WAIT_MIN)
    max_h = float(os.environ.get("AGENT_LIMIT_WAIT_MAX_H", "") or DEFAULT_LIMIT_WAIT_MAX_H)
    waited = 0.0
    transient_attempts = 0
    timeout_attempts = 0
    limit_attempt = 0

    async def on_transient(message: str) -> None:
        nonlocal transient_attempts
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

    async def on_limit(message: str) -> None:
        nonlocal limit_attempt, waited
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

    while True:
        try:
            return await coro_fn()
        except AgentTimeout as exc:
            timeout_attempts += 1
            if timeout_attempts >= MAX_TIMEOUT_ATTEMPTS:
                raise
            print(f"== [{label}] {exc} - one fresh attempt", flush=True)
        except AgentTransient as exc:
            await on_transient(str(exc))
        except AgentRateLimited as exc:
            await on_limit(str(exc))
        except AgentError:
            raise  # unavailable / incomplete / exhausted: terminal by design
        except Exception as exc:  # foreign SDK/CLI exception: string fallback
            kind = classify_error_message(str(exc))
            if kind == "transient":
                await on_transient(str(exc))
            elif kind == "limit":
                await on_limit(str(exc))
            else:
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


DEFAULT_BACKEND = "openai"
MODEL_ENV = "ECA_PP_AGENT_MODEL"
_DEFAULT_MODEL = {
    "claude": "claude-sonnet-5",
    "deepseek": "doubao-seed-2-1-turbo-260628",
    "openai": "doubao-seed-2-1-turbo-260628",
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
        elif backend == "openai":
            from ._harness_openai import check_available as check
        else:
            raise AgentUnavailable(
                f"unknown HARNESS backend {backend!r} "
                "(expected 'deepseek', 'openai', or 'claude')"
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
        elif backend == "openai":
            from ._harness_openai import run_agent as run
        else:
            raise AgentUnavailable(
                f"unknown HARNESS backend {backend!r} "
                "(expected 'deepseek', 'openai', or 'claude')"
            )
    except ImportError as exc:
        raise AgentUnavailable(
            f"HARNESS={backend} dependencies are missing: {exc}"
        ) from None

    resolved_model = default_model(model)
    if resolved_model.lower().startswith("claude") and backend != "claude":
        # Batch scripts used to ship --model claude-sonnet-5 without
        # HARNESS=claude; Doubao would reject the name and the whole run
        # silently degraded. Make the mismatch explicit instead.
        raise AgentUnavailable(
            f"model {resolved_model!r} needs HARNESS=claude (current HARNESS="
            f"{backend!r}); export HARNESS=claude or choose a {backend} model"
        )

    wall = wall_seconds()

    async def attempt() -> AgentRunResult:
        coro = run(
            tools=tools,
            submit_tool=submit_tool,
            prompt=prompt,
            system_prompt=system_prompt,
            cwd=cwd,
            model=resolved_model,
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
