"""Validated JSON interface shared by all eca-pp agent call sites.

The backend is selected by ``HARNESS`` through :mod:`eca_pp.harness`.
Answers are accepted only through a caller-named submit tool, so OpenAI,
DSH, and Claude share the validation and retry contract.
"""

from __future__ import annotations

import json
import os
from collections.abc import Callable

from . import harness

AgentUnavailable = harness.AgentUnavailable
AgentTimeout = harness.AgentTimeout
DEFAULT_MODEL = harness._DEFAULT_MODEL[harness.DEFAULT_BACKEND]
MODEL_ENV = harness.MODEL_ENV


def backend_name() -> str:
    return harness.backend_name()


def model_name(model: str | None = None) -> str:
    return harness.default_model(model)


def check_available() -> None:
    harness.check_available()


def ask_json(
    *,
    system_prompt: str,
    message: str,
    cwd: str,
    submit_tool: str,
    schema: dict,
    validate: Callable[[dict], dict | None],
    allowed_builtin: tuple[str, ...] = (),
    max_turns: int = 6,
    model: str | None = None,
    label: str = "agent",
) -> tuple[dict, str | None, list[dict], dict]:
    """Run one session and return a validated JSON object plus audit data.

    ``validate`` may mutate the object and return ``None``, or return a
    replacement object. Raising ``ValueError`` rejects the tool call and
    sends the reason back to the model so it can correct and resubmit within
    the same session.
    """
    import anyio

    os.makedirs(cwd, exist_ok=True)

    async def submit(args: dict) -> dict:
        try:
            payload = json.loads(args["response_json"])
            if not isinstance(payload, dict):
                raise ValueError("response must be a JSON object")
            normalized = validate(payload)
            if normalized is not None:
                payload = normalized
        except (KeyError, json.JSONDecodeError, TypeError, ValueError) as exc:
            return {
                "content": [{
                    "type": "text",
                    "text": f"invalid response, fix and resubmit: {exc}",
                }],
                "is_error": True,
            }
        return {
            "content": [{"type": "text", "text": "response accepted"}],
            "is_error": False,
            "_submitted": payload,
        }

    tool = harness.ToolSpec(
        name=submit_tool,
        description=(
            "Submit the final decision as response_json, a JSON string matching "
            "this schema:\n" + json.dumps(schema, ensure_ascii=False, indent=1)
        ),
        input_schema={"response_json": str},
        handler=submit,
    )
    prompt = (
        message
        + f"\n\nFinish by calling {submit_tool}. Put the complete JSON object in "
          "its response_json argument. Do not finish with only a prose answer."
    )

    async def run():
        return await harness.run_agent(
            tools=[tool],
            submit_tool=submit_tool,
            prompt=prompt,
            system_prompt=system_prompt,
            cwd=cwd,
            model=model_name(model),
            max_turns=max_turns,
            allowed_builtin=allowed_builtin,
            label=label,
            max_buffer_size=32 * 1024 * 1024,
        )

    try:
        result = anyio.run(run)
    except harness.AgentError as exc:
        raise AgentUnavailable(str(exc)) from exc
    except Exception as exc:  # backend SDK exceptions share one caller contract
        raise AgentUnavailable(f"agent exchange failed: {exc}") from exc
    assert result.submitted is not None
    return result.submitted, result.transcript_text, result.tools_used, result.usage
