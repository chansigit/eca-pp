"""Decision policies for identify-columns: who chooses the next candidate and
when to conclude. Everything agent-related lives HERE and nowhere else —
``ClaudeAgentPolicy`` (Agent SDK) with ``HeuristicPolicy`` as the
deterministic baseline and test double. Any policy failure raises
:class:`PolicyUnavailable`; the caller degrades deterministically.
"""

from __future__ import annotations

import json
import os


class PolicyUnavailable(Exception):
    pass


class HeuristicPolicy:
    """Deterministic bottom-up walker — the audit baseline and test double."""

    def decide(self, state: dict) -> dict:
        tried = {t["batch_col"] for t in state["trials"]}
        for t in state["trials"]:
            if t["verdict"] == "correction_unnecessary":
                return {"action": "conclude_unnecessary", "candidate": t["batch_col"],
                        "reason": "pre-integration iLISI already high",
                        "cell_type": state["best_cell_type"]}
            if t["verdict"] == "adopted":
                return {"action": "adopt", "candidate": t["batch_col"],
                        "reason": "converged with iLISI gain and cLISI preserved",
                        "cell_type": state["best_cell_type"]}
        for c in state["candidates"]["batch"]:
            if not c["excluded"] and c["label"] not in tried:
                return {"action": "probe", "candidate": c["label"],
                        "reason": f"next viable candidate ({c['class']}, "
                                  f"{c['n_groups']} groups, bottom-up order)"}
        if not any(not c["excluded"] for c in state["candidates"]["batch"]):
            return {"action": "conclude_no_batch",
                    "reason": "no viable grouping column exists",
                    "cell_type": state["best_cell_type"]}
        return {"action": "give_up",
                "reason": "all viable candidates probed, none qualified"}


class ClaudeAgentPolicy:
    """Agent SDK-backed policy: the model reads the profile, trial metrics and
    UMAP panels (via its Read tool) and returns one structured decision per
    round. Construction or any exchange failure raises PolicyUnavailable —
    the caller degrades deterministically."""

    def __init__(self, outdir: str, model: str | None = None):
        os.makedirs(outdir, exist_ok=True)  # the SDK needs an existing cwd
        # Auth: an explicit API key, or the Claude Code CLI's stored
        # credentials (the SDK spawns the CLI, which can use either).
        cli_creds = os.path.expanduser("~/.claude/.credentials.json")
        if not os.environ.get("ANTHROPIC_API_KEY") and not os.path.isfile(cli_creds):
            raise PolicyUnavailable(
                "no ANTHROPIC_API_KEY and no Claude CLI credentials")
        try:
            from claude_agent_sdk import ClaudeAgentOptions, ClaudeSDKClient
        except ImportError as exc:
            raise PolicyUnavailable(f"claude-agent-sdk not installed: {exc}")
        import shutil

        # Prefer an external `claude` (env override, then PATH): the SDK's
        # bundled native binary needs a newer glibc than old cluster OSes
        # ship, while the npm-installed JS CLI runs anywhere node runs.
        cli_path = os.environ.get("ECASTEPS_CLAUDE_CLI") or shutil.which("claude")
        self._client_cls = ClaudeSDKClient
        # model=None → the claude CLI's own default; the model actually used
        # is captured from the reply stream into usage["model"] either way.
        # max_buffer_size: the state message (obs profile) for wide datasets
        # can exceed the SDK's 1 MiB default decode buffer.
        self._options = ClaudeAgentOptions(
            system_prompt=PROMPT, allowed_tools=["Read"], max_turns=6,
            cwd=outdir, permission_mode="default", cli_path=cli_path,
            model=model, max_buffer_size=32 * 1024 * 1024)

    def _ask(self, message: str) -> tuple[str, list, dict]:
        """One decision = one self-contained agent session (the full state is
        resent every round, so no cross-round session memory is needed — and
        every decision stays independently reproducible). Returns the reply
        text, the tool calls the agent made (e.g. Reading UMAP panels), and
        the session's token/cost usage."""
        import anyio

        async def go():
            async with self._client_cls(options=self._options) as client:
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
                            tools.append({"tool": name,
                                          "target": inp.get("file_path", "")})
                    if hasattr(msg, "total_cost_usd"):  # final ResultMessage
                        u = getattr(msg, "usage", None) or {}
                        get = (u.get if isinstance(u, dict)
                               else lambda k, d=None: getattr(u, k, d))
                        usage.update({
                            "cost_usd": getattr(msg, "total_cost_usd", None),
                            "input_tokens": get("input_tokens"),
                            "output_tokens": get("output_tokens"),
                            "cache_creation_tokens":
                                get("cache_creation_input_tokens"),
                            "cache_read_tokens": get("cache_read_input_tokens"),
                            "num_turns": getattr(msg, "num_turns", None)})
                return "".join(chunks), tools, usage

        try:
            return anyio.run(go)
        except Exception as exc:  # noqa: BLE001
            raise PolicyUnavailable(f"agent exchange failed: {exc}")

    def decide(self, state: dict) -> dict:
        message = (
            "Current state (JSON):\n```json\n"
            + json.dumps(state, ensure_ascii=False, default=str)
            + "\n```\nTrial UMAP panels are PNG files in the working directory "
              "(paths in trials[].umap); Read them if helpful.\n"
              "Reply with EXACTLY one fenced json block:\n"
              '{"action": "probe|adopt|conclude_unnecessary|conclude_no_batch'
              '|give_up", "candidate": "<label or null>", '
              '"cell_type": "<label or null>", "reason": "<one sentence>"}')
        reply, tools_used, usage = self._ask(message)
        try:
            payload = reply.split("```json", 1)[1].split("```", 1)[0]
            decision = json.loads(payload)
            assert decision.get("action") in (
                "probe", "adopt", "conclude_unnecessary", "conclude_no_batch",
                "give_up")
            decision["tools_used"] = tools_used
            decision["raw_reply"] = reply  # the agent's full reply text
            decision["usage"] = usage      # tokens + cost for this round
            return decision
        except Exception as exc:  # noqa: BLE001
            raise PolicyUnavailable(f"unparseable agent decision: {exc}")


PROMPT = """\
You identify two roles among the obs columns of a standardized scRNA-seq
dataset: the BATCH column (for integration) and the CELL TYPE column. You
work for an atlas-building pipeline: the purpose of integration is to align
cell identities across experiments.

Evidence provided each round: a three-layer profile (per-column stats with
sampled values and per-value cell counts; group-size health; a nesting/
equivalence graph among grouping columns), candidate lists with pre-check
results, and the metrics and UMAP panel of every probe trial so far.

Doctrine:
1. Classify grouping columns first: technical (lane/channel/library/run/
   pool/hash), donor, experimental condition (disease/treatment/timepoint/
   genotype), annotation, QC numeric, identifier.
2. Batch candidates: ANY grouping across which cell identities should be
   aligned - technical, donor, AND experimental condition (atlas setting:
   merging condition-driven expression shifts is the goal; count data stay
   untouched for downstream differential analysis). When you adopt a
   condition-like column, state this consequence in your reason.
3. Never batch: annotation columns, QC numeric columns, per-cell-unique
   identifiers, constant columns.
4. Nested groupings: prefer the finest viable technical level (correcting
   at library level already aligns across conditions when libraries nest
   within conditions). If a level is pathological or probes poorly, move
   one level up.
5. Orthogonal groupings: choose ONE column - the one that probes better;
   record the other's existence in your reason.
6. Decide from probe metrics first (iLISI gain, cLISI preservation,
   convergence; thresholds are provided). Use the UMAP image only as
   supporting evidence or a veto, and state the reason when vetoing.
   If pre-integration iLISI is already high, conclude
   "conclude_unnecessary". If no viable grouping exists at all, conclude
   "conclude_no_batch". If nothing qualifies after the probe budget,
   "give_up" rather than guessing.

Every reply must be EXACTLY one fenced json block:
{"action": "probe|adopt|conclude_unnecessary|conclude_no_batch|give_up",
 "candidate": "<batch candidate label>",
 "cell_type": "<cell-type column label or null>",
 "reason": "<one sentence naming the evidence>"}
"candidate" is REQUIRED for probe, adopt, and conclude_unnecessary — it names
the batch candidate the action applies to (for conclude_unnecessary: the
identified batch column whose correction is unnecessary). Use null only for
conclude_no_batch / give_up.
"""
