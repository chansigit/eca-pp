"""Decision policies for identify-columns: who chooses the next candidate and
when to conclude. Everything agent-related lives HERE and nowhere else —
``AgentPolicy`` (with ``ClaudeAgentPolicy`` retained as a compatibility alias)
with ``HeuristicPolicy`` as the
deterministic baseline and test double. Any policy failure raises
:class:`PolicyUnavailable`; the caller degrades deterministically.
"""

from __future__ import annotations

import json


class PolicyUnavailable(Exception):
    def __init__(self, message: str, *, kind: str = "error"):
        super().__init__(message)
        self.kind = kind


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
        active_tier = state.get("active_batch_tier")
        if state.get("probes_left", 0) <= 0:
            return {"action": "give_up", "candidate": None,
                    "reason": "probe budget exhausted without a qualifying batch",
                    "cell_type": state["best_cell_type"]}
        for c in state["candidates"]["batch"]:
            if (not c["excluded"] and not c.get("equivalent_to")
                    and c["label"] not in tried):
                if active_tier and c.get("tier") != active_tier:
                    continue
                return {"action": "probe", "candidate": c["label"],
                        "cell_type": state["best_cell_type"],
                        "reason": f"next viable candidate ({c['class']}, "
                                  f"{c['n_groups']} groups, {c.get('tier', 'primary')} tier)"}
        if not any(not c["excluded"] for c in state["candidates"]["batch"]):
            return {"action": "conclude_no_batch",
                    "reason": "no viable grouping column exists",
                    "cell_type": state["best_cell_type"]}
        return {"action": "give_up", "candidate": None,
                "reason": "all viable candidates probed, none qualified",
                "cell_type": state["best_cell_type"]}


class AgentPolicy:
    """Harness-backed policy (OpenAI by default, DSH/Claude when selected).

    The historical ``ClaudeAgentPolicy`` name remains an alias for source
    compatibility. Construction or any exchange failure raises
    PolicyUnavailable and the caller degrades deterministically.
    """

    def __init__(self, outdir: str, model: str | None = None):
        from eca_pp import agent
        try:
            agent.check_available()
            self._outdir = outdir
            self._model = model
        except agent.AgentUnavailable as exc:
            kind = "timeout" if isinstance(exc.__cause__, agent.AgentTimeout) else "error"
            raise PolicyUnavailable(str(exc), kind=kind) from exc

    def _ask(self, state: dict, message: str) -> tuple[dict, str | None, list, dict]:
        """One decision = one self-contained agent session (the full state is
        resent every round, so no cross-round session memory is needed — and
        every decision stays independently reproducible). Returns the reply
        text, the tool calls the agent made, and the session's token/cost
        usage."""
        from eca_pp import agent

        def validate(decision: dict) -> dict:
            missing = [key for key in ("action", "candidate", "cell_type", "reason")
                       if key not in decision]
            if missing:
                raise ValueError(f"missing field(s): {missing}")
            action = decision["action"]
            if action not in ACTIONS:
                raise ValueError(f"unknown action {action!r}")
            if not isinstance(decision["reason"], str) or not decision["reason"].strip():
                raise ValueError("reason must be a non-empty sentence")
            cell_type = decision["cell_type"]
            valid_cell_types = {
                item["label"] for item in state["candidates"]["cell_type"]
                if item.get("class") == "annotation"
            }
            if cell_type is not None and cell_type not in valid_cell_types:
                raise ValueError(
                    f"cell_type must be an author annotation or null; got {cell_type!r}"
                )
            candidate = decision["candidate"]
            eligible = set(state.get("eligible_batch_candidates", []))
            if action == "probe":
                if state.get("probes_left", 0) <= 0:
                    raise ValueError("probe budget is exhausted")
                if candidate not in eligible:
                    raise ValueError(
                        f"probe candidate must be currently eligible: {sorted(eligible)}"
                    )
            elif action in ("adopt", "conclude_unnecessary"):
                expected = "adopted" if action == "adopt" else "correction_unnecessary"
                accepted = {trial["batch_col"] for trial in state["trials"]
                            if trial["verdict"] == expected}
                if candidate not in accepted:
                    raise ValueError(
                        f"{action} candidate must have verdict {expected!r}: {sorted(accepted)}"
                    )
            elif candidate is not None:
                raise ValueError(f"candidate must be null for {action}")
            if action in ("conclude_no_batch", "give_up") and eligible \
                    and state.get("probes_left", 0) > 0:
                raise ValueError(
                    f"eligible candidates remain and must be probed first: {sorted(eligible)}"
                )
            return decision

        try:
            return agent.ask_json(
                system_prompt=PROMPT,
                message=message,
                cwd=self._outdir,
                submit_tool="submit_column_decision",
                schema=DECISION_SCHEMA,
                validate=validate,
                # The complete decision state is embedded above. Keep the
                # model-facing surface to the validated submit tool only.
                allowed_builtin=(),
                max_turns=6,
                model=self._model,
                label="identify columns",
            )
        except agent.AgentUnavailable as exc:
            kind = "timeout" if isinstance(exc.__cause__, agent.AgentTimeout) else "error"
            raise PolicyUnavailable(str(exc), kind=kind) from exc

    def decide(self, state: dict) -> dict:
        message = (
            "Current state (JSON):\n```json\n"
            + json.dumps(state, ensure_ascii=False, default=str)
            + "\n```\n"
              "Set \"cell_type\" in EVERY reply (probe included): it is the "
              "author annotation to report. candidates.cell_type is ranked; "
              "best_cell_type is the current default. A constant annotation "
              "is valid output but the program will omit it from cLISI.\n"
              "Only propose a batch candidate whose tier equals "
              "active_batch_tier; primary technical/donor candidates must be "
              "exhausted before fallback biological/unknown candidates. "
              "A clear primary probe may be finalized locally by the metric "
              "fast path, so do not promise later probes in your reason; "
              "equivalent groupings need not all be probed.")
        decision, transcript, tools_used, usage = self._ask(state, message)
        decision["tools_used"] = tools_used
        decision["raw_reply"] = transcript or json.dumps(decision, ensure_ascii=False)
        decision["usage"] = usage
        return decision


ACTIONS = ("probe", "adopt", "conclude_unnecessary", "conclude_no_batch",
           "give_up")

DECISION_SCHEMA = {
    "type": "object",
    "properties": {
        "action": {"enum": list(ACTIONS)},
        "candidate": {"type": ["string", "null"]},
        "cell_type": {"type": ["string", "null"]},
        "reason": {"type": "string"},
    },
    "required": ["action", "candidate", "cell_type", "reason"],
}


PROMPT = """\
You identify two roles among the obs columns of a standardized scRNA-seq
dataset: the BATCH column (for integration) and the CELL TYPE column. You
work for an atlas-building pipeline: the purpose of integration is to align
cell identities across experiments.

Evidence provided each round: a three-layer profile (per-column stats with
sampled values and per-value cell counts; group-size health; a nesting/
equivalence graph among grouping columns), candidate lists with pre-check
results, and the metrics of every probe trial so far.

Doctrine:
1. Classify grouping columns first: technical (lane/channel/library/run/
   pool/hash), donor, experimental condition (disease/treatment/timepoint/
   genotype), annotation, QC numeric, identifier.
2. Batch candidates have two strict tiers. PRIMARY = technical and donor/
   sample factors, including technical structure derived from barcodes. Probe
   every viable primary candidate before considering FALLBACK = experimental
   condition or unknown grouping. Never jump to fallback while an untried
   primary candidate remains. A fallback is acceptable only when no primary
   candidate qualifies and its probe preserves biological structure; state
   this biological-risk consequence in the reason.
3. Never batch: annotation columns, QC numeric columns, per-cell-unique
   identifiers, constant columns.
4. Nested groupings: prefer the finest viable technical level (correcting
   at library level already aligns across conditions when libraries nest
   within conditions). If a level is pathological or probes poorly, move
   one level up.
5. Orthogonal groupings: choose ONE column - the one that probes better;
   record the other's existence in your reason.
6. Decide from probe metrics (iLISI gain, cLISI preservation and convergence;
   thresholds are provided).
   If pre-integration iLISI is already high, conclude
   "conclude_unnecessary". If no viable grouping exists at all, conclude
   "conclude_no_batch". If nothing qualifies after the probe budget,
   "give_up" rather than guessing.
   The program may finalize a clear primary trial locally with a conservative
   metric fast path. Do not promise that equivalent candidates will be probed
   later; equivalent groupings need not all be probed.
7. Cell type column = the AUTHOR'S cell-type annotation (biological names
   such as "T cell", "hepatocyte", ontology terms), NOT an algorithmic
   clustering (leiden / louvain / seurat_clusters / numeric cluster IDs).
   candidates.cell_type is pre-ranked (class "annotation" before
   "cluster"; text labels before bare integers) and best_cell_type is
   the current default. Override it when the sampled values show the
   default is wrong (e.g. the "annotation" column holds integers while
   another annotation column holds cell-type names).
   A constant author annotation is still a valid output, although it cannot
   be used for cLISI. If only clusters exist, report null; the probe will make
   its own pseudo-labels. When several annotation columns exist (e.g. coarse and fine), prefer
   the one with recognizable cell-type names at a usable granularity and
   mention the other in your reason. If no author annotation exists, set null;
   that is a successful unattended outcome, not an error.
8. Exhausted or ambiguous batch evidence is also a successful unattended
   outcome: use give_up and let the pipeline report batch=null with warnings.

Finish by submitting this JSON object through the provided submit tool:
{"action": "probe|adopt|conclude_unnecessary|conclude_no_batch|give_up",
 "candidate": "<batch candidate label>",
 "cell_type": "<cell-type column label or null>",
 "reason": "<one sentence naming the evidence>"}
"candidate" is REQUIRED for probe, adopt, and conclude_unnecessary — it names
the batch candidate the action applies to (for conclude_unnecessary: the
identified batch column whose correction is unnecessary). Use null only for
conclude_no_batch / give_up.
"cell_type" is REQUIRED in every reply, probe included: the trials use it
as the cLISI label column, so a wrong choice corrupts the probe metrics.
"""


# Public compatibility for callers written before the harness became
# backend-neutral. New code should use AgentPolicy.
ClaudeAgentPolicy = AgentPolicy
