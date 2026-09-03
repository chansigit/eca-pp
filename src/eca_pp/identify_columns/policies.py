"""Decision policies for identify-columns: who chooses the next candidate and
when to conclude. Everything agent-related lives HERE and nowhere else —
``ClaudeAgentPolicy`` (Agent SDK) with ``HeuristicPolicy`` as the
deterministic baseline and test double. Any policy failure raises
:class:`PolicyUnavailable`; the caller degrades deterministically.
"""

from __future__ import annotations

import json


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
    """Agent SDK-backed policy (via :mod:`eca_pp.agent`): the model reads the
    profile, trial metrics and UMAP panels (via its Read tool) and returns one
    structured decision per round. Construction or any exchange failure
    raises PolicyUnavailable — the caller degrades deterministically."""

    def __init__(self, outdir: str, model: str | None = None):
        from eca_pp import agent
        try:
            agent.check_available()
            self._options = agent.make_options(
                system_prompt=PROMPT, cwd=outdir, allowed_tools=["Read"],
                max_turns=6, model=model)
        except agent.AgentUnavailable as exc:
            raise PolicyUnavailable(str(exc))

    def _ask(self, message: str) -> tuple[str, list, dict]:
        """One decision = one self-contained agent session (the full state is
        resent every round, so no cross-round session memory is needed — and
        every decision stays independently reproducible). Returns the reply
        text, the tool calls the agent made (e.g. Reading UMAP panels), and
        the session's token/cost usage."""
        from eca_pp import agent
        try:
            return agent.ask(self._options, message)
        except agent.AgentUnavailable as exc:
            raise PolicyUnavailable(str(exc))

    def decide(self, state: dict) -> dict:
        message = (
            "Current state (JSON):\n```json\n"
            + json.dumps(state, ensure_ascii=False, default=str)
            + "\n```\nTrial UMAP panels are PNG files in the working directory "
              "(paths in trials[].umap); Read them if helpful.\n"
              "Set \"cell_type\" in EVERY reply (probe included): it is the "
              "cLISI label column for the trials. candidates.cell_type is "
              "ranked; best_cell_type is the current default.\n"
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
7. Cell type column = the AUTHOR'S cell-type annotation (biological names
   such as "T cell", "hepatocyte", ontology terms), NOT an algorithmic
   clustering (leiden / louvain / seurat_clusters / numeric cluster IDs).
   candidates.cell_type is pre-ranked (class "annotation" before
   "cluster"; text labels before bare integers) and best_cell_type is
   the current default. Override it when the sampled values show the
   default is wrong (e.g. the "annotation" column holds integers while
   another column holds cell-type names). Fall back to a cluster column
   ONLY when no annotation column exists; if none is usable, set null.
   When several annotation columns exist (e.g. coarse and fine), prefer
   the one with recognizable cell-type names at a usable granularity and
   mention the other in your reason.

Every reply must be EXACTLY one fenced json block:
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
