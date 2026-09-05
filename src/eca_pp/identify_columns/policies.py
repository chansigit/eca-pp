"""Column classifiers for identify-columns: ONE call that reads every obs
column's value-count table and names the batch column(s) and the cell-type
column. The host then verifies the batch choice with integration probes.

``AgentClassifier`` asks the configured harness once and validates the
submission; ``HeuristicClassifier`` is the deterministic fallback (name
heuristics) used when no model is available or the call fails. Both return
the same ``classification`` dict::

    {"batch_ranked": [{"column": label, "class": cls, "reason": str}, ...],
     "cell_type": label | None, "cell_type_reason": str,
     "columns": {column: cls, ...}, "notes": str}
"""

from __future__ import annotations

import json

CLASSES = ("technical", "donor", "condition", "annotation", "cluster", "state",
           "qc_numeric", "identifier", "constant", "other")
MAX_BATCH_RANKED = 3


class PolicyUnavailable(Exception):
    def __init__(self, message: str, *, kind: str = "error"):
        super().__init__(message)
        self.kind = kind


class HeuristicClassifier:
    """Deterministic fallback: the host's name-heuristic candidate ranking."""

    def classify(self, state: dict) -> dict:
        ranked = []
        for c in state["candidates"]["batch"]:
            if c["excluded"] or c.get("equivalent_to"):
                continue
            ranked.append({"column": c["label"], "class": c["class"],
                           "reason": f"name heuristic: {c['class']} grouping "
                                     f"with {c['n_groups']} groups"})
            if len(ranked) >= MAX_BATCH_RANKED:
                break
        annotations = [c for c in state["candidates"]["cell_type"]
                       if c["class"] == "annotation"]
        cell_type = annotations[0]["label"] if annotations else None
        return {"batch_ranked": ranked, "cell_type": cell_type,
                "cell_type_reason": ("name heuristic: annotation-named column"
                                     if cell_type else "no annotation-named column"),
                "columns": {c["label"]: c["class"]
                            for c in state["candidates"]["batch"]},
                "notes": "deterministic name heuristics (no model)"}


class AgentClassifier:
    """One harness call over the obs profile; validated submit tool only."""

    def __init__(self, outdir: str, model: str | None = None):
        from eca_pp import agent
        try:
            agent.check_available()
            self._outdir = outdir
            self._model = model
        except agent.AgentUnavailable as exc:
            kind = "timeout" if isinstance(exc.__cause__, agent.AgentTimeout) else "error"
            raise PolicyUnavailable(str(exc), kind=kind) from exc

    def classify(self, state: dict) -> dict:
        from eca_pp import agent

        allowed_batch = {c["label"]: c for c in state["candidates"]["batch"]
                         if not c["excluded"]}
        excluded_batch = {c["label"]: c["note"] for c in state["candidates"]["batch"]
                          if c["excluded"]}
        columns = {e["column"]: e for e in state["profile"]["columns"]}
        heuristic_class = state["heuristic_class"]

        def validate(answer: dict) -> dict:
            missing = [k for k in ("batch_ranked", "cell_type", "cell_type_reason")
                       if k not in answer]
            if missing:
                raise ValueError(f"missing field(s): {missing}")
            ranked = answer["batch_ranked"]
            if not isinstance(ranked, list) or len(ranked) > MAX_BATCH_RANKED:
                raise ValueError(f"batch_ranked must be a list of at most "
                                 f"{MAX_BATCH_RANKED} entries")
            seen = set()
            for item in ranked:
                if not isinstance(item, dict) or "column" not in item:
                    raise ValueError("each batch_ranked entry needs a column")
                col = item["column"]
                if col in seen:
                    raise ValueError(f"batch_ranked repeats {col!r}")
                seen.add(col)
                if col in excluded_batch:
                    raise ValueError(
                        f"{col!r} was excluded before probing: {excluded_batch[col]}")
                if col not in allowed_batch:
                    raise ValueError(
                        f"{col!r} is not a probeable grouping column; choose from "
                        f"{sorted(allowed_batch)}")
                if item.get("class") not in CLASSES:
                    item["class"] = allowed_batch[col]["class"]
                if not str(item.get("reason", "")).strip():
                    raise ValueError(f"give a reason for {col!r}")
            ct = answer["cell_type"]
            if ct is not None:
                entry = columns.get(ct)
                if entry is None:
                    raise ValueError(f"cell_type {ct!r} is not an obs column")
                if entry["is_per_cell_unique"] or entry["n_unique"] < 1:
                    raise ValueError(f"cell_type {ct!r} is a per-cell identifier")
                if heuristic_class.get(ct) == "cluster":
                    raise ValueError(
                        f"{ct!r} holds algorithmic cluster IDs, not the author's "
                        "annotation; pick a column with cell-type names or null")
                if not str(answer["cell_type_reason"]).strip():
                    raise ValueError("cell_type_reason must quote the values")
            cols = answer.get("columns") or {}
            answer["columns"] = {k: v for k, v in cols.items()
                                 if k in columns and v in CLASSES}
            answer.setdefault("notes", "")
            return answer

        message = ("Dataset obs profile (JSON):\n```json\n"
                   + json.dumps(state["evidence"], ensure_ascii=False, default=str)
                   + "\n```")
        try:
            answer, transcript, tools_used, usage = agent.ask_json(
                system_prompt=PROMPT,
                message=message,
                cwd=self._outdir,
                submit_tool="submit_column_classification",
                schema=SCHEMA,
                validate=validate,
                allowed_builtin=(),
                max_turns=6,
                model=self._model,
                label="identify columns",
            )
        except agent.AgentUnavailable as exc:
            kind = "timeout" if isinstance(exc.__cause__, agent.AgentTimeout) else "error"
            raise PolicyUnavailable(str(exc), kind=kind) from exc
        answer["tools_used"] = tools_used
        answer["raw_reply"] = transcript or json.dumps(answer, ensure_ascii=False)
        answer["usage"] = usage
        return answer


SCHEMA = {
    "type": "object",
    "properties": {
        "batch_ranked": {
            "type": "array", "maxItems": MAX_BATCH_RANKED,
            "items": {"type": "object",
                      "properties": {"column": {"type": "string"},
                                     "class": {"enum": list(CLASSES)},
                                     "reason": {"type": "string"}},
                      "required": ["column", "reason"]}},
        "cell_type": {"type": ["string", "null"]},
        "cell_type_reason": {"type": "string"},
        "columns": {"type": "object",
                    "additionalProperties": {"enum": list(CLASSES)}},
        "notes": {"type": "string"},
    },
    "required": ["batch_ranked", "cell_type", "cell_type_reason"],
}


PROMPT = """\
You are given the obs (cell metadata) profile of a standardized scRNA-seq
dataset: for every column its name, dtype, number of distinct values,
missing fraction, and a value -> cell-count table (all values when there are
at most 50, else the 50 most frequent), plus nesting/equivalence relations
between grouping columns and derived candidates (barcode prefix/suffix,
two-column composites). Read the VALUES of every column — names are hints,
values are the truth — and answer two questions in one submission.

1. BATCH column(s), ranked, at most 3. The program will run a small Harmony
   integration trial on each in order and keep the first one that qualifies
   (clear iLISI gain with cell-type structure preserved, or "already mixed").
   - Prefer technical factors (lane/channel/library/run/pool/10x well),
     then donor/sample/animal, then experimental condition. Among nested
     technical levels prefer the finest one whose groups are not mostly tiny.
   - A column nested inside a cell-type-like column, or whose values look like
     "<batch>-<cell type>" (e.g. "ABM2-ILC2P.4"), is batch x annotation:
     use the coarser technical column instead.
   - Never a batch: annotation columns, QC numbers, per-cell identifiers,
     constants, cluster IDs, and per-cell biological STATES such as cell-cycle
     phase, activation/stress state, or doublet/QC bins — correcting on them
     would erase biology. Only columns listed as probeable are allowed, but
     "probeable" only means the program can run a trial on it, not that it
     is a batch: return an EMPTY list rather than a state or annotation
     column when no sample/technical structure exists.
   - An empty list means no plausible batch structure exists in obs.
2. CELL TYPE column: the AUTHOR'S cell-type annotation, judged from values
   (lineage names, ontology terms, abbreviations such as proB, CDP, ILC2P,
   "1:CDP-like"), whatever the column is called (ann0608, ImmGen_refine,
   labels_v2 ...). Never an algorithmic clustering (leiden/louvain/
   seurat_clusters/bare integers). When several exist prefer the one with
   recognizable names at usable granularity and mention the others. null if
   none exists.

Also classify each grouping column (technical/donor/condition/annotation/
cluster/state/qc_numeric/identifier/constant/other) in "columns".

Submit exactly this JSON through the provided tool:
{"batch_ranked": [{"column": "<name>", "class": "<class>", "reason": "<why, citing values>"}],
 "cell_type": "<column or null>", "cell_type_reason": "<quote 2-3 values>",
 "columns": {"<column>": "<class>"}, "notes": "<anything else worth recording>"}
"""
