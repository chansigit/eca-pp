"""F4a — species resolution ladder (spec §5.2).

    T0  --species CODE          the driver's explicit call; always wins
    T1  stangene.infer_species  deterministic (ID prefixes, mito styles,
                                reference symbol-inventory overlap)
    T2  single-shot LLM         --llm only; validated submit tool; any failure
                                falls through — never retried, never looped
    T3  unresolved              caller blocks with exit 3 and the evidence

The LLM tier is one tool-less harness session via :mod:`eca_pp.agent` (DSH by
default, Claude when selected), whose answer is validated against stangene's
supported species before adoption.
"""

from __future__ import annotations

import logging
import math
from dataclasses import dataclass, field

import stangene

CODE_BY_SPECIES = {
    "human": "hs", "mouse": "mm", "rat": "rn", "zebrafish": "dr",
    "fruit_fly": "dm", "c_elegans": "ce", "cynomolgus": "cyno",
    "rhesus": "rhesus", "marmoset": "marmoset", "mouse_lemur": "lemur",
}

_LLM_SAMPLE = 300  # gene names shown to the LLM
log = logging.getLogger(__name__)


@dataclass
class SpeciesResolution:
    resolved: str | None
    code: str | None
    source: str | None  # cli | inferred | llm | None
    confidence: float
    evidence: dict = field(default_factory=dict)

    def as_dict(self) -> dict:
        return {"resolved": self.resolved, "code": self.code, "source": self.source,
                "confidence": self.confidence, "evidence": self.evidence}


def _gene_ids(adata):
    if "gene_ids" in adata.var.columns:
        return [str(x) for x in adata.var["gene_ids"]]
    return None


def _t1_infer(adata) -> dict:
    return stangene.infer_species(list(adata.var_names), gene_ids=_gene_ids(adata))


def resolve(adata, *, cli_species: str | None = None, llm: bool = False) -> SpeciesResolution:
    """Walk the ladder. ``resolved is None`` means T3: the caller must block.

    Raises ``ValueError`` for an unknown ``--species`` code (a driver mistake,
    not a data problem).
    """
    if cli_species:  # T0
        canon = stangene.resolve_species(cli_species)
        return SpeciesResolution(canon, CODE_BY_SPECIES.get(canon), "cli", 1.0, {})

    t1 = _t1_infer(adata)  # T1
    if t1["species"]:
        return SpeciesResolution(t1["species"], CODE_BY_SPECIES.get(t1["species"]),
                                 "inferred", float(t1["confidence"]), t1["evidence"])

    if llm:  # T2 — explicitly enabled only
        guess = _llm_infer([str(s) for s in adata.var_names[:_LLM_SAMPLE]],
                           t1["evidence"])
        if guess is not None:
            sp, conf = guess[0], guess[1]
            if len(guess) > 2 and guess[2]:  # token usage of the single call
                t1["evidence"]["llm_usage"] = guess[2]
            return SpeciesResolution(sp, CODE_BY_SPECIES.get(sp), "llm", conf,
                                     t1["evidence"])

    return SpeciesResolution(None, None, None, 0.0, t1["evidence"])  # T3


def _llm_infer(symbols_sample: list[str], evidence: dict):
    """One harness session; ``(canonical_species, confidence, usage)`` or
    ``None`` on any failure. Deterministic T3 is always the fallback."""
    import tempfile

    from eca_pp import agent

    supported = ", ".join(sorted(CODE_BY_SPECIES))
    system = (
        "You identify the species of a single-cell RNA-seq dataset from its "
        f"gene identifiers. Answer with one of: {supported}. If the names are "
        "genuinely uninformative, still pick the most likely species but give "
        "a low confidence (< 0.5). Submit exactly this object: "
        '{"species": "<one of the supported names>", "confidence": <0..1>, '
        '"reason": "<one sentence>"}.')
    user = (
        f"Deterministic inference was inconclusive. Its evidence:\n{evidence}\n\n"
        f"A sample of the dataset's feature names:\n{symbols_sample}\n\n"
        "Which species is this dataset from?")
    try:
        agent.check_available()

        def validate(payload: dict) -> dict:
            missing = [key for key in ("species", "confidence", "reason")
                       if key not in payload]
            if missing:
                raise ValueError(f"missing field(s): {missing}")
            payload["species"] = stangene.resolve_species(str(payload["species"]))
            raw_confidence = payload["confidence"]
            if isinstance(raw_confidence, bool):
                raise TypeError("confidence must be a number between 0 and 1")
            confidence = float(raw_confidence)
            if not math.isfinite(confidence) or not 0.0 <= confidence <= 1.0:
                raise ValueError("confidence must be a finite number between 0 and 1")
            payload["confidence"] = confidence
            if not isinstance(payload["reason"], str) or not payload["reason"].strip():
                raise ValueError("reason must be a non-empty sentence")
            return payload

        with tempfile.TemporaryDirectory(prefix="eca-pp-species-") as cwd:
            parsed, _reply, _tools, usage = agent.ask_json(
                system_prompt=system,
                message=user,
                cwd=cwd,
                submit_tool="submit_species",
                schema={
                    "type": "object",
                    "properties": {
                        "species": {"enum": sorted(CODE_BY_SPECIES)},
                        "confidence": {"type": "number", "minimum": 0, "maximum": 1},
                        "reason": {"type": "string"},
                    },
                    "required": ["species", "confidence", "reason"],
                },
                validate=validate,
                allowed_builtin=(),
                max_turns=3,
                label="infer species",
            )
        canon = parsed["species"]
        conf = parsed["confidence"]
        if conf < 0.5:  # the model itself is unsure — let T3 block instead
            log.warning("species LLM unresolved: confidence %.3f below 0.5", conf)
            return None
        usage = {"model": usage.get("model") or agent.model_name(),
                 "cost_usd": usage.get("cost_usd"),
                 "input_tokens": usage.get("input_tokens"),
                 "output_tokens": usage.get("output_tokens"),
                 "backend": usage.get("backend")}
        return canon, conf, usage
    except Exception:  # Any failure falls through to deterministic T3.
        # Preserve diagnostic context without weakening the unresolved fallback.
        log.exception("species LLM failed; falling back to unresolved species")
        return None
