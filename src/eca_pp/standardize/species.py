"""F4a — species resolution ladder (spec §5.2).

    T0  --species CODE          the driver's explicit call; always wins
    T1  stangene.infer_species  deterministic (ID prefixes, mito styles,
                                reference symbol-inventory overlap)
    T2  single-shot LLM         --llm only; structured output; any failure
                                falls through — never retried, never looped
    T3  unresolved              caller blocks with exit 3 and the evidence

The LLM tier is one single-turn, tool-less Agent SDK exchange via
:mod:`eca_pp.agent` (same model and credentials as identify-columns, NOT an
agent loop), whose answer is validated against stangene's supported species
before adoption.
"""

from __future__ import annotations

import os
from dataclasses import dataclass, field

import stangene

CODE_BY_SPECIES = {
    "human": "hs", "mouse": "mm", "rat": "rn", "zebrafish": "dr",
    "fruit_fly": "dm", "c_elegans": "ce", "cynomolgus": "cyno",
    "rhesus": "rhesus", "marmoset": "marmoset", "mouse_lemur": "lemur",
}

_LLM_SAMPLE = 300  # gene names shown to the LLM


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
    """One single-turn Agent SDK exchange; ``(canonical_species, confidence,
    usage)`` or ``None`` on ANY failure (SDK missing, no credentials, exchange
    failed, unparsable or unsupported answer, low confidence). Deterministic
    T3 is the fallback — this tier never loops."""
    import json
    import tempfile

    from eca_pp import agent

    supported = ", ".join(sorted(CODE_BY_SPECIES))
    system = (
        "You identify the species of a single-cell RNA-seq dataset from its "
        f"gene identifiers. Answer with one of: {supported}. If the names are "
        "genuinely uninformative, still pick the most likely species but give "
        "a low confidence (< 0.5). Reply with EXACTLY one fenced json block: "
        '{"species": "<one of the supported names>", "confidence": <0..1>, '
        '"reason": "<one sentence>"} and nothing else.')
    user = (
        f"Deterministic inference was inconclusive. Its evidence:\n{evidence}\n\n"
        f"A sample of the dataset's feature names:\n{symbols_sample}\n\n"
        "Which species is this dataset from?")
    try:
        agent.check_available()
        with tempfile.TemporaryDirectory(prefix="eca-pp-species-") as cwd:
            options = agent.make_options(system_prompt=system, cwd=cwd,
                                         allowed_tools=(), max_turns=1)
            reply, _tools, usage = agent.ask(options, user)
        parsed = json.loads(agent.extract_json(reply))
        canon = stangene.resolve_species(str(parsed["species"]))
        conf = max(0.0, min(1.0, float(parsed["confidence"])))
        if conf < 0.5:  # the model itself is unsure — let T3 block instead
            return None
        usage = {"model": usage.get("model") or agent.model_name(),
                 "cost_usd": usage.get("cost_usd"),
                 "input_tokens": usage.get("input_tokens"),
                 "output_tokens": usage.get("output_tokens"),
                 "billing_url": "https://console.anthropic.com/settings/usage"}
        return canon, conf, usage
    except Exception:  # noqa: BLE001 - any failure -> deterministic T3
        return None
