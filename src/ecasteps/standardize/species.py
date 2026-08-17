"""F4a — species resolution ladder (spec §5.2).

    T0  --species CODE          the driver's explicit call; always wins
    T1  stangene.infer_species  deterministic (ID prefixes, mito styles,
                                reference symbol-inventory overlap)
    T2  single-shot LLM         --llm only; structured output; any failure
                                falls through — never retried, never looped
    T3  unresolved              caller blocks with exit 3 and the evidence

The LLM tier is one ``messages.parse`` call (anthropic SDK, NOT an agent loop),
whose answer is validated against stangene's supported species before adoption.
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

DEFAULT_LLM_MODEL = "claude-opus-5"
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
            sp, conf = guess
            return SpeciesResolution(sp, CODE_BY_SPECIES.get(sp), "llm", conf,
                                     t1["evidence"])

    return SpeciesResolution(None, None, None, 0.0, t1["evidence"])  # T3


def _llm_infer(symbols_sample: list[str], evidence: dict):
    """One structured-output LLM call; ``(canonical_species, confidence)`` or
    ``None`` on ANY failure (not installed, no key, API error, unsupported
    answer). Deterministic T3 is the fallback — this tier never retries."""
    try:
        import anthropic
        from pydantic import BaseModel
    except Exception:  # noqa: BLE001 - [llm] extra not installed
        return None

    class SpeciesGuess(BaseModel):
        species: str
        confidence: float
        reason: str

    supported = ", ".join(sorted(CODE_BY_SPECIES))
    system = (
        "You identify the species of a single-cell RNA-seq dataset from its "
        f"gene identifiers. Answer with one of: {supported}. If the names are "
        "genuinely uninformative, still pick the most likely species but give "
        "a low confidence (< 0.5).")
    user = (
        f"Deterministic inference was inconclusive. Its evidence:\n{evidence}\n\n"
        f"A sample of the dataset's feature names:\n{symbols_sample}\n\n"
        "Which species is this dataset from?")
    try:
        client = anthropic.Anthropic()
        resp = client.messages.parse(
            model=os.environ.get("ECASTEPS_LLM_MODEL", DEFAULT_LLM_MODEL),
            max_tokens=1024, system=system,
            messages=[{"role": "user", "content": user}],
            output_format=SpeciesGuess)
        parsed = getattr(resp, "parsed_output", None)
        if parsed is None:
            return None
        canon = stangene.resolve_species(parsed.species)
        conf = max(0.0, min(1.0, float(parsed.confidence)))
        if conf < 0.5:  # the model itself is unsure — let T3 block instead
            return None
        return canon, conf
    except Exception:  # noqa: BLE001 - any failure -> deterministic T3
        return None
