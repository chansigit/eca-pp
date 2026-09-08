"""Build evals/gold/cases.json from finished production runs on Oak.

A case stores the *path* to a run's result.json plus the expected answer; the evidence itself is
rebuilt at eval time by cli.build_evidence(profile, candidates), so a case stays small and always
reflects the current evidence format.

Gold is deliberately tiny — see README. Everything here is either human-corrected or a failure
mode characterised in the issue tracker / git log / memory.
"""

from __future__ import annotations

import json
import sys
from pathlib import Path

OAK = Path("/home/users/chensj16/oak/data/sc")
OUT = Path(__file__).resolve().parent / "gold" / "cases.json"

# mouse-pansci: 12 organs, one obs schema, unambiguous truth. `batch` is a real technical column
# (215 sequencing sub-libraries) that carries no correctable effect; Age_group/Genotype are the
# biology the atlas exists to measure. 4 organs needed a human correction (result.json.orig kept
# the machine verdict), 8 got it right natively — same question, so they are matched controls.
PANSCI_HARD = ["BAT", "gWAT", "iWAT", "liver"]
PANSCI_CONTROL = ["brain", "colon", "duodenum", "heart", "ileum", "jejunum", "muscle", "stomach"]
PANSCI_EXPECT = {
    # what the CLASSIFIER should say (track 1) — not what the pipeline should adopt (track 2)
    "cell_type": "Main_cell_type",
    "column_class": {"batch": "technical", "Age_group": "condition", "Genotype": "condition"},
    "final_batch": None,  # track 2 only: the verdict the pipeline should reach
}

# Failure modes with a named cause. Paths are resolved leniently: a case that cannot be located
# on disk is reported and skipped rather than silently dropped.
CHARACTERISED = [
    # v05-comparison is the post-fix rerun, so it is the one that reflects the intended answer
    {"id": "3ca/Aynaud2020_CellLines",
     "glob": "3ca/eca-pp/v05-comparison-*/Aynaud2020_othermodels_CellLines/identify_columns/result.json",
     "expect": {"column_class": {"cell_cycle_phase": "per_cell_state"}},
     "why": "50cd3da: per-cell state (cell-cycle phase) was ranked as batch because it was the "
            "only probeable column"},
    {"id": "abm-ilcp/ann0608", "glob": "abm-ilcp/eca-pp/identify_columns-v3/result.json",
     "expect": {"cell_type": "ann0608"},
     "why": "1b40a71/38cc972: v0.4 rule prefilter made the real annotation column unselectable"},
    {"id": "tabula-muris-drop/Mammary_Gland",
     "glob": "tabula-muris-drop/eca-pp/*Mammary*/identify_columns/result.json",
     "expect": {"cell_type": "cell_ontology_class"},
     "why": "tracker CAVEATS: agent returned explicit null and cli.py did not fall back; "
            "cell_ontology_class is the obvious column"},
    {"id": "tabula-muris-facs/Kidney",
     "glob": "tabula-muris-facs/eca-pp/*Kidney*/identify_columns/result.json",
     "expect": {"final_batch": None},
     "why": "empty-string artefact: 1-63% of cells have blank tissue/mouse.id, the profiler treats "
            "blank as a real group and 'annotated vs unannotated' looks like a batch"},
]


def pansci_cases():
    out = []
    for organ in PANSCI_HARD + PANSCI_CONTROL:
        p = OAK / "mouse-pansci" / "eca-pp" / organ / "identify_columns" / "result.json"
        if not p.exists():
            print(f"  MISSING {p}", file=sys.stderr)
            continue
        hard = organ in PANSCI_HARD
        out.append({
            "id": f"mouse-pansci/{organ}", "result": str(p),
            "kind": "hard" if hard else "control",
            "expect": PANSCI_EXPECT,
            "why": ("human-corrected: machine adopted a condition column as batch "
                    "(result.json.orig keeps the original verdict)" if hard else
                    "same obs schema as the hard cases; model got it right natively"),
        })
    return out


def characterised_cases():
    out = []
    for spec in CHARACTERISED:
        hits = sorted(OAK.glob(spec["glob"]))
        if not hits:
            print(f"  UNRESOLVED {spec['id']}: no match for {spec['glob']}", file=sys.stderr)
            continue
        out.append({"id": spec["id"], "result": str(hits[0]), "kind": "hard",
                    "expect": spec["expect"], "why": spec["why"]})
        if len(hits) > 1:
            print(f"  note: {spec['id']} matched {len(hits)} paths, took {hits[0].name}",
                  file=sys.stderr)
    return out


def main():
    cases = pansci_cases() + characterised_cases()
    OUT.parent.mkdir(parents=True, exist_ok=True)
    OUT.write_text(json.dumps({"cases": cases}, indent=2) + "\n")
    kinds = {}
    for c in cases:
        kinds[c["kind"]] = kinds.get(c["kind"], 0) + 1
    print(f"wrote {OUT} — {len(cases)} cases {kinds}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
