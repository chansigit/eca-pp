"""Track 2: replay recorded trial metrics through the adoption rule. No LLM, no data, seconds.

Every finished run stores its trials with the metrics AND the class the classifier assigned, so the
per-trial verdict (`cli.trial_verdict`) can be re-derived offline and checked against the gold
verdict. Issue #1 lived here: the metrics alone adopted `condition` columns. Both rules are
reported so the fix stays measurable:

  metrics-only  the pre-#1 rule (class ignored)
  current       trial_verdict() as shipped in the src under test

  python evals/run_track2.py
  ECA_PP_SRC=/path/to/other/checkout/src python evals/run_track2.py   # score another tree
"""

from __future__ import annotations

import json
import os
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
sys.path.insert(0, os.environ.get("ECA_PP_SRC") or str(HERE.parent / "src"))


def replay(trials: list[dict], *, class_aware: bool) -> str | None:
    """The adoption loop over recorded trials: first candidate whose verdict adopts wins."""
    from eca_pp.identify_columns.cli import trial_verdict

    for t in trials:
        cls = t.get("class") if class_aware else None
        if trial_verdict(t.get("metrics") or {}, cls) in ("adopted", "correction_unnecessary"):
            return t.get("batch_col")
    return None


def main() -> int:
    cases = json.loads((HERE / "gold" / "cases.json").read_text())["cases"]
    cases = [c for c in cases if "final_batch" in c["expect"]]
    print(f"{len(cases)} cases carry a final_batch expectation\n")
    print(f"{'case':34} {'gold':10} {'metrics-only':13} {'current':12} verdict")
    base_ok = cur_ok = 0
    for c in cases:
        j = json.loads(Path(c["result"]).read_text())
        # .orig preserves the machine verdict on human-corrected runs — that is the honest input
        orig = Path(c["result"] + ".orig")
        trials = json.loads(orig.read_text())["trials"] if orig.exists() else j.get("trials", [])
        gold = c["expect"]["final_batch"]
        base, cur = replay(trials, class_aware=False), replay(trials, class_aware=True)
        base_ok += base == gold
        cur_ok += cur == gold
        flag = ("ok" if cur == gold else "WRONG") + ("  (class guard changed it)" if base != cur else "")
        print(f"{c['id']:34} {str(gold):10} {str(base):13} {str(cur):12} {flag}")
    print(f"\nmetrics-only rule {base_ok}/{len(cases)}   <- before issue #1")
    print(f"current rule      {cur_ok}/{len(cases)}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
