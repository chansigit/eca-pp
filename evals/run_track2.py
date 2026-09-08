"""Track 2: replay recorded trial metrics through the adoption rule. No LLM, no data, seconds.

Every finished run stores its trials with the metrics AND the class the classifier assigned, so the
adoption decision can be re-derived offline and checked against the gold verdict. This is where
issue #1 lives: `qualifies()` (cli.py) looks only at metrics, and `class` is consulted solely to
emit a `biological_batch_fallback` warning — after which the candidate is adopted anyway.

Also reports what a class-aware rule would do, which is the concrete measurement behind issue #1's
proposed fix ("never adopt a candidate the classifier called `condition`").

  python evals/run_track2.py
"""

from __future__ import annotations

import json
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parent / "src"))


def replay(trials: list[dict], class_aware: bool) -> str | None:
    """The current adoption loop, optionally with issue #1's proposed guard."""
    from eca_pp.identify_columns.cli import correction_unnecessary, qualifies

    for t in trials:
        m = t.get("metrics") or {}
        if not (qualifies(m) or correction_unnecessary(m)):
            continue
        if class_aware and t.get("class") in ("condition", "other"):
            continue  # proposed: a biological grouping is never a batch
        return t.get("batch_col")
    return None


def main() -> int:
    cases = json.loads((HERE / "gold" / "cases.json").read_text())["cases"]
    cases = [c for c in cases if "final_batch" in c["expect"]]
    print(f"{len(cases)} cases carry a final_batch expectation\n")
    print(f"{'case':34} {'gold':10} {'current':12} {'class-aware':12} verdict")
    now_ok = aware_ok = 0
    for c in cases:
        j = json.loads(Path(c["result"]).read_text())
        # .orig preserves the machine verdict on human-corrected runs — that is the honest input
        orig = Path(c["result"] + ".orig")
        trials = json.loads(orig.read_text())["trials"] if orig.exists() else j.get("trials", [])
        gold = c["expect"]["final_batch"]
        cur, aware = replay(trials, False), replay(trials, True)
        now_ok += cur == gold
        aware_ok += aware == gold
        flag = "ok" if cur == gold else ("FIXED by class-aware" if aware == gold else "still wrong")
        print(f"{c['id']:34} {str(gold):10} {str(cur):12} {str(aware):12} {flag}")
    print(f"\ncurrent rule     {now_ok}/{len(cases)}")
    print(f"class-aware rule {aware_ok}/{len(cases)}   <- issue #1's proposed guard")
    return 0


if __name__ == "__main__":
    sys.exit(main())
