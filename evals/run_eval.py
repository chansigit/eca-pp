"""Replay gold cases through the classifier at whatever HARNESS/model is configured, and score.

Track 1 only — this scores the model's *classification*, not the pipeline's adoption logic.
See README for why the two must not be merged.

  HARNESS=openai python evals/run_eval.py
  HARNESS=claude ECA_PP_AGENT_MODEL=claude-sonnet-5 python evals/run_eval.py --tag sonnet5
  python evals/run_eval.py --compare
"""

from __future__ import annotations

import argparse
import json
import os
import sys
import tempfile
import time
from pathlib import Path

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parent / "src"))
RUNS = HERE / "runs"


def load_state(result_path: str) -> dict:
    """Rebuild exactly what the classifier saw, from a finished run's result.json.

    Mirrors cli.py's own state construction (see `state = {...}` there) so the replay feeds the
    model the same bytes production did — if that shape ever changes, this must follow.
    """
    from eca_pp.identify_columns.cli import build_evidence, classify_column

    j = json.loads(Path(result_path).read_text())
    profile, candidates = j["profile"], j["candidates"]
    return {"profile": profile, "candidates": candidates,
            "evidence": build_evidence(profile, candidates),
            "heuristic_class": {e["column"]: classify_column(e) for e in profile["columns"]}}


def score(answer: dict, expect: dict) -> dict:
    """Per-case checks. Only the assertions a case actually carries are counted."""
    checks = {}
    if "cell_type" in expect:
        checks["cell_type"] = (answer.get("cell_type") == expect["cell_type"])
    for col, want in (expect.get("column_class") or {}).items():
        got = (answer.get("columns") or {}).get(col)
        if got is None:  # some backends only report the ranked list
            got = next((r.get("class") for r in answer.get("batch_ranked", [])
                        if r.get("column") == col), None)
        checks[f"class:{col}"] = (got == want)
        # a condition column ranked #1 as batch is the issue-#1 shape, worth flagging separately
        if want == "condition":
            ranked = answer.get("batch_ranked") or []
            checks[f"not_top_batch:{col}"] = not (ranked and ranked[0].get("column") == col)
    return checks


def run_case(case: dict, model: str | None) -> dict:
    from eca_pp.identify_columns.policies import AgentClassifier

    t0 = time.time()
    with tempfile.TemporaryDirectory() as td:
        clf = AgentClassifier(outdir=td, model=model)
        answer = clf.classify(load_state(case["result"]))
    return {"id": case["id"], "kind": case["kind"], "elapsed_s": round(time.time() - t0, 1),
            "answer": answer, "checks": score(answer, case["expect"])}


def summarise(records: list[dict]) -> dict:
    passed = sum(sum(r["checks"].values()) for r in records)
    total = sum(len(r["checks"]) for r in records)
    by_kind = {}
    for r in records:
        k = by_kind.setdefault(r["kind"], [0, 0])
        k[0] += sum(r["checks"].values())
        k[1] += len(r["checks"])
    return {"passed": passed, "total": total,
            "rate": round(passed / total, 3) if total else None,
            "by_kind": {k: f"{v[0]}/{v[1]}" for k, v in by_kind.items()}}


def compare() -> int:
    rows = []
    for p in sorted(RUNS.glob("*.json")):
        d = json.loads(p.read_text())
        rows.append((p.stem, d["harness"], d["model"], d["summary"]))
    if not rows:
        print(f"no runs in {RUNS}")
        return 1
    print(f"{'tag':22} {'harness':9} {'model':34} {'score':>9}  by kind")
    for tag, h, m, s in rows:
        print(f"{tag:22} {h:9} {str(m):34} {s['passed']:>4}/{s['total']:<4} {s['by_kind']}")
    return 0


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--tag", help="name for this run (default: harness+model)")
    ap.add_argument("--only", help="substring filter on case id")
    ap.add_argument("--compare", action="store_true", help="table of previous runs, no LLM calls")
    args = ap.parse_args()
    if args.compare:
        return compare()

    from eca_pp.harness import backend_name, default_model

    harness, model = backend_name(), default_model(None)
    cases = json.loads((HERE / "gold" / "cases.json").read_text())["cases"]
    if args.only:
        cases = [c for c in cases if args.only in c["id"]]
    # A case whose only assertion is `final_batch` is about the adoption logic, not the model.
    # Scoring it here would silently award it 0/0; keep it visible and out of the LLM budget.
    track2 = [c for c in cases if not score({}, c["expect"])]
    cases = [c for c in cases if score({}, c["expect"])]
    print(f"HARNESS={harness} model={model} — {len(cases)} track-1 cases"
          + (f", skipping {len(track2)} track-2-only ({', '.join(c['id'] for c in track2)})"
             if track2 else ""), flush=True)

    records = []
    for c in cases:
        try:
            r = run_case(c, os.environ.get("ECA_PP_AGENT_MODEL"))
        except Exception as exc:  # a backend that dies IS a benchmark result
            r = {"id": c["id"], "kind": c["kind"], "error": f"{type(exc).__name__}: {exc}",
                 "checks": {k: False for k in score({}, c["expect"])}}
        records.append(r)
        marks = "".join("." if v else "X" for v in r["checks"].values())
        print(f"  {r['id']:34} {marks:8} {r.get('error', '')}", flush=True)

    summary = summarise(records)
    tag = args.tag or f"{harness}-{model}".replace("/", "_")
    RUNS.mkdir(exist_ok=True)
    (RUNS / f"{tag}.json").write_text(json.dumps(
        {"harness": harness, "model": model, "summary": summary, "records": records}, indent=2))
    print(f"\n{summary['passed']}/{summary['total']} checks passed "
          f"({summary['rate']}) by kind {summary['by_kind']} -> runs/{tag}.json")
    return 0


if __name__ == "__main__":
    sys.exit(main())
