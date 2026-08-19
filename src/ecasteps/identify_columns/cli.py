"""identify-columns — the project's first declared agent-implemented step
(identify-columns spec). To its caller it is a plain CLI (h5ad in →
result.json out, exit codes 0/3/2/1); internally an agent chooses which batch
candidate to probe next and when to conclude, while every tool invocation —
profiling and integration-probe trials — is a deterministic CLI, recorded
round by round in ``trials``.

Degraded mode (no Agent SDK / no API key / --no-probe): the deterministic
profile + heuristic ranking are still produced, status=needs_review, exit 3.
"""

from __future__ import annotations

import argparse
import json
import logging
import os
import shutil
import sys
import time

from ecasteps.identify_columns import obsprofile
from ecasteps.identify_columns.policies import (
    ClaudeAgentPolicy,
    HeuristicPolicy,  # noqa: F401 - part of the public policy surface
    PolicyUnavailable,
)
from ecasteps.probe import cli as probe
from ecasteps.core.colspec import write_values_tsv
from ecasteps.core.result import EXIT_BLOCKED, EXIT_ERROR, EXIT_OK, \
    new_result, write_result

log = logging.getLogger("ecasteps.identify_columns")

# --- decision thresholds (v0.3 defaults; all recorded in result.json) --------
PATHOLOGICAL_TINY_GROUP_FRAC = 0.5   # majority of groups tiny -> column is out
TINY_NOTE_CELL_FRAC = 0.05           # below this, tiny groups only get a note
PRE_MIXED_ILISI = 0.8                # pre-iLISI at/above -> correction unnecessary
PRE_MIXED_PCR = 0.05                 # ...together with PC-regression R2 below
ILISI_GAIN_MIN = 0.05                # normalized iLISI gain required to adopt
CLISI_DROP_TOL = 0.05                # tolerated normalized cLISI drop
CELLS_PER_BATCH = 50                 # adaptive sampling: expected cells/batch
N_CELLS_FLOOR, N_CELLS_CAP = 5000, 30000
MAX_PROBES = 6

TECH_TOKENS = ("lane", "channel", "library", "batch", "run", "pool", "hash",
               "chip", "well", "plate", "flowcell", "kit", "10x", "lib",
               "gem", "seq")
DONOR_TOKENS = ("donor", "mouse", "patient", "subject", "individual", "animal",
                "sample", "specimen", "rep")
CONDITION_TOKENS = ("condition", "disease", "treatment", "treat", "stim",
                    "genotype", "timepoint", "time", "stage", "diet", "dose",
                    "age", "sex", "group")
ANNOTATION_TOKENS = ("celltype", "cell_type", "annotation", "cluster",
                     "louvain", "leiden", "ontology", "class", "lineage",
                     "subtype", "cell_label")


def build_parser() -> argparse.ArgumentParser:
    p = argparse.ArgumentParser(
        prog="ecasteps-identify-columns",
        description="Identify the batch column and cell-type column of a "
                    "standardized h5ad, verified by small-scale integration "
                    "trials; writes OUTDIR/result.json (+ batch.tsv for a "
                    "derived batch column).")
    p.add_argument("src", help="standardized .h5ad")
    p.add_argument("-o", "--outdir", required=True)
    p.add_argument("--max-probes", type=int, default=MAX_PROBES)
    p.add_argument("--n-cells", type=int, default=None,
                   help="probe subsample size (default: adaptive, "
                        "clamp(50×max_batches, 5000, 30000))")
    p.add_argument("--no-probe", action="store_true",
                   help="profile + ranking only (degraded mode, exit 3)")
    p.add_argument("--seed", type=int, default=0)
    p.add_argument("--model", default=os.environ.get("ECASTEPS_AGENT_MODEL"),
                   help="agent model ID (e.g. claude-sonnet-5); default: "
                        "$ECASTEPS_AGENT_MODEL, else the claude CLI's "
                        "default model")
    return p


# ------------------------------------------------------------ classification

def _norm(name: str) -> str:
    return name.lower().replace("_", "").replace(".", "").replace(" ", "")


def classify_column(entry: dict) -> str:
    """technical | donor | condition | annotation | qc_numeric | identifier |
    constant | other — doctrine §4.1; structural classes win over names."""
    if entry["is_constant"]:
        return "constant"
    if entry["is_per_cell_unique"]:
        return "identifier"
    if entry["dtype"] == "float":
        return "qc_numeric"
    n = _norm(entry["column"])
    for tokens, label in ((ANNOTATION_TOKENS, "annotation"),
                          (TECH_TOKENS, "technical"),
                          (DONOR_TOKENS, "donor"),
                          (CONDITION_TOKENS, "condition")):
        if any(t in n for t in tokens):
            return label
    return "other"


def _pathology(entry: dict) -> str | None:
    gs = entry.get("group_sizes")
    if gs is None:
        return "not a grouping column"
    if gs["tiny_group_frac"] > PATHOLOGICAL_TINY_GROUP_FRAC:
        return (f"pathological: {gs['n_tiny']}/{gs['n_groups']} groups are "
                f"tiny (<{obsprofile.TINY_GROUP_CELLS} cells)")
    return None


def build_candidates(profile: dict) -> dict:
    """{'batch': [...], 'cell_type': [...]} — classified, pre-checked,
    bottom-up ordered (finest viable technical level first)."""
    batch, cell_type = [], []
    class_of = {e["column"]: classify_column(e) for e in profile["columns"]}
    for e in profile["columns"]:
        cls = class_of[e["column"]]
        if cls == "annotation" and obsprofile.is_grouping_candidate(e):
            cell_type.append({"label": e["column"], "kind": "existing",
                              "class": cls, "n_groups": e["n_unique"]})
        if cls in ("technical", "donor", "condition", "other") \
                and obsprofile.is_grouping_candidate(e):
            note = _pathology(e)
            gs = e["group_sizes"]
            cand = {"label": e["column"], "kind": "existing", "class": cls,
                    "n_groups": gs["n_groups"],
                    "tiny_cell_frac": gs["tiny_cell_frac"],
                    "excluded": bool(note), "note": note or ""}
            if not note and 0 < gs["tiny_cell_frac"] <= TINY_NOTE_CELL_FRAC:
                cand["note"] = (f"{gs['n_tiny']} tiny groups "
                                f"({gs['tiny_cell_frac']:.1%} of cells)")
            batch.append(cand)
    for d in profile["derived"]:
        # A composite built on an annotation column would fold biology into
        # the batch — structurally disqualified (doctrine §4.3).
        if d["kind"] == "composite":
            parts = d["label"].split(":", 1)[1].split("+")
            if any(class_of.get(p) == "annotation" for p in parts):
                continue
        batch.append({"label": d["label"], "kind": d["kind"], "class": "derived",
                      "n_groups": d["n_groups"], "excluded": False, "note": ""})
    order = {"technical": 0, "donor": 1, "condition": 2, "other": 3, "derived": 4}
    batch.sort(key=lambda c: (c["excluded"], order[c["class"]], -c["n_groups"]))
    return {"batch": batch, "cell_type": cell_type}


# ------------------------------------------------------- trial verdict rules

def qualifies(m: dict) -> bool:
    return bool(m.get("harmony_converged")
                and m.get("ilisi_norm_post") is not None
                and m["ilisi_norm_post"] - m["ilisi_norm_pre"] >= ILISI_GAIN_MIN
                and m["clisi_norm_post"] >= m["clisi_norm_pre"] - CLISI_DROP_TOL)


def correction_unnecessary(m: dict) -> bool:
    return bool(m["ilisi_norm_pre"] >= PRE_MIXED_ILISI
                and m.get("pc_regression_r2", 1.0) <= PRE_MIXED_PCR)


# ------------------------------------------------------------------- flow

def _adaptive_n_cells(candidates: dict, override: int | None) -> int:
    if override:
        return override
    viable = [c["n_groups"] for c in candidates["batch"] if not c["excluded"]]
    want = CELLS_PER_BATCH * max(viable, default=1)
    return max(N_CELLS_FLOOR, min(N_CELLS_CAP, want))


def _candidate_spec(adata, cand: dict, outdir: str) -> str:
    """obs column name, or a materialized TSV for a derived candidate."""
    if cand["kind"] == "existing":
        return cand["label"]
    values = obsprofile.derive_values(adata, cand["label"])
    safe = cand["label"].replace(":", "_").replace("+", "_").replace("|", "_")
    path = os.path.join(outdir, "candidates", f"{safe}.tsv")
    return write_values_tsv(path, adata.obs_names, values)


def _run_trial(args, adata, cand: dict, n_cells: int, trial_no: int,
               cell_type_spec: str | None, outdir: str) -> dict:
    spec = _candidate_spec(adata, cand, outdir)
    trial_dir = os.path.join(outdir, f"trial_{trial_no}")
    argv = [args.src, "-o", trial_dir, "--batch-col", spec,
            "--n-cells", str(n_cells), "--seed", str(args.seed)]
    if cell_type_spec:
        argv += ["--cell-type-col", cell_type_spec]
    code = probe.main(argv)
    with open(os.path.join(trial_dir, "result.json")) as fh:
        pr = json.load(fh)
    umap_rel = ""
    src_umap = os.path.join(trial_dir, probe.UMAP_FILENAME)
    if os.path.exists(src_umap):
        umap_rel = f"trial_{trial_no}_umap.png"
        shutil.copyfile(src_umap, os.path.join(outdir, umap_rel))
    m = pr["metrics"]
    if code != 0:
        verdict = "rejected"
    elif correction_unnecessary(m):
        verdict = "correction_unnecessary"
    elif qualifies(m):
        verdict = "adopted"
    else:
        verdict = "rejected"
    return {"batch_col": cand["label"], "spec": spec, "exit_code": code,
            "metrics": {k: m.get(k) for k in
                        ("ilisi_pre", "ilisi_post", "ilisi_norm_pre",
                         "ilisi_norm_post", "clisi_norm_pre", "clisi_norm_post",
                         "clisi_labels", "harmony_converged", "n_batches",
                         "n_batches_sampled", "pc_regression_r2")},
            "umap": umap_rel, "verdict": verdict, "reason": ""}


def _best_cell_type(candidates: dict) -> str | None:
    ct = candidates["cell_type"]
    return ct[0]["label"] if ct else None


def _billing_url() -> str:
    """Where the caller can review the account-level total spend."""
    if os.environ.get("ANTHROPIC_API_KEY"):
        return "https://console.anthropic.com/settings/usage"
    return "https://claude.ai/settings/usage"  # subscription (CLI credentials)


def _tally_llm(res: dict, usage: dict | None) -> None:
    """Aggregate per-round agent usage into metrics.llm (tokens + dollars)."""
    if not usage or not any(v is not None for v in usage.values()):
        return
    llm = res["metrics"].setdefault("llm", {
        "calls": 0, "models": [], "input_tokens": 0, "output_tokens": 0,
        "cache_creation_tokens": 0, "cache_read_tokens": 0,
        "cost_usd": 0.0, "cost_complete": True, "billing_url": _billing_url()})
    llm["calls"] += 1
    model = usage.get("model")
    if model and model not in llm["models"]:
        llm["models"].append(model)
    for k in ("input_tokens", "output_tokens",
              "cache_creation_tokens", "cache_read_tokens"):
        if usage.get(k):
            llm[k] += usage[k]
    if usage.get("cost_usd") is not None:
        llm["cost_usd"] = round(llm["cost_usd"] + usage["cost_usd"], 6)
    else:  # provider did not report a per-call cost (e.g. subscription auth)
        llm["cost_complete"] = False


def _run(args, res: dict, policy) -> int:
    import anndata as ad

    timings = res["metrics"].setdefault("timings", {})
    t0 = time.perf_counter()
    adata = ad.read_h5ad(args.src)
    profile = obsprofile.profile_obs(adata)
    res["profile"] = profile
    candidates = build_candidates(profile)
    res["candidates"] = candidates
    res["thresholds"] = {
        "ilisi_gain_min": ILISI_GAIN_MIN, "clisi_drop_tol": CLISI_DROP_TOL,
        "pre_mixed_ilisi": PRE_MIXED_ILISI, "pre_mixed_pcr": PRE_MIXED_PCR,
        "pathological_tiny_group_frac": PATHOLOGICAL_TINY_GROUP_FRAC}
    timings["profile"] = round(time.perf_counter() - t0, 3)

    best_ct = _best_cell_type(candidates)
    if args.no_probe or policy is None:
        res["columns"] = {"batch": None, "cell_type": _ct_block(best_ct)}
        res["status"] = "needs_review"
        res["reasons"].append(
            "degraded mode (--no-probe or agent unavailable): profile and "
            "ranking produced, no trials run — confirm the columns and re-run, "
            "or consume the ranking directly")
        return EXIT_BLOCKED

    n_cells = _adaptive_n_cells(candidates, args.n_cells)
    res["metrics"]["probe_n_cells"] = n_cells
    trials = res["trials"] = []
    decisions = res["decisions"] = []  # full audit trail, incl. agent tool use
    by_label = {c["label"]: c for c in candidates["batch"]}

    while True:
        state = {"profile": profile, "candidates": candidates, "trials": trials,
                 "thresholds": res["thresholds"], "best_cell_type": best_ct,
                 "probes_left": args.max_probes - len(trials)}
        decision = policy.decide(state)
        action = decision.get("action")
        decisions.append({"action": action,
                          "candidate": decision.get("candidate"),
                          "reason": decision.get("reason", ""),
                          "tools_used": decision.get("tools_used", []),
                          "raw_reply": decision.get("raw_reply", ""),
                          "usage": decision.get("usage", {})})
        _tally_llm(res, decision.get("usage"))
        log.info("policy: %s %s — %s", action, decision.get("candidate"),
                 decision.get("reason"))

        if action == "probe":
            cand = by_label.get(decision.get("candidate"))
            if cand is None or cand["excluded"] or len(trials) >= args.max_probes:
                res["reasons"].append(
                    f"policy proposed invalid/exhausted probe "
                    f"({decision.get('candidate')!r}) — stopping")
                return _blocked(res, best_ct)
            t0 = time.perf_counter()
            trial = _run_trial(args, adata, cand, n_cells, len(trials) + 1,
                               _ct_spec(best_ct), args.outdir)
            trial["reason"] = decision.get("reason", "")
            trials.append(trial)
            timings[f"trial_{len(trials)}"] = round(time.perf_counter() - t0, 3)
            continue

        chosen = decision.get("candidate")
        ct_choice = decision.get("cell_type", best_ct)
        if action in ("adopt", "conclude_unnecessary"):
            if chosen is None:  # tolerate an omitted candidate when unambiguous
                want = ("correction_unnecessary" if action == "conclude_unnecessary"
                        else "adopted")
                matching = [t for t in trials if t["verdict"] == want]
                if len(matching) == 1:
                    chosen = matching[0]["batch_col"]
            trial = next((t for t in trials if t["batch_col"] == chosen), None)
            cand = by_label.get(chosen)
            if trial is None or cand is None:
                res["reasons"].append(
                    f"policy concluded on unknown/unprobed candidate "
                    f"{chosen!r} — stopping")
                return _blocked(res, best_ct)
            value, kind = chosen, "existing"
            if cand["kind"] != "existing":
                value = os.path.join(args.outdir, "batch.tsv")
                shutil.copyfile(trial["spec"], value)
                kind = "derived"
            res["columns"] = {
                "batch": {"value": value, "kind": kind,
                          "correction": ("recommended" if action == "adopt"
                                         else "unnecessary"),
                          "confidence": 0.9,
                          "evidence": decision.get("reason", "")},
                "cell_type": _ct_block(ct_choice)}
            res["status"] = "ok"
            return EXIT_OK
        if action == "conclude_no_batch":
            res["columns"] = {"batch": None, "cell_type": _ct_block(ct_choice)}
            res["columns"]["batch_evidence"] = decision.get("reason", "")
            res["status"] = "ok"
            return EXIT_OK
        # give_up or anything else
        res["reasons"].append(decision.get("reason", "no qualifying candidate"))
        return _blocked(res, best_ct)


def _ct_spec(label: str | None) -> str | None:
    return label


def _ct_block(label: str | None) -> dict | None:
    if label is None:
        return None
    return {"value": label, "kind": "existing", "confidence": 0.7,
            "evidence": "annotation-classified column (name/value heuristics); "
                        "cLISI in trials corroborates"}


def _blocked(res: dict, best_ct: str | None) -> int:
    res["columns"] = {"batch": None, "cell_type": _ct_block(best_ct)}
    res["status"] = "needs_review"
    return EXIT_BLOCKED


def main(argv=None, *, policy="auto") -> int:
    args = build_parser().parse_args(argv)
    logging.basicConfig(level=logging.INFO, stream=sys.stderr,
                        format="%(levelname)s %(name)s: %(message)s")
    params = {"max_probes": args.max_probes, "n_cells": args.n_cells,
              "no_probe": args.no_probe, "seed": args.seed,
              "model": args.model}
    res = new_result("identify_columns", os.path.abspath(args.src), params)

    if policy == "auto":
        try:
            policy = ClaudeAgentPolicy(args.outdir, model=args.model)
        except PolicyUnavailable as exc:
            log.warning("agent unavailable (%s) — degraded mode", exc)
            res["reasons"].append(f"agent unavailable: {exc}")
            policy = None

    t0 = time.perf_counter()
    try:
        code = _run(args, res, policy)
    except PolicyUnavailable as exc:
        res["status"] = "needs_review"
        res["reasons"].append(f"agent failed mid-run: {exc}")
        res.setdefault("columns", {"batch": None, "cell_type": None})
        code = EXIT_BLOCKED
    except Exception as exc:  # noqa: BLE001
        log.exception("unexpected error")
        res["status"] = "error"
        res["reasons"].append(f"{type(exc).__name__}: {exc}")
        code = EXIT_ERROR
    res["metrics"].setdefault("timings", {})["total"] = \
        round(time.perf_counter() - t0, 3)
    res["exit_code"] = code
    try:
        write_result(args.outdir, res)
    except Exception:  # noqa: BLE001
        log.exception("failed to write result.json")
        if code == EXIT_OK:
            code = EXIT_ERROR
    return code


def cli() -> None:
    raise SystemExit(main())


if __name__ == "__main__":
    cli()
