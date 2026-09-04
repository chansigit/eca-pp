"""identify-columns — the project's first declared agent-implemented step
(identify-columns spec). To its caller it is a plain CLI (h5ad in →
result.json out); internally an agent chooses which batch
candidate to probe next and when to conclude, while every tool invocation —
profiling and integration-probe trials — is a deterministic CLI, recorded
round by round in ``trials``.

No selected harness/credential: deterministic policy takes over. --no-probe leaves
batch null. Both finish successfully with structured warnings.
"""

from __future__ import annotations

import argparse
import json
import logging
import os
import re
import sys
import time

from eca_pp.core.atomic_io import copyfile_atomic
from eca_pp.core.colspec import write_values_tsv
from eca_pp.core.result import (
    EXIT_ERROR,
    EXIT_OK,
    EXIT_REJECTED,
    new_result,
    write_result,
)
from eca_pp.identify_columns import obsprofile
from eca_pp.identify_columns.policies import (
    AgentPolicy,
    HeuristicPolicy,
    PolicyUnavailable,
)
from eca_pp.probe import cli as probe

log = logging.getLogger("eca_pp.identify_columns")

# --- decision thresholds (v0.4 defaults; all recorded in result.json) --------
PATHOLOGICAL_TINY_GROUP_FRAC = 0.5   # majority of groups tiny -> column is out
TINY_NOTE_CELL_FRAC = 0.05           # below this, tiny groups only get a note
PRE_MIXED_ILISI = 0.8                # pre-iLISI at/above -> correction unnecessary
PRE_MIXED_PCR = 0.05                 # ...together with PC-regression R2 below
ILISI_GAIN_MIN = 0.05                # normalized iLISI gain required to adopt
CLISI_DROP_TOL = 0.05                # annotated cLISI drop tolerance
PSEUDO_CLISI_DROP_TOL = 0.15         # pseudo-labels are weaker evidence
CLEAR_GAIN_FACTOR = 1.5               # fast path requires metric headroom
CLEAR_CLISI_TOL_FACTOR = 0.5          # preserve half the cLISI budget
PRE_MIXED_ILISI_HEADROOM = 0.05       # fast path likewise for pre-mixed data
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
ANNOTATION_TOKENS = ("celltype", "annotation", "ontology", "class",
                     "lineage", "subtype", "celllabel")
ANNOTATION_EXACT = frozenset({"ct", "celltype", "celltypes", "celllabel",
                              "celllabels"})
# "ann" is a common annotation shorthand but only as a whole affix — as a
# bare substring it would swallow e.g. "channel".
ANNOTATION_AFFIX = re.compile(r"(?i)(^|[_.\s])ann([_.\s]|$)")
# Algorithmic cluster IDs: a valid cLISI label set only as a LAST resort —
# never mistaken for the author's cell-type annotation when one exists.
CLUSTER_TOKENS = ("cluster", "louvain", "leiden")


def build_parser() -> argparse.ArgumentParser:
    p = argparse.ArgumentParser(
        prog="eca-pp-identify-columns",
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
                   help="profile + ranking only; batch=null with warning")
    p.add_argument("--seed", type=int, default=0)
    p.add_argument("--model", default=None,
                   help="agent model ID; default: $ECA_PP_AGENT_MODEL, else "
                        "backend-specific (Doubao for deepseek/openai, Claude for claude)")
    return p


# ------------------------------------------------------------ classification

def _norm(name: str) -> str:
    return name.lower().replace("_", "").replace(".", "").replace(" ", "")


def classify_column(entry: dict) -> str:
    """technical | donor | condition | annotation | cluster | qc_numeric |
    identifier | constant | other — doctrine §4.1. A recognized annotation
    name wins over constant structure because a single-cell-type dataset still
    has valid author metadata. Annotation is also checked before cluster so
    "cluster_annotation" is an annotation, while "seurat_clusters" and
    "leiden" are clusters."""
    n = _norm(entry["column"])
    # Integer-valued QC measurements remain measurements, even if their
    # storage dtype also permits discrete group IDs.
    if (n.startswith(("pctcounts", "percent", "nfeature", "ncount", "ngenes"))
            or n in {"totalcounts", "doubletscore", "scrubletscore", "pctmt", "pcthb"}):
        return "qc_numeric"
    # A named author annotation remains useful metadata even when a dataset
    # contains only one cell type.  It simply cannot serve as a cLISI label.
    if (n in ANNOTATION_EXACT or ANNOTATION_AFFIX.search(entry["column"])
            or any(t in n for t in ANNOTATION_TOKENS)):
        return "annotation"
    if entry["is_constant"]:
        return "constant"
    if entry["is_per_cell_unique"]:
        return "identifier"
    if entry["dtype"] == "float":
        return "qc_numeric"
    for tokens, label in ((CLUSTER_TOKENS, "cluster"),
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


def _numeric_labels(entry: dict) -> bool:
    """All sampled values are bare integers (cluster IDs, not names)."""
    keys = [k for k in entry.get("examples", {}) if k != "<NA>"]
    return bool(keys) and all(k.lstrip("-").isdigit() for k in keys)


CELL_TYPE_CLASS_ORDER = {"annotation": 0, "cluster": 1}


def _batch_tier(cls: str) -> str:
    """Primary technical evidence first; biological/unknown labels later."""
    return "primary" if cls in ("technical", "donor") else "fallback"


def build_candidates(profile: dict) -> dict:
    """{'batch': [...], 'cell_type': [...]} — classified, pre-checked.
    batch: bottom-up ordered (finest viable technical level first).
    cell_type: ranked — named annotation columns first, algorithmic cluster
    columns last, numeric-ID columns after text-labelled ones within a
    class, obs order as the final tiebreak. ``cell_type[0]`` is the
    default the deterministic spine uses; the agent may pick any other."""
    batch, cell_type = [], []
    class_of = {e["column"]: classify_column(e) for e in profile["columns"]}
    for e in profile["columns"]:
        cls = class_of[e["column"]]
        cell_type_eligible = (
            cls == "annotation"
            and not e["is_per_cell_unique"]
            and 1 <= e["n_unique"] <= obsprofile.MAX_GROUPING_CARD)
        cluster_probe_only = (
            cls == "cluster" and obsprofile.is_grouping_candidate(e))
        if cell_type_eligible or cluster_probe_only:
            numeric = _numeric_labels(e)
            usable_for_clisi = e["n_unique"] >= 2
            note = ("algorithmic cluster IDs — probe support only; never "
                    "reported as the author's cell type" if cls == "cluster" else "")
            if numeric and cls == "annotation":
                note = "annotation-named but values are bare integers"
            if not usable_for_clisi:
                note = ((note + "; ") if note else "") + \
                    "constant annotation — valid metadata, unavailable for cLISI"
            cell_type.append({"label": e["column"], "kind": "existing",
                              "class": cls, "n_groups": e["n_unique"],
                              "numeric_labels": numeric,
                              "output_eligible": cls == "annotation",
                              "usable_for_clisi": usable_for_clisi,
                              "note": note})
        if cls in ("technical", "donor", "condition", "other") \
                and obsprofile.is_grouping_candidate(e):
            note = _pathology(e)
            gs = e["group_sizes"]
            cand = {"label": e["column"], "kind": "existing", "class": cls,
                    "tier": _batch_tier(cls),
                    "n_groups": gs["n_groups"],
                    "missing_frac": e["missing_frac"],
                    "tiny_cell_frac": gs["tiny_cell_frac"],
                    "excluded": bool(note), "note": note or ""}
            if not note and 0 < gs["tiny_cell_frac"] <= TINY_NOTE_CELL_FRAC:
                cand["note"] = (f"{gs['n_tiny']} tiny groups "
                                f"({gs['tiny_cell_frac']:.1%} of cells)")
            batch.append(cand)
    for d in profile["derived"]:
        # Combining columns cannot make an ineligible QC/annotation/ID column
        # into a batch factor (doctrine §4.3).
        if d["kind"] == "composite":
            parts = d["label"].split(":", 1)[1].split("+")
            if any(class_of.get(p) not in ("technical", "donor", "condition", "other")
                   for p in parts):
                continue
        tier = "primary" if d["kind"] == "barcode" else "fallback"
        if d["kind"] == "composite":
            parts = d["label"].split(":", 1)[1].split("+")
            if all(class_of.get(p) in ("technical", "donor") for p in parts):
                tier = "primary"
        note = _pathology(d)
        batch.append({"label": d["label"], "kind": d["kind"], "class": "derived",
                      "tier": tier, "n_groups": d["n_groups"],
                      "missing_frac": d.get("missing_frac", 0.0),
                      "tiny_cell_frac": d["group_sizes"]["tiny_cell_frac"],
                      "_equivalent_with": d.get("equivalent_with", []),
                      "excluded": bool(note), "note": note or ""})
    order = {"technical": 0, "donor": 1, "derived": 2,
             "condition": 3, "other": 4}
    batch.sort(key=lambda c: (c["excluded"], c["tier"] != "primary",
                              order[c["class"]], -c["n_groups"]))

    # Collapse equivalent partitions before any expensive probe. Existing
    # columns use the profile relation graph; derived candidates additionally
    # carry exact partition matches computed while their values are available.
    parent = {c["label"]: c["label"] for c in batch}

    def find(label):
        while parent[label] != label:
            parent[label] = parent[parent[label]]
            label = parent[label]
        return label

    def union(left, right):
        if left not in parent or right not in parent:
            return
        left_root, right_root = find(left), find(right)
        if left_root != right_root:
            parent[right_root] = left_root

    for relation in profile.get("relations", []):
        if relation["kind"] == "equivalent":
            union(relation["finer"], relation["coarser"])
    for candidate in batch:
        for existing in candidate.pop("_equivalent_with", []):
            union(candidate["label"], existing)

    representative = {}
    for candidate in batch:
        if candidate["excluded"]:
            continue
        root = find(candidate["label"])
        if root in representative:
            candidate["equivalent_to"] = representative[root]
        else:
            representative[root] = candidate["label"]
    cell_type.sort(key=lambda c: (CELL_TYPE_CLASS_ORDER[c["class"]],
                                  c["numeric_labels"]))  # stable: obs order
    return {"batch": batch, "cell_type": cell_type}


# ------------------------------------------------------- trial verdict rules

def qualifies(m: dict) -> bool:
    clisi_tol = (PSEUDO_CLISI_DROP_TOL
                 if m.get("clisi_labels") == "pseudo" else CLISI_DROP_TOL)
    return bool(m.get("harmony_converged")
                and m.get("ilisi_norm_post") is not None
                and m["ilisi_norm_post"] - m["ilisi_norm_pre"] >= ILISI_GAIN_MIN
                and m["clisi_norm_post"] >= m["clisi_norm_pre"] - clisi_tol)


def correction_unnecessary(m: dict) -> bool:
    return bool(m["ilisi_norm_pre"] >= PRE_MIXED_ILISI
                and m.get("pc_regression_r2", 1.0) <= PRE_MIXED_PCR)


def clear_metric_decision(trial: dict, candidate: dict) -> dict | None:
    """Conclude locally only when a primary trial has ample metric headroom.

    The first agent round still supplies the semantic column choice. Fallback
    candidates, missing batch labels, rejected probes, and borderline metrics
    continue to a second agent round.
    """
    if candidate.get("tier") != "primary" or candidate.get("missing_frac", 0):
        return None
    m = trial["metrics"]
    verdict = trial.get("verdict")
    if verdict == "adopted":
        gain = m["ilisi_norm_post"] - m["ilisi_norm_pre"]
        drop = m["clisi_norm_pre"] - m["clisi_norm_post"]
        tolerance = (PSEUDO_CLISI_DROP_TOL
                     if m.get("clisi_labels") == "pseudo" else CLISI_DROP_TOL)
        if gain < ILISI_GAIN_MIN * CLEAR_GAIN_FACTOR \
                or drop > tolerance * CLEAR_CLISI_TOL_FACTOR:
            return None
        return {
            "action": "adopt", "candidate": trial["batch_col"],
            "cell_type": trial.get("cell_type_col"),
            "reason": (
                f"{trial.get('reason', '').rstrip('.')} — metric fast path: "
                "primary candidate, Harmony converged, "
                f"normalized iLISI gain={gain:.4f} and cLISI drop={drop:.4f} "
                "clear the guarded thresholds"
            ),
        }
    if verdict == "correction_unnecessary":
        if m["ilisi_norm_pre"] < PRE_MIXED_ILISI + PRE_MIXED_ILISI_HEADROOM \
                or m.get("pc_regression_r2", 1.0) > PRE_MIXED_PCR / 2:
            return None
        return {
            "action": "conclude_unnecessary", "candidate": trial["batch_col"],
            "cell_type": trial.get("cell_type_col"),
            "reason": (
                f"{trial.get('reason', '').rstrip('.')} — metric fast path: "
                "primary candidate is already well mixed "
                f"(normalized pre-iLISI={m['ilisi_norm_pre']:.4f}, "
                f"PC regression R2={m.get('pc_regression_r2', 0):.4f})"
            ),
        }
    return None


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
    m = pr["metrics"]
    if code == EXIT_ERROR or pr.get("status") == "error":
        reasons = "; ".join(map(str, pr.get("reasons") or []))
        detail = reasons or "no diagnostic reason was recorded"
        raise RuntimeError(
            f"integration probe failed for {cand['label']!r}: {detail}"
        )
    if code not in (EXIT_OK, EXIT_REJECTED):
        raise RuntimeError(
            f"integration probe returned unexpected exit code {code} "
            f"for {cand['label']!r}"
        )
    if code == EXIT_REJECTED:
        verdict = "rejected"
    elif correction_unnecessary(m):
        verdict = "correction_unnecessary"
    elif qualifies(m):
        verdict = "adopted"
    else:
        verdict = "rejected"
    return {"batch_col": cand["label"], "spec": spec,
            "cell_type_col": cell_type_spec, "exit_code": code,
            "metrics": {k: m.get(k) for k in
                        ("ilisi_pre", "ilisi_post", "ilisi_norm_pre",
                         "ilisi_norm_post", "clisi_norm_pre", "clisi_norm_post",
                         "clisi_labels", "pseudo_label_graph",
                         "cell_type_missing_sampled",
                         "cell_type_coverage_sampled", "n_cells_clisi",
                         "harmony_converged", "n_batches",
                         "n_batches_sampled", "pc_regression_r2", "timings")},
            "verdict": verdict, "reason": ""}


def _best_cell_type(candidates: dict) -> str | None:
    """Top-ranked author annotation; clusters are probe support only."""
    ct = [c for c in candidates["cell_type"] if c["output_eligible"]]
    return ct[0]["label"] if ct else None


def _policy_cell_type(decision: dict, candidates: dict,
                      current: str | None, res: dict) -> str | None:
    """The agent's cell-type choice if it names a listed candidate, else
    ``current``. Applied to EVERY decision (probe included) so the trials'
    cLISI labels follow the agent's choice rather than the heuristic default."""
    choice = decision.get("cell_type")
    if choice and any(c["label"] == choice and c["output_eligible"]
                      for c in candidates["cell_type"]):
        return choice
    if choice:
        _warn(res, "invalid_cell_type_choice",
              "policy proposed a missing or non-author-annotation cell-type column; ignored",
              candidate=choice)
    return current


def _warn(res: dict, code: str, message: str, **details) -> None:
    """Append one structured, non-blocking warning."""
    warning = {"code": code, "message": message}
    if details:
        warning["details"] = details
    if not any(w.get("code") == code and w.get("details") == warning.get("details")
               for w in res.setdefault("warnings", [])):
        res["warnings"].append(warning)


def _active_batch_tier(candidates: dict, trials: list) -> str | None:
    if any(t["verdict"] in ("adopted", "correction_unnecessary")
           for t in trials):
        return None
    tried = {t["batch_col"] for t in trials}
    viable = [c for c in candidates["batch"]
              if not c["excluded"] and not c.get("equivalent_to")
              and c["label"] not in tried]
    if any(c["tier"] == "primary" for c in viable):
        return "primary"
    if any(c["tier"] == "fallback" for c in viable):
        return "fallback"
    return None


def _billing_url() -> str:
    """Where the caller can review the account-level total spend."""
    from eca_pp import agent

    if agent.backend_name() in ("deepseek", "openai"):
        return "https://console.volcengine.com/ark"
    if os.environ.get("ANTHROPIC_API_KEY"):
        return "https://console.anthropic.com/settings/usage"
    return "https://claude.ai/settings/usage"  # subscription (CLI credentials)


def _llm_metrics(res: dict) -> dict:
    llm = res["metrics"].setdefault("llm", {
        "calls": 0, "models": [], "input_tokens": 0, "output_tokens": 0,
        "reasoning_tokens": 0,
        "cache_creation_tokens": 0, "cache_read_tokens": 0,
        "cost_usd": 0.0, "cost_complete": True, "billing_url": _billing_url(),
        "successful_calls": 0, "failed_calls": 0, "timeout_calls": 0,
        "failed_seconds": 0.0, "failures": []})
    for key, default in (("successful_calls", 0), ("failed_calls", 0),
                         ("timeout_calls", 0), ("failed_seconds", 0.0),
                         ("failures", [])):
        llm.setdefault(key, default)
    return llm


def _tally_llm(res: dict, usage: dict | None) -> None:
    """Aggregate one successful agent attempt into metrics.llm."""
    if not usage or not any(v is not None for v in usage.values()):
        return
    llm = _llm_metrics(res)
    llm["calls"] += 1
    llm["successful_calls"] += 1
    model = usage.get("model")
    if model and model not in llm["models"]:
        llm["models"].append(model)
    for k in ("input_tokens", "output_tokens", "reasoning_tokens",
              "cache_creation_tokens", "cache_read_tokens"):
        if usage.get(k):
            llm[k] += usage[k]
    if usage.get("cost_usd") is not None:
        llm["cost_usd"] = round(llm["cost_usd"] + usage["cost_usd"], 6)
    else:  # provider did not report a per-call cost (e.g. subscription auth)
        llm["cost_complete"] = False


def _tally_llm_failure(res: dict, *, error: Exception, elapsed: float,
                       model: str, backend: str) -> None:
    """Record a failed attempt even though it produced no decision/usage."""
    llm = _llm_metrics(res)
    message = str(error)
    timed_out = getattr(error, "kind", "error") == "timeout"
    llm["calls"] += 1
    llm["failed_calls"] += 1
    llm["timeout_calls"] += int(timed_out)
    llm["failed_seconds"] = round(llm["failed_seconds"] + elapsed, 3)
    if model and model not in llm["models"]:
        llm["models"].append(model)
    llm["cost_complete"] = False
    llm["failures"].append({
        "backend": backend, "model": model,
        "kind": "timeout" if timed_out else "error",
        "elapsed_seconds": round(elapsed, 3), "error": message,
    })


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
        "pseudo_clisi_drop_tol": PSEUDO_CLISI_DROP_TOL,
        "pre_mixed_ilisi": PRE_MIXED_ILISI, "pre_mixed_pcr": PRE_MIXED_PCR,
        "pathological_tiny_group_frac": PATHOLOGICAL_TINY_GROUP_FRAC}
    timings["profile"] = round(time.perf_counter() - t0, 3)

    best_ct = _best_cell_type(candidates)
    if args.no_probe:
        res["columns"] = {"batch": None,
                          "cell_type": _ct_block(best_ct, candidates)}
        res["status"] = "ok"
        _warn(res, "probe_disabled",
              "probe disabled: batch left null; profile and cell-type inference produced")
        _warn_if_null_cell_type(res, best_ct, candidates)
        return EXIT_OK
    if policy is None:
        policy = HeuristicPolicy()
        _warn(res, "agent_unavailable",
              "agent unavailable: continued with deterministic policy")

    n_cells = _adaptive_n_cells(candidates, args.n_cells)
    res["metrics"]["probe_n_cells"] = n_cells
    trials = res["trials"] = []
    decisions = res["decisions"] = []  # full audit trail, incl. agent tool use
    by_label = {c["label"]: c for c in candidates["batch"]}

    pending_decision = None
    while True:
        if len(trials) >= args.max_probes and not any(
            trial["verdict"] in ("adopted", "correction_unnecessary")
            for trial in trials
        ):
            return _null_batch_ok(
                res, best_ct, candidates,
                "probe budget exhausted without a qualifying batch",
                code="batch_evidence_insufficient",
            )
        active_tier = _active_batch_tier(candidates, trials)
        state = {"profile": profile, "candidates": candidates, "trials": trials,
                 "thresholds": res["thresholds"], "best_cell_type": best_ct,
                 "active_batch_tier": active_tier,
                 "eligible_batch_candidates": [
                     c["label"] for c in candidates["batch"]
                     if not c["excluded"] and not c.get("equivalent_to")
                     and c["tier"] == active_tier
                     and c["label"] not in {t["batch_col"] for t in trials}],
                 "probes_left": args.max_probes - len(trials)}
        if pending_decision is not None:
            decision = pending_decision
            pending_decision = None
            source = "metric_fast_path"
            res["metrics"]["metric_fast_path"] = True
        else:
            source = "agent" if isinstance(policy, AgentPolicy) else "deterministic"
            decision_no = len(decisions) + 1
            t0 = time.perf_counter()
            try:
                decision = policy.decide(state)
            except PolicyUnavailable as exc:
                failed_elapsed = time.perf_counter() - t0
                if source == "agent":
                    from eca_pp import agent
                    _tally_llm_failure(
                        res, error=exc, elapsed=failed_elapsed,
                        model=agent.model_name(args.model),
                        backend=agent.backend_name(),
                    )
                _warn(res, "agent_failed",
                      "agent failed during the decision loop; continued deterministically",
                      error=str(exc))
                policy = HeuristicPolicy()
                source = "deterministic_fallback"
                decision = policy.decide(state)
            timings[f"decision_{decision_no}"] = round(
                time.perf_counter() - t0, 3)
        action = decision.get("action")
        decisions.append({"action": action,
                          "candidate": decision.get("candidate"),
                          "cell_type": decision.get("cell_type"),
                          "reason": decision.get("reason", ""),
                          "source": source,
                          "tools_used": decision.get("tools_used", []),
                          "raw_reply": decision.get("raw_reply", ""),
                          "usage": decision.get("usage", {})})
        _tally_llm(res, decision.get("usage"))
        log.info("policy: %s %s — %s", action, decision.get("candidate"),
                 decision.get("reason"))
        chosen_ct = _policy_cell_type(decision, candidates, best_ct, res)
        if chosen_ct != best_ct:
            log.info("policy: cell_type %r -> %r", best_ct, chosen_ct)
            best_ct = chosen_ct

        if action == "probe":
            cand = by_label.get(decision.get("candidate"))
            if (cand is None or cand["excluded"] or cand.get("equivalent_to")
                    or cand["tier"] != active_tier
                    or len(trials) >= args.max_probes):
                reason = ("policy proposed an invalid, out-of-tier, or "
                          "exhausted probe")
                _warn(res, "invalid_policy_decision",
                      reason + "; continued deterministically when budget allowed",
                      candidate=decision.get("candidate"), active_tier=active_tier)
                if len(trials) >= args.max_probes:
                    return _null_batch_ok(res, best_ct, candidates,
                                          "probe budget exhausted without a qualifying batch",
                                          code="batch_evidence_insufficient")
                policy = HeuristicPolicy()
                continue
            t0 = time.perf_counter()
            trial = _run_trial(args, adata, cand, n_cells, len(trials) + 1,
                               _ct_spec(best_ct, candidates), args.outdir)
            trial["reason"] = decision.get("reason", "")
            trials.append(trial)
            timings[f"trial_{len(trials)}"] = round(time.perf_counter() - t0, 3)
            if isinstance(policy, AgentPolicy):
                pending_decision = clear_metric_decision(trial, cand)
            continue

        chosen = decision.get("candidate")
        if action in ("adopt", "conclude_unnecessary"):
            if chosen is None:  # tolerate an omitted candidate when unambiguous
                want = ("correction_unnecessary" if action == "conclude_unnecessary"
                        else "adopted")
                matching = [t for t in trials if t["verdict"] == want]
                if len(matching) == 1:
                    chosen = matching[0]["batch_col"]
            trial = next((t for t in trials if t["batch_col"] == chosen), None)
            cand = by_label.get(chosen)
            expected = ("adopted" if action == "adopt"
                        else "correction_unnecessary")
            if trial is None or cand is None or trial["verdict"] != expected:
                reason = (f"policy concluded on unknown, unprobed, or "
                          f"non-qualifying candidate {chosen!r}")
                _warn(res, "invalid_policy_decision", reason,
                      candidate=chosen, action=action)
                policy = HeuristicPolicy()
                if any(t["verdict"] in ("adopted", "correction_unnecessary")
                       for t in trials):
                    continue
                if len(trials) < args.max_probes and active_tier is not None:
                    continue
                return _null_batch_ok(res, best_ct, candidates, reason,
                                      code="batch_evidence_insufficient")
            value, kind = chosen, "existing"
            if cand["kind"] != "existing":
                value = os.path.join(args.outdir, "batch.tsv")
                copyfile_atomic(trial["spec"], value)
                kind = "derived"
            res["columns"] = {
                "batch": {"value": value, "kind": kind,
                          "correction": ("recommended" if action == "adopt"
                                         else "unnecessary"),
                          "confidence": 0.9,
                          "evidence": decision.get("reason", "")},
                "cell_type": _ct_block(best_ct, candidates, trials)}
            res["status"] = "ok"
            if cand["tier"] == "fallback":
                _warn(res, "biological_batch_fallback",
                      "no technical/donor candidate qualified; selected a fallback batch",
                      candidate=chosen, candidate_class=cand["class"])
            if cand.get("missing_frac", 0):
                _warn(res, "selected_batch_has_missing_values",
                      "selected batch column contains missing values",
                      candidate=chosen, missing_frac=cand["missing_frac"])
            _warn_if_null_cell_type(res, best_ct, candidates)
            return EXIT_OK
        qualifying_trials = [t for t in trials if t["verdict"] in
                             ("adopted", "correction_unnecessary")]
        if action in ("conclude_no_batch", "give_up") and qualifying_trials:
            _warn(res, "invalid_policy_decision",
                  "policy ignored a qualifying trial; continued deterministically",
                  qualifying=[t["batch_col"] for t in qualifying_trials])
            policy = HeuristicPolicy()
            continue
        if action == "conclude_no_batch":
            if active_tier is not None and len(trials) < args.max_probes:
                _warn(res, "invalid_policy_decision",
                      "policy concluded no batch before viable candidates were exhausted",
                      active_tier=active_tier)
                policy = HeuristicPolicy()
                continue
            if active_tier is not None:
                return _null_batch_ok(
                    res, best_ct, candidates,
                    "probe budget exhausted before viable candidates were exhausted",
                    code="batch_evidence_insufficient")
            res["columns"] = {"batch": None,
                              "cell_type": _ct_block(best_ct, candidates)}
            res["columns"]["batch_evidence"] = decision.get("reason", "")
            res["status"] = "ok"
            _warn_if_null_cell_type(res, best_ct, candidates)
            return EXIT_OK
        # give_up or anything else
        if active_tier is not None and len(trials) < args.max_probes:
            _warn(res, "invalid_policy_decision",
                  "policy gave up before viable candidates were exhausted; continued deterministically",
                  active_tier=active_tier)
            policy = HeuristicPolicy()
            continue
        reason = decision.get("reason", "no qualifying candidate")
        return _null_batch_ok(res, best_ct, candidates, reason,
                              code="batch_evidence_insufficient")


def _ct_spec(label: str | None, candidates: dict) -> str | None:
    cand = next((c for c in candidates["cell_type"] if c["label"] == label), None)
    return label if cand and cand["usable_for_clisi"] else None


CELL_TYPE_CONFIDENCE = {"annotation": 0.7, "cluster": 0.4}


def _ct_block(label: str | None, candidates: dict,
              trials: list | None = None) -> dict | None:
    if label is None:
        return None
    cand = next((c for c in candidates["cell_type"] if c["label"] == label),
                None)
    cls = cand["class"] if cand else "annotation"
    evidence = (f"{cls}-classified column (name/value heuristics)"
                if cls == "annotation" else
                "algorithmic cluster IDs — no author annotation column found")
    if cand and cand["numeric_labels"]:
        evidence += "; values are bare integers"
    if cand and not cand["usable_for_clisi"]:
        evidence += "; constant annotation, not used for cLISI"
    others = [c["label"] for c in candidates["cell_type"] if c["label"] != label]
    if others:
        evidence += f"; other candidates: {', '.join(others)}"
    if trials and any(
        t.get("cell_type_col") == label
        and t.get("metrics", {}).get("clisi_labels") == "annotated"
        for t in trials
    ):
        evidence += "; used as cLISI labels in the trials"
    return {"value": label, "kind": "existing",
            "confidence": CELL_TYPE_CONFIDENCE.get(cls, 0.5),
            "evidence": evidence}


def _warn_if_null_cell_type(res: dict, best_ct: str | None,
                            candidates: dict) -> None:
    if best_ct is not None:
        return
    clusters = [c["label"] for c in candidates["cell_type"]
                if c["class"] == "cluster"]
    if clusters:
        _warn(res, "cell_type_not_found",
              "no author cell-type annotation found; cluster columns are not reported as cell type",
              cluster_columns=clusters)
    else:
        _warn(res, "cell_type_not_found",
              "no author cell-type annotation found; cell_type is null")


def _null_batch_ok(res: dict, best_ct: str | None, candidates: dict,
                   reason: str, *, code: str) -> int:
    res["columns"] = {"batch": None,
                      "cell_type": _ct_block(best_ct, candidates)}
    res["columns"]["batch_evidence"] = reason
    res["status"] = "ok"
    _warn(res, code, reason)
    _warn_if_null_cell_type(res, best_ct, candidates)
    return EXIT_OK


def main(argv=None, *, policy="auto") -> int:
    args = build_parser().parse_args(argv)
    logging.basicConfig(level=logging.INFO, stream=sys.stderr,
                        format="%(levelname)s %(name)s: %(message)s")
    params = {"max_probes": args.max_probes, "n_cells": args.n_cells,
              "no_probe": args.no_probe, "seed": args.seed,
              "model": args.model}
    res = new_result("identify_columns", os.path.abspath(args.src), params)
    res["warnings"] = []

    if policy == "auto" and args.no_probe:
        policy = None
    elif policy == "auto":
        try:
            policy = AgentPolicy(args.outdir, model=args.model)
        except PolicyUnavailable as exc:
            log.warning("agent unavailable (%s) — deterministic fallback", exc)
            _warn(res, "agent_unavailable",
                  "agent unavailable: continued with deterministic policy",
                  error=str(exc))
            policy = HeuristicPolicy()

    t0 = time.perf_counter()
    try:
        from eca_pp.core.run_outputs import archive_outputs

        archive_outputs(args.outdir, "identify_columns", args.src)
        code = _run(args, res, policy)
    except PolicyUnavailable as exc:
        # Defensive guard: normal mid-loop failures are already converted to
        # deterministic decisions inside _run.
        res["status"] = "ok"
        _warn(res, "agent_failed", "agent failed; output degraded to null columns",
              error=str(exc))
        res.setdefault("columns", {"batch": None, "cell_type": None})
        code = EXIT_OK
    except Exception as exc:
        log.exception("unexpected error")
        res["status"] = "error"
        res["reasons"].append(f"{type(exc).__name__}: {exc}")
        code = EXIT_ERROR
    res["metrics"].setdefault("timings", {})["total"] = \
        round(time.perf_counter() - t0, 3)
    res["exit_code"] = code
    try:
        write_result(args.outdir, res)
    except Exception:
        log.exception("failed to write result.json")
        if code == EXIT_OK:
            code = EXIT_ERROR
    return code


def cli() -> None:
    raise SystemExit(main())


if __name__ == "__main__":
    cli()
