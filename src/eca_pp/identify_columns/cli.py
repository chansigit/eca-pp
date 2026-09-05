"""identify-columns — name the batch column and the cell-type column of a
standardized h5ad (identify-columns spec).

Flow: profile obs (deterministic) → ONE classification call in which the
model reads every column's value-count table and ranks up to three batch
candidates plus the author's cell-type column → integration probes verify
the ranked batch candidates in order (at most ``--max-probes``) → result.json.
Name heuristics only feed the prompt, the fallback when no model is
available, and the "probeable" allow-list.
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
    AgentClassifier,
    HeuristicClassifier,
    PolicyUnavailable,
)
from eca_pp.probe import cli as probe

log = logging.getLogger("eca_pp.identify_columns")

# --- decision thresholds (recorded in result.json) ---------------------------
PATHOLOGICAL_TINY_GROUP_FRAC = 0.5   # majority of groups tiny -> column is out
TINY_NOTE_CELL_FRAC = 0.05           # below this, tiny groups only get a note
PRE_MIXED_ILISI = 0.8                # pre-iLISI at/above -> correction unnecessary
PRE_MIXED_PCR = 0.05                 # ...together with PC-regression R2 below
ILISI_GAIN_MIN = 0.05                # normalized iLISI gain required to adopt
CLISI_DROP_TOL = 0.05                # annotated cLISI drop tolerance
PSEUDO_CLISI_DROP_TOL = 0.15         # pseudo-labels are weaker evidence
CELLS_PER_BATCH = 50                 # adaptive sampling: expected cells/batch
N_CELLS_FLOOR, N_CELLS_CAP = 5000, 30000
MAX_PROBES = 2                       # the classifier ranks; probes verify in order

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
# "ann" is a common annotation shorthand as a whole affix, optionally followed
# by a date/version tag (ann0608, ann_v2) — as a bare substring it would
# swallow e.g. "channel".
ANNOTATION_AFFIX = re.compile(r"(?i)(^|[_.\s])ann(?:\d+|_?v?\d+)?([_.\s]|$)")
# Algorithmic cluster IDs: never the author's cell-type annotation.
CLUSTER_TOKENS = ("cluster", "louvain", "leiden")
# Per-cell biological states (cell-cycle phase, ...): never a batch factor.
STATE_TOKENS = ("cellcycle", "phase", "cyclestate")


def build_parser() -> argparse.ArgumentParser:
    p = argparse.ArgumentParser(
        prog="eca-pp-identify-columns",
        description="Identify the batch column and cell-type column of a "
                    "standardized h5ad: one model classification over the "
                    "obs profile, verified by small integration trials; "
                    "writes OUTDIR/result.json (+ batch.tsv for a derived "
                    "batch column).")
    p.add_argument("src", help="standardized .h5ad")
    p.add_argument("-o", "--outdir", required=True)
    p.add_argument("--max-probes", type=int, default=MAX_PROBES,
                   help=f"verify at most this many ranked batch candidates "
                        f"(default {MAX_PROBES})")
    p.add_argument("--n-cells", type=int, default=None,
                   help="probe subsample size (default: adaptive, "
                        "clamp(50×max_batches, 5000, 30000))")
    p.add_argument("--no-probe", action="store_true",
                   help="profile + classification only; batch=null with warning")
    p.add_argument("--seed", type=int, default=0)
    p.add_argument("--model", default=None,
                   help="agent model ID; default: $ECA_PP_AGENT_MODEL, else "
                        "backend-specific (Doubao for deepseek/openai, Claude for claude)")
    return p


# ------------------------------------------------------------ classification

def _norm(name: str) -> str:
    return name.lower().replace("_", "").replace(".", "").replace(" ", "")


def classify_column(entry: dict) -> str:
    """Name heuristic: technical | donor | condition | annotation | cluster |
    state | qc_numeric | identifier | constant | other. A hint for the model
    and the deterministic fallback; it also defines the probeable set (state /
    annotation / cluster / qc / identifier columns are never probed)."""
    n = _norm(entry["column"])
    if (n.startswith(("pctcounts", "percent", "nfeature", "ncount", "ngenes"))
            or n in {"totalcounts", "doubletscore", "scrubletscore", "pctmt", "pcthb"}):
        return "qc_numeric"
    if (n in ANNOTATION_EXACT or ANNOTATION_AFFIX.search(entry["column"])
            or any(t in n for t in ANNOTATION_TOKENS)):
        return "annotation"
    if entry["is_constant"]:
        return "constant"
    if entry["is_per_cell_unique"]:
        return "identifier"
    if entry["dtype"] == "float":
        return "qc_numeric"
    if any(t in n for t in STATE_TOKENS):
        return "state"
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


CELL_TYPE_CLASS_ORDER = {"annotation": 0, "other": 1, "cluster": 2}


def build_candidates(profile: dict) -> dict:
    """{'batch': [...], 'cell_type': [...]} from name heuristics + group health.

    batch: every grouping column / derived candidate with its heuristic class,
    pathology exclusion, nesting parents and equivalence collapse — the
    probeable allow-list and the fallback ranking (technical → donor → derived
    → condition → other, more groups first).
    cell_type: annotation-named columns, unplaced text columns (class
    "other"), and cluster columns (probe support only)."""
    batch, cell_type = [], []
    class_of = {e["column"]: classify_column(e) for e in profile["columns"]}
    for e in profile["columns"]:
        cls = class_of[e["column"]]
        text_like = (e["dtype"] in ("categorical", "string")
                     and not e["is_per_cell_unique"]
                     and 1 <= e["n_unique"] <= obsprofile.MAX_GROUPING_CARD)
        if cls == "annotation" and not e["is_per_cell_unique"] \
                and 1 <= e["n_unique"] <= obsprofile.MAX_GROUPING_CARD:
            ct_cls = "annotation"
        elif cls == "cluster" and obsprofile.is_grouping_candidate(e):
            ct_cls = "cluster"
        elif cls in ("other", "condition") and text_like and not _numeric_labels(e):
            ct_cls = "other"
        else:
            ct_cls = None
        if ct_cls:
            numeric = _numeric_labels(e)
            usable = e["n_unique"] >= 2
            note = {"cluster": "algorithmic cluster IDs — probe support only; "
                               "never reported as the author's cell type",
                    "other": "text labels not recognized as an annotation by "
                             "name; judged from the sampled values",
                    "annotation": ""}[ct_cls]
            if numeric and ct_cls == "annotation":
                note = "annotation-named but values are bare integers"
            if not usable:
                note = ((note + "; ") if note else "") + \
                    "constant annotation — valid metadata, unavailable for cLISI"
            cell_type.append({"label": e["column"], "kind": "existing",
                              "class": ct_cls, "n_groups": e["n_unique"],
                              "numeric_labels": numeric,
                              "usable_for_clisi": usable, "note": note})
        if cls in ("technical", "donor", "condition", "other") \
                and obsprofile.is_grouping_candidate(e):
            note = _pathology(e)
            gs = e["group_sizes"]
            cand = {"label": e["column"], "kind": "existing", "class": cls,
                    "n_groups": gs["n_groups"],
                    "missing_frac": e["missing_frac"],
                    "tiny_cell_frac": gs["tiny_cell_frac"],
                    "excluded": bool(note), "note": note or ""}
            if not note and 0 < gs["tiny_cell_frac"] <= TINY_NOTE_CELL_FRAC:
                cand["note"] = (f"{gs['n_tiny']} tiny groups "
                                f"({gs['tiny_cell_frac']:.1%} of cells)")
            batch.append(cand)
    for d in profile["derived"]:
        if d["kind"] == "composite":
            # Combining columns cannot turn a QC/annotation/ID column into a
            # batch factor.
            parts = d["label"].split(":", 1)[1].split("+")
            if any(class_of.get(p) not in ("technical", "donor", "condition", "other")
                   for p in parts):
                continue
        note = _pathology(d)
        batch.append({"label": d["label"], "kind": d["kind"], "class": "derived",
                      "n_groups": d["n_groups"],
                      "missing_frac": d.get("missing_frac", 0.0),
                      "tiny_cell_frac": d["group_sizes"]["tiny_cell_frac"],
                      "_equivalent_with": d.get("equivalent_with", []),
                      "excluded": bool(note), "note": note or ""})
    order = {"technical": 0, "donor": 1, "derived": 2, "condition": 3, "other": 4}
    batch.sort(key=lambda c: (c["excluded"], order[c["class"]], -c["n_groups"]))

    # Equivalent partitions collapse onto the first-listed representative.
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

    nested_within: dict[str, list[str]] = {}
    for relation in profile.get("relations", []):
        if relation["kind"] == "nested":
            nested_within.setdefault(relation["finer"], []).append(relation["coarser"])
    for candidate in batch:
        parents = nested_within.get(candidate["label"], [])
        if parents:
            candidate["nested_within"] = [
                {"column": parent, "class": class_of.get(parent, "derived")}
                for parent in parents]
    cell_type.sort(key=lambda c: (CELL_TYPE_CLASS_ORDER[c["class"]],
                                  c["numeric_labels"]))  # stable: obs order
    return {"batch": batch, "cell_type": cell_type}


def build_evidence(profile: dict, candidates: dict) -> dict:
    """The compact table the classifier reads: one row per column with its
    value-count table, plus relations, derived candidates and pathologies."""
    heuristic = {e["column"]: classify_column(e) for e in profile["columns"]}
    probeable = {c["label"]: c for c in candidates["batch"]}
    rows = []
    for e in profile["columns"]:
        row = {"column": e["column"], "dtype": e["dtype"],
               "n_unique": e["n_unique"], "missing_frac": e["missing_frac"],
               "heuristic_class": heuristic[e["column"]],
               "value_counts": e["examples"]}
        if e["is_per_cell_unique"]:
            row["per_cell_unique"] = True
        if e["is_constant"]:
            row["constant"] = True
        gs = e.get("group_sizes")
        if gs:
            row["group_sizes"] = {k: gs[k] for k in ("n_groups", "min", "median",
                                                     "max", "n_tiny")}
        cand = probeable.get(e["column"])
        if cand:
            row["probeable_as_batch"] = not cand["excluded"]
            if cand["note"]:
                row["note"] = cand["note"]
            if cand.get("nested_within"):
                row["nested_within"] = cand["nested_within"]
            if cand.get("equivalent_to"):
                row["equivalent_to"] = cand["equivalent_to"]
        rows.append(row)
    derived = []
    for c in candidates["batch"]:
        if c["kind"] != "existing":
            derived.append({"label": c["label"], "kind": c["kind"],
                            "n_groups": c["n_groups"],
                            "probeable_as_batch": not c["excluded"],
                            "note": c["note"],
                            "equivalent_to": c.get("equivalent_to")})
    return {"n_obs": profile["n_obs"], "columns": rows,
            "relations": profile.get("relations", []),
            "derived_candidates": derived,
            "probeable_batch_columns": sorted(
                c["label"] for c in candidates["batch"] if not c["excluded"])}


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
            "verdict": verdict,
            "probe_reasons": pr.get("reasons") or [],
            "reason": ""}


def _warn(res: dict, code: str, message: str, **details) -> None:
    """Append one structured, non-blocking warning."""
    warning = {"code": code, "message": message}
    if details:
        warning["details"] = details
    if not any(w.get("code") == code and w.get("details") == warning.get("details")
               for w in res.setdefault("warnings", [])):
        res["warnings"].append(warning)


def _billing_url() -> str:
    """Where the caller can review the account-level total spend."""
    from eca_pp import agent

    if agent.backend_name() in ("deepseek", "openai"):
        return "https://console.volcengine.com/ark"
    if os.environ.get("ANTHROPIC_API_KEY"):
        return "https://console.anthropic.com/settings/usage"
    return "https://claude.ai/settings/usage"  # subscription (CLI credentials)


def _llm_metrics(res: dict) -> dict:
    return res["metrics"].setdefault("llm", {
        "calls": 0, "models": [], "input_tokens": 0, "output_tokens": 0,
        "reasoning_tokens": 0,
        "cache_creation_tokens": 0, "cache_read_tokens": 0,
        "cost_usd": 0.0, "cost_complete": True, "billing_url": _billing_url(),
        "successful_calls": 0, "failed_calls": 0, "timeout_calls": 0,
        "failed_seconds": 0.0, "failures": []})


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
    """Record a failed attempt even though it produced no answer/usage."""
    llm = _llm_metrics(res)
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
        "elapsed_seconds": round(elapsed, 3), "error": str(error),
    })


def _classify(args, res: dict, classifier, state: dict) -> tuple[dict, str]:
    """One classification; a failing agent degrades to the name heuristics."""
    source = "agent" if isinstance(classifier, AgentClassifier) else "deterministic"
    t0 = time.perf_counter()
    try:
        answer = classifier.classify(state)
    except PolicyUnavailable as exc:
        elapsed = time.perf_counter() - t0
        if source == "agent":
            from eca_pp import agent
            _tally_llm_failure(res, error=exc, elapsed=elapsed,
                               model=agent.model_name(args.model),
                               backend=agent.backend_name())
        _warn(res, "agent_failed",
              "agent classification failed; continued with name heuristics",
              error=str(exc))
        answer = HeuristicClassifier().classify(state)
        source = "deterministic_fallback"
    res["metrics"]["timings"]["classification"] = round(time.perf_counter() - t0, 3)
    _tally_llm(res, answer.get("usage"))
    return answer, source


def _cell_type_entry(label: str | None, profile: dict, candidates: dict) -> dict | None:
    """Per-column facts for the chosen cell type (None when label is None or
    not an obs column)."""
    if label is None:
        return None
    entry = next((e for e in profile["columns"] if e["column"] == label), None)
    if entry is None or entry["is_per_cell_unique"]:
        return None
    cand = next((c for c in candidates["cell_type"] if c["label"] == label), None)
    return {"label": label, "n_groups": entry["n_unique"],
            "usable_for_clisi": entry["n_unique"] >= 2,
            "numeric_labels": _numeric_labels(entry),
            "heuristic_class": cand["class"] if cand else classify_column(entry)}


def _ct_block(ct: dict | None, classification: dict, source: str,
              candidates: dict, trials: list) -> dict | None:
    if ct is None:
        return None
    evidence = classification.get("cell_type_reason") or ""
    if ct["numeric_labels"]:
        evidence += "; values are bare integers"
    if not ct["usable_for_clisi"]:
        evidence += "; constant annotation, not used for cLISI"
    others = [c["label"] for c in candidates["cell_type"]
              if c["label"] != ct["label"] and c["class"] != "cluster"]
    if others:
        evidence += f"; other candidates: {', '.join(others)}"
    if any(t.get("cell_type_col") == ct["label"]
           and t.get("metrics", {}).get("clisi_labels") == "annotated"
           for t in trials):
        evidence += "; used as cLISI labels in the trials"
    return {"value": ct["label"], "kind": "existing",
            "confidence": 0.8 if source == "agent" else 0.6,
            "evidence": evidence.strip("; ")}


def _warn_if_null_cell_type(res: dict, ct: dict | None, candidates: dict) -> None:
    if ct is not None:
        return
    clusters = [c["label"] for c in candidates["cell_type"] if c["class"] == "cluster"]
    if clusters:
        _warn(res, "cell_type_not_found",
              "no author cell-type annotation found; cluster columns are not reported as cell type",
              cluster_columns=clusters)
    else:
        _warn(res, "cell_type_not_found",
              "no author cell-type annotation found; cell_type is null")


def _finish(res: dict, batch_block: dict | None, ct_block: dict | None,
            batch_evidence: str | None = None) -> int:
    res["columns"] = {"batch": batch_block, "cell_type": ct_block}
    if batch_block is None and batch_evidence:
        res["columns"]["batch_evidence"] = batch_evidence
    res["status"] = "ok"
    return EXIT_OK


def _run(args, res: dict, classifier) -> int:
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

    # ① one classification over the whole profile
    if classifier is None:
        classifier = HeuristicClassifier()
        _warn(res, "agent_unavailable",
              "agent unavailable: continued with name heuristics")
    state = {"profile": profile, "candidates": candidates,
             "evidence": build_evidence(profile, candidates),
             "heuristic_class": {e["column"]: classify_column(e)
                                 for e in profile["columns"]}}
    answer, source = _classify(args, res, classifier, state)
    classification = {
        "source": source,
        "batch_ranked": answer.get("batch_ranked") or [],
        "cell_type": answer.get("cell_type"),
        "cell_type_reason": answer.get("cell_type_reason", ""),
        "columns": answer.get("columns") or {},
        "notes": answer.get("notes", ""),
        "tools_used": answer.get("tools_used", []),
        "raw_reply": answer.get("raw_reply", ""),
        "usage": answer.get("usage", {}),
    }
    res["classification"] = classification
    # Audit-trail shape kept for downstream summaries (one entry per call).
    res["decisions"] = [{
        "action": "classify",
        "candidate": (classification["batch_ranked"][0]["column"]
                      if classification["batch_ranked"] else None),
        "cell_type": classification["cell_type"], "source": source,
        "reason": classification["notes"] or classification["cell_type_reason"],
        "tools_used": classification["tools_used"],
        "raw_reply": classification["raw_reply"], "usage": classification["usage"]}]
    log.info("classification (%s): batch ranked %s; cell_type %r", source,
             [b["column"] for b in classification["batch_ranked"]],
             classification["cell_type"])

    # ② cell type: the classifier's column, sanity-checked against the profile
    ct = _cell_type_entry(classification["cell_type"], profile, candidates)
    if classification["cell_type"] is not None and ct is None:
        _warn(res, "invalid_cell_type_choice",
              "classifier named a column that is missing or per-cell unique; ignored",
              candidate=classification["cell_type"])
    if ct is not None and ct["heuristic_class"] != "annotation":
        _warn(res, "cell_type_identified_from_values",
              "cell-type column was not recognized by name; identified from its values",
              candidate=ct["label"], reason=classification["cell_type_reason"])
    _warn_if_null_cell_type(res, ct, candidates)
    ct_spec = ct["label"] if ct and ct["usable_for_clisi"] else None

    # ③ verify the ranked batch candidates in order
    trials = res["trials"] = []
    by_label = {c["label"]: c for c in candidates["batch"]}
    ranked = [b for b in classification["batch_ranked"]
              if b.get("column") in by_label and not by_label[b["column"]]["excluded"]]
    dropped = [b.get("column") for b in classification["batch_ranked"]
               if b not in ranked]
    if dropped:
        _warn(res, "invalid_batch_choice",
              "classifier ranked columns that are not probeable; skipped",
              candidates=dropped)
    if args.no_probe:
        _warn(res, "probe_disabled",
              "probe disabled: batch left null; profile and classification produced")
        return _finish(res, None, _ct_block(ct, classification, source, candidates, trials),
                       "; ".join(f"{b['column']}: {b['reason']}" for b in ranked) or None)
    if not ranked:
        _warn(res, "no_batch_candidate",
              "classifier found no plausible batch structure in obs")
        return _finish(res, None, _ct_block(ct, classification, source, candidates, trials),
                       classification["notes"] or "no plausible batch column")
    if adata.n_obs < probe.MIN_CELLS:
        reason = (f"dataset has {adata.n_obs} cells (< {probe.MIN_CELLS}); too small "
                  f"for integration trials, so the ranked batch candidates "
                  f"{[b['column'] for b in ranked]} were not probed")
        _warn(res, "dataset_too_small_to_probe", reason)
        return _finish(res, None, _ct_block(ct, classification, source, candidates, trials), reason)

    n_cells = _adaptive_n_cells(candidates, args.n_cells)
    res["metrics"]["probe_n_cells"] = n_cells
    verdicts = []
    for b in ranked[:max(args.max_probes, 0)]:
        cand = by_label[b["column"]]
        t1 = time.perf_counter()
        trial = _run_trial(args, adata, cand, n_cells, len(trials) + 1,
                           ct_spec, args.outdir)
        trial["reason"] = b.get("reason", "")
        trial["class"] = b.get("class") or cand["class"]
        trials.append(trial)
        timings[f"trial_{len(trials)}"] = round(time.perf_counter() - t1, 3)
        log.info("trial %d: %s -> %s", len(trials), cand["label"], trial["verdict"])
        verdicts.append(f"{cand['label']}: {trial['verdict']}")
        if trial["verdict"] not in ("adopted", "correction_unnecessary"):
            continue
        value, kind = cand["label"], "existing"
        if cand["kind"] != "existing":
            value = os.path.join(args.outdir, "batch.tsv")
            copyfile_atomic(trial["spec"], value)
            kind = "derived"
        m = trial["metrics"]
        evidence = (f"{b.get('reason', '')}".rstrip(".") + " — probe: " + (
            f"normalized iLISI {m['ilisi_norm_pre']} → {m['ilisi_norm_post']}, "
            f"cLISI {m['clisi_norm_pre']} → {m['clisi_norm_post']} ({m['clisi_labels']} labels)"
            if trial["verdict"] == "adopted" else
            f"already mixed (normalized pre-iLISI {m['ilisi_norm_pre']}, "
            f"PC regression R2 {m.get('pc_regression_r2')})"))
        batch_block = {"value": value, "kind": kind,
                       "correction": ("recommended" if trial["verdict"] == "adopted"
                                      else "unnecessary"),
                       "confidence": 0.9, "evidence": evidence}
        if trial["class"] in ("condition", "other"):
            _warn(res, "biological_batch_fallback",
                  "selected batch is an experimental condition or unclassified grouping, "
                  "not a technical/donor factor",
                  candidate=cand["label"], candidate_class=trial["class"])
        if cand.get("missing_frac", 0):
            _warn(res, "selected_batch_has_missing_values",
                  "selected batch column contains missing values",
                  candidate=cand["label"], missing_frac=cand["missing_frac"])
        return _finish(res, batch_block,
                       _ct_block(ct, classification, source, candidates, trials))
    untried = [b["column"] for b in ranked[len(trials):]]
    reason = "no ranked batch candidate qualified in the probes (" + "; ".join(verdicts) + ")"
    if untried:
        reason += f"; not probed within --max-probes: {untried}"
    _warn(res, "batch_evidence_insufficient", reason)
    return _finish(res, None, _ct_block(ct, classification, source, candidates, trials), reason)


def main(argv=None, *, classifier="auto") -> int:
    args = build_parser().parse_args(argv)
    logging.basicConfig(level=logging.INFO, stream=sys.stderr,
                        format="%(levelname)s %(name)s: %(message)s")
    params = {"max_probes": args.max_probes, "n_cells": args.n_cells,
              "no_probe": args.no_probe, "seed": args.seed,
              "model": args.model}
    res = new_result("identify_columns", os.path.abspath(args.src), params)
    res["warnings"] = []

    if classifier == "auto":
        try:
            classifier = AgentClassifier(args.outdir, model=args.model)
        except PolicyUnavailable as exc:
            log.warning("agent unavailable (%s) — name heuristics", exc)
            _warn(res, "agent_unavailable",
                  "agent unavailable: continued with name heuristics",
                  error=str(exc))
            classifier = HeuristicClassifier()

    t0 = time.perf_counter()
    try:
        from eca_pp.core.run_outputs import archive_outputs

        archive_outputs(args.outdir, "identify_columns", args.src)
        code = _run(args, res, classifier)
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
