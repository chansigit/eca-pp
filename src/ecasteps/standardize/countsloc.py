"""Counts location — the three-layer defence around ``stancounts.get_counts``. Spec §5.

stancounts' whitelist (counts/count/raw_counts/umi/…) misses oddly-named counts
layers. Rather than guessing from names (or asking an LLM), this module:

1. runs a deterministic **census** of every layer (integer-ness, sparsity, max), and
   when stancounts falls back to log1p *recovery* while unrecognized integer layers
   exist, tries to **prove** a candidate is the true counts via the consistency check
   ``log1p(normalize_total(candidate)) ≈ target`` — a decisive, reproducible test;
2. honours an explicit ``--counts-layer`` override (the driver's decision);
3. reports ambiguity as a **blocked** outcome and doubt as **needs_review**, with the
   full census in result.json so a human or agent can decide without reopening the
   file.

No LLM anywhere in this module.
"""

from __future__ import annotations

from dataclasses import dataclass, field

import numpy as np
import scipy.sparse as sp

import stancounts
from stancounts.counts import DEFAULT_EXCLUDE_LAYERS, DEFAULT_PREFER_LAYERS

from ecasteps.standardize.qc import is_integer_matrix

PREFER = set(DEFAULT_PREFER_LAYERS)
EXCLUDE = set(DEFAULT_EXCLUDE_LAYERS)


@dataclass
class Resolution:
    """Outcome of counts location. ``outcome``: ok | blocked | unavailable."""

    outcome: str = "ok"
    counts: object = None
    source: str = ""
    name_recognized: bool = True
    adopted_by: str = "stancounts"  # stancounts | consistency_check | override
    counts_integer: bool | None = None
    needs_review: list = field(default_factory=list)
    blocked: list = field(default_factory=list)
    census: list = field(default_factory=list)
    x_normalization: dict = field(default_factory=dict)


def _row_dense(M, i) -> np.ndarray:
    row = M[i]
    if sp.issparse(row):
        return np.asarray(row.todense(), dtype=np.float64).ravel()
    return np.asarray(row, dtype=np.float64).ravel()


def consistent_with(candidate, lognorm_target, *, n_rows: int = 50, seed: int = 0,
                    ratio_tol: float = 0.02, structure_tol: float = 1e-3,
                    min_pass_frac: float = 0.98) -> bool:
    """Prove ``candidate`` is the counts matrix behind ``lognorm_target``.

    If target = log1p(candidate / rowsum * T), then per cell: expm1(target_row) has
    the same nonzero set as candidate_row, and their ratio is a single per-row scale
    factor. Both properties are checked on sampled rows — deterministic and cheap.
    """
    if lognorm_target is None or candidate.shape != lognorm_target.shape:
        return False
    n = candidate.shape[0]
    rng = np.random.default_rng(seed)
    rows = np.arange(n) if n <= n_rows else rng.choice(n, n_rows, replace=False)
    passed = 0
    for i in rows:
        c = _row_dense(candidate, i)
        x = np.expm1(_row_dense(lognorm_target, i))
        cm, xm = c > 0, x > 1e-9
        union = int((cm | xm).sum())
        if union == 0:  # empty cell in both — vacuously consistent
            passed += 1
            continue
        if int((cm ^ xm).sum()) > structure_tol * union + 0.5:
            continue  # nonzero structures differ
        both = cm & xm
        if not both.any():
            continue
        r = x[both] / c[both]
        med = float(np.median(r))
        if med > 0 and float(np.max(np.abs(r - med))) <= ratio_tol * med:
            passed += 1
    return passed >= min_pass_frac * len(rows)


def _matrix_stats(M, *, n_sample: int = 200, seed: int = 0) -> dict:
    """dtype / integer-ness / sparsity / max for one matrix (sampled where a full
    pass would be costly on dense data)."""
    if sp.issparse(M):
        size = int(M.shape[0]) * int(M.shape[1])
        nnz = int(M.nnz)
        mx = float(M.data.max()) if nnz else 0.0
        sparsity = 1.0 - (nnz / size if size else 0.0)
    else:
        A = np.asarray(M)
        rng = np.random.RandomState(seed)
        idx = (np.arange(A.shape[0]) if A.shape[0] <= n_sample
               else rng.choice(A.shape[0], n_sample, replace=False))
        sub = A[idx]
        sparsity = float((sub == 0).mean()) if sub.size else 1.0
        mx = float(sub.max()) if sub.size else 0.0
    return {"dtype": str(M.dtype), "is_integer": is_integer_matrix(M),
            "sparsity": round(float(sparsity), 4), "max": mx}


def _census_entry(census: list, name: str) -> dict:
    return next(e for e in census if e["name"] == name)


def _unrecognized_integer_layers(census: list) -> list[str]:
    return [e["name"] for e in census
            if e["is_integer"] and e["name"] not in PREFER and e["name"] not in EXCLUDE]


def resolve(adata, *, counts_layer: str | None = None) -> Resolution:
    """Locate the counts matrix for ``adata`` per the three-layer defence."""
    res = Resolution()
    X = adata.X

    for name in adata.layers:
        if name is None:  # anndata >=0.13 exposes X as layers[None]; X is
            continue      # covered by x_normalization, not the layer census
        res.census.append({"name": name, **_matrix_stats(adata.layers[name]),
                           "consistent_with_X": None})
    try:
        res.x_normalization = dict(stancounts.detect_normalization(X)) if X is not None else {}
    except Exception:  # noqa: BLE001 - the detection is informational only
        res.x_normalization = {}

    # Layer 2 of the defence: explicit override — the driver has decided.
    if counts_layer:
        if counts_layer not in adata.layers:
            res.outcome = "blocked"
            res.blocked.append(
                f"--counts-layer {counts_layer!r} not found; layers present: "
                f"{[n for n in adata.layers if n is not None]}")
            return res
        M = adata.layers[counts_layer]
        if not is_integer_matrix(M):
            res.outcome = "blocked"
            res.blocked.append(
                f"designated layer {counts_layer!r} failed the integer check "
                f"(see the layer census)")
            return res
        res.counts, res.source = M, f"layer:{counts_layer}"
        res.adopted_by = "override"
        res.name_recognized = counts_layer in PREFER
        res.counts_integer = True
        return res

    # Primary path: stancounts.
    try:
        got = stancounts.get_counts(adata)
    except stancounts.CountsUnavailable as exc:
        cands = _unrecognized_integer_layers(res.census)
        if cands:
            res.outcome = "blocked"
            res.blocked.append(
                f"stancounts found no counts ({exc}), but unrecognized integer "
                f"layer(s) exist: {cands} — decide with --counts-layer")
        else:
            res.outcome = "unavailable"
            res.blocked.append(str(exc))
        return res

    src = got["source"]
    if not src.startswith("recovered"):
        # whitelist layer / integer X / raw — unambiguous.
        res.counts, res.source = got["counts"], src
        res.counts_integer = True
        return res

    # Layer 1 of the defence: recovery happened — is a pristine counts layer being
    # overlooked? Prove candidates against the matrix recovery would have reversed.
    cands = _unrecognized_integer_layers(res.census)
    target = (X if src == "recovered"
              else adata.layers.get(src.split("recovered:layer:", 1)[1]))
    proven = []
    for name in cands:
        ok = consistent_with(adata.layers[name], target)
        _census_entry(res.census, name)["consistent_with_X"] = bool(ok)
        if ok:
            proven.append(name)

    if len(proven) == 1:
        name = proven[0]
        res.counts, res.source = adata.layers[name], f"layer:{name}"
        res.adopted_by = "consistency_check"
        res.name_recognized = False
        res.counts_integer = True
        return res
    if len(proven) > 1:
        res.outcome = "blocked"
        res.blocked.append(
            f"multiple integer layers are consistent with the log-normalized "
            f"matrix: {proven} — pick one with --counts-layer")
        return res

    # No candidate proven: accept the recovery, flagging doubt when candidates exist.
    res.counts, res.source = got["counts"], src
    res.counts_integer = is_integer_matrix(res.counts)
    if cands:
        res.needs_review.append(
            f"counts were RECOVERED by reversing log1p while unrecognized integer "
            f"layer(s) {cands} exist (consistency check failed for all) — verify, "
            f"or re-run with --counts-layer")
    return res
