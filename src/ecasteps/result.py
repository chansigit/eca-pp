"""result.json + exit-code convention — the formal interface every step exposes to
its driver (a shell script, Snakemake, an agent, or a human). Spec §7.

Exit codes:
    0  ok (including non-blocking needs_review — output produced, flags recorded)
    1  unexpected error (I/O, memory, bug)      → retry may help
    2  permanent data problem                   → skip, log, never retry
    3  blocked, needs a driver decision         → read the evidence, re-run with params
"""

from __future__ import annotations

import json
import os
from datetime import datetime, timezone

import numpy as np

from ecasteps import __version__
from ecasteps.atomic_io import write_bytes_atomic

SCHEMA_VERSION = 2  # v0.2: adds species / harmonization / qc blocks + "output"

EXIT_OK = 0
EXIT_ERROR = 1
EXIT_REJECTED = 2
EXIT_BLOCKED = 3

RESULT_FILENAME = "result.json"


def new_result(step: str, src: str, params: dict) -> dict:
    """The result.json skeleton. ``status`` starts as ``error`` so a crash that
    somehow skips the outcome-setting path is never reported as success."""
    return {
        "schema_version": SCHEMA_VERSION,
        "step": step,
        "step_version": __version__,
        "status": "error",  # ok | rejected | needs_review | error
        "reasons": [],
        "rejected_at": None,  # input | pre_gate | final_gate | counts_recovery | null
        "exit_code": None,
        "src": src,
        "params": params,
        "output": None,  # OUTDIR/standardized.h5ad once written (F7)
        "species": None,  # F4a: {resolved, code, source, confidence, evidence}
        "metrics": {},
        "layers": [],
        "finished_at": None,
    }


def _json_default(o):
    if isinstance(o, np.integer):
        return int(o)
    if isinstance(o, np.floating):
        return float(o)
    if isinstance(o, np.bool_):
        return bool(o)
    if isinstance(o, np.ndarray):
        return o.tolist()
    return str(o)


def write_result(outdir: str, payload: dict) -> str:
    """Atomically write ``outdir/result.json``; returns the path."""
    payload["finished_at"] = datetime.now(timezone.utc).isoformat(timespec="seconds")
    path = os.path.join(outdir, RESULT_FILENAME)
    data = json.dumps(payload, indent=2, ensure_ascii=False, default=_json_default)
    write_bytes_atomic(path, data.encode("utf-8"))
    return path
