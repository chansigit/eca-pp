"""standardize v0.2 — validate, gate, locate counts, resolve species, harmonize
genes, compute QC, write the standard form.

Flow (spec §4):

    ①  F1       no load: file exists → is_hdf5 → has obs/var structure
    ②  F3-pre   no load: n_cells peeked from HDF5 metadata (millisecond fast-reject)
    ③  load     anndata.read_h5ad
    ④  F3-pre   provisional genes gate on X's nonzero structure (skipped when
                untrustworthy: scaled X, missing X)
    ⑤  F2       counts location & recovery (countsloc: stancounts + 3-layer defence)
    ⑥  F3-final authoritative re-check of both gates on the true counts
    ⑦  F4a      species resolution ladder (T0 --species → T1 infer → T2 --llm → T3 block)
    ⑧  F4       harmonize gene names; drop unmappable features (default); re-gate
    ⑨  F5       authoritative QC obs columns on the final gene space
    ⑩  F7       build (counts layer + lognorm X + provenance) → atomic standardized.h5ad
    ⑪          atomically write result.json; exit per spec §7

Produces ``OUTDIR/standardized.h5ad`` + ``OUTDIR/result.json``. report.md /
qc.png and metacols (F6/F8) join in v0.3; the CLI signature only grows.
"""

from __future__ import annotations

import argparse
import logging
import os
import sys
import time

import anndata as ad
import h5py

from eca_pp import __version__
from eca_pp.standardize import build, countsloc, harmonize
from eca_pp.standardize import species as species_ladder
from eca_pp.standardize.qc import apply_qc, count_n_genes_detected, has_negative
from eca_pp.core.result import (
    EXIT_BLOCKED,
    EXIT_ERROR,
    EXIT_OK,
    EXIT_REJECTED,
    new_result,
    write_result,
)

log = logging.getLogger("eca_pp.standardize")


class _Stop(Exception):
    """Terminal non-ok outcome; carries everything result.json needs."""

    def __init__(self, exit_code: int, status: str, rejected_at: str | None, reasons):
        super().__init__("; ".join(reasons))
        self.exit_code = exit_code
        self.status = status
        self.rejected_at = rejected_at
        self.reasons = list(reasons)


class _Timer:
    """Per-stage wall-time laps, recorded live into result.json's
    ``metrics.timings`` — partial timings survive any early stop."""

    def __init__(self, sink: dict):
        self._t = time.perf_counter()
        self.laps = sink

    def lap(self, label: str) -> None:
        now = time.perf_counter()
        self.laps[label] = round(now - self._t, 3)
        self._t = now


def build_parser() -> argparse.ArgumentParser:
    p = argparse.ArgumentParser(
        prog="eca-pp-standardize",
        description="Standardize one single-cell h5ad sample: validate, hard-gate "
                    "QC, locate counts, resolve species, harmonize gene names, "
                    "compute QC; writes OUTDIR/standardized.h5ad + result.json.")
    p.add_argument("src", help="source .h5ad")
    p.add_argument("-o", "--outdir", required=True, help="output directory")
    p.add_argument("--min-cells", type=int, default=100,
                   help="reject samples below this many cells (default 100)")
    p.add_argument("--min-genes", type=int, default=5000,
                   help="reject samples below this many detected genes (default 5000)")
    p.add_argument("--counts-layer", default=None, metavar="NAME",
                   help="use this layer as counts (skips inference)")
    p.add_argument("--no-gate", action="store_true",
                   help="disable the hard QC gates (metrics still recorded)")
    p.add_argument("--species", default=None, metavar="CODE",
                   help="species (skips inference): hs=human, mm=mouse, rn=rat, "
                        "dr=zebrafish, dm=fruit_fly, ce=c_elegans, "
                        "cyno=cynomolgus, rhesus, marmoset, lemur=mouse_lemur; "
                        "full names accepted too")
    p.add_argument("--llm", action="store_true",
                   help="allow one LLM call as a species-inference fallback "
                        "(off by default; needs ANTHROPIC_API_KEY)")
    p.add_argument("--keep-unmapped", action="store_true",
                   help="keep unmappable features under their original names "
                        "instead of dropping them (default: drop)")
    return p


def _peek_n_obs(f: h5py.File) -> int | None:
    """n_obs from HDF5 structure without loading any matrix; None if undeterminable."""
    x = f.get("X")
    if isinstance(x, h5py.Dataset):
        return int(x.shape[0])
    if isinstance(x, h5py.Group):
        shape = x.attrs.get("shape")
        if shape is not None:
            return int(shape[0])
    obs = f.get("obs")
    if isinstance(obs, h5py.Group):
        idx = obs.attrs.get("_index")
        if isinstance(idx, bytes):
            idx = idx.decode()
        if isinstance(idx, str):
            ds = obs.get(idx)
            if isinstance(ds, h5py.Dataset):
                return int(ds.shape[0])
    return None


def _pre_gate_matrix(adata, counts_layer: str | None):
    """(matrix, trusted, exact) for the provisional genes gate.

    A designated counts layer is exact. Otherwise X is used — its nonzero structure
    equals the counts' under log1p/normalize (both preserve zeros) — unless X is
    missing or carries negatives (scaled data: structure not trustworthy).
    """
    if counts_layer and counts_layer in adata.layers:
        return adata.layers[counts_layer], True, True
    X = adata.X
    if X is None or has_negative(X):
        return None, False, False
    return X, True, False


def _run(args, res: dict) -> int:
    src = res["src"]
    gates = not args.no_gate
    timer = _Timer(res["metrics"].setdefault("timings", {}))

    # ① F1 — input validation, nothing loaded.
    if not os.path.isfile(src):
        raise _Stop(EXIT_REJECTED, "rejected", "input", [f"no such file: {src}"])
    if not h5py.is_hdf5(src):
        raise _Stop(EXIT_REJECTED, "rejected", "input",
                    [f"not an HDF5 file (signature check failed): {src}"])
    with h5py.File(src, "r") as f:
        if "obs" not in f or "var" not in f:
            raise _Stop(EXIT_REJECTED, "rejected", "input",
                        ["HDF5 but not an h5ad: missing obs/var groups"])
        n_obs = _peek_n_obs(f)

    # ② F3-pre — cells gate from metadata (fast reject before any load).
    if n_obs is not None:
        res["metrics"]["n_cells"] = n_obs
        if gates and n_obs < args.min_cells:
            raise _Stop(EXIT_REJECTED, "rejected", "pre_gate",
                        [f"n_cells {n_obs} < min_cells {args.min_cells}"])
    timer.lap("f1_validate")

    # ③ load.
    log.info("loading %s", src)
    adata = ad.read_h5ad(src)
    res["metrics"]["n_cells"] = int(adata.n_obs)
    res["metrics"]["n_vars"] = int(adata.n_vars)
    if gates and adata.n_obs < args.min_cells:  # authoritative (peek may have missed)
        raise _Stop(EXIT_REJECTED, "rejected", "pre_gate",
                    [f"n_cells {adata.n_obs} < min_cells {args.min_cells}"])

    # ④ F3-pre — provisional genes gate.
    pm, trusted, exact = _pre_gate_matrix(adata, args.counts_layer)
    if gates and trusted:
        nd = count_n_genes_detected(pm)
        res["metrics"]["n_genes_detected"] = nd
        if nd < args.min_genes:
            note = "" if exact else " (assessed on X's nonzero structure)"
            raise _Stop(EXIT_REJECTED, "rejected", "pre_gate",
                        [f"n_genes_detected {nd} < min_genes {args.min_genes}{note}"])
    timer.lap("load")

    # ⑤ F2 — counts location & recovery.
    log.info("locating counts ...")
    loc = countsloc.resolve(adata, counts_layer=args.counts_layer)
    res["layers"] = loc.census
    res["metrics"]["x_normalization"] = loc.x_normalization
    if loc.outcome == "blocked":
        raise _Stop(EXIT_BLOCKED, "needs_review", None, loc.blocked)
    if loc.outcome == "unavailable":
        raise _Stop(EXIT_REJECTED, "rejected", "counts_recovery", loc.blocked)
    res["metrics"].update({
        "counts_source": loc.source,
        "counts_integer": loc.counts_integer,
        "counts_name_recognized": loc.name_recognized,
        "counts_adopted_by": loc.adopted_by,
    })
    log.info("counts: %s (via %s)", loc.source, loc.adopted_by)
    timer.lap("f2_counts")

    # ⑥ F3-final — authoritative gates on the true counts.
    n_cells = int(loc.counts.shape[0])
    nd = count_n_genes_detected(loc.counts)
    res["metrics"]["n_cells"] = n_cells
    res["metrics"]["n_genes_detected"] = nd
    if gates:
        reasons = []
        if n_cells < args.min_cells:
            reasons.append(f"n_cells {n_cells} < min_cells {args.min_cells}")
        if nd < args.min_genes:
            reasons.append(f"n_genes_detected {nd} < min_genes {args.min_genes}")
        if reasons:
            raise _Stop(EXIT_REJECTED, "rejected", "final_gate", reasons)
    review = list(loc.needs_review)
    timer.lap("f3_gates")

    # Attach counts as the canonical layer NOW so F4's gene dropping subsets it
    # together with X and every other layer (F7's first half).
    build.attach_counts(adata, loc.counts, loc.source)

    # ⑦ F4a — species resolution ladder.
    try:
        spec = species_ladder.resolve(adata, cli_species=args.species, llm=args.llm)
    except ValueError as exc:  # unknown --species code: a driver mistake, decide again
        raise _Stop(EXIT_BLOCKED, "needs_review", None, [str(exc)])
    res["species"] = spec.as_dict()
    if spec.resolved is None:
        hint = ("LLM fallback also failed — " if args.llm else "") + \
               "decide and re-run with --species CODE"
        raise _Stop(EXIT_BLOCKED, "needs_review", None,
                    [f"species could not be resolved deterministically; {hint} "
                     f"(evidence in result.json species.evidence)"])
    log.info("species: %s (%s, confidence %.2f)",
             spec.resolved, spec.source, spec.confidence)
    timer.lap("f4a_species")

    # ⑧ F4 — harmonize gene names; drop unmappable features (default); re-gate.
    adata, hstats, hreasons = harmonize.harmonize_genes(
        adata, spec.resolved, keep_unmapped=args.keep_unmapped)
    res["metrics"]["harmonization"] = hstats
    res["metrics"]["n_vars"] = int(adata.n_vars)
    review.extend(hreasons)
    nd = count_n_genes_detected(adata.layers["counts"])
    res["metrics"]["n_genes_detected"] = nd
    if gates and nd < args.min_genes:
        raise _Stop(EXIT_REJECTED, "rejected", "final_gate",
                    [f"n_genes_detected {nd} < min_genes {args.min_genes} "
                     f"after dropping unmappable features"])
    log.info("harmonized: %d genes kept, %d genes dropped", hstats["genes_kept"],
             sum(hstats["genes_dropped"].values()))
    timer.lap("f4_harmonize")

    # ⑨ F5 — authoritative QC obs columns on the final gene space. Zero mt/hb
    # matches are NORMAL (pre-filtered matrices, species without RBC hemoglobin)
    # — recorded in metrics.qc, logged, but never flagged for review.
    qstats = apply_qc(adata, spec.resolved)
    res["metrics"]["qc"] = qstats
    if qstats["n_mt_genes"] == 0:
        log.info("0 mitochondrial genes matched — pct_counts_mt is all zeros")
    if qstats["n_hb_genes"] == 0:
        log.info("0 hemoglobin genes matched — pct_counts_hb is all zeros")
    timer.lap("f5_qc")

    # ⑩ F7 — build the standard form and write it atomically.
    raw_dropped = build.finalize(adata, {
        "step_version": __version__,
        "species": spec.resolved,
        "species_code": spec.code or "",
        "species_source": spec.source or "",
        "counts_source": loc.source,
        "counts_adopted_by": loc.adopted_by,
    })
    if raw_dropped:
        res["metrics"]["raw_dropped"] = True
    res["output"] = build.write_h5ad(adata, args.outdir)
    timer.lap("f7_build_write")
    log.info("wrote %s", res["output"])

    if review:
        res["status"] = "needs_review"
        res["reasons"].extend(review)
    else:
        res["status"] = "ok"
    return EXIT_OK


def main(argv=None) -> int:
    args = build_parser().parse_args(argv)
    logging.basicConfig(level=logging.INFO, stream=sys.stderr,
                        format="%(levelname)s %(name)s: %(message)s")
    params = {"min_cells": args.min_cells, "min_genes": args.min_genes,
              "counts_layer": args.counts_layer, "no_gate": args.no_gate,
              "species": args.species, "llm": args.llm,
              "keep_unmapped": args.keep_unmapped}
    res = new_result("standardize", os.path.abspath(args.src), params)

    t0 = time.perf_counter()
    try:
        code = _run(args, res)
    except _Stop as stop:
        res["status"] = stop.status
        res["rejected_at"] = stop.rejected_at
        res["reasons"].extend(stop.reasons)
        code = stop.exit_code
        log.warning("%s: %s", stop.status, "; ".join(stop.reasons))
    except Exception as exc:  # noqa: BLE001 - anything else is exit 1
        log.exception("unexpected error")
        res["status"] = "error"
        res["reasons"].append(f"{type(exc).__name__}: {exc}")
        code = EXIT_ERROR

    res["metrics"].setdefault("timings", {})["total"] = \
        round(time.perf_counter() - t0, 3)
    res["exit_code"] = code
    try:
        path = write_result(args.outdir, res)
        log.info("wrote %s (status=%s, exit=%d)", path, res["status"], code)
    except Exception:  # noqa: BLE001 - a success we can't record is not a success
        log.exception("failed to write result.json")
        if code == EXIT_OK:
            code = EXIT_ERROR
    return code


def cli() -> None:
    raise SystemExit(main())


if __name__ == "__main__":
    cli()
