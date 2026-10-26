"""subset — a cell subset of a standardized product, written in the same standard
form (spec: docs/subset-spec.md).

The source is a product of standardize or subset (its OUTDIR, or the
standardized.h5ad in it, with result.json beside; status ok or needs_review).
A cell is kept when every ``--obs COL=VALUE`` matches and its counts pass
``--min-genes`` / ``--min-counts`` (strictly greater, on ``layers["counts"]``).

Flow:

    ①  source  result.json ready, integer counts layer (I1), species in uns
    ②  filter  per-cell drop reasons → OUTDIR/dropped_cells.csv.gz
    ③  QC      the four QC columns recomputed on the kept cells (I4)
    ④  gates   the source's hard gates re-checked on the subset (I7)
    ⑤  build   X ← lognorm(counts) (I2); obsp dropped (the full product's
               graph) and with --counts-only the other layers, obsm, varm,
               both right after the load
    ⑥          atomic standardized.h5ad + result.json (step "subset",
               derived_from: source, its sha256, step and filters)

A hand-made subset is how tabula-sapiens-clean ended with raw counts in X
(2026-09-07); this step is the supported way to derive one product from another.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import logging
import os
import sys
import time

import anndata as ad
import numpy as np
import pandas as pd
import scipy.sparse as sp

from eca_pp import __version__
from eca_pp.core.atomic_io import atomic_write
from eca_pp.core.layers import layer_names
from eca_pp.core.result import (
    EXIT_ERROR,
    EXIT_OK,
    EXIT_REJECTED,
    RESULT_FILENAME,
    new_result,
    write_result,
)
from eca_pp.standardize import build
from eca_pp.standardize.cli import _Stop, _Timer
from eca_pp.standardize.qc import QC_COLS, apply_qc, count_n_genes_detected, is_integer_matrix

log = logging.getLogger("eca_pp.subset")

DROPPED_FILENAME = "dropped_cells.csv.gz"
READY = ("ok", "needs_review")


def build_parser() -> argparse.ArgumentParser:
    p = argparse.ArgumentParser(
        prog="eca-pp-subset",
        description="Subset the cells of a standardized product and write the "
                    "subset in the same standard form: OUTDIR/standardized.h5ad "
                    "+ result.json + dropped_cells.csv.gz.")
    p.add_argument("src", help="a standardized product: its OUTDIR or its standardized.h5ad")
    p.add_argument("-o", "--outdir", required=True, help="output directory")
    p.add_argument("--obs", action="append", default=[], metavar="COL=VALUE",
                   help="keep cells whose obs COL equals VALUE (repeatable; all must match)")
    p.add_argument("--min-genes", type=int, default=None, metavar="N",
                   help="keep cells with more than N detected genes in the counts layer")
    p.add_argument("--min-counts", type=float, default=None, metavar="N",
                   help="keep cells with more than N total counts")
    p.add_argument("--counts-only", action="store_true",
                   help="keep only the counts layer; drop the other layers, obsm and varm")
    return p


def _source(src: str) -> tuple[str, str]:
    """(standardized.h5ad, result.json) of a product directory or of its h5ad."""
    if os.path.isdir(src):
        return os.path.join(src, build.H5AD_FILENAME), os.path.join(src, RESULT_FILENAME)
    return src, os.path.join(os.path.dirname(src) or ".", RESULT_FILENAME)


def _filters(pairs) -> list[tuple[str, str]]:
    out = []
    for pair in pairs:
        col, sep, value = pair.partition("=")
        if not sep or not col:
            raise _Stop(EXIT_REJECTED, "rejected", "input", [f"--obs needs COL=VALUE, got {pair!r}"])
        out.append((col, value))
    return out


def _sha256(path: str) -> str:
    h = hashlib.sha256()
    with open(path, "rb") as fh:
        for block in iter(lambda: fh.read(1 << 23), b""):
            h.update(block)
    return h.hexdigest()


def _read(h5: str, counts_only: bool):
    """The source. With --counts-only only obs, var, uns and layers['counts'] are read: X is recomputed and the
    rest is dropped anyway, and a tabula-sapiens organ holds four more matrices of its size (12 GB on disk)."""
    if not counts_only:
        return ad.read_h5ad(h5)
    import h5py
    try:
        from anndata.io import read_elem
    except ImportError:  # anndata < 0.11
        from anndata.experimental import read_elem
    with h5py.File(h5, "r") as f:
        if "counts" not in f.get("layers", {}):
            return ad.read_h5ad(h5)  # ① rejects it
        return ad.AnnData(obs=read_elem(f["obs"]), var=read_elem(f["var"]),
                          uns=read_elem(f["uns"]) if "uns" in f else {},
                          layers={"counts": read_elem(f["layers/counts"])})


def drop_reasons(adata, filters, min_genes=None, min_counts=None) -> pd.Series:
    """Per cell, why it is not kept ('' = kept): every failed criterion, joined with '+'."""
    parts = []
    for col, value in filters:
        if col not in adata.obs.columns:
            raise _Stop(EXIT_REJECTED, "rejected", "input", [f"no obs column {col!r}"])
        parts.append(np.where(adata.obs[col].astype(str).to_numpy() == value, "", f"obs {col}!={value}"))
    if min_genes is not None or min_counts is not None:
        C = sp.csr_matrix(adata.layers["counts"], copy=True)
        C.eliminate_zeros()  # stored zeros are not detected genes (as qc.count_n_genes_detected)
        if min_genes is not None:
            parts.append(np.where(np.diff(C.indptr) > min_genes, "", f"n_genes<={min_genes}"))
        if min_counts is not None:
            total = np.asarray(C.sum(axis=1)).ravel()
            parts.append(np.where(total > min_counts, "", f"n_counts<={min_counts:g}"))
    reasons = ["+".join(r for r in row if r) for row in zip(*parts)] if parts else [""] * adata.n_obs
    return pd.Series(reasons, index=adata.obs_names, dtype=object)


def _run(args, res: dict) -> int:
    from stancounts import detect_normalization

    timer = _Timer(res["metrics"].setdefault("timings", {}))
    filters = _filters(args.obs)

    # ① the source must be a ready product of the standard form.
    h5, rj = _source(args.src)
    if not (os.path.isfile(h5) and os.path.isfile(rj)):
        raise _Stop(EXIT_REJECTED, "rejected", "input",
                    [f"not a standardized product ({build.H5AD_FILENAME} and {RESULT_FILENAME}): {args.src}"])
    with open(rj) as fh:
        source = json.load(fh)
    if source.get("status") not in READY:
        raise _Stop(EXIT_REJECTED, "rejected", "input",
                    [f"source status {source.get('status')!r} is not ok or needs_review"])
    log.info("loading %s", h5)
    adata = _read(h5, args.counts_only)
    n_source = int(adata.n_obs)
    if "counts" not in layer_names(adata) or not is_integer_matrix(adata.layers["counts"]):
        raise _Stop(EXIT_REJECTED, "rejected", "input", ["source has no integer layers['counts'] (contract I1)"])
    provenance = dict(adata.uns.get("eca_pp_standardize") or {})
    if not provenance.get("species"):
        raise _Stop(EXIT_REJECTED, "rejected", "input", ["source lacks uns['eca_pp_standardize']['species']"])
    # The full product's graphs do not describe the subset; drop what is not kept before
    # the subset copy (a tabula-sapiens organ carries dense scale_data of ~10 GB).
    for key in list(adata.obsp.keys()):
        del adata.obsp[key]
    if args.counts_only:
        for name in layer_names(adata):
            if name != "counts":
                del adata.layers[name]
        for mapping in (adata.obsm, adata.varm):
            for key in list(mapping.keys()):
                del mapping[key]
    res["species"] = source.get("species")
    res["derived_from"] = {"src": os.path.abspath(h5), "src_sha256": _sha256(h5),
                           "src_step": source.get("step"), "src_step_version": source.get("step_version"),
                           "filters": {"obs": dict(filters), "min_genes": args.min_genes,
                                       "min_counts": args.min_counts},
                           "n_source_cells": n_source}
    timer.lap("load")

    # ② filter; every dropped cell keeps its reason.
    reasons = drop_reasons(adata, filters, args.min_genes, args.min_counts)
    keep = (reasons == "").to_numpy()
    dropped = reasons[~keep].rename_axis("cell_id").rename("reason").reset_index()
    with atomic_write(os.path.join(args.outdir, DROPPED_FILENAME)) as tmp:
        dropped.to_csv(tmp, index=False, compression="gzip")
    res["derived_from"]["n_dropped"] = int(len(dropped))
    adata = adata[keep].copy()
    timer.lap("filter")

    # ③ QC on the kept cells; drop the source's values first, which apply_qc would back up as __original.
    for col in QC_COLS:
        if col in adata.obs.columns:
            del adata.obs[col]
    res["metrics"]["qc"] = apply_qc(adata, provenance["species"])

    # ④ the source's hard gates, on the subset.
    n_cells, nd = int(adata.n_obs), count_n_genes_detected(adata.layers["counts"])
    res["metrics"].update(n_cells=n_cells, n_vars=int(adata.n_vars), n_genes_detected=nd, counts_integer=True)
    gate = source.get("params") or {}
    if not gate.get("no_gate"):
        min_cells, min_genes = gate.get("min_cells", 100), gate.get("min_genes", 5000)
        failed = ([f"n_cells {n_cells} < min_cells {min_cells}"] if n_cells < min_cells else []) + \
                 ([f"n_genes_detected {nd} < min_genes {min_genes}"] if nd < min_genes else [])
        if failed:
            raise _Stop(EXIT_REJECTED, "rejected", "final_gate", failed)
    timer.lap("qc_gates")

    # ⑤ build.
    build.finalize(adata, provenance)
    adata.uns["eca_pp_subset"] = {"step_version": __version__, "src": res["derived_from"]["src"],
                                  "filters": json.dumps(res["derived_from"]["filters"])}
    res["metrics"]["x_normalization"] = detect_normalization(adata.X)
    res["output"] = build.write_h5ad(adata, args.outdir)
    timer.lap("build_write")
    log.info("wrote %s (%d of %d cells)", res["output"], n_cells, n_source)

    res["reasons"] = list(source.get("reasons") or [])  # the source's review notes still hold
    res["status"] = source["status"]
    return EXIT_OK


def main(argv=None) -> int:
    args = build_parser().parse_args(argv)
    logging.basicConfig(level=logging.INFO, stream=sys.stderr,
                        format="%(levelname)s %(name)s: %(message)s")
    params = {"obs": args.obs, "min_genes": args.min_genes, "min_counts": args.min_counts,
              "counts_only": args.counts_only}
    res = new_result("subset", os.path.abspath(args.src), params)

    t0 = time.perf_counter()
    try:
        from eca_pp.core.run_outputs import archive_outputs

        archive_outputs(args.outdir, "subset", _source(args.src)[0])
        code = _run(args, res)
    except _Stop as stop:
        res["status"] = stop.status
        res["rejected_at"] = stop.rejected_at
        res["reasons"] = list(stop.reasons) + res["reasons"]
        code = stop.exit_code
        log.warning("%s: %s", stop.status, "; ".join(stop.reasons))
    except Exception as exc:  # noqa: BLE001 - anything else is exit 1
        log.exception("unexpected error")
        res["status"] = "error"
        res["reasons"].append(f"{type(exc).__name__}: {exc}")
        code = EXIT_ERROR

    res["metrics"].setdefault("timings", {})["total"] = round(time.perf_counter() - t0, 3)
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
