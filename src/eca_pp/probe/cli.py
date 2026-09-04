"""integration-probe — the deterministic small-scale integration trial
(identify-columns spec §5).

A seeded group-preserving subsample → HVG → PCA → Harmony → iLISI/cLISI.
A non-converging or crashing Harmony run is a LEGITIMATE observation (status
ok, ``harmony_converged: false``) — it is exactly the pathology this tool
exists to detect. Heavy dependencies (scanpy/harmonypy) are imported
lazily and ship in the ``[probe]`` extra.
"""

from __future__ import annotations

import argparse
import logging
import os
import sys
import time

import numpy as np
import pandas as pd

from eca_pp.core.colspec import ColumnSpecError, resolve_spec
from eca_pp.core.result import EXIT_ERROR, EXIT_OK, EXIT_REJECTED, new_result, write_result

log = logging.getLogger("eca_pp.probe")

DEFAULT_N_CELLS = 5000
DEFAULT_N_HVG = 2000
N_PCS = 30
LISI_PERPLEXITY = 30
LEIDEN_RESOLUTION = 1.0
PSEUDO_LABEL_GRAPH = "knn_gauss_sklearn"
MIN_CELLS = 300
MIN_SAMPLED_PER_BATCH = 5


def build_parser() -> argparse.ArgumentParser:
    p = argparse.ArgumentParser(
        prog="eca-pp-integration-probe",
        description="Small-scale integration trial: subsample, run Harmony on "
                    "the given batch column and report iLISI/cLISI. "
                    "Diagnoses a batch-column choice; not a production "
                    "integration.")
    p.add_argument("src", help="standardized .h5ad")
    p.add_argument("-o", "--outdir", required=True)
    p.add_argument("--batch-col", required=True, metavar="SPEC",
                   help="obs column name, or a TSV value file (cell_id<TAB>value)")
    p.add_argument("--cell-type-col", default=None, metavar="SPEC",
                   help="optional cell-type labels (column or TSV); without it "
                        "cLISI uses Leiden pseudo-labels")
    p.add_argument("--n-cells", type=int, default=DEFAULT_N_CELLS)
    p.add_argument("--n-hvg", type=int, default=DEFAULT_N_HVG)
    p.add_argument("--seed", type=int, default=0)
    return p


class _Stop(Exception):
    def __init__(self, exit_code: int, reasons: list):
        super().__init__("; ".join(reasons))
        self.exit_code = exit_code
        self.reasons = reasons


def _group_stats(values: pd.Series) -> dict:
    sizes = values.value_counts(dropna=True)
    if not len(sizes):
        return {"n_batches": 0, "min_cells": 0,
                "median_cells": 0.0, "max_cells": 0}
    return {"n_batches": int(len(sizes)), "min_cells": int(sizes.min()),
            "median_cells": float(sizes.median()), "max_cells": int(sizes.max())}


def _stratified_indices(batch: pd.Series, n_keep: int, seed: int) -> np.ndarray:
    """Seeded sample that preserves every batch group when the budget allows.

    Reserve up to ``MIN_SAMPLED_PER_BATCH`` cells per group, then fill the
    remaining budget uniformly from all unselected cells.  This prevents a
    valid rare barcode-derived batch from disappearing by chance while still
    keeping the bulk of the sample representative of the full dataset.
    """
    n = len(batch)
    if n_keep >= n:
        return np.arange(n)
    rng = np.random.default_rng(seed)
    groups = []
    arr = batch.astype("string").to_numpy()
    for value in sorted(pd.unique(arr), key=str):
        groups.append(np.flatnonzero(arr == value))
    reserve = min(MIN_SAMPLED_PER_BATCH, n_keep // len(groups))
    selected = []
    if reserve:
        for idx in groups:
            selected.extend(rng.choice(idx, min(reserve, len(idx)),
                                       replace=False).tolist())
    selected_arr = np.asarray(selected, dtype=int)
    left = n_keep - len(selected_arr)
    if left:
        available = np.setdiff1d(np.arange(n), selected_arr, assume_unique=False)
        selected_arr = np.concatenate(
            [selected_arr, rng.choice(available, left, replace=False)])
    return np.sort(selected_arr)


def _lisi_median(embedding: np.ndarray, labels: pd.Series) -> float:
    from harmonypy import compute_lisi

    meta = pd.DataFrame({"label": labels.to_numpy()})
    perp = min(LISI_PERPLEXITY, max(5, embedding.shape[0] // 4))
    return float(np.median(compute_lisi(embedding, meta, ["label"], perp)))


def _norm_ilisi(v: float, n: int) -> float:
    return round((v - 1) / (n - 1), 4) if n > 1 else 0.0


def _norm_clisi(v: float, n: int) -> float:
    return round((n - v) / (n - 1), 4) if n > 1 else 1.0


def _pc_regression_r2(pcs: np.ndarray, var_ratio: np.ndarray, batch: pd.Series) -> float:
    """Variance-ratio-weighted one-way ANOVA R^2 of batch on each PC."""
    codes = pd.Categorical(batch.to_numpy()).codes
    r2 = []
    for i in range(pcs.shape[1]):
        x = pcs[:, i]
        gm = x.mean()
        ss_tot = ((x - gm) ** 2).sum()
        if ss_tot <= 0:
            r2.append(0.0)
            continue
        df = pd.DataFrame({"x": x, "g": codes}).groupby("g")["x"]
        ss_between = (df.count() * (df.mean() - gm) ** 2).sum()
        r2.append(float(ss_between / ss_tot))
    w = var_ratio[: pcs.shape[1]]
    return round(float((np.asarray(r2) * w).sum() / w.sum()), 4)


def _run_harmony(pcs: np.ndarray, batch: pd.Series, seed: int):
    """(embedding, converged, note) — crash/non-convergence is an observation."""
    import harmonypy

    records: list[str] = []

    class _Capture(logging.Handler):
        def emit(self, rec):  # noqa: D401
            records.append(rec.getMessage())

    hlog = logging.getLogger("harmonypy")
    cap = _Capture()
    hlog.addHandler(cap)
    try:
        meta = pd.DataFrame({"batch": batch.to_numpy()})
        try:
            ho = harmonypy.run_harmony(pcs, meta, ["batch"],
                                       max_iter_harmony=20, random_state=seed)
        except TypeError:  # older harmonypy without random_state
            ho = harmonypy.run_harmony(pcs, meta, ["batch"], max_iter_harmony=20)
        Z = ho.Z_corr
        if hasattr(Z, "cpu"):  # harmonypy >=0.2 may hand back a torch tensor
            Z = Z.cpu().numpy()
        Z = np.asarray(Z)
        # classic harmonypy returns (n_pcs, n_cells); torch builds (n_cells, n_pcs)
        post = Z if Z.shape[0] == len(batch) else Z.T
        if post.shape[0] != len(batch):
            raise ValueError(f"unexpected Z_corr shape {Z.shape}")
        converged = any("Converged" in m for m in records)
        note = "" if converged else "harmony stopped without reporting convergence"
        return post, converged, note
    except Exception as exc:  # noqa: BLE001 - the crash IS the signal
        return None, False, f"harmony raised {type(exc).__name__}: {exc}"
    finally:
        hlog.removeHandler(cap)


def _pseudo_labels(adata, seed: int) -> pd.Series:
    """Leiden communities on the PRE-integration kNN graph (spec §5)."""
    import scanpy as sc

    # At the probe cap (5k cells), exact sklearn neighbors plus Gaussian graph
    # weights avoid the per-process Numba compilation used by UMAP weights.
    sc.pp.neighbors(adata, use_rep="X_pca", random_state=seed,
                    transformer="sklearn", method="gauss")
    try:
        sc.tl.leiden(adata, resolution=LEIDEN_RESOLUTION, random_state=seed,
                     key_added="_pseudo_ct", flavor="leidenalg")
    except ImportError as exc:
        raise _Stop(EXIT_ERROR, [f"leiden unavailable: {exc} — install the "
                                 f"[probe] extra (leidenalg)"])
    return adata.obs["_pseudo_ct"].astype("string")


def _run(args, res: dict) -> int:
    import anndata as ad
    import scanpy as sc

    timings = res["metrics"].setdefault("timings", {})
    t0 = time.perf_counter()

    adata = ad.read_h5ad(args.src)
    if adata.n_obs < MIN_CELLS:
        raise _Stop(EXIT_REJECTED, [f"only {adata.n_obs} cells (< {MIN_CELLS}) "
                                    f"— too small to probe"])
    try:
        batch_full, _ = resolve_spec(adata, args.batch_col)
    except ColumnSpecError as exc:
        raise _Stop(EXIT_REJECTED, [str(exc)])
    full_stats = _group_stats(batch_full)
    batch_missing = batch_full.isna()
    full_stats["n_batch_missing"] = int(batch_missing.sum())
    full_stats["batch_missing_frac"] = round(float(batch_missing.mean()), 4)
    res["metrics"].update(full_stats)
    if full_stats["n_batches"] < 2:
        raise _Stop(EXIT_REJECTED, [f"batch column {args.batch_col!r} has a "
                                    f"single group — nothing to correct"])
    ct_full, ct_kind = (None, "pseudo")
    if args.cell_type_col:
        try:
            ct_full, ct_kind = resolve_spec(adata, args.cell_type_col)[0], "annotated"
        except ColumnSpecError as exc:
            raise _Stop(EXIT_REJECTED, [str(exc)])
    if batch_missing.any():
        adata = adata[~batch_missing.to_numpy()].copy()
        batch_full = batch_full.loc[~batch_missing]
        if ct_full is not None:
            ct_full = ct_full.loc[~batch_missing]
        res["reasons"].append(
            f"excluded {int(batch_missing.sum())} cells with missing batch labels from the probe")
    if adata.n_obs < MIN_CELLS:
        raise _Stop(EXIT_REJECTED, [
            f"only {adata.n_obs} cells have known batch labels (< {MIN_CELLS})"])
    timings["load"] = round(time.perf_counter() - t0, 3); t0 = time.perf_counter()

    # One seeded, group-preserving subsample for this candidate.
    n_keep = min(args.n_cells, adata.n_obs)
    idx = _stratified_indices(batch_full, n_keep, args.seed)
    A = adata[idx].copy()
    batch = batch_full.iloc[idx]
    n_b = int(batch.nunique())
    res["metrics"]["n_cells_sampled"] = n_keep
    res["metrics"]["n_batches_sampled"] = n_b
    res["metrics"]["sampling"] = "stratified_min_per_batch"
    if n_b < 2:
        raise _Stop(EXIT_REJECTED, ["subsample retained a single batch — "
                                    "raise --n-cells"])

    sc.pp.highly_variable_genes(A, n_top_genes=min(args.n_hvg, A.n_vars - 1))
    A = A[:, A.var.highly_variable].copy()
    sc.pp.scale(A, max_value=10)
    n_comps = min(N_PCS, A.n_obs - 1, A.n_vars - 1)
    sc.tl.pca(A, n_comps=n_comps, svd_solver="arpack", random_state=args.seed)
    pre = A.obsm["X_pca"]
    timings["hvg_pca"] = round(time.perf_counter() - t0, 3); t0 = time.perf_counter()

    ct = (ct_full.iloc[idx].fillna("<missing>")
          if ct_full is not None else _pseudo_labels(A, args.seed))
    n_t = int(ct.nunique())
    timings["labels"] = round(time.perf_counter() - t0, 3); t0 = time.perf_counter()

    post, converged, note = _run_harmony(pre, batch, args.seed)
    timings["harmony"] = round(time.perf_counter() - t0, 3); t0 = time.perf_counter()

    m = res["metrics"]
    m["harmony_converged"] = converged
    if note:
        res["reasons"].append(note)
    m["pc_regression_r2"] = _pc_regression_r2(
        pre, A.uns["pca"]["variance_ratio"], batch)
    ilisi_pre = _lisi_median(pre, batch)
    clisi_pre = _lisi_median(pre, ct)
    m.update({"ilisi_pre": round(ilisi_pre, 3),
              "ilisi_norm_pre": _norm_ilisi(ilisi_pre, n_b),
              "clisi_norm_pre": _norm_clisi(clisi_pre, n_t),
              "clisi_labels": ct_kind, "n_cell_types": n_t})
    if ct_kind == "pseudo":
        m["pseudo_label_graph"] = PSEUDO_LABEL_GRAPH
    if post is not None:
        ilisi_post = _lisi_median(post, batch)
        clisi_post = _lisi_median(post, ct)
        m.update({"ilisi_post": round(ilisi_post, 3),
                  "ilisi_norm_post": _norm_ilisi(ilisi_post, n_b),
                  "clisi_norm_post": _norm_clisi(clisi_post, n_t)})
        timings["lisi"] = round(time.perf_counter() - t0, 3)
    else:
        m.update({"ilisi_post": None, "ilisi_norm_post": None,
                  "clisi_norm_post": None})
        timings["lisi"] = round(time.perf_counter() - t0, 3)

    res["status"] = "ok"
    return EXIT_OK


def main(argv=None) -> int:
    args = build_parser().parse_args(argv)
    logging.basicConfig(level=logging.INFO, stream=sys.stderr,
                        format="%(levelname)s %(name)s: %(message)s")
    params = {"batch_col": args.batch_col, "cell_type_col": args.cell_type_col,
              "n_cells": args.n_cells, "n_hvg": args.n_hvg, "seed": args.seed}
    res = new_result("integration_probe", os.path.abspath(args.src), params)

    t0 = time.perf_counter()
    try:
        code = _run(args, res)
    except _Stop as stop:
        res["status"] = "rejected"
        res["rejected_at"] = "input"
        res["reasons"].extend(stop.reasons)
        code = stop.exit_code
        log.warning("rejected: %s", "; ".join(stop.reasons))
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
