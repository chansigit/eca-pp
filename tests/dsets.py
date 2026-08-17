"""Shared synthetic-dataset builders for the acceptance tests.

Gene names/IDs come from stangene's bundled references (real approved symbols +
Ensembl IDs), so from v0.2 on the F4a species ladder resolves and F4 maps —
letting the v0.1 counts-location assertions run through the full pipeline.
"""

from __future__ import annotations

import json

import anndata as ad
import numpy as np
import pandas as pd
import scipy.sparse as sp

from stangene.references import load_reference

from ecasteps.standardize import main

RNG = np.random.default_rng(0)
N, G = 300, 6000  # default gates: min_cells=100, min_genes=5000 — this passes both

_REF_CACHE: dict = {}


def ref_genes(species: str = "human", n: int = G):
    """(symbols, ensembl_ids): n unique named+ID'd genes from the bundled
    reference, guaranteed to include the species' mito/hb genes (so F5 has
    nonzero masks and happy-path runs stay status=ok)."""
    key = (species, n)
    if key not in _REF_CACHE:
        import pandas as pd
        import stangene

        gt = load_reference(species)["gene_table"]
        rows = gt[gt["symbol"].notna() & gt["ensembl_id"].notna()]
        rows = rows.drop_duplicates("symbol").drop_duplicates("ensembl_id")
        special = (stangene.mito_mask(rows["symbol"], species)
                   | stangene.hb_mask(rows["symbol"], species))
        picked = pd.concat([rows[~special].head(n - int(special.sum())),
                            rows[special]])
        assert len(picked) == n, f"reference for {species} has < {n} usable genes"
        _REF_CACHE[key] = (list(picked["symbol"].astype(str)),
                           list(picked["ensembl_id"].astype(str)))
    return _REF_CACHE[key]


def make_counts(n: int = N, g: int = G, lam: float = 0.5) -> np.ndarray:
    """Dense integer counts with (virtually) every gene detected."""
    return RNG.poisson(lam, size=(n, g)).astype(np.float32)


def lognorm(counts: np.ndarray) -> sp.csr_matrix:
    """log1p(normalize_total(counts, 1e4)) — what scanpy would produce."""
    C = counts.astype(np.float64)
    s = C.sum(axis=1, keepdims=True)
    s[s == 0] = 1.0
    return sp.csr_matrix(np.log1p(C / s * 1e4).astype(np.float32))


def write_h5ad(path, X, layers=None, *, species: str = "human",
               var_names=None, gene_ids=None, with_ids: bool = True, obs=None):
    """Write a synthetic h5ad whose features are real reference genes.

    ``var_names``/``gene_ids`` override the reference-derived defaults;
    ``with_ids=False`` omits the gene_ids column (symbol-only dataset).
    """
    g = X.shape[1]
    if var_names is None:
        var_names, ids = ref_genes(species, g)
        if gene_ids is None:
            gene_ids = ids
    A = ad.AnnData(X=X, layers=layers or {})
    A.var_names = pd.Index(list(var_names), dtype=object)
    if with_ids and gene_ids is not None:
        A.var["gene_ids"] = list(gene_ids)
    if obs:
        for k, v in obs.items():
            A.obs[k] = v
    A.write_h5ad(path)
    return path


def run_cli(tmp_path, src, *extra) -> tuple[int, dict]:
    out = tmp_path / "out"
    code = main([str(src), "-o", str(out), *map(str, extra)])
    res = json.loads((out / "result.json").read_text())
    assert res["exit_code"] == code
    return code, res
