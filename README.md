# ecasteps

Data curation tools for the Open Cell Atlas pipeline. Each step is a
command-line tool — **h5ad in → files out + `result.json`** — that runs
offline by default and reports through exit codes, so steps compose cleanly
into scripts and larger workflows.

**`ecasteps-standardize`** turns a single-sample `.h5ad` of unknown
provenance into a standardized form the rest of the pipeline can trust.

## What standardize does

1. **Validates** the input and fast-rejects samples below hard QC thresholds
   (minimum cells / detected genes — checked from HDF5 metadata before the
   matrix is even loaded).
2. **Locates the raw counts**: recognized layer names, integer X, or log1p
   reversal — and can *prove* an oddly-named candidate layer is the true
   counts via a deterministic consistency check against X.
3. **Resolves the species**: explicit `--species` flag → deterministic
   inference (stable-ID prefixes, mitochondrial naming styles, symbol overlap
   with bundled references) → optional single LLM call (`--llm`) → otherwise
   stops and asks the caller to decide.
4. **Harmonizes gene names** to canonical symbols (via stangene), dropping
   unmappable features by default (`--keep-unmapped` to keep them).
5. **Computes per-cell QC columns** — `pct_counts_mt`, `pct_counts_hb`,
   `total_counts`, `n_genes_by_counts` — with species-aware mito/hemoglobin
   gene sets. Columns of the same name brought by the data are preserved under
   a `__original` suffix and overwritten.
6. **Writes the outputs atomically**: `standardized.h5ad` (integer counts
   layer + log-normalized X + full provenance in `var`/`uns`) and
   `result.json` recording every decision, per-stage timings included —
   written on failure too.

## Usage

```bash
ecasteps-standardize SRC.h5ad -o OUTDIR \
    [--species hs] [--llm] \
    [--min-cells 100] [--min-genes 5000] \
    [--counts-layer NAME] [--no-gate] [--keep-unmapped]
```

`--species` accepts a short code or full name: `hs`/`human`, `mm`/`mouse`,
`rn`/`rat`, `dr`/`zebrafish`, `dm`/`fruit_fly`, `ce`/`c_elegans`,
`cyno`/`cynomolgus`, `rhesus`, `marmoset`, `lemur`/`mouse_lemur`.
An unrecognized value exits 3 with the supported list in `result.json`.

Samples below `--min-cells` / `--min-genes` (cells; detected genes) are
rejected with exit 2 and no h5ad is produced. `--no-gate` disables both
thresholds — useful for rare, deliberately small samples — while the
metrics are still computed and recorded in `result.json`.

By default, features that cannot be mapped to a canonical gene — unknown
names, ambiguous old symbols, and non-gene features such as ERCC spike-ins —
are removed from the matrix (counted per category under
`harmonization.genes_dropped` in `result.json`; cells are never removed).
`--keep-unmapped` keeps them instead, under their original names.

On Stanford Sherlock, `bash run.sh standardize ...` wraps the same command
with the cluster environment set up (compute nodes only).

## Exit codes — the caller's contract

| code | meaning | caller's action |
|---|---|---|
| 0 | success (`result.json` may carry non-blocking review notes) | use the outputs |
| 2 | permanent data problem (too few cells/genes, no recoverable counts, not an h5ad) | skip this sample; retrying cannot help |
| 3 | blocked on a decision (ambiguous counts layer, unresolvable species) | read the evidence in `result.json`, re-run with `--counts-layer` / `--species` |
| 1 | unexpected error | retry / investigate |

`result.json` is always written — status, reasons, species evidence, a
per-layer census, gene-drop statistics (`genes_kept` / `genes_dropped` —
features, never cells), QC gene-set hit counts, and wall-time per stage — so
reviewing a run never requires reopening the h5ad.

## Install

Requires Python ≥ 3.10. The stancounts and stangene dependencies are not on
PyPI; install them from their source checkouts first:

```bash
pip install /path/to/stancounts /path/to/stangene
pip install .
```

Optional extras: `.[llm]` for the LLM species fallback, `.[test]` for pytest.

## Container

The full test suite passes inside a stock `python:3.12-slim` image with
everything pip-installed from local sources. On Sherlock this is exercised
via Apptainer:

```bash
apptainer pull python312-slim.sif docker://python:3.12-slim
bash scripts/test-in-container.sh
```

## Layout

```
src/ecasteps/
  standardize.py   CLI + pipeline flow
  countsloc.py     counts location: layer census + consistency proof (stancounts)
  species.py       species resolution ladder
  harmonize.py     gene-name harmonization, default drop of unmappables
  qc.py            gate metrics + per-cell QC columns
  build.py         standard-form assembly + atomic h5ad write
  result.py        result.json schema + exit codes
  atomic_io.py     atomic writes
tests/             acceptance tests (run natively or in the container)
run.sh             environment bootstrap for Stanford Sherlock
```

## Docs

- [`docs/tutorial.md`](docs/tutorial.md) — hands-on walkthrough on a real
  dataset (Tabula Muris droplet, 12 organs), with per-stage timings. Chinese.
- [`docs/standardize-spec.md`](docs/standardize-spec.md) — requirements and
  design document.
