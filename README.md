# ecasteps

Standalone curation steps for the Open Cell Atlas pipeline. Each step is a
plain CLI — **h5ad in → files out + `result.json`** — with no workflow
framework, no daemon, and no network access by default, so it can be driven by
anything: a shell loop, Snakemake, an agent, or a human.

Currently shipped: **`ecasteps-standardize`**, which turns a single-sample
`.h5ad` of unknown provenance into a standardized form the rest of the
pipeline can trust.

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
4. **Harmonizes gene names** to canonical symbols (via
   [stangene](https://github.com/chansigit)), dropping unmappable features by
   default (`--keep-unmapped` to keep them).
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

On Stanford Sherlock, use the bootstrap wrapper instead (compute nodes only):

```bash
bash run.sh standardize SRC.h5ad -o OUTDIR [...]
bash run.sh test tests -q
```

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

```bash
pip install <stancounts-src> <stangene-src> .
```

Requires Python ≥ 3.10. The `stancounts` and `stangene` dependencies are
installed from their source directories for now (not yet on PyPI). Optional
extras: `.[llm]` for the LLM species fallback, `.[test]` for pytest.

## Container

No cluster dependency: the full test suite passes inside a stock
`python:3.12-slim` image with everything pip-installed from local sources.
On Sherlock this is exercised via Apptainer:

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
  result.py        result.json schema + exit codes (shared by future steps)
  atomic_io.py     atomic writes
tests/             32 acceptance tests (run natively or in the container)
run.sh             Sherlock env bootstrap (the only cluster-specific file)
```

## Docs

- [`docs/tutorial.md`](docs/tutorial.md) — hands-on walkthrough on a real
  dataset (Tabula Muris droplet, 12 organs), with per-stage timings. Chinese.
- [`docs/standardize-spec.md`](docs/standardize-spec.md) — requirements &
  design document (the internal F1–F8 numbering and phase plan live there).
