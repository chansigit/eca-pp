# ecasteps

Mining scRNA-seq datasets from the published literature is slow, manual work.
Public datasets are pervasively non-standard: raw counts hide in arbitrarily
named layers or arrive already log-transformed, the species goes unstated,
gene identifiers mix outdated and ambiguous symbols with spike-ins, and
metadata follows no convention. Curating datasets at atlas scale means
resolving the same ambiguities over and over by hand.

**ecasteps** automates this curation for the Open Cell Atlas pipeline. Each
step is a command-line tool — **h5ad in → files out + `result.json`** — that
resolves everything it can deterministically, records the evidence behind
every decision, and escalates only the genuinely ambiguous cases (exit
code 3, optionally with LLM assistance) for a human or an agent to settle.
Steps run offline by default and report through exit codes, so they compose
cleanly into scripts and larger workflows.

The first step, **`ecasteps-standardize`**, turns a single-sample `.h5ad` of
unknown provenance into a standardized form the rest of the pipeline can trust.

## Quick start

```bash
# stancounts and stangene are not on PyPI yet — install from source checkouts
pip install /path/to/stancounts /path/to/stangene
pip install .

ecasteps-standardize SRC.h5ad -o out/sample1
echo $?          # 0 = success; see the exit-code table below
```

On Stanford Sherlock, `bash run.sh standardize ...` wraps the same command
with the cluster environment set up (compute nodes only). Python ≥ 3.10.
Optional extras: `.[llm]` (LLM species fallback), `.[test]` (pytest).

## What you get

| file | contents |
|---|---|
| `OUTDIR/standardized.h5ad` | `layers["counts"]` (integer raw counts) · `X` = log-normalized counts (float32) · canonical gene symbols in `var_names` with the original names and full mapping provenance in `var` · per-cell QC columns in `obs` (`pct_counts_mt`, `pct_counts_hb`, `total_counts`, `n_genes_by_counts`) · run provenance in `uns` |
| `OUTDIR/result.json` | every decision and its evidence: status + reasons, species evidence, a per-layer census, gene-drop statistics (`genes_kept`/`genes_dropped` — features, never cells), QC gene-set hit counts, and wall-time per stage. **Written on failure too** — reviewing a run never requires reopening the h5ad. |

Both files are written atomically: a crash never leaves a torn output.

## Exit codes — the caller's contract

`0` success · `2` permanent data problem, don't retry · `3` blocked on a
decision, re-run with a flag · `1` unexpected error.

<details>
<summary>Full table</summary>

| code | meaning | caller's action |
|---|---|---|
| 0 | success (`result.json` may carry non-blocking review notes) | use the outputs |
| 2 | permanent data problem (too few cells/genes, no recoverable counts, not an h5ad) | skip this sample; retrying cannot help |
| 3 | blocked on a decision (ambiguous counts layer, unresolvable species) | read the evidence in `result.json`, re-run with `--counts-layer` / `--species` |
| 1 | unexpected error | retry / investigate |

</details>

## Options

All flags are optional; with none given, everything is inferred and default
gates apply.

<details>
<summary>Flag reference</summary>

| flag | default | effect |
|---|---|---|
| `--species CODE` | inferred | pin the species, skipping inference. Accepts a short code or full name: `hs`/`human`, `mm`/`mouse`, `rn`/`rat`, `dr`/`zebrafish`, `dm`/`fruit_fly`, `ce`/`c_elegans`, `cyno`/`cynomolgus`, `rhesus`, `marmoset`, `lemur`/`mouse_lemur`. An unrecognized value exits 3 with the supported list in `result.json`. |
| `--counts-layer NAME` | inferred | pin which layer holds the raw counts, skipping inference |
| `--min-cells N` / `--min-genes N` | 100 / 5000 | hard QC gates: samples below either threshold are rejected (exit 2, no h5ad) |
| `--no-gate` | off | disable both gates (e.g. for rare, deliberately small samples); metrics are still computed and recorded |
| `--keep-unmapped` | off | keep features that cannot be mapped to a canonical gene (unknown names, ambiguous old symbols, spike-ins) under their original names, instead of dropping them. Drops are counted per category in `result.json`; cells are never removed either way. |
| `--llm` | off | allow one LLM call as a species-inference fallback (needs `ANTHROPIC_API_KEY`); any failure falls back to exit 3 |

</details>

## How it works

1. **Validate & fast-reject** — input checks and the cell-count gate run on
   HDF5 metadata, before the matrix is loaded.
2. **Locate the raw counts** — recognized layer names, integer X, or log1p
   reversal; an oddly-named candidate layer is *proven* to be the true counts
   by a deterministic consistency check against X.
3. **Resolve the species** — explicit flag → deterministic inference
   (stable-ID prefixes, mitochondrial naming styles, symbol overlap with
   bundled references) → optional single LLM call → otherwise stop and ask.
4. **Harmonize gene names** to canonical symbols (via stangene); unmappable
   features are dropped by default.
5. **Compute QC columns** with species-aware mito/hemoglobin gene sets;
   same-named columns brought by the data survive under a `__original` suffix.
6. **Write atomically** — the standard-form h5ad plus `result.json`.

## Testing

```bash
pip install .[test] && pytest tests -q       # any machine
bash run.sh test tests -q                    # Sherlock compute node
```

The same suite also passes inside a stock `python:3.12-slim` container with
everything pip-installed from local sources — on Sherlock via Apptainer:

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
