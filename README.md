# ecasteps

Standalone, orchestrator-free curation steps for the Open Cell Atlas pipeline —
the framework-independent successor to eca-prefect-v2's step chain. Each step is a
plain CLI: **h5ad in → files out + `result.json`**, drivable by anything (a shell
loop, Snakemake, an agent, or a human).

Spec: [`docs/standardize-spec.md`](docs/standardize-spec.md) · Status: **v0.2** —
standardize implements F1 (input validation) → F3 (hard QC gates, two-gate design)
→ F2 (counts location & recovery) → F4a (species ladder) → F4 (gene harmonization,
unmappable features dropped by default) → F5 (authoritative QC obs columns) →
F7 (atomic `standardized.h5ad`). v0.3 adds F6 (metacols, rank-only) and F8 (reports).

## Run (Sherlock)

All cluster fixups live in `run.sh`; always use a compute node:

```bash
bash run.sh standardize SRC.h5ad -o OUTDIR [--species hs] [--llm] \
    [--min-cells 100] [--min-genes 5000] \
    [--counts-layer NAME] [--no-gate] [--keep-unmapped]
bash run.sh test tests/ -q
```

Outputs: `OUTDIR/standardized.h5ad` (counts layer + lognorm X + harmonized
var_names + QC obs columns) and `OUTDIR/result.json` (always written, even on
failure).

## Run (anywhere else)

```bash
pip install -e .          # needs stancounts + stangene installed/reachable
ecasteps-standardize SRC.h5ad -o OUTDIR
```

## Container (portability proof)

No cluster dependency: the acceptance suite passes inside stock `python:3.12-slim`
with everything pip-installed from local sources (verified on Sherlock via Apptainer):

```bash
apptainer pull python312-slim.sif docker://python:3.12-slim   # unset APPTAINER_DOCKER_* first
bash scripts/test-in-container.sh                              # venv-in-container + pytest
```

## Exit codes (the driver contract)

| code | meaning | driver action |
|---|---|---|
| 0 | ok (may carry non-blocking `needs_review` flags) | next step |
| 2 | permanent data problem (QC reject — incl. post-drop, no counts, not an h5ad) | skip, never retry |
| 3 | blocked — needs a decision (counts ambiguity, unresolvable species) | read `result.json` evidence, re-run with `--counts-layer` / `--species` |
| 1 | unexpected error | retry / investigate |

`OUTDIR/result.json` carries status, reasons, metrics, and a full per-layer census
(integer-ness, sparsity, consistency-with-X) so review never requires reopening the
h5ad.

## Layout

```
src/ecasteps/
  standardize.py   CLI + flow ①–⑪ (F1 → F3 → F2 → F4a → F4 → F5 → F7)
  countsloc.py     counts location: stancounts + census + consistency proof
  species.py       F4a ladder: --species → stangene.infer_species → --llm → block
  harmonize.py     F4: stangene mapping applied; unmappable dropped by default
  qc.py            gate metrics + F5 authoritative QC obs columns
  build.py         F7: counts layer + lognorm X (scipy) + atomic h5ad write
  result.py        result.json schema + exit codes (shared by future steps)
  atomic_io.py     atomic writes
tests/             acceptance tests (spec §9); dsets.py builds real-gene datasets
run.sh             Sherlock env bootstrap (the only cluster-specific file)
```
