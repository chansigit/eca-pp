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

Two steps are available:

- **`ecasteps-standardize`** turns a single-sample `.h5ad` of unknown
  provenance into a standardized form the rest of the pipeline can trust.
- **`ecasteps-identify-columns`** identifies the batch column (for
  integration) and the cell-type column among the obs columns — an
  agent-driven step that verifies its candidates by running small-scale
  integration trials (**`ecasteps-integration-probe`**, also usable
  standalone) and records every decision with its evidence.

## Quick start

```bash
pip install git+https://github.com/chansigit/stancounts \
            git+https://github.com/chansigit/stangene
pip install ".[probe,agent]"   # from a checkout of this repo

ecasteps-standardize SRC.h5ad -o out/sample1
ecasteps-identify-columns out/sample1/standardized.h5ad -o out/sample1_columns
echo $?                        # 0 = success; see the exit-code table below
```

Python ≥ 3.10. Extras: `[probe]` (scanpy/harmonypy stack, needed by
identify-columns and the probe), `[agent]` (Claude Agent SDK for
identify-columns; authenticates via `ANTHROPIC_API_KEY` or the Claude Code
CLI's stored credentials), `[llm]` (LLM species fallback for standardize),
`[test]` (pytest). On Stanford Sherlock, `bash run.sh <tool> ...` wraps the
same commands with the cluster environment set up (compute nodes only).

## What you get

| tool | outputs |
|---|---|
| standardize | `standardized.h5ad` — integer counts layer · log-normalized `X` · canonical gene symbols with full mapping provenance in `var` · authoritative QC columns in `obs` (`pct_counts_mt`, `pct_counts_hb`, `total_counts`, `n_genes_by_counts`) · run provenance in `uns` |
| identify-columns | the verdict in `result.json → columns` (batch column + whether correction is even needed, cell-type column, each with confidence and evidence) · `batch.tsv` when the batch is a derived column (barcode/composite) · one UMAP panel per trial · a full audit trail (`decisions` with the agent's per-round reasoning and tool use, `trials` with iLISI/cLISI metrics) |
| every tool | `result.json` — every decision and its evidence, per-stage wall times, **written on failure too**. All writes are atomic: a crash never leaves a torn output. |

## Exit codes — the caller's contract

`0` success · `2` permanent data problem, don't retry · `3` blocked on a
decision, re-run with a flag or decide from the evidence · `1` unexpected error.

<details>
<summary>Full table</summary>

| code | meaning | caller's action |
|---|---|---|
| 0 | success (`result.json` may carry non-blocking review notes) | use the outputs |
| 2 | permanent data problem (too few cells/genes, no recoverable counts, not an h5ad) | skip this sample; retrying cannot help |
| 3 | blocked on a decision (ambiguous counts layer, unresolvable species, undecidable batch column) | read the evidence in `result.json`; re-run with `--counts-layer` / `--species`, or settle the column choice yourself |
| 1 | unexpected error | retry / investigate |

</details>

## Options

All flags are optional; with none given, everything is inferred and defaults
apply.

<details>
<summary>ecasteps-standardize flags</summary>

| flag | default | effect |
|---|---|---|
| `--species CODE` | inferred | pin the species, skipping inference. Accepts a short code or full name: `hs`/`human`, `mm`/`mouse`, `rn`/`rat`, `dr`/`zebrafish`, `dm`/`fruit_fly`, `ce`/`c_elegans`, `cyno`/`cynomolgus`, `rhesus`, `marmoset`, `lemur`/`mouse_lemur`. An unrecognized value exits 3 with the supported list in `result.json`. |
| `--counts-layer NAME` | inferred | pin which layer holds the raw counts, skipping inference |
| `--min-cells N` / `--min-genes N` | 100 / 5000 | hard QC gates: samples below either threshold are rejected (exit 2, no h5ad) |
| `--no-gate` | off | disable both gates (e.g. for rare, deliberately small samples); metrics are still computed and recorded |
| `--keep-unmapped` | off | keep features that cannot be mapped to a canonical gene (unknown names, ambiguous old symbols, spike-ins) under their original names, instead of dropping them. Drops are counted per category in `result.json`; cells are never removed either way. |
| `--llm` | off | allow one LLM call as a species-inference fallback (needs `ANTHROPIC_API_KEY`); any failure falls back to exit 3 |

</details>

<details>
<summary>ecasteps-identify-columns flags</summary>

| flag | default | effect |
|---|---|---|
| `--max-probes N` | 6 | budget of integration trials the agent may run |
| `--n-cells N` | adaptive | probe subsample size; default `clamp(50 × max_batches, 5000, 30000)` |
| `--no-probe` | off | profile + candidate ranking only, no trials (degraded mode, exit 3) |
| `--seed N` | 0 | sampling / integration seed; same input + seed → same trial metrics |

Without Agent SDK credentials the step degrades the same way as
`--no-probe`. `ECASTEPS_CLAUDE_CLI` can point at a specific `claude`
executable. Downstream tools accept the identified batch as
`--batch-col <obs column or batch.tsv path>`.

</details>

## How it works

**standardize** — validate & fast-reject from HDF5 metadata → locate the raw
counts (recognized layers, integer X, or log1p reversal, with a deterministic
consistency proof for oddly-named layers) → resolve the species (flag →
deterministic inference → optional LLM → stop and ask) → harmonize gene names
(unmappable features dropped by default) → compute species-aware QC columns →
write atomically.

**identify-columns** — profile every obs column (stats + sampled values,
group-size health, a nesting/equivalence graph) and enumerate derived
candidates (barcode prefixes, composites) → classify candidates
(technical / donor / condition are eligible batches; annotation, QC and
identifier columns never are) → an agent picks candidates bottom-up and
verifies each with a small-scale Harmony trial (iLISI mixing gain, cLISI
structure preservation, convergence, UMAP panel) → concludes one of four
verdicts: batch + correction recommended, batch + correction unnecessary,
no batch structure, or undecidable (exit 3). Every round is recorded.

## Docs

- [`docs/tutorial.md`](docs/tutorial.md) — hands-on walkthrough on a real
  dataset (Tabula Muris droplet, 12 organs), with per-stage timings. Chinese.
- [`docs/standardize-spec.md`](docs/standardize-spec.md) — standardize
  requirements and design.
- [`docs/identify-columns-spec.md`](docs/identify-columns-spec.md) —
  identify-columns requirements and design, including the agent doctrine
  and the full prompt.
