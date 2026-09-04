# eca-pp

Mining scRNA-seq datasets from the published literature is slow, manual work.
Public datasets are pervasively non-standard: raw counts hide in arbitrarily
named layers or arrive already log-transformed, the species goes unstated,
gene identifiers mix outdated and ambiguous symbols with spike-ins, and
metadata follows no convention. Curating datasets at atlas scale means
resolving the same ambiguities over and over by hand.

**eca-pp** automates this curation for the Open Cell Atlas pipeline. Each
step is a command-line tool — **h5ad in → files out + `result.json`** — that
resolves everything it can deterministically, records the evidence behind
every decision, and finishes unattended whenever a safe null result is
possible. Ambiguity is preserved as structured warnings instead of forcing a
risky guess.
Steps run offline by default and report through exit codes, so they compose
cleanly into scripts and larger workflows.

Two steps are available:

- **`eca-pp-standardize`** turns a single-sample `.h5ad` of unknown
  provenance into a standardized form the rest of the pipeline can trust.
- **`eca-pp-identify-columns`** identifies the batch column (for
  integration) and the cell-type column among the obs columns — an
  agent-driven step that verifies its candidates by running small-scale
  integration trials (**`eca-pp-integration-probe`**, also usable
  standalone) and records every decision with its evidence.

## Quick start

```bash
pip install git+https://github.com/chansigit/stancounts \
            git+https://github.com/chansigit/stangene
pip install ".[probe,agent]"   # from a checkout of this repo

eca-pp-standardize SRC.h5ad -o out/sample1
eca-pp-identify-columns out/sample1/standardized.h5ad -o out/sample1_columns
echo $?                        # 0 = success; see the exit-code table below
```

Python ≥ 3.10. Extras: `[probe]` (scanpy/harmonypy stack), `[agent]`
(DeepSeek Harness plus its MCP bridge; `[llm]` is an alias), `[claude]`
(optional Claude Agent SDK fallback), and `[test]` (pytest). Agent calls use
`HARNESS=deepseek` by default, driving Doubao through Volcengine Ark; set
`HARNESS=claude` to use Claude instead. Models default to
`doubao-seed-2-1-turbo-260628` / `claude-sonnet-5` according to the backend
and can be overridden with `--model` / `ECA_PP_AGENT_MODEL`. On Stanford
Sherlock, `bash run.sh <tool> ...` wraps the same commands with the cluster
environment set up (compute nodes only).

## What you get

| tool | outputs |
|---|---|
| standardize | `standardized.h5ad` — integer counts layer · log-normalized `X` · canonical gene symbols with full mapping provenance in `var` · authoritative QC columns in `obs` (`pct_counts_mt`, `pct_counts_hb`, `total_counts`, `n_genes_by_counts`) · run provenance in `uns` |
| identify-columns | the verdict in `result.json → columns` (batch column + whether correction is even needed, cell-type column, each with confidence and evidence) · `batch.tsv` when the batch is a derived column (barcode/composite) · one UMAP panel per trial · a full audit trail (`decisions` with the agent's per-round reasoning, tool use, and token/cost usage — totals in `metrics.llm` with a billing URL, `trials` with iLISI/cLISI metrics) |
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
| 3 | blocked on a decision (for example an ambiguous counts layer or unresolvable species) | read the evidence in `result.json`; re-run with an explicit override such as `--counts-layer` / `--species` |
| 1 | unexpected error | retry / investigate |

</details>

## Options

All flags are optional; with none given, everything is inferred and defaults
apply.

<details>
<summary>eca-pp-standardize flags</summary>

| flag | default | effect |
|---|---|---|
| `--species CODE` | inferred | pin the species, skipping inference. Accepts a short code or full name: `hs`/`human`, `mm`/`mouse`, `rn`/`rat`, `dr`/`zebrafish`, `dm`/`fruit_fly`, `ce`/`c_elegans`, `cyno`/`cynomolgus`, `rhesus`, `marmoset`, `lemur`/`mouse_lemur`. An unrecognized value exits 3 with the supported list in `result.json`. |
| `--counts-layer NAME` | inferred | pin which layer holds the raw counts, skipping inference |
| `--min-cells N` / `--min-genes N` | 100 / 5000 | hard QC gates: samples below either threshold are rejected (exit 2, no h5ad) |
| `--no-gate` | off | disable both gates (e.g. for rare, deliberately small samples); metrics are still computed and recorded |
| `--keep-unmapped` | off | keep features that cannot be mapped to a canonical gene (unknown names, ambiguous old symbols, spike-ins) under their original names, instead of dropping them. Drops are counted per category in `result.json`; cells are never removed either way. |
| `--llm` | off | allow one tool-less harness session as a species-inference fallback (same backend and model as identify-columns); any failure falls back to exit 3 |

</details>

<details>
<summary>eca-pp-identify-columns flags</summary>

| flag | default | effect |
|---|---|---|
| `--max-probes N` | 6 | budget of integration trials the agent may run |
| `--n-cells N` | adaptive | probe subsample size; default `clamp(50 × max_batches, 5000, 30000)` |
| `--no-probe` | off | profile + candidate ranking only; batch is safely left null with a structured warning (exit 0) |
| `--seed N` | 0 | sampling / integration seed; same input + seed → same trial metrics |
| `--model ID` | backend default | agent model; also settable via `ECA_PP_AGENT_MODEL`. The model actually used is recorded per round in `result.json` and summarized in `metrics.llm.models`. |

Without selected-backend credentials, or if the agent fails during a run, the step
continues with its deterministic policy and records that fallback in
`warnings`. DSH uses `ARK_API_KEY`; `DSH_BIN` may point at a source-built dsh
CLI (`apps/cli/lib/bin.js`) and otherwise defaults to
`$SCRATCH/tools/deepseek-harness-src/apps/cli/lib/bin.js`. For the fallback,
`ECA_PP_CLAUDE_CLI` can point at a specific `claude` executable.
Downstream tools accept the identified batch as
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
candidates (barcode prefixes, composites) → classify candidates (technical /
donor candidates are exhausted before condition or unknown fallbacks;
annotation, cluster, QC and identifier columns never are batches; only author
annotations are reported as cell type) → an agent picks candidates bottom-up and
verifies each with a small-scale Harmony trial (iLISI mixing gain, cLISI
structure preservation, convergence, UMAP panel) → concludes one of four
verdicts: batch + correction recommended, batch + correction unnecessary, or
batch null with structured evidence. Missing cell type is likewise a valid
null result. Every round is recorded.

## Docs

- [`docs/tutorial.md`](docs/tutorial.md) — hands-on walkthrough on a real
  dataset (Tabula Muris droplet, 12 organs), with per-stage timings. Chinese.
- [`docs/standardize-spec.md`](docs/standardize-spec.md) — standardize
  requirements and design.
- [`docs/identify-columns-spec.md`](docs/identify-columns-spec.md) —
  identify-columns requirements and design, including the agent doctrine
  and the full prompt.
