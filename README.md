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
Standardization runs offline unless `--llm` is enabled. Column identification
uses the selected model when credentials are available and otherwise follows
a deterministic policy. Exit codes let both steps compose into larger workflows.

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

Harness and model are independent choices. For example, both of these use the
OpenAI Agents SDK, while selecting different Doubao models:

```bash
HARNESS=openai eca-pp-identify-columns standardized.h5ad -o out-turbo \
  --model doubao-seed-2-1-turbo-260628
HARNESS=openai eca-pp-identify-columns standardized.h5ad -o out-pro \
  --model doubao-seed-2-1-pro-260628
```

The same `--model` switch also works with `HARNESS=deepseek`; use
`ECA_PP_AGENT_MODEL` to select a model for a whole campaign.

Python ≥ 3.10. Extras: `[probe]` (scanpy/harmonypy stack), `[agent]`
(default OpenAI Agents SDK plus the DSH fallback), `[llm]` (legacy DSH alias),
`[claude]` (optional Claude Agent SDK fallback), `[openai]` (OpenAI-only
lightweight agent extra), and `[test]` (pytest). Agent calls use
`HARNESS=openai` by default, driving Doubao through Volcengine Ark. Set
`HARNESS=deepseek` to use DSH, or `HARNESS=claude` to use Claude.
Models default to `doubao-seed-2-1-turbo-260628` for both Doubao backends and
`claude-sonnet-5` for Claude, and can be overridden with `--model` /
`ECA_PP_AGENT_MODEL`. On Stanford
Sherlock, `bash run.sh <tool> ...` wraps the same commands with the cluster
environment set up (compute nodes only).

## What you get

| tool | outputs |
|---|---|
| standardize | `standardized.h5ad` — integer counts layer · log-normalized `X` · canonical gene symbols with full mapping provenance in `var` · authoritative QC columns in `obs` (`pct_counts_mt`, `pct_counts_hb`, `total_counts`, `n_genes_by_counts`) · run provenance in `uns` |
| identify-columns | the verdict in `result.json → columns` (batch column + whether correction is even needed, cell-type column, each with confidence and evidence) · `batch.tsv` when the batch is a derived column (barcode/composite) · a full audit trail (`decisions` with the agent's per-round reasoning, tool use, and token/cost usage — totals in `metrics.llm` with a billing URL, `trials` with iLISI/cLISI metrics) |
| every tool | `result.json` — every decision and its evidence, per-stage wall times, **written on failure too** when the destination is writable. Formal result files are published atomically; runtime logs are incremental. |

Reusing an output directory moves known previous step outputs into
`.history/<step>-<unique suffix>/` before running. This includes `result.json`,
`standardized.h5ad` for standardize, and `batch.tsv`, `candidates/`, and `trial_N/`
for identify-columns. Other files are retained. Use a different output directory
if the input overlaps those previous outputs. Read the current `result.json`
before consuming output files; publication of multiple files is not one transaction.

Existing QC and gene-mapping metadata are preserved in `name__original[_N]`
backup columns. Equal backups within that field's backup family share the shortest
existing name; distinct values use the next free name. Unrelated author columns
are retained even when they contain equal values.

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

The source path and `-o/--outdir` are required. Other options use the defaults below.

<details>
<summary>eca-pp-standardize flags</summary>

| flag | default | effect |
|---|---|---|
| `--species CODE` | inferred | pin the species, skipping inference. Accepts a short code or full name: `hs`/`human`, `mm`/`mouse`, `rn`/`rat`, `dr`/`zebrafish`, `dm`/`fruit_fly`, `ce`/`c_elegans`, `cyno`/`cynomolgus`, `rhesus`, `marmoset`, `lemur`/`mouse_lemur`. An unrecognized value exits 3 with the supported list in `result.json`. |
| `--counts-layer NAME` | inferred | pin which layer holds the raw counts, skipping inference |
| `--min-cells N` / `--min-genes N` | 100 / 5000 | hard QC gates: samples below either threshold are rejected (exit 2, no h5ad) |
| `--no-gate` | off | disable both gates (e.g. for rare, deliberately small samples); metrics are still computed and recorded |
| `--keep-unmapped` | off | keep features that cannot be mapped to a canonical gene (unknown names, ambiguous old symbols, spike-ins) under their original names, instead of dropping them. Drops are counted per category in `result.json`; cells are never removed either way. |
| `--llm` | off | allow a harness session with a validated `submit_species` tool when deterministic inference is inconclusive; uses `HARNESS` and the model environment settings. Unresolved inference exits 3 and logs the failure cause. |

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
`warnings`. Both Doubao backends use `ARK_API_KEY`; `HARNESS=openai` additionally
needs `pip install ".[openai]"`. Its Responses backend defaults to minimal reasoning,
server-side response chaining, serialized function calls, and up to two continuation
nudges when the model ends without submitting. Override these with
`OPENAI_AGENTS_REASONING_EFFORT`, `OPENAI_AGENTS_SERVER_STATE=0`, and
`OPENAI_AGENTS_MAX_NUDGES`. DSH's `DSH_BIN` may point at a source-built dsh
CLI (`apps/cli/lib/bin.js`) and otherwise defaults to
`$SCRATCH/tools/deepseek-harness-src/apps/cli/lib/bin.js`. For the fallback,
`ECA_PP_CLAUDE_CLI` can point at a specific `claude` executable.
Downstream tools accept the identified batch as
`--batch-col <obs column or batch.tsv path>`.

All three backends receive the structured decision state in the prompt and expose
only the validated submit tool for this step. Backend failures trigger the local
deterministic policy; they do not automatically switch harness or model. A probe
technical failure (exit 1 or an unexpected exit) makes identify-columns fail with
exit 1. Only probe exit 2 is treated as candidate rejection.

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
annotations are reported as cell type; equivalent existing/derived partitions
share one canonical probe) → an agent picks candidates bottom-up and
verifies each with a small-scale Harmony trial (iLISI mixing gain, cLISI
structure preservation on non-missing annotation cells, convergence); clear primary-candidate results finish
through a guarded metric fast path, while borderline and fallback results return
to the agent → concludes one of four
verdicts: batch + correction recommended, batch + correction unnecessary,
no batch structure, or insufficient evidence. The latter two return batch null
with structured evidence. Missing cell type is likewise a valid
null result. Every round is recorded.

## Docs

- [`docs/architecture.html`](docs/architecture.html) — interactive architecture;
  its editable source is `docs/architecture.archify.json`.
- [`docs/tutorial.md`](docs/tutorial.md) — hands-on walkthrough on a real
  dataset (Tabula Muris droplet, 12 organs), with per-stage timings. Chinese.
- [`docs/standardize-spec.md`](docs/standardize-spec.md) — standardize
  requirements and design.
- [`docs/identify-columns-spec.md`](docs/identify-columns-spec.md) —
  identify-columns requirements and design, including the agent doctrine
  and the full prompt.
