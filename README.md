# ECA-PP

**Standardized preprocessing of real-world single-cell data for Ensemble Cell Atlas.**

ECA-PP is a component of the **Ensemble Cell Atlas (ECA)** ecosystem. It turns
heterogeneous public single-cell RNA sequencing datasets into a consistent,
documented starting point for **ECA-RSI (Recursive Self-Improvement)**, ECA's
automated data governance system. Handling routine preprocessing here reduces
the burden on ECA-RSI, leaving it to focus on further quality assessment,
annotation, and iterative refinement.

## Why public data needs more than format conversion

Mining published single-cell data means dealing with the choices made by many
different authors. Expression matrices may arrive as text tables, sparse matrix
files, Seurat objects, or H5AD files, with metadata in separate supplements.
Even after conversion, the same-looking dataset can mean very different things:

- The expression matrix may contain raw counts, normalized values, or scaled
  values; the original counts may be tucked away in another layer.
- Gene identifiers may mix Ensembl IDs, old symbols, aliases, and non-gene features.
- Sample, sequencing batch, donor, and cell-type columns may use unfamiliar names,
  duplicate one another, or contain missing values.
- Existing quality-control measurements may follow different definitions or have
  been calculated before genes were filtered out.

A file that opens successfully is not necessarily ready for scientific reuse.
Repeating these checks manually—or asking an agent to invent a new preprocessing
script for every study—makes consistent treatment across datasets difficult.

## Where ECA-PP fits

General-purpose agents can usually handle downloading files, unpacking archives,
and scripting conversions into `.h5ad`, the AnnData format used by Scanpy.
**ECA-PP starts at H5AD.** From that point onward, decisions about counts, gene
identity, QC, and batch structure need shared, domain-specific rules that are
applied consistently across studies.

| Stage | Responsibility |
|---|---|
| Data acquisition and conversion | Upstream tools or agents gather the published files and create an H5AD dataset. |
| Standardized preprocessing | ECA-PP checks the data, standardizes expression and genes, evaluates metadata, and records unresolved issues. |
| Iterative data governance | ECA-RSI uses the prepared data and evidence for further assessment and refinement. |

ECA-PP can also be used independently of ECA-RSI in your own analysis pipeline.

## What makes it useful?

- **Consistent scientific rules.** Counts validation, gene mapping, and QC follow
  a shared implementation. QC and normalized expression use the same final gene
  set, making their definitions consistent across datasets.
- **Decisions checked against the data.** Batch candidates are evaluated through
  small integration trials. Technical and donor factors take priority over
  biological conditions; model suggestions must pass programmatic validation.
- **Automation with explicit uncertainty.** Built-in rules keep column
  identification moving when a model is unavailable. Insufficient evidence can
  produce an empty selection with an explanation; essential unresolved choices
  stop preprocessing for review.
- **Traceable changes.** The source file stays unchanged. Outputs preserve author
  metadata, back up replaced fields, and record gene-mapping changes, decisions,
  and trial results so downstream systems can inspect the evidence.

## What does it do?

eca-pp currently provides two steps:

1. **Standardize the data.** Find raw counts (the original expression measurements),
   recover them from supported transformed data when possible, identify the
   species, standardize gene names, and calculate quality-control measurements.
2. **Identify useful metadata.** Find the existing cell-type annotation and a
   suitable batch column, such as a sequencing run or sample identifier. Small
   integration trials help assess whether batch correction would be useful and
   whether it would preserve cell-type structure.

For example, a study might store counts in an unfamiliar layer and describe its
cells using columns named `channel`, `mouse.id`, and `cell_ontology_class`.
eca-pp can locate the counts, produce consistent gene names and QC measurements,
and evaluate which metadata columns to use downstream. Its report explains the
selection—including when correction appears unnecessary.

By default, features that cannot be mapped to a canonical gene are removed from
the output; this can be changed with `--keep-unmapped`.

## Is it right for my data?

Use eca-pp when you have an H5AD dataset and want a consistent starting point for
further analysis. It runs on CPUs and can be used for one dataset or called from
a larger processing pipeline.

The current release prepares data and evaluates metadata. It does not assign
new biological cell-type labels, remove individual low-quality cells or doublets,
or produce a final integrated atlas. A successful run still needs to be interpreted
in the context of your study.

When evidence is insufficient, eca-pp can leave the batch or cell-type result
empty (`null`) and explain why. If a choice is required to standardize the data,
such as resolving an ambiguous species, it stops with the evidence needed to
continue.

## Try it

Use Python 3.10 or newer, preferably in a dedicated environment:

```bash
git clone https://github.com/chansigit/eca-pp.git
cd eca-pp
pip install git+https://github.com/chansigit/stancounts \
            git+https://github.com/chansigit/stangene
pip install ".[probe,openai]"
```

Start by standardizing your file:

```bash
eca-pp-standardize your-data.h5ad -o results/standardize
```

Read `results/standardize/result.json`. If standardization succeeds, use the
generated file to identify batch and cell-type columns:

```bash
eca-pp-identify-columns results/standardize/standardized.h5ad \
  -o results/columns
```

Standardization runs locally by default. For AI-assisted column identification,
set `ARK_API_KEY` in your environment before the second command. The default uses
Doubao Turbo through the OpenAI Agents SDK. Without model credentials, column
identification continues using built-in rules and integration trials.
See the [tutorial](docs/tutorial.md#8-识别批次列--细胞类型列identify-columns)
for configuration and model choices.

The default size checks require at least 100 cells and 5,000 detected genes.
For smaller datasets or other input-specific settings, consult the
[tutorial's options](docs/tutorial.md#4-常用参数) or run either command with `--help`.

## What will I get?

| Output | What it gives you |
|---|---|
| `results/standardize/standardized.h5ad` | Counts, normalized expression, standardized gene names, and QC measurements for downstream analysis. |
| `results/standardize/result.json` | Input checks, counts source, species, changes made, and any issues requiring review. |
| `results/columns/result.json` | Selected batch and cell-type columns, whether correction is recommended, and supporting evidence. |
| `results/columns/batch.tsv` | Batch labels when the selected grouping is derived from barcodes or a combination of columns. Created only when needed. |

Start with each step's `result.json`: it records the outcome and reasons.
An empty selection means no suitable column was identified; a selected batch with
`correction: "unnecessary"` means the evidence did not support correcting it.
Execution errors are reported separately. The
[tutorial](docs/tutorial.md#3-退出码--下一步动作) explains how to respond to each outcome.

## Learn more

**Using eca-pp:** the [hands-on tutorial](docs/tutorial.md) covers a real mouse
dataset, result interpretation, common options, model setup, and reruns.
The tutorial is currently in Chinese.

**Understanding or developing eca-pp:** see the
[standardization specification](docs/standardize-spec.md) and
[column-identification specification](docs/identify-columns-spec.md) for methods,
interfaces, and tests. The [interactive architecture diagram](docs/architecture.html)
can be downloaded and opened in a browser.

For questions, unexpected results, or feature requests, please
[open an issue](https://github.com/chansigit/eca-pp/issues).
