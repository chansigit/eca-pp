# eca-pp

**Prepare published single-cell RNA sequencing data for reuse.**

Public datasets often need considerable cleanup before they can be combined or
reanalyzed. Gene names differ between studies, raw expression counts can be hard
to find, and columns describing samples and cell types follow no common naming
convention.

eca-pp helps automate that preparation for the Open Cell Atlas pipeline. Give it
an `.h5ad` file—the data format used by Scanpy—and it produces standardized data
and a record of what it found, changed, or could not resolve.

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

The source file stays unchanged. The output retains author metadata, backs up
fields that need replacement, and records gene-mapping changes. By default,
features that cannot be mapped to a canonical gene are removed from the output;
this can be changed with `--keep-unmapped`.

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
