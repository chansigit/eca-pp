# identify-columns agent benchmark

A replay benchmark for the `identify-columns` classifier: feed a model the exact evidence a
past production run showed it, and score its answer against a curated gold set.

## Why replay works

`cli.build_evidence(profile, candidates)` is a pure function, and every finished run stores both
`profile` and `candidates` in its `result.json`. So a case can be replayed with **one LLM call** —
no h5ad, no scanpy, no Harmony, no GPU. A full sweep over the gold set costs minutes, which is what
makes "is model X good enough?" a question you can actually answer per model.

## Two tracks — do not merge them

The failure modes on record split into two kinds, and conflating them produces the wrong
conclusion ("we paid for a better model and the score did not move").

**Track 1 — classification (model-sensitive).** Does the model label the columns correctly?
Cluster IDs are not cell types; a cell-cycle phase is not a batch candidate; an aging atlas's
`Age_group` is a `condition`, not a technical factor. Swapping `HARNESS`/`--model` changes this
score. This is the model benchmark.

**Track 2 — adoption logic (model-independent).** Given a *correct* classification plus the trial
metrics, does the pipeline reach the right verdict? Issue #1 lives entirely here: the classifier
labelled `Genotype` as `class: "condition"` and said in its own reason that it was only a fallback —
and the adoption path took it anyway, because nothing consults that field. A better model cannot fix
this; only code can. Track 2 needs no LLM at all and belongs in the normal test suite.

**Scoring consequence:** for the mouse-pansci cases the gold answer is *not* "the model should output
batch=null". The model is judged on whether it classes `batch` as `technical` and `Age_group` /
`Genotype` as `condition`. Whether those get adopted is Track 2's question.

## Gold set

Ground truth is scarce and must stay honest. Only these sources qualify:

| source | n | what it gives |
| --- | --- | --- |
| mouse-pansci, 4 organs with `result.json.orig` | 4 | machine verdict vs human correction, same dataset |
| mouse-pansci, the other 8 organs | 8 | same obs schema, model got it right natively — controls |
| characterised bugs from git log / memory | 4 | CellLines `cell_cycle_phase`, abm-ilcp `ann0608`, TM-drop `cell_type=null`, TM-FACS empty-string artefact |

**Not gold:** the other ~271 finished runs are unreviewed machine output. Scoring against them
measures agreement with a past model, not correctness. The 105-dataset 3CA campaign is worse than
neutral — it was produced by v0.4 and memory records it "likely carries similar errors".

Growing the set means a human labelling pass, not more scraping.

## Usage

```bash
python evals/build_gold.py                 # refresh evals/gold/cases.json from Oak
HARNESS=openai python evals/run_eval.py    # score the current default model
HARNESS=claude ECA_PP_AGENT_MODEL=claude-sonnet-5 python evals/run_eval.py --tag claude-sonnet-5
python evals/run_eval.py --compare runs/   # table of every model that has been scored
```

Each run writes `runs/<tag>.json` with per-case answers, so a regression can be diffed rather than
re-argued.
