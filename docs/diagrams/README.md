# Project workflow diagram

The reader-facing workflow lives in `docs/workflow.archify.json`. It was checked
against ECA-PP source at commit `796e915c2d1680e511232482876d791744338f63`, particularly
`standardize/cli.py`, `identify_columns/cli.py`, and their result contracts.
The ECA-RSI handoff denotes the intended ecosystem workflow, not an automatic
ECA-RSI launch by either CLI.

From the repository root, with the archify skill installed:

```bash
node /path/to/archify/bin/archify.mjs validate workflow docs/workflow.archify.json --quality showcase --json
node /path/to/archify/bin/archify.mjs deliver workflow docs/workflow.archify.json docs/workflow.html --quality showcase --json
python3 docs/diagrams/export_svg.py
node /path/to/archify/bin/archify.mjs visual-check docs/workflow.html --json
```

The SVG export helper is adapted from OSP's `docs/diagrams/export_svg.py`.
It extracts the delivered SVG and resolves its theme colors for GitHub images;
the interactive HTML is kept intact. Light/dark previews intentionally contain
only the diagram; explanatory cards are available in the HTML.

## Delivery record

- Diagram type: workflow (schema v2).
- Specification SHA-256: `4e98e4501fd105915765a78da8b14500a6b4f2278a96c5ba779ce547eadf7975` (6,057 bytes).
- HTML SHA-256: `d2a0a8d444f5505bef2e4b83c942aa943d4883e231996db1d0f6e419e547f1ee` (716,309 bytes).
- Validation: 9/9 showcase checks, 0 errors, 0 warnings.
- Browser evidence: skipped; Chrome/Chromium unavailable in the authoring environment.
- Visual review: both exported SVGs inspected through CairoSVG rasterization at 1,920 pixels wide. Text and routes are readable; this does not verify HTML interactions or viewport containment.
- Visual correction rounds: 1 (clarified legend wording).

Regenerate the record when the diagram changes. Browser-check sidecars are local
diagnostics and need not be committed.
