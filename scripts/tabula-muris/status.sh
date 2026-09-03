#!/bin/bash
# One line per sample: exit codes + key verdicts pulled from result.json.
#   bash scripts/tabula-muris/status.sh        # all datasets under $SCRATCH/eca-pp-runs
set -uo pipefail
python3 - "$SCRATCH/eca-pp-runs" <<'PY'
import glob, json, os, re, sys
base = sys.argv[1]
hdr = ("ds", "organ", "std", "idc", "cells", "species", "counts", "batch", "corr", "cell_type", "probes", "cost$", "wall_s")
rows = []
for dsdir in sorted(glob.glob(f"{base}/*/")):
    ds = os.path.basename(dsdir.rstrip("/")).replace("tabula-", "")
    for d in sorted(glob.glob(f"{dsdir}*/")):
        organ = os.path.basename(d.rstrip("/"))
        st = os.path.join(d, "status.txt")
        if not os.path.isfile(st):
            continue
        txt = open(st).read()
        g = lambda k: (re.findall(rf"^{k} exit=(\d+)$", txt, re.M) or ["-"])[-1]
        r1 = r2 = None
        try: r1 = json.load(open(os.path.join(d, "standardize/result.json")))
        except Exception: pass
        try: r2 = json.load(open(os.path.join(d, "identify_columns/result.json")))
        except Exception: pass
        cells = sp = cnt = "-"
        if r1:
            cells = r1.get("metrics", {}).get("n_cells", "-")
            sp = (r1.get("species") or {}).get("code", "-")
            cnt = (r1.get("counts") or {}).get("source") or r1.get("metrics", {}).get("counts_source", "-")
        batch = corr = ct = probes = cost = wall = "-"
        if r2:
            c = r2.get("columns") or {}
            b = c.get("batch") or {}
            batch = b.get("value", "-"); corr = b.get("correction", "-")
            ct = (c.get("cell_type") or {}).get("value", "-")
            probes = len(r2.get("trials") or [])
            llm = r2.get("metrics", {}).get("llm") or {}
            cost = f"{llm.get('cost_usd', 0):.2f}"
            wall = f"{r2.get('metrics', {}).get('timings', {}).get('total', 0):.0f}"
        rows.append((ds, organ, g("standardize"), g("identify-columns"), cells, sp, cnt,
                     batch, corr, ct, probes, cost, wall))
w = [max(len(str(x)) for x in col) for col in zip(hdr, *rows)] if rows else [len(h) for h in hdr]
for r in (hdr, *rows):
    print("  ".join(str(x).ljust(n) for x, n in zip(r, w)).rstrip())
if rows:
    tot = sum(float(r[11]) for r in rows if r[11] != "-")
    print(f"\n{len(rows)} samples, agent cost total ${tot:.2f}")
PY
