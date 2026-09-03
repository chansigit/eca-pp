"""Replace empty-string values in batch-type obs columns with "missing".

    python fill_missing.py FILE.h5ad [...]   # given files, in place (atomic rewrite)
    python fill_missing.py                   # every standardized.h5ad under $SCRATCH/eca-pp-runs/

Records what was changed in uns["eca_pp_standardize"]["obs_empty_filled"] and prints one
line per file. Run on a compute node (bash run.sh python scripts/tabula-muris/fill_missing.py ...).
"""
import glob, os, sys
import anndata as ad
import pandas as pd

BATCH_COLS = ("tissue", "subtissue", "channel", "mouse.id", "plate.barcode", "mouse.sex",
              "sex", "age", "method", "donor", "assay", "mouse_id", "sample_id")
FILL = "missing"

files = sys.argv[1:] or sorted(glob.glob(
    os.path.join(os.environ["SCRATCH"], "eca-pp-runs", "*", "*", "standardize", "standardized.h5ad")))

for f in files:
    tag = "/".join(f.split("/")[-4:-2])
    A = ad.read_h5ad(f)
    filled = {}
    for c in BATCH_COLS:
        if c not in A.obs:
            continue
        s = A.obs[c]
        if not (isinstance(s.dtype, pd.CategoricalDtype) or s.dtype == object or str(s.dtype) == "string"):
            continue
        mask = s.astype("string").fillna("") == ""
        n = int(mask.sum())
        if n == 0:
            continue
        if isinstance(s.dtype, pd.CategoricalDtype):
            s = s.cat.add_categories([FILL]) if FILL not in s.cat.categories else s
            s = s.where(~mask, FILL).cat.remove_unused_categories()
        else:
            s = s.astype(object).where(~mask, FILL).astype("category")
        A.obs[c] = s
        filled[c] = n
    if not filled:
        print(f"{tag:28s} no empty strings in batch columns")
        continue
    prov = A.uns.setdefault("eca_pp_standardize", {})
    prov["obs_empty_filled"] = {"value": FILL, "columns": {k: v for k, v in filled.items()}}
    tmp = f + ".tmp"
    A.write_h5ad(tmp)
    os.replace(tmp, f)
    print(f"{tag:28s} filled {filled}")
