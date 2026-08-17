#!/bin/bash
# Sherlock-only environment bootstrap for ecasteps — ALL cluster-specific fixups
# live here, none in the Python code. On any other machine, skip this file and
# `pip install -e .` into a normal environment instead.
#
# Fixups (same three as eca-prefect-v2, minus .testdeps which v0.1 doesn't need):
#   - unset PYTHONPATH    : drop Lmod's py3.12 numpy/h5py that shadow the venv
#   - ml load hdf5/1.14.4 : dl2025's h5py needs libhdf5.so.310 at runtime
#   - dl2025 venv python  : prebuilt shared env (anndata, scipy, stancounts, pytest)
#
# ALWAYS run on a compute node, never the login node:
#   bash run.sh standardize SRC.h5ad -o OUTDIR [--min-cells N ...]
#   bash run.sh test [pytest args]
#   bash run.sh python script.py
set -euo pipefail

REPO="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
DL="${ECA_DL_VENV:-/scratch/users/chensj16/venvs/dl2025/.venv}"
# Local stangene source must SHADOW any older stangene installed in dl2025
# (v0.2 needs stangene>=0.5 with infer_species).
STANGENE_SRC="${STANGENE_SRC:-/home/users/chensj16/s/projects/stangene/src}"

unset PYTHONPATH || true
ml load hdf5/1.14.4 2>/dev/null || true
export PYTHONPATH="$REPO/src:$STANGENE_SRC"

cmd="${1:-}"
shift || true
case "$cmd" in
  standardize) exec "$DL/bin/python" -m ecasteps.standardize "$@" ;;
  test)        cd "$REPO"
               exec "$DL/bin/python" -m pytest -p no:cacheprovider -o addopts="" "$@" ;;
  python)      exec "$DL/bin/python" "$@" ;;
  *) echo "usage: bash run.sh {standardize|test|python} [args...]" >&2; exit 64 ;;
esac
