#!/bin/bash
# Sherlock-only environment bootstrap for ecasteps — ALL cluster-specific fixups
# live here, none in the Python code. On any other machine, skip this file and
# `pip install -e .` into a normal environment instead.
#
# Fixups:
#   - unset PYTHONPATH    : drop Lmod's py3.12 numpy/h5py that shadow the venv
#   - ml load hdf5/1.14.4 : dl2025's h5py needs libhdf5.so.310 at runtime
#   - dl2025 venv python  : prebuilt shared env (anndata, scipy, stancounts, pytest)
#
# Agent SDK note: claude-agent-sdk's BUNDLED claude binary needs glibc >= 2.25
# (CentOS 7 has 2.17). identify-columns therefore prefers the npm-installed
# `claude` on PATH (newer anyway), overridable via ECASTEPS_CLAUDE_CLI.
# Alternative when no npm CLI exists:
#   ml load system polyfill-glibc
#   polyfill-glibc --target-glibc=2.17 <site-packages>/claude_agent_sdk/_bundled/claude
# (re-patch after every claude-agent-sdk upgrade).
#
# ALWAYS run on a compute node, never the login node:
#   bash run.sh standardize SRC.h5ad -o OUTDIR [--min-cells N ...]
#   bash run.sh identify-columns STD.h5ad -o OUTDIR [--max-probes N ...]
#   bash run.sh integration-probe STD.h5ad --batch-col COL -o OUTDIR
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
  standardize)       exec "$DL/bin/python" -m ecasteps.standardize "$@" ;;
  identify-columns)  exec "$DL/bin/python" -m ecasteps.identify_columns "$@" ;;
  integration-probe) exec "$DL/bin/python" -m ecasteps.probe "$@" ;;
  test)        cd "$REPO"
               exec "$DL/bin/python" -m pytest -p no:cacheprovider -o addopts="" "$@" ;;
  python)      exec "$DL/bin/python" "$@" ;;
  *) echo "usage: bash run.sh {standardize|identify-columns|integration-probe|test|python} [args...]" >&2; exit 64 ;;
esac
