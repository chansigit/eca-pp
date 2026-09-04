#!/bin/bash
# Sherlock-only environment bootstrap for eca-pp — ALL cluster-specific fixups
# live here, none in the Python code. On any other machine, skip this file and
# `pip install -e .` into a normal environment instead.
#
# Fixups:
#   - unset PYTHONPATH    : drop Lmod's py3.12 numpy/h5py that shadow the venv
#   - ml load hdf5/1.14.4 : dl2025's h5py needs libhdf5.so.310 at runtime
#   - dl2025 venv python  : prebuilt shared env (anndata, scipy, stancounts, pytest)
#
# Agent harness note: HARNESS=deepseek is the default and uses the source-built
# dsh CLI at $SCRATCH/tools/deepseek-harness-src/apps/cli/lib/bin.js (override
# with DSH_BIN). Set HARNESS=openai for the OpenAI Agents SDK + Doubao comparison,
# or HARNESS=claude for Claude Agent SDK; ECA_PP_CLAUDE_CLI may point at an
# npm-installed Claude CLI on old glibc hosts.
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
# Agent SDK initialize handshake: the npm `claude` CLI cold-starts slowly on a
# compute node (node + NFS-backed ~/.claude + plugins), so allow 3 min instead
# of the SDK's 60 s default (value in ms; eca_pp.agent also retries transients).
export CLAUDE_CODE_STREAM_CLOSE_TIMEOUT="${CLAUDE_CODE_STREAM_CLOSE_TIMEOUT:-180000}"

cmd="${1:-}"
shift || true
case "$cmd" in
  standardize)       exec "$DL/bin/python" -m eca_pp.standardize "$@" ;;
  identify-columns)  exec "$DL/bin/python" -m eca_pp.identify_columns "$@" ;;
  integration-probe) exec "$DL/bin/python" -m eca_pp.probe "$@" ;;
  test)        cd "$REPO"
               exec "$DL/bin/python" -m pytest -p no:cacheprovider -o addopts="" "$@" ;;
  python)      exec "$DL/bin/python" "$@" ;;
  *) echo "usage: bash run.sh {standardize|identify-columns|integration-probe|test|python} [args...]" >&2; exit 64 ;;
esac
