#!/bin/bash
# Sherlock-only environment bootstrap for eca-pp — ALL cluster-specific fixups
# live here, none in the Python code. On any other machine, skip this file and
# `pip install -e .` into a normal environment instead.
#
# Fixups:
#   - unset PYTHONPATH    : drop Lmod's py3.12 numpy/h5py that shadow the venv
#   - eca-pp-ct python    : container-only interpreter (python312-slim.sif, glibc 2.41)
#                           so official manylinux wheels work on this glibc-2.17 host.
#                           stancounts/stangene/eca-pp are editable-installed inside it,
#                           so no STANGENE_SRC shadowing and no hdf5 module are needed.
#                           Set ECA_PP_PYTHON=/path/to/python to override (e.g. a plain
#                           venv on another machine); the old dl2025 native route is
#                           ECA_PP_PYTHON=/scratch/users/chensj16/venvs/dl2025/.venv/bin/python
#                           plus `ml load hdf5/1.14.4` and STANGENE_SRC on PYTHONPATH.
#
# Agent harness note: HARNESS=openai is the default and uses the OpenAI Agents
# SDK with Doubao Turbo and medium reasoning. Set HARNESS=deepseek for the
# source-built dsh CLI at $SCRATCH/tools/deepseek-harness-src/apps/cli/lib/bin.js
# (override with DSH_BIN), or HARNESS=claude for Claude Agent SDK;
# ECA_PP_CLAUDE_CLI may point at an npm-installed Claude CLI on old glibc hosts.
#
# ALWAYS run on a compute node, never the login node:
#   bash run.sh standardize SRC.h5ad -o OUTDIR [--min-cells N ...]
#   bash run.sh identify-columns STD.h5ad -o OUTDIR [--max-probes N ...]
#   bash run.sh integration-probe STD.h5ad --batch-col COL -o OUTDIR
#   bash run.sh test [pytest args]
#   bash run.sh python script.py
set -euo pipefail

REPO="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PY="${ECA_PP_PYTHON:-/scratch/users/chensj16/venvs/eca-pp-ct/python}"

unset PYTHONPATH || true
# The container wrapper forwards PYTHONPATH via APPTAINERENV_PYTHONPATH, so worktree
# overrides still work; $REPO/src keeps a dirty checkout ahead of the editable install.
export PYTHONPATH="$REPO/src"
# Agent SDK initialize handshake: the npm `claude` CLI cold-starts slowly on a
# compute node (node + NFS-backed ~/.claude + plugins), so allow 3 min instead
# of the SDK's 60 s default (value in ms; eca_pp.agent also retries transients).
export CLAUDE_CODE_STREAM_CLOSE_TIMEOUT="${CLAUDE_CODE_STREAM_CLOSE_TIMEOUT:-180000}"

cmd="${1:-}"
shift || true
case "$cmd" in
  standardize)       exec "$PY" -m eca_pp.standardize "$@" ;;
  identify-columns)  exec "$PY" -m eca_pp.identify_columns "$@" ;;
  integration-probe) exec "$PY" -m eca_pp.probe "$@" ;;
  test)        cd "$REPO"
               exec "$PY" -m pytest -p no:cacheprovider -o addopts="" "$@" ;;
  python)      exec "$PY" "$@" ;;
  *) echo "usage: bash run.sh {standardize|identify-columns|integration-probe|test|python} [args...]" >&2; exit 64 ;;
esac
