#!/bin/bash
# Portability proof: inside a stock python:3.12-slim container (no dl2025 venv, no
# Lmod modules, no cluster fixups), build a venv from LOCAL SOURCE ONLY — stancounts
# + ecasteps via `pip install <dir>` — and run the same acceptance suite run.sh runs
# natively. Third-party deps (numpy/scipy/h5py/anndata) come from PyPI wheels.
#
#   bash scripts/test-in-container.sh [image.sif]
set -euo pipefail

REPO="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
SIF="${1:-$REPO/python312-slim.sif}"
STANCOUNTS_SRC="${STANCOUNTS_SRC:-/home/users/chensj16/s/projects/stancounts}"
STANGENE_SRC="${STANGENE_SRC:-/home/users/chensj16/s/projects/stangene}"
PIP_CACHE="${PIP_CACHE:-$SCRATCH/pip-cache}"

# --cleanenv: no host env leaks in (PYTHONPATH etc.) — the point is purity.
exec apptainer exec --cleanenv \
    -B /scratch,/home/users \
    --env REPO="$REPO",STANCOUNTS_SRC="$STANCOUNTS_SRC",STANGENE_SRC="$STANGENE_SRC",PIP_CACHE="$PIP_CACHE" \
    "$SIF" bash -c '
  set -euo pipefail
  echo "[container] $(python --version) @ $(cat /etc/os-release | grep PRETTY | cut -d= -f2)"
  python -m venv --clear "$REPO/.ctr-venv"
  V="$REPO/.ctr-venv/bin"
  echo "[container] pip install from local sources ..."
  "$V/pip" install -q --cache-dir "$PIP_CACHE" \
      "$STANCOUNTS_SRC" "$STANGENE_SRC" "$REPO" pytest
  cd "$REPO"
  echo "[container] running acceptance suite ..."
  exec "$V/python" -m pytest tests -q
'
