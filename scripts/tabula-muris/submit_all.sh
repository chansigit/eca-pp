#!/bin/bash
# Submit every generated job of the given dataset(s) whose run has not already finished with exit 0.
#   bash scripts/tabula-muris/submit_all.sh tabula-muris-senis-drop            # skip finished ones
#   bash scripts/tabula-muris/submit_all.sh --force tabula-sapiens             # resubmit everything
#   bash scripts/tabula-muris/submit_all.sh tabula-muris-drop tabula-muris-facs
set -euo pipefail
OAK_BASE=/home/users/chensj16/oak/data/sc
force=0; sets=()
for a in "$@"; do case "$a" in --force) force=1 ;; *) sets+=("$a") ;; esac; done
[ ${#sets[@]} -gt 0 ] || { echo "usage: submit_all.sh [--force] <dataset-dir> [...]" >&2; exit 64; }

for ds in "${sets[@]}"; do
  short="${ds#tabula-}"
  for job in "$OAK_BASE/$ds"/eca-pp/jobs/*.sbatch; do
    organ="$(basename "$job" .sbatch)"
    st="$SCRATCH/eca-pp-runs/$ds/$organ/status.txt"
    if [ $force -eq 0 ] && [ -f "$st" ] && grep -q '^end=.* exit=0$' "$st"; then
      echo "skip  $ds/$organ (done)"; continue
    fi
    if squeue --me -h -o %j | grep -qx "ecapp-$short-$organ"; then
      echo "skip  $ds/$organ (queued/running)"; continue
    fi
    echo -n "submit $ds/$organ -> "; sbatch "$job"
  done
done
