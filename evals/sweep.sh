#!/bin/bash
# Sequential model sweep. Serial on purpose: Ark rate limits have already cost us
# two stages today, and concurrent runs make a 429 look like a model weakness.
cd /scratch/users/chensj16/worktrees/eca-pp-agent-eval
PY=/scratch/users/chensj16/venvs/eca-pp-ct/python
export PYTHONPATH=$PWD/src

while pgrep -f "run_eval.py --tag doubao-turbo" >/dev/null; do sleep 30; done

run() {  # tag harness model
  echo "===== $1 ($2 / $3) ====="
  HARNESS=$2 ECA_PP_AGENT_MODEL=$3 $PY evals/run_eval.py --tag "$1" 2>&1 \
    | grep -vE "^== \[identify columns\] agent:" | tail -22
}
run doubao-pro    openai doubao-seed-2-1-pro-260628
run claude-sonnet-5 claude claude-sonnet-5
run claude-opus-5   claude claude-opus-5
echo "===== 汇总 ====="
$PY evals/run_eval.py --compare
