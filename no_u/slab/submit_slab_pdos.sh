#!/usr/bin/env bash
set -euo pipefail

ACCOUNT="def-peslherb"
TASKS_PER_NODE=64
MEM="140G"
TIME="15:00:00"

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PY_SCRIPT="$SCRIPT_DIR/run_slab_pdos.py"

SLABS=(TiN VN ScN NbN ZrN)

for slab in "${SLABS[@]}"; do
  sbatch \
    --account="$ACCOUNT" \
    --nodes=1 \
    --ntasks="$TASKS_PER_NODE" \
    --cpus-per-task=1 \
    --mem="$MEM" \
    --time="$TIME" \
    --job-name="pdos_${slab}" \
    --output="$SCRIPT_DIR/${slab}_pdos_%j.out" \
    --error="$SCRIPT_DIR/${slab}_pdos_%j.err" \
    --wrap="cd '$SCRIPT_DIR' && python '$PY_SCRIPT' --structure '$slab' --run"
done
