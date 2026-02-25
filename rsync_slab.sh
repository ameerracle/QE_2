#!/bin/bash
set -euo pipefail

REMOTE_USER="anizami"
REMOTE_HOST="rorqual.alliancecan.ca"
REMOTE_DIR="/lustre10/scratch/anizami/QE_2/no_u/slab/input_slab_144"
LOCAL_DIR="/home/ameer_ubuntu/Git_projects/QE_2/no_u/slab/output_slab_144"

# Dry run mode: set DRY_RUN=1 to preview without transferring
DRY_RUN="${DRY_RUN:-0}"

REMOTE_SRC="${REMOTE_DIR%/}/"
LOCAL_DST="${LOCAL_DIR%/}/"

# Base rsync options
RSYNC_OPTS=(-avzP)
[[ "${DRY_RUN}" == "1" ]] && RSYNC_OPTS+=(-n)

# Create local directory
mkdir -p "${LOCAL_DST}"

# Transfer only .out files from *_slab_relax subdirectories
echo "Pulling .out files from *_slab_relax directories on ${REMOTE_HOST}:${REMOTE_SRC} → ${LOCAL_DST}"
[[ "${DRY_RUN}" == "1" ]] && echo "[DRY RUN]"

rsync "${RSYNC_OPTS[@]}" \
    --include='*_slab_relax/' \
    --include='*_slab_relax/*.out' \
    --exclude='*/' \
    --exclude='*' \
    "${REMOTE_USER}@${REMOTE_HOST}:${REMOTE_SRC}" \
    "${LOCAL_DST}"

echo "Done. Results in ${LOCAL_DST}"