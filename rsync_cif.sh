#!/bin/bash
set -euo pipefail

REMOTE_USER="anizami"
REMOTE_HOST="rorqual.alliancecan.ca"
REMOTE_DIR="/lustre10/scratch/anizami/QE_2/cif/atomic_cif"
LOCAL_DIR="/home/ameer_ubuntu/Git_projects/QE_2/cif/atomic_cif"

# fallback if remote root is not there
if ssh "${REMOTE_USER}@${REMOTE_HOST}" "[ ! -d '${REMOTE_DIR}' ]"; then
  REMOTE_DIR="/lustre10/scratch/anizami/QE_2/no_u/cif/atomic_cif"
fi

# Dry run mode: set DRY_RUN=1 to preview without transferring
DRY_RUN="${DRY_RUN:-0}"

REMOTE_SRC="${REMOTE_DIR%/}/"
LOCAL_DST="${LOCAL_DIR%/}/"

# Base rsync options
RSYNC_OPTS=(-avzP)
[[ "${DRY_RUN}" == "1" ]] && RSYNC_OPTS+=(-n)

# Create local directory
mkdir -p "${LOCAL_DST}"

# Transfer only .out files from all metal subdirectories recursively
echo "Pulling .out files from all directories on ${REMOTE_HOST}:${REMOTE_SRC} → ${LOCAL_DST}"
[[ "${DRY_RUN}" == "1" ]] && echo "[DRY RUN]"

rsync "${RSYNC_OPTS[@]}" \
    --include='*/' \
    --include='*.out' \
    --exclude='*' \
    "${REMOTE_USER}@${REMOTE_HOST}:${REMOTE_SRC}" \
    "${LOCAL_DST}"

echo "Done. Results in ${LOCAL_DST}"