#!/bin/bash
set -euo pipefail

REMOTE_USER="anizami"
REMOTE_HOST="rorqual.alliancecan.ca"
REMOTE_DIR="${REMOTE_DIR:-/lustre10/scratch/anizami/QE_2/no_u/NEB_Li}"
LOCAL_DIR="${LOCAL_DIR:-/home/ameer_ubuntu/Git_projects/QE_2/no_u/NEB_Li}"

# Dry run mode: set DRY_RUN=1 to preview without transferring
DRY_RUN="${DRY_RUN:-0}"

REMOTE_SRC="${REMOTE_DIR%/}/"
LOCAL_DST="${LOCAL_DIR%/}/"

# Base rsync options
RSYNC_OPTS=(-avzP)
[[ "${DRY_RUN}" == "1" ]] && RSYNC_OPTS+=(-n)

# Create local directory
mkdir -p "${LOCAL_DST}"

# Pull only the requested Li relaxation folders and their .out files
echo "Pulling .out files from Li subfolders on ${REMOTE_HOST}:${REMOTE_SRC} → ${LOCAL_DST}"
[[ "${DRY_RUN}" == "1" ]] && echo "[DRY RUN]"

rsync "${RSYNC_OPTS[@]}" \
    --include='*_NEB_Li/' \
    --include='*_NEB_Li/*.out' \
    --include='*_NEB_Li/*.path' \
    --include='*_NEB_Li/*.axsf' \
    --include='*_NEB_Li/*.xyz' \
    --exclude='*/' \
    --exclude='*' \
    "${REMOTE_USER}@${REMOTE_HOST}:${REMOTE_SRC}" \
    "${LOCAL_DST}"

echo "Done. Results in ${LOCAL_DST}"