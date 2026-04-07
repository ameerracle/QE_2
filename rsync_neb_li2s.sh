#!/bin/bash
set -euo pipefail

REMOTE_USER="anizami"
REMOTE_HOST="rorqual.alliancecan.ca"
REMOTE_DIR="${REMOTE_DIR:-/lustre10/scratch/anizami/QE_2/no_u/NEB_Li2S}"
LOCAL_DIR="${LOCAL_DIR:-/home/ameer_ubuntu/Git_projects/QE_2/no_u/NEB_Li2S}"

# Dry run mode: set DRY_RUN=1 to preview without transferring
DRY_RUN="${DRY_RUN:-0}"

# Keep SSH non-interactive and limit auth/reconnect attempts.
RSYNC_SSH="ssh -o BatchMode=yes -o PreferredAuthentications=publickey -o KbdInteractiveAuthentication=no -o PasswordAuthentication=no -o NumberOfPasswordPrompts=1 -o ConnectionAttempts=1 -o StrictHostKeyChecking=accept-new"

REMOTE_SRC="${REMOTE_DIR%/}/"
LOCAL_DST="${LOCAL_DIR%/}/"

# Base rsync options
RSYNC_OPTS=(-avzP)
[[ "${DRY_RUN}" == "1" ]] && RSYNC_OPTS+=(-n)

# Create local directory
mkdir -p "${LOCAL_DST}"

echo "Pulling .out files from metal-specific subdirectories on ${REMOTE_HOST}:${REMOTE_SRC} → ${LOCAL_DST}"
[[ "${DRY_RUN}" == "1" ]] && echo "[DRY RUN]"

rsync "${RSYNC_OPTS[@]}" \
    --include='*/' \
    --include='*.out' \
    --exclude='*' \
    "${REMOTE_USER}@${REMOTE_HOST}:${REMOTE_SRC}" \
    "${LOCAL_DST}"

echo "Done. Results in ${LOCAL_DST}"