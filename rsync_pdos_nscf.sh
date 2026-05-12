#!/bin/bash
set -euo pipefail

REMOTE_USER="anizami"
REMOTE_HOST="rorqual.alliancecan.ca"
REMOTE_DIR="${REMOTE_DIR:-/scratch/anizami/QE_2/no_u/PDOS}"
LOCAL_DIR="${LOCAL_DIR:-/home/ameer_ubuntu/Git_projects/QE_2/no_u/PDOS}"

# Dry run mode: set DRY_RUN=1 to preview without transferring.
DRY_RUN="${DRY_RUN:-0}"

REMOTE_SRC="${REMOTE_DIR%/}/"
LOCAL_DST="${LOCAL_DIR%/}/"

# Keep SSH interactive enough for Alliance MFA while still preferring key auth.
RSYNC_SSH="ssh -o BatchMode=no -o PreferredAuthentications=publickey,keyboard-interactive -o KbdInteractiveAuthentication=yes -o PasswordAuthentication=no -o ServerAliveInterval=30 -o ServerAliveCountMax=4 -o StrictHostKeyChecking=accept-new"

# Base rsync options.
RSYNC_OPTS=(-avz --info=progress2)
[[ "${DRY_RUN}" == "1" ]] && RSYNC_OPTS+=(-n)

# Create the local destination folder.
mkdir -p "${LOCAL_DST}"

echo "Pulling *_nscf.out files from *_slab directories on ${REMOTE_HOST}:${REMOTE_SRC} -> ${LOCAL_DST}"
[[ "${DRY_RUN}" == "1" ]] && echo "[DRY RUN]"

REMOTE_FILE_LIST=$(
    ssh -o BatchMode=no \
        -o PreferredAuthentications=publickey,keyboard-interactive \
        -o KbdInteractiveAuthentication=yes \
        -o PasswordAuthentication=no \
        -o ServerAliveInterval=30 \
        -o ServerAliveCountMax=4 \
        -o StrictHostKeyChecking=accept-new \
        "${REMOTE_USER}@${REMOTE_HOST}" \
        "find '${REMOTE_SRC%/}' -mindepth 2 -maxdepth 2 -type f -name '*_nscf.out' | sed 's#^${REMOTE_SRC%/}/##'" \
)

printf '%s\n' "${REMOTE_FILE_LIST}" | rsync -e "${RSYNC_SSH}" "${RSYNC_OPTS[@]}" \
    --files-from=- \
    --relative \
    "${REMOTE_USER}@${REMOTE_HOST}:${REMOTE_SRC}" \
    "${LOCAL_DST}"

echo "Done. Results in ${LOCAL_DST}"
