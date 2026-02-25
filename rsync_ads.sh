#!/bin/bash
set -euo pipefail

REMOTE_USER="anizami"
REMOTE_HOST="narval.alliancecan.ca"
# allow overriding with environment variables
REMOTE_DIR="${REMOTE_DIR:-/lustre07/scratch/anizami/QE_2/adsorbates}"
LOCAL_DIR="${LOCAL_DIR:-/home/ameer_ubuntu/Git_projects/QE_2/adsorbates}"

# Dry run mode: set DRY_RUN=1 to preview without transferring
DRY_RUN="${DRY_RUN:-0}"

REMOTE_SRC="${REMOTE_DIR%/}/"
LOCAL_DST="${LOCAL_DIR%/}/"

# Base rsync options
RSYNC_OPTS=(-avzP)
[[ "${DRY_RUN}" == "1" ]] && RSYNC_OPTS+=(-n)

# verify remote directory exists before attempting transfer
if ! ssh -q "${REMOTE_USER}@${REMOTE_HOST}" test -d "${REMOTE_DIR}"; then
    echo "ERROR: remote directory '${REMOTE_DIR}' not found on ${REMOTE_HOST}" >&2
    exit 1
fi

# Create local directory
mkdir -p "${LOCAL_DST}"

# default adsorbate list; can override with ADSORBATES_ENV (comma-separated)
ADSORBATES=("Li2S" "Li2S2" "Li2S4" "Li2S6" "Li2S8" "S8")
if [[ -n "${ADSORBATES_ENV:-}" ]]; then
    IFS=',' read -r -a ADSORBATES <<<"$ADSORBATES_ENV"
fi

# build include patterns for each adsorbate
INCLUDE_PATTERNS=()
for a in "${ADSORBATES[@]}"; do
    INCLUDE_PATTERNS+=("${a}/" "${a}/*.pwo")
done
EXCLUDE_PATTERNS=("*/" "*")

echo "Pulling adsorbate folders (${ADSORBATES[*]}) from ${REMOTE_HOST}:${REMOTE_SRC} → ${LOCAL_DST}"
[[ "${DRY_RUN}" == "1" ]] && echo "[DRY RUN]"

# assemble rsync arguments to avoid issues with line breaks
RSYNC_ARGS=("${RSYNC_OPTS[@]}")
for p in "${INCLUDE_PATTERNS[@]}"; do
    RSYNC_ARGS+=("--include=$p")
done
for p in "${EXCLUDE_PATTERNS[@]}"; do
    RSYNC_ARGS+=("--exclude=$p")
done

rsync "${RSYNC_ARGS[@]}" \
    "${REMOTE_USER}@${REMOTE_HOST}:${REMOTE_SRC}" \
    "${LOCAL_DST}"

echo "Done. Results in ${LOCAL_DST}"