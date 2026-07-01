#!/usr/bin/env bash
set -euo pipefail

CONFIG="${CONFIG:-$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)/config/config.sh}"
source "${CONFIG}"

ARCHIVE_DIR="${REPO_DIR}/data/archived_inputs"
BULK_ARCHIVE="${ARCHIVE_DIR}/bulk"
PRED_ARCHIVE="${ARCHIVE_DIR}/predictions"

echo "RESULT_ROOT=${RESULT_ROOT}"
echo "RUN_TAG=${RUN_TAG}"
echo "MODEL=${MODEL}"
echo "PROJECTS=${PROJECTS}"

mkdir -p "${RESULT_ROOT}"

for PROJECT in ${PROJECTS}; do
  echo "============================================================"
  echo "PROJECT=${PROJECT}"

  RUN_DIR="${RESULT_ROOT}/${PROJECT}/${RUN_TAG}"
  BULK_DIR="${RUN_DIR}/bulk"
  PRED_DIR="${RUN_DIR}/pred"

  mkdir -p "${BULK_DIR}" "${PRED_DIR}"

  SRC_BULK="${BULK_ARCHIVE}/${PROJECT}/bulk.tsv"

  if [[ ! -f "${SRC_BULK}" ]]; then
    echo "ERROR: archived bulk not found: ${SRC_BULK}" >&2
    exit 1
  fi

  cp -f "${SRC_BULK}" "${BULK_DIR}/bulk.tsv"

  FILES=(
    "TCGA_${PROJECT}_Scaden_${MODEL}_predictions.txt"
    "TCGA_${PROJECT}_Scaden_${MODEL}_sample_by_state.csv"
    "TCGA_${PROJECT}_Scaden_${MODEL}_state_by_sample.csv"
    "TCGA_${PROJECT}_Scaden_${MODEL}_state_QC.csv"
  )

  for f in "${FILES[@]}"; do
    SRC="${PRED_ARCHIVE}/${PROJECT}/${f}"
    DST="${PRED_DIR}/${f}"

    if [[ ! -f "${SRC}" ]]; then
      echo "ERROR: archived prediction not found: ${SRC}" >&2
      exit 1
    fi

    cp -f "${SRC}" "${DST}"
  done

  echo "Copied archived inputs to ${RUN_DIR}"
done

echo "Finished archived input setup."
