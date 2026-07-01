#!/usr/bin/env bash
set -euo pipefail

CONFIG="${CONFIG:-$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)/config/config.sh}"
source "${CONFIG}"

PROJECTS_CSV="$(projects_csv)"
WRAPPER="${CORE_DIR}/run_TME_projection.sh"

echo "RESULT_ROOT=${RESULT_ROOT}"
echo "RUN_TAG=${RUN_TAG}"
echo "MODEL=${MODEL}"
echo "PROJ_SUBDIR=${PROJ_SUBDIR}"
echo "PROJECTS=${PROJECTS}"
echo "PROJECTS_CSV=${PROJECTS_CSV}"

if [[ ! -f "${WRAPPER}" ]]; then
  echo "ERROR: projection wrapper not found: ${WRAPPER}" >&2
  exit 1
fi

for PROJECT in ${PROJECTS}; do
  PRED="${RESULT_ROOT}/${PROJECT}/${RUN_TAG}/pred/TCGA_${PROJECT}_Scaden_${MODEL}_predictions.txt"
  if [[ ! -f "${PRED}" ]]; then
    echo "ERROR: missing prediction file: ${PRED}" >&2
    echo "Run scripts/00_use_archived_inputs.sh first." >&2
    exit 1
  fi
done

# Remove old long log directories before running.
rm -rf "${RESULT_ROOT}/log"
find "${RESULT_ROOT}" -maxdepth 1 -type d -name "module_projection_km_batch_logs*" -exec rm -rf {} +

WD_DIR="${WD}" \
BASE_DIR="${RESULT_ROOT}" \
RESULT_ROOT="${RESULT_ROOT}" \
RUN_TAG="${RUN_TAG}" \
MODEL="${MODEL}" \
PROJECTS="${PROJECTS}" \
PROJECTS_CSV="${PROJECTS_CSV}" \
TME_KEEP="${TME_KEEP}" \
PROJECTION_SUBDIR="${PROJ_SUBDIR}" \
PROJ_SUBDIR="${PROJ_SUBDIR}" \
DOGMA_ABUNDANCE_FILE="${DOGMA_ABUNDANCE_FILE}" \
DOGMA_TME_ASSIGNMENT_FILE="${DOGMA_TME_ASSIGNMENT_FILE}" \
DOGMA_COUNTERPART_FILE="${DOGMA_COUNTERPART_FILE}" \
bash "${WRAPPER}"

# Core scripts may still write module_signature_template_projection_noTME2.
# Normalize each project output to the short directory name.
for PROJECT in ${PROJECTS}; do
  RUN_DIR="${RESULT_ROOT}/${PROJECT}/${RUN_TAG}"
  SHORT_DIR="${RUN_DIR}/${PROJ_SUBDIR}"
  OLD_DIR="${RUN_DIR}/module_signature_template_projection_noTME2"

  if [[ ! -f "${SHORT_DIR}/integrated_TME_assignment.tsv" && -f "${OLD_DIR}/integrated_TME_assignment.tsv" ]]; then
    echo "Copying old projection directory to short name for ${PROJECT}"
    rm -rf "${SHORT_DIR}"
    rsync -a "${OLD_DIR}/" "${SHORT_DIR}/"
  fi

  if [[ ! -f "${SHORT_DIR}/integrated_TME_assignment.tsv" ]]; then
    echo "ERROR: missing assignment file: ${SHORT_DIR}/integrated_TME_assignment.tsv" >&2
    exit 1
  fi

  # Remove the old noTME2 directory from clean outputs.
  if [[ -d "${OLD_DIR}" && "${OLD_DIR}" != "${SHORT_DIR}" ]]; then
    rm -rf "${OLD_DIR}"
  fi

  echo "OK: ${PROJECT} assignment"
done

# Normalize projection batch log directory.
LOG_DIR="${RESULT_ROOT}/log"
mkdir -p "${LOG_DIR}"

while IFS= read -r d; do
  [[ -d "${d}" ]] || continue
  echo "Moving projection log directory:"
  echo "  ${d}"
  echo "  -> ${LOG_DIR}"
  rsync -a "${d}/" "${LOG_DIR}/"
  rm -rf "${d}"
done < <(find "${RESULT_ROOT}" -maxdepth 1 -type d -name "module_projection_km_batch_logs*" | sort)

echo "Final projection log directory:"
echo "  ${LOG_DIR}"

echo "Finished TME projection."
