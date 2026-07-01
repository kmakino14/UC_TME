#!/usr/bin/env bash
set -uo pipefail

# Run module/signature/template projection and integrated-TME KM analysis
# for TCGA pan-cancer projects, excluding TME2/TLS-like from projection targets.
# Projection targets by default: TME1,TME3,TME4,TME5,TME6

PROJECTS=(
  BLCA
  BRCA
  CESC
  COAD
  ESCA
  GBM
  HNSC
  KIRC
  KIRP
  LGG
  LIHC
  LUAD
  LUSC
  OV
  PRAD
  SARC
  SKCM
  STAD
  THCA
  UCEC
)

RUN_TAG="${RUN_TAG:-v15k_ranktemplate_tme2_tme3_tme5_recovery_ns40000_steps40000}"
MODEL="${MODEL:-m256}"
TME_KEEP="${TME_KEEP:-TME1,TME3,TME4,TME5,TME6}"
PROJ_SUBDIR="${PROJ_SUBDIR:-module_signature_template_projection_noTME2}"

WD_DIR="${WD_DIR:-/home/k-makino/wd_home/UC_DOGMA_reseq}"
CODE_DIR="${CODE_DIR:-/home/k-makino/code/UC_DOGMA_reseq}"
BASE_DIR="${BASE_DIR:-${WD_DIR}/scaden/TCGA_Pancancer_TME}"

DOGMA_ABUNDANCE_FILE="${DOGMA_ABUNDANCE_FILE:-${WD_DIR}/ECOTYPE/cell_abundance_w_epithelial.txt}"
DOGMA_TME_ASSIGNMENT_FILE="${DOGMA_TME_ASSIGNMENT_FILE:-${WD_DIR}/ECOTYPE/TME_assignment_DOGMA.txt}"

SCRIPT_PROJECTION="${SCRIPT_PROJECTION:-${CODE_DIR}/260621_scaden_11_tme_module_signature_template_projection_noTME2.R}"
SCRIPT_KM="${SCRIPT_KM:-${CODE_DIR}/260619_scaden_12_integrated_TME_assignment_KM_BLCA.R}"

LOG_ROOT="${LOG_ROOT:-${BASE_DIR}/module_projection_km_batch_logs_v15k_panproject_noTME2}"
mkdir -p "${LOG_ROOT}"

STATUS_TSV="${LOG_ROOT}/batch_status_v15k_panproject_noTME2_$(date +%Y%m%d_%H%M%S).tsv"
printf "PROJECT\tRUN_TAG\tTME_KEEP\tprojection_status\tkm_status\tmessage\n" > "${STATUS_TSV}"

message() {
  echo "[$(date '+%Y-%m-%d %H:%M:%S')] $*"
}

require_file() { [[ -f "$1" ]] || { echo "ERROR: file not found: $1" >&2; exit 1; }; }

message "RUN_TAG=${RUN_TAG}"
message "MODEL=${MODEL}"
message "TME_KEEP=${TME_KEEP}"
message "PROJ_SUBDIR=${PROJ_SUBDIR}"
message "BASE_DIR=${BASE_DIR}"
message "SCRIPT_PROJECTION=${SCRIPT_PROJECTION}"
message "SCRIPT_KM=${SCRIPT_KM}"
message "LOG_ROOT=${LOG_ROOT}"
message "STATUS_TSV=${STATUS_TSV}"

require_file "${SCRIPT_PROJECTION}"
require_file "${SCRIPT_KM}"
require_file "${DOGMA_ABUNDANCE_FILE}"
require_file "${DOGMA_TME_ASSIGNMENT_FILE}"

for PROJECT in "${PROJECTS[@]}"; do
  message "============================================================"
  message "PROJECT=${PROJECT}"

  RUN_DIR="${BASE_DIR}/${PROJECT}/${RUN_TAG}"
  PRED_FILE="${RUN_DIR}/pred/TCGA_${PROJECT}_Scaden_${MODEL}_predictions.txt"
  BULK_FILE="${RUN_DIR}/bulk/bulk.tsv"
  PROJ_DIR="${RUN_DIR}/${PROJ_SUBDIR}"
  ASSIGN_FILE="${PROJ_DIR}/integrated_TME_assignment.tsv"
  KM_DIR="${PROJ_DIR}/KM_integrated_TME"
  KM_PDF="${KM_DIR}/TCGA_${PROJECT}_integrated_TME_OS_KM.pdf"

  LOG_PROJECTION="${LOG_ROOT}/${PROJECT}_${RUN_TAG}.noTME2.projection.log"
  LOG_KM="${LOG_ROOT}/${PROJECT}_${RUN_TAG}.noTME2.km.log"

  if [[ ! -d "${RUN_DIR}" ]]; then
    message "SKIP: run directory not found: ${RUN_DIR}"
    printf "%s\t%s\t%s\tSKIP\tSKIP\trun_dir_not_found\n" "${PROJECT}" "${RUN_TAG}" "${TME_KEEP}" >> "${STATUS_TSV}"
    continue
  fi

  if [[ ! -f "${PRED_FILE}" ]]; then
    message "SKIP: prediction file not found: ${PRED_FILE}"
    printf "%s\t%s\t%s\tSKIP\tSKIP\tprediction_file_not_found\n" "${PROJECT}" "${RUN_TAG}" "${TME_KEEP}" >> "${STATUS_TSV}"
    continue
  fi

  if [[ ! -f "${BULK_FILE}" ]]; then
    message "SKIP: bulk file not found: ${BULK_FILE}"
    printf "%s\t%s\t%s\tSKIP\tSKIP\tbulk_file_not_found\n" "${PROJECT}" "${RUN_TAG}" "${TME_KEEP}" >> "${STATUS_TSV}"
    continue
  fi

  message "Running no-TME2 module/signature/template projection..."
  (
    cd "${CODE_DIR}" || exit 1
    PROJECT="${PROJECT}" \
    RUN_TAG="${RUN_TAG}" \
    MODEL="${MODEL}" \
    BASE_DIR="${BASE_DIR}" \
    PRED_FILE="${PRED_FILE}" \
    BULK_FILE="${BULK_FILE}" \
    OUT_DIR="${PROJ_DIR}" \
    TME_KEEP="${TME_KEEP}" \
    DOGMA_ABUNDANCE_FILE="${DOGMA_ABUNDANCE_FILE}" \
    DOGMA_TME_ASSIGNMENT_FILE="${DOGMA_TME_ASSIGNMENT_FILE}" \
    Rscript "${SCRIPT_PROJECTION}"
  ) > "${LOG_PROJECTION}" 2>&1

  if [[ $? -ne 0 ]]; then
    message "FAILED projection: ${PROJECT}. See ${LOG_PROJECTION}"
    printf "%s\t%s\t%s\tFAIL\tSKIP\tprojection_failed\n" "${PROJECT}" "${RUN_TAG}" "${TME_KEEP}" >> "${STATUS_TSV}"
    continue
  fi

  if [[ ! -f "${ASSIGN_FILE}" ]]; then
    message "FAILED projection: assignment file not created: ${ASSIGN_FILE}"
    printf "%s\t%s\t%s\tFAIL\tSKIP\tassignment_file_not_created\n" "${PROJECT}" "${RUN_TAG}" "${TME_KEEP}" >> "${STATUS_TSV}"
    continue
  fi

  message "Running no-TME2 integrated-TME KM analysis..."
  (
    cd "${CODE_DIR}" || exit 1
    PROJECT="${PROJECT}" \
    RUN_TAG="${RUN_TAG}" \
    MODEL="${MODEL}" \
    ASSIGNMENT_FILE="${ASSIGN_FILE}" \
    OUT_DIR="${KM_DIR}" \
    Rscript "${SCRIPT_KM}"
  ) > "${LOG_KM}" 2>&1

  if [[ $? -ne 0 ]]; then
    message "FAILED KM: ${PROJECT}. See ${LOG_KM}"
    printf "%s\t%s\t%s\tOK\tFAIL\tkm_failed\n" "${PROJECT}" "${RUN_TAG}" "${TME_KEEP}" >> "${STATUS_TSV}"
    continue
  fi

  if [[ ! -f "${KM_PDF}" ]]; then
    message "WARNING: KM finished but expected PDF not found: ${KM_PDF}"
    printf "%s\t%s\t%s\tOK\tWARN\tkm_pdf_not_found\n" "${PROJECT}" "${RUN_TAG}" "${TME_KEEP}" >> "${STATUS_TSV}"
    continue
  fi

  message "DONE: ${PROJECT}. PDF: ${KM_PDF}"
  printf "%s\t%s\t%s\tOK\tOK\tdone\n" "${PROJECT}" "${RUN_TAG}" "${TME_KEEP}" >> "${STATUS_TSV}"
done

message "============================================================"
message "Batch finished. Status: ${STATUS_TSV}"
cat "${STATUS_TSV}"
