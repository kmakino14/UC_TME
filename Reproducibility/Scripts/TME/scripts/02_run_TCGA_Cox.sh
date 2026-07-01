#!/usr/bin/env bash
set -euo pipefail

CONFIG="${CONFIG:-$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)/config/config.sh}"
source "${CONFIG}"

R_SCRIPT="${CORE_DIR}/TCGA_multivariable_Cox.R"
PROJECTS_CSV="$(projects_csv)"

echo "RESULT_ROOT=${RESULT_ROOT}"
echo "RUN_TAG=${RUN_TAG}"
echo "PROJ_SUBDIR=${PROJ_SUBDIR}"
echo "COX_OUTDIR=${COX_OUTDIR}"
echo "PROJECTS_CSV=${PROJECTS_CSV}"

if [[ ! -f "${R_SCRIPT}" ]]; then
  echo "ERROR: Cox R script not found: ${R_SCRIPT}" >&2
  exit 1
fi

for PROJECT in ${PROJECTS}; do
  ASSIGN="${RESULT_ROOT}/${PROJECT}/${RUN_TAG}/${PROJ_SUBDIR}/integrated_TME_assignment.tsv"
  if [[ ! -f "${ASSIGN}" ]]; then
    echo "ERROR: missing assignment file: ${ASSIGN}" >&2
    exit 1
  fi
done

# Clean both the intended short output and the old long output.
LONG_COX_OUTDIR="${RESULT_ROOT}/TCGA_multivariable_Cox_TME_noTME2_gender_cancertype_${RUN_TAG}"

rm -rf "${COX_OUTDIR}"
rm -rf "${LONG_COX_OUTDIR}"

WD_DIR="${WD}" \
BASE_DIR="${RESULT_ROOT}" \
RESULT_ROOT="${RESULT_ROOT}" \
RUN_TAG="${RUN_TAG}" \
MODEL="${MODEL}" \
PROJECTS="${PROJECTS_CSV}" \
TME_KEEP="${TME_KEEP}" \
PROJECTION_SUBDIR="${PROJ_SUBDIR}" \
PROJ_SUBDIR="${PROJ_SUBDIR}" \
CLINICAL_FILE="${CLINICAL_FILE}" \
CDR_FILE="${CDR_FILE}" \
OUT_DIR="${COX_OUTDIR}" \
Rscript "${R_SCRIPT}"

# The original R script may ignore OUT_DIR and create the historical long folder.
# For the clean GitHub workflow, normalize it to TCGA_multivariable_Cox.
if [[ -d "${LONG_COX_OUTDIR}" ]]; then
  echo "Normalizing Cox output directory:"
  echo "  ${LONG_COX_OUTDIR}"
  echo "  -> ${COX_OUTDIR}"

  rm -rf "${COX_OUTDIR}"
  mv "${LONG_COX_OUTDIR}" "${COX_OUTDIR}"
fi

if [[ ! -d "${COX_OUTDIR}" ]]; then
  echo "ERROR: Cox output dir not found after normalization: ${COX_OUTDIR}" >&2
  exit 1
fi

echo "Final Cox output directory:"
echo "  ${COX_OUTDIR}"

echo "Cox output files:"
find "${COX_OUTDIR}" -maxdepth 1 -type f | sort

echo "Finished TCGA Cox."
