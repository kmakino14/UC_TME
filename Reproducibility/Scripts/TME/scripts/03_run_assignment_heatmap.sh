#!/usr/bin/env bash
set -euo pipefail

CONFIG="${CONFIG:-$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)/config/config.sh}"
source "${CONFIG}"

WRAPPER="${CORE_DIR}/run_assignment_heatmap.sh"
PROJECTS_CSV="$(projects_csv)"

echo "RESULT_ROOT=${RESULT_ROOT}"
echo "RUN_TAG=${RUN_TAG}"
echo "PROJ_SUBDIR=${PROJ_SUBDIR}"
echo "HEATMAP_OUTDIR=${HEATMAP_OUTDIR}"
echo "PROJECTS_CSV=${PROJECTS_CSV}"

if [[ ! -f "${WRAPPER}" ]]; then
  echo "ERROR: heatmap wrapper not found: ${WRAPPER}" >&2
  exit 1
fi

for PROJECT in ${PROJECTS}; do
  ASSIGN="${RESULT_ROOT}/${PROJECT}/${RUN_TAG}/${PROJ_SUBDIR}/integrated_TME_assignment.tsv"
  if [[ ! -f "${ASSIGN}" ]]; then
    echo "ERROR: missing assignment file: ${ASSIGN}" >&2
    exit 1
  fi
done

LONG_HEATMAP_OUTDIR="${RESULT_ROOT}/TCGA_noTME2_assignment_evidence_heatmap_${RUN_TAG}"

rm -rf "${HEATMAP_OUTDIR}"
rm -rf "${LONG_HEATMAP_OUTDIR}"

WD_DIR="${WD}" \
BASE_DIR="${RESULT_ROOT}" \
RESULT_ROOT="${RESULT_ROOT}" \
RUN_TAG="${RUN_TAG}" \
MODEL="${MODEL}" \
PROJECTS="${PROJECTS_CSV}" \
TME_KEEP="${TME_KEEP}" \
PROJECTION_SUBDIR="${PROJ_SUBDIR}" \
PROJ_SUBDIR="${PROJ_SUBDIR}" \
OUT_DIR="${HEATMAP_OUTDIR}" \
bash "${WRAPPER}"

# Normalize historical long output directory.
if [[ -d "${LONG_HEATMAP_OUTDIR}" ]]; then
  echo "Normalizing heatmap output directory:"
  echo "  ${LONG_HEATMAP_OUTDIR}"
  echo "  -> ${HEATMAP_OUTDIR}"
  rm -rf "${HEATMAP_OUTDIR}"
  mv "${LONG_HEATMAP_OUTDIR}" "${HEATMAP_OUTDIR}"
fi

if [[ ! -d "${HEATMAP_OUTDIR}" ]]; then
  echo "ERROR: heatmap output dir not found after normalization: ${HEATMAP_OUTDIR}" >&2
  exit 1
fi

# Remove noTME2 from output file names.
find "${HEATMAP_OUTDIR}" -depth -name "*noTME2*" | while IFS= read -r f; do
  dir="$(dirname "${f}")"
  base="$(basename "${f}")"
  newbase="${base//_noTME2/}"
  newbase="${newbase//noTME2_/}"
  newbase="${newbase//noTME2/}"
  newbase="$(echo "${newbase}" | sed 's/__*/_/g; s/_\././g; s/^_//; s/_$//')"
  new="${dir}/${newbase}"

  if [[ "${f}" != "${new}" ]]; then
    echo "Renaming:"
    echo "  ${f}"
    echo "  -> ${new}"
    mv "${f}" "${new}"
  fi
done

echo "Final heatmap output directory:"
echo "  ${HEATMAP_OUTDIR}"

find "${HEATMAP_OUTDIR}" -maxdepth 1 -type f | sort

echo "Finished assignment evidence heatmap."
