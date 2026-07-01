#!/usr/bin/env bash
set -euo pipefail

export REPO_DIR="${REPO_DIR:-$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)}"

# Local working root.
export WD="${WD:-/home/k-makino/wd_home/UC_DOGMA_reseq}"

# Output root.
export RESULT_ROOT="${RESULT_ROOT:-${WD}/Reproducibility/Results/TME}"

# Archived prediction run.
export RUN_TAG="${RUN_TAG:-archived_m256}"
export MODEL="${MODEL:-m256}"

# Space-separated project list for bash loops.
export PROJECTS="${PROJECTS:-BLCA BRCA COAD HNSC KIRC LUAD LUSC OV SKCM STAD UCEC GBM SARC ESCA LIHC KIRP CESC THCA LGG PRAD}"

# Final assignment categories.
export TME_KEEP="${TME_KEEP:-TME1,TME3,TME4,TME5,TME6}"

# Short output names.
export PROJ_SUBDIR="${PROJ_SUBDIR:-module_signature_template_projection}"
export COX_OUTDIR="${COX_OUTDIR:-${RESULT_ROOT}/TCGA_multivariable_Cox}"
export HEATMAP_OUTDIR="${HEATMAP_OUTDIR:-${RESULT_ROOT}/TCGA_assignment_evidence_heatmap}"

# Clinical file.
export CLINICAL_FILE="${CLINICAL_FILE:-/home/k-makino/TCGA/Pancancer/GDCdata/clinical_TCGA_CDR_all_cancers.tsv}"
export CDR_FILE="${CDR_FILE:-${CLINICAL_FILE}}"

# DOGMA/TME resources.
export DOGMA_ABUNDANCE_FILE="${DOGMA_ABUNDANCE_FILE:-${WD}/ECOTYPE/cell_abundance_w_epithelial.txt}"
export DOGMA_TME_ASSIGNMENT_FILE="${DOGMA_TME_ASSIGNMENT_FILE:-${WD}/ECOTYPE/TME_assignment_DOGMA.txt}"
export DOGMA_COUNTERPART_FILE="${DOGMA_COUNTERPART_FILE:-${WD}/ECOTYPE/TME_counterpart.txt}"

export CORE_DIR="${CORE_DIR:-${REPO_DIR}/scripts/core}"
export SCRIPT_DIR="${SCRIPT_DIR:-${REPO_DIR}/scripts}"

projects_csv() {
  echo "${PROJECTS}" | tr ' ' ',' | sed 's/,,*/,/g; s/^,//; s/,$//'
}
