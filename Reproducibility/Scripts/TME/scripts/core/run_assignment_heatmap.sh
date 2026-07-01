#!/usr/bin/env bash
set -euo pipefail

CODE_DIR="${CODE_DIR:-/home/k-makino/code/UC_DOGMA_reseq}"
RUN_TAG="${RUN_TAG:-v15k_ranktemplate_tme2_tme3_tme5_recovery_ns40000_steps40000}"
PROJECTION_SUBDIR="${PROJECTION_SUBDIR:-module_signature_template_projection_noTME2}"

cd "${CODE_DIR}"

RUN_TAG="${RUN_TAG}" \
PROJECTION_SUBDIR="${PROJECTION_SUBDIR}" \
ORDER_WITHIN_TME="${ORDER_WITHIN_TME:-margin}" \
CLUSTER_ROWS="${CLUSTER_ROWS:-1}" \
SHOW_ROW_DEND="${SHOW_ROW_DEND:-0}" \
ORDER_ROW_CLUSTERS_BY_TME="${ORDER_ROW_CLUSTERS_BY_TME:-1}" \
SPLIT_ROWS_BY_TME_CLUSTER="${SPLIT_ROWS_BY_TME_CLUSTER:-1}" \
ROW_CLUSTER_K="${ROW_CLUSTER_K:-0}" \
ROW_CLUSTER_METHOD="${ROW_CLUSTER_METHOD:-ward.D2}" \
ROW_CLUSTER_ORDER_TME="${ROW_CLUSTER_ORDER_TME:-TME1,TME3,TME4,TME5,TME6}" \
FORCE_RASTER="${FORCE_RASTER:-1}" \
Rscript "${CODE_DIR}/260621_tcga_noTME2_assignment_evidence_heatmap_panproject_v8_tmeordered_rowclusters.R"
