#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(data.table)
  library(dplyr)
  library(stringr)
  library(tibble)
  library(ggplot2)
})

## ============================================================
## TCGA pan-cancer no-TME2 assignment evidence heatmap v8 TME-ordered row clusters
##
## Purpose:
##   Visualize the matrix used as evidence for assigning TCGA bulk samples
##   to TME1, TME3, TME4, TME5, or TME6.
##
## Input per cancer type:
##   <BASE_DIR>/<PROJECT>/<RUN_TAG>/<PROJECTION_SUBDIR>/
##     integrated_TME_assignment.tsv
##     integrated_TME_scores_ranknorm.tsv
##     scaden_module_scores_ranknorm.tsv
##     bulk_gene_signature_scores_ranknorm.tsv
##     TCGA_to_DOGMA_TME_template_pearson.tsv
##
## Output:
##   1) Pan-cancer patient x feature matrices
##   2) One combined heatmap with all TCGA samples mixed and split by TME
##   3) Per-TME heatmaps and matrices
## ============================================================

get_env <- function(x, default = NULL) {
  v <- Sys.getenv(x, unset = NA_character_)
  if (is.na(v) || v == "") return(default)
  v
}

BASE_DIR <- get_env("BASE_DIR", "/home/k-makino/wd_home/UC_DOGMA_reseq/scaden/TCGA_Pancancer_TME")
RUN_TAG <- get_env("RUN_TAG", "v15k_ranktemplate_tme2_tme3_tme5_recovery_ns40000_steps40000")
PROJECTION_SUBDIR <- get_env("PROJECTION_SUBDIR", "module_signature_template_projection_noTME2")
PROJECTS <- strsplit(get_env("PROJECTS", "BLCA,BRCA,CESC,COAD,ESCA,GBM,HNSC,KIRC,KIRP,LGG,LIHC,LUAD,LUSC,OV,PRAD,SARC,SKCM,STAD,THCA,UCEC"), ",")[[1]] |> trimws()
TME_KEEP <- strsplit(get_env("TME_KEEP", "TME1,TME3,TME4,TME5,TME6"), ",")[[1]] |> trimws()

## Which feature groups to include in the evidence heatmap.
INCLUDE_INTEGRATED <- as.logical(as.integer(get_env("INCLUDE_INTEGRATED", "1")))
INCLUDE_SCADEN_MODULES <- as.logical(as.integer(get_env("INCLUDE_SCADEN_MODULES", "1")))
INCLUDE_BULK_SIGNATURES <- as.logical(as.integer(get_env("INCLUDE_BULK_SIGNATURES", "1")))
INCLUDE_TEMPLATE <- as.logical(as.integer(get_env("INCLUDE_TEMPLATE", "1")))

## Limit columns if you need a quick preview. 0 means all samples.
MAX_SAMPLES_PER_TME <- as.integer(get_env("MAX_SAMPLES_PER_TME", "0"))
ORDER_WITHIN_TME <- get_env("ORDER_WITHIN_TME", "margin")
HEATMAP_ZLIM <- as.numeric(get_env("HEATMAP_ZLIM", "2.5"))
## By default, keep zero-variance features and draw them as z=0.
## In some rank/template outputs, all selected evidence features can be constant
## after pan-cancer merging; dropping them produced empty heatmaps.
DROP_ZERO_VARIANCE_FEATURES <- as.logical(as.integer(get_env("DROP_ZERO_VARIANCE_FEATURES", "0")))
MIN_FINITE_PER_FEATURE <- as.integer(get_env("MIN_FINITE_PER_FEATURE", "5"))
RASTER_COLUMN_THRESHOLD <- as.integer(get_env("RASTER_COLUMN_THRESHOLD", "800"))
FORCE_RASTER <- as.logical(as.integer(get_env("FORCE_RASTER", "0")))

## Row direction in ComplexHeatmap corresponds to evidence variables/features
## because the plotted matrix is transposed as features x samples.
CLUSTER_ROWS <- as.logical(as.integer(get_env("CLUSTER_ROWS", "1")))
SHOW_ROW_DEND <- as.logical(as.integer(get_env("SHOW_ROW_DEND", "1")))
SHOW_FEATURE_GROUP_ROW_ANNOTATION <- as.logical(as.integer(get_env("SHOW_FEATURE_GROUP_ROW_ANNOTATION", "1")))

## v8: cluster features first, then order feature clusters by the TME in which
## each cluster is most enriched.  This puts TME1-like feature clusters at the
## top, then TME3-like, TME4-like, TME5-like, and TME6-like.
ORDER_ROW_CLUSTERS_BY_TME <- as.logical(as.integer(get_env("ORDER_ROW_CLUSTERS_BY_TME", "1")))
SPLIT_ROWS_BY_TME_CLUSTER <- as.logical(as.integer(get_env("SPLIT_ROWS_BY_TME_CLUSTER", "1")))
ROW_CLUSTER_K <- as.integer(get_env("ROW_CLUSTER_K", "0"))
ROW_CLUSTER_METHOD <- get_env("ROW_CLUSTER_METHOD", "ward.D2")
ROW_CLUSTER_ORDER_TME <- strsplit(get_env("ROW_CLUSTER_ORDER_TME", paste(TME_KEEP, collapse = ",")), ",")[[1]] |> trimws()

OUT_DIR <- get_env(
  "OUT_DIR",
  file.path(BASE_DIR, paste0("TCGA_noTME2_assignment_evidence_heatmap_", RUN_TAG))
)
dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)

message("BASE_DIR:          ", BASE_DIR)
message("RUN_TAG:           ", RUN_TAG)
message("PROJECTION_SUBDIR: ", PROJECTION_SUBDIR)
message("PROJECTS:          ", paste(PROJECTS, collapse = ", "))
message("TME_KEEP:          ", paste(TME_KEEP, collapse = ", "))
message("OUT_DIR:           ", OUT_DIR)
message("MAX_SAMPLES_PER_TME: ", MAX_SAMPLES_PER_TME)
message("DROP_ZERO_VARIANCE_FEATURES: ", DROP_ZERO_VARIANCE_FEATURES)
message("MIN_FINITE_PER_FEATURE: ", MIN_FINITE_PER_FEATURE)
message("RASTER_COLUMN_THRESHOLD: ", RASTER_COLUMN_THRESHOLD)
message("FORCE_RASTER: ", FORCE_RASTER)
message("ORDER_WITHIN_TME: ", ORDER_WITHIN_TME)
message("CLUSTER_ROWS: ", CLUSTER_ROWS)
message("SHOW_ROW_DEND: ", SHOW_ROW_DEND)
message("SHOW_FEATURE_GROUP_ROW_ANNOTATION: ", SHOW_FEATURE_GROUP_ROW_ANNOTATION)
message("ORDER_ROW_CLUSTERS_BY_TME: ", ORDER_ROW_CLUSTERS_BY_TME)
message("SPLIT_ROWS_BY_TME_CLUSTER: ", SPLIT_ROWS_BY_TME_CLUSTER)
message("ROW_CLUSTER_K: ", ROW_CLUSTER_K)
message("ROW_CLUSTER_METHOD: ", ROW_CLUSTER_METHOD)
message("ROW_CLUSTER_ORDER_TME: ", paste(ROW_CLUSTER_ORDER_TME, collapse = ", "))
message("==============================")

## ----------------------------
## Helpers
## ----------------------------

normalize_sample_id <- function(x) {
  x <- as.character(x)
  x <- trimws(x)
  x <- gsub("\\.", "-", x)
  is_tcga <- grepl("^TCGA-[A-Za-z0-9]{2}-[A-Za-z0-9]{4}-[0-9]{2}", x)
  x[is_tcga] <- substr(x[is_tcga], 1, 15)
  x
}

clean_feature_name <- function(x) {
  x <- as.character(x)
  x <- gsub("^TME1_desert_epithelial$", "Score_TME1_desert_epithelial", x)
  x <- gsub("^TME2_TLS_active$", "Score_TME2_TLS_active", x)
  x <- gsub("^TME3_naive_memory$", "Score_TME3_naive_memory", x)
  x <- gsub("^TME4_dysfunctional$", "Score_TME4_dysfunctional", x)
  x <- gsub("^TME5_cytotoxic_Th17$", "Score_TME5_cytotoxic_Th17", x)
  x <- gsub("^TME6_stroma$", "Score_TME6_stroma", x)
  x
}

clean_tme <- function(x) {
  x <- as.character(x)
  x <- trimws(x)
  x <- gsub("_", "", x)
  x <- toupper(x)
  ifelse(grepl("TME[0-9]", x), sub(".*(TME[0-9]).*", "\\1", x), x)
}

first_existing_ci <- function(x, candidates) {
  lx <- tolower(x)
  lc <- tolower(candidates)
  idx <- match(lc, lx)
  idx <- idx[!is.na(idx)]
  if (length(idx) == 0) return(NA_character_)
  x[idx[[1]]]
}

read_table_with_sample <- function(path, sample_candidates = c("sample", "Sample", "ID", "sample_id", "patient_id", "barcode")) {
  if (!file.exists(path)) stop("file not found: ", path)
  dt <- fread(path, data.table = FALSE, check.names = FALSE)
  cn <- colnames(dt)
  sample_col <- first_existing_ci(cn, sample_candidates)
  if (is.na(sample_col)) sample_col <- cn[[1]]
  row_ids <- normalize_sample_id(dt[[sample_col]])
  df <- dt[, setdiff(cn, sample_col), drop = FALSE]
  for (cc in colnames(df)) {
    if (!is.numeric(df[[cc]])) {
      suppressWarnings(df[[cc]] <- as.numeric(df[[cc]]))
    }
  }
  rownames(df) <- row_ids
  if (anyDuplicated(rownames(df))) {
    ## Average duplicate aliquots if any.
    m <- as.matrix(df)
    storage.mode(m) <- "numeric"
    m <- rowsum(m, group = rownames(df), reorder = FALSE)
    n <- as.numeric(table(factor(rownames(df), levels = rownames(m))))
    m <- sweep(m, 1, n, "/")
    df <- as.data.frame(m, check.names = FALSE)
  }
  df
}

read_assignment <- function(path, project) {
  if (!file.exists(path)) stop("assignment file not found: ", path)
  dt <- fread(path, data.table = FALSE, check.names = FALSE)
  cn <- colnames(dt)
  sample_col <- first_existing_ci(cn, c("sample", "Sample", "ID", "sample_id", "patient_id", "barcode"))
  if (is.na(sample_col)) sample_col <- cn[[1]]
  tme_col <- first_existing_ci(cn, c("predicted_TME", "TME", "TME_assignment", "assigned_TME", "best_TME", "cluster", "group"))
  if (is.na(tme_col)) {
    cand <- cn[grepl("TME", cn, ignore.case = TRUE)]
    if (length(cand) > 0) tme_col <- cand[[1]]
  }
  if (is.na(tme_col)) stop("Could not detect TME column in: ", path)
  margin_col <- first_existing_ci(cn, c("margin", "score_margin", "delta", "best_margin"))
  max_col <- first_existing_ci(cn, c("max_score", "best_score"))
  second_col <- first_existing_ci(cn, c("second_score", "runnerup_score"))
  out <- data.frame(
    sample = normalize_sample_id(dt[[sample_col]]),
    cancer_type = project,
    predicted_TME = clean_tme(dt[[tme_col]]),
    max_score = if (!is.na(max_col)) suppressWarnings(as.numeric(dt[[max_col]])) else NA_real_,
    second_score = if (!is.na(second_col)) suppressWarnings(as.numeric(dt[[second_col]])) else NA_real_,
    margin = if (!is.na(margin_col)) suppressWarnings(as.numeric(dt[[margin_col]])) else NA_real_,
    stringsAsFactors = FALSE
  )
  out <- out %>%
    filter(!is.na(sample), sample != "", predicted_TME %in% TME_KEEP) %>%
    distinct(sample, .keep_all = TRUE)
  out
}

zscale_columns <- function(df) {
  m <- as.matrix(df)
  storage.mode(m) <- "numeric"
  mu <- colMeans(m, na.rm = TRUE)
  sdv <- apply(m, 2, sd, na.rm = TRUE)
  sdv[!is.finite(sdv) | sdv <= 1e-8] <- 1
  z <- sweep(m, 2, mu, "-")
  z <- sweep(z, 2, sdv, "/")
  z[!is.finite(z)] <- 0
  z
}

cap_matrix <- function(m, zlim = 2.5) {
  m <- as.matrix(m)
  m[m > zlim] <- zlim
  m[m < -zlim] <- -zlim
  m
}

make_feature_groups <- function(feature_names) {
  g <- rep("Other", length(feature_names))
  g[grepl("^Score_TME", feature_names)] <- "Integrated TME score"
  g[grepl("^Template", feature_names)] <- "DOGMA template similarity"
  g[grepl("^Sig_", feature_names)] <- "Bulk gene signature"
  g[grepl("^(Epithelial|TLS_B|Naive_memory|Cytotoxic_Th17|Dysfunctional|MDSC|Myeloid|Stroma|NK_ILC|Minor_T|CD8_|MDSC_)", feature_names)] <- "Scaden module"
  factor(g, levels = c("Integrated TME score", "Scaden module", "Bulk gene signature", "DOGMA template similarity", "Other"))
}

order_rows_by_tme_characteristic_clusters <- function(mat, ann, out_prefix = NULL) {
  ## mat: features x samples, already z-scaled/capped.
  ## ann: sample annotation ordered to colnames(mat).
  out <- list(
    mat = mat,
    row_dend = FALSE,
    row_cluster = factor(rep("all_features", nrow(mat)), levels = "all_features"),
    row_cluster_tme = factor(rep("Other", nrow(mat)), levels = c(ROW_CLUSTER_ORDER_TME, "Other")),
    cluster_info = data.frame()
  )
  if (!ORDER_ROW_CLUSTERS_BY_TME || !CLUSTER_ROWS || nrow(mat) < 2) {
    return(out)
  }

  row_sd <- apply(mat, 1, sd, na.rm = TRUE)
  variable_rows <- is.finite(row_sd) & row_sd > 1e-8
  if (sum(variable_rows) < 2) {
    return(out)
  }

  ## Cluster all rows.  Constant rows are allowed; dist() is still finite after
  ## non-finite values have been replaced by 0 upstream.
  row_dist <- stats::dist(mat)
  hc <- stats::hclust(row_dist, method = ROW_CLUSTER_METHOD)
  leaf_order <- rownames(mat)[hc$order]

  k <- ROW_CLUSTER_K
  if (is.na(k) || k <= 1) {
    ## For ~30-40 evidence features, this usually gives 6-8 interpretable
    ## feature clusters.  Users can override with ROW_CLUSTER_K.
    k <- min(nrow(mat), max(length(TME_KEEP), ceiling(sqrt(nrow(mat)))))
  }
  k <- max(1, min(k, nrow(mat)))
  cl <- stats::cutree(hc, k = k)
  names(cl) <- rownames(mat)

  tme_levels <- ROW_CLUSTER_ORDER_TME[ROW_CLUSTER_ORDER_TME %in% as.character(unique(ann$predicted_TME))]
  if (length(tme_levels) == 0) {
    tme_levels <- as.character(unique(ann$predicted_TME))
  }
  ann_tme <- factor(as.character(ann$predicted_TME), levels = tme_levels)

  tme_mean <- vapply(tme_levels, function(tm) {
    idx <- which(ann_tme == tm)
    if (length(idx) == 0) return(rep(NA_real_, nrow(mat)))
    rowMeans(mat[, idx, drop = FALSE], na.rm = TRUE)
  }, numeric(nrow(mat)))
  rownames(tme_mean) <- rownames(mat)

  cluster_ids <- sort(unique(cl))
  cl_info <- lapply(cluster_ids, function(cid) {
    rows <- names(cl)[cl == cid]
    cm <- colMeans(tme_mean[rows, , drop = FALSE], na.rm = TRUE)
    if (all(!is.finite(cm))) {
      best_tme <- "Other"
      best_score <- NA_real_
    } else {
      best_tme <- names(which.max(cm))
      best_score <- unname(max(cm, na.rm = TRUE))
    }
    data.frame(
      row_cluster_id = cid,
      row_cluster = paste0(best_tme, "_cluster", sprintf("%02d", cid)),
      characteristic_TME = best_tme,
      mean_enrichment = best_score,
      n_features = length(rows),
      features = paste(rows, collapse = ";"),
      stringsAsFactors = FALSE
    )
  }) |> bind_rows()
  cl_info$tme_rank <- match(cl_info$characteristic_TME, ROW_CLUSTER_ORDER_TME)
  cl_info$tme_rank[is.na(cl_info$tme_rank)] <- length(ROW_CLUSTER_ORDER_TME) + 1
  cl_info <- cl_info |> arrange(tme_rank, desc(mean_enrichment), row_cluster_id)

  ordered_rows <- unlist(lapply(cl_info$row_cluster_id, function(cid) {
    leaf_order[cl[leaf_order] == cid]
  }), use.names = FALSE)
  ordered_rows <- ordered_rows[ordered_rows %in% rownames(mat)]

  mat2 <- mat[ordered_rows, , drop = FALSE]
  cl2 <- cl[ordered_rows]
  cl_info_by_id <- cl_info[match(cl2, cl_info$row_cluster_id), , drop = FALSE]
  row_cluster <- factor(cl_info_by_id$row_cluster, levels = unique(cl_info$row_cluster))
  row_cluster_tme <- factor(cl_info_by_id$characteristic_TME, levels = c(ROW_CLUSTER_ORDER_TME, "Other"))

  if (!is.null(out_prefix)) {
    fwrite(cl_info |> select(-tme_rank), paste0(out_prefix, "_row_cluster_summary.tsv"), sep = "	")
    fwrite(
      data.frame(
        feature = rownames(mat2),
        row_cluster = as.character(row_cluster),
        characteristic_TME = as.character(row_cluster_tme),
        stringsAsFactors = FALSE
      ),
      paste0(out_prefix, "_feature_order.tsv"),
      sep = "	"
    )
  }

  list(
    mat = mat2,
    row_dend = FALSE,
    row_cluster = row_cluster,
    row_cluster_tme = row_cluster_tme,
    cluster_info = cl_info
  )
}

## ----------------------------
## Load all projects
## ----------------------------

assignment_list <- list()
feature_list <- list()
load_status <- list()

for (project in PROJECTS) {
  message("==============================")
  message("PROJECT: ", project)
  proj_dir <- file.path(BASE_DIR, project, RUN_TAG, PROJECTION_SUBDIR)
  assign_file <- file.path(proj_dir, "integrated_TME_assignment.tsv")
  integrated_file <- file.path(proj_dir, "integrated_TME_scores_ranknorm.tsv")
  scaden_file <- file.path(proj_dir, "scaden_module_scores_ranknorm.tsv")
  sig_file <- file.path(proj_dir, "bulk_gene_signature_scores_ranknorm.tsv")
  template_file <- file.path(proj_dir, "TCGA_to_DOGMA_TME_template_pearson.tsv")

  if (!file.exists(assign_file)) {
    warning("Skip ", project, ": missing assignment file: ", assign_file)
    load_status[[project]] <- data.frame(project = project, status = "skip", message = "assignment_missing", stringsAsFactors = FALSE)
    next
  }

  a <- read_assignment(assign_file, project)
  if (nrow(a) == 0) {
    warning("Skip ", project, ": no assignment rows after TME_KEEP filter")
    load_status[[project]] <- data.frame(project = project, status = "skip", message = "no_assignment_after_filter", stringsAsFactors = FALSE)
    next
  }

  mats <- list()

  if (INCLUDE_INTEGRATED) {
    if (!file.exists(integrated_file)) stop("missing integrated score file for ", project, ": ", integrated_file)
    x <- read_table_with_sample(integrated_file)
    colnames(x) <- clean_feature_name(colnames(x))
    score_to_tme <- ifelse(grepl("^Score_TME[0-9]", colnames(x)), sub("^Score_(TME[0-9]).*", "\\1", colnames(x)), NA_character_)
    keep_cols <- is.na(score_to_tme) | score_to_tme %in% TME_KEEP
    x <- x[, keep_cols, drop = FALSE]
    mats[["integrated"]] <- x
  }

  if (INCLUDE_SCADEN_MODULES) {
    if (!file.exists(scaden_file)) stop("missing scaden module file for ", project, ": ", scaden_file)
    x <- read_table_with_sample(scaden_file)
    mats[["scaden"]] <- x
  }

  if (INCLUDE_BULK_SIGNATURES) {
    if (!file.exists(sig_file)) stop("missing bulk signature file for ", project, ": ", sig_file)
    x <- read_table_with_sample(sig_file)
    mats[["signature"]] <- x
  }

  if (INCLUDE_TEMPLATE) {
    if (!file.exists(template_file)) stop("missing template similarity file for ", project, ": ", template_file)
    x <- read_table_with_sample(template_file)
    ## Keep only no-TME2 template similarity columns.
    keep_cols <- colnames(x)[sub("^TemplatePearson_(TME[0-9]).*", "\\1", colnames(x)) %in% TME_KEEP]
    x <- x[, keep_cols, drop = FALSE]
    mats[["template"]] <- x
  }

  common <- Reduce(intersect, c(list(a$sample), lapply(mats, rownames)))
  if (length(common) == 0) {
    warning("Skip ", project, ": no common samples between assignment and evidence matrices")
    load_status[[project]] <- data.frame(project = project, status = "skip", message = "no_common_samples", stringsAsFactors = FALSE)
    next
  }

  a <- a[match(common, a$sample), , drop = FALSE]
  ## Use a project-prefixed key for pan-cancer merging.  Direct rbind() of
  ## a named list of data.frames can silently prefix row names with project
  ## names and break downstream matching by TCGA barcode.
  a$sample_key <- paste(project, a$sample, sep = "__")

  feat <- do.call(cbind, lapply(mats, function(x) x[common, , drop = FALSE]))
  ## Make column names unique in case of accidental overlap.
  colnames(feat) <- make.unique(colnames(feat), sep = "__")
  rownames(feat) <- a$sample_key

  assignment_list[[project]] <- a
  feature_list[[project]] <- feat
  load_status[[project]] <- data.frame(
    project = project,
    status = "ok",
    message = "loaded",
    n_assignment = nrow(a),
    n_features = ncol(feat),
    stringsAsFactors = FALSE
  )
  message("Loaded samples: ", nrow(a), "; features: ", ncol(feat))
}

status_tbl <- bind_rows(load_status)
fwrite(status_tbl, file.path(OUT_DIR, "TCGA_noTME2_heatmap_load_status.tsv"), sep = "\t")

assignment <- bind_rows(assignment_list)
## Combine feature matrices without allowing R to prefix row names with project
## names. feature_list rows are already keyed by assignment$sample_key.
feature_raw <- do.call(rbind, unname(feature_list))

if (nrow(assignment) == 0 || nrow(feature_raw) == 0) {
  stop("No samples loaded. Check PROJECTION_SUBDIR/RUN_TAG and projection outputs.")
}

if (!"sample_key" %in% colnames(assignment)) {
  assignment$sample_key <- paste(assignment$cancer_type, assignment$sample, sep = "__")
}

missing_keys <- setdiff(assignment$sample_key, rownames(feature_raw))
if (length(missing_keys) > 0) {
  stop("Feature matrix is missing ", length(missing_keys), " sample_key rows. Example: ", paste(head(missing_keys, 5), collapse = ", "))
}

## Ensure same order and force all evidence columns to numeric after pan-cancer merge.
feature_raw <- feature_raw[assignment$sample_key, , drop = FALSE]
rownames(feature_raw) <- assignment$sample_key
for (cc in colnames(feature_raw)) {
  if (!is.numeric(feature_raw[[cc]])) {
    suppressWarnings(feature_raw[[cc]] <- as.numeric(feature_raw[[cc]]))
  }
}

## Remove only features that are effectively all missing by default.
## Do NOT drop zero-variance features unless explicitly requested: if every evidence
## feature is constant after merging, dropping by variance creates an empty matrix and
## an empty PDF. Constant features are retained and z-scaled to 0.
feature_summary <- data.frame(
  feature = colnames(feature_raw),
  n_finite = vapply(feature_raw, function(x) sum(is.finite(as.numeric(x))), numeric(1)),
  sd = vapply(feature_raw, function(x) sd(as.numeric(x), na.rm = TRUE), numeric(1)),
  stringsAsFactors = FALSE
)
feature_summary$sd[!is.finite(feature_summary$sd)] <- 0
fwrite(feature_summary, file.path(OUT_DIR, "TCGA_noTME2_feature_numeric_summary_before_filter.tsv"), sep = "	")
if (sum(feature_summary$n_finite > 0) == 0) {
  debug_n <- min(20, nrow(assignment), nrow(feature_raw))
  fwrite(
    data.frame(
      assignment_sample_key_head = head(assignment$sample_key, debug_n),
      feature_rowname_head = head(rownames(feature_raw), debug_n),
      stringsAsFactors = FALSE
    ),
    file.path(OUT_DIR, "TCGA_noTME2_empty_numeric_debug_sample_keys.tsv"),
    sep = "	"
  )
}

if (DROP_ZERO_VARIANCE_FEATURES) {
  keep_feature <- feature_summary$n_finite >= MIN_FINITE_PER_FEATURE & feature_summary$sd > 1e-8
} else {
  keep_feature <- feature_summary$n_finite >= MIN_FINITE_PER_FEATURE
}

if (!any(keep_feature)) {
  warning("No feature passed the default finite/variance filter. Keeping all features with at least one finite value.")
  keep_feature <- feature_summary$n_finite > 0
}
if (!any(keep_feature)) {
  stop("No numeric evidence feature remained after filtering. Check input files and sample orientation.")
}

feature_raw <- feature_raw[, keep_feature, drop = FALSE]
message("Evidence features kept for heatmap: ", ncol(feature_raw), " / ", nrow(feature_summary))

feature_z <- zscale_columns(feature_raw)
feature_z_cap <- cap_matrix(feature_z, HEATMAP_ZLIM)

feature_group <- make_feature_groups(colnames(feature_z_cap))
feature_group_tbl <- data.frame(
  feature = colnames(feature_z_cap),
  feature_group = as.character(feature_group),
  stringsAsFactors = FALSE
)

## Order samples by assigned TME.  All cancer types are mixed in the same heatmap,
## with cancer type shown as column annotation rather than separated into facets.
assignment$predicted_TME <- factor(as.character(assignment$predicted_TME), levels = TME_KEEP)
if (ORDER_WITHIN_TME == "margin") {
  ord <- order(assignment$predicted_TME, -assignment$margin, assignment$cancer_type, assignment$sample, na.last = TRUE)
} else if (ORDER_WITHIN_TME == "cancer") {
  ord <- order(assignment$predicted_TME, assignment$cancer_type, assignment$sample, na.last = TRUE)
} else {
  ord <- order(assignment$predicted_TME, assignment$cancer_type, -assignment$margin, assignment$sample, na.last = TRUE)
}
assignment <- assignment[ord, , drop = FALSE]

## IMPORTANT: after pan-cancer merging, rownames(feature_raw/feature_z_cap) are
## project-prefixed sample_key values (e.g. BLCA__TCGA-...).  Using assignment$sample
## can be out-of-range / all-NA because TCGA sample IDs are duplicated across contexts
## and are not the row identifiers anymore.
feature_raw <- feature_raw[assignment$sample_key, , drop = FALSE]
feature_z_cap <- feature_z_cap[assignment$sample_key, , drop = FALSE]

if (!identical(rownames(feature_raw), assignment$sample_key) || !identical(rownames(feature_z_cap), assignment$sample_key)) {
  stop("Internal ordering error: feature matrix rows are not aligned to assignment$sample_key after sorting.")
}

## Optional downsample for preview heatmap only. Full matrices are always saved.
plot_samples <- assignment$sample_key
if (MAX_SAMPLES_PER_TME > 0) {
  plot_samples <- assignment %>%
    group_by(predicted_TME) %>%
    arrange(desc(margin), .by_group = TRUE) %>%
    slice_head(n = MAX_SAMPLES_PER_TME) %>%
    ungroup() %>%
    pull(sample_key)
}

plot_assignment <- assignment[match(plot_samples, assignment$sample_key), , drop = FALSE]
plot_mat <- feature_z_cap[plot_samples, , drop = FALSE]

## ----------------------------
## Save matrices
## ----------------------------

fwrite(assignment, file.path(OUT_DIR, "TCGA_noTME2_all_projects_assignment_for_heatmap.tsv"), sep = "\t")
fwrite(data.frame(sample = rownames(feature_raw), feature_raw, check.names = FALSE),
       file.path(OUT_DIR, "TCGA_noTME2_all_projects_patient_by_feature_matrix_raw.tsv"), sep = "\t")
fwrite(data.frame(sample = rownames(feature_z_cap), feature_z_cap, check.names = FALSE),
       file.path(OUT_DIR, "TCGA_noTME2_all_projects_patient_by_feature_matrix_z_capped.tsv"), sep = "\t")
fwrite(feature_group_tbl, file.path(OUT_DIR, "TCGA_noTME2_feature_groups.tsv"), sep = "\t")

counts_tbl <- assignment %>%
  count(cancer_type, predicted_TME, name = "n") %>%
  arrange(predicted_TME, cancer_type)
fwrite(counts_tbl, file.path(OUT_DIR, "TCGA_noTME2_assignment_counts_by_cancer.tsv"), sep = "\t")

counts_tme <- assignment %>% count(predicted_TME, name = "n")
fwrite(counts_tme, file.path(OUT_DIR, "TCGA_noTME2_assignment_counts_by_TME.tsv"), sep = "\t")

for (tm in TME_KEEP) {
  smp <- assignment$sample_key[assignment$predicted_TME == tm]
  if (length(smp) == 0) next
  ann_sub <- assignment[match(smp, assignment$sample_key), , drop = FALSE]
  meta_sub <- ann_sub[, c("sample_key", "sample", "cancer_type", "predicted_TME", "margin"), drop = FALSE]
  fwrite(ann_sub,
         file.path(OUT_DIR, paste0("TCGA_noTME2_", tm, "_assignment.tsv")), sep = "	")
  fwrite(data.frame(meta_sub, feature_raw[smp, , drop = FALSE], check.names = FALSE),
         file.path(OUT_DIR, paste0("TCGA_noTME2_", tm, "_patient_by_feature_matrix_raw.tsv")), sep = "	")
  fwrite(data.frame(meta_sub, feature_z_cap[smp, , drop = FALSE], check.names = FALSE),
         file.path(OUT_DIR, paste0("TCGA_noTME2_", tm, "_patient_by_feature_matrix_z_capped.tsv")), sep = "	")
}

## ----------------------------
## Heatmap functions
## ----------------------------

plot_complex_heatmap <- function(mat_sample_by_feature, ann, out_pdf, title,
                                 split_by_tme = TRUE, width = 16, height = 9) {
  message("  heatmap matrix: ", nrow(mat_sample_by_feature), " samples x ", ncol(mat_sample_by_feature), " features -> ", out_pdf)
  if (nrow(mat_sample_by_feature) == 0 || ncol(mat_sample_by_feature) == 0) {
    warning("Skip heatmap because matrix has zero rows or columns: ", out_pdf)
    return(invisible(NULL))
  }

  ## ComplexHeatmap draws features as rows and samples as columns.
  mat <- t(as.matrix(mat_sample_by_feature))
  storage.mode(mat) <- "numeric"
  mat[!is.finite(mat)] <- 0

  ## Match sample annotations to heatmap columns. In pan-cancer plots,
  ## columns are keyed by sample_key; sample alone is not safe after rbind.
  key_col <- if ("sample_key" %in% colnames(ann)) "sample_key" else "sample"
  ann <- ann[match(colnames(mat), ann[[key_col]]), , drop = FALSE]
  if (any(is.na(ann[[key_col]]))) {
    stop("Could not match all heatmap columns to annotation rows for: ", out_pdf)
  }

  ## Define feature groups directly from the current matrix. This avoids relying on
  ## a global feature_group vector and prevents NULL row_order errors when the
  ## selected matrix contains only a subset of features.
  fg <- make_feature_groups(rownames(mat))
  fg_chr <- as.character(fg)
  fg_chr[is.na(fg_chr) | fg_chr == ""] <- "Other"
  fg <- factor(fg_chr, levels = levels(make_feature_groups(character(0))))
  fg <- droplevels(fg)

  ## Row direction corresponds to evidence variables/features.
  ## v8: first cluster feature rows, then order the feature clusters by the TME
  ## in which each cluster is most enriched.  This gives rows ordered as
  ## TME1-like clusters, TME3-like clusters, TME4-like clusters, etc.
  row_cluster <- factor(rep("all_features", nrow(mat)), levels = "all_features")
  row_cluster_tme <- factor(rep("Other", nrow(mat)), levels = c(ROW_CLUSTER_ORDER_TME, "Other"))
  row_dend <- FALSE

  out_prefix <- sub("\\.pdf$", "", out_pdf)
  if (ORDER_ROW_CLUSTERS_BY_TME) {
    ord_res <- order_rows_by_tme_characteristic_clusters(mat, ann, out_prefix = out_prefix)
    mat <- ord_res$mat
    row_cluster <- ord_res$row_cluster
    row_cluster_tme <- ord_res$row_cluster_tme
    row_dend <- FALSE
  } else if (CLUSTER_ROWS && nrow(mat) >= 2) {
    row_sd <- apply(mat, 1, sd, na.rm = TRUE)
    if (sum(is.finite(row_sd) & row_sd > 1e-8) >= 2) {
      row_dist <- stats::dist(mat)
      row_hclust <- stats::hclust(row_dist, method = ROW_CLUSTER_METHOD)
      row_dend <- stats::as.dendrogram(row_hclust)
    } else {
      row_ord <- order(fg, seq_along(fg), na.last = TRUE)
      mat <- mat[row_ord, , drop = FALSE]
      fg <- droplevels(fg[row_ord])
      row_dend <- FALSE
    }
  } else {
    row_ord <- order(fg, seq_along(fg), na.last = TRUE)
    mat <- mat[row_ord, , drop = FALSE]
    fg <- droplevels(fg[row_ord])
    row_dend <- FALSE
  }

  ## Recompute feature groups after possible row reordering.
  fg <- make_feature_groups(rownames(mat)) |> droplevels()

  if (requireNamespace("ComplexHeatmap", quietly = TRUE) && requireNamespace("circlize", quietly = TRUE)) {
    ## Do not assign to ComplexHeatmap::ht_opt$message directly; assignment to a
    ## namespace-qualified object can fail with: object 'ComplexHeatmap' not found.
    ## getExportedValue() returns the exported GlobalOptions object safely.
    try({
      ht_opt_obj <- getExportedValue("ComplexHeatmap", "ht_opt")
      ht_opt_obj$message <- FALSE
    }, silent = TRUE)

    tme_levels <- TME_KEEP[TME_KEEP %in% as.character(unique(ann$predicted_TME))]
    ann$predicted_TME <- factor(as.character(ann$predicted_TME), levels = tme_levels)
    ann$cancer_type <- factor(ann$cancer_type, levels = sort(unique(ann$cancer_type)))

    tme_palette <- c("#D9D9D9", "#80B1D3", "#FB8072", "#FDB462", "#B3DE69")
    tme_cols <- setNames(tme_palette[seq_along(tme_levels)], tme_levels)
    cancer_levels <- levels(ann$cancer_type)
    cancer_cols <- setNames(grDevices::hcl.colors(length(cancer_levels), palette = "Dark 3"), cancer_levels)

    ## Numeric margin annotation is omitted when all values are missing.
    if ("margin" %in% colnames(ann) && any(is.finite(ann$margin))) {
      ha <- ComplexHeatmap::HeatmapAnnotation(
        TME = ann$predicted_TME,
        cancer = ann$cancer_type,
        margin = ann$margin,
        col = list(TME = tme_cols, cancer = cancer_cols),
        annotation_name_side = "left",
        show_annotation_name = TRUE
      )
    } else {
      ha <- ComplexHeatmap::HeatmapAnnotation(
        TME = ann$predicted_TME,
        cancer = ann$cancer_type,
        col = list(TME = tme_cols, cancer = cancer_cols),
        annotation_name_side = "left",
        show_annotation_name = TRUE
      )
    }

    row_ha <- NULL
    row_ann_args <- list(show_annotation_name = TRUE, annotation_name_side = "top")
    row_ann_cols <- list()

    if (ORDER_ROW_CLUSTERS_BY_TME && length(row_cluster_tme) == nrow(mat)) {
      row_ann_args$ClusterTME <- row_cluster_tme
      cluster_tme_cols <- c(tme_cols, Other = "#BDBDBD")
      row_ann_cols$ClusterTME <- cluster_tme_cols[names(cluster_tme_cols) %in% levels(droplevels(row_cluster_tme))]
    }
    if (SHOW_FEATURE_GROUP_ROW_ANNOTATION && length(unique(fg)) >= 2) {
      fg_levels <- levels(droplevels(fg))
      row_ann_args$FeatureGroup <- fg
      row_ann_cols$FeatureGroup <- setNames(grDevices::hcl.colors(length(fg_levels), palette = "Set 3"), fg_levels)
    }
    if (length(row_ann_args) > 2) {
      row_ann_args$col <- row_ann_cols
      row_ha <- do.call(ComplexHeatmap::rowAnnotation, row_ann_args)
    }

    ## Columns are samples. Keep cluster_columns=FALSE so that, within each TME,
    ## samples remain ordered by margin rather than by cancer type or clustering.
    row_split <- NULL
    if (ORDER_ROW_CLUSTERS_BY_TME && SPLIT_ROWS_BY_TME_CLUSTER && length(row_cluster) == nrow(mat) && length(unique(row_cluster)) >= 2) {
      row_split <- row_cluster
    }
    column_split <- NULL
    if (split_by_tme && length(tme_levels) >= 2) {
      column_split <- ann$predicted_TME
    }

    ht <- ComplexHeatmap::Heatmap(
      mat,
      name = "z",
      col = circlize::colorRamp2(c(-HEATMAP_ZLIM, 0, HEATMAP_ZLIM), c("#2166AC", "#F7F7F7", "#B2182B")),
      top_annotation = ha,
      row_split = row_split,
      column_split = column_split,
      cluster_rows = if (ORDER_ROW_CLUSTERS_BY_TME) FALSE else row_dend,
      show_row_dend = SHOW_ROW_DEND && !ORDER_ROW_CLUSTERS_BY_TME && !identical(row_dend, FALSE),
      left_annotation = row_ha,
      cluster_columns = FALSE,
      cluster_row_slices = FALSE,
      cluster_column_slices = FALSE,
      show_column_names = FALSE,
      show_row_names = TRUE,
      row_names_gp = grid::gpar(fontsize = 7),
      column_title = title,
      use_raster = FORCE_RASTER || ncol(mat) > RASTER_COLUMN_THRESHOLD,
      raster_device = "png",
      raster_quality = 2,
      heatmap_legend_param = list(title = paste0("z capped\n±", HEATMAP_ZLIM))
    )

    pdf(out_pdf, width = width, height = height, useDingbats = FALSE)
    ComplexHeatmap::draw(ht)
    dev.off()

    png(sub("\\.pdf$", ".png", out_pdf), width = width, height = height, units = "in", res = 220)
    ComplexHeatmap::draw(ht)
    dev.off()

  } else if (requireNamespace("pheatmap", quietly = TRUE)) {
    annot <- data.frame(
      TME = ann$predicted_TME,
      cancer_type = ann$cancer_type,
      row.names = if ("sample_key" %in% colnames(ann)) ann$sample_key else ann$sample
    )
    pdf(out_pdf, width = width, height = height, useDingbats = FALSE)
    pheatmap::pheatmap(
      mat,
      color = colorRampPalette(c("#2166AC", "#F7F7F7", "#B2182B"))(101),
      cluster_rows = CLUSTER_ROWS,
      cluster_cols = FALSE,
      show_colnames = FALSE,
      show_rownames = TRUE,
      annotation_col = annot,
      fontsize_row = 7,
      main = title,
      border_color = NA
    )
    dev.off()
  } else {
    message("Neither ComplexHeatmap nor pheatmap is installed. Skipping heatmap: ", out_pdf)
  }
}

plot_integrated_only <- function(mat_sample_by_feature, ann, out_pdf, title) {
  cols <- colnames(mat_sample_by_feature)[grepl("^Score_TME", colnames(mat_sample_by_feature))]
  if (length(cols) < 2) return(invisible(NULL))
  submat <- mat_sample_by_feature[, cols, drop = FALSE]
  plot_complex_heatmap(submat, ann, out_pdf, title, split_by_tme = TRUE, width = 15, height = 4.8)
}

## ----------------------------
## Make heatmaps
## ----------------------------

message("Plotting combined evidence heatmap...")
plot_complex_heatmap(
  plot_mat,
  plot_assignment,
  file.path(OUT_DIR, "TCGA_noTME2_all_projects_assignment_evidence_heatmap_by_TME.pdf"),
  title = "TCGA pan-cancer no-TME2 assignment evidence matrix",
  split_by_tme = TRUE,
  width = ifelse(nrow(plot_assignment) > 2000, 24, 18),
  height = 10
)

message("Plotting integrated TME-score-only heatmap...")
plot_integrated_only(
  plot_mat,
  plot_assignment,
  file.path(OUT_DIR, "TCGA_noTME2_all_projects_integrated_TME_scores_heatmap_by_TME.pdf"),
  title = "TCGA pan-cancer integrated TME scores used for assignment"
)

for (tm in TME_KEEP) {
  smp <- assignment$sample_key[assignment$predicted_TME == tm]
  if (length(smp) < 2) next
  if (MAX_SAMPLES_PER_TME > 0 && length(smp) > MAX_SAMPLES_PER_TME) {
    smp <- assignment %>%
      filter(predicted_TME == tm) %>%
      arrange(desc(margin)) %>%
      slice_head(n = MAX_SAMPLES_PER_TME) %>%
      pull(sample_key)
  }
  ann_tm <- assignment[match(smp, assignment$sample_key), , drop = FALSE]
  mat_tm <- feature_z_cap[smp, , drop = FALSE]
  message("Plotting per-TME heatmap: ", tm, " n=", length(smp))
  plot_complex_heatmap(
    mat_tm,
    ann_tm,
    file.path(OUT_DIR, paste0("TCGA_noTME2_", tm, "_assignment_evidence_heatmap.pdf")),
    title = paste0("TCGA pan-cancer ", tm, " assignment evidence matrix"),
    split_by_tme = FALSE,
    width = ifelse(length(smp) > 1200, 20, 14),
    height = 9
  )
}

message("==============================")
message("Done.")
message("Output directory: ", OUT_DIR)
message("Key outputs:")
message("  ", file.path(OUT_DIR, "TCGA_noTME2_all_projects_patient_by_feature_matrix_raw.tsv"))
message("  ", file.path(OUT_DIR, "TCGA_noTME2_all_projects_patient_by_feature_matrix_z_capped.tsv"))
message("  ", file.path(OUT_DIR, "TCGA_noTME2_all_projects_assignment_evidence_heatmap_by_TME.pdf"))
message("  ", file.path(OUT_DIR, "TCGA_noTME2_all_projects_integrated_TME_scores_heatmap_by_TME.pdf"))
