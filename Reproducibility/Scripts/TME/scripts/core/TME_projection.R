#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(data.table)
  library(dplyr)
})

## ============================================================
## TCGA TME module / gene signature / DOGMA-template projection
##
## Goal:
##   Use Scaden as a broad-module deconvolution input, not as a
##   single-cell-state oracle. Combine:
##     1) Scaden-derived robust module scores
##     2) TCGA bulk gene signature scores
##     3) DOGMA TME template projection
##   Then rank-normalize the resulting module/TME scores.
## ============================================================

get_env <- function(x, default = NULL) {
  v <- Sys.getenv(x, unset = NA_character_)
  if (is.na(v) || v == "") return(default)
  v
}

PROJECT <- get_env("PROJECT", "BLCA")
RUN_TAG <- get_env("RUN_TAG", "v15j_hard_zero_dys_rank_orthogonal_ns40000_steps40000")
MODEL <- get_env("MODEL", "m256")

BASE_DIR <- get_env("SCADEN_BASE_DIR", "/home/k-makino/wd_home/UC_DOGMA_reseq/scaden/TCGA_Pancancer_TME")
RUN_DIR <- file.path(BASE_DIR, PROJECT, RUN_TAG)

PRED_FILE <- get_env(
  "PRED_FILE",
  file.path(RUN_DIR, "pred", paste0("TCGA_", PROJECT, "_Scaden_", MODEL, "_predictions.txt"))
)
if (!file.exists(PRED_FILE)) {
  alt <- file.path(RUN_DIR, "pred", paste0("TCGA_", PROJECT, "_Scaden_", MODEL, "_sample_by_state.csv"))
  if (file.exists(alt)) PRED_FILE <- alt
}

BULK_FILE <- get_env("BULK_FILE", file.path(RUN_DIR, "bulk", "bulk.tsv"))
DOGMA_ABUNDANCE_FILE <- get_env(
  "DOGMA_ABUNDANCE_FILE",
  "/home/k-makino/wd_home/UC_DOGMA_reseq/ECOTYPE/cell_abundance_w_epithelial.txt"
)
DOGMA_TME_ASSIGNMENT_FILE <- get_env(
  "DOGMA_TME_ASSIGNMENT_FILE",
  "/home/k-makino/wd_home/UC_DOGMA_reseq/ECOTYPE/TME_assignment_DOGMA.txt"
)

## Projection target TME classes.
## Default excludes TME2 because TME2/TLS-like is considered a BCG-induced,
## treatment-specific state and may not be suitable for pan-cancer / external
## pretreatment cohort projection.
TME_KEEP <- unlist(strsplit(get_env("TME_KEEP", "TME1,TME3,TME4,TME5,TME6"), ","))
TME_KEEP <- trimws(TME_KEEP)
TME_KEEP <- TME_KEEP[nzchar(TME_KEEP)]
if (length(TME_KEEP) < 2) stop("TME_KEEP must contain at least two TME labels.")

OUT_DIR <- get_env("OUT_DIR", file.path(RUN_DIR, "module_signature_template_projection_noTME2"))
dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)

message("PROJECT: ", PROJECT)
message("RUN_TAG: ", RUN_TAG)
message("PRED_FILE: ", PRED_FILE)
message("BULK_FILE: ", BULK_FILE)
message("DOGMA_ABUNDANCE_FILE: ", DOGMA_ABUNDANCE_FILE)
message("DOGMA_TME_ASSIGNMENT_FILE: ", DOGMA_TME_ASSIGNMENT_FILE)
message("TME_KEEP: ", paste(TME_KEEP, collapse = ","))
message("OUT_DIR: ", OUT_DIR)

stopifnot(file.exists(PRED_FILE))
stopifnot(file.exists(BULK_FILE))
stopifnot(file.exists(DOGMA_ABUNDANCE_FILE))
stopifnot(file.exists(DOGMA_TME_ASSIGNMENT_FILE))

## -----------------------------
## Utility functions
## -----------------------------

clean_state <- function(x) {
  x <- as.character(x)
  x <- gsub("-", "_", x)
  x <- gsub("\\.", "_", x)
  x <- gsub(" ", "_", x)
  x
}

rank_norm_vec <- function(x) {
  x <- as.numeric(x)
  ok <- is.finite(x)
  y <- rep(NA_real_, length(x))
  if (sum(ok) <= 1) {
    y[ok] <- 0
    return(y)
  }
  r <- rank(x[ok], ties.method = "average", na.last = "keep")
  y[ok] <- qnorm((r - 0.5) / length(r))
  y
}

rank_norm_rows <- function(mat) {
  mat <- as.matrix(mat)
  out <- t(apply(mat, 1, rank_norm_vec))
  rownames(out) <- rownames(mat)
  colnames(out) <- colnames(mat)
  out
}

cohort_robust_z <- function(mat, eps = 1e-8, log_transform = TRUE) {
  mat <- as.matrix(mat)
  x <- mat
  if (log_transform) x <- log10(pmax(x, 0) + eps)
  med <- apply(x, 2, median, na.rm = TRUE)
  madv <- apply(x, 2, mad, na.rm = TRUE)
  madv[!is.finite(madv) | madv < 1e-8] <- 1
  z <- sweep(x, 2, med, "-")
  z <- sweep(z, 2, madv, "/")
  z
}


normalize_tcga_sample_id <- function(x) {
  x <- as.character(x)
  ## Collapse aliquot/vial-level TCGA barcodes to tumor sample-type level, e.g.
  ## TCGA-BL-A0C8-01A / TCGA-BL-A0C8-01B -> TCGA-BL-A0C8-01.
  is_tcga <- grepl("^TCGA-[A-Za-z0-9]{2}-[A-Za-z0-9]{4}-[0-9]{2}", x)
  x[is_tcga] <- substr(x[is_tcga], 1, 15)
  x
}

collapse_matrix_by_id <- function(m, ids, label = "matrix") {
  m <- as.matrix(m)
  ids <- normalize_tcga_sample_id(ids)
  if (length(ids) != nrow(m)) stop("ids length does not match nrow for ", label)
  if (anyDuplicated(ids)) {
    dup_ids <- unique(ids[duplicated(ids)])
    message("Collapsing duplicated sample IDs in ", label, ": ", length(dup_ids), " duplicated IDs")
    rs <- rowsum(m, group = ids, reorder = FALSE)
    n <- as.numeric(table(factor(ids, levels = rownames(rs))))
    rs <- sweep(rs, 1, n, "/")
    return(rs)
  }
  rownames(m) <- ids
  m
}

collapse_df_by_id <- function(df, ids, label = "data.frame") {
  m <- as.matrix(df)
  storage.mode(m) <- "numeric"
  m <- collapse_matrix_by_id(m, ids, label = label)
  out <- as.data.frame(m, check.names = FALSE)
  rownames(out) <- rownames(m)
  out
}

safe_sum <- function(df, states, weights = NULL) {
  states <- clean_state(states)
  states <- intersect(states, colnames(df))
  if (length(states) == 0) return(rep(0, nrow(df)))
  X <- as.matrix(df[, states, drop = FALSE])
  if (is.null(weights)) {
    rowSums(X, na.rm = TRUE)
  } else {
    w <- weights[states]
    w[is.na(w)] <- 1
    as.numeric(X %*% as.numeric(w))
  }
}

read_scaden_prediction <- function(path) {
  dt <- fread(path)
  df <- as.data.frame(dt)
  ## fread on tab files with empty first header often calls first column V1.
  row_ids <- normalize_tcga_sample_id(df[[1]])
  df <- df[, -1, drop = FALSE]
  colnames(df) <- clean_state(colnames(df))
  for (cc in colnames(df)) df[[cc]] <- as.numeric(df[[cc]])
  df <- collapse_df_by_id(df, row_ids, label = "Scaden prediction")
  df
}

read_bulk_matrix <- function(path) {
  dt <- fread(path)
  df <- as.data.frame(dt, check.names = FALSE)
  first_col <- df[[1]]
  df <- df[, -1, drop = FALSE]
  rownames(df) <- make.unique(as.character(first_col))
  m <- as.matrix(df)
  storage.mode(m) <- "numeric"

  ## Guess orientation and return samples x genes.
  ## TCGA bulk.tsv usually has genes x TCGA-samples.
  ## External cohorts may have non-TCGA sample names, so TCGA barcode detection
  ## is insufficient.  In that case, use signature-gene enrichment in rownames
  ## versus colnames to detect genes x samples.
  sample_like_rows <- mean(grepl("^TCGA", rownames(m)))
  sample_like_cols <- mean(grepl("^TCGA", colnames(m)))

  marker_genes_for_orientation <- unique(toupper(c(
    "MS4A1", "CD79A", "CD79B", "CD19", "BANK1", "MZB1", "JCHAIN",
    "CD3D", "CD3E", "CD8A", "CD8B", "CD4", "CCR7", "IL7R", "SELL",
    "TCF7", "LEF1", "NKG7", "PRF1", "GZMB", "GZMH", "GNLY", "IFNG",
    "PDCD1", "LAG3", "HAVCR2", "TIGIT", "TOX", "CXCL13", "FOXP3",
    "IL2RA", "CTLA4", "COL1A1", "COL1A2", "COL3A1", "DCN", "LUM",
    "ACTA2", "TAGLN", "RGS5", "MCAM", "PDGFRB", "LST1", "C1QA",
    "C1QB", "APOE", "TREM2", "KRT5", "KRT14", "KRT20", "EPCAM"
  )))
  row_gene_hits <- sum(toupper(rownames(m)) %in% marker_genes_for_orientation)
  col_gene_hits <- sum(toupper(colnames(m)) %in% marker_genes_for_orientation)

  transpose_to_samples_by_genes <- FALSE
  if (is.finite(sample_like_cols) && sample_like_cols > sample_like_rows) {
    transpose_to_samples_by_genes <- TRUE
  } else if (row_gene_hits > col_gene_hits) {
    transpose_to_samples_by_genes <- TRUE
  }

  if (transpose_to_samples_by_genes) {
    message("Detected bulk matrix as genes x samples; transposing to samples x genes for signature scoring.")
    m <- t(m)
  } else {
    message("Detected bulk matrix as samples x genes; using as-is for signature scoring.")
  }

  ## Collapse aliquot/vial duplicates after standardizing to sample-type barcode.
  m <- collapse_matrix_by_id(m, rownames(m), label = "bulk expression")
  m
}

logcpm <- function(counts) {
  counts <- as.matrix(counts)
  counts[counts < 0 | !is.finite(counts)] <- 0
  lib <- rowSums(counts)
  lib[lib <= 0] <- 1
  log2(sweep(counts, 1, lib, "/") * 1e6 + 1)
}

score_signature <- function(expr_log, genes) {
  genes <- intersect(unique(genes), colnames(expr_log))
  if (length(genes) < 2) {
    return(rep(NA_real_, nrow(expr_log)))
  }
  X <- expr_log[, genes, drop = FALSE]
  Z <- scale(X)
  rowMeans(Z, na.rm = TRUE)
}

parse_tme_assignment <- function(path) {
  txt <- readLines(path, warn = FALSE)
  ## Case 1: R-style named vector: BC_029 = "TME1"
  m <- regexec("([A-Za-z0-9_\\-]+)\\s*=\\s*['\"]?(TME[0-9]+)['\"]?", txt)
  rr <- regmatches(txt, m)
  rr <- rr[lengths(rr) == 3]
  if (length(rr) > 0) {
    out <- data.frame(
      sample = vapply(rr, `[`, character(1), 2),
      TME = vapply(rr, `[`, character(1), 3),
      stringsAsFactors = FALSE
    )
    return(out)
  }
  ## Case 2: table.
  dt <- fread(path)
  df <- as.data.frame(dt)
  if (ncol(df) < 2) stop("Could not parse TME assignment file: ", path)
  tme_col <- which(sapply(df, function(x) any(grepl("^TME[0-9]+$", as.character(x)))))
  if (length(tme_col) == 0) stop("Could not find TME column in: ", path)
  tme_col <- tme_col[1]
  sample_col <- setdiff(seq_len(ncol(df)), tme_col)[1]
  out <- data.frame(
    sample = as.character(df[[sample_col]]),
    TME = as.character(df[[tme_col]]),
    stringsAsFactors = FALSE
  )
  out <- out[grepl("^TME[0-9]+$", out$TME), ]
  out
}

read_dogma_abundance <- function(path) {
  dt <- fread(path)
  df <- as.data.frame(dt)
  row_ids <- df[[1]]
  df <- df[, -1, drop = FALSE]
  rownames(df) <- clean_state(row_ids)
  for (cc in colnames(df)) df[[cc]] <- as.numeric(df[[cc]])
  m <- as.matrix(df)
  ## Output samples x states.
  m <- t(m)
  colnames(m) <- clean_state(colnames(m))
  m
}

## -----------------------------
## Modules
## -----------------------------

module_list <- list(
  Epithelial = c("Epithelial"),
  TLS_B = c("Atypical_B", "B_memory", "B_naive", "GC_B", "Plasma", "CD4_Tfh_like"),
  Naive_memory = c("CD4_Tn", "CD4_Tcm", "CD4_Tsen", "CD8_Tn", "CD8_Tcm", "Treg_naive", "B_naive", "NK_CD56_CD49a_Lo"),
  Cytotoxic_Th17 = c("CD8_Tem", "CD4_Th17"),
  Dysfunctional = c("CD8_Tex", "Treg_effector", "CD4_Tph_like"),
  MDSC_myeloid_suppressive = c("MDSC_like"),
  Myeloid_TAM = c("Mono", "TAM_FOLR2", "TAM_TREM2", "cDC2", "cDC1", "preDC", "mregDC", "pDC"),
  Stroma_CAF = c("matCAF", "iCAF", "vSMC", "contCAP", "Endothelial"),
  NK_ILC = c("NK_CD56_dim", "NK_CD56_CD49a_Hi", "NK_CD56_CD49a_Lo", "ILC3"),
  Minor_T = c("MAIT", "CD4_CTL", "CD8_Temra", "CD8_Trm")
)

compute_scaden_modules <- function(theta) {
  theta <- as.data.frame(theta)
  colnames(theta) <- clean_state(colnames(theta))
  out <- data.frame(row.names = rownames(theta))
  for (nm in names(module_list)) {
    out[[nm]] <- safe_sum(theta, module_list[[nm]])
  }
  ## Add key ratios. These are often more stable than raw fine-state fractions.
  eps <- 1e-8
  t_total <- safe_sum(theta, c("CD4_Tn", "CD4_Tcm", "CD4_Th17", "CD4_Tfh_like", "CD4_Tph_like", "CD4_Tsen", "CD8_Tn", "CD8_Tcm", "CD8_Tem", "CD8_Tex", "CD8_Trm", "CD8_Temra", "Treg_naive", "Treg_effector", "MAIT", "CD4_CTL"))
  cd8_total <- safe_sum(theta, c("CD8_Tn", "CD8_Tcm", "CD8_Tem", "CD8_Tex", "CD8_Trm", "CD8_Temra"))
  myeloid_total <- safe_sum(theta, c("Mono", "TAM_FOLR2", "TAM_TREM2", "MDSC_like", "cDC1", "cDC2", "preDC", "mregDC", "pDC"))
  out$CD8_Tem_over_CD8 <- safe_sum(theta, "CD8_Tem") / (cd8_total + eps)
  out$CD8_Tex_over_CD8 <- safe_sum(theta, "CD8_Tex") / (cd8_total + eps)
  out$Dysfunctional_over_T <- out$Dysfunctional / (t_total + eps)
  out$MDSC_over_myeloid <- safe_sum(theta, "MDSC_like") / (myeloid_total + eps)
  out
}

## -----------------------------
## Gene signatures
## -----------------------------

signature_genes <- list(
  Sig_cytotoxic = c("NKG7", "PRF1", "GZMB", "GZMH", "GNLY", "IFNG", "CCL5", "CXCL9", "CXCL10"),
  Sig_exhaustion = c("PDCD1", "LAG3", "HAVCR2", "TIGIT", "TOX", "ENTPD1", "CXCL13"),
  Sig_Treg = c("FOXP3", "IL2RA", "CTLA4", "IKZF2", "TNFRSF18", "TNFRSF4"),
  Sig_Tph_Tfh = c("CXCL13", "PDCD1", "ICOS", "MAF", "TOX", "BCL6"),
  Sig_TLS_B = c("MS4A1", "CD79A", "CD79B", "CD19", "CD37", "BANK1", "MZB1", "JCHAIN", "LTB", "CXCL13"),
  Sig_naive_T = c("CCR7", "IL7R", "SELL", "TCF7", "LEF1", "LTB"),
  Sig_CAF_stroma = c("COL1A1", "COL1A2", "COL3A1", "DCN", "LUM", "COL6A1", "ACTA2", "TAGLN", "RGS5", "MCAM", "PDGFRB"),
  Sig_TAM_myeloid = c("LST1", "C1QA", "C1QB", "C1QC", "APOE", "FCGR3A", "TYROBP", "AIF1", "TREM2"),
  Sig_epithelial = c("EPCAM", "KRT8", "KRT18", "KRT19", "KRT7"),
  Sig_EMT_pEMT = c("VIM", "FN1", "ITGA5", "ZEB1", "ZEB2", "SNAI2", "TWIST1")
)

compute_signature_scores <- function(bulk_file) {
  counts <- read_bulk_matrix(bulk_file)
  expr <- logcpm(counts)
  genes <- colnames(expr)
  colnames(expr) <- make.unique(toupper(colnames(expr)))
  signature_genes_upper <- lapply(signature_genes, toupper)
  out <- data.frame(row.names = rownames(expr))
  for (nm in names(signature_genes_upper)) {
    out[[nm]] <- score_signature(expr, signature_genes_upper[[nm]])
  }
  out
}

## -----------------------------
## Template projection
## -----------------------------

pearson_sim <- function(x, y) {
  x <- as.numeric(x); y <- as.numeric(y)
  ok <- is.finite(x) & is.finite(y)
  if (sum(ok) < 3) return(NA_real_)
  suppressWarnings(cor(x[ok], y[ok], method = "pearson"))
}

cosine_sim <- function(x, y) {
  x <- as.numeric(x); y <- as.numeric(y)
  ok <- is.finite(x) & is.finite(y)
  if (sum(ok) < 3) return(NA_real_)
  den <- sqrt(sum(x[ok]^2)) * sqrt(sum(y[ok]^2))
  if (!is.finite(den) || den <= 0) return(NA_real_)
  sum(x[ok] * y[ok]) / den
}

make_template_projection <- function(tcga_module_rn, dogma_module_rn, dogma_tme) {
  dogma_tme <- dogma_tme[dogma_tme$sample %in% rownames(dogma_module_rn), ]
  common_modules <- intersect(colnames(tcga_module_rn), colnames(dogma_module_rn))
  stopifnot(length(common_modules) >= 3)
  templates <- lapply(sort(unique(dogma_tme$TME)), function(tme) {
    smp <- dogma_tme$sample[dogma_tme$TME == tme]
    colMeans(dogma_module_rn[smp, common_modules, drop = FALSE], na.rm = TRUE)
  })
  names(templates) <- sort(unique(dogma_tme$TME))
  template_mat <- do.call(rbind, templates)
  colnames(template_mat) <- common_modules

  pear <- matrix(NA_real_, nrow = nrow(tcga_module_rn), ncol = nrow(template_mat),
                 dimnames = list(rownames(tcga_module_rn), paste0("TemplatePearson_", rownames(template_mat))))
  cos <- matrix(NA_real_, nrow = nrow(tcga_module_rn), ncol = nrow(template_mat),
                dimnames = list(rownames(tcga_module_rn), paste0("TemplateCosine_", rownames(template_mat))))
  for (i in seq_len(nrow(tcga_module_rn))) {
    x <- tcga_module_rn[i, common_modules]
    for (j in seq_len(nrow(template_mat))) {
      y <- template_mat[j, common_modules]
      pear[i, j] <- pearson_sim(x, y)
      cos[i, j] <- cosine_sim(x, y)
    }
  }
  list(template_mat = template_mat, pearson = pear, cosine = cos)
}

## -----------------------------
## Integrated TME scores
## -----------------------------

zcols <- function(df) {
  df <- as.data.frame(df)
  out <- as.data.frame(scale(df))
  out[!is.finite(as.matrix(out))] <- 0
  rownames(out) <- rownames(df)
  out
}

build_integrated_scores <- function(sc_mod_z, sig_z, template_pearson) {
  common <- Reduce(intersect, list(rownames(sc_mod_z), rownames(sig_z), rownames(template_pearson)))
  sc <- sc_mod_z[common, , drop = FALSE]
  sg <- sig_z[common, , drop = FALSE]
  tp <- template_pearson[common, , drop = FALSE]
  val <- function(M, nm) if (nm %in% colnames(M)) M[, nm] else rep(0, nrow(M))

  out <- data.frame(row.names = common)
  out$TME1_desert_epithelial <-
    0.8 * val(sc, "Epithelial") + 0.8 * val(sg, "Sig_epithelial") -
    0.4 * val(sc, "TLS_B") - 0.4 * val(sc, "Dysfunctional") - 0.4 * val(sc, "Cytotoxic_Th17") - 0.4 * val(sc, "Stroma_CAF")

  out$TME2_TLS_active <-
    0.8 * val(sc, "TLS_B") + 0.7 * val(sg, "Sig_TLS_B") +
    0.4 * val(sc, "Cytotoxic_Th17") + 0.4 * val(sg, "Sig_cytotoxic") -
    0.7 * val(sc, "Dysfunctional") - 0.5 * val(sc, "Stroma_CAF")

  out$TME3_naive_memory <-
    0.9 * val(sc, "Naive_memory") + 0.8 * val(sg, "Sig_naive_T") +
    0.2 * val(sc, "TLS_B") -
    0.7 * val(sc, "Dysfunctional") - 0.5 * val(sc, "Cytotoxic_Th17") - 0.4 * val(sc, "Stroma_CAF")

  out$TME4_dysfunctional <-
    0.9 * val(sc, "Dysfunctional") + 0.6 * val(sc, "Dysfunctional_over_T") +
    0.7 * val(sg, "Sig_exhaustion") + 0.7 * val(sg, "Sig_Treg") + 0.5 * val(sg, "Sig_Tph_Tfh") -
    0.5 * val(sc, "Cytotoxic_Th17") - 0.4 * val(sc, "Stroma_CAF")

  out$TME5_cytotoxic_Th17 <-
    0.9 * val(sc, "Cytotoxic_Th17") + 0.7 * val(sc, "CD8_Tem_over_CD8") +
    0.8 * val(sg, "Sig_cytotoxic") -
    0.8 * val(sc, "Dysfunctional") - 0.5 * val(sg, "Sig_Treg") - 0.4 * val(sc, "Stroma_CAF")

  out$TME6_stroma <-
    1.0 * val(sc, "Stroma_CAF") + 0.9 * val(sg, "Sig_CAF_stroma") +
    0.2 * val(sc, "Myeloid_TAM") + 0.2 * val(sg, "Sig_TAM_myeloid") -
    0.7 * val(sc, "Dysfunctional") - 0.5 * val(sc, "Cytotoxic_Th17")

  ## Add template similarity as an optional stabilizer.
  for (tme in paste0("TME", 1:6)) {
    cn <- paste0("TemplatePearson_", tme)
    if (cn %in% colnames(tp)) {
      idx <- switch(tme,
                    TME1 = "TME1_desert_epithelial",
                    TME2 = "TME2_TLS_active",
                    TME3 = "TME3_naive_memory",
                    TME4 = "TME4_dysfunctional",
                    TME5 = "TME5_cytotoxic_Th17",
                    TME6 = "TME6_stroma")
      out[[idx]] <- out[[idx]] + 0.5 * val(tp, cn)
    }
  }
  out
}

## -----------------------------
## Main
## -----------------------------

theta <- read_scaden_prediction(PRED_FILE)
scaden_modules_raw <- compute_scaden_modules(theta)
scaden_modules_z <- cohort_robust_z(scaden_modules_raw, log_transform = TRUE)
scaden_modules_rn <- rank_norm_rows(scaden_modules_z)

sig_scores_raw <- compute_signature_scores(BULK_FILE)
common_sig <- intersect(rownames(scaden_modules_raw), rownames(sig_scores_raw))
sig_scores_raw <- sig_scores_raw[common_sig, , drop = FALSE]
sig_scores_z <- zcols(sig_scores_raw)
sig_scores_rn <- rank_norm_rows(sig_scores_z)

## DOGMA templates using the same Scaden-derived module definitions.
dogma_ab <- read_dogma_abundance(DOGMA_ABUNDANCE_FILE)
dogma_tme <- parse_tme_assignment(DOGMA_TME_ASSIGNMENT_FILE)
rownames(dogma_ab) <- normalize_tcga_sample_id(rownames(dogma_ab))
dogma_tme$sample <- normalize_tcga_sample_id(dogma_tme$sample)

## Remove TME classes that should not be used as projection targets.
## This removes them both from the DOGMA template construction and from the
## final integrated score columns used for max-score assignment.
dogma_tme_original <- dogma_tme
dogma_tme <- dogma_tme[dogma_tme$TME %in% TME_KEEP, , drop = FALSE]
if (nrow(dogma_tme) == 0) stop("No DOGMA samples remain after filtering by TME_KEEP: ", paste(TME_KEEP, collapse = ","))
missing_keep <- setdiff(TME_KEEP, unique(dogma_tme$TME))
if (length(missing_keep) > 0) {
  warning("Requested TME_KEEP labels absent from DOGMA templates: ", paste(missing_keep, collapse = ","))
}
message("DOGMA template counts after TME_KEEP filtering:")
print(table(dogma_tme$TME))

dogma_mod_raw <- compute_scaden_modules(as.data.frame(dogma_ab))
dogma_mod_z <- cohort_robust_z(dogma_mod_raw, log_transform = TRUE)
dogma_mod_rn <- rank_norm_rows(dogma_mod_z)
proj <- make_template_projection(scaden_modules_rn, dogma_mod_rn, dogma_tme)
template_pearson <- proj$pearson
template_cosine <- proj$cosine

template_pearson_z <- zcols(as.data.frame(template_pearson))

common_all <- Reduce(intersect, list(rownames(scaden_modules_rn), rownames(sig_scores_rn), rownames(template_pearson_z)))
integrated_raw <- build_integrated_scores(
  zcols(as.data.frame(scaden_modules_z[common_all, , drop = FALSE])),
  zcols(as.data.frame(sig_scores_z[common_all, , drop = FALSE])),
  zcols(as.data.frame(template_pearson[common_all, , drop = FALSE]))
)

## Keep only the requested projection targets.
score_to_tme <- sub("^TME([0-9]).*", "TME\\1", colnames(integrated_raw))
keep_score_cols <- colnames(integrated_raw)[score_to_tme %in% TME_KEEP]
if (length(keep_score_cols) < 2) {
  stop("Fewer than two integrated score columns remain after TME_KEEP filtering: ",
       paste(TME_KEEP, collapse = ","))
}
integrated_raw_all6 <- integrated_raw
integrated_raw <- integrated_raw[, keep_score_cols, drop = FALSE]
integrated_rn <- rank_norm_rows(integrated_raw)

predicted_TME <- colnames(integrated_raw)[max.col(as.matrix(integrated_raw), ties.method = "first")]
predicted_TME <- sub("^TME([0-9]).*", "TME\\1", predicted_TME)
assignment <- data.frame(
  sample = rownames(integrated_raw),
  predicted_TME = predicted_TME,
  max_score = apply(integrated_raw, 1, max, na.rm = TRUE),
  second_score = apply(integrated_raw, 1, function(x) sort(x, decreasing = TRUE)[2]),
  stringsAsFactors = FALSE
)
assignment$margin <- assignment$max_score - assignment$second_score

## Write outputs.
fwrite(data.frame(TME = names(table(dogma_tme_original$TME)), n_original = as.integer(table(dogma_tme_original$TME))),
       file.path(OUT_DIR, "DOGMA_TME_template_counts_original.tsv"), sep = "\t")
fwrite(data.frame(TME = names(table(dogma_tme$TME)), n_used = as.integer(table(dogma_tme$TME))),
       file.path(OUT_DIR, "DOGMA_TME_template_counts_used.tsv"), sep = "\t")
fwrite(data.frame(TME_KEEP = TME_KEEP), file.path(OUT_DIR, "projection_TME_KEEP.tsv"), sep = "\t")

fwrite(data.frame(sample = rownames(scaden_modules_raw), scaden_modules_raw), file.path(OUT_DIR, "scaden_module_scores_raw.tsv"), sep = "\t")
fwrite(data.frame(sample = rownames(scaden_modules_z), scaden_modules_z), file.path(OUT_DIR, "scaden_module_scores_cohort_z.tsv"), sep = "\t")
fwrite(data.frame(sample = rownames(scaden_modules_rn), scaden_modules_rn), file.path(OUT_DIR, "scaden_module_scores_ranknorm.tsv"), sep = "\t")

fwrite(data.frame(sample = rownames(sig_scores_raw), sig_scores_raw), file.path(OUT_DIR, "bulk_gene_signature_scores_raw.tsv"), sep = "\t")
fwrite(data.frame(sample = rownames(sig_scores_z), sig_scores_z), file.path(OUT_DIR, "bulk_gene_signature_scores_z.tsv"), sep = "\t")
fwrite(data.frame(sample = rownames(sig_scores_rn), sig_scores_rn), file.path(OUT_DIR, "bulk_gene_signature_scores_ranknorm.tsv"), sep = "\t")

fwrite(data.frame(TME = rownames(proj$template_mat), proj$template_mat), file.path(OUT_DIR, "DOGMA_TME_module_templates_ranknorm.tsv"), sep = "\t")
fwrite(data.frame(sample = rownames(template_pearson), template_pearson), file.path(OUT_DIR, "TCGA_to_DOGMA_TME_template_pearson.tsv"), sep = "\t")
fwrite(data.frame(sample = rownames(template_cosine), template_cosine), file.path(OUT_DIR, "TCGA_to_DOGMA_TME_template_cosine.tsv"), sep = "\t")

fwrite(data.frame(sample = rownames(integrated_raw_all6), integrated_raw_all6), file.path(OUT_DIR, "integrated_TME_scores_raw_all6_before_TME_KEEP_filter.tsv"), sep = "\t")
fwrite(data.frame(sample = rownames(integrated_raw), integrated_raw), file.path(OUT_DIR, "integrated_TME_scores_raw.tsv"), sep = "\t")
fwrite(data.frame(sample = rownames(integrated_rn), integrated_rn), file.path(OUT_DIR, "integrated_TME_scores_ranknorm.tsv"), sep = "\t")
fwrite(assignment, file.path(OUT_DIR, "integrated_TME_assignment.tsv"), sep = "\t")
fwrite(as.data.frame(table(assignment$predicted_TME)), file.path(OUT_DIR, "integrated_TME_assignment_counts.tsv"), sep = "\t")

## Optional heatmaps.
plot_heatmap <- function(mat, out_pdf, title) {
  mat <- as.matrix(mat)
  if (requireNamespace("ComplexHeatmap", quietly = TRUE) && requireNamespace("circlize", quietly = TRUE)) {
    pdf(out_pdf, width = 12, height = 6)
    ht <- ComplexHeatmap::Heatmap(
      t(mat),
      name = "ranknorm",
      col = circlize::colorRamp2(c(-2, 0, 2), c("#2166AC", "#F7F7F7", "#B2182B")),
      cluster_rows = TRUE,
      cluster_columns = TRUE,
      show_column_names = FALSE,
      column_title = title
    )
    print(ht)
    dev.off()
  } else if (requireNamespace("pheatmap", quietly = TRUE)) {
    pdf(out_pdf, width = 12, height = 6)
    pheatmap::pheatmap(t(mat), show_colnames = FALSE, main = title)
    dev.off()
  } else {
    message("Neither ComplexHeatmap nor pheatmap is installed; skip heatmap: ", out_pdf)
  }
}

plot_heatmap(scaden_modules_rn[common_all, , drop = FALSE], file.path(OUT_DIR, "heatmap_scaden_modules_ranknorm.pdf"), "Scaden robust module rank-normalized heatmap")
plot_heatmap(sig_scores_rn[common_all, , drop = FALSE], file.path(OUT_DIR, "heatmap_gene_signatures_ranknorm.pdf"), "Bulk gene signature rank-normalized heatmap")
plot_heatmap(template_pearson_z[common_all, , drop = FALSE], file.path(OUT_DIR, "heatmap_DOGMA_template_similarity_z.pdf"), "DOGMA TME template similarity")
plot_heatmap(integrated_rn, file.path(OUT_DIR, "heatmap_integrated_TME_scores_ranknorm.pdf"), "Integrated TME score rank-normalized heatmap")

message("Done. Output written to: ", OUT_DIR)
message("Key files:")
message("  ", file.path(OUT_DIR, "integrated_TME_scores_ranknorm.tsv"))
message("  ", file.path(OUT_DIR, "integrated_TME_assignment.tsv"))
message("  ", file.path(OUT_DIR, "heatmap_integrated_TME_scores_ranknorm.pdf"))
