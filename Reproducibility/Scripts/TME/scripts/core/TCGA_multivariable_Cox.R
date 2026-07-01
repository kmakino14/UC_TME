#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(data.table)
  library(dplyr)
  library(tibble)
  library(stringr)
  library(survival)
  library(ggplot2)
})

## ============================================================
## TCGA no-TME2 multivariable Cox analysis
##   1) within each cancer type: Surv(OS.time, OS) ~ TME + gender
##   2) pan-cancer integrated:  Surv(OS.time, OS) ~ TME + cancer_type + gender
##
## TME2/TLS-like is excluded by default; TME5 is used as the reference level.
## ============================================================

BASE_DIR <- Sys.getenv("BASE_DIR", unset = "/home/k-makino/wd_home/UC_DOGMA_reseq/scaden/TCGA_Pancancer_TME")
RUN_TAG <- Sys.getenv("RUN_TAG", unset = "v15k_ranktemplate_tme2_tme3_tme5_recovery_ns40000_steps40000")
CDR_FILE <- Sys.getenv("TCGA_CDR_FILE", unset = "/home/k-makino/TCGA/Pancancer/GDCdata/clinical_TCGA_CDR_all_cancers.tsv")
REF_TME <- Sys.getenv("REF_TME", unset = "TME5")
PROJECTION_SUBDIR <- Sys.getenv("PROJECTION_SUBDIR", unset = "module_signature_template_projection_noTME2")
TME_LEVELS <- strsplit(Sys.getenv("TME_LEVELS", unset = "TME1,TME3,TME4,TME5,TME6"), ",")[[1]]
TME_LEVELS <- trimws(TME_LEVELS)
TME_LEVELS <- TME_LEVELS[TME_LEVELS != ""]
if (!REF_TME %in% TME_LEVELS) {
  stop("REF_TME must be one of TME_LEVELS. REF_TME=", REF_TME, "; TME_LEVELS=", paste(TME_LEVELS, collapse = ","))
}
PROJECTS <- strsplit(Sys.getenv("PROJECTS", unset = "BLCA,BRCA,CESC,COAD,ESCA,GBM,HNSC,KIRC,KIRP,LGG,LIHC,LUAD,LUSC,OV,PRAD,SARC,SKCM,STAD,THCA,UCEC"), ",")[[1]]
PROJECTS <- trimws(PROJECTS)
MIN_EVENTS <- as.integer(Sys.getenv("MIN_EVENTS", unset = "5"))
MIN_N_PER_MODEL <- as.integer(Sys.getenv("MIN_N_PER_MODEL", unset = "20"))

OUT_DIR <- file.path(BASE_DIR, paste0("TCGA_multivariable_Cox_TME_noTME2_gender_cancertype_", RUN_TAG))
dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)

message("BASE_DIR: ", BASE_DIR)
message("RUN_TAG:  ", RUN_TAG)
message("CDR_FILE: ", CDR_FILE)
message("REF_TME:  ", REF_TME)
message("TME_LEVELS: ", paste(TME_LEVELS, collapse = ", "))
message("PROJECTION_SUBDIR: ", PROJECTION_SUBDIR)
message("OUT_DIR:  ", OUT_DIR)
message("PROJECTS: ", paste(PROJECTS, collapse = ", "))
message("==============================")

## ----------------------------
## Helpers
## ----------------------------
first_existing <- function(x, candidates) {
  hit <- candidates[candidates %in% x]
  if (length(hit) == 0) return(NA_character_)
  hit[[1]]
}

first_existing_ci <- function(x, candidates) {
  lx <- tolower(x)
  lc <- tolower(candidates)
  idx <- match(lc, lx)
  idx <- idx[!is.na(idx)]
  if (length(idx) == 0) return(NA_character_)
  x[idx[[1]]]
}

clean_barcode <- function(x) {
  x <- as.character(x)
  x <- trimws(x)
  x <- gsub("\\.", "-", x)
  x
}

tcga_patient_id <- function(x) {
  x <- clean_barcode(x)
  out <- ifelse(grepl("^TCGA-[A-Z0-9]{2}-[A-Z0-9]{4}", x), substr(x, 1, 12), x)
  out
}

tcga_sample15 <- function(x) {
  x <- clean_barcode(x)
  out <- ifelse(grepl("^TCGA-[A-Z0-9]{2}-[A-Z0-9]{4}-[0-9A-Z]{2}", x), substr(x, 1, 15), x)
  out
}

norm_project <- function(x) {
  x <- toupper(as.character(x))
  x <- gsub("^TCGA[-_]", "", x)
  x <- gsub("[^A-Z0-9]", "", x)
  x
}

clean_gender <- function(x) {
  x <- trimws(tolower(as.character(x)))
  x[x %in% c("male", "m", "1")] <- "Male"
  x[x %in% c("female", "f", "0")] <- "Female"
  x[!(x %in% c("Male", "Female"))] <- NA_character_
  factor(x, levels = c("Female", "Male"))
}

clean_tme <- function(x) {
  x <- as.character(x)
  x <- trimws(x)
  x <- gsub("_", "", x)
  x <- toupper(x)
  x <- ifelse(grepl("TME[1-6]", x), sub(".*(TME[1-6]).*", "\\1", x), x)
  x
}

pick_assignment_columns <- function(df) {
  cn <- colnames(df)
  sample_col <- first_existing_ci(cn, c("sample", "Sample", "ID", "sample_id", "patient_id", "barcode", "ID_full"))
  if (is.na(sample_col)) sample_col <- cn[[1]]
  tme_col <- first_existing_ci(cn, c("predicted_TME", "TME", "TME_assignment", "assigned_TME", "best_TME", "cluster", "group"))
  if (is.na(tme_col)) {
    cand <- cn[grepl("TME", cn, ignore.case = TRUE)]
    if (length(cand) > 0) tme_col <- cand[[1]]
  }
  if (is.na(tme_col)) stop("Could not detect TME assignment column in integrated_TME_assignment.tsv")
  margin_col <- first_existing_ci(cn, c("margin", "score_margin", "delta", "best_margin"))
  list(sample_col = sample_col, tme_col = tme_col, margin_col = margin_col)
}

read_assignment <- function(project) {
  f <- file.path(BASE_DIR, project, RUN_TAG, PROJECTION_SUBDIR, "integrated_TME_assignment.tsv")
  if (!file.exists(f)) stop("assignment file not found: ", f)
  x <- fread(f, data.table = FALSE, check.names = FALSE)
  cols <- pick_assignment_columns(x)
  out <- x %>%
    transmute(
      sample = clean_barcode(.data[[cols$sample_col]]),
      patient_id = tcga_patient_id(.data[[cols$sample_col]]),
      sample15 = tcga_sample15(.data[[cols$sample_col]]),
      TME = clean_tme(.data[[cols$tme_col]]),
      margin = if (!is.na(cols$margin_col)) suppressWarnings(as.numeric(.data[[cols$margin_col]])) else NA_real_
    ) %>%
    filter(!is.na(sample), !is.na(patient_id), !is.na(TME), TME %in% TME_LEVELS)
  ## collapse multiple aliquots/samples per patient; keep highest confidence if margin exists
  out <- out %>%
    mutate(margin2 = ifelse(is.na(margin), -Inf, margin)) %>%
    arrange(patient_id, desc(margin2), sample) %>%
    group_by(patient_id) %>%
    slice(1) %>%
    ungroup() %>%
    select(-margin2)
  out
}

read_tcga_cdr <- function() {
  if (!file.exists(CDR_FILE)) stop("TCGA CDR file not found: ", CDR_FILE)
  cdr <- fread(CDR_FILE, data.table = FALSE, check.names = FALSE)
  cn <- colnames(cdr)
  os_col <- first_existing_ci(cn, c("OS", "OS_event", "os"))
  time_col <- first_existing_ci(cn, c("OS.time", "OS_time", "os.time", "os_time"))
  gender_col <- first_existing_ci(cn, c("gender", "sex", "Gender", "Sex"))
  project_col <- first_existing_ci(cn, c("cancer_type", "type", "project", "project_id", "Project"))
  patient_col <- first_existing_ci(cn, c("patient_id", "bcr_patient_barcode", "Patient_Identifier", "case_submitter_id", "submitter_id", "barcode"))

  ## If no explicit patient column is detected, use the first column when it looks like a TCGA barcode.
  if (is.na(patient_col)) {
    tcga_like_cols <- cn[vapply(cdr, function(v) mean(grepl("^TCGA-", as.character(v))) > 0.1, numeric(1))]
    if (length(tcga_like_cols) > 0) patient_col <- tcga_like_cols[[1]]
  }
  if (is.na(patient_col)) stop("Could not detect patient barcode column in CDR file.")
  if (is.na(os_col)) stop("Could not detect OS event column in CDR file.")
  if (is.na(time_col)) stop("Could not detect OS time column in CDR file.")
  if (is.na(gender_col)) {
    warning("Could not detect gender/sex column in CDR file. Gender adjustment will be skipped where unavailable.")
    cdr$.gender_tmp <- NA_character_
    gender_col <- ".gender_tmp"
  }
  if (is.na(project_col)) {
    cdr$.project_tmp <- NA_character_
    project_col <- ".project_tmp"
  }

  out <- cdr %>%
    mutate(
      patient_id = tcga_patient_id(.data[[patient_col]]),
      cancer_type = norm_project(.data[[project_col]]),
      OS_event = suppressWarnings(as.numeric(.data[[os_col]])),
      OS_time = suppressWarnings(as.numeric(.data[[time_col]])),
      gender = clean_gender(.data[[gender_col]])
    ) %>%
    filter(!is.na(patient_id), patient_id != "", !is.na(OS_time), OS_time > 0, !is.na(OS_event)) %>%
    mutate(OS_event = ifelse(OS_event > 0, 1, 0)) %>%
    distinct(patient_id, .keep_all = TRUE) %>%
    select(patient_id, cancer_type, OS_time, OS_event, gender, everything())
  attr(out, "columns") <- list(patient_col = patient_col, project_col = project_col, os_col = os_col, time_col = time_col, gender_col = gender_col)
  out
}

add_reference_rows <- function(coef_df, present_tmes, model_label, cohort, endpoint, adjustment_note) {
  all_tmes <- TME_LEVELS
  rows <- lapply(all_tmes, function(tm) {
    if (tm == REF_TME) {
      data.frame(
        model_label = model_label,
        cohort = cohort,
        endpoint = endpoint,
        TME = tm,
        reference = TRUE,
        HR = 1,
        HR_lower95 = 1,
        HR_upper95 = 1,
        coef = 0,
        se = NA_real_,
        z = NA_real_,
        p = NA_real_,
        n = NA_integer_,
        events = NA_integer_,
        adjustment_note = adjustment_note,
        note = "reference",
        stringsAsFactors = FALSE
      )
    } else if (tm %in% coef_df$TME) {
      coef_df[coef_df$TME == tm, , drop = FALSE]
    } else {
      data.frame(
        model_label = model_label,
        cohort = cohort,
        endpoint = endpoint,
        TME = tm,
        reference = FALSE,
        HR = NA_real_,
        HR_lower95 = NA_real_,
        HR_upper95 = NA_real_,
        coef = NA_real_,
        se = NA_real_,
        z = NA_real_,
        p = NA_real_,
        n = NA_integer_,
        events = NA_integer_,
        adjustment_note = adjustment_note,
        note = ifelse(tm %in% present_tmes, "not_estimable", "absent"),
        stringsAsFactors = FALSE
      )
    }
  })
  bind_rows(rows)
}

extract_tme_terms <- function(fit, model_label, cohort, endpoint, dat, adjustment_note) {
  sm <- summary(fit)
  co <- as.data.frame(sm$coefficients)
  ci <- as.data.frame(sm$conf.int)
  if (nrow(co) == 0) {
    coef_df <- data.frame()
  } else {
    co$term <- rownames(co)
    ci$term <- rownames(ci)
    p_col <- grep("Pr\\(>\\|z\\|\\)", colnames(co), value = TRUE)
    if (length(p_col) == 0) p_col <- grep("Pr", colnames(co), value = TRUE)[1]
    tme_terms <- co$term[grepl("^TME_factor", co$term)]
    coef_df <- lapply(tme_terms, function(term) {
      tm <- sub("^TME_factor", "", term)
      data.frame(
        model_label = model_label,
        cohort = cohort,
        endpoint = endpoint,
        TME = tm,
        reference = FALSE,
        HR = ci[term, "exp(coef)"],
        HR_lower95 = ci[term, "lower .95"],
        HR_upper95 = ci[term, "upper .95"],
        coef = co[term, "coef"],
        se = co[term, "se(coef)"],
        z = co[term, "z"],
        p = co[term, p_col],
        n = sum(dat$TME == tm, na.rm = TRUE),
        events = sum(dat$OS_event[dat$TME == tm] == 1, na.rm = TRUE),
        adjustment_note = adjustment_note,
        note = "OK",
        stringsAsFactors = FALSE
      )
    }) %>% bind_rows()
  }
  add_reference_rows(coef_df, present_tmes = sort(unique(dat$TME)), model_label, cohort, endpoint, adjustment_note)
}

fit_cox_model <- function(dat, model_label, cohort, endpoint, adjust_cancer = FALSE, adjust_gender = TRUE) {
  dat <- dat %>%
    filter(!is.na(OS_time), OS_time > 0, !is.na(OS_event), !is.na(TME)) %>%
    mutate(
      TME = clean_tme(TME),
      TME_factor = factor(TME, levels = c(REF_TME, setdiff(TME_LEVELS, REF_TME))),
      cancer_type_factor = factor(cancer_type),
      gender = droplevels(clean_gender(gender))
    ) %>%
    filter(TME %in% TME_LEVELS)

  if (!REF_TME %in% dat$TME) stop("Reference ", REF_TME, " is absent in ", cohort)
  if (nrow(dat) < MIN_N_PER_MODEL) stop("Too few matched rows: ", nrow(dat))
  if (sum(dat$OS_event == 1, na.rm = TRUE) < MIN_EVENTS) stop("Too few events: ", sum(dat$OS_event == 1, na.rm = TRUE))
  if (n_distinct(dat$TME) < 2) stop("Only one TME group present.")

  rhs <- c("TME_factor")
  notes <- c()
  if (adjust_cancer) {
    if (n_distinct(dat$cancer_type_factor) >= 2) {
      rhs <- c(rhs, "cancer_type_factor")
      notes <- c(notes, "cancer_type_adjusted")
    } else {
      notes <- c(notes, "cancer_type_not_adjusted_one_level")
    }
  }
  if (adjust_gender) {
    if (sum(!is.na(dat$gender)) >= MIN_N_PER_MODEL && n_distinct(dat$gender[!is.na(dat$gender)]) >= 2) {
      rhs <- c(rhs, "gender")
      notes <- c(notes, "gender_adjusted")
      dat_model <- dat %>% filter(!is.na(gender))
    } else {
      notes <- c(notes, "gender_not_adjusted_one_or_missing_level")
      dat_model <- dat
    }
  } else {
    dat_model <- dat
    notes <- c(notes, "gender_not_requested")
  }
  if (nrow(dat_model) < MIN_N_PER_MODEL) stop("Too few rows after covariate filtering: ", nrow(dat_model))
  if (sum(dat_model$OS_event == 1, na.rm = TRUE) < MIN_EVENTS) stop("Too few events after covariate filtering: ", sum(dat_model$OS_event == 1, na.rm = TRUE))

  form <- as.formula(paste("Surv(OS_time, OS_event) ~", paste(rhs, collapse = " + ")))
  fit <- coxph(form, data = dat_model, ties = "efron")
  res <- extract_tme_terms(fit, model_label, cohort, endpoint, dat_model, paste(notes, collapse = ";"))
  list(result = res, fit = fit, data = dat_model, formula = deparse(form))
}

plot_forest <- function(df, out_pdf, title, facet_by_cohort = FALSE) {
  pdat <- df %>%
    mutate(
      TME = factor(TME, levels = rev(TME_LEVELS)),
      HR_plot = HR,
      HR_lower_plot = HR_lower95,
      HR_upper_plot = HR_upper95,
      significant = ifelse(!is.na(p_BH_global) & p_BH_global < 0.05, "BH < 0.05", "n.s./NA")
    ) %>%
    filter(!is.na(HR_plot), HR_plot > 0)
  if (nrow(pdat) == 0) return(invisible(NULL))

  g <- ggplot(pdat, aes(x = HR_plot, y = TME)) +
    geom_vline(xintercept = 1, linetype = "dashed", linewidth = 0.4) +
    geom_errorbarh(aes(xmin = HR_lower_plot, xmax = HR_upper_plot), height = 0.18, linewidth = 0.35, na.rm = TRUE) +
    geom_point(aes(shape = reference, alpha = significant), size = 2.3, na.rm = TRUE) +
    scale_x_log10() +
    scale_shape_manual(values = c(`FALSE` = 16, `TRUE` = 18), guide = "none") +
    scale_alpha_manual(values = c("BH < 0.05" = 1, "n.s./NA" = 0.55), guide = "none") +
    labs(title = title, x = paste0("Hazard ratio relative to ", REF_TME, " (log scale)"), y = NULL) +
    theme_bw(base_size = 10) +
    theme(
      plot.title = element_text(face = "bold", hjust = 0.5),
      panel.grid.minor = element_blank()
    )
  if (facet_by_cohort) {
    g <- g + facet_wrap(~ cohort, ncol = 3)
  }
  ggsave(out_pdf, g, width = ifelse(facet_by_cohort, 13, 6.5), height = ifelse(facet_by_cohort, 11, 4.8), useDingbats = FALSE)
  invisible(g)
}

## ----------------------------
## Load data
## ----------------------------
cdr <- read_tcga_cdr()
message("CDR columns used:")
print(attr(cdr, "columns"))

assignment_list <- list()
input_list <- list()
diagnostics <- list()

for (project in PROJECTS) {
  message("==============================")
  message("Loading assignment: ", project)
  a <- tryCatch(read_assignment(project), error = function(e) {
    message("ERROR reading assignment for ", project, ": ", conditionMessage(e))
    NULL
  })
  if (is.null(a)) next

  cdr_project <- cdr %>% filter(cancer_type == norm_project(project))
  if (nrow(cdr_project) == 0) {
    ## Fallback: use all CDR and rely on barcode matching.
    cdr_project <- cdr
  }
  dat <- a %>%
    inner_join(cdr_project %>% select(patient_id, clinical_cancer_type = cancer_type, OS_time, OS_event, gender), by = "patient_id") %>%
    mutate(cancer_type = project, endpoint = "OS")

  assignment_list[[project]] <- a
  input_list[[project]] <- dat
  diagnostics[[project]] <- data.frame(
    cohort = project,
    n_assignment_patients = nrow(a),
    n_cdr_project_rows = nrow(cdr_project),
    n_matched_rows = nrow(dat),
    n_events = sum(dat$OS_event == 1, na.rm = TRUE),
    n_TME1 = sum(dat$TME == "TME1", na.rm = TRUE),
    n_TME3 = sum(dat$TME == "TME3", na.rm = TRUE),
    n_TME4 = sum(dat$TME == "TME4", na.rm = TRUE),
    n_TME5 = sum(dat$TME == "TME5", na.rm = TRUE),
    n_TME6 = sum(dat$TME == "TME6", na.rm = TRUE),
    n_gender_levels = n_distinct(dat$gender[!is.na(dat$gender)]),
    example_assignment_id = paste(head(a$patient_id, 3), collapse = ";"),
    example_matched_id = paste(head(dat$patient_id, 3), collapse = ";"),
    stringsAsFactors = FALSE
  )
}

all_inputs <- bind_rows(input_list)
all_diagnostics <- bind_rows(diagnostics)
fwrite(all_diagnostics, file.path(OUT_DIR, "TCGA_TME_multivariable_cox_matching_diagnostics.tsv"), sep = "\t")
fwrite(all_inputs, file.path(OUT_DIR, "TCGA_TME_multivariable_cox_matched_input.tsv"), sep = "\t")

if (nrow(all_inputs) == 0) stop("No matched TCGA samples across all projects.")

## ----------------------------
## Per-cancer models
## ----------------------------
per_cancer_results <- list()
per_cancer_formula <- list()

for (project in PROJECTS) {
  message("==============================")
  message("Fitting per-cancer Cox model: ", project)
  dat <- input_list[[project]]
  if (is.null(dat) || nrow(dat) == 0) {
    per_cancer_results[[project]] <- data.frame(
      model_label = "within_cancer_gender_adjusted",
      cohort = project,
      endpoint = "OS",
      TME = TME_LEVELS,
      reference = TME_LEVELS == REF_TME,
      HR = NA_real_, HR_lower95 = NA_real_, HR_upper95 = NA_real_,
      coef = NA_real_, se = NA_real_, z = NA_real_, p = NA_real_,
      n = NA_integer_, events = NA_integer_,
      adjustment_note = NA_character_,
      note = "no_matched_samples",
      stringsAsFactors = FALSE
    )
    next
  }
  ans <- tryCatch(
    fit_cox_model(dat, "within_cancer_gender_adjusted", project, "OS", adjust_cancer = FALSE, adjust_gender = TRUE),
    error = function(e) e
  )
  if (inherits(ans, "error")) {
    message("ERROR in ", project, ": ", conditionMessage(ans))
    per_cancer_results[[project]] <- data.frame(
      model_label = "within_cancer_gender_adjusted",
      cohort = project,
      endpoint = "OS",
      TME = TME_LEVELS,
      reference = TME_LEVELS == REF_TME,
      HR = NA_real_, HR_lower95 = NA_real_, HR_upper95 = NA_real_,
      coef = NA_real_, se = NA_real_, z = NA_real_, p = NA_real_,
      n = NA_integer_, events = NA_integer_,
      adjustment_note = NA_character_,
      note = paste0("ERROR: ", conditionMessage(ans)),
      stringsAsFactors = FALSE
    )
  } else {
    per_cancer_results[[project]] <- ans$result
    per_cancer_formula[[project]] <- ans$formula
  }
}

per_cancer_tbl <- bind_rows(per_cancer_results) %>%
  group_by(cohort) %>%
  mutate(p_BH_within_cohort = p.adjust(p, method = "BH"),
         p_Holm_within_cohort = p.adjust(p, method = "holm")) %>%
  ungroup() %>%
  mutate(
    p_BH_global = p.adjust(p, method = "BH"),
    p_Holm_global = p.adjust(p, method = "holm"),
    direction_poor_vs_ref = ifelse(!is.na(HR) & HR > 1 & TME != REF_TME, TRUE, FALSE)
  )

fwrite(per_cancer_tbl, file.path(OUT_DIR, "TCGA_within_cancer_TME_HR_gender_adjusted.tsv"), sep = "\t")

## ----------------------------
## Pan-cancer integrated model
## ----------------------------
message("==============================")
message("Fitting pan-cancer Cox model: TME + cancer_type + gender")
pan_ans <- tryCatch(
  fit_cox_model(all_inputs, "pan_cancer_cancertype_gender_adjusted", "TCGA_pan_cancer", "OS", adjust_cancer = TRUE, adjust_gender = TRUE),
  error = function(e) e
)

if (inherits(pan_ans, "error")) {
  message("ERROR in pan-cancer model: ", conditionMessage(pan_ans))
  pan_tbl <- data.frame(
    model_label = "pan_cancer_cancertype_gender_adjusted",
    cohort = "TCGA_pan_cancer",
    endpoint = "OS",
    TME = TME_LEVELS,
    reference = TME_LEVELS == REF_TME,
    HR = NA_real_, HR_lower95 = NA_real_, HR_upper95 = NA_real_,
    coef = NA_real_, se = NA_real_, z = NA_real_, p = NA_real_,
    n = NA_integer_, events = NA_integer_,
    adjustment_note = NA_character_,
    note = paste0("ERROR: ", conditionMessage(pan_ans)),
    stringsAsFactors = FALSE
  )
} else {
  pan_tbl <- pan_ans$result
  writeLines(pan_ans$formula, con = file.path(OUT_DIR, "TCGA_pan_cancer_Cox_formula.txt"))
}

pan_tbl <- pan_tbl %>%
  mutate(
    p_BH_global = p.adjust(p, method = "BH"),
    p_Holm_global = p.adjust(p, method = "holm"),
    direction_poor_vs_ref = ifelse(!is.na(HR) & HR > 1 & TME != REF_TME, TRUE, FALSE)
  )
fwrite(pan_tbl, file.path(OUT_DIR, "TCGA_pan_cancer_TME_HR_cancertype_gender_adjusted.tsv"), sep = "\t")

## ----------------------------
## Counts and plots
## ----------------------------
counts_tbl <- all_inputs %>%
  count(cancer_type, TME, gender, name = "n") %>%
  arrange(cancer_type, TME, gender)
fwrite(counts_tbl, file.path(OUT_DIR, "TCGA_TME_gender_counts.tsv"), sep = "\t")

plot_forest(
  pan_tbl,
  file.path(OUT_DIR, "TCGA_pan_cancer_TME_HR_cancertype_gender_adjusted_forest.pdf"),
  title = paste0("TCGA pan-cancer Cox model: TME + cancer type + gender (ref = ", REF_TME, ")"),
  facet_by_cohort = FALSE
)

plot_forest(
  per_cancer_tbl,
  file.path(OUT_DIR, "TCGA_within_cancer_TME_HR_gender_adjusted_forest.pdf"),
  title = paste0("TCGA within-cancer Cox models: TME + gender (ref = ", REF_TME, ")"),
  facet_by_cohort = TRUE
)

## Combined export for convenience
combined_tbl <- bind_rows(per_cancer_tbl, pan_tbl)
fwrite(combined_tbl, file.path(OUT_DIR, "TCGA_TME_multivariable_Cox_all_results.tsv"), sep = "\t")

message("==============================")
message("Done.")
message("Output directory: ", OUT_DIR)
message("Main outputs:")
message("  ", file.path(OUT_DIR, "TCGA_within_cancer_TME_HR_gender_adjusted.tsv"))
message("  ", file.path(OUT_DIR, "TCGA_pan_cancer_TME_HR_cancertype_gender_adjusted.tsv"))
message("  ", file.path(OUT_DIR, "TCGA_within_cancer_TME_HR_gender_adjusted_forest.pdf"))
message("  ", file.path(OUT_DIR, "TCGA_pan_cancer_TME_HR_cancertype_gender_adjusted_forest.pdf"))
