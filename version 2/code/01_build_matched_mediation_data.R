#!/usr/bin/env Rscript

options(stringsAsFactors = FALSE)

if (!requireNamespace("R.matlab", quietly = TRUE)) {
  stop("Package 'R.matlab' is required. Install it with install.packages('R.matlab').")
}

project_root <- normalizePath(getwd(), mustWork = TRUE)
shape_dir <- file.path(project_root, "data", "ShapeMAdata")
spls_dir <- file.path(project_root, "data", "SPLS-Coxdata")
out_dir <- file.path(project_root, "data")

required_files <- c(
  file.path(shape_dir, "ADNI_CCseg_707subjects_SRVF.mat"),
  file.path(shape_dir, "ADNI_CCseg_707subjects_info_v2.csv"),
  file.path(shape_dir, "snps_correct_negative_pc5.mat"),
  file.path(shape_dir, "top_snp_4_correct_negative.mat"),
  file.path(shape_dir, "x_design_correct_negative.mat"),
  file.path(spls_dir, "clinical.dat"),
  file.path(spls_dir, "ccdata.csv")
)
missing_files <- required_files[!file.exists(required_files)]
if (length(missing_files) > 0) {
  stop("Missing required files:\n", paste(missing_files, collapse = "\n"))
}

info <- read.csv(
  file.path(shape_dir, "ADNI_CCseg_707subjects_info_v2.csv"),
  check.names = FALSE
)
if (!"RID" %in% names(info)) stop("The ShapeMA information table has no RID column.")
info$RID <- as.integer(info$RID)
if (anyNA(info$RID) || anyDuplicated(info$RID)) {
  stop("The ShapeMA information table has missing or duplicated RID values.")
}
if (!"PTEDUCAT" %in% names(info)) {
  stop("The ShapeMA information table has no PTEDUCAT column.")
}
education_years_all <- as.numeric(info$PTEDUCAT)
if (anyNA(education_years_all) || any(!is.finite(education_years_all))) {
  stop("PTEDUCAT must be complete and numeric before it is added to W.")
}

read_shape_mat <- suppressWarnings(R.matlab::readMat(
  file.path(shape_dir, "ADNI_CCseg_707subjects_SRVF.mat")
))
CC_q <- read_shape_mat$CC.q
if (is.null(CC_q) || !identical(dim(CC_q), c(2L, 707L, 100L))) {
  stop("Unexpected CC_q dimensions; expected 2 x 707 x 100.")
}

read_snp_mat <- suppressWarnings(R.matlab::readMat(
  file.path(shape_dir, "top_snp_4_correct_negative.mat")
))
X_all <- read_snp_mat$top.snp.4
snp_names <- as.character(read_snp_mat$top.snp.4.name[, 1])
if (is.null(X_all) || !identical(dim(X_all), c(707L, 4L))) {
  stop("Unexpected top_snp_4 dimensions; expected 707 x 4.")
}
colnames(X_all) <- snp_names

read_pc5_mat <- suppressWarnings(R.matlab::readMat(
  file.path(shape_dir, "snps_correct_negative_pc5.mat")
))
W_pc5_all <- read_pc5_mat$x.design
if (is.null(W_pc5_all) || !identical(dim(W_pc5_all), c(707L, 10L))) {
  stop("Unexpected 5-PC design dimensions; expected 707 x 10.")
}
colnames(W_pc5_all) <- c(
  "intercept", "male", "age_z", "handedness", "APOE4",
  paste0("PC", 1:5)
)
W_pc5_all <- cbind(W_pc5_all, education_years = education_years_all)

read_pc2_mat <- suppressWarnings(R.matlab::readMat(
  file.path(shape_dir, "x_design_correct_negative.mat")
))
W_pc2_all <- read_pc2_mat$x.design
if (is.null(W_pc2_all) || !identical(dim(W_pc2_all), c(707L, 7L))) {
  stop("Unexpected 2-PC design dimensions; expected 707 x 7.")
}
colnames(W_pc2_all) <- c("intercept", "male", "age_z", "handedness", "APOE4", "PC1", "PC2")
W_pc2_all <- cbind(W_pc2_all, education_years = education_years_all)

spls_env <- new.env(parent = emptyenv())
load(file.path(spls_dir, "clinical.dat"), envir = spls_env)
if (!exists("scoredata", envir = spls_env, inherits = FALSE)) {
  stop("clinical.dat does not contain the expected object 'scoredata'.")
}
scoredata <- get("scoredata", envir = spls_env, inherits = FALSE)
if (!is.matrix(scoredata) || ncol(scoredata) < 3L) {
  stop("Unexpected scoredata structure in clinical.dat.")
}

# The SPLS-Cox repository removes this row because the corresponding CC shape
# row is missing; this reproduces the repository's real-data analysis cohort.
drop_row <- 331L
if (nrow(scoredata) != 373L) {
  stop("Unexpected scoredata row count; expected 373 before the repository's row-331 removal.")
}
scoredata_clean <- scoredata[-drop_row, , drop = FALSE]
if (nrow(scoredata_clean) != nrow(read.csv(file.path(spls_dir, "ccdata.csv"), header = FALSE))) {
  stop("SPLS-Cox clinical and shape row counts do not agree after row removal.")
}

survival_data <- data.frame(
  spls_row = seq_len(nrow(scoredata_clean)),
  RID = as.integer(scoredata_clean[, 3]),
  time = as.numeric(scoredata_clean[, 1]),
  event = as.integer(scoredata_clean[, 2])
)
if (anyNA(survival_data$RID) || anyDuplicated(survival_data$RID)) {
  stop("The SPLS-Cox survival data has missing or duplicated RID values.")
}
if (!all(survival_data$event %in% c(0L, 1L)) || any(!is.finite(survival_data$time)) || any(survival_data$time <= 0)) {
  stop("Invalid time/event values in SPLS-Cox clinical.dat.")
}

current_row <- match(survival_data$RID, info$RID)
matched <- !is.na(current_row)
if (sum(matched) == 0L) stop("No RID values could be matched.")

survival_matched <- survival_data[matched, , drop = FALSE]
current_row_matched <- current_row[matched]
subject <- info[current_row_matched, , drop = FALSE]
subject$source_shape_row <- current_row_matched
subject$source_spls_row <- survival_matched$spls_row
subject$time <- survival_matched$time
subject$event <- survival_matched$event

M <- aperm(CC_q[, current_row_matched, , drop = FALSE], c(2, 3, 1))
X <- X_all[current_row_matched, , drop = FALSE]
W_pc5 <- W_pc5_all[current_row_matched, , drop = FALSE]
W_pc2 <- W_pc2_all[current_row_matched, , drop = FALSE]

if (!identical(as.integer(subject$RID), as.integer(survival_matched$RID))) {
  stop("Internal RID order check failed after matching.")
}
if (!identical(dim(M), c(nrow(subject), 100L, 2L))) {
  stop("Unexpected matched mediator dimensions.")
}

model_data <- list(
  subject = subject,
  X = X,
  M = M,
  W_pc5 = W_pc5,
  W_pc2 = W_pc2,
  snp_names = snp_names,
  shape_grid = seq(0, 1, length.out = 100),
  source = list(
    shapema_dir = shape_dir,
    spls_cox_dir = spls_dir,
    shape_file = "ADNI_CCseg_707subjects_SRVF.mat",
    survival_file = "clinical.dat",
    matching_key = "RID",
    spls_removed_row = drop_row,
    education_variable = "PTEDUCAT (years)"
  )
)

saveRDS(model_data, file.path(out_dir, "mediation_model_matched.rds"))
save(model_data, file = file.path(out_dir, "mediation_model_matched.RData"))
write.csv(subject, file.path(out_dir, "mediation_model_matched_subjects.csv"), row.names = FALSE)

audit <- merge(
  survival_data,
  data.frame(RID = info$RID, shapeMA_row = seq_len(nrow(info))),
  by = "RID", all.x = TRUE, sort = FALSE
)
audit$matched <- !is.na(audit$shapeMA_row)
write.csv(audit, file.path(out_dir, "mediation_model_RID_match_audit.csv"), row.names = FALSE)

cat("Matched subjects:", nrow(subject), "\n")
cat("Events:", sum(subject$event == 1L), " Censored:", sum(subject$event == 0L), "\n")
cat("Mediator dimensions (subject x grid x coordinate):", paste(dim(M), collapse = " x "), "\n")
cat("Outputs:\n")
cat(file.path(out_dir, "mediation_model_matched.rds"), "\n")
cat(file.path(out_dir, "mediation_model_matched.RData"), "\n")
cat(file.path(out_dir, "mediation_model_matched_subjects.csv"), "\n")
cat(file.path(out_dir, "mediation_model_RID_match_audit.csv"), "\n")
