#!/usr/bin/env Rscript

script_dir <- function() {
  args <- commandArgs(trailingOnly = FALSE)
  file_arg <- grep("^--file=", args, value = TRUE)
  if (length(file_arg) == 1L) {
    return(dirname(normalizePath(sub("^--file=", "", file_arg))))
  }
  normalizePath(getwd())
}

code_dir <- script_dir()
project_dir <- dirname(code_dir)
matched_file <- file.path(project_dir, "data", "mediation_model_matched.rds")
spls_file <- file.path(project_dir, "data", "SPLS-Coxdata", "clinical.dat")
output_dir <- file.path(project_dir, "data")

if (!file.exists(matched_file)) stop("Matched data not found: ", matched_file)
if (!file.exists(spls_file)) stop("SPLS clinical data not found: ", spls_file)

matched <- readRDS(matched_file)
spls_env <- new.env(parent = emptyenv())
load(spls_file, envir = spls_env)
if (!exists("scoredata", envir = spls_env, inherits = FALSE)) {
  stop("clinical.dat does not contain scoredata.")
}
scoredata <- get("scoredata", envir = spls_env, inherits = FALSE)
if (!is.matrix(scoredata) || ncol(scoredata) != 15L || nrow(scoredata) != 373L) {
  stop("Unexpected scoredata dimensions; expected 373 x 15.")
}

scoredata_clean <- scoredata[-331L, , drop = FALSE]
subject <- as.data.frame(matched$subject)
if (!all(c("source_spls_row", "RID", "time", "event") %in% names(subject))) {
  stop("Matched subject data lacks source_spls_row, RID, time, or event.")
}
if (anyNA(subject$source_spls_row) || anyDuplicated(subject$source_spls_row)) {
  stop("source_spls_row must be complete and unique.")
}
if (any(subject$source_spls_row < 1L | subject$source_spls_row > nrow(scoredata_clean))) {
  stop("source_spls_row is outside the cleaned SPLS clinical data.")
}

spls_matched <- scoredata_clean[subject$source_spls_row, , drop = FALSE]
if (!identical(as.integer(spls_matched[, 3]), as.integer(subject$RID))) {
  stop("SPLS clinical rows do not agree with subject RIDs.")
}
if (!identical(as.numeric(spls_matched[, 1]), as.numeric(subject$time)) ||
    !identical(as.integer(spls_matched[, 2]), as.integer(subject$event))) {
  stop("SPLS time/event values do not agree with matched subject data.")
}

spls_names <- c("SPLS_time", "SPLS_event", "SPLS_RID", paste0("SPLS_clinical_col", sprintf("%02d", 4:15)))
spls_clinical_raw <- as.data.frame(spls_matched, check.names = FALSE)
names(spls_clinical_raw) <- spls_names

analysis_frame <- subject
for (snp in colnames(matched$X)) {
  analysis_frame[[paste0("X_", snp)]] <- as.numeric(matched$X[, snp])
}
for (variable in colnames(matched$W_pc2)) {
  analysis_frame[[paste0("W_pc2_", variable)]] <- as.numeric(matched$W_pc2[, variable])
}
for (variable in colnames(matched$W_pc5)) {
  analysis_frame[[paste0("W_pc5_", variable)]] <- as.numeric(matched$W_pc5[, variable])
}
analysis_frame <- cbind(analysis_frame, spls_clinical_raw)

candidate_covariates <- c(
  "AGE", "PTGENDER", "PTEDUCAT", "PTETHCAT", "PTRACCAT", "PTMARRY", "SITE",
  "APOE4", "PTHAND", "ICV", "ADAS11", "ADAS13", "MMSE", "CDRSB", "FAQ",
  "FDG", "ABETA", "PTAU", "TAU", "W_pc2_male", "W_pc2_age_z",
  "W_pc2_handedness", "W_pc2_APOE4", paste0("W_pc2_PC", 1:2),
  "W_pc2_education_years",
  "W_pc5_male", "W_pc5_age_z", "W_pc5_handedness", "W_pc5_APOE4",
  paste0("W_pc5_PC", 1:5), "W_pc5_education_years"
)
candidate_covariates <- candidate_covariates[candidate_covariates %in% names(analysis_frame)]

variable_catalog <- data.frame(
  variable = character(),
  role = character(),
  source = character(),
  notes = character(),
  stringsAsFactors = FALSE
)
add_catalog <- function(variable, role, source, notes) {
  data.frame(variable = variable, role = role, source = source, notes = notes, stringsAsFactors = FALSE)
}

variable_catalog <- rbind(
  variable_catalog,
  add_catalog("RID", "identifier", "ShapeMA info / SPLS clinical", "RID matching key"),
  add_catalog("time", "survival outcome", "SPLS clinical.dat", "Observed survival time"),
  add_catalog("event", "survival outcome", "SPLS clinical.dat", "1=AD conversion event; 0=censored"),
  add_catalog("DX.bl", "eligibility", "ShapeMA info", "Baseline diagnosis; LMCI is the intended MCI cohort"),
  add_catalog("M / CC_q", "functional mediator", "ShapeMA SRVF MAT", "333 x 100 x 2; one multivariate shape mediator"),
  add_catalog("X_*", "exposure", "ShapeMA SNP MAT", "Four SNPs; additive dosage 0/1/2"),
  add_catalog("W_pc2_*", "main adjustment set", "ShapeMA design MAT + ShapeMA info", "Includes PTEDUCAT as education_years"),
  add_catalog("W_pc5_*", "sensitivity adjustment set", "ShapeMA design MAT + ShapeMA info", "Includes PTEDUCAT as education_years"),
  add_catalog("PTEDUCAT", "adjustment variable", "ShapeMA info", "Education years; included in W_pc2 and W_pc5"),
  add_catalog("SITE", "candidate center effect", "ShapeMA info", "56 recruitment sites; consider stratification or center effect"),
  add_catalog("ICV", "candidate size adjustment", "ShapeMA info", "Intracranial volume"),
  add_catalog("ADAS13/MMSE/CDRSB/FAQ", "clinical candidate", "ShapeMA info", "Baseline cognition/function; causal role requires prespecification"),
  add_catalog("FDG/ABETA/PTAU/TAU", "biomarker candidate", "ShapeMA info", "Disease-related biomarkers; missingness and downstream status require review"),
  add_catalog("SPLS_clinical_col04:15", "raw clinical source fields", "SPLS clinical.dat", "Repository provides no column names; preserved by original position")
)

complete_data <- matched
complete_data$subject <- subject
complete_data$analysis_frame <- analysis_frame
complete_data$spls_clinical_raw <- spls_clinical_raw
complete_data$candidate_covariates <- analysis_frame[, candidate_covariates, drop = FALSE]
complete_data$candidate_covariate_names <- candidate_covariates
complete_data$variable_catalog <- variable_catalog
complete_data$source$complete_dataset_script <- "05_build_complete_analysis_dataset.R"
complete_data$source$spls_clinical_columns <- "SPLS_clinical_col04 through SPLS_clinical_col15 retain original unnamed positions"

saveRDS(complete_data, file.path(output_dir, "mediation_model_complete.rds"))
save(complete_data, file = file.path(output_dir, "mediation_model_complete.RData"))
write.csv(analysis_frame, file.path(output_dir, "mediation_model_complete_analysis_frame.csv"), row.names = FALSE)
write.csv(variable_catalog, file.path(output_dir, "mediation_model_complete_variable_catalog.csv"), row.names = FALSE)

cat("Complete analysis dataset created.\n")
cat("Subjects:", nrow(analysis_frame), "\n")
cat("Analysis-frame columns:", ncol(analysis_frame), "\n")
cat("Candidate covariates:", length(candidate_covariates), "\n")
cat("Outputs:\n")
cat(file.path(output_dir, "mediation_model_complete.rds"), "\n")
cat(file.path(output_dir, "mediation_model_complete.RData"), "\n")
cat(file.path(output_dir, "mediation_model_complete_analysis_frame.csv"), "\n")
cat(file.path(output_dir, "mediation_model_complete_variable_catalog.csv"), "\n")
