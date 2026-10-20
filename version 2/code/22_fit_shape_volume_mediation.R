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
penalized_script <- file.path(code_dir, "16_fit_penalized_functional_cox.R")
if (!file.exists(penalized_script)) {
  stop("Penalized functional Cox script not found: ", penalized_script)
}
# Also sources 02_fit_four_snp_mediation_models.R.
source(penalized_script, local = FALSE)

output_dir <- file.path(project_dir, "results", "shape_volume_mediation")
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

# Second mediator: total corpus callosum volume (sum of the five FreeSurfer
# subregions) divided by intracranial volume, standardized in the analysis
# sample. The standardization is a fixed linear transformation, so indirect
# effects do not depend on it.
volume_file <- file.path(project_dir, "data", "SPLS-Coxdata", "CCinfo.csv")
if (!file.exists(volume_file)) {
  stop("Corpus callosum volume file not found: ", volume_file)
}
cc_regions <- c(
  "CorpusCallosumAnterior", "CorpusCallosumMidAnterior", "CorpusCallosumCentral",
  "CorpusCallosumMidPosterior", "CorpusCallosumPosterior"
)
volume_table <- utils::read.csv(volume_file)
volume_rows <- match(as.integer(subject$RID), as.integer(volume_table$RID))
if (anyNA(volume_rows)) {
  stop("Corpus callosum volumes are missing for some analysis subjects.")
}
volume_table <- volume_table[volume_rows, , drop = FALSE]
if (anyNA(volume_table[, c(cc_regions, "ICV")])) {
  stop("Corpus callosum volumes or ICV contain missing values.")
}
volume_ratio <- rowSums(volume_table[, cc_regions]) / volume_table$ICV
volume_mediator <- cbind(
  cc_volume_icv_z = (volume_ratio - mean(volume_ratio)) / stats::sd(volume_ratio)
)

# Two mediators in one model:
#   M_j(s) = alpha_j(s) X + lambda_j(s)' W + e_j(s)        (shape, as before)
#   V      = a X + c' W + e                                 (volume, linear)
#   h(t)   = h0(t) exp[tau X + kappa' W + sum_j int M_j beta_j + b V]
# Indirect effects: shape = sum_j int alpha_j beta_j, volume = a b.
fit_shape_volume_replicate <- function(indices, replicate_id) {
  b_subject <- subject[indices, , drop = FALSE]
  b_X <- X[indices, , drop = FALSE]
  b_M <- M[indices, , , drop = FALSE]
  b_W <- W[indices, , drop = FALSE]
  b_volume <- volume_mediator[indices, , drop = FALSE]
  rows <- lapply(snp_names, function(snp) {
    x <- b_X[, snp]
    alpha <- fit_shape_layer(x, b_W, b_M, grid)$alpha
    volume_design <- cbind(SNP = as.numeric(x), b_W)
    snp_to_volume <- unname(qr.coef(qr(volume_design), b_volume)["SNP", 1L])
    fit_warnings <- character()
    cox_fit <- withCallingHandlers(
      fit_penalized_functional_cox(
        x, b_W, b_subject, b_M, grid, main_basis_size,
        scalar_mediators = b_volume
      ),
      warning = function(w) {
        fit_warnings <<- c(fit_warnings, conditionMessage(w))
        invokeRestart("muffleWarning")
      }
    )
    volume_to_hazard <- unname(cox_fit$mediator_coefficients[1L])
    shape <- penalized_effects(alpha, cox_fit$beta, cox_fit$grid_weights)
    indirect_volume <- snp_to_volume * volume_to_hazard
    indirect_total <- shape[["indirect_log_HR"]] + indirect_volume
    data.frame(
      replicate = replicate_id,
      SNP = snp,
      direct_log_HR = cox_fit$tau,
      indirect_shape_log_HR_coordinate_1 = shape[["indirect_log_HR_coordinate_1"]],
      indirect_shape_log_HR_coordinate_2 = shape[["indirect_log_HR_coordinate_2"]],
      indirect_shape_log_HR = shape[["indirect_log_HR"]],
      indirect_volume_log_HR = indirect_volume,
      indirect_total_log_HR = indirect_total,
      total_log_HR = cox_fit$tau + indirect_total,
      snp_to_volume = snp_to_volume,
      volume_to_log_hazard = volume_to_hazard,
      edf_beta_coordinate_1 = cox_fit$edf[1L],
      edf_beta_coordinate_2 = cox_fit$edf[2L],
      cox_warning = paste(unique(fit_warnings), collapse = " | ")
    )
  })
  do.call(rbind, rows)
}

run_shape_volume_point_analysis <- function() {
point_data <- fit_shape_volume_replicate(seq_len(n), 0L)

# Path-specific tests on the original data, for description only.
path_rows <- lapply(snp_names, function(snp) {
  x <- X[, snp]
  volume_fit <- stats::lm(volume_mediator[, 1L] ~ x + W[, -1L])
  cox_fit <- fit_penalized_functional_cox(
    x, W, subject, M, grid, main_basis_size, scalar_mediators = volume_mediator
  )
  parametric <- summary(cox_fit$fit)$p.table
  data.frame(
    SNP = snp,
    snp_to_volume_t = summary(volume_fit)$coefficients["x", "t value"],
    volume_to_log_hazard_z = parametric["scalarcc_volume_icv_z", "z value"],
    volume_to_log_hazard_p = parametric["scalarcc_volume_icv_z", "Pr(>|z|)"]
  )
})
point_data <- merge(point_data, do.call(rbind, path_rows), by = "SNP", sort = FALSE)
for (effect in c("direct", "indirect_shape", "indirect_volume", "indirect_total", "total")) {
  point_data[[paste0(effect, "_HR")]] <- exp(point_data[[paste0(effect, "_log_HR")]])
}

mediator_table <- data.frame(
  RID = subject$RID,
  cc_total_volume = rowSums(volume_table[, cc_regions]),
  ICV = volume_table$ICV,
  cc_volume_icv_ratio = volume_ratio,
  cc_volume_icv_z = volume_mediator[, 1L]
)

utils::write.csv(
  point_data,
  file.path(output_dir, "shape_volume_effect_summary.csv"),
  row.names = FALSE
)
utils::write.csv(
  mediator_table,
  file.path(output_dir, "cc_volume_mediator.csv"),
  row.names = FALSE
)

message("Completed. Results written to: ", output_dir)
print(
  point_data[, c(
    "SNP", "direct_HR", "indirect_shape_HR", "indirect_volume_HR",
    "indirect_total_HR", "total_HR", "snp_to_volume_t",
    "volume_to_log_hazard_z", "edf_beta_coordinate_1", "edf_beta_coordinate_2"
  )],
  row.names = FALSE, digits = 4
)
invisible(point_data)
}

if (sys.nframe() == 0L) {
  run_shape_volume_point_analysis()
}
