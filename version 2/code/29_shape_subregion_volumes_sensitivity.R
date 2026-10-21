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
bootstrap_script <- file.path(code_dir, "23_bootstrap_shape_volume_mediation.R")
if (!file.exists(bootstrap_script)) {
  stop("Shape and volume bootstrap script not found: ", bootstrap_script)
}
# Also sources 22, 16, and 02, and reads the same arguments: replicates, seed, cores.
source(bootstrap_script, local = FALSE)

# Sensitivity analysis: replace total CC volume / ICV with the five subregion
# volumes, each divided by ICV and standardized, alongside the shape mediator.
#   V_r  = a_r X + c_r' W + e_r,  r = 1..5
#   h(t) = h0(t) exp[tau X + kappa' W + sum_j int M_j beta_j + sum_r b_r V_r]
# Indirect through volume = sum_r a_r b_r; each a_r b_r is also reported.
region_names <- c(
  anterior = "CorpusCallosumAnterior",
  mid_anterior = "CorpusCallosumMidAnterior",
  central = "CorpusCallosumCentral",
  mid_posterior = "CorpusCallosumMidPosterior",
  posterior = "CorpusCallosumPosterior"
)
volume_mediator <- scale(as.matrix(volume_table[, region_names]) / volume_table$ICV)
colnames(volume_mediator) <- names(region_names)

fit_shape_volume_replicate <- function(indices, replicate_id) {
  b_subject <- subject[indices, , drop = FALSE]
  b_X <- X[indices, , drop = FALSE]
  b_M <- M[indices, , , drop = FALSE]
  b_W <- W[indices, , drop = FALSE]
  b_volume <- volume_mediator[indices, , drop = FALSE]
  rows <- lapply(snp_names, function(snp) {
    x <- b_X[, snp]
    alpha <- fit_shape_layer(x, b_W, b_M, grid)$alpha
    snp_to_volume <- qr.coef(qr(cbind(SNP = as.numeric(x), b_W)), b_volume)["SNP", ]
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
    region_indirect <- snp_to_volume * cox_fit$mediator_coefficients[colnames(b_volume)]
    shape <- penalized_effects(alpha, cox_fit$beta, cox_fit$grid_weights)
    indirect_volume <- sum(region_indirect)
    indirect_total <- shape[["indirect_log_HR"]] + indirect_volume
    data.frame(
      replicate = replicate_id,
      SNP = snp,
      direct_log_HR = cox_fit$tau,
      indirect_shape_log_HR_coordinate_1 = shape[["indirect_log_HR_coordinate_1"]],
      indirect_shape_log_HR_coordinate_2 = shape[["indirect_log_HR_coordinate_2"]],
      indirect_shape_log_HR = shape[["indirect_log_HR"]],
      indirect_volume_log_HR = indirect_volume,
      as.list(stats::setNames(region_indirect, paste0("indirect_", names(region_names), "_log_HR"))),
      indirect_total_log_HR = indirect_total,
      total_log_HR = cox_fit$tau + indirect_total,
      edf_beta_coordinate_1 = cox_fit$edf[1L],
      edf_beta_coordinate_2 = cox_fit$edf[2L],
      cox_warning = paste(unique(fit_warnings), collapse = " | ")
    )
  })
  do.call(rbind, rows)
}

effect_names <- c(
  "direct_log_HR",
  "indirect_shape_log_HR_coordinate_1",
  "indirect_shape_log_HR_coordinate_2",
  "indirect_shape_log_HR",
  "indirect_volume_log_HR",
  paste0("indirect_", names(region_names), "_log_HR"),
  "indirect_total_log_HR",
  "total_log_HR"
)

output_dir <- file.path(project_dir, "results", "shape_subregion_volumes_bootstrap")
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

if (sys.nframe() == 0L) {
  run_shape_volume_bootstrap()
}
