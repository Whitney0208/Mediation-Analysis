#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  if (!requireNamespace("mgcv", quietly = TRUE)) {
    stop("Package 'mgcv' is required.")
  }
  if (!requireNamespace("survival", quietly = TRUE)) {
    stop("Package 'survival' is required.")
  }
})

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
model_script <- file.path(code_dir, "02_fit_four_snp_mediation_models.R")
if (!file.exists(model_script)) {
  stop("Model script not found: ", model_script)
}
source(model_script, local = FALSE)

output_dir <- file.path(project_dir, "results", "penalized_functional_cox")
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

# Main analysis: cubic P-splines (second-order difference penalty) with 20 basis
# functions per coordinate; 10 and 30 are sensitivity checks.
main_basis_size <- 20L
basis_sizes <- c(10L, 20L, 30L)

# Second layer in the style of S-GMAS (refund::pfr with lf terms), adapted to a
# Cox outcome:
#   h(t) = h0(t) exp[tau X + kappa' W + sum_j integral M_j(s) beta_j(s) ds].
# Each beta_j(s) is a penalized spline; the integral uses trapezoid weights and
# the smoothing parameters are chosen by REML.
fit_penalized_functional_cox <- function(
  x,
  W_matrix,
  subject_data,
  shape_array,
  shape_grid,
  basis_size = main_basis_size
) {
  n_subjects <- dim(shape_array)[1L]
  p_grid <- length(shape_grid)
  weights_grid <- trapezoid_weights(shape_grid)
  W_cox <- W_matrix[, setdiff(colnames(W_matrix), "intercept"), drop = FALSE]
  grid_matrix <- matrix(shape_grid, n_subjects, p_grid, byrow = TRUE)
  cox_data <- list(
    time = as.numeric(subject_data$time),
    event = as.integer(subject_data$event),
    scalar = cbind(SNP = as.numeric(x), W_cox),
    s1 = grid_matrix,
    s2 = grid_matrix,
    m1 = sweep(shape_array[, , 1L], 2L, weights_grid, FUN = "*"),
    m2 = sweep(shape_array[, , 2L], 2L, weights_grid, FUN = "*")
  )
  fit <- mgcv::gam(
    time ~ scalar +
      s(s1, by = m1, bs = "ps", k = basis_size) +
      s(s2, by = m2, bs = "ps", k = basis_size),
    family = mgcv::cox.ph(),
    weights = event,
    data = cox_data,
    method = "REML"
  )
  coefficients <- stats::coef(fit)
  if (any(!is.finite(coefficients))) {
    stop("The penalized functional Cox model did not converge.")
  }

  smooth_parts <- lapply(fit$smooth, function(smooth_term) {
    columns <- smooth_term$first.para:smooth_term$last.para
    grid_data <- data.frame(shape_grid, 1)
    names(grid_data) <- c(smooth_term$term, smooth_term$by)
    basis_at_grid <- mgcv::PredictMat(smooth_term, grid_data)
    list(
      columns = columns,
      basis_at_grid = basis_at_grid,
      beta = as.numeric(basis_at_grid %*% coefficients[columns]),
      edf = sum(fit$edf[columns]),
      penalty = smooth_term$S[[1L]]
    )
  })

  list(
    fit = fit,
    tau = unname(coefficients["scalarSNP"]),
    beta = cbind(
      coordinate_1 = smooth_parts[[1L]]$beta,
      coordinate_2 = smooth_parts[[2L]]$beta
    ),
    edf = vapply(smooth_parts, `[[`, numeric(1L), "edf"),
    smoothing_parameters = unname(fit$sp),
    smooth_parts = smooth_parts,
    grid_weights = weights_grid,
    basis_size = basis_size
  )
}

# Rewrites a fitted model as an equivalent coxph ridge problem, which accepts
# subject weights. Each spline coefficient vector is split into the penalty null
# space (unpenalized columns) and the penalized space rescaled so that the
# penalty becomes a plain ridge with weight equal to the REML smoothing parameter.
penalized_cox_ridge_form <- function(penalized_fit) {
  model_matrix <- stats::model.matrix(penalized_fit$fit)
  lapply(seq_along(penalized_fit$smooth_parts), function(index) {
    part <- penalized_fit$smooth_parts[[index]]
    decomposition <- eigen(part$penalty, symmetric = TRUE)
    penalized <- decomposition$values > max(decomposition$values) * 1e-8
    to_penalized <- decomposition$vectors[, penalized, drop = FALSE] %*%
      diag(1 / sqrt(decomposition$values[penalized]))
    to_null <- decomposition$vectors[, !penalized, drop = FALSE]
    smooth_columns <- model_matrix[, part$columns, drop = FALSE]
    list(
      penalized_columns = smooth_columns %*% to_penalized,
      null_columns = smooth_columns %*% to_null,
      to_penalized = to_penalized,
      to_null = to_null,
      basis_at_grid = part$basis_at_grid,
      smoothing_parameter = penalized_fit$smoothing_parameters[index]
    )
  })
}

# Penalized Cox fit with subject weights at fixed smoothing parameters.
fit_weighted_penalized_cox <- function(
  x,
  W_matrix,
  subject_data,
  ridge_form,
  subject_weights
) {
  W_cox <- W_matrix[, setdiff(colnames(W_matrix), "intercept"), drop = FALSE]
  cox_data <- data.frame(
    time = as.numeric(subject_data$time),
    event = as.integer(subject_data$event)
  )
  cox_data$scalar <- cbind(SNP = as.numeric(x), W_cox)
  cox_data$null_1 <- ridge_form[[1L]]$null_columns
  cox_data$null_2 <- ridge_form[[2L]]$null_columns
  cox_data$penalized_1 <- ridge_form[[1L]]$penalized_columns
  cox_data$penalized_2 <- ridge_form[[2L]]$penalized_columns
  sp_1 <- ridge_form[[1L]]$smoothing_parameter
  sp_2 <- ridge_form[[2L]]$smoothing_parameter

  fit_warnings <- character()
  fit <- withCallingHandlers(
    survival::coxph(
      survival::Surv(time, event) ~ scalar + null_1 + null_2 +
        survival::ridge(penalized_1, theta = sp_1, scale = FALSE) +
        survival::ridge(penalized_2, theta = sp_2, scale = FALSE),
      data = cox_data,
      weights = subject_weights,
      ties = "breslow"
    ),
    warning = function(w) {
      fit_warnings <<- c(fit_warnings, conditionMessage(w))
      invokeRestart("muffleWarning")
    }
  )
  coefficients <- stats::coef(fit)
  if (any(!is.finite(coefficients))) {
    stop("The weighted penalized Cox model did not converge.")
  }

  beta <- vapply(seq_along(ridge_form), function(index) {
    form <- ridge_form[[index]]
    penalized_coef <- coefficients[grep(sprintf("penalized_%d", index), names(coefficients))]
    null_coef <- coefficients[grep(sprintf("^null_%d", index), names(coefficients))]
    spline_coef <- form$to_penalized %*% penalized_coef + form$to_null %*% null_coef
    as.numeric(form$basis_at_grid %*% spline_coef)
  }, numeric(nrow(ridge_form[[1L]]$basis_at_grid)))
  colnames(beta) <- c("coordinate_1", "coordinate_2")

  list(
    tau = unname(coefficients["scalarSNP"]),
    beta = beta,
    warnings = unique(fit_warnings)
  )
}

penalized_effects <- function(alpha, beta, weights_grid) {
  indirect_1 <- sum(weights_grid * alpha[, 1L] * beta[, 1L])
  indirect_2 <- sum(weights_grid * alpha[, 2L] * beta[, 2L])
  c(
    indirect_log_HR_coordinate_1 = indirect_1,
    indirect_log_HR_coordinate_2 = indirect_2,
    indirect_log_HR = indirect_1 + indirect_2
  )
}

run_penalized_point_analysis <- function() {
effect_rows <- list()
function_rows <- list()
for (snp in snp_names) {
  x <- X[, snp]
  alpha <- fit_shape_layer(x, W, M, grid)$alpha
  for (basis_size in basis_sizes) {
    message("Fitting ", snp, " with ", basis_size, " basis functions per coordinate")
    cox_fit <- fit_penalized_functional_cox(x, W, subject, M, grid, basis_size)
    effects <- penalized_effects(alpha, cox_fit$beta, cox_fit$grid_weights)
    direct <- cox_fit$tau
    effect_rows[[length(effect_rows) + 1L]] <- data.frame(
      SNP = snp,
      basis_size = basis_size,
      main_analysis = basis_size == main_basis_size,
      edf_beta_coordinate_1 = cox_fit$edf[1L],
      edf_beta_coordinate_2 = cox_fit$edf[2L],
      smoothing_parameter_coordinate_1 = cox_fit$smoothing_parameters[1L],
      smoothing_parameter_coordinate_2 = cox_fit$smoothing_parameters[2L],
      direct_log_HR = direct,
      as.list(effects),
      total_log_HR = direct + effects[["indirect_log_HR"]],
      direct_HR = exp(direct),
      indirect_HR = exp(effects[["indirect_log_HR"]]),
      total_HR = exp(direct + effects[["indirect_log_HR"]])
    )
    if (basis_size == main_basis_size) {
      function_rows[[length(function_rows) + 1L]] <- data.frame(
        SNP = snp,
        grid = grid,
        alpha_coordinate_1 = alpha[, 1L],
        alpha_coordinate_2 = alpha[, 2L],
        beta_coordinate_1 = cox_fit$beta[, 1L],
        beta_coordinate_2 = cox_fit$beta[, 2L]
      )
    }
  }
}
effect_summary <- do.call(rbind, effect_rows)
function_table <- do.call(rbind, function_rows)

utils::write.csv(
  effect_summary,
  file.path(output_dir, "penalized_functional_cox_effect_summary.csv"),
  row.names = FALSE
)
utils::write.csv(
  function_table,
  file.path(output_dir, sprintf("penalized_functional_cox_functions_k%d.csv", main_basis_size)),
  row.names = FALSE
)

message("Completed. Results written to: ", output_dir)
print(
  effect_summary[, c(
    "SNP", "basis_size", "edf_beta_coordinate_1", "edf_beta_coordinate_2",
    "indirect_log_HR_coordinate_1", "indirect_log_HR_coordinate_2",
    "direct_HR", "indirect_HR", "total_HR"
  )],
  row.names = FALSE, digits = 4
)
invisible(effect_summary)
}

if (sys.nframe() == 0L) {
  run_penalized_point_analysis()
}
