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
wild_script <- file.path(code_dir, "07_wild_bootstrap_four_snp_mediation.R")
for (required_script in c(penalized_script, wild_script)) {
  if (!file.exists(required_script)) {
    stop("Required script not found: ", required_script)
  }
}
# Also sources 02_fit_four_snp_mediation_models.R.
source(penalized_script, local = FALSE)

# Reuse the multiplier weights and weighted first layer of the FLCRM wild
# bootstrap. Script 07 runs its bootstrap at top level, so only these two
# function definitions are evaluated from it.
for (script_expression in parse(wild_script)) {
  if (
    is.call(script_expression) &&
    identical(script_expression[[1L]], as.name("<-")) &&
    is.name(script_expression[[2L]]) &&
    as.character(script_expression[[2L]]) %in%
      c("normalize_multiplier_weights", "fit_shape_layer_weighted")
  ) {
    eval(script_expression, envir = globalenv())
  }
}

output_dir <- file.path(project_dir, "results", "penalized_functional_cox_wild_bootstrap")
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

args <- commandArgs(trailingOnly = TRUE)
n_bootstrap <- if (length(args) >= 1L) as.integer(args[1L]) else 10000L
seed <- if (length(args) >= 2L) as.integer(args[2L]) else 20260920L
if (is.na(n_bootstrap) || n_bootstrap < 1L) {
  stop("The number of wild bootstrap replicates must be a positive integer.")
}
if (is.na(seed)) {
  stop("The random seed must be an integer.")
}

effect_names <- c(
  "direct_log_HR",
  "indirect_log_HR_coordinate_1",
  "indirect_log_HR_coordinate_2",
  "indirect_log_HR",
  "total_log_HR"
)

run_penalized_wild_bootstrap <- function() {
# The smoothing parameters are selected once by REML on the observed data and
# held fixed in every replicate; the spline coefficients are re-estimated.
main_fits <- lapply(
  stats::setNames(snp_names, snp_names),
  function(snp) {
    fit <- fit_penalized_functional_cox(X[, snp], W, subject, M, grid, main_basis_size)
    list(fit = fit, ridge_form = penalized_cox_ridge_form(fit))
  }
)

point_rows <- lapply(snp_names, function(snp) {
  fit <- main_fits[[snp]]$fit
  alpha <- fit_shape_layer(X[, snp], W, M, grid)$alpha
  effects <- penalized_effects(alpha, fit$beta, fit$grid_weights)
  check <- fit_weighted_penalized_cox(
    X[, snp], W, subject, main_fits[[snp]]$ridge_form, rep(1, n)
  )
  if (max(abs(check$beta - fit$beta)) > 1e-5 || abs(check$tau - fit$tau) > 1e-5) {
    stop("The ridge form does not reproduce the REML fit for ", snp, ".")
  }
  data.frame(
    SNP = snp,
    direct_log_HR = fit$tau,
    as.list(effects),
    total_log_HR = fit$tau + effects[["indirect_log_HR"]]
  )
})
point_data <- do.call(rbind, point_rows)

set.seed(seed)
bootstrap_rows <- vector("list", n_bootstrap)
error_rows <- list()
for (replicate_id in seq_len(n_bootstrap)) {
  subject_weights <- normalize_multiplier_weights(stats::rexp(n, rate = 1), n)
  replicate_result <- tryCatch(
    do.call(rbind, lapply(snp_names, function(snp) {
      alpha <- fit_shape_layer_weighted(X[, snp], W, M, grid, subject_weights)$alpha
      cox_fit <- fit_weighted_penalized_cox(
        X[, snp], W, subject, main_fits[[snp]]$ridge_form, subject_weights
      )
      effects <- penalized_effects(alpha, cox_fit$beta, trapezoid_weights(grid))
      data.frame(
        replicate = replicate_id,
        SNP = snp,
        direct_log_HR = cox_fit$tau,
        as.list(effects),
        total_log_HR = cox_fit$tau + effects[["indirect_log_HR"]],
        cox_warning = paste(cox_fit$warnings, collapse = " | ")
      )
    })),
    error = function(e) e
  )
  if (inherits(replicate_result, "error")) {
    error_rows[[length(error_rows) + 1L]] <- data.frame(
      replicate = replicate_id,
      error = conditionMessage(replicate_result)
    )
  } else {
    bootstrap_rows[[replicate_id]] <- replicate_result
  }
  if (replicate_id %% 500L == 0L || replicate_id == n_bootstrap) {
    message("Completed wild bootstrap replicate ", replicate_id, "/", n_bootstrap)
  }
}

bootstrap_data <- do.call(rbind, bootstrap_rows[!vapply(bootstrap_rows, is.null, logical(1))])
if (is.null(bootstrap_data)) {
  stop("All wild bootstrap replicates failed.")
}
errors <- if (length(error_rows) > 0L) {
  do.call(rbind, error_rows)
} else {
  data.frame(replicate = integer(), error = character())
}
warning_data <- bootstrap_data[
  nzchar(bootstrap_data$cox_warning),
  c("replicate", "SNP", "cox_warning"),
  drop = FALSE
]
valid_data <- bootstrap_data[!nzchar(bootstrap_data$cox_warning), , drop = FALSE]

summary_rows <- list()
for (snp in snp_names) {
  for (effect in effect_names) {
    values <- valid_data[valid_data$SNP == snp, effect]
    point <- point_data[point_data$SNP == snp, effect]
    percentile_quantiles <- stats::quantile(values, probs = c(0.025, 0.975), names = FALSE)
    centered_quantiles <- stats::quantile(values - point, probs = c(0.025, 0.975), names = FALSE)
    basic_ci <- c(point - centered_quantiles[2L], point - centered_quantiles[1L])
    summary_rows[[length(summary_rows) + 1L]] <- data.frame(
      SNP = snp,
      effect = effect,
      point_log_HR = point,
      wild_mean_log_HR = mean(values),
      wild_bias_log_HR = mean(values) - point,
      wild_SE_log_HR = stats::sd(values),
      basic_CI_2.5_log_HR = basic_ci[1L],
      basic_CI_97.5_log_HR = basic_ci[2L],
      percentile_CI_2.5_log_HR = percentile_quantiles[1L],
      percentile_CI_97.5_log_HR = percentile_quantiles[2L],
      point_HR = exp(point),
      basic_CI_2.5_HR = exp(basic_ci[1L]),
      basic_CI_97.5_HR = exp(basic_ci[2L]),
      percentile_CI_2.5_HR = exp(percentile_quantiles[1L]),
      percentile_CI_97.5_HR = exp(percentile_quantiles[2L]),
      successful_replicates = length(values)
    )
  }
}
wild_summary <- do.call(rbind, summary_rows)

smoothing_table <- do.call(rbind, lapply(snp_names, function(snp) {
  fit <- main_fits[[snp]]$fit
  data.frame(
    SNP = snp,
    basis_size = fit$basis_size,
    smoothing_parameter_coordinate_1 = fit$smoothing_parameters[1L],
    smoothing_parameter_coordinate_2 = fit$smoothing_parameters[2L],
    edf_beta_coordinate_1 = fit$edf[1L],
    edf_beta_coordinate_2 = fit$edf[2L]
  )
}))

file_suffix <- sprintf("B%d_seed%d_k%d", n_bootstrap, seed, main_basis_size)
utils::write.csv(
  bootstrap_data,
  file.path(output_dir, paste0("wild_bootstrap_estimates_", file_suffix, ".csv")),
  row.names = FALSE
)
utils::write.csv(
  wild_summary,
  file.path(output_dir, paste0("wild_bootstrap_summary_", file_suffix, ".csv")),
  row.names = FALSE
)
utils::write.csv(
  errors,
  file.path(output_dir, paste0("wild_bootstrap_errors_", file_suffix, ".csv")),
  row.names = FALSE
)
utils::write.csv(
  warning_data,
  file.path(output_dir, paste0("wild_bootstrap_warnings_", file_suffix, ".csv")),
  row.names = FALSE
)
utils::write.csv(
  smoothing_table,
  file.path(output_dir, paste0("fixed_smoothing_parameters_", file_suffix, ".csv")),
  row.names = FALSE
)

message("Wild bootstrap completed. Results written to: ", output_dir)
print(
  wild_summary[
    wild_summary$effect %in% c("direct_log_HR", "indirect_log_HR", "total_log_HR"),
    c("SNP", "effect", "point_HR", "basic_CI_2.5_HR", "basic_CI_97.5_HR",
      "percentile_CI_2.5_HR", "percentile_CI_97.5_HR", "wild_SE_log_HR",
      "successful_replicates")
  ],
  row.names = FALSE, digits = 4
)
invisible(wild_summary)
}

if (sys.nframe() == 0L) {
  run_penalized_wild_bootstrap()
}
