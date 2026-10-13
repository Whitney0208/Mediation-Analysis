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
bca_script <- file.path(code_dir, "19_bca_penalized_functional_cox.R")
for (required_script in c(penalized_script, bca_script)) {
  if (!file.exists(required_script)) {
    stop("Required script not found: ", required_script)
  }
}
# Also sources 02_fit_four_snp_mediation_models.R.
source(penalized_script, local = FALSE)

# Reuse the BCa formula of script 19 (same as script 08). Script 19 reads the
# pairs-bootstrap output at top level, so only bca_interval is evaluated from it.
for (script_expression in parse(bca_script)) {
  if (
    is.call(script_expression) &&
    identical(script_expression[[1L]], as.name("<-")) &&
    is.name(script_expression[[2L]]) &&
    as.character(script_expression[[2L]]) == "bca_interval"
  ) {
    eval(script_expression, envir = globalenv())
  }
}

wild_dir <- file.path(project_dir, "results", "penalized_functional_cox_wild_bootstrap")
output_dir <- wild_dir

args <- commandArgs(trailingOnly = TRUE)
n_bootstrap <- if (length(args) >= 1L) as.integer(args[1L]) else 10000L
seed <- if (length(args) >= 2L) as.integer(args[2L]) else 20260920L
if (is.na(n_bootstrap) || is.na(seed)) {
  stop("Replicates and seed must be integers.")
}
file_suffix <- sprintf("B%d_seed%d_k%d", n_bootstrap, seed, main_basis_size)
wild_file <- file.path(wild_dir, paste0("wild_bootstrap_estimates_", file_suffix, ".csv"))
if (!file.exists(wild_file)) {
  stop("Run 17_wild_bootstrap_penalized_functional_cox.R first: ", wild_file)
}

effect_names <- c(
  "direct_log_HR",
  "indirect_log_HR_coordinate_1",
  "indirect_log_HR_coordinate_2",
  "indirect_log_HR",
  "total_log_HR"
)

# Leave-one-subject-out fit with the smoothing parameters held at their REML
# values from the full data, matching the wild bootstrap.
fit_fixed_smoothing_jackknife <- function(snp, ridge_form, left_out) {
  keep <- setdiff(seq_len(n), left_out)
  reduced_form <- lapply(ridge_form, function(form) {
    form$penalized_columns <- form$penalized_columns[keep, , drop = FALSE]
    form$null_columns <- form$null_columns[keep, , drop = FALSE]
    form
  })
  alpha <- fit_shape_layer(X[keep, snp], W[keep, , drop = FALSE], M[keep, , , drop = FALSE], grid)$alpha
  cox_fit <- fit_weighted_penalized_cox(
    X[keep, snp], W[keep, , drop = FALSE], subject[keep, , drop = FALSE],
    reduced_form, rep(1, length(keep))
  )
  effects <- penalized_effects(alpha, cox_fit$beta, trapezoid_weights(grid))
  data.frame(
    left_out = left_out,
    SNP = snp,
    direct_log_HR = cox_fit$tau,
    as.list(effects),
    total_log_HR = cox_fit$tau + effects[["indirect_log_HR"]],
    cox_warning = paste(cox_fit$warnings, collapse = " | ")
  )
}

run_wild_bca <- function() {
wild_data <- utils::read.csv(wild_file, colClasses = c(cox_warning = "character"))
# read.csv turns the empty warning fields into NA.
wild_data$cox_warning[is.na(wild_data$cox_warning)] <- ""
valid_data <- wild_data[!nzchar(wild_data$cox_warning), , drop = FALSE]

summary_rows <- list()
jackknife_rows <- list()
for (snp in snp_names) {
  message("Jackknife for ", snp)
  full_fit <- fit_penalized_functional_cox(X[, snp], W, subject, M, grid, main_basis_size)
  ridge_form <- penalized_cox_ridge_form(full_fit)
  alpha <- fit_shape_layer(X[, snp], W, M, grid)$alpha
  effects <- penalized_effects(alpha, full_fit$beta, full_fit$grid_weights)
  point <- c(
    direct_log_HR = full_fit$tau,
    effects,
    total_log_HR = full_fit$tau + effects[["indirect_log_HR"]]
  )
  jackknife <- do.call(rbind, lapply(seq_len(n), function(left_out) {
    fit_fixed_smoothing_jackknife(snp, ridge_form, left_out)
  }))
  if (any(nzchar(jackknife$cox_warning))) {
    warning("Some jackknife fits for ", snp, " produced warnings.")
  }
  jackknife_rows[[snp]] <- jackknife

  for (effect in effect_names) {
    values <- valid_data[valid_data$SNP == snp, effect]
    bca <- bca_interval(values, point[[effect]], jackknife[[effect]])
    percentile_quantiles <- stats::quantile(values, probs = c(0.025, 0.975), names = FALSE)
    basic_ci <- 2 * point[[effect]] - rev(percentile_quantiles)
    summary_rows[[length(summary_rows) + 1L]] <- data.frame(
      SNP = snp,
      effect = effect,
      point_log_HR = point[[effect]],
      wild_SE_log_HR = stats::sd(values),
      bias_correction_z0 = bca$bias_correction,
      acceleration = bca$acceleration,
      BCa_lower_probability = bca$adjusted_probabilities[1L],
      BCa_upper_probability = bca$adjusted_probabilities[2L],
      point_HR = exp(point[[effect]]),
      basic_CI_2.5_HR = exp(basic_ci[1L]),
      basic_CI_97.5_HR = exp(basic_ci[2L]),
      percentile_CI_2.5_HR = exp(percentile_quantiles[1L]),
      percentile_CI_97.5_HR = exp(percentile_quantiles[2L]),
      BCa_CI_2.5_HR = exp(bca$interval[1L]),
      BCa_CI_97.5_HR = exp(bca$interval[2L]),
      successful_replicates = length(values)
    )
  }
}
bca_summary <- do.call(rbind, summary_rows)
jackknife_data <- do.call(rbind, jackknife_rows)

utils::write.csv(
  jackknife_data,
  file.path(output_dir, sprintf("fixed_smoothing_jackknife_n%d_k%d.csv", n, main_basis_size)),
  row.names = FALSE
)
utils::write.csv(
  bca_summary,
  file.path(output_dir, paste0("wild_bootstrap_bca_summary_", file_suffix, ".csv")),
  row.names = FALSE
)

message("Wild BCa completed. Results written to: ", output_dir)
print(
  bca_summary[
    bca_summary$effect %in% c("direct_log_HR", "indirect_log_HR", "total_log_HR"),
    c("SNP", "effect", "point_HR", "basic_CI_2.5_HR", "basic_CI_97.5_HR",
      "percentile_CI_2.5_HR", "percentile_CI_97.5_HR", "BCa_CI_2.5_HR", "BCa_CI_97.5_HR")
  ],
  row.names = FALSE, digits = 4
)
invisible(bca_summary)
}

if (sys.nframe() == 0L) {
  run_wild_bca()
}
