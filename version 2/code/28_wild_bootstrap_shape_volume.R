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
model_script <- file.path(code_dir, "22_fit_shape_volume_mediation.R")
wild_script <- file.path(code_dir, "07_wild_bootstrap_four_snp_mediation.R")
bca_script <- file.path(code_dir, "19_bca_penalized_functional_cox.R")
for (required_script in c(model_script, wild_script, bca_script)) {
  if (!file.exists(required_script)) {
    stop("Required script not found: ", required_script)
  }
}
# Also sources 16_fit_penalized_functional_cox.R and 02_fit_four_snp_mediation_models.R.
source(model_script, local = FALSE)

# Reuse the multiplier weights and weighted shape layer of script 07 and the
# BCa formula of script 19. Both scripts run analyses at top level, so only
# these function definitions are evaluated from them.
reused_functions <- c("normalize_multiplier_weights", "fit_shape_layer_weighted", "bca_interval")
for (script_file in c(wild_script, bca_script)) {
  for (script_expression in parse(script_file)) {
    if (
      is.call(script_expression) &&
      identical(script_expression[[1L]], as.name("<-")) &&
      is.name(script_expression[[2L]]) &&
      as.character(script_expression[[2L]]) %in% reused_functions
    ) {
      eval(script_expression, envir = globalenv())
    }
  }
}

output_dir <- file.path(project_dir, "results", "shape_volume_wild_bootstrap")
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

args <- commandArgs(trailingOnly = TRUE)
n_bootstrap <- if (length(args) >= 1L) as.integer(args[1L]) else 10000L
seed <- if (length(args) >= 2L) as.integer(args[2L]) else 20260920L
if (is.na(n_bootstrap) || n_bootstrap < 1L || is.na(seed)) {
  stop("Replicates and seed must be integers.")
}

effect_names <- c(
  "direct_log_HR",
  "indirect_shape_log_HR_coordinate_1",
  "indirect_shape_log_HR_coordinate_2",
  "indirect_shape_log_HR",
  "indirect_volume_log_HR",
  "indirect_total_log_HR",
  "total_log_HR"
)

# One weighted fit of both layers for one SNP. The smoothing parameters stay
# at their full-data REML values (ridge form of the penalized Cox model).
fit_weighted_shape_volume <- function(snp, ridge_form, rows, subject_weights) {
  x <- X[rows, snp]
  b_W <- W[rows, , drop = FALSE]
  b_volume <- volume_mediator[rows, , drop = FALSE]
  alpha <- fit_shape_layer_weighted(x, b_W, M[rows, , , drop = FALSE], grid, subject_weights)$alpha
  sqrt_weights <- sqrt(subject_weights / mean(subject_weights))
  volume_design <- cbind(SNP = as.numeric(x), b_W) * sqrt_weights
  snp_to_volume <- unname(qr.coef(qr(volume_design), b_volume * sqrt_weights)["SNP", 1L])
  reduced_form <- lapply(ridge_form, function(form) {
    form$penalized_columns <- form$penalized_columns[rows, , drop = FALSE]
    form$null_columns <- form$null_columns[rows, , drop = FALSE]
    form
  })
  cox_fit <- fit_weighted_penalized_cox(
    x, b_W, subject[rows, , drop = FALSE], reduced_form, subject_weights,
    scalar_mediators = b_volume
  )
  shape <- penalized_effects(alpha, cox_fit$beta, trapezoid_weights(grid))
  indirect_volume <- snp_to_volume * unname(cox_fit$mediator_coefficients[1L])
  indirect_total <- shape[["indirect_log_HR"]] + indirect_volume
  data.frame(
    SNP = snp,
    direct_log_HR = cox_fit$tau,
    indirect_shape_log_HR_coordinate_1 = shape[["indirect_log_HR_coordinate_1"]],
    indirect_shape_log_HR_coordinate_2 = shape[["indirect_log_HR_coordinate_2"]],
    indirect_shape_log_HR = shape[["indirect_log_HR"]],
    indirect_volume_log_HR = indirect_volume,
    indirect_total_log_HR = indirect_total,
    total_log_HR = cox_fit$tau + indirect_total,
    cox_warning = paste(cox_fit$warnings, collapse = " | ")
  )
}

run_shape_volume_wild_bootstrap <- function() {
ridge_forms <- lapply(stats::setNames(snp_names, snp_names), function(snp) {
  fit <- fit_penalized_functional_cox(
    X[, snp], W, subject, M, grid, main_basis_size, scalar_mediators = volume_mediator
  )
  penalized_cox_ridge_form(fit)
})

# Point estimates from the REML fit; the unit-weight ridge refit must match them.
point_data <- fit_shape_volume_replicate(seq_len(n), 0L)
for (snp in snp_names) {
  check <- fit_weighted_shape_volume(snp, ridge_forms[[snp]], seq_len(n), rep(1, n))
  reml <- point_data[point_data$SNP == snp, ]
  if (max(abs(unlist(check[effect_names]) - unlist(reml[effect_names]))) > 1e-5) {
    stop("The ridge form does not reproduce the REML fit for ", snp, ".")
  }
}

set.seed(seed)
bootstrap_rows <- vector("list", n_bootstrap)
error_rows <- list()
for (replicate_id in seq_len(n_bootstrap)) {
  subject_weights <- normalize_multiplier_weights(stats::rexp(n, rate = 1), n)
  replicate_result <- tryCatch(
    do.call(rbind, lapply(snp_names, function(snp) {
      fit_weighted_shape_volume(snp, ridge_forms[[snp]], seq_len(n), subject_weights)
    })),
    error = function(e) e
  )
  if (inherits(replicate_result, "error")) {
    error_rows[[length(error_rows) + 1L]] <- data.frame(
      replicate = replicate_id,
      error = conditionMessage(replicate_result)
    )
  } else {
    bootstrap_rows[[replicate_id]] <- cbind(replicate = replicate_id, replicate_result)
  }
  if (replicate_id %% 1000L == 0L || replicate_id == n_bootstrap) {
    message("Completed wild bootstrap replicate ", replicate_id, "/", n_bootstrap)
  }
}
bootstrap_data <- do.call(rbind, bootstrap_rows[!vapply(bootstrap_rows, is.null, logical(1))])
valid_data <- bootstrap_data[!nzchar(bootstrap_data$cox_warning), , drop = FALSE]

message("Running leave-one-subject-out jackknife with fixed smoothing parameters")
jackknife_data <- do.call(rbind, lapply(seq_len(n), function(left_out) {
  rows <- setdiff(seq_len(n), left_out)
  cbind(left_out = left_out, do.call(rbind, lapply(snp_names, function(snp) {
    fit_weighted_shape_volume(snp, ridge_forms[[snp]], rows, rep(1, length(rows)))
  })))
}))

summary_rows <- list()
for (snp in snp_names) {
  for (effect in effect_names) {
    values <- valid_data[valid_data$SNP == snp, effect]
    point <- point_data[point_data$SNP == snp, effect]
    bca <- bca_interval(values, point, jackknife_data[jackknife_data$SNP == snp, effect])
    percentile_quantiles <- stats::quantile(values, probs = c(0.025, 0.975), names = FALSE)
    basic_ci <- 2 * point - rev(percentile_quantiles)
    summary_rows[[length(summary_rows) + 1L]] <- data.frame(
      SNP = snp,
      effect = effect,
      point_log_HR = point,
      wild_SE_log_HR = stats::sd(values),
      bias_correction_z0 = bca$bias_correction,
      acceleration = bca$acceleration,
      point_HR = exp(point),
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
wild_summary <- do.call(rbind, summary_rows)
errors <- if (length(error_rows) > 0L) {
  do.call(rbind, error_rows)
} else {
  data.frame(replicate = integer(), error = character())
}

file_suffix <- sprintf("B%d_seed%d_k%d", n_bootstrap, seed, main_basis_size)
utils::write.csv(wild_summary, file.path(output_dir, paste0("wild_bootstrap_summary_", file_suffix, ".csv")), row.names = FALSE)
utils::write.csv(bootstrap_data, file.path(output_dir, paste0("wild_bootstrap_estimates_", file_suffix, ".csv")), row.names = FALSE)
utils::write.csv(jackknife_data, file.path(output_dir, sprintf("fixed_smoothing_jackknife_n%d_k%d.csv", n, main_basis_size)), row.names = FALSE)
utils::write.csv(errors, file.path(output_dir, paste0("wild_bootstrap_errors_", file_suffix, ".csv")), row.names = FALSE)
message(
  "Wild bootstrap completed: ", nrow(errors), " failed replicates, ",
  sum(nzchar(bootstrap_data$cox_warning)), " warning-flagged fits. Results written to: ", output_dir
)
print(
  wild_summary[
    wild_summary$effect %in% c("indirect_shape_log_HR", "indirect_volume_log_HR", "indirect_total_log_HR"),
    c("SNP", "effect", "point_HR", "basic_CI_2.5_HR", "basic_CI_97.5_HR",
      "percentile_CI_2.5_HR", "percentile_CI_97.5_HR", "BCa_CI_2.5_HR", "BCa_CI_97.5_HR")
  ],
  row.names = FALSE, digits = 4
)
invisible(wild_summary)
}

if (sys.nframe() == 0L) {
  run_shape_volume_wild_bootstrap()
}
