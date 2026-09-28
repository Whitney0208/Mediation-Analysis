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
model_script <- file.path(code_dir, "02_fit_four_snp_mediation_models.R")
point_file <- file.path(
  project_dir, "results", "four_snp_point_estimates", "four_snp_model_results.rds"
)
output_dir <- file.path(project_dir, "results", "four_snp_wild_bootstrap")
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

if (!file.exists(model_script)) {
  stop("Model script not found: ", model_script)
}
if (!file.exists(point_file)) {
  stop("Point-estimate result not found: ", point_file)
}

source(model_script, local = FALSE)
point_result <- readRDS(point_file)
output_dir <- file.path(project_dir, "results", "four_snp_wild_bootstrap")
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

normalize_multiplier_weights <- function(weights, expected_length) {
  weights <- as.numeric(weights)
  if (
    length(weights) != expected_length || any(!is.finite(weights)) ||
    any(weights <= 0)
  ) {
    stop("Multiplier weights must be finite, positive, and subject-specific.")
  }
  weights / mean(weights)
}

fit_shape_layer_weighted <- function(
  x,
  W_matrix,
  shape_array,
  shape_grid,
  subject_weights
) {
  p_grid <- length(shape_grid)
  design <- cbind(SNP = as.numeric(x), W_matrix)
  subject_weights <- normalize_multiplier_weights(subject_weights, nrow(design))
  sqrt_weights <- sqrt(subject_weights)
  weighted_design <- design * sqrt_weights
  design_qr <- qr(weighted_design)
  if (design_qr$rank < ncol(design)) {
    stop("The weighted first-stage design matrix is rank deficient.")
  }

  shape_matrix <- cbind(shape_array[, , 1L], shape_array[, , 2L])
  weighted_shape <- shape_matrix * sqrt_weights
  raw_coef <- qr.coef(design_qr, weighted_shape)
  rownames(raw_coef) <- colnames(design)
  smoothed <- smooth_shape_coefficients(raw_coef, shape_grid, p_grid)

  raw_array <- array(
    NA_real_,
    dim = c(nrow(raw_coef), p_grid, 2L),
    dimnames = dimnames(smoothed$coefficients)
  )
  raw_array[, , 1L] <- raw_coef[, seq_len(p_grid), drop = FALSE]
  raw_array[, , 2L] <- raw_coef[, p_grid + seq_len(p_grid), drop = FALSE]

  list(
    design_names = colnames(design),
    coefficients_raw = raw_array,
    coefficients_smoothed = smoothed$coefficients,
    smoothing = smoothed$smoothing,
    alpha = smoothed$coefficients["SNP", , , drop = TRUE]
  )
}

fit_joint_fpca_weighted <- function(
  shape_array,
  shape_grid,
  subject_weights,
  variance_threshold = 0.85,
  n_components = NULL
) {
  p_grid <- length(shape_grid)
  subject_weights <- normalize_multiplier_weights(subject_weights, dim(shape_array)[1L])
  weights_grid <- trapezoid_weights(shape_grid)
  weights_flat <- rep(weights_grid, times = 2L)
  shape_matrix <- cbind(shape_array[, , 1L], shape_array[, , 2L])
  shape_mean <- colSums(shape_matrix * subject_weights) / sum(subject_weights)
  centered <- sweep(shape_matrix, 2L, shape_mean, FUN = "-")
  spatially_weighted <- sweep(centered, 2L, sqrt(weights_flat), FUN = "*")
  decomposition <- svd(spatially_weighted * sqrt(subject_weights), nu = 0L)

  effective_denominator <- sum(subject_weights) -
    sum(subject_weights^2) / sum(subject_weights)
  eigenvalues <- decomposition$d^2 / effective_denominator
  variance_fraction <- eigenvalues / sum(eigenvalues)
  cumulative_variance <- cumsum(variance_fraction)
  if (is.null(n_components)) {
    n_components <- which(cumulative_variance >= variance_threshold)[1L]
  } else {
    n_components <- as.integer(n_components)
    if (
      length(n_components) != 1L || is.na(n_components) ||
      n_components < 1L || n_components > ncol(decomposition$v)
    ) {
      stop("n_components is outside the available weighted FPCA range.")
    }
  }

  weighted_eigenvectors <- decomposition$v[, seq_len(n_components), drop = FALSE]
  eigenfunctions <- sweep(
    weighted_eigenvectors,
    1L,
    sqrt(weights_flat),
    FUN = "/"
  )
  # Subject weights estimate the basis; scores remain projections of observed shapes.
  scores <- spatially_weighted %*% weighted_eigenvectors
  colnames(scores) <- sprintf("FPC%02d", seq_len(n_components))

  list(
    scores = scores,
    eigenfunctions = eigenfunctions,
    eigenvalues = eigenvalues,
    variance_fraction = variance_fraction,
    cumulative_variance = cumulative_variance,
    n_components = n_components,
    retained_variance = cumulative_variance[n_components],
    shape_mean = shape_mean,
    grid_weights = weights_grid,
    flat_weights = weights_flat,
    p_grid = p_grid
  )
}

fit_cox_layer_weighted <- function(
  x,
  W_matrix,
  subject_data,
  fpca_fit,
  subject_weights
) {
  subject_weights <- normalize_multiplier_weights(subject_weights, nrow(W_matrix))
  W_cox <- W_matrix[, setdiff(colnames(W_matrix), "intercept"), drop = FALSE]
  cox_data <- data.frame(
    time = as.numeric(subject_data$time),
    event = as.integer(subject_data$event),
    SNP = as.numeric(x),
    W_cox,
    fpca_fit$scores,
    check.names = FALSE
  )
  fit_warnings <- character()
  fit <- withCallingHandlers(
    survival::coxph(
      survival::Surv(time, event) ~ .,
      data = cox_data,
      weights = subject_weights,
      ties = "breslow",
      robust = FALSE,
      x = TRUE,
      y = TRUE,
      model = TRUE
    ),
    warning = function(w) {
      fit_warnings <<- c(fit_warnings, conditionMessage(w))
      invokeRestart("muffleWarning")
    }
  )
  if (any(!is.finite(stats::coef(fit)))) {
    stop("The weighted Cox model did not converge to finite coefficients.")
  }

  fpc_names <- colnames(fpca_fit$scores)
  theta <- stats::coef(fit)[fpc_names]
  beta_flat <- as.numeric(fpca_fit$eigenfunctions %*% theta)
  p_grid <- fpca_fit$p_grid
  beta <- cbind(
    coordinate_1 = beta_flat[seq_len(p_grid)],
    coordinate_2 = beta_flat[p_grid + seq_len(p_grid)]
  )

  list(
    tau = unname(stats::coef(fit)["SNP"]),
    beta = beta,
    concordance = unname(summary(fit)$concordance[1L]),
    warnings = unique(fit_warnings)
  )
}

calculate_effects <- function(shape_fit, cox_fit, fpca_fit) {
  alpha <- shape_fit$alpha
  beta <- cox_fit$beta
  indirect_1 <- sum(fpca_fit$grid_weights * alpha[, 1L] * beta[, 1L])
  indirect_2 <- sum(fpca_fit$grid_weights * alpha[, 2L] * beta[, 2L])
  indirect <- indirect_1 + indirect_2
  direct <- cox_fit$tau
  c(
    direct_log_HR = direct,
    indirect_log_HR_coordinate_1 = indirect_1,
    indirect_log_HR_coordinate_2 = indirect_2,
    indirect_log_HR = indirect,
    total_log_HR = direct + indirect
  )
}

fit_one_wild_bootstrap <- function(subject_weights, replicate_id) {
  subject_weights <- normalize_multiplier_weights(subject_weights, n)
  weighted_fpca <- fit_joint_fpca_weighted(
    M,
    grid,
    subject_weights,
    variance_threshold = 0.85
  )
  effective_sample_size <- sum(subject_weights)^2 / sum(subject_weights^2)
  replicate_rows <- vector("list", length(snp_names))

  for (snp_index in seq_along(snp_names)) {
    snp <- snp_names[snp_index]
    shape_fit <- fit_shape_layer_weighted(
      X[, snp], W, M, grid, subject_weights
    )
    cox_fit <- fit_cox_layer_weighted(
      X[, snp], W, subject, weighted_fpca, subject_weights
    )
    effects <- calculate_effects(shape_fit, cox_fit, weighted_fpca)
    replicate_rows[[snp_index]] <- data.frame(
      replicate = replicate_id,
      SNP = snp,
      as.list(effects),
      functional_components = weighted_fpca$n_components,
      functional_variance_retained = weighted_fpca$retained_variance,
      multiplier_effective_sample_size = effective_sample_size,
      multiplier_min = min(subject_weights),
      multiplier_max = max(subject_weights),
      cox_warning = paste(cox_fit$warnings, collapse = " | "),
      check.names = FALSE
    )
  }

  replicate_rows
}

bootstrap_rows <- vector("list", n_bootstrap * length(snp_names))
error_rows <- list()
row_counter <- 1L
error_counter <- 1L
set.seed(seed)

checkpoint_file <- file.path(
  output_dir,
  sprintf("wild_bootstrap_checkpoint_B%d_seed%d.rds", n_bootstrap, seed)
)

for (replicate_id in seq_len(n_bootstrap)) {
  subject_weights <- stats::rexp(n, rate = 1)
  replicate_result <- tryCatch(
    fit_one_wild_bootstrap(subject_weights, replicate_id),
    error = function(e) e
  )

  if (inherits(replicate_result, "error")) {
    error_rows[[error_counter]] <- data.frame(
      replicate = replicate_id,
      error = conditionMessage(replicate_result)
    )
    error_counter <- error_counter + 1L
  } else {
    for (result_row in replicate_result) {
      bootstrap_rows[[row_counter]] <- result_row
      row_counter <- row_counter + 1L
    }
  }

  if (replicate_id %% 25L == 0L || replicate_id == n_bootstrap) {
    completed_rows <- bootstrap_rows[!vapply(bootstrap_rows, is.null, logical(1))]
    completed_data <- if (length(completed_rows) > 0L) {
      do.call(rbind, completed_rows)
    } else {
      data.frame()
    }
    saveRDS(
      list(
        completed_replicates = replicate_id,
        requested_replicates = n_bootstrap,
        seed = seed,
        estimates = completed_data,
        errors = error_rows
      ),
      checkpoint_file
    )
    message("Completed wild bootstrap replicate ", replicate_id, "/", n_bootstrap)
  }
}

bootstrap_rows <- bootstrap_rows[!vapply(bootstrap_rows, is.null, logical(1))]
if (length(bootstrap_rows) == 0L) {
  stop("All wild bootstrap replicates failed.")
}
bootstrap_data <- do.call(rbind, bootstrap_rows)
errors <- if (length(error_rows) > 0L) {
  do.call(rbind, error_rows)
} else {
  data.frame(replicate = integer(), error = character())
}
warnings <- bootstrap_data[
  nzchar(bootstrap_data$cox_warning),
  c("replicate", "SNP", "cox_warning"),
  drop = FALSE
]
valid_data <- bootstrap_data[!nzchar(bootstrap_data$cox_warning), , drop = FALSE]

summarize_effect <- function(snp, effect) {
  values <- valid_data[valid_data$SNP == snp, effect]
  point <- point_result$effects[point_result$effects$SNP == snp, effect]
  percentile_quantiles <- stats::quantile(
    values,
    probs = c(0.025, 0.975),
    names = FALSE
  )
  centered_quantiles <- stats::quantile(
    values - point,
    probs = c(0.025, 0.975),
    names = FALSE
  )
  basic_ci <- c(
    point - centered_quantiles[2L],
    point - centered_quantiles[1L]
  )
  data.frame(
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

summary_rows <- list()
summary_counter <- 1L
for (snp in snp_names) {
  for (effect in effect_names) {
    summary_rows[[summary_counter]] <- summarize_effect(snp, effect)
    summary_counter <- summary_counter + 1L
  }
}
bootstrap_summary <- do.call(rbind, summary_rows)

bootstrap_data$direct_HR <- exp(bootstrap_data$direct_log_HR)
bootstrap_data$indirect_HR <- exp(bootstrap_data$indirect_log_HR)
bootstrap_data$total_HR <- exp(bootstrap_data$total_log_HR)

component_summary <- aggregate(
  functional_components ~ replicate,
  data = bootstrap_data,
  FUN = function(x) x[1L]
)
weight_diagnostics <- unique(bootstrap_data[c(
  "replicate",
  "multiplier_effective_sample_size",
  "multiplier_min",
  "multiplier_max"
)])

final_result <- list(
  summary = bootstrap_summary,
  estimates = bootstrap_data,
  errors = errors,
  warnings = warnings,
  weight_diagnostics = weight_diagnostics,
  fpca_component_distribution = table(component_summary$functional_components),
  settings = list(
    requested_replicates = n_bootstrap,
    successful_replicates = length(unique(bootstrap_data$replicate)),
    failed_replicates = nrow(errors),
    warning_free_model_fits_by_snp = table(valid_data$SNP),
    seed = seed,
    resampling_method = "subject-level exponential multiplier wild bootstrap",
    multiplier_distribution = "independent Exp(1), normalized to mean 1",
    confidence_interval_primary = "95% centered basic multiplier interval",
    confidence_interval_secondary = "95% weighted-bootstrap percentile interval",
    same_subject_weights_across_shape_fpca_cox_and_snp_models = TRUE,
    weighted_fpca_recomputed_each_replicate = TRUE,
    observed_subjects_and_risk_sets_retained = TRUE
  )
)

file_suffix <- sprintf("B%d_seed%d", n_bootstrap, seed)
utils::write.csv(
  bootstrap_data,
  file.path(output_dir, paste0("wild_bootstrap_estimates_", file_suffix, ".csv")),
  row.names = FALSE
)
utils::write.csv(
  bootstrap_summary,
  file.path(output_dir, paste0("wild_bootstrap_summary_", file_suffix, ".csv")),
  row.names = FALSE
)
utils::write.csv(
  errors,
  file.path(output_dir, paste0("wild_bootstrap_errors_", file_suffix, ".csv")),
  row.names = FALSE
)
utils::write.csv(
  warnings,
  file.path(output_dir, paste0("wild_bootstrap_warnings_", file_suffix, ".csv")),
  row.names = FALSE
)
utils::write.csv(
  weight_diagnostics,
  file.path(output_dir, paste0("wild_bootstrap_weight_diagnostics_", file_suffix, ".csv")),
  row.names = FALSE
)
saveRDS(
  final_result,
  file.path(output_dir, paste0("wild_bootstrap_results_", file_suffix, ".rds")),
  compress = "xz"
)

if (file.exists(checkpoint_file)) {
  unlink(checkpoint_file)
}

message("Wild bootstrap completed. Results written to: ", output_dir)
print(
  bootstrap_summary[
    bootstrap_summary$effect %in% c(
      "direct_log_HR", "indirect_log_HR", "total_log_HR"
    ),
  ],
  row.names = FALSE
)
