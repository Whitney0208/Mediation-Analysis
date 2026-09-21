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
output_dir <- file.path(project_dir, "results", "four_snp_bca_bootstrap")
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

if (!file.exists(model_script)) {
  stop("Model script not found: ", model_script)
}
source(model_script, local = FALSE)
output_dir <- file.path(project_dir, "results", "four_snp_bca_bootstrap")
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

args <- commandArgs(trailingOnly = TRUE)
n_bootstrap <- if (length(args) >= 1L) as.integer(args[1L]) else 5000L
seed <- if (length(args) >= 2L) as.integer(args[2L]) else 20260921L
fixed_k <- if (length(args) >= 3L) as.integer(args[3L]) else 18L
if (is.na(n_bootstrap) || n_bootstrap < 1L) {
  stop("The number of bootstrap replicates must be a positive integer.")
}
if (is.na(seed)) {
  stop("The random seed must be an integer.")
}
if (is.na(fixed_k) || fixed_k < 1L) {
  stop("The fixed FPCA dimension must be a positive integer.")
}

effect_names <- c(
  "direct_log_HR",
  "indirect_log_HR_coordinate_1",
  "indirect_log_HR_coordinate_2",
  "indirect_log_HR",
  "total_log_HR"
)

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

fit_indices <- function(indices, iteration_id, iteration_name) {
  selected_subject <- subject[indices, , drop = FALSE]
  selected_X <- X[indices, , drop = FALSE]
  selected_M <- M[indices, , , drop = FALSE]
  selected_W <- W[indices, , drop = FALSE]
  selected_fpca <- fit_joint_fpca(
    selected_M,
    grid,
    n_components = fixed_k
  )
  shape_fits <- setNames(
    lapply(
      snp_names,
      function(snp) {
        fit_shape_layer(selected_X[, snp], selected_W, selected_M, grid)
      }
    ),
    snp_names
  )
  rows <- vector("list", length(snp_names))

  for (snp_index in seq_along(snp_names)) {
    snp <- snp_names[snp_index]
    cox_fit <- fit_cox_layer(
      selected_X[, snp],
      selected_W,
      selected_subject,
      selected_fpca
    )
    effects <- calculate_effects(shape_fits[[snp]], cox_fit, selected_fpca)
    rows[[snp_index]] <- data.frame(
      iteration = iteration_id,
      iteration_type = iteration_name,
      SNP = snp,
      as.list(effects),
      functional_components = selected_fpca$n_components,
      functional_variance_retained = selected_fpca$retained_variance,
      cox_warning = paste(cox_fit$warnings, collapse = " | "),
      check.names = FALSE
    )
  }

  rows
}

message("Fitting the fixed-K point estimate")
point_rows <- fit_indices(seq_len(n), 0L, "point")
point_estimates <- do.call(rbind, point_rows)
if (any(nzchar(point_estimates$cox_warning))) {
  stop("The fixed-K point estimate generated a Cox warning.")
}

bootstrap_rows <- vector("list", n_bootstrap * length(snp_names))
bootstrap_error_rows <- list()
bootstrap_row_counter <- 1L
bootstrap_error_counter <- 1L
set.seed(seed)

bootstrap_checkpoint_file <- file.path(
  output_dir,
  sprintf("bca_bootstrap_checkpoint_B%d_seed%d_K%d.rds", n_bootstrap, seed, fixed_k)
)

for (replicate_id in seq_len(n_bootstrap)) {
  indices <- sample.int(n, size = n, replace = TRUE)
  replicate_result <- tryCatch(
    fit_indices(indices, replicate_id, "bootstrap"),
    error = function(e) e
  )

  if (inherits(replicate_result, "error")) {
    bootstrap_error_rows[[bootstrap_error_counter]] <- data.frame(
      replicate = replicate_id,
      error = conditionMessage(replicate_result)
    )
    bootstrap_error_counter <- bootstrap_error_counter + 1L
  } else {
    for (result_row in replicate_result) {
      bootstrap_rows[[bootstrap_row_counter]] <- result_row
      bootstrap_row_counter <- bootstrap_row_counter + 1L
    }
  }

  if (replicate_id %% 100L == 0L || replicate_id == n_bootstrap) {
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
        fixed_k = fixed_k,
        estimates = completed_data,
        errors = bootstrap_error_rows
      ),
      bootstrap_checkpoint_file
    )
    message("Completed BCa bootstrap replicate ", replicate_id, "/", n_bootstrap)
  }
}

bootstrap_rows <- bootstrap_rows[!vapply(bootstrap_rows, is.null, logical(1))]
if (length(bootstrap_rows) == 0L) {
  stop("All BCa bootstrap replicates failed.")
}
bootstrap_data <- do.call(rbind, bootstrap_rows)
bootstrap_errors <- if (length(bootstrap_error_rows) > 0L) {
  do.call(rbind, bootstrap_error_rows)
} else {
  data.frame(replicate = integer(), error = character())
}
bootstrap_warnings <- bootstrap_data[
  nzchar(bootstrap_data$cox_warning),
  c("iteration", "SNP", "cox_warning"),
  drop = FALSE
]
names(bootstrap_warnings)[1L] <- "replicate"
valid_bootstrap_data <- bootstrap_data[
  !nzchar(bootstrap_data$cox_warning),
  ,
  drop = FALSE
]

message("Running leave-one-subject-out jackknife")
jackknife_rows <- vector("list", n * length(snp_names))
jackknife_error_rows <- list()
jackknife_row_counter <- 1L
jackknife_error_counter <- 1L
jackknife_checkpoint_file <- file.path(
  output_dir,
  sprintf("bca_jackknife_checkpoint_n%d_K%d.rds", n, fixed_k)
)

for (left_out in seq_len(n)) {
  indices <- setdiff(seq_len(n), left_out)
  jackknife_result <- tryCatch(
    fit_indices(indices, left_out, "jackknife"),
    error = function(e) e
  )

  if (inherits(jackknife_result, "error")) {
    jackknife_error_rows[[jackknife_error_counter]] <- data.frame(
      left_out = left_out,
      RID = subject$RID[left_out],
      error = conditionMessage(jackknife_result)
    )
    jackknife_error_counter <- jackknife_error_counter + 1L
  } else {
    for (result_row in jackknife_result) {
      result_row$left_out_RID <- subject$RID[left_out]
      jackknife_rows[[jackknife_row_counter]] <- result_row
      jackknife_row_counter <- jackknife_row_counter + 1L
    }
  }

  if (left_out %% 25L == 0L || left_out == n) {
    completed_rows <- jackknife_rows[!vapply(jackknife_rows, is.null, logical(1))]
    completed_data <- if (length(completed_rows) > 0L) {
      do.call(rbind, completed_rows)
    } else {
      data.frame()
    }
    saveRDS(
      list(
        completed_leave_one_out_fits = left_out,
        requested_leave_one_out_fits = n,
        fixed_k = fixed_k,
        estimates = completed_data,
        errors = jackknife_error_rows
      ),
      jackknife_checkpoint_file
    )
    message("Completed jackknife fit ", left_out, "/", n)
  }
}

jackknife_rows <- jackknife_rows[!vapply(jackknife_rows, is.null, logical(1))]
jackknife_data <- if (length(jackknife_rows) > 0L) {
  do.call(rbind, jackknife_rows)
} else {
  stop("All jackknife fits failed.")
}
jackknife_errors <- if (length(jackknife_error_rows) > 0L) {
  do.call(rbind, jackknife_error_rows)
} else {
  data.frame(left_out = integer(), RID = integer(), error = character())
}
if (nrow(jackknife_errors) > 0L) {
  stop("BCa acceleration requires all leave-one-subject-out fits to succeed.")
}
jackknife_warnings <- jackknife_data[
  nzchar(jackknife_data$cox_warning),
  c("iteration", "left_out_RID", "SNP", "cox_warning"),
  drop = FALSE
]
names(jackknife_warnings)[1L] <- "left_out"

calculate_bca_summary <- function(snp, effect) {
  values <- valid_bootstrap_data[valid_bootstrap_data$SNP == snp, effect]
  point <- point_estimates[point_estimates$SNP == snp, effect]
  jackknife_values <- jackknife_data[jackknife_data$SNP == snp, effect]
  if (length(jackknife_values) != n || any(!is.finite(jackknife_values))) {
    stop("Incomplete jackknife values for ", snp, " and ", effect, ".")
  }

  proportion_less <- (
    sum(values < point) + 0.5 * sum(values == point)
  ) / length(values)
  probability_bound <- 1 / (2 * length(values))
  proportion_less <- min(
    max(proportion_less, probability_bound),
    1 - probability_bound
  )
  bias_correction <- stats::qnorm(proportion_less)

  jackknife_mean <- mean(jackknife_values)
  jackknife_influence <- jackknife_mean - jackknife_values
  acceleration_denominator <- 6 * sum(jackknife_influence^2)^(3 / 2)
  if (!is.finite(acceleration_denominator) || acceleration_denominator == 0) {
    stop("BCa acceleration is undefined for ", snp, " and ", effect, ".")
  }
  acceleration <- sum(jackknife_influence^3) / acceleration_denominator

  nominal_probabilities <- c(0.025, 0.975)
  nominal_z <- stats::qnorm(nominal_probabilities)
  adjusted_probabilities <- stats::pnorm(
    bias_correction +
      (bias_correction + nominal_z) /
        (1 - acceleration * (bias_correction + nominal_z))
  )
  quantile_bound <- 1 / (length(values) + 1)
  adjusted_probabilities <- pmin(
    pmax(adjusted_probabilities, quantile_bound),
    1 - quantile_bound
  )

  bca_quantiles <- stats::quantile(
    values,
    probs = adjusted_probabilities,
    names = FALSE,
    type = 7
  )
  percentile_quantiles <- stats::quantile(
    values,
    probs = nominal_probabilities,
    names = FALSE,
    type = 7
  )

  data.frame(
    SNP = snp,
    effect = effect,
    point_log_HR = point,
    bootstrap_mean_log_HR = mean(values),
    bootstrap_bias_log_HR = mean(values) - point,
    bootstrap_SE_log_HR = stats::sd(values),
    bias_correction_z0 = bias_correction,
    jackknife_acceleration = acceleration,
    adjusted_lower_probability = adjusted_probabilities[1L],
    adjusted_upper_probability = adjusted_probabilities[2L],
    BCa_CI_2.5_log_HR = bca_quantiles[1L],
    BCa_CI_97.5_log_HR = bca_quantiles[2L],
    percentile_CI_2.5_log_HR = percentile_quantiles[1L],
    percentile_CI_97.5_log_HR = percentile_quantiles[2L],
    point_HR = exp(point),
    BCa_CI_2.5_HR = exp(bca_quantiles[1L]),
    BCa_CI_97.5_HR = exp(bca_quantiles[2L]),
    percentile_CI_2.5_HR = exp(percentile_quantiles[1L]),
    percentile_CI_97.5_HR = exp(percentile_quantiles[2L]),
    successful_bootstrap_replicates = length(values),
    jackknife_replicates = length(jackknife_values)
  )
}

summary_rows <- list()
summary_counter <- 1L
for (snp in snp_names) {
  for (effect in effect_names) {
    summary_rows[[summary_counter]] <- calculate_bca_summary(snp, effect)
    summary_counter <- summary_counter + 1L
  }
}
bca_summary <- do.call(rbind, summary_rows)

bootstrap_data$direct_HR <- exp(bootstrap_data$direct_log_HR)
bootstrap_data$indirect_HR <- exp(bootstrap_data$indirect_log_HR)
bootstrap_data$total_HR <- exp(bootstrap_data$total_log_HR)

file_suffix <- sprintf("B%d_seed%d_K%d", n_bootstrap, seed, fixed_k)
utils::write.csv(
  point_estimates,
  file.path(output_dir, paste0("bca_point_estimates_", file_suffix, ".csv")),
  row.names = FALSE
)
utils::write.csv(
  bca_summary,
  file.path(output_dir, paste0("bca_summary_", file_suffix, ".csv")),
  row.names = FALSE
)
utils::write.csv(
  bootstrap_data,
  file.path(output_dir, paste0("bca_bootstrap_estimates_", file_suffix, ".csv")),
  row.names = FALSE
)
utils::write.csv(
  jackknife_data,
  file.path(output_dir, paste0("bca_jackknife_estimates_", file_suffix, ".csv")),
  row.names = FALSE
)
utils::write.csv(
  bootstrap_errors,
  file.path(output_dir, paste0("bca_bootstrap_errors_", file_suffix, ".csv")),
  row.names = FALSE
)
utils::write.csv(
  bootstrap_warnings,
  file.path(output_dir, paste0("bca_bootstrap_warnings_", file_suffix, ".csv")),
  row.names = FALSE
)
utils::write.csv(
  jackknife_errors,
  file.path(output_dir, paste0("bca_jackknife_errors_", file_suffix, ".csv")),
  row.names = FALSE
)
utils::write.csv(
  jackknife_warnings,
  file.path(output_dir, paste0("bca_jackknife_warnings_", file_suffix, ".csv")),
  row.names = FALSE
)

final_result <- list(
  point_estimates = point_estimates,
  summary = bca_summary,
  bootstrap_estimates = bootstrap_data,
  jackknife_estimates = jackknife_data,
  bootstrap_errors = bootstrap_errors,
  bootstrap_warnings = bootstrap_warnings,
  jackknife_errors = jackknife_errors,
  jackknife_warnings = jackknife_warnings,
  settings = list(
    requested_bootstrap_replicates = n_bootstrap,
    completed_bootstrap_replicates = length(unique(bootstrap_data$iteration)),
    failed_bootstrap_replicates = nrow(bootstrap_errors),
    warning_free_bootstrap_fits_by_snp = table(valid_bootstrap_data$SNP),
    completed_jackknife_replicates = length(unique(jackknife_data$iteration)),
    seed = seed,
    fixed_fpca_components = fixed_k,
    resampling_unit = "subject",
    confidence_interval = "95% bias-corrected and accelerated bootstrap",
    quantile_type = 7L,
    fpca_basis_recomputed_each_bootstrap_and_jackknife_fit = TRUE,
    four_snp_models_share_resampled_subject_indices = TRUE
  )
)
saveRDS(
  final_result,
  file.path(output_dir, paste0("bca_results_", file_suffix, ".rds")),
  compress = "xz"
)

if (file.exists(bootstrap_checkpoint_file)) {
  unlink(bootstrap_checkpoint_file)
}
if (file.exists(jackknife_checkpoint_file)) {
  unlink(jackknife_checkpoint_file)
}

message("BCa bootstrap completed. Results written to: ", output_dir)
print(
  bca_summary[
    bca_summary$effect %in% c(
      "direct_log_HR", "indirect_log_HR", "total_log_HR"
    ),
  ],
  row.names = FALSE
)
