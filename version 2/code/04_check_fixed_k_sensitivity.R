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
if (!file.exists(model_script)) {
  stop("Model script not found: ", model_script)
}
source(model_script, local = FALSE)

code_dir <- file.path(project_dir, "code")
output_dir <- file.path(project_dir, "results", "fixed_k_sensitivity")
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

args <- commandArgs(trailingOnly = TRUE)
n_bootstrap <- if (length(args) >= 1L) as.integer(args[1L]) else 1000L
seed <- if (length(args) >= 2L) as.integer(args[2L]) else 20260905L
if (is.na(n_bootstrap) || n_bootstrap < 1L) {
  stop("The number of bootstrap replicates must be a positive integer.")
}
if (is.na(seed)) {
  stop("The random seed must be an integer.")
}

k_values <- c(8L, 10L, 12L, 15L, 18L)
effect_names <- c(
  "direct_log_HR",
  "indirect_log_HR_coordinate_1",
  "indirect_log_HR_coordinate_2",
  "indirect_log_HR",
  "total_log_HR"
)

subset_fpca <- function(fpca_fit, k) {
  if (k > ncol(fpca_fit$scores)) {
    stop("Requested K exceeds the fitted FPCA dimension.")
  }
  result <- fpca_fit
  result$scores <- fpca_fit$scores[, seq_len(k), drop = FALSE]
  result$eigenfunctions <- fpca_fit$eigenfunctions[, seq_len(k), drop = FALSE]
  result$n_components <- k
  result$retained_variance <- fpca_fit$cumulative_variance[k]
  result
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

fit_point_estimates <- function() {
  fpca_max <- fit_joint_fpca(M, grid, n_components = max(k_values))
  shape_fits <- setNames(
    lapply(snp_names, function(snp) fit_shape_layer(X[, snp], W, M, grid)),
    snp_names
  )
  rows <- vector("list", length(k_values) * length(snp_names))
  counter <- 1L

  for (k in k_values) {
    fpca_k <- subset_fpca(fpca_max, k)
    for (snp in snp_names) {
      cox_fit <- fit_cox_layer(X[, snp], W, subject, fpca_k)
      effects <- calculate_effects(shape_fits[[snp]], cox_fit, fpca_k)
      rows[[counter]] <- data.frame(
        K = k,
        SNP = snp,
        as.list(effects),
        direct_HR = exp(effects["direct_log_HR"]),
        indirect_HR = exp(effects["indirect_log_HR"]),
        total_HR = exp(effects["total_log_HR"]),
        cox_concordance = cox_fit$concordance,
        functional_variance_retained = fpca_k$retained_variance,
        cox_warning = paste(cox_fit$warnings, collapse = ""),
        check.names = FALSE
      )
      counter <- counter + 1L
    }
  }
  do.call(rbind, rows)
}

fit_bootstrap_replicate <- function(indices, replicate_id) {
  b_subject <- subject[indices, , drop = FALSE]
  b_X <- X[indices, , drop = FALSE]
  b_M <- M[indices, , , drop = FALSE]
  b_W <- W[indices, , drop = FALSE]
  b_fpca_max <- fit_joint_fpca(
    b_M,
    grid,
    n_components = max(k_values)
  )
  b_shape_fits <- setNames(
    lapply(
      snp_names,
      function(snp) fit_shape_layer(b_X[, snp], b_W, b_M, grid)
    ),
    snp_names
  )
  rows <- vector("list", length(k_values) * length(snp_names))
  counter <- 1L

  for (k in k_values) {
    b_fpca_k <- subset_fpca(b_fpca_max, k)
    for (snp in snp_names) {
      cox_fit <- fit_cox_layer(b_X[, snp], b_W, b_subject, b_fpca_k)
      effects <- calculate_effects(b_shape_fits[[snp]], cox_fit, b_fpca_k)
      rows[[counter]] <- data.frame(
        replicate = replicate_id,
        K = k,
        SNP = snp,
        as.list(effects),
        functional_variance_retained = b_fpca_k$retained_variance,
        cox_warning = paste(cox_fit$warnings, collapse = " | "),
        check.names = FALSE
      )
      counter <- counter + 1L
    }
  }
  rows
}

message("Fitting point estimates for fixed K values")
point_estimates <- fit_point_estimates()

bootstrap_rows <- vector(
  "list",
  n_bootstrap * length(k_values) * length(snp_names)
)
error_rows <- list()
row_counter <- 1L
error_counter <- 1L
set.seed(seed)

checkpoint_file <- file.path(
  output_dir,
  sprintf("fixed_k_checkpoint_B%d_seed%d.rds", n_bootstrap, seed)
)

for (replicate_id in seq_len(n_bootstrap)) {
  indices <- sample.int(n, size = n, replace = TRUE)
  replicate_result <- tryCatch(
    fit_bootstrap_replicate(indices, replicate_id),
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
    message("Completed fixed-K bootstrap replicate ", replicate_id, "/", n_bootstrap)
  }
}

bootstrap_rows <- bootstrap_rows[!vapply(bootstrap_rows, is.null, logical(1))]
if (length(bootstrap_rows) == 0L) {
  stop("All fixed-K bootstrap replicates failed.")
}
bootstrap_data <- do.call(rbind, bootstrap_rows)
errors <- if (length(error_rows) > 0L) {
  do.call(rbind, error_rows)
} else {
  data.frame(replicate = integer(), error = character())
}
warnings <- bootstrap_data[
  nzchar(bootstrap_data$cox_warning),
  c("replicate", "K", "SNP", "cox_warning"),
  drop = FALSE
]
valid_data <- bootstrap_data[!nzchar(bootstrap_data$cox_warning), , drop = FALSE]

summarize_effect <- function(k, snp, effect) {
  selected <- valid_data$K == k & valid_data$SNP == snp
  values <- valid_data[selected, effect]
  point <- point_estimates[
    point_estimates$K == k & point_estimates$SNP == snp,
    effect
  ]
  quantiles <- stats::quantile(values, c(0.025, 0.975), names = FALSE)
  data.frame(
    K = k,
    SNP = snp,
    effect = effect,
    point_log_HR = point,
    bootstrap_mean_log_HR = mean(values),
    bootstrap_bias_log_HR = mean(values) - point,
    bootstrap_SE_log_HR = stats::sd(values),
    CI_2.5_log_HR = quantiles[1L],
    CI_97.5_log_HR = quantiles[2L],
    point_HR = exp(point),
    CI_2.5_HR = exp(quantiles[1L]),
    CI_97.5_HR = exp(quantiles[2L]),
    successful_replicates = length(values)
  )
}

summary_rows <- list()
summary_counter <- 1L
for (k in k_values) {
  for (snp in snp_names) {
    for (effect in effect_names) {
      summary_rows[[summary_counter]] <- summarize_effect(k, snp, effect)
      summary_counter <- summary_counter + 1L
    }
  }
}
bootstrap_summary <- do.call(rbind, summary_rows)

bootstrap_data$direct_HR <- exp(bootstrap_data$direct_log_HR)
bootstrap_data$indirect_HR <- exp(bootstrap_data$indirect_log_HR)
bootstrap_data$total_HR <- exp(bootstrap_data$total_log_HR)

file_suffix <- sprintf("B%d_seed%d", n_bootstrap, seed)
utils::write.csv(
  point_estimates,
  file.path(output_dir, "fixed_k_point_estimates.csv"),
  row.names = FALSE
)
utils::write.csv(
  bootstrap_summary,
  file.path(output_dir, paste0("fixed_k_bootstrap_summary_", file_suffix, ".csv")),
  row.names = FALSE
)
utils::write.csv(
  bootstrap_data,
  file.path(output_dir, paste0("fixed_k_bootstrap_estimates_", file_suffix, ".csv")),
  row.names = FALSE
)
utils::write.csv(
  warnings,
  file.path(output_dir, paste0("fixed_k_bootstrap_warnings_", file_suffix, ".csv")),
  row.names = FALSE
)
utils::write.csv(
  errors,
  file.path(output_dir, paste0("fixed_k_bootstrap_errors_", file_suffix, ".csv")),
  row.names = FALSE
)

final_result <- list(
  point_estimates = point_estimates,
  bootstrap_summary = bootstrap_summary,
  bootstrap_estimates = bootstrap_data,
  warnings = warnings,
  errors = errors,
  settings = list(
    K_values = k_values,
    requested_replicates = n_bootstrap,
    completed_replicates = length(unique(bootstrap_data$replicate)),
    seed = seed,
    resampling_unit = "subject",
    confidence_interval = "95% percentile bootstrap",
    fpca_basis_recomputed_each_replicate = TRUE,
    K_fixed_within_each_model = TRUE,
    paired_resampling_across_K_and_SNP = TRUE
  )
)
saveRDS(
  final_result,
  file.path(output_dir, paste0("fixed_k_sensitivity_", file_suffix, ".rds")),
  compress = "xz"
)

if (file.exists(checkpoint_file)) {
  unlink(checkpoint_file)
}

message("Fixed-K sensitivity analysis completed: ", output_dir)
print(
  bootstrap_summary[bootstrap_summary$effect == "indirect_log_HR", ],
  row.names = FALSE
)
