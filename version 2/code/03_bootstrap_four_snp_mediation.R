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
output_dir <- file.path(project_dir, "results", "four_snp_bootstrap")
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

if (!file.exists(model_script)) {
  stop("Model script not found: ", model_script)
}
if (!file.exists(point_file)) {
  stop("Point-estimate result not found: ", point_file)
}

source(model_script, local = FALSE)
point_result <- readRDS(point_file)
output_dir <- file.path(project_dir, "results", "four_snp_bootstrap")
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

args <- commandArgs(trailingOnly = TRUE)
n_bootstrap <- if (length(args) >= 1L) as.integer(args[1L]) else 100L
seed <- if (length(args) >= 2L) as.integer(args[2L]) else 20260905L
if (is.na(n_bootstrap) || n_bootstrap < 1L) {
  stop("The number of bootstrap replicates must be a positive integer.")
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

bootstrap_rows <- vector("list", n_bootstrap * length(snp_names))
error_rows <- list()
row_counter <- 1L
error_counter <- 1L
set.seed(seed)

checkpoint_file <- file.path(
  output_dir,
  sprintf("bootstrap_checkpoint_B%d_seed%d.rds", n_bootstrap, seed)
)

fit_one_bootstrap <- function(indices, replicate_id) {
  b_subject <- subject[indices, , drop = FALSE]
  b_X <- X[indices, , drop = FALSE]
  b_M <- M[indices, , , drop = FALSE]
  b_W <- W[indices, , drop = FALSE]
  b_fpca <- fit_joint_fpca(b_M, grid, variance_threshold = 0.85)
  replicate_rows <- vector("list", length(snp_names))

  for (snp_index in seq_along(snp_names)) {
    snp <- snp_names[snp_index]
    x <- b_X[, snp]
    shape_fit <- fit_shape_layer(x, b_W, b_M, grid)
    cox_fit <- fit_cox_layer(x, b_W, b_subject, b_fpca)
    alpha <- shape_fit$alpha
    beta <- cox_fit$beta
    indirect_1 <- sum(b_fpca$grid_weights * alpha[, 1L] * beta[, 1L])
    indirect_2 <- sum(b_fpca$grid_weights * alpha[, 2L] * beta[, 2L])
    indirect <- indirect_1 + indirect_2
    direct <- cox_fit$tau

    replicate_rows[[snp_index]] <- data.frame(
      replicate = replicate_id,
      SNP = snp,
      direct_log_HR = direct,
      indirect_log_HR_coordinate_1 = indirect_1,
      indirect_log_HR_coordinate_2 = indirect_2,
      indirect_log_HR = indirect,
      total_log_HR = direct + indirect,
      functional_components = b_fpca$n_components,
      functional_variance_retained = b_fpca$retained_variance,
      cox_warning = paste(cox_fit$warnings, collapse = " | ")
    )
  }

  replicate_rows
}

for (replicate_id in seq_len(n_bootstrap)) {
  indices <- sample.int(n, size = n, replace = TRUE)
  replicate_result <- tryCatch(
    fit_one_bootstrap(indices, replicate_id),
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
    checkpoint <- list(
      completed_replicates = replicate_id,
      requested_replicates = n_bootstrap,
      seed = seed,
      estimates = completed_data,
      errors = error_rows
    )
    saveRDS(checkpoint, checkpoint_file)
    message("Completed bootstrap replicate ", replicate_id, "/", n_bootstrap)
  }
}

bootstrap_rows <- bootstrap_rows[!vapply(bootstrap_rows, is.null, logical(1))]
bootstrap_data <- if (length(bootstrap_rows) > 0L) {
  do.call(rbind, bootstrap_rows)
} else {
  stop("All bootstrap replicates failed.")
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
valid_bootstrap_data <- bootstrap_data[!nzchar(bootstrap_data$cox_warning), , drop = FALSE]

summarize_effect <- function(snp, effect) {
  values <- valid_bootstrap_data[valid_bootstrap_data$SNP == snp, effect]
  point <- point_result$effects[point_result$effects$SNP == snp, effect]
  quantiles <- stats::quantile(values, probs = c(0.025, 0.975), names = FALSE)
  data.frame(
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

final_result <- list(
  summary = bootstrap_summary,
  estimates = bootstrap_data,
  errors = errors,
  warnings = warning_data,
  fpca_component_distribution = table(component_summary$functional_components),
  settings = list(
    requested_replicates = n_bootstrap,
    successful_replicates = length(unique(bootstrap_data$replicate)),
    failed_replicates = nrow(errors),
    warning_free_model_fits_by_snp = table(valid_bootstrap_data$SNP),
    seed = seed,
    resampling_unit = "subject",
    confidence_interval = "95% percentile bootstrap",
    fpca_recomputed_each_replicate = TRUE,
    four_snp_models_share_resampled_subject_indices = TRUE
  )
)

file_suffix <- sprintf("B%d_seed%d", n_bootstrap, seed)
utils::write.csv(
  bootstrap_data,
  file.path(output_dir, paste0("bootstrap_estimates_", file_suffix, ".csv")),
  row.names = FALSE
)
utils::write.csv(
  bootstrap_summary,
  file.path(output_dir, paste0("bootstrap_summary_", file_suffix, ".csv")),
  row.names = FALSE
)
utils::write.csv(
  errors,
  file.path(output_dir, paste0("bootstrap_errors_", file_suffix, ".csv")),
  row.names = FALSE
)
utils::write.csv(
  warning_data,
  file.path(output_dir, paste0("bootstrap_warnings_", file_suffix, ".csv")),
  row.names = FALSE
)
saveRDS(
  final_result,
  file.path(output_dir, paste0("bootstrap_results_", file_suffix, ".rds")),
  compress = "xz"
)

if (file.exists(checkpoint_file)) {
  unlink(checkpoint_file)
}

message("Bootstrap completed. Results written to: ", output_dir)
print(
  bootstrap_summary[
    bootstrap_summary$effect %in% c("direct_log_HR", "indirect_log_HR", "total_log_HR"),
  ],
  row.names = FALSE
)
