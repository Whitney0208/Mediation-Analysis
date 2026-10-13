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

output_dir <- file.path(project_dir, "results", "penalized_functional_cox_bootstrap")
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

args <- commandArgs(trailingOnly = TRUE)
n_bootstrap <- if (length(args) >= 1L) as.integer(args[1L]) else 1000L
seed <- if (length(args) >= 2L) as.integer(args[2L]) else 20260905L
n_cores <- if (length(args) >= 3L) as.integer(args[3L]) else 1L
if (is.na(n_bootstrap) || n_bootstrap < 1L) {
  stop("The number of bootstrap replicates must be a positive integer.")
}
if (is.na(seed)) {
  stop("The random seed must be an integer.")
}
if (is.na(n_cores) || n_cores < 1L) {
  stop("The number of cores must be a positive integer.")
}

effect_names <- c(
  "direct_log_HR",
  "indirect_log_HR_coordinate_1",
  "indirect_log_HR_coordinate_2",
  "indirect_log_HR",
  "total_log_HR"
)

# Refits both layers on one set of subject rows. The smoothing parameters of the
# second layer are reselected by REML in every fit.
fit_penalized_replicate <- function(indices, replicate_id) {
  b_subject <- subject[indices, , drop = FALSE]
  b_X <- X[indices, , drop = FALSE]
  b_M <- M[indices, , , drop = FALSE]
  b_W <- W[indices, , drop = FALSE]
  rows <- lapply(snp_names, function(snp) {
    alpha <- fit_shape_layer(b_X[, snp], b_W, b_M, grid)$alpha
    fit_warnings <- character()
    cox_fit <- withCallingHandlers(
      fit_penalized_functional_cox(b_X[, snp], b_W, b_subject, b_M, grid, main_basis_size),
      warning = function(w) {
        fit_warnings <<- c(fit_warnings, conditionMessage(w))
        invokeRestart("muffleWarning")
      }
    )
    effects <- penalized_effects(alpha, cox_fit$beta, cox_fit$grid_weights)
    data.frame(
      replicate = replicate_id,
      SNP = snp,
      direct_log_HR = cox_fit$tau,
      as.list(effects),
      total_log_HR = cox_fit$tau + effects[["indirect_log_HR"]],
      edf_beta_coordinate_1 = cox_fit$edf[1L],
      edf_beta_coordinate_2 = cox_fit$edf[2L],
      cox_warning = paste(unique(fit_warnings), collapse = " | ")
    )
  })
  do.call(rbind, rows)
}

run_penalized_bootstrap <- function() {
point_data <- fit_penalized_replicate(seq_len(n), 0L)

# Indices are drawn in sequence exactly as in 03_bootstrap_four_snp_mediation.R,
# so replicate b uses the same subjects as replicate b of the main bootstrap.
set.seed(seed)
replicate_indices <- lapply(seq_len(n_bootstrap), function(replicate_id) {
  sample.int(n, size = n, replace = TRUE)
})

file_suffix <- sprintf("B%d_seed%d_k%d", n_bootstrap, seed, main_basis_size)
checkpoint_file <- file.path(output_dir, paste0("bootstrap_checkpoint_", file_suffix, ".rds"))
results <- if (file.exists(checkpoint_file)) readRDS(checkpoint_file) else list()
chunk_size <- 500L
chunk_starts <- seq(1L, n_bootstrap, by = chunk_size)
for (chunk_start in chunk_starts) {
  chunk <- chunk_start:min(chunk_start + chunk_size - 1L, n_bootstrap)
  if (all(as.character(chunk) %in% names(results))) {
    next
  }
  chunk_results <- parallel::mclapply(
    chunk,
    function(replicate_id) {
      tryCatch(
        fit_penalized_replicate(replicate_indices[[replicate_id]], replicate_id),
        error = function(e) e
      )
    },
    mc.cores = n_cores
  )
  names(chunk_results) <- as.character(chunk)
  results[names(chunk_results)] <- chunk_results
  saveRDS(results, checkpoint_file)
  message("Completed bootstrap replicate ", max(chunk), "/", n_bootstrap)
}

failed <- vapply(results, function(result) inherits(result, "error") || !is.data.frame(result), logical(1))
errors <- if (any(failed)) {
  data.frame(
    replicate = as.integer(names(results)[failed]),
    error = vapply(results[failed], function(result) {
      if (inherits(result, "error")) conditionMessage(result) else "Replicate returned no result."
    }, character(1))
  )
} else {
  data.frame(replicate = integer(), error = character())
}
bootstrap_data <- do.call(rbind, results[!failed])
if (is.null(bootstrap_data)) {
  stop("All bootstrap replicates failed.")
}
bootstrap_data <- bootstrap_data[order(bootstrap_data$replicate), , drop = FALSE]
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
    basic_ci <- 2 * point - rev(percentile_quantiles)
    summary_rows[[length(summary_rows) + 1L]] <- data.frame(
      SNP = snp,
      effect = effect,
      point_log_HR = point,
      bootstrap_mean_log_HR = mean(values),
      bootstrap_bias_log_HR = mean(values) - point,
      bootstrap_SE_log_HR = stats::sd(values),
      CI_2.5_log_HR = percentile_quantiles[1L],
      CI_97.5_log_HR = percentile_quantiles[2L],
      basic_CI_2.5_log_HR = basic_ci[1L],
      basic_CI_97.5_log_HR = basic_ci[2L],
      point_HR = exp(point),
      CI_2.5_HR = exp(percentile_quantiles[1L]),
      CI_97.5_HR = exp(percentile_quantiles[2L]),
      basic_CI_2.5_HR = exp(basic_ci[1L]),
      basic_CI_97.5_HR = exp(basic_ci[2L]),
      successful_replicates = length(values)
    )
  }
}
bootstrap_summary <- do.call(rbind, summary_rows)

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
if (file.exists(checkpoint_file)) {
  unlink(checkpoint_file)
}

message("Bootstrap completed. Results written to: ", output_dir)
print(
  bootstrap_summary[
    bootstrap_summary$effect %in% c("direct_log_HR", "indirect_log_HR", "total_log_HR"),
    c("SNP", "effect", "point_HR", "CI_2.5_HR", "CI_97.5_HR",
      "basic_CI_2.5_HR", "basic_CI_97.5_HR", "bootstrap_SE_log_HR",
      "successful_replicates")
  ],
  row.names = FALSE, digits = 4
)
invisible(bootstrap_summary)
}

if (sys.nframe() == 0L) {
  run_penalized_bootstrap()
}
