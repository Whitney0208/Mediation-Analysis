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
bca_script <- file.path(code_dir, "19_bca_penalized_functional_cox.R")
for (required_script in c(model_script, bca_script)) {
  if (!file.exists(required_script)) {
    stop("Required script not found: ", required_script)
  }
}
# Also sources 16_fit_penalized_functional_cox.R and 02_fit_four_snp_mediation_models.R.
source(model_script, local = FALSE)

# Reuse the BCa formula of script 19 (same as script 08). Script 19 reads the
# shape-only bootstrap output at top level, so only bca_interval is evaluated.
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

output_dir <- file.path(project_dir, "results", "shape_volume_mediation_bootstrap")
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

args <- commandArgs(trailingOnly = TRUE)
n_bootstrap <- if (length(args) >= 1L) as.integer(args[1L]) else 1000L
seed <- if (length(args) >= 2L) as.integer(args[2L]) else 20260905L
n_cores <- if (length(args) >= 3L) as.integer(args[3L]) else 1L
if (is.na(n_bootstrap) || n_bootstrap < 1L || is.na(seed) || is.na(n_cores) || n_cores < 1L) {
  stop("Replicates, seed, and cores must be positive integers.")
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

run_shape_volume_bootstrap <- function() {
point_data <- fit_shape_volume_replicate(seq_len(n), 0L)

# Indices are drawn in sequence exactly as in 03_bootstrap_four_snp_mediation.R
# and 18_bootstrap_penalized_functional_cox.R.
set.seed(seed)
replicate_indices <- lapply(seq_len(n_bootstrap), function(replicate_id) {
  sample.int(n, size = n, replace = TRUE)
})

file_suffix <- sprintf("B%d_seed%d_k%d", n_bootstrap, seed, main_basis_size)
checkpoint_file <- file.path(output_dir, paste0("bootstrap_checkpoint_", file_suffix, ".rds"))
results <- if (file.exists(checkpoint_file)) readRDS(checkpoint_file) else list()
chunk_size <- 500L
for (chunk_start in seq(1L, n_bootstrap, by = chunk_size)) {
  chunk <- chunk_start:min(chunk_start + chunk_size - 1L, n_bootstrap)
  if (all(as.character(chunk) %in% names(results))) {
    next
  }
  chunk_results <- parallel::mclapply(
    chunk,
    function(replicate_id) {
      tryCatch(
        fit_shape_volume_replicate(replicate_indices[[replicate_id]], replicate_id),
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

failed <- vapply(results, function(result) !is.data.frame(result), logical(1))
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
valid_data <- bootstrap_data[!nzchar(bootstrap_data$cox_warning), , drop = FALSE]

message("Running leave-one-subject-out jackknife")
jackknife_results <- parallel::mclapply(
  seq_len(n),
  function(left_out) {
    tryCatch(
      fit_shape_volume_replicate(setdiff(seq_len(n), left_out), left_out),
      error = function(e) e
    )
  },
  mc.cores = n_cores
)
if (any(vapply(jackknife_results, inherits, logical(1), what = "error"))) {
  stop("BCa acceleration requires all leave-one-subject-out fits to succeed.")
}
jackknife_data <- do.call(rbind, jackknife_results)
names(jackknife_data)[names(jackknife_data) == "replicate"] <- "left_out"

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
      bootstrap_mean_log_HR = mean(values),
      bootstrap_SE_log_HR = stats::sd(values),
      bias_correction_z0 = bca$bias_correction,
      acceleration = bca$acceleration,
      point_HR = exp(point),
      BCa_CI_2.5_HR = exp(bca$interval[1L]),
      BCa_CI_97.5_HR = exp(bca$interval[2L]),
      percentile_CI_2.5_HR = exp(percentile_quantiles[1L]),
      percentile_CI_97.5_HR = exp(percentile_quantiles[2L]),
      basic_CI_2.5_HR = exp(basic_ci[1L]),
      basic_CI_97.5_HR = exp(basic_ci[2L]),
      successful_replicates = length(values)
    )
  }
}
bootstrap_summary <- do.call(rbind, summary_rows)
warning_data <- bootstrap_data[
  nzchar(bootstrap_data$cox_warning),
  c("replicate", "SNP", "cox_warning"),
  drop = FALSE
]

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
  jackknife_data,
  file.path(output_dir, sprintf("jackknife_estimates_n%d_k%d.csv", n, main_basis_size)),
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
    bootstrap_summary$effect %in% c(
      "direct_log_HR", "indirect_shape_log_HR", "indirect_volume_log_HR",
      "indirect_total_log_HR", "total_log_HR"
    ),
    c("SNP", "effect", "point_HR", "BCa_CI_2.5_HR", "BCa_CI_97.5_HR",
      "percentile_CI_2.5_HR", "percentile_CI_97.5_HR", "successful_replicates")
  ],
  row.names = FALSE, digits = 4
)
invisible(bootstrap_summary)
}

if (sys.nframe() == 0L) {
  run_shape_volume_bootstrap()
}
