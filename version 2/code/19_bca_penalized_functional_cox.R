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
bootstrap_script <- file.path(code_dir, "18_bootstrap_penalized_functional_cox.R")
if (!file.exists(bootstrap_script)) {
  stop("Penalized bootstrap script not found: ", bootstrap_script)
}
# Also sources 16_fit_penalized_functional_cox.R and 02_fit_four_snp_mediation_models.R,
# and reads the same arguments: replicates, seed, cores.
source(bootstrap_script, local = FALSE)

bootstrap_dir <- file.path(project_dir, "results", "penalized_functional_cox_bootstrap")
output_dir <- file.path(project_dir, "results", "penalized_functional_cox_bca")
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

file_suffix <- sprintf("B%d_seed%d_k%d", n_bootstrap, seed, main_basis_size)
bootstrap_file <- file.path(bootstrap_dir, paste0("bootstrap_estimates_", file_suffix, ".csv"))
if (!file.exists(bootstrap_file)) {
  stop("Run 18_bootstrap_penalized_functional_cox.R first: ", bootstrap_file)
}

# BCa interval with the same bias correction and jackknife acceleration as
# 08_bca_bootstrap_four_snp_mediation.R.
bca_interval <- function(values, point, jackknife_values) {
  proportion_less <- (sum(values < point) + 0.5 * sum(values == point)) / length(values)
  probability_bound <- 1 / (2 * length(values))
  proportion_less <- min(max(proportion_less, probability_bound), 1 - probability_bound)
  bias_correction <- stats::qnorm(proportion_less)

  jackknife_influence <- mean(jackknife_values) - jackknife_values
  acceleration_denominator <- 6 * sum(jackknife_influence^2)^(3 / 2)
  if (!is.finite(acceleration_denominator) || acceleration_denominator == 0) {
    stop("BCa acceleration is undefined.")
  }
  acceleration <- sum(jackknife_influence^3) / acceleration_denominator

  nominal_z <- stats::qnorm(c(0.025, 0.975))
  adjusted_probabilities <- stats::pnorm(
    bias_correction +
      (bias_correction + nominal_z) / (1 - acceleration * (bias_correction + nominal_z))
  )
  quantile_bound <- 1 / (length(values) + 1)
  adjusted_probabilities <- pmin(pmax(adjusted_probabilities, quantile_bound), 1 - quantile_bound)
  list(
    bias_correction = bias_correction,
    acceleration = acceleration,
    adjusted_probabilities = adjusted_probabilities,
    interval = stats::quantile(values, probs = adjusted_probabilities, names = FALSE, type = 7)
  )
}

run_penalized_bca <- function() {
bootstrap_data <- utils::read.csv(bootstrap_file, colClasses = c(cox_warning = "character"))
# read.csv turns the empty warning fields into NA.
bootstrap_data$cox_warning[is.na(bootstrap_data$cox_warning)] <- ""
valid_data <- bootstrap_data[!nzchar(bootstrap_data$cox_warning), , drop = FALSE]
point_data <- fit_penalized_replicate(seq_len(n), 0L)

message("Running leave-one-subject-out jackknife")
jackknife_results <- parallel::mclapply(
  seq_len(n),
  function(left_out) {
    tryCatch(
      fit_penalized_replicate(setdiff(seq_len(n), left_out), left_out),
      error = function(e) e
    )
  },
  mc.cores = n_cores
)
failed <- vapply(jackknife_results, inherits, logical(1), what = "error")
if (any(failed)) {
  stop("BCa acceleration requires all leave-one-subject-out fits to succeed.")
}
jackknife_data <- do.call(rbind, jackknife_results)
names(jackknife_data)[names(jackknife_data) == "replicate"] <- "left_out"
if (any(nzchar(jackknife_data$cox_warning))) {
  warning("Some jackknife fits produced warnings; see the jackknife estimates file.")
}

summary_rows <- list()
for (snp in snp_names) {
  for (effect in effect_names) {
    values <- valid_data[valid_data$SNP == snp, effect]
    point <- point_data[point_data$SNP == snp, effect]
    jackknife_values <- jackknife_data[jackknife_data$SNP == snp, effect]
    bca <- bca_interval(values, point, jackknife_values)
    percentile_quantiles <- stats::quantile(values, probs = c(0.025, 0.975), names = FALSE)
    basic_ci <- 2 * point - rev(percentile_quantiles)
    summary_rows[[length(summary_rows) + 1L]] <- data.frame(
      SNP = snp,
      effect = effect,
      point_log_HR = point,
      bootstrap_SE_log_HR = stats::sd(values),
      bias_correction_z0 = bca$bias_correction,
      acceleration = bca$acceleration,
      BCa_lower_probability = bca$adjusted_probabilities[1L],
      BCa_upper_probability = bca$adjusted_probabilities[2L],
      BCa_CI_2.5_log_HR = bca$interval[1L],
      BCa_CI_97.5_log_HR = bca$interval[2L],
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
bca_summary <- do.call(rbind, summary_rows)

utils::write.csv(
  jackknife_data,
  file.path(output_dir, sprintf("jackknife_estimates_n%d_k%d.csv", n, main_basis_size)),
  row.names = FALSE
)
utils::write.csv(
  bca_summary,
  file.path(output_dir, paste0("bca_summary_", file_suffix, ".csv")),
  row.names = FALSE
)

message("BCa completed. Results written to: ", output_dir)
print(
  bca_summary[
    bca_summary$effect %in% c("direct_log_HR", "indirect_log_HR", "total_log_HR"),
    c("SNP", "effect", "point_HR", "BCa_CI_2.5_HR", "BCa_CI_97.5_HR",
      "percentile_CI_2.5_HR", "percentile_CI_97.5_HR", "bias_correction_z0", "acceleration")
  ],
  row.names = FALSE, digits = 4
)
invisible(bca_summary)
}

if (sys.nframe() == 0L) {
  run_penalized_bca()
}
