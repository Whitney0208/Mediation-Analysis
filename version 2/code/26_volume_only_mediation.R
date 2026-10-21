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
bca_script <- file.path(code_dir, "19_bca_penalized_functional_cox.R")
for (required_script in c(model_script, bca_script)) {
  if (!file.exists(required_script)) {
    stop("Required script not found: ", required_script)
  }
}
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

output_dir <- file.path(project_dir, "results", "volume_only_mediation")
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

args <- commandArgs(trailingOnly = TRUE)
n_bootstrap <- if (length(args) >= 1L) as.integer(args[1L]) else 1000L
seed <- if (length(args) >= 2L) as.integer(args[2L]) else 20260905L
n_cores <- if (length(args) >= 3L) as.integer(args[3L]) else 1L
if (is.na(n_bootstrap) || n_bootstrap < 1L || is.na(seed) || is.na(n_cores) || n_cores < 1L) {
  stop("Replicates, seed, and cores must be positive integers.")
}

# Mediators: the five FreeSurfer CC subregion volumes, each divided by ICV and
# standardized in the analysis sample. This mirrors the volumetric mediators of
# S-GMAS, with a Cox second layer because the outcome is time to AD conversion.
volume_file <- file.path(project_dir, "data", "SPLS-Coxdata", "CCinfo.csv")
cc_regions <- c(
  anterior = "CorpusCallosumAnterior",
  mid_anterior = "CorpusCallosumMidAnterior",
  central = "CorpusCallosumCentral",
  mid_posterior = "CorpusCallosumMidPosterior",
  posterior = "CorpusCallosumPosterior"
)
volume_table <- utils::read.csv(volume_file)
volume_rows <- match(as.integer(subject$RID), as.integer(volume_table$RID))
if (anyNA(volume_rows) || anyNA(volume_table[volume_rows, c(cc_regions, "ICV")])) {
  stop("Corpus callosum subregion volumes are missing for some analysis subjects.")
}
volume_table <- volume_table[volume_rows, , drop = FALSE]
region_ratio <- as.matrix(volume_table[, cc_regions]) / volume_table$ICV
region_mediators <- scale(region_ratio)
colnames(region_mediators) <- names(cc_regions)
W_cox <- W[, setdiff(colnames(W), "intercept"), drop = FALSE]

# First layer: V_r = a_r X + c_r' W + e_r for each region.
# Joint model:      h(t) = h0(t) exp[tau X + kappa' W + sum_r b_r V_r],
#                   indirect = sum_r a_r b_r.
# Single-region:    one model per region with only V_r, indirect = a_r b_r.
fit_volume_replicate <- function(indices, replicate_id) {
  b_volume <- region_mediators[indices, , drop = FALSE]
  b_W <- W[indices, , drop = FALSE]
  b_W_cox <- W_cox[indices, , drop = FALSE]
  b_time <- as.numeric(subject$time[indices])
  b_event <- as.integer(subject$event[indices])
  rows <- lapply(snp_names, function(snp) {
    x <- as.numeric(X[indices, snp])
    snp_to_volume <- qr.coef(qr(cbind(SNP = x, b_W)), b_volume)["SNP", ]
    fit_warnings <- character()
    fit_cox <- function(mediators) {
      withCallingHandlers(
        survival::coxph(
          survival::Surv(b_time, b_event) ~ x + b_W_cox + mediators,
          ties = "breslow"
        ),
        warning = function(w) {
          fit_warnings <<- c(fit_warnings, conditionMessage(w))
          invokeRestart("muffleWarning")
        }
      )
    }
    joint_fit <- fit_cox(b_volume)
    joint_b <- stats::coef(joint_fit)[paste0("mediators", colnames(b_volume))]
    joint_direct <- unname(stats::coef(joint_fit)["x"])
    joint_indirect <- sum(snp_to_volume * joint_b)
    single_indirect <- vapply(colnames(b_volume), function(region) {
      single_fit <- fit_cox(b_volume[, region, drop = FALSE])
      single_b <- stats::coef(single_fit)[grep("^mediators", names(stats::coef(single_fit)))]
      unname(snp_to_volume[region] * single_b)
    }, numeric(1L))
    data.frame(
      replicate = replicate_id,
      SNP = snp,
      direct_log_HR = joint_direct,
      indirect_joint_log_HR = joint_indirect,
      total_log_HR = joint_direct + joint_indirect,
      as.list(stats::setNames(single_indirect, paste0("indirect_", names(single_indirect), "_log_HR"))),
      cox_warning = paste(unique(fit_warnings), collapse = " | ")
    )
  })
  do.call(rbind, rows)
}

effect_names <- c(
  "direct_log_HR", "indirect_joint_log_HR", "total_log_HR",
  paste0("indirect_", names(cc_regions), "_log_HR")
)

run_volume_only <- function() {
point_data <- fit_volume_replicate(seq_len(n), 0L)

# Indices are drawn exactly as in 18_bootstrap_penalized_functional_cox.R.
set.seed(seed)
replicate_indices <- lapply(seq_len(n_bootstrap), function(replicate_id) {
  sample.int(n, size = n, replace = TRUE)
})
results <- parallel::mclapply(
  seq_len(n_bootstrap),
  function(replicate_id) {
    tryCatch(fit_volume_replicate(replicate_indices[[replicate_id]], replicate_id), error = function(e) e)
  },
  mc.cores = n_cores
)
failed <- vapply(results, function(result) !is.data.frame(result), logical(1))
bootstrap_data <- do.call(rbind, results[!failed])
valid_data <- bootstrap_data[!nzchar(bootstrap_data$cox_warning), , drop = FALSE]

jackknife_data <- do.call(rbind, parallel::mclapply(
  seq_len(n),
  function(left_out) fit_volume_replicate(setdiff(seq_len(n), left_out), left_out),
  mc.cores = n_cores
))

summary_rows <- list()
for (snp in snp_names) {
  for (effect in effect_names) {
    values <- valid_data[valid_data$SNP == snp, effect]
    point <- point_data[point_data$SNP == snp, effect]
    bca <- bca_interval(values, point, jackknife_data[jackknife_data$SNP == snp, effect])
    percentile_quantiles <- stats::quantile(values, probs = c(0.025, 0.975), names = FALSE)
    summary_rows[[length(summary_rows) + 1L]] <- data.frame(
      SNP = snp,
      effect = effect,
      point_log_HR = point,
      bootstrap_SE_log_HR = stats::sd(values),
      bias_correction_z0 = bca$bias_correction,
      acceleration = bca$acceleration,
      point_HR = exp(point),
      BCa_CI_2.5_HR = exp(bca$interval[1L]),
      BCa_CI_97.5_HR = exp(bca$interval[2L]),
      percentile_CI_2.5_HR = exp(percentile_quantiles[1L]),
      percentile_CI_97.5_HR = exp(percentile_quantiles[2L]),
      successful_replicates = length(values)
    )
  }
}
bootstrap_summary <- do.call(rbind, summary_rows)

file_suffix <- sprintf("B%d_seed%d", n_bootstrap, seed)
utils::write.csv(point_data, file.path(output_dir, "volume_only_point_estimates.csv"), row.names = FALSE)
utils::write.csv(bootstrap_summary, file.path(output_dir, paste0("bootstrap_summary_", file_suffix, ".csv")), row.names = FALSE)
utils::write.csv(bootstrap_data, file.path(output_dir, paste0("bootstrap_estimates_", file_suffix, ".csv")), row.names = FALSE)
utils::write.csv(jackknife_data, file.path(output_dir, sprintf("jackknife_estimates_n%d.csv", n)), row.names = FALSE)
message(
  "Volume-only mediation completed: ", sum(failed), " failed replicates, ",
  sum(nzchar(bootstrap_data$cox_warning)), " warning-flagged fits. Results written to: ", output_dir
)
print(
  bootstrap_summary[, c("SNP", "effect", "point_HR", "BCa_CI_2.5_HR", "BCa_CI_97.5_HR", "successful_replicates")],
  row.names = FALSE, digits = 4
)
invisible(bootstrap_summary)
}

if (sys.nframe() == 0L) {
  run_volume_only()
}
