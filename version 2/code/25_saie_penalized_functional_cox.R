#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    stop("Package 'ggplot2' is required.")
  }
})

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

output_dir <- file.path(project_dir, "results", "penalized_functional_cox_saie")
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

args <- commandArgs(trailingOnly = TRUE)
n_bootstrap <- if (length(args) >= 1L) as.integer(args[1L]) else 1000L
seed <- if (length(args) >= 2L) as.integer(args[2L]) else 20260905L
n_cores <- if (length(args) >= 3L) as.integer(args[3L]) else 1L
if (is.na(n_bootstrap) || n_bootstrap < 1L || is.na(seed) || is.na(n_cores) || n_cores < 1L) {
  stop("Replicates, seed, and cores must be positive integers.")
}

# Spatial indirect effect: SAIE(s) = sum_j alpha_j(s) beta_j(s). Its trapezoid
# integral over s equals the indirect log hazard ratio.
fit_saie_replicate <- function(indices) {
  b_subject <- subject[indices, , drop = FALSE]
  b_X <- X[indices, , drop = FALSE]
  b_M <- M[indices, , , drop = FALSE]
  b_W <- W[indices, , drop = FALSE]
  curves <- vapply(snp_names, function(snp) {
    alpha <- fit_shape_layer(b_X[, snp], b_W, b_M, grid)$alpha
    fit_warnings <- character()
    cox_fit <- withCallingHandlers(
      fit_penalized_functional_cox(b_X[, snp], b_W, b_subject, b_M, grid, main_basis_size),
      warning = function(w) {
        fit_warnings <<- c(fit_warnings, conditionMessage(w))
        invokeRestart("muffleWarning")
      }
    )
    if (length(fit_warnings) > 0L) {
      return(rep(NA_real_, length(grid)))
    }
    rowSums(alpha * cox_fit$beta)
  }, numeric(length(grid)))
  curves
}

run_saie <- function() {
weights_grid <- trapezoid_weights(grid)
point_curves <- fit_saie_replicate(seq_len(n))

# Indices are drawn exactly as in 18_bootstrap_penalized_functional_cox.R, so
# replicate b uses the same subjects as replicate b of the main pairs bootstrap.
set.seed(seed)
replicate_indices <- lapply(seq_len(n_bootstrap), function(replicate_id) {
  sample.int(n, size = n, replace = TRUE)
})

checkpoint_file <- file.path(output_dir, sprintf("saie_checkpoint_B%d_seed%d.rds", n_bootstrap, seed))
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
      tryCatch(fit_saie_replicate(replicate_indices[[replicate_id]]), error = function(e) e)
    },
    mc.cores = n_cores
  )
  names(chunk_results) <- as.character(chunk)
  results[names(chunk_results)] <- chunk_results
  saveRDS(results, checkpoint_file)
  message("Completed SAIE bootstrap replicate ", max(chunk), "/", n_bootstrap)
}
failed <- vapply(results, function(result) !is.matrix(result), logical(1))
curves <- simplify2array(results[!failed])  # grid x SNP x replicate

# Consistency check: integrals of the replicate curves must reproduce the
# indirect effects stored by script 18 when the same replicates are available.
pairs_file <- file.path(
  project_dir, "results", "penalized_functional_cox_bootstrap",
  sprintf("bootstrap_estimates_B%d_seed%d_k%d.csv", n_bootstrap, seed, main_basis_size)
)
consistency <- NA_real_
if (file.exists(pairs_file)) {
  pairs_data <- utils::read.csv(pairs_file)
  check_ids <- intersect(head(as.integer(names(results)[!failed]), 50L), pairs_data$replicate)
  differences <- unlist(lapply(check_ids, function(replicate_id) {
    integral <- colSums(weights_grid * results[[as.character(replicate_id)]])
    stored <- pairs_data[pairs_data$replicate == replicate_id, ]
    integral[stored$SNP] - stored$indirect_log_HR
  }))
  consistency <- max(abs(differences), na.rm = TRUE)
  if (consistency > 1e-6) {
    stop("SAIE integrals do not reproduce the stored indirect effects (max difference ", consistency, ").")
  }
}

summary_rows <- list()
band_rows <- list()
for (snp_index in seq_along(snp_names)) {
  snp <- snp_names[snp_index]
  point <- point_curves[, snp_index]
  replicate_curves <- curves[, snp_index, ]
  replicate_curves <- replicate_curves[, colSums(is.na(replicate_curves)) == 0L, drop = FALSE]
  pointwise_sd <- apply(replicate_curves, 1L, stats::sd)
  pointwise_quantiles <- apply(replicate_curves, 1L, stats::quantile, probs = c(0.025, 0.975))
  # Simultaneous band: 95th percentile of the maximum standardized deviation.
  max_deviation <- apply(abs(replicate_curves - point) / pointwise_sd, 2L, max)
  critical_value <- stats::quantile(max_deviation, 0.95, names = FALSE)
  band_lower <- point - critical_value * pointwise_sd
  band_upper <- point + critical_value * pointwise_sd
  summary_rows[[snp]] <- data.frame(
    SNP = snp,
    grid = grid,
    saie = point,
    bootstrap_sd = pointwise_sd,
    pointwise_lower = pointwise_quantiles[1L, ],
    pointwise_upper = pointwise_quantiles[2L, ],
    simultaneous_lower = band_lower,
    simultaneous_upper = band_upper,
    simultaneous_excludes_zero = band_lower > 0 | band_upper < 0,
    pointwise_excludes_zero = pointwise_quantiles[1L, ] > 0 | pointwise_quantiles[2L, ] < 0
  )
  band_rows[[snp]] <- data.frame(
    SNP = snp,
    indirect_log_HR = sum(weights_grid * point),
    simultaneous_critical_value = critical_value,
    grid_points_band_excludes_zero = sum(band_lower > 0 | band_upper < 0),
    grid_points_pointwise_excludes_zero = sum(pointwise_quantiles[1L, ] > 0 | pointwise_quantiles[2L, ] < 0),
    valid_replicates = ncol(replicate_curves)
  )
}
saie_table <- do.call(rbind, summary_rows)
band_table <- do.call(rbind, band_rows)
band_table$failed_replicates <- sum(failed)
band_table$max_integral_difference_vs_script18 <- consistency

file_suffix <- sprintf("B%d_seed%d_k%d", n_bootstrap, seed, main_basis_size)
utils::write.csv(saie_table, file.path(output_dir, paste0("saie_summary_", file_suffix, ".csv")), row.names = FALSE)
utils::write.csv(band_table, file.path(output_dir, paste0("saie_band_summary_", file_suffix, ".csv")), row.names = FALSE)
if (file.exists(checkpoint_file)) {
  unlink(checkpoint_file)
}

# Figure: SAIE(s) with pointwise and simultaneous 95% bands.
surface_color <- "#fcfcfb"
ink_primary <- "#0b0b0b"
ink_secondary <- "#52514e"
ink_muted <- "#898781"
saie_table$SNP <- factor(saie_table$SNP, levels = snp_names)
highlight <- saie_table[saie_table$simultaneous_excludes_zero, ]
figure <- ggplot2::ggplot(saie_table, ggplot2::aes(x = grid)) +
  ggplot2::geom_ribbon(
    ggplot2::aes(ymin = simultaneous_lower, ymax = simultaneous_upper),
    fill = "#2a78d6", alpha = 0.12
  ) +
  ggplot2::geom_ribbon(
    ggplot2::aes(ymin = pointwise_lower, ymax = pointwise_upper),
    fill = "#2a78d6", alpha = 0.22
  ) +
  ggplot2::geom_hline(yintercept = 0, colour = "#c3c2b7", linewidth = 0.3) +
  ggplot2::geom_line(ggplot2::aes(y = saie), colour = "#2a78d6", linewidth = 0.8) +
  ggplot2::geom_point(
    data = highlight, ggplot2::aes(y = saie),
    colour = "#e34948", size = 1.2
  ) +
  ggplot2::facet_wrap(~SNP, ncol = 2L) +
  ggplot2::scale_x_continuous(breaks = seq(0, 1, by = 0.25)) +
  ggplot2::labs(
    title = "Spatial indirect effect along the corpus callosum contour",
    subtitle = paste0(
      "SAIE(s) = alpha_1(s) beta_1(s) + alpha_2(s) beta_2(s); its integral over s is the indirect log HR.\n",
      "Dark band: pointwise 95% interval. Light band: simultaneous 95% band. ",
      "Red points: positions where the simultaneous band excludes 0."
    ),
    x = "Position along the normalized contour (s)",
    y = "SAIE(s), log hazard ratio per unit s",
    caption = sprintf(
      "Pairs bootstrap, %s replicates, smoothing parameters reselected by REML in every replicate.",
      format(n_bootstrap, big.mark = ",")
    )
  ) +
  ggplot2::theme_minimal(base_size = 11) +
  ggplot2::theme(
    plot.background = ggplot2::element_rect(fill = surface_color, colour = NA),
    plot.margin = ggplot2::margin(16, 20, 12, 16),
    plot.title = ggplot2::element_text(face = "bold", size = 15, colour = ink_primary),
    plot.title.position = "plot",
    plot.subtitle = ggplot2::element_text(size = 10, colour = ink_secondary, margin = ggplot2::margin(t = 2, b = 8)),
    plot.caption = ggplot2::element_text(size = 8.5, colour = ink_muted, hjust = 0, margin = ggplot2::margin(t = 10)),
    plot.caption.position = "plot",
    panel.grid.major = ggplot2::element_line(colour = "#e1e0d9", linewidth = 0.25),
    panel.grid.minor = ggplot2::element_blank(),
    panel.spacing.x = grid::unit(16, "pt"),
    panel.spacing.y = grid::unit(12, "pt"),
    strip.text = ggplot2::element_text(size = 10, face = "bold", colour = ink_primary, hjust = 0),
    axis.text = ggplot2::element_text(size = 8, colour = ink_muted),
    axis.title = ggplot2::element_text(size = 10, colour = ink_secondary)
  )
ggplot2::ggsave(
  file.path(output_dir, paste0("saie_", file_suffix, ".png")),
  figure, width = 11, height = 7.5, dpi = 200, bg = surface_color
)

message("SAIE completed. Results written to: ", output_dir)
print(band_table, row.names = FALSE, digits = 4)
invisible(list(saie = saie_table, bands = band_table))
}

if (sys.nframe() == 0L) {
  run_saie()
}
