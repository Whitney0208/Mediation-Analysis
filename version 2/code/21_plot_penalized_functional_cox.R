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
result_dir <- file.path(project_dir, "results", "penalized_functional_cox")
function_file <- file.path(result_dir, "penalized_functional_cox_functions_k20.csv")
if (!file.exists(function_file)) {
  stop("Run 16_fit_penalized_functional_cox.R first: ", function_file)
}

functions <- utils::read.csv(function_file)
snp_levels <- unique(functions$SNP)
column_labels <- c(
  alpha = "First layer: alpha_j(s)",
  beta = "Second layer: beta_j(s)",
  product = "Indirect-effect integrand: alpha_j(s) beta_j(s)"
)
plot_data <- do.call(rbind, lapply(1:2, function(coordinate) {
  alpha <- functions[[sprintf("alpha_coordinate_%d", coordinate)]]
  beta <- functions[[sprintf("beta_coordinate_%d", coordinate)]]
  do.call(rbind, lapply(names(column_labels), function(quantity) {
    data.frame(
      SNP = factor(functions$SNP, levels = snp_levels),
      grid = functions$grid,
      panel = sprintf("%s, coordinate %d", column_labels[[quantity]], coordinate),
      value = switch(quantity, alpha = alpha, beta = beta, product = alpha * beta)
    )
  }))
}))
plot_data$panel <- factor(
  plot_data$panel,
  levels = as.vector(t(outer(
    unname(column_labels), 1:2,
    function(label, coordinate) sprintf("%s, coordinate %d", label, coordinate)
  )))[c(1, 3, 5, 2, 4, 6)]
)

surface_color <- "#fcfcfb"
ink_primary <- "#0b0b0b"
ink_secondary <- "#52514e"
ink_muted <- "#898781"
snp_colors <- stats::setNames(c("#2a78d6", "#eb6834", "#1baf7a", "#eda100"), snp_levels)

figure <- ggplot2::ggplot(
  plot_data,
  ggplot2::aes(x = grid, y = value, colour = SNP)
) +
  ggplot2::geom_hline(yintercept = 0, colour = "#c3c2b7", linewidth = 0.3) +
  ggplot2::geom_line(linewidth = 0.7, lineend = "round", linejoin = "round") +
  ggplot2::facet_wrap(~panel, ncol = 3L, scales = "free_y") +
  ggplot2::scale_colour_manual(values = snp_colors, name = NULL) +
  ggplot2::scale_x_continuous(breaks = seq(0, 1, by = 0.25)) +
  ggplot2::guides(colour = ggplot2::guide_legend(override.aes = list(linewidth = 1.2))) +
  ggplot2::labs(
    title = "Penalized functional Cox model: coefficient functions by SNP model",
    subtitle = paste0(
      "Each SNP has its own first-layer alpha_j(s) and second-layer beta_j(s). ",
      "The area under the right-hand curves, summed over both\ncoordinates, ",
      "is the indirect log hazard ratio. Point estimates only; confidence bands are not shown."
    ),
    x = "Position along the normalized contour",
    y = NULL,
    caption = paste0(
      "beta_j(s): 20 cubic P-spline basis functions per coordinate, second-order ",
      "difference penalty, smoothing parameters selected by REML."
    )
  ) +
  ggplot2::theme_minimal(base_size = 11) +
  ggplot2::theme(
    plot.background = ggplot2::element_rect(fill = surface_color, colour = NA),
    plot.margin = ggplot2::margin(16, 20, 12, 16),
    plot.title = ggplot2::element_text(face = "bold", size = 15, colour = ink_primary),
    plot.title.position = "plot",
    plot.subtitle = ggplot2::element_text(
      size = 10.5, colour = ink_secondary, margin = ggplot2::margin(t = 2, b = 6)
    ),
    plot.caption = ggplot2::element_text(
      size = 8.5, colour = ink_muted, hjust = 0, margin = ggplot2::margin(t = 10)
    ),
    plot.caption.position = "plot",
    legend.position = "top",
    legend.justification = "left",
    legend.margin = ggplot2::margin(0, 0, 0, 0),
    legend.text = ggplot2::element_text(size = 10, colour = ink_secondary),
    legend.key.width = grid::unit(18, "pt"),
    panel.grid.major = ggplot2::element_line(colour = "#e1e0d9", linewidth = 0.25),
    panel.grid.minor = ggplot2::element_blank(),
    panel.spacing.x = grid::unit(16, "pt"),
    panel.spacing.y = grid::unit(12, "pt"),
    strip.text = ggplot2::element_text(
      size = 9, face = "bold", colour = ink_primary, hjust = 0
    ),
    axis.text = ggplot2::element_text(size = 8, colour = ink_muted),
    axis.title = ggplot2::element_text(size = 10, colour = ink_secondary),
    axis.title.x = ggplot2::element_text(margin = ggplot2::margin(t = 8))
  )

output_file <- file.path(result_dir, "penalized_functional_cox_functions_k20.png")
ggplot2::ggsave(
  output_file,
  figure,
  width = 13, height = 7.2, dpi = 200, bg = surface_color
)

message("Completed. Figure written to: ", output_file)
