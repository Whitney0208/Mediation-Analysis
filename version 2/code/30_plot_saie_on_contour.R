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
model_script <- file.path(code_dir, "02_fit_four_snp_mediation_models.R")
source(model_script, local = FALSE)

saie_dir <- file.path(project_dir, "results", "penalized_functional_cox_saie")
saie_file <- file.path(saie_dir, "saie_summary_B10000_seed20260905_k20.csv")
if (!file.exists(saie_file)) {
  stop("Run 25_saie_penalized_functional_cox.R first: ", saie_file)
}
saie <- utils::read.csv(saie_file)

# Mean contour rebuilt from the sample-mean SRVF, as in 12_plot_fpc_shape.R:
# the contour is the running integral of q(s) |q(s)|, centered and rotated
# 90 degrees clockwise so that the long axis is horizontal.
p_grid <- length(grid)
mean_srvf <- cbind(colMeans(M[, , 1L]), colMeans(M[, , 2L]))
speed <- sqrt(rowSums(mean_srvf^2))
velocity <- mean_srvf * speed
step <- diff(grid)
contour_x <- c(0, cumsum(step * (velocity[-1L, 1L] + velocity[-p_grid, 1L]) / 2))
contour_y <- c(0, cumsum(step * (velocity[-1L, 2L] + velocity[-p_grid, 2L]) / 2))
contour <- data.frame(
  grid = grid,
  x = contour_y - mean(contour_y),
  y = -(contour_x - mean(contour_x))
)

# Attach contour coordinates by grid order; the SAIE file stores the same
# 100 grid points per SNP, so matching on floating-point grid values is avoided.
snp_levels <- unique(saie$SNP)
plot_data <- do.call(rbind, lapply(snp_levels, function(snp) {
  rows <- saie[saie$SNP == snp, ]
  rows <- rows[order(rows$grid), ]
  if (nrow(rows) != p_grid) {
    stop("Expected ", p_grid, " SAIE grid points for ", snp, ".")
  }
  cbind(rows, contour[, c("x", "y")])
}))
plot_data$SNP <- factor(plot_data$SNP, levels = snp_levels)
marker_positions <- c(0, 0.25, 0.5, 0.75)
markers <- contour[vapply(marker_positions, function(position) which.min(abs(grid - position)), integer(1L)), ]
markers$label <- c("s = 0", "s = 0.25", "s = 0.5", "s = 0.75")
# Fixed offsets keep labels off the contour: s = 0.25 lies on the lower
# boundary and s = 0.75 on the upper boundary of the same region.
offset <- 0.02
markers$label_x <- markers$x + c(-offset, 0, offset, 0)
markers$label_y <- markers$y + c(-offset, -offset, -offset, offset)
peak_rows <- do.call(rbind, lapply(split(plot_data, plot_data$SNP), function(rows) rows[which.max(abs(rows$saie)), ]))

surface_color <- "#fcfcfb"
ink_primary <- "#0b0b0b"
ink_secondary <- "#52514e"
ink_muted <- "#898781"
limit <- max(abs(plot_data$saie))

figure <- ggplot2::ggplot(plot_data, ggplot2::aes(x = x, y = y)) +
  ggplot2::geom_polygon(data = contour, fill = ink_secondary, alpha = 0.06, colour = NA) +
  ggplot2::geom_path(ggplot2::aes(colour = saie), linewidth = 2.2, lineend = "round", linejoin = "round") +
  ggplot2::geom_point(
    data = peak_rows, shape = 21, size = 3.4, stroke = 1.2,
    fill = NA, colour = ink_primary
  ) +
  ggplot2::geom_point(
    data = markers, ggplot2::aes(x = x, y = y), inherit.aes = FALSE,
    size = 1.6, colour = ink_muted
  ) +
  ggplot2::geom_text(
    data = markers, ggplot2::aes(x = label_x, y = label_y, label = label),
    inherit.aes = FALSE, size = 3, colour = ink_muted
  ) +
  ggplot2::facet_wrap(~SNP, ncol = 2L) +
  ggplot2::coord_equal() +
  ggplot2::scale_colour_gradient2(
    low = "#2a78d6", mid = "#c3c2b7", high = "#e34948", midpoint = 0,
    limits = c(-limit, limit), name = "SAIE(s)"
  ) +
  ggplot2::labs(
    title = "Spatial indirect effect on the mean corpus callosum contour",
    subtitle = paste0(
      "Colour: SAIE(s) point estimate (red raises, blue lowers the hazard of AD conversion). ",
      "Circle: position of the largest |SAIE|.\nNo position had a simultaneous or pointwise 95% ",
      "bootstrap interval excluding 0."
    ),
    caption = paste0(
      "Mean contour rebuilt from the sample-mean SRVF and rotated 90 degrees clockwise from the source ",
      "landmark orientation.\nGrey dots mark contour positions s; the anatomical labels of the two ends ",
      "should be confirmed against the source landmarks."
    )
  ) +
  ggplot2::theme_void(base_size = 11) +
  ggplot2::theme(
    plot.background = ggplot2::element_rect(fill = surface_color, colour = NA),
    plot.margin = ggplot2::margin(16, 20, 12, 16),
    plot.title = ggplot2::element_text(face = "bold", size = 15, colour = ink_primary),
    plot.title.position = "plot",
    plot.subtitle = ggplot2::element_text(size = 10, colour = ink_secondary, margin = ggplot2::margin(t = 4, b = 8)),
    plot.caption = ggplot2::element_text(size = 8.5, colour = ink_muted, hjust = 0, margin = ggplot2::margin(t = 10)),
    plot.caption.position = "plot",
    strip.text = ggplot2::element_text(size = 10, face = "bold", colour = ink_primary, hjust = 0, margin = ggplot2::margin(b = 4)),
    panel.spacing = grid::unit(18, "pt"),
    legend.title = ggplot2::element_text(size = 9, colour = ink_secondary),
    legend.text = ggplot2::element_text(size = 8, colour = ink_muted)
  )

output_file <- file.path(saie_dir, "saie_on_mean_contour.png")
ggplot2::ggsave(output_file, figure, width = 11, height = 7, dpi = 200, bg = surface_color)
message("Completed. Figure written to: ", output_file)
