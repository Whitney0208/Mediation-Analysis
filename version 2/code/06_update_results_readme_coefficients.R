#!/usr/bin/env Rscript

script_args <- commandArgs(trailingOnly = FALSE)
script_file <- grep("^--file=", script_args, value = TRUE)
if (length(script_file) != 1L) {
  stop("Run this script with Rscript.")
}
project_dir <- dirname(dirname(normalizePath(sub("^--file=", "", script_file))))
results_dir <- file.path(project_dir, "results")
coefficient_dir <- file.path(results_dir, "four_snp_point_estimates")
readme_path <- file.path(results_dir, "README.md")
effect_data <- read.csv(file.path(coefficient_dir, "four_snp_effect_summary.csv"))
readme <- readLines(readme_path, warn = FALSE)
begin_marker <- "<!-- BEGIN FOUR-MODEL COEFFICIENTS -->"
end_marker <- "<!-- END FOUR-MODEL COEFFICIENTS -->"
begin <- which(readme == begin_marker)
end <- which(readme == end_marker)
if (length(begin) != 1L || length(end) != 1L || begin >= end) {
  stop("The results README must contain one ordered pair of coefficient markers.")
}

number <- function(x) sprintf("%.6f", x)
predictors <- c(
  "SNP", "intercept", "male", "age_z", "handedness", "APOE4",
  "PC1", "PC2", "education_years"
)
scalar_terms <- c("SNP", predictors[-c(1L, 2L)])
fpc_terms <- sprintf("FPC%02d", seq_len(18L))
table_row <- function(cells) paste0("| ", paste(cells, collapse = " | "), " |")

render_model <- function(snp) {
  shape_path <- file.path(coefficient_dir, paste0(snp, "_shape_coefficients.csv"))
  cox_path <- file.path(coefficient_dir, paste0(snp, "_cox_coefficients.csv"))
  function_path <- file.path(coefficient_dir, paste0(snp, "_functional_parameters.csv"))
  shape <- read.csv(shape_path)
  cox <- read.csv(cox_path)
  functional <- read.csv(function_path)
  if (
    nrow(shape) != 900L || nrow(cox) != 26L || nrow(functional) != 100L ||
    !identical(unique(shape$predictor), predictors) ||
    !identical(as.character(cox$term), c(scalar_terms, fpc_terms)) ||
    !identical(as.character(unique(shape$SNP_model)), snp) ||
    !identical(as.character(unique(cox$SNP_model)), snp)
  ) {
    stop("Coefficient dimensions or term order changed for ", snp)
  }

  scalar <- cox[seq_along(scalar_terms), , drop = FALSE]
  fpc <- cox[-seq_along(scalar_terms), , drop = FALSE]
  scalar_table <- c(
    table_row(c("Term", "Log-hazard coefficient", "HR per unit")),
    "|---|---:|---:|",
    vapply(seq_len(nrow(scalar)), function(i) {
      table_row(c(
        paste0("`", scalar$term[i], "`"),
        number(scalar$coefficient[i]), number(scalar$HR_per_unit[i])
      ))
    }, character(1))
  )
  fpc_table <- c(
    table_row(c("Shape FPC", "Log-hazard coefficient")),
    "|---|---:|",
    vapply(seq_len(nrow(fpc)), function(i) {
      table_row(c(paste0("`", fpc$term[i], "`"), number(fpc$coefficient[i])))
    }, character(1))
  )

  coordinate_table <- function(coordinate) {
    column <- paste0("coordinate_", coordinate)
    beta_column <- paste0("beta_coordinate_", coordinate)
    values <- vapply(predictors, function(predictor) {
      rows <- shape[shape$predictor == predictor, , drop = FALSE]
      if (nrow(rows) != 100L || !identical(rows$grid, functional$grid)) {
        stop("Shape grid mismatch for ", snp, ": ", predictor)
      }
      rows[[column]]
    }, numeric(100L))
    if (!isTRUE(all.equal(
      values[, "SNP"], functional[[paste0("alpha_coordinate_", coordinate)]],
      tolerance = 1e-12
    ))) {
      stop("SNP coefficient mismatch for ", snp, " coordinate ", coordinate)
    }
    c(
      table_row(c("Grid s", paste0("`", predictors, "`"), "M: beta(s)")),
      paste0("|---:", paste(rep("|---:", length(predictors) + 1L), collapse = ""), "|"),
      vapply(seq_len(100L), function(i) {
        table_row(c(
          number(functional$grid[i]), number(values[i, ]),
          number(functional[[beta_column]][i])
        ))
      }, character(1))
    )
  }

  c(
    paste0("### ", snp),
    "",
    "Second-layer Cox coefficients for the SNP and confounders:",
    "",
    scalar_table,
    "",
    "Second-layer coefficients of the 18 joint shape FPC scores:",
    "",
    fpc_table,
    "",
    "First-layer coefficients and reconstructed shape-mediator coefficient functions:",
    "",
    "<details>",
    "<summary>Coordinate 1: all 100 contour positions</summary>",
    "",
    coordinate_table(1L),
    "",
    "</details>",
    "",
    "<details>",
    "<summary>Coordinate 2: all 100 contour positions</summary>",
    "",
    coordinate_table(2L),
    "",
    "</details>",
    "",
    paste0(
      "Full precision: [shape coefficients](four_snp_point_estimates/", snp,
      "_shape_coefficients.csv), [Cox coefficients](four_snp_point_estimates/", snp,
      "_cox_coefficients.csv), [shape beta and indirect integrands](four_snp_point_estimates/", snp,
      "_functional_parameters.csv)."
    ),
    ""
  )
}

models <- unlist(lapply(as.character(effect_data$SNP), render_model), use.names = FALSE)
updated <- c(readme[seq_len(begin)], "", models, readme[end:length(readme)])
writeLines(updated, readme_path, useBytes = TRUE)
cat("Updated coefficients for", nrow(effect_data), "models in", readme_path, "\n")
