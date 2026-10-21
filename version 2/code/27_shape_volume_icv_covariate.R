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
bootstrap_script <- file.path(code_dir, "23_bootstrap_shape_volume_mediation.R")
if (!file.exists(bootstrap_script)) {
  stop("Shape and volume bootstrap script not found: ", bootstrap_script)
}
# Also sources 22, 16, and 02, and reads the same arguments: replicates, seed, cores.
source(bootstrap_script, local = FALSE)

# Sensitivity analysis for the ICV handling of the shape + volume model.
# Main version: mediator = total CC volume / ICV.
# This version: mediator = total CC volume (standardized), and standardized ICV
# is added to the adjustment set of every equation (both layers).
raw_volume <- rowSums(volume_table[, cc_regions])
volume_mediator <- cbind(cc_volume_z = (raw_volume - mean(raw_volume)) / stats::sd(raw_volume))
icv_z <- (volume_table$ICV - mean(volume_table$ICV)) / stats::sd(volume_table$ICV)
W <- cbind(W, ICV_z = icv_z)

output_dir <- file.path(project_dir, "results", "shape_volume_icv_covariate_bootstrap")
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

if (sys.nframe() == 0L) {
  run_shape_volume_bootstrap()
}
