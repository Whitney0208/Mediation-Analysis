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

output_dir <- file.path(project_dir, "results", "penalized_functional_cox_diagnostics")
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

# The fitted functional term sum_j integral M_j(s) beta_j(s) ds is summarized
# as one shape score per subject. Refitting the Cox model with that score as a
# covariate allows standard Schoenfeld-residual and interaction tests. These
# tests condition on the estimated beta_j(s).
shape_score <- function(cox_fit, shape_array) {
  rowSums(vapply(seq_len(2L), function(coordinate) {
    as.numeric(
      sweep(shape_array[, , coordinate], 2L, cox_fit$grid_weights, FUN = "*") %*%
        cox_fit$beta[, coordinate]
    )
  }, numeric(dim(shape_array)[1L])))
}

run_assumption_checks <- function() {
W_cox <- W[, setdiff(colnames(W), "intercept"), drop = FALSE]
ph_rows <- list()
interaction_rows <- list()
stratified_rows <- list()

for (snp in snp_names) {
  message("Checking ", snp)
  cox_fit <- fit_penalized_functional_cox(X[, snp], W, subject, M, grid, main_basis_size)
  cox_data <- data.frame(
    time = as.numeric(subject$time),
    event = as.integer(subject$event),
    SNP = as.numeric(X[, snp]),
    W_cox,
    shape_score = shape_score(cox_fit, M),
    check.names = FALSE
  )
  base_formula <- stats::as.formula(paste(
    "survival::Surv(time, event) ~",
    paste(c("SNP", colnames(W_cox), "shape_score"), collapse = " + ")
  ))
  base_fit <- survival::coxph(base_formula, data = cox_data, ties = "breslow")

  # Proportional hazards: scaled Schoenfeld residuals against KM-transformed time.
  zph <- survival::cox.zph(base_fit, transform = "km")
  ph_rows[[snp]] <- data.frame(
    SNP_model = snp,
    term = rownames(zph$table),
    chisq = zph$table[, "chisq"],
    df = zph$table[, "df"],
    p_value = zph$table[, "p"],
    row.names = NULL
  )

  # Exposure-mediator interaction: SNP x shape score, likelihood ratio test.
  interaction_fit <- survival::coxph(
    stats::update(base_formula, . ~ . + SNP:shape_score),
    data = cox_data,
    ties = "breslow"
  )
  lr_statistic <- 2 * (interaction_fit$loglik[2L] - base_fit$loglik[2L])
  interaction_summary <- summary(interaction_fit)$coefficients["SNP:shape_score", ]
  interaction_rows[[snp]] <- data.frame(
    SNP_model = snp,
    shape_score_coefficient = unname(stats::coef(base_fit)["shape_score"]),
    interaction_coefficient = unname(interaction_summary["coef"]),
    interaction_se = unname(interaction_summary["se(coef)"]),
    likelihood_ratio_chisq = lr_statistic,
    likelihood_ratio_p = stats::pchisq(lr_statistic, df = 1L, lower.tail = FALSE)
  )

  # APOE4 violated proportional hazards in the first run of these checks, so
  # refit with a separate baseline hazard for each APOE4 genotype (0/1/2).
  # The first layer still adjusts for APOE4 linearly.
  alpha <- fit_shape_layer(X[, snp], W, M, grid)$alpha
  stratified_fit <- fit_penalized_functional_cox(
    X[, snp], W[, colnames(W) != "APOE4", drop = FALSE], subject, M, grid,
    main_basis_size, strata = W[, "APOE4"]
  )
  main_effects <- penalized_effects(alpha, cox_fit$beta, cox_fit$grid_weights)
  stratified_effects <- penalized_effects(alpha, stratified_fit$beta, stratified_fit$grid_weights)
  stratified_rows[[snp]] <- data.frame(
    SNP_model = snp,
    model = c("main (APOE4 covariate)", "APOE4-stratified baseline hazard"),
    edf_beta_coordinate_1 = c(cox_fit$edf[1L], stratified_fit$edf[1L]),
    edf_beta_coordinate_2 = c(cox_fit$edf[2L], stratified_fit$edf[2L]),
    direct_HR = exp(c(cox_fit$tau, stratified_fit$tau)),
    indirect_HR = exp(c(main_effects[["indirect_log_HR"]], stratified_effects[["indirect_log_HR"]])),
    total_HR = exp(c(
      cox_fit$tau + main_effects[["indirect_log_HR"]],
      stratified_fit$tau + stratified_effects[["indirect_log_HR"]]
    ))
  )
}

ph_table <- do.call(rbind, ph_rows)
interaction_table <- do.call(rbind, interaction_rows)
stratified_table <- do.call(rbind, stratified_rows)

# Basis dimension: mgcv::k.check is not available for the cox.ph family, so the
# check uses the refits with 10, 20, and 30 basis functions from script 16.
basis_file <- file.path(project_dir, "results", "penalized_functional_cox", "penalized_functional_cox_effect_summary.csv")
basis_table <- utils::read.csv(basis_file)[, c(
  "SNP", "basis_size", "edf_beta_coordinate_1", "edf_beta_coordinate_2",
  "direct_HR", "indirect_HR", "total_HR"
)]

utils::write.csv(ph_table, file.path(output_dir, "proportional_hazards_tests.csv"), row.names = FALSE)
utils::write.csv(interaction_table, file.path(output_dir, "snp_shape_interaction_tests.csv"), row.names = FALSE)
utils::write.csv(basis_table, file.path(output_dir, "basis_dimension_check.csv"), row.names = FALSE)
utils::write.csv(stratified_table, file.path(output_dir, "apoe4_stratified_sensitivity.csv"), row.names = FALSE)

message("Completed. Results written to: ", output_dir)
print(ph_table, row.names = FALSE, digits = 3)
print(interaction_table, row.names = FALSE, digits = 3)
print(stratified_table, row.names = FALSE, digits = 4)
invisible(list(
  ph = ph_table, interaction = interaction_table,
  stratified = stratified_table, basis = basis_table
))
}

if (sys.nframe() == 0L) {
  run_assumption_checks()
}
