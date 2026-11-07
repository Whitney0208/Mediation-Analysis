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
source(penalized_script, local = FALSE)
if (!requireNamespace("parallel", quietly = TRUE)) {
  stop("Package 'parallel' is required.")
}

args <- commandArgs(trailingOnly = TRUE)
permutations <- if (length(args) >= 1L) as.integer(args[1L]) else 10000L
seed <- if (length(args) >= 2L) as.integer(args[2L]) else 20261102L
cores <- if (length(args) >= 3L) as.integer(args[3L]) else 6L

output_dir <- file.path(project_dir, "results", "penalized_functional_cox_path_tests")
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

# Global tests of the two paths of the main model, for each SNP.
#   SNP -> shape:    H0: alpha_1(s) = alpha_2(s) = 0 for all s, given W.
#   shape -> hazard: H0: beta_1(s) = beta_2(s) = 0 for all s, given SNP and W.
# Both are Freedman-Lane permutation tests. Subjects are permuted as whole
# curves, so dependence along the contour and between coordinates is kept.
# The same permutations are used for every SNP.
set.seed(seed)
permutation_table <- replicate(permutations, sample.int(n))

weights_flat <- rep(trapezoid_weights(grid), times = 2L)
shape_matrix <- cbind(M[, , 1L], M[, , 2L])

pointwise_t <- function(x_residual, shape_residual, residual_df) {
  x_sum <- sum(x_residual^2)
  slope <- drop(crossprod(x_residual, shape_residual)) / x_sum
  rss <- colSums(shape_residual^2) - slope^2 * x_sum
  slope / sqrt(rss / residual_df / x_sum)
}

# Integrated squared t over both coordinates, and its maximum counterpart.
shape_path_statistics <- function(t_values) {
  c(integrated = sum(weights_flat * t_values^2), maximum = max(abs(t_values)))
}

test_shape_path <- function(x) {
  W_qr <- qr(W)
  x_residual <- qr.resid(W_qr, x)
  shape_residual <- qr.resid(W_qr, shape_matrix)
  residual_df <- n - W_qr$rank - 1L
  t_values <- pointwise_t(x_residual, shape_residual, residual_df)
  observed <- shape_path_statistics(t_values)
  permuted <- apply(permutation_table, 2L, function(order) {
    permuted_residual <- qr.resid(W_qr, shape_residual[order, , drop = FALSE])
    shape_path_statistics(pointwise_t(x_residual, permuted_residual, residual_df))
  })
  largest <- which.max(abs(t_values))
  p_grid <- length(grid)
  data.frame(
    shape_integrated_t2 = observed[["integrated"]],
    shape_p_integrated = (1 + sum(permuted["integrated", ] >= observed[["integrated"]])) /
      (permutations + 1),
    shape_max_abs_t = observed[["maximum"]],
    shape_p_maximum = (1 + sum(permuted["maximum", ] >= observed[["maximum"]])) /
      (permutations + 1),
    shape_max_coordinate = if (largest <= p_grid) 1L else 2L,
    shape_max_grid = grid[(largest - 1L) %% p_grid + 1L],
    points_abs_t_above_1.96 = sum(abs(t_values) > 1.96)
  )
}

partial_log_likelihood <- function(linear_predictor) {
  survival::coxph(
    survival::Surv(subject$time, subject$event) ~ offset(linear_predictor),
    ties = "breslow"
  )$loglik[1L]
}

# Likelihood ratio of the REML fit against the Cox model without shape. Its null
# distribution comes from the permutations, not from a chi-squared limit.
hazard_path_statistic <- function(x, shape_array, null_log_likelihood) {
  fit <- fit_penalized_functional_cox(x, W, subject, shape_array, grid)
  2 * (partial_log_likelihood(as.numeric(stats::predict(fit$fit, type = "link"))) -
    null_log_likelihood)
}

test_hazard_path <- function(x) {
  W_cox <- W[, setdiff(colnames(W), "intercept"), drop = FALSE]
  null_fit <- survival::coxph(
    survival::Surv(subject$time, subject$event) ~ x + W_cox,
    ties = "breslow"
  )
  null_log_likelihood <- null_fit$loglik[2L]
  observed_fit <- fit_penalized_functional_cox(x, W, subject, M, grid)
  observed <- hazard_path_statistic(x, M, null_log_likelihood)
  approximate <- summary(observed_fit$fit)$s.table

  design_qr <- qr(cbind(SNP = as.numeric(x), W))
  fitted_shape <- qr.fitted(design_qr, shape_matrix)
  residual_shape <- qr.resid(design_qr, shape_matrix)
  chunks <- split(seq_len(permutations), cut(seq_len(permutations), cores, labels = FALSE))
  permuted <- unlist(parallel::mclapply(chunks, function(replicates) {
    vapply(replicates, function(replicate_id) {
      shape_star <- fitted_shape + residual_shape[permutation_table[, replicate_id], , drop = FALSE]
      tryCatch(
        suppressWarnings(hazard_path_statistic(
          x, array(shape_star, dim = dim(M)), null_log_likelihood
        )),
        error = function(e) NA_real_
      )
    }, numeric(1L))
  }, mc.cores = cores, mc.preschedule = FALSE), use.names = FALSE)
  completed <- sum(is.finite(permuted))

  data.frame(
    hazard_lr = observed,
    hazard_p_permutation = (1 + sum(permuted >= observed, na.rm = TRUE)) / (completed + 1),
    hazard_permutations_completed = completed,
    edf_beta_coordinate_1 = approximate[1L, "edf"],
    edf_beta_coordinate_2 = approximate[2L, "edf"],
    p_approximate_coordinate_1 = approximate[1L, "p-value"],
    p_approximate_coordinate_2 = approximate[2L, "p-value"]
  )
}

rows <- lapply(snp_names, function(snp) {
  message("Testing ", snp)
  x <- X[, snp]
  cbind(SNP = snp, test_shape_path(x), test_hazard_path(x))
})
path_tests <- do.call(rbind, rows)
path_tests$shape_p_integrated_holm <- stats::p.adjust(path_tests$shape_p_integrated, "holm")
path_tests$shape_p_maximum_holm <- stats::p.adjust(path_tests$shape_p_maximum, "holm")

output_file <- file.path(
  output_dir,
  sprintf("path_tests_B%d_seed%d_k%d.csv", permutations, seed, main_basis_size)
)
utils::write.csv(path_tests, output_file, row.names = FALSE)
message("Completed. Results written to: ", output_file)
print(path_tests, row.names = FALSE, digits = 4)
