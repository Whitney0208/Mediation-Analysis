#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  if (!requireNamespace("survival", quietly = TRUE)) {
    stop("Package 'survival' is required.")
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
input_file <- file.path(project_dir, "data", "mediation_model_matched.rds")
output_dir <- file.path(project_dir, "results", "four_snp_point_estimates")
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

if (!file.exists(input_file)) {
  stop("Input file not found: ", input_file)
}

dat <- readRDS(input_file)
required_objects <- c("subject", "X", "M", "W_pc2", "snp_names", "shape_grid")
missing_objects <- setdiff(required_objects, names(dat))
if (length(missing_objects) > 0L) {
  stop("Missing objects in input: ", paste(missing_objects, collapse = ", "))
}

subject <- as.data.frame(dat$subject)
X <- as.matrix(dat$X)
M <- dat$M
W <- as.matrix(dat$W_pc2)
grid <- as.numeric(dat$shape_grid)
snp_names <- as.character(dat$snp_names)

n <- nrow(X)
if (!identical(nrow(subject), n) || !identical(dim(M)[1], n) || !identical(nrow(W), n)) {
  stop("Subject counts do not agree across subject, X, M, and W_pc2.")
}
if (length(dim(M)) != 3L || dim(M)[3] != 2L) {
  stop("M must have dimensions subject x grid x 2 coordinates.")
}
if (dim(M)[2] != length(grid)) {
  stop("The M grid dimension does not agree with shape_grid.")
}
if (!all(c("RID", "time", "event") %in% names(subject))) {
  stop("subject must contain RID, time, and event.")
}
if (!all(snp_names %in% colnames(X))) {
  stop("Not all snp_names are present in X.")
}
if (!"intercept" %in% colnames(W)) {
  stop("W_pc2 must contain an intercept column.")
}

complete <- complete.cases(subject[, c("RID", "time", "event")], X, W) &
  apply(M, 1L, function(z) all(is.finite(z)))
if (!all(complete)) {
  message("Removing ", sum(!complete), " subjects with incomplete model data.")
  subject <- subject[complete, , drop = FALSE]
  X <- X[complete, , drop = FALSE]
  M <- M[complete, , , drop = FALSE]
  W <- W[complete, , drop = FALSE]
  n <- nrow(X)
}

if (!all(subject$event %in% c(0, 1))) {
  stop("event must be coded 0/1.")
}
if (any(subject$time <= 0)) {
  stop("All survival times must be positive.")
}

trapezoid_weights <- function(x) {
  if (length(x) < 2L || any(diff(x) <= 0)) {
    stop("shape_grid must be strictly increasing.")
  }
  dx <- diff(x)
  c(dx[1] / 2, (dx[-1] + dx[-length(dx)]) / 2, dx[length(dx)] / 2)
}

smooth_shape_coefficients <- function(raw_coef, shape_grid, p_grid) {
  n_coef <- nrow(raw_coef)
  smooth_coef <- array(
    NA_real_,
    dim = c(n_coef, p_grid, 2L),
    dimnames = list(rownames(raw_coef), NULL, c("coordinate_1", "coordinate_2"))
  )
  smoothing <- vector("list", n_coef * 2L)
  names(smoothing) <- as.vector(outer(
    rownames(raw_coef), c("coordinate_1", "coordinate_2"), paste, sep = ":"
  ))

  counter <- 1L
  for (coordinate in seq_len(2L)) {
    cols <- ((coordinate - 1L) * p_grid + 1L):(coordinate * p_grid)
    for (coefficient in seq_len(n_coef)) {
      fit <- stats::smooth.spline(shape_grid, raw_coef[coefficient, cols])
      smooth_coef[coefficient, , coordinate] <-
        stats::predict(fit, x = shape_grid)$y
      smoothing[[counter]] <- c(spar = fit$spar, df = fit$df)
      counter <- counter + 1L
    }
  }

  list(coefficients = smooth_coef, smoothing = smoothing)
}

fit_shape_layer <- function(x, W_matrix, shape_array, shape_grid) {
  p_grid <- length(shape_grid)
  design <- cbind(SNP = as.numeric(x), W_matrix)
  design_qr <- qr(design)
  if (design_qr$rank < ncol(design)) {
    stop("The first-stage design matrix is rank deficient.")
  }

  shape_matrix <- cbind(shape_array[, , 1L], shape_array[, , 2L])
  raw_coef <- qr.coef(design_qr, shape_matrix)
  rownames(raw_coef) <- colnames(design)
  smoothed <- smooth_shape_coefficients(raw_coef, shape_grid, p_grid)

  raw_array <- array(
    NA_real_,
    dim = c(nrow(raw_coef), p_grid, 2L),
    dimnames = dimnames(smoothed$coefficients)
  )
  raw_array[, , 1L] <- raw_coef[, seq_len(p_grid), drop = FALSE]
  raw_array[, , 2L] <- raw_coef[, p_grid + seq_len(p_grid), drop = FALSE]

  fitted_matrix <- design %*% raw_coef
  residual_matrix <- shape_matrix - fitted_matrix

  list(
    design_names = colnames(design),
    coefficients_raw = raw_array,
    coefficients_smoothed = smoothed$coefficients,
    smoothing = smoothed$smoothing,
    alpha = smoothed$coefficients["SNP", , , drop = TRUE],
    residual_sd = apply(residual_matrix, 2L, stats::sd)
  )
}

fit_joint_fpca <- function(
  shape_array,
  shape_grid,
  variance_threshold = 0.85,
  n_components = NULL
) {
  p_grid <- length(shape_grid)
  weights_grid <- trapezoid_weights(shape_grid)
  weights_flat <- rep(weights_grid, times = 2L)
  shape_matrix <- cbind(shape_array[, , 1L], shape_array[, , 2L])
  shape_mean <- colMeans(shape_matrix)
  centered <- sweep(shape_matrix, 2L, shape_mean, FUN = "-")
  weighted <- sweep(centered, 2L, sqrt(weights_flat), FUN = "*")
  decomposition <- svd(weighted, nu = 0L)
  eigenvalues <- decomposition$d^2 / (nrow(shape_matrix) - 1L)
  variance_fraction <- eigenvalues / sum(eigenvalues)
  cumulative_variance <- cumsum(variance_fraction)
  if (is.null(n_components)) {
    n_components <- which(cumulative_variance >= variance_threshold)[1L]
  } else {
    n_components <- as.integer(n_components)
    if (
      length(n_components) != 1L || is.na(n_components) ||
      n_components < 1L || n_components > ncol(decomposition$v)
    ) {
      stop("n_components is outside the available FPCA range.")
    }
  }

  weighted_eigenvectors <- decomposition$v[, seq_len(n_components), drop = FALSE]
  eigenfunctions <- sweep(
    weighted_eigenvectors,
    1L,
    sqrt(weights_flat),
    FUN = "/"
  )
  scores <- weighted %*% weighted_eigenvectors
  colnames(scores) <- sprintf("FPC%02d", seq_len(n_components))

  list(
    scores = scores,
    eigenfunctions = eigenfunctions,
    eigenvalues = eigenvalues,
    variance_fraction = variance_fraction,
    cumulative_variance = cumulative_variance,
    n_components = n_components,
    retained_variance = cumulative_variance[n_components],
    shape_mean = shape_mean,
    grid_weights = weights_grid,
    flat_weights = weights_flat,
    p_grid = p_grid
  )
}

fit_cox_layer <- function(x, W_matrix, subject_data, fpca_fit) {
  W_cox <- W_matrix[, setdiff(colnames(W_matrix), "intercept"), drop = FALSE]
  cox_data <- data.frame(
    time = as.numeric(subject_data$time),
    event = as.integer(subject_data$event),
    SNP = as.numeric(x),
    W_cox,
    fpca_fit$scores,
    check.names = FALSE
  )
  fit_warnings <- character()
  fit <- withCallingHandlers(
    survival::coxph(
      survival::Surv(time, event) ~ .,
      data = cox_data,
      ties = "breslow",
      x = TRUE,
      y = TRUE,
      model = TRUE
    ),
    warning = function(w) {
      fit_warnings <<- c(fit_warnings, conditionMessage(w))
      invokeRestart("muffleWarning")
    }
  )
  if (any(!is.finite(stats::coef(fit)))) {
    stop("The Cox model did not converge to finite coefficients.")
  }

  fpc_names <- colnames(fpca_fit$scores)
  theta <- stats::coef(fit)[fpc_names]
  beta_flat <- as.numeric(fpca_fit$eigenfunctions %*% theta)
  p_grid <- fpca_fit$p_grid
  beta <- cbind(
    coordinate_1 = beta_flat[seq_len(p_grid)],
    coordinate_2 = beta_flat[p_grid + seq_len(p_grid)]
  )
  fit_summary <- summary(fit)

  list(
    fit = fit,
    tau = unname(stats::coef(fit)["SNP"]),
    kappa = stats::coef(fit)[colnames(W_cox)],
    fpc_coefficients = theta,
    beta = beta,
    baseline_hazard = survival::basehaz(fit, centered = FALSE),
    concordance = unname(fit_summary$concordance[1L]),
    warnings = unique(fit_warnings)
  )
}

run_point_analysis <- function() {
fpca <- fit_joint_fpca(M, grid, variance_threshold = 0.85)
model_results <- setNames(vector("list", length(snp_names)), snp_names)
effect_rows <- vector("list", length(snp_names))
scalar_rows <- vector("list", length(snp_names))
shape_coefficient_rows <- vector("list", length(snp_names))
cox_coefficient_rows <- vector("list", length(snp_names))

for (index in seq_along(snp_names)) {
  snp <- snp_names[index]
  message("Fitting ", snp, " (", index, "/", length(snp_names), ")")
  x <- X[, snp]
  if (!all(x %in% c(0, 1, 2))) {
    stop("SNP ", snp, " is not coded as additive dosage 0/1/2.")
  }

  shape_fit <- fit_shape_layer(x, W, M, grid)
  cox_fit <- fit_cox_layer(x, W, subject, fpca)
  alpha <- shape_fit$alpha
  beta <- cox_fit$beta
  spatial_indirect <- alpha[, 1L] * beta[, 1L] + alpha[, 2L] * beta[, 2L]
  indirect_coordinate_1 <- sum(fpca$grid_weights * alpha[, 1L] * beta[, 1L])
  indirect_coordinate_2 <- sum(fpca$grid_weights * alpha[, 2L] * beta[, 2L])
  indirect_log_hr <- indirect_coordinate_1 + indirect_coordinate_2
  direct_log_hr <- cox_fit$tau
  total_log_hr <- direct_log_hr + indirect_log_hr

  effect_rows[[index]] <- data.frame(
    SNP = snp,
    n = n,
    events = sum(subject$event),
    effect_contrast_alleles = 1,
    direct_log_HR = direct_log_hr,
    direct_HR = exp(direct_log_hr),
    indirect_log_HR_coordinate_1 = indirect_coordinate_1,
    indirect_log_HR_coordinate_2 = indirect_coordinate_2,
    indirect_log_HR = indirect_log_hr,
    indirect_HR = exp(indirect_log_hr),
    total_log_HR = total_log_hr,
    total_HR = exp(total_log_hr),
    cox_concordance = cox_fit$concordance,
    functional_components = fpca$n_components,
    functional_variance_retained = fpca$retained_variance
  )

  all_scalar <- stats::coef(cox_fit$fit)[c("SNP", setdiff(colnames(W), "intercept"))]
  scalar_rows[[index]] <- data.frame(
    SNP_model = snp,
    parameter = names(all_scalar),
    log_HR_coefficient = unname(all_scalar),
    HR = exp(unname(all_scalar))
  )

  shape_coefficients <- shape_fit$coefficients_smoothed
  shape_coefficients_raw <- shape_fit$coefficients_raw
  shape_coefficient_table <- do.call(
    rbind,
    lapply(seq_len(nrow(shape_coefficients)), function(parameter_index) {
      data.frame(
        SNP_model = snp,
        predictor = rownames(shape_coefficients)[parameter_index],
        grid = grid,
        coordinate_1 = shape_coefficients[parameter_index, , 1L],
        coordinate_2 = shape_coefficients[parameter_index, , 2L],
        coordinate_1_raw = shape_coefficients_raw[parameter_index, , 1L],
        coordinate_2_raw = shape_coefficients_raw[parameter_index, , 2L],
        row.names = NULL
      )
    })
  )
  shape_coefficient_rows[[index]] <- shape_coefficient_table
  utils::write.csv(
    shape_coefficient_table,
    file.path(output_dir, paste0(snp, "_shape_coefficients.csv")),
    row.names = FALSE
  )

  cox_coefficients <- stats::coef(cox_fit$fit)
  cox_coefficient_table <- data.frame(
    SNP_model = snp,
    term = names(cox_coefficients),
    coefficient = unname(cox_coefficients),
    HR_per_unit = exp(unname(cox_coefficients)),
    row.names = NULL
  )
  cox_coefficient_rows[[index]] <- cox_coefficient_table
  utils::write.csv(
    cox_coefficient_table,
    file.path(output_dir, paste0(snp, "_cox_coefficients.csv")),
    row.names = FALSE
  )

  function_table <- data.frame(
    grid = grid,
    alpha_coordinate_1 = alpha[, 1L],
    alpha_coordinate_2 = alpha[, 2L],
    beta_coordinate_1 = beta[, 1L],
    beta_coordinate_2 = beta[, 2L],
    indirect_integrand_coordinate_1 = alpha[, 1L] * beta[, 1L],
    indirect_integrand_coordinate_2 = alpha[, 2L] * beta[, 2L],
    indirect_integrand_total = spatial_indirect
  )
  utils::write.csv(
    function_table,
    file.path(output_dir, paste0(snp, "_functional_parameters.csv")),
    row.names = FALSE
  )

  model_results[[snp]] <- list(
    snp = snp,
    dosage_counts = table(factor(x, levels = 0:2)),
    shape_layer = shape_fit,
    cox_layer = cox_fit,
    effects = effect_rows[[index]],
    shape_coefficients = shape_coefficient_table,
    cox_coefficients = cox_coefficient_table,
    functional_parameters = function_table
  )
}

effect_summary <- do.call(rbind, effect_rows)
scalar_parameters <- do.call(rbind, scalar_rows)
shape_coefficients_all <- do.call(rbind, shape_coefficient_rows)
cox_coefficients_all <- do.call(rbind, cox_coefficient_rows)

analysis_result <- list(
  effects = effect_summary,
  scalar_parameters = scalar_parameters,
  models = model_results,
  fpca = fpca,
  settings = list(
    input_file = input_file,
    subjects = n,
    events = sum(subject$event),
    confounders = setdiff(colnames(W), "intercept"),
    exposure_coding = "additive SNP dosage 0/1/2; effect contrast is one allele",
    first_layer = "pointwise multivariate OLS with GCV smoothing splines for coefficient functions",
    second_layer = "Cox FLCRM using joint weighted FPCA of the two shape coordinates",
    fpca_variance_threshold = 0.85,
    fpca_components = fpca$n_components,
    fpca_variance_retained = fpca$retained_variance,
    cox_ties = "breslow",
    inference = "point estimates only; no bootstrap or sensitivity analysis"
  )
)

utils::write.csv(
  effect_summary,
  file.path(output_dir, "four_snp_effect_summary.csv"),
  row.names = FALSE
)
utils::write.csv(
  scalar_parameters,
  file.path(output_dir, "four_snp_scalar_parameters.csv"),
  row.names = FALSE
)
utils::write.csv(
  shape_coefficients_all,
  file.path(output_dir, "four_snp_shape_coefficients.csv"),
  row.names = FALSE
)
utils::write.csv(
  cox_coefficients_all,
  file.path(output_dir, "four_snp_cox_coefficients.csv"),
  row.names = FALSE
)
saveRDS(
  analysis_result,
  file.path(output_dir, "four_snp_model_results.rds"),
  compress = "xz"
)

message("Completed. Results written to: ", output_dir)
print(effect_summary, row.names = FALSE)
invisible(analysis_result)
}

if (sys.nframe() == 0L) {
  run_point_analysis()
}
