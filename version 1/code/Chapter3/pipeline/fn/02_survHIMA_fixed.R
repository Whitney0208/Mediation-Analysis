safe_screen_beta <- function(mediator, X, Z, OT, status) {
  dat <- data.frame(
    OT = OT,
    status = status,
    mediator = mediator,
    exposure = as.numeric(X)
  )
  if (!is.null(Z)) {
    dat <- cbind(dat, as.data.frame(Z))
  }

  fit <- tryCatch(
    suppressWarnings(
      survival::coxph(
        survival::Surv(OT, status) ~ .,
        data = dat,
        singular.ok = TRUE,
        ties = "efron",
        model = FALSE,
        x = FALSE,
        y = FALSE
      )
    ),
    error = function(e) NULL
  )

  if (is.null(fit)) {
    return(NA_real_)
  }

  coef_vec <- stats::coef(fit)
  if (is.null(coef_vec) || !("mediator" %in% names(coef_vec))) {
    return(NA_real_)
  }

  beta <- unname(coef_vec[["mediator"]])
  if (!is.finite(beta)) {
    return(NA_real_)
  }
  beta
}

safe_alpha_fit <- function(mediator, X, Z) {
  dat <- data.frame(
    mediator = mediator,
    exposure = as.numeric(X)
  )
  if (!is.null(Z)) {
    dat <- cbind(dat, as.data.frame(Z))
  }

  fit <- tryCatch(
    suppressWarnings(stats::lm(mediator ~ ., data = dat)),
    error = function(e) NULL
  )
  if (is.null(fit)) {
    return(c(est = NA_real_, se = NA_real_, p = NA_real_))
  }

  coef_tab <- summary(fit)$coefficients
  if (!("exposure" %in% rownames(coef_tab))) {
    return(c(est = NA_real_, se = NA_real_, p = NA_real_))
  }

  c(
    est = unname(coef_tab["exposure", "Estimate"]),
    se = unname(coef_tab["exposure", "Std. Error"]),
    p = unname(coef_tab["exposure", "Pr(>|t|)"])
  )
}

survHIMA_fixed <- function(X, Z, M, OT, status, FDRcut = 0.3, scale = TRUE, verbose = FALSE) {
  X <- matrix(X, ncol = 1)
  M <- as.matrix(M)

  mediator_names <- colnames(M)
  if (is.null(mediator_names)) {
    mediator_names <- paste0("Mediator", seq_len(ncol(M)))
  }

  nuisance_idx <- if (is.null(Z)) integer(0) else seq_len(ncol(as.matrix(Z)))
  MZ <- if (is.null(Z)) cbind(M, X) else cbind(M, as.matrix(Z), X)
  if (scale) {
    MZ <- scale(MZ)
  }

  n <- nrow(M)
  p <- ncol(M)
  q <- length(nuisance_idx)

  message("Step 1: Sure Independent Screening ...", "     (", format(Sys.time(), "%X"), ")")

  beta_SIS <- vapply(
    seq_len(p),
    function(i) safe_screen_beta(M[, i], X = X[, 1], Z = Z, OT = OT, status = status),
    numeric(1)
  )
  alpha_screen <- t(vapply(
    seq_len(p),
    function(i) safe_alpha_fit(M[, i], X = X[, 1], Z = Z),
    numeric(3)
  ))

  alpha_SIS <- alpha_screen[, "est"]
  valid_screen <- which(is.finite(alpha_SIS) & is.finite(beta_SIS))
  if (length(valid_screen) == 0) {
    message("No mediators survived the screening step.", "     (", format(Sys.time(), "%X"), ")")
    return(empty_survhima_result())
  }

  d_0 <- max(1, round(n / log(n)))
  screen_score <- abs(alpha_SIS[valid_screen] * beta_SIS[valid_screen])
  keep_n <- min(length(valid_screen), d_0)
  ID_SIS <- valid_screen[order(screen_score, decreasing = TRUE)[seq_len(keep_n)]]

  d <- length(ID_SIS)
  if (verbose) {
    message("        ", d, " mediators selected from the screening.")
  }

  message("Step 2: De-biased Lasso estimates ...", "     (", format(Sys.time(), "%X"), ")")

  P_beta_SIS <- rep(NA_real_, d)
  beta_DLASSO_SIS_est <- rep(NA_real_, d)
  beta_DLASSO_SIS_SE <- rep(NA_real_, d)
  MZ_SIS <- MZ[, c(ID_SIS, (p + 1):(p + q + 1)), drop = FALSE]
  MZ_SIS_1 <- MZ_SIS[, 1]

  for (i in seq_len(d)) {
    V <- MZ_SIS
    V[, 1] <- V[, i]
    V[, i] <- MZ_SIS_1

    ldpe_res <- tryCatch(
      LDPE_func(ID = 1, X = V, OT = OT, status = status),
      error = function(e) c(est = NA_real_, se = NA_real_)
    )

    beta_est <- as.numeric(ldpe_res["est"])
    beta_se <- as.numeric(ldpe_res["se"])
    if (is.finite(beta_est) && is.finite(beta_se) && beta_se > 0) {
      z_score <- abs(beta_est) / beta_se
      P_beta_SIS[i] <- 2 * (1 - stats::pnorm(z_score))
      beta_DLASSO_SIS_est[i] <- beta_est
      beta_DLASSO_SIS_SE[i] <- beta_se
    }
  }

  alpha_SIS_est <- rep(NA_real_, d)
  alpha_SIS_SE <- rep(NA_real_, d)
  P_alpha_SIS <- rep(NA_real_, d)
  for (i in seq_len(d)) {
    alpha_fit <- safe_alpha_fit(M[, ID_SIS[i]], X = X[, 1], Z = Z)
    alpha_SIS_est[i] <- alpha_fit["est"]
    alpha_SIS_SE[i] <- alpha_fit["se"]
    P_alpha_SIS[i] <- alpha_fit["p"]
  }

  message("Step 3: Multiple-testing procedure ...", "     (", format(Sys.time(), "%X"), ")")

  PA <- cbind(P_alpha_SIS, P_beta_SIS)
  valid_joint <- which(stats::complete.cases(PA))
  if (length(valid_joint) == 0) {
    message("No mediators survived the inference step.", "     (", format(Sys.time(), "%X"), ")")
    return(empty_survhima_result())
  }

  PA_valid <- PA[valid_joint, , drop = FALSE]
  P_value <- apply(PA_valid, 1, max)
  N0 <- nrow(PA_valid) * ncol(PA_valid)
  input_pvalues <- PA_valid + matrix(stats::runif(N0, 0, 1e-10), nrow(PA_valid), 2)

  fdrcut <- tryCatch(
    {
      nullprop <- HIMA:::null_estimation(input_pvalues, lambda = 0.5)
      HDMT::fdr_est(
        nullprop$alpha00,
        nullprop$alpha01,
        nullprop$alpha10,
        nullprop$alpha1,
        nullprop$alpha2,
        input_pvalues,
        exact = 0
      )
    },
    error = function(e) stats::p.adjust(P_value, method = "BH")
  )

  ID_fdr_local <- which(fdrcut <= FDRcut)
  if (length(ID_fdr_local) == 0) {
    message("No mediators passed the multiple-testing threshold.", "     (", format(Sys.time(), "%X"), ")")
    return(empty_survhima_result())
  }

  selected_local <- valid_joint[ID_fdr_local]
  selected_global <- ID_SIS[selected_local]

  out_result <- data.frame(
    ID = mediator_names[selected_global],
    alpha = alpha_SIS_est[selected_local],
    alpha_se = alpha_SIS_SE[selected_local],
    beta = beta_DLASSO_SIS_est[selected_local],
    beta_se = beta_DLASSO_SIS_SE[selected_local],
    p.joint = P_value[ID_fdr_local],
    stringsAsFactors = FALSE
  )

  message("Done!", "     (", format(Sys.time(), "%X"), ")")
  out_result
}

survHIMA_fixed_diagnostics <- function(X, Z, M, OT, status, scale = TRUE, verbose = FALSE) {
  X <- matrix(X, ncol = 1)
  M <- as.matrix(M)

  mediator_names <- colnames(M)
  if (is.null(mediator_names)) {
    mediator_names <- paste0("Mediator", seq_len(ncol(M)))
  }

  nuisance_idx <- if (is.null(Z)) integer(0) else seq_len(ncol(as.matrix(Z)))
  MZ <- if (is.null(Z)) cbind(M, X) else cbind(M, as.matrix(Z), X)
  if (scale) {
    MZ <- scale(MZ)
  }

  n <- nrow(M)
  p <- ncol(M)
  q <- length(nuisance_idx)

  message("Step 1: Sure Independent Screening ...", "     (", format(Sys.time(), "%X"), ")")

  beta_SIS <- vapply(
    seq_len(p),
    function(i) safe_screen_beta(M[, i], X = X[, 1], Z = Z, OT = OT, status = status),
    numeric(1)
  )
  alpha_screen <- t(vapply(
    seq_len(p),
    function(i) safe_alpha_fit(M[, i], X = X[, 1], Z = Z),
    numeric(3)
  ))

  alpha_SIS <- alpha_screen[, "est"]
  valid_screen <- which(is.finite(alpha_SIS) & is.finite(beta_SIS))
  if (length(valid_screen) == 0) {
    message("No mediators survived the screening step.", "     (", format(Sys.time(), "%X"), ")")
    return(list(
      selected = empty_survhima_result(),
      candidates = empty_survhima_result(),
      screening_n = 0L
    ))
  }

  d_0 <- max(1, round(n / log(n)))
  screen_score <- abs(alpha_SIS[valid_screen] * beta_SIS[valid_screen])
  keep_n <- min(length(valid_screen), d_0)
  ID_SIS <- valid_screen[order(screen_score, decreasing = TRUE)[seq_len(keep_n)]]

  d <- length(ID_SIS)
  if (verbose) {
    message("        ", d, " mediators selected from the screening.")
  }

  message("Step 2: De-biased Lasso estimates ...", "     (", format(Sys.time(), "%X"), ")")

  P_beta_SIS <- rep(NA_real_, d)
  beta_DLASSO_SIS_est <- rep(NA_real_, d)
  beta_DLASSO_SIS_SE <- rep(NA_real_, d)
  MZ_SIS <- MZ[, c(ID_SIS, (p + 1):(p + q + 1)), drop = FALSE]
  MZ_SIS_1 <- MZ_SIS[, 1]

  for (i in seq_len(d)) {
    V <- MZ_SIS
    V[, 1] <- V[, i]
    V[, i] <- MZ_SIS_1

    ldpe_res <- tryCatch(
      LDPE_func(ID = 1, X = V, OT = OT, status = status),
      error = function(e) c(est = NA_real_, se = NA_real_)
    )

    beta_est <- as.numeric(ldpe_res["est"])
    beta_se <- as.numeric(ldpe_res["se"])
    if (is.finite(beta_est) && is.finite(beta_se) && beta_se > 0) {
      z_score <- abs(beta_est) / beta_se
      P_beta_SIS[i] <- 2 * (1 - stats::pnorm(z_score))
      beta_DLASSO_SIS_est[i] <- beta_est
      beta_DLASSO_SIS_SE[i] <- beta_se
    }
  }

  alpha_SIS_est <- rep(NA_real_, d)
  alpha_SIS_SE <- rep(NA_real_, d)
  P_alpha_SIS <- rep(NA_real_, d)
  for (i in seq_len(d)) {
    alpha_fit <- safe_alpha_fit(M[, ID_SIS[i]], X = X[, 1], Z = Z)
    alpha_SIS_est[i] <- alpha_fit["est"]
    alpha_SIS_SE[i] <- alpha_fit["se"]
    P_alpha_SIS[i] <- alpha_fit["p"]
  }

  message("Step 3: Multiple-testing procedure ...", "     (", format(Sys.time(), "%X"), ")")

  PA <- cbind(P_alpha_SIS, P_beta_SIS)
  valid_joint <- which(stats::complete.cases(PA))
  if (length(valid_joint) == 0) {
    message("No mediators survived the inference step.", "     (", format(Sys.time(), "%X"), ")")
    return(list(
      selected = empty_survhima_result(),
      candidates = empty_survhima_result(),
      screening_n = d
    ))
  }

  PA_valid <- PA[valid_joint, , drop = FALSE]
  P_value <- apply(PA_valid, 1, max)
  N0 <- nrow(PA_valid) * ncol(PA_valid)
  input_pvalues <- PA_valid + matrix(stats::runif(N0, 0, 1e-10), nrow(PA_valid), 2)

  fdr_values <- tryCatch(
    {
      nullprop <- HIMA:::null_estimation(input_pvalues, lambda = 0.5)
      HDMT::fdr_est(
        nullprop$alpha00,
        nullprop$alpha01,
        nullprop$alpha10,
        nullprop$alpha1,
        nullprop$alpha2,
        input_pvalues,
        exact = 0
      )
    },
    error = function(e) stats::p.adjust(P_value, method = "BH")
  )

  selected_local <- valid_joint
  selected_global <- ID_SIS[selected_local]
  candidates <- data.frame(
    ID = mediator_names[selected_global],
    alpha = alpha_SIS_est[selected_local],
    alpha_se = alpha_SIS_SE[selected_local],
    beta = beta_DLASSO_SIS_est[selected_local],
    beta_se = beta_DLASSO_SIS_SE[selected_local],
    p.joint = P_value,
    fdr_value = fdr_values,
    stringsAsFactors = FALSE
  )
  candidates <- candidates[order(candidates$fdr_value, candidates$p.joint), , drop = FALSE]

  list(
    selected = candidates,
    candidates = candidates,
    screening_n = d
  )
}
