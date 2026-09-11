hima2_survival_fixed <- function(formula, data.pheno, data.M, FDRcut = 0.3, scale = TRUE, verbose = FALSE) {
  lhs_vars <- all.vars(formula[[2]])
  rhs_vars <- attr(stats::terms(formula), "term.labels")

  if (length(lhs_vars) != 2) {
    stop("Use a survival formula of the form Surv(Time, Status) ~ Exposure + covariates.")
  }
  if (length(rhs_vars) < 1) {
    stop("The formula must include an exposure variable on the right-hand side.")
  }

  time_var <- lhs_vars[1]
  status_var <- lhs_vars[2]
  exposure_var <- rhs_vars[1]
  covariate_vars <- rhs_vars[-1]

  Z <- if (length(covariate_vars) > 0) data.pheno[, covariate_vars, drop = FALSE] else NULL
  res <- survHIMA_fixed(
    X = data.pheno[[exposure_var]],
    Z = Z,
    M = data.M,
    OT = data.pheno[[time_var]],
    status = data.pheno[[status_var]],
    FDRcut = FDRcut,
    scale = scale,
    verbose = verbose
  )

  data.frame(
    ID = res$ID,
    alpha = res$alpha,
    alpha_se = res$alpha_se,
    beta = res$beta,
    beta_se = res$beta_se,
    p.joint = res$p.joint,
    stringsAsFactors = FALSE
  )
}
