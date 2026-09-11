empirical_two_sided_p <- function(values, null_value = 0) {
  valid_values <- values[is.finite(values)]
  n_valid <- length(valid_values)
  if (n_valid == 0) {
    return(NA_real_)
  }

  p_val <- 2 * min(
    mean(valid_values >= null_value),
    mean(valid_values <= null_value)
  )
  p_val <- max(p_val, 1 / n_valid)
  min(p_val, 1)
}

summarize_bootstrap_metric <- function(metric_name, estimate, draws) {
  valid_draws <- draws[is.finite(draws)]
  n_valid <- length(valid_draws)
  if (n_valid == 0) {
    return(data.frame(
      metric = metric_name,
      estimate = estimate,
      bootstrap_se = NA_real_,
      ci_low = NA_real_,
      ci_high = NA_real_,
      p_value = NA_real_,
      n_valid = 0L,
      stringsAsFactors = FALSE
    ))
  }

  data.frame(
    metric = metric_name,
    estimate = estimate,
    bootstrap_se = stats::sd(valid_draws),
    ci_low = as.numeric(stats::quantile(valid_draws, 0.025, names = FALSE)),
    ci_high = as.numeric(stats::quantile(valid_draws, 0.975, names = FALSE)),
    p_value = empirical_two_sided_p(valid_draws, null_value = 0),
    n_valid = n_valid,
    stringsAsFactors = FALSE
  )
}

fit_mediator_models <- function(basis_scores, sample_data, gene_name, covariate_names) {
  lambda_rows <- vector("list", length = ncol(basis_scores))
  mediator_summary_rows <- vector("list", length = ncol(basis_scores))
  fit_list <- vector("list", length = ncol(basis_scores))

  for (k in seq_len(ncol(basis_scores))) {
    dat <- data.frame(
      Response = basis_scores[, k],
      Exposure = sample_data[[gene_name]],
      sample_data[, covariate_names, drop = FALSE],
      stringsAsFactors = FALSE
    )

    fit <- stats::lm(Response ~ Exposure + ., data = dat)
    fit_list[[k]] <- fit

    coef_tab <- summary(fit)$coefficients
    if (!("Exposure" %in% rownames(coef_tab))) {
      stop("Exposure coefficient missing in mediator model for ", gene_name, ", basis ", k, ".")
    }

    lambda_rows[[k]] <- data.frame(
      Gene = gene_name,
      Basis = colnames(basis_scores)[k],
      basis_index = k,
      parameter = "lambda",
      estimate = unname(coef_tab["Exposure", "Estimate"]),
      std_error = unname(coef_tab["Exposure", "Std. Error"]),
      statistic = unname(coef_tab["Exposure", "t value"]),
      p_value = unname(coef_tab["Exposure", "Pr(>|t|)"]),
      stringsAsFactors = FALSE
    )

    mediator_summary_rows[[k]] <- data.frame(
      Gene = gene_name,
      Basis = colnames(basis_scores)[k],
      basis_index = k,
      term = rownames(coef_tab),
      estimate = coef_tab[, "Estimate"],
      std_error = coef_tab[, "Std. Error"],
      statistic = coef_tab[, "t value"],
      p_value = coef_tab[, "Pr(>|t|)"],
      stringsAsFactors = FALSE
    )
  }

  list(
    lambda_table = do.call(rbind, lambda_rows),
    mediator_model_summary = do.call(rbind, mediator_summary_rows),
    fits = fit_list
  )
}

fit_survival_model <- function(basis_scores, sample_data, gene_name, covariate_names) {
  dat <- data.frame(
    Time = sample_data$Time,
    Status = sample_data$Status,
    Exposure = sample_data[[gene_name]],
    sample_data[, covariate_names, drop = FALSE],
    basis_scores,
    stringsAsFactors = FALSE
  )

  rhs_terms <- c("Exposure", covariate_names, colnames(basis_scores))
  cox_formula <- stats::as.formula(
    paste("survival::Surv(Time, Status) ~", paste(rhs_terms, collapse = " + "))
  )

  fit <- survival::coxph(
    cox_formula,
    data = dat,
    ties = "efron",
    model = TRUE,
    x = TRUE
  )

  coef_tab <- summary(fit)$coefficients
  if (!("Exposure" %in% rownames(coef_tab))) {
    stop("Exposure coefficient missing in Cox model for ", gene_name, ".")
  }

  alpha_row <- data.frame(
    Gene = gene_name,
    parameter = "alpha",
    estimate = unname(coef_tab["Exposure", "coef"]),
    std_error = unname(coef_tab["Exposure", "se(coef)"]),
    statistic = unname(coef_tab["Exposure", "z"]),
    p_value = unname(coef_tab["Exposure", "Pr(>|z|)"]),
    stringsAsFactors = FALSE
  )

  basis_idx <- match(colnames(basis_scores), rownames(coef_tab))
  if (any(is.na(basis_idx))) {
    stop("Basis coefficients missing in Cox model for ", gene_name, ".")
  }

  beta_table <- data.frame(
    Gene = gene_name,
    Basis = colnames(basis_scores),
    basis_index = seq_len(ncol(basis_scores)),
    parameter = "beta",
    estimate = coef_tab[basis_idx, "coef"],
    std_error = coef_tab[basis_idx, "se(coef)"],
    statistic = coef_tab[basis_idx, "z"],
    p_value = coef_tab[basis_idx, "Pr(>|z|)"],
    stringsAsFactors = FALSE
  )

  cox_model_summary <- data.frame(
    Gene = gene_name,
    term = rownames(coef_tab),
    estimate = coef_tab[, "coef"],
    hazard_ratio = coef_tab[, "exp(coef)"],
    std_error = coef_tab[, "se(coef)"],
    statistic = coef_tab[, "z"],
    p_value = coef_tab[, "Pr(>|z|)"],
    stringsAsFactors = FALSE
  )

  list(
    fit = fit,
    alpha_row = alpha_row,
    beta_table = beta_table,
    cox_model_summary = cox_model_summary
  )
}

compute_gene_effects <- function(gene_name, alpha_row, lambda_table, beta_table) {
  beta_aligned <- beta_table[match(lambda_table$Basis, beta_table$Basis), , drop = FALSE]
  if (any(is.na(beta_aligned$estimate))) {
    stop("Unable to align lambda and beta coefficients for ", gene_name, ".")
  }

  direct_effect <- alpha_row$estimate[1]
  indirect_effect <- sum(lambda_table$estimate * beta_aligned$estimate)
  total_effect <- direct_effect + indirect_effect

  data.frame(
    Gene = gene_name,
    direct_effect = direct_effect,
    direct_effect_hr = exp(direct_effect),
    indirect_effect = indirect_effect,
    indirect_effect_hr = exp(indirect_effect),
    total_effect = total_effect,
    total_effect_hr = exp(total_effect),
    stringsAsFactors = FALSE
  )
}

validate_gene_fit <- function(gene_name, alpha_row, lambda_table, beta_table, gene_effects) {
  numeric_blocks <- list(
    alpha = alpha_row$estimate,
    lambda = lambda_table$estimate,
    beta = beta_table$estimate,
    direct_effect = gene_effects$direct_effect,
    indirect_effect = gene_effects$indirect_effect,
    total_effect = gene_effects$total_effect
  )

  has_bad_values <- vapply(
    numeric_blocks,
    function(x) any(!is.finite(x)),
    logical(1)
  )
  if (any(has_bad_values)) {
    stop(
      "Non-finite estimate encountered for ", gene_name, ": ",
      paste(names(has_bad_values)[has_bad_values], collapse = ", ")
    )
  }
}

fit_functional_gene_model <- function(curve_matrix, sample_data, basis_obj, gene_name, covariate_names) {
  basis_scores <- project_curves_to_basis(curve_matrix, basis_obj)
  mediator_fit <- fit_mediator_models(
    basis_scores = basis_scores,
    sample_data = sample_data,
    gene_name = gene_name,
    covariate_names = covariate_names
  )
  survival_fit <- fit_survival_model(
    basis_scores = basis_scores,
    sample_data = sample_data,
    gene_name = gene_name,
    covariate_names = covariate_names
  )

  basis_coefficients <- rbind(
    mediator_fit$lambda_table,
    survival_fit$beta_table
  )
  gene_effects <- compute_gene_effects(
    gene_name = gene_name,
    alpha_row = survival_fit$alpha_row,
    lambda_table = mediator_fit$lambda_table,
    beta_table = survival_fit$beta_table
  )
  validate_gene_fit(
    gene_name = gene_name,
    alpha_row = survival_fit$alpha_row,
    lambda_table = mediator_fit$lambda_table,
    beta_table = survival_fit$beta_table,
    gene_effects = gene_effects
  )

  list(
    basis_scores = basis_scores,
    lambda_table = mediator_fit$lambda_table,
    beta_table = survival_fit$beta_table,
    alpha_row = survival_fit$alpha_row,
    basis_coefficients = basis_coefficients,
    gene_effects = gene_effects,
    mediator_model_summary = mediator_fit$mediator_model_summary,
    cox_model_summary = survival_fit$cox_model_summary
  )
}

run_gene_bootstrap <- function(curve_matrix,
                               sample_data,
                               basis_obj,
                               gene_name,
                               covariate_names,
                               bootstrap_reps = 100,
                               random_seed = 123) {
  set.seed(random_seed)

  n_subjects <- nrow(sample_data)
  bootstrap_rows <- vector("list", bootstrap_reps)

  for (b in seq_len(bootstrap_reps)) {
    sample_idx <- sample.int(n_subjects, size = n_subjects, replace = TRUE)

    bootstrap_rows[[b]] <- tryCatch(
      {
        fit <- suppressWarnings(
          fit_functional_gene_model(
            curve_matrix = curve_matrix[sample_idx, , drop = FALSE],
            sample_data = sample_data[sample_idx, , drop = FALSE],
            basis_obj = basis_obj,
            gene_name = gene_name,
            covariate_names = covariate_names
          )
        )

        data.frame(
          replicate = b,
          direct_effect = fit$gene_effects$direct_effect[1],
          indirect_effect = fit$gene_effects$indirect_effect[1],
          total_effect = fit$gene_effects$total_effect[1],
          status = "ok",
          message = NA_character_,
          stringsAsFactors = FALSE
        )
      },
      error = function(e) {
        data.frame(
          replicate = b,
          direct_effect = NA_real_,
          indirect_effect = NA_real_,
          total_effect = NA_real_,
          status = "failed",
          message = conditionMessage(e),
          stringsAsFactors = FALSE
        )
      }
    )
  }

  do.call(rbind, bootstrap_rows)
}

summarize_gene_bootstrap <- function(gene_name, full_sample_effects, bootstrap_draws) {
  direct_summary <- summarize_bootstrap_metric(
    metric_name = "direct_effect",
    estimate = full_sample_effects$direct_effect[1],
    draws = bootstrap_draws$direct_effect
  )
  indirect_summary <- summarize_bootstrap_metric(
    metric_name = "indirect_effect",
    estimate = full_sample_effects$indirect_effect[1],
    draws = bootstrap_draws$indirect_effect
  )
  total_summary <- summarize_bootstrap_metric(
    metric_name = "total_effect",
    estimate = full_sample_effects$total_effect[1],
    draws = bootstrap_draws$total_effect
  )

  bootstrap_summary <- rbind(direct_summary, indirect_summary, total_summary)
  bootstrap_summary$Gene <- gene_name
  bootstrap_summary$n_attempted <- nrow(bootstrap_draws)
  bootstrap_summary$n_failed <- sum(bootstrap_draws$status != "ok")
  bootstrap_summary
}

assemble_gene_level_summary <- function(full_sample_effects, bootstrap_summary) {
  get_metric_col <- function(metric_name, column_name) {
    bootstrap_summary[bootstrap_summary$metric == metric_name, column_name][1]
  }

  data.frame(
    Gene = full_sample_effects$Gene[1],
    direct_effect = full_sample_effects$direct_effect[1],
    direct_effect_hr = full_sample_effects$direct_effect_hr[1],
    direct_bootstrap_se = get_metric_col("direct_effect", "bootstrap_se"),
    direct_ci_low = get_metric_col("direct_effect", "ci_low"),
    direct_ci_high = get_metric_col("direct_effect", "ci_high"),
    direct_p_value = get_metric_col("direct_effect", "p_value"),
    indirect_effect = full_sample_effects$indirect_effect[1],
    indirect_effect_hr = full_sample_effects$indirect_effect_hr[1],
    indirect_bootstrap_se = get_metric_col("indirect_effect", "bootstrap_se"),
    indirect_ci_low = get_metric_col("indirect_effect", "ci_low"),
    indirect_ci_high = get_metric_col("indirect_effect", "ci_high"),
    indirect_p_value = get_metric_col("indirect_effect", "p_value"),
    total_effect = full_sample_effects$total_effect[1],
    total_effect_hr = full_sample_effects$total_effect_hr[1],
    total_bootstrap_se = get_metric_col("total_effect", "bootstrap_se"),
    total_ci_low = get_metric_col("total_effect", "ci_low"),
    total_ci_high = get_metric_col("total_effect", "ci_high"),
    total_p_value = get_metric_col("total_effect", "p_value"),
    bootstrap_attempted = get_metric_col("direct_effect", "n_attempted"),
    bootstrap_failed = get_metric_col("direct_effect", "n_failed"),
    bootstrap_valid = get_metric_col("direct_effect", "n_valid"),
    stringsAsFactors = FALSE
  )
}
