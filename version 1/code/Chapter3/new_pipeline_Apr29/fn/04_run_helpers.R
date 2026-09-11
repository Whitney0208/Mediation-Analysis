prepare_run_directories <- function(run_root) {
  combined_dir <- file.path(run_root, "combined")
  per_gene_dir <- file.path(run_root, "per_gene")
  diagnostics_dir <- file.path(run_root, "diagnostics")

  dir.create(combined_dir, recursive = TRUE, showWarnings = FALSE)
  dir.create(per_gene_dir, recursive = TRUE, showWarnings = FALSE)
  dir.create(diagnostics_dir, recursive = TRUE, showWarnings = FALSE)

  list(
    run_root = run_root,
    combined_dir = combined_dir,
    per_gene_dir = per_gene_dir,
    diagnostics_dir = diagnostics_dir
  )
}

write_basis_diagnostics <- function(chapter3, basis_obj, basis_scores, diagnostics_dir) {
  utils::write.csv(
    make_basis_grid_table(basis_obj),
    file = file.path(diagnostics_dir, "basis_grid.csv"),
    row.names = FALSE
  )
  save_basis_plot(basis_obj, file.path(diagnostics_dir, "basis_plot.png"))
  utils::write.csv(
    data.frame(
      SampleID = chapter3$sample_data$SampleID,
      basis_scores,
      check.names = FALSE
    ),
    file = file.path(diagnostics_dir, "basis_scores.csv"),
    row.names = FALSE
  )
}

run_functional_mediation_analysis <- function(chapter3,
                                              gene_names,
                                              basis_k,
                                              bootstrap_reps,
                                              random_seed,
                                              save_bootstrap_draws,
                                              run_root,
                                              output_label,
                                              write_combined_basis_summary = TRUE) {
  run_dirs <- prepare_run_directories(run_root)

  message("Building cubic spline basis with K = ", basis_k, "...")
  basis_obj <- build_spline_basis(chapter3$grid, basis_k = basis_k)
  orthogonality_check <- verify_basis_orthonormality(basis_obj)
  message("Basis max orthonormality deviation: ", signif(orthogonality_check$max_abs_deviation, 4))

  full_basis_scores <- project_curves_to_basis(chapter3$mediator_matrix, basis_obj)
  write_basis_diagnostics(
    chapter3 = chapter3,
    basis_obj = basis_obj,
    basis_scores = full_basis_scores,
    diagnostics_dir = run_dirs$diagnostics_dir
  )

  gene_summary_rows <- list()
  basis_summary_rows <- list()

  for (gene_idx in seq_along(gene_names)) {
    gene_name <- gene_names[gene_idx]
    message("Running functional mediation model for ", gene_name, "...")

    gene_dir <- file.path(run_dirs$per_gene_dir, gene_name)
    dir.create(gene_dir, recursive = TRUE, showWarnings = FALSE)

    full_fit <- fit_functional_gene_model(
      curve_matrix = chapter3$mediator_matrix,
      sample_data = chapter3$sample_data,
      basis_obj = basis_obj,
      gene_name = gene_name,
      covariate_names = chapter3$covariate_names
    )

    bootstrap_draws <- run_gene_bootstrap(
      curve_matrix = chapter3$mediator_matrix,
      sample_data = chapter3$sample_data,
      basis_obj = basis_obj,
      gene_name = gene_name,
      covariate_names = chapter3$covariate_names,
      bootstrap_reps = bootstrap_reps,
      random_seed = random_seed + gene_idx * 1000L
    )

    bootstrap_summary <- summarize_gene_bootstrap(
      gene_name = gene_name,
      full_sample_effects = full_fit$gene_effects,
      bootstrap_draws = bootstrap_draws
    )
    gene_level_summary <- assemble_gene_level_summary(
      full_sample_effects = full_fit$gene_effects,
      bootstrap_summary = bootstrap_summary
    )

    basis_coefficients <- full_fit$basis_coefficients
    basis_coefficients$hazard_ratio <- ifelse(
      basis_coefficients$parameter == "beta",
      exp(basis_coefficients$estimate),
      NA_real_
    )

    utils::write.csv(gene_level_summary, file = file.path(gene_dir, "gene_level_effects.csv"), row.names = FALSE)
    utils::write.csv(basis_coefficients, file = file.path(gene_dir, "basis_coefficients.csv"), row.names = FALSE)
    if (isTRUE(save_bootstrap_draws)) {
      utils::write.csv(bootstrap_draws, file = file.path(gene_dir, "bootstrap_effects.csv"), row.names = FALSE)
    } else {
      utils::write.csv(bootstrap_summary, file = file.path(gene_dir, "bootstrap_effects.csv"), row.names = FALSE)
    }
    utils::write.csv(full_fit$cox_model_summary, file = file.path(gene_dir, "cox_model_summary.csv"), row.names = FALSE)
    utils::write.csv(full_fit$mediator_model_summary, file = file.path(gene_dir, "mediator_model_summary.csv"), row.names = FALSE)

    gene_summary_rows[[gene_name]] <- gene_level_summary
    basis_summary_rows[[gene_name]] <- basis_coefficients
  }

  combined_gene_summary <- do.call(rbind, gene_summary_rows)
  combined_basis_summary <- do.call(rbind, basis_summary_rows)

  run_config <- data.frame(
    output_label = output_label,
    basis_k = basis_k,
    bootstrap_reps = bootstrap_reps,
    random_seed = random_seed,
    save_bootstrap_draws = save_bootstrap_draws,
    gene_names = paste(gene_names, collapse = ","),
    n_samples = nrow(chapter3$sample_data),
    n_grid = ncol(chapter3$mediator_matrix),
    removed_n = chapter3$removed_n,
    orthonormality_max_abs_deviation = orthogonality_check$max_abs_deviation,
    stringsAsFactors = FALSE
  )

  utils::write.csv(combined_gene_summary, file = file.path(run_dirs$combined_dir, "gene_level_summary.csv"), row.names = FALSE)
  if (isTRUE(write_combined_basis_summary)) {
    utils::write.csv(combined_basis_summary, file = file.path(run_dirs$combined_dir, "basis_level_summary.csv"), row.names = FALSE)
  }
  utils::write.csv(run_config, file = file.path(run_dirs$combined_dir, "run_config.csv"), row.names = FALSE)

  list(
    combined_gene_summary = combined_gene_summary,
    combined_basis_summary = combined_basis_summary,
    run_config = run_config,
    orthogonality_check = orthogonality_check,
    run_dirs = run_dirs
  )
}

format_ie_direction <- function(values) {
  direction_codes <- ifelse(
    values > 0, "+",
    ifelse(values < 0, "-", "0")
  )
  paste(direction_codes, collapse = "")
}

summarize_indirect_effect_stability <- function(gene_level_sensitivity_summary) {
  genes <- unique(gene_level_sensitivity_summary$Gene)
  stability_rows <- vector("list", length(genes))

  for (i in seq_along(genes)) {
    gene_name <- genes[i]
    gene_df <- gene_level_sensitivity_summary[
      gene_level_sensitivity_summary$Gene == gene_name,
      ,
      drop = FALSE
    ]
    gene_df <- gene_df[order(gene_df$basis_k), , drop = FALSE]

    ie_significant <- gene_df$indirect_p_value < 0.05
    stability_rows[[i]] <- data.frame(
      Gene = gene_name,
      basis_k_values = paste(gene_df$basis_k, collapse = ","),
      ie_direction_pattern = format_ie_direction(gene_df$indirect_effect),
      ie_min = min(gene_df$indirect_effect),
      ie_max = max(gene_df$indirect_effect),
      ie_range = diff(range(gene_df$indirect_effect)),
      any_ie_significant = any(ie_significant, na.rm = TRUE),
      significant_basis_k = if (any(ie_significant, na.rm = TRUE)) {
        paste(gene_df$basis_k[ie_significant], collapse = ",")
      } else {
        ""
      },
      significance_stability = if (!any(ie_significant, na.rm = TRUE)) {
        "none_significant"
      } else if (all(ie_significant, na.rm = TRUE)) {
        "significant_all_k"
      } else if (sum(ie_significant, na.rm = TRUE) == 1L) {
        "significant_one_k_only"
      } else {
        "significant_multiple_not_all"
      },
      stringsAsFactors = FALSE
    )
  }

  do.call(rbind, stability_rows)
}
