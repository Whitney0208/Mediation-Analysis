required_pkgs <- c("survival", "splines")
missing_pkgs <- required_pkgs[!vapply(required_pkgs, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing_pkgs) > 0) {
  stop("Please install required package(s) first: ", paste(missing_pkgs, collapse = ", "))
}

invisible(lapply(required_pkgs, library, character.only = TRUE))

base_dir <- "/Users/wangwanbing/Desktop/Chao Lab/Mediation Analysis/wwanbing/code/Chapter3/new_pipeline_Apr29"
fn_dir <- file.path(base_dir, "fn")
source(file.path(fn_dir, "01_data_helpers.R"))
source(file.path(fn_dir, "02_basis_helpers.R"))
source(file.path(fn_dir, "03_model_helpers.R"))
source(file.path(fn_dir, "04_run_helpers.R"))

input_dir <- "/Users/wangwanbing/Desktop/Chao Lab/Mediation Analysis/wwanbing/input_data/Chpater3"
output_dir <- "/Users/wangwanbing/Desktop/Chao Lab/Mediation Analysis/wwanbing/output/Apr29"
gene_names <- c("DDIT3", "EGFR", "KIT", "MDM4", "PDGFRA", "PIK3CA", "PTEN")
basis_k_values <- c(3L, 5L, 7L, 9L)
bootstrap_reps <- 100L
random_seed <- 123L
save_bootstrap_draws <- TRUE

run_root <- file.path(output_dir, "functional_spline_k_sensitivity")
combined_dir <- file.path(run_root, "combined")
dir.create(combined_dir, recursive = TRUE, showWarnings = FALSE)

message("Loading Chapter 3 functional mediation inputs...")
chapter3 <- load_chapter3_functional_data(input_dir)
gene_names <- intersect(gene_names, chapter3$gene_names)
if (length(gene_names) == 0) {
  stop("None of the requested gene names were found in the input data.")
}

message("Loaded ", nrow(chapter3$sample_data), " complete samples out of ", chapter3$original_n, ".")
if (chapter3$removed_n > 0) {
  message("Dropped ", chapter3$removed_n, " sample(s) with missing values.")
}
message("Mediator matrix dimensions: ", nrow(chapter3$mediator_matrix), " x ", ncol(chapter3$mediator_matrix))

sensitivity_rows <- list()
run_config_rows <- list()

for (basis_idx in seq_along(basis_k_values)) {
  basis_k <- basis_k_values[basis_idx]
  run_label <- paste0("k", basis_k)
  message("Running basis-k sensitivity for K = ", basis_k, "...")

  run_result <- run_functional_mediation_analysis(
    chapter3 = chapter3,
    gene_names = gene_names,
    basis_k = basis_k,
    bootstrap_reps = bootstrap_reps,
    random_seed = random_seed,
    save_bootstrap_draws = save_bootstrap_draws,
    run_root = file.path(run_root, run_label),
    output_label = run_label,
    write_combined_basis_summary = FALSE
  )

  gene_summary <- run_result$combined_gene_summary
  gene_summary$basis_k <- basis_k
  sensitivity_rows[[run_label]] <- gene_summary

  cfg <- run_result$run_config
  cfg$input_dir <- input_dir
  cfg$output_dir <- output_dir
  cfg$output_subfolder <- file.path("functional_spline_k_sensitivity", run_label)
  run_config_rows[[run_label]] <- cfg
}

gene_level_sensitivity_summary <- do.call(rbind, sensitivity_rows)
gene_level_sensitivity_summary <- gene_level_sensitivity_summary[
  order(gene_level_sensitivity_summary$Gene, gene_level_sensitivity_summary$basis_k),
  ,
  drop = FALSE
]

indirect_effect_stability <- summarize_indirect_effect_stability(gene_level_sensitivity_summary)
run_config_summary <- do.call(rbind, run_config_rows)

utils::write.csv(
  gene_level_sensitivity_summary,
  file = file.path(combined_dir, "gene_level_sensitivity_summary.csv"),
  row.names = FALSE
)
utils::write.csv(
  indirect_effect_stability,
  file = file.path(combined_dir, "indirect_effect_stability.csv"),
  row.names = FALSE
)
utils::write.csv(
  run_config_summary,
  file = file.path(combined_dir, "run_config.csv"),
  row.names = FALSE
)

message("Done.")
message("Outputs written to: ", run_root)
