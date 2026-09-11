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
basis_k <- 5L
bootstrap_reps <- 100L
random_seed <- 123L
save_bootstrap_draws <- TRUE

run_root <- file.path(output_dir, "functional_spline_k5_bootstrap")

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

run_result <- run_functional_mediation_analysis(
  chapter3 = chapter3,
  gene_names = gene_names,
  basis_k = basis_k,
  bootstrap_reps = bootstrap_reps,
  random_seed = random_seed,
  save_bootstrap_draws = save_bootstrap_draws,
  run_root = run_root,
  output_label = "functional_spline_k5_bootstrap",
  write_combined_basis_summary = TRUE
)

message("Done.")
message("Outputs written to: ", run_result$run_dirs$run_root)
