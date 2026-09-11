required_pkgs <- c("HIMA", "survival", "qvalue", "HDMT", "glmnet", "ncvreg")
missing_pkgs <- required_pkgs[!vapply(required_pkgs, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing_pkgs) > 0) {
  stop("Please install required package(s) first: ", paste(missing_pkgs, collapse = ", "))
}

invisible(lapply(required_pkgs, library, character.only = TRUE))

base_dir <- "/Users/wangwanbing/Desktop/Chao Lab/Mediation Analysis/wwanbing/code/Chapter3/pipeline"
fn_dir <- file.path(base_dir, "fn")
source(file.path(fn_dir, "03_pipeline_helpers.R"))
source(file.path(fn_dir, "01_ldpe_helpers.R"))
source(file.path(fn_dir, "02_survHIMA_fixed.R"))
source(file.path(fn_dir, "04_hima2_survival_fixed.R"))

input_dir <- "/Users/wangwanbing/Desktop/Chao Lab/Mediation Analysis/wwanbing/input_data/Chpater3"
random_seed <- as.integer(Sys.getenv("CHAPTER3_RANDOM_SEED", unset = "123"))
if (is.na(random_seed)) {
  stop("CHAPTER3_RANDOM_SEED must be an integer.")
}
set.seed(random_seed)

thresholds <- c(0.30, 0.35, 0.40, 0.45, 0.50, 0.60, 0.70, 0.80, 0.90, 1.00)
output_dir <- "/Users/wangwanbing/Desktop/Chao Lab/Mediation Analysis/wwanbing/output/Chapter3/pipeline/pooled_model/fdr_sensitivity"
if (!dir.exists(output_dir)) {
  dir.create(output_dir, recursive = TRUE)
}

chapter3 <- load_chapter3_data(input_dir)
pooled <- build_pooled_data(
  sample_data = chapter3$sample_data,
  mediator_matrix = chapter3$mediator_matrix,
  snp_names = chapter3$snp_names,
  covariate_names = chapter3$covariate_names
)

message("Running pooled FDR sensitivity analysis.")
message("Using random seed = ", random_seed, ".")
message("Pooling ", pooled$n_snp, " SNP blocks in this order: ", paste(chapter3$snp_names, collapse = ", "))

diag_res <- survHIMA_fixed_diagnostics(
  X = pooled$Exposure,
  Z = pooled$Covariates,
  M = pooled$Mediator,
  OT = pooled$Time,
  status = pooled$Status,
  scale = TRUE,
  verbose = TRUE
)

candidates <- diag_res$candidates
utils::write.csv(
  candidates,
  file = file.path(output_dir, "pooled_survHIMA_fixed_candidate_table.csv"),
  row.names = FALSE
)

summary_rows <- lapply(
  thresholds,
  function(thr) {
    data.frame(
      FDR = thr,
      n_hits = sum(candidates$fdr_value <= thr, na.rm = TRUE),
      stringsAsFactors = FALSE
    )
  }
)
summary_df <- do.call(rbind, summary_rows)
utils::write.csv(
  summary_df,
  file = file.path(output_dir, "pooled_survHIMA_fixed_fdr_summary.csv"),
  row.names = FALSE
)

first_threshold <- summary_df$FDR[summary_df$n_hits > 0][1]
if (is.na(first_threshold)) {
  message("No mediators were selected for any tested FDR threshold up to 1.00.")
} else {
  message("First tested FDR threshold with at least one mediator: ", first_threshold)
}
message("Candidate table written to: ", file.path(output_dir, "pooled_survHIMA_fixed_candidate_table.csv"))
message("Summary written to: ", file.path(output_dir, "pooled_survHIMA_fixed_fdr_summary.csv"))
message("Note: pooled hima2_survival_fixed shares the same survHIMA_fixed core in this pipeline, so this sensitivity summary applies to both pooled methods.")
