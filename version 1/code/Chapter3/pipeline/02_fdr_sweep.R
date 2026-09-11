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

input_dir <- "/Users/wangwanbing/Desktop/Chao Lab/Mediation Analysis/wwanbing/input_data/Chpater3"
output_dir <- "/Users/wangwanbing/Desktop/Chao Lab/Mediation Analysis/wwanbing/output/Chapter3/pipeline"
if (!dir.exists(output_dir)) {
  dir.create(output_dir, recursive = TRUE)
}
random_seed <- as.integer(Sys.getenv("CHAPTER3_RANDOM_SEED", unset = "123"))
if (is.na(random_seed)) {
  stop("CHAPTER3_RANDOM_SEED must be an integer.")
}
set.seed(random_seed)

thresholds <- c(0.35, 0.40, 0.45, 0.50, 0.60, 0.70, 0.80, 0.90, 1.00)
chapter3 <- load_chapter3_data(input_dir)
sweep_rows <- list()
first_positive <- NULL

message("Using random seed = ", random_seed, ".")
for (thr in thresholds) {
  message("Evaluating FDR = ", thr, "...")
  total_hits <- 0L

  for (snp_name in chapter3$snp_names) {
    pheno <- build_snp_pheno(
      sample_data = chapter3$sample_data,
      snp_name = snp_name,
      covariate_names = chapter3$covariate_names
    )

    result <- survHIMA_fixed(
      X = pheno$Exposure,
      Z = pheno[, chapter3$covariate_names, drop = FALSE],
      M = chapter3$mediator_matrix,
      OT = pheno$Time,
      status = pheno$Status,
      FDRcut = thr,
      scale = TRUE,
      verbose = FALSE
    )

    hit_n <- nrow(result)
    total_hits <- total_hits + hit_n
    sweep_rows[[length(sweep_rows) + 1L]] <- data.frame(
      FDR = thr,
      SNP = snp_name,
      n_hits = hit_n,
      stringsAsFactors = FALSE
    )
  }

  if (is.null(first_positive) && total_hits > 0) {
    first_positive <- thr
  }
}

sweep_summary <- do.call(rbind, sweep_rows)
utils::write.csv(
  sweep_summary,
  file = file.path(output_dir, "survHIMA_fixed_fdr_sweep_summary.csv"),
  row.names = FALSE
)

if (is.null(first_positive)) {
  message("No mediators were selected for any tested FDR threshold up to 1.00.")
} else {
  message("First tested FDR threshold with at least one mediator: ", first_positive)
}
