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
fdr_cutoff <- as.numeric(Sys.getenv("CHAPTER3_FDR_CUTOFF", unset = "0.30"))
random_seed <- as.integer(Sys.getenv("CHAPTER3_RANDOM_SEED", unset = "123"))
if (!is.finite(fdr_cutoff) || fdr_cutoff <= 0 || fdr_cutoff > 1) {
  stop("CHAPTER3_FDR_CUTOFF must be a numeric value in (0, 1].")
}
if (is.na(random_seed)) {
  stop("CHAPTER3_RANDOM_SEED must be an integer.")
}
set.seed(random_seed)

fdr_tag <- gsub("\\.", "_", sprintf("fdr_%0.2f", fdr_cutoff))
output_dir <- file.path(
  "/Users/wangwanbing/Desktop/Chao Lab/Mediation Analysis/wwanbing/output/Chapter3/pipeline/pooled_model",
  fdr_tag
)
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

message("Running pooled SNP mediation model.")
message("Using FDR cutoff = ", fdr_cutoff, ".")
message("Using random seed = ", random_seed, ".")
message("Loaded ", nrow(chapter3$sample_data), " complete samples out of ", chapter3$original_n, ".")
if (chapter3$removed_n > 0) {
  message("Dropped ", chapter3$removed_n, " samples with missing values before modeling.")
}
message("Pooling ", pooled$n_snp, " SNP blocks in this order: ", paste(chapter3$snp_names, collapse = ", "))
message(
  "Pooled dimensions: Exposure=", length(pooled$Exposure),
  ", Time/Status=", length(pooled$Time),
  ", Covariates=", nrow(pooled$Covariates), "x", ncol(pooled$Covariates),
  ", Mediator=", nrow(pooled$Mediator), "x", ncol(pooled$Mediator)
)

results <- list()
for (method_name in c("survHIMA_fixed", "hima2_survival_fixed")) {
  message("Running pooled ", method_name, "...")

  result <- if (identical(method_name, "survHIMA_fixed")) {
    run_model_safe(
      function() {
        survHIMA_fixed(
          X = pooled$Exposure,
          Z = pooled$Covariates,
          M = pooled$Mediator,
          OT = pooled$Time,
          status = pooled$Status,
          FDRcut = fdr_cutoff,
          scale = TRUE,
          verbose = TRUE
        )
      },
      snp_name = "pooled",
      method_name = method_name
    )
  } else {
    pooled_pheno <- data.frame(
      Time = pooled$Time,
      Status = pooled$Status,
      Exposure = pooled$Exposure,
      pooled$Covariates,
      stringsAsFactors = FALSE
    )
    run_model_safe(
      function() {
        hima2_survival_fixed(
          survival::Surv(Time, Status) ~ Exposure + Age + Gender + Classical + Mesenchymal + Neural + Proneural,
          data.pheno = pooled_pheno,
          data.M = pooled$Mediator,
          FDRcut = fdr_cutoff,
          scale = TRUE,
          verbose = FALSE
        )
      },
      snp_name = "pooled",
      method_name = method_name
    )
  }

  annotated <- annotate_pooled_result(result, method_name = method_name)
  csv_file <- file.path(output_dir, paste0("pooled_", method_name, ".csv"))
  png_file <- file.path(output_dir, paste0("pooled_", method_name, "_ide.png"))

  utils::write.csv(annotated, file = csv_file, row.names = FALSE)
  save_indirect_effect_plot(
    annotated,
    out_file = png_file,
    title_text = paste("pooled", method_name)
  )

  message("Selected mediators for ", method_name, ": ", nrow(annotated))
  results[[method_name]] <- annotated
}

message("Done.")
message("Outputs written to: ", output_dir)
