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
  "/Users/wangwanbing/Desktop/Chao Lab/Mediation Analysis/wwanbing/output/Chapter3/pipeline",
  fdr_tag
)
if (!dir.exists(output_dir)) {
  dir.create(output_dir, recursive = TRUE)
}

config <- list(
  fdr_cutoff = fdr_cutoff,
  scale_matrices = TRUE
)

chapter3 <- load_chapter3_data(input_dir)
message("Using FDR cutoff = ", config$fdr_cutoff, ".")
message("Using random seed = ", random_seed, ".")
message("Loaded ", nrow(chapter3$sample_data), " complete samples out of ", chapter3$original_n, ".")
if (chapter3$removed_n > 0) {
  message("Dropped ", chapter3$removed_n, " samples with missing values before modeling.")
}

all_results <- list()

for (method_name in c("survHIMA_fixed", "hima2_survival_fixed")) {
  method_dir <- file.path(output_dir, method_name)
  if (!dir.exists(method_dir)) {
    dir.create(method_dir, recursive = TRUE)
  }

  method_results <- list()

  for (snp_name in chapter3$snp_names) {
    message("Running ", method_name, " for ", snp_name, "...")

    pheno <- build_snp_pheno(
      sample_data = chapter3$sample_data,
      snp_name = snp_name,
      covariate_names = chapter3$covariate_names
    )

    result <- if (identical(method_name, "survHIMA_fixed")) {
      run_model_safe(
        function() {
          survHIMA_fixed(
            X = pheno$Exposure,
            Z = pheno[, chapter3$covariate_names, drop = FALSE],
            M = chapter3$mediator_matrix,
            OT = pheno$Time,
            status = pheno$Status,
            FDRcut = config$fdr_cutoff,
            scale = config$scale_matrices,
            verbose = TRUE
          )
        },
        snp_name = snp_name,
        method_name = method_name
      )
    } else {
      run_model_safe(
        function() {
          hima2_survival_fixed(
            survival::Surv(Time, Status) ~ Exposure + Age + Gender + Classical + Mesenchymal + Neural + Proneural,
            data.pheno = pheno,
            data.M = chapter3$mediator_matrix,
            FDRcut = config$fdr_cutoff,
            scale = config$scale_matrices,
            verbose = FALSE
          )
        },
        snp_name = snp_name,
        method_name = method_name
      )
    }

    annotated <- annotate_result(result, snp_name = snp_name, method_name = method_name)
    method_results[[snp_name]] <- annotated

    utils::write.csv(
      annotated,
      file = file.path(method_dir, paste0(snp_name, ".csv")),
      row.names = FALSE
    )
    save_indirect_effect_plot(
      annotated,
      out_file = file.path(method_dir, paste0(snp_name, "_ide.png")),
      title_text = paste(method_name, "-", snp_name)
    )
  }

  combined <- do.call(rbind, method_results)
  if (is.null(combined)) {
    combined <- annotate_result(empty_survhima_result(), snp_name = character(0), method_name = method_name)
  }
  utils::write.csv(
    combined,
    file = file.path(output_dir, paste0(method_name, "_combined.csv")),
    row.names = FALSE
  )
  all_results[[method_name]] <- combined
}

message("Done.")
message("Outputs written to: ", output_dir)
