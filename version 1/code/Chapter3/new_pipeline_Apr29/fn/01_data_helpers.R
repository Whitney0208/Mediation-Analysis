load_chapter3_functional_data <- function(input_dir, snp_cols = 14:20) {
  mediator_raw <- readRDS(file.path(input_dir, "LQD_transformed_pixel.RDS"))
  mediator_matrix <- t(mediator_raw)
  mediator_matrix <- as.matrix(mediator_matrix)
  mode(mediator_matrix) <- "numeric"
  colnames(mediator_matrix) <- sprintf("Grid%04d", seq_len(ncol(mediator_matrix)))

  cov <- read.csv(file.path(input_dir, "Q3_Covariates.csv"))

  status_text <- cov[["Corr..Overall.survival.status"]]
  status <- ifelse(
    status_text == "DECEASED",
    1L,
    ifelse(status_text == "LIVING", 0L, NA_integer_)
  )
  gender <- ifelse(
    cov$Gender == "MALE",
    1L,
    ifelse(cov$Gender == "FEMALE", 0L, NA_integer_)
  )

  sample_data <- data.frame(
    SampleID = cov$Case.ID,
    Time = cov$Survival.Overall,
    Status = status,
    Age = cov$Age,
    Gender = gender,
    Classical = cov$Classical,
    Mesenchymal = cov$Mesenchymal,
    Neural = cov$Neural,
    Proneural = cov$Proneural,
    cov[, snp_cols, drop = FALSE],
    stringsAsFactors = FALSE
  )

  complete_idx <- complete.cases(sample_data) & complete.cases(mediator_matrix)
  removed_n <- sum(!complete_idx)

  list(
    sample_data = sample_data[complete_idx, , drop = FALSE],
    mediator_matrix = mediator_matrix[complete_idx, , drop = FALSE],
    covariate_names = c("Age", "Gender", "Classical", "Mesenchymal", "Neural", "Proneural"),
    gene_names = names(cov)[snp_cols],
    grid = seq(0, 1, length.out = ncol(mediator_matrix)),
    removed_n = removed_n,
    original_n = nrow(sample_data)
  )
}
