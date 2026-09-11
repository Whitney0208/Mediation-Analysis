empty_survhima_result <- function() {
  data.frame(
    ID = character(0),
    alpha = numeric(0),
    alpha_se = numeric(0),
    beta = numeric(0),
    beta_se = numeric(0),
    p.joint = numeric(0),
    stringsAsFactors = FALSE
  )
}

load_chapter3_data <- function(input_dir, snp_cols = 14:20) {
  mediator_raw <- readRDS(file.path(input_dir, "LQD_transformed_pixel.RDS"))
  mediator_matrix <- t(mediator_raw)
  mediator_matrix <- as.matrix(mediator_matrix)
  colnames(mediator_matrix) <- paste0("Mediator", seq_len(ncol(mediator_matrix)))

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
    snp_names = names(cov)[snp_cols],
    covariate_names = c("Age", "Gender", "Classical", "Mesenchymal", "Neural", "Proneural"),
    removed_n = removed_n,
    original_n = nrow(sample_data)
  )
}

build_snp_pheno <- function(sample_data, snp_name, covariate_names) {
  data.frame(
    Time = sample_data$Time,
    Status = sample_data$Status,
    Exposure = sample_data[[snp_name]],
    sample_data[, covariate_names, drop = FALSE],
    stringsAsFactors = FALSE
  )
}

build_pooled_data <- function(sample_data, mediator_matrix, snp_names, covariate_names) {
  n_samples <- nrow(sample_data)
  n_snp <- length(snp_names)

  pooled_mediator <- do.call(rbind, replicate(n_snp, mediator_matrix, simplify = FALSE))
  pooled_time <- rep(sample_data$Time, times = n_snp)
  pooled_status <- rep(sample_data$Status, times = n_snp)
  pooled_exposure <- unlist(sample_data[, snp_names, drop = FALSE], use.names = FALSE)
  pooled_covariates <- sample_data[
    rep(seq_len(n_samples), times = n_snp),
    covariate_names,
    drop = FALSE
  ]
  pooled_block <- rep(snp_names, each = n_samples)

  rownames(pooled_mediator) <- paste0("PooledSample", seq_len(nrow(pooled_mediator)))
  rownames(pooled_covariates) <- paste0("PooledSample", seq_len(nrow(pooled_covariates)))

  list(
    Exposure = pooled_exposure,
    Time = pooled_time,
    Status = pooled_status,
    Covariates = pooled_covariates,
    Mediator = pooled_mediator,
    SNPBlock = pooled_block,
    n_samples = n_samples,
    n_snp = n_snp
  )
}

annotate_result <- function(result, snp_name, method_name) {
  if (nrow(result) == 0) {
    result$SNP <- character(0)
    result$Method <- character(0)
    result$ide <- numeric(0)
    return(result)
  }

  result$SNP <- snp_name
  result$Method <- method_name
  result$ide <- result$alpha * result$beta
  result[, c("SNP", "Method", "ID", "alpha", "alpha_se", "beta", "beta_se", "p.joint", "ide")]
}

annotate_pooled_result <- function(result, method_name) {
  if (nrow(result) == 0) {
    result$Model <- character(0)
    result$Method <- character(0)
    result$ide <- numeric(0)
    return(result)
  }

  result$Model <- "pooled"
  result$Method <- method_name
  result$ide <- result$alpha * result$beta
  result[, c("Model", "Method", "ID", "alpha", "alpha_se", "beta", "beta_se", "p.joint", "ide")]
}

save_indirect_effect_plot <- function(result, out_file, title_text) {
  grDevices::png(filename = out_file, width = 1600, height = 900, res = 150)
  on.exit(grDevices::dev.off(), add = TRUE)

  if (nrow(result) == 0) {
    graphics::plot.new()
    graphics::text(0.5, 0.5, "No mediators passed the threshold")
    graphics::title(main = title_text)
    return(invisible(NULL))
  }

  ide <- result$ide
  y_lim <- range(c(0, ide), finite = TRUE)
  if (diff(y_lim) == 0) {
    y_lim <- y_lim + c(-0.5, 0.5)
  } else {
    pad <- 0.1 * diff(y_lim)
    y_lim <- y_lim + c(-pad, pad)
  }

  graphics::barplot(
    ide,
    names.arg = result$ID,
    las = 2,
    main = title_text,
    xlab = "Significant Mediators",
    ylab = "Indirect Effect Estimates",
    ylim = y_lim,
    cex.names = 0.6
  )
}

run_model_safe <- function(fun, snp_name, method_name) {
  tryCatch(
    fun(),
    error = function(e) {
      warning(method_name, " failed for ", snp_name, ": ", conditionMessage(e))
      empty_survhima_result()
    }
  )
}
