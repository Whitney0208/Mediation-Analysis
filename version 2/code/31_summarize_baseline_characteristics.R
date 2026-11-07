#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  if (!requireNamespace("survival", quietly = TRUE)) {
    stop("Package 'survival' is required.")
  }
})

script_dir <- function() {
  args <- commandArgs(trailingOnly = FALSE)
  file_arg <- grep("^--file=", args, value = TRUE)
  if (length(file_arg) == 1L) {
    return(dirname(normalizePath(sub("^--file=", "", file_arg))))
  }
  normalizePath(getwd())
}

code_dir <- script_dir()
project_dir <- dirname(code_dir)
model_script <- file.path(code_dir, "02_fit_four_snp_mediation_models.R")
source(model_script, local = FALSE)

output_dir <- file.path(project_dir, "results", "baseline_characteristics")
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

cc_regions <- c(
  "CorpusCallosumAnterior", "CorpusCallosumMidAnterior", "CorpusCallosumCentral",
  "CorpusCallosumMidPosterior", "CorpusCallosumPosterior"
)
volume_table <- utils::read.csv(file.path(project_dir, "data", "SPLS-Coxdata", "CCinfo.csv"))
volume_rows <- match(as.integer(subject$RID), as.integer(volume_table$RID))
if (anyNA(volume_rows)) {
  stop("Corpus callosum volumes are missing for some analysis subjects.")
}
volume_table <- volume_table[volume_rows, , drop = FALSE]
cc_volume <- rowSums(volume_table[, cc_regions])
icv_cm3 <- volume_table$ICV / 1000

mean_sd <- function(values) {
  missing_count <- sum(is.na(values))
  text <- sprintf("%.1f (%.1f)", mean(values, na.rm = TRUE), stats::sd(values, na.rm = TRUE))
  if (missing_count > 0L) sprintf("%s [%d missing]", text, missing_count) else text
}
count_percent <- function(flag) {
  sprintf("%d (%.1f)", sum(flag, na.rm = TRUE), 100 * mean(flag, na.rm = TRUE))
}

groups <- list(
  Overall = rep(TRUE, nrow(subject)),
  Converted = subject$event == 1L,
  `Did not convert` = subject$event == 0L
)

continuous_rows <- list(
  `Age, years, mean (SD)` = subject$AGE,
  `Education, years, mean (SD)` = subject$PTEDUCAT,
  `MMSE, mean (SD)` = as.numeric(subject$MMSE),
  `ADAS-Cog 13, mean (SD)` = as.numeric(subject$ADAS13),
  `CDR-SB, mean (SD)` = as.numeric(subject$CDRSB),
  `FAQ, mean (SD)` = as.numeric(subject$FAQ),
  `Corpus callosum volume, mm3, mean (SD)` = cc_volume,
  `Intracranial volume, cm3, mean (SD)` = icv_cm3
)
categorical_rows <- list(
  `Male, n (%)` = subject$PTGENDER == "Male",
  `Right-handed, n (%)` = subject$PTHAND == 1,
  `White, n (%)` = subject$PTRACCAT == "White",
  `Not Hispanic or Latino, n (%)` = subject$PTETHCAT == "Not Hisp/Latino",
  `APOE4 alleles = 0, n (%)` = subject$APOE4 == 0,
  `APOE4 alleles = 1, n (%)` = subject$APOE4 == 1,
  `APOE4 alleles = 2, n (%)` = subject$APOE4 == 2
)

baseline_table <- do.call(rbind, c(
  list(data.frame(
    Characteristic = "Subjects, n",
    t(vapply(groups, function(group) as.character(sum(group)), character(1L))),
    check.names = FALSE
  )),
  lapply(names(continuous_rows), function(label) {
    data.frame(
      Characteristic = label,
      t(vapply(groups, function(group) mean_sd(continuous_rows[[label]][group]), character(1L))),
      check.names = FALSE
    )
  }),
  lapply(names(categorical_rows), function(label) {
    data.frame(
      Characteristic = label,
      t(vapply(groups, function(group) count_percent(categorical_rows[[label]][group]), character(1L))),
      check.names = FALSE
    )
  })
))
utils::write.csv(
  baseline_table,
  file.path(output_dir, "baseline_characteristics.csv"),
  row.names = FALSE
)

# Hardy-Weinberg exact test (Wigginton, Cutler and Abecasis, 2005).
hwe_exact_p <- function(genotype_counts) {
  hom_1 <- genotype_counts[[1L]]
  het <- genotype_counts[[2L]]
  hom_2 <- genotype_counts[[3L]]
  subjects <- hom_1 + het + hom_2
  rare_alleles <- 2L * min(hom_1, hom_2) + het
  het_values <- seq(rare_alleles %% 2L, rare_alleles, by = 2L)
  log_probability <- vapply(het_values, function(h) {
    rare_hom <- (rare_alleles - h) / 2
    common_hom <- subjects - h - rare_hom
    h * log(2) + lfactorial(subjects) - lfactorial(rare_hom) - lfactorial(h) -
      lfactorial(common_hom) + lfactorial(rare_alleles) +
      lfactorial(2L * subjects - rare_alleles) - lfactorial(2L * subjects)
  }, numeric(1L))
  probability <- exp(log_probability)
  probability <- probability / sum(probability)
  observed <- probability[match(het, het_values)]
  min(1, sum(probability[probability <= observed * (1 + 1e-9)]))
}

genotype_table <- do.call(rbind, lapply(snp_names, function(snp) {
  dosage <- X[, snp]
  counts <- as.integer(table(factor(dosage, levels = 0:2)))
  converted <- as.integer(table(factor(dosage[subject$event == 1L], levels = 0:2)))
  not_converted <- as.integer(table(factor(dosage[subject$event == 0L], levels = 0:2)))
  data.frame(
    SNP = snp,
    dosage_0 = counts[1L],
    dosage_1 = counts[2L],
    dosage_2 = counts[3L],
    coded_allele_frequency = mean(dosage) / 2,
    hwe_exact_p = hwe_exact_p(counts),
    converted_0_1_2 = paste(converted, collapse = "/"),
    not_converted_0_1_2 = paste(not_converted, collapse = "/")
  )
}))
utils::write.csv(
  genotype_table,
  file.path(output_dir, "snp_genotype_summary.csv"),
  row.names = FALSE
)

# Correlation of allele dosages; its square is the usual r2 measure of LD.
dosage_correlation <- stats::cor(X[, snp_names])
utils::write.csv(
  dosage_correlation^2,
  file.path(output_dir, "snp_dosage_r_squared.csv")
)
utils::write.csv(
  dosage_correlation,
  file.path(output_dir, "snp_dosage_correlation.csv")
)

survival_object <- survival::Surv(subject$time, subject$event)
conversion_fit <- survival::survfit(survival_object ~ 1)
follow_up_fit <- survival::survfit(survival::Surv(subject$time, 1L - subject$event) ~ 1)
landmarks <- c(365, 730, 1095)
conversion_at <- summary(conversion_fit, times = landmarks)
follow_up_table <- data.frame(
  quantity = c(
    "Observed time, median",
    "Observed time, 25th percentile",
    "Observed time, 75th percentile",
    "Observed time, minimum",
    "Observed time, maximum",
    "Median follow-up, reverse Kaplan-Meier",
    "Median time to conversion, Kaplan-Meier",
    "Median time to conversion, lower 95% limit",
    "Median time to conversion, upper 95% limit",
    sprintf("Cumulative conversion at time %d", landmarks)
  ),
  value = c(
    stats::median(subject$time),
    unname(stats::quantile(subject$time, 0.25)),
    unname(stats::quantile(subject$time, 0.75)),
    min(subject$time),
    max(subject$time),
    unname(summary(follow_up_fit)$table["median"]),
    unname(summary(conversion_fit)$table["median"]),
    unname(summary(conversion_fit)$table["0.95LCL"]),
    unname(summary(conversion_fit)$table["0.95UCL"]),
    1 - conversion_at$surv
  )
)
utils::write.csv(
  follow_up_table,
  file.path(output_dir, "follow_up_summary.csv"),
  row.names = FALSE
)

message("Completed. Results written to: ", output_dir)
print(baseline_table, row.names = FALSE)
print(genotype_table, row.names = FALSE, digits = 3)
print(round(dosage_correlation^2, 3))
print(follow_up_table, row.names = FALSE, digits = 4)
print(table(subject$PTHAND, W[, "handedness"], useNA = "ifany"))
