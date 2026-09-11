required_pkgs <- c("HIMA", "survival", "qvalue", "HDMT")
missing_pkgs <- required_pkgs[!vapply(required_pkgs, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing_pkgs) > 0) {
  stop("Please install required package(s) first: ", paste(missing_pkgs, collapse = ", "))
}

invisible(lapply(required_pkgs, library, character.only = TRUE))

# This script intentionally keeps the original Chapter 3 modeling logic
# and only adapts local input/output paths plus source() calls.
fn_dir <- "/Users/wangwanbing/Desktop/Chao Lab/Mediation Analysis/wwanbing/code/Chapter3/fn"
source(file.path(fn_dir, "1.R"))
source(file.path(fn_dir, "3.R"))
source(file.path(fn_dir, "7.R"))

input_dir <- "/Users/wangwanbing/Desktop/Chao Lab/Mediation Analysis/wwanbing/input_data/Chpater3"
output_dir <- "/Users/wangwanbing/Desktop/Chao Lab/Mediation Analysis/wwanbing/output/Chapter3"
if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)

r = readRDS(file.path(input_dir, "LQD_transformed_pixel.RDS"))
dim(r)
lqd = as.matrix(t(r))
dim(lqd)
rownames(lqd) = c(paste("S",1:64, sep=""))
colnames(lqd) = c(paste("P",1:1024, sep=""))
colnames(lqd) = c(paste("Mediator",1:1024, sep=""))

cov = read.csv(file.path(input_dir, "Q3_Covariates.csv"))
dim(cov)

snp_cols <- 14:20
snp_names <- names(cov)[snp_cols]
gender = ifelse(cov$Gender == "MALE",1,0)
covariates = data.frame(
  Age = cov$Age,
  Gender = gender,
  Classical = cov$Classical,
  Mesenchymal = cov$Mesenchymal,
  Neural = cov$Neural,
  Proneural = cov$Proneural
)
surv = cov$Survival.Overall
surv_status = cov[,7]
status = ifelse(surv_status == "DECEASED",T,F)
status = ifelse(surv_status == "DECEASED",0,1)
id = cov[,2]

sample_level_dat <- data.frame(
  Time = surv,
  Status = status,
  covariates,
  cov[, snp_cols, drop = FALSE]
)
complete_idx <- complete.cases(sample_level_dat) & complete.cases(lqd)
removed_n <- sum(!complete_idx)
if (removed_n > 0) {
  message("Dropping ", removed_n, " sample(s) with missing values before model fitting.")
}

lqd = lqd[complete_idx, , drop = FALSE]
sample_level_dat = sample_level_dat[complete_idx, , drop = FALSE]
n_samples <- nrow(sample_level_dat)
n_snp <- length(snp_names)

lqd2 = do.call(rbind, replicate(n_snp, lqd, simplify = FALSE))
rownames(lqd2) = paste("S", seq_len(nrow(lqd2)), sep = "")

snps = unlist(sample_level_dat[, snp_names, drop = FALSE], use.names = FALSE)
surv2 = rep(sample_level_dat$Time, times = n_snp)
status2 = rep(sample_level_dat$Status, times = n_snp)
cova = sample_level_dat[rep(seq_len(n_samples), times = n_snp),
                        c("Age", "Gender", "Classical", "Mesenchymal", "Neural", "Proneural"),
                        drop = FALSE]
d1 = data.frame(snps,status2,surv2,cova)
colnames(d1) = c("SNPs","Status","Time","Age","Gender","Classical","Mesenchymal","Neural","Proneural")

message("Running original survHIMA branch...")
sur = survHIMA(X = d1$SNPs,
               Z = d1[,c("Age","Gender","Classical","Mesenchymal","Neural","Proneural")],
               M = lqd2,
               OT = d1$Time,
               status = d1$Status,
               FDRcut = 0.35, verbose = TRUE, scale = FALSE)
ide = sur$alpha * sur$beta
result = cbind(sur,ide)
write.csv(result, file.path(output_dir, "chapter3_original_paths_only_survHIMA.csv"))
png(
  filename = file.path(output_dir, "chapter3_original_paths_only_survHIMA_ide.png"),
  width = 1600,
  height = 900,
  res = 150
)
if (nrow(result) == 0) {
  plot.new()
  text(0.5, 0.5, "No mediators passed the survHIMA threshold")
} else {
  barplot(ide, xlab = "Significant Mediators", ylab = "Indirect Effect Estimates", ylim = c(-1, 1))
}
dev.off()

message("Running original hima2 branch...")
sur = hima2(Surv(Status, Time) ~ SNPs + Gender + Age,
            data.pheno = d1,
            data.M = lqd2,
            outcome.family = "survival",
            mediator.family = "gaussian",
            penalty = "DBlasso",
            scale = TRUE,
            verbose = FALSE)
ide = sur$alpha * sur$beta
result = cbind(sur,ide)
write.csv(result, file.path(output_dir, "chapter3_original_paths_only_hima2.csv"))
png(
  filename = file.path(output_dir, "chapter3_original_paths_only_hima2_ide.png"),
  width = 1600,
  height = 900,
  res = 150
)
barplot(ide, xlab = "Significant Mediators", ylab = "Indirect Effect Estimates",ylim = c(-2,2))
dev.off()

message("Done.")
message("Outputs written to: ", output_dir)
