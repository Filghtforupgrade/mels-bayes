#!/usr/bin/env Rscript
# A4: Bayes Factor — bridge_sampler 两步法 (Schad et al. 2022, Dora R L164-171)
options(width=200); options(warn=1)
get_project_root <- function() {
  tryCatch({
    p <- sys.frame(1)$ofile
    if (!is.null(p)) return(normalizePath(file.path(dirname(p), "..", "..")))
  }, error = function(e) {})
  tryCatch({
    args <- commandArgs(trailingOnly = FALSE)
    f <- sub("^--file=", "", args[grep("^--file=", args)])
    if (length(f) > 0 && file.exists(f)) return(normalizePath(file.path(dirname(f), "..", "..")))
  }, error = function(e) {})
  getwd()
}
project_root <- get_project_root()
SCRIPT_DIR    <- file.path(project_root, "scripts", "pipeline")
MODULE_DIR <- file.path(project_root, "scripts", "modules")
source(file.path(MODULE_DIR, "00_setup.R"))
source(file.path(MODULE_DIR, "01_data.R"))
source(file.path(MODULE_DIR, "02_priors.R"))

cat("\n═══ A4: Bayes Factor ═══\n"); t0 <- Sys.time()

# Load full model
B1_MELS <- readRDS(file.path(fits_dir, "B1_MELS.rds"))

# ---- Null-sigma: keep mu health, remove sigma health ----
B1_MELS_null_sigma_form <- brms::bf(
  mood  ~ health_session_c + health_session_between_c +
          mood_lag1_day_c + study_day_c + session_numeric + (1 | user_id),
  sigma ~ 1
)
# ---- Null-full: remove all health ----
B1_MELS_null_full_form <- brms::bf(
  mood  ~ 1 + mood_lag1_day_c + study_day_c + session_numeric + (1 | user_id),
  sigma ~ 1
)

# Simplified priors for null models (no sigma predictors)
priors_null <- c(
  brms::prior(normal(0, 1),       class = "b"),
  brms::prior(normal(3, 2),       class = "Intercept"),
  brms::prior(student_t(3, 0, 1), class = "sd"),
  brms::prior(student_t(3, 0, 1), class = "Intercept", dpar = "sigma")
)

cat("Fitting null-sigma...\n")
B1_MELS_null_sigma <- fit_or_load("B1_MELS_null_sigma", B1_MELS_null_sigma_form,
  sample_lagged, priors_null, iter=6000L, warmup=3000L, threads=2L, sample_prior="yes")

cat("Fitting null-full...\n")
B1_MELS_null_full <- fit_or_load("B1_MELS_null_full", B1_MELS_null_full_form,
  sample_lagged, priors_null, iter=6000L, warmup=3000L, threads=2L, sample_prior="yes")

# ---- bridge_sampler ----
cat("Computing bridge_sampler...\n")
bs_full       <- bridgesampling::bridge_sampler(B1_MELS, silent=TRUE)
bs_null_sigma <- bridgesampling::bridge_sampler(B1_MELS_null_sigma, silent=TRUE)
bs_null_full  <- bridgesampling::bridge_sampler(B1_MELS_null_full, silent=TRUE)

# 注意：bs$cv 不存在 → 抛 "differing number of rows" 异常。
# 正确接口：bridgesampling::error_measures(bs)$percentage 返回百分比字符串如 "1.23%"；
# 解析为数值用于可靠性门控（< 5% reliable）。同时 CSV schema 与 reliable 字段保持程序化生成，
# cv_note 由 error_measures() 直接计算，无需手工占位。
cv_pct <- function(bs) {
  em <- tryCatch(bridgesampling::error_measures(bs), error = function(e) NULL)
  if (is.null(em) || is.null(em$percentage)) return(NA_real_)
  as.numeric(sub("%", "", em$percentage)) / 100
}
cv_full       <- cv_pct(bs_full)
cv_null_sigma <- cv_pct(bs_null_sigma)
cv_null_full  <- cv_pct(bs_null_full)
cat(sprintf("Bridge sampler CV:\n  Full: %.2f%%\n  Null-sigma: %.2f%%\n  Null-full: %.2f%%\n",
            100*cv_full, 100*cv_null_sigma, 100*cv_null_full))

# Bayes Factors — 从 logml 直接计算，避免 bf()$bf 兼容问题
log_bf_sigma <- bs_full$logml - bs_null_sigma$logml
log_bf_full  <- bs_full$logml - bs_null_full$logml
log_bf_mu    <- bs_null_sigma$logml - bs_null_full$logml

bf_results <- data.frame(
  comparison = c("Full vs Null-sigma (sigma effect)",
                 "Full vs Null-full (overall health effect)",
                 "Null-sigma vs Null-full (mu effect only)"),
  BF10 = exp(c(log_bf_sigma, log_bf_full, log_bf_mu)),
  log_BF10 = c(log_bf_sigma, log_bf_full, log_bf_mu),
  cv_model1 = c(cv_full, cv_full, cv_null_sigma),
  cv_model2 = c(cv_null_sigma, cv_null_full, cv_null_full)
)
bf_results$reliable <- !is.na(bf_results$cv_model1) & !is.na(bf_results$cv_model2) &
  bf_results$cv_model1 < 0.05 & bf_results$cv_model2 < 0.05
print(bf_results, row.names=FALSE, digits=3)
write.csv(bf_results, file.path(output_dir, "bayes_factor.csv"), row.names=FALSE)

md_path <- file.path(tables_dir, "table_12_bayes_factors.md")
md_lines <- c(
  "# Table S12: Bayes Factors (bridge sampling)", "",
  "| Comparison | BF₁₀ | log(BF₁₀) | Reliable? |",
  "|:-----------|-----:|----------:|:---------:|"
)
for(i in seq_len(nrow(bf_results))) {
  r <- bf_results[i,]
  md_lines <- c(md_lines, sprintf("| %s | %.2f | %.2f | %s |",
              r$comparison, r$BF10, r$log_BF10,
              ifelse(r$reliable, "✓ yes", "⚠ cv > 5%")))
}
writeLines(md_lines, md_path)
elapsed <- difftime(Sys.time(), t0, units="mins")
cat(sprintf("═══ A4 完成 — %.1f min ═══\n", as.numeric(elapsed)))
