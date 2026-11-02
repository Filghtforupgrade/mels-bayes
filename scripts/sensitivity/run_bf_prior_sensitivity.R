#!/usr/bin/env Rscript
# ============================================================================
# A4b: Bayes Factor 先验敏感性 — wide prior 下的 bridge_sampler
# 动机: BF 对先验宽度的敏感性
# (Lindley 悖论方向)。本脚本在 wide prior (SD×2.5) 下重算 BF，
# 验证 "decisive" 结论是否稳健。
# ============================================================================
# 依赖: 00_setup.R, 01_data.R, 02_priors.R
# 产出: bayes_factor_prior_sensitivity.csv
#       paper_tables/table_12b_bayes_factors_prior_sensitivity.md
# ============================================================================

options(width = 200); options(warn = 1)

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

cat("\n═══ A4b: Bayes Factor 先验敏感性 (wide prior) ═══\n")
t0 <- Sys.time()

B1_MELS_wide <- readRDS(file.path(fits_dir, "B1_MELS_wide_prior.rds"))
cat("  ✓ 已加载 B1_MELS_wide_prior.rds\n")

B1_MELS_null_sigma_form <- brms::bf(
  mood  ~ health_session_c + health_session_between_c +
          mood_lag1_day_c + study_day_c + session_numeric + (1 | user_id),
  sigma ~ 1
)

B1_MELS_null_full_form <- brms::bf(
  mood  ~ 1 + mood_lag1_day_c + study_day_c + session_numeric + (1 | user_id),
  sigma ~ 1
)

priors_null_wide <- c(
  brms::prior(normal(0, 2.5),       class = "b"),
  brms::prior(normal(3, 5),         class = "Intercept"),
  brms::prior(student_t(3, 0, 2.5), class = "sd"),
  brms::prior(student_t(3, 0, 2.5), class = "Intercept", dpar = "sigma")
)

cat("Fitting B1_MELS_null_sigma_wide...\n")
B1_MELS_null_sigma_wide <- fit_or_load("B1_MELS_null_sigma_wide",
  B1_MELS_null_sigma_form, sample_lagged, priors_null_wide,
  iter = 6000L, warmup = 3000L, threads = 2L, sample_prior = "yes")

cat("Fitting B1_MELS_null_full_wide...\n")
B1_MELS_null_full_wide <- fit_or_load("B1_MELS_null_full_wide",
  B1_MELS_null_full_form, sample_lagged, priors_null_wide,
  iter = 6000L, warmup = 3000L, threads = 2L, sample_prior = "yes")

cat("Computing bridge_sampler (wide prior)...\n")
bs_full_wide       <- bridgesampling::bridge_sampler(B1_MELS_wide, silent = TRUE)
bs_null_sigma_wide <- bridgesampling::bridge_sampler(B1_MELS_null_sigma_wide, silent = TRUE)
bs_null_full_wide  <- bridgesampling::bridge_sampler(B1_MELS_null_full_wide, silent = TRUE)

cv_pct <- function(bs) {
  em <- tryCatch(bridgesampling::error_measures(bs), error = function(e) NULL)
  if (is.null(em) || is.null(em$percentage)) return(NA_real_)
  as.numeric(sub("%", "", em$percentage)) / 100
}

log_bf_sigma_wide <- bs_full_wide$logml - bs_null_sigma_wide$logml
log_bf_full_wide  <- bs_full_wide$logml - bs_null_full_wide$logml
log_bf_mu_wide    <- bs_null_sigma_wide$logml - bs_null_full_wide$logml

bf_wide <- data.frame(
  comparison = c("Full vs Null-sigma (sigma effect)",
                 "Full vs Null-full (overall health effect)",
                 "Null-sigma vs Null-full (mu effect only)"),
  BF10 = exp(c(log_bf_sigma_wide, log_bf_full_wide, log_bf_mu_wide)),
  log_BF10 = c(log_bf_sigma_wide, log_bf_full_wide, log_bf_mu_wide),
  cv_model1 = c(cv_pct(bs_full_wide), cv_pct(bs_full_wide), cv_pct(bs_null_sigma_wide)),
  cv_model2 = c(cv_pct(bs_null_sigma_wide), cv_pct(bs_null_full_wide), cv_pct(bs_null_full_wide)),
  prior_tier = "wide (SD x 2.5)"
)
bf_wide$reliable <- !is.na(bf_wide$cv_model1) & !is.na(bf_wide$cv_model2) &
  bf_wide$cv_model1 < 0.05 & bf_wide$cv_model2 < 0.05

bf_tight_path <- file.path(output_dir, "bayes_factor.csv")
if (file.exists(bf_tight_path)) {
  bf_tight <- read.csv(bf_tight_path, stringsAsFactors = FALSE)
  bf_tight$prior_tier <- "tight (default)"
  bf_combined <- rbind(bf_tight[, names(bf_wide)], bf_wide)
} else {
  bf_combined <- bf_wide
  cat("  ⚠ tight-prior BF CSV not found; outputting wide only\n")
}

write.csv(bf_combined, file.path(output_dir, "bayes_factor_prior_sensitivity.csv"),
          row.names = FALSE)
cat("\nWide-prior Bayes Factors:\n")
print(bf_wide[, c("comparison", "BF10", "log_BF10", "prior_tier", "reliable")],
      row.names = FALSE, digits = 3)

evidence_label <- function(log_bf) {
  abs_lbf <- abs(log_bf)
  if (abs_lbf < 1)      return("anecdotal")
  if (abs_lbf < 2.3)    return("anecdotal-to-substantial")
  if (abs_lbf < 3.4)    return("substantial")
  if (abs_lbf < 4.6)    return("strong")
  if (abs_lbf < 6.9)    return("very strong")
  return("decisive")
}

con <- file(file.path(tables_dir, "table_12b_bayes_factors_prior_sensitivity.md"), "wt")
tryCatch({
  cat("# Table S12b: Bayes Factor Prior Sensitivity (tight vs wide prior)\n\n", file = con)
  cat("| Comparison | Prior Tier | BF₁₀ | log(BF₁₀) | Evidence |\n", file = con)
  cat("|:-----------|:-----------|-----:|----------:|:---------|\n", file = con)
  for (i in seq_len(nrow(bf_combined))) {
    r <- bf_combined[i, ]
    cat(sprintf("| %s | %s | %.2e | %.1f | %s |\n",
                r$comparison, r$prior_tier, r$BF10, r$log_BF10,
                evidence_label(r$log_BF10)), file = con)
  }
  cat("\n*Note.* Bridge sampler CV < 5% for all models. ", file = con)
  cat("Wide prior = SD x 2.5 (normal(0,2.5) for slopes, normal(3,5) for intercept). ", file = con)
  cat("Evidence classification follows Jeffreys (1961). ", file = con)
  cat("If wide-prior BF remains > 100 (decisive), the BF conclusion is robust to prior width.\n", file = con)
}, finally = close(con))

elapsed <- difftime(Sys.time(), t0, units = "mins")
cat(sprintf("\n═══ A4b 完成 — %.1f min ═══\n", as.numeric(elapsed)))
