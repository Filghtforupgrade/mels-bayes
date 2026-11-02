#!/usr/bin/env Rscript
# ============================================================================
# A2:  缺失数据报告 + 机制评估
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

cat("\n═══ A2: 缺失数据报告 + 机制评估 ═══\n")

# ---- 1. 逐变量缺失统计 ----
key_vars <- c("mood", "mood_lag1_day_c", "mood_prev_session_c",
              "health_session_c", "health_session_between_c",
              "unhealthy_session_c", "unhealthy_session_between_c",
              "balance_session_c",
              "study_day_c", "study_day_squared", "is_weekend", "session_numeric",
              "health_x_prev_mood_c")
l2_vars <- c("pre_phq9_c", "pre_control_c", "person_std_mood_c")

report_rows <- list()
for (sample_name in c("sample_B", "sample_A4", "sample_lagged")) {
  df <- get(sample_name)
  for (v in c(key_vars, l2_vars)) {
    if (!(v %in% names(df))) next
    n_total <- nrow(df)
    n_miss <- sum(is.na(df[[v]]))
    report_rows[[length(report_rows) + 1]] <- data.frame(
      sample = sample_name, variable = v,
      level = ifelse(v %in% l2_vars, "L2_person", "L1_observation"),
      n_total = n_total, n_missing = n_miss,
      pct_missing = round(100 * n_miss / n_total, 2)
    )
  }
}
report <- do.call(rbind, report_rows)
rownames(report) <- NULL
write.csv(report, file.path(output_dir, "missingness_report.csv"), row.names = FALSE)

cat("\n--- 缺失统计汇总 ---\n")
print(report[report$pct_missing > 0, ], row.names = FALSE)

# L2 缺失影响
n_B  <- nlevels(sample_B$user_id)
n_A4 <- nlevels(sample_A4$user_id)
n_L  <- nlevels(sample_lagged$user_id)
cat(sprintf("\n--- L2 缺失导致的人员丢弃 ---\n"))
cat(sprintf("  sample_B (完整):     %d 人\n", n_B))
cat(sprintf("  sample_A4 (有基线):  %d 人\n", n_A4))
cat(sprintf("  sample_lagged:       %d 人\n", n_L))
cat(sprintf("  因基线缺失丢弃:      %d 人 (%.1f%%)\n", n_B - n_A4, 100*(n_B-n_A4)/n_B))

# ---- 2. 缺失机制评估 ----
cat("\n--- 缺失机制评估 ---\n")
suppressPackageStartupMessages(library(naniar))

# Person-level dataframe
person_df <- sample_B %>%
  group_by(user_id) %>%
  summarise(
    person_mean_mood = mean(mood, na.rm = TRUE),
    person_mean_health = mean(health_session_c, na.rm = TRUE),
    person_std_mood = sd(mood, na.rm = TRUE),
    n_sessions = n(),
    n_days = n_distinct(study_day_c),
    .groups = "drop"
  ) %>%
  left_join(
    sample_A4 %>% distinct(user_id, pre_phq9_c, pre_control_c),
    by = "user_id"
  )

# Little's MCAR test
mcar_test_vars <- c("pre_phq9_c", "pre_control_c",
                     "person_mean_mood", "person_mean_health",
                     "person_std_mood", "n_sessions", "n_days")
mcar_df <- person_df[, mcar_test_vars]
mcar_result <- tryCatch(
  naniar::mcar_test(mcar_df),
  error = function(e) list(statistic = NA, df = NA, p.value = NA)
)

cat(sprintf("  Little's MCAR test: χ² = %s, df = %s, p = %s\n",
            if(is.na(mcar_result$statistic)) "NA" else sprintf("%.2f", mcar_result$statistic),
            if(is.na(mcar_result$df)) "NA" else as.character(mcar_result$df),
            if(is.na(mcar_result$p.value)) "NA" else sprintf("%.4f", mcar_result$p.value)))

mc_conclusion <- if (!is.na(mcar_result$p.value) && mcar_result$p.value > 0.05) {
  "MCAR not rejected — MI assumptions are strongest"
} else {
  "MCAR rejected or NA — assuming MAR; discuss MNAR in Limitations"
}
cat(sprintf("  Conclusion: %s\n", mc_conclusion))

# 缺失指示变量 × 观测变量的相关
miss_indicators <- data.frame(
  miss_phq = as.integer(is.na(person_df$pre_phq9_c)),
  miss_pc  = as.integer(is.na(person_df$pre_control_c))
)
obs_vars <- person_df[, c("person_mean_mood", "person_mean_health",
                           "person_std_mood", "n_sessions", "n_days")]
cor_matrix <- cor(cbind(miss_indicators, obs_vars), use = "pairwise.complete.obs")
cat("\n  缺失指示 × 观测变量的相关矩阵:\n")
print(round(cor_matrix, 3))

# 输入机制报告
mechanism_report <- data.frame(
  test = "Little_MCAR",
  statistic = if(is.na(mcar_result$statistic)) NA else mcar_result$statistic,
  df = if(is.na(mcar_result$df)) NA else mcar_result$df,
  p_value = if(is.na(mcar_result$p.value)) NA else mcar_result$p.value,
  conclusion = mc_conclusion
)
write.csv(mechanism_report,
          file.path(output_dir, "missingness_mechanism_assessment.csv"),
          row.names = FALSE)

cat(sprintf("\n  ✓ missingness_report.csv\n"))
cat(sprintf("  ✓ missingness_mechanism_assessment.csv\n"))
cat("═══ A2 完成 ═══\n")
