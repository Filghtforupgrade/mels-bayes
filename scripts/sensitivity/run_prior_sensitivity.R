#!/usr/bin/env Rscript
# ============================================================================
# A1:  先验敏感性分析 — 3 档先验 × 2 模型
# 对标 Dora et al. (2024) Example 2 的 weakly informative → flat 对比
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

cat("\n═══ A1: 先验敏感性分析 — 3 档先验 × 2 模型 ═══\n")
cat(sprintf("开始时间: %s\n", format(Sys.time(), "%Y-%m-%d %H:%M:%S")))
t0 <- Sys.time()

# ---- 2. 加载  tight fits ----
B1_MELS_tight  <- readRDS(file.path(fits_dir, "B1_MELS.rds"))
B1_RS_PHQ_tight <- readRDS(file.path(fits_dir, "B1_RS_PHQ.rds"))
cat("  ✓ 已加载 B1_MELS + B1_RS_PHQ tight fits\n")

# ---- 3. 定义 3 档先验 ----
# 档 1 (tight): 已有 priors_mels / priors_rs —  默认

# 档 2 (wide): SD 放大 2.5 倍
priors_mels_wide <- c(
  brms::prior(normal(0, 2.5),         class = "b"),
  brms::prior(normal(0, 2.5),         class = "b", dpar = "sigma"),
  brms::prior(normal(3, 5),           class = "Intercept"),
  brms::prior(student_t(3, 0, 2.5),   class = "Intercept", dpar = "sigma"),
  brms::prior(student_t(3, 0, 2.5),   class = "sd")
)
priors_rs_wide <- c(
  brms::prior(normal(0, 2.5),         class = "b"),
  brms::prior(normal(3, 5),           class = "Intercept"),
  brms::prior(student_t(3, 0, 2.5),   class = "sd"),
  brms::prior(student_t(3, 0, 6.25),   class = "sigma"),
  brms::prior(lkj(2),                 class = "cor")
)

# 档 3 (flat): 近似无信息（Dora Example 2 对标）— 仅对 MELS
priors_mels_flat <- c(
  brms::prior(normal(0, 10),          class = "b"),
  brms::prior(normal(0, 10),          class = "b", dpar = "sigma"),
  brms::prior(normal(3, 10),          class = "Intercept"),
  brms::prior(student_t(3, 0, 10),    class = "Intercept", dpar = "sigma"),
  brms::prior(student_t(3, 0, 10),    class = "sd")
)
cat("  ✓ 3 档先验已定义\n")

# ---- 4. 拟合 ----
B1_MELS_form <- brms::bf(
  mood  ~ health_session_c + health_session_between_c +
          mood_lag1_day_c + study_day_c + session_numeric +
          (1 | user_id),
  sigma ~ health_session_c + health_session_between_c
)

cat("\n--- 拟合 B1_MELS_wide_prior ---\n")
B1_MELS_wide <- fit_or_load("B1_MELS_wide_prior", B1_MELS_form,
  sample_lagged, priors_mels_wide, iter = 6000L, warmup = 3000L,
  threads = 2L, sample_prior = "yes")

cat("\n--- 拟合 B1_MELS_flat_prior ---\n")
B1_MELS_flat <- fit_or_load("B1_MELS_flat_prior", B1_MELS_form,
  sample_lagged, priors_mels_flat, iter = 6000L, warmup = 3000L,
  threads = 2L, sample_prior = "yes")

# B1_RS_PHQ: 定义公式（避免 source 04b 触发重跑）
B1_RS_PHQ_form <- brms::bf(
  mood ~ health_session_c * pre_phq9_c +
         health_session_between_c +
         mood_prev_session_c + study_day_c + session_numeric +
         (1 + health_session_c | user_id)
)

cat("\n--- 拟合 B1_RS_PHQ_wide_prior ---\n")
B1_RS_PHQ_wide <- fit_or_load("B1_RS_PHQ_wide_prior", B1_RS_PHQ_form,
  sample_A4, priors_rs_wide, iter = 4000L, warmup = 2000L,
  threads = 2L, sample_prior = "yes")

# ---- 5. 对比函数 ----
compare_priors <- function(fits_named, params) {
  result <- data.frame()
  for (p in params) {
    row <- data.frame(parameter = p)
    for (nm in names(fits_named)) {
      s <- extract_posterior_summary(fits_named[[nm]])
      r <- s[s$variable == p, ]
      if (nrow(r) == 0) next
      draws <- brms::as_draws_df(fits_named[[nm]])
      if (p %in% names(draws)) {
        # sigma 参数的"有意义"方向是负值（压缩波动）
        if (grepl("^b_sigma_", p)) {
          p_dir <- mean(draws[[p]] < 0, na.rm = TRUE)
        } else {
          p_dir <- mean(draws[[p]] > 0, na.rm = TRUE)
        }
      } else { p_dir <- NA }
      row[[paste0(nm, "_median")]] <- r$median
      row[[paste0(nm, "_q2.5")]]   <- r$q2.5
      row[[paste0(nm, "_q97.5")]]  <- r$q97.5
      row[[paste0(nm, "_rhat")]]   <- r$rhat
      row[[paste0(nm, "_P_dir")]]  <- p_dir
      row[[paste0(nm, "_excludes_zero")]] <- (r$q2.5 > 0 | r$q97.5 < 0)
    }
    result <- rbind(result, row)
  }
  # 计算结论翻转标志（tight vs wide）
  if (all(c("tight_excludes_zero", "wide_excludes_zero") %in% names(result))) {
    result$conclusion_flipped_tight_wide <- 
      result$tight_excludes_zero != result$wide_excludes_zero
  }
  if (all(c("tight_excludes_zero", "flat_excludes_zero") %in% names(result))) {
    result$conclusion_flipped_tight_flat <- 
      result$tight_excludes_zero != result$flat_excludes_zero
  }
  result
}

# ---- 6. 跑对比 ----
core_params_mels <- c(
  "b_health_session_c", "b_health_session_between_c",
  "b_mood_lag1_day_c", "b_study_day_c", "b_session_numeric",
  "b_sigma_health_session_c", "b_sigma_health_session_between_c",
  "sd_user_id__Intercept"
)

cat("\n--- 对比 B1_MELS: tight vs wide vs flat ---\n")
cmp_mels <- compare_priors(
  list(tight = B1_MELS_tight, wide = B1_MELS_wide, flat = B1_MELS_flat),
  core_params_mels)
print(cmp_mels[, c("parameter", "tight_median", "wide_median", "flat_median",
                    "tight_excludes_zero", "wide_excludes_zero", "flat_excludes_zero",
                    "conclusion_flipped_tight_wide", "conclusion_flipped_tight_flat")],
      row.names = FALSE, digits = 3)

cat("\n--- 对比 B1_RS_PHQ: tight vs wide ---\n")
cmp_rs <- compare_priors(
  list(tight = B1_RS_PHQ_tight, wide = B1_RS_PHQ_wide),
  c("b_health_session_c", "b_health_session_c:pre_phq9_c",
    "b_mood_prev_session_c", "sd_user_id__health_session_c"))
print(cmp_rs, row.names = FALSE, digits = 3)

# ---- 7. 写入 CSV ----
write.csv(cmp_mels, file.path(output_dir, "prior_sensitivity_B1_MELS.csv"), row.names = FALSE)
write.csv(cmp_rs,   file.path(output_dir, "prior_sensitivity_B1_RS_PHQ.csv"), row.names = FALSE)
cat("\n  ✓ CSV 已保存\n")

# ---- 8. 森林图（3 色：tight=steelblue, wide=darkorange, flat=grey40）----
library(ggplot2)

make_forest <- function(cmp_df, file_name, title_str) {
  # 重构为 long format
  prefix_map <- list()
  for (nm in c("tight", "wide", "flat")) {
    if (paste0(nm, "_median") %in% names(cmp_df)) {
      prefix_map[[nm]] <- nm
    }
  }
  if (length(prefix_map) < 2) return(invisible(NULL))
  
  long_rows <- list()
  for (nm in names(prefix_map)) {
    long_rows[[nm]] <- data.frame(
      parameter = cmp_df$parameter,
      prior = nm,
      median = cmp_df[[paste0(nm, "_median")]],
      q2.5   = cmp_df[[paste0(nm, "_q2.5")]],
      q97.5  = cmp_df[[paste0(nm, "_q97.5")]]
    )
  }
  long <- do.call(rbind, long_rows)
  long$parameter <- factor(long$parameter, levels = rev(unique(long$parameter)))
  
  color_map <- c(tight = "steelblue", wide = "darkorange", flat = "grey40")
  colors_use <- color_map[names(prefix_map)]
  
  p <- ggplot(long, aes(x = median, y = parameter, color = prior)) +
    geom_vline(xintercept = 0, linetype = "dashed", color = "grey50") +
    geom_errorbarh(aes(xmin = q2.5, xmax = q97.5),
                   position = position_dodge(width = 0.5), height = 0.2) +
    geom_point(position = position_dodge(width = 0.5), size = 2) +
    scale_color_manual(values = colors_use) +
    labs(title = title_str,
         subtitle = paste0("Comparison across ", length(prefix_map), " prior tiers"),
         x = "Posterior median (95% CrI)", y = NULL, color = "Prior") +
    theme_minimal(base_size = 11) +
    theme(legend.position = "top")
  
  ggsave(file.path(figures_dir, file_name), plot = p, width = 10, height = 7, dpi = 150)
  cat(sprintf("  ✓ %s\n", file_name))
}

make_forest(cmp_mels, "fig_prior_sensitivity.png",
            "Prior Sensitivity: B1_MELS (tight vs wide vs flat)")
make_forest(cmp_rs, "fig_prior_sensitivity_B1_RS_PHQ.png",
            "Prior Sensitivity: B1_RS_PHQ (tight vs wide)")

# ---- 9. MD 表 ----
write_prior_sensitivity_md <- function(cmp_df, out_name, title) {
  out_path <- file.path(tables_dir, out_name)
  con <- file(out_path, "wt")
  tryCatch({
    cat(sprintf("# %s\n\n", title), file = con)
    cat("Posterior median and 95% CrI for each prior tier.\n", file = con, append = TRUE)
    cat("**Prior tiers**: tight =  default; wide = SDs × 2.5; flat ≈ uninformative (normal(0,10)).\n\n", file = con, append = TRUE)
    cat("| Parameter | Tight median [95% CrI] | Wide median [95% CrI] | Flat median [95% CrI] | Conclusion flipped? |\n", file = con, append = TRUE)
    cat("|:----------|:-----------------------|:----------------------|:----------------------|:-------------------:|\n", file = con, append = TRUE)

    for (i in seq_len(nrow(cmp_df))) {
      r <- cmp_df[i, ]
      tight_str <- sprintf("%.3f [%.3f, %.3f]", r$tight_median, r$tight_q2.5, r$tight_q97.5)
      wide_str  <- if ("wide_median" %in% names(r) && !is.na(r$wide_median))
        sprintf("%.3f [%.3f, %.3f]", r$wide_median, r$wide_q2.5, r$wide_q97.5) else "—"
      flat_str  <- if ("flat_median" %in% names(r) && !is.na(r$flat_median))
        sprintf("%.3f [%.3f, %.3f]", r$flat_median, r$flat_q2.5, r$flat_q97.5) else "—"
      flipped   <- if ("conclusion_flipped_tight_wide" %in% names(r) && !is.na(r$conclusion_flipped_tight_wide))
        ifelse(r$conclusion_flipped_tight_wide, "⚠ YES", "no") else "—"
      cat(sprintf("| %s | %s | %s | %s | %s |\n",
                  r$parameter, tight_str, wide_str, flat_str, flipped),
          file = con, append = TRUE)
    }
  }, finally = close(con))
  cat(sprintf("  ✓ %s\n", out_name))
}

write_prior_sensitivity_md(cmp_mels, "table_11_prior_sensitivity.md",
                           "Table S11: Prior Sensitivity — B1_MELS")

# ---- 10. 总结 ----
elapsed <- difftime(Sys.time(), t0, units = "mins")
n_flipped_mels <- sum(cmp_mels$conclusion_flipped_tight_wide, na.rm = TRUE)

cat(sprintf("\n═══ A1 完成 — 耗时 %.1f 分钟 ═══\n", as.numeric(elapsed)))
cat(sprintf("  B1_MELS: %d 参数, %d 个结论翻转 (tight→wide)\n",
            nrow(cmp_mels), n_flipped_mels))
cat(sprintf("  B1_RS_PHQ: %d 参数, %d 个结论翻转\n",
            nrow(cmp_rs), sum(cmp_rs$conclusion_flipped_tight_wide, na.rm = TRUE)))

if (n_flipped_mels == 0) {
  cat("\n  ✓✓ 先验不影响实质性结论 — 核心发现对 prior 稳健\n")
} else {
  cat(sprintf("\n  ⚠ %d 个参数在 widen prior 后结论翻转，需在 Limitations 讨论\n", n_flipped_mels))
}
