# ============================================================================
# 模块 07b (): 跨范式对照表 + MD 表格
# ============================================================================
# 依赖: 03b/04b/05b/06b
# 产出: 4 张 paper_tables/ 目录下的 MD 表
#       1) table_01_mels_summary.md           — MELS mu+sigma 后验 (B1/R3)
#       2) table_02_bayes_rs_summary.md       — RS 模型固定+随机效应
#       3) table_03_brms_vs_lme4_comparison.md — brms vs lme4 跨范式对照
#       4) table_04_item_mels_summary.md      — 条目级 MELS 后验 (4 specs × 8 items)
#          ↑ 由 03_dry_runs/regenerate_item_mels_artifacts.R 在 master pipeline
#            阶段 9 生成（依赖 run_item_mels_1.R 在阶段 8 完成的 4 个条目级拟合）。
#            本模块 (阶段 6) 仅生成 Tables 1-3。
# ============================================================================

cat("\n--- 模块07b (): MD 表格生成 + 跨范式对照 ---\n")

# ============================================================================
# 7Ba. 加载  lme4 拟合（用于对照）
# ============================================================================
#  输出在 results_dir 下，但其 lme4 拟合是 RDS 还是 list 待确认
# 备用方案：重新拟合  的 B1/B2/B3/R3 lme4 模型（用与  完全相同的公式）
cat("  重新拟合  lme4 基线模型（B1/B2/B3/R3）以保对照同步...\n")
suppressPackageStartupMessages({
  library(lme4)
  library(lmerTest)
})

shared_controls <- "mood_prev_session_c + study_day_c + session_numeric"

# 用 lmerTest::lmer 替代 lme4::lmer，summary() 返回 Satterthwaite p 值
# （供 07b Table 3 的 lme4_p 列使用；compare_to_lme4(include_p=TRUE) 才会读取）
lme4_fits <- list()
lme4_fits$B1 <- lmerTest::lmer(
  as.formula(paste("mood ~ health_session_c + health_session_between_c +",
                   shared_controls, "+ (1 | user_id)")),
  data = sample_B, REML = TRUE)
#  新增（新发现 η）: B1_MELS 用 sample_lagged + mood_lag1_day_c，
# 与 lme4_fits$B1（sample_B + mood_prev_session_c）的 AR 控制项不一致。
# 创建 B1_lagged 专用对照基线，确保 B1_MELS vs lme4 公平对照。
lme4_fits$B1_lagged <- lmerTest::lmer(
  mood ~ health_session_c + health_session_between_c +
         mood_lag1_day_c + study_day_c + session_numeric +
         (1 | user_id),
  data = sample_lagged, REML = TRUE)
lme4_fits$B2 <- lmerTest::lmer(
  as.formula(paste("mood ~ unhealthy_session_c + unhealthy_session_between_c +",
                   shared_controls, "+ (1 | user_id)")),
  data = sample_B, REML = TRUE)
lme4_fits$B3 <- lmerTest::lmer(
  as.formula(paste("mood ~ balance_session_c +",
                   shared_controls, "+ (1 | user_id)")),
  data = sample_B, REML = TRUE)
lme4_fits$R3 <- lmerTest::lmer(
  mood ~ mood_lag1_day_c + health_session_c + health_x_prev_mood_c +
         study_day_c + study_day_squared + is_weekend + session_numeric +
         (1 | user_id),
  data = sample_lagged, REML = TRUE)
cat("  ✓ 5 个 lme4 基线模型已拟合（B1 / B1_lagged / B2 / B3 / R3，使用 lmerTest::lmer 以提供 Satterthwaite p）\n")

# ============================================================================
# 7Bb. Table 1: MELS 后验摘要
# ============================================================================
write_mels_table <- function() {
  out_path <- file.path(tables_dir, "table_01_mels_summary.md")
  con <- file(out_path, "wt")
  on.exit(close(con), add = TRUE)

  cat("# Table 1 (): MELS 模型后验摘要\n\n", file = con)
  cat("每个 MELS 模型的 mu 与 sigma 子方程参数后验中位数 + 95% CrI。\n",
      file = con, append = TRUE)
  cat("**关键发现**：sigma 子方程的健康行为系数 — 若负且 95% CrI 不含 0，\n",
      file = con, append = TRUE)
  cat("说明健康行为不仅升高 mood，还压缩个人内 mood 波动。\n\n",
      file = con, append = TRUE)

  for (mname in c("B1_MELS", "R3_MELS")) {
    if (is.null(brms_fits[[mname]])) next
    cat(sprintf("\n## %s\n\n", mname), file = con, append = TRUE)
    s <- posterior_summaries[[mname]]
    # mu 子方程
    mu_rows <- s[grepl("^b_(?!sigma_)", s$variable, perl = TRUE), ]
    cat("\n**mu 子方程（位置）**:\n\n", file = con, append = TRUE)
    cat("| 参数 | Median | SD | 95% CrI | R-hat | ESS bulk |\n",
        file = con, append = TRUE)
    cat("|:-----|-------:|---:|:--------|------:|---------:|\n",
        file = con, append = TRUE)
    for (i in seq_len(nrow(mu_rows))) {
      cat(sprintf("| %s | %.3f | %.3f | [%.3f, %.3f] | %.3f | %.0f |\n",
                  sub("^b_", "", mu_rows$variable[i]),
                  mu_rows$median[i], mu_rows$sd[i],
                  mu_rows$q2.5[i], mu_rows$q97.5[i],
                  mu_rows$rhat[i], mu_rows$ess_bulk[i]),
          file = con, append = TRUE)
    }

    # sigma 子方程
    sigma_rows <- s[grepl("^b_sigma_", s$variable), ]
    if (nrow(sigma_rows) > 0) {
      cat("\n**sigma 子方程（尺度，log 链接）**:\n\n", file = con, append = TRUE)
      cat("| 参数 | Median | SD | 95% CrI | exp(Median) | exp(95% CrI) | R-hat | ESS bulk |\n",
          file = con, append = TRUE)
      cat("|:-----|-------:|---:|:--------|------------:|:-------------|------:|---------:|\n",
          file = con, append = TRUE)
      for (i in seq_len(nrow(sigma_rows))) {
        cat(sprintf("| %s | %.3f | %.3f | [%.3f, %.3f] | %.3f | [%.3f, %.3f] | %.3f | %.0f |\n",
                    sub("^b_sigma_", "", sigma_rows$variable[i]),
                    sigma_rows$median[i], sigma_rows$sd[i],
                    sigma_rows$q2.5[i], sigma_rows$q97.5[i],
                    exp(sigma_rows$median[i]),
                    exp(sigma_rows$q2.5[i]), exp(sigma_rows$q97.5[i]),
                    sigma_rows$rhat[i], sigma_rows$ess_bulk[i]),
            file = con, append = TRUE)
      }
      cat("\n*exp(Median) = 每单位预测因子对残差 SD 的乘性效应 *\n",
          file = con, append = TRUE)
    }
  }
  cat(sprintf("  ✓ Table 1 已保存: %s\n", out_path))
}
write_mels_table()

# ============================================================================
# 7Bc. Table 2: 贝叶斯随机斜率后验摘要
# ============================================================================
write_rs_table <- function() {
  out_path <- file.path(tables_dir, "table_02_bayes_rs_summary.md")
  con <- file(out_path, "wt")
  on.exit(close(con), add = TRUE)

  cat("# Table 2 (): 贝叶斯随机斜率模型后验摘要\n\n", file = con)
  cat("**核心发现**：lme4 A1 在边界单一性失败 (gradient=0.009)；\n",
      file = con, append = TRUE)
  cat("贝叶斯估计给出可信的 sd(slope) 与 cor(intercept,slope) 后验区间，\n",
      file = con, append = TRUE)
  cat("从而支持后续 person-specific slope 分析与跨层调节检验。\n\n",
      file = con, append = TRUE)

  rs_models <- c("B1_RS", "B1_RS_PHQ", "B1_RS_PC", "B1_RS_sigmaM")
  for (mname in rs_models) {
    if (is.null(brms_fits[[mname]])) next
    cat(sprintf("\n## %s\n\n", mname), file = con, append = TRUE)

    s <- posterior_summaries[[mname]]
    # 固定效应
    fix_rows <- s[grepl("^b_", s$variable), ]
    cat("**固定效应**:\n\n", file = con, append = TRUE)
    cat("| 参数 | Median | SD | 95% CrI | R-hat | ESS bulk |\n",
        file = con, append = TRUE)
    cat("|:-----|-------:|---:|:--------|------:|---------:|\n",
        file = con, append = TRUE)
    for (i in seq_len(nrow(fix_rows))) {
      cat(sprintf("| %s | %.3f | %.3f | [%.3f, %.3f] | %.3f | %.0f |\n",
                  sub("^b_", "", fix_rows$variable[i]),
                  fix_rows$median[i], fix_rows$sd[i],
                  fix_rows$q2.5[i], fix_rows$q97.5[i],
                  fix_rows$rhat[i], fix_rows$ess_bulk[i]),
          file = con, append = TRUE)
    }

    # 随机效应（sd_user_id__* 与 cor_user_id__*）
    re_rows <- s[grepl("^(sd_|cor_)", s$variable), ]
    if (nrow(re_rows) > 0) {
      cat("\n**随机效应**:\n\n", file = con, append = TRUE)
      cat("| 参数 | Median | SD | 95% CrI | R-hat | ESS bulk |\n",
          file = con, append = TRUE)
      cat("|:-----|-------:|---:|:--------|------:|---------:|\n",
          file = con, append = TRUE)
      for (i in seq_len(nrow(re_rows))) {
        cat(sprintf("| %s | %.3f | %.3f | [%.3f, %.3f] | %.3f | %.0f |\n",
                    re_rows$variable[i],
                    re_rows$median[i], re_rows$sd[i],
                    re_rows$q2.5[i], re_rows$q97.5[i],
                    re_rows$rhat[i], re_rows$ess_bulk[i]),
            file = con, append = TRUE)
      }
    }
  }

  # B1_RS 个体斜率分布摘要
  if (exists("person_slopes_B1") && nrow(person_slopes_B1) > 0) {
    cat("\n## B1_RS 个体斜率分布（rescue 实质性证据）\n\n",
        file = con, append = TRUE)
    n_total <- nrow(person_slopes_B1)
    n_pos <- sum(person_slopes_B1$credibly_positive, na.rm = TRUE)
    n_neg <- sum(person_slopes_B1$credibly_negative, na.rm = TRUE)
    cat(sprintf("- 个体总数: %d\n", n_total), file = con, append = TRUE)
    cat(sprintf("- 95%% CrI 完全 > 0: %d (%.1f%%)\n",
                n_pos, 100 * n_pos / n_total), file = con, append = TRUE)
    cat(sprintf("- 95%% CrI 完全 < 0: %d (%.1f%%)\n",
                n_neg, 100 * n_neg / n_total), file = con, append = TRUE)
    cat(sprintf("- 个体斜率范围: [%.3f, %.3f]\n",
                min(person_slopes_B1$person_slope, na.rm = TRUE),
                max(person_slopes_B1$person_slope, na.rm = TRUE)),
        file = con, append = TRUE)
    cat(sprintf("- 个体斜率 IQR: [%.3f, %.3f]\n",
                quantile(person_slopes_B1$person_slope, 0.25, na.rm = TRUE),
                quantile(person_slopes_B1$person_slope, 0.75, na.rm = TRUE)),
        file = con, append = TRUE)
  }
  cat(sprintf("  ✓ Table 2 已保存: %s\n", out_path))
}
write_rs_table()

# ============================================================================
# 7Bd. Table 3: brms vs lme4 跨范式对照（最关键）
# ============================================================================
write_comparison_table <- function() {
  out_path <- file.path(tables_dir, "table_03_brms_vs_lme4_comparison.md")
  con <- file(out_path, "wt")
  on.exit(close(con), add = TRUE)

  cat("# Table 3 (): brms 后验中位数 vs lme4 REML β 跨范式对照\n\n",
      file = con)
  cat("**对比协议（先验承诺）**：mu 方程的 brms posterior median 应在 lme4 REML β\n",
      file = con, append = TRUE)
  cat("的 ≈1 SE 内；若 |Δ|/SE_lme4 > 1.0，必须解释（先验影响 vs 随机效应收缩传导）。\n",
      file = con, append = TRUE)
  cat("**pp_same_sign**：后验中与 lme4 估计同号的 draw 比例。\n\n",
      file = con, append = TRUE)

  # 对应关系：brms_fit ↔ lme4_fit ↔ 共同 terms
  # B1_MELS 用 mood_lag1_day_c，对照需用 lme4_fits$B1_lagged
  # 而非 lme4_fits$B1（后者用 mood_prev_session_c，AR 项不一致）
  pairs <- list(
    list(brms = "B1_MELS", lme4 = "B1_lagged",
         terms = c("health_session_c", "health_session_between_c",
                   "mood_lag1_day_c", "study_day_c", "session_numeric")),
    list(brms = "R3_MELS", lme4 = "R3",
         terms = c("mood_lag1_day_c", "health_session_c", "health_x_prev_mood_c",
                   "study_day_c", "study_day_squared", "is_weekend",
                   "session_numeric")),
    list(brms = "B1_RS", lme4 = "B1",
         terms = c("health_session_c", "health_session_between_c",
                   "mood_prev_session_c", "study_day_c", "session_numeric"))
  )

  all_comparison <- list()
  for (pr in pairs) {
    if (is.null(brms_fits[[pr$brms]]) || is.null(lme4_fits[[pr$lme4]])) next
    cat(sprintf("\n## %s (brms) vs %s (lme4)\n\n", pr$brms, pr$lme4),
        file = con, append = TRUE)
    cmp <- compare_to_lme4(brms_fits[[pr$brms]], lme4_fits[[pr$lme4]],
                            pr$terms, include_p = TRUE)
    if (nrow(cmp) == 0) {
      cat("  *(无可对照参数)*\n", file = con, append = TRUE)
      next
    }
    cmp$brms_model <- pr$brms
    cmp$lme4_model <- pr$lme4
    all_comparison[[paste0(pr$brms, "_vs_", pr$lme4)]] <- cmp

    cat("| Term | lme4 β | lme4 SE | lme4 p | brms median | brms 95% CrI | Δ | Δ/SE | pp_same_sign |\n",
        file = con, append = TRUE)
    cat("|:-----|-------:|--------:|-------:|------------:|:-------------|--:|-----:|-------------:|\n",
        file = con, append = TRUE)
    for (i in seq_len(nrow(cmp))) {
      flag <- if (abs(cmp$delta_over_se[i]) > 1.0) " ⚠" else ""
      cat(sprintf("| %s | %.3f | %.3f | %.3f | %.3f | [%.3f, %.3f] | %.3f | %.2f%s | %.3f |\n",
                  cmp$term[i], cmp$lme4_beta[i], cmp$lme4_se[i],
                  cmp$lme4_p[i], cmp$brms_median[i],
                  cmp$brms_q2.5[i], cmp$brms_q97.5[i],
                  cmp$delta[i], cmp$delta_over_se[i], flag,
                  cmp$pp_same_sign[i]),
          file = con, append = TRUE)
    }
  }

  cat(sprintf("  ✓ Table 3 已保存: %s\n", out_path))

  # 写出 CSV 副本
  if (length(all_comparison) > 0) {
    cmb <- do.call(rbind, all_comparison)
    rownames(cmb) <- NULL
    write.csv(cmb, file.path(output_dir, "comparison_brms_vs_lme4.csv"),
              row.names = FALSE)
    comparison_brms_vs_lme4 <<- all_comparison

    # 触发警告：任何 |Δ|/SE > 1.0
    n_warn_terms <- sum(abs(cmb$delta_over_se) > 1.0, na.rm = TRUE)
    if (n_warn_terms > 0) {
      cat(sprintf("\n  ⚠ 跨范式比较门：%d 个 term 的 |Δ|/SE > 1.0，需人工复核\n",
                  n_warn_terms))
      print(cmb[abs(cmb$delta_over_se) > 1.0,
                c("brms_model", "term", "lme4_beta", "brms_median",
                  "delta_over_se")],
            row.names = FALSE, digits = 3)
    } else {
      cat("\n  ✓ 跨范式比较门通过：所有 |Δ|/SE ≤ 1.0\n")
    }
  }
}
write_comparison_table()

cat("\n  ✓ 模块07b 3 张 MD 表 + 对照 CSV 完成\n")
