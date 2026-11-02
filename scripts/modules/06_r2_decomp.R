# ============================================================================
# 模块 06b (): bayes_R2 + 五分量方差分解
# ============================================================================
# 依赖: 03b/04b/05b
# 产出: bayes_r2_decomposition_.csv
#       bayes_r2_total_.csv（每模型 bayes_R2 的后验中位数+CrI）
# ============================================================================
# 设计决策（计划文件第 2.06 节）:
#   - bayes_R2() 给总 R²
#   - 自定义 bayes_r2_decomp() 给五分量（与 r2mlm 平行）
#   - 仅对 PASS/WARN 模型运行（FAIL 模型方差分解无意义）
# ============================================================================
# var_map 口径：MELS 模型 = sample_lagged + mood_lag1_day_c。
#   后果：bayes_r2_decomp() 的 fname_short %in% within_vars 匹配失败，
#         AR 控制项 b_mood_lag1_day_c 的方差贡献被遗漏，f1_within 被低估
#   修正：B1/B2/B3_MELS 的 data/within 列同步到 sample_lagged + mood_lag1_day_c
# ============================================================================

cat("\n--- 模块06b (): bayes_R2 + 五分量方差分解 ---\n")

# ============================================================================
# 6Ba. 各模型的 within / between 变量映射
# ============================================================================
# 与 r2mlm 一致的变量分类：L1 within / L2 between
# : B1_MELS / R3_MELS 的 data/within 列已对齐到 03b 的实际拟合规格
model_var_map <- list(
  B1_MELS = list(
    data = sample_lagged,
    within  = c("health_session_c", "mood_lag1_day_c",
                "study_day_c", "session_numeric"),
    between = c("health_session_between_c")
  ),
  R3_MELS = list(
    data = sample_lagged,
    within  = c("mood_lag1_day_c", "health_session_c", "health_x_prev_mood_c",
                "study_day_c", "study_day_squared", "is_weekend",
                "session_numeric"),
    between = character(0)
  ),
  B1_RS = list(
    data = sample_B,
    within  = c("health_session_c", "mood_prev_session_c",
                "study_day_c", "session_numeric"),
    between = c("health_session_between_c")
  ),
  B1_RS_PHQ = list(
    data = sample_A4,
    within  = c("health_session_c", "mood_prev_session_c",
                "study_day_c", "session_numeric",
                "health_session_c:pre_phq9_c"),
    between = c("health_session_between_c", "pre_phq9_c")
  ),
  B1_RS_PC = list(
    data = sample_A4,
    within  = c("health_session_c", "mood_prev_session_c",
                "study_day_c", "session_numeric",
                "health_session_c:pre_control_c"),
    between = c("health_session_between_c", "pre_control_c")
  ),
  B1_RS_sigmaM = list(
    data = sample_A4,
    within  = c("health_session_c", "mood_prev_session_c",
                "study_day_c", "session_numeric",
                "health_session_c:person_std_mood_c"),
    between = c("health_session_between_c", "person_std_mood_c")
  )
)

# ============================================================================
# 6Bb. 跑 bayes_R2 + 五分量分解
# ============================================================================
r2_total_rows <- list()
r2_decomp_rows <- list()

for (mname in names(brms_fits)) {
  # 跳过 FAIL 模型
  status_row <- convergence[[mname]]
  if (!is.null(status_row) && status_row$status == "FAIL") {
    cat(sprintf("  [%s] FAIL — 跳过方差分解\n", mname))
    next
  }

  cat(sprintf("\n  [%s] 计算 bayes_R2 + 五分量分解...\n", mname))
  fit <- brms_fits[[mname]]
  vmap <- model_var_map[[mname]]

  # === 总 bayes_R2 ===
  tryCatch({
    r2_draws <- brms::bayes_R2(fit, summary = FALSE)
    r2_total <- if (is.matrix(r2_draws)) {
      data.frame(
        model = mname,
        target = colnames(r2_draws),
        median = apply(r2_draws, 2, median),
        q2.5   = apply(r2_draws, 2, quantile, 0.025),
        q97.5  = apply(r2_draws, 2, quantile, 0.975)
      )
    } else {
      data.frame(
        model = mname,
        target = "mu",
        median = median(r2_draws),
        q2.5   = quantile(r2_draws, 0.025),
        q97.5  = quantile(r2_draws, 0.975)
      )
    }
    r2_total_rows[[mname]] <- r2_total
    cat(sprintf("    bayes_R2 (mu): median=%.3f, 95%% CrI=[%.3f, %.3f]\n",
                r2_total$median[1], r2_total$q2.5[1], r2_total$q97.5[1]))
  }, error = function(e) {
    cat(sprintf("    ⚠ bayes_R2 失败: %s\n", e$message))
  })

  # === 五分量分解 ===
  # 五分量分解的 m_random 仅含随机截距方差
  #   (tau_int^2)，不含随机斜率方差。对 (1 + slope | user_id) 模型，随机斜率对
  #   总方差的贡献被系统性漏算，分解不完整。决策：随机斜率模型
  #   只报告总 Bayes R²，不做五分量分解；仅随机截距模型 (B1_MELS/R3_MELS) 分解。
  vc_chk <- tryCatch(brms::VarCorr(fit, summary = TRUE), error = function(e) NULL)
  has_random_slope <- !is.null(vc_chk) && !is.null(vc_chk$user_id) &&
    nrow(vc_chk$user_id$sd) > 1L   # >1 行 = 截距之外还有随机斜率
  if (has_random_slope) {
    cat(sprintf("    [%s] 含随机斜率 — 跳过五分量分解（F3：仅报告总 Bayes R²）\n",
                mname))
    next
  }

  tryCatch({
    decomp <- bayes_r2_decomp(
      fit,
      data = vmap$data,
      within_vars = vmap$within,
      between_vars = vmap$between,
      n_draws = 500L  # 500 draws 已足够稳定，节省时间
    )
    decomp$model <- mname
    r2_decomp_rows[[mname]] <- decomp[, c("model", "component",
                                           "median", "q2.5", "q97.5")]
    bayes_r2_results[[mname]] <- decomp
    cat(sprintf("    五分量分解完成: f1=%.3f, f2=%.3f, m=%.3f, resid=%.3f, fvm=%.3f\n",
                decomp$median[decomp$component == "f1_within"],
                decomp$median[decomp$component == "f2_between"],
                decomp$median[decomp$component == "m_random"],
                decomp$median[decomp$component == "sigma2_resid"],
                decomp$median[decomp$component == "fvm_total"]))
  }, error = function(e) {
    cat(sprintf("    ⚠ 五分量分解失败: %s\n", e$message))
  })
}

# ============================================================================
# 6Bc. 保存汇总 CSV
# ============================================================================
if (length(r2_total_rows) > 0) {
  r2_total_combined <- do.call(rbind, r2_total_rows)
  rownames(r2_total_combined) <- NULL
  write.csv(r2_total_combined,
            file.path(output_dir, "bayes_r2_total_.csv"),
            row.names = FALSE)
  cat(sprintf("\n  ✓ bayes_R2 总值已保存\n"))
}

if (length(r2_decomp_rows) > 0) {
  r2_decomp_combined <- do.call(rbind, r2_decomp_rows)
  rownames(r2_decomp_combined) <- NULL
  write.csv(r2_decomp_combined,
            file.path(output_dir, "bayes_r2_decomposition_.csv"),
            row.names = FALSE)
  cat(sprintf("  ✓ 五分量分解已保存\n"))
  cat("\n  方差分解汇总（component × model）:\n")
  pivot_df <- tidyr::pivot_wider(
    r2_decomp_combined[, c("model", "component", "median")],
    names_from = "component", values_from = "median"
  )
  print(pivot_df, row.names = FALSE, digits = 3)
}

cat("\n  ✓ 模块06b 方差分解完成\n")
