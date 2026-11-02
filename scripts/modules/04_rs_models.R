# ============================================================================
# 模块 04b (): 贝叶斯随机斜率模型 — 4 个 brms 拟合
# ============================================================================
# 依赖: 00b/01b/02b
# 产出: 4 个 brmsfit 对象 + 4 个 posterior CSV
#       person-specific slopes（B1_RS）→ 06b 进一步分析
# ============================================================================
# 模型清单说明：unhealthy/balance 与 health 的共线性在随机斜率模型中同样
#   存在（敏感度互为镜像/线性变换，与 B1_RS 不构成独立证据），仅拟合 B1 系列。
# 设计决策（计划文件第 2.04 节）:
#   - lme4 A1 boundary failure (gradient=0.009) 的贝叶斯救援
#   - 弱信息先验对 sd 设 student_t(3,0,2.5) → 不会被钉到 0
#   - 跨层交互模型用 adapt_delta=0.99, max_treedepth=15（更高难度）
#   - 4 个模型：B1_RS, B1_RS×PHQ, B1_RS×PC, B1_RS×σ^M
#     · 主调节项跑在 sample_A4（pre_phq9_c, pre_control_c, person_std_mood_c
#       仅在 A4 完整无缺）
#     · B1_RS 跑在 sample_B
# ============================================================================

cat("\n--- 模块04b (): 贝叶斯随机斜率模型 (B1_RS + 3 调节版) ---\n")

# ============================================================================
# 4Ba. B1_RS — 健康行为随机斜率（headline rescue）
# ============================================================================
# 核心问题：lme4 说"健康行为对 mood 的斜率没有跨人差异"，
#           贝叶斯能否给出可信的斜率方差区间？
B1_RS_form <- brms::bf(
  mood ~ health_session_c + health_session_between_c +
         mood_prev_session_c + study_day_c + session_numeric +
         (1 + health_session_c | user_id)
)

brms_fits$B1_RS <- fit_or_load(
  model_name = "B1_RS",
  formula = B1_RS_form,
  data = sample_B,
  prior = priors_rs,
  adapt_delta = 0.95,
  iter = 4000L,
  warmup = 2000L,
  threads = 2L
)
posterior_summaries$B1_RS <- extract_posterior_summary(brms_fits$B1_RS)
write.csv(posterior_summaries$B1_RS,
          file.path(posterior_dir, "B1_RS_posterior.csv"),
          row.names = FALSE)

# ============================================================================
# 4Bb. B1_RS_PHQ — 健康行为 × PHQ-9 跨层调节
# ============================================================================
# A1 lme4 失败 → 无法估计；本模型直接给出后验中位数 + 95% CrI
B1_RS_PHQ_form <- brms::bf(
  mood ~ health_session_c * pre_phq9_c +
         health_session_between_c +
         mood_prev_session_c + study_day_c + session_numeric +
         (1 + health_session_c | user_id)
)

brms_fits$B1_RS_PHQ <- fit_or_load(
  model_name = "B1_RS_PHQ",
  formula = B1_RS_PHQ_form,
  data = sample_A4,
  prior = priors_rs,
  adapt_delta = 0.99,
  max_treedepth = 15L,
  iter = 4000L,
  warmup = 2000L,
  threads = 2L
)
posterior_summaries$B1_RS_PHQ <- extract_posterior_summary(brms_fits$B1_RS_PHQ)
write.csv(posterior_summaries$B1_RS_PHQ,
          file.path(posterior_dir, "B1_RS_PHQ_posterior.csv"),
          row.names = FALSE)

# ============================================================================
# 4Bc. B1_RS_PC — 健康行为 × 感知控制 跨层调节
# ============================================================================
B1_RS_PC_form <- brms::bf(
  mood ~ health_session_c * pre_control_c +
         health_session_between_c +
         mood_prev_session_c + study_day_c + session_numeric +
         (1 + health_session_c | user_id)
)

brms_fits$B1_RS_PC <- fit_or_load(
  model_name = "B1_RS_PC",
  formula = B1_RS_PC_form,
  data = sample_A4,
  prior = priors_rs,
  adapt_delta = 0.99,
  max_treedepth = 15L,
  iter = 4000L,
  warmup = 2000L,
  threads = 2L
)
posterior_summaries$B1_RS_PC <- extract_posterior_summary(brms_fits$B1_RS_PC)
write.csv(posterior_summaries$B1_RS_PC,
          file.path(posterior_dir, "B1_RS_PC_posterior.csv"),
          row.names = FALSE)

# ============================================================================
# 4Bd. B1_RS_sigmaM — 健康行为 × 心境波动性 跨层调节
# ============================================================================
# 实质问题：心境本就更不稳定的人，是否对健康行为更敏感（更陡的斜率）？
B1_RS_sigmaM_form <- brms::bf(
  mood ~ health_session_c * person_std_mood_c +
         health_session_between_c +
         mood_prev_session_c + study_day_c + session_numeric +
         (1 + health_session_c | user_id)
)

brms_fits$B1_RS_sigmaM <- fit_or_load(
  model_name = "B1_RS_sigmaM",
  formula = B1_RS_sigmaM_form,
  data = sample_A4,
  prior = priors_rs,
  adapt_delta = 0.99,
  max_treedepth = 15L,
  iter = 4000L,
  warmup = 2000L,
  threads = 2L
)
posterior_summaries$B1_RS_sigmaM <- extract_posterior_summary(brms_fits$B1_RS_sigmaM)
write.csv(posterior_summaries$B1_RS_sigmaM,
          file.path(posterior_dir, "B1_RS_sigmaM_posterior.csv"),
          row.names = FALSE)

# ============================================================================
# 4Be. B2_RS / B3_RS — 已在  model-spec revision 中移除
# ============================================================================
# 与 B2_MELS / B3_MELS 平行删除：unhealthy 与 balance 在 L1 与 health 共线性问题
# 在随机斜率模型中同样存在（"个体对 unhealthy 的敏感度"与"个体对 health 的敏感度"
# 是同一个体差异的镜像/线性变换）。
# 历史 posterior CSV 仍保留在 04_results/posterior_summaries/ 作为 audit。

# ============================================================================
# 4Bg. 提取 B1_RS 个人级斜率（rescue 的实质性证据）
# ============================================================================
# coef(fit)$user_id[, , "health_session_c"] → 250 × 4 矩阵 (Estimate/SE/q2.5/q97.5)
extract_person_slopes <- function(fit, slope_var = "health_session_c") {
  re <- brms::ranef(fit, summary = TRUE)
  if (is.null(re$user_id)) return(NULL)
  # ranef 返回三维数组 [user, summary, effect]
  arr <- re$user_id
  if (!(slope_var %in% dimnames(arr)[[3]])) return(NULL)

  # 全局固定效应斜率
  fix_med <- brms::fixef(brms_fits$B1_RS)[slope_var, "Estimate"]

  data.frame(
    user_id = dimnames(arr)[[1]],
    slope_deviation = arr[, "Estimate", slope_var],
    se = arr[, "Est.Error", slope_var],
    q2.5 = arr[, "Q2.5", slope_var],
    q97.5 = arr[, "Q97.5", slope_var],
    # 个体绝对斜率 = 固定效应 + 个体偏离
    person_slope = fix_med + arr[, "Estimate", slope_var],
    # CrI 是否完全大于 0（积极效应）
    credibly_positive = arr[, "Q2.5", slope_var] + fix_med > 0,
    # CrI 是否完全小于 0
    credibly_negative = arr[, "Q97.5", slope_var] + fix_med < 0
  )
}

person_slopes_B1 <- extract_person_slopes(brms_fits$B1_RS, "health_session_c")
write.csv(person_slopes_B1,
          file.path(posterior_dir, "B1_RS_person_slopes.csv"),
          row.names = FALSE)

n_credibly_pos <- sum(person_slopes_B1$credibly_positive, na.rm = TRUE)
n_credibly_neg <- sum(person_slopes_B1$credibly_negative, na.rm = TRUE)
n_total <- nrow(person_slopes_B1)
cat(sprintf("\n  B1_RS 个体斜率分布:\n"))
cat(sprintf("    可信为正 (95%% CrI > 0): %d / %d (%.1f%%)\n",
            n_credibly_pos, n_total, 100 * n_credibly_pos / n_total))
cat(sprintf("    可信为负 (95%% CrI < 0): %d / %d (%.1f%%)\n",
            n_credibly_neg, n_total, 100 * n_credibly_neg / n_total))
cat(sprintf("    个体斜率范围: [%.3f, %.3f]\n",
            min(person_slopes_B1$person_slope, na.rm = TRUE),
            max(person_slopes_B1$person_slope, na.rm = TRUE)))

cat("\n  ✓ 模块04b 4 个贝叶斯随机斜率模型完成（ 删除 B2_RS / B3_RS；4 个模型均用 threads=2L within-chain 并行）\n")
