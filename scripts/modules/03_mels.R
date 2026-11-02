# ============================================================================
# 模块 03b (): MELS 模型 — B1 / R3 location-scale 扩展
# ============================================================================
# 依赖: 00b/01b/02b
# 产出: 2 个 brmsfit 对象（存入 brms_fits 列表 + .rds 文件）
#       2 个 posterior_summary CSV
#       MELS sigma 子方程系数（mels_sigma_results）
# ============================================================================
# 模型清单说明：unhealthy_session_c 与 health_session_c 100% 共线（恒等式
#   health_noon + unhealthy_noon ≡ 4）；balance_session_c 与 health_session_c
#   r = +0.920（线性变换），均不携带独立信息，故本模块仅拟合 B1。
#   条目级分解见 scripts/validation/run_item_mels.R（8 个 CWC 行为条目的独立贡献）。
# 方案 B（定稿；方案 A 仍 FAIL: B1_MELS z_2[1,221] R̂=1.58 ESS=7）:
#   方案 A 移除 |p| 后, sigma 子方程的 (1|user_id) 仍难识别。
#   方案 B: 完全去掉 sigma 子方程的 RE,只保留固定效应。
#   失去: 个体异方差差异 (sd_user_id__sigma_Intercept) 估计
#   保留: MELS 核心 — 行为对残差 SD 的固定效应 b_sigma_health_session_c
#         (方案 A 中此参数 R̂=1.10 ESS=111, 已基本收敛)
# 方案 A：保留 mu 端 (1|user_id)，不足以解决 sigma 端问题
# 基线规格（保留）:
#   - AR 控制项: mood_prev_session_c → mood_lag1_day_c
#   - sigma 截距 / sd 先验已在 02b 收紧到 student_t(3,0,1)
# 设计决策保留:
#   - bf(mu ~ ..., sigma ~ ...) 双子方程
#   - R3 控制集匹配  R3 真实规格
# ============================================================================

cat("\n--- 模块03b (): MELS 模型 (B1 + R3) ---\n")

# ============================================================================
# 3Ba. B1_MELS — 健康行为 MELS
# ============================================================================
# 实质问题：当个体在某 session 比自己平时多做一项健康行为时，
#   (a) mood 是否升高（mu 子方程，与  B1 同问题）
#   (b) mood 残差方差是否压缩（sigma 子方程， 完全无法回答）
B1_MELS_form <- brms::bf(
  mood  ~ health_session_c + health_session_between_c +
          mood_lag1_day_c + study_day_c + session_numeric +
          (1 | user_id),
  sigma ~ health_session_c + health_session_between_c
)

brms_fits$B1_MELS <- fit_or_load(
  model_name = "B1_MELS",
  formula = B1_MELS_form,
  data = sample_lagged,
  prior = priors_mels,
  adapt_delta = 0.95,
  max_treedepth = 12L,
  iter = 6000L,
  warmup = 3000L,
  threads = 2L,
  sample_prior = "yes"
)
posterior_summaries$B1_MELS <- extract_posterior_summary(brms_fits$B1_MELS)
write.csv(posterior_summaries$B1_MELS,
          file.path(posterior_dir, "B1_MELS_posterior.csv"),
          row.names = FALSE)

# ============================================================================
# 3Bb. B2_MELS / B3_MELS — 已在  model-spec revision 中移除
# ============================================================================
# 历史 posterior CSV 仍保留在 04_results/posterior_summaries/ 作为 audit。
# 替代分析（条目级 MELS）：03_dry_runs/run_item_mels_1.R

# ============================================================================
# 3Bd. R3_MELS — 健康行为 × 前次 mood 缓冲效应 MELS
# ============================================================================
# 注意：R3 控制集 =  R3 真实规格（不是 mood_prev_session_c）
#   mood_lag1_day_c + study_day_c + study_day_squared + is_weekend + session_numeric
# 用 sample_lagged 数据集（删除首次测量）
R3_MELS_form <- brms::bf(
  mood  ~ mood_lag1_day_c + health_session_c + health_x_prev_mood_c +
          study_day_c + study_day_squared + is_weekend + session_numeric +
          (1 | user_id),
  sigma ~ health_session_c + health_x_prev_mood_c
)

brms_fits$R3_MELS <- fit_or_load(
  model_name = "R3_MELS",
  formula = R3_MELS_form,
  data = sample_lagged,
  prior = priors_mels,
  adapt_delta = 0.95,
  max_treedepth = 12L,
  iter = 6000L,
  warmup = 3000L,
  threads = 2L,
  sample_prior = "yes"
)
posterior_summaries$R3_MELS <- extract_posterior_summary(brms_fits$R3_MELS)
write.csv(posterior_summaries$R3_MELS,
          file.path(posterior_dir, "R3_MELS_posterior.csv"),
          row.names = FALSE)

# ============================================================================
# 3Be. 提取 sigma 子方程的关键系数（MELS 实质性发现）
# ============================================================================
# 把 b_sigma_health_session_c 等参数单独抽出 → 用于 forest 图与报告主表
extract_sigma_coefs <- function(fit, model_name) {
  s <- extract_posterior_summary(fit)
  sigma_rows <- s[grepl("^b_sigma_", s$variable), ]
  if (nrow(sigma_rows) == 0) return(NULL)
  data.frame(
    model = model_name,
    parameter = sub("^b_sigma_", "", sigma_rows$variable),
    median = sigma_rows$median,
    sd = sigma_rows$sd,
    q2.5 = sigma_rows$q2.5,
    q97.5 = sigma_rows$q97.5,
    rhat = sigma_rows$rhat,
    ess_bulk = sigma_rows$ess_bulk,
    # 解释性指标：exp(median) = 每单位预测因子对残差 SD 的乘性效应
    multiplicative_effect_on_sd = exp(sigma_rows$median),
    multiplicative_q2.5 = exp(sigma_rows$q2.5),
    multiplicative_q97.5 = exp(sigma_rows$q97.5)
  )
}

mels_sigma_results$B1_MELS <- extract_sigma_coefs(brms_fits$B1_MELS, "B1_MELS")
mels_sigma_results$R3_MELS <- extract_sigma_coefs(brms_fits$R3_MELS, "R3_MELS")

mels_sigma_combined <- do.call(rbind, mels_sigma_results)
write.csv(mels_sigma_combined,
          file.path(posterior_dir, "mels_sigma_coefficients.csv"),
          row.names = FALSE)

cat("\n  MELS sigma 子方程系数（实质性发现）:\n")
print(mels_sigma_combined, row.names = FALSE, digits = 3)

cat("\n  ✓ 模块03b 2 个 MELS 模型完成（B1 + R3；B2/B3 已在  移除）\n")
