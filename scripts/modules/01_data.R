# ============================================================================
# 模块 01b (): 数据导入 — 复用  prep + 不变量断言
# ============================================================================
# 依赖: 00b_setup_brms.R (路径 + 容器)
# 产出: sample_B, sample_A4, sample_lagged（与  完全同源）
#       4 项不变量已验证； lme4 helper 已 source（用于 07b 对照）
# ============================================================================
# 设计决策（计划文件第 4 节）:
#   - 优先加载  已生成的 checkpoint，避免重跑 350 行数据准备
#   - 缺失则触发  00_config + 01_data_prep（自动写检查点）
#   - 4 项不变量：人数=250, 行数>8000, 关键变量齐全, person_std_mood_c 在 A4
# ============================================================================

cat("\n--- 模块01b (): 数据导入与不变量断言 ---\n")

# ============================================================================
# 1Ba. 优先 checkpoint 加载，否则触发  prep
# ============================================================================
ckpt_path <- file.path(checkpoint_dir, "checkpoint_01_data.RData")

if (file.exists(ckpt_path)) {
  cat(sprintf("  从 checkpoint 加载: %s\n", basename(ckpt_path)))
  load(ckpt_path)
  cat("    ✓ session_df / sample_B / sample_A4 / sample_lagged 已加载\n")
} else {
  cat("  checkpoint 不存在，触发  数据准备...\n")
  #  00_config 依赖一些脚本路径变量；为避免覆盖  自身路径，
  # 临时 source  模块到独立环境后只导出 data frame
  source(file.path(module_dir, "00_config.R"), local = TRUE)
  source(file.path(module_dir, "01_data_prep.R"))
  cat("    ✓  prep 完成，checkpoint 已写入\n")
}

# 始终 source  helper（用于 07b 模块的 lme4 对照工具）
source(file.path(module_dir, "02_helper_funcs.R"), local = TRUE)
cat("  ✓  02_helper_funcs.R 已 source（仅供对照）\n")

# ============================================================================
# 1Ba2. L2 调节变量标准化
# ============================================================================
#  data_prep 只对 pre_phq9_c / pre_control_c 做了中心化（未除以 SD）：
#   实测 SD(pre_phq9_c)=4.34, SD(pre_control_c)=9.68，而 person_std_mood_c SD≈1.14。
# 后果：(1) 三个跨层交互系数尺度不可比；(2) SESOI 统一阈值 0.02/0.10 对不同
#   尺度变量失效；(3) 与未标准化 health 相乘放大 R² 分解的列错配 bug。
# 修复：在  分支内（不动  文件）对 sample_A4 重新标准化为 SD=1。
#   仅在 sample_A4 上做——sample_B 含 6 个缺失者，scale() 会因 NA 污染整列。
#   person_std_mood_c 已是 SD≈1，保持不变。
for (.v in c("pre_phq9_c", "pre_control_c")) {
  if (.v %in% names(sample_A4)) {
    sample_A4[[.v]] <- as.numeric(scale(sample_A4[[.v]]))
  }
}
rm(.v)
cat(sprintf("  ✓ F2: sample_A4 的 pre_phq9_c / pre_control_c 已标准化 (SD: %.3f / %.3f → 1)\n",
            1, 1))

# ============================================================================
# 1Bb. 不变量断言 1: 样本规模
# ============================================================================
cat("\n  不变量断言:\n")
n_user_B <- nlevels(sample_B$user_id)
n_obs_B  <- nrow(sample_B)
n_user_A4 <- nlevels(sample_A4$user_id)
n_obs_A4 <- nrow(sample_A4)
n_user_L <- nlevels(sample_lagged$user_id)
n_obs_L  <- nrow(sample_lagged)

stopifnot("sample_B 人数应在 [200, 300]" = n_user_B >= 200L && n_user_B <= 300L)
stopifnot("sample_B 观测数应 > 8000" = n_obs_B > 8000L)
stopifnot("sample_A4 人数应在 [180, 300]" = n_user_A4 >= 180L && n_user_A4 <= 300L)
stopifnot("sample_lagged 观测数应 > 6000" = n_obs_L > 6000L)

cat(sprintf("    ✓ Sample B:      %d 人, %d 观测\n", n_user_B, n_obs_B))
cat(sprintf("    ✓ Sample A4:     %d 人, %d 观测\n", n_user_A4, n_obs_A4))
cat(sprintf("    ✓ Sample Lagged: %d 人, %d 观测\n", n_user_L, n_obs_L))

# ============================================================================
# 1Bc. 不变量断言 2: 关键  预测变量齐全
# ============================================================================
required_in_B <- c(
  "mood", "user_id",
  "health_session_c", "health_session_between_c",
  "unhealthy_session_c", "unhealthy_session_between_c",
  "balance_session_c",
  "mood_prev_session_c", "mood_lag1_day_c",
  "health_x_prev_mood_c",
  "study_day_c", "study_day_squared", "is_weekend", "session_numeric"
)
missing_in_B <- setdiff(required_in_B, names(sample_B))
if (length(missing_in_B) > 0) {
  stop(sprintf("sample_B 缺失变量: %s", paste(missing_in_B, collapse = ", ")))
}
cat(sprintf("    ✓ sample_B 包含所有 %d 个  必需变量\n",
            length(required_in_B)))

required_in_A4 <- c(required_in_B, "pre_phq9_c", "pre_control_c",
                    "person_std_mood_c")
missing_in_A4 <- setdiff(required_in_A4, names(sample_A4))
if (length(missing_in_A4) > 0) {
  stop(sprintf("sample_A4 缺失变量: %s", paste(missing_in_A4, collapse = ", ")))
}
cat(sprintf("    ✓ sample_A4 包含所有 %d 个  必需变量（含 σ^M）\n",
            length(required_in_A4)))

# ============================================================================
# 1Bd. 不变量断言 3: 关键变量描述统计 — 报告 Rule 13 要求详细输出
# ============================================================================
cat("\n  关键预测变量描述统计（sample_B / sample_A4）:\n")
safe_mean <- function(x) {
  if (is.numeric(x)) mean(x, na.rm = TRUE) else NA_real_
}
safe_sd <- function(x) {
  if (is.numeric(x)) sd(x, na.rm = TRUE) else NA_real_
}

desc_B <- data.frame(
  variable = required_in_B,
  M  = sapply(required_in_B, function(v) safe_mean(sample_B[[v]])),
  SD = sapply(required_in_B, function(v) safe_sd(sample_B[[v]])),
  N_valid = sapply(required_in_B, function(v) sum(!is.na(sample_B[[v]]))),
  N_NA = sapply(required_in_B, function(v) sum(is.na(sample_B[[v]])))
)
print(desc_B, row.names = FALSE, digits = 4)

# A4 独有变量
extra_A4 <- setdiff(required_in_A4, required_in_B)
desc_A4 <- data.frame(
  variable = extra_A4,
  M  = sapply(extra_A4, function(v) safe_mean(sample_A4[[v]])),
  SD = sapply(extra_A4, function(v) safe_sd(sample_A4[[v]])),
  N_valid = sapply(extra_A4, function(v) sum(!is.na(sample_A4[[v]]))),
  N_NA = sapply(extra_A4, function(v) sum(is.na(sample_A4[[v]])))
)
cat("\n  Sample A4 独有变量:\n")
print(desc_A4, row.names = FALSE, digits = 4)

# 描述统计写入 logs（供报告引用）
desc_path <- file.path(logs_dir, "desc_predictors_.csv")
write.csv(rbind(
  cbind(sample = "B",  desc_B),
  cbind(sample = "A4", desc_A4)
), desc_path, row.names = FALSE)
cat(sprintf("  ✓ 描述统计已保存: %s\n", desc_path))

# ============================================================================
# 1Be. 不变量断言 4: 因变量 mood 在合法范围 [-10, +10]
# ============================================================================
mood_range <- range(sample_B$mood, na.rm = TRUE)
stopifnot("mood 应在 [-10, +10] 闭区间" =
            mood_range[1] >= -10 && mood_range[2] <= 10)
cat(sprintf("    ✓ mood 范围: [%.1f, %.1f]\n", mood_range[1], mood_range[2]))

# ============================================================================
# 1Bf. 不变量断言 5 (): mood_lag1_day_c 在 sample_lagged 中无 NA
# ============================================================================
# 该变量是 03b MELS 模型的核心 AR 控制项。若含 NA，brms 静默删除对应行，
# 拟合样本量与上方报告的 n_obs_L 会不一致，且用户无法察觉。
n_na_lag <- sum(is.na(sample_lagged$mood_lag1_day_c))
stopifnot("sample_lagged 的 mood_lag1_day_c 必须无 NA" = n_na_lag == 0)
cat(sprintf("    ✓ mood_lag1_day_c 在 sample_lagged 中无 NA (n_NA=0)\n"))

cat("\n  ✓ 模块01b 数据导入与不变量断言完成\n")
