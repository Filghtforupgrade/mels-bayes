# ============================================================================
# 模块 01: 数据准备 — 数据加载 + 变量构建(CWC/CGM/Lag) + 样本定义
# ============================================================================
# 依赖: 00_config.R (project_root, data_dir, output_dir, result containers)
# 产出: session_df, person_vars, sample_B, sample_A4, sample_lagged
#       noon_behaviors, evening_behaviors, all_behaviors
# ============================================================================

cat("\n--- 模块01: 数据加载与变量构建 ---\n")

# ============================================================================
# 1A. 数据加载
# ============================================================================
cat("\n【数据加载】\n")

session_raw <- read.csv(file.path(data_dir, "ema_model_input_session.csv"), stringsAsFactors = FALSE)
person_raw <- read.csv(file.path(data_dir, "ema_model_input_person.csv"), stringsAsFactors = FALSE)
day_raw <- read.csv(file.path(data_dir, "ema_model_input_day.csv"), stringsAsFactors = FALSE)

session_df <- session_raw %>%
  filter(!is.na(mood)) %>%
  mutate(user_id = as.character(user_id))

person_vars <- person_raw %>%
  select(user_id, pre_phq9_total, pre_perceived_control_total,
         pre_emotion_controllability_total,
         included_sample_C, has_posttest, delta_phq9,
         overall_completion_rate,
         included_sample_B,
         pre_seasonality_total,
         commit_emotion, commit_routine, commit_selfcare,
         post_phq9_total) %>%
  mutate(user_id = as.character(user_id))

session_df <- session_df %>%
  left_join(person_vars, by = "user_id") %>%
  filter(!is.na(included_sample_B)) %>%
  mutate(user_id = as.factor(user_id))

cat(sprintf("  Session数据: %d人, %d观测\n", nlevels(session_df$user_id), nrow(session_df)))
cat(sprintf("  Person数据: %d人\n", nrow(person_vars)))

# ============================================================================
# 1B. 构建共享变量
# ============================================================================
cat("\n【构建共享变量】\n")

# --- 行为分类（: Session-Aware定义）---
# 午间行为 (session_noon=1): sleep_late, wake_early, breakfast, aspiration
# 晚间行为 (session_noon=0): walk, exercise, dinner, breathing
noon_behaviors <- c('sleep_late', 'wake_early', 'breakfast', 'aspiration')
evening_behaviors <- c('walk', 'exercise', 'dinner', 'breathing')
all_behaviors <- c(noon_behaviors, evening_behaviors)

# ★ : 不做NA→0填充。行为变量在设计性缺失的session保持NA，
# 仅对实际测量的session项目求和，消除跨session测量偏差。

# --- Session-Aware组合变量 ---
# 午间健康行为: 夜卧(否)+早起(是)+早餐(是)+立志(是) = 0-4
# 晚间健康行为: 散步(是)+锻炼(是)+晚餐(是)+呼吸练习(是) = 0-4
session_df <- session_df %>%
  mutate(
    health_noon  = (1 - sleep_late) + wake_early + breakfast + aspiration,
    health_evening = walk + exercise + dinner + breathing,
    unhealthy_noon = sleep_late + (1 - wake_early) + (1 - breakfast) + (1 - aspiration),
    unhealthy_evening = (1 - dinner)
  )

# 统一为session-aware变量
session_df <- session_df %>%
  mutate(
    health_session = ifelse(session_noon == 1, health_noon, health_evening),
    unhealthy_session = ifelse(session_noon == 1, unhealthy_noon, unhealthy_evening)
  )

# 行为平衡指数（午间/晚间分别计算比例）
session_df <- session_df %>%
  mutate(
    balance_noon = health_noon / 4 - unhealthy_noon / 4,
    balance_evening = health_evening / 4 - unhealthy_evening / 1,
    balance_session = ifelse(session_noon == 1, balance_noon, balance_evening)
  )

# 高健康/有不健康行为（基于session-aware组合，午间和晚间最大均为4）
session_df <- session_df %>%
  mutate(
    high_health_behavior = ifelse(health_session >= 3, 1, 0),
    has_unhealthy_behavior = ifelse(unhealthy_session >= 1, 1, 0)
  )

# --- Session-Specific CWC（组内中心化）---
# 按用户×session分别中心化，确保"偏离平常"相对于相同时段的自己
session_df <- session_df %>%
  group_by(user_id) %>%
  mutate(
    health_noon_c  = health_noon - mean(health_noon[session_noon == 1], na.rm = TRUE),
    health_evening_c = health_evening - mean(health_evening[session_noon == 0], na.rm = TRUE),
    unhealthy_noon_c = unhealthy_noon - mean(unhealthy_noon[session_noon == 1], na.rm = TRUE),
    unhealthy_evening_c = unhealthy_evening - mean(unhealthy_evening[session_noon == 0], na.rm = TRUE),
    balance_noon_c = balance_noon - mean(balance_noon[session_noon == 1], na.rm = TRUE),
    balance_evening_c = balance_evening - mean(balance_evening[session_noon == 0], na.rm = TRUE)
  ) %>%
  ungroup()

# 统一为session-aware CWC变量
session_df <- session_df %>%
  mutate(
    health_session_c = ifelse(session_noon == 1, health_noon_c, health_evening_c),
    unhealthy_session_c = ifelse(session_noon == 1, unhealthy_noon_c, unhealthy_evening_c),
    balance_session_c = ifelse(session_noon == 1, balance_noon_c, balance_evening_c)
  )

# 高健康/不健康行为的CWC版本
session_df <- session_df %>%
  group_by(user_id) %>%
  mutate(
    high_health_behavior_c = high_health_behavior - mean(high_health_behavior, na.rm = TRUE),
    has_unhealthy_behavior_c = has_unhealthy_behavior - mean(has_unhealthy_behavior, na.rm = TRUE)
  ) %>%
  ungroup()

# 单个行为CWC化（仅基于有效session的person-mean，NA保持NA）
for (beh in all_behaviors) {
  if (beh %in% names(session_df)) {
    beh_c_name <- paste0(beh, "_c")
    session_df[[beh_c_name]] <- ave(session_df[[beh]], session_df$user_id,
                                     FUN = function(x) x - mean(x, na.rm = TRUE))
  }
}

# positive_session: 积极行为 = 健康行为总和（用于B12/R5，原positive_sum已融入health_session）
session_df <- session_df %>%
  mutate(positive_session = health_session)

session_df <- session_df %>%
  group_by(user_id) %>%
  mutate(positive_session_c = positive_session - mean(positive_session, na.rm = TRUE)) %>%
  ungroup()

# --- 情绪滞后变量CWC化 ---
cat("\n  ★情绪滞后变量CWC化\n")

session_df <- session_df %>%
  group_by(user_id) %>%
  mutate(person_mean_mood = mean(mood, na.rm = TRUE)) %>%
  ungroup()

if ("mood_lag1_day" %in% names(session_df)) {
  session_df <- session_df %>%
    mutate(mood_lag1_day_c = mood_lag1_day - person_mean_mood)
  cat("    mood_lag1_day_c 已创建（个人均值中心化）\n")
} else {
  session_df <- session_df %>%
    arrange(user_id, study_day_c, session_numeric) %>%
    group_by(user_id) %>%
    mutate(
      mood_lag1_day = lag(mood, 1),
      mood_lag1_day_c = mood_lag1_day - person_mean_mood
    ) %>%
    ungroup()
  cat("    mood_lag1_day 和 mood_lag1_day_c 已创建\n")
}

# --- : mood_prev_session_c 情绪基线控制 ---
# 对标Lischetzke et al. (2021) morning calmness控制策略
cat("\n  ★: mood_prev_session_c 情绪基线控制\n")

if ("mood_prev_session" %in% names(session_df)) {
  session_df <- session_df %>%
    mutate(mood_prev_session_c = mood_prev_session - person_mean_mood)
} else {
  session_df <- session_df %>%
    arrange(user_id, study_day_c, session_numeric) %>%
    group_by(user_id) %>%
    mutate(
      mood_prev_session = lag(mood, 1),
      mood_prev_session_c = mood_prev_session - person_mean_mood
    ) %>%
    ungroup()
}

# 首次测量标记（无前次情绪）
session_df <- session_df %>%
  group_by(user_id) %>%
  mutate(
    is_first_session = ifelse(is.na(mood_prev_session), 1, 0),
    mood_prev_session_c = ifelse(is.na(mood_prev_session_c), 0, mood_prev_session_c)
  ) %>%
  ungroup()

n_first <- sum(session_df$is_first_session)
cat(sprintf("    mood_prev_session_c 已创建（CWC，首次测量=%d设为0+标记）\n", n_first))

# --- 个体间均值（L2 Person-Mean）---
session_df <- session_df %>%
  group_by(user_id) %>%
  mutate(
    health_session_between = mean(health_session, na.rm = TRUE),
    unhealthy_session_between = mean(unhealthy_session, na.rm = TRUE)
  ) %>%
  ungroup() %>%
  mutate(
    health_session_between_c = health_session_between - mean(health_session_between, na.rm = TRUE),
    unhealthy_session_between_c = unhealthy_session_between - mean(unhealthy_session_between, na.rm = TRUE)
  )

# --- 时间变量 ---
if ("day_of_study_centered" %in% names(session_df)) {
  session_df$study_day_c <- ifelse(is.na(session_df$day_of_study_centered), 0, session_df$day_of_study_centered)
} else if ("study_day_c" %in% names(session_df)) {
  session_df$study_day_c <- ifelse(is.na(session_df$study_day_c), 0, session_df$study_day_c)
} else {
  session_df <- session_df %>%
    group_by(user_id) %>%
    mutate(study_day_c = row_number() - 10) %>%
    ungroup()
}

session_df$study_day_squared <- session_df$study_day_c ^ 2

if (!"is_weekend" %in% names(session_df)) {
  if ("weekday" %in% names(session_df)) {
    session_df$weekday <- as.numeric(session_df$weekday)
    session_df$is_weekend <- ifelse(session_df$weekday >= 5, 1, 0)
  } else {
    session_df$is_weekend <- 0
  }
}

if (!"session_numeric" %in% names(session_df)) {
  if ("session" %in% names(session_df)) {
    session_df$session_numeric <- ifelse(session_df$session == "evening", 1, 0)
  } else if ("session_noon" %in% names(session_df)) {
    session_df$session_numeric <- session_df$session_noon
  } else {
    session_df$session_numeric <- 0
  }
}

# --- 个体间总均值中心化 (CGM) ---
session_df <- session_df %>%
  mutate(
    pre_phq9_c = pre_phq9_total - mean(pre_phq9_total, na.rm = TRUE),
    pre_control_c = pre_perceived_control_total - mean(pre_perceived_control_total, na.rm = TRUE),
    pre_emotion_c = pre_emotion_controllability_total - mean(pre_emotion_controllability_total, na.rm = TRUE),
    pre_seasonality_c = pre_seasonality_total - mean(pre_seasonality_total, na.rm = TRUE)
  )

# --- 情绪调节交互项 (阻断效应) ---
session_df <- session_df %>%
  mutate(
    health_x_prev_mood_c = health_session_c * mood_lag1_day_c,
    unhealthy_x_prev_mood_c = unhealthy_session_c * mood_lag1_day_c,
    positive_x_prev_mood_c = positive_session_c * mood_lag1_day_c
  )

for (beh in all_behaviors) {
  beh_c_var <- paste0(beh, "_c")
  if (beh_c_var %in% names(session_df)) {
    session_df[[paste0(beh, "_x_prev_mood_c")]] <- session_df[[beh_c_var]] * session_df$mood_lag1_day_c
  }
}

# --- 兼容保留变量（模型不使用）---
session_df <- session_df %>%
  mutate(
    health_sum = health_session,
    unhealthy_sum = unhealthy_session,
    health_behavior_sum = health_session,
    unhealthy_behavior_sum = unhealthy_session,
    health_behavior_ratio = health_session / 4,
    unhealthy_behavior_ratio = unhealthy_session / 4,
    balance_c = balance_session_c,
    health_sum_x_day = health_session_c * study_day_c,
    unhealthy_sum_x_day = unhealthy_session_c * study_day_c,
    balance_x_day = balance_session_c * study_day_c,
    health_sum_x_phq9 = health_session_c * pre_phq9_c,
    unhealthy_sum_x_phq9 = unhealthy_session_c * pre_phq9_c,
    balance_x_phq9 = balance_session_c * pre_phq9_c,
    positive_behavior_sum = health_session,
    positive_sum = health_session,
    health_sum_c = health_session_c,
    unhealthy_sum_c = unhealthy_session_c,
    positive_sum_c = positive_session_c
  )

session_df <- session_df %>%
  group_by(user_id) %>%
  mutate(positive_behavior_c = positive_behavior_sum - mean(positive_behavior_sum, na.rm = TRUE)) %>%
  ungroup()

# --- 交叉时滞变量 (Behavior_{t-1} -> Emotion_t) ---
session_df <- session_df %>%
  arrange(user_id, study_day_c, session_numeric) %>%
  group_by(user_id) %>%
  mutate(
    health_lag1 = lag(health_session_c, 1),
    unhealthy_lag1 = lag(unhealthy_session_c, 1),
    positive_lag1 = lag(positive_session_c, 1)
  ) %>%
  ungroup()

session_df <- session_df %>% arrange(user_id, study_day_c, session_numeric)
for (beh in all_behaviors) {
  beh_c_var <- paste0(beh, "_c")
  beh_lag1_name <- paste0(beh, "_lag1")
  if (beh_c_var %in% names(session_df)) {
    session_df[[beh_lag1_name]] <- ave(session_df[[beh_c_var]], session_df$user_id, 
                                        FUN = function(x) c(NA, head(x, -1)))
  }
}

cat("  ✓ 共享变量构建完成（: Session-Aware组合 + Session-Specific CWC + mood_prev_session_c）\n")

# ============================================================================
# 1C. 样本定义
# ============================================================================
sample_B <- session_df %>% filter(included_sample_B == 1) %>% droplevels()
sample_A4 <- sample_B %>% filter(!is.na(pre_phq9_c), !is.na(pre_control_c)) %>% droplevels()
sample_lagged <- sample_B %>% filter(!is.na(health_lag1), !is.na(mood_lag1_day_c)) %>% droplevels()

cat(sprintf("\nSample B: %d 人, %d 观测\n", nlevels(sample_B$user_id), nrow(sample_B)))
cat(sprintf("Sample A4: %d 人, %d 观测\n", nlevels(sample_A4$user_id), nrow(sample_A4)))
cat(sprintf("Sample Lagged: %d 人, %d 观测\n", nlevels(sample_lagged$user_id), nrow(sample_lagged)))

# ============================================================================
# 1D. : 预计算 person_std_mood_c (情绪波动性Level-2变量)
# ============================================================================
cat("\n  ★: 预计算情绪波动性 (person_std_mood_c)\n")

sample_A4 <- sample_A4 %>%
  group_by(user_id) %>%
  mutate(person_std_mood_val = sd(mood, na.rm = TRUE)) %>%
  ungroup() %>%
  mutate(person_std_mood_c = person_std_mood_val - mean(person_std_mood_val, na.rm = TRUE))

cat(sprintf("    person_std_mood: M=%.4f SD=%.4f NA=%d\n",
            mean(sample_A4$person_std_mood_val, na.rm = TRUE),
            sd(sample_A4$person_std_mood_val, na.rm = TRUE),
            sum(is.na(sample_A4$person_std_mood_c))))

# ============================================================================
# 1E. 保存检查点
# ============================================================================
save(session_df, person_vars, sample_B, sample_A4, sample_lagged,
     noon_behaviors, evening_behaviors, all_behaviors,
     session_raw,
     file = file.path(checkpoint_dir, "checkpoint_01_data.RData"))
cat("  ✓ 检查点已保存: checkpoint_01_data.RData\n")

cat("\n  ✓ 模块01 数据准备完成\n")
