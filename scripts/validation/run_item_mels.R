#!/usr/bin/env Rscript
# ============================================================================
#  条目级 MELS 细化分析 — 午间 + 晚间，主/稳健性双控制集
# ============================================================================
# Revision tag:  (aligned with manuscript + sensitivity suite, 2026-05-28)
#   - : manuscript rewrite; figure regeneration (A9/A11)
#   - : prior sensitivity, MI sensitivity, conditional effects
#   -  (prior synced with main pipeline; see ledger R3-10)
#   - priors_item 与 02_modules/02_priors.R::priors_mels 完全一致：
#       Intercept       ← normal(3, 2)         (R3-01)
#       Intercept σ-dpar← student_t(3, 0, 1)   (Plan B 维持)
#       sd              ← student_t(3, 0, 1)   (R3-02)
#       b / b σ-dpar    ← normal(0, 1)
#   - 替代  的 B2_MELS 与 B3_MELS（已确认 L1 共线性接近恒等式，R2-01/02）
#   - 作为 B1_MELS 的细化分析：识别 8 个行为中哪几个驱动了组合效应
#
# 模型矩阵 (4 个 brmsfit):
#   Item_MELS_Noon_Main      4 午间条目 + [mood_lag1_day_c, study_day_c, health_session_between_c]
#   Item_MELS_Noon_Robust    Noon_Main + mood_prev_session_c
#   Item_MELS_Evening_Main   4 晚间条目 + 同主控制集
#   Item_MELS_Evening_Robust Evening_Main + mood_prev_session_c
#
# 编码方向（统一为"做了=正值，预期与 mood 正相关"）:
#   sleep_ontime_c = -sleep_late_c    （按时睡=+1）
#   其余 7 个条目保持原 CWC（做了=+1）
#
# 输出目录命名: item_mels_3_<TIMESTAMP>/（ 起的命名约定；
# 旧的  /  子目录保留为 audit trail）
#
# 用法:
#   cd /path/to/EMA_MELS
#   Rscript scripts/validation/run_item_mels.R
#
#   # 干跑（chains=2 iter=500 verify formula only，约 2 min/model）
#   ITEM_MELS_MODE=dry Rscript ...
#
#   # 选择性运行（其他用环境变量已存的 rds）
#   ITEM_MELS_MODELS=Item_MELS_Noon_Main Rscript ...
# ============================================================================

options(width = 200)
options(warn = 1)

# ============================================================================
# 0. 路径与环境
# ============================================================================
current_file <- commandArgs(trailingOnly = FALSE)
script_path <- sub("^--file=", "", current_file[grep("^--file=", current_file)])
if (exists("project_root", inherits = TRUE) && dir.exists(get("project_root", inherits = TRUE))) {
  # 被 source() 调用时：父脚本已定义 project_root，直接复用
  SCRIPT_DIR <- file.path(project_root, "scripts", "validation")
} else if (length(script_path) > 0 && file.exists(script_path)) {
  SCRIPT_DIR <- dirname(normalizePath(script_path))
} else {
  SCRIPT_DIR <- getwd()
}

#  项目根目录（脚本在 scripts/validation/ 下）
if (!exists("project_root", inherits = TRUE) || !dir.exists(get("project_root", inherits = TRUE))) {
  project_root <- normalizePath(file.path(SCRIPT_DIR, "..", ".."))
}

suppressPackageStartupMessages({
  library(brms)
  library(cmdstanr)
  library(posterior)
  library(bayesplot)
  library(ggplot2)
  library(dplyr)
})

MODE <- Sys.getenv("ITEM_MELS_MODE", "production")
DRY  <- identical(MODE, "dry")

MODELS_ENV <- Sys.getenv("ITEM_MELS_MODELS", "")
SELECT_MODELS <- if (nchar(MODELS_ENV) > 0) {
  trimws(strsplit(MODELS_ENV, ",")[[1]])
} else {
  c("Item_MELS_Noon_Main",
    "Item_MELS_Noon_Robust",
    "Item_MELS_Evening_Main",
    "Item_MELS_Evening_Robust")
}

# 时间戳输出目录（项目规则：自含归档）
TS <- format(Sys.time(), "%Y%m%d_%H%M%S")
output_dir <- file.path(project_root, "results", sprintf("item_mels_3_%s", TS))
fits_dir   <- file.path(output_dir, "fits")
logs_dir   <- file.path(output_dir, "logs")
posterior_dir <- file.path(output_dir, "posterior_summaries")
figures_dir   <- file.path(output_dir, "figures")
for (d in c(output_dir, fits_dir, logs_dir, posterior_dir, figures_dir)) {
  dir.create(d, showWarnings = FALSE, recursive = TRUE)
}

cat("\n==============================================================\n")
cat("   条目级 MELS 细化分析\n")
cat("==============================================================\n")
cat(sprintf("  开始时间    : %s\n", format(Sys.time(), "%Y-%m-%d %H:%M:%S")))
cat(sprintf("  模式        : %s%s\n", MODE, if (DRY) " (chains=2, iter=500)" else " (chains=4, iter=6000)"))
cat(sprintf("  待跑模型    : %s\n", paste(SELECT_MODELS, collapse = ", ")))
cat(sprintf("  输出目录    : %s\n", output_dir))

# ============================================================================
# 1. 数据加载与编码方向统一
# ============================================================================
ckpt_path <- file.path(project_root, "data", "checkpoints", "checkpoint_01_data.RData")
if (!file.exists(ckpt_path)) stop("checkpoint 不存在: ", ckpt_path)
load(ckpt_path)

if (!"sample_lagged" %in% ls()) stop("sample_lagged 未在 checkpoint 中找到")

# sleep_ontime_c = -sleep_late_c（按时睡=+1，与其余 7 条目方向一致）
sample_lagged <- sample_lagged %>%
  mutate(sleep_ontime_c = -sleep_late_c)

sample_noon <- sample_lagged %>% filter(session_noon == 1)
sample_eve  <- sample_lagged %>% filter(session_noon == 0)

# drop unused user_id levels（factor 残留会让 brms 警告）
sample_noon <- sample_noon %>% mutate(user_id = droplevels(user_id))
sample_eve  <- sample_eve  %>% mutate(user_id = droplevels(user_id))

cat(sprintf("\n  数据 (sample_lagged):\n"))
cat(sprintf("    午间样本: %d 人, %d 观测\n", nlevels(sample_noon$user_id), nrow(sample_noon)))
cat(sprintf("    晚间样本: %d 人, %d 观测\n", nlevels(sample_eve$user_id),  nrow(sample_eve)))

# 验证条目方向（sleep_ontime_c 与 mood 应弱正相关）
noon_items <- c("sleep_ontime_c", "wake_early_c", "breakfast_c", "aspiration_c")
eve_items  <- c("walk_c", "exercise_c", "dinner_c", "breathing_c")

cat("\n  L1 条目与 mood 的零阶相关（CWC 已移除个体均值差异）:\n")
for (it in noon_items) {
  r <- cor(sample_noon[[it]], sample_noon$mood, use = "complete.obs")
  cat(sprintf("    [午间] %-16s r(mood) = %+.4f  SD=%.3f\n",
              it, r, sd(sample_noon[[it]], na.rm = TRUE)))
}
for (it in eve_items) {
  r <- cor(sample_eve[[it]], sample_eve$mood, use = "complete.obs")
  cat(sprintf("    [晚间] %-16s r(mood) = %+.4f  SD=%.3f\n",
              it, r, sd(sample_eve[[it]], na.rm = TRUE)))
}

# ============================================================================
# 2. 先验（ 与 02b priors_mels 完全一致）
# ============================================================================
# A: Intercept normal(3, 2) 校准到 mood 量表中点 M≈3.31
# B: sd student_t(3, 0, 1) 折中收紧
# C: σ_Intercept student_t(3, 0, 1) Plan B 收敛保护维持
priors_item <- c(
  brms::prior(normal(0, 1),         class = "b"),
  brms::prior(normal(0, 1),         class = "b",         dpar = "sigma"),
  brms::prior(normal(3, 2),         class = "Intercept"),                # A
  brms::prior(student_t(3, 0, 1),   class = "Intercept", dpar = "sigma"),# C
  brms::prior(student_t(3, 0, 1),   class = "sd")                        # B
)

# ============================================================================
# 3. 公式定义
# ============================================================================
# 主控制集
form_noon_main <- brms::bf(
  mood  ~ sleep_ontime_c + wake_early_c + breakfast_c + aspiration_c +
          mood_lag1_day_c + study_day_c + health_session_between_c +
          (1 | user_id),
  sigma ~ sleep_ontime_c + wake_early_c + breakfast_c + aspiration_c
)

form_noon_robust <- brms::bf(
  mood  ~ sleep_ontime_c + wake_early_c + breakfast_c + aspiration_c +
          mood_lag1_day_c + mood_prev_session_c + study_day_c +
          health_session_between_c +
          (1 | user_id),
  sigma ~ sleep_ontime_c + wake_early_c + breakfast_c + aspiration_c
)

form_eve_main <- brms::bf(
  mood  ~ walk_c + exercise_c + dinner_c + breathing_c +
          mood_lag1_day_c + study_day_c + health_session_between_c +
          (1 | user_id),
  sigma ~ walk_c + exercise_c + dinner_c + breathing_c
)

form_eve_robust <- brms::bf(
  mood  ~ walk_c + exercise_c + dinner_c + breathing_c +
          mood_lag1_day_c + mood_prev_session_c + study_day_c +
          health_session_between_c +
          (1 | user_id),
  sigma ~ walk_c + exercise_c + dinner_c + breathing_c
)

# ============================================================================
# 4. 拟合配置（dry vs production）
# ============================================================================
if (DRY) {
  n_chains <- 2L; n_iter <- 500L;  n_warmup <- 250L
} else {
  n_chains <- 4L; n_iter <- 6000L; n_warmup <- 3000L
}
mcmc_ctrl <- list(adapt_delta = 0.95, max_treedepth = 12L)
GLOBAL_SEED <- 20260519L

fit_one <- function(name, formula, data) {
  rds_path <- file.path(fits_dir, paste0(name, ".rds"))
  log_path <- file.path(logs_dir, paste0(name, "_cmdstan.log"))

  if (file.exists(rds_path)) {
    cat(sprintf("\n  [%s] 加载已有 .rds\n", name))
    return(readRDS(rds_path))
  }

  cat(sprintf("\n  [%s] 拟合中（chains=%d, iter=%d, warmup=%d, threads=2）...\n",
              name, n_chains, n_iter, n_warmup))
  t0 <- Sys.time()

  con <- file(log_path, open = "wt")
  sink(con, type = "message")
  on.exit({ sink(NULL, type = "message"); close(con) }, add = TRUE)

  fit <- tryCatch(
    brms::brm(
      formula = formula,
      data    = data,
      family  = brms::brmsfamily("gaussian"),
      prior   = priors_item,
      chains  = n_chains,
      iter    = n_iter,
      warmup  = n_warmup,
      cores   = n_chains,
      seed    = GLOBAL_SEED,
      control = mcmc_ctrl,
      backend = "cmdstanr",
      save_pars = brms::save_pars(all = TRUE),
      refresh = if (DRY) 100 else 200,
      threads = brms::threading(2L),
      stan_model_args = list(cpp_options = list(stan_threads = TRUE)),
      sample_prior = if (identical(name, "Item_MELS_Noon_Main")) "yes" else "no"
    ),
    error = function(e) {
      cat(sprintf("  [%s] ✗ 拟合错误: %s\n", name, conditionMessage(e)))
      NULL
    }
  )

  if (!is.null(fit)) {
    saveRDS(fit, rds_path, compress = "xz")
    elapsed <- difftime(Sys.time(), t0, units = "mins")
    cat(sprintf("  [%s] ✓ 完成，耗时 %.2f 分钟\n", name, as.numeric(elapsed)))
  }
  fit
}

# ============================================================================
# 5. 拟合循环
# ============================================================================
model_specs <- list(
  Item_MELS_Noon_Main      = list(form = form_noon_main,   data = sample_noon),
  Item_MELS_Noon_Robust    = list(form = form_noon_robust, data = sample_noon),
  Item_MELS_Evening_Main   = list(form = form_eve_main,    data = sample_eve),
  Item_MELS_Evening_Robust = list(form = form_eve_robust,  data = sample_eve)
)

fits <- list()
for (mname in SELECT_MODELS) {
  if (!mname %in% names(model_specs)) {
    cat(sprintf("\n  ⚠ 未知模型: %s（跳过）\n", mname))
    next
  }
  spec <- model_specs[[mname]]
  fits[[mname]] <- fit_one(mname, spec$form, spec$data)
}

if (length(Filter(Negate(is.null), fits)) == 0) {
  cat("\n  ✗ 无任何模型成功拟合，终止\n")
  # 被 source() 调用时不 quit()（避免杀死父脚本），改用 stop()
  if (interactive() || sys.nframe() > 0L) {
    stop("Item-level MELS: 无任何模型成功拟合")
  } else {
    quit(status = 1)
  }
}

# ============================================================================
# 6. 收敛诊断
# ============================================================================
cat("\n--- 收敛诊断 ---\n")

summarise_fit <- function(fit) {
  posterior::summarise_draws(
    posterior::as_draws_df(fit),
    mean, median, sd,
    ~ quantile(.x, probs = c(0.025, 0.975), na.rm = TRUE),
    posterior::default_convergence_measures()
  )
}

diag_summary <- function(fit, s, name) {
  np <- bayesplot::nuts_params(fit)
  div <- sum(np$Value[np$Parameter == "divergent__"], na.rm = TRUE)
  max_rhat <- max(s$rhat, na.rm = TRUE)
  min_ess  <- min(s$ess_bulk, na.rm = TRUE)
  worst_rhat_var <- s$variable[which.max(s$rhat)]
  worst_ess_var  <- s$variable[which.min(s$ess_bulk)]

  cat(sprintf("\n  %s:\n", name))
  cat(sprintf("    max R-hat  = %.4f (%s)\n", max_rhat, worst_rhat_var))
  cat(sprintf("    min ESS    = %.0f (%s)\n", min_ess, worst_ess_var))
  cat(sprintf("    divergent  = %d / %d (%.2f%%)\n",
              div, brms::ndraws(fit), 100 * div / brms::ndraws(fit)))
  invisible(data.frame(
    model = name,
    max_rhat = max_rhat,
    min_ess  = min_ess,
    worst_rhat_var = worst_rhat_var,
    worst_ess_var  = worst_ess_var,
    n_divergent = div,
    n_draws = brms::ndraws(fit),
    pct_divergent = 100 * div / brms::ndraws(fit)
  ))
}

posterior_list <- list()
conv_rows <- list()
for (mname in names(fits)) {
  fit <- fits[[mname]]
  if (is.null(fit)) next
  s <- summarise_fit(fit)
  posterior_list[[mname]] <- s
  write.csv(s, file.path(posterior_dir, paste0(mname, "_posterior.csv")),
            row.names = FALSE)
  conv_rows[[mname]] <- diag_summary(fit, s, mname)
}

conv_df <- do.call(rbind, conv_rows)
write.csv(conv_df, file.path(output_dir, "convergence_summary.csv"), row.names = FALSE)

# ============================================================================
# 7. 提取条目级效应（mu + sigma）
# ============================================================================
extract_item_effects <- function(s, items, model_label) {
  out <- data.frame(
    model  = model_label,
    item   = items,
    stringsAsFactors = FALSE
  )
  pick <- function(prefix, item, col) {
    r <- s[s$variable == paste0(prefix, item), ]
    if (nrow(r) == 0) NA_real_ else as.numeric(r[[col]])
  }
  out$mu_median     <- sapply(items, function(x) pick("b_", x, "median"))
  out$mu_q2.5       <- sapply(items, function(x) pick("b_", x, "2.5%"))
  out$mu_q97.5      <- sapply(items, function(x) pick("b_", x, "97.5%"))
  out$mu_rhat       <- sapply(items, function(x) pick("b_", x, "rhat"))
  out$mu_ess_bulk   <- sapply(items, function(x) pick("b_", x, "ess_bulk"))

  out$sig_log_median <- sapply(items, function(x) pick("b_sigma_", x, "median"))
  out$sig_log_q2.5   <- sapply(items, function(x) pick("b_sigma_", x, "2.5%"))
  out$sig_log_q97.5  <- sapply(items, function(x) pick("b_sigma_", x, "97.5%"))
  out$sig_rhat       <- sapply(items, function(x) pick("b_sigma_", x, "rhat"))
  out$sig_ess_bulk   <- sapply(items, function(x) pick("b_sigma_", x, "ess_bulk"))

  out$sig_mult_median <- exp(out$sig_log_median)
  out$sig_mult_q2.5   <- exp(out$sig_log_q2.5)
  out$sig_mult_q97.5  <- exp(out$sig_log_q97.5)
  out$sig_compression_pct <- (1 - out$sig_mult_median) * 100

  # CrI 是否排除 0 → 显著性指示
  out$mu_excludes_zero  <- (out$mu_q2.5 > 0) | (out$mu_q97.5 < 0)
  out$sig_excludes_zero <- (out$sig_log_q2.5 > 0) | (out$sig_log_q97.5 < 0)

  rownames(out) <- NULL
  out
}

items_by_model <- list(
  Item_MELS_Noon_Main      = noon_items,
  Item_MELS_Noon_Robust    = noon_items,
  Item_MELS_Evening_Main   = eve_items,
  Item_MELS_Evening_Robust = eve_items
)

effect_rows <- list()
for (mname in names(posterior_list)) {
  effect_rows[[mname]] <- extract_item_effects(
    posterior_list[[mname]],
    items_by_model[[mname]],
    mname
  )
}
all_effects <- do.call(rbind, effect_rows)
rownames(all_effects) <- NULL

write.csv(all_effects, file.path(output_dir, "item_mels_effects.csv"),
          row.names = FALSE)

cat("\n========================================\n")
cat("  条目级效应总表（mu + sigma）\n")
cat("========================================\n")
print(all_effects[, c("model", "item", "mu_median", "mu_q2.5", "mu_q97.5",
                       "mu_excludes_zero",
                       "sig_log_median", "sig_mult_median",
                       "sig_compression_pct", "sig_excludes_zero")],
      row.names = FALSE, digits = 4)

# ============================================================================
# 8. 可视化（forest plot）
# ============================================================================
all_effects$label <- paste0(all_effects$model, ": ", all_effects$item)
all_effects$session <- ifelse(grepl("Noon", all_effects$model), "Noon", "Evening")
all_effects$spec    <- ifelse(grepl("Robust", all_effects$model), "Robust", "Main")

# mu forest
p_mu <- ggplot(all_effects, aes(x = mu_median, y = reorder(label, mu_median),
                                 color = spec, shape = session)) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "grey40") +
  geom_errorbarh(aes(xmin = mu_q2.5, xmax = mu_q97.5), height = 0.2) +
  geom_point(size = 2.5) +
  scale_color_manual(values = c("Main" = "steelblue", "Robust" = "darkorange")) +
  labs(title = "Item-level MELS: mu (location of mood)",
       subtitle = "Posterior median + 95% CrI; positive = behavior raises mood",
       x = "Posterior median", y = NULL) +
  theme_minimal(base_size = 11) +
  theme(legend.position = "bottom")
ggsave(file.path(figures_dir, "fig_item_mels_mu.png"),
       p_mu, width = 10, height = 6, dpi = 120)

# sigma forest
all_effects$sig_dir <- ifelse(all_effects$sig_log_median < 0, "compression", "inflation")
p_sig <- ggplot(all_effects,
                aes(x = sig_log_median, y = reorder(label, sig_log_median),
                    color = sig_dir, shape = session)) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "grey40") +
  geom_errorbarh(aes(xmin = sig_log_q2.5, xmax = sig_log_q97.5), height = 0.2) +
  geom_point(size = 2.5) +
  scale_color_manual(values = c("compression" = "steelblue", "inflation" = "tomato")) +
  labs(title = "Item-level MELS: sigma (scale of mood, log-link)",
       subtitle = "Posterior median + 95% CrI; negative = compresses variability",
       x = "Posterior median (log-SD)", y = NULL) +
  theme_minimal(base_size = 11) +
  theme(legend.position = "bottom")
ggsave(file.path(figures_dir, "fig_item_mels_sigma.png"),
       p_sig, width = 10, height = 6, dpi = 120)

cat(sprintf("\n  ✓ 效应表与图保存至 %s\n", output_dir))

# ============================================================================
# 9. README 自动生成（claim → source 映射）
# ============================================================================
readme_path <- file.path(output_dir, "README.md")

n_noon_obs <- nrow(sample_noon); n_noon_user <- nlevels(sample_noon$user_id)
n_eve_obs  <- nrow(sample_eve);  n_eve_user  <- nlevels(sample_eve$user_id)

cat("#  Item-level MELS Results\n\n", file = readme_path)
cat(sprintf("**Run timestamp**: %s\n", TS), file = readme_path, append = TRUE)
cat(sprintf("**Mode**: %s\n", MODE), file = readme_path, append = TRUE)
cat("**Revision tag**:  (aligned with manuscript + sensitivity suite)\n",
    file = readme_path, append = TRUE)
cat(sprintf("**Models fitted**: %s\n\n", paste(names(fits), collapse = ", ")),
    file = readme_path, append = TRUE)

cat("## 数据\n\n", file = readme_path, append = TRUE)
cat(sprintf("- 数据源：`data/checkpoints/checkpoint_01_data.RData` → `sample_lagged`\n"),
    file = readme_path, append = TRUE)
cat(sprintf("- 午间子样本：%d 人，%d 观测\n", n_noon_user, n_noon_obs),
    file = readme_path, append = TRUE)
cat(sprintf("- 晚间子样本：%d 人，%d 观测\n\n", n_eve_user, n_eve_obs),
    file = readme_path, append = TRUE)

cat("## 模型规格\n\n", file = readme_path, append = TRUE)
cat("**编码方向统一**：`sleep_ontime_c = -sleep_late_c`，所有条目正系数=促进 mood\n\n",
    file = readme_path, append = TRUE)
cat("| 模型 | 时段 | 条目 | 额外控制 |\n", file = readme_path, append = TRUE)
cat("|------|------|------|---------|\n", file = readme_path, append = TRUE)
cat("| Item_MELS_Noon_Main      | 午间 | 4 项 | mood_lag1_day_c + study_day_c + health_session_between_c |\n",
    file = readme_path, append = TRUE)
cat("| Item_MELS_Noon_Robust    | 午间 | 4 项 | 主控制 + mood_prev_session_c |\n",
    file = readme_path, append = TRUE)
cat("| Item_MELS_Evening_Main   | 晚间 | 4 项 | 同主控制集 |\n",
    file = readme_path, append = TRUE)
cat("| Item_MELS_Evening_Robust | 晚间 | 4 项 | 主控制 + mood_prev_session_c |\n\n",
    file = readme_path, append = TRUE)

cat(sprintf("**采样**: chains=%d × iter=%d (warmup=%d), adapt_delta=%.2f, max_treedepth=%d, seed=%d\n\n",
    n_chains, n_iter, n_warmup, mcmc_ctrl$adapt_delta, mcmc_ctrl$max_treedepth, GLOBAL_SEED),
    file = readme_path, append = TRUE)

cat("## 输出文件\n\n", file = readme_path, append = TRUE)
cat("| 文件 | 用途 |\n", file = readme_path, append = TRUE)
cat("|------|------|\n", file = readme_path, append = TRUE)
cat("| `fits/<model>.rds` | brmsfit 对象（xz 压缩） |\n",
    file = readme_path, append = TRUE)
cat("| `logs/<model>_cmdstan.log` | cmdstan 原始日志 |\n",
    file = readme_path, append = TRUE)
cat("| `posterior_summaries/<model>_posterior.csv` | 全参数后验摘要 |\n",
    file = readme_path, append = TRUE)
cat("| `item_mels_effects.csv` | 8 条目 × mu/sigma 效应汇总（论文主表来源） |\n",
    file = readme_path, append = TRUE)
cat("| `convergence_summary.csv` | 4 个模型的 max R-hat / min ESS / divergent |\n",
    file = readme_path, append = TRUE)
cat("| `figures/fig_item_mels_mu.png` | mu forest plot（mood 位置效应） |\n",
    file = readme_path, append = TRUE)
cat("| `figures/fig_item_mels_sigma.png` | sigma forest plot（mood 尺度效应） |\n\n",
    file = readme_path, append = TRUE)

cat("## 关键统计量来源映射（验证用）\n\n", file = readme_path, append = TRUE)
cat("所有报告的 mu/sigma 系数、CrI、R-hat、ESS 均来源于 `posterior_summaries/<model>_posterior.csv`\n",
    file = readme_path, append = TRUE)
cat("（由 `posterior::summarise_draws()` 直接生成，未硬编码）。\n\n",
    file = readme_path, append = TRUE)
cat("`item_mels_effects.csv` 中：\n", file = readme_path, append = TRUE)
cat("- `mu_median` 来自 `b_<item>` 行的 `median` 列\n",
    file = readme_path, append = TRUE)
cat("- `sig_log_median` 来自 `b_sigma_<item>` 行的 `median` 列\n",
    file = readme_path, append = TRUE)
cat("- `sig_mult_median = exp(sig_log_median)`，对应每单位条目带来的残差 SD 乘性变化\n",
    file = readme_path, append = TRUE)
cat("- `sig_compression_pct = (1 - sig_mult_median) * 100`，正值=压缩 mood 波动\n\n",
    file = readme_path, append = TRUE)

cat("---\n", file = readme_path, append = TRUE)
cat(sprintf("*生成于 %s by `run_item_mels.R`*\n",
            format(Sys.time(), "%Y-%m-%d %H:%M:%S")),
    file = readme_path, append = TRUE)

cat(sprintf("  ✓ README 写入: %s\n", readme_path))

# 完成标记（项目规则：worker 完成信号）
done_marker <- file.path(output_dir, ".done")
writeLines(format(Sys.time(), "%Y-%m-%d %H:%M:%S"), done_marker)

cat(sprintf("\n  ✓ 全部完成。完成时间: %s\n",
            format(Sys.time(), "%Y-%m-%d %H:%M:%S")))
cat(sprintf("\n  输出根目录: %s\n", output_dir))
