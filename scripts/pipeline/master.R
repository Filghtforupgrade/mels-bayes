#!/usr/bin/env Rscript
# -*- coding: utf-8 -*-
# ============================================================================
# EMA  贝叶斯分析 — 统一主控脚本
# ============================================================================
# Revision tag: v10.0.0 (unified pipeline, 2026-06-02)
#
# 目的: 单一入口覆盖全部分析管线——核心建模、补充材料、敏感性分析、端到端验证。
#
# 模块列表 (scripts/modules/):
#   00_setup.R               — cmdstanr/brms 环境检测 + 输出目录 + 容器初始化
#   01_data.R              — 优先复用 checkpoint；不变量断言
#   02_priors.R           — 三套先验 + fit_or_load + bayes_r2_decomp + compare_to_lme4
#   03_mels.R              — B1 / R3 MELS（2 个 brms 拟合）
#   04_rs_models.R          — 4 个贝叶斯随机斜率（B1_RS + 3 跨层调节）
#   05_diagnostics.R              — R-hat/ESS/divergent + pp_check (PASS/WARN/FAIL)
#   06_r2_decomp.R          — bayes_R2 + 五分量方差分解
#   06_rs_decomp.R  — RS 模型六分量方差分解
#   06_loo.R           — LOO 模型比较
#   07_tables.R           — 3 张 MD 表 + brms vs lme4 跨范式对照
#   07_supp_tables.R — 6 张补充表格（Dora 2024 对标）
#   08_figures_main.R            — 6 张主图
#   08_figures_paper.R       — 4 张论文图（Dora 2024 对标）
#   08_figures_student_zh.R     — 学生演讲中文图
#
# 敏感性分析 (scripts/sensitivity/):
#   run_prior_sensitivity.R   — 3 档先验 × 2 模型
#   run_bayes_factor.R        — BF (tight prior)
#   run_bf_prior_sensitivity.R — BF (wide prior)
#   run_mi_sensitivity.R      — MI M=10
#   run_conditional_effects.R — 条件效应图
#   run_sesoi_posterior.R     — SESOI 后验概率
#   run_plan_ab_trace.R — Plan A vs B trace
#
# 用法:
#   Rscript master.R                           # 全量运行 (Stage 0-10)
#   Rscript master.R --stage 0-8               # 仅核心管线 (到论文图)
#   Rscript master.R --stage 9                 # 仅敏感性分析
#   Rscript master.R --stage 0-4               # 仅建模+诊断
#   Rscript master.R --dry-run                 # 干跑模式
#
# 环境变量:
#   FORCE_REFIT="B1_MELS,B1_RS"            # 强制重拟合
#   DRY_RUN=TRUE                                # 干跑模式
# ============================================================================

options(width = 200)
options(warn = 1)

# ----- 命令行参数解析 -----
args <- commandArgs(trailingOnly = TRUE)
stage_start <- 0L
stage_end   <- 10L
dry_run     <- FALSE

i <- 1L
while (i <= length(args)) {
  if (args[i] == "--stage") {
    if (i + 1L <= length(args)) {
      parts <- strsplit(args[i + 1L], "-")[[1]]
      if (length(parts) == 2L) {
        stage_start <- as.integer(parts[1])
        stage_end   <- as.integer(parts[2])
      } else if (length(parts) == 1L) {
        stage_start <- as.integer(parts[1])
        stage_end   <- stage_start
      }
      i <- i + 1L
    }
  } else if (args[i] == "--dry-run") {
    dry_run <- TRUE
  }
  i <- i + 1L
}

if (dry_run || isTRUE(as.logical(Sys.getenv("DRY_RUN", "FALSE")))) {
  dry_run <- TRUE
  cat("\n  ⚠ DRY-RUN MODE: using chains=2, iter=500\n")
}

# ----- 获取脚本所在目录 -----
current_file <- commandArgs(trailingOnly = FALSE)
script_path <- sub("^--file=", "", current_file[grep("^--file=", current_file)])
if (length(script_path) > 0 && file.exists(script_path)) {
  SCRIPT_DIR <- dirname(script_path)
} else {
  SCRIPT_DIR <- getwd()
}
MODULE_DIR <- file.path(SCRIPT_DIR, "..", "modules")
SENSITIVITY_DIR <- file.path(SCRIPT_DIR, "..", "sensitivity")
VALIDATION_DIR  <- file.path(SCRIPT_DIR, "..", "validation")

cat("\n")
cat("╔", rep("═", 78), "╗\n", sep = "")
cat("║  EMA  贝叶斯分析 — 统一主控脚本                                ║\n")
cat(sprintf("║  阶段范围: %d–%d%-50s║\n", stage_start, stage_end, if (dry_run) " [DRY-RUN]" else ""))
cat("╚", rep("═", 78), "╝\n", sep = "")

v10_t0 <- Sys.time()

# ============================================================================
# Stage 0: 基础设施层
# ============================================================================
run_stage <- function(num, label) {
  cat(sprintf("\n╔══ Stage %d: %s %s╗\n", num, label,
              paste(rep("═", max(0, 70 - nchar(label))), collapse = "")))
}

if (stage_start <= 0L && stage_end >= 0L) {
  run_stage(0L, "基础设施层")
  source(file.path(MODULE_DIR, "00_setup.R"))
  cat("  ✓ 00_setup.R\n")
}

# ============================================================================
# Stage 1: 数据准备层
# ============================================================================
if (stage_start <= 1L && stage_end >= 1L) {
  run_stage(1L, "数据准备层（复用 checkpoint）")
  source(file.path(MODULE_DIR, "01_data.R"))
  cat("  ✓ 01_data.R\n")
}

# ============================================================================
# Stage 2: 先验与辅助函数层
# ============================================================================
if (stage_start <= 2L && stage_end >= 2L) {
  run_stage(2L, "先验与辅助函数层")
  source(file.path(MODULE_DIR, "02_priors.R"))
  cat("  ✓ 02_priors.R\n")
}

# ============================================================================
# Stage 3: 贝叶斯建模层
# ============================================================================
if (stage_start <= 3L && stage_end >= 3L) {
  run_stage(3L, "贝叶斯建模层 — 6 个 brms 拟合")

  cat("\n--- 子阶段 3a: 2 个 MELS 模型 (B1 + R3) ---\n")
  source(file.path(MODULE_DIR, "03_mels.R"))
  cat("  ✓ 03_mels.R\n")

  cat("\n--- 子阶段 3b: 4 个贝叶斯随机斜率模型 (B1_RS + 3 调节版) ---\n")
  source(file.path(MODULE_DIR, "04_rs_models.R"))
  cat("  ✓ 04_rs_models.R\n")
}

# ============================================================================
# Stage 4: 收敛诊断层
# ============================================================================
if (stage_start <= 4L && stage_end >= 4L) {
  run_stage(4L, "收敛诊断层")
  source(file.path(MODULE_DIR, "05_diagnostics.R"))
  cat("  ✓ 05_diagnostics.R\n")
}

# ============================================================================
# Stage 5: 方差分解层
# ============================================================================
if (stage_start <= 5L && stage_end >= 5L) {
  run_stage(5L, "方差分解层")

  # 5a: 五分量分解 + bayes_R2（随机截距模型）
  cat("\n--- 子阶段 5a: bayes_R2 + 五分量方差分解 (06b) ---\n")
  source(file.path(MODULE_DIR, "06_r2_decomp.R"))
  cat("  ✓ 06_r2_decomp.R\n")

  # 5b: 六分量分解（随机斜率模型）
  cat("\n--- 子阶段 5b: RS 模型六分量方差分解 (06d) ---\n")
  source(file.path(MODULE_DIR, "06_rs_decomp.R"))
  cat("  ✓ 06_rs_decomp.R\n")
}

# ============================================================================
# Stage 6: LOO 模型比较
# ============================================================================
if (stage_start <= 6L && stage_end >= 6L) {
  run_stage(6L, "LOO 模型比较")
  source(file.path(MODULE_DIR, "06_loo.R"))
  cat("  ✓ 06_loo.R\n")
}

# ============================================================================
# Stage 7: 表格输出层
# ============================================================================
if (stage_start <= 7L && stage_end >= 7L) {
  run_stage(7L, "表格输出层")

  # 7a: 3 张主表 + 跨范式对照
  cat("\n--- 子阶段 7a: 主表格 + brms vs lme4 跨范式对照 (07b) ---\n")
  source(file.path(MODULE_DIR, "07_tables.R"))
  cat("  ✓ 07_tables.R\n")

  # 7b: 6 张补充表格
  cat("\n--- 子阶段 7b: 6 张补充表格 (07c) ---\n")
  source(file.path(MODULE_DIR, "07_supp_tables.R"))
  cat("  ✓ 07_supp_tables.R\n")
}

# ============================================================================
# Stage 8: 可视化层
# ============================================================================
if (stage_start <= 8L && stage_end >= 8L) {
  run_stage(8L, "可视化层")

  # 8a: 6 张主图
  cat("\n--- 子阶段 8a: 6 张主图 (08b) ---\n")
  source(file.path(MODULE_DIR, "08_figures_main.R"))
  cat("  ✓ 08_figures_main.R\n")

  # 8b: 4 张论文图
  cat("\n--- 子阶段 8b: 4 张论文图 (08c) ---\n")
  source(file.path(MODULE_DIR, "08_figures_paper.R"))
  cat("  ✓ 08_figures_paper.R\n")

  # 8c: 学生演讲中文图
  cat("\n--- 子阶段 8c: 学生演讲中文图 (08c student_zh) ---\n")
  source(file.path(MODULE_DIR, "08_figures_student_zh.R"))
  cat("  ✓ 08_figures_student_zh.R\n")
}

# ============================================================================
# Stage 9: 敏感性分析套件
# ============================================================================
if (stage_start <= 9L && stage_end >= 9L) {
  run_stage(9L, "敏感性分析套件（7 项）")

  # 确保主模型 fits 已加载（如跳过 Stage 3 则从 .rds 加载）
  if (!exists("brms_fits") || length(brms_fits) == 0L) {
    cat("\n  加载主模型 .rds 到 brms_fits 容器...\n")
    brms_fits <- list()
    main_models <- c("B1_MELS", "R3_MELS", "B1_RS",
                     "B1_RS_PHQ", "B1_RS_PC", "B1_RS_sigmaM")
    for (m in main_models) {
      rds_path <- file.path(fits_dir, paste0(m, ".rds"))
      if (file.exists(rds_path)) {
        brms_fits[[m]] <- readRDS(rds_path)
        cat(sprintf("    ✓ %s\n", m))
      }
    }
  }

  # 9a: 先验敏感性 — 3 档先验 × 2 模型
  cat("\n--- 子阶段 9a: 先验敏感性 (tight/wide/flat) ---\n")
  source(file.path(SENSITIVITY_DIR, "run_prior_sensitivity.R"))

  # 9b: Bayes Factor (tight prior)
  cat("\n--- 子阶段 9b: Bayes Factor (tight prior) ---\n")
  source(file.path(SENSITIVITY_DIR, "run_bayes_factor.R"))

  # 9c: Bayes Factor 先验敏感性 (wide prior)
  cat("\n--- 子阶段 9c: Bayes Factor 先验敏感性 (wide prior) ---\n")
  source(file.path(SENSITIVITY_DIR, "run_bf_prior_sensitivity.R"))

  # 9d: 多重插补稳健性
  cat("\n--- 子阶段 9d: 多重插补稳健性 (M=10) ---\n")
  source(file.path(SENSITIVITY_DIR, "run_mi_sensitivity.R"))

  # 9e: 条件效应图
  cat("\n--- 子阶段 9e: 条件效应图 ---\n")
  source(file.path(SENSITIVITY_DIR, "run_conditional_effects.R"))

  # 9f: SESOI 后验概率
  cat("\n--- 子阶段 9f: SESOI 后验概率 ---\n")
  source(file.path(SENSITIVITY_DIR, "run_sesoi_posterior.R"))

  # 9g: Plan A vs B trace
  cat("\n--- 子阶段 9g: Plan A vs B trace 对比 ---\n")
  source(file.path(SENSITIVITY_DIR, "run_plan_ab_trace.R"))

  cat("\n  ✓ 全部 7 项敏感性分析完成\n")
}

# ============================================================================
# Stage 10: 端到端验证
# ============================================================================
if (stage_start <= 10L && stage_end >= 10L) {
  run_stage(10L, "端到端验证")
  source(file.path(VALIDATION_DIR, "e2e_test.R"))
  cat("  ✓ e2e_test.R\n")
}

# ============================================================================
# 完成 + 总结
# ============================================================================
v10_elapsed <- difftime(Sys.time(), v10_t0, units = "mins")

cat("\n")
cat("╔", rep("═", 78), "╗\n", sep = "")
cat("║   统一主控执行完毕！                                            ║\n")
cat("╚", rep("═", 78), "╝\n", sep = "")
cat(sprintf("\n  总耗时: %.1f 分钟\n", as.numeric(v10_elapsed)))

if (exists("brms_fits")) {
  cat(sprintf("  brms 拟合数: %d\n", length(brms_fits)))
}
if (exists("output_dir")) {
  cat(sprintf("  输出目录: %s\n", output_dir))
  if (exists("fits_dir") && dir.exists(fits_dir)) {
    cat(sprintf("    fits/             %d 个 .rds\n",
                length(list.files(fits_dir, pattern = "\\.rds$"))))
  }
  if (exists("posterior_dir") && dir.exists(posterior_dir)) {
    cat(sprintf("    posterior_summaries/  %d 个 CSV\n",
                length(list.files(posterior_dir, pattern = "\\.csv$"))))
  }
  if (exists("figures_dir") && dir.exists(figures_dir)) {
    cat(sprintf("    figures/          %d 个 PNG\n",
                length(list.files(figures_dir, pattern = "\\.png$", recursive = TRUE))))
  }
  if (exists("tables_dir") && dir.exists(tables_dir)) {
    cat(sprintf("    tables/           %d 个文件\n",
                length(list.files(tables_dir))))
  }
}

# 写 session info
if (exists("logs_dir")) {
  sink(file.path(logs_dir, "session_info_v10.txt"))
  cat(sprintf(" 运行时间: %s\n", Sys.time()))
  cat(sprintf("阶段范围: %d–%d\n", stage_start, stage_end))
  cat(sprintf("干跑模式: %s\n", if (dry_run) "TRUE" else "FALSE"))
  cat(sprintf("总耗时: %.1f 分钟\n\n", as.numeric(v10_elapsed)))
  print(sessionInfo())
  sink()
  cat("  ✓ session_info 已写入 logs/\n")
}

cat("\n  下一步: 检查 paper_tables/ 中的 MD 表 + results/ 中的 CSV/PNG 产出\n\n")
