# ============================================================================
# 模块 00: 基础设施 — 包加载 + 路径配置 + 全局结果容器
# ============================================================================
# 依赖: 无
# 产出: project_root, data_dir, output_dir, table_dir, 全部结果容器
# ============================================================================

cat("\n--- 模块00: 包加载与路径配置 ---\n")

# ============================================================================
# 0A. 包加载
# ============================================================================
suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(lme4)
  library(lmerTest)
  library(car)        # 用于VIF计算
  library(performance) # 用于ICC计算
  library(ggplot2)    # 用于可视化
  library(r2mlm)       # 多层模型R²方差分解（Shaw et al., 2023, BRM）
})

suppressPackageStartupMessages({
  tryCatch({
    library(simr)
    cat("  ✓ simr包已加载 (事后功效分析)\n")
  }, error = function(e) {
    cat("  ⚠ simr包未安装，功效分析将使用解析公式替代\n")
    cat("    安装命令: install.packages('simr')\n")
  })
})

suppressPackageStartupMessages({
  tryCatch({ library(modelsummary); cat("  ✓ modelsummary包已加载 (论文表格)\n") },
    error = function(e) cat("  ⚠ modelsummary未安装, install.packages('modelsummary')\n"))
  tryCatch({ library(broom.mixed); cat("  ✓ broom.mixed包已加载 (lmer支持)\n") },
    error = function(e) cat("  ⚠ broom.mixed未安装, install.packages('broom.mixed')\n"))
  tryCatch({ library(gtsummary); cat("  ✓ gtsummary包已加载 (基线表/回归表)\n") },
    error = function(e) cat("  ⚠ gtsummary未安装, install.packages('gtsummary')\n"))
  tryCatch({ library(sjPlot); cat("  ✓ sjPlot包已加载 (APA多层模型表)\n") },
    error = function(e) cat("  ⚠ sjPlot未安装, install.packages('sjPlot')\n"))
  tryCatch({ library(flextable); cat("  ✓ flextable包已加载 (Word表格引擎)\n") },
    error = function(e) cat("  ⚠ flextable未安装, install.packages('flextable')\n"))
  tryCatch({ library(gt); cat("  ✓ gt包已加载 (HTML表格引擎)\n") },
    error = function(e) cat("  ⚠ gt未安装, install.packages('gt')\n"))
  tryCatch({ library(kableExtra); cat("  ✓ kableExtra包已加载 (LaTeX表格引擎)\n") },
    error = function(e) cat("  ⚠ kableExtra未安装, install.packages('kableExtra')\n"))
})

# ============================================================================
# 0B. 路径配置
# ============================================================================
# 从master脚本的SCRIPT_DIR向上推算项目根目录
# master: project_root/code_archive/mini/ema_modeling_v7_3_master.R
# module: project_root/code_archive/mini/modules/00_config.R
# 使用master中已定义的 SCRIPT_DIR (= mini/) 进行计算
if (!exists("SCRIPT_DIR")) {
  SCRIPT_DIR <- getwd()
}
mini_dir <- normalizePath(SCRIPT_DIR, mustWork = FALSE)        # scripts/modules/fallback/
project_root <- normalizePath(file.path(mini_dir, "..", ".."), mustWork = FALSE)  # EMA_MELS/
module_dir <- file.path(mini_dir, "modules")                    # scripts/modules/fallback/

# 敏感性分析: 支持环境变量 EMA_DATA_DIR 动态注入数据目录
# 用法: EMA_DATA_DIR="cleaned_data_d5_s10" Rscript ema_modeling_master.R
# 主分析默认: data/processed/cleaned/d7_s14（≥7天 / ≥14 sessions，236人）
env_data_dir <- Sys.getenv("EMA_DATA_DIR", unset = "")
if (nchar(env_data_dir) > 0) {
  data_dir <- ifelse(
    startsWith(env_data_dir, "/"),
    env_data_dir,
    file.path(project_root, env_data_dir)
  )
  cat(sprintf("  ★ 环境变量注入数据目录: %s\n", data_dir))
} else {
  data_dir <- file.path(project_root, "data", "processed", "cleaned", "d7_s14")
}
output_dir <- file.path(project_root, "data", "checkpoints", "modeling_results")
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

# 敏感性分析: 支持环境变量 EMA_OUTPUT_DIR 动态注入输出目录
env_output_dir <- Sys.getenv("EMA_OUTPUT_DIR", unset = "")
if (nchar(env_output_dir) > 0) {
  output_dir <- ifelse(
    startsWith(env_output_dir, "/"),
    env_output_dir,
    file.path(project_root, env_output_dir)
  )
  dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)
  cat(sprintf("  ★ 环境变量注入输出目录: %s\n", output_dir))
}

table_dir <- file.path(output_dir, "paper_tables")
dir.create(table_dir, showWarnings = FALSE, recursive = TRUE)

checkpoint_dir <- file.path(project_root, "data", "checkpoints")
dir.create(checkpoint_dir, showWarnings = FALSE, recursive = TRUE)

cat(sprintf("  项目根目录: %s\n", project_root))
cat(sprintf("  数据目录:   %s\n", data_dir))
cat(sprintf("  输出目录:   %s\n", output_dir))
cat(sprintf("  检查点目录: %s\n", checkpoint_dir))

# ============================================================================
# 0C. 全局结果容器
# ============================================================================
all_results <- list()
all_diagnostics <- list()
all_icc_results <- list()
all_individual_slopes <- data.frame()
all_delta_r2_results <- data.frame()
all_r2mlm_results <- list()
all_robust_tests <- list()

# : 新增结果容器
all_convergence_diag <- list()
all_power_results <- list()
all_missing_patterns <- list()
all_r2mlm_verification <- list()
all_b4_full_results <- list()
all_model_store <- list()

cat("\n  ✓ 模块00 基础设施初始化完成\n")
cat(sprintf("    输出目录: %s\n", output_dir))
cat(sprintf("    表格目录: %s\n", table_dir))
