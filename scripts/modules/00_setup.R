# ============================================================================
# 模块 00b (): brms / cmdstanr 环境检测 + 输出目录 + 全局容器
# ============================================================================
# 依赖: 无（ 主线起点）
# 产出: 全局路径变量 (mini_dir, output_dir, fits_dir, ...)
#       brms 后端配置, 全局结果容器
# ============================================================================
# 设计决策（来自  计划）:
#   - 不自动安装 cmdstanr，缺失即报错并打印安装命令（规则9）
#   - 后端 cmdstanr，mc.cores = detectCores() - 1
#   - 4 个输出子目录：fits/posterior_summaries/figures/paper_tables/logs
# ============================================================================

cat("\n--- 模块00b (): brms / cmdstanr 环境与输出配置 ---\n")

# ============================================================================
# 0Ba. 包加载
# ============================================================================
suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(ggplot2)
})

# brms 与 cmdstanr 是核心依赖：缺失即报错
required_pkgs <- c("brms", "cmdstanr", "posterior", "bayesplot",
                   "loo", "tidybayes")
missing_pkgs <- required_pkgs[!sapply(required_pkgs, requireNamespace, quietly = TRUE)]
if (length(missing_pkgs) > 0) {
  cat("\n  ✗ 缺失以下 R 包：\n")
  for (p in missing_pkgs) cat(sprintf("      - %s\n", p))
  cat("\n  安装指引：\n")
  cat("    install.packages(c('brms', 'posterior', 'bayesplot', 'loo', 'tidybayes'))\n")
  cat("    install.packages('cmdstanr',\n")
  cat("      repos = c('https://stan-dev.r-universe.dev', getOption('repos')))\n")
  cat("    cmdstanr::check_cmdstan_toolchain(fix = TRUE)\n")
  cat("    cmdstanr::install_cmdstan()\n")
  stop(" 缺失必需包，已中止。手动安装后重跑。")
}

suppressPackageStartupMessages({
  library(brms)
  library(cmdstanr)
  library(posterior)
  library(bayesplot)
  library(loo)
})

# 验证 cmdstan 工具链
tryCatch({
  cmdstanr::check_cmdstan_toolchain(fix = FALSE, quiet = TRUE)
  cmdstan_ver <- cmdstanr::cmdstan_version()
  cat(sprintf("  ✓ cmdstan 已检测到，版本 %s\n", cmdstan_ver))
}, error = function(e) {
  cat("\n  ✗ cmdstan 未正确安装或工具链不完整\n")
  cat("    运行: cmdstanr::install_cmdstan()\n")
  stop(e$message)
})

# ============================================================================
# 0Bb. 路径配置（: 读路径指向 legacy/ 副本，写路径指向 results/）
# ============================================================================
if (!exists("SCRIPT_DIR")) SCRIPT_DIR <- getwd()
project_root <- normalizePath(
  file.path(SCRIPT_DIR, "..", ".."),   # scripts/pipeline → EMA_MELS/
  mustWork = FALSE
)

# ---- 数据准备回退模块路径 ----
module_dir  <- file.path(project_root, "scripts", "modules", "fallback")
data_dir    <- file.path(project_root, "data", "processed", "cleaned", "d7_s14")
results_dir <- file.path(project_root, "data", "checkpoints", "modeling_results")

# checkpoint 目录（单目录；历史多版本 checkpoint 已合并）
checkpoint_dir <- file.path(project_root, "data", "checkpoints")
if (!file.exists(file.path(checkpoint_dir, "checkpoint_01_data.RData"))) {
  checkpoint_dir <- file.path(project_root, "data", "checkpoints")
}

# ----   输出路径 ----
output_dir   <- file.path(project_root, "results")
fits_dir         <- file.path(output_dir, "fits")
posterior_dir    <- file.path(output_dir, "posterior_summaries")
figures_dir      <- file.path(output_dir, "figures")
tables_dir   <- file.path(output_dir, "tables")
logs_dir     <- file.path(output_dir, "logs")

for (d in c(output_dir, fits_dir, posterior_dir, figures_dir,
            tables_dir, logs_dir)) {
  dir.create(d, showWarnings = FALSE, recursive = TRUE)
}

cat(sprintf("  项目根:        %s\n", project_root))
cat(sprintf("   模块:     %s\n", module_dir))
cat(sprintf("   结果:     %s\n", results_dir))
cat(sprintf("   输出: %s\n", output_dir))

# ============================================================================
# 0Bc. brms 全局选项
# ============================================================================
options(
  brms.backend  = "cmdstanr",
  brms.file_refit = "on_change",
  mc.cores = max(1L, parallel::detectCores() - 1L)
)
cat(sprintf("  brms 后端:    %s\n", getOption("brms.backend")))
cat(sprintf("  并行核心数:   %d\n", getOption("mc.cores")))

# ============================================================================
# 0Bd. 全局结果容器
# ============================================================================
brms_fits <- list()
posterior_summaries <- list()
bayes_r2_results <- list()
convergence <- list()
mels_sigma_results <- list()
comparison_brms_vs_lme4 <- list()
ppcheck_results <- list()

# ============================================================================
# 0Be. 全局随机种子（可重现性）
# ============================================================================
GLOBAL_SEED <- 20260511L
set.seed(GLOBAL_SEED)
cat(sprintf("  全局种子:     %d\n", GLOBAL_SEED))

V914_TAG <- "_decomp_fix_L2_standardized"
cat(sprintf("  版本标记:     %s\n", V914_TAG))

# ============================================================================
# 0Bf. 选择性强制重拟合（环境变量控制）
# ============================================================================
# 用法: FORCE_REFIT="B1_MELS,B1_RS" Rscript ema_modeling__master.R
force_refit_env <- Sys.getenv("FORCE_REFIT", "")
FORCE_REFIT <- if (nchar(force_refit_env) > 0) {
  trimws(strsplit(force_refit_env, ",")[[1]])
} else {
  character(0)
}
if (length(FORCE_REFIT) > 0) {
  cat(sprintf("  强制重拟合:   %s\n", paste(FORCE_REFIT, collapse = ", ")))
}

# ============================================================================
# 0Bg. 记录 session_info
# ============================================================================
session_info_path <- file.path(logs_dir, "session_info_.txt")
sink(session_info_path)
cat("=====  Session Info =====\n")
cat(sprintf("Run timestamp: %s\n", format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z")))
cat(sprintf("R version: %s\n", R.version.string))
cat(sprintf("brms version: %s\n", as.character(packageVersion("brms"))))
cat(sprintf("cmdstanr version: %s\n", as.character(packageVersion("cmdstanr"))))
cat(sprintf("cmdstan version: %s\n", cmdstan_ver))
cat(sprintf("posterior version: %s\n", as.character(packageVersion("posterior"))))
cat("\n--- sessionInfo() ---\n")
print(sessionInfo())
sink()
cat(sprintf("  ✓ session_info 已记录: %s\n", session_info_path))

cat("\n  ✓ 模块00b 环境与输出配置完成\n")
