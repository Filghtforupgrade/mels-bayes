#!/usr/bin/env Rscript
cat("=== EMA_MELS 冒烟测试 ===\n")
cat("时间:", format(Sys.time()), "\n\n")

# 1. 检查R版本
cat("1. R版本检查\n")
cat("R版本:", R.version.string, "\n")
cat("平台:", R.version$platform, "\n\n")

# 2. 检查关键包
cat("2. 关键包检查\n")
required_packages <- c("brms", "cmdstanr", "bayesplot", "tidybayes", "ggplot2", "dplyr", "lme4", "jsonlite")
missing_packages <- c()

for (pkg in required_packages) {
  if (requireNamespace(pkg, quietly = TRUE)) {
    pkg_version <- packageVersion(pkg)
    cat(sprintf("  ✓ %s (版本: %s)\n", pkg, pkg_version))
  } else {
    cat(sprintf("  ✗ %s (未安装)\n", pkg))
    missing_packages <- c(missing_packages, pkg)
  }
}

if (length(missing_packages) > 0) {
  cat("\n缺少以下包:", paste(missing_packages, collapse = ", "), "\n")
} else {
  cat("\n所有关键包已安装\n")
}

# 3. 检查cmdstanr配置
cat("\n3. cmdstanr配置检查\n")
if (requireNamespace("cmdstanr", quietly = TRUE)) {
  tryCatch({
    cmdstan_path <- cmdstanr::cmdstan_path()
    cat("  ✓ cmdstan路径:", cmdstan_path, "\n")
    cat("  ✓ cmdstan版本:", cmdstanr::cmdstan_version(), "\n")
  }, error = function(e) {
    cat("  ✗ cmdstanr配置错误:", conditionMessage(e), "\n")
  })
}

# 4. 检查项目结构
cat("\n4. 项目结构检查\n")

# 确定项目根目录（不依赖here包）
get_project_root <- function() {
  # 尝试从脚本位置推导
  tryCatch({
    script_path <- sys.frame(1)$ofile
    if (!is.null(script_path)) {
      script_dir <- normalizePath(dirname(script_path))
      # 从scripts/validation向上两级到项目根目录
      return(file.path(script_dir, "..", ".."))
    }
  }, error = function(e) {})
  
  # 回退到命令行参数
  tryCatch({
    args <- commandArgs(trailingOnly = FALSE)
    file_arg <- grep("^--file=", args, value = TRUE)
    if (length(file_arg) > 0) {
      script_path <- sub("^--file=", "", file_arg[1])
      script_dir <- normalizePath(dirname(script_path))
      return(file.path(script_dir, "..", ".."))
    }
  }, error = function(e) {})
  
  # 最后回退到当前工作目录
  return(getwd())
}

project_root <- normalizePath(get_project_root())
cat("项目根目录:", project_root, "\n")

# 检查关键目录
key_dirs <- c("scripts/modules", "data/processed", "results/fits", "results/logs")
for (dir in key_dirs) {
  full_path <- file.path(project_root, dir)
  if (dir.exists(full_path)) {
    file_count <- length(list.files(full_path, recursive = TRUE))
    cat(sprintf("  ✓ %s (文件数: %d)\n", dir, file_count))
  } else {
    cat(sprintf("  ✗ %s (目录不存在)\n", dir))
  }
}

# 5. 尝试加载核心模块
cat("\n5. 核心模块加载测试\n")
module_files <- c("00_setup.R", "02_priors.R")
modules_dir <- file.path(project_root, "scripts/modules")

# 设置SCRIPT_DIR为scripts/pipeline，这是00_setup.R期望的位置
SCRIPT_DIR <- file.path(project_root, "scripts", "pipeline")

for (module in module_files) {
  module_path <- file.path(modules_dir, module)
  if (file.exists(module_path)) {
    tryCatch({
      source(module_path, local = TRUE)
      cat(sprintf("  ✓ %s 加载成功\n", module))
    }, error = function(e) {
      cat(sprintf("  ✗ %s 加载失败: %s\n", module, conditionMessage(e)))
    })
  } else {
    cat(sprintf("  ✗ %s 文件不存在\n", module))
  }
}

# 6. 检查审计日志目录
cat("\n6. 审计日志检查\n")
logs_dir <- file.path(project_root, "results/logs")
if (dir.exists(logs_dir)) {
  log_files <- list.files(logs_dir, pattern = "\\.jsonl$")
  cat(sprintf("  ✓ 日志目录存在，包含 %d 个JSONL文件\n", length(log_files)))
  if (length(log_files) > 0) {
    cat("  最新日志文件:\n")
    for (log_file in head(log_files, 3)) {
      cat(sprintf("    - %s\n", log_file))
    }
  }
} else {
  cat("  ✗ 日志目录不存在\n")
}

cat("\n=== 冒烟测试完成 ===\n")
cat("如果所有检查都通过，说明环境基本正常。\n")
cat("如果需要完整测试，请运行主分析脚本。\n")
