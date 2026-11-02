#!/usr/bin/env Rscript
options(width = 200)
options(warn = 1)

get_project_root <- function() {
  tryCatch({
    p <- sys.frame(1)$ofile
    if (!is.null(p)) return(normalizePath(file.path(dirname(p), "..", "..")))
  }, error = function(e) {})
  tryCatch({
    args <- commandArgs(trailingOnly = FALSE)
    f <- sub("^--file=", "", args[grep("^--file=", args)])
    if (length(f) > 0 && file.exists(f)) return(normalizePath(file.path(dirname(f), "..", "..")))
  }, error = function(e) {})
  getwd()
}
project_root <- get_project_root()
PIPELINE_DIR <- file.path(project_root, "scripts", "pipeline")
MODULE_DIR   <- file.path(project_root, "scripts", "modules")
RESULTS_DIR  <- file.path(project_root, "results")

results <- list()

cat("\n")
cat("╔", rep("═", 76), "╗\n", sep = "")
cat("║  EMA_MELS  — 全管线端到端测试                                 ║\n")
cat("╚", rep("═", 76), "╝\n", sep = "")
cat(sprintf("  项目根目录: %s\n", project_root))
cat(sprintf("  开始时间:   %s\n\n", format(Sys.time(), "%Y-%m-%d %H:%M:%S")))

# ============================================================================
# E2E-1: 主流水线 master.R (00b→08b)
# ============================================================================
cat("── E2E-1: 主流水线 master.R ──────────────────────────────\n")
e2e1_t0 <- Sys.time()
e2e1_status <- tryCatch({
  source(file.path(PIPELINE_DIR, "master.R"))
  "PASS"
}, error = function(e) {
  cat(sprintf("  ✗ FAIL: %s\n", conditionMessage(e)))
  "FAIL"
})
e2e1_elapsed <- as.numeric(difftime(Sys.time(), e2e1_t0, units = "mins"))
results[["E2E-1 master.R"]] <- list(status = e2e1_status, elapsed = e2e1_elapsed)
cat(sprintf("  → %s (%.1f min)\n\n", e2e1_status, e2e1_elapsed))

# ============================================================================
# E2E-2: 补充阶段 — master.R --stage 6-8 (06c/07c/08c，原 supplementary.R 功能)
# ============================================================================
cat("── E2E-2: 补充阶段 (master.R --stage 6-8) ────────────────\n")
e2e2_t0 <- Sys.time()
if (!exists("SCRIPT_DIR")) SCRIPT_DIR <- PIPELINE_DIR
if (!exists("MODULE_DIR")) MODULE_DIR <- file.path(project_root, "scripts", "modules")
e2e2_status <- tryCatch({
  # 
  # Stage 6 = 06c LOO, Stage 7 = 07b+07c tables, Stage 8 = 08b+08c figures
  capture.output(source(file.path(MODULE_DIR, "06_loo.R")), type = "message")
  capture.output(source(file.path(MODULE_DIR, "07_tables.R")), type = "message")
  capture.output(source(file.path(MODULE_DIR, "07_supp_tables.R")), type = "message")
  capture.output(source(file.path(MODULE_DIR, "08_figures_main.R")), type = "message")
  capture.output(source(file.path(MODULE_DIR, "08_figures_paper.R")), type = "message")
  capture.output(source(file.path(MODULE_DIR, "08_figures_student_zh.R")), type = "message")
  "PASS"
}, error = function(e) {
  cat(sprintf("  ✗ FAIL: %s\n", conditionMessage(e)))
  "FAIL"
})
e2e2_elapsed <- as.numeric(difftime(Sys.time(), e2e2_t0, units = "mins"))
results[["E2E-2 suppl. pipeline (06c-08c)"]] <- list(status = e2e2_status, elapsed = e2e2_elapsed)
cat(sprintf("  → %s (%.1f min)\n\n", e2e2_status, e2e2_elapsed))

# ============================================================================
# E2E-2b: RS 六分量方差分解 06d
# ============================================================================
cat("── E2E-2b: 06d RS 方差分解 ──────────────────────────────\n")
e2e2b_t0 <- Sys.time()
e2e2b_status <- tryCatch({
  conv_csv <- file.path(project_root, "results", "convergence_diagnostics_.csv")
  if (file.exists(conv_csv)) {
    convergence <- read.csv(conv_csv, stringsAsFactors = FALSE)
    convergence <- split(convergence, convergence$model)
  } else {
    convergence <- list()
  }
  source(file.path(MODULE_DIR, "06_rs_decomp.R"))
  "PASS"
}, error = function(e) {
  cat(sprintf("  ✗ FAIL: %s\n", conditionMessage(e)))
  "FAIL"
})
e2e2b_elapsed <- as.numeric(difftime(Sys.time(), e2e2b_t0, units = "mins"))
results[["E2E-2b 06d RS decomp"]] <- list(status = e2e2b_status, elapsed = e2e2b_elapsed)
cat(sprintf("  → %s (%.1f min)\n\n", e2e2b_status, e2e2b_elapsed))

# ============================================================================
# E2E-3: 敏感性脚本路径验证（快速模式：仅 source 00b/01b/02b）
# ============================================================================
cat("── E2E-3: 敏感性脚本路径验证 ───────────────────────────\n")

sensitivity_scripts <- c(
  "run_bayes_factor.R",
  "run_prior_sensitivity.R",
  "run_mi_sensitivity.R",
  "run_conditional_effects.R",
  "run_sesoi_posterior.R",
  "run_bf_prior_sensitivity.R",
  "run_plan_ab_trace.R"
)

sensitivity_dir <- file.path(project_root, "scripts", "sensitivity")

for (script in sensitivity_scripts) {
  script_path <- file.path(sensitivity_dir, script)
  if (!file.exists(script_path)) {
    results[[paste0("E2E-3 ", script)]] <- list(status = "SKIP (not found)", elapsed = 0)
    cat(sprintf("  → SKIP %s (文件不存在)\n", script))
    next
  }

  e2e3_status <- tryCatch({
    parse(file = script_path)

    lines <- readLines(script_path)
    mod_ref <- grep("MODULE_DIR.*file\\.path", lines, value = TRUE)
    fit_ref <- grep("(fits_dir|figures_dir|tables_dir|checkpoint_dir|output_dir)", lines, value = TRUE)

    has_old_path <- any(grepl("EMA_habit", lines))
    if (has_old_path) {
      cat(sprintf("  ✗ %s: 仍含硬编码旧路径\n", script))
      "FAIL"
    } else {
      "PASS"
    }
  }, error = function(e) {
    cat(sprintf("  ✗ %s: 语法错误 — %s\n", script, conditionMessage(e)))
    "FAIL"
  })

  results[[paste0("E2E-3 ", script)]] <- list(status = e2e3_status, elapsed = 0)
  cat(sprintf("  → %s %s\n", e2e3_status, script))
}
cat("\n")

# ============================================================================
# E2E-4: 产出清单验证
# ============================================================================
cat("── E2E-4: 产出清单验证 ────────────────────────────────\n")

check_artifacts <- function() {
  fails <- 0
  passes <- 0

  fits_dir <- file.path(RESULTS_DIR, "fits")
  n_fits <- length(list.files(fits_dir, pattern = "\\.rds$", recursive = TRUE))
  if (n_fits >= 10) {
    cat(sprintf("  ✓ results/fits/: %d .rds 文件\n", n_fits))
    passes <- passes + 1
  } else {
    cat(sprintf("  ✗ results/fits/: 仅 %d .rds 文件（预期 >= 10）\n", n_fits))
    fails <- fails + 1
  }

  tables_dir <- file.path(RESULTS_DIR, "tables")
  n_csv <- length(list.files(tables_dir, pattern = "\\.csv$"))
  n_md  <- length(list.files(tables_dir, pattern = "\\.md$"))
  if (n_csv >= 5 && n_md >= 3) {
    cat(sprintf("  ✓ results/tables/: %d CSV + %d MD\n", n_csv, n_md))
    passes <- passes + 1
  } else {
    cat(sprintf("  ✗ results/tables/: %d CSV + %d MD（预期 >= 5 CSV, >= 3 MD）\n", n_csv, n_md))
    fails <- fails + 1
  }

  figures_dir <- file.path(RESULTS_DIR, "figures")
  n_png <- length(list.files(figures_dir, pattern = "\\.png$", recursive = TRUE))
  if (n_png >= 30) {
    cat(sprintf("  ✓ results/figures/: %d PNG 文件\n", n_png))
    passes <- passes + 1
  } else {
    cat(sprintf("  ✗ results/figures/: 仅 %d PNG 文件（预期 >= 30）\n", n_png))
    fails <- fails + 1
  }

  posterior_dir <- file.path(RESULTS_DIR, "posterior_summaries")
  n_post <- length(list.files(posterior_dir, pattern = "\\.csv$"))
  if (n_post >= 5) {
    cat(sprintf("  ✓ results/posterior_summaries/: %d CSV\n", n_post))
    passes <- passes + 1
  } else {
    cat(sprintf("  ✗ results/posterior_summaries/: 仅 %d CSV（预期 >= 5）\n", n_post))
    fails <- fails + 1
  }

  cat("\n  ──  特定验证 ──\n")

  conv_path <- file.path(tables_dir, "convergence_diagnostics_.csv")
  if (file.exists(conv_path)) {
    conv_df <- read.csv(conv_path, stringsAsFactors = FALSE)
    n_pass_conv <- sum(conv_df$status == "PASS")
    n_warn_conv <- sum(conv_df$status == "WARN")
    if (n_pass_conv == 6 && n_warn_conv == 0) {
      cat(sprintf("  ✓  收敛: %d PASS / %d WARN\n", n_pass_conv, n_warn_conv))
      passes <- passes + 1
    } else {
      cat(sprintf("  ✗  收敛: %d PASS / %d WARN（预期 6/0）\n", n_pass_conv, n_warn_conv))
      fails <- fails + 1
    }
  } else {
    cat("  ✗  收敛 CSV 不存在\n")
    fails <- fails + 1
  }

  rs_decomp_path <- file.path(tables_dir, "bayes_r2_decomposition_rs.csv")
  if (file.exists(rs_decomp_path)) {
    rs_df <- read.csv(rs_decomp_path, stringsAsFactors = FALSE)
    rs_models <- unique(rs_df$model)
    sum_ok <- TRUE
    for (m in rs_models) {
      sub <- rs_df[rs_df$model == m, ]
      comp_sum <- sum(sub$median[sub$component != "fvm_total"])
      if (abs(comp_sum - 1.0) > 0.01) {
        cat(sprintf("  ✗  RS 方差分解 %s: sum = %.4f (预期 ≈1.0)\n", m, comp_sum))
        sum_ok <- FALSE
      }
    }
    if (sum_ok) {
      cat(sprintf("  ✓  RS 方差分解: %d 模型, sum ≈ 1.0\n", length(rs_models)))
      passes <- passes + 1
    } else {
      fails <- fails + 1
    }
  } else {
    cat("  ✗  RS 方差分解 CSV 不存在\n")
    fails <- fails + 1
  }

  list(passes = passes, fails = fails)
}

artifacts <- check_artifacts()
results[["E2E-4 artifacts"]] <- list(
  status = if (artifacts$fails == 0) "PASS" else paste0("FAIL (", artifacts$fails, " checks)"),
  elapsed = 0
)
cat("\n")

# ============================================================================
# 汇总报告
# ============================================================================
cat("╔", rep("═", 76), "╗\n", sep = "")
cat("║  端到端测试汇总                                                    ║\n")
cat("╚", rep("═", 76), "╝\n", sep = "")

n_pass <- 0
n_fail <- 0
n_skip <- 0

for (name in names(results)) {
  r <- results[[name]]
  status_mark <- if (grepl("^PASS", r$status)) {
    n_pass <<- n_pass + 1; "✓"
  } else if (grepl("^FAIL", r$status)) {
    n_fail <<- n_fail + 1; "✗"
  } else {
    n_skip <<- n_skip + 1; "—"
  }
  cat(sprintf("  %s %-50s %s\n", status_mark, name, r$status))
  if (r$elapsed > 0) {
    cat(sprintf("    (%.1f min)\n", r$elapsed))
  }
}

cat("\n")
cat(sprintf("  总计: %d PASS / %d FAIL / %d SKIP\n", n_pass, n_fail, n_skip))
cat(sprintf("  结束时间: %s\n", format(Sys.time(), "%Y-%m-%d %H:%M:%S")))

if (n_fail > 0) {
  cat("\n  ⚠ 存在失败项，请检查上述 FAIL 输出。\n")
  quit(status = 1)
} else {
  cat("\n  ✓ 全部测试通过！\n")
}
