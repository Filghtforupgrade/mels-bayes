# ============================================================================
# 模块 06c (): LOO 模型比较 — 信息准则视角
# ============================================================================
# 依赖: 03b/04b 拟合完成（brms_fits 包含 6 个主模型）
#       可选：item-level fits（如已加载到 brms_fits 或 fits 容器中）
# 产出: 04_results/loo_results.csv         （各模型 elpd_loo 等）
#       04_results/loo_comparison.csv      （loo_compare 结果）
#       04_results/loo_pareto_k.csv        （高 k 观测预警）
#       04_results/logs/loo_warnings.txt  （Pareto k > 0.7 详细列表）
# ============================================================================
# 设计决策（ plan 阶段 4）:
#   - 对所有 PASS/WARN 主模型计算 loo()
#   - 关键比较对：
#       (1) B1_MELS vs B1_RS          — MELS 双子方程 vs RS 单方程
#       (2) R3_MELS vs B1_MELS         — 加 prev_mood × health 交互
#       (3) B1_RS_PHQ / PC / sigmaM    — 跨层调节对预测的提升
#   - 对每个 loo() 检查 Pareto k：
#       k > 0.7 → 自动用 moment_match=TRUE 重算
#       仍 > 0.7 → 记录到 warning 文件，注明该模型 elpd 估计不可靠
#   - 不要求 item-level 与主线模型比较（数据不同，elpd 无可比性）
# 注意:
#   - 6 个主模型 loo 约 1-3 分钟（每个 fit 24-48k draws，2.4k obs）
#   - 高 k 在 MELS 模型中常见（sigma 子方程让 cross-validation 难度提升）
# ============================================================================

cat("\n--- 模块06c (): LOO 模型比较 ---\n")

suppressPackageStartupMessages({
  library(brms)
  library(loo)
})

# ============================================================================
# 6Ca. 收集所有可用模型
# ============================================================================
if (!exists("brms_fits") || length(brms_fits) == 0) {
  stop("brms_fits 容器为空。请先 source 03b + 04b 加载 6 个主模型。")
}

# 仅对 PASS/WARN 模型跑 LOO（FAIL 模型 loo 不可靠）
loo_target_models <- if (exists("convergence")) {
  names(brms_fits)[
    vapply(names(brms_fits), function(m) {
      row <- convergence[[m]]
      !is.null(row) && row$status %in% c("PASS", "WARN")
    }, FUN.VALUE = logical(1))
  ]
} else {
  names(brms_fits)
}
cat(sprintf("  对 %d 个模型计算 LOO: %s\n",
            length(loo_target_models),
            paste(loo_target_models, collapse = ", ")))

# ============================================================================
# 6Cb. compute_loo_safe() — Pareto k 自动 moment_match 处理
# ============================================================================
compute_loo_safe <- function(fit, model_name, k_threshold = 0.7) {
  cat(sprintf("\n  [%s] 计算 loo()...\n", model_name))
  t0 <- Sys.time()
  loo_obj <- tryCatch(
    brms::loo(fit, cores = 4L),
    error = function(e) {
      cat(sprintf("    ✗ loo() 失败: %s\n", e$message))
      NULL
    }
  )
  if (is.null(loo_obj)) return(NULL)

  pk <- loo_obj$diagnostics$pareto_k
  n_high_k <- sum(pk > k_threshold, na.rm = TRUE)
  n_very_high <- sum(pk > 1.0, na.rm = TRUE)
  cat(sprintf("    elpd_loo = %.2f (SE = %.2f), Pareto k > %.1f: %d / %d (%.2f%%)\n",
              loo_obj$estimates["elpd_loo", "Estimate"],
              loo_obj$estimates["elpd_loo", "SE"],
              k_threshold, n_high_k, length(pk),
              100 * n_high_k / length(pk)))

  # 如有较多高 k，启用 moment_match
  if (n_high_k > 0 && n_high_k <= 50) {
    cat(sprintf("    ⟳ Pareto k 高，启用 moment_match = TRUE 重算...\n"))
    loo_obj <- tryCatch(
      brms::loo(fit, moment_match = TRUE, cores = 4L),
      error = function(e) {
        cat(sprintf("    ⚠ moment_match 失败: %s（保留原 loo）\n", e$message))
        loo_obj
      }
    )
    pk <- loo_obj$diagnostics$pareto_k
    n_high_k_post <- sum(pk > k_threshold, na.rm = TRUE)
    cat(sprintf("    moment_match 后 Pareto k > %.1f: %d / %d\n",
                k_threshold, n_high_k_post, length(pk)))
  }
  elapsed <- difftime(Sys.time(), t0, units = "secs")
  cat(sprintf("    ✓ 完成（%.1f 秒）\n", as.numeric(elapsed)))

  list(loo = loo_obj, n_high_k = n_high_k, n_very_high = n_very_high,
       max_k = max(pk, na.rm = TRUE), n_obs = length(pk))
}

loo_results <- list()
for (mname in loo_target_models) {
  res <- compute_loo_safe(brms_fits[[mname]], mname)
  if (!is.null(res)) loo_results[[mname]] <- res
}

# ============================================================================
# 6Cc. 汇总 loo_results 为 CSV
# ============================================================================
if (length(loo_results) > 0) {
  loo_summary <- do.call(rbind, lapply(names(loo_results), function(m) {
    obj <- loo_results[[m]]$loo
    est <- obj$estimates
    data.frame(
      model = m,
      elpd_loo = est["elpd_loo", "Estimate"],
      elpd_loo_se = est["elpd_loo", "SE"],
      p_loo = est["p_loo", "Estimate"],
      p_loo_se = est["p_loo", "SE"],
      looic = est["looic", "Estimate"],
      looic_se = est["looic", "SE"],
      n_high_k = loo_results[[m]]$n_high_k,
      n_very_high_k = loo_results[[m]]$n_very_high,
      max_pareto_k = loo_results[[m]]$max_k,
      n_obs = loo_results[[m]]$n_obs,
      stringsAsFactors = FALSE
    )
  }))
  rownames(loo_summary) <- NULL
  write.csv(loo_summary,
            file.path(output_dir, "loo_results.csv"),
            row.names = FALSE)
  cat("\n  ✓ loo_results.csv 写入\n")
  print(loo_summary, row.names = FALSE, digits = 3)

  # ============================================================================
  # 6Cd. loo_compare() 关键对比
  # ============================================================================
  # 注: loo_compare 要求模型在相同数据上拟合，故主模型间比较应分组
  # 组 1: 共用 sample_lagged 的模型 (B1_MELS, R3_MELS)
  # 组 2: 共用 sample_B 的模型 (B1_RS)
  # 组 3: 共用 sample_A4 的模型 (B1_RS_PHQ, B1_RS_PC, B1_RS_sigmaM)
  cat("\n  关键 LOO 比较（按数据集分组）:\n")

  compare_groups <- list(
    "sample_lagged" = c("B1_MELS", "R3_MELS"),
    "sample_A4_crosslevel" = c("B1_RS_PHQ", "B1_RS_PC", "B1_RS_sigmaM")
  )

  all_compare_rows <- list()
  for (grp in names(compare_groups)) {
    models_in_grp <- intersect(compare_groups[[grp]], names(loo_results))
    if (length(models_in_grp) < 2) {
      cat(sprintf("    [%s] 跳过（仅 %d 个模型可用）\n",
                  grp, length(models_in_grp)))
      next
    }
    cat(sprintf("\n    === %s: %s ===\n",
                grp, paste(models_in_grp, collapse = " vs ")))
    loo_list <- lapply(models_in_grp, function(m) loo_results[[m]]$loo)
    names(loo_list) <- models_in_grp
    cmp <- tryCatch(loo::loo_compare(loo_list), error = function(e) {
      cat(sprintf("    ⚠ loo_compare 失败: %s\n", e$message))
      NULL
    })
    if (!is.null(cmp)) {
      print(cmp, digits = 3)
      cmp_df <- as.data.frame(cmp)
      cmp_df$model <- rownames(cmp)
      cmp_df$group <- grp
      all_compare_rows[[grp]] <- cmp_df[, c("group", "model",
                                             setdiff(colnames(cmp_df),
                                                      c("group", "model")))]
    }
  }

  if (length(all_compare_rows) > 0) {
    cmp_combined <- do.call(rbind, all_compare_rows)
    rownames(cmp_combined) <- NULL
    write.csv(cmp_combined,
              file.path(output_dir, "loo_comparison.csv"),
              row.names = FALSE)
    cat("\n  ✓ loo_comparison.csv 写入\n")
  }

  # ============================================================================
  # 6Ce. Pareto k 高警告记录
  # ============================================================================
  warn_path <- file.path(logs_dir, "loo_warnings.txt")
  con <- file(warn_path, "wt")
  on.exit(close(con), add = TRUE)
  cat(" LOO Pareto k > 0.7 警告记录\n", file = con)
  cat(sprintf("生成时间: %s\n", format(Sys.time(), "%Y-%m-%d %H:%M:%S")),
      file = con, append = TRUE)
  cat(strrep("-", 70), "\n", sep = "", file = con, append = TRUE)
  for (m in names(loo_results)) {
    pk <- loo_results[[m]]$loo$diagnostics$pareto_k
    high_idx <- which(pk > 0.7)
    cat(sprintf("\n[%s] n_obs = %d, n_high_k = %d (max k = %.3f)\n",
                m, length(pk), length(high_idx),
                max(pk, na.rm = TRUE)),
        file = con, append = TRUE)
    if (length(high_idx) > 0 && length(high_idx) <= 20) {
      cat("  高 k 观测索引: ", paste(high_idx, collapse = ", "), "\n",
          file = con, append = TRUE)
    } else if (length(high_idx) > 20) {
      cat(sprintf("  高 k 观测超过 20 个，前 20: %s\n",
                  paste(head(high_idx, 20), collapse = ", ")),
          file = con, append = TRUE)
    }
  }
  cat(sprintf("  ✓ Pareto k 警告记录: %s\n", warn_path))
}

cat("\n  ✓ 模块06c LOO 模型比较完成 ()\n")
