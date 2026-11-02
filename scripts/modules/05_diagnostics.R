# ============================================================================
# 模块 05b (): 收敛诊断 + 后验预测检查 + 先验预测检查 (gating)
# ============================================================================
# 依赖: 03b/04b（6 个 brms_fits 已就绪）
# 产出: convergence_diagnostics_.csv
#       每个拟合的 trace / pp_check / energy PNG（≈ 18 个图）
#       MELS 专属：sd_grouped pp_check（验证 sigma 子方程）
#        新增: prior_pp_check/ 子目录（6 张先验预测检查图）
#                    *_ppcheck_sd_by_health.png（B1_MELS / R3_MELS 分组 PPC）
# ============================================================================
# 设计决策（计划文件第 2.05 节 +  §2.6 项 D/J）:
#   - PASS:  R-hat<1.01, ESS>400, divergent<1%
#   - WARN:  R-hat∈[1.01,1.05] 或 ESS∈[100,400] 或 divergent∈[1%,5%]
#   - FAIL:  超出 WARN → 自动建议 iter=8000, adapt_delta=0.99 重跑
#   - 不静默接受 WARN（规则 5）
#   - : 在 5Ba 之前增加 prior PPC gating（PASS ≥ 0.70 / FAIL < 0.10 才阻断）
# ============================================================================

cat("\n--- 模块05b (): 收敛诊断 + 后验预测检查 + 先验预测检查 ---\n")

# ============================================================================
# 5B0. 先验预测检查 (gating， 新增，建议文档项 D)
# ============================================================================
# 在批量诊断 (5Bb) 之前验证 6 个模型的先验生成 y 大致落在 mood 量表 [-10, 10] 内。
#  三级判定（详见 02b prior_predictive_check 注释）：
#   - PASS  (≥ 0.70): 先验合理
#   - WARN  (0.10 ≤ pct_in_range < 0.70): 偏宽或 predictor 尺度大；不阻断
#       MELS log-link sigma 与 R3_MELS study_day_squared 系实证 WARN 来源
#   - FAIL  (< 0.10): catastrophic — 先验明显病态，停止流程
cat("\n  [Prior PPC gating] 6 个模型先验预测检查（PASS ≥ 0.70 / FAIL < 0.10）...\n")
prior_ppc_dir <- file.path(figures_dir, "prior_pp_check")
dir.create(prior_ppc_dir, showWarnings = FALSE, recursive = TRUE)

ppc_configs <- list(
  B1_MELS      = list(form = B1_MELS_form,      data = sample_lagged, prior = priors_mels),
  R3_MELS      = list(form = R3_MELS_form,      data = sample_lagged, prior = priors_mels),
  B1_RS        = list(form = B1_RS_form,        data = sample_B,      prior = priors_rs),
  B1_RS_PHQ    = list(form = B1_RS_PHQ_form,    data = sample_A4,     prior = priors_rs),
  B1_RS_PC     = list(form = B1_RS_PC_form,     data = sample_A4,     prior = priors_rs),
  B1_RS_sigmaM = list(form = B1_RS_sigmaM_form, data = sample_A4,     prior = priors_rs)
)

prior_ppc_rows <- list()
for (mname in names(ppc_configs)) {
  cfg <- ppc_configs[[mname]]
  cat(sprintf("    [%s] prior PPC ... ", mname))
  ppc <- prior_predictive_check(cfg$form, cfg$data, cfg$prior, n_draws = 500L)
  prior_ppc_rows[[mname]] <- data.frame(
    model = mname,
    pct_in_range = ppc$pct_in_range,
    median_mean = ppc$median_mean,
    median_sd = ppc$median_sd,
    status = ppc$status,
    pass = ppc$pass
  )
  if (ppc$status == "FAIL") {
    cat(sprintf("FAIL (pct_in_range = %.1f%% < 10%%, catastrophic)\n",
                100 * ppc$pct_in_range))
    stop(sprintf("Prior PPC gating FAIL: %s — 先验明显病态，停止流程", mname))
  }
  cat(sprintf("%s (pct_in_range = %.1f%%, median_mean = %.2f, median_sd = %.2f)\n",
              ppc$status, 100 * ppc$pct_in_range,
              ppc$median_mean, ppc$median_sd))
  tryCatch({
    p_prior <- brms::pp_check(ppc$fit_prior, type = "dens_overlay", ndraws = 50) +
      ggplot2::ggtitle(sprintf("Prior PPC: %s [%s]", mname, ppc$status))
    ggplot2::ggsave(file.path(prior_ppc_dir, paste0(mname, "_prior_ppc.png")),
                    plot = p_prior, width = 8, height = 5, dpi = 100)
  }, error = function(e) {
    cat(sprintf("    ⚠ prior PPC plot 失败: %s\n", e$message))
  })
}
prior_ppc_combined <- do.call(rbind, prior_ppc_rows)
rownames(prior_ppc_combined) <- NULL
write.csv(prior_ppc_combined,
          file.path(output_dir, "prior_ppc_summary_.csv"),
          row.names = FALSE)
n_pass_ppc <- sum(prior_ppc_combined$status == "PASS")
n_warn_ppc <- sum(prior_ppc_combined$status == "WARN")
cat(sprintf("  ✓ Prior PPC 完成：%d PASS, %d WARN, 0 FAIL（catastrophic stop 未触发）\n",
            n_pass_ppc, n_warn_ppc))

# ============================================================================
# 5Ba. 单模型诊断函数
# ============================================================================
diagnose_brms <- function(fit, model_name, save_plots = TRUE) {
  cat(sprintf("\n  [%s] 诊断中...\n", model_name))

  # 后验摘要
  s <- extract_posterior_summary(fit)

  # R-hat / ESS 总体指标
  max_rhat <- max(s$rhat, na.rm = TRUE)
  min_ess_bulk <- min(s$ess_bulk, na.rm = TRUE)
  min_ess_tail <- min(s$ess_tail, na.rm = TRUE)

  # 散度过渡（divergent transitions）
  np <- bayesplot::nuts_params(fit)
  divergent_count <- sum(np$Value[np$Parameter == "divergent__"], na.rm = TRUE)
  total_post_draws <- brms::ndraws(fit)
  divergent_pct <- 100 * divergent_count / total_post_draws

  # 状态判定
  status <- if (max_rhat < 1.01 && min_ess_bulk > 400 && divergent_pct < 1) {
    "PASS"
  } else if (max_rhat < 1.05 && min_ess_bulk >= 100 && divergent_pct < 5) {
    "WARN"
  } else {
    "FAIL"
  }

  # 找出最差的参数（若有 R-hat 超标）
  worst_rhat_param <- s$variable[which.max(s$rhat)]
  worst_ess_param  <- s$variable[which.min(s$ess_bulk)]

  cat(sprintf("    max R-hat = %.4f (%s) | min ESS = %.0f (%s) | divergent = %d (%.2f%%) | %s\n",
              max_rhat, worst_rhat_param,
              min_ess_bulk, worst_ess_param,
              divergent_count, divergent_pct, status))

  if (status != "PASS") {
    cat(sprintf("    ⚠ 建议: 重跑 iter=8000, warmup=4000, adapt_delta=0.99\n"))
  }

  # ============================================================================
  # 图表（仅 PASS/WARN 模型生成；FAIL 模型也生成以便诊断）
  # ============================================================================
  if (save_plots) {
    # Trace plot — top 6 关键参数
    key_params <- c(grep("^b_", brms::variables(fit), value = TRUE)[1:6])
    key_params <- key_params[!is.na(key_params)]
    if (length(key_params) > 0) {
      p_trace <- bayesplot::mcmc_trace(fit, pars = key_params)
      ggplot2::ggsave(file.path(figures_dir, paste0(model_name, "_trace.png")),
                      plot = p_trace, width = 10, height = 6, dpi = 100)
    }

    # Energy diagnostic
    tryCatch({
      p_energy <- bayesplot::mcmc_nuts_energy(np)
      ggplot2::ggsave(file.path(figures_dir, paste0(model_name, "_energy.png")),
                      plot = p_energy, width = 8, height = 4, dpi = 100)
    }, error = function(e) cat(sprintf("    ⚠ energy plot 失败: %s\n", e$message)))

    # pp_check 密度叠加
    tryCatch({
      p_dens <- brms::pp_check(fit, type = "dens_overlay", ndraws = 100)
      ggplot2::ggsave(file.path(figures_dir, paste0(model_name, "_ppcheck_dens.png")),
                      plot = p_dens, width = 8, height = 5, dpi = 100)
    }, error = function(e) cat(sprintf("    ⚠ pp_check dens 失败: %s\n", e$message)))

    # MELS 专属：grouped sd 检查（模型 implied 个人内 SD vs 实测）
    is_mels <- grepl("MELS", model_name)
    if (is_mels) {
      tryCatch({
        p_sd <- brms::pp_check(fit, type = "stat_grouped",
                                group = "user_id", stat = "sd", ndraws = 50)
        ggplot2::ggsave(file.path(figures_dir,
                                   paste0(model_name, "_ppcheck_sd_grouped.png")),
                        plot = p_sd, width = 10, height = 6, dpi = 100)
      }, error = function(e) cat(sprintf("    ⚠ pp_check sd_grouped 失败: %s\n",
                                          e$message)))
    }
  }

  # 返回单行诊断结果
  data.frame(
    model = model_name,
    max_rhat = max_rhat,
    min_ess_bulk = min_ess_bulk,
    min_ess_tail = min_ess_tail,
    divergent_count = divergent_count,
    divergent_pct = divergent_pct,
    total_post_draws = total_post_draws,
    worst_rhat_param = worst_rhat_param,
    worst_ess_param = worst_ess_param,
    status = status
  )
}

# ============================================================================
# 5Bb. 对所有 6 个 brms 拟合执行诊断（ 起 B2/B3 已删，从 10 缩至 6）
# ============================================================================
diag_rows <- list()
for (mname in names(brms_fits)) {
  diag_rows[[mname]] <- diagnose_brms(brms_fits[[mname]], mname, save_plots = TRUE)
  convergence[[mname]] <- diag_rows[[mname]]
}

convergence_combined <- do.call(rbind, diag_rows)
write.csv(convergence_combined,
          file.path(output_dir, "convergence_diagnostics_.csv"),
          row.names = FALSE)

cat("\n  收敛诊断汇总:\n")
print(convergence_combined[, c("model", "max_rhat", "min_ess_bulk",
                                "divergent_pct", "status")],
      row.names = FALSE, digits = 3)

# ============================================================================
# 5Bc. 收敛门：FAIL 即终止主管线，强制人工干预
# ============================================================================
n_fail <- sum(convergence_combined$status == "FAIL")
n_warn <- sum(convergence_combined$status == "WARN")
n_pass <- sum(convergence_combined$status == "PASS")

cat(sprintf("\n  收敛状态: PASS=%d, WARN=%d, FAIL=%d\n", n_pass, n_warn, n_fail))

if (n_fail > 0) {
  cat("\n  ✗ FAIL 模型清单:\n")
  for (m in convergence_combined$model[convergence_combined$status == "FAIL"]) {
    cat(sprintf("    - %s\n", m))
  }
  cat("\n  请删除对应 .rds 后用更高 iter 重跑：\n")
  cat("    FORCE_REFIT='<model_name>' Rscript ema_modeling__master.R\n")
  # v2 修复: FAIL 模型链未混合，下游 R²/表格/图全部无意义。
  # 改为 stop()，阻止 06b/07b/08b 在坏链上算出"看起来很合理"的假数字。
  # 若确需在 FAIL 链上继续（例如单纯做诊断画图），手动 source 06b/07b 即可。
  stop(sprintf(
    " 有 %d 个模型未通过收敛门 — 已阻断下游模块。请删除对应 .rds 后用 FORCE_REFIT 重跑。",
    n_fail))
}

if (n_warn > 0) {
  cat("\n  ⚠ WARN 模型清单（规则 5：不静默接受）:\n")
  for (m in convergence_combined$model[convergence_combined$status == "WARN"]) {
    diag_row <- convergence_combined[convergence_combined$model == m, ]
    cat(sprintf("    - %s: max R-hat=%.4f (%s), min ESS=%.0f (%s), divergent=%.2f%%\n",
                m, diag_row$max_rhat, diag_row$worst_rhat_param,
                diag_row$min_ess_bulk, diag_row$worst_ess_param,
                diag_row$divergent_pct))
  }
}

# ============================================================================
# 5Bd. 分组 PPC ( 新增，建议文档项 J)
# ============================================================================
# 对 B1_MELS / R3_MELS 各出一张 stat_grouped (group = high/low health_session_c)
# 验证 "高 health session 的 mood SD 是否真的更低"——这是 MELS 实质性发现的
# 模型 implied 后验预测视觉化。
cat("\n  [Grouped PPC] B1_MELS / R3_MELS 按 health_level 分组 SD 检查...\n")
# 构造 health_level 分组向量（health_level 不在原拟合 data 中，
# brms::pp_check 即使传 newdata 也只在 fit$data 查 group 列名，
# 故改用 bayesplot::ppc_stat_grouped 并直接传 group 向量）
health_level_vec <- ifelse(sample_lagged$health_session_c > 0, "high", "low")

for (mname in c("B1_MELS", "R3_MELS")) {
  if (is.null(brms_fits[[mname]])) {
    cat(sprintf("    [%s] 跳过（fit 不存在）\n", mname))
    next
  }
  tryCatch({
    y_rep <- brms::posterior_predict(brms_fits[[mname]], ndraws = 100)
    p_grouped <- bayesplot::ppc_stat_grouped(
      y = sample_lagged$mood,
      yrep = y_rep,
      group = health_level_vec,
      stat = "sd"
    ) + ggplot2::ggtitle(sprintf("%s: Posterior SD by health_level", mname))
    out_path <- file.path(figures_dir,
                          paste0(mname, "_ppcheck_sd_by_health.png"))
    ggplot2::ggsave(out_path, plot = p_grouped,
                    width = 8, height = 5, dpi = 120)
    cat(sprintf("    [%s] ✓ %s\n", mname, basename(out_path)))
  }, error = function(e) {
    cat(sprintf("    [%s] ⚠ grouped PPC 失败: %s\n", mname, e$message))
  })
}

cat("\n  ✓ 模块05b 收敛诊断 + 先验/分组 PPC 完成 ()\n")
