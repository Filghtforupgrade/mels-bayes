# ============================================================================
# 模块 08c (): 4 张新论文图 — 对标 Dora et al. (2024) AMPPS tutorial
# ============================================================================
# 依赖: 03b/04b/05b 拟合完成（特别是 B1_MELS, R3_MELS 用 sample_prior="yes" 重跑）
# 产出: 4 张新 PNG + 1 张合并 PNG，存入 figures_dir 与 06_paper/figures/
#   1) fig_prior_distributions.png     — Dora Fig.2 对标（先验理论曲线）
#   2) fig_mu_posterior_densities.png  — Dora Fig.6/11 对标（mu 后验密度）
#   3) fig_prior_posterior_overlay.png — Dora Fig.8 对标（先验-后验对比）
#   4) fig_item_mels_ppc.png           — 顶刊要求（item-level PPC）
#   5) fig_item_mels_combined.png      — 优化既有 4a/4b 为单图分面
# ============================================================================
#  设计决策:
#   - 全部输出 dpi=150（A4 paper standard）；尺寸 8x6 in（双栏稿件适配）
#   - 同时复制到 06_paper/figures/ 供 LaTeX/Markdown 编译
#   - 标签全英文（投稿前最终化）
# ============================================================================

cat("\n--- 模块08c (): 4 张新论文图（Dora 2024 对标）---\n")

suppressPackageStartupMessages({
  library(ggplot2)
  library(bayesplot)
  library(patchwork)
  library(dplyr)
})

# 论文图输出到 figures_dir（paper/figures 通过符号链接指向此处）
save_dual <- function(p, basename, width, height, dpi = 150) {
  out <- file.path(figures_dir, basename)
  ggplot2::ggsave(out, plot = p, width = width, height = height,
                  dpi = dpi, bg = "white")
  cat(sprintf("  ✓ %s (%dx%d in @ %d dpi)\n", basename, width, height, dpi))
}

# ============================================================================
# 8Ca. Fig S1 — 先验分布理论曲线（Dora 2024 Fig.2 对标）
# ============================================================================
fig_prior_distributions <- function() {
  out_name <- "fig_prior_distributions.png"

  # student_t(df, mu, sigma) 密度（dt 是标准 t，需位置-尺度变换）
  dt_loc_scale <- function(x, df, mu, sigma) {
    dt((x - mu) / sigma, df = df) / sigma
  }
  # LKJ correlation density（对于 2x2 矩阵，rho ~ Beta((eta-1+K/2), (eta-1+K/2))
  # 但更准确地，LKJ(eta) 对 K=2 给出 rho 的密度 ∝ (1-rho^2)^(eta-1)
  dlkj_2x2 <- function(rho, eta) {
    norm_c <- gamma(eta) / (sqrt(pi) * gamma(eta - 0.5))
    norm_c * (1 - rho^2)^(eta - 1)
  }

  prior_grid <- list(
    list(label = "normal(0, 1)\nfixed slopes (b)",
         x_range = c(-4, 4),
         fun = function(x) dnorm(x, 0, 1),
         color = "#1f77b4"),
    list(label = "normal(3, 2)\nIntercept",
         x_range = c(-3, 9),
         fun = function(x) dnorm(x, 3, 2),
         color = "#ff7f0e"),
    list(label = "student_t(3, 0, 1)\nsd, sigma_Intercept",
         x_range = c(-6, 6),
         fun = function(x) dt_loc_scale(x, 3, 0, 1),
         color = "#2ca02c"),
    list(label = "student_t(3, 0, 2.5)\nsigma (non-MELS)",
         x_range = c(-10, 10),
         fun = function(x) dt_loc_scale(x, 3, 0, 2.5),
         color = "#d62728"),
    list(label = "LKJ(2)\nRE correlation",
         x_range = c(-1, 1),
         fun = function(x) dlkj_2x2(x, 2),
         color = "#9467bd")
  )

  panels <- lapply(prior_grid, function(cfg) {
    df <- data.frame(x = seq(cfg$x_range[1], cfg$x_range[2], length.out = 400))
    df$y <- cfg$fun(df$x)
    ggplot2::ggplot(df, ggplot2::aes(x = x, y = y)) +
      ggplot2::geom_area(fill = cfg$color, alpha = 0.4) +
      ggplot2::geom_line(color = cfg$color, linewidth = 1) +
      ggplot2::labs(title = cfg$label, x = "Parameter value", y = "Density") +
      ggplot2::theme_minimal(base_size = 10) +
      ggplot2::theme(plot.title = ggplot2::element_text(size = 9, lineheight = 1.0))
  })

  p_combined <- patchwork::wrap_plots(panels, ncol = 3) +
    patchwork::plot_annotation(
      title = "Prior distributions used in the Bayesian MELS / RS models",
      subtitle = NULL
    )
  save_dual(p_combined, out_name, width = 11, height = 6)
}
tryCatch(fig_prior_distributions(), error = function(e) {
  cat(sprintf("  ✗ Fig S1 失败: %s\n", e$message))
})

# ============================================================================
# 8Cb. Fig S2 — mu 后验密度图（Dora 2024 Fig.6/11 对标）
# ============================================================================
fig_mu_posterior_densities <- function() {
  out_name <- "fig_mu_posterior_densities.png"
  fit <- brms_fits$B1_MELS
  if (is.null(fit)) {
    cat("  ⚠ B1_MELS not available\n")
    return(invisible(NULL))
  }
  # 选取 mu 子方程的核心参数 + sigma 子方程的 health 项
  candidate <- c("b_health_session_c", "b_health_session_between_c",
                 "b_mood_lag1_day_c", "b_study_day_c",
                 "b_sigma_health_session_c", "b_sigma_health_session_between_c")
  pars <- intersect(candidate, brms::variables(fit))
  if (length(pars) == 0) {
    cat("  ⚠ no key params found in B1_MELS\n")
    return(invisible(NULL))
  }
  p <- bayesplot::mcmc_areas(fit, pars = pars, prob = 0.95, prob_outer = 1) +
    ggplot2::ggtitle("B1_MELS: posterior densities (95% CrI shaded)") +
    ggplot2::labs(subtitle = "Black bar = posterior median; shaded = 95% CrI.") +
    ggplot2::theme_minimal(base_size = 11)
  save_dual(p, out_name, width = 9, height = 7)
}
tryCatch(fig_mu_posterior_densities(), error = function(e) {
  cat(sprintf("  ✗ Fig S2 失败: %s\n", e$message))
})

# ============================================================================
# 8Cc. Fig S3 — 先验-后验对比图（Dora 2024 Fig.8 对标）
# ============================================================================
# 要求: B1_MELS 必须用 sample_prior = "yes" 拟合（ 已启用）
fig_prior_posterior_overlay <- function() {
  out_name <- "fig_prior_posterior_overlay.png"
  fit <- brms_fits$B1_MELS
  if (is.null(fit)) {
    cat("  ⚠ B1_MELS not available\n")
    return(invisible(NULL))
  }
  draws <- brms::as_draws_df(fit)
  prior_cols <- grep("^prior_", colnames(draws), value = TRUE)
  if (length(prior_cols) == 0) {
    cat("  ⚠ sample_prior='yes' draws not found. Re-fit with  settings.\n")
    return(invisible(NULL))
  }
  # 关注的参数对：(posterior name, prior name)
  pair_specs <- list(
    list(post = "b_health_session_c",       prior = "prior_b",
         label = "b_health_session_c"),
    list(post = "b_sigma_health_session_c", prior = "prior_b_sigma",
         label = "b_sigma_health_session_c (log-SD scale)"),
    list(post = "b_mood_lag1_day_c",        prior = "prior_b",
         label = "b_mood_lag1_day_c"),
    list(post = "Intercept",                prior = "prior_Intercept",
         label = "Intercept (mu)")
  )

  panels <- list()
  for (sp in pair_specs) {
    if (!(sp$post %in% colnames(draws)) ||
        !(sp$prior %in% colnames(draws))) next
    df_post  <- data.frame(value = draws[[sp$post]],  type = "Posterior")
    df_prior <- data.frame(value = draws[[sp$prior]], type = "Prior")
    df <- rbind(df_post, df_prior)
    df$type <- factor(df$type, levels = c("Prior", "Posterior"))

    p <- ggplot2::ggplot(df, ggplot2::aes(x = value, fill = type, color = type)) +
      ggplot2::geom_density(alpha = 0.35, linewidth = 0.8) +
      ggplot2::scale_fill_manual(values = c("Prior" = "grey60",
                                              "Posterior" = "#1f77b4")) +
      ggplot2::scale_color_manual(values = c("Prior" = "grey30",
                                               "Posterior" = "#1f3a6d")) +
      ggplot2::labs(title = sp$label, x = NULL, y = "Density", fill = NULL, color = NULL) +
      ggplot2::theme_minimal(base_size = 10) +
      ggplot2::theme(legend.position = "bottom")
    panels[[sp$post]] <- p
  }

  if (length(panels) == 0) {
    cat("  ⚠ no prior-posterior pairs available\n")
    return(invisible(NULL))
  }
  # 先在每个 panel 设置 legend，避免 patchwork `&` 兼容问题
  panels <- lapply(panels, function(p)
    p + ggplot2::theme(legend.position = "bottom"))
  p_combined <- patchwork::wrap_plots(panels, ncol = 2) +
    patchwork::plot_layout(guides = "collect") +
    patchwork::plot_annotation(
      title = "B1_MELS: prior vs posterior densities",
      subtitle = "Dora et al. (2024) Fig.8 analog. Data dominates priors for within-person effects; intercept prior is informative."
    )
  save_dual(p_combined, out_name, width = 10, height = 7)
}
tryCatch(fig_prior_posterior_overlay(), error = function(e) {
  cat(sprintf("  ✗ Fig S3 失败: %s\n", e$message))
})

# ============================================================================
# 8Cd. Fig S4 — Item-level MELS 的 PPC（顶刊要求）
# ============================================================================
fig_item_mels_ppc <- function() {
  out_name <- "fig_item_mels_ppc.png"
  # 寻找最新的 item_mels 输出目录
  modeling_root <- file.path(project_root, "legacy", "bayesian_analysis_",
                              "modeling_results_")
  if (!dir.exists(modeling_root)) {
    cat(sprintf("  ⚠ %s 不存在\n", modeling_root))
    return(invisible(NULL))
  }
  item_dirs <- list.dirs(modeling_root, recursive = FALSE)
  # 取最新的 item_mels_<N>_<TS> 目录
  cand_dirs <- item_dirs[grepl("item_mels_[34]_\\d{8}_\\d{6}$",
                                basename(item_dirs))]
  if (length(cand_dirs) == 0) {
        cand_dirs <- item_dirs[grepl("item_mels_\\d_\\d{8}_\\d{6}$",
                                  basename(item_dirs))]
  }
  if (length(cand_dirs) == 0) {
    cat("  ⚠ 未找到 item_mels_* 目录\n")
    return(invisible(NULL))
  }
  latest_dir <- cand_dirs[which.max(file.mtime(cand_dirs))]
  cat(sprintf("    使用 item-level 目录: %s\n", latest_dir))

  item_names <- c("Item_MELS_Noon_Main", "Item_MELS_Evening_Main",
                  "Item_MELS_Noon_Robust", "Item_MELS_Evening_Robust")
  panels <- list()
  for (im in item_names) {
    rds_path <- file.path(latest_dir, "fits", paste0(im, ".rds"))
    if (!file.exists(rds_path)) next
    fit_i <- readRDS(rds_path)
    p <- tryCatch(
      brms::pp_check(fit_i, type = "dens_overlay", ndraws = 80) +
        ggplot2::ggtitle(im) +
        ggplot2::theme_minimal(base_size = 9) +
        ggplot2::theme(plot.title = ggplot2::element_text(size = 9)),
      error = function(e) {
        cat(sprintf("    ⚠ pp_check 失败 (%s): %s\n", im, e$message))
        NULL
      }
    )
    if (!is.null(p)) panels[[im]] <- p
  }
  if (length(panels) == 0) {
    cat("  ⚠ no item-level PPCs generated\n")
    return(invisible(NULL))
  }
  panels <- lapply(panels, function(p)
    p + ggplot2::theme(legend.position = "bottom"))
  p_combined <- patchwork::wrap_plots(panels, ncol = 2) +
    patchwork::plot_layout(guides = "collect") +
    patchwork::plot_annotation(
      title = "Item-level MELS: posterior predictive checks (density overlay)",
      subtitle = "Each panel = y (dark) vs y_rep (80 draws, light). Good fit = light curves envelop dark curve."
    )
  save_dual(p_combined, out_name, width = 11, height = 8)
}
tryCatch(fig_item_mels_ppc(), error = function(e) {
  cat(sprintf("  ✗ Fig S4 失败: %s\n", e$message))
})

# ============================================================================
# 8Ce. 优化 — 合并既有 fig_item_mels_mu + fig_item_mels_sigma 为单图分面
# ============================================================================
# 注: 原图是 ggplot 对象但已 saveRDS 为 PNG，无法直接重组。
# 此处通过重新读取 item_mels_effects.csv 重画为分面 forest。
fig_item_mels_combined <- function() {
  out_name <- "fig_item_mels_combined.png"
  modeling_root <- file.path(project_root, "legacy", "bayesian_analysis_",
                              "modeling_results_")
  item_dirs <- list.dirs(modeling_root, recursive = FALSE)
  cand_dirs <- item_dirs[grepl("item_mels_[134]_\\d{8}_\\d{6}$",
                                basename(item_dirs))]
  if (length(cand_dirs) == 0) {
    cat("  ⚠ no item-level dir for combined figure\n")
    return(invisible(NULL))
  }
  latest_dir <- cand_dirs[which.max(file.mtime(cand_dirs))]
  csv_path <- file.path(latest_dir, "item_mels_effects.csv")
  if (!file.exists(csv_path)) {
    cat(sprintf("  ⚠ %s 不存在\n", csv_path))
    return(invisible(NULL))
  }
  df <- read.csv(csv_path, stringsAsFactors = FALSE)
  if (nrow(df) == 0) return(invisible(NULL))

  # df 期望列：model, item, mu_median, mu_q2.5, mu_q97.5,
  #                       sig_log_median, sig_log_q2.5, sig_log_q97.5
  required_cols <- c("model", "item", "mu_median", "mu_q2.5", "mu_q97.5",
                     "sig_log_median", "sig_log_q2.5", "sig_log_q97.5")
  missing_cols <- setdiff(required_cols, colnames(df))
  if (length(missing_cols) > 0) {
    cat(sprintf("  ⚠ item_mels_effects.csv 缺列: %s\n",
                paste(missing_cols, collapse = ", ")))
    return(invisible(NULL))
  }

  df_long <- rbind(
    data.frame(model = df$model, item = df$item,
               estimate = df$mu_median, low = df$mu_q2.5, high = df$mu_q97.5,
               submodel = "mu (location)"),
    data.frame(model = df$model, item = df$item,
               estimate = df$sig_log_median, low = df$sig_log_q2.5,
               high = df$sig_log_q97.5,
               submodel = "sigma (log-SD scale)")
  )
  df_long$item_model <- paste(df_long$item, df_long$model, sep = " | ")
  df_long$item_model <- factor(df_long$item_model,
                                 levels = unique(df_long$item_model))

  p <- ggplot2::ggplot(df_long,
                        ggplot2::aes(x = estimate, y = item_model)) +
    ggplot2::geom_vline(xintercept = 0, linetype = "dashed", color = "grey50") +
    ggplot2::geom_errorbarh(ggplot2::aes(xmin = low, xmax = high),
                             height = 0.2, color = "steelblue") +
    ggplot2::geom_point(size = 1.6, color = "steelblue") +
    ggplot2::facet_wrap(~ submodel, ncol = 2, scales = "free_x") +
    ggplot2::labs(
      title = "Item-level MELS: posterior medians + 95% CrI",
      subtitle = "mu (left) and sigma (right) facets. Heterogeneity is the central finding (aspiration = pure sigma-compressor).",
      x = "Posterior median", y = NULL
    ) +
    ggplot2::theme_minimal(base_size = 9) +
    ggplot2::theme(strip.text = ggplot2::element_text(face = "bold"))
  save_dual(p, out_name, width = 12, height = 7)
}
tryCatch(fig_item_mels_combined(), error = function(e) {
  cat(sprintf("  ✗ Fig combined 失败: %s\n", e$message))
})

cat("\n  ✓ 模块08c  新论文图完成\n")
