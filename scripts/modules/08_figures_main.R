# ============================================================================
# 模块 08b (): 6 张主图（论文级别）
# ============================================================================
# 依赖: 03b/04b/05b/06b/07b
# 产出: 6 张主图，存入 figures_dir
#   1) fig_trace_grid.png            — 关键参数 trace
#   2) fig_ppcheck_mels.png          — MELS 4 模型 dens overlay
#   3) fig_sigma_effects.png         — MELS sigma forest（核心发现）
#   4) fig_person_slopes_caterpillar.png — B1_RS 个体斜率
#   5) fig_crosslevel_interactions.png   — 跨层调节简单斜率
#   6) fig_brms_vs_lme4_forest.png   — 跨范式对照 forest
# ============================================================================

cat("\n--- 模块08b (): 6 张主图 ---\n")

suppressPackageStartupMessages({
  library(ggplot2)
  library(bayesplot)
})

# ============================================================================
# 8Ba. Fig 1: trace_grid — 6 个关键参数 trace
# ============================================================================
fig_trace_grid <- function() {
  out_path <- file.path(figures_dir, "fig_trace_grid.png")
  # 选择 B1_MELS 的关键 mu/sigma 参数
  fit <- brms_fits$B1_MELS
  if (is.null(fit)) return(invisible(NULL))
  key_params <- c("b_health_session_c", "b_health_session_between_c",
                  "b_mood_prev_session_c",
                  "b_sigma_health_session_c",
                  "sd_user_id__Intercept",
                  "sd_user_id__sigma_Intercept")
  available <- intersect(key_params, brms::variables(fit))
  if (length(available) == 0) return(invisible(NULL))
  p <- bayesplot::mcmc_trace(fit, pars = available) +
    ggplot2::ggtitle("B1_MELS — Trace plots of key parameters") +
    #  F1: mcmc_trace 在本环境默认深色底，期刊要求浅色 → 显式 theme_bw
    ggplot2::theme_bw() +
    ggplot2::theme(text = ggplot2::element_text(size = 14),
                   axis.text = ggplot2::element_text(size = 12),
                   plot.title = ggplot2::element_text(size = 16))
  ggplot2::ggsave(out_path, plot = p, width = 12, height = 8, dpi = 120)
  cat(sprintf("  ✓ Fig 1: %s\n", basename(out_path)))
}
fig_trace_grid()

# ============================================================================
# 8Bb. Fig 2: ppcheck_mels — 主 MELS + 条目级 MELS 的 dens overlay
# ============================================================================
#  A11: 扩展为 6-panel，加入 4 个 item-level MELS 的 PPC
fig_ppcheck_mels <- function() {
  out_path <- file.path(figures_dir, "fig_ppcheck_mels.png")
  mels_models <- c("B1_MELS", "R3_MELS")
  available <- mels_models[mels_models %in% names(brms_fits)]
  if (length(available) == 0) return(invisible(NULL))

  plots <- lapply(available, function(m) {
    brms::pp_check(brms_fits[[m]], type = "dens_overlay", ndraws = 100) +
      ggplot2::ggtitle(m)
  })
  
  #  A11: 尝试加载条目级 MELS 的 PPC
  item_root <- file.path(project_root, "legacy", "bayesian_analysis_",
                        "modeling_results_")
  item_dirs <- list.dirs(item_root, recursive = FALSE)
  cand_dirs <- item_dirs[grepl("item_mels_[34]_\\d{8}_\\d{6}$", basename(item_dirs))]
  if (length(cand_dirs) > 0) {
    latest_dir <- cand_dirs[which.max(file.mtime(cand_dirs))]
    item_models <- c("Item_MELS_Noon_Main", "Item_MELS_Evening_Main",
                     "Item_MELS_Noon_Robust", "Item_MELS_Evening_Robust")
    for (im in item_models) {
      rds_path <- file.path(latest_dir, "fits", paste0(im, ".rds"))
      if (file.exists(rds_path)) {
        fit_i <- readRDS(rds_path)
        p_i <- tryCatch(
          brms::pp_check(fit_i, type = "dens_overlay", ndraws = 60) +
            ggplot2::ggtitle(gsub("Item_MELS_", "", im)) +
            ggplot2::theme_minimal(base_size = 9) +
            ggplot2::theme(plot.title = ggplot2::element_text(size = 9)),
          error = function(e) NULL)
        if (!is.null(p_i)) plots[[length(plots) + 1]] <- p_i
      }
    }
    cat(sprintf("  ✓ added %d item-level PPC panels\n",
                length(plots) - length(available)))
  }
  
  # 用 patchwork 拼图（如缺失则降级到逐图保存）
  if (requireNamespace("patchwork", quietly = TRUE)) {
    library(patchwork)
    ncol <- if (length(plots) >= 6) 3 else min(length(plots), 3)
    p <- Reduce(`+`, plots) + plot_layout(ncol = ncol)
    ggplot2::ggsave(out_path, plot = p, width = 14, height = 10, dpi = 120)
    cat(sprintf("  ✓ Fig 2: %s (%d panels)\n", basename(out_path), length(plots)))
  } else {
    cat("  ⚠ patchwork 未安装，逐图保存 ppcheck\n")
    for (i in seq_along(plots)) {
      ggplot2::ggsave(file.path(figures_dir,
                                 sprintf("fig_ppcheck_%s.png", available[i])),
                      plot = plots[[i]], width = 8, height = 5, dpi = 120)
    }
  }
}
fig_ppcheck_mels()

# ============================================================================
# 8Bc. Fig 3: sigma_effects — MELS sigma 子方程系数 forest（核心发现）
# ============================================================================
fig_sigma_effects <- function() {
  out_path <- file.path(figures_dir, "fig_sigma_effects.png")
  if (!exists("mels_sigma_results") || length(mels_sigma_results) == 0) {
    return(invisible(NULL))
  }
  df <- do.call(rbind, mels_sigma_results)
  if (nrow(df) == 0) return(invisible(NULL))

  df$model_param <- paste0(df$model, " | ", df$parameter)
  df$model_param <- factor(df$model_param, levels = rev(unique(df$model_param)))
  
  #  A9: 按 L1/L2/interaction 着色
  df$predictor_type <- dplyr::case_when(
    grepl("_session_c$", df$parameter) & !grepl("_between_c$", df$parameter) ~ "L1 within-person",
    grepl("_between_c$", df$parameter) ~ "L2 between-person",
    grepl("_x_|:", df$parameter) ~ "cross-level interaction",
    TRUE ~ "other"
  )
  col_map <- c("L1 within-person" = "#2E8B57", "L2 between-person" = "#C0392B",
               "cross-level interaction" = "#8E44AD", "other" = "grey50")

  p <- ggplot2::ggplot(df,
                        ggplot2::aes(x = median, y = model_param, color = predictor_type)) +
    ggplot2::geom_vline(xintercept = 0, linetype = "dashed",
                         color = "grey40") +
    ggplot2::geom_errorbarh(ggplot2::aes(xmin = q2.5, xmax = q97.5),
                             height = 0.2) +
    ggplot2::geom_point(size = 3) +
    ggplot2::scale_color_manual(values = col_map, name = "Predictor Level") +
    ggplot2::labs(
      title = "MELS sigma-submodel coefficients (log-SD link)",
      subtitle = "Negative values with 95% CrI excluding 0 indicate compression of residual mood variability",
      x = "Posterior Median (log-SD scale)",
      y = NULL
    ) +
    ggplot2::theme_minimal(base_size = 13)
  ggplot2::ggsave(out_path, plot = p, width = 10, height = 6, dpi = 120)
  cat(sprintf("  ✓ Fig 3: %s\n", basename(out_path)))
}
fig_sigma_effects()

# ============================================================================
# 8Bd. Fig 4: B1_RS 个体斜率 caterpillar
# ============================================================================
fig_person_slopes <- function() {
  out_path <- file.path(figures_dir, "fig_person_slopes_caterpillar.png")
  if (!exists("person_slopes_B1") || nrow(person_slopes_B1) == 0) {
    return(invisible(NULL))
  }
  df <- person_slopes_B1[order(person_slopes_B1$person_slope), ]
  df$rank <- seq_len(nrow(df))
  pop_mean <- brms::fixef(brms_fits$B1_RS)["health_session_c", "Estimate"]

  df$color_cat <- with(df,
    ifelse(credibly_positive, "credibly +",
      ifelse(credibly_negative, "credibly −", "uncertain")))

  p <- ggplot2::ggplot(df, ggplot2::aes(x = rank, y = person_slope,
                                          color = color_cat)) +
    ggplot2::geom_hline(yintercept = 0, linetype = "dashed",
                         color = "grey40") +
    ggplot2::geom_hline(yintercept = pop_mean, linetype = "solid",
                         color = "darkred", linewidth = 0.5) +
    ggplot2::geom_errorbar(ggplot2::aes(ymin = pop_mean + q2.5,
                                          ymax = pop_mean + q97.5),
                            alpha = 0.5, width = 0) +
    ggplot2::geom_point(size = 0.8) +
    ggplot2::scale_color_manual(values = c("credibly +" = "#1b7837",
                                            "uncertain" = "grey60",
                                            "credibly −" = "#762a83")) +
    ggplot2::labs(
      title = "B1_RS: person-specific health-to-mood slopes (sorted by slope)",
      subtitle = sprintf("Red line = population average (%.3f); colored = 95%% CrI entirely on one side of zero", pop_mean),
      x = "Person rank (sorted by slope)",
      y = "Person-specific slope on mood",
      color = "95% CrI category"
    ) +
    ggplot2::theme_minimal(base_size = 12) +
    ggplot2::theme(legend.position = "top")
  ggplot2::ggsave(out_path, plot = p, width = 12, height = 6, dpi = 120)
  cat(sprintf("  ✓ Fig 4: %s\n", basename(out_path)))
}
fig_person_slopes()

# ============================================================================
# 8Be. Fig 5: 跨层调节 simple slopes（B1_RS_PHQ / PC / sigmaM）
# ============================================================================
fig_crosslevel_interactions <- function() {
  out_path <- file.path(figures_dir, "fig_crosslevel_interactions.png")
  models <- list(
    PHQ = list(fit = brms_fits$B1_RS_PHQ,
               mod_var = "pre_phq9_c",
               data = sample_A4),
    PC = list(fit = brms_fits$B1_RS_PC,
              mod_var = "pre_control_c",
              data = sample_A4),
    sigmaM = list(fit = brms_fits$B1_RS_sigmaM,
                  mod_var = "person_std_mood_c",
                  data = sample_A4)
  )
  models <- models[!sapply(models, function(x) is.null(x$fit))]
  if (length(models) == 0) return(invisible(NULL))

  panels <- list()
  for (mname in names(models)) {
    m <- models[[mname]]
    # 简单斜率：在 mod_var 的 -1SD/Mean/+1SD 三个水平上预测 mood ~ health_session_c
    sd_mod <- sd(m$data[[m$mod_var]], na.rm = TRUE)
    mod_levels <- c(-sd_mod, 0, sd_mod)
    health_grid <- seq(quantile(m$data$health_session_c, 0.05, na.rm = TRUE),
                        quantile(m$data$health_session_c, 0.95, na.rm = TRUE),
                        length.out = 30)

    newdata <- expand.grid(
      health_session_c = health_grid,
      mod_lvl = mod_levels
    )
    newdata[[m$mod_var]] <- newdata$mod_lvl
    newdata$health_session_between_c <- 0
    newdata$mood_prev_session_c <- 0
    newdata$study_day_c <- 0
    newdata$session_numeric <- 0

    pred <- brms::posterior_epred(m$fit, newdata = newdata,
                                    re_formula = NA, ndraws = 500)
    newdata$pred_median <- apply(pred, 2, median)
    newdata$pred_q2.5 <- apply(pred, 2, quantile, 0.025)
    newdata$pred_q97.5 <- apply(pred, 2, quantile, 0.975)
    newdata$mod_label <- factor(newdata$mod_lvl,
      labels = c("Low (-1 SD)", "Mean", "High (+1 SD)"))

    p <- ggplot2::ggplot(newdata,
                          ggplot2::aes(x = health_session_c,
                                        y = pred_median,
                                        color = mod_label,
                                        fill = mod_label)) +
      ggplot2::geom_ribbon(ggplot2::aes(ymin = pred_q2.5, ymax = pred_q97.5),
                            alpha = 0.15, color = NA) +
      ggplot2::geom_line(linewidth = 1) +
      ggplot2::labs(
        title = sprintf("B1_RS × %s", mname),
        x = "health_session_c (within-person)",
        y = "Predicted mood",
        color = mname, fill = mname
      ) +
      ggplot2::theme_minimal(base_size = 11) +
      ggplot2::theme(legend.position = "bottom")
    panels[[mname]] <- p
  }

  if (requireNamespace("patchwork", quietly = TRUE)) {
    library(patchwork)
    p_all <- Reduce(`+`, panels) + plot_layout(ncol = length(panels))
    ggplot2::ggsave(out_path, plot = p_all,
                     width = 5 * length(panels), height = 5, dpi = 120)
    cat(sprintf("  ✓ Fig 5: %s\n", basename(out_path)))
  } else {
    for (mname in names(panels)) {
      ggplot2::ggsave(file.path(figures_dir,
                                 sprintf("fig_simple_slopes_%s.png", mname)),
                      plot = panels[[mname]], width = 6, height = 5, dpi = 120)
    }
  }
}
fig_crosslevel_interactions()

# ============================================================================
# 8Bf. Fig 6: brms vs lme4 对照 forest
# ============================================================================
fig_brms_vs_lme4_forest <- function() {
  out_path <- file.path(figures_dir, "fig_brms_vs_lme4_forest.png")
  cmp_csv <- file.path(output_dir, "comparison_brms_vs_lme4.csv")
  if (!file.exists(cmp_csv)) return(invisible(NULL))
  cmb <- read.csv(cmp_csv)
  # 只取 MELS 模型对比，避免图过密（：仅 B1 + R3）
  cmb <- cmb[cmb$brms_model %in% c("B1_MELS", "R3_MELS"), ]
  if (nrow(cmb) == 0) return(invisible(NULL))

  long_df <- rbind(
    data.frame(
      pipeline = "lme4 REML",
      model = cmb$brms_model,
      term = cmb$term,
      estimate = cmb$lme4_beta,
      lower = cmb$lme4_beta - 1.96 * cmb$lme4_se,
      upper = cmb$lme4_beta + 1.96 * cmb$lme4_se
    ),
    data.frame(
      pipeline = "brms posterior",
      model = cmb$brms_model,
      term = cmb$term,
      estimate = cmb$brms_median,
      lower = cmb$brms_q2.5,
      upper = cmb$brms_q97.5
    )
  )
  long_df$model_term <- paste0(long_df$model, " | ", long_df$term)
  long_df$model_term <- factor(long_df$model_term,
                                 levels = rev(unique(long_df$model_term)))

  p <- ggplot2::ggplot(long_df,
                        ggplot2::aes(x = estimate, y = model_term,
                                      color = pipeline)) +
    ggplot2::geom_vline(xintercept = 0, linetype = "dashed",
                         color = "grey40") +
    ggplot2::geom_errorbarh(ggplot2::aes(xmin = lower, xmax = upper),
                             position = ggplot2::position_dodge(width = 0.5),
                             height = 0.2) +
    ggplot2::geom_point(position = ggplot2::position_dodge(width = 0.5),
                         size = 2) +
    ggplot2::scale_color_manual(values = c("lme4 REML" = "darkorange",
                                            "brms posterior" = "steelblue")) +
    ggplot2::labs(
      title = "brms vs lme4 cross-paradigm comparison (B1_MELS and R3_MELS)",
      x = "Estimate (95% CI/CrI)", y = NULL, color = NULL
    ) +
    ggplot2::theme_minimal(base_size = 11) +
    ggplot2::theme(legend.position = "top")
  ggplot2::ggsave(out_path, plot = p, width = 12, height = 10, dpi = 120)
  cat(sprintf("  ✓ Fig 6: %s\n", basename(out_path)))
}
fig_brms_vs_lme4_forest()

cat("\n  ✓ 模块08b 6 张主图完成\n")
