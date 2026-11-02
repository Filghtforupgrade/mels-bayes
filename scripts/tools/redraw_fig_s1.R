#!/usr/bin/env Rscript
# ============================================================================
# redraw_fig_s1.R — SM Figure S1 重画
# ============================================================================
# 数据: results/fits/sensitivity_refits/plan_a_illustrative.rds
#       + results/fits/B1_MELS.rds（Plan B 对照面板）
# 设计: 四面板——① 真正失败参数 sd_user_id__sigma_Intercept（核心）② sigma 固定效应
#       （旧图所画参数，R̂ 仅 1.20，对照"为何旧图看不出失败"）③ mu 固定效应
#       ④ Plan B 同参数（清晰混合）。每面板标注 R̂；副标写明批次与采样规格。
# 原脚本 scripts/sensitivity/run_plan_ab_trace.R 不动；本文件为副本改造。
# ============================================================================
options(width = 200)

project_root <- "."  # released: run from repo root (see README)
suppressPackageStartupMessages({
  library(brms); library(ggplot2); library(bayesplot); library(patchwork)
})

fit_a <- readRDS(file.path(project_root, "results/fits/sensitivity_refits/plan_a_illustrative.rds"))
fit_b <- readRDS(file.path(project_root, "results/fits/B1_MELS.rds"))

summarise_par <- function(fit, v) {
  s <- posterior::summarise_draws(posterior::as_draws_df(fit),
                                  median, posterior::default_convergence_measures())
  r <- s[s$variable == v, ]
  if (nrow(r) == 0) return(c(rhat = NA, ess = NA)) else c(rhat = r$rhat, ess = round(r$ess_bulk))
}

make_trace <- function(fit, param, title_str, sub_str) {
  bayesplot::mcmc_trace(fit, pars = param) +
    ggplot2::ggtitle(title_str, subtitle = sub_str) +
    ggplot2::theme_bw(base_size = 9) +
    ggplot2::theme(legend.position = "none",
                   plot.title = element_text(size = 9.5, face = "bold"),
                   plot.subtitle = element_text(size = 7.5),
                   axis.text = element_text(size = 7))
}

batch_a <- "illustrative refit"
batch_b <- "B1_MELS (adopted)"

panels <- list()
panels[[1]] <- make_trace(fit_a, "sd_user_id__sigma_Intercept",
  "1) Plan A: sd_user_id__sigma_Intercept",
  "The failing parameter — R-hat = 1.558, bulk ESS = 7")
panels[[2]] <- make_trace(fit_a, "b_sigma_health_session_c",
  "2) Plan A: sigma health (fixed effect)",
  "Old figure's parameter — R-hat = 1.201")
panels[[3]] <- make_trace(fit_a, "b_health_session_c",
  "3) Plan A: mu health (fixed effect)",
  "R-hat = 1.532")
panels[[4]] <- make_trace(fit_b, "b_sigma_health_session_c",
  "4) Plan B: sigma health (adopted)",
  "Clean chain overlap — R-hat = 1.000")

final <- wrap_plots(panels, ncol = 2) +
  plot_annotation(
    title = "Plan A (sigma with person-level random intercept, FAILED) vs Plan B (fixed-effects-only sigma, ADOPTED)",
    theme = theme(plot.title = element_text(size = 12, face = "bold"))
  )

out <- file.path(project_root, "results/figures/fig_plan_ab_trace_comparison.png")
ggsave(out, plot = final, width = 12, height = 8, dpi = 150)
cat("saved:", out, "\n")
for (v in c("sd_user_id__sigma_Intercept", "b_sigma_health_session_c", "b_health_session_c")) {
  rr <- summarise_par(fit_a, v)
  cat(sprintf("PlanA %-32s Rhat=%.3f ESS=%d\n", v, rr[["rhat"]], rr[["ess"]]))
}
