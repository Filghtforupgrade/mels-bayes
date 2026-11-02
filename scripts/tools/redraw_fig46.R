#!/usr/bin/env Rscript
# ============================================================================
# redraw_fig46.R — Figure 4 分面 + Figure 6 面板标签
# ============================================================================
# Figure 4: sigma 子方程系数森林图改分面（predictor level 三个小面 + other），
#           free x 轴——修复 Intercept(≈0.92) 把轴拉到 0–1 导致小效应 CrI 不可辨的问题。
#           数据源 results/posterior_summaries/mels_sigma_coefficients.csv（管线权威输出）。
# Figure 6: conditional effects 三面板，面板标题补全 moderator 名称、x 轴补构念定义。
# 原脚本（08_figures_main.R / run_conditional_effects.R）不动；本文件为副本改造。
# ============================================================================
options(width = 200)

project_root <- "."  # released: run from repo root (see README)
SCRIPT_DIR   <- file.path(project_root, "scripts", "pipeline")
MODULE_DIR   <- file.path(project_root, "scripts", "modules")
source(file.path(MODULE_DIR, "00_setup.R"))
source(file.path(MODULE_DIR, "01_data.R"))
source(file.path(MODULE_DIR, "02_priors.R"))

suppressPackageStartupMessages({
  library(ggplot2); library(patchwork)
})

# ---------------- Figure 4: 分面 sigma forest ----------------
df <- read.csv(file.path(project_root, "results/posterior_summaries/mels_sigma_coefficients.csv"))
# 读者友好标签 + 分面内仅显示实际参数（free y 去空白行）
model_lab <- c("B1_MELS" = "B1 (primary MELS)", "R3_MELS" = "R3 (extension MELS)")
param_lab <- c(
  "health_session_c"         = "Within-person health behavior",
  "health_session_between_c" = "Between-person health behavior (person mean)",
  "health_x_prev_mood_c"     = "Health × prior-mood interaction",
  "Intercept"                = "Intercept"
)
df$model_param <- paste0(model_lab[df$model], ": ", param_lab[df$parameter])
df$predictor_type <- dplyr::case_when(
  grepl("_session_c$", df$parameter) & !grepl("_between_c$", df$parameter) ~ "L1 within-person",
  grepl("_between_c$", df$parameter) ~ "L2 between-person",
  grepl("_x_|:", df$parameter) ~ "cross-level interaction",
  TRUE ~ "other (intercepts)"
)
lv <- c("L1 within-person", "L2 between-person", "cross-level interaction", "other (intercepts)")
df$predictor_type <- factor(df$predictor_type, levels = lv)
df <- df[order(df$predictor_type, -abs(df$median)), ]
df$model_param <- factor(df$model_param, levels = rev(unique(df$model_param)))

col_map <- c("L1 within-person" = "#2E8B57", "L2 between-person" = "#C0392B",
             "cross-level interaction" = "#8E44AD", "other (intercepts)" = "grey50")

p4 <- ggplot(df, aes(x = median, y = model_param)) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "grey40") +
  geom_errorbarh(aes(xmin = q2.5, xmax = q97.5, color = predictor_type), height = 0.2) +
  geom_point(aes(color = predictor_type), size = 2.6) +
  facet_wrap(~predictor_type, scales = "free", ncol = 2) +
  scale_color_manual(values = col_map, guide = "none") +
  labs(
    title = "MELS sigma-submodel coefficients (log-SD link), faceted by predictor level",
    subtitle = "Each panel shows only the parameters estimated at that level,\nwith panel-specific axes so that small effect sizes remain readable; dashed line = 0",
    x = "Posterior median (log-SD scale; panel-specific axis)",
    y = NULL) +
  theme_bw(base_size = 10) +
  theme(strip.text = element_text(face = "bold", size = 9))

for (ext in c("png", "svg")) {
  out <- file.path(project_root, "results/figures", paste0("fig_sigma_effects_faceted.", ext))
  if (ext == "png") {
    ggplot2::ggsave(out, plot = p4, width = 10, height = 6, dpi = 300, bg = "white")
  } else {
    ggplot2::ggsave(out, plot = p4, width = 10, height = 6, bg = "white")
  }
  file.copy(out, file.path(project_root, "submission/bmc_medicine/figures_export",
                           paste0("figure4.", ext)), overwrite = TRUE)
  cat("saved:", out, "+ figures_export copy\n")
}

# ---------------- Figure 6: 面板标题与轴标签补全 ----------------
B1_RS_PHQ    <- readRDS(file.path(fits_dir, "B1_RS_PHQ.rds"))
B1_RS_PC     <- readRDS(file.path(fits_dir, "B1_RS_PC.rds"))
B1_RS_sigmaM <- readRDS(file.path(fits_dir, "B1_RS_sigmaM.rds"))

sd_phq <- sd(sample_A4$pre_phq9_c, na.rm = TRUE)
sd_pc  <- sd(sample_A4$pre_control_c, na.rm = TRUE)
sd_sm  <- sd(sample_A4$person_std_mood_c, na.rm = TRUE)

configs <- list(
  list(fit = B1_RS_PHQ,    mod = "pre_phq9_c",         sd = sd_phq, title = "Moderator: PHQ-9 (baseline)"),
  list(fit = B1_RS_PC,     mod = "pre_control_c",      sd = sd_pc,  title = "Moderator: Perceived control (baseline)"),
  list(fit = B1_RS_sigmaM, mod = "person_std_mood_c",  sd = sd_sm,  title = "Moderator: EMA-derived mood variability")
)

plots <- list()
for (cfg in configs) {
  ce <- conditional_effects(cfg$fit,
    effects = paste0("health_session_c:", cfg$mod),
    int_conditions = setNames(list(c(-1, 0, 1) * cfg$sd), cfg$mod))
  d <- ce[[1]]
  d$mod_label <- factor(d[[cfg$mod]], labels = c("Low (-1 SD)", "Mean", "High (+1 SD)"))
  p <- ggplot(d, aes(x = health_session_c, y = estimate__, color = mod_label, fill = mod_label)) +
    geom_ribbon(aes(ymin = lower__, ymax = upper__), alpha = 0.15, color = NA) +
    geom_line(linewidth = 1) +
    labs(title = cfg$title,
         x = "Within-person health behavior (session-level, CWC)",
         y = "Predicted mood") +
    theme_minimal(base_size = 10) +
    theme(legend.position = "bottom", legend.title = element_blank(),
          plot.title = element_text(size = 11))
  plots[[length(plots) + 1]] <- p
}

p6 <- wrap_plots(plots, ncol = 3) +
  plot_annotation(title = "Cross-level moderation of the within-person health-mood association")

for (ext in c("png", "svg")) {
  out <- file.path(project_root, "results/figures", paste0("fig_conditional_effects.", ext))
  if (ext == "png") {
    ggplot2::ggsave(out, plot = p6, width = 15, height = 5, dpi = 150, bg = "white")
  } else {
    ggplot2::ggsave(out, plot = p6, width = 15, height = 5, bg = "white")
  }
  file.copy(out, file.path(project_root, "submission/bmc_medicine/figures_export",
                           paste0("figure6.", ext)), overwrite = TRUE)
  cat("saved:", out, "+ figures_export copy\n")
}
cat("done\n")
