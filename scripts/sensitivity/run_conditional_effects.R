#!/usr/bin/env Rscript
# ============================================================================
# A5:  conditional_effects() + 个体斜率图（对标 Dora Fig.9/12）
#  B1: §3 绘图段已与 run_conditional_effects_popavg.R 统一
#   （re_formula = NA；标题/副标逐字一致）；两脚本互为镜像，改一处须同步另一处
# ============================================================================
options(width = 200); options(warn = 1)

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
SCRIPT_DIR    <- file.path(project_root, "scripts", "pipeline")
MODULE_DIR <- file.path(project_root, "scripts", "modules")

source(file.path(MODULE_DIR, "00_setup.R"))
source(file.path(MODULE_DIR, "01_data.R"))
source(file.path(MODULE_DIR, "02_priors.R"))

suppressPackageStartupMessages({
  library(ggplot2)
  library(patchwork)
  library(scales)
})

cat("\n═══ A5: conditional_effects + 个体斜率图 ═══\n")
t0 <- Sys.time()

# ---- 1. 加载所有 RS 模型 ----
B1_RS_PHQ    <- readRDS(file.path(fits_dir, "B1_RS_PHQ.rds"))
B1_RS_PC     <- readRDS(file.path(fits_dir, "B1_RS_PC.rds"))
B1_RS_sigmaM <- readRDS(file.path(fits_dir, "B1_RS_sigmaM.rds"))
B1_RS        <- readRDS(file.path(fits_dir, "B1_RS.rds"))
cat("  ✓ 4 models loaded\n")

# ---- 2. conditional_effects 3-panel ----
sd_phq  <- sd(sample_A4$pre_phq9_c, na.rm = TRUE)
sd_pc   <- sd(sample_A4$pre_control_c, na.rm = TRUE)
sd_sm   <- sd(sample_A4$person_std_mood_c, na.rm = TRUE)

configs <- list(
  list(fit = B1_RS_PHQ, mod = "pre_phq9_c", sd = sd_phq, title = "× PHQ-9"),
  list(fit = B1_RS_PC,  mod = "pre_control_c", sd = sd_pc, title = "× Perceived Control"),
  list(fit = B1_RS_sigmaM, mod = "person_std_mood_c", sd = sd_sm, title = "× mood variability")
)

plots <- list()
for (cfg in configs) {
  ce <- conditional_effects(
    cfg$fit,
    effects = paste0("health_session_c:", cfg$mod),
    int_conditions = setNames(list(c(-1, 0, 1) * cfg$sd), cfg$mod)
  )
  d <- ce[[1]]
  d$mod_label <- factor(d[[cfg$mod]], labels = c("Low (-1 SD)", "Mean", "High (+1 SD)"))
  p <- ggplot(d, aes(x = health_session_c, y = estimate__,
                      color = mod_label, fill = mod_label)) +
    geom_ribbon(aes(ymin = lower__, ymax = upper__), alpha = 0.15, color = NA) +
    geom_line(linewidth = 1) +
    labs(title = cfg$title, x = "Health behavior (CWC)", y = "Predicted Mood") +
    theme_minimal(base_size = 10) +
    theme(legend.position = "bottom", legend.title = element_blank())
  plots[[length(plots) + 1]] <- p
}

p_all <- wrap_plots(plots, ncol = 3) +
  plot_annotation(
    title = "Cross-Level Moderation: Conditional Effects"
  )

out_ce <- file.path(figures_dir, "fig_conditional_effects.png")
ggsave(out_ce, plot = p_all, width = 15, height = 5, dpi = 150, bg = "white")
cat(sprintf("  ✓ %s\n", basename(out_ce)))

# ---- 3. 个体斜率图（Dora Fig.9 对标； B1：本节与 run_conditional_effects_popavg.R 逐字一致） ----
set.seed(1)
sample_ids <- sample(unique(B1_RS$data$user_id), 12)

# 固定效应总体平均带（re_formula = NA，不含 group-level 效应；与 popavg 脚本一致）
re <- conditional_effects(B1_RS, effects = "health_session_c",
                          re_formula = NA)
rd <- re$health_session_c

p_slopes <- ggplot(rd, aes(x = effect1__, y = estimate__)) +
  geom_ribbon(aes(ymin = lower__, ymax = upper__), alpha = 0.15, fill = "steelblue") +
  geom_smooth(
    data = B1_RS$data[B1_RS$data$user_id %in% sample_ids, ],
    aes(x = health_session_c, y = mood, group = user_id),
    method = "lm", se = FALSE, linewidth = 0.3,
    color = alpha("grey40", 0.6)
  ) +
  geom_line(color = "black", linewidth = 1.2) +
  labs(
    title = "Health-to-Mood Association for 12 Participants",
    subtitle = paste0("Grey lines: descriptive within-person OLS fits; ",
                      "black line: population-average slope with 95% credible band (fixed effects only)"),
    x = "Health Behavior (within-person, CWC)",
    y = "Mood (-10 to +10)"
  ) +
  theme_minimal(base_size = 11)

out_slopes <- file.path(figures_dir, "fig_individual_slopes.png")
ggsave(out_slopes, plot = p_slopes, width = 10, height = 7, dpi = 150, bg = "white")
cat(sprintf("  ✓ %s\n", basename(out_slopes)))

elapsed <- difftime(Sys.time(), t0, units = "mins")
cat(sprintf("\n═══ A5 完成 — 耗时 %.1f 分钟 ═══\n", as.numeric(elapsed)))
