#!/usr/bin/env Rscript
# ============================================================================
#  B1: Figure 5 population-average 重绘
# 复制自 run_conditional_effects.R 第 3 节（个体斜率图）。
#  B1: 原脚本 §3 已并入本脚本口径（re_formula = NA + 同标题/副标），
#   两脚本互为镜像，改一处须同步另一处。差异仅剩输出路径（本脚本直出 figures_export）。
# 差异：
#   1. re_formula = ~ 1 + health_session_c | user_id  →  re_formula = NA
#      （可信带只含固定效应不确定性，与图注 population-average 一致）
#   2. 标题去 "B1_RS:" 模型名前缀，改内容式标题
#   3. 灰线保留（描述性 within-person OLS），副标标注其性质
#   4. 直接导出 submission/bmc_medicine/figures_export/figure5.{png,svg}
#      （不写 results/figures/，既有图零触碰）
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

suppressPackageStartupMessages({
  library(ggplot2)
})

cat("\n═══ Figure 5 population-average 重绘 ═══\n")
t0 <- Sys.time()

B1_RS <- readRDS(file.path(fits_dir, "B1_RS.rds"))
cat("  ✓ B1_RS loaded\n")

set.seed(1)
sample_ids <- sample(unique(B1_RS$data$user_id), 12)

# 固定效应总体平均带（不含 group-level 效应）
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

export_dir <- file.path(project_root, "submission", "bmc_medicine", "figures_export")
out_png <- file.path(export_dir, "figure5.png")
out_svg <- file.path(export_dir, "figure5.svg")
ggsave(out_png, plot = p_slopes, width = 10, height = 7, dpi = 240, bg = "white")
ggsave(out_svg, plot = p_slopes, width = 10, height = 7, bg = "white")
cat(sprintf("  ✓ %s\n", basename(out_png)))
cat(sprintf("  ✓ %s\n", basename(out_svg)))

elapsed <- difftime(Sys.time(), t0, units = "mins")
cat(sprintf("\n═══ 完成 — 耗时 %.1f 分钟 ═══\n", as.numeric(elapsed)))
