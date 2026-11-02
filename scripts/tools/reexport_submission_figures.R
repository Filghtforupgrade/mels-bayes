#!/usr/bin/env Rscript
# ============================================================================
# 主图重导出（BMC Medicine submission-ready）
# ============================================================================
# 策略：临时接管 ggplot2::ggsave ——
#   1) 所有输出重定向 results/figures_reexport_tmp/（results/figures 既有文件零触碰）
#   2) 每图额外输出 _hr.png（dpi 缩放至 ≥2400×1600 目标）与同名 .svg
#   3) 源脚本（08_figures_main / 08_figures_paper / run_conditional_effects）原样 source，零修改
# 产出：submission/bmc_medicine/figures_export/figure{1..6}.{png,svg}
# 图号映射（执行时双重确认后落 mapping 表）：
#   Fig1 = fig_prior_distributions   (Methods §2.5 先验校准语境)
#   Fig2 = fig_trace_grid             (Methods §2.7 收敛语境)
#   Fig3 = fig_ppcheck_mels           (Methods §2.8 PPC density overlay)
#   Fig4 = fig_sigma_effects          (Results §3.3)
#   Fig5 = fig_individual_slopes     (Results §3.4)
#   Fig6 = fig_conditional_effects   (Results §3.4)
# ============================================================================
options(warn = 1)

project_root <- normalizePath(getwd())
tmp_dir   <- file.path(project_root, "results", "figures_reexport_tmp")
export_dir <- file.path(project_root, "submission", "bmc_medicine", "figures_export")
dir.create(tmp_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(export_dir, recursive = TRUE, showWarnings = FALSE)

# ---- setup（提供 figures_dir / fits_dir / project_root 等全局）+ 数据/先验 ----
# SCRIPT_DIR 必须预定义：setup 以 SCRIPT_DIR/../.. 推项目根，回退 getwd() 会错指家目录
SCRIPT_DIR <- file.path(project_root, "scripts", "pipeline")
source(file.path(project_root, "scripts", "modules", "00_setup.R"))
source(file.path(project_root, "scripts", "modules", "01_data.R"))
source(file.path(project_root, "scripts", "modules", "02_priors.R"))
# 03_mels：经 fit_or_load 直接 readRDS 装载 B1/R3_MELS（磁盘存在不重跑 MCMC），
# 填充 brms_fits 与 mels_sigma_results（Fig2/3/4 的依赖）
source(file.path(project_root, "scripts", "modules", "03_mels.R"))

# ---- ggsave 接管 ----
ori_ggsave <- ggplot2::ggsave
hr_ggsave <- function(filename, plot = ggplot2::last_plot(), ..., dpi = 300) {
  args <- list(...)
  bn <- basename(filename)
  base_noext <- sub("\\.png$", "", bn)
  # 原分辨率（重定向至 tmp）
  do.call(ori_ggsave, c(list(filename = file.path(tmp_dir, bn), plot = plot), args))
  w <- args$width; h <- args$height
  if (is.null(w) || is.null(h)) return(invisible(NULL))
  scale <- max(1, 2400 / (w * dpi), 1600 / (h * dpi))
  dpi_hr <- round(dpi * scale)
  args_hr <- args; args_hr$units <- args$units %||% "in"
  # 高分辨率 PNG
  ok_png <- tryCatch({
    do.call(ori_ggsave, c(list(filename = file.path(tmp_dir, paste0(base_noext, "_hr.png")),
                            plot = plot), args_hr, list(dpi = dpi_hr))); TRUE
  }, error = function(e) { cat(sprintf("    [hr-png 失败] %s: %s\n", bn, conditionMessage(e))); FALSE })
  # SVG 矢量
  ok_svg <- tryCatch({
    do.call(ori_ggsave, c(list(filename = file.path(tmp_dir, paste0(base_noext, ".svg")),
                            plot = plot, device = svglite::svglite), args_hr)); TRUE
  }, error = function(e) { cat(sprintf("    [svg 失败] %s: %s\n", bn, conditionMessage(e))); FALSE })
  invisible(list(png_hr = ok_png, svg = ok_svg))
}
`%||%` <- function(a, b) if (is.null(a)) b else a
assignInNamespace("ggsave", hr_ggsave, ns = "ggplot2")
# 全局环境遮蔽：拦截裸 ggsave(...) 调用（解析顺序 globalenv 先于 package:ggplot2；

ggsave <- hr_ggsave

# ---- 依次 source 三个源脚本（不修改；单个失败不阻塞其余）----
src_scripts <- c(
  main = file.path(project_root, "scripts", "modules", "08_figures_main.R"),
  paper = file.path(project_root, "scripts", "modules", "08_figures_paper.R"),
  cond = file.path(project_root, "scripts", "sensitivity", "run_conditional_effects.R")
)
for (nm in names(src_scripts)) {
  cat(sprintf("\n=== source [%s]: %s ===\n", nm, basename(src_scripts[[nm]])))
  tryCatch(source(src_scripts[[nm]], encoding = "UTF-8", echo = FALSE),
           error = function(e) cat(sprintf("  !! [%s] 失败: %s\n", nm, conditionMessage(e))))
}

# 恢复原 ggsave
assignInNamespace("ggsave", ori_ggsave, ns = "ggplot2")

# ---- 提升 6 张目标图到 figures_export ----
fig_map <- c(
  figure1 = "fig_prior_distributions",
  figure2 = "fig_trace_grid",
  figure3 = "fig_ppcheck_mels",
  figure4 = "fig_sigma_effects",
  figure5 = "fig_individual_slopes",
  figure6 = "fig_conditional_effects"
)
cat("\n=== 提升目标图 ===\n")
for (i in seq_along(fig_map)) {
  out_name <- names(fig_map)[i]
  src_hr <- file.path(tmp_dir, paste0(fig_map[[i]], "_hr.png"))
  src_svg <- file.path(tmp_dir, paste0(fig_map[[i]], ".svg"))
  dst_png <- file.path(export_dir, paste0(out_name, ".png"))
  dst_svg <- file.path(export_dir, paste0(out_name, ".svg"))
  if (file.exists(src_hr)) { file.copy(src_hr, dst_png, overwrite = TRUE); cat(sprintf("  %s.png <- %s_hr.png\n", out_name, fig_map[[i]])) }
  else cat(sprintf("  !! %s.png 缺失（%s）\n", out_name, basename(src_hr)))
  if (file.exists(src_svg)) { file.copy(src_svg, dst_svg, overwrite = TRUE); cat(sprintf("  %s.svg <- %s.svg\n", out_name, fig_map[[i]])) }
  else cat(sprintf("  !! %s.svg 缺失\n", out_name))
}
cat("\n完成。产物目录:", export_dir, "\n")
