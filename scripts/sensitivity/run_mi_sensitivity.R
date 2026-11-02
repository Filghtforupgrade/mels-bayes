#!/usr/bin/env Rscript
# ============================================================================
# A14:  多重插补稳健性检查 — 对标 Dora et al. (2024) Example 3
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

cat("\n═══ A14: 多重插补稳健性检查 — M=10 ═══\n")
t0 <- Sys.time()
suppressPackageStartupMessages({ library(mice); library(ggplot2) })

# ---- 阶段 1: L2 插补 ----
cat("\n--- 阶段 1: mice L2 插补 (m=10, maxit=20) ---\n")
person_df <- sample_B %>%
  group_by(user_id) %>%
  summarise(
    person_mean_mood   = mean(mood, na.rm = TRUE),
    person_mean_health = mean(health_session_c, na.rm = TRUE),
    person_std_mood    = sd(mood, na.rm = TRUE),
    n_sessions = n(), n_days = n_distinct(study_day_c),
    .groups = "drop"
  ) %>%
  left_join(sample_A4 %>% distinct(user_id, pre_phq9_c, pre_control_c), by = "user_id")

imp_vars <- c("pre_phq9_c", "pre_control_c",
              "person_mean_mood", "person_mean_health",
              "person_std_mood", "n_sessions", "n_days")
cat(sprintf("  Person-level data: %d rows, %d with missing PHQ9/PC\n",
            nrow(person_df), sum(is.na(person_df$pre_phq9_c))))

set.seed(20260525)
imp <- mice(person_df[, imp_vars], m = 10, maxit = 20, method = "pmm", printFlag = FALSE)

# 诊断：密度叠加图
png(file.path(figures_dir, "fig_mi_density.png"), width = 8, height = 5, units = "in", res = 150)
densityplot(imp)
dev.off()
cat("  ✓ MI density plot saved\n")

# 生成 10 套插补后的 sample_A4 数据集
imputed_samples <- lapply(1:10, function(i) {
  completed_person <- complete(imp, i)
  completed_person$user_id <- person_df$user_id
  # 去掉 sample_B 中已有的缺失列，避免 left_join 产生 .x/.y 后缀
  sb <- sample_B[, setdiff(names(sample_B), c("pre_phq9_c", "pre_control_c"))]
  sb <- sb %>%
    left_join(completed_person %>% select(user_id, pre_phq9_c, pre_control_c), by = "user_id")
  # 删除此处的二次 scale()。
  # 01b 已把 sample_A4 的 pre_phq9_c/pre_control_c 标准化为 SD=1，person_df
  # 由 sample_A4 取值、mice 在该标准化尺度上插补，故插补结果已是 SD≈1 尺度，
  # 与 listwise 主线 (sample_A4) 同尺度。原 scale() 会按"含插补的全样本"重标准化，
  # 引入 listwise(person-weighted) vs MI 的尺度错配——正是 F4 要消除的问题。
  sb
})
cat(sprintf("  ✓ %d imputed datasets generated (%d rows each)\n", length(imputed_samples), nrow(imputed_samples[[1]])))

# person_std_mood_c 仅存在于 sample_A4，需手动计算并 join
# mood 无缺失 → person_std_mood 不依赖 MI，直接从 person_df 取
person_std_mood_c_vec <- as.numeric(scale(person_df$person_std_mood))
pstd_df <- data.frame(user_id = person_df$user_id,
                       person_std_mood_c = person_std_mood_c_vec)
imputed_samples <- lapply(imputed_samples, function(sb) {
  merge(sb, pstd_df, by = "user_id", all.x = TRUE)
})
cat(sprintf("  ✓ person_std_mood_c joined (%d NA)\n",
            sum(is.na(imputed_samples[[1]]$person_std_mood_c))))

# ---- 阶段 2: brm_multiple ----
cat("\n--- 阶段 2: brm_multiple on 3 models ---\n")

# 定义模型公式（本地，避免 source 04b 触发重跑）
B1_RS_PHQ_form <- brms::bf(
  mood ~ health_session_c * pre_phq9_c + health_session_between_c +
         mood_prev_session_c + study_day_c + session_numeric +
         (1 + health_session_c | user_id)
)
B1_RS_PC_form <- brms::bf(
  mood ~ health_session_c * pre_control_c + health_session_between_c +
         mood_prev_session_c + study_day_c + session_numeric +
         (1 + health_session_c | user_id)
)
B1_RS_sigmaM_form <- brms::bf(
  mood ~ health_session_c * person_std_mood_c + health_session_between_c +
         mood_prev_session_c + study_day_c + session_numeric +
         (1 + health_session_c | user_id)
)

fit_mi_manual <- function(form, fname, imputed_data) {
  out_rds <- file.path(fits_dir, paste0(fname, "_mi.rds"))
  if (file.exists(out_rds)) {
    cat(sprintf("  [%s] loaded from .rds\n", fname))
    return(readRDS(out_rds))
  }
  cat(sprintf("  [%s] fitting across 10 imputed datasets ...\n", fname))
  t1 <- Sys.time()
  #  增强：每个 set 单独 checkpoint 到 set_k.rds，
  # 这样 process 中途死掉时，下次启动可跳过已完成的 set。
  ckpt_dir <- file.path(fits_dir, paste0(fname, "_mi_checkpoints"))
  dir.create(ckpt_dir, showWarnings = FALSE, recursive = TRUE)

  # 逐套拟合 brms（避免 brm_multiple 的后端兼容问题）
  fit_list <- list()
  for (k in seq_along(imputed_data)) {
    set_rds <- file.path(ckpt_dir, sprintf("set_%02d.rds", k))
    if (file.exists(set_rds)) {
      cat(sprintf("    set %d/10... [loaded from checkpoint]\n", k))
      fit_list[[k]] <- readRDS(set_rds)
      next
    }
    cat(sprintf("    set %d/10...", k))
    fit_k <- brms::brm(
      form, data = imputed_data[[k]], family = gaussian(), prior = priors_rs,
      chains = 2, iter = 2000, warmup = 1000,
      sample_prior = "yes", save_pars = save_pars(all = TRUE),
      control = list(adapt_delta = 0.99, max_treedepth = 15),
      backend = "cmdstanr", refresh = 0, silent = 2
    )
    saveRDS(fit_k, set_rds, compress = "xz")   # per-set checkpoint
    fit_list[[k]] <- fit_k
    cat(" done\n")
  }
  # 合并后验 draws（按 Dora 建议）
  all_draws <- lapply(fit_list, as_draws_df)
  combined_draws <- do.call(rbind, all_draws)
  # 必须在 fit_list[[1]]$fit <- NULL 之前
  # 计算 rhats_list，否则后续 extract_posterior_summary(fit_list[[1]]) 会调用
  # as_draws_df(brmsfit)，而 $fit 已被置空 → "is.stanfit(x) is not TRUE"。
  rhats_list <- lapply(fit_list, function(f) {
    s <- extract_posterior_summary(f); setNames(s$rhat, s$variable)
  })
  fit_list[[1]]$fit <- NULL  # 移除 Stan 对象节省空间
  fit_list[[1]]$combined_draws <- combined_draws
  fit_list[[1]]$rhats_list <- rhats_list
  saveRDS(fit_list[[1]], out_rds, compress = "xz")
  elapsed <- difftime(Sys.time(), t1, units = "mins")
  cat(sprintf("  [%s] ✓ %.1f min\n", fname, as.numeric(elapsed)))
  fit_list[[1]]
}

B1_RS_PHQ_mi    <- fit_mi_manual(B1_RS_PHQ_form,    "B1_RS_PHQ",    imputed_samples)
B1_RS_PC_mi     <- fit_mi_manual(B1_RS_PC_form,     "B1_RS_PC",     imputed_samples)
B1_RS_sigmaM_mi <- fit_mi_manual(B1_RS_sigmaM_form, "B1_RS_sigmaM", imputed_samples)

# ---- 阶段 3: Rhat 逐套检查 (Dora p.14 假阳性) ----
cat("\n--- 阶段 3: Rhat 逐集检查 ---\n")
for (mname in c("B1_RS_PHQ", "B1_RS_PC", "B1_RS_sigmaM")) {
  fit_obj <- get(paste0(mname, "_mi"))
  if (!is.null(fit_obj$rhats_list)) {
    # 合并为矩阵：参数行 × 插补集列
    rhat_mat <- do.call(rbind, lapply(fit_obj$rhats_list, function(x) x[names(fit_obj$rhats_list[[1]])]))
    rhat_max <- apply(rhat_mat, 2, max, na.rm = TRUE)
    cat(sprintf("  %s: max R-hat per set = %s\n",
                mname, paste(round(rhat_max, 3), collapse = ", ")))
    if (any(rhat_max > 1.01)) warning(sprintf("  ⚠ %s: some sets Rhat > 1.01", mname))
  }
}

# ---- 阶段 4: listwise vs MI 对比 ----
cat("\n--- 阶段 4: listwise vs MI 对比 ---\n")
B1_RS_PHQ_lw    <- readRDS(file.path(fits_dir, "B1_RS_PHQ.rds"))
B1_RS_PC_lw     <- readRDS(file.path(fits_dir, "B1_RS_PC.rds"))
B1_RS_sigmaM_lw <- readRDS(file.path(fits_dir, "B1_RS_sigmaM.rds"))

compare_mi <- function(fit_lw, fit_mi, params, model_name) {
  result <- data.frame()
  for (p in params) {
    s_lw <- tryCatch(fixef(fit_lw)[p, ], error = function(e) NULL)
    # MI: 使用已保存的 combined_draws
    draws_mi <- fit_mi$combined_draws
    if (is.null(draws_mi)) next
    b_p <- paste0("b_", p)
    if (is.null(s_lw) || !(b_p %in% names(draws_mi))) next
    draws_lw <- as_draws_df(fit_lw)
    mi_median <- median(draws_mi[[b_p]], na.rm = TRUE)
    mi_q2.5   <- quantile(draws_mi[[b_p]], 0.025, na.rm = TRUE)
    mi_q97.5  <- quantile(draws_mi[[b_p]], 0.975, na.rm = TRUE)
    P_lw <- if(b_p %in% names(draws_lw)) mean(draws_lw[[b_p]] > 0, na.rm = TRUE) else NA
    P_mi <- mean(draws_mi[[b_p]] > 0, na.rm = TRUE)
    result <- rbind(result, data.frame(
      model = model_name, parameter = p,
      listwise_median = s_lw["Estimate"], listwise_q2.5 = s_lw["Q2.5"], listwise_q97.5 = s_lw["Q97.5"],
      mi_median = mi_median, mi_q2.5 = mi_q2.5, mi_q97.5 = mi_q97.5,
      listwise_P_pos = P_lw, mi_P_pos = P_mi,
      delta_median = mi_median - s_lw["Estimate"]
    ))
  }
  result$conclusion_consistent <- with(result,
    (listwise_q2.5 > 0 | listwise_q97.5 < 0) == (mi_q2.5 > 0 | mi_q97.5 < 0))
  result
}

cmp_phq    <- compare_mi(B1_RS_PHQ_lw,    B1_RS_PHQ_mi,
                          c("health_session_c", "health_session_c:pre_phq9_c"), "B1_RS_PHQ")
cmp_pc     <- compare_mi(B1_RS_PC_lw,     B1_RS_PC_mi,
                          c("health_session_c", "health_session_c:pre_control_c"), "B1_RS_PC")
cmp_sigmaM <- compare_mi(B1_RS_sigmaM_lw, B1_RS_sigmaM_mi,
                          c("health_session_c", "health_session_c:person_std_mood_c"), "B1_RS_sigmaM")

mi_comparison <- rbind(cmp_phq, cmp_pc, cmp_sigmaM)
print(mi_comparison[, c("model","parameter","listwise_median","mi_median","delta_median","conclusion_consistent")],
      row.names = FALSE, digits = 3)

write.csv(mi_comparison, file.path(output_dir, "mi_comparison.csv"), row.names = FALSE)

out_md <- file.path(tables_dir, "table_13_mi_sensitivity.md")
md_lines <- c(
  "# Table S13: Multiple Imputation Sensitivity (M=10)", "",
  "| Model | Parameter | Listwise median [CrI] | MI median [CrI] | Δ | Conclusion consistent? |",
  "|:------|:----------|:---------------------|:----------------|--:|:----------------------:|"
)
for (i in seq_len(nrow(mi_comparison))) {
  r <- mi_comparison[i, ]
  lw_str <- sprintf("%.3f [%.3f, %.3f]", r$listwise_median, r$listwise_q2.5, r$listwise_q97.5)
  mi_str <- sprintf("%.3f [%.3f, %.3f]", r$mi_median, r$mi_q2.5, r$mi_q97.5)
  md_lines <- c(md_lines, sprintf("| %s | %s | %s | %s | %.3f | %s |",
              r$model, r$parameter, lw_str, mi_str, r$delta_median,
              ifelse(r$conclusion_consistent, "✓ yes", "⚠ no")))
}
writeLines(md_lines, out_md)
cat(sprintf("\n\n  ✓ %s\n", basename(out_md)))

# ---- Forest plot: listwise vs MI ----
if (nrow(mi_comparison) > 0) {
  long <- rbind(
    data.frame(method = "Listwise", parameter = mi_comparison$parameter,
               model = mi_comparison$model,
               median = mi_comparison$listwise_median,
               q2.5 = mi_comparison$listwise_q2.5,
               q97.5 = mi_comparison$listwise_q97.5),
    data.frame(method = "MI (M=10)", parameter = mi_comparison$parameter,
               model = mi_comparison$model,
               median = mi_comparison$mi_median,
               q2.5 = mi_comparison$mi_q2.5,
               q97.5 = mi_comparison$mi_q97.5)
  )
  long$label <- paste0(long$model, "\n", long$parameter)
  long$label <- factor(long$label, levels = rev(unique(long$label)))
  p_forest <- ggplot(long, aes(x = median, y = label, color = method)) +
    geom_vline(xintercept = 0, linetype = "dashed", color = "grey50") +
    geom_errorbarh(aes(xmin = q2.5, xmax = q97.5),
                   position = position_dodge(width = 0.5), height = 0.3) +
    geom_point(position = position_dodge(width = 0.5), size = 2) +
    scale_color_manual(values = c("Listwise" = "steelblue", "MI (M=10)" = "darkorange")) +
    labs(title = "Listwise vs Multiple Imputation",
         x = "Posterior median (95% CrI)", y = NULL, color = "Method") +
    theme_minimal(base_size = 11) + theme(legend.position = "top")
  ggsave(file.path(figures_dir, "fig_mi_comparison.png"),
         plot = p_forest, width = 10, height = 5, dpi = 150, bg = "white")
  cat("  ✓ fig_mi_comparison.png\n")
}

elapsed <- difftime(Sys.time(), t0, units = "mins")
cat(sprintf("\n═══ A14 完成 — 耗时 %.1f 分钟 ═══\n", as.numeric(elapsed)))
