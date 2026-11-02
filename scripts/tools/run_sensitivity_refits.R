#!/usr/bin/env Rscript
# ============================================================================
# run_sensitivity_refits.R 
# ============================================================================
# 子命令:
#   f1  Plan A illustrative 重跑（sigma 随机截距版，adapt_delta 0.99/td 15，
#       run_plan_ab_trace.R 原值 0.95/12）
#       → results/fits/sensitivity_refits/plan_a_illustrative.rds + t5_plan_a_batches.md
#   f2  item-level prior-only 拟合 ×4（含 Noon_Main 统一口径）
#       → 4 个 rds + t7_prior_ppc.md（PASS/WARN/FAIL 沿用房内 ≥0.70/≥0.10 阈值）
#   f3  B1_RS_PHQ flat refit（16,000/8,000 对齐主拟合终版；priors_rs_flat 全类 SD≈10，
#       与 B1_MELS_flat_prior 同约定）→ b1_rs_phq_flat_prior.rds + t6_flat_s6.md/csv
#   f4  administered-session 完成率（纯计算，三条 stopifnot）→ t12_administered_rate.md
#
# 本脚本不修改任何现有拟合脚本与稿件；产物仅入 results/fits/sensitivity_refits/ 与
# config_audit_log.jsonl 追加登记。
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
if (!file.exists(file.path(project_root, "scripts", "modules", "00_setup.R"))) {
  project_root <- "."  # released fallback: run from repo root  # 兜底：后台任务 CWD 不确定时用绝对路径
}
# 00_setup.R 依赖：SCRIPT_DIR 未定义时回退 getwd() 并上跳两级（scripts/pipeline 假设），
# 从 tools/ 调用会得到错误根目录 —— 必须在 source 前显式定义
SCRIPT_DIR <- file.path(project_root, "scripts", "pipeline")
MODULE_DIR <- file.path(project_root, "scripts", "modules")

source(file.path(MODULE_DIR, "00_setup.R"))
source(file.path(MODULE_DIR, "01_data.R"))
source(file.path(MODULE_DIR, "02_priors.R"))

suppressPackageStartupMessages({
  library(brms)
})

R5B_DIR  <- file.path(project_root, "results", "fits", "sensitivity_refits")
R5B_LOG  <- file.path(R5B_DIR, "logs")
RPT_DIR  <- file.path(project_root, "docs", "reports", "sensitivity")
AUDIT_LOG <- file.path(project_root, "results", "logs", "config_audit_log.jsonl")
for (d in c(R5B_DIR, R5B_LOG, RPT_DIR)) dir.create(d, showWarnings = FALSE, recursive = TRUE)

append_audit <- function(model, params, change, note) {
  ts <- format(Sys.time(), "%Y%m%d_%H%M%S")
  if (requireNamespace("jsonlite", quietly = TRUE)) {
    line <- jsonlite::toJSON(
      list(timestamp = ts, framework = "brms", model = model,
           params = params, change = change, note = note),
      auto_unbox = TRUE)
  } else {
    line <- sprintf('{"timestamp": "%s", "framework": "brms", "model": "%s", "note": "%s"}', ts, model, note)
  }
  cat(as.character(line), "\n", file = AUDIT_LOG, append = TRUE)
  cat("  [config_audit] 已登记:", model, "\n")
}

sink_log <- function(name) {
  con <- file(file.path(R5B_LOG, paste0(name, "_console.log")), open = "wt")
  sink(con, type = "output"); sink(con, type = "message")
  con
}

cmd <- commandArgs(trailingOnly = TRUE)
if (length(cmd) == 0) stop("用法: Rscript run_sensitivity_refits.R {f1|f2|f3|f4}")

# ============================================================================
# F1 — Plan A illustrative 重跑（T5）
# ============================================================================
if (identical(cmd, "f1")) {
  cat("═══ F1: Plan A illustrative 重跑（adapt_delta 0.99 / treedepth 15, seed 20260511）═══\n")
  t0 <- Sys.time()

  PlanA_form <- brms::bf(
    mood  ~ health_session_c + health_session_between_c +
            mood_lag1_day_c + study_day_c + session_numeric +
            (1 | user_id),
    sigma ~ health_session_c + health_session_between_c +
            (1 | user_id)
  )

  fit <- brms::brm(
    PlanA_form, data = sample_lagged, family = brmsfamily("gaussian"),
    prior = priors_mels, chains = 4L, iter = 2000L, warmup = 1000L,
    cores = min(4L, getOption("mc.cores")), seed = 20260511,
    control = list(adapt_delta = 0.99, max_treedepth = 15L),
    backend = "cmdstanr", save_pars = save_pars(all = TRUE),
    refresh = 0
  )
  saveRDS(fit, file.path(R5B_DIR, "plan_a_illustrative.rds"))
  cat("  ✓ rds 落盘: plan_a_illustrative.rds\n")

  s <- extract_posterior_summary(fit)
  pick <- function(v) {
    r <- s[s$variable == v, ]
    if (nrow(r) == 0) return(c(rhat = NA, ess = NA)) else c(rhat = r$rhat, ess = r$ess_bulk)
  }
  p_fail  <- pick("sd_user_id__sigma_Intercept")
  p_disp  <- pick("b_sigma_health_session_c")
  p_mu    <- pick("b_health_session_c")
  np  <- brms::nuts_params(fit)
  div <- sum(subset(np, Parameter == "divergent__")$Value)
  n_draws <- 4 * 1000

  out <- data.frame(
    quantity = c("max R-hat (all params)", "min bulk ESS (all params)",
                 "sd_user_id__sigma_Intercept R-hat", "sd_user_id__sigma_Intercept ESS_bulk",
                 "b_sigma_health_session_c R-hat", "b_health_session_c R-hat",
                 "divergent transitions", "post-warmup draws"),
    value = c(max(s$rhat, na.rm = TRUE), min(s$ess_bulk, na.rm = TRUE),
              p_fail[["rhat"]], p_fail[["ess"]],
              p_disp[["rhat"]], p_mu[["rhat"]],
              div, n_draws)
  )
  write.csv(out, file.path(R5B_DIR, "t5_diagnostics.csv"), row.names = FALSE)
  print(out)

  append_audit("PlanA_illustrative_sensitivity_refits",
    list(chains = 4, iter = 2000, warmup = 1000, adapt_delta = 0.99,
         max_treedepth = 15, seed = 20260511),
    list(fit = list(old = "none (illustrative refit, first on-disk save)",
                    new = "results/fits/sensitivity_refits/plan_a_illustrative.rds")),
    "Plan A sigma-RE illustrative refit; adapt_delta 0.99/td 15, aligned with the primary B1_RS_PHQ final rerun spec (run_plan_ab_trace.R original 0.95/12); expected divergent geometry is the diagnostic product, not a failure")

  cat(sprintf("═══ F1 完成 — 耗时 %.1f 分钟 ═══\n",
              as.numeric(difftime(Sys.time(), t0, units = "mins"))))
}

# ============================================================================
# F2 — item-level prior-only 拟合 ×4（T7）
# ============================================================================
if (identical(cmd, "f2")) {
  cat("═══ F2: item-level prior predictive checks（sample_prior=only, seed 20260519）═══\n")
  t0 <- Sys.time()

  ckpt <- file.path(project_root, "data", "checkpoints", "checkpoint_01_data.RData")
  load(ckpt)
  sample_lagged$sleep_ontime_c <- -sample_lagged$sleep_late_c
  sample_noon <- droplevels(subset(sample_lagged, session_noon == 1))
  sample_eve  <- droplevels(subset(sample_lagged, session_noon == 0))
  cat(sprintf("  样本: noon %d obs / %d 人, eve %d obs / %d 人\n",
              nrow(sample_noon), nlevels(sample_noon$user_id),
              nrow(sample_eve), nlevels(sample_eve$user_id)))

  forms <- list(
    Item_MELS_Noon_Main = list(
      form = bf(mood ~ sleep_ontime_c + wake_early_c + breakfast_c + aspiration_c +
                       mood_lag1_day_c + study_day_c + health_session_between_c + (1 | user_id),
                sigma ~ sleep_ontime_c + wake_early_c + breakfast_c + aspiration_c),
      data = sample_noon, rds = "item_noon_main_prior_only.rds"),
    Item_MELS_Noon_Robust = list(
      form = bf(mood ~ sleep_ontime_c + wake_early_c + breakfast_c + aspiration_c +
                       mood_lag1_day_c + mood_prev_session_c + study_day_c + health_session_between_c + (1 | user_id),
                sigma ~ sleep_ontime_c + wake_early_c + breakfast_c + aspiration_c),
      data = sample_noon, rds = "item_noon_robust_prior_only.rds"),
    Item_MELS_Evening_Main = list(
      form = bf(mood ~ walk_c + exercise_c + dinner_c + breathing_c +
                       mood_lag1_day_c + study_day_c + health_session_between_c + (1 | user_id),
                sigma ~ walk_c + exercise_c + dinner_c + breathing_c),
      data = sample_eve, rds = "item_eve_main_prior_only.rds"),
    Item_MELS_Evening_Robust = list(
      form = bf(mood ~ walk_c + exercise_c + dinner_c + breathing_c +
                       mood_lag1_day_c + mood_prev_session_c + study_day_c + health_session_between_c + (1 | user_id),
                sigma ~ walk_c + exercise_c + dinner_c + breathing_c),
      data = sample_eve, rds = "item_eve_robust_prior_only.rds")
  )

  rows <- list()
  for (nm in names(forms)) {
    cfg <- forms[[nm]]
    cat(sprintf("\n  [%s] prior-only 拟合...\n", nm))
    fit <- brms::brm(
      cfg$form, data = cfg$data, family = brmsfamily("gaussian"),
      prior = priors_mels, sample_prior = "only",
      chains = 4L, iter = 2000L, warmup = 1000L,
      cores = min(4L, getOption("mc.cores")), seed = 20260519,
      control = list(adapt_delta = 0.95, max_treedepth = 12L),
      backend = "cmdstanr", save_pars = save_pars(all = TRUE), refresh = 0
    )
    saveRDS(fit, file.path(R5B_DIR, cfg$rds))

    s <- extract_posterior_summary(fit)
    max_rhat <- max(s$rhat, na.rm = TRUE)
    y_prior <- brms::posterior_predict(fit, ndraws = 500)
    pct_in <- mean(y_prior >= -10 & y_prior <= 10, na.rm = TRUE)
    status <- if (pct_in >= 0.70) "PASS" else if (pct_in >= 0.10) "WARN" else "FAIL"
    rows[[nm]] <- data.frame(
      model = nm, pct_in_range = round(pct_in, 4),
      median_mean = round(median(rowMeans(y_prior, na.rm = TRUE)), 2),
      median_sd = round(median(apply(y_prior, 1, sd, na.rm = TRUE)), 2),
      max_rhat_priorfit = round(max_rhat, 4), status = status)
    print(rows[[nm]])

    append_audit(paste0(nm, "_prior_only"),
      list(chains = 4, iter = 2000, warmup = 1000, adapt_delta = 0.95,
           max_treedepth = 12, seed = 20260519, sample_prior = "only"),
      list(fit = list(old = "none", new = paste0("results/fits/sensitivity_refits/", cfg$rds))),
      paste0("Item-level prior predictive check; metric follows modules/02_priors.R prior_predictive_check (pct of prior-predictive mood in [-10,10]; PASS>=0.70 WARN>=0.10)"))
  }
  tbl <- do.call(rbind, rows)
  write.csv(tbl, file.path(R5B_DIR, "t7_prior_ppc.csv"), row.names = FALSE)
  cat("\n  ✓ t7_prior_ppc.csv 落盘\n")
  cat(sprintf("═══ F2 完成 — 耗时 %.1f 分钟 ═══\n",
              as.numeric(difftime(Sys.time(), t0, units = "mins"))))
}

# ============================================================================
# F3 — B1_RS_PHQ flat refit（T6；16,000/8,000 对齐主拟合终版）
# ============================================================================
if (identical(cmd, "f3")) {
  cat("═══ F3: B1_RS_PHQ flat-prior refit（4×16,000/8,000, ad 0.99, seed 20260511）═══\n")
  t0 <- Sys.time()

  # flat 档与 B1_MELS_flat_prior 同约定：全类 SD≈10（b/Intercept/sd/sigma），cor 保持 lkj(2)
  priors_rs_flat <- c(
    brms::prior(normal(0, 10),        class = "b"),
    brms::prior(normal(3, 10),        class = "Intercept"),
    brms::prior(student_t(3, 0, 10),  class = "sd"),
    brms::prior(student_t(3, 0, 10),  class = "sigma"),
    brms::prior(lkj(2),               class = "cor")
  )

  B1_RS_PHQ_form <- bf(
    mood ~ health_session_c * pre_phq9_c +
           health_session_between_c +
           mood_prev_session_c + study_day_c + session_numeric +
           (1 + health_session_c | user_id)
  )

  fit <- brms::brm(
    B1_RS_PHQ_form, data = sample_A4, family = brmsfamily("gaussian"),
    prior = priors_rs_flat,
    chains = 4L, iter = 16000L, warmup = 8000L,
    cores = min(4L, getOption("mc.cores")), seed = 20260511,
    control = list(adapt_delta = 0.99, max_treedepth = 15L),
    backend = "cmdstanr", save_pars = save_pars(all = TRUE),
    threads = brms::threading(2L),
    stan_model_args = list(cpp_options = list(stan_threads = TRUE)),
    sample_prior = "yes", refresh = 0
  )
  saveRDS(fit, file.path(R5B_DIR, "b1_rs_phq_flat_prior.rds"))
  cat("  ✓ rds 落盘: b1_rs_phq_flat_prior.rds\n")

  s <- extract_posterior_summary(fit)
  np  <- brms::nuts_params(fit)
  div <- sum(subset(np, Parameter == "divergent__")$Value)
  diag <- data.frame(
    quantity = c("max R-hat", "min bulk ESS", "divergent transitions", "post-warmup draws"),
    value = c(max(s$rhat, na.rm = TRUE), min(s$ess_bulk, na.rm = TRUE), div, 4 * 8000))
  write.csv(diag, file.path(R5B_DIR, "t6_diagnostics.csv"), row.names = FALSE)
  print(diag)

  focal <- c("b_health_session_c", "b_health_session_c:pre_phq9_c",
             "b_mood_prev_session_c", "sd_user_id__health_session_c")
  fs <- s[s$variable %in% focal, ]
  fs$median <- round(fs$median, 3)
  fs$ci <- sprintf("%.3f, %.3f", fs$ci_low, fs$ci_high)
  flat_out <- data.frame(parameter = fs$variable,
                         flat_median = fs$median, flat_95crI = fs$ci,
                         flat_excludes_zero = fs$ci_low > 0 | fs$ci_high < 0)
  write.csv(flat_out, file.path(R5B_DIR, "t6_flat_s6.csv"), row.names = FALSE)
  print(flat_out)

  append_audit("B1_RS_PHQ_flat_prior",
    list(chains = 4, iter = 16000, warmup = 8000, adapt_delta = 0.99,
         max_treedepth = 15, threads = 2, seed = 20260511,
         prior_b = "normal(0,10)", prior_sd = "student_t(3,0,10)"),
    list(prior = list(old = "priors_rs (tight: b normal(0,1), sd student_t(3,0,1))",
                      new = "flat (SD approx 10, same convention as B1_MELS_flat_prior)")),
    "Flat-prior refit to complete the three-tier ladder for the 4 B1_RS_PHQ focal parameters; sampler spec aligned with primary B1_RS_PHQ final rerun (16000/8000)")

  cat(sprintf("═══ F3 完成 — 耗时 %.1f 分钟 ═══\n",
              as.numeric(difftime(Sys.time(), t0, units = "mins"))))
}

# ============================================================================
# F4 — administered-session 完成率（T12；纯计算）
# ============================================================================
if (identical(cmd, "f4")) {
  cat("═══ F4: administered-session 完成率 ═══\n")
  ckpt <- file.path(project_root, "data", "checkpoints", "checkpoint_01_data.RData")
  load(ckpt)

  df <- session_df[, c("user_id", "date", "mood")]
  df <- df[!is.na(df$mood), ]
  df$date <- as.Date(df$date)
  stopifnot(nrow(df) == 8309)

  per <- do.call(rbind, lapply(split(df, df$user_id), function(d) {
    f <- min(d$date); l <- max(d$date)
    data.frame(user_id = as.character(d$user_id[1]),
               completed = nrow(d), span_days = as.numeric(l - f),
               first = as.character(f), last = as.character(l))
  }))
  per$administered <- pmin((per$span_days + 1) * 2, 48)

  stopifnot(all(per$administered <= 48))                       # ①
  stopifnot(all(per$completed <= per$administered))             # ②
  # ③ 提案原文断言 max(span)=21 与数据不符（实测 max=23，中位 21；24 个完成 44 次者
  #   23 人 span=21→(21+1)×2=44、1 人 span=22）。改为数据锚定断言：完成上限与 S11 的
  #   max 44 一致，且完成数不超过首末活跃日推得的不封顶槽位。分母定义未改（D18 登记）。
  stopifnot(max(per$completed) == 44)
  stopifnot(all(per$completed <= (per$span_days + 1) * 2))

  total_admin <- sum(per$administered)
  rate_admin <- 8309 / total_admin
  rate_proto <- 8309 / 11328

  out <- data.frame(
    metric = c("completed sessions (total)", "participants",
               "protocol denominator (236 x 48)", "protocol completion rate",
               "administered denominator (sum)", "administered completion rate",
               "administered slots per person: median [min,max]",
               "span days per person: median [min,max]"),
    value = c(8309, nrow(per), 11328, sprintf("%.1f%%", 100 * rate_proto),
              total_admin, sprintf("%.1f%%", 100 * rate_admin),
              sprintf("%d [%d, %d]", median(per$administered), min(per$administered), max(per$administered)),
              sprintf("%d [%d, %d]", median(per$span_days), min(per$span_days), max(per$span_days))))
  write.csv(out, file.path(R5B_DIR, "t12_administered_rate.csv"), row.names = FALSE)
  write.csv(per, file.path(R5B_DIR, "t12_per_person.csv"), row.names = FALSE)
  print(out)
  cat("  ✓ 三条 stopifnot 全过；t12 csv 落盘\n")
}

cat("\n[run_sensitivity_refits] 完成:", paste(cmd, collapse = ","), "\n")
