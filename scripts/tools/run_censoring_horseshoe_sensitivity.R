#!/usr/bin/env Rscript
# ============================================================================
# Censoring & horseshoe sensitivity batch
#   F1 b1_censored        B1_MELS 规格 + 双侧删失（−10 左删失 / +10 右删失）
#   F2 item_noon_censored 午间条目 Main 规格 + 删失
#   F3 item_eve_censored  晚间条目 Main 规格 + 删失
#   F4 item_noon_horseshoe 午间条目 Main 规格（未删失）+ sigma 子模型 horseshoe
#   F5 item_eve_horseshoe  晚间条目 Main 规格（未删失）+ sigma 子模型 horseshoe
#
# 规格卡：
#   共同数据源 checkpoint_01_data.RData::sample_lagged（6,836 obs / 236 人）
#   与既有 B1_MELS / Item_MELS_*_Main 同源同编码（sleep_ontime_c = -sleep_late_c）
#   先验 = priors_mels（Intercept N(3,2)；b N(0,1)；σ截距 t(3,0,1)；sd t(3,0,1)）
#   F4/F5 仅改 sigma 子模型 class="b" 先验 → horseshoe(df=1)
#   采样 4 chains × iter 6000 / warmup 3000；adapt_delta 0.95（发散>0 升 0.99 并记录）
#   max_treedepth 15（F1）/ 12（F2–F5，同原模型）；seed 42 / 20261026 / 20261027 / 20260519 / 20260519
#
# 用法：
#   R5_MODE=dry   Rscript scripts/tools/run_censoring_horseshoe_sensitivity.R   # 30 min 内验语法
#   R5_STAGE=fit  Rscript ...（默认，断点续跑：已有 .rds 跳过）
#   R5_STAGE=report Rscript ...（只出报告与 PPC 图，读取 .rds）
#   R5_MODELS=b1_censored,item_noon_censored ...（子集）
# ============================================================================
options(width = 200); options(warn = 1)

get_root <- function() {
  a <- commandArgs(trailingOnly = FALSE)
  f <- sub("^--file=", "", a[grep("^--file=", a)])
  if (length(f) > 0 && file.exists(f)) return(normalizePath(file.path(dirname(f), "..", "..")))
  getwd()
}
project_root <- get_root()

suppressPackageStartupMessages({
  library(brms); library(cmdstanr); library(posterior)
  library(bayesplot); library(ggplot2); library(dplyr); library(jsonlite)
})

MODE   <- Sys.getenv("R5_MODE", "production")   # production | dry
STAGE  <- Sys.getenv("R5_STAGE", "fit")         # fit | report | both
DRY    <- identical(MODE, "dry")
ALL_NAMES <- c("b1_censored", "item_noon_censored", "item_eve_censored",
               "item_noon_horseshoe", "item_eve_horseshoe")
SEL_ENV <- Sys.getenv("R5_MODELS", "")
SELECT  <- if (DRY) c("b1_censored", "item_noon_horseshoe") else
           if (nchar(SEL_ENV) > 0) trimws(strsplit(SEL_ENV, ",")[[1]]) else ALL_NAMES

sensitivity_main_dir <- file.path(project_root, "results", "fits", "sensitivity_main")
sensitivity_main_log <- file.path(sensitivity_main_dir, "logs")
ppc_dir    <- file.path(project_root, "docs", "reports", "sensitivity", "ppc")
for (d in c(sensitivity_main_dir, sensitivity_main_log, ppc_dir)) dir.create(d, showWarnings = FALSE, recursive = TRUE)

cat("\n═══ censoring & horseshoe sensitivity ═══\n")
cat(sprintf("  mode=%s stage=%s models=%s\n", MODE, STAGE, paste(SELECT, collapse = ", ")))

# ---------------------------------------------------------------------------
# 1. 数据（只读 checkpoint；与 run_item_mels.R 同一构造路径）
# ---------------------------------------------------------------------------
e <- new.env()
load(file.path(project_root, "data", "checkpoints", "checkpoint_01_data.RData"), envir = e)
sample_lagged <- e$sample_lagged %>%
  mutate(sleep_ontime_c = -sleep_late_c) %>%
  mutate(cens_mood = ifelse(mood <= -10, "left",
                     ifelse(mood >=  10, "right", "none")))
sample_noon <- sample_lagged %>% filter(session_noon == 1) %>% mutate(user_id = droplevels(user_id))
sample_eve  <- sample_lagged %>% filter(session_noon == 0) %>% mutate(user_id = droplevels(user_id))
cat(sprintf("  data: lagged %d | noon %d | eve %d | cens left %d / right %d\n",
            nrow(sample_lagged), nrow(sample_noon), nrow(sample_eve),
            sum(sample_lagged$cens_mood == "left"), sum(sample_lagged$cens_mood == "right")))

# ---------------------------------------------------------------------------
# 2. 规格卡（固定）
# ---------------------------------------------------------------------------
priors_mels <- c(
  brms::prior(normal(0, 1),         class = "b"),
  brms::prior(normal(0, 1),         class = "b",         dpar = "sigma"),
  brms::prior(normal(3, 2),         class = "Intercept"),
  brms::prior(student_t(3, 0, 1),   class = "Intercept", dpar = "sigma"),
  brms::prior(student_t(3, 0, 1),   class = "sd")
)
priors_hs_sigma <- c(
  brms::prior(normal(0, 1),         class = "b"),
  brms::prior(normal(3, 2),         class = "Intercept"),
  brms::prior(student_t(3, 0, 1),   class = "Intercept", dpar = "sigma"),
  brms::prior(student_t(3, 0, 1),   class = "sd"),
  brms::prior(horseshoe(df = 1, par_ratio = NULL), class = "b", dpar = "sigma")
)
# 备选（仅当 dpar 层 horseshoe 不被支持）：全部 population-level b 施加 horseshoe
priors_hs_all <- c(
  brms::prior(horseshoe(df = 1, par_ratio = NULL), class = "b"),
  brms::prior(normal(3, 2),         class = "Intercept"),
  brms::prior(student_t(3, 0, 1),   class = "Intercept", dpar = "sigma"),
  brms::prior(student_t(3, 0, 1),   class = "sd")
)

form_b1 <- brms::bf(
  mood | cens(cens_mood) ~ health_session_c + health_session_between_c +
    mood_lag1_day_c + study_day_c + session_numeric + (1 | user_id),
  sigma ~ health_session_c + health_session_between_c
)
form_noon <- brms::bf(
  mood | cens(cens_mood) ~ sleep_ontime_c + wake_early_c + breakfast_c + aspiration_c +
    mood_lag1_day_c + study_day_c + health_session_between_c + (1 | user_id),
  sigma ~ sleep_ontime_c + wake_early_c + breakfast_c + aspiration_c
)
form_eve <- brms::bf(
  mood | cens(cens_mood) ~ walk_c + exercise_c + dinner_c + breathing_c +
    mood_lag1_day_c + study_day_c + health_session_between_c + (1 | user_id),
  sigma ~ walk_c + exercise_c + dinner_c + breathing_c
)
form_noon_unc <- brms::bf(
  mood ~ sleep_ontime_c + wake_early_c + breakfast_c + aspiration_c +
    mood_lag1_day_c + study_day_c + health_session_between_c + (1 | user_id),
  sigma ~ sleep_ontime_c + wake_early_c + breakfast_c + aspiration_c
)
form_eve_unc <- brms::bf(
  mood ~ walk_c + exercise_c + dinner_c + breathing_c +
    mood_lag1_day_c + study_day_c + health_session_between_c + (1 | user_id),
  sigma ~ walk_c + exercise_c + dinner_c + breathing_c
)

SPECS <- list(
  b1_censored = list(form = form_b1,   data = sample_lagged, prior = priors_mels,
                     family = "gaussian", cens = TRUE,  seed = 42L,       mt = 15L),
  item_noon_censored  = list(form = form_noon, data = sample_noon, prior = priors_mels,
                     family = "gaussian", cens = TRUE,  seed = 20261026L, mt = 12L),
  item_eve_censored   = list(form = form_eve,  data = sample_eve,  prior = priors_mels,
                     family = "gaussian", cens = TRUE,  seed = 20261027L, mt = 12L),
  item_noon_horseshoe = list(form = form_noon_unc, data = sample_noon, prior = priors_hs_sigma,
                     family = "gaussian", cens = FALSE, seed = 20260519L, mt = 12L,
                     prior_label = "horseshoe(df=1) on sigma-submodel b"),
  item_eve_horseshoe  = list(form = form_eve_unc,  data = sample_eve,  prior = priors_hs_sigma,
                     family = "gaussian", cens = FALSE, seed = 20260519L, mt = 12L,
                     prior_label = "horseshoe(df=1) on sigma-submodel b")
)

n_chains <- if (DRY) 2L else 4L
n_iter   <- if (nchar(Sys.getenv("R5_ITER", "")) > 0) as.integer(Sys.getenv("R5_ITER")) else
            if (DRY) 500L else 6000L
n_warmup <- if (nchar(Sys.getenv("R5_WARMUP", "")) > 0) as.integer(Sys.getenv("R5_WARMUP")) else
            if (DRY) 250L else 3000L
ad_start <- if (nchar(Sys.getenv("R5_AD", "")) > 0) as.numeric(Sys.getenv("R5_AD")) else 0.95

# ---------------------------------------------------------------------------
# 3. 审计日志（与 fit_or_load 同格式；results/logs/config_audit_log.jsonl）
# ---------------------------------------------------------------------------
audit_file <- file.path(project_root, "results", "logs", "config_audit_log.jsonl")
audit_append <- function(model, params) {
  entry <- list(timestamp = format(Sys.time(), "%Y%m%d_%H%M%S"),
                framework = "brms", model = paste0("sensitivity_main/", model), params = params)
  cat(toJSON(entry, auto_unbox = TRUE), "\n", file = audit_file, append = file.exists(audit_file))
  cat(sprintf("  [audit] %s @ %s\n", model, entry$timestamp))
}

# ---------------------------------------------------------------------------
# 4. 收敛门
# ---------------------------------------------------------------------------
convergence_gate <- function(fit) {
  s <- posterior::summarise_draws(posterior::as_draws_df(fit),
                                  posterior::default_convergence_measures())
  s <- s[!startsWith(s$variable, "lp"), ]
  np <- brms::nuts_params(fit)
  div <- sum(subset(np, Parameter == "divergent__")$Value > 0)
  list(max_rhat = max(s$rhat, na.rm = TRUE),
       min_ess_bulk = min(s$ess_bulk, na.rm = TRUE),
       min_ess_tail = min(s$ess_tail, na.rm = TRUE),
       divergent = div,
       table = s)
}

# ---------------------------------------------------------------------------
# 5. 拟合（断点续跑 + 发散升 0.99 重跑）
# ---------------------------------------------------------------------------
fit_one <- function(name) {
  sp <- SPECS[[name]]
  rds_path <- if (DRY) file.path(sensitivity_main_dir, paste0("dry_", name, ".rds")) else
                       file.path(sensitivity_main_dir, paste0(name, ".rds"))
  log_path <- file.path(sensitivity_main_log, paste0(name, "_cmdstan.log"))
  if (file.exists(rds_path)) { cat(sprintf("\n[%s] 已有 .rds，跳过\n", name)); return(readRDS(rds_path)) }

  run_fit <- function(ad) {
    cat(sprintf("\n[%s] 拟合中 chains=%d iter=%d warmup=%d adapt_delta=%.2f seed=%d ...\n",
                name, n_chains, n_iter, n_warmup, ad, sp$seed))
    t0 <- Sys.time()
    con <- file(log_path, open = "wt")
    sink(con, type = "message")
    fit <- tryCatch(
      brms::brm(sp$form, data = sp$data, family = brms::brmsfamily("gaussian"),
                prior = sp$prior, chains = n_chains, iter = n_iter, warmup = n_warmup,
                cores = n_chains, seed = sp$seed,
                control = list(adapt_delta = ad, max_treedepth = sp$mt),
                backend = "cmdstanr", save_pars = brms::save_pars(all = TRUE),
                refresh = if (DRY) 100 else 500,
                threads = brms::threading(2L),
                stan_model_args = list(cpp_options = list(stan_threads = TRUE))),
      error = function(e) { cat(sprintf("[%s] ERROR: %s\n", name, conditionMessage(e))); NULL },
      finally = { sink(NULL, type = "message"); close(con) })
    if (!is.null(fit)) attr(fit, "elapsed_min") <- as.numeric(difftime(Sys.time(), t0, units = "mins"))
    fit
  }

  if (!DRY) audit_append(name, list(chains = n_chains, iter = n_iter, warmup = n_warmup,
                          adapt_delta = ad_start, max_treedepth = sp$mt, seed = sp$seed,
                          family = if (sp$cens) "gaussian (two-sided censored at -10/+10)" else "gaussian",
                          prior = sp$prior_label %||% "priors_mels (as pipeline)"))
  fit <- run_fit(ad_start)
  if (is.null(fit)) return(NULL)

  gate <- convergence_gate(fit)
  cat(sprintf("  [%s] gate: Rhat=%.4f ess_bulk=%.0f div=%d\n",
              name, gate$max_rhat, gate$min_ess_bulk, gate$divergent))
  if (!DRY && gate$divergent > 0 && ad_start < 0.99) {
    cat(sprintf("  [%s] divergent=%d > 0 → adapt_delta 0.95 → 0.99 重跑（记录）\n", name, gate$divergent))
    saveRDS(fit, file.path(sensitivity_main_dir, paste0(name, "_ad95_div.rds")), compress = "xz")
    audit_append(name, list(chains = n_chains, iter = n_iter, warmup = n_warmup,
                            adapt_delta = 0.99, max_treedepth = sp$mt, seed = sp$seed,
                            note = "escalated from 0.95 due to divergent transitions"))
    fit2 <- run_fit(0.99)
    if (!is.null(fit2)) fit <- fit2
  }
  if (!is.null(fit)) saveRDS(fit, rds_path, compress = "xz")
  fit
}
`%||%` <- function(a, b) if (is.null(a)) b else a

if (STAGE %in% c("fit", "both")) {
  fits <- list()
  for (nm in SELECT) {
    f <- fit_one(nm)
    if (!is.null(f)) {
      g <- convergence_gate(f)
      fits[[nm]] <- list(fit = f, gate = g,
                         pass = !DRY && g$max_rhat < 1.01 && g$min_ess_bulk > 400 && g$divergent == 0)
      cat(sprintf("[%s] PASS=%s (Rhat=%.4f ess_bulk=%.0f div=%d)\n", nm,
                  fits[[nm]]$pass, g$max_rhat, g$min_ess_bulk, g$divergent))
    } else {
      cat(sprintf("[%s] FIT FAILED\n", nm))
    }
  }
  # dry 模式 PPC 管道验证（censored fit 的 posterior_predict 可用性）
  if (DRY && !is.null(fits$b1_censored)) {
    yr <- tryCatch(posterior_predict(fits$b1_censored$fit, ndraws = 20),
                   error = function(e) { cat("PPC ERROR:", conditionMessage(e), "\n"); NULL })
    cat(sprintf("[dry] PPC pipeline: %s | yrep dim %s\n",
                if (is.null(yr)) "FAIL" else "OK",
                if (is.null(yr)) "-" else paste(dim(yr), collapse = "x")))
  }
}

# ---------------------------------------------------------------------------
# 6. 报告（R5_STAGE=report 可独立运行）
# ---------------------------------------------------------------------------
summarize_par <- function(fit, v) {
  dr <- posterior::as_draws_df(fit)[[v]]
  if (is.null(dr)) return(rep(NA_real_, 3))
  c(median(dr), quantile(dr, .025, na.rm = TRUE), quantile(dr, .975, na.rm = TRUE))
}
fmt3 <- function(x) sprintf("%.3f [%.3f, %.3f]", x[1], x[2], x[3])

if (STAGE %in% c("report", "both")) {
  cat("\n═══ 报告生成 ═══\n")
  item_rds_dir <- file.path(project_root, "legacy", "bayesian_analysis_",
                            "modeling_results_", "item_mels_fits")
  B1_orig    <- readRDS(file.path(project_root, "results", "fits", "B1_MELS.rds"))
  Noon_orig  <- readRDS(file.path(item_rds_dir, "Item_MELS_Noon_Main.rds"))
  Eve_orig   <- readRDS(file.path(item_rds_dir, "Item_MELS_Evening_Main.rds"))
  F1 <- readRDS(file.path(sensitivity_main_dir, "b1_censored.rds"))
  F2 <- readRDS(file.path(sensitivity_main_dir, "item_noon_censored.rds"))
  F3 <- readRDS(file.path(sensitivity_main_dir, "item_eve_censored.rds"))
  F4 <- readRDS(file.path(sensitivity_main_dir, "item_noon_horseshoe.rds"))
  F5 <- readRDS(file.path(sensitivity_main_dir, "item_eve_horseshoe.rds"))

  gate_all <- function(fit) {
    g <- convergence_gate(fit)
    sprintf("R̂ max %.4f, bulk ESS min %.0f, divergent %d", g$max_rhat, g$min_ess_bulk, g$divergent)
  }

  # ---- censored_comparison.md ----
  rows <- list(
    c("composite mu (health_session_c)",
      fmt3(summarize_par(B1_orig, "b_health_session_c")),
      fmt3(summarize_par(F1, "b_health_session_c"))),
    c("composite log-sigma (health_session_c)",
      fmt3(summarize_par(B1_orig, "b_sigma_health_session_c")),
      fmt3(summarize_par(F1, "b_sigma_health_session_c"))),
    c("composite log-sigma (between)",
      fmt3(summarize_par(B1_orig, "b_sigma_health_session_between_c")),
      fmt3(summarize_par(F1, "b_sigma_health_session_between_c"))),
    c("aspiration mu",      fmt3(summarize_par(Noon_orig, "b_aspiration_c")), fmt3(summarize_par(F2, "b_aspiration_c"))),
    c("aspiration log-sigma", fmt3(summarize_par(Noon_orig, "b_sigma_aspiration_c")), fmt3(summarize_par(F2, "b_sigma_aspiration_c"))),
    c("exercise mu",        fmt3(summarize_par(Eve_orig, "b_exercise_c")),   fmt3(summarize_par(F3, "b_exercise_c"))),
    c("exercise log-sigma", fmt3(summarize_par(Eve_orig, "b_sigma_exercise_c")), fmt3(summarize_par(F3, "b_sigma_exercise_c")))
  )
  noon_items <- c("sleep_ontime_c", "wake_early_c", "breakfast_c", "aspiration_c")
  eve_items  <- c("walk_c", "exercise_c", "dinner_c", "breathing_c")

  direction_verdict <- function(orig_str, new_str) {
    nums_of <- function(s) as.numeric(regmatches(s, gregexpr("-?[0-9]+\\.?[0-9]*", s))[[1]])
    o <- nums_of(orig_str); n <- nums_of(new_str)
    if (length(o) < 3 || length(n) < 3 || is.na(n[1])) return("n/a")
    if (sign(o[1]) == sign(n[1]) && n[1] != 0) {
      excl0 <- (n[2] > 0 | n[3] < 0)
      return(if (excl0) "方向保持（CrI 不含 0）" else "方向保持（CrI 含 0，强度减弱）")
    }
    "方向改变 ⚠"
  }
  hdr <- c("| 对比项 | 原正态模型 (median [95% CrI]) | censored 模型 | 判定 |",
           "|---|---|---|---|")
  body <- ""
  key_rows <- rows[1:7]
  for (i in seq_along(key_rows)) {
    r <- key_rows[[i]]
    body <- paste0(body, "| ", r[1], " | ", r[2], " | ", r[3], " | ",
                   direction_verdict(r[2], r[3]), " |\n")
  }
  all_items <- c(paste0("noon ", noon_items), paste0("evening ", eve_items))
  full_body <- "\n\n### 全八条目 mu 与 log-sigma（censored vs 原正态）\n\n| 条目 | 子模型 | 原 (Main) | censored |\n|---|---|---|---|\n"
  for (it in noon_items) for (dp in c("mu", "sigma")) {
    v <- if (dp == "mu") paste0("b_", it) else paste0("b_sigma_", it)
    full_body <- paste0(full_body, "| noon ", it, " | ", dp, " | ",
                        fmt3(summarize_par(Noon_orig, v)), " | ", fmt3(summarize_par(F2, v)), " |\n")
  }
  for (it in eve_items) for (dp in c("mu", "sigma")) {
    v <- if (dp == "mu") paste0("b_", it) else paste0("b_sigma_", it)
    full_body <- paste0(full_body, "| evening ", it, " | ", dp, " | ",
                        fmt3(summarize_par(Eve_orig, v)), " | ", fmt3(summarize_par(F3, v)), " |\n")
  }
  conv_body <- paste0("\n\n### 收敛门槛（R̂ < 1.01, bulk ESS > 400, divergent = 0）\n\n",
    "| fit | 诊断 |\n|---|---|\n",
    "| F1 b1_censored | ", gate_all(F1), " |\n",
    "| F2 item_noon_censored | ", gate_all(F2), " |\n",
    "| F3 item_eve_censored | ", gate_all(F3), " |\n",
    "| F4 item_noon_horseshoe | ", gate_all(F4), " |\n",
    "| F5 item_eve_horseshoe | ", gate_all(F5), " |\n")
  md <- c("# Censored sensitivity comparison", "",
          "双侧删失 Gaussian location-scale：mood = −10 左删失、+10 右删失（Tobit 边界似然），",
          "规格与 B1_MELS / Item_MELS_*_Main 完全一致（仅 likelihood 改变）；先验 priors_mels。",
          "seed：F1 42（iter 8000/4000）/ F2 20261026 / F3 20261027；",
          "adapt_delta 0.95（censored 三 fit 无发散，未升级）。",
          "F5 item_eve_horseshoe 未过 div=0 门（0.999 下 2/12000 残余发散，全部尝试留档），",
          "其系数仅作披露性参考，详见 horseshoe_shrinkage.md。",
          "", hdr, body, full_body, conv_body,
          "", sprintf("生成 %s；脚本 scripts/tools/run_censoring_horseshoe_sensitivity.R", format(Sys.time())))
  writeLines(md, file.path(project_root, "docs", "reports", "sensitivity", "censored_comparison.md"))
  cat("  ✓ censored_comparison.md\n")

  # ---- horseshoe_shrinkage.md ----
  hs_body <- "| 时段 | 条目 | 原后验中位数 [CrI] | horseshoe 后中位数 [CrI] | 是否收缩向 0 | CrI 是否含 0 |\n|---|---|---|---|---|---|\n"
  for (it in noon_items) {
    o <- summarize_par(Noon_orig, paste0("b_sigma_", it)); n <- summarize_par(F4, paste0("b_sigma_", it))
    shrunk <- if (is.na(n[1])) "n/a" else if (abs(n[1]) < abs(o[1])) sprintf("是（|%.3f|→|%.3f|）", o[1], n[1]) else "否"
    excl <- if (is.na(n[1])) "n/a" else if (n[2] > 0 | n[3] < 0) "不含 0" else "含 0"
    hs_body <- paste0(hs_body, "| 午间 | ", it, " | ", fmt3(o), " | ", fmt3(n), " | ", shrunk, " | ", excl, " |\n")
  }
  for (it in eve_items) {
    o <- summarize_par(Eve_orig, paste0("b_sigma_", it)); n <- summarize_par(F5, paste0("b_sigma_", it))
    shrunk <- if (is.na(n[1])) "n/a" else if (abs(n[1]) < abs(o[1])) sprintf("是（|%.3f|→|%.3f|）", o[1], n[1]) else "否"
    excl <- if (is.na(n[1])) "n/a" else if (n[2] > 0 | n[3] < 0) "不含 0" else "含 0"
    hs_body <- paste0(hs_body, "| 晚间 | ", it, " | ", fmt3(o), " | ", fmt3(n), " | ", shrunk, " | ", excl, " |\n")
  }
  hs_md <- c("# Horseshoe shrinkage sensitivity", "",
             "对午/晚条目模型 sigma 子模型的 4 个系数施加 horseshoe(df = 1) 收缩先验",
             "（mu 子模型保持 normal(0,1)；规格与原 Main 模型一致，seed 20260519）。",
             "判定关注点：aspiration 的 −0.365 与 exercise 的 +0.121 在收缩后是否保持方向与不含 0。",
             "",
             "**门状态**：F4 午间 = PASS（adapt_delta 0.999；0.95/0.99 尝试各留档 `_ad95_div.rds` / `_ad99_div16.rds`）；",
             "F5 晚间 = **未过 div=0 门**（failure_type: convergence_error）：0.95→200 div、0.99→4 div、",
             "0.999→2 div、0.999+iter 8000→16 div，全部尝试留档（`item_eve_horseshoe_ad*.rds`）；",
             "下表 F5 列为 0.999 + 4×6000 尝试（残余 2/12000 发散）的披露性参考值，R̂/ESS 达标但不可宣称完全通过。",
             "", hs_body, "",
             paste0("F4 收敛：", gate_all(F4), "；F5 收敛：", gate_all(F5)),
             "", sprintf("生成 %s；脚本 scripts/tools/run_censoring_horseshoe_sensitivity.R", format(Sys.time())))
  writeLines(hs_md, file.path(project_root, "docs", "reports", "sensitivity", "horseshoe_shrinkage.md"))
  cat("  ✓ horseshoe_shrinkage.md\n")

  # ---- PPC 图（白底，DECISIONS 2026-10-26-B 标准）----
  make_ppc <- function(fit, data, name, label) {
    yr <- posterior_predict(fit, ndraws = 100)
    p_ge10 <- apply(yr >= 10, 1, mean); p_le10m <- apply(yr <= -10, 1, mean)
    obs_ge10 <- mean(data$mood >= 10); obs_le <- mean(data$mood <= -10)
    df <- data.frame(x = as.vector(yr[sample.int(nrow(yr), min(100, nrow(yr))), ]))
    p <- ggplot(df, aes(x = x)) +
      geom_density(color = "steelblue", alpha = 0.1, fill = "steelblue") +
      geom_density(data = data, aes(x = mood), color = "black", inherit.aes = FALSE,
                   linetype = "dashed") +
      geom_vline(xintercept = c(-10, 10), color = "grey40", linetype = "dotted") +
      annotate("text", x = 0, y = Inf, vjust = 2,
               label = sprintf("P(y_rep ≥ +10) = %.1f%% (obs %.1f%%)\nP(y_rep ≤ −10) = %.1f%% (obs %.1f%%)",
                               100 * mean(p_ge10), 100 * obs_ge10, 100 * mean(p_le10m), 100 * obs_le)) +
      labs(title = paste0("PPC density overlay (100 draws): ", label),
           subtitle = "dashed = observed; dotted = scale bounds; latent predictive draws",
           x = "mood", y = "density") +
      theme_minimal(base_size = 11)
    ggsave(file.path(ppc_dir, paste0(name, "_ppc.png")), plot = p,
           width = 8, height = 5, dpi = 150, bg = "white")
    cat(sprintf("  ✓ ppc/%s_ppc.png (P≥+10 %.1f%% vs obs %.1f%% | P≤−10 %.1f%% vs obs %.1f%%)\n",
                name, 100 * mean(p_ge10), 100 * obs_ge10, 100 * mean(p_le10m), 100 * obs_le))
  }
  make_ppc(F1, sample_lagged, "f1_b1_censored", "F1 composite censored")
  make_ppc(F2, sample_noon,   "f2_noon_censored", "F2 noon items censored")
  make_ppc(F3, sample_eve,    "f3_eve_censored", "F3 evening items censored")
  make_ppc(F4, sample_noon,   "f4_noon_horseshoe", "F4 noon horseshoe")
  make_ppc(F5, sample_eve,    "f5_eve_horseshoe", "F5 evening horseshoe")
  cat("\n═══ 报告完成 ═══\n")
}
