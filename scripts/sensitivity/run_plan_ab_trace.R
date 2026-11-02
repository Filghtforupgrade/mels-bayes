#!/usr/bin/env Rscript
# ============================================================================
# A3: Plan A vs Plan B trace 对比图 — Dora Fig.3 对标
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

cat("\n═══ A3: Plan A vs Plan B trace 对比 ═══\n")
t0 <- Sys.time()
suppressPackageStartupMessages({
  library(ggplot2); library(bayesplot); library(patchwork)
})

# ---- 1. Plan A: 故意失败版 (sigma 含随机截距) ----
PlanA_form <- brms::bf(
  mood  ~ health_session_c + health_session_between_c +
          mood_lag1_day_c + study_day_c + session_numeric +
          (1 | user_id),
  sigma ~ health_session_c + health_session_between_c +
          (1 | user_id)   # ← 这会导致发散
)

cat("--- 拟合 Plan A (sigma 含 RE，预期不收敛) ---\n")
set.seed(20260511)
B1_MELS_PlanA <- brms::brm(
  PlanA_form, data = sample_lagged, family = brmsfamily("gaussian"),
  prior = priors_mels, chains = 4L, iter = 2000L, warmup = 1000L,
  cores = min(4L, getOption("mc.cores")), seed = 20260511,
  control = list(adapt_delta = 0.95, max_treedepth = 12L),
  backend = "cmdstanr", save_pars = save_pars(all = TRUE),
  refresh = 0
)
planA_s <- extract_posterior_summary(B1_MELS_PlanA)
planA_rhat <- max(planA_s$rhat, na.rm = TRUE)
planA_ess  <- min(planA_s$ess_bulk, na.rm = TRUE)
planA_np   <- brms::nuts_params(B1_MELS_PlanA)
planA_ndiv <- sum(subset(planA_np, Parameter == "divergent__")$Value)
cat(sprintf("  Plan A: max R-hat = %.3f, min ESS = %.0f (%d divergent)\n", planA_rhat, planA_ess, planA_ndiv))

# ---- 2. 加载 Plan B ( 主跑) ----
B1_MELS_PlanB <- readRDS(file.path(fits_dir, "B1_MELS.rds"))
planB_s <- extract_posterior_summary(B1_MELS_PlanB)
planB_rhat <- max(planB_s$rhat, na.rm = TRUE)
planB_ess  <- min(planB_s$ess_bulk, na.rm = TRUE)
cat(sprintf("  Plan B: max R-hat = %.3f, min ESS = %.0f\n", planB_rhat, planB_ess))

# ---- 3. 4-panel trace grid ----
make_trace <- function(fit, param, title_str, plan_s) {
  r <- plan_s[plan_s$variable == param, ]
  rhat_v <- if(nrow(r) > 0) r$rhat else NA
  p <- bayesplot::mcmc_trace(fit, pars = param) +
    ggplot2::ggtitle(sprintf("%s (R-hat=%.3f)", title_str, rhat_v)) +
    ggplot2::theme_minimal(base_size = 9) +
    ggplot2::theme(legend.position = "none")
  p
}

panels <- list()

pa_sigma <- "b_sigma_health_session_c"
if (pa_sigma %in% brms::variables(B1_MELS_PlanA)) {
  panels[[1]] <- make_trace(B1_MELS_PlanA, pa_sigma,
                             "Plan A: sigma_health (FAILED)", planA_s)
}
pb_sigma <- "b_sigma_health_session_c"
if (pb_sigma %in% brms::variables(B1_MELS_PlanB)) {
  panels[[2]] <- make_trace(B1_MELS_PlanB, pb_sigma,
                             "Plan B: sigma_health (PASSED)", planB_s)
}
pa_mu <- "b_health_session_c"
if (pa_mu %in% brms::variables(B1_MELS_PlanA)) {
  panels[[3]] <- make_trace(B1_MELS_PlanA, pa_mu,
                             "Plan A: health mu", planA_s)
}
pb_mu <- "b_health_session_c"
if (pb_mu %in% brms::variables(B1_MELS_PlanB)) {
  panels[[4]] <- make_trace(B1_MELS_PlanB, pb_mu,
                             "Plan B: health mu", planB_s)
}

if (length(panels) >= 2) {
  final <- wrap_plots(panels, ncol = 2) +
    plot_annotation(
      title = "Plan A (sigma with random intercept, FAILED) vs Plan B (fixed-only sigma, PASSED)",
      subtitle = sprintf("Plan A: 4 chains × 2000 iter; Plan B: 4 chains × 6000 iter"),
      theme = theme(plot.title = element_text(size = 13, face = "bold"))
    )
  
  out_path <- file.path(figures_dir, "fig_plan_ab_trace_comparison.png")
  ggsave(out_path, plot = final, width = 12, height = 8, dpi = 150)
  cat(sprintf("\n  ✓ %s\n", basename(out_path)))
}

elapsed <- difftime(Sys.time(), t0, units = "mins")
cat(sprintf("═══ A3 完成 — 耗时 %.1f 分钟 ═══\n", as.numeric(elapsed)))
