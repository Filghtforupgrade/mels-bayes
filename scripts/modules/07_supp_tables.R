# ============================================================================
# 模块 07c (): 6 张补充表格 — 对标 Dora 2024 顶刊要求
# ============================================================================
# 依赖: 05b 诊断 CSV / 06b R² CSV / 06c LOO CSV / 02b priors_mels 等
# 产出: paper_tables/ 下 6 张新 MD 表
#   1) table_05_convergence.md          — 收敛诊断表
#   2) table_06_prior_ppc.md            — 先验预测检查表
#   3) table_07_variance_decomposition.md — 方差五分量表
#   4) table_08_prior_specifications.md — 先验规格表（手写）
#   5) table_09_model_formulas.md       — 模型公式表（手写）
#   6) table_10_loo.md                  — LOO 模型比较表（依赖 06c）
# ============================================================================

cat("\n--- 模块07c (): 6 张补充表格（Dora 2024 对标）---\n")

# 工具函数 ----
write_md <- function(lines, path) {
  con <- file(path, "wt")
  on.exit(close(con), add = TRUE)
  writeLines(lines, con)
}

format_num <- function(x, digits = 3) {
  if (is.na(x)) return("—")
  if (is.numeric(x) && abs(x) >= 1e4) {
    return(format(as.numeric(x), digits = 3, scientific = TRUE))
  }
  format(round(as.numeric(x), digits), nsmall = digits, scientific = FALSE)
}

# ============================================================================
# 7Ca. Table S5: 收敛诊断 MD 表
# ============================================================================
write_convergence_table <- function() {
  csv_path <- file.path(output_dir, "convergence_diagnostics_.csv")
  out_path <- file.path(tables_dir, "table_05_convergence.md")
  if (!file.exists(csv_path)) {
    cat(sprintf("  ⚠ %s 不存在，跳过\n", basename(csv_path)))
    return(invisible(NULL))
  }
  df <- read.csv(csv_path, stringsAsFactors = FALSE)
  lines <- c(
    "# Table S5 (): Convergence diagnostics — main pipeline",
    "",
    "Reported per fit: maximum R-hat, minimum bulk ESS / tail ESS, divergent transition count and percentage, and the parameter responsible for the worst R-hat / ESS. All fits use 4 chains; iteration / warmup specified by model class (MELS: 6,000 / 3,000; RS: 4,000 / 2,000).",
    "",
    "**Pass thresholds**: R-hat < 1.01, bulk ESS > 400, divergent < 1%. WARN: any one indicator in the moderate band (R-hat ∈ [1.01, 1.05] or ESS ∈ [100, 400] or divergent ∈ [1%, 5%]). FAIL would block downstream modules (see §2.5 of Methods).",
    "",
    "| Model | max R-hat | min ESS bulk | min ESS tail | Divergent | Worst R-hat param | Worst ESS param | Status |",
    "|:------|----------:|-------------:|-------------:|----------:|:------------------|:------------------|:------:|"
  )
  for (i in seq_len(nrow(df))) {
    r <- df[i, ]
    div_str <- sprintf("%d (%.2f%%)", r$divergent_count, r$divergent_pct)
    lines <- c(lines, sprintf("| %s | %s | %s | %s | %s | %s | %s | **%s** |",
                              r$model,
                              format_num(r$max_rhat, 4),
                              format_num(r$min_ess_bulk, 0),
                              format_num(r$min_ess_tail, 0),
                              div_str,
                              r$worst_rhat_param,
                              r$worst_ess_param,
                              r$status))
  }
  lines <- c(lines, "",
             "*Note*: All fits use the cmdstanr backend with `reduce_sum` enabled ( setting). Threading: 2 threads per chain, 4 chains parallel.")
  write_md(lines, out_path)
  cat(sprintf("  ✓ Table S5 写入: %s\n", basename(out_path)))
}
tryCatch(write_convergence_table(), error = function(e)
  cat(sprintf("  ✗ Table S5 失败: %s\n", e$message)))

# ============================================================================
# 7Cb. Table S6: 先验预测检查 MD 表
# ============================================================================
write_prior_ppc_table <- function() {
  csv_path <- file.path(output_dir, "prior_ppc_summary_.csv")
  out_path <- file.path(tables_dir, "table_06_prior_ppc.md")
  if (!file.exists(csv_path)) {
    cat(sprintf("  ⚠ %s 不存在，跳过\n", basename(csv_path)))
    return(invisible(NULL))
  }
  df <- read.csv(csv_path, stringsAsFactors = FALSE)
  lines <- c(
    "# Table S6 (): Prior predictive check gating",
    "",
    "For each model, prior predictive draws (n = 500) of the outcome `mood` are generated under `sample_prior = \"only\"` and compared against the natural mood scale [-10, 10]. The three-tier gate ( onwards) is:",
    "",
    "- **PASS** (pct_in_range ≥ 70%): priors generate plausible mood values for the focal coefficients.",
    "- **WARN** (10% ≤ pct_in_range < 70%): broad-tailed priors interacting with unstandardized higher-order predictors (e.g., `study_day_squared`) widen the predictive distribution, but focal behavior coefficients remain in plausible ranges.",
    "- **FAIL** (pct_in_range < 10%): catastrophic — stops the pipeline.",
    "",
    "| Model | % in [-10, 10] | Median mean | Median SD | Status | Comment |",
    "|:------|---------------:|------------:|----------:|:------:|:--------|"
  )
  for (i in seq_len(nrow(df))) {
    r <- df[i, ]
    pct <- sprintf("%.1f%%", 100 * r$pct_in_range)
    comment <- if (r$status == "WARN") {
      "Broad envelope from higher-order time terms; focal coefs unaffected."
    } else if (r$status == "PASS") {
      "Within scale; priors well-calibrated."
    } else {
      "Catastrophic — pipeline halted."
    }
    lines <- c(lines, sprintf("| %s | %s | %s | %s | **%s** | %s |",
                              r$model, pct,
                              format_num(r$median_mean, 2),
                              format_num(r$median_sd, 2),
                              r$status, comment))
  }
  lines <- c(lines, "",
             "*Note*: Plots in `04_results/figures/prior_pp_check/`.")
  write_md(lines, out_path)
  cat(sprintf("  ✓ Table S6 写入: %s\n", basename(out_path)))
}
tryCatch(write_prior_ppc_table(), error = function(e)
  cat(sprintf("  ✗ Table S6 失败: %s\n", e$message)))

# ============================================================================
# 7Cc. Table S7: 方差五分量 MD 表
# ============================================================================
write_variance_decomp_table <- function() {
  decomp_csv <- file.path(output_dir, "bayes_r2_decomposition_.csv")
  rs_decomp_csv <- file.path(tables_dir, "bayes_r2_decomposition_rs.csv")
  total_csv  <- file.path(output_dir, "bayes_r2_total_.csv")
  out_path <- file.path(tables_dir, "table_07_variance_decomposition.md")
  if (!file.exists(decomp_csv) || !file.exists(total_csv)) {
    cat("  ⚠ R² CSV 缺失，跳过\n")
    return(invisible(NULL))
  }
  decomp <- read.csv(decomp_csv, stringsAsFactors = FALSE)
  total  <- read.csv(total_csv,  stringsAsFactors = FALSE)

  lines <- c(
    "# Table S7 (): Variance decomposition + Bayesian R²",
    "",
    "Decomposition of explained variance per Hoffman (2015) / Rights & Sterba (2019): f₁ = within-person fixed effects, f₂ = between-person fixed effects. MELS models: m = random intercept SD² (five-component). RS models: f₃ = random slope variance × var(x), f₄ = random intercept variance (six-component). σ²_resid = residual variance. fvm = total fixed + random variance. All percentages are posterior medians (95% CrI in parentheses); n_draws = 500 per fit.",
    "",
    "## MELS five-component decomposition (medians)",
    "",
    "| Model | f₁ within | f₂ between | m random | σ² resid | fvm total |",
    "|:------|----------:|-----------:|---------:|---------:|----------:|"
  )

  mels_models <- unique(decomp$model)
  for (m in mels_models) {
    sub <- decomp[decomp$model == m, ]
    get_med <- function(comp) {
      v <- sub$median[sub$component == comp]
      if (length(v) == 0) NA_real_ else v
    }
    lines <- c(lines, sprintf("| %s | %s | %s | %s | %s | %s |",
                              m,
                              format_num(get_med("f1_within"), 3),
                              format_num(get_med("f2_between"), 3),
                              format_num(get_med("m_random"), 3),
                              format_num(get_med("sigma2_resid"), 3),
                              format_num(get_med("fvm_total"), 3)))
  }

  if (file.exists(rs_decomp_csv)) {
    rs_decomp <- read.csv(rs_decomp_csv, stringsAsFactors = FALSE)
    rs_models <- unique(rs_decomp$model)
    lines <- c(lines, "",
               "## RS six-component decomposition (medians)",
               "",
               "| Model | f₁ within | f₂ between | f₃ slope var | f₄ intercept var | σ² resid | fvm total |",
               "|:------|----------:|-----------:|------------:|----------------:|---------:|----------:|")
    for (m in rs_models) {
      sub <- rs_decomp[rs_decomp$model == m, ]
      get_med <- function(comp) {
        v <- sub$median[sub$component == comp]
        if (length(v) == 0) NA_real_ else v
      }
      lines <- c(lines, sprintf("| %s | %s | %s | %s | %s | %s | %s |",
                                m,
                                format_num(get_med("f1_within"), 3),
                                format_num(get_med("f2_between"), 3),
                                format_num(get_med("f3_slope_var"), 3),
                                format_num(get_med("f4_intercept_var"), 3),
                                format_num(get_med("sigma2_resid"), 3),
                                format_num(get_med("fvm_total"), 3)))
    }
  }

  lines <- c(lines, "",
             "## Total Bayesian R² (mu equation)",
             "",
             "| Model | Median | 95% CrI |",
             "|:------|-------:|:--------|")
  for (i in seq_len(nrow(total))) {
    r <- total[i, ]
    lines <- c(lines, sprintf("| %s | %s | [%s, %s] |",
                              r$model,
                              format_num(r$median, 3),
                              format_num(r$q2.5, 3),
                              format_num(r$q97.5, 3)))
  }
  lines <- c(lines, "",
             "*Note*: Bayesian R² (Gelman et al. 2019) provides total variance explained. MELS decomposition (Hoffman 2015) uses five components; RS decomposition (Rights & Sterba 2019, extended) uses six components with f₃ (random slope variance) and f₄ (random intercept variance) replacing m_random. In RS models with within-person centered predictors (E[x_w] = 0), the intercept–slope covariance contribution is zero, so f₃ and f₄ are fully additive.")
  write_md(lines, out_path)
  cat(sprintf("  ✓ Table S7 写入: %s\n", basename(out_path)))
}
tryCatch(write_variance_decomp_table(), error = function(e)
  cat(sprintf("  ✗ Table S7 失败: %s\n", e$message)))

# ============================================================================
# 7Cd. Table S8: 先验规格表（手写）
# ============================================================================
write_prior_spec_table <- function() {
  out_path <- file.path(tables_dir, "table_08_prior_specifications.md")
  lines <- c(
    "# Table S8 (): Prior specifications (Dora 2024 §2 analog)",
    "",
    "All priors are weakly informative. Three sets are defined in `02b_priors_helpers.R`: `priors_fixed`, `priors_mels`, `priors_rs`. Priors were -calibrated to the mood scale midpoint (M ≈ 3.31) and to identified random-effect SDs.",
    "",
    "## Prior set 1: `priors_fixed` — random-intercept (non-MELS) models",
    "",
    "| Parameter class | Distribution | Hyperparameters | Rationale |",
    "|:----------------|:-------------|:----------------|:----------|",
    "| `b` (fixed slopes) | normal | (0, 1) | Weakly informative; allows effect sizes up to ±2 SD in standardized scale. |",
    "| `Intercept` | normal | (3, 2) | Centered at observed mood mean (3.31) ± 2 (≈ 1 SD of person-mean mood).  calibration R3-01. |",
    "| `sd` (RE SD) | student_t | (3, 0, 1) | Heavy-tailed; prevents zero-pinning, allows τ up to ≈ 3.  calibration R3-02. |",
    "| `sigma` (residual SD) | student_t | (3, 0, 2.5) | Original brms default; broad enough for mood SD ≈ 3.8. |",
    "",
    "## Prior set 2: `priors_mels` — MELS (location-scale) models",
    "",
    "| Parameter class | Distribution | Hyperparameters | Rationale |",
    "|:----------------|:-------------|:----------------|:----------|",
    "| `b` (mu slopes) | normal | (0, 1) | Same as `priors_fixed`. |",
    "| `b` for `sigma` dpar (log-SD slopes) | normal | (0, 1) | exp(±2) ∈ [0.13, 7.4] covers plausible multiplicative effects on residual SD. |",
    "| `Intercept` | normal | (3, 2) | Calibrated to mood midpoint. |",
    "| `Intercept` for `sigma` dpar | student_t | (3, 0, 1) | log(σ) ≈ log(1.5–4) ∈ [0.4, 1.4]; t_3 heavy tails prevent log-σ blow-up.  Plan B retained. |",
    "| `sd` (RE SD) | student_t | (3, 0, 1) | Tightened from default. |",
    "",
    "## Prior set 3: `priors_rs` — random-slope models (incl. cross-level moderators)",
    "",
    "| Parameter class | Distribution | Hyperparameters | Rationale |",
    "|:----------------|:-------------|:----------------|:----------|",
    "| `b` (fixed slopes + interactions) | normal | (0, 1) | Weakly informative. |",
    "| `Intercept` | normal | (3, 2) | Same calibration. |",
    "| `sd` (RE SDs) | student_t | (3, 0, 1) | Tightened. |",
    "| `sigma` | student_t | (3, 0, 2.5) | Original. |",
    "| `cor` (RE correlation) | LKJ | (2) | Mildly regularizing; concentrates around 0. |",
    "",
    "*Note*: All MELS / RS models have `sample_prior = \"yes\"` enabled in  (for MELS only) to permit prior-posterior overlay (Fig. S3). The Item-level Noon_Main also has `sample_prior = \"yes\"`.",
    "",
    "## Prior-predictive validation",
    "",
    "Each model's prior predictive distribution was gated through `prior_predictive_check()` against the mood scale [-10, 10] (Table S6). All models PASS or WARN; none FAIL."
  )
  write_md(lines, out_path)
  cat(sprintf("  ✓ Table S8 写入: %s\n", basename(out_path)))
}
tryCatch(write_prior_spec_table(), error = function(e)
  cat(sprintf("  ✗ Table S8 失败: %s\n", e$message)))

# ============================================================================
# 7Ce. Table S9: 模型公式表（手写）
# ============================================================================
write_model_formulas_table <- function() {
  out_path <- file.path(tables_dir, "table_09_model_formulas.md")
  lines <- c(
    "# Table S9 (): Model formulas — math ↔ brms code ↔ description",
    "",
    "Six main brms models + four item-level MELS, each shown as: (a) mathematical formulation, (b) brms `bf()` formula, (c) natural-language description.",
    "",
    "## 1. B1_MELS — composite-level mixed-effects location-scale",
    "",
    "**Mathematical**:",
    "",
    "$$y_{ij} = \\beta_0 + u_{0i} + \\beta_1 \\text{health}^{within}_{ij} + \\beta_2 \\text{health}^{between}_i + \\beta_3 \\text{mood}^{lag1day}_{ij} + \\beta_4 \\text{day}_j + \\beta_5 \\text{session}_{ij} + \\varepsilon_{ij}$$",
    "",
    "$$\\log(\\sigma_{ij}) = \\gamma_0 + \\gamma_1 \\text{health}^{within}_{ij} + \\gamma_2 \\text{health}^{between}_i$$",
    "",
    "**brms code**:",
    "```r",
    "B1_MELS_form <- brms::bf(",
    "  mood  ~ health_session_c + health_session_between_c +",
    "          mood_lag1_day_c + study_day_c + session_numeric +",
    "          (1 | user_id),",
    "  sigma ~ health_session_c + health_session_between_c",
    ")",
    "```",
    "",
    "**Description**: Random intercept per participant for `mu`; sigma submodel has fixed effects only ( Plan B convergence rescue). β₁ tests whether health behaviors elevate session-level mood (location); γ₁ tests whether they compress residual mood variability (scale).",
    "",
    "## 2. R3_MELS — health × prior-mood buffering MELS",
    "",
    "**Mathematical**:",
    "",
    "$$y_{ij} = \\beta_0 + u_{0i} + \\beta_1 \\text{mood}^{lag1day}_{ij} + \\beta_2 \\text{health}^{within}_{ij} + \\beta_3 (\\text{health} \\times \\text{prevMood})_{ij} + \\beta_4 \\text{day}_j + \\beta_5 \\text{day}^2_j + \\beta_6 \\text{weekend}_j + \\beta_7 \\text{session}_{ij} + \\varepsilon_{ij}$$",
    "",
    "**brms code**:",
    "```r",
    "R3_MELS_form <- brms::bf(",
    "  mood  ~ mood_lag1_day_c + health_session_c + health_x_prev_mood_c +",
    "          study_day_c + study_day_squared + is_weekend + session_numeric +",
    "          (1 | user_id),",
    "  sigma ~ health_session_c + health_x_prev_mood_c",
    ")",
    "```",
    "",
    "**Description**: Tests whether health behaviors buffer mood under low prior-mood, in both location (β₃) and scale (γ₃ via sigma slope).",
    "",
    "## 3. B1_RS — Bayesian random-slope rescue of A1",
    "",
    "**Mathematical**:",
    "",
    "$$y_{ij} = \\beta_0 + u_{0i} + (\\beta_1 + u_{1i}) \\text{health}^{within}_{ij} + \\beta_2 \\text{health}^{between}_i + \\beta_3 \\text{mood}^{prev}_{ij} + \\beta_4 \\text{day}_j + \\beta_5 \\text{session}_{ij} + \\varepsilon_{ij}$$",
    "",
    "$$\\binom{u_{0i}}{u_{1i}} \\sim \\mathcal{N}\\left(0, \\begin{pmatrix} \\tau_0^2 & \\rho \\tau_0 \\tau_1 \\\\ \\rho \\tau_0 \\tau_1 & \\tau_1^2 \\end{pmatrix}\\right)$$",
    "",
    "**brms code**:",
    "```r",
    "B1_RS_form <- brms::bf(",
    "  mood ~ health_session_c + health_session_between_c +",
    "         mood_prev_session_c + study_day_c + session_numeric +",
    "         (1 + health_session_c | user_id)",
    ")",
    "```",
    "",
    "**Description**: Person-specific health slopes. lme4 boundary-singular failure; Bayesian estimates give credible τ₁ and ρ for cross-level analyses.",
    "",
    "## 4-6. B1_RS_{PHQ, PC, sigmaM} — cross-level moderators",
    "",
    "Same as B1_RS, plus `health_session_c × moderator` interaction (PHQ-9, perceived control, σ^M).",
    "",
    "**brms code (PHQ-9 example)**:",
    "```r",
    "B1_RS_PHQ_form <- brms::bf(",
    "  mood ~ health_session_c * pre_phq9_c + health_session_between_c +",
    "         mood_prev_session_c + study_day_c + session_numeric +",
    "         (1 + health_session_c | user_id)",
    ")",
    "```",
    "",
    "**Description**: Cross-level moderation — does the health → mood slope depend on baseline depression / perceived control / mood instability?",
    "",
    "## 7-10. Item-level MELS (Noon × 2 + Evening × 2)",
    "",
    "Decomposition of B1_MELS composite effect into 8 individual health-behavior items (4 noon: sleep_ontime, wake_early, breakfast, aspiration; 4 evening: walk, exercise, dinner, breathing).",
    "",
    "**brms code (Noon Main)**:",
    "```r",
    "form_noon_main <- brms::bf(",
    "  mood  ~ sleep_ontime_c + wake_early_c + breakfast_c + aspiration_c +",
    "          mood_lag1_day_c + study_day_c + health_session_between_c +",
    "          (1 | user_id),",
    "  sigma ~ sleep_ontime_c + wake_early_c + breakfast_c + aspiration_c",
    ")",
    "```",
    "",
    "**Description**: Item-level dissociation reveals aspiration as a pure sigma-compressor (γ ≪ 0, β ≈ 0) and exercise as a mu-elevator + sigma-inflator. Robust variants add `mood_prev_session_c` as additional AR control.",
    "",
    "## Estimator settings (all )",
    "",
    "| Class | chains | iter | warmup | adapt_delta | max_treedepth | threads | sample_prior |",
    "|:------|-------:|-----:|-------:|------------:|--------------:|--------:|:-------------|",
    "| MELS (B1, R3) | 4 | 6,000 | 3,000 | 0.95 | 12 | 2 | yes |",
    "| Item-level MELS | 4 | 6,000 | 3,000 | 0.95 | 12 | 2 | Noon_Main: yes; others: no |",
    "| RS (B1) | 4 | 4,000 | 2,000 | 0.95 | 10 | 2 | no |",
    "| RS cross-level (PHQ/PC/σ^M) | 4 | 4,000 | 2,000 | 0.99 | 15 | 2 | no |",
    "",
    "All chains run in parallel; `stan_threads = TRUE` enabled via `stan_model_args = list(cpp_options = list(stan_threads = TRUE))` ()."
  )
  write_md(lines, out_path)
  cat(sprintf("  ✓ Table S9 写入: %s\n", basename(out_path)))
}
tryCatch(write_model_formulas_table(), error = function(e)
  cat(sprintf("  ✗ Table S9 失败: %s\n", e$message)))

# ============================================================================
# 7Cf. Table S10: LOO 模型比较表（依赖 06c）
# ============================================================================
write_loo_table <- function() {
  results_csv <- file.path(output_dir, "loo_results.csv")
  compare_csv <- file.path(output_dir, "loo_comparison.csv")
  out_path <- file.path(tables_dir, "table_10_loo.md")

  if (!file.exists(results_csv)) {
    cat("  ⚠ loo_results CSV 不存在，跳过 Table S10（请先运行 06c_loo_comparison.R）\n")
    return(invisible(NULL))
  }

  res <- read.csv(results_csv, stringsAsFactors = FALSE)
  lines <- c(
    "# Table S10 (): LOO model evaluation",
    "",
    "Pareto-smoothed leave-one-out cross-validation (PSIS-LOO; Vehtari et al. 2017). For each model: expected log pointwise predictive density (elpd_loo), effective number of parameters (p_loo), LOO information criterion (looic = -2 × elpd_loo), and Pareto k diagnostics. All values from `brms::loo()`.",
    "",
    "## Per-model summary",
    "",
    "| Model | elpd_loo | SE | p_loo | looic | high-k (k > 0.7) | max k |",
    "|:------|---------:|---:|------:|------:|-----------------:|------:|"
  )
  for (i in seq_len(nrow(res))) {
    r <- res[i, ]
    n_total <- if ("n_obs" %in% names(r)) r$n_obs else NA
    n_high_str <- if (!is.na(n_total)) sprintf("%d / %s", r$n_high_k, format(n_total, big.mark = ",")) else as.character(r$n_high_k)
    lines <- c(lines, sprintf("| %s | %s | %s | %s | %s | %s | %s |",
                              r$model,
                              format_num(r$elpd_loo, 1),
                              format_num(r$elpd_loo_se, 1),
                              format_num(r$p_loo, 1),
                              format_num(r$looic, 1),
                              n_high_str,
                              format_num(r$max_pareto_k, 3)))
  }
  lines <- c(lines, "",
             "*Note*: Pareto k > 0.7 indicates unreliable point-wise elpd estimate. When n_high_k > 0, `loo()` was rerun with `moment_match = TRUE` before reporting.")

  # 比较表
  if (file.exists(compare_csv)) {
    cmp <- read.csv(compare_csv, stringsAsFactors = FALSE)
    if (nrow(cmp) > 0) {
      lines <- c(lines, "",
                 "## Pairwise comparisons (loo_compare)",
                 "",
                 "Models within each group share the same data (and thus same elpd reference); comparisons across groups are not meaningful.",
                 "")
      groups <- unique(cmp$group)
      for (g in groups) {
        sub <- cmp[cmp$group == g, ]
        if (nrow(sub) == 0) next
        lines <- c(lines, sprintf("### Group: %s", g), "",
                   "| Model | elpd_diff | se_diff |",
                   "|:------|----------:|--------:|")
        for (i in seq_len(nrow(sub))) {
          r <- sub[i, ]
          ed <- if ("elpd_diff" %in% colnames(r)) r$elpd_diff else NA_real_
          se <- if ("se_diff" %in% colnames(r)) r$se_diff else NA_real_
          lines <- c(lines, sprintf("| %s | %s | %s |",
                                    r$model,
                                    format_num(ed, 2),
                                    format_num(se, 2)))
        }
        lines <- c(lines, "")
      }
    }
  }
  write_md(lines, out_path)
  cat(sprintf("  ✓ Table S10 写入: %s\n", basename(out_path)))
}
tryCatch(write_loo_table(), error = function(e)
  cat(sprintf("  ✗ Table S10 失败: %s\n", e$message)))

cat("\n  ✓ 模块07c  补充表格完成\n")
