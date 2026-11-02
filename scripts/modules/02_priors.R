# ============================================================================
# 模块 02b (): 先验定义 + brms 辅助函数
# ============================================================================
# 依赖: 00b_setup_brms.R, 01b_data_import.R
# 产出: priors_fixed / priors_mels / priors_rs（三套先验对象）
#       fit_or_load() 断点续跑包装器（支持多框架审计）
#       extract_posterior_summary() 后验摘要
#       bayes_r2_decomp() 五分量分解
#       compare_to_lme4() 跨范式对照
# ============================================================================
# 设计决策:
#   - 弱信息先验：fixed normal(0,1); SDs student_t(3,0,2.5); cor LKJ(2)
#   - MELS sigma 截距用更紧的 student_t(3,0,1) 防止 log-scale 后验飘出
#   - .rds 用 xz 压缩；force_refit 通过 FORCE_REFIT 全局向量控制
#   - : 多框架 CONFIG + 自动审计追踪
# ============================================================================

cat("\n--- 模块02b (): 先验定义 + brms 辅助函数 ---\n")

# ============================================================================
# 2B0. 多框架配置注册表（ 新增）
# ============================================================================
# 修改此处的任何值都必须在提案中明确说明
# 支持多框架：brms（贝叶斯）/ lme4（频率学派）
# ============================================================================
CONFIG <- list(
  brms = list(
    iter            = 6000L,
    warmup          = 3000L,
    chains          = 4L,
    adapt_delta     = 0.95,
    max_treedepth   = 15L,
    prior_intercept = "normal(3, 2)",
    prior_coef      = "normal(0, 1)",
    prior_sd        = "student_t(3, 0, 1)",
    seed            = 42L
  ),
  lme4 = list(
    reml      = TRUE,
    optimizer = "bobyqa",
    control   = list(maxit = 10000L, optCtrl = list()),
    seed      = 42L
  ),
  config_version = "v10.1.0",
  last_modified  = "2026-05-31"
)

cat(sprintf("  ✓ 多框架 CONFIG 已加载 (版本: %s)\n", CONFIG$config_version))

# ============================================================================
# 2Ba. 三套先验对象
# ============================================================================
#  (2026-05-19) 先验升级（建议文档 A+B+C）:
#   - Intercept: student_t(3,0,2.5) → normal(3,2)  （A: 校准到 mood 量表中点 M≈3.31）
#   - sd:        student_t(3,0,2.5) → student_t(3,0,1)  （B: 折中收紧，保留重尾）
#   - σ_Intercept: student_t(3,0,1) 维持  （C: Plan B 收敛保护不变）
# 不采纳: sd → normal(0,0.5) 过激；σ_Intercept → normal(0,1.5) 失去保护
# 先验规格变更记录: results/logs/config_audit_log.jsonl
#
# 适用于：所有非 MELS 的随机截距 brms 模型（无 RS、无 sigma 方程）
priors_fixed <- c(
  brms::prior(normal(0, 1),         class = "b"),
  brms::prior(normal(3, 2),         class = "Intercept"),   # A
  brms::prior(student_t(3, 0, 1),   class = "sd"),          # B
  brms::prior(student_t(3, 0, 2.5), class = "sigma")
)

# 适用于：MELS（有 sigma 子方程）模型 — sigma 截距用更紧先验防溢出
# sigma 子方程的 sd（随机截距 SD）单独收紧到 student_t(3,0,1)
# 方案 A: 移除 |p| 管道, 删除 lkj(2) cor 先验
# 方案 B: sigma 子方程已无 (1|user_id),
#   故删除 class="sd" dpar="sigma" 先验(否则 brms 报 "no parameter")
priors_mels <- c(
  brms::prior(normal(0, 1),         class = "b"),
  brms::prior(normal(0, 1),         class = "b",         dpar = "sigma"),
  brms::prior(normal(3, 2),         class = "Intercept"),                # A
  brms::prior(student_t(3, 0, 1),   class = "Intercept", dpar = "sigma"),# C 维持
  brms::prior(student_t(3, 0, 1),   class = "sd")                        # B
)

# 适用于：贝叶斯随机斜率（含 1+slope|user_id 与跨层交互）
priors_rs <- c(
  brms::prior(normal(0, 1),         class = "b"),
  brms::prior(normal(3, 2),         class = "Intercept"),  # A
  brms::prior(student_t(3, 0, 1),   class = "sd"),         # B
  brms::prior(student_t(3, 0, 2.5), class = "sigma"),
  brms::prior(lkj(2),               class = "cor")
)

cat("  ✓ 三套先验已定义 (): priors_fixed / priors_mels / priors_rs\n")

# ============================================================================
# 2Bb. fit_or_load() — 断点续跑包装器（: 多框架 + 自动审计）
# ============================================================================
# 行为:
#   1. 检查 fits/<model_name>.rds 是否存在
#   2. 存在且 model_name 不在 FORCE_REFIT 中 → 直接 readRDS()
#   3. 否则根据 framework 参数调用对应拟合函数
#   4. 自动记录审计日志到 results/logs/config_audit_log.jsonl
#   5. 自动检测同一模型历史参数变更，变更时控制台打印警告
# 用法:
#   fit <- fit_or_load("B1_MELS", form, data = sample_B,
#                      prior = priors_mels, framework = "brms", ...)
#   fit <- fit_or_load("B1_MLM", form, data = sample_B,
#                      framework = "lme4", reml = TRUE, ...)
fit_or_load <- function(model_name,
                        formula,
                        data,
                        prior = NULL,
                        framework = "brms",
                        chains = 4L,
                        iter = 2500L,
                        warmup = 1000L,
                        adapt_delta = 0.95,
                        max_treedepth = 12L,
                        threads = NULL,
                        stan_threads = TRUE,
                        family = brms::brmsfamily("gaussian"),
                        seed = GLOBAL_SEED,
                        sample_prior = "no",
                        reml = TRUE,
                        optimizer = "bobyqa",
                        control_maxit = 10000L,
                        ...) {
  rds_path <- file.path(fits_dir, paste0(model_name, ".rds"))
  log_path <- file.path(logs_dir, paste0(model_name, "_cmdstan.log"))

  refit_needed <- (!file.exists(rds_path)) || (model_name %in% FORCE_REFIT)

  if (!refit_needed) {
    cat(sprintf("  [%s] 加载已有 .rds (跳过拟合)\n", model_name))
    fit <- readRDS(rds_path)
    return(fit)
  }

  if (!framework %in% c("brms", "lme4")) {
    stop(sprintf("未知框架: %s。支持的框架: brms, lme4", framework))
  }

  default_params <- CONFIG[[framework]]

  if (framework == "brms") {
    actual_params <- list(
      iter          = if (missing(iter)) default_params$iter else iter,
      warmup        = if (missing(warmup)) default_params$warmup else warmup,
      chains        = chains,
      adapt_delta   = adapt_delta,
      max_treedepth = max_treedepth,
      seed          = seed
    )
  } else {
    actual_params <- list(
      reml          = reml,
      optimizer     = optimizer,
      control       = list(maxit = control_maxit, optCtrl = list()),
      seed          = seed
    )
  }

  audit_entry <- list(
    timestamp = format(Sys.time(), "%Y%m%d_%H%M%S"),
    framework = framework,
    model     = model_name,
    params    = actual_params
  )

  audit_file <- file.path(output_dir, "logs", "config_audit_log.jsonl")
  dir.create(dirname(audit_file), showWarnings = FALSE, recursive = TRUE)

  if (file.exists(audit_file)) {
    all_lines <- readLines(audit_file, warn = FALSE)
    if (length(all_lines) > 0) {
      all_records <- lapply(all_lines, function(l) {
        tryCatch(jsonlite::fromJSON(l), error = function(e) NULL)
      })
      all_records <- all_records[!sapply(all_records, is.null)]
      same_model <- all_records[sapply(all_records, function(r) 
        r$framework == framework && r$model == model_name)]
      
      if (length(same_model) >= 1) {
        last_record <- same_model[[length(same_model)]]
        diff_params <- list()
        param_names <- union(names(actual_params), names(last_record$params))
        for (p in param_names) {
          v1 <- actual_params[[p]]
          v2 <- last_record$params[[p]]
          if (!identical(v1, v2)) {
            diff_params[[p]] <- list(old = v2, new = v1)
          }
        }
        
        if (length(diff_params) > 0) {
          warning_msg <- sprintf(
            "\n╔══════════════════════════════════════════════════════════════╗\n║  ⚠️  %s/%s 参数变更 detected                              ║\n╚══════════════════════════════════════════════════════════════╝\n  上次运行: %s\n  变更项: %s\n  变更详情见 results/logs/config_audit_log.jsonl\n",
            framework, model_name,
            last_record$timestamp,
            paste(names(diff_params), collapse = ", "),
            model_name
          )
          warning(warning_msg)
        }
      }
    }
  }

  tryCatch({
    json_line <- jsonlite::toJSON(audit_entry, auto_unbox = TRUE)
    cat(json_line, "\n", file = audit_file, append = file.exists(audit_file))
    cat(sprintf("  审计记录: %s/%s @ %s\n", framework, model_name, audit_entry$timestamp))
  }, error = function(e) {
    warning(sprintf("审计日志写入失败: %s", e$message))
  })

  if (framework == "brms") {
    thread_info <- if (!is.null(threads)) sprintf(", threads=%d", threads) else ""
    cat(sprintf("  [%s] 拟合中 ... (chains=%d, iter=%d, warmup=%d, adapt_delta=%.2f%s)\n",
                model_name, actual_params$chains, actual_params$iter, actual_params$warmup, 
                actual_params$adapt_delta, thread_info))
    t0 <- Sys.time()

    con <- file(log_path, open = "wt")
    sink(con, type = "message")
    on.exit({
      sink(NULL, type = "message")
      close(con)
    }, add = TRUE)

    fit <- brms::brm(
      formula = formula,
      data = data,
      family = family,
      prior = prior,
      chains = actual_params$chains,
      iter = actual_params$iter,
      warmup = actual_params$warmup,
      cores = min(actual_params$chains, getOption("mc.cores")),
      seed = actual_params$seed,
      control = list(adapt_delta = actual_params$adapt_delta,
                     max_treedepth = actual_params$max_treedepth),
      backend = "cmdstanr",
      save_pars = brms::save_pars(all = TRUE),
      refresh = 200,
      threads = if (!is.null(threads)) brms::threading(threads) else NULL,
      stan_model_args = if (isTRUE(stan_threads)) list(cpp_options = list(stan_threads = TRUE)) else list(),
      sample_prior = sample_prior,
      ...
    )

    saveRDS(fit, rds_path, compress = "xz")
    elapsed <- difftime(Sys.time(), t0, units = "mins")
    cat(sprintf("  [%s] ✓ 拟合完成，耗时 %.2f 分钟，已保存 %s\n",
                model_name, as.numeric(elapsed), basename(rds_path)))
  } else if (framework == "lme4") {
    cat(sprintf("  [%s] 拟合中 ... (reml=%s, optimizer=%s)\n",
                model_name, actual_params$reml, actual_params$optimizer))
    t0 <- Sys.time()

    fit <- lme4::lmer(
      formula = formula,
      data = data,
      REML = actual_params$reml,
      control = lme4::lmerControl(
        optimizer = actual_params$optimizer,
        optCtrl = actual_params$control$optCtrl
      ),
      ...
    )

    saveRDS(fit, rds_path, compress = "xz")
    elapsed <- difftime(Sys.time(), t0, units = "mins")
    cat(sprintf("  [%s] ✓ 拟合完成，耗时 %.2f 分钟，已保存 %s\n",
                model_name, as.numeric(elapsed), basename(rds_path)))
  }

  fit
}

# ============================================================================
# 2Bc. extract_posterior_summary() — 后验摘要表
# ============================================================================
# 返回: data.frame，每行一个参数，列含 mean/median/sd/q2.5/q97.5/rhat/ess_bulk/ess_tail
# 用法: write.csv(extract_posterior_summary(fit), "...csv")
extract_posterior_summary <- function(fit) {
  draws <- posterior::as_draws_df(fit)
  s <- posterior::summarise_draws(
    draws,
    mean, median, sd,
    ~ quantile(.x, probs = c(0.025, 0.975), na.rm = TRUE),
    posterior::default_convergence_measures()
  )
  # summarise_draws 默认列名 2.5%, 97.5% 不便引用 → 改名
  names(s)[names(s) == "2.5%"]  <- "q2.5"
  names(s)[names(s) == "97.5%"] <- "q97.5"
  s
}

# ============================================================================
# 2Bc2. prior_predictive_check() —  新增（建议文档 D 项）
# ============================================================================
# 在拟合前对先验做预测检查，验证先验生成的 y 大致落在 mood 量表 [-10, 10] 内。
# 用作 05b 的 gating 步骤； 设三级判定:
#   - PASS  (pct_in_range ≥ 0.70): 先验合理，无须调整
#   - WARN  (0.10 ≤ pct_in_range < 0.70): 先验偏宽或 predictor 尺度大；不阻断流程
#                                          但记录到 CSV 供 + 决策
#   - FAIL  (pct_in_range < 0.10): catastrophic — 先验明显病态，停止流程
# 阈值选择理由：
#   - 0.70 适用于固定 sigma 模型（RS 系列）
#   - MELS log-link sigma 子方程下 student_t(3,0,1) σ_Intercept 给 sigma ∈ [0.08, 12]
#     使 y_prior 自然超出 [-10, 10]；约 0.76 是该模型族的实证下限
#   - R3_MELS 含 study_day_squared (range ~0-100) 与 normal(0,1) β 相乘，
#     先验 linear predictor 范围极大；约 0.29 是该模型的实证下限
#   - 0.10 catastrophic 阈值保证全失控的先验仍被截获
prior_predictive_check <- function(formula, data, prior,
                                   n_draws = 500L,
                                   y_range = c(-10, 10),
                                   pct_threshold_warn = 0.70,
                                   pct_threshold_fail = 0.10,
                                   chains = 2L,
                                   iter = 1000L,
                                   warmup = 500L,
                                   seed = GLOBAL_SEED) {
  fit_prior <- brms::brm(
    formula = formula,
    data = data,
    family = brms::brmsfamily("gaussian"),
    prior = prior,
    sample_prior = "only",
    chains = chains,
    iter = iter,
    warmup = warmup,
    cores = min(chains, getOption("mc.cores")),
    seed = seed,
    backend = "cmdstanr",
    refresh = 0,
    save_pars = brms::save_pars(all = TRUE)
  )
  y_prior <- brms::posterior_predict(fit_prior, ndraws = n_draws)
  pct_in_range <- mean(y_prior >= y_range[1] & y_prior <= y_range[2],
                       na.rm = TRUE)
  status <- if (pct_in_range >= pct_threshold_warn) {
    "PASS"
  } else if (pct_in_range >= pct_threshold_fail) {
    "WARN"
  } else {
    "FAIL"
  }
  list(
    pct_in_range  = pct_in_range,
    median_mean   = median(rowMeans(y_prior, na.rm = TRUE)),
    median_sd     = median(apply(y_prior, 1, sd, na.rm = TRUE)),
    pass          = status != "FAIL",
    status        = status,
    fit_prior     = fit_prior
  )
}

# ============================================================================
# 2Bd. bayes_r2_decomp() — 五分量方差分解（与 r2mlm 平行）
# ============================================================================
# 输入: brmsfit, within_vars/between_vars 字符向量（以 sample 中变量名表示）
# 输出: data.frame with cols [component, median, q2.5, q97.5]
#   component ∈ {f1_within, f2_between, m_random, sigma2_resid, fvm_total}
# 算法:
#   逐 draw 计算：
#     mu_w = X_within %*% beta_within_draw
#     mu_b = Z_between %*% gamma_between_draw
#     var(mu_w), var(mu_b), tau_intercept^2, sigma^2
#   占比 = 各分量 / 总和
#   原 MELS 分支用 posterior_predict 在循环内重抽样，且 var(y - yhat)
#   包含随机截距 + 残差噪声 + sigma 异质性的总和，与其他四分量重复计入。
#   修正：循环外一次性 posterior_epred(dpar="sigma") 得 sigma_per_row,
#         var_resid_draw = mean(sigma_i^2) over rows，与"在 MELS 框架下
#         的期望残差方差"概念一致，且不再双重计入 random intercept。
#         附加收益：MELS 单模型分解从 ~10 分钟降到 ~30 秒。
bayes_r2_decomp <- function(fit, data, within_vars, between_vars = character(0),
                            n_draws = 1000L) {
  # 抽样 n_draws 个后验样本（默认 1000，全部 8000 太慢）
  total_draws <- brms::ndraws(fit)
  draw_ids <- sample.int(total_draws, size = min(n_draws, total_draws))

  # 提取所有 draw 的 fixed effect coefficient 矩阵
  fix_draws <- as.matrix(brms::as_draws_matrix(fit, variable = "^b_",
                                                regex = TRUE))[draw_ids, , drop = FALSE]

  # 提取每 draw 的 sigma （residual SD）
  is_mels <- !("sigma" %in% brms::variables(fit))
  sigma_draws <- if (!is_mels) {
    as.numeric(brms::as_draws_matrix(fit, variable = "sigma")[draw_ids])
  } else {
    rep(NA_real_, length(draw_ids))
  }

  # MELS 分支：循环外一次性预测 per-row sigma，避免 var(y - yhat) 双重计入
  mels_var_resid_draws <- if (is_mels) {
    # dpar="sigma" 返回 ndraws × N 矩阵；re_formula = NA 排除随机效应
    # （beta_v2.2 起 sigma 子方程已无 (1|user_id)，此参数为安全显式声明）
    sigma_per_row <- brms::posterior_epred(
      fit, dpar = "sigma",
      ndraws = length(draw_ids),
      draw_ids = draw_ids,
      re_formula = NA
    )
    rowMeans(sigma_per_row^2, na.rm = TRUE)
  } else {
    rep(NA_real_, length(draw_ids))
  }

  # 提取 random intercept SD（user_id 截距）
  vc <- brms::VarCorr(fit, summary = FALSE)
  tau_int_draws <- if (!is.null(vc$user_id)) {
    as.numeric(vc$user_id$sd[draw_ids, "Intercept"])
  } else {
    rep(0, length(draw_ids))
  }

  # 构造 within / between 设计矩阵
  # 用 model.matrix() 处理交互项（如 "health_session_c:pre_phq9_c"）
  # 这些列在 data 中不作为字面量存在，必须由 model.matrix 自动生成。
  fnames <- colnames(fix_draws)
  fname_short <- sub("^b_", "", fnames)

  # 构造完整设计矩阵（含交互项）
  # 用正则剥离随机效应项（lme4::nobars 不识别 brms 风格的
  # (1 + slope | user_id) 语法）。这样可在 stats::model.matrix 中得到 fixed
  # design matrix（含交互项展开）。
  # 注意: stats::formula(brmsfit) 返回 brmsformula 对象（class "brmsformula"），
  # 需先取 $formula slot 才是 R formula。
  strip_re_terms <- function(form) {
    s <- paste(deparse(form), collapse = " ")
    # 反复移除 ( ... | ... ) 直到无随机效应残留
    while (grepl("\\([^()]*\\|[^()]*\\)", s)) {
      s <- gsub("\\s*\\+?\\s*\\([^()]*\\|[^()]*\\)\\s*", "", s)
    }
    # 清理可能残留的多余加号
    s <- gsub("\\+\\s*$", "", s)
    s <- gsub("~\\s*\\+", "~", s)
    stats::as.formula(s)
  }

  bf_obj <- stats::formula(fit)
  # brms 返回的可能是 brmsformula 对象或普通 formula；统一访问 $formula
  raw_form <- if (inherits(bf_obj, "brmsformula") && !is.null(bf_obj$formula)) {
    bf_obj$formula
  } else {
    bf_obj
  }
  fixed_form <- tryCatch(
    strip_re_terms(raw_form),
    error = function(e) raw_form
  )
  mm <- tryCatch(
    stats::model.matrix(stats::delete.response(stats::terms(fixed_form)),
                        data = data),
    error = function(e) NULL
  )
  mm_colnames <- if (!is.null(mm)) colnames(mm) else character(0)

  # 模糊匹配：var_name 在 mm 列名中可能以 "var1:var2" 或 "var2:var1" 出现
  match_design_col <- function(var_name, available) {
    if (var_name %in% available) return(var_name)
    # 冒号形式：尝试反序
    if (grepl(":", var_name, fixed = TRUE)) {
      parts <- strsplit(var_name, ":", fixed = TRUE)[[1]]
      reversed <- paste(rev(parts), collapse = ":")
      if (reversed %in% available) return(reversed)
    }
    NA_character_
  }

  # 注意：不能按"位置"配对系数与设计矩阵列 ——
  #   beta 列序来自 which(fname_short %in% vars)（= brms 内部存储顺序），
  #   X   列序来自 vapply(vars, ...)（= vars 传入顺序），两者不一致时
  #   系数被乘到错误的预测列上（实测 B1_RS_PC 的 var_f2 由 1.66 膨胀到 181）。
  # 改为"按名对齐"：对每个 var 同时解析它的系数列与设计列，二者以相同顺序入栈，
  #   保证 beta[k] 永远乘以同一预测因子的设计列。
  build_aligned <- function(vars) {
    if (length(vars) == 0) return(NULL)
    beta_idx <- integer(0)
    Xcols    <- list()
    for (v in vars) {
      coef_match <- match_design_col(v, fname_short)   # 系数侧（可反序匹配交互项）
      if (is.na(coef_match)) next                       # 无对应后验系数 → 跳过
      xcol <- NULL
      if (!is.null(mm)) {
        dm <- match_design_col(v, mm_colnames)          # 设计矩阵侧
        if (!is.na(dm)) xcol <- mm[, dm]
      }
      if (is.null(xcol) && v %in% colnames(data)) xcol <- data[[v]]
      if (is.null(xcol)) next                           # 设计列无法解析 → 跳过
      beta_idx <- c(beta_idx, which(fname_short == coef_match)[1])
      Xcols[[length(Xcols) + 1]] <- as.numeric(xcol)
    }
    if (length(beta_idx) == 0) return(NULL)
    X <- do.call(cbind, Xcols)
    X[is.na(X)] <- 0
    list(beta_idx = beta_idx, X = X)   # beta_idx[k] 与 X[, k] 同序，逐 var 配对
  }

  aln_w <- build_aligned(within_vars)
  aln_b <- build_aligned(between_vars)

  result_df <- data.frame()

  for (s in seq_along(draw_ids)) {
    # within / between：按名对齐后，beta 行与 X 列严格同序
    var_f1 <- if (!is.null(aln_w)) {
      beta_w <- fix_draws[s, aln_w$beta_idx, drop = FALSE]
      as.numeric(var(as.numeric(aln_w$X %*% t(beta_w))))
    } else 0

    var_f2 <- if (!is.null(aln_b)) {
      beta_b <- fix_draws[s, aln_b$beta_idx, drop = FALSE]
      as.numeric(var(as.numeric(aln_b$X %*% t(beta_b))))
    } else 0

    var_m <- tau_int_draws[s]^2
    var_resid <- if (is_mels) {
      mels_var_resid_draws[s]
    } else {
      sigma_draws[s]^2
    }

    var_total <- var_f1 + var_f2 + var_m + var_resid

    result_df <- rbind(result_df, data.frame(
      draw = s,
      f1_within    = var_f1 / var_total,
      f2_between   = var_f2 / var_total,
      m_random     = var_m  / var_total,
      sigma2_resid = var_resid / var_total,
      fvm_total    = (var_f1 + var_f2 + var_m) / var_total
    ))
  }

  # 汇总成四分位
  comp_cols <- c("f1_within", "f2_between", "m_random",
                 "sigma2_resid", "fvm_total")
  summary_df <- do.call(rbind, lapply(comp_cols, function(cc) {
    data.frame(
      component = cc,
      median = median(result_df[[cc]]),
      q2.5   = quantile(result_df[[cc]], 0.025),
      q97.5  = quantile(result_df[[cc]], 0.975)
    )
  }))
  rownames(summary_df) <- NULL
  summary_df
}

# ============================================================================
# 2Be. compare_to_lme4() — brms vs lme4 跨范式对照
# ============================================================================
# 输入:
#   brms_fit  — brmsfit 对象
#   lme4_fit  — lmerMod 对象（ 已存的）
#   terms     — 字符向量，要对照的 fixed effect 名（lme4 与 brms 命名应一致）
#   include_p — 是否在输出中包含 lme4_p 列（ 默认 FALSE：项 δ，
#               brms 端无对应概念，列易误导读者）
# 输出: data.frame，每行一个 term，含 lme4 β/SE [+ p]、brms median/sd/q2.5/q97.5、Δ、pp_sign
compare_to_lme4 <- function(brms_fit, lme4_fit, terms, include_p = FALSE) {
  # lme4 端
  lme4_coef <- summary(lme4_fit)$coefficients
  # brms 端
  brms_summary <- extract_posterior_summary(brms_fit)
  brms_b_short <- sub("^b_", "", brms_summary$variable)

  out <- data.frame()
  for (tm in terms) {
    if (!(tm %in% rownames(lme4_coef))) {
      cat(sprintf("    ⚠ lme4 中无 term '%s'，跳过\n", tm))
      next
    }
    if (!(tm %in% brms_b_short)) {
      cat(sprintf("    ⚠ brms 中无 term '%s'，跳过\n", tm))
      next
    }
    lme4_row <- lme4_coef[tm, ]
    brms_idx <- which(brms_b_short == tm)
    brms_row <- brms_summary[brms_idx, ]

    # p-value 列名在 lmerTest::summary 下叫 "Pr(>|t|)"（仅在 include_p=TRUE 时使用）
    p_val <- if (include_p) {
      pcol <- if ("Pr(>|t|)" %in% colnames(lme4_coef)) "Pr(>|t|)" else NA
      if (!is.na(pcol)) lme4_row[pcol] else NA_real_
    } else {
      NULL
    }

    # 后验同号概率
    draws_b <- as.numeric(brms::as_draws_matrix(brms_fit,
                                                variable = paste0("b_", tm)))
    pp_same_sign <- if (lme4_row["Estimate"] >= 0) {
      mean(draws_b > 0)
    } else {
      mean(draws_b < 0)
    }

    delta <- as.numeric(brms_row$median) - as.numeric(lme4_row["Estimate"])
    delta_over_se <- delta / as.numeric(lme4_row["Std. Error"])

    row_data <- data.frame(
      term = tm,
      lme4_beta   = as.numeric(lme4_row["Estimate"]),
      lme4_se     = as.numeric(lme4_row["Std. Error"]),
      brms_median = as.numeric(brms_row$median),
      brms_sd     = as.numeric(brms_row$sd),
      brms_q2.5   = as.numeric(brms_row$q2.5),
      brms_q97.5  = as.numeric(brms_row$q97.5),
      delta       = delta,
      delta_over_se = delta_over_se,
      pp_same_sign  = pp_same_sign,
      brms_rhat   = as.numeric(brms_row$rhat),
      brms_ess_bulk = as.numeric(brms_row$ess_bulk)
    )
    if (include_p) {
      row_data$lme4_p <- as.numeric(p_val)
    }
    out <- rbind(out, row_data)
  }
  out
}

cat("  ✓ 5 个辅助函数已定义: fit_or_load / extract_posterior_summary / prior_predictive_check / bayes_r2_decomp / compare_to_lme4\n")

cat("\n  ✓ 模块02b 先验与辅助函数加载完成 ()\n")
