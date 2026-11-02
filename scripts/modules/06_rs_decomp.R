# ============================================================================
# 模块 06d (): RS 模型六分量方差分解
# ============================================================================
# 依赖: 03b/04b/05b（brms_fits + convergence 已就位）
# 产出: bayes_r2_decomposition_rs.csv
# 方法: Rights & Sterba (2019) 框架，适配 brms 随机斜率模型
# 关键: health_session_c 为 within-person centered (E[x_w] = 0)
#       随机截距-斜率协方差贡献 = 2 * rho_01 * tau_0 * tau_1 * E[x_w] = 0
#       截距与斜率方差完全可加
# ============================================================================
# 六分量定义:
#   f1_within      = var(X_within * beta_within) / var_total
#   f2_between     = var(X_between * beta_between) / var_total
#   f3_slope_var   = tau_1^2 * var(health_session_c) / var_total
#   f4_intercept_var = tau_0^2 / var_total
#   sigma2_resid   = sigma^2 / var_total
#   fvm_total      = 1 - sigma2_resid
# ============================================================================

cat("\n--- 模块06d (): RS 模型六分量方差分解 ---\n")

rs_models <- c("B1_RS", "B1_RS_PHQ", "B1_RS_PC", "B1_RS_sigmaM")

rs_var_map <- list(
  B1_RS = list(
    within  = c("health_session_c", "mood_prev_session_c",
                "study_day_c", "session_numeric"),
    between = c("health_session_between_c")
  ),
  B1_RS_PHQ = list(
    within  = c("health_session_c", "mood_prev_session_c",
                "study_day_c", "session_numeric",
                "health_session_c:pre_phq9_c"),
    between = c("health_session_between_c", "pre_phq9_c")
  ),
  B1_RS_PC = list(
    within  = c("health_session_c", "mood_prev_session_c",
                "study_day_c", "session_numeric",
                "health_session_c:pre_control_c"),
    between = c("health_session_between_c", "pre_control_c")
  ),
  B1_RS_sigmaM = list(
    within  = c("health_session_c", "mood_prev_session_c",
                "study_day_c", "session_numeric",
                "health_session_c:person_std_mood_c"),
    between = c("health_session_between_c", "person_std_mood_c")
  )
)

rs_r2_decomp_6comp <- function(fit, data, within_vars, between_vars,
                                slope_var = "health_session_c", n_draws = 500L) {
  total_draws <- brms::ndraws(fit)
  draw_ids <- sample.int(total_draws, size = min(n_draws, total_draws))

  fix_draws <- as.matrix(brms::as_draws_matrix(fit, variable = "^b_",
                                                regex = TRUE))[draw_ids, , drop = FALSE]

  vc <- brms::VarCorr(fit, summary = FALSE)
  tau_int_draws <- as.numeric(vc$user_id$sd[draw_ids, "Intercept"])
  tau_slope_draws <- as.numeric(vc$user_id$sd[draw_ids, slope_var])
  var_slope_pred <- var(data[[slope_var]], na.rm = TRUE)

  var_f3 <- tau_slope_draws^2 * var_slope_pred
  var_f4 <- tau_int_draws^2

  sigma_draws <- as.numeric(brms::as_draws_matrix(fit, variable = "sigma")[draw_ids])
  var_resid <- sigma_draws^2

  fnames <- colnames(fix_draws)
  fname_short <- sub("^b_", "", fnames)

  strip_re_terms <- function(form) {
    s <- paste(deparse(form), collapse = " ")
    while (grepl("\\([^()]*\\|[^()]*\\)", s)) {
      s <- gsub("\\s*\\+?\\s*\\([^()]*\\|[^()]*\\)\\s*", "", s)
    }
    s <- gsub("\\+\\s*$", "", s)
    s <- gsub("~\\s*\\+", "~", s)
    stats::as.formula(s)
  }

  bf_obj <- stats::formula(fit)
  raw_form <- if (inherits(bf_obj, "brmsformula") && !is.null(bf_obj$formula)) {
    bf_obj$formula
  } else {
    bf_obj
  }
  fixed_form <- tryCatch(strip_re_terms(raw_form), error = function(e) raw_form)
  mm <- tryCatch(
    stats::model.matrix(stats::delete.response(stats::terms(fixed_form)), data = data),
    error = function(e) NULL
  )
  mm_colnames <- if (!is.null(mm)) colnames(mm) else character(0)

  match_design_col <- function(var_name, available) {
    if (var_name %in% available) return(var_name)
    if (grepl(":", var_name, fixed = TRUE)) {
      parts <- strsplit(var_name, ":", fixed = TRUE)[[1]]
      reversed <- paste(rev(parts), collapse = ":")
      if (reversed %in% available) return(reversed)
    }
    NA_character_
  }

  build_aligned <- function(vars) {
    if (length(vars) == 0) return(NULL)
    beta_idx <- integer(0)
    Xcols <- list()
    for (v in vars) {
      coef_match <- match_design_col(v, fname_short)
      if (is.na(coef_match)) next
      xcol <- NULL
      if (!is.null(mm)) {
        dm <- match_design_col(v, mm_colnames)
        if (!is.na(dm)) xcol <- mm[, dm]
      }
      if (is.null(xcol) && v %in% colnames(data)) xcol <- data[[v]]
      if (is.null(xcol)) next
      beta_idx <- c(beta_idx, which(fname_short == coef_match)[1])
      Xcols[[length(Xcols) + 1]] <- as.numeric(xcol)
    }
    if (length(beta_idx) == 0) return(NULL)
    X <- do.call(cbind, Xcols)
    X[is.na(X)] <- 0
    list(beta_idx = beta_idx, X = X)
  }

  aln_w <- build_aligned(within_vars)
  aln_b <- build_aligned(between_vars)

  n_s <- length(draw_ids)
  prop_f1 <- numeric(n_s)
  prop_f2 <- numeric(n_s)
  prop_f3 <- numeric(n_s)
  prop_f4 <- numeric(n_s)
  prop_resid <- numeric(n_s)

  for (s in seq_len(n_s)) {
    var_f1 <- if (!is.null(aln_w)) {
      beta_w <- fix_draws[s, aln_w$beta_idx, drop = FALSE]
      as.numeric(var(as.numeric(aln_w$X %*% t(beta_w))))
    } else 0

    var_f2 <- if (!is.null(aln_b)) {
      beta_b <- fix_draws[s, aln_b$beta_idx, drop = FALSE]
      as.numeric(var(as.numeric(aln_b$X %*% t(beta_b))))
    } else 0

    v_f3 <- var_f3[s]
    v_f4 <- var_f4[s]
    v_resid <- var_resid[s]

    var_total <- var_f1 + var_f2 + v_f3 + v_f4 + v_resid

    prop_f1[s] <- var_f1 / var_total
    prop_f2[s] <- var_f2 / var_total
    prop_f3[s] <- v_f3 / var_total
    prop_f4[s] <- v_f4 / var_total
    prop_resid[s] <- v_resid / var_total
  }

  comp_names <- c("f1_within", "f2_between", "f3_slope_var",
                  "f4_intercept_var", "sigma2_resid", "fvm_total")
  prop_fvm <- 1 - prop_resid

  result_rows <- list()
  for (comp_idx in seq_along(comp_names)) {
    cname <- comp_names[comp_idx]
    vals <- switch(cname,
                  f1_within = prop_f1,
                  f2_between = prop_f2,
                  f3_slope_var = prop_f3,
                  f4_intercept_var = prop_f4,
                  sigma2_resid = prop_resid,
                  fvm_total = prop_fvm)
    result_rows[[comp_idx]] <- data.frame(
      component = cname,
      median = median(vals),
      q2.5   = quantile(vals, 0.025),
      q97.5  = quantile(vals, 0.975),
      stringsAsFactors = FALSE
    )
  }

  do.call(rbind, result_rows)
}

# ============================================================================
# Execute for all 4 RS models
# ============================================================================
r2_rs_rows <- list()

for (mname in rs_models) {
  if (is.null(brms_fits[[mname]])) {
    cat(sprintf("  [%s] fit not in brms_fits, skipping\n", mname))
    next
  }

  status_row <- convergence[[mname]]
  if (!is.null(status_row) && status_row$status == "FAIL") {
    cat(sprintf("  [%s] FAIL — skipping\n", mname))
    next
  }

  cat(sprintf("\n  [%s] Computing 6-component decomposition...\n", mname))

  fit <- brms_fits[[mname]]
  vmap <- rs_var_map[[mname]]

  fit_data <- tryCatch(
    brms::model.frame(fit),
    error = function(e) {
      if (mname %in% c("B1_RS_PHQ", "B1_RS_PC", "B1_RS_sigmaM")) sample_A4
      else sample_B
    }
  )

  tryCatch({
    decomp <- rs_r2_decomp_6comp(
      fit,
      data = fit_data,
      within_vars = vmap$within,
      between_vars = vmap$between,
      slope_var = "health_session_c",
      n_draws = 500L
    )
    decomp$model <- mname
    r2_rs_rows[[mname]] <- decomp[, c("model", "component", "median", "q2.5", "q97.5")]

    cat(sprintf("    f1_within=%.3f, f2_between=%.3f, f3_slope=%.3f, f4_intercept=%.3f, resid=%.3f, fvm=%.3f\n",
                decomp$median[decomp$component == "f1_within"],
                decomp$median[decomp$component == "f2_between"],
                decomp$median[decomp$component == "f3_slope_var"],
                decomp$median[decomp$component == "f4_intercept_var"],
                decomp$median[decomp$component == "sigma2_resid"],
                decomp$median[decomp$component == "fvm_total"]))

    # Self-check: sum ≈ 1
    sum_check <- sum(decomp$median[decomp$component != "fvm_total"])
    cat(sprintf("    Sum check: %.6f (should be ~1.0)\n", sum_check))

    # Cross-check against bayes_R2
    r2_total <- tryCatch(brms::bayes_R2(fit), error = function(e) NULL)
    if (!is.null(r2_total) && "R2" %in% colnames(r2_total)) {
      r2_med <- median(r2_total[, "R2"])
      fvm_med <- decomp$median[decomp$component == "fvm_total"]
      cat(sprintf("    bayes_R2 median=%.3f, fvm_total=%.3f, diff=%.3f\n",
                  r2_med, fvm_med, abs(r2_med - fvm_med)))
    } else {
      cat("    bayes_R2 cross-check skipped (internal error)\n")
    }
  }, error = function(e) {
    cat(sprintf("    ⚠ Decomposition failed: %s\n", e$message))
  })
}

# ============================================================================
# Save
# ============================================================================
if (length(r2_rs_rows) > 0) {
  r2_rs_combined <- do.call(rbind, r2_rs_rows)
  rownames(r2_rs_combined) <- NULL
  out_path <- file.path(tables_dir, "bayes_r2_decomposition_rs.csv")
  write.csv(r2_rs_combined, out_path, row.names = FALSE)
  cat(sprintf("\n  ✓ RS 6-component decomposition saved: %s\n", out_path))

  cat("\n  RS Variance Decomposition Summary:\n")
  pivot_df <- tidyr::pivot_wider(
    r2_rs_combined[, c("model", "component", "median")],
    names_from = "component", values_from = "median"
  )
  print(pivot_df, row.names = FALSE, digits = 3)
} else {
  cat("\n  ✗ No RS models decomposed\n")
}

cat("\n  ✓ 模块06d RS 方差分解完成\n")
