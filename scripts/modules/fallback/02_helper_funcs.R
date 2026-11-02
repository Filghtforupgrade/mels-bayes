# ============================================================================
# 模块 02: 辅助函数 — 建模诊断 + r2mlm验证 + 收敛诊断 + 缺失分析
# ============================================================================
# 依赖: 00_config.R + 01_data_prep.R
# 产出: extract_coef_v71, get_model_series, print_model_diag,
#       compute_r2mlm_flagship, extended_convergence_check,
#       generate_convergence_report, verify_r2mlm_output,
#       analyze_missing_patterns
# ============================================================================

cat("\n--- 模块02: 辅助函数定义 ---\n")

# ============================================================================
# 2A. 系数提取
# ============================================================================
extract_coef_v71 <- function(model, model_name, vars_of_interest = NULL) {
  s <- coef(summary(model))
  
  if (is.null(vars_of_interest)) {
    vars_of_interest <- rownames(s)[-1]
  }
  
  results_list <- list()
  
  for (v in vars_of_interest) {
    if (v %in% rownames(s)) {
      beta <- s[v, 1]
      se <- s[v, 2]
      if (ncol(s) >= 5) {
        pv <- s[v, 5]
      } else {
        pv <- 2 * (1 - pt(abs(s[v, 4]), df = model@devcomp$dims["nmp"]))
      }
      
      ci_l <- beta - 1.96 * se
      ci_u <- beta + 1.96 * se
      
      sig <- if (pv < 0.001) "***" else if (pv < 0.01) "**" else if (pv < 0.05) "*" else if (pv < 0.10) "†" else ""
      
      results_list[[v]] <- data.frame(
        model = model_name,
        series = get_model_series(model_name),
        predictor = v,
        beta = beta,
        se = se,
        p_value = pv,
        ci_lower = ci_l,
        ci_upper = ci_u,
        significance = sig,
        stringsAsFactors = FALSE
      )
    }
  }
  
  return(do.call(rbind, results_list))
}

# ============================================================================
# 2B. 模型系列推断
# ============================================================================
get_model_series <- function(name) {
  if (startsWith(name, "ICC") || startsWith(name, "Null")) {
    return("ICC_null_model")
  } else if (startsWith(name, "G")) {
    return("G_growth_trajectory")
  } else if (startsWith(name, "B")) {
    return("B_behavior_combination")
  } else if (startsWith(name, "R")) {
    return("R_emotion_regulation")
  } else if (startsWith(name, "CL")) {
    return("CL_cross_lagged")
  } else if (startsWith(name, "S")) {
    return("S_random_slope")
  } else if (startsWith(name, "RE")) {
    return("RE_reverse_prediction")
  } else if (startsWith(name, "A")) {
    return("A_diagnostic")
  } else {
    return("other")
  }
}

# ============================================================================
# 2C. 模型诊断打印
# ============================================================================
print_model_diag <- function(model, model_name, calculate_vif = TRUE) {
  diag <- data.frame(
    model = model_name,
    converged = ifelse(is.null(model@optinfo$conv$code) || model@optinfo$conv$code == 0, TRUE, FALSE),
    n_obs = nobs(model),
    aic = AIC(model),
    bic = BIC(model),
    llf = as.numeric(logLik(model)),
    max_vif = NA,
    stringsAsFactors = FALSE
  )
  
  if (calculate_vif && length(fixef(model)) > 3) {
    tryCatch({
      vif_vals <- car::vif(model)
      if (is.numeric(vif_vals)) {
        diag$max_vif <- max(vif_vals)
      }
    }, error = function(e) {})
  }
  
  cat(sprintf("\n--- %s ---\n", model_name))
  cat(sprintf("  收敛: %s\n", ifelse(diag$converged, "✅", "❌")))
  cat(sprintf("  观测: %d\n", diag$n_obs))
  cat(sprintf("  AIC: %.2f, BIC: %.2f\n", diag$aic, diag$bic))
  
  if (!is.na(diag$max_vif) && diag$max_vif > 1) {
    cat(sprintf("  最大VIF: %.2f %s\n", diag$max_vif, ifelse(diag$max_vif < 5, "✅", "⚠️")))
  }
  
  return(diag)
}

# ============================================================================
# 2D. r2mlm旗舰模型R²分解
# ============================================================================
compute_r2mlm_flagship <- function(model, model_name, within_vars, between_vars = NULL,
                                    random_vars = NULL, data_obj = NULL, verbose = TRUE) {
  if (is.null(data_obj)) data_obj <- model@frame
  
  result_entry <- list(
    model = model_name,
    success = FALSE,
    f1_total = NA, f1_within = NA,
    f2_total = NA, f2_between = NA,
    v_total = NA, v_within = NA,
    m_total = NA, m_between = NA,
    sigma2 = NA,
    fvm_total = NA,
    note = ""
  )
  
  tryCatch({
    gamma_w <- fixef(model)[within_vars]
    gamma_w <- gamma_w[!is.na(gamma_w)]
    
    gamma_b <- fixef(model)["(Intercept)"]
    if (!is.null(between_vars) && length(between_vars) > 0) {
      gamma_b <- c(gamma_b, fixef(model)[between_vars])
    }
    
    Tau <- as.matrix(VarCorr(model)$user_id)
    if (length(Tau) == 0 || nrow(Tau) == 0) {
      Tau <- matrix(attr(VarCorr(model)$user_id, "stddev")^2, 1, 1)
    }
    
    sigma2 <- attr(VarCorr(model), "sc")^2
    
    if (is.null(between_vars) || length(between_vars) == 0) {
      between_vars_input <- NULL
    } else {
      between_vars_input <- between_vars
    }
    
    if (is.null(random_vars) || length(random_vars) == 0) {
      random_vars_input <- NULL
    } else {
      random_vars_input <- random_vars
    }
    
    r2_out <- r2mlm_manual(
      data = data_obj,
      within_covs = within_vars,
      between_covs = between_vars_input,
      random_covs = random_vars_input,
      gamma_w = as.numeric(gamma_w),
      gamma_b = as.numeric(gamma_b),
      Tau = Tau,
      sigma2 = sigma2,
      has_intercept = TRUE,
      clustermeancentered = TRUE,
      bargraph = FALSE
    )
    
    if (is.list(r2_out) && "R2s" %in% names(r2_out)) {
      r2s <- r2_out$R2s
      result_entry$success <- TRUE
      result_entry$f1_total  <- as.numeric(r2s["f1", "total"])
      result_entry$f1_within <- as.numeric(r2s["f1", "within"])
      result_entry$f2_total  <- as.numeric(r2s["f2", "total"])
      result_entry$f2_between<- as.numeric(r2s["f2", "between"])
      result_entry$v_total   <- as.numeric(r2s["v", "total"])
      result_entry$v_within  <- as.numeric(r2s["v", "within"])
      result_entry$m_total   <- as.numeric(r2s["m", "total"])
      result_entry$m_between <- as.numeric(r2s["m", "between"])
      result_entry$sigma2    <- 1.0 - as.numeric(r2s["fvm", "total"])
      result_entry$fvm_total <- as.numeric(r2s["fvm", "total"])
    }
  }, error = function(e) {
    result_entry$note <- paste("r2mlm失败:", e$message)
  })
  
  if (verbose) {
    if (result_entry$success) {
      cat(sprintf("  ★r2mlm R²[%s]: f1_tot=%.3f f1_within=%.3f m_tot=%.3f m_between=%.3f\n",
                  model_name, result_entry$f1_total, result_entry$f1_within,
                  result_entry$m_total, result_entry$m_between))
    } else {
      cat(sprintf("  ⚠ r2mlm[%s]: %s\n", model_name, result_entry$note))
    }
  }
  
  return(result_entry)
}

# ============================================================================
# 2E.  M-5: 扩展收敛诊断函数
# ============================================================================
extended_convergence_check <- function(model, model_name) {
  diag_entry <- data.frame(
    model = model_name,
    converged = NA,
    conv_code = NA,
    max_gradient = NA,
    singular_fit = NA,
    hessian_pd = NA,
    n_obs = NA,
    boundary_warning = NA,
    message = "",
    stringsAsFactors = FALSE
  )
  
  tryCatch({
    diag_entry$n_obs <- nobs(model)
    
    conv_code <- model@optinfo$conv$opt
    diag_entry$conv_code <- ifelse(is.null(conv_code), 0, conv_code)
    diag_entry$converged <- ifelse(is.null(conv_code) || conv_code == 0, TRUE, FALSE)
    
    if (!is.null(model@optinfo$derivs) && length(model@optinfo$derivs$gradient) > 0) {
      diag_entry$max_gradient <- max(abs(model@optinfo$derivs$gradient))
    }
    
    diag_entry$singular_fit <- isSingular(model, tol = 1e-4)
    
    if (!is.null(model@optinfo$derivs$Hessian) && length(model@optinfo$derivs$Hessian) > 0) {
      hess_eig <- tryCatch(eigen(model@optinfo$derivs$Hessian)$values, error = function(e) NULL)
      if (!is.null(hess_eig)) {
        diag_entry$hessian_pd <- !any(hess_eig <= 0)
        if (any(hess_eig <= 0)) {
          diag_entry$boundary_warning <- TRUE
        }
      }
    }
    
    if (!diag_entry$converged) {
      diag_entry$message <- paste("收敛代码:", diag_entry$conv_code)
      if (!is.na(diag_entry$max_gradient) && diag_entry$max_gradient > 0.001) {
        diag_entry$message <- paste(diag_entry$message, 
          sprintf("| 最大梯度: %.6f", diag_entry$max_gradient))
      }
    }
    if (diag_entry$singular_fit) {
      diag_entry$message <- paste(diag_entry$message, "| ⚠边界拟合(singular fit)")
    }
    if (!is.na(diag_entry$hessian_pd) && !diag_entry$hessian_pd) {
      diag_entry$message <- paste(diag_entry$message, "| ⚠Hessian非正定")
    }
    
  }, error = function(e) {
    diag_entry$message <- paste("诊断提取失败:", e$message)
  })
  
  return(diag_entry)
}

# ============================================================================
# 2F.  M-5: 批量收敛诊断报告生成
# ============================================================================
generate_convergence_report <- function(model_list) {
  report <- do.call(rbind, lapply(names(model_list), function(mname) {
    extended_convergence_check(model_list[[mname]], mname)
  }))
  
  n_total <- nrow(report)
  n_converged <- sum(report$converged, na.rm = TRUE)
  n_singular <- sum(report$singular_fit, na.rm = TRUE)
  n_warning <- sum(!report$converged | report$singular_fit | 
                   (!is.na(report$hessian_pd) & !report$hessian_pd), na.rm = TRUE)
  
  cat(sprintf("\n  ★收敛诊断报告: %d个模型\n", n_total))
  cat(sprintf("    成功收敛: %d/%d\n", n_converged, n_total))
  cat(sprintf("    边界拟合警告(singular): %d\n", n_singular))
  cat(sprintf("    任何警告: %d\n", n_warning))
  
  if (n_warning > 0) {
    warn_models <- report %>%
      filter(!converged | singular_fit | (!is.na(hessian_pd) & !hessian_pd))
    cat("    ⚠ 警告模型:\n")
    for (i in 1:nrow(warn_models)) {
      cat(sprintf("      - %s: %s\n", warn_models$model[i], warn_models$message[i]))
    }
  }
  
  return(report)
}

# ============================================================================
# 2G.  M-6: r2mlm输出验证
# ============================================================================
verify_r2mlm_output <- function(r2mlm_entry, model_name, tolerance = 0.001) {
  verify <- data.frame(
    model = model_name,
    f1_total = r2mlm_entry$f1_total,
    f2_total = r2mlm_entry$f2_total,
    m_total  = r2mlm_entry$m_total,
    v_total  = r2mlm_entry$v_total,
    sigma2   = r2mlm_entry$sigma2,
    fvm_total = r2mlm_entry$fvm_total,
    sum_check = NA,
    sigma2_match = NA,
    fvm_match = NA,
    all_pass = NA,
    note = "",
    stringsAsFactors = FALSE
  )
  
  if (r2mlm_entry$success) {
    components <- c(r2mlm_entry$f1_total, r2mlm_entry$f2_total, 
                    r2mlm_entry$v_total, r2mlm_entry$m_total)
    components_clean <- components[!is.na(components)]
    
    if (length(components_clean) == 4) {
      fvm_computed <- sum(components_clean)
      verify$sum_check <- fvm_computed + r2mlm_entry$sigma2
      verify$sigma2_match <- abs(verify$sum_check - 1.0) < tolerance
      
      if (!is.na(fvm_computed) && !is.na(r2mlm_entry$fvm_total)) {
        verify$fvm_match <- abs(fvm_computed - r2mlm_entry$fvm_total) < tolerance
      }
      
      if (!verify$sigma2_match) {
        verify$note <- sprintf("f1+f2+v+m+sigma2=%.4f (偏离1.0达%.6f)", 
                               verify$sum_check, verify$sum_check - 1.0)
      }
      if (!is.na(verify$fvm_match) && !verify$fvm_match) {
        verify$note <- paste(verify$note, sprintf("| fvm不匹配(手动=%.4f, 包=%.4f)", 
                             fvm_computed, r2mlm_entry$fvm_total))
      }
      
      verify$all_pass <- verify$sigma2_match && 
                         (is.na(verify$fvm_match) || verify$fvm_match)
    } else {
      verify$note <- "分量不完整，跳过验证"
    }
  } else {
    verify$note <- "r2mlm计算失败"
  }
  
  return(verify)
}

# ============================================================================
# 2H.  M-4: 缺失数据模式分析
# ============================================================================
analyze_missing_patterns <- function(df, behavior_cols, session_raw = NULL) {
  cat("\n★ M-4: 缺失数据模式分析 (EMA设计结构审计)\n")
  cat("  关键发现: EMA将8个行为分时段收集,非缺失问题\n\n")
  
  missing_report <- list()
  
  if (!is.null(session_raw)) {
    raw_with_mood <- session_raw %>% filter(!is.na(mood))
    
    cat("  ★EMA设计结构: 行为分时段收集\n")
    cat("  ─────────────────────────────────────────────\n")
    cat("  noon时段收集: breakfast, sleep_late, wake_early, aspiration\n")
    cat("  evening时段收集: walk, exercise, dinner, breathing\n\n")
    
    cat("  (a) 各时段session数及行为覆盖:\n")
    
    behavior_info <- data.frame()
    for (s in unique(raw_with_mood$session)) {
      s_data <- raw_with_mood %>% filter(session == s)
      n_s <- nrow(s_data)
      cat(sprintf("\n  --- %s (n=%d sessions) ---\n", s, n_s))
      
      for (b in behavior_cols) {
        n_na <- sum(is.na(s_data[[b]]))
        n_asked <- n_s - n_na
        n_done <- sum(s_data[[b]] > 0, na.rm = TRUE)
        n_not_done <- sum(s_data[[b]] == 0, na.rm = TRUE)
        
        if (n_na == n_s) {
          cat(sprintf("    %s: 该时段未询问 (NA=%d/%d)\n", b, n_na, n_s))
        } else {
          exec_rate <- 100 * n_done / n_asked
          cat(sprintf("    %s: 执行=%d/%d(%.1f%%), 未执行=%d\n", 
                      b, n_done, n_asked, exec_rate, n_not_done))
          behavior_info <- rbind(behavior_info, data.frame(
            behavior = b,
            session = s,
            n_asked = n_asked,
            n_done = n_done,
            n_not_done = n_not_done,
            exec_rate = exec_rate,
            stringsAsFactors = FALSE
          ))
        }
      }
    }
    
    cat(sprintf("\n  ─────────────────────────────────────────────\n"))
    cat(sprintf("  (b) 汇总: 各行为在\"被询问的时段\"中的执行率\n\n"))
    behavior_summary <- behavior_info %>%
      group_by(behavior) %>%
      summarise(
        n_asked_total = sum(n_asked),
        n_done_total = sum(n_done),
        exec_rate = 100 * n_done_total / n_asked_total,
        .groups = 'drop'
      ) %>%
      arrange(desc(exec_rate))
    for (i in 1:nrow(behavior_summary)) {
      cat(sprintf("    %-14s 执行率=%5.1f%% (%d/%d)\n",
                  behavior_summary$behavior[i],
                  behavior_summary$exec_rate[i],
                  behavior_summary$n_done_total[i],
                  behavior_summary$n_asked_total[i]))
    }
    
    missing_report[["behavior_exec_rates"]] <- behavior_summary
    missing_report[["behavior_by_session"]] <- behavior_info
    
    cat(sprintf("\n  (c) 参与者在被询问时至少执行1个行为:\n"))
    raw_with_mood <- raw_with_mood %>%
      mutate(
        noon_behaviors_done = ifelse(session == "noon",
          rowSums(select(., intersect(c("breakfast","sleep_late","wake_early","aspiration"), names(.))) > 0, na.rm = TRUE), NA),
        evening_behaviors_done = ifelse(session == "evening",
          rowSums(select(., intersect(c("walk","exercise","dinner","breathing"), names(.))) > 0, na.rm = TRUE), NA)
      )
    cat(sprintf("    noon 至少1行为: %.1f%%\n",
                100 * mean(raw_with_mood$noon_behaviors_done >= 1, na.rm = TRUE)))
    cat(sprintf("    evening 至少1行为: %.1f%%\n",
                100 * mean(raw_with_mood$evening_behaviors_done >= 1, na.rm = TRUE)))
    
    cat(sprintf("\n  (d) 补充说明:\n"))
    cat("    item_missing_n 反映了EMA设计结构(每时段固定询问4个行为),\n")
    cat("    并非参与者不响应。本次EMA未收集行为质量/持续时间/强度等维度。\n")
    cat("    每个行为在其被询问时段内的执行率见上表,代表了真实的\n")
    cat("    \"参与者执行了该行为\"的比例(而非\"缺失响应\")。\n")
  } else {
    cat("  ⚠ session_raw未传入,跳过结构分析\n")
    missing_report <- list()
  }
  
  return(missing_report)
}

cat("\n  ✓ 模块02 辅助函数定义完成\n")
cat(sprintf("    extract_coef_v71 | get_model_series | print_model_diag\n"))
cat(sprintf("    compute_r2mlm_flagship | extended_convergence_check\n"))
cat(sprintf("    generate_convergence_report | verify_r2mlm_output\n"))
cat(sprintf("    analyze_missing_patterns\n"))
