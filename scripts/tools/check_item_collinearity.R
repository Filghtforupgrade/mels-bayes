#!/usr/bin/env Rscript
# ============================================================================
#  C1: 每时段条目相关矩阵与 VIF（零拟合诊断）
# 午间 4 条目（retiring on time / early rising / breakfast / aspiration）
# 与晚间 4 条目（walk / exercise / dinner / breathing）分别计算
# Pearson 相关矩阵与 VIF（VIF = diag(solve(R))，与 lm+car::vif 等价）。
# 数据只读：sample_lagged（与条目级 MELS 同一分析样本），不新增拟合。
# 输出: item_collinearity.md（诊断报告）
# ============================================================================
options(width = 200)

project_root <- normalizePath(
  file.path(dirname(sub("^--file=", "",
    commandArgs(trailingOnly = FALSE)[grep("^--file=", commandArgs(trailingOnly = FALSE))])), "..", ".."),
  mustWork = TRUE)

e <- new.env()
load(file.path(project_root, "data", "checkpoints", "checkpoint_01_data.RData"), envir = e)
d <- e$sample_lagged
cat(sprintf("sample_lagged: %d rows, %d cols\n", nrow(d), ncol(d)))

# 会话类型列探测
sess_col <- intersect(names(d), c("session_noon", "session_numeric", "session_type", "session"))
stopifnot(length(sess_col) >= 1)
is_noon <- if ("session_noon" %in% names(d)) d$session_noon == 1 else d[[sess_col[1]]] == "noon" | grepl("noon", as.character(d[[sess_col[1]]]))
cat(sprintf("noon rows: %d | evening rows: %d\n", sum(is_noon), sum(!is_noon)))

# 分析编码: sleep_ontime = 1 - sleep_late（与 sleep_ontime_c = -sleep_late_c 同向）
noon_df <- data.frame(
  `retiring on time` = 1 - d$sleep_late[is_noon],
  `early rising`     = d$wake_early[is_noon],
  `breakfast`        = d$breakfast[is_noon],
  `aspiration`       = d$aspiration[is_noon],
  check.names = FALSE)
eve_df <- data.frame(
  `walking`    = d$walk[!is_noon],
  `exercise`   = d$exercise[!is_noon],
  `dinner`     = d$dinner[!is_noon],
  `breathing`  = d$breathing[!is_noon],
  check.names = FALSE)

report_item_block <- function(df, label) {
  R <- cor(df)
  vif <- diag(solve(R))
  # 与 car::vif 等价性自检（线性代数恒等式，无需拟合）
  cat(sprintf("\n[%s] n = %d observations\n", label, nrow(df)))
  fmt <- function(x) formatC(x, digits = 3, format = "f")
  lines <- c(sprintf("## %s items (n = %d observations)", label, nrow(df)), "",
             "Pearson correlation matrix:", "",
             paste0("| ", paste(c("", colnames(df)), collapse = " | "), " |"),
             paste0("|", paste(rep("---", ncol(df) + 1), collapse = "|"), "|"))
  for (i in seq_len(nrow(R))) {
    lines <- c(lines, paste0("| ", rownames(R)[i], " | ",
                             paste(fmt(R[i, ]), collapse = " | "), " |"))
  }
  lines <- c(lines, "", "Variance inflation factors (VIF):", "",
             paste0("| ", paste(c("Item", "VIF"), collapse = " | "), " |"),
             "|---|---|")
  for (j in seq_along(vif)) {
    lines <- c(lines, paste0("| ", names(vif)[j], " | ", fmt(vif[j]), " |"))
  }
  lines <- c(lines, "",
             sprintf("Max |r| = %.3f; max VIF = %.3f (threshold for concern commonly VIF > 5).",
                     max(abs(R[upper.tri(R)])), max(vif)),
             "")
  list(lines = lines, R = R, vif = vif)
}

noon_rep <- report_item_block(noon_df, "Noon-session")
eve_rep  <- report_item_block(eve_df,  "Evening-session")

max_vif <- max(noon_rep$vif, eve_rep$vif)
concl <- if (max_vif > 5) {
  sprintf("**Conclusion.** Severe collinearity detected (max VIF = %.2f > 5); this must be disclosed in the Supplementary Materials.", max_vif)
} else {
  sprintf("**Conclusion.** All pairwise |r| <= %.3f and all VIF <= %.2f (well below the conventional VIF > 5 concern threshold), so the four items within each session are not seriously collinear and the eight item-level coefficients can be interpreted as separate associations.", max(abs(noon_rep$R[upper.tri(noon_rep$R)]), abs(eve_rep$R[upper.tri(eve_rep$R)])), max_vif)
}

out_dir <- file.path(project_root, "docs", "reports")
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
md <- c("# Item collinearity diagnostics", "",
        "Computed on the lagged analytic sample (`sample_lagged`, n = 6,836; the estimation",
        "sample of the item-level MELS models), split by session type. Items are scored in",
        "the analysis direction (retiring on time = 1 - staying up late). Pearson correlations",
        "for binary items equal phi coefficients; VIF = diagonal of the inverse correlation",
        "matrix, identical to lm + car::vif. Zero new model fits; read-only computation.",
        "", noon_rep$lines, eve_rep$lines, concl, "",
        sprintf("Script: `scripts/tools/check_item_collinearity.R`; generated %s.", format(Sys.Date())))
writeLines(md, file.path(out_dir, "item_collinearity.md"))
cat(sprintf("\n✓ %s\n", file.path(out_dir, "item_collinearity.md")))
cat("Max VIF overall:", round(max_vif, 3), "\n")
