#!/usr/bin/env Rscript
# A15: SESOI 后验概率 — 对标 Dora R L352-355
options(width=200); options(warn=1)
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
cat("\n═══ A15: SESOI 后验概率 ═══\n")

compute_sesoi <- function(fit, param, minimal, moderate, model_name) {
  draws <- posterior::as_draws_df(fit)
  if(!(param %in% names(draws))) return(NULL)
  vals <- draws[[param]]
  n <- length(vals)
  is_sigma <- grepl("^b_sigma_", param)
  if(is_sigma) {
    p_min <- mean(vals < -minimal)
    p_mod <- mean(vals < -moderate)
  } else {
    p_min <- mean(vals > minimal)
    p_mod <- mean(vals > moderate)
  }
  data.frame(model=model_name, param=param, median=median(vals),
             q2.5=quantile(vals,.025), q97.5=quantile(vals,.975),
             P_positive=mean(vals>0), P_meaningful_min=round(p_min,4), P_meaningful_mod=round(p_mod,4),
             SESOI_min=minimal, SESOI_mod=moderate, n_draws=n)
}

# Load fits
B1_MELS <- readRDS(file.path(fits_dir,"B1_MELS.rds"))
B1_RS_PHQ <- readRDS(file.path(fits_dir,"B1_RS_PHQ.rds"))
B1_RS_PC <- readRDS(file.path(fits_dir,"B1_RS_PC.rds"))
B1_RS_sigmaM <- readRDS(file.path(fits_dir,"B1_RS_sigmaM.rds"))

sesoi_all <- rbind(
  compute_sesoi(B1_MELS, "b_health_session_c", 0.05, 0.50, "B1_MELS"),
  compute_sesoi(B1_MELS, "b_sigma_health_session_c", 0.03, 0.10, "B1_MELS"),
  compute_sesoi(B1_RS_PHQ, "b_health_session_c:pre_phq9_c", 0.02, 0.10, "B1_RS_PHQ"),
  compute_sesoi(B1_RS_PC, "b_health_session_c:pre_control_c", 0.02, 0.10, "B1_RS_PC"),
  compute_sesoi(B1_RS_sigmaM, "b_health_session_c:person_std_mood_c", 0.02, 0.10, "B1_RS_sigmaM")
)

write.csv(sesoi_all[,c("model","param","median","q2.5","q97.5","P_positive","P_meaningful_min","P_meaningful_mod","SESOI_min","SESOI_mod")],
          file.path(output_dir,"posterior_probability_sesoi.csv"), row.names=FALSE)

cat("\n--- SESOI 后验概率 ---\n")
print(sesoi_all[,c("param","median","q2.5","q97.5","P_meaningful_min","P_meaningful_mod")], row.names=FALSE, digits=3)

# MD table
write_sesoi_md <- function(df, path) {
  con <- file(path, "wt")
  tryCatch({
    cat("# Table S14: Posterior Probability of Meaningful Effects (SESOI)\n\n", file=con)
    cat("| Model | Parameter | Median [95% CrI] | P(>SESOI minimal) | P(>SESOI moderate) | SESOI minimal | SESOI moderate |\n", file=con)
    cat("|:------|:----------|:----------------|:-----------------:|:------------------:|:-------------:|:---------------:|\n", file=con)
    for(i in seq_len(nrow(df))) {
      r <- df[i,]
      cat(sprintf("| %s | %s | %.3f [%.3f, %.3f] | %.3f | %.3f | %.3f | %.3f |\n",
                  r$model, r$param, r$median, r$q2.5, r$q97.5,
                  r$P_meaningful_min, r$P_meaningful_mod, r$SESOI_min, r$SESOI_mod), file=con)
    }
  }, finally = close(con))
}
write_sesoi_md(sesoi_all, file.path(tables_dir,"table_14_sesoi.md"))
cat("\n═══ A15 完成 ═══\n")
