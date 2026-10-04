#!/usr/bin/env Rscript
# =============================================================================
# nDT upper-bound sensitivity ladder — engine version (corrected density)
# =============================================================================
# WHY.  Methods has to justify the nDT upper bound of 3.0 s, and the justification
# currently on file predates the 2026-08-07 transition-density fix.  Two facts make
# the bound a live reviewer target:
#   * at 1.5 s the ceiling censored ~74% of participants, which is why the suite
#     moved to 3.0 (`ndt-bound-forensic-and-venue`);
#   * at 3.0 s the ceiling is NOT gone -- 14.2% still pile up (tolerance 0.025-0.03;
#     8.3% at 0.001, 16.0% at 0.05), asymmetrically by group (HC 10.2%, DEP 19.9%),
#     and the nDT group difference does not survive dropping them (p .013 -> .245).
# The question this answers is therefore NOT "does tau die past 2.0" (that was the
# old, pre-fix narrative and it is wrong: tau survives the ceiling-drop test at
# p = .024).  It is: **is 3.0 a defensible choice, or does the tau result move with
# the bound?**  Everything the manuscript claims about tau should be flat across
# bounds once the pile-up clears.
#
# Refits the MAIN 4j (mean-first) model at a ladder of nDT upper bounds and reports
# per bound: ceiling pile-up (three tolerances), tau group difference, ANCOVA with
# and without median RT, and r(tau, median RT).
#
#   Rscript scripts/tDDM/nDT_bound_ladder_engine.R
#   Rscript scripts/tDDM/nDT_bound_ladder_engine.R --bounds 1.5,2.0,3.0 --cores 8
# =============================================================================
suppressPackageStartupMessages({ library(DEoptim); library(parallel) })
args <- commandArgs(trailingOnly = TRUE)
ga <- function(f, d) { if (f %in% args) { i <- which(args == f); if (i < length(args)) return(args[i+1]) }; d }

BOUNDS  <- as.numeric(strsplit(ga("--bounds", "1.5,2.0,2.5,3.0,3.5"), ",")[[1]])
DE_ITER <- as.integer(ga("--de-iter", 200))     # matched to refit_corrected_likelihood.R
DE_NP   <- as.integer(ga("--de-np", 80))
NC      <- as.integer(ga("--cores", max(1, parallel::detectCores() - 2)))
# 2026-09-23: keep the BEST-NLL run of R repeats, matching the canonical procedure
# (docs/PAPER/Memo_Analyses_0919.md S4).
RESTARTS <- as.integer(ga("--restarts", 1))

sd_ <- tryCatch(dirname(sys.frame(1)$ofile), error = function(e) {
  ca <- commandArgs(FALSE); fa <- sub("^--file=", "", ca[grep("^--file=", ca)])
  if (length(fa)) dirname(fa) else getwd() })
source(file.path(sd_, "..", "..", "config.R"))
source(file.path(MODEL_DIR, "tDDM", "tddm_engine.R"))

S1 <- "mean_sig"; S2 <- "contrast_sig"            # mean-first, the manuscript model
LL <- mk_ll_tau(S1, S2, trans = get_transition("image"))

data <- tddm_load_study4(DATA_DIR)
subjects <- unique(data$subject)
medRT <- tapply(data$RT_ddm, data$subject, median)

cat(sprintf("=== nDT-bound ladder (4j mean-first, CORRECTED density) ===\n"))
cat(sprintf("  bounds = {%s} | %d subjects | DE %d/%d | restarts = %d | cores = %d\n\n",
            paste(BOUNDS, collapse = ", "), length(subjects), DE_ITER, DE_NP, RESTARTS, NC))

ctrl <- DEoptim.control(itermax = DE_ITER, NP = DE_NP, trace = FALSE)   # NEVER set VTR

fit_one <- function(s, ub) {
  td <- data[data$subject == s, ]
  if (nrow(td) < 5) return(NULL)
  f <- NULL
  for (k in seq_len(RESTARTS)) {
    fk <- tryCatch(DEoptim(LL, lower = c(-3, -3, 0.5, 0.01, -2),
                           upper = c(3, 3, 4, ub, 3), control = ctrl, td = td),
                   error = function(e) NULL)
    if (!is.null(fk) && (is.null(f) || fk$optim$bestval < f$optim$bestval)) f <- fk
  }
  if (is.null(f)) return(NULL)
  e <- f$optim$bestmem
  data.frame(subject = s, group = td$group[1], n_trials = nrow(td),
             negLogLik = f$optim$bestval, de_iter = f$optim$iter,
             d_mean = e[1], d_contrast = e[2], thres = e[3], nDT = e[4], tau = max(e[5], 0), tau_raw = e[5],   # 2026-09-28: tau<=0 -> 0
             stringsAsFactors = FALSE)
}

ts <- format(Sys.time(), "%Y%m%d_%H%M%S")
od <- file.path(RESULTS_DIR, "tDDM", sprintf("nDT_bound_ladder_engine_%s", ts))
dir.create(od, recursive = TRUE, showWarnings = FALSE)

t0 <- Sys.time(); S <- list()
for (ub in BOUNDS) {
  cat(sprintf("--- nDT bound %.1f ---\n", ub)); flush.console()
  res <- do.call(rbind, Filter(is.data.frame, mclapply(subjects, fit_one, ub = ub, mc.cores = NC)))
  stopifnot(min(res$de_iter) > 0)                                        # VTR guard
  write.csv(res, file.path(od, sprintf("fits_4j_nDT%.1f.csv", ub)), row.names = FALSE)

  res$medRT <- as.numeric(medRT[as.character(res$subject)])
  ceil <- function(tol) 100 * mean(res$nDT >= ub - tol)
  hc <- res$tau[res$group == "HC"]; dp <- res$tau[res$group == "DEP"]
  tt <- t.test(dp, hc)
  m1 <- summary(lm(tau ~ group, data = res))$coefficients["groupHC", ]
  m2 <- summary(lm(tau ~ group + medRT, data = res))$coefficients["groupHC", ]
  S[[length(S) + 1]] <- data.frame(
    bound = ub, nDT_median = median(res$nDT),
    ceiling_001 = ceil(0.001), ceiling_01 = ceil(0.01), ceiling_05 = ceil(0.05),
    ceiling_HC_05 = 100 * mean(res$nDT[res$group == "HC"]  >= ub - 0.05),
    ceiling_DEP_05 = 100 * mean(res$nDT[res$group == "DEP"] >= ub - 0.05),
    dtau = mean(dp) - mean(hc), p_dtau = tt$p.value,
    beta_grp = m1[1], p_grp = m1[4],
    beta_grp_medRT = m2[1], p_grp_medRT = m2[4],
    r_tau_medRT = cor(res$tau, res$medRT), tau_mean = mean(res$tau))
  x <- S[[length(S)]]
  cat(sprintf("  nDT med %.2f | ceiling %.0f%% (HC %.0f / DEP %.0f) | Dtau %+.3f p=%.4f | b(grp|medRT) %+.3f p=%.3f | r(tau,medRT) %.2f  [%.1f min]\n",
              x$nDT_median, x$ceiling_05, x$ceiling_HC_05, x$ceiling_DEP_05,
              x$dtau, x$p_dtau, x$beta_grp_medRT, x$p_grp_medRT, x$r_tau_medRT,
              as.numeric(difftime(Sys.time(), t0, units = "mins"))))
  flush.console()
}
S <- do.call(rbind, S)
write.csv(S, file.path(od, "ladder_summary.csv"), row.names = FALSE)
cat("\n=== LADDER SUMMARY ===\n")
cat(sprintf("%-7s %-9s %-22s %-18s %-18s %-8s\n",
            "bound", "nDT_med", "ceiling%@.05 (HC/DEP)", "Dtau (p)", "b(grp|medRT) (p)", "r(t,RT)"))
for (i in seq_len(nrow(S))) cat(sprintf("%-7.1f %-9.2f %-22s %+.3f (%.4f)   %+.3f (%.3f)     %.2f\n",
  S$bound[i], S$nDT_median[i],
  sprintf("%.0f (%.0f/%.0f)", S$ceiling_05[i], S$ceiling_HC_05[i], S$ceiling_DEP_05[i]),
  S$dtau[i], S$p_dtau[i], S$beta_grp_medRT[i], S$p_grp_medRT[i], S$r_tau_medRT[i]))
cat(sprintf("\nSaved: %s\nTotal: %.1f min\n", od, as.numeric(difftime(Sys.time(), t0, units = "mins"))))
