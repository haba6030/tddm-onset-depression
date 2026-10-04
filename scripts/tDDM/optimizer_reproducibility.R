#!/usr/bin/env Rscript
# =============================================================================
# Optimizer reproducibility of the 4j fit  (G4, 2026-08-20)
# =============================================================================
# WHY.  Two independent runs of the SAME canonical fit (mean-first 4j, corrected
# density, nDT bound 3.0, DE itermax 200 / NP 80) returned different group onset
# differences: +0.325 (corrected_likelihood_20260807_222938) and +0.234 (the
# nDT-bound ladder at 3.0).  The two runs are tied on fit -- mean NLL difference
# +0.010, median 0.000, 144 subjects better in one and 145 in the other -- so this
# is not optimizer quality.  Across the two runs the other parameters are stable
# (r = 0.989 d_mean, 0.994 threshold, 0.998 nDT) while tau is not (r = 0.758).
# The tau likelihood is flat for many participants, so equally good optima sit at
# different tau.  The reported point estimate therefore carries run-to-run
# variability that has never been quantified.
#
# This script fits the canonical model K times with different RNG streams and
# reports the distribution of the group onset difference across runs.  With
# --restarts R > 1 each participant is fitted R times per run and the best fit is
# kept, which is the candidate remedy.
#
#   Rscript scripts/tDDM/optimizer_reproducibility.R --reps 10 --cores 12
#   Rscript scripts/tDDM/optimizer_reproducibility.R --reps 5 --restarts 3
# =============================================================================
suppressPackageStartupMessages({ library(DEoptim); library(parallel) })
args <- commandArgs(trailingOnly = TRUE)
ga <- function(f, d) { if (f %in% args) { i <- which(args == f); if (i < length(args)) return(args[i+1]) }; d }
N_REP    <- as.integer(ga("--reps", 10))
RESTARTS <- as.integer(ga("--restarts", 1))
ITERMAX  <- as.integer(ga("--itermax", 200))
NP       <- as.integer(ga("--np", 80))
NC       <- as.integer(ga("--cores", max(1, parallel::detectCores() - 2)))
SEED0    <- as.integer(ga("--seed", 20260820))

sd_ <- tryCatch(dirname(sys.frame(1)$ofile), error = function(e) {
  ca <- commandArgs(FALSE); fa <- sub("^--file=", "", ca[grep("^--file=", ca)])
  if (length(fa)) dirname(fa) else getwd() })
source(file.path(sd_, "..", "..", "config.R"))
source(file.path(MODEL_DIR, "tDDM", "tddm_engine.R"))

# Orientation of the two signals. Default "mean" reproduces the canonical fit
# exactly; "contrast" runs the alternative ordering reported as a robustness check.
ORIENT <- match.arg(ga("--orientation", "mean"), c("mean", "contrast"))
S1 <- if (ORIENT == "mean") "mean_sig" else "contrast_sig"
S2 <- if (ORIENT == "mean") "contrast_sig" else "mean_sig"
LL <- mk_ll_tau(S1, S2, trans = get_transition("image"))
LO <- c(-3, -3, 0.5, 0.01, -2); UP <- c(3, 3, 4, 3, 3)
ctrl <- DEoptim.control(itermax = ITERMAX, NP = NP, trace = FALSE)   # NEVER set VTR

data <- tddm_load_study4(DATA_DIR)
subjects <- unique(data$subject)
cat(sprintf("=== optimizer reproducibility of the canonical 4j fit ===\n"))
cat(sprintf("  %d subjects | orientation = %s-first | reps = %d | restarts = %d | DE %d/%d | cores = %d\n\n",
            length(subjects), ORIENT, N_REP, RESTARTS, ITERMAX, NP, NC))

fit_one <- function(s) {
  td <- data[data$subject == s, ]
  if (nrow(td) < 5) return(NULL)
  best <- NULL
  for (k in seq_len(RESTARTS)) {
    f <- tryCatch(DEoptim(LL, LO, UP, ctrl, td = td), error = function(e) NULL)
    if (!is.null(f) && (is.null(best) || f$optim$bestval < best$optim$bestval)) best <- f
  }
  if (is.null(best)) return(NULL)
  e <- best$optim$bestmem
  data.frame(subject = s, group = td$group[1], n_trials = nrow(td),
             negLogLik = best$optim$bestval, de_iter = best$optim$iter,
             # e[1] is the LEADING signal weight and e[2] the joining one, so the
             # labels follow the orientation rather than the position.
             d_mean     = if (ORIENT == "mean") e[1] else e[2],
             d_contrast = if (ORIENT == "mean") e[2] else e[1],
             thres = e[3], nDT = e[4], tau = e[5],
             stringsAsFactors = FALSE)
}

ts <- format(Sys.time(), "%Y%m%d_%H%M%S")
od <- file.path(RESULTS_DIR, "tDDM", sprintf("optimizer_reproducibility_%s_r%d_%s",
                                             ORIENT, RESTARTS, ts))
dir.create(od, recursive = TRUE, showWarnings = FALSE)
RNGkind("L'Ecuyer-CMRG")

t0 <- Sys.time(); S <- list(); TAU <- NULL
for (r in seq_len(N_REP)) {
  set.seed(SEED0 + r)
  res <- do.call(rbind, Filter(is.data.frame, mclapply(subjects, fit_one, mc.cores = NC)))
  stopifnot(min(res$de_iter) > 0)                                     # VTR guard
  write.csv(res, file.path(od, sprintf("fits_rep%02d.csv", r)), row.names = FALSE)
  hc <- res$tau[res$group == "HC"]; dp <- res$tau[res$group == "DEP"]
  tt <- t.test(dp, hc)
  gd <- function(p) mean(res[[p]][res$group == "DEP"]) - mean(res[[p]][res$group == "HC"])
  S[[r]] <- data.frame(rep = r, seed = SEED0 + r, dtau = mean(dp) - mean(hc), p_dtau = tt$p.value,
                       tau_HC = mean(hc), tau_DEP = mean(dp), mean_NLL = mean(res$negLogLik),
                       d_dmean = gd("d_mean"), d_thres = gd("thres"), d_nDT = gd("nDT"),
                       d_dcontrast = gd("d_contrast"))
  TAU <- if (is.null(TAU)) res[, c("subject", "tau")] else merge(TAU, res[, c("subject", "tau")], by = "subject")
  cat(sprintf("  rep %2d/%d | Dtau %+.3f (p=%.4f) | tau HC %.3f DEP %.3f | mean NLL %.4f | %.1f min\n",
              r, N_REP, S[[r]]$dtau, S[[r]]$p_dtau, S[[r]]$tau_HC, S[[r]]$tau_DEP,
              S[[r]]$mean_NLL, as.numeric(difftime(Sys.time(), t0, units = "mins"))))
  flush.console()
}
S <- do.call(rbind, S); write.csv(S, file.path(od, "per_rep.csv"), row.names = FALSE)
names(TAU)[-1] <- paste0("tau_rep", seq_len(ncol(TAU) - 1))
write.csv(TAU, file.path(od, "tau_by_rep.csv"), row.names = FALSE)

cat("\n=== ACROSS-RUN VARIABILITY (same data, same settings, different RNG) ===\n")
cat(sprintf("  group Dtau : mean %+.3f | sd %.3f | range [%+.3f, %+.3f]\n",
            mean(S$dtau), sd(S$dtau), min(S$dtau), max(S$dtau)))
cat(sprintf("  p(Dtau)    : min %.4f | max %.4f | %d of %d below .05\n",
            min(S$p_dtau), max(S$p_dtau), sum(S$p_dtau < .05), nrow(S)))
cat(sprintf("  mean NLL   : sd %.5f  (a tie across runs means this is not optimizer quality)\n", sd(S$mean_NLL)))
cat(sprintf("\n  other group differences, sd across runs: d_mean %.4f | thres %.4f | nDT %.4f | d_contrast %.4f\n",
            sd(S$d_dmean), sd(S$d_thres), sd(S$d_nDT), sd(S$d_dcontrast)))
M <- as.matrix(TAU[, -1]); cr <- cor(M)
cat(sprintf("\n  pairwise r between runs for individual tau: mean %.3f | min %.3f\n",
            mean(cr[upper.tri(cr)]), min(cr[upper.tri(cr)])))
cat(sprintf("  within-subject sd of tau across runs: median %.3f\n", median(apply(M, 1, sd))))
cat(sprintf("\nSaved: %s\nTotal: %.1f min\n", od, as.numeric(difftime(Sys.time(), t0, units = "mins"))))
