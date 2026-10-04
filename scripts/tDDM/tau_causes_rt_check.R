#!/usr/bin/env Rscript
# =============================================================================
# Is median RT DOWNSTREAM of tau?  (the premise of the "over-control" reading)
# =============================================================================
# Controlling tau's group effect for median RT drops it from beta = +0.325
# (p = .004) to +0.174 (p = .110).  Two readings:
#
#   (a) medRT is a CONFOUND -> the control is appropriate and the onset effect
#       does not survive it.
#   (b) medRT is partly a CONSEQUENCE of tau (a later onset means the second
#       signal arrives later, so decisions take longer) -> conditioning on it
#       removes real signal, i.e. over-control.
#
# (b) is only defensible if tau ACTUALLY moves RT.  That is a pure forward
# question: simulate the same subjects at different tau with EVERYTHING else
# held fixed and measure d(medRT)/d(tau).  No fitting is involved, so this costs
# seconds, not hours.
#
# The number that matters is not just "is the path non-zero" but "how much of
# the OBSERVED r(tau_hat, medRT) = +0.337 can the causal path account for?"  The
# rest would be estimator artefact (tau_hat reading RT geometry), which the
# over-control argument cannot excuse.
#
# Usage: Rscript scripts/tDDM/tau_causes_rt_check.R [--reps 5] [--cores 4]
# =============================================================================
suppressPackageStartupMessages({ library(parallel) })
args <- commandArgs(trailingOnly = TRUE)
ga <- function(f, d) { i <- which(args == f); if (length(i) == 1 && i < length(args)) as.numeric(args[i + 1]) else d }
N_REP <- ga("--reps", 5)
NC    <- ga("--cores", 4)

sd_ <- tryCatch(dirname(sys.frame(1)$ofile), error = function(e) {
  ca <- commandArgs(FALSE); fa <- sub("^--file=", "", ca[grep("^--file=", ca)])
  if (length(fa)) dirname(fa) else getwd() })
source(file.path(sd_, "..", "..", "config.R"))
source(file.path(MODEL_DIR, "tDDM", "tddm_engine.R"))

# forward simulator, mean-first, UNRESCALED space (matches the fitting convention)
# 2026-09-28: simulator replaced by the shared bridge-corrected tddm_simulate_block()
# (models/tDDM/tddm_engine.R). The old Euler copy lengthened mean DT by 17-39 ms.
sim1 <- function(d1, d2, th, nDT, tau, s1, s2, sg = TDDM_SIGMA, dt = .001, mx = 10)
  tddm_simulate_trial(d1, d2, th, nDT, tau, s1, s2, sg = sg, dt = dt, mx = mx)$rt

dat  <- tddm_load_study4(DATA_DIR)
fits <- read.csv(file.path(PROJECT_ROOT, "results/tDDM/canonical_4j_accurate_bestof10/fits_4j_accurate_best10.csv"),   # 2026-09-28
                 stringsAsFactors = FALSE)
subj <- list()
for (i in seq_len(nrow(fits))) {
  tr <- dat[dat$subject == fits$subject[i], c("mean_sig", "contrast_sig")]
  tr <- tr[complete.cases(tr), ]
  if (nrow(tr) < 5) next
  subj[[length(subj) + 1]] <- list(d1 = fits$d_mean[i], d2 = fits$d_contrast[i],
                                   th = fits$thres[i], nDT = fits$nDT[i],
                                   grp = fits$group[i], tau = fits$tau[i],
                                   s1 = tr$mean_sig, s2 = tr$contrast_sig)
}
cat(sprintf("=== Does tau move median RT?  (forward simulation, no fitting) ===\n"))
cat(sprintf("  %d subjects, %d reps, tau held at each subject's fitted value + offset\n\n",
            length(subj), N_REP))

# each subject simulated at its own fitted tau PLUS a set of offsets; everything
# else (drift weights, threshold, nDT, that subject's own signals) held fixed
OFFSETS <- c(0, 0.25, 0.50, 1.00)
med_rt <- function(s, off) {
  n <- length(s$s1)
  median(vapply(seq_len(n), function(j)
    sim1(s$d1, s$d2, s$th, s$nDT, max(s$tau + off, 0), s$s1[j], s$s2[j]), 0))
}

set.seed(20260811)
rows <- list()
for (rep in 1:N_REP) {
  for (off in OFFSETS) {
    mr <- unlist(mclapply(subj, function(s) med_rt(s, off), mc.cores = NC))
    rows[[length(rows) + 1]] <- data.frame(rep = rep, offset = off,
                                           med_rt = mean(mr, na.rm = TRUE))
  }
  cat(sprintf("  rep %d/%d done\n", rep, N_REP)); flush.console()
}
res <- do.call(rbind, rows)

cat("\n--- mean of per-subject median RT, by tau offset ---\n")
# aggregate() with a multi-value FUN returns a matrix column; index it directly
mm <- tapply(res$med_rt, res$offset, mean)
ss <- tapply(res$med_rt, res$offset, sd)
for (o in names(mm))
  cat(sprintf("  tau + %.2f s  ->  median RT %.4f s  (sd over reps %.4f)\n",
              as.numeric(o), mm[[o]], ss[[o]]))

fit <- lm(med_rt ~ offset, data = res)
b  <- coef(fit)[2]; se <- summary(fit)$coefficients[2, 2]
cat(sprintf("\n  d(median RT)/d(tau) = %+.4f  [95%% CI %+.4f, %+.4f]\n",
            b, b - 1.96 * se, b + 1.96 * se))

# how much of the observed r(tau_hat, medRT) could this causal path produce?
sd_tau <- sd(fits$tau)
dat2 <- dat[dat$subject %in% fits$subject, ]
sd_medrt <- sd(tapply(dat2$RT_ddm, dat2$subject, median))
r_implied <- b * sd_tau / sd_medrt
cat(sprintf("\n  SD(tau_hat) = %.3f s, SD(median RT across subjects) = %.3f s\n", sd_tau, sd_medrt))
cat(sprintf("  => causal path alone implies r(tau, medRT) = %+.3f\n", r_implied))
cat(sprintf("     OBSERVED r(tau_hat, medRT) = +0.337\n"))
cat(sprintf("     causal share = %.0f%% of the observed correlation\n", 100 * r_implied / 0.337))
cat("\n  If this share is small, most of the tau-RT link is the ESTIMATOR reading RT\n")
cat("  geometry, and the over-control argument does not rescue the group effect.\n")

ts <- format(Sys.time(), "%Y%m%d_%H%M%S")
od <- file.path(RESULTS_DIR, "tDDM", sprintf("tau_causes_rt_%s", ts))
dir.create(od, recursive = TRUE, showWarnings = FALSE)
write.csv(res, file.path(od, "per_rep.csv"), row.names = FALSE)
cat(sprintf("\nOutput: %s\n", od))
