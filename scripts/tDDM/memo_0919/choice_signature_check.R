#!/usr/bin/env Rscript
# =============================================================================
# Does the time-contrast model reproduce the model-free group x contrast
# interaction on ACCEPTANCE, and through which parameter?
# =============================================================================
# Model-free (Results :152, memo_analyses.R m1): the contrast signal lowers
# acceptance in HC (beta -0.81) and not in DEP (+0.26), interaction +1.07,
# p < .001. The fitted contrast WEIGHT does not differ between groups (p = .76).
#
# Test: simulate choice + RT for every participant's own trials from their own
# best-run 4j fits, fit the SAME logistic mixed model, and compare. Three arms:
#   fitted    : fitted parameters as estimated
#   tau_equal : every participant's tau set to the pooled median tau
#               (removes the group onset difference, keeps everything else)
#   dc_perm   : contrast weights permuted across participants, ignoring group
#               (removes any group difference in the weight, keeps its spread)
# If `fitted` reproduces the interaction and `tau_equal` loses it, the model
# carries the group difference in contrast sensitivity through the onset.
#
# Usage: Rscript scripts/tDDM/memo_0919/choice_signature_check.R [--reps 20] [--cores 8]
# =============================================================================
suppressPackageStartupMessages({ library(parallel); library(lme4) })
if (!file.exists("config.R")) stop("run from the repository root")
source("config.R"); source(file.path(MODEL_DIR, "tDDM", "tddm_engine.R"))

args <- commandArgs(trailingOnly = TRUE)
ga <- function(f, d) { i <- which(args == f)
  if (length(i) == 1 && i < length(args)) as.numeric(args[i + 1]) else d }
N_REP <- ga("--reps", 20); NC <- ga("--cores", 8)
FITS <- "results/tDDM/canonical_4j_accurate_bestof10/fits_4j_accurate_best10.csv"   # 2026-09-28: accurate-integral canonical

# 2026-09-28: simulator replaced by the shared bridge-corrected tddm_simulate_block()
# (models/tDDM/tddm_engine.R). The old Euler copy lengthened mean DT by 17-39 ms.
sim_block <- function(dm, dc, th, ndt, tau, ms, cs, sg = TDDM_SIGMA, dt = 0.001, tmax = 10)
  tddm_simulate_block(dm, dc, th, ndt, tau, ms, cs, sg = sg, dt = dt, tmax = tmax)

dat  <- tddm_load_study4(DATA_DIR)
fits <- read.csv(FITS, stringsAsFactors = FALSE)
dat  <- dat[dat$subject %in% fits$subject, ]
dat$DEP <- as.numeric(dat$group == "DEP"); dat$acc <- as.numeric(dat$choice == 0)

TERMS <- c("mean_sig", "contrast_sig", "mean_sig:DEP", "contrast_sig:DEP")
run_glmer <- function(df) {
  m <- glmer(acc ~ (mean_sig + contrast_sig) * DEP + (1 + mean_sig | subject), df,
             family = binomial, control = glmerControl(optimizer = "bobyqa"))
  b <- fixef(m)[TERMS]
  c(b, DEP_contrast_slope = unname(b["contrast_sig"] + b["contrast_sig:DEP"]))
}

obs <- run_glmer(dat)
cat("=== OBSERVED ===\n"); print(round(obs, 3))

tau_med <- median(fits$tau)
arms <- c("fitted", "tau_equal", "dc_perm")
set.seed(20260927); seeds <- sample.int(1e6, N_REP)
one <- function(seed, arm) {
  set.seed(seed)
  f <- fits
  if (arm == "tau_equal") f$tau <- tau_med
  if (arm == "dc_perm")   f$d_contrast <- sample(f$d_contrast)
  out <- do.call(rbind, lapply(seq_len(nrow(f)), function(i) {
    tr <- dat[dat$subject == f$subject[i], ]
    s  <- sim_block(f$d_mean[i], f$d_contrast[i], f$thres[i], f$nDT[i], f$tau[i],
                    tr$mean_sig, tr$contrast_sig)
    data.frame(subject = tr$subject, DEP = tr$DEP, mean_sig = tr$mean_sig,
               contrast_sig = tr$contrast_sig, acc = s$acc, RT = s$RT)
  }))
  out <- out[!is.na(out$acc) & out$RT > 0.15 & out$RT <= 10, ]
  tryCatch(run_glmer(out), error = function(e) NULL)
}

od <- file.path(RESULTS_DIR, "tDDM", sprintf("choice_signature_check_%s", format(Sys.time(), "%Y%m%d_%H%M%S")))
dir.create(od, recursive = TRUE, showWarnings = FALSE)
summ <- list()
for (arm in arms) {
  res <- mclapply(seeds, one, arm = arm, mc.cores = NC)
  res <- do.call(rbind, res[!vapply(res, is.null, TRUE)])
  write.csv(res, file.path(od, sprintf("sim_%s.csv", arm)), row.names = FALSE)
  cat(sprintf("\n=== SIMULATED: %s (%d reps) ===\n", arm, nrow(res)))
  for (k in colnames(res)) {
    v <- res[, k]
    cat(sprintf("  %-20s mean %+.3f  95%% [%+.3f, %+.3f]  observed %+.3f\n",
                k, mean(v), quantile(v, .025), quantile(v, .975), obs[[k]]))
    summ[[length(summ) + 1]] <- data.frame(arm = arm, term = k, sim_mean = mean(v),
      lo = quantile(v, .025), hi = quantile(v, .975), observed = obs[[k]])
  }
}
write.csv(do.call(rbind, summ), file.path(od, "summary.csv"), row.names = FALSE)
cat(sprintf("\nOutput: %s\n", od))
