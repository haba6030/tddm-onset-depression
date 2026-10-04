#!/usr/bin/env Rscript
# =============================================================================
# Does 4j REPRODUCE the observed RT signature of mean and contrast?
# =============================================================================
# Lab-meeting item #4 established, model-free, that the two axes act differently
# on RT (mixed model on the raw data, scripts/tDDM/labfeedback_2026_08.py):
#
#   mean      -> SIGNED effect on log RT     (higher common value, faster)  = DRIFT
#   contrast  -> no signed effect, but a positive |contrast| effect (slower) = CONFLICT
#
# That regression touches no tDDM machinery, so the transition-density fix does
# not invalidate it.  What the fix DOES change is whether the model can be asked
# the question at all: under the shipped density only 44% of trials ever crossed
# tau, so the contrast channel was inert in most of the likelihood.  Corrected,
# 95% of trials enter phase 2 -- the contrast term is finally doing work, and its
# RT signature becomes a testable prediction.
#
# The test: simulate each participant's own trials from their own fitted 4j
# parameters, run the SAME mixed model on the simulated RTs, and compare the four
# coefficients to the empirical ones.
#
# The prediction that matters:  4j enters contrast as a SIGNED drift term
# (v = d_mean*mean + d_contrast*contrast after tau).  In a diffusion process RT
# falls with |v|, so a larger |contrast| should on average make responses FASTER.
# The data say the opposite (slower).  If the simulation confirms that sign flip,
# 4j misspecifies the contrast channel -- a model-adequacy finding, not a nuisance.
#
# Usage: Rscript scripts/tDDM/rt_signature_model_check.R [--reps 20] [--cores 8]
# =============================================================================
suppressPackageStartupMessages({ library(parallel); library(lme4) })

args <- commandArgs(trailingOnly = TRUE)
ga <- function(f, d) { i <- which(args == f)
  if (length(i) == 1 && i < length(args)) as.numeric(args[i + 1]) else d }
N_REP <- ga("--reps", 20)
NC    <- ga("--cores", 8)

sd_ <- tryCatch(dirname(sys.frame(1)$ofile), error = function(e) {
  ca <- commandArgs(FALSE); fa <- sub("^--file=", "", ca[grep("^--file=", ca)])
  if (length(fa)) dirname(fa) else getwd() })
source(file.path(sd_, "..", "..", "config.R"))
source(file.path(MODEL_DIR, "tDDM", "tddm_engine.R"))

# 2026-09-26: canonical fits are the highest-likelihood of ten optimizer runs.
# Pass --fits to point at another fits file.
ga2 <- function(f, d) { if (f %in% args) { i <- which(args == f); if (i < length(args)) return(args[i + 1]) }; d }
FITS <- file.path(PROJECT_ROOT, ga2("--fits",
  "results/tDDM/canonical_4j_accurate_bestof10/fits_4j_accurate_best10.csv"))

# ---- vectorised two-phase simulator ----------------------------------------
# One path per trial, all trials stepped together.  Mean-first: drift is
# d_mean*mean before tau, plus d_contrast*contrast after.  Bias fixed at a/2
# (the 5-parameter 4j has no free bias).
# 2026-09-28: simulator replaced by the shared bridge-corrected tddm_simulate_block()
# (models/tDDM/tddm_engine.R). The old Euler copy lengthened mean DT by 17-39 ms.
sim_block <- function(dm, dc, th, ndt, tau, ms, cs, sg = TDDM_SIGMA, dt = 0.001, tmax = 10)
  tddm_simulate_block(dm, dc, th, ndt, tau, ms, cs, sg = sg, dt = dt, tmax = tmax)$RT

# ---- data + fits ------------------------------------------------------------
dat  <- tddm_load_study4(DATA_DIR)
fits <- read.csv(FITS, stringsAsFactors = FALSE)

subj <- list()
for (i in seq_len(nrow(fits))) {
  tr <- dat[dat$subject == fits$subject[i], c("mean_sig", "contrast_sig", "RT_ddm")]
  tr <- tr[complete.cases(tr), ]
  if (nrow(tr) < 5) next
  subj[[length(subj) + 1]] <- list(
    sub = fits$subject[i], grp = fits$group[i],
    dm = fits$d_mean[i], dc = fits$d_contrast[i], th = fits$thres[i],
    ndt = fits$nDT[i], tau = fits$tau[i],
    ms = tr$mean_sig, cs = tr$contrast_sig, rt = tr$RT_ddm)
}
cat(sprintf("=== RT signature: does 4j reproduce it? ===\n"))
cat(sprintf("  %d subjects, %d trials, %d reps\n\n",
            length(subj), sum(vapply(subj, function(s) length(s$ms), 0L)), N_REP))

TERMS <- c("mean_sig", "abs_mean", "contrast_sig", "abs_contrast")
run_lmm <- function(df) {
  df$logRT <- log(df$RT)
  df$abs_mean <- abs(df$mean_sig); df$abs_contrast <- abs(df$contrast_sig)
  m <- lmer(logRT ~ mean_sig + abs_mean + contrast_sig + abs_contrast + (1 | sub),
            data = df, REML = FALSE,
            control = lmerControl(calc.derivs = FALSE))
  setNames(fixef(m)[TERMS], TERMS)
}

# ---- observed ---------------------------------------------------------------
obs <- do.call(rbind, lapply(subj, function(s)
  data.frame(sub = s$sub, mean_sig = s$ms, contrast_sig = s$cs, RT = s$rt)))
b_obs <- run_lmm(obs)

cat("--- OBSERVED (raw behaviour, model-free) ---\n")
for (k in TERMS) cat(sprintf("  %-14s %+.4f\n", k, b_obs[[k]]))

# ---- simulated --------------------------------------------------------------
set.seed(20260814)
seeds <- sample.int(1e6, N_REP)
sim_one <- function(sd_i) {
  set.seed(sd_i)
  df <- do.call(rbind, lapply(subj, function(s)
    data.frame(sub = s$sub, mean_sig = s$ms, contrast_sig = s$cs,
               RT = sim_block(s$dm, s$dc, s$th, s$ndt, s$tau, s$ms, s$cs))))
  df <- df[df$RT > 0.15 & df$RT <= 10, ]
  run_lmm(df)
}
res <- mclapply(seeds, sim_one, mc.cores = NC)
res <- do.call(rbind, res[!vapply(res, is.null, TRUE)])

cat(sprintf("\n--- SIMULATED from 4j fits (%d reps) ---\n", nrow(res)))
cat("  term            model beta      95%% CI over reps     observed    match?\n")
verdict <- list()
for (k in TERMS) {
  v  <- res[, k]
  m  <- mean(v); lo <- quantile(v, .025); hi <- quantile(v, .975)
  ok <- if (b_obs[[k]] >= lo && b_obs[[k]] <= hi) "yes" else
        if (sign(m) != sign(b_obs[[k]]))          "SIGN FLIP" else "outside"
  verdict[[k]] <- ok
  cat(sprintf("  %-14s %+.4f      [%+.4f, %+.4f]    %+.4f     %s\n",
              k, m, lo, hi, b_obs[[k]], ok))
}

cat("\n--- verdict ---\n")
if (identical(verdict$abs_contrast, "SIGN FLIP"))
  cat("  |contrast|: model and data DISAGREE IN SIGN.\n",
      "  4j enters contrast as a signed drift, so larger |contrast| speeds\n",
      "  responses; the data slow down. The contrast channel is misspecified\n",
      "  as drift -- it behaves like a CONFLICT term.\n") else
  cat(sprintf("  |contrast|: %s\n", verdict$abs_contrast))

ts <- format(Sys.time(), "%Y%m%d_%H%M%S")
od <- file.path(RESULTS_DIR, "tDDM", sprintf("rt_signature_check_%s", ts))
dir.create(od, recursive = TRUE, showWarnings = FALSE)
write.csv(res, file.path(od, "sim_coefficients.csv"), row.names = FALSE)
write.csv(data.frame(term = TERMS, observed = as.numeric(b_obs[TERMS]),
                     model_mean = colMeans(res[, TERMS]),
                     verdict = unlist(verdict[TERMS])),
          file.path(od, "comparison.csv"), row.names = FALSE)

# quintile curves for the overlay figure
qbin <- function(df, col) {
  q <- cut(df[[col]], breaks = quantile(df[[col]], seq(0, 1, .2)),
           include.lowest = TRUE, labels = FALSE)
  data.frame(x = tapply(df[[col]], q, mean), rt = tapply(df$RT, q, median))
}
set.seed(1)
simdf <- do.call(rbind, lapply(subj, function(s)
  data.frame(sub = s$sub, mean_sig = s$ms, contrast_sig = s$cs,
             RT = sim_block(s$dm, s$dc, s$th, s$ndt, s$tau, s$ms, s$cs))))
simdf$abs_contrast <- abs(simdf$contrast_sig)
obs$abs_contrast <- abs(obs$contrast_sig)

# 2026-08-20: trial-level join for fig4_rt_signature.py (observed vs simulated RT
# on the same trials). obs and simdf are built from the same subject loop in the
# same order, so they align row by row.
stopifnot(nrow(obs) == nrow(simdf), all(obs$sub == simdf$sub))
grp <- vapply(subj, function(s) s$grp, "")[match(obs$sub, vapply(subj, function(s) s$sub, obs$sub[1]))]
write.csv(data.frame(sub = obs$sub, group = grp,
                     mean_sig = obs$mean_sig, contrast_sig = obs$contrast_sig,
                     rt_obs = obs$RT, rt_sim = simdf$RT),
          file.path(od, "trial_level.csv"), row.names = FALSE)
for (v in c("mean_sig", "contrast_sig", "abs_contrast")) {
  o <- qbin(obs, v); m <- qbin(simdf, v)
  write.csv(data.frame(x = o$x, rt_obs = o$rt, rt_model = m$rt),
            file.path(od, sprintf("quintile_%s.csv", v)), row.names = FALSE)
}
cat(sprintf("\nOutput: %s\n", od))
