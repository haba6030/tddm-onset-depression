#!/usr/bin/env Rscript
# =============================================================================
# NULL-tau generative control, re-run under the CORRECTED transition density
# =============================================================================
# The 2026-06-28 version of this test (scripts/tDDM/nulltau_generative_sim.R)
# concluded that the estimator MANUFACTURES the observed group onset gap: with a
# single common tau for every subject (true Dtau = 0) it still recovered
# Dtau_hat = -0.193, i.e. 114% of the observed -0.169, in 30/30 reps.
#
# That test was run with the SHIPPED transition density, which we have since
# shown to be misspecified (mass +70% at tau = 2 s, +133% at a low threshold;
# models/tDDM/tddm_engine.R).  Since the bias it induces grows with tau and
# shrinks with the threshold -- exactly the quantities that differ between
# groups -- the manufactured gap may be an artefact OF THE ARTEFACT TEST.
#
# This re-run changes ONE thing: the fitting likelihood uses the method-of-images
# transition density.  Everything else (generating parameters, common tau, DE
# settings, orientation) is held at the 2026-06-28 values so the comparison is
# clean.
#
#   recovered Dtau_hat ~ -0.19  -> artefact confirmed even under a correct likelihood
#   recovered Dtau_hat ~ 0      -> the artefact was a symptom of the bad density,
#                                  and the observed group gap survives
#
# Usage:
#   Rscript scripts/tDDM/nulltau_corrected.R --reps 8 --cores 4
#   Rscript scripts/tDDM/nulltau_corrected.R --trans pipeline   # reproduce the old result
# =============================================================================
suppressPackageStartupMessages({ library(DEoptim); library(parallel) })
args <- commandArgs(trailingOnly = TRUE)
ga <- function(f, d) { if (f %in% args) { i <- which(args == f); if (i < length(args)) return(args[i+1]) }; d }
N_REP   <- as.integer(ga("--reps", 8))
DE_ITER <- as.integer(ga("--de-iter", 60))    # 2026-06-28 setting, kept for comparability
DE_NP   <- as.integer(ga("--de-np", 40))
NC      <- as.integer(ga("--cores", 4))
TRANS   <- match.arg(ga("--trans", "image"), c("image", "pipeline"))
# 2026-08-19 additions
#   --fits  generating parameters.  The default reproduces the 2026-06-28 design
#           (shipped-density fits) so the published -0.150 stays reproducible;
#           pass the corrected fits when the question is about the CURRENT model.
#   --wgap  group difference in the STARTING POINT (w_HC = .5 + gap/2,
#           w_DEP = .5 - gap/2), generated but NOT fitted (the likelihood fixes
#           w = 0.5).  This is the second negative control: how much onset gap is
#           manufactured by an unmodelled starting-point difference?  The data
#           bound such a difference at |dw| <= 0.016 (the whole HC-DEP accept-rate
#           difference, -0.0159, p = .53, attributed to the starting point at zero
#           drift), so 0.05 and 0.10 are 3x and 6x the empirical ceiling.
FITS_ARG <- ga("--fits", "")
WGAP     <- as.numeric(ga("--wgap", 0))
# 2026-09-23: match the observed procedure, which now keeps the BEST-NLL run of R
# repeats of the tau fit (docs/PAPER/Memo_Analyses_0919.md S4, TODO-1).
RESTARTS <- as.integer(ga("--restarts", 1))

sd_ <- tryCatch(dirname(sys.frame(1)$ofile), error = function(e) {
  ca <- commandArgs(FALSE); fa <- sub("^--file=", "", ca[grep("^--file=", ca)])
  if (length(fa)) dirname(fa) else getwd() })
source(file.path(sd_, "..", "..", "config.R"))
source(file.path(MODEL_DIR, "tDDM", "tddm_engine.R"))
TDDM_NGRID <- as.integer(ga("--ng", 25))      # 2026-06-28 setting

TRANS_FUN <- get_transition(TRANS)
# MEAN-first (the manuscript model): mean leads, contrast joins at tau
LL <- mk_ll_tau("s1", "s2", trans = TRANS_FUN)   # s1 = mean, s2 = contrast (named below)

# 2026-09-27: body replaced by the shared bridge-corrected simulator (models/tDDM/tddm_engine.R,
# tddm_simulate_trial). The old Euler copy lengthened mean DT by 17-39 ms vs the likelihood.
simulate_ddm <- function(d1, d2, th, nDT, tau, s1, s2, w = 0.5, sg = TDDM_SIGMA, dt = .001, mx = 10)
  tddm_simulate_trial(d1, d2, th, nDT, tau, s1, s2, w = w, sg = sg, dt = dt, mx = mx)

cat("=== NULL-tau generative control, CORRECTED density ===\n")
cat(sprintf("  transition = %s | reps = %d | DE itermax=%d NP=%d ngrid=%d | cores = %d\n",
            TRANS, N_REP, DE_ITER, DE_NP, TDDM_NGRID, NC))

dat <- tddm_load_study4(DATA_DIR)
# generating parameters: the same mean-first fits the 2026-06-28 test used, so the
# only thing that changes between the two runs is the FITTING likelihood
FITS_PATH <- if (nzchar(FITS_ARG)) FITS_ARG else
  file.path(RESULTS_DIR, "tDDM", "canonical_4j_accurate_bestof10", "fits_4j_accurate_best10.csv")
fits <- read.csv(FITS_PATH, stringsAsFactors = FALSE)
cat(sprintf("  generating parameters: %s\n", FITS_PATH))
if (WGAP != 0) cat(sprintf("  starting-point gap: w_HC = %.3f, w_DEP = %.3f (fitted with w fixed at 0.5)\n",
                           0.5 + WGAP/2, 0.5 - WGAP/2))
# The common tau must sit at the OPERATING POINT of the likelihood being tested.
# The shipped-density fits give tau ~ 1.80, but the corrected density puts the
# real-data tau at ~0.72 -- and the induced bias depends on how many trials
# actually cross into phase 2, which is a function of tau.  Testing at 1.80 while
# the corrected model operates at 0.72 tests the wrong regime.
TAU_COMMON <- as.numeric(ga("--tau-common", NA_real_))
if (!is.finite(TAU_COMMON)) {
  if (TRANS == "image")
    stop("--tau-common is REQUIRED with --trans image.\n",
         "  The generating tau must sit at the OPERATING POINT of the density under test.\n",
         "  The fits file read here was produced with the SHIPPED density (mean tau ", sprintf("%.3f", mean(fits$tau)), "),\n",
         "  but the corrected density puts the real-data tau near 0.72. Testing at 1.80 while the\n",
         "  model operates at 0.72 tests the wrong regime and flipped this conclusion once already.\n",
         "  Pass e.g. --tau-common 0.720 (mean of results/tDDM/corrected_likelihood_*/fits_4j_image.csv).")
  TAU_COMMON <- mean(fits$tau)      # shipped-density default, only for --trans pipeline
}
obs_gap <- mean(fits$tau[fits$group == "HC"]) - mean(fits$tau[fits$group == "DEP"])
cat(sprintf("  observed Dtau (HC-DEP) = %+.3f | common tau assigned to all = %.3f (true Dtau = 0)\n\n",
            obs_gap, TAU_COMMON))

subj <- list()
for (i in seq_len(nrow(fits))) {
  sid <- fits$subject[i]
  tr <- dat[dat$subject == sid, c("mean_sig", "contrast_sig")]
  tr <- tr[complete.cases(tr), ]
  if (nrow(tr) < 5) next
  subj[[length(subj) + 1]] <- list(d1 = fits$d_mean[i], d2 = fits$d_contrast[i],
                                   th = fits$thres[i], nDT = fits$nDT[i],
                                   grp = fits$group[i], s1 = tr$mean_sig, s2 = tr$contrast_sig,
                                   w = if (fits$group[i] == "HC") 0.5 + WGAP/2 else 0.5 - WGAP/2)
}
NSUB <- as.integer(ga("--nsub", length(subj)))
if (NSUB < length(subj)) subj <- subj[sample(seq_along(subj), NSUB)]
cat(sprintf("  subjects = %d\n\n", length(subj)))

ctrl <- DEoptim.control(itermax = DE_ITER, NP = DE_NP, trace = FALSE)   # no VTR
# 2026-08-19: return ALL five estimates, not just tau_hat.  The generator already
# carries the real per-group differences in drift / threshold / nDT, so the same
# run answers "does the estimator manufacture a group difference in THOSE too?" --
# previously only tau was checked, which left three of the four reported group
# effects without a generative control.
fit1 <- function(s, tau_true) {
  n <- length(s$s1); rt <- numeric(n); ch <- integer(n)
  for (j in 1:n) { r <- simulate_ddm(s$d1, s$d2, s$th, s$nDT, tau_true, s$s1[j], s$s2[j], w = s$w)
                   rt[j] <- r$rt; ch[j] <- r$ch }
  td <- data.frame(s1 = s$s1, s2 = s$s2, RT_ddm = rt, choice_ddm = ch)
  td <- td[td$RT_ddm > 0.15 & td$RT_ddm <= 10, ]
  if (nrow(td) < 5) return(NULL)
  f <- NULL
  for (k in seq_len(RESTARTS)) {
    fk <- tryCatch(DEoptim(LL, c(-3, -3, 0.5, 0.01, -2), c(3, 3, 4, 3, 3), ctrl, td = td),
                   error = function(e) NULL)
    if (!is.null(fk) && (is.null(f) || fk$optim$bestval < f$optim$bestval)) f <- fk
  }
  if (is.null(f)) return(NULL)
  e <- f$optim$bestmem
  # 2026-09-28: every tau <= 0 gives the contrast model (identical likelihood), so tau is
  # clamped at 0 before any summary, matching the manuscript's canonical fits.
  data.frame(tau_hat = max(e[5], 0), tau_hat_raw = e[5], d_mean_hat = e[1], d_contrast_hat = e[2],
             thres_hat = e[3], nDT_hat = e[4],
             d_mean_true = s$d1, d_contrast_true = s$d2,
             thres_true = s$th, nDT_true = s$nDT,
             grp = s$grp, de_iter = f$optim$iter)
}

RNGkind("L'Ecuyer-CMRG")   # 2026-09-27: forked workers inherit reproducible streams
set.seed(20260807)
ts <- format(Sys.time(), "%Y%m%d_%H%M%S")
od <- file.path(RESULTS_DIR, "tDDM", sprintf("nulltau_corrected_%s_%s", TRANS, ts))
dir.create(od, recursive = TRUE, showWarnings = FALSE)
pf <- file.path(od, "per_rep.csv"); t0 <- Sys.time(); rows <- list()
for (rep in 1:N_REP) {
  # 2026-09-28: mclapply seeds each worker from the parent's L'Ecuyer stream WITHOUT advancing it,
  # so every call needs its own set.seed(); otherwise all reps reproduce the same simulated data.
  set.seed(20260807 + rep)
  r <- mclapply(subj, function(s) fit1(s, TAU_COMMON), mc.cores = NC, mc.preschedule = FALSE)
  r <- do.call(rbind, Filter(Negate(is.null), r))
  if (is.null(r) || nrow(r) < 10) { cat(sprintf("  rep %d SKIP\n", rep)); next }
  stopifnot(min(r$de_iter) > 0)                                  # VTR guard
  gd <- function(v) mean(v[r$grp == "HC"]) - mean(v[r$grp == "DEP"])   # HC - DEP
  dh <- gd(r$tau_hat)
  row <- data.frame(rep = rep, recovered_dtau = dh,
                    tauhat_HC = mean(r$tau_hat[r$grp == "HC"]),
                    tauhat_DEP = mean(r$tau_hat[r$grp == "DEP"]), n = nrow(r),
                    # recovered vs TRUE group difference, per parameter
                    d_mean_rec     = gd(r$d_mean_hat),     d_mean_true     = gd(r$d_mean_true),
                    d_contrast_rec = gd(r$d_contrast_hat), d_contrast_true = gd(r$d_contrast_true),
                    thres_rec      = gd(r$thres_hat),      thres_true      = gd(r$thres_true),
                    nDT_rec        = gd(r$nDT_hat),        nDT_true        = gd(r$nDT_true),
                    wgap = WGAP)
  write.csv(r, file.path(od, sprintf("subject_level_rep%02d.csv", rep)), row.names = FALSE)
  rows[[length(rows) + 1]] <- row
  h <- !file.exists(pf); write.table(row, pf, sep = ",", row.names = FALSE, col.names = h, append = !h)
  cat(sprintf("  rep %d/%d: recovered Dtau_hat = %+.3f  (true 0; observed %+.3f)  [%.1f min]\n",
              rep, N_REP, dh, obs_gap, as.numeric(difftime(Sys.time(), t0, units = "mins"))))
  flush.console()
}
res <- do.call(rbind, rows)
m <- mean(res$recovered_dtau); s <- sd(res$recovered_dtau); n <- nrow(res)
ci <- m + c(-1.96, 1.96) * s / sqrt(n)
cat("\n============================================================\n")
cat(sprintf("  transition density : %s\n", TRANS))
cat(sprintf("  mean recovered Dtau_hat = %+.3f  (sd %.3f, 95%% CI [%+.3f, %+.3f], n_rep = %d)\n", m, s, ci[1], ci[2], n))
cat(sprintf("  observed Dtau           = %+.3f\n", obs_gap))
cat(sprintf("  manufactured share      = %.0f%% of the observed gap\n", 100 * m / obs_gap))
cat(sprintf("  reps reproducing the observed sign: %d/%d\n", sum(sign(res$recovered_dtau) == sign(obs_gap)), n))

# ---- group-difference recovery for EVERY parameter (2026-08-19) -------------
cat("\n  --- group-difference recovery, all parameters (HC - DEP) ---\n")
cat(sprintf("  %-12s %10s %10s %10s\n", "parameter", "TRUE", "recovered", "bias"))
for (pn in c("d_mean", "d_contrast", "thres", "nDT")) {
  tv <- mean(res[[paste0(pn, "_true")]]); rv <- mean(res[[paste0(pn, "_rec")]])
  cat(sprintf("  %-12s %+10.3f %+10.3f %+10.3f\n", pn, tv, rv, rv - tv))
}
cat(sprintf("  %-12s %+10.3f %+10.3f %+10.3f   <- true 0 by construction\n", "tau", 0, m, m))
cat("\n  2026-06-28 result with the SHIPPED density: -0.193, i.e. 114% of observed, 30/30 reps\n")
cat("============================================================\n")
cat(sprintf("\nOutput: %s\n", od))
