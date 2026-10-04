#!/usr/bin/env Rscript
# Manuscript numbers from the stage-6 controls (accurate quadrature, C++ likelihood,
# bridge-corrected simulator, per-call seeds). Same formulas as the grid-era scripts
# (inversion_bestrun.R, null_correlations_overcontrol.py). Observed values come from the
# accurate canonical fits with tau <= 0 set to 0.
#   Rscript scripts/tDDM/memo_0919/stage6_summary_accurate.R
if (!file.exists("config.R")) stop("run from the repository root")
source("config.R"); source(file.path(MODEL_DIR, "tDDM", "tddm_engine.R"))
R <- file.path(RESULTS_DIR, "tDDM")
OD <- file.path(R, "canonical_4j_accurate_bestof10"); sink(file.path(OD, "stage6_summary.txt"), split = TRUE)
f <- read.csv(file.path(OD, "fits_4j_accurate_best10.csv")); f$tau <- pmax(f$tau, 0)
g <- f$group; tt <- t.test(f$tau[g == "DEP"], f$tau[g == "HC"])
OBS <- unname(diff(rev(tt$estimate))); OBS <- mean(f$tau[g == "DEP"]) - mean(f$tau[g == "HC"]); OBS_SE <- unname(tt$stderr)
cat(sprintf("observed dtau (DEP - HC, tau<=0 -> 0) = %.3f s, SE %.3f\n\n", OBS, OBS_SE))

cat("== positive control, 20 trials, 4 levels ==\n")
pp <- read.csv(file.path(R, "positive_control_tau_20260928_231637", "per_point.csv"))
m <- lm(recovered_gap ~ true_gap, pp); V <- vcov(m); ci <- confint(m)
b0 <- coef(m)[1]; b1 <- coef(m)[2]; s0 <- sqrt(V[1, 1]); s1 <- sqrt(V[2, 2]); rho <- V[1, 2] / (s0 * s1)
cat(sprintf("  %d points | sensitivity %.3f [%.3f, %.3f] (SE %.3f) | bias %.3f (SE %.3f)\n", nrow(pp), b1, ci[2, 1], ci[2, 2], s1, b0, s0))
set.seed(20260926); N <- 200000; z1 <- rnorm(N); z2 <- rnorm(N)
d0 <- b0 + s0 * z1; d1 <- b1 + s1 * (rho * z1 + sqrt(1 - rho^2) * z2)
tv <- (rnorm(N, OBS, OBS_SE) - d0) / d1; tv <- tv[is.finite(tv) & d1 > 0]
q <- quantile(tv, c(.025, .5, .975))
cat(sprintf("  latent dtau: median %+.3f s | 95%% [%+.3f, %+.3f] | P(>0) = %.3f\n", q[2], q[1], q[3], mean(tv > 0)))
cat("== positive control, 80 trials, 2 levels ==\n")
p8 <- read.csv(file.path(R, "positive_control_tau_20260929_012719", "per_point.csv"))
m8 <- lm(recovered_gap ~ true_gap, p8); c8 <- summary(m8)$coefficients; ci8 <- confint(m8)
cat(sprintf("  %d points | sensitivity %.3f [%.3f, %.3f] | bias %.3f (SE %.3f)\n\n", nrow(p8), c8[2, 1], ci8[2, 1], ci8[2, 2], c8[1, 1], c8[1, 2]))

cat("== null onset (true tau common) ==\n")
ND <- file.path(R, "nulltau_corrected_image_20260928_211818")
pr <- read.csv(file.path(ND, "per_rep.csv")); nul <- -pr$recovered_dtau          # file is HC - DEP
nm <- mean(nul); nse <- sd(nul) / sqrt(length(nul))
res <- OBS - nm; rse <- sqrt(OBS_SE^2 + nse^2); rp <- 2 * pnorm(-abs(res / rse))
cat(sprintf("  %d reps | recovered %.3f s (sd %.3f) = %.0f%% of observed | residual %.3f s, SE %.3f, p = %.2f\n",
            length(nul), nm, sd(nul), 100 * nm / OBS, res, rse, rp))
bias <- c(d_mean = mean(pr$d_mean_rec - pr$d_mean_true), thres = mean(pr$thres_rec - pr$thres_true),
          d_contrast = mean(pr$d_contrast_rec - pr$d_contrast_true))
cat(sprintf("  bias of group differences (HC - DEP, recovered - true): drift %+.3f, threshold %+.3f, contrast %+.3f\n",
            bias["d_mean"], bias["thres"], bias["d_contrast"]))
dat <- tddm_load_study4(DATA_DIR)
ces <- sapply(f$subject, function(s) dat$cesScore[dat$subject == s][1])
rs <- sapply(sort(Sys.glob(file.path(ND, "subject_level_rep*.csv"))), function(p) {
  s <- read.csv(p); stopifnot(isTRUE(all.equal(s$d_mean_true, f$d_mean)))     # rows follow the fits file
  c(thres = cor(s$tau_hat, s$thres_hat), ces = cor(s$tau_hat, ces)) })
cat(sprintf("  r(tau_hat, thres_hat) mean %.2f | r(tau_hat, CES-D) mean %.3f (over %d reps)\n",
            mean(rs["thres", ]), mean(rs["ces", ]), ncol(rs)))
o <- sapply(c("thres", "medRT"), function(v) NA)
x <- merge(f, do.call(rbind, lapply(split(dat, dat$subject), function(s) data.frame(subject = s$subject[1], ces = s$cesScore[1]))))
cat(sprintf("  observed: r(tau, thres) %.2f | r(tau, CES-D) %.3f\n", cor(x$tau, x$thres), cor(x$tau, x$ces)))
for (gp in c("0.05", "0.10")) {
  d <- read.csv(file.path(R, if (gp == "0.05") "nulltau_corrected_image_20260928_222147" else "nulltau_corrected_image_20260928_224927", "per_rep.csv"))
  v <- -d$recovered_dtau
  cat(sprintf("  starting-point gap %s: recovered %.3f s (sd %.3f, %d reps)\n", gp, mean(v), sd(v), length(v)))
}

cat("\n== individual recovery ==\n")
rc <- read.csv(file.path(R, "recovery_redesign_tau_20260928_192519", "recovery_summary.csv"))
for (i in seq_len(nrow(rc))) cat(sprintf("  %d trials: r = %.2f [%.2f, %.2f] | error SD %.2f vs true SD %.2f | bias %+.2f s\n",
  rc$n_trials[i], rc$r_tau[i], rc$r_lo[i], rc$r_hi[i], rc$sd_err_tau[i], rc$sd_true_tau[i], rc$bias[i]))
cat("\n== non-decision-time bound ladder ==\n")
print(read.csv(file.path(R, "nDT_bound_ladder_engine_20260928_192521", "ladder_summary.csv")), row.names = FALSE)
sink()
