#!/usr/bin/env Rscript
# ==============================================================================
# RECOVERY REDESIGN — onset tau (4j mean/contrast) @ nDT bound 3.0
# ==============================================================================
# Answers lab feedback #4: "scatter has spread but correlation is low — why?
#   are the all-positive estimates a bias? why are there negative recoveries
#   when true tau is all positive? is it a bounds difference?"
#
# The 2026-06-28 run (results/tDDM/parameter_recovery_4j_nDT3_20260628_012335)
# reported r(tau) = 0.083 but was NOT a fair test:
#   (a) truth was drawn U(1.30, 2.40) => SD_true = 0.30, while the estimator's
#       error SD is ~1.5 s. Even a PERFECT estimator caps at r ~ 0.20 there.
#       => range restriction, not (only) estimator failure.
#   (b) fit bound on tau was (-2, 3) while truth was strictly positive, so noise
#       spills into negative tau -- and tau <= 0 is a DEGENERATE BRANCH of the
#       likelihood (variant4j_mean_contrast.R:63-64): it collapses to a
#       single-phase DDM driven by drift2 only, so d1 drops out of the
#       likelihood entirely. 26-30% of recoveries landed there.
#
# This script runs a 2x2 factorial that separates those two causes, and reports
# the ATTAINABLE ceiling r for every cell so the observed r is interpretable.
#
#   ARM              tau fit bound      tau truth range     purpose
#   ---------------  -----------------  ------------------  ---------------------
#   legacy           (-2, 3)            U(1.30, 2.40)       reproduce r ~ .08
#   boundfix         ( 0, 3)            U(1.30, 2.40)       isolate degenerate branch
#   widetruth        (-2, 3)            U(0.30, 2.70)       isolate range restriction
#   redesign         ( 0, 3)            U(0.30, 2.70)       fair test of identifiability
#
# Each (arm x trial count x subject) is simulated+refit REPS times with different
# noise seeds, giving a direct Monte-Carlo SE(tau_hat) per subject -- the honest
# quantity to report instead of, or alongside, a correlation.
#
# Usage:
#   Rscript scripts/tDDM/recovery_redesign_tau.R                     # full
#   Rscript scripts/tDDM/recovery_redesign_tau.R --fast              # smoke test
#   Rscript scripts/tDDM/recovery_redesign_tau.R --N 100 --reps 4 --cores 10
# ==============================================================================

suppressPackageStartupMessages({ library(DEoptim); library(parallel) })

args <- commandArgs(trailingOnly = TRUE)
FAST_MODE <- any(args %in% c("--fast", "--dry"))
get_arg <- function(flag, default) {
  if (flag %in% args) { i <- which(args == flag); if (i < length(args)) return(args[i + 1]) }
  default
}

script_dir <- tryCatch(dirname(sys.frame(1)$ofile), error = function(e) {
  ca <- commandArgs(trailingOnly = FALSE); fa <- sub("^--file=", "", ca[grep("^--file=", ca)])
  if (length(fa) > 0) dirname(fa) else getwd()
})
source(file.path(script_dir, "..", "..", "config.R"))
RESULT_ROOT <- file.path(RESULTS_DIR, "tDDM"); dir.create(RESULT_ROOT, showWarnings = FALSE, recursive = TRUE)

if (FAST_MODE) {
  TRIAL_COUNTS <- c(20, 80); N_SUBJECTS <- 12; REPS <- 2; ITERMAX <- 80; NP <- 40
} else {
  TRIAL_COUNTS <- c(20, 40, 80, 160, 320); N_SUBJECTS <- 100; REPS <- 4; ITERMAX <- 200; NP <- 80
}
TRIAL_COUNTS <- as.integer(strsplit(get_arg("--trials", paste(TRIAL_COUNTS, collapse = ",")), ",")[[1]])
N_SUBJECTS <- as.integer(get_arg("--N", N_SUBJECTS))
REPS       <- as.integer(get_arg("--reps", REPS))
ITERMAX    <- as.integer(get_arg("--itermax", ITERMAX))
NP         <- as.integer(get_arg("--np", NP))
N_CORES    <- as.integer(get_arg("--cores", max(1, detectCores() - 2)))
# 2026-09-26: keep the BEST-NLL run of R repeats, matching the canonical procedure
# (docs/PAPER/Memo_Analyses_0919.md S4).
RESTARTS   <- as.integer(get_arg("--restarts", 1))
NDT_UB     <- 3.0


# ---- arms -------------------------------------------------------------------
ARMS <- list(
  legacy    = list(tau_lb = -2, tau_truth = c(1.30, 2.40), label = "legacy (bound -2..3, truth 1.3-2.4)"),
  boundfix  = list(tau_lb =  0, tau_truth = c(1.30, 2.40), label = "bound-fix (bound 0..3, truth 1.3-2.4)"),
  widetruth = list(tau_lb = -2, tau_truth = c(0.30, 2.70), label = "wide-truth (bound -2..3, truth 0.3-2.7)"),
  redesign  = list(tau_lb =  0, tau_truth = c(0.30, 2.70), label = "REDESIGN (bound 0..3, truth 0.3-2.7)")
)

cat("=== RECOVERY REDESIGN: onset tau, 4j mean/contrast, nDT bound 3.0 ===\n")
cat(sprintf("  N=%d subjects x REPS=%d x trials={%s} x %d arms = %d fits\n",
            N_SUBJECTS, REPS, paste(TRIAL_COUNTS, collapse = ","), length(ARMS),
            N_SUBJECTS * REPS * length(TRIAL_COUNTS) * length(ARMS)))
cat(sprintf("  itermax=%d NP=%d restarts=%d cores=%d\n\n", ITERMAX, NP, RESTARTS, N_CORES))

# ----------------------------------------------- engine (SHARED) --------------
# 2026-08-10: this script used to carry its own verbatim copy of the Navarro-Fuss
# density, the split-at-tau transition density and the piecewise likelihood.  Four
# separate defects were traced to such copies drifting apart, so the definitions
# now live in ONE place.  Regression-checked with
#   Rscript scripts/tDDM/test_engine_equivalence.R --trans image|pipeline
source(file.path(MODEL_DIR, "tDDM", "tddm_engine.R"))
TRANS_MODE <- get_arg("--trans", "image")  # default was "pipeline" until 2026-09-27
cat(sprintf("transition density: %s\n", TRANS_MODE))

SIGMA   <- TDDM_SIGMA       # kept: simulate_piecewise_ddm() below works in the
EPS     <- TDDM_EPS         # UNRESCALED space and needs sigma explicitly
PENALTY <- TDDM_PENALTY
N_GRID  <- TDDM_CANON$ngrid

# td columns here are signal1 / signal2 / RT_ddm / choice_ddm, which is exactly
# what mk_ll_tau() expects once the two signal columns are named.
ll_4j_fit <- mk_ll_tau("signal1", "signal2", trans = get_transition(TRANS_MODE))

# 2026-09-27: body replaced by the shared bridge-corrected simulator (models/tDDM/tddm_engine.R,
# tddm_simulate_trial). The old Euler copy lengthened mean DT by 17-39 ms vs the likelihood.
simulate_piecewise_ddm <- function(d1, d2, thres, nDT, tau, signal1, signal2, sigma = SIGMA, dt = 0.001, max_time = 10) {
  r <- tddm_simulate_trial(d1, d2, thres, nDT, tau, signal1, signal2, sg = sigma, dt = dt, mx = max_time)
  list(rt = r$rt, choice = r$ch)
}

# ----------------------------------------------- data + signal pool ------------
cat("Step 1: signal pool from real data...\n")
files <- list.files(DATA_DIR, pattern = "^study4_.*\\.csv$", full.names = TRUE)
real_data <- do.call(rbind, lapply(files, function(p) {
  df <- read.csv(p, stringsAsFactors = FALSE)
  df$group <- toupper(sub("\\.csv$", "", sub("^study4_", "", basename(p)))); df }))
real_data <- real_data[is.finite(real_data$choice_RT) & real_data$choice_RT > 0.15, ]
real_data$mean_sig     <- (real_data$RPE_scale + real_data$EPE_X1_scale) / 2
real_data$contrast_sig <- (real_data$RPE_scale - real_data$EPE_X1_scale) / 2
signal_pool <- real_data[, c("mean_sig", "contrast_sig")]
names(signal_pool) <- c("s1", "s2")
signal_pool <- signal_pool[complete.cases(signal_pool), ]
cat(sprintf("  pool=%d pairs; r(mean,contrast)=%.3f; Var(mean)=%.3f Var(contrast)=%.3f (ratio %.1fx)\n\n",
            nrow(signal_pool), cor(signal_pool$s1, signal_pool$s2),
            var(signal_pool$s1), var(signal_pool$s2), var(signal_pool$s1) / var(signal_pool$s2)))

# ----------------------------------------------- truth -------------------------
# d1/d2/thres/nDT are IDENTICAL across arms (same seed) so arms differ only in
# what the feedback asks about: the tau truth range and the tau fit bound.
set.seed(20260804)
base_truth <- data.frame(
  id    = 1:N_SUBJECTS,
  d1    = runif(N_SUBJECTS, -0.50, 1.75),
  d2    = runif(N_SUBJECTS, -2.00, 2.00),
  thres = runif(N_SUBJECTS,  1.10, 2.70),
  nDT   = runif(N_SUBJECTS,  1.20, 3.00))
u_tau <- runif(N_SUBJECTS)   # common quantiles -> same rank order of tau in every arm

# NOTE 2026-08-04: do NOT set VTR here.  DEoptim stops as soon as bestval <= VTR,
# and the NLLs here are ~30-300, so the VTR = 1e5 used by the older recovery
# scripts made DEoptim return the best of its INITIAL RANDOM POPULATION after
# 0 generations -- i.e. no optimisation at all.  Verified: on a smooth 5-D test
# problem, VTR=1e5 gives bestval 29.5 in 0 iterations vs the true optimum 25.0
# reached in 200 iterations with the default VTR = -Inf.
fit_ctrl <- DEoptim.control(itermax = ITERMAX, NP = NP, trace = FALSE)

recover_one <- function(tp, n_trials, tau_lb, seed) {
  set.seed(seed)
  idx <- sample(seq_len(nrow(signal_pool)), n_trials, replace = TRUE)
  ts <- signal_pool[idx, ]
  rt <- numeric(n_trials); ch <- integer(n_trials)
  for (j in seq_len(n_trials)) {
    sim <- simulate_piecewise_ddm(tp$d1, tp$d2, tp$thres, tp$nDT, tp$tau, ts$s1[j], ts$s2[j])
    rt[j] <- sim$rt; ch[j] <- sim$choice
  }
  td <- data.frame(signal1 = ts$s1, signal2 = ts$s2, RT_ddm = rt, choice_ddm = ch)
  td <- td[td$RT_ddm > 0.15 & td$RT_ddm <= 10, ]
  if (nrow(td) < 5) return(NULL)
  lo <- c(-3, -3, 0.5, 0.01, tau_lb); up <- c(3, 3, 4.0, NDT_UB, 3)
  fit <- NULL
  for (k in seq_len(RESTARTS)) {
    fk <- tryCatch(DEoptim(ll_4j_fit, lo, up, control = fit_ctrl, td = td), error = function(e) NULL)
    if (!is.null(fk) && (is.null(fit) || fk$optim$bestval < fit$optim$bestval)) fit <- fk
  }
  if (is.null(fit)) return(NULL)
  est <- fit$optim$bestmem
  de_iter <- fit$optim$iter        # guard against a silent early stop (see VTR note)
  # identifiability window: tau is only informed by trials whose decision time
  # exceeds tau (those are the only ones where phase 2 is ever entered)
  dt_true  <- td$RT_ddm - tp$nDT
  frac_two <- mean(dt_true > tp$tau)          # trials informing d2 / tau under TRUTH
  window   <- max(dt_true)                    # largest identifiable tau
  data.frame(
    subject = tp$id, n_trials = n_trials, rep = tp$rep, n_valid = nrow(td),
    true_d1 = tp$d1, est_d1 = est[1], true_d2 = tp$d2, est_d2 = est[2],
    true_thres = tp$thres, est_thres = est[3], true_nDT = tp$nDT, est_nDT = est[4],
    true_tau = tp$tau, est_tau = est[5],
    frac_two_phase = frac_two, ident_window = window,
    nll = fit$optim$bestval, de_iter = de_iter,
    at_tau_lb = as.integer(est[5] <= tau_lb + 0.01),
    at_tau_ub = as.integer(est[5] >= 3 - 0.01),
    degenerate = as.integer(est[5] <= 0),
    accept_rate = mean(ch), stringsAsFactors = FALSE)
}

# ----------------------------------------------- run ---------------------------
cat("Step 2: simulate + recover across arms...\n")
t0 <- Sys.time(); all_rows <- list()
ARM_SEL <- strsplit(get_arg("--arms", paste(names(ARMS), collapse=",")), ",")[[1]]
ARMS <- ARMS[ARM_SEL]
for (arm_name in names(ARMS)) {
  arm <- ARMS[[arm_name]]
  tau_truth <- arm$tau_truth[1] + u_tau * (arm$tau_truth[2] - arm$tau_truth[1])
  truth <- base_truth; truth$tau <- tau_truth
  for (n_trials in TRIAL_COUNTS) {
    jobs <- expand.grid(si = seq_len(N_SUBJECTS), rp = seq_len(REPS))
    res <- mclapply(seq_len(nrow(jobs)), function(k) {
      si <- jobs$si[k]; rp <- jobs$rp[k]
      tp <- truth[si, ]; tp$rep <- rp
      seed <- 1e6 * which(names(ARMS) == arm_name) + 1e4 * match(n_trials, TRIAL_COUNTS) + 100 * si + rp
      recover_one(tp, n_trials, arm$tau_lb, seed)
    }, mc.cores = N_CORES, mc.preschedule = FALSE)
    res <- Filter(function(z) is.data.frame(z), res)
    if (length(res)) {
      df <- do.call(rbind, res); df$arm <- arm_name
      all_rows[[length(all_rows) + 1]] <- df
    }
    cat(sprintf("  [%-9s] n_trials=%-4d: %d/%d fits (%.1f min)\n", arm_name, n_trials,
                length(res), nrow(jobs), as.numeric(difftime(Sys.time(), t0, units = "mins"))))
    flush.console()
  }
}
results_df <- do.call(rbind, all_rows)
cat(sprintf("\nTotal: %.1f min, %d fits\n", as.numeric(difftime(Sys.time(), t0, units = "mins")), nrow(results_df)))
cat(sprintf("DEoptim generations actually run: median %.0f, min %.0f (0 = the VTR early-stop bug)\n\n",
            median(results_df$de_iter), min(results_df$de_iter)))
if (min(results_df$de_iter) == 0)
  cat("*** WARNING: some fits ran 0 generations -- check DEoptim.control for VTR ***\n\n")

# ----------------------------------------------- metrics -----------------------
fisher_ci <- function(r, n) {
  if (!is.finite(r) || n < 4) return(c(NA, NA))
  z <- atanh(max(min(r, 0.9999), -0.9999)); se <- 1 / sqrt(n - 3)
  tanh(c(z - 1.96 * se, z + 1.96 * se))
}

summ <- list()
for (arm_name in names(ARMS)) for (nt in TRIAL_COUNTS) {
  s <- results_df[results_df$arm == arm_name & results_df$n_trials == nt, ]
  if (nrow(s) < 5) next
  # per-subject aggregation across reps (mean tau_hat, MC SE of tau_hat)
  agg <- do.call(rbind, lapply(split(s, s$subject), function(g) data.frame(
    subject = g$subject[1], true_tau = g$true_tau[1],
    mean_est = mean(g$est_tau), mc_se = if (nrow(g) > 1) sd(g$est_tau) else NA_real_)))
  sd_true <- sd(s$true_tau); sd_err <- sd(s$est_tau - s$true_tau)
  r_obs   <- cor(s$true_tau, s$est_tau)
  ci      <- fisher_ci(r_obs, nrow(s))
  # attainable r given THIS truth spread and THIS estimator noise
  r_ceiling <- sd_true / sqrt(sd_true^2 + sd_err^2)
  summ[[length(summ) + 1]] <- data.frame(
    arm = arm_name, n_trials = nt, n_fits = nrow(s), n_subjects = nrow(agg),
    r_tau = r_obs, r_lo = ci[1], r_hi = ci[2], r_ceiling = r_ceiling,
    r_subjmean = cor(agg$true_tau, agg$mean_est),
    slope = unname(coef(lm(s$est_tau ~ s$true_tau))[2]),
    bias = mean(s$est_tau - s$true_tau), rmse = sqrt(mean((s$est_tau - s$true_tau)^2)),
    sd_true_tau = sd_true, sd_err_tau = sd_err,
    mc_se_tau = median(agg$mc_se, na.rm = TRUE),
    pct_negative = 100 * mean(s$est_tau <= 0),
    pct_degenerate = 100 * mean(s$degenerate),
    pct_at_lb = 100 * mean(s$at_tau_lb), pct_at_ub = 100 * mean(s$at_tau_ub),
    pct_tau_above_window = 100 * mean(s$true_tau > s$ident_window),
    median_frac_two_phase = median(s$frac_two_phase),
    r_d1 = cor(s$true_d1, s$est_d1), r_d2 = cor(s$true_d2, s$est_d2),
    r_thres = cor(s$true_thres, s$est_thres), r_nDT = cor(s$true_nDT, s$est_nDT),
    stringsAsFactors = FALSE)
}
summary_df <- do.call(rbind, summ)

cat("=========================================================================\n")
cat("  TAU RECOVERY BY ARM  (r_ceiling = best attainable r given truth spread\n")
cat("                        and this estimator's error SD)\n")
cat("=========================================================================\n\n")
cat(sprintf("%-10s %5s %7s %7s %8s %8s %8s %7s %7s\n",
            "arm", "ntrl", "r_obs", "r_ceil", "SD_true", "SD_err", "MC_SE", "%neg", "%2phase"))
for (i in seq_len(nrow(summary_df))) with(summary_df[i, ],
  cat(sprintf("%-10s %5d %7.3f %7.3f %8.3f %8.3f %8.3f %6.0f%% %6.0f%%\n",
              arm, n_trials, r_tau, r_ceiling, sd_true_tau, sd_err_tau, mc_se_tau,
              pct_negative, 100 * median_frac_two_phase)))
cat("\n--- other parameters (r) ---\n")
cat(sprintf("%-10s %5s %7s %7s %7s %7s\n", "arm", "ntrl", "d_mean", "d_cont", "thres", "nDT"))
for (i in seq_len(nrow(summary_df))) with(summary_df[i, ],
  cat(sprintf("%-10s %5d %7.3f %7.3f %7.3f %7.3f\n", arm, n_trials, r_d1, r_d2, r_thres, r_nDT)))
cat("\n")

# ----------------------------------------------- save --------------------------
tstamp <- format(Sys.time(), "%Y%m%d_%H%M%S")
outdir <- file.path(RESULT_ROOT, sprintf("recovery_redesign_tau_%s", tstamp))
dir.create(outdir, recursive = TRUE)
write.csv(results_df, file.path(outdir, "recovery_fits.csv"), row.names = FALSE)
write.csv(summary_df, file.path(outdir, "recovery_summary.csv"), row.names = FALSE)

# per-subject MC SE table (the headline number to report instead of r)
mcse_tbl <- do.call(rbind, lapply(split(results_df, list(results_df$arm, results_df$n_trials), drop = TRUE),
  function(s) {
    agg <- do.call(rbind, lapply(split(s, s$subject), function(g)
      data.frame(mc_se = if (nrow(g) > 1) sd(g$est_tau) else NA_real_)))
    data.frame(arm = s$arm[1], n_trials = s$n_trials[1],
               mc_se_median = median(agg$mc_se, na.rm = TRUE),
               mc_se_mean = mean(agg$mc_se, na.rm = TRUE),
               sd_true_between = sd(s$true_tau))
  }))
write.csv(mcse_tbl, file.path(outdir, "tau_mc_se.csv"), row.names = FALSE)

# ----------------------------------------------- figures -----------------------
if (requireNamespace("ggplot2", quietly = TRUE)) {
  library(ggplot2)
  arm_lab <- setNames(sapply(ARMS, `[[`, "label"), names(ARMS))
  results_df$arm_f <- factor(results_df$arm, levels = names(ARMS), labels = arm_lab)
  summary_df$arm_f <- factor(summary_df$arm, levels = names(ARMS), labels = arm_lab)

  p1 <- ggplot(results_df[results_df$n_trials %in% c(20, 80, 320), ], aes(true_tau, est_tau)) +
    geom_hline(yintercept = 0, colour = "grey70", linewidth = 0.3) +
    geom_point(aes(colour = degenerate == 1), alpha = 0.4, size = 1.1) +
    scale_colour_manual(values = c(`FALSE` = "grey25", `TRUE` = "firebrick"),
                        labels = c("two-phase fit", "degenerate (tau<=0)"), name = NULL) +
    geom_abline(slope = 1, intercept = 0, linetype = "dashed", colour = "red") +
    geom_smooth(method = "lm", se = FALSE, colour = "blue", linewidth = 0.6) +
    facet_grid(arm_f ~ n_trials, labeller = labeller(n_trials = label_both)) +
    labs(title = "Onset tau recovery: bound and truth-range decomposed",
         subtitle = "red dashed = identity; blue = OLS. Red points = degenerate branch (tau<=0 collapses the model to drift2-only)",
         x = "True tau (s)", y = "Estimated tau (s)") +
    theme_minimal(base_size = 10) + theme(legend.position = "top")
  ggsave(file.path(outdir, "fig_tau_scatter_by_arm.png"), p1, width = 11, height = 9.5, dpi = 200)

  long_r <- rbind(
    data.frame(arm_f = summary_df$arm_f, n_trials = summary_df$n_trials,
               value = summary_df$r_tau, kind = "observed r"),
    data.frame(arm_f = summary_df$arm_f, n_trials = summary_df$n_trials,
               value = summary_df$r_ceiling, kind = "ceiling r (attainable)"))
  p2 <- ggplot(long_r, aes(factor(n_trials), value, fill = kind)) +
    geom_col(position = position_dodge(0.8), width = 0.7) +
    geom_hline(yintercept = 0, linewidth = 0.3) +
    facet_wrap(~ arm_f, nrow = 2) +
    scale_fill_manual(values = c("observed r" = "#3b6ea5", "ceiling r (attainable)" = "#c9c9c9"), name = NULL) +
    labs(title = "Observed vs attainable tau recovery correlation",
         subtitle = "grey = the best r ANY estimator could reach given this truth spread and this error SD",
         x = "Trials per subject", y = "r(true tau, estimated tau)") +
    theme_minimal(base_size = 11) + theme(legend.position = "top")
  ggsave(file.path(outdir, "fig_r_vs_ceiling.png"), p2, width = 10, height = 6.5, dpi = 200)

  p3 <- ggplot(summary_df, aes(n_trials, mc_se_tau, colour = arm_f)) +
    geom_line(linewidth = 0.8) + geom_point(size = 2) +
    geom_hline(aes(yintercept = sd_true_tau, linetype = "between-subject SD of true tau"),
               data = summary_df[!duplicated(summary_df$arm), ], colour = "grey40") +
    scale_x_log10(breaks = TRIAL_COUNTS) +
    scale_linetype_manual(values = c("between-subject SD of true tau" = "dashed"), name = NULL) +
    labs(title = "Per-subject Monte-Carlo SE of tau_hat vs between-subject signal",
         subtitle = "individual tau is recoverable only where SE falls below the between-subject SD",
         x = "Trials per subject (log)", y = "SE(tau_hat), seconds", colour = NULL) +
    theme_minimal(base_size = 11) + theme(legend.position = "right")
  ggsave(file.path(outdir, "fig_tau_se_vs_trials.png"), p3, width = 9.5, height = 5, dpi = 200)

  p4 <- ggplot(summary_df, aes(factor(n_trials), pct_negative, fill = arm_f)) +
    geom_col(position = position_dodge(0.8), width = 0.7) +
    labs(title = "Share of fits landing in the degenerate branch (tau_hat <= 0)",
         subtitle = "tau<=0 is NOT 'emotion-first': the likelihood collapses to drift2-only and d_mean drops out",
         x = "Trials per subject", y = "% of fits with tau_hat <= 0", fill = NULL) +
    theme_minimal(base_size = 11) + theme(legend.position = "top")
  ggsave(file.path(outdir, "fig_degenerate_rate.png"), p4, width = 9.5, height = 5.5, dpi = 200)
  cat("  figures saved\n")
}

cat(sprintf("\nOutput: %s\n", outdir))
