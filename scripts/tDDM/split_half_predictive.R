#!/usr/bin/env Rscript
# ==============================================================================
# SPLIT-HALF PREDICTIVE VALIDATION  --  Variants 4a, 4i, 4j
# ==============================================================================
# Question (user): randomly split the data in half, fit the model on one half,
#   and see whether it RELIABLY PREDICTS the behavioral results in the 2nd half.
#
# WHAT "PREDICT THE BEHAVIORAL RESULTS" MEANS HERE (two complementary tests):
#
#   [PRIMARY] EFFECT REPLICATION.  Fit BOTH halves independently (subject-level,
#     stratified). For each model, take its signature GROUP effect and ask: does
#     it appear in BOTH halves with the SAME sign?  Reliability = fraction of K
#     random splits where the two halves agree.  We also report the effect SIZE
#     (Cohen's d) in each half so a null cannot masquerade as ~50% sign-agreement.
#       - 4a            : d_mean   (expect HC > DEP)
#       - 4i_epe_first  : tau      (expect HC < DEP)   [memory: 1.646 vs 2.001]
#       - 4i_rpe_first  : tau      (expect HC < DEP)
#       - 4j_mean_first : tau      (no a-priori; replication = halves agree)
#       - 4j_contrast_first : tau  (no a-priori)
#
#   [SECONDARY] HELD-OUT PREDICTION.  Apply the FIT-half group-mean params to the
#     TEST-half trials and score them:
#       - held-out NLL / trial  (all models; = cross-validated model comparison)
#       - 4a ONLY: analytical accept-rate & median-RT prediction vs observed
#         (the most literal "predict behavioral results"; this is what exposed
#          4a's structural accept-rate miss, which NLL hides because RT density
#          dominates the likelihood).
#
# WHY THIS DESIGN (see git/advisor notes):
#   * 20 trials/subject -> trial-level individual splits (10/10) are unsound
#     (recovery r=0.47 @ N=20). We split SUBJECTS; each fit-half subject keeps
#     all 20 trials and is fit with the primary individual pipeline, then
#     averaged to group means (pooling trials attenuates drift).
#   * 4i orthogonalization (PCA / Gram-Schmidt) is FIT TO DATA -> computed on the
#     fit-half only for held-out scoring (no leakage). For the independent
#     replication fits, each half uses its own orthogonalization (self-contained).
#   * Piecewise likelihood (4i/4j) is ~30x the 4a per-trial cost -> mclapply over
#     subjects + a reduced DEoptim budget (tau is stable across budgets; group
#     means average out individual-fit noise via CLT).
#
# Usage:
#   Rscript scripts/tDDM/split_half_predictive.R --time-only      # time ONE split
#   Rscript scripts/tDDM/split_half_predictive.R --K 20           # 20 splits
#   Rscript scripts/tDDM/split_half_predictive.R --K 30 --models 4a,4i_epe_first
#   Rscript scripts/tDDM/split_half_predictive.R --fast           # K=3 smoke test
# ==============================================================================

suppressPackageStartupMessages({
  library(DEoptim)
  library(parallel)
})

# ------------------------------------------------------------------ args ------
args <- commandArgs(trailingOnly = TRUE)
FAST_MODE <- any(args %in% c("--fast", "--dry"))
TIME_ONLY <- any(args %in% c("--time-only"))
get_arg <- function(flag, default) {
  i <- which(args == flag); if (length(i) == 1 && i < length(args)) args[i + 1] else default
}
K_SPLITS <- as.integer(get_arg("--K", if (FAST_MODE) 3 else 30))
SEED     <- as.integer(get_arg("--seed", 20260605))
N_CORES  <- max(1, as.integer(get_arg("--cores", detectCores() - 2)))

# Reduced DEoptim budget for the per-subject split fits. VALIDATED (same-seed
# split): itermax=30/NP=30 vs 80/45 preserves the tau replication SIGN for all 5
# models (all HC<DEP, agree=1) at 3.4x speed (6.9 vs 23.3 min/split). Magnitudes
# wobble within single-split noise; group means over K splits average it out.
# Bump --itermax/--np for a final effect-SIZE run if precise Delta-tau is needed.
IT_PW <- as.integer(get_arg("--itermax", 30))   # piecewise (4i/4j) itermax
NP_PW <- as.integer(get_arg("--np", 30))        # piecewise NP
# 2026-09-23: keep the BEST-NLL run of R repeats of each piecewise fit, matching
# the canonical procedure (docs/PAPER/Memo_Analyses_0919.md S4).  The static
# models converge in one run and are left alone.
RESTARTS <- as.integer(get_arg("--restarts", 1))
# 2026-08-10: default raised 1.5 -> 3.0. The canonical suite moved to nDT <= 3.0
# on 2026-06-26 because the 1.5 bound censored 74% of subjects at the ceiling and
# the resulting "nDT is NS" was an artefact of that. Leaving this at 1.5 silently
# ran the 2026-08-10 split-half under the OLD bound while every other analysis
# used 3.0.
NDT_UB <- as.numeric(get_arg("--ndt-ub", 3.0))  # nDT upper bound (canonical = 3.0)
# --trans pipeline (DEFAULT, historical) | image (corrected absorbing transition)
TRANS_MODE <- get_arg("--trans", "image")  # default was "pipeline" until 2026-09-27
cat(sprintf("transition density: %s\n", TRANS_MODE))

MODELS_ARG <- get_arg("--models",
                      "4a,4i_epe_first,4i_rpe_first,4j_mean_first,4j_contrast_first")
SEL_MODELS <- strsplit(MODELS_ARG, ",")[[1]]

# ------------------------------------------------------------- paths ----------
script_dir <- tryCatch(dirname(sys.frame(1)$ofile), error = function(e) {
  ca <- commandArgs(trailingOnly = FALSE); fa <- sub("^--file=", "", ca[grep("^--file=", ca)])
  if (length(fa) > 0) dirname(fa) else getwd()
})
source(file.path(script_dir, "..", "..", "config.R"))
RESULT_ROOT <- file.path(RESULTS_DIR, "tDDM"); dir.create(RESULT_ROOT, FALSE, TRUE)

cat("=== Split-Half Predictive Validation (4a / 4i / 4j) ===\n")
cat(sprintf("  K=%d splits | seed=%d | cores=%d | piecewise budget itermax=%d NP=%d restarts=%d\n",
            K_SPLITS, SEED, N_CORES, IT_PW, NP_PW, RESTARTS))
cat(sprintf("  nDT UB=%.2f | transition density = %s\n", NDT_UB, TRANS_MODE))
cat(sprintf("  Models: %s\n\n", paste(SEL_MODELS, collapse = ", ")))

# ==============================================================================
# Wiener density + piecewise engine (shared with variant4i/4j)
# ==============================================================================
# ---- engine (SHARED) ---------------------------------------------------------
# 2026-08-10: replaced this script's private copies of the Navarro-Fuss density,
# the split-at-tau transition density and the piecewise likelihood with the one
# shared definition.  Regression-checked with
#   Rscript scripts/tDDM/test_engine_equivalence.R --trans image|pipeline
source(file.path(MODEL_DIR, "tDDM", "tddm_engine.R"))
SIGMA <- TDDM_SIGMA; EPS <- TDDM_EPS; PENALTY <- TDDM_PENALTY
N_GRID <- TDDM_CANON$ngrid
# Resolve the density ONCE. ll_piecewise_trial() takes `trans` as its 9th
# argument and DEFAULTS to the image density, so every call below passes it
# explicitly -- otherwise a --trans pipeline run would silently use the
# corrected density wherever the argument was omitted.
TRANS_FUN <- get_transition(TRANS_MODE)

# RESTORED 2026-08-10: these two live in this script (not the engine) and were
# deleted by accident when the private engine block was swapped out -- the block
# boundary ran past them. Recovered verbatim from commit c5e2cab. They are used
# by the 4a behavioural prediction (analytical accept-rate and median RT).
choice_prob_upper <- function(v, a, w = 0.5) {
  if (abs(v) < 1e-10) return(w)
  va <- v * a
  if (abs(va) > 500) return(if (va > 0) 1.0 else 0.0)
  max(0.001, min(0.999, (exp(-2 * va * (1 - w)) - 1) / (exp(-2 * va) - 1)))
}
predict_rt_median <- function(v, a, w, t_max = 12, n_grid = 400) {
  tg <- seq(0.01, t_max, length.out = n_grid); dt <- tg[2] - tg[1]
  dens <- vapply(tg, function(t) dwiener_nf(t, v, a, w), numeric(1))
  cdf <- cumsum(dens) * dt; tot <- cdf[n_grid]
  if (tot < 1e-3) return(NA_real_)
  tg[which(cdf / tot >= 0.5)[1]]
}

ll_4a <- function(p, td) {                                   # d_mean, thres, nDT
  bs <- 2 * p[2] / SIGMA; w <- 0.5; if (bs <= 0.01 || p[3] <= 0) return(PENALTY)
  nll <- 0
  for (i in 1:nrow(td)) {
    drift <- p[1] * td$mean_sig[i] / SIGMA; dt <- td$RT_ddm[i] - p[3]
    if (dt <= 0.001) { nll <- nll - log(EPS); next }
    dens <- if (td$choice_ddm[i] == 1) dwiener_nf(dt, -drift, bs, 1 - w) else dwiener_nf(dt, drift, bs, w)
    nll <- nll - log(max(dens, EPS))
  }
  if (!is.finite(nll)) PENALTY else nll
}
# Generic piecewise: phase-1 signal s1, phase-2 adds signal s2 (column names).
# This is exactly the shared engine's mk_ll_tau().
make_ll_pw <- function(s1col, s2col) mk_ll_tau(s1col, s2col, trans = TRANS_FUN)

# STATIC two-signal (red-team #2): both signals from t=0, tau FIXED at 0. This is
# the missing nested control between 4a {mean} and 4j {mean,contrast,tau}: it adds
# the second drift component WITHOUT a free onset, so comparing its held-out NLL to
# 4j isolates whether the predictive gain is from TIMING (tau) or merely from the
# extra drift signal. 4 params: d1, d2, thres, nDT.
make_ll_static2 <- function(s1col, s2col) function(p, td) {  # d1, d2, thres, nDT (tau=0)
  bs <- 2 * p[3] / SIGMA; w <- 0.5; if (bs <= 0.01 || p[4] <= 0) return(PENALTY)
  nll <- 0
  for (i in 1:nrow(td)) {
    drift2 <- (p[1] * td[[s1col]][i] + p[2] * td[[s2col]][i]) / SIGMA
    nll <- nll - log(max(ll_piecewise_trial(drift2, drift2, 0, bs, w, p[4],
                                            td$RT_ddm[i], td$choice_ddm[i], TRANS_FUN), EPS))
  }
  if (!is.finite(nll)) PENALTY else nll
}

PW_LO <- c(-3, -3, 0.5, 0.01, -2); PW_UP <- c(3, 3, 4, NDT_UB, 3)
MODEL_REGISTRY <- list(
  "4a" = list(name = "4a (mean, tau=0)", ll = ll_4a,
              lower = c(-3, 0.5, 0.01), upper = c(3, 4, NDT_UB),
              pnames = c("d_mean", "thres", "nDT"), k = 3,
              sig_param = "d_mean", expect = "HC>DEP", piecewise = FALSE),
  "4i_epe_first" = list(name = "4i EPE-first (EPE->RPE_orth)",
              ll = make_ll_pw("EPE_X1_scale", "RPE_orth"),
              lower = PW_LO, upper = PW_UP,
              pnames = c("d_EPE", "d_RPE_orth", "thres", "nDT", "tau"), k = 5,
              sig_param = "tau", expect = "HC<DEP", piecewise = TRUE),
  "4i_rpe_first" = list(name = "4i RPE-first (RPE->EPE_orth)",
              ll = make_ll_pw("RPE_scale", "EPE_orth"),
              lower = PW_LO, upper = PW_UP,
              pnames = c("d_RPE", "d_EPE_orth", "thres", "nDT", "tau"), k = 5,
              sig_param = "tau", expect = "HC<DEP", piecewise = TRUE),
  "4j_mean_first" = list(name = "4j mean-first (mean->contrast)",
              ll = make_ll_pw("mean_sig", "contrast_sig"),
              lower = PW_LO, upper = PW_UP,
              pnames = c("d_mean", "d_contrast", "thres", "nDT", "tau"), k = 5,
              sig_param = "tau", expect = "either", piecewise = TRUE),
  "4j_contrast_first" = list(name = "4j contrast-first (contrast->mean)",
              ll = make_ll_pw("contrast_sig", "mean_sig"),
              lower = PW_LO, upper = PW_UP,
              pnames = c("d_contrast", "d_mean", "thres", "nDT", "tau"), k = 5,
              sig_param = "tau", expect = "either", piecewise = TRUE),
  # red-team #2: nested static two-signal control (tau=0, both drifts). Fit at the
  # SAME piecewise budget as 4j so the ONLY difference is the free onset tau ->
  # held-out NLL(4c) vs NLL(4j) isolates the predictive value of TIMING.
  "4c_mean_static" = list(name = "4c static (mean+contrast, tau=0)",
              ll = make_ll_static2("mean_sig", "contrast_sig"),
              lower = c(-3, -3, 0.5, 0.01), upper = c(3, 3, 4, NDT_UB),
              pnames = c("d_mean", "d_contrast", "thres", "nDT"), k = 4,
              sig_param = "d_mean", expect = "HC>DEP", piecewise = TRUE)
)
SEL_MODELS <- SEL_MODELS[SEL_MODELS %in% names(MODEL_REGISTRY)]
stopifnot(length(SEL_MODELS) > 0)

# ==============================================================================
# Data + fixed regressors
# ==============================================================================
cat("Step 1: Loading data...\n")
files <- list.files(DATA_DIR, pattern = "^study4_.*\\.csv$", full.names = TRUE)
data <- do.call(rbind, lapply(files, function(p) {
  df <- read.csv(p, stringsAsFactors = FALSE)
  df$group <- toupper(sub("\\.csv$", "", sub("^study4_", "", basename(p)))); df
}))
if ("sub" %in% names(data) && !"subject" %in% names(data)) names(data)[names(data) == "sub"] <- "subject"
data <- data[is.finite(data$choice_RT) & data$choice_RT > 0.15 & data$choice_RT <= 10, ]
data$choice_ddm   <- 1 - data$choice
data$RT_ddm       <- data$choice_RT
data$mean_sig     <- (data$RPE_scale + data$EPE_X1_scale) / 2   # 4a/4j common signal
data$contrast_sig <- (data$RPE_scale - data$EPE_X1_scale) / 2   # 4j differentiation
data$uid <- paste(data$group, data$subject, sep = "_")

subs_HC  <- unique(data$uid[data$group == "HC"])
subs_DEP <- unique(data$uid[data$group == "DEP"])
cat(sprintf("  HC: %d subjects, DEP: %d subjects, %d trials\n\n",
            length(subs_HC), length(subs_DEP), nrow(data)))

# Orthogonalization fit to a data set, applied to (possibly another) data set.
compute_ortho <- function(df) {
  list(pca = prcomp(df[, c("RPE_scale", "EPE_X1_scale")], center = TRUE, scale. = TRUE),
       b_re = coef(lm(RPE_scale ~ EPE_X1_scale, data = df))[2],
       b_er = coef(lm(EPE_X1_scale ~ RPE_scale, data = df))[2])
}
apply_ortho <- function(df, o) {
  pcs <- predict(o$pca, df[, c("RPE_scale", "EPE_X1_scale")])
  df$PC1 <- pcs[, 1]; df$PC2 <- pcs[, 2]
  df$RPE_orth <- df$RPE_scale - o$b_re * df$EPE_X1_scale
  df$EPE_orth <- df$EPE_X1_scale - o$b_er * df$RPE_scale
  df
}
NEED_ORTHO <- any(vapply(SEL_MODELS, function(m) MODEL_REGISTRY[[m]]$piecewise &&
                           grepl("4i", m), logical(1)))

# ==============================================================================
# Fit all subjects of a data frame for one model -> per-subject param matrix
# ==============================================================================
fit_individuals <- function(dat, model) {
  uids <- unique(dat$uid)
  rows <- mclapply(uids, function(u) {
    sd <- dat[dat$uid == u, ]; if (nrow(sd) < 5) return(NULL)
    it <- if (model$piecewise) IT_PW else 120
    np <- if (model$piecewise) NP_PW else 50
    nrun <- if (model$piecewise) RESTARTS else 1
    fit <- NULL
    for (k in seq_len(nrun)) {
      fk <- tryCatch(suppressWarnings(DEoptim(model$ll, model$lower, model$upper,
        control = DEoptim.control(itermax = it, NP = np, trace = FALSE), td = sd)),
        error = function(e) NULL)
      if (!is.null(fk) && (is.null(fit) || fk$optim$bestval < fit$optim$bestval)) fit <- fk
    }
    if (is.null(fit)) return(NULL)
    # 2026-09-28: every tau <= 0 gives the contrast model (identical likelihood), so tau is
    # clamped at 0 before any summary, matching the manuscript's canonical fits.
    est <- fit$optim$bestmem
    if ("tau" %in% model$pnames) est[model$pnames == "tau"] <- max(est[model$pnames == "tau"], 0)
    setNames(c(est, fit$optim$bestval, nrow(sd),
               if (sd$group[1] == "HC") 1 else 0),
             c(model$pnames, "negLogLik", "n_trials", "isHC"))
  }, mc.cores = N_CORES, mc.preschedule = FALSE)  # dynamic load-balance (no straggler)
  rows <- Filter(Negate(is.null), rows)
  if (length(rows) == 0) return(NULL)
  as.data.frame(do.call(rbind, rows))
}
group_means <- function(fitmat, model) {                     # named param vector per group
  list(
    HC  = vapply(model$pnames, function(p) mean(fitmat[[p]][fitmat$isHC == 1]), numeric(1)),
    DEP = vapply(model$pnames, function(p) mean(fitmat[[p]][fitmat$isHC == 0]), numeric(1))
  )
}
cohens_d <- function(a, b) {
  if (length(a) < 2 || length(b) < 2) return(NA_real_)
  (mean(a) - mean(b)) / sqrt((var(a) + var(b)) / 2)
}

# ==============================================================================
# One split: fit BOTH halves, compute replication + held-out prediction
# ==============================================================================
run_one_split <- function(k) {
  fit_uids  <- c(sample(subs_HC,  floor(length(subs_HC)  / 2)),
                 sample(subs_DEP, floor(length(subs_DEP) / 2)))
  fit_raw  <- data[data$uid %in%  fit_uids, ]
  test_raw <- data[!(data$uid %in% fit_uids), ]

  # Orthogonalization: each half self-contained for independent replication fits;
  # fit-half transform applied to test for held-out scoring (no leakage).
  if (NEED_ORTHO) {
    o_fit  <- compute_ortho(fit_raw); o_test <- compute_ortho(test_raw)
    fit_d   <- apply_ortho(fit_raw,  o_fit)
    test_d  <- apply_ortho(test_raw, o_test)   # independent fit of test half
    test_hd <- apply_ortho(test_raw, o_fit)    # held-out scoring with fit transform
  } else {
    fit_d <- fit_raw; test_d <- test_raw; test_hd <- test_raw
  }

  rows <- list()
  for (mid in SEL_MODELS) {
    m <- MODEL_REGISTRY[[mid]]
    f_fit  <- fit_individuals(fit_d,  m)
    f_test <- fit_individuals(test_d, m)
    if (is.null(f_fit) || is.null(f_test)) next

    sp <- m$sig_param
    sig_fit_HC  <- f_fit[[sp]][f_fit$isHC == 1];  sig_fit_DEP  <- f_fit[[sp]][f_fit$isHC == 0]
    sig_test_HC <- f_test[[sp]][f_test$isHC == 1]; sig_test_DEP <- f_test[[sp]][f_test$isHC == 0]
    delta_fit  <- mean(sig_fit_HC)  - mean(sig_fit_DEP)
    delta_test <- mean(sig_test_HC) - mean(sig_test_DEP)
    agree <- as.integer(sign(delta_fit) == sign(delta_test) && delta_fit != 0)

    gm <- group_means(f_fit, m)   # fit-half group-mean params for held-out scoring

    # held-out NLL/trial on test, using each group's fit-half mean params
    th_HC  <- test_hd[test_hd$group == "HC", ];  th_DEP <- test_hd[test_hd$group == "DEP", ]
    nll_HC_tot <- m$ll(gm$HC, th_HC); nll_DEP_tot <- m$ll(gm$DEP, th_DEP)  # compute once
    oos_nll_HC  <- nll_HC_tot  / nrow(th_HC)
    oos_nll_DEP <- nll_DEP_tot / nrow(th_DEP)

    row <- data.frame(
      split = k, model = mid, sig_param = sp,
      delta_fit = delta_fit, delta_test = delta_test,
      d_fit  = cohens_d(sig_fit_HC,  sig_fit_DEP),
      d_test = cohens_d(sig_test_HC, sig_test_DEP),
      agree = agree,
      sig_fitHC = mean(sig_fit_HC), sig_fitDEP = mean(sig_fit_DEP),
      sig_testHC = mean(sig_test_HC), sig_testDEP = mean(sig_test_DEP),
      oos_nll_HC = oos_nll_HC, oos_nll_DEP = oos_nll_DEP,
      oos_nll = (nll_HC_tot + nll_DEP_tot) / (nrow(th_HC) + nrow(th_DEP)),
      stringsAsFactors = FALSE
    )

    # 4a-only: analytical accept-rate & RT-median prediction vs observed (test)
    if (mid == "4a") {
      for (g in c("HC", "DEP")) {
        par <- gm[[g]]; bs <- 2 * par["thres"] / SIGMA
        tg <- test_hd[test_hd$group == g, ]
        drift <- par["d_mean"] * tg$mean_sig / SIGMA
        p_up <- vapply(drift, function(v) choice_prob_upper(v, bs, 0.5), numeric(1))
        mdr <- par["d_mean"] * mean(tg$mean_sig) / SIGMA
        prt <- (if (choice_prob_upper(mdr, bs, 0.5) >= 0.5)
                  predict_rt_median(-mdr, bs, 0.5) else predict_rt_median(mdr, bs, 0.5)) + par["nDT"]
        row[[paste0("accept_obs_", g)]]  <- mean(tg$choice_ddm == 1)
        row[[paste0("accept_pred_", g)]] <- mean(p_up)
        row[[paste0("rt_obs_", g)]]      <- median(tg$RT_ddm)
        row[[paste0("rt_pred_", g)]]     <- prt
      }
    }
    rows[[mid]] <- row
  }
  if (length(rows) == 0) return(NULL)
  # bind allowing differing columns (4a has extra)
  all_cols <- unique(unlist(lapply(rows, names)))
  do.call(rbind, lapply(rows, function(r) { r[setdiff(all_cols, names(r))] <- NA; r[all_cols] }))
}

# ==============================================================================
# --time-only: run ONE split, report wall-clock, project to K
# ==============================================================================
set.seed(SEED)
if (TIME_ONLY) {
  cat("Timing ONE full split (both halves, all selected models, parallel)...\n")
  t0 <- Sys.time(); r1 <- run_one_split(1)
  el <- as.numeric(difftime(Sys.time(), t0, units = "secs"))
  cat(sprintf("  -> %.1f s/split\n\n", el))
  cat("PROJECTED TOTAL RUNTIME:\n")
  for (kk in c(10, 20, 30, 50, 100)) cat(sprintf("  K=%-4d : ~%.1f min\n", kk, kk * el / 60))
  cat("\nPer-model signature (this one split, fit half vs test half):\n")
  for (i in seq_len(nrow(r1))) cat(sprintf("  %-20s %s: fitΔ=%.3f testΔ=%.3f agree=%d\n",
      r1$model[i], r1$sig_param[i], r1$delta_fit[i], r1$delta_test[i], r1$agree[i]))
  quit(save = "no")
}

# ==============================================================================
# Run K splits
# ==============================================================================
cat(sprintf("Step 2: Running %d random split-half iterations...\n", K_SPLITS))
# Create output dir up front so we can save INCREMENTALLY (long run -> monitor /
# survive a crash). per_split CSV is rewritten after every split.
ts <- format(Sys.time(), "%Y%m%d_%H%M%S")
outdir <- file.path(RESULT_ROOT, sprintf("split_half_predictive_%s_%s", TRANS_MODE, ts)); dir.create(outdir, recursive = TRUE)
per_split_path <- file.path(outdir, "split_half_per_split.csv")
cat(sprintf("  Output dir (incremental): %s\n", outdir)); flush.console()

t_start <- Sys.time()
all_rows <- list()
for (k in seq_len(K_SPLITS)) {
  all_rows[[k]] <- tryCatch(run_one_split(k),
                            error = function(e) { cat(sprintf("  split %d ERROR: %s\n", k, e$message)); NULL })
  # incremental save of everything completed so far
  done <- all_rows[!vapply(all_rows, is.null, logical(1))]
  if (length(done) > 0) write.csv(do.call(rbind, done), per_split_path, row.names = FALSE)
  el <- as.numeric(difftime(Sys.time(), t_start, units = "secs"))
  cat(sprintf("  split %3d/%3d  (%.1f min elapsed, ~%.1f s/split, ETA %.1f min)\n",
              k, K_SPLITS, el / 60, el / k, (el / k) * (K_SPLITS - k) / 60)); flush.console()
}
res <- do.call(rbind, all_rows[!vapply(all_rows, is.null, logical(1))])
total_min <- as.numeric(difftime(Sys.time(), t_start, units = "secs")) / 60

# ==============================================================================
# Summarize
# ==============================================================================
m_ <- function(x) mean(x, na.rm = TRUE); s_ <- function(x) sd(x, na.rm = TRUE)
cat(sprintf("\nStep 3: Summary across splits (%.1f min total)\n", total_min))

cat("\n=== [PRIMARY] EFFECT REPLICATION (both halves agree on signature direction) ===\n")
cat(sprintf("  %-20s %-9s %8s %8s %8s %8s %10s\n",
            "model", "sig", "fitΔ", "testΔ", "d_fit", "d_test", "RELIAB."))
for (mid in SEL_MODELS) {
  r <- res[res$model == mid, ]; if (nrow(r) == 0) next
  cat(sprintf("  %-20s %-9s %8.3f %8.3f %8.2f %8.2f %9.1f%%\n",
              mid, r$sig_param[1], m_(r$delta_fit), m_(r$delta_test),
              m_(r$d_fit), m_(r$d_test), 100 * m_(r$agree)))
}
cat("  (RELIAB. = % of splits where fit-half and test-half effects share sign;\n")
cat("   read together with d_test: high reliability + |d_test|>~0.2 = genuine replication.)\n")

cat("\n=== [SECONDARY] HELD-OUT PREDICTION (fit-half params scoring test-half) ===\n")
cat(sprintf("  %-20s %12s %12s %12s\n", "model", "NLL/trial", "NLL_HC", "NLL_DEP"))
for (mid in SEL_MODELS) {
  r <- res[res$model == mid, ]; if (nrow(r) == 0) next
  cat(sprintf("  %-20s %12.3f %12.3f %12.3f\n",
              mid, m_(r$oos_nll), m_(r$oos_nll_HC), m_(r$oos_nll_DEP)))
}
if ("4a" %in% SEL_MODELS && "accept_obs_HC" %in% names(res)) {
  r <- res[res$model == "4a", ]
  cat("\n  4a behavioral prediction (test half):\n")
  for (g in c("HC", "DEP")) {
    cat(sprintf("    %-3s accept obs=%.3f pred=%.3f | RT median obs=%.2fs pred=%.2fs\n", g,
                m_(r[[paste0("accept_obs_", g)]]), m_(r[[paste0("accept_pred_", g)]]),
                m_(r[[paste0("rt_obs_", g)]]),     m_(r[[paste0("rt_pred_", g)]])))
  }
  cat("    (note: 4a fixes bias=0.5 with no drift intercept -> accept-rate is\n")
  cat("     structurally ~0.5 regardless of fit; predictive validity rests on RT.)\n")
}

# ==============================================================================
# Save
# ==============================================================================
# outdir created before the K-loop; per_split.csv already written incrementally.
write.csv(res, per_split_path, row.names = FALSE)
summ <- do.call(rbind, lapply(SEL_MODELS, function(mid) {
  r <- res[res$model == mid, ]; if (nrow(r) == 0) return(NULL)
  data.frame(model = mid, sig_param = r$sig_param[1], K = nrow(r),
             delta_fit = m_(r$delta_fit), delta_test = m_(r$delta_test),
             d_fit = m_(r$d_fit), d_test = m_(r$d_test),
             reliability_agree = m_(r$agree),
             oos_nll = m_(r$oos_nll))
}))
write.csv(summ, file.path(outdir, "split_half_summary.csv"), row.names = FALSE)
cat(sprintf("\nResults saved to: %s\nDone!\n", outdir))
