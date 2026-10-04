#!/usr/bin/env Rscript
# =============================================================================
# Re-establish the tau model under a CORRECTED transition density
# =============================================================================
# The paper's framing rests on the BIC advantage of the tau model.  On the
# shipped likelihood that advantage is large (dBIC = -12.1/subject, 94% of 351
# subjects at dBIC < -2).  But the split-at-tau branch -- the branch only the tau
# model can use -- carries an inflated density (mass +70% at tau = 2 s, +133% at
# a low threshold; models/tDDM/tddm_engine.R).  An inflated density mechanically
# lowers NLL, so the BIC advantage has to be re-established with the corrected
# transition density before it can be defended in review.
#
# Fits, on identical data with a converged optimiser (NO VTR):
#   4c  static two-signal, K = 4          -- never touches the transition density
#   4j  contrast-first tau model, K = 5   -- SHIPPED transition density
#   4j  contrast-first tau model, K = 5   -- IMAGE (corrected) transition density
#
# Reports dBIC(4j - 4c) under both densities, plus the group contrasts.
#
# Usage:
#   Rscript scripts/tDDM/refit_corrected_likelihood.R --smoke      # 12 subjects
#   Rscript scripts/tDDM/refit_corrected_likelihood.R --cores 10
# =============================================================================
suppressPackageStartupMessages({ library(DEoptim); library(parallel) })
args <- commandArgs(trailingOnly = TRUE)
get_arg <- function(f, d) { if (f %in% args) { i <- which(args == f); if (i < length(args)) return(args[i+1]) }; d }
SMOKE   <- any(args == "--smoke")
N_CORES <- as.integer(get_arg("--cores", max(1, detectCores() - 2)))
ITERMAX <- as.integer(get_arg("--itermax", 200))
NP      <- as.integer(get_arg("--np", 80))
# 2026-08-07: the paper's narrative ("DEP integrates the reward-emotion contrast
# later") is MEAN-FIRST -- the common signal leads and the contrast joins at tau.
# The negative controls (nulltau_generative_sim.R, positive_control_tau.R) were
# already built mean-first; the file named fits_4j_contrast_first.csv is the odd
# one out.  Default to mean-first so the whole chain matches the manuscript.
ORIENT  <- match.arg(get_arg("--orientation", "mean"), c("mean", "contrast"))
ARMS    <- strsplit(get_arg("--arms", "4c_static,4j_pipeline,4j_image"), ",")[[1]]

script_dir <- tryCatch(dirname(sys.frame(1)$ofile), error = function(e) {
  ca <- commandArgs(trailingOnly = FALSE); fa <- sub("^--file=", "", ca[grep("^--file=", ca)])
  if (length(fa) > 0) dirname(fa) else getwd() })
source(file.path(script_dir, "..", "..", "config.R"))
source(file.path(MODEL_DIR, "tDDM", "tddm_engine.R"))

data <- tddm_load_study4(DATA_DIR)
subjects <- unique(data$subject)
if (SMOKE) subjects <- subjects[1:12]
cat(sprintf("=== refit under corrected transition density ===\n"))
cat(sprintf("  %d subjects, %d trials, itermax=%d NP=%d, cores=%d\n\n",
            length(subjects), sum(data$subject %in% subjects), ITERMAX, NP, N_CORES))

# mean-first : s1 = mean_sig leads, s2 = contrast_sig joins at tau  (paper model)
# contrast-first: the reverse
S1 <- if (ORIENT == "mean") "mean_sig" else "contrast_sig"
S2 <- if (ORIENT == "mean") "contrast_sig" else "mean_sig"
cat(sprintf("  orientation = %s-first  (s1 = %s leads, s2 = %s joins at tau)\n", ORIENT, S1, S2))
cat(sprintf("  arms = %s\n\n", paste(ARMS, collapse = ", ")))

MODELS <- list(
  `4a_static1`  = list(ll = mk_ll_static1(S1),                           K = 3,
                       lo = c(-3, 0.5, 0.01),        up = c(3, 4, 3)),
  `4c_static`   = list(ll = mk_ll_static(S1, S2),                        K = 4,
                       lo = c(-3, -3, 0.5, 0.01),     up = c(3, 3, 4, 3)),
  `4j_pipeline` = list(ll = mk_ll_tau(S1, S2, trans = p_trans_pipeline), K = 5,
                       lo = c(-3, -3, 0.5, 0.01, -2), up = c(3, 3, 4, 3, 3)),
  `4j_image`    = list(ll = mk_ll_tau(S1, S2, trans = p_trans_image),    K = 5,
                       lo = c(-3, -3, 0.5, 0.01, -2), up = c(3, 3, 4, 3, 3)))
MODELS <- MODELS[ARMS]

ctrl <- DEoptim.control(itermax = ITERMAX, NP = NP, trace = FALSE)   # NEVER set VTR

fit_one <- function(s, m) {
  td <- data[data$subject == s, ]
  if (nrow(td) < 5) return(NULL)
  f <- tryCatch(DEoptim(m$ll, lower = m$lo, upper = m$up, control = ctrl, td = td),
                error = function(e) NULL)
  if (is.null(f)) return(NULL)
  e <- f$optim$bestmem; n <- nrow(td)
  # e[1] is the weight on the LEADING signal, e[2] on the one that joins at tau
  # k = 3 (baseline) has no join weight, so its parameter vector is shifted
  one_sig <- m$K == 3
  d_join  <- if (one_sig) NA_real_ else e[2]
  thres_  <- if (one_sig) e[2] else e[3]
  nDT_    <- if (one_sig) e[3] else e[4]
  out <- data.frame(subject = s, group = td$group[1], n_trials = n,
                    negLogLik = f$optim$bestval,
                    BIC = 2 * f$optim$bestval + m$K * log(n),
                    de_iter = f$optim$iter,
                    d_lead = e[1], d_join = d_join, thres = thres_, nDT = nDT_,
                    stringsAsFactors = FALSE)
  out$d_mean     <- if (ORIENT == "mean") e[1] else d_join
  out$d_contrast <- if (ORIENT == "mean") d_join else e[1]
  out$tau <- if (m$K == 5) e[5] else NA_real_
  out
}

fits <- list(); t0 <- Sys.time()
for (nm in names(MODELS)) {
  res <- mclapply(subjects, fit_one, m = MODELS[[nm]], mc.cores = N_CORES)
  res <- Filter(is.data.frame, res)
  fits[[nm]] <- do.call(rbind, res)
  cat(sprintf("  [%-12s] %d/%d fits, %.1f min elapsed\n", nm, length(res), length(subjects),
              as.numeric(difftime(Sys.time(), t0, units = "mins"))))
  flush.console()
}
stopifnot(all(vapply(fits, function(f) min(f$de_iter) > 0, TRUE)))   # VTR guard

ts <- format(Sys.time(), "%Y%m%d_%H%M%S")
outdir <- file.path(RESULTS_DIR, "tDDM", sprintf("corrected_likelihood_%s", ts))
dir.create(outdir, recursive = TRUE)
for (nm in names(fits)) write.csv(fits[[nm]], file.path(outdir, sprintf("fits_%s.csv", nm)), row.names = FALSE)

# ------------------------------------------------------- BIC comparison -----
# 4c may come from this run or from an earlier one (the static model is
# orientation-invariant: drift = d1*s1 + d2*s2 is symmetric in the two signals,
# so only the parameter LABELS swap).  Same for the shipped-density 4j arm, which
# refit_gateB_sensitivity.R already produced mean-first with identical settings
# (itermax 200, NP 80, N_GRID 30, nDT bound 3.0, no VTR).
REUSE_4C <- get_arg("--reuse-4c", "")
REUSE_PIPE <- get_arg("--reuse-pipeline", "")
if (nzchar(REUSE_4C) && is.null(fits[["4c_static"]])) {
  z <- read.csv(REUSE_4C); fits[["4c_static"]] <- z
  cat(sprintf("  reused 4c fits: %s (%d rows)\n", REUSE_4C, nrow(z)))
}
if (nzchar(REUSE_PIPE) && is.null(fits[["4j_pipeline"]])) {
  z <- read.csv(REUSE_PIPE)
  if (!"BIC" %in% names(z)) z$BIC <- 2 * z$negLogLik + 5 * log(z$n_trials)
  fits[["4j_pipeline"]] <- z
  cat(sprintf("  reused 4j shipped-density fits: %s (%d rows)\n", REUSE_PIPE, nrow(z)))
}

cmp <- NULL
if (!is.null(fits[["4c_static"]])) {
  cmp <- fits[["4c_static"]][, c("subject", "group", "n_trials", "BIC", "negLogLik")]
  names(cmp)[4:5] <- c("BIC_4c", "NLL_4c")
  for (nm in intersect(c("4j_pipeline", "4j_image"), names(fits))) {
    x <- fits[[nm]][, c("subject", "BIC", "negLogLik", "tau")]
    names(x)[2:4] <- paste0(c("BIC_", "NLL_", "tau_"), sub("4j_", "", nm))
    cmp <- merge(cmp, x, by = "subject")
  }
  write.csv(cmp, file.path(outdir, "bic_comparison.csv"), row.names = FALSE)
}

cat("\n=========================================================================\n")
cat("  DOES THE TAU MODEL STILL WIN ON BIC UNDER THE CORRECTED DENSITY?\n")
cat("=========================================================================\n\n")
cat(sprintf("%-24s %10s %10s %10s %10s\n", "transition density", "mean dBIC", "sum dBIC", "% dBIC<-2", "mean dNLL"))
for (tag in c("pipeline", "image")) {
  if (is.null(cmp) || !paste0("BIC_", tag) %in% names(cmp)) next
  d <- cmp[[paste0("BIC_", tag)]] - cmp$BIC_4c
  n <- cmp[[paste0("NLL_", tag)]] - cmp$NLL_4c
  cat(sprintf("%-24s %+10.2f %+10.0f %9.0f%% %+10.3f\n",
              if (tag == "pipeline") "shipped (inflated)" else "image (corrected)",
              mean(d), sum(d), 100 * mean(d < -2), mean(n)))
}
if (!is.null(cmp)) cat(sprintf("\n  (negative dBIC favours the tau model; BIC penalty for K=5 vs 4 is log(n) = %.2f)\n",
            log(median(cmp$n_trials))))

# ------------------------------------------------------- group contrasts ----
cat("\n--- group contrasts (HC - DEP) ---\n")
cat(sprintf("%-14s %-12s %8s %8s %8s %8s\n", "model", "param", "HC", "DEP", "t", "p"))
for (nm in names(fits)) {
  f <- fits[[nm]]
  for (p in c("d_contrast", "d_mean", "thres", "nDT", "tau")) {
    v <- f[[p]]; if (all(is.na(v))) next
    hc <- v[f$group == "HC"]; dp <- v[f$group == "DEP"]
    if (length(hc) < 3 || length(dp) < 3) {
      cat(sprintf("%-14s %-12s %8.3f %8.3f %8s %8s\n", nm, p,
                  if (length(hc)) mean(hc) else NA, if (length(dp)) mean(dp) else NA, "-", "-"))
      next
    }
    tt <- t.test(hc, dp)
    cat(sprintf("%-14s %-12s %8.3f %8.3f %+8.2f %8.4f\n", nm, p, mean(hc), mean(dp),
                tt$statistic, tt$p.value))
  }
}

cat("\n--- tau distribution ---\n")
for (tag in c("pipeline", "image")) {
  if (is.null(cmp) || !paste0("tau_", tag) %in% names(cmp)) next
  v <- cmp[[paste0("tau_", tag)]]
  cat(sprintf("  %-10s mean %.3f  sd %.3f  median %.3f  %% <= 0: %.0f%%  %% at bound 3: %.0f%%\n",
              tag, mean(v), sd(v), median(v), 100 * mean(v <= 0), 100 * mean(v >= 2.99)))
}
cat(sprintf("\nOutput: %s\n", outdir))
