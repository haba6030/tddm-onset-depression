# Does the simulator used by the generative controls produce the distribution that the
# accurate likelihood assigns? Compares the old per-script Euler copy (dt = 1 ms, no
# crossing correction) with tddm_simulate_trial() (bridge-corrected, step cut at tau).
# 2026-09-27.
#   Rscript scripts/tDDM/diag_simulator_vs_likelihood.R [--n 20000] [--cores 6]
suppressPackageStartupMessages(library(parallel))
if (!file.exists("config.R")) stop("run from the repository root")
source("config.R"); source(file.path(MODEL_DIR, "tDDM", "tddm_engine.R"))
args <- commandArgs(trailingOnly = TRUE)
ga <- function(f, d) if (f %in% args) args[which(args == f) + 1] else d
N <- as.integer(ga("--n", 20000)); NC <- as.integer(ga("--cores", 6))
src <- readLines("scripts/tDDM/diag_quadrature_spectral.R")
eval(parse(text = src[grep("^acc <- function", src):(grep("^grid <- function", src) - 1)])); fallback_n <- 0

old_euler <- function(d1, d2, th, nDT, tau, s1, s2, sg = TDDM_SIGMA, dt = .001, mx = 10) {
  a <- 2 * th; x <- a / 2; t <- 0; v1 <- d1 * s1; v2 <- d1 * s1 + d2 * s2
  while (t < mx) {
    v <- if (tau > 0 && t < tau) v1 else v2
    x <- x + v * dt + sg * sqrt(dt) * rnorm(1); t <- t + dt
    if (x >= a) return(list(rt = t + nDT, ch = 1)); if (x <= 0) return(list(rt = t + nDT, ch = 0))
  }
  list(rt = mx + nDT, ch = as.integer(x > a / 2))
}
cases <- data.frame(d1 = c(0.8, -0.5, 1.2, 0.3), d2 = c(1.0, 1.1, -1.6, 0.0),
                    thres = c(1.9, 2.1, 1.5, 2.6), tau = c(0.6, 0.9, 0.3, 0))
qs <- c(.1, .5, .9)
RNGkind("L'Ecuyer-CMRG"); set.seed(20260927)
for (k in seq_len(nrow(cases))) {
  p <- cases[k, ]; a <- 2 * p$thres / TDDM_SIGMA; v1 <- p$d1 / TDDM_SIGMA
  v2 <- (p$d1 + p$d2) / TDDM_SIGMA   # simulator phase-2 drift = d1*s1 + d2*s2 with s1 = s2 = 1
  tt <- seq(0.0005, 10, by = 0.001)
  lik <- function(t, ch) if (p$tau <= 0) (if (ch == 1) dwiener_nf(t, -v2, a, .5) else dwiener_nf(t, v2, a, .5)) else acc(v1, v2, p$tau, a, t, ch)
  f1 <- sapply(tt, lik, ch = 1); f0 <- sapply(tt, lik, ch = 0)
  # both sides conditioned on DT < 10 s (the simulator censors at 10 s)
  P1 <- sum(f1) / sum(f1 + f0); m <- sum(tt * (f1 + f0)) / sum(f1 + f0)
  F1 <- cumsum(f1) / sum(f1); q_lik <- sapply(qs, function(q) tt[which(F1 >= q)[1]])
  cat(sprintf("case %d (d1 %.1f, d2 %.1f, thres %.1f, tau %.1f) | likelihood P(upper) %.4f, mean DT %.4f, upper DT q10/50/90 %s\n",
              k, p$d1, p$d2, p$thres, p$tau, P1, m, paste(sprintf("%.3f", q_lik), collapse = "/")))
  for (lab in c("old Euler", "bridge")) {
    fn <- if (lab == "old Euler") old_euler else tddm_simulate_trial
    s <- do.call(rbind, mclapply(1:NC, function(c) do.call(rbind, lapply(seq_len(ceiling(N / NC)), function(i) {
      r <- fn(p$d1, p$d2, p$thres, 0, p$tau, 1, 1); c(r$rt, r$ch) })), mc.cores = NC))
    s <- s[s[, 1] < 10, , drop = FALSE]
    se <- sd(s[, 1]) / sqrt(nrow(s))
    q_sim <- quantile(s[s[, 2] == 1, 1], qs)
    cat(sprintf("  %-9s P(upper) %.4f (%+.4f, SE %.4f) | mean DT %+.1f ms (SE %.1f) | quantiles %s ms\n", lab,
                mean(s[, 2]), mean(s[, 2]) - P1, sqrt(P1 * (1 - P1) / nrow(s)), 1000 * (mean(s[, 1]) - m), 1000 * se,
                paste(sprintf("%+.0f", 1000 * (q_sim - q_lik)), collapse = "/")))
  }
}
