# Independent reference for the split-at-tau trial likelihood (2026-09-27).
#
# Both phases are solved by the sine-eigenfunction expansion of the Fokker-Planck
# equation (unit diffusion, absorbing at 0 and a). Phase 1 gives p(x, tau) as a sine
# series; phase 2 projects it onto the phase-2 eigenfunctions with a closed-form
# overlap integral I_jk, and the first-passage density is the boundary flux. The
# construction uses no method of images, no Navarro-Fuss series and no quadrature,
# so it is independent of both the engine grid and the accurate integral.
#
#   V0  spectral vs closed-form Wiener density (drift1 == drift2) and vs the 4 MC cases
#   V1  1000 random configurations, rem from 1e-7 to 3 s: grid30 / accurate vs spectral
#   V2  every phase-2 trial of all 351 participants at the grid-canonical parameters
#       and at the accurate-integral rep01 parameters
#   J1  grid refinement: does a finer fixed grid converge (ngrid 30 .. 3000)?
#
#   Rscript scripts/tDDM/diag_quadrature_spectral.R
if (!file.exists("config.R")) stop("run from the repository root")
source("config.R"); source(file.path(MODEL_DIR, "tDDM", "tddm_engine.R"))
R <- file.path(RESULTS_DIR, "tDDM")
OD <- file.path(R, paste0("diag_quadrature_spectral_", format(Sys.time(), "%Y%m%d_%H%M%S")))
dir.create(OD); sink(file.path(OD, "summary.txt"), split = TRUE)

# ---- spectral first-passage density after the drift change -------------------
# lower boundary = choice_ddm 0, upper = choice_ddm 1 (engine convention)
spec_dens <- function(v1, v2, tau, a, t, ch, x0 = a / 2) {
  if (t < tau) {                                   # single phase: pure sine series
    K <- min(60000, ceiling(a * sqrt(80 / t) / pi) + 50); k <- 1:K
    lam <- k^2 * pi^2 / (2 * a^2)
    s <- sum(k * sin(k * pi * x0 / a) * exp(-lam * t) * (if (ch == 1) (-1)^(k + 1) else 1))
    pre <- if (ch == 1) exp(v1 * (a - x0) - v1^2 * t / 2) else exp(-v1 * x0 - v1^2 * t / 2)
    return(pre * pi / a^2 * s)
  }
  s <- t - tau
  J <- min(4000, ceiling(a * sqrt(80 / tau) / pi) + 50); j <- 1:J
  K <- min(60000, ceiling(a * sqrt(80 / max(s, 1e-12)) / pi) + 50); k <- 1:K
  c <- v1 - v2; al <- j * pi / a; be <- k * pi / a
  sg <- outer(j, k, function(p, q) (-1)^(p - q))
  Cm <- function(g2, same) {                        # int_0^a e^{cy} cos(g y) dy, g*a = m*pi
    out <- c * (sg * exp(c * a) - 1) / (c^2 + g2)
    if (!is.null(same)) out[same] <- if (abs(c) < 1e-12) a else expm1(c * a) / c
    out
  }
  dmin <- outer(al, be, "-")^2; dplus <- outer(al, be, "+")^2
  same <- which(outer(j, k, "=="))
  if (abs(c) < 1e-12) { I <- matrix(0, J, K); I[same] <- a / 2 }
  else I <- 0.5 * (Cm(dmin, same) - Cm(dplus, NULL))
  wj <- sin(al * x0) * exp(-j^2 * pi^2 * tau / (2 * a^2))
  b <- (2 / a)^2 * exp(-v1 * x0 - v1^2 * tau / 2) * as.vector(wj %*% I)
  dk <- b * be * exp(-v2^2 * s / 2 - k^2 * pi^2 * s / (2 * a^2))
  if (ch == 0) 0.5 * sum(dk) else 0.5 * exp(v2 * a) * sum(dk * (-1)^(k + 1))
}

# ---- the two quadratures under test ------------------------------------------
fallback_n <- 0
acc <- function(v1, v2, tau, a, t, ch) {
  if (t < tau) return(if (ch == 1) dwiener_nf(t, -v1, a, .5) else dwiener_nf(t, v1, a, .5))
  rem <- t - tau
  fx <- function(x) p_trans_image(x, tau, a / 2, v1, a) *
    (if (ch == 1) dwiener_nf_vec(rem, -v2, a, 1 - x / a) else dwiener_nf_vec(rem, v2, a, x / a))
  h <- min(a / 4, 20 * sqrt(rem)); br <- unique(c(0, h, a - h, a))
  tryCatch(sum(vapply(seq_len(length(br) - 1), function(k)   # identical to diag_quadrature_refit.R
    integrate(fx, br[k], br[k + 1], subdivisions = 1000L, rel.tol = 1e-6)$value, 0)),
    error = function(e) { fallback_n <<- fallback_n + 1
      xg <- seq(0, a, length.out = 4001); y <- fx(xg); sum((y[-1] + y[-length(y)]) / 2) * diff(xg)[1] })
}
grid <- function(v1, v2, tau, a, t, ch, n = 30)
  ll_piecewise_trial(v1, v2, tau, a, .5, 0, t, ch, p_trans_image, ngrid = n)
rel <- function(x, ref) abs(x - ref) / ref

# ---- V0: sanity of the reference itself --------------------------------------
cat("== V0a: spectral vs closed-form Wiener density, drift1 == drift2 (200 configs) ==\n")
set.seed(10)
V0 <- t(replicate(200, {
  v <- runif(1, -3, 3); a <- runif(1, 0.7, 5.7); tau <- runif(1, 0.01, 2)
  rem <- 10^runif(1, -6, 0.5); ch <- rbinom(1, 1, .5)
  ex <- if (ch == 1) dwiener_nf(tau + rem, -v, a, .5) else dwiener_nf(tau + rem, v, a, .5)
  c(ex = ex, sp = spec_dens(v, v, tau, a, tau + rem, ch))
}))
k <- V0[, "ex"] > 1e-8
cat(sprintf("  max rel err %.1e (n = %d with density > 1e-8)\n", max(rel(V0[k, "sp"], V0[k, "ex"])), sum(k)))

cat("== V0b: spectral vs Monte Carlo (400k paths, 20 ms windows), 4 canonical trials ==\n")
mc <- read.csv(file.path(R, "diag_quadrature_refit_20260927_151946", "mc_cases.csv"))
mcv <- list(c(.1467, .1359, .1105), c(.1604, .2869, .6006), c(.1230, .1444, .1844), c(.9036, 1.0018, 1.1161))
for (i in seq_len(nrow(mc))) {
  cc <- mc[i, ]; wins <- list(c(-.02, 0), c(0, .02), c(.04, .06))
  sp <- sapply(wins, function(w) {                 # window average of the spectral density
    tt <- cc$dt + seq(w[1], w[2], length.out = 41); tt <- tt[tt > 0]
    mean(sapply(tt, function(t) spec_dens(cc$d1, cc$d2, cc$tau, cc$bsep, t, 0))) })
  cat(sprintf("  subject %3d | MC %s | spectral %s | grid30 at dt %.3f\n", cc$subject,
              paste(sprintf("%.3f", mcv[[i]]), collapse = " "), paste(sprintf("%.3f", sp), collapse = " "), cc$grid30))
}

# ---- V1: random configurations ------------------------------------------------
cat("== V1: 1000 random configurations, drift1 != drift2, rem 1e-7..3 s ==\n")
set.seed(11)
V1 <- as.data.frame(t(replicate(1000, {
  v1 <- runif(1, -4, 4); v2 <- runif(1, -8, 8); a <- runif(1, 0.72, 5.7); tau <- 10^runif(1, -2, 0.5)
  rem <- 10^runif(1, -7, 0.5); ch <- rbinom(1, 1, .5); t <- tau + rem
  c(v1 = v1, v2 = v2, a = a, tau = tau, rem = rem, ch = ch,
    sp = spec_dens(v1, v2, tau, a, t, ch), acc = acc(v1, v2, tau, a, t, ch), g30 = grid(v1, v2, tau, a, t, ch))
})))
V1 <- V1[V1$sp > 1e-6, ]
bands <- cut(V1$rem, c(0, 1e-5, 1e-4, 1e-3, 1e-2, 0.1, Inf),
             labels = c("<10us", "10-100us", "0.1-1ms", "1-10ms", "10-100ms", ">100ms"))
tab <- do.call(rbind, lapply(split(V1, bands), function(d) data.frame(
  n = nrow(d), acc_median = median(rel(d$acc, d$sp)), acc_max = max(rel(d$acc, d$sp)),
  grid_median_ratio = median(d$g30 / d$sp), grid_max_ratio = max(d$g30 / d$sp), grid_min_ratio = min(d$g30 / d$sp))))
print(signif(tab, 3))
cat(sprintf("  accurate integral: fallback used %d times\n", fallback_n))
write.csv(V1, file.path(OD, "V1_random_configs.csv"), row.names = FALSE)

# ---- V2: every phase-2 trial of the real data ---------------------------------
data <- tddm_load_study4(DATA_DIR)
on_data <- function(fits, lab) {
  rows <- list()
  for (i in seq_len(nrow(fits))) {
    f <- fits[i, ]; td <- data[data$subject == f$subject, ]
    a <- 2 * f$thres / TDDM_SIGMA
    for (r in seq_len(nrow(td))) {
      t <- td$RT_ddm[r] - f$nDT
      if (f$tau <= 0 || t <= 0.001 || t < f$tau) next
      v1 <- f$d_mean * td$mean_sig[r] / TDDM_SIGMA
      v2 <- (f$d_mean * td$mean_sig[r] + f$d_contrast * td$contrast_sig[r]) / TDDM_SIGMA
      ch <- td$choice_ddm[r]
      rows[[length(rows) + 1]] <- data.frame(subject = f$subject, rem = t - f$tau,
        sp = spec_dens(v1, v2, f$tau, a, t, ch), acc = acc(v1, v2, f$tau, a, t, ch), g30 = grid(v1, v2, f$tau, a, t, ch))
    }
  }
  d <- do.call(rbind, rows)
  d$l_acc <- log(pmax(d$acc, TDDM_EPS)) - log(pmax(d$sp, TDDM_EPS))
  d$l_g30 <- log(pmax(d$g30, TDDM_EPS)) - log(pmax(d$sp, TDDM_EPS))
  ps <- aggregate(cbind(l_acc, l_g30) ~ subject, d, sum)
  mx <- aggregate(l_g30 ~ subject, d, max)
  cat(sprintf("== V2 [%s]: %d phase-2 trials, %d participants ==\n", lab, nrow(d), nrow(ps)))
  cat(sprintf("  per trial  |log acc - log spectral|: median %.1e, max %.1e\n", median(abs(d$l_acc)), max(abs(d$l_acc))))
  cat(sprintf("  per trial  log grid30 - log spectral: median %+.4f, max %+.3f (ratio %.2f)\n",
              median(d$l_g30), max(d$l_g30), exp(max(d$l_g30))))
  cat(sprintf("  per participant NLL(true) - NLL(grid30): mean %.3f, median %.3f | accurate: max |diff| %.1e\n",
              mean(ps$l_g30), median(ps$l_g30), max(abs(ps$l_acc))))
  cat(sprintf("  share of that gap carried by the single largest trial: median %.0f%%\n",
              100 * median(mx$l_g30 / ps$l_g30[match(mx$subject, ps$subject)])))
  write.csv(d, file.path(OD, sprintf("V2_trials_%s.csv", lab)), row.names = FALSE)
  invisible(d)
}
gf <- file.path(R, "canonical_4j_bestof10_20260923", "fits_4j_image_best10.csv")  # absent from the public repository
fa <- read.csv(file.path(R, "diag_quadrature_refit_20260927_151946", "fits_4j_accurate_rep01seed.csv"))
if (file.exists(gf)) { fallback_n <- 0; on_data(read.csv(gf), "grid_canonical"); cat(sprintf("  fallback used %d times\n", fallback_n)) }
fallback_n <- 0; on_data(fa, "accurate_rep01"); cat(sprintf("  fallback used %d times\n", fallback_n))

# ---- J1: does a finer fixed grid converge? ------------------------------------
cat("== J1: worst grid/true ratio over rem in [1e-7, 0.1] s, 50 random (v1, v2, a, tau) ==\n")
set.seed(12)
cfg <- replicate(50, c(v1 = runif(1, -3, 3), v2 = runif(1, -6, 6), a = runif(1, 1, 5), tau = runif(1, .1, 2)))
rems <- 10^seq(-7, -1, length.out = 61)
for (n in c(30, 100, 300, 1000, 3000)) {
  worst <- apply(cfg, 2, function(p) {
    sp <- sapply(rems, function(r) spec_dens(p["v1"], p["v2"], p["tau"], p["a"], p["tau"] + r, 0))
    g <- sapply(rems, function(r) grid(p["v1"], p["v2"], p["tau"], p["a"], p["tau"] + r, 0, n))
    k <- sp > 1e-6; c(max(g[k] / sp[k]), rems[k][which.max(g[k] / sp[k])]) })
  cat(sprintf("  ngrid %4d | worst ratio median %.2f, max %.2f | at rem median %.1e s\n",
              n, median(worst[1, ]), max(worst[1, ]), median(worst[2, ])))
}
sink(); cat("written:", OD, "\n")
