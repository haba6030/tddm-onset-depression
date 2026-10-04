# =============================================================================
# Shared analytical tDDM engine  (single source of truth)
# =============================================================================
# Every fitting/recovery script should `source()` this instead of pasting its own
# copy.  Before 2026-08-07 the Navarro-Fuss density and the split-at-tau
# transition density were duplicated verbatim in 14 non-archive scripts, so a
# defect in one was a defect in all of them.
#
# Conventions (unchanged from the pipeline):
#   * sigma = 1.4 is hard-coded in the original C++ engine, so the analytical
#     likelihood works in the RESCALED space: drift' = drift / SIGMA,
#     boundary'  = 2 * thres / SIGMA, w = 0.5.  In that space the diffusion
#     coefficient is 1.
#   * dwiener_nf() returns the LOWER-boundary (y = 0) density.
#     upper -> dwiener_nf(t, -v, a, 1 - w) ;  lower -> dwiener_nf(t, v, a, w)
#   * acceptance (choice == 0) maps to the UPPER boundary, i.e. choice_ddm = 1.
#
# TRANSITION DENSITY (the 2026-08-07 fix)
# ---------------------------------------
# A trial whose decision time exceeds tau is split; the likelihood integrates
# over the accumulator position x at time tau:
#     L = INT p_transition(x, tau | x0, drift1) * f_wiener(dt - tau | drift2, x/a) dx
# p_transition must be the density of being at x at time tau AND not having been
# absorbed yet.  Two implementations:
#
#   "pipeline"  dnorm(x, x0 + v*t, SIGMA*sqrt(t))       <- as shipped; WRONG twice
#               (i)  uses SIGMA in the rescaled space, where the diffusion
#                    coefficient is 1  -> spread inflated by 40%
#               (ii) unconditional Gaussian: ignores absorption, so paths that
#                    already produced a response keep contributing
#   "image"     Girsanov x method of images, absorbing at 0 and a; a
#               SUB-probability density integrating to P(alive at tau)
#
# Verified against 150k-path simulation (scripts/tDDM/verify_transition_density.R):
#   tau=2.0, a=2.57 : true survival 0.272 | pipeline mass 0.464 (+70%) | image 0.265
#   tau=1.0, a=1.57 : true survival 0.177 | pipeline mass 0.413 (+133%) | image 0.165
#   shape L1 error  : pipeline 0.17-0.24  | image 0.015-0.027
#
# Keep "pipeline" available so historical results can be reproduced exactly.
# =============================================================================

TDDM_SIGMA   <- 1.4
TDDM_EPS     <- 1e-10
TDDM_PENALTY <- 1e12
TDDM_NGRID   <- 30      # used only when TDDM_QUAD == "grid"
# 2026-09-28: the split-at-tau integral is computed by adaptive quadrature with the
# two boundary layers split out. The old fixed 30-point grid on [0.01, bsep - 0.01]
# overstated the trial likelihood up to ~15x when dt - tau was tiny, and the optimizer
# parked tau at a trial's decision time (docs/PAPER/Quadrature_Validation_260927.md).
# Set TDDM_QUAD <- "grid" only to reproduce results produced before 2026-09-28.
TDDM_QUAD    <- "accurate"
TDDM_NIMAGE  <- 4

# ---- canonical settings (single source of truth) -----------------------------
# Every error found on 2026-08-10 had the same shape: a script's DEFAULT still
# pointed at a superseded setting and the caller did not override it.
#   E1  split_half_predictive.R kept nDT_UB = 1.5 after the suite moved to 3.0
#   E2  compute_tau_SE.R auto-picked the newest fits dir (a 20-subject pilot)
#   E3  end_to_end injected a truth taken from SHIPPED-density fits while
#       fitting with the corrected density
#   E7  the x-integration grid was 25 in three scripts and 30 in the others
# Scripts should read their defaults from here rather than hard-coding them.
TDDM_CANON <- list(
  ndt_ub    = 3.0,        # 1.5 censored 74% of subjects at the ceiling
  ngrid     = 30,         # x-grid for the split-at-tau integral (TDDM_QUAD = "grid" only)
  quad      = "accurate", # split-at-tau quadrature (2026-09-28)
  orientation = "mean",   # manuscript model: mean leads, contrast joins at tau
  rt_filter = c(0.15, 10),
  trans_default = "image",
  # tau operating point of each density -- negative/positive controls MUST
  # generate at the operating point of the density they are testing
  tau_operating_point = c(image = 0.720, pipeline = 1.800),
  # fits that correspond to each density (mean-first, N = 351, nDT <= 3.0)
  fits = c(image    = "results/tDDM/corrected_likelihood_20260807_222938/fits_4j_image.csv",
           pipeline = "results/tDDM/refit_gateB_20260624_172551/fits_4j_widenDT.csv")
)

# ---- Navarro-Fuss first-passage density -------------------------------------
# NOTE: the shipped version had `if (abs(term) < 1e-12) break` inside the
# large-time loop.  With w = 0.5 the k = 2 term is exactly 0 (sin(k*pi/2) = 0 for
# even k), so the break fired at k = 2 and truncated the series after one term,
# overstating the density by up to 5.8% near t = a^2/pi^2.  Impact on the fits
# was verified negligible (parameter r = 0.998-1.000, group p-values unchanged;
# results/tDDM/engine_break_bug_20260804_041524) because almost no trial lands in
# the affected window, but the break is removed here.
dwiener_large <- function(t, v, a, w, K = 7) {
  if (t <= 0 || a <= 0) return(0)
  prefix <- (pi / a^2) * exp(-v * a * w - v^2 * t / 2); s <- 0
  for (k in 1:K) s <- s + k * sin(k * pi * w) * exp(-k^2 * pi^2 * t / (2 * a^2))
  prefix * s
}
dwiener_small <- function(t, v, a, w, K = 7) {
  if (t <= 0 || a <= 0) return(0)
  prefix <- (a / sqrt(2 * pi * t^3)) * exp(-v * a * w - v^2 * t / 2); s <- 0
  for (k in (-K):K) { wk <- w + 2 * k; s <- s + wk * exp(-wk^2 * a^2 / (2 * t)) }
  prefix * s
}
dwiener_nf <- function(t, v, a, w) {
  if (t <= 0 || a <= 0) return(0)
  if (t < a^2 / (pi^2)) dwiener_small(t, v, a, w) else dwiener_large(t, v, a, w)
}

# ---- vectorised Wiener density over w -----------------------------------------
# The split-at-tau integral evaluates the phase-2 density at every grid point.
# t (= dt_obs - tau), v and a are IDENTICAL across grid points -- only w = x/a
# varies -- so the small/large-time branch is uniform and the whole grid can be
# done in one vectorised call instead of `ngrid` scalar ones.  Same arithmetic,
# verified by scripts/tDDM/test_engine_equivalence.R to < 1e-10.
dwiener_nf_vec <- function(t, v, a, w) {
  if (t <= 0 || a <= 0) return(rep(0, length(w)))
  if (t < a^2 / (pi^2)) {
    prefix <- (a / sqrt(2 * pi * t^3)) * exp(-v * a * w - v^2 * t / 2)
    s <- 0
    for (k in (-7):7) { wk <- w + 2 * k; s <- s + wk * exp(-wk^2 * a^2 / (2 * t)) }
  } else {
    prefix <- (pi / a^2) * exp(-v * a * w - v^2 * t / 2)
    s <- 0
    for (k in 1:7) s <- s + k * sin(k * pi * w) * exp(-k^2 * pi^2 * t / (2 * a^2))
  }
  prefix * s
}

# ---- transition densities ----------------------------------------------------
p_trans_pipeline <- function(x, t, x0, v, a) {
  if (t <= 0) return(rep(0, length(x)))
  dnorm(x, x0 + v * t, TDDM_SIGMA * sqrt(t))
}
p_trans_image <- function(x, t, x0, v, a, nmax = TDDM_NIMAGE) {
  if (t <= 0) return(rep(0, length(x)))
  s <- 0
  for (n in (-nmax):nmax)
    s <- s + exp(-(x - x0 - 2 * n * a)^2 / (2 * t)) - exp(-(x + x0 - 2 * n * a)^2 / (2 * t))
  s[s <= 0] <- 0
  # log-space prefactor so a large |v| cannot overflow before the exponential decay
  lp <- v * (x - x0) - v^2 * t / 2 - 0.5 * log(2 * pi * t)
  val <- exp(lp) * s
  val[!is.finite(val) | val < 0] <- 0
  val
}
get_transition <- function(which = c("image", "pipeline")) {
  switch(match.arg(which), image = p_trans_image, pipeline = p_trans_pipeline)
}

# ---- piecewise (two-phase) trial likelihood ---------------------------------
# tau <= 0            -> drift2 from t = 0, i.e. the static two-signal model
# dt_obs < tau        -> drift1 only; the phase-2-only weight is uninformed here
# otherwise           -> split at tau and integrate over the crossing state
# `trans` is deliberately REQUIRED (no default). It used to default to
# p_trans_image, which meant any caller that forgot the argument silently used
# the corrected density even in a --trans pipeline run -- a wrong answer rather
# than an error. Found 2026-08-10 when the Layer-1 harness called this kernel
# directly and reported a 4.1-unit mismatch in pipeline mode.
# Resolve the density once with get_transition(TRANS_MODE) and pass it through.
ll_piecewise_trial <- function(drift1, drift2, tau, bsep, w, nDT, obs_rt, obs_choice,
                               trans, ngrid = TDDM_NGRID, quad = TDDM_QUAD) {
  if (missing(trans)) stop("ll_piecewise_trial(): `trans` is required -- pass get_transition(\"image\"|\"pipeline\")")
  dt_obs <- obs_rt - nDT
  if (dt_obs <= 0.001) return(1e-12)
  if (dt_obs < tau || tau <= 0) {
    eff <- if (tau <= 0) drift2 else drift1
    d <- if (obs_choice == 1) dwiener_nf(dt_obs, -eff, bsep, 1 - w)
         else                 dwiener_nf(dt_obs,  eff, bsep, w)
    return(max(d, 1e-12))
  }
  rem <- dt_obs - tau
  p_from <- function(x) if (obs_choice == 1) dwiener_nf_vec(rem, -drift2, bsep, 1 - x / bsep)
                        else                 dwiener_nf_vec(rem,  drift2, bsep, x / bsep)
  if (quad == "grid") {                  # pre-2026-09-28 behaviour, kept for reproduction
    x_grid <- seq(0.01, bsep - 0.01, length.out = ngrid); dx <- diff(x_grid)[1]
    return(max(sum(trans(x_grid, tau, bsep * w, drift1, bsep) * p_from(x_grid)) * dx, 1e-12))
  }
  # accurate: adaptive quadrature with the two boundary layers (width ~ sqrt(rem) in
  # sigma-rescaled units) integrated separately from the interior. Validated against an
  # independent eigenfunction solution (scripts/tDDM/diag_quadrature_spectral.R).
  fx <- function(x) trans(x, tau, bsep * w, drift1, bsep) * p_from(x)
  h  <- min(bsep / 4, 20 * sqrt(rem))
  br <- unique(c(0, h, bsep - h, bsep))
  v <- tryCatch(sum(vapply(seq_len(length(br) - 1), function(k)
         integrate(fx, br[k], br[k + 1], subdivisions = 1000L, rel.tol = 1e-6)$value, 0)),
       error = function(e) {             # fallback: dense trapezoid including the endpoints
         xg <- seq(0, bsep, length.out = 4001); y <- fx(xg)
         sum((y[-1] + y[-length(y)]) / 2) * diff(xg)[1] })
  max(v, 1e-12)
}

# ---- model likelihoods -------------------------------------------------------
# Two-signal tDDM with a free onset.  `s1` leads; `s2` joins at tau.
# Canonical 4j is MEAN-first: s1 = mean_sig, s2 = contrast_sig (contrast-first is the
# orientation robustness check).
mk_ll_tau <- function(s1, s2, trans = p_trans_image) {
  use_cpp <- isTRUE(TDDM_CPP) && TDDM_QUAD == "accurate" && identical(trans, p_trans_image)
  if (use_cpp) return(function(params, td)          # C++ path: same algorithm, verified to 1e-10
    nll_tau_cpp(params, td[[s1]], td[[s2]], td$RT_ddm, as.integer(td$choice_ddm), TDDM_SIGMA, TDDM_EPS))
  function(params, td) {
    d1 <- params[1]; d2 <- params[2]; thres <- params[3]; nDT <- params[4]; tau <- params[5]
    bsep <- 2 * thres / TDDM_SIGMA; w <- 0.5
    if (bsep <= 0.01 || nDT <= 0) return(TDDM_PENALTY)
    nll <- 0
    for (i in 1:nrow(td)) {
      drift1 <- d1 * td[[s1]][i] / TDDM_SIGMA
      drift2 <- (d1 * td[[s1]][i] + d2 * td[[s2]][i]) / TDDM_SIGMA
      nll <- nll - log(max(ll_piecewise_trial(drift1, drift2, tau, bsep, w, nDT,
                                              td$RT_ddm[i], td$choice_ddm[i], trans),
                           TDDM_EPS))
    }
    if (!is.finite(nll)) TDDM_PENALTY else nll
  }
}
# Static two-signal model (tau-free).  Never touches the transition density.
mk_ll_static <- function(s1, s2) {
  function(params, td) {
    d1 <- params[1]; d2 <- params[2]; thres <- params[3]; nDT <- params[4]
    bsep <- 2 * thres / TDDM_SIGMA; w <- 0.5
    if (bsep <= 0.01 || nDT <= 0) return(TDDM_PENALTY)
    nll <- 0
    for (i in 1:nrow(td)) {
      drift <- (d1 * td[[s1]][i] + d2 * td[[s2]][i]) / TDDM_SIGMA
      dt <- td$RT_ddm[i] - nDT
      if (dt <= 0.001) { nll <- nll - log(TDDM_EPS); next }
      dens <- if (td$choice_ddm[i] == 1) dwiener_nf(dt, -drift, bsep, 1 - w)
              else                       dwiener_nf(dt,  drift, bsep, w)
      nll <- nll - log(max(dens, TDDM_EPS))
    }
    if (!is.finite(nll)) TDDM_PENALTY else nll
  }
}

# Single-signal static model (the manuscript's BASELINE, k = 3).  Only the
# leading signal drives drift; there is no contrast weight and no onset.  Added
# 2026-08-19: the model space needs it for the baseline -> contrast comparison
# (the surviving structural statement "adding the contrast weight without an
# onset buys nothing"), and every existing 4a fit is in results/archive/ from
# before the VTR fix, so none of them are usable.
mk_ll_static1 <- function(s1) {
  function(params, td) {
    d1 <- params[1]; thres <- params[2]; nDT <- params[3]
    bsep <- 2 * thres / TDDM_SIGMA; w <- 0.5
    if (bsep <= 0.01 || nDT <= 0) return(TDDM_PENALTY)
    nll <- 0
    for (i in 1:nrow(td)) {
      drift <- (d1 * td[[s1]][i]) / TDDM_SIGMA
      dt <- td$RT_ddm[i] - nDT
      if (dt <= 0.001) { nll <- nll - log(TDDM_EPS); next }
      dens <- if (td$choice_ddm[i] == 1) dwiener_nf(dt, -drift, bsep, 1 - w)
              else                       dwiener_nf(dt,  drift, bsep, w)
      nll <- nll - log(max(dens, TDDM_EPS))
    }
    if (!is.finite(nll)) TDDM_PENALTY else nll
  }
}

# ---- C++ likelihood (2026-09-28) -------------------------------------------------
# mk_ll_tau() uses nll_tau_cpp() from tddm_split.cpp when TDDM_CPP is TRUE, the
# quadrature is "accurate" and the transition density is the method of images. The C++
# code calls the same QUADPACK routine as R integrate() with the same breakpoints,
# tolerances and fallback; a full DEoptim fit is ~10x faster and returns the identical
# optimum (scripts/tDDM/verify_cpp_likelihood.R). Set TDDM_CPP <- FALSE before sourcing
# this file to force the R path.
if (!exists("TDDM_CPP")) TDDM_CPP <- TRUE
if (isTRUE(TDDM_CPP)) {
  .tddm_dir <- if (exists("MODEL_DIR")) file.path(MODEL_DIR, "tDDM") else file.path("models", "tDDM")
  TDDM_CPP <- tryCatch({
    Rcpp::sourceCpp(file.path(.tddm_dir, "tddm_split.cpp"), cacheDir = file.path(.tddm_dir, ".rcpp_cache"))
    TRUE
  }, error = function(e) { warning("tddm_split.cpp failed to compile; using the R likelihood: ",
                                   conditionMessage(e)); FALSE })
}

# ---- forward simulation of one two-phase trial --------------------------------
# Single source for every generative control (recovery, positive control, null-tau,
# PBCM). Unrescaled space: boundaries 0 and 2*th, start 2*th*w, noise sd `sg`.
# Added 2026-09-27. The per-script Euler copies (dt = 1 ms, a crossing counted
# only when the step ENDS outside) lengthened the mean decision time by 19-44 ms
# relative to the likelihood (scripts/tDDM/diag_simulator_vs_likelihood.R).
# Two corrections: (1) a Brownian-bridge test catches paths that cross and return
# within a step, P = exp(-2 d0 d1 / (sg^2 h)); (2) the step is cut at tau so the
# drift switches exactly at the onset. A hit is timed at the step midpoint.
tddm_simulate_trial <- function(d1, d2, th, nDT, tau, s1, s2, w = 0.5, sg = TDDM_SIGMA,
                                dt = 0.001, mx = 10) {
  a <- 2 * th; x <- a * w; t <- 0
  v1 <- d1 * s1; v2 <- d1 * s1 + d2 * s2
  while (t < mx) {
    ph1 <- tau > 0 && t < tau
    h <- if (ph1) min(dt, tau - t) else dt
    xn <- x + (if (ph1) v1 else v2) * h + sg * sqrt(h) * rnorm(1)
    if (xn >= a) return(list(rt = t + h / 2 + nDT, ch = 1))
    if (xn <= 0) return(list(rt = t + h / 2 + nDT, ch = 0))
    k <- 2 / (sg^2 * h)
    if (runif(1) < exp(-k * x * xn))             return(list(rt = t + h / 2 + nDT, ch = 0))
    if (runif(1) < exp(-k * (a - x) * (a - xn))) return(list(rt = t + h / 2 + nDT, ch = 1))
    x <- xn; t <- t + h
  }
  list(rt = mx + nDT, ch = as.integer(x > a / 2))
}

# Vectorised version of tddm_simulate_trial() for one participant's whole trial set
# (same corrections: bridge crossing test, step cut at tau, hit at the step midpoint).
# Returns RT (non-decision time added; paths alive at tmax are censored at tmax + ndt)
# and acc (1 = upper = acceptance, NA if censored).
tddm_simulate_block <- function(dm, dc, th, ndt, tau, ms, cs, sg = TDDM_SIGMA,
                                dt = 0.001, tmax = 10) {
  n <- length(ms); a <- 2 * th; x <- rep(a / 2, n)
  rt <- rep(NA_real_, n); acc <- rep(NA_integer_, n)
  v1 <- dm * ms; v2 <- dm * ms + dc * cs
  live <- rep(TRUE, n); t <- 0
  while (t < tmax && any(live)) {
    ph1 <- tau > 0 && t < tau
    h <- if (ph1) min(dt, tau - t) else dt
    i <- which(live); x0 <- x[i]
    xn <- x0 + (if (ph1) v1[i] else v2[i]) * h + sg * sqrt(h) * rnorm(length(i))
    k <- 2 / (sg^2 * h)
    up <- xn >= a; lo <- xn <= 0
    ins <- !(up | lo)
    lo[ins] <- runif(sum(ins)) < exp(-k * x0[ins] * xn[ins])
    ins <- !(up | lo)
    up[ins] <- runif(sum(ins)) < exp(-k * (a - x0[ins]) * (a - xn[ins]))
    if (any(up)) { rt[i[up]] <- t + h / 2 + ndt; acc[i[up]] <- 1L }
    if (any(lo)) { rt[i[lo]] <- t + h / 2 + ndt; acc[i[lo]] <- 0L }
    live[i[up | lo]] <- FALSE
    x[i] <- xn; t <- t + h
  }
  rt[is.na(rt)] <- tmax + ndt
  data.frame(RT = rt, acc = acc)
}

# ---- data preparation (identical filter/coding to the canonical pipeline) ----
tddm_load_study4 <- function(data_dir) {
  files <- list.files(data_dir, pattern = "^study4_.*[.]csv$", full.names = TRUE)
  d <- do.call(rbind, lapply(files, function(p) {
    df <- read.csv(p, stringsAsFactors = FALSE)
    df$group <- toupper(sub("[.]csv$", "", sub("^study4_", "", basename(p)))); df }))
  if ("sub" %in% names(d) && !"subject" %in% names(d))
    names(d)[names(d) == "sub"] <- "subject"
  d <- d[is.finite(d$choice_RT) & d$choice_RT > 0.15 & d$choice_RT <= 10, ]
  d$choice_ddm   <- 1 - d$choice                 # accept (choice == 0) -> upper
  d$RT_ddm       <- d$choice_RT
  d$mean_sig     <- (d$RPE_scale + d$EPE_X1_scale) / 2
  d$contrast_sig <- (d$RPE_scale - d$EPE_X1_scale) / 2
  d[complete.cases(d[, c("mean_sig", "contrast_sig", "RT_ddm", "choice_ddm")]), ]
}
