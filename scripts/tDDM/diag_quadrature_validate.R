# Validation of the split-at-tau trial likelihood, grid (engine, 30 points) vs accurate
# (adaptive quadrature with the boundary layers split out). 2026-09-27.
#  Test 1 (identity): with drift1 == drift2 the split integral must equal the unsplit
#          Wiener density (Chapman-Kolmogorov), known in closed form.
#  Test 2 (normalisation): with drift1 != drift2 the defective densities at the two
#          boundaries, integrated over all t, must sum to 1.
if (!file.exists("config.R")) stop("run from the repository root")
source("config.R"); source(file.path(MODEL_DIR, "tDDM", "tddm_engine.R"))
acc <- function(d1, d2, tau, bsep, dt, ch) {
  if (dt < tau) return(if (ch == 1) dwiener_nf(dt, -d1, bsep, .5) else dwiener_nf(dt, d1, bsep, .5))
  rem <- dt - tau
  fx <- function(x) p_trans_image(x, tau, bsep / 2, d1, bsep) *
    (if (ch == 1) dwiener_nf_vec(rem, -d2, bsep, 1 - x / bsep) else dwiener_nf_vec(rem, d2, bsep, x / bsep))
  h <- min(bsep / 4, 20 * sqrt(rem)); br <- unique(c(0, h, bsep - h, bsep))
  sum(vapply(seq_len(length(br) - 1), function(k)
    integrate(fx, br[k], br[k + 1], subdivisions = 2000L, rel.tol = 1e-8)$value, 0))
}
grid <- function(d1, d2, tau, bsep, dt, ch) ll_piecewise_trial(d1, d2, tau, bsep, .5, 0, dt, ch, p_trans_image)

set.seed(1)
cat("== Test 1: identity, 300 random configurations ==\n")
R1 <- t(replicate(300, {
  v <- runif(1, -2, 2); bsep <- runif(1, 0.8, 4); tau <- runif(1, 0.05, 2)
  rem <- 10^runif(1, -5, 0.3); ch <- rbinom(1, 1, .5)
  ex <- if (ch == 1) dwiener_nf(tau + rem, -v, bsep, .5) else dwiener_nf(tau + rem, v, bsep, .5)
  c(rem = rem, ex = ex, g = grid(v, v, tau, bsep, tau + rem, ch), a = acc(v, v, tau, bsep, tau + rem, ch))
}))
R1 <- as.data.frame(R1); R1 <- R1[R1$ex > 1e-6, ]
re <- function(x) abs(x - R1$ex) / R1$ex
for (lab in c("rem < 0.01 s", "rem >= 0.01 s")) {
  k <- if (lab == "rem < 0.01 s") R1$rem < .01 else R1$rem >= .01
  cat(sprintf("  %-14s n=%3d | accurate max rel err %.1e | grid median %.1e, max %.2f\n",
              lab, sum(k), max(re(R1$a)[k]), median(re(R1$g)[k]), max(re(R1$g)[k])))
}

cat("== Test 2: normalisation with drift1 != drift2, 12 random configurations ==\n")
tot <- function(f, d1, d2, tau, bsep) {
  # integrate over t on a grid that is dense just after tau, both boundaries
  tt <- sort(unique(c(seq(0.002, tau, length.out = 400), tau + 10^seq(-6, log10(12), length.out = 1500))))
  y <- sapply(tt, function(t) f(d1, d2, tau, bsep, t, 0) + f(d1, d2, tau, bsep, t, 1))
  sum((y[-1] + y[-length(y)]) / 2 * diff(tt))
}
set.seed(2)
R2 <- t(replicate(12, {
  d1 <- runif(1, -1.5, 1.5); d2 <- runif(1, -2, 2); bsep <- runif(1, 1, 3); tau <- runif(1, 0.1, 1.5)
  c(d1 = d1, d2 = d2, bsep = bsep, tau = tau, accurate = tot(acc, d1, d2, tau, bsep), grid = tot(grid, d1, d2, tau, bsep))
}))
print(round(R2, 4))
cat(sprintf("  total probability: accurate %.4f-%.4f | grid %.4f-%.4f (should be 1)\n",
            min(R2[, "accurate"]), max(R2[, "accurate"]), min(R2[, "grid"]), max(R2[, "grid"])))
