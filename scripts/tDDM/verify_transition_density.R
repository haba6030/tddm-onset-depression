#!/usr/bin/env Rscript
# =============================================================================
# Is the transition density formula correct?  Check it against simulation.
# =============================================================================
# The piecewise likelihood needs p(x at time tau, path not yet absorbed).
# Three candidates:
#   pipeline  dnorm(x, x0 + v*tau, SIGMA*sqrt(tau))     <- as shipped
#   gauss     dnorm(x, x0 + v*tau, sqrt(tau))           <- sigma scale fixed
#   image     Girsanov x method of images                <- absorbing boundaries
#
# Ground truth: simulate the rescaled diffusion (unit variance, boundaries 0/a),
# keep only paths still alive at tau, and histogram their positions.  A correct
# formula must match BOTH the shape AND the total mass (= survival probability).
#
# Usage: Rscript scripts/tDDM/verify_transition_density.R [--nsim 200000]
# =============================================================================
args <- commandArgs(trailingOnly = TRUE)
get_arg <- function(f, d) { if (f %in% args) { i <- which(args == f); if (i < length(args)) return(args[i+1]) }; d }
NSIM <- as.integer(get_arg("--nsim", 200000))
script_dir <- tryCatch(dirname(sys.frame(1)$ofile), error = function(e) {
  ca <- commandArgs(trailingOnly = FALSE); fa <- sub("^--file=", "", ca[grep("^--file=", ca)])
  if (length(fa) > 0) dirname(fa) else getwd() })
source(file.path(script_dir, "..", "..", "config.R"))
SIGMA <- 1.4

trans_pipeline <- function(x, t, x0, v, a) dnorm(x, x0 + v * t, SIGMA * sqrt(t))
trans_gauss    <- function(x, t, x0, v, a) dnorm(x, x0 + v * t, sqrt(t))
trans_image <- function(x, t, x0, v, a, nmax = 4) {
  pref <- exp(v * (x - x0) - v^2 * t / 2) / sqrt(2 * pi * t)
  s <- 0
  for (n in (-nmax):nmax)
    s <- s + exp(-(x - x0 - 2 * n * a)^2 / (2 * t)) - exp(-(x + x0 - 2 * n * a)^2 / (2 * t))
  pmax(pref * s, 0)
}

# ground truth: Euler simulation of the RESCALED process (unit diffusion)
sim_alive <- function(t_target, x0, v, a, nsim = NSIM, dt = 0.0005) {
  x <- rep(x0, nsim); alive <- rep(TRUE, nsim)
  nstep <- round(t_target / dt)
  sd_step <- sqrt(dt)
  for (k in seq_len(nstep)) {
    idx <- which(alive)
    if (!length(idx)) break
    x[idx] <- x[idx] + v * dt + rnorm(length(idx), 0, sd_step)
    alive[idx][x[idx] <= 0 | x[idx] >= a] <- FALSE
  }
  x[alive]
}

# representative operating points from the canonical fits:
# thres 1.1-2.7 -> a = 2*thres/SIGMA in 1.57-3.86 ; drift/SIGMA typically |.|<1
cases <- data.frame(
  a   = c(2.57, 2.57, 2.57, 3.86, 1.57),
  tau = c(0.50, 1.00, 2.00, 1.00, 1.00),
  v   = c(0.30, 0.30, 0.30, 0.30, 0.30))

cat("=== transition density vs simulation (rescaled space, unit diffusion) ===\n")
cat(sprintf("  %d paths per case, Euler dt=5e-4\n\n", NSIM))
cat(sprintf("%5s %5s %6s | %9s %9s %9s %9s | %s\n",
            "a", "tau", "v", "P(alive)", "mass_pipe", "mass_gauss", "mass_image", "shape error (L1)"))

out <- list()
for (i in seq_len(nrow(cases))) {
  a <- cases$a[i]; tau <- cases$tau[i]; v <- cases$v[i]; x0 <- a / 2
  xs <- sim_alive(tau, x0, v, a)
  p_alive <- length(xs) / NSIM
  grid <- seq(0.01, a - 0.01, length.out = 200); dx <- diff(grid)[1]
  # empirical density on the same grid, normalised so it integrates to P(alive)
  edges <- c(grid - dx / 2, grid[length(grid)] + dx / 2)
  h <- as.numeric(table(cut(xs, breaks = edges, include.lowest = TRUE)))
  emp <- h / NSIM / dx
  dens <- list(pipe  = trans_pipeline(grid, tau, x0, v, a),
               gauss = trans_gauss(grid, tau, x0, v, a),
               image = trans_image(grid, tau, x0, v, a))
  mass <- vapply(dens, function(z) sum(z) * dx, 0)
  l1   <- vapply(dens, function(z) sum(abs(z - emp)) * dx, 0)
  cat(sprintf("%5.2f %5.2f %6.2f | %9.3f %9.3f %9.3f %9.3f | pipe %.3f  gauss %.3f  image %.3f\n",
              a, tau, v, p_alive, mass[["pipe"]], mass[["gauss"]], mass[["image"]],
              l1[["pipe"]], l1[["gauss"]], l1[["image"]]))
  out[[i]] <- data.frame(a = a, tau = tau, v = v, p_alive = p_alive,
                         mass_pipeline = mass[["pipe"]], mass_gauss = mass[["gauss"]],
                         mass_image = mass[["image"]], l1_pipeline = l1[["pipe"]],
                         l1_gauss = l1[["gauss"]], l1_image = l1[["image"]])
}
res <- do.call(rbind, out)
ts <- format(Sys.time(), "%Y%m%d_%H%M%S")
outdir <- file.path(RESULTS_DIR, "tDDM", sprintf("transition_density_check_%s", ts))
dir.create(outdir, recursive = TRUE)
# 2026-08-20: record the simulation settings alongside the numbers.  The
# 2026-08-07 run wrote only the comparison table, so the path count it used
# could not be recovered afterwards and the manuscript could not state it.
res$nsim <- NSIM; res$dt <- 5e-4
write.csv(res, file.path(outdir, "density_check.csv"), row.names = FALSE)
write.csv(data.frame(nsim = NSIM, euler_dt = 5e-4, n_cases = nrow(res),
                     mass_err_pct_min = min(abs(100 * (res$mass_image - res$p_alive) / res$p_alive)),
                     mass_err_pct_max = max(abs(100 * (res$mass_image - res$p_alive) / res$p_alive)),
                     l1_min = min(res$l1_image), l1_max = max(res$l1_image),
                     run_utc = format(Sys.time(), "%Y-%m-%d %H:%M:%S", tz = "UTC")),
          file.path(outdir, "settings_and_summary.csv"), row.names = FALSE)

cat("\n  P(alive)   = simulated survival probability at tau — the mass a correct density must have\n")
cat("  mass_*     = integral of each formula over (0, a)\n")
cat("  shape L1   = integral |formula - empirical| ; 0 = perfect\n")
cat(sprintf("\nOutput: %s\n", outdir))
