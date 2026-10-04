#!/usr/bin/env Rscript
# =============================================================================
# HC vs DEP group statistics for a tDDM fits file
# =============================================================================
# Reports, per parameter: group means, Welch t-test, Cohen's d (pooled SD) with
# a 95% CI, and the share of estimates sitting on a bound.  Written so the
# manuscript's baseline (4c) and main (4j) tables come from one code path.
#
#   Rscript scripts/tDDM/group_stats.R --fits <path/to/fits.csv> [--label 4c]
# =============================================================================
args <- commandArgs(trailingOnly = TRUE)
ga <- function(f, d) { if (f %in% args) { i <- which(args == f); if (i < length(args)) return(args[i+1]) }; d }
FITS  <- ga("--fits", "")
LABEL <- ga("--label", basename(dirname(FITS)))
if (FITS == "") stop("--fits is required")

f <- read.csv(FITS, stringsAsFactors = FALSE)
PARAMS <- intersect(c("d_mean", "d_contrast", "thres", "nDT", "tau"), names(f))
# k = 3 fits carry an all-NA d_contrast, and tau-free fits an all-NA tau; drop
# any column with no finite values rather than letting t.test() abort the run.
PARAMS <- PARAMS[vapply(PARAMS, function(p) sum(is.finite(f[[p]])) > 2, TRUE)]
# bounds used by refit_corrected_likelihood.R / fit_4c_static_full.R
BOUNDS <- list(d_mean = c(-3, 3), d_contrast = c(-3, 3), thres = c(0.5, 4),
               nDT = c(0.01, 3), tau = c(-2, 3))

cohen_d <- function(x, y) {
  nx <- length(x); ny <- length(y)
  s <- sqrt(((nx - 1) * var(x) + (ny - 1) * var(y)) / (nx + ny - 2))
  d <- (mean(x) - mean(y)) / s
  se <- sqrt((nx + ny) / (nx * ny) + d^2 / (2 * (nx + ny - 2)))
  c(d = d, lo = d - 1.96 * se, hi = d + 1.96 * se)
}

cat(sprintf("=== group statistics: %s ===\n", LABEL))
cat(sprintf("  %s\n", FITS))
cat(sprintf("  N = %d  (HC %d / DEP %d), median %d trials\n\n",
            nrow(f), sum(f$group == "HC"), sum(f$group == "DEP"), median(f$n_trials)))
cat(sprintf("  %-11s %9s %9s %9s %11s %8s %-18s %7s\n",
            "param", "HC", "DEP", "DEP-HC", "t", "p", "d [95% CI]", "%bound"))
for (p in PARAMS) {
  hc <- f[[p]][f$group == "HC"]; dp <- f[[p]][f$group == "DEP"]
  tt <- t.test(dp, hc)              # DEP - HC
  d  <- cohen_d(dp, hc)
  b  <- BOUNDS[[p]]
  pb <- if (is.null(b)) NA else 100 * mean(abs(f[[p]] - b[1]) < 1e-3 | abs(f[[p]] - b[2]) < 1e-3)
  cat(sprintf("  %-11s %9.3f %9.3f %+9.3f %11.3f %8.4f  %+.3f [%+.3f,%+.3f] %6.1f%%\n",
              p, mean(hc), mean(dp), mean(dp) - mean(hc), tt$statistic, tt$p.value,
              d["d"], d["lo"], d["hi"], pb))
}
cat("\n  (t and d are DEP - HC; positive d = DEP higher. %bound = estimates at a box constraint.)\n")
