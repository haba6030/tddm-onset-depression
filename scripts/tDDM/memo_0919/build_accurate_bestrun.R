# Best-NLL run of the 10 accurate-integral 4j refits (diag_quadrature_refit_*), compared
# with the grid-quadrature canonical best-of-10. Seeds 20260821-30 match the grid set.
# Run from the repository root:  Rscript scripts/tDDM/memo_0919/build_accurate_bestrun.R
if (!file.exists("config.R")) stop("run from the repository root")
source("config.R"); source(file.path(MODEL_DIR, "tDDM", "tddm_engine.R"))
R <- file.path(RESULTS_DIR, "tDDM")
dirs <- sort(Sys.glob(file.path(R, "diag_quadrature_refit_2026*")))
reps <- lapply(dirs, function(d) read.csv(file.path(d, "fits_4j_accurate_rep01seed.csv")))
reps <- lapply(reps, function(r) r[order(r$subject), ])
stopifnot(all(vapply(reps, function(r) identical(r$subject, reps[[1]]$subject), TRUE)))
cat(sprintf("%d accurate-integral runs\n", length(reps)))

NLL <- sapply(reps, `[[`, "negLogLik"); bi <- apply(NLL, 1, which.min)
out <- reps[[1]][, c("subject", "group", "n_trials")]
for (p in c("d_mean", "d_contrast", "thres", "nDT", "tau", "negLogLik", "min_gap"))
  out[[p]] <- sapply(seq_along(bi), function(i) reps[[bi[i]]][[p]][i])
out$BIC <- 2 * out$negLogLik + 5 * log(out$n_trials)
out$best_run <- bi
out$nll_range <- apply(NLL, 1, max) - apply(NLL, 1, min)
od <- file.path(R, "canonical_4j_accurate_bestof10")
dir.create(od, showWarnings = FALSE)
write.csv(out, file.path(od, "fits_4j_accurate_best10.csv"), row.names = FALSE)

# grid-quadrature canonical, for comparison only; absent from the public repository
gf <- file.path(R, "canonical_4j_bestof10_20260923", "fits_4j_image_best10.csv")
has_grid <- file.exists(gf); if (has_grid) grid <- read.csv(gf)
c4 <- read.csv(file.path(R, "corrected_likelihood_20260819_005659", "fits_4c_static.csv"))
a4 <- read.csv(file.path(R, "corrected_likelihood_20260819_012638", "fits_4a_static1.csv"))
if (has_grid) { grid <- grid[order(grid$subject), ]; stopifnot(identical(out$subject, grid$subject)) }
c4 <- c4[order(c4$subject), ]; a4 <- a4[order(a4$subject), ]
stopifnot(identical(out$subject, c4$subject), identical(out$subject, a4$subject))

g <- out$group
wt <- function(x, lab) {
  a <- x[g == "DEP"]; b <- x[g == "HC"]; t <- t.test(a, b)
  n1 <- length(a); n2 <- length(b)
  sp <- sqrt(((n1 - 1) * var(a) + (n2 - 1) * var(b)) / (n1 + n2 - 2)); d <- (mean(a) - mean(b)) / sp
  se <- sqrt((n1 + n2) / (n1 * n2) + d^2 / (2 * (n1 + n2)))
  cat(sprintf("  %-22s HC %.3f DEP %.3f | DEP-HC %+.3f t(%.0f)=%.2f p=%.4f d=%.2f [%.2f, %.2f]\n", lab,
              mean(b), mean(a), mean(a) - mean(b), t$parameter, t$statistic, t$p.value, d, d - 1.96 * se, d + 1.96 * se))
}
sink(file.path(od, "summary.txt"), split = TRUE)
cat("== accurate integral, best of", length(reps), "runs ==\n")
for (p in c("d_mean", "d_contrast", "thres", "nDT", "tau")) wt(out[[p]], p)
if (has_grid) {
  cat("== grid canonical, best of 10 runs ==\n")
  wt(grid$tau, "tau (grid)")
  cat(sprintf("r(tau accurate, tau grid) = %.3f | tau<=0 accurate %.1f%% vs grid %.1f%%\n",
              cor(out$tau, grid$tau), 100 * mean(out$tau <= 0), 100 * mean(grid$tau <= 0)))
}
pos <- out$tau > 0
cat(sprintf("tau>0 within 1 ms of a trial DT: %.1f%%\n", 100 * mean(out$min_gap[pos] < 1e-3, na.rm = TRUE)))
p <- out$tau > 0; cat("  excluding tau<=0: "); gg <- g; g <- g[p]; wt(out$tau[p], "tau>0"); g <- gg
w <- wilcox.test(out$tau[g == "DEP"], out$tau[g == "HC"]); cat(sprintf("  Wilcoxon p = %.4f\n", w$p.value))
cat(sprintf("NLL range across runs: median %.3f\n", median(out$nll_range)))
for (cmp in list(list("4j - 4c", out$BIC - c4$BIC), list("4j - 4a", out$BIC - a4$BIC), list("4c - 4a", c4$BIC - a4$BIC)))
  cat(sprintf("dBIC(%s): mean %+.3f | <0 %.1f%% | <-2 %.1f%%\n", cmp[[1]], mean(cmp[[2]]),
              100 * mean(cmp[[2]] < 0), 100 * mean(cmp[[2]] < -2)))
sink()
cat("written:", file.path(od, "fits_4j_accurate_best10.csv"), "\n")
