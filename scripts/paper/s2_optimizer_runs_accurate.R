#!/usr/bin/env Rscript
# Supplement S2 numbers: agreement among the 10 accurate-integral runs of the main model
# (diag_quadrature_refit_*, seeds 20260821-30). tau <= 0 is set to 0 (single-phase model).
#   Rscript scripts/paper/s2_optimizer_runs_accurate.R
if (!file.exists("config.R")) stop("run from the repository root")
source("config.R")
R <- file.path(RESULTS_DIR, "tDDM")
dirs <- sort(Sys.glob(file.path(R, "diag_quadrature_refit_2026*")))
reps <- lapply(dirs, function(d) { r <- read.csv(file.path(d, "fits_4j_accurate_rep01seed.csv")); r[order(r$subject), ] })
stopifnot(length(reps) == 10, all(vapply(reps, function(r) identical(r$subject, reps[[1]]$subject), TRUE)))
NLL <- sapply(reps, `[[`, "negLogLik"); TAU <- pmax(sapply(reps, `[[`, "tau"), 0)
best <- apply(NLL, 1, min); near <- NLL - best <= 0.01
cat(sprintf("all 10 runs within 0.01 NLL of the best: %.1f%%\n", 100 * mean(rowSums(near) == 10)))
agree <- sapply(seq_len(nrow(NLL)), function(i) diff(range(TAU[i, near[i, ]])) <= 0.1)
cat(sprintf("onset within 0.1 s among runs within 0.01 of the best: %.1f%%\n", 100 * mean(agree)))
g <- reps[[1]]$group
per <- t(apply(TAU, 2, function(t) { x <- t.test(t[g == "DEP"], t[g == "HC"])
  c(dtau = mean(t[g == "DEP"]) - mean(t[g == "HC"]), p = x$p.value) }))
print(round(per, 4))
cat(sprintf("group onset difference across runs: %+.3f to %+.3f s, p from %.3f to %.3f\n",
            min(per[, "dtau"]), max(per[, "dtau"]), min(per[, "p"]), max(per[, "p"])))
