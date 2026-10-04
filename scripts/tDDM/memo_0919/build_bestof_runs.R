#!/usr/bin/env Rscript
# Canonical fits = for each participant, the run with the lowest negative log-likelihood
# among repeated fits. Works on the output folder of optimizer_reproducibility.R
# (fits_rep01.csv ... fits_rep10.csv) for either signal ordering. tau is stored as fitted;
# analyses set tau <= 0 to 0 (every such value gives the contrast model).
#   Rscript scripts/tDDM/memo_0919/build_bestof_runs.R --runs-dir <dir> --out <file.csv>
args <- commandArgs(trailingOnly = TRUE)
ga <- function(f) { i <- which(args == f); if (length(i) && i < length(args)) args[i + 1] else stop("missing ", f) }
D <- ga("--runs-dir"); OUT <- ga("--out")
files <- sort(Sys.glob(file.path(D, "fits_rep*.csv"))); stopifnot(length(files) > 0)
reps <- lapply(files, function(p) { r <- read.csv(p); r[order(r$subject), ] })
stopifnot(all(vapply(reps, function(r) identical(r$subject, reps[[1]]$subject), TRUE)))
NLL <- sapply(reps, `[[`, "negLogLik"); bi <- apply(NLL, 1, which.min)
out <- reps[[1]][, c("subject", "group", "n_trials")]
for (p in c("d_mean", "d_contrast", "thres", "nDT", "tau", "negLogLik"))
  out[[p]] <- sapply(seq_along(bi), function(i) reps[[bi[i]]][[p]][i])
out$BIC <- 2 * out$negLogLik + 5 * log(out$n_trials)
out$best_run <- bi
out$nll_range <- apply(NLL, 1, max) - apply(NLL, 1, min)
dir.create(dirname(OUT), showWarnings = FALSE, recursive = TRUE)
write.csv(out, OUT, row.names = FALSE)
g <- out$group; tau <- pmax(out$tau, 0); t <- t.test(tau[g == "DEP"], tau[g == "HC"])
cat(sprintf("%d runs | N = %d | tau (<= 0 set to 0) HC %.3f DEP %.3f | dtau %+.3f p = %.4f | median NLL range %.4f\nwritten: %s\n",
            length(reps), nrow(out), mean(tau[g == "HC"]), mean(tau[g == "DEP"]),
            mean(tau[g == "DEP"]) - mean(tau[g == "HC"]), t$p.value, median(out$nll_range), OUT))
