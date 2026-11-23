#!/usr/bin/env Rscript
# ==============================================================================
# PARAMETER RECOVERY -- baseline model (4a: mean-signal drift, threshold, nDT; tau = 0)
# ==============================================================================
# The manuscript's primary group comparison uses the baseline model, but Table S4
# reports recovery only for the time-contrast model (recovery_redesign_tau.R). This
# script recovers the baseline parameters with the SAME design as Table S4 so the
# two can be compared row by row:
#   * truth: d_mean ~ U(-0.50, 1.75), thres ~ U(1.10, 2.70), nDT ~ U(1.20, 3.00),
#     drawn with the same seed (20260804) and in the same order as recovery_redesign_tau.R
#     (d2 is drawn and discarded so thres/nDT get identical values)
#   * trials: signals resampled from the observed (mean, contrast) pool
#   * generator: shared forward simulator tddm_simulate_trial() with d2 = 0, tau = 0
#   * fit: shared engine likelihood mk_ll_static1(), canonical 4a bounds
#     lo = (-3, 0.5, 0.01), up = (3, 4, 3), DEoptim itermax 200 / NP 80, best of RESTARTS
# Usage:
#   Rscript scripts/tDDM/recovery_baseline.R --N 50 --reps 3 --restarts 3 --trials 20,80
#   Rscript scripts/tDDM/recovery_baseline.R --fast
# ==============================================================================
suppressPackageStartupMessages({ library(DEoptim); library(parallel) })
args <- commandArgs(trailingOnly = TRUE)
FAST <- any(args %in% "--fast")
get_arg <- function(flag, default) { i <- which(args == flag); if (length(i) == 1 && i < length(args)) args[i + 1] else default }
script_dir <- tryCatch(dirname(sys.frame(1)$ofile), error = function(e) {
  ca <- commandArgs(trailingOnly = FALSE); fa <- sub("^--file=", "", ca[grep("^--file=", ca)])
  if (length(fa) > 0) dirname(fa) else getwd() })
source(file.path(script_dir, "..", "..", "config.R"))
source(file.path(MODEL_DIR, "tDDM", "tddm_engine.R"))

TRIALS   <- as.integer(strsplit(get_arg("--trials", "20,80"), ",")[[1]])
N_SUBJ   <- as.integer(get_arg("--N", if (FAST) 6 else 50))
REPS     <- as.integer(get_arg("--reps", if (FAST) 1 else 3))
RESTARTS <- as.integer(get_arg("--restarts", if (FAST) 1 else 3))
ITERMAX  <- as.integer(get_arg("--itermax", if (FAST) 40 else 200))
NP       <- as.integer(get_arg("--np", if (FAST) 30 else 80))
N_CORES  <- as.integer(get_arg("--cores", max(1, detectCores() - 2)))
LO <- c(-3, 0.5, 0.01); UP <- c(3, 4, 3)

ts <- format(Sys.time(), "%Y%m%d_%H%M%S")
outdir <- file.path(RESULTS_DIR, "tDDM", sprintf("recovery_baseline_%s%s", ts, if (FAST) "_fast" else ""))
if (dir.exists(outdir)) stop("output directory exists: ", outdir)
dir.create(outdir, recursive = TRUE)
info <- file.path(outdir, "run_info.txt")
writeLines(c(paste("command: Rscript scripts/tDDM/recovery_baseline.R", paste(args, collapse = " ")),
             paste("commit:", tryCatch(system("git rev-parse --short HEAD", intern = TRUE), error = function(e) NA)),
             paste("dirty:", length(tryCatch(system("git status --porcelain scripts models", intern = TRUE), error = function(e) NULL)) > 0),
             paste("host:", Sys.info()[["nodename"]]), paste("R:", R.version.string),
             paste("DEoptim:", as.character(packageVersion("DEoptim"))),
             sprintf("settings: trials=%s N=%d reps=%d restarts=%d itermax=%d NP=%d cores=%d lo=(%s) up=(%s)",
                     paste(TRIALS, collapse = ","), N_SUBJ, REPS, RESTARTS, ITERMAX, NP, N_CORES,
                     paste(LO, collapse = ","), paste(UP, collapse = ",")),
             paste("start:", format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z"))), info)

# signal pool: identical construction to recovery_redesign_tau.R
files <- list.files(DATA_DIR, pattern = "^study4_.*\\.csv$", full.names = TRUE)
rd <- do.call(rbind, lapply(files, read.csv, stringsAsFactors = FALSE))
rd <- rd[is.finite(rd$choice_RT) & rd$choice_RT > 0.15, ]
pool <- data.frame(s1 = (rd$RPE_scale + rd$EPE_X1_scale) / 2, s2 = (rd$RPE_scale - rd$EPE_X1_scale) / 2)
pool <- pool[complete.cases(pool), ]

set.seed(20260804)                      # same draws as recovery_redesign_tau.R
d1 <- runif(N_SUBJ, -0.50, 1.75); invisible(runif(N_SUBJ, -2, 2))
truth <- data.frame(id = 1:N_SUBJ, d1 = d1, thres = runif(N_SUBJ, 1.10, 2.70), nDT = runif(N_SUBJ, 1.20, 3.00))

ll <- mk_ll_static1("signal1")
ctrl <- DEoptim.control(itermax = ITERMAX, NP = NP, trace = FALSE)   # never set VTR
one <- function(tp, nt, rep, seed) {
  set.seed(seed)
  s <- pool[sample(nrow(pool), nt, replace = TRUE), ]
  sim <- lapply(s$s1, function(x) tddm_simulate_trial(tp$d1, 0, tp$thres, tp$nDT, 0, x, 0))
  td <- data.frame(signal1 = s$s1, RT_ddm = vapply(sim, `[[`, 0, "rt"), choice_ddm = vapply(sim, `[[`, 0, "ch"))
  td <- td[td$RT_ddm > 0.15 & td$RT_ddm <= 10, ]
  if (nrow(td) < 5) return(NULL)
  fit <- NULL
  for (k in seq_len(RESTARTS)) {
    fk <- tryCatch(DEoptim(ll, LO, UP, control = ctrl, td = td), error = function(e) NULL)
    if (!is.null(fk) && (is.null(fit) || fk$optim$bestval < fit$optim$bestval)) fit <- fk
  }
  if (is.null(fit)) return(NULL)
  e <- fit$optim$bestmem
  data.frame(subject = tp$id, n_trials = nt, rep = rep, n_valid = nrow(td),
             true_d_mean = tp$d1, est_d_mean = e[1], true_thres = tp$thres, est_thres = e[2],
             true_nDT = tp$nDT, est_nDT = e[3], nll = fit$optim$bestval, de_iter = fit$optim$iter,
             accept_rate = mean(td$choice_ddm))
}

t0 <- Sys.time(); out <- list()
for (nt in TRIALS) {
  jobs <- expand.grid(si = seq_len(N_SUBJ), rp = seq_len(REPS))
  r <- mclapply(seq_len(nrow(jobs)), function(k) one(truth[jobs$si[k], ], nt, jobs$rp[k],
                5e6 + 1e4 * match(nt, TRIALS) + 100 * jobs$si[k] + jobs$rp[k]),
                mc.cores = N_CORES, mc.preschedule = FALSE)
  r <- Filter(is.data.frame, r); out[[length(out) + 1]] <- do.call(rbind, r)
  cat(sprintf("  trials=%d: %d/%d fits (%.1f min)\n", nt, length(r), nrow(jobs),
              as.numeric(difftime(Sys.time(), t0, units = "mins")))); flush.console()
}
fits <- do.call(rbind, out)
write.csv(fits, file.path(outdir, "recovery_fits.csv"), row.names = FALSE)

ci <- function(r, n) tanh(atanh(r) + c(-1.96, 1.96) / sqrt(n - 3))
summ <- do.call(rbind, lapply(TRIALS, function(nt) { s <- fits[fits$n_trials == nt, ]
  do.call(rbind, lapply(c("d_mean", "thres", "nDT"), function(p) {
    tr <- s[[paste0("true_", p)]]; es <- s[[paste0("est_", p)]]; r <- cor(tr, es); cc <- ci(r, nrow(s))
    data.frame(n_trials = nt, param = p, n_fits = nrow(s), r = r, r_lo = cc[1], r_hi = cc[2],
               bias = mean(es - tr), sd_err = sd(es - tr), sd_true = sd(tr),
               slope = unname(coef(lm(es ~ tr))[2]),
               pct_at_bound = 100 * mean(abs(es - LO[match(p, c("d_mean", "thres", "nDT"))]) < 0.05 |
                                         abs(es - UP[match(p, c("d_mean", "thres", "nDT"))]) < 0.05)) })) }))
write.csv(summ, file.path(outdir, "recovery_summary.csv"), row.names = FALSE)
print(summ, digits = 3)
cat(sprintf("DEoptim generations: median %.0f, min %.0f\n", median(fits$de_iter), min(fits$de_iter)))
cat(paste("end:", format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z"), "| exit 0\n"), file = info, append = TRUE)
cat("saved:", outdir, "\n")
