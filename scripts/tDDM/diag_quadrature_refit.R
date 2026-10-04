# Diagnostic (2026-09-27): does the 30-point split-at-tau quadrature drive the 4j fits?
#
# ll_piecewise_trial() integrates p_reach(x) * p_from(x | rem) over a fixed 30-point
# grid on [0.01, bsep - 0.01]. When rem = dt_obs - tau is tiny, p_from is a boundary
# layer of width ~sqrt(rem) and the grid overweights it (up to ~7x per trial). In the
# canonical best-of-10 fits, 97.7% of participants with tau > 0 place tau within 1 ms
# of one trial's decision time (chance 0.3%).
#
# This script refits the mean-first 4j ONCE with an accurate integral (adaptive
# quadrature with the two boundary layers split out) and the same settings and seed
# as optimizer_reproducibility_r1_20260820_011554/fits_rep01.csv, then compares.
# The shared engine is not modified.
#
#   Rscript scripts/tDDM/diag_quadrature_refit.R [--cores 12] [--seed 20260821] [--nsub 0]
suppressPackageStartupMessages({ library(DEoptim); library(parallel) })
args <- commandArgs(trailingOnly = TRUE)
ga <- function(f, d) { if (f %in% args) { i <- which(args == f); if (i < length(args)) return(args[i+1]) }; d }
NC   <- as.integer(ga("--cores", 12))
SEED <- as.integer(ga("--seed", 20260821))       # = rep01 of the canonical 10-run set
NSUB <- as.integer(ga("--nsub", 0))              # 0 = all participants
if (!file.exists("config.R")) stop("run from the repository root")
source("config.R"); source(file.path(MODEL_DIR, "tDDM", "tddm_engine.R"))

# accurate split-at-tau trial likelihood (same model, different quadrature)
ll_piecewise_trial_acc <- function(drift1, drift2, tau, bsep, w, nDT, obs_rt, obs_choice) {
  dt_obs <- obs_rt - nDT
  if (dt_obs <= 0.001) return(1e-12)
  if (dt_obs < tau || tau <= 0) {
    eff <- if (tau <= 0) drift2 else drift1
    d <- if (obs_choice == 1) dwiener_nf(dt_obs, -eff, bsep, 1 - w) else dwiener_nf(dt_obs, eff, bsep, w)
    return(max(d, 1e-12))
  }
  rem <- dt_obs - tau
  fx <- function(x) p_trans_image(x, tau, bsep * w, drift1, bsep) *
    (if (obs_choice == 1) dwiener_nf_vec(rem, -drift2, bsep, 1 - x / bsep)
     else                 dwiener_nf_vec(rem,  drift2, bsep, x / bsep))
  h  <- min(bsep / 4, 20 * sqrt(rem))              # boundary-layer width
  br <- unique(c(0, h, bsep - h, bsep))
  v <- tryCatch(sum(vapply(seq_len(length(br) - 1), function(k)
         integrate(fx, br[k], br[k + 1], subdivisions = 1000L, rel.tol = 1e-6)$value, 0)),
       error = function(e) {                        # fallback: dense grid incl. endpoints
         xg <- seq(0, bsep, length.out = 4001); y <- fx(xg)
         sum((y[-1] + y[-length(y)]) / 2) * diff(xg)[1] })
  max(v, 1e-12)
}
mk_ll_tau_acc <- function(s1, s2) function(params, td) {
  d1 <- params[1]; d2 <- params[2]; thres <- params[3]; nDT <- params[4]; tau <- params[5]
  bsep <- 2 * thres / TDDM_SIGMA
  if (bsep <= 0.01 || nDT <= 0) return(TDDM_PENALTY)
  nll <- 0
  for (i in seq_len(nrow(td))) {
    drift1 <- d1 * td[[s1]][i] / TDDM_SIGMA
    drift2 <- (d1 * td[[s1]][i] + d2 * td[[s2]][i]) / TDDM_SIGMA
    nll <- nll - log(max(ll_piecewise_trial_acc(drift1, drift2, tau, bsep, 0.5, nDT,
                                                td$RT_ddm[i], td$choice_ddm[i]), TDDM_EPS))
  }
  if (!is.finite(nll)) TDDM_PENALTY else nll
}
LLacc <- mk_ll_tau_acc("mean_sig", "contrast_sig")
LLgrid <- mk_ll_tau("mean_sig", "contrast_sig", trans = get_transition("image"))
LO <- c(-3, -3, 0.5, 0.01, -2); UP <- c(3, 3, 4, 3, 3)
ctrl <- DEoptim.control(itermax = 200, NP = 80, trace = FALSE)   # NEVER set VTR

data <- tddm_load_study4(DATA_DIR)
subjects <- unique(data$subject)
if (NSUB > 0) { set.seed(1); subjects <- sort(sample(subjects, NSUB)) }
cat(sprintf("=== quadrature diagnostic: 4j mean-first, accurate integral | %d subjects | seed %d | cores %d ===\n",
            length(subjects), SEED, NC))

fit_one <- function(s) {
  td <- data[data$subject == s, ]
  f <- tryCatch(DEoptim(LLacc, LO, UP, ctrl, td = td), error = function(e) NULL)
  if (is.null(f)) return(NULL)
  e <- f$optim$bestmem
  dt <- td$RT_ddm - e[4]; g <- dt[dt >= e[5]] - e[5]
  data.frame(subject = s, group = td$group[1], n_trials = nrow(td),
             negLogLik = f$optim$bestval, nll_grid_at_est = LLgrid(e, td), de_iter = f$optim$iter,
             d_mean = e[1], d_contrast = e[2], thres = e[3], nDT = e[4], tau = e[5],
             min_gap = if (length(g)) min(g) else NA, stringsAsFactors = FALSE)
}
RNGkind("L'Ecuyer-CMRG"); set.seed(SEED)
t0 <- Sys.time()
res <- do.call(rbind, Filter(is.data.frame, mclapply(subjects, fit_one, mc.cores = NC)))
stopifnot(min(res$de_iter) > 0)
res$BIC <- 2 * res$negLogLik + 5 * log(res$n_trials)
od <- file.path(RESULTS_DIR, "tDDM", sprintf("diag_quadrature_refit_%s", format(Sys.time(), "%Y%m%d_%H%M%S")))
dir.create(od, recursive = TRUE, showWarnings = FALSE)
write.csv(res, file.path(od, "fits_4j_accurate_rep01seed.csv"), row.names = FALSE)

# comparison with the grid fits (same seed) and the tau-free 4c (unaffected: no split)
# (the grid run is absent from the public repository; the comparison is then skipped)
gf <- file.path(RESULTS_DIR, "tDDM", "optimizer_reproducibility_r1_20260820_011554", "fits_rep01.csv")
if (!file.exists(gf)) { cat("Saved:", od, "(grid comparison skipped)\n"); quit(save = "no") }
g1 <- read.csv(gf)
c4 <- read.csv(file.path(RESULTS_DIR, "tDDM", "corrected_likelihood_20260819_005659", "fits_4c_static.csv"))
m <- merge(merge(res, g1[, c("subject", "tau", "negLogLik")], by = "subject", suffixes = c("", "_grid")),
           c4[, c("subject", "BIC")], by = "subject", suffixes = c("", "_4c"))
grp <- function(x, g) { t <- t.test(x[g == "DEP"], x[g == "HC"])
  sprintf("HC %.3f DEP %.3f | DEP-HC %+.3f t(%.0f)=%.2f p=%.4f", mean(x[g == "HC"]), mean(x[g == "DEP"]),
          mean(x[g == "DEP"]) - mean(x[g == "HC"]), t$parameter, t$statistic, t$p.value) }
sink(file.path(od, "summary.txt"), split = TRUE)
cat(sprintf("N = %d | %.1f min\n", nrow(m), as.numeric(difftime(Sys.time(), t0, units = "mins"))))
cat("tau (accurate) :", grp(m$tau, m$group), "\n")
cat("tau (grid rep01):", grp(m$tau_grid, m$group), "\n")
for (p in c("d_mean", "thres", "nDT", "d_contrast")) cat(sprintf("%-15s: %s\n", p, grp(m[[p]], m$group)))
pos <- m$tau > 0
cat(sprintf("tau<=0: %.1f%% | tau>0 within 1 ms of a trial DT: %.1f%% (grid canonical 97.7%%, chance ~0.3%%)\n",
            100 * mean(!pos), 100 * mean(m$min_gap[pos] < 1e-3, na.rm = TRUE)))
cat(sprintf("r(tau accurate, tau grid) = %.3f\n", cor(m$tau, m$tau_grid)))
db <- m$BIC - m$BIC_4c
cat(sprintf("dBIC(4j accurate - 4c): mean %+.3f | <0 %.1f%% | <-2 %.1f%%\n", mean(db), 100 * mean(db < 0), 100 * mean(db < -2)))
cat(sprintf("grid NLL at the accurate estimates minus accurate NLL: median %+.3f\n",
            median(m$nll_grid_at_est - m$negLogLik)))
sink()
cat("Saved:", od, "\n")
