# Verifies that the C++ likelihood (models/tDDM/tddm_split.cpp, used by mk_ll_tau() when
# TDDM_CPP is TRUE) reproduces the R accurate path. 2026-09-28.
#   V1  stored canonical NLLs (mean-first and contrast-first, N = 351 each)
#   V2  3,000 random parameter points over the fitting box, plus near-zero tau
#   V3  simulated trial sets of the kind the controls fit
#   V4  full DEoptim fits (200 generations, NP 80) with the same seed: identical optimum
#   V5  the same evaluations inside forked mclapply workers
#   Rscript scripts/tDDM/verify_cpp_likelihood.R
suppressPackageStartupMessages({ library(DEoptim); library(parallel) })
if (!file.exists("config.R")) stop("run from the repository root")
source("config.R"); source(file.path(MODEL_DIR, "tDDM", "tddm_engine.R"))
stopifnot(isTRUE(TDDM_CPP))
R <- file.path(RESULTS_DIR, "tDDM"); data <- tddm_load_study4(DATA_DIR)
mk_both <- function(s1, s2) {
  cpp <- mk_ll_tau(s1, s2, trans = p_trans_image)
  TDDM_CPP <<- FALSE; r <- mk_ll_tau(s1, s2, trans = p_trans_image); TDDM_CPP <<- TRUE
  list(cpp = cpp, r = r)
}
MF <- mk_both("mean_sig", "contrast_sig"); CF <- mk_both("contrast_sig", "mean_sig")
stopifnot(!identical(body(MF$cpp), body(MF$r)))          # the two paths really differ
cols <- c("d_mean", "d_contrast", "thres", "nDT", "tau")
tdof <- function(s) data[data$subject == s, ]
cat("== V1: stored canonical NLLs ==\n")
for (lab in c("mean-first", "contrast-first")) {
  f <- read.csv(file.path(R, "canonical_4j_accurate_bestof10",
                          if (lab == "mean-first") "fits_4j_accurate_best10.csv" else "fits_4j_accurate_contrast_first_best10.csv"))
  LL <- if (lab == "mean-first") MF$cpp else CF$cpp
  tcol <- if ("tau_raw" %in% names(f)) "tau_raw" else "tau"  # the stored NLL belongs to the unclamped estimate
  d <- sapply(seq_len(nrow(f)), function(i) {
    p <- unlist(f[i, cols]); p[5] <- f[[tcol]][i]
    if (lab == "contrast-first") p[1:2] <- p[2:1]      # position 1 = leading signal's weight
    LL(p, tdof(f$subject[i])) - f$negLogLik[i] })
  cat(sprintf("  %-15s N = %d | max |C++ - stored| = %.1e\n", lab, nrow(f), max(abs(d))))
}
cat("== V2: random and near-zero-tau points, C++ vs R ==\n")
set.seed(20260928)
LO <- c(-3, -3, 0.5, 0.01, -2); UP <- c(3, 3, 4, 3, 3)
P <- t(replicate(3000, runif(5, LO, UP)))
P[1:300, 5] <- 10^runif(300, -6, -1)                         # near-zero onsets
subs <- sample(unique(data$subject), 3000, replace = TRUE)
d2 <- sapply(seq_len(nrow(P)), function(i) { td <- tdof(subs[i]); MF$cpp(P[i, ], td) - MF$r(P[i, ], td) })
cat(sprintf("  3000 points | max |diff| %.1e | median %.1e | > 1e-9: %d\n", max(abs(d2)), median(abs(d2)), sum(abs(d2) > 1e-9)))
cat("== V3: simulated trial sets (controls) ==\n")
d3 <- sapply(1:300, function(i) {
  p <- runif(5, c(-1, -2, 1.1, 1.2, 0.1), c(1.8, 2, 2.7, 3, 2.7))
  s <- data[sample(nrow(data), 20), ]
  sim <- tddm_simulate_block(p[1], p[2], p[3], p[4], p[5], s$mean_sig, s$contrast_sig)
  td <- data.frame(mean_sig = s$mean_sig, contrast_sig = s$contrast_sig, RT_ddm = sim$RT, choice_ddm = sim$acc)
  td <- td[!is.na(td$choice_ddm) & td$RT_ddm > 0.15 & td$RT_ddm <= 10, ]
  q <- runif(5, LO, UP); MF$cpp(q, td) - MF$r(q, td) })
cat(sprintf("  300 simulated sets | max |diff| %.1e\n", max(abs(d3))))
cat("== V4: full DEoptim fits, same seed ==\n")
ctl <- DEoptim.control(itermax = 200, NP = 80, trace = FALSE)
for (s in sample(unique(data$subject), 4)) {
  td <- tdof(s)
  set.seed(s); t1 <- system.time(a <- DEoptim(MF$r, LO, UP, ctl, td = td))[3]
  set.seed(s); t2 <- system.time(b <- DEoptim(MF$cpp, LO, UP, ctl, td = td))[3]
  pa <- a$optim$bestmem; pb <- b$optim$bestmem; pa[5] <- max(pa[5], 0); pb[5] <- max(pb[5], 0)   # tau <= 0 is one model
  cat(sprintf("  subject %3d | NLL R %.10f  C++ %.10f | max |par diff| (tau<=0 -> 0) %.1e | %.1f s vs %.1f s\n",
              s, a$optim$bestval, b$optim$bestval, max(abs(pa - pb)), t1, t2))
}
cat("== V5: forked workers ==\n")
seq_v <- sapply(1:40, function(i) MF$cpp(P[i, ], tdof(subs[i])))
par_v <- unlist(mclapply(1:40, function(i) MF$cpp(P[i, ], tdof(subs[i])), mc.cores = 4))
cat(sprintf("  40 evaluations in 4 forked workers | max |parallel - sequential| %.1e\n", max(abs(par_v - seq_v))))
