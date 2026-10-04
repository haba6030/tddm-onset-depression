#!/usr/bin/env Rscript
# Every Results number that comes from the 4j fits file, recomputed on the
# accurate-integral best-of-10 fits (canonical_4j_accurate_bestof10). tau <= 0 is
# set to 0: every tau <= 0 gives the same single-phase model (identical NLL to 4c),
# so the value below 0 carries no information. Same statistics as reanalysis_bestrun.R.
#   Rscript scripts/tDDM/memo_0919/reanalysis_accurate.R
suppressPackageStartupMessages({ library(BayesFactor) })
if (!file.exists("config.R")) stop("run from the repository root")
source("config.R"); source(file.path(MODEL_DIR, "tDDM", "tddm_engine.R"))
R <- file.path(RESULTS_DIR, "tDDM"); OD <- file.path(R, "canonical_4j_accurate_bestof10")
f <- read.csv(file.path(OD, "fits_4j_accurate_best10.csv"))
f$tau_raw <- f$tau; f$tau <- pmax(f$tau, 0)
g <- f$group
sink(file.path(OD, "reanalysis_accurate.txt"), split = TRUE)
wt <- function(x, g, lab) {
  a <- x[g == "DEP"]; b <- x[g == "HC"]; t <- t.test(a, b)
  n1 <- length(a); n2 <- length(b)
  sp <- sqrt(((n1 - 1) * var(a) + (n2 - 1) * var(b)) / (n1 + n2 - 2)); d <- (mean(a) - mean(b)) / sp
  se <- sqrt((n1 + n2) / (n1 * n2) + d^2 / (2 * (n1 + n2)))
  cat(sprintf("  %-26s HC=%.3f DEP=%.3f diff=%+.3f SE=%.3f CI[%+.3f, %+.3f] t(%.0f)=%.2f p=%.4f d=%.2f [%.2f, %.2f]\n",
      lab, mean(b), mean(a), mean(a) - mean(b), unname(t$stderr), t$conf.int[1], t$conf.int[2],
      t$parameter, abs(t$statistic), t$p.value, d, d - 1.96 * se, d + 1.96 * se))
  invisible(t)
}
cat(sprintf("N = %d (HC %d, DEP %d)\n", nrow(f), sum(g == "HC"), sum(g == "DEP")))
cat("=== 1. group comparisons ===\n")
for (p in c("d_mean", "d_contrast", "thres", "nDT", "tau")) wt(f[[p]], g, p)

cat("\n=== 2. boundary rates ===\n")
for (tol in c(0.05, 0.01)) {
  dc <- abs(abs(f$d_contrast) - 3) < tol; nd <- abs(f$nDT - 3) < tol
  cat(sprintf("  tol %.2f: d_contrast at bound %.1f%% | nDT at ceiling %.1f%% (HC %.1f%% / DEP %.1f%%)\n",
      tol, 100 * mean(dc), 100 * mean(nd), 100 * mean(nd[g == "HC"]), 100 * mean(nd[g == "DEP"])))
}
cat(sprintf("  tau <= 0 (set to 0): %.1f%% (n = %d; HC %.1f%%, DEP %.1f%%)\n", 100 * mean(f$tau_raw <= 0),
    sum(f$tau_raw <= 0), 100 * mean(f$tau_raw[g == "HC"] <= 0), 100 * mean(f$tau_raw[g == "DEP"] <= 0)))

cat("\n=== 3. tau: robustness ===\n")
keep <- f$tau_raw > 0
cat(sprintf("  excluding tau<=0 (n=%d dropped, n=%d kept):\n", sum(!keep), sum(keep))); wt(f$tau[keep], g[keep], "tau")
w <- wilcox.test(f$tau[g == "DEP"], f$tau[g == "HC"], conf.int = TRUE)
nh <- sum(g == "HC"); nd <- sum(g == "DEP"); z <- qnorm(w$p.value / 2)
cat(sprintf("  Wilcoxon W=%.0f p=%.4f HL shift=%+.3f CI[%+.2f, %+.2f] r=%.2f\n", w$statistic, w$p.value,
    w$estimate, w$conf.int[1], w$conf.int[2], abs(z) / sqrt(nh + nd)))
bd <- abs(abs(f$d_contrast) - 3) < 0.05
cat(sprintf("  contrast weight inside bounds (n=%d):\n", sum(!bd))); wt(f$tau[!bd], g[!bd], "tau")
cat(sprintf("  contrast weight at bound (n=%d):\n", sum(bd)));     wt(f$tau[bd], g[bd], "tau")
d <- data.frame(tau = f$tau, DEP = as.numeric(g == "DEP"), bnd = as.numeric(bd))
m <- summary(lm(tau ~ DEP * bnd, d)); k <- m$coefficients["DEP:bnd", ]
bf <- extractBF(lmBF(tau ~ DEP + bnd + DEP:bnd, data = d) / lmBF(tau ~ DEP + bnd, data = d))$bf
cat(sprintf("  interaction b=%+.3f SE=%.3f t(%d)=%.2f p=%.3f | BF01=%.2f\n", k[1], k[2], m$df[2], k[3], k[4], 1 / bf))

cat("\n=== 4. CES-D, age, gender, correlations ===\n")
dat <- tddm_load_study4(DATA_DIR)
ps <- do.call(rbind, lapply(split(dat, dat$subject), function(s) data.frame(
  subject = s$subject[1], ces = s$cesScore[1], age = s$age[1], gender = s$gender[1], medRT = median(s$RT_ddm))))
x <- merge(f, ps, by = "subject"); x$DEP <- as.numeric(x$group == "DEP")
ctp <- function(a, b, lab) { ct <- cor.test(a, b)
  cat(sprintf("  %-14s r(%d)=%.3f p=%.4f CI[%.2f, %.2f]\n", lab, ct$parameter, ct$estimate, ct$p.value, ct$conf.int[1], ct$conf.int[2])) }
ctp(x$ces, x$tau, "CES-D ~ tau"); ctp(x$age, x$tau, "age ~ tau")
ctp(x$tau, x$thres, "tau ~ thres"); ctp(x$tau, x$medRT, "tau ~ medRT"); ctp(x$thres, x$medRT, "thres ~ medRT")
m <- lm(tau ~ DEP + age + factor(gender), x); k <- coef(summary(m))["DEP", ]; ci <- confint(m)["DEP", ]
cat(sprintf("  tau ~ group + age + gender: b=%+.3f SE=%.3f t(%d)=%.2f p=%.4f CI[%.3f, %.3f]\n",
    k[1], k[2], m$df.residual, k[3], k[4], ci[1], ci[2]))

cat("\n=== 5. ANCOVA ===\n")
for (fm in list(tau ~ DEP, tau ~ DEP + thres, tau ~ DEP + thres + nDT, tau ~ DEP + thres + nDT + d_mean,
                tau ~ DEP + thres + nDT + d_mean + medRT)) {
  m <- lm(fm, x); k <- coef(summary(m))["DEP", ]; ci <- confint(m)["DEP", ]
  r0 <- summary(update(m, . ~ . - DEP))$r.squared
  cat(sprintf("  %-38s b=%+.3f SE=%.3f t(%d)=%.2f p=%.4f CI[%.3f, %.3f] unique R2=%.1f%%\n", deparse(fm[[3]]),
      k[1], k[2], m$df.residual, k[3], k[4], ci[1], ci[2], 100 * (summary(m)$r.squared - r0)))
}

cat("\n=== 6. BIC ladder ===\n")
c4 <- read.csv(file.path(R, "corrected_likelihood_20260819_005659", "fits_4c_static.csv"))
a4 <- read.csv(file.path(R, "corrected_likelihood_20260819_012638", "fits_4a_static1.csv"))
c4 <- c4[match(f$subject, c4$subject), ]; a4 <- a4[match(f$subject, a4$subject), ]
stopifnot(all(f$subject == c4$subject), all(f$subject == a4$subject))
for (cmp in list(list("4j - 4a", f$BIC - a4$BIC), list("4c - 4a", c4$BIC - a4$BIC), list("4j - 4c", f$BIC - c4$BIC)))
  cat(sprintf("  dBIC(%s): mean %+.3f | median %+.3f | <0 %.1f%% | <-2 %.1f%% | >+2 %.1f%%\n", cmp[[1]], mean(cmp[[2]]),
      median(cmp[[2]]), 100 * mean(cmp[[2]] < 0), 100 * mean(cmp[[2]] < -2), 100 * mean(cmp[[2]] > 2)))
sink()
