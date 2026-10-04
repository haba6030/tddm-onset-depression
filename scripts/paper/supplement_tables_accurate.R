#!/usr/bin/env Rscript
# Numbers for Supplementary Tables S3-S10 under the accurate-integral canonical fits
# (tau <= 0 set to 0) and the stage-6 controls. 2026-09-30.
#   Rscript scripts/paper/supplement_tables_accurate.R
suppressPackageStartupMessages(library(BayesFactor))
if (!file.exists("config.R")) stop("run from the repository root")
source("config.R"); source(file.path(MODEL_DIR, "tDDM", "tddm_engine.R"))
R <- file.path(RESULTS_DIR, "tDDM"); OD <- file.path(R, "canonical_4j_accurate_bestof10")
sink(file.path(OD, "supplement_tables.txt"), split = TRUE)
welch <- function(x, g) {                     # DEP - HC, Cohen's d with pooled SD and its CI
  a <- x[g == "DEP"]; b <- x[g == "HC"]; t <- t.test(a, b); n1 <- length(a); n2 <- length(b)
  sp <- sqrt(((n1 - 1) * var(a) + (n2 - 1) * var(b)) / (n1 + n2 - 2)); d <- (mean(a) - mean(b)) / sp
  se <- sqrt((n1 + n2) / (n1 * n2) + d^2 / (2 * (n1 + n2)))
  c(diff = mean(a) - mean(b), t = unname(abs(t$statistic)), df = unname(t$parameter), p = t$p.value,
    d = d, dlo = d - 1.96 * se, dhi = d + 1.96 * se)
}
f <- read.csv(file.path(OD, "fits_4j_accurate_best10.csv")); f$tau <- pmax(f$tau, 0)
dat <- tddm_load_study4(DATA_DIR)
ps <- do.call(rbind, lapply(split(dat, dat$subject), function(s) data.frame(
  subject = s$subject[1], ces = s$cesScore[1], age = s$age[1], gender = s$gender[1])))
x <- merge(f, ps, by = "subject"); P <- c("d_mean", "d_contrast", "thres", "nDT", "tau")

cat("== Table S3: non-decision-time bound ladder (one fit per bound, 3 restarts; tau <= 0 -> 0) ==\n")
LD <- file.path(R, "nDT_bound_ladder_engine_20260928_192521")
for (ub in c(1.5, 2.0, 2.5, 3.0, 3.5)) {
  l <- read.csv(file.path(LD, sprintf("fits_4j_nDT%.1f.csv", ub))); w <- welch(l$tau, l$group)
  ce <- function(k) 100 * mean(l$nDT[k] >= ub - 0.05)
  cat(sprintf("  %.1f | median nDT %.3f | ceiling %.1f (HC %.1f, DEP %.1f) | dtau %+.3f | t(%.0f) = %.2f | p = %.4f | d = %.2f [%.2f, %.2f]\n",
              ub, median(l$nDT), ce(TRUE), ce(l$group == "HC"), ce(l$group == "DEP"), w["diff"], w["df"], w["t"], w["p"], w["d"], w["dlo"], w["dhi"]))
}

cat("\n== Table S4: individual recovery ==\n")
rc <- read.csv(file.path(R, "recovery_redesign_tau_20260928_192519", "recovery_fits.csv"))
rs <- read.csv(file.path(R, "recovery_redesign_tau_20260928_192519", "recovery_summary.csv"))
for (i in seq_len(nrow(rs))) with(rs[i, ], cat(sprintf(
  "  %d trials | r: drift %.3f, contrast %.3f, threshold %.3f, nDT %.3f, onset %+.3f [%.3f, %.3f] | error SD %.3f | bias %+.3f | slope %+.3f | MC SE of tau (median over participants) %.3f\n",
  n_trials, r_d1, r_d2, r_thres, r_nDT, r_tau, r_lo, r_hi, sd_err_tau, bias, slope, mc_se_tau)))
for (nt in c(20, 80)) {
  s <- rc[rc$n_trials == nt, ]; rr <- sapply(split(s, s$rep), function(g) cor(g$true_tau, g$est_tau))
  b <- cut(s$true_tau, c(0.3, 0.9, 1.5, 2.1, 2.7), include.lowest = TRUE)
  cat(sprintf("  %d trials | SD of onset r across the %d simulations %.3f | bias by true onset: %s\n", nt, length(rr), sd(rr),
              paste(sprintf("%s %+.2f", levels(b), tapply(s$est_tau - s$true_tau, b, mean)), collapse = ", ")))
}
nl <- read.csv(file.path(R, "nulltau_corrected_image_20260928_211818", "per_rep.csv"))
pc <- read.csv(file.path(R, "positive_control_tau_20260928_231637", "per_point.csv")); p0 <- pc[pc$true_gap == 0, ]
cat(sprintf("  bias at the fitted onset, same true onset in both groups: no-onset simulation (tau 0.891) HC %+.3f, DEP %+.3f | positive control at zero gap (tau 0.804) HC %+.3f, DEP %+.3f\n",
            mean(nl$tauhat_HC) - 0.891, mean(nl$tauhat_DEP) - 0.891, mean(p0$tauhat_HC) - 0.804, mean(p0$tauhat_DEP) - 0.804))

cat("\n== Table S6: negative control, DEP - HC (per_rep stores HC - DEP) ==\n")
for (p in c("d_mean", "d_contrast", "thres", "nDT"))
  cat(sprintf("  %-10s generating %+.3f | recovered %+.3f | difference %+.3f\n", p, -mean(nl[[paste0(p, "_true")]]),
              -mean(nl[[paste0(p, "_rec")]]), -mean(nl[[paste0(p, "_rec")]] - nl[[paste0(p, "_true")]])))
cat(sprintf("  onset      generating 0 | recovered %+.3f (SD %.3f, %d reps)\n", -mean(nl$recovered_dtau), sd(nl$recovered_dtau), nrow(nl)))
for (gp in c("222147", "224927")) { d <- read.csv(file.path(R, paste0("nulltau_corrected_image_20260928_", gp), "per_rep.csv"))
  cat(sprintf("  starting-point gap %.2f: recovered %+.3f (SD %.3f, %d reps)\n", d$wgap[1], -mean(d$recovered_dtau), sd(d$recovered_dtau), nrow(d))) }

cat("\n== Table S7: positive control ==\n")
for (k in list(list("20 trials", pc), list("80 trials", read.csv(file.path(R, "positive_control_tau_20260929_012719", "per_point.csv"))))) {
  d <- k[[2]]; a <- aggregate(recovered_gap ~ true_gap, d, function(v) c(m = mean(v), s = sd(v), n = length(v)))
  m <- lm(recovered_gap ~ true_gap, d); ci <- confint(m)
  cat(sprintf("  %s: %s | slope %.3f [%.3f, %.3f], intercept %.3f (SE %.3f)\n", k[[1]],
              paste(sprintf("%.3f -> %.3f (SD %.3f, n %d)", a$true_gap, a$recovered_gap[, "m"], a$recovered_gap[, "s"], a$recovered_gap[, "n"]), collapse = "; "),
              coef(m)[2], ci[2, 1], ci[2, 2], coef(m)[1], sqrt(vcov(m)[1, 1])))
}

cat("\n== Table S8: response-time signature ==\n")
print(read.csv(file.path(R, "rt_signature_check_20260928_005349", "comparison.csv")), row.names = FALSE)

cat("\n== contrast-weight subgroups (interaction) ==\n")
bd <- abs(abs(f$d_contrast) - 3) < 0.05
d <- data.frame(tau = f$tau, DEP = as.numeric(f$group == "DEP"), bnd = as.numeric(bd))
m <- lm(tau ~ DEP * bnd, d); k <- summary(m)$coefficients["DEP:bnd", ]; ci <- confint(m)["DEP:bnd", ]
bf <- extractBF(lmBF(tau ~ DEP + bnd + DEP:bnd, data = d) / lmBF(tau ~ DEP + bnd, data = d))$bf
cat(sprintf("  inside %d, at bound %d | interaction (at bound minus inside) %+.3f SE %.3f t(%d) = %.2f p = %.3f 95%% CI [%.2f, %.2f] BF01 = %.2f\n",
            sum(!bd), sum(bd), k[1], k[2], m$df.residual, k[3], k[4], ci[1], ci[2], 1 / bf))

cat("\n== Table S9: correlations with CES-D and age ==\n")
for (p in P) { a <- cor.test(x[[p]], x$ces); b <- cor.test(x[[p]], x$age)
  cat(sprintf("  %-10s CES-D r %+.3f p %.4f [%+.3f, %+.3f] | age r %+.3f p %.4f\n", p, a$estimate, a$p.value, a$conf.int[1], a$conf.int[2], b$estimate, b$p.value)) }

cat("\n== Table S10: group contrast before and after age and gender (standardized parameter) ==\n")
x$DEP <- as.numeric(x$group == "DEP")
for (p in P) { w <- welch(x[[p]], x$group); z <- (x[[p]] - mean(x[[p]])) / sd(x[[p]])
  m <- lm(z ~ DEP + age + factor(gender), x); k <- summary(m)$coefficients["DEP", ]; ci <- confint(m)["DEP", ]
  cat(sprintf("  %-10s d %+.3f t(%.0f) = %.2f p %.4f [%+.2f, %+.2f] | beta %+.3f t(%d) = %.2f p %.4f [%+.3f, %+.3f]\n",
              p, w["d"], w["df"], w["t"], w["p"], w["dlo"], w["dhi"], k[1], m$df.residual, abs(k[3]), k[4], ci[1], ci[2])) }
m <- lm(tau ~ DEP + age + factor(gender), x); cat(sprintf("  onset coefficient on its original scale: %.3f s\n", coef(m)["DEP"]))
sink()
