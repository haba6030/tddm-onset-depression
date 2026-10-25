#!/usr/bin/env Rscript
# Numbers added by the 2026-10-19 reframe (baseline model primary, onset exploratory):
#   1. baseline-model group comparisons (Results, Fig. 3a-b)
#   2. baseline rows of Table S9 (CES-D and age) and Table S10 (age and gender adjustment)
#   3. lowest-BIC model by group and its chi-square test (Results "Model comparison")
#   4. BIC of the two signal orderings of the time-contrast model (Results, exploratory section)
#   Rscript scripts/paper/baseline_primary_accurate.R
# Output: results/tDDM/canonical_4j_accurate_bestof10/baseline_primary.txt
if (!file.exists("config.R")) stop("run from the repository root")
source("config.R"); source(file.path(MODEL_DIR, "tDDM", "tddm_engine.R"))
R <- file.path(RESULTS_DIR, "tDDM"); OD <- file.path(R, "canonical_4j_accurate_bestof10")
sink(file.path(OD, "baseline_primary.txt"), split = TRUE)
welch <- function(x, g) {                     # DEP - HC, Cohen's d with pooled SD and its CI
  a <- x[g == "DEP"]; b <- x[g == "HC"]; t <- t.test(a, b); n1 <- length(a); n2 <- length(b)
  sp <- sqrt(((n1 - 1) * var(a) + (n2 - 1) * var(b)) / (n1 + n2 - 2)); d <- (mean(a) - mean(b)) / sp
  se <- sqrt((n1 + n2) / (n1 * n2) + d^2 / (2 * (n1 + n2)))
  c(hc = mean(b), dep = mean(a), t = unname(abs(t$statistic)), df = unname(t$parameter), p = t$p.value,
    d = d, dlo = d - 1.96 * se, dhi = d + 1.96 * se)
}
f4a <- read.csv(file.path(R, "corrected_likelihood_20260819_012638", "fits_4a_static1.csv"))
f4c <- read.csv(file.path(R, "corrected_likelihood_20260819_005659", "fits_4c_static.csv"))
f4j <- read.csv(file.path(OD, "fits_4j_accurate_best10.csv"))
f4k <- read.csv(file.path(OD, "fits_4j_accurate_contrast_first_best10.csv"))
srt <- function(d) d[order(d$subject), ]      # the fits files store participants in different orders
f4a <- srt(f4a); f4c <- srt(f4c); f4j <- srt(f4j); f4k <- srt(f4k)
stopifnot(nrow(f4a) == 351, all(f4a$subject == f4j$subject), all(f4c$subject == f4j$subject),
          all(f4k$subject == f4j$subject), all(f4a$group == f4j$group))
dat <- tddm_load_study4(DATA_DIR)
ps <- do.call(rbind, lapply(split(dat, dat$subject), function(s) data.frame(
  subject = s$subject[1], ces = s$cesScore[1], age = s$age[1], gender = s$gender[1])))
x <- merge(f4a, ps, by = "subject"); x$DEP <- as.numeric(x$group == "DEP")
P <- c(d_mean = "mean-signal drift", thres = "threshold")
cat(sprintf("N = %d (HC %d, DEP %d)\n", nrow(x), sum(x$DEP == 0), sum(x$DEP == 1)))

cat("\n== 1. baseline model (4a): group comparisons ==\n")
for (p in names(P)) { w <- welch(x[[p]], x$group)
  cat(sprintf("  %-18s HC = %.3f DEP = %.3f t(%.0f) = %.2f p = %.4f d = %+.2f [%+.2f, %+.2f]\n",
              P[p], w["hc"], w["dep"], w["df"], w["t"], w["p"], w["d"], w["dlo"], w["dhi"])) }

cat("\n== 2a. Table S9, baseline rows: correlations with CES-D and age ==\n")
for (p in names(P)) { a <- cor.test(x[[p]], x$ces); b <- cor.test(x[[p]], x$age)
  cat(sprintf("  %-18s CES-D r(%d) = %+.3f p %.4f [%+.3f, %+.3f] | age r %+.3f p %.4f\n", P[p], a$parameter,
              a$estimate, a$p.value, a$conf.int[1], a$conf.int[2], b$estimate, b$p.value)) }

cat("\n== 2b. Table S10, baseline rows: group contrast before and after age and gender (standardized) ==\n")
for (p in names(P)) { w <- welch(x[[p]], x$group); z <- (x[[p]] - mean(x[[p]])) / sd(x[[p]])
  m <- lm(z ~ DEP + age + factor(gender), x); k <- summary(m)$coefficients["DEP", ]; ci <- confint(m)["DEP", ]
  cat(sprintf("  %-18s d %+.3f t(%.0f) = %.2f p %.4f [%+.2f, %+.2f] | beta %+.3f t(%d) = %.2f p %.4f [%+.3f, %+.3f]\n",
              P[p], w["d"], w["df"], w["t"], w["p"], w["dlo"], w["dhi"], k[1], m$df.residual, abs(k[3]), k[4], ci[1], ci[2])) }

cat("\n== 3. lowest-BIC model by group (baseline 4a, contrast 4c, time-contrast 4j) ==\n")
B <- cbind(baseline = f4a$BIC, contrast = f4c$BIC, time_contrast = f4j$BIC)
best <- factor(colnames(B)[apply(B, 1, which.min)], levels = colnames(B))
tb <- table(group = factor(f4j$group, levels = c("HC", "DEP")), best = best)
print(tb); cat("  row percentages\n"); print(round(100 * prop.table(tb, 1), 1))
ch <- suppressWarnings(chisq.test(tb))
cat(sprintf("  chi-square(%d) = %.2f p = %.4f | Fisher exact p = %.4f\n",
            ch$parameter, ch$statistic, ch$p.value, fisher.test(tb)$p.value))
d <- f4j$BIC - f4a$BIC
cat(sprintf("  dBIC(4j - 4a): mean %+.3f | time-contrast ahead (dBIC < 0) for %.1f%%\n", mean(d), 100 * mean(d < 0)))

cat("\n== 4. signal ordering: BIC(mean first) - BIC(contrast first), best of 10 runs each ==\n")
d <- f4j$BIC - f4k$BIC; t <- t.test(d)
cat(sprintf("  mean %+.3f | t(%.0f) = %.2f p = %.2g | contrast first lower for %.1f%%\n",
            mean(d), t$parameter, t$statistic, t$p.value, 100 * mean(d > 0)))
sink()
