# Model-free behaviour reported in Results (acceptance, RT, CES-D, trial-level choice and
# RT models). Section A of memo_analyses.R, re-pointed at the best-of-10 canonical fits
# (used only to define the 351 participants). Run from the repository root. 2026-09-27
suppressPackageStartupMessages({ library(lme4) })
if (!file.exists("config.R")) stop("run from the repository root")
source("config.R"); source(file.path(MODEL_DIR, "tDDM", "tddm_engine.R"))
R <- "results/tDDM/"
wt <- function(x, g, lab) {  # Welch HC vs DEP, DEP-HC
  t <- t.test(x[g == "DEP"], x[g == "HC"])
  sp <- sqrt((var(x[g=="HC"]) + var(x[g=="DEP"]))/2)
  cat(sprintf("  %-28s HC=%.3f DEP=%.3f diff=%+.3f t(%.0f)=%.2f p=%.4f d=%+.2f\n", lab,
      mean(x[g=="HC"]), mean(x[g=="DEP"]), mean(x[g=="DEP"])-mean(x[g=="HC"]),
      t$parameter, t$statistic, t$p.value, (mean(x[g=="DEP"])-mean(x[g=="HC"]))/sp))
}
fit4j <- read.csv(paste0(R, "canonical_4j_accurate_bestof10/fits_4j_accurate_best10.csv"))
g <- fit4j$group
cat("N =", nrow(fit4j), table(g), "\n")
## ---------------- A. model-free ----------------
cat("\n==== A. Model-free behaviour ====\n")
d <- tddm_load_study4(DATA_DIR)
d <- d[d$subject %in% fit4j$subject, ]
ps <- do.call(rbind, lapply(split(d, d$subject), function(s) data.frame(
  subject = s$subject[1], group = s$group[1], n = nrow(s), acc = mean(s$choice == 0),
  medRT = median(s$RT_ddm), meanRT = mean(s$RT_ddm), ces = s$cesScore[1],
  medRT_acc = median(s$RT_ddm[s$choice == 0]), medRT_rej = median(s$RT_ddm[s$choice == 1]))))
cat(sprintf("  trials/participant: median %d, range %d-%d\n", median(ps$n), min(ps$n), max(ps$n)))
wt(ps$acc, ps$group, "acceptance rate")
wt(ps$medRT, ps$group, "median RT (s)")
wt(ps$meanRT, ps$group, "mean RT (s)")
wt(log(ps$medRT), ps$group, "log median RT")
w <- wilcox.test(medRT ~ group, ps); cat(sprintf("  Wilcoxon median RT p=%.4f\n", w$p.value))
cat(sprintf("  CES-D ~ medRT r=%.3f p=%.4f ; CES-D ~ acc r=%.3f p=%.4f\n",
    cor(ps$ces, ps$medRT), cor.test(ps$ces, ps$medRT)$p.value, cor(ps$ces, ps$acc), cor.test(ps$ces, ps$acc)$p.value))
# trial-level: choice sensitivity to the two signals, by group
d$DEP <- as.numeric(d$group == "DEP"); d$acc <- as.numeric(d$choice == 0)
m1 <- glmer(acc ~ (mean_sig + contrast_sig) * DEP + (1 + mean_sig | subject), d, family = binomial,
            control = glmerControl(optimizer = "bobyqa"))
cat("  glmer acc ~ (mean + contrast) * DEP:\n"); print(round(coef(summary(m1)), 4))
# DEP simple slopes: refit with HC as the indicator so the main effects are the DEP slopes
d$HC <- 1 - d$DEP
m1h <- glmer(acc ~ (mean_sig + contrast_sig) * HC + (1 + mean_sig | subject), d, family = binomial,
             control = glmerControl(optimizer = "bobyqa"))
cat("  DEP simple slopes (HC-indicator refit):\n"); print(round(coef(summary(m1h))[c("mean_sig", "contrast_sig"), ], 4))
m2 <- lmer(log(RT_ddm) ~ (mean_sig + contrast_sig) * DEP + (1 | subject), d)
cat("  lmer logRT ~ (mean + contrast) * DEP:\n"); print(round(coef(summary(m2)), 4))
d$absmean <- abs(d$mean_sig); d$abscon <- abs(d$contrast_sig)
m3 <- lmer(log(RT_ddm) ~ (absmean + abscon) * DEP + (1 | subject), d)
cat("  lmer logRT ~ (|mean| + |contrast|) * DEP:\n"); print(round(coef(summary(m3)), 4))

