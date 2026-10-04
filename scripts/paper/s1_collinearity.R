#!/usr/bin/env Rscript
# ==============================================================================
# Paper numbers: Supplement S1 collinearity of the two prediction errors
# ==============================================================================
# "Across all trials and participants the two prediction errors correlate at
#  r = .757. Within participants the median correlation is .854 (IQR .741 to
#  .918), which gives a median variance inflation factor of 3.68 (IQR 2.22 to
#  6.37) for the two-signal model. The factor exceeds 5 for 36.5% of
#  participants and 10 for 11.7%."
#
# Trials: engine loader tddm_load_study4() (0.15 < RT <= 10 s, complete cases),
# restricted to the 351 fitted participants. Signals: RPE_scale, EPE_X1_scale.
# VIF with two predictors = 1 / (1 - r^2), r = within-participant correlation.
# Quantiles use R's default (type 7).
#
# Usage (from repo root):  Rscript scripts/paper/s1_collinearity.R
# Output: results/paper_numbers/s1_collinearity.csv, s1_collinearity_per_subject.csv
# ==============================================================================
source("config.R")
source(file.path(MODEL_DIR, "tDDM", "tddm_engine.R"))
OUT <- file.path(RESULTS_DIR, "paper_numbers"); dir.create(OUT, showWarnings = FALSE, recursive = TRUE)

fits <- read.csv(file.path(RESULTS_DIR, "tDDM", "canonical_4j_accurate_bestof10", "fits_4j_accurate_best10.csv"))
d_all <- tddm_load_study4(DATA_DIR)
d <- d_all[d_all$subject %in% fits$subject, ]
cat(sprintf("Retained trials: %d from %d participants (loader returned %d participants)\n",
            nrow(d), length(unique(d$subject)), length(unique(d_all$subject))))

r_pooled <- cor(d$RPE_scale, d$EPE_X1_scale)
cat(sprintf("Pooled r(RPE, EPE_X1), retained trials: %.3f\n", r_pooled))

ps <- do.call(rbind, lapply(split(d, d$subject), function(s)
  data.frame(subject = s$subject[1], group = s$group[1], n = nrow(s),
             r = cor(s$RPE_scale, s$EPE_X1_scale))))
ps$vif <- 1 / (1 - ps$r^2)
n_na <- sum(!is.finite(ps$r))
if (n_na) cat("Participants with undefined r (constant signal):", n_na, "\n")
ok <- is.finite(ps$r)
qr <- quantile(ps$r[ok], c(.25, .5, .75)); qv <- quantile(ps$vif[ok], c(.25, .5, .75))
cat(sprintf("Within-participant r: median %.3f (IQR %.3f to %.3f), N = %d\n", qr[2], qr[1], qr[3], sum(ok)))
cat(sprintf("VIF:                  median %.2f (IQR %.2f to %.2f)\n", qv[2], qv[1], qv[3]))
cat(sprintf("VIF > 5: %.1f%%   VIF > 10: %.1f%%\n", 100 * mean(ps$vif[ok] > 5), 100 * mean(ps$vif[ok] > 10)))

res <- data.frame(
  quantity = c("r_pooled", "r_within_median", "r_within_q25", "r_within_q75",
               "vif_median", "vif_q25", "vif_q75", "pct_vif_gt5", "pct_vif_gt10", "n_participants", "n_trials"),
  value = c(r_pooled, qr[2], qr[1], qr[3], qv[2], qv[1], qv[3],
            100 * mean(ps$vif[ok] > 5), 100 * mean(ps$vif[ok] > 10), sum(ok), nrow(d)))
write.csv(res, file.path(OUT, "s1_collinearity.csv"), row.names = FALSE)
write.csv(ps, file.path(OUT, "s1_collinearity_per_subject.csv"), row.names = FALSE)
cat("Saved:", file.path(OUT, "s1_collinearity.csv"), "\n")
