#!/usr/bin/env Rscript
# ==============================================================================
# CES-D CONTINUOUS CORRELATION with tDDM parameters
# ==============================================================================
# Robustness companion to the binary HC/DEP group comparison: instead of
# splitting at the CES-D cutoff (16), treat depressive symptom load as a
# CONTINUOUS predictor (CES-D 0-48) across all subjects. Answers:
#   (a) Does each fitted tDDM parameter correlate with CES-D? (raw)
#   (b) Does that correlation survive controlling median RT? (the RT-confound
#       channel that manufactures the onset tau effect)
#   (c) Are the continuous results consistent with the binary group contrasts?
#
# Reads per-subject tDDM fits + per-subject CES-D + median choice_RT, merges,
# and reports for each parameter: r(param, CES), raw beta_CES + p, and the
# medRT-adjusted beta_CES + p. Also prints the binary DEP-HC contrast for
# reference. Writes a tidy CSV to results/tDDM/.
#
# NOTE on median-RT control: for tau, medRT is a CONFOUND (tau is non-identified;
# its apparent CES association runs through RT censoring). For threshold/nDT,
# medRT is a MECHANICAL CONSEQUENCE (higher threshold -> longer RT), so adjusting
# for it is OVER-CONTROL and does not by itself imply artifact. Distinguishing
# artifact from over-controlled-real requires the generative positive control
# (scripts/tDDM/positive_control_tau*.R), not this regression.
#
# Usage:
#   Rscript scripts/tDDM/ces_continuous_correlation.R
#   Rscript scripts/tDDM/ces_continuous_correlation.R --fits results/tDDM/<dir>/fits.csv
# ==============================================================================
args <- commandArgs(trailingOnly = TRUE)
ga_s <- function(f, d) { if (f %in% args) { i <- which(args == f); if (i < length(args)) return(args[i + 1]) }; d }

script_dir <- tryCatch(dirname(sys.frame(1)$ofile), error = function(e) {
  ca <- commandArgs(FALSE); fa <- sub("^--file=", "", ca[grep("^--file=", ca)])
  if (length(fa)) dirname(fa) else getwd()
})
source(file.path(script_dir, "..", "..", "config.R"))

# 2026-08-14: default switched to the CORRECTED-transition-density fits. The old
# default (refit_gateB_20260624_172551/fits_4j_widenDT.csv) was produced with the
# misspecified split-at-tau density (no absorption, sigma-scaled SD) and its tau
# estimates are inflated ~2.5x (mean 1.800 -> 0.720). See 260808_LABMEETING_FEEDBACK.md.
FITS_REL   <- ga_s("--fits", "results/tDDM/canonical_4j_accurate_bestof10/fits_4j_accurate_best10.csv")
OUT_DIR    <- ga_s("--outdir", "results/tDDM/cohort_demographics")
PARAMS     <- c("tau", "thres", "nDT", "d_mean", "d_contrast")  # tDDM params to test
RT_LO <- 0.15; RT_HI <- 10                                       # choice_RT filter (matches fitting)

# ---- load + merge ------------------------------------------------------------
fits <- read.csv(file.path(PROJECT_ROOT, FITS_REL), stringsAsFactors = FALSE)
fits$tau <- pmax(fits$tau, 0)   # 2026-09-30: every tau <= 0 is the contrast model
hc <- read.csv(file.path(DATA_DIR, "study4_HC.csv")); hc$grp <- "HC"
dep <- read.csv(file.path(DATA_DIR, "study4_dep.csv")); dep$grp <- "DEP"
all <- rbind(hc, dep)

agg <- do.call(rbind, lapply(split(all, all$sub), function(d) {
  rt <- d$choice_RT[is.finite(d$choice_RT) & d$choice_RT > RT_LO & d$choice_RT <= RT_HI]
  data.frame(sub = d$sub[1], grp = d$grp[1], ces = d$cesScore[1],
             ces_within_sd = sd(d$cesScore, na.rm = TRUE), medRT = median(rt))
}))
m <- merge(fits, agg, by.x = "subject", by.y = "sub")
HC <- m$grp == "HC"; DEP <- m$grp == "DEP"
pz <- function(p) if (is.na(p)) "NA" else if (p < .001) "<.001" else sprintf("%.3f", p)

cat(sprintf("Fits: %s\n", FITS_REL))
cat(sprintf("N=%d | CES-D within-subject SD max=%.3f (0 = constant per subject, OK)\n",
            nrow(m), max(agg$ces_within_sd, na.rm = TRUE)))
cat(sprintf("CES-D: mean=%.1f sd=%.1f range=[%.0f,%.0f] | HC mean=%.1f (max=%.0f)  DEP mean=%.1f (min=%.0f)\n",
            mean(m$ces), sd(m$ces), min(m$ces), max(m$ces),
            mean(m$ces[HC]), max(m$ces[HC]), mean(m$ces[DEP]), min(m$ces[DEP])))
cat(sprintf("cor(CES-D, medRT) = %.3f   [the RT-confound channel]\n\n", cor(m$ces, m$medRT)))

# ---- per-parameter continuous + adjusted + binary ----------------------------
rows <- lapply(PARAMS, function(y) {
  Y <- m[[y]]
  r    <- cor(Y, m$ces)
  c1   <- summary(lm(Y ~ ces, m))$coefficients["ces", ]
  c2   <- summary(lm(Y ~ ces + medRT, m))$coefficients["ces", ]
  tt   <- t.test(Y[DEP], Y[HC])
  dcoh <- (mean(Y[DEP]) - mean(Y[HC])) / sd(Y)
  data.frame(parameter = y,
             r_ces = r,
             beta_ces_raw = c1[1], p_ces_raw = c1[4],
             beta_ces_adjRT = c2[1], p_ces_adjRT = c2[4],
             binary_delta_DEPminusHC = mean(Y[DEP]) - mean(Y[HC]),
             binary_cohen_d = dcoh, binary_p = tt$p.value,
             row.names = NULL)
})
res <- do.call(rbind, rows)

cat("=== Continuous CES-D (raw and median-RT-adjusted) vs binary DEP-HC ===\n")
for (i in seq_len(nrow(res))) with(res[i, ], cat(sprintf(
  "%-11s r=%+.3f | raw b=%+.4f p=%-6s | +medRT b=%+.4f p=%-6s || binary d=%+.2f p=%s\n",
  parameter, r_ces, beta_ces_raw, pz(p_ces_raw), beta_ces_adjRT, pz(p_ces_adjRT),
  binary_cohen_d, pz(binary_p))))

# ---- write ------------------------------------------------------------------
ts <- format(Sys.time(), "%Y%m%d_%H%M%S")
outdir <- file.path(PROJECT_ROOT, OUT_DIR); dir.create(outdir, showWarnings = FALSE, recursive = TRUE)
out <- file.path(outdir, sprintf("ces_continuous_correlation_%s.csv", ts))
write.csv(res, out, row.names = FALSE)
cat(sprintf("\nSaved: %s\n", out))
