#!/usr/bin/env Rscript
# ==============================================================================
# DEMOGRAPHICS TABLE — HC vs DEP
# ==============================================================================
# Companion to the group comparison: are the two CES-D-defined groups matched on
# age and gender, and how do sample-level task quantities (usable trials, median
# choice RT) compare? Reviewers ask for this whenever groups are defined by a
# self-report cutoff rather than randomisation.
#
# Reports per group: n, age (mean/SD, range), gender counts + % female, CES-D
# (mean/SD, range), usable trials after the fitting RT filter, and median RT.
# Tests: Welch t for continuous, chi-square for gender.
#
# Usage:
#   Rscript scripts/tDDM/demographics_table.R
#   Rscript scripts/tDDM/demographics_table.R --fits results/tDDM/<dir>/fits.csv
# ==============================================================================
args <- commandArgs(trailingOnly = TRUE)
ga_s <- function(f, d) { if (f %in% args) { i <- which(args == f); if (i < length(args)) return(args[i + 1]) }; d }

script_dir <- tryCatch(dirname(sys.frame(1)$ofile), error = function(e) {
  ca <- commandArgs(FALSE); fa <- sub("^--file=", "", ca[grep("^--file=", ca)])
  if (length(fa)) dirname(fa) else getwd()
})
source(file.path(script_dir, "..", "..", "config.R"))

RT_LO <- 0.15; RT_HI <- 10   # must match the fitting filter
OUT_DIR <- ga_s("--outdir", "results/tDDM/cohort_demographics")

hc <- read.csv(file.path(DATA_DIR, "study4_HC.csv")); hc$grp <- "HC"
dep <- read.csv(file.path(DATA_DIR, "study4_dep.csv")); dep$grp <- "DEP"
all <- rbind(hc, dep)

# one row per subject
subj <- do.call(rbind, lapply(split(all, all$sub), function(d) {
  keep <- is.finite(d$choice_RT) & d$choice_RT > RT_LO & d$choice_RT <= RT_HI
  data.frame(sub = d$sub[1], grp = d$grp[1],
             age = d$age[1], gender = as.character(d$gender[1]),
             ces = d$cesScore[1],
             n_trials_raw = nrow(d), n_trials_used = sum(keep),
             medRT = median(d$choice_RT[keep]),
             stringsAsFactors = FALSE)
}))
HC <- subj$grp == "HC"; DEP <- subj$grp == "DEP"

pz <- function(p) if (is.na(p)) "NA" else if (p < .001) "<.001" else sprintf("%.3f", p)
ms <- function(x) sprintf("%.2f (%.2f)", mean(x, na.rm = TRUE), sd(x, na.rm = TRUE))
rg <- function(x) sprintf("[%g, %g]", min(x, na.rm = TRUE), max(x, na.rm = TRUE))

cat(sprintf("\nN = %d  (HC = %d, DEP = %d)   CES-D cutoff = 16\n\n",
            nrow(subj), sum(HC), sum(DEP)))
cat(sprintf("Gender codes present: %s\n",
            paste(sort(unique(subj$gender)), collapse = ", ")))

rows <- list()
cont <- c(age = "Age (years)", ces = "CES-D total",
          n_trials_used = "Usable trials", medRT = "Median choice RT (s)")
for (v in names(cont)) {
  x <- subj[[v]]
  tt <- t.test(x[DEP], x[HC])
  d <- (mean(x[DEP], na.rm = TRUE) - mean(x[HC], na.rm = TRUE)) / sd(x, na.rm = TRUE)
  rows[[v]] <- data.frame(
    variable = cont[[v]],
    HC = paste(ms(x[HC]), rg(x[HC])), DEP = paste(ms(x[DEP]), rg(x[DEP])),
    test = sprintf("t(%.0f) = %.2f", tt$parameter, tt$statistic),
    p = pz(tt$p.value), cohen_d = sprintf("%+.2f", d), row.names = NULL)
}

# gender: chi-square over observed categories
gt <- table(subj$grp, subj$gender)
cs <- suppressWarnings(chisq.test(gt))
fem <- grep("^(f|female|2)$", colnames(gt), ignore.case = TRUE, value = TRUE)
pct <- function(g) if (length(fem)) sprintf("%d/%d female (%.1f%%)",
    sum(gt[g, fem]), sum(gt[g, ]), 100 * sum(gt[g, fem]) / sum(gt[g, ])) else
    paste(sprintf("%s=%d", colnames(gt), gt[g, ]), collapse = ", ")
rows[["gender"]] <- data.frame(
  variable = "Gender", HC = pct("HC"), DEP = pct("DEP"),
  test = sprintf("chi2(%d) = %.2f", cs$parameter, cs$statistic),
  p = pz(cs$p.value), cohen_d = "-", row.names = NULL)

out <- do.call(rbind, rows[c("age", "gender", "ces", "n_trials_used", "medRT")])
print(out, row.names = FALSE)

cat("\nGender x group contingency:\n"); print(gt)

# ---- do the group effects survive adjusting for age and gender? --------------
# The groups are NOT age-matched, so every group contrast must be re-checked with
# age (and gender) in the model. Note the direction: age correlates POSITIVELY
# with threshold and nDT while DEP is YOUNGER, so the unadjusted contrast is if
# anything conservative -- adjustment should strengthen, not weaken, the effects.
# 2026-08-14: default switched to the CORRECTED-transition-density fits (see
# 260808_LABMEETING_FEEDBACK.md). The demographic table itself does not depend on
# the fits; only this adjustment block does.
FITS_REL <- ga_s("--fits", "results/tDDM/canonical_4j_accurate_bestof10/fits_4j_accurate_best10.csv")
fits_path <- file.path(PROJECT_ROOT, FITS_REL)
if (file.exists(fits_path)) {
  fits <- read.csv(fits_path)
  m <- merge(fits, subj, by.x = "subject", by.y = "sub")
  m$grpDEP <- as.numeric(m$grp == "DEP")
  cat(sprintf("\n=== Group contrasts adjusted for age + gender (%s) ===\n", FITS_REL))
  cat("param       r(age)   p      | unadj d    p      | adj b(group)  p\n")
  adj <- do.call(rbind, lapply(c("tau", "thres", "nDT", "d_mean", "d_contrast"),
    function(v) {
      Y <- m[[v]]
      ra <- cor.test(Y, m$age)
      tt <- t.test(Y[m$grpDEP == 1], Y[m$grpDEP == 0])
      d <- (mean(Y[m$grpDEP == 1]) - mean(Y[m$grpDEP == 0])) / sd(Y)
      co <- summary(lm(Y ~ grpDEP + age + gender, m))$coefficients["grpDEP", ]
      cat(sprintf("%-11s %+.3f  %-6s | %+.2f      %-6s | %+.4f      %-6s\n",
                  v, ra$estimate, pz(ra$p.value), d, pz(tt$p.value), co[1], pz(co[4])))
      data.frame(parameter = v, r_age = ra$estimate, p_age = ra$p.value,
                 unadj_d = d, unadj_p = tt$p.value,
                 adj_beta_group = co[1], adj_p_group = co[4], row.names = NULL)
    }))
  cat("\nNOTE: adjustment does not bear on whether the onset contrast is timing-\n",
      "specific -- that is settled by the generative controls, not by covariates.\n", sep = "")
} else {
  adj <- NULL
  cat(sprintf("\n[skip] fits not found: %s\n", fits_path))
}

stamp <- format(Sys.time(), "%Y%m%d_%H%M%S")
outdir <- file.path(PROJECT_ROOT, OUT_DIR); dir.create(outdir, showWarnings = FALSE, recursive = TRUE)
f <- file.path(outdir, sprintf("demographics_table_%s.csv", stamp))
write.csv(out, f, row.names = FALSE)
cat(sprintf("\nSaved: %s\n", f))
if (!is.null(adj)) {
  f2 <- file.path(outdir, sprintf("demographics_adjusted_%s.csv", stamp))
  write.csv(adj, f2, row.names = FALSE)
  cat(sprintf("Saved: %s\n", f2))
}
