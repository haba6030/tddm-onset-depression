suppressPackageStartupMessages({})
if (!file.exists("config.R")) stop("run from the repository root"); source("config.R"); source(file.path(MODEL_DIR,"tDDM","tddm_engine.R"))
f <- read.csv("results/tDDM/canonical_4j_accurate_bestof10/fits_4j_accurate_best10.csv")
d <- tddm_load_study4(DATA_DIR); d <- d[d$subject %in% f$subject, ]
od <- "results/tDDM/memo_analyses_0919/starting_point_sensitivity_20260927"; dir.create(od, recursive = TRUE, showWarnings = FALSE)
write.csv(d[, c("subject","group","mean_sig","contrast_sig","choice","choice_ddm","RT_ddm")], file.path(od, "trials.csv"), row.names=FALSE)
cat("TDDM_SIGMA", TDDM_SIGMA, "\n"); print(table(d$choice, d$choice_ddm))
ps <- aggregate(choice ~ subject + group, d, function(x) mean(x==0))
t <- t.test(ps$choice[ps$group=="HC"], ps$choice[ps$group=="DEP"]); print(t)
