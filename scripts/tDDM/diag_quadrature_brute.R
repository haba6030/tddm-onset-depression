# Adjudicates the V1 configurations where the accurate integral and the spectral
# reference disagree (> 1e-4 relative), using a 1.2e6-point trapezoid of the same
# integrand with log-spaced layers at both walls. Run after diag_quadrature_spectral.R.
#   Rscript scripts/tDDM/diag_quadrature_brute.R
if (!file.exists("config.R")) stop("run from the repository root")
source("config.R"); source(file.path(MODEL_DIR, "tDDM", "tddm_engine.R"))
OD <- tail(sort(Sys.glob("results/tDDM/diag_quadrature_spectral_*")), 1)
V <- read.csv(file.path(OD, "V1_random_configs.csv"))
V$ca <- abs(V$v1-V$v2)*V$a; V$re <- abs(V$acc-V$sp)/V$sp
b <- V[V$re > 1e-4, ]; b <- b[order(b$ca), ]
brute <- function(v1, v2, tau, a, rem, ch) {
  # 1e6-point trapezoid plus a 1e5-point log-spaced layer at each wall
  e <- 10^seq(-9, log10(a/2), length.out = 1e5)
  xg <- sort(unique(c(seq(0, a, length.out = 1e6), e, a - e)))
  y <- p_trans_image(xg, tau, a/2, v1, a) *
    (if (ch == 1) dwiener_nf_vec(rem, -v2, a, 1 - xg/a) else dwiener_nf_vec(rem, v2, a, xg/a))
  sum((y[-1] + y[-length(y)]) / 2 * diff(xg))
}
b$brute <- mapply(brute, b$v1, b$v2, b$tau, b$a, b$rem, b$ch)
b$re_acc_brute <- abs(b$acc - b$brute)/b$brute; b$re_sp_brute <- abs(b$sp - b$brute)/b$brute
print(signif(b[, c("v1","v2","a","tau","rem","ca","sp","acc","brute","re_acc_brute","re_sp_brute")], 3))
write.csv(b, file.path(OD, "V1b_brute_adjudication.csv"), row.names = FALSE)
cat(sprintf("n = %d | accurate vs brute max rel err %.1e | spectral vs brute max rel err %.1e\n", nrow(b), max(b$re_acc_brute), max(b$re_sp_brute)))
