#!/usr/bin/env bash
# Refit everything from the raw data, in the order the paper's results depend on.
# About 14 hours of wall-clock time on 8 cores with the C++ likelihood. Each step writes
# a new timestamped directory under results/tDDM/.
set -euo pipefail
cd "$(dirname "$0")"
C=${CORES:-8}
[ -f data/study4_HC.csv ] || python data/fetch_data.py
stamp () { ls -dt results/tDDM/"$1"_* | head -1; }

# main model: 10 runs per signal ordering, keep the lowest-NLL run per participant (~20 min each)
Rscript scripts/tDDM/optimizer_reproducibility.R --reps 10 --cores "$C"
Rscript scripts/tDDM/memo_0919/build_bestof_runs.R --runs-dir "$(stamp optimizer_reproducibility_mean_r1)" \
  --out results/tDDM/refit_bestof10/fits_4j_mean_first_best10.csv
FITS=results/tDDM/refit_bestof10/fits_4j_mean_first_best10.csv

# baselines (<1 h each)
Rscript scripts/tDDM/refit_corrected_likelihood.R --arms 4a_static1
Rscript scripts/tDDM/refit_corrected_likelihood.R --arms 4c_static

# simulation controls, contrast-first refit and split-half (~13 h on 8 cores).
# The runner's operating points (--tau-common 0.891, --tau-base 0.804) are the pooled and
# HC mean onsets of the paper's fits with tau <= 0 set to 0; recompute them for new fits.
FITS="$FITS" CORES="$C" bash scripts/tDDM/memo_0919/run_stage6_accurate_controls.sh
Rscript scripts/tDDM/memo_0919/build_bestof_runs.R --runs-dir "$(stamp optimizer_reproducibility_contrast_r1)" \
  --out results/tDDM/refit_bestof10/fits_4j_contrast_first_best10.csv

# checks of the likelihood and the simulator (minutes each, except the brute-force check)
Rscript scripts/tDDM/verify_transition_density.R
Rscript scripts/tDDM/verify_cpp_likelihood.R
Rscript scripts/tDDM/diag_simulator_vs_likelihood.R
Rscript scripts/tDDM/diag_quadrature_spectral.R
Rscript scripts/tDDM/diag_quadrature_brute.R
echo "Fitting done. The summary and figure scripts read the shipped run directories;"
echo "point them at the new ones (README section 4) and run run_fast.sh."
